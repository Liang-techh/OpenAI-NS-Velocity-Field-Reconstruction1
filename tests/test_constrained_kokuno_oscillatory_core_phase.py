from math import pi

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_complete_curl import KokunoCompleteCurlCorrection
from openai_ns_reconstruction.kokuno_core_composite_checkpoint import KokunoCoreCompositeCandidate
from openai_ns_reconstruction.kokuno_leading_core_series import KokunoLeadingCoreSeriesCandidate
from openai_ns_reconstruction.kokuno_oscillatory_core_phase import (
    AMPLITUDE,
    PHASES,
    _sample_core_points,
    _select_training_phase,
    evaluate_phase_grid,
    generate_core_phase_report,
)


def test_phase_grid_is_bounded_and_changes_public_velocity():
    assert len(PHASES) == 8
    assert all(-pi <= phase <= pi for phase in PHASES)
    assert 0.0 in PHASES

    leading = KokunoLeadingCoreSeriesCandidate()
    point = np.array([[0.071, -0.043, 0.087]])
    zero_phase = KokunoCompleteCurlCorrection(amplitude=AMPLITUDE, phase=0.0)
    quarter_phase = KokunoCompleteCurlCorrection(amplitude=AMPLITUDE, phase=pi / 2.0)
    u0 = leading.velocity(point[:, 0], point[:, 1], point[:, 2], 0.5) + zero_phase.at_points(point, 0.5)
    u1 = leading.velocity(point[:, 0], point[:, 1], point[:, 2], 0.5) + quarter_phase.at_points(point, 0.5)
    assert np.isfinite(u0).all() and np.isfinite(u1).all()
    assert not np.allclose(u0, u1)


def test_small_phase_screen_runs_with_frozen_amplitude():
    leading = KokunoLeadingCoreSeriesCandidate()
    points = _sample_core_points(9172881, 6)
    rows = evaluate_phase_grid(
        leading,
        points,
        phases=(0.0, pi / 2.0),
        times=(0.5,),
        steps=(0.01,),
    )
    assert len(rows) == 2
    assert {row["phase"] for row in rows} == {0.0, pi / 2.0}
    assert all(row["momentum_rms"] > 0.0 for row in rows)
    assert all(np.isfinite(row["momentum_max"]) for row in rows)

    with pytest.raises(ValueError, match="freezes amplitude"):
        evaluate_phase_grid(
            leading,
            points,
            phases=(0.0,),
            times=(0.5,),
            steps=(0.01,),
            template=KokunoCompleteCurlCorrection(amplitude=0.25),
        )


def test_training_selector_is_deterministic_and_not_holdout_aware():
    aggregates = {
        "a": {"phase": -pi / 2.0, "mean_momentum_rms": 2.0},
        "b": {"phase": pi / 2.0, "mean_momentum_rms": 1.0},
        "c": {"phase": 0.0, "mean_momentum_rms": 1.0},
    }
    assert _select_training_phase(aggregates) == 0.0


def test_report_separates_train_holdout_and_roundtrips_candidate(tmp_path):
    report = generate_core_phase_report(
        output_dir=tmp_path,
        train_seed=9172881,
        holdout_seed=9172891,
        train_count=8,
        holdout_count=10,
    )
    assert report["screen"]["training_seed"] != report["screen"]["holdout_seed"]
    assert report["training"]["selected_phase"] in PHASES
    assert report["holdout"]["mean_selected_phase_momentum_rms"] > 0.0
    assert report["holdout"]["mean_phase0_momentum_rms"] > 0.0
    assert np.isfinite(report["holdout"]["selected_over_phase0"])
    assert report["truth_boundary"]["formal_full_domain_pde_gate_assessed"] is False
    assert report["truth_boundary"]["pde_validated"] is False
    assert report["truth_boundary"]["paper_exact"] is False

    replay = KokunoCoreCompositeCandidate.load_json(tmp_path / "selected_phase_core_candidate.json")
    assert replay.sha256 == report["candidate"]["sha256"]
    grid = replay.grid(
        np.array([-0.05, 0.0, 0.05]),
        np.array([-0.05, 0.0, 0.05]),
        np.array([-0.05, 0.0, 0.05]),
        np.array([0.4, 0.6]),
    )
    assert grid.shape == (2, 3, 3, 3, 3)
    assert np.isfinite(grid).all()
    assert np.max(np.abs(grid)) > 0.0


def test_train_holdout_seed_alias_is_rejected(tmp_path):
    with pytest.raises(ValueError, match="must differ"):
        generate_core_phase_report(
            output_dir=tmp_path,
            train_seed=123,
            holdout_seed=123,
            train_count=8,
            holdout_count=8,
        )
