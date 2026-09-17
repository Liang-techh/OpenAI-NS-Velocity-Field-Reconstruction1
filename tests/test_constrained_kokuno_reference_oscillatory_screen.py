from __future__ import annotations

import numpy as np

from openai_ns_reconstruction.kokuno_complete_curl import KokunoCompleteCurlCorrection
from openai_ns_reconstruction.kokuno_reference_continuation import (
    KokunoReferenceContinuationCandidate,
)
from openai_ns_reconstruction.kokuno_reference_oscillatory_screen import (
    AMPLITUDE,
    PHASE,
    _sample_points,
    run_screen,
)


def test_post_core_stratum_is_new_reference_continuation_domain() -> None:
    leading = KokunoReferenceContinuationCandidate()
    samples = _sample_points(count=6, seed=9172981)
    points = samples["post_core"]
    time = 0.375

    coordinates = leading.coordinates(points[:, 0], points[:, 1], points[:, 2], time)
    assert np.min(coordinates["X"]) > 4.1 / leading.Lambda

    with np.testing.assert_raises((ValueError, RuntimeError)):
        leading.core.at_points(points, time)

    values = leading.at_points(points, time)
    assert values.shape == (6, 3)
    assert np.all(np.isfinite(values))
    assert np.max(np.linalg.norm(values, axis=1)) > 0.0


def test_frozen_complete_curl_is_divergence_free_and_composable() -> None:
    leading = KokunoReferenceContinuationCandidate()
    correction = KokunoCompleteCurlCorrection(amplitude=AMPLITUDE, phase=PHASE)
    points = _sample_points(count=6, seed=9172983)["inner"]
    time = 0.5

    base = leading.at_points(points, time)
    delta = correction.at_points(points, time)
    total = base + delta
    divergence = correction.divergence(
        points[:, 0], points[:, 1], points[:, 2], time
    )

    assert base.shape == delta.shape == total.shape == (6, 3)
    assert np.all(np.isfinite(total))
    assert np.max(np.abs(divergence)) < 1.0e-12
    assert np.max(np.linalg.norm(delta, axis=1)) > 0.0


def test_screen_is_fail_closed_and_performs_no_selection() -> None:
    report = run_screen(
        count=4,
        seed=9172985,
        times=(0.5,),
        steps=(0.01,),
    )

    assert report["task_id"] == "K2-OSC-005"
    assert report["validation"]["post_core_stratum_old_core_rejected"] is True
    assert report["validation"]["post_core_X_min"] > report["validation"]["old_core_nominal_X_limit"]
    assert report["oscillation"]["amplitude"] == AMPLITUDE
    assert report["oscillation"]["phase"] == PHASE
    assert report["oscillation"]["parameter_selection_performed"] is False
    assert report["routing"]["routing_changed_by_this_screen"] is False
    assert report["routing"]["candidate_promoted"] is False
    assert report["fixed_gates"]["changed"] is False
    assert report["truth_boundary"]["formal_full_domain_pde_gate_assessed"] is False
    assert report["truth_boundary"]["pde_validated"] is False
    assert report["truth_boundary"]["paper_exact"] is False
    assert report["truth_boundary"]["openai_field_identified"] is False
    assert len(report["rows"]) == 2
    for row in report["rows"]:
        assert np.isfinite(row["residual_rms_ratio"])
        assert row["residual_rms_ratio"] > 0.0
        assert row["analytic_correction_divergence_max"] < 1.0e-12
