import json

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_core_composite_checkpoint import (
    DIAGNOSTIC_OSCILLATORY_AMPLITUDE,
    PDE_THRESHOLD,
    KokunoCoreCompositeCandidate,
    generate_core_composite_checkpoint,
)


def test_core_composite_velocity_round_trip_and_grid(tmp_path):
    candidate = KokunoCoreCompositeCandidate()
    points = np.array([
        [0.0, 0.0, 0.0],
        [0.08, -0.04, 0.06],
        [-0.06, 0.07, -0.05],
    ])
    before = candidate.at_points(points, 0.5)
    path = candidate.save_json(tmp_path / "candidate.json")
    replay = KokunoCoreCompositeCandidate.load_json(path)
    after = replay.at_points(points, 0.5)
    np.testing.assert_allclose(after, before, rtol=0.0, atol=0.0)
    assert replay.sha256 == candidate.sha256
    assert replay.oscillatory.amplitude == DIAGNOSTIC_OSCILLATORY_AMPLITUDE

    grid = replay.grid([-0.04, 0.04], [-0.04, 0.04], [-0.04, 0.04], [0.5])
    assert grid.shape == (1, 2, 2, 2, 3)
    assert np.isfinite(grid).all()


def test_core_composite_artifact_fails_closed_on_truth_laundering(tmp_path):
    candidate = KokunoCoreCompositeCandidate()
    path = candidate.save_json(tmp_path / "candidate.json")
    payload = json.loads(path.read_text())
    payload["truth_boundary"]["pde_validated"] = True
    path.write_text(json.dumps(payload))
    with pytest.raises(ValueError):
        KokunoCoreCompositeCandidate.load_json(path)


def test_integration_checkpoint_keeps_full_gate_unassessed(tmp_path):
    report = generate_core_composite_checkpoint(
        output_dir=tmp_path / "checkpoint", seed=9172861, point_count=4
    )
    assert report["registered_gate"]["held_out_normalized_ns_residual_threshold"] == PDE_THRESHOLD
    assert report["registered_gate"]["assessed"] is False
    assert report["registered_gate"]["passed"] is False
    assert report["stage_state"]["leading_core_ready"] is True
    assert report["stage_state"]["leading_ready"] is False
    assert report["stage_state"]["oscillatory_ready"] is True
    assert report["stage_state"]["correction_ready"] is False
    assert report["stage_state"]["core_velocity_export_ready"] is True
    assert report["stage_state"]["velocity_export_ready"] is False
    assert report["stage_state"]["pde_validated"] is False
    assert report["agent3_bridge_cycle_evidence"]["mean_correction_applied_to_core_composite"] is False
    assert report["sampling"]["forcing"].startswith("zero diagnostic")
    assert len(report["baseline_vs_core_composite"]["rows"]) == 9
    for row in report["baseline_vs_core_composite"]["rows"]:
        assert np.isfinite(row["momentum_rms_ratio"])
        assert row["leading_only"]["momentum_rms"] > 0.0
        assert row["leading_plus_oscillatory"]["momentum_rms"] > 0.0
