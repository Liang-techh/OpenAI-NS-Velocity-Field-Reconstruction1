import json

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_reference_composite_checkpoint import (
    REGISTERED_DIVERGENCE_THRESHOLD,
    REGISTERED_PDE_THRESHOLD,
    KokunoReferenceCompositeCandidate,
    build_checkpoint,
    write_checkpoint,
)


def test_reference_composite_is_exact_public_sum_and_roundtrips(tmp_path):
    candidate = KokunoReferenceCompositeCandidate()
    points = np.array(
        [[0.09, 0.04, 0.01], [0.29, -0.13, -0.02], [0.62, 0.14, 0.03]],
        dtype=float,
    )
    time = 0.5
    expected = candidate.leading.at_points(points, time) + candidate.correction.at_points(
        points, time
    )
    np.testing.assert_allclose(candidate.at_points(points, time), expected, rtol=0.0, atol=0.0)

    path = candidate.save_json(tmp_path / "candidate.json")
    restored = KokunoReferenceCompositeCandidate.load_json(path)
    assert restored.sha256 == candidate.sha256
    np.testing.assert_allclose(restored.at_points(points, time), expected, rtol=0.0, atol=0.0)

    grid = candidate.grid([-0.1, 0.1], [-0.08, 0.08], [-0.02, 0.02], [0.5])
    assert grid.shape == (1, 2, 2, 2, 3)
    assert np.isfinite(grid).all()


def test_checkpoint_binds_all_four_agents_and_fails_closed():
    checkpoint = build_checkpoint()
    evidence = checkpoint["upstream_evidence"]
    assert {entry["pr"] for entry in evidence.values()} == {251, 252, 253, 254}
    assert (
        evidence["agent1_reference_continuation"]["candidate_sha256"]
        == evidence["agent4_independent_reference_audit"]["candidate_sha256"]
    )
    assert checkpoint["fixed_gates"]["held_out_normalized_full_momentum_residual"] == REGISTERED_PDE_THRESHOLD
    assert checkpoint["fixed_gates"]["divergence_max"] == REGISTERED_DIVERGENCE_THRESHOLD
    assert checkpoint["fixed_gates"]["changed"] is False

    findings = checkpoint["integration_findings"]
    assert findings["retained_oscillation_improves_inner_mean_rms"] is True
    assert findings["retained_oscillation_materially_repairs_post_core_seam"] is False
    assert findings["agent4_large_sampled_transition_and_post_core_divergence"] is True
    assert findings["agent4_large_sampled_momentum_obstruction"] is True
    assert findings["agent3_rank1_signed_covariance_realizability_ready"] is True
    assert findings["agent3_public_velocity_correction_materialized"] is False
    assert findings["reference_stage_eligible_for_global_promotion"] is False

    states = checkpoint["states"]
    assert states["reference_composite_velocity_export_ready"] is True
    assert states["leading_ready"] is False
    assert states["correction_ready"] is False
    assert states["velocity_export_ready"] is False
    assert states["formal_full_domain_pde_gate_assessed"] is False
    assert states["pde_validated"] is False


def test_candidate_truth_boundary_and_sha_are_fail_closed(tmp_path):
    candidate = KokunoReferenceCompositeCandidate()
    path = candidate.save_json(tmp_path / "candidate.json")
    payload = json.loads(path.read_text())
    payload["truth_boundary"]["pde_validated"] = True
    path.write_text(json.dumps(payload))
    with pytest.raises(ValueError, match="metadata changed"):
        KokunoReferenceCompositeCandidate.load_json(path)


def test_one_command_checkpoint_writes_replayable_bundle(tmp_path):
    checkpoint = write_checkpoint(tmp_path)
    assert (tmp_path / "reference_composite_candidate.json").is_file()
    assert (tmp_path / "integration_checkpoint.json").is_file()
    assert checkpoint["delivery_smoke"]["finite"] is True
    assert checkpoint["delivery_smoke"]["nontrivial"] is True
    assert checkpoint["delivery_smoke"]["save_load_same_sha"] is True
    assert checkpoint["delivery_smoke"]["used_for_pde_acceptance"] is False
    stored = json.loads((tmp_path / "integration_checkpoint.json").read_text())
    assert stored["candidate_sha256"] == checkpoint["candidate_sha256"]
    assert stored["checkpoint_sha256"] == checkpoint["checkpoint_sha256"]
