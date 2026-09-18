import copy

import pytest

from openai_ns_reconstruction.kokuno_rational_rectangle_routing_checkpoint import (
    build_checkpoint,
    checkpoint_sha256,
    validate_checkpoint,
)


def test_checkpoint_is_deterministic_and_fail_closed():
    first = build_checkpoint()
    second = build_checkpoint()
    assert first == second
    validate_checkpoint(first)
    assert first["checkpoint_sha256"] == checkpoint_sha256(first)
    assert first["states"]["source_rational_rectangle_separation_contract_executable"]
    assert first["states"]["autonomous_rational_rectangle_witness_ready"]
    assert not first["states"]["source_rectangle_parameters_recovered"]
    assert not first["states"]["oscillatory_ready"]
    assert not first["states"]["correction_ready"]
    assert not first["states"]["candidate_artifact_instantiated"]
    assert not first["states"]["pde_validated"]


def test_agent2_receipt_preserves_autonomous_truth_boundary():
    payload = build_checkpoint()
    a2 = payload["upstream"]["agent2"]
    assert a2["pr"] == 439
    assert a2["head"] == "e5bc36691d16be09cab602885fc49cb99f3d2f71"
    assert a2["dedicated_status"] == "success"
    assert a2["standard_status"] == "success"
    assert a2["source_rational_rectangle_separation_contract_executable"]
    assert a2["autonomous_rational_witness_executable"]
    assert not a2["source_rectangle_centers_recovered"]
    assert not a2["source_rectangle_radius_r0_recovered"]
    assert not a2["source_color_count_recovered"]
    assert not a2["source_center_denominator_recovered"]
    assert not a2["actual_positive_order_background_bound"]
    assert not a2["actual_auxiliary_torus_mode_family_bound"]
    assert not a2["public_q_scaled_by_beta_total_xyz_t_velocity_ready"]


def test_pending_agent3_is_not_promoted_and_old_rank_blocker_remains():
    payload = build_checkpoint()
    stable = payload["upstream"]["agent3_stable"]
    pending = payload["upstream"]["agent3_latest_pending"]
    assert stable["required_rank_two_nodes"] == 25
    assert stable["rank_two_nodes"] == 0
    assert not stable["upstream_bounded_inverse_passed"]
    assert not stable["finite_correction_cycle_rerun_allowed"]
    assert pending["pr"] == 440
    assert pending["dedicated_status_at_routing_audit"] == "in_progress"
    assert pending["standard_status_at_routing_audit"] == "in_progress"
    assert not pending["promoted_to_stable_receipt"]


def test_formal_gates_and_st006_baseline_are_unchanged():
    payload = build_checkpoint()
    gates = payload["formal_gates"]
    baseline = payload["baseline_vs_kokuno"]["st006"]
    assert gates["held_out_normalized_momentum_max"] == 1.0e-3
    assert gates["held_out_normalized_momentum_l2"] == 1.0e-3
    assert gates["held_out_divergence_max"] == 1.0e-5
    assert gates["held_out_divergence_l2"] == 1.0e-5
    assert not gates["changed_this_round"]
    assert baseline["momentum_sampled_max"] == 0.1082289305112118
    assert baseline["momentum_volume_l2"] == 0.10758432876230622
    assert payload["baseline_vs_kokuno"]["kokuno_current_comparable_full_domain_receipt"] is None


def test_validator_rejects_truth_or_gate_promotion():
    payload = build_checkpoint()
    bad = copy.deepcopy(payload)
    bad["states"]["pde_validated"] = True
    bad["checkpoint_sha256"] = checkpoint_sha256(bad)
    with pytest.raises(ValueError, match="fail-closed"):
        validate_checkpoint(bad)

    bad = copy.deepcopy(payload)
    bad["upstream"]["agent2"]["source_rectangle_centers_recovered"] = True
    bad["checkpoint_sha256"] = checkpoint_sha256(bad)
    with pytest.raises(ValueError, match="hidden source rectangle"):
        validate_checkpoint(bad)

    bad = copy.deepcopy(payload)
    bad["formal_gates"]["held_out_normalized_momentum_max"] = 2.0e-3
    bad["checkpoint_sha256"] = checkpoint_sha256(bad)
    with pytest.raises(ValueError, match="momentum gate changed"):
        validate_checkpoint(bad)
