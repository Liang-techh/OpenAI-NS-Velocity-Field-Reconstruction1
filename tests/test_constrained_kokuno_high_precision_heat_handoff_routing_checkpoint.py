from __future__ import annotations

import copy
import json

import pytest

from openai_ns_reconstruction.kokuno_high_precision_heat_handoff_routing_checkpoint import (
    FIXED_GATES,
    STATES,
    ST006_REFERENCE,
    build_checkpoint,
    load_checkpoint,
    validate_checkpoint,
    write_bundle,
)


@pytest.fixture(scope="module")
def checkpoint() -> dict:
    return build_checkpoint()


def test_high_precision_rebuild_removes_old_dominant_channel_blocker(checkpoint: dict) -> None:
    old = checkpoint["previous_agent4_rejection"]
    new = checkpoint["agent1_high_precision_handoff"]
    assert old["scientific_result"] == "REJECT_LEGACY_FLOAT64_MOMENT_TENSOR"
    assert old["legacy_max_I_sub_relative_error"] > 1.0e-10
    assert new["dominant_I_sub_guard_passed"] is True
    assert new["new_I_sub_relative_error"] < 1.0e-10
    assert new["new_I_sub_relative_error"] < new["legacy_I_sub_relative_error"] * 1.0e-3
    assert new["improvement_factor"] > 1000.0
    assert new["construction_max_relative_residual"] < 1.0e-40
    assert new["target_sign"] == [1, -1, 1]


def test_tiny_channels_remain_fail_closed_pending_high_precision_independent_audit(checkpoint: dict) -> None:
    handoff = checkpoint["agent1_high_precision_handoff"]
    assert handoff["tiny_Cp_S_independently_high_precision_certified"] is False
    assert handoff["independent_high_precision_moment_audit_completed"] is False
    assert handoff["continuous_source_moment_compensation_certified"] is False
    assert checkpoint["states"]["independent_high_precision_moment_audit_completed"] is False
    assert checkpoint["states"]["continuous_source_moment_compensation_certified"] is False
    assert checkpoint["states"]["heat_compensation_completed"] is False
    assert checkpoint["states"]["leading_ready"] is False
    assert checkpoint["routing"]["promote_repaired_I2_to_completed_heat_compensation_now"] is False


def test_sibling_routes_are_bound_without_laundering(checkpoint: dict) -> None:
    a2 = checkpoint["upstream"]["agent2"]
    a3 = checkpoint["upstream"]["agent3"]
    assert a2["homogeneous_pulse_spatial_sensitivity_ready"] is True
    assert a2["actual_positive_order_source_path_instantiated"] is False
    assert a2["public_xyz_t_velocity_correction_materialized"] is False
    assert a2["genuinely_independent_second_covariance_column_ready"] is False
    assert a3["spacetime_screen_times"] == [0.375, 0.5, 0.625]
    assert a3["spacetime_screen_z"] == [0.06, 0.08, 0.10]
    assert a3["actual_second_public_covariance_column_available"] is False
    assert a3["finite_cycle_rerun_allowed"] is False
    assert checkpoint["routing"]["rerun_one_or_duplicate_column_finite_cycle"] is False


def test_st006_and_formal_gates_remain_separate(checkpoint: dict) -> None:
    assert checkpoint["fixed_gates"] == FIXED_GATES
    assert checkpoint["st006_reference"] == ST006_REFERENCE
    rows = checkpoint["baseline_vs_kokuno"]
    assert rows[0]["momentum_sampled_max"] == ST006_REFERENCE["momentum_sampled_max"]
    assert rows[0]["volume_l2"] == ST006_REFERENCE["volume_l2"]
    assert rows[1]["directly_comparable_to_ST006"] is False
    assert checkpoint["states"] == STATES
    assert checkpoint["states"]["formal_full_domain_pde_gate_assessed"] is False
    assert checkpoint["states"]["pde_validated"] is False


def test_checkpoint_rejects_truth_laundering(checkpoint: dict) -> None:
    mutations = []

    item = copy.deepcopy(checkpoint)
    item["states"]["leading_ready"] = True
    mutations.append(item)

    item = copy.deepcopy(checkpoint)
    item["agent1_high_precision_handoff"]["tiny_Cp_S_independently_high_precision_certified"] = True
    mutations.append(item)

    item = copy.deepcopy(checkpoint)
    item["baseline_vs_kokuno"][1]["directly_comparable_to_ST006"] = True
    mutations.append(item)

    item = copy.deepcopy(checkpoint)
    item["routing"]["run_formal_full_domain_gate_now"] = True
    mutations.append(item)

    item = copy.deepcopy(checkpoint)
    item["truth_boundary"]["threshold_relaxed"] = True
    mutations.append(item)

    for mutated in mutations:
        with pytest.raises(ValueError, match="SHA mismatch"):
            validate_checkpoint(mutated)


def test_bundle_round_trip(tmp_path) -> None:
    checkpoint = write_bundle(tmp_path)
    path = tmp_path / "high_precision_heat_handoff_routing_checkpoint.json"
    assert load_checkpoint(path) == checkpoint
    raw = json.loads(path.read_text(encoding="utf-8"))
    assert raw["truth_boundary"]["ordinary_precision_I_sub_crosscheck_is_not_full_high_precision_certification"] is True
    assert raw["truth_boundary"]["tiny_Cp_S_channels_still_pending_independent_high_precision_audit"] is True
    assert raw["truth_boundary"]["free_residual_defined_forcing_used"] is False
    assert raw["truth_boundary"]["threshold_relaxed"] is False
    assert raw["routing"]["run_formal_full_domain_gate_now"] is False
