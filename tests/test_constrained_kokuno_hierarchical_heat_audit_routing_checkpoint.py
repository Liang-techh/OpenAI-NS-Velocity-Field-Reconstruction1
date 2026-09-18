from __future__ import annotations

import copy
import json

import pytest

from openai_ns_reconstruction.kokuno_hierarchical_heat_audit_routing_checkpoint import (
    FIXED_GATES,
    STATES,
    ST006_REFERENCE,
    bind_hierarchical_heat_audit,
    build_checkpoint,
    load_checkpoint,
    validate_checkpoint,
    write_bundle,
)


@pytest.fixture(scope="module")
def audit() -> dict:
    return bind_hierarchical_heat_audit()


@pytest.fixture(scope="module")
def checkpoint(audit: dict) -> dict:
    return build_checkpoint(audit)


def test_checkpoint_preserves_independent_scientific_rejection(checkpoint: dict) -> None:
    audit = checkpoint["agent4_independent_audit"]
    assert audit["scientific_result"] == "REJECT_CONTINUOUS_MOMENT_CLOSURE"
    assert audit["local_audit_completed"] is False
    assert audit["summary"]["independent_continuous_moment_closure"] is False
    assert audit["summary"]["frozen_discrete_map_only"] is True
    assert audit["local_guards"]["dominant_I_sub_channel_relative_error_le_1e-10"] is False
    assert checkpoint["fixed_gates"] == FIXED_GATES
    assert checkpoint["st006_reference"] == ST006_REFERENCE
    assert checkpoint["states"] == STATES
    assert checkpoint["states"]["pde_validated"] is False


def test_discrete_decimal_closure_is_not_laundered_into_continuous_closure(checkpoint: dict) -> None:
    construction = checkpoint["agent1_discrete_hierarchical_receipt"]
    gap = checkpoint["construction_vs_independent_precision_gap"]
    assert construction["precision_digits"] >= 480
    assert construction["target_sign"] == [1, -1, 1]
    assert construction["max_relative_residual"] < 1.0e-40
    assert construction["continuous_source_moment_compensation_certified"] is False
    assert gap["independent_I_sub_relative_error"] > 1.0e-10
    assert gap["guard_excess_orders"] > 4.0
    assert gap["minimum_Cp_log10_relative_error"] > 50.0
    assert gap["minimum_S_log10_relative_error"] > 50.0
    assert checkpoint["routing"]["promote_current_hierarchical_repair_to_heat_compensation"] is False


def test_sibling_routes_remain_fail_closed(checkpoint: dict) -> None:
    a2 = checkpoint["upstream"]["agent2"]
    a3 = checkpoint["upstream"]["agent3"]
    assert a2["source_homogeneous_primary_pulse_ready"] is True
    assert a2["caller_supplied_source_path_geometry_still_required"] is True
    assert a2["public_xyz_t_velocity_correction_materialized"] is False
    assert a2["genuinely_independent_second_covariance_column_ready"] is False
    assert a3["node_time_requirements_total"] == 75
    assert a3["rank_two_requirements_passed"] == 0
    assert a3["finite_cycle_rerun_allowed"] is False
    assert checkpoint["routing"]["rerun_one_or_duplicate_column_finite_cycle"] is False


def test_checkpoint_rejects_truth_laundering(checkpoint: dict) -> None:
    mutated = copy.deepcopy(checkpoint)
    mutated["states"]["pde_validated"] = True
    with pytest.raises(ValueError, match="SHA mismatch"):
        validate_checkpoint(mutated)

    mutated = copy.deepcopy(checkpoint)
    mutated["baseline_vs_kokuno"][1]["directly_comparable_to_ST006"] = True
    with pytest.raises(ValueError, match="SHA mismatch"):
        validate_checkpoint(mutated)

    mutated = copy.deepcopy(checkpoint)
    mutated["agent4_independent_audit"]["scientific_result"] = "PASS"
    with pytest.raises(ValueError, match="SHA mismatch"):
        validate_checkpoint(mutated)


def test_bundle_round_trip(tmp_path) -> None:
    checkpoint = write_bundle(tmp_path)
    path = tmp_path / "hierarchical_heat_audit_routing_checkpoint.json"
    assert load_checkpoint(path) == checkpoint
    raw = json.loads(path.read_text(encoding="utf-8"))
    assert raw["truth_boundary"]["scientific_rejection_preserved_in_green_CI"] is True
    assert raw["truth_boundary"]["free_residual_defined_forcing_used"] is False
    assert raw["truth_boundary"]["threshold_relaxed"] is False
    assert raw["routing"]["run_formal_full_domain_gate_now"] is False
