from __future__ import annotations

import copy
import json

import pytest

from openai_ns_reconstruction.kokuno_signed_log_heat_routing_checkpoint import (
    FIXED_GATES,
    STATES,
    ST006_REFERENCE,
    bind_independent_heat_audit,
    build_checkpoint,
    load_checkpoint,
    validate_checkpoint,
    write_bundle,
)


@pytest.fixture(scope="module")
def audit() -> dict:
    return bind_independent_heat_audit()


@pytest.fixture(scope="module")
def checkpoint(audit: dict) -> dict:
    return build_checkpoint(audit)


def test_checkpoint_binds_independent_heat_preflight_fail_closed(checkpoint: dict) -> None:
    assert checkpoint["independent_signed_log_heat_audit"]["structural_preflight_passed"] is True
    assert all(checkpoint["independent_signed_log_heat_audit"]["guards"].values())
    assert checkpoint["fixed_gates"] == FIXED_GATES
    assert checkpoint["st006_reference"] == ST006_REFERENCE
    assert checkpoint["baseline_vs_kokuno"][1]["directly_comparable_to_ST006"] is False
    assert checkpoint["states"] == STATES
    assert checkpoint["states"]["pde_validated"] is False
    assert checkpoint["states"]["formal_full_domain_pde_gate_assessed"] is False


def test_default_precision_frontier_preserves_nonzero_underflow_channels(checkpoint: dict) -> None:
    frontier = checkpoint["default_eta_0p2_precision_frontier"]
    assert frontier["sign"] == [1, -1, 1]
    assert frontier["underflow_channels"] == ["C_p", "S"]
    assert frontier["channel_span_decades"] > 400.0
    assert frontier["fully_float64_materializable"] is False
    assert frontier["existing_float64_three_bump_inverse_faithful_for_default"] is False


def test_sibling_lane_routing_stays_blocked_until_physical_inputs_exist(checkpoint: dict) -> None:
    a2 = checkpoint["upstream"]["agent2"]
    a3 = checkpoint["upstream"]["agent3"]
    assert a2["source_real_conjugate_pair_kernel_ready"] is True
    assert a2["actual_source_pulse_background_instantiated"] is False
    assert a2["public_xyz_t_velocity_correction_materialized"] is False
    assert a2["genuinely_independent_second_covariance_column_ready"] is False
    assert a3["velocity_column_covariance_adapter_executable"] is True
    assert a3["nodes_requiring_second_direction"] == 25
    assert a3["rank_two_required_nodes"] == 0
    assert a3["finite_cycle_rerun_allowed"] is False
    assert checkpoint["routing"]["rerun_duplicate_or_one_column_correction_cycle"] is False


def test_checkpoint_rejects_truth_laundering(checkpoint: dict) -> None:
    mutated = copy.deepcopy(checkpoint)
    mutated["states"]["pde_validated"] = True
    with pytest.raises(ValueError, match="SHA mismatch"):
        validate_checkpoint(mutated)

    mutated = copy.deepcopy(checkpoint)
    mutated["baseline_vs_kokuno"][1]["directly_comparable_to_ST006"] = True
    with pytest.raises(ValueError, match="SHA mismatch"):
        validate_checkpoint(mutated)


def test_bundle_round_trip(tmp_path) -> None:
    checkpoint = write_bundle(tmp_path)
    path = tmp_path / "signed_log_heat_routing_checkpoint.json"
    assert load_checkpoint(path) == checkpoint
    raw = json.loads(path.read_text(encoding="utf-8"))
    assert raw["typed_heat_target_component"]["velocity_api_ready"] is False
    assert raw["truth_boundary"]["free_residual_defined_forcing_used"] is False
    assert raw["truth_boundary"]["threshold_relaxed"] is False
