from __future__ import annotations

import copy

import pytest

from openai_ns_reconstruction.kokuno_segmented_frontier_routing_checkpoint import (
    FIXED_GATES,
    ST006_REFERENCE,
    STATES,
    _sha,
    build_checkpoint,
    validate_checkpoint,
)


def _audit_fixture() -> dict:
    finest = {
        "step": 0.001,
        "sample_count": 4,
        "max_vector_residual": 0.12,
        "sample_l2_vector_residual": 0.08,
        "component_rms": [0.01, 0.02, 0.07],
        "divergence_max_abs": 2.0e-6,
        "divergence_sample_l2": 8.0e-7,
    }
    middle = {**finest, "step": 0.002, "max_vector_residual": 0.121}
    coarse = {**finest, "step": 0.004, "max_vector_residual": 0.124}
    region = {
        "resolution_ladder": [coarse, middle, finest],
        "fine_to_previous_relative_change": {
            "sample_l2_vector_residual": 0.0,
            "max_vector_residual": 0.008264462809917363,
        },
        "finest_local_residual_gate_met": False,
        "finest_local_divergence_gate_met": True,
    }
    return {
        "task_id": "KOKUNO-A4-SEGMENTED-PUBLIC-CONTRACT-AUDIT-016",
        "base": {
            "pr": 345,
            "head": "57df8375f134277ebea509847f6b319aa84acfd7",
        },
        "seed": 9173071,
        "fixed_project_contract": {
            "nu": 0.01,
            "normalized_momentum_max_gate": 1.0e-3,
            "normalized_momentum_l2_gate": 1.0e-3,
            "divergence_max_gate": 1.0e-5,
            "divergence_l2_gate": 1.0e-5,
            "thresholds_changed": False,
        },
        "independent_operator": {
            "reference_velocity_input": "public only",
            "reference_pressure_input": "public only",
            "forcing": "zero local diagnostic",
            "spatial_and_time_derivatives": "independent centered Cartesian FD2",
            "steps": [0.004, 0.002, 0.001],
            "training_or_construction_derivatives_read": False,
            "training_loss_read": False,
        },
        "reference_stage": {"off_grid": copy.deepcopy(region), "axis_near": copy.deepcopy(region)},
        "mutation_calibration": {
            "declared_mutation": 0.02,
            "mean_divergence_shift": 0.02,
            "max_divergence_shift_error": 1.0e-12,
            "mean_pressure_gradient_x_shift": 0.02,
            "max_pressure_gradient_shift_error": 1.0e-12,
        },
        "i2_public_contract": {
            "sample_count": 4,
            "public_float64_equals_uncorrected_base_all_probes": True,
            "public_float64_equals_uncorrected_base_fraction": 1.0,
            "router_vs_uncorrected_base_max_abs": 0.0,
            "public_vs_independent_source_base_max_relative": 1.0e-14,
            "public_speed_min": 1.0e-33,
            "public_speed_max": 2.0e-33,
            "public_float64_heat_repair_visible": False,
            "module_level_after_heat_repair_attribution_available": False,
        },
        "preregistered_local_guards": {
            "public_vs_independent_I2_source_base_relative_le_5e-10": True,
            "mutation_shift_error_le_1e-8": True,
            "all_local_implementation_guards_passed": True,
        },
        "truth_boundary": {
            "formal_full_domain_pde_gate_assessed": False,
            "pde_validated": False,
        },
    }


def test_checkpoint_binds_segmented_router_without_global_promotion() -> None:
    checkpoint = build_checkpoint(_audit_fixture())
    assert validate_checkpoint(checkpoint) == checkpoint
    assert checkpoint["fixed_gates"] == FIXED_GATES
    assert checkpoint["st006_reference"] == ST006_REFERENCE
    assert checkpoint["states"]["segmented_leading_velocity_router_ready"] is True
    assert checkpoint["states"]["leading_ready"] is False
    assert checkpoint["states"]["correction_ready"] is False
    assert checkpoint["states"]["velocity_export_ready"] is False
    assert checkpoint["states"]["formal_full_domain_pde_gate_assessed"] is False
    assert checkpoint["states"]["pde_validated"] is False
    assert checkpoint["states"]["i2_float64_public_repair_observable_at_agent4_probes"] is False
    assert checkpoint["segmented_leading_contract"]["coverage"]["global_coverage_complete"] is False
    assert checkpoint["segmented_leading_contract"]["coverage"]["unreconstructed_log_X_gaps"]
    assert checkpoint["baseline_vs_kokuno"][1]["directly_comparable_to_ST006"] is False


def test_checkpoint_keeps_agent2_and_agent3_sibling_boundaries() -> None:
    checkpoint = build_checkpoint(_audit_fixture())
    assert checkpoint["upstream"]["agent2"]["consumed_in_executable_ancestry"] is False
    assert checkpoint["upstream"]["agent2"]["public_source_oscillatory_xyz_t_velocity_ready"] is False
    assert checkpoint["upstream"]["agent3"]["consumed_in_executable_ancestry"] is False
    assert checkpoint["upstream"]["agent3"]["radial_target_stable_for_future_column_screen"] is True
    assert checkpoint["upstream"]["agent3"]["finite_correction_cycle_rerun_allowed"] is False
    assert checkpoint["upstream"]["agent3"]["duplicate_control_rank2_required_nodes"] == 0
    assert checkpoint["upstream"]["agent3"]["duplicate_control_required_nodes"] == 99


@pytest.mark.parametrize(
    ("path", "value"),
    [
        (("states", "pde_validated"), True),
        (("states", "leading_ready"), True),
        (("states", "velocity_export_ready"), True),
        (("truth_boundary", "threshold_relaxed"), True),
        (("truth_boundary", "free_residual_defined_forcing_used"), True),
        (("segmented_leading_contract", "coverage", "global_coverage_complete"), True),
    ],
)
def test_checkpoint_rejects_truth_boundary_mutations(path: tuple[str, ...], value: object) -> None:
    checkpoint = build_checkpoint(_audit_fixture())
    mutated = copy.deepcopy(checkpoint)
    target = mutated
    for key in path[:-1]:
        target = target[key]
    target[path[-1]] = value
    unsigned = dict(mutated)
    unsigned.pop("checkpoint_sha256")
    mutated["checkpoint_sha256"] = _sha(unsigned)
    with pytest.raises(ValueError):
        validate_checkpoint(mutated)


def test_checkpoint_rejects_threshold_mutation_even_with_rehashed_payload() -> None:
    checkpoint = build_checkpoint(_audit_fixture())
    mutated = copy.deepcopy(checkpoint)
    mutated["fixed_gates"]["held_out_normalized_full_momentum_max"] = 2.0e-3
    unsigned = dict(mutated)
    unsigned.pop("checkpoint_sha256")
    mutated["checkpoint_sha256"] = _sha(unsigned)
    with pytest.raises(ValueError, match="fixed gates changed"):
        validate_checkpoint(mutated)


def test_checkpoint_rejects_agent4_local_guard_failure() -> None:
    audit = _audit_fixture()
    audit["preregistered_local_guards"]["all_local_implementation_guards_passed"] = False
    with pytest.raises(ValueError, match="local implementation guards"):
        build_checkpoint(audit)


def test_static_state_contract_remains_fail_closed() -> None:
    assert STATES["leading_ready"] is False
    assert STATES["public_source_oscillatory_xyz_t_velocity_ready"] is False
    assert STATES["genuinely_independent_second_covariance_column_ready"] is False
    assert STATES["correction_ready"] is False
    assert STATES["complete_kokuno_composite_velocity_ready"] is False
    assert STATES["velocity_export_ready"] is False
    assert STATES["formal_full_domain_pde_gate_assessed"] is False
    assert STATES["pde_validated"] is False
