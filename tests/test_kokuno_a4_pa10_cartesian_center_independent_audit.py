from __future__ import annotations

import inspect

import pytest

import openai_ns_reconstruction.kokuno_a4_pa10_cartesian_center_independent_audit as audit_module
from openai_ns_reconstruction.kokuno_a4_pa10_cartesian_center_independent_audit import (
    DIVERGENCE_L2_GATE,
    DIVERGENCE_MAX_GATE,
    FD_STEPS,
    Q_RELATIVE_MAX_GATE,
    VELOCITY_RELATIVE_MAX_GATE,
    run_independent_audit,
)


@pytest.fixture(scope="module")
def audit():
    return run_independent_audit().payload


def test_independent_cartesian_center_audit_passes_frozen_guards(audit):
    assert audit["passed"] is True
    assert audit["failed_guards"] == []
    measurements = audit["measurements"]
    assert measurements["q_relative_max"] <= Q_RELATIVE_MAX_GATE
    assert measurements["X_relative_max"] <= Q_RELATIVE_MAX_GATE
    assert measurements["eta_relative_max"] <= Q_RELATIVE_MAX_GATE
    assert measurements["velocity_relative_max"] <= VELOCITY_RELATIVE_MAX_GATE
    assert measurements["velocity_rms"] > 0.0
    assert measurements["velocity_max_abs"] > 0.0
    assert measurements["axis_transverse_max_abs"] == 0.0

    ladder = measurements["divergence_by_step"]
    assert list(float(key) for key in ladder) == list(FD_STEPS)
    finest = ladder[f"{FD_STEPS[-1]:.6f}"]
    assert finest["max_abs"] <= DIVERGENCE_MAX_GATE
    assert finest["l2"] <= DIVERGENCE_L2_GATE


def test_validator_does_not_reuse_candidate_coordinate_helpers():
    source = inspect.getsource(audit_module)
    forbidden = (
        ".solve_q(",
        ".similarity_coordinates(",
        ".cartesian_from_similarity(",
        ".report(",
    )
    for token in forbidden:
        assert token not in source
    assert "safeguarded Newton" in source
    assert "_fd4_divergence" in source


def test_mutation_controls_are_pre_registered_and_detected(audit):
    mutations = audit["mutation_controls"]
    assert (
        mutations["public_velocity_times_0p999_relative_max"]
        > VELOCITY_RELATIVE_MAX_GATE
    )
    assert (
        mutations["independent_q_times_0p999_coordinate_error"]
        > Q_RELATIVE_MAX_GATE
    )
    assert audit["guards"]["scaled_velocity_mutation_detected"] is True
    assert audit["guards"]["q_coordinate_mutation_detected"] is True


def test_three_resolution_divergence_protocol_is_frozen(audit):
    protocol = audit["protocol"]
    assert protocol["fd4_spatial_steps"] == [0.004, 0.002, 0.001]
    assert protocol["divergence_max_gate"] == 1e-5
    assert protocol["divergence_l2_gate"] == 1e-5
    assert protocol["final_project_momentum_gate"] == 1e-3
    assert protocol["final_project_divergence_gate"] == 1e-5
    assert protocol["free_residual_defined_forcing_forbidden"] is True
    assert audit["guards"]["fd4_three_level_resolution_stable"] is True


def test_truth_boundary_remains_inner_center_only(audit):
    truth = audit["truth_boundary"]
    assert truth["inner_cartesian_center_velocity_independently_audited"] is True
    assert truth["agent4_independent_coordinate_inverse_used"] is True
    assert truth["agent4_public_velocity_only_fd4_divergence_used"] is True
    assert truth["source_center_is_final_corrected_fixed_point"] is False
    assert truth["fixed_point_correction_materialized"] is False
    assert truth["global_cartesian_spacetime_leading_velocity_materialized"] is False
    assert truth["outer_join_localization_materialized"] is False
    assert truth["matched_global_pressure_materialized"] is False
    assert truth["complete_restricted_forcing_materialized"] is False
    assert truth["complete_kokuno_composite_velocity"] is False
    assert truth["leading_only_ns_residual_assessed"] is False
    assert truth["leading_plus_oscillatory_ns_residual_assessed"] is False
    assert truth["after_correction_ns_residual_assessed"] is False
    assert truth["heldout_ns_momentum_residual_assessed"] is False
    assert truth["formal_full_domain_pde_gate_assessed"] is False
    assert truth["pde_validated"] is False
