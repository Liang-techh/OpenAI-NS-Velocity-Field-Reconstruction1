from __future__ import annotations

import inspect

import pytest

from openai_ns_reconstruction.kokuno_current_i2_nonlinear_mean_attribution import (
    AGENT1_LEADING_HEAD,
    AGENT1_LEADING_PR,
    AGENT1_LEADING_SOURCE_BLOB,
    AGENT2_COMPOSITE_HEAD,
    AGENT2_COMPOSITE_PR,
    AGENT2_COMPOSITE_SOURCE_BLOB,
    AGENT2_DIFFERENTIAL_HEAD,
    AGENT2_DIFFERENTIAL_PR,
    AGENT2_DIFFERENTIAL_SOURCE_BLOB,
    FINAL_NORMALIZED_DIVERGENCE_GATE,
    FINAL_NORMALIZED_MOMENTUM_GATE,
    PARENT_AGENT3_HEAD,
    PARENT_AGENT3_PR,
    materialize_current_i2_nonlinear_mean_attribution,
    truth_boundary,
)


def test_public_surface_has_no_scientific_tuning_controls() -> None:
    assert tuple(inspect.signature(materialize_current_i2_nonlinear_mean_attribution).parameters) == (
        "backend",
        "radius",
        "z",
        "t",
    )
    truth = truth_boundary()
    assert truth["public_parameters"] == ("backend", "radius", "z", "t")
    assert truth["forbidden_scientific_controls_exposed"] is False
    assert truth["final_normalized_momentum_gate"] == FINAL_NORMALIZED_MOMENTUM_GATE == 1.0e-3
    assert truth["final_normalized_divergence_gate"] == FINAL_NORMALIZED_DIVERGENCE_GATE == 1.0e-5


def test_exact_current_i2_lineage_is_pinned() -> None:
    truth = truth_boundary()
    assert truth["parent_agent3_pr"] == PARENT_AGENT3_PR == 1063
    assert truth["parent_agent3_head"] == PARENT_AGENT3_HEAD == "e6ac8bed0ef775232e18e81e87f6a7572dcba8fb"
    assert truth["agent2_composite_pr"] == AGENT2_COMPOSITE_PR == 1071
    assert truth["agent2_composite_head"] == AGENT2_COMPOSITE_HEAD == "48d69e37e78e7f7f0e4e9936f28ff7974719288d"
    assert truth["agent2_composite_source_blob"] == AGENT2_COMPOSITE_SOURCE_BLOB == "be68248aca97173a15474324f63efff2c9ffce56"
    assert truth["agent2_differential_pr"] == AGENT2_DIFFERENTIAL_PR == 960
    assert truth["agent2_differential_head"] == AGENT2_DIFFERENTIAL_HEAD == "6d2fb1f701a34f783dca15a267ae2ce0734ba741"
    assert truth["agent2_differential_source_blob"] == AGENT2_DIFFERENTIAL_SOURCE_BLOB == "12df3bf6c949baeebf609b366f2973e7f515a157"
    assert truth["agent1_leading_pr"] == AGENT1_LEADING_PR == 1061
    assert truth["agent1_leading_head"] == AGENT1_LEADING_HEAD == "9c2abbc69ba1ef8c8a9e8d1699a466a75d0293c3"
    assert truth["agent1_leading_source_blob"] == AGENT1_LEADING_SOURCE_BLOB == "bfeff163d304c4cf1fb84eb8ff1cdfcf430870b4"


def test_truth_boundary_stays_scoped_to_i2_mean_attribution() -> None:
    truth = truth_boundary()
    assert truth["current_I2_leading_plus_oscillatory_identity_consumed"] is True
    assert truth["current_I2_mixed_nonlinear_mean_materialized"] is True
    assert truth["current_I2_quadratic_oscillatory_mean_materialized"] is True
    assert truth["current_I2_aggregate_nonlinear_mean_materialized"] is True
    assert truth["public_float64_I2_delta_verified"] is False
    assert truth["agent2_oscillatory_curl_or_jacobian_reimplemented"] is False

    for key in (
        "current_I3_overlay_consumed",
        "current_I4_overlay_consumed",
        "outer_global_velocity_consumed",
        "radial_inverse_performed_in_this_increment",
        "compact_radial_stress_materialized_in_this_increment",
        "radial_force_materialized_in_this_increment",
        "scoped_nonlinear_mean_authorized_as_complete_ns_defect",
        "scoped_nonlinear_mean_authorized_as_correction_target",
        "matched_global_pressure_materialized",
        "restricted_forcing_materialized",
        "complete_ns_defect",
        "current_real_ns_correction_velocity_materialized",
        "real_candidate_finite_correction_cycle_run",
        "heldout_normalized_ns_residual_assessed",
        "residual_reduction_claimed",
        "same_protocol_comparable_to_st006",
        "paper_exact",
        "openai_field_identified",
        "blowup_proof_claimed",
        "pde_validated",
    ):
        assert truth[key] is False


def test_materializer_rejects_untyped_backend_before_any_residual_work() -> None:
    with pytest.raises(TypeError, match="ExactCurrentI2NonlinearBackend"):
        materialize_current_i2_nonlinear_mean_attribution(object(), 0.3, 0.0, 0.5)
