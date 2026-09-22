from __future__ import annotations

import inspect

import pytest

from openai_ns_reconstruction.kokuno_current_i4_nonlinear_mean_attribution import (
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
    SOURCE_MEAN_DEFECT_FORMULAS,
    SOURCE_MEAN_INTERVAL_FORMULA,
    materialize_current_i4_nonlinear_mean_attribution,
    truth_boundary,
)


def test_public_surface_has_no_scientific_tuning_controls() -> None:
    assert tuple(inspect.signature(materialize_current_i4_nonlinear_mean_attribution).parameters) == (
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


def test_exact_current_i4_lineage_is_pinned() -> None:
    truth = truth_boundary()
    assert truth["parent_agent3_pr"] == PARENT_AGENT3_PR == 1090
    assert truth["parent_agent3_head"] == PARENT_AGENT3_HEAD == "5e8512392df683375fabd132f35577e1b94773ff"
    assert truth["agent2_composite_pr"] == AGENT2_COMPOSITE_PR == 1080
    assert truth["agent2_composite_head"] == AGENT2_COMPOSITE_HEAD == "c40d8ddecd2971544a6e07dab093436b423cf326"
    assert truth["agent2_composite_source_blob"] == AGENT2_COMPOSITE_SOURCE_BLOB == "2a0a5aa5966b02da856bdcf51940f3c186042802"
    assert truth["agent2_differential_pr"] == AGENT2_DIFFERENTIAL_PR == 960
    assert truth["agent2_differential_head"] == AGENT2_DIFFERENTIAL_HEAD == "6d2fb1f701a34f783dca15a267ae2ce0734ba741"
    assert truth["agent2_differential_source_blob"] == AGENT2_DIFFERENTIAL_SOURCE_BLOB == "12df3bf6c949baeebf609b366f2973e7f515a157"
    assert truth["agent1_leading_pr"] == AGENT1_LEADING_PR == 1079
    assert truth["agent1_leading_head"] == AGENT1_LEADING_HEAD == "b06742ca6e189499192ede3cce40f62cdc1e35ca"
    assert truth["agent1_leading_source_blob"] == AGENT1_LEADING_SOURCE_BLOB == "6f04ce0a856b44430402576dad88438da90d1ebb"


def test_truth_boundary_records_source_i4_role_without_promoting_correction() -> None:
    truth = truth_boundary()
    assert truth["source_mean_interval_formula"] == SOURCE_MEAN_INTERVAL_FORMULA == "RF40a: I_mean = I4"
    assert truth["source_mean_defect_formulas"] == SOURCE_MEAN_DEFECT_FORMULAS == "RF30--RF49"
    assert truth["current_I4_leading_plus_oscillatory_identity_consumed"] is True
    assert truth["current_I4_mixed_nonlinear_mean_materialized"] is True
    assert truth["current_I4_quadratic_oscillatory_mean_materialized"] is True
    assert truth["current_I4_aggregate_nonlinear_mean_materialized"] is True
    assert truth["source_I4_reserved_mean_correction_interval_recorded"] is True
    assert truth["agent2_oscillatory_curl_or_jacobian_reimplemented"] is False

    for key in (
        "source_I3_positive_order_correction_materialized",
        "source_I4_five_row_mean_correction_materialized",
        "current_I4_correction_velocity_materialized",
        "post_I4_or_pulse_leading_identity_consumed",
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


def test_materializer_rejects_untyped_backend_before_any_defect_work() -> None:
    with pytest.raises(TypeError, match="ExactCurrentI4NonlinearBackend"):
        materialize_current_i4_nonlinear_mean_attribution(object(), 0.3, 0.0, 0.5)
