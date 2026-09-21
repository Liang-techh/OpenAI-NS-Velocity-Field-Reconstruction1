from __future__ import annotations

import inspect

from openai_ns_reconstruction.kokuno_complete_ns_defect_cycle_admission import (
    PARENT_AGENT3_HEAD,
    PARENT_AGENT3_PR,
    PARENT_AGENT3_SOURCE_BLOB,
    REQUIRED_EVIDENCE,
    materialize_current_complete_ns_defect_cycle_admission,
    truth_boundary,
)


def test_public_admission_api_accepts_no_scientific_payload() -> None:
    signature = inspect.signature(materialize_current_complete_ns_defect_cycle_admission)
    assert tuple(signature.parameters) == ()
    boundary = truth_boundary()
    assert boundary["forbidden_public_parameters_absent"] is True
    assert boundary["caller_supplied_surrogate_defect_allowed"] is False
    assert boundary["caller_supplied_residual_arrays_allowed"] is False
    assert boundary["caller_supplied_forcing_allowed"] is False
    assert boundary["residual_defined_forcing_allowed"] is False


def test_current_x4_lineage_fails_closed_before_finite_cycle() -> None:
    witness = materialize_current_complete_ns_defect_cycle_admission()
    receipt = witness.to_receipt()

    assert receipt["parent_agent3_pr"] == PARENT_AGENT3_PR == 1045
    assert receipt["parent_agent3_head"] == PARENT_AGENT3_HEAD
    assert receipt["parent_agent3_source_blob"] == PARENT_AGENT3_SOURCE_BLOB
    assert receipt["finite_cycle_admitted"] is False
    assert tuple(receipt["missing_required_evidence"]) == REQUIRED_EVIDENCE
    assert receipt["fixed_scientific_gates"] == {
        "normalized_momentum": 1.0e-3,
        "normalized_divergence": 1.0e-5,
    }

    parent = receipt["parent_truth_boundary"]
    for key in (
        "outer_global_leading_velocity_materialized",
        "matched_global_pressure_materialized",
        "restricted_forcing_included",
        "complete_ns_defect",
        "current_real_ns_correction_velocity_materialized",
        "real_candidate_finite_correction_cycle_run",
        "heldout_normalized_ns_residual_assessed",
    ):
        assert parent[key] is False, key

    boundary = receipt["truth_boundary"]
    for key in (
        "complete_ns_defect",
        "complete_defect_identity_bound_to_global_candidate",
        "matched_global_pressure_materialized",
        "preregistered_restricted_forcing_materialized",
        "forcing_proven_not_residual_defined",
        "held_in_held_out_sample_sets_disjoint",
        "independent_second_covariance_column_passed",
        "second_column_bounded_inverse_preflight_passed",
        "current_real_ns_correction_velocity_materialized",
        "correction_cycle_gain_guard_accepted",
        "finite_correction_cycle_admitted",
        "real_candidate_finite_correction_cycle_run",
        "heldout_normalized_ns_residual_assessed",
        "residual_reduction_claimed",
        "same_protocol_comparable_to_st006",
        "pde_validated",
        "blowup_proved",
    ):
        assert boundary[key] is False, key

    assert boundary["same_samples_allowed_for_held_in_and_held_out"] is False
    assert boundary["scoped_mean_stress_force_counts_as_complete_defect"] is False
    assert boundary["historical_gain_guard_counts_as_current_acceptance"] is False
    assert boundary["historical_second_column_preflight_counts_as_current_evidence"] is False
    assert boundary["final_normalized_momentum_gate"] == 1.0e-3
    assert boundary["final_normalized_divergence_gate"] == 1.0e-5
