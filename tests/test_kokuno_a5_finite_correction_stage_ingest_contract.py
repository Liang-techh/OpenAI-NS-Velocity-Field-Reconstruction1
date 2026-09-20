from __future__ import annotations

import copy

import pytest

from openai_ns_reconstruction.kokuno_a5_finite_correction_stage_ingest_contract import (
    AGENT3_DEDICATED_RUN,
    AGENT3_HEAD,
    AGENT3_SOURCE_BLOB,
    AGENT3_TESTS_RUN,
    AGENT3_WORKFLOW_BLOB,
    FINAL_NORMALIZED_DIVERGENCE_GATE,
    FINAL_NORMALIZED_MOMENTUM_GATE,
    PARENT_A5_HEAD,
    deterministic_finite_correction_stage_ingest_contract,
    validate_finite_correction_stage_ingest_contract,
)
from openai_ns_reconstruction.kokuno_a5_strict_inner_vorticity_artifact_ingest_contract import (
    deterministic_strict_inner_vorticity_artifact_ingest_contract,
)

EXACT_HEAD = "0" * 40


def _payload():
    return deterministic_finite_correction_stage_ingest_contract(exact_head=EXACT_HEAD)


def test_agent3_stage_is_registered_but_not_a_real_cycle() -> None:
    payload = _payload()
    validate_finite_correction_stage_ingest_contract(payload)

    assert payload["parent_a5"]["pr"] == 884
    assert payload["parent_a5"]["head"] == PARENT_A5_HEAD
    a3 = payload["agent3_binding"]
    assert a3["pr"] == 888
    assert a3["head"] == AGENT3_HEAD
    assert a3["source_blob_sha"] == AGENT3_SOURCE_BLOB
    assert a3["workflow_blob_sha"] == AGENT3_WORKFLOW_BLOB
    assert a3["dedicated_run"] == AGENT3_DEDICATED_RUN == 35530175857
    assert a3["tests_run"] == AGENT3_TESTS_RUN == 35530175820
    assert a3["module"] == "openai_ns_reconstruction.kokuno_finite_correction_stage"
    assert a3["public_api"] == (
        "run_finite_correction_stage(backend,candidate,correction,held_in,held_out)"
    )
    assert a3["parent_pr"] == 882
    assert a3["parent_head"] == "29b385398d0c7636ed8d4049821a534d80728b33"

    evidence = payload["evidence"]
    assert evidence["parent_vorticity_artifact_ingest_admitted"] is False
    assert evidence["agent3_exact_head_ci_conclusion"] is None
    assert evidence["agent3_typed_stage_contract_present"] is True
    assert evidence["agent3_mechanics_receipt_is_candidate_residual_evidence"] is False
    assert evidence["agent3_real_candidate_finite_correction_cycle_run"] is False
    assert evidence["agent4_dedicated_finite_correction_stage_audit_present"] is False
    assert evidence["agent4_finite_correction_stage_audit_conclusion"] is None
    assert evidence["finite_correction_stage_contract_ingest_admitted"] is False
    assert evidence["real_correction_scientifically_admitted"] is False

    handoff = payload["finite_correction_stage_handoff"]
    assert handoff["registered"] is True
    assert handoff["status"] == "registered_unresolved"
    assert handoff["construction_authority"] == "Agent3#888"
    assert handoff["candidate_and_correction_provenance_typed"] is True
    assert handoff["disjoint_held_in_held_out_enforced"] is True
    assert handoff["held_out_fit_leakage_rejected"] is True
    assert handoff["same_residual_protocol_before_after_enforced"] is True
    assert handoff["complete_ns_residual_required"] is True
    assert handoff["matched_pressure_required"] is True
    assert handoff["restricted_forcing_required"] is True
    assert handoff["restricted_forcing_preregistered_required"] is True
    assert handoff["residual_as_forcing_shortcut_forbidden"] is True
    assert handoff["correction_divergence_gate"] == 1e-5
    assert handoff["mechanics_only_stage_receipt_is_scientific_candidate_evidence"] is False
    assert handoff["analytic_kokuno_gain_available"] is False
    assert handoff["dedicated_agent4_stage_audit_available"] is False
    assert handoff["usable_for_real_cycle_now"] is False
    assert handoff["usable_as_st006_comparison_now"] is False
    assert handoff["usable_for_final_independent_pde_validation_now"] is False

    assert payload["ingest_status"]["typed_finite_correction_stage_contract_registered"] is True
    assert payload["ingest_status"]["finite_correction_stage_contract_ingest_admitted"] is False
    assert payload["ingest_status"]["real_finite_correction_cycle_admitted"] is False
    assert payload["artifact_status"]["real_agent3_correction_velocity_materialized"] is False
    assert payload["artifact_status"]["real_finite_correction_cycle_receipt_materialized"] is False
    assert payload["stage_state"] == {
        "leading_ready": False,
        "oscillatory_ready": True,
        "correction_ready": False,
        "velocity_export_ready": False,
        "pde_validated": False,
    }


def test_protocol_freezes_complete_defect_and_antileakage_requirements() -> None:
    protocol = _payload()["agent3_binding"]["protocol"]
    assert protocol["candidate_admission_requires"] == [
        "complete_ns_defect",
        "corrected_global_leading_join_complete",
        "matched_pressure_gradient_included",
        "restricted_forcing_included",
        "restricted_forcing_preregistered",
        "residual_as_forcing_shortcut_used=false",
    ]
    assert "fit_sample_ids_subset_of_held_in" in protocol["correction_provenance_requires"]
    assert "fit_sample_ids_disjoint_from_held_out" in protocol["correction_provenance_requires"]
    assert protocol["held_in_held_out_partition_ids_must_differ"] is True
    assert protocol["held_in_held_out_sample_ids_must_be_disjoint"] is True
    assert protocol["backend_recomputes_complete_ns_residual_before_and_after"] is True
    assert protocol["same_residual_protocol_sha_required_before_and_after"] is True
    assert protocol["incomplete_ns_residual_rejected"] is True
    assert protocol["reported_momentum_metrics"] == [
        "normalized_momentum_sample_max",
        "normalized_momentum_grid_l2",
        "normalized_momentum_volume_l2",
    ]
    assert protocol["reported_divergence_metrics"] == [
        "normalized_divergence_sample_max",
        "normalized_divergence_grid_l2",
        "normalized_divergence_volume_l2",
    ]
    assert protocol["correction_divergence_gate"] == FINAL_NORMALIZED_DIVERGENCE_GATE
    assert protocol["held_in_momentum_growth_rejected"] is True
    assert protocol["held_out_momentum_growth_rejected"] is True
    assert protocol["held_in_divergence_growth_rejected"] is True
    assert protocol["held_out_divergence_growth_rejected"] is True
    assert protocol["caller_supplied_residual_or_gain_allowed"] is False
    assert protocol["caller_supplied_pressure_or_forcing_allowed"] is False
    assert protocol["caller_supplied_scientific_threshold_allowed"] is False
    assert protocol["analytic_kokuno_gain_bindings_complete"] is False
    assert protocol["analytic_kokuno_gain_available"] is False
    assert protocol["observed_gain_is_analytic_kokuno_bound"] is False
    assert protocol["mechanics_only_receipt_is_candidate_residual_evidence"] is False
    assert protocol["mechanics_only_receipt_can_set_pde_validated"] is False


def test_current_source_and_a4_context_remain_fail_closed() -> None:
    payload = _payload()
    current = payload["current_source_admission"]
    assert current["complete_ns_defect"] is False
    assert current["pressure_gradient_included"] is False
    assert current["restricted_forcing_included"] is False
    assert current["scoped_transport_stress_authorized_as_correction_target"] is False
    assert current["admitted_to_real_finite_correction_cycle"] is False
    assert current["surrogate_defect_substitution_allowed"] is False
    assert current["residual_as_forcing_shortcut_allowed"] is False
    assert current["real_candidate_finite_correction_cycle_run"] is False
    assert current["heldout_normalized_ns_residual_assessed"] is False
    assert current["residual_reduction_claimed"] is False
    assert current["same_protocol_comparable_to_st006"] is False

    a4 = payload["agent4_context"]
    assert a4["latest_relevant_candidate_audit_pr"] == 889
    assert a4["latest_relevant_candidate_audit_head"] == (
        "574fd5f7d67e5fbf905ef91a483a7e9e8124698d"
    )
    assert a4["latest_relevant_candidate_audit_dedicated_run"] == 35530344739
    assert a4["latest_relevant_candidate_audit_tests_run"] == 35530344662
    assert a4["dedicated_finite_correction_stage_audit_present"] is False
    assert a4["finite_correction_stage_audit_conclusion"] is None
    assert a4["latest_a4_scope_is_finite_correction_stage_audit"] is False
    assert a4["latest_a4_scope_is_final_pde_validation"] is False


def test_parent_candidate_api_baseline_and_gates_are_unchanged() -> None:
    payload = _payload()
    parent = deterministic_strict_inner_vorticity_artifact_ingest_contract(
        exact_head=PARENT_A5_HEAD
    )
    assert payload["candidate_api_handoff"] == parent["candidate_api_handoff"]
    assert payload["candidate_api_handoff"]["pressure"] is None
    assert payload["candidate_api_handoff"]["forcing"] is None
    assert payload["candidate_api_handoff"]["complete_candidate_api_ready"] is False
    assert payload["baseline_vs_kokuno"] == parent["baseline_vs_kokuno"]
    assert payload["final_project_gates_unchanged"] == parent["final_project_gates_unchanged"]
    assert FINAL_NORMALIZED_MOMENTUM_GATE == 1e-3
    assert FINAL_NORMALIZED_DIVERGENCE_GATE == 1e-5

    truth = payload["truth_boundary"]
    assert truth["agent3_mechanics_receipt_laundered_as_candidate_evidence"] is False
    assert truth["agent4_stage_audit_invented"] is False
    assert truth["complete_ns_defect_invented"] is False
    assert truth["matched_pressure_invented"] is False
    assert truth["restricted_forcing_invented"] is False
    assert truth["correction_velocity_invented"] is False
    assert truth["outer_global_join_invented"] is False
    assert truth["heldout_sample_leakage_allowed"] is False
    assert truth["residual_protocol_drift_allowed"] is False
    assert truth["free_residual_defined_forcing_allowed"] is False
    assert truth["threshold_relaxed"] is False
    assert truth["pde_validated"] is False


@pytest.mark.parametrize(
    "mutate",
    [
        lambda p: p["evidence"].__setitem__("agent3_exact_head_ci_conclusion", "success"),
        lambda p: p["agent3_binding"].__setitem__("head", "1" * 40),
        lambda p: p["agent3_binding"].__setitem__("source_blob_sha", "2" * 40),
        lambda p: p["agent3_binding"].__setitem__("workflow_blob_sha", "3" * 40),
        lambda p: p["agent3_binding"]["protocol"].__setitem__("correction_divergence_gate", 2e-5),
        lambda p: p["agent3_binding"]["protocol"].__setitem__("held_out_momentum_growth_rejected", False),
        lambda p: p["agent3_binding"]["protocol"].__setitem__("caller_supplied_pressure_or_forcing_allowed", True),
        lambda p: p["agent4_context"].__setitem__("dedicated_finite_correction_stage_audit_present", True),
        lambda p: p["current_source_admission"].__setitem__("complete_ns_defect", True),
        lambda p: p["current_source_admission"].__setitem__("admitted_to_real_finite_correction_cycle", True),
        lambda p: p["current_source_admission"].__setitem__("surrogate_defect_substitution_allowed", True),
        lambda p: p["finite_correction_stage_handoff"].__setitem__("usable_for_real_cycle_now", True),
        lambda p: p["finite_correction_stage_handoff"].__setitem__("usable_as_st006_comparison_now", True),
        lambda p: p["artifact_status"].__setitem__("real_agent3_correction_velocity_materialized", True),
        lambda p: p["ingest_status"].__setitem__("real_finite_correction_cycle_admitted", True),
        lambda p: p["stage_state"].__setitem__("correction_ready", True),
        lambda p: p["stage_state"].__setitem__("pde_validated", True),
        lambda p: p["truth_boundary"].__setitem__("heldout_sample_leakage_allowed", True),
        lambda p: p["truth_boundary"].__setitem__("residual_protocol_drift_allowed", True),
        lambda p: p["truth_boundary"].__setitem__("free_residual_defined_forcing_allowed", True),
        lambda p: p["truth_boundary"].__setitem__("threshold_relaxed", True),
    ],
)
def test_truth_mutations_fail_closed(mutate) -> None:
    payload = copy.deepcopy(_payload())
    mutate(payload)
    with pytest.raises(ValueError, match="drifted"):
        validate_finite_correction_stage_ingest_contract(payload)


def test_invalid_exact_head_fails_closed() -> None:
    payload = deterministic_finite_correction_stage_ingest_contract(exact_head=None)
    payload["exact_head"] = "not-a-commit"
    with pytest.raises(ValueError, match="40-character"):
        validate_finite_correction_stage_ingest_contract(payload)
