from __future__ import annotations

from copy import deepcopy
from dataclasses import replace

import pytest

from openai_ns_reconstruction.kokuno_a5_source_leading_admission_firewall_v5 import (
    A4_WSTAR_R1_RECEIPT,
    FINAL_PROJECT_GATES,
    FROZEN_EVIDENCE,
    SCHEMA,
    ST006_BASELINE,
    derive_state,
    deterministic_receipt,
    validate_receipt,
)


def test_exact_head_wstar_receipt_admits_only_wstar_increment():
    receipt = deterministic_receipt(exact_head="test-head")
    assert receipt["schema"] == SCHEMA
    assert receipt["exact_head"] == "test-head"

    new = receipt["new_independent_receipt"]
    assert new == A4_WSTAR_R1_RECEIPT
    assert new["workflow_conclusion"] == "success"
    assert new["head"] == "279de4c794935bcb3331627046be3c96b865bac7"
    assert new["workflow_id"] == 35477137066
    assert new["artifact_id"] == 10596193830
    assert new["artifact_digest"] == (
        "sha256:08b514794efa085297ce9ab14d6e8eebc0e9d737b1a753e3ffd2dfd03baaa2b7"
    )
    assert new["report_sha256"] == (
        "76488aa1440e9651ded65480a529eac7b1fed86ea299e698eb5e6d8c4a315c6e"
    )
    assert new["receipt_sha256"] == (
        "f0816da8a9c1264ef6822ed99034bfed9a15e6a5d1200fae1d9b378426da60ec"
    )
    assert new["failed_guards"] == []
    assert new["independent_preflight_passed"] is True

    state = receipt["derived_state"]
    assert state["source_r2_all_ordinary_slots_independently_admitted"] is True
    assert state["source_r1_zeta_ordinary_slot_independently_admitted"] is True
    assert state["source_r1_averaged_eta_ordinary_slot_independently_admitted"] is True
    assert state["source_r1_wstar_logradial_ordinary_slot_independently_admitted"] is True
    assert state["source_r1_averaged_logradial_ordinary_slot_independently_admitted"] is False
    assert state["source_r1_hstar_phi_eta_ordinary_slot_independently_admitted"] is False
    assert state["source_r1_du_phi_eta_ordinary_slot_independently_admitted"] is False
    assert state["source_r1_all_ordinary_slots_independently_admitted"] is False
    assert state["source_full_post_J2_R1_machine_bound"] is False
    assert state["leading_ready"] is False
    assert state["correction_ready"] is False
    assert state["velocity_export_ready"] is False
    assert state["pde_validated"] is False

    assert receipt["final_project_gates_unchanged"] == FINAL_PROJECT_GATES
    assert receipt["baseline_vs_kokuno"]["retained_repository_baseline"] == ST006_BASELINE
    assert (
        receipt["baseline_vs_kokuno"][
            "kokuno_same_protocol_full_candidate_residual_available"
        ]
        is False
    )
    validate_receipt(receipt)


def test_removing_wstar_independent_admission_fails_closed():
    evidence = dict(FROZEN_EVIDENCE)
    evidence["r1_wstar_logradial_ordinary"] = replace(
        evidence["r1_wstar_logradial_ordinary"],
        workflow_conclusion_at_freeze="failure",
        independently_admitted=False,
    )
    state = derive_state(evidence)
    assert state["source_r1_zeta_ordinary_slot_independently_admitted"] is True
    assert state["source_r1_averaged_eta_ordinary_slot_independently_admitted"] is True
    assert state["source_r1_wstar_logradial_ordinary_slot_independently_admitted"] is False
    assert state["source_r1_all_ordinary_slots_independently_admitted"] is False
    assert state["leading_ready"] is False
    assert state["pde_validated"] is False


def test_11_of_11_agent1_self_ledger_does_not_promote_full_r1():
    receipt = deterministic_receipt()
    progress = receipt["upstream_agent1_r1_progress"]
    assert progress["latest_pr"] == 751
    assert progress["ordinary_r1_source_bound_slots"] == 11
    assert progress["ordinary_r1_total_slots"] == 11
    assert progress["ordinary_r1_remaining_slots"] == 0
    assert progress["remaining_slot_names"] == []
    assert progress["agent1_full_post_J2_R1_self_envelope_materialized"] is True
    assert progress["agent1_self_certificate_is_not_independent_admission"] is True
    assert progress["agent1_dedicated_workflow_conclusion_at_freeze"] == "queued"

    state = receipt["derived_state"]
    assert state["source_r1_all_ordinary_slots_independently_admitted"] is False
    assert state["source_full_post_J2_R1_machine_bound"] is False
    assert state["leading_ready"] is False


def test_unresolved_a4_r1_audits_remain_fail_closed():
    receipt = deterministic_receipt()
    pending = receipt["pending_agent4_r1_audits"]
    assert pending["averaged_logradial"]["workflow_id"] == 35479780213
    assert pending["averaged_logradial"]["workflow_conclusion_at_freeze"] == "queued"
    assert pending["hstar_phi_eta"]["workflow_id"] == 35482164500
    assert pending["hstar_phi_eta"]["workflow_conclusion_at_freeze"] == "queued"
    assert pending["du_phi_eta"]["workflow_id"] == 35484899251
    assert pending["du_phi_eta"]["workflow_conclusion_at_freeze"] == "queued"
    assert pending["earlier_algebraic_r1_slots"][
        "independent_agent4_coverage_complete_at_freeze"
    ] is False

    parent_ci = receipt["parent_a5_ci_at_freeze"]
    assert parent_ci["dedicated_workflow_id"] == 35482609891
    assert parent_ci["dedicated_workflow_conclusion_at_freeze"] == "queued"
    assert parent_ci["parent_queued_ci_treated_as_success"] is False


def test_validator_rejects_broader_state_promotion_or_receipt_mutation():
    receipt = deterministic_receipt(exact_head="test-head")

    promoted = deepcopy(receipt)
    promoted["derived_state"]["leading_ready"] = True
    with pytest.raises(RuntimeError, match="fail-closed states were promoted"):
        validate_receipt(promoted)

    fake_full_r1 = deepcopy(receipt)
    fake_full_r1["derived_state"][
        "source_r1_all_ordinary_slots_independently_admitted"
    ] = True
    with pytest.raises(RuntimeError, match="fail-closed states were promoted"):
        validate_receipt(fake_full_r1)

    mutated = deepcopy(receipt)
    mutated["new_independent_receipt"]["workflow_id"] = 1
    with pytest.raises(RuntimeError, match="receipt changed"):
        validate_receipt(mutated)

    fake_residual = deepcopy(receipt)
    fake_residual["baseline_vs_kokuno"][
        "kokuno_same_protocol_full_candidate_residual_available"
    ] = True
    with pytest.raises(RuntimeError, match="nonexistent full Kokuno residual"):
        validate_receipt(fake_residual)

    parent_launder = deepcopy(receipt)
    parent_launder["parent_a5_ci_at_freeze"]["parent_queued_ci_treated_as_success"] = True
    with pytest.raises(RuntimeError, match="queued parent A5 CI was promoted"):
        validate_receipt(parent_launder)
