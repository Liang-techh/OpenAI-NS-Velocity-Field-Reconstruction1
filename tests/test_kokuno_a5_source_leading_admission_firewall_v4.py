from __future__ import annotations

from copy import deepcopy
from dataclasses import replace

import pytest

from openai_ns_reconstruction.kokuno_a5_source_leading_admission_firewall_v4 import (
    A4_AVERAGED_ETA_R1_RECEIPT,
    FINAL_PROJECT_GATES,
    FROZEN_EVIDENCE,
    SCHEMA,
    ST006_BASELINE,
    derive_state,
    deterministic_receipt,
    validate_receipt,
)


def test_frozen_receipt_admits_only_exact_head_r1_averaged_eta_increment():
    receipt = deterministic_receipt(exact_head="test-head")
    assert receipt["schema"] == SCHEMA
    assert receipt["exact_head"] == "test-head"

    new = receipt["new_independent_receipt"]
    assert new == A4_AVERAGED_ETA_R1_RECEIPT
    assert new["workflow_event"] == "push"
    assert new["workflow_conclusion"] == "success"
    assert new["head"] == "c4515a3c7803b7fa28d06f5b546cb943804dcd35"
    assert new["workflow_id"] == 35474229697
    assert new["artifact_id"] == 10595470641
    assert new["artifact_digest"] == (
        "sha256:65907488f0cab2fe431e6e0786152ee13e5656a227009e7b1d05df2bee6df158"
    )

    state = receipt["derived_state"]
    assert state["source_r2_all_ordinary_slots_independently_admitted"] is True
    assert state["source_r1_zeta_ordinary_slot_independently_admitted"] is True
    assert state["source_r1_averaged_eta_ordinary_slot_independently_admitted"] is True
    assert state["source_r1_wstar_logradial_ordinary_slot_independently_admitted"] is False
    assert state["source_r1_averaged_logradial_ordinary_slot_independently_admitted"] is False
    assert state["source_r1_all_ordinary_slots_independently_admitted"] is False
    assert state["source_full_post_J2_R1_machine_bound"] is False
    assert state["leading_ready"] is False
    assert state["correction_ready"] is False
    assert state["velocity_export_ready"] is False
    assert state["pde_validated"] is False

    assert receipt["final_project_gates_unchanged"] == FINAL_PROJECT_GATES
    assert (
        receipt["baseline_vs_kokuno"]["retained_repository_baseline"]
        == ST006_BASELINE
    )
    assert (
        receipt["baseline_vs_kokuno"][
            "kokuno_same_protocol_full_candidate_residual_available"
        ]
        is False
    )
    validate_receipt(receipt)


def test_removing_r1_averaged_eta_independent_admission_fails_closed():
    evidence = dict(FROZEN_EVIDENCE)
    evidence["r1_averaged_eta_ordinary"] = replace(
        evidence["r1_averaged_eta_ordinary"],
        workflow_conclusion_at_freeze="failure",
        independently_admitted=False,
    )
    state = derive_state(evidence)
    assert state["source_r1_zeta_ordinary_slot_independently_admitted"] is True
    assert state["source_r1_averaged_eta_ordinary_slot_independently_admitted"] is False
    assert state["source_r1_all_ordinary_slots_independently_admitted"] is False
    assert state["leading_ready"] is False
    assert state["pde_validated"] is False


def test_two_green_r1_slots_cannot_promote_full_r1_or_leading():
    state = derive_state()
    forbidden = (
        "source_r1_all_ordinary_slots_independently_admitted",
        "source_full_post_J2_R1_machine_bound",
        "source_operator_MK_machine_bound",
        "source_B0_Tsh_machine_bound",
        "selected_PA16_handoff_allowed",
        "global_matched_pressure_materialized",
        "global_leading_velocity_materialized",
        "leading_ready",
        "correction_ready",
        "velocity_export_ready",
        "formal_heldout_pde_gate_assessed",
        "pde_validated",
    )
    assert all(state[name] is False for name in forbidden)


def test_pending_a4_and_agent1_progress_remain_fail_closed():
    receipt = deterministic_receipt()
    pending = receipt["pending_agent4_r1_audits"]
    assert pending["wstar_logradial"]["workflow_conclusion_at_freeze"] == "queued"
    assert pending["averaged_logradial"][
        "independent_agent4_audit_available_at_freeze"
    ] is False

    progress = receipt["upstream_agent1_r1_progress"]
    assert progress["latest_pr"] == 733
    assert progress["ordinary_r1_source_bound_slots"] == 9
    assert progress["ordinary_r1_total_slots"] == 11
    assert progress["ordinary_r1_remaining_slots"] == 2
    assert progress["remaining_slot_names"] == [
        "H_* Phi_eta",
        "Lambda^-1 d u Phi_eta",
    ]
    assert progress["agent1_self_certificate_is_not_independent_admission"] is True
    assert receipt["derived_state"]["source_r1_all_ordinary_slots_independently_admitted"] is False


def test_validator_rejects_broader_state_promotion_or_receipt_mutation():
    receipt = deterministic_receipt(exact_head="test-head")

    promoted = deepcopy(receipt)
    promoted["derived_state"]["leading_ready"] = True
    with pytest.raises(RuntimeError, match="fail-closed states were promoted"):
        validate_receipt(promoted)

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

    fake_full_r1 = deepcopy(receipt)
    fake_full_r1["derived_state"]["source_r1_all_ordinary_slots_independently_admitted"] = True
    with pytest.raises(RuntimeError, match="fail-closed states were promoted"):
        validate_receipt(fake_full_r1)
