from __future__ import annotations

from copy import deepcopy
from dataclasses import replace

import pytest

from openai_ns_reconstruction.kokuno_a5_source_leading_admission_firewall_v3 import (
    A4_ZETA_R1_RECEIPT,
    FINAL_PROJECT_GATES,
    FROZEN_EVIDENCE,
    SCHEMA,
    ST006_BASELINE,
    derive_state,
    deterministic_receipt,
    validate_receipt,
)


def test_frozen_receipt_admits_only_exact_head_r1_zeta_increment():
    receipt = deterministic_receipt(exact_head="test-head")
    assert receipt["schema"] == SCHEMA
    assert receipt["exact_head"] == "test-head"

    new = receipt["new_independent_receipt"]
    assert new == A4_ZETA_R1_RECEIPT
    assert new["workflow_event"] == "push"
    assert new["workflow_conclusion"] == "success"
    assert new["head"] == "a6616ec2cdd21948b1e45a1cccc27e5cf3206a46"
    assert new["artifact_id"] == 10594044299

    state = receipt["derived_state"]
    assert state["source_r2_all_ordinary_slots_independently_admitted"] is True
    assert state["source_r1_zeta_ordinary_slot_independently_admitted"] is True
    assert state["source_r1_averaged_eta_ordinary_slot_independently_admitted"] is False
    assert state["source_r1_wstar_logradial_ordinary_slot_independently_admitted"] is False
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


def test_removing_r1_zeta_independent_admission_fails_closed():
    evidence = dict(FROZEN_EVIDENCE)
    evidence["r1_zeta_ordinary"] = replace(
        evidence["r1_zeta_ordinary"],
        workflow_conclusion_at_freeze="failure",
        independently_admitted=False,
    )
    state = derive_state(evidence)
    assert state["source_r2_all_ordinary_slots_independently_admitted"] is True
    assert state["source_r1_zeta_ordinary_slot_independently_admitted"] is False
    assert state["source_r1_all_ordinary_slots_independently_admitted"] is False
    assert state["leading_ready"] is False
    assert state["pde_validated"] is False


def test_one_green_r1_slot_cannot_promote_full_r1_or_leading():
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


def test_pending_a4_audits_remain_pending_and_agent1_progress_is_not_admission():
    receipt = deterministic_receipt()
    pending = receipt["pending_agent4_r1_audits"]
    assert pending["averaged_eta"]["workflow_conclusion_at_freeze"] == "queued"
    assert pending["wstar_logradial"]["workflow_conclusion_at_freeze"] == "queued"

    progress = receipt["upstream_agent1_r1_progress"]
    assert progress["ordinary_r1_source_bound_slots"] == 8
    assert progress["ordinary_r1_total_slots"] == 11
    assert progress["ordinary_r1_remaining_slots"] == 3
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
