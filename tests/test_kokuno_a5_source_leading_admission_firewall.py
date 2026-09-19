from __future__ import annotations

from dataclasses import replace

import pytest

from openai_ns_reconstruction.kokuno_a5_source_leading_admission_firewall import (
    FINAL_PROJECT_GATES,
    FROZEN_EVIDENCE,
    SCHEMA,
    derive_state,
    deterministic_receipt,
    validate_receipt,
)


def test_frozen_receipt_admits_only_green_independent_prerequisites():
    receipt = deterministic_receipt(exact_head="test-head")
    assert receipt["schema"] == SCHEMA
    assert receipt["exact_head"] == "test-head"
    evidence = receipt["evidence"]

    assert evidence["source_axis"]["independently_admitted"] is True
    assert evidence["source_axis"]["workflow_conclusion_at_freeze"] == "success"
    assert evidence["pressure_r2_ordinary"]["independently_admitted"] is True
    assert evidence["source_u_ball"]["independently_admitted"] is True

    assert evidence["algebraic_r2_ordinary"]["independently_admitted"] is False
    assert evidence["algebraic_r2_ordinary"]["workflow_conclusion_at_freeze"] == "queued"
    assert evidence["derivative_r2_ordinary"]["independently_admitted"] is False
    assert evidence["derivative_r2_ordinary"]["workflow_conclusion_at_freeze"] == "queued"

    validate_receipt(receipt)


def test_one_pending_r2_subset_cannot_promote_all_ordinary_r2():
    evidence = dict(FROZEN_EVIDENCE)
    evidence["algebraic_r2_ordinary"] = replace(
        evidence["algebraic_r2_ordinary"], independently_admitted=True
    )
    state = derive_state(evidence)
    assert state["source_r2_algebraic_ordinary_slots_independently_admitted"] is True
    assert state["source_r2_derivative_ordinary_slots_independently_admitted"] is False
    assert state["source_r2_all_ordinary_slots_independently_admitted"] is False
    assert state["leading_ready"] is False


def test_even_synthetic_full_ordinary_r2_does_not_launder_full_r2_or_leading():
    evidence = dict(FROZEN_EVIDENCE)
    evidence["algebraic_r2_ordinary"] = replace(
        evidence["algebraic_r2_ordinary"], independently_admitted=True
    )
    evidence["derivative_r2_ordinary"] = replace(
        evidence["derivative_r2_ordinary"], independently_admitted=True
    )
    state = derive_state(evidence)
    assert state["source_r2_all_ordinary_slots_independently_admitted"] is True
    assert state["source_full_post_J1_R2_independently_admitted"] is False
    assert state["all_R1_ordinary_slots_source_bound"] is False
    assert state["source_operator_MK_machine_bound"] is False
    assert state["leading_ready"] is False
    assert state["pde_validated"] is False


def test_missing_admission_evidence_fails_closed():
    evidence = dict(FROZEN_EVIDENCE)
    del evidence["source_u_ball"]
    with pytest.raises(ValueError, match="missing required admission evidence"):
        derive_state(evidence)


def test_fixed_project_gates_and_forcing_policy_are_unchanged():
    receipt = deterministic_receipt()
    assert receipt["final_project_gates_unchanged"] == {
        "normalized_momentum_max": 1.0e-3,
        "normalized_momentum_L2": 1.0e-3,
        "divergence_max": 1.0e-5,
        "divergence_L2": 1.0e-5,
    }
    assert receipt["final_project_gates_unchanged"] == FINAL_PROJECT_GATES
    truth = receipt["truth_boundary"]
    assert truth["free_residual_defined_forcing_allowed"] is False
    assert truth["threshold_relaxed"] is False
    assert truth["queued_workflow_treated_as_success"] is False
    assert truth["agent1_self_certificate_treated_as_independent_admission"] is False
    assert truth["kokuno_replay_used_as_final_independent_pde_validation"] is False
    assert truth["pde_validated"] is False


def test_receipt_hash_is_deterministic_and_exact_head_bound():
    a = deterministic_receipt(exact_head="abc")
    b = deterministic_receipt(exact_head="abc")
    c = deterministic_receipt(exact_head="def")
    assert a["receipt_sha256"] == b["receipt_sha256"]
    assert a["receipt_sha256"] != c["receipt_sha256"]
