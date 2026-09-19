from __future__ import annotations

from dataclasses import replace
from copy import deepcopy

import pytest

from openai_ns_reconstruction.kokuno_a5_source_leading_admission_firewall_v2 import (
    A4_DERIVATIVE_R2_RECEIPT,
    FINAL_PROJECT_GATES,
    FROZEN_EVIDENCE,
    SCHEMA,
    derive_state,
    deterministic_receipt,
    validate_receipt,
)


def test_exact_agent4_derivative_r2_pass_is_the_only_new_admission():
    receipt = deterministic_receipt(exact_head="test-head")
    assert receipt["schema"] == SCHEMA
    assert receipt["exact_head"] == "test-head"
    assert receipt["new_independent_receipt"] == A4_DERIVATIVE_R2_RECEIPT

    evidence = receipt["evidence"]
    derivative = evidence["derivative_r2_ordinary"]
    assert derivative["producer"] == "Agent 4"
    assert derivative["pr"] == 704
    assert derivative["head"] == "a96bef02c0cb5c2fdc292d4dc340e5eb464f8d42"
    assert derivative["workflow_id"] == 35468474402
    assert derivative["workflow_conclusion_at_freeze"] == "success"
    assert derivative["independently_admitted"] is True

    algebraic = evidence["algebraic_r2_ordinary"]
    assert algebraic["pr"] == 695
    assert algebraic["independently_admitted"] is False

    state = receipt["derived_state"]
    assert state["source_r2_derivative_ordinary_slots_independently_admitted"] is True
    assert state["source_r2_algebraic_ordinary_slots_independently_admitted"] is False
    assert state["source_r2_all_ordinary_slots_independently_admitted"] is False
    assert state["source_full_post_J1_R2_independently_admitted"] is False
    assert state["leading_ready"] is False
    assert state["correction_ready"] is False
    assert state["velocity_export_ready"] is False
    assert state["pde_validated"] is False
    validate_receipt(receipt)


def test_future_algebraic_pass_would_only_complete_ordinary_subset_not_full_r2():
    evidence = dict(FROZEN_EVIDENCE)
    evidence["algebraic_r2_ordinary"] = replace(
        evidence["algebraic_r2_ordinary"],
        workflow_conclusion_at_freeze="success",
        independently_admitted=True,
    )
    state = derive_state(evidence)
    assert state["source_r2_all_ordinary_slots_independently_admitted"] is True
    assert state["source_full_post_J1_R2_independently_admitted"] is False
    assert state["source_operator_MK_machine_bound"] is False
    assert state["leading_ready"] is False
    assert state["pde_validated"] is False


def test_derivative_admission_cannot_be_inferred_from_agent1_self_certificate():
    evidence = dict(FROZEN_EVIDENCE)
    evidence["derivative_r2_ordinary"] = replace(
        evidence["derivative_r2_ordinary"],
        producer="Agent 1",
        pr=693,
        workflow_id=35464666115,
        workflow_conclusion_at_freeze="success",
        independently_admitted=False,
    )
    state = derive_state(evidence)
    assert state["source_r2_derivative_ordinary_slots_independently_admitted"] is False
    assert state["source_r2_all_ordinary_slots_independently_admitted"] is False
    assert state["leading_ready"] is False


def test_receipt_rejects_mutated_independent_provenance():
    receipt = deterministic_receipt()
    bad = deepcopy(receipt)
    bad["new_independent_receipt"]["artifact_digest"] = "sha256:" + "0" * 64
    with pytest.raises(RuntimeError, match="receipt changed"):
        validate_receipt(bad)


def test_fixed_project_gates_and_no_free_forcing_remain_unchanged():
    receipt = deterministic_receipt()
    assert receipt["final_project_gates_unchanged"] == {
        "normalized_momentum_max": 1.0e-3,
        "normalized_momentum_L2": 1.0e-3,
        "divergence_max": 1.0e-5,
        "divergence_L2": 1.0e-5,
    }
    assert receipt["final_project_gates_unchanged"] == FINAL_PROJECT_GATES
    truth = receipt["truth_boundary"]
    assert truth["single_promotion"] == (
        "source_r2_derivative_ordinary_slots_independently_admitted"
    )
    assert truth["free_residual_defined_forcing_allowed"] is False
    assert truth["threshold_relaxed"] is False
    assert truth["agent1_self_certificate_treated_as_independent_admission"] is False
    assert truth["kokuno_replay_used_as_final_independent_pde_validation"] is False
    assert truth["pde_validated"] is False


def test_receipt_hash_is_deterministic_and_exact_head_bound():
    a = deterministic_receipt(exact_head="abc")
    b = deterministic_receipt(exact_head="abc")
    c = deterministic_receipt(exact_head="def")
    assert a["receipt_sha256"] == b["receipt_sha256"]
    assert a["receipt_sha256"] != c["receipt_sha256"]
