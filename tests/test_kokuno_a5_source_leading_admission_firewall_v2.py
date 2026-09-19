from __future__ import annotations

from dataclasses import replace
from copy import deepcopy

import pytest

from openai_ns_reconstruction.kokuno_a5_source_leading_admission_firewall_v2 import (
    A4_ALGEBRAIC_R2_RECEIPT,
    A4_DERIVATIVE_R2_RECEIPT,
    FINAL_PROJECT_GATES,
    FROZEN_EVIDENCE,
    SCHEMA,
    derive_state,
    deterministic_receipt,
    validate_receipt,
)


def test_exact_agent4_audits_complete_only_the_ordinary_r2_subset():
    receipt = deterministic_receipt(exact_head="test-head")
    assert receipt["schema"] == SCHEMA
    assert receipt["exact_head"] == "test-head"
    assert receipt["new_independent_receipts"] == {
        "algebraic_r2_ordinary": A4_ALGEBRAIC_R2_RECEIPT,
        "derivative_r2_ordinary": A4_DERIVATIVE_R2_RECEIPT,
    }

    evidence = receipt["evidence"]
    algebraic = evidence["algebraic_r2_ordinary"]
    assert algebraic["producer"] == "Agent 4"
    assert algebraic["pr"] == 695
    assert algebraic["head"] == "c38112a8d0172825119ef4f174b3b57af3a369ea"
    assert algebraic["workflow_id"] == 35465298813
    assert algebraic["workflow_conclusion_at_freeze"] == "success"
    assert algebraic["independently_admitted"] is True

    derivative = evidence["derivative_r2_ordinary"]
    assert derivative["producer"] == "Agent 4"
    assert derivative["pr"] == 704
    assert derivative["head"] == "a96bef02c0cb5c2fdc292d4dc340e5eb464f8d42"
    assert derivative["workflow_id"] == 35468474402
    assert derivative["workflow_conclusion_at_freeze"] == "success"
    assert derivative["independently_admitted"] is True

    state = receipt["derived_state"]
    assert state["source_r2_algebraic_ordinary_slots_independently_admitted"] is True
    assert state["source_r2_derivative_ordinary_slots_independently_admitted"] is True
    assert state["source_r2_all_ordinary_slots_independently_admitted"] is True
    assert state["source_full_post_J1_R2_independently_admitted"] is False
    assert state["all_R1_ordinary_slots_source_bound"] is False
    assert state["leading_ready"] is False
    assert state["correction_ready"] is False
    assert state["velocity_export_ready"] is False
    assert state["pde_validated"] is False
    validate_receipt(receipt)


def test_removing_either_independent_subset_fails_closed_for_all_ordinary_r2():
    for key in ("algebraic_r2_ordinary", "derivative_r2_ordinary"):
        evidence = dict(FROZEN_EVIDENCE)
        evidence[key] = replace(evidence[key], independently_admitted=False)
        state = derive_state(evidence)
        assert state["source_r2_all_ordinary_slots_independently_admitted"] is False
        assert state["source_full_post_J1_R2_independently_admitted"] is False
        assert state["leading_ready"] is False


def test_complete_ordinary_subset_cannot_launder_full_r2_or_leading():
    state = derive_state()
    assert state["source_r2_all_ordinary_slots_independently_admitted"] is True
    assert state["source_full_post_J1_R2_independently_admitted"] is False
    assert state["source_full_post_J2_R1_machine_bound"] is False
    assert state["source_operator_MK_machine_bound"] is False
    assert state["source_B0_Tsh_machine_bound"] is False
    assert state["global_matched_pressure_materialized"] is False
    assert state["global_leading_velocity_materialized"] is False
    assert state["leading_ready"] is False
    assert state["pde_validated"] is False


def test_agent1_self_certificate_cannot_replace_either_agent4_receipt():
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
    bad["new_independent_receipts"]["algebraic_r2_ordinary"]["artifact_digest"] = (
        "sha256:" + "0" * 64
    )
    with pytest.raises(RuntimeError, match="receipts changed"):
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
    assert truth["single_integration_increment"] == (
        "admit independently audited ordinary-R2 subsets"
    )
    assert truth["ordinary_r2_subset_is_not_full_post_J1_R2"] is True
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
