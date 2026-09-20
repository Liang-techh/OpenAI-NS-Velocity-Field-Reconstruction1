from __future__ import annotations

import copy

import pytest

from openai_ns_reconstruction import kokuno_a5_canonical_quadrature_evidence_ingest_contract as c


HEAD = "0123456789abcdef0123456789abcdef01234567"


def fresh() -> dict:
    return c.build_contract(exact_head=HEAD)


def resign(value: dict) -> dict:
    value["contract_sha256"] = c._sha256(c._without_sha(value))
    return value


def expect_reject(mutator) -> None:
    value = copy.deepcopy(fresh())
    mutator(value)
    resign(value)
    with pytest.raises(ValueError):
        c.validate_contract(value)


def test_registered_contract_is_fail_closed_and_frozen() -> None:
    value = fresh()
    c.validate_contract(value)
    a4 = value["agent4_canonical_quadrature_handoff"]
    seam = value["admission_seam"]
    assert value["parent_a5"]["head"] == "21ff4a86638a1d8e179e25ed25886a68a34b4e6b"
    assert a4["head"] == "065810d94a37ade8d81acd2c7f622c1ff8cc6bdd"
    assert a4["source_blob"] == "4a14a2e2bd6c06280d992802414109444d41e84b"
    assert a4["workflow_blob"] == "d8dc0e4bc971e94f5b886cc57ffc3a29bade58e6"
    assert a4["quadrature_orders_per_axis"] == [24, 48, 96]
    assert a4["real_canonical_quadrature_run_completed"] is False
    assert a4["real_receipt_available"] is False
    assert a4["scientific_admission"] is False
    assert seam["typed_receipt_required"] is True
    assert seam["current_evidence_admitted"] is False
    assert value["readiness"]["pde_validated"] is False


def test_parent_and_a4_provenance_drift_rejected() -> None:
    expect_reject(lambda x: x["parent_a5"].update(head="0" * 40))
    expect_reject(lambda x: x["agent4_canonical_quadrature_handoff"].update(head="f" * 40))
    expect_reject(lambda x: x["agent4_canonical_quadrature_handoff"].update(source_blob="f" * 40))
    expect_reject(lambda x: x["agent4_canonical_quadrature_handoff"].update(workflow_blob="f" * 40))


def test_queued_evidence_cannot_be_laundered_into_admission() -> None:
    expect_reject(lambda x: x["agent4_canonical_quadrature_handoff"].update(scientific_admission=True))
    expect_reject(lambda x: x["agent4_canonical_quadrature_handoff"].update(real_receipt_available=True))
    expect_reject(lambda x: x["agent4_canonical_quadrature_handoff"].update(real_canonical_quadrature_run_completed=True))
    expect_reject(lambda x: x["admission_seam"].update(current_evidence_admitted=True))
    expect_reject(lambda x: x["readiness"].update(pde_validated=True))


def test_identity_and_digest_requirements_cannot_be_weakened() -> None:
    for key in (
        "candidate_identity_must_match",
        "stage_identity_must_match",
        "physical_contract_identity_must_match",
        "validator_operator_identity_must_match",
        "heldout_protocol_identity_must_match",
        "receipt_digest_must_verify",
    ):
        expect_reject(lambda x, key=key: x["admission_seam"].update({key: False}))


def test_quadrature_ladder_and_norm_firewall_cannot_be_weakened() -> None:
    expect_reject(lambda x: x["agent4_canonical_quadrature_handoff"].update(quadrature_orders_per_axis=[24, 48]))
    expect_reject(lambda x: x["agent4_canonical_quadrature_handoff"].update(quadrature_stability_gate_48_to_96=0.2))
    expect_reject(lambda x: x["agent4_canonical_quadrature_handoff"].update(caller_sample_weights_are_canonical_quadrature=True))
    expect_reject(lambda x: x["agent4_canonical_quadrature_handoff"].update(bare_boolean_authorizes_pde_promotion=True))
    expect_reject(lambda x: x["admission_seam"].update(sampled_rms_cannot_be_relabelled_canonical_volume_l2=False))
    expect_reject(lambda x: x["admission_seam"].update(monte_carlo_weights_cannot_replace_canonical_quadrature=False))


def test_final_gates_and_missing_scientific_surfaces_stay_locked() -> None:
    expect_reject(lambda x: x["final_gate"].update(momentum_volume_l2=1.001e-3))
    expect_reject(lambda x: x["final_gate"].update(divergence_volume_l2=1.001e-5))
    for key in (
        "global_join_present",
        "matched_pressure_present",
        "restricted_forcing_present",
        "restricted_forcing_preregistered",
        "real_correction_velocity_present",
        "complete_candidate_api_ready",
        "current_candidate_eligible_for_full_ns_validation",
        "same_protocol_st006_comparison_available_now",
    ):
        expect_reject(lambda x, key=key: x["admission_seam"].update({key: True}))


def test_free_forcing_and_paper_exact_relabel_are_rejected() -> None:
    expect_reject(lambda x: x["admission_seam"].update(residual_defined_forcing_forbidden=False))
    expect_reject(lambda x: x["truth_boundary"].update(residual_defined_forcing_forbidden=False))
    expect_reject(lambda x: x["truth_boundary"].update(kokuno_reconstruction_is_paper_exact=True))
    expect_reject(lambda x: x["truth_boundary"].update(threshold_relaxation_allowed=True))


def test_contract_digest_detects_unsigned_mutation() -> None:
    value = fresh()
    value["readiness"]["correction_ready"] = True
    with pytest.raises(ValueError, match="readiness drift|contract SHA mismatch"):
        c.validate_contract(value)
