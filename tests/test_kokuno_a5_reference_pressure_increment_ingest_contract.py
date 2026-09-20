from __future__ import annotations

import copy

import pytest

from openai_ns_reconstruction.kokuno_a5_reference_pressure_increment_ingest_contract import (
    AGENT1_PRESSURE_INCREMENT,
    AGENT4_PRESSURE_INCREMENT_AUDIT,
    FINAL_GATE,
    PARENT_A5,
    READINESS,
    ST006_BASELINE,
    build_contract,
    validate_contract,
)


def _contract():
    return build_contract(exact_head="f" * 40)


def test_registered_pressure_increment_and_independent_audit_are_scoped() -> None:
    contract = _contract()
    validate_contract(contract)

    assert contract["parent_a5"] == PARENT_A5
    a1 = contract["agent1_pressure_increment"]
    assert a1["head"] == AGENT1_PRESSURE_INCREMENT["head"]
    assert a1["source_blob"] == AGENT1_PRESSURE_INCREMENT["source_blob"]
    assert a1["surface"] == "C_p_reference(X,eta)=Pi_r(X,eta)-Pi_0(eta)"
    assert a1["source_determined_pressure_increment"] is True
    assert a1["absolute_axis_pressure_Pi0_materialized"] is False
    assert a1["absolute_reference_pressure_materialized"] is False
    assert a1["matched_global_pressure_materialized"] is False
    assert a1["cartesian_pressure_gradient_materialized"] is False

    a4 = contract["agent4_pressure_increment_audit"]
    assert a4["head"] == AGENT4_PRESSURE_INCREMENT_AUDIT["head"]
    assert a4["source_blob"] == AGENT4_PRESSURE_INCREMENT_AUDIT["source_blob"]
    assert a4["test_blob"] == AGENT4_PRESSURE_INCREMENT_AUDIT["test_blob"]
    assert a4["dedicated_workflow_present"] is False
    assert a4["saved_reloaded_public_surface_only"] is True
    assert a4["implementation_distinct_from_agent1_production_derivatives"] is True
    assert a4["reference_pressure_increment_independently_audited"] is False
    assert a4["cartesian_pressure_gradient_assessed"] is False
    assert a4["complete_ns_residual_assessed"] is False


def test_ci_is_recorded_as_unresolved_not_admitted() -> None:
    contract = _contract()
    a1_ci = contract["agent1_pressure_increment"]["observed_ci"]
    a4_ci = contract["agent4_pressure_increment_audit"]["observed_ci"]
    assert a1_ci["dedicated"] == {"run_id": 35539131381, "status": "queued", "conclusion": None}
    assert a1_ci["repository_tests"] == {"run_id": 35539131366, "status": "queued", "conclusion": None}
    assert a4_ci["repository_tests"] == {"run_id": 35539845126, "status": "queued", "conclusion": None}
    assert contract["agent1_pressure_increment"]["scientific_admission"] is False
    assert contract["agent4_pressure_increment_audit"]["scientific_admission"] is False
    assert contract["pressure_ingest_seam"]["queued_or_unresolved_ci_is_pass"] is False


def test_readiness_baseline_and_final_gates_do_not_move() -> None:
    contract = _contract()
    assert contract["readiness"] == READINESS
    assert contract["final_gate"] == FINAL_GATE
    assert contract["baseline"]["st006"] == ST006_BASELINE
    assert contract["readiness"] == {
        "leading_ready": False,
        "oscillatory_ready": True,
        "correction_ready": False,
        "velocity_export_ready": False,
        "pde_validated": False,
    }
    assert contract["pressure_ingest_seam"]["complete_candidate_pressure_api_ready"] is False
    assert contract["pressure_ingest_seam"]["current_candidate_eligible_for_full_ns_validation"] is False
    assert contract["pressure_ingest_seam"]["same_protocol_st006_comparison_available_now"] is False


@pytest.mark.parametrize(
    ("path", "value"),
    [
        (("agent1_pressure_increment", "absolute_axis_pressure_Pi0_materialized"), True),
        (("agent1_pressure_increment", "matched_global_pressure_materialized"), True),
        (("agent4_pressure_increment_audit", "reference_pressure_increment_independently_audited"), True),
        (("agent4_pressure_increment_audit", "scientific_admission"), True),
        (("pressure_ingest_seam", "pressure_increment_is_absolute_pressure"), True),
        (("pressure_ingest_seam", "pressure_increment_authorized_as_matched_global_pressure"), True),
        (("pressure_ingest_seam", "pressure_increment_authorized_as_cartesian_grad_p"), True),
        (("pressure_ingest_seam", "residual_defined_forcing_forbidden"), False),
        (("readiness", "leading_ready"), True),
        (("readiness", "pde_validated"), True),
        (("truth_boundary", "matched_pressure_present"), True),
        (("truth_boundary", "same_protocol_st006_improvement_claimed"), True),
        (("truth_boundary", "paper_exact"), True),
    ],
)
def test_truth_promotions_fail_closed(path: tuple[str, str], value: object) -> None:
    mutated = copy.deepcopy(_contract())
    mutated[path[0]][path[1]] = value
    with pytest.raises(ValueError):
        validate_contract(mutated)


def test_provenance_and_gate_mutations_fail_closed() -> None:
    for mutate in (
        lambda c: c["agent1_pressure_increment"].__setitem__("head", "0" * 40),
        lambda c: c["agent4_pressure_increment_audit"].__setitem__("source_blob", "0" * 40),
        lambda c: c["final_gate"].__setitem__("momentum_volume_l2", 1.001e-3),
        lambda c: c["final_gate"].__setitem__("divergence_volume_l2", 1.001e-5),
        lambda c: c["baseline"]["st006"].__setitem__("momentum_volume_l2", 0.1),
    ):
        mutated = copy.deepcopy(_contract())
        mutate(mutated)
        with pytest.raises(ValueError):
            validate_contract(mutated)


def test_checksum_tamper_is_rejected() -> None:
    mutated = copy.deepcopy(_contract())
    mutated["contract_sha256"] = "0" * 64
    with pytest.raises(ValueError, match="checksum"):
        validate_contract(mutated)
