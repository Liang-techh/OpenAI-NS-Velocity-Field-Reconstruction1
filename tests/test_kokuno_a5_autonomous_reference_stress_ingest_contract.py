from __future__ import annotations

import copy

import pytest

from openai_ns_reconstruction.kokuno_a5_autonomous_reference_stress_ingest_contract import (
    AGENT1_AUTONOMOUS_REFERENCE_STRESS,
    AGENT4_AUTONOMOUS_REFERENCE_STRESS_AUDIT,
    FINAL_GATE,
    PARENT_A5,
    READINESS,
    ST006_BASELINE,
    _sha256,
    _without_sha,
    build_contract,
    validate_contract,
)


EXACT_HEAD = "0123456789abcdef0123456789abcdef01234567"


def _reseal(payload: dict) -> dict:
    out = copy.deepcopy(payload)
    out["contract_sha256"] = _sha256(_without_sha(out))
    return out


def test_build_contract_registers_only_autonomous_reference_stress_seam() -> None:
    c = build_contract(exact_head=EXACT_HEAD)
    assert c["exact_head"] == EXACT_HEAD
    assert c["parent_a5"] == PARENT_A5

    a1 = c["agent1_autonomous_reference_stress"]
    assert a1["head"] == AGENT1_AUTONOMOUS_REFERENCE_STRESS["head"]
    assert a1["source_blob"] == "44159b63c440ee9764261dad41fb3f9d3fec7b0d"
    assert a1["autonomous_pressure_reference_stress_pair_callable"] is True
    assert a1["source_prepared_appendixA_Pi0_materialized"] is False
    assert a1["source_prepared_full_reference_stress_pair_materialized"] is False
    assert a1["matched_global_pressure_materialized"] is False
    assert a1["scientific_admission"] is False

    a4 = c["agent4_autonomous_reference_stress_audit"]
    assert a4["head"] == AGENT4_AUTONOMOUS_REFERENCE_STRESS_AUDIT["head"]
    assert a4["source_blob"] == "4707b4de00d79f7ae9db6531ff625bece4ead2df"
    assert a4["saved_reloaded_public_values_only"] is True
    assert a4["implementation_distinct_from_agent1_radial_derivative"] is True
    assert a4["autonomous_reference_stress_independently_audited"] is False
    assert a4["complete_ns_residual_assessed"] is False

    seam = c["reference_stress_ingest_seam"]
    assert seam["autonomous_reference_stress_registered"] is True
    assert seam["independent_audit_registered"] is True
    assert seam["registration_is_scientific_admission"] is False
    assert seam["queued_or_unresolved_ci_is_pass"] is False
    assert seam["autonomous_pressure_realization_is_source_prepared_Pi0"] is False
    assert seam["actual_fixed_kappa_continuation_materialized"] is False
    assert seam["outer_global_leading_velocity_materialized"] is False
    assert seam["complete_candidate_api_ready"] is False
    assert seam["current_candidate_eligible_for_full_ns_validation"] is False
    assert seam["residual_defined_forcing_forbidden"] is True

    assert c["readiness"] == READINESS
    assert c["final_gate"] == FINAL_GATE
    assert c["baseline"] == {"st006": ST006_BASELINE}
    validate_contract(c)


def test_frozen_agent4_protocol_is_not_a_pde_gate() -> None:
    c = build_contract(exact_head=EXACT_HEAD)
    protocol = c["agent4_autonomous_reference_stress_audit"]["protocol"]
    assert protocol["seed"] == 9173591
    assert protocol["scale_aware_fd6_steps"] == [2.0e-3, 1.0e-3, 5.0e-4]
    assert protocol["Ns_relative_rms_gate"] == 5.0e-3
    assert protocol["Ns_relative_sampled_max_gate"] == 2.0e-2
    assert protocol["refinement_ratio_gate"] == 12.0
    assert c["final_gate"]["momentum_volume_l2"] == 1.0e-3
    assert c["final_gate"]["divergence_volume_l2"] == 1.0e-5
    assert c["truth_boundary"]["heldout_complete_ns_residual_assessed"] is False
    assert c["truth_boundary"]["pde_validated"] is False


@pytest.mark.parametrize(
    ("section", "key", "value"),
    [
        ("agent1_autonomous_reference_stress", "source_prepared_appendixA_Pi0_materialized", True),
        ("agent1_autonomous_reference_stress", "source_prepared_full_reference_stress_pair_materialized", True),
        ("agent1_autonomous_reference_stress", "matched_global_pressure_materialized", True),
        ("agent4_autonomous_reference_stress_audit", "autonomous_reference_stress_independently_audited", True),
        ("agent4_autonomous_reference_stress_audit", "matched_cartesian_pressure_gradient_assessed", True),
        ("agent4_autonomous_reference_stress_audit", "complete_ns_residual_assessed", True),
        ("reference_stress_ingest_seam", "queued_or_unresolved_ci_is_pass", True),
        ("reference_stress_ingest_seam", "autonomous_pressure_realization_is_source_prepared_Pi0", True),
        ("reference_stress_ingest_seam", "actual_fixed_kappa_continuation_materialized", True),
        ("reference_stress_ingest_seam", "outer_global_leading_velocity_materialized", True),
        ("reference_stress_ingest_seam", "complete_candidate_api_ready", True),
        ("reference_stress_ingest_seam", "current_candidate_eligible_for_full_ns_validation", True),
        ("truth_boundary", "matched_pressure_present", True),
        ("truth_boundary", "complete_full_ns_candidate_present", True),
        ("truth_boundary", "same_protocol_st006_improvement_claimed", True),
        ("truth_boundary", "pde_validated", True),
    ],
)
def test_rejects_checksum_valid_scientific_promotion(section: str, key: str, value: object) -> None:
    c = build_contract(exact_head=EXACT_HEAD)
    c[section][key] = value
    c = _reseal(c)
    with pytest.raises(ValueError):
        validate_contract(c)


def test_rejects_gate_readiness_and_baseline_drift_even_when_resealed() -> None:
    mutations = []

    c = build_contract(exact_head=EXACT_HEAD)
    c["final_gate"]["momentum_volume_l2"] = 1.001e-3
    mutations.append(c)

    c = build_contract(exact_head=EXACT_HEAD)
    c["final_gate"]["divergence_volume_l2"] = 1.001e-5
    mutations.append(c)

    c = build_contract(exact_head=EXACT_HEAD)
    c["readiness"]["leading_ready"] = True
    mutations.append(c)

    c = build_contract(exact_head=EXACT_HEAD)
    c["readiness"]["correction_ready"] = True
    mutations.append(c)

    c = build_contract(exact_head=EXACT_HEAD)
    c["baseline"]["st006"]["momentum_volume_l2"] = 0.1
    mutations.append(c)

    for mutated in mutations:
        with pytest.raises(ValueError):
            validate_contract(_reseal(mutated))


def test_rejects_upstream_provenance_or_ci_snapshot_drift() -> None:
    mutations = []

    c = build_contract(exact_head=EXACT_HEAD)
    c["agent1_autonomous_reference_stress"]["head"] = "f" * 40
    mutations.append(c)

    c = build_contract(exact_head=EXACT_HEAD)
    c["agent4_autonomous_reference_stress_audit"]["source_blob"] = "e" * 40
    mutations.append(c)

    c = build_contract(exact_head=EXACT_HEAD)
    c["agent1_autonomous_reference_stress"]["observed_ci_at_registration"]["dedicated"]["status"] = "success"
    mutations.append(c)

    c = build_contract(exact_head=EXACT_HEAD)
    c["agent4_autonomous_reference_stress_audit"]["observed_ci_at_registration"]["repository_tests"]["conclusion"] = "success"
    mutations.append(c)

    for mutated in mutations:
        with pytest.raises(ValueError):
            validate_contract(_reseal(mutated))


def test_rejects_free_residual_defined_forcing_promotion() -> None:
    c = build_contract(exact_head=EXACT_HEAD)
    c["reference_stress_ingest_seam"]["residual_defined_forcing_forbidden"] = False
    with pytest.raises(ValueError):
        validate_contract(_reseal(c))
