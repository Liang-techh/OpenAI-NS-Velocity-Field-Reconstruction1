import copy

import pytest

from openai_ns_reconstruction.kokuno_a5_current_cartesian_leading_ingest_contract import (
    AGENT1_CARTESIAN_LEADING,
    AGENT4_CARTESIAN_DIVERGENCE_AUDIT,
    FINAL_GATE,
    PARENT_A5,
    READINESS,
    ST006_BASELINE,
    TRUTH_BOUNDARY,
    _sha256,
    build_contract,
    validate_contract,
)


HEAD = "1" * 40


def _rehash(payload):
    payload = copy.deepcopy(payload)
    payload.pop("contract_sha256", None)
    payload["contract_sha256"] = _sha256(payload)
    return payload


def test_contract_is_checksum_bound_and_valid():
    payload = build_contract(HEAD)
    assert validate_contract(payload) == []
    assert payload["parent_a5"]["head"] == PARENT_A5["head"]
    assert payload["agent1_cartesian_leading"]["head"] == AGENT1_CARTESIAN_LEADING["head"]
    assert payload["agent4_cartesian_divergence_audit"]["head"] == AGENT4_CARTESIAN_DIVERGENCE_AUDIT["head"]
    assert payload["st006_baseline"] == ST006_BASELINE
    assert payload["final_gate"] == FINAL_GATE
    assert payload["readiness"] == READINESS
    assert payload["truth_boundary"] == TRUTH_BOUNDARY
    assert build_contract(HEAD)["contract_sha256"] == payload["contract_sha256"]


def test_exact_head_must_be_lowercase_sha40():
    with pytest.raises(ValueError):
        build_contract("not-a-sha")
    with pytest.raises(ValueError):
        build_contract("A" * 40)


def test_rejects_a1_a4_lineage_drift():
    payload = build_contract(HEAD)
    payload["agent4_cartesian_divergence_audit"]["audited_agent1_head"] = "2" * 40
    errors = validate_contract(_rehash(payload))
    assert "a1_a4_lineage_mismatch" in errors


def test_rejects_velocity_beyond_xh_or_global_leading_laundering():
    payload = build_contract(HEAD)
    payload["truth_boundary"]["velocity_beyond_xh_materialized"] = True
    payload["truth_boundary"]["global_cartesian_leading_velocity_materialized"] = True
    errors = validate_contract(_rehash(payload))
    assert "global_leading_laundering" in errors


def test_rejects_scoped_divergence_as_canonical_or_momentum_evidence():
    payload = build_contract(HEAD)
    payload["agent4_cartesian_divergence_audit"]["canonical_whole_domain_divergence_l2"] = True
    payload["agent4_cartesian_divergence_audit"]["momentum_residual_assessed"] = True
    payload["agent4_cartesian_divergence_audit"]["complete_ns_residual_audit"] = True
    errors = validate_contract(_rehash(payload))
    assert "canonical_norm_laundering" in errors
    assert "a4_scope_laundering" in errors


def test_rejects_queued_evidence_promotion():
    payload = build_contract(HEAD)
    payload["truth_boundary"]["current_cartesian_leading_a4_divergence_audit_admitted"] = True
    payload["truth_boundary"]["current_cartesian_leading_divergence_scoped_assessed"] = True
    payload["truth_boundary"]["upstream_ci_admitted_as_pass"] = True
    payload["truth_boundary"]["scientific_admission"] = True
    errors = validate_contract(_rehash(payload))
    assert "queued_evidence_laundering" in errors
    assert "scientific_admission_laundering" in errors


def test_rejects_unmaterialized_cross_lane_composition():
    payload = build_contract(HEAD)
    payload["agent2_sibling"]["coupled_into_current_cartesian_leading"] = True
    payload["truth_boundary"]["leading_plus_oscillatory_cartesian_velocity_materialized"] = True
    payload["agent3_sibling"]["real_ns_correction_velocity_materialized"] = True
    payload["truth_boundary"]["real_agent3_ns_correction_velocity_materialized"] = True
    errors = validate_contract(_rehash(payload))
    assert "oscillatory_composition_laundering" in errors
    assert "correction_velocity_laundering" in errors


def test_rejects_pressure_forcing_and_complete_defect_laundering():
    payload = build_contract(HEAD)
    payload["truth_boundary"]["matched_cartesian_pressure_materialized"] = True
    payload["truth_boundary"]["restricted_forcing_materialized"] = True
    payload["agent3_sibling"]["discrepancy_from_complete_ns_defect"] = True
    payload["truth_boundary"]["discrepancy_from_complete_ns_defect"] = True
    payload["agent3_sibling"]["authorized_for_gain_gated_ns_stage"] = True
    payload["truth_boundary"]["authorized_for_gain_gated_ns_stage"] = True
    errors = validate_contract(_rehash(payload))
    assert "pressure_forcing_laundering" in errors
    assert "complete_defect_laundering" in errors
    assert "gain_stage_laundering" in errors


def test_rejects_threshold_free_forcing_and_quadrature_drift():
    payload = build_contract(HEAD)
    payload["frozen_science"]["residual_defined_free_forcing_forbidden"] = False
    payload["final_gate"]["normalized_momentum_volume_l2"] = 2.0e-3
    payload["final_gate"]["divergence_volume_l2"] = 2.0e-5
    payload["final_gate"]["canonical_volume_quadrature_ladder"] = [24, 48]
    errors = validate_contract(_rehash(payload))
    assert "free_forcing_firewall_weakened" in errors
    assert "momentum_gate_drift" in errors
    assert "divergence_gate_drift" in errors
    assert "canonical_quadrature_drift" in errors


def test_rejects_readiness_or_pde_promotion():
    for key in ("leading_ready", "correction_ready", "velocity_export_ready", "pde_validated"):
        payload = build_contract(HEAD)
        payload["readiness"][key] = True
        assert "readiness_premature_promotion" in validate_contract(_rehash(payload))


def test_rejects_st006_or_paper_exact_drift():
    payload = build_contract(HEAD)
    payload["st006_baseline"]["momentum_volume_l2"] = 0.01
    payload["truth_boundary"]["paper_exact"] = True
    errors = validate_contract(_rehash(payload))
    assert "st006_baseline_drift" in errors
    assert "paper_exact_laundering" in errors
