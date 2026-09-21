from __future__ import annotations

import copy

import pytest

from openai_ns_reconstruction.kokuno_a5_current_tsh_certificate_ingest_contract import (
    AGENT1_CURRENT_TSH,
    AGENT2_SIBLING_STATUS,
    AGENT3_CURRENT_XI_HANDOFF,
    AGENT4_CURRENT_TSH_AUDIT,
    FINAL_GATE,
    PARENT_A5,
    READINESS,
    ST006_BASELINE,
    TRUTH_BOUNDARY,
    _sha256,
    build_contract,
    validate_contract,
)


DUMMY_HEAD = "a" * 40


def _reseal(payload: dict) -> dict:
    payload = copy.deepcopy(payload)
    payload.pop("contract_sha256", None)
    payload["contract_sha256"] = _sha256(payload)
    return payload


def test_current_tsh_ingest_contract_pins_exact_lineage_and_keeps_states_closed() -> None:
    payload = build_contract(DUMMY_HEAD)
    assert validate_contract(payload) == []
    assert payload["parent_a5"] == PARENT_A5
    assert payload["agent1_current_tsh"] == AGENT1_CURRENT_TSH
    assert payload["agent4_current_tsh_audit"] == AGENT4_CURRENT_TSH_AUDIT
    assert payload["agent3_current_xi_handoff"] == AGENT3_CURRENT_XI_HANDOFF
    assert payload["agent2_sibling_status"] == AGENT2_SIBLING_STATUS
    assert payload["final_gate"] == FINAL_GATE
    assert payload["st006_baseline"] == ST006_BASELINE
    assert payload["readiness"] == READINESS
    assert payload["truth_boundary"] == TRUTH_BOUNDARY

    assert payload["agent1_current_tsh"]["candidate_side_numerical_T_sh_instantiated"] is True
    assert payload["agent1_current_tsh"]["source_T_sh_lower_bound_verified"] is False
    assert payload["agent4_current_tsh_audit"]["audit_protocol_registered"] is True
    assert payload["agent4_current_tsh_audit"]["audit_admitted"] is False
    assert payload["agent4_current_tsh_audit"]["candidate_side_route_feasibility"] == (
        "unresolved_pending_exact_head_receipt"
    )
    assert payload["agent3_current_xi_handoff"]["pa16_repair_applied"] is False
    assert payload["agent3_current_xi_handoff"]["authorized_for_gain_gated_ns_stage"] is False
    assert payload["agent2_sibling_status"]["creates_correction_velocity"] is False

    assert payload["readiness"] == {
        "leading_ready": False,
        "oscillatory_ready": True,
        "correction_ready": False,
        "velocity_export_ready": False,
        "pde_validated": False,
    }


def test_current_tsh_ingest_does_not_launder_queued_ci_into_admission() -> None:
    payload = build_contract(DUMMY_HEAD)
    for section in (
        "parent_a5",
        "agent1_current_tsh",
        "agent4_current_tsh_audit",
        "agent3_current_xi_handoff",
        "agent2_sibling_status",
    ):
        for run in payload[section]["observed_ci"].values():
            assert run["status"] == "queued"
            assert run["conclusion"] is None
    assert payload["truth_boundary"]["upstream_ci_admitted_as_pass"] is False
    assert payload["truth_boundary"]["agent4_current_tsh_independent_audit_admitted"] is False
    assert payload["truth_boundary"]["scientific_admission"] is False


def test_current_tsh_ingest_preserves_frozen_pde_gate_and_baseline() -> None:
    payload = build_contract(DUMMY_HEAD)
    assert payload["final_gate"]["normalized_momentum_sampled_max"] == pytest.approx(1.0e-3)
    assert payload["final_gate"]["normalized_momentum_volume_l2"] == pytest.approx(1.0e-3)
    assert payload["final_gate"]["divergence_sampled_max"] == pytest.approx(1.0e-5)
    assert payload["final_gate"]["divergence_volume_l2"] == pytest.approx(1.0e-5)
    assert payload["frozen_science"]["residual_defined_free_forcing_forbidden"] is True
    assert payload["st006_baseline"]["momentum_sampled_max"] == pytest.approx(
        0.1082289305112118
    )
    assert payload["st006_baseline"]["momentum_volume_l2"] == pytest.approx(
        0.10758432876230622
    )


@pytest.mark.parametrize(
    ("section", "key", "value", "expected_error"),
    [
        ("truth_boundary", "source_T_sh_lower_bound_verified", True, "truth_boundary_drift"),
        ("truth_boundary", "pa16_repair_applied", True, "truth_boundary_drift"),
        ("truth_boundary", "actual_inner_to_outer_global_join_materialized", True, "truth_boundary_drift"),
        ("truth_boundary", "matched_cartesian_pressure_materialized", True, "truth_boundary_drift"),
        ("truth_boundary", "restricted_forcing_materialized", True, "truth_boundary_drift"),
        ("truth_boundary", "real_agent3_correction_velocity_materialized", True, "truth_boundary_drift"),
        ("truth_boundary", "scientific_admission", True, "truth_boundary_drift"),
        ("readiness", "leading_ready", True, "readiness_drift"),
        ("readiness", "correction_ready", True, "readiness_drift"),
        ("readiness", "velocity_export_ready", True, "readiness_drift"),
        ("readiness", "pde_validated", True, "readiness_drift"),
        ("final_gate", "normalized_momentum_volume_l2", 1.0e-2, "final_gate_drift"),
        ("final_gate", "divergence_volume_l2", 1.0e-4, "final_gate_drift"),
        ("frozen_science", "residual_defined_free_forcing_forbidden", False, "frozen_science_drift"),
    ],
)
def test_current_tsh_ingest_rejects_truth_and_gate_laundering(
    section: str, key: str, value: object, expected_error: str
) -> None:
    payload = build_contract(DUMMY_HEAD)
    payload[section][key] = value
    payload = _reseal(payload)
    errors = validate_contract(payload)
    assert expected_error in errors


def test_current_tsh_ingest_rejects_upstream_identity_drift_even_with_new_checksum() -> None:
    payload = build_contract(DUMMY_HEAD)
    payload["agent1_current_tsh"]["head"] = "b" * 40
    payload = _reseal(payload)
    assert "agent1_current_tsh_drift" in validate_contract(payload)

    payload = build_contract(DUMMY_HEAD)
    payload["agent4_current_tsh_audit"]["protocol"]["seed"] += 1
    payload = _reseal(payload)
    assert "agent4_current_tsh_audit_drift" in validate_contract(payload)


def test_current_tsh_ingest_rejects_premature_a4_admission() -> None:
    payload = build_contract(DUMMY_HEAD)
    payload["agent4_current_tsh_audit"]["audit_admitted"] = True
    payload["truth_boundary"]["agent4_current_tsh_independent_audit_admitted"] = True
    payload["truth_boundary"]["upstream_ci_admitted_as_pass"] = True
    payload = _reseal(payload)
    errors = validate_contract(payload)
    assert "agent4_current_tsh_audit_drift" in errors
    assert "truth_boundary_drift" in errors


def test_current_tsh_ingest_checksum_is_fail_closed() -> None:
    payload = build_contract(DUMMY_HEAD)
    payload["pipeline_position"]["next_shortest_blocker"] = "pretend complete"
    assert "contract_sha256_mismatch" in validate_contract(payload)


def test_current_tsh_ingest_requires_real_sha_shape() -> None:
    with pytest.raises(ValueError, match="lowercase 40-hex SHA"):
        build_contract("not-a-sha")
