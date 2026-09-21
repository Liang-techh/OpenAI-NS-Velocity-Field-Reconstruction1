from __future__ import annotations

import copy

from openai_ns_reconstruction import kokuno_a5_current_pa16_profile_application_ingest_contract as c

HEAD = "0123456789abcdef0123456789abcdef01234567"


def _resign(payload: dict) -> dict:
    out = copy.deepcopy(payload)
    out.pop("contract_sha256", None)
    out["contract_sha256"] = c._sha256(out)
    return out


def test_frozen_contract_validates() -> None:
    payload = c.build_contract(HEAD)
    assert c.validate_contract(payload) == []
    assert payload["truth_boundary"]["current_pa16_source_profile_application_materialized"] is True
    assert payload["truth_boundary"]["five_moment_repair_closed"] is False
    assert payload["readiness"] == c.READINESS


def test_checksum_mutation_fails() -> None:
    payload = c.build_contract(HEAD)
    payload["pipeline_position"]["stage"] = "laundered"
    assert "contract_sha256_mismatch" in c.validate_contract(payload)


def test_parent_and_agent_heads_are_exact() -> None:
    payload = c.build_contract(HEAD)
    assert payload["parent_a5"]["head"] == "77d29056ac960fa09fdc899ce596dce68f842cd1"
    assert payload["agent1_joined_profile"]["head"] == "61b8730fba5623c5acfbd3ec42d54f7cc4f45748"
    assert payload["agent3_pa16_application"]["head"] == "670f7b71d72735d441dcf48f2391c271ee5c163c"
    assert payload["agent4_joined_profile_audit"]["head"] == "ac0d099777e9269719302850bcc0c567b139ca88"
    assert payload["agent2_sibling"]["head"] == "6d2fb1f701a34f783dca15a267ae2ce0734ba741"


def test_blob_pins_are_frozen() -> None:
    payload = c.build_contract(HEAD)
    assert payload["agent1_joined_profile"]["source_blob"] == "e38b9b79378aaf97a17226e41da75b477f26ae23"
    assert payload["agent3_pa16_application"]["source_blob"] == "cb3dcd25e0c55a7935daee1b77ffb1741320f83b"
    assert payload["agent4_joined_profile_audit"]["source_blob"] == "f3c02e1ee28e1b4f4b162be486953d1eca60a37b"
    assert payload["agent2_sibling"]["source_blob"] == "12df3bf6c949baeebf609b366f2973e7f515a157"


def test_a4_scope_cannot_be_laundered_into_repaired_profile_audit() -> None:
    payload = c.build_contract(HEAD)
    payload["truth_boundary"]["current_repaired_profile_independent_a4_audit_available"] = True
    payload = _resign(payload)
    errors = c.validate_contract(payload)
    assert "truth_boundary_drift" in errors
    assert "a4_scope_laundering" in errors


def test_a4_target_cannot_be_changed_to_agent3_coefficients() -> None:
    payload = c.build_contract(HEAD)
    payload["agent4_joined_profile_audit"]["audits_agent3_repair_coefficients"] = True
    payload = _resign(payload)
    errors = c.validate_contract(payload)
    assert "agent4_joined_profile_audit_drift" in errors
    assert "a4_scope_laundering" in errors


def test_source_profile_application_is_not_five_moment_closure() -> None:
    payload = c.build_contract(HEAD)
    payload["truth_boundary"]["five_moment_repair_closed"] = True
    payload = _resign(payload)
    errors = c.validate_contract(payload)
    assert "five_moment_closure_laundering" in errors


def test_source_profile_application_is_not_ns_correction() -> None:
    payload = c.build_contract(HEAD)
    payload["truth_boundary"]["source_profile_repair_is_complete_ns_correction"] = True
    payload = _resign(payload)
    assert "ns_correction_laundering" in c.validate_contract(payload)


def test_cartesian_global_leading_cannot_be_promoted() -> None:
    payload = c.build_contract(HEAD)
    payload["truth_boundary"]["corrected_global_cartesian_leading_velocity_materialized"] = True
    payload = _resign(payload)
    assert "forbidden_truth_promotion:corrected_global_cartesian_leading_velocity_materialized" in c.validate_contract(payload)


def test_pressure_forcing_and_complete_ns_remain_false() -> None:
    for key in (
        "matched_cartesian_pressure_materialized",
        "cartesian_matched_pressure_gradient_materialized",
        "restricted_forcing_materialized",
        "same_protocol_full_ns_residual_available",
    ):
        payload = c.build_contract(HEAD)
        payload["truth_boundary"][key] = True
        payload = _resign(payload)
        assert f"forbidden_truth_promotion:{key}" in c.validate_contract(payload)


def test_residual_defined_free_forcing_cannot_be_enabled() -> None:
    payload = c.build_contract(HEAD)
    payload["frozen_science"]["residual_defined_free_forcing_forbidden"] = False
    payload = _resign(payload)
    errors = c.validate_contract(payload)
    assert "frozen_science_drift" in errors
    assert "free_forcing_firewall_weakened" in errors


def test_final_gate_cannot_be_relaxed() -> None:
    payload = c.build_contract(HEAD)
    payload["final_gate"]["normalized_momentum_volume_l2"] = 2.0e-3
    payload = _resign(payload)
    errors = c.validate_contract(payload)
    assert "final_gate_drift" in errors
    assert "final_gate_changed" in errors


def test_st006_baseline_cannot_be_rewritten() -> None:
    payload = c.build_contract(HEAD)
    payload["st006_baseline"]["momentum_volume_l2"] = 1.0e-3
    payload = _resign(payload)
    errors = c.validate_contract(payload)
    assert "st006_baseline_drift" in errors
    assert "st006_baseline_changed" in errors


def test_queued_a4_evidence_cannot_be_admitted() -> None:
    payload = c.build_contract(HEAD)
    payload["truth_boundary"]["current_joined_profile_independent_a4_audit_admitted"] = True
    payload = _resign(payload)
    errors = c.validate_contract(payload)
    assert "queued_a4_evidence_cannot_be_admitted" in errors


def test_core_readiness_cannot_be_promoted() -> None:
    payload = c.build_contract(HEAD)
    payload["readiness"]["leading_ready"] = True
    payload = _resign(payload)
    errors = c.validate_contract(payload)
    assert "readiness_drift" in errors
    assert "readiness_promotion_forbidden" in errors


def test_coefficient_consistency_number_is_not_pde_evidence() -> None:
    payload = c.build_contract(HEAD)
    a3 = payload["agent3_pa16_application"]
    assert a3["coefficient_inf_abs_error"] == 1.8126505079421038e-13
    assert a3["coefficient_error_is_pde_residual"] is False
    assert a3["independently_certifies_post_application_source_moments"] is False


def test_invalid_exact_head_fails_closed() -> None:
    try:
        c.build_contract("not-a-sha")
    except ValueError:
        pass
    else:
        raise AssertionError("invalid exact head was accepted")
