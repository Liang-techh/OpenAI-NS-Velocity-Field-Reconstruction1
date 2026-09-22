from __future__ import annotations

import copy

import pytest

from openai_ns_reconstruction import kokuno_a5_holdprefix_multiharmonic_rf49_audit_routing as a5


def _rehash(registration: dict) -> dict:
    registration["registration_sha256"] = a5.registration_sha256(registration)
    return registration


def test_registration_locks_fresh_frontiers_and_core_truth() -> None:
    reg = a5.build_registration()
    frontiers = reg["frontiers"]
    truth = reg["truth_boundary"]

    assert frontiers["leading"]["pr"] == 1148
    assert frontiers["leading"]["bounded_velocity_xyzt_materialized"] is True
    assert frontiers["leading"]["relative_swirl_bumps_materialized"] is False
    assert frontiers["postpulse_leading_validator"]["pr"] == 1143
    assert frontiers["postpulse_leading_validator"]["audited_agent1_pr"] == 1140
    assert frontiers["postpulse_leading_validator"]["does_not_audit_agent1_pr"] == 1148

    assert frontiers["bounded_multiharmonic_diagnostic"]["pr"] == 1149
    assert frontiers["bounded_multiharmonic_diagnostic"]["source_provider_self_contained"] is False
    assert frontiers["bounded_multiharmonic_diagnostic"]["complete_ns_residual_assessed"] is False

    assert frontiers["rf49_gain_firewall"]["pr"] == 1150
    assert frontiers["rf49_gain_firewall"]["gain_role"] == (
        "asymptotic-class-exponent-gain-not-raw-residual-contraction"
    )
    assert frontiers["rf49_gain_firewall"]["current_i4_source_chart_backend_materialized"] is False

    assert truth["source_relative_swirl_bumps_materialized"] is False
    assert truth["source_terminal_hold_after_eta_flattening_materialized"] is False
    assert truth["matching_latest_holdprefix_agent2_project_composite_materialized"] is False
    assert truth["matching_latest_holdprefix_agent4_audit_present"] is False
    assert truth["complete_identity_bound_ns_defect_materialized"] is False
    assert truth["heldout_complete_ns_residual_assessed"] is False
    assert truth["pde_validated"] is False

    assert reg["core_state"]["leading_ready"] is False
    assert reg["core_state"]["oscillatory_ready"] is True
    assert reg["core_state"]["correction_ready"] is False
    assert reg["core_state"]["velocity_export_ready"] is False
    assert reg["core_state"]["pde_validated"] is False


def test_registration_roundtrip_is_digest_bound(tmp_path) -> None:
    path = tmp_path / "registration.json"
    written = a5.write_registration(path)
    loaded = a5.load_registration(path)
    assert loaded == written
    assert loaded["registration_sha256"] == a5.registration_sha256(loaded)


def test_rejects_newer_a1_audit_evidence_transfer() -> None:
    reg = copy.deepcopy(a5.build_registration())
    reg["frontiers"]["postpulse_leading_validator"]["audited_agent1_pr"] = 1148
    _rehash(reg)
    with pytest.raises(ValueError, match="frontier identity drift"):
        a5.validate_registration(reg)


def test_rejects_rf49_gain_relabel_as_norm_contraction() -> None:
    reg = copy.deepcopy(a5.build_registration())
    reg["frontiers"]["rf49_gain_firewall"]["gain_role"] = "raw-residual-contraction-factor"
    _rehash(reg)
    with pytest.raises(ValueError, match="frontier identity drift"):
        a5.validate_registration(reg)


def test_rejects_truth_promotion_even_with_rehashed_receipt() -> None:
    reg = copy.deepcopy(a5.build_registration())
    reg["truth_boundary"]["pde_validated"] = True
    _rehash(reg)
    with pytest.raises(ValueError, match="truth boundary promotion/drift"):
        a5.validate_registration(reg)


def test_rejects_core_readiness_promotion() -> None:
    reg = copy.deepcopy(a5.build_registration())
    reg["core_state"]["correction_ready"] = True
    _rehash(reg)
    with pytest.raises(ValueError, match="core readiness promotion/drift"):
        a5.validate_registration(reg)


def test_rejects_project_gate_relaxation() -> None:
    reg = copy.deepcopy(a5.build_registration())
    reg["frozen_gates"]["normalized_momentum_sampled_max"] = 2.0e-3
    _rehash(reg)
    with pytest.raises(ValueError, match="frozen project gates changed"):
        a5.validate_registration(reg)


def test_identity_firewall_is_fail_closed() -> None:
    firewall = a5.build_identity_firewall()
    assert firewall["agent1_1148_hold_prefix_not_consumed_by_xi11_agent2_1117"] is True
    assert firewall["agent4_1143_audit_of_agent1_1140_not_evidence_for_agent1_1148"] is True
    assert firewall["agent2_1149_provider_diagnostic_not_self_contained_or_full_ns_evidence"] is True
    assert firewall["agent3_1150_source_exponent_gain_not_raw_residual_contraction"] is True
    assert firewall["agent3_1150_mechanics_not_current_i4_candidate_gain_evidence"] is True
    assert firewall["cross_identity_evidence_transfer_allowed"] is False
