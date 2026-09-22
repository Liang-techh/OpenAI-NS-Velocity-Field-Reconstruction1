from __future__ import annotations

import copy
import json

import pytest

from openai_ns_reconstruction import kokuno_a5_pulse_mj_physical_harmonic_rf31_routing as routing
from openai_ns_reconstruction import kokuno_a5_xi11_composite_i4_force_routing as parent


def _rehash(reg: dict) -> dict:
    reg["registration_sha256"] = routing.registration_sha256(reg)
    return reg


def test_parent_a5_exact_contract_is_preserved() -> None:
    assert parent.SCHEMA == routing.PARENT_A5["schema"]
    assert parent.TASK == routing.PARENT_A5["task"]
    assert routing.PARENT_A5["pr"] == 1120
    assert routing.PARENT_A5["head"] == "337ad149f6227fd054e496f4f38e48216b2ff2f4"
    assert routing.PARENT_A5["source_blob"] == "9559ea30dcc6bea7c4d8179b555e4884108c240a"


def test_new_upstream_frontiers_are_exact_and_non_promoting() -> None:
    reg = routing.build_registration()
    frontiers = reg["frontiers"]

    a1 = frontiers["leading_end_bookkeeping"]
    assert a1["pr"] == 1124
    assert a1["current_lineage_J_entry_materialized"] is True
    assert a1["current_J_assumed_zero"] is False
    assert a1["cartesian_end_compensation_composed"] is False
    assert a1["terminal_global_leading_velocity_materialized"] is False

    a2_source = frontiers["source_oscillatory_physicalization_sibling"]
    assert a2_source["pr"] == 1125
    assert a2_source["source_physical_Q_scaling_applied"] is True
    assert a2_source["self_contained_velocity_xyzt_provider"] is False

    composite = frontiers["latest_self_contained_project_composite"]
    assert composite["pr"] == 1117
    assert composite["head"] == "27741d9c0a27262f7fabf61eebaa2fbd507e9f03"
    assert composite["matching_agent4_audit_present"] is False

    a3 = frontiers["rf30_rf31_typed_executor_sibling"]
    assert a3["pr"] == 1126
    assert a3["rf30_rf31_formula_executor_materialized"] is True
    assert a3["current_i4_source_fixed_q_backend_materialized"] is False
    assert a3["current_i4_candidate_defect_evidence"] is False
    assert a3["correction_coefficients_solved"] is False
    assert a3["correction_applied"] is False
    assert a3["nonlinear_remainder_recomputed"] is False


def test_radial_force_a4_claim_is_pending_not_evidence() -> None:
    reg = routing.build_registration()
    pending = reg["frontiers"]["current_i4_radial_force_validator_pending"]
    assert pending["task"] == "K4-VAL-116"
    assert pending["status"] == "claimed-pending-delivery"
    assert pending["exact_base_pr"] == 1118
    assert pending["exact_base_head"] == "922f7aa10460ded44af212d313eceff33a2ac647"
    assert pending["delivered_pr"] is None
    assert pending["delivered_head"] is None
    assert pending["present"] is False
    assert pending["registered"] is False
    assert pending["scientifically_admitted"] is False

    stress = reg["frontiers"]["current_i4_stress_validator"]
    assert stress["pr"] == 1119
    assert stress["present"] is True
    assert stress["registered"] is True
    assert stress["scientifically_admitted"] is False
    assert stress["audits_radial_force"] is False


def test_truth_boundary_keeps_only_real_narrow_advances_true() -> None:
    truth = routing.build_registration()["truth_boundary"]
    assert truth["current_lineage_J_entry_materialized"] is True
    assert truth["source_physicalized_localized_harmonic_materialized"] is True
    assert truth["current_i4_radial_force_materialized"] is True
    assert truth["rf30_rf31_formula_executor_materialized"] is True
    assert truth["agent4_matching_current_i4_radial_stress_audit_present"] is True
    assert truth["agent4_current_i4_radial_force_audit_claimed"] is True

    for key in (
        "current_cartesian_end_compensation_composed",
        "source_exact_main_pulse_amplitude_materialized",
        "terminal_global_leading_velocity_materialized",
        "matching_agent4_xi11_composite_audit_present",
        "agent4_matching_current_i4_radial_force_audit_present",
        "agent4_matching_current_i4_radial_force_audit_registered",
        "agent4_matching_current_i4_radial_force_audit_admitted",
        "current_i4_source_fixed_q_backend_materialized",
        "current_i4_rf30_defect_materialized",
        "current_i4_rf31_five_row_system_materialized",
        "correction_coefficients_solved",
        "correction_applied",
        "rf44_rf49_nonlinear_remainder_recomputed",
        "scoped_current_i4_force_authorized_as_complete_ns_correction_target",
        "matched_cartesian_pressure_gradient_materialized",
        "preregistered_restricted_forcing_materialized",
        "restricted_forcing_proved_not_residual_defined",
        "complete_identity_bound_ns_defect_materialized",
        "cartesian_correction_velocity_materialized",
        "finite_correction_cycle_run",
        "heldout_complete_ns_residual_assessed",
        "canonical_whole_domain_admission_run",
        "same_protocol_comparable_to_st006",
        "scientific_admission",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
        "pde_validated",
    ):
        assert truth[key] is False, key


def test_core_states_and_fixed_scientific_gates_do_not_move() -> None:
    reg = routing.build_registration()
    assert reg["core_state"] == {
        "leading_ready": False,
        "oscillatory_ready": True,
        "correction_ready": False,
        "velocity_export_ready": False,
        "pde_validated": False,
    }
    gates = reg["frozen_gates"]
    assert gates["normalized_momentum_sampled_max"] == 1e-3
    assert gates["normalized_momentum_volume_l2"] == 1e-3
    assert gates["normalized_divergence_sampled_max"] == 1e-5
    assert gates["normalized_divergence_volume_l2"] == 1e-5
    assert gates["canonical_quadrature"] == [24, 48, 96]
    assert gates["residual_defined_free_forcing_allowed"] is False
    assert gates["st006_momentum_sampled_max"] == 0.1082289305112118
    assert gates["st006_momentum_volume_l2"] == 0.10758432876230622


def test_registration_is_deterministic_and_roundtrips(tmp_path) -> None:
    first = routing.build_registration()
    second = routing.build_registration()
    assert first == second
    assert len(first["registration_sha256"]) == 64

    path = tmp_path / "routing.json"
    written = routing.write_registration(path)
    loaded = routing.load_registration(path)
    assert loaded == written == first
    assert json.loads(path.read_text()) == first


def test_digest_and_truth_promotion_mutations_fail_closed() -> None:
    reg = routing.build_registration()
    corrupt = copy.deepcopy(reg)
    corrupt["truth_boundary"]["pde_validated"] = True
    with pytest.raises(ValueError, match="digest"):
        routing.validate_registration(corrupt)

    promoted = copy.deepcopy(reg)
    promoted["truth_boundary"]["pde_validated"] = True
    _rehash(promoted)
    with pytest.raises(ValueError, match="truth boundary"):
        routing.validate_registration(promoted)


def test_pending_a4_claim_cannot_be_rehashed_into_delivered_evidence() -> None:
    reg = routing.build_registration()
    promoted = copy.deepcopy(reg)
    pending = promoted["frontiers"]["current_i4_radial_force_validator_pending"]
    pending["delivered_pr"] = 9999
    pending["delivered_head"] = "a" * 40
    pending["present"] = True
    pending["registered"] = True
    _rehash(promoted)
    with pytest.raises(ValueError, match="pending A4"):
        routing.validate_registration(promoted)


def test_identity_evidence_laundering_fails_closed() -> None:
    reg = routing.build_registration()
    mutated = copy.deepcopy(reg)
    mutated["identity_firewall"]["cross_identity_evidence_transfer_allowed"] = True
    _rehash(mutated)
    with pytest.raises(ValueError, match="identity firewall"):
        routing.validate_registration(mutated)


def test_gate_relaxation_fails_closed_even_with_fresh_digest() -> None:
    reg = routing.build_registration()
    relaxed = copy.deepcopy(reg)
    relaxed["frozen_gates"]["normalized_momentum_sampled_max"] = 2e-3
    _rehash(relaxed)
    with pytest.raises(ValueError, match="frozen project gates"):
        routing.validate_registration(relaxed)
