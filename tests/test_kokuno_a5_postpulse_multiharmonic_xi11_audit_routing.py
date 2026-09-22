from __future__ import annotations

import copy
import json

import pytest

from openai_ns_reconstruction import kokuno_a5_xi13_axissafe_rf39_routing as parent
from openai_ns_reconstruction import kokuno_a5_postpulse_multiharmonic_xi11_audit_routing as routing


def _rehash(reg: dict) -> dict:
    reg["registration_sha256"] = routing.registration_sha256(reg)
    return reg


def test_exact_parent_a5_contract_is_preserved() -> None:
    assert parent.SCHEMA == routing.PARENT_A5["schema"]
    assert parent.TASK == routing.PARENT_A5["task"]
    assert routing.PARENT_A5["pr"] == 1135
    assert routing.PARENT_A5["head"] == "a6976c654926c3c96364c66bb18eac461688bde7"
    assert routing.PARENT_A5["source_blob"] == "a04996365234f33e81172a2d6a4f6bee8a7da993"


def test_a1_postpulse_eta_flattening_is_registered_without_global_promotion() -> None:
    reg = routing.build_registration()
    a1 = reg["frontiers"]["leading"]
    assert a1["pr"] == 1140
    assert a1["head"] == "c5442c11165d7f0976889b826bf58dce38a3d3dd"
    assert a1["source_blob"] == "32cbfc00406aa8d47904e2492ad8916730541d28"
    assert a1["autonomous_T_f"] == 100.0
    assert a1["source_exact_T_f_recovered"] is False
    assert a1["save_load_materialized"] is True
    assert a1["bounded_velocity_xyzt_materialized"] is True
    assert a1["terminal_global_leading_velocity_materialized"] is False


def test_a2_bounded_multiharmonic_family_remains_provider_driven() -> None:
    reg = routing.build_registration()
    a2 = reg["frontiers"]["bounded_multiharmonic_sibling"]
    assert a2["pr"] == 1141
    assert a2["head"] == "4241b0b4fcdc36e750798133101baaf803ba8c44"
    assert a2["source_blob"] == "06bc24722b4219419e91cd2d5b08cb1294665f24"
    assert a2["max_terms"] == 4
    assert a2["amplitude_l1_cap"] == 1.0
    assert a2["nontriviality_floor"] == 1e-12
    assert a2["complete_curl_family_materialized"] is True
    assert a2["provider_driven_velocity_xyzt_materialized"] is True
    assert a2["source_exact_amplitude_phase_recovered"] is False
    assert a2["project_domain_source_input_provider_materialized"] is False
    assert a2["self_contained_velocity_xyzt_provider"] is False


def test_xi11_self_contained_composite_now_has_matching_a4_audit_registered() -> None:
    reg = routing.build_registration()
    composite = reg["frontiers"]["latest_self_contained_project_composite"]
    audit = reg["frontiers"]["xi11_composite_validator"]
    assert composite["pr"] == 1117
    assert composite["public_pulse_endpoint_xi"] == 11.0
    assert composite["matching_agent4_audit_present"] is True
    assert composite["matching_agent4_audit_pr"] == 1136
    assert composite["matching_agent4_audit_scoped_gate_passed"] is None
    assert composite["matching_postpulse_agent2_composite"] is False

    assert audit["pr"] == 1136
    assert audit["head"] == "03269e7ba2053e0bc35813089577d274d9ddf8f9"
    assert audit["source_blob"] == "bd605bfb650b7eadbebb8d860e9720b434540ba3"
    assert audit["audited_agent2_pr"] == 1117
    assert audit["audited_agent1_pr"] == 1107
    assert audit["fd_steps"] == [0.02, 0.01, 0.005]
    assert audit["seed"] == 9173891
    assert audit["present"] is True
    assert audit["registered"] is True
    assert audit["scoped_gate_passed"] is None
    assert audit["scientifically_admitted"] is False
    assert audit["complete_ns_momentum_residual_assessed"] is False
    assert audit["pde_validated"] is False
    assert audit["actions_status_at_registration"] == "queued"


def test_a3_rf39_mechanics_and_older_i4_force_path_remain_separate() -> None:
    reg = routing.build_registration()
    a3 = reg["frontiers"]["rf34_rf39_mechanics_sibling"]
    assert a3["pr"] == 1134
    assert a3["current_i4_source_chart_backend_materialized"] is False
    assert a3["current_i4_rf30_defect_materialized"] is False
    assert a3["current_i4_rf34_rf39_correction_materialized"] is False
    assert a3["candidate_residual_evidence"] is False

    force = reg["frontiers"]["current_i4_real_radial_force"]
    validator = reg["frontiers"]["current_i4_radial_force_validator"]
    assert force["pr"] == 1118
    assert force["radial_force_materialized"] is True
    assert force["complete_ns_correction_authorized"] is False
    assert validator["pr"] == 1127
    assert validator["scientifically_admitted"] is False


def test_truth_boundary_marks_only_narrow_real_advances_true() -> None:
    truth = routing.build_registration()["truth_boundary"]
    for key in (
        "postpulse_eta_flattening_leading_materialized",
        "postpulse_eta_flattening_save_load_materialized",
        "bounded_multiharmonic_complete_curl_family_materialized",
        "bounded_multiharmonic_provider_driven_velocity_xyzt_materialized",
        "latest_self_contained_project_composite_is_agent2_1117",
        "agent4_matching_xi11_composite_audit_present",
        "agent4_matching_xi11_composite_audit_registered",
        "rf34_rf39_compact_mechanics_materialized",
        "current_i4_radial_force_materialized",
        "agent4_matching_current_i4_radial_force_audit_present",
    ):
        assert truth[key] is True, key

    assert truth["agent4_matching_xi11_composite_audit_scoped_gate_passed"] is None

    for key in (
        "source_exact_postpulse_T_f_recovered",
        "terminal_global_leading_velocity_materialized",
        "matching_postpulse_agent2_project_composite_materialized",
        "matching_postpulse_agent4_audit_present",
        "bounded_multiharmonic_self_contained_velocity_xyzt_provider",
        "source_exact_multiharmonic_amplitude_phase_recovered",
        "project_domain_source_input_provider_materialized",
        "agent4_matching_xi11_composite_audit_admitted",
        "current_i4_source_chart_backend_materialized",
        "current_i4_rf30_defect_materialized",
        "current_i4_rf34_rf39_correction_materialized",
        "cartesian_correction_velocity_materialized",
        "agent4_matching_current_i4_radial_force_audit_admitted",
        "matched_cartesian_pressure_gradient_materialized",
        "preregistered_restricted_forcing_materialized",
        "restricted_forcing_proved_not_residual_defined",
        "complete_candidate_api_ready",
        "complete_identity_bound_ns_defect_materialized",
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


def test_core_states_and_frozen_scientific_gates_do_not_move() -> None:
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


def test_truth_promotion_fails_closed_even_with_fresh_digest() -> None:
    reg = routing.build_registration()
    promoted = copy.deepcopy(reg)
    promoted["truth_boundary"]["pde_validated"] = True
    _rehash(promoted)
    with pytest.raises(ValueError, match="truth boundary"):
        routing.validate_registration(promoted)


def test_cross_identity_evidence_transfer_fails_closed() -> None:
    reg = routing.build_registration()
    mutated = copy.deepcopy(reg)
    mutated["identity_firewall"]["cross_identity_evidence_transfer_allowed"] = True
    _rehash(mutated)
    with pytest.raises(ValueError, match="identity firewall"):
        routing.validate_registration(mutated)


def test_xi11_audit_cannot_be_rehashed_into_scoped_pass_or_full_ns() -> None:
    reg = routing.build_registration()
    promoted = copy.deepcopy(reg)
    promoted["frontiers"]["xi11_composite_validator"]["scoped_gate_passed"] = True
    _rehash(promoted)
    with pytest.raises(ValueError, match="frontier identity"):
        routing.validate_registration(promoted)

    promoted = copy.deepcopy(reg)
    promoted["frontiers"]["xi11_composite_validator"]["complete_ns_momentum_residual_assessed"] = True
    _rehash(promoted)
    with pytest.raises(ValueError, match="frontier identity"):
        routing.validate_registration(promoted)


def test_provider_family_cannot_be_rehashed_into_self_contained_candidate() -> None:
    reg = routing.build_registration()
    promoted = copy.deepcopy(reg)
    promoted["frontiers"]["bounded_multiharmonic_sibling"]["self_contained_velocity_xyzt_provider"] = True
    _rehash(promoted)
    with pytest.raises(ValueError, match="frontier identity"):
        routing.validate_registration(promoted)


def test_gate_relaxation_fails_closed_even_with_fresh_digest() -> None:
    reg = routing.build_registration()
    relaxed = copy.deepcopy(reg)
    relaxed["frozen_gates"]["normalized_momentum_sampled_max"] = 2e-3
    _rehash(relaxed)
    with pytest.raises(ValueError, match="frozen project gates"):
        routing.validate_registration(relaxed)
