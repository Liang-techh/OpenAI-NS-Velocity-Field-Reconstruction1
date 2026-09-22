from __future__ import annotations

import copy
import json

import pytest

from openai_ns_reconstruction import kokuno_a5_pulse_mj_physical_harmonic_rf31_routing as parent
from openai_ns_reconstruction import kokuno_a5_xi13_axissafe_rf39_routing as routing


def _rehash(reg: dict) -> dict:
    reg["registration_sha256"] = routing.registration_sha256(reg)
    return reg


def test_exact_parent_a5_contract_is_preserved() -> None:
    assert parent.SCHEMA == routing.PARENT_A5["schema"]
    assert parent.TASK == routing.PARENT_A5["task"]
    assert routing.PARENT_A5["pr"] == 1128
    assert routing.PARENT_A5["head"] == "f69f5842ea77d0837fc47ed9036b2e7cf799d101"
    assert routing.PARENT_A5["source_blob"] == "4a0d14b06301b8d2db6f5925676f29eea38e3c52"


def test_fresh_a1_xi13_leading_is_registered_without_global_promotion() -> None:
    reg = routing.build_registration()
    a1 = reg["frontiers"]["leading"]
    assert a1["pr"] == 1133
    assert a1["head"] == "3b4af71c4ab547bc11f9f3e320afa4b62c25b930"
    assert a1["source_blob"] == "4b8fad0781d74c934b6cc6fe7805a11d16807a92"
    assert a1["cartesian_end_compensation_composed"] is True
    assert a1["public_pulse_endpoint_xi"] == 13.0
    assert a1["save_load_materialized"] is True
    assert a1["source_exact_amplitude_materialized"] is False
    assert a1["source_exact_hidden_bump_recovered"] is False
    assert a1["terminal_global_leading_velocity_materialized"] is False


def test_a2_axis_safe_provider_is_not_laundered_into_self_contained_composite() -> None:
    reg = routing.build_registration()
    a2 = reg["frontiers"]["source_axis_safe_harmonic_sibling"]
    assert a2["pr"] == 1132
    assert a2["source_blob"] == "9067cb68dc49e08a3da7aaa429366f11d2073cdd"
    assert a2["provider_driven_velocity_xyzt_materialized"] is True
    assert a2["axis_zero_core_fail_closed"] is True
    assert a2["project_domain_source_input_provider_materialized"] is False
    assert a2["global_axis_safe_source_velocity_materialized"] is False
    assert a2["self_contained_velocity_xyzt_provider"] is False

    composite = reg["frontiers"]["latest_self_contained_project_composite"]
    assert composite["pr"] == 1117
    assert composite["public_pulse_endpoint_xi"] == 11.0
    assert composite["matching_xi13_agent2_composite"] is False
    assert composite["matching_agent4_audit_present"] is False


def test_a3_rf39_realization_remains_mechanics_only_without_real_rf30_backend() -> None:
    reg = routing.build_registration()
    a3 = reg["frontiers"]["rf34_rf39_mechanics_sibling"]
    assert a3["pr"] == 1134
    assert a3["head"] == "ecbc7e3702e366a5f355800dd662b11766273a41"
    assert a3["source_blob"] == "392975eafd3787c91e2dc82ad9b735567513271d"
    assert a3["rf34_rf39_compact_realization_executable"] is True
    assert a3["autonomous_pointwise_bump_realization"] is True
    assert a3["source_exact_pointwise_bump_recovered"] is False
    assert a3["current_i4_source_chart_backend_materialized"] is False
    assert a3["current_i4_rf30_defect_materialized"] is False
    assert a3["current_i4_rf34_rf39_correction_materialized"] is False
    assert a3["candidate_residual_evidence"] is False
    assert a3["cartesian_correction_velocity_materialized"] is False


def test_real_current_i4_force_and_a4_scoped_audit_stay_separate_from_newer_frontiers() -> None:
    reg = routing.build_registration()
    force = reg["frontiers"]["current_i4_real_radial_force"]
    assert force["pr"] == 1118
    assert force["radial_force_materialized"] is True
    assert force["complete_ns_correction_authorized"] is False

    audit = reg["frontiers"]["current_i4_radial_force_validator"]
    assert audit["pr"] == 1127
    assert audit["head"] == "add922982c2f3a4cf5a72d84910817e63f05c487"
    assert audit["source_blob"] == "c254b248e2703eaec1c4f757657ca40b9e211c7d"
    assert audit["audited_agent3_pr"] == 1118
    assert audit["scoped_gate_passed"] is None
    assert audit["scientifically_admitted"] is False
    assert audit["authorizes_complete_ns_correction"] is False
    assert audit["actions_status_at_registration"] == "queued"


def test_truth_boundary_marks_only_narrow_real_advances_true() -> None:
    truth = routing.build_registration()["truth_boundary"]
    for key in (
        "current_cartesian_leading_through_xi13_materialized",
        "current_cartesian_end_compensation_composed",
        "provider_driven_axis_safe_source_harmonic_velocity_materialized",
        "latest_self_contained_project_composite_is_agent2_1117",
        "rf34_rf39_compact_mechanics_materialized",
        "current_i4_radial_force_materialized",
        "agent4_matching_current_i4_radial_force_audit_present",
        "agent4_matching_current_i4_radial_force_audit_registered",
    ):
        assert truth[key] is True, key

    assert truth["agent4_matching_current_i4_radial_force_audit_scoped_gate_passed"] is None

    for key in (
        "source_exact_main_pulse_amplitude_materialized",
        "source_exact_hidden_end_bump_recovered",
        "terminal_global_leading_velocity_materialized",
        "matching_xi13_agent2_project_composite_materialized",
        "matching_xi13_agent4_audit_present",
        "project_domain_source_input_provider_materialized",
        "global_axis_safe_source_velocity_materialized",
        "source_harmonic_self_contained_velocity_xyzt_provider",
        "current_i4_source_chart_backend_materialized",
        "current_i4_rf30_defect_materialized",
        "current_i4_rf34_rf39_correction_materialized",
        "cartesian_correction_velocity_materialized",
        "agent4_matching_current_i4_radial_force_audit_admitted",
        "scoped_current_i4_force_authorized_as_complete_ns_correction_target",
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


def test_rf39_mechanics_cannot_be_rehashed_into_candidate_evidence() -> None:
    reg = routing.build_registration()
    promoted = copy.deepcopy(reg)
    promoted["frontiers"]["rf34_rf39_mechanics_sibling"]["current_i4_rf34_rf39_correction_materialized"] = True
    _rehash(promoted)
    with pytest.raises(ValueError, match="frontier identity"):
        routing.validate_registration(promoted)


def test_queued_a4_audit_cannot_be_rehashed_into_scoped_pass() -> None:
    reg = routing.build_registration()
    promoted = copy.deepcopy(reg)
    promoted["frontiers"]["current_i4_radial_force_validator"]["scoped_gate_passed"] = True
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
