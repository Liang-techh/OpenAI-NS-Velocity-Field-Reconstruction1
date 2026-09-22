from __future__ import annotations

import copy
import json

import pytest

from openai_ns_reconstruction import kokuno_a5_current_i4_mean_audit_pulse_routing as a5


def test_registration_binds_exact_current_i4_mean_and_audit_without_promotion() -> None:
    receipt = a5.build_registration()

    assert receipt["schema"] == a5.SCHEMA
    assert receipt["task"] == a5.TASK
    assert receipt["parent_a5"]["head"] == "666d5dad83dee8a9313c3a59bf102b10bb2917cc"

    a1 = receipt["agent1_current_main_pulse"]
    assert a1["head"] == "2a6c58ff6bf6d1ba6b875bca0e6abf9d1118e156"
    assert a1["current_cartesian_main_pulse_prefix_materialized"] is True
    assert a1["source_exact_amplitude_root_materialized"] is False
    assert a1["full_source_xi_11_current_cartesian_materialized"] is False
    assert a1["terminal_global_leading_velocity_materialized"] is False
    assert a1["matching_agent2_project_cartesian_composite_materialized"] is False

    a2 = receipt["agent2_source_shell_edge"]
    assert a2["head"] == "3f0988ff016f41126ec546681cba696c5fa8001d"
    assert a2["source_shell_edge_sqrt_zeta_materialized"] is True
    assert a2["complete_curl_localization_retains_support_gradient_terms"] is True
    assert a2["project_domain_coordinate_map_applied"] is False
    assert a2["self_contained_velocity_xyzt_provider"] is False
    assert a2["consumed_by_agent1_1100_identity"] is False

    a3 = receipt["agent3_current_i4_nonlinear_mean"]
    assert a3["head"] == "6c9a4e7b0dec4b3669fb1dc60297b18f10163e0b"
    assert a3["agent2_composite_head"] == "c40d8ddecd2971544a6e07dab093436b423cf326"
    assert a3["candidate_stage"] == "current-I4"
    assert a3["nonlinear_mean_attribution_materialized"] is True
    assert a3["source_i4_five_row_mean_correction_materialized"] is False
    assert a3["current_i4_compact_radial_stress_materialized"] is False
    assert a3["current_i4_radial_force_materialized"] is False
    assert a3["authorized_as_complete_ns_correction_target"] is False

    a4 = receipt["agent4_current_i4_nonlinear_mean_audit"]
    assert a4["head"] == "d3999c0e4e0e3483d4199f009e936060c42848c8"
    assert a4["audited_agent3_head"] == a3["head"]
    assert a4["agent2_composite_head"] == a3["agent2_composite_head"]
    assert a4["implementation_distinct"] is True
    assert a4["registered"] is True
    assert a4["scientific_admission"] is False
    assert a4["complete_ns_residual_evidence"] is False


def test_frontier_and_readiness_states_remain_fail_closed() -> None:
    receipt = a5.build_registration()
    truth = receipt["truth_boundary"]

    assert truth["leading_frontier_stage"] == "active-main-pulse-leading-only-finite-X-prefix"
    assert truth["latest_self_contained_composite_frontier_stage"] == "current-I4"
    assert truth["correction_frontier_stage"] == "current-I4-nonlinear-mean"
    assert truth["current_i4_nonlinear_mean_attribution_materialized"] is True
    assert truth["agent4_matching_current_i4_nonlinear_mean_audit_registered"] is True
    assert truth["agent4_matching_current_i4_nonlinear_mean_audit_admitted"] is False
    assert truth["current_i4_compact_radial_stress_materialized"] is False
    assert truth["current_i4_radial_force_materialized"] is False
    assert truth["older_current_i2_stress_force_evidence_transferred_to_current_i4"] is False
    assert truth["older_current_i2_stress_force_evidence_transferred_to_main_pulse"] is False
    assert truth["matched_cartesian_pressure_materialized"] is False
    assert truth["preregistered_restricted_forcing_materialized"] is False
    assert truth["complete_identity_bound_ns_defect_materialized"] is False
    assert truth["finite_correction_cycle_admitted"] is False
    assert truth["heldout_normalized_full_ns_residual_assessed"] is False
    assert truth["scientific_admission"] is False
    assert truth["pde_validated"] is False

    assert receipt["readiness"] == {
        "leading_ready": False,
        "oscillatory_ready": True,
        "correction_ready": False,
        "velocity_export_ready": False,
        "pde_validated": False,
    }


def test_fixed_project_gates_and_baseline_are_inherited_unchanged() -> None:
    receipt = a5.build_registration()
    gates = receipt["final_gates"]

    assert gates["normalized_momentum_sampled_max"] == 1.0e-3
    assert gates["normalized_momentum_volume_l2"] == 1.0e-3
    assert gates["normalized_divergence_sampled_max"] == 1.0e-5
    assert gates["normalized_divergence_volume_l2"] == 1.0e-5
    assert gates["canonical_quadrature"] == [24, 48, 96]
    assert gates["residual_defined_free_forcing_allowed"] is False

    baseline = receipt["st006_baseline"]
    assert baseline == a5.parent_a5.ST006_BASELINE


def test_deterministic_json_roundtrip(tmp_path) -> None:
    first = a5.build_registration()
    second = a5.build_registration()
    assert first == second
    assert first["digest"] == second["digest"]

    path = tmp_path / "routing.json"
    a5.save_registration(path)
    loaded = a5.load_registration(path)
    assert loaded == first
    assert json.loads(path.read_text(encoding="utf-8")) == first


@pytest.mark.parametrize(
    ("block", "key", "value"),
    [
        ("agent1_current_main_pulse", "source_exact_amplitude_root_materialized", True),
        ("agent1_current_main_pulse", "terminal_global_leading_velocity_materialized", True),
        ("agent2_source_shell_edge", "self_contained_velocity_xyzt_provider", True),
        ("agent2_source_shell_edge", "consumed_by_agent1_1100_identity", True),
        ("agent3_current_i4_nonlinear_mean", "current_i4_compact_radial_stress_materialized", True),
        ("agent3_current_i4_nonlinear_mean", "authorized_as_complete_ns_correction_target", True),
        ("agent4_current_i4_nonlinear_mean_audit", "scientific_admission", True),
        ("agent4_current_i4_nonlinear_mean_audit", "complete_ns_residual_evidence", True),
    ],
)
def test_evidence_promotion_mutations_fail_closed(block: str, key: str, value: object) -> None:
    receipt = a5.build_registration()
    mutated = copy.deepcopy(receipt)
    mutated[block][key] = value
    with pytest.raises(ValueError):
        a5.validate_registration(mutated)


@pytest.mark.parametrize(
    "key",
    [
        "agent4_matching_current_i4_nonlinear_mean_audit_admitted",
        "current_i4_compact_radial_stress_materialized",
        "current_i4_radial_force_materialized",
        "older_current_i2_stress_force_evidence_transferred_to_current_i4",
        "older_current_i2_stress_force_evidence_transferred_to_main_pulse",
        "matched_cartesian_pressure_materialized",
        "preregistered_restricted_forcing_materialized",
        "complete_identity_bound_ns_defect_materialized",
        "finite_correction_cycle_admitted",
        "heldout_normalized_full_ns_residual_assessed",
        "pde_validated",
    ],
)
def test_truth_boundary_promotions_fail_closed(key: str) -> None:
    receipt = a5.build_registration()
    mutated = copy.deepcopy(receipt)
    mutated["truth_boundary"][key] = True
    with pytest.raises(ValueError):
        a5.validate_registration(mutated)


def test_gate_relaxation_and_free_forcing_fail_closed() -> None:
    receipt = a5.build_registration()

    relaxed = copy.deepcopy(receipt)
    relaxed["final_gates"]["normalized_momentum_sampled_max"] = 2.0e-3
    with pytest.raises(ValueError):
        a5.validate_registration(relaxed)

    free_force = copy.deepcopy(receipt)
    free_force["final_gates"]["residual_defined_free_forcing_allowed"] = True
    with pytest.raises(ValueError):
        a5.validate_registration(free_force)


def test_exact_identity_mutation_fails_closed() -> None:
    receipt = a5.build_registration()
    mutated = copy.deepcopy(receipt)
    mutated["agent4_current_i4_nonlinear_mean_audit"]["audited_agent3_head"] = "0" * 40
    with pytest.raises(ValueError):
        a5.validate_registration(mutated)
