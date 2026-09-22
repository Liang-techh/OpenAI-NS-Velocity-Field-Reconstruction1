from __future__ import annotations

import copy
import json

import pytest

from openai_ns_reconstruction import kokuno_a5_xi11_pulse_i4_stress_routing as a5


def test_registration_harvests_exact_a1_a2_a3_a4_without_identity_laundering() -> None:
    receipt = a5.build_registration()
    assert receipt["schema"] == a5.SCHEMA
    assert receipt["task"] == a5.TASK
    assert receipt["parent_a5"]["head"] == "68535271cdbb79cead9431109b5311d6047cfd0e"

    a1 = receipt["agent1_xi11_main_pulse"]
    assert a1["head"] == "45da043dd2b4cd067f005a72c7e21fd0f2bcf309"
    assert a1["full_source_xi_11_current_cartesian_materialized"] is True
    assert a1["source_exact_amplitude_root_materialized"] is False
    assert a1["source_pulse_end_mj_corrections_materialized"] is False
    assert a1["terminal_global_leading_velocity_materialized"] is False

    a2 = receipt["agent2_main_pulse_prefix_composite"]
    assert a2["head"] == "9f221dd57ef4b5e2d2e80e24c8cf531e991a5387"
    assert a2["consumed_agent1_head"] == a1["parent_head"]
    assert a2["newer_agent1_head"] == a1["head"]
    assert a2["newer_agent1_1107_consumed"] is False
    assert a2["self_contained_project_cartesian_velocity"] is True
    assert a2["full_source_xi_11_current_cartesian_materialized"] is False

    a3 = receipt["agent3_current_i4_radial_stress"]
    assert a3["head"] == "49590ff311fef1dec4bd850b3013fd989485c87a"
    assert a3["parent_head"] == a5.parent_a5.AGENT3_CURRENT_I4_NONLINEAR_MEAN["head"]
    assert a3["current_i4_compact_radial_stress_materialized"] is True
    assert a3["radial_force_partial_z_sigma_1_materialized"] is False
    assert a3["matching_agent4_stress_audit_present"] is False
    assert a3["authorized_as_complete_ns_correction_target"] is False

    a4 = receipt["agent4_main_pulse_prefix_divergence_audit"]
    assert a4["head"] == "8a9fb9ed1a67ec5ad7d9d856be71ee23467654bf"
    assert a4["audited_agent2_head"] == a2["head"]
    assert a4["audited_agent1_head"] == a2["consumed_agent1_head"]
    assert a4["newer_agent1_1107_consumed"] is False
    assert a4["implementation_distinct"] is True
    assert a4["registered"] is True
    assert a4["scientific_admission"] is False
    assert a4["complete_ns_residual_evidence"] is False


def test_frontiers_and_core_readiness_remain_fail_closed() -> None:
    receipt = a5.build_registration()
    truth = receipt["truth_boundary"]

    assert truth["leading_frontier_stage"] == "current-main-pulse-leading-only-logX-through-xi11"
    assert truth["latest_self_contained_composite_frontier_stage"] == "current-main-pulse-finite-X-prefix"
    assert truth["correction_frontier_stage"] == "current-I4-compact-radial-stress"
    assert truth["agent1_xi11_logx_main_pulse_leading_materialized"] is True
    assert truth["main_pulse_prefix_agent2_project_composite_materialized"] is True
    assert truth["matching_xi11_agent2_project_composite_materialized"] is False
    assert truth["agent4_matching_main_pulse_prefix_divergence_audit_registered"] is True
    assert truth["agent4_matching_main_pulse_prefix_divergence_audit_admitted"] is False
    assert truth["current_i4_compact_radial_stress_materialized"] is True
    assert truth["current_i4_radial_force_materialized"] is False
    assert truth["agent4_matching_current_i4_radial_stress_audit_present"] is False
    assert truth["complete_identity_bound_ns_defect_materialized"] is False
    assert truth["pde_validated"] is False

    assert receipt["readiness"] == {
        "leading_ready": False,
        "oscillatory_ready": True,
        "correction_ready": False,
        "velocity_export_ready": False,
        "pde_validated": False,
    }


def test_fixed_project_gates_and_st006_are_inherited_unchanged() -> None:
    receipt = a5.build_registration()
    gates = receipt["final_gates"]
    assert gates["normalized_momentum_sampled_max"] == 1.0e-3
    assert gates["normalized_momentum_volume_l2"] == 1.0e-3
    assert gates["normalized_divergence_sampled_max"] == 1.0e-5
    assert gates["normalized_divergence_volume_l2"] == 1.0e-5
    assert gates["canonical_quadrature"] == [24, 48, 96]
    assert gates["residual_defined_free_forcing_allowed"] is False
    assert receipt["st006_baseline"] == a5.parent_a5.ST006_BASELINE


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
        ("agent1_xi11_main_pulse", "source_exact_amplitude_root_materialized", True),
        ("agent1_xi11_main_pulse", "terminal_global_leading_velocity_materialized", True),
        ("agent1_xi11_main_pulse", "matching_agent2_composite_consuming_this_exact_identity_materialized", True),
        ("agent2_main_pulse_prefix_composite", "newer_agent1_1107_consumed", True),
        ("agent2_main_pulse_prefix_composite", "full_source_xi_11_current_cartesian_materialized", True),
        ("agent3_current_i4_radial_stress", "radial_force_partial_z_sigma_1_materialized", True),
        ("agent3_current_i4_radial_stress", "matching_agent4_stress_audit_present", True),
        ("agent3_current_i4_radial_stress", "authorized_as_complete_ns_correction_target", True),
        ("agent4_main_pulse_prefix_divergence_audit", "scientific_admission", True),
        ("agent4_main_pulse_prefix_divergence_audit", "complete_ns_residual_evidence", True),
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
        "source_exact_main_pulse_amplitude_materialized",
        "source_pulse_end_mj_corrections_materialized",
        "terminal_global_leading_completion_materialized",
        "main_pulse_prefix_agent2_composite_consumes_newer_agent1_1107",
        "matching_xi11_agent2_project_composite_materialized",
        "agent4_matching_main_pulse_prefix_divergence_audit_admitted",
        "agent4_matching_main_pulse_prefix_divergence_audit_is_complete_ns_evidence",
        "agent4_main_pulse_prefix_audit_transferred_to_xi11_identity",
        "current_i4_radial_force_materialized",
        "agent4_matching_current_i4_radial_stress_audit_present",
        "source_i4_five_row_mean_correction_materialized",
        "current_i4_stress_authorized_as_complete_ns_correction_target",
        "current_i4_stress_evidence_transferred_to_main_pulse_prefix",
        "current_i4_stress_evidence_transferred_to_xi11_pulse",
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
    mutated["agent4_main_pulse_prefix_divergence_audit"]["audited_agent2_head"] = "0" * 40
    with pytest.raises(ValueError):
        a5.validate_registration(mutated)
