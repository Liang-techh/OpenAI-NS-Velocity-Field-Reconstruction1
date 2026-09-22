from __future__ import annotations

import json
from pathlib import Path

import pytest

from openai_ns_reconstruction import (
    kokuno_a5_scalar_refutation_runtime_rebind_routing as r,
)


def test_exact_parent_and_upstream_heads_are_frozen() -> None:
    assert r.PARENT_A5_EXACT_HEAD == "8ed2641571bbce4340252e1809550959f57acbec"
    assert r.AGENT1_EXACT_HEAD == "7e987ff0a3bea743fe3f0fe9d02d5b3c8c593ea3"
    assert r.AGENT4_EXACT_HEAD == "d38aa3dacc205deacf785308b616c24ffb4eb135"
    assert r.SELECTED_AGENT2_EXACT_HEAD == "3afc18d93334b6d371856cce28f90c6bff40bbc0"
    assert r.DUPLICATE_AGENT3_EXACT_HEAD == "20120093abbd981e2d2cb5d95093cea1eebb28ef"


def test_common_scalar_witness_is_one_way_and_result_unresolved() -> None:
    a1 = r.leading_route()["agent1"]
    assert a1["finite_unique_root_witness_surface_present"] is True
    assert a1["positive_witness_would_refute_one_common_scalar_bridge"] is True
    assert a1["negative_witness_would_establish_common_scalar_bridge"] is False
    assert a1["actual_refutation_result_ingested"] is False
    assert a1["common_scalar_bridge_refuted"] is None
    assert a1["full_eta_common_scalar_bridge_established"] is False
    assert a1["eta_dependent_cartesian_bridge_materialized"] is False


def test_matching_a4_surface_is_present_but_not_ingested() -> None:
    a4 = r.leading_route()["agent4"]
    assert a4["audits_agent1_exact_head"] == r.AGENT1_EXACT_HEAD
    assert a4["matching_independent_audit_surface_present"] is True
    assert a4["scientific_result_ingested"] is False
    assert a4["scoped_gate_passed"] is None
    assert a4["audit_is_complete_ns_validation"] is False
    assert a4["pde_validated"] is False


def test_representation_firewall_rejects_finite_to_continuum_promotion() -> None:
    firewall = r.leading_route()["representation_firewall"]
    assert firewall
    assert not any(firewall.values())


def test_single_downstream_correction_route_is_selected() -> None:
    route = r.correction_route()
    selected = route["selected_route"]
    duplicate = route["duplicate_nonselected_route"]
    topo = route["topology_firewall"]
    assert selected["parent_pr"] == 1234
    assert selected["agent2_pr"] == 1243
    assert selected["selected_as_single_downstream_unified_route"] is True
    assert duplicate["agent3_pr"] == 1236
    assert duplicate["selected_for_downstream"] is False
    assert duplicate["evidence_unioned_into_selected_route"] is False
    assert topo["selected_route_is_agent2_1234_to_1243"] is True
    assert topo["duplicate_agent3_1236_may_be_silently_unioned"] is False
    assert topo["duplicate_capability_count"] == 1


def test_runtime_rebind_surface_does_not_claim_execution_before_ci() -> None:
    selected = r.correction_route()["selected_route"]
    assert selected["exact_a2_1080_external_snapshot_declared"] is True
    assert selected["exact_a1_1079_external_snapshot_declared"] is True
    assert selected["runtime_rebind_contract_surface_present"] is True
    assert selected["exact_head_ci_resolved"] is False
    assert selected["runtime_rebind_result_ingested"] is False
    assert selected["exact_backend_rebind_promoted"] is False
    assert selected["rf30_to_rf39_candidate_execution_recorded"] is False
    assert selected["rf44_post_update_state_materialized"] is False
    assert selected["cartesian_delta_u_materialized"] is False
    assert selected["finite_correction_cycle_run"] is False


def test_latest_self_contained_candidate_does_not_move() -> None:
    assert r.latest_self_contained_candidate() == {
        "agent2_pr": 1117,
        "agent1_pr": 1107,
        "stage": "xi=11",
        "self_contained_leading_plus_oscillatory_velocity_xyzt": True,
        "superseded_by_agent1_1242": False,
        "superseded_by_agent2_1243": False,
        "superseded_by_agent4_1244": False,
    }


def test_core_states_remain_fail_closed() -> None:
    assert r.core_states() == {
        "leading_ready": False,
        "oscillatory_ready": True,
        "correction_ready": False,
        "velocity_export_ready": False,
        "pde_validated": False,
    }


def test_validation_protocol_is_unchanged() -> None:
    p = r.fixed_validation_protocol()
    assert p["st006_momentum_sampled_max"] == pytest.approx(0.1082289305112118)
    assert p["st006_momentum_volume_l2"] == pytest.approx(0.10758432876230622)
    assert p["normalized_momentum_sampled_max_gate"] == pytest.approx(1e-3)
    assert p["normalized_momentum_volume_l2_gate"] == pytest.approx(1e-3)
    assert p["normalized_divergence_sampled_max_gate"] == pytest.approx(1e-5)
    assert p["normalized_divergence_volume_l2_gate"] == pytest.approx(1e-5)
    assert p["canonical_quadrature"] == [24, 48, 96]
    assert p["residual_defined_free_forcing_forbidden"] is True
    assert p["posthoc_threshold_relaxation_forbidden"] is True
    assert p["kokuno_replay_is_independent_final_validation"] is False
    assert p["green_software_ci_is_pde_validation"] is False
    assert p["queued_or_running_ci_is_pass"] is False


def test_closed_and_open_blockers_are_distinguished() -> None:
    b = r.remaining_blockers()
    assert b["a1_common_scalar_refutation_surface"] is False
    assert b["matching_a4_refutation_audit_surface"] is False
    assert b["single_correction_route_selected"] is False
    assert b["exact_runtime_rebind_contract_surface"] is False
    assert b["a4_refutation_audit_result_resolved"] is True
    assert b["full_eta_bridge_geometry"] is True
    assert b["exact_runtime_rebind_ci_and_result_resolved"] is True
    assert b["real_candidate_rf30_to_rf39_execution"] is True
    assert b["complete_identity_bound_ns_defect"] is True
    assert b["agent4_heldout_canonical_complete_ns_gate"] is True


def test_shortest_closure_preserves_conditional_bridge_routing() -> None:
    closure = "\n".join(r.shortest_closure())
    assert "If the finite disjoint-root witness refutes" in closure
    assert "full-eta target totality" in closure
    assert "exact #1243 runtime-rebind CI" in closure
    assert "RF30->RF39" in closure
    assert "[24,48,96]" in closure
    assert "momentum <=1e-3" in closure
    assert "divergence <=1e-5" in closure


def test_build_artifact_validates() -> None:
    artifact = r.build_artifact()
    r.validate_artifact(artifact)
    assert artifact["truth_boundary"]["complete_ns_residual_assessed"] is False
    assert artifact["truth_boundary"]["pde_validated"] is False


@pytest.mark.parametrize(
    ("path", "value"),
    [
        (["core_states", "leading_ready"], True),
        (["core_states", "correction_ready"], True),
        (["core_states", "velocity_export_ready"], True),
        (["core_states", "pde_validated"], True),
        (["truth_boundary", "complete_ns_residual_assessed"], True),
        (["truth_boundary", "scientifically_admitted"], True),
        (["truth_boundary", "pde_validated"], True),
        (["leading_route", "agent1", "actual_refutation_result_ingested"], True),
        (["leading_route", "agent1", "common_scalar_bridge_refuted"], True),
        (["leading_route", "agent4", "scoped_gate_passed"], True),
        (["correction_route", "selected_route", "exact_backend_rebind_promoted"], True),
        (["correction_route", "selected_route", "rf30_to_rf39_candidate_execution_recorded"], True),
        (["correction_route", "duplicate_nonselected_route", "evidence_unioned_into_selected_route"], True),
        (["fixed_validation_protocol", "normalized_momentum_sampled_max_gate"], 2e-3),
    ],
)
def test_truth_promotions_fail_closed(path: list[str], value: object) -> None:
    r.assert_mutation_fails(path, value)


def test_artifact_round_trip(tmp_path: Path) -> None:
    path = r.write_artifact(tmp_path / "routing.json")
    loaded = json.loads(path.read_text())
    assert loaded == r.build_artifact()
    r.validate_artifact(loaded)
