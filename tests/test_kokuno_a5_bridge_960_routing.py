from __future__ import annotations

import json
from pathlib import Path

import pytest

from openai_ns_reconstruction import kokuno_a5_bridge_960_routing as r


def test_exact_parent_and_upstream_heads_are_frozen() -> None:
    assert r.PARENT_A5_EXACT_HEAD == "9784da1596a71d7d60dc4759f7571cb33f17d904"
    assert r.AGENT1_EXACT_HEAD == "6ae44fdac17483b1fd042dca6f5bdf06701ae068"
    assert r.AGENT4_EXACT_HEAD == "ae26fda9c373f7cfb08a43a0ad48977ab13bcc71"
    assert r.AGENT2_EXACT_HEAD == "f90bdb9c232fb16df881dfbcac28be1cd5c29780"
    assert r.AGENT3_EXACT_HEAD == "20120093abbd981e2d2cb5d95093cea1eebb28ef"


def test_agent1_bridge_transport_is_not_bridge_geometry() -> None:
    a1 = r.leading_route()["agent1"]
    assert a1["current_l_minus_h_qs_transport_materialized"] is True
    assert a1["full_eta_common_scalar_bridge_length_established"] is False
    assert a1["eta_dependent_cartesian_bridge_geometry_materialized"] is False
    assert a1["current_cartesian_terminal_multiplier_composed"] is False


def test_matching_agent4_audit_is_scoped_and_unresolved() -> None:
    a4 = r.leading_route()["agent4"]
    assert a4["audits_agent1_exact_head"] == r.AGENT1_EXACT_HEAD
    assert a4["matching_independent_audit_surface_present"] is True
    assert a4["scientific_result_ingested"] is False
    assert a4["scoped_gate_passed"] is None
    assert a4["audit_is_complete_ns_validation"] is False
    assert a4["pde_validated"] is False


def test_eta_representation_firewall_is_fail_closed() -> None:
    fw = r.leading_route()["representation_firewall"]
    assert fw
    assert not any(fw.values())


def test_correction_routes_have_same_logical_capability() -> None:
    routes = r.correction_routes()
    a2 = routes["agent2_route"]
    a3 = routes["agent3_route"]
    keys = [
        "rf44_prestate_blob",
        "exact_backend_firewall_blob",
        "a2_960_differential_blob",
        "a2_960_axis_blob",
        "a2_960_vorticity_blob",
        "rf44_prestate_present",
        "exact_backend_firewall_present",
        "exact_a2_960_runtime_present",
        "exact_a2_1080_runtime_present",
        "real_exact_backend_bind_executed",
    ]
    for key in keys:
        assert a2[key] == a3[key]


def test_duplicate_routes_are_not_double_counted_or_unioned() -> None:
    topo = r.correction_routes()["topology"]
    assert topo["agent2_1234_and_agent3_1236_are_diverged_siblings"] is True
    assert topo["agent2_is_ancestor_of_agent3"] is False
    assert topo["agent3_is_ancestor_of_agent2"] is False
    assert topo["old_1224_1226_unified_runtime_seam_structurally_closed"] is True
    assert topo["logical_capability_count"] == 1
    assert topo["double_count_duplicate_routes"] is False
    assert topo["silent_cross_branch_evidence_union_allowed"] is False
    assert topo["downstream_must_consume_exactly_one_unified_route"] is True
    assert topo["canonical_downstream_route_selected_here"] is False


def test_a2_1080_still_blocks_real_backend_execution() -> None:
    for route_name in ("agent2_route", "agent3_route"):
        route = r.correction_routes()[route_name]
        assert route["exact_a2_1080_runtime_present"] is False
        assert route["real_exact_backend_bind_executed"] is False
        assert route["real_rf30_to_rf39_candidate_execution_recorded"] is False
        assert route["rf44_post_update_state_materialized"] is False
        assert route["cartesian_delta_u_materialized"] is False
        assert route["finite_correction_cycle_run"] is False


def test_latest_self_contained_candidate_does_not_move() -> None:
    candidate = r.latest_self_contained_candidate()
    assert candidate == {
        "agent2_pr": 1117,
        "agent1_pr": 1107,
        "stage": "xi=11",
        "self_contained_leading_plus_oscillatory_velocity_xyzt": True,
        "superseded_by_agent1_1233": False,
        "superseded_by_agent2_1234": False,
        "superseded_by_agent3_1236": False,
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
    assert p["green_ci_is_pde_validation"] is False
    assert p["queued_or_running_ci_is_pass"] is False


def test_closed_and_open_blockers_are_distinguished() -> None:
    b = r.remaining_blockers()
    assert b["actual_current_release2_qs"] is False
    assert b["current_l_minus_h_qs_transport"] is False
    assert b["matching_agent4_bridge_audit_surface"] is False
    assert b["rf44_prestate_firewall_a2_960_colocation"] is False
    assert b["duplicate_unified_960_route_deduplication"] is False
    assert b["agent4_bridge_scoped_result_resolved"] is True
    assert b["full_eta_bridge_geometry"] is True
    assert b["exact_a2_1080_runtime_on_one_unified_route"] is True
    assert b["real_exact_backend_bind"] is True
    assert b["complete_identity_bound_ns_defect"] is True
    assert b["agent4_heldout_canonical_complete_ns_gate"] is True


def test_shortest_closure_points_to_geometry_and_1080() -> None:
    closure = "\n".join(r.shortest_closure())
    assert "full-eta bridge geometry" in closure
    assert "exact A2 #1080" in closure
    assert "exactly one" in closure
    assert "[24,48,96]" in closure
    assert "momentum <=1e-3" in closure
    assert "divergence <=1e-5" in closure


def test_build_artifact_validates() -> None:
    artifact = r.build_artifact()
    r.validate_artifact(artifact)
    assert artifact["truth_boundary"]["pde_validated"] is False
    assert artifact["truth_boundary"]["complete_ns_residual_assessed"] is False


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
        (["correction_routes", "topology", "double_count_duplicate_routes"], True),
        (["correction_routes", "topology", "silent_cross_branch_evidence_union_allowed"], True),
        (["leading_route", "agent1", "full_eta_common_scalar_bridge_length_established"], True),
        (["leading_route", "agent4", "scoped_gate_passed"], True),
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
