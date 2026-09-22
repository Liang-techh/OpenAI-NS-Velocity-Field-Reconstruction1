from __future__ import annotations

import copy
import json
from pathlib import Path

import pytest

from openai_ns_reconstruction.kokuno_a5_current_qs_rf44_960_routing import (
    AGENT1_EXACT_HEAD,
    AGENT1_PR,
    AGENT2_EXACT_HEAD,
    AGENT2_PR,
    AGENT3_EXACT_HEAD,
    AGENT3_PR,
    AGENT4_EXACT_HEAD,
    AGENT4_PR,
    CANONICAL_QUADRATURE,
    CORRECTION_COMMON_PARENT_HEAD,
    FINAL_NORMALIZED_DIVERGENCE_GATE,
    FINAL_NORMALIZED_MOMENTUM_GATE,
    KokunoA5RoutingError,
    PARENT_A5_EXACT_HEAD,
    PARENT_A5_PR,
    ST006_MOMENTUM_SAMPLED_MAX,
    ST006_MOMENTUM_VOLUME_L2,
    TASK,
    build_artifact,
    validate_artifact,
    write_artifact,
)


def test_fresh_frontier_and_parent_are_frozen() -> None:
    a = build_artifact()
    assert a["task"] == TASK
    assert a["parent_a5"] == {"pr": PARENT_A5_PR, "exact_head": PARENT_A5_EXACT_HEAD}
    routes = a["upstream_routes"]
    assert (routes["agent1"]["pr"], routes["agent1"]["exact_head"]) == (AGENT1_PR, AGENT1_EXACT_HEAD)
    assert (routes["agent2"]["pr"], routes["agent2"]["exact_head"]) == (AGENT2_PR, AGENT2_EXACT_HEAD)
    assert (routes["agent3"]["pr"], routes["agent3"]["exact_head"]) == (AGENT3_PR, AGENT3_EXACT_HEAD)
    assert (routes["agent4"]["pr"], routes["agent4"]["exact_head"]) == (AGENT4_PR, AGENT4_EXACT_HEAD)


def test_actual_current_qs_closes_only_scalar_state_blocker() -> None:
    a = build_artifact()
    a1 = a["upstream_routes"]["agent1"]
    blockers = a["remaining_blockers"]
    assert a1["current_q_s_release2_endpoint_materialized"] is True
    assert a1["source_ideal_q_s_relabelled_as_current"] is False
    assert a1["current_l_minus_h_matching_bridge_materialized"] is False
    assert a1["current_cartesian_terminal_multiplier_composed"] is False
    assert blockers["a1_current_release2_end_q_s"] is False
    assert blockers["a1_current_l_minus_h_matching_bridge_to_qp"] is True


def test_source_ideal_a4_audit_does_not_transfer_to_current_qs() -> None:
    a = build_artifact()
    a4 = a["upstream_routes"]["agent4"]
    topology = a["branch_topology"]
    assert a4["source_ideal_q_s_schedule_audited"] is True
    assert a4["current_candidate_q_s_audited"] is False
    assert a4["matching_audit_for_agent1_1225_present"] is False
    assert topology["a1_1225_and_a4_1223_are_diverged_siblings"] is True
    assert topology["a4_1223_is_ancestor_of_a1_1225"] is False
    assert topology["a1_1225_is_ancestor_of_a4_1223"] is False
    assert a["remaining_blockers"]["matching_a4_actual_current_qs_audit"] is True


def test_a2_restack_closes_old_sibling_seam_but_not_runtime_bind() -> None:
    route = build_artifact()["upstream_routes"]["agent2"]
    assert route["rf44_pre_update_state_materialized"] is True
    assert route["exact_backend_firewall_present"] is True
    assert route["concrete_a2_960_runtime_present"] is False
    assert route["concrete_a2_1080_runtime_present"] is False
    assert route["real_exact_backend_bind_executed"] is False
    assert route["rf44_post_update_state_materialized"] is False
    assert route["cartesian_delta_u_materialized"] is False


def test_a3_restores_only_960_half_on_parallel_lineage() -> None:
    route = build_artifact()["upstream_routes"]["agent3"]
    topology = build_artifact()["branch_topology"]
    assert route["exact_backend_firewall_present"] is True
    assert route["concrete_a2_960_runtime_present"] is True
    assert route["concrete_a2_1080_runtime_present"] is False
    assert route["rf44_pre_update_state_from_agent2_1215_present"] is False
    assert route["real_exact_backend_bind_executed"] is False
    assert route["real_rf30_to_rf39_candidate_execution_recorded"] is False
    assert topology["correction_common_parent_head"] == CORRECTION_COMMON_PARENT_HEAD
    assert topology["a2_1224_and_a3_1226_are_diverged_descendants"] is True
    assert topology["a2_1224_is_ancestor_of_a3_1226"] is False
    assert topology["a3_1226_is_ancestor_of_a2_1224"] is False
    assert topology["silent_cross_branch_evidence_union_allowed"] is False
    assert topology["single_stack_with_rf44_prestate_firewall_and_a2_960_runtime_materialized"] is False


def test_latest_self_contained_candidate_does_not_move() -> None:
    c = build_artifact()["latest_self_contained_candidate"]
    assert c["agent2_pr"] == 1117
    assert c["agent1_pr"] == 1107
    assert c["stage"] == "xi=11"
    assert c["self_contained_leading_plus_oscillatory_velocity_xyzt"] is True
    assert c["superseded_by_agent1_1225"] is False
    assert c["superseded_by_agent2_1224"] is False
    assert c["superseded_by_agent3_1226"] is False


def test_core_state_and_validation_contract_remain_frozen() -> None:
    a = build_artifact()
    assert a["core_states"] == {
        "leading_ready": False,
        "oscillatory_ready": True,
        "correction_ready": False,
        "velocity_export_ready": False,
        "pde_validated": False,
    }
    p = a["fixed_validation_protocol"]
    assert p["st006_momentum_sampled_max"] == ST006_MOMENTUM_SAMPLED_MAX
    assert p["st006_momentum_volume_l2"] == ST006_MOMENTUM_VOLUME_L2
    assert p["normalized_momentum_sampled_max_gate"] == FINAL_NORMALIZED_MOMENTUM_GATE
    assert p["normalized_momentum_volume_l2_gate"] == FINAL_NORMALIZED_MOMENTUM_GATE
    assert p["normalized_divergence_sampled_max_gate"] == FINAL_NORMALIZED_DIVERGENCE_GATE
    assert p["normalized_divergence_volume_l2_gate"] == FINAL_NORMALIZED_DIVERGENCE_GATE
    assert p["canonical_quadrature"] == CANONICAL_QUADRATURE
    assert p["residual_defined_free_forcing_forbidden"] is True
    assert p["posthoc_threshold_relaxation_forbidden"] is True
    assert p["green_ci_is_pde_validation"] is False


def test_evidence_firewall_is_all_false() -> None:
    fw = build_artifact()["evidence_firewall"]
    assert fw
    assert all(value is False for value in fw.values())


@pytest.mark.parametrize(
    ("path", "value"),
    [
        (("core_states", "pde_validated"), True),
        (("core_states", "correction_ready"), True),
        (("branch_topology", "silent_cross_branch_evidence_union_allowed"), True),
        (("branch_topology", "single_stack_with_rf44_prestate_firewall_and_a2_960_runtime_materialized"), True),
        (("upstream_routes", "agent1", "current_cartesian_terminal_multiplier_composed"), True),
        (("upstream_routes", "agent2", "real_exact_backend_bind_executed"), True),
        (("upstream_routes", "agent2", "cartesian_delta_u_materialized"), True),
        (("upstream_routes", "agent3", "concrete_a2_1080_runtime_present"), True),
        (("upstream_routes", "agent3", "real_rf30_to_rf39_candidate_execution_recorded"), True),
        (("upstream_routes", "agent4", "current_candidate_q_s_audited"), True),
        (("upstream_routes", "agent4", "matching_audit_for_agent1_1225_present"), True),
        (("truth_boundary", "scientifically_admitted"), True),
        (("truth_boundary", "pde_validated"), True),
    ],
)
def test_false_promotions_fail_closed(path: tuple[str, ...], value: object) -> None:
    a = copy.deepcopy(build_artifact())
    cursor = a
    for key in path[:-1]:
        cursor = cursor[key]
    cursor[path[-1]] = value
    with pytest.raises(KokunoA5RoutingError):
        validate_artifact(a)


def test_gate_and_forcing_mutations_fail_closed() -> None:
    a = copy.deepcopy(build_artifact())
    a["fixed_validation_protocol"]["normalized_momentum_sampled_max_gate"] = 2e-3
    with pytest.raises(KokunoA5RoutingError):
        validate_artifact(a)

    a = copy.deepcopy(build_artifact())
    a["fixed_validation_protocol"]["residual_defined_free_forcing_forbidden"] = False
    with pytest.raises(KokunoA5RoutingError):
        validate_artifact(a)


def test_blocker_truth_cannot_be_silently_promoted() -> None:
    a = copy.deepcopy(build_artifact())
    a["remaining_blockers"]["matching_a4_actual_current_qs_audit"] = False
    with pytest.raises(KokunoA5RoutingError):
        validate_artifact(a)

    a = copy.deepcopy(build_artifact())
    a["remaining_blockers"]["a2_a3_1224_1226_unified_runtime_lineage"] = False
    with pytest.raises(KokunoA5RoutingError):
        validate_artifact(a)

    a = copy.deepcopy(build_artifact())
    a["remaining_blockers"]["a1_current_release2_end_q_s"] = True
    with pytest.raises(KokunoA5RoutingError):
        validate_artifact(a)


def test_deterministic_artifact_roundtrip(tmp_path: Path) -> None:
    p1 = write_artifact(tmp_path / "a.json")
    p2 = write_artifact(tmp_path / "b.json")
    assert p1.read_bytes() == p2.read_bytes()
    loaded = json.loads(p1.read_text())
    validate_artifact(loaded)
    assert loaded == build_artifact()
