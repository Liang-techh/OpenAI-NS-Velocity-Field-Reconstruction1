from __future__ import annotations

import copy
import json
from pathlib import Path

import pytest

from openai_ns_reconstruction.kokuno_a5_terminal_rf44_backend_routing import (
    AGENT1_EXACT_HEAD,
    AGENT1_PR,
    AGENT2_EXACT_HEAD,
    AGENT2_PR,
    AGENT2_PROVIDER_SOURCE_BLOB,
    AGENT3_COMMON_PARENT_EXACT_HEAD,
    AGENT3_EXACT_HEAD,
    AGENT3_PR,
    AGENT4_EXACT_HEAD,
    AGENT4_PR,
    CANONICAL_QUADRATURE,
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


def test_fresh_exact_frontier_and_parent_are_frozen() -> None:
    a = build_artifact()
    assert a["task"] == TASK
    assert a["parent_a5"] == {"pr": PARENT_A5_PR, "exact_head": PARENT_A5_EXACT_HEAD}
    routes = a["upstream_routes"]
    assert (routes["agent1"]["pr"], routes["agent1"]["exact_head"]) == (AGENT1_PR, AGENT1_EXACT_HEAD)
    assert (routes["agent2"]["pr"], routes["agent2"]["exact_head"]) == (AGENT2_PR, AGENT2_EXACT_HEAD)
    assert (routes["agent3"]["pr"], routes["agent3"]["exact_head"]) == (AGENT3_PR, AGENT3_EXACT_HEAD)
    assert (routes["agent4"]["pr"], routes["agent4"]["exact_head"]) == (AGENT4_PR, AGENT4_EXACT_HEAD)


def test_terminal_target_and_matching_a4_surface_close_only_scoped_structure() -> None:
    a = build_artifact()
    a1 = a["upstream_routes"]["agent1"]
    a4 = a["upstream_routes"]["agent4"]
    blockers = a["remaining_blockers"]
    assert a1["public_terminal_target_api_materialized"] is True
    assert a1["current_release2_end_q_s_materialized"] is False
    assert a1["l_minus_h_matching_bridge_materialized"] is False
    assert a1["cartesian_terminal_multiplier_composed"] is False
    assert a4["matching_terminal_target_audit_present"] is True
    assert a4["audits_agent1_pr"] == AGENT1_PR
    assert a4["scoped_gate_passed"] is None
    assert a4["cartesian_leading_divergence_assessed_here"] is False
    assert blockers["a1_terminal_qp_target_algebra"] is False
    assert blockers["matching_a4_terminal_target_audit_surface"] is False
    assert blockers["a1_current_release2_end_q_s"] is True
    assert blockers["a1_l_minus_h_matching_bridge_to_qp"] is True


def test_agent2_rf44_prestate_is_not_postupdate_or_cartesian_correction() -> None:
    route = build_artifact()["upstream_routes"]["agent2"]
    assert route["raw_auxiliary_t2_provider_source_blob"] == AGENT2_PROVIDER_SOURCE_BLOB
    assert route["rf44_pre_update_state_materialized"] is True
    assert route["rf44_pre_state_reduces_to_authenticated_rf30_on_frozen_premean_state"] is True
    assert route["rf44_post_update_state_materialized"] is False
    assert route["rf44_correction_increment_materialized"] is False
    assert route["cartesian_correction_velocity_materialized"] is False
    assert route["repository_candidate_scientific_evidence"] is False


def test_agent3_firewall_records_real_runtime_dependency_blocker() -> None:
    route = build_artifact()["upstream_routes"]["agent3"]
    assert route["exact_runtime_rebind_required"] is True
    assert route["pinned_agent2_1198_provider_available_in_parent_lineage"] is True
    assert route["concrete_a2_1080_runtime_module_present_in_exact_stack"] is False
    assert route["concrete_a2_960_differential_module_present_in_exact_stack"] is False
    assert route["real_exact_backend_rebind_executed"] is False
    assert route["real_current_candidate_compact_correction_numerically_executed"] is False
    assert route["direct_constructed_backend_fixture_is_scientific_evidence"] is False
    assert route["parent_1210_promotion_without_rebind_is_sufficient_scientific_evidence"] is False
    assert route["repository_candidate_correction_scientifically_admitted"] is False


def test_a2_a3_sibling_topology_prevents_silent_evidence_union() -> None:
    topology = build_artifact()["sibling_topology"]
    assert topology["common_parent_exact_head"] == AGENT3_COMMON_PARENT_EXACT_HEAD
    assert topology["agent2_and_agent3_are_sibling_descendants_of_1210"] is True
    assert topology["agent2_1215_is_ancestor_of_agent3_1214"] is False
    assert topology["agent3_1214_is_ancestor_of_agent2_1215"] is False
    assert topology["silent_cross_branch_evidence_union_allowed"] is False
    assert topology["unified_stack_containing_agent2_1215_and_agent3_1214_materialized"] is False


def test_latest_self_contained_candidate_does_not_move() -> None:
    c = build_artifact()["latest_self_contained_candidate"]
    assert c["agent2_pr"] == 1117
    assert c["agent1_pr"] == 1107
    assert c["stage"] == "xi=11"
    assert c["self_contained_leading_plus_oscillatory_velocity_xyzt"] is True
    assert c["superseded_by_agent1_1213_as_unified_candidate"] is False
    assert c["superseded_by_agent2_1215_as_unified_candidate"] is False
    assert c["superseded_by_agent3_1214_as_corrected_candidate"] is False


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
    assert p["ci_success_is_scientific_admission"] is False


def test_evidence_firewall_is_all_false() -> None:
    firewall = build_artifact()["evidence_firewall"]
    assert firewall
    assert all(value is False for value in firewall.values())


@pytest.mark.parametrize(
    ("path", "value"),
    [
        (("core_states", "pde_validated"), True),
        (("core_states", "correction_ready"), True),
        (("sibling_topology", "silent_cross_branch_evidence_union_allowed"), True),
        (("sibling_topology", "unified_stack_containing_agent2_1215_and_agent3_1214_materialized"), True),
        (("upstream_routes", "agent1", "cartesian_terminal_multiplier_composed"), True),
        (("upstream_routes", "agent2", "rf44_post_update_state_materialized"), True),
        (("upstream_routes", "agent2", "cartesian_correction_velocity_materialized"), True),
        (("upstream_routes", "agent3", "real_exact_backend_rebind_executed"), True),
        (("upstream_routes", "agent3", "repository_candidate_correction_scientifically_admitted"), True),
        (("upstream_routes", "agent4", "scoped_gate_passed"), True),
        (("truth_boundary", "scientifically_admitted"), True),
        (("truth_boundary", "complete_ns_residual_assessed"), True),
    ],
)
def test_false_promotion_mutations_fail_closed(path: tuple[str, ...], value: object) -> None:
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


def test_reopening_closed_or_closing_open_blockers_fails() -> None:
    a = copy.deepcopy(build_artifact())
    a["remaining_blockers"]["agent3_checksum_pin_of_exact_agent2_1198_provider_blob"] = True
    with pytest.raises(KokunoA5RoutingError):
        validate_artifact(a)

    a = copy.deepcopy(build_artifact())
    a["remaining_blockers"]["agent3_real_exact_backend_rebind"] = False
    with pytest.raises(KokunoA5RoutingError):
        validate_artifact(a)


def test_deterministic_artifact_roundtrip(tmp_path: Path) -> None:
    p1 = write_artifact(tmp_path / "a.json")
    p2 = write_artifact(tmp_path / "b.json")
    assert p1.read_bytes() == p2.read_bytes()
    loaded = json.loads(p1.read_text())
    validate_artifact(loaded)
    assert loaded == build_artifact()
