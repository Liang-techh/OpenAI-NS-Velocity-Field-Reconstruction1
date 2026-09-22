from __future__ import annotations

import copy
import json
from pathlib import Path

import pytest

from openai_ns_reconstruction.kokuno_a5_release2_auxt2_derivative_rf49_routing import (
    AGENT1_EXACT_HEAD,
    AGENT1_PR,
    AGENT2_EXACT_HEAD,
    AGENT2_PR,
    AGENT2_PROVIDER_SOURCE_BLOB,
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


def test_release2_and_derivative_audit_close_only_structural_availability() -> None:
    a = build_artifact()
    routes = a["upstream_routes"]
    blockers = a["remaining_blockers"]

    assert routes["agent1"]["source_l_minus1_to_minus_h_transition_materialized"] is True
    assert routes["agent1"]["source_terminal_multiplier_materialized"] is False
    assert routes["agent1"]["outer_global_leading_velocity_materialized"] is False

    assert routes["agent2"]["current_i4_raw_auxiliary_t2_provider_materialized"] is True
    assert routes["agent2"]["three_resolution_slow_derivative_diagnostic_materialized"] is True
    assert routes["agent2"]["provider_parent_source_blob"] == AGENT2_PROVIDER_SOURCE_BLOB
    assert routes["agent2"]["provider_checksum_pinned_by_agent3"] is False
    assert routes["agent2"]["rf30_repository_candidate_state_authorized"] is False

    assert blockers["a1_minus1_to_minus_h_transition"] is False
    assert blockers["a2_three_resolution_auxiliary_derivative_audit_missing"] is False
    assert blockers["agent3_checksum_pin_of_exact_agent2_1198_provider_blob"] is True


def test_agent3_unpinned_mechanics_remain_non_scientific() -> None:
    route = build_artifact()["upstream_routes"]["agent3"]
    assert route["pinned_same_identity_repository_provider_blob"] is None
    assert route["sees_agent2_1198_provider_as_pinned"] is False
    assert route["sees_agent2_1205_derivative_audit_as_scientific_authorization"] is False
    assert route["repository_candidate_rf30_evidence"] is False
    assert route["repository_candidate_rf44_rf49_evidence"] is False
    assert route["cartesian_correction_velocity_materialized"] is False
    assert route["finite_correction_cycle_run"] is False


def test_older_agent4_audit_is_not_transferred_to_release2() -> None:
    route = build_artifact()["upstream_routes"]["agent4"]
    assert route["audits_agent1_pr"] == 1196
    assert route["audits_latest_agent1_release2_pr"] is False
    assert route["matching_release2_audit_present"] is False
    assert route["scoped_gate_passed"] is None
    assert route["complete_ns_residual_assessed"] is False
    assert route["pde_validated"] is False


def test_latest_self_contained_candidate_does_not_move() -> None:
    c = build_artifact()["latest_self_contained_candidate"]
    assert c["agent2_pr"] == 1117
    assert c["agent1_pr"] == 1107
    assert c["stage"] == "xi=11"
    assert c["self_contained_leading_plus_oscillatory_velocity_xyzt"] is True
    assert c["superseded_by_agent1_1204_as_unified_candidate"] is False
    assert c["superseded_by_agent2_1205_as_unified_candidate"] is False
    assert c["superseded_by_agent3_1197_as_corrected_candidate"] is False


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
        (("upstream_routes", "agent2", "provider_checksum_pinned_by_agent3"), True),
        (("upstream_routes", "agent2", "rf30_repository_candidate_state_authorized"), True),
        (("upstream_routes", "agent3", "repository_candidate_rf30_evidence"), True),
        (("upstream_routes", "agent3", "cartesian_correction_velocity_materialized"), True),
        (("upstream_routes", "agent4", "matching_release2_audit_present"), True),
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


def test_deterministic_artifact_roundtrip(tmp_path: Path) -> None:
    p1 = write_artifact(tmp_path / "a.json")
    p2 = write_artifact(tmp_path / "b.json")
    assert p1.read_bytes() == p2.read_bytes()
    loaded = json.loads(p1.read_text())
    validate_artifact(loaded)
    assert loaded == build_artifact()
