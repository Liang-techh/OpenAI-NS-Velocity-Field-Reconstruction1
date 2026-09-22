from __future__ import annotations

import copy
import json

import pytest

from openai_ns_reconstruction.kokuno_a5_decimal_total_oriented_rf30_audit_routing import (
    AGENT1_EXACT_HEAD,
    AGENT2_EXACT_HEAD,
    AGENT3_EXACT_HEAD,
    AGENT4_EXACT_HEAD,
    LATEST_SELF_CONTAINED_A1_PR,
    LATEST_SELF_CONTAINED_A2_PR,
    PARENT_A5_EXACT_HEAD,
    KokunoA5RoutingError,
    build_artifact,
    dump_artifact,
    validate_artifact,
)


def test_exact_routes_and_latest_self_contained_identity() -> None:
    artifact = build_artifact()
    assert artifact["parent_a5"] == {
        "pr": 1173,
        "exact_head": PARENT_A5_EXACT_HEAD,
    }
    routes = artifact["upstream_routes"]
    assert routes["agent1"]["exact_head"] == AGENT1_EXACT_HEAD
    assert routes["agent2"]["exact_head"] == AGENT2_EXACT_HEAD
    assert routes["agent3"]["exact_head"] == AGENT3_EXACT_HEAD
    assert routes["agent4"]["exact_head"] == AGENT4_EXACT_HEAD

    candidate = artifact["latest_self_contained_candidate"]
    assert candidate["agent2_pr"] == LATEST_SELF_CONTAINED_A2_PR == 1117
    assert candidate["agent1_pr"] == LATEST_SELF_CONTAINED_A1_PR == 1107
    assert candidate["stage"] == "xi=11"
    assert candidate["self_contained_leading_plus_oscillatory_velocity_xyzt"] is True
    assert candidate["superseded_by_agent1_1179_as_unified_candidate"] is False
    assert candidate["superseded_by_agent2_1180_as_unified_candidate"] is False


def test_new_positive_facts_are_scoped_and_fail_closed() -> None:
    artifact = build_artifact()
    a1 = artifact["upstream_routes"]["agent1"]
    assert a1["fixed_decimal_significant_digits"] == 96
    assert a1["current_cartesian_relative_swirl_composed_precision_qualified"] is True
    assert a1["binary64_total_relative_swirl_sum_is_resolved"] is False
    assert a1["terminal_exterior_global_leading_materialized"] is False
    assert a1["matched_pressure_materialized"] is False

    a2 = artifact["upstream_routes"]["agent2"]
    assert a2["provider_driven"] is True
    assert a2["self_contained_velocity_xyzt_provider"] is False
    assert a2["source_exact_orientation_recovered"] is False
    assert a2["complete_ns_residual_assessed"] is False

    a3 = artifact["upstream_routes"]["agent3"]
    assert a3["typed_rf30_mechanics_materialized"] is True
    assert a3["pinned_same_identity_repository_provider_blob"] is None
    assert a3["repository_candidate_rf30_defect_scientifically_admitted"] is False
    assert a3["rf31_five_row_system_materialized"] is False
    assert a3["correction_applied"] is False

    a4 = artifact["upstream_routes"]["agent4"]
    assert a4["audits_agent1_exact_head"] == AGENT1_EXACT_HEAD
    assert a4["save_load_public_velocity_only_scientific_path"] is True
    assert a4["scoped_only"] is True
    assert a4["scoped_gate_passed"] is None
    assert a4["canonical_absolute_fd_transfer_ready"] is False
    assert a4["scientifically_admitted"] is False
    assert a4["pde_validated"] is False


def test_core_states_and_frozen_gates_remain_exact() -> None:
    artifact = build_artifact()
    assert artifact["core_states"] == {
        "leading_ready": False,
        "oscillatory_ready": True,
        "correction_ready": False,
        "velocity_export_ready": False,
        "pde_validated": False,
    }
    protocol = artifact["fixed_validation_protocol"]
    assert protocol["st006_momentum_sampled_max"] == 0.1082289305112118
    assert protocol["st006_momentum_volume_l2"] == 0.10758432876230622
    assert protocol["normalized_momentum_sampled_max_gate"] == 1.0e-3
    assert protocol["normalized_momentum_volume_l2_gate"] == 1.0e-3
    assert protocol["normalized_divergence_sampled_max_gate"] == 1.0e-5
    assert protocol["normalized_divergence_volume_l2_gate"] == 1.0e-5
    assert protocol["canonical_quadrature"] == [24, 48, 96]
    assert protocol["residual_defined_free_forcing_forbidden"] is True
    assert protocol["posthoc_threshold_relaxation_forbidden"] is True
    assert protocol["ci_success_is_scientific_admission"] is False


def test_firewall_rejects_cross_identity_and_scientific_promotions() -> None:
    artifact = build_artifact()

    promoted = copy.deepcopy(artifact)
    promoted["core_states"]["pde_validated"] = True
    with pytest.raises(KokunoA5RoutingError):
        validate_artifact(promoted)

    promoted = copy.deepcopy(artifact)
    promoted["upstream_routes"]["agent4"]["pde_validated"] = True
    with pytest.raises(KokunoA5RoutingError):
        validate_artifact(promoted)

    promoted = copy.deepcopy(artifact)
    promoted["upstream_routes"]["agent3"]["pinned_same_identity_repository_provider_blob"] = "fabricated"
    with pytest.raises(KokunoA5RoutingError):
        validate_artifact(promoted)

    promoted = copy.deepcopy(artifact)
    promoted["latest_self_contained_candidate"]["superseded_by_agent2_1180_as_unified_candidate"] = True
    with pytest.raises(KokunoA5RoutingError):
        validate_artifact(promoted)

    promoted = copy.deepcopy(artifact)
    promoted["fixed_validation_protocol"]["normalized_momentum_sampled_max_gate"] = 2.0e-3
    with pytest.raises(KokunoA5RoutingError):
        validate_artifact(promoted)

    promoted = copy.deepcopy(artifact)
    promoted["fixed_validation_protocol"]["residual_defined_free_forcing_forbidden"] = False
    with pytest.raises(KokunoA5RoutingError):
        validate_artifact(promoted)


def test_deterministic_json_roundtrip(tmp_path) -> None:
    first = dump_artifact(tmp_path / "first.json")
    second = dump_artifact(tmp_path / "second.json")
    assert first.read_bytes() == second.read_bytes()
    parsed = json.loads(first.read_text(encoding="utf-8"))
    validate_artifact(parsed)
    assert all(parsed["remaining_blockers"].values())
    assert not any(parsed["truth_boundary"].values())
