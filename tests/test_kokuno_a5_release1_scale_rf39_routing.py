from __future__ import annotations

import copy
import json

import pytest

from openai_ns_reconstruction.kokuno_a5_release1_scale_rf39_routing import (
    AGENT1_EXACT_HEAD,
    AGENT2_EXACT_HEAD,
    AGENT3_EXACT_HEAD,
    AGENT4_AUDITED_A1_HEAD,
    AGENT4_EXACT_HEAD,
    PARENT_A5_EXACT_HEAD,
    KokunoA5RoutingError,
    build_artifact,
    dump_artifact,
    validate_artifact,
)


def test_exact_sibling_routes_and_identity_firewall() -> None:
    artifact = build_artifact()
    assert artifact["parent_a5"] == {"pr": 1183, "exact_head": PARENT_A5_EXACT_HEAD}
    routes = artifact["upstream_routes"]
    assert routes["agent1"]["exact_head"] == AGENT1_EXACT_HEAD
    assert routes["agent2"]["exact_head"] == AGENT2_EXACT_HEAD
    assert routes["agent3"]["exact_head"] == AGENT3_EXACT_HEAD
    assert routes["agent4"]["exact_head"] == AGENT4_EXACT_HEAD
    assert routes["agent4"]["audits_agent1_exact_head"] == AGENT4_AUDITED_A1_HEAD
    assert routes["agent4"]["audits_agent1_exact_head"] != AGENT1_EXACT_HEAD
    assert routes["agent4"]["audits_latest_agent1_release1"] is False


def test_fresh_positive_facts_remain_scoped() -> None:
    artifact = build_artifact()
    a1 = artifact["upstream_routes"]["agent1"]
    assert a1["post_relative_swirl_release1_materialized"] is True
    assert a1["fixed_decimal_significant_digits"] == 96
    assert a1["source_l_minus1_hold_materialized"] is False
    assert a1["source_l_minus1_to_minus_h_transition_materialized"] is False
    assert a1["source_terminal_multiplier_materialized"] is False
    assert a1["source_exterior_heat_replacement_materialized"] is False
    assert a1["outer_global_leading_velocity_materialized"] is False
    assert a1["matched_global_pressure_materialized"] is False

    a2 = artifact["upstream_routes"]["agent2"]
    assert a2["bounded_isotropic_scale_materialized"] is True
    assert a2["scale_bounds"] == [0.5, 2.0]
    assert a2["provider_driven"] is True
    assert a2["self_contained_velocity_xyzt_provider"] is False
    assert a2["source_exact_scale_recovered"] is False
    assert a2["same_identity_auxiliary_t2_provider_checksum_pinned"] is False

    a3 = artifact["upstream_routes"]["agent3"]
    assert a3["typed_rf30_to_rf34_rf39_chain_executable"] is True
    assert a3["rf31_basis_frozen_before_defect"] is True
    assert a3["pinned_same_identity_repository_provider_blob"] is None
    assert a3["repository_candidate_scientific_correction_materialized"] is False
    assert a3["correction_applied_to_candidate"] is False
    assert a3["cartesian_correction_velocity_materialized"] is False
    assert a3["finite_correction_cycle_run"] is False


def test_latest_unified_candidate_and_core_states_are_not_promoted() -> None:
    artifact = build_artifact()
    candidate = artifact["latest_self_contained_candidate"]
    assert candidate["agent2_pr"] == 1117
    assert candidate["agent1_pr"] == 1107
    assert candidate["stage"] == "xi=11"
    assert candidate["self_contained_leading_plus_oscillatory_velocity_xyzt"] is True
    assert candidate["superseded_by_agent1_1188_as_unified_candidate"] is False
    assert candidate["superseded_by_agent2_1189_as_unified_candidate"] is False
    assert candidate["superseded_by_agent3_1190_as_corrected_candidate"] is False
    assert artifact["core_states"] == {
        "leading_ready": False,
        "oscillatory_ready": True,
        "correction_ready": False,
        "velocity_export_ready": False,
        "pde_validated": False,
    }


def test_frozen_st006_and_project_gates_are_exact() -> None:
    protocol = build_artifact()["fixed_validation_protocol"]
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


def test_mutations_fail_closed() -> None:
    artifact = build_artifact()

    bad = copy.deepcopy(artifact)
    bad["upstream_routes"]["agent1"]["outer_global_leading_velocity_materialized"] = True
    with pytest.raises(KokunoA5RoutingError):
        validate_artifact(bad)

    bad = copy.deepcopy(artifact)
    bad["upstream_routes"]["agent2"]["self_contained_velocity_xyzt_provider"] = True
    with pytest.raises(KokunoA5RoutingError):
        validate_artifact(bad)

    bad = copy.deepcopy(artifact)
    bad["upstream_routes"]["agent3"]["pinned_same_identity_repository_provider_blob"] = "fabricated"
    with pytest.raises(KokunoA5RoutingError):
        validate_artifact(bad)

    bad = copy.deepcopy(artifact)
    bad["upstream_routes"]["agent3"]["correction_applied_to_candidate"] = True
    with pytest.raises(KokunoA5RoutingError):
        validate_artifact(bad)

    bad = copy.deepcopy(artifact)
    bad["upstream_routes"]["agent4"]["audits_latest_agent1_release1"] = True
    with pytest.raises(KokunoA5RoutingError):
        validate_artifact(bad)

    bad = copy.deepcopy(artifact)
    bad["latest_self_contained_candidate"]["superseded_by_agent2_1189_as_unified_candidate"] = True
    with pytest.raises(KokunoA5RoutingError):
        validate_artifact(bad)

    bad = copy.deepcopy(artifact)
    bad["core_states"]["pde_validated"] = True
    with pytest.raises(KokunoA5RoutingError):
        validate_artifact(bad)

    bad = copy.deepcopy(artifact)
    bad["fixed_validation_protocol"]["normalized_momentum_sampled_max_gate"] = 2.0e-3
    with pytest.raises(KokunoA5RoutingError):
        validate_artifact(bad)

    bad = copy.deepcopy(artifact)
    bad["fixed_validation_protocol"]["residual_defined_free_forcing_forbidden"] = False
    with pytest.raises(KokunoA5RoutingError):
        validate_artifact(bad)


def test_deterministic_artifact_roundtrip(tmp_path) -> None:
    first = dump_artifact(tmp_path / "first.json")
    second = dump_artifact(tmp_path / "second.json")
    assert first.read_bytes() == second.read_bytes()
    parsed = json.loads(first.read_text(encoding="utf-8"))
    validate_artifact(parsed)
    assert all(parsed["remaining_blockers"].values())
    assert not any(parsed["evidence_firewall"].values())
    assert not any(parsed["truth_boundary"].values())
