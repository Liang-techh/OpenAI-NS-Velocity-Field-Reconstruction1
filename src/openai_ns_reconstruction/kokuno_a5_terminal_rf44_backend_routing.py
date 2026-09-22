"""Fail-closed Kokuno Agent-5 routing artifact for terminal/RF44/backend seams.

This module is integration/provenance glue only. It freshness-binds the latest
Agent-1--4 frontier, records sibling-branch topology, and keeps complete-NS/PDE
states fail-closed. It does not reimplement upstream mathematics.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

SCHEMA_VERSION = "kokuno-a5-terminal-rf44-backend-routing-v1"
TASK = "KOKUNO-A5-TERMINAL-RF44-BACKEND-ROUTING-128"

PARENT_A5_PR = 1208
PARENT_A5_EXACT_HEAD = "78b055f96458f9c029ca8ef9243f75f00c3341cd"

AGENT1_PR = 1213
AGENT1_EXACT_HEAD = "fb81b8e90272caa9e8bbd075f05b94dbb1195310"
AGENT1_SOURCE_BLOB = "fd3236607fba9c3060a1dc3058d2cec7272969d6"
AGENT1_PATH = "src/openai_ns_reconstruction/kokuno_public_terminal_multiplier_target.py"

AGENT2_PR = 1215
AGENT2_EXACT_HEAD = "cab837f6358c8d29e3d470061879a58b3f7e1454"
AGENT2_SOURCE_BLOB = "0d8bb6471a4b7825b1b4b428d3b2c38d207729fd"
AGENT2_PATH = "src/openai_ns_reconstruction/kokuno_current_i4_rf44_prestate_provider.py"

AGENT3_COMMON_PARENT_PR = 1210
AGENT3_COMMON_PARENT_EXACT_HEAD = "e75940c33127b0725ddc703eaae97590158f2a56"
AGENT3_COMMON_PARENT_SOURCE_BLOB = "b99b3bee2b2dad6003c11acae7946c9be53347fb"

AGENT3_PR = 1214
AGENT3_EXACT_HEAD = "25f2793091cbf8dc2d8f966b941c32eaff63d876"
AGENT3_SOURCE_BLOB = "0152a85f7fbaf9c32c24c9d5c7d0d11166e3b157"
AGENT3_PATH = "src/openai_ns_reconstruction/kokuno_current_i4_exact_backend_correction_firewall.py"

AGENT4_PR = 1216
AGENT4_EXACT_HEAD = "7cb6ca4c5f771bc7c1ae216730852eaa96d83bea"
AGENT4_SOURCE_BLOB = "8a82b1c1b6bf8d155c746254f5f58123654f10fe"
AGENT4_PATH = "src/openai_ns_reconstruction/kokuno_a4_terminal_target_independent_audit.py"

AGENT2_PROVIDER_PR = 1198
AGENT2_PROVIDER_SOURCE_BLOB = "96168ac6583ad5e71aa044bd558c6feeadf051f4"

LATEST_SELF_CONTAINED_A2_PR = 1117
LATEST_SELF_CONTAINED_A2_HEAD = "27741d9c0a27262f7fabf61eebaa2fbd507e9f03"
LATEST_SELF_CONTAINED_A1_PR = 1107
LATEST_SELF_CONTAINED_STAGE = "xi=11"

ST006_MOMENTUM_SAMPLED_MAX = 0.1082289305112118
ST006_MOMENTUM_VOLUME_L2 = 0.10758432876230622
FINAL_NORMALIZED_MOMENTUM_GATE = 1.0e-3
FINAL_NORMALIZED_DIVERGENCE_GATE = 1.0e-5
CANONICAL_QUADRATURE = [24, 48, 96]


class KokunoA5RoutingError(RuntimeError):
    """Raised when routing identity or truth-state data is promoted incorrectly."""


def sibling_topology() -> dict[str, Any]:
    return {
        "common_parent_pr": AGENT3_COMMON_PARENT_PR,
        "common_parent_exact_head": AGENT3_COMMON_PARENT_EXACT_HEAD,
        "common_parent_source_blob": AGENT3_COMMON_PARENT_SOURCE_BLOB,
        "agent2_pr": AGENT2_PR,
        "agent2_exact_head": AGENT2_EXACT_HEAD,
        "agent3_pr": AGENT3_PR,
        "agent3_exact_head": AGENT3_EXACT_HEAD,
        "agent2_and_agent3_are_sibling_descendants_of_1210": True,
        "agent2_1215_is_ancestor_of_agent3_1214": False,
        "agent3_1214_is_ancestor_of_agent2_1215": False,
        "silent_cross_branch_evidence_union_allowed": False,
        "unified_stack_containing_agent2_1215_and_agent3_1214_materialized": False,
    }


def upstream_routes() -> dict[str, dict[str, Any]]:
    return {
        "agent1": {
            "pr": AGENT1_PR,
            "exact_head": AGENT1_EXACT_HEAD,
            "source_blob": AGENT1_SOURCE_BLOB,
            "path": AGENT1_PATH,
            "delivery": "public terminal-multiplier target algebra and deterministic Q_p target",
            "release2_cartesian_parent_pr": 1204,
            "public_terminal_target_api_materialized": True,
            "repository_autonomous_c_o": 1.0 / 64.0,
            "current_autonomous_h": 0.005,
            "q_p_approx": 5.746912488789217e-4,
            "current_release2_end_q_s_materialized": False,
            "l_minus_h_matching_bridge_materialized": False,
            "cartesian_terminal_multiplier_composed": False,
            "exterior_heat_replacement_materialized": False,
            "outer_global_leading_velocity_materialized": False,
            "matched_global_pressure_materialized": False,
            "cross_language_global_velocity_export_ready": False,
            "heldout_complete_ns_residual_assessed": False,
            "exact_head_ci_conclusion_ingested": False,
        },
        "agent2": {
            "pr": AGENT2_PR,
            "exact_head": AGENT2_EXACT_HEAD,
            "source_blob": AGENT2_SOURCE_BLOB,
            "path": AGENT2_PATH,
            "delivery": "same-identity current-I4 RF44 pre-update state provider",
            "raw_auxiliary_t2_provider_pr": AGENT2_PROVIDER_PR,
            "raw_auxiliary_t2_provider_source_blob": AGENT2_PROVIDER_SOURCE_BLOB,
            "rf44_pre_update_state_materialized": True,
            "rf44_pre_state_reduces_to_authenticated_rf30_on_frozen_premean_state": True,
            "rf44_extension_has_separate_semantic_identity": True,
            "rf44_post_update_state_materialized": False,
            "rf44_correction_increment_materialized": False,
            "cartesian_correction_velocity_materialized": False,
            "finite_correction_cycle_run": False,
            "repository_candidate_scientific_evidence": False,
            "complete_ns_residual_assessed": False,
            "exact_head_ci_conclusion_ingested": False,
        },
        "agent3": {
            "pr": AGENT3_PR,
            "exact_head": AGENT3_EXACT_HEAD,
            "source_blob": AGENT3_SOURCE_BLOB,
            "path": AGENT3_PATH,
            "delivery": "exact-runtime rebind firewall before current-I4 compact correction evidence",
            "parent_pr": AGENT3_COMMON_PARENT_PR,
            "pinned_agent2_1198_provider_available_in_parent_lineage": True,
            "rf30_to_rf39_parent_mechanics_available": True,
            "exact_runtime_rebind_required": True,
            "concrete_a2_1080_runtime_module_present_in_exact_stack": False,
            "concrete_a2_960_differential_module_present_in_exact_stack": False,
            "real_exact_backend_rebind_executed": False,
            "real_current_candidate_compact_correction_numerically_executed": False,
            "direct_constructed_backend_fixture_is_scientific_evidence": False,
            "parent_1210_promotion_without_rebind_is_sufficient_scientific_evidence": False,
            "repository_candidate_correction_scientifically_admitted": False,
            "cartesian_correction_velocity_materialized": False,
            "rf44_rf49_postupdate_remainder_materialized_on_real_candidate": False,
            "finite_correction_cycle_run": False,
            "heldout_complete_ns_residual_assessed": False,
            "exact_head_ci_conclusion_ingested": False,
        },
        "agent4": {
            "pr": AGENT4_PR,
            "exact_head": AGENT4_EXACT_HEAD,
            "source_blob": AGENT4_SOURCE_BLOB,
            "path": AGENT4_PATH,
            "delivery": "implementation-distinct independent audit of exact A1 #1213 terminal Q_p target algebra",
            "audits_agent1_pr": AGENT1_PR,
            "matching_terminal_target_audit_present": True,
            "scoped_only": True,
            "scoped_gate_passed": None,
            "cartesian_leading_divergence_assessed_here": False,
            "complete_ns_residual_assessed": False,
            "scientifically_admitted": False,
            "pde_validated": False,
            "exact_head_ci_conclusion_ingested": False,
        },
    }


def latest_self_contained_candidate() -> dict[str, Any]:
    return {
        "agent2_pr": LATEST_SELF_CONTAINED_A2_PR,
        "agent2_exact_head": LATEST_SELF_CONTAINED_A2_HEAD,
        "agent1_pr": LATEST_SELF_CONTAINED_A1_PR,
        "stage": LATEST_SELF_CONTAINED_STAGE,
        "self_contained_leading_plus_oscillatory_velocity_xyzt": True,
        "superseded_by_agent1_1213_as_unified_candidate": False,
        "superseded_by_agent2_1215_as_unified_candidate": False,
        "superseded_by_agent3_1214_as_corrected_candidate": False,
    }


def core_states() -> dict[str, bool]:
    return {
        "leading_ready": False,
        "oscillatory_ready": True,
        "correction_ready": False,
        "velocity_export_ready": False,
        "pde_validated": False,
    }


def fixed_validation_protocol() -> dict[str, Any]:
    return {
        "st006_momentum_sampled_max": ST006_MOMENTUM_SAMPLED_MAX,
        "st006_momentum_volume_l2": ST006_MOMENTUM_VOLUME_L2,
        "normalized_momentum_sampled_max_gate": FINAL_NORMALIZED_MOMENTUM_GATE,
        "normalized_momentum_volume_l2_gate": FINAL_NORMALIZED_MOMENTUM_GATE,
        "normalized_divergence_sampled_max_gate": FINAL_NORMALIZED_DIVERGENCE_GATE,
        "normalized_divergence_volume_l2_gate": FINAL_NORMALIZED_DIVERGENCE_GATE,
        "canonical_quadrature": list(CANONICAL_QUADRATURE),
        "residual_defined_free_forcing_forbidden": True,
        "posthoc_threshold_relaxation_forbidden": True,
        "ci_success_is_scientific_admission": False,
        "scoped_operator_consistency_is_complete_ns_validation": False,
    }


def evidence_firewall() -> dict[str, bool]:
    return {
        "agent1_1213_target_algebra_is_cartesian_terminal_velocity": False,
        "agent4_1216_target_audit_is_cartesian_divergence_audit": False,
        "agent4_1216_target_audit_is_complete_ns_validation": False,
        "agent2_1215_prestate_is_postupdate_state": False,
        "agent2_1215_prestate_is_cartesian_correction": False,
        "agent2_1215_and_agent3_1214_sibling_evidence_can_be_silently_unioned": False,
        "agent3_1210_parent_promotion_bypasses_1214_exact_backend_firewall": False,
        "agent3_direct_constructor_fixture_is_real_candidate_execution": False,
        "queued_or_running_ci_means_pass": False,
        "kokuno_replay_is_independent_final_validation": False,
    }


def remaining_blockers() -> dict[str, bool]:
    """True means the blocker remains unresolved."""
    return {
        "a1_terminal_qp_target_algebra": False,
        "matching_a4_terminal_target_audit_surface": False,
        "a1_current_release2_end_q_s": True,
        "a1_l_minus_h_matching_bridge_to_qp": True,
        "a1_cartesian_terminal_multiplier": True,
        "a1_exterior_heat_and_global_leading": True,
        "matched_cartesian_pressure_and_grad_p": True,
        "cross_language_global_velocity_export_representation": True,
        "matching_global_self_contained_a2_composite": True,
        "raw_same_identity_auxiliary_t2_provider_missing": False,
        "agent3_checksum_pin_of_exact_agent2_1198_provider_blob": False,
        "rf44_pre_update_state_provider": False,
        "agent3_exact_a2_1080_runtime_dependency_in_stack": True,
        "agent3_exact_a2_960_differential_dependency_in_stack": True,
        "agent3_real_exact_backend_rebind": True,
        "a2_1215_a3_1214_unified_same_identity_stack": True,
        "rf44_post_update_state_and_increment_provider": True,
        "scientific_candidate_specific_rf30_to_rf49_correction": True,
        "cartesian_correction_velocity": True,
        "rf44_rf49_repository_candidate_remainder": True,
        "real_finite_correction_cycle": True,
        "preregistered_restricted_non_residual_defined_forcing": True,
        "complete_identity_bound_ns_defect": True,
        "matching_agent4_global_corrected_audit": True,
        "agent4_heldout_canonical_complete_ns_gate": True,
    }


def shortest_closure() -> list[str]:
    return [
        "A1: compute the exact current release2-end Q_s, then materialize the source-prescribed l=-h matching bridge until Q_s reaches the already-materialized Q_p before composing the terminal multiplier.",
        "A3/integration: make the exact A2 #1080 current-I4 runtime and A2 #960 differential runtime available unchanged/checksum-bound in the #1214 lineage, force ExactCurrentI4NonlinearBackend.bind, and record a real candidate execution of the #1210 RF30->RF39 path.",
        "A2/A3 integration: restack the new #1215 RF44 pre-state extension with the #1214 exact-backend-firewall lineage; sibling descendants of #1210 cannot be silently unioned. Then materialize RF44 post-update state/correction increment and Cartesian delta-u.",
        "A1/A2: finish terminal/exterior/global leading plus matched pressure and deterministic Python/MATLAB-safe export; only then build the matching self-contained global leading+oscillatory candidate.",
        "Integration/A4: bind preregistered restricted non-residual-defined forcing to the same corrected identity, run a real finite correction cycle, form the complete NS defect, and independently enforce held-out/canonical [24,48,96] momentum <=1e-3 and divergence <=1e-5 before pde_validated=true.",
    ]


def build_artifact() -> dict[str, Any]:
    artifact = {
        "schema": SCHEMA_VERSION,
        "task": TASK,
        "parent_a5": {"pr": PARENT_A5_PR, "exact_head": PARENT_A5_EXACT_HEAD},
        "sibling_topology": sibling_topology(),
        "upstream_routes": upstream_routes(),
        "latest_self_contained_candidate": latest_self_contained_candidate(),
        "core_states": core_states(),
        "fixed_validation_protocol": fixed_validation_protocol(),
        "evidence_firewall": evidence_firewall(),
        "remaining_blockers": remaining_blockers(),
        "shortest_closure": shortest_closure(),
        "truth_boundary": {
            "paper_exact": False,
            "openai_field_identified": False,
            "complete_ns_residual_assessed": False,
            "same_protocol_st006_comparison_available": False,
            "scientifically_admitted": False,
            "pde_validated": False,
        },
    }
    validate_artifact(artifact)
    return artifact


def validate_artifact(artifact: dict[str, Any]) -> None:
    if artifact.get("schema") != SCHEMA_VERSION:
        raise KokunoA5RoutingError("schema drift")
    if artifact.get("parent_a5") != {"pr": PARENT_A5_PR, "exact_head": PARENT_A5_EXACT_HEAD}:
        raise KokunoA5RoutingError("parent A5 identity drift")

    routes = artifact.get("upstream_routes", {})
    expected = {
        "agent1": (AGENT1_EXACT_HEAD, AGENT1_SOURCE_BLOB),
        "agent2": (AGENT2_EXACT_HEAD, AGENT2_SOURCE_BLOB),
        "agent3": (AGENT3_EXACT_HEAD, AGENT3_SOURCE_BLOB),
        "agent4": (AGENT4_EXACT_HEAD, AGENT4_SOURCE_BLOB),
    }
    if set(routes) != set(expected):
        raise KokunoA5RoutingError("upstream route set drift")
    for lane, (head, blob) in expected.items():
        route = routes[lane]
        if route.get("exact_head") != head or route.get("source_blob") != blob:
            raise KokunoA5RoutingError(f"{lane} identity/blob drift")

    topology = artifact.get("sibling_topology", {})
    if topology.get("common_parent_exact_head") != AGENT3_COMMON_PARENT_EXACT_HEAD:
        raise KokunoA5RoutingError("A2/A3 common parent drift")
    if not topology.get("agent2_and_agent3_are_sibling_descendants_of_1210"):
        raise KokunoA5RoutingError("A2/A3 sibling topology lost")
    if topology.get("agent2_1215_is_ancestor_of_agent3_1214") or topology.get("agent3_1214_is_ancestor_of_agent2_1215"):
        raise KokunoA5RoutingError("A2/A3 sibling branches falsely linearized")
    if topology.get("silent_cross_branch_evidence_union_allowed") or topology.get("unified_stack_containing_agent2_1215_and_agent3_1214_materialized"):
        raise KokunoA5RoutingError("cross-branch evidence silently unioned")

    a1, a2, a3, a4 = (routes[k] for k in ("agent1", "agent2", "agent3", "agent4"))
    if not a1.get("public_terminal_target_api_materialized"):
        raise KokunoA5RoutingError("A1 terminal target delivery lost")
    for key in (
        "current_release2_end_q_s_materialized",
        "l_minus_h_matching_bridge_materialized",
        "cartesian_terminal_multiplier_composed",
        "exterior_heat_replacement_materialized",
        "outer_global_leading_velocity_materialized",
        "matched_global_pressure_materialized",
        "cross_language_global_velocity_export_ready",
        "heldout_complete_ns_residual_assessed",
    ):
        if a1.get(key):
            raise KokunoA5RoutingError(f"A1 downstream stage falsely promoted: {key}")

    if a2.get("raw_auxiliary_t2_provider_source_blob") != AGENT2_PROVIDER_SOURCE_BLOB:
        raise KokunoA5RoutingError("A2 raw provider identity drift")
    if not a2.get("rf44_pre_update_state_materialized"):
        raise KokunoA5RoutingError("A2 RF44 pre-state delivery lost")
    for key in (
        "rf44_post_update_state_materialized",
        "rf44_correction_increment_materialized",
        "cartesian_correction_velocity_materialized",
        "finite_correction_cycle_run",
        "repository_candidate_scientific_evidence",
        "complete_ns_residual_assessed",
    ):
        if a2.get(key):
            raise KokunoA5RoutingError(f"A2 pre-state falsely promoted: {key}")

    if not a3.get("exact_runtime_rebind_required"):
        raise KokunoA5RoutingError("A3 exact-backend firewall lost")
    if a3.get("concrete_a2_1080_runtime_module_present_in_exact_stack") or a3.get("concrete_a2_960_differential_module_present_in_exact_stack"):
        raise KokunoA5RoutingError("A3 missing runtime dependency falsely materialized")
    for key in (
        "real_exact_backend_rebind_executed",
        "real_current_candidate_compact_correction_numerically_executed",
        "direct_constructed_backend_fixture_is_scientific_evidence",
        "parent_1210_promotion_without_rebind_is_sufficient_scientific_evidence",
        "repository_candidate_correction_scientifically_admitted",
        "cartesian_correction_velocity_materialized",
        "rf44_rf49_postupdate_remainder_materialized_on_real_candidate",
        "finite_correction_cycle_run",
        "heldout_complete_ns_residual_assessed",
    ):
        if a3.get(key):
            raise KokunoA5RoutingError(f"A3 firewall bypass/promotion: {key}")

    if a4.get("audits_agent1_pr") != AGENT1_PR or not a4.get("matching_terminal_target_audit_present"):
        raise KokunoA5RoutingError("A4 terminal-target identity drift")
    if a4.get("scoped_gate_passed") is not None:
        raise KokunoA5RoutingError("unresolved A4 scoped result pre-ingested")
    if a4.get("cartesian_leading_divergence_assessed_here") or a4.get("complete_ns_residual_assessed") or a4.get("scientifically_admitted") or a4.get("pde_validated"):
        raise KokunoA5RoutingError("A4 scoped target audit falsely promoted")

    candidate = artifact.get("latest_self_contained_candidate", {})
    if candidate.get("agent2_pr") != LATEST_SELF_CONTAINED_A2_PR or candidate.get("stage") != LATEST_SELF_CONTAINED_STAGE:
        raise KokunoA5RoutingError("latest self-contained candidate identity drift")
    for key in (
        "superseded_by_agent1_1213_as_unified_candidate",
        "superseded_by_agent2_1215_as_unified_candidate",
        "superseded_by_agent3_1214_as_corrected_candidate",
    ):
        if candidate.get(key):
            raise KokunoA5RoutingError("non-unified sibling falsely superseded project candidate")

    states = artifact.get("core_states", {})
    if states != {
        "leading_ready": False,
        "oscillatory_ready": True,
        "correction_ready": False,
        "velocity_export_ready": False,
        "pde_validated": False,
    }:
        raise KokunoA5RoutingError("core state promotion/drift")

    protocol = artifact.get("fixed_validation_protocol", {})
    if protocol.get("st006_momentum_sampled_max") != ST006_MOMENTUM_SAMPLED_MAX:
        raise KokunoA5RoutingError("ST006 max baseline drift")
    if protocol.get("st006_momentum_volume_l2") != ST006_MOMENTUM_VOLUME_L2:
        raise KokunoA5RoutingError("ST006 L2 baseline drift")
    if protocol.get("normalized_momentum_sampled_max_gate") != FINAL_NORMALIZED_MOMENTUM_GATE:
        raise KokunoA5RoutingError("momentum max gate changed")
    if protocol.get("normalized_momentum_volume_l2_gate") != FINAL_NORMALIZED_MOMENTUM_GATE:
        raise KokunoA5RoutingError("momentum L2 gate changed")
    if protocol.get("normalized_divergence_sampled_max_gate") != FINAL_NORMALIZED_DIVERGENCE_GATE:
        raise KokunoA5RoutingError("divergence max gate changed")
    if protocol.get("normalized_divergence_volume_l2_gate") != FINAL_NORMALIZED_DIVERGENCE_GATE:
        raise KokunoA5RoutingError("divergence L2 gate changed")
    if protocol.get("canonical_quadrature") != CANONICAL_QUADRATURE:
        raise KokunoA5RoutingError("canonical quadrature changed")
    if not protocol.get("residual_defined_free_forcing_forbidden") or not protocol.get("posthoc_threshold_relaxation_forbidden"):
        raise KokunoA5RoutingError("validation/forcing firewall weakened")

    firewall = artifact.get("evidence_firewall", {})
    if any(firewall.values()):
        raise KokunoA5RoutingError("one or more forbidden evidence transfers became true")

    blockers = artifact.get("remaining_blockers", {})
    for closed in (
        "a1_terminal_qp_target_algebra",
        "matching_a4_terminal_target_audit_surface",
        "raw_same_identity_auxiliary_t2_provider_missing",
        "agent3_checksum_pin_of_exact_agent2_1198_provider_blob",
        "rf44_pre_update_state_provider",
    ):
        if blockers.get(closed) is not False:
            raise KokunoA5RoutingError(f"closed structural blocker reopened: {closed}")
    for open_key in (
        "a1_current_release2_end_q_s",
        "a1_l_minus_h_matching_bridge_to_qp",
        "a1_cartesian_terminal_multiplier",
        "a1_exterior_heat_and_global_leading",
        "matched_cartesian_pressure_and_grad_p",
        "cross_language_global_velocity_export_representation",
        "matching_global_self_contained_a2_composite",
        "agent3_exact_a2_1080_runtime_dependency_in_stack",
        "agent3_exact_a2_960_differential_dependency_in_stack",
        "agent3_real_exact_backend_rebind",
        "a2_1215_a3_1214_unified_same_identity_stack",
        "rf44_post_update_state_and_increment_provider",
        "scientific_candidate_specific_rf30_to_rf49_correction",
        "cartesian_correction_velocity",
        "rf44_rf49_repository_candidate_remainder",
        "real_finite_correction_cycle",
        "preregistered_restricted_non_residual_defined_forcing",
        "complete_identity_bound_ns_defect",
        "matching_agent4_global_corrected_audit",
        "agent4_heldout_canonical_complete_ns_gate",
    ):
        if blockers.get(open_key) is not True:
            raise KokunoA5RoutingError(f"open blocker falsely closed: {open_key}")

    truth = artifact.get("truth_boundary", {})
    if any(truth.values()):
        raise KokunoA5RoutingError("truth boundary promoted")


def write_artifact(path: str | Path) -> Path:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(build_artifact(), indent=2, sort_keys=True) + "\n")
    return target


def _main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    write_artifact(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(_main())
