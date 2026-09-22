"""Fail-closed Kokuno Agent-5 routing artifact for current-Qs/RF44/#960 seams.

Integration/provenance glue only.  It binds the fresh Agent-1--4 frontier,
records the two new branch/evidence separations, and keeps complete-NS/PDE
states fail-closed.  It does not reimplement upstream mathematics.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

SCHEMA_VERSION = "kokuno-a5-current-qs-rf44-960-routing-v1"
TASK = "KOKUNO-A5-CURRENT-QS-RF44-960-ROUTING-129"

PARENT_A5_PR = 1218
PARENT_A5_EXACT_HEAD = "01ed8bf5418166f1aaa0ef164381d9df2340743e"

A1_COMMON_PARENT_PR = 1217
A1_COMMON_PARENT_HEAD = "3978560078104bff7c9d5556bd9e658e4874c963"
AGENT1_PR = 1225
AGENT1_EXACT_HEAD = "8503f2ede2dd604392198b31ffdeb79321822a30"
AGENT1_SOURCE_BLOB = "7c6b7e6a95813db9422fff9194ea373faddc0895"

AGENT4_PR = 1223
AGENT4_EXACT_HEAD = "47d1f19bcea8cdbeaa4ae836eb191a76d22260ae"
AGENT4_SOURCE_BLOB = "8dbab20ec11dbd0048a5b5f4c7754b8a204ba535"

CORRECTION_COMMON_PARENT_PR = 1210
CORRECTION_COMMON_PARENT_HEAD = "e75940c33127b0725ddc703eaae97590158f2a56"
AGENT2_PR = 1224
AGENT2_EXACT_HEAD = "107cac08050fddeab778416c1caf019be0f5795b"
AGENT2_RF44_PRESTATE_BLOB = "0d8bb6471a4b7825b1b4b428d3b2c38d207729fd"
AGENT2_FIREWALL_BLOB = "0152a85f7fbaf9c32c24c9d5c7d0d11166e3b157"
AGENT2_RAW_AUX_T2_PROVIDER_BLOB = "96168ac6583ad5e71aa044bd558c6feeadf051f4"

AGENT3_PR = 1226
AGENT3_EXACT_HEAD = "571abc47c98001886cdcc6042832635e7bd93f31"
AGENT3_A2_960_DIFFERENTIAL_BLOB = "12df3bf6c949baeebf609b366f2973e7f515a157"
AGENT3_A2_960_AXIS_BLOB = "599baa190ec742d0e532e5517078534124e68d6b"
AGENT3_A2_960_VORTICITY_BLOB = "4ae525bad1e4b83c9dcabc1a97d3931562d099a5"

LATEST_SELF_CONTAINED_A2_PR = 1117
LATEST_SELF_CONTAINED_A1_PR = 1107
LATEST_SELF_CONTAINED_STAGE = "xi=11"

ST006_MOMENTUM_SAMPLED_MAX = 0.1082289305112118
ST006_MOMENTUM_VOLUME_L2 = 0.10758432876230622
FINAL_NORMALIZED_MOMENTUM_GATE = 1.0e-3
FINAL_NORMALIZED_DIVERGENCE_GATE = 1.0e-5
CANONICAL_QUADRATURE = [24, 48, 96]


class KokunoA5RoutingError(RuntimeError):
    """Raised when a routing/provenance or truth-boundary claim drifts."""


def branch_topology() -> dict[str, Any]:
    return {
        "leading_common_parent_pr": A1_COMMON_PARENT_PR,
        "leading_common_parent_head": A1_COMMON_PARENT_HEAD,
        "a1_1225_and_a4_1223_are_diverged_siblings": True,
        "a4_1223_is_ancestor_of_a1_1225": False,
        "a1_1225_is_ancestor_of_a4_1223": False,
        "correction_common_parent_pr": CORRECTION_COMMON_PARENT_PR,
        "correction_common_parent_head": CORRECTION_COMMON_PARENT_HEAD,
        "a2_1224_and_a3_1226_are_diverged_descendants": True,
        "a2_1224_is_ancestor_of_a3_1226": False,
        "a3_1226_is_ancestor_of_a2_1224": False,
        "silent_cross_branch_evidence_union_allowed": False,
        "single_stack_with_rf44_prestate_firewall_and_a2_960_runtime_materialized": False,
    }


def upstream_routes() -> dict[str, dict[str, Any]]:
    return {
        "agent1": {
            "pr": AGENT1_PR,
            "exact_head": AGENT1_EXACT_HEAD,
            "source_blob": AGENT1_SOURCE_BLOB,
            "delivery": "actual current release2-end Q_s from the general radial identity",
            "current_q_s_release2_endpoint_materialized": True,
            "source_ideal_q_s_relabelled_as_current": False,
            "current_l_minus_h_matching_bridge_materialized": False,
            "current_cartesian_terminal_multiplier_composed": False,
            "outer_global_leading_velocity_materialized": False,
            "matched_global_pressure_materialized": False,
            "complete_ns_residual_assessed": False,
        },
        "agent2": {
            "pr": AGENT2_PR,
            "exact_head": AGENT2_EXACT_HEAD,
            "rf44_prestate_blob": AGENT2_RF44_PRESTATE_BLOB,
            "firewall_blob": AGENT2_FIREWALL_BLOB,
            "raw_auxiliary_t2_provider_blob": AGENT2_RAW_AUX_T2_PROVIDER_BLOB,
            "delivery": "RF44 pre-state and exact-backend firewall co-resident in one ancestry",
            "rf44_pre_update_state_materialized": True,
            "exact_backend_firewall_present": True,
            "concrete_a2_960_runtime_present": False,
            "concrete_a2_1080_runtime_present": False,
            "real_exact_backend_bind_executed": False,
            "rf44_post_update_state_materialized": False,
            "cartesian_delta_u_materialized": False,
            "repository_candidate_correction_scientifically_admitted": False,
        },
        "agent3": {
            "pr": AGENT3_PR,
            "exact_head": AGENT3_EXACT_HEAD,
            "a2_960_differential_blob": AGENT3_A2_960_DIFFERENTIAL_BLOB,
            "a2_960_axis_blob": AGENT3_A2_960_AXIS_BLOB,
            "a2_960_vorticity_blob": AGENT3_A2_960_VORTICITY_BLOB,
            "delivery": "exact A2 #960 differential runtime restored on the #1214 lineage",
            "exact_backend_firewall_present": True,
            "concrete_a2_960_runtime_present": True,
            "concrete_a2_1080_runtime_present": False,
            "rf44_pre_update_state_from_agent2_1215_present": False,
            "real_exact_backend_bind_executed": False,
            "real_rf30_to_rf39_candidate_execution_recorded": False,
            "cartesian_delta_u_materialized": False,
            "finite_correction_cycle_run": False,
        },
        "agent4": {
            "pr": AGENT4_PR,
            "exact_head": AGENT4_EXACT_HEAD,
            "source_blob": AGENT4_SOURCE_BLOB,
            "delivery": "implementation-distinct audit of source-ideal exterior Q_s schedule only",
            "audits_agent1_pr": A1_COMMON_PARENT_PR,
            "source_ideal_q_s_schedule_audited": True,
            "current_candidate_q_s_audited": False,
            "matching_audit_for_agent1_1225_present": False,
            "complete_ns_residual_assessed": False,
            "scientifically_admitted": False,
            "pde_validated": False,
        },
    }


def latest_self_contained_candidate() -> dict[str, Any]:
    return {
        "agent2_pr": LATEST_SELF_CONTAINED_A2_PR,
        "agent1_pr": LATEST_SELF_CONTAINED_A1_PR,
        "stage": LATEST_SELF_CONTAINED_STAGE,
        "self_contained_leading_plus_oscillatory_velocity_xyzt": True,
        "superseded_by_agent1_1225": False,
        "superseded_by_agent2_1224": False,
        "superseded_by_agent3_1226": False,
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
        "green_ci_is_pde_validation": False,
    }


def evidence_firewall() -> dict[str, bool]:
    return {
        "source_ideal_a4_audit_transfers_to_current_a1_qs": False,
        "actual_current_qs_scalar_is_cartesian_terminal_velocity": False,
        "a2_1224_and_a3_1226_evidence_can_be_silently_unioned": False,
        "a2_1224_firewall_plus_prestate_implies_a2_960_runtime": False,
        "a3_1226_a2_960_runtime_implies_rf44_prestate": False,
        "a2_960_runtime_without_a2_1080_is_full_exact_backend": False,
        "scoped_qs_audit_is_complete_ns_validation": False,
        "kokuno_replay_is_independent_final_validation": False,
        "queued_or_running_ci_means_pass": False,
    }


def remaining_blockers() -> dict[str, bool]:
    """True means the blocker is still unresolved."""
    return {
        "a1_current_release2_end_q_s": False,
        "matching_a4_actual_current_qs_audit": True,
        "a1_current_l_minus_h_matching_bridge_to_qp": True,
        "a1_cartesian_terminal_multiplier": True,
        "a1_exterior_heat_and_global_leading": True,
        "matched_cartesian_pressure_and_grad_p": True,
        "matching_global_self_contained_a2_composite": True,
        "a2_rf44_prestate_plus_firewall_same_lineage": False,
        "a2_a3_1224_1226_unified_runtime_lineage": True,
        "exact_a2_960_runtime_available_somewhere": False,
        "exact_a2_1080_composite_runtime_in_unified_stack": True,
        "real_exact_backend_bind": True,
        "real_candidate_rf30_to_rf39_execution": True,
        "rf44_post_update_state_and_increment": True,
        "cartesian_correction_delta_u": True,
        "real_finite_correction_cycle": True,
        "preregistered_restricted_non_residual_defined_forcing": True,
        "complete_identity_bound_ns_defect": True,
        "agent4_heldout_canonical_complete_ns_gate": True,
    }


def shortest_closure() -> list[str]:
    return [
        "A1/A4: use the now-materialized actual release2-end Q_s to build the current l=-h bridge to Q_p, then obtain an implementation-distinct A4 audit of that exact current identity before Cartesian terminal promotion.",
        "A2/A3 integration: explicitly restack #1224 with #1226 so RF44 pre-state + exact-backend firewall + exact A2 #960 runtime coexist in one lineage; do not cross-union sibling evidence.",
        "A3: restore/checksum-bind the remaining exact A2 #1080 current-I4 composite runtime, force ExactCurrentI4NonlinearBackend.bind, and record real candidate RF30->RF39 execution.",
        "A2/A3: materialize RF44 post-update/increment, Cartesian delta-u, nonlinear remainder, and a real finite correction cycle.",
        "A1/A2/A4: finish terminal/exterior/global leading plus matched pressure and deterministic Python/MATLAB export; build the matching self-contained global leading+oscillatory candidate, bind preregistered restricted forcing, form complete NS defect, then independently enforce held-out/canonical [24,48,96] momentum <=1e-3 and divergence <=1e-5.",
    ]


def build_artifact() -> dict[str, Any]:
    artifact = {
        "schema": SCHEMA_VERSION,
        "task": TASK,
        "parent_a5": {"pr": PARENT_A5_PR, "exact_head": PARENT_A5_EXACT_HEAD},
        "branch_topology": branch_topology(),
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
            "same_protocol_st006_improvement_assessed": False,
            "scientifically_admitted": False,
            "pde_validated": False,
        },
    }
    validate_artifact(artifact)
    return artifact


def validate_artifact(a: dict[str, Any]) -> None:
    if a.get("schema") != SCHEMA_VERSION or a.get("task") != TASK:
        raise KokunoA5RoutingError("schema/task drift")
    if a.get("parent_a5") != {"pr": PARENT_A5_PR, "exact_head": PARENT_A5_EXACT_HEAD}:
        raise KokunoA5RoutingError("parent A5 identity drift")

    if a.get("upstream_routes") != upstream_routes():
        raise KokunoA5RoutingError("upstream route/truth drift")
    if a.get("branch_topology") != branch_topology():
        raise KokunoA5RoutingError("branch topology drift")
    if a.get("latest_self_contained_candidate") != latest_self_contained_candidate():
        raise KokunoA5RoutingError("candidate identity drift")
    if a.get("core_states") != core_states():
        raise KokunoA5RoutingError("core readiness promotion/drift")
    if a.get("fixed_validation_protocol") != fixed_validation_protocol():
        raise KokunoA5RoutingError("validation protocol drift")
    if a.get("evidence_firewall") != evidence_firewall():
        raise KokunoA5RoutingError("evidence firewall drift")
    if a.get("remaining_blockers") != remaining_blockers():
        raise KokunoA5RoutingError("blocker state drift")
    if a.get("shortest_closure") != shortest_closure():
        raise KokunoA5RoutingError("closure plan drift")

    expected_truth = {
        "paper_exact": False,
        "openai_field_identified": False,
        "complete_ns_residual_assessed": False,
        "same_protocol_st006_improvement_assessed": False,
        "scientifically_admitted": False,
        "pde_validated": False,
    }
    if a.get("truth_boundary") != expected_truth:
        raise KokunoA5RoutingError("scientific truth boundary drift")


def write_artifact(path: str | Path) -> Path:
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(build_artifact(), indent=2, sort_keys=True) + "\n")
    return out


def _main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    write_artifact(args.output)


if __name__ == "__main__":
    _main()
