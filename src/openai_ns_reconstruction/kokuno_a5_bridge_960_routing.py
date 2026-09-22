"""Fail-closed Kokuno Agent-5 routing checkpoint for the bridge/#960 frontier.

This module is integration/provenance glue only.  It does not reimplement
Agent-1--4 mathematics and it cannot promote a scoped Kokuno replay into the
project's independent full-NS admission gate.
"""
from __future__ import annotations

import argparse
import json
from copy import deepcopy
from pathlib import Path
from typing import Any

SCHEMA_VERSION = "kokuno-a5-bridge-960-routing-v1"
TASK = "KOKUNO-A5-BRIDGE-960-ROUTING-130"

PARENT_A5_PR = 1228
PARENT_A5_EXACT_HEAD = "9784da1596a71d7d60dc4759f7571cb33f17d904"

AGENT1_PR = 1233
AGENT1_EXACT_HEAD = "6ae44fdac17483b1fd042dca6f5bdf06701ae068"
AGENT1_SOURCE_BLOB = "0db8baa2514ddf82b7f0a35c4fcb2fec66d21623"

AGENT4_PR = 1237
AGENT4_EXACT_HEAD = "ae26fda9c373f7cfb08a43a0ad48977ab13bcc71"
AGENT4_SOURCE_BLOB = "5dd00173615d1ee5c6b3e26cfb729196399c5667"

CORRECTION_COMMON_PARENT_PR = 1224
CORRECTION_COMMON_PARENT_HEAD = "107cac08050fddeab778416c1caf019be0f5795b"
AGENT2_PR = 1234
AGENT2_EXACT_HEAD = "f90bdb9c232fb16df881dfbcac28be1cd5c29780"
AGENT3_PR = 1236
AGENT3_EXACT_HEAD = "20120093abbd981e2d2cb5d95093cea1eebb28ef"

RF44_PRESTATE_BLOB = "0d8bb6471a4b7825b1b4b428d3b2c38d207729fd"
EXACT_BACKEND_FIREWALL_BLOB = "0152a85f7fbaf9c32c24c9d5c7d0d11166e3b157"
A2_960_DIFFERENTIAL_BLOB = "12df3bf6c949baeebf609b366f2973e7f515a157"
A2_960_AXIS_BLOB = "599baa190ec742d0e532e5517078534124e68d6b"
A2_960_VORTICITY_BLOB = "4ae525bad1e4b83c9dcabc1a97d3931562d099a5"
A2_1080_RUNTIME_PRESENT = False

LATEST_SELF_CONTAINED_A2_PR = 1117
LATEST_SELF_CONTAINED_A1_PR = 1107
LATEST_SELF_CONTAINED_STAGE = "xi=11"

ST006_MOMENTUM_SAMPLED_MAX = 0.1082289305112118
ST006_MOMENTUM_VOLUME_L2 = 0.10758432876230622
FINAL_NORMALIZED_MOMENTUM_GATE = 1.0e-3
FINAL_NORMALIZED_DIVERGENCE_GATE = 1.0e-5
CANONICAL_QUADRATURE = [24, 48, 96]


class KokunoA5RoutingError(RuntimeError):
    """Raised when routing, provenance, or scientific truth is promoted."""


def leading_route() -> dict[str, Any]:
    return {
        "agent1": {
            "pr": AGENT1_PR,
            "exact_head": AGENT1_EXACT_HEAD,
            "source_blob": AGENT1_SOURCE_BLOB,
            "delivery": "eta-resolved current inhomogeneous l=-h Q_s transport",
            "current_release2_qs_materialized": True,
            "current_l_minus_h_qs_transport_materialized": True,
            "pointwise_qp_first_hit_diagnostics_materialized": True,
            "full_eta_common_scalar_bridge_length_established": False,
            "eta_dependent_cartesian_bridge_geometry_materialized": False,
            "current_cartesian_terminal_multiplier_composed": False,
            "outer_global_leading_velocity_materialized": False,
            "matched_pressure_materialized": False,
        },
        "agent4": {
            "pr": AGENT4_PR,
            "exact_head": AGENT4_EXACT_HEAD,
            "source_blob": AGENT4_SOURCE_BLOB,
            "audits_agent1_pr": AGENT1_PR,
            "audits_agent1_exact_head": AGENT1_EXACT_HEAD,
            "matching_independent_audit_surface_present": True,
            "scientific_result_ingested": False,
            "scoped_gate_passed": None,
            "audit_is_complete_ns_validation": False,
            "pde_validated": False,
        },
        "representation_firewall": {
            "finite_eta_probes_define_one_scalar_bridge_length": False,
            "target_time_min_max_mean_rms_define_bridge_geometry": False,
            "representative_eta_defines_bridge_geometry": False,
            "source_ideal_matching_length_transfers_to_current_candidate": False,
            "pointwise_target_time_diagnostic_is_cartesian_terminal_surface": False,
        },
    }


def correction_routes() -> dict[str, Any]:
    shared_capability = {
        "rf44_prestate_blob": RF44_PRESTATE_BLOB,
        "exact_backend_firewall_blob": EXACT_BACKEND_FIREWALL_BLOB,
        "a2_960_differential_blob": A2_960_DIFFERENTIAL_BLOB,
        "a2_960_axis_blob": A2_960_AXIS_BLOB,
        "a2_960_vorticity_blob": A2_960_VORTICITY_BLOB,
        "rf44_prestate_present": True,
        "exact_backend_firewall_present": True,
        "exact_a2_960_runtime_present": True,
        "exact_a2_1080_runtime_present": A2_1080_RUNTIME_PRESENT,
        "real_exact_backend_bind_executed": False,
        "real_rf30_to_rf39_candidate_execution_recorded": False,
        "rf44_post_update_state_materialized": False,
        "cartesian_delta_u_materialized": False,
        "finite_correction_cycle_run": False,
    }
    return {
        "common_parent": {
            "pr": CORRECTION_COMMON_PARENT_PR,
            "exact_head": CORRECTION_COMMON_PARENT_HEAD,
        },
        "agent2_route": {
            "pr": AGENT2_PR,
            "exact_head": AGENT2_EXACT_HEAD,
            **shared_capability,
        },
        "agent3_route": {
            "pr": AGENT3_PR,
            "exact_head": AGENT3_EXACT_HEAD,
            **shared_capability,
        },
        "topology": {
            "agent2_1234_and_agent3_1236_are_diverged_siblings": True,
            "agent2_is_ancestor_of_agent3": False,
            "agent3_is_ancestor_of_agent2": False,
            "old_1224_1226_unified_runtime_seam_structurally_closed": True,
            "logical_capability_count": 1,
            "double_count_duplicate_routes": False,
            "silent_cross_branch_evidence_union_allowed": False,
            "downstream_must_consume_exactly_one_unified_route": True,
            "canonical_downstream_route_selected_here": False,
        },
    }


def latest_self_contained_candidate() -> dict[str, Any]:
    return {
        "agent2_pr": LATEST_SELF_CONTAINED_A2_PR,
        "agent1_pr": LATEST_SELF_CONTAINED_A1_PR,
        "stage": LATEST_SELF_CONTAINED_STAGE,
        "self_contained_leading_plus_oscillatory_velocity_xyzt": True,
        "superseded_by_agent1_1233": False,
        "superseded_by_agent2_1234": False,
        "superseded_by_agent3_1236": False,
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
        "kokuno_replay_is_independent_final_validation": False,
        "green_ci_is_pde_validation": False,
        "queued_or_running_ci_is_pass": False,
    }


def remaining_blockers() -> dict[str, bool]:
    """True means unresolved; False means structurally closed."""
    return {
        "actual_current_release2_qs": False,
        "current_l_minus_h_qs_transport": False,
        "matching_agent4_bridge_audit_surface": False,
        "agent4_bridge_scoped_result_resolved": True,
        "full_eta_bridge_geometry": True,
        "cartesian_terminal_multiplier": True,
        "exterior_global_leading": True,
        "matched_cartesian_pressure_and_grad_p": True,
        "deterministic_global_python_matlab_export": True,
        "rf44_prestate_firewall_a2_960_colocation": False,
        "duplicate_unified_960_route_deduplication": False,
        "exact_a2_1080_runtime_on_one_unified_route": True,
        "real_exact_backend_bind": True,
        "real_candidate_rf30_to_rf39_execution": True,
        "rf44_post_update_and_increment": True,
        "cartesian_correction_delta_u": True,
        "real_finite_correction_cycle": True,
        "matching_global_self_contained_a2_composite": True,
        "preregistered_restricted_non_residual_defined_forcing": True,
        "complete_identity_bound_ns_defect": True,
        "agent4_heldout_canonical_complete_ns_gate": True,
    }


def shortest_closure() -> list[str]:
    return [
        "A1/A4: resolve the exact #1237 scoped audit, then establish mathematically valid full-eta bridge geometry for exact #1233: one common scalar length only if justified over the full applicable eta interval, otherwise an explicit smooth eta-dependent Cartesian representation with a new save/load identity; only then compose the terminal stage.",
        "A3/integration: choose exactly one of the duplicate unified #960 routes (#1234 or #1236), restore/checksum-bind exact A2 #1080 and its dependency closure on that route, force ExactCurrentI4NonlinearBackend.bind, and record real RF30->RF39 candidate execution. Do not union sibling evidence or count the duplicated capability twice.",
        "A2/A3: materialize RF44 post-update/increment, Cartesian delta-u, nonlinear remainder, and a real finite correction cycle without held-out tuning.",
        "A1/A2: finish terminal/exterior/global leading, matched pressure/grad-p and deterministic Python/MATLAB export; then build a matching self-contained global leading+oscillatory candidate on the same semantic identity.",
        "Integration/A4: bind preregistered restricted non-residual-defined forcing, form the complete NS defect, and independently enforce held-out/canonical [24,48,96] momentum <=1e-3 and divergence <=1e-5 before any pde_validated=true promotion.",
    ]


def build_artifact() -> dict[str, Any]:
    artifact = {
        "schema": SCHEMA_VERSION,
        "task": TASK,
        "parent_a5": {"pr": PARENT_A5_PR, "exact_head": PARENT_A5_EXACT_HEAD},
        "leading_route": leading_route(),
        "correction_routes": correction_routes(),
        "latest_self_contained_candidate": latest_self_contained_candidate(),
        "core_states": core_states(),
        "fixed_validation_protocol": fixed_validation_protocol(),
        "remaining_blockers": remaining_blockers(),
        "shortest_closure": shortest_closure(),
        "truth_boundary": {
            "terminal_or_global_cartesian_velocity_materialized": False,
            "complete_candidate_api_materialized": False,
            "complete_ns_residual_assessed": False,
            "same_protocol_st006_improvement_assessed": False,
            "scientifically_admitted": False,
            "paper_exact": False,
            "openai_field_identified": False,
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
    if a.get("leading_route") != leading_route():
        raise KokunoA5RoutingError("leading route/truth drift")
    if a.get("correction_routes") != correction_routes():
        raise KokunoA5RoutingError("correction route/topology drift")
    if a.get("latest_self_contained_candidate") != latest_self_contained_candidate():
        raise KokunoA5RoutingError("self-contained candidate identity drift")
    if a.get("core_states") != core_states():
        raise KokunoA5RoutingError("readiness promotion/drift")
    if a.get("fixed_validation_protocol") != fixed_validation_protocol():
        raise KokunoA5RoutingError("validation protocol drift")
    if a.get("remaining_blockers") != remaining_blockers():
        raise KokunoA5RoutingError("blocker state drift")
    if a.get("shortest_closure") != shortest_closure():
        raise KokunoA5RoutingError("shortest-closure drift")

    expected_truth = {
        "terminal_or_global_cartesian_velocity_materialized": False,
        "complete_candidate_api_materialized": False,
        "complete_ns_residual_assessed": False,
        "same_protocol_st006_improvement_assessed": False,
        "scientifically_admitted": False,
        "paper_exact": False,
        "openai_field_identified": False,
        "pde_validated": False,
    }
    if a.get("truth_boundary") != expected_truth:
        raise KokunoA5RoutingError("scientific truth-boundary drift")


def assert_mutation_fails(path: list[str], value: Any) -> None:
    """Small utility used by regressions/CI to prove fail-closed mutations."""
    mutated = deepcopy(build_artifact())
    cursor: Any = mutated
    for key in path[:-1]:
        cursor = cursor[key]
    cursor[path[-1]] = value
    try:
        validate_artifact(mutated)
    except KokunoA5RoutingError:
        return
    raise AssertionError(f"mutation unexpectedly validated: {path}")


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
