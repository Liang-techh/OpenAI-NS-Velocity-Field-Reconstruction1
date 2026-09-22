"""Fail-closed Kokuno Agent-5 routing checkpoint for scalar-bridge/refutation and runtime-rebind frontier.

Integration/provenance glue only.  This module does not reimplement Agent-1--4
mathematics and cannot promote scoped Kokuno evidence into the project's
independent complete-Navier--Stokes admission gate.
"""
from __future__ import annotations

import argparse
from copy import deepcopy
import json
from pathlib import Path
from typing import Any

SCHEMA_VERSION = "kokuno-a5-scalar-refutation-runtime-rebind-routing-v1"
TASK = "KOKUNO-A5-SCALAR-REFUTATION-RUNTIME-REBIND-ROUTING-131"

PARENT_A5_PR = 1238
PARENT_A5_EXACT_HEAD = "8ed2641571bbce4340252e1809550959f57acbec"

AGENT1_PR = 1242
AGENT1_EXACT_HEAD = "7e987ff0a3bea743fe3f0fe9d02d5b3c8c593ea3"
AGENT1_SOURCE_BLOB = "409031a36425befefa7038ebda8ca848cdde239a"

AGENT4_PR = 1244
AGENT4_EXACT_HEAD = "d38aa3dacc205deacf785308b616c24ffb4eb135"
AGENT4_SOURCE_BLOB = "359d141d5406ec12e18315fc618ef8bb96cd696c"

SELECTED_CORRECTION_PARENT_PR = 1234
SELECTED_CORRECTION_PARENT_HEAD = "f90bdb9c232fb16df881dfbcac28be1cd5c29780"
SELECTED_AGENT2_PR = 1243
SELECTED_AGENT2_EXACT_HEAD = "3afc18d93334b6d371856cce28f90c6bff40bbc0"
SELECTED_AGENT2_SOURCE_BLOB = "4faeea46ac15b96fb030957e4f8be73c23612ff9"

DUPLICATE_AGENT3_PR = 1236
DUPLICATE_AGENT3_EXACT_HEAD = "20120093abbd981e2d2cb5d95093cea1eebb28ef"

A2_1080_HEAD = "c40d8ddecd2971544a6e07dab093436b423cf326"
A2_1080_COMPOSITE_BLOB = "2a0a5aa5966b02da856bdcf51940f3c186042802"
A1_1079_HEAD = "b06742ca6e189499192ede3cce40f62cdc1e35ca"
A1_1079_LEADING_BLOB = "6f04ce0a856b44430402576dad88438da90d1ebb"
A2_960_DIFFERENTIAL_BLOB = "12df3bf6c949baeebf609b366f2973e7f515a157"

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
            "delivery": "one-way common-scalar bridge refutation certificate",
            "finite_unique_root_witness_surface_present": True,
            "positive_witness_would_refute_one_common_scalar_bridge": True,
            "negative_witness_would_establish_common_scalar_bridge": False,
            "actual_refutation_result_ingested": False,
            "common_scalar_bridge_refuted": None,
            "full_eta_target_totality_established": False,
            "full_eta_root_uniqueness_established": False,
            "full_eta_transversality_established": False,
            "smooth_eta_target_time_map_established": False,
            "full_eta_common_scalar_bridge_established": False,
            "eta_dependent_cartesian_bridge_materialized": False,
            "terminal_cartesian_stage_composed": False,
            "global_leading_velocity_materialized": False,
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
            "finite_absence_of_refutation_establishes_scalar_bridge": False,
            "finite_refutation_establishes_smooth_eta_dependent_bridge": False,
            "finite_root_samples_establish_full_eta_totality": False,
            "finite_root_samples_establish_global_uniqueness": False,
            "mean_min_max_rms_target_time_defines_bridge_geometry": False,
            "representative_eta_defines_bridge_geometry": False,
            "terminal_composition_authorized_before_continuum_geometry": False,
        },
    }


def correction_route() -> dict[str, Any]:
    return {
        "selected_route": {
            "parent_pr": SELECTED_CORRECTION_PARENT_PR,
            "parent_exact_head": SELECTED_CORRECTION_PARENT_HEAD,
            "agent2_pr": SELECTED_AGENT2_PR,
            "exact_head": SELECTED_AGENT2_EXACT_HEAD,
            "source_blob": SELECTED_AGENT2_SOURCE_BLOB,
            "selected_as_single_downstream_unified_route": True,
            "exact_a2_1080_external_snapshot_declared": True,
            "a2_1080_head": A2_1080_HEAD,
            "a2_1080_composite_blob": A2_1080_COMPOSITE_BLOB,
            "exact_a1_1079_external_snapshot_declared": True,
            "a1_1079_head": A1_1079_HEAD,
            "a1_1079_leading_blob": A1_1079_LEADING_BLOB,
            "local_exact_a2_960_blob": A2_960_DIFFERENTIAL_BLOB,
            "runtime_rebind_contract_surface_present": True,
            "exact_head_ci_resolved": False,
            "runtime_rebind_result_ingested": False,
            "exact_backend_rebind_promoted": False,
            "rf30_to_rf39_candidate_execution_recorded": False,
            "rf44_post_update_state_materialized": False,
            "cartesian_delta_u_materialized": False,
            "finite_correction_cycle_run": False,
        },
        "duplicate_nonselected_route": {
            "agent3_pr": DUPLICATE_AGENT3_PR,
            "exact_head": DUPLICATE_AGENT3_EXACT_HEAD,
            "logical_capability_already_counted_on_selected_parent_route": True,
            "selected_for_downstream": False,
            "evidence_unioned_into_selected_route": False,
        },
        "topology_firewall": {
            "selected_route_is_agent2_1234_to_1243": True,
            "duplicate_agent3_1236_may_be_silently_unioned": False,
            "duplicate_capability_count": 1,
            "runtime_authentication_is_rf30_rf39_execution": False,
            "runtime_authentication_is_pde_validation": False,
        },
    }


def latest_self_contained_candidate() -> dict[str, Any]:
    return {
        "agent2_pr": LATEST_SELF_CONTAINED_A2_PR,
        "agent1_pr": LATEST_SELF_CONTAINED_A1_PR,
        "stage": LATEST_SELF_CONTAINED_STAGE,
        "self_contained_leading_plus_oscillatory_velocity_xyzt": True,
        "superseded_by_agent1_1242": False,
        "superseded_by_agent2_1243": False,
        "superseded_by_agent4_1244": False,
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
        "green_software_ci_is_pde_validation": False,
        "queued_or_running_ci_is_pass": False,
    }


def remaining_blockers() -> dict[str, bool]:
    """True means unresolved; False means the structural availability seam is closed."""
    return {
        "a1_common_scalar_refutation_surface": False,
        "matching_a4_refutation_audit_surface": False,
        "a4_refutation_audit_result_resolved": True,
        "full_eta_bridge_geometry": True,
        "terminal_cartesian_stage": True,
        "exterior_global_leading": True,
        "matched_cartesian_pressure_and_grad_p": True,
        "deterministic_global_python_matlab_export": True,
        "single_correction_route_selected": False,
        "exact_runtime_rebind_contract_surface": False,
        "exact_runtime_rebind_ci_and_result_resolved": True,
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
        "A1/A4: resolve exact #1244 against exact #1242. If the finite disjoint-root witness refutes one common scalar bridge, abandon scalarization and construct an explicit smooth eta-dependent Cartesian matching surface with a new save/load identity; if it does not refute, still establish full-eta target totality, intended-root uniqueness/transversality and continuum geometry before any scalar bridge or terminal composition.",
        "A2/A3: let exact #1243 runtime-rebind CI resolve on the selected #1234 route; only after the exact #1080/#1079/#960 bind is authenticated may Agent 3 execute the real RF30->RF39 correction on that same ancestry. Runtime authentication alone is not correction execution.",
        "A2/A3: materialize RF44 post-update/increment, Cartesian delta-u, nonlinear remainder and a real finite correction cycle without held-out tuning.",
        "A1/A2: after legal bridge geometry, finish terminal/exterior/global leading, matched pressure/grad-p and deterministic Python/MATLAB export, then build a matching self-contained global leading+oscillatory candidate on the same semantic identity.",
        "Integration/A4: bind preregistered restricted non-residual-defined forcing, form the complete NS defect, and independently enforce held-out/canonical [24,48,96] momentum <=1e-3 and divergence <=1e-5 before any pde_validated=true promotion.",
    ]


def build_artifact() -> dict[str, Any]:
    artifact = {
        "schema": SCHEMA_VERSION,
        "task": TASK,
        "parent_a5": {"pr": PARENT_A5_PR, "exact_head": PARENT_A5_EXACT_HEAD},
        "leading_route": leading_route(),
        "correction_route": correction_route(),
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
    expected = {
        "leading_route": leading_route(),
        "correction_route": correction_route(),
        "latest_self_contained_candidate": latest_self_contained_candidate(),
        "core_states": core_states(),
        "fixed_validation_protocol": fixed_validation_protocol(),
        "remaining_blockers": remaining_blockers(),
        "shortest_closure": shortest_closure(),
    }
    for key, value in expected.items():
        if a.get(key) != value:
            raise KokunoA5RoutingError(f"{key} drift")

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
