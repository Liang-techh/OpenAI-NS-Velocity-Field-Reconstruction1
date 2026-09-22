"""Fail-closed Kokuno Agent-5 routing artifact for the release2/AuxT2/RF49 seam.

This module is integration/provenance glue only.  It records the exact fresh
Agent-1--4 frontier and the shortest remaining closure.  It does not
reimplement any upstream mathematics and it deliberately refuses to convert
scoped numerical/CI evidence into complete-NS or PDE validation.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

SCHEMA_VERSION = "kokuno-a5-release2-auxt2-derivative-rf49-routing-v1"
TASK = "KOKUNO-A5-RELEASE2-AUXT2-DERIVATIVE-RF49-ROUTING-127"

PARENT_A5_PR = 1200
PARENT_A5_EXACT_HEAD = "51dab8d613283a339c13e8633a88d79b4c85f2a2"

AGENT1_PR = 1204
AGENT1_EXACT_HEAD = "18e75e6e2df43206147db34788c68ee6a7fa0f37"
AGENT1_SOURCE_BLOB = "101d8da6b3deff2322e5c2efaea278fefa7d19e3"
AGENT1_PATH = "src/openai_ns_reconstruction/kokuno_pa16_current_cartesian_postswirl_release2.py"

AGENT2_PR = 1205
AGENT2_EXACT_HEAD = "addc5dd3e6c61bc9d483e015071b4f8f4573324d"
AGENT2_SOURCE_BLOB = "33d3a37706a9436ce6f37f924332c7e07307e5c5"
AGENT2_PATH = "src/openai_ns_reconstruction/kokuno_current_i4_auxiliary_t2_derivative_diagnostics.py"

AGENT2_PROVIDER_PR = 1198
AGENT2_PROVIDER_EXACT_HEAD = "5765c2b3bab7482df514efa87f9b6bba48e04b7f"
AGENT2_PROVIDER_SOURCE_BLOB = "96168ac6583ad5e71aa044bd558c6feeadf051f4"

AGENT3_PR = 1197
AGENT3_EXACT_HEAD = "aa6b1a77bee404e1e852835fe048a7365b58f379"
AGENT3_SOURCE_BLOB = "8f42375022a9fcd0aeb24eb07542fc6fd298b62d"
AGENT3_PATH = "src/openai_ns_reconstruction/kokuno_current_i4_rf44_rf49_handoff.py"

AGENT4_PR = 1199
AGENT4_EXACT_HEAD = "184a628da361c9d233df5011afb7a68e2e776680"
AGENT4_SOURCE_BLOB = "773a9bc7009deff6218c06d46c45ee5455467786"
AGENT4_PATH = "src/openai_ns_reconstruction/kokuno_a4_lminus1_hold_independent_audit.py"

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


def upstream_routes() -> dict[str, dict[str, Any]]:
    return {
        "agent1": {
            "pr": AGENT1_PR,
            "exact_head": AGENT1_EXACT_HEAD,
            "source_blob": AGENT1_SOURCE_BLOB,
            "path": AGENT1_PATH,
            "delivery": "current Cartesian leading through the public post-swirl -1 -> -h release2 transition",
            "source_l_minus1_hold_materialized": True,
            "source_l_minus1_to_minus_h_transition_materialized": True,
            "current_cartesian_l_minus1_to_minus_h_composed": True,
            "current_autonomous_h": 0.005,
            "decimal_output_significant_digits": 96,
            "end_to_end_96digit_arithmetic_materialized": False,
            "source_terminal_multiplier_materialized": False,
            "source_exterior_heat_replacement_materialized": False,
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
            "delivery": "three-resolution independent slow-derivative audit of exact current-I4 raw auxiliary-T2 provider",
            "provider_parent_pr": AGENT2_PROVIDER_PR,
            "provider_parent_exact_head": AGENT2_PROVIDER_EXACT_HEAD,
            "provider_parent_source_blob": AGENT2_PROVIDER_SOURCE_BLOB,
            "current_i4_raw_auxiliary_t2_provider_materialized": True,
            "three_resolution_slow_derivative_diagnostic_materialized": True,
            "audit_relative_steps": [2.0**-10, 2.0**-11, 2.0**-12],
            "audit_recomputes_derivatives_from_raw_wave": True,
            "production_private_fd4_helper_used_by_audit_operator": False,
            "source_exact_auxiliary_phase_assignment_recovered": False,
            "source_exact_slow_derivatives_claimed": False,
            "agent3_provider_pin_modified_here": False,
            "provider_checksum_pinned_by_agent3": False,
            "rf30_repository_candidate_state_authorized": False,
            "complete_ns_residual_assessed": False,
            "exact_head_ci_conclusion_ingested": False,
        },
        "agent3": {
            "pr": AGENT3_PR,
            "exact_head": AGENT3_EXACT_HEAD,
            "source_blob": AGENT3_SOURCE_BLOB,
            "path": AGENT3_PATH,
            "delivery": "identity-preserving RF30/RF34-RF39 to RF44-RF49 recompute handoff mechanics",
            "typed_rf30_to_rf39_chain_available": True,
            "rf44_rf49_postupdate_handoff_executable": True,
            "pinned_same_identity_repository_provider_blob": None,
            "sees_agent2_1198_provider_as_pinned": False,
            "sees_agent2_1205_derivative_audit_as_scientific_authorization": False,
            "repository_candidate_rf30_evidence": False,
            "repository_candidate_rf44_rf49_evidence": False,
            "cartesian_correction_velocity_materialized": False,
            "finite_correction_cycle_run": False,
            "heldout_complete_ns_residual_assessed": False,
            "exact_head_ci_conclusion_ingested": False,
        },
        "agent4": {
            "pr": AGENT4_PR,
            "exact_head": AGENT4_EXACT_HEAD,
            "source_blob": AGENT4_SOURCE_BLOB,
            "path": AGENT4_PATH,
            "delivery": "implementation-distinct save/load FD4 divergence audit of exact A1 #1196 l=-1 hold",
            "audits_agent1_pr": 1196,
            "audits_latest_agent1_release2_pr": False,
            "matching_release2_audit_present": False,
            "scoped_only": True,
            "scoped_gate_passed": None,
            "canonical_absolute_fd_transfer_ready": False,
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
        "superseded_by_agent1_1204_as_unified_candidate": False,
        "superseded_by_agent2_1205_as_unified_candidate": False,
        "superseded_by_agent3_1197_as_corrected_candidate": False,
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
        "agent1_1204_leading_stage_plus_agent2_1205_provider_audit_form_unified_candidate": False,
        "agent2_1205_derivative_audit_authorizes_agent3_rf30": False,
        "agent2_1198_provider_existence_alone_authorizes_agent3_rf30": False,
        "agent3_1197_unpinned_handoff_creates_scientific_correction": False,
        "agent3_1197_mechanics_authorizes_cartesian_delta_u": False,
        "agent4_1199_scoped_audit_validates_agent1_1204": False,
        "agent4_1199_scoped_audit_validates_agent2_1205": False,
        "agent4_1199_scoped_audit_validates_agent3_1197": False,
        "agent4_1199_scoped_audit_validates_complete_ns": False,
        "queued_or_running_ci_means_pass": False,
        "kokuno_replay_is_independent_final_validation": False,
    }


def remaining_blockers() -> dict[str, bool]:
    """True means the blocker remains unresolved."""
    return {
        "a1_l_minus1_hold": False,
        "a1_minus1_to_minus_h_transition": False,
        "a1_terminal_multiplier": True,
        "a1_exterior_heat_and_global_leading": True,
        "matched_cartesian_pressure_and_grad_p": True,
        "cross_language_global_velocity_export_representation": True,
        "matching_global_self_contained_a2_composite": True,
        "raw_same_identity_auxiliary_t2_provider_missing": False,
        "a2_three_resolution_auxiliary_derivative_audit_missing": False,
        "agent3_checksum_pin_of_exact_agent2_1198_provider_blob": True,
        "scientific_candidate_specific_rf30_to_rf49_correction": True,
        "cartesian_correction_velocity": True,
        "rf44_rf49_repository_candidate_remainder": True,
        "preregistered_restricted_non_residual_defined_forcing": True,
        "complete_identity_bound_ns_defect": True,
        "real_finite_correction_cycle": True,
        "matching_agent4_release2_audit": True,
        "matching_agent4_global_corrected_audit": True,
        "agent4_heldout_canonical_complete_ns_gate": True,
    }


def shortest_closure() -> list[str]:
    return [
        "A3: independently checksum-pin exact A2 #1198 provider blob 96168ac6583ad5e71aa044bd558c6feeadf051f4 and replay the current-I4 RF30->RF49 chain; A2 #1205 is numerical derivative audit evidence, not authorization.",
        "A1: continue exact #1204 through the source terminal multiplier and exterior/global leading; materialize matched pressure and a deterministic Python/MATLAB-safe export representation.",
        "A3: after provider authorization, materialize the correction as Cartesian delta-u, recompute the nonlinear remainder, and run a real finite correction cycle without held-out tuning.",
        "A2: once global A1 exists, build the matching self-contained global leading+oscillatory candidate rather than transferring current-I4 evidence across identities.",
        "Integration/A4: bind matched pressure and preregistered restricted non-residual-defined forcing to one corrected semantic identity, form the complete NS defect, then independently enforce held-out/canonical [24,48,96] momentum <=1e-3 and divergence <=1e-5 before pde_validated=true.",
    ]


def build_artifact() -> dict[str, Any]:
    artifact = {
        "schema": SCHEMA_VERSION,
        "task": TASK,
        "parent_a5": {"pr": PARENT_A5_PR, "exact_head": PARENT_A5_EXACT_HEAD},
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

    a1, a2, a3, a4 = (routes[k] for k in ("agent1", "agent2", "agent3", "agent4"))
    if not a1.get("source_l_minus1_to_minus_h_transition_materialized"):
        raise KokunoA5RoutingError("A1 release2 delivery lost")
    if a1.get("source_terminal_multiplier_materialized") or a1.get("outer_global_leading_velocity_materialized"):
        raise KokunoA5RoutingError("A1 terminal/global stages falsely promoted")
    if a1.get("end_to_end_96digit_arithmetic_materialized"):
        raise KokunoA5RoutingError("A1 Decimal output mislabelled as full precision arithmetic")

    if a2.get("provider_parent_source_blob") != AGENT2_PROVIDER_SOURCE_BLOB:
        raise KokunoA5RoutingError("A2 provider parent blob drift")
    if not a2.get("three_resolution_slow_derivative_diagnostic_materialized"):
        raise KokunoA5RoutingError("A2 derivative audit delivery lost")
    if a2.get("production_private_fd4_helper_used_by_audit_operator"):
        raise KokunoA5RoutingError("A2 audit no longer implementation-distinct")
    if a2.get("provider_checksum_pinned_by_agent3") or a2.get("rf30_repository_candidate_state_authorized"):
        raise KokunoA5RoutingError("A2 audit/provider availability falsely promoted to A3 authorization")

    if not a3.get("rf44_rf49_postupdate_handoff_executable"):
        raise KokunoA5RoutingError("A3 RF44-RF49 handoff delivery lost")
    if a3.get("pinned_same_identity_repository_provider_blob") is not None:
        raise KokunoA5RoutingError("A3 provider pin falsely inferred")
    if a3.get("repository_candidate_rf30_evidence") or a3.get("repository_candidate_rf44_rf49_evidence"):
        raise KokunoA5RoutingError("A3 mechanics falsely promoted to repository-candidate evidence")
    if a3.get("cartesian_correction_velocity_materialized") or a3.get("finite_correction_cycle_run"):
        raise KokunoA5RoutingError("A3 correction/cycle falsely promoted")

    if a4.get("audits_latest_agent1_release2_pr") or a4.get("matching_release2_audit_present"):
        raise KokunoA5RoutingError("older A4 audit falsely transferred to A1 #1204")
    if a4.get("scoped_gate_passed") is not None:
        raise KokunoA5RoutingError("unresolved A4 CI result pre-ingested")
    if a4.get("complete_ns_residual_assessed") or a4.get("scientifically_admitted") or a4.get("pde_validated"):
        raise KokunoA5RoutingError("A4 scoped evidence falsely promoted")

    candidate = artifact.get("latest_self_contained_candidate", {})
    if candidate.get("agent2_pr") != LATEST_SELF_CONTAINED_A2_PR or candidate.get("stage") != LATEST_SELF_CONTAINED_STAGE:
        raise KokunoA5RoutingError("latest self-contained candidate identity drift")
    for key in (
        "superseded_by_agent1_1204_as_unified_candidate",
        "superseded_by_agent2_1205_as_unified_candidate",
        "superseded_by_agent3_1197_as_corrected_candidate",
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
    if protocol.get("normalized_momentum_sampled_max_gate") != FINAL_NORMALIZED_MOMENTUM_GATE:
        raise KokunoA5RoutingError("momentum gate changed")
    if protocol.get("normalized_momentum_volume_l2_gate") != FINAL_NORMALIZED_MOMENTUM_GATE:
        raise KokunoA5RoutingError("momentum L2 gate changed")
    if protocol.get("normalized_divergence_sampled_max_gate") != FINAL_NORMALIZED_DIVERGENCE_GATE:
        raise KokunoA5RoutingError("divergence gate changed")
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
        "a1_l_minus1_hold",
        "a1_minus1_to_minus_h_transition",
        "raw_same_identity_auxiliary_t2_provider_missing",
        "a2_three_resolution_auxiliary_derivative_audit_missing",
    ):
        if blockers.get(closed) is not False:
            raise KokunoA5RoutingError(f"closed structural blocker reopened: {closed}")
    for open_key in (
        "a1_terminal_multiplier",
        "a1_exterior_heat_and_global_leading",
        "agent3_checksum_pin_of_exact_agent2_1198_provider_blob",
        "cartesian_correction_velocity",
        "complete_identity_bound_ns_defect",
        "real_finite_correction_cycle",
        "agent4_heldout_canonical_complete_ns_gate",
    ):
        if blockers.get(open_key) is not True:
            raise KokunoA5RoutingError(f"remaining blocker falsely closed: {open_key}")

    truth = artifact.get("truth_boundary", {})
    if any(bool(truth.get(k)) for k in (
        "paper_exact",
        "openai_field_identified",
        "complete_ns_residual_assessed",
        "same_protocol_st006_comparison_available",
        "scientifically_admitted",
        "pde_validated",
    )):
        raise KokunoA5RoutingError("scientific truth boundary falsely promoted")


def write_artifact(path: str | Path) -> Path:
    artifact = build_artifact()
    out = Path(path)
    out.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", default="kokuno_a5_release2_auxt2_derivative_rf49_routing.json")
    args = parser.parse_args()
    path = write_artifact(args.output)
    print(path)
    print(json.dumps(build_artifact(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
