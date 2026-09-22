"""Kokuno Agent-5 routing artifact for the Decimal-total / oriented / RF30 seam.

This module is deliberately integration-only.  It records immutable upstream
identities and fail-closed truth boundaries; it does not reimplement Agent-1--4
mathematics and it cannot promote scoped CI or replay evidence to PDE validity.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

SCHEMA_VERSION = "kokuno-a5-decimal-total-oriented-rf30-routing-v1"
TASK = "KOKUNO-A5-DECIMAL-TOTAL-ORIENTED-RF30-ROUTING-124"

PARENT_A5_PR = 1173
PARENT_A5_EXACT_HEAD = "ae5b344ecdedbd532a2d6739e69d12460d76616b"

AGENT1_PR = 1179
AGENT1_EXACT_HEAD = "8c5c5b6d55a68de8285ed1dd0ebb55128f3078a4"
AGENT1_PATH = "src/openai_ns_reconstruction/kokuno_pa16_current_cartesian_relative_swirl_decimal_composed.py"

AGENT2_PR = 1180
AGENT2_EXACT_HEAD = "3bdf0d6a33a629c4262875b146e1378d550023d0"
AGENT2_PATH = "src/openai_ns_reconstruction/kokuno_source_oriented_multiharmonic_diagnostics.py"

AGENT3_PR = 1181
AGENT3_EXACT_HEAD = "d157d2c718d7cb27935abb86229f7a2f0e8fbda3"
AGENT3_PATH = "src/openai_ns_reconstruction/kokuno_current_i4_rf30_typed_defect.py"

AGENT4_PR = 1182
AGENT4_EXACT_HEAD = "5e12ed7dc55e0740a40f9a210e23e98cd71707e2"
AGENT4_PATH = "src/openai_ns_reconstruction/kokuno_a4_relative_swirl_decimal_total_independent_audit.py"

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
    """Raised when the integration artifact attempts an invalid promotion."""


def upstream_routes() -> dict[str, dict[str, Any]]:
    """Return the exact upstream identities consumed only as routing evidence."""
    return {
        "agent1": {
            "pr": AGENT1_PR,
            "exact_head": AGENT1_EXACT_HEAD,
            "path": AGENT1_PATH,
            "delivery": "precision-qualified Decimal total relative-swirl composition",
            "fixed_decimal_significant_digits": 96,
            "current_cartesian_relative_swirl_composed_precision_qualified": True,
            "binary64_total_relative_swirl_sum_is_resolved": False,
            "terminal_exterior_global_leading_materialized": False,
            "matched_pressure_materialized": False,
        },
        "agent2": {
            "pr": AGENT2_PR,
            "exact_head": AGENT2_EXACT_HEAD,
            "path": AGENT2_PATH,
            "delivery": "three-resolution oriented complete-curl public-velocity diagnostic",
            "provider_driven": True,
            "self_contained_velocity_xyzt_provider": False,
            "source_exact_orientation_recovered": False,
            "complete_ns_residual_assessed": False,
        },
        "agent3": {
            "pr": AGENT3_PR,
            "exact_head": AGENT3_EXACT_HEAD,
            "path": AGENT3_PATH,
            "delivery": "typed current-I4 RF30 defect mechanics from auxiliary-Haar state",
            "typed_rf30_mechanics_materialized": True,
            "pinned_same_identity_repository_provider_blob": None,
            "repository_candidate_rf30_defect_scientifically_admitted": False,
            "rf31_five_row_system_materialized": False,
            "correction_applied": False,
            "heldout_ns_residual_assessed": False,
        },
        "agent4": {
            "pr": AGENT4_PR,
            "exact_head": AGENT4_EXACT_HEAD,
            "path": AGENT4_PATH,
            "audits_agent1_exact_head": AGENT1_EXACT_HEAD,
            "delivery": "implementation-distinct Decimal-total representation/divergence audit",
            "save_load_public_velocity_only_scientific_path": True,
            "scoped_only": True,
            "scoped_gate_passed": None,
            "canonical_absolute_fd_transfer_ready": False,
            "complete_ns_residual_assessed": False,
            "scientifically_admitted": False,
            "pde_validated": False,
        },
    }


def latest_self_contained_candidate() -> dict[str, Any]:
    """Keep newer provider/representation lanes separate from the last unified candidate."""
    return {
        "agent2_pr": LATEST_SELF_CONTAINED_A2_PR,
        "agent2_exact_head": LATEST_SELF_CONTAINED_A2_HEAD,
        "agent1_pr": LATEST_SELF_CONTAINED_A1_PR,
        "stage": LATEST_SELF_CONTAINED_STAGE,
        "self_contained_leading_plus_oscillatory_velocity_xyzt": True,
        "superseded_by_agent1_1179_as_unified_candidate": False,
        "superseded_by_agent2_1180_as_unified_candidate": False,
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
        "visual_success_is_pde_validation": False,
        "ci_success_is_scientific_admission": False,
    }


def evidence_firewall() -> dict[str, bool]:
    return {
        "agent4_1182_may_validate_agent2_1180": False,
        "agent4_1182_may_validate_agent3_1181": False,
        "agent4_1182_may_validate_xi11_composite_complete_ns": False,
        "agent4_1182_may_set_pde_validated": False,
        "agent2_1180_diagnostic_creates_self_contained_candidate": False,
        "agent3_1181_unpinned_provider_creates_scientific_rf30_defect": False,
        "agent3_1181_mechanics_authorizes_applied_correction": False,
        "queued_or_running_ci_means_pass": False,
    }


def remaining_blockers() -> dict[str, bool]:
    """True means blocker remains unresolved."""
    return {
        "terminal_exterior_global_leading": True,
        "cross_language_global_velocity_export_representation": True,
        "matching_global_self_contained_a2_composite": True,
        "checksum_pinned_same_identity_raw_auxiliary_t2_provider": True,
        "scientific_candidate_specific_rf30_defect": True,
        "rf31_rf34_rf49_applied_correction_and_remainder": True,
        "matched_cartesian_pressure_and_grad_p": True,
        "preregistered_restricted_non_residual_defined_forcing": True,
        "complete_identity_bound_ns_defect": True,
        "real_finite_correction_cycle": True,
        "agent4_heldout_canonical_complete_ns_gate": True,
    }


def shortest_closure() -> list[str]:
    return [
        "A1: extend the exact precision-qualified Decimal total identity through terminal/exterior/global leading and matched pressure while preserving a deterministic cross-language export representation.",
        "A2: checksum-pin a same-identity raw auxiliary-T2 provider for A3 and, after global A1 exists, build the matching self-contained global leading+oscillatory candidate.",
        "A3: with the pinned provider, materialize scientific candidate-specific RF30, then RF31/RF34-RF49 applied correction and recomputed nonlinear remainder.",
        "Integration: bind preregistered restricted non-residual-defined forcing and matched pressure to the same candidate identity, form the complete NS defect, and run the real finite correction cycle.",
        "A4: independently validate the final held-out complete NS candidate and enforce the unchanged canonical [24,48,96] momentum <=1e-3 and divergence <=1e-5 gates.",
    ]


def build_artifact() -> dict[str, Any]:
    artifact: dict[str, Any] = {
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
    parent = artifact.get("parent_a5", {})
    if parent != {"pr": PARENT_A5_PR, "exact_head": PARENT_A5_EXACT_HEAD}:
        raise KokunoA5RoutingError("parent A5 identity drift")

    routes = artifact.get("upstream_routes", {})
    expected_heads = {
        "agent1": AGENT1_EXACT_HEAD,
        "agent2": AGENT2_EXACT_HEAD,
        "agent3": AGENT3_EXACT_HEAD,
        "agent4": AGENT4_EXACT_HEAD,
    }
    if set(routes) != set(expected_heads):
        raise KokunoA5RoutingError("upstream route set drift")
    for lane, head in expected_heads.items():
        if routes[lane].get("exact_head") != head:
            raise KokunoA5RoutingError(f"{lane} exact-head drift")

    a1, a2, a3, a4 = (routes[k] for k in ("agent1", "agent2", "agent3", "agent4"))
    if not a1.get("current_cartesian_relative_swirl_composed_precision_qualified"):
        raise KokunoA5RoutingError("A1 Decimal composition was lost")
    if a1.get("binary64_total_relative_swirl_sum_is_resolved"):
        raise KokunoA5RoutingError("binary64 total was falsely promoted")
    if not a2.get("provider_driven") or a2.get("self_contained_velocity_xyzt_provider"):
        raise KokunoA5RoutingError("A2 provider truth boundary drift")
    if not a3.get("typed_rf30_mechanics_materialized"):
        raise KokunoA5RoutingError("A3 typed RF30 mechanics missing")
    if a3.get("pinned_same_identity_repository_provider_blob") is not None:
        raise KokunoA5RoutingError("A3 provider was promoted without a pinned sibling")
    if a3.get("repository_candidate_rf30_defect_scientifically_admitted") or a3.get("correction_applied"):
        raise KokunoA5RoutingError("A3 scientific correction was falsely promoted")
    if a4.get("audits_agent1_exact_head") != AGENT1_EXACT_HEAD:
        raise KokunoA5RoutingError("A4 audit is detached from exact A1 identity")
    if a4.get("scoped_gate_passed") is not None:
        raise KokunoA5RoutingError("A4 scoped numerical outcome is not ingested by this routing artifact")
    if not a4.get("scoped_only") or a4.get("scientifically_admitted") or a4.get("pde_validated"):
        raise KokunoA5RoutingError("A4 scoped evidence was over-promoted")

    candidate = artifact.get("latest_self_contained_candidate", {})
    if candidate.get("agent2_pr") != LATEST_SELF_CONTAINED_A2_PR or candidate.get("agent1_pr") != LATEST_SELF_CONTAINED_A1_PR:
        raise KokunoA5RoutingError("self-contained candidate identity drift")
    if candidate.get("superseded_by_agent1_1179_as_unified_candidate") or candidate.get("superseded_by_agent2_1180_as_unified_candidate"):
        raise KokunoA5RoutingError("non-unified sibling was promoted to unified candidate")

    if artifact.get("core_states") != {
        "leading_ready": False,
        "oscillatory_ready": True,
        "correction_ready": False,
        "velocity_export_ready": False,
        "pde_validated": False,
    }:
        raise KokunoA5RoutingError("core-state promotion detected")

    protocol = artifact.get("fixed_validation_protocol", {})
    if protocol.get("st006_momentum_sampled_max") != ST006_MOMENTUM_SAMPLED_MAX:
        raise KokunoA5RoutingError("ST006 max drift")
    if protocol.get("st006_momentum_volume_l2") != ST006_MOMENTUM_VOLUME_L2:
        raise KokunoA5RoutingError("ST006 L2 drift")
    if protocol.get("normalized_momentum_sampled_max_gate") != 1.0e-3 or protocol.get("normalized_momentum_volume_l2_gate") != 1.0e-3:
        raise KokunoA5RoutingError("momentum gate drift")
    if protocol.get("normalized_divergence_sampled_max_gate") != 1.0e-5 or protocol.get("normalized_divergence_volume_l2_gate") != 1.0e-5:
        raise KokunoA5RoutingError("divergence gate drift")
    if protocol.get("canonical_quadrature") != [24, 48, 96]:
        raise KokunoA5RoutingError("canonical quadrature drift")
    if not protocol.get("residual_defined_free_forcing_forbidden"):
        raise KokunoA5RoutingError("free-forcing prohibition weakened")

    firewall = artifact.get("evidence_firewall", {})
    if any(firewall.values()):
        raise KokunoA5RoutingError("cross-identity/scoped evidence promotion detected")
    if not all(artifact.get("remaining_blockers", {}).values()):
        raise KokunoA5RoutingError("blocker was cleared without same-identity evidence")
    truth = artifact.get("truth_boundary", {})
    if any(truth.values()):
        raise KokunoA5RoutingError("scientific truth boundary was promoted")


def dump_artifact(path: str | Path) -> Path:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(build_artifact(), sort_keys=True, indent=2) + "\n", encoding="utf-8")
    return target


__all__ = [
    "SCHEMA_VERSION",
    "TASK",
    "PARENT_A5_EXACT_HEAD",
    "build_artifact",
    "validate_artifact",
    "dump_artifact",
]
