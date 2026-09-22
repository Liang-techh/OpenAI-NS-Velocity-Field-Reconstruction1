"""Fail-closed Kokuno Agent-5 routing artifact for the release1/scale/RF39 seam.

Integration only: this module records exact sibling identities, construction
frontiers, the frozen validation contract, and the shortest remaining closure.
It does not reimplement Agent-1--4 mathematics and cannot promote scoped replay,
CI, or provider-driven mechanics to PDE validation.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

SCHEMA_VERSION = "kokuno-a5-release1-scale-rf39-routing-v1"
TASK = "KOKUNO-A5-RELEASE1-SCALE-RF39-ROUTING-125"

PARENT_A5_PR = 1183
PARENT_A5_EXACT_HEAD = "9f4bc6b415e2e0334ccdc1b681b83539972cee52"

AGENT1_PR = 1188
AGENT1_EXACT_HEAD = "f5b5b8b532e0802e181bfead5e748cd4008906bf"
AGENT1_SOURCE_BLOB = "98be6d0d31a382ad7638318cec670d317bcc5cce"
AGENT1_PATH = "src/openai_ns_reconstruction/kokuno_pa16_current_cartesian_postswirl_release1.py"

AGENT2_PR = 1189
AGENT2_EXACT_HEAD = "7e9d8ce853f43f96eecde9343dcbe0c3138a6042"
AGENT2_SOURCE_BLOB = "26b59059b2341f7a6a4d1f985fb719df007b4d19"
AGENT2_PATH = "src/openai_ns_reconstruction/kokuno_source_bounded_scale_oriented_multiharmonic_velocity.py"

AGENT3_PR = 1190
AGENT3_EXACT_HEAD = "21b0ade966fda3b7344273c7bff2efe7e7b9a93e"
AGENT3_SOURCE_BLOB = "fe08f93f48aeb2d4055b12889325a6ae5e05a5a6"
AGENT3_PATH = "src/openai_ns_reconstruction/kokuno_current_i4_rf34_rf39_correction.py"

AGENT4_PR = 1182
AGENT4_EXACT_HEAD = "5e12ed7dc55e0740a40f9a210e23e98cd71707e2"
AGENT4_SOURCE_BLOB = "9fbc39f4e9740a504e277a9cb758a144eb219c62"
AGENT4_PATH = "src/openai_ns_reconstruction/kokuno_a4_relative_swirl_decimal_total_independent_audit.py"
AGENT4_AUDITED_A1_PR = 1179
AGENT4_AUDITED_A1_HEAD = "8c5c5b6d55a68de8285ed1dd0ebb55128f3078a4"

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
    """Raised when routing/provenance or truth-state promotion fails closed."""


def upstream_routes() -> dict[str, dict[str, Any]]:
    return {
        "agent1": {
            "pr": AGENT1_PR,
            "exact_head": AGENT1_EXACT_HEAD,
            "source_blob": AGENT1_SOURCE_BLOB,
            "path": AGENT1_PATH,
            "delivery": "first post-relative-swirl release interval with Decimal Cartesian velocity",
            "post_relative_swirl_release1_materialized": True,
            "fixed_decimal_significant_digits": 96,
            "source_l_minus1_hold_materialized": False,
            "source_l_minus1_to_minus_h_transition_materialized": False,
            "source_terminal_multiplier_materialized": False,
            "source_exterior_heat_replacement_materialized": False,
            "outer_global_leading_velocity_materialized": False,
            "matched_global_pressure_materialized": False,
            "cross_language_global_velocity_export_ready": False,
            "heldout_ns_residual_assessed": False,
        },
        "agent2": {
            "pr": AGENT2_PR,
            "exact_head": AGENT2_EXACT_HEAD,
            "source_blob": AGENT2_SOURCE_BLOB,
            "path": AGENT2_PATH,
            "delivery": "bounded isotropic physical scale on oriented complete-curl family",
            "bounded_isotropic_scale_materialized": True,
            "scale_bounds": [0.5, 2.0],
            "provider_driven": True,
            "self_contained_velocity_xyzt_provider": False,
            "source_exact_scale_recovered": False,
            "same_identity_auxiliary_t2_provider_checksum_pinned": False,
            "complete_ns_residual_assessed": False,
        },
        "agent3": {
            "pr": AGENT3_PR,
            "exact_head": AGENT3_EXACT_HEAD,
            "source_blob": AGENT3_SOURCE_BLOB,
            "path": AGENT3_PATH,
            "delivery": "typed current-I4 RF30 -> RF31/RF34-RF39 mechanics adapter",
            "typed_rf30_to_rf34_rf39_chain_executable": True,
            "rf31_basis_frozen_before_defect": True,
            "pinned_same_identity_repository_provider_blob": None,
            "repository_candidate_scientific_correction_materialized": False,
            "correction_applied_to_candidate": False,
            "cartesian_correction_velocity_materialized": False,
            "rf44_rf49_nonlinear_remainder_recomputed": False,
            "finite_correction_cycle_run": False,
            "heldout_ns_residual_assessed": False,
        },
        "agent4": {
            "pr": AGENT4_PR,
            "exact_head": AGENT4_EXACT_HEAD,
            "source_blob": AGENT4_SOURCE_BLOB,
            "path": AGENT4_PATH,
            "audits_agent1_pr": AGENT4_AUDITED_A1_PR,
            "audits_agent1_exact_head": AGENT4_AUDITED_A1_HEAD,
            "delivery": "independent Decimal-total representation/divergence audit of A1 #1179",
            "save_load_public_velocity_only_scientific_path": True,
            "scoped_only": True,
            "scoped_gate_passed": None,
            "canonical_absolute_fd_transfer_ready": False,
            "audits_latest_agent1_release1": False,
            "audits_agent2_1189": False,
            "audits_agent3_1190": False,
            "complete_ns_residual_assessed": False,
            "scientifically_admitted": False,
            "pde_validated": False,
        },
    }


def latest_self_contained_candidate() -> dict[str, Any]:
    return {
        "agent2_pr": LATEST_SELF_CONTAINED_A2_PR,
        "agent2_exact_head": LATEST_SELF_CONTAINED_A2_HEAD,
        "agent1_pr": LATEST_SELF_CONTAINED_A1_PR,
        "stage": LATEST_SELF_CONTAINED_STAGE,
        "self_contained_leading_plus_oscillatory_velocity_xyzt": True,
        "superseded_by_agent1_1188_as_unified_candidate": False,
        "superseded_by_agent2_1189_as_unified_candidate": False,
        "superseded_by_agent3_1190_as_corrected_candidate": False,
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
        "agent4_1182_may_validate_agent1_1188": False,
        "agent4_1182_may_validate_agent2_1189": False,
        "agent4_1182_may_validate_agent3_1190": False,
        "agent4_1182_may_validate_complete_ns": False,
        "agent4_1182_may_set_pde_validated": False,
        "agent2_1189_provider_family_creates_self_contained_candidate": False,
        "agent3_1190_unpinned_provider_creates_scientific_correction": False,
        "agent3_1190_mechanics_authorizes_cartesian_delta_u": False,
        "queued_or_running_ci_means_pass": False,
    }


def remaining_blockers() -> dict[str, bool]:
    """True means the blocker remains unresolved."""
    return {
        "a1_l_minus1_hold_and_minus1_to_minus_h_transition": True,
        "a1_terminal_multiplier_exterior_heat_and_global_leading": True,
        "cross_language_global_velocity_export_representation": True,
        "matched_cartesian_pressure_and_grad_p": True,
        "matching_global_self_contained_a2_composite": True,
        "checksum_pinned_same_identity_raw_auxiliary_t2_provider": True,
        "scientific_candidate_specific_rf30_to_rf39_correction": True,
        "cartesian_correction_velocity": True,
        "rf44_rf49_recomputed_remainder": True,
        "preregistered_restricted_non_residual_defined_forcing": True,
        "complete_identity_bound_ns_defect": True,
        "real_finite_correction_cycle": True,
        "matching_agent4_latest_global_composite_audit": True,
        "agent4_heldout_canonical_complete_ns_gate": True,
    }


def shortest_closure() -> list[str]:
    return [
        "A1: continue exact #1188 through the l=-1 hold, -1->-h transition, terminal multiplier and exterior/global leading while preserving the precision-qualified representation; materialize matched pressure and a deterministic Python/MATLAB export representation.",
        "A2: provide and checksum-pin a same-identity raw auxiliary-T2 provider for A3; after global A1 exists, build the matching self-contained global leading+oscillatory candidate using the bounded complete-curl family.",
        "A3: consume that pinned provider through the already-executable #1190 RF30->RF39 chain, apply the correction as Cartesian delta-u, and recompute RF44-RF49/nonlinear remainder under one candidate identity.",
        "Integration: bind matched pressure and preregistered restricted non-residual-defined forcing to that same identity, form the complete NS defect, and run the real finite correction cycle.",
        "A4: independently audit the matching global/corrected identity and enforce the unchanged held-out canonical [24,48,96] momentum <=1e-3 and divergence <=1e-5 gates before any pde_validated promotion.",
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
        if routes[lane].get("exact_head") != head or routes[lane].get("source_blob") != blob:
            raise KokunoA5RoutingError(f"{lane} identity/blob drift")

    a1, a2, a3, a4 = (routes[k] for k in ("agent1", "agent2", "agent3", "agent4"))
    if not a1.get("post_relative_swirl_release1_materialized") or a1.get("outer_global_leading_velocity_materialized"):
        raise KokunoA5RoutingError("A1 release/global boundary drift")
    if a1.get("source_l_minus1_hold_materialized") or a1.get("source_exterior_heat_replacement_materialized"):
        raise KokunoA5RoutingError("A1 later stages falsely promoted")
    if not a2.get("bounded_isotropic_scale_materialized") or not a2.get("provider_driven"):
        raise KokunoA5RoutingError("A2 bounded-scale delivery lost")
    if a2.get("self_contained_velocity_xyzt_provider") or a2.get("same_identity_auxiliary_t2_provider_checksum_pinned"):
        raise KokunoA5RoutingError("A2 provider state falsely promoted")
    if not a3.get("typed_rf30_to_rf34_rf39_chain_executable"):
        raise KokunoA5RoutingError("A3 RF30->RF39 mechanics lost")
    if a3.get("pinned_same_identity_repository_provider_blob") is not None:
        raise KokunoA5RoutingError("A3 provider falsely pinned")
    if a3.get("repository_candidate_scientific_correction_materialized") or a3.get("correction_applied_to_candidate"):
        raise KokunoA5RoutingError("A3 scientific correction falsely promoted")
    if a4.get("audits_agent1_exact_head") != AGENT4_AUDITED_A1_HEAD or a4.get("audits_latest_agent1_release1"):
        raise KokunoA5RoutingError("A4 evidence transferred to latest A1 identity")
    if a4.get("scoped_gate_passed") is not None or a4.get("scientifically_admitted") or a4.get("pde_validated"):
        raise KokunoA5RoutingError("A4 scoped evidence over-promoted")

    candidate = artifact.get("latest_self_contained_candidate", {})
    if candidate.get("agent2_pr") != 1117 or candidate.get("agent1_pr") != 1107 or candidate.get("stage") != "xi=11":
        raise KokunoA5RoutingError("latest self-contained candidate identity drift")
    if any(candidate.get(k) for k in (
        "superseded_by_agent1_1188_as_unified_candidate",
        "superseded_by_agent2_1189_as_unified_candidate",
        "superseded_by_agent3_1190_as_corrected_candidate",
    )):
        raise KokunoA5RoutingError("non-unified sibling falsely promoted")

    if artifact.get("core_states") != core_states():
        raise KokunoA5RoutingError("core-state promotion detected")
    protocol = artifact.get("fixed_validation_protocol", {})
    if protocol.get("st006_momentum_sampled_max") != ST006_MOMENTUM_SAMPLED_MAX or protocol.get("st006_momentum_volume_l2") != ST006_MOMENTUM_VOLUME_L2:
        raise KokunoA5RoutingError("ST006 baseline drift")
    if protocol.get("normalized_momentum_sampled_max_gate") != 1.0e-3 or protocol.get("normalized_momentum_volume_l2_gate") != 1.0e-3:
        raise KokunoA5RoutingError("momentum gate drift")
    if protocol.get("normalized_divergence_sampled_max_gate") != 1.0e-5 or protocol.get("normalized_divergence_volume_l2_gate") != 1.0e-5:
        raise KokunoA5RoutingError("divergence gate drift")
    if protocol.get("canonical_quadrature") != [24, 48, 96]:
        raise KokunoA5RoutingError("canonical quadrature drift")
    if not protocol.get("residual_defined_free_forcing_forbidden") or not protocol.get("posthoc_threshold_relaxation_forbidden"):
        raise KokunoA5RoutingError("scientific firewall weakened")
    if any(artifact.get("evidence_firewall", {}).values()):
        raise KokunoA5RoutingError("cross-identity evidence transfer detected")
    if not all(artifact.get("remaining_blockers", {}).values()):
        raise KokunoA5RoutingError("blocker cleared without same-identity evidence")
    if any(artifact.get("truth_boundary", {}).values()):
        raise KokunoA5RoutingError("scientific truth boundary promoted")


def dump_artifact(path: str | Path) -> Path:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(build_artifact(), sort_keys=True, indent=2) + "\n", encoding="utf-8")
    return target


__all__ = [
    "SCHEMA_VERSION", "TASK", "PARENT_A5_EXACT_HEAD",
    "AGENT1_EXACT_HEAD", "AGENT2_EXACT_HEAD", "AGENT3_EXACT_HEAD", "AGENT4_EXACT_HEAD",
    "KokunoA5RoutingError", "build_artifact", "validate_artifact", "dump_artifact",
]
