"""Fail-closed Kokuno Agent-5 routing artifact for the l=-1/AuxT2/RF49 seam.

This module is integration/provenance glue only. It records exact sibling
identities and the shortest remaining closure without reimplementing Agent-1--4
mathematics or converting scoped/queued evidence into PDE validation.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

SCHEMA_VERSION = "kokuno-a5-lminus1-auxt2-rf49-audit-routing-v1"
TASK = "KOKUNO-A5-LMINUS1-AUXT2-RF49-AUDIT-ROUTING-126"

PARENT_A5_PR = 1192
PARENT_A5_EXACT_HEAD = "18f9df62154dd5230eede98eb3634a74560b9e97"

AGENT1_PR = 1196
AGENT1_EXACT_HEAD = "e0cf691884fd371c477de5ae55a18cdabbd72e81"
AGENT1_SOURCE_BLOB = "466c87d26524a6f8b2931da9d501a2795d55db08"
AGENT1_PATH = "src/openai_ns_reconstruction/kokuno_pa16_current_cartesian_postswirl_lminus1_hold.py"

AGENT2_PR = 1198
AGENT2_EXACT_HEAD = "5765c2b3bab7482df514efa87f9b6bba48e04b7f"
AGENT2_SOURCE_BLOB = "96168ac6583ad5e71aa044bd558c6feeadf051f4"
AGENT2_PATH = "src/openai_ns_reconstruction/kokuno_current_i4_auxiliary_t2_provider.py"

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
    """Raised when a routing identity or truth boundary is promoted incorrectly."""


def upstream_routes() -> dict[str, dict[str, Any]]:
    return {
        "agent1": {
            "pr": AGENT1_PR,
            "exact_head": AGENT1_EXACT_HEAD,
            "source_blob": AGENT1_SOURCE_BLOB,
            "path": AGENT1_PATH,
            "delivery": "current Cartesian leading through the public post-swirl l=-1 hold",
            "post_relative_swirl_release1_materialized": True,
            "source_l_minus1_hold_materialized": True,
            "current_autonomous_h": 0.005,
            "hold_length_4log1overh": 21.193269466192145,
            "decimal_output_significant_digits": 96,
            "end_to_end_96digit_arithmetic_materialized": False,
            "source_l_minus1_to_minus_h_transition_materialized": False,
            "source_terminal_multiplier_materialized": False,
            "source_exterior_heat_replacement_materialized": False,
            "outer_global_leading_velocity_materialized": False,
            "matched_global_pressure_materialized": False,
            "cross_language_global_velocity_export_ready": False,
            "heldout_complete_ns_residual_assessed": False,
        },
        "agent2": {
            "pr": AGENT2_PR,
            "exact_head": AGENT2_EXACT_HEAD,
            "source_blob": AGENT2_SOURCE_BLOB,
            "path": AGENT2_PATH,
            "delivery": "same-identity current-I4 raw auxiliary-T2 wave provider",
            "current_i4_raw_auxiliary_t2_provider_materialized": True,
            "recomputed_from_actual_frozen_candidate": True,
            "surrogate_or_preaveraged_covariance_used": False,
            "heldout_data_used": False,
            "residual_as_forcing_used": False,
            "repository_autonomous_phase_lift": True,
            "source_exact_auxiliary_phase_assignment_recovered": False,
            "provider_blob_available_for_agent3_pin": True,
            "provider_checksum_pinned_by_agent3": False,
            "rf30_repository_candidate_state_authorized": False,
            "complete_ns_residual_assessed": False,
        },
        "agent3": {
            "pr": AGENT3_PR,
            "exact_head": AGENT3_EXACT_HEAD,
            "source_blob": AGENT3_SOURCE_BLOB,
            "path": AGENT3_PATH,
            "delivery": "identity-preserving RF34-RF39 to RF44-RF49 recompute handoff mechanics",
            "typed_rf30_to_rf39_chain_available": True,
            "rf44_rf49_postupdate_handoff_executable": True,
            "pinned_same_identity_repository_provider_blob": None,
            "sees_agent2_1198_provider_as_pinned": False,
            "repository_candidate_rf30_evidence": False,
            "repository_candidate_rf44_rf49_evidence": False,
            "cartesian_correction_velocity_materialized": False,
            "finite_correction_cycle_run": False,
            "heldout_complete_ns_residual_assessed": False,
        },
        "agent4": {
            "pr": AGENT4_PR,
            "exact_head": AGENT4_EXACT_HEAD,
            "source_blob": AGENT4_SOURCE_BLOB,
            "path": AGENT4_PATH,
            "delivery": "implementation-distinct save/load FD4 divergence audit of exact A1 #1196 l=-1 hold",
            "audits_agent1_pr": AGENT1_PR,
            "audits_agent1_exact_head": AGENT1_EXACT_HEAD,
            "save_load_public_velocity_only_scientific_path": True,
            "scoped_relative_fd_ladder": [2.0e-6, 1.0e-6, 5.0e-7],
            "canonical_absolute_fd_ladder": [0.02, 0.01, 0.005],
            "scoped_only": True,
            "scoped_gate_passed": None,
            "canonical_absolute_fd_transfer_ready": False,
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
        "superseded_by_agent1_1196_as_unified_candidate": False,
        "superseded_by_agent2_1198_as_unified_candidate": False,
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
        "scoped_divergence_is_complete_ns_validation": False,
    }


def evidence_firewall() -> dict[str, bool]:
    return {
        "agent2_1198_provider_existence_alone_authorizes_agent3_rf30": False,
        "agent3_1197_unpinned_handoff_creates_scientific_correction": False,
        "agent3_1197_mechanics_authorizes_cartesian_delta_u": False,
        "agent4_1199_scoped_audit_validates_agent2_1198": False,
        "agent4_1199_scoped_audit_validates_agent3_1197": False,
        "agent4_1199_scoped_audit_validates_complete_ns": False,
        "agent4_1199_scoped_audit_may_set_pde_validated": False,
        "queued_or_running_ci_means_pass": False,
    }


def remaining_blockers() -> dict[str, bool]:
    """True means the blocker remains unresolved."""
    return {
        "a1_l_minus1_hold": False,
        "a1_minus1_to_minus_h_transition": True,
        "a1_terminal_multiplier_exterior_heat_and_global_leading": True,
        "matched_cartesian_pressure_and_grad_p": True,
        "cross_language_global_velocity_export_representation": True,
        "matching_global_self_contained_a2_composite": True,
        "raw_same_identity_auxiliary_t2_provider_missing": False,
        "agent3_checksum_pin_of_exact_agent2_1198_provider_blob": True,
        "scientific_candidate_specific_rf30_to_rf49_correction": True,
        "cartesian_correction_velocity": True,
        "rf44_rf49_repository_candidate_remainder": True,
        "preregistered_restricted_non_residual_defined_forcing": True,
        "complete_identity_bound_ns_defect": True,
        "real_finite_correction_cycle": True,
        "matching_agent4_global_corrected_audit": True,
        "agent4_heldout_canonical_complete_ns_gate": True,
    }


def shortest_closure() -> list[str]:
    return [
        "A3: independently consume and checksum-pin exact A2 #1198 provider blob 96168ac6583ad5e71aa044bd558c6feeadf051f4, then replay the existing current-I4 RF30->RF49 chain on that exact identity; do not promote to Cartesian correction until the repository-candidate receipt is authorized.",
        "A1: continue exact #1196 through the public -1->-h transition, terminal multiplier and exterior/global leading; materialize matched pressure and a deterministic Python/MATLAB-safe export representation.",
        "A3: after provider authorization, materialize the correction as Cartesian delta-u, recompute the nonlinear remainder, and run an actual finite correction cycle without held-out tuning.",
        "A2: once global A1 exists, build a matching self-contained global leading+oscillatory candidate instead of transferring current-I4 provider evidence across identities.",
        "Integration/A4: bind matched pressure and preregistered restricted non-residual-defined forcing to one corrected semantic identity, form the complete NS defect, then enforce the unchanged held-out/canonical [24,48,96] momentum <=1e-3 and divergence <=1e-5 gates before any pde_validated promotion.",
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
    if not a1.get("source_l_minus1_hold_materialized"):
        raise KokunoA5RoutingError("A1 l=-1 hold delivery lost")
    if a1.get("source_l_minus1_to_minus_h_transition_materialized") or a1.get("outer_global_leading_velocity_materialized"):
        raise KokunoA5RoutingError("A1 later/global stages falsely promoted")
    if a1.get("end_to_end_96digit_arithmetic_materialized"):
        raise KokunoA5RoutingError("A1 output encoding mislabelled as construction arithmetic")

    if not a2.get("current_i4_raw_auxiliary_t2_provider_materialized") or not a2.get("provider_blob_available_for_agent3_pin"):
        raise KokunoA5RoutingError("A2 raw provider delivery lost")
    if a2.get("provider_checksum_pinned_by_agent3") or a2.get("rf30_repository_candidate_state_authorized"):
        raise KokunoA5RoutingError("A2 provider availability falsely promoted to A3 authorization")
    if a2.get("surrogate_or_preaveraged_covariance_used") or a2.get("residual_as_forcing_used"):
        raise KokunoA5RoutingError("A2 provider provenance firewall violated")

    if not a3.get("rf44_rf49_postupdate_handoff_executable"):
        raise KokunoA5RoutingError("A3 RF44-RF49 handoff delivery lost")
    if a3.get("pinned_same_identity_repository_provider_blob") is not None or a3.get("sees_agent2_1198_provider_as_pinned"):
        raise KokunoA5RoutingError("A3 provider pin falsely inferred")
    if a3.get("repository_candidate_rf30_evidence") or a3.get("repository_candidate_rf44_rf49_evidence"):
        raise KokunoA5RoutingError("A3 mechanics falsely promoted to repository-candidate evidence")
    if a3.get("cartesian_correction_velocity_materialized") or a3.get("finite_correction_cycle_run"):
        raise KokunoA5RoutingError("A3 correction/cycle falsely promoted")

    if a4.get("audits_agent1_exact_head") != AGENT1_EXACT_HEAD:
        raise KokunoA5RoutingError("A4 audit is not bound to exact A1 #1196")
    if a4.get("scoped_gate_passed") is not None or a4.get("scientifically_admitted") or a4.get("pde_validated"):
        raise KokunoA5RoutingError("A4 unresolved scoped audit over-promoted")

    candidate = artifact.get("latest_self_contained_candidate", {})
    if candidate.get("agent2_pr") != LATEST_SELF_CONTAINED_A2_PR or candidate.get("agent1_pr") != LATEST_SELF_CONTAINED_A1_PR or candidate.get("stage") != LATEST_SELF_CONTAINED_STAGE:
        raise KokunoA5RoutingError("latest self-contained candidate identity drift")
    for key in (
        "superseded_by_agent1_1196_as_unified_candidate",
        "superseded_by_agent2_1198_as_unified_candidate",
        "superseded_by_agent3_1197_as_corrected_candidate",
    ):
        if candidate.get(key):
            raise KokunoA5RoutingError("cross-identity candidate promotion")

    if artifact.get("core_states") != core_states():
        raise KokunoA5RoutingError("core state promotion")
    if any(artifact.get("evidence_firewall", {}).values()):
        raise KokunoA5RoutingError("evidence firewall opened")

    blockers = artifact.get("remaining_blockers", {})
    if blockers.get("a1_l_minus1_hold") is not False or blockers.get("raw_same_identity_auxiliary_t2_provider_missing") is not False:
        raise KokunoA5RoutingError("newly closed blocker not registered")
    for key in (
        "a1_minus1_to_minus_h_transition",
        "agent3_checksum_pin_of_exact_agent2_1198_provider_blob",
        "cartesian_correction_velocity",
        "complete_identity_bound_ns_defect",
        "agent4_heldout_canonical_complete_ns_gate",
    ):
        if blockers.get(key) is not True:
            raise KokunoA5RoutingError(f"remaining blocker falsely closed: {key}")

    protocol = artifact.get("fixed_validation_protocol", {})
    if protocol.get("normalized_momentum_sampled_max_gate") != FINAL_NORMALIZED_MOMENTUM_GATE:
        raise KokunoA5RoutingError("momentum gate drift")
    if protocol.get("normalized_divergence_sampled_max_gate") != FINAL_NORMALIZED_DIVERGENCE_GATE:
        raise KokunoA5RoutingError("divergence gate drift")
    if protocol.get("canonical_quadrature") != CANONICAL_QUADRATURE:
        raise KokunoA5RoutingError("quadrature drift")
    if not protocol.get("residual_defined_free_forcing_forbidden") or not protocol.get("posthoc_threshold_relaxation_forbidden"):
        raise KokunoA5RoutingError("validation firewall weakened")

    truth = artifact.get("truth_boundary", {})
    if any(truth.values()):
        raise KokunoA5RoutingError("scientific truth falsely promoted")


def write_artifact(path: str | Path) -> Path:
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(build_artifact(), sort_keys=True, indent=2) + "\n", encoding="utf-8")
    return out


def main() -> None:
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", default="kokuno_a5_lminus1_auxt2_rf49_audit_routing.json")
    args = parser.parse_args()
    write_artifact(args.output)


if __name__ == "__main__":
    main()
