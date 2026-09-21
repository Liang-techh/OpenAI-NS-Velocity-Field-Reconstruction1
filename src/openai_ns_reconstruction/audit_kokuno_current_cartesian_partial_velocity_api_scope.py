"""Fail-closed CR002 audit for the current Kokuno Cartesian-leading API scope.

Agent-1 #965 now exposes a real vectorized, save/load capable
``velocity(x,y,z,t)->[...,3]`` surface, but only on points whose governed
source-coordinate image satisfies ``X <= X_h``.  Agent-5 #968 correctly keeps
Kokuno ``velocity_export_ready`` false.  This audit machine-locks the missing
representation distinction: callable/serializable/identity-bound is not the
same thing as a total project-domain velocity provider.

The deterministic partial-provider witness in this module is autonomous
mechanics only.  It is not Kokuno/OpenAI numerical evidence.
"""
from __future__ import annotations

import copy
import json
from pathlib import Path
from typing import Any, Mapping

from . import kokuno_a5_current_cartesian_leading_ingest_contract as parent_a5


ROOT = Path(__file__).resolve().parents[2]
SCOPE_PATH = ROOT / "configs" / "kokuno_current_cartesian_partial_velocity_api_scope.json"
CONSTRAINTS_PATH = ROOT / "configs" / "constraints.json"
PROJECT_STATUS_PATH = ROOT / "project_status.json"

EXPECTED_PARENT_PR = 968
EXPECTED_PARENT_HEAD = "5558d0859d142914ae8921a4286f56077fe2232c"
EXPECTED_PARENT_SOURCE_BLOB = "a3456aa45f8372359a4c9b53c8ecab317291f2ca"
EXPECTED_A1_PR = 965
EXPECTED_A1_HEAD = "7dc09c59584a6b1bfaabcbb8df00283d9f9a40c1"
EXPECTED_A1_SOURCE_BLOB = "ca8b80b0451be1a8f31deaf620f542d8c2e92c0c"
EXPECTED_CONSTRAINTS_BLOB = "6c559e42895a606e2ef025ade4cb448966d75814"


class PartialVelocityAPIScopeError(RuntimeError):
    """Raised by the CLI verifier when the governed scope drifts."""


def _load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise PartialVelocityAPIScopeError(f"{path} must contain a JSON object")
    return value


def _toy_partial_velocity(s: float) -> tuple[float, float, float]:
    """Autonomous mechanics witness: callable on a strict subdomain only."""
    s = float(s)
    if s > 1.0:
        raise ValueError("toy provider is intentionally undefined for s>1")
    return (s, 0.0, 0.0)


def _mechanics_witness_ok() -> bool:
    if _toy_partial_velocity(0.5) != (0.5, 0.0, 0.0):
        return False
    try:
        _toy_partial_velocity(2.0)
    except ValueError:
        return True
    return False


def audit_scope(
    scope: Mapping[str, Any] | None = None,
    constraints: Mapping[str, Any] | None = None,
    project_status: Mapping[str, Any] | None = None,
) -> list[str]:
    """Return fail-closed drift labels; an empty list means governance consistency."""
    scope = copy.deepcopy(dict(scope)) if scope is not None else _load_json(SCOPE_PATH)
    constraints = (
        copy.deepcopy(dict(constraints)) if constraints is not None else _load_json(CONSTRAINTS_PATH)
    )
    project_status = (
        copy.deepcopy(dict(project_status))
        if project_status is not None
        else _load_json(PROJECT_STATUS_PATH)
    )
    errors: list[str] = []

    if scope.get("schema") != "cr002-kokuno-current-cartesian-partial-api-scope-v1":
        errors.append("scope_schema_drift")
    if scope.get("task") != "CR002-KOKUNO-CURRENT-CARTESIAN-PARTIAL-API-TOTALITY-097":
        errors.append("scope_task_drift")

    parent = scope.get("parent", {})
    if parent.get("pr") != EXPECTED_PARENT_PR:
        errors.append("parent_pr_drift")
    if parent.get("exact_head") != EXPECTED_PARENT_HEAD:
        errors.append("parent_head_drift")
    if parent.get("source_blob") != EXPECTED_PARENT_SOURCE_BLOB:
        errors.append("parent_source_blob_drift")

    upstream = scope.get("upstream_cartesian_leading", {})
    if upstream.get("pr") != EXPECTED_A1_PR:
        errors.append("a1_pr_drift")
    if upstream.get("exact_head") != EXPECTED_A1_HEAD:
        errors.append("a1_head_drift")
    if upstream.get("source_blob") != EXPECTED_A1_SOURCE_BLOB:
        errors.append("a1_source_blob_drift")
    if upstream.get("public_api") != ["velocity", "save_configuration", "load_configuration"]:
        errors.append("a1_public_api_drift")

    # Bind the exact Agent-5 registration semantics rather than inferring from prose.
    a1 = parent_a5.AGENT1_CARTESIAN_LEADING
    readiness = parent_a5.READINESS
    truth = parent_a5.TRUTH_BOUNDARY
    if a1.get("pr") != EXPECTED_A1_PR or a1.get("head") != EXPECTED_A1_HEAD:
        errors.append("parent_a5_a1_identity_drift")
    if a1.get("source_blob") != EXPECTED_A1_SOURCE_BLOB:
        errors.append("parent_a5_a1_blob_drift")
    if a1.get("public_api") != ["velocity", "save_configuration", "load_configuration"]:
        errors.append("parent_a5_public_api_drift")
    if a1.get("vectorized_cartesian_velocity") is not True:
        errors.append("parent_a5_vectorized_callable_lost")
    if a1.get("candidate_save_reload") is not True:
        errors.append("parent_a5_save_reload_lost")
    if a1.get("velocity_materialized_through_xh") is not True:
        errors.append("parent_a5_through_xh_materialization_lost")
    if a1.get("velocity_beyond_xh_materialized") is not False:
        errors.append("parent_a5_beyond_xh_promotion")
    if a1.get("global_compact_supported_leading_velocity") is not False:
        errors.append("parent_a5_global_compact_promotion")
    if readiness.get("velocity_export_ready") is not False:
        errors.append("parent_a5_velocity_export_promotion")
    if truth.get("complete_velocity_candidate_materialized") is not False:
        errors.append("parent_a5_complete_velocity_promotion")
    if truth.get("complete_candidate_api_ready") is not False:
        errors.append("parent_a5_complete_api_promotion")
    if truth.get("global_cartesian_leading_velocity_materialized") is not False:
        errors.append("parent_a5_global_velocity_promotion")
    if truth.get("velocity_beyond_xh_materialized") is not False:
        errors.append("parent_a5_truth_beyond_xh_promotion")

    delivery = scope.get("delivery_contract", {})
    required_true = (
        "cartesian_velocity_callable_through_current_xh",
        "cartesian_velocity_vectorized_through_current_xh",
        "configuration_save_load_materialized",
        "semantic_identity_materialized",
        "calls_beyond_current_xh_fail_closed",
        "provider_domain_is_partial_relative_to_global_project_delivery",
        "canonical_eq45_velocity_export_ready_independent",
    )
    required_false = (
        "global_or_r3_total_velocity_materialized",
        "velocity_beyond_current_xh_materialized",
        "complete_candidate_api_ready",
        "kokuno_velocity_export_ready",
        "callable_implies_global_api_totality",
        "save_load_implies_global_api_totality",
        "semantic_identity_implies_global_api_totality",
        "pde_pending_blocks_callable_velocity_delivery",
    )
    for key in required_true:
        if delivery.get(key) is not True:
            errors.append(f"delivery_required_true_{key}")
    for key in required_false:
        if delivery.get(key) is not False:
            errors.append(f"delivery_forbidden_promotion_{key}")

    provenance = scope.get("provenance", {})
    if set(provenance) != {
        "user_requirement",
        "public_source_fact",
        "autonomous_design",
        "pending_unknown",
    }:
        errors.append("four_way_provenance_partition_drift")
    witness = scope.get("mechanics_witness", {})
    if witness.get("classification") != "autonomous_design":
        errors.append("mechanics_witness_source_laundering")
    if witness.get("not_public_source_or_candidate_evidence") is not True:
        errors.append("mechanics_witness_evidence_laundering")
    if not _mechanics_witness_ok():
        errors.append("mechanics_partial_provider_witness_failed")

    truth_scope = scope.get("truth_boundary", {})
    if truth_scope.get("current_cartesian_leading_velocity_materialized_through_xh") is not True:
        errors.append("scope_through_xh_truth_lost")
    if truth_scope.get("current_cartesian_leading_save_load_materialized") is not True:
        errors.append("scope_save_load_truth_lost")
    forbidden_truth_promotions = (
        "global_cartesian_leading_velocity_materialized",
        "velocity_beyond_xh_materialized",
        "global_compact_supported_leading_velocity_materialized",
        "leading_plus_oscillatory_cartesian_velocity_materialized",
        "complete_velocity_candidate_materialized",
        "matched_cartesian_pressure_materialized",
        "restricted_forcing_materialized",
        "canonical_whole_domain_divergence_l2_assessed",
        "same_protocol_full_ns_residual_available",
        "visual_correspondence_verified",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    )
    for key in forbidden_truth_promotions:
        if truth_scope.get(key) is not False:
            errors.append(f"truth_boundary_promotion_{key}")

    # Canonical CR001 remains immutable and is not replaced by this scoped API audit.
    lock = scope.get("cr001_lock", {})
    if lock.get("constraints_blob") != EXPECTED_CONSTRAINTS_BLOB:
        errors.append("constraints_blob_drift")
    domain = constraints.get("domain", {})
    forcing = constraints.get("forcing", {})
    nontriviality = constraints.get("nontriviality", {})
    validation = constraints.get("validation", {})
    thresholds = validation.get("thresholds", {})
    if constraints.get("nu") != 0.01 or lock.get("nu") != 0.01:
        errors.append("cr001_nu_drift")
    if domain.get("physical") != "R^3" or lock.get("physical_domain") != "R^3":
        errors.append("cr001_physical_domain_drift")
    if domain.get("evaluation_box") != [[-2, 2], [-2, 2], [-2, 2]]:
        errors.append("cr001_evaluation_box_drift")
    if lock.get("evaluation_box") != [[-2, 2], [-2, 2], [-2, 2]]:
        errors.append("scope_evaluation_box_drift")
    if domain.get("support") != "r < 2 and abs(z) < 2" or lock.get("support") != domain.get("support"):
        errors.append("cr001_support_drift")
    if domain.get("time_interval") != [0.25, 0.75] or lock.get("time_interval") != [0.25, 0.75]:
        errors.append("cr001_time_interval_drift")
    if forcing.get("mode") != "restricted_two_parameter_family":
        errors.append("cr001_forcing_mode_drift")
    restriction = str(forcing.get("restriction", ""))
    if "No residual-dependent basis or pointwise free force" not in restriction:
        errors.append("cr001_free_force_firewall_lost")
    if nontriviality.get("reference_energy") != 1.0:
        errors.append("cr001_reference_energy_drift")
    if nontriviality.get("reference_energy_abs_tolerance") != 0.001:
        errors.append("cr001_energy_tolerance_drift")
    if "reject collapsed candidates" not in str(nontriviality.get("enforcement", "")):
        errors.append("cr001_collapse_firewall_lost")
    if validation.get("seed") != 914027 or validation.get("held_out_points") != 4096:
        errors.append("cr001_validation_split_drift")
    if validation.get("derivative_steps") != [0.02, 0.01, 0.005]:
        errors.append("cr001_derivative_ladder_drift")
    if validation.get("quadrature_orders_per_axis") != [24, 48, 96]:
        errors.append("cr001_quadrature_ladder_drift")
    if thresholds.get("pde_residual_max") != 0.001 or thresholds.get("pde_residual_L2") != 0.001:
        errors.append("cr001_momentum_gate_drift")
    if thresholds.get("divergence_max") != 1.0e-5 or thresholds.get("divergence_L2") != 1.0e-5:
        errors.append("cr001_divergence_gate_drift")
    if lock.get("residual_defined_free_forcing_forbidden") is not True:
        errors.append("scope_free_force_firewall_lost")
    if lock.get("candidate_collapse_forbidden") is not True:
        errors.append("scope_candidate_collapse_firewall_lost")
    if lock.get("post_hoc_threshold_relaxation_forbidden") is not True:
        errors.append("scope_threshold_relaxation_firewall_lost")

    # Keep canonical callable delivery independent from this incomplete Kokuno sibling.
    states = project_status.get("states", {})
    if project_status.get("candidate_family") != "eq45_supported_velocity_candidate_v1":
        errors.append("canonical_delivery_family_drift")
    if states.get("velocity_export_ready") is not True:
        errors.append("canonical_eq45_velocity_export_ready_lost")
    for key in ("visual_correspondence_verified", "pde_validated", "paper_exact", "openai_field_identified"):
        if states.get(key) is not False:
            errors.append(f"canonical_delivery_truth_promotion_{key}")

    return errors


def verify_scope() -> None:
    errors = audit_scope()
    if errors:
        raise PartialVelocityAPIScopeError("; ".join(errors))


def _main() -> None:
    verify_scope()
    print("CR002 partial Cartesian velocity API scope: PASS")


if __name__ == "__main__":
    _main()
