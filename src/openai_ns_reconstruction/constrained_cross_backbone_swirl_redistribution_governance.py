"""Fail-closed CR002 governance for parent-conditioned swirl redistribution.

This module changes no velocity values.  It prevents a parent-specific
first-order energy-balanced redistribution coordinate from being treated as one
portable physical parameter across ST048-S, ST050R-C, or future backbones.
"""
from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
from typing import Any, Mapping

TASK_ID = "CR002-CROSS-BACKBONE-ENERGY-NEUTRAL-REDISTRIBUTION-GOVERNANCE-068"
LIVE_ROUTE = "materialize_integrated_axial_cap_poloidal_child_then_fresh_validate"
CANONICAL_FAMILY = "eq45_supported_velocity_candidate_v1"
CANONICAL_API = "openai_ns_reconstruction.eq45_supported_delivery:velocity"
_ALLOWED_CLASSES = {
    "user_requirement",
    "public_source_fact",
    "autonomous_design",
    "pending_unknown",
}
_FALSE_STATES = (
    "visualization_ready",
    "visual_correspondence_verified",
    "pde_validated",
    "paper_exact",
    "openai_field_identified",
)


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _load_json(path: Path) -> dict[str, Any]:
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise ValueError(f"expected JSON object in {path}")
    return data


def load_scope() -> dict[str, Any]:
    return _load_json(_repo_root() / "configs" / "cross_backbone_swirl_redistribution_governance.json")


def load_constraints() -> dict[str, Any]:
    return _load_json(_repo_root() / "configs" / "constraints.json")


def load_project_status() -> dict[str, Any]:
    return _load_json(_repo_root() / "project_status.json")


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _require_equal(actual: Any, expected: Any, label: str) -> None:
    _require(actual == expected, f"{label} drifted: {actual!r} != {expected!r}")


def _require_close(actual: Any, expected: Any, label: str, tol: float = 1e-12) -> None:
    _require(abs(float(actual) - float(expected)) <= tol, f"{label} drifted: {actual!r} != {expected!r}")


def _audit_cr001(snapshot: Mapping[str, Any], constraints: Mapping[str, Any]) -> None:
    domain = constraints.get("domain", {})
    forcing = constraints.get("forcing", {})
    nontriviality = constraints.get("nontriviality", {})
    validation = constraints.get("validation", {})
    thresholds = validation.get("thresholds", {})

    _require_close(snapshot.get("nu"), constraints.get("nu"), "nu")
    _require_equal(snapshot.get("physical_domain"), domain.get("physical"), "physical domain")
    _require_equal(snapshot.get("evaluation_box"), domain.get("evaluation_box"), "evaluation box")
    _require_equal(snapshot.get("support"), domain.get("support"), "support")
    _require_equal(snapshot.get("time_interval"), domain.get("time_interval"), "time interval")
    _require_equal(snapshot.get("force_mode"), forcing.get("mode"), "force mode")
    _require_equal(snapshot.get("force_bounds"), forcing.get("parameters"), "force bounds")
    _require_close(snapshot.get("reference_energy"), nontriviality.get("reference_energy"), "reference energy")
    _require_close(
        snapshot.get("reference_energy_abs_tolerance"),
        nontriviality.get("reference_energy_abs_tolerance"),
        "reference energy tolerance",
    )
    _require_equal(snapshot.get("validation_seed"), validation.get("seed"), "validation seed")
    _require_equal(snapshot.get("held_out_points"), validation.get("held_out_points"), "held-out points")
    _require_equal(snapshot.get("validation_times"), validation.get("times"), "validation times")
    _require_equal(snapshot.get("derivative_steps"), validation.get("derivative_steps"), "derivative steps")
    for key in ("divergence_max", "divergence_L2", "pde_residual_max", "pde_residual_L2"):
        _require_close(snapshot.get(key), thresholds.get(key), key)

    restriction = str(forcing.get("restriction", ""))
    _require("No residual-dependent basis" in restriction, "residual-dependent forcing ban missing")
    _require("pointwise free force" in restriction, "pointwise free force ban missing")


def audit_cross_backbone_swirl_redistribution_governance(
    scope: Mapping[str, Any] | None = None,
    constraints: Mapping[str, Any] | None = None,
    project_status: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Audit parent-conditioned transfer semantics without promoting a child."""

    scope = deepcopy(dict(load_scope() if scope is None else scope))
    constraints = deepcopy(dict(load_constraints() if constraints is None else constraints))
    project_status = deepcopy(dict(load_project_status() if project_status is None else project_status))

    _require_equal(scope.get("schema_version"), 1, "schema version")
    _require_equal(scope.get("task_id"), TASK_ID, "task id")
    _require_equal(scope.get("status"), "governance_only_no_candidate_promotion", "status")

    _require_equal(set(scope.get("allowed_source_classes", [])), _ALLOWED_CLASSES, "source classes")
    source = scope.get("source_classification", {})
    _require(bool(source), "source classification missing")
    for label, classification in source.items():
        _require(classification in _ALLOWED_CLASSES, f"invalid source class for {label}: {classification!r}")
    _require_equal(source.get("callable_saveable_velocity_delivery"), "user_requirement", "delivery class")
    _require_equal(
        source.get("public_openai_velocity_field_visualization_as_observable_target"),
        "public_source_fact",
        "public target class",
    )
    for key in (
        "st048s_parent_and_temporal_piola_challenger",
        "st050rc_parent_challenger",
        "inner_outer_radial_windows",
        "parent_conditioned_energy_balance_alpha",
        "redistribution_gain_grid_and_clean_crossing",
        "common_reference_energy_normalization",
        "frozen_material_path_protocol",
    ):
        _require_equal(source.get(key), "autonomous_design", f"{key} class")
    for key in (
        "openai_hidden_radial_swirl_redistribution",
        "openai_hidden_coordinate_time_camera_alignment",
        "openai_hidden_numeric_velocity_field",
    ):
        _require_equal(source.get(key), "pending_unknown", f"{key} class")

    snapshots = scope.get("evidence_snapshots", {})
    st048 = snapshots.get("st048s", {})
    st050 = snapshots.get("st050rc", {})
    _require_equal(st048.get("redistribution_pr"), 447, "ST048 redistribution PR")
    _require_equal(st048.get("material_path_replay_pr"), 455, "ST048 material-path PR")
    _require_equal(st050.get("redistribution_pr"), 458, "ST050R-C redistribution PR")
    _require_equal(st050.get("material_path_replay_pr"), 466, "ST050R-C material-path PR")
    _require_equal(st048.get("inner_window"), st050.get("inner_window"), "shared inner window")
    _require_equal(st048.get("outer_window"), st050.get("outer_window"), "shared outer window")
    _require_equal(st048.get("gain_grid"), st050.get("gain_grid"), "shared diagnostic gain grid")
    _require_close(st048.get("first_clean_diagnostic_gain"), 0.025, "ST048 clean gain")
    _require_close(st050.get("first_clean_diagnostic_gain"), 0.025, "ST050 clean gain")
    _require_close(st048.get("balance_alpha"), 2.415185613025038, "ST048 alpha")
    _require_close(st050.get("balance_alpha"), 2.520520814687742, "ST050 alpha")
    _require(
        abs(float(st050.get("balance_alpha")) - float(st048.get("balance_alpha"))) > 0.1,
        "parent-conditioned alpha difference unexpectedly disappeared",
    )
    _require_equal(st048.get("material_path_evidence_status"), "successful_parent_local_diagnostic_only", "ST048 path status")
    _require_close(st048.get("material_path_mean_turns_relative_change"), 0.0150902, "ST048 mean-turn change", 1e-9)
    _require_close(st048.get("material_path_max_turns_relative_change"), 0.0221848, "ST048 max-turn change", 1e-9)

    _require_equal(st050.get("material_path_dedicated_run"), 35383631171, "ST050 dedicated run")
    _require_equal(
        st050.get("material_path_dedicated_run_status"),
        "execution_failed_before_numerical_receipt",
        "ST050 dedicated status",
    )
    _require_equal(st050.get("material_path_failure_class"), "ModuleNotFoundError", "ST050 failure class")
    _require_equal(st050.get("material_path_failure_missing_module"), "validate", "ST050 missing module")
    _require_equal(st050.get("material_path_standard_tests_status"), "success", "ST050 standard CI status")
    _require_equal(st050.get("material_path_research_publication_status"), "success", "ST050 publication CI status")
    _require_equal(st050.get("material_path_numerical_receipt_exists_at_this_head"), False, "ST050 numerical receipt")

    identity = scope.get("representation_identity", {})
    for key in (
        "alpha_is_part_of_transform_identity",
        "reference_energy_normalization_is_part_of_transform_identity",
        "parent_change_requires_alpha_recomputation_for_first_order_balance",
        "parent_change_requires_fresh_normalization",
        "materialized_child_requires_new_representation_identity",
        "materialized_child_requires_new_candidate_sha256",
    ):
        _require_equal(identity.get(key), True, key)
    for key in (
        "same_numeric_kappa_with_different_alpha_is_same_transform",
        "same_numeric_kappa_is_portable_production_parameter_across_parent_families",
        "same_windows_and_gain_grid_make_cross_parent_children_single_factor_ablations",
    ):
        _require_equal(identity.get(key), False, key)
    _require("alpha_parent" in str(identity.get("profile_definition", "")), "profile must expose parent-conditioned alpha")

    comparability = scope.get("cross_backbone_comparability", {})
    for key in (
        "within_parent_child_vs_parent_comparison_allowed_under_frozen_protocol",
        "parent_local_clean_crossing_is_capacity_evidence",
        "same_clean_crossing_on_two_parents_is_robustness_evidence",
        "valid_single_factor_backbone_ablation_requires_frozen_transform_profile_or_explicit_confound_report",
    ):
        _require_equal(comparability.get(key), True, key)
    for key in (
        "same_clean_crossing_on_two_parents_selects_production_gain",
        "same_clean_crossing_on_two_parents_recovers_openai_hidden_parameter",
        "cross_parent_child_difference_attributable_only_to_backbone_when_alpha_rebalanced",
        "cross_parent_child_difference_attributable_only_to_backbone_when_normalization_recomputed",
        "pathline_improvement_is_visual_correspondence_certificate",
        "pathline_improvement_is_pde_evidence",
        "eulerian_proxy_improvement_is_pde_evidence",
        "rank_or_novelty_is_pde_leverage_proof",
    ):
        _require_equal(comparability.get(key), False, key)

    ci = scope.get("ci_evidence_semantics", {})
    for key in (
        "green_standard_ci_may_replace_missing_dedicated_scientific_receipt",
        "focused_unit_tests_may_replace_failed_exact_source_replay",
        "failed_dedicated_replay_before_numerical_receipt_is_scientific_negative_result",
    ):
        _require_equal(ci.get(key), False, key)
    for key in (
        "failed_dedicated_replay_before_numerical_receipt_is_execution_or_portability_blocker",
        "numerical_material_path_claim_requires_successful_exact_source_replay_and_receipt",
    ):
        _require_equal(ci.get(key), True, key)

    promotion = scope.get("promotion_contract", {})
    for key in (
        "diagnostic_gain_may_become_production_bound_without_explicit_registration",
        "diagnostic_gain_may_become_production_value_without_explicit_selection",
        "parent_pressure_or_forcing_receipts_transfer_through_redistribution",
        "parent_pde_receipts_transfer_through_redistribution",
        "parent_material_paths_transfer_through_redistribution",
    ):
        _require_equal(promotion.get(key), False, key)
    for key in (
        "fresh_axis_support_core_divergence_checks_required",
        "fresh_material_paths_required_for_material_path_claim",
        "fresh_compatible_pressure_and_restricted_force_required_for_pde_claim",
        "fresh_unused_held_out_full_momentum_and_divergence_required_for_pde_promotion",
        "no_residual_defined_free_force",
        "no_zero_or_amplitude_collapse_shortcut",
    ):
        _require_equal(promotion.get(key), True, key)

    _audit_cr001(scope.get("canonical_cr001_snapshot", {}), constraints)

    routing = scope.get("routing_scope", {})
    _require_equal(routing.get("live_route_must_remain"), LIVE_ROUTE, "scope live route")
    _require_equal(routing.get("canonical_velocity_family_must_remain"), CANONICAL_FAMILY, "scope family")
    _require_equal(routing.get("canonical_velocity_api_must_remain"), CANONICAL_API, "scope API")
    _require_equal(project_status.get("active_scientific_route"), LIVE_ROUTE, "project live route")
    _require_equal(project_status.get("candidate_family"), CANONICAL_FAMILY, "project family")
    _require_equal(project_status.get("velocity_api"), CANONICAL_API, "project API")

    states = project_status.get("states", {})
    _require_equal(states.get("velocity_export_ready"), True, "velocity export state")
    for key in _FALSE_STATES:
        _require_equal(states.get(key), False, f"project state {key}")

    return {
        "task_id": TASK_ID,
        "audit": "pass",
        "same_numeric_gain_is_portable_across_parents": False,
        "st048_alpha": st048.get("balance_alpha"),
        "st050rc_alpha": st050.get("balance_alpha"),
        "st050rc_material_path_receipt_exists_at_audited_head": False,
        "canonical_velocity_unchanged": True,
        "cr001_unchanged": True,
    }


def main() -> int:
    print(json.dumps(audit_cross_backbone_swirl_redistribution_governance(), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
