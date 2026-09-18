"""Fail-closed CR002 governance for the integrated axial-cap poloidal route.

This module changes no velocity value. It audits provenance, nonlinear morphology
semantics, the canonical CR001 contract, and the live materialization-routing
boundary around the integrated axial-cap capacity direction.
"""
from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
from typing import Any, Mapping

TASK_ID = "CR002-BIPOLAR-AXIAL-CAP-POLOIDAL-REPRESENTATION-SCOPE-038"
MODE_NAME = "AXIAL_CAP_BANDED_C4_ODD_Z_POLOIDAL"
_ALLOWED_CLASSES = {
    "user_requirement",
    "public_source_fact",
    "autonomous_design",
    "pending_unknown",
}
_FALSE_TRUTH_STATES = (
    "axial_cap_poloidal_nonzero_child_materialized",
    "axial_cap_poloidal_coefficient_selected",
    "candidate_selection_resolved",
    "visualization_ready",
    "visual_correspondence_verified",
    "pde_validated",
    "paper_exact",
    "openai_field_identified",
    "blowup_proved",
)
_PROMOTION_REQUIREMENTS = [
    "explicit_autonomous_coefficient_bound",
    "explicit_nonzero_coefficient_value",
    "new_representation_family_identity",
    "new_candidate_sha256",
    "fresh_full_candidate_validation",
]


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _load_json(path: Path) -> dict[str, Any]:
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise ValueError(f"expected JSON object in {path}")
    return data


def load_scope() -> dict[str, Any]:
    return _load_json(_repo_root() / "configs" / "bipolar_axial_cap_poloidal_representation_scope.json")


def load_constraints() -> dict[str, Any]:
    return _load_json(_repo_root() / "configs" / "constraints.json")


def load_project_status() -> dict[str, Any]:
    return _load_json(_repo_root() / "project_status.json")


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _require_equal(actual: Any, expected: Any, label: str) -> None:
    _require(actual == expected, f"{label} drifted: {actual!r} != {expected!r}")


def audit_axial_cap_poloidal_scope(
    scope: Mapping[str, Any] | None = None,
    constraints: Mapping[str, Any] | None = None,
    project_status: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Audit integrated axial-cap evidence without promoting it into a candidate."""

    scope = deepcopy(dict(load_scope() if scope is None else scope))
    constraints = deepcopy(dict(load_constraints() if constraints is None else constraints))
    project_status = deepcopy(dict(load_project_status() if project_status is None else project_status))

    _require_equal(scope.get("task_id"), TASK_ID, "scope task id")
    upstream = scope.get("governed_upstream", {})
    _require_equal(upstream.get("pr"), 287, "upstream PR")
    _require_equal(upstream.get("task_id"), "CR003-BIPOLAR-AXIAL-CAP-POLOIDAL-CAPACITY-043", "upstream task")
    _require_equal(upstream.get("head_sha"), "080936d6a2ad82ce5e9a69b27e8eba97550ad911", "upstream head")
    _require_equal(upstream.get("mode_name"), MODE_NAME, "mode name")
    _require_equal(
        scope.get("representation_kind"),
        "autonomous_physical_space_additive_curl_correction_capacity_only",
        "representation kind",
    )

    allowed = set(scope.get("allowed_source_classes", ()))
    _require_equal(allowed, _ALLOWED_CLASSES, "source-class vocabulary")
    source = scope.get("source_classification", {})
    _require(source and set(source.values()) <= _ALLOWED_CLASSES, "invalid source classification")
    _require_equal(source.get("callable_velocity_delivery"), "user_requirement", "velocity delivery class")
    _require_equal(source.get("eq45_backbone"), "public_source_fact", "Eq45 backbone class")
    for key in (
        "axial_cap_band_and_envelope",
        "cap_band_0p55_to_0p95",
        "diagnostic_trial_amplitude_0p25",
    ):
        _require_equal(source.get(key), "autonomous_design", f"source class {key}")
    _require_equal(source.get("openai_hidden_numerical_profile"), "pending_unknown", "hidden profile class")

    capacity = scope.get("capacity_evidence_only", {})
    _require_equal(capacity.get("public_velocity_rank"), 5, "public velocity rank")
    _require(capacity.get("normalized_condition_number", 0.0) > 1.0, "condition number missing")
    _require(0.9 < capacity.get("novelty_outside_previous_four_direction_span", 0.0) <= 1.0, "novelty drifted")
    _require(0.0 < capacity.get("cartesian_fd_divergence_max_abs", 1.0) < 1e-7, "basis divergence evidence drifted")
    _require_equal(capacity.get("diagnostic_coefficients"), [-0.25, 0.25], "diagnostic coefficients")
    for key in (
        "diagnostic_coefficients_select_value",
        "diagnostic_coefficients_define_materialization_bound",
        "inherited_profile_limit_is_materialization_bound",
        "basis_divergence_is_registered_full_candidate_divergence",
        "public_velocity_rank_is_momentum_residual_jacobian_rank",
        "morphology_capacity_is_visual_correspondence",
        "morphology_capacity_is_pde_improvement",
    ):
        _require_equal(capacity.get(key), False, key)

    nonlinear = scope.get("nonlinear_morphology_semantics", {})
    _require(nonlinear.get("diagnostic_abs_coefficient_q90_over_Zp", 0.0) > nonlinear.get("baseline_q90_over_Zp", 1.0), "q90 capacity missing")
    _require(nonlinear.get("diagnostic_abs_coefficient_q99_over_Zp", 0.0) > nonlinear.get("baseline_q99_over_Zp", 1.0), "q99 capacity missing")
    _require_equal(nonlinear.get("plus_minus_q90_span"), 0.0, "q90 plus/minus span")
    _require_equal(nonlinear.get("plus_minus_q99_span"), 0.0, "q99 plus/minus span")
    _require_equal(nonlinear.get("sign_symmetric_reach_change_observed"), True, "sign-symmetric reach observation")
    _require_equal(nonlinear.get("reach_change_is_even_quadratic_enstrophy_effect"), True, "quadratic enstrophy interpretation")
    for key in (
        "zero_plus_minus_span_means_zero_capacity",
        "reach_change_selects_coefficient_sign",
        "reach_change_establishes_linear_sensitivity",
        "reach_change_establishes_openai_tip_geometry",
        "reach_change_establishes_visual_correspondence",
    ):
        _require_equal(nonlinear.get(key), False, key)

    routing = scope.get("integration_routing_scope", {})
    _require_equal(
        routing.get("scope_kind"),
        "active_child_materialization_route_with_prematerialization_evidence",
        "routing scope kind",
    )
    _require_equal(routing.get("live_integrated_mode"), MODE_NAME, "live integrated mode")
    _require_equal(
        routing.get("live_route"),
        project_status.get("active_scientific_route"),
        "live scientific route",
    )
    _require_equal(routing.get("live_mode_nonzero_child_materialized"), False, "live materialization state")
    _require_equal(
        routing.get("axial_cap_capacity_is_integrated_prematerialization_evidence"),
        True,
        "axial-cap integration state",
    )
    _require_equal(
        routing.get("additional_capacity_growth_allowed_before_live_materialization_screen"),
        False,
        "capacity-growth guard",
    )
    _require_equal(
        routing.get("route_guard"),
        "axial_cap_prematerialization_evidence_does_not_itself_constitute_production_materialization_promotion_or_pde_validation",
        "route semantic guard",
    )
    _require_equal(routing.get("promotion_requires"), _PROMOTION_REQUIREMENTS, "promotion requirements")

    materialization = scope.get("future_nonzero_materialization", {})
    for key in (
        "explicit_autonomous_coefficient_bound_required",
        "explicit_coefficient_value_required",
        "new_representation_family_identity_required",
        "new_candidate_sha_required",
        "recheck_reference_energy",
        "recheck_validation_time_energy",
        "recheck_core_rotation_sign_nontriviality",
        "recheck_bipolar_radial_axial_direction_parity",
        "recheck_axis_support_regularity",
        "fresh_full_per_component_momentum_divergence_required",
        "model_selection_data_may_not_be_reused_as_acceptance_data",
    ):
        _require_equal(materialization.get(key), True, key)

    snap = scope.get("canonical_cr001_snapshot", {})
    domain = constraints.get("domain", {})
    forcing = constraints.get("forcing", {})
    nontriviality = constraints.get("nontriviality", {})
    validation = constraints.get("validation", {})
    thresholds = validation.get("thresholds", {})
    expected_pairs = (
        (snap.get("nu"), constraints.get("nu"), "nu"),
        (snap.get("physical_domain"), domain.get("physical"), "physical domain"),
        (snap.get("evaluation_box"), domain.get("evaluation_box"), "evaluation box"),
        (snap.get("support"), domain.get("support"), "support"),
        (snap.get("time_interval"), domain.get("time_interval"), "time interval"),
        (snap.get("force_mode"), forcing.get("mode"), "force mode"),
        (snap.get("force_bounds"), forcing.get("parameters"), "force bounds"),
        (snap.get("reference_energy"), nontriviality.get("reference_energy"), "reference energy"),
        (snap.get("reference_energy_abs_tolerance"), nontriviality.get("reference_energy_abs_tolerance"), "reference energy tolerance"),
        (
            snap.get("validation_energy_range"),
            [nontriviality.get("minimum_energy_each_validation_time"), nontriviality.get("maximum_energy_each_validation_time")],
            "validation energy range",
        ),
        (snap.get("validation_seed"), validation.get("seed"), "validation seed"),
        (snap.get("held_out_points"), validation.get("held_out_points"), "held-out points"),
        (snap.get("derivative_steps"), validation.get("derivative_steps"), "derivative steps"),
        (snap.get("divergence_max"), thresholds.get("divergence_max"), "divergence max"),
        (snap.get("divergence_L2"), thresholds.get("divergence_L2"), "divergence L2"),
        (snap.get("pde_residual_max"), thresholds.get("pde_residual_max"), "PDE max"),
        (snap.get("pde_residual_L2"), thresholds.get("pde_residual_L2"), "PDE L2"),
    )
    for expected, actual, label in expected_pairs:
        _require_equal(actual, expected, label)
    _require("No residual-dependent basis" in forcing.get("restriction", ""), "free-force restriction weakened")

    latest = project_status.get("latest_integrated_axial_cap_capacity", {})
    _require_equal(latest.get("mode"), MODE_NAME, "project integrated axial-cap mode")
    _require("capacity_integrated_nonzero_child_not_materialized" in latest.get("status", ""), "axial-cap unexpectedly materialized")
    _require_equal(latest.get("production_coefficient_bound_selected"), False, "integrated bound selection")
    _require_equal(latest.get("production_coefficient_value_selected"), False, "integrated value selection")
    _require_equal(latest.get("candidate_sha_created"), False, "integrated candidate SHA state")
    _require_equal(latest.get("fresh_full_momentum_validated"), False, "integrated full-momentum state")

    next_task = project_status.get("next_integration_task", "")
    for token in (
        MODE_NAME,
        "explicit autonomous coefficient bound/value",
        "new representation identity/candidate SHA",
        "fresh full-momentum/divergence validation",
        "do not add another capacity basis first",
    ):
        _require(token in next_task, f"live next integration task lost token: {token}")
    avoided = set(project_status.get("integration_policy", {}).get("avoid", ()))
    _require(
        "additional_capacity_basis_growth_before_a_current_integrated_direction_is_materialized_and_screened" in avoided,
        "live capacity-growth block missing",
    )

    truth = scope.get("truth_boundary", {})
    _require_equal(truth.get("velocity_export_ready_may_remain_true_for_existing_candidate"), True, "export independence")
    for key in _FALSE_TRUTH_STATES:
        _require_equal(truth.get(key), False, key)

    states = project_status.get("states", {})
    _require_equal(states.get("velocity_export_ready"), True, "current export state")
    for key in (
        "visualization_ready",
        "pde_validated",
        "physical_support_validated",
        "visual_correspondence_verified",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    ):
        _require_equal(states.get(key), False, f"project state {key}")

    return {
        "task_id": TASK_ID,
        "status": "governance_pass",
        "mode_name": MODE_NAME,
        "upstream_pr": 287,
        "source_classification": source,
        "canonical_thresholds_unchanged": True,
        "diagnostic_amplitude_selected": False,
        "nonlinear_reach_is_capacity_only": True,
        "axial_cap_capacity_status": "integrated_prematerialization",
        "live_next_integration_task_preserved": True,
        "velocity_export_ready": True,
        "visual_correspondence_verified": False,
        "pde_validated": False,
    }


def main() -> None:
    print(json.dumps(audit_axial_cap_poloidal_scope(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
