"""Fail-closed CR002 governance for the open ST048-S post-Piola swirl-gain screen.

This module changes no velocity values. It keeps PR #428's constant azimuthal
capacity coordinate separate from production selection, immutable candidate
identity, the live route, visual correspondence, and Navier-Stokes validation.
"""
from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
from typing import Any, Mapping

TASK_ID = "CR002-ST048S-POST-PIOLA-SWIRL-GAIN-GOVERNANCE-065"
UPSTREAM_PR = 428
UPSTREAM_HEAD = "b2b4888ece83a9f862435ce1adb847c94ad3ad7d"
PARENT_PR = 418
PARENT_HEAD = "ae10c3ccacd16e1f5e0054502f1dc5a6f68320b5"
LIVE_ROUTE = "materialize_integrated_axial_cap_poloidal_child_then_fresh_validate"
LIVE_MODE = "AXIAL_CAP_BANDED_C4_ODD_Z_POLOIDAL"
_ALLOWED_CLASSES = {
    "user_requirement",
    "public_source_fact",
    "autonomous_design",
    "pending_unknown",
}
_FALSE_PROJECT_STATES = (
    "visualization_ready",
    "visual_correspondence_verified",
    "pde_validated",
    "paper_exact",
    "openai_field_identified",
    "blowup_proved",
)


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _load_json(path: Path) -> dict[str, Any]:
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise ValueError(f"expected JSON object in {path}")
    return data


def load_scope() -> dict[str, Any]:
    return _load_json(_repo_root() / "configs" / "st048s_post_piola_swirl_gain_governance.json")


def load_constraints() -> dict[str, Any]:
    return _load_json(_repo_root() / "configs" / "constraints.json")


def load_project_status() -> dict[str, Any]:
    return _load_json(_repo_root() / "project_status.json")


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _require_equal(actual: Any, expected: Any, label: str) -> None:
    _require(actual == expected, f"{label} drifted: {actual!r} != {expected!r}")


def _require_close(actual: float, expected: float, label: str, tol: float = 1e-12) -> None:
    _require(abs(float(actual) - expected) <= tol, f"{label} drifted: {actual!r} != {expected!r}")


def audit_st048s_swirl_gain_governance(
    scope: Mapping[str, Any] | None = None,
    constraints: Mapping[str, Any] | None = None,
    project_status: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Audit PR #428 semantics without selecting or promoting its child."""

    scope = deepcopy(dict(load_scope() if scope is None else scope))
    constraints = deepcopy(dict(load_constraints() if constraints is None else constraints))
    project_status = deepcopy(dict(load_project_status() if project_status is None else project_status))

    _require_equal(scope.get("task_id"), TASK_ID, "task id")
    upstream = scope.get("governed_upstream", {})
    _require_equal(upstream.get("pr"), UPSTREAM_PR, "upstream PR")
    _require_equal(upstream.get("head_sha"), UPSTREAM_HEAD, "upstream head")
    _require_equal(upstream.get("status"), "open_unintegrated_stacked_experimental_evidence", "upstream status")
    _require_equal(upstream.get("parent_pr"), PARENT_PR, "parent PR")
    _require_equal(upstream.get("parent_head_sha"), PARENT_HEAD, "parent head")
    _require_equal(upstream.get("candidate_materialized"), False, "upstream materialization")

    representation = scope.get("representation", {})
    _require_equal(
        representation.get("parent_temporal_piola_schedule"),
        "beta_gamma(t)=0.075+0.025*(4*t-2)^2",
        "temporal Piola schedule",
    )
    _require_equal(representation.get("swirl_transform"), "u_theta -> (1+kappa)*u_theta", "swirl transform")
    _require_equal(representation.get("screen_kappa_grid"), [0.0, 0.025, 0.05, 0.075, 0.1, 0.15], "kappa grid")
    _require_equal(representation.get("smallest_preregistered_clean_winding_crossing_kappa"), 0.05, "diagnostic kappa crossing")
    _require_close(representation.get("diagnostic_crossing_multiplier", 0.0), 1.05, "diagnostic multiplier")
    _require_close(
        representation.get("diagnostic_crossing_common_reference_energy_scale", 0.0),
        0.9726627718,
        "diagnostic common scale",
    )
    for key in (
        "kappa_constant_in_time_in_screen",
        "common_energy_scale_constant_in_time_per_kappa",
        "positive_finite_constant_swirl_gain_preserves_axisymmetric_divergence_structure",
        "positive_finite_constant_swirl_gain_preserves_parent_u_theta_O_r_axis_regular_order",
    ):
        _require_equal(representation.get(key), True, key)
    for key in (
        "navier_stokes_invariant_under_swirl_gain",
        "navier_stokes_invariant_under_common_energy_scale",
    ):
        _require_equal(representation.get(key), False, key)

    derivative = scope.get("derivative_and_axis_contract", {})
    _require_equal(derivative.get("constant_kappa_adds_kappa_prime_u_theta_term"), False, "constant kappa derivative")
    for key in (
        "temporal_piola_map_motion_terms_still_required",
        "future_time_dependent_kappa_requires_kappa_prime_u_theta_term",
        "future_time_dependent_scale_requires_scale_prime_u_term",
        "cylindrical_to_cartesian_export_must_use_regular_axis_extension",
        "near_axis_finite_cartesian_velocity_required",
        "fresh_numerical_divergence_recheck_required_after_materialization",
    ):
        _require_equal(derivative.get(key), True, key)
    for key in (
        "parent_pressure_force_momentum_receipt_transfer_allowed",
        "parent_pde_validation_transfer_allowed",
    ):
        _require_equal(derivative.get(key), False, key)

    allowed = set(scope.get("allowed_source_classes", ()))
    _require_equal(allowed, _ALLOWED_CLASSES, "source vocabulary")
    source = scope.get("source_classification", {})
    _require(source and set(source.values()) <= _ALLOWED_CLASSES, "invalid source classification")
    _require_equal(source.get("callable_saveable_velocity_delivery"), "user_requirement", "delivery source class")
    _require_equal(source.get("public_openai_visualization_as_observable_target"), "public_source_fact", "public target source class")
    for key in (
        "st048s_parent_candidate",
        "temporal_piola_family",
        "beta_base_and_gamma_0p025_schedule",
        "post_piola_swirl_gain_transform",
        "kappa_screen_grid",
        "kappa_0p05_clean_winding_crossing",
        "positive_common_reference_energy_scale",
        "three_column_sensitivity_screen",
    ):
        _require_equal(source.get(key), "autonomous_design", f"source class {key}")
    for key in (
        "openai_hidden_swirl_gain_or_coefficient",
        "openai_hidden_coordinate_time_camera_alignment",
    ):
        _require_equal(source.get(key), "pending_unknown", f"source class {key}")

    evidence = scope.get("screen_evidence_only", {})
    _require_equal(evidence.get("kappa_0p05_q90_moves_33_grid_bin"), False, "q90 grid-bin evidence")
    _require_equal(evidence.get("kappa_0p05_q99_moves_33_grid_bin"), False, "q99 grid-bin evidence")
    _require(evidence.get("kappa_0p05_independent_cartesian_fd_divergence_max", 1.0) < 1e-8, "diagnostic divergence evidence drifted")
    _require_equal(evidence.get("three_column_sensitivity_rank"), 3, "sensitivity rank")
    _require(evidence.get("three_column_sensitivity_condition_number", 1e9) < 4.0, "sensitivity condition evidence drifted")
    for key in (
        "selects_production_kappa",
        "screen_grid_defines_production_kappa_bound",
        "winding_proxy_crossing_is_visual_correspondence",
        "rank_three_sensitivity_is_pde_leverage_proof",
        "divergence_support_energy_core_preflights_are_full_pde_acceptance",
    ):
        _require_equal(evidence.get(key), False, key)

    routing = scope.get("lineage_and_routing_scope", {})
    _require_equal(routing.get("live_route_must_remain"), LIVE_ROUTE, "scope live route")
    _require_equal(routing.get("live_integrated_mode_must_remain"), LIVE_MODE, "scope live mode")
    _require_equal(routing.get("swirl_gain_evidence_is_open_unintegrated_stacked"), True, "stacked evidence state")
    _require_equal(routing.get("branch_tip_is_scientific_identity"), False, "branch-tip identity")
    _require_equal(
        routing.get("exact_parent_commit_representation_identity_and_candidate_sha_required_before_promotion"),
        True,
        "immutable parent lineage guard",
    )
    _require_equal(routing.get("downstream_green_ci_promotes_parent_or_child"), False, "CI promotion guard")
    _require_equal(routing.get("may_silently_replace_canonical_velocity_or_live_route"), False, "silent route replacement")
    _require_equal(project_status.get("active_scientific_route"), LIVE_ROUTE, "project live route")
    latest = project_status.get("latest_integrated_axial_cap_capacity", {})
    _require_equal(latest.get("mode"), LIVE_MODE, "project live mode")
    _require("not_materialized" in latest.get("status", ""), "live axial-cap unexpectedly materialized")
    _require_equal(latest.get("candidate_sha_created"), False, "live candidate SHA state")

    materialization = scope.get("future_materialization", {})
    for key in (
        "explicit_autonomous_beta_gamma_kappa_bounds_required",
        "explicit_nonzero_parameter_values_required",
        "explicit_normalization_rule_required",
        "new_representation_identity_required",
        "new_candidate_sha_required",
        "save_load_replay_required",
        "analytic_or_independently_checked_temporal_chain_rule_required",
        "axis_regular_cartesian_export_check_required",
        "rebuild_or_refit_restricted_pressure_force_required",
        "recheck_energy_core_axis_support_required",
        "fresh_full_momentum_divergence_required",
        "unused_acceptance_data_required",
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
        (snap.get("validation_seed"), validation.get("seed"), "validation seed"),
        (snap.get("held_out_points"), validation.get("held_out_points"), "held-out points"),
        (snap.get("validation_times"), validation.get("times"), "validation times"),
        (snap.get("derivative_steps"), validation.get("derivative_steps"), "derivative steps"),
        (snap.get("divergence_max"), thresholds.get("divergence_max"), "divergence max"),
        (snap.get("divergence_L2"), thresholds.get("divergence_L2"), "divergence L2"),
        (snap.get("pde_residual_max"), thresholds.get("pde_residual_max"), "PDE max"),
        (snap.get("pde_residual_L2"), thresholds.get("pde_residual_L2"), "PDE L2"),
    )
    for expected, actual, label in expected_pairs:
        _require_equal(actual, expected, label)
    _require("No residual-dependent basis" in forcing.get("restriction", ""), "free-force restriction weakened")

    truth = scope.get("truth_boundary", {})
    _require_equal(truth.get("existing_velocity_export_ready_may_remain_true"), True, "export independence")
    for key in (
        "swirl_gain_child_materialized",
        "production_kappa_selected",
        "visualization_ready",
        "visual_correspondence_verified",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    ):
        _require_equal(truth.get(key), False, key)
    states = project_status.get("states", {})
    _require_equal(states.get("velocity_export_ready"), True, "project export state")
    for key in _FALSE_PROJECT_STATES:
        _require_equal(states.get(key), False, f"project state {key}")

    return {
        "task_id": TASK_ID,
        "status": "governance_pass",
        "upstream_pr": UPSTREAM_PR,
        "upstream_status": "open_unintegrated_stacked_experimental_evidence",
        "kappa_0p05_is_diagnostic_only": True,
        "constant_kappa_adds_no_kappa_prime_term": True,
        "temporal_piola_chain_terms_still_required": True,
        "axis_regular_extension_required": True,
        "parent_pde_evidence_transfer_allowed": False,
        "live_route_preserved": True,
        "canonical_thresholds_unchanged": True,
        "velocity_export_ready": True,
        "visual_correspondence_verified": False,
        "pde_validated": False,
    }


def main() -> None:
    print(json.dumps(audit_st048s_swirl_gain_governance(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
