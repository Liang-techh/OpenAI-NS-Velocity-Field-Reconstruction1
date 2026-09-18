"""Fail-closed CR002 governance for the open ST048-S temporal Piola screen.

This module changes no velocity values.  It makes the time-dependent coordinate
map's u_t chain rule explicit and keeps PR #418's morphology/capacity evidence
separate from production selection, the live scientific route, and PDE claims.
"""
from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
from typing import Any, Mapping

TASK_ID = "CR002-ST048S-TEMPORAL-PIOLA-GOVERNANCE-064"
UPSTREAM_PR = 418
UPSTREAM_HEAD = "ae10c3ccacd16e1f5e0054502f1dc5a6f68320b5"
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
    return _load_json(_repo_root() / "configs" / "st048s_temporal_piola_warp_governance.json")


def load_constraints() -> dict[str, Any]:
    return _load_json(_repo_root() / "configs" / "constraints.json")


def load_project_status() -> dict[str, Any]:
    return _load_json(_repo_root() / "project_status.json")


def beta_gamma(t: float, gamma: float) -> float:
    return 0.075 + gamma * (4.0 * t - 2.0) ** 2


def beta_gamma_prime(t: float, gamma: float) -> float:
    return 8.0 * gamma * (4.0 * t - 2.0)


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _require_equal(actual: Any, expected: Any, label: str) -> None:
    _require(actual == expected, f"{label} drifted: {actual!r} != {expected!r}")


def audit_st048s_temporal_piola_governance(
    scope: Mapping[str, Any] | None = None,
    constraints: Mapping[str, Any] | None = None,
    project_status: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Audit PR #418 semantics without selecting or promoting its child."""

    scope = deepcopy(dict(load_scope() if scope is None else scope))
    constraints = deepcopy(dict(load_constraints() if constraints is None else constraints))
    project_status = deepcopy(dict(load_project_status() if project_status is None else project_status))

    _require_equal(scope.get("task_id"), TASK_ID, "task id")
    upstream = scope.get("governed_upstream", {})
    _require_equal(upstream.get("pr"), UPSTREAM_PR, "upstream PR")
    _require_equal(upstream.get("head_sha"), UPSTREAM_HEAD, "upstream head")
    _require_equal(upstream.get("status"), "open_unintegrated_stacked_experimental_evidence", "upstream status")
    _require_equal(upstream.get("parent_pr"), 407, "upstream parent PR")
    _require_equal(upstream.get("candidate_materialized"), False, "upstream materialization")

    representation = scope.get("representation", {})
    _require("H(z,t)" in representation.get("map", ""), "time-dependent map missing")
    _require("H_z" in representation.get("velocity", ""), "Piola velocity missing")
    _require_equal(representation.get("beta_schedule"), "beta_gamma(t)=0.075+gamma*(4*t-2)^2", "beta schedule")
    _require_equal(representation.get("beta_prime"), "8*gamma*(4*t-2)", "beta prime")
    _require_equal(representation.get("beta_second"), "32*gamma", "beta second")
    _require_equal(representation.get("screen_gamma_grid"), [0.0, 0.025, 0.05, 0.075, 0.125], "gamma grid")
    expected_beta = [beta_gamma(t, 0.025) for t in (0.25, 0.5, 0.75)]
    expected_prime = [beta_gamma_prime(t, 0.025) for t in (0.25, 0.5, 0.75)]
    _require_equal(representation.get("gamma_0p025_beta_at_t025_t05_t075"), expected_beta, "gamma crossing beta values")
    _require_equal(representation.get("gamma_0p025_beta_prime_at_t025_t05_t075"), expected_prime, "gamma crossing beta-prime values")
    _require_equal(representation.get("screen_common_energy_scale_is_constant_in_time_per_gamma"), True, "screen scale time dependence")
    _require_equal(representation.get("navier_stokes_invariant_under_time_dependent_warp"), False, "temporal warp NS invariance")
    _require_equal(representation.get("navier_stokes_invariant_under_common_energy_scale"), False, "common scale NS invariance")

    derivative = scope.get("time_derivative_contract", {})
    _require("H_zt*u_i" in derivative.get("horizontal_formula", ""), "horizontal H_zt chain term missing")
    _require("H_t*u_i,z" in derivative.get("horizontal_formula", ""), "horizontal map-motion chain term missing")
    _require("H_t*u_z,z" in derivative.get("axial_formula", ""), "axial map-motion chain term missing")
    for key in (
        "H_t_contains_beta_prime",
        "H_zt_contains_beta_prime",
        "future_time_dependent_scale_requires_scale_prime_term",
        "endpoint_time_derivative_must_respect_registered_time_window",
    ):
        _require_equal(derivative.get(key), True, key)
    for key in (
        "omitting_map_motion_terms_allowed",
        "parent_static_warp_u_t_receipt_transfer_allowed",
        "parent_pressure_force_momentum_receipt_transfer_allowed",
        "constant_screen_scale_adds_scale_prime_term",
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
        "piola_warp_family",
        "beta_base_0p075",
        "quadratic_even_time_schedule",
        "gamma_screen_grid",
        "gamma_0p025_clean_crossing",
        "positive_common_reference_energy_scale",
        "sensitivity_rank_screen",
    ):
        _require_equal(source.get(key), "autonomous_design", f"source class {key}")
    for key in (
        "openai_hidden_time_dependent_coordinate_map",
        "openai_hidden_frame_time_camera_alignment",
    ):
        _require_equal(source.get(key), "pending_unknown", f"source class {key}")

    evidence = scope.get("screen_evidence_only", {})
    _require_equal(evidence.get("smallest_preregistered_clean_endpoint_crossing_gamma"), 0.025, "diagnostic gamma crossing")
    _require(abs(evidence.get("gamma_0p025_common_scale", 0.0) - 0.998791257501046) < 1e-15, "common scale evidence drifted")
    _require_equal(evidence.get("gamma_0p025_q99_moves"), False, "q99 evidence")
    _require_equal(evidence.get("sensitivity_rank"), 2, "sensitivity rank")
    _require(evidence.get("sensitivity_condition_number", 1e9) < 4.0, "sensitivity conditioning evidence drifted")
    for key in (
        "selects_production_gamma",
        "screen_grid_defines_production_gamma_bound",
        "morphology_crossing_is_visual_correspondence",
        "rank_two_sensitivity_is_pde_leverage_proof",
        "energy_support_core_divergence_preflights_are_full_pde_acceptance",
    ):
        _require_equal(evidence.get(key), False, key)

    routing = scope.get("integration_routing_scope", {})
    _require_equal(routing.get("live_route_must_remain"), LIVE_ROUTE, "scope live route")
    _require_equal(routing.get("live_integrated_mode_must_remain"), LIVE_MODE, "scope live mode")
    _require_equal(routing.get("temporal_piola_evidence_is_open_unintegrated_stacked"), True, "stacked status")
    _require_equal(routing.get("may_silently_replace_live_route"), False, "silent route replacement")
    _require_equal(routing.get("explicit_routing_decision_required_before_materialization"), True, "routing decision guard")
    _require_equal(routing.get("upstream_parent_lineage_must_be_resolved_before_promotion"), True, "lineage guard")
    _require_equal(project_status.get("active_scientific_route"), LIVE_ROUTE, "project live route")
    latest = project_status.get("latest_integrated_axial_cap_capacity", {})
    _require_equal(latest.get("mode"), LIVE_MODE, "project live mode")
    _require("not_materialized" in latest.get("status", ""), "live axial-cap unexpectedly materialized")
    _require_equal(latest.get("candidate_sha_created"), False, "live candidate SHA state")

    materialization = scope.get("future_materialization", {})
    for key in (
        "explicit_autonomous_beta_and_gamma_bounds_required",
        "explicit_nonzero_beta_gamma_values_required",
        "explicit_normalization_rule_required",
        "new_representation_identity_required",
        "new_candidate_sha_required",
        "save_load_replay_required",
        "analytic_or_independently_checked_u_t_chain_rule_required",
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
        "temporal_piola_child_materialized",
        "production_gamma_selected",
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
        "gamma_0p025_is_diagnostic_only": True,
        "beta_prime_endpoint_abs_at_gamma_0p025": 0.2,
        "u_t_map_motion_terms_required": True,
        "parent_pde_evidence_transfer_allowed": False,
        "live_route_preserved": True,
        "canonical_thresholds_unchanged": True,
        "velocity_export_ready": True,
        "visual_correspondence_verified": False,
        "pde_validated": False,
    }


def main() -> None:
    print(json.dumps(audit_st048s_temporal_piola_governance(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
