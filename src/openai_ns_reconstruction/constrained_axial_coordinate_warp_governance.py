"""Fail-closed CR002 governance for open axial-coordinate-warp evidence.

This module changes no velocity value. It keeps PR #367's autonomous Piola
coordinate-warp capacity evidence separate from the currently integrated
axial-cap materialization route, public-source facts, and PDE/visual claims.
"""
from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
from typing import Any, Mapping

TASK_ID = "CR002-AXIAL-COORDINATE-WARP-GOVERNANCE-057"
UPSTREAM_TASK_ID = "CR003-AXIAL-COORDINATE-WARP-SCREEN-056"
UPSTREAM_PR = 367
UPSTREAM_HEAD = "25797f27619dcdf8798e0128b764d6f4eccba7f7"
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
    "physical_support_validated",
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
    return _load_json(_repo_root() / "configs" / "axial_coordinate_warp_governance.json")


def load_constraints() -> dict[str, Any]:
    return _load_json(_repo_root() / "configs" / "constraints.json")


def load_project_status() -> dict[str, Any]:
    return _load_json(_repo_root() / "project_status.json")


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _require_equal(actual: Any, expected: Any, label: str) -> None:
    _require(actual == expected, f"{label} drifted: {actual!r} != {expected!r}")


def audit_axial_coordinate_warp_governance(
    scope: Mapping[str, Any] | None = None,
    constraints: Mapping[str, Any] | None = None,
    project_status: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Audit warp capacity semantics without promoting the open experiment."""

    scope = deepcopy(dict(load_scope() if scope is None else scope))
    constraints = deepcopy(dict(load_constraints() if constraints is None else constraints))
    project_status = deepcopy(dict(load_project_status() if project_status is None else project_status))

    _require_equal(scope.get("task_id"), TASK_ID, "task id")
    upstream = scope.get("governed_upstream", {})
    _require_equal(upstream.get("pr"), UPSTREAM_PR, "upstream PR")
    _require_equal(upstream.get("task_id"), UPSTREAM_TASK_ID, "upstream task")
    _require_equal(upstream.get("head_sha"), UPSTREAM_HEAD, "upstream head")
    _require_equal(upstream.get("status"), "open_unintegrated_experimental_evidence", "upstream status")
    _require_equal(upstream.get("candidate_materialized"), False, "upstream materialization")

    representation = scope.get("representation", {})
    _require_equal(
        representation.get("kind"),
        "autonomous_support_preserving_piola_coordinate_warp_capacity_screen",
        "representation kind",
    )
    _require("h_beta(z)" in representation.get("map_definition", ""), "warp map missing")
    _require("h'(z)" in representation.get("velocity_definition", ""), "Piola velocity map missing")
    _require("div(u_beta)=h'(z)" in representation.get("divergence_identity", ""), "divergence identity missing")
    _require_equal(representation.get("screen_beta_is_time_independent"), True, "screen beta time dependence")
    _require_equal(representation.get("spatial_derivatives_change_under_warp"), True, "spatial derivative semantics")
    _require_equal(representation.get("navier_stokes_invariant_under_warp"), False, "warp NS invariance")
    _require_equal(
        representation.get("navier_stokes_invariant_under_common_energy_scale"),
        False,
        "common-scale NS invariance",
    )

    allowed = set(scope.get("allowed_source_classes", ()))
    _require_equal(allowed, _ALLOWED_CLASSES, "source vocabulary")
    source = scope.get("source_classification", {})
    _require(source and set(source.values()) <= _ALLOWED_CLASSES, "invalid source classification")
    _require_equal(source.get("callable_saveable_velocity_delivery"), "user_requirement", "delivery class")
    _require_equal(source.get("eq45_backbone"), "public_source_fact", "Eq45 class")
    for key in (
        "piola_coordinate_warp_family",
        "quartic_collar_fade_power_4",
        "screen_beta_grid_0_to_0p5",
        "beta_0p40_diagnostic_crossing",
        "positive_common_energy_renormalization",
    ):
        _require_equal(source.get(key), "autonomous_design", f"source class {key}")
    for key in (
        "openai_hidden_coordinate_or_aspect_map",
        "openai_hidden_frame_time_camera_alignment",
    ):
        _require_equal(source.get(key), "pending_unknown", f"source class {key}")

    evidence = scope.get("screen_evidence_only", {})
    _require_equal(evidence.get("beta_grid"), [0.0, 0.1, 0.2, 0.3, 0.4, 0.5], "beta grid")
    _require_equal(evidence.get("smallest_sampled_preflight_passing_multi_time_q90_mover"), 0.4, "diagnostic crossing")
    _require(0.0 < evidence.get("beta_0p40_raw_reference_energy", 0.0) < 1.0, "raw energy evidence missing")
    _require(evidence.get("beta_0p40_common_velocity_scale", 0.0) > 1.0, "common scale evidence missing")
    _require(abs(evidence.get("beta_0p40_jointly_renormalized_reference_energy", 0.0) - 1.0) < 1e-12, "renormalized energy drifted")
    _require_equal(evidence.get("beta_0p40_q90_gain_over_Zp_by_time"), [0.15234375, 0.15234375, 0.076171875], "q90 evidence")
    _require_equal(evidence.get("beta_0p40_q99_gain_over_Zp_by_time"), [0.15234375, 0.15234375, 0.076171875], "q99 evidence")
    _require(evidence.get("beta_0p40_coordinate_jacobian_min", 0.0) > 0.0, "positive Jacobian evidence missing")
    _require(0.0 < evidence.get("beta_0p40_cartesian_fd_divergence_max_abs", 1.0) < 1e-7, "FD divergence evidence drifted")
    for key in (
        "beta_0p40_selects_production_value",
        "screen_beta_grid_defines_production_bound",
        "positive_jacobian_establishes_openai_coordinate_map",
        "support_endpoint_preservation_is_physical_support_validation",
        "piola_identity_is_registered_full_candidate_divergence_acceptance",
        "fd_divergence_is_full_pde_acceptance",
        "target_free_morphology_is_visual_correspondence",
        "lower_collar_or_radial_cost_is_candidate_superiority",
        "common_scale_is_bookkeeping_only",
        "parent_pressure_force_or_pde_evidence_transfer_allowed",
    ):
        _require_equal(evidence.get(key), False, key)

    derivatives = scope.get("derivative_semantics", {})
    _require_equal(derivatives.get("current_screen_beta_constant_in_time"), True, "current beta schedule")
    _require_equal(derivatives.get("current_screen_adds_beta_prime_term_to_u_t"), False, "current u_t semantics")
    _require_equal(derivatives.get("warped_spatial_chain_rule_must_be_revalidated"), True, "spatial revalidation")
    _require_equal(derivatives.get("future_time_dependent_beta_requires_fresh_u_t_chain_rule"), True, "future u_t guard")
    _require_equal(derivatives.get("future_time_dependent_beta_may_inherit_static_warp_pde_evidence"), False, "future PDE inheritance")

    routing = scope.get("integration_routing_scope", {})
    _require_equal(routing.get("live_route_must_remain"), LIVE_ROUTE, "scope live route")
    _require_equal(routing.get("live_integrated_mode_must_remain"), LIVE_MODE, "scope live mode")
    _require_equal(routing.get("coordinate_warp_evidence_is_open_unintegrated"), True, "warp integration state")
    _require_equal(routing.get("coordinate_warp_may_silently_replace_live_route"), False, "silent route replacement")
    _require_equal(routing.get("explicit_routing_decision_required_before_warp_materialization"), True, "routing decision guard")
    _require_equal(
        routing.get("additional_capacity_basis_growth_allowed_before_live_materialization_screen"),
        False,
        "capacity-growth guard",
    )
    _require_equal(project_status.get("active_scientific_route"), LIVE_ROUTE, "project live route")
    latest = project_status.get("latest_integrated_axial_cap_capacity", {})
    _require_equal(latest.get("mode"), LIVE_MODE, "project live mode")
    _require("not_materialized" in latest.get("status", ""), "live child unexpectedly materialized")
    _require_equal(latest.get("production_coefficient_bound_selected"), False, "live bound selection")
    _require_equal(latest.get("production_coefficient_value_selected"), False, "live value selection")
    _require_equal(latest.get("candidate_sha_created"), False, "live SHA state")
    next_task = project_status.get("next_integration_task", "")
    for token in (
        LIVE_MODE,
        "explicit autonomous coefficient bound/value",
        "new representation identity/candidate SHA",
        "fresh full-momentum/divergence validation",
        "do not add another capacity basis first",
    ):
        _require(token in next_task, f"live next-integration token missing: {token}")
    avoided = set(project_status.get("integration_policy", {}).get("avoid", ()))
    _require(
        "additional_capacity_basis_growth_before_a_current_integrated_direction_is_materialized_and_screened" in avoided,
        "live basis-growth block missing",
    )

    materialization = scope.get("future_warp_materialization", {})
    for key in (
        "explicit_autonomous_beta_bound_required",
        "explicit_nonzero_beta_value_required",
        "explicit_normalization_rule_required",
        "new_representation_family_identity_required",
        "new_candidate_sha_required",
        "candidate_save_load_replay_required",
        "recheck_reference_energy",
        "recheck_validation_time_energy",
        "recheck_core_sign_nontriviality",
        "recheck_axis_support_regularity",
        "rebuild_or_refit_restricted_pressure_force_contract",
        "fresh_full_per_component_momentum_divergence_required",
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
        (snap.get("validation_energy_range"), [nontriviality.get("minimum_energy_each_validation_time"), nontriviality.get("maximum_energy_each_validation_time")], "validation energy range"),
        (snap.get("validation_seed"), validation.get("seed"), "validation seed"),
        (snap.get("held_out_points"), validation.get("held_out_points"), "held-out points"),
        (snap.get("validation_times"), validation.get("times"), "validation times"),
        (snap.get("derivative_steps"), validation.get("derivative_steps"), "derivative steps"),
        (snap.get("divergence_max"), thresholds.get("divergence_max"), "divergence max"),
        (snap.get("divergence_L2"), thresholds.get("divergence_L2"), "divergence L2"),
        (snap.get("pde_residual_max"), thresholds.get("pde_residual_max"), "PDE max"),
        (snap.get("pde_residual_L2"), thresholds.get("pde_residual_L2"), "PDE L2"),
        (snap.get("boundary_velocity_max"), thresholds.get("boundary_velocity_max"), "boundary velocity"),
        (snap.get("boundary_pressure_max"), thresholds.get("boundary_pressure_max"), "boundary pressure"),
    )
    for expected, actual, label in expected_pairs:
        _require_equal(actual, expected, label)
    _require("No residual-dependent basis" in forcing.get("restriction", ""), "free-force restriction weakened")

    truth = scope.get("truth_boundary", {})
    _require_equal(truth.get("existing_velocity_export_ready_may_remain_true"), True, "export independence")
    for key in (
        "coordinate_warp_child_materialized",
        "coordinate_warp_parameter_selected",
        "candidate_selection_resolved",
        "physical_support_validated",
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
        "upstream_status": "open_unintegrated_experimental_evidence",
        "beta_0p40_is_diagnostic_only": True,
        "warp_is_navier_stokes_invariance": False,
        "parent_pde_evidence_transfer_allowed": False,
        "live_route_preserved": True,
        "canonical_thresholds_unchanged": True,
        "velocity_export_ready": True,
        "visual_correspondence_verified": False,
        "pde_validated": False,
    }


def main() -> None:
    print(json.dumps(audit_axial_coordinate_warp_governance(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
