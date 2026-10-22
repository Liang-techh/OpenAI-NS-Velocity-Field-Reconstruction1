"""Fail-closed CR002 governance for the unintegrated axial-cap energy envelope.

This module changes no velocity, coefficient, pressure, forcing, sample, norm, or
threshold.  It records how Agent-7 PR #327 may be interpreted without promoting
its fixed-parent energy/morphology screen into a production coefficient, a PDE
result, or an OpenAI-field identity claim.
"""
from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
from typing import Any, Mapping

TASK_ID = "CR002-BIPOLAR-AXIAL-CAP-POLOIDAL-ENERGY-ENVELOPE-SCOPE-040"
MODE_NAME = "AXIAL_CAP_BANDED_C4_ODD_Z_POLOIDAL"
_ALLOWED_CLASSES = {
    "user_requirement",
    "public_source_fact",
    "autonomous_design",
    "pending_unknown",
}
_FALSE_TRUTH_STATES = (
    "axial_cap_poloidal_nonzero_child_materialized",
    "axial_cap_poloidal_coefficient_bound_selected",
    "axial_cap_poloidal_coefficient_selected",
    "candidate_selection_resolved",
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
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"expected JSON object in {path}")
    return value


def load_scope() -> dict[str, Any]:
    return _load_json(_repo_root() / "configs" / "bipolar_axial_cap_poloidal_energy_envelope_scope.json")


def load_constraints() -> dict[str, Any]:
    return _load_json(_repo_root() / "configs" / "constraints.json")


def load_project_status() -> dict[str, Any]:
    return _load_json(_repo_root() / "project_status.json")


def load_prior_representation_scope() -> dict[str, Any]:
    return _load_json(_repo_root() / "configs" / "bipolar_axial_cap_poloidal_representation_scope.json")


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _require_equal(actual: Any, expected: Any, label: str) -> None:
    _require(actual == expected, f"{label} drifted: {actual!r} != {expected!r}")


def _audit_cr001_snapshot(scope: Mapping[str, Any], constraints: Mapping[str, Any]) -> None:
    snapshot = scope.get("canonical_cr001_snapshot", {})
    domain = constraints["domain"]
    forcing = constraints["forcing"]
    nontriviality = constraints["nontriviality"]
    validation = constraints["validation"]
    thresholds = validation["thresholds"]

    _require_equal(snapshot.get("nu"), constraints["nu"], "nu")
    _require_equal(snapshot.get("physical_domain"), domain["physical"], "physical domain")
    _require_equal(snapshot.get("evaluation_box"), domain["evaluation_box"], "evaluation box")
    _require_equal(snapshot.get("support"), domain["support"], "support")
    _require_equal(snapshot.get("time_interval"), domain["time_interval"], "time interval")
    _require_equal(snapshot.get("force_mode"), forcing["mode"], "force mode")
    _require_equal(snapshot.get("force_bounds"), forcing["parameters"], "force bounds")
    _require_equal(snapshot.get("reference_energy"), nontriviality["reference_energy"], "reference energy")
    _require_equal(
        snapshot.get("reference_energy_abs_tolerance"),
        nontriviality["reference_energy_abs_tolerance"],
        "reference energy tolerance",
    )
    _require_equal(
        snapshot.get("validation_energy_range"),
        [
            nontriviality["minimum_energy_each_validation_time"],
            nontriviality["maximum_energy_each_validation_time"],
        ],
        "validation energy range",
    )
    _require_equal(snapshot.get("validation_seed"), validation["seed"], "validation seed")
    _require_equal(snapshot.get("held_out_points"), validation["held_out_points"], "held-out points")
    _require_equal(snapshot.get("derivative_steps"), validation["derivative_steps"], "derivative steps")
    for key in ("divergence_max", "divergence_L2", "pde_residual_max", "pde_residual_L2"):
        _require_equal(snapshot.get(key), thresholds[key], key)
    _require("No residual-dependent basis" in forcing["restriction"], "free residual-dependent forcing became allowed")


def audit_axial_cap_energy_envelope_scope(
    scope: Mapping[str, Any] | None = None,
    constraints: Mapping[str, Any] | None = None,
    project_status: Mapping[str, Any] | None = None,
    prior_representation_scope: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Audit the experimental axial-cap energy-envelope interpretation contract."""

    scope = deepcopy(dict(load_scope() if scope is None else scope))
    constraints = deepcopy(dict(load_constraints() if constraints is None else constraints))
    project_status = deepcopy(dict(load_project_status() if project_status is None else project_status))
    prior = deepcopy(dict(load_prior_representation_scope() if prior_representation_scope is None else prior_representation_scope))

    _require_equal(scope.get("task_id"), TASK_ID, "scope task id")
    upstream = scope.get("governed_upstream", {})
    _require_equal(upstream.get("pr"), 327, "upstream PR")
    _require_equal(upstream.get("task_id"), "CR003-BIPOLAR-AXIAL-CAP-ENERGY-ENVELOPE-046", "upstream task")
    _require_equal(upstream.get("head_sha"), "d4e740649bbd0825b0d8aaa336776fd54ee8f707", "upstream head")
    _require_equal(upstream.get("mode_name"), MODE_NAME, "mode name")
    _require_equal(
        upstream.get("source_candidate_sha256"),
        "7875214ad65f5e1ba4f7c5e362217f05ef47bc3d893e8641e8e40a85fc4e7609",
        "source candidate",
    )
    _require_equal(
        scope.get("representation_kind"),
        "experimental_unintegrated_axial_cap_fixed_parent_energy_envelope_audit_only",
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
        "zero_centered_energy_envelope_method",
        "quadrature_orders_96_192",
        "morphology_grid_size_33",
    ):
        _require_equal(source.get(key), "autonomous_design", f"source class {key}")
    _require_equal(source.get("openai_hidden_numerical_profile"), "pending_unknown", "hidden profile class")
    _require_equal(scope.get("evidence_relation_labels_are_source_classes"), False, "evidence/source separation")
    _require_equal(set(scope.get("evidence_relation_labels", {}).values()), {"project_derived_measurement"}, "evidence relation vocabulary")

    energy = scope.get("energy_envelope_evidence_only", {})
    _require_equal(energy.get("reference_time"), 0.25, "reference time")
    _require(0.999 <= float(energy.get("reference_baseline_energy", 0.0)) <= 1.001, "baseline energy evidence missing")
    _require(float(energy.get("energy_linear_coefficient", 0.0)) < 0.0, "axial-cap linear energy coefficient sign drifted")
    _require(float(energy.get("energy_quadratic_coefficient", 0.0)) > 0.2, "axial-cap quadratic energy coefficient missing")
    half_width = float(energy.get("largest_symmetric_zero_centered_diagnostic_half_width", 0.0))
    _require(0.06 < half_width < 0.07, "diagnostic fixed-parent envelope drifted")
    _require_equal(energy.get("quadrature_orders"), [96, 192], "quadrature orders")
    _require(0.0 <= float(energy.get("maximum_relative_quadratic_coefficient_change", 1.0)) < 1e-3, "energy refinement drifted")
    _require_equal(energy.get("previous_capacity_trial_abs_coefficient"), 0.25, "capacity trial")
    _require(float(energy.get("previous_trial_over_diagnostic_half_width_ratio", 0.0)) > 3.7, "capacity/envelope ratio drifted")
    _require_equal(energy.get("capacity_trial_inside_fixed_parent_energy_envelope"), False, "capacity-trial energy status")
    for key in (
        "diagnostic_envelope_is_selected_materialization_bound",
        "diagnostic_envelope_selects_coefficient",
        "diagnostic_envelope_is_complete_feasible_coefficient_interval",
        "diagnostic_exact_normalization_root_selects_coefficient",
        "fixed_parent_energy_admissibility_transfers_through_joint_renormalization",
        "validation_time_energy_range_is_full_candidate_acceptance",
        "quadrature_refinement_is_pde_validation",
        "previous_capacity_trial_invalidated_as_capacity_evidence",
    ):
        _require_equal(energy.get(key), False, key)

    morphology = scope.get("morphology_semantics", {})
    _require_equal(morphology.get("time"), 0.5, "morphology time")
    _require_equal(morphology.get("grid_size"), 33, "morphology grid size")
    _require_equal(morphology.get("baseline_q90_over_Zp"), morphology.get("energy_envelope_max_q90_over_Zp"), "q90 fixed-parent envelope reach")
    _require_equal(morphology.get("baseline_q99_over_Zp"), morphology.get("energy_envelope_max_q99_over_Zp"), "q99 fixed-parent envelope reach")
    _require_equal(morphology.get("energy_envelope_q90_moved_on_this_grid"), False, "q90 movement state")
    _require_equal(morphology.get("energy_envelope_q99_moved_on_this_grid"), False, "q99 movement state")
    _require(0.0 < float(morphology.get("maximum_axial_rms_relative_gain_approx", 0.0)) < 0.05, "axial RMS evidence drifted")
    _require(0.0 < float(morphology.get("maximum_outer_tail_fraction_remains_small_approx", 1.0)) < 0.01, "tail-fraction scale drifted")
    for key in (
        "zero_q90_q99_movement_proves_no_pde_leverage",
        "zero_q90_q99_movement_proves_mode_useless",
        "relative_tail_gain_proves_gross_tip_extension",
        "relative_tail_gain_establishes_visual_correspondence",
        "morphology_screen_is_momentum_residual_jacobian",
        "morphology_screen_selects_coefficient_or_sign",
    ):
        _require_equal(morphology.get(key), False, key)

    routing = scope.get("integration_routing_scope", {})
    _require_equal(routing.get("live_integrated_mode"), "COMPACT_C4_ODD_Z_POLOIDAL", "live integrated mode")
    _require_equal(routing.get("live_mode_nonzero_child_materialized"), False, "live materialization state")
    _require_equal(routing.get("live_next_integration_task_must_remain_materialize_existing_compact_poloidal_child"), True, "live next-task preservation")
    _require_equal(routing.get("axial_cap_energy_envelope_is_experimental_unintegrated"), True, "axial-cap integration state")
    _require_equal(routing.get("axial_cap_energy_envelope_may_replace_live_next_task_without_explicit_integration_decision"), False, "routing promotion guard")
    _require_equal(routing.get("axial_cap_joint_renormalization_may_be_tested_as_later_pre_materialization_experiment"), True, "later renorm experiment allowance")
    _require_equal(routing.get("axial_cap_joint_renormalization_is_candidate_promotion"), False, "renorm promotion state")
    _require_equal(routing.get("additional_axial_cap_or_higher_odd_q_basis_growth_supported_by_this_result"), False, "basis-growth guard")
    _require_equal(routing.get("reference_energy_gate_may_be_relaxed_for_tip_reach"), False, "energy threshold guard")

    prior_routing = prior.get("integration_routing_scope", {})
    _require_equal(prior_routing.get("axial_cap_capacity_is_experimental_unintegrated"), True, "prior axial-cap integration state")
    _require_equal(prior_routing.get("live_next_integration_task_must_remain_materialize_existing_compact_poloidal_child"), True, "prior live next-task state")
    _require(MODE_NAME in prior.get("governed_upstream", {}).get("mode_name", ""), "prior representation scope no longer binds axial-cap mode")

    _require_equal(project_status.get("states", {}).get("velocity_export_ready"), True, "existing velocity export state")
    _require_equal(project_status.get("states", {}).get("pde_validated"), False, "project PDE state")
    _require("COMPACT_C4_ODD_Z_POLOIDAL" in project_status.get("next_integration_task", ""), "project next task no longer matches live compact-poloidal materialization route")

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
        "joint_renormalization_requires_fresh_nonlinear_pde_validation",
        "model_selection_data_may_not_be_reused_as_acceptance_data",
    ):
        _require_equal(materialization.get(key), True, key)

    _audit_cr001_snapshot(scope, constraints)

    truth = scope.get("truth_boundary", {})
    _require_equal(truth.get("velocity_export_ready_may_remain_true_for_existing_candidate"), True, "velocity export independence")
    for key in _FALSE_TRUTH_STATES:
        _require_equal(truth.get(key), False, f"truth state {key}")

    return {
        "task_id": TASK_ID,
        "governed_upstream_pr": 327,
        "mode_name": MODE_NAME,
        "fixed_parent_diagnostic_half_width": half_width,
        "gross_q90_q99_reach_moved_inside_envelope": False,
        "live_next_task_preserved": True,
        "canonical_cr001_unchanged": True,
        "truth_boundary_preserved": True,
    }


if __name__ == "__main__":
    print(json.dumps(audit_axial_cap_energy_envelope_scope(), indent=2, sort_keys=True))
