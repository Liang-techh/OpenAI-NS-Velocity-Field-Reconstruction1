"""Fail-closed CR002 governance for the compact-poloidal energy envelope.

This module changes no velocity, coefficient, force, pressure, sample, norm, or
threshold. It keeps Agent-7 PR #309's fixed-parent energy/morphology audit from
being promoted into a production coefficient, a PDE result, or an OpenAI-field
identity claim.
"""
from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
from typing import Any, Mapping

TASK_ID = "CR002-BIPOLAR-COMPACT-POLOIDAL-ENERGY-ENVELOPE-SCOPE-039"
MODE_NAME = "COMPACT_C4_ODD_Z_POLOIDAL"
_ALLOWED_CLASSES = {
    "user_requirement",
    "public_source_fact",
    "autonomous_design",
    "pending_unknown",
}
_FALSE_TRUTH_STATES = (
    "compact_poloidal_nonzero_child_materialized",
    "compact_poloidal_coefficient_bound_selected",
    "compact_poloidal_coefficient_selected",
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
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise ValueError(f"expected JSON object in {path}")
    return data


def load_scope() -> dict[str, Any]:
    return _load_json(_repo_root() / "configs" / "bipolar_compact_poloidal_energy_envelope_scope.json")


def load_constraints() -> dict[str, Any]:
    return _load_json(_repo_root() / "configs" / "constraints.json")


def load_project_status() -> dict[str, Any]:
    return _load_json(_repo_root() / "project_status.json")


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _require_equal(actual: Any, expected: Any, label: str) -> None:
    _require(actual == expected, f"{label} drifted: {actual!r} != {expected!r}")


def audit_compact_poloidal_energy_envelope_scope(
    scope: Mapping[str, Any] | None = None,
    constraints: Mapping[str, Any] | None = None,
    project_status: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Audit fixed-parent energy-envelope evidence without selecting a child."""

    scope = deepcopy(dict(load_scope() if scope is None else scope))
    constraints = deepcopy(dict(load_constraints() if constraints is None else constraints))
    project_status = deepcopy(dict(load_project_status() if project_status is None else project_status))

    _require_equal(scope.get("task_id"), TASK_ID, "scope task id")
    upstream = scope.get("governed_upstream", {})
    _require_equal(upstream.get("pr"), 309, "upstream PR")
    _require_equal(
        upstream.get("task_id"),
        "CR003-BIPOLAR-COMPACT-POLOIDAL-ENERGY-ENVELOPE-044",
        "upstream task",
    )
    _require_equal(
        upstream.get("head_sha"),
        "1f8be88f9b6f367298092f14ad077f26259cc4ef",
        "upstream head",
    )
    _require_equal(upstream.get("mode_name"), MODE_NAME, "mode name")
    _require_equal(
        upstream.get("source_candidate_sha256"),
        "7875214ad65f5e1ba4f7c5e362217f05ef47bc3d893e8641e8e40a85fc4e7609",
        "upstream source candidate",
    )
    _require_equal(
        scope.get("representation_kind"),
        "integrated_autonomous_compact_poloidal_response_energy_envelope_audit_only",
        "representation kind",
    )

    allowed = set(scope.get("allowed_source_classes", ()))
    _require_equal(allowed, _ALLOWED_CLASSES, "source-class vocabulary")
    source = scope.get("source_classification", {})
    _require(source and set(source.values()) <= _ALLOWED_CLASSES, "invalid source classification")
    _require_equal(source.get("callable_velocity_delivery"), "user_requirement", "velocity delivery class")
    _require_equal(source.get("eq45_backbone"), "public_source_fact", "Eq45 backbone class")
    for key in (
        "compact_poloidal_mode",
        "zero_centered_energy_envelope_method",
        "quadrature_orders_96_192",
        "morphology_grid_size_25",
    ):
        _require_equal(source.get(key), "autonomous_design", f"source class {key}")
    _require_equal(source.get("openai_hidden_numerical_profile"), "pending_unknown", "hidden profile class")
    _require_equal(scope.get("evidence_relation_labels_are_source_classes"), False, "evidence/source separation")
    evidence_labels = set(scope.get("evidence_relation_labels", {}).values())
    _require_equal(evidence_labels, {"project_derived_measurement"}, "evidence relation vocabulary")

    energy = scope.get("energy_envelope_evidence_only", {})
    _require_equal(energy.get("reference_time"), 0.25, "reference time")
    _require(0.999 <= energy.get("reference_baseline_energy", 0.0) <= 1.001, "baseline energy evidence missing")
    _require(energy.get("energy_linear_coefficient", 0.0) != 0.0, "linear energy coefficient missing")
    _require(energy.get("energy_quadratic_coefficient", 0.0) > 0.0, "quadratic energy coefficient missing")
    half_width = float(energy.get("largest_symmetric_zero_centered_diagnostic_half_width", 0.0))
    _require(0.0 < half_width < 0.25, "diagnostic energy envelope drifted")
    _require_equal(energy.get("quadrature_orders"), [96, 192], "quadrature orders")
    _require(0.0 <= energy.get("maximum_relative_quadratic_coefficient_change", 1.0) < 1e-6, "energy refinement drifted")
    _require_equal(energy.get("previous_capacity_trial_abs_coefficient"), 0.25, "previous capacity amplitude")
    _require(energy.get("previous_trial_over_diagnostic_half_width_ratio", 0.0) > 9.0, "capacity/envelope ratio drifted")
    for key in (
        "diagnostic_envelope_is_selected_materialization_bound",
        "diagnostic_envelope_selects_coefficient",
        "diagnostic_envelope_is_complete_feasible_coefficient_interval",
        "diagnostic_exact_normalization_root_selects_coefficient",
        "fixed_parent_energy_admissibility_transfers_through_joint_renormalization",
        "validation_time_energy_range_is_full_candidate_acceptance",
        "quadrature_refinement_is_pde_validation",
        "previous_capacity_trial_invalidated_by_energy_envelope",
    ):
        _require_equal(energy.get(key), False, key)

    morphology = scope.get("morphology_semantics", {})
    _require_equal(morphology.get("time"), 0.5, "morphology time")
    _require_equal(morphology.get("grid_size"), 25, "morphology grid")
    _require(morphology.get("axial_vorticity_rms_span_over_Zp", 1.0) < 0.001, "axial morphology span drifted")
    _require(morphology.get("radial_vorticity_rms_span_over_Rp", 1.0) < 0.001, "radial morphology span drifted")
    _require_equal(morphology.get("axial_q90_span_over_Zp"), 0.0, "q90 span")
    _require_equal(morphology.get("axial_q99_span_over_Zp"), 0.0, "q99 span")
    for key in (
        "small_span_proves_no_pde_leverage",
        "zero_q90_q99_span_proves_mode_useless",
        "morphology_span_is_momentum_residual_jacobian",
        "morphology_span_is_visual_correspondence",
        "morphology_span_selects_coefficient_or_sign",
    ):
        _require_equal(morphology.get(key), False, key)

    routing = scope.get("routing_scope", {})
    _require_equal(routing.get("live_integrated_mode"), MODE_NAME, "live integrated mode")
    _require_equal(routing.get("live_mode_nonzero_child_materialized"), False, "live materialization state")
    _require_equal(routing.get("energy_envelope_may_inform_future_autonomous_bound_design"), True, "bound-design evidence use")
    _require_equal(routing.get("energy_envelope_may_automatically_define_materialization_bound"), False, "automatic bound promotion")
    _require_equal(routing.get("exact_energy_root_may_automatically_define_materialized_child"), False, "root promotion")
    _require_equal(routing.get("larger_move_requires_explicit_joint_materialization_or_renormalization_contract"), True, "joint materialization contract")
    _require_equal(routing.get("reference_energy_gate_may_be_relaxed_for_larger_move"), False, "energy-gate relaxation")
    _require_equal(routing.get("live_next_task_remains_materialize_existing_compact_poloidal_child"), True, "live next task")
    _require_equal(routing.get("additional_basis_growth_allowed_before_materialization_screen"), False, "basis-growth guard")

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

    latest = project_status.get("latest_integrated_poloidal_capacity", {})
    _require_equal(latest.get("mode"), MODE_NAME, "project integrated poloidal mode")
    _require("nonzero_child_not_materialized" in latest.get("status", ""), "integrated mode unexpectedly materialized")
    _require_equal(latest.get("coefficient_bound_selected"), False, "integrated bound selection")
    _require_equal(latest.get("coefficient_value_selected"), False, "integrated value selection")
    _require_equal(latest.get("candidate_selection_resolved"), False, "integrated candidate selection")

    next_task = project_status.get("next_integration_task", "")
    for token in (
        MODE_NAME,
        "explicit autonomous coefficient bound/value",
        "new representation identity/candidate SHA",
        "fresh full-momentum/divergence validation",
    ):
        _require(token in next_task, f"live next integration task lost token: {token}")
    avoided = set(project_status.get("integration_policy", {}).get("avoid", ()))
    _require(
        "additional_capacity_basis_growth_before_a_current_integrated_direction_is_materialized_and_screened" in avoided,
        "live capacity-growth block missing",
    )
    _require(
        "coefficient_selection_from_capacity_or_morphology_alone" in avoided,
        "capacity-to-coefficient promotion block missing",
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
        "upstream_pr": 309,
        "source_classification": source,
        "canonical_thresholds_unchanged": True,
        "diagnostic_energy_envelope_selected_as_bound": False,
        "normalization_root_selected_as_coefficient": False,
        "joint_renormalization_requires_new_candidate_validation": True,
        "small_morphology_span_is_capacity_evidence_only": True,
        "live_next_integration_task_preserved": True,
        "velocity_export_ready": True,
        "visual_correspondence_verified": False,
        "pde_validated": False,
    }


def main() -> None:
    print(json.dumps(audit_compact_poloidal_energy_envelope_scope(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
