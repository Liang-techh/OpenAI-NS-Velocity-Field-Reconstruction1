"""Fail-closed CR002 governance for the compact bipolar poloidal capacity mode.

This module changes no velocity value.  It audits the inference boundary around
``COMPACT_C4_ODD_Z_POLOIDAL`` after that capacity experiment was integrated.
"""
from __future__ import annotations

from copy import deepcopy
import inspect
import json
from pathlib import Path
from typing import Any, Mapping

from . import constrained_eq45_bipolar_compact_poloidal_capacity as upstream

TASK_ID = "CR002-BIPOLAR-COMPACT-POLOIDAL-REPRESENTATION-SCOPE-036"
_ALLOWED_CLASSES = {
    "user_requirement",
    "public_source_fact",
    "autonomous_design",
    "pending_unknown",
}
_FALSE_TRUTH_STATES = (
    "compact_poloidal_nonzero_child_materialized",
    "compact_poloidal_coefficient_selected",
    "candidate_selection_resolved",
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
        raise ValueError(f"expected object in {path}")
    return data


def load_scope() -> dict[str, Any]:
    return _load_json(_repo_root() / "configs" / "bipolar_compact_poloidal_representation_scope.json")


def load_constraints() -> dict[str, Any]:
    return _load_json(_repo_root() / "configs" / "constraints.json")


def load_project_status() -> dict[str, Any]:
    return _load_json(_repo_root() / "project_status.json")


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _require_equal(actual: Any, expected: Any, label: str) -> None:
    _require(actual == expected, f"{label} drifted: {actual!r} != {expected!r}")


def audit_compact_poloidal_scope(
    scope: Mapping[str, Any] | None = None,
    constraints: Mapping[str, Any] | None = None,
    project_status: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Audit the compact-poloidal representation and claim boundary.

    Optional mappings exist so regressions can mutate one semantic at a time
    without modifying repository files.
    """

    scope = deepcopy(dict(load_scope() if scope is None else scope))
    constraints = deepcopy(dict(load_constraints() if constraints is None else constraints))
    project_status = deepcopy(dict(load_project_status() if project_status is None else project_status))

    _require_equal(scope.get("task_id"), TASK_ID, "scope task id")
    _require_equal(
        scope.get("governed_upstream_task"),
        "CR003-BIPOLAR-COMPACT-POLOIDAL-AXIAL-CAPACITY-041",
        "governed upstream task",
    )
    _require_equal(scope.get("mode_name"), "COMPACT_C4_ODD_Z_POLOIDAL", "mode name")
    _require_equal(
        scope.get("representation_kind"),
        "autonomous_physical_space_additive_curl_correction",
        "representation kind",
    )

    allowed = set(scope.get("allowed_source_classes", ()))
    _require_equal(allowed, _ALLOWED_CLASSES, "source-class vocabulary")
    source = scope.get("source_classification", {})
    _require(source and set(source.values()) <= _ALLOWED_CLASSES, "invalid source classification")
    _require_equal(source.get("callable_velocity_delivery"), "user_requirement", "velocity delivery class")
    _require_equal(source.get("eq45_backbone"), "public_source_fact", "Eq45 backbone class")
    _require_equal(source.get("compact_poloidal_shape"), "autonomous_design", "compact shape class")
    _require_equal(source.get("plateau_geometry_reuse"), "autonomous_design", "plateau reuse class")
    _require_equal(source.get("openai_hidden_numerical_profile"), "pending_unknown", "hidden profile class")

    capacity = scope.get("capacity_evidence_only", {})
    _require_equal(capacity.get("public_velocity_rank"), 3, "capacity rank")
    _require(capacity.get("normalized_condition_number", 0.0) > 1.0, "condition number missing")
    _require(0.0 < capacity.get("novelty_outside_phi01_phi03_span", 0.0) < 1.0, "novelty out of range")
    for key in (
        "diagnostic_coefficients_select_value",
        "diagnostic_coefficients_define_materialization_bound",
        "inherited_profile_coefficient_limit_is_materialization_bound",
        "zero_audited_flank_response_is_support_validation",
        "analytic_or_fd_divergence_response_is_registered_divergence_acceptance",
        "morphology_leverage_is_visual_correspondence",
        "capacity_rank_is_pde_leverage",
    ):
        _require_equal(capacity.get(key), False, key)
    _require_equal(capacity.get("diagnostic_coefficients"), [-0.25, 0.25], "diagnostic coefficients")

    semantics = scope.get("shape_semantics", {})
    for key in (
        "axis_regular_by_construction",
        "radial_component_even_in_z",
        "axial_component_odd_in_z",
        "swirl_component_zero_for_basis_response",
        "parity_compatible_with_current_bipolar_class",
        "shape_semantics_are_not_public_source_facts",
    ):
        _require_equal(semantics.get(key), True, key)
    _require_equal(semantics.get("positive_coefficient_midplane_radial_direction"), "inward", "radial sign")
    _require_equal(
        semantics.get("positive_coefficient_axial_parity"),
        "opposite_sign_above_below_midplane",
        "axial parity",
    )

    morphology = scope.get("morphology_scope", {})
    _require_equal(morphology.get("bulk_axial_rms_leverage_observed"), True, "bulk axial leverage")
    _require_equal(morphology.get("radial_rms_leverage_also_nonzero"), True, "radial collateral")
    _require_equal(morphology.get("diagnostic_q90_q99_axial_reach_moved"), False, "q90/q99 reach")
    _require_equal(morphology.get("tip_geometry_solution_claimed"), False, "tip-geometry claim")

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
        (
            snap.get("reference_energy_abs_tolerance"),
            nontriviality.get("reference_energy_abs_tolerance"),
            "reference energy tolerance",
        ),
        (
            snap.get("validation_energy_range"),
            [
                nontriviality.get("minimum_energy_each_validation_time"),
                nontriviality.get("maximum_energy_each_validation_time"),
            ],
            "validation energy range",
        ),
        (snap.get("validation_seed"), validation.get("seed"), "validation seed"),
        (snap.get("held_out_points"), validation.get("held_out_points"), "held-out count"),
        (snap.get("derivative_steps"), validation.get("derivative_steps"), "derivative ladder"),
        (snap.get("divergence_max"), thresholds.get("divergence_max"), "divergence max"),
        (snap.get("divergence_L2"), thresholds.get("divergence_L2"), "divergence L2"),
        (snap.get("pde_residual_max"), thresholds.get("pde_residual_max"), "PDE max"),
        (snap.get("pde_residual_L2"), thresholds.get("pde_residual_L2"), "PDE L2"),
    )
    for expected, actual, label in expected_pairs:
        _require_equal(actual, expected, label)
    _require("No residual-dependent basis" in forcing.get("restriction", ""), "free-force restriction weakened")

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

    _require_equal(upstream.TASK_ID, "CR003-BIPOLAR-COMPACT-POLOIDAL-AXIAL-CAPACITY-041", "upstream task id")
    _require_equal(upstream.COMPACT_POLOIDAL, "COMPACT_C4_ODD_Z_POLOIDAL", "upstream mode id")
    for key in (
        "velocity_changed",
        "candidate_artifact_changed",
        "new_basis_materialized",
        "compact_poloidal_coefficient_selected",
        "perturbed_candidate_pde_residual_evaluated",
        "force_or_pressure_changed",
        "visualization_ready",
        "visual_correspondence_verified",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
    ):
        _require_equal(upstream.TRUTH_BOUNDARY.get(key), False, f"upstream truth {key}")
    source_text = inspect.getsource(upstream._compact_velocity)
    _require(
        "implementation guard for the diagnostic, not a selected bound" in source_text,
        "upstream coefficient-limit semantics drifted",
    )

    return {
        "task_id": TASK_ID,
        "status": "governance_pass",
        "mode_name": scope["mode_name"],
        "source_classification": source,
        "canonical_thresholds_unchanged": True,
        "nonzero_child_materialized": False,
        "coefficient_selected": False,
        "pde_validated": False,
        "visual_correspondence_verified": False,
    }


def main() -> None:
    print(json.dumps(audit_compact_poloidal_scope(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
