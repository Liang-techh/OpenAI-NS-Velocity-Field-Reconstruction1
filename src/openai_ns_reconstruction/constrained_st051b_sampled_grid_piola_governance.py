"""Fail-closed CR002 governance for sampled-grid Piola semantics.

The exact Piola divergence identity is a continuum statement about a sufficiently
smooth callable parent.  PR #499 instead applies the same spatial map through a
trilinear interpolant of a frozen 33^3 sampled grid.  That object is useful
render/capacity evidence, but it is a distinct numerical representation and may
not inherit continuum divergence, CR001 energy, PDE, visual-correspondence, or
candidate-identity claims.

This module changes no velocity, candidate, pressure, forcing, route, threshold,
or upstream implementation.
"""
from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
from typing import Any, Mapping


TASK_ID = "CR002-SAMPLED-GRID-PIOLA-SEMANTICS-083"
_INTEGRATION_HEAD = "fec97a1f25e914069e528dac0852d9a46244f159"
_PR499_HEAD = "cf5c678117678abe39704f9ee1c66180a0409916"
_PR499_RUN = 35400911257
_PR499_ARTIFACT = 10570795164
_PR499_ARTIFACT_DIGEST = (
    "sha256:ec99470af7028e0004877def1065d11e4ee66c7b234c95283b20e923f2235f6d"
)
_CANONICAL_FAMILY = "eq45_supported_velocity_candidate_v1"
_CANONICAL_SHA = "2fdff812c22131d56eac1b7d6e455187ff3207500c6385aa9abf282a8e3d7b1d"
_CANONICAL_API = "openai_ns_reconstruction.eq45_supported_delivery:velocity"
_ACTIVE_ROUTE = "materialize_integrated_axial_cap_poloidal_child_then_fresh_validate"
_ALLOWED_CLASSES = {
    "user_requirement",
    "public_source_fact",
    "autonomous_design",
    "pending_unknown",
}


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _load_json(path: Path) -> dict[str, Any]:
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise ValueError(f"expected object in {path}")
    return data


def load_contract() -> dict[str, Any]:
    return _load_json(_repo_root() / "configs" / "st051b_sampled_grid_piola_governance.json")


def load_constraints() -> dict[str, Any]:
    return _load_json(_repo_root() / "configs" / "constraints.json")


def load_project_status() -> dict[str, Any]:
    return _load_json(_repo_root() / "project_status.json")


def load_axial_warp_governance() -> dict[str, Any]:
    return _load_json(_repo_root() / "configs" / "axial_coordinate_warp_governance.json")


def load_grid_identity_contract() -> dict[str, Any]:
    return _load_json(_repo_root() / "configs" / "velocity_grid_delivery_identity_contract.json")


def load_render_scope_contract() -> dict[str, Any]:
    return _load_json(_repo_root() / "configs" / "render_readiness_scope_contract.json")


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _equal(actual: Any, expected: Any, label: str) -> None:
    _require(actual == expected, f"{label} drifted: {actual!r} != {expected!r}")


def _false(value: Any, label: str) -> None:
    _equal(value, False, label)


def _true(value: Any, label: str) -> None:
    _equal(value, True, label)


def _audit_cr001(contract: Mapping[str, Any], constraints: Mapping[str, Any]) -> None:
    frozen = contract["canonical_cr001"]
    domain = constraints["domain"]
    forcing = constraints["forcing"]
    nontrivial = constraints["nontriviality"]
    validation = constraints["validation"]
    thresholds = validation["thresholds"]

    _equal(frozen["nu"], constraints["nu"], "nu")
    _equal(frozen["physical_domain"], domain["physical"], "physical domain")
    _equal(frozen["evaluation_box"], domain["evaluation_box"], "evaluation box")
    _equal(frozen["support"], domain["support"], "support")
    _equal(frozen["time_interval"], domain["time_interval"], "time interval")
    _equal(frozen["force_mode"], forcing["mode"], "force mode")
    _equal(frozen["force_bounds"], forcing["parameters"], "force bounds")
    _equal(frozen["reference_energy"], nontrivial["reference_energy"], "reference energy")
    _equal(
        frozen["reference_energy_abs_tolerance"],
        nontrivial["reference_energy_abs_tolerance"],
        "reference energy tolerance",
    )
    _equal(frozen["validation_seed"], validation["seed"], "validation seed")
    _equal(frozen["held_out_points"], validation["held_out_points"], "held-out points")
    _equal(frozen["validation_times"], validation["times"], "validation times")
    _equal(frozen["derivative_steps"], validation["derivative_steps"], "derivative steps")
    _equal(frozen["divergence_max"], thresholds["divergence_max"], "divergence max")
    _equal(frozen["divergence_L2"], thresholds["divergence_L2"], "divergence L2")
    _equal(frozen["pde_residual_max"], thresholds["pde_residual_max"], "momentum max")
    _equal(frozen["pde_residual_L2"], thresholds["pde_residual_L2"], "momentum L2")
    _require("No residual-dependent basis or pointwise free force" in forcing["restriction"], "free-force guard drifted")


def audit_sampled_grid_piola_governance(
    contract: Mapping[str, Any] | None = None,
    constraints: Mapping[str, Any] | None = None,
    project_status: Mapping[str, Any] | None = None,
    axial_warp_governance: Mapping[str, Any] | None = None,
    grid_identity_contract: Mapping[str, Any] | None = None,
    render_scope_contract: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Audit the sampled/interpolated Piola evidence boundary against live contracts."""

    contract = deepcopy(dict(load_contract() if contract is None else contract))
    constraints = deepcopy(dict(load_constraints() if constraints is None else constraints))
    project_status = deepcopy(dict(load_project_status() if project_status is None else project_status))
    axial = deepcopy(
        dict(load_axial_warp_governance() if axial_warp_governance is None else axial_warp_governance)
    )
    grid = deepcopy(
        dict(load_grid_identity_contract() if grid_identity_contract is None else grid_identity_contract)
    )
    render = deepcopy(
        dict(load_render_scope_contract() if render_scope_contract is None else render_scope_contract)
    )

    _equal(contract.get("schema_version"), 1, "contract schema")
    _equal(contract.get("task_id"), TASK_ID, "task id")

    snapshot = contract["snapshot"]
    _equal(snapshot["active_integration_branch"], "codex/cr001-constraints", "integration branch")
    _equal(snapshot["active_integration_head"], _INTEGRATION_HEAD, "integration head")
    _equal(snapshot["observed_pr"], 499, "observed PR")
    _equal(snapshot["observed_pr_head"], _PR499_HEAD, "observed PR head")
    _equal(snapshot["observed_pr_state"], "open_unintegrated_capacity_evidence", "PR state")
    _equal(snapshot["dedicated_workflow_run"], _PR499_RUN, "dedicated workflow run")
    _equal(snapshot["dedicated_workflow_conclusion"], "success", "dedicated workflow conclusion")
    _equal(snapshot["receipt_artifact_id"], _PR499_ARTIFACT, "receipt artifact")
    _equal(snapshot["receipt_artifact_digest"], _PR499_ARTIFACT_DIGEST, "receipt digest")

    classes = set(contract["allowed_source_classes"])
    _equal(classes, _ALLOWED_CLASSES, "allowed source classes")
    for key, value in contract["source_classification"].items():
        _require(value in classes, f"unknown source class for {key}: {value}")
    _equal(contract["source_classification"]["callable_saveable_velocity_x_y_z_t"], "user_requirement", "velocity source class")
    _equal(contract["source_classification"]["public_openai_visualization_observables"], "public_source_fact", "public visualization class")
    _equal(contract["source_classification"]["static_axial_piola_map_beta_0p075"], "autonomous_design", "Piola source class")
    _equal(contract["source_classification"]["openai_hidden_coordinate_map"], "pending_unknown", "hidden map class")

    layers = contract["representation_layers"]
    continuum = layers["continuum_piola_operator"]
    sampled_parent = layers["sampled_parent"]
    sampled_child = layers["sampled_piola_child"]
    coarse = layers["nested_coarse_view"]

    _equal(continuum["beta"], 0.075, "screen beta")
    _false(continuum["beta_time_dependent"], "beta time dependence")
    _equal(continuum["identity_scope"], "sufficiently_smooth_callable_parent", "continuum identity scope")
    _true(continuum["identity_is_analytic_consequence_of_autonomous_map"], "continuum identity status")
    _equal(sampled_parent["fine_resolution"], 33, "fine resolution")
    _equal(sampled_parent["times"], [0.25, 0.5, 0.75], "sample times")
    _equal(sampled_parent["interpolation"], "trilinear", "interpolation")
    _false(sampled_parent["is_underlying_callable_parent"], "sampled grid callable-parent identity")
    _false(sampled_parent["is_cr001_held_out_sample"], "sampled grid acceptance-sample identity")
    _true(sampled_child["new_numerical_representation"], "sampled Piola representation identity")
    _false(sampled_child["same_object_as_P_beta_of_underlying_callable_parent"], "sampled/continuum object equivalence")
    _false(sampled_child["beta_alone_identifies_child"], "beta-only child identity")
    _false(sampled_child["continuum_divergence_identity_is_discrete_divergence_certificate"], "continuum-to-discrete divergence transfer")
    _false(sampled_child["grid_level_divergence_preflight_if_run_is_pde_acceptance"], "grid divergence/PDE transfer")
    required_identity_inputs = {
        "source_grid_identity",
        "grid_resolution_and_nodes",
        "sample_times",
        "interpolation_kernel_and_boundary_behavior",
        "beta",
        "sampled_energy_normalization_rule",
    }
    _equal(set(sampled_child["identity_depends_on"]), required_identity_inputs, "sampled child identity inputs")
    _equal(coarse["resolution"], 17, "coarse resolution")
    _false(coarse["independent_resample_from_callable_parent"], "coarse independent resampling")
    _false(coarse["continuum_convergence_study"], "coarse continuum convergence")

    # Bind the new seam to the already-integrated generic governance instead of
    # silently weakening either contract.
    _equal(
        axial["representation"]["divergence_identity"],
        continuum["divergence_identity"],
        "continuum Piola identity",
    )
    _true(axial["representation"]["navier_stokes_invariant_under_warp"] is False, "warp non-invariance")
    _false(grid["identity_layers"]["sampled_grid"]["sampled_grid_is_continuum_velocity_api"], "grid/continuum API distinction")
    _true(
        grid["identity_layers"]["sampled_grid"]["off_grid_or_unsampled_time_interpolation_requires_explicit_new_numerical_contract"],
        "interpolation contract requirement",
    )
    _false(grid["grid_evidence_scope"]["grid_derivatives_are_cr001_pde_acceptance_evidence"], "grid derivative PDE evidence")
    _equal(render["integrated_pr490"]["interpolation"], "trilinear_on_sampled_grid", "render interpolation")
    _false(render["integrated_pr490"]["sampled_vorticity_may_be_claimed_as_pde_acceptance_evidence"], "render vorticity PDE evidence")

    evidence = contract["observed_pr499_evidence"]
    _require(evidence["source_roundtrip_max_abs"] <= 5e-13, "sampled source roundtrip guard failed")
    _true(evidence["source_roundtrip_proves_only_frozen_sampled_transform_replay"], "roundtrip scope")
    _false(evidence["source_roundtrip_recovers_underlying_analytic_parent_off_grid"], "roundtrip off-grid recovery")
    _equal(evidence["fine_common_reference_energy_scale"], 1.0035328035496631, "fine sampled energy scale")
    _equal(evidence["coarse_common_reference_energy_scale"], 1.006988221666327, "coarse sampled energy scale")
    _require(
        abs(evidence["fine_common_reference_energy_scale"] - evidence["coarse_common_reference_energy_scale"]) > 1e-3,
        "resolution-specific energy-scale distinction disappeared",
    )
    _true(evidence["sampled_energy_scale_is_resolution_specific"], "sampled energy scale scope")
    _false(evidence["sampled_energy_restoration_is_cr001_reference_energy_validation"], "sampled energy/CR001 energy transfer")
    _true(evidence["useful_axial_capacity_transfer"], "capacity result")
    _equal(evidence["response_rank"], 2, "response rank")
    _false(evidence["capacity_pass_selects_production_beta"], "capacity production-beta selection")
    _false(evidence["capacity_pass_selects_production_candidate"], "capacity production-candidate selection")
    _false(evidence["capacity_pass_is_visual_correspondence"], "capacity visual correspondence")
    _false(evidence["capacity_pass_is_pde_validation"], "capacity PDE validation")

    for inference in contract["forbidden_inferences"]:
        _false(inference["allowed"], f"forbidden inference {inference['from']} -> {inference['to']}")

    future = contract["future_materialization_requirements"]
    for key in (
        "explicit_production_beta_bound_and_value_required",
        "new_representation_identity_required",
        "new_candidate_sha_required",
        "candidate_save_load_replay_required",
        "fresh_reference_and_validation_energy_required",
        "fresh_axis_near_axis_support_and_boundary_checks_required",
        "fresh_independent_cartesian_divergence_required_before_divergence_claim",
        "rebuild_or_refit_compatible_pressure_and_restricted_force_required_before_PDE_claim",
        "fresh_full_momentum_and_divergence_on_unused_CR001_acceptance_data_required_before_PDE_promotion",
        "public_reference_comparison_and_independent_resolution_stability_required_before_visual_correspondence",
    ):
        _true(future[key], f"future requirement {key}")

    _audit_cr001(contract, constraints)

    delivery = contract["canonical_delivery"]
    _equal(delivery["candidate_family"], _CANONICAL_FAMILY, "canonical family")
    _equal(delivery["candidate_sha256"], _CANONICAL_SHA, "canonical SHA")
    _equal(delivery["velocity_api"], _CANONICAL_API, "canonical velocity API")
    _true(delivery["velocity_export_ready"], "canonical velocity export")
    _equal(delivery["active_scientific_route"], _ACTIVE_ROUTE, "active scientific route")
    _false(delivery["pr499_changes_canonical_candidate"], "PR499 canonical candidate mutation")
    _false(delivery["pr499_changes_canonical_velocity_api"], "PR499 canonical API mutation")
    _false(delivery["pr499_changes_active_scientific_route"], "PR499 route mutation")

    _equal(project_status["candidate_family"], _CANONICAL_FAMILY, "project-status family")
    _equal(project_status["candidate_sha256"], _CANONICAL_SHA, "project-status SHA")
    _equal(project_status["velocity_api"], _CANONICAL_API, "project-status API")
    _equal(project_status["active_scientific_route"], _ACTIVE_ROUTE, "project-status route")
    _true(project_status["states"]["velocity_export_ready"], "project-status velocity export")
    for state in ("visualization_ready", "visual_correspondence_verified", "pde_validated", "paper_exact", "openai_field_identified", "blowup_proved"):
        _false(project_status["states"][state], f"project-status {state}")

    truth = contract["truth_boundary"]
    _true(truth["canonical_velocity_export_ready"], "truth velocity export")
    for key in (
        "pr499_sampled_piola_child_materialized_as_callable_candidate",
        "production_beta_selected",
        "production_candidate_selected",
        "visualization_ready",
        "visual_correspondence_verified",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    ):
        _false(truth[key], f"truth boundary {key}")

    for key, value in contract["mutation_scope"].items():
        _false(value, f"mutation scope {key}")

    return {
        "task_id": TASK_ID,
        "audit": "passed",
        "observed_pr": 499,
        "observed_pr_head": _PR499_HEAD,
        "sampled_piola_is_distinct_numerical_representation": True,
        "continuum_divergence_transfer_allowed": False,
        "sampled_energy_is_cr001_energy_acceptance": False,
        "nested_17_is_continuum_convergence": False,
        "canonical_velocity_export_ready": True,
        "pde_validated": False,
    }


def main() -> int:
    print(json.dumps(audit_sampled_grid_piola_governance(), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
