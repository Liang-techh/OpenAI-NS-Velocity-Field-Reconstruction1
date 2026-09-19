"""Fail-closed CR002 audit for PR #568 visual-fingerprint resolution semantics.

The governed experiment compares the exact PR #559 compensated ST052-M child on
25^3 and 33^3 fixed grids. Its fine-grid magnitude guards pass, while the
preregistered cross-resolution sign guard fails. This module preserves that
mixed result without laundering it into production selection, visual
correspondence, continuum convergence, or PDE validation.

It changes no velocity, pressure, forcing, threshold, candidate, or canonical
delivery state.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Mapping

TASK_ID = "CR002-ST052M-VISUAL-FINGERPRINT-RESOLUTION-SCOPE-089"
INTEGRATION_HEAD = "8ff26526e687a4c7178e7996018af09671a49b90"
PARENT_HEAD = "39b106ad8cb8df2064cabead3a12682089575e74"
PATH_HEAD = "196bcc5445039a2b3533128daa5c20cbb151cd58"
VISUAL_HEAD = "85f2f0940d80757313651d2ab50c6ef89dc22914"
CANONICAL_FAMILY = "eq45_supported_velocity_candidate_v1"
CANONICAL_SHA = "2fdff812c22131d56eac1b7d6e455187ff3207500c6385aa9abf282a8e3d7b1d"
CANONICAL_API = "openai_ns_reconstruction.eq45_supported_delivery:velocity"
ALLOWED_CLASSES = {"user_requirement", "public_source_fact", "autonomous_design", "pending_unknown"}
REGISTERED_STEPS = [0.02, 0.01, 0.005]
REGISTERED_QUADRATURE = [24, 48, 96]


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _load_json(path: Path) -> dict[str, Any]:
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise ValueError(f"expected JSON object in {path}")
    return data


def load_contract() -> dict[str, Any]:
    return _load_json(_repo_root() / "configs" / "st052m_visual_fingerprint_resolution_governance.json")


def load_constraints() -> dict[str, Any]:
    return _load_json(_repo_root() / "configs" / "constraints.json")


def load_delivery() -> dict[str, Any]:
    return _load_json(_repo_root() / "configs" / "velocity_delivery_contract.json")


def load_project_status() -> dict[str, Any]:
    return _load_json(_repo_root() / "project_status.json")


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _equal(actual: Any, expected: Any, label: str) -> None:
    _require(actual == expected, f"{label}: expected {expected!r}, got {actual!r}")


def _close(actual: float, expected: float, label: str, atol: float = 1e-12) -> None:
    _require(abs(float(actual) - float(expected)) <= atol, f"{label}: expected {expected!r}, got {actual!r}")


def audit_contract(
    contract: Mapping[str, Any],
    constraints: Mapping[str, Any],
    delivery: Mapping[str, Any],
    project_status: Mapping[str, Any],
) -> dict[str, Any]:
    _equal(contract["schema_version"], 1, "schema_version")
    _equal(contract["task_id"], TASK_ID, "task_id")

    snap = contract["snapshot"]
    _equal(snap["active_integration_branch"], "codex/cr001-constraints", "integration branch")
    _equal(snap["active_integration_head"], INTEGRATION_HEAD, "integration head")
    _equal(snap["governed_parent_pr"], 559, "governed parent PR")
    _equal(snap["governed_parent_head"], PARENT_HEAD, "governed parent head")
    _equal(snap["material_path_pr"], 565, "material path PR")
    _equal(snap["material_path_head"], PATH_HEAD, "material path head")
    _equal(snap["visual_fingerprint_pr"], 568, "visual fingerprint PR")
    _equal(snap["visual_fingerprint_head"], VISUAL_HEAD, "visual fingerprint head")
    _equal(snap["visual_fingerprint_standard_run"], 35423257712, "standard run")
    _equal(snap["visual_fingerprint_dedicated_run"], 35423257771, "dedicated run")
    _equal(snap["visual_fingerprint_replay_run"], 35423257732, "replay run")
    _equal(snap["visual_fingerprint_artifact_id"], 10578815385, "artifact id")
    _equal(
        snap["visual_fingerprint_artifact_digest"],
        "sha256:a82ea7d94e956f1213e2f1fb6b08e940e975eeeabdb9b3c2d249662986ea075c",
        "artifact digest",
    )

    _equal(set(contract["allowed_source_classes"]), ALLOWED_CLASSES, "allowed source classes")
    classes = contract["source_classification"]
    _equal(set(classes.values()), ALLOWED_CLASSES, "all four source classes represented")
    _equal(classes["callable_saveable_velocity_x_y_z_t"], "user_requirement", "delivery source class")
    _equal(classes["public_openai_velocity_visualization_observables"], "public_source_fact", "public visual source class")
    for key in (
        "st052m_kappa_0p05_taper_local_swirl_child",
        "grid_sizes_times_quantile_fd2_metrics_and_thresholds",
        "preregistered_cross_resolution_rule",
    ):
        _equal(classes[key], "autonomous_design", f"source class for {key}")
    for key in ("openai_hidden_numerical_velocity_or_morphology_thresholds", "openai_hidden_time_camera_seed_mapping"):
        _equal(classes[key], "pending_unknown", f"source class for {key}")

    protocol = contract["protocol"]
    _equal(protocol["child_identity"], "exact PR559 compensated child", "child identity")
    _equal(protocol["grid_sizes"], [25, 33], "grid sizes")
    _equal(protocol["times"], [0.25, 0.50, 0.75], "visual times")
    _equal(protocol["vorticity_operator"], "second_order_cartesian_curl", "vorticity operator")
    _close(protocol["top_vorticity_fraction"], 0.015, "top-vorticity fraction", atol=0.0)
    _require(protocol["openai_image_or_numeric_target_used"] is False, "target-free screen cannot gain an OpenAI numeric target")
    _require(protocol["visual_acceptance_threshold_against_openai_used"] is False, "target-free screen cannot become an OpenAI acceptance test")

    rule = contract["preregistered_rule"]
    _equal(rule["fine_grid"], 33, "fine grid")
    _equal(rule["coarse_grid"], 25, "coarse grid")
    _close(rule["fine_top_axial_gain_min"], 0.04, "fine axial threshold", atol=0.0)
    _close(rule["fine_top_radial_gain_max"], -0.015, "fine radial threshold", atol=0.0)
    _close(rule["fine_tip_radial_gain_max"], -0.007, "fine tip threshold", atol=0.0)
    _close(rule["fine_central_radial_abs_change_max"], 0.005, "fine central threshold", atol=0.0)
    _close(rule["fine_normalized_taper_axial_gain_retention_min"], 0.70, "retention threshold", atol=0.0)
    _require(rule["desired_sign_agreement_required_on_both_grids"] is True, "cross-resolution sign rule must remain preregistered")
    _require(rule["aggregate_clean_verdict"] is False, "scientific FAIL cannot be relabelled PASS")
    _equal(rule["failed_guard"], "cross_resolution_top_cloud_radial_sign_agreement", "failed guard")

    fine = contract["fine_grid_evidence"]
    _equal(fine["grid_size"], 33, "fine evidence grid")
    _equal(fine["top_axial_relative_changes"], [0.092626, 0.093798, 0.084810], "fine top axial changes")
    _equal(fine["top_radial_relative_changes"], [-0.030890, -0.035790, -0.055940], "fine top radial changes")
    _equal(fine["tip_radial_relative_changes"], [-0.012168, -0.012133, -0.012354], "fine tip changes")
    _equal(fine["central_radial_relative_changes"], [0.000270, 0.000134, 0.000087], "fine central changes")
    _equal(fine["normalized_taper_axial_gain_retention"], [1.0092, 1.0380, 1.0201], "fine axial retention")
    _require(fine["all_frozen_magnitude_guards_pass"] is True, "fine-grid magnitude result must remain recorded")

    coarse = contract["coarse_grid_evidence"]
    _equal(coarse["grid_size"], 25, "coarse evidence grid")
    _equal(coarse["top_axial_relative_changes"], [0.027411, 0.023478, 0.015450], "coarse top axial changes")
    _equal(coarse["top_radial_relative_changes"], [-0.008391, -0.003799, 0.004750], "coarse top radial changes")
    _equal(coarse["tip_radial_relative_changes"], [-0.015264, -0.015125, -0.015508], "coarse tip changes")
    _equal(coarse["central_radial_relative_changes"], [0.000977, 0.000461, 0.000237], "coarse central changes")
    _require(coarse["top_cloud_radial_sign_flips_at_t_0p75"] is True, "observed coarse-grid sign flip must remain explicit")
    _require(coarse["top_radial_relative_changes"][-1] > 0.0, "coarse t=.75 top-cloud radial metric must retain observed positive sign")
    _require(fine["top_radial_relative_changes"][-1] < 0.0, "fine t=.75 top-cloud radial metric must retain observed negative sign")

    resolution = contract["resolution_semantics"]
    for key in (
        "fine_grid_pass_can_override_preregistered_aggregate_fail",
        "two_grid_comparison_is_continuum_convergence_certificate",
        "two_grid_comparison_is_resolution_stability_certificate_for_top_cloud_radial_metric",
        "top_cloud_radial_tightening_has_same_desired_sign_on_both_grids_all_times",
        "stable_tip_band_sign_is_visual_correspondence_evidence",
        "stable_tip_band_sign_is_pde_evidence",
    ):
        _require(resolution[key] is False, f"forbidden resolution/evidence promotion became true: {key}")
    for key in (
        "thresholded_top_1p5pct_cloud_is_discretization_sensitive",
        "tip_band_radial_thinning_has_same_desired_sign_on_both_grids_all_times",
        "top_cloud_axial_gain_has_same_desired_sign_on_both_grids_all_times",
        "stable_tip_band_sign_is_descriptive_morphology_evidence",
        "changing_grid_quantile_metric_threshold_or_sign_rule_after_result_requires_new_experiment_identity",
    ):
        _require(resolution[key] is True, f"required resolution boundary lost: {key}")

    science = contract["scientific_scope"]
    _equal(science["scientific_verdict"], "FAIL", "scientific verdict")
    _require(science["failure_is_in_preregistered_resolution_guard"] is True, "failure scope must remain explicit")
    _require(science["fine_grid_morphology_signal_remains_descriptive"] is True, "fine-grid descriptive evidence must not be erased")
    for key in (
        "production_candidate_selected",
        "production_visual_fingerprint_selected",
        "visualization_ready",
        "visual_correspondence_verified",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    ):
        _require(science[key] is False, f"scientific state overpromoted: {key}")

    child = contract["child_delivery_scope"]
    for key in (
        "versioned_child_materialized",
        "child_candidate_id_assigned",
        "child_candidate_sha256_assigned",
        "save_load_contract_exists",
        "unified_velocity_api_registered",
        "velocity_export_ready",
        "visual_fingerprint_receipt_can_replace_canonical_delivery",
    ):
        _require(child[key] is False, f"unmaterialized child delivery overpromoted: {key}")

    cr = contract["canonical_cr001"]
    _close(cr["nu"], constraints["nu"], "nu", atol=0.0)
    _equal(cr["physical_domain"], constraints["domain"]["physical"], "physical domain")
    _equal(cr["evaluation_box"], constraints["domain"]["evaluation_box"], "evaluation box")
    _equal(cr["support"], constraints["domain"]["support"], "support")
    _equal(cr["time_interval"], constraints["domain"]["time_interval"], "time interval")
    _equal(cr["force_mode"], constraints["forcing"]["mode"], "force mode")
    _equal(cr["force_bounds"], constraints["forcing"]["parameters"], "force bounds")
    _close(cr["reference_energy"], constraints["nontriviality"]["reference_energy"], "reference energy", atol=0.0)
    _close(cr["reference_energy_abs_tolerance"], constraints["nontriviality"]["reference_energy_abs_tolerance"], "energy tolerance", atol=0.0)
    validation = constraints["validation"]
    _equal(cr["validation_seed"], validation["seed"], "validation seed")
    _equal(cr["held_out_points"], validation["held_out_points"], "held-out points")
    _equal(cr["validation_times"], validation["times"], "validation times")
    _equal(cr["derivative_steps"], REGISTERED_STEPS, "registered derivative steps")
    _equal(cr["derivative_steps"], validation["derivative_steps"], "constraint derivative steps")
    _equal(cr["quadrature_orders_per_axis"], REGISTERED_QUADRATURE, "registered energy ladder")
    _equal(cr["quadrature_orders_per_axis"], validation["quadrature_orders_per_axis"], "constraint energy ladder")
    _close(cr["divergence_max"], validation["thresholds"]["divergence_max"], "divergence max", atol=0.0)
    _close(cr["divergence_L2"], validation["thresholds"]["divergence_L2"], "divergence L2", atol=0.0)
    _close(cr["pde_residual_max"], validation["thresholds"]["pde_residual_max"], "momentum max", atol=0.0)
    _close(cr["pde_residual_L2"], validation["thresholds"]["pde_residual_L2"], "momentum L2", atol=0.0)
    _require(cr["free_or_residual_defined_force_allowed"] is False, "free residual force cannot be enabled")
    _require(cr["collapsed_velocity_success_allowed"] is False, "u->0 shortcut cannot be enabled")
    _require(cr["threshold_relaxation_under_same_experiment_allowed"] is False, "threshold relaxation cannot be enabled")

    canonical = contract["canonical_delivery"]
    primary = delivery["primary_deliverable"]
    _equal(canonical["candidate_family"], CANONICAL_FAMILY, "canonical family")
    _equal(canonical["candidate_sha256"], CANONICAL_SHA, "canonical SHA")
    _equal(canonical["velocity_api"], CANONICAL_API, "canonical API")
    _equal(primary["candidate_family"], CANONICAL_FAMILY, "delivery family")
    _equal(primary["candidate_sha256"], CANONICAL_SHA, "delivery SHA")
    _equal(primary["api"], CANONICAL_API, "delivery API")
    _equal(project_status["candidate_family"], CANONICAL_FAMILY, "project-status family")
    _equal(project_status["candidate_sha256"], CANONICAL_SHA, "project-status SHA")
    _equal(project_status["velocity_api"], CANONICAL_API, "project-status API")
    _require(canonical["velocity_export_ready"] is True, "canonical velocity must remain export-ready")
    _require(canonical["visual_fingerprint_child_replaces_canonical_candidate"] is False, "visual receipt cannot replace canonical candidate")
    _require(canonical["visual_fingerprint_child_replaces_canonical_velocity_api"] is False, "visual receipt cannot replace canonical API")

    states = contract["independent_truth_states"]
    _require(states["canonical_velocity_export_ready"] is True, "canonical velocity readiness must remain true")
    _require(states["visual_fingerprint_child_velocity_export_ready"] is False, "experimental child cannot inherit canonical readiness")
    for key in ("visualization_ready", "visual_correspondence_verified", "pde_validated", "paper_exact", "openai_field_identified", "blowup_proved"):
        _require(states[key] is False, f"independent truth state overpromoted: {key}")
    _require(delivery["claim_status"]["velocity_export_ready"] is True, "delivery contract lost canonical readiness")
    _require(delivery["claim_status"]["visual_correspondence_verified"] is False, "delivery contract visual correspondence unexpectedly promoted")
    _require(delivery["claim_status"]["pde_validated"] is False, "delivery contract PDE unexpectedly promoted")
    _require(project_status["states"]["velocity_export_ready"] is True, "project status lost canonical readiness")
    _require(project_status["states"]["visual_correspondence_verified"] is False, "project status visual correspondence unexpectedly promoted")
    _require(project_status["states"]["pde_validated"] is False, "project status PDE unexpectedly promoted")

    future = contract["future_research_rules"]
    for key in (
        "production_promotion_requires_versioned_child_identity_and_save_load",
        "visual_correspondence_requires_public_reference_comparison_and_resolution_stability",
        "continuum_claim_requires_more_than_two_fixed_grids_and_thresholded_cloud_sign_checks",
        "pde_promotion_requires_exact_child_pressure_restricted_force_and_full_registered_validation",
        "failed_cross_resolution_guard_may_not_be_relabelled_pass_by_reporting_only_33_grid",
    ):
        _require(future[key] is True, f"future governance rule lost: {key}")

    return {
        "status": "PASS",
        "task_id": TASK_ID,
        "governed_pr": 568,
        "scientific_verdict": "FAIL",
        "fine_grid_magnitude_guards_pass": True,
        "cross_resolution_sign_guard_pass": False,
        "stable_tip_band_sign": True,
        "top_cloud_radial_resolution_stable": False,
        "continuum_convergence_certified": False,
        "visual_correspondence_verified": False,
        "child_velocity_export_ready": False,
        "canonical_velocity_export_ready": True,
        "pde_validated": False,
        "paper_exact": False,
    }


def audit_repository() -> dict[str, Any]:
    return audit_contract(load_contract(), load_constraints(), load_delivery(), load_project_status())
