"""Fail-closed CR002 audit for the ST052-M local-swirl energy-root child.

PR #559 closes one autonomous order-64 reference-energy equation for a fixed
localized taper plus a compact shoulder swirl multiplier. PR #565 then replays
that exact representation on the frozen 48-path protocol. This module keeps the
energy root, zero-change path receipt, delivery state, visual state, and PDE
state separate and candidate-specific.

It changes no velocity, pressure, forcing, threshold, candidate, or canonical
delivery state.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Mapping

TASK_ID = "CR002-ST052M-LOCAL-SWIRL-ENERGY-ROOT-SCOPE-088"
INTEGRATION_HEAD = "253ad1f9042e4eb5eebe6cb36e94c467a6a03d0e"
PARENT_HEAD = "779ffca71066e2864496d37de55a7aafc45d6f57"
TAPER_HEAD = "093c7171cd61c6bd439afa30b2da69598a02d182"
FAILED_COMP_HEAD = "62e8c170427d5d830d7f897ba31768e0fc4ce56a"
ROOT_HEAD = "39b106ad8cb8df2064cabead3a12682089575e74"
PATH_HEAD = "196bcc5445039a2b3533128daa5c20cbb151cd58"
CANONICAL_FAMILY = "eq45_supported_velocity_candidate_v1"
CANONICAL_SHA = "2fdff812c22131d56eac1b7d6e455187ff3207500c6385aa9abf282a8e3d7b1d"
CANONICAL_API = "openai_ns_reconstruction.eq45_supported_delivery:velocity"
ALLOWED_CLASSES = {"user_requirement", "public_source_fact", "autonomous_design", "pending_unknown"}
REGISTERED_TIMES = [0.25, 0.3125, 0.4375, 0.5625, 0.6875, 0.75]
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
    return _load_json(_repo_root() / "configs" / "st052m_local_swirl_energy_root_governance.json")


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
    _equal(snap["redistributed_parent_pr"], 528, "redistributed parent PR")
    _equal(snap["redistributed_parent_head"], PARENT_HEAD, "redistributed parent head")
    _equal(snap["localized_taper_pr"], 544, "localized taper PR")
    _equal(snap["localized_taper_head"], TAPER_HEAD, "localized taper head")
    _equal(snap["failed_piola_compensation_pr"], 549, "failed compensation PR")
    _equal(snap["failed_piola_compensation_head"], FAILED_COMP_HEAD, "failed compensation head")
    _equal(snap["local_swirl_energy_pr"], 559, "local swirl root PR")
    _equal(snap["local_swirl_energy_head"], ROOT_HEAD, "local swirl root head")
    _equal(snap["local_swirl_energy_standard_run"], 35420587077, "root standard run")
    _equal(snap["local_swirl_energy_dedicated_run"], 35420587120, "root dedicated run")
    _equal(snap["local_swirl_energy_replay_run"], 35420587103, "root replay run")
    _equal(snap["material_path_pr"], 565, "material-path PR")
    _equal(snap["material_path_head"], PATH_HEAD, "material-path head")
    _equal(snap["material_path_standard_run"], 35422700284, "path standard run")
    _equal(snap["material_path_dedicated_run"], 35422700348, "path dedicated run")
    _equal(snap["material_path_publication_run"], 35422700425, "path publication run")
    _equal(snap["material_path_artifact_id"], 10578310175, "path artifact id")
    _equal(snap["material_path_artifact_digest"], "sha256:740eb9f60e1dbfe3aa58cff79ed46c8572a4ff430554f7e6da041831423f91d2", "path artifact digest")
    _equal(snap["material_path_report_sha256"], "9e82a2743cb28da9fdcb94ff2a38f8d4e3420bd071bd8182f8dbda75f8f2d736", "path report sha")

    _equal(set(contract["allowed_source_classes"]), ALLOWED_CLASSES, "allowed source classes")
    classes = contract["source_classification"]
    _equal(set(classes.values()), ALLOWED_CLASSES, "all source classes represented")
    for key in (
        "st052m_and_frozen_kappa_0p05_redistribution",
        "localized_piola_tau_0p05_and_tip_window",
        "local_swirl_shoulder_window_and_energy_root_rule",
        "beta_0p08837490297155456",
        "order_64_energy_root_screen",
        "frozen_48_path_protocol",
    ):
        _equal(classes[key], "autonomous_design", f"source class for {key}")
    for key in (
        "openai_hidden_numerical_velocity_or_compensation_map",
        "openai_hidden_parameters_time_camera_seed_mapping",
    ):
        _equal(classes[key], "pending_unknown", f"source class for {key}")

    root = contract["root_identity"]
    _equal(root["task_id"], "CR003-ST052M-LOCAL-SWIRL-ENERGY-082", "root task id")
    _equal(root["parent_representation"], "ST052-M+frozen-redistribution-.05", "root parent")
    _close(root["redistribution_alpha"], 2.520520814687742, "redistribution alpha", atol=2e-15)
    _close(root["redistribution_gain"], 0.05, "redistribution gain", atol=0.0)
    _equal(root["redistribution_inner_window"], [0.30, 1.05], "inner window")
    _equal(root["redistribution_outer_window"], [0.95, 1.85], "outer window")
    _close(root["taper_tau"], 0.05, "taper tau", atol=0.0)
    _equal(root["tip_window_abs_z_over_2"], [0.50, 0.82], "tip window")
    _equal(root["shoulder_window_abs_z_over_2"], [0.36, 0.49], "shoulder window")
    _equal(root["beta_bracket"], [0.0, 0.9], "beta bracket")
    _close(root["beta"], 0.08837490297155456, "energy beta", atol=2e-15)
    _equal(root["energy_quadrature_order"], 64, "energy quadrature")
    _close(root["post_transform_common_scale"], 1.0, "post transform scale", atol=0.0)
    _require(root["parameter_scan_performed"] is False, "energy root cannot become a parameter scan")
    _require(root["energy_only_dependent_root_solve"] is True, "beta must remain energy-only dependent")
    _equal(root["independently_tuned_new_coefficients"], 0, "independently tuned coefficients")
    _require(float(root["energy_relative_error_max_reported"]) <= 2.0e-12, "reported energy root closure drifted")
    _require(float(root["minimum_swirl_factor_reported"]) > 0.0, "screened swirl multiplier must remain positive")
    _require(root["root_exists_in_frozen_bracket"] is True, "successful frozen energy root must remain recorded")
    _require(root["clean_local_swirl_energy_capacity"] is True, "capacity result must remain recorded")

    semantics = contract["root_semantics"]
    _require(semantics["beta_is_conditioned_on_exact_parent_taper_window_quadrature_and_energy_target"] is True, "beta conditioning must remain explicit")
    for key in (
        "beta_is_portable_across_changed_parent_or_transform_identity",
        "beta_is_a_production_coefficient",
        "beta_is_an_openai_hidden_parameter",
        "energy_root_is_a_visual_fit",
        "energy_root_is_a_pde_fit",
        "reusing_beta_without_rechecking_energy_on_changed_representation_allowed",
    ):
        _require(semantics[key] is False, f"forbidden beta inference became true: {key}")
    _require(semantics["changing_parent_taper_window_quadrature_or_energy_rule_requires_new_experiment_identity"] is True, "changed root definition requires new experiment identity")

    rep = contract["representation_scope"]
    for key in (
        "axisymmetric_theta_only_multiplier",
        "shoulder_multiplier_is_exactly_inactive_in_registered_central_band",
        "post_transform_common_scale_is_unity",
        "central_velocity_pointwise_identity_under_registered_central_probe",
        "minimum_reported_multiplier_is_positive",
        "axis_swirl_regularity_order_preserved_if_parent_is_axis_regular",
        "continuous_axisymmetric_divergence_structure_preserved_by_theta_only_z_multiplier",
    ):
        _require(rep[key] is True, f"representation property lost: {key}")
    for key in (
        "central_identity_implies_entire_velocity_field_identity",
        "swirl_sign_reversal_introduced_by_screened_multiplier",
        "theta_only_divergence_structure_is_navier_stokes_invariance",
        "parent_pressure_receipt_may_transfer",
        "parent_restricted_force_receipt_may_transfer",
        "parent_momentum_receipt_may_transfer",
        "custom_cartesian_fd_divergence_is_formal_cr001_acceptance",
    ):
        _require(rep[key] is False, f"forbidden representation/PDE inference became true: {key}")

    euler = contract["eulerian_capacity_evidence"]
    _require(euler["tip_radial_vorticity_rms_reduction_observed"] is True, "tip morphology evidence must remain recorded")
    _require(euler["central_radial_vorticity_change_registered_as_zero"] is True, "central morphology identity must remain recorded")
    _require(euler["global_axial_vorticity_rms_gain_observed"] is True, "axial morphology evidence must remain recorded")
    _equal(euler["response_rank"], 4, "response rank")
    _require(float(euler["response_condition_number"]) < 2.0, "response condition unexpectedly degraded in governance receipt")
    _require(euler["target_free_morphology_only"] is True, "morphology evidence must remain target-free")
    for key in ("public_openai_numeric_target_used", "can_promote_visual_correspondence_verified", "can_promote_pde_validated"):
        _require(euler[key] is False, f"Eulerian evidence overpromoted: {key}")

    paths = contract["material_path_evidence"]
    _require(paths["exact_child_replayed"] is True, "exact-child path replay must remain recorded")
    _equal(paths["frozen_path_count"], 48, "path count")
    _equal(paths["frozen_axial_pair_count"], 24, "axial pair count")
    _equal(paths["time_interval"], [0.25, 0.75], "path time interval")
    _equal(paths["seed_radii"], [0.6, 0.9, 1.2], "seed radii")
    _equal(paths["seed_z"], [-0.3, 0.3], "seed z")
    _equal(paths["azimuth_count_per_radius_z"], 8, "azimuth count")
    _close(paths["rtol"], 1e-9, "path rtol", atol=0.0)
    _close(paths["atol"], 1e-11, "path atol", atol=0.0)
    _close(paths["max_step"], 0.01, "path max step", atol=0.0)
    for left, right, label in (
        ("mean_absolute_turns_control", "mean_absolute_turns_child", "mean turns"),
        ("maximum_absolute_turns_control", "maximum_absolute_turns_child", "max turns"),
        ("contraction_magnitude_control", "contraction_magnitude_child", "contraction"),
        ("mean_pair_axial_separation_control", "mean_pair_axial_separation_child", "pair separation"),
    ):
        _close(paths[left], paths[right], f"reported zero-change {label}", atol=0.0)
    _require(paths["all_reported_aggregate_relative_changes_round_to_zero_percent"] is True, "zero-change path receipt must remain recorded")
    _equal(paths["inward_paths_control"], 48, "control inward paths")
    _equal(paths["inward_paths_child"], 48, "child inward paths")
    for key in (
        "zero_change_on_frozen_paths_proves_global_trajectory_identity",
        "zero_change_on_frozen_paths_proves_velocity_field_identity",
        "zero_change_on_frozen_paths_negates_eulerian_morphology_change",
        "material_path_evidence_is_pde_evidence",
    ):
        _require(paths[key] is False, f"frozen path evidence overgeneralized: {key}")

    energy = contract["formal_cr001_energy_distinction"]
    _equal(energy["root_screen_quadrature_order"], 64, "root screen quadrature")
    _equal(energy["registered_quadrature_orders_per_axis"], REGISTERED_QUADRATURE, "registered energy ladder")
    _close(energy["registered_reference_time"], 0.25, "reference time", atol=0.0)
    _close(energy["registered_reference_energy"], 1.0, "reference energy", atol=0.0)
    _close(energy["registered_reference_energy_abs_tolerance"], 0.001, "reference energy tolerance", atol=0.0)
    _require(energy["order_64_root_screen_is_registered_reference_energy_acceptance"] is False, "order-64 root cannot replace CR001 energy acceptance")
    _require(energy["future_materialized_child_requires_registered_energy_ladder"] is True, "future child must rerun registered energy ladder")

    child = contract["child_delivery_scope"]
    for key in (
        "versioned_child_materialized",
        "child_candidate_id_assigned",
        "child_candidate_sha256_assigned",
        "save_load_contract_exists",
        "unified_velocity_api_registered",
        "velocity_export_ready",
        "successful_energy_root_or_path_replay_can_replace_canonical_delivery_without_promotion",
    ):
        _require(child[key] is False, f"unpromoted child delivery state became true: {key}")

    cr = contract["canonical_cr001"]
    _close(cr["nu"], constraints["nu"], "nu", atol=0.0)
    _equal(cr["physical_domain"], constraints["domain"]["physical"], "physical domain")
    _equal(cr["evaluation_box"], constraints["domain"]["evaluation_box"], "evaluation box")
    _equal(cr["support"], constraints["domain"]["support"], "support")
    _equal(cr["time_interval"], constraints["domain"]["time_interval"], "time interval")
    _equal(cr["force_mode"], constraints["forcing"]["mode"], "force mode")
    _equal(cr["force_bounds"], constraints["forcing"]["parameters"], "force bounds")
    _close(cr["reference_energy"], constraints["nontriviality"]["reference_energy"], "CR001 energy", atol=0.0)
    _close(cr["reference_energy_abs_tolerance"], constraints["nontriviality"]["reference_energy_abs_tolerance"], "CR001 energy tolerance", atol=0.0)
    val = constraints["validation"]
    _equal(cr["validation_seed"], val["seed"], "validation seed")
    _equal(cr["held_out_points"], val["held_out_points"], "held-out points")
    _equal(cr["validation_times"], REGISTERED_TIMES, "validation times")
    _equal(val["times"], REGISTERED_TIMES, "live validation times")
    _equal(cr["derivative_steps"], REGISTERED_STEPS, "derivative steps")
    _equal(val["derivative_steps"], REGISTERED_STEPS, "live derivative steps")
    _equal(cr["quadrature_orders_per_axis"], REGISTERED_QUADRATURE, "quadrature ladder")
    _equal(val["quadrature_orders_per_axis"], REGISTERED_QUADRATURE, "live quadrature ladder")
    thresholds = val["thresholds"]
    _close(cr["divergence_max"], thresholds["divergence_max"], "divergence max", atol=0.0)
    _close(cr["divergence_L2"], thresholds["divergence_L2"], "divergence L2", atol=0.0)
    _close(cr["pde_residual_max"], thresholds["pde_residual_max"], "momentum max", atol=0.0)
    _close(cr["pde_residual_L2"], thresholds["pde_residual_L2"], "momentum L2", atol=0.0)
    for key in (
        "free_or_residual_defined_force_allowed",
        "collapsed_velocity_success_allowed",
        "threshold_relaxation_under_same_experiment_allowed",
    ):
        _require(cr[key] is False, f"CR001 prohibition weakened: {key}")

    canonical = contract["canonical_delivery"]
    primary = delivery["primary_deliverable"]
    _equal(canonical["candidate_family"], CANONICAL_FAMILY, "contract canonical family")
    _equal(canonical["candidate_sha256"], CANONICAL_SHA, "contract canonical SHA")
    _equal(canonical["velocity_api"], CANONICAL_API, "contract canonical API")
    _equal(primary["candidate_family"], CANONICAL_FAMILY, "delivery family")
    _equal(primary["candidate_sha256"], CANONICAL_SHA, "delivery SHA")
    _equal(primary["api"], CANONICAL_API, "delivery API")
    _equal(project_status["candidate_family"], CANONICAL_FAMILY, "project family")
    _equal(project_status["candidate_sha256"], CANONICAL_SHA, "project SHA")
    _equal(project_status["velocity_api"], CANONICAL_API, "project API")
    _require(canonical["velocity_export_ready"] is True, "canonical velocity delivery must remain ready")
    _require(canonical["local_swirl_energy_child_replaces_canonical_candidate"] is False, "unpromoted child cannot replace canonical candidate")
    _require(canonical["local_swirl_energy_child_replaces_canonical_velocity_api"] is False, "unpromoted child cannot replace canonical API")

    states = contract["independent_truth_states"]
    _require(states["canonical_velocity_export_ready"] is True, "canonical export readiness must remain true")
    for key in (
        "local_swirl_energy_child_velocity_export_ready",
        "visualization_ready",
        "visual_correspondence_verified",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    ):
        _require(states[key] is False, f"independent truth state became true: {key}")
    live_states = project_status["states"]
    _require(live_states["velocity_export_ready"] is True, "live canonical export readiness changed")
    for key in ("visualization_ready", "visual_correspondence_verified", "pde_validated", "paper_exact", "openai_field_identified", "blowup_proved"):
        _require(live_states[key] is False, f"live project truth state unexpectedly true: {key}")

    future = contract["future_research_rules"]
    for key, value in future.items():
        _require(value is True, f"future research guard weakened: {key}")

    forbidden = contract["forbidden_inferences"]
    _require(len(forbidden) >= 8, "forbidden inference registry unexpectedly shortened")
    _require(all(item.get("allowed") is False for item in forbidden), "all registered forbidden inferences must stay disallowed")

    return {
        "task_id": TASK_ID,
        "status": "PASS",
        "energy_root_found": True,
        "beta": root["beta"],
        "frozen_material_path_aggregate_change": "reported_zero_at_receipt_precision",
        "child_velocity_export_ready": False,
        "canonical_velocity_export_ready": True,
        "pde_validated": False,
        "visual_correspondence_verified": False,
        "paper_exact": False,
    }


def audit_repository() -> dict[str, Any]:
    return audit_contract(load_contract(), load_constraints(), load_delivery(), load_project_status())


def main() -> int:
    print(json.dumps(audit_repository(), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
