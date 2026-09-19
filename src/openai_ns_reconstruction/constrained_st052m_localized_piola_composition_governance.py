"""Fail-closed CR002 audit for the PR #544 localized Piola composition.

PR #544 applies one fixed callable contravariant Piola taper to the diagnostic
ST052-M + kappa=.05 redistribution parent.  This audit keeps the coordinate map,
common energy scale, sampled morphology/divergence screen, parent material-path
receipt, parent PDE receipt, delivery identity, and scientific claim states
separate.

It changes no velocity, pressure, force, threshold, validation sample, or
canonical delivery state.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Mapping

TASK_ID = "CR002-ST052M-LOCALIZED-PIOLA-COMPOSITION-087"
_INTEGRATION_HEAD = "b9246717cccafa2ea6595f9393d4effdf67d2854"
_TAPER_HEAD = "093c7171cd61c6bd439afa30b2da69598a02d182"
_TAPER_BASE_HEAD = "779ffca71066e2864496d37de55a7aafc45d6f57"
_PARENT_GOV_TASK = "CR002-ST052M-REDISTRIBUTION-EVIDENCE-COMPOSITION-086"
_CANONICAL_FAMILY = "eq45_supported_velocity_candidate_v1"
_CANONICAL_SHA = "2fdff812c22131d56eac1b7d6e455187ff3207500c6385aa9abf282a8e3d7b1d"
_CANONICAL_API = "openai_ns_reconstruction.eq45_supported_delivery:velocity"
_ALLOWED_CLASSES = {"user_requirement", "public_source_fact", "autonomous_design", "pending_unknown"}
_REGISTERED_TIMES = [0.25, 0.3125, 0.4375, 0.5625, 0.6875, 0.75]
_REGISTERED_STEPS = [0.02, 0.01, 0.005]
_REGISTERED_QUADRATURE = [24, 48, 96]


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _load_json(path: Path) -> dict[str, Any]:
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise ValueError(f"expected JSON object in {path}")
    return data


def load_contract() -> dict[str, Any]:
    return _load_json(_repo_root() / "configs" / "st052m_localized_piola_composition_governance.json")


def load_parent_governance() -> dict[str, Any]:
    return _load_json(_repo_root() / "configs" / "st052m_redistribution_evidence_composition_governance.json")


def load_constraints() -> dict[str, Any]:
    return _load_json(_repo_root() / "configs" / "constraints.json")


def load_project_status() -> dict[str, Any]:
    return _load_json(_repo_root() / "project_status.json")


def load_velocity_delivery_contract() -> dict[str, Any]:
    return _load_json(_repo_root() / "configs" / "velocity_delivery_contract.json")


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
    parent: Mapping[str, Any],
    project_status: Mapping[str, Any],
    delivery: Mapping[str, Any],
) -> dict[str, Any]:
    _equal(contract["schema_version"], 1, "schema_version")
    _equal(contract["task_id"], TASK_ID, "task_id")

    snap = contract["snapshot"]
    _equal(snap["active_integration_branch"], "codex/cr001-constraints", "integration branch")
    _equal(snap["active_integration_head"], _INTEGRATION_HEAD, "integration head")
    _equal(snap["parent_governance_task"], _PARENT_GOV_TASK, "parent governance task")
    _equal(snap["taper_pr"], 544, "taper PR")
    _equal(snap["taper_pr_head"], _TAPER_HEAD, "taper head")
    _equal(snap["taper_pr_base_head"], _TAPER_BASE_HEAD, "taper base head")
    _equal(snap["preregister_issue"], 542, "preregister issue")
    _equal(snap["dedicated_run"], 35416797383, "dedicated CI")
    _equal(snap["standard_run"], 35416797305, "standard CI")
    _equal(snap["inherited_st052_replay_run"], 35416797276, "ST052 replay CI")
    _equal(snap["artifact_id"], 10576196732, "artifact id")
    _equal(
        snap["artifact_digest"],
        "sha256:1ff9209f249ca049ed11a647ab2383296bfbcd28c577fc5e312e0c38960a5057",
        "artifact digest",
    )

    _equal(set(contract["allowed_source_classes"]), _ALLOWED_CLASSES, "allowed source classes")
    _equal(set(contract["source_classification"].values()), _ALLOWED_CLASSES, "source classes represented")

    _equal(parent["task_id"], _PARENT_GOV_TASK, "loaded parent governance task")
    _require(parent["child_material_path_evidence"]["positive_material_path_transfer_under_frozen_protocol"] is True,
             "governed parent must retain its positive material-path receipt")
    _require(parent["child_delivery_identity"]["materialized_versioned_candidate_exists"] is False,
             "redistributed diagnostic parent must remain unmaterialized")
    _require(parent["parent_pde_evidence"]["may_transfer_to_transformed_child"] is False,
             "ST052-M PDE receipt must not transfer even to the redistributed parent")

    ident = contract["composition_identity"]
    _equal(ident["grandparent"], "ST052-M", "grandparent")
    _equal(ident["diagnostic_parent"], "ST052-M plus frozen kappa=.05 redistribution", "diagnostic parent")
    for key in (
        "diagnostic_parent_versioned_save_load_candidate",
        "diagnostic_child_versioned_save_load_candidate",
        "diagnostic_child_candidate_id_assigned",
        "diagnostic_child_candidate_sha256_assigned",
        "diagnostic_child_unified_velocity_api_registered",
        "diagnostic_child_canonical_repository_candidate",
        "diagnostic_child_production_selected",
    ):
        _require(ident[key] is False, f"unmaterialized composition state became true: {key}")

    taper = contract["frozen_taper"]
    _close(taper["tau"], 0.05, "tau", atol=0.0)
    _close(taper["axial_support"], 2.0, "axial support", atol=0.0)
    _equal(taper["window_abs_z_over_2"], [0.50, 0.82], "taper window")
    _equal(taper["q_range"], [0.0, 1.0], "q range")
    expected_a = [1.0, 1.0 + taper["tau"]]
    _equal(taper["a_range"], expected_a, "a range")
    _require(min(taper["a_range"]) > 0.0 and taper["orientation_preserving"] is True,
             "localized map must remain orientation preserving")
    _require(taper["map_identity_for_abs_z_over_2_le_0p50"] is True,
             "central coordinate-map identity must remain explicit")
    _require(taper["map_identity_near_axial_support_endpoints"] is True,
             "endpoint coordinate-map identity must remain explicit")
    _close(taper["taper_energy_scale_reported"], 0.9929528816, "reported taper scale", atol=5e-10)
    _require(0.0 < taper["taper_energy_scale_reported"] < 1.0, "taper scale must remain positive and non-unit")
    _require(taper["parameter_scan_performed"] is False, "PR544 must remain a fixed taper screen")
    _equal(taper["new_velocity_degree_count"], 1, "new velocity degree count")

    rep = contract["representation_audit"]
    for key in (
        "transform_is_applied_to_callable_parent_not_trilinear_sample_grid",
        "continuous_piola_identity_conditionally_preserves_divergence",
        "common_constant_scale_preserves_divergence_identity",
        "smooth_map_introduces_no_new_cartesian_axis_singularity",
        "a_ge_1_and_z_unchanged_preserve_parent_zero_extension_outside_registered_support",
        "central_final_velocity_is_positive_common_scale_times_parent",
        "central_direction_and_normalized_geometry_can_be_unchanged_under_positive_scale",
    ):
        _require(rep[key] is True, f"representation invariant must remain true: {key}")
    for key in (
        "nonzero_support_region_is_proved_identical_to_parent",
        "coordinate_map_identity_implies_final_velocity_identity_after_common_scale",
        "central_velocity_amplitude_is_exactly_unchanged",
        "piola_divergence_identity_is_navier_stokes_invariance",
        "common_energy_scale_is_navier_stokes_invariance",
        "pressure_can_be_inherited_without_rebuild",
        "restricted_force_fit_can_be_inherited_without_rebuild",
        "momentum_residual_can_be_inherited_from_parent",
    ):
        _require(rep[key] is False, f"forbidden representation inference became true: {key}")

    sampled = contract["sampled_taper_evidence"]
    _equal(sampled["times"], [0.25, 0.50, 0.75], "sampled times")
    _equal(sampled["energy_quadrature_order"], 64, "screen energy quadrature")
    _require(sampled["energy_rescaling_restores_sampled_parent_reference_energy"] is True,
             "screen must retain its sampled energy-rescaling receipt")
    _require(sampled["registered_CR001_energy_quadrature_ladder_replayed"] is False,
             "order-64 screen cannot become registered energy acceptance")
    _equal(sampled["central_radial_vorticity_rms_relative"], [0.0, 0.0, 0.0], "central morphology rows")
    _close(sampled["central_angular_rate_retention_r_0p6"], taper["taper_energy_scale_reported"],
           "central r=.6 retention vs common scale", atol=2e-6)
    _close(sampled["central_angular_rate_retention_r_0p9"], taper["taper_energy_scale_reported"],
           "central r=.9 retention vs common scale", atol=2e-6)
    _close(sampled["sampled_support_max_abs"], 0.0, "sampled support max", atol=0.0)
    _close(sampled["cartesian_fd_divergence_max"], 2.613e-9, "custom FD divergence", atol=5e-13)
    _close(sampled["central_coordinate_map_identity_error"], 0.0, "central map identity", atol=0.0)
    _equal(sampled["response_rank"], 3, "response rank")
    _close(sampled["response_condition_number"], 1.109075, "response condition", atol=5e-7)
    _close(sampled["response_max_pairwise_abs_cosine"], 0.084803, "response cosine", atol=5e-7)
    _require(sampled["clean_localized_taper_capacity"] is True, "clean capacity receipt must remain true")
    for key in (
        "formal_CR001_divergence_acceptance_replayed",
        "formal_CR001_full_momentum_replayed",
        "continuum_morphology_certificate",
        "visual_correspondence_verified",
        "pde_validated",
    ):
        _require(sampled[key] is False, f"sampled evidence was over-promoted: {key}")

    ownership = contract["evidence_ownership_and_noninheritance"]
    for key in (
        "st052m_parent_pde_receipt_belongs_to_taper_child",
        "redistributed_parent_material_paths_belong_to_taper_child",
        "redistributed_parent_positive_material_path_transfer_is_taper_material_path_evidence",
        "taper_material_paths_integrated",
        "taper_render_or_public_reference_comparison_completed",
        "rank3_response_selects_production_tau",
        "central_map_identity_proves_central_velocity_identity",
        "sampled_order64_energy_scale_proves_CR001_energy_acceptance",
        "custom_fd_divergence_below_1e5_proves_CR001_divergence_acceptance",
        "target_free_morphology_pass_proves_OpenAI_correspondence",
        "green_ci_is_scientific_acceptance",
    ):
        _require(ownership[key] is False, f"forbidden evidence inheritance became true: {key}")

    routing = contract["routing"]
    _require(routing["localized_taper_is_expression_capacity_evidence"] is True,
             "localized taper must remain capacity evidence")
    _require(routing["localized_taper_supersedes_governed_kappa_0p05_materialization_priority"] is False,
             "capacity screen cannot supersede governed parent handoff")
    _require(routing["combined_taper_child_ready_for_pressure_or_restricted_force_reconstruction"] is False,
             "PDE reconstruction cannot be promoted before taper selection")
    _require(routing["combined_taper_child_ready_for_production_selection"] is False,
             "production selection cannot be promoted")
    _require(routing["independent_material_path_or_render_evidence_required_before_taper_selection"] is True,
             "selection must require an independent downstream lane")
    _require(routing["pde_failure_blocks_truth_bounded_diagnostic_callable_or_future_save_load"] is False,
             "PDE failure must not block truth-bounded velocity delivery")

    cr = contract["canonical_cr001"]
    _close(cr["nu"], constraints["nu"], "nu", atol=0.0)
    _equal(cr["physical_domain"], constraints["domain"]["physical"], "physical domain")
    _equal(cr["evaluation_box"], constraints["domain"]["evaluation_box"], "evaluation box")
    _equal(cr["support"], constraints["domain"]["support"], "support")
    _equal(cr["time_interval"], constraints["domain"]["time_interval"], "time interval")
    _equal(cr["force_mode"], constraints["forcing"]["mode"], "force mode")
    _equal(cr["force_bounds"], constraints["forcing"]["parameters"], "force bounds")
    _close(cr["reference_energy"], constraints["nontriviality"]["reference_energy"], "reference energy", atol=0.0)
    _close(cr["reference_energy_abs_tolerance"], constraints["nontriviality"]["reference_energy_abs_tolerance"],
           "reference energy tolerance", atol=0.0)
    val = constraints["validation"]
    _equal(cr["validation_seed"], val["seed"], "validation seed")
    _equal(cr["held_out_points"], val["held_out_points"], "held-out points")
    _equal(cr["validation_times"], _REGISTERED_TIMES, "registered validation times")
    _equal(val["times"], _REGISTERED_TIMES, "constraints validation times")
    _equal(cr["energy_quadrature_orders_per_axis"], _REGISTERED_QUADRATURE, "registered quadrature ladder")
    _equal(val["quadrature_orders_per_axis"], _REGISTERED_QUADRATURE, "constraints quadrature ladder")
    _equal(cr["derivative_steps"], _REGISTERED_STEPS, "registered derivative ladder")
    _equal(val["derivative_steps"], _REGISTERED_STEPS, "constraints derivative ladder")
    _close(cr["divergence_max"], val["thresholds"]["divergence_max"], "divergence max", atol=0.0)
    _close(cr["divergence_L2"], val["thresholds"]["divergence_L2"], "divergence L2", atol=0.0)
    _close(cr["pde_residual_max"], val["thresholds"]["pde_residual_max"], "momentum max", atol=0.0)
    _close(cr["pde_residual_L2"], val["thresholds"]["pde_residual_L2"], "momentum L2", atol=0.0)
    _require(cr["free_or_residual_defined_force_allowed"] is False, "free force must remain forbidden")
    _require(cr["collapsed_velocity_success_allowed"] is False, "amplitude collapse must remain forbidden")
    _require(cr["threshold_relaxation_under_same_experiment_allowed"] is False, "threshold relaxation must remain forbidden")

    canon = contract["canonical_delivery"]
    _equal(canon["candidate_family"], _CANONICAL_FAMILY, "canonical family")
    _equal(canon["candidate_sha256"], _CANONICAL_SHA, "canonical SHA")
    _equal(canon["velocity_api"], _CANONICAL_API, "canonical API")
    _require(canon["velocity_export_ready"] is True, "canonical callable delivery must stay ready")
    _require(canon["localized_taper_replaces_canonical_candidate"] is False, "taper cannot replace canonical candidate")
    _require(canon["localized_taper_replaces_canonical_velocity_api"] is False, "taper cannot replace canonical API")
    _equal(project_status["candidate_family"], _CANONICAL_FAMILY, "project-status family")
    _equal(project_status["candidate_sha256"], _CANONICAL_SHA, "project-status SHA")
    _equal(project_status["velocity_api"], _CANONICAL_API, "project-status API")
    _equal(delivery["primary_deliverable"]["candidate_family"], _CANONICAL_FAMILY, "delivery family")
    _equal(delivery["primary_deliverable"]["candidate_sha256"], _CANONICAL_SHA, "delivery SHA")
    _equal(delivery["primary_deliverable"]["api"], _CANONICAL_API, "delivery API")

    states = contract["independent_truth_states"]
    _require(states["velocity_export_ready"] is True, "canonical export readiness must stay true")
    _require(project_status["states"]["velocity_export_ready"] is True, "project-status export readiness must stay true")
    _require(delivery["claim_status"]["velocity_export_ready"] is True, "delivery export readiness must stay true")
    for key in (
        "localized_taper_child_velocity_export_ready",
        "visualization_ready",
        "visual_correspondence_verified",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    ):
        _require(states[key] is False, f"independent truth state became true: {key}")
    for key in ("visualization_ready", "pde_validated", "visual_correspondence_verified", "paper_exact", "openai_field_identified", "blowup_proved"):
        _require(project_status["states"][key] is False, f"project status scientific state became true: {key}")
        _require(delivery["claim_status"][key] is False, f"delivery scientific state became true: {key}")

    req = contract["future_selection_requirements"]
    for key in (
        "assign_new_candidate_id_and_sha256_if_materialized",
        "create_reproducible_save_load_recipe_if_materialized",
        "register_explicit_velocity_x_y_z_t_api_or_adapter_if_materialized",
        "recheck_axis_near_axis_support_core_and_registered_reference_energy",
        "rerun_material_paths_on_exact_taper_child_before_material_path_claim",
        "run_public_reference_render_comparison_and_resolution_stability_before_visual_correspondence",
        "rebuild_or_refit_compatible_pressure_before_pde_claim",
        "rebuild_or_refit_only_existing_restricted_two_parameter_force_family",
        "do_not_use_free_or_residual_defined_force",
        "run_fresh_full_momentum_and_divergence_validation_on_exact_materialized_child",
        "run_registered_CR001_acceptance_contract_before_formal_pde_promotion",
    ):
        _require(req[key] is True, f"future selection gate must remain required: {key}")

    return {
        "task_id": TASK_ID,
        "status": "pass",
        "taper_child_materialized": False,
        "central_coordinate_map_identity": True,
        "central_final_velocity_identity": False,
        "taper_material_paths_integrated": False,
        "pde_validated": False,
        "visual_correspondence_verified": False,
        "canonical_velocity_unchanged": True,
    }


def audit_current_repository() -> dict[str, Any]:
    return audit_contract(
        load_contract(),
        load_constraints(),
        load_parent_governance(),
        load_project_status(),
        load_velocity_delivery_contract(),
    )


if __name__ == "__main__":
    print(json.dumps(audit_current_repository(), indent=2, sort_keys=True))
