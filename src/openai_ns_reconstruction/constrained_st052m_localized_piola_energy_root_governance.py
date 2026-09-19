"""Fail-closed CR002 audit for the ST052-M localized-Piola energy-root screen.

PR #549 is intentionally a scientific negative result with green execution: the
preregistered shoulder-compensation family has no energy-neutral root on its
frozen alpha bracket.  This module keeps that statement narrow.  In particular,
it prevents evidence from the scaled #544 taper or its #546 material-path replay
from being assigned to a nonexistent #549 child, and it keeps the order-64 root
screen separate from the registered CR001 energy/validation contract.

This module changes no velocity, pressure, forcing, threshold, candidate, or
canonical delivery state.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Mapping

TASK_ID = "CR002-ST052M-LOCALIZED-PIOLA-ENERGY-ROOT-SCOPE-087"
INTEGRATION_HEAD = "b9246717cccafa2ea6595f9393d4effdf67d2854"
TAPER_HEAD = "093c7171cd61c6bd439afa30b2da69598a02d182"
PATH_HEAD = "1bc3f912e65e88c002002dbce47d597933706e78"
ROOT_HEAD = "62e8c170427d5d830d7f897ba31768e0fc4ce56a"
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
    return _load_json(_repo_root() / "configs" / "st052m_localized_piola_energy_root_governance.json")


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
    _equal(snap["scaled_taper_pr"], 544, "scaled taper PR")
    _equal(snap["scaled_taper_head"], TAPER_HEAD, "scaled taper head")
    _equal(snap["scaled_taper_material_path_pr"], 546, "material-path PR")
    _equal(snap["scaled_taper_material_path_head"], PATH_HEAD, "material-path head")
    _equal(snap["energy_root_pr"], 549, "energy-root PR")
    _equal(snap["energy_root_head"], ROOT_HEAD, "energy-root head")
    _equal(snap["energy_root_standard_run"], 35418075526, "energy-root standard run")
    _equal(snap["energy_root_dedicated_run"], 35418075593, "energy-root dedicated run")

    _equal(set(contract["allowed_source_classes"]), ALLOWED_CLASSES, "allowed source classes")
    classes = contract["source_classification"]
    _equal(set(classes.values()), ALLOWED_CLASSES, "all source classes represented")
    for key in (
        "st052m_and_frozen_kappa_0p05_redistribution",
        "localized_piola_tau_0p05_and_tip_window",
        "shoulder_window_and_signed_compensation_family",
        "alpha_energy_bracket_0_to_4",
        "order_64_energy_root_screen",
    ):
        _equal(classes[key], "autonomous_design", f"source class for {key}")
    for key in (
        "openai_hidden_numerical_velocity_or_compensation_map",
        "openai_hidden_parameters_time_camera_seed_mapping",
    ):
        _equal(classes[key], "pending_unknown", f"source class for {key}")

    scaled = contract["scaled_taper_evidence"]
    _close(scaled["tau"], 0.05, "scaled taper tau", atol=0.0)
    _equal(scaled["tip_window_abs_z_over_2"], [0.50, 0.82], "scaled taper window")
    _close(scaled["post_piola_common_scale"], 0.9929528816318394, "scaled taper normalization", atol=2e-12)
    _require(scaled["central_coordinate_map_identity"] is True, "central map identity must remain recorded")
    _require(scaled["central_velocity_pointwise_identity"] is False, "global scale prevents pointwise central velocity identity")
    _require(scaled["material_path_replay_exists"] is True, "scaled child path replay must remain recorded")
    _close(scaled["material_path_mean_turns_relative_change"], -0.007387, "mean turns change", atol=5e-7)
    _close(scaled["material_path_contraction_magnitude_relative_change"], -0.007091, "contraction change", atol=5e-7)
    _close(scaled["material_path_mean_pair_axial_separation_relative_change"], -0.009409, "pair separation change", atol=5e-7)
    _require(scaled["material_path_evidence_transfers_to_hypothetical_unscaled_energy_root_child"] is False, "scaled-child path evidence cannot transfer")
    _require(scaled["eulerian_morphology_evidence_transfers_to_hypothetical_unscaled_energy_root_child"] is False, "scaled-child morphology cannot transfer")

    root = contract["energy_root_experiment"]
    _equal(root["task_id"], "CR003-ST052M-ENERGY-NEUTRAL-TAPER-081", "root task id")
    _close(root["tau"], 0.05, "root tau", atol=0.0)
    _equal(root["tip_window_abs_z_over_2"], [0.50, 0.82], "root tip window")
    _equal(root["shoulder_window_abs_z_over_2"], [0.36, 0.49], "root shoulder window")
    _equal(root["alpha_E_bracket"], [0.0, 4.0], "alpha bracket")
    _equal(root["energy_quadrature_order"], 64, "root quadrature order")
    _require(root["post_piola_common_scale_allowed"] is False, "root family cannot silently restore energy by common scaling")
    _require(root["parameter_scan_performed"] is False, "root experiment must remain fixed-family evidence")
    _require(root["energy_only_dependent_root_solve"] is True, "alpha_E may depend only on frozen energy equation")
    _close(root["control_energy"], 1.0000000000014384, "control energy", atol=2e-12)
    _close(root["endpoint_energy_residual_alpha_0"], 0.014244634697393677, "alpha=0 energy residual", atol=2e-12)
    _close(root["endpoint_energy_residual_alpha_4"], 0.10321550869365947, "alpha=4 energy residual", atol=2e-12)
    _require(root["endpoint_energy_residual_alpha_0"] > 0.0 and root["endpoint_energy_residual_alpha_4"] > 0.0, "both frozen bracket endpoints must retain positive energy residual")
    for key in (
        "root_exists_in_preregistered_bracket",
        "alpha_E_selected",
        "clean_energy_neutral_taper_capacity",
        "morphology_evaluated_for_energy_neutral_child",
        "material_paths_integrated_for_energy_neutral_child",
        "versioned_child_materialized",
        "child_candidate_id_assigned",
        "child_candidate_sha256_assigned",
        "save_load_contract_exists",
        "unified_velocity_api_registered",
        "velocity_export_ready",
    ):
        _require(root[key] is False, f"nonexistent root child state became true: {key}")

    failure = contract["scientific_failure_scope"]
    _require(failure["workflow_execution_success"] is True, "successful experiment execution must remain recorded")
    _require(failure["scientific_energy_neutral_crossing_found"] is False, "scientific crossing must remain false")
    _require(failure["no_root_is_negative_capacity_evidence_for_frozen_family"] is True, "narrow negative result must remain recorded")
    for key in (
        "no_root_proves_no_energy_neutral_localized_piola_transform_exists",
        "no_root_proves_no_other_divergence_preserving_compensation_can_work",
        "no_root_identifies_any_openai_hidden_parameter",
        "green_ci_is_scientific_pass",
        "same_experiment_may_relax_or_expand_bracket_after_observing_failure",
    ):
        _require(failure[key] is False, f"forbidden failure inference became true: {key}")
    _require(failure["changing_sign_window_bracket_or_compensation_family_after_result_requires_new_experiment_identity"] is True, "post-result redesign must get new experiment identity")

    rep = contract["representation_and_pde_scope"]
    _require(rep["continuous_callable_piola_structure_may_preserve_divergence_under_its_mathematical_assumptions"] is True, "Piola structural scope should remain explicit")
    for key in (
        "piola_divergence_structure_is_navier_stokes_invariance",
        "common_energy_scaling_is_navier_stokes_invariance",
        "parent_pressure_receipt_may_transfer",
        "parent_restricted_force_receipt_may_transfer",
        "parent_momentum_receipt_may_transfer",
        "custom_or_local_fd_divergence_is_formal_cr001_acceptance",
        "target_free_morphology_or_render_can_promote_pde_validated",
        "target_free_morphology_or_render_can_promote_visual_correspondence_verified",
    ):
        _require(rep[key] is False, f"forbidden representation/PDE inference became true: {key}")

    energy = contract["formal_cr001_energy_distinction"]
    _equal(energy["root_screen_quadrature_order"], 64, "root screen quadrature")
    _equal(energy["registered_quadrature_orders_per_axis"], REGISTERED_QUADRATURE, "registered energy ladder")
    _close(energy["registered_reference_time"], 0.25, "reference time", atol=0.0)
    _close(energy["registered_reference_energy"], 1.0, "reference energy", atol=0.0)
    _close(energy["registered_reference_energy_abs_tolerance"], 0.001, "reference energy tolerance", atol=0.0)
    _require(energy["order_64_root_screen_is_registered_reference_energy_acceptance"] is False, "order-64 screen cannot become CR001 acceptance")
    _require(energy["future_materialized_energy_root_child_requires_registered_energy_ladder"] is True, "future child must rerun registered energy ladder")

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
    _require(canonical["energy_root_experiment_replaces_canonical_candidate"] is False, "failed root experiment cannot replace canonical candidate")
    _require(canonical["energy_root_experiment_replaces_canonical_velocity_api"] is False, "failed root experiment cannot replace canonical API")

    states = contract["independent_truth_states"]
    _require(states["velocity_export_ready"] is True, "canonical velocity_export_ready must remain true")
    for key in (
        "energy_root_child_velocity_export_ready",
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
    for key in (
        "a_new_compensation_family_must_preregister_new_identity_before_evaluation",
        "a_future_root_requires_new_candidate_id_and_sha_before_delivery_claim",
        "a_future_root_requires_reproducible_save_load_and_velocity_adapter",
        "a_future_root_requires_fresh_axis_support_core_and_registered_energy_checks",
        "a_future_root_requires_fresh_material_paths_before_path_claims",
        "a_future_root_requires_fresh_visual_diagnostics_before_visual_claims",
        "a_future_root_requires_compatible_pressure_and_existing_restricted_force_rebuild_before_pde_claim",
        "a_future_root_requires_fresh_full_momentum_and_divergence_on_exact_child_before_pde_claim",
        "pde_failure_does_not_block_truth_bounded_callable_velocity_delivery",
    ):
        _require(future[key] is True, f"future research guard weakened: {key}")

    forbidden = contract["forbidden_inferences"]
    _require(len(forbidden) >= 8, "forbidden inference registry unexpectedly shortened")
    _require(all(item.get("allowed") is False for item in forbidden), "all registered forbidden inferences must stay disallowed")

    return {
        "task_id": TASK_ID,
        "status": "PASS",
        "energy_root_found": False,
        "scientific_result": "negative_capacity_evidence_for_frozen_family_only",
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
