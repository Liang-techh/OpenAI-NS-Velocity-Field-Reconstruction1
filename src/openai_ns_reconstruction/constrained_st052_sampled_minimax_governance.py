"""Fail-closed CR002 governance for ST052-M sampled-minimax evidence.

PR #508 freezes a replayable experimental child after finite sampled optimization
and reports a real reduction in fresh sampled momentum maximum together with a
small worsening in fixed-time volume-L2.  This module prevents those diagnostics
from replacing the preregistered CR001 acceptance sample or from promoting PDE,
visual, paper-exact, OpenAI-field, or canonical-delivery state.

No velocity, candidate coefficient, pressure, forcing, threshold, validation
sample, or active scientific route is changed here.
"""
from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
from typing import Any, Mapping

TASK_ID = "CR002-ST052-SAMPLED-MINIMAX-VALIDATION-SCOPE-085"
_INTEGRATION_HEAD = "3c028688b3ff3d2cece5f734a98f392725585493"
_PR508_HEAD = "b3b8bfdbe1077f9ec967d158602951997d81e17d"
_PR508_BASE = "4b784f1b8457af2ead49295631d834d4e882000b"
_PR508_TEST_RUN = 35404426881
_PR508_REPLAY_RUN = 35404426929
_PR508_ARTIFACT = 10572166706
_PR508_ARTIFACT_DIGEST = "sha256:4f379950d9aa9267e03feda49abdf9808635d490822fcefdafaf8bd1f9311239"
_ST052_SHA = "e078e753fab38ebfa0284d28ba64d26cb8538b700a3afc7705668849c43e12da"
_ST051B_SHA = "0071106ef10a5d77b620b942becc00b58c8a6765fd5fb7bc4dcd195ac65c970d"
_MODIFIER_SHA = "4d8b6e92f47151a7eed112985779a58fdf1188bd30c35c3a59263910629f1be9"
_CANONICAL_FAMILY = "eq45_supported_velocity_candidate_v1"
_CANONICAL_SHA = "2fdff812c22131d56eac1b7d6e455187ff3207500c6385aa9abf282a8e3d7b1d"
_CANONICAL_API = "openai_ns_reconstruction.eq45_supported_delivery:velocity"
_ALLOWED_CLASSES = {"user_requirement", "public_source_fact", "autonomous_design", "pending_unknown"}
_REGISTERED_TIMES = [0.25, 0.3125, 0.4375, 0.5625, 0.6875, 0.75]


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _load_json(path: Path) -> dict[str, Any]:
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise ValueError(f"expected object in {path}")
    return data


def load_contract() -> dict[str, Any]:
    return _load_json(_repo_root() / "configs" / "st052_sampled_minimax_validation_governance.json")


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
    _require(actual == expected, f"{label} drifted: {actual!r} != {expected!r}")


def _true(value: Any, label: str) -> None:
    _equal(value, True, label)


def _false(value: Any, label: str) -> None:
    _equal(value, False, label)


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
    _equal(frozen["validation_times"], _REGISTERED_TIMES, "registered six validation times")
    _equal(frozen["derivative_steps"], validation["derivative_steps"], "derivative steps")
    _equal(frozen["divergence_max"], thresholds["divergence_max"], "divergence max")
    _equal(frozen["divergence_L2"], thresholds["divergence_L2"], "divergence L2")
    _equal(frozen["pde_residual_max"], thresholds["pde_residual_max"], "momentum max threshold")
    _equal(frozen["pde_residual_L2"], thresholds["pde_residual_L2"], "momentum L2 threshold")
    _true(frozen["formal_momentum_gate_requires_max_and_L2"], "joint momentum acceptance")
    _false(frozen["free_or_residual_defined_force_allowed"], "free-force allowance")
    _false(frozen["collapsed_velocity_success_allowed"], "collapsed-velocity shortcut")
    _false(frozen["threshold_relaxation_under_same_experiment_allowed"], "threshold relaxation")
    _require(
        "No residual-dependent basis or pointwise free force" in forcing["restriction"],
        "free-force guard drifted",
    )
    _require("reject collapsed candidates" in nontrivial["enforcement"], "nontriviality guard drifted")


def audit_st052_sampled_minimax_governance(
    contract: Mapping[str, Any] | None = None,
    constraints: Mapping[str, Any] | None = None,
    project_status: Mapping[str, Any] | None = None,
    delivery_contract: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Audit ST052-M evidence boundaries against live CR001/delivery truth state."""
    contract = deepcopy(dict(load_contract() if contract is None else contract))
    constraints = deepcopy(dict(load_constraints() if constraints is None else constraints))
    project_status = deepcopy(dict(load_project_status() if project_status is None else project_status))
    delivery_contract = deepcopy(
        dict(load_velocity_delivery_contract() if delivery_contract is None else delivery_contract)
    )

    _equal(contract.get("schema_version"), 1, "contract schema")
    _equal(contract.get("task_id"), TASK_ID, "task id")

    snapshot = contract["snapshot"]
    _equal(snapshot["active_integration_branch"], "codex/cr001-constraints", "integration branch")
    _equal(snapshot["active_integration_head"], _INTEGRATION_HEAD, "integration head")
    _equal(snapshot["observed_pr"], 508, "observed PR")
    _equal(snapshot["observed_pr_head"], _PR508_HEAD, "PR508 head")
    _equal(snapshot["observed_pr_base"], _PR508_BASE, "PR508 base")
    _equal(snapshot["observed_pr_state"], "open_draft_unintegrated_experimental_candidate", "PR508 state")
    _equal(snapshot["standard_tests_run"], _PR508_TEST_RUN, "PR508 standard run")
    _equal(snapshot["standard_tests_conclusion"], "success", "PR508 standard run conclusion")
    _equal(snapshot["dedicated_replay_run"], _PR508_REPLAY_RUN, "PR508 replay run")
    _equal(snapshot["dedicated_replay_conclusion"], "success", "PR508 replay conclusion")
    _equal(snapshot["dedicated_artifact_id"], _PR508_ARTIFACT, "PR508 replay artifact")
    _equal(snapshot["dedicated_artifact_digest"], _PR508_ARTIFACT_DIGEST, "PR508 replay digest")

    classes = set(contract["allowed_source_classes"])
    _equal(classes, _ALLOWED_CLASSES, "allowed source classes")
    for key, value in contract["source_classification"].items():
        _require(value in classes, f"unknown source class for {key}: {value}")
    _equal(
        contract["source_classification"]["callable_saveable_velocity_x_y_z_t"],
        "user_requirement",
        "velocity objective source class",
    )
    _equal(
        contract["source_classification"]["public_openai_visualization_observables"],
        "public_source_fact",
        "public visualization source class",
    )
    _equal(
        contract["source_classification"]["st052_sampled_minimax_objective_and_training_pool"],
        "autonomous_design",
        "ST052 optimization source class",
    )
    _equal(
        contract["source_classification"]["openai_hidden_numerical_velocity_and_coefficients"],
        "pending_unknown",
        "hidden OpenAI field source class",
    )

    candidate = contract["candidate_identity"]
    _equal(candidate["candidate_id"], "ST052-M", "candidate id")
    _equal(candidate["parent_id"], "ST051-B", "parent id")
    _equal(candidate["parent_sha256"], _ST051B_SHA, "parent SHA")
    _equal(candidate["raw_candidate_sha256"], _ST052_SHA, "ST052 SHA")
    _equal(candidate["modifier_sha256"], _MODIFIER_SHA, "modifier SHA")
    _true(candidate["frozen_before_fresh_validation"], "pre-validation freeze")
    _true(candidate["replay_recipe_present"], "replay recipe")
    _true(candidate["experimental_materialized_child"], "materialized experimental child")
    _false(candidate["canonical_repository_candidate"], "canonical ST052 promotion")
    _false(candidate["production_candidate_selected"], "production ST052 selection")
    _false(candidate["paper_exact"], "ST052 paper exactness")
    _false(candidate["openai_field_identified"], "ST052 OpenAI identity")

    selection = contract["optimization_and_selection"]
    _equal(selection["training_seed"], 9175250, "training seed")
    _equal(selection["training_pool_points"], 23079, "training pool size")
    _equal(selection["active_set_sizes"], [528, 833], "active set sizes")
    _equal(selection["optimizer_stages"], 2, "optimizer stage count")
    _false(selection["optimizer_stages_converged"], "optimizer convergence claim")
    _true(selection["optimizer_stages_hit_wall_budget"], "optimizer wall-stop record")
    _true(selection["pre_holdout_feasibility_restoration_allowed"], "pre-holdout repair policy")
    _true(selection["pre_holdout_feasibility_restoration_used"], "pre-holdout repair used")
    _equal(selection["restoration_corrections_used"], 1, "repair correction count")
    _require(selection["restoration_min_constraint"] >= -1e-7, "repair feasibility tolerance failed")
    _true(selection["restored_checkpoint_is_finite_training_feasible"], "finite training feasibility")
    _false(selection["restored_checkpoint_is_local_or_global_optimum"], "optimizer optimality laundering")
    _false(selection["sampled_epigraph_is_continuum_bound"], "sampled epigraph continuum promotion")
    _false(selection["training_pool_max_is_continuum_bound"], "training max continuum promotion")
    _false(
        selection["post_holdout_retuning_under_same_experiment_identity_allowed"],
        "post-holdout retuning",
    )
    _false(
        selection["fresh_validation_samples_may_be_recycled_for_model_selection_and_still_called_held_out"],
        "held-out sample recycling",
    )

    observed = contract["observed_fresh_validation"]
    _equal(observed["validation_points_per_seed"], 4096, "fresh validation point count")
    _equal(observed["seeds"], [9175291, 9175292], "fresh validation seeds")
    _require(914027 not in observed["seeds"], "fresh seeds unexpectedly include registered CR001 seed")
    _equal(len(observed["rows"]), 2, "fresh validation row count")
    for row in observed["rows"]:
        _require(row["child_momentum_max"] < row["parent_momentum_max"], "sampled max did not improve")
        _require(row["child_momentum_L2"] > row["parent_momentum_L2"], "sampled L2 tradeoff disappeared")
        _require(row["child_momentum_max"] > 0.001, "child sampled max no longer fails registered gate")
        _require(row["child_momentum_L2"] > 0.001, "child sampled L2 no longer fails registered gate")
    _true(observed["sampled_momentum_max_improves_on_both_fresh_seeds"], "fresh sampled max result")
    _false(observed["sampled_momentum_L2_improves_on_both_fresh_seeds"], "fresh sampled L2 result")
    _false(observed["both_registered_momentum_thresholds_pass"], "joint momentum gate promotion")
    _false(observed["pde_validated"], "ST052 PDE promotion")
    _true(observed["fresh_seed_results_are_additional_generalization_evidence"], "fresh evidence scope")
    _false(
        observed["fresh_seed_results_are_the_preregistered_CR001_acceptance_sample"],
        "fresh-seed/registered-sample identity",
    )
    _false(observed["registered_CR001_seed_replayed_by_pr508"], "registered seed replay claim")
    _false(
        observed["same_validator_shape_or_point_count_may_replace_registered_seed"],
        "registered validation replacement",
    )
    _false(observed["unqualified_candidate_improved_claim_allowed"], "unqualified improvement claim")

    sampled = contract["sampled_evidence_scope"]
    for key in (
        "active_set_epigraph_is_acceptance_evidence",
        "full_training_pool_max_is_acceptance_evidence",
        "independent_random_or_edge_grid_max_is_continuum_upper_bound",
        "finite_core_structure_fractions_are_global_support_proofs",
        "DOP853_self_refinement_is_independent_integrator_validation",
        "trajectory_preservation_is_visual_correspondence",
        "vorticity_moment_change_is_public_target_match",
        "green_ci_is_scientific_acceptance",
    ):
        _false(sampled[key], f"sampled evidence promotion {key}")

    for inference in contract["forbidden_inferences"]:
        _false(inference["allowed"], f"forbidden inference {inference['from']} -> {inference['to']}")

    _audit_cr001(contract, constraints)

    delivery = contract["canonical_delivery"]
    _equal(delivery["candidate_family"], _CANONICAL_FAMILY, "canonical family")
    _equal(delivery["candidate_sha256"], _CANONICAL_SHA, "canonical SHA")
    _equal(delivery["velocity_api"], _CANONICAL_API, "canonical velocity API")
    _true(delivery["velocity_export_ready"], "canonical velocity export")
    _false(delivery["st052_replaces_canonical_candidate"], "ST052 canonical replacement")
    _false(delivery["st052_replaces_canonical_velocity_api"], "ST052 API replacement")
    _false(delivery["st052_fresh_validation_promotes_canonical_pde_state"], "canonical PDE promotion")

    _equal(project_status["candidate_family"], _CANONICAL_FAMILY, "project-status family")
    _equal(project_status["candidate_sha256"], _CANONICAL_SHA, "project-status SHA")
    _equal(project_status["velocity_api"], _CANONICAL_API, "project-status velocity API")
    _true(project_status["states"]["velocity_export_ready"], "project-status velocity export")
    for key in (
        "visualization_ready",
        "visual_correspondence_verified",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    ):
        _false(project_status["states"][key], f"project-status {key}")

    primary = delivery_contract["primary_deliverable"]
    _equal(primary["candidate_family"], _CANONICAL_FAMILY, "delivery-contract family")
    _equal(primary["candidate_sha256"], _CANONICAL_SHA, "delivery-contract SHA")
    _equal(primary["api"], _CANONICAL_API, "delivery-contract API")
    _equal(primary["api_signature"], "velocity(x,y,z,t)->[u,v,w]", "delivery API signature")
    _true(delivery_contract["claim_status"]["velocity_export_ready"], "delivery velocity export")
    for key in (
        "visualization_ready",
        "visual_correspondence_verified",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    ):
        _false(delivery_contract["claim_status"][key], f"delivery-contract {key}")

    future = contract["future_promotion_requirements"]
    for key in (
        "preserve_st052_identity_and_recipe_hashes",
        "no_retuning_after_acceptance_data_are_read_without_new_experiment_identity",
        "run_preregistered_seed_914027_4096_point_acceptance_contract_before_formal_CR001_promotion",
        "use_exact_registered_six_times_and_derivative_ladder",
        "require_momentum_max_and_L2_both_at_or_below_1e-3",
        "retain_divergence_and_other_hard_gates",
        "retain_restricted_forcing_and_no_free_f_equals_R_shortcut",
        "canonical_delivery_change_requires_explicit_versioned_identity_api_and_save_load_update",
        "visual_correspondence_requires_public_reference_comparison_and_independent_resolution_stability",
        "pde_failure_does_not_block_truth_bounded_candidate_save_load_or_visualization_work",
    ):
        _true(future[key], f"future requirement {key}")

    truth = contract["independent_truth_states"]
    _true(truth["velocity_export_ready"], "canonical delivery readiness")
    for key in (
        "visualization_ready",
        "visual_correspondence_verified",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    ):
        _false(truth[key], f"independent truth state {key}")

    return {
        "audit": "passed",
        "task_id": TASK_ID,
        "observed_pr": 508,
        "st052_candidate_sha256": _ST052_SHA,
        "sampled_momentum_max_improved_on_two_fresh_seeds": True,
        "sampled_momentum_L2_improved_on_two_fresh_seeds": False,
        "registered_cr001_seed_replayed": False,
        "pde_validated": False,
        "st052_is_canonical_candidate": False,
        "canonical_velocity_export_ready": True,
    }


if __name__ == "__main__":
    print(json.dumps(audit_st052_sampled_minimax_governance(), indent=2, sort_keys=True))
