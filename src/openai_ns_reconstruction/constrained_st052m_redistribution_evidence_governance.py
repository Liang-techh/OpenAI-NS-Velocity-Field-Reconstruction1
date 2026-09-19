"""Fail-closed CR002 audit for ST052-M redistribution evidence composition.

The ST052-M parent has its own sampled momentum receipt.  PR #528 applies a
frozen swirl-only redistribution and records Eulerian representation evidence;
PR #533 replays the same transform on material paths.  Those evidence lanes are
candidate-specific and must not be composed into a child PDE, visual, or
production claim before the transformed child receives a versioned identity and
fresh validation.

This module changes no velocity, coefficient, pressure, force, threshold, sample,
or canonical delivery state.
"""
from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
from typing import Any, Mapping

TASK_ID = "CR002-ST052M-REDISTRIBUTION-EVIDENCE-COMPOSITION-086"
_INTEGRATION_HEAD = "f18adf51953fc378c872904111f0d27deaaa7ccc"
_PARENT_HEAD = "b3b8bfdbe1077f9ec967d158602951997d81e17d"
_EULERIAN_HEAD = "779ffca71066e2864496d37de55a7aafc45d6f57"
_PATH_HEAD = "8e34021f08e3e778c5f5a07f5639c2e3cde3a5e1"
_PARENT_SHA = "e078e753fab38ebfa0284d28ba64d26cb8538b700a3afc7705668849c43e12da"
_PARENT_MODIFIER_SHA = "4d8b6e92f47151a7eed112985779a58fdf1188bd30c35c3a59263910629f1be9"
_CANONICAL_FAMILY = "eq45_supported_velocity_candidate_v1"
_CANONICAL_SHA = "2fdff812c22131d56eac1b7d6e455187ff3207500c6385aa9abf282a8e3d7b1d"
_CANONICAL_API = "openai_ns_reconstruction.eq45_supported_delivery:velocity"
_ALLOWED_CLASSES = {"user_requirement", "public_source_fact", "autonomous_design", "pending_unknown"}
_REGISTERED_TIMES = [0.25, 0.3125, 0.4375, 0.5625, 0.6875, 0.75]
_REGISTERED_STEPS = [0.02, 0.01, 0.005]


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _load_json(path: Path) -> dict[str, Any]:
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise ValueError(f"expected JSON object in {path}")
    return data


def load_contract() -> dict[str, Any]:
    return _load_json(_repo_root() / "configs" / "st052m_redistribution_evidence_composition_governance.json")


def load_parent_governance() -> dict[str, Any]:
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
    _equal(snap["parent_pr"], 508, "parent PR")
    _equal(snap["parent_pr_head"], _PARENT_HEAD, "parent head")
    _equal(snap["eulerian_transform_pr"], 528, "Eulerian PR")
    _equal(snap["eulerian_transform_head"], _EULERIAN_HEAD, "Eulerian head")
    _equal(snap["material_path_pr"], 533, "material-path PR")
    _equal(snap["material_path_head"], _PATH_HEAD, "material-path head")
    _equal(snap["eulerian_standard_run"], 35412489248, "Eulerian standard CI")
    _equal(snap["eulerian_dedicated_run"], 35412489260, "Eulerian dedicated CI")
    _equal(snap["eulerian_parent_replay_run"], 35412489236, "parent replay CI")
    _equal(snap["material_path_standard_run"], 35414322398, "path standard CI")
    _equal(snap["material_path_dedicated_run"], 35414322426, "path dedicated CI")
    _equal(snap["material_path_publication_run"], 35414322455, "path publication CI")
    _equal(
        snap["material_path_report_sha256"],
        "9ab3d96af9fbe15fcf147188a390b53e2c1af1df81804e2f9e8a150daddcecce",
        "material-path report SHA",
    )

    _equal(set(contract["allowed_source_classes"]), _ALLOWED_CLASSES, "allowed source classes")
    _equal(set(contract["source_classification"].values()), _ALLOWED_CLASSES, "source classes represented")

    ident = contract["parent_identity"]
    _equal(ident["candidate_id"], "ST052-M", "parent candidate id")
    _equal(ident["raw_candidate_sha256"], _PARENT_SHA, "parent candidate SHA")
    _equal(ident["modifier_sha256"], _PARENT_MODIFIER_SHA, "parent modifier SHA")
    _require(ident["experimental_materialized_child"] is True, "ST052-M must remain a materialized experiment")
    _require(ident["canonical_repository_candidate"] is False, "ST052-M cannot become canonical here")
    _require(ident["production_candidate_selected"] is False, "production candidate cannot be selected here")

    parent_ident = parent["candidate_identity"]
    _equal(parent_ident["candidate_id"], ident["candidate_id"], "parent governance candidate id")
    _equal(parent_ident["raw_candidate_sha256"], ident["raw_candidate_sha256"], "parent governance candidate SHA")
    _equal(parent_ident["modifier_sha256"], ident["modifier_sha256"], "parent governance modifier SHA")

    transform = contract["frozen_transform"]
    _close(transform["alpha"], 2.520520814687742, "alpha", atol=1e-15)
    _close(transform["kappa"], 0.05, "kappa", atol=0.0)
    _equal(transform["inner_window"], [0.30, 1.05], "inner window")
    _equal(transform["outer_window"], [0.95, 1.85], "outer window")
    _close(transform["exact_eulerian_normalization_receipt"], 1.0032534663681094, "Eulerian normalization", atol=2e-12)
    _close(transform["independent_material_path_normalization_replay"], 1.003253466368113, "path normalization", atol=2e-12)
    lower = 1.0 - transform["kappa"] * transform["alpha"]
    upper = 1.0 + transform["kappa"]
    _close(transform["multiplier_global_lower_bound_from_profile_range"], lower, "multiplier lower bound", atol=1e-15)
    _close(transform["multiplier_global_upper_bound_from_profile_range"], upper, "multiplier upper bound", atol=1e-15)
    _require(lower > 0.0 and transform["multiplier_positive_everywhere"] is True, "swirl multiplier must stay positive")
    _require(transform["identity_near_axis_below_r_0p30"] is True, "transform must be identity near axis")
    _equal(transform["new_velocity_degree_count"], 1, "velocity degree count")
    _require(transform["parameter_scan_performed_in_pr528"] is False, "PR528 must remain fixed-gain evidence")

    rep = contract["representation_audit"]
    for key in (
        "axisymmetric_parent_assumed_by_transform",
        "h_is_smooth_compact_radial_profile",
        "near_axis_parent_swirl_order_is_unchanged",
        "finite_positive_multiplier_preserves_parent_swirl_sign",
        "pointwise_multiplier_preserves_parent_zero_support",
        "common_positive_scale_preserves_parent_zero_support",
        "axisymmetric_swirl_only_multiplier_adds_no_theta_divergence_term",
        "common_constant_scale_preserves_divergence_identity",
    ):
        _require(rep[key] is True, f"representation invariant must remain true: {key}")
    for key in (
        "analytic_divergence_structure_preservation_is_full_pde_preservation",
        "pressure_can_be_inherited_without_rebuild",
        "restricted_force_fit_can_be_inherited_without_rebuild",
        "momentum_residual_can_be_inherited_from_parent",
        "reference_energy_normalization_is_navier_stokes_invariance",
    ):
        _require(rep[key] is False, f"forbidden representation inference became true: {key}")

    pde = contract["parent_pde_evidence"]
    _equal(pde["belongs_to_candidate"], "ST052-M", "PDE evidence owner")
    _equal(pde["fresh_validation_seeds"], [9175291, 9175292], "fresh validation seeds")
    _require(pde["sampled_max_improved"] is True, "parent max improvement must remain recorded")
    _require(pde["fixed_time_volume_L2_improved"] is False, "parent L2 regression must remain visible")
    _require(pde["formal_1e3_momentum_gates_pass"] is False, "formal momentum gate cannot pass")
    _require(pde["registered_CR001_seed_914027_replayed"] is False, "noncanonical fresh seeds cannot replace seed 914027")
    _require(pde["may_transfer_to_transformed_child"] is False, "parent PDE evidence cannot transfer to transformed child")
    _equal(parent["observed_fresh_validation"]["seeds"], pde["fresh_validation_seeds"], "parent evidence seeds")
    _require(parent["observed_fresh_validation"]["pde_validated"] is False, "parent must remain PDE-unvalidated")

    eul = contract["child_eulerian_evidence"]
    _close(eul["normalization"], transform["exact_eulerian_normalization_receipt"], "child normalization", atol=2e-12)
    _close(eul["divergence_fd_max"], 2.497e-9, "Eulerian FD divergence", atol=5e-13)
    _require(eul["clean_expression_capacity_transfer"] is True, "clean Eulerian transfer must remain recorded")
    _require(eul["full_registered_divergence_acceptance_replayed"] is False, "local FD screen cannot become formal divergence acceptance")
    _require(eul["full_momentum_replayed"] is False, "child momentum must remain unevaluated")
    _require(eul["pde_validated"] is False, "Eulerian screen cannot validate PDE")

    paths = contract["child_material_path_evidence"]
    _equal(paths["path_count"], 48, "path count")
    _equal(paths["paired_material_line_count"], 24, "pair count")
    _close(paths["mean_absolute_turns_relative_change"], 0.019601, "mean winding change", atol=5e-7)
    _close(paths["maximum_absolute_turns_relative_change"], 0.044819, "max winding change", atol=5e-7)
    _close(paths["radial_contraction_magnitude_relative_change"], 0.003274, "contraction change", atol=5e-7)
    _close(paths["mean_pair_axial_separation_change_relative"], 0.004353, "pair separation change", atol=5e-7)
    _equal(paths["inward_path_count_parent"], 48, "parent inward paths")
    _equal(paths["inward_path_count_child"], 48, "child inward paths")
    _require(paths["positive_material_path_transfer_under_frozen_protocol"] is True, "positive frozen path transfer must remain recorded")
    _require(paths["visual_correspondence_verified"] is False, "path replay is not visual correspondence")
    _require(paths["pde_evidence"] is False, "path replay is not PDE evidence")

    child = contract["child_delivery_identity"]
    for key in (
        "materialized_versioned_candidate_exists",
        "candidate_id_assigned",
        "candidate_sha256_assigned",
        "save_load_contract_exists",
        "unified_velocity_api_registered",
        "sampled_grid_export_identity_registered",
        "canonical_repository_candidate",
        "production_candidate_selected",
        "diagnostic_callable_is_delivery_identity",
    ):
        _require(child[key] is False, f"unmaterialized child state became true: {key}")
    _require(child["diagnostic_callable_may_be_used_for_truth_bounded_replay"] is True, "diagnostic replay must stay allowed")

    composition = contract["evidence_composition_rules"]
    for key in (
        "parent_pde_plus_child_eulerian_equals_child_pde",
        "parent_pde_plus_child_material_paths_equals_child_pde",
        "child_eulerian_plus_child_material_paths_equals_visual_correspondence",
        "clean_fd_divergence_plus_parent_momentum_equals_child_full_acceptance",
        "positive_material_path_transfer_selects_production_candidate",
        "materialization_priority_is_scientific_promotion",
        "green_ci_is_scientific_acceptance",
    ):
        _require(composition[key] is False, f"forbidden evidence composition became true: {key}")
    _require(composition["positive_material_path_transfer_supports_materialization_priority"] is True, "routing priority should remain explicit")

    cr = contract["canonical_cr001"]
    _close(cr["nu"], constraints["nu"], "nu", atol=0.0)
    _equal(cr["physical_domain"], constraints["domain"]["physical"], "physical domain")
    _equal(cr["evaluation_box"], constraints["domain"]["evaluation_box"], "evaluation box")
    _equal(cr["support"], constraints["domain"]["support"], "support")
    _equal(cr["time_interval"], constraints["domain"]["time_interval"], "time interval")
    _equal(cr["force_mode"], constraints["forcing"]["mode"], "force mode")
    _equal(cr["force_bounds"], constraints["forcing"]["parameters"], "force bounds")
    _close(cr["reference_energy"], constraints["nontriviality"]["reference_energy"], "reference energy", atol=0.0)
    _close(cr["reference_energy_abs_tolerance"], constraints["nontriviality"]["reference_energy_abs_tolerance"], "reference energy tolerance", atol=0.0)
    val = constraints["validation"]
    _equal(cr["validation_seed"], val["seed"], "validation seed")
    _equal(cr["held_out_points"], val["held_out_points"], "held-out points")
    _equal(cr["validation_times"], _REGISTERED_TIMES, "registered times")
    _equal(cr["validation_times"], val["times"], "constraint times")
    _equal(cr["derivative_steps"], _REGISTERED_STEPS, "registered derivative steps")
    _equal(cr["derivative_steps"], val["derivative_steps"], "constraint derivative steps")
    _close(cr["divergence_max"], val["thresholds"]["divergence_max"], "divergence max", atol=0.0)
    _close(cr["divergence_L2"], val["thresholds"]["divergence_L2"], "divergence L2", atol=0.0)
    _close(cr["pde_residual_max"], val["thresholds"]["pde_residual_max"], "momentum max", atol=0.0)
    _close(cr["pde_residual_L2"], val["thresholds"]["pde_residual_L2"], "momentum L2", atol=0.0)
    _require(cr["free_or_residual_defined_force_allowed"] is False, "free residual-defined forcing cannot be enabled")
    _require(cr["collapsed_velocity_success_allowed"] is False, "velocity collapse cannot be enabled")
    _require(cr["threshold_relaxation_under_same_experiment_allowed"] is False, "threshold relaxation cannot be enabled")

    canonical = contract["canonical_delivery"]
    _equal(canonical["candidate_family"], _CANONICAL_FAMILY, "canonical family")
    _equal(canonical["candidate_sha256"], _CANONICAL_SHA, "canonical SHA")
    _equal(canonical["velocity_api"], _CANONICAL_API, "canonical API")
    _equal(project_status["candidate_family"], _CANONICAL_FAMILY, "project-status family")
    _equal(project_status["candidate_sha256"], _CANONICAL_SHA, "project-status SHA")
    _equal(project_status["velocity_api"], _CANONICAL_API, "project-status API")
    primary = delivery["primary_deliverable"]
    _equal(primary["candidate_family"], _CANONICAL_FAMILY, "delivery family")
    _equal(primary["candidate_sha256"], _CANONICAL_SHA, "delivery SHA")
    _equal(primary["api"], _CANONICAL_API, "delivery API")
    _require(canonical["velocity_export_ready"] is True, "canonical velocity must stay export ready")
    _require(canonical["transformed_st052m_replaces_canonical_candidate"] is False, "transformed child cannot replace canonical candidate")
    _require(canonical["transformed_st052m_replaces_canonical_velocity_api"] is False, "transformed child cannot replace canonical API")

    truth = contract["independent_truth_states"]
    _require(truth["velocity_export_ready"] is True, "canonical velocity_export_ready must remain true")
    for key in (
        "transformed_child_velocity_export_ready",
        "visualization_ready",
        "visual_correspondence_verified",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    ):
        _require(truth[key] is False, f"independent truth state cannot be promoted: {key}")
    _equal(project_status["states"]["velocity_export_ready"], truth["velocity_export_ready"], "project velocity-export state")
    for key in ("visualization_ready", "visual_correspondence_verified", "pde_validated", "paper_exact", "openai_field_identified", "blowup_proved"):
        _equal(project_status["states"][key], truth[key], f"project truth state {key}")

    future = contract["future_materialization_requirements"]
    for key, required in future.items():
        _require(required is True, f"future materialization requirement disabled: {key}")

    for row in contract["forbidden_inferences"]:
        _require(row["allowed"] is False, f"forbidden inference enabled: {row['from']} -> {row['to']}")

    return {
        "task_id": TASK_ID,
        "status": "pass",
        "parent_candidate": "ST052-M",
        "transformed_child_materialized": False,
        "parent_pde_receipt_transfers": False,
        "positive_material_path_transfer": True,
        "materialization_priority_supported": True,
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


def audited_copy_with_mutation(path: tuple[str, ...], value: Any) -> dict[str, Any]:
    """Test helper: mutate one contract leaf and run the same fail-closed audit."""
    contract = deepcopy(load_contract())
    cursor: dict[str, Any] = contract
    for key in path[:-1]:
        cursor = cursor[key]
    cursor[path[-1]] = value
    return audit_contract(
        contract,
        load_constraints(),
        load_parent_governance(),
        load_project_status(),
        load_velocity_delivery_contract(),
    )


if __name__ == "__main__":
    print(json.dumps(audit_current_repository(), indent=2, sort_keys=True))
