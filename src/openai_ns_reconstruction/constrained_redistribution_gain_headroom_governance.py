"""Fail-closed CR002 governance for redistribution gain-headroom evidence.

A higher preregistered gain can establish representation capacity without
inheriting material-path, sampled-grid, pressure/forcing, PDE, visual, or
identity evidence obtained at a lower gain.  This module changes no numerical
velocity, candidate, route, or scientific threshold.
"""
from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
from typing import Any, Mapping

TASK_ID = "CR002-REDISTRIBUTION-GAIN-HEADROOM-080"
_ALLOWED_CLASSES = {
    "user_requirement",
    "public_source_fact",
    "autonomous_design",
    "pending_unknown",
}
_INTEGRATION_HEAD = "c0f20e0712f4347698485a2de1e3cb980b5ecfb1"
_HEADROOM_PR = 480
_HEADROOM_HEAD = "187157d377b14bae56415648ac100342f2b8bbb0"
_PATH_PR = 471
_PATH_HEAD = "ad028f512346ad875e6e865018c5d234b1526ea1"
_CANONICAL_SHA = "2fdff812c22131d56eac1b7d6e455187ff3207500c6385aa9abf282a8e3d7b1d"
_CANONICAL_API = "openai_ns_reconstruction.eq45_supported_delivery:velocity"
_LIVE_ROUTE = "materialize_integrated_axial_cap_poloidal_child_then_fresh_validate"


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _load_json(path: Path) -> dict[str, Any]:
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise ValueError(f"expected object in {path}")
    return data


def load_contract() -> dict[str, Any]:
    return _load_json(_repo_root() / "configs" / "redistribution_gain_headroom_governance.json")


def load_cross_backbone_contract() -> dict[str, Any]:
    return _load_json(_repo_root() / "configs" / "cross_backbone_swirl_redistribution_governance.json")


def load_delivery_state_contract() -> dict[str, Any]:
    return _load_json(_repo_root() / "configs" / "delivery_state_contract.json")


def load_constraints() -> dict[str, Any]:
    return _load_json(_repo_root() / "configs" / "constraints.json")


def load_project_status() -> dict[str, Any]:
    return _load_json(_repo_root() / "project_status.json")


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _require_equal(actual: Any, expected: Any, label: str) -> None:
    _require(actual == expected, f"{label} drifted: {actual!r} != {expected!r}")


def audit_redistribution_gain_headroom(
    contract: Mapping[str, Any] | None = None,
    cross_backbone_contract: Mapping[str, Any] | None = None,
    delivery_state_contract: Mapping[str, Any] | None = None,
    constraints: Mapping[str, Any] | None = None,
    project_status: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Audit gain-scoped evidence and unchanged repository contracts."""

    contract = deepcopy(dict(load_contract() if contract is None else contract))
    cross = deepcopy(
        dict(load_cross_backbone_contract() if cross_backbone_contract is None else cross_backbone_contract)
    )
    delivery = deepcopy(
        dict(load_delivery_state_contract() if delivery_state_contract is None else delivery_state_contract)
    )
    constraints = deepcopy(dict(load_constraints() if constraints is None else constraints))
    project = deepcopy(dict(load_project_status() if project_status is None else project_status))

    _require_equal(contract.get("schema_version"), 1, "contract schema")
    _require_equal(contract.get("task_id"), TASK_ID, "task id")
    _require_equal(contract.get("status"), "governance_only_no_candidate_promotion", "governance status")

    snapshot = contract.get("snapshot", {})
    _require_equal(snapshot.get("active_integration_branch"), "codex/cr001-constraints", "integration branch")
    _require_equal(snapshot.get("active_integration_head"), _INTEGRATION_HEAD, "integration head")
    _require_equal(snapshot.get("observed_headroom_pr"), _HEADROOM_PR, "headroom PR")
    _require_equal(snapshot.get("observed_headroom_head"), _HEADROOM_HEAD, "headroom head")
    _require_equal(snapshot.get("observed_material_path_pr"), _PATH_PR, "material-path PR")
    _require_equal(snapshot.get("observed_material_path_head"), _PATH_HEAD, "material-path head")
    _require_equal(snapshot.get("observed_sampled_grid_pr"), 476, "sampled-grid PR")

    _require_equal(set(contract.get("allowed_source_classes", ())), _ALLOWED_CLASSES, "source-class vocabulary")
    _require_equal(
        set(delivery.get("classification_vocabulary", {})),
        _ALLOWED_CLASSES,
        "delivery-state source vocabulary",
    )
    source = contract.get("source_classification", {})
    _require(source and set(source.values()) <= _ALLOWED_CLASSES, "invalid source classification")
    _require_equal(source.get("callable_saveable_velocity_delivery"), "user_requirement", "callable delivery class")
    _require_equal(
        source.get("public_openai_velocity_visualization_observables"),
        "public_source_fact",
        "public visualization class",
    )
    for key in (
        "st051b_parent",
        "frozen_redistribution_profile_alpha_windows",
        "gain_grid_and_headroom_thresholds",
        "common_reference_energy_normalization",
        "material_path_protocol",
        "sampled_grid_protocol",
        "cr001_validation_protocol",
    ):
        _require_equal(source.get(key), "autonomous_design", f"source class {key}")
    for key in (
        "openai_hidden_swirl_redistribution_gain",
        "openai_hidden_numerical_velocity",
        "openai_hidden_frame_time_camera_mapping",
    ):
        _require_equal(source.get(key), "pending_unknown", f"source class {key}")

    identity = contract.get("fixed_transform_identity", {})
    cross_st050 = cross.get("evidence_snapshots", {}).get("st050rc", {})
    _require_equal(identity.get("parent"), "ST051-B", "redistribution parent")
    _require_equal(identity.get("parent_pr"), 460, "parent PR")
    _require_equal(identity.get("parent_head"), "4b784f1b8457af2ead49295631d834d4e882000b", "parent head")
    _require_equal(identity.get("alpha"), cross_st050.get("balance_alpha"), "frozen alpha identity")
    _require_equal(identity.get("inner_window"), cross_st050.get("inner_window"), "inner window identity")
    _require_equal(identity.get("outer_window"), cross_st050.get("outer_window"), "outer window identity")
    _require_equal(identity.get("profile_rebalanced_on_st051b"), False, "ST051-B rebalance status")
    _require_equal(identity.get("gain_is_part_of_materialized_child_identity"), True, "gain identity rule")
    _require_equal(identity.get("changing_gain_requires_new_child_identity"), True, "new child identity rule")
    _require_equal(identity.get("changing_gain_requires_new_candidate_sha256"), True, "new candidate SHA rule")

    by_gain = contract.get("evidence_by_gain", {})
    low = by_gain.get("gain_0_025", {})
    high = by_gain.get("gain_0_05", {})
    _require_equal(low.get("gain"), 0.025, "lower evidence gain")
    _require_equal(low.get("material_path_evidence_exists"), True, "lower-gain material-path evidence")
    _require_equal(low.get("material_path_pr"), _PATH_PR, "lower-gain material-path PR")
    _require_equal(low.get("sampled_grid_export_pr"), 476, "lower-gain sampled-grid PR")
    _require_equal(low.get("sampled_grid_evidence_is_scoped_to_this_gain"), True, "sampled-grid gain scope")
    _require_equal(low.get("production_gain_selected"), False, "lower-gain production status")
    _require_equal(low.get("visual_correspondence_verified"), False, "lower-gain visual truth")
    _require_equal(low.get("pde_validated"), False, "lower-gain PDE truth")

    _require_equal(high.get("gain"), 0.05, "higher evidence gain")
    _require_equal(high.get("headroom_pr"), _HEADROOM_PR, "higher-gain headroom PR")
    _require_equal(high.get("preregistered_gain_grid"), [0.0, 0.015, 0.025, 0.035, 0.05], "headroom gain grid")
    _require_equal(high.get("first_clean_higher_gain_crossing"), True, "headroom crossing status")
    _require_equal(high.get("common_reference_energy_normalization"), 1.002495582074808, "headroom normalization")
    _require_equal(high.get("cartesian_fd_divergence_max"), 1.5592e-9, "headroom divergence diagnostic")
    for key in (
        "material_path_evidence_exists",
        "sampled_grid_export_evidence_exists",
        "pressure_forcing_refit_exists",
        "fresh_full_momentum_receipt_exists",
        "production_gain_selected",
        "visual_correspondence_verified",
        "pde_validated",
    ):
        _require_equal(high.get(key), False, f"higher-gain truth {key}")

    semantics = contract.get("gain_evidence_semantics", {})
    _require_equal(semantics.get("clean_headroom_crossing_is_capacity_evidence"), True, "capacity-evidence rule")
    for key in (
        "clean_headroom_crossing_selects_production_gain",
        "largest_preregistered_clean_gain_is_production_upper_bound",
        "same_profile_at_different_gain_is_same_materialized_child",
        "evidence_may_transfer_between_gain_values_without_fresh_replay",
        "lower_gain_material_paths_may_be_linearly_extrapolated_to_higher_gain",
        "lower_gain_material_paths_may_be_inherited_by_higher_gain",
        "lower_gain_sampled_grid_may_be_relabelled_as_higher_gain",
        "parent_or_lower_gain_pressure_forcing_may_be_inherited",
        "parent_or_lower_gain_pde_receipt_may_be_inherited",
        "sampled_q90_q99_stability_is_continuum_morphology_certificate",
        "eulerian_proxy_improvement_is_visual_correspondence_certificate",
        "eulerian_proxy_improvement_is_pde_evidence",
        "clean_crossing_recovers_openai_hidden_gain",
    ):
        _require_equal(semantics.get(key), False, f"gain evidence rule {key}")

    cross_promotion = cross.get("promotion_contract", {})
    _require_equal(
        cross_promotion.get("diagnostic_gain_may_become_production_value_without_explicit_selection"),
        False,
        "existing diagnostic-to-production guard",
    )
    _require_equal(
        cross_promotion.get("parent_material_paths_transfer_through_redistribution"),
        False,
        "existing material-path transfer guard",
    )
    _require_equal(
        cross_promotion.get("parent_pde_receipts_transfer_through_redistribution"),
        False,
        "existing PDE transfer guard",
    )

    requirements = contract.get("materialization_requirements_for_gain_0_05", {})
    for key in (
        "explicit_production_bound_and_value_registration",
        "new_representation_identity",
        "new_candidate_sha256",
        "deterministic_save_load_replay",
        "fresh_reference_energy_check",
        "fresh_axis_near_axis_support_core_checks",
        "fresh_cartesian_divergence_check",
        "fresh_material_paths_required_before_material_path_claim",
        "fresh_sampled_grid_required_before_grid_export_claim",
        "fresh_compatible_pressure_and_restricted_force_required_before_pde_claim",
        "fresh_unused_held_out_full_momentum_and_divergence_required_before_pde_promotion",
        "no_residual_defined_free_force",
        "no_zero_or_amplitude_collapse_shortcut",
    ):
        _require_equal(requirements.get(key), True, f"materialization requirement {key}")

    domain = constraints.get("domain", {})
    forcing = constraints.get("forcing", {})
    nontriviality = constraints.get("nontriviality", {})
    validation = constraints.get("validation", {})
    thresholds = validation.get("thresholds", {})
    cr001 = contract.get("canonical_cr001", {})
    expected_pairs = (
        (cr001.get("nu"), constraints.get("nu"), "nu"),
        (cr001.get("physical_domain"), domain.get("physical"), "physical domain"),
        (cr001.get("evaluation_box"), domain.get("evaluation_box"), "evaluation box"),
        (cr001.get("support"), domain.get("support"), "support"),
        (cr001.get("time_interval"), domain.get("time_interval"), "time interval"),
        (cr001.get("force_mode"), forcing.get("mode"), "force mode"),
        (cr001.get("force_bounds"), forcing.get("parameters"), "force bounds"),
        (cr001.get("reference_energy"), nontriviality.get("reference_energy"), "reference energy"),
        (cr001.get("reference_energy_abs_tolerance"), nontriviality.get("reference_energy_abs_tolerance"), "energy tolerance"),
        (cr001.get("validation_seed"), validation.get("seed"), "validation seed"),
        (cr001.get("held_out_points"), validation.get("held_out_points"), "held-out points"),
        (cr001.get("validation_times"), validation.get("times"), "validation times"),
        (cr001.get("derivative_steps"), validation.get("derivative_steps"), "derivative ladder"),
        (cr001.get("divergence_max"), thresholds.get("divergence_max"), "divergence max"),
        (cr001.get("divergence_L2"), thresholds.get("divergence_L2"), "divergence L2"),
        (cr001.get("pde_residual_max"), thresholds.get("pde_residual_max"), "PDE max"),
        (cr001.get("pde_residual_L2"), thresholds.get("pde_residual_L2"), "PDE L2"),
    )
    for expected, actual, label in expected_pairs:
        _require_equal(actual, expected, label)
    _require("No residual-dependent basis" in forcing.get("restriction", ""), "free-force restriction weakened")

    canonical = contract.get("canonical_delivery_and_truth", {})
    states = project.get("states", {})
    _require_equal(canonical.get("candidate_family"), project.get("candidate_family"), "canonical family")
    _require_equal(canonical.get("candidate_sha256"), _CANONICAL_SHA, "canonical SHA")
    _require_equal(canonical.get("candidate_sha256"), project.get("candidate_sha256"), "project candidate SHA")
    _require_equal(canonical.get("velocity_api"), _CANONICAL_API, "canonical velocity API")
    _require_equal(canonical.get("velocity_api"), project.get("velocity_api"), "project velocity API")
    _require_equal(canonical.get("active_scientific_route"), _LIVE_ROUTE, "live scientific route")
    _require_equal(canonical.get("active_scientific_route"), project.get("active_scientific_route"), "project scientific route")
    _require_equal(canonical.get("velocity_export_ready"), True, "canonical export truth")
    _require_equal(states.get("velocity_export_ready"), True, "project export truth")
    for key in (
        "visualization_ready",
        "visual_correspondence_verified",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    ):
        _require_equal(canonical.get(key), False, f"contract truth {key}")
        _require_equal(states.get(key), False, f"project truth {key}")

    mutation = contract.get("mutation_scope", {})
    for key in (
        "velocity_changed",
        "redistribution_gain_changed",
        "candidate_changed",
        "pressure_changed",
        "forcing_changed",
        "optimizer_changed",
        "validation_sample_changed",
        "norm_definition_changed",
        "threshold_changed",
        "active_route_changed",
        "canonical_velocity_api_changed",
    ):
        _require_equal(mutation.get(key), False, f"mutation scope {key}")

    return {
        "task_id": TASK_ID,
        "status": "governance_pass",
        "headroom_pr": _HEADROOM_PR,
        "lower_gain_with_material_paths": 0.025,
        "higher_gain_capacity_crossing": 0.05,
        "higher_gain_production_selected": False,
        "higher_gain_material_path_evidence_exists": False,
        "lower_gain_evidence_transfer_allowed": False,
        "higher_gain_requires_new_identity_and_sha": True,
        "canonical_velocity_export_ready": True,
        "canonical_thresholds_unchanged": True,
        "visual_correspondence_verified": False,
        "pde_validated": False,
    }


def main() -> None:
    print(json.dumps(audit_redistribution_gain_headroom(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
