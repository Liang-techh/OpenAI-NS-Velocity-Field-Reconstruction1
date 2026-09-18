"""Fail-closed CR002 governance for ST048-S material-path evidence.

This module changes no velocity values. It keeps Eulerian fixed-probe proxies,
Lagrangian pathline diagnostics, visual correspondence, candidate identity, and
Navier-Stokes acceptance as separate evidence layers.
"""
from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
from typing import Any, Mapping

TASK_ID = "CR002-ST048S-MATERIAL-PATH-EVIDENCE-GOVERNANCE-066"
MATERIAL_PATH_PR = 435
MATERIAL_PATH_HEAD = "6699a698c9fb2f0a0986e7aefc5fe0ee431cd5ad"
SWIRL_PARENT_PR = 428
SWIRL_PARENT_HEAD = "b2b4888ece83a9f862435ce1adb847c94ad3ad7d"
RADIAL_RING_PR = 437
RADIAL_RING_HEAD = "6d485adb3b9c65e9719bd247fd5cc1ddfa9afb20"
LIVE_ROUTE = "materialize_integrated_axial_cap_poloidal_child_then_fresh_validate"
CANONICAL_FAMILY = "eq45_supported_velocity_candidate_v1"
CANONICAL_API = "openai_ns_reconstruction.eq45_supported_delivery:velocity"

_ALLOWED_CLASSES = {
    "user_requirement",
    "public_source_fact",
    "autonomous_design",
    "pending_unknown",
}
_FALSE_PROJECT_STATES = (
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
    return _load_json(_repo_root() / "configs" / "st048s_material_path_evidence_governance.json")


def load_constraints() -> dict[str, Any]:
    return _load_json(_repo_root() / "configs" / "constraints.json")


def load_project_status() -> dict[str, Any]:
    return _load_json(_repo_root() / "project_status.json")


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _require_equal(actual: Any, expected: Any, label: str) -> None:
    _require(actual == expected, f"{label} drifted: {actual!r} != {expected!r}")


def _require_close(actual: Any, expected: float, label: str, tol: float = 1e-12) -> None:
    _require(abs(float(actual) - expected) <= tol, f"{label} drifted: {actual!r} != {expected!r}")


def _audit_cr001(snapshot: Mapping[str, Any], constraints: Mapping[str, Any]) -> None:
    domain = constraints.get("domain", {})
    forcing = constraints.get("forcing", {})
    nontriviality = constraints.get("nontriviality", {})
    validation = constraints.get("validation", {})
    thresholds = validation.get("thresholds", {})

    _require_close(snapshot.get("nu"), constraints.get("nu"), "nu")
    _require_equal(snapshot.get("physical_domain"), domain.get("physical"), "physical domain")
    _require_equal(snapshot.get("evaluation_box"), domain.get("evaluation_box"), "evaluation box")
    _require_equal(snapshot.get("support"), domain.get("support"), "support")
    _require_equal(snapshot.get("time_interval"), domain.get("time_interval"), "time interval")
    _require_equal(snapshot.get("force_mode"), forcing.get("mode"), "force mode")
    _require_equal(snapshot.get("force_bounds"), forcing.get("parameters"), "force bounds")
    _require_close(snapshot.get("reference_energy"), nontriviality.get("reference_energy"), "reference energy")
    _require_close(
        snapshot.get("reference_energy_abs_tolerance"),
        nontriviality.get("reference_energy_abs_tolerance"),
        "reference energy tolerance",
    )
    _require_equal(snapshot.get("validation_seed"), validation.get("seed"), "validation seed")
    _require_equal(snapshot.get("held_out_points"), validation.get("held_out_points"), "held-out points")
    _require_equal(snapshot.get("validation_times"), validation.get("times"), "validation times")
    _require_equal(snapshot.get("derivative_steps"), validation.get("derivative_steps"), "derivative steps")
    for key in ("divergence_max", "divergence_L2", "pde_residual_max", "pde_residual_L2"):
        _require_close(snapshot.get(key), thresholds.get(key), key)

    restriction = str(forcing.get("restriction", ""))
    _require("No residual-dependent basis" in restriction, "residual-dependent forcing ban missing")
    _require("pointwise free force" in restriction, "pointwise free force ban missing")


def audit_material_path_evidence_governance(
    scope: Mapping[str, Any] | None = None,
    constraints: Mapping[str, Any] | None = None,
    project_status: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Audit material-path evidence without selecting or promoting a child."""

    scope = deepcopy(dict(load_scope() if scope is None else scope))
    constraints = deepcopy(dict(load_constraints() if constraints is None else constraints))
    project_status = deepcopy(dict(load_project_status() if project_status is None else project_status))

    _require_equal(scope.get("schema_version"), 1, "schema version")
    _require_equal(scope.get("task_id"), TASK_ID, "task id")

    governed = scope.get("governed_evidence", {})
    _require_equal(governed.get("material_path_pr"), MATERIAL_PATH_PR, "material-path PR")
    _require_equal(governed.get("material_path_head_sha"), MATERIAL_PATH_HEAD, "material-path head")
    _require_equal(
        governed.get("material_path_status"),
        "open_unintegrated_diagnostic_evidence",
        "material-path status",
    )
    _require_equal(governed.get("swirl_gain_parent_pr"), SWIRL_PARENT_PR, "swirl parent PR")
    _require_equal(governed.get("swirl_gain_parent_head_sha"), SWIRL_PARENT_HEAD, "swirl parent head")
    _require_equal(governed.get("radial_ring_pr"), RADIAL_RING_PR, "radial-ring PR")
    _require_equal(governed.get("radial_ring_head_sha"), RADIAL_RING_HEAD, "radial-ring head")
    _require_equal(
        governed.get("radial_ring_status"),
        "open_unintegrated_proxy_screen_evidence",
        "radial-ring status",
    )
    _require_equal(governed.get("candidate_materialized"), False, "candidate materialization")

    allowed = scope.get("allowed_source_classes", [])
    _require_equal(set(allowed), _ALLOWED_CLASSES, "allowed source classes")
    source = scope.get("source_classification", {})
    _require(bool(source), "source classification missing")
    for label, classification in source.items():
        _require(classification in _ALLOWED_CLASSES, f"invalid source class for {label}: {classification!r}")
    _require_equal(source.get("callable_saveable_velocity_delivery"), "user_requirement", "delivery source class")
    _require_equal(
        source.get("public_openai_velocity_field_visualization_as_observable_target"),
        "public_source_fact",
        "public visualization class",
    )
    for key in (
        "st048s_parent_candidate",
        "temporal_piola_schedule",
        "kappa_0p05_diagnostic_swirl_gain",
        "material_path_seed_set_and_integrator_contract",
        "material_path_winding_contraction_and_pair_separation_metrics",
        "radial_ring_profile_and_gain_grid",
        "fixed_probe_angular_rate_proxy",
    ):
        _require_equal(source.get(key), "autonomous_design", f"{key} source class")
    for key in (
        "openai_hidden_particle_seeds_or_material_trajectories",
        "openai_hidden_camera_frame_time_alignment",
        "openai_visualization_contains_tracked_material_particles_or_pathlines",
    ):
        _require_equal(source.get(key), "pending_unknown", f"{key} source class")

    paths = scope.get("material_path_contract", {})
    _require_equal(paths.get("definition"), "dX/dt = velocity(X(t), t)", "pathline definition")
    _require_equal(paths.get("time_interval"), [0.25, 0.75], "pathline time interval")
    _require_equal(paths.get("initial_radii"), [0.6, 0.9, 1.2], "pathline radii")
    _require_equal(paths.get("initial_z"), [-0.3, 0.3], "pathline z seeds")
    _require_equal(paths.get("azimuth_count"), 8, "azimuth count")
    _require_equal(paths.get("path_count"), 48, "path count")
    _require_equal(paths.get("axial_pair_count"), 24, "axial pair count")
    _require_equal(paths.get("output_time_samples"), 33, "output sample count")
    _require_equal(paths.get("integrator"), "scipy.solve_ivp:DOP853", "integrator")
    _require_close(paths.get("rtol"), 1e-9, "pathline rtol")
    _require_close(paths.get("atol"), 1e-11, "pathline atol")
    _require_close(paths.get("max_step"), 0.01, "pathline max step")
    for key in (
        "same_run_comparisons_require_same_initial_seed_set",
        "same_run_comparisons_require_same_time_interval",
        "same_run_comparisons_require_same_integrator_and_tolerances",
        "same_run_comparisons_require_immutable_field_lineage",
        "changing_any_contract_item_creates_a_new_diagnostic_protocol",
    ):
        _require_equal(paths.get(key), True, key)
    _require_equal(paths.get("successful_ode_integration_is_pde_validation"), False, "ODE/PDE separation")

    semantics = scope.get("observable_semantics", {})
    _require("frozen time" in str(semantics.get("instantaneous_streamline_definition", "")), "streamline definition weakened")
    _require("dX/dt=velocity" in str(semantics.get("material_path_definition", "")), "material-path definition weakened")
    _require("fixed spatial probes" in str(semantics.get("fixed_probe_eulerian_proxy_definition", "")), "proxy definition weakened")
    for key in (
        "streamlines_and_material_paths_are_interchangeable_for_time_dependent_fields",
        "fixed_probe_eulerian_proxy_is_material_path_evidence",
        "eulerian_vorticity_morphology_is_material_path_evidence",
        "material_path_winding_is_visual_correspondence",
        "material_path_winding_is_pde_evidence",
        "material_path_contraction_is_pde_evidence",
        "material_path_pair_separation_is_blowup_evidence",
        "positive_common_energy_rescaling_preserves_fixed_physical_time_material_paths_in_general",
        "coordinate_warp_preserves_material_paths_without_fresh_integration",
    ):
        _require_equal(semantics.get(key), False, key)

    receipt = scope.get("pr435_diagnostic_receipt", {})
    _require_close(receipt.get("kappa"), 0.05, "PR435 diagnostic kappa")
    _require_close(receipt.get("common_reference_energy_scale"), 0.972662771822918, "PR435 energy scale")
    _require_close(receipt.get("mean_absolute_turns_relative_change"), 0.0189884, "mean winding change", 1e-9)
    _require_close(receipt.get("maximum_absolute_turns_relative_change"), 0.0196838, "max winding change", 1e-9)
    _require_close(receipt.get("radial_contraction_magnitude_relative_change"), -0.0259322, "contraction change", 1e-9)
    _require_close(receipt.get("mean_axial_pair_separation_relative_change"), -0.037589, "pair separation change", 1e-9)
    _require_equal(receipt.get("inward_paths"), 48, "inward path count")
    _require_equal(receipt.get("total_paths"), 48, "total path count")
    _require_equal(receipt.get("axial_grow_pairs"), 16, "growing axial pair count")
    _require_equal(receipt.get("axial_shrink_pairs"), 8, "shrinking axial pair count")
    _require_equal(receipt.get("interpretation"), "diagnostic_tradeoff_only", "PR435 interpretation")
    for key in ("production_kappa_selected", "visual_correspondence_verified", "pde_validated"):
        _require_equal(receipt.get(key), False, f"PR435 {key}")

    ring = scope.get("pr437_proxy_receipt", {})
    _require_equal(ring.get("smallest_preregistered_clean_radial_swirl_replacement_crossing"), None, "ring crossing")
    _require_equal(ring.get("material_path_replay_run"), False, "ring material-path replay")
    for key in (
        "production_ring_gain_selected",
        "no_proxy_crossing_proves_no_material_path_improvement",
        "proxy_failure_promotes_global_kappa_to_production",
        "independent_rank_four_sensitivity_proves_pde_leverage",
    ):
        _require_equal(ring.get(key), False, f"ring {key}")

    promotion = scope.get("promotion_contract", {})
    for key in (
        "material_path_metrics_may_select_production_parameters_without_other_governed_gates",
        "proxy_metrics_may_select_production_parameters_without_material_path_replay_when_pathlines_are_the_claim",
        "material_path_or_proxy_improvement_may_replace_canonical_velocity",
        "material_path_or_proxy_improvement_may_replace_live_scientific_route",
    ):
        _require_equal(promotion.get(key), False, key)
    for key in (
        "new_materialized_child_requires_new_representation_identity",
        "new_materialized_child_requires_new_candidate_sha",
        "save_load_replay_required",
        "fresh_material_path_replay_required_if_material_path_claim_is_made",
        "fresh_pressure_restricted_force_and_full_pde_validation_required_for_pde_claim",
        "unused_acceptance_data_required_for_formal_pde_promotion",
    ):
        _require_equal(promotion.get(key), True, key)

    _audit_cr001(scope.get("canonical_cr001_snapshot", {}), constraints)

    routing = scope.get("routing_scope", {})
    _require_equal(routing.get("live_route_must_remain"), LIVE_ROUTE, "scope live route")
    _require_equal(routing.get("canonical_velocity_family_must_remain"), CANONICAL_FAMILY, "scope family")
    _require_equal(routing.get("canonical_velocity_api_must_remain"), CANONICAL_API, "scope API")
    _require_equal(routing.get("velocity_export_ready_may_remain_true"), True, "velocity export independence")
    _require_equal(project_status.get("active_scientific_route"), LIVE_ROUTE, "project live route")
    _require_equal(project_status.get("candidate_family"), CANONICAL_FAMILY, "project candidate family")
    _require_equal(project_status.get("velocity_api"), CANONICAL_API, "project velocity API")

    states = project_status.get("states", {})
    _require_equal(states.get("velocity_export_ready"), True, "existing velocity export state")
    for key in _FALSE_PROJECT_STATES:
        _require_equal(states.get(key), False, f"project state {key}")

    truth = scope.get("truth_boundary", {})
    _require_equal(truth.get("material_path_protocol_registered"), True, "material-path protocol registration")
    for key in (
        "material_path_child_materialized",
        "production_kappa_selected",
        "production_ring_gain_selected",
        "visualization_ready",
        "visual_correspondence_verified",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    ):
        _require_equal(truth.get(key), False, f"truth boundary {key}")

    return {
        "task_id": TASK_ID,
        "status": "pass",
        "material_path_protocol_registered": True,
        "canonical_velocity_unchanged": True,
        "live_route_unchanged": True,
        "visual_correspondence_verified": False,
        "pde_validated": False,
    }


def main() -> None:
    print(json.dumps(audit_material_path_evidence_governance(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
