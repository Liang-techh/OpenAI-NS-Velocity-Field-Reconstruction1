"""Fail-closed CR002 governance for ST048-S swirl redistribution evidence.

This module changes no velocity values. It separates first-order reference-time
energy orthogonality, finite-gain normalization, proxy/pathline evidence, and
Navier--Stokes acceptance.
"""
from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
from typing import Any, Mapping

TASK_ID = "CR002-ST048S-ENERGY-NEUTRAL-SWIRL-REDISTRIBUTION-GOVERNANCE-067"
REDISTRIBUTION_PR = 447
REDISTRIBUTION_HEAD = "527e20e7eaa506d144c725d2c90c468da6aa5d34"
RADIUS_ATTRIBUTION_PR = 445
RADIUS_ATTRIBUTION_HEAD = "8d04a0dae7d2f5ff0bed60ce7b1c1aaaba58c46f"
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
    return _load_json(_repo_root() / "configs" / "st048s_energy_neutral_swirl_redistribution_governance.json")


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


def audit_energy_neutral_swirl_redistribution_governance(
    scope: Mapping[str, Any] | None = None,
    constraints: Mapping[str, Any] | None = None,
    project_status: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Audit the open redistribution capacity evidence without promoting it."""

    scope = deepcopy(dict(load_scope() if scope is None else scope))
    constraints = deepcopy(dict(load_constraints() if constraints is None else constraints))
    project_status = deepcopy(dict(load_project_status() if project_status is None else project_status))

    _require_equal(scope.get("schema_version"), 1, "schema version")
    _require_equal(scope.get("task_id"), TASK_ID, "task id")

    governed = scope.get("governed_evidence", {})
    _require_equal(governed.get("redistribution_pr"), REDISTRIBUTION_PR, "redistribution PR")
    _require_equal(governed.get("redistribution_head_sha"), REDISTRIBUTION_HEAD, "redistribution head")
    _require_equal(
        governed.get("redistribution_status"),
        "open_unintegrated_capacity_evidence",
        "redistribution status",
    )
    _require_equal(governed.get("radius_attribution_pr"), RADIUS_ATTRIBUTION_PR, "radius attribution PR")
    _require_equal(governed.get("radius_attribution_head_sha"), RADIUS_ATTRIBUTION_HEAD, "radius attribution head")
    _require_equal(
        governed.get("radius_attribution_status"),
        "open_unintegrated_diagnostic_evidence",
        "radius attribution status",
    )
    _require_equal(governed.get("material_path_governance_pr"), 444, "material-path governance PR")
    _require_equal(governed.get("st048s_parent_status"), "derived_challenger_not_canonical", "ST048-S parent status")
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
        "st006_radius_deficit_attribution_protocol",
        "inner_and_outer_bump_profiles",
        "reference_time_balance_coefficient_alpha",
        "redistribution_gain_grid",
        "kappa_redist_0p025_clean_crossing",
        "common_reference_energy_renormalization",
        "rank_and_novelty_diagnostic",
    ):
        _require_equal(source.get(key), "autonomous_design", f"{key} source class")
    for key in (
        "openai_hidden_radial_winding_distribution",
        "openai_hidden_swirl_redistribution_profile",
        "openai_hidden_particle_seeds_or_material_trajectories",
        "openai_hidden_camera_frame_time_alignment",
    ):
        _require_equal(source.get(key), "pending_unknown", f"{key} source class")

    definition = scope.get("redistribution_definition", {})
    _require("kappa_redist*h(r)" in str(definition.get("formula", "")), "redistribution formula weakened")
    _require_equal(definition.get("inner_support_interval"), [0.30, 1.05], "inner support interval")
    _require_equal(definition.get("outer_support_interval"), [0.95, 1.85], "outer support interval")
    _require("t=0.25" in str(definition.get("alpha_definition", "")), "alpha reference time missing")
    _require_equal(definition.get("diagnostic_gain_grid"), [0.0, 0.015, 0.025, 0.035, 0.05], "gain grid")
    _require_close(definition.get("first_clean_crossing"), 0.025, "diagnostic crossing")
    _require_equal(definition.get("production_coefficient_selected"), False, "production coefficient selection")
    _require_equal(definition.get("production_coefficient_bound_selected"), False, "production bound selection")
    _require_equal(
        definition.get("replacement_for_global_kappa_not_additive_fourth_degree"),
        True,
        "replacement routing semantics",
    )

    energy = scope.get("energy_semantics", {})
    _require_close(energy.get("reference_time"), 0.25, "energy reference time")
    _require_equal(
        energy.get("first_order_energy_derivative_at_zero_gain_is_zero_by_registered_quadrature"),
        True,
        "first-order energy statement",
    )
    for key in (
        "first_order_energy_neutral_means_exact_energy_invariant_for_finite_gain",
        "first_order_energy_neutral_means_energy_neutral_at_all_times",
        "balance_coefficient_is_parent_independent",
        "balance_coefficient_is_source_recovered",
        "positive_common_renormalization_is_navier_stokes_invariance",
        "parent_pressure_forcing_or_pde_receipt_transfers_through_redistribution",
        "parent_material_paths_transfer_through_redistribution",
    ):
        _require_equal(energy.get(key), False, key)
    for key in (
        "changing_parent_or_reference_time_requires_recomputing_balance",
        "finite_gain_rows_still_require_exact_reference_energy_check",
    ):
        _require_equal(energy.get(key), True, key)

    receipt = scope.get("pr447_diagnostic_receipt", {})
    _require_close(receipt.get("inner_swirl_energy_moment"), 0.24996768654689747, "inner energy moment")
    _require_close(receipt.get("outer_swirl_energy_moment"), 0.10349833370935456, "outer energy moment")
    _require_close(receipt.get("alpha"), 2.415185613025038, "alpha")
    _require_close(receipt.get("balanced_moment_defect"), 2.7755575615628914e-17, "balanced moment defect", 1e-15)
    _require_close(receipt.get("kinetic_energy_first_derivative"), 5.551115123125783e-17, "energy derivative", 1e-15)
    _require_close(receipt.get("alpha_refinement_relative_change_48_to_64"), 0.000512514, "alpha refinement", 1e-12)
    _require_close(receipt.get("maximum_energy_quadrature_refinement_change"), 0.0000129936, "energy refinement", 1e-12)
    _require_close(receipt.get("h_at_seed_radius_0p6"), 0.959189, "h(.6)", 1e-9)
    _require_close(receipt.get("h_at_seed_radius_0p9"), 0.569783, "h(.9)", 1e-9)
    _require_close(receipt.get("h_at_seed_radius_1p2"), -1.888197, "h(1.2)", 1e-9)
    _require_close(receipt.get("kappa_redist"), 0.025, "diagnostic redistribution gain")
    _require_close(receipt.get("common_reference_energy_scale"), 0.9985868857, "reference energy scale", 1e-10)
    _require_equal(
        receipt.get("seed_radius_angular_rate_percent_change"),
        {"r_0p6": 2.3770, "r_0p9": 1.4037, "r_1p2": -4.7400},
        "seed-radius angular-rate changes",
    )
    _require_equal(receipt.get("q90_q99_bins_changed"), False, "q90/q99 status")
    _require_equal(receipt.get("held_out_full_momentum_evaluated"), False, "held-out momentum evaluation")
    _require_equal(receipt.get("material_path_replay_run"), False, "material-path replay")
    _require_equal(receipt.get("visual_correspondence_verified"), False, "PR447 visual correspondence")
    _require_equal(receipt.get("pde_validated"), False, "PR447 PDE state")
    _require_equal(receipt.get("interpretation"), "representation_capacity_crossing_only", "PR447 interpretation")

    radius = scope.get("radius_attribution_semantics", {})
    _require_equal(radius.get("st006_is_openai_visual_truth"), False, "ST006/OpenAI truth separation")
    _require_equal(radius.get("st006_is_retained_repository_pde_benchmark"), True, "ST006 benchmark role")
    _require_equal(
        radius.get("radius_deficit_attribution_identifies_openai_hidden_radial_target"),
        False,
        "radius attribution/OpenAI separation",
    )
    _require_equal(radius.get("radius_deficit_attribution_may_motivate_autonomous_profile_design"), True, "routing use")
    _require_close(radius.get("r_0p6_net_deficit_contribution_percent"), 89.75, "r=.6 deficit contribution")
    _require_close(radius.get("r_0p9_net_deficit_contribution_percent"), 21.35, "r=.9 deficit contribution")
    _require_close(radius.get("r_1p2_net_deficit_contribution_percent"), -11.10, "r=1.2 deficit contribution")

    axis = scope.get("representation_and_axis_contract", {})
    _require_equal(axis.get("redistribution_changes_only_axisymmetric_azimuthal_component_before_common_scale"), True, "component scope")
    _require_equal(axis.get("finite_radial_multiplier_preserves_parent_u_theta_order_O_r_if_multiplier_is_bounded"), True, "axis order")
    _require_equal(axis.get("analytic_axis_regular_order_is_cartesian_axis_validation"), False, "axis validation separation")
    for key in (
        "materialized_child_requires_multiplier_finite_and_positive_on_support",
        "materialized_child_requires_fresh_axis_and_near_axis_finiteness_check",
        "materialized_child_requires_fresh_support_and_core_sign_check",
        "materialized_child_requires_fresh_independent_divergence_check",
    ):
        _require_equal(axis.get(key), True, key)

    sensitivity = scope.get("sensitivity_semantics", {})
    _require_equal(sensitivity.get("reported_rank"), 4, "reported rank")
    _require_equal(sensitivity.get("reported_singular_values"), [1.37132, 0.99141, 0.95470, 0.47449], "singular values")
    _require_close(sensitivity.get("reported_condition_number"), 2.8901, "condition number", 1e-9)
    _require_close(sensitivity.get("redistribution_novelty_fraction_outside_prior_span"), 0.9826, "novelty fraction", 1e-9)
    _require_close(sensitivity.get("cosine_with_global_kappa"), 0.0592, "global-kappa cosine", 1e-9)
    for key in (
        "rank_four_proves_pde_leverage",
        "rank_four_selects_a_production_coefficient",
        "rank_four_requires_adding_a_fourth_production_degree",
    ):
        _require_equal(sensitivity.get(key), False, key)

    promotion = scope.get("promotion_contract", {})
    for key in (
        "diagnostic_crossing_may_be_called_production_coefficient",
        "first_order_energy_neutrality_may_bypass_nontriviality_or_energy_checks",
        "fixed_probe_angular_rate_gain_may_be_called_material_path_gain",
        "st006_radius_attribution_may_be_called_openai_correspondence",
        "capacity_or_rank_evidence_may_replace_canonical_velocity",
        "capacity_or_rank_evidence_may_replace_live_scientific_route",
    ):
        _require_equal(promotion.get(key), False, key)
    for key in (
        "new_materialized_child_requires_new_representation_identity",
        "new_materialized_child_requires_new_candidate_sha",
        "save_load_replay_required",
        "fresh_material_path_replay_required_for_material_path_claim",
        "fresh_pressure_and_restricted_forcing_rebuild_required_for_pde_claim",
        "fresh_unused_held_out_full_momentum_and_divergence_required_for_pde_claim",
        "formal_thresholds_may_not_be_relaxed",
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
    _require_equal(truth.get("redistribution_capacity_registered"), True, "capacity registration")
    for key in (
        "redistribution_child_materialized",
        "production_redistribution_gain_selected",
        "production_redistribution_bound_selected",
        "material_path_gain_verified",
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
        "first_order_energy_neutrality_only": True,
        "diagnostic_crossing_only": True,
        "material_path_gain_verified": False,
        "canonical_velocity_unchanged": True,
        "live_route_unchanged": True,
        "visual_correspondence_verified": False,
        "pde_validated": False,
    }


if __name__ == "__main__":
    print(json.dumps(audit_energy_neutral_swirl_redistribution_governance(), indent=2, sort_keys=True))
