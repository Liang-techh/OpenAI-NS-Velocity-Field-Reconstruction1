"""Fail-closed governance for the ST052-M linear temporal activation experiment.

This module changes no velocity, coefficient, pressure, force, optimizer, threshold,
or route.  It keeps PR #587's time-dependent representation distinct from the
static PR #559 child and from the canonical Eq45 delivery, with special attention
to endpoint identity, time-derivative chain-rule terms, energy scope, and evidence
transfer.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
CONTRACT_PATH = ROOT / "configs" / "st052m_linear_temporal_activation_governance_contract.json"


class GovernanceError(RuntimeError):
    """Raised when the governed temporal-representation boundary is violated."""


def _load(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise GovernanceError(message)


def audit(
    root: Path = ROOT,
    *,
    contract_path: Path | None = None,
) -> dict[str, Any]:
    """Audit the frozen PR #587 evidence and current live CR001/delivery truth."""

    contract = _load(contract_path or (root / "configs" / CONTRACT_PATH.name))
    status = _load(root / "project_status.json")
    constraints = _load(root / "configs" / "constraints.json")
    delivery = _load(root / "configs" / "delivery_state_contract.json")
    prior_stage = _load(root / "configs" / "st052m_delivery_stage_accounting_contract.json")

    expected_classes = {
        "user_requirement",
        "public_source_fact",
        "autonomous_design",
        "pending_unknown",
    }
    _require(set(contract["allowed_source_classes"]) == expected_classes, "source vocabulary drift")
    _require(set(contract["source_classification"].values()) <= expected_classes, "invalid source class")
    _require(
        contract["source_classification"]["linear_activation_law_g_of_t"] == "autonomous_design",
        "autonomous time law relabelled as source fact",
    )
    _require(
        contract["source_classification"]["openai_hidden_time_parameterization"] == "pending_unknown",
        "hidden OpenAI time law claimed known",
    )

    snap = contract["snapshot"]
    _require(snap["active_integration_branch"] == "codex/cr001-constraints", "integration branch drift")
    _require(
        snap["active_integration_head"] == "97d1649da5913402003d3e62dacdf4c88208310d",
        "temporal-governance base snapshot drift",
    )
    _require(snap["temporal_pr"] == 587, "temporal PR drift")
    _require(
        snap["temporal_head"] == "0b93819095f6c8576a7571bdc2d2fbef4154944d",
        "temporal head drift",
    )
    _require(snap["static_endpoint_pr"] == 559, "static endpoint PR drift")
    _require(
        snap["static_endpoint_head"] == "39b106ad8cb8df2064cabead3a12682089575e74",
        "static endpoint head drift",
    )
    _require(snap["dedicated_workflow_run"] == 35428790564, "dedicated workflow identity drift")
    _require(snap["standard_tests_run"] == 35428790464, "standard workflow identity drift")
    _require(snap["st052_replay_run"] == 35428790494, "replay workflow identity drift")
    _require(snap["all_three_workflow_conclusions"] == "success", "PR587 exact-head workflows not green")

    temporal = contract["temporal_representation"]
    _require(temporal["time_interval"] == [0.25, 0.75], "temporal interval drift")
    _require(temporal["activation_formula"] == "g(t)=2*(t-0.25)", "activation law drift")
    _require(temporal["activation_derivative"] == 2.0, "activation derivative drift")
    _require(temporal["redistribution_kappa"] == 0.05, "redistribution gain drift")
    _require(temporal["taper_tau_formula"] == "tau(t)=0.05*g(t)", "tau law drift")
    _require(temporal["taper_tau_derivative"] == 0.1, "tau time derivative drift")
    _require(
        temporal["shoulder_beta_formula"] == "beta(t)=0.08837490297155456*g(t)",
        "beta law drift",
    )
    _require(
        math.isclose(temporal["shoulder_beta_derivative"], 2.0 * 0.08837490297155456, rel_tol=0.0, abs_tol=1e-16),
        "beta time derivative inconsistent with frozen activation",
    )
    _require(temporal["post_transform_common_scale"] == 1.0, "post-transform common scale introduced")
    _require(temporal["parameter_scan_performed"] is False, "temporal experiment relabelled as scan")
    _require(
        temporal["velocity_is_linear_interpolation_of_endpoint_velocities"] is False,
        "nonlinear transform composition mislabelled as velocity interpolation",
    )
    _require(temporal["time_law_is_part_of_candidate_identity"] is True, "time law dropped from identity")
    for key in (
        "materialized_on_live_constrained_ancestry",
        "versioned_candidate_id_assigned",
        "versioned_candidate_sha_assigned",
        "save_load_contract_registered",
        "unified_velocity_api_registered",
        "velocity_export_ready",
        "production_candidate_selected",
    ):
        _require(temporal[key] is False, f"temporal child delivery state over-promoted: {key}")

    endpoints = contract["endpoint_identities"]
    _require(endpoints["t025_activation"] == 0.0, "t=.25 endpoint activation drift")
    _require(endpoints["t075_activation"] == 1.0, "t=.75 endpoint activation drift")
    _require(
        endpoints["t025_identity"] == "untapered_st052m_kappa005_redistributed_control",
        "t=.25 endpoint identity drift",
    )
    _require(
        endpoints["t075_identity"] == "static_pr559_local_swirl_energy_child",
        "t=.75 endpoint identity drift",
    )
    _require(endpoints["endpoint_identity_implies_interior_identity"] is False, "endpoint laundering into interior identity")
    _require(
        endpoints["endpoint_identity_implies_material_path_identity"] is False,
        "endpoint equality laundered into path-history identity",
    )
    _require(
        endpoints["whole_interval_evidence_requires_temporal_child_specific_replay"] is True,
        "whole-interval evidence no longer requires temporal replay",
    )

    derivative = contract["time_derivative_semantics"]
    _require(derivative["activation_has_nonzero_time_derivative"] is True, "nonzero activation derivative lost")
    _require(
        derivative["extra_transform_time_derivative_terms_present_in_general"] is True,
        "chain-rule transform term erased",
    )
    for key in (
        "static_child_time_derivative_receipt_transfers",
        "parent_pressure_transfers",
        "parent_restricted_forcing_fit_transfers",
        "parent_or_static_momentum_residual_receipt_transfers",
    ):
        _require(derivative[key] is False, f"PDE evidence transfer over-promoted: {key}")
    _require(
        derivative["fresh_exact_child_time_derivative_required_for_pde_promotion"] is True,
        "fresh temporal derivative requirement weakened",
    )
    _require(
        derivative["fresh_compatible_pressure_and_restricted_forcing_required_for_pde_promotion"] is True,
        "fresh pressure/force rebuild requirement weakened",
    )

    energy = contract["energy_semantics"]
    _require(energy["reference_time"] == 0.25, "reference-energy time drift")
    _require(energy["reference_energy_relative_error_reported"] == 0.0, "reported endpoint energy identity drift")
    _require(energy["reference_time_identity_is_valid_nontriviality_evidence"] is True, "valid reference-time evidence discarded")
    _require(
        energy["reference_time_identity_implies_all_validation_time_energy_acceptance"] is False,
        "reference-time energy laundered into all-time energy acceptance",
    )
    _require(
        energy["order64_descriptive_energy_replaces_cr001_24_48_96_ladder"] is False,
        "order-64 diagnostic replaced CR001 energy ladder",
    )
    _require(energy["full_cr001_energy_checks_required_before_scientific_promotion"] is True, "full energy checks weakened")

    evidence = contract["pr587_evidence_scope"]
    _require(evidence["clean_linear_temporal_activation_capacity"] is True, "positive capacity receipt lost")
    _require(evidence["off_midplane_path_count"] == 48, "path count drift")
    _require(evidence["path_solver"] == "DOP853", "path solver drift")
    _require(evidence["late_field_matches_static_pr559_on_frozen_grid_abs_max"] == 0.0, "late endpoint equality lost")
    _require(evidence["sampled_cartesian_fd_divergence_max"] < 1e-5, "reported sampled divergence screen drift")
    _require(evidence["sampled_divergence_is_cr001_heldout_acceptance"] is False, "sampled screen laundered into CR001 acceptance")
    _require(evidence["public_openai_numeric_image_target_used"] is False, "OpenAI image target invented")
    _require(evidence["visual_acceptance_threshold_used"] is False, "visual acceptance threshold invented")
    _require(evidence["heldout_full_momentum_evaluated"] is False, "unrun momentum gate promoted")
    _require(evidence["counts_as_visual_correspondence_verified"] is False, "target-free capacity promoted to visual correspondence")
    _require(evidence["counts_as_pde_validated"] is False, "capacity receipt promoted to PDE validation")

    transfer = contract["static_evidence_transfer"]
    for key in (
        "static_pr559_full_interval_material_paths_transfer_to_temporal_child",
        "static_pr583_three_time_render_receipt_transfers_wholesale",
        "static_pr576_three_time_morphology_receipt_transfers_wholesale",
    ):
        _require(transfer[key] is False, f"static evidence transferred wholesale: {key}")
    _require(transfer["t075_same_time_pointwise_static_identity_can_be_replayed_as_equivalence_evidence"] is True,
             "legitimate t=.75 same-time equivalence evidence blocked")
    _require(transfer["t025_control_identity_can_be_replayed_as_equivalence_evidence"] is True,
             "legitimate t=.25 control equivalence evidence blocked")
    _require(transfer["interior_time_t050_requires_temporal_child_specific_evidence"] is True,
             "interior-time evidence requirement weakened")
    _require(transfer["material_paths_require_temporal_history_specific_evidence"] is True,
             "path-history evidence requirement weakened")

    identity = contract["delivery_identity_requirements"]
    _require(identity["automatic_promotion_from_pr587_allowed"] is False, "PR587 automatically promoted")
    _require(identity["candidate_identity_must_include_activation_formula"] is True, "activation law omitted from candidate identity")
    _require(identity["candidate_identity_must_include_static_endpoint_parameters"] is True, "endpoint parameters omitted from identity")
    _require(identity["same_endpoints_alone_are_sufficient_identity"] is False, "endpoints alone made full identity")
    _require(identity["same_static_pr559_parameters_alone_are_sufficient_identity"] is False, "static params alone made temporal identity")
    _require(identity["pde_failure_or_pending_blocks_callable_velocity_delivery"] is False, "PDE state improperly blocks velocity delivery")
    required = set(identity["required_for_materialization"])
    _require(
        {
            "exact_functional_binding_to_pr587_time_law",
            "versioned_candidate_id_and_sha",
            "deterministic_save_load_replay",
            "unified_velocity_xyzt_api_replay",
            "finite_output_and_support_checks",
        } <= required,
        "temporal materialization prerequisites weakened",
    )

    # The live route still points to materializing the static #559 child.  PR #587
    # is capacity evidence only and must not silently replace that route.
    _require(
        status["delivery_shortest_route"] == prior_stage["snapshot"]["project_status_delivery_shortest_route"],
        "live delivery route silently changed by temporal capacity experiment",
    )
    latest = status["latest_st052m_visual_candidate_evidence"]
    _require(latest["exact_child_source_pr"] == 559, "live source-child identity silently switched")
    _require(latest["production_candidate_selected"] is False, "live experimental child production-selected")
    _require(latest["velocity_export_ready"] is False, "live experimental child export-promoted")
    _require(latest["visual_correspondence_verified"] is False, "live visual correspondence promoted")
    _require(latest["pde_validated"] is False, "live PDE state promoted")

    canonical = contract["canonical_delivery"]
    _require(status["candidate_family"] == canonical["candidate_family"], "canonical family drift")
    _require(status["candidate_sha256"] == canonical["candidate_sha256"], "canonical candidate SHA drift")
    _require(status["velocity_api"] == canonical["velocity_api"], "canonical velocity API drift")
    _require(status["states"]["velocity_export_ready"] is canonical["velocity_export_ready"] is True,
             "canonical export readiness drift")
    for key in ("visualization_ready", "visual_correspondence_verified", "pde_validated", "paper_exact", "openai_field_identified"):
        _require(status["states"][key] is canonical[key] is False, f"canonical truth state promoted: {key}")

    _require(
        delivery["states"]["velocity_export_ready"]["positive_evidence"]
        == ["versioned_candidate_artifact", "unified_velocity_api", "finite_velocity_output", "reproducible_export"],
        "velocity export evidence contract drift",
    )
    _require(
        "pde_validated" in delivery["states"]["velocity_export_ready"]["does_not_imply"],
        "velocity export improperly coupled to PDE validation",
    )

    cr = contract["canonical_cr001"]
    _require(constraints["nu"] == cr["nu"], "nu drift")
    _require(constraints["domain"]["physical"] == cr["physical_domain"], "physical domain drift")
    _require(constraints["domain"]["evaluation_box"] == cr["evaluation_box"], "evaluation box drift")
    _require(constraints["domain"]["support"] == cr["support"], "support drift")
    _require(constraints["domain"]["time_interval"] == cr["time_interval"], "time interval drift")
    _require(constraints["forcing"]["mode"] == cr["force_mode"], "forcing family drift")
    _require(constraints["forcing"]["parameters"] == cr["force_bounds"], "forcing bounds drift")
    nontriviality = constraints["nontriviality"]
    for key in (
        "reference_energy",
        "reference_energy_abs_tolerance",
        "minimum_energy_each_validation_time",
        "maximum_energy_each_validation_time",
    ):
        _require(nontriviality[key] == cr[key], f"nontriviality contract drift: {key}")
    validation = constraints["validation"]
    _require(validation["seed"] == cr["validation_seed"], "validation seed drift")
    _require(validation["held_out_points"] == cr["held_out_points"], "held-out count drift")
    _require(validation["times"] == cr["validation_times"], "validation times drift")
    _require(validation["derivative_steps"] == cr["derivative_steps"], "derivative ladder drift")
    _require(validation["quadrature_orders_per_axis"] == cr["quadrature_orders_per_axis"], "energy quadrature ladder drift")
    thresholds = validation["thresholds"]
    for key in ("divergence_max", "divergence_L2", "pde_residual_max", "pde_residual_L2"):
        _require(thresholds[key] == cr[key], f"CR001 threshold drift: {key}")

    truth = contract["truth_boundary"]
    _require(truth["canonical_velocity_export_ready"] is True, "canonical export state corrupted")
    for key in (
        "temporal_child_velocity_export_ready",
        "production_candidate_selected",
        "visualization_ready",
        "visual_correspondence_verified",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    ):
        _require(truth[key] is False, f"truth-state promotion in temporal contract: {key}")

    _require(not any(contract["mutation_scope"].values()), "governance increment mutates scientific state")

    return {
        "task_id": contract["task_id"],
        "status": "pass",
        "temporal_pr_head": snap["temporal_head"],
        "activation_law_bound_to_identity": True,
        "endpoint_identity_scoped": True,
        "fresh_time_derivative_required_for_pde": True,
        "static_evidence_transfer_fail_closed": True,
        "temporal_child_velocity_export_ready": False,
        "canonical_candidate_unchanged": True,
        "cr001_unchanged": True,
    }


def main() -> None:
    print(json.dumps(audit(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
