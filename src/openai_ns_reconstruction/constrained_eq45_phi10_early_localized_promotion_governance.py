"""Fail-closed promotion governance for the early-localized supported Phi(1,0) field.

This module changes claim semantics only.  The quadratic temporal candidate is a
callable, serializable visualization/research artifact.  Its target-free morphology
screen, fixed/fresh pressure-free PDE diagnostics, green CI, or exportability do not
establish public visual correspondence, formal PDE validity, paper exactness, or
OpenAI-field identity.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Mapping

from .constrained_eq45_supported_phi10_early_localized_delivery_capsule import (
    delivery_capsule,
)


SCHEMA = "eq45_phi10_early_localized_promotion_gate_v1"
TASK_ID = "CR002-EQ45-PHI10-EARLY-LOCALIZED-PROMOTION-GATE-023"
EXPECTED_BASE_SHA256 = "2fdff812c22131d56eac1b7d6e455187ff3207500c6385aa9abf282a8e3d7b1d"
EXPECTED_CANDIDATE_SHA256 = "5b08abbfdce0e96d39a5fc394ceedea0d0ff549cb921ff50e84004355c0e2bfe"

_EXPECTED_CLASSIFICATION = {
    "velocity_visual_delivery": "user_requirement",
    "eq45_backbone": "public_source_fact",
    "phi10_early_localized_temporal_shape": "autonomous_design",
    "target_free_internal_morphology": "autonomous_design",
    "public_visual_correspondence": "pending_unknown",
    "pde_validity": "pending_unknown",
    "candidate_specific_energy_acceptance": "pending_unknown",
}

_EXPECTED_STATE = {
    "callable_serializable": True,
    "velocity_export_ready": True,
    "visualization_candidate_only": True,
    "canonical_velocity_changed": False,
    "screened_temporal_shape_materialized": True,
    "production_temporal_shape_promoted": False,
    "physical_support_connection_implemented": True,
    "physical_support_validated": False,
    "visualization_ready": False,
    "visual_correspondence_verified": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}

_REQUIRED_VIS_READY = {
    "callable_velocity",
    "supported_visualization_entry_point",
    "candidate_specific_resolution_sanity",
    "whole_domain_temporal_morphology_review_resolved",
    "three_dimensional_render_review",
}
_REQUIRED_VIS_CORR = {
    "public_reference_only",
    "candidate_specific_independent_visual_diagnostics",
    "candidate_specific_resolution_stability",
    "explicit_public_observable_geometry_and_time_comparison",
}
_REQUIRED_PDE = {
    "independent_validator",
    "preregistered_thresholds",
    "held_out_validation",
}
_FALSE_POLICY_KEYS = {
    "target_free_selection_implies_public_visual_correspondence",
    "early_endpoint_morphology_improvement_implies_visualization_ready",
    "static_midpoint_or_endpoint_identity_implies_intermediate_visual_acceptance",
    "intermediate_resolution_stability_implies_public_visual_correspondence",
    "fixed_probe_pressure_free_diagnostic_implies_pde_validation",
    "fresh_seed_pressure_free_diagnostic_implies_pde_validation",
    "green_ci_implies_visual_or_pde_acceptance",
    "velocity_export_implies_visualization_ready",
    "visual_similarity_implies_pde_validation",
    "visual_correspondence_implies_paper_exact",
    "pde_failure_blocks_velocity_export",
    "candidate_energy_failure_blocks_velocity_export",
    "open_sibling_evidence_is_integrated_acceptance",
    "production_temporal_shape_is_promoted",
    "canonical_velocity_is_changed",
}


def _read_json(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"{path} must contain a JSON object")
    return payload


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def audit_phi10_early_localized_promotion_gate(
    contract: Mapping[str, Any], *, repo_root: str | Path
) -> dict[str, Any]:
    """Audit source, constraint, evidence, and promotion semantics for the field."""
    root = Path(repo_root)
    _require(contract.get("schema") == SCHEMA, "unexpected early-localized promotion schema")
    _require(contract.get("task_id") == TASK_ID, "early-localized promotion task identity drifted")

    delivery = _read_json(root / "configs" / "delivery_state_contract.json")
    constraints = _read_json(root / "configs" / "constraints.json")
    vocabulary = delivery.get("classification_vocabulary")
    _require(isinstance(vocabulary, Mapping), "delivery classification vocabulary is missing")

    classification = contract.get("classification")
    _require(classification == _EXPECTED_CLASSIFICATION, "early-localized source classification drifted")
    _require(
        set(classification.values()).issubset(set(vocabulary)),
        "early-localized gate uses a noncanonical source class",
    )

    capsule = delivery_capsule()
    candidate = contract.get("candidate")
    _require(isinstance(candidate, Mapping), "candidate section is missing")
    _require(
        candidate.get("source_capsule_task") == capsule.get("task_id"),
        "source delivery capsule task drifted",
    )
    _require(
        candidate.get("base_supported_sha256")
        == capsule.get("base_supported_sha256")
        == EXPECTED_BASE_SHA256,
        "supported parent identity drifted",
    )
    _require(
        candidate.get("candidate_sha256")
        == capsule.get("candidate_sha256")
        == EXPECTED_CANDIDATE_SHA256,
        "early-localized candidate identity drifted",
    )
    _require(candidate.get("temporal_mode") == ["phi", 1, 0], "temporal mode drifted")
    _require(float(candidate.get("early_delta")) == -1.4, "early_delta drifted")
    _require(
        candidate.get("schedule") == "0.5*early_delta*tau*(tau-1)",
        "quadratic temporal schedule drifted",
    )
    _require(
        candidate.get("coefficient_start_mid_intermediate_end") == [-1.7, -0.3, -0.125, -0.3],
        "registered temporal coefficients drifted",
    )
    _require(candidate.get("coefficient_bound") == [-4.0, 4.0], "coefficient bound drifted")
    _require(candidate.get("time_interval") == [0.25, 0.75], "candidate time interval drifted")
    _require(
        candidate.get("spatial_box") == [[-2.0, 2.0], [-2.0, 2.0], [-2.0, 2.0]],
        "candidate spatial box drifted",
    )

    mode = capsule.get("temporal_mode")
    _require(isinstance(mode, Mapping), "live capsule temporal mode is missing")
    _require(mode.get("family") == "phi" and mode.get("index") == [1, 0], "live capsule mode drifted")
    _require(float(mode.get("early_delta")) == -1.4, "live capsule early_delta drifted")
    _require(mode.get("schedule") == candidate.get("schedule"), "live capsule schedule drifted")
    _require(mode.get("coefficient_bound") == candidate.get("coefficient_bound"), "live coefficient bound drifted")

    capsule_delivery = capsule.get("delivery")
    _require(isinstance(capsule_delivery, Mapping), "live capsule delivery section is missing")
    _require(capsule_delivery.get("time_interval") == candidate.get("time_interval"), "capsule time interval drifted")
    _require(capsule_delivery.get("spatial_box") == candidate.get("spatial_box"), "capsule spatial box drifted")
    _require(
        set(capsule_delivery.get("velocity_interfaces", []))
        == {"velocity", "at_points", "velocity_xyz", "grid"},
        "public velocity interface inventory drifted",
    )
    _require(
        set(capsule_delivery.get("serialization_interfaces", [])) == {"save_json", "load_json"},
        "serialization interface inventory drifted",
    )
    _require(capsule_delivery.get("component_order") == ["u", "v", "w"], "component order drifted")

    domain = constraints.get("domain")
    forcing = constraints.get("forcing")
    nontriviality = constraints.get("nontriviality")
    validation = constraints.get("validation")
    _require(float(constraints.get("nu")) == 0.01, "CR001 viscosity drifted")
    _require(isinstance(domain, Mapping), "CR001 domain contract is missing")
    _require(domain.get("physical") == "R^3", "CR001 physical domain drifted")
    _require(domain.get("evaluation_box") == candidate.get("spatial_box"), "delivery box drifted from CR001")
    _require(domain.get("time_interval") == candidate.get("time_interval"), "time interval drifted from CR001")
    _require(domain.get("support") == "r < 2 and abs(z) < 2", "CR001 support drifted")
    _require(
        domain.get("boundary") == "smooth zero extension of velocity and pressure outside support",
        "CR001 boundary convention drifted",
    )
    _require(isinstance(forcing, Mapping), "CR001 forcing contract is missing")
    _require(forcing.get("mode") == "restricted_two_parameter_family", "restricted forcing family drifted")
    _require(
        "No residual-dependent basis or pointwise free force" in str(forcing.get("restriction")),
        "free/residual-defined force restriction drifted",
    )
    _require(isinstance(nontriviality, Mapping), "CR001 nontriviality contract is missing")
    _require(float(nontriviality.get("reference_time")) == 0.25, "energy reference time drifted")
    _require(float(nontriviality.get("reference_energy")) == 1.0, "energy reference value drifted")
    _require(float(nontriviality.get("reference_energy_abs_tolerance")) == 0.001, "energy tolerance drifted")
    _require(isinstance(validation, Mapping), "CR001 validation contract is missing")
    _require(validation.get("held_out_points") == 4096, "held-out validation size drifted")
    _require(validation.get("derivative_steps") == [0.02, 0.01, 0.005], "derivative ladder drifted")
    thresholds = validation.get("thresholds")
    _require(isinstance(thresholds, Mapping), "CR001 validation thresholds are missing")
    _require(thresholds.get("divergence_max") == 1e-5, "divergence max threshold drifted")
    _require(thresholds.get("divergence_L2") == 1e-5, "divergence L2 threshold drifted")
    _require(thresholds.get("pde_residual_max") == 0.001, "PDE max threshold drifted")
    _require(thresholds.get("pde_residual_L2") == 0.001, "PDE L2 threshold drifted")

    current_state = contract.get("current_state")
    _require(current_state == _EXPECTED_STATE, "early-localized current state drifted")
    _require(capsule.get("status") == _EXPECTED_STATE, "live delivery capsule truth state drifted")

    evidence = contract.get("evidence_scope")
    _require(isinstance(evidence, Mapping), "evidence_scope is missing")
    selection = evidence.get("target_free_selection")
    _require(isinstance(selection, Mapping), "target-free selection evidence is missing")
    _require(selection.get("source_pr") == 150, "target-free selection source drifted")
    _require(selection.get("status") == "integrated_in_capsule_ancestry", "selection ancestry drifted")
    _require(selection.get("may_support_visualization_trial_selection") is True, "trial-selection evidence was disabled")
    _require(selection.get("may_support_public_visual_correspondence") is False, "target-free selection cannot prove public correspondence")
    _require(selection.get("may_support_pde_validation") is False, "target-free selection cannot prove PDE validity")
    capsule_selection = capsule.get("selection_evidence")
    _require(isinstance(capsule_selection, Mapping), "capsule selection evidence is missing")
    _require(capsule_selection.get("capacity_pr") == 150, "capsule selection source drifted")
    _require(capsule_selection.get("whole_domain_morphology_pr") == 146, "capsule morphology provenance drifted")
    _require(capsule_selection.get("public_openai_reference_used_for_selection") is False, "target-free selection was relabeled as public fitting")
    _require(capsule_selection.get("visual_correspondence_established") is False, "internal morphology cannot establish public correspondence")

    early = evidence.get("early_endpoint_whole_domain_morphology")
    _require(isinstance(early, Mapping), "early-endpoint morphology evidence is missing")
    _require(early.get("source_pr") == 146, "early-endpoint morphology source drifted")
    _require(early.get("status") == "bound_as_numeric_provenance_via_capsule", "early-endpoint provenance status drifted")
    _require(float(early.get("early_radial_q99_trial")) < float(early.get("early_radial_q99_baseline")), "registered early radial improvement drifted")
    _require(
        float(early.get("early_radial_outer_vorticity2_fraction_trial"))
        < float(early.get("early_radial_outer_vorticity2_fraction_baseline")),
        "registered early collar-energy improvement drifted",
    )
    _require(early.get("may_support_internal_morphology_review") is True, "early morphology evidence was disabled")
    _require(early.get("may_support_public_visual_correspondence") is False, "internal morphology cannot prove public correspondence")

    intermediate = evidence.get("intermediate_time_morphology")
    _require(isinstance(intermediate, Mapping), "intermediate morphology evidence is missing")
    _require(intermediate.get("source_pr") == 154, "intermediate morphology source drifted")
    _require(intermediate.get("source_head") == "385c4ae7c4cda43874cb451aa5a5d2b8966e6dfb", "intermediate morphology head drifted")
    _require(intermediate.get("source_ci_run") == 35212652157 and intermediate.get("source_ci_success") is True, "intermediate morphology CI evidence drifted")
    _require(intermediate.get("status") == "open_unconsumed_sibling_evidence", "PR #154 must remain unconsumed in this ancestry")
    _require(float(intermediate.get("time")) == 0.625, "intermediate morphology time drifted")
    _require(intermediate.get("finest_resolution") == 81, "intermediate morphology resolution drifted")
    _require(intermediate.get("may_trigger_followup_review") is True, "intermediate morphology may guide follow-up")
    _require(intermediate.get("may_promote_current_ancestry_state") is False, "open sibling evidence cannot promote this ancestry")
    _require(intermediate.get("reported_visualization_ready") is False, "PR #154 did not establish visualization readiness")
    _require(intermediate.get("reported_visual_correspondence_verified") is False, "PR #154 did not establish public correspondence")

    fixed_pde = evidence.get("fixed_holdout_pde_tradeoff")
    fresh_pde = evidence.get("fresh_seed_pde_generalization")
    _require(isinstance(fixed_pde, Mapping) and isinstance(fresh_pde, Mapping), "PDE evidence scope is missing")
    _require(fixed_pde.get("source_pr") == 152, "fixed-holdout PDE source drifted")
    _require(float(fixed_pde.get("early_localized_vs_static_fractional_change")) > 0.0, "fixed-holdout warning must remain a worsening")
    _require(fixed_pde.get("formal_pde_gate_assessed") is False, "fixed pressure-free diagnostic is not the formal PDE gate")
    _require(fixed_pde.get("may_support_pde_validation") is False, "fixed pressure-free diagnostic cannot promote PDE validity")
    _require(fresh_pde.get("source_pr") == 153, "fresh-seed PDE source drifted")
    _require(fresh_pde.get("source_head") == "73429a605c750707cbcc5e6e26ddcc5511894ac4", "fresh-seed PDE head drifted")
    _require(fresh_pde.get("source_ci_run") == 35211476178 and fresh_pde.get("source_ci_success") is True, "fresh-seed PDE CI evidence drifted")
    _require(fresh_pde.get("validation_seeds") == [914457, 914471, 914483], "fresh-seed set drifted")
    _require(fresh_pde.get("status") == "integrated_provenance_without_bound_numeric_summary", "fresh-seed evidence scope drifted")
    _require(fresh_pde.get("numeric_summary_bound_in_capsule") is False, "unbound fresh-seed numerics cannot be invented")
    _require(fresh_pde.get("formal_pde_gate_assessed") is False, "fresh-seed diagnostic is not the formal PDE gate")
    _require(fresh_pde.get("may_support_pde_validation") is False, "fresh-seed diagnostic cannot promote PDE validity")
    capsule_pde = capsule.get("pde_evidence")
    _require(isinstance(capsule_pde, Mapping), "capsule PDE evidence is missing")
    _require(capsule_pde.get("fresh_seed_numeric_summary_bound_here") is False, "capsule must keep fresh-seed numeric summary unbound")
    _require(capsule_pde.get("formal_pde_gate_assessed") is False, "capsule diagnostics are not formal PDE acceptance")
    _require(capsule_pde.get("pde_validated") is False, "capsule cannot over-promote PDE state")

    energy = evidence.get("candidate_energy")
    _require(isinstance(energy, Mapping), "candidate energy evidence is missing")
    _require(energy.get("status") == "pending_candidate_specific_revalidation", "candidate energy must remain pending")
    _require(energy.get("may_inherit_parent_energy_acceptance") is False, "temporal child cannot inherit parent energy acceptance")
    _require(energy.get("failure_blocks_velocity_export") is False, "candidate energy failure cannot block velocity export")
    capsule_pending = capsule.get("pending")
    _require(isinstance(capsule_pending, Mapping), "capsule pending state is missing")
    _require(capsule_pending.get("candidate_specific_energy_acceptance") is True, "capsule energy acceptance must remain pending")
    _require(capsule_pending.get("public_observable_comparison") is True, "public comparison must remain pending")
    _require(capsule_pending.get("three_dimensional_render_review") is True, "3D render review must remain pending")
    _require(capsule_pending.get("formal_pde_gate") is True, "formal PDE gate must remain pending")
    _require(capsule_pending.get("canonical_promotion") is False, "candidate must not be canonically promoted")

    requirements = contract.get("promotion_requirements")
    _require(isinstance(requirements, Mapping), "promotion requirements are missing")
    _require(set(requirements.get("visualization_ready", [])) == _REQUIRED_VIS_READY, "visualization-ready requirements drifted")
    _require(set(requirements.get("visual_correspondence_verified", [])) == _REQUIRED_VIS_CORR, "visual-correspondence requirements drifted")
    _require(set(requirements.get("pde_validated", [])) == _REQUIRED_PDE, "PDE-validation requirements drifted")

    global_states = delivery.get("states")
    _require(isinstance(global_states, Mapping), "global delivery states are missing")
    for state_name, required in (
        ("visualization_ready", {"callable_velocity", "supported_visualization_entry_point", "resolution_sanity_check"}),
        ("visual_correspondence_verified", {"public_reference_only", "independent_visual_diagnostics", "resolution_stability"}),
        ("pde_validated", _REQUIRED_PDE),
    ):
        state = global_states.get(state_name)
        _require(isinstance(state, Mapping), f"global {state_name} state is missing")
        _require(required.issubset(set(state.get("positive_evidence", []))), f"early-localized gate cannot weaken global {state_name} evidence")

    policy = contract.get("policy")
    _require(isinstance(policy, Mapping), "promotion policy is missing")
    _require(set(policy) == _FALSE_POLICY_KEYS, "promotion policy key set drifted")
    for key in _FALSE_POLICY_KEYS:
        _require(policy.get(key) is False, f"forbidden promotion inference enabled: {key}")

    boundary = contract.get("truth_boundary")
    _require(isinstance(boundary, Mapping), "truth_boundary is missing")
    for key in (
        "velocity_changed_by_this_gate",
        "temporal_schedule_changed_by_this_gate",
        "support_transform_changed",
        "forcing_family_changed",
        "thresholds_changed",
        "candidate_promoted",
    ):
        _require(boundary.get(key) is False, f"governance-only increment drifted: {key}")

    return {
        "contract_pass": True,
        "candidate_sha256": capsule["candidate_sha256"],
        "base_supported_sha256": capsule["base_supported_sha256"],
        "velocity_export_ready": True,
        "visualization_candidate_only": True,
        "visualization_ready": False,
        "visual_correspondence_verified": False,
        "pde_validated": False,
        "paper_exact": False,
        "openai_field_identified": False,
        "candidate_energy_status": "pending_candidate_specific_revalidation",
        "intermediate_morphology_status": "open_unconsumed_sibling_evidence",
        "pde_evidence_scope": "pressure_free_diagnostics_not_formal_gate",
        "next_promotion_gate": "consume_intermediate_morphology_then_candidate_energy_and_public_observable_3d_review",
    }


def audit_phi10_early_localized_promotion_gate_file(
    path: str | Path, *, repo_root: str | Path
) -> dict[str, Any]:
    return audit_phi10_early_localized_promotion_gate(_read_json(Path(path)), repo_root=repo_root)
