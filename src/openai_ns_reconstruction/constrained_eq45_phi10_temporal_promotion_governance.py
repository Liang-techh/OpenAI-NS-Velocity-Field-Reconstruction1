"""Fail-closed promotion governance for the selected supported Phi(1,0) temporal trial.

This module governs claim semantics only.  The target-free slope-1.4 morphology
trial is a callable/exportable visualization candidate, but neither its internal
morphology improvement nor mixed pressure-free PDE diagnostics establish public
visual correspondence, formal PDE validity, paper exactness, or OpenAI-field
identity.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Mapping

from .constrained_eq45_supported_phi10_delivery_capsule import delivery_capsule


SCHEMA = "eq45_phi10_temporal_promotion_gate_v1"
TASK_ID = "CR002-EQ45-PHI10-TEMPORAL-PROMOTION-GATE-022"
EXPECTED_BASE_SHA256 = "2fdff812c22131d56eac1b7d6e455187ff3207500c6385aa9abf282a8e3d7b1d"
_EXPECTED_CLASSIFICATION = {
    "velocity_visual_delivery": "user_requirement",
    "eq45_backbone": "public_source_fact",
    "phi10_temporal_mode_and_slope_selection": "autonomous_design",
    "target_free_internal_morphology_metric": "autonomous_design",
    "public_visual_correspondence": "pending_unknown",
    "pde_validity": "pending_unknown",
}
_EXPECTED_STATE = {
    "velocity_export_ready": True,
    "visualization_candidate_only": True,
    "visualization_ready": False,
    "visual_correspondence_verified": False,
    "physical_support_validated": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "canonical_velocity_changed": False,
    "production_slope_promoted": False,
}
_REQUIRED_VIS_READY = {
    "callable_velocity",
    "supported_visualization_entry_point",
    "candidate_specific_resolution_sanity",
    "whole_domain_temporal_morphology_review_resolved",
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
    "early_morphology_improvement_implies_visualization_ready",
    "whole_domain_resolution_stability_implies_public_visual_correspondence",
    "mixed_seed_pde_changes_imply_pde_validation",
    "restricted_force_projection_implies_pde_validation",
    "green_ci_implies_visual_or_pde_acceptance",
    "velocity_export_implies_visualization_ready",
    "visual_similarity_implies_pde_validation",
    "visual_correspondence_implies_paper_exact",
    "pde_failure_blocks_velocity_export",
    "open_sibling_evidence_is_integrated_acceptance",
    "production_slope_is_promoted",
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


def audit_phi10_temporal_promotion_gate(
    contract: Mapping[str, Any], *, repo_root: str | Path
) -> dict[str, Any]:
    """Audit promotion semantics for the selected temporal visualization trial."""
    root = Path(repo_root)
    _require(contract.get("schema") == SCHEMA, "unexpected Phi10 promotion schema")
    _require(contract.get("task_id") == TASK_ID, "Phi10 promotion task identity drifted")

    delivery = _read_json(root / "configs" / "delivery_state_contract.json")
    constraints = _read_json(root / "configs" / "constraints.json")
    vocabulary = delivery.get("classification_vocabulary")
    _require(isinstance(vocabulary, Mapping), "delivery classification vocabulary is missing")

    classification = contract.get("classification")
    _require(classification == _EXPECTED_CLASSIFICATION, "Phi10 source classification drifted")
    _require(
        set(classification.values()).issubset(set(vocabulary)),
        "Phi10 promotion gate uses a noncanonical source class",
    )

    capsule = delivery_capsule()
    candidate = contract.get("candidate")
    _require(isinstance(candidate, Mapping), "candidate section is missing")
    _require(
        candidate.get("source_capsule_task") == capsule.get("task_id"),
        "source capsule task drifted",
    )
    _require(
        candidate.get("base_supported_sha256") == capsule.get("base_supported_sha256") == EXPECTED_BASE_SHA256,
        "supported parent identity drifted",
    )
    _require(candidate.get("temporal_mode") == ["phi", 1, 0], "temporal mode drifted")
    _require(float(candidate.get("slope")) == 1.4, "selected temporal slope drifted")
    _require(
        candidate.get("coefficient_start_mid_end") == [-1.7, -0.3, 1.1],
        "temporal endpoint coefficients drifted",
    )
    _require(candidate.get("coefficient_bound") == [-4.0, 4.0], "coefficient bound drifted")
    _require(candidate.get("time_interval") == [0.25, 0.75], "candidate time interval drifted")
    _require(
        candidate.get("spatial_box") == [[-2.0, 2.0], [-2.0, 2.0], [-2.0, 2.0]],
        "candidate spatial box drifted",
    )

    domain = constraints.get("domain")
    forcing = constraints.get("forcing")
    validation = constraints.get("validation")
    _require(isinstance(domain, Mapping), "CR001 domain contract is missing")
    _require(isinstance(forcing, Mapping), "CR001 forcing contract is missing")
    _require(isinstance(validation, Mapping), "CR001 validation contract is missing")
    _require(domain.get("evaluation_box") == candidate.get("spatial_box"), "delivery box drifted from CR001")
    _require(domain.get("time_interval") == candidate.get("time_interval"), "time interval drifted from CR001")
    _require(forcing.get("mode") == "restricted_two_parameter_family", "restricted forcing family drifted")
    _require(
        "No residual-dependent basis or pointwise free force" in str(forcing.get("restriction")),
        "free/residual-defined force restriction drifted",
    )
    thresholds = validation.get("thresholds")
    _require(isinstance(thresholds, Mapping), "CR001 validation thresholds are missing")
    _require(thresholds.get("pde_residual_max") == 0.001, "PDE max threshold drifted")
    _require(thresholds.get("pde_residual_L2") == 0.001, "PDE L2 threshold drifted")

    current_state = contract.get("current_state")
    _require(current_state == _EXPECTED_STATE, "Phi10 current state drifted")
    _require(capsule.get("status") == _EXPECTED_STATE, "live Phi10 capsule truth state drifted")

    evidence = contract.get("evidence_scope")
    _require(isinstance(evidence, Mapping), "evidence_scope is missing")
    selection = evidence.get("target_free_selection")
    _require(isinstance(selection, Mapping), "target-free selection evidence is missing")
    _require(selection.get("source_pr") == 143, "target-free selection source drifted")
    _require(selection.get("status") == "integrated_in_capsule_ancestry", "selection ancestry drifted")
    _require(selection.get("may_support_visualization_trial_selection") is True, "trial-selection evidence was disabled")
    _require(selection.get("may_support_public_visual_correspondence") is False, "target-free screen cannot prove public correspondence")
    _require(selection.get("may_support_pde_validation") is False, "target-free screen cannot prove PDE validity")
    _require(capsule["selection_evidence"].get("openai_image_used_as_fit_target") is False, "trial was not selected against a public image")

    fixed_pde = evidence.get("fixed_holdout_pde_tradeoff")
    seed_pde = evidence.get("fresh_seed_pde_generalization")
    _require(isinstance(fixed_pde, Mapping) and isinstance(seed_pde, Mapping), "PDE evidence scope is missing")
    _require(fixed_pde.get("source_pr") == 144, "fixed-holdout PDE source drifted")
    _require(seed_pde.get("source_pr") == 145, "seed-generalization PDE source drifted")
    _require(float(fixed_pde.get("finest_zero_force_fractional_change")) > 0.0, "fixed holdout must retain the observed worsening")
    _require(fixed_pde.get("formal_pde_gate_assessed") is False, "fixed-probe diagnostic is not the formal PDE gate")
    changes = seed_pde.get("fresh_seed_fractional_changes")
    _require(isinstance(changes, list) and any(x < 0 for x in changes) and any(x > 0 for x in changes), "fresh-seed PDE tradeoff must remain mixed")
    _require(seed_pde.get("uniform_generalization") is False, "mixed seed result cannot be called uniform")
    _require(seed_pde.get("formal_pde_gate_assessed") is False, "seed diagnostic is not the formal PDE gate")
    _require(fixed_pde.get("may_support_pde_validation") is False, "fixed holdout cannot promote PDE validation")
    _require(seed_pde.get("may_support_pde_validation") is False, "mixed seed evidence cannot promote PDE validation")

    whole_domain = evidence.get("whole_domain_morphology")
    _require(isinstance(whole_domain, Mapping), "whole-domain morphology evidence is missing")
    _require(whole_domain.get("source_pr") == 146, "whole-domain morphology source drifted")
    _require(whole_domain.get("status") == "open_unconsumed_sibling_evidence", "PR #146 must remain unconsumed in this ancestry")
    _require(whole_domain.get("may_trigger_followup_review") is True, "whole-domain evidence may guide follow-up")
    _require(whole_domain.get("may_promote_current_ancestry_state") is False, "open sibling evidence cannot promote this ancestry")
    _require(whole_domain.get("reported_visualization_ready") is False, "PR #146 did not claim visualization readiness")
    _require(whole_domain.get("reported_visual_correspondence_verified") is False, "PR #146 did not claim public correspondence")
    _require(capsule["pending"].get("whole_domain_morphology_acceptance") == "PR146_open_unconsumed", "capsule sibling status drifted")

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
        _require(required.issubset(set(state.get("positive_evidence", []))), f"Phi10 gate cannot weaken global {state_name} evidence")

    policy = contract.get("policy")
    _require(isinstance(policy, Mapping), "promotion policy is missing")
    _require(set(policy) == _FALSE_POLICY_KEYS, "promotion policy key set drifted")
    for key in _FALSE_POLICY_KEYS:
        _require(policy.get(key) is False, f"forbidden promotion inference enabled: {key}")

    boundary = contract.get("truth_boundary")
    _require(isinstance(boundary, Mapping), "truth_boundary is missing")
    for key in (
        "velocity_changed_by_this_gate",
        "temporal_slope_changed_by_this_gate",
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
        "whole_domain_morphology_status": "open_unconsumed_sibling_evidence",
        "pde_tradeoff": "mixed_seed_sensitive_not_formal_gate",
        "next_promotion_gate": "integrate_whole_domain_morphology_then_public_observable_3d_comparison",
    }


def audit_phi10_temporal_promotion_gate_file(
    path: str | Path, *, repo_root: str | Path
) -> dict[str, Any]:
    return audit_phi10_temporal_promotion_gate(_read_json(Path(path)), repo_root=repo_root)
