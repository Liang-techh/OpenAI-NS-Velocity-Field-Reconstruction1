"""Fail-closed promotion governance for the cubic-localized supported Phi(1,0) field.

The cubic candidate is callable and serializable, but exact snapshot equality with its
static supported parent at selected times is not temporal-derivative identity.  This
module governs that distinction without changing the velocity field or any scientific
threshold.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Mapping

from .constrained_eq45_supported_phi10_cubic_delivery_capsule import delivery_capsule


SCHEMA = "eq45_phi10_cubic_promotion_gate_v1"
TASK_ID = "CR002-EQ45-PHI10-CUBIC-PROMOTION-GATE-024"
EXPECTED_BASE_SHA256 = "2fdff812c22131d56eac1b7d6e455187ff3207500c6385aa9abf282a8e3d7b1d"
EXPECTED_CANDIDATE_SHA256 = "310fc2ad1bf3d721860369b9a58c7599202182964158708867691f4445f21576"
EXPECTED_STATIC_RETURN_TIMES = [0.5, 0.625, 0.75]
EXPECTED_OFF_KEYFRAME_TIMES = [0.375, 0.6875]

_EXPECTED_CLASSIFICATION = {
    "velocity_visual_delivery": "user_requirement",
    "eq45_backbone": "public_source_fact",
    "phi10_cubic_temporal_shape": "autonomous_design",
    "snapshot_identity_observation": "autonomous_design",
    "temporal_derivative_validation": "autonomous_design",
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
    "spatial_basis_grown": False,
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
    "snapshot_identity_implies_time_derivative_identity",
    "snapshot_identity_allows_parent_pde_inheritance",
    "snapshot_identity_allows_parent_energy_inheritance",
    "snapshot_identity_allows_parent_visual_inheritance",
    "snapshot_identity_allows_parent_canonical_inheritance",
    "target_free_selection_implies_public_visual_correspondence",
    "temporal_derivative_convergence_implies_pde_validation",
    "restricted_force_diagnostic_implies_pde_validation",
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


def audit_phi10_cubic_promotion_gate(
    contract: Mapping[str, Any], *, repo_root: str | Path
) -> dict[str, Any]:
    """Audit source, constraint, identity-inheritance and promotion semantics."""
    root = Path(repo_root)
    _require(contract.get("schema") == SCHEMA, "unexpected cubic promotion schema")
    _require(contract.get("task_id") == TASK_ID, "cubic promotion task identity drifted")

    delivery = _read_json(root / "configs" / "delivery_state_contract.json")
    constraints = _read_json(root / "configs" / "constraints.json")
    vocabulary = delivery.get("classification_vocabulary")
    _require(isinstance(vocabulary, Mapping), "delivery classification vocabulary is missing")

    classification = contract.get("classification")
    _require(classification == _EXPECTED_CLASSIFICATION, "cubic source classification drifted")
    _require(
        set(classification.values()).issubset(set(vocabulary)),
        "cubic gate uses a noncanonical source class",
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
        "cubic candidate identity drifted",
    )
    _require(candidate.get("temporal_mode") == ["phi", 1, 0], "temporal mode drifted")
    _require(float(candidate.get("early_delta")) == -1.4, "early_delta drifted")
    _require(
        candidate.get("schedule") == "early_delta*tau*(tau-0.5)*(tau-1)/(-3)",
        "cubic temporal schedule drifted",
    )
    _require(candidate.get("coefficient_bound") == [-4.0, 4.0], "coefficient bound drifted")
    _require(candidate.get("static_return_times") == EXPECTED_STATIC_RETURN_TIMES, "static-return times drifted")
    _require(candidate.get("off_keyframe_nontrivial_times") == EXPECTED_OFF_KEYFRAME_TIMES, "off-keyframe times drifted")
    _require(candidate.get("time_interval") == [0.25, 0.75], "candidate time interval drifted")
    _require(
        candidate.get("spatial_box") == [[-2.0, 2.0], [-2.0, 2.0], [-2.0, 2.0]],
        "candidate spatial box drifted",
    )

    representation = capsule.get("representation_evidence")
    temporal = capsule.get("temporal_derivative_evidence")
    _require(isinstance(representation, Mapping), "capsule representation evidence is missing")
    _require(isinstance(temporal, Mapping), "capsule temporal-derivative evidence is missing")
    _require(representation.get("capacity_pr") == 158, "cubic selection provenance drifted")
    _require(representation.get("materialization_pr") == 159, "cubic materialization provenance drifted")
    _require(representation.get("static_return_times") == EXPECTED_STATIC_RETURN_TIMES, "capsule static-return times drifted")
    _require(representation.get("off_keyframe_nontrivial_times") == EXPECTED_OFF_KEYFRAME_TIMES, "capsule off-keyframe times drifted")
    _require(representation.get("public_openai_reference_used_for_selection") is False, "target-free selection was relabeled as public fitting")
    _require(representation.get("visual_correspondence_established") is False, "internal representation evidence cannot establish public correspondence")
    _require(temporal.get("audit_pr") == 161, "temporal derivative audit provenance drifted")
    _require(temporal.get("snapshot_identity_not_derivative_identity") is True, "snapshot/derivative distinction was lost")
    _require(temporal.get("formal_pde_gate_assessed") is False, "temporal derivative audit is not the formal PDE gate")
    _require(temporal.get("pde_validated") is False, "temporal convergence cannot promote PDE validity")
    orders = [float(x) for x in temporal.get("observed_rms_orders", [])]
    _require(len(orders) == 2 and min(orders) > 1.9, "registered temporal derivative convergence drifted")

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
    _require(validation.get("derivative_steps") == [0.02, 0.01, 0.005], "spatial derivative ladder drifted")
    thresholds = validation.get("thresholds")
    _require(isinstance(thresholds, Mapping), "CR001 validation thresholds are missing")
    _require(thresholds.get("divergence_max") == 1e-5, "divergence max threshold drifted")
    _require(thresholds.get("divergence_L2") == 1e-5, "divergence L2 threshold drifted")
    _require(thresholds.get("pde_residual_max") == 0.001, "PDE max threshold drifted")
    _require(thresholds.get("pde_residual_L2") == 0.001, "PDE L2 threshold drifted")

    current_state = contract.get("current_state")
    _require(current_state == _EXPECTED_STATE, "cubic current state drifted")
    _require(capsule.get("status") == _EXPECTED_STATE, "live cubic capsule truth state drifted")

    evidence = contract.get("evidence_scope")
    _require(isinstance(evidence, Mapping), "evidence_scope is missing")
    selection = evidence.get("target_free_cubic_selection")
    _require(isinstance(selection, Mapping), "target-free cubic selection evidence is missing")
    _require(selection.get("source_pr") == 158, "target-free selection source drifted")
    _require(selection.get("status") == "integrated_in_capsule_ancestry", "selection ancestry drifted")
    _require(selection.get("public_openai_reference_used_for_selection") is False, "target-free selection was relabeled")
    _require(selection.get("may_support_visualization_trial_selection") is True, "trial-selection evidence was disabled")
    _require(selection.get("may_support_public_visual_correspondence") is False, "target-free selection cannot prove public correspondence")
    _require(selection.get("may_support_pde_validation") is False, "target-free selection cannot prove PDE validity")

    derivative = evidence.get("public_output_temporal_derivative")
    _require(isinstance(derivative, Mapping), "public-output temporal derivative evidence is missing")
    _require(derivative.get("source_pr") == 161, "temporal derivative source drifted")
    _require(derivative.get("source_head") == temporal.get("audit_exact_head"), "temporal derivative head drifted")
    _require(derivative.get("probe_seed") == temporal.get("probe_seed") == 914509, "temporal derivative seed drifted")
    _require(derivative.get("probe_count") == temporal.get("probe_count") == 24, "temporal derivative probe count drifted")
    _require(derivative.get("second_order_steps") == temporal.get("second_order_steps"), "temporal derivative ladder drifted")
    _require(float(derivative.get("fourth_order_reference_step")) == float(temporal.get("fourth_order_reference_step")), "temporal derivative reference step drifted")
    _require(derivative.get("snapshot_identity_not_derivative_identity") is True, "snapshot/derivative evidence drifted")
    _require(derivative.get("formal_pde_gate_assessed") is False, "temporal derivative diagnostic is not formal PDE acceptance")
    _require(derivative.get("may_support_public_velocity_derivative_confidence") is True, "valid derivative convergence evidence was disabled")
    _require(derivative.get("may_support_pde_validation") is False, "derivative convergence alone cannot promote PDE validity")

    keyframes = evidence.get("static_return_keyframes")
    _require(isinstance(keyframes, Mapping), "static-return keyframe evidence is missing")
    _require(keyframes.get("times") == EXPECTED_STATIC_RETURN_TIMES, "governed keyframe times drifted")
    _require(keyframes.get("pointwise_velocity_equals_static_parent") is True, "static-return snapshot identity drifted")
    for key in (
        "may_inherit_parent_time_derivative_evidence",
        "may_inherit_parent_pde_acceptance",
        "may_inherit_parent_energy_acceptance",
        "may_inherit_parent_visual_acceptance",
        "may_inherit_parent_canonical_status",
    ):
        _require(keyframes.get(key) is False, f"snapshot identity cannot enable {key}")

    sibling = evidence.get("restricted_force_sibling")
    _require(isinstance(sibling, Mapping), "restricted-force sibling evidence is missing")
    _require(sibling.get("source_pr") == 160, "restricted-force sibling provenance drifted")
    _require(sibling.get("status") == "open_unconsumed_sibling_evidence", "force sibling must remain unconsumed")
    _require(sibling.get("numeric_summary_bound_here") is False, "unconsumed sibling numerics cannot be bound here")
    _require(sibling.get("formal_pde_gate_assessed") is False, "restricted-force diagnostic is not formal PDE acceptance")
    _require(sibling.get("may_promote_pde_validated") is False, "sibling evidence cannot promote PDE validity")

    energy = evidence.get("candidate_energy")
    _require(isinstance(energy, Mapping), "candidate energy scope is missing")
    _require(energy.get("status") == "pending_candidate_specific_revalidation", "candidate energy status drifted")
    _require(energy.get("may_inherit_parent_energy_acceptance") is False, "time-dependent child cannot inherit parent energy acceptance")
    _require(energy.get("failure_blocks_velocity_export") is False, "energy failure cannot block velocity export")

    requirements = contract.get("promotion_requirements")
    _require(isinstance(requirements, Mapping), "promotion requirements are missing")
    _require(set(requirements.get("visualization_ready", [])) == _REQUIRED_VIS_READY, "visualization-ready evidence contract drifted")
    _require(set(requirements.get("visual_correspondence_verified", [])) == _REQUIRED_VIS_CORR, "public visual correspondence evidence contract drifted")
    _require(set(requirements.get("pde_validated", [])) == _REQUIRED_PDE, "PDE validation evidence contract drifted")

    policy = contract.get("policy")
    _require(isinstance(policy, Mapping), "policy is missing")
    _require(_FALSE_POLICY_KEYS.issubset(policy), "required fail-closed policy keys are missing")
    for key in _FALSE_POLICY_KEYS:
        _require(policy.get(key) is False, f"forbidden promotion inference enabled: {key}")

    truth = contract.get("truth_boundary")
    _require(isinstance(truth, Mapping), "truth boundary is missing")
    for key in (
        "velocity_changed_by_this_gate",
        "temporal_schedule_changed_by_this_gate",
        "support_transform_changed",
        "forcing_family_changed",
        "thresholds_changed",
        "candidate_promoted",
    ):
        _require(truth.get(key) is False, f"governance gate must not change/promote: {key}")

    return {
        "contract_pass": True,
        "candidate_sha256": EXPECTED_CANDIDATE_SHA256,
        "base_supported_sha256": EXPECTED_BASE_SHA256,
        "velocity_export_ready": True,
        "visualization_candidate_only": True,
        "visualization_ready": False,
        "visual_correspondence_verified": False,
        "pde_validated": False,
        "paper_exact": False,
        "openai_field_identified": False,
        "static_return_times": list(EXPECTED_STATIC_RETURN_TIMES),
        "snapshot_identity_not_derivative_identity": True,
        "candidate_energy_status": energy["status"],
        "force_sibling_status": sibling["status"],
    }
