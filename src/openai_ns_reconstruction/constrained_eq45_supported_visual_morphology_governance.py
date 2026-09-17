"""Fail-closed visualization promotion governance for the supported Eq45 child.

A physical support connection is a non-identity change to ``[u,v,w]``.  Its
exact zero extension, inner identity plateau, or green CI therefore cannot make
that transformed child visualization-ready by themselves.  This auditor keeps
engineering delivery, visual morphology evidence, public-reference
correspondence, and PDE validation as independent states.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Mapping

from .eq45_supported_delivery import default_field


SCHEMA = "eq45_supported_visual_morphology_promotion_gate_v1"
_EXPECTED_CLASSIFICATION = {
    "velocity_visual_delivery": "user_requirement",
    "eq45_backbone": "public_source_fact",
    "physical_support_transform": "autonomous_design",
    "support_transform_visual_neutrality": "pending_unknown",
    "public_visual_correspondence": "pending_unknown",
}
_EXPECTED_CHILD_SHA = "2fdff812c22131d56eac1b7d6e455187ff3207500c6385aa9abf282a8e3d7b1d"
_EXPECTED_PARENT_SHA = "48f1845fcfd71ec95495c98bb7aac3fca4653e748a8856a5421234e9601525a7"
_EXPECTED_CURRENT_STATE = {
    "callable_serializable": True,
    "velocity_export_ready": True,
    "physical_support_connection_implemented": True,
    "physical_support_validated": False,
    "visualization_ready": False,
    "visual_correspondence_verified": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}
_REQUIRED_VISUALIZATION_READY = {
    "callable_velocity",
    "supported_visualization_entry_point",
    "candidate_specific_resolution_sanity",
    "support_transform_morphology_review_resolved",
}
_REQUIRED_VISUAL_CORRESPONDENCE = {
    "public_reference_only",
    "candidate_specific_independent_visual_diagnostics",
    "candidate_specific_resolution_stability",
    "support_transform_morphology_review_resolved",
}
_FALSE_POLICY_KEYS = {
    "support_connection_implies_visualization_ready",
    "exact_zero_outside_support_implies_visualization_ready",
    "identity_plateau_implies_whole_domain_morphology_neutral",
    "pointwise_collar_stability_implies_whole_domain_morphology_neutral",
    "resolution_stability_implies_visual_correspondence",
    "visual_similarity_implies_pde_validation",
    "visual_correspondence_implies_paper_exact",
    "pde_failure_blocks_velocity_export",
    "open_sibling_visual_evidence_is_integrated_acceptance",
    "unresolved_support_morphology_allows_visualization_ready_promotion",
    "unresolved_public_comparison_allows_visual_correspondence_promotion",
}


def _read_json(path: Path) -> dict[str, Any]:
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise ValueError(f"{path} must contain a JSON object")
    return data


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def audit_eq45_supported_visual_morphology_gate(
    contract: Mapping[str, Any], *, repo_root: str | Path
) -> dict[str, Any]:
    """Audit support-transform visual-state promotion without rendering.

    Passing this audit means only that the repository keeps the state/evidence
    boundaries truthful.  It is not positive visual correspondence, PDE, or
    paper-exact evidence.
    """

    root = Path(repo_root)
    _require(contract.get("schema") == SCHEMA, "unexpected visual morphology gate schema")

    delivery_state = _read_json(root / "configs" / "delivery_state_contract.json")
    constraints = _read_json(root / "configs" / "constraints.json")
    vocabulary = delivery_state.get("classification_vocabulary")
    _require(isinstance(vocabulary, Mapping), "delivery state vocabulary is missing")

    classification = contract.get("classification")
    _require(classification == _EXPECTED_CLASSIFICATION, "visual morphology source classification drifted")
    _require(
        set(classification.values()).issubset(vocabulary),
        "visual morphology gate uses noncanonical source classification",
    )

    field = default_field()
    metadata = field.metadata()
    candidate = contract.get("candidate")
    _require(isinstance(candidate, Mapping), "candidate section is missing")
    _require(candidate.get("family") == metadata.get("family"), "supported candidate family drifted")
    _require(candidate.get("child_sha256") == field.sha256 == _EXPECTED_CHILD_SHA, "supported child identity drifted")
    _require(
        candidate.get("parent_sha256") == metadata.get("parent_sha256") == _EXPECTED_PARENT_SHA,
        "supported parent identity drifted",
    )
    _require(
        candidate.get("public_velocity_api") == metadata.get("entrypoint"),
        "supported public velocity API drifted",
    )

    domain = constraints.get("domain")
    _require(isinstance(domain, Mapping), "constraints.domain is missing")
    _require(candidate.get("physical_support") == domain.get("support"), "physical support drifted from CR001")
    _require(candidate.get("time_interval") == domain.get("time_interval"), "time interval drifted from CR001")
    _require(metadata.get("time_interval") == domain.get("time_interval"), "live supported field time interval drifted")

    current_state = contract.get("current_state")
    _require(current_state == _EXPECTED_CURRENT_STATE, "supported visual current-state contract drifted")
    truth = metadata.get("truth_boundary")
    _require(isinstance(truth, Mapping), "supported field truth boundary is missing")
    for key, expected in _EXPECTED_CURRENT_STATE.items():
        _require(truth.get(key) is expected, f"live supported field truth state drifted: {key}")

    evidence = contract.get("current_visual_evidence")
    _require(isinstance(evidence, Mapping), "current_visual_evidence is missing")
    _require(
        evidence.get("support_transform_morphology_review") == "pending_unresolved",
        "support-transform morphology review must stay unresolved until candidate-specific evidence is integrated",
    )
    _require(
        evidence.get("candidate_specific_resolution_sanity") == "pending_unresolved",
        "candidate-specific visual resolution sanity is not integrated in this ancestry",
    )
    _require(
        evidence.get("public_reference_comparison") == "pending_unresolved",
        "public-reference comparison is not integrated in this ancestry",
    )
    sibling = evidence.get("external_sibling_signal")
    _require(isinstance(sibling, Mapping), "external sibling visual signal is missing")
    _require(sibling.get("source_pr") == 139, "external sibling signal no longer points to PR #139")
    _require(
        sibling.get("status") == "open_unconsumed_sibling_evidence",
        "open sibling evidence must remain explicitly unconsumed in this ancestry",
    )
    _require(sibling.get("may_trigger_followup_review") is True, "sibling diagnostic may guide follow-up review")
    _require(
        sibling.get("may_promote_current_ancestry_state") is False,
        "open sibling evidence cannot promote current-ancestry state",
    )

    requirements = contract.get("candidate_specific_promotion_requirements")
    _require(isinstance(requirements, Mapping), "candidate-specific promotion requirements are missing")
    _require(
        set(requirements.get("visualization_ready", [])) == _REQUIRED_VISUALIZATION_READY,
        "supported-child visualization-ready requirements drifted",
    )
    _require(
        set(requirements.get("visual_correspondence_verified", [])) == _REQUIRED_VISUAL_CORRESPONDENCE,
        "supported-child visual-correspondence requirements drifted",
    )

    states = delivery_state.get("states")
    _require(isinstance(states, Mapping), "delivery state definitions are missing")
    global_vis = states.get("visualization_ready")
    global_corr = states.get("visual_correspondence_verified")
    _require(isinstance(global_vis, Mapping), "global visualization_ready definition is missing")
    _require(isinstance(global_corr, Mapping), "global visual_correspondence definition is missing")
    _require(
        {"callable_velocity", "supported_visualization_entry_point", "resolution_sanity_check"}.issubset(
            set(global_vis.get("positive_evidence", []))
        ),
        "candidate gate cannot weaken global visualization-ready evidence",
    )
    _require(
        {"public_reference_only", "independent_visual_diagnostics", "resolution_stability"}.issubset(
            set(global_corr.get("positive_evidence", []))
        ),
        "candidate gate cannot weaken global visual-correspondence evidence",
    )

    policy = contract.get("policy")
    _require(isinstance(policy, Mapping), "visual morphology policy is missing")
    _require(set(policy) == _FALSE_POLICY_KEYS, "visual morphology policy key set drifted")
    for key in _FALSE_POLICY_KEYS:
        _require(policy.get(key) is False, f"forbidden visual promotion inference enabled: {key}")

    boundary = contract.get("truth_boundary")
    _require(isinstance(boundary, Mapping), "governance truth boundary is missing")
    for key in ("velocity_changed", "taper_changed", "forcing_changed", "thresholds_changed", "candidate_promoted"):
        _require(boundary.get(key) is False, f"governance-only increment drifted: {key}")

    return {
        "contract_pass": True,
        "candidate_sha256": field.sha256,
        "velocity_export_ready": True,
        "physical_support_connection_implemented": True,
        "support_transform_morphology_review": "pending_unresolved",
        "visualization_ready": False,
        "visual_correspondence_verified": False,
        "pde_validated": False,
        "paper_exact": False,
        "openai_field_identified": False,
        "next_visual_gate": "candidate_specific_support_transform_morphology_review",
    }


def audit_eq45_supported_visual_morphology_gate_file(
    path: str | Path, *, repo_root: str | Path
) -> dict[str, Any]:
    return audit_eq45_supported_visual_morphology_gate(
        _read_json(Path(path)), repo_root=repo_root
    )
