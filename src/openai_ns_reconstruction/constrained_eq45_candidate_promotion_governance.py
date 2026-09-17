"""Fail-closed governance for promoting an optimized Eq45 candidate.

This module separates a useful callable/exportable research artifact from a
canonical visualization-candidate decision.  Optimizer convergence, residual
reduction, serialization, and green CI are not visual-correspondence evidence.
PDE validation remains an independent claim and is not a prerequisite for
saving or rendering a visualization candidate.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Mapping

from .constrained_eq45_candidate import Eq45VelocityCandidate

SCHEMA = "eq45_candidate_promotion_gate_v1"
TASK_ID = "CR002-EQ45-CANONICAL-PROMOTION-GATE-017"
CANONICAL_VOCABULARY = {
    "user_requirement",
    "public_source_fact",
    "autonomous_design",
    "pending_unknown",
}
REQUIRED_FALSE_CLAIMS = {
    "physical_support_validated",
    "visualization_ready",
    "visual_correspondence_verified",
    "pde_validated",
    "paper_exact",
    "openai_field_identified",
    "blowup_proved",
}
REQUIRED_PROMOTION_EVIDENCE = [
    "independent_public_visual_comparison",
    "independent_visual_resolution_stability",
    "explicit_truth_boundary_review",
]


def _load(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"expected JSON object: {path}")
    return value


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def audit_eq45_candidate_promotion(
    contract: Mapping[str, Any], *, repo_root: str | Path
) -> dict[str, Any]:
    """Audit the current optimized Eq45 artifact's promotion status."""
    root = Path(repo_root)
    _require(contract.get("schema") == SCHEMA, "unsupported promotion-gate schema")
    _require(contract.get("task_id") == TASK_ID, "unexpected promotion-gate task id")

    vocab_path = root / str(contract.get("classification_vocabulary_contract", ""))
    _require(vocab_path.is_file(), "classification vocabulary contract is missing")
    delivery_states = _load(vocab_path)
    vocabulary = delivery_states.get("classification_vocabulary")
    _require(isinstance(vocabulary, Mapping), "delivery-state vocabulary is missing")
    _require(set(vocabulary) == CANONICAL_VOCABULARY, "canonical source vocabulary drifted")

    classifications = contract.get("source_classification")
    _require(isinstance(classifications, Mapping), "source_classification is missing")
    _require(
        set(classifications.values()).issubset(CANONICAL_VOCABULARY),
        "promotion contract uses a noncanonical source classification",
    )
    _require(
        classifications.get("optimized_candidate_matches_public_visualization") == "pending_unknown",
        "visual correspondence must remain pending without public comparison evidence",
    )
    _require(
        classifications.get("optimized_candidate_is_exact_openai_field") == "pending_unknown",
        "OpenAI-field identity must remain pending without public identity evidence",
    )

    inputs = contract.get("inputs")
    evidence = contract.get("current_evidence")
    policy = contract.get("promotion_policy")
    _require(isinstance(inputs, Mapping), "promotion inputs are missing")
    _require(isinstance(evidence, Mapping), "current promotion evidence is missing")
    _require(isinstance(policy, Mapping), "promotion policy is missing")

    seed_path = root / str(inputs.get("canonical_seed_artifact", ""))
    optimized_path = root / str(inputs.get("optimized_artifact", ""))
    receipt_path = root / str(inputs.get("optimization_receipt", ""))
    _require(seed_path.is_file(), "canonical Eq45 seed is missing")
    _require(optimized_path.is_file(), "optimized Eq45 artifact is missing")
    _require(receipt_path.is_file(), "optimization receipt is missing")

    seed = Eq45VelocityCandidate.load_json(seed_path)
    optimized = Eq45VelocityCandidate.load_json(optimized_path)
    receipt = _load(receipt_path)
    _require(seed.sha256 == inputs.get("canonical_seed_sha256"), "canonical seed SHA drifted")
    _require(optimized.sha256 == inputs.get("optimized_sha256"), "optimized candidate SHA drifted")
    _require(seed.sha256 != optimized.sha256, "optimized artifact unexpectedly equals the seed")

    fit = receipt.get("fit")
    holdout = receipt.get("holdout")
    optimized_meta = receipt.get("optimized_candidate")
    receipt_status = receipt.get("status")
    velocity_change = receipt.get("velocity_change_on_holdout_probes")
    _require(isinstance(fit, Mapping), "optimization fit section is missing")
    _require(isinstance(holdout, list) and holdout, "optimization holdout evidence is missing")
    _require(isinstance(optimized_meta, Mapping), "optimized candidate receipt metadata is missing")
    _require(isinstance(receipt_status, Mapping), "optimization status is missing")
    _require(isinstance(velocity_change, Mapping), "velocity-change receipt is missing")

    _require(fit.get("optimizer_success") is True, "optimizer convergence receipt drifted")
    _require(optimized_meta.get("canonical_seed_promoted") is False, "optimizer receipt promoted the canonical seed")
    _require(optimized_meta.get("sha256") == optimized.sha256, "optimizer/artifact identity mismatch")
    _require(optimized_meta.get("roundtrip_serializable") is True, "optimized candidate lost round-trip evidence")
    finest = holdout[-1]
    _require(isinstance(finest, Mapping), "finest holdout row is malformed")
    _require(float(finest.get("improvement_vs_seed", 0.0)) > 0.0, "held-out obstruction did not improve")

    recorded_delta = float(velocity_change.get("relative_delta_rms"))
    _require(recorded_delta > 0.0, "optimized velocity did not materially differ from the seed")
    _require(
        abs(recorded_delta - float(evidence.get("optimized_velocity_changed_relative_rms"))) <= 1e-15,
        "promotion contract velocity-change evidence drifted from optimization receipt",
    )
    _require(evidence.get("optimizer_converged") is True, "optimizer evidence drifted")
    _require(evidence.get("heldout_pressure_curl_obstruction_improved") is True, "holdout evidence drifted")
    _require(evidence.get("artifact_roundtrip_serializable") is True, "serialization evidence drifted")
    _require(evidence.get("candidate_local_velocity_export_ready") is True, "candidate export readiness drifted")

    truth = optimized.to_dict().get("truth_boundary")
    _require(isinstance(truth, Mapping), "optimized artifact truth boundary is missing")
    required_false = contract.get("required_false_scientific_claims")
    _require(isinstance(required_false, list), "required false claim set is missing")
    _require(set(required_false) == REQUIRED_FALSE_CLAIMS, "required false scientific claims drifted")
    for key in REQUIRED_FALSE_CLAIMS:
        _require(truth.get(key) is False, f"optimized artifact illegally promotes {key}")
    _require(truth.get("velocity_export_ready") is True, "optimized artifact lost velocity export readiness")
    _require(truth.get("callable_serializable") is True, "optimized artifact lost callable/serialization readiness")

    for key in ("pde_validated", "visual_correspondence_verified", "paper_exact", "openai_field_identified"):
        _require(receipt_status.get(key) is False, f"optimization receipt illegally promotes {key}")
    _require(receipt_status.get("visualization_ready_promoted") is False, "optimizer receipt promoted visualization readiness")

    _require(policy.get("diagnostic_candidate_may_be_saved_loaded_exported_and_rendered") is True, "diagnostic use was incorrectly blocked")
    _require(policy.get("pde_validation_is_required_for_callable_export") is False, "PDE validation cannot gate callable export")
    _require(policy.get("pde_validation_is_required_for_visualization_candidate_status") is False, "PDE validation cannot gate visualization-candidate status")
    _require(policy.get("pde_validation_is_required_for_pde_validated_status") is True, "PDE status lost its independent evidence gate")
    _require(policy.get("optimizer_convergence_is_promotion_evidence") is False, "optimizer convergence cannot promote a candidate")
    _require(policy.get("residual_reduction_is_visual_correspondence_evidence") is False, "residual reduction cannot establish visual correspondence")
    _require(policy.get("serialization_or_ci_green_is_visual_correspondence_evidence") is False, "serialization/CI cannot establish visual correspondence")
    _require(policy.get("required_before_canonical_visual_promotion") == REQUIRED_PROMOTION_EVIDENCE, "canonical visual-promotion evidence set drifted")

    evidence_ready = all(evidence.get(key) is True for key in REQUIRED_PROMOTION_EVIDENCE)
    promoted = policy.get("canonical_visual_candidate_promoted")
    _require(isinstance(promoted, bool), "canonical promotion state must be boolean")
    if promoted and not evidence_ready:
        raise ValueError("canonical visual promotion requires independent visual evidence")
    _require(promoted is False, "current optimized Eq45 artifact must remain non-canonical")
    _require(evidence.get("independent_public_visual_comparison") == "pending_unknown", "public visual comparison is not yet established")
    _require(evidence.get("independent_visual_resolution_stability") == "pending_unknown", "optimized visual-resolution stability is not yet established")
    _require(evidence.get("explicit_truth_boundary_review") is True, "truth-boundary review evidence is missing")

    return {
        "contract_pass": True,
        "canonical_seed_sha256": seed.sha256,
        "optimized_candidate_sha256": optimized.sha256,
        "optimized_velocity_relative_rms_change": recorded_delta,
        "candidate_local_velocity_export_ready": True,
        "diagnostic_rendering_allowed": True,
        "canonical_visual_candidate_promoted": False,
        "canonical_visual_promotion_eligible": evidence_ready,
        "visual_correspondence_verified": False,
        "pde_validated": False,
        "paper_exact": False,
        "openai_field_identified": False,
    }


def audit_eq45_candidate_promotion_file(
    path: str | Path, *, repo_root: str | Path
) -> dict[str, Any]:
    return audit_eq45_candidate_promotion(_load(Path(path)), repo_root=repo_root)
