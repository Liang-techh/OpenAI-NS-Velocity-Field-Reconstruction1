import copy
import json
from pathlib import Path

import pytest

from openai_ns_reconstruction.constrained_eq45_candidate_promotion_governance import (
    audit_eq45_candidate_promotion,
    audit_eq45_candidate_promotion_file,
)


ROOT = Path(__file__).resolve().parents[1]
CONTRACT = ROOT / "configs" / "eq45_candidate_promotion_gate.json"


def _contract() -> dict:
    return json.loads(CONTRACT.read_text(encoding="utf-8"))


def test_repository_optimized_candidate_remains_truth_bounded_and_exportable() -> None:
    result = audit_eq45_candidate_promotion_file(CONTRACT, repo_root=ROOT)

    assert result["contract_pass"] is True
    assert result["candidate_local_velocity_export_ready"] is True
    assert result["diagnostic_rendering_allowed"] is True
    assert result["optimized_velocity_relative_rms_change"] == pytest.approx(0.200973334054043)
    assert result["canonical_visual_candidate_promoted"] is False
    assert result["canonical_visual_promotion_eligible"] is False
    assert result["visual_correspondence_verified"] is False
    assert result["pde_validated"] is False
    assert result["paper_exact"] is False
    assert result["openai_field_identified"] is False


def test_optimizer_or_residual_success_cannot_promote_without_visual_evidence() -> None:
    mutated = copy.deepcopy(_contract())
    mutated["promotion_policy"]["canonical_visual_candidate_promoted"] = True

    with pytest.raises(ValueError, match="requires independent visual evidence"):
        audit_eq45_candidate_promotion(mutated, repo_root=ROOT)

    mutated = copy.deepcopy(_contract())
    mutated["promotion_policy"]["residual_reduction_is_visual_correspondence_evidence"] = True
    with pytest.raises(ValueError, match="residual reduction cannot establish visual correspondence"):
        audit_eq45_candidate_promotion(mutated, repo_root=ROOT)


def test_pde_failure_cannot_become_export_or_visualization_candidate_blocker() -> None:
    mutated = copy.deepcopy(_contract())
    mutated["promotion_policy"]["pde_validation_is_required_for_callable_export"] = True
    with pytest.raises(ValueError, match="cannot gate callable export"):
        audit_eq45_candidate_promotion(mutated, repo_root=ROOT)

    mutated = copy.deepcopy(_contract())
    mutated["promotion_policy"]["pde_validation_is_required_for_visualization_candidate_status"] = True
    with pytest.raises(ValueError, match="cannot gate visualization-candidate status"):
        audit_eq45_candidate_promotion(mutated, repo_root=ROOT)


def test_pending_visual_or_openai_identity_cannot_be_reclassified_as_known() -> None:
    mutated = copy.deepcopy(_contract())
    mutated["source_classification"]["optimized_candidate_matches_public_visualization"] = "public_source_fact"
    with pytest.raises(ValueError, match="visual correspondence must remain pending"):
        audit_eq45_candidate_promotion(mutated, repo_root=ROOT)

    mutated = copy.deepcopy(_contract())
    mutated["source_classification"]["optimized_candidate_is_exact_openai_field"] = "autonomous_design"
    with pytest.raises(ValueError, match="OpenAI-field identity must remain pending"):
        audit_eq45_candidate_promotion(mutated, repo_root=ROOT)
