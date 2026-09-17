import copy
import json
from pathlib import Path

import pytest

from openai_ns_reconstruction.constrained_eq45_supported_visual_morphology_governance import (
    audit_eq45_supported_visual_morphology_gate,
    audit_eq45_supported_visual_morphology_gate_file,
)


ROOT = Path(__file__).resolve().parents[1]
CONTRACT = ROOT / "configs" / "eq45_supported_visual_morphology_gate.json"


def _contract() -> dict:
    return json.loads(CONTRACT.read_text(encoding="utf-8"))


def test_repository_supported_visual_morphology_gate_is_truthful() -> None:
    result = audit_eq45_supported_visual_morphology_gate_file(
        CONTRACT, repo_root=ROOT
    )

    assert result["contract_pass"] is True
    assert result["velocity_export_ready"] is True
    assert result["physical_support_connection_implemented"] is True
    assert result["support_transform_morphology_review"] == "pending_unresolved"
    assert result["visualization_ready"] is False
    assert result["visual_correspondence_verified"] is False
    assert result["pde_validated"] is False
    assert result["paper_exact"] is False


def test_support_connection_cannot_promote_visualization_ready() -> None:
    mutated = copy.deepcopy(_contract())
    mutated["policy"]["support_connection_implies_visualization_ready"] = True

    with pytest.raises(ValueError, match="forbidden visual promotion inference"):
        audit_eq45_supported_visual_morphology_gate(mutated, repo_root=ROOT)


def test_identity_plateau_cannot_claim_whole_domain_visual_neutrality() -> None:
    mutated = copy.deepcopy(_contract())
    mutated["policy"]["identity_plateau_implies_whole_domain_morphology_neutral"] = True

    with pytest.raises(ValueError, match="forbidden visual promotion inference"):
        audit_eq45_supported_visual_morphology_gate(mutated, repo_root=ROOT)


def test_pointwise_collar_stability_cannot_claim_whole_domain_neutrality() -> None:
    mutated = copy.deepcopy(_contract())
    mutated["policy"]["pointwise_collar_stability_implies_whole_domain_morphology_neutral"] = True

    with pytest.raises(ValueError, match="forbidden visual promotion inference"):
        audit_eq45_supported_visual_morphology_gate(mutated, repo_root=ROOT)


def test_unresolved_morphology_cannot_be_relabelled_as_visualization_ready() -> None:
    mutated = copy.deepcopy(_contract())
    mutated["current_state"]["visualization_ready"] = True

    with pytest.raises(ValueError, match="current-state contract drifted"):
        audit_eq45_supported_visual_morphology_gate(mutated, repo_root=ROOT)


def test_unconsumed_sibling_evidence_cannot_promote_current_ancestry() -> None:
    mutated = copy.deepcopy(_contract())
    mutated["current_visual_evidence"]["external_sibling_signal"][
        "may_promote_current_ancestry_state"
    ] = True

    with pytest.raises(ValueError, match="cannot promote current-ancestry state"):
        audit_eq45_supported_visual_morphology_gate(mutated, repo_root=ROOT)


def test_visual_similarity_cannot_promote_pde_validation() -> None:
    mutated = copy.deepcopy(_contract())
    mutated["policy"]["visual_similarity_implies_pde_validation"] = True

    with pytest.raises(ValueError, match="forbidden visual promotion inference"):
        audit_eq45_supported_visual_morphology_gate(mutated, repo_root=ROOT)


def test_pde_failure_cannot_block_velocity_export() -> None:
    mutated = copy.deepcopy(_contract())
    mutated["policy"]["pde_failure_blocks_velocity_export"] = True

    with pytest.raises(ValueError, match="forbidden visual promotion inference"):
        audit_eq45_supported_visual_morphology_gate(mutated, repo_root=ROOT)


def test_support_transform_cannot_be_reclassified_as_public_source_fact() -> None:
    mutated = copy.deepcopy(_contract())
    mutated["classification"]["physical_support_transform"] = "public_source_fact"

    with pytest.raises(ValueError, match="source classification drifted"):
        audit_eq45_supported_visual_morphology_gate(mutated, repo_root=ROOT)
