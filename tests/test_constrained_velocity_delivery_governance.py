from __future__ import annotations

import copy
import json
from pathlib import Path

import pytest

from openai_ns_reconstruction.constrained_velocity_delivery_governance import (
    audit_velocity_delivery_contract,
    audit_velocity_delivery_contract_file,
)


ROOT = Path(".")


def _contract() -> dict:
    return json.loads(Path("configs/velocity_delivery_contract.json").read_text(encoding="utf-8"))


def test_repository_velocity_delivery_contract_is_truthful() -> None:
    result = audit_velocity_delivery_contract_file(
        "configs/velocity_delivery_contract.json", repo_root=ROOT
    )
    assert result["contract_pass"] is True
    assert result["velocity_export_ready"] is True
    assert result["visualization_ready"] is False
    assert result["visual_correspondence_verified"] is False
    assert result["pde_validated"] is False
    assert result["candidate_family"] == "eq45_supported_velocity_candidate_v1"
    assert result["velocity_api"] == "openai_ns_reconstruction.eq45_supported_delivery:velocity"
    assert set(result["source_classes"]) == {
        "autonomous_design",
        "pending_unknown",
        "public_source_fact",
        "user_requirement",
    }


def test_visual_correspondence_cannot_be_promoted_without_evidence() -> None:
    mutated = copy.deepcopy(_contract())
    mutated["claim_status"]["visual_correspondence_verified"] = True
    with pytest.raises(ValueError, match="time mapping|unsupported claim promotion"):
        audit_velocity_delivery_contract(mutated, repo_root=ROOT)


def test_pde_and_paper_claims_fail_closed() -> None:
    mutated = copy.deepcopy(_contract())
    mutated["claim_status"]["pde_validated"] = True
    with pytest.raises(ValueError, match="PDE-valid status lacks independent evidence"):
        audit_velocity_delivery_contract(mutated, repo_root=ROOT)

    mutated = copy.deepcopy(_contract())
    mutated["claim_status"]["paper_exact"] = True
    with pytest.raises(ValueError, match="unsupported claim promotion: paper_exact"):
        audit_velocity_delivery_contract(mutated, repo_root=ROOT)


def test_domain_or_source_classification_drift_fails_closed() -> None:
    mutated = copy.deepcopy(_contract())
    mutated["primary_deliverable"]["time_interval"] = [0.0, 1.0]
    with pytest.raises(ValueError, match="time interval drifted"):
        audit_velocity_delivery_contract(mutated, repo_root=ROOT)

    mutated = copy.deepcopy(_contract())
    mutated["source_classification"][0]["classification"] = "public_source_fact"
    with pytest.raises(ValueError, match="coverage is incomplete"):
        audit_velocity_delivery_contract(mutated, repo_root=ROOT)


def test_legacy_surface_cannot_be_relabelled_canonical() -> None:
    mutated = copy.deepcopy(_contract())
    mutated["legacy_compatibility_delivery"]["canonical"] = True
    with pytest.raises(ValueError, match="legacy delivery cannot be canonical"):
        audit_velocity_delivery_contract(mutated, repo_root=ROOT)

    mutated = copy.deepcopy(_contract())
    mutated["primary_deliverable"]["api"] = "openai_ns_reconstruction.velocity_components:velocity"
    with pytest.raises(ValueError, match="unexpected primary velocity API"):
        audit_velocity_delivery_contract(mutated, repo_root=ROOT)
