from __future__ import annotations

import copy
import json
from pathlib import Path

import pytest

from openai_ns_reconstruction.constrained_velocity_delivery_governance import (
    audit_velocity_delivery_contract,
    audit_velocity_delivery_contract_file,
)


def _fixture_repo(tmp_path: Path) -> tuple[Path, dict]:
    (tmp_path / "configs").mkdir()
    (tmp_path / "docs").mkdir()
    (tmp_path / "src/openai_ns_reconstruction/data").mkdir(parents=True)
    for name in ("PROJECT_GOAL.md", "VISUAL_TARGET.md", "VELOCITY_API.md"):
        (tmp_path / "docs" / name).write_text("fixture\n", encoding="utf-8")
    (tmp_path / "src/openai_ns_reconstruction/data/velocity_candidate.json").write_text(
        json.dumps({"family": "coupled_velocity_v1", "status": "candidate"}), encoding="utf-8"
    )
    (tmp_path / "configs/constraints_coupled.json").write_text(
        json.dumps({"domain": {"time_interval": [0.25, 0.75], "support": "r < 2 and abs(z) < 2"}}),
        encoding="utf-8",
    )
    contract = json.loads(Path("configs/velocity_delivery_contract.json").read_text(encoding="utf-8"))
    return tmp_path, contract


def test_repository_velocity_delivery_contract_is_truthful() -> None:
    result = audit_velocity_delivery_contract_file(
        "configs/velocity_delivery_contract.json", repo_root=Path(".")
    )
    assert result["contract_pass"] is True
    assert result["velocity_export_ready"] is True
    assert result["visual_correspondence_verified"] is False
    assert result["pde_validated"] is False
    assert result["candidate_family"] == "coupled_velocity_v1"
    assert set(result["source_classes"]) == {
        "autonomous_design",
        "pending_unknown",
        "public_source_fact",
        "user_requirement",
    }


def test_visual_correspondence_cannot_be_promoted_without_evidence(tmp_path: Path) -> None:
    root, contract = _fixture_repo(tmp_path)
    mutated = copy.deepcopy(contract)
    mutated["claim_status"]["openai_correspondence_verified"] = True
    with pytest.raises(ValueError, match="time mapping|unsupported claim promotion"):
        audit_velocity_delivery_contract(mutated, repo_root=root)


def test_pde_and_paper_claims_fail_closed(tmp_path: Path) -> None:
    root, contract = _fixture_repo(tmp_path)
    mutated = copy.deepcopy(contract)
    mutated["claim_status"]["pde_validated"] = True
    with pytest.raises(ValueError, match="PDE-valid status lacks independent evidence"):
        audit_velocity_delivery_contract(mutated, repo_root=root)

    mutated = copy.deepcopy(contract)
    mutated["claim_status"]["paper_exact"] = True
    with pytest.raises(ValueError, match="unsupported claim promotion: paper_exact"):
        audit_velocity_delivery_contract(mutated, repo_root=root)


def test_domain_or_source_classification_drift_fails_closed(tmp_path: Path) -> None:
    root, contract = _fixture_repo(tmp_path)
    mutated = copy.deepcopy(contract)
    mutated["primary_deliverable"]["time_interval"] = [0.0, 1.0]
    with pytest.raises(ValueError, match="time interval drifted"):
        audit_velocity_delivery_contract(mutated, repo_root=root)

    mutated = copy.deepcopy(contract)
    mutated["source_classification"][0]["classification"] = "public_source"
    with pytest.raises(ValueError, match="invalid source classification"):
        audit_velocity_delivery_contract(mutated, repo_root=root)
