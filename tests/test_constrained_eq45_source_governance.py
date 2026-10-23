from __future__ import annotations

import copy
import json
from pathlib import Path

import pytest

from openai_ns_reconstruction.constrained_eq45_source_governance import (
    audit_eq45_source_contract,
    audit_eq45_source_contract_file,
)


def _fixture_repo(tmp_path: Path) -> tuple[Path, dict]:
    (tmp_path / "docs").mkdir()
    for name in ("VELOCITY_FORMULAS.md", "PROJECT_GOAL.md"):
        (tmp_path / "docs" / name).write_text("fixture\n", encoding="utf-8")
    contract = json.loads(Path("configs/eq45_source_contract.json").read_text(encoding="utf-8"))
    return tmp_path, contract


def test_repository_eq45_source_contract_is_truthful() -> None:
    result = audit_eq45_source_contract_file(
        "configs/eq45_source_contract.json", repo_root=Path(".")
    )
    assert result["contract_pass"] is True
    assert result["q_relation"] == "q-z^2*q^(2*h)=1-t"
    assert result["q_exponent_sign"] == "positive"
    assert result["coordinate_h_open_interval"] == [0.0, 0.5]
    assert result["leading_field_only"] is True
    assert result["pde_validated"] is False
    assert result["paper_exact"] is False
    assert set(result["source_classes"]) == {
        "autonomous_design",
        "pending_unknown",
        "public_source_fact",
        "user_requirement",
    }


def test_negative_q_exponent_drift_fails_closed(tmp_path: Path) -> None:
    root, contract = _fixture_repo(tmp_path)
    mutated = copy.deepcopy(contract)
    mutated["coordinate_contract"]["implicit_q_relation"] = "q-z^2*q^(-2*h)=1-t"
    mutated["coordinate_contract"]["q_exponent_sign"] = "negative"
    with pytest.raises(ValueError, match="coordinate drift: implicit_q_relation|coordinate drift: q_exponent_sign"):
        audit_eq45_source_contract(mutated, repo_root=root)


def test_source_domain_drift_fails_closed(tmp_path: Path) -> None:
    root, contract = _fixture_repo(tmp_path)
    mutated = copy.deepcopy(contract)
    mutated["coordinate_contract"]["h_open_interval"] = [0.0, 1.0]
    with pytest.raises(ValueError, match="h-domain drift"):
        audit_eq45_source_contract(mutated, repo_root=root)


def test_leading_formula_cannot_promote_stronger_claims(tmp_path: Path) -> None:
    root, contract = _fixture_repo(tmp_path)
    mutated = copy.deepcopy(contract)
    mutated["claim_status"]["paper_exact"] = True
    with pytest.raises(ValueError, match="unsupported claim promotion: paper_exact"):
        audit_eq45_source_contract(mutated, repo_root=root)

    mutated = copy.deepcopy(contract)
    mutated["velocity_contract"]["leading_field_only"] = False
    with pytest.raises(ValueError, match="leading-field only"):
        audit_eq45_source_contract(mutated, repo_root=root)


def test_legacy_source_class_alias_fails_closed(tmp_path: Path) -> None:
    root, contract = _fixture_repo(tmp_path)
    mutated = copy.deepcopy(contract)
    mutated["source_classification"][0]["classification"] = "public_source"
    with pytest.raises(ValueError, match="invalid source classification"):
        audit_eq45_source_contract(mutated, repo_root=root)


def test_pending_profile_data_cannot_be_reclassified_as_public_fact(tmp_path: Path) -> None:
    root, contract = _fixture_repo(tmp_path)
    mutated = copy.deepcopy(contract)
    mutated["source_classification"][3]["classification"] = "public_source_fact"
    with pytest.raises(ValueError, match="source classification coverage is incomplete|source classification drift"):
        audit_eq45_source_contract(mutated, repo_root=root)
