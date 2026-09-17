import copy
import json
from pathlib import Path

import pytest

from openai_ns_reconstruction.constrained_eq45_export_scope_governance import (
    audit_eq45_export_scope_contract,
    audit_eq45_export_scope_contract_file,
)


ROOT = Path(__file__).resolve().parents[1]
CONTRACT = ROOT / "configs" / "eq45_export_scope_contract.json"


def _load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _copy_inputs(tmp_path: Path, contract: dict) -> None:
    paths = [
        contract["package_unified_velocity"]["delivery_contract"],
        contract["state_semantics"]["delivery_state_contract"],
        contract["eq45_candidate"]["artifact"],
    ]
    for relative in paths:
        source = ROOT / relative
        target = tmp_path / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(source.read_text(encoding="utf-8"), encoding="utf-8")


def test_repository_eq45_export_scope_contract_is_truthful() -> None:
    result = audit_eq45_export_scope_contract_file(CONTRACT, repo_root=ROOT)

    assert result["contract_pass"] is True
    assert result["eq45_scope"] == "candidate_local"
    assert result["eq45_callable_serializable"] is True
    assert result["eq45_velocity_export_ready"] is True
    assert result["eq45_is_package_unified_default"] is False
    assert result["package_candidate_family"] == "coupled_velocity_v1"
    assert result["pde_validated"] is False
    assert result["paper_exact"] is False


def test_candidate_local_readiness_cannot_claim_package_default() -> None:
    contract = _load(CONTRACT)
    mutated = copy.deepcopy(contract)
    mutated["eq45_candidate"]["is_package_unified_default"] = True

    with pytest.raises(ValueError, match="package unified/default"):
        audit_eq45_export_scope_contract(mutated, repo_root=ROOT)


def test_candidate_local_readiness_cannot_promote_package_binding() -> None:
    contract = _load(CONTRACT)
    mutated = copy.deepcopy(contract)
    mutated["state_semantics"]["eq45_local_readiness_is_evidence_for_package_unified_binding"] = True

    with pytest.raises(ValueError, match="must not promote"):
        audit_eq45_export_scope_contract(mutated, repo_root=ROOT)


def test_pde_failure_does_not_block_candidate_local_export() -> None:
    contract = _load(CONTRACT)
    mutated = copy.deepcopy(contract)
    mutated["state_semantics"]["pde_failure_blocks_candidate_local_export"] = True

    with pytest.raises(ValueError, match="must not block"):
        audit_eq45_export_scope_contract(mutated, repo_root=ROOT)


def test_eq45_scientific_claim_promotion_fails_closed(tmp_path: Path) -> None:
    contract = _load(CONTRACT)
    _copy_inputs(tmp_path, contract)
    candidate_path = tmp_path / contract["eq45_candidate"]["artifact"]
    candidate = _load(candidate_path)
    candidate["truth_boundary"]["pde_validated"] = True
    candidate_path.write_text(json.dumps(candidate), encoding="utf-8")

    with pytest.raises(ValueError, match="unsupported Eq45 claim promotion: pde_validated"):
        audit_eq45_export_scope_contract(contract, repo_root=tmp_path)
