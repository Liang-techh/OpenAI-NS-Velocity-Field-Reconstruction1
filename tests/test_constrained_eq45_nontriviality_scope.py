import copy
import json
from pathlib import Path

import pytest

from openai_ns_reconstruction.constrained_eq45_nontriviality_scope import (
    audit_eq45_nontriviality_scope,
    audit_eq45_nontriviality_scope_file,
)


ROOT = Path(__file__).resolve().parents[1]
CONTRACT = ROOT / "configs" / "eq45_nontriviality_scope_contract.json"


def _load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _copy_inputs(tmp_path: Path, contract: dict) -> None:
    paths = [
        contract["eq45_candidate"]["artifact"],
        contract["cr001_nontriviality"]["constraints"],
    ]
    for relative in paths:
        source = ROOT / relative
        target = tmp_path / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(source.read_text(encoding="utf-8"), encoding="utf-8")


def test_repository_eq45_nontriviality_scope_is_truthful() -> None:
    result = audit_eq45_nontriviality_scope_file(CONTRACT, repo_root=ROOT)

    assert result["contract_pass"] is True
    assert result["eq45_coefficient_norm"] > result["eq45_coefficient_guard_min_exclusive"]
    assert result["eq45_activity_speed_max"] > 0.0
    assert result["eq45_activity_speed_rms"] > 0.0
    assert result["eq45_cr001_energy_status"] == "pending_unknown"
    assert result["candidate_local_export_blocked_by_energy_status"] is False
    assert result["velocity_export_ready"] is True
    assert result["pde_validated"] is False
    assert result["paper_exact"] is False


@pytest.mark.parametrize(
    "key",
    [
        "eq45_coefficient_guard_is_cr001_energy_evidence",
        "eq45_activity_probe_is_cr001_energy_evidence",
    ],
)
def test_candidate_local_nonzero_evidence_cannot_become_cr001_energy_evidence(key: str) -> None:
    contract = _load(CONTRACT)
    mutated = copy.deepcopy(contract)
    mutated["state_semantics"][key] = True

    with pytest.raises(ValueError, match="cannot become CR001 energy evidence"):
        audit_eq45_nontriviality_scope(mutated, repo_root=ROOT)


def test_pending_cr001_energy_cannot_be_promoted_without_measurement() -> None:
    contract = _load(CONTRACT)
    mutated = copy.deepcopy(contract)
    mutated["state_semantics"]["eq45_cr001_energy_status"] = "pass"

    with pytest.raises(ValueError, match="must remain pending_unknown"):
        audit_eq45_nontriviality_scope(mutated, repo_root=ROOT)


def test_unassessed_cr001_energy_cannot_block_candidate_local_export() -> None:
    contract = _load(CONTRACT)
    mutated = copy.deepcopy(contract)
    mutated["eq45_candidate"]["candidate_local_export_requires_cr001_energy_normalization"] = True

    with pytest.raises(ValueError, match="must not become a candidate-local export gate"):
        audit_eq45_nontriviality_scope(mutated, repo_root=ROOT)


def test_cr001_energy_contract_drift_fails_closed() -> None:
    contract = _load(CONTRACT)
    mutated = copy.deepcopy(contract)
    mutated["cr001_nontriviality"]["reference_energy"] = 0.5

    with pytest.raises(ValueError, match="reference_energy"):
        audit_eq45_nontriviality_scope(mutated, repo_root=ROOT)


def test_eq45_scientific_claim_promotion_fails_closed(tmp_path: Path) -> None:
    contract = _load(CONTRACT)
    _copy_inputs(tmp_path, contract)
    candidate_path = tmp_path / contract["eq45_candidate"]["artifact"]
    candidate = _load(candidate_path)
    candidate["truth_boundary"]["pde_validated"] = True
    candidate_path.write_text(json.dumps(candidate), encoding="utf-8")

    with pytest.raises(ValueError, match="unsupported Eq45 claim promotion: pde_validated"):
        audit_eq45_nontriviality_scope(contract, repo_root=tmp_path)
