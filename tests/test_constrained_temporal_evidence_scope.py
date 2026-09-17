from __future__ import annotations

import copy
import json
from pathlib import Path

import pytest

from openai_ns_reconstruction.constrained_temporal_evidence_scope import (
    audit_temporal_evidence_scope,
    audit_temporal_evidence_scope_file,
)


def _contract() -> dict:
    return json.loads(Path("configs/temporal_evidence_scope_contract.json").read_text(encoding="utf-8"))


def test_repository_temporal_evidence_scope_is_truthful() -> None:
    result = audit_temporal_evidence_scope_file(
        "configs/temporal_evidence_scope_contract.json", repo_root=Path(".")
    )
    assert result["contract_pass"] is True
    assert result["pde_time_interval"] == [0.25, 0.75]
    assert result["delivery_time_interval"] == [0.25, 0.75]
    assert result["visualization_reference_times"] == [0.25, 0.5, 0.75]
    assert result["visual_to_pde_transfer_allowed"] is False
    assert result["pde_failure_blocks_velocity_delivery"] is False
    assert result["pde_validated"] is False


def test_pde_time_scope_cannot_drift_from_preregistration() -> None:
    mutated = copy.deepcopy(_contract())
    mutated["scopes"]["pde_validation"]["time_interval"] = [0.2, 0.8]
    with pytest.raises(ValueError, match="PDE time interval drifted"):
        audit_temporal_evidence_scope(mutated, repo_root=Path("."))

    mutated = copy.deepcopy(_contract())
    mutated["scopes"]["pde_validation"]["validation_times"][-1] = 0.74
    with pytest.raises(ValueError, match="PDE validation times drifted"):
        audit_temporal_evidence_scope(mutated, repo_root=Path("."))


def test_visualization_times_stay_autonomous_and_inside_callable_domain() -> None:
    mutated = copy.deepcopy(_contract())
    mutated["scopes"]["visualization_reference"]["classification"] = "pending"
    with pytest.raises(ValueError, match="visualization reference times are autonomous"):
        audit_temporal_evidence_scope(mutated, repo_root=Path("."))

    mutated = copy.deepcopy(_contract())
    mutated["scopes"]["visualization_reference"]["times"] = [0.25, 0.5, 0.8]
    with pytest.raises(ValueError, match="outside callable velocity interval"):
        audit_temporal_evidence_scope(mutated, repo_root=Path("."))


def test_evidence_transfer_and_claim_promotion_fail_closed() -> None:
    mutated = copy.deepcopy(_contract())
    mutated["evidence_transfer"]["visualization_snapshot_to_pde_validation"] = True
    with pytest.raises(ValueError, match="cannot be transferred"):
        audit_temporal_evidence_scope(mutated, repo_root=Path("."))

    mutated = copy.deepcopy(_contract())
    mutated["scopes"]["public_visual_time_mapping"]["status"] = "verified"
    with pytest.raises(ValueError, match="not independently verified"):
        audit_temporal_evidence_scope(mutated, repo_root=Path("."))

    mutated = copy.deepcopy(_contract())
    mutated["claim_status"]["pde_validated"] = True
    with pytest.raises(ValueError, match="unsupported claim promotion: pde_validated"):
        audit_temporal_evidence_scope(mutated, repo_root=Path("."))
