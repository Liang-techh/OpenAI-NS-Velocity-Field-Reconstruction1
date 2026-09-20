from __future__ import annotations

from copy import deepcopy
from pathlib import Path

import pytest

from openai_ns_reconstruction import constrained_st052_runtime_dependency_acceptance as acceptance


ROOT = Path(__file__).resolve().parents[1]
EXPECTED_RUNTIME_SHA256 = "a0a4682d62677911758e6f9240439c133ee1e3787d3801b87c25de6dcd32f89c"


def _contract() -> dict:
    return acceptance._load_contract(ROOT)


def _audit_mutated(monkeypatch: pytest.MonkeyPatch, contract: dict) -> None:
    monkeypatch.setattr(acceptance, "_load_contract", lambda _repo_root: contract)
    acceptance.audit(ROOT)


def test_live_runtime_dependency_acceptance_receipt() -> None:
    receipt = acceptance.audit(ROOT)
    assert receipt["task_id"] == "A8-DELIVERY-01"
    assert receipt["candidate_id"] == "ST052-M-linear-temporal-child-v1"
    assert receipt["exact_source_runtime_identity_sha256"] == EXPECTED_RUNTIME_SHA256
    assert receipt["external_runtime_dependency_accepted_for_visualization_delivery"] is True
    assert receipt["standalone_package_parent_runtime_ready"] is False
    assert receipt["velocity_export_ready"] is False
    assert receipt["visualization_ready"] is False
    assert receipt["pde_validated"] is False
    assert receipt["next_task"] == "A8-DELIVERY-02"


def test_source_runtime_identity_drift_fails_closed(monkeypatch: pytest.MonkeyPatch) -> None:
    contract = deepcopy(_contract())
    contract["required_source_runtime"]["source_head"] = "0" * 40
    with pytest.raises(ValueError, match="source_head"):
        _audit_mutated(monkeypatch, contract)


def test_standalone_runtime_promotion_fails_closed(monkeypatch: pytest.MonkeyPatch) -> None:
    contract = deepcopy(_contract())
    contract["runtime_policy"]["standalone_package_parent_runtime_ready"] = True
    with pytest.raises(ValueError, match="standalone package"):
        _audit_mutated(monkeypatch, contract)


def test_portable_runtime_identity_laundering_fails_closed(monkeypatch: pytest.MonkeyPatch) -> None:
    contract = deepcopy(_contract())
    contract["runtime_policy"]["portable_bitwise_identity_across_dependency_builds_claimed"] = True
    with pytest.raises(ValueError, match="portable bitwise"):
        _audit_mutated(monkeypatch, contract)


def test_runtime_receipt_requirement_cannot_be_removed(monkeypatch: pytest.MonkeyPatch) -> None:
    contract = deepcopy(_contract())
    contract["runtime_policy"]["integration_smoke_must_record_resolved_python_numpy_scipy_sympy_versions"] = False
    with pytest.raises(ValueError, match="record the resolved numerical runtime"):
        _audit_mutated(monkeypatch, contract)


def test_scientific_state_promotion_fails_closed(monkeypatch: pytest.MonkeyPatch) -> None:
    base_contract = _contract()
    for key in ("velocity_export_ready", "visualization_ready", "pde_validated", "openai_field_identified"):
        contract = deepcopy(base_contract)
        contract["truth_boundary"][key] = True
        with monkeypatch.context() as scoped:
            with pytest.raises(ValueError, match=key):
                _audit_mutated(scoped, contract)


def test_runtime_acceptance_unblocks_only_next_delivery_smoke(monkeypatch: pytest.MonkeyPatch) -> None:
    contract = deepcopy(_contract())
    contract["unblocks"] = ["A8-DELIVERY-02", "PDE-01"]
    with pytest.raises(ValueError, match="unblock only"):
        _audit_mutated(monkeypatch, contract)
