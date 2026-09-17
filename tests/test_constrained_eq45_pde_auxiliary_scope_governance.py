from __future__ import annotations

import copy
import json
from pathlib import Path

import pytest

from openai_ns_reconstruction.constrained_eq45_pde_auxiliary_scope_governance import (
    audit_documents,
    audit_eq45_pde_auxiliary_scope,
)


ROOT = Path(__file__).resolve().parents[1]


def _read(relative: str):
    return json.loads((ROOT / relative).read_text(encoding="utf-8"))


def _documents():
    return (
        _read("configs/eq45_pde_auxiliary_scope_contract.json"),
        _read("configs/constraints.json"),
        _read("artifacts/constrained/eq45_velocity_candidate_seed.json"),
        _read("configs/delivery_state_contract.json"),
    )


def test_current_scope_keeps_export_and_pde_auxiliaries_independent():
    report = audit_eq45_pde_auxiliary_scope(ROOT)
    assert report["velocity_export_ready"] is True
    assert report["pressure_binding"] == "pending_unknown"
    assert report["forcing_binding"] == "pending_unknown"
    assert report["registered_forcing_mode"] == "restricted_two_parameter_family"
    assert report["residual_defined_forcing_allowed"] is False
    assert report["pde_validated"] is False
    assert report["velocity_changed"] is False


def test_residual_defined_or_unrestricted_forcing_drift_fails_closed():
    contract, constraints, candidate, delivery = _documents()
    mutated = copy.deepcopy(constraints)
    mutated["forcing"]["mode"] = "residual_defined"
    with pytest.raises(ValueError, match="forcing family drifted"):
        audit_documents(contract, mutated, candidate, delivery)

    mutated = copy.deepcopy(constraints)
    mutated["forcing"]["restriction"] = "Any pointwise force may be fitted."
    with pytest.raises(ValueError, match="residual-dependent"):
        audit_documents(contract, mutated, candidate, delivery)


def test_pending_auxiliary_cannot_become_export_blocker_or_pde_success():
    contract, constraints, candidate, delivery = _documents()
    mutated_contract = copy.deepcopy(contract)
    mutated_contract["policy"]["pending_auxiliary_blocks_velocity_export"] = True
    with pytest.raises(ValueError, match="must not block"):
        audit_documents(mutated_contract, constraints, candidate, delivery)

    mutated_candidate = copy.deepcopy(candidate)
    mutated_candidate["truth_boundary"]["pde_validated"] = True
    with pytest.raises(ValueError, match="cannot be true"):
        audit_documents(contract, constraints, mutated_candidate, delivery)

    mutated_contract = copy.deepcopy(contract)
    mutated_contract["policy"]["pressure_or_forcing_solver_success_implies_pde_validation"] = True
    with pytest.raises(ValueError, match="solver success"):
        audit_documents(mutated_contract, constraints, candidate, delivery)

    mutated_contract = copy.deepcopy(contract)
    mutated_contract["classification"]["eq45_pressure_binding"] = "public_source_fact"
    with pytest.raises(ValueError, match="classification drifted"):
        audit_documents(mutated_contract, constraints, candidate, delivery)


def test_frozen_velocity_artifact_must_not_smuggle_pde_auxiliaries():
    contract, constraints, candidate, delivery = _documents()
    mutated = copy.deepcopy(candidate)
    mutated["forcing"] = {"mode": "residual_defined"}
    with pytest.raises(ValueError, match="unexpectedly embeds"):
        audit_documents(contract, constraints, mutated, delivery)
