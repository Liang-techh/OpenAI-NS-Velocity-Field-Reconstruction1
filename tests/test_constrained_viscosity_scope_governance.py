from __future__ import annotations

import copy
import json
from pathlib import Path

import pytest

from openai_ns_reconstruction.constrained_viscosity_scope_governance import (
    audit_repository,
    audit_viscosity_scope,
)


def _contract() -> dict:
    return json.loads(
        (
            Path(__file__).parents[1]
            / "configs"
            / "function_first_viscosity_scope.json"
        ).read_text(encoding="utf-8")
    )


def _problem() -> dict:
    return {"nu": 0.01}


def _profile() -> dict:
    return {"metadata": {"viscosity": 1.0}}


def test_current_repository_viscosity_scope_is_explicit() -> None:
    result = audit_repository(Path(__file__).parents[1])
    assert result["problem_nu"] == pytest.approx(0.01)
    assert result["profile_equation_nu"] == pytest.approx(1.0)
    assert result["same_viscosity"] is False
    assert result["direct_pde_evidence_transfer_allowed"] is False
    assert result["mapping_status"] == "pending_unknown"
    assert result["function_first_profile_callable"] is True
    assert result["function_first_profile_pde_contract_compatible"] is False
    assert result["pde_validated"] is False


def test_viscosity_scope_mutations_fail_closed() -> None:
    contract = _contract()

    changed_problem = _problem()
    changed_problem["nu"] = 1.0
    with pytest.raises(ValueError, match="problem viscosity drifted"):
        audit_viscosity_scope(changed_problem, _profile(), contract)

    changed_profile = _profile()
    changed_profile["metadata"]["viscosity"] = 0.01
    with pytest.raises(ValueError, match="profile viscosity metadata drifted"):
        audit_viscosity_scope(_problem(), changed_profile, contract)

    free_transfer = copy.deepcopy(contract)
    free_transfer["compatibility"]["direct_pde_evidence_transfer_allowed"] = True
    with pytest.raises(ValueError, match="direct PDE evidence transfer is forbidden"):
        audit_viscosity_scope(_problem(), _profile(), free_transfer)


def test_claim_and_source_promotions_fail_closed() -> None:
    contract = _contract()

    sourced_nu = copy.deepcopy(contract)
    sourced_nu["active_problem"]["classification"] = "public_source_fact"
    with pytest.raises(ValueError, match="CR001 nu=0.01 choice"):
        audit_viscosity_scope(_problem(), _profile(), sourced_nu)

    mapped = copy.deepcopy(contract)
    mapped["compatibility"]["mapping_status"] = "verified"
    with pytest.raises(ValueError, match="pending_unknown"):
        audit_viscosity_scope(_problem(), _profile(), mapped)

    promoted = copy.deepcopy(contract)
    promoted["claim_status"]["pde_validated"] = True
    with pytest.raises(ValueError, match="pde_validated"):
        audit_viscosity_scope(_problem(), _profile(), promoted)
