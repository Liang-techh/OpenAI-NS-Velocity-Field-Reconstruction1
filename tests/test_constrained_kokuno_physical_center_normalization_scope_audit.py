from __future__ import annotations

import copy
import json

import pytest

from openai_ns_reconstruction.constrained_kokuno_physical_center_normalization_scope_audit import (
    CONTRACT_RELATIVE_PATH,
    CONSTRAINTS_RELATIVE_PATH,
    _default_a5_payload,
    _repo_root,
    audit_normalization_scope,
    validate_normalization_scope_contract,
)


def _load_contract() -> dict:
    return json.loads((_repo_root() / CONTRACT_RELATIVE_PATH).read_text())


def _load_constraints() -> dict:
    return json.loads((_repo_root() / CONSTRAINTS_RELATIVE_PATH).read_text())


def test_live_normalization_scope_audit_passes() -> None:
    receipt = audit_normalization_scope()
    assert receipt["audit_passed"] is True
    assert receipt["source_C_role_separated_from_CR001_energy_normalization"] is True
    assert receipt["source_complex_C_normalization_certified"] is False
    assert receipt["physical_center_profile_ingest_admitted"] is False
    assert receipt["kokuno_velocity_export_ready"] is False
    assert receipt["pde_validated"] is False


def test_rejects_source_C_as_posthoc_pde_amplitude_knob() -> None:
    contract = _load_contract()
    constraints = _load_constraints()
    mutated = copy.deepcopy(contract)
    mutated["source_C_role"]["may_be_tuned_posthoc_to_reduce_NS_residual"] = True
    with pytest.raises(ValueError, match="may_be_tuned_posthoc_to_reduce_NS_residual"):
        validate_normalization_scope_contract(mutated, constraints)


def test_rejects_energy_normalization_as_source_C_certificate() -> None:
    contract = _load_contract()
    constraints = _load_constraints()
    mutated = copy.deepcopy(contract)
    mutated["cr001_energy_normalization_role"]["certifies_source_complex_C_condition"] = True
    with pytest.raises(ValueError, match="certifies_source_complex_C_condition"):
        validate_normalization_scope_contract(mutated, constraints)


def test_rejects_cr001_momentum_threshold_relaxation() -> None:
    contract = _load_contract()
    constraints = _load_constraints()
    mutated = copy.deepcopy(constraints)
    mutated["validation"]["thresholds"]["pde_residual_max"] = 0.002
    with pytest.raises(ValueError, match="CR001 momentum_max"):
        validate_normalization_scope_contract(contract, mutated)


def test_rejects_a5_velocity_export_promotion() -> None:
    contract = _load_contract()
    constraints = _load_constraints()
    a5_payload = _default_a5_payload()
    mutated = copy.deepcopy(a5_payload)
    mutated["truth_boundary"]["velocity_export_ready"] = True
    with pytest.raises(ValueError, match="velocity_export_ready"):
        validate_normalization_scope_contract(contract, constraints, a5_payload=mutated)


def test_rejects_configured_C_as_public_source_data() -> None:
    contract = _load_contract()
    constraints = _load_constraints()
    mutated = copy.deepcopy(contract)
    mutated["classification"]["autonomous_design"] = [
        item
        for item in mutated["classification"]["autonomous_design"]
        if "configured numerical C" not in item
    ]
    with pytest.raises(ValueError, match="configured C must remain autonomous"):
        validate_normalization_scope_contract(mutated, constraints)
