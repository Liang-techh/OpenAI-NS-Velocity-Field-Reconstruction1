from __future__ import annotations

import copy
import json
from pathlib import Path

import pytest

from openai_ns_reconstruction.constrained_st052m_grid_representation_audit import (
    DEFAULT_CONTRACT_PATH,
    GridRepresentationAuditError,
    audit_grid_representation_contract,
    load_contract,
    sampled_grid_nonuniqueness_witness,
)


def _write_mutation(tmp_path: Path, mutate) -> Path:
    contract = copy.deepcopy(load_contract())
    mutate(contract)
    out = tmp_path / "mutated_grid_representation_contract.json"
    out.write_text(json.dumps(contract, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return out


def test_integrated_grid_representation_contract_audits_cleanly() -> None:
    report = audit_grid_representation_contract(DEFAULT_CONTRACT_PATH)
    assert report["grid_shape_xyz"] == [33, 33, 33]
    assert report["times"] == [0.25, 0.375, 0.5, 0.625, 0.75]
    assert report["callable_grid_node_replay_atol"] == 5e-12
    assert report["sampled_grid_is_complete_continuous_velocity_representation"] is False
    assert report["underlying_runtime_callable_delivery_blocked_by_this_audit"] is False


def test_finite_grid_does_not_identify_off_grid_continuous_field() -> None:
    witness = sampled_grid_nonuniqueness_witness()
    assert witness["grid_spacing"] == pytest.approx(0.125)
    assert witness["max_abs_at_serialized_nodes"] < 1e-12
    assert witness["abs_at_first_half_grid_point"] > 1.0 - 1e-12


def test_rejects_laundering_grid_into_complete_continuous_velocity(tmp_path: Path) -> None:
    path = _write_mutation(
        tmp_path,
        lambda contract: contract["truth_boundary"].__setitem__(
            "sampled_grid_is_complete_continuous_velocity_representation", True
        ),
    )
    with pytest.raises(GridRepresentationAuditError, match="sampled-grid claim must remain false"):
        audit_grid_representation_contract(path)


def test_rejects_unversioned_interpolation_claim(tmp_path: Path) -> None:
    path = _write_mutation(
        tmp_path,
        lambda contract: contract["truth_boundary"].__setitem__(
            "off_grid_interpolant_specified_and_identity_bound", True
        ),
    )
    with pytest.raises(GridRepresentationAuditError, match="sampled-grid claim must remain false"):
        audit_grid_representation_contract(path)


def test_rejects_pde_transfer_from_grid_node_replay(tmp_path: Path) -> None:
    path = _write_mutation(
        tmp_path,
        lambda contract: contract["truth_boundary"].__setitem__(
            "pde_validation_transfers_to_grid_interpolant", True
        ),
    )
    with pytest.raises(GridRepresentationAuditError, match="sampled-grid claim must remain false"):
        audit_grid_representation_contract(path)


def test_rejects_relabelling_engineering_replay_tolerance_as_scientific(tmp_path: Path) -> None:
    path = _write_mutation(
        tmp_path,
        lambda contract: contract["representation"].__setitem__(
            "callable_grid_node_replay_tolerance_role", "pde_acceptance_threshold"
        ),
    )
    with pytest.raises(GridRepresentationAuditError, match="must not be relabelled"):
        audit_grid_representation_contract(path)
