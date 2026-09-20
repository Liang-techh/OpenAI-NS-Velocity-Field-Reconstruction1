from __future__ import annotations

import json
from pathlib import Path

import pytest

from openai_ns_reconstruction.audit_kokuno_strict_inner_transport_mean_nonlinear_scope import (
    AuditError,
    CONTRACT_REL,
    CONSTRAINTS_REL,
    TARGET_REL,
    audit_repo,
)


REPO_ROOT = Path(__file__).resolve().parents[1]


def _copy_inputs(tmp_path: Path) -> Path:
    for rel in (CONTRACT_REL, CONSTRAINTS_REL, TARGET_REL):
        target = tmp_path / rel
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes((REPO_ROOT / rel).read_bytes())
    return tmp_path


def _read_contract(root: Path) -> dict:
    return json.loads((root / CONTRACT_REL).read_text())


def _write_contract(root: Path, data: dict) -> None:
    (root / CONTRACT_REL).write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")


def _read_constraints(root: Path) -> dict:
    return json.loads((root / CONSTRAINTS_REL).read_text())


def _write_constraints(root: Path, data: dict) -> None:
    (root / CONSTRAINTS_REL).write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")


def test_baseline_contract_passes_and_witnesses_noncommutation() -> None:
    report = audit_repo(REPO_ROOT)
    assert report["ok"] is True
    witness = report["manufactured_witness"]
    assert witness["sampled_mean_velocity_cylindrical"] == pytest.approx([0.0, 0.0, 0.0], abs=1e-12)
    assert witness["sampled_mean_advection_cylindrical"] == pytest.approx([-0.128, 0.0, 0.0], abs=1e-12)
    assert witness["advection_of_mean_velocity_cylindrical"] == [0.0, 0.0, 0.0]
    assert witness["noncommutation_radial_gap"] == pytest.approx(-0.128, abs=1e-15)
    assert report["scientific_status_promoted"] is False


def test_rejects_mean_advection_to_advection_of_mean_promotion(tmp_path: Path) -> None:
    root = _copy_inputs(tmp_path)
    contract = _read_contract(root)
    contract["nonlinear_mean_boundary"]["mean_of_advection_equals_advection_of_mean"] = True
    _write_contract(root, contract)
    with pytest.raises(AuditError, match="mean_of_advection_equals_advection_of_mean"):
        audit_repo(root)


def test_rejects_mean_transport_to_transport_of_mean_promotion(tmp_path: Path) -> None:
    root = _copy_inputs(tmp_path)
    contract = _read_contract(root)
    contract["nonlinear_mean_boundary"]["mean_of_transport_equals_transport_of_mean"] = True
    _write_contract(root, contract)
    with pytest.raises(AuditError, match="mean_of_transport_equals_transport_of_mean"):
        audit_repo(root)


def test_rejects_closed_mean_pde_promotion(tmp_path: Path) -> None:
    root = _copy_inputs(tmp_path)
    contract = _read_contract(root)
    contract["nonlinear_mean_boundary"]["closed_mean_velocity_pde_established"] = True
    _write_contract(root, contract)
    with pytest.raises(AuditError, match="closed_mean_velocity_pde_established"):
        audit_repo(root)


def test_rejects_pde_validation_promotion(tmp_path: Path) -> None:
    root = _copy_inputs(tmp_path)
    contract = _read_contract(root)
    contract["truth_boundary"]["pde_validated"] = True
    _write_contract(root, contract)
    with pytest.raises(AuditError, match="pde_validated"):
        audit_repo(root)


def test_rejects_st006_comparability_promotion(tmp_path: Path) -> None:
    root = _copy_inputs(tmp_path)
    contract = _read_contract(root)
    contract["truth_boundary"]["st006_same_protocol_comparable"] = True
    _write_contract(root, contract)
    with pytest.raises(AuditError, match="st006_same_protocol_comparable"):
        audit_repo(root)


def test_rejects_momentum_threshold_relaxation_in_canonical_constraints(tmp_path: Path) -> None:
    root = _copy_inputs(tmp_path)
    constraints = _read_constraints(root)
    constraints["validation"]["thresholds"]["pde_residual_max"] = 0.002
    _write_constraints(root, constraints)
    with pytest.raises(AuditError, match="pde_residual_max"):
        audit_repo(root)


def test_rejects_momentum_threshold_relaxation_in_snapshot(tmp_path: Path) -> None:
    root = _copy_inputs(tmp_path)
    contract = _read_contract(root)
    contract["cr001_snapshot"]["pde_residual_L2"] = 0.002
    _write_contract(root, contract)
    with pytest.raises(AuditError, match="pde_residual_L2"):
        audit_repo(root)


def test_rejects_free_residual_defined_force_permission(tmp_path: Path) -> None:
    root = _copy_inputs(tmp_path)
    contract = _read_contract(root)
    contract["cr001_snapshot"]["residual_defined_pointwise_free_force_allowed"] = True
    _write_contract(root, contract)
    with pytest.raises(AuditError, match="free residual force"):
        audit_repo(root)


def test_rejects_amplitude_collapse_permission(tmp_path: Path) -> None:
    root = _copy_inputs(tmp_path)
    contract = _read_contract(root)
    contract["cr001_snapshot"]["amplitude_collapse_allowed"] = True
    _write_contract(root, contract)
    with pytest.raises(AuditError, match="amplitude collapse"):
        audit_repo(root)


def test_rejects_post_hoc_threshold_relaxation_permission(tmp_path: Path) -> None:
    root = _copy_inputs(tmp_path)
    contract = _read_contract(root)
    contract["cr001_snapshot"]["post_hoc_threshold_relaxation_allowed"] = True
    _write_contract(root, contract)
    with pytest.raises(AuditError, match="threshold relaxation"):
        audit_repo(root)


def test_rejects_autonomous_angular_protocol_laundered_as_public_source(tmp_path: Path) -> None:
    root = _copy_inputs(tmp_path)
    contract = _read_contract(root)
    entry = "the angular ladder [32,64,128]"
    contract["classification"]["autonomous_design"].remove(entry)
    contract["classification"]["public_source_fact"].append(entry)
    _write_contract(root, contract)
    with pytest.raises(AuditError, match="angular ladder must remain autonomous_design"):
        audit_repo(root)


def test_rejects_pde_pending_as_velocity_delivery_blocker(tmp_path: Path) -> None:
    root = _copy_inputs(tmp_path)
    contract = _read_contract(root)
    contract["truth_boundary"]["pde_pending_blocks_callable_velocity_delivery"] = True
    _write_contract(root, contract)
    with pytest.raises(AuditError, match="pde_pending_blocks_callable_velocity_delivery"):
        audit_repo(root)


def test_rejects_target_source_identity_drift(tmp_path: Path) -> None:
    root = _copy_inputs(tmp_path)
    source_path = root / TARGET_REL
    source = source_path.read_text().replace(
        'TASK = "KOKUNO-A3-STRICT-INNER-TRANSPORT-MEAN-072"',
        'TASK = "KOKUNO-A3-STRICT-INNER-TRANSPORT-MEAN-DRIFTED"',
    )
    source_path.write_text(source)
    with pytest.raises(AuditError, match="target source anchor missing"):
        audit_repo(root)


def test_rejects_manufactured_witness_promoted_to_source_fact(tmp_path: Path) -> None:
    root = _copy_inputs(tmp_path)
    contract = _read_contract(root)
    contract["manufactured_annular_witness"]["candidate_or_source_fact_claimed"] = True
    _write_contract(root, contract)
    with pytest.raises(AuditError, match="candidate/source fact"):
        audit_repo(root)
