from __future__ import annotations

import copy
import json
from pathlib import Path

import pytest

from openai_ns_reconstruction.constrained_st052m_temporal_serialization_governance import (
    EXPECTED_TEMPORAL_SPEC_SHA256,
    audit_contract,
    run_audit,
)
from openai_ns_reconstruction.st052_linear_temporal_transform import (
    DEFAULT_TEMPORAL_SPEC,
)


ROOT = Path(__file__).resolve().parents[1]
CONTRACT_PATH = ROOT / "configs" / "st052m_temporal_serialization_governance_contract.json"


def _contract() -> dict:
    return json.loads(CONTRACT_PATH.read_text(encoding="utf-8"))


def test_live_temporal_serialization_governance_passes() -> None:
    receipt = run_audit(root=ROOT)
    assert receipt == {
        "status": "pass",
        "temporal_transform_spec_sha256": EXPECTED_TEMPORAL_SPEC_SHA256,
        "transform_spec_save_load_ready": True,
        "complete_candidate_save_load_ready": False,
        "experimental_velocity_export_ready": False,
        "canonical_velocity_export_ready": True,
        "pde_validated": False,
    }
    assert DEFAULT_TEMPORAL_SPEC.sha256() == EXPECTED_TEMPORAL_SPEC_SHA256


def test_rejects_transform_sha_as_complete_candidate_sha() -> None:
    contract = _contract()
    contract["identity_semantics"]["transform_spec_sha_is_complete_candidate_sha"] = True
    with pytest.raises(ValueError, match="identity laundering"):
        audit_contract(contract, root=ROOT)


def test_rejects_spec_roundtrip_as_complete_candidate_save_load() -> None:
    contract = _contract()
    contract["identity_semantics"]["spec_save_load_is_complete_candidate_save_load"] = True
    with pytest.raises(ValueError, match="identity laundering"):
        audit_contract(contract, root=ROOT)


def test_rejects_parent_metadata_as_authentication() -> None:
    contract = _contract()
    contract["identity_semantics"]["parent_metadata_strings_authenticate_parent_callable"] = True
    with pytest.raises(ValueError, match="identity laundering"):
        audit_contract(contract, root=ROOT)


def test_rejects_premature_temporal_export_promotion() -> None:
    contract = _contract()
    contract["live_temporal_transform"]["experimental_velocity_export_ready"] = True
    with pytest.raises(ValueError, match="premature temporal delivery promotion"):
        audit_contract(contract, root=ROOT)


def test_rejects_missing_complete_child_prerequisite() -> None:
    contract = _contract()
    contract["requirements_before_temporal_child_export_ready"].remove(
        "assign immutable parent artifact or candidate SHA and bind it at runtime"
    )
    with pytest.raises(ValueError, match="prerequisite removed"):
        audit_contract(contract, root=ROOT)


def test_rejects_canonical_delivery_relabelling() -> None:
    contract = _contract()
    contract["canonical_delivery"]["candidate_family"] = "st052m_temporal_v1"
    with pytest.raises(ValueError, match="canonical family relabelled"):
        audit_contract(contract, root=ROOT)


def test_rejects_pde_truth_promotion() -> None:
    contract = _contract()
    contract["canonical_delivery"]["pde_validated"] = True
    with pytest.raises(ValueError, match="canonical truth state promoted"):
        audit_contract(contract, root=ROOT)


def test_rejects_cr001_threshold_laundering() -> None:
    contract = _contract()
    contract["cr001_nonmutation"]["pde_residual_max"] = 0.01
    with pytest.raises(ValueError, match="CR001 threshold drift"):
        audit_contract(contract, root=ROOT)


def test_rejects_cr001_validation_sample_drift() -> None:
    contract = _contract()
    contract["cr001_nonmutation"]["held_out_points"] = 2048
    with pytest.raises(ValueError, match="CR001 sample drift"):
        audit_contract(contract, root=ROOT)


def test_rejects_source_classification_drift() -> None:
    contract = _contract()
    mutated = copy.deepcopy(contract)
    mutated["source_classification"][2]["classification"] = "public_source_fact"
    with pytest.raises(ValueError, match="source classification drift"):
        audit_contract(mutated, root=ROOT)
