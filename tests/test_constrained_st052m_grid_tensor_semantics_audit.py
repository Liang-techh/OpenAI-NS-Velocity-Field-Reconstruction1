from __future__ import annotations

import json
from pathlib import Path
import shutil

import pytest

from openai_ns_reconstruction.constrained_st052m_grid_tensor_semantics_audit import (
    CONTRACT_FILENAME,
    audit_grid_tensor_semantics,
    axis_permutation_witness,
)


ROOT = Path(__file__).resolve().parents[1]


def _copy_audit_inputs(tmp_path: Path) -> Path:
    (tmp_path / "configs").mkdir(parents=True)
    (tmp_path / "src" / "openai_ns_reconstruction").mkdir(parents=True)
    shutil.copy2(ROOT / "configs" / CONTRACT_FILENAME, tmp_path / "configs" / CONTRACT_FILENAME)
    shutil.copy2(ROOT / "configs" / "constraints.json", tmp_path / "configs" / "constraints.json")
    shutil.copy2(
        ROOT / "src" / "openai_ns_reconstruction" / "st052_grid_export.py",
        tmp_path / "src" / "openai_ns_reconstruction" / "st052_grid_export.py",
    )
    return tmp_path


def test_live_grid_tensor_semantics_audit_passes() -> None:
    report = audit_grid_tensor_semantics(ROOT)
    assert report["live_export_contract_checked"] is True
    assert report["cr001_checked_unchanged"] is True
    witness = report["axis_permutation_witness"]
    assert witness["shape"] == [33, 33, 33]
    assert witness["coordinate_vectors_equal"] is True
    assert witness["absolute_difference"] > 0.0
    assert witness["shape_only_semantics_sufficient"] is False


def test_equal_cubic_shape_does_not_bind_axis_semantics() -> None:
    witness = axis_permutation_witness()
    assert witness["shape"] == [33, 33, 33]
    assert witness["original_value"] != witness["permuted_value"]


@pytest.mark.parametrize(
    "claim",
    [
        "shape_only_consumer_axis_semantics_verified",
        "matlab_octave_consumer_layout_metadata_enforced",
        "component_basis_machine_bound_in_consumer",
        "axis_permutation_equivalent_to_same_physical_field",
        "derived_grid_vorticity_semantics_verified",
        "visual_correspondence_verified",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    ],
)
def test_unverified_semantic_or_scientific_claim_promotion_fails_closed(
    tmp_path: Path, claim: str
) -> None:
    root = _copy_audit_inputs(tmp_path)
    contract_path = root / "configs" / CONTRACT_FILENAME
    contract = json.loads(contract_path.read_text(encoding="utf-8"))
    contract["verified_state"][claim] = True
    contract_path.write_text(json.dumps(contract, indent=2) + "\n", encoding="utf-8")
    with pytest.raises(AssertionError, match="promoted"):
        audit_grid_tensor_semantics(root)


def test_cr001_threshold_drift_fails_closed(tmp_path: Path) -> None:
    root = _copy_audit_inputs(tmp_path)
    constraints_path = root / "configs" / "constraints.json"
    constraints = json.loads(constraints_path.read_text(encoding="utf-8"))
    constraints["validation"]["thresholds"]["pde_residual_max"] = 0.01
    constraints_path.write_text(json.dumps(constraints, indent=2) + "\n", encoding="utf-8")
    with pytest.raises(AssertionError, match="pde_residual_max drifted"):
        audit_grid_tensor_semantics(root)


def test_free_force_prohibition_drift_fails_closed(tmp_path: Path) -> None:
    root = _copy_audit_inputs(tmp_path)
    constraints_path = root / "configs" / "constraints.json"
    constraints = json.loads(constraints_path.read_text(encoding="utf-8"))
    constraints["forcing"]["restriction"] = "Only a,c may be fitted."
    constraints_path.write_text(json.dumps(constraints, indent=2) + "\n", encoding="utf-8")
    with pytest.raises(AssertionError, match="free forcing prohibition drifted"):
        audit_grid_tensor_semantics(root)
