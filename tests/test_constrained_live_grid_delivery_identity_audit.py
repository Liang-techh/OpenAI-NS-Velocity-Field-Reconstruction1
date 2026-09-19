from __future__ import annotations

import json
from pathlib import Path
import shutil
from typing import Callable

import pytest

from openai_ns_reconstruction.constrained_live_grid_delivery_identity_audit import (
    CONTRACT_FILENAME,
    audit_live_grid_delivery_identity,
)


ROOT = Path(__file__).resolve().parents[1]


def _copy_audit_inputs(tmp_path: Path) -> Path:
    (tmp_path / "configs").mkdir(parents=True)
    (tmp_path / "src" / "openai_ns_reconstruction").mkdir(parents=True)
    (tmp_path / "scripts").mkdir(parents=True)
    shutil.copy2(ROOT / "configs" / CONTRACT_FILENAME, tmp_path / "configs" / CONTRACT_FILENAME)
    shutil.copy2(ROOT / "configs" / "constraints.json", tmp_path / "configs" / "constraints.json")
    shutil.copy2(
        ROOT / "src" / "openai_ns_reconstruction" / "st052_grid_export.py",
        tmp_path / "src" / "openai_ns_reconstruction" / "st052_grid_export.py",
    )
    shutil.copy2(
        ROOT / "scripts" / "smoke_st052_grid_octave.m",
        tmp_path / "scripts" / "smoke_st052_grid_octave.m",
    )
    return tmp_path


def _mutate_contract(tmp_path: Path, mutator: Callable[[dict], None]) -> Path:
    root = _copy_audit_inputs(tmp_path)
    path = root / "configs" / CONTRACT_FILENAME
    obj = json.loads(path.read_text(encoding="utf-8"))
    mutator(obj)
    path.write_text(json.dumps(obj, indent=2) + "\n", encoding="utf-8")
    return root


def test_live_grid_delivery_identity_audit_passes() -> None:
    report = audit_live_grid_delivery_identity(ROOT)
    assert report["live_st052_grid_delivery_checked"] is True
    assert report["octave_software_scope_checked"] is True
    assert report["continuous_representation_boundary_checked"] is True
    assert report["cr001_checked_unchanged"] is True
    assert report["truth_boundary"]["visual_correspondence_verified"] is False
    assert report["truth_boundary"]["pde_validated"] is False


def test_rejects_stale_or_rewritten_grid_times(tmp_path: Path) -> None:
    root = _mutate_contract(
        tmp_path,
        lambda obj: obj["current_st052_grid_delivery"].__setitem__(
            "times", [0.25, 0.5, 0.75]
        ),
    )
    with pytest.raises(AssertionError, match="grid times drifted"):
        audit_live_grid_delivery_identity(root)


@pytest.mark.parametrize(
    "key",
    [
        "sampled_grid_is_continuum_velocity_api",
        "node_equality_implies_continuous_function_identity",
        "off_grid_interpolant_specified",
        "off_grid_error_bound_established",
        "continuous_derivative_equivalence_established",
        "grid_derivatives_are_cr001_pde_acceptance_evidence",
        "grid_nodes_may_replace_cr001_held_out_sample",
        "five_grid_times_may_replace_six_cr001_validation_times",
        "visual_interpolation_may_be_called_exact_candidate_without_contract",
        "sampled_grid_may_be_called_exact_openai_field",
    ],
)
def test_rejects_finite_grid_claim_laundering(tmp_path: Path, key: str) -> None:
    root = _mutate_contract(
        tmp_path,
        lambda obj: obj["representation_scope"].__setitem__(key, True),
    )
    with pytest.raises(AssertionError, match="unverified finite-grid inference promoted"):
        audit_live_grid_delivery_identity(root)


@pytest.mark.parametrize(
    "key",
    [
        "actual_matlab_runtime_executed",
        "layout_metadata_field_enforced",
        "whole_candidate_identity_metadata_enforced",
        "grid_payload_checksum_metadata_enforced",
        "component_basis_machine_bound_from_metadata",
        "derived_grid_vorticity_scientific_semantics_verified",
    ],
)
def test_rejects_octave_smoke_as_stronger_semantic_evidence(
    tmp_path: Path, key: str
) -> None:
    root = _mutate_contract(
        tmp_path,
        lambda obj: obj["octave_consumer"].__setitem__(key, True),
    )
    with pytest.raises(AssertionError, match="unverified Octave/metadata claim promoted"):
        audit_live_grid_delivery_identity(root)


@pytest.mark.parametrize(
    "key",
    [
        "st052_velocity_export_ready",
        "visualization_ready",
        "visual_correspondence_verified",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    ],
)
def test_rejects_truth_state_promotion_without_evidence(tmp_path: Path, key: str) -> None:
    root = _mutate_contract(
        tmp_path,
        lambda obj: obj["truth_boundary"].__setitem__(key, True),
    )
    with pytest.raises(AssertionError, match="truth state promoted without evidence"):
        audit_live_grid_delivery_identity(root)


def test_rejects_engineering_replay_tolerance_relabelled_as_scientific_gate(
    tmp_path: Path,
) -> None:
    root = _mutate_contract(
        tmp_path,
        lambda obj: obj["current_st052_grid_delivery"].__setitem__(
            "callable_node_replay_tolerance_scope", "pde_acceptance_threshold"
        ),
    )
    with pytest.raises(AssertionError, match="callable node replay tolerance scope drifted"):
        audit_live_grid_delivery_identity(root)


def test_rejects_cr001_threshold_drift(tmp_path: Path) -> None:
    root = _copy_audit_inputs(tmp_path)
    path = root / "configs" / "constraints.json"
    obj = json.loads(path.read_text(encoding="utf-8"))
    obj["validation"]["thresholds"]["pde_residual_max"] = 0.01
    path.write_text(json.dumps(obj, indent=2) + "\n", encoding="utf-8")
    with pytest.raises(AssertionError, match="pde_residual_max drifted"):
        audit_live_grid_delivery_identity(root)


def test_rejects_octave_script_metadata_behavior_change_without_reaudit(tmp_path: Path) -> None:
    root = _copy_audit_inputs(tmp_path)
    path = root / "scripts" / "smoke_st052_grid_octave.m"
    path.write_text(
        path.read_text(encoding="utf-8")
        + "\n% future behavior: validate array_layout metadata\n",
        encoding="utf-8",
    )
    with pytest.raises(AssertionError, match="metadata behavior changed"):
        audit_live_grid_delivery_identity(root)
