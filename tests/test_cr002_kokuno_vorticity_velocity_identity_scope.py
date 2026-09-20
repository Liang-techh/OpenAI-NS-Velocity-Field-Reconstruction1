from __future__ import annotations

import json
import shutil
from pathlib import Path

import pytest

from openai_ns_reconstruction.audit_kokuno_vorticity_velocity_identity_scope import (
    CONFIG_REL,
    CONSTRAINTS_REL,
    UPSTREAM_REL,
    audit,
    mechanics_witness,
)


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def _copy_fragment(tmp_path: Path) -> Path:
    root = _repo_root()
    for rel in (CONFIG_REL, CONSTRAINTS_REL, UPSTREAM_REL):
        target = tmp_path / rel
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(root / rel, target)
    return tmp_path


def _rewrite_json(path: Path, mutator) -> None:
    payload = json.loads(path.read_text(encoding="utf-8"))
    mutator(payload)
    path.write_text(json.dumps(payload, sort_keys=True, indent=2) + "\n", encoding="utf-8")


def test_scope_audit_passes_on_frozen_upstream_and_cr001() -> None:
    receipt = audit(repo_root=_repo_root())
    assert receipt["status"] == "pass"
    assert receipt["scientific_admission_changed"] is False
    assert receipt["candidate_bytes_changed"] is False
    assert receipt["thresholds_changed"] is False


def test_same_vorticity_and_divergence_do_not_identify_local_velocity() -> None:
    witness = mechanics_witness()
    assert witness["max_divergence_difference"] <= 1.0e-8
    assert witness["max_curl_difference"] <= 1.0e-8
    assert witness["max_expected_curl_error"] <= 1.0e-8
    assert witness["min_velocity_difference"] >= 0.299999999


def test_rejects_vorticity_sha_promoted_to_velocity_identity(tmp_path: Path) -> None:
    repo = _copy_fragment(tmp_path)
    _rewrite_json(
        repo / CONFIG_REL,
        lambda payload: payload["identity_boundary"].__setitem__(
            "vorticity_sha256_is_velocity_candidate_identity", True
        ),
    )
    with pytest.raises(AssertionError, match="forbidden identity promotion"):
        audit(repo_root=repo)


def test_rejects_local_curl_divergence_promoted_to_unique_velocity(tmp_path: Path) -> None:
    repo = _copy_fragment(tmp_path)
    _rewrite_json(
        repo / CONFIG_REL,
        lambda payload: payload["identity_boundary"].__setitem__(
            "matching_vorticity_and_divergence_locally_identifies_velocity_without_global_boundary_data",
            True,
        ),
    )
    with pytest.raises(AssertionError, match="forbidden identity promotion"):
        audit(repo_root=repo)


def test_rejects_pde_or_visual_promotion_from_scoped_curl_audit(tmp_path: Path) -> None:
    repo = _copy_fragment(tmp_path)
    _rewrite_json(
        repo / CONFIG_REL,
        lambda payload: payload["route_state"].update(
            {
                "kokuno_visual_correspondence_verified": True,
                "kokuno_pde_validated": True,
            }
        ),
    )
    with pytest.raises(AssertionError, match="route truth state promoted"):
        audit(repo_root=repo)


def test_rejects_cr001_momentum_threshold_relaxation(tmp_path: Path) -> None:
    repo = _copy_fragment(tmp_path)
    constraints_path = repo / CONSTRAINTS_REL
    _rewrite_json(
        constraints_path,
        lambda payload: payload["validation"]["thresholds"].__setitem__(
            "pde_residual_max", 2.0e-3
        ),
    )
    with pytest.raises(AssertionError, match="canonical CR001 blob drift|CR001 threshold drift"):
        audit(repo_root=repo)


def test_rejects_provenance_class_laundering(tmp_path: Path) -> None:
    repo = _copy_fragment(tmp_path)
    _rewrite_json(
        repo / CONFIG_REL,
        lambda payload: payload["provenance_classes"].pop("pending_unknown"),
    )
    with pytest.raises(AssertionError, match="four-way provenance classification drift"):
        audit(repo_root=repo)


def test_scoped_kokuno_audit_does_not_block_canonical_eq45_delivery(tmp_path: Path) -> None:
    repo = _copy_fragment(tmp_path)
    _rewrite_json(
        repo / CONFIG_REL,
        lambda payload: payload["route_state"].__setitem__(
            "canonical_eq45_velocity_export_ready_unaffected", False
        ),
    )
    with pytest.raises(AssertionError, match="incorrectly blocks canonical Eq45 delivery"):
        audit(repo_root=repo)
