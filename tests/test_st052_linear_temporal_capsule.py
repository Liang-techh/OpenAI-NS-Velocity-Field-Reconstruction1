from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

from openai_ns_reconstruction.st052_linear_temporal_capsule import (
    PARENT_CANDIDATE_FILENAME,
    TASK_ID,
    WHOLE_MANIFEST_FILENAME,
    build_bundle,
    verify_bundle,
)
from openai_ns_reconstruction.st052_parent_capsule import SOURCE_HEAD


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _inputs(tmp_path: Path):
    candidate = tmp_path / "candidate.json"
    validation = tmp_path / "validation.json"
    parent_manifest = tmp_path / "parent_manifest.json"
    candidate.write_text('{"schema":"test-parent","pde_validated":false}\n')
    validation.write_text(
        json.dumps(
            {
                "pde_validated": False,
                "all_numeric_gates_pass": False,
                "gates": {"momentum_max": False, "momentum_L2": False},
            },
            sort_keys=True,
        )
        + "\n"
    )
    parent_manifest.write_text(
        json.dumps(
            {
                "candidate_id": "ST052-M",
                "source": {"head": SOURCE_HEAD},
                "replay_identity": {"sha256": "1" * 64},
                "materialized_candidate": {"sha256": _sha(candidate)},
                "validation": {"sha256": _sha(validation)},
                "truth_boundary": {
                    "pde_validated": False,
                    "source_correspondence_verified": False,
                },
            },
            sort_keys=True,
        )
        + "\n"
    )
    return candidate, validation, parent_manifest


def test_bundle_identity_is_deterministic(tmp_path: Path):
    candidate, validation, parent_manifest = _inputs(tmp_path)
    a = build_bundle(
        parent_candidate=candidate,
        parent_validation=validation,
        parent_manifest=parent_manifest,
        out_dir=tmp_path / "a",
    )
    b = build_bundle(
        parent_candidate=candidate,
        parent_validation=validation,
        parent_manifest=parent_manifest,
        out_dir=tmp_path / "b",
    )
    assert a["task_id"] == TASK_ID
    assert a["whole_candidate_identity_sha256"] == b["whole_candidate_identity_sha256"]
    assert verify_bundle(tmp_path / "a")["whole_candidate_identity_sha256"] == a["whole_candidate_identity_sha256"]


def test_bundle_rejects_parent_candidate_not_bound_by_parent_manifest(tmp_path: Path):
    candidate, validation, parent_manifest = _inputs(tmp_path)
    candidate.write_text('{"schema":"tampered-before-bundle"}\n')
    with pytest.raises(ValueError, match="parent candidate bytes"):
        build_bundle(
            parent_candidate=candidate,
            parent_validation=validation,
            parent_manifest=parent_manifest,
            out_dir=tmp_path / "out",
        )


def test_bundle_verifier_rejects_file_tamper(tmp_path: Path):
    candidate, validation, parent_manifest = _inputs(tmp_path)
    build_bundle(
        parent_candidate=candidate,
        parent_validation=validation,
        parent_manifest=parent_manifest,
        out_dir=tmp_path / "bundle",
    )
    target = tmp_path / "bundle" / PARENT_CANDIDATE_FILENAME
    target.write_bytes(target.read_bytes() + b" ")
    with pytest.raises(ValueError, match="file checksum mismatch"):
        verify_bundle(tmp_path / "bundle")


def test_bundle_verifier_rejects_manifest_identity_tamper(tmp_path: Path):
    candidate, validation, parent_manifest = _inputs(tmp_path)
    build_bundle(
        parent_candidate=candidate,
        parent_validation=validation,
        parent_manifest=parent_manifest,
        out_dir=tmp_path / "bundle",
    )
    manifest_path = tmp_path / "bundle" / WHOLE_MANIFEST_FILENAME
    obj = json.loads(manifest_path.read_text())
    obj["whole_candidate_identity_sha256"] = "0" * 64
    manifest_path.write_text(json.dumps(obj, sort_keys=True) + "\n")
    with pytest.raises(ValueError, match="identity checksum mismatch"):
        verify_bundle(tmp_path / "bundle")


def test_bundle_requires_expected_parent_momentum_rejection(tmp_path: Path):
    candidate, validation, parent_manifest = _inputs(tmp_path)
    obj = json.loads(validation.read_text())
    obj["gates"]["momentum_max"] = True
    validation.write_text(json.dumps(obj, sort_keys=True) + "\n")
    pm = json.loads(parent_manifest.read_text())
    pm["validation"]["sha256"] = _sha(validation)
    parent_manifest.write_text(json.dumps(pm, sort_keys=True) + "\n")
    with pytest.raises(ValueError, match="momentum rejection"):
        build_bundle(
            parent_candidate=candidate,
            parent_validation=validation,
            parent_manifest=parent_manifest,
            out_dir=tmp_path / "out",
        )
