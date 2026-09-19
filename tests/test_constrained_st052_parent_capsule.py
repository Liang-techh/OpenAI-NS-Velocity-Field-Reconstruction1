from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

from openai_ns_reconstruction import st052_parent_capsule as capsule


def _git_blob_sha1(payload: bytes) -> str:
    return hashlib.sha1(f"blob {len(payload)}\0".encode("ascii") + payload).hexdigest()


def _fixture(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    source = tmp_path / "source"
    recipe_path = source / capsule.SOURCE_RECIPE_PATH
    recipe_path.parent.mkdir(parents=True)
    recipe = {
        "parent_id": "ST051-B",
        "parent_sha256": "1" * 64,
        "modifiers_sha256": "2" * 64,
        "reference": {
            "points": [[0.1, 0.0, 0.2]],
            "times": [0.5],
            "velocity": [[1.0, 2.0, 3.0]],
            "pressure": [4.0],
        },
        "pde_validated": False,
        "source_correspondence_verified": False,
    }
    payload = json.dumps(recipe, sort_keys=True).encode()
    recipe_path.write_bytes(payload)
    monkeypatch.setattr(capsule, "SOURCE_RECIPE_GIT_BLOB_SHA1", _git_blob_sha1(payload))

    candidate = tmp_path / "candidate.json"
    candidate.write_text('{"synthetic":"candidate"}\n')
    validation = tmp_path / "validation.json"
    validation.write_text(json.dumps({
        "pde_validated": False,
        "all_numeric_gates_pass": False,
        "gates": {"momentum_max": False, "momentum_L2": False},
    }))
    return source, candidate, validation


def test_replay_identity_is_source_recipe_identity_not_candidate_raw_sha(tmp_path, monkeypatch):
    source, candidate, validation = _fixture(tmp_path, monkeypatch)
    first = capsule.build_manifest(source, candidate, validation)
    first_replay = first["replay_identity"]["sha256"]
    first_file = first["materialized_candidate"]["sha256"]

    candidate.write_text('{"synthetic":"same-field-different-metadata"}\n')
    second = capsule.build_manifest(source, candidate, validation)

    assert second["replay_identity"]["sha256"] == first_replay
    assert second["materialized_candidate"]["sha256"] != first_file
    assert second["materialized_candidate"]["historical_raw_sha_identity_claimed"] is False
    assert second["truth_boundary"]["whole_child_save_load_ready"] is False


def test_rejects_pde_validated_recipe(tmp_path, monkeypatch):
    source, candidate, validation = _fixture(tmp_path, monkeypatch)
    recipe_path = source / capsule.SOURCE_RECIPE_PATH
    recipe = json.loads(recipe_path.read_text())
    recipe["pde_validated"] = True
    payload = json.dumps(recipe, sort_keys=True).encode()
    recipe_path.write_bytes(payload)
    monkeypatch.setattr(capsule, "SOURCE_RECIPE_GIT_BLOB_SHA1", _git_blob_sha1(payload))
    with pytest.raises(ValueError, match="recipe.pde_validated"):
        capsule.build_manifest(source, candidate, validation)


def test_rejects_validation_without_expected_momentum_failure(tmp_path, monkeypatch):
    source, candidate, validation = _fixture(tmp_path, monkeypatch)
    report = json.loads(validation.read_text())
    report["gates"]["momentum_max"] = True
    validation.write_text(json.dumps(report))
    with pytest.raises(ValueError, match="momentum rejection"):
        capsule.build_manifest(source, candidate, validation)


def test_rejects_wrong_source_recipe_blob(tmp_path, monkeypatch):
    source, candidate, validation = _fixture(tmp_path, monkeypatch)
    monkeypatch.setattr(capsule, "SOURCE_RECIPE_GIT_BLOB_SHA1", "0" * 40)
    with pytest.raises(ValueError, match="source recipe blob mismatch"):
        capsule.build_manifest(source, candidate, validation)
