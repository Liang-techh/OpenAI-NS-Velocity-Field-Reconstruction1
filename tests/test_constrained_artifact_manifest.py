import json
from pathlib import Path

from openai_ns_reconstruction.constrained_artifact_manifest import git_blob_sha1, validate_manifest


def _write(path: Path, payload):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def _fixture(tmp_path: Path):
    candidate = tmp_path / "artifacts/candidate.json"
    training = tmp_path / "artifacts/training.json"
    validation = tmp_path / "artifacts/validation.json"
    config = tmp_path / "configs/constraints.json"
    _write(candidate, {"schema_version": 1, "status": "candidate", "family": "demo", "parameters": {"a": 1.0}})
    _write(training, {"status": "training_only_not_validated", "seed": 11})
    _write(validation, {
        "status": "failed_validation", "candidate": "artifacts/candidate.json", "seed": 22,
        "rows": [
            {"time": .25, "step": .02, "residual_sampled_max": .2, "divergence_sampled_max": 2e-4},
            {"time": .25, "step": .01, "residual_sampled_max": .19, "divergence_sampled_max": 2e-5},
        ],
    })
    _write(config, {
        "optimization": {"seed": 11},
        "validation": {"seed": 22, "thresholds": {"pde_residual_max": 1e-3, "divergence_max": 1e-5}},
    })
    manifest = tmp_path / "manifest.json"
    payload = {
        "schema_version": 1,
        "artifact_kind": "constrained_reconstruction_reproducibility_manifest",
        "candidate_family": "demo",
        "validation_status": "failed_validation",
        "repository_commit": "deadbeef",
        "truth_boundary": {"paper_exact": False, "full_blowup_proof": False},
        "config_binding": {"status": "exact"},
        "seeds": {"training": 11, "validation": 22},
        "files": {
            "candidate": {"path": "artifacts/candidate.json", "git_blob_sha1": git_blob_sha1(candidate)},
            "training": {"path": "artifacts/training.json", "git_blob_sha1": git_blob_sha1(training)},
            "validation": {"path": "artifacts/validation.json", "git_blob_sha1": git_blob_sha1(validation)},
            "config": {"path": "configs/constraints.json", "git_blob_sha1": git_blob_sha1(config)},
        },
        "observed_finest_validation": {"step": .01, "pde_residual_sampled_max": .19, "divergence_sampled_max": 2e-5},
        "constraint_snapshot": {
            "pde_residual_max": {"status": "fail"},
            "divergence_max": {"status": "fail"},
            "structure": {"status": "pending"},
        },
        "acceptance_ready": False,
    }
    _write(manifest, payload)
    return manifest, payload


def test_failed_candidate_can_have_consistent_reproducibility_manifest(tmp_path):
    manifest, _ = _fixture(tmp_path)
    assert validate_manifest(manifest, tmp_path) == []


def test_hash_mutation_is_detected(tmp_path):
    manifest, _ = _fixture(tmp_path)
    path = tmp_path / "artifacts/candidate.json"
    data = json.loads(path.read_text())
    data["parameters"]["a"] = 2.0
    _write(path, data)
    errors = validate_manifest(manifest, tmp_path)
    assert any("hash mismatch" in error for error in errors)


def test_seed_reuse_and_false_acceptance_fail_closed(tmp_path):
    manifest, payload = _fixture(tmp_path)
    validation = tmp_path / "artifacts/validation.json"
    val = json.loads(validation.read_text())
    val["seed"] = 11
    _write(validation, val)
    payload["files"]["validation"]["git_blob_sha1"] = git_blob_sha1(validation)
    payload["seeds"]["validation"] = 11
    payload["acceptance_ready"] = True
    _write(manifest, payload)
    errors = validate_manifest(manifest, tmp_path)
    assert "training and validation seeds must differ" in errors
    assert "acceptance_ready requires passed_validation" in errors
    assert "acceptance_ready requires every recorded constraint to pass" in errors


def test_checked_in_localized_swirl_manifest_is_consistent():
    repo_root = Path(__file__).resolve().parents[1]
    manifest = repo_root / "artifacts/constrained/localized_swirl/reproducibility_manifest.json"
    required = [
        manifest,
        repo_root / "artifacts/constrained/localized_swirl/candidate.json",
        repo_root / "artifacts/constrained/localized_swirl/training.json",
        repo_root / "artifacts/constrained/localized_swirl/validation.json",
        repo_root / "configs/constraints_localized_swirl.json",
    ]
    if not all(path.exists() for path in required):
        import pytest
        pytest.skip("checked-in constrained artifact bundle not present in isolated unit-test fixture")
    assert validate_manifest(manifest, repo_root) == []
