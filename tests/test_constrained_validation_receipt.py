import hashlib
import json
from pathlib import Path

import numpy as np
import pytest

from openai_ns_reconstruction.constrained_validation_receipt import (
    build_validation_receipt,
)


class FakeField:
    def __init__(self, candidate_path):
        self.path = Path(candidate_path)
        self.sha256 = hashlib.sha256(self.path.read_bytes()).hexdigest()

    def metadata(self):
        return {"family": "fake_v1", "candidate_sha256": self.sha256}

    def at_points(self, points, time):
        p = np.asarray(points, dtype=float)
        return p * (1.0 + float(time))


def make_case(tmp_path, *, status="failed_validation"):
    repo = tmp_path / "repo"
    candidate = repo / "artifacts" / "candidate.json"
    candidate.parent.mkdir(parents=True)
    candidate.write_text('{"family":"fake_v1","a":1}\n', encoding="utf-8")
    validation = repo / "artifacts" / "validation.json"
    validation.write_text(
        json.dumps({
            "status": status,
            "candidate": "artifacts/candidate.json",
            "seed": 914027,
            "points": 128,
            "scope": "held-out sampled evidence",
        }) + "\n",
        encoding="utf-8",
    )
    return repo, candidate, validation, FakeField(candidate)


def test_receipt_is_deterministic_and_identity_bound(tmp_path):
    repo, candidate, validation, field = make_case(tmp_path)
    a = build_validation_receipt(field, validation, repo_root=repo)
    b = build_validation_receipt(field, validation, repo_root=repo)
    assert a == b
    assert a["candidate_sha256"] == hashlib.sha256(candidate.read_bytes()).hexdigest()
    assert a["validation_seed"] == 914027
    assert a["validation_point_count"] == 128
    assert a["truth_boundary"]["public_velocity_identity_bound"] is True


def test_candidate_drift_fails_closed(tmp_path):
    repo, candidate, validation, field = make_case(tmp_path)
    candidate.write_text('{"family":"fake_v1","a":2}\n', encoding="utf-8")
    with pytest.raises(ValueError, match="does not match public evaluator"):
        build_validation_receipt(field, validation, repo_root=repo)


def test_sampled_pass_never_promotes_pde_validated(tmp_path):
    repo, _, validation, field = make_case(
        tmp_path, status="sampled_pde_thresholds_passed"
    )
    report = build_validation_receipt(field, validation, repo_root=repo)
    assert report["validation_status"] == "sampled_pde_thresholds_passed"
    assert report["truth_boundary"]["pde_validated"] is False
    assert report["truth_boundary"]["visualization_ready"] == "not_assessed_here"


def test_candidate_path_cannot_escape_repo(tmp_path):
    repo, _, validation, field = make_case(tmp_path)
    outside = tmp_path / "outside.json"
    outside.write_text('{"family":"fake_v1"}\n', encoding="utf-8")
    data = json.loads(validation.read_text())
    data["candidate"] = "../outside.json"
    validation.write_text(json.dumps(data))
    with pytest.raises(ValueError, match="leaves repository root"):
        build_validation_receipt(field, validation, repo_root=repo)


def test_checked_in_validation_receipt_binds_to_public_velocity():
    from openai_ns_reconstruction.velocity_components import VelocityField

    repo_root = Path(__file__).resolve().parents[1]
    validation = repo_root / "artifacts" / "constrained" / "coupled_joint" / "validation.json"
    report = build_validation_receipt(
        VelocityField(), validation, repo_root=repo_root
    )
    assert report["candidate_sha256"] == "1ab8793073c69ebef59c4fdbd1888e181b1a6bf3aadd6a018eb7b21748dcb65f"
    assert report["validation_status"] == "failed_validation"
    assert report["validation_seed"] == 914027
    assert report["validation_point_count"] == 4096
    assert report["truth_boundary"]["pde_validated"] is False
