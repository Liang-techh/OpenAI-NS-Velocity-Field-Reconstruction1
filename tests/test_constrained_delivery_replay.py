import json
from pathlib import Path
import shutil

import numpy as np
import pytest

from openai_ns_reconstruction.constrained_delivery_replay import replay_current_delivery


REQUIRED = (
    "artifacts/constrained/coupled_joint/candidate.json",
    "artifacts/constrained/coupled_joint/training.json",
    "artifacts/constrained/coupled_joint/validation.json",
    "configs/constraints.json",
    "src/openai_ns_reconstruction/data/velocity_candidate.json",
)


def _copy_replay_inputs(source_root: Path, target_root: Path) -> None:
    for relative in REQUIRED:
        source = source_root / relative
        target = target_root / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, target)


def test_checked_in_delivery_replays_through_public_velocity_and_independent_operator():
    repo_root = Path(__file__).resolve().parents[1]
    report = replay_current_delivery(repo_root)

    assert report["replay_pass"] is True
    assert report["candidate_family"] == "coupled_velocity_v1"
    assert report["candidate_sha256"] == "1ab8793073c69ebef59c4fdbd1888e181b1a6bf3aadd6a018eb7b21748dcb65f"
    assert report["training_seed"] == 20260916
    assert report["validation_seed"] == 914027
    assert report["replay_seed"] not in {report["training_seed"], report["validation_seed"]}
    assert report["validation_status"] == "failed_validation"
    assert report["validation_points"] == 4096
    assert report["grid_shape"] == [3, 5, 5, 5, 3]
    assert report["probe_speed_rms"] > 0.0
    assert len(report["public_probe_sha256"]) == 64
    assert len(report["grid_values_sha256"]) == 64
    assert all(np.isfinite(value) for value in report["pde_operator_smoke"].values())
    assert report["delivery_state"] == {
        "velocity_export_ready": True,
        "visualization_ready": False,
        "pde_validated": False,
        "paper_exact": False,
        "openai_correspondence_verified": False,
    }


def test_replay_fails_closed_if_packaged_candidate_drifted(tmp_path):
    repo_root = Path(__file__).resolve().parents[1]
    _copy_replay_inputs(repo_root, tmp_path)
    packaged = tmp_path / "src/openai_ns_reconstruction/data/velocity_candidate.json"
    packaged.write_bytes(packaged.read_bytes() + b"\n")

    with pytest.raises(ValueError, match="selected and packaged candidate bytes differ"):
        replay_current_delivery(tmp_path)


def test_replay_fails_closed_on_training_validation_seed_reuse(tmp_path):
    repo_root = Path(__file__).resolve().parents[1]
    _copy_replay_inputs(repo_root, tmp_path)
    validation_path = tmp_path / "artifacts/constrained/coupled_joint/validation.json"
    validation = json.loads(validation_path.read_text(encoding="utf-8"))
    validation["seed"] = 20260916
    validation_path.write_text(json.dumps(validation), encoding="utf-8")

    with pytest.raises(ValueError, match="training and validation seeds must be distinct"):
        replay_current_delivery(tmp_path)
