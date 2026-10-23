"""Truth-bounded replay of the currently packaged constrained velocity delivery.

This module is intentionally small: it checks that the selected optimizer output,
the packaged public velocity artifact, the independent validation metadata, and
the public grid evaluator still describe one coherent frozen candidate.  It also
runs a small independent PDE-operator smoke through the public velocity API.

A successful replay means the current delivery chain is reproducible in the
runtime executing this command.  It does not promote visualization
correspondence, PDE validation, paper exactness, or a blow-up claim.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np

from .constrained_coupled import CoupledCandidate
from .constrained_force import RestrictedForce
from .constrained_validation import residual, sampled_norms
from .velocity_components import VelocityField


_SELECTED = Path("artifacts/constrained/coupled_joint/candidate.json")
_TRAINING = Path("artifacts/constrained/coupled_joint/training.json")
_VALIDATION = Path("artifacts/constrained/coupled_joint/validation.json")
_CONFIG = Path("configs/constraints.json")
_PACKAGED = Path("src/openai_ns_reconstruction/data/velocity_candidate.json")
_REPLAY_SEED = 914061


def _read_json(path: Path) -> dict[str, Any]:
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise ValueError(f"expected JSON object: {path}")
    return data


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _array_sha256(*arrays: np.ndarray) -> str:
    digest = hashlib.sha256()
    for array in arrays:
        stable = np.ascontiguousarray(np.asarray(array, dtype="<f8"))
        digest.update(stable.tobytes())
    return digest.hexdigest()


def replay_current_delivery(repo_root: str | Path | None = None) -> dict[str, Any]:
    """Replay the frozen public velocity chain without changing the candidate.

    The replay seed is deliberately distinct from both training and preregistered
    validation seeds.  The PDE calculation is a finite smoke only; the stored
    held-out validation artifact remains authoritative for the current failure
    status.
    """
    root = (
        Path(repo_root).resolve()
        if repo_root is not None
        else Path(__file__).resolve().parents[2]
    )
    paths = {
        "selected_candidate": root / _SELECTED,
        "packaged_candidate": root / _PACKAGED,
        "training": root / _TRAINING,
        "validation": root / _VALIDATION,
        "config": root / _CONFIG,
    }
    missing = [name for name, path in paths.items() if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"missing replay inputs: {', '.join(missing)}")

    selected_bytes = paths["selected_candidate"].read_bytes()
    packaged_bytes = paths["packaged_candidate"].read_bytes()
    if selected_bytes != packaged_bytes:
        raise ValueError("selected and packaged candidate bytes differ")

    selected_payload = _read_json(paths["selected_candidate"])
    training = _read_json(paths["training"])
    validation = _read_json(paths["validation"])
    config = _read_json(paths["config"])

    if selected_payload.get("family") != "coupled_velocity_v1":
        raise ValueError("unexpected selected candidate family")
    if validation.get("candidate") != _SELECTED.as_posix():
        raise ValueError("validation artifact does not point at selected candidate")
    if training.get("force") != validation.get("force"):
        raise ValueError("training and validation force metadata differ")

    training_seed = int(training["seed"])
    validation_seed = int(validation["seed"])
    if training_seed == validation_seed:
        raise ValueError("training and validation seeds must be distinct")
    if int(config["optimization"]["seed"]) != training_seed:
        raise ValueError("training seed does not match preregistered optimization seed")
    if int(config["validation"]["seed"]) != validation_seed:
        raise ValueError("validation seed does not match preregistered validation seed")
    if _REPLAY_SEED in {training_seed, validation_seed}:
        raise ValueError("replay seed must be distinct from training and validation")

    field = VelocityField(paths["packaged_candidate"])
    metadata = field.metadata()
    candidate_sha256 = _sha256(paths["selected_candidate"])
    if field.sha256 != candidate_sha256:
        raise ValueError("public velocity candidate SHA256 mismatch")
    if metadata.get("candidate_sha256") != candidate_sha256:
        raise ValueError("public velocity metadata candidate SHA256 mismatch")
    if metadata.get("family") != selected_payload["family"]:
        raise ValueError("public velocity family does not match selected candidate")

    box = np.asarray(config["domain"]["evaluation_box"], dtype=float)
    if box.shape != (3, 2) or not np.all(np.isfinite(box)) or not np.all(box[:, 1] > box[:, 0]):
        raise ValueError("invalid configured evaluation box")
    time_bounds = tuple(float(value) for value in config["domain"]["time_interval"])
    if len(time_bounds) != 2 or not time_bounds[0] < time_bounds[1]:
        raise ValueError("invalid configured time interval")

    rng = np.random.default_rng(_REPLAY_SEED)
    probe_lo = 0.5 * box[:, 0]
    probe_hi = 0.5 * box[:, 1]
    probe_points = rng.uniform(probe_lo, probe_hi, size=(32, 3))
    probe_time = float(sum(time_bounds) / 2.0)
    probe_velocity = np.asarray(field.at_points(probe_points, probe_time), dtype=float)
    if probe_velocity.shape != (32, 3) or not np.all(np.isfinite(probe_velocity)):
        raise ValueError("public velocity probe returned malformed values")
    probe_speed_rms = float(np.sqrt(np.mean(np.sum(probe_velocity * probe_velocity, axis=1))))
    if not probe_speed_rms > 0.0:
        raise ValueError("public velocity probe is numerically trivial")

    candidate = CoupledCandidate.load(paths["selected_candidate"])
    force = RestrictedForce(**validation["force"])
    pde_result = residual(
        lambda points, time: field.at_points(points, time),
        candidate.pressure,
        force,
        probe_points,
        probe_time,
        nu=float(config["nu"]),
        step=0.01,
        time_bounds=time_bounds,
    )
    volume = float(np.prod(box[:, 1] - box[:, 0]))
    pde_smoke = sampled_norms(pde_result, volume)
    if not all(np.isfinite(value) for value in pde_smoke.values()):
        raise ValueError("independent public-velocity PDE smoke returned non-finite metrics")

    axes = np.linspace(box[:, 0], box[:, 1], 5, axis=1)
    times = np.asarray([time_bounds[0], probe_time, time_bounds[1]], dtype=float)
    grid = np.asarray(field.grid(axes[0], axes[1], axes[2], times), dtype=float)
    if grid.shape != (3, 5, 5, 5, 3) or not np.all(np.isfinite(grid)):
        raise ValueError("public velocity grid replay returned malformed values")

    return {
        "schema": "constrained_delivery_replay_v1",
        "replay_pass": True,
        "candidate_family": selected_payload["family"],
        "candidate_sha256": candidate_sha256,
        "training_seed": training_seed,
        "validation_seed": validation_seed,
        "replay_seed": _REPLAY_SEED,
        "validation_status": validation.get("status"),
        "validation_points": int(validation["points"]),
        "probe_time": probe_time,
        "probe_points": int(probe_points.shape[0]),
        "probe_speed_rms": probe_speed_rms,
        "public_probe_sha256": _array_sha256(probe_points, probe_velocity),
        "pde_operator_smoke": pde_smoke,
        "reference_times": [float(value) for value in times],
        "grid_shape": list(grid.shape),
        "grid_values_sha256": _array_sha256(axes, times, grid),
        "delivery_state": {
            "velocity_export_ready": True,
            "visualization_ready": False,
            "pde_validated": False,
            "paper_exact": False,
            "openai_correspondence_verified": False,
        },
        "truth_boundary": (
            "fresh runtime replay of one frozen candidate through the public velocity "
            "API and independent finite-difference smoke only; the stored held-out "
            "validation status is not upgraded"
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = replay_current_delivery(args.repo_root)
    text = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding="utf-8")
    print(text, end="")


if __name__ == "__main__":
    main()
