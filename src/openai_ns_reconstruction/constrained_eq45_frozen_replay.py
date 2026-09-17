"""Deterministic replay receipt for the first frozen Eq. (4.5) velocity artifact.

This is a CR011 reproducibility layer, not a scientific acceptance test.  It
loads the checked Eq45 candidate artifact, consumes velocity only through its
public ``at_points`` interface, and fingerprints fixed-seed off-grid probes plus
small reference grids.  The receipt keeps delivery, visualization, support and
PDE states separate.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np

from .constrained_eq45_candidate import Eq45VelocityCandidate


SCHEMA = "eq45_frozen_replay_receipt_v1"
EXPECTED_CANDIDATE_SHA256 = "48f1845fcfd71ec95495c98bb7aac3fca4653e748a8856a5421234e9601525a7"
DEFAULT_PROBE_SEED = 914117
DEFAULT_PROBE_COUNT = 64
DEFAULT_PROBE_HALF_WIDTH = 0.75
DEFAULT_GRID_SIZE = 5
REFERENCE_TIMES = (0.25, 0.50, 0.75)


def _fingerprint_arrays(*arrays: np.ndarray) -> str:
    digest = hashlib.sha256()
    for value in arrays:
        array = np.ascontiguousarray(np.asarray(value, dtype="<f8"))
        digest.update(str(array.shape).encode("ascii"))
        digest.update(array.tobytes(order="C"))
    return digest.hexdigest()


def replay_eq45_candidate(
    candidate_path: str | Path,
    *,
    expected_candidate_sha256: str = EXPECTED_CANDIDATE_SHA256,
    probe_seed: int = DEFAULT_PROBE_SEED,
    probe_count: int = DEFAULT_PROBE_COUNT,
    probe_half_width: float = DEFAULT_PROBE_HALF_WIDTH,
    grid_size: int = DEFAULT_GRID_SIZE,
) -> dict[str, Any]:
    """Replay one frozen candidate through its public velocity interface.

    The default expected hash intentionally pins this receipt to the first
    serialized Eq45 seed.  A deliberate candidate update must therefore update
    the replay contract rather than silently inheriting old evidence.
    """
    path = Path(candidate_path)
    candidate = Eq45VelocityCandidate.load_json(path)

    if not isinstance(expected_candidate_sha256, str) or len(expected_candidate_sha256) != 64:
        raise ValueError("expected_candidate_sha256 must be a 64-character hex digest")
    try:
        int(expected_candidate_sha256, 16)
    except ValueError as exc:
        raise ValueError("expected_candidate_sha256 must be hexadecimal") from exc
    if candidate.sha256 != expected_candidate_sha256:
        raise ValueError("candidate SHA256 does not match frozen replay contract")

    if not isinstance(probe_seed, int):
        raise ValueError("probe_seed must be an integer")
    if not isinstance(probe_count, int) or probe_count < 8:
        raise ValueError("probe_count must be an integer >= 8")
    if not np.isfinite(probe_half_width) or not (0.0 < probe_half_width <= 1.0):
        raise ValueError("probe_half_width must lie in (0, 1]")
    if not isinstance(grid_size, int) or grid_size < 3 or grid_size % 2 == 0:
        raise ValueError("grid_size must be an odd integer >= 3")

    rng = np.random.default_rng(probe_seed)
    points = rng.uniform(-probe_half_width, probe_half_width, size=(probe_count, 3))
    times = rng.uniform(candidate.time_start, candidate.time_end, size=probe_count)
    velocity = np.asarray(candidate.at_points(points, times), dtype=float)
    if velocity.shape != points.shape or not np.all(np.isfinite(velocity)):
        raise RuntimeError("public Eq45 velocity replay returned malformed/nonfinite values")

    speed_rms = float(np.sqrt(np.mean(np.sum(velocity * velocity, axis=-1))))
    if not np.isfinite(speed_rms) or speed_rms <= 1e-12:
        raise RuntimeError("frozen Eq45 replay is numerically inactive")

    axis = np.linspace(-probe_half_width, probe_half_width, grid_size)
    xx, yy, zz = np.meshgrid(axis, axis, axis, indexing="ij")
    grid_points = np.stack((xx, yy, zz), axis=-1).reshape(-1, 3)
    grid_values = []
    for time in REFERENCE_TIMES:
        if not (candidate.time_start <= time <= candidate.time_end):
            raise RuntimeError("reference time fell outside candidate delivery interval")
        values = np.asarray(candidate.at_points(grid_points, time), dtype=float)
        if values.shape != grid_points.shape or not np.all(np.isfinite(values)):
            raise RuntimeError("public Eq45 grid replay returned malformed/nonfinite values")
        grid_values.append(values)
    grid_velocity = np.stack(grid_values, axis=0)

    receipt = {
        "schema": SCHEMA,
        "candidate": {
            "path": path.as_posix(),
            "sha256": candidate.sha256,
            "evaluator": "Eq45VelocityCandidate.at_points",
            "time_interval": [float(candidate.time_start), float(candidate.time_end)],
        },
        "replay": {
            "probe_seed": probe_seed,
            "probe_count": probe_count,
            "probe_half_width": float(probe_half_width),
            "probe_fingerprint_sha256": _fingerprint_arrays(points, times, velocity),
            "probe_speed_rms": speed_rms,
            "grid_size": grid_size,
            "reference_times": list(REFERENCE_TIMES),
            "grid_fingerprint_sha256": _fingerprint_arrays(axis, np.asarray(REFERENCE_TIMES), grid_velocity),
        },
        "states": {
            "velocity_export_ready": True,
            "visualization_ready": False,
            "physical_support_validated": False,
            "pde_validated": False,
            "visual_correspondence_verified": False,
            "paper_exact": False,
            "openai_field_identified": False,
            "blowup_proved": False,
        },
        "claim_scope": "reproducible_velocity_replay_only",
    }
    return receipt


def save_replay_receipt(receipt: dict[str, Any], path: str | Path) -> None:
    if receipt.get("schema") != SCHEMA:
        raise ValueError("receipt schema does not match Eq45 frozen replay contract")
    states = receipt.get("states")
    if not isinstance(states, dict):
        raise ValueError("receipt states are missing")
    forbidden_true = (
        "visualization_ready",
        "physical_support_validated",
        "pde_validated",
        "visual_correspondence_verified",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    )
    if any(states.get(name) is not False for name in forbidden_true):
        raise ValueError("replay receipt attempted to promote an unsupported scientific state")
    if states.get("velocity_export_ready") is not True:
        raise ValueError("replay receipt lost the candidate-local export-ready state")
    Path(path).write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
