"""Bind held-out validation evidence to the exact public velocity artifact.

This is provenance/integration evidence only. It does not recompute the PDE,
relax thresholds, or infer PDE/visualization acceptance from a validator status.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np


DEFAULT_PROBE_POINTS = np.array(
    [
        [0.0, 0.0, 0.0],
        [0.25, -0.20, 0.30],
        [-0.45, 0.35, -0.55],
        [0.70, 0.10, 0.80],
    ],
    dtype=float,
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _resolve_inside(root: Path, relative: str) -> Path:
    root = root.resolve()
    path = (root / relative).resolve()
    try:
        path.relative_to(root)
    except ValueError as exc:
        raise ValueError("validation candidate path leaves repository root") from exc
    return path


def build_validation_receipt(
    field,
    validation_path,
    *,
    repo_root,
    probe_points=DEFAULT_PROBE_POINTS,
    probe_time: float = 0.5,
) -> dict:
    """Return a deterministic receipt tying validation evidence to ``field``.

    The validation artifact must name its candidate path. That exact file is
    hashed and required to equal the candidate loaded by the public evaluator.
    A small public-API velocity probe is also fingerprinted so downstream
    reports can detect accidental evaluator/candidate drift.

    The receipt deliberately keeps acceptance claims fail-closed: a sampled
    validator status is recorded as evidence but never promoted to
    ``pde_validated=True``.
    """
    validation_path = Path(validation_path)
    repo_root = Path(repo_root)
    _require(validation_path.is_file(), "validation artifact is missing")
    _require(np.isfinite(probe_time), "probe_time must be finite")

    validation = json.loads(validation_path.read_text(encoding="utf-8"))
    _require(isinstance(validation, dict), "validation artifact must contain an object")
    candidate_ref = validation.get("candidate")
    _require(isinstance(candidate_ref, str) and candidate_ref, "validation missing candidate path")
    _require(isinstance(validation.get("status"), str), "validation missing status")
    _require(isinstance(validation.get("seed"), int), "validation missing integer seed")
    _require(
        isinstance(validation.get("points"), int) and validation["points"] > 0,
        "validation missing positive point count",
    )

    candidate_path = _resolve_inside(repo_root, candidate_ref)
    _require(candidate_path.is_file(), "validation candidate file is missing")
    candidate_sha256 = _sha256(candidate_path)

    field_sha256 = getattr(field, "sha256", None)
    _require(
        isinstance(field_sha256, str) and len(field_sha256) == 64,
        "public evaluator must expose candidate sha256",
    )
    _require(
        candidate_sha256 == field_sha256,
        "validation candidate does not match public evaluator candidate",
    )

    metadata = field.metadata()
    _require(isinstance(metadata, dict), "public evaluator metadata must be a mapping")
    _require(
        metadata.get("candidate_sha256") == field_sha256,
        "public evaluator metadata candidate hash mismatch",
    )

    points = np.asarray(probe_points, dtype=float)
    _require(
        points.ndim == 2 and points.shape[1] == 3 and points.shape[0] > 0,
        "probe_points must have shape (N,3)",
    )
    _require(np.all(np.isfinite(points)), "probe_points must be finite")
    velocity = np.asarray(field.at_points(points, float(probe_time)), dtype=float)
    _require(velocity.shape == points.shape, "public evaluator returned wrong probe shape")
    _require(np.all(np.isfinite(velocity)), "public evaluator returned nonfinite probe values")

    probe_hasher = hashlib.sha256()
    probe_hasher.update(np.ascontiguousarray(points, dtype="<f8").tobytes())
    probe_hasher.update(np.asarray([probe_time], dtype="<f8").tobytes())
    probe_hasher.update(np.ascontiguousarray(velocity, dtype="<f8").tobytes())

    return {
        "schema_version": 1,
        "candidate_family": metadata.get("family"),
        "candidate_sha256": candidate_sha256,
        "public_evaluator": f"{type(field).__module__}.{type(field).__qualname__}.at_points",
        "velocity_probe_sha256": probe_hasher.hexdigest(),
        "probe_time": float(probe_time),
        "probe_point_count": int(points.shape[0]),
        "validation_artifact_sha256": _sha256(validation_path),
        "validation_status": validation["status"],
        "validation_seed": int(validation["seed"]),
        "validation_point_count": int(validation["points"]),
        "validation_scope": validation.get("scope"),
        "truth_boundary": {
            "validation_artifact_identity_bound": True,
            "public_velocity_identity_bound": True,
            "velocity_export_ready": "not_assessed_here",
            "visualization_ready": "not_assessed_here",
            "pde_validated": False,
            "paper_exact": False,
            "openai_field_identified": False,
            "note": (
                "This receipt binds existing sampled validation evidence to the "
                "exact public velocity artifact; it does not upgrade acceptance."
            ),
        },
    }
