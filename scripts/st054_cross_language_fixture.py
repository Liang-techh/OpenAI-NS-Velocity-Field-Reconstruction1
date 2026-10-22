#!/usr/bin/env python3
"""Generate a frozen off-grid ST054 Python velocity fixture for native MATLAB parity.

This is software/cross-language evidence only.  It does not evaluate Navier--Stokes
acceptance or visual correspondence and does not fit any field parameter.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from openai_ns_reconstruction import available_st054_models, load_st054
from openai_ns_reconstruction.st054_snapshot import ST054_REFERENCE_ATOL

SEED = 9174091
POINT_COUNT = 96
TIMES = np.asarray([0.271, 0.389, 0.503, 0.617, 0.739], dtype=float)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _offgrid_points() -> np.ndarray:
    rng = np.random.default_rng(SEED)
    accepted: list[np.ndarray] = []
    while sum(chunk.shape[0] for chunk in accepted) < POINT_COUNT:
        trial = rng.uniform([-1.55, -1.55, -1.70], [1.55, 1.55, 1.70], size=(256, 3))
        radius = np.hypot(trial[:, 0], trial[:, 1])
        mask = (radius > 0.075) & (radius < 1.70) & (np.abs(trial[:, 2]) < 1.70)
        accepted.append(trial[mask])
    points = np.concatenate(accepted, axis=0)[:POINT_COUNT]
    if not np.isfinite(points).all():
        raise RuntimeError("Generated non-finite cross-language points")
    return points


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    points = _offgrid_points()
    files: dict[str, dict[str, object]] = {}
    mat_sha: str | None = None
    for model_id in available_st054_models():
        field = load_st054(model_id)
        if mat_sha is None:
            mat_sha = field.mat_sha256
        elif field.mat_sha256 != mat_sha:
            raise RuntimeError("ST054 models unexpectedly resolve to different MAT payloads")
        rows: list[np.ndarray] = []
        for time in TIMES:
            velocity = np.asarray(field.velocity(points[:, 0], points[:, 1], points[:, 2], float(time)))
            if velocity.shape != (POINT_COUNT, 3) or not np.isfinite(velocity).all():
                raise RuntimeError(f"Malformed Python velocity output for {model_id}")
            rows.append(
                np.column_stack(
                    (
                        points,
                        np.full(POINT_COUNT, float(time)),
                        velocity,
                    )
                )
            )
        table = np.concatenate(rows, axis=0)
        filename = model_id.lower().replace("-", "_") + ".csv"
        path = args.output / filename
        np.savetxt(path, table, delimiter=",", fmt="%.17g")
        files[model_id] = {
            "file": filename,
            "rows": int(table.shape[0]),
            "sha256": _sha256(path),
        }

    receipt = {
        "schema": "st054_python_matlab_offgrid_parity_fixture_v1",
        "seed": SEED,
        "points_per_time": POINT_COUNT,
        "times": TIMES.tolist(),
        "models": list(available_st054_models()),
        "columns": ["x", "y", "z", "t", "u", "v", "w"],
        "engineering_atol": ST054_REFERENCE_ATOL,
        "engineering_scope": "cross_language_velocity_parity_only_not_pde_or_visual_acceptance",
        "source_mat_sha256": mat_sha,
        "files": files,
        "pde_validated": False,
        "visual_correspondence_verified": False,
        "openai_field_identified": False,
    }
    manifest = args.output / "fixture_manifest.json"
    manifest.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(receipt, sort_keys=True))


if __name__ == "__main__":
    main()
