#!/usr/bin/env python3
"""Export and verify the CR-A9-071 ST052 whole-child visualization grid."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from openai_ns_reconstruction.st052_grid_export import (
    DEFAULT_GRID_POINTS,
    DEFAULT_TIMES,
    MANIFEST_FILENAME,
    NPZ_FILENAME,
    export_velocity_grid,
    verify_velocity_grid_export,
)
from openai_ns_reconstruction.st052_linear_temporal_capsule import load_bundle_runtime


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--source-root", type=Path, required=True)
    p.add_argument("--bundle-dir", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--report", type=Path, required=True)
    p.add_argument("--grid-points", type=int, default=DEFAULT_GRID_POINTS)
    a = p.parse_args()

    candidate = load_bundle_runtime(a.bundle_dir, exact_source_root=a.source_root)
    manifest = export_velocity_grid(
        candidate,
        a.out,
        grid_points=a.grid_points,
        times=DEFAULT_TIMES,
    )
    checked = verify_velocity_grid_export(a.out)

    # Independent direct-node replay against the loaded runtime.  This is not a
    # visual acceptance test; it only checks that the serialized grid contains
    # the values returned by the authenticated whole-candidate callable.
    with np.load(a.out / NPZ_FILENAME, allow_pickle=False) as data:
        x = np.asarray(data["x"], dtype=float)
        y = np.asarray(data["y"], dtype=float)
        z = np.asarray(data["z"], dtype=float)
        times = np.asarray(data["t"], dtype=float)
        components = [np.asarray(data[name], dtype=float) for name in ("u", "v", "w")]
    node_indices = (0, len(x) // 4, len(x) // 2, (3 * len(x)) // 4, len(x) - 1)
    max_error = 0.0
    comparisons = 0
    for ti, time in enumerate(times):
        for ix in node_indices:
            iy = (2 * ix + ti) % len(y)
            iz = (3 * ix + 2 * ti) % len(z)
            actual = np.asarray(candidate.velocity(x[ix], y[iy], z[iz], float(time)), dtype=float).reshape(3)
            saved = np.array([c[ti, ix, iy, iz] for c in components], dtype=float)
            max_error = max(max_error, float(np.max(np.abs(actual - saved))))
            comparisons += 1
    if not np.isfinite(max_error) or max_error != 0.0:
        raise SystemExit(f"serialized grid differs from callable at frozen nodes: {max_error:.17g}")

    report = {
        "task_id": "CR-A9-071",
        "whole_candidate_identity_sha256": candidate.identity_sha256,
        "candidate_id": candidate.candidate_id,
        "grid_payload_sha256": manifest["grid_payload_sha256"],
        "grid": manifest["grid"],
        "files": manifest["files"],
        "manifest_sha256": _sha(a.out / MANIFEST_FILENAME),
        "npz_mat_roundtrip": {
            "exact_equal": checked["npz_mat_exact_equal"],
            "npz_payload_sha256": checked["npz_payload_sha256"],
            "mat_payload_sha256": checked["mat_payload_sha256"],
        },
        "callable_grid_node_replay": {
            "comparisons": comparisons,
            "max_component_error": max_error,
            "passed": True,
        },
        "external_method": manifest["external_method"],
        "execution_environment": manifest["execution_environment"],
        "matlab_compatibility": manifest["matlab_compatibility"],
        "truth_boundary": manifest["truth_boundary"],
    }
    a.report.parent.mkdir(parents=True, exist_ok=True)
    a.report.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
