#!/usr/bin/env python3
"""Build and replay the CR-A9-070 runtime-bound ST052-M temporal capsule."""
from __future__ import annotations

import argparse
import hashlib
import importlib
import json
from pathlib import Path
import sys

import numpy as np

from openai_ns_reconstruction.st052_linear_temporal_capsule import (
    build_bundle,
    load_bundle_runtime,
    verify_bundle,
)
from openai_ns_reconstruction.st052_linear_temporal_transform import (
    DEFAULT_TEMPORAL_SPEC,
    St052LinearTemporalAdapter,
)
from openai_ns_reconstruction.st052_parent_capsule import write_manifest
from openai_ns_reconstruction.st052_source_runtime_identity import (
    authenticate_source_runtime,
    source_runtime_identity_sha256,
)


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _exact_reconstruct(source_root: Path):
    st052 = source_root.resolve() / "experiments/root_st052"
    sys.path.insert(0, str(st052))
    try:
        sys.modules.pop("replay_st052", None)
        replay = importlib.import_module("replay_st052")
        family, raw = replay.reconstruct()
    finally:
        if sys.path and sys.path[0] == str(st052):
            sys.path.pop(0)
    return family, raw


def _base(family, raw):
    def evaluate(points: np.ndarray, time: float) -> np.ndarray:
        points = np.asarray(points, dtype=float)
        times = np.full(len(points), float(time), dtype=float)
        velocity, _ = family.fields(raw, points, times)
        return np.asarray(velocity, dtype=float)
    return evaluate


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--source-root", type=Path, required=True)
    p.add_argument("--parent-dir", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--report", type=Path, required=True)
    a = p.parse_args()

    runtime_receipt = authenticate_source_runtime(a.source_root)
    candidate = a.parent_dir / "candidate.json"
    validation = a.parent_dir / "validation.json"
    parent_manifest = a.parent_dir / "parent_replay_manifest.json"
    write_manifest(a.source_root, candidate, validation, parent_manifest)

    bundle = a.out / "ST052-M-linear-temporal-child-v1"
    manifest = build_bundle(
        parent_candidate=candidate,
        parent_validation=validation,
        parent_manifest=parent_manifest,
        out_dir=bundle,
    )
    verified = verify_bundle(bundle)
    loaded = load_bundle_runtime(bundle, exact_source_root=a.source_root)

    # Authenticate again before the independent direct-source comparator is
    # imported, so both sides of the parity check share the frozen #508 runtime.
    authenticate_source_runtime(a.source_root)
    family, raw = _exact_reconstruct(a.source_root)
    direct = St052LinearTemporalAdapter(
        _base(family, raw),
        parent_candidate_id=DEFAULT_TEMPORAL_SPEC.parent_candidate_id,
        parent_source_head=DEFAULT_TEMPORAL_SPEC.parent_source_head,
        spec=DEFAULT_TEMPORAL_SPEC,
    )

    rng = np.random.default_rng(9069)
    points = rng.uniform((-1.75, -1.75, -1.75), (1.75, 1.75, 1.75), size=(96, 3))
    times = (0.25, 0.375, 0.5, 0.625, 0.75)
    errors = []
    for time in times:
        expected = direct.velocity(points[:, 0], points[:, 1], points[:, 2], time)
        actual = loaded.velocity(points[:, 0], points[:, 1], points[:, 2], time)
        errors.append(float(np.max(np.abs(actual - expected))))
    max_error = max(errors)
    if not np.isfinite(max_error) or max_error > 5.0e-12:
        raise SystemExit(f"whole-child save/load parity failed: {max_error:.17g}")

    validation_obj = json.loads(validation.read_text())
    fine = [row for row in validation_obj.get("spatial_refinement", []) if row.get("space_step") == 0.005]
    momentum_max = max((float(row["momentum_max"]) for row in fine), default=None)
    momentum_l2 = max((float(row["momentum_L2"]) for row in fine), default=None)
    runtime_id = source_runtime_identity_sha256()
    report = {
        "task_id": "CR-A9-070",
        "whole_candidate_identity_sha256": manifest["whole_candidate_identity_sha256"],
        "bundle_manifest_reverified": verified["whole_candidate_identity_sha256"] == manifest["whole_candidate_identity_sha256"],
        "parent_candidate_sha256": _sha(candidate),
        "parent_validation_sha256": _sha(validation),
        "parent_replay_manifest_sha256": _sha(parent_manifest),
        "temporal_transform_spec_sha256": DEFAULT_TEMPORAL_SPEC.sha256(),
        "exact_source_runtime": {
            "identity_sha256": runtime_id,
            "manifest_identity_matches": manifest["runtime"]["exact_source_runtime_identity_sha256"] == runtime_id,
            "authenticated_head": runtime_receipt["authenticated_head"],
            "authenticated_tree": runtime_receipt["authenticated_tree"],
            "tracked_worktree_clean": runtime_receipt["tracked_worktree_clean"],
            "untracked_python_sources_absent": runtime_receipt["untracked_python_sources_absent"],
        },
        "parity": {
            "seed": 9069,
            "points_per_time": int(len(points)),
            "times": list(times),
            "max_component_error_by_time": errors,
            "max_component_error": max_error,
            "threshold": 5.0e-12,
            "passed": True,
        },
        "parent_scientific_receipt": {
            "pde_validated": False,
            "all_numeric_gates_pass": False,
            "momentum_max_gate": validation_obj["gates"]["momentum_max"],
            "momentum_L2_gate": validation_obj["gates"]["momentum_L2"],
            "fine_grid_momentum_max": momentum_max,
            "fine_grid_momentum_L2": momentum_l2,
            "transferred_to_temporal_child": False,
        },
        "truth_boundary": manifest["truth_boundary"],
    }
    a.report.parent.mkdir(parents=True, exist_ok=True)
    a.report.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
