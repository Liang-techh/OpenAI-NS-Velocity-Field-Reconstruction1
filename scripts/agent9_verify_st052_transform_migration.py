"""Exact-source parity check for the package-level ST052 #559 transform kernel."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

from openai_ns_reconstruction.st052_local_swirl_transform import (
    DEFAULT_SPEC,
    St052LocalSwirlAdapter,
    TRUTH_BOUNDARY,
)

TASK_ID = "CR-A9-065"
SOURCE_PR = 559
SOURCE_HEAD = "39b106ad8cb8df2064cabead3a12682089575e74"
SOURCE_PARENT_HEAD = "b3b8bfdbe1077f9ec967d158602951997d81e17d"
SOURCE_REDIS_HEAD = "779ffca71066e2864496d37de55a7aafc45d6f57"


def run(source_root: Path, out: Path) -> dict:
    exp = source_root / "experiments" / "root_st052"
    if not exp.is_dir():
        raise FileNotFoundError(f"missing exact source checkout: {exp}")
    sys.path.insert(0, str(exp))
    sys.path.insert(0, str(source_root / "experiments" / "root_st051"))

    import agent7_st052m_frozen_redistribution as base
    import agent7_st052m_local_swirl_energy as comp
    import replay_st052

    field, raw = replay_st052.reconstruct()
    redistribution_scale, _, _ = base.child_scale(field, raw)
    energy = comp.solve_energy_beta(field, raw, redistribution_scale)
    if not energy["root_exists"]:
        raise RuntimeError("exact #559 energy root disappeared")
    beta = float(energy["beta"])
    if abs(redistribution_scale - DEFAULT_SPEC.redistribution_scale) > 2.0e-12:
        raise RuntimeError("redistribution scale identity drift")
    if abs(beta - DEFAULT_SPEC.shoulder_beta) > 2.0e-10:
        raise RuntimeError("shoulder beta identity drift")

    def parent(points, time):
        points = np.asarray(points, dtype=float)
        u, _ = field.fields(raw, points, np.full(len(points), float(time)))
        return np.asarray(u, dtype=float)

    migrated = St052LocalSwirlAdapter(
        parent,
        parent_candidate_id="ST052-M",
        parent_source_head=SOURCE_PARENT_HEAD,
    )

    rng = np.random.default_rng(9176065)
    max_abs = 0.0
    rms_acc = []
    per_time = []
    for time in (0.25, 0.50, 0.75):
        radius = rng.uniform(0.0, 1.95, size=256)
        angle = rng.uniform(0.0, 2.0 * np.pi, size=256)
        z = rng.uniform(-1.95, 1.95, size=256)
        points = np.column_stack((radius * np.cos(angle), radius * np.sin(angle), z))
        expected = comp.compensated_velocity(
            field,
            raw,
            points,
            time,
            redistribution_scale=redistribution_scale,
            beta=beta,
        )
        actual = migrated.points(points, time)
        diff = actual - expected
        this_max = float(np.max(np.abs(diff)))
        this_rms = float(np.sqrt(np.mean(diff * diff)))
        max_abs = max(max_abs, this_max)
        rms_acc.append(this_rms)
        per_time.append({"time": time, "max_abs": this_max, "rms": this_rms})

    parity_tolerance = 5.0e-12
    exact_source_parity = bool(max_abs <= parity_tolerance)
    report = {
        "task_id": TASK_ID,
        "source": {
            "repo": "Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1",
            "source_pr": SOURCE_PR,
            "source_head": SOURCE_HEAD,
            "parent_head": SOURCE_PARENT_HEAD,
            "redistribution_head": SOURCE_REDIS_HEAD,
            "license": "repository-native internal migration; no third-party code copied",
            "migration_class": "direct_internal_migration",
            "migration_scope": "exact deterministic #559 post-processing transform only",
        },
        "transform_spec_sha256": DEFAULT_SPEC.sha256(),
        "source_recomputed": {
            "redistribution_scale": float(redistribution_scale),
            "shoulder_beta": beta,
            "post_transform_common_scale": 1.0,
        },
        "parity": {
            "seed": 9176065,
            "points_per_time": 256,
            "times": [0.25, 0.5, 0.75],
            "max_abs": max_abs,
            "rms_max": float(max(rms_acc)),
            "per_time": per_time,
            "tolerance": parity_tolerance,
            "exact_source_parity": exact_source_parity,
        },
        "direct_contribution": (
            "moves the exact #559 taper/redistribution/shoulder-swirl transform from a stacked "
            "experiment script into a package-level vectorized x/y/z/t adapter with checksumable spec"
        ),
        "remaining_limitations": [
            "the checksum-bound ST052-M parent evaluator itself is not materialized by this increment",
            "this is not a complete candidate save/load capsule until the parent evaluator is bound",
            "pressure and restricted forcing were not rebuilt",
            "fresh held-out CR001 divergence/momentum validation was not run",
            "no OpenAI image-derived numeric target or visual acceptance threshold was used",
        ],
        **TRUTH_BOUNDARY,
    }
    if not exact_source_parity:
        raise RuntimeError(f"exact-source parity failed: max_abs={max_abs}")
    out.parent.mkdir(parents=True, exist_ok=True)
    text = json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n"
    out.write_text(text, encoding="utf-8")
    print(f"exact_source_parity={exact_source_parity}")
    print(f"max_abs={max_abs:.17g}")
    print(f"report_sha256={hashlib.sha256(text.encode()).hexdigest()}")
    return report


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--source-root", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    args = p.parse_args()
    run(args.source_root, args.out)


if __name__ == "__main__":
    main()
