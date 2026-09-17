"""Cross-seed generalization audit for the frozen Eq. (4.5) vorticity check.

This module deliberately reuses the independent pressure-free vorticity-equation
operator from ``constrained_eq45_vorticity_residual`` and changes only the held-out
probe draw.  Every run reloads the serialized Eq45 candidate, and that upstream
audit obtains velocity values only through public ``at_points(...)->[u,v,w]``.

The purpose is to test whether the already observed O(10) vorticity obstruction is
specific to one deterministic validation seed.  It is not a new PDE acceptance
criterion, optimizer, pressure fit, force fit, or visual-correspondence test.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
import json
from pathlib import Path
from typing import Iterable

import numpy as np

from .constrained_eq45_vorticity_residual import audit_frozen_eq45


DEFAULT_SEEDS = (914211, 914223, 914237, 914251, 914269)


@dataclass(frozen=True)
class SeedGeneralizationRow:
    seed: int
    zero_force_max: float
    zero_force_rms: float
    zero_force_relative_rms: float
    restricted_force_max: float
    restricted_force_rms: float
    restricted_force_relative_rms: float
    force_rms_reduction_fraction: float

    def as_dict(self) -> dict[str, float | int]:
        return asdict(self)


def _validated_seeds(seeds: Iterable[int]) -> tuple[int, ...]:
    try:
        values = tuple(seeds)
    except TypeError as exc:
        raise ValueError("seeds must be an iterable of integers") from exc
    if len(values) < 3:
        raise ValueError("at least three independent held-out seeds are required")
    if any(not isinstance(seed, (int, np.integer)) for seed in values):
        raise ValueError("all seeds must be integers")
    values = tuple(int(seed) for seed in values)
    if len(set(values)) != len(values):
        raise ValueError("held-out seeds must be unique")
    if any(seed < 0 for seed in values):
        raise ValueError("held-out seeds must be nonnegative")
    return values


def _summary(values: list[float]) -> dict[str, float]:
    array = np.asarray(values, dtype=float)
    if array.ndim != 1 or len(array) == 0 or not np.all(np.isfinite(array)):
        raise ValueError("summary values must be a nonempty finite vector")
    mean = float(np.mean(array))
    std = float(np.std(array))
    return {
        "min": float(np.min(array)),
        "median": float(np.median(array)),
        "mean": mean,
        "max": float(np.max(array)),
        "population_std": std,
        "coefficient_of_variation": std / max(abs(mean), np.finfo(float).tiny),
    }


def audit_eq45_vorticity_seed_generalization(
    *,
    seeds: Iterable[int] = DEFAULT_SEEDS,
    sample_count: int = 12,
    candidate_path: str | Path | None = None,
    constraints_path: str | Path | None = None,
) -> dict:
    """Repeat the frozen Eq45 vorticity audit on independent held-out draws.

    ``audit_frozen_eq45`` reloads the candidate on every seed, uses the registered
    three spatial derivative levels, and routes all velocity samples through
    ``Eq45VelocityCandidate.at_points``.  This wrapper records only the finest-level
    result from each draw and summarizes cross-seed spread.
    """
    checked_seeds = _validated_seeds(seeds)
    if not isinstance(sample_count, int) or sample_count < 4:
        raise ValueError("sample_count must be an integer >=4")

    rows: list[SeedGeneralizationRow] = []
    candidate_sha256: str | None = None
    spatial_steps: tuple[float, ...] | None = None
    sample_box = None
    sample_time_range = None
    nu = None
    fixed_time_step = None
    force_choice = None

    for seed in checked_seeds:
        report = audit_frozen_eq45(
            candidate_path=candidate_path,
            constraints_path=constraints_path,
            sample_count=sample_count,
            seed=seed,
        )
        current_sha = str(report["candidate_sha256"])
        if candidate_sha256 is None:
            candidate_sha256 = current_sha
        elif current_sha != candidate_sha256:
            raise RuntimeError("candidate identity changed across held-out seed runs")

        current_steps = tuple(float(step) for step in report["spatial_steps"])
        if spatial_steps is None:
            spatial_steps = current_steps
        elif current_steps != spatial_steps:
            raise RuntimeError("registered derivative levels changed across seed runs")

        if len(report["rows"]) != len(current_steps):
            raise RuntimeError("vorticity audit row count does not match derivative levels")
        finest = report["rows"][-1]
        if float(finest["spatial_step"]) != current_steps[-1]:
            raise RuntimeError("vorticity audit finest row does not match finest derivative step")

        if sample_box is None:
            sample_box = report["sample_box"]
            sample_time_range = report["sample_time_range"]
            nu = float(report["nu"])
            fixed_time_step = float(report["fixed_time_step"])
            force_choice = dict(report["restricted_force_comparison"])

        rows.append(
            SeedGeneralizationRow(
                seed=seed,
                zero_force_max=float(finest["zero_force_max"]),
                zero_force_rms=float(finest["zero_force_rms"]),
                zero_force_relative_rms=float(finest["zero_force_relative_rms"]),
                restricted_force_max=float(finest["restricted_force_max"]),
                restricted_force_rms=float(finest["restricted_force_rms"]),
                restricted_force_relative_rms=float(finest["restricted_force_relative_rms"]),
                force_rms_reduction_fraction=float(finest["force_rms_reduction_fraction"]),
            )
        )

    assert candidate_sha256 is not None and spatial_steps is not None
    return {
        "schema": "eq45_vorticity_seed_generalization_v1",
        "candidate_sha256": candidate_sha256,
        "validation_seeds": list(checked_seeds),
        "seed_count": len(checked_seeds),
        "sample_count_per_seed": sample_count,
        "sample_box": sample_box,
        "sample_time_range": sample_time_range,
        "nu": nu,
        "spatial_steps": list(spatial_steps),
        "reported_spatial_step": spatial_steps[-1],
        "fixed_time_step": fixed_time_step,
        "restricted_force_comparison": force_choice,
        "rows": [row.as_dict() for row in rows],
        "cross_seed_summary": {
            "zero_force_max": _summary([row.zero_force_max for row in rows]),
            "zero_force_rms": _summary([row.zero_force_rms for row in rows]),
            "zero_force_relative_rms": _summary(
                [row.zero_force_relative_rms for row in rows]
            ),
            "restricted_force_max": _summary(
                [row.restricted_force_max for row in rows]
            ),
            "restricted_force_rms": _summary(
                [row.restricted_force_rms for row in rows]
            ),
            "restricted_force_relative_rms": _summary(
                [row.restricted_force_relative_rms for row in rows]
            ),
            "force_rms_reduction_fraction": _summary(
                [row.force_rms_reduction_fraction for row in rows]
            ),
        },
        "interpretation": (
            "held-out probe-seed robustness of the existing pressure-free vorticity "
            "necessary-condition audit; no new PDE acceptance threshold"
        ),
        "velocity_changed": False,
        "pressure_fitted": False,
        "forcing_fitted_in_this_audit": False,
        "training_holdout_separate": True,
        "visualization_ready": False,
        "visual_correspondence_verified": False,
        "pde_validated": False,
        "paper_exact": False,
        "openai_field_identified": False,
        "blowup_proved": False,
    }


def main() -> None:
    print(json.dumps(audit_eq45_vorticity_seed_generalization(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
