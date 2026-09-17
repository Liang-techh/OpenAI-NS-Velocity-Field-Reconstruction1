"""Independent vorticity-equation cross-check for the bounded Eq45 profile update.

This validator consumes two checked artifacts: the frozen Eq45 seed candidate and
Agent 2's bounded profile/force optimization report.  It reconstructs exactly the
reported diagnostic candidate, verifies its serialized identity, and then evaluates
an independently implemented pressure-free vorticity equation on new held-out probe
draws.  No pressure, force, or profile parameter is fitted here.

The check is a necessary-condition diagnostic only.  Residual reduction does not
promote the field to PDE-valid, visualization-verified, paper-exact, or an identified
OpenAI hidden field.
"""
from __future__ import annotations

from dataclasses import asdict
import json
from pathlib import Path
from typing import Iterable

import numpy as np

from .constrained_eq45_candidate import Eq45VelocityCandidate
from .constrained_eq45_profile_force_optimization import candidate_with_profile_pair
from .constrained_eq45_vorticity_residual import audit_vorticity_equation
from .constrained_force import RestrictedForce

DEFAULT_VALIDATION_SEEDS = (914307, 914319, 914333)
DEFAULT_SAMPLE_COUNT = 12
DEFAULT_TIME_STEP = 0.0025
SAMPLE_BOX = (-0.22, 0.22)
TIME_RANGE = (0.40, 0.60)
# Seeds already used by the training/holdout optimization or earlier Agent-3 audits.
DISALLOWED_VALIDATION_SEEDS = {
    20260917,
    914117,
    914131,
    914211,
    914223,
    914237,
    914251,
    914269,
}


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
    overlap = set(values) & DISALLOWED_VALIDATION_SEEDS
    if overlap:
        raise ValueError(f"validation seeds overlap prior fit/holdout draws: {sorted(overlap)}")
    return values


def _summary(values: Iterable[float]) -> dict[str, float]:
    array = np.asarray(tuple(values), dtype=float)
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


def _probe_cloud(seed: int, sample_count: int) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    points = rng.uniform(SAMPLE_BOX[0], SAMPLE_BOX[1], size=(sample_count, 3))
    times = rng.uniform(TIME_RANGE[0], TIME_RANGE[1], size=sample_count)
    return points, times


def audit_optimized_eq45_vorticity(
    *,
    seeds: Iterable[int] = DEFAULT_VALIDATION_SEEDS,
    sample_count: int = DEFAULT_SAMPLE_COUNT,
    candidate_path: str | Path | None = None,
    optimization_path: str | Path | None = None,
    constraints_path: str | Path | None = None,
) -> dict:
    """Cross-check the frozen optimized profile with a different PDE identity.

    Every field evaluation is delegated to ``Eq45VelocityCandidate.at_points`` by
    ``audit_vorticity_equation``.  This function only freezes inputs and summarizes
    the independent held-out results.
    """
    checked_seeds = _validated_seeds(seeds)
    if not isinstance(sample_count, int) or sample_count < 4:
        raise ValueError("sample_count must be an integer >=4")

    root = Path(__file__).resolve().parents[2]
    candidate_path = Path(
        candidate_path or root / "artifacts/constrained/eq45_velocity_candidate_seed.json"
    )
    optimization_path = Path(
        optimization_path or root / "artifacts/constrained/eq45_profile_force_optimization.json"
    )
    constraints_path = Path(constraints_path or root / "configs/constraints.json")

    baseline = Eq45VelocityCandidate.load_json(candidate_path)
    optimization = json.loads(optimization_path.read_text(encoding="utf-8"))
    constraints = json.loads(constraints_path.read_text(encoding="utf-8"))

    fit = optimization["fit"]
    seed_fit = optimization["seed"]
    optimized = candidate_with_profile_pair(
        baseline,
        float(fit["Phi_02"]),
        float(fit["F_02"]),
    )
    roundtrip = Eq45VelocityCandidate.from_dict(optimized.to_dict())
    expected_sha = str(optimization["optimized_candidate"]["sha256"])
    if optimized.sha256 != expected_sha or roundtrip.sha256 != expected_sha:
        raise RuntimeError("optimized candidate identity does not match the checked optimization artifact")

    nu = float(constraints["nu"])
    spatial_steps = tuple(float(x) for x in constraints["validation"]["derivative_steps"])
    if len(spatial_steps) < 3:
        raise RuntimeError("registered validation ladder must contain at least three levels")

    baseline_force = RestrictedForce(
        a=float(seed_fit["force_a"]),
        c=float(seed_fit["force_c"]),
    )
    optimized_force = RestrictedForce(
        a=float(fit["force_a"]),
        c=float(fit["force_c"]),
    )

    seed_reports: list[dict] = []
    for seed in checked_seeds:
        points, times = _probe_cloud(seed, sample_count)
        baseline_rows = audit_vorticity_equation(
            baseline,
            points,
            times,
            spatial_steps=spatial_steps,
            time_step=DEFAULT_TIME_STEP,
            nu=nu,
            restricted_force=baseline_force,
        )
        optimized_rows = audit_vorticity_equation(
            roundtrip,
            points,
            times,
            spatial_steps=spatial_steps,
            time_step=DEFAULT_TIME_STEP,
            nu=nu,
            restricted_force=optimized_force,
        )
        if len(baseline_rows) != len(spatial_steps) or len(optimized_rows) != len(spatial_steps):
            raise RuntimeError("vorticity audit row count changed unexpectedly")

        levels: list[dict] = []
        for baseline_row, optimized_row in zip(baseline_rows, optimized_rows):
            if baseline_row.spatial_step != optimized_row.spatial_step:
                raise RuntimeError("baseline and optimized derivative ladders differ")
            zero_reduction = 1.0 - optimized_row.zero_force_rms / max(
                baseline_row.zero_force_rms, np.finfo(float).tiny
            )
            restricted_reduction = 1.0 - optimized_row.restricted_force_rms / max(
                baseline_row.restricted_force_rms, np.finfo(float).tiny
            )
            levels.append(
                {
                    "spatial_step": baseline_row.spatial_step,
                    "baseline": asdict(baseline_row),
                    "optimized": asdict(optimized_row),
                    "profile_only_zero_force_rms_reduction_fraction": float(zero_reduction),
                    "profile_plus_frozen_force_rms_reduction_fraction": float(
                        restricted_reduction
                    ),
                }
            )
        seed_reports.append({"seed": seed, "levels": levels})

    finest = [report["levels"][-1] for report in seed_reports]
    zero_reductions = [
        row["profile_only_zero_force_rms_reduction_fraction"] for row in finest
    ]
    restricted_reductions = [
        row["profile_plus_frozen_force_rms_reduction_fraction"] for row in finest
    ]
    return {
        "schema": "eq45_optimized_profile_vorticity_crosscheck_v1",
        "task_id": "CR009-EQ45-OPTIMIZED-PROFILE-VORTICITY-CROSSCHECK-015",
        "baseline_candidate_sha256": baseline.sha256,
        "optimized_candidate_sha256": optimized.sha256,
        "optimized_candidate_roundtrip_equal": roundtrip.sha256 == optimized.sha256,
        "source_optimization_task_id": str(optimization["task_id"]),
        "profile_parameters": {
            "Phi_02": float(fit["Phi_02"]),
            "F_02": float(fit["F_02"]),
        },
        "frozen_force_choices": {
            "baseline": {"a": baseline_force.a, "c": baseline_force.c},
            "optimized": {"a": optimized_force.a, "c": optimized_force.c},
            "fitted_in_this_audit": False,
        },
        "validation_seeds": list(checked_seeds),
        "sample_count_per_seed": sample_count,
        "sample_box": [[SAMPLE_BOX[0], SAMPLE_BOX[1]]] * 3,
        "sample_time_range": [TIME_RANGE[0], TIME_RANGE[1]],
        "nu": nu,
        "spatial_steps": list(spatial_steps),
        "fixed_time_step": DEFAULT_TIME_STEP,
        "seed_reports": seed_reports,
        "finest_level_summary": {
            "baseline_zero_force_rms": _summary(
                row["baseline"]["zero_force_rms"] for row in finest
            ),
            "optimized_zero_force_rms": _summary(
                row["optimized"]["zero_force_rms"] for row in finest
            ),
            "profile_only_zero_force_rms_reduction_fraction": _summary(zero_reductions),
            "baseline_restricted_force_rms": _summary(
                row["baseline"]["restricted_force_rms"] for row in finest
            ),
            "optimized_restricted_force_rms": _summary(
                row["optimized"]["restricted_force_rms"] for row in finest
            ),
            "profile_plus_frozen_force_rms_reduction_fraction": _summary(
                restricted_reductions
            ),
        },
        "interpretation": (
            "independent pressure-free vorticity-equation held-out cross-check of the "
            "bounded profile update; reductions are validation evidence only and not a "
            "Navier-Stokes acceptance result"
        ),
        "velocity_changed_by_source_optimization": True,
        "velocity_changed_by_this_validator": False,
        "training_or_force_fit_performed_here": False,
        "visualization_ready": False,
        "visual_correspondence_verified": False,
        "pde_validated": False,
        "paper_exact": False,
        "openai_field_identified": False,
        "blowup_proved": False,
    }


def main() -> None:
    print(json.dumps(audit_optimized_eq45_vorticity(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
