"""Held-out vorticity-equation robustness of the optimized Eq45 profile.

This is a validation-side perturbation audit.  It starts from the checked bounded
Eq45 profile update, makes only small feasible *inward* changes to the two profile
coefficients that were optimized in CR005, and evaluates those changed public
velocity fields with the independent pressure-free vorticity-equation operator.
The restricted force is frozen at the optimized CR005 value for every perturbation;
nothing is refit on validation probes.

The result measures local robustness/sensitivity only.  It cannot establish
Navier--Stokes validity, visual correspondence, paper exactness, or identification
of an OpenAI hidden field.
"""
from __future__ import annotations

from dataclasses import asdict
import json
from pathlib import Path
from typing import Iterable, Mapping

import numpy as np

from .constrained_eq45_candidate import Eq45VelocityCandidate
from .constrained_eq45_optimized_vorticity_crosscheck import (
    DEFAULT_TIME_STEP,
    DISALLOWED_VALIDATION_SEEDS as PRIOR_DISALLOWED_SEEDS,
    SAMPLE_BOX,
    TIME_RANGE,
)
from .constrained_eq45_profile_force_optimization import candidate_with_profile_pair
from .constrained_eq45_vorticity_residual import audit_vorticity_equation
from .constrained_force import RestrictedForce

DEFAULT_VALIDATION_SEEDS = (914351, 914363, 914377)
DEFAULT_SAMPLE_COUNT = 8
DEFAULT_AMPLITUDES = (0.02, 0.04)
DIRECTIONS: Mapping[str, tuple[float, float]] = {
    "phi_inward": (-1.0, 0.0),
    "F_inward": (0.0, -1.0),
    "joint_inward": (-2.0**-0.5, -2.0**-0.5),
}
DISALLOWED_VALIDATION_SEEDS = set(PRIOR_DISALLOWED_SEEDS) | {914307, 914319, 914333}


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


def _validated_amplitudes(amplitudes: Iterable[float]) -> tuple[float, ...]:
    try:
        values = tuple(float(value) for value in amplitudes)
    except (TypeError, ValueError) as exc:
        raise ValueError("amplitudes must be finite positive numbers") from exc
    if len(values) < 2:
        raise ValueError("at least two perturbation amplitudes are required")
    if any(not np.isfinite(value) or value <= 0.0 for value in values):
        raise ValueError("amplitudes must be finite and positive")
    if any(next_value <= value for value, next_value in zip(values, values[1:])):
        raise ValueError("amplitudes must be strictly increasing")
    return values


def _probe_cloud(seed: int, sample_count: int) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    points = rng.uniform(SAMPLE_BOX[0], SAMPLE_BOX[1], size=(sample_count, 3))
    times = rng.uniform(TIME_RANGE[0], TIME_RANGE[1], size=sample_count)
    return points, times


def _vector_rms(values: np.ndarray) -> float:
    array = np.asarray(values, dtype=float)
    if array.ndim < 1 or array.shape[-1] != 3 or not np.all(np.isfinite(array)):
        raise ValueError("values must be a finite vector field with final dimension 3")
    return float(np.sqrt(np.mean(np.sum(array * array, axis=-1))))


def _velocity_response(
    baseline: Eq45VelocityCandidate,
    trial: Eq45VelocityCandidate,
    points: np.ndarray,
    times: np.ndarray,
) -> dict[str, float]:
    base_velocity = np.asarray(baseline.at_points(points, times), dtype=float)
    trial_velocity = np.asarray(trial.at_points(points, times), dtype=float)
    if base_velocity.shape != points.shape or trial_velocity.shape != points.shape:
        raise RuntimeError("public velocity evaluator returned an unexpected shape")
    delta = trial_velocity - base_velocity
    baseline_rms = _vector_rms(base_velocity)
    delta_rms = _vector_rms(delta)
    delta_norms = np.linalg.norm(delta, axis=-1)
    return {
        "baseline_velocity_rms": baseline_rms,
        "delta_velocity_rms": delta_rms,
        "delta_velocity_relative_rms": delta_rms / max(baseline_rms, np.finfo(float).tiny),
        "delta_velocity_max": float(np.max(delta_norms)),
    }


def _summary(values: Iterable[float]) -> dict[str, float]:
    array = np.asarray(tuple(values), dtype=float)
    if array.ndim != 1 or len(array) == 0 or not np.all(np.isfinite(array)):
        raise ValueError("summary values must be a nonempty finite vector")
    return {
        "min": float(np.min(array)),
        "median": float(np.median(array)),
        "mean": float(np.mean(array)),
        "max": float(np.max(array)),
        "population_std": float(np.std(array)),
    }


def audit_optimized_eq45_parameter_robustness(
    *,
    seeds: Iterable[int] = DEFAULT_VALIDATION_SEEDS,
    amplitudes: Iterable[float] = DEFAULT_AMPLITUDES,
    sample_count: int = DEFAULT_SAMPLE_COUNT,
    candidate_path: str | Path | None = None,
    optimization_path: str | Path | None = None,
    constraints_path: str | Path | None = None,
) -> dict:
    """Measure independent PDE sensitivity around the bounded optimized profile."""
    checked_seeds = _validated_seeds(seeds)
    checked_amplitudes = _validated_amplitudes(amplitudes)
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

    seed_candidate = Eq45VelocityCandidate.load_json(candidate_path)
    optimization = json.loads(optimization_path.read_text(encoding="utf-8"))
    constraints = json.loads(constraints_path.read_text(encoding="utf-8"))
    if optimization.get("task_id") != "CR005-EQ45-PROFILE-FORCE-OPT-015":
        raise RuntimeError("optimization artifact task identity changed")

    fit = optimization["fit"]
    optimized = candidate_with_profile_pair(
        seed_candidate,
        float(fit["Phi_02"]),
        float(fit["F_02"]),
    )
    expected_sha = str(optimization["optimized_candidate"]["sha256"])
    if optimized.sha256 != expected_sha:
        raise RuntimeError("optimized candidate identity does not match the checked receipt")

    coefficient_limit = float(optimized.profile_basis.coefficient_limit)
    if not np.isfinite(coefficient_limit) or coefficient_limit <= 0.0:
        raise RuntimeError("candidate coefficient limit is invalid")
    baseline_parameters = np.array([float(fit["Phi_02"]), float(fit["F_02"])])
    frozen_force = RestrictedForce(a=float(fit["force_a"]), c=float(fit["force_c"]))

    nu = float(constraints["nu"])
    spatial_steps = tuple(float(value) for value in constraints["validation"]["derivative_steps"])
    if len(spatial_steps) < 3:
        raise RuntimeError("registered validation ladder must contain at least three levels")
    if any(next_step >= step for step, next_step in zip(spatial_steps, spatial_steps[1:])):
        raise RuntimeError("registered derivative steps must be strictly decreasing")

    perturbations: list[tuple[str, float, np.ndarray, Eq45VelocityCandidate]] = []
    for direction_name, raw_direction in DIRECTIONS.items():
        direction = np.asarray(raw_direction, dtype=float)
        if not np.isclose(np.linalg.norm(direction), 1.0, rtol=0.0, atol=1e-14):
            raise RuntimeError("perturbation directions must have unit Euclidean norm")
        for amplitude in checked_amplitudes:
            parameters = baseline_parameters + amplitude * direction
            if np.any(parameters < -coefficient_limit) or np.any(parameters > coefficient_limit):
                raise ValueError(
                    f"perturbation {direction_name}@{amplitude} leaves registered coefficient bounds"
                )
            trial = candidate_with_profile_pair(optimized, parameters[0], parameters[1])
            perturbations.append((direction_name, amplitude, parameters, trial))

    seed_reports: list[dict] = []
    for seed in checked_seeds:
        points, times = _probe_cloud(seed, sample_count)
        baseline_rows = audit_vorticity_equation(
            optimized,
            points,
            times,
            spatial_steps=spatial_steps,
            time_step=DEFAULT_TIME_STEP,
            nu=nu,
            restricted_force=frozen_force,
        )
        if len(baseline_rows) != len(spatial_steps):
            raise RuntimeError("baseline vorticity audit row count changed unexpectedly")

        trial_reports: list[dict] = []
        for direction_name, amplitude, parameters, trial in perturbations:
            trial_rows = audit_vorticity_equation(
                trial,
                points,
                times,
                spatial_steps=spatial_steps,
                time_step=DEFAULT_TIME_STEP,
                nu=nu,
                restricted_force=frozen_force,
            )
            if len(trial_rows) != len(spatial_steps):
                raise RuntimeError("perturbed vorticity audit row count changed unexpectedly")
            levels: list[dict] = []
            for base_row, trial_row in zip(baseline_rows, trial_rows):
                if base_row.spatial_step != trial_row.spatial_step:
                    raise RuntimeError("baseline and perturbation derivative ladders differ")
                levels.append(
                    {
                        "spatial_step": float(base_row.spatial_step),
                        "baseline": asdict(base_row),
                        "perturbed": asdict(trial_row),
                        "zero_force_rms_relative_change": float(
                            trial_row.zero_force_rms / max(base_row.zero_force_rms, np.finfo(float).tiny)
                            - 1.0
                        ),
                        "frozen_force_rms_relative_change": float(
                            trial_row.restricted_force_rms
                            / max(base_row.restricted_force_rms, np.finfo(float).tiny)
                            - 1.0
                        ),
                    }
                )
            trial_reports.append(
                {
                    "direction": direction_name,
                    "amplitude": float(amplitude),
                    "Phi_02": float(parameters[0]),
                    "F_02": float(parameters[1]),
                    "candidate_sha256": trial.sha256,
                    "velocity_response": _velocity_response(optimized, trial, points, times),
                    "levels": levels,
                }
            )
        seed_reports.append(
            {
                "seed": int(seed),
                "baseline_levels": [asdict(row) for row in baseline_rows],
                "perturbations": trial_reports,
            }
        )

    finest_changes_zero: list[float] = []
    finest_changes_frozen: list[float] = []
    finest_baseline_zero: list[float] = []
    finest_baseline_frozen: list[float] = []
    velocity_relative_changes: list[float] = []
    for seed_report in seed_reports:
        finest_baseline_zero.append(seed_report["baseline_levels"][-1]["zero_force_rms"])
        finest_baseline_frozen.append(seed_report["baseline_levels"][-1]["restricted_force_rms"])
        for trial_report in seed_report["perturbations"]:
            finest = trial_report["levels"][-1]
            finest_changes_zero.append(finest["zero_force_rms_relative_change"])
            finest_changes_frozen.append(finest["frozen_force_rms_relative_change"])
            velocity_relative_changes.append(
                trial_report["velocity_response"]["delta_velocity_relative_rms"]
            )

    return {
        "schema": "eq45_optimized_parameter_robustness_v1",
        "task_id": "CR009-EQ45-OPTIMIZED-PARAM-PERTURB-VORTICITY-016",
        "source_optimization_task_id": str(optimization["task_id"]),
        "optimized_candidate_sha256": optimized.sha256,
        "baseline_profile_parameters": {
            "Phi_02": float(baseline_parameters[0]),
            "F_02": float(baseline_parameters[1]),
        },
        "coefficient_bound": [-coefficient_limit, coefficient_limit],
        "perturbation_directions": {
            key: [float(value) for value in direction] for key, direction in DIRECTIONS.items()
        },
        "perturbation_amplitudes": list(checked_amplitudes),
        "validation_seeds": list(checked_seeds),
        "sample_count_per_seed": sample_count,
        "sample_box": [[SAMPLE_BOX[0], SAMPLE_BOX[1]]] * 3,
        "sample_time_range": [TIME_RANGE[0], TIME_RANGE[1]],
        "nu": nu,
        "spatial_steps": list(spatial_steps),
        "fixed_time_step": DEFAULT_TIME_STEP,
        "frozen_force": {
            "a": frozen_force.a,
            "c": frozen_force.c,
            "refit_on_validation": False,
        },
        "seed_reports": seed_reports,
        "finest_level_summary": {
            "optimized_zero_force_rms_across_seeds": _summary(finest_baseline_zero),
            "optimized_frozen_force_rms_across_seeds": _summary(finest_baseline_frozen),
            "perturbed_zero_force_rms_relative_change": _summary(finest_changes_zero),
            "perturbed_frozen_force_rms_relative_change": _summary(finest_changes_frozen),
            "public_velocity_relative_rms_change": _summary(velocity_relative_changes),
        },
        "interpretation": (
            "local held-out parameter-robustness audit around the bounded optimized Eq45 profile; "
            "a small residual change under perturbation is robustness evidence only, while a large "
            "change diagnoses local fragility; neither outcome is PDE acceptance"
        ),
        "training_performed_here": False,
        "pressure_fitted_here": False,
        "force_refit_here": False,
        "velocity_changed_by_audit": False,
        "visualization_candidate_only": True,
        "visualization_ready": False,
        "visual_correspondence_verified": False,
        "pde_validated": False,
        "paper_exact": False,
        "openai_field_identified": False,
        "blowup_proved": False,
    }


def main() -> None:
    print(json.dumps(audit_optimized_eq45_parameter_robustness(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
