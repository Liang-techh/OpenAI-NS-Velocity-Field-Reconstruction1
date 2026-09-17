"""Held-out seed generalization for the supported Phi(1,0) temporal trial.

The morphology-selected affine slope and the restricted-force coefficients are
frozen before this audit.  Fresh deterministic off-grid probe sets are then used
to ask whether the pressure-free vorticity-equation trade-off observed on the
single fixed holdout in the upstream cross-check is a sampling accident.

Every velocity value is obtained through the public ``at_points`` interface.  No
profile coefficient, temporal slope, pressure, force coefficient, derivative
step, viscosity, or acceptance threshold is fitted on these validation seeds.
The report is therefore independent generalization evidence for a visualization
trial, not a full momentum/PDE acceptance test.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
import json
from pathlib import Path
from typing import Iterable

import numpy as np

from .constrained_eq45_supported_force_projection import _curl_force
from .constrained_eq45_supported_phi10_force_crosscheck import (
    FROZEN_SLOPE,
    compare_supported_phi10_temporal_force,
)
from .constrained_eq45_supported_phi10_temporal_replay import build_phi10_temporal_trial
from .constrained_eq45_supported_vorticity import vorticity_equation_terms
from .constrained_force import RestrictedForce

DEFAULT_SEEDS = (914401, 914419, 914437)
DEFAULT_POINTS_PER_REGION = 2
REGIONS = ("plateau", "radial_collar", "axial_collar", "corner_collar")


@dataclass(frozen=True)
class FinestSeedRow:
    seed: int
    static_zero_rms: float
    temporal_zero_rms: float
    zero_fractional_change: float
    static_forced_rms: float
    temporal_forced_rms: float
    forced_fractional_change: float
    velocity_relative_rms_change: float

    def as_dict(self) -> dict[str, float | int]:
        return asdict(self)


def _validated_seeds(seeds: Iterable[int]) -> tuple[int, ...]:
    try:
        values = tuple(seeds)
    except TypeError as exc:
        raise ValueError("seeds must be an iterable of integers") from exc
    if len(values) < 3:
        raise ValueError("at least three held-out seeds are required")
    if any(not isinstance(seed, (int, np.integer)) for seed in values):
        raise ValueError("all seeds must be integers")
    checked = tuple(int(seed) for seed in values)
    if any(seed < 0 for seed in checked):
        raise ValueError("held-out seeds must be nonnegative")
    if len(set(checked)) != len(checked):
        raise ValueError("held-out seeds must be unique")
    return checked


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


def _draw_region(
    rng: np.random.Generator,
    count: int,
    *,
    radial_range: tuple[float, float],
    axial_range: tuple[float, float],
    axial_signed: bool,
) -> np.ndarray:
    radius = rng.uniform(radial_range[0], radial_range[1], size=count)
    angle = rng.uniform(0.0, 2.0 * np.pi, size=count)
    z = rng.uniform(axial_range[0], axial_range[1], size=count)
    if axial_signed:
        z *= rng.choice(np.asarray([-1.0, 1.0]), size=count)
    return np.column_stack((radius * np.cos(angle), radius * np.sin(angle), z))


def held_out_probe_set(seed: int, points_per_region: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Generate a deterministic stratified off-grid validation set.

    The collar samples remain at least 0.22 away from the physical support faces,
    leaving margin for the nested centered stencils at the coarsest registered
    spatial step.  Times are also separated from the delivery interval endpoints.
    """
    if not isinstance(seed, (int, np.integer)) or int(seed) < 0:
        raise ValueError("seed must be a nonnegative integer")
    if not isinstance(points_per_region, int) or points_per_region < 1:
        raise ValueError("points_per_region must be a positive integer")
    rng = np.random.default_rng(int(seed))

    plateau = _draw_region(
        rng,
        points_per_region,
        radial_range=(0.30, 1.20),
        axial_range=(-1.20, 1.20),
        axial_signed=False,
    )
    radial = _draw_region(
        rng,
        points_per_region,
        radial_range=(1.65, 1.78),
        axial_range=(-1.20, 1.20),
        axial_signed=False,
    )
    axial = _draw_region(
        rng,
        points_per_region,
        radial_range=(0.30, 1.20),
        axial_range=(1.65, 1.78),
        axial_signed=True,
    )
    corner = _draw_region(
        rng,
        points_per_region,
        radial_range=(1.65, 1.78),
        axial_range=(1.65, 1.78),
        axial_signed=True,
    )
    points = np.vstack((plateau, radial, axial, corner))
    labels = np.repeat(np.asarray(REGIONS, dtype=object), points_per_region)
    times = rng.uniform(0.34, 0.66, size=len(points))
    return points, times, labels


def _term_scale(terms: dict[str, np.ndarray]) -> float:
    squared = np.zeros(len(np.asarray(terms["residual"])), dtype=float)
    for key in ("omega_t", "advection", "stretching", "diffusion"):
        value = np.asarray(terms[key], dtype=float)
        squared += np.sum(value * value, axis=1)
    return float(np.sqrt(np.mean(squared)))


def _metrics(residual: np.ndarray, scale: float) -> dict[str, float]:
    value = np.asarray(residual, dtype=float)
    if value.ndim != 2 or value.shape[1] != 3 or not np.all(np.isfinite(value)):
        raise ValueError("residual must be a finite (n,3) array")
    norms = np.linalg.norm(value, axis=1)
    rms = float(np.sqrt(np.mean(norms * norms)))
    return {
        "max": float(np.max(norms)),
        "rms": rms,
        "term_normalized_rms": rms / max(float(scale), np.finfo(float).tiny),
    }


def _region_metrics(
    residual: np.ndarray,
    labels: np.ndarray,
    scale: float,
) -> dict[str, dict[str, float]]:
    return {
        region: _metrics(np.asarray(residual)[labels == region], scale)
        for region in REGIONS
    }


def _velocity_relative_rms(static, temporal, points: np.ndarray, times: np.ndarray) -> float:
    base = np.asarray(static.at_points(points, times), dtype=float)
    changed = np.asarray(temporal.at_points(points, times), dtype=float)
    numerator = float(np.sqrt(np.mean(np.sum((changed - base) ** 2, axis=1))))
    denominator = float(np.sqrt(np.mean(np.sum(base**2, axis=1))))
    return numerator / max(denominator, np.finfo(float).tiny)


def audit_supported_phi10_temporal_seed_generalization(
    *,
    seeds: Iterable[int] = DEFAULT_SEEDS,
    points_per_region: int = DEFAULT_POINTS_PER_REGION,
    constraints_path: str | Path | None = None,
) -> dict:
    """Generalize the frozen slope-1.4 trade-off over fresh held-out draws."""
    checked_seeds = _validated_seeds(seeds)
    if not isinstance(points_per_region, int) or points_per_region < 1:
        raise ValueError("points_per_region must be a positive integer")

    root = Path(__file__).resolve().parents[2]
    constraints_path = Path(constraints_path or root / "configs/constraints.json")
    constraints = json.loads(constraints_path.read_text(encoding="utf-8"))
    nu = float(constraints["nu"])
    steps = tuple(float(step) for step in constraints["validation"]["derivative_steps"])
    if len(steps) < 3 or any(right >= left for left, right in zip(steps, steps[1:])):
        raise ValueError("registered derivative ladder must contain >=3 decreasing levels")

    frozen = compare_supported_phi10_temporal_force(
        slope=FROZEN_SLOPE,
        constraints_path=constraints_path,
    )
    trial = build_phi10_temporal_trial(slope=FROZEN_SLOPE)
    static = trial.base
    if static.sha256 != frozen["static_supported_sha256"]:
        raise RuntimeError("static supported candidate identity drifted from the frozen cross-check")
    if trial.sha256 != frozen["temporal_trial_sha256"]:
        raise RuntimeError("temporal trial identity drifted from the frozen cross-check")

    static_fit = frozen["static_projection"]["fit"]
    temporal_fit = frozen["temporal_projection"]["fit"]
    static_force = RestrictedForce(a=float(static_fit["a"]), c=float(static_fit["c"]))
    temporal_force = RestrictedForce(a=float(temporal_fit["a"]), c=float(temporal_fit["c"]))
    time_step = float(temporal_fit["time_step"])

    seed_reports: list[dict[str, object]] = []
    finest_rows: list[FinestSeedRow] = []
    for seed in checked_seeds:
        points, times, labels = held_out_probe_set(seed, points_per_region)
        velocity_change = _velocity_relative_rms(static, trial, points, times)
        level_rows: list[dict[str, object]] = []
        for step in steps:
            static_terms = vorticity_equation_terms(
                static,
                points,
                times,
                spatial_step=step,
                time_step=time_step,
                nu=nu,
            )
            temporal_terms = vorticity_equation_terms(
                trial,
                points,
                times,
                spatial_step=step,
                time_step=time_step,
                nu=nu,
            )
            static_zero = np.asarray(static_terms["residual"], dtype=float)
            temporal_zero = np.asarray(temporal_terms["residual"], dtype=float)
            static_force_curl = _curl_force(static_force, points, times, step)
            temporal_force_curl = _curl_force(temporal_force, points, times, step)
            static_after = static_zero - static_force_curl
            temporal_after = temporal_zero - temporal_force_curl
            static_scale = _term_scale(static_terms)
            temporal_scale = _term_scale(temporal_terms)
            static_zero_metrics = _metrics(static_zero, static_scale)
            temporal_zero_metrics = _metrics(temporal_zero, temporal_scale)
            static_forced_metrics = _metrics(static_after, static_scale)
            temporal_forced_metrics = _metrics(temporal_after, temporal_scale)
            level_rows.append(
                {
                    "spatial_step": step,
                    "static_zero": static_zero_metrics,
                    "temporal_zero": temporal_zero_metrics,
                    "zero_fractional_change": (
                        temporal_zero_metrics["rms"] - static_zero_metrics["rms"]
                    ) / max(static_zero_metrics["rms"], np.finfo(float).tiny),
                    "static_frozen_force": static_forced_metrics,
                    "temporal_frozen_force": temporal_forced_metrics,
                    "forced_fractional_change": (
                        temporal_forced_metrics["rms"] - static_forced_metrics["rms"]
                    ) / max(static_forced_metrics["rms"], np.finfo(float).tiny),
                    "by_region": {
                        "static_zero": _region_metrics(static_zero, labels, static_scale),
                        "temporal_zero": _region_metrics(temporal_zero, labels, temporal_scale),
                        "static_frozen_force": _region_metrics(static_after, labels, static_scale),
                        "temporal_frozen_force": _region_metrics(temporal_after, labels, temporal_scale),
                    },
                }
            )
        finest = level_rows[-1]
        finest_rows.append(
            FinestSeedRow(
                seed=int(seed),
                static_zero_rms=float(finest["static_zero"]["rms"]),
                temporal_zero_rms=float(finest["temporal_zero"]["rms"]),
                zero_fractional_change=float(finest["zero_fractional_change"]),
                static_forced_rms=float(finest["static_frozen_force"]["rms"]),
                temporal_forced_rms=float(finest["temporal_frozen_force"]["rms"]),
                forced_fractional_change=float(finest["forced_fractional_change"]),
                velocity_relative_rms_change=float(velocity_change),
            )
        )
        seed_reports.append(
            {
                "seed": int(seed),
                "probe_count": int(len(points)),
                "points_per_region": int(points_per_region),
                "velocity_relative_rms_change": float(velocity_change),
                "levels": level_rows,
            }
        )

    zero_changes = [row.zero_fractional_change for row in finest_rows]
    forced_changes = [row.forced_fractional_change for row in finest_rows]
    velocity_changes = [row.velocity_relative_rms_change for row in finest_rows]
    return {
        "schema": "eq45_supported_phi10_temporal_seed_generalization_v1",
        "claim_scope": "fresh_seed_generalization_of_frozen_visualization_trial_pde_tradeoff",
        "validation_seeds": list(checked_seeds),
        "seed_count": len(checked_seeds),
        "points_per_region": int(points_per_region),
        "probe_count_per_seed": int(points_per_region * len(REGIONS)),
        "regions": list(REGIONS),
        "sampling_contract": {
            "plateau": "0.30<=r<=1.20, |z|<=1.20",
            "radial_collar": "1.65<=r<=1.78, |z|<=1.20",
            "axial_collar": "0.30<=r<=1.20, 1.65<=|z|<=1.78",
            "corner_collar": "1.65<=r<=1.78, 1.65<=|z|<=1.78",
            "time_range": [0.34, 0.66],
        },
        "nu": nu,
        "spatial_steps": list(steps),
        "fixed_time_step": time_step,
        "slope": FROZEN_SLOPE,
        "mode": {"family": "phi", "index": [1, 0]},
        "static_supported_sha256": static.sha256,
        "temporal_trial_sha256": trial.sha256,
        "slope_frozen_before_validation": True,
        "force_coefficients_frozen_before_validation": True,
        "static_frozen_force": {"a": static_force.a, "c": static_force.c},
        "temporal_frozen_force": {"a": temporal_force.a, "c": temporal_force.c},
        "force_fitted_on_generalization_seeds": False,
        "pressure_fitted": False,
        "new_force_direction_added": False,
        "velocity_access": "public_at_points_only",
        "seed_reports": seed_reports,
        "finest_rows": [row.as_dict() for row in finest_rows],
        "cross_seed_summary": {
            "zero_fractional_change": _summary(zero_changes),
            "forced_fractional_change": _summary(forced_changes),
            "velocity_relative_rms_change": _summary(velocity_changes),
            "zero_improved_seed_fraction": float(np.mean(np.asarray(zero_changes) < 0.0)),
            "forced_improved_seed_fraction": float(np.mean(np.asarray(forced_changes) < 0.0)),
        },
        "registered_pde_thresholds": {
            "max": float(constraints["validation"]["thresholds"]["pde_residual_max"]),
            "L2": float(constraints["validation"]["thresholds"]["pde_residual_L2"]),
        },
        "formal_pde_gate_assessed": False,
        "reason_formal_pde_gate_unassessed": (
            "fresh fixed probes evaluate a pressure-free vorticity necessary condition, not the "
            "preregistered volume-weighted full momentum residual/global maximum"
        ),
        "visualization_candidate_only": True,
        "visualization_ready": False,
        "visual_correspondence_verified": False,
        "pde_validated": False,
        "paper_exact": False,
        "openai_field_identified": False,
        "blowup_proved": False,
    }


def main() -> None:
    print(json.dumps(audit_supported_phi10_temporal_seed_generalization(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
