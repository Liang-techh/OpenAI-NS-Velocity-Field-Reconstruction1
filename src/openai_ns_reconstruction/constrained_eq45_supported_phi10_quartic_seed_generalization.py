"""Fresh-seed PDE generalization for the derivative-balanced quartic Phi(1,0) trial.

The quartic schedule and the static/cubic/quartic RestrictedForce(a,c) projections
are reproduced from the already-frozen CR005 cross-check before this audit sees any
fresh seed.  The three candidate objects are then serialized and reloaded; every
production velocity value used below is obtained only through public ``at_points``.

This module changes no velocity coefficient, force coefficient, pressure model,
derivative step, viscosity, sampling threshold, or acceptance threshold.  It asks
one CR009 question: does the small fixed-holdout quartic-vs-cubic vorticity-residual
advantage generalize to new stratified off-grid draws, or is it sampling-sensitive?
The result is a pressure-free necessary-condition diagnostic for a visualization
candidate, not the preregistered full-momentum PDE acceptance gate.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Iterable

import numpy as np

from .constrained_eq45_candidate import Eq45VelocityCandidate
from .constrained_eq45_supported_candidate import Eq45SupportedVelocityCandidate
from .constrained_eq45_supported_force_projection import _curl_force
from .constrained_eq45_supported_phi10_cubic_temporal_mode import (
    Eq45SupportedPhi10CubicLocalizedTemporalCandidate,
)
from .constrained_eq45_supported_phi10_quartic_force_crosscheck import (
    FROZEN_EARLY_DELTA,
    FROZEN_NULLSPACE_COEFFICIENT,
    compare_quartic_balanced_phi10_force,
)
from .constrained_eq45_supported_phi10_quartic_temporal_mode import (
    Eq45SupportedPhi10QuarticDerivativeBalancedTemporalCandidate,
)
from .constrained_eq45_supported_vorticity import vorticity_equation_terms
from .constrained_force import RestrictedForce

DEFAULT_SEEDS = (914531, 914547, 914563)
DEFAULT_POINTS_PER_REGION = 2
REGIONS = ("plateau", "radial_collar", "axial_collar", "corner_collar")


@dataclass(frozen=True)
class FinestSeedRow:
    seed: int
    static_zero_rms: float
    cubic_zero_rms: float
    quartic_zero_rms: float
    quartic_vs_static_zero_fractional_change: float
    quartic_vs_cubic_zero_fractional_change: float
    static_forced_rms: float
    cubic_forced_rms: float
    quartic_forced_rms: float
    quartic_vs_static_forced_fractional_change: float
    quartic_vs_cubic_forced_fractional_change: float
    quartic_vs_static_velocity_relative_rms: float
    quartic_vs_cubic_velocity_relative_rms: float
    quartic_vs_cubic_zero_sign_consistent_across_levels: bool
    quartic_vs_cubic_forced_sign_consistent_across_levels: bool

    def as_dict(self) -> dict[str, float | int | bool]:
        return asdict(self)


def _validated_seeds(seeds: Iterable[int]) -> tuple[int, ...]:
    try:
        raw = tuple(seeds)
    except TypeError as exc:
        raise ValueError("seeds must be an iterable of integers") from exc
    if len(raw) < 3:
        raise ValueError("at least three fresh held-out seeds are required")
    if any(not isinstance(seed, (int, np.integer)) for seed in raw):
        raise ValueError("all seeds must be integers")
    checked = tuple(int(seed) for seed in raw)
    if any(seed < 0 for seed in checked):
        raise ValueError("seeds must be nonnegative")
    if len(set(checked)) != len(checked):
        raise ValueError("fresh held-out seeds must be unique")
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
    signed_axial: bool,
) -> np.ndarray:
    radius = rng.uniform(radial_range[0], radial_range[1], size=count)
    angle = rng.uniform(0.0, 2.0 * np.pi, size=count)
    z = rng.uniform(axial_range[0], axial_range[1], size=count)
    if signed_axial:
        z *= rng.choice(np.asarray([-1.0, 1.0]), size=count)
    return np.column_stack((radius * np.cos(angle), radius * np.sin(angle), z))


def held_out_probe_set(
    seed: int, points_per_region: int
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Generate one deterministic off-grid draw across plateau and taper collars."""
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
        signed_axial=False,
    )
    radial = _draw_region(
        rng,
        points_per_region,
        radial_range=(1.65, 1.78),
        axial_range=(-1.20, 1.20),
        signed_axial=False,
    )
    axial = _draw_region(
        rng,
        points_per_region,
        radial_range=(0.30, 1.20),
        axial_range=(1.65, 1.78),
        signed_axial=True,
    )
    corner = _draw_region(
        rng,
        points_per_region,
        radial_range=(1.65, 1.78),
        axial_range=(1.65, 1.78),
        signed_axial=True,
    )
    points = np.vstack((plateau, radial, axial, corner))
    labels = np.repeat(np.asarray(REGIONS, dtype=object), points_per_region)
    # Stay away from temporal endpoints so all candidates use the same centered
    # time-difference stencil in the inherited independent vorticity operator.
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
    if value.ndim != 2 or value.shape[1] != 3 or len(value) == 0:
        raise ValueError("residual must have shape (n,3) with n>0")
    if not np.all(np.isfinite(value)):
        raise FloatingPointError("residual must be finite")
    norms = np.linalg.norm(value, axis=1)
    rms = float(np.sqrt(np.mean(norms * norms)))
    return {
        "max": float(np.max(norms)),
        "rms": rms,
        "term_normalized_rms": rms / max(float(scale), np.finfo(float).tiny),
    }


def _region_metrics(
    residual: np.ndarray, labels: np.ndarray, scale: float
) -> dict[str, dict[str, float]]:
    return {
        region: _metrics(np.asarray(residual)[labels == region], scale)
        for region in REGIONS
    }


def _fractional_change(new: float, old: float) -> float:
    return float((new - old) / max(abs(old), np.finfo(float).tiny))


def _velocity_relative_rms(first, second, points: np.ndarray, times: np.ndarray) -> float:
    baseline = np.asarray(first.at_points(points, times), dtype=float)
    changed = np.asarray(second.at_points(points, times), dtype=float)
    numerator = float(np.sqrt(np.mean(np.sum((changed - baseline) ** 2, axis=1))))
    denominator = float(np.sqrt(np.mean(np.sum(baseline * baseline, axis=1))))
    return numerator / max(denominator, np.finfo(float).tiny)


def _same_nonzero_sign(values: Iterable[float]) -> bool:
    array = np.asarray(tuple(values), dtype=float)
    if array.ndim != 1 or len(array) == 0 or not np.all(np.isfinite(array)):
        raise ValueError("sign values must be a nonempty finite vector")
    signs = np.sign(array)
    return bool(signs[0] != 0.0 and np.all(signs == signs[0]))


def _serialized_candidates():
    """Construct then reload all three public velocity objects before validation."""
    constructed_static = Eq45SupportedVelocityCandidate(parent=Eq45VelocityCandidate.seed())
    constructed_cubic = Eq45SupportedPhi10CubicLocalizedTemporalCandidate(
        base=constructed_static, early_delta=FROZEN_EARLY_DELTA
    )
    constructed_quartic = Eq45SupportedPhi10QuarticDerivativeBalancedTemporalCandidate(
        base=constructed_static, early_delta=FROZEN_EARLY_DELTA
    )
    with TemporaryDirectory() as tmp:
        root = Path(tmp)
        static_path = root / "static.json"
        cubic_path = root / "cubic.json"
        quartic_path = root / "quartic.json"
        constructed_static.save_json(static_path)
        constructed_cubic.save_json(cubic_path)
        constructed_quartic.save_json(quartic_path)
        static = Eq45SupportedVelocityCandidate.load_json(static_path)
        cubic = Eq45SupportedPhi10CubicLocalizedTemporalCandidate.load_json(cubic_path)
        quartic = Eq45SupportedPhi10QuarticDerivativeBalancedTemporalCandidate.load_json(
            quartic_path
        )
    if static.sha256 != constructed_static.sha256:
        raise AssertionError("static supported serialization changed identity")
    if cubic.sha256 != constructed_cubic.sha256:
        raise AssertionError("cubic serialization changed identity")
    if quartic.sha256 != constructed_quartic.sha256:
        raise AssertionError("quartic serialization changed identity")
    return static, cubic, quartic


def audit_quartic_phi10_seed_generalization(
    *,
    seeds: Iterable[int] = DEFAULT_SEEDS,
    points_per_region: int = DEFAULT_POINTS_PER_REGION,
    constraints_path: str | Path | None = None,
) -> dict:
    """Audit frozen quartic-vs-cubic PDE ordering on fresh deterministic seeds."""
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

    # Reproduce the already-frozen force cross-check before touching fresh seeds.
    frozen = compare_quartic_balanced_phi10_force(constraints_path=constraints_path)
    static, cubic, quartic = _serialized_candidates()
    if static.sha256 != frozen["static_supported_sha256"]:
        raise RuntimeError("static identity drifted from frozen quartic cross-check")
    if cubic.sha256 != frozen["cubic_trial_sha256"]:
        raise RuntimeError("cubic identity drifted from frozen quartic cross-check")
    if quartic.sha256 != frozen["quartic_trial_sha256"]:
        raise RuntimeError("quartic identity drifted from frozen quartic cross-check")
    if quartic.nullspace_coefficient != FROZEN_NULLSPACE_COEFFICIENT:
        raise RuntimeError("quartic target-free balancing coefficient drifted")

    static_fit = frozen["static_projection"]["fit"]
    cubic_fit = frozen["cubic_projection"]["fit"]
    quartic_fit = frozen["quartic_projection"]["fit"]
    static_force = RestrictedForce(a=float(static_fit["a"]), c=float(static_fit["c"]))
    cubic_force = RestrictedForce(a=float(cubic_fit["a"]), c=float(cubic_fit["c"]))
    quartic_force = RestrictedForce(a=float(quartic_fit["a"]), c=float(quartic_fit["c"]))
    time_step = float(quartic_fit["time_step"])
    if float(static_fit["time_step"]) != time_step or float(cubic_fit["time_step"]) != time_step:
        raise RuntimeError("frozen force cross-check time-step contract drifted")

    seed_reports: list[dict[str, object]] = []
    finest_rows: list[FinestSeedRow] = []
    for seed in checked_seeds:
        points, times, labels = held_out_probe_set(seed, points_per_region)
        quartic_vs_static_velocity = _velocity_relative_rms(static, quartic, points, times)
        quartic_vs_cubic_velocity = _velocity_relative_rms(cubic, quartic, points, times)
        level_rows: list[dict[str, object]] = []
        for step in steps:
            candidates = {
                "static": (static, static_force),
                "cubic": (cubic, cubic_force),
                "quartic": (quartic, quartic_force),
            }
            values: dict[str, dict[str, object]] = {}
            for name, (candidate, force) in candidates.items():
                terms = vorticity_equation_terms(
                    candidate,
                    points,
                    times,
                    spatial_step=step,
                    time_step=time_step,
                    nu=nu,
                )
                zero = np.asarray(terms["residual"], dtype=float)
                forced = zero - _curl_force(force, points, times, step)
                scale = _term_scale(terms)
                values[name] = {
                    "zero": _metrics(zero, scale),
                    "frozen_force": _metrics(forced, scale),
                    "by_region_zero": _region_metrics(zero, labels, scale),
                    "by_region_frozen_force": _region_metrics(forced, labels, scale),
                }

            static_zero = float(values["static"]["zero"]["rms"])
            cubic_zero = float(values["cubic"]["zero"]["rms"])
            quartic_zero = float(values["quartic"]["zero"]["rms"])
            static_forced = float(values["static"]["frozen_force"]["rms"])
            cubic_forced = float(values["cubic"]["frozen_force"]["rms"])
            quartic_forced = float(values["quartic"]["frozen_force"]["rms"])
            level_rows.append(
                {
                    "spatial_step": step,
                    "static": values["static"],
                    "cubic": values["cubic"],
                    "quartic": values["quartic"],
                    "quartic_vs_static_zero_fractional_change": _fractional_change(
                        quartic_zero, static_zero
                    ),
                    "quartic_vs_cubic_zero_fractional_change": _fractional_change(
                        quartic_zero, cubic_zero
                    ),
                    "quartic_vs_static_forced_fractional_change": _fractional_change(
                        quartic_forced, static_forced
                    ),
                    "quartic_vs_cubic_forced_fractional_change": _fractional_change(
                        quartic_forced, cubic_forced
                    ),
                }
            )

        zero_vs_cubic = [
            float(row["quartic_vs_cubic_zero_fractional_change"]) for row in level_rows
        ]
        forced_vs_cubic = [
            float(row["quartic_vs_cubic_forced_fractional_change"])
            for row in level_rows
        ]
        finest = level_rows[-1]
        finest_rows.append(
            FinestSeedRow(
                seed=int(seed),
                static_zero_rms=float(finest["static"]["zero"]["rms"]),
                cubic_zero_rms=float(finest["cubic"]["zero"]["rms"]),
                quartic_zero_rms=float(finest["quartic"]["zero"]["rms"]),
                quartic_vs_static_zero_fractional_change=float(
                    finest["quartic_vs_static_zero_fractional_change"]
                ),
                quartic_vs_cubic_zero_fractional_change=float(
                    finest["quartic_vs_cubic_zero_fractional_change"]
                ),
                static_forced_rms=float(finest["static"]["frozen_force"]["rms"]),
                cubic_forced_rms=float(finest["cubic"]["frozen_force"]["rms"]),
                quartic_forced_rms=float(finest["quartic"]["frozen_force"]["rms"]),
                quartic_vs_static_forced_fractional_change=float(
                    finest["quartic_vs_static_forced_fractional_change"]
                ),
                quartic_vs_cubic_forced_fractional_change=float(
                    finest["quartic_vs_cubic_forced_fractional_change"]
                ),
                quartic_vs_static_velocity_relative_rms=float(quartic_vs_static_velocity),
                quartic_vs_cubic_velocity_relative_rms=float(quartic_vs_cubic_velocity),
                quartic_vs_cubic_zero_sign_consistent_across_levels=_same_nonzero_sign(
                    zero_vs_cubic
                ),
                quartic_vs_cubic_forced_sign_consistent_across_levels=_same_nonzero_sign(
                    forced_vs_cubic
                ),
            )
        )
        seed_reports.append(
            {
                "seed": int(seed),
                "probe_count": int(len(points)),
                "points_per_region": int(points_per_region),
                "quartic_vs_static_velocity_relative_rms": float(quartic_vs_static_velocity),
                "quartic_vs_cubic_velocity_relative_rms": float(quartic_vs_cubic_velocity),
                "levels": level_rows,
            }
        )

    zero_vs_static = [row.quartic_vs_static_zero_fractional_change for row in finest_rows]
    zero_vs_cubic = [row.quartic_vs_cubic_zero_fractional_change for row in finest_rows]
    forced_vs_static = [row.quartic_vs_static_forced_fractional_change for row in finest_rows]
    forced_vs_cubic = [row.quartic_vs_cubic_forced_fractional_change for row in finest_rows]
    return {
        "schema": "eq45_supported_phi10_quartic_seed_generalization_v1",
        "task_id": "CR009-EQ45-SUPPORTED-PHI10-QUARTIC-SEED-GENERALIZATION-028",
        "claim_scope": "fresh_seed_generalization_of_frozen_quartic_vs_cubic_pde_tradeoff",
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
        "early_delta": FROZEN_EARLY_DELTA,
        "nullspace_coefficient": FROZEN_NULLSPACE_COEFFICIENT,
        "mode": {"family": "phi", "index": [1, 0]},
        "static_supported_sha256": static.sha256,
        "cubic_trial_sha256": cubic.sha256,
        "quartic_trial_sha256": quartic.sha256,
        "candidates_serialized_and_reloaded_before_validation": True,
        "schedule_frozen_before_validation": True,
        "force_coefficients_frozen_before_validation": True,
        "static_frozen_force": {"a": static_force.a, "c": static_force.c},
        "cubic_frozen_force": {"a": cubic_force.a, "c": cubic_force.c},
        "quartic_frozen_force": {"a": quartic_force.a, "c": quartic_force.c},
        "fresh_seed_force_refit": False,
        "pressure_fitted": False,
        "new_force_direction_added": False,
        "residual_defined_force_allowed": False,
        "seed_reports": seed_reports,
        "finest_rows": [row.as_dict() for row in finest_rows],
        "finest_cross_seed_summary": {
            "quartic_vs_static_zero_fractional_change": _summary(zero_vs_static),
            "quartic_vs_cubic_zero_fractional_change": _summary(zero_vs_cubic),
            "quartic_vs_static_forced_fractional_change": _summary(forced_vs_static),
            "quartic_vs_cubic_forced_fractional_change": _summary(forced_vs_cubic),
            "quartic_beats_cubic_zero_seed_count": int(
                sum(value < 0.0 for value in zero_vs_cubic)
            ),
            "quartic_beats_cubic_forced_seed_count": int(
                sum(value < 0.0 for value in forced_vs_cubic)
            ),
        },
        "formal_pde_gate_assessed": False,
        "reason_formal_pde_gate_unassessed": (
            "fresh fixed-probe pressure-free vorticity RMS is not the preregistered "
            "volume-weighted full-momentum L2/global maximum and pressure is not bound"
        ),
        "visualization_candidate_only": True,
        "canonical_velocity_changed": False,
        "visualization_ready": False,
        "visual_correspondence_verified": False,
        "pde_validated": False,
        "paper_exact": False,
        "openai_field_identified": False,
        "blowup_proved": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = audit_quartic_phi10_seed_generalization()
    text = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding="utf-8")
    print(text, end="")


if __name__ == "__main__":
    main()
