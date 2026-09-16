"""Independent PDE audit that consumes the frozen public velocity interface.

The existing constrained validator is independent of the optimizer, but candidate
validation historically called ``candidate.velocity`` directly.  This module
adds a user-facing check that obtains every velocity value through
``VelocityField.at_points`` (or an equivalent public field object), while
pressure and the preregistered force remain explicit inputs.

A passing finite-sample report is sampled evidence only.  It never promotes
``pde_validated`` and it does not establish visual correspondence or paper
identity.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Iterable, Mapping

import numpy as np

from .constrained_force import RestrictedForce
from .constrained_validation import residual, sampled_norms
from .velocity_components import DEFAULT_CANDIDATE, VelocityField


_REQUIRED_THRESHOLDS = (
    "pde_residual_max",
    "pde_residual_L2",
    "divergence_max",
    "divergence_L2",
)


def _validate_thresholds(thresholds: Mapping[str, float]) -> dict[str, float]:
    if not isinstance(thresholds, Mapping):
        raise TypeError("thresholds must be a mapping")
    result = {}
    for name in _REQUIRED_THRESHOLDS:
        if name not in thresholds:
            raise ValueError(f"missing threshold {name}")
        value = float(thresholds[name])
        if not np.isfinite(value) or value <= 0.0:
            raise ValueError(f"threshold {name} must be positive and finite")
        result[name] = value
    return result


def _validate_sampling(points, times, steps, volume):
    points = np.asarray(points, dtype=float)
    if points.ndim != 2 or points.shape[1] != 3 or points.shape[0] < 1:
        raise ValueError("points must have shape (N,3) with N>=1")
    if not np.all(np.isfinite(points)):
        raise ValueError("points must be finite")

    times = np.asarray(tuple(times), dtype=float)
    if times.ndim != 1 or times.size == 0 or not np.all(np.isfinite(times)):
        raise ValueError("times must be a nonempty finite 1D sequence")

    steps = np.asarray(tuple(steps), dtype=float)
    if steps.ndim != 1 or steps.size < 3 or not np.all(np.isfinite(steps)):
        raise ValueError("at least three finite derivative steps are required")
    if np.any(steps <= 0.0) or np.unique(steps).size != steps.size:
        raise ValueError("derivative steps must be positive and distinct")
    steps = np.sort(steps)[::-1]

    volume = float(volume)
    if not np.isfinite(volume) or volume <= 0.0:
        raise ValueError("volume must be positive and finite")
    return points, times, steps, volume


def _public_velocity(field, points, time):
    value = np.asarray(field.at_points(points, float(time)), dtype=float)
    points = np.asarray(points, dtype=float)
    if value.shape != points.shape or not np.all(np.isfinite(value)):
        raise ValueError("public velocity interface returned malformed/nonfinite data")
    return value


def _component_metrics(momentum):
    momentum = np.asarray(momentum, dtype=float)
    if momentum.ndim != 2 or momentum.shape[1] != 3 or not np.all(np.isfinite(momentum)):
        raise ValueError("momentum residual must have finite shape (N,3)")
    component_rms = np.sqrt(np.mean(momentum * momentum, axis=0))
    component_max = np.max(np.abs(momentum), axis=0)
    return {
        "momentum_component_rms": [float(v) for v in component_rms],
        "momentum_component_max_abs": [float(v) for v in component_max],
    }


def audit_public_velocity_refinement(
    field,
    pressure,
    force,
    points,
    times: Iterable[float],
    steps: Iterable[float],
    *,
    nu: float,
    time_bounds,
    volume: float,
    thresholds: Mapping[str, float],
):
    """Evaluate independent residuals using only the public velocity callable.

    ``pressure`` and ``force`` remain explicit callables.  The velocity values
    supplied to the independent finite-difference operator are obtained only
    from ``field.at_points``.  Ratios to the preregistered thresholds provide a
    normalized error report without redefining or relaxing those thresholds.
    """
    points, times, steps, volume = _validate_sampling(points, times, steps, volume)
    thresholds = _validate_thresholds(thresholds)
    nu = float(nu)
    time_bounds = tuple(float(v) for v in time_bounds)
    if len(time_bounds) != 2 or not np.all(np.isfinite(time_bounds)):
        raise ValueError("time_bounds must contain two finite values")
    if nu <= 0.0 or time_bounds[1] <= time_bounds[0]:
        raise ValueError("positive viscosity and increasing time bounds required")
    if np.any((times < time_bounds[0]) | (times > time_bounds[1])):
        raise ValueError("audit times must stay inside time_bounds")

    velocity = lambda p, t: _public_velocity(field, p, t)
    rows = []
    for time in times:
        for step in steps:
            result = residual(
                velocity,
                pressure,
                force,
                points,
                float(time),
                nu=nu,
                step=float(step),
                time_bounds=time_bounds,
            )
            norms = sampled_norms(result, volume)
            row = {
                "time": float(time),
                "step": float(step),
                **norms,
                **_component_metrics(result["momentum"]),
                "residual_max_over_threshold": float(
                    norms["residual_sampled_max"] / thresholds["pde_residual_max"]
                ),
                "residual_L2_over_threshold": float(
                    norms["residual_L2_estimate"] / thresholds["pde_residual_L2"]
                ),
                "divergence_max_over_threshold": float(
                    norms["divergence_sampled_max"] / thresholds["divergence_max"]
                ),
                "divergence_L2_over_threshold": float(
                    norms["divergence_L2_estimate"] / thresholds["divergence_L2"]
                ),
            }
            rows.append(row)

    finest = float(np.min(steps))
    finest_rows = [row for row in rows if row["step"] == finest]
    sampled_pass = all(
        row["residual_sampled_max"] <= thresholds["pde_residual_max"]
        and row["residual_L2_estimate"] <= thresholds["pde_residual_L2"]
        and row["divergence_sampled_max"] <= thresholds["divergence_max"]
        and row["divergence_L2_estimate"] <= thresholds["divergence_L2"]
        for row in finest_rows
    )
    candidate_sha256 = getattr(field, "sha256", None)
    if candidate_sha256 is not None and not isinstance(candidate_sha256, str):
        raise ValueError("field.sha256 must be a string when present")

    return {
        "candidate_sha256": candidate_sha256,
        "point_count": int(points.shape[0]),
        "times": [float(v) for v in times],
        "derivative_steps": [float(v) for v in steps],
        "finest_step": finest,
        "rows": rows,
        "sampled_finest_thresholds_passed": bool(sampled_pass),
        "velocity_source": "public field.at_points only",
        "pde_validated": False,
        "truth_boundary": (
            "finite held-out/public-interface PDE evidence only; a sampled pass does not by itself "
            "establish PDE validation, visual correspondence, paper identity, or blow-up"
        ),
    }


def audit_packaged_candidate_smoke(
    candidate_path=DEFAULT_CANDIDATE,
    *,
    config_path: str | Path | None = None,
    training_path: str | Path | None = None,
    seed: int = 914031,
    point_count: int = 16,
    times=(0.3125, 0.6875),
):
    """Small reproducible public-interface audit of a serialized candidate.

    This is intentionally a smoke-sized held-out check, not a replacement for
    the preregistered 4096-point validation artifact.
    """
    root = Path(__file__).resolve().parents[2]
    config_path = Path(config_path) if config_path is not None else root / "configs" / "constraints.json"
    training_path = (
        Path(training_path)
        if training_path is not None
        else root / "artifacts" / "constrained" / "coupled_joint" / "training.json"
    )
    config = json.loads(config_path.read_text(encoding="utf-8"))
    training = json.loads(training_path.read_text(encoding="utf-8"))
    if isinstance(seed, (bool, np.bool_)) or int(seed) != seed:
        raise TypeError("seed must be an integer")
    seed = int(seed)
    if seed == int(config["optimization"]["seed"]):
        raise ValueError("public audit seed must remain distinct from the training seed")
    if isinstance(point_count, (bool, np.bool_)) or int(point_count) != point_count:
        raise TypeError("point_count must be an integer")
    point_count = int(point_count)
    if point_count < 8:
        raise ValueError("point_count must be >=8")

    box = np.asarray(config["domain"]["evaluation_box"], dtype=float)
    rng = np.random.default_rng(seed)
    points = rng.uniform(box[:, 0], box[:, 1], size=(point_count, 3))
    field = VelocityField(candidate_path)
    force = RestrictedForce(**training["force"])
    report = audit_public_velocity_refinement(
        field,
        field.candidate.pressure,
        force,
        points,
        times,
        config["validation"]["derivative_steps"],
        nu=config["nu"],
        time_bounds=config["domain"]["time_interval"],
        volume=float(np.prod(box[:, 1] - box[:, 0])),
        thresholds=config["validation"]["thresholds"],
    )
    report.update(
        {
            "candidate_path": str(Path(candidate_path)),
            "config_path": str(config_path),
            "training_path": str(training_path),
            "seed": seed,
            "training_seed": int(config["optimization"]["seed"]),
            "validation_seed": int(config["validation"]["seed"]),
            "scope": "small held-out public-interface smoke; not the full preregistered validation sample",
        }
    )
    return report


def main():
    print(json.dumps(audit_packaged_candidate_smoke(), indent=2))


if __name__ == "__main__":
    main()
