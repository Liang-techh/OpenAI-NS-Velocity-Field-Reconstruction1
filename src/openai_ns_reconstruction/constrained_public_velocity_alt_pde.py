"""Alternate finite-difference PDE cross-check through the public velocity API.

This module intentionally implements a second-order Cartesian spatial operator
separately from ``constrained_validation``, whose production validator uses
fourth-order spatial stencils. The two operators can be compared on the same
frozen public velocity artifact. Pressure and forcing remain explicit inputs.

The report is numerical cross-check evidence only. It never promotes a field
to ``pde_validated`` and does not establish visual correspondence, paper
identity, or blow-up.
"""
from __future__ import annotations

from collections.abc import Iterable, Mapping

import numpy as np

from .constrained_validation import residual as fourth_order_residual


_REQUIRED_THRESHOLDS = (
    "pde_residual_max",
    "pde_residual_L2",
    "divergence_max",
    "divergence_L2",
)


def _public_velocity(field, points, time):
    points = np.asarray(points, dtype=float)
    value = np.asarray(field.at_points(points, float(time)), dtype=float)
    if value.shape != points.shape or not np.all(np.isfinite(value)):
        raise ValueError("public velocity interface returned malformed/nonfinite data")
    return value


def _vector_call(name, fn, points, time):
    value = np.asarray(fn(points, float(time)), dtype=float)
    if value.shape != points.shape or not np.all(np.isfinite(value)):
        raise ValueError(f"{name} must return finite shape (N,3)")
    return value


def _scalar_call(name, fn, points, time):
    value = np.asarray(fn(points, float(time)), dtype=float)
    if value.shape == (len(points), 1):
        value = value[:, 0]
    if value.shape != (len(points),) or not np.all(np.isfinite(value)):
        raise ValueError(f"{name} must return finite shape (N,)")
    return value


def second_order_public_residual(
    field,
    pressure,
    force,
    points,
    time,
    *,
    nu=0.01,
    step=0.005,
    time_bounds=(0.25, 0.75),
):
    """Evaluate NS residual with an independently coded second-order stencil."""
    points = np.asarray(points, dtype=float)
    if points.ndim != 2 or points.shape[1] != 3 or points.shape[0] == 0:
        raise ValueError("points must have finite shape (N,3) with N>=1")
    if not np.all(np.isfinite(points)):
        raise ValueError("points must have finite shape (N,3) with N>=1")
    time = float(time)
    nu = float(nu)
    step = float(step)
    bounds = tuple(float(value) for value in time_bounds)
    if len(bounds) != 2 or not np.all(np.isfinite(bounds)):
        raise ValueError("time_bounds must contain two finite values")
    lo, hi = bounds
    if not np.isfinite([time, nu, step]).all() or nu <= 0.0 or step <= 0.0:
        raise ValueError("time, viscosity and step must be finite with nu/step positive")
    if hi <= lo or not lo <= time <= hi or 2.0 * step > hi - lo:
        raise ValueError("invalid time domain or time stencil width")

    velocity = lambda p, t: _public_velocity(field, p, t)
    u = velocity(points, time)
    jac = np.empty((len(points), 3, 3), dtype=float)
    lap = np.zeros_like(u)
    gradp = np.empty_like(u)

    for axis in range(3):
        shift = np.zeros(3, dtype=float)
        shift[axis] = step
        minus = points - shift
        plus = points + shift
        um = velocity(minus, time)
        up = velocity(plus, time)
        jac[:, :, axis] = (up - um) / (2.0 * step)
        lap += (up - 2.0 * u + um) / (step * step)
        pm = _scalar_call("pressure", pressure, minus, time)
        pp = _scalar_call("pressure", pressure, plus, time)
        gradp[:, axis] = (pp - pm) / (2.0 * step)

    if time - step < lo:
        ut = (
            -3.0 * u
            + 4.0 * velocity(points, time + step)
            - velocity(points, time + 2.0 * step)
        ) / (2.0 * step)
    elif time + step > hi:
        ut = (
            3.0 * u
            - 4.0 * velocity(points, time - step)
            + velocity(points, time - 2.0 * step)
        ) / (2.0 * step)
    else:
        ut = (velocity(points, time + step) - velocity(points, time - step)) / (2.0 * step)

    forcing = _vector_call("force", force, points, time)
    momentum = ut + np.einsum("nij,nj->ni", jac, u) + gradp - nu * lap - forcing
    divergence = np.trace(jac, axis1=1, axis2=2)
    if not np.all(np.isfinite(momentum)) or not np.all(np.isfinite(divergence)):
        raise ValueError("alternate residual produced nonfinite values")
    return {"momentum": momentum, "divergence": divergence}


def _norms(result, volume):
    volume = float(volume)
    if not np.isfinite(volume) or volume <= 0.0:
        raise ValueError("volume must be positive and finite")
    momentum = np.asarray(result["momentum"], dtype=float)
    divergence = np.asarray(result["divergence"], dtype=float)
    if momentum.ndim != 2 or momentum.shape[1] != 3 or divergence.shape != (len(momentum),):
        raise ValueError("malformed residual arrays")
    if not np.all(np.isfinite(momentum)) or not np.all(np.isfinite(divergence)):
        raise ValueError("residual arrays must be finite")
    mag = np.linalg.norm(momentum, axis=1)
    div = np.abs(divergence)
    return {
        "residual_sampled_max": float(np.max(mag)),
        "residual_L2_estimate": float(np.sqrt(volume * np.mean(mag * mag))),
        "divergence_sampled_max": float(np.max(div)),
        "divergence_L2_estimate": float(np.sqrt(volume * np.mean(div * div))),
    }


def _thresholds(values: Mapping[str, float]) -> dict[str, float]:
    if not isinstance(values, Mapping):
        raise TypeError("thresholds must be a mapping")
    result = {}
    for name in _REQUIRED_THRESHOLDS:
        if name not in values:
            raise ValueError(f"missing threshold {name}")
        value = float(values[name])
        if not np.isfinite(value) or value <= 0.0:
            raise ValueError(f"threshold {name} must be positive and finite")
        result[name] = value
    return result


def _sampling(points, times: Iterable[float], steps: Iterable[float]):
    points = np.asarray(points, dtype=float)
    if points.ndim != 2 or points.shape[1] != 3 or points.shape[0] == 0:
        raise ValueError("points must have shape (N,3) with N>=1")
    if not np.all(np.isfinite(points)):
        raise ValueError("points must be finite")
    times = np.asarray(tuple(times), dtype=float)
    steps = np.asarray(tuple(steps), dtype=float)
    if times.ndim != 1 or times.size == 0 or not np.all(np.isfinite(times)):
        raise ValueError("times must be a nonempty finite sequence")
    if steps.ndim != 1 or steps.size < 3 or not np.all(np.isfinite(steps)):
        raise ValueError("at least three finite derivative steps are required")
    if np.any(steps <= 0.0) or np.unique(steps).size != steps.size:
        raise ValueError("derivative steps must be positive and distinct")
    return points, times, np.sort(steps)[::-1]


def audit_public_velocity_operator_crosscheck(
    field,
    pressure,
    force,
    points,
    times,
    steps,
    *,
    nu,
    time_bounds,
    volume,
    thresholds,
):
    """Compare independent second- and fourth-order operators on one public field."""
    points, times, steps = _sampling(points, times, steps)
    thresholds = _thresholds(thresholds)
    bounds = tuple(float(value) for value in time_bounds)
    if len(bounds) != 2 or not np.all(np.isfinite(bounds)) or bounds[1] <= bounds[0]:
        raise ValueError("time_bounds must be finite and increasing")
    if np.any((times < bounds[0]) | (times > bounds[1])):
        raise ValueError("audit times must stay inside time_bounds")
    nu = float(nu)
    volume = float(volume)
    if not np.isfinite(nu) or nu <= 0.0:
        raise ValueError("nu must be positive and finite")
    if not np.isfinite(volume) or volume <= 0.0:
        raise ValueError("volume must be positive and finite")

    velocity = lambda p, t: _public_velocity(field, p, t)
    rows = []
    for time in times:
        for step in steps:
            alternate = second_order_public_residual(
                field,
                pressure,
                force,
                points,
                float(time),
                nu=nu,
                step=float(step),
                time_bounds=bounds,
            )
            reference = fourth_order_residual(
                velocity,
                pressure,
                force,
                points,
                float(time),
                nu=nu,
                step=float(step),
                time_bounds=bounds,
            )
            alt_norms = _norms(alternate, volume)
            ref_norms = _norms(reference, volume)
            momentum_delta = alternate["momentum"] - reference["momentum"]
            divergence_delta = alternate["divergence"] - reference["divergence"]
            momentum_delta_mag = np.linalg.norm(momentum_delta, axis=1)
            rows.append(
                {
                    "time": float(time),
                    "step": float(step),
                    "second_order": alt_norms,
                    "fourth_order": ref_norms,
                    "operator_delta_residual_rms": float(np.sqrt(np.mean(momentum_delta_mag**2))),
                    "operator_delta_residual_max": float(np.max(momentum_delta_mag)),
                    "operator_delta_divergence_rms": float(np.sqrt(np.mean(divergence_delta**2))),
                    "operator_delta_divergence_max": float(np.max(np.abs(divergence_delta))),
                    "second_order_residual_max_over_threshold": float(
                        alt_norms["residual_sampled_max"] / thresholds["pde_residual_max"]
                    ),
                    "second_order_residual_L2_over_threshold": float(
                        alt_norms["residual_L2_estimate"] / thresholds["pde_residual_L2"]
                    ),
                    "second_order_divergence_max_over_threshold": float(
                        alt_norms["divergence_sampled_max"] / thresholds["divergence_max"]
                    ),
                    "second_order_divergence_L2_over_threshold": float(
                        alt_norms["divergence_L2_estimate"] / thresholds["divergence_L2"]
                    ),
                    "fourth_order_residual_max_over_threshold": float(
                        ref_norms["residual_sampled_max"] / thresholds["pde_residual_max"]
                    ),
                    "fourth_order_divergence_max_over_threshold": float(
                        ref_norms["divergence_sampled_max"] / thresholds["divergence_max"]
                    ),
                }
            )

    candidate_sha256 = getattr(field, "sha256", None)
    if candidate_sha256 is not None and not isinstance(candidate_sha256, str):
        raise ValueError("field.sha256 must be a string when present")
    return {
        "candidate_sha256": candidate_sha256,
        "velocity_source": "public field.at_points only",
        "point_count": int(len(points)),
        "times": [float(value) for value in times],
        "derivative_steps": [float(value) for value in steps],
        "finest_step": float(np.min(steps)),
        "rows": rows,
        "operators": {
            "alternate": "independently coded second-order Cartesian spatial finite differences",
            "reference": "existing optimizer-independent fourth-order Cartesian spatial finite differences",
        },
        "thresholds_unchanged": thresholds,
        "pde_validated": False,
        "visualization_candidate_only_allowed": True,
        "truth_boundary": (
            "cross-operator finite-sample evidence only; operator agreement or a sampled threshold pass "
            "does not establish PDE validation, visual correspondence, paper identity, or blow-up"
        ),
    }
