"""Time-step-only PDE refinement through the public velocity API.

This diagnostic keeps the fourth-order Cartesian spatial stencil fixed while
varying only the second-order time-difference step.  It is deliberately
separate from training loss and from the existing validator whose single
``step`` changes spatial and temporal differentiation together.

The report is finite-sample convergence evidence only.  It never promotes a
field to PDE-valid, visually matched, paper-exact, or singular.
"""
from __future__ import annotations

from collections.abc import Iterable, Mapping

import numpy as np


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


def _validate_common(points, time, nu, spatial_step, time_bounds):
    points = np.asarray(points, dtype=float)
    if points.ndim != 2 or points.shape[1] != 3 or points.shape[0] == 0:
        raise ValueError("points must have shape (N,3) with N>=1")
    if not np.all(np.isfinite(points)):
        raise ValueError("points must be finite")
    time = float(time)
    nu = float(nu)
    spatial_step = float(spatial_step)
    bounds = tuple(float(value) for value in time_bounds)
    if len(bounds) != 2 or not np.all(np.isfinite(bounds)) or bounds[1] <= bounds[0]:
        raise ValueError("time_bounds must be finite and increasing")
    if not bounds[0] <= time <= bounds[1]:
        raise ValueError("time must lie inside time_bounds")
    if not np.isfinite(nu) or nu <= 0.0:
        raise ValueError("nu must be positive and finite")
    if not np.isfinite(spatial_step) or spatial_step <= 0.0:
        raise ValueError("spatial_step must be positive and finite")
    return points, time, nu, spatial_step, bounds


def _fixed_spatial_terms(field, pressure, force, points, time, *, nu, spatial_step):
    velocity = lambda p, t: _public_velocity(field, p, t)
    u = velocity(points, time)
    jac = np.empty((len(points), 3, 3), dtype=float)
    lap = np.zeros_like(u)
    gradp = np.empty_like(u)

    for axis in range(3):
        shift = np.zeros(3, dtype=float)
        shift[axis] = spatial_step
        um2, um1, up1, up2 = [
            velocity(points + multiplier * shift, time)
            for multiplier in (-2.0, -1.0, 1.0, 2.0)
        ]
        jac[:, :, axis] = (um2 - 8.0 * um1 + 8.0 * up1 - up2) / (12.0 * spatial_step)
        lap += (-up2 + 16.0 * up1 - 30.0 * u + 16.0 * um1 - um2) / (
            12.0 * spatial_step * spatial_step
        )
        pm2, pm1, pp1, pp2 = [
            _scalar_call("pressure", pressure, points + multiplier * shift, time)
            for multiplier in (-2.0, -1.0, 1.0, 2.0)
        ]
        gradp[:, axis] = (pm2 - 8.0 * pm1 + 8.0 * pp1 - pp2) / (12.0 * spatial_step)

    forcing = _vector_call("force", force, points, time)
    base = np.einsum("nij,nj->ni", jac, u) + gradp - nu * lap - forcing
    divergence = np.trace(jac, axis1=1, axis2=2)
    if not np.all(np.isfinite(base)) or not np.all(np.isfinite(divergence)):
        raise ValueError("fixed spatial operator produced nonfinite values")
    return u, base, divergence


def _time_derivative(field, points, time, step, bounds, u=None):
    step = float(step)
    if not np.isfinite(step) or step <= 0.0:
        raise ValueError("time step must be positive and finite")
    lo, hi = bounds
    if 2.0 * step > hi - lo:
        raise ValueError("time step is too large for the declared interval")
    velocity = lambda p, t: _public_velocity(field, p, t)
    if u is None:
        u = velocity(points, time)
    if time - step < lo:
        if time + 2.0 * step > hi:
            raise ValueError("forward time stencil leaves declared interval")
        return (-3.0 * u + 4.0 * velocity(points, time + step) - velocity(points, time + 2.0 * step)) / (
            2.0 * step
        )
    if time + step > hi:
        if time - 2.0 * step < lo:
            raise ValueError("backward time stencil leaves declared interval")
        return (3.0 * u - 4.0 * velocity(points, time - step) + velocity(points, time - 2.0 * step)) / (
            2.0 * step
        )
    return (velocity(points, time + step) - velocity(points, time - step)) / (2.0 * step)


def temporal_refinement_residual(
    field,
    pressure,
    force,
    points,
    time,
    *,
    nu=0.01,
    spatial_step=0.005,
    time_step=0.005,
    time_bounds=(0.25, 0.75),
):
    """Evaluate residual while keeping spatial and temporal steps independent."""
    points, time, nu, spatial_step, bounds = _validate_common(
        points, time, nu, spatial_step, time_bounds
    )
    u, base, divergence = _fixed_spatial_terms(
        field,
        pressure,
        force,
        points,
        time,
        nu=nu,
        spatial_step=spatial_step,
    )
    ut = _time_derivative(field, points, time, time_step, bounds, u=u)
    momentum = ut + base
    if not np.all(np.isfinite(momentum)):
        raise ValueError("temporal refinement residual produced nonfinite values")
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
    magnitude = np.linalg.norm(momentum, axis=1)
    return {
        "residual_sampled_max": float(np.max(magnitude)),
        "residual_L2_estimate": float(np.sqrt(volume * np.mean(magnitude * magnitude))),
        "divergence_sampled_max": float(np.max(np.abs(divergence))),
        "divergence_L2_estimate": float(np.sqrt(volume * np.mean(divergence * divergence))),
        "component_rms": [float(value) for value in np.sqrt(np.mean(momentum * momentum, axis=0))],
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


def audit_public_velocity_time_refinement(
    field,
    pressure,
    force,
    points,
    times: Iterable[float],
    time_steps: Iterable[float],
    *,
    nu,
    spatial_step,
    time_bounds,
    volume,
    thresholds,
):
    """Vary only the time-difference step on one frozen public velocity field."""
    points = np.asarray(points, dtype=float)
    if points.ndim != 2 or points.shape[1] != 3 or points.shape[0] == 0 or not np.all(np.isfinite(points)):
        raise ValueError("points must be a finite (N,3) array with N>=1")
    times = np.asarray(tuple(times), dtype=float)
    steps = np.asarray(tuple(time_steps), dtype=float)
    if times.ndim != 1 or times.size == 0 or not np.all(np.isfinite(times)):
        raise ValueError("times must be a nonempty finite sequence")
    if steps.ndim != 1 or steps.size < 3 or not np.all(np.isfinite(steps)):
        raise ValueError("at least three finite time steps are required")
    if np.any(steps <= 0.0) or np.unique(steps).size != steps.size:
        raise ValueError("time steps must be positive and distinct")
    steps = np.sort(steps)[::-1]
    thresholds = _thresholds(thresholds)
    bounds = tuple(float(value) for value in time_bounds)
    if len(bounds) != 2 or not np.all(np.isfinite(bounds)) or bounds[1] <= bounds[0]:
        raise ValueError("time_bounds must be finite and increasing")
    if np.any((times < bounds[0]) | (times > bounds[1])):
        raise ValueError("audit times must stay inside time_bounds")
    nu = float(nu)
    spatial_step = float(spatial_step)
    volume = float(volume)
    if not np.isfinite(nu) or nu <= 0.0:
        raise ValueError("nu must be positive and finite")
    if not np.isfinite(spatial_step) or spatial_step <= 0.0:
        raise ValueError("spatial_step must be positive and finite")
    if not np.isfinite(volume) or volume <= 0.0:
        raise ValueError("volume must be positive and finite")

    rows = []
    by_time = []
    for time in times:
        u, base, divergence = _fixed_spatial_terms(
            field,
            pressure,
            force,
            points,
            float(time),
            nu=nu,
            spatial_step=spatial_step,
        )
        momenta = []
        time_rows = []
        for step in steps:
            ut = _time_derivative(field, points, float(time), float(step), bounds, u=u)
            result = {"momentum": ut + base, "divergence": divergence}
            norms = _norms(result, volume)
            row = {
                "time": float(time),
                "time_step": float(step),
                **norms,
                "residual_max_over_threshold": float(norms["residual_sampled_max"] / thresholds["pde_residual_max"]),
                "residual_L2_over_threshold": float(norms["residual_L2_estimate"] / thresholds["pde_residual_L2"]),
                "divergence_max_over_threshold": float(norms["divergence_sampled_max"] / thresholds["divergence_max"]),
                "divergence_L2_over_threshold": float(norms["divergence_L2_estimate"] / thresholds["divergence_L2"]),
            }
            rows.append(row)
            time_rows.append(row)
            momenta.append(result["momentum"])

        adjacent = []
        for index in range(len(steps) - 1):
            delta = momenta[index] - momenta[index + 1]
            magnitude = np.linalg.norm(delta, axis=1)
            adjacent.append(
                {
                    "coarse_step": float(steps[index]),
                    "fine_step": float(steps[index + 1]),
                    "momentum_delta_rms": float(np.sqrt(np.mean(magnitude * magnitude))),
                    "momentum_delta_max": float(np.max(magnitude)),
                }
            )
        by_time.append(
            {
                "time": float(time),
                "rows": time_rows,
                "adjacent_time_step_deltas": adjacent,
                "finest_to_coarsest_residual_max_ratio": float(
                    time_rows[-1]["residual_sampled_max"] / time_rows[0]["residual_sampled_max"]
                ) if time_rows[0]["residual_sampled_max"] > 0.0 else None,
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
        "fixed_spatial_step": spatial_step,
        "time_steps": [float(value) for value in steps],
        "rows": rows,
        "time_refinement": by_time,
        "thresholds_unchanged": thresholds,
        "pde_validated": False,
        "visualization_candidate_only_allowed": True,
        "truth_boundary": (
            "time-step-only finite-sample convergence evidence; a stable visual field may still be a "
            "visualization candidate only, and neither refinement nor threshold crossing establishes "
            "paper identity, full PDE validation, singularity, or blow-up"
        ),
    }
