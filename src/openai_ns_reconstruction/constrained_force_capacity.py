"""Capacity diagnostic for the preregistered two-parameter forcing family.

This module never constructs a residual-defined force.  It asks only how much of
a *fixed* momentum residual lies in the span of the already-preregistered
RestrictedForce(a, c) family with its original coefficient bounds.

A fit performed on the same residual samples is a capacity ceiling and is not
independent validation or evidence that the PDE has passed.
"""
from __future__ import annotations

from typing import Any

import numpy as np
from scipy.optimize import lsq_linear

from .constrained_force import RestrictedForce
from .constrained_validation import residual


def _finite_points(points: Any) -> np.ndarray:
    x = np.asarray(points, dtype=float)
    if x.ndim != 2 or x.shape[1] != 3 or not np.all(np.isfinite(x)):
        raise ValueError("points must be a finite (N,3) array")
    if len(x) == 0:
        raise ValueError("points must be nonempty")
    return x


def _finite_residual(base_residual: Any, n: int) -> np.ndarray:
    r = np.asarray(base_residual, dtype=float)
    if r.shape != (n, 3) or not np.all(np.isfinite(r)):
        raise ValueError("base_residual must be a finite (N,3) array")
    return r


def fit_restricted_force_capacity(
    base_residual: Any,
    points: Any,
    time: float,
    *,
    lower: float = 0.0,
    upper: float = 10.0,
) -> dict[str, Any]:
    """Best bounded residual reduction available to the frozen force family.

    ``base_residual`` is the momentum residual with ``a=c=0`` while velocity
    and pressure are held fixed.  Because RestrictedForce is exactly linear in
    ``a`` and ``c``, the bounded least-squares solve below is the exact capacity
    ceiling of that predeclared two-dimensional forcing family on these samples.

    The returned ``after_*`` metrics are deliberately labelled capacity-only:
    the same samples are used to fit and score the two force coefficients.
    """
    x = _finite_points(points)
    r0 = _finite_residual(base_residual, len(x))
    if not np.isfinite(time):
        raise ValueError("time must be finite")
    if not np.isfinite([lower, upper]).all() or lower < 0 or upper <= lower:
        raise ValueError("forcing bounds must be finite with 0 <= lower < upper")
    if lower < 0 or upper > 10:
        raise ValueError("diagnostic bounds may not exceed the preregistered [0,10] force bounds")

    fa = RestrictedForce(1.0, 0.0)(x, time)
    fc = RestrictedForce(0.0, 1.0)(x, time)
    matrix = np.column_stack((fa.ravel(), fc.ravel()))
    singular_values = np.linalg.svd(matrix, compute_uv=False)
    if singular_values.size:
        tol = max(matrix.shape) * np.finfo(float).eps * singular_values[0]
        rank = int(np.sum(singular_values > tol))
    else:
        rank = 0
    condition = (
        float(singular_values[0] / singular_values[rank - 1])
        if rank == 2 and singular_values[rank - 1] > 0
        else None
    )

    fit = lsq_linear(
        matrix,
        r0.ravel(),
        bounds=([lower, lower], [upper, upper]),
        tol=1e-12,
        max_iter=200,
    )
    if not fit.success or not np.all(np.isfinite(fit.x)):
        raise RuntimeError("bounded restricted-force capacity solve failed: " + fit.message)

    fitted_force = (matrix @ fit.x).reshape(r0.shape)
    remaining = r0 - fitted_force
    before_vec = np.linalg.norm(r0, axis=1)
    after_vec = np.linalg.norm(remaining, axis=1)
    before_rms = float(np.sqrt(np.mean(np.sum(r0 * r0, axis=1))))
    after_rms = float(np.sqrt(np.mean(np.sum(remaining * remaining, axis=1))))
    ratio_floor = 100 * np.finfo(float).eps
    recoverable_fraction = (
        float(max(0.0, min(1.0, 1.0 - after_rms / before_rms)))
        if before_rms > ratio_floor
        else None
    )
    bound_tol = 1e-8 * max(1.0, abs(lower), abs(upper))
    labels = ("a", "c")
    active_bounds: list[str] = []
    for label, value in zip(labels, fit.x):
        if abs(value - lower) <= bound_tol:
            active_bounds.append(f"{label}:lower")
        elif abs(value - upper) <= bound_tol:
            active_bounds.append(f"{label}:upper")

    return {
        "scope": "same-sample restricted-force capacity ceiling; not independent validation",
        "forcing_family": "RestrictedForce(a,c) preregistered before optimization",
        "coefficient_bounds": {"a": [float(lower), float(upper)], "c": [float(lower), float(upper)]},
        "coefficients": {"a": float(fit.x[0]), "c": float(fit.x[1])},
        "active_bounds": active_bounds,
        "design_rank": rank,
        "design_condition_number": condition,
        "design_singular_values": [float(v) for v in singular_values],
        "before_rms": before_rms,
        "before_max": float(before_vec.max()),
        "capacity_rms": after_rms,
        "capacity_max": float(after_vec.max()),
        "recoverable_fraction": recoverable_fraction,
        "point_count": int(len(x)),
        "time": float(time),
        "pde_validated": False,
    }


def diagnose_fixed_candidate_force_capacity(
    velocity,
    pressure,
    points: Any,
    time: float,
    *,
    nu: float = 0.01,
    step: float = 0.005,
    time_bounds: tuple[float, float] = (0.25, 0.75),
) -> dict[str, Any]:
    """Compute the unforced independent residual, then its restricted-force ceiling.

    Velocity and pressure callables are never mutated or fitted.  The force fit
    sees the diagnostic samples, so this function must not be used to claim an
    independently validated Navier--Stokes residual.
    """
    x = _finite_points(points)
    base = residual(
        velocity,
        pressure,
        RestrictedForce(0.0, 0.0),
        x,
        time,
        nu=nu,
        step=step,
        time_bounds=time_bounds,
    )["momentum"]
    result = fit_restricted_force_capacity(base, x, time)
    result.update(
        {
            "operator": "independent constrained_validation.residual with fourth-order spatial / second-order time stencils",
            "velocity_and_pressure_frozen": True,
            "force_fit_uses_diagnostic_samples": True,
            "nu": float(nu),
            "step": float(step),
        }
    )
    return result
