"""Derivative contract and independent finite-difference oracle for CR003 candidates.

The production candidate may later provide analytic or automatic derivatives through
``velocity_derivatives``. This module deliberately keeps a separate numerical
reference implementation so derivative bugs can be detected without reusing the
candidate's own formulas.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol, runtime_checkable

import numpy as np


@dataclass(frozen=True)
class VelocityDerivativeBundle:
    """Velocity plus derivatives needed by an incompressible NS residual.

    Conventions:
      gradient[..., i, j] = d u_i / d x_j
      hessian[..., i, j, k] = d^2 u_i / (d x_j d x_k)
    """

    value: np.ndarray
    time: np.ndarray
    gradient: np.ndarray
    hessian: np.ndarray
    laplacian: np.ndarray

    @property
    def divergence(self) -> np.ndarray:
        return np.trace(self.gradient, axis1=-2, axis2=-1)


@runtime_checkable
class SupportsVelocity(Protocol):
    def velocity(self, points, time):
        ...


@runtime_checkable
class SupportsAnalyticVelocityDerivatives(SupportsVelocity, Protocol):
    def velocity_derivatives(self, points, time) -> VelocityDerivativeBundle:
        ...


def _positive_step(value: float, name: str) -> float:
    if isinstance(value, (bool, np.bool_)):
        raise TypeError(f"{name} must be a positive finite real number")
    try:
        value = float(value)
    except (TypeError, ValueError) as exc:
        raise TypeError(f"{name} must be a positive finite real number") from exc
    if not np.isfinite(value) or value <= 0:
        raise ValueError(f"{name} must be positive and finite")
    return value


def _prepare_points_time(points, time) -> tuple[np.ndarray, np.ndarray, tuple[int, ...]]:
    points = np.asarray(points, dtype=float)
    time = np.asarray(time, dtype=float)
    if points.ndim < 1 or points.shape[-1] != 3:
        raise ValueError("points must have shape (..., 3)")
    if not np.all(np.isfinite(points)) or not np.all(np.isfinite(time)):
        raise ValueError("points and time must be finite")
    try:
        shape = np.broadcast_shapes(points.shape[:-1], time.shape)
    except ValueError as exc:
        raise ValueError("time must broadcast against points[..., 0]") from exc
    return points, time, shape


def _eval_velocity(
    field: SupportsVelocity,
    points: np.ndarray,
    time: np.ndarray,
    expected_shape: tuple[int, ...],
) -> np.ndarray:
    value = np.asarray(field.velocity(points, time), dtype=float)
    expected = expected_shape + (3,)
    if value.shape != expected:
        raise ValueError(f"velocity returned shape {value.shape}, expected {expected}")
    if not np.all(np.isfinite(value)):
        raise ValueError("velocity returned non-finite values")
    return value


def reference_velocity_derivatives(
    field: SupportsVelocity,
    points,
    time,
    *,
    spatial_step: float = 1.0e-4,
    time_step: float = 1.0e-5,
) -> VelocityDerivativeBundle:
    """Compute an independent central-difference derivative bundle.

    This is a validation oracle, not a claim that finite differences are the final
    production derivative implementation. Callers must choose points/time far
    enough from any chart or time-domain boundary for the symmetric stencil.
    """

    if not hasattr(field, "velocity") or not callable(field.velocity):
        raise TypeError("field must provide a callable velocity(points, time)")
    h = _positive_step(spatial_step, "spatial_step")
    ht = _positive_step(time_step, "time_step")
    points, time, shape = _prepare_points_time(points, time)
    base = _eval_velocity(field, points, time, shape)

    gradient = np.empty(shape + (3, 3), dtype=float)
    hessian = np.empty(shape + (3, 3, 3), dtype=float)

    for j in range(3):
        delta = np.zeros_like(points, dtype=float)
        delta[..., j] = h
        plus = _eval_velocity(field, points + delta, time, shape)
        minus = _eval_velocity(field, points - delta, time, shape)
        gradient[..., :, j] = (plus - minus) / (2.0 * h)
        hessian[..., :, j, j] = (plus - 2.0 * base + minus) / (h * h)

    for j in range(3):
        for k in range(j + 1, 3):
            dj = np.zeros_like(points, dtype=float)
            dk = np.zeros_like(points, dtype=float)
            dj[..., j] = h
            dk[..., k] = h
            pp = _eval_velocity(field, points + dj + dk, time, shape)
            pm = _eval_velocity(field, points + dj - dk, time, shape)
            mp = _eval_velocity(field, points - dj + dk, time, shape)
            mm = _eval_velocity(field, points - dj - dk, time, shape)
            mixed = (pp - pm - mp + mm) / (4.0 * h * h)
            hessian[..., :, j, k] = mixed
            hessian[..., :, k, j] = mixed

    plus_t = _eval_velocity(field, points, time + ht, shape)
    minus_t = _eval_velocity(field, points, time - ht, shape)
    time_derivative = (plus_t - minus_t) / (2.0 * ht)
    laplacian = sum(hessian[..., :, j, j] for j in range(3))

    return VelocityDerivativeBundle(
        value=base,
        time=time_derivative,
        gradient=gradient,
        hessian=hessian,
        laplacian=laplacian,
    )
