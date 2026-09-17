"""Pressure-independent compatibility audit for a frozen velocity candidate.

For zero forcing, define the pressure-free momentum term

    M = u_t + (u . grad)u - nu * laplacian(u).

If a scalar pressure can close Navier--Stokes, then M = -grad(p), so curl(M)=0.
This module measures that necessary condition from the public velocity interface.
A small measured curl is not sufficient evidence for pressure existence or PDE validity.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any, Iterable

import numpy as np

Array = np.ndarray
VelocityLike = Any


@dataclass(frozen=True)
class PressureCompatibilityLevel:
    spatial_step: float
    time_step: float
    momentum_rms: float
    momentum_max: float
    curl_momentum_rms: float
    curl_momentum_max: float
    normalized_curl_rms: float


@dataclass(frozen=True)
class PressureCompatibilityReport:
    nu: float
    forcing_mode: str
    pressure_model: str
    characteristic_length: float
    point_count: int
    levels: tuple[PressureCompatibilityLevel, ...]
    velocity_changed: bool = False
    forcing_fitted: bool = False
    pressure_fitted: bool = False
    pde_validated: bool = False

    def to_dict(self) -> dict[str, Any]:
        payload = asdict(self)
        payload["levels"] = [asdict(level) for level in self.levels]
        payload["interpretation"] = (
            "curl(M)=0 is necessary for pressure-only zero-forcing closure; "
            "small numerical curl is not sufficient for a global pressure, boundary "
            "compatibility, or PDE validation"
        )
        return payload


def _as_points(points: Array) -> Array:
    arr = np.asarray(points, dtype=float)
    if arr.ndim != 2 or arr.shape[1] != 3 or arr.shape[0] == 0:
        raise ValueError("points must have shape (N, 3) with N > 0")
    if not np.all(np.isfinite(arr)):
        raise ValueError("points must be finite")
    return arr


def _as_times(times: float | Array, count: int) -> Array:
    arr = np.asarray(times, dtype=float)
    if arr.ndim == 0:
        arr = np.full(count, float(arr))
    else:
        arr = np.broadcast_to(arr, (count,)).astype(float, copy=False)
    if not np.all(np.isfinite(arr)):
        raise ValueError("times must be finite")
    return arr


def _evaluate_velocity(source: VelocityLike, points: Array, times: Array) -> Array:
    if hasattr(source, "velocity"):
        values = source.velocity(points, times)
    elif callable(source):
        values = source(points, times)
    else:
        raise TypeError("velocity source must be callable or expose velocity(points, time)")
    values = np.asarray(values, dtype=float)
    if values.shape != points.shape:
        raise ValueError(f"velocity output must have shape {points.shape}, got {values.shape}")
    if not np.all(np.isfinite(values)):
        raise ValueError("velocity output must be finite")
    return values


def pressure_free_momentum(
    velocity: VelocityLike,
    points: Array,
    times: float | Array,
    *,
    nu: float = 0.01,
    spatial_step: float = 0.01,
    time_step: float | None = None,
) -> Array:
    """Return M = u_t + (u.grad)u - nu*laplacian(u) by centered differences."""
    pts = _as_points(points)
    ts = _as_times(times, pts.shape[0])
    h = float(spatial_step)
    ht = h if time_step is None else float(time_step)
    viscosity = float(nu)
    if not np.isfinite(viscosity) or viscosity < 0.0:
        raise ValueError("nu must be finite and nonnegative")
    if not np.isfinite(h) or h <= 0.0 or not np.isfinite(ht) or ht <= 0.0:
        raise ValueError("finite-difference steps must be positive and finite")

    center = _evaluate_velocity(velocity, pts, ts)
    ddt = (
        _evaluate_velocity(velocity, pts, ts + ht)
        - _evaluate_velocity(velocity, pts, ts - ht)
    ) / (2.0 * ht)

    convective = np.zeros_like(center)
    laplacian = np.zeros_like(center)
    for axis in range(3):
        shift = np.zeros(3)
        shift[axis] = h
        plus = _evaluate_velocity(velocity, pts + shift, ts)
        minus = _evaluate_velocity(velocity, pts - shift, ts)
        derivative = (plus - minus) / (2.0 * h)
        convective += center[:, axis, None] * derivative
        laplacian += (plus - 2.0 * center + minus) / (h * h)

    return ddt + convective - viscosity * laplacian


def curl_pressure_free_momentum(
    velocity: VelocityLike,
    points: Array,
    times: float | Array,
    *,
    nu: float = 0.01,
    spatial_step: float = 0.01,
    time_step: float | None = None,
) -> tuple[Array, Array]:
    """Return (M, curl(M)) using an independent outer centered-difference curl."""
    pts = _as_points(points)
    ts = _as_times(times, pts.shape[0])
    h = float(spatial_step)
    center = pressure_free_momentum(
        velocity, pts, ts, nu=nu, spatial_step=h, time_step=time_step
    )

    derivatives: list[Array] = []
    for axis in range(3):
        shift = np.zeros(3)
        shift[axis] = h
        plus = pressure_free_momentum(
            velocity, pts + shift, ts, nu=nu, spatial_step=h, time_step=time_step
        )
        minus = pressure_free_momentum(
            velocity, pts - shift, ts, nu=nu, spatial_step=h, time_step=time_step
        )
        derivatives.append((plus - minus) / (2.0 * h))

    d_dx, d_dy, d_dz = derivatives
    curl = np.stack(
        (
            d_dy[:, 2] - d_dz[:, 1],
            d_dz[:, 0] - d_dx[:, 2],
            d_dx[:, 1] - d_dy[:, 0],
        ),
        axis=-1,
    )
    return center, curl


def audit_pressure_compatibility(
    velocity: VelocityLike,
    points: Array,
    times: float | Array,
    *,
    nu: float = 0.01,
    derivative_steps: Iterable[float] = (0.02, 0.01, 0.005),
    characteristic_length: float = 1.0,
) -> PressureCompatibilityReport:
    """Audit a frozen public velocity under the zero-forcing pressure necessity."""
    pts = _as_points(points)
    ts = _as_times(times, pts.shape[0])
    length = float(characteristic_length)
    if not np.isfinite(length) or length <= 0.0:
        raise ValueError("characteristic_length must be positive and finite")

    steps = tuple(float(step) for step in derivative_steps)
    if len(steps) < 2 or any((not np.isfinite(step) or step <= 0.0) for step in steps):
        raise ValueError("derivative_steps must contain at least two positive finite levels")
    if any(b >= a for a, b in zip(steps, steps[1:])):
        raise ValueError("derivative_steps must be strictly decreasing")

    levels: list[PressureCompatibilityLevel] = []
    for step in steps:
        momentum, curl_momentum = curl_pressure_free_momentum(
            velocity,
            pts,
            ts,
            nu=nu,
            spatial_step=step,
            time_step=step,
        )
        momentum_norm = np.linalg.norm(momentum, axis=-1)
        curl_norm = np.linalg.norm(curl_momentum, axis=-1)
        momentum_rms = float(np.sqrt(np.mean(momentum_norm * momentum_norm)))
        curl_rms = float(np.sqrt(np.mean(curl_norm * curl_norm)))
        levels.append(
            PressureCompatibilityLevel(
                spatial_step=step,
                time_step=step,
                momentum_rms=momentum_rms,
                momentum_max=float(np.max(momentum_norm)),
                curl_momentum_rms=curl_rms,
                curl_momentum_max=float(np.max(curl_norm)),
                normalized_curl_rms=float(length * curl_rms / max(momentum_rms, 1e-15)),
            )
        )

    return PressureCompatibilityReport(
        nu=float(nu),
        forcing_mode="zero_forcing_subcase",
        pressure_model="scalar_pressure_gradient_only",
        characteristic_length=length,
        point_count=pts.shape[0],
        levels=tuple(levels),
    )
