"""Near-axis scaling diagnostic for callable Eq. (4.5) velocity candidates.

This module consumes only a public ``velocity(points, time) -> [...,3]`` callable.
It is intended to detect numerical/representation axis artifacts that could make
streamlines or vorticity tubes look grid-dependent.  It does not validate the
Navier--Stokes equations, physical support, visual correspondence, or blow-up.
"""
from __future__ import annotations

from typing import Callable, Iterable
import math

import numpy as np


Velocity = Callable[[np.ndarray, float], np.ndarray]


def _finite_tuple(name: str, values: Iterable[float]) -> tuple[float, ...]:
    result = tuple(float(value) for value in values)
    if not result or not all(math.isfinite(value) for value in result):
        raise ValueError(f"{name} must contain finite values")
    return result


def audit_eq45_near_axis_scaling(
    velocity: Velocity,
    *,
    times: Iterable[float] = (0.25, 0.5, 0.75),
    z_levels: Iterable[float] = (0.0, 0.15, -0.15),
    radii: Iterable[float] = (0.08, 0.04, 0.02, 0.01),
    azimuth_count: int = 8,
) -> dict:
    """Measure local axis scaling and azimuthal consistency.

    For a regular smooth axisymmetric field, transverse speed should approach
    ``O(r)`` and the axial component should differ from its axis value by
    ``O(r^2)``.  The returned observed orders are finite-radius diagnostics, not
    universal acceptance thresholds.
    """
    if not callable(velocity):
        raise TypeError("velocity must be callable")
    times_t = _finite_tuple("times", times)
    z_t = _finite_tuple("z_levels", z_levels)
    radii_t = _finite_tuple("radii", radii)
    if len(radii_t) < 3:
        raise ValueError("at least three radii are required")
    if any(radius <= 0.0 for radius in radii_t):
        raise ValueError("radii must be positive")
    if any(next_radius >= radius for radius, next_radius in zip(radii_t, radii_t[1:])):
        raise ValueError("radii must be strictly decreasing")
    if not isinstance(azimuth_count, (int, np.integer)) or azimuth_count < 4:
        raise ValueError("azimuth_count must be an integer >= 4")

    angles = 2.0 * np.pi * np.arange(azimuth_count, dtype=float) / azimuth_count
    rows = []
    all_axis_transverse = []
    all_transverse_spreads = []
    all_axial_spreads = []

    for time in times_t:
        for z in z_t:
            axis_point = np.array([[0.0, 0.0, z]], dtype=float)
            axis_velocity = np.asarray(velocity(axis_point, time), dtype=float)
            if axis_velocity.shape != (1, 3) or not np.all(np.isfinite(axis_velocity)):
                raise ValueError("velocity must return finite arrays matching points.shape")
            axis_velocity = axis_velocity[0]
            all_axis_transverse.append(float(np.linalg.norm(axis_velocity[:2])))

            transverse_means = []
            axial_departures = []
            max_transverse_spread = 0.0
            max_axial_spread = 0.0
            for radius in radii_t:
                points = np.stack(
                    (
                        radius * np.cos(angles),
                        radius * np.sin(angles),
                        np.full_like(angles, z),
                    ),
                    axis=-1,
                )
                values = np.asarray(velocity(points, time), dtype=float)
                if values.shape != points.shape or not np.all(np.isfinite(values)):
                    raise ValueError("velocity must return finite arrays matching points.shape")

                transverse = np.linalg.norm(values[:, :2], axis=1)
                transverse_over_r = transverse / radius
                axial = values[:, 2]
                transverse_means.append(float(np.mean(transverse)))
                axial_departures.append(float(abs(np.mean(axial) - axis_velocity[2])))
                max_transverse_spread = max(
                    max_transverse_spread, float(np.ptp(transverse_over_r))
                )
                max_axial_spread = max(max_axial_spread, float(np.ptp(axial)))

            transverse_means = np.asarray(transverse_means, dtype=float)
            axial_departures = np.asarray(axial_departures, dtype=float)
            if np.any(transverse_means <= 0.0):
                raise ValueError("transverse speed must be positive on sampled rings")
            if np.any(axial_departures <= 0.0):
                raise ValueError("axial departure must be positive on sampled rings")

            log_radius_ratio = np.log(
                np.asarray(radii_t[:-1], dtype=float) / np.asarray(radii_t[1:], dtype=float)
            )
            transverse_orders = np.log(transverse_means[:-1] / transverse_means[1:]) / log_radius_ratio
            axial_orders = np.log(axial_departures[:-1] / axial_departures[1:]) / log_radius_ratio

            rows.append(
                {
                    "time": time,
                    "z": z,
                    "axis_velocity": axis_velocity.tolist(),
                    "transverse_orders": transverse_orders.tolist(),
                    "axial_departure_orders": axial_orders.tolist(),
                    "max_transverse_over_r_azimuthal_spread": max_transverse_spread,
                    "max_axial_azimuthal_spread": max_axial_spread,
                }
            )
            all_transverse_spreads.append(max_transverse_spread)
            all_axial_spreads.append(max_axial_spread)

    finest_transverse = [row["transverse_orders"][-1] for row in rows]
    finest_axial = [row["axial_departure_orders"][-1] for row in rows]
    return {
        "sampling": {
            "times": list(times_t),
            "z_levels": list(z_t),
            "radii": list(radii_t),
            "azimuth_count": int(azimuth_count),
        },
        "rows": rows,
        "summary": {
            "min_finest_transverse_order": float(min(finest_transverse)),
            "max_finest_transverse_order": float(max(finest_transverse)),
            "min_finest_axial_departure_order": float(min(finest_axial)),
            "max_finest_axial_departure_order": float(max(finest_axial)),
            "max_transverse_over_r_azimuthal_spread": float(max(all_transverse_spreads)),
            "max_axial_azimuthal_spread": float(max(all_axial_spreads)),
            "max_axis_transverse_speed": float(max(all_axis_transverse)),
        },
    }
