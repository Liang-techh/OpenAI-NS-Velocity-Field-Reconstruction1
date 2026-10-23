"""Canonical exterior pressure gauge for the constrained candidate contract.

The CR001 pressure convention is ``p=0`` outside the cylindrical compact
support. In incompressible momentum equations, adding a spatially constant
function of time to pressure does not change ``grad(p)``. This module removes
only that gauge freedom. It deliberately refuses to repair a spatially
varying exterior pressure tail.
"""
from dataclasses import asdict, dataclass
from typing import Callable

import numpy as np

PressureCallable = Callable[[np.ndarray, np.ndarray], np.ndarray]


def _points(value, *, name):
    array = np.asarray(value, dtype=float)
    if array.ndim == 1:
        if array.shape != (3,):
            raise ValueError(f"{name} must have shape (3,) or (N,3)")
        array = array[None, :]
    if array.ndim != 2 or array.shape[1] != 3 or array.shape[0] == 0:
        raise ValueError(f"{name} must have shape (N,3) with N>=1")
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must be finite")
    return array


def _times(value, count):
    array = np.asarray(value, dtype=float)
    if array.ndim == 0:
        array = np.full(count, float(array))
    else:
        try:
            array = np.broadcast_to(array, (count,)).astype(float, copy=False)
        except ValueError as exc:
            raise ValueError("time must be scalar or broadcast to one value per point") from exc
    if not np.all(np.isfinite(array)):
        raise ValueError("time must be finite")
    return array


def _pressure_values(pressure, points, times):
    values = np.asarray(pressure(points, times), dtype=float)
    if values.shape != (len(points),):
        raise ValueError("pressure callable must return shape (N,)")
    if not np.all(np.isfinite(values)):
        raise ValueError("pressure callable returned nonfinite values")
    return values


@dataclass(frozen=True)
class ExteriorPressureGaugeReceipt:
    """Audit record for a canonical exterior pressure gauge."""

    reference_point: tuple[float, float, float]
    audit_times: tuple[float, ...]
    reference_offsets: tuple[float, ...]
    max_exterior_spatial_deviation: float
    max_abs_gauged_exterior: float
    probe_count: int
    tolerance: float
    radial_support: float
    axial_half_support: float
    gauge: str = "p=0 outside compact support"
    gradient_invariant_by_construction: bool = True
    velocity_changed: bool = False
    forcing_changed: bool = False
    pde_validated: bool = False

    def as_dict(self):
        return asdict(self)


@dataclass(frozen=True)
class ExteriorGaugedPressure:
    """Pressure wrapper subtracting ``p(reference_point,t)`` at every query."""

    base_pressure: PressureCallable
    reference_point: tuple[float, float, float]

    def __call__(self, points, time):
        points_array = _points(points, name="points")
        times = _times(time, len(points_array))
        reference = np.broadcast_to(
            np.asarray(self.reference_point, dtype=float), points_array.shape
        ).copy()
        values = _pressure_values(self.base_pressure, points_array, times)
        offsets = _pressure_values(self.base_pressure, reference, times)
        return values - offsets


def canonicalize_exterior_pressure_gauge(
    pressure: PressureCallable,
    exterior_probe_points,
    audit_times,
    *,
    radial_support: float = 2.0,
    axial_half_support: float = 2.0,
    tolerance: float = 1e-10,
):
    """Return a pressure with the CR001 exterior gauge and an audit receipt.

    The supplied probes must all lie outside the strict support
    ``r < radial_support`` and ``|z| < axial_half_support``. At every audit
    time the exterior pressure must be spatially constant to ``tolerance``;
    otherwise an additive pressure gauge cannot enforce the contract and this
    routine fails closed. The constant may depend on time, which is harmless
    to the spatial pressure gradient.
    """
    if not callable(pressure):
        raise ValueError("pressure must be callable")
    if not np.isfinite(radial_support) or radial_support <= 0:
        raise ValueError("radial_support must be finite and positive")
    if not np.isfinite(axial_half_support) or axial_half_support <= 0:
        raise ValueError("axial_half_support must be finite and positive")
    if not np.isfinite(tolerance) or tolerance < 0:
        raise ValueError("tolerance must be finite and nonnegative")

    probes = _points(exterior_probe_points, name="exterior_probe_points")
    radius = np.hypot(probes[:, 0], probes[:, 1])
    outside = (radius >= radial_support) | (np.abs(probes[:, 2]) >= axial_half_support)
    if not np.all(outside):
        raise ValueError("all pressure-gauge probes must lie outside compact support")

    times = np.asarray(audit_times, dtype=float)
    if times.ndim == 0:
        times = times.reshape(1)
    if times.ndim != 1 or len(times) == 0 or not np.all(np.isfinite(times)):
        raise ValueError("audit_times must be a nonempty finite one-dimensional array")

    reference_point = tuple(float(v) for v in probes[0])
    max_deviation = 0.0
    max_after = 0.0
    offsets = []
    for time in times:
        time_vector = np.full(len(probes), float(time))
        values = _pressure_values(pressure, probes, time_vector)
        offset = float(values[0])
        deviation = float(np.max(np.abs(values - offset)))
        if deviation > tolerance:
            raise ValueError(
                "exterior pressure is spatially nonconstant; an additive gauge "
                "cannot repair the boundary/support violation"
            )
        max_deviation = max(max_deviation, deviation)
        max_after = max(max_after, deviation)
        offsets.append(offset)

    gauged = ExteriorGaugedPressure(pressure, reference_point)
    receipt = ExteriorPressureGaugeReceipt(
        reference_point=reference_point,
        audit_times=tuple(float(v) for v in times),
        reference_offsets=tuple(offsets),
        max_exterior_spatial_deviation=max_deviation,
        max_abs_gauged_exterior=max_after,
        probe_count=len(probes),
        tolerance=float(tolerance),
        radial_support=float(radial_support),
        axial_half_support=float(axial_half_support),
    )
    return gauged, receipt
