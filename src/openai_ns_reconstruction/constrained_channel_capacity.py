"""Small residual-channel probes for the current axisymmetric candidate family.

The purpose is diagnostic: distinguish directions that can affect the dominant
axial residual from directions that are structurally blind to it. These probes
are not production candidate parameters and never modify pressure or forcing.
"""
from __future__ import annotations

import numpy as np

from .cutoffs import standard_cutoff, standard_cutoff_derivative

PEAK_RADIUS = 1.0153968253968253
PEAK_Z = 0.6879365079365078


def _points(points):
    value = np.asarray(points, dtype=float)
    if value.ndim < 1 or value.shape[-1] != 3 or not np.all(np.isfinite(value)):
        raise ValueError("points must be finite with trailing dimension 3")
    return value


def _time_factor(time, shape):
    value = np.asarray(time, dtype=float)
    if not np.all(np.isfinite(value)):
        raise ValueError("time must be finite")
    try:
        return np.broadcast_to(2.0 * (value - 0.25), shape)
    except ValueError as exc:
        raise ValueError("time is not broadcast-compatible with points") from exc


def _cutoff(values):
    values = np.asarray(values, dtype=float)
    flat = np.fromiter((standard_cutoff(float(v)) for v in values.ravel()), dtype=float, count=values.size)
    return flat.reshape(values.shape)


def _cutoff_derivative(values):
    values = np.asarray(values, dtype=float)
    flat = np.fromiter((standard_cutoff_derivative(float(v)) for v in values.ravel()), dtype=float, count=values.size)
    return flat.reshape(values.shape)


def pure_swirl_probe(points, time):
    """Compact axisymmetric pure-swirl mutation centered near the axial peak."""
    value = _points(points)
    x, y, z = np.moveaxis(value, -1, 0)
    r2 = x * x + y * y
    envelope = _cutoff(r2 / 4.0) * _cutoff(z * z / 4.0)
    local = np.exp(-((r2 - PEAK_RADIUS**2) / 0.35) ** 2 - ((z - PEAK_Z) / 0.35) ** 2)
    factor = _time_factor(time, r2.shape) * envelope * local
    return np.stack((-y * factor, x * factor, np.zeros_like(factor)), axis=-1)


def poloidal_curl_probe(points, time):
    """Compact divergence-free poloidal mutation from ``curl((-y q,x q,0))``."""
    value = _points(points)
    x, y, z = np.moveaxis(value, -1, 0)
    r2 = x * x + y * y
    t = _time_factor(time, r2.shape)
    sr = r2 / 4.0
    sz = z * z / 4.0
    cr = _cutoff(sr)
    cz = _cutoff(sz)
    dcr_dr2 = _cutoff_derivative(sr) / 4.0
    dcz_dz = _cutoff_derivative(sz) * z / 2.0
    dr = r2 - PEAK_RADIUS**2
    dz = z - PEAK_Z
    gaussian = np.exp(-(dr / 0.35) ** 2 - (dz / 0.35) ** 2)
    dgaussian_dr2 = gaussian * (-2.0 * dr / 0.35**2)
    dgaussian_dz = gaussian * (-2.0 * dz / 0.35**2)
    q = t * cr * cz * gaussian
    dq_dr2 = t * cz * (dcr_dr2 * gaussian + cr * dgaussian_dr2)
    dq_dz = t * cr * (dcz_dz * gaussian + cz * dgaussian_dz)
    radial_factor = -dq_dz
    ux = x * radial_factor
    uy = y * radial_factor
    uz = 2.0 * q + 2.0 * r2 * dq_dr2
    return np.stack((ux, uy, uz), axis=-1)


def perturbed_velocity(base_velocity, probe, amplitude):
    """Return a velocity callable with one fixed diagnostic coefficient."""
    amplitude = float(amplitude)
    if not np.isfinite(amplitude):
        raise ValueError("amplitude must be finite")

    def velocity(points, time):
        base = np.asarray(base_velocity(points, time), dtype=float)
        delta = np.asarray(probe(points, time), dtype=float)
        if base.shape != delta.shape or not np.all(np.isfinite(base)) or not np.all(np.isfinite(delta)):
            raise ValueError("base velocity and probe must return matching finite values")
        return base + amplitude * delta

    return velocity


__all__ = ["PEAK_RADIUS", "PEAK_Z", "perturbed_velocity", "poloidal_curl_probe", "pure_swirl_probe"]
