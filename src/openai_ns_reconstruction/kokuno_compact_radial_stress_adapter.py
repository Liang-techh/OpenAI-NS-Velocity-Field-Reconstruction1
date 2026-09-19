"""Integration-only pure adapter for the admitted compact radial-stress operator.

Agent 5 needs the already-admitted Agent-3 radial inverse without importing the
oscillatory-candidate-specific transitive dependency graph owned by Agent 2.
The numerical helper bodies below are copied exactly from
``kokuno_actual_oscillatory_mean_stress.py`` at Agent-3 #677 head
``56024553b981833a283f98648c8599011d6f38ab``.  No coefficient, exponent,
quadrature rule, support rule, tolerance, or sign is retuned here.

This file is glue, not a new mathematical lane and not a new validation of the
operator.  Independent radial-operator admission remains the earlier A3/A4
#598/#607 evidence.  Final PDE validation remains separate.
"""
from __future__ import annotations

import math
from typing import Any

import numpy as np

UPSTREAM_AGENT3_PR = 677
UPSTREAM_AGENT3_HEAD = "56024553b981833a283f98648c8599011d6f38ab"
UPSTREAM_OPERATOR_MODULE = "kokuno_actual_oscillatory_mean_stress.py"


def _scalar_rms(values: np.ndarray) -> float:
    array = np.asarray(values, dtype=float)
    return float(np.sqrt(np.mean(array * array)))


def _cumulative_trapezoid(values: np.ndarray, coordinates: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    coordinates = np.asarray(coordinates, dtype=float)
    if values.ndim != 1 or coordinates.ndim != 1 or values.shape != coordinates.shape:
        raise ValueError("cumulative trapezoid requires matching one-dimensional arrays")
    if coordinates.size < 3 or np.any(np.diff(coordinates) <= 0.0):
        raise ValueError("coordinates must be strictly increasing with at least three nodes")
    out = np.zeros_like(values)
    increments = 0.5 * (values[1:] + values[:-1]) * np.diff(coordinates)
    out[1:] = np.cumsum(increments)
    return out


def _compact_cos8_bump(radii: np.ndarray, center: float, halfwidth: float) -> np.ndarray:
    radii = np.asarray(radii, dtype=float)
    s = (radii - center) / halfwidth
    out = np.zeros_like(s)
    mask = np.abs(s) < 1.0
    if np.any(mask):
        out[mask] = np.cos(0.5 * math.pi * s[mask]) ** 8
    return out


def _compact_radial_stress(
    radii: np.ndarray,
    source: np.ndarray,
    *,
    exponent: int,
    bump_center: float,
    bump_halfwidth: float,
) -> dict[str, Any]:
    """Apply the discrete moment-complement form of Kokuno's compact inverse."""
    r = np.asarray(radii, dtype=float)
    f = np.asarray(source, dtype=float)
    if r.ndim != 1 or f.shape != r.shape or r.size < 9:
        raise ValueError("radial stress requires matching one-dimensional arrays with >=9 nodes")
    if exponent not in (1, 2):
        raise ValueError("only the source axial e=1 and angular e=2 channels are supported")
    if np.any(r <= 0.0) or not np.all(np.isfinite(f)):
        raise ValueError("radii must be positive and source finite")

    weight = r ** exponent
    moment = float(_cumulative_trapezoid(weight * f, r)[-1])
    bump_raw = _compact_cos8_bump(r, bump_center, bump_halfwidth)
    bump_weighted_integral = float(_cumulative_trapezoid(weight * bump_raw, r)[-1])
    if not math.isfinite(bump_weighted_integral) or bump_weighted_integral <= 0.0:
        raise RuntimeError("compact bump lost positive weighted normalization")
    bump = bump_raw / bump_weighted_integral
    complement = f - bump * moment
    complement_moment = float(_cumulative_trapezoid(weight * complement, r)[-1])
    primitive = _cumulative_trapezoid(weight * complement, r)
    stress = -primitive / weight

    # This is the exact algebraic RHS of (d_r+e/r)sigma_e=-F+b_e M_e.
    reconstructed_force = -f + bump * moment
    return {
        "exponent": exponent,
        "weighted_moment": moment,
        "bump_weighted_integral": float(
            _cumulative_trapezoid(weight * bump, r)[-1]
        ),
        "moment_complement_weighted_moment": complement_moment,
        "stress": stress,
        "stress_rms": _scalar_rms(stress),
        "stress_max_abs": float(np.max(np.abs(stress))),
        "stress_inner_edge": float(stress[0]),
        "stress_outer_edge": float(stress[-1]),
        "reconstructed_force_rms": _scalar_rms(reconstructed_force),
    }
