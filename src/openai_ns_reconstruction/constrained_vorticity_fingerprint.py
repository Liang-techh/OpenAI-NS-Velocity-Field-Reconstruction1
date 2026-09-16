"""Candidate-agnostic vorticity morphology fingerprint for visualization diagnostics.

This module samples only a public velocity(points, time) callable. The metrics are
finite-grid visualization diagnostics, not Navier--Stokes validation and not
evidence of a singularity.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Callable, Iterable

import numpy as np


@dataclass(frozen=True)
class VorticityFingerprint:
    grid_size: int
    time: float
    omega_peak: float
    omega_rms: float
    enstrophy_r50: float
    enstrophy_r90: float
    enstrophy_absz50: float
    enstrophy_absz90: float
    half_peak_volume_fraction: float

    def to_dict(self):
        return asdict(self)


def _validate_box(box):
    arr = np.asarray(box, dtype=float)
    if arr.shape != (3, 2) or not np.all(np.isfinite(arr)):
        raise ValueError("box must be finite with shape (3,2)")
    if np.any(arr[:, 1] <= arr[:, 0]):
        raise ValueError("box upper bounds must exceed lower bounds")
    return arr


def _weighted_quantile(values, weights, q):
    values = np.asarray(values, dtype=float).ravel()
    weights = np.asarray(weights, dtype=float).ravel()
    if values.shape != weights.shape or values.size == 0:
        raise ValueError("values/weights must be nonempty and shape-matched")
    if not np.all(np.isfinite(values)) or not np.all(np.isfinite(weights)):
        raise ValueError("values/weights must be finite")
    if np.any(weights < 0) or not 0.0 <= q <= 1.0:
        raise ValueError("invalid weights or quantile")
    total = float(np.sum(weights))
    if total <= 0:
        raise ValueError("positive total weight required")
    order = np.argsort(values, kind="mergesort")
    sv = values[order]
    sw = weights[order]
    cdf = np.cumsum(sw)
    idx = int(np.searchsorted(cdf, q * total, side="left"))
    idx = min(idx, sv.size - 1)
    return float(sv[idx])


def diagnose_vorticity_fingerprint(
    velocity: Callable,
    time: float,
    *,
    box=((-2.0, 2.0), (-2.0, 2.0), (-2.0, 2.0)),
    grid_size: int = 25,
) -> VorticityFingerprint:
    """Compute a finite-grid vorticity morphology fingerprint.

    The curl uses second-order finite differences on a uniform Cartesian grid.
    Enstrophy-weighted radial/axial quantiles measure where vorticity is
    geometrically concentrated. They are visualization diagnostics only.
    """
    box = _validate_box(box)
    if not np.isfinite(time):
        raise ValueError("time must be finite")
    if not isinstance(grid_size, int) or grid_size < 5:
        raise ValueError("grid_size must be an integer >=5")

    axes = [np.linspace(lo, hi, grid_size) for lo, hi in box]
    x, y, z = np.meshgrid(*axes, indexing="ij")
    points = np.stack((x, y, z), axis=-1)
    vel = np.asarray(velocity(points, time), dtype=float)
    if vel.shape != points.shape or not np.all(np.isfinite(vel)):
        raise ValueError("velocity must return finite values with shape (...,3)")

    dx, dy, dz = (axis[1] - axis[0] for axis in axes)
    u, v, w = np.moveaxis(vel, -1, 0)
    du_dx, du_dy, du_dz = np.gradient(u, dx, dy, dz, edge_order=2)
    dv_dx, dv_dy, dv_dz = np.gradient(v, dx, dy, dz, edge_order=2)
    dw_dx, dw_dy, dw_dz = np.gradient(w, dx, dy, dz, edge_order=2)
    del du_dx, dv_dy, dw_dz

    omega_x = dw_dy - dv_dz
    omega_y = du_dz - dw_dx
    omega_z = dv_dx - du_dy
    omega_mag = np.sqrt(omega_x * omega_x + omega_y * omega_y + omega_z * omega_z)
    enstrophy = omega_mag * omega_mag
    total = float(np.sum(enstrophy))
    if not np.isfinite(total) or total <= 0:
        raise ValueError("vorticity fingerprint requires positive finite enstrophy")

    radius = np.sqrt(x * x + y * y)
    absz = np.abs(z)
    peak = float(np.max(omega_mag))
    rms = float(np.sqrt(np.mean(enstrophy)))
    half_peak_fraction = float(np.mean(omega_mag >= 0.5 * peak))

    return VorticityFingerprint(
        grid_size=grid_size,
        time=float(time),
        omega_peak=peak,
        omega_rms=rms,
        enstrophy_r50=_weighted_quantile(radius, enstrophy, 0.5),
        enstrophy_r90=_weighted_quantile(radius, enstrophy, 0.9),
        enstrophy_absz50=_weighted_quantile(absz, enstrophy, 0.5),
        enstrophy_absz90=_weighted_quantile(absz, enstrophy, 0.9),
        half_peak_volume_fraction=half_peak_fraction,
    )


_METRICS = (
    "omega_peak",
    "omega_rms",
    "enstrophy_r50",
    "enstrophy_r90",
    "enstrophy_absz50",
    "enstrophy_absz90",
    "half_peak_volume_fraction",
)


def audit_vorticity_resolution(
    velocity: Callable,
    time: float,
    *,
    box=((-2.0, 2.0), (-2.0, 2.0), (-2.0, 2.0)),
    grid_sizes: Iterable[int] = (17, 25, 33),
):
    """Return three-or-more resolution levels plus relative changes to finest.

    No universal pass threshold is imposed. Consumers should compare the
    reported changes with their visualization tolerance and keep that choice
    separate from PDE acceptance.
    """
    sizes = tuple(grid_sizes)
    if len(sizes) < 3 or any(not isinstance(n, int) or n < 5 for n in sizes):
        raise ValueError("grid_sizes must contain at least three integers >=5")
    if any(b <= a for a, b in zip(sizes, sizes[1:])):
        raise ValueError("grid_sizes must be strictly increasing")

    levels = [
        diagnose_vorticity_fingerprint(velocity, time, box=box, grid_size=n)
        for n in sizes
    ]
    finest = levels[-1]
    changes = {}
    for metric in _METRICS:
        ref = float(getattr(finest, metric))
        denom = max(abs(ref), np.finfo(float).eps)
        changes[metric] = [
            float(abs(float(getattr(level, metric)) - ref) / denom)
            for level in levels[:-1]
        ]
    return {
        "time": float(time),
        "box": np.asarray(box, dtype=float).tolist(),
        "levels": [level.to_dict() for level in levels],
        "relative_change_to_finest": changes,
        "truth_boundary": (
            "finite-grid vorticity morphology only; not PDE validation, "
            "OpenAI-field identification, or blow-up evidence"
        ),
    }
