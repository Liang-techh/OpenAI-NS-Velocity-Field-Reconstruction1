"""Resolution-aware cylindrical velocity morphology diagnostics.

These diagnostics consume only a public ``velocity(points, time)`` callable and
measure how much kinetic energy lies in azimuthal (swirl) versus poloidal
(radial/axial) motion. They are visualization diagnostics only: they do not
validate Navier--Stokes, identify the OpenAI numerical field, or prove blow-up.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Callable, Iterable

import numpy as np


@dataclass(frozen=True)
class CylindricalVelocityFingerprint:
    grid_size: int
    time: float
    speed2_integral: float
    swirl_energy_fraction: float
    poloidal_energy_fraction: float
    radial_fraction_of_poloidal: float | None
    swirl_to_poloidal_ratio: float | None
    pitch_deg_q10: float
    pitch_deg_q50: float
    pitch_deg_q90: float
    swirl_sign_coherence: float
    axis_transverse_energy_fraction: float

    def to_dict(self):
        return asdict(self)


def _validate_box(box):
    arr = np.asarray(box, dtype=float)
    if arr.shape != (3, 2) or not np.all(np.isfinite(arr)):
        raise ValueError("box must be finite with shape (3,2)")
    if np.any(arr[:, 1] <= arr[:, 0]):
        raise ValueError("box upper bounds must exceed lower bounds")
    if not (arr[0, 0] <= 0.0 <= arr[0, 1] and arr[1, 0] <= 0.0 <= arr[1, 1]):
        raise ValueError("x/y box must contain the symmetry axis")
    return arr


def _weighted_quantile(values, weights, q):
    values = np.asarray(values, dtype=float).ravel()
    weights = np.asarray(weights, dtype=float).ravel()
    if values.shape != weights.shape or values.size == 0:
        raise ValueError("values/weights must be nonempty and shape-matched")
    if not np.all(np.isfinite(values)) or not np.all(np.isfinite(weights)):
        raise ValueError("values/weights must be finite")
    if np.any(weights < 0.0) or not 0.0 <= q <= 1.0:
        raise ValueError("invalid weights or quantile")
    total = float(np.sum(weights))
    if total <= 0.0:
        raise ValueError("positive total weight required")
    order = np.argsort(values, kind="mergesort")
    sv = values[order]
    sw = weights[order]
    cdf = np.cumsum(sw)
    idx = int(np.searchsorted(cdf, q * total, side="left"))
    return float(sv[min(idx, sv.size - 1)])


def diagnose_cylindrical_velocity_fingerprint(
    velocity: Callable,
    time: float,
    *,
    box=((-2.0, 2.0), (-2.0, 2.0), (-2.0, 2.0)),
    grid_size: int = 25,
) -> CylindricalVelocityFingerprint:
    """Measure swirl/poloidal balance and streamline-pitch morphology.

    The field is sampled on a uniform Cartesian grid and decomposed into local
    cylindrical components away from the axis. Axis transverse energy is
    reported separately because the radial/azimuthal basis is undefined at
    r=0. The continuum axis has zero volume, so excluding it from the split is
    appropriate only when that separately reported fraction is small.
    """
    box = _validate_box(box)
    if not np.isfinite(time):
        raise ValueError("time must be finite")
    if not isinstance(grid_size, int) or grid_size < 5 or grid_size % 2 == 0:
        raise ValueError("grid_size must be an odd integer >=5")

    axes = [np.linspace(lo, hi, grid_size) for lo, hi in box]
    x, y, z = np.meshgrid(*axes, indexing="ij")
    points = np.stack((x, y, z), axis=-1)
    vel = np.asarray(velocity(points, time), dtype=float)
    if vel.shape != points.shape or not np.all(np.isfinite(vel)):
        raise ValueError("velocity must return finite values with shape (...,3)")

    u, v, w = np.moveaxis(vel, -1, 0)
    r = np.sqrt(x * x + y * y)
    speed2 = u * u + v * v + w * w

    one_d_weights = []
    for axis_values in axes:
        h = float(axis_values[1] - axis_values[0])
        weights_1d = np.full(grid_size, h, dtype=float)
        weights_1d[0] *= 0.5
        weights_1d[-1] *= 0.5
        one_d_weights.append(weights_1d)
    cell_weights = (
        one_d_weights[0][:, None, None]
        * one_d_weights[1][None, :, None]
        * one_d_weights[2][None, None, :]
    )
    total = float(np.sum(cell_weights * speed2))
    if not np.isfinite(total) or total <= 0.0:
        raise ValueError("fingerprint requires positive finite kinetic energy")

    step = min(
        float(abs(axes[0][1] - axes[0][0])),
        float(abs(axes[1][1] - axes[1][0])),
    )
    axis_tol = max(np.finfo(float).eps * 32.0, step * 1.0e-12)
    off = r > axis_tol
    axis = ~off

    ur = np.zeros_like(u)
    utheta = np.zeros_like(u)
    ur[off] = (x[off] * u[off] + y[off] * v[off]) / r[off]
    utheta[off] = (-y[off] * u[off] + x[off] * v[off]) / r[off]

    swirl_e = utheta * utheta
    radial_e = ur * ur
    axial_e = w * w
    poloidal_e = radial_e + axial_e

    split_total = float(np.sum(cell_weights * (swirl_e + poloidal_e)))
    if split_total <= 0.0:
        raise ValueError("off-axis cylindrical decomposition has zero energy")
    swirl_sum = float(np.sum(cell_weights * swirl_e))
    poloidal_sum = float(np.sum(cell_weights * poloidal_e))
    swirl_fraction = swirl_sum / split_total
    poloidal_fraction = poloidal_sum / split_total

    floor = np.finfo(float).eps * max(split_total, 1.0) * 64.0
    ratio = None if poloidal_sum <= floor else swirl_sum / poloidal_sum
    radial_fraction = (
        None
        if poloidal_sum <= floor
        else float(np.sum(cell_weights * radial_e)) / poloidal_sum
    )

    pitch = np.degrees(np.arctan2(np.sqrt(poloidal_e[off]), np.abs(utheta[off])))
    pitch_weights = (cell_weights * (swirl_e + poloidal_e))[off]
    if float(np.sum(pitch_weights)) <= 0.0:
        raise ValueError("positive off-axis energy required for pitch statistics")

    abs_swirl = np.abs(utheta[off])
    swirl_linear_weights = cell_weights[off]
    abs_swirl_sum = float(np.sum(swirl_linear_weights * abs_swirl))
    coherence = (
        0.0
        if abs_swirl_sum <= floor
        else abs(float(np.sum(swirl_linear_weights * utheta[off]))) / abs_swirl_sum
    )

    axis_transverse = float(
        np.sum((cell_weights * (u * u + v * v))[axis])
    ) / total

    return CylindricalVelocityFingerprint(
        grid_size=grid_size,
        time=float(time),
        speed2_integral=total,
        swirl_energy_fraction=float(swirl_fraction),
        poloidal_energy_fraction=float(poloidal_fraction),
        radial_fraction_of_poloidal=(
            None if radial_fraction is None else float(radial_fraction)
        ),
        swirl_to_poloidal_ratio=None if ratio is None else float(ratio),
        pitch_deg_q10=_weighted_quantile(pitch, pitch_weights, 0.10),
        pitch_deg_q50=_weighted_quantile(pitch, pitch_weights, 0.50),
        pitch_deg_q90=_weighted_quantile(pitch, pitch_weights, 0.90),
        swirl_sign_coherence=float(coherence),
        axis_transverse_energy_fraction=float(axis_transverse),
    )


_METRICS = (
    "swirl_energy_fraction",
    "poloidal_energy_fraction",
    "pitch_deg_q10",
    "pitch_deg_q50",
    "pitch_deg_q90",
    "swirl_sign_coherence",
    "axis_transverse_energy_fraction",
)


def audit_cylindrical_velocity_resolution(
    velocity: Callable,
    time: float,
    *,
    box=((-2.0, 2.0), (-2.0, 2.0), (-2.0, 2.0)),
    grid_sizes: Iterable[int] = (17, 25, 33),
):
    """Return at least three resolutions and absolute changes to the finest.

    Fractions and pitch angles are compared by absolute change, avoiding
    unstable relative errors near physically meaningful zeros. No universal
    pass threshold is imposed.
    """
    sizes = tuple(grid_sizes)
    if len(sizes) < 3:
        raise ValueError("grid_sizes must contain at least three levels")
    if any(not isinstance(n, int) or n < 5 or n % 2 == 0 for n in sizes):
        raise ValueError("grid_sizes must contain odd integers >=5")
    if any(b <= a for a, b in zip(sizes, sizes[1:])):
        raise ValueError("grid_sizes must be strictly increasing")

    levels = [
        diagnose_cylindrical_velocity_fingerprint(
            velocity, time, box=box, grid_size=n
        )
        for n in sizes
    ]
    finest = levels[-1]
    changes = {}
    for metric in _METRICS:
        ref = float(getattr(finest, metric))
        changes[metric] = [
            float(abs(float(getattr(level, metric)) - ref))
            for level in levels[:-1]
        ]

    return {
        "time": float(time),
        "box": np.asarray(box, dtype=float).tolist(),
        "levels": [level.to_dict() for level in levels],
        "absolute_change_to_finest": changes,
        "truth_boundary": (
            "finite-grid cylindrical velocity morphology only; not PDE "
            "validation, OpenAI-field identification, or blow-up evidence"
        ),
    }
