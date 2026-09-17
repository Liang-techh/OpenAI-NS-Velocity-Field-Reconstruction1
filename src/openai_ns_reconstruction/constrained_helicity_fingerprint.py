"""Truth-bounded helicity/alignment morphology diagnostic for callable 3-D velocity fields."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Iterable

import numpy as np


VelocityCallable = Callable[[np.ndarray, float], np.ndarray]


@dataclass(frozen=True)
class HelicityFingerprint:
    time: float
    grid_size: int
    box_half_width: float
    spacing: float
    active_fraction: float
    helicity_integral: float
    helicity_rms: float
    signed_alignment: float
    absolute_alignment: float
    positive_alignment_weight_fraction: float
    negative_alignment_weight_fraction: float
    alignment_q10: float
    alignment_q50: float
    alignment_q90: float
    divergence_rms: float
    provenance: str
    velocity_changed: bool = False
    visualization_ready: bool = False
    visual_correspondence_verified: bool = False
    pde_validated: bool = False
    paper_exact: bool = False
    openai_field_identified: bool = False
    blowup_proved: bool = False


@dataclass(frozen=True)
class HelicityResolutionAudit:
    time: float
    grid_sizes: tuple[int, ...]
    fingerprints: tuple[HelicityFingerprint, ...]
    finest_grid_size: int
    signed_alignment_abs_delta_to_finest: tuple[float, ...]
    absolute_alignment_abs_delta_to_finest: tuple[float, ...]
    helicity_integral_rel_delta_to_finest: tuple[float, ...]
    provenance: str
    visualization_ready: bool = False
    visual_correspondence_verified: bool = False
    pde_validated: bool = False
    paper_exact: bool = False
    openai_field_identified: bool = False
    blowup_proved: bool = False


def _validate_common(time: float, box_half_width: float, grid_size: int, provenance: str) -> None:
    if not np.isfinite(time):
        raise ValueError("time must be finite")
    if not np.isfinite(box_half_width) or box_half_width <= 0:
        raise ValueError("box_half_width must be finite and positive")
    if isinstance(grid_size, bool) or int(grid_size) != grid_size or grid_size < 7 or grid_size % 2 == 0:
        raise ValueError("grid_size must be an odd integer >= 7")
    if not isinstance(provenance, str) or not provenance.strip():
        raise ValueError("provenance must be a non-empty string")


def _weighted_quantile(values: np.ndarray, weights: np.ndarray, quantiles: Iterable[float]) -> np.ndarray:
    values = np.asarray(values, dtype=float).reshape(-1)
    weights = np.asarray(weights, dtype=float).reshape(-1)
    quantiles = np.asarray(tuple(quantiles), dtype=float)
    if values.size == 0 or values.shape != weights.shape:
        raise ValueError("weighted quantile inputs must be non-empty and shape-matched")
    if not np.isfinite(values).all() or not np.isfinite(weights).all() or np.any(weights < 0):
        raise ValueError("weighted quantile inputs must be finite with nonnegative weights")
    total = float(np.sum(weights))
    if total <= 0:
        raise ValueError("weighted quantile total weight must be positive")
    order = np.argsort(values, kind="mergesort")
    v = values[order]
    w = weights[order]
    centers = (np.cumsum(w) - 0.5 * w) / total
    return np.interp(quantiles, centers, v, left=v[0], right=v[-1])


def _evaluate_velocity(velocity: VelocityCallable, points: np.ndarray, time: float) -> np.ndarray:
    raw = np.asarray(velocity(points, time), dtype=float)
    if raw.shape != points.shape:
        raise ValueError("velocity(points,time) must return shape (N,3)")
    if not np.isfinite(raw).all():
        raise ValueError("velocity output must be finite")
    return raw


def diagnose_helicity_fingerprint(
    velocity: VelocityCallable,
    time: float,
    *,
    box_half_width: float = 1.5,
    grid_size: int = 17,
    activity_floor_relative: float = 1e-12,
    provenance: str,
) -> HelicityFingerprint:
    """Measure velocity-vorticity alignment on one fixed Cartesian visualization grid.

    This is a visualization/morphology diagnostic, not a PDE validator.  Vorticity is
    reconstructed independently from public velocity samples using centered second-order
    finite differences.  No registration, camera fit, component rescaling, or field
    modification is performed.
    """
    _validate_common(time, box_half_width, grid_size, provenance)
    if not np.isfinite(activity_floor_relative) or not (0 <= activity_floor_relative < 1):
        raise ValueError("activity_floor_relative must lie in [0,1)")

    axis = np.linspace(-box_half_width, box_half_width, grid_size, dtype=float)
    h = float(axis[1] - axis[0])
    xx, yy, zz = np.meshgrid(axis, axis, axis, indexing="ij")
    points = np.stack((xx, yy, zz), axis=-1).reshape(-1, 3)
    sampled = _evaluate_velocity(velocity, points, float(time)).reshape(grid_size, grid_size, grid_size, 3)

    center = sampled[1:-1, 1:-1, 1:-1]
    dudx = (sampled[2:, 1:-1, 1:-1] - sampled[:-2, 1:-1, 1:-1]) / (2.0 * h)
    dudy = (sampled[1:-1, 2:, 1:-1] - sampled[1:-1, :-2, 1:-1]) / (2.0 * h)
    dudz = (sampled[1:-1, 1:-1, 2:] - sampled[1:-1, 1:-1, :-2]) / (2.0 * h)

    omega = np.empty_like(center)
    omega[..., 0] = dudy[..., 2] - dudz[..., 1]
    omega[..., 1] = dudz[..., 0] - dudx[..., 2]
    omega[..., 2] = dudx[..., 1] - dudy[..., 0]
    divergence = dudx[..., 0] + dudy[..., 1] + dudz[..., 2]

    speed = np.linalg.norm(center, axis=-1)
    omega_mag = np.linalg.norm(omega, axis=-1)
    max_speed = float(np.max(speed))
    max_omega = float(np.max(omega_mag))
    if max_speed <= 0:
        raise ValueError("sampled interior velocity is identically zero")
    if max_omega <= 0:
        raise ValueError("sampled interior vorticity is identically zero")

    active = (speed > activity_floor_relative * max_speed) & (omega_mag > activity_floor_relative * max_omega)
    if not np.any(active):
        raise ValueError("no velocity-vorticity samples survive the declared activity floor")

    helicity = np.einsum("...i,...i->...", center, omega)
    denom = speed * omega_mag
    active_h = helicity[active]
    active_w = denom[active]
    total_w = float(np.sum(active_w))
    if not np.isfinite(total_w) or total_w <= 0:
        raise ValueError("active velocity-vorticity weight is degenerate")

    alignment = np.clip(active_h / active_w, -1.0, 1.0)
    q10, q50, q90 = _weighted_quantile(alignment, active_w, (0.10, 0.50, 0.90))
    positive_weight = float(np.sum(active_w[alignment > 0]))
    negative_weight = float(np.sum(active_w[alignment < 0]))

    return HelicityFingerprint(
        time=float(time),
        grid_size=int(grid_size),
        box_half_width=float(box_half_width),
        spacing=h,
        active_fraction=float(np.mean(active)),
        helicity_integral=float(np.sum(helicity) * h**3),
        helicity_rms=float(np.sqrt(np.mean(helicity**2))),
        signed_alignment=float(np.sum(active_h) / total_w),
        absolute_alignment=float(np.sum(np.abs(active_h)) / total_w),
        positive_alignment_weight_fraction=positive_weight / total_w,
        negative_alignment_weight_fraction=negative_weight / total_w,
        alignment_q10=float(q10),
        alignment_q50=float(q50),
        alignment_q90=float(q90),
        divergence_rms=float(np.sqrt(np.mean(divergence**2))),
        provenance=provenance.strip(),
    )


def audit_helicity_resolution(
    velocity: VelocityCallable,
    time: float,
    *,
    grid_sizes: Iterable[int] = (17, 25, 33),
    box_half_width: float = 1.5,
    activity_floor_relative: float = 1e-12,
    provenance: str,
) -> HelicityResolutionAudit:
    """Repeat the helicity fingerprint on >=3 declared grids without a pass threshold."""
    sizes = tuple(grid_sizes)
    if len(sizes) < 3:
        raise ValueError("at least three grid sizes are required")
    for size in sizes:
        _validate_common(time, box_half_width, size, provenance)
    if any(b <= a for a, b in zip(sizes, sizes[1:])):
        raise ValueError("grid_sizes must be strictly increasing")

    fps = tuple(
        diagnose_helicity_fingerprint(
            velocity,
            time,
            box_half_width=box_half_width,
            grid_size=size,
            activity_floor_relative=activity_floor_relative,
            provenance=provenance,
        )
        for size in sizes
    )
    finest = fps[-1]
    denom = max(abs(finest.helicity_integral), np.finfo(float).tiny)
    return HelicityResolutionAudit(
        time=float(time),
        grid_sizes=tuple(int(v) for v in sizes),
        fingerprints=fps,
        finest_grid_size=int(sizes[-1]),
        signed_alignment_abs_delta_to_finest=tuple(abs(fp.signed_alignment - finest.signed_alignment) for fp in fps),
        absolute_alignment_abs_delta_to_finest=tuple(abs(fp.absolute_alignment - finest.absolute_alignment) for fp in fps),
        helicity_integral_rel_delta_to_finest=tuple(abs(fp.helicity_integral - finest.helicity_integral) / denom for fp in fps),
        provenance=provenance.strip(),
    )
