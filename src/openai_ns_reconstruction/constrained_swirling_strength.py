"""Truth-bounded 3-D swirling-strength morphology diagnostic.

This module consumes only a public ``velocity(points, time) -> [u,v,w]``
callable.  It reconstructs the Cartesian velocity-gradient tensor on a fixed
uniform grid and measures the imaginary part of its complex eigenvalue pair.
The result is a visualization/morphology fingerprint, not PDE acceptance.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Iterable

import numpy as np


VelocityCallable = Callable[[np.ndarray, float], np.ndarray]


@dataclass(frozen=True)
class SwirlingStrengthFingerprint:
    time: float
    grid_size: int
    box_half_width: float
    spacing: float
    velocity_rms_speed: float
    swirling_strength_rms: float
    swirling_strength_max: float
    swirling_strength_integral: float
    swirling_detected: bool
    active_fraction: float
    swirl_weighted_centroid: tuple[float, float, float] | None
    swirl_weighted_radius_rms: float | None
    swirl_weighted_abs_z_mean: float | None
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
class SwirlingStrengthResolutionAudit:
    time: float
    grid_sizes: tuple[int, ...]
    fingerprints: tuple[SwirlingStrengthFingerprint, ...]
    finest_grid_size: int
    swirling_rms_abs_delta_to_finest: tuple[float, ...]
    swirling_max_abs_delta_to_finest: tuple[float, ...]
    active_fraction_abs_delta_to_finest: tuple[float, ...]
    radius_rms_abs_delta_to_finest: tuple[float | None, ...]
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


def _evaluate_velocity(velocity: VelocityCallable, points: np.ndarray, time: float) -> np.ndarray:
    raw = np.asarray(velocity(points, time), dtype=float)
    if raw.shape != points.shape:
        raise ValueError("velocity(points,time) must return shape (N,3)")
    if not np.isfinite(raw).all():
        raise ValueError("velocity output must be finite")
    return raw


def _swirling_from_gradient(gradient: np.ndarray) -> np.ndarray:
    """Return lambda_ci for real 3x3 velocity-gradient tensors.

    Tiny imaginary parts no larger than a fixed machine-precision guard are
    zeroed.  This guard is purely numerical and is never used as a candidate
    acceptance threshold.
    """
    flat = np.asarray(gradient, dtype=float).reshape(-1, 3, 3)
    if not np.isfinite(flat).all():
        raise ValueError("velocity gradient must be finite")
    eigenvalues = np.linalg.eigvals(flat)
    imag_max = np.max(np.abs(eigenvalues.imag), axis=1)
    scale = np.maximum(1.0, np.linalg.norm(flat, axis=(1, 2)))
    numerical_guard = 64.0 * np.finfo(float).eps * scale
    imag_max[imag_max <= numerical_guard] = 0.0
    return imag_max.reshape(gradient.shape[:-2])


def diagnose_swirling_strength(
    velocity: VelocityCallable,
    time: float,
    *,
    box_half_width: float = 1.5,
    grid_size: int = 17,
    provenance: str,
) -> SwirlingStrengthFingerprint:
    """Measure intrinsic local swirling strength on one fixed Cartesian grid.

    ``lambda_ci`` is the magnitude of the imaginary part of the complex
    conjugate eigenvalue pair of ``grad(u)``.  No image threshold, registration,
    field rescaling, or candidate modification is performed.
    """
    _validate_common(time, box_half_width, grid_size, provenance)

    axis = np.linspace(-box_half_width, box_half_width, grid_size, dtype=float)
    spacing = float(axis[1] - axis[0])
    xx, yy, zz = np.meshgrid(axis, axis, axis, indexing="ij")
    points = np.stack((xx, yy, zz), axis=-1).reshape(-1, 3)
    sampled = _evaluate_velocity(velocity, points, float(time)).reshape(grid_size, grid_size, grid_size, 3)

    center = sampled[1:-1, 1:-1, 1:-1]
    speed_sq = np.sum(center * center, axis=-1)
    velocity_rms_speed = float(np.sqrt(np.mean(speed_sq)))
    if velocity_rms_speed == 0.0:
        raise ValueError("exact-zero interior velocity is not a valid morphology candidate")

    d_dx = (sampled[2:, 1:-1, 1:-1] - sampled[:-2, 1:-1, 1:-1]) / (2.0 * spacing)
    d_dy = (sampled[1:-1, 2:, 1:-1] - sampled[1:-1, :-2, 1:-1]) / (2.0 * spacing)
    d_dz = (sampled[1:-1, 1:-1, 2:] - sampled[1:-1, 1:-1, :-2]) / (2.0 * spacing)

    gradient = np.empty(center.shape[:-1] + (3, 3), dtype=float)
    gradient[..., :, 0] = d_dx
    gradient[..., :, 1] = d_dy
    gradient[..., :, 2] = d_dz

    swirling = _swirling_from_gradient(gradient)
    swirling_strength_rms = float(np.sqrt(np.mean(swirling * swirling)))
    swirling_strength_max = float(np.max(swirling))
    swirling_strength_integral = float(np.sum(swirling) * spacing**3)
    active = swirling > 0.0
    active_fraction = float(np.mean(active))
    swirling_detected = bool(np.any(active))

    divergence = gradient[..., 0, 0] + gradient[..., 1, 1] + gradient[..., 2, 2]
    divergence_rms = float(np.sqrt(np.mean(divergence * divergence)))

    centroid: tuple[float, float, float] | None
    radius_rms: float | None
    abs_z_mean: float | None
    total_weight = float(np.sum(swirling))
    if total_weight > 0.0:
        interior_axis = axis[1:-1]
        xi, yi, zi = np.meshgrid(interior_axis, interior_axis, interior_axis, indexing="ij")
        cx = float(np.sum(swirling * xi) / total_weight)
        cy = float(np.sum(swirling * yi) / total_weight)
        cz = float(np.sum(swirling * zi) / total_weight)
        centroid = (cx, cy, cz)
        radius_rms = float(np.sqrt(np.sum(swirling * (xi * xi + yi * yi)) / total_weight))
        abs_z_mean = float(np.sum(swirling * np.abs(zi)) / total_weight)
    else:
        centroid = None
        radius_rms = None
        abs_z_mean = None

    return SwirlingStrengthFingerprint(
        time=float(time),
        grid_size=int(grid_size),
        box_half_width=float(box_half_width),
        spacing=spacing,
        velocity_rms_speed=velocity_rms_speed,
        swirling_strength_rms=swirling_strength_rms,
        swirling_strength_max=swirling_strength_max,
        swirling_strength_integral=swirling_strength_integral,
        swirling_detected=swirling_detected,
        active_fraction=active_fraction,
        swirl_weighted_centroid=centroid,
        swirl_weighted_radius_rms=radius_rms,
        swirl_weighted_abs_z_mean=abs_z_mean,
        divergence_rms=divergence_rms,
        provenance=provenance.strip(),
    )


def audit_swirling_strength_resolution(
    velocity: VelocityCallable,
    time: float,
    grid_sizes: Iterable[int],
    *,
    box_half_width: float = 1.5,
    provenance: str,
) -> SwirlingStrengthResolutionAudit:
    """Report grid sensitivity against the finest caller-declared grid.

    This function intentionally reports discrepancies only.  It defines no
    post-hoc pass/fail threshold for visualization correspondence.
    """
    sizes = tuple(grid_sizes)
    if len(sizes) < 3:
        raise ValueError("grid_sizes must contain at least three resolutions")
    for size in sizes:
        _validate_common(time, box_half_width, size, provenance)
    if tuple(sorted(set(sizes))) != sizes:
        raise ValueError("grid_sizes must be strictly increasing with no duplicates")

    fingerprints = tuple(
        diagnose_swirling_strength(
            velocity,
            time,
            box_half_width=box_half_width,
            grid_size=size,
            provenance=f"{provenance.strip()}; grid_size={size}",
        )
        for size in sizes
    )
    finest = fingerprints[-1]

    def radius_delta(item: SwirlingStrengthFingerprint) -> float | None:
        if item.swirl_weighted_radius_rms is None or finest.swirl_weighted_radius_rms is None:
            return None
        return abs(item.swirl_weighted_radius_rms - finest.swirl_weighted_radius_rms)

    return SwirlingStrengthResolutionAudit(
        time=float(time),
        grid_sizes=tuple(int(s) for s in sizes),
        fingerprints=fingerprints,
        finest_grid_size=int(sizes[-1]),
        swirling_rms_abs_delta_to_finest=tuple(
            abs(item.swirling_strength_rms - finest.swirling_strength_rms) for item in fingerprints
        ),
        swirling_max_abs_delta_to_finest=tuple(
            abs(item.swirling_strength_max - finest.swirling_strength_max) for item in fingerprints
        ),
        active_fraction_abs_delta_to_finest=tuple(
            abs(item.active_fraction - finest.active_fraction) for item in fingerprints
        ),
        radius_rms_abs_delta_to_finest=tuple(radius_delta(item) for item in fingerprints),
        provenance=provenance.strip(),
    )
