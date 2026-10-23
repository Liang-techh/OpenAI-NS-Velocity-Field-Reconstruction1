"""Truth-bounded axial vortex-core probe based on parallel-vector flow geometry.

This module consumes only a public ``velocity(points, time) -> [u,v,w]``
callable. It uses a fixed Cartesian grid and the local velocity gradient to
construct a visualization-only axial probe trace. The construction is inspired
by parallel-vector vortex-core methods, but it is deliberately not a full
Sujudi-Haimes/VTK core-line extractor and is never PDE acceptance evidence.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np


VelocityCallable = Callable[[np.ndarray, float], np.ndarray]


@dataclass(frozen=True)
class ParallelVectorCoreTrace:
    time: float
    grid_size: int
    box_half_width: float
    spacing: float
    trace_points: tuple[tuple[float, float, float], ...]
    trace_swirl_strength: tuple[float, ...]
    trace_parallel_residual: tuple[float, ...]
    trace_local_parallel_sine: tuple[float | None, ...]
    trace_point_count: int
    valid_slice_fraction: float
    axial_span: float
    transverse_radius_rms: float | None
    transverse_step_rms: float | None
    trace_swirl_strength_median: float | None
    trace_parallel_residual_median: float | None
    trace_parallel_residual_max: float | None
    trace_local_parallel_sine_median: float | None
    velocity_rms_speed: float
    acceleration_rms: float
    swirling_strength_rms: float
    swirling_strength_max: float
    divergence_rms: float
    provenance: str
    method_scope: str = "axial_visualization_probe_only"
    velocity_changed: bool = False
    visualization_ready: bool = False
    visual_correspondence_verified: bool = False
    pde_validated: bool = False
    paper_exact: bool = False
    openai_field_identified: bool = False
    blowup_proved: bool = False


def _validate_inputs(time: float, box_half_width: float, grid_size: int, provenance: str) -> None:
    if not np.isfinite(time) or not (0.25 <= float(time) <= 0.75):
        raise ValueError("time must be finite and inside the governed [0.25,0.75] interval")
    if not np.isfinite(box_half_width) or not (0.0 < float(box_half_width) <= 2.0):
        raise ValueError(
            "box_half_width must be finite, positive, and no larger than the governed box half-width 2"
        )
    if (
        isinstance(grid_size, bool)
        or int(grid_size) != grid_size
        or int(grid_size) < 7
        or int(grid_size) % 2 == 0
    ):
        raise ValueError("grid_size must be an odd integer >= 7")
    if not isinstance(provenance, str) or not provenance.strip():
        raise ValueError("provenance must be a non-empty string")


def _evaluate_velocity(velocity: VelocityCallable, points: np.ndarray, time: float) -> np.ndarray:
    raw = np.asarray(velocity(points, float(time)), dtype=float)
    if raw.shape != points.shape:
        raise ValueError("velocity(points,time) must return shape (N,3)")
    if not np.isfinite(raw).all():
        raise ValueError("velocity output must be finite")
    return raw


def _swirling_strength(gradient: np.ndarray) -> np.ndarray:
    flat = np.asarray(gradient, dtype=float).reshape(-1, 3, 3)
    if not np.isfinite(flat).all():
        raise ValueError("velocity gradient must be finite")
    eigvals = np.linalg.eigvals(flat)
    lam_ci = np.max(np.abs(eigvals.imag), axis=1)
    scale = np.maximum(1.0, np.linalg.norm(flat, axis=(1, 2)))
    numerical_guard = 64.0 * np.finfo(float).eps * scale
    lam_ci[lam_ci <= numerical_guard] = 0.0
    return lam_ci.reshape(gradient.shape[:-2])


def diagnose_parallel_vector_core_trace(
    velocity: VelocityCallable,
    time: float,
    *,
    box_half_width: float = 1.5,
    grid_size: int = 17,
    provenance: str,
) -> ParallelVectorCoreTrace:
    """Return a deterministic axial probe trace for local swirling geometry.

    At each interior z-slice, the routine finds the numerically strongest local
    swirling-strength level. Within machine-precision ties on that level, it
    chooses the point with the smallest dimensionless parallel-vector residual
    ``|u x ((grad u)u)| / (u_rms * a_rms)``. This makes no visual pass/fail
    decision and does not fit any camera, hidden time, candidate coefficient,
    pressure, or forcing.

    The returned trace is only an axial visualization probe. It is not a full
    cell-intersection vortex-core line and cannot establish PDE validity or
    correspondence to any hidden OpenAI field.
    """
    _validate_inputs(time, box_half_width, grid_size, provenance)
    n = int(grid_size)
    axis = np.linspace(-float(box_half_width), float(box_half_width), n, dtype=float)
    spacing = float(axis[1] - axis[0])
    xx, yy, zz = np.meshgrid(axis, axis, axis, indexing="ij")
    points = np.stack((xx, yy, zz), axis=-1).reshape(-1, 3)
    sampled = _evaluate_velocity(velocity, points, float(time)).reshape(n, n, n, 3)

    center = sampled[1:-1, 1:-1, 1:-1]
    speed_sq = np.sum(center * center, axis=-1)
    velocity_rms_speed = float(np.sqrt(np.mean(speed_sq)))
    if velocity_rms_speed == 0.0:
        raise ValueError("exact-zero interior velocity is not a valid visualization candidate")

    d_dx = (sampled[2:, 1:-1, 1:-1] - sampled[:-2, 1:-1, 1:-1]) / (2.0 * spacing)
    d_dy = (sampled[1:-1, 2:, 1:-1] - sampled[1:-1, :-2, 1:-1]) / (2.0 * spacing)
    d_dz = (sampled[1:-1, 1:-1, 2:] - sampled[1:-1, 1:-1, :-2]) / (2.0 * spacing)

    gradient = np.empty(center.shape[:-1] + (3, 3), dtype=float)
    gradient[..., :, 0] = d_dx
    gradient[..., :, 1] = d_dy
    gradient[..., :, 2] = d_dz

    swirling = _swirling_strength(gradient)
    swirling_strength_rms = float(np.sqrt(np.mean(swirling * swirling)))
    swirling_strength_max = float(np.max(swirling))

    acceleration = np.einsum("...ij,...j->...i", gradient, center)
    acceleration_rms = float(np.sqrt(np.mean(np.sum(acceleration * acceleration, axis=-1))))
    cross_norm = np.linalg.norm(np.cross(center, acceleration), axis=-1)

    global_parallel_scale = velocity_rms_speed * acceleration_rms
    if global_parallel_scale > 0.0:
        parallel_residual = cross_norm / global_parallel_scale
    else:
        parallel_residual = np.zeros_like(cross_norm)

    local_denom = np.linalg.norm(center, axis=-1) * np.linalg.norm(acceleration, axis=-1)
    local_sine = np.full_like(local_denom, np.nan)
    local_guard = 128.0 * np.finfo(float).eps * max(1.0, global_parallel_scale)
    local_active = local_denom > local_guard
    local_sine[local_active] = np.minimum(
        1.0, cross_norm[local_active] / local_denom[local_active]
    )

    divergence = gradient[..., 0, 0] + gradient[..., 1, 1] + gradient[..., 2, 2]
    divergence_rms = float(np.sqrt(np.mean(divergence * divergence)))

    interior_axis = axis[1:-1]
    trace_points: list[tuple[float, float, float]] = []
    trace_swirl: list[float] = []
    trace_residual: list[float] = []
    trace_sine: list[float | None] = []

    eps = np.finfo(float).eps
    for k, z in enumerate(interior_axis):
        lam_slice = swirling[:, :, k]
        max_lam = float(np.max(lam_slice))
        if max_lam <= 0.0:
            continue

        tie_tol = 128.0 * eps * max(1.0, max_lam)
        tied = (max_lam - lam_slice) <= tie_tol
        ranked = np.where(tied, parallel_residual[:, :, k], np.inf)
        flat_index = int(np.argmin(ranked))
        i, j = np.unravel_index(flat_index, ranked.shape)

        point = (float(interior_axis[i]), float(interior_axis[j]), float(z))
        trace_points.append(point)
        trace_swirl.append(float(lam_slice[i, j]))
        trace_residual.append(float(parallel_residual[i, j, k]))
        sine_value = (
            float(local_sine[i, j, k]) if np.isfinite(local_sine[i, j, k]) else None
        )
        trace_sine.append(sine_value)

    count = len(trace_points)
    valid_slice_fraction = float(count / (n - 2))
    if count:
        trace_array = np.asarray(trace_points, dtype=float)
        axial_span = float(trace_array[-1, 2] - trace_array[0, 2]) if count > 1 else 0.0
        transverse_radius_rms = float(
            np.sqrt(np.mean(trace_array[:, 0] ** 2 + trace_array[:, 1] ** 2))
        )
        if count > 1:
            transverse_steps = np.diff(trace_array[:, :2], axis=0)
            transverse_step_rms = float(
                np.sqrt(np.mean(np.sum(transverse_steps * transverse_steps, axis=1)))
            )
        else:
            transverse_step_rms = None
        swirl_median = float(np.median(np.asarray(trace_swirl)))
        residual_median = float(np.median(np.asarray(trace_residual)))
        residual_max = float(np.max(np.asarray(trace_residual)))
        finite_sines = np.asarray(
            [value for value in trace_sine if value is not None], dtype=float
        )
        sine_median = float(np.median(finite_sines)) if finite_sines.size else None
    else:
        axial_span = 0.0
        transverse_radius_rms = None
        transverse_step_rms = None
        swirl_median = None
        residual_median = None
        residual_max = None
        sine_median = None

    return ParallelVectorCoreTrace(
        time=float(time),
        grid_size=n,
        box_half_width=float(box_half_width),
        spacing=spacing,
        trace_points=tuple(trace_points),
        trace_swirl_strength=tuple(trace_swirl),
        trace_parallel_residual=tuple(trace_residual),
        trace_local_parallel_sine=tuple(trace_sine),
        trace_point_count=count,
        valid_slice_fraction=valid_slice_fraction,
        axial_span=axial_span,
        transverse_radius_rms=transverse_radius_rms,
        transverse_step_rms=transverse_step_rms,
        trace_swirl_strength_median=swirl_median,
        trace_parallel_residual_median=residual_median,
        trace_parallel_residual_max=residual_max,
        trace_local_parallel_sine_median=sine_median,
        velocity_rms_speed=velocity_rms_speed,
        acceleration_rms=acceleration_rms,
        swirling_strength_rms=swirling_strength_rms,
        swirling_strength_max=swirling_strength_max,
        divergence_rms=divergence_rms,
        provenance=provenance.strip(),
    )
