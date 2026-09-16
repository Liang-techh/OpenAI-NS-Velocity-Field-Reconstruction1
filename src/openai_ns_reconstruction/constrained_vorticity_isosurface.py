from __future__ import annotations

from collections import deque
from dataclasses import dataclass
from typing import Callable, Iterable

import numpy as np

VelocityFn = Callable[[np.ndarray, float], np.ndarray]


@dataclass(frozen=True)
class VorticityIsosurfaceFingerprint:
    time: float
    grid_size: int
    box_half_width: float
    iso_fraction: float
    omega_peak: float
    iso_level: float
    active_voxel_count: int
    active_volume_fraction: float
    component_count: int
    largest_component_fraction: float
    largest_component_volume_fraction: float
    largest_radial_extent: float
    largest_axial_span: float
    largest_axial_half_extent: float
    tip_to_midplane_area_ratio: float | None
    claim_scope: str = "visualization_morphology_only"
    visualization_ready: bool = False
    pde_validated: bool = False
    paper_exact: bool = False
    openai_field_identified: bool = False
    blowup_proved: bool = False


@dataclass(frozen=True)
class VorticityIsosurfaceResolutionAudit:
    levels: tuple[VorticityIsosurfaceFingerprint, ...]
    component_count_stable: bool
    max_largest_component_fraction_delta: float
    max_radial_extent_relative_delta: float
    max_axial_span_relative_delta: float
    max_tip_ratio_absolute_delta: float | None
    claim_scope: str = "visualization_resolution_only"
    visualization_ready: bool = False
    pde_validated: bool = False
    paper_exact: bool = False
    openai_field_identified: bool = False
    blowup_proved: bool = False


def _validate_velocity(values: np.ndarray, expected: int) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    if values.shape != (expected, 3):
        raise ValueError(f"velocity must return shape ({expected}, 3), got {values.shape}")
    if not np.all(np.isfinite(values)):
        raise ValueError("velocity returned non-finite values")
    return values


def _curl_magnitude(velocity: VelocityFn, time: float, grid_size: int, box_half_width: float):
    if not isinstance(grid_size, (int, np.integer)) or grid_size < 7 or grid_size % 2 == 0:
        raise ValueError("grid_size must be an odd integer >= 7")
    if not np.isfinite(time):
        raise ValueError("time must be finite")
    if not np.isfinite(box_half_width) or box_half_width <= 0:
        raise ValueError("box_half_width must be finite and positive")

    axis = np.linspace(-box_half_width, box_half_width, grid_size, dtype=float)
    h = float(axis[1] - axis[0])
    x, y, z = np.meshgrid(axis, axis, axis, indexing="ij")
    points = np.stack((x, y, z), axis=-1).reshape(-1, 3)
    values = _validate_velocity(velocity(points, float(time)), points.shape[0])
    field = values.reshape(grid_size, grid_size, grid_size, 3)

    # Centered second-order derivatives. Morphology is evaluated only where all
    # derivatives use a centered stencil, so boundary one-sided artifacts do not
    # define the isosurface fingerprint.
    core = (slice(1, -1), slice(1, -1), slice(1, -1))
    du_dy = (field[1:-1, 2:, 1:-1, 0] - field[1:-1, :-2, 1:-1, 0]) / (2.0 * h)
    du_dz = (field[1:-1, 1:-1, 2:, 0] - field[1:-1, 1:-1, :-2, 0]) / (2.0 * h)
    dv_dx = (field[2:, 1:-1, 1:-1, 1] - field[:-2, 1:-1, 1:-1, 1]) / (2.0 * h)
    dv_dz = (field[1:-1, 1:-1, 2:, 1] - field[1:-1, 1:-1, :-2, 1]) / (2.0 * h)
    dw_dx = (field[2:, 1:-1, 1:-1, 2] - field[:-2, 1:-1, 1:-1, 2]) / (2.0 * h)
    dw_dy = (field[1:-1, 2:, 1:-1, 2] - field[1:-1, :-2, 1:-1, 2]) / (2.0 * h)

    omega_x = dw_dy - dv_dz
    omega_y = du_dz - dw_dx
    omega_z = dv_dx - du_dy
    omega = np.sqrt(omega_x * omega_x + omega_y * omega_y + omega_z * omega_z)
    xc = x[core]
    yc = y[core]
    zc = z[core]
    return omega, xc, yc, zc, h


def _label_6(mask: np.ndarray) -> tuple[np.ndarray, list[int]]:
    if mask.ndim != 3:
        raise ValueError("mask must be three-dimensional")
    labels = np.zeros(mask.shape, dtype=np.int32)
    sizes: list[int] = []
    label = 0
    nx, ny, nz = mask.shape
    for start in np.argwhere(mask):
        i, j, k = (int(start[0]), int(start[1]), int(start[2]))
        if labels[i, j, k] != 0:
            continue
        label += 1
        q = deque([(i, j, k)])
        labels[i, j, k] = label
        count = 0
        while q:
            a, b, c = q.popleft()
            count += 1
            for da, db, dc in ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)):
                aa, bb, cc = a + da, b + db, c + dc
                if 0 <= aa < nx and 0 <= bb < ny and 0 <= cc < nz:
                    if mask[aa, bb, cc] and labels[aa, bb, cc] == 0:
                        labels[aa, bb, cc] = label
                        q.append((aa, bb, cc))
        sizes.append(count)
    return labels, sizes


def diagnose_vorticity_isosurface(
    velocity: VelocityFn,
    time: float,
    *,
    grid_size: int = 33,
    box_half_width: float = 1.5,
    iso_fraction: float = 0.5,
    activity_floor: float = 1e-12,
) -> VorticityIsosurfaceFingerprint:
    """Measure connected high-vorticity morphology on a direct velocity grid.

    ``iso_fraction`` is an autonomous visualization diagnostic, not a physical or
    PDE acceptance threshold. The returned component geometry is a finite-grid
    fingerprint of the superlevel set ``|curl u| >= iso_fraction * max|curl u|``.
    """
    if not np.isfinite(iso_fraction) or not (0.0 < iso_fraction < 1.0):
        raise ValueError("iso_fraction must lie strictly between 0 and 1")
    if not np.isfinite(activity_floor) or activity_floor <= 0:
        raise ValueError("activity_floor must be finite and positive")

    omega, x, y, z, h = _curl_magnitude(velocity, time, grid_size, box_half_width)
    peak = float(np.max(omega))
    if not np.isfinite(peak) or peak <= activity_floor:
        raise ValueError("vorticity is numerically inactive at this resolution")
    level = iso_fraction * peak
    mask = omega >= level
    active = int(np.count_nonzero(mask))
    if active == 0:
        raise ValueError("isosurface superlevel set is empty")

    labels, sizes = _label_6(mask)
    if not sizes:
        raise ValueError("isosurface component labeling failed")
    largest_label = int(np.argmax(sizes)) + 1
    largest = labels == largest_label
    largest_count = int(sizes[largest_label - 1])

    r = np.sqrt(x * x + y * y)
    r_extent = float(np.max(r[largest]))
    z_vals = z[largest]
    z_min = float(np.min(z_vals))
    z_max = float(np.max(z_vals))
    axial_span = z_max - z_min
    axial_half_extent = max(abs(z_min), abs(z_max))

    # A simple, scale-free taper diagnostic: compare the maximum active cross-
    # sectional area in the outer 20% of the largest component's |z| extent with
    # the maximum area in its central 20%. A smaller ratio means a narrower tip.
    # The 20% split is an autonomous diagnostic choice, not a source-derived fact.
    tip_ratio: float | None
    if axial_half_extent <= 2.0 * h:
        tip_ratio = None
    else:
        abs_z_axis = np.abs(z[0, 0, :])
        plane_counts = np.count_nonzero(largest, axis=(0, 1)).astype(float)
        center_sel = abs_z_axis <= 0.2 * axial_half_extent
        tip_sel = (abs_z_axis >= 0.8 * axial_half_extent) & (abs_z_axis <= axial_half_extent + 0.5 * h)
        center_max = float(np.max(plane_counts[center_sel])) if np.any(center_sel) else 0.0
        tip_max = float(np.max(plane_counts[tip_sel])) if np.any(tip_sel) else 0.0
        tip_ratio = None if center_max <= 0.0 else tip_max / center_max

    interior_voxels = int(mask.size)
    return VorticityIsosurfaceFingerprint(
        time=float(time),
        grid_size=int(grid_size),
        box_half_width=float(box_half_width),
        iso_fraction=float(iso_fraction),
        omega_peak=peak,
        iso_level=float(level),
        active_voxel_count=active,
        active_volume_fraction=active / interior_voxels,
        component_count=len(sizes),
        largest_component_fraction=largest_count / active,
        largest_component_volume_fraction=largest_count / interior_voxels,
        largest_radial_extent=r_extent,
        largest_axial_span=float(axial_span),
        largest_axial_half_extent=float(axial_half_extent),
        tip_to_midplane_area_ratio=tip_ratio,
    )


def _relative_delta(value: float, reference: float) -> float:
    scale = max(abs(reference), np.finfo(float).tiny)
    return abs(value - reference) / scale


def audit_vorticity_isosurface_resolution(
    velocity: VelocityFn,
    time: float,
    *,
    grid_sizes: Iterable[int] = (17, 25, 33),
    box_half_width: float = 1.5,
    iso_fraction: float = 0.5,
    activity_floor: float = 1e-12,
) -> VorticityIsosurfaceResolutionAudit:
    sizes = tuple(int(n) for n in grid_sizes)
    if len(sizes) < 3 or tuple(sorted(set(sizes))) != sizes:
        raise ValueError("grid_sizes must contain at least three strictly increasing values")
    levels = tuple(
        diagnose_vorticity_isosurface(
            velocity,
            time,
            grid_size=n,
            box_half_width=box_half_width,
            iso_fraction=iso_fraction,
            activity_floor=activity_floor,
        )
        for n in sizes
    )
    finest = levels[-1]
    component_stable = all(level.component_count == finest.component_count for level in levels[:-1])
    largest_delta = max(abs(level.largest_component_fraction - finest.largest_component_fraction) for level in levels[:-1])
    radial_delta = max(_relative_delta(level.largest_radial_extent, finest.largest_radial_extent) for level in levels[:-1])
    axial_delta = max(_relative_delta(level.largest_axial_span, finest.largest_axial_span) for level in levels[:-1])
    tip_values = [level.tip_to_midplane_area_ratio for level in levels]
    if finest.tip_to_midplane_area_ratio is None or any(value is None for value in tip_values[:-1]):
        tip_delta = None
    else:
        tip_delta = max(abs(float(value) - float(finest.tip_to_midplane_area_ratio)) for value in tip_values[:-1])
    return VorticityIsosurfaceResolutionAudit(
        levels=levels,
        component_count_stable=component_stable,
        max_largest_component_fraction_delta=float(largest_delta),
        max_radial_extent_relative_delta=float(radial_delta),
        max_axial_span_relative_delta=float(axial_delta),
        max_tip_ratio_absolute_delta=None if tip_delta is None else float(tip_delta),
    )
