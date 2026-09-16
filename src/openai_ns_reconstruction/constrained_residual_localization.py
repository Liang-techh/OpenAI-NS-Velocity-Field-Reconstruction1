from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Iterable

import numpy as np


ResidualFn = Callable[[np.ndarray, float], np.ndarray]


@dataclass(frozen=True)
class ResidualLocalizationSnapshot:
    time: float
    grid_size: int
    total_rms: float
    total_max: float
    cylindrical_component_rms: dict[str, float]
    cylindrical_component_max: dict[str, float]
    cylindrical_component_energy_fraction: dict[str, float]
    region_energy_fraction: dict[str, float]
    region_component_energy_fraction: dict[str, dict[str, float]]
    region_thresholds: dict[str, float]
    sampled_volume: float
    claim_scope: str = "residual_localization_only"
    pde_validated: bool = False
    paper_exact: bool = False


@dataclass(frozen=True)
class ResidualLocalizationResolutionAudit:
    snapshots: tuple[ResidualLocalizationSnapshot, ...]
    max_region_fraction_delta_to_finest: tuple[float, ...]
    max_component_fraction_delta_to_finest: tuple[float, ...]
    relative_total_rms_delta_to_finest: tuple[float, ...]
    claim_scope: str = "residual_localization_resolution_only"
    pde_validated: bool = False
    paper_exact: bool = False


def _validate_scalar(name: str, value: float, *, positive: bool = False) -> float:
    value = float(value)
    if not np.isfinite(value) or (positive and value <= 0.0):
        qualifier = "positive finite" if positive else "finite"
        raise ValueError(f"{name} must be {qualifier}")
    return value


def _validate_fraction(name: str, value: float) -> float:
    value = _validate_scalar(name, value)
    if not (0.0 < value < 1.0):
        raise ValueError(f"{name} must lie strictly between 0 and 1")
    return value


def _cell_centers(half_width: float, count: int) -> np.ndarray:
    edges = np.linspace(-half_width, half_width, count + 1, dtype=float)
    return 0.5 * (edges[:-1] + edges[1:])


def diagnose_residual_localization(
    residual: ResidualFn,
    *,
    time: float,
    radial_half_width: float,
    axial_half_height: float,
    grid_size: int,
    core_radius_fraction: float = 0.35,
    collar_start_fraction: float = 0.75,
    tip_start_fraction: float = 0.75,
) -> ResidualLocalizationSnapshot:
    """Localize an already-independent momentum residual on a finite box.

    The callback must return Cartesian momentum residual samples with shape (N, 3).
    This routine does not differentiate velocity, refit pressure/forcing, or decide PDE
    acceptance. It only decomposes supplied residual values by cylindrical component and
    by four non-overlapping spatial regions.
    """

    time = _validate_scalar("time", time)
    radial_half_width = _validate_scalar("radial_half_width", radial_half_width, positive=True)
    axial_half_height = _validate_scalar("axial_half_height", axial_half_height, positive=True)
    if isinstance(grid_size, bool) or not isinstance(grid_size, (int, np.integer)) or grid_size < 4:
        raise ValueError("grid_size must be an integer >= 4")
    grid_size = int(grid_size)
    core_radius_fraction = _validate_fraction("core_radius_fraction", core_radius_fraction)
    collar_start_fraction = _validate_fraction("collar_start_fraction", collar_start_fraction)
    tip_start_fraction = _validate_fraction("tip_start_fraction", tip_start_fraction)
    if core_radius_fraction >= collar_start_fraction:
        raise ValueError("core_radius_fraction must be smaller than collar_start_fraction")

    x = _cell_centers(radial_half_width, grid_size)
    y = _cell_centers(radial_half_width, grid_size)
    z = _cell_centers(axial_half_height, grid_size)
    xx, yy, zz = np.meshgrid(x, y, z, indexing="ij")
    points = np.column_stack((xx.ravel(), yy.ravel(), zz.ravel()))

    values = np.asarray(residual(points, time), dtype=float)
    if values.shape != points.shape:
        raise ValueError(f"residual must return shape {points.shape}, got {values.shape}")
    if not np.all(np.isfinite(values)):
        raise ValueError("residual returned non-finite values")

    r = np.hypot(points[:, 0], points[:, 1])
    nonaxis = r > 0.0
    erx = np.zeros_like(r)
    ery = np.zeros_like(r)
    etx = np.zeros_like(r)
    ety = np.zeros_like(r)
    erx[nonaxis] = points[nonaxis, 0] / r[nonaxis]
    ery[nonaxis] = points[nonaxis, 1] / r[nonaxis]
    etx[nonaxis] = -points[nonaxis, 1] / r[nonaxis]
    ety[nonaxis] = points[nonaxis, 0] / r[nonaxis]

    radial = values[:, 0] * erx + values[:, 1] * ery
    azimuthal = values[:, 0] * etx + values[:, 1] * ety
    axial = values[:, 2]
    cyl = np.column_stack((radial, azimuthal, axial))
    component_names = ("radial", "azimuthal", "axial")

    sq = np.sum(values * values, axis=1)
    total_mean_sq = float(np.mean(sq))
    if not np.isfinite(total_mean_sq) or total_mean_sq <= np.finfo(float).tiny:
        raise ValueError("residual is numerically zero on the sampled grid")

    component_sq = cyl * cyl
    component_energy = np.mean(component_sq, axis=0)
    component_fraction = component_energy / total_mean_sq

    rho = r / radial_half_width
    zeta = np.abs(points[:, 2]) / axial_half_height
    tip = zeta >= tip_start_fraction
    collar = (~tip) & (rho >= collar_start_fraction)
    core = (~tip) & (~collar) & (rho <= core_radius_fraction)
    mid = ~(tip | collar | core)
    regions = {
        "core": core,
        "mid": mid,
        "radial_collar": collar,
        "axial_tips": tip,
    }

    region_energy_fraction: dict[str, float] = {}
    region_component_energy_fraction: dict[str, dict[str, float]] = {}
    total_sq_sum = float(np.sum(sq))
    component_sq_sum = np.sum(component_sq, axis=0)
    for name, mask in regions.items():
        region_energy_fraction[name] = float(np.sum(sq[mask]) / total_sq_sum)
        per_component: dict[str, float] = {}
        for j, component_name in enumerate(component_names):
            denom = float(component_sq_sum[j])
            per_component[component_name] = (
                float(np.sum(component_sq[mask, j]) / denom)
                if denom > np.finfo(float).tiny
                else 0.0
            )
        region_component_energy_fraction[name] = per_component

    cell_volume = (
        (2.0 * radial_half_width / grid_size)
        * (2.0 * radial_half_width / grid_size)
        * (2.0 * axial_half_height / grid_size)
    )
    sampled_volume = float(points.shape[0] * cell_volume)

    return ResidualLocalizationSnapshot(
        time=time,
        grid_size=grid_size,
        total_rms=float(np.sqrt(total_mean_sq)),
        total_max=float(np.max(np.linalg.norm(values, axis=1))),
        cylindrical_component_rms={
            name: float(np.sqrt(np.mean(component_sq[:, j])))
            for j, name in enumerate(component_names)
        },
        cylindrical_component_max={
            name: float(np.max(np.abs(cyl[:, j])))
            for j, name in enumerate(component_names)
        },
        cylindrical_component_energy_fraction={
            name: float(component_fraction[j])
            for j, name in enumerate(component_names)
        },
        region_energy_fraction=region_energy_fraction,
        region_component_energy_fraction=region_component_energy_fraction,
        region_thresholds={
            "core_radius_fraction": core_radius_fraction,
            "collar_start_fraction": collar_start_fraction,
            "tip_start_fraction": tip_start_fraction,
        },
        sampled_volume=sampled_volume,
    )


def audit_residual_localization_resolution(
    residual: ResidualFn,
    *,
    time: float,
    radial_half_width: float,
    axial_half_height: float,
    grid_sizes: Iterable[int],
    core_radius_fraction: float = 0.35,
    collar_start_fraction: float = 0.75,
    tip_start_fraction: float = 0.75,
) -> ResidualLocalizationResolutionAudit:
    sizes = tuple(grid_sizes)
    if len(sizes) < 3:
        raise ValueError("at least three grid sizes are required")
    if any(
        isinstance(n, bool) or not isinstance(n, (int, np.integer)) or int(n) < 4
        for n in sizes
    ):
        raise ValueError("grid sizes must be integers >= 4")
    sizes = tuple(int(n) for n in sizes)
    if any(b <= a for a, b in zip(sizes, sizes[1:])):
        raise ValueError("grid sizes must be strictly increasing")

    snapshots = tuple(
        diagnose_residual_localization(
            residual,
            time=time,
            radial_half_width=radial_half_width,
            axial_half_height=axial_half_height,
            grid_size=n,
            core_radius_fraction=core_radius_fraction,
            collar_start_fraction=collar_start_fraction,
            tip_start_fraction=tip_start_fraction,
        )
        for n in sizes
    )
    finest = snapshots[-1]
    region_names = tuple(finest.region_energy_fraction)
    component_names = tuple(finest.cylindrical_component_energy_fraction)

    region_deltas = []
    component_deltas = []
    rms_deltas = []
    for snap in snapshots:
        region_deltas.append(
            max(
                abs(
                    snap.region_energy_fraction[name]
                    - finest.region_energy_fraction[name]
                )
                for name in region_names
            )
        )
        component_deltas.append(
            max(
                abs(
                    snap.cylindrical_component_energy_fraction[name]
                    - finest.cylindrical_component_energy_fraction[name]
                )
                for name in component_names
            )
        )
        rms_deltas.append(abs(snap.total_rms - finest.total_rms) / finest.total_rms)

    return ResidualLocalizationResolutionAudit(
        snapshots=snapshots,
        max_region_fraction_delta_to_finest=tuple(float(x) for x in region_deltas),
        max_component_fraction_delta_to_finest=tuple(float(x) for x in component_deltas),
        relative_total_rms_delta_to_finest=tuple(float(x) for x in rms_deltas),
    )
