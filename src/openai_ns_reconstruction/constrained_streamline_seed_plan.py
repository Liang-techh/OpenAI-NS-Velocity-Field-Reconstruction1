"""Truth-bounded, activity-aware streamline seed planning.

This module is visualization support only. It does not alter the candidate,
pressure, forcing, residual, validation thresholds, or any scientific claim.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np
from scipy.spatial.distance import pdist
from scipy.stats import qmc

VelocityCallable = Callable[[np.ndarray, float], np.ndarray]


@dataclass(frozen=True)
class StreamlineSeedPlan:
    seeds: np.ndarray
    speeds: np.ndarray
    visualization_seed: int
    requested_count: int
    poisson_radius_unit: float
    min_normalized_separation: float
    activity_floor: float
    support_radius: float
    support_half_height: float
    support_margin_fraction: float
    claim_scope: str = "visualization_seed_selection_only"
    visualization_ready: bool = False
    visual_correspondence_verified: bool = False
    pde_validated: bool = False
    paper_exact: bool = False
    openai_field_identified: bool = False
    blowup_proved: bool = False


def _velocity_values(velocity: VelocityCallable, points: np.ndarray, time: float) -> np.ndarray:
    values = np.asarray(velocity(points, float(time)), dtype=float)
    if values.shape != points.shape:
        raise ValueError("velocity must return shape (n, 3)")
    if not np.all(np.isfinite(values)):
        raise ValueError("velocity returned non-finite values")
    return values


def _poisson_engine(radius: float, seed: int) -> qmc.PoissonDisk:
    # SciPy >=1.15 prefers rng=; the repository still supports SciPy 1.10,
    # where seed= is the compatible public spelling.
    try:
        return qmc.PoissonDisk(
            d=3,
            radius=radius,
            optimization=None,
            rng=np.random.default_rng(seed),
        )
    except TypeError:
        return qmc.PoissonDisk(
            d=3,
            radius=radius,
            optimization=None,
            seed=seed,
        )


def plan_activity_aware_streamline_seeds(
    velocity: VelocityCallable,
    time: float,
    *,
    count: int,
    support_radius: float,
    support_half_height: float,
    visualization_seed: int = 914113,
    support_margin_fraction: float = 0.95,
    relative_activity_floor: float = 0.02,
    pool_factor: int = 5,
    max_attempts: int = 8,
    initial_poisson_radius_unit: float | None = None,
    min_peak_speed: float = 1.0e-12,
) -> StreamlineSeedPlan:
    """Return exactly ``count`` well-spaced active streamline seeds.

    Poisson-disk sampling is performed in a normalized cube. Points are then
    restricted to a normalized cylinder and mapped into the caller-declared
    physical support. Among the well-spaced candidate points, seeds below a
    declared relative speed floor are rejected; the fastest remaining points
    are retained. This is an autonomous visualization choice, not a PDE or
    OpenAI-source acceptance rule.
    """

    if not isinstance(count, (int, np.integer)) or count < 1:
        raise ValueError("count must be a positive integer")
    if not np.isfinite(time):
        raise ValueError("time must be finite")
    if not np.isfinite(support_radius) or support_radius <= 0:
        raise ValueError("support_radius must be positive and finite")
    if not np.isfinite(support_half_height) or support_half_height <= 0:
        raise ValueError("support_half_height must be positive and finite")
    if not isinstance(visualization_seed, (int, np.integer)):
        raise ValueError("visualization_seed must be an integer")
    if not (0.0 < support_margin_fraction < 1.0):
        raise ValueError("support_margin_fraction must lie in (0, 1)")
    if not (0.0 <= relative_activity_floor < 1.0):
        raise ValueError("relative_activity_floor must lie in [0, 1)")
    if not isinstance(pool_factor, (int, np.integer)) or pool_factor < 2:
        raise ValueError("pool_factor must be an integer >= 2")
    if not isinstance(max_attempts, (int, np.integer)) or max_attempts < 1:
        raise ValueError("max_attempts must be a positive integer")
    if not np.isfinite(min_peak_speed) or min_peak_speed <= 0:
        raise ValueError("min_peak_speed must be positive and finite")

    if initial_poisson_radius_unit is None:
        radius0 = 0.12 * (200.0 / float(count)) ** (1.0 / 3.0)
        radius0 = float(np.clip(radius0, 0.025, 0.22))
    else:
        radius0 = float(initial_poisson_radius_unit)
        if not np.isfinite(radius0) or not (0.0 < radius0 < 0.5):
            raise ValueError("initial_poisson_radius_unit must lie in (0, 0.5)")

    pool_target = max(int(count) * int(pool_factor), int(count) + 64)
    last_active = 0
    last_peak = 0.0

    for attempt in range(int(max_attempts)):
        radius = radius0 * (0.85 ** attempt)
        engine = _poisson_engine(radius, int(visualization_seed) + attempt)
        unit = np.asarray(engine.random(pool_target), dtype=float)
        if unit.ndim != 2 or unit.shape[1] != 3 or not np.all(np.isfinite(unit)):
            raise RuntimeError("PoissonDisk returned malformed samples")

        ux = 2.0 * unit[:, 0] - 1.0
        uy = 2.0 * unit[:, 1] - 1.0
        uz = 2.0 * unit[:, 2] - 1.0
        inside = ux * ux + uy * uy <= 1.0
        normalized = np.column_stack((ux[inside], uy[inside], uz[inside]))
        if normalized.shape[0] < count:
            continue

        margin = float(support_margin_fraction)
        points = np.empty_like(normalized)
        points[:, 0] = normalized[:, 0] * support_radius * margin
        points[:, 1] = normalized[:, 1] * support_radius * margin
        points[:, 2] = normalized[:, 2] * support_half_height * margin

        values = _velocity_values(velocity, points, float(time))
        speeds = np.linalg.norm(values, axis=1)
        peak = float(np.max(speeds)) if speeds.size else 0.0
        last_peak = peak
        if not np.isfinite(peak) or peak < min_peak_speed:
            raise ValueError("velocity is numerically inactive at the seed-planning probes")

        activity_floor = max(min_peak_speed, relative_activity_floor * peak)
        active_indices = np.flatnonzero(speeds >= activity_floor)
        last_active = int(active_indices.size)
        if active_indices.size < count:
            continue

        # Stable descending speed order keeps the choice deterministic while
        # preserving the Poisson-disk separation inherited from the pool.
        local_order = np.argsort(-speeds[active_indices], kind="stable")
        chosen = active_indices[local_order[:count]]
        selected_points = np.array(points[chosen], copy=True)
        selected_speeds = np.array(speeds[chosen], copy=True)

        normalized_selected = np.empty_like(selected_points)
        normalized_selected[:, 0] = selected_points[:, 0] / (support_radius * margin)
        normalized_selected[:, 1] = selected_points[:, 1] / (support_radius * margin)
        normalized_selected[:, 2] = selected_points[:, 2] / (support_half_height * margin)
        min_sep = float(np.min(pdist(normalized_selected))) if count > 1 else float("inf")

        selected_points.setflags(write=False)
        selected_speeds.setflags(write=False)
        return StreamlineSeedPlan(
            seeds=selected_points,
            speeds=selected_speeds,
            visualization_seed=int(visualization_seed),
            requested_count=int(count),
            poisson_radius_unit=float(radius),
            min_normalized_separation=min_sep,
            activity_floor=float(activity_floor),
            support_radius=float(support_radius),
            support_half_height=float(support_half_height),
            support_margin_fraction=margin,
        )

    raise ValueError(
        "could not obtain the requested active seed count under the declared "
        f"settings; last_active={last_active}, last_peak_speed={last_peak:.6g}"
    )
