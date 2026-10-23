"""Independent finite-domain energy concentration diagnostic for axisymmetric fields.

The diagnostic intentionally consumes only a ``velocity(points, time)`` callable and
performs its own Gauss--Legendre quadrature in a meridional half-plane. It assumes
axisymmetry solely for the cylindrical reduction ``dV = 2*pi*r dr dz``; callers must
establish that representation property separately.

The output is a finite-window numerical diagnostic. Stable concentration radii over
several quadrature orders do not imply finite-time blow-up or singularity formation.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Iterable
import math
import numpy as np

from .quadrature import unit_rule

Velocity = Callable[[np.ndarray, float], np.ndarray]


@dataclass(frozen=True)
class EnergyConcentrationRow:
    order: int
    scale: float
    energy: float
    fraction_of_full_support_energy: float


@dataclass(frozen=True)
class EnergyConcentrationReport:
    time: float
    support_radius: float
    support_half_height: float
    levels: tuple[float, ...]
    orders: tuple[int, ...]
    rows: tuple[EnergyConcentrationRow, ...]
    concentration_scales: dict[str, tuple[float, ...]]
    concentration_scale_spread: dict[str, float]
    finest_total_energy: float
    finest_fraction_profile: tuple[float, ...]
    max_adjacent_order_fraction_delta: float
    scope: str


def _finite_positive(name: str, value: float) -> float:
    value = float(value)
    if not math.isfinite(value) or value <= 0.0:
        raise ValueError(f"{name} must be finite and positive")
    return value


def _normalise_levels(levels: Iterable[float]) -> tuple[float, ...]:
    values = tuple(float(v) for v in levels)
    if not values:
        raise ValueError("levels must be nonempty")
    if any((not math.isfinite(v)) or v <= 0.0 or v > 1.0 for v in values):
        raise ValueError("levels must lie in (0, 1]")
    if any(b <= a for a, b in zip(values, values[1:])):
        raise ValueError("levels must be strictly increasing")
    if values[-1] != 1.0:
        values = values + (1.0,)
    return values


def _normalise_orders(orders: Iterable[int]) -> tuple[int, ...]:
    values = tuple(int(v) for v in orders)
    if len(values) < 2:
        raise ValueError("at least two quadrature orders are required")
    if any(v < 2 or v > 2048 for v in values):
        raise ValueError("quadrature orders must lie in [2, 2048]")
    if any(b <= a for a, b in zip(values, values[1:])):
        raise ValueError("quadrature orders must be strictly increasing")
    return values


def _scaled_cylinder_energy(
    velocity: Velocity,
    time: float,
    *,
    support_radius: float,
    support_half_height: float,
    scale: float,
    order: int,
) -> float:
    nodes, weights = unit_rule(order)
    radial_extent = support_radius * scale
    axial_extent = support_half_height * scale
    r = radial_extent * nodes
    z = axial_extent * (2.0 * nodes - 1.0)
    rr, zz = np.meshgrid(r, z, indexing="ij")
    points = np.stack((rr, np.zeros_like(rr), zz), axis=-1)
    values = np.asarray(velocity(points, time), dtype=float)
    if values.shape != points.shape:
        raise ValueError("velocity(points, time) must return an array with points.shape")
    if not np.all(np.isfinite(values)):
        raise ValueError("velocity returned nonfinite values")
    density = 0.5 * np.sum(values * values, axis=-1)
    weighted = weights[:, None] * weights[None, :] * rr * density
    # r = R*s*xi and z = Z*s*(2*eta-1), so dr*dz = 2*R*Z*s^2 dxi deta.
    return float(4.0 * math.pi * support_radius * support_half_height * scale * scale * np.sum(weighted))


def _crossing_scale(levels: tuple[float, ...], fractions: np.ndarray, target: float) -> float:
    if not 0.0 < target < 1.0:
        raise ValueError("target fraction must lie in (0, 1)")
    idx = int(np.searchsorted(fractions, target, side="left"))
    if idx == 0:
        return float(levels[0] * target / fractions[0]) if fractions[0] > 0.0 else float("nan")
    if idx >= len(levels):
        return float("nan")
    f0, f1 = float(fractions[idx - 1]), float(fractions[idx])
    q0, q1 = levels[idx - 1], levels[idx]
    if f1 <= f0:
        return float(q1)
    return float(q0 + (target - f0) * (q1 - q0) / (f1 - f0))


def diagnose_axisymmetric_energy_concentration(
    velocity: Velocity,
    time: float,
    *,
    support_radius: float,
    support_half_height: float,
    levels: Iterable[float] = (0.25, 0.5, 0.75, 1.0),
    orders: Iterable[int] = (16, 32, 64),
    target_fractions: Iterable[float] = (0.5, 0.9),
) -> EnergyConcentrationReport:
    """Measure energy concentration in homothetic cylinders.

    A scale ``q`` integrates over ``0 <= r <= q*support_radius`` and
    ``|z| <= q*support_half_height``. Fractions are always normalised by the
    full-support energy computed at the *same* quadrature order, keeping the
    resolution comparison separate from any training objective.
    """
    if not callable(velocity):
        raise TypeError("velocity must be callable")
    time = float(time)
    if not math.isfinite(time):
        raise ValueError("time must be finite")
    radius = _finite_positive("support_radius", support_radius)
    half_height = _finite_positive("support_half_height", support_half_height)
    levels_t = _normalise_levels(levels)
    orders_t = _normalise_orders(orders)
    targets = tuple(float(v) for v in target_fractions)
    if not targets or any((not math.isfinite(v)) or not 0.0 < v < 1.0 for v in targets):
        raise ValueError("target_fractions must be finite values in (0, 1)")
    if len(set(targets)) != len(targets):
        raise ValueError("target_fractions must be unique")

    rows: list[EnergyConcentrationRow] = []
    per_order_profiles: list[np.ndarray] = []
    scales_by_target: dict[str, list[float]] = {f"q{int(round(100*t))}": [] for t in targets}
    totals: list[float] = []
    for order in orders_t:
        energies = np.array([
            _scaled_cylinder_energy(
                velocity,
                time,
                support_radius=radius,
                support_half_height=half_height,
                scale=q,
                order=order,
            )
            for q in levels_t
        ], dtype=float)
        total = float(energies[-1])
        if not math.isfinite(total) or total <= 0.0:
            raise ValueError("full-support kinetic energy must be finite and positive")
        fractions = energies / total
        if np.any(np.diff(fractions) < -1e-12):
            raise ValueError("energy fractions are not numerically monotone")
        fractions = np.maximum.accumulate(fractions)
        fractions[-1] = 1.0
        totals.append(total)
        per_order_profiles.append(fractions)
        for q, energy, fraction in zip(levels_t, energies, fractions):
            rows.append(EnergyConcentrationRow(order, q, float(energy), float(fraction)))
        for target in targets:
            scales_by_target[f"q{int(round(100*target))}"].append(
                _crossing_scale(levels_t, fractions, target)
            )

    adjacent = [
        float(np.max(np.abs(per_order_profiles[i + 1] - per_order_profiles[i])))
        for i in range(len(per_order_profiles) - 1)
    ]
    scales_tuple = {key: tuple(values) for key, values in scales_by_target.items()}
    spreads = {
        key: float(np.max(values) - np.min(values)) if np.all(np.isfinite(values)) else float("nan")
        for key, values in scales_tuple.items()
    }
    return EnergyConcentrationReport(
        time=time,
        support_radius=radius,
        support_half_height=half_height,
        levels=levels_t,
        orders=orders_t,
        rows=tuple(rows),
        concentration_scales=scales_tuple,
        concentration_scale_spread=spreads,
        finest_total_energy=totals[-1],
        finest_fraction_profile=tuple(float(v) for v in per_order_profiles[-1]),
        max_adjacent_order_fraction_delta=max(adjacent),
        scope=(
            "finite-window axisymmetric kinetic-energy concentration diagnostic; "
            "quadrature stability or decreasing concentration scale is not evidence of blow-up"
        ),
    )
