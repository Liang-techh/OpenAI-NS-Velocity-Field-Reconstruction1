"""Finite-domain kinetic-energy truncation diagnostics for a frozen velocity field.

The diagnostic consumes only a ``velocity(points, time) -> (..., 3)`` callable.
It integrates kinetic energy over nested Cartesian cubes with independent
Gauss--Legendre quadrature orders.  The outermost cube is only a finite
reference domain: energy not seen beyond it is explicitly *unresolved*, so a
small outer-shell fraction is not promoted to a rigorous tail bound.

This is geometry/visualization-support evidence only.  It does not evaluate a
Navier--Stokes residual and does not imply PDE validation, visual
correspondence, paper exactness, or blow-up.
"""
from __future__ import annotations

from typing import Any, Callable, Iterable, Sequence

import numpy as np
from numpy.polynomial.legendre import leggauss


VelocityFn = Callable[[np.ndarray, float], np.ndarray]


def _validate_levels(name: str, values: Iterable[float], *, integer: bool) -> tuple:
    raw = tuple(values)
    if len(raw) < 3:
        raise ValueError(f"{name} must contain at least three levels")
    if integer:
        numeric = tuple(float(v) for v in raw)
        if any(
            (not np.isfinite(v)) or v < 2.0 or v != float(int(v))
            for v in numeric
        ):
            raise ValueError(f"{name} entries must be finite integers >= 2")
        result = tuple(int(v) for v in numeric)
    else:
        result = tuple(float(v) for v in raw)
        if any((not np.isfinite(v)) or v <= 0.0 for v in result):
            raise ValueError(f"{name} entries must be positive and finite")
    if any(a >= b for a, b in zip(result, result[1:])):
        raise ValueError(f"{name} must be strictly increasing")
    return result


def _evaluate_velocity(
    velocity: VelocityFn,
    points: np.ndarray,
    time: float,
) -> np.ndarray:
    raw = np.asarray(velocity(points, float(time)), dtype=float)
    if raw.shape != points.shape:
        raise ValueError("velocity must return shape (n, 3)")
    if not np.all(np.isfinite(raw)):
        raise ValueError("velocity returned non-finite values")
    return raw


def integrate_cube_kinetic_energy(
    velocity: VelocityFn,
    time: float,
    half_width: float,
    quadrature_order: int,
    *,
    chunk_size: int = 65536,
) -> float:
    """Return ``1/2 ∫_{[-L,L]^3} |u|^2 dx`` by tensor Gauss quadrature."""
    if not np.isfinite(time):
        raise ValueError("time must be finite")
    half_width = float(half_width)
    if not np.isfinite(half_width) or half_width <= 0.0:
        raise ValueError("half_width must be positive and finite")
    if int(quadrature_order) != quadrature_order or quadrature_order < 2:
        raise ValueError("quadrature_order must be an integer >= 2")
    if int(chunk_size) != chunk_size or chunk_size <= 0:
        raise ValueError("chunk_size must be a positive integer")

    order = int(quadrature_order)
    nodes, weights = leggauss(order)
    coords = half_width * nodes
    scaled_weights = half_width * weights

    x, y, z = np.meshgrid(coords, coords, coords, indexing="ij")
    wx, wy, wz = np.meshgrid(
        scaled_weights, scaled_weights, scaled_weights, indexing="ij"
    )
    points = np.stack((x.ravel(), y.ravel(), z.ravel()), axis=-1)
    volume_weights = (wx * wy * wz).ravel()

    total = 0.0
    for start in range(0, points.shape[0], int(chunk_size)):
        stop = min(start + int(chunk_size), points.shape[0])
        velocity_values = _evaluate_velocity(velocity, points[start:stop], float(time))
        speed_sq = np.einsum("ij,ij->i", velocity_values, velocity_values)
        total += 0.5 * float(np.dot(volume_weights[start:stop], speed_sq))

    if not np.isfinite(total) or total < 0.0:
        raise ValueError("kinetic-energy quadrature produced an invalid result")
    return total


def diagnose_finite_domain_energy_truncation(
    velocity: VelocityFn,
    time: float,
    *,
    half_widths: Sequence[float] = (1.0, 1.5, 2.0),
    quadrature_orders: Sequence[int] = (12, 24, 48),
    chunk_size: int = 65536,
) -> dict[str, Any]:
    """Audit finite-window energy capture across domain size and quadrature.

    The result intentionally distinguishes two effects:
    - ``quadrature_sensitivity``: integration error at one fixed cube size;
    - ``domain_shell_fractions``: energy newly exposed by enlarging the cube.

    The largest cube is not treated as all of R^3.  Therefore
    ``tail_beyond_largest_cube_resolved`` is always false unless another
    analytic argument outside this diagnostic supplies such a bound.
    """
    if not np.isfinite(time):
        raise ValueError("time must be finite")
    widths = _validate_levels("half_widths", half_widths, integer=False)
    orders = _validate_levels("quadrature_orders", quadrature_orders, integer=True)

    energy_rows: list[dict[str, Any]] = []
    energies_by_order: dict[int, np.ndarray] = {}

    for order in orders:
        energies = np.asarray(
            [
                integrate_cube_kinetic_energy(
                    velocity,
                    float(time),
                    width,
                    order,
                    chunk_size=chunk_size,
                )
                for width in widths
            ],
            dtype=float,
        )
        energies_by_order[order] = energies
        reference = float(energies[-1])
        if reference <= np.finfo(float).tiny:
            raise ValueError("outermost finite-domain energy is numerically zero")

        captured = energies / reference
        shell = np.empty_like(energies)
        shell[0] = energies[0] / reference
        shell[1:] = np.diff(energies) / reference
        scale = max(1.0, reference)
        monotonic = bool(np.all(np.diff(energies) >= -1e-12 * scale))

        energy_rows.append(
            {
                "quadrature_order": int(order),
                "energies": [float(v) for v in energies],
                "captured_fractions_of_largest_cube": [float(v) for v in captured],
                "domain_shell_fractions_of_largest_cube": [float(v) for v in shell],
                "outermost_added_shell_fraction": float(shell[-1]),
                "energy_monotone_with_domain": monotonic,
            }
        )

    finest = energies_by_order[orders[-1]]
    quadrature_sensitivity = []
    for width_index, width in enumerate(widths):
        finest_value = float(finest[width_index])
        entries = []
        for order in orders:
            value = float(energies_by_order[order][width_index])
            absolute = abs(value - finest_value)
            relative = None if finest_value <= np.finfo(float).tiny else absolute / finest_value
            entries.append(
                {
                    "quadrature_order": int(order),
                    "energy": value,
                    "absolute_delta_to_finest": float(absolute),
                    "relative_delta_to_finest": None if relative is None else float(relative),
                }
            )
        quadrature_sensitivity.append(
            {
                "half_width": float(width),
                "entries": entries,
            }
        )

    finest_reference = float(finest[-1])
    finest_captured = finest / finest_reference
    finest_shell = np.empty_like(finest)
    finest_shell[0] = finest[0] / finest_reference
    finest_shell[1:] = np.diff(finest) / finest_reference

    return {
        "claim_scope": "finite_domain_energy_truncation_observation_only",
        "time": float(time),
        "half_widths": [float(v) for v in widths],
        "quadrature_orders": [int(v) for v in orders],
        "energy_definition": "0.5 * integral_[cube] |u|^2 dx",
        "energy_rows": energy_rows,
        "finest_order": int(orders[-1]),
        "finest_energies": [float(v) for v in finest],
        "finest_captured_fractions_of_largest_cube": [
            float(v) for v in finest_captured
        ],
        "finest_domain_shell_fractions_of_largest_cube": [
            float(v) for v in finest_shell
        ],
        "finest_outermost_added_shell_fraction": float(finest_shell[-1]),
        "quadrature_sensitivity": quadrature_sensitivity,
        "tail_beyond_largest_cube_resolved": False,
        "tail_boundary_statement": (
            "The largest cube is a finite reference window only. Energy outside "
            "that cube is unresolved here and requires a larger-domain sweep or "
            "an independent analytic support/tail bound."
        ),
        "visualization_ready": False,
        "pde_validated": False,
        "paper_exact": False,
        "blowup_proved": False,
    }
