"""Independent checks for ``SmoothExtendedAxial``.

The checks use split Gauss--Legendre rules at the saved geometry and bump
supports.  They include selected Chebyshev nodes and independent off-grid Z
values, core preservation, exact zero primitive tail, the same-grid eta
derivative, and quadrature refinement.  Passing these checks only validates
the finite numerical adapter and its stated gauge.
"""

from __future__ import annotations

import json
from pathlib import Path
import sys
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
LOCAL = ROOT / "experiments" / "root_st073"
for path in (SRC, LOCAL):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from lei_ren_part1_smooth_extended_axial import (  # noqa: E402
    SmoothExtendedAxial,
)


NODE_FRACTIONS = (0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0)
OFFGRID_Z = (-0.975, -0.91, -0.77, -0.61, -0.43, -0.17,
             0.06, 0.31, 0.48, 0.73, 0.93, 0.975)
MOMENT_ORDERS = (96, 128)
MOMENT_THRESHOLDS = {
    "mass_max_abs": 2.0e-10,
    "mixed_max_abs": 1.0e-7,
    "quadratic_max_abs": 1.0e-7,
    "quadrature_refinement_max_abs": 2.0e-5,
}
CORE_THRESHOLD = 2.0e-12
TAIL_THRESHOLD = 0.0
ETA_DERIVATIVE_STEP = 1.0e-6
ETA_DERIVATIVE_THRESHOLD = 2.0e-5
AVERAGE_THRESHOLD = 2.0e-8


def _integral_moments(matched: SmoothExtendedAxial, z: float, order: int) -> dict[str, float]:
    radii, weights = matched._radial_rule(float(z), order)
    f_values = np.asarray([matched.F(float(radius), z) for radius in radii], dtype=float)
    u_values = np.asarray([matched.U(float(radius), z) for radius in radii], dtype=float)
    return {
        "mass": float(weights @ u_values),
        "mixed": float(weights @ (2.0 * radii * f_values * u_values)),
        "quadratic": float(weights @ (u_values * u_values - radii * f_values * f_values)),
    }


def _primitive_integral(matched: SmoothExtendedAxial, radius: float, z: float, order: int = 128) -> float:
    if radius == 0.0:
        return 0.0
    clipped = min(float(radius), matched.Rmax)
    # Reuse the same split geometry by integrating the public U on a truncated
    # rule assembled locally.  This is an independent check of average_U.
    nodes, weights = np.polynomial.legendre.leggauss(order)
    edges = [0.0, clipped, matched.R_core, matched.R_cutoff]
    edges.extend(edge for pair in matched.bump_supports for edge in pair)
    edges = sorted({value for value in edges if 0.0 <= value <= clipped})
    total = 0.0
    for left, right in zip(edges[:-1], edges[1:]):
        if right <= left:
            continue
        width = right - left
        points = left + width * (nodes + 1.0) / 2.0
        radial_weights = width * weights / 2.0
        total += float(
            radial_weights
            @ np.asarray([matched.U(float(point), z) for point in points], dtype=float)
        )
    return float(total / clipped)


def _max_abs(values: list[float]) -> float:
    return float(max((abs(float(value)) for value in values), default=0.0))


def run() -> dict[str, Any]:
    matched = SmoothExtendedAxial()
    node_indices = tuple(
        int(round(fraction * (matched.z_nodes - 1)))
        for fraction in NODE_FRACTIONS
    )
    node_z = [float(matched.z_grid.eta[index]) for index in node_indices]
    all_z = tuple(node_z) + tuple(float(value) for value in OFFGRID_Z)
    moment_rows: list[dict[str, Any]] = []
    residuals = {"mass": [], "mixed": [], "quadratic": [], "refinement": []}

    for z, location in [(value, "node") for value in node_z] + [
        (float(value), "offgrid") for value in OFFGRID_Z
    ]:
        row: dict[str, Any] = {"Z": z, "location": location, "orders": {}}
        target = float(matched.quadratic_target(z))
        for order in MOMENT_ORDERS:
            moments = _integral_moments(matched, z, order)
            residual = {
                "mass": moments["mass"],
                "mixed": moments["mixed"],
                "quadratic": moments["quadratic"] - target,
            }
            row["orders"][str(order)] = {
                "moments": moments,
                "target": target,
                "residual": residual,
            }
        refinement = {
            key: row["orders"]["128"]["moments"][key]
            - row["orders"]["96"]["moments"][key]
            for key in ("mass", "mixed", "quadratic")
        }
        row["refinement_128_minus_96"] = refinement
        independent = row["orders"]["128"]["residual"]
        for key in ("mass", "mixed", "quadratic"):
            residuals[key].append(float(independent[key]))
        residuals["refinement"].extend(float(value) for value in refinement.values())
        moment_rows.append(row)

    moment_max = {key: _max_abs(values) for key, values in residuals.items()}

    # The near-axis branch delegates U and U_eta to the shared 257-node core.
    core_errors: list[float] = []
    for z in OFFGRID_Z[:6]:
        for radius in (0.0, 0.001, 0.004, 0.0049):
            core_errors.extend(
                (
                    matched.U(radius, z) - float(matched.core.U(radius, z)),
                    matched.dU_deta(radius, z) - float(matched.core.dU_deta(radius, z)),
                    matched.average_U(radius, z)
                    - float(matched.core.radial_average_U(radius, z)),
                    matched.average_dU_deta(radius, z)
                    - float(matched.core.radial_average_dU_deta(radius, z)),
                )
            )
    core_max = _max_abs(core_errors)

    # The primitive and its eta derivative are explicitly zero at and beyond Ra.
    tail_errors: list[float] = []
    for z in OFFGRID_Z[:6]:
        for radius in (matched.Rmax, matched.Rmax * 1.01, matched.Rmax * 2.0):
            tail_errors.extend(
                (
                    matched.U(radius, z),
                    matched.dU_deta(radius, z),
                    matched.average_U(radius, z),
                    matched.average_dU_deta(radius, z),
                )
            )
    tail_max = _max_abs(tail_errors)

    # Verify the public primitive independently at representative radii.
    average_errors: list[float] = []
    for z in OFFGRID_Z[::2]:
        for radius in (0.01, 0.015, 0.05, 0.16 * matched.Rmax, 0.52 * matched.Rmax,
                       0.93 * matched.Rmax):
            average_errors.append(
                matched.average_U(radius, z)
                - _primitive_integral(matched, radius, z)
            )
    average_max = _max_abs(average_errors)

    # Same-grid eta interpolation derivative, checked against an independent
    # scalar central difference away from the endpoints and flat joins.
    derivative_errors: list[float] = []
    derivative_rows: list[dict[str, float]] = []
    derivative_radii = (
        0.012,
        sum(matched.bump_supports[0]) / 2.0,
        sum(matched.bump_supports[1]) / 2.0,
        sum(matched.bump_supports[2]) / 2.0,
        sum(matched.bump_supports[3]) / 2.0,
    )
    for z in OFFGRID_Z[1:-1:2]:
        for radius in derivative_radii:
            h = ETA_DERIVATIVE_STEP
            finite_difference = (
                matched.U(radius, z + h) - matched.U(radius, z - h)
            ) / (2.0 * h)
            analytic = matched.dU_deta(radius, z)
            error = float(finite_difference - analytic)
            derivative_errors.append(abs(error))
            derivative_rows.append(
                {
                    "R": float(radius),
                    "Z": float(z),
                    "finite_difference": float(finite_difference),
                    "analytic": float(analytic),
                    "error": error,
                }
            )
    derivative_max = _max_abs(derivative_errors)

    average_derivative_errors: list[float] = []
    for z in OFFGRID_Z[1:-1:2]:
        for radius in (0.01, 0.16 * matched.Rmax, 0.52 * matched.Rmax,
                       0.93 * matched.Rmax):
            h = ETA_DERIVATIVE_STEP
            finite_difference = (
                matched.average_U(radius, z + h)
                - matched.average_U(radius, z - h)
            ) / (2.0 * h)
            average_derivative_errors.append(
                abs(float(finite_difference - matched.average_dU_deta(radius, z)))
            )
    average_derivative_max = _max_abs(average_derivative_errors)

    metadata = matched.metadata()
    checks_passed = bool(
        moment_max["mass"] <= MOMENT_THRESHOLDS["mass_max_abs"]
        and moment_max["mixed"] <= MOMENT_THRESHOLDS["mixed_max_abs"]
        and moment_max["quadratic"] <= MOMENT_THRESHOLDS["quadratic_max_abs"]
        and moment_max["refinement"] <= MOMENT_THRESHOLDS["quadrature_refinement_max_abs"]
        and core_max <= CORE_THRESHOLD
        and tail_max <= TAIL_THRESHOLD
        and average_max <= AVERAGE_THRESHOLD
        and derivative_max <= ETA_DERIVATIVE_THRESHOLD
        and average_derivative_max <= ETA_DERIVATIVE_THRESHOLD
        and metadata["branch_orientation"]["all_positive_at_nodes"]
    )

    # Keep a compact replay receipt so a later process can inspect the saved
    # grid and coefficients without parsing the implementation internals.
    replay = {
        "z_grid": np.asarray(matched.z_grid.eta, dtype=float).tolist(),
        "coefficient_grid": np.asarray(matched._coefficient_grid, dtype=float).tolist(),
        "coefficient_eta_grid": np.asarray(matched._coefficient_eta_grid, dtype=float).tolist(),
        "quadratic_target_grid": np.asarray(matched._quadratic_target_grid, dtype=float).tolist(),
        "base_mass_grid": np.asarray(matched._base_mass_grid, dtype=float).tolist(),
        "bump_mixed_grid": np.asarray(matched._bump_mixed_grid, dtype=float).tolist(),
    }
    report: dict[str, Any] = {
        "status": "completed" if checks_passed else "failed",
        "checks_passed": checks_passed,
        "source": metadata["source"],
        "source_version": metadata["source_version"],
        "profile": "SmoothExtendedAxial(load_extended_profile())",
        "metadata": metadata,
        "thresholds": {
            **MOMENT_THRESHOLDS,
            "core_max_abs": CORE_THRESHOLD,
            "tail_max_abs": TAIL_THRESHOLD,
            "average_primitive_max_abs": AVERAGE_THRESHOLD,
            "eta_derivative_max_abs": ETA_DERIVATIVE_THRESHOLD,
        },
        "checked_node_indices": list(node_indices),
        "checked_node_Z": node_z,
        "offgrid_Z": list(OFFGRID_Z),
        "moment_rows": moment_rows,
        "moment_max_abs": moment_max,
        "core_preservation_max_abs": core_max,
        "zero_tail_max_abs": tail_max,
        "average_primitive_max_abs": average_max,
        "eta_derivative_fd_step": ETA_DERIVATIVE_STEP,
        "eta_derivative_fd_max_abs": derivative_max,
        "average_eta_derivative_fd_max_abs": average_derivative_max,
        "eta_derivative_rows": derivative_rows,
        "replay": replay,
        "five_moment_closure": False,
        "pde_validated": False,
        "whole_cone_validated": False,
        "scope": (
            "Finite smooth Z-dependent axial repair of the shared-pressure "
            "ExtendedSwirl candidate. Independent off-grid moment checks and "
            "primitive/eta-derivative checks are numerical only; no complete "
            "Part I field, PDE, stress cone, global energy, or five-moment "
            "certificate is claimed."
        ),
    }
    destination = Path(__file__).with_suffix(".json")
    destination.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "status": report["status"],
                "checks_passed": report["checks_passed"],
                "moment_max_abs": moment_max,
                "core_max_abs": core_max,
                "tail_max_abs": tail_max,
                "average_max_abs": average_max,
                "eta_derivative_max_abs": derivative_max,
                "average_eta_derivative_max_abs": average_derivative_max,
            }
        ),
        flush=True,
    )
    if not checks_passed:
        raise AssertionError(
            "smooth extended axial checks failed; see "
            f"{destination.name} for the full receipt"
        )
    return report


if __name__ == "__main__":
    run()
