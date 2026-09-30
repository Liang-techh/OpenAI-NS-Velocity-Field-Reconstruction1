"""Independent checks for the finite Part I axial-moment repair.

The moment integrals below use fresh Gauss--Legendre rules of orders 96 and
128, split at every cutoff and bump edge.  They do not call the object's
cached radial primitives.  The receipt intentionally does not claim the five
moment closure, stress cone, PDE, energy, forcing, or recursive construction.
"""

from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
LOCAL = Path(__file__).resolve().parent
for path in (SRC, LOCAL):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from lei_ren_part1_axial_match import (  # noqa: E402
    AxialMatchedProfile,
    BUMP_SUPPORTS,
    R_CORE,
    R_JOIN,
    SOURCE,
    SOURCE_VERSION,
    load_core,
)


Z_OFFGRID = (-0.43, -0.17, 0.06, 0.31, 0.48)
RADIAL_EDGES = (0.0, R_CORE, BUMP_SUPPORTS[0][0], BUMP_SUPPORTS[0][1],
                BUMP_SUPPORTS[1][0], BUMP_SUPPORTS[1][1], R_JOIN)


def _integral(function, *, order: int) -> float:
    nodes, weights = np.polynomial.legendre.leggauss(order)
    total = 0.0
    for left, right in zip(RADIAL_EDGES[:-1], RADIAL_EDGES[1:]):
        midpoint = (left + right) / 2.0
        half = (right - left) / 2.0
        points = midpoint + half * nodes
        values = np.asarray([function(float(point)) for point in points], dtype=float)
        total += half * float(weights @ values)
    return float(total)


def _max_abs(values) -> float:
    return float(np.max(np.abs(np.asarray(values, dtype=float))))


def run() -> dict[str, object]:
    matcher = AxialMatchedProfile(*load_core())
    metadata = matcher.metadata()
    conditioning = metadata["matrix_conditioning"]
    thresholds = {
        "conditioning_relative_min_lower": 1.0e-8,
        "nodal_axial_moment_max": 5.0e-13,
        "nodal_mixed_moment_max": 5.0e-13,
        "offgrid_axial_moment_max": 2.0e-10,
        "offgrid_mixed_moment_max": 2.0e-8,
        "eta_derivative_fd_max": 2.0e-5,
        "core_unchanged_max": 2.0e-13,
    }
    if conditioning["relative_min_abs"] <= thresholds["conditioning_relative_min_lower"]:
        raise AssertionError("axial moment matrix conditioning is below the declared bound")
    if metadata["nodal_moment_residuals"]["max_abs_integral_U"] > thresholds["nodal_axial_moment_max"]:
        raise AssertionError("nodal axial moment identity failed")
    if metadata["nodal_moment_residuals"]["max_abs_integral_2RFU"] > thresholds["nodal_mixed_moment_max"]:
        raise AssertionError("nodal mixed moment identity failed")

    independent = []
    for eta in Z_OFFGRID:
        row = {"Z": eta, "orders": {}}
        for order in (96, 128):
            m0 = _integral(lambda radius: matcher.U(radius, eta), order=order)
            m1 = _integral(
                lambda radius: 2.0 * radius * matcher.F(radius, eta) * matcher.U(radius, eta),
                order=order,
            )
            row["orders"][str(order)] = {
                "integral_U": m0,
                "integral_2RFU": m1,
            }
        row["refinement_128_minus_96"] = {
            "integral_U": row["orders"]["128"]["integral_U"] - row["orders"]["96"]["integral_U"],
            "integral_2RFU": row["orders"]["128"]["integral_2RFU"] - row["orders"]["96"]["integral_2RFU"],
        }
        independent.append(row)

    offgrid_m0 = [row["orders"]["128"]["integral_U"] for row in independent]
    offgrid_m1 = [row["orders"]["128"]["integral_2RFU"] for row in independent]
    offgrid_m0_max = _max_abs(offgrid_m0)
    offgrid_m1_max = _max_abs(offgrid_m1)
    if offgrid_m0_max > thresholds["offgrid_axial_moment_max"]:
        raise AssertionError("off-grid axial moment identity failed")
    if offgrid_m1_max > thresholds["offgrid_mixed_moment_max"]:
        raise AssertionError("off-grid mixed moment identity failed")

    # Verify the repaired axial profile keeps the saved core exactly through
    # R=.05 and is identically zero, including radial averages, outside .2.
    branch_z = np.asarray((*Z_OFFGRID, 0.0), dtype=float)
    core, joined = load_core()
    core_errors = []
    for radius in (0.0, 0.01, 0.049):
        for eta in branch_z:
            core_errors.extend((
                matcher.U(radius, eta) - float(core.U(radius, eta)),
                matcher.dU_deta(radius, eta) - float(core.dU_deta(radius, eta)),
                matcher.average_U(radius, eta) - float(core.radial_average_U(radius, eta)),
                matcher.average_dU_deta(radius, eta) - float(core.radial_average_dU_deta(radius, eta)),
            ))
    core_unchanged_max = _max_abs(core_errors)
    if core_unchanged_max > thresholds["core_unchanged_max"]:
        raise AssertionError("saved core branch changed below R=.05")

    exterior_rows = []
    for eta in branch_z:
        exterior_rows.append({
            "Z": float(eta),
            "U_at_Rjoin": matcher.U(R_JOIN, float(eta)),
            "average_U_at_Rjoin": matcher.average_U(R_JOIN, float(eta)),
            "average_dU_deta_at_Rjoin": matcher.average_dU_deta(R_JOIN, float(eta)),
            "U_exterior": matcher.U(0.25, float(eta)),
            "average_U_exterior": matcher.average_U(0.25, float(eta)),
            "average_dU_deta_exterior": matcher.average_dU_deta(0.25, float(eta)),
            "F_exterior_minus_joined": matcher.F(0.25, float(eta)) - float(joined.F(0.25, float(eta))),
        })
    exterior_max = _max_abs([
        value
        for row in exterior_rows
        for key, value in row.items()
        if key != "Z"
    ])
    if exterior_max != 0.0:
        raise AssertionError("the exterior axial branch or F handoff is not exact")

    # The eta derivative uses the cached Chebyshev derivative for c1 and the
    # analytic identity c2=-B-c1.  Compare it against an independent central
    # difference of the scalar U evaluator away from the eta endpoints.
    derivative_rows = []
    derivative_errors = []
    step = 1.0e-6
    for radius in (0.025, 0.075, 0.14, 0.19, 0.25):
        for eta in Z_OFFGRID:
            finite_difference = (matcher.U(radius, eta + step) - matcher.U(radius, eta - step)) / (2.0 * step)
            error = finite_difference - matcher.dU_deta(radius, eta)
            derivative_errors.append(abs(float(error)))
            derivative_rows.append({
                "R": radius,
                "Z": eta,
                "finite_difference": float(finite_difference),
                "analytic_interpolated": float(matcher.dU_deta(radius, eta)),
                "error": float(error),
            })
    derivative_max = float(max(derivative_errors))
    if derivative_max > thresholds["eta_derivative_fd_max"]:
        raise AssertionError("eta derivative check failed")

    report = {
        "status": "completed",
        "candidate_id": "lr1-axial-moment-match-001",
        "source": SOURCE,
        "source_version": SOURCE_VERSION,
        "source_section": "Section 2.5, equations (2.21)-(2.22)",
        "profile": "AxialMatchedProfile(load_core())",
        "parameters": {
            "R_core": R_CORE,
            "R_join": R_JOIN,
            "bump_supports": [list(item) for item in BUMP_SUPPORTS],
            "quadrature_order_constructor": matcher.quadrature_order,
            "independent_quadrature_orders": [96, 128],
            "offgrid_Z": list(Z_OFFGRID),
        },
        "thresholds": thresholds,
        "matrix_conditioning": conditioning,
        "nodal_moment_residuals": metadata["nodal_moment_residuals"],
        "independent_offgrid_moments": independent,
        "independent_offgrid_max_abs_order128": {
            "integral_U": offgrid_m0_max,
            "integral_2RFU": offgrid_m1_max,
        },
        "core_unchanged_max_abs": core_unchanged_max,
        "exterior_branch": exterior_rows,
        "exterior_branch_max_abs": exterior_max,
        "eta_derivative_fd_step": step,
        "eta_derivative_fd_max_abs": derivative_max,
        "eta_derivative_rows": derivative_rows,
        "pressure_owner": "JoinedOuterPressure.pressure unchanged",
        "five_moment_closure": False,
        "admissible_stress": False,
        "finite_energy_certified": False,
        "pde_validated": False,
        "scale_recursion_established": False,
        "scope": (
            "Finite two-axial-moment repair of the saved degree-four core. "
            "Independent Gauss checks cover off-grid eta values, exact core and "
            "exterior branches, and eta derivatives. No Part I O.1-O.8 outer "
            "construction, collar, all-five-moment closure, stress cone, PDE, "
            "energy, forcing, or recursion claim."
        ),
        "implementation_metadata": metadata,
    }
    destination = Path(__file__).with_suffix(".json")
    destination.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": report["status"],
        "offgrid_max_abs": report["independent_offgrid_max_abs_order128"],
        "eta_derivative_fd_max_abs": derivative_max,
        "matrix_conditioning": conditioning,
    }))
    return report


if __name__ == "__main__":
    run()
