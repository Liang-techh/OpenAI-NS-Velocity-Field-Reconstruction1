"""Analytic and refinement checks for the finite Lei--Ren moment helper.

The checks use low-degree polynomial profiles whose five integrals are known
from coefficient-wise antiderivatives.  They exercise vector broadcasting,
pressure consistency, finite-input validation, and quadrature refinement.
They do not certify a Navier--Stokes field or any infinite exterior.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from lei_ren_part1_moments import (
    MOMENT_NAMES,
    SOURCE,
    SOURCE_SECTION,
    SOURCE_VERSION,
    cumulative_five_moments,
)


def _integrate_polynomial(coefficients: np.ndarray, rmax: np.ndarray) -> np.ndarray:
    """Integrate R-polynomial coefficients from 0 to each terminal radius."""

    result = np.zeros_like(rmax, dtype=float)
    for power, coefficient in enumerate(coefficients):
        result = result + coefficient * rmax ** (power + 1) / (power + 1)
    return result


def _poly_mul(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    """Multiply R-polynomials whose coefficient axis is the first axis."""

    product = np.zeros((left.shape[0] + right.shape[0] - 1,) + left.shape[1:])
    for i, left_coefficient in enumerate(left):
        for j, right_coefficient in enumerate(right):
            product[i + j] = product[i + j] + left_coefficient * right_coefficient
    return product


def _exact_moments(rmax: np.ndarray, z: np.ndarray) -> dict[str, np.ndarray]:
    """Return exact moments for the manufactured polynomial profiles."""

    # F = 1 + (Z/2) R + R^2/4; Uz = (-3/10 + Z/5) + 2 R^2/5.
    f = np.asarray([np.ones_like(z), z / 2.0, np.full_like(z, 1.0 / 4.0)])
    uz = np.asarray([
        -3.0 / 10.0 + z / 5.0,
        np.zeros_like(z),
        np.full_like(z, 2.0 / 5.0),
    ])
    angular = np.pad(2.0 * f, ((1, 0), (0, 0)))
    mixed = np.pad(2.0 * _poly_mul(f, uz), ((1, 0), (0, 0)))
    quadratic = np.pad(_poly_mul(uz, uz), ((0, 1), (0, 0))) - np.pad(
        _poly_mul(f, f), ((1, 0), (0, 0))
    )
    pressure = _poly_mul(f, f)
    return {
        "angular": _integrate_polynomial(angular, rmax),
        "axial": _integrate_polynomial(uz, rmax),
        "mixed": _integrate_polynomial(mixed, rmax),
        "quadratic": _integrate_polynomial(quadratic, rmax),
        "pressure": _integrate_polynomial(pressure, rmax),
    }


def _max_abs_difference(left: dict[str, np.ndarray],
                        right: dict[str, np.ndarray]) -> dict[str, float]:
    return {
        name: float(np.max(np.abs(np.asarray(left[name]) - np.asarray(right[name]))))
        for name in MOMENT_NAMES
    }


def run() -> dict[str, object]:
    rmax = np.asarray([0.5, 1.1, 1.7], dtype=float)
    z = np.asarray([-0.4, 0.15, 0.8], dtype=float)

    def F(rho, eta):
        return 1.0 + eta * rho / 2.0 + rho * rho / 4.0

    def Uz(rho, eta):
        return -3.0 / 10.0 + eta / 5.0 + 2.0 * rho * rho / 5.0

    expected = _exact_moments(rmax, z)
    coarse = cumulative_five_moments(F, Uz, rmax, z, n=2)
    resolved = cumulative_five_moments(F, Uz, rmax, z, n=8)
    exact_error = _max_abs_difference(resolved, expected)
    coarse_error = _max_abs_difference(coarse, expected)
    assert max(exact_error.values()) < 2.0e-13
    assert max(coarse_error.values()) > 1.0e-8

    # A second shape pair exercises NumPy broadcasting (Rmax shape (2, 1),
    # Z shape (3,)) and the optional pressure identity.
    broadcast_rmax = np.asarray([[0.3], [0.9]])
    broadcast_z = np.asarray([-0.6, 0.0, 0.6])
    broadcast = cumulative_five_moments(
        F,
        Uz,
        broadcast_rmax,
        broadcast_z,
        n=8,
        P0=lambda eta: 0.25 - eta / 7.0,
        P=lambda radius, eta: (0.25 - eta / 7.0)
        + cumulative_five_moments(F, Uz, radius, eta, n=8)["pressure"],
    )
    assert broadcast["angular"].shape == (2, 3)
    pressure_defect_max = float(np.max(np.abs(broadcast["pressure_defect"])))
    assert pressure_defect_max < 2.0e-13

    # Explicitly verify fail-closed validation for finite input and interval
    # requirements.  The exact error text is intentionally not part of the
    # interface, so these checks only require the documented exception type.
    validation_checks: dict[str, str] = {}
    for label, thunk in {
        "negative_Rmax": lambda: cumulative_five_moments(F, Uz, -1.0, 0.0),
        "nonfinite_Z": lambda: cumulative_five_moments(F, Uz, 1.0, np.nan),
        "nonpositive_n": lambda: cumulative_five_moments(F, Uz, 1.0, 0.0, n=0),
        "incomplete_pressure": lambda: cumulative_five_moments(
            F, Uz, 1.0, 0.0, P0=0.0
        ),
    }.items():
        try:
            thunk()
        except (TypeError, ValueError):
            validation_checks[label] = "passed"
        else:
            raise AssertionError(f"validation did not reject {label}")

    report = {
        "status": "completed",
        "source": SOURCE,
        "source_version": SOURCE_VERSION,
        "source_section": SOURCE_SECTION,
        "quadrature_orders": [2, 8],
        "manufactured_profile": {
            "F": "1 + (Z/2) R + R^2/4",
            "Uz": "-3/10 + Z/5 + 2 R^2/5",
            "Rmax": rmax.tolist(),
            "Z": z.tolist(),
        },
        "exact_error_at_n8": exact_error,
        "coarse_error_at_n2": coarse_error,
        "broadcast_shape": list(broadcast["angular"].shape),
        "pressure_defect_max_abs": pressure_defect_max,
        "validation_checks": validation_checks,
        "pde_validated": False,
        "scale_recursion_established": False,
        "scope": (
            "Finite-interval Gauss--Legendre checks for the five cumulative "
            "moments only; no infinite tail, heat-exterior certification, "
            "PDE validation, or scale-recursion claim."
        ),
    }
    output = Path(__file__).with_suffix(".json")
    output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report))
    return report


if __name__ == "__main__":
    run()
