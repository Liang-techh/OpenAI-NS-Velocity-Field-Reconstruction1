"""Independent finite checks for the Lei--Ren Part I heat collar.

The receipt checks the supplied heat boundary, the radial pressure identity,
the Section 3 stress ODE, the normalized cone inequalities, and the finite
Laguerre and Taylor numerical choices.  It is deliberately a local collar
check: it does not certify moment closure, the PDE, the full outer intervals,
or a global flatness construction.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
LOCAL = Path(__file__).resolve().parent
for path in (SRC, LOCAL):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from lei_ren_part1_heat_collar import (  # noqa: E402
    HeatCollar,
    SOURCE,
    SOURCE_SECTIONS,
    SOURCE_VERSION,
)
from openai_ns_reconstruction.heat_exterior import heat_factor  # noqa: E402


Z_CHECKS = (-0.5, 0.0, 0.5)
Y_ODE = (0.2, 0.3, 0.45)


def _relative_error(observed: float, expected: float) -> float:
    return abs(observed - expected) / max(abs(expected), 1.0e-30)


def _pressure_identity(collar: HeatCollar) -> dict[str, object]:
    rows: list[dict[str, float]] = []
    for y in Y_ODE:
        radius = collar.R_at_y(y)
        dr = radius * 1.0e-5
        for z in Z_CHECKS:
            observed = (collar.P(radius + dr, z) - collar.P(radius - dr, z)) / (2.0 * dr)
            expected = collar.F(radius, z) ** 2
            rows.append(
                {
                    "y": y,
                    "Z": z,
                    "finite_difference": observed,
                    "F_squared": expected,
                    "absolute_error": abs(observed - expected),
                    "relative_error": _relative_error(observed, expected),
                }
            )
    return {
        "rows": rows,
        "max_absolute_error": max(row["absolute_error"] for row in rows),
        "max_relative_error": max(row["relative_error"] for row in rows),
    }


def _stress_ode(collar: HeatCollar) -> dict[str, object]:
    rows: list[dict[str, float]] = []
    for y in Y_ODE:
        step = 2.0e-5
        for z in Z_CHECKS:
            r_lo = collar.R_at_y(y - step)
            r_hi = collar.R_at_y(y + step)
            dtheta_dy = (collar.I_theta(r_hi, z) - collar.I_theta(r_lo, z)) / (2.0 * step)
            dz_dy = (collar.I_z(r_hi, z) - collar.I_z(r_lo, z)) / (2.0 * step)
            theta_residual = -dtheta_dy + collar.I_theta(collar.R_at_y(y), z) - collar.N_theta(collar.R_at_y(y), z)
            z_residual = -dz_dy + 0.5 * collar.I_z(collar.R_at_y(y), z) - collar.N_z(collar.R_at_y(y), z)
            rows.append(
                {
                    "y": y,
                    "Z": z,
                    "theta_residual": theta_residual,
                    "z_residual": z_residual,
                    "theta_scale": max(abs(collar.N_theta(collar.R_at_y(y), z)), 1.0e-30),
                    "z_scale": max(abs(collar.N_z(collar.R_at_y(y), z)), 1.0e-30),
                }
            )
    return {
        "rows": rows,
        "max_theta_absolute_residual": max(abs(row["theta_residual"]) for row in rows),
        "max_z_absolute_residual": max(abs(row["z_residual"]) for row in rows),
    }


def _heat_taylor(collar: HeatCollar) -> dict[str, object]:
    rows: list[dict[str, float | str | None]] = []
    errors: list[float] = []
    derivative_errors: list[float] = []
    for xi in (0.0, 1.0e-4, 5.0e-4, 1.0e-3, 2.0e-3, 3.0e-3):
        jet = collar.heat_jet(xi)
        exact = heat_factor(xi, collar.h)
        exact_derivative = heat_factor(xi, collar.h, derivative=1)
        value_error = abs(float(jet["value"]) - exact)
        derivative_error = abs(float(jet["derivative"]) - exact_derivative)
        errors.append(value_error)
        derivative_errors.append(derivative_error)
        rows.append(
            {
                "xi": xi,
                "method": jet["method"],
                "value_error": value_error,
                "derivative_error": derivative_error,
                "value_error_bound": jet["value_error_bound"],
                "derivative_error_bound": jet["derivative_error_bound"],
            }
        )
    return {
        "rows": rows,
        "max_value_error": max(errors),
        "max_derivative_error": max(derivative_errors),
    }


def _laguerre_comparison(collar: HeatCollar) -> dict[str, object]:
    coarse = HeatCollar(
        h=collar.h,
        c_inf=collar.c_inf,
        R_b=collar.R_b,
        ell=collar.ell,
        epsilon=collar.epsilon,
        laguerre_order=32,
    )
    fine = HeatCollar(
        h=collar.h,
        c_inf=collar.c_inf,
        R_b=collar.R_b,
        ell=collar.ell,
        epsilon=collar.epsilon,
        laguerre_order=64,
    )
    rows: list[dict[str, float]] = []
    differences: dict[str, list[float]] = {
        "btheta": [],
        "bz": [],
        "Tz_over_Ttheta": [],
        "delta_P_over_flat": [],
        "delta_P_Z_over_flat": [],
    }
    for y in (1.0e-4, 1.0e-2, 0.2, 0.3, 0.45):
        for z in Z_CHECKS:
            radius = collar.R_at_y(y)
            row: dict[str, float] = {"y": y, "Z": z}
            coarse_data = coarse.normalized(radius, z)
            fine_data = fine.normalized(radius, z)
            for key in ("btheta", "bz", "Tz_over_Ttheta"):
                difference = abs(float(coarse_data[key]) - float(fine_data[key]))
                differences[key].append(difference)
                row[f"{key}_difference"] = difference
            for key, method in (
                ("delta_P_over_flat", "_delta_P_over_flat"),
                ("delta_P_Z_over_flat", "_delta_P_Z_over_flat"),
            ):
                coarse_value = float(getattr(coarse, method)(y, z))
                fine_value = float(getattr(fine, method)(y, z))
                difference = abs(coarse_value - fine_value)
                differences[key].append(difference)
                row[f"{key}_difference"] = difference
            rows.append(row)
    return {
        "rows": rows,
        "max_absolute_difference": {
            key: max(values) for key, values in differences.items()
        },
    }


def _cone_grid(collar: HeatCollar) -> dict[str, object]:
    rows: list[dict[str, float]] = []
    y_values = np.concatenate(([0.0], np.geomspace(1.0e-5, collar.ell, 9)))
    z_values = np.linspace(-1.0, 1.0, 17)
    for y in y_values:
        radius = collar.R_at_y(float(y))
        for z in z_values:
            data = collar.normalized(radius, float(z))
            rows.append(
                {
                    "y": float(y),
                    "Z": float(z),
                    "btheta": float(data["btheta"]),
                    "bz": float(data["bz"]),
                    "Tz_over_Ttheta": float(data["Tz_over_Ttheta"]),
                    "directional_margin": float(data["directional_margin"]),
                    "kappa_margin": float(data["kappa_margin"]),
                }
            )
    return {
        "rows": rows,
        "min_btheta": min(row["btheta"] for row in rows),
        "max_abs_bz": max(abs(row["bz"]) for row in rows),
        "max_abs_ratio": max(abs(row["Tz_over_Ttheta"]) for row in rows),
        "min_directional_margin": min(row["directional_margin"] for row in rows),
        "min_kappa_margin": min(row["kappa_margin"] for row in rows),
        "all_btheta_positive": all(row["btheta"] > 0.0 for row in rows),
        "all_directional_margin_positive": all(row["directional_margin"] > 0.0 for row in rows),
        "all_kappa_margin_positive": all(row["kappa_margin"] > 0.0 for row in rows),
    }


def _edge_limits(collar: HeatCollar) -> dict[str, object]:
    rows: list[dict[str, float]] = []
    errors_theta: list[float] = []
    errors_z: list[float] = []
    for y in (1.0e-5, 3.0e-5, 1.0e-4, 3.0e-4, 1.0e-3, 3.0e-3, 1.0e-2):
        for z in (-1.0, -0.5, 0.0, 0.5, 1.0):
            radius = collar.R_at_y(y)
            endpoint = collar.normalized(collar.R_b, z)
            data = collar.normalized(radius, z)
            theta_error = abs(float(data["btheta"]) - float(endpoint["btheta"]))
            z_error = abs(float(data["bz"]) - float(endpoint["bz"]))
            errors_theta.append(theta_error)
            errors_z.append(z_error)
            rows.append(
                {
                    "y": y,
                    "Z": z,
                    "btheta_error_to_y0": theta_error,
                    "bz_error_to_y0": z_error,
                }
            )
    # The endpoint formulas freeze the heat factor at R_b.  At y>0 the heat
    # factor itself changes by O(y), so the finite-distance discrepancy is
    # tested in the asymptotic band y <= 5e-5 and reported as an O(y) slope
    # over the full geometric grid below.
    near_rows = [row for row in rows if row["y"] <= 5.0e-5]
    return {
        "rows": rows,
        "edge_y_threshold": 5.0e-5,
        "max_btheta_error": max(row["btheta_error_to_y0"] for row in near_rows),
        "max_bz_error": max(row["bz_error_to_y0"] for row in near_rows),
        "max_btheta_error_over_y": max(
            row["btheta_error_to_y0"] / row["y"] for row in rows
        ),
        "max_bz_error_over_y": max(
            row["bz_error_to_y0"] / row["y"] for row in rows
        ),
        "all_errors_linear_in_y": bool(
            max(row["btheta_error_to_y0"] / row["y"] for row in rows) < 2.0e-8
            and max(row["bz_error_to_y0"] / row["y"] for row in rows) < 5.0e-10
        ),
        "stress_at_y0_exactly_zero": bool(
            all(collar.T_theta(collar.R_b, z) == 0.0 for z in (-1.0, 0.0, 1.0))
            and all(collar.T_z(collar.R_b, z) == 0.0 for z in (-1.0, 0.0, 1.0))
        ),
    }


def run() -> dict[str, object]:
    collar = HeatCollar()
    pressure = _pressure_identity(collar)
    stress = _stress_ode(collar)
    heat = _heat_taylor(collar)
    laguerre = _laguerre_comparison(collar)
    cone = _cone_grid(collar)
    edge = _edge_limits(collar)

    thresholds = {
        "pressure_identity_max_relative": 2.0e-7,
        "stress_theta_max_absolute": 1.0e-10,
        "stress_z_max_absolute": 1.0e-10,
        "heat_value_max_absolute": 1.0e-14,
        "heat_derivative_max_absolute": 1.0e-13,
        "laguerre_btheta_max_absolute": 1.0e-14,
        "laguerre_bz_max_absolute": 1.0e-16,
        "laguerre_ratio_max_absolute": 1.0e-13,
        "laguerre_pressure_max_absolute": 1.0e-14,
        "laguerre_pressure_z_max_absolute": 1.0e-14,
        "edge_btheta_max_absolute": 2.0e-12,
        "edge_bz_max_absolute": 2.0e-14,
    }
    checks_passed = bool(
        pressure["max_relative_error"] < thresholds["pressure_identity_max_relative"]
        and stress["max_theta_absolute_residual"] < thresholds["stress_theta_max_absolute"]
        and stress["max_z_absolute_residual"] < thresholds["stress_z_max_absolute"]
        and heat["max_value_error"] < thresholds["heat_value_max_absolute"]
        and heat["max_derivative_error"] < thresholds["heat_derivative_max_absolute"]
        and laguerre["max_absolute_difference"]["btheta"] < thresholds["laguerre_btheta_max_absolute"]
        and laguerre["max_absolute_difference"]["bz"] < thresholds["laguerre_bz_max_absolute"]
        and laguerre["max_absolute_difference"]["Tz_over_Ttheta"] < thresholds["laguerre_ratio_max_absolute"]
        and laguerre["max_absolute_difference"]["delta_P_over_flat"] < thresholds["laguerre_pressure_max_absolute"]
        and laguerre["max_absolute_difference"]["delta_P_Z_over_flat"] < thresholds["laguerre_pressure_z_max_absolute"]
        and edge["max_btheta_error"] < thresholds["edge_btheta_max_absolute"]
        and edge["max_bz_error"] < thresholds["edge_bz_max_absolute"]
        and edge["all_errors_linear_in_y"]
        and cone["all_btheta_positive"]
        and cone["all_directional_margin_positive"]
        and cone["all_kappa_margin_positive"]
        and edge["stress_at_y0_exactly_zero"]
    )
    report: dict[str, object] = {
        "source": SOURCE,
        "source_version": SOURCE_VERSION,
        "source_sections": SOURCE_SECTIONS,
        "status": "completed finite numerical collar checks",
        "parameters": collar.metadata(),
        "pressure_radial_identity": pressure,
        "stress_ode": stress,
        "heat_taylor_vs_source": heat,
        "laguerre_32_vs_64": laguerre,
        "normalized_cone_grid": cone,
        "normalized_edge_limits": edge,
        "thresholds": thresholds,
        "scope": {
            "global_moment_matching": False,
            "full_outer_O1_O8": False,
            "axis_core_compatibility": False,
            "PDE_validation": False,
            "flatness_certification": False,
        },
        "checks_passed": checks_passed,
    }
    output = Path(__file__).with_suffix(".json")
    output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    summary = {
        key: value
        for key, value in report.items()
        if key not in {"parameters", "pressure_radial_identity", "stress_ode", "heat_taylor_vs_source", "laguerre_32_vs_64", "normalized_cone_grid", "normalized_edge_limits"}
    }
    summary["pressure_identity_max_relative"] = pressure["max_relative_error"]
    summary["stress_theta_max_absolute"] = stress["max_theta_absolute_residual"]
    summary["stress_z_max_absolute"] = stress["max_z_absolute_residual"]
    summary["cone_min_directional_margin"] = cone["min_directional_margin"]
    summary["cone_min_kappa_margin"] = cone["min_kappa_margin"]
    print(json.dumps(summary, indent=2))
    if not checks_passed:
        raise AssertionError("heat collar checks failed")
    return report


if __name__ == "__main__":
    run()
