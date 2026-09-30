"""Checks for actual SmoothExtendedAxial profile-derived stress.

The stress and five moments are evaluated from one actual ``F/U`` profile.
The heat/collar target moments used by earlier receipts are intentionally not
inserted here.  Cone booleans and tail residual comparisons are reported as
diagnostics; this file does not turn a partial candidate into a whole-cone or
five-moment acceptance claim.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
LOCAL = ROOT / "experiments" / "root_st073"
for path in (SRC, LOCAL):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from lei_ren_part1_extended_profile_stress import (  # noqa: E402
    ExtendedProfileStressAdapter,
    SOURCE,
    SOURCE_SECTIONS,
    SOURCE_VERSION,
)
from lei_ren_part1_smooth_extended_axial import SmoothExtendedAxial  # noqa: E402


RADII = (0.002, 0.03, 1.0, 100.0, 300.0, 600.0, 1000.0, 1500.0, 2300.0)
Z_VALUES = (0.0, 0.15, -0.3)
INTERIOR_Z = (-0.5, -0.15, 0.15, 0.5)
MOMENT_KEYS = ("theta", "z", "theta_z", "z_theta", "p")
STRESS_KEYS = ("I_theta", "I_z", "S_theta", "S_z", "T_theta", "T_z")


def _cone_data(data: dict[str, object]) -> dict[str, float | bool | None]:
    F = float(data["F"])
    S_theta = float(data["S_theta"])
    S_z = float(data["S_z"])
    T_theta = float(data["T_theta"])
    T_z = float(data["T_z"])
    shear_sq = S_theta * S_theta + S_z * S_z
    if abs(F * S_theta) <= 1.0e-30 or shear_sq <= 1.0e-30:
        return {
            "kappa": None,
            "T_dot_S": None,
            "T_dot_S_perp": None,
            "cone_margin": None,
            "kappa_gt_2": False,
            "T_dot_S_negative": False,
            "directional_cone": False,
        }
    kappa = -shear_sq / (F * S_theta)
    dot_s = T_theta * S_theta + T_z * S_z
    dot_s_perp = T_theta * (-S_z) + T_z * S_theta
    cone_margin = 2.0 * dot_s * dot_s - (kappa - 2.0) * dot_s_perp * dot_s_perp
    return {
        "kappa": kappa,
        "T_dot_S": dot_s,
        "T_dot_S_perp": dot_s_perp,
        "cone_margin": cone_margin,
        "kappa_gt_2": bool(kappa > 2.0),
        "T_dot_S_negative": bool(dot_s < 0.0),
        "directional_cone": bool(cone_margin > 0.0),
    }


def _point_rows(stress: object, adapter: ExtendedProfileStressAdapter) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for radius in RADII:
        for z in Z_VALUES:
            data = stress.evaluate(radius, z)
            cone = _cone_data(data)
            row: dict[str, object] = {
                "R": radius,
                "Z": z,
                "F": float(data["F"]),
                "U": float(data["U"]),
                "P": float(data["P"]),
                "U_r": float(data["U_r"]),
                "kappa": cone["kappa"],
                "T_dot_S": cone["T_dot_S"],
                "T_dot_S_perp": cone["T_dot_S_perp"],
                "cone_margin": cone["cone_margin"],
                "kappa_gt_2": cone["kappa_gt_2"],
                "T_dot_S_negative": cone["T_dot_S_negative"],
                "directional_cone": cone["directional_cone"],
            }
            for key in STRESS_KEYS:
                row[key] = float(data[key])
            row["moments"] = {
                key: float(data["moments"][key]) for key in MOMENT_KEYS
            }
            row["moments_Z"] = {
                key: float(data["moments_Z"][key]) for key in MOMENT_KEYS
            }
            rows.append(row)
    return rows


def _quadrature_refinement(
    stress32: object,
    stress64: object,
    rows64: list[dict[str, object]],
) -> dict[str, object]:
    rows: list[dict[str, object]] = []
    moment_errors: dict[str, list[float]] = {key: [] for key in MOMENT_KEYS}
    stress_errors: dict[str, list[float]] = {key: [] for key in STRESS_KEYS}
    for row64 in rows64:
        radius = float(row64["R"])
        z = float(row64["Z"])
        coarse = stress32.evaluate(radius, z)
        fine = stress64.evaluate(radius, z)
        row: dict[str, object] = {"R": radius, "Z": z}
        for key in MOMENT_KEYS:
            error = abs(float(coarse["moments"][key]) - float(fine["moments"][key]))
            moment_errors[key].append(error)
            row[f"moment_{key}_abs_difference"] = error
        for key in STRESS_KEYS:
            error = abs(float(coarse[key]) - float(fine[key]))
            stress_errors[key].append(error)
            row[f"stress_{key}_abs_difference"] = error
        rows.append(row)
    three_moment_errors = [
        value
        for key in ("theta", "theta_z", "z_theta")
        for value in moment_errors[key]
    ]
    return {
        "rows": rows,
        "max_abs_moment_difference": {
            key: max(values) for key, values in moment_errors.items()
        },
        "max_abs_stress_difference": {
            key: max(values) for key, values in stress_errors.items()
        },
        "three_moment_refinement_support_max_abs": max(three_moment_errors),
    }


def _independent_derivative_checks(adapter: ExtendedProfileStressAdapter) -> dict[str, object]:
    rows: list[dict[str, float | str]] = []
    pressure_errors: list[float] = []
    moment_errors: list[float] = []
    pressure_step_errors: list[float] = []
    moment_step_errors: dict[str, list[float]] = {key: [] for key in MOMENT_KEYS}
    for radius in RADII:
        for z in INTERIOR_Z:
            supplied_p = adapter.P_Z(radius, z)
            supplied = adapter.moments_Z(radius, z)
            fd_values: dict[float, dict[str, float]] = {}
            fd_pressure: dict[float, float] = {}
            for step in (2.0e-4, 1.0e-4):
                lower_p = adapter.P(radius, z - step)
                upper_p = adapter.P(radius, z + step)
                fd_pressure[step] = (upper_p - lower_p) / (2.0 * step)
                lower = adapter.moments(radius, z - step)
                upper = adapter.moments(radius, z + step)
                fd_values[step] = {
                    key: (upper[key] - lower[key]) / (2.0 * step)
                    for key in MOMENT_KEYS
                }
            p_error = abs(fd_pressure[1.0e-4] - supplied_p)
            pressure_errors.append(p_error)
            pressure_step_errors.append(
                abs(fd_pressure[1.0e-4] - fd_pressure[2.0e-4])
            )
            moment_fd_errors = {
                key: abs(fd_values[1.0e-4][key] - supplied[key])
                for key in MOMENT_KEYS
            }
            for key in MOMENT_KEYS:
                moment_errors.append(moment_fd_errors[key])
                moment_step_errors[key].append(
                    abs(fd_values[1.0e-4][key] - fd_values[2.0e-4][key])
                )
            rows.append(
                {
                    "R": radius,
                    "Z": z,
                    "P_Z_fd_error": p_error,
                    "moment_Z_fd_max_error": max(moment_fd_errors.values()),
                    "P_Z_fd_step_halving_difference": pressure_step_errors[-1],
                    "moment_Z_fd_step_halving_max_difference": max(
                        moment_step_errors[key][-1] for key in MOMENT_KEYS
                    ),
                }
            )
    return {
        "rows": rows,
        "difference_scheme": "independent centered second-order FD at steps 2e-4 and 1e-4; interior nonzero Z only",
        "pressure_Z_fd_max_abs_error": max(pressure_errors),
        "moment_Z_fd_max_abs_error": max(moment_errors),
        "pressure_Z_fd_step_halving_max_abs_difference": max(pressure_step_errors),
        "moment_Z_fd_step_halving_max_abs_difference": {
            key: max(values) for key, values in moment_step_errors.items()
        },
    }


def _tail_comparison(rows: list[dict[str, object]], refinement: dict[str, object]) -> dict[str, object]:
    tail = [row for row in rows if float(row["R"]) > 1242.1747910914733]
    tail_components = [
        abs(float(row["T_theta"]))
        for row in tail
    ] + [
        abs(float(row["T_z"]))
        for row in tail
    ]
    stress_support = [
        float(value)
        for key, value in refinement["max_abs_stress_difference"].items()
        if key in ("T_theta", "T_z")
    ]
    moment_support = float(refinement["three_moment_refinement_support_max_abs"])
    return {
        "tail_rows": tail,
        "tail_max_abs_T_component": max(tail_components),
        "n32_n64_stress_support_max_abs": max(stress_support),
        "n32_n64_three_moment_support_max_abs": moment_support,
        "tail_below_stress_refinement_support": bool(max(tail_components) <= max(stress_support)),
        "tail_below_three_moment_refinement_support": bool(max(tail_components) <= moment_support),
        "comparison_note": "Supports are numerical n32/n64 refinement ranges; no physical acceptance threshold is inferred from them.",
    }


def run() -> dict[str, object]:
    matched = SmoothExtendedAxial()
    adapter64 = ExtendedProfileStressAdapter(matched, radial_quadrature_order=64)
    adapter32 = ExtendedProfileStressAdapter(matched, radial_quadrature_order=32)
    stress64 = adapter64.profile
    stress32 = adapter32.profile
    rows64 = _point_rows(stress64, adapter64)
    refinement = _quadrature_refinement(stress32, stress64, rows64)
    derivatives = _independent_derivative_checks(adapter64)
    tail = _tail_comparison(rows64, refinement)

    finite = all(
        math.isfinite(float(row[key]))
        for row in rows64
        for key in ("F", "U", "P", "I_theta", "I_z", "S_theta", "S_z", "T_theta", "T_z")
    )
    adapter_checks_passed = bool(
        finite
        and refinement["max_abs_moment_difference"]["theta"] < 2.0e-6
        and refinement["max_abs_moment_difference"]["theta_z"] < 2.0e-6
        and refinement["max_abs_moment_difference"]["z_theta"] < 2.0e-6
        and derivatives["pressure_Z_fd_max_abs_error"] < 2.0e-6
    )
    cone_pass_count = sum(
        bool(row["kappa_gt_2"] and row["T_dot_S_negative"] and row["directional_cone"])
        for row in rows64
    )
    report: dict[str, object] = {
        "status": "completed actual-profile stress diagnostics",
        "source": SOURCE,
        "source_version": SOURCE_VERSION,
        "source_sections": SOURCE_SECTIONS,
        "parameters": adapter64.metadata(),
        "historical_run_metadata": {
            "original_candidate_metadata_unrecorded": False,
            "note": (
                "These numerical rows were regenerated with the current candidate; "
                "the complete matched.metadata() subobject is embedded under "
                "parameters.matched_metadata. The prior interpolated-coefficient "
                "rows are preserved separately in the historical receipt."
            ),
        },
        "points": {"R": list(RADII), "Z": list(Z_VALUES)},
        "actual_profile_rows_order64": rows64,
        "quadrature_refinement_n32_vs_n64": refinement,
        "independent_derivative_checks": derivatives,
        "tail_stress_comparison": tail,
        "cone_diagnostic": {
            "tested_point_count": len(rows64),
            "passing_point_count": cone_pass_count,
            "all_tested_points_pass": bool(cone_pass_count == len(rows64)),
            "whole_cone_claim": False,
        },
        "adapter_checks_passed": adapter_checks_passed,
        "scope": {
            "actual_F_U_moments_used": True,
            "heat_or_collar_moment_targets_substituted": False,
            "full_five_moment_closure": False,
            "whole_cone_validation": False,
            "pde_validation": False,
            "flatness_certification": False,
        },
    }
    output = Path(__file__).with_suffix(".json")
    output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "adapter_checks_passed": adapter_checks_passed,
        "cone_pass_count": cone_pass_count,
        "point_count": len(rows64),
        "tail_max_abs_T": tail["tail_max_abs_T_component"],
        "tail_below_stress_refinement_support": tail["tail_below_stress_refinement_support"],
        "tail_below_three_moment_refinement_support": tail["tail_below_three_moment_refinement_support"],
    }, indent=2))
    return report


if __name__ == "__main__":
    run()
