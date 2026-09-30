"""Independent checks for the Section 3 profile-derived stress evaluator.

The checks use two supplied-profile receipts:

* the exact heat exterior with its finite-radius angular and quadratic tail
  targets; and
* the implemented inward heat collar, with the actual exterior-target
  moments transported to each test radius.

The collar comparison checks the new moment formulas against the existing
ODE-integrated ``HeatCollar`` stress.  It does not claim the missing core,
global five-moment closure, PDE, or the complete outer construction.
"""

from __future__ import annotations

import json
import math
from functools import lru_cache
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
LOCAL = Path(__file__).resolve().parent
for path in (SRC, LOCAL):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from lei_ren_part1_exterior_targets import exterior_targets  # noqa: E402
from lei_ren_part1_heat_collar import HeatCollar  # noqa: E402
from lei_ren_part1_heat_moments import heat_tail_moments  # noqa: E402
from lei_ren_part1_profile_stress import (  # noqa: E402
    MOMENT_KEYS,
    ProfileStress,
    SOURCE,
    SOURCE_SECTIONS,
    SOURCE_VERSION,
)


Z_POINTS = (-0.5, 0.0, 0.5)
HEAT_RADII = (100.0, 500.0, 2048.0)
COLLAR_Y = (0.0, 0.12, 0.25, 0.4, 0.5)


class HeatTailProfile:
    """Heat boundary data with the backward tail targets used in Part I."""

    def __init__(self, collar: HeatCollar) -> None:
        self.collar = collar
        self.heat = collar.heat
        self.delta = collar.delta

    def F(self, R: float, Z: float) -> float:
        return self.collar.F_heat(R, Z)

    def F_R(self, R: float, Z: float) -> float:
        return self.collar.F_heat_R(R, Z)

    def F_eta(self, R: float, Z: float) -> float:
        return self.collar.F_heat_Z(R, Z)

    def U(self, R: float, Z: float) -> float:
        del R, Z
        return 0.0

    def U_R(self, R: float, Z: float) -> float:
        del R, Z
        return 0.0

    def U_eta(self, R: float, Z: float) -> float:
        del R, Z
        return 0.0

    @lru_cache(maxsize=4096)
    def _tail(self, R: float, Z: float) -> dict[str, float]:
        return heat_tail_moments(R, Z, self.heat, series_order=12)

    def moments(self, R: float, Z: float) -> dict[str, float]:
        tail = self._tail(float(R), float(Z))
        # The tail routine's quadratic_tail is the negative future integral;
        # the cumulative inward target is therefore -quadratic_tail.  Mp is
        # represented relative to its terminal pressure normalization; P is
        # supplied independently below and is the quantity used in (3.17).
        return {
            "theta": float(tail["angular_target"]),
            "z": 0.0,
            "theta_z": 0.0,
            "z_theta": float(-tail["quadratic_tail"]),
            "p": float(self.P(R, Z)),
        }

    def moments_Z(self, R: float, Z: float) -> dict[str, float]:
        R = float(R)
        Z = float(Z)
        step = 2.0e-5
        lower = self._tail(R, Z - step)
        upper = self._tail(R, Z + step)
        return {
            "theta": float((upper["angular_target"] - lower["angular_target"]) / (2.0 * step)),
            "z": 0.0,
            "theta_z": 0.0,
            "z_theta": float(
                -(
                    upper["quadratic_tail"] - lower["quadratic_tail"]
                )
                / (2.0 * step)
            ),
            "p": float(self.P_Z(R, Z)),
        }

    def P(self, R: float, Z: float) -> float:
        return self.collar.P_heat(R, Z)

    def P_Z(self, R: float, Z: float) -> float:
        return self.collar.P_heat_Z(R, Z)


class CollarTargetProfile:
    """Implemented collar with its actual backward exterior moment targets."""

    def __init__(self, collar: HeatCollar) -> None:
        self.collar = collar
        self.heat = collar.heat
        self.delta = collar.delta
        self.R_b = collar.R_b

    def F(self, R: float, Z: float) -> float:
        return self.collar.F(R, Z)

    def F_R(self, R: float, Z: float) -> float:
        return self.collar.F_R(R, Z)

    def F_eta(self, R: float, Z: float) -> float:
        return self.collar.F_Z(R, Z)

    def U(self, R: float, Z: float) -> float:
        del R, Z
        return 0.0

    def U_R(self, R: float, Z: float) -> float:
        del R, Z
        return 0.0

    def U_eta(self, R: float, Z: float) -> float:
        del R, Z
        return 0.0

    @lru_cache(maxsize=4096)
    def _target(self, R: float, Z: float) -> dict[str, float]:
        result = exterior_targets(
            self.collar.F,
            self.heat,
            float(R),
            self.R_b,
            float(Z),
            n=64,
            P0=0.0,
        )
        return {key: float(value) for key, value in result.items() if isinstance(value, (int, float))}

    def moments(self, R: float, Z: float) -> dict[str, float]:
        target = self._target(float(R), float(Z))
        return {
            "theta": target["angular_target"],
            "z": target["axial_target"],
            "theta_z": target["mixed_target"],
            "z_theta": target["quadratic_target"],
            "p": target["pressure_moment_target"],
        }

    def moments_Z(self, R: float, Z: float) -> dict[str, float]:
        R = float(R)
        Z = float(Z)
        step = 2.0e-5
        lower = self.moments(R, Z - step)
        upper = self.moments(R, Z + step)
        return {
            key: float((upper[key] - lower[key]) / (2.0 * step))
            for key in MOMENT_KEYS
        }

    def P(self, R: float, Z: float) -> float:
        return self.collar.P(R, Z)

    def P_Z(self, R: float, Z: float) -> float:
        return self.collar.P_Z(R, Z)


def _max_abs(rows: list[dict[str, float]], key: str) -> float:
    return max(abs(float(row[key])) for row in rows)


def _heat_checks(stress: ProfileStress, profile: HeatTailProfile) -> dict[str, object]:
    rows: list[dict[str, float]] = []
    moment_errors: list[float] = []
    pressure_errors: list[float] = []
    stress_rows: list[dict[str, float]] = []
    for R in HEAT_RADII:
        for Z in Z_POINTS:
            data = stress.evaluate(R, Z)
            tail = profile._tail(float(R), float(Z))
            moment_errors.extend(
                [
                    abs(data["moments"]["z"]),
                    abs(data["moments"]["theta_z"]),
                    abs(data["moments"]["z_theta"] + tail["quadratic_tail"]),
                    abs(data["moments"]["theta"] - tail["angular_target"]),
                ]
            )
            pressure_errors.append(abs(float(data["P"]) - profile.P(R, Z)))
            rows.append(
                {
                    "R": R,
                    "Z": Z,
                    "T_theta": float(data["T_theta"]),
                    "T_z": float(data["T_z"]),
                    "I_theta": float(data["I_theta"]),
                    "I_z": float(data["I_z"]),
                    "S_theta": float(data["S_theta"]),
                    "S_z": float(data["S_z"]),
                }
            )
            # The supplied moment derivatives are central finite differences;
            # retain this radial-equation result as numerical evidence only.
            step = R * 2.0e-5
            left = stress.evaluate(R - step, Z)
            right = stress.evaluate(R + step, Z)
            stress_rows.append(
                {
                    "R": R,
                    "Z": Z,
                    "theta_ode_residual": R * (right["I_theta"] - left["I_theta"]) / (2.0 * step) + data["I_theta"] - data["N_theta"],
                    "z_ode_residual": R * (right["I_z"] - left["I_z"]) / (2.0 * step) + 0.5 * data["I_z"] - data["N_z"],
                }
            )
    return {
        "rows": rows,
        "moment_identity_max_absolute_error": max(moment_errors),
        "pressure_heat_max_absolute_error": max(pressure_errors),
        "max_abs_T_theta": _max_abs(rows, "T_theta"),
        "max_abs_T_z": _max_abs(rows, "T_z"),
        "radial_ode_fd_rows": stress_rows,
        "radial_ode_fd_max_theta_absolute": _max_abs(stress_rows, "theta_ode_residual"),
        "radial_ode_fd_max_z_absolute": _max_abs(stress_rows, "z_ode_residual"),
        "moment_origin_note": "heat angular/quadratic values are the supplied backward tail targets; no axis core is inferred",
    }


def _collar_checks(stress: ProfileStress, profile: CollarTargetProfile) -> dict[str, object]:
    rows: list[dict[str, float]] = []
    moment_errors: list[float] = []
    stress_errors: dict[str, list[float]] = {
        "I_theta": [],
        "I_z": [],
        "S_theta": [],
        "S_z": [],
        "T_theta": [],
        "T_z": [],
    }
    for y in COLLAR_Y:
        R = profile.R_b * math.exp(-y)
        for Z in Z_POINTS:
            data = stress.evaluate(R, Z)
            target = profile._target(R, Z)
            expected = {
                "I_theta": profile.collar.I_theta(R, Z),
                "I_z": profile.collar.I_z(R, Z),
                "S_theta": profile.collar.S_theta(R, Z),
                "S_z": profile.collar.S_z(R, Z),
                "T_theta": profile.collar.T_theta(R, Z),
                "T_z": profile.collar.T_z(R, Z),
            }
            for key in stress_errors:
                stress_errors[key].append(abs(float(data[key]) - expected[key]))
            moment_errors.extend(
                [
                    abs(float(data["moments"]["theta"]) - target["angular_target"]),
                    abs(float(data["moments"]["z_theta"]) - target["quadratic_target"]),
                    abs(float(data["moments"]["p"]) - target["pressure_moment_target"]),
                ]
            )
            rows.append(
                {
                    "y": y,
                    "R": R,
                    "Z": Z,
                    "profile_T_theta": float(data["T_theta"]),
                    "collar_T_theta": expected["T_theta"],
                    "profile_T_z": float(data["T_z"]),
                    "collar_T_z": expected["T_z"],
                }
            )
    return {
        "rows": rows,
        "moment_target_max_absolute_error": max(moment_errors),
        "stress_match_max_absolute_error": {
            key: max(values) for key, values in stress_errors.items()
        },
        "stress_match_max_relative_error": {
            key: max(
                abs(values[index])
                / max(abs(float(rows[index]["profile_T_theta"])), 1.0e-30)
                for index in range(len(values))
            )
            if key == "T_theta"
            else None
            for key, values in stress_errors.items()
        },
        "moment_source": "exterior_targets(collar.F, heat, R, R_b, Z, P0=0)",
    }


def _derivative_checks(
    stress_heat: ProfileStress,
    heat_profile: HeatTailProfile,
    stress_collar: ProfileStress,
    collar_profile: CollarTargetProfile,
) -> dict[str, object]:
    rows: list[dict[str, float | str]] = []
    p_errors: list[float] = []
    m_errors: list[float] = []
    for label, stress, profile, radii in (
        ("heat", stress_heat, heat_profile, HEAT_RADII),
        ("collar", stress_collar, collar_profile, tuple(collar_profile.R_b * math.exp(-y) for y in COLLAR_Y[1:])),
    ):
        for R in radii:
            for Z in Z_POINTS:
                step = 1.0e-4
                p_fd = (profile.P(R, Z + step) - profile.P(R, Z - step)) / (2.0 * step)
                p_error = abs(p_fd - profile.P_Z(R, Z))
                p_errors.append(p_error)
                m_fd = {}
                lower = profile.moments(R, Z - step)
                upper = profile.moments(R, Z + step)
                supplied = profile.moments_Z(R, Z)
                for key in MOMENT_KEYS:
                    m_fd[key] = (upper[key] - lower[key]) / (2.0 * step)
                    m_errors.append(abs(m_fd[key] - supplied[key]))
                rows.append(
                    {
                        "profile": label,
                        "R": R,
                        "Z": Z,
                        "P_Z_fd_absolute_error": p_error,
                        "moment_Z_fd_max_absolute_error": max(
                            abs(m_fd[key] - supplied[key]) for key in MOMENT_KEYS
                        ),
                    }
                )
    return {
        "rows": rows,
        "pressure_Z_fd_max_absolute_error": max(p_errors),
        "moment_Z_fd_max_absolute_error": max(m_errors),
        "difference_method": "independent central finite differences with step 1e-4",
    }


def run() -> dict[str, object]:
    collar = HeatCollar()
    heat_profile = HeatTailProfile(collar)
    collar_profile = CollarTargetProfile(collar)
    heat_stress = ProfileStress(heat_profile)
    collar_stress = ProfileStress(collar_profile)

    heat = _heat_checks(heat_stress, heat_profile)
    collar_result = _collar_checks(collar_stress, collar_profile)
    derivatives = _derivative_checks(
        heat_stress,
        heat_profile,
        collar_stress,
        collar_profile,
    )

    thresholds = {
        "heat_moment_identity_max_absolute": 2.0e-12,
        "heat_pressure_max_absolute": 2.0e-12,
        "heat_T_theta_max_absolute": 2.0e-9,
        "heat_T_z_max_absolute": 2.0e-9,
        "collar_moment_target_max_absolute": 2.0e-11,
        "collar_I_theta_match_max_absolute": 2.0e-8,
        "collar_I_z_match_max_absolute": 2.0e-8,
        "collar_S_theta_match_max_absolute": 2.0e-12,
        "collar_S_z_match_max_absolute": 2.0e-12,
        "collar_T_theta_match_max_absolute": 2.0e-8,
        "collar_T_z_match_max_absolute": 2.0e-8,
        "pressure_Z_fd_max_absolute": 2.0e-8,
        "moment_Z_fd_max_absolute": 2.0e-6,
    }
    collar_errors = collar_result["stress_match_max_absolute_error"]
    checks_passed = bool(
        heat["moment_identity_max_absolute_error"] < thresholds["heat_moment_identity_max_absolute"]
        and heat["pressure_heat_max_absolute_error"] < thresholds["heat_pressure_max_absolute"]
        and heat["max_abs_T_theta"] < thresholds["heat_T_theta_max_absolute"]
        and heat["max_abs_T_z"] < thresholds["heat_T_z_max_absolute"]
        and collar_result["moment_target_max_absolute_error"] < thresholds["collar_moment_target_max_absolute"]
        and collar_errors["I_theta"] < thresholds["collar_I_theta_match_max_absolute"]
        and collar_errors["I_z"] < thresholds["collar_I_z_match_max_absolute"]
        and collar_errors["S_theta"] < thresholds["collar_S_theta_match_max_absolute"]
        and collar_errors["S_z"] < thresholds["collar_S_z_match_max_absolute"]
        and collar_errors["T_theta"] < thresholds["collar_T_theta_match_max_absolute"]
        and collar_errors["T_z"] < thresholds["collar_T_z_match_max_absolute"]
        and derivatives["pressure_Z_fd_max_absolute_error"] < thresholds["pressure_Z_fd_max_absolute"]
        and derivatives["moment_Z_fd_max_absolute_error"] < thresholds["moment_Z_fd_max_absolute"]
    )
    report: dict[str, object] = {
        "source": SOURCE,
        "source_version": SOURCE_VERSION,
        "source_sections": SOURCE_SECTIONS,
        "status": "completed finite profile-derived stress checks",
        "parameters": collar.metadata(),
        "api_contract": {
            "profile_methods": [
                "F(R,Z)", "U(R,Z)", "F_R(R,Z)", "F_eta(R,Z)",
                "U_R(R,Z)", "U_eta(R,Z)", "moments(R,Z)",
                "moments_Z(R,Z)", "P(R,Z)", "P_Z(R,Z)",
            ],
            "moment_keys": MOMENT_KEYS,
            "stress_formulas": "3.9, 3.12, 3.16, 3.17, 3.18",
        },
        "exact_heat": heat,
        "collar_target_match": collar_result,
        "independent_Z_derivatives": derivatives,
        "thresholds": thresholds,
        "scope": {
            "full_five_moment_closure": False,
            "axis_core_match": False,
            "full_outer_construction": False,
            "PDE_validation": False,
            "flatness_certification": False,
        },
        "checks_passed": checks_passed,
    }
    output = Path(__file__).with_suffix(".json")
    output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    summary = {
        "checks_passed": checks_passed,
        "heat_max_abs_T_theta": heat["max_abs_T_theta"],
        "heat_max_abs_T_z": heat["max_abs_T_z"],
        "collar_stress_match_max": collar_result["stress_match_max_absolute_error"],
        "pressure_Z_fd_max_absolute_error": derivatives["pressure_Z_fd_max_absolute_error"],
        "moment_Z_fd_max_absolute_error": derivatives["moment_Z_fd_max_absolute_error"],
    }
    print(json.dumps(summary, indent=2))
    if not checks_passed:
        raise AssertionError("profile-derived stress checks failed")
    return report


if __name__ == "__main__":
    run()
