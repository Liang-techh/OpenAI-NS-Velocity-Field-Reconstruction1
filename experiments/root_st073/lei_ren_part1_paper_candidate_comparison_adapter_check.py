"""Smoke receipt for the nominal CandidateMPComparison adapter."""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import mpmath as mp


HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from lei_ren_part1_paper_candidate_comparison_adapter import CandidateMPComparison
from lei_ren_part1_paper_candidate_comparison_jets import CandidateComparisonJets
from lei_ren_part1_paper_candidate_exit_analytic import candidate_core


def relative(left: mp.mpf, right: mp.mpf) -> mp.mpf:
    scale = max(abs(left), abs(right), mp.mpf("1e-300000"))
    return abs(left - right) / scale


def run() -> dict[str, object]:
    calc = CandidateComparisonJets(32)
    adapter = CandidateMPComparison(calc, candidate_core(calc))
    with mp.workdps(adapter.precision + 40):
        z = ".3"
        endpoint = adapter.evaluate(2 * adapter.h_b, z)
        frozen = adapter.frozen_driver_coefficients(z)
        rows = []
        max_driver_error = mp.mpf(0)
        max_moment_error = mp.mpf(0)
        max_pressure_error = mp.mpf(0)
        for radius in (mp.mpf("1"), mp.mpf("100")):
            y = mp.log(radius / adapter.R_a)
            value = adapter.evaluate(y, z)
            D_expected = (
                frozen["D"][1] * radius
                + frozen["D"][0]
                + frozen["D"][-1] / radius
            )
            I_expected = (
                frozen["I_z"][3] * radius ** mp.mpf("1.5")
                + frozen["I_z"][1] * radius ** mp.mpf(".5")
                + frozen["I_z"][-1] * radius ** mp.mpf("-.5")
            )
            driver_error = max(
                relative(value["D"], D_expected),
                relative(value["I_z"], I_expected),
            )
            delta_R = radius - endpoint["R"]
            delta_R2 = radius * radius - endpoint["R"] * endpoint["R"]
            expected_moments = {
                "theta": endpoint["moments"]["theta"] + endpoint["F"] * delta_R2,
                "z": endpoint["moments"]["z"] + endpoint["Uz"] * delta_R,
                "theta_z": endpoint["moments"]["theta_z"]
                + endpoint["F"] * endpoint["Uz"] * delta_R2,
                "z_theta": endpoint["moments"]["z_theta"]
                + endpoint["Uz"] ** 2 * delta_R
                - endpoint["F"] ** 2 * delta_R2 / 2,
                "p": endpoint["moments"]["p"] + endpoint["F"] ** 2 * delta_R,
            }
            moment_error = max(
                relative(value["moments"][name], expected)
                for name, expected in expected_moments.items()
            )
            pressure_expected = endpoint["P"] + endpoint["F"] ** 2 * delta_R
            pressure_error = relative(value["P"], pressure_expected)
            max_driver_error = max(max_driver_error, driver_error)
            max_moment_error = max(max_moment_error, moment_error)
            max_pressure_error = max(max_pressure_error, pressure_error)
            rows.append(
                {
                    "R": mp.nstr(radius, 30),
                    "driver_relative_error": mp.nstr(driver_error, 30),
                    "moment_relative_error": mp.nstr(moment_error, 30),
                    "pressure_relative_error": mp.nstr(pressure_error, 30),
                }
            )
        tolerance = mp.mpf("1e-8")
        if max_driver_error > tolerance or max_moment_error > tolerance or max_pressure_error > tolerance:
            raise AssertionError(
                "frozen adapter consistency exceeded tolerance: "
                f"{max_driver_error}, {max_moment_error}, {max_pressure_error}"
            )
        result = {
            "branch": adapter.branch,
            "center": mp.nstr(adapter.center, 40),
            "Lambda": mp.nstr(adapter.Lambda, 40),
            "radial_degree": 124,
            "state_sha256": calc.state_hash,
            "rows": rows,
            "max_driver_relative_error": mp.nstr(max_driver_error, 40),
            "max_moment_relative_error": mp.nstr(max_moment_error, 40),
            "max_pressure_relative_error": mp.nstr(max_pressure_error, 40),
            "frozen_coefficients_reused": True,
            "amplitude_factor_cancelled_algebraically": True,
            "pressure_reset_or_fit": False,
            "nominal_center_only": True,
            "transition_discretization_enclosed": False,
            "whole_axis_certified": False,
            "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        }
    Path(__file__).with_suffix(".json").write_text(
        json.dumps(result, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(result, indent=2))
    return result


if __name__ == "__main__":
    run()
