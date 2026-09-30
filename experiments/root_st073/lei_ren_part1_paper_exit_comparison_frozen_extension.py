"""Focused receipt for the frozen Section 9.23 comparison extension.

This runner deliberately keeps the finite-core receipt untouched.  It checks
the exact frozen moment propagation and the closed D/I_z power formulas at
physical radii 1, 100, and 110.  The zero radial shear in the frozen branch
is recorded as a strict source-cone failure; no cone certification is made.
"""

from __future__ import annotations

import json
from pathlib import Path
import sys
from typing import Any

import mpmath as mp


HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from lei_ren_part1_paper_exit_comparison import (  # noqa: E402
    MOMENT_KEYS,
    _nstr,
    _relative,
    _signed_log,
    build_comparison,
)


def _moment_expectation(endpoint: dict[str, Any], R: mp.mpf) -> tuple[dict[str, Any], dict[str, Any]]:
    rb = endpoint["R"]
    dr = R - rb
    dr2 = R * R - rb * rb
    f = endpoint["F"]
    fz = endpoint["FZ"]
    u = endpoint["Uz"]
    uz = endpoint["UZ"]
    moments = endpoint["moments"]
    moments_z = endpoint["momentsZ"]
    expected = {
        "theta": moments["theta"] + f * dr2,
        "z": moments["z"] + u * dr,
        "theta_z": moments["theta_z"] + f * u * dr2,
        "z_theta": moments["z_theta"] + u * u * dr - f * f * dr2 / 2,
        "p": moments["p"] + f * f * dr,
    }
    expected_z = {
        "theta": moments_z["theta"] + fz * dr2,
        "z": moments_z["z"] + uz * dr,
        "theta_z": moments_z["theta_z"] + (fz * u + f * uz) * dr2,
        "z_theta": moments_z["z_theta"] + 2 * u * uz * dr - f * fz * dr2,
        "p": moments_z["p"] + 2 * f * fz * dr,
    }
    return expected, expected_z


def run_checks() -> dict[str, Any]:
    with mp.workdps(160):
        comparison = build_comparison(precision=160, degree=18, Lambda="1e36", h_b=".005")
        Z = ".3"
        endpoint = comparison.evaluate(2 * comparison.hb, Z)
        coefficients = comparison.frozen_driver_coefficients(Z)
        rows = []
        max_moment_error = mp.mpf(0)
        max_moment_z_error = mp.mpf(0)
        max_D_error = mp.mpf(0)
        max_I_z_error = mp.mpf(0)
        cone_failures: dict[str, int] = {
            "D_not_positive": 0,
            "Utheta_not_positive": 0,
            "S_theta_not_strictly_negative": 0,
            "kappa_not_defined": 0,
        }
        for R_text in ("1", "100", "110"):
            R = mp.mpf(R_text)
            y = mp.log(R / comparison.R_a)
            value = comparison.evaluate(y, Z)
            expected, expected_z = _moment_expectation(endpoint, R)
            moment_error = max(
                _relative(value["moments"][key], expected[key]) for key in MOMENT_KEYS
            )
            moment_z_error = max(
                _relative(value["momentsZ"][key], expected_z[key]) for key in MOMENT_KEYS
            )
            driver_D = (
                coefficients["D"][1] * R
                + coefficients["D"][0]
                + coefficients["D"][-1] / R
            )
            driver_I_z = (
                coefficients["I_z"][3] * R ** mp.mpf("1.5")
                + coefficients["I_z"][1] * R ** mp.mpf(".5")
                + coefficients["I_z"][-1] * R ** mp.mpf("-.5")
            )
            D_error = _relative(value["D"], driver_D)
            I_z_error = _relative(value["I_z"], driver_I_z)
            max_moment_error = max(max_moment_error, moment_error)
            max_moment_z_error = max(max_moment_z_error, moment_z_error)
            max_D_error = max(max_D_error, D_error)
            max_I_z_error = max(max_I_z_error, I_z_error)
            if value["D"] <= 0:
                cone_failures["D_not_positive"] += 1
            if value["F"] <= 0:
                cone_failures["Utheta_not_positive"] += 1
            if value["S_theta"] >= 0:
                cone_failures["S_theta_not_strictly_negative"] += 1
            if value["S_theta"] == 0:
                cone_failures["kappa_not_defined"] += 1
            rows.append(
                {
                    "R": R_text,
                    "y": _nstr(y, 40),
                    "core_scope": value["core_scope"],
                    "core_reference_R": _nstr(value["core_reference_R"], 40),
                    "F": _signed_log(value["F"]),
                    "Uz": _signed_log(value["Uz"]),
                    "D": _signed_log(value["D"]),
                    "E": _signed_log(value["E"]),
                    "I_z": _signed_log(value["I_z"]),
                    "S_theta": _signed_log(value["S_theta"]),
                    "T_theta": _signed_log(value["T_theta"]),
                    "T_z": _signed_log(value["T_z"]),
                    "FR_zero": value["FR"] == 0,
                    "UR_zero": value["UR"] == 0,
                    "moment_relative_error": _nstr(moment_error, 30),
                    "moment_Z_relative_error": _nstr(moment_z_error, 30),
                    "D_driver_relative_error": _nstr(D_error, 30),
                    "I_z_driver_relative_error": _nstr(I_z_error, 30),
                    "D_positive": value["D"] > 0,
                    "S_theta_strictly_negative": value["S_theta"] < 0,
                }
            )
        return {
            "metadata": {
                "source": "https://arxiv.org/html/2609.35406v1",
                "source_version": "2609.35406v1",
                "equation": "Section 9.23 / source (9.23)",
                "profile": "Section923Comparison frozen endpoint continuation",
                "Lambda": _nstr(comparison.Lambda, 50),
                "radial_degree": comparison.bundle.get("radial_degree"),
                "h_b": _nstr(comparison.hb, 50),
                "R_a": _nstr(comparison.R_a, 50),
                "R_b": _nstr(coefficients["R_b"], 50),
                "Z": Z,
                "precision": comparison.precision,
                "core_scope": "core queried only at R_b; all R > R_b use endpoint frozen fields",
                "D_formula": "D(R)=D[1] R+D[0]+D[-1] R^(-1)",
                "I_z_formula": "I_z(R)=I_z[3] R^(3/2)+I_z[1] R^(1/2)+I_z[-1] R^(-1/2)",
                "I_z_coefficients": {
                    str(key): _nstr(value, 50) for key, value in coefficients["I_z"].items()
                },
                "D_coefficients": {
                    str(key): _nstr(value, 50) for key, value in coefficients["D"].items()
                },
            },
            "rows": rows,
            "max_frozen_moment_relative_error": _nstr(max_moment_error, 30),
            "max_frozen_moment_Z_relative_error": _nstr(max_moment_z_error, 30),
            "max_D_driver_relative_error": _nstr(max_D_error, 30),
            "max_I_z_driver_relative_error": _nstr(max_I_z_error, 30),
            "source_cone_raw_failures": cone_failures,
            "source_cone_certified": False,
            "scope": (
                "Frozen comparison stress and moment identities only. "
                "The frozen branch has zero radial shear, so strict source-cone "
                "conditions are recorded as failures rather than certified."
            ),
        }


if __name__ == "__main__":
    receipt = run_checks()
    output = Path(__file__).with_suffix(".json")
    output.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "rows": len(receipt["rows"]),
                "max_D_driver_relative_error": receipt["max_D_driver_relative_error"],
                "max_I_z_driver_relative_error": receipt["max_I_z_driver_relative_error"],
                "source_cone_raw_failures": receipt["source_cone_raw_failures"],
            }
        ),
        flush=True,
    )
