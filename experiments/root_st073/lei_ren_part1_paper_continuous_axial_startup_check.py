"""Focused checks for the continuous axial-pulse startup endpoint.

The startup primitive is evaluated with the same ``sigma_pair`` source as
``ContinuousAxialPulse.value_jet``.  This diagnostic compares the transformed
positive primitive with direct quadrature where the latter is resolvable,
checks logarithmic derivatives with endpoint-scaled steps, and evaluates the
full pulse-energy atom.  It does not install coefficients or certify the
quadrature globally.
"""

from __future__ import annotations

import json
from pathlib import Path
import sys

import mpmath as mp


HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from lei_ren_part1_paper_axial_pulse import pulse_K_p  # noqa: E402
from lei_ren_part1_paper_continuous_axial_pulse import (  # noqa: E402
    ContinuousAxialPulse,
)
from lei_ren_part1_paper_continuous_pulse_energy import (  # noqa: E402
    _continuous_parts,
)


OUT_JSON = Path(__file__).with_suffix(".json")
ENERGY_JSON = HERE / "lei_ren_part1_paper_continuous_pulse_energy.json"


def _n(value: mp.mpf, digits: int = 80) -> str:
    return mp.nstr(mp.mpf(value), digits)


def _relative(left: mp.mpf, right: mp.mpf) -> mp.mpf:
    scale = max(abs(left), abs(right))
    return abs(left - right) / scale if scale else abs(left - right)


def _direct_startup_primitive(pulse: ContinuousAxialPulse, q: mp.mpf) -> mp.mpf:
    """Independent original-coordinate quadrature for moderate ``q``."""

    return mp.quad(lambda u: pulse.sigma_pair(u)[0], [0, q]) / 50


def _stored_legacy_kp() -> mp.mpf | None:
    if not ENERGY_JSON.exists():
        return None
    with ENERGY_JSON.open(encoding="utf-8") as handle:
        receipt = json.load(handle)
    value = receipt.get("legacy_float_comparison", {}).get("pulse_K_p_order_192")
    if isinstance(value, dict) and value.get("arbitrary_exponent_value") is not None:
        return mp.mpf(value["arbitrary_exponent_value"])
    return None


def run() -> dict[str, object]:
    precision = 150
    with mp.workdps(precision + 30):
        pulse = ContinuousAxialPulse(precision=precision)

        startup_rows = []
        for xi in (mp.mpf(".01"), mp.mpf(".001"), mp.mpf(".00001")):
            jet = pulse.value_jet(xi)
            h = xi**4
            left = pulse.value_jet(xi - h)["value"]
            right = pulse.value_jet(xi + h)["value"]
            numerical_log_derivative = (
                mp.log(right) - mp.log(left)
            ) / (2 * h)
            analytic_log_derivative = jet["derivative"] / jet["value"]
            startup_rows.append(
                {
                    "xi": _n(xi),
                    "step": _n(h),
                    "value": _n(jet["value"]),
                    "log10_value": _n(mp.log10(jet["value"])),
                    "derivative": _n(jet["derivative"]),
                    "positive": bool(jet["value"] > 0),
                    "log_derivative": _n(analytic_log_derivative),
                    "log_derivative_finite_difference": _n(
                        numerical_log_derivative
                    ),
                    "log_derivative_relative_error": _n(
                        abs(numerical_log_derivative / analytic_log_derivative - 1)
                    ),
                }
            )

        baseline_rows = []
        for xi in (mp.mpf(".005"), mp.mpf(".01")):
            q = 50 * xi
            direct = _direct_startup_primitive(pulse, q)
            transformed = pulse._startup_primitive(q)
            baseline_rows.append(
                {
                    "xi": _n(xi),
                    "q": _n(q),
                    "direct_primitive": _n(direct),
                    "transformed_primitive": _n(transformed),
                    "relative_error": _n(_relative(direct, transformed)),
                }
            )

        continuity_rows = []
        epsilon = mp.mpf("1e-20")
        for boundary in (mp.mpf(".01"), mp.mpf(".02")):
            left = pulse.value_jet(boundary - epsilon)["value"]
            center = pulse.value_jet(boundary)["value"]
            right = pulse.value_jet(boundary + epsilon)["value"]
            continuity_rows.append(
                {
                    "boundary": _n(boundary),
                    "epsilon": _n(epsilon),
                    "left": _n(left),
                    "center": _n(center),
                    "right": _n(right),
                    "left_right_relative_error": _n(_relative(left, right)),
                    "left_center_relative_error": _n(_relative(left, center)),
                    "right_center_relative_error": _n(_relative(right, center)),
                }
            )

        support_zero = {
            "xi_zero": pulse.value_jet(mp.mpf("0"))["value"] == 0,
            "xi_eleven": pulse.value_jet(mp.mpf("11"))["value"] == 0,
        }

        parts = _continuous_parts(precision=60, quadrature_order=48)
        legacy_kp = _stored_legacy_kp()
        if legacy_kp is None:
            legacy_kp = mp.mpf(str(pulse_K_p(quadrature_order=192)))
        energy = {
            "startup": _n(parts["startup"]),
            "plateau": _n(parts["plateau"]),
            "cutoff": _n(parts["cutoff"]),
            "total": _n(parts["total"]),
            "stored_legacy_Kp": _n(legacy_kp),
            "difference_from_stored_legacy_Kp": _n(parts["total"] - legacy_kp),
            "relative_difference_from_stored_legacy_Kp": _n(
                _relative(parts["total"], legacy_kp)
            ),
            "working_precision": 60,
            "quadrature_order": 48,
            "comparison_scope": (
                "The stored Kp is a legacy float receipt; this is a nominal MP "
                "comparison and not an enclosure."
            ),
        }

        derivative_ok = all(
            mp.mpf(row["log_derivative_relative_error"]) < mp.mpf("1e-10")
            for row in startup_rows
        )
        baseline_ok = all(
            mp.mpf(row["relative_error"]) < mp.mpf("1e-90")
            for row in baseline_rows
        )
        continuity_ok = all(
            mp.mpf(row["left_right_relative_error"]) < mp.mpf("1e-15")
            for row in continuity_rows
        )
        report: dict[str, object] = {
            "kind": "continuous_axial_startup_check",
            "source": "ContinuousAxialPulse.value_jet",
            "startup_rows": startup_rows,
            "independent_direct_quadrature_rows": baseline_rows,
            "continuity_rows": continuity_rows,
            "support_zero": support_zero,
            "energy": energy,
            "checks": {
                "startup_values_positive": all(row["positive"] for row in startup_rows),
                "log_derivative_checks_pass": derivative_ok,
                "direct_quadrature_checks_pass": baseline_ok,
                "continuity_checks_pass": continuity_ok,
                "support_endpoint_checks_pass": all(support_zero.values()),
                "energy_evaluated": mp.isfinite(parts["total"]),
            },
            "scope": (
                "Endpoint-scaled positive startup primitive using the shared sigma "
                "definition. Numerical quadrature error, coefficient closure, "
                "global installation, and finite-energy certification remain open."
            ),
        }
    OUT_JSON.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    return report


if __name__ == "__main__":
    result = run()
    print(
        json.dumps(
            {
                "checks": result["checks"],
                "energy_total": result["energy"]["total"],
                "energy_relative_difference_from_stored_legacy_Kp": result[
                    "energy"
                ]["relative_difference_from_stored_legacy_Kp"],
            },
            indent=2,
        )
    )
