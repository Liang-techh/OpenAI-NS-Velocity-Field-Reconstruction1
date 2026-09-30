"""Independent resolved and extreme-scale checks for higher flat-tail orders.

The resolved rows compare each ``B`` derivative against direct ``mp.quad``
of its differentiated integrand.  The extreme rows only compare against the
leading flat-tail saddle formula; their absolute log discrepancy is expected
to be large because the leading log itself is of order ``10**101``.  No row
is an enclosure or a global field claim.
"""

from __future__ import annotations

import importlib
import json
from pathlib import Path
from typing import Any

import mpmath as mp

try:
    from .lei_ren_part1_paper_axial_primitive import _sigma_mp
except (ImportError, ValueError):
    from lei_ren_part1_paper_axial_primitive import _sigma_mp


PRECISION = 100
T_RESOLVED = mp.mpf(400)
K_RESOLVED = mp.mpf(".2")
M_RESOLVED = mp.mpf(2)
B_RESOLVED = mp.mpf(-4)
ORDERS = (1, 2, 3)
CUTS = tuple(mp.mpf(value) for value in (0, 10, 25, 50, 100, 200, 300, 400))
ERROR_THRESHOLD = mp.mpf("1e-14")
EXTREME_SCALED_THRESHOLD = mp.mpf("1e-50")


def _quad_split(function: Any) -> mp.mpf:
    return mp.fsum(
        mp.quad(function, [left, right])
        for left, right in zip(CUTS, CUTS[1:])
    )


def _direct(order: int) -> mp.mpf:
    def integrand(ell: mp.mpf) -> mp.mpf:
        q = _sigma_mp(ell / T_RESOLVED)
        return (
            mp.exp(-K_RESOLVED * ell)
            * (M_RESOLVED * q) ** order
            * mp.exp(M_RESOLVED * B_RESOLVED * q)
        )

    return _quad_split(integrand)


def _leading_extreme(order: int) -> tuple[mp.mpf, mp.mpf]:
    k = mp.mpf(".2")
    m = mp.mpf(2)
    T = mp.mpf("4e152")
    mode = (2 * order * T * T / k) ** (mp.mpf(1) / 3)
    log_integral = (
        -mp.mpf("1.5") * k * mode
        + order * (1 + mp.log(m))
        + mp.mpf(".5") * mp.log(2 * mp.pi * mode / (3 * k))
    )
    return mode, log_integral


def run_fixture(*, precision: int = PRECISION) -> dict[str, Any]:
    precision = max(80, int(precision))
    with mp.workdps(precision):
        module = importlib.import_module(
            "lei_ren_part1_paper_flat_shape_defect"
        )
        rows: list[dict[str, Any]] = []
        maximum_resolved_log_error = mp.mpf(0)
        maximum_resolved_relative_error = mp.mpf(0)
        for order in ORDERS:
            expected = _direct(order)
            result = module.evaluate_flat_shape_derivative(
                K_RESOLVED,
                M_RESOLVED,
                B_RESOLVED,
                T_RESOLVED,
                order,
                precision=precision,
                order=32,
                window=24,
            )
            log_error = abs(result["log_abs"] - mp.log(expected))
            relative_error = abs(result["value"] - expected) / expected
            maximum_resolved_log_error = max(
                maximum_resolved_log_error, log_error
            )
            maximum_resolved_relative_error = max(
                maximum_resolved_relative_error, relative_error
            )
            rows.append(
                {
                    "kind": "resolved",
                    "order": order,
                    "expected": mp.nstr(expected, 35),
                    "actual": mp.nstr(result["value"], 35),
                    "log_error": mp.nstr(log_error, 25),
                    "relative_error": mp.nstr(relative_error, 25),
                    "mode": mp.nstr(result["mode"], 25),
                    "width": mp.nstr(result["width"], 25),
                    "positive_sign": result["sign"] == 1,
                }
            )

        if maximum_resolved_log_error >= ERROR_THRESHOLD:
            raise AssertionError(
                f"resolved log derivative mismatch: {maximum_resolved_log_error}"
            )
        if maximum_resolved_relative_error >= ERROR_THRESHOLD:
            raise AssertionError(
                f"resolved derivative mismatch: {maximum_resolved_relative_error}"
            )

        extreme_rows: list[dict[str, Any]] = []
        maximum_extreme_mode_error = mp.mpf(0)
        maximum_extreme_scaled_log_error = mp.mpf(0)
        for order in ORDERS:
            mode_expected, log_expected = _leading_extreme(order)
            result = module.evaluate_flat_shape_derivative(
                ".2",
                "2",
                "-8.62e36",
                "4e152",
                order,
                precision=260,
                order=32,
                window=24,
            )
            mode_error = abs(result["mode"] / mode_expected - 1)
            scaled_log_error = abs(result["log_abs"] - log_expected) / max(
                mp.mpf(1), abs(log_expected)
            )
            maximum_extreme_mode_error = max(maximum_extreme_mode_error, mode_error)
            maximum_extreme_scaled_log_error = max(
                maximum_extreme_scaled_log_error, scaled_log_error
            )
            extreme_rows.append(
                {
                    "kind": "extreme_leading_saddle",
                    "order": order,
                    "mode": mp.nstr(result["mode"], 35),
                    "leading_mode": mp.nstr(mode_expected, 35),
                    "mode_relative_error": mp.nstr(mode_error, 25),
                    "scaled_log_error": mp.nstr(scaled_log_error, 25),
                    "absolute_log_error": mp.nstr(
                        abs(result["log_abs"] - log_expected), 25
                    ),
                }
            )

        if maximum_extreme_mode_error >= EXTREME_SCALED_THRESHOLD:
            raise AssertionError(
                f"extreme saddle mode mismatch: {maximum_extreme_mode_error}"
            )
        if maximum_extreme_scaled_log_error >= EXTREME_SCALED_THRESHOLD:
            raise AssertionError(
                "extreme leading-log mismatch: "
                f"{maximum_extreme_scaled_log_error}"
            )

        report = {
            "passed": True,
            "resolved_precision": precision,
            "resolved_T": mp.nstr(T_RESOLVED, 25),
            "resolved_k": mp.nstr(K_RESOLVED, 25),
            "resolved_m": mp.nstr(M_RESOLVED, 25),
            "resolved_B": mp.nstr(B_RESOLVED, 25),
            "orders": list(ORDERS),
            "cuts": [mp.nstr(value, 20) for value in CUTS],
            "maximum_resolved_log_error": mp.nstr(
                maximum_resolved_log_error, 25
            ),
            "maximum_resolved_relative_error": mp.nstr(
                maximum_resolved_relative_error, 25
            ),
            "maximum_extreme_mode_relative_error": mp.nstr(
                maximum_extreme_mode_error, 25
            ),
            "maximum_extreme_scaled_log_error": mp.nstr(
                maximum_extreme_scaled_log_error, 25
            ),
            "resolved_rows": rows,
            "extreme_rows": extreme_rows,
            "leading_mode_formula": "(2*n*T^2/k)^(1/3)",
            "leading_log_formula": (
                "-1.5*k*mode + n*(1+log(m)) "
                "+ .5*log(2*pi*mode/(3*k))"
            ),
            "independent_quadrature": "mp.quad over resolved cuts",
            "quadrature_error_enclosed": False,
            "saddle_remainder_enclosed": False,
            "global_field_installed": False,
        }
        Path(__file__).with_suffix(".json").write_text(
            json.dumps(report, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        return report


def fixture(**kwargs: Any) -> dict[str, Any]:
    return run_fixture(**kwargs)


def main() -> None:
    print(json.dumps(run_fixture(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()


__all__ = ["fixture", "run_fixture"]
