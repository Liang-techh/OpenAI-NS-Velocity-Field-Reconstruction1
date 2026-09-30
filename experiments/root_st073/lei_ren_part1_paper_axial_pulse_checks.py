"""Focused checks for the fixed Part I Section 7.5 axial pulse."""

from __future__ import annotations

import json
import math
from pathlib import Path
import sys

import numpy as np
from scipy.integrate import quad


ROOT = Path(__file__).resolve().parents[2]
LOCAL = Path(__file__).resolve().parent
if str(LOCAL) not in sys.path:
    sys.path.insert(0, str(LOCAL))

from lei_ren_part1_paper_axial_pulse import (  # noqa: E402
    SOURCE,
    SOURCE_EQUATION,
    SOURCE_KP_BOUND,
    SOURCE_SECTION,
    SOURCE_VERSION,
    PaperAxialPulse,
    pulse_derivative,
    pulse_K_p,
    pulse_primitive,
    pulse_primitive_derivative,
    pulse_value,
    smooth_step,
)


REPORT_PATH = LOCAL / "lei_ren_part1_paper_axial_pulse_checks.json"


def _fd(function, x: float, step: float) -> float:
    return (function(x + step) - function(x - step)) / (2.0 * step)


def _piecewise_quad() -> tuple[float, float]:
    intervals = ((0.0, 0.02), (0.02, 10.0), (10.0, 10.02), (10.02, 11.0))
    values = []
    errors = []
    for left, right in intervals:
        value, error = quad(
            lambda xi: math.exp(-2.0 * xi) * pulse_value(xi) ** 2,
            left,
            right,
            epsabs=2.0e-13,
            epsrel=2.0e-13,
            limit=300,
        )
        values.append(value)
        errors.append(error)
    return float(sum(values)), float(sum(errors))


def main() -> int:
    pulse = PaperAxialPulse()
    kp64 = pulse_K_p(quadrature_order=64)
    kp128 = pulse_K_p(quadrature_order=128)
    kp256 = pulse_K_p(quadrature_order=256)
    kp_quad, kp_quad_error = _piecewise_quad()

    primitive_points = (0.001, 0.01, 0.019, 0.03, 0.5)
    primitive_rows = []
    for point in primitive_points:
        step = min(1.0e-5, point / 10.0)
        observed = _fd(pulse_primitive, point, step)
        expected = pulse_primitive_derivative(point)
        primitive_rows.append(
            {
                "xi": point,
                "finite_difference": observed,
                "sigma_50xi": expected,
                "absolute_error": abs(observed - expected),
            }
        )

    derivative_points = (0.005, 0.015, 0.5, 9.5, 10.005, 10.5)
    derivative_rows = []
    for point in derivative_points:
        step = 1.0e-5
        observed = _fd(pulse_value, point, step)
        expected = pulse_derivative(point)
        derivative_rows.append(
            {
                "xi": point,
                "finite_difference": observed,
                "analytic_derivative": expected,
                "absolute_error": abs(observed - expected),
            }
        )

    support_checks = {
        "value_at_negative": pulse_value(-1.0),
        "value_at_zero": pulse_value(0.0),
        "value_at_eleven": pulse_value(11.0),
        "value_at_thirteen": pulse_value(13.0),
        "derivative_at_negative": pulse_derivative(-1.0),
        "derivative_at_eleven": pulse_derivative(11.0),
    }
    stage_checks = {
        "mu": 1.0e-28,
        "xi": 5.0,
        "value_at_stage": pulse.value_at_stage(5.0 / 1.0e-28, 1.0e-28),
        "value_direct": pulse.value(5.0),
        "derivative_at_stage": pulse.derivative_at_stage(5.0 / 1.0e-28, 1.0e-28),
        "expected_stage_derivative": 1.0e-28 * pulse.derivative(5.0),
    }

    max_primitive_error = max(row["absolute_error"] for row in primitive_rows)
    max_derivative_error = max(row["absolute_error"] for row in derivative_rows)
    kp_refinement = max(abs(kp64 - kp128), abs(kp128 - kp256), abs(kp256 - kp_quad))
    bound_passed = SOURCE_KP_BOUND[0] < kp_quad < SOURCE_KP_BOUND[1]
    report = {
        "source": SOURCE,
        "source_version": SOURCE_VERSION,
        "source_section": SOURCE_SECTION,
        "source_equation": SOURCE_EQUATION,
        "support": [0.0, 11.0],
        "K_p": {
            "gauss_legendre_64": kp64,
            "gauss_legendre_128": kp128,
            "gauss_legendre_256": kp256,
            "piecewise_scipy_quad": kp_quad,
            "piecewise_quad_error_estimate": kp_quad_error,
            "max_refinement_difference": kp_refinement,
            "source_bound": list(SOURCE_KP_BOUND),
            "source_bound_passed": bound_passed,
        },
        "primitive_derivative": {
            "rows": primitive_rows,
            "max_absolute_error": max_primitive_error,
        },
        "pulse_derivative": {
            "rows": derivative_rows,
            "max_absolute_error": max_derivative_error,
        },
        "support_checks": support_checks,
        "stage_local_checks": stage_checks,
        "scope": {
            "fixed_pulse_only": True,
            "actual_affine_rhs_bound": False,
            "full_axial_closure": False,
            "corrected_angular_energy_input_open": True,
            "absolute_large_radius_coordinates_constructed": False,
            "open_dependencies": [
                "Compute normalized incoming moment and pulse RHS data from the actual profile.",
                "Supply the corrected angular tail to the (7.34) energy target.",
            ],
        },
    }
    report["checks_passed"] = bool(
        bound_passed
        and max_primitive_error < 1.0e-7
        # The largest absolute derivative is about 83 on the terminal
        # transition; the observed 6.5e-7 absolute FD error is below 1e-8
        # relative there and is retained in the receipt.
        and max_derivative_error < 1.0e-5
        and max(abs(support_checks[key]) for key in ("value_at_negative", "value_at_zero", "value_at_eleven", "value_at_thirteen")) < 1.0e-14
        and abs(stage_checks["value_at_stage"] - stage_checks["value_direct"]) < 1.0e-14
        and abs(stage_checks["derivative_at_stage"] - stage_checks["expected_stage_derivative"]) < 1.0e-35
    )
    REPORT_PATH.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "report": str(REPORT_PATH),
        "checks_passed": report["checks_passed"],
        "K_p": kp_quad,
        "K_p_refinement": kp_refinement,
        "max_primitive_derivative_error": max_primitive_error,
        "max_pulse_derivative_error": max_derivative_error,
    }, indent=2))
    return 0 if report["checks_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
