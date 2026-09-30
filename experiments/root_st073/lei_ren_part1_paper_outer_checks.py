"""Focused checks for the source-pinned Section 6 outer schedule.

The checks exercise the explicit reference and terminal schedule only.  They
do not certify pressure, moments, PDE residuals, an admissible cone, or the
paper's scale recursion.
"""

from __future__ import annotations

from decimal import Decimal, localcontext
import json
import math
from pathlib import Path
import sys
from typing import Any


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SRC = ROOT / "src"
for path in (SRC, HERE):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from lei_ren_part1_paper_outer import (  # noqa: E402
    PaperOuterSchedule,
    SOURCE,
    SOURCE_SECTION,
    SOURCE_VERSION,
)


DEMO_INPUTS = {
    "logPstar": "0",
    "logRref": "-230",
    "delta": "0.002",
    "Md": "0.5",
    "c_mu": "0.5",
    "c_delta": "0.25",
    "c_epsilon": "0.01",
}

# This second input set keeps the theorem-style smallness flags visible while
# still allowing the check to probe a genuinely enormous Decimal checkpoint.
SOURCE_INPUTS = {
    "logPstar": "74",
    "logRref": "10",
    "delta": "1e-32",
    "Md": "0.5",
    "c_mu": "0.001",
    "c_delta": "0.001",
    "c_epsilon": "0.01",
}


def _float(value: Any) -> float:
    return float(value)


def _close(value: Any, target: Any, tolerance: float) -> tuple[bool, float]:
    error = abs(_float(value) - _float(target))
    return error <= tolerance, error


def _finite_schedule_validation() -> dict[str, Any]:
    schedule = PaperOuterSchedule(**DEMO_INPUTS)
    failures: list[str] = []
    for name, call in (
        ("negative_R", lambda: schedule.at_radius(-1, 0.0)),
        ("nan_R", lambda: schedule.at_radius(float("nan"), 0.0)),
        ("out_of_range_Z", lambda: schedule.at_log_radius(0, 1.01)),
        (
            "negative_waiting_length",
            lambda: PaperOuterSchedule(**DEMO_INPUTS, waiting_length="-1"),
        ),
    ):
        try:
            call()
        except (TypeError, ValueError, OverflowError):
            continue
        failures.append(name)
    return {"passed": not failures, "unexpectedly_accepted": failures}


def run() -> dict[str, Any]:
    schedule = PaperOuterSchedule(**DEMO_INPUTS)
    source_schedule = PaperOuterSchedule(**SOURCE_INPUTS)
    checks: dict[str, Any] = {}

    # Reference branch: U_theta is R^0.1 and F is R^-0.4.
    reference_rows = [
        schedule.at_log_radius(schedule.logRref + Decimal(str(y)), 0.3)
        for y in (-4, -2, 0)
    ]
    u_ok = max(abs(_float(row["logarithmic_slope"]) - 0.1) for row in reference_rows)
    f_ok = max(abs(_float(row["logF_slope"]) + 0.4) for row in reference_rows)
    checks["reference_power_laws"] = {
        "max_Utheta_log_slope_error": u_ok,
        "max_F_log_slope_error": f_ok,
        "passed": u_ok < 1e-14 and f_ok < 1e-14,
    }

    # Source slope schedule at exact checkpoints.  At y=1 the reference
    # transition is complete; y_d is already beyond the axial turnoff; y_f is
    # after the z-flattening stage; and the terminal branch has slope
    # -(1+delta)/2 before its heat factor radial derivative is introduced.
    stage_points = {
        "slope_transition_ref": (Decimal(1), -0.5),
        "axial_turnoff": (schedule.y_d + Decimal(1), -0.5 - float(schedule.mu)),
        "z_flatten_end": (schedule.y_f, -0.5 - float(schedule.mu)),
        "terminal_start": (schedule.y_tail, -(1 + float(schedule.delta)) / 2),
    }
    stage_rows: dict[str, Any] = {}
    stage_passed = True
    for name, (y, expected) in stage_points.items():
        row = schedule.at_log_radius(schedule.logRref + y, 0.0)
        ok, error = _close(row["logarithmic_slope"], expected, 2e-12)
        stage_rows[name] = {
            "stage": row["stage"],
            "observed_slope": str(row["logarithmic_slope"]),
            "expected_slope": expected,
            "absolute_error": error,
            "passed": ok,
        }
        stage_passed = stage_passed and ok
    checks["stage_slopes"] = {"rows": stage_rows, "passed": stage_passed}

    exact_heat_y = schedule.y_b + Decimal(1)
    exact_heat_log_radius = schedule.logRref + exact_heat_y
    exact_heat_row = schedule.at_log_radius(exact_heat_log_radius, 0.0)
    exact_heat = schedule._heat_info(exact_heat_log_radius, 0.0)
    exact_heat_expected = (
        -(1 + float(schedule.delta)) / 2
        - float(exact_heat["xi"]) * float(exact_heat["H_prime"])
        / float(exact_heat["H"])
    )
    exact_heat_ok, exact_heat_error = _close(
        exact_heat_row["logarithmic_slope"], exact_heat_expected, 2e-12
    )
    stage_rows["exact_heat"] = {
        "stage": exact_heat_row["stage"],
        "observed_slope": str(exact_heat_row["logarithmic_slope"]),
        "expected_slope": exact_heat_expected,
        "absolute_error": exact_heat_error,
        "passed": exact_heat_ok,
    }
    checks["stage_slopes"]["passed"] = stage_passed and exact_heat_ok

    # The paper replaces (1+Z^2)^-1 by the constant 1/2 after the flat
    # z-transition.  Check both the value factor and its Z derivative.
    flatten_rows = []
    flatten_passed = True
    log_a = schedule._log_A(schedule.y_f)
    for z in (0.0, 0.5, 1.0):
        row = schedule.at_log_radius(schedule.logRref + schedule.y_f, z)
        factor = math.exp(float(row["log_angular_amplitude"] - log_a))
        factor_error = abs(factor - 0.5)
        derivative = abs(float(row["dlogU_dZ"]))
        row_ok = factor_error < 1e-14 and derivative < 1e-14
        flatten_rows.append(
            {
                "Z": z,
                "factor": factor,
                "factor_error": factor_error,
                "abs_dlogU_dZ": derivative,
                "passed": row_ok,
            }
        )
        flatten_passed = flatten_passed and row_ok
    checks["z_flattening"] = {"rows": flatten_rows, "passed": flatten_passed}

    # Uz=4Z through the reference/early schedule and zero after the cutoff.
    uz_before = schedule.at_log_radius(schedule.logRref + Decimal(1), 0.4)["Uz"]
    uz_after = schedule.at_log_radius(
        schedule.logRref + schedule.y_d + Decimal(1), 0.4
    )["Uz"]
    checks["axial_turnoff"] = {
        "Uz_before": uz_before,
        "Uz_expected_before": 1.6,
        "Uz_after": uz_after,
        "passed": abs(uz_before - 1.6) < 1e-14 and abs(uz_after) < 1e-14,
    }

    # At y_tail the c_inf normalization makes the pre-tail value A/2.  At y_b
    # the interpolation has become the exact heat factor; compare both value
    # and d/d log R slope against that factor.
    z = 0.3
    tail_row = schedule.at_log_radius(schedule.logRref + schedule.y_tail, z)
    tail_expected = schedule._log_A(schedule.y_tail) - Decimal(str(math.log(2.0)))
    tail_error = abs(float(tail_row["log_angular_amplitude"] - tail_expected))
    join_log_radius = schedule.logRref + schedule.y_b
    join_row = schedule.at_log_radius(join_log_radius, z)
    heat = schedule._heat_info(join_log_radius, z)
    a = (Decimal(1) + schedule.delta) / Decimal(2)
    exact_heat_log = (
        schedule._log_c_inf
        - a * join_log_radius
        + Decimal(str(math.log(float(heat["H"]))))
    )
    join_value_error = abs(float(join_row["log_angular_amplitude"] - exact_heat_log))
    join_slope_expected = -float(a) - float(heat["xi"]) * float(heat["H_prime"]) / float(heat["H"])
    join_slope_error = abs(float(join_row["logarithmic_slope"]) - join_slope_expected)
    checks["heat_join"] = {
        "tail_value_log_error": tail_error,
        "R_b_value_log_error": join_value_error,
        "R_b_slope_error": join_slope_error,
        "heat_method": join_row["heat_method"],
        "declared_heat_limit_bound": join_row["heat_limit_bound"],
        "passed": tail_error < 1e-12 and join_value_error < 1e-12 and join_slope_error < 1e-12,
    }

    # Waiting only moves the terminal checkpoint through a constant-slope
    # segment, so c_inf is invariant.  This is a schedule identity, not a
    # completed pressure matching result.
    waiting_schedules = [
        PaperOuterSchedule(**DEMO_INPUTS, waiting_length=str(waiting))
        for waiting in (0, 0.7, 10)
    ]
    log_c_values = [str(item._log_c_inf) for item in waiting_schedules]
    checks["waiting_length_invariance"] = {
        "log_c_inf_values": log_c_values,
        "max_absolute_log_difference": max(
            abs(float(item._log_c_inf - waiting_schedules[0]._log_c_inf))
            for item in waiting_schedules
        ),
        "passed": all(item._log_c_inf == waiting_schedules[0]._log_c_inf for item in waiting_schedules),
    }

    # Large-log-Pstar schedule: checkpoint offsets remain exact Decimal
    # differences even when the physical radii cannot be represented as float.
    with localcontext() as context:
        context.prec = source_schedule.decimal_precision + 16
        yf_difference = source_schedule.y_f - source_schedule.y_v
        yb_difference = source_schedule.y_b - source_schedule.y_tail
    checks["large_schedule_precision"] = {
        "decimal_precision": source_schedule.decimal_precision,
        "y_v_adjusted_exponent": source_schedule.y_v.adjusted(),
        "y_f_minus_y_v": str(yf_difference),
        "y_b_minus_y_tail": str(yb_difference),
        "c_inf_numeric_available": source_schedule.c_inf_decimal is not None,
        "passed": yf_difference == Decimal(100) and yb_difference == Decimal(3),
    }

    checks["finite_input_validation"] = _finite_schedule_validation()

    passed = all(value["passed"] for value in checks.values())
    report = {
        "source": SOURCE,
        "source_version": SOURCE_VERSION,
        "source_section": SOURCE_SECTION,
        "api": "PaperOuterSchedule + PaperOuterProfile",
        "demo_inputs": DEMO_INPUTS,
        "source_precision_inputs": SOURCE_INPUTS,
        "schedule_metadata": schedule.metadata(),
        "source_precision_metadata": source_schedule.metadata(),
        "checks": checks,
        "checks_passed": passed,
        "scope": (
            "Section 6 explicit source outer velocity schedule and Decimal log-radius "
            "checkpoint arithmetic only. No pressure, five-moment closure, PDE, "
            "cone, admissibility, or scale-recursion certification."
        ),
        "uncertified_assumptions": [
            "Absolute source constants are supplied numerical inputs.",
            "The moderate demo intentionally does not satisfy all theorem smallness flags.",
            "H=1 is used only with an explicit bound when a heat argument is outside float range.",
        ],
    }
    return report


if __name__ == "__main__":
    output = run()
    output_path = Path(__file__).with_suffix(".json")
    output_path.write_text(json.dumps(output, indent=2, default=str) + "\n", encoding="utf-8")
    print(json.dumps({"checks_passed": output["checks_passed"], "checks": output["checks"]}, indent=2, default=str))
    if not output["checks_passed"]:
        raise SystemExit(1)
