"""Directed accepted-schedule axis-pressure Taylor jets at ``Z = 0.3``.

The receipt combines the accepted positive true-mass bounds into the analytic
fixed-beta terms and a Cauchy enclosure for the variable-beta flatten term.
It deliberately contains no finite pressure-parameter truncation and makes no
claim about omitted source or parameter errors.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext

from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode


HERE = Path(__file__).resolve().parent
FIXED_NAME = "lei_ren_part1_paper_coherent_uniform_fixed_beta_error.json"
COMPLETE_NAME = "lei_ren_part1_paper_coherent_uniform_preheat_error_budget.json"
OLD_AXIS_NAME = "lei_ren_part1_paper_axis_pressure_jets.json"
DEFAULT_LENGTH = 160
DEFAULT_PRECISION = 260
CENTER_Z = "0.3"
CAUCHY_RADIUS = "0.25"


def _interval_from_exact(ctx: MPIntervalContext, row: dict[str, Any]):
    """Restore an encoded interval with directed MP arithmetic."""

    return ctx.mpf(
        [
            mp.make_mpf(tuple(row["lower_exact_mpf_tuple"])),
            mp.make_mpf(tuple(row["upper_exact_mpf_tuple"])),
        ]
    )


def _scalar_upper_interval(ctx: MPIntervalContext, row: dict[str, Any]):
    """Restore an encoded nonnegative scalar as the interval [0, value]."""

    value = mp.make_mpf(tuple(row["exact_mpf_tuple"]))
    return ctx.mpf([mp.mpf(0), value])


def _endpoints(value):
    """Return ordinary MP endpoints without importing the checker helper."""

    lo, hi = value._mpi_
    return mp.make_mpf(tuple(lo)), mp.make_mpf(tuple(hi))


def _symmetric_from_upper(ctx: MPIntervalContext, value):
    """Return [-upper(value), upper(value)] for a nonnegative interval."""

    _, hi = _endpoints(value)
    return ctx.mpf([-hi, hi])


def _q_inverse_squared_coefficients(ctx: MPIntervalContext, center, length: int):
    """Taylor coefficients of (1 + Z^2)^(-2) around ``center``.

    If x = Z - center, q = a + b x + x^2.  The ODE q f' + 2 q' f = 0
    gives a stable directed recurrence for ordinary Taylor coefficients.
    """

    a = 1 + center * center
    b = 2 * center
    coefficients = [1 / (a * a)]
    previous = ctx.mpf(0)
    for n in range(length):
        current = coefficients[n]
        next_coefficient = -(
            b * (n + 2) * current + (n + 3) * previous
        ) / (a * (n + 1))
        coefficients.append(next_coefficient)
        previous = current
    return coefficients


def _independent_q_checks(ctx: MPIntervalContext, center, coefficients, count: int):
    """Compare the first coefficients with independent high-precision mp.diff."""

    checks = []
    with mp.workdps(ctx.dps + 60):
        center_plain = mp.mpf("0.3")

        def q_inverse_squared(z):
            return (1 + z * z) ** (-2)

        tolerance = mp.mpf("1e-220")
        all_contained = True
        for order in range(count):
            reference = mp.diff(q_inverse_squared, center_plain, order) / mp.factorial(order)
            lo, hi = _endpoints(coefficients[order])
            contained = lo - tolerance <= reference <= hi + tolerance
            all_contained = all_contained and contained
            checks.append(
                {
                    "order": order,
                    "reference": reference,
                    "recurrence_interval": coefficients[order],
                    "contained_with_tolerance": contained,
                    "absolute_midpoint_error": (
                        (lo + hi) / 2 - reference
                    ).__abs__(),
                }
            )
    return checks, all_contained


def _old_axis_consistency(base: Path, normalized_center_interval):
    """Compare against the old nominal axis-jet value when available."""

    path = base / OLD_AXIS_NAME
    if not path.exists():
        return {
            "available": False,
            "scope": "No old axis-pressure receipt was present.",
        }
    try:
        receipt = json.loads(path.read_text(encoding="utf-8"))
        candidates = [
            item
            for item in receipt.get("rows", [])
            if abs(float(item.get("Z", 99)) - 0.3) <= 1e-12
        ]
        row = candidates[0] if candidates else None
    except Exception:
        # Keep this diagnostic optional and nonfatal.  The accepted-mass jet is
        # the authoritative output of this module.
        row = None
    if row is None:
        return {
            "available": False,
            "scope": "Old axis-pressure receipt had no Z = 0.3 row.",
        }
    nominal = mp.mpf(row["full_P0"]["arbitrary_exponent_value"])
    lo, hi = _endpoints(normalized_center_interval)
    return {
        "available": True,
        "old_receipt": OLD_AXIS_NAME,
        "old_nominal_normalized_full_P0": nominal,
        "new_center_normalized_interval": normalized_center_interval,
        "nominal_inside_new_interval": lo <= nominal <= hi,
        "scope": "Nominal consistency only; old receipt is not an interval certificate.",
    }


def run(length: int = DEFAULT_LENGTH, precision: int = DEFAULT_PRECISION):
    if int(length) != length or length < 0:
        raise ValueError("length must be a nonnegative integer")
    if int(precision) != precision or precision < 80:
        raise ValueError("precision must be an integer at least 80")
    length = int(length)
    precision = int(precision)

    fixed_path = HERE / FIXED_NAME
    complete_path = HERE / COMPLETE_NAME
    fixed = json.loads(fixed_path.read_text(encoding="utf-8"))
    complete = json.loads(complete_path.read_text(encoding="utf-8"))
    if len(fixed.get("stages", {})) != 13:
        raise AssertionError("fixed-beta receipt must contain 13 stages")
    if not complete.get("all_pressure_stages_included"):
        raise AssertionError("complete pressure receipt does not include all stages")
    accepted_sha = fixed.get("accepted_schedule_sha256")
    if accepted_sha != complete.get("accepted_schedule_sha256"):
        raise AssertionError("accepted schedule hashes differ")
    if not fixed.get("accepted_source_alignment_certified"):
        raise AssertionError("fixed-beta source alignment is not certified")
    if not complete.get("accepted_source_alignment_certified"):
        raise AssertionError("complete pressure source alignment is not certified")

    ctx = MPIntervalContext()
    ctx.dps = precision
    with mp.workdps(precision + 40):
        mass0 = ctx.mpf(0)
        mass2 = ctx.mpf(0)
        stage_counts = {"beta0": 0, "beta2": 0}
        stage_mass_upper = {}
        for name, row in fixed["stages"].items():
            beta = int(row["beta"])
            if beta not in (0, 2):
                raise AssertionError(f"unexpected fixed beta {beta!r} in {name}")
            mass = _interval_from_exact(ctx, row["true_mass_interval"])
            lo, _ = _endpoints(mass)
            if lo < 0:
                raise AssertionError(f"negative accepted true mass in {name}")
            stage_mass_upper[name] = mass
            stage_counts[f"beta{beta}"] += 1
            if beta == 0:
                mass0 += mass
            else:
                mass2 += mass

        flatten_row = complete["flatten_true_mass_upper"]
        flatten_mass = _scalar_upper_interval(ctx, flatten_row)
        _, flatten_upper = _endpoints(flatten_mass)

        center = ctx.mpf(CENTER_Z)
        cauchy_radius = ctx.mpf(CAUCHY_RADIUS)
        q_lower = 1 - (abs(center) + cauchy_radius) ** 2
        if q_lower <= 0:
            raise AssertionError("flatten Cauchy disk reaches the q pole")

        q_coefficients = _q_inverse_squared_coefficients(ctx, center, length)
        q_checks, q_checks_pass = _independent_q_checks(
            ctx, center, q_coefficients, min(length + 1, 8)
        )
        if not q_checks_pass:
            raise AssertionError("independent fixed-beta q coefficient check failed")

        pstar_squared = ctx.exp(28)
        rows = []
        for order in range(length + 1):
            fixed_beta2 = mass2 * q_coefficients[order]
            fixed_beta0 = mass0 if order == 0 else ctx.mpf(0)
            if order == 0:
                flatten = flatten_mass
            else:
                flatten = _symmetric_from_upper(
                    ctx, flatten_mass * (q_lower ** (-2)) * (4 ** order)
                )
            normalized_total = fixed_beta2 + fixed_beta0 + flatten
            physical_total = -pstar_squared * normalized_total
            rows.append(
                {
                    "order": order,
                    "fixed_beta2_normalized": fixed_beta2,
                    "fixed_beta0_normalized": fixed_beta0,
                    "flatten_normalized": flatten,
                    "total_normalized_pressure_coefficient": -normalized_total,
                    "physical_pressure_coefficient": physical_total,
                }
            )

        center_interval = rows[0]["total_normalized_pressure_coefficient"]
        old_consistency = _old_axis_consistency(HERE, center_interval)

        input_paths = [fixed_path, complete_path]
        report = {
            "center_Z": center,
            "Taylor_length": length,
            "precision": precision,
            "Cauchy_radius": cauchy_radius,
            "flatten_q_modulus_lower": q_lower,
            "flatten_coefficient_bound_formula": "M_flat * q_lower^(-2) * 4^k for k > 0; [0, M_flat] at k = 0",
            "fixed_beta2_formula": "M_beta2 * (1 + Z^2)^(-2)",
            "fixed_beta0_formula": "M_beta0",
            "physical_pressure_formula": "-exp(28) * (fixed_beta2 + fixed_beta0 + flatten)",
            "pressure_units": "physical_P",
            "normalized_pressure_units": "P0_over_Pstar_squared",
            "Pstar_squared": pstar_squared,
            "fixed_beta0_mass_upper_normalized": mass0,
            "fixed_beta2_mass_upper_normalized": mass2,
            "flatten_true_mass_upper_normalized": flatten_mass,
            "stage_counts": stage_counts,
            "stage_count_total": len(fixed["stages"]) + 1,
            "stage_mass_intervals": stage_mass_upper,
            "q_inverse_squared_coefficients": q_coefficients,
            "fixed_q_independent_checks": q_checks,
            "fixed_q_independent_checks_pass": q_checks_pass,
            "old_axis_pressure_consistency": old_consistency,
            "Taylor_rows": rows,
            "input_hashes": {
                path.name: hashlib.sha256(path.read_bytes()).hexdigest()
                for path in input_paths
            },
            "accepted_schedule_sha256": accepted_sha,
            "all_14_true_pressure_stages_included": True,
            "finite_pressure_quadrature_approximation_used": False,
            "pressure_parameter_order_truncated": False,
            "original_parameter_errors_enclosed": False,
            "core_source_error_enclosed": False,
            "adapter_evaluation_roundoff_enclosed": False,
            "five_defect_interval_closure": False,
            "accepted_stored_schedule_scope": True,
            "source_error_scope": "accepted stored schedule only; source and parameter uncertainties remain outside this jet",
            "no_shared_core_or_matching_claim": True,
        }
        output = HERE / Path(__file__).with_suffix(".json").name
        output.write_text(
            json.dumps(encode(report), indent=2) + "\n", encoding="utf-8"
        )
        print(
            "pressure axis jets length",
            length,
            "q checks",
            q_checks_pass,
            "old consistency",
            old_consistency.get("nominal_inside_new_interval"),
            flush=True,
        )
        return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--length", type=int, default=DEFAULT_LENGTH)
    parser.add_argument("--precision", type=int, default=DEFAULT_PRECISION)
    args = parser.parse_args()
    run(args.length, args.precision)
