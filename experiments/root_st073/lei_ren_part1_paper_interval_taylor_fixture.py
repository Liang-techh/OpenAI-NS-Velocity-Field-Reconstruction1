"""Independent checks for the bounded interval Taylor algebra."""

import json
from pathlib import Path

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext

from lei_ren_part1_paper_interval_taylor import (
    IntervalTaylor,
    constant,
    integrate_symmetric,
    variable,
)


def bounds(value):
    return tuple(mp.make_mpf(item) for item in value._mpi_)


def reference_exp_rational(x):
    return mp.exp(x) / (1 + x * x)


def reference_poly_exp(x):
    return x * x * mp.exp(x)


def main():
    with mp.workdps(110):
        ctx = MPIntervalContext()
        ctx.dps = 95
        order = 8
        x0 = mp.mpf("0.3")
        x = variable(ctx, ctx.mpf("0.3"), order)
        rational = x.exp() / (constant(ctx, 1, order) + x * x)
        derivative_checks = []
        for degree, coefficient in enumerate(rational.coefficients):
            exact = mp.diff(reference_exp_rational, x0, degree) / mp.factorial(degree)
            lower, upper = bounds(coefficient)
            assert lower <= exact <= upper, (degree, lower, exact, upper)
            derivative_checks.append({
                "degree": degree,
                "lower": mp.nstr(lower, 35),
                "upper": mp.nstr(upper, 35),
                "reference": mp.nstr(exact, 35),
            })

        radius_text = "0.25"
        radius = ctx.mpf(radius_text)
        integral_checks = []
        cell = ctx.mpf(["-0.25", "0.25"])
        for degree in (5, 6):
            center_x = variable(ctx, ctx.mpf(0), degree)
            cell_x = variable(ctx, cell, degree)
            center_function = center_x * center_x * center_x.exp()
            cell_function = cell_x * cell_x * cell_x.exp()
            enclosure = integrate_symmetric(center_function, cell_function, radius)
            lower, upper = bounds(enclosure)
            exact = mp.quad(reference_poly_exp, [-mp.mpf(radius_text), mp.mpf(radius_text)])
            assert lower <= exact <= upper, (degree, lower, exact, upper)
            integral_checks.append({
                "order": degree,
                "lower": mp.nstr(lower, 35),
                "upper": mp.nstr(upper, 35),
                "reference": mp.nstr(exact, 35),
            })

        # A cell-valued C0 variable produces interval coefficient bounds.
        cell_variable = variable(ctx, cell, 4)
        cell_rational = cell_variable.exp() / (constant(ctx, 1, 4) + cell_variable * cell_variable)
        cell_coefficient_bounds = [
            [mp.nstr(item, 30) for item in bounds(coefficient)]
            for coefficient in cell_rational.coefficients
        ]

        result = {
            "all_checks_passed": True,
            "directed_interval_arithmetic": True,
            "ordinary_taylor_coefficients": True,
            "source_certificate": False,
            "derivative_checks": derivative_checks,
            "symmetric_integral_checks": integral_checks,
            "cell_variable_coefficients": cell_coefficient_bounds,
        }
    output = Path(__file__).with_suffix(".json")
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"all_checks_passed": True, "output": str(output)}))


if __name__ == "__main__":
    main()
