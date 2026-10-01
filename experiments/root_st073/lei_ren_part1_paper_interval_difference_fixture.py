"""Independent checks for nominal/perturbation interval arithmetic."""

import json
from pathlib import Path

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext

from lei_ren_part1_paper_interval_difference import IntervalDifference


def bounds(value):
    return tuple(mp.make_mpf(item) for item in value._mpi_)


def main():
    with mp.workdps(100):
        ctx = MPIntervalContext()
        ctx.dps = 80
        x = IntervalDifference(
            ctx,
            ctx.mpf(["1.8", "2.2"]),
            ctx.mpf(["-0.04", "0.04"]),
        )
        y = IntervalDifference(
            ctx,
            ctx.mpf(["0.9", "1.1"]),
            ctx.mpf(["-0.03", "0.03"]),
        )
        expression = ((x * y + 1) / (y + 2)) ** 3
        reciprocal = x / y
        samples = []
        for x0 in (mp.mpf("1.85"), mp.mpf("2.0"), mp.mpf("2.15")):
            for dx in (mp.mpf("-0.03"), mp.mpf("0.02")):
                for y0 in (mp.mpf(".92"), mp.mpf("1.0"), mp.mpf("1.08")):
                    for dy in (mp.mpf("-0.02"), mp.mpf(".025")):
                        nominal = ((x0 * y0 + 1) / (y0 + 2)) ** 3
                        actual = (((x0 + dx) * (y0 + dy) + 1) / (y0 + dy + 2)) ** 3
                        target_difference = actual - nominal
                        lower, upper = bounds(expression.nominal)
                        d_lower, d_upper = bounds(expression.difference)
                        assert lower <= nominal <= upper
                        assert d_lower <= target_difference <= d_upper
                        r_nominal = x0 / y0
                        r_actual = (x0 + dx) / (y0 + dy)
                        r_lower, r_upper = bounds(reciprocal.nominal)
                        rd_lower, rd_upper = bounds(reciprocal.difference)
                        assert r_lower <= r_nominal <= r_upper
                        assert rd_lower <= r_actual - r_nominal <= rd_upper
                        samples.append({
                            "x0": mp.nstr(x0, 8),
                            "dx": mp.nstr(dx, 8),
                            "y0": mp.nstr(y0, 8),
                            "dy": mp.nstr(dy, 8),
                        })

        zero = IntervalDifference(ctx, ctx.mpf(["1.1", "1.2"]))
        zero_expression = ((zero * zero + 1) / (zero + 2)) ** -3
        zlo, zhi = bounds(zero_expression.difference)
        assert zlo == 0 and zhi == 0

        converter = IntervalDifference.converter(ctx)
        converted = converter(ctx.mpf("1.25"))
        clo, chi = bounds(converted.difference)
        assert clo == 0 and chi == 0

        denominator_rejection = False
        try:
            x / IntervalDifference(ctx, ctx.mpf(["-.1", ".1"]))
        except ZeroDivisionError:
            denominator_rejection = True
        assert denominator_rejection

        # The one-argument converter is compatible with the existing radial
        # recurrence and preserves an exactly zero perturbation component.
        from lei_ren_part1_paper_core_recursion import core_coefficients

        zero_atom = converter(0)
        count = 5
        formal = core_coefficients(
            ".2",
            ".01",
            F0_Z_taylor=[converter("1.2"), converter(".1")] + [zero_atom] * (count - 2),
            U0_Z_taylor=[converter(".4"), converter(".05")] + [zero_atom] * (count - 2),
            P0_Z_taylor=[converter("-1"), converter(".02")] + [zero_atom] * (count - 2),
            radial_degree=2,
            precision=80,
            scalar_converter=converter,
        )
        formal_zero_components = 0
        for rows in (formal["F"], formal["Uz"], formal["P"]):
            for row in rows:
                for coefficient in row:
                    if isinstance(coefficient, IntervalDifference):
                        lo, hi = bounds(coefficient.difference)
                        assert lo == 0 and hi == 0
                        formal_zero_components += 1

        result = {
            "all_checks_passed": True,
            "directed_interval_arithmetic": True,
            "independent_resolved_samples": len(samples),
            "zero_difference_expression_exact": True,
            "core_converter_zero_components": formal_zero_components,
            "denominator_crossing_rejected": denominator_rejection,
            "source_certificate": False,
        }
    output = Path(__file__).with_suffix(".json")
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"all_checks_passed": True, "output": str(output)}))


if __name__ == "__main__":
    main()
