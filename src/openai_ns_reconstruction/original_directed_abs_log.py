"""Absolute-value logarithmic caps without scalar-context endpoint rounding.

This is a bound operator, never a selector of a signed function value.
Only exact MPF tuples are compared or have their sign bit removed. All
logarithmic arithmetic then stays in the caller's interval context.
"""
from __future__ import annotations

from mpmath import make_mpf
from mpmath.ctx_iv import MPIntervalContext
from mpmath.libmp import mpf_cmp


def _finite_tuple(row):
    return (isinstance(row, tuple) and len(row) == 4
            and all(type(value) is int for value in row)
            and row[0] in (0, 1) and row[1] >= 0
            and row[3] == row[1].bit_length())


def abs_upper(context, value):
    """Enclose the exact maximum endpoint magnitude, including exact zero."""
    if not isinstance(context, MPIntervalContext) or not hasattr(value, "_mpi_"):
        raise TypeError("An interval context and a finite interval are required")
    lo, hi = value._mpi_
    if not (_finite_tuple(lo) and _finite_tuple(hi)) or mpf_cmp(lo, hi) > 0:
        raise ValueError("Finite ordered exact interval endpoints are required")
    lo_abs, hi_abs = (0, *lo[1:]), (0, *hi[1:])
    magnitude = lo_abs if mpf_cmp(lo_abs, hi_abs) >= 0 else hi_abs
    exact = make_mpf(magnitude)
    return context.mpf([exact, exact])


def abs_log_upper(context, value):
    """Return a directed logarithmic absolute cap; None means exact zero."""
    if value is None:
        return None
    magnitude = abs_upper(context, value)
    if magnitude._mpi_[1][1] == 0:
        return None
    return context.ln(magnitude)


def fixture_bound_check():
    """New rational and signed-interval witnesses, independent of global mp.prec."""
    from mpmath import mp
    from fractions import Fraction

    previous = mp.prec
    tested = 0
    try:
        for scalar_precision in (37, 53, 113):
            mp.prec = scalar_precision
            context = MPIntervalContext()
            context.prec = 1600
            for fraction in (Fraction(1, 3), Fraction(7, 10), Fraction(-7, 10),
                             Fraction(4, 5), Fraction(10**100, 3), Fraction(1, 10**100)):
                value = context.mpf(fraction.numerator) / fraction.denominator
                upper = abs_upper(context, value)
                log_cap = abs_log_upper(context, value)
                for endpoint in value._mpi_:
                    magnitude = (0, *endpoint[1:])
                    if mpf_cmp(upper._mpi_[1], magnitude) < 0:
                        raise AssertionError("Directed magnitude underbounds an input endpoint")
                reference = context.ln(upper)
                if mpf_cmp(log_cap._mpi_[1], reference._mpi_[1]) < 0:
                    raise AssertionError("Directed logarithm loses its upper endpoint")
                tested += 1
            crossing = context.mpf([-context.mpf(7)/10, context.mpf(1)/3])
            bound = abs_upper(context, crossing)
            if bound.a < (context.mpf(7)/10).a:
                raise AssertionError("Signed interval magnitude is incorrect")
            if abs_log_upper(context, context.mpf(0)) is not None:
                raise AssertionError("Exact zero sentinel changed")
    finally:
        mp.prec = previous
    return dict(passed=True, rational_cases=tested, scalar_precisions=[37, 53, 113],
                interval_precision_bits=1600, signed_crossing_checked=True,
                exact_zero_checked=True, scalar_abs_called=False,
                signed_function_value_selected=False, accepted_self_checks_replayed=False)


__all__ = ["abs_upper", "abs_log_upper", "fixture_bound_check"]
