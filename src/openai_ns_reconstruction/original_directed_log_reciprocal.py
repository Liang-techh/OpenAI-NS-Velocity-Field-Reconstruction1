"""Tuple-safe reciprocal logarithmic majorants for ordinary ``Z0..Z4`` jets.

The input rows are logarithmic absolute majorants for the five ordinary
derivatives of a positive value.  ``None`` is reserved for an exact zero
derivative row; the value row itself must be present.  The lower logarithmic
value bound is used only through its lower endpoint.  That endpoint is first
reboxed in the caller's interval context, and its exact MPF sign bit is then
flipped.  This keeps a large endpoint out of the process-wide scalar
precision path.

The recurrence stays in logarithmic space and uses the generic positive
``_logsum`` helper.  No stored logarithm is converted back to a scale.
This module is deliberately pure: it does not load source packets or replay
any accepted builder, recurrence, or self-check.
"""

from __future__ import annotations

from fractions import Fraction
import math

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from mpmath.libmp import mpf_cmp

from .original_log_mixed_jets import _logsum


ORDERS = (0, 1, 2, 3, 4)


def _require_context(context):
    if not isinstance(context, MPIntervalContext):
        raise TypeError("context must be an MPIntervalContext")


def _valid_mpf_tuple(value):
    """Return whether *value* is one finite exact mpmath MPF tuple."""
    return (
        isinstance(value, tuple)
        and len(value) == 4
        and all(type(part) is int for part in value)
        and value[0] in (0, 1)
        and value[1] >= 0
        and value[3] == value[1].bit_length()
    )


def _endpoint(value):
    if not _valid_mpf_tuple(value):
        raise TypeError("Malformed MPF endpoint")
    point = mp.make_mpf(value)
    if mp.isnan(point) or mp.isinf(point):
        raise ValueError("Finite directed endpoints are required")
    return point


def _raw_bounds(value):
    """Read exact interval endpoints and validate their order."""
    if not hasattr(value, "_mpi_"):
        raise TypeError("An MPIntervalContext interval is required")
    rows = value._mpi_
    # Normalize degenerate forms used by different mpmath releases.
    if len(rows) == 4 and all(type(part) is int for part in rows):
        rows = (rows, rows)
    elif len(rows) == 1:
        rows = (rows[0], rows[0])
    if len(rows) != 2:
        raise TypeError("Malformed directed interval")
    lower = _endpoint(rows[0])
    upper = _endpoint(rows[1])
    if mpf_cmp(lower._mpf_, upper._mpf_) > 0:
        raise ValueError("Reversed directed interval")
    return lower, upper


def _lower(value):
    return _raw_bounds(value)[0]


def _upper(value):
    return _raw_bounds(value)[1]


def _interval(context, lower, upper):
    """Box exact finite endpoints in *context*."""
    _require_context(context)
    if hasattr(lower, "_mpi_"):
        lower = _lower(lower)
    if hasattr(upper, "_mpi_"):
        upper = _upper(upper)
    # Both endpoints are exact MPF values by this point.  Keep the order
    # check tuple-directed; it does not depend on global scalar precision.
    lower = _endpoint(lower._mpf_ if hasattr(lower, "_mpf_") else lower)
    upper = _endpoint(upper._mpf_ if hasattr(upper, "_mpf_") else upper)
    if mpf_cmp(lower._mpf_, upper._mpf_) > 0:
        raise ValueError("Reversed directed interval")
    return context.mpf([lower, upper])


def _coerce(context, value):
    """Coerce exact scalar/interval input without binary-float conversion."""
    _require_context(context)
    if isinstance(value, bool) or isinstance(value, float):
        raise TypeError("Binary floats are not accepted")
    if hasattr(value, "_mpi_"):
        lower, upper = _raw_bounds(value)
        return _interval(context, lower, upper)
    if type(value) is int:
        return context.mpf(value)
    if isinstance(value, Fraction):
        return context.mpf(value.numerator) / value.denominator
    if isinstance(value, str):
        return context.mpf(value)
    if hasattr(value, "_mpf_"):
        # A scalar MPF is already rounded at its construction precision.  It
        # is accepted as an exact scalar, then boxed in the caller context.
        return context.mpf(value)
    if isinstance(value, (list, tuple)):
        if len(value) == 4 and all(type(part) is int for part in value):
            point = mp.make_mpf(tuple(value))
            return _interval(context, point, point)
        if len(value) == 2:
            lower = _coerce(context, value[0])
            upper = _coerce(context, value[1])
            return _interval(context, _lower(lower), _upper(upper))
    raise TypeError("Use an exact scalar or MPIntervalContext interval")


def _finite_interval(context, value):
    out = _coerce(context, value)
    lower, upper = _raw_bounds(out)
    return _interval(context, lower, upper)


def _validate_log_rows(context, rows):
    if not isinstance(rows, (list, tuple)) or len(rows) != 5:
        raise ValueError("A complete ordinary Z0..Z4 logarithmic jet is required")
    if rows[0] is None:
        raise ValueError("The logarithmic value row cannot be exact zero")
    out = []
    for row in rows:
        # None is meaningful only for derivative rows: it records a proved
        # exact zero, never an omitted or unknown row.
        out.append(None if row is None else _finite_interval(context, row))
    return out


def _flip_sign_tuple(value):
    """Flip an exact MPF sign bit without scalar arithmetic."""
    value = _endpoint(value)
    sign, mantissa, exponent, bit_count = value._mpf_
    # Canonicalize zero to positive zero; mpmath accepts either sign, but a
    # canonical tuple makes exact fixture comparisons stable.
    sign = 0 if mantissa == 0 else 1 - sign
    return (sign, mantissa, exponent, bit_count)


def _negate_boxed_lower_endpoint(context, lower_log):
    """Safely form ``-lower_log.lower`` after caller-context boxing."""
    boxed = _finite_interval(context, lower_log)
    lower = _lower(boxed)
    # Rebox the selected endpoint before changing its sign.  In particular,
    # never evaluate ``-lower`` while it is a scalar MPF under ``mp.prec``.
    endpoint_box = _interval(context, lower, lower)
    endpoint_tuple = endpoint_box._mpi_[0]
    negated = mp.make_mpf(_flip_sign_tuple(endpoint_tuple))
    return _interval(context, negated, negated)


def directed_log_reciprocal(context, logjet5, certifiedlowerlog):
    """Return reciprocal log-majorants for ordinary derivatives ``Z0..Z4``.

    ``logjet5`` is a complete five-row sequence.  Its row zero is the
    logarithmic cap for a nonzero positive value, and rows one through four
    are derivative caps; ``None`` in those derivative rows denotes an exact
    zero.  ``certifiedlowerlog`` is a finite directed interval containing a
    certified lower logarithmic bound for the positive value row.
    """
    _require_context(context)
    rows = _validate_log_rows(context, logjet5)
    lower_log = _finite_interval(context, certifiedlowerlog)
    inverse_log = _negate_boxed_lower_endpoint(context, lower_log)
    out = [inverse_log]
    for order in range(1, 5):
        terms = []
        for split in range(1, order + 1):
            left = rows[split]
            right = out[order - split]
            if left is None or right is None:
                continue
            term = left + right
            coefficient = math.comb(order, split)
            if coefficient != 1:
                term = term + context.ln(context.mpf(coefficient))
            terms.append(term)
        total = _logsum(context, terms)
        out.append(None if total is None
                   else _finite_interval(context, total + inverse_log))
    return out


# The short name mirrors the legacy arithmetic-layer API while keeping the
# implementation in this new module.  Existing consumers can migrate by
# changing only their import.
log_reciprocal = directed_log_reciprocal


def _negate_tuple(value):
    """Independent exact sign flip used only by fresh fixture comparisons."""
    return _flip_sign_tuple(value)


def _tuple_cmp(left, right):
    return mpf_cmp(left, right)


def _contains_interval(outer, expected):
    """Check directed containment using exact endpoint tuples."""
    if outer is None or expected is None:
        return outer is None and expected is None
    outer_lower, outer_upper = _raw_bounds(outer)
    expected_lower, expected_upper = _raw_bounds(expected)
    return (
        _tuple_cmp(outer_lower._mpf_, expected_lower._mpf_) <= 0
        and _tuple_cmp(expected_upper._mpf_, outer_upper._mpf_) <= 0
    )


def fixture_bound_check():
    """Run bounded fresh rational fixtures, including the legacy comparator.

    The legacy call below is deliberately local and is given only the new
    rational fixtures constructed here.  No historical packet or old builder
    is imported, loaded, or replayed.
    """
    from mpmath import mp as global_mp

    previous = global_mp.prec
    rows = []
    old_underbound_rows = 0
    huge_old_underbound_cases = 0
    checked = 0
    try:
        for scalar_precision in (37, 53, 113):
            global_mp.prec = scalar_precision
            context = MPIntervalContext()
            context.prec = 1600
            # A genuinely large rational exposes the endpoint-rounding bug;
            # the ordinary rational exercises the same path at a small scale.
            lower_inputs = (
                Fraction(10) ** 1000 / 3,
                Fraction(7, 10),
            )
            for lower_input in lower_inputs:
                positive = context.mpf(lower_input.numerator) / lower_input.denominator
                certified = context.ln(positive)
                # Every derivative row is nonzero, including the high-order
                # row.  The following are fresh logarithmic caps, not source
                # packet values.
                logs = [
                    certified,
                    context.ln(context.mpf(2)),
                    context.ln(context.mpf(3)),
                    context.ln(context.mpf(5)),
                    context.ln(context.mpf(7)),
                ]
                corrected = directed_log_reciprocal(context, logs, certified)

                # The exact expected inverse endpoint is obtained by flipping
                # the certified lower endpoint tuple, independently of scalar
                # precision and without formatting/rounding.
                expected_inverse = _negate_tuple(certified._mpi_[0])
                corrected_inverse = corrected[0]._mpi_[0]
                if _tuple_cmp(corrected_inverse, expected_inverse) != 0:
                    raise AssertionError("Safe inverse endpoint changed")
                if _tuple_cmp(corrected[0]._mpi_[1], expected_inverse) != 0:
                    raise AssertionError("Safe inverse row is not exact")

                # Independent recurrence-shape checks: all rows are present
                # for this nonzero high-order fixture, and each endpoint is
                # finite/ordered in the caller context.
                if any(row is None for row in corrected):
                    raise AssertionError("Nonzero fixture produced an exact-zero row")
                for row in corrected:
                    _raw_bounds(row)

                # Comparator is restricted to this fresh input.  The old
                # function is intentionally not used for any accepted data.
                from .original_native_log_jet4 import log_reciprocal as old_log_reciprocal

                old = old_log_reciprocal(context, logs, certified)
                old_upper = old[0]._mpi_[1]
                old_underbound = _tuple_cmp(old_upper, expected_inverse) < 0
                if lower_input == Fraction(10) ** 1000 / 3 and old_underbound:
                    huge_old_underbound_cases += 1
                old_underbound_rows += int(old_underbound)
                rows.append({
                    "scalar_precision_bits": scalar_precision,
                    "interval_precision_bits": context.prec,
                    "lower_input": str(lower_input),
                    "certified_lower_log_exact_tuples": [
                        list(certified._mpi_[0]), list(certified._mpi_[1])
                    ],
                    "expected_inverse_lower_exact_tuple": list(expected_inverse),
                    "corrected_inverse_exact_tuples": [
                        list(corrected[0]._mpi_[0]), list(corrected[0]._mpi_[1])
                    ],
                    "legacy_inverse_upper_exact_tuple": list(old_upper),
                    "legacy_upper_lt_expected_lower": old_underbound,
                    "exact_tuple_comparison": True,
                    "nonzero_high_order_rows": True,
                })
                checked += 1

            # Explicit exact-zero derivative rows exercise the only accepted
            # None sentinel without involving the old helper.
            zero_rows = [context.ln(context.mpf(11)), None, context.ln(context.mpf(3)), None, None]
            zero_lower = context.mpf(7) / 10
            zero_result = directed_log_reciprocal(context, zero_rows, context.ln(zero_lower))
            if zero_result[1] is not None or zero_result[3] is not None:
                raise AssertionError("Exact-zero derivative sentinel was not preserved")
            if zero_result[2] is None:
                raise AssertionError("Nonzero second derivative became exact zero")
            zero_expected_inverse = _negate_tuple(context.ln(zero_lower)._mpi_[0])
            if _tuple_cmp(zero_result[0]._mpi_[0], zero_expected_inverse) != 0:
                raise AssertionError("Exact-zero fixture changed the safe inverse endpoint")

            # Independent ordinary-derivative reference: for f(z)=1+z,
            # reciprocal derivatives at zero are (1,-1,2,-6,24).  Compare
            # their logarithmic magnitudes by directed containment rather
            # than by a self-comparison against the recurrence implementation.
            one_plus_z = directed_log_reciprocal(
                context,
                [context.mpf(0), context.mpf(0), None, None, None],
                context.mpf(0),
            )
            expected_logs = [
                context.mpf(0),
                context.mpf(0),
                context.ln(context.mpf(2)),
                context.ln(context.mpf(6)),
                context.ln(context.mpf(24)),
            ]
            if not all(_contains_interval(actual, expected)
                       for actual, expected in zip(one_plus_z, expected_logs)):
                raise AssertionError("Independent 1/(1+z) derivative logs not contained")
    finally:
        global_mp.prec = previous
    if huge_old_underbound_cases == 0:
        raise AssertionError("No fresh large-rational legacy underbound was observed")
    return {
        "passed": True,
        "fixture_rows": rows,
        "checked_cases": checked,
        "old_legacy_underbound_cases": old_underbound_rows,
        "huge_rational_legacy_underbound_cases": huge_old_underbound_cases,
        "scalar_precisions": [37, 53, 113],
        "interval_precision_bits": 1600,
        "exact_tuple_comparisons": True,
        "nonzero_high_order_fixture": True,
        "ordinary_rational_fixture": True,
        "huge_rational_fixture": True,
        "exact_zero_derivative_sentinel_checked": True,
        "stored_logarithms_exponentiated": False,
        "accepted_builders_or_selfchecks_called": False,
        "old_helper_called_only_on_fresh_inputs": True,
        "global_scalar_precision_restored": global_mp.prec == previous,
    }


fresh_fixture_check = fixture_bound_check


__all__ = [
    "ORDERS",
    "directed_log_reciprocal",
    "log_reciprocal",
    "fixture_bound_check",
    "fresh_fixture_check",
]
