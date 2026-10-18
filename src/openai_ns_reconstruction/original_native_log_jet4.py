"""Directed ordinary ``Z0..Z4`` arithmetic for native source providers.

The five entries in a value jet are ordinary derivatives, in order.  Value
jets contain directed intervals and use ``c.mpf(0)`` for an exact zero.  A
logarithmic jet contains directed intervals for logarithmic absolute caps and
uses ``None`` only for a proved exact-zero row.  The module deliberately stays
at the arithmetic layer: it has no source loading, numerical selection, or
field/pressure construction.

All logarithmic combinations use the exponent-free ``_logsum`` helper.  In
particular, this module never exponentiates a stored logarithm, so it does not
materialize a possibly enormous scale merely to form a cap.
"""

from __future__ import annotations

from fractions import Fraction
import math

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext

from .original_log_mixed_jets import _logsum


ORDERS = (0, 1, 2, 3, 4)
LOG_DERIVATIVE_ORDERS = (0, 1, 2, 3)


def _require_context(c):
    if not isinstance(c, MPIntervalContext):
        raise TypeError("c must be an MPIntervalContext")


def _mpf_endpoint(value):
    """Return one exact mpmath endpoint from an interval-record component."""
    try:
        return mp.make_mpf(tuple(value))
    except (TypeError, ValueError, IndexError) as exc:
        raise TypeError("Malformed MPIntervalContext interval") from exc


def _raw_bounds(value):
    """Read interval endpoints without using midpoint or endpoint strings."""
    if not hasattr(value, "_mpi_"):
        raise TypeError("An MPIntervalContext interval is required")
    rows = value._mpi_
    # A few mpmath versions expose a degenerate interval as one raw endpoint
    # (or expose the endpoint tuple itself).  Normalize those forms while
    # retaining the exact binary endpoints.
    if len(rows) == 4 and all(type(part) is int for part in rows):
        rows = (rows, rows)
    elif len(rows) == 1:
        rows = (rows[0], rows[0])
    if len(rows) != 2:
        raise TypeError("Malformed MPIntervalContext interval")
    lower = _mpf_endpoint(rows[0])
    upper = _mpf_endpoint(rows[1])
    if any(mp.isnan(endpoint) or mp.isinf(endpoint)
           for endpoint in (lower, upper)):
        raise ValueError("Finite directed intervals are required")
    if lower > upper:
        raise ValueError("Reversed directed interval")
    return lower, upper


def _lower(value):
    return _raw_bounds(value)[0]


def _upper(value):
    return _raw_bounds(value)[1]


def _interval(c, lower, upper):
    """Construct a finite directed interval from exact mpmath endpoints."""
    _require_context(c)
    if lower > upper:
        raise ValueError("Reversed directed interval")
    if any(mp.isnan(endpoint) or mp.isinf(endpoint)
           for endpoint in (lower, upper)):
        raise ValueError("Finite directed intervals are required")
    return c.mpf([lower, upper])


def _scalar(c, value):
    """Coerce an exact scalar or directed interval, rejecting binary floats."""
    _require_context(c)
    if isinstance(value, bool) or isinstance(value, float):
        raise TypeError("Binary floats are not accepted")
    if hasattr(value, "_mpi_"):
        lower, upper = _raw_bounds(value)
        return _interval(c, lower, upper)
    if type(value) is int:
        return c.mpf(value)
    if isinstance(value, Fraction):
        return c.mpf(value.numerator) / value.denominator
    if isinstance(value, str):
        return c.mpf(value)
    if hasattr(value, "_mpf_"):
        # An mpmath scalar has an exact binary representation already.  It is
        # converted through the interval context before entering arithmetic.
        return c.mpf(value)
    if isinstance(value, (list, tuple)) and len(value) == 2:
        left = _scalar(c, value[0])
        right = _scalar(c, value[1])
        return _interval(c, _lower(left), _upper(right))
    raise TypeError("Use an exact scalar or MPIntervalContext interval")


def _finite_interval(c, value):
    """Coerce and normalize one scalar/interval to the caller's context."""
    out = _scalar(c, value)
    # _scalar already checks the raw endpoints.  Reading them once more makes
    # the returned context independent of an interval from another context.
    return _interval(c, _lower(out), _upper(out))


def _value_rows(c, rows):
    """Validate exactly five ordinary signed-value derivative rows."""
    _require_context(c)
    if not isinstance(rows, (list, tuple)) or len(rows) != 5:
        raise ValueError("A value jet must provide exactly five rows")
    out = []
    for order, row in enumerate(rows):
        if row is None:
            raise TypeError(
                "Value jet row %d must use c.mpf(0), not None" % order
            )
        out.append(_finite_interval(c, row))
    return out


def _log_rows(c, rows):
    """Validate exactly five logarithmic absolute-cap rows."""
    _require_context(c)
    if not isinstance(rows, (list, tuple)) or len(rows) != 5:
        raise ValueError("A logarithmic jet must provide exactly five rows")
    out = []
    for row in rows:
        out.append(None if row is None else _finite_interval(c, row))
    return out


def _log_rows4(c, rows):
    """Validate the four nonconstant logarithmic derivative rows."""
    _require_context(c)
    if not isinstance(rows, (list, tuple)) or len(rows) != 4:
        raise ValueError(
            "The logarithmic derivative input must provide exactly four rows"
        )
    return [None if row is None else _finite_interval(c, row) for row in rows]


def _sum_values(c, terms):
    total = c.mpf(0)
    for term in terms:
        total = total + term
    return total


def _log_integer(c, value):
    if type(value) is not int or value <= 0:
        raise ValueError("A positive integer logarithm is required")
    return c.ln(c.mpf(value))


def _weighted_log_term(c, left, right, coefficient):
    term = left + right
    if coefficient != 1:
        term = term + _log_integer(c, coefficient)
    return term


def abs_log(c, interval):
    """Return a logarithmic cap for the absolute value of ``interval``.

    ``None`` is returned only when the supplied interval is the exact zero.
    """
    _require_context(c)
    if interval is None:
        return None
    value = _finite_interval(c, interval)
    lower, upper = _raw_bounds(value)
    magnitude = max(abs(lower), abs(upper))
    if magnitude == 0:
        return None
    return c.ln(_interval(c, magnitude, magnitude))


def log_add(c, *terms):
    """Add logarithmic absolute caps without exponentiating them."""
    _require_context(c)
    normalized = [
        _finite_interval(c, term) for term in terms if term is not None
    ]
    return _logsum(c, normalized)


def value_product(c, left5, right5):
    """Return the ordinary derivative jet of a product through order four."""
    left = _value_rows(c, left5)
    right = _value_rows(c, right5)
    out = []
    for order in ORDERS:
        terms = (
            c.mpf(math.comb(order, split))
            * left[split]
            * right[order - split]
            for split in range(order + 1)
        )
        out.append(_sum_values(c, terms))
    return out


def value_reciprocal(c, positive5):
    """Return the ordinary derivative jet of ``1 / f``.

    The value row must have a strictly positive lower endpoint.  The
    recurrence is the differentiated identity ``f * (1/f) = 1`` and hence
    preserves signed interval information from the supplied value jet.
    """
    value = _value_rows(c, positive5)
    if _lower(value[0]) <= 0:
        raise ValueError("The reciprocal value row must have a positive lower endpoint")
    out = [c.mpf(1) / value[0]]
    for order in range(1, 5):
        total = c.mpf(0)
        for split in range(1, order + 1):
            total = total + (
                c.mpf(math.comb(order, split))
                * value[split]
                * out[order - split]
            )
        out.append(-out[0] * total)
    return out


def value_log_derivatives(c, positive5):
    """Return ``[log(f0), (log f)'(0), ..., (log f)''''(0)]``.

    The derivative rows follow the ordinary recurrence obtained from
    ``f * (log f)' = f'``.  No logarithm of a stored magnitude cap is used;
    the input is a positive signed-value jet.
    """
    value = _value_rows(c, positive5)
    if _lower(value[0]) <= 0:
        raise ValueError("The logarithm value row must have a positive lower endpoint")
    out = [c.ln(value[0])]
    for order in range(1, 5):
        remainder = value[order]
        # Differentiate f*(log f)' exactly order-1 times.  The j=0 term is
        # f0*(log f)^(order), which is the row being solved for.
        for split in range(1, order):
            remainder = remainder - (
                c.mpf(math.comb(order - 1, split))
                * value[split]
                * out[order - split]
            )
        out.append(remainder / value[0])
    return out


def log_product(c, left5, right5):
    """Return logarithmic absolute derivative caps for an ordinary product."""
    left = _log_rows(c, left5)
    right = _log_rows(c, right5)
    out = []
    for order in ORDERS:
        terms = [
            _weighted_log_term(
                c,
                left[split],
                right[order - split],
                math.comb(order, split),
            )
            for split in range(order + 1)
            if left[split] is not None and right[order - split] is not None
        ]
        out.append(_logsum(c, terms))
    return out


def log_reciprocal(c, logjet5, log_positive_lower):
    """Return logarithmic absolute caps for the ordinary reciprocal jet.

    ``log_positive_lower`` is a directed interval containing a certified lower
    logarithmic bound for the positive value row.  Its lower endpoint is used
    for the reciprocal cap, so the returned row is safe even when that input
    bound itself has directed width.
    """
    logjet = _log_rows(c, logjet5)
    if logjet[0] is None:
        raise ValueError("The reciprocal logarithmic value row cannot be exact zero")
    lower_log = _finite_interval(c, log_positive_lower)
    inverse_log = _interval(c, -_lower(lower_log), -_lower(lower_log))
    out = [inverse_log]
    for order in range(1, 5):
        terms = [
            _weighted_log_term(
                c,
                logjet[split],
                out[order - split],
                math.comb(order, split),
            )
            for split in range(1, order + 1)
            if logjet[split] is not None and out[order - split] is not None
        ]
        total = _logsum(c, terms)
        out.append(None if total is None else total + inverse_log)
    return out


def exp_relative(c, ell_derivative_logs4):
    """Return log caps for ``|D^k exp(ell)| / exp(ell)`` through order four.

    The input rows are logarithmic absolute caps for ``ell'`` through
    ``ell''''``.  Bell's ordinary derivative recurrence is evaluated entirely
    in log space, with ``None`` propagating only through exact-zero terms.
    """
    ell = _log_rows4(c, ell_derivative_logs4)
    out = [c.mpf(0)]
    for order in range(1, 5):
        terms = [
            _weighted_log_term(
                c,
                ell[derivative_index],
                out[order - 1 - derivative_index],
                math.comb(order - 1, derivative_index),
            )
            for derivative_index in range(order)
            if ell[derivative_index] is not None
            and out[order - 1 - derivative_index] is not None
        ]
        out.append(_logsum(c, terms))
    return out


def log_scale(c, jet5, logscalar):
    """Add one logarithmic scalar cap to every nonzero row of ``jet5``."""
    jet = _log_rows(c, jet5)
    if logscalar is None:
        return [None for _ in ORDERS]
    scalar = _finite_interval(c, logscalar)
    return [None if row is None else row + scalar for row in jet]


def _contains(c, outer, expected):
    """Check interval containment without selecting a midpoint."""
    if expected is None:
        return outer is None
    if outer is None:
        return False
    expected_interval = _finite_interval(c, expected)
    return (
        _lower(outer) <= _lower(expected_interval)
        and _upper(expected_interval) <= _upper(outer)
    )


def _all_contains(c, outer_rows, expected_rows):
    if len(outer_rows) != len(expected_rows):
        return False
    return all(_contains(c, outer, expected)
               for outer, expected in zip(outer_rows, expected_rows))


def self_check():
    """Run pure rational-point recurrence checks and return their receipt."""
    c = MPIntervalContext()
    c.prec = 180

    # All rows describe ordinary derivatives at the rational point z=0.
    f = [c.mpf(1), c.mpf(1), c.mpf(0), c.mpf(0), c.mpf(0)]
    inv = value_reciprocal(c, f)
    logf = value_log_derivatives(c, f)
    product = value_product(c, f, inv)
    f_logs = [c.mpf(0), c.mpf(0), None, None, None]
    log_product_rows = log_product(c, f_logs, f_logs)
    log_inverse_rows = log_reciprocal(c, f_logs, c.mpf(0))
    exp_z2_rows = exp_relative(c, [None, c.ln(c.mpf(2)), None, None])

    checks = {
        "value_reciprocal_one_plus_z": _all_contains(
            c, inv, [1, -1, 2, -6, 24]
        ),
        "value_log_derivatives_one_plus_z": _all_contains(
            c, logf, [0, 1, -1, 2, -6]
        ),
        "leibniz_product_with_reciprocal": _all_contains(
            c, product, [1, 0, 0, 0, 0]
        ),
        "log_product_one_plus_z_squared": _all_contains(
            c,
            log_product_rows,
            [c.mpf(0), c.ln(c.mpf(2)), c.ln(c.mpf(2)), None, None],
        ),
        "log_reciprocal_one_plus_z": _all_contains(
            c,
            log_inverse_rows,
            [c.mpf(0), c.mpf(0), c.ln(c.mpf(2)),
             c.ln(c.mpf(6)), c.ln(c.mpf(24))],
        ),
        "exp_z_squared_bell_rows": _all_contains(
            c,
            exp_z2_rows,
            [c.mpf(0), None, c.ln(c.mpf(2)), None, c.ln(c.mpf(12))],
        ),
        "absolute_log_interval": _contains(
            c, abs_log(c, c.mpf([-3, 4])), c.ln(c.mpf(4))
        ),
        "log_add_without_exponentiation": _contains(
            c, log_add(c, c.ln(c.mpf(2)), c.ln(c.mpf(3))), c.ln(c.mpf(5))
        ),
    }
    result = {
        "passed": all(checks.values()),
        "checks": checks,
        "ordinary_orders": list(ORDERS),
        "point": "z=0",
        "directed_interval_context": True,
        "source_values_selected": False,
        "stored_logarithms_exponentiated": False,
        "signed_values_returned_from_log_caps": False,
    }
    if not result["passed"]:
        raise ArithmeticError("original_native_log_jet4 self-check failed")
    return result


__all__ = [
    "ORDERS",
    "LOG_DERIVATIVE_ORDERS",
    "abs_log",
    "log_add",
    "value_product",
    "value_reciprocal",
    "value_log_derivatives",
    "log_product",
    "log_reciprocal",
    "exp_relative",
    "log_scale",
    "self_check",
]
