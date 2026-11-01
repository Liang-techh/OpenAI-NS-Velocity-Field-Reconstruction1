"""Bounded directed addition of finite raw mpmath dyadic intervals.

The wire representation is the finite four-integer mpmath ``mpf`` tuple
``(sign, mantissa, exponent, bitcount)``.  Its value is
``(-1)**sign * mantissa * 2**exponent``.  This module deliberately keeps
that representation all the way through the public API:

* endpoints and intervals are tuples, never scalar values or ``Fraction``s;
* exponents are retained as integers and are never used as shift counts;
* every endpoint is first rounded in the requested signed direction to a
  bounded precision, with the discarded part retained as a sticky bit;
* addition aligns only bounded mantissas.  If an exponent gap is too large
  for alignment, the discarded nonzero term selects the adjacent ulp in the
  correct direction instead of being silently dropped;
* NaN/Inf encodings, malformed finite tuples, reversed intervals and resource
  overflows are rejected.

The result is an outward interval conditional on the caller's finite input
enclosures.  No midpoint, point value, source/provider identity or old
arithmetic helper is used here.  A tuple of tuples is JSON-safe under the
standard Python JSON encoder (it is emitted as nested arrays).
"""

from __future__ import annotations

import operator
from typing import Any


SCHEMA = "original_bounded_dyadic_sum_v1"
SOURCE_ID = "synthetic:bounded-dyadic-sum-v1"
OUTPUT_SOURCE_ID = "conditional:directed-finite-dyadic-interval-sum"

# These are implementation limits, not mathematical assumptions.  The
# mantissa limit is independent of the working precision so that a caller may
# supply a full finite mpf and ask for a smaller directed enclosure.
DEFAULT_PRECISION_BITS = 256
MAX_PRECISION_BITS = 4096
MAX_INPUT_MANTISSA_BITS = 65_536
MAX_EXPONENT_BITS = 256
MAX_TERMS = 100_000

# JSON-safe static receipt for callers that need to retain the arithmetic
# policy alongside the tuple result.  It is deliberately separate from the
# API result: the public operations accept and return only raw tuple data.
ROUNDING_RECEIPT = {
    "schema": SCHEMA + "_receipt",
    "operation": "directed finite dyadic interval addition",
    "endpoint_rounding": "signed floor/ceil with sticky remainder",
    "large_exponent_path": "bounded mantissa alignment plus adjacent ulp",
    "accumulation": "sequential outward interval sum",
    "point_value_selected": False,
    "special_values_rejected": True,
    "exponent_shift_by_magnitude": False,
    "max_precision_bits": MAX_PRECISION_BITS,
    "max_input_mantissa_bits": MAX_INPUT_MANTISSA_BITS,
    "max_exponent_bits": MAX_EXPONENT_BITS,
    "max_terms": MAX_TERMS,
}

Endpoint = tuple[int, int, int, int]
Interval = tuple[Endpoint, Endpoint]
_ZERO: Endpoint = (0, 0, 0, 0)


class DyadicSumInputError(ValueError):
    """Raised for malformed, special or reversed dyadic input."""


class DyadicSumResourceError(DyadicSumInputError):
    """Raised before an operation would exceed an explicit resource bound."""


def _index(value: Any, label: str) -> int:
    """Read one integer field without accepting bool or scalar coercion."""

    if isinstance(value, bool):
        raise DyadicSumInputError(f"{label} must be an integer, not bool")
    try:
        result = operator.index(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise DyadicSumInputError(f"{label} must support operator.index") from exc
    if isinstance(result, bool) or type(result) is not int:
        raise DyadicSumInputError(f"{label} did not produce a plain integer")
    return result


def _validate_precision(precision: Any) -> int:
    value = _index(precision, "precision")
    if value <= 0:
        raise DyadicSumInputError("precision must be positive")
    if value > MAX_PRECISION_BITS:
        raise DyadicSumResourceError(
            "precision exceeds MAX_PRECISION_BITS; no unbounded mantissa is allowed"
        )
    return value


def _safe_exponent(value: int, label: str) -> int:
    """Keep exponent fields bounded without using the exponent as a shift."""

    if type(value) is not int:
        raise DyadicSumResourceError(f"{label} is not a plain integer")
    if abs(value).bit_length() > MAX_EXPONENT_BITS:
        raise DyadicSumResourceError(
            f"{label} exceeds the {MAX_EXPONENT_BITS}-bit exponent resource limit"
        )
    return value


def _raw_endpoint(value: Any, label: str) -> tuple[Any, Any, Any, Any]:
    """Require exactly a tuple of four fields, never a list or scalar."""

    if not isinstance(value, tuple) or len(value) != 4:
        raise DyadicSumInputError(
            f"{label} must be a tuple of exactly four integer fields"
        )
    return value[0], value[1], value[2], value[3]


def _canonical_magnitude(mantissa: int, exponent: int, label: str) -> Endpoint:
    """Return a positive canonical endpoint with a bounded local shift."""

    if type(mantissa) is not int or mantissa <= 0:
        raise DyadicSumInputError(f"{label} must have a positive mantissa")
    exponent = _safe_exponent(exponent, f"{label}.exponent")
    trailing = (mantissa & -mantissa).bit_length() - 1
    if trailing:
        mantissa >>= trailing
        exponent = _safe_exponent(
            exponent + trailing, f"{label}.normalized exponent"
        )
    return (0, mantissa, exponent, mantissa.bit_length())


def _normalize_endpoint(value: Any, label: str) -> Endpoint:
    """Validate one finite raw mpf tuple and strip only mantissa-local zeros."""

    raw_sign, raw_mantissa, raw_exponent, raw_bitcount = _raw_endpoint(value, label)
    sign = _index(raw_sign, f"{label}.sign")
    mantissa = _index(raw_mantissa, f"{label}.mantissa")
    exponent = _index(raw_exponent, f"{label}.exponent")
    bitcount = _index(raw_bitcount, f"{label}.bitcount")
    _safe_exponent(exponent, f"{label}.exponent")

    if sign not in (0, 1):
        raise DyadicSumInputError(f"{label}.sign must be 0 or 1")
    if mantissa < 0:
        raise DyadicSumInputError(f"{label}.mantissa must be nonnegative")
    if bitcount < 0:
        # mpmath uses negative bit counts for NaN/Inf sentinels.
        raise DyadicSumInputError(f"{label}.bitcount is a special-value sentinel")

    if mantissa == 0:
        # Only the ordinary finite zero is admitted.  Signed zero is
        # canonicalized so it cannot create a fake negative lower endpoint.
        if exponent != 0 or bitcount != 0:
            raise DyadicSumInputError(
                f"{label} is a non-finite or nonstandard zero-mantissa encoding"
            )
        return _ZERO

    actual_bits = mantissa.bit_length()
    if bitcount == 0 or bitcount != actual_bits:
        raise DyadicSumInputError(
            f"{label}.bitcount must equal mantissa.bit_length()"
        )
    if actual_bits > MAX_INPUT_MANTISSA_BITS:
        raise DyadicSumResourceError(
            f"{label}.mantissa exceeds MAX_INPUT_MANTISSA_BITS"
        )

    trailing = (mantissa & -mantissa).bit_length() - 1
    if trailing:
        mantissa >>= trailing
        exponent = _safe_exponent(
            exponent + trailing, f"{label}.normalized exponent"
        )
        bitcount -= trailing
    # The tuple is a signed finite dyadic.  Keep the sign only after zero
    # canonicalization; all arithmetic below assumes a positive magnitude.
    return (sign, mantissa, exponent, bitcount)


def _is_zero(value: Endpoint) -> bool:
    return value[1] == 0


def _compare_abs(left: Endpoint, right: Endpoint) -> int:
    """Compare magnitudes with shifts bounded by the mantissa sizes."""

    left_zero = _is_zero(left)
    right_zero = _is_zero(right)
    if left_zero or right_zero:
        if left_zero and right_zero:
            return 0
        return -1 if left_zero else 1

    left_top = left[2] + left[3]
    right_top = right[2] + right[3]
    if left_top != right_top:
        return 1 if left_top > right_top else -1

    # Equal top positions imply that this shift is a mantissa-bitcount
    # difference, not an astronomical exponent gap.
    shift = left[2] - right[2]
    if abs(shift) > max(left[3], right[3]):
        raise DyadicSumResourceError("internal comparison shift escaped its bound")
    if shift >= 0:
        lm = left[1] << shift
        rm = right[1]
    else:
        lm = left[1]
        rm = right[1] << (-shift)
    return (lm > rm) - (lm < rm)


def _compare_signed(left: Endpoint, right: Endpoint) -> int:
    """Compare signed finite canonical endpoints."""

    left_zero = _is_zero(left)
    right_zero = _is_zero(right)
    if left_zero or right_zero:
        if left_zero and right_zero:
            return 0
        if left_zero:
            return 1 if right[0] else -1
        return -1 if left[0] else 1
    if left[0] != right[0]:
        return -1 if left[0] else 1
    magnitude = _compare_abs(left, right)
    return magnitude if left[0] == 0 else -magnitude


def _signed_magnitude(sign: int, magnitude: Endpoint) -> Endpoint:
    if _is_zero(magnitude):
        return _ZERO
    return (sign, magnitude[1], magnitude[2], magnitude[3])


def _round_magnitude(
    mantissa: int,
    exponent: int,
    precision: int,
    upward: bool,
    *,
    label: str,
) -> Endpoint:
    """Round a positive dyadic magnitude with a finite sticky remainder."""

    if mantissa <= 0 or type(mantissa) is not int:
        raise DyadicSumInputError(f"{label} must have a positive integer mantissa")
    exponent = _safe_exponent(exponent, f"{label}.exponent")
    bits = mantissa.bit_length()
    if bits <= precision:
        return _canonical_magnitude(mantissa, exponent, label)

    drop = bits - precision
    # ``drop`` is a mantissa-local count.  It is bounded by the input or
    # alignment resource cap; it is never an exponent difference.
    quotient = mantissa >> drop
    remainder_nonzero = (mantissa & ((1 << drop) - 1)) != 0
    if upward and remainder_nonzero:
        quotient += 1
    rounded_exponent = _safe_exponent(
        exponent + drop, f"{label}.rounded exponent"
    )
    return _canonical_magnitude(quotient, rounded_exponent, label)


def _round_signed(value: Endpoint, precision: int, direction: int) -> Endpoint:
    """Round one signed value toward -infinity or +infinity."""

    if direction not in (-1, 1):
        raise DyadicSumInputError("internal direction must be -1 or 1")
    if _is_zero(value):
        return _ZERO

    # For a negative value, floor(-x) = -ceil(x), while ceil(-x) = -floor(x).
    upward_magnitude = (direction == 1 and value[0] == 0) or (
        direction == -1 and value[0] == 1
    )
    magnitude = _round_magnitude(
        value[1],
        value[2],
        precision,
        upward_magnitude,
        label="directed endpoint",
    )
    return _signed_magnitude(value[0], magnitude)


def _neighbor_magnitude(
    dominant: Endpoint,
    precision: int,
    delta: int,
) -> Endpoint:
    """Return dominant magnitude plus or minus one p-bit ulp.

    This helper is used only when a nonzero term lies below the p-bit grid of
    the dominant term.  The mantissa shift is at most ``precision`` bits.
    """

    if delta not in (-1, 1):
        raise DyadicSumInputError("internal neighbor delta must be -1 or 1")
    if _is_zero(dominant):
        raise DyadicSumInputError("cannot form an ulp neighbor of zero")
    # Immediately below an exact power of two the p-bit spacing is half the
    # spacing above it.  Use the preceding binade, including at p=1.
    if delta == -1 and dominant[1] == 1:
        predecessor = (1 << precision) - 1
        exponent = _safe_exponent(dominant[2] - precision, "predecessor exponent")
        return _canonical_magnitude(predecessor, exponent, "predecessor magnitude")
    top = dominant[2] + dominant[3]
    ulp_exponent = _safe_exponent(top - precision, "neighbor ulp exponent")
    local_shift = dominant[2] - ulp_exponent
    if local_shift < 0 or local_shift > precision:
        raise DyadicSumResourceError("internal neighbor shift escaped its bound")
    scaled = dominant[1] << local_shift
    neighbor = scaled + delta
    if neighbor <= 0:
        raise DyadicSumResourceError("neighbor would leave the positive magnitude range")
    return _canonical_magnitude(neighbor, ulp_exponent, "neighbor magnitude")


def _add_directed(left: Endpoint, right: Endpoint, precision: int, direction: int) -> Endpoint:
    """Round the exact sum of two already p-bit endpoints directionally."""

    if _is_zero(left):
        return _round_signed(right, precision, direction)
    if _is_zero(right):
        return _round_signed(left, precision, direction)

    # A bounded alignment is enough whenever cancellation or a carry could
    # affect the retained p bits.  The exponent gap is compared as an integer
    # only; it is never passed to a shift operation.
    exponent_gap = abs(left[2] - right[2])
    alignment_limit = 2 * precision + 2
    if exponent_gap <= alignment_limit:
        common_exponent = min(left[2], right[2])
        left_mantissa = left[1] << (left[2] - common_exponent)
        right_mantissa = right[1] << (right[2] - common_exponent)
        if left[0] == right[0]:
            magnitude = left_mantissa + right_mantissa
            upward_magnitude = (
                direction == 1 if left[0] == 0 else direction == -1
            )
            return _signed_magnitude(
                left[0],
                _round_magnitude(
                    magnitude,
                    common_exponent,
                    precision,
                    upward_magnitude,
                    label="aligned sum",
                ),
            )

        if left_mantissa == right_mantissa:
            return _ZERO
        if left_mantissa > right_mantissa:
            sign = left[0]
            magnitude = left_mantissa - right_mantissa
        else:
            sign = right[0]
            magnitude = right_mantissa - left_mantissa
        upward_magnitude = direction == 1 if sign == 0 else direction == -1
        return _signed_magnitude(
            sign,
            _round_magnitude(
                magnitude,
                common_exponent,
                precision,
                upward_magnitude,
                label="aligned difference",
            ),
        )

    # The exponent gap is now more than 2p+2.  Since both operands have at
    # most p bits, the smaller operand is strictly below one p-bit ulp of the
    # larger one.  Its nonzero sticky bit still changes exactly one directed
    # neighbor: same signs use [dominant, dominant + ulp], opposite signs use
    # [dominant - ulp, dominant].  This is the critical cancellation-safe
    # path for a negative microscopic term next to a huge positive term.
    if _compare_abs(left, right) >= 0:
        dominant, tiny = left, right
    else:
        dominant, tiny = right, left
    if _is_zero(tiny):
        return _round_signed(dominant, precision, direction)

    top_gap = abs((dominant[2] + dominant[3]) - (tiny[2] + tiny[3]))
    if top_gap <= precision:
        # This should be unreachable after the exponent-gap test.  Keep the
        # guard explicit so a future representation change cannot drop a term.
        raise DyadicSumResourceError("large-gap addition lost its sticky separation")

    if dominant[0] == tiny[0]:
        lower_magnitude = _canonical_magnitude(
            dominant[1], dominant[2], "dominant magnitude"
        )
        upper_magnitude = _neighbor_magnitude(dominant, precision, 1)
    else:
        lower_magnitude = _neighbor_magnitude(dominant, precision, -1)
        upper_magnitude = _canonical_magnitude(
            dominant[1], dominant[2], "dominant magnitude"
        )

    if (direction == -1 and dominant[0] == 0) or (
        direction == 1 and dominant[0] == 1
    ):
        chosen = lower_magnitude
    else:
        chosen = upper_magnitude
    return _signed_magnitude(dominant[0], chosen)


def _validate_interval(value: Any, label: str) -> Interval:
    if not isinstance(value, tuple) or len(value) != 2:
        raise DyadicSumInputError(
            f"{label} must be a tuple of exactly two endpoint tuples"
        )
    lower = _normalize_endpoint(value[0], f"{label}.lower")
    upper = _normalize_endpoint(value[1], f"{label}.upper")
    if _compare_signed(lower, upper) > 0:
        raise DyadicSumInputError(f"{label} endpoints are reversed")
    return (lower, upper)


def _rounded_interval(value: Interval, precision: int) -> Interval:
    lower = _round_signed(value[0], precision, -1)
    upper = _round_signed(value[1], precision, 1)
    if _compare_signed(lower, upper) > 0:
        raise DyadicSumResourceError("directed endpoint rounding reversed an interval")
    return (lower, upper)


def _sum_validated(left: Interval, right: Interval, precision: int) -> Interval:
    left_r = _rounded_interval(left, precision)
    right_r = _rounded_interval(right, precision)
    lower = _add_directed(left_r[0], right_r[0], precision, -1)
    upper = _add_directed(left_r[1], right_r[1], precision, 1)
    if _compare_signed(lower, upper) > 0:
        raise DyadicSumResourceError("outward addition produced a reversed interval")
    return (lower, upper)


def add_interval(
    left: Any,
    right: Any,
    precision: int = DEFAULT_PRECISION_BITS,
) -> Interval:
    """Return the outward finite-dyadic sum of two closed intervals.

    ``precision`` is a number of significant binary bits, not a decimal
    precision.  Each input endpoint is directed-rounded before the endpoint
    sums are formed, so every discarded bit contributes through a sticky
    floor/ceil decision.  The returned nested tuple contains only integer
    fields and can be passed directly to ``json.dumps``.
    """

    working_precision = _validate_precision(precision)
    left_i = _validate_interval(left, "left")
    right_i = _validate_interval(right, "right")
    return _sum_validated(left_i, right_i, working_precision)


def sum_intervals(
    intervals: Any,
    precision: int = DEFAULT_PRECISION_BITS,
) -> Interval:
    """Return an outward sum of a tuple of closed finite-dyadic intervals.

    The outer collection is intentionally required to be a tuple.  Empty
    input denotes the exact zero interval.  Terms are accumulated one at a
    time with directed rounding; no term is discarded merely because its
    exponent is far from the running sum.  ``MAX_TERMS`` is checked before
    any unbounded iteration can occur.
    """

    working_precision = _validate_precision(precision)
    if not isinstance(intervals, tuple):
        raise DyadicSumInputError("intervals must be a tuple of intervals")
    if len(intervals) > MAX_TERMS:
        raise DyadicSumResourceError("interval count exceeds MAX_TERMS")
    if not intervals:
        return (_ZERO, _ZERO)

    accumulated = _validate_interval(intervals[0], "intervals[0]")
    for index in range(1, len(intervals)):
        current = _validate_interval(intervals[index], f"intervals[{index}]")
        accumulated = _sum_validated(accumulated, current, working_precision)
    return _rounded_interval(accumulated, working_precision)


__all__ = [
    "SCHEMA",
    "SOURCE_ID",
    "OUTPUT_SOURCE_ID",
    "DEFAULT_PRECISION_BITS",
    "MAX_PRECISION_BITS",
    "MAX_INPUT_MANTISSA_BITS",
    "MAX_EXPONENT_BITS",
    "MAX_TERMS",
    "ROUNDING_RECEIPT",
    "Endpoint",
    "Interval",
    "DyadicSumInputError",
    "DyadicSumResourceError",
    "add_interval",
    "sum_intervals",
]
