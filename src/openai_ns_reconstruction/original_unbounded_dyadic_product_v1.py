"""Conditional exact products of finite raw mpmath dyadic intervals.

The input wire format is the finite mpmath ``mpf`` endpoint tuple
``[sign, mantissa, exponent, bitcount]``.  Its value is
``(-1)**sign * mantissa * 2**exponent``.  This leaf deliberately stays at
that representation boundary:

* exponents are carried as Python integers and are never used as shift sizes;
* comparison uses the signed top bit ``exponent + bitcount`` and only shifts
  a mantissa after equality of the top bits has bounded the shift;
* multiplication multiplies mantissas and adds exponents, then canonicalizes
  trailing zero bits;
* interval multiplication selects the minimum and maximum of its four exact
  endpoint products.

The result is a finite raw interval product conditional on the caller's
finite endpoint enclosures.  No point is selected, and no source/provider
identity is inferred or certified here.
"""

from __future__ import annotations

import operator
from typing import Any, Callable


SCHEMA = "original_unbounded_dyadic_product_v1"
SOURCE_ID = "synthetic:unbounded-dyadic-product-v1"
LEFT_SOURCE_ID = "caller:left-finite-dyadic-enclosure"
RIGHT_SOURCE_ID = "caller:right-finite-dyadic-enclosure"
OUTPUT_SOURCE_ID = "conditional:finite-dyadic-interval-product"

# Exponents do not drive memory allocation, but their own integer size is
# bounded for resource safety; this still admits magnitudes near 2**256.
MAX_MANTISSA_BITS = 4096
MAX_PRODUCT_MANTISSA_BITS = 2 * MAX_MANTISSA_BITS
ABSOLUTE_MAX_MANTISSA_BITS = 65_536
MAX_EXPONENT_BITS = 256

Endpoint = list[int]
Interval = list[Endpoint]


class DyadicProductInputError(ValueError):
    """Raised when a raw endpoint or interval is not a finite dyadic."""


def _index(value: Any, label: str) -> int:
    """Read an exact integer field without accepting bool or scalar coercion."""

    if isinstance(value, bool):
        raise DyadicProductInputError(f"{label} must be an integer, not bool")
    try:
        result = operator.index(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise DyadicProductInputError(f"{label} must support operator.index") from exc
    if isinstance(result, bool) or type(result) is not int:
        raise DyadicProductInputError(f"{label} did not produce a plain integer")
    return result


def _validate_limit(max_mantissa_bits: Any) -> int:
    limit = _index(max_mantissa_bits, "max_mantissa_bits")
    if limit <= 0:
        raise DyadicProductInputError("max_mantissa_bits must be positive")
    if limit > ABSOLUTE_MAX_MANTISSA_BITS:
        raise DyadicProductInputError(
            "max_mantissa_bits exceeds the finite implementation resource cap"
        )
    return limit


def _raw_endpoint(value: Any, label: str) -> tuple[Any, Any, Any, Any]:
    """Require exactly four raw fields; do not accept scalar/short encodings."""

    if not isinstance(value, (tuple, list)):
        raise DyadicProductInputError(
            f"{label} must be a tuple or list of exactly four fields"
        )
    if len(value) != 4:
        raise DyadicProductInputError(f"{label} must contain exactly four fields")
    return value[0], value[1], value[2], value[3]


def _normalize_fields(
    raw_sign: Any,
    raw_mantissa: Any,
    raw_exponent: Any,
    raw_bitcount: Any,
    *,
    label: str,
    max_mantissa_bits: int,
) -> Endpoint:
    """Validate and canonicalize one finite raw endpoint.

    A zero is finite only in the ordinary ``exponent == bitcount == 0``
    encoding.  mpmath's NaN/Inf encodings use a zero mantissa with a nonzero
    exponent and/or a negative bitcount, and are rejected before arithmetic.
    """

    sign = _index(raw_sign, f"{label}.sign")
    mantissa = _index(raw_mantissa, f"{label}.mantissa")
    exponent = _index(raw_exponent, f"{label}.exponent")
    bitcount = _index(raw_bitcount, f"{label}.bitcount")
    if abs(exponent).bit_length() > MAX_EXPONENT_BITS:
        raise DyadicProductInputError(f"{label}.exponent exceeds resource limit")

    if sign not in (0, 1):
        raise DyadicProductInputError(f"{label}.sign must be 0 or 1")
    if mantissa < 0:
        raise DyadicProductInputError(f"{label}.mantissa must be nonnegative")
    if bitcount < 0:
        raise DyadicProductInputError(
            f"{label}.bitcount must be nonnegative for a finite endpoint"
        )

    if mantissa == 0:
        # Signed zero is finite and canonicalized.  Every nonstandard zero
        # tuple is a special value or an unknown encoding, never ordinary zero.
        if exponent != 0 or bitcount != 0:
            raise DyadicProductInputError(
                f"{label} is a non-finite or unknown zero-mantissa encoding"
            )
        return [0, 0, 0, 0]

    actual_bitcount = mantissa.bit_length()
    if bitcount == 0 or bitcount != actual_bitcount:
        raise DyadicProductInputError(
            f"{label}.bitcount must equal mantissa.bit_length()"
        )
    if actual_bitcount > max_mantissa_bits:
        raise DyadicProductInputError(
            f"{label}.mantissa exceeds max_mantissa_bits"
        )

    # mpmath canonicalization strips trailing zero bits.  This shift is by a
    # mantissa-local count, never by the caller's astronomical exponent.
    trailing = (mantissa & -mantissa).bit_length() - 1
    if trailing:
        mantissa >>= trailing
        exponent += trailing
        bitcount -= trailing
        if abs(exponent).bit_length() > MAX_EXPONENT_BITS:
            raise DyadicProductInputError(f"{label}.normalized exponent exceeds resource limit")
    return [sign, mantissa, exponent, bitcount]


def normalize_endpoint(
    value: Any,
    max_mantissa_bits: int = MAX_MANTISSA_BITS,
) -> Endpoint:
    """Validate and return one canonical finite raw endpoint."""

    limit = _validate_limit(max_mantissa_bits)
    raw = _raw_endpoint(value, "endpoint")
    return _normalize_fields(*raw, label="endpoint", max_mantissa_bits=limit)


def _is_zero(value: Endpoint) -> bool:
    return value[1] == 0


def _compare_normalized(left: Endpoint, right: Endpoint) -> int:
    """Compare two canonical finite endpoints without exponent-sized shifts."""

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

    # The highest occupied binary position is exponent + bitcount.  This
    # comparison remains tiny even when either exponent is 10**17 or larger.
    left_top = left[2] + left[3]
    right_top = right[2] + right[3]
    if left_top != right_top:
        magnitude = 1 if left_top > right_top else -1
    else:
        # Equality of top bits implies shift == right.bitcount-left.bitcount,
        # hence this is bounded by the mantissa resource limits.  Never shift
        # by an exponent difference unless this proof has already held.
        shift = left[2] - right[2]
        if abs(shift) > max(left[3], right[3]):
            raise DyadicProductInputError("internal unbounded comparison shift")
        if shift >= 0:
            left_mantissa = left[1] << shift
            right_mantissa = right[1]
        else:
            left_mantissa = left[1]
            right_mantissa = right[1] << (-shift)
        magnitude = (left_mantissa > right_mantissa) - (
            left_mantissa < right_mantissa
        )

    return magnitude if left[0] == 0 else -magnitude


def compare_dyadics(
    left: Any,
    right: Any,
    max_mantissa_bits: int = MAX_MANTISSA_BITS,
) -> int:
    """Return ``-1``, ``0`` or ``1`` for two finite raw dyadics."""

    limit = _validate_limit(max_mantissa_bits)
    left_n = _normalize_fields(
        *_raw_endpoint(left, "left"),
        label="left",
        max_mantissa_bits=limit,
    )
    right_n = _normalize_fields(
        *_raw_endpoint(right, "right"),
        label="right",
        max_mantissa_bits=limit,
    )
    return _compare_normalized(left_n, right_n)


def _multiply_normalized(
    left: Endpoint,
    right: Endpoint,
    *,
    max_mantissa_bits: int,
) -> Endpoint:
    """Multiply canonical endpoints exactly, with no exponent-sized shift."""

    if _is_zero(left) or _is_zero(right):
        return [0, 0, 0, 0]
    sign = left[0] ^ right[0]
    mantissa = left[1] * right[1]
    exponent = left[2] + right[2]
    product_limit = 2 * max_mantissa_bits
    if mantissa.bit_length() > product_limit:
        raise DyadicProductInputError("exact product mantissa exceeds implementation limit")
    return _normalize_fields(
        sign,
        mantissa,
        exponent,
        mantissa.bit_length(),
        label="product",
        max_mantissa_bits=product_limit,
    )


def multiply_endpoint(
    left: Any,
    right: Any,
    max_mantissa_bits: int = MAX_MANTISSA_BITS,
) -> Endpoint:
    """Return one exact signed product of two finite raw endpoints."""

    limit = _validate_limit(max_mantissa_bits)
    left_n = _normalize_fields(
        *_raw_endpoint(left, "left"),
        label="left",
        max_mantissa_bits=limit,
    )
    right_n = _normalize_fields(
        *_raw_endpoint(right, "right"),
        label="right",
        max_mantissa_bits=limit,
    )
    return _multiply_normalized(left_n, right_n, max_mantissa_bits=limit)


def _validate_interval(value: Any, *, label: str, max_mantissa_bits: int) -> Interval:
    if not isinstance(value, (tuple, list)):
        raise DyadicProductInputError(
            f"{label} must be a two-endpoint raw interval"
        )
    if len(value) != 2:
        raise DyadicProductInputError(f"{label} must contain exactly two endpoints")
    lower = _normalize_fields(
        *_raw_endpoint(value[0], f"{label}.lower"),
        label=f"{label}.lower",
        max_mantissa_bits=max_mantissa_bits,
    )
    upper = _normalize_fields(
        *_raw_endpoint(value[1], f"{label}.upper"),
        label=f"{label}.upper",
        max_mantissa_bits=max_mantissa_bits,
    )
    if _compare_normalized(lower, upper) > 0:
        raise DyadicProductInputError(f"{label} endpoints are reversed")
    return [lower, upper]


def interval_product(
    left: Any,
    right: Any,
    max_mantissa_bits: int = MAX_MANTISSA_BITS,
) -> Interval:
    """Return the exact finite interval product conditional on two enclosures.

    The four products are exact endpoint dyadics.  The lower and upper result
    are selected solely by finite signed comparison; no midpoint, source
    value, addition, exponentiation, integration, or rounding is performed.
    """

    limit = _validate_limit(max_mantissa_bits)
    left_i = _validate_interval(left, label="left", max_mantissa_bits=limit)
    right_i = _validate_interval(right, label="right", max_mantissa_bits=limit)
    products = [
        _multiply_normalized(a, b, max_mantissa_bits=limit)
        for a in left_i
        for b in right_i
    ]
    lower = products[0]
    upper = products[0]
    for candidate in products[1:]:
        if _compare_normalized(candidate, lower) < 0:
            lower = candidate
        if _compare_normalized(candidate, upper) > 0:
            upper = candidate
    return [lower, upper]


# Focused aliases make the operation easy to discover without creating a
# second arithmetic implementation.
def multiply_interval(
    left: Any,
    right: Any,
    max_mantissa_bits: int = MAX_MANTISSA_BITS,
) -> Interval:
    return interval_product(left, right, max_mantissa_bits)


def interval_mul(
    left: Any,
    right: Any,
    max_mantissa_bits: int = MAX_MANTISSA_BITS,
) -> Interval:
    return interval_product(left, right, max_mantissa_bits)


def multiply_intervals(
    left: Any,
    right: Any,
    max_mantissa_bits: int = MAX_MANTISSA_BITS,
) -> Interval:
    return interval_product(left, right, max_mantissa_bits)


def conditional_interval_product(
    left: Any,
    right: Any,
    max_mantissa_bits: int = MAX_MANTISSA_BITS,
) -> Interval:
    return interval_product(left, right, max_mantissa_bits)


def make_kernel(
    max_mantissa_bits: int = MAX_MANTISSA_BITS,
    *,
    left_source_id: Any = LEFT_SOURCE_ID,
    right_source_id: Any = RIGHT_SOURCE_ID,
) -> Callable[[Any, Any], Interval]:
    """Build a bounded kernel while retaining independent caller source IDs."""

    limit = _validate_limit(max_mantissa_bits)
    if not isinstance(left_source_id, str) or not left_source_id:
        raise DyadicProductInputError("left_source_id must be a nonempty string")
    if not isinstance(right_source_id, str) or not right_source_id:
        raise DyadicProductInputError("right_source_id must be a nonempty string")

    def kernel(left: Any, right: Any) -> Interval:
        return interval_product(left, right, limit)

    kernel.__name__ = "conditional_unbounded_dyadic_interval_product"
    kernel.__qualname__ = kernel.__name__
    kernel.max_mantissa_bits = limit  # type: ignore[attr-defined]
    kernel.source_ids = {"left": left_source_id, "right": right_source_id}  # type: ignore[attr-defined]
    kernel.receipt = {  # type: ignore[attr-defined]
        "schema": SCHEMA + "_receipt",
        "operation": "exact signed finite raw-dyadic interval product",
        "left_source_id": left_source_id,
        "right_source_id": right_source_id,
        "output_source_id": OUTPUT_SOURCE_ID,
        "source_ids_separate": left_source_id != right_source_id,
        "conditional_on_input_enclosures": True,
        "input_endpoints_finite_only": True,
        "point_value_selected": False,
        "source_identity_verified": False,
        "original_source_certified": False,
        "max_mantissa_bits": limit,
    }
    return kernel


make_interval_product_kernel = make_kernel


__all__ = [
    "SCHEMA",
    "SOURCE_ID",
    "LEFT_SOURCE_ID",
    "RIGHT_SOURCE_ID",
    "OUTPUT_SOURCE_ID",
    "MAX_MANTISSA_BITS",
    "MAX_PRODUCT_MANTISSA_BITS",
    "MAX_EXPONENT_BITS",
    "DyadicProductInputError",
    "normalize_endpoint",
    "compare_dyadics",
    "multiply_endpoint",
    "interval_product",
    "multiply_interval",
    "interval_mul",
    "multiply_intervals",
    "conditional_interval_product",
    "make_kernel",
    "make_interval_product_kernel",
]
