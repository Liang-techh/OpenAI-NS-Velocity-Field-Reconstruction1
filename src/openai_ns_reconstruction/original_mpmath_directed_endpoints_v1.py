"""Strict exact conversion of mpmath 1.3.x interval endpoints.

``mpmath.iv`` stores a directed interval on a private ``_mpi_`` member.
For the 1.3.x representation each endpoint is a four-tuple
``(sign, mantissa, exponent, bitcount)`` and denotes

``(-1)**sign * mantissa * 2**exponent``.

This module is intentionally only a representation boundary. It does not
import mpmath, inspect a provider, inspect a source graph, or infer a rounded
endpoint from a scalar. Only the explicitly recognised 1.3.0 spelling of the
1.3.x contract is accepted; callers can pass the returned callable directly to
``original_q4_fraction_interval_v2.extract_component_interval``.
"""

from __future__ import annotations

import operator
from fractions import Fraction
from typing import Any, Callable


SCHEMA = "original_mpmath_directed_endpoints_v1"
RECEIPT_SCHEMA = SCHEMA + "_receipt"
MPMATH_VERSION_FAMILY = "1.3.x"
SUPPORTED_VERSION = "1.3.0"
SUPPORTED_VERSIONS = (SUPPORTED_VERSION,)

# The default is deliberately finite. It is large enough for ordinary
# mpmath 1.3.x precision settings while keeping an accidental untrusted input
# from requesting an unbounded Python integer shift.
DEFAULT_MAX_BITS = 4096
_ABSOLUTE_MAX_BITS = 1_000_000
_MISSING = object()


def _index(value: Any, label: str) -> int:
    """Return an exact index value, rejecting bool and scalar coercions."""

    # ``bool`` is an ``int`` subclass and also has an index implementation;
    # neither True nor False is a valid field in this wire representation.
    if isinstance(value, bool):
        raise TypeError(f"{label} must be an integer index, not bool")
    try:
        result = operator.index(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise TypeError(f"{label} must support operator.index") from exc
    if isinstance(result, bool):
        raise TypeError(f"{label} must be an integer index, not bool")
    # CPython's operator.index currently returns ``int`` even for an int-like
    # object. Keep this check explicit so a non-conforming implementation
    # cannot smuggle another scalar type into arithmetic below.
    if type(result) is not int:
        raise TypeError(f"{label} did not produce a plain integer index")
    return result


def _validate_version(version: Any) -> str:
    """Accept only the known mpmath 1.3.x contract spelling."""

    if not isinstance(version, str):
        raise TypeError("version must be the exact string '1.3.0'")
    if version not in SUPPORTED_VERSIONS:
        raise ValueError(
            "unsupported mpmath version; only the recognised 1.3.x contract "
            "('1.3.0') is accepted"
        )
    return version


def _validate_max_bits(max_bits: Any) -> int:
    """Validate the finite resource budget before decoding any endpoint."""

    bits = _index(max_bits, "max_bits")
    if bits <= 0:
        raise ValueError("max_bits must be positive")
    if bits > _ABSOLUTE_MAX_BITS:
        raise ValueError("max_bits exceeds the implementation resource cap")
    return bits


def _raw_endpoint(value: Any, label: str) -> tuple[Any, Any, Any, Any]:
    """Validate the structural four-field endpoint without numeric coercion."""

    if not isinstance(value, (tuple, list)):
        raise TypeError(f"{label} must be a tuple or list of exactly four fields")
    if len(value) != 4:
        raise ValueError(f"{label} must contain exactly four fields")
    # Indexing after the exact length check avoids accepting a scalar iterable,
    # and avoids traversing an arbitrary unbounded iterator.
    return value[0], value[1], value[2], value[3]


def _endpoint_fraction(value: Any, *, max_bits: int, label: str) -> Fraction:
    """Decode one finite 1.3.x endpoint to an exact :class:`Fraction`."""

    raw_sign, raw_mantissa, raw_exponent, raw_bitcount = _raw_endpoint(value, label)
    sign = _index(raw_sign, f"{label}.sign")
    mantissa = _index(raw_mantissa, f"{label}.mantissa")
    exponent = _index(raw_exponent, f"{label}.exponent")
    bitcount = _index(raw_bitcount, f"{label}.bitcount")

    if sign not in (0, 1):
        raise ValueError(f"{label}.sign must be 0 or 1")
    if mantissa < 0:
        raise ValueError(f"{label}.mantissa must be nonnegative")
    if bitcount < 0:
        # This also rejects mpmath's negative bit-count sentinels for NaN and
        # infinities; they are never allowed to collapse into zero.
        raise ValueError(f"{label}.bitcount must be nonnegative")

    if mantissa == 0:
        # The sole accepted zero encoding. Every other zero-mantissa code is
        # a special value or an unknown encoding, and must not silently become
        # Fraction(0).
        if exponent == 0 and bitcount == 0:
            return Fraction(0, 1)
        raise ValueError(f"{label} is a non-finite or unknown zero-mantissa code")

    actual_bitcount = mantissa.bit_length()
    if bitcount != actual_bitcount:
        raise ValueError(
            f"{label}.bitcount does not match mantissa.bit_length()"
        )
    if actual_bitcount > max_bits:
        raise ValueError(f"{label}.mantissa exceeds max_bits")

    # Reject an exponent whose magnitude could make a shift large before
    # calculating numerator_bits or denominator_bits. bit_length() is a
    # bounded-size inspection; the negation below happens only after the
    # exponent has passed this guard.
    max_bits_bit_length = max_bits.bit_length()
    if exponent >= 0:
        if exponent.bit_length() > max_bits_bit_length or exponent > max_bits:
            raise ValueError(f"{label}.exponent exceeds max_bits")
        exponent_magnitude = exponent
    else:
        if exponent.bit_length() > max_bits_bit_length:
            raise ValueError(f"{label}.exponent exceeds max_bits")
        exponent_magnitude = -exponent
        if exponent_magnitude > max_bits:
            raise ValueError(f"{label}.exponent exceeds max_bits")

    # These are the actual bit lengths of the unreduced numerator and
    # denominator that the exact construction would allocate. Check them
    # before any left shift.
    if exponent >= 0:
        numerator_bits = actual_bitcount + exponent
        denominator_bits = 1
    else:
        numerator_bits = actual_bitcount
        denominator_bits = exponent_magnitude + 1
    if numerator_bits > max_bits or denominator_bits > max_bits:
        raise ValueError(f"{label} numerator/denominator exceeds max_bits")

    if exponent >= 0:
        numerator = mantissa << exponent
        denominator = 1
    else:
        numerator = mantissa
        denominator = 1 << exponent_magnitude
    result = Fraction(numerator, denominator)
    return -result if sign else result


def _convert_validated(value: Any, *, max_bits: int) -> tuple[Fraction, Fraction]:
    """Convert an already version- and budget-validated interval object."""

    raw = getattr(value, "_mpi_", _MISSING)
    if raw is _MISSING:
        raise TypeError("value must expose an mpmath 1.3.x _mpi_ member")
    if not isinstance(raw, (tuple, list)):
        raise TypeError("value._mpi_ must be a tuple or list of exactly two endpoints")
    if len(raw) != 2:
        raise ValueError("value._mpi_ must contain exactly two endpoints")

    lower = _endpoint_fraction(raw[0], max_bits=max_bits, label="lower")
    upper = _endpoint_fraction(raw[1], max_bits=max_bits, label="upper")
    if lower > upper:
        # Directed endpoint order is evidence supplied by the caller. Never
        # repair a reversed enclosure by swapping its endpoints.
        raise ValueError("directed lower endpoint exceeds upper endpoint")
    return lower, upper


def _make_receipt(version: str, max_bits: int) -> dict[str, Any]:
    """Create a caller-retainable scope receipt for one converter."""

    return {
        "schema": RECEIPT_SCHEMA,
        "conversion_semantics": (
            "mpmath 1.3.x private interval._mpi_ four-field endpoint "
            "contract only"
        ),
        "version_family": MPMATH_VERSION_FAMILY,
        "version": version,
        "version_supplied_by_caller": True,
        "runtime_version_verified": False,
        "runtime_version_observed": None,
        "source_provenance_unknown": True,
        "max_bits": max_bits,
        "conversion_contract_recognized": True,
        "provider_verified": False,
        "source_verified": False,
        "rounding_verified": False,
        "actual_provider_verified": False,
        "source_identity_verified": False,
        "endpoint_rounding_verified": False,
        "actual_provider_called": False,
        "source_provider_called": False,
        "actual_value_selected": False,
        "value_selected": False,
        "source_values_selected": False,
        "limitations": [
            "The caller supplied the version; the installed runtime version is not verified.",
            "The converter authenticates only the recognised 1.3.x representation contract.",
            "The actual provider, source identity, and endpoint rounding are unverified.",
            "No numerical value is selected from the enclosure.",
        ],
    }


def conversion_receipt(
    version: str = SUPPORTED_VERSION,
    max_bits: int = DEFAULT_MAX_BITS,
) -> dict[str, Any]:
    """Return a fresh metadata receipt for a recognised converter contract."""

    checked_version = _validate_version(version)
    checked_bits = _validate_max_bits(max_bits)
    return _make_receipt(checked_version, checked_bits)


def make_converter(
    version: str = SUPPORTED_VERSION,
    max_bits: int = DEFAULT_MAX_BITS,
) -> Callable[[Any], tuple[Fraction, Fraction]]:
    """Build a one-argument exact endpoint converter.

    The returned callable is suitable for
    ``extract_component_interval(component, converter)``. Its ``receipt``,
    ``metadata_receipt`` and ``conversion_receipt`` attributes all point to a
    caller-retainable dictionary describing the limited representation scope.
    """

    checked_version = _validate_version(version)
    checked_bits = _validate_max_bits(max_bits)
    receipt = _make_receipt(checked_version, checked_bits)

    def converter(value: Any) -> tuple[Fraction, Fraction]:
        return _convert_validated(value, max_bits=checked_bits)

    converter.__name__ = "convert_mpmath_1_3_x_interval"
    converter.__qualname__ = converter.__name__
    converter.version = checked_version  # type: ignore[attr-defined]
    converter.max_bits = checked_bits  # type: ignore[attr-defined]
    converter.receipt = receipt  # type: ignore[attr-defined]
    converter.metadata_receipt = receipt  # type: ignore[attr-defined]
    converter.conversion_receipt = receipt  # type: ignore[attr-defined]
    return converter


def convert(
    value: Any,
    *,
    version: str = SUPPORTED_VERSION,
    max_bits: int = DEFAULT_MAX_BITS,
) -> tuple[Fraction, Fraction]:
    """Convert one mpmath 1.3.x interval object to Fraction endpoints."""

    _validate_version(version)
    checked_bits = _validate_max_bits(max_bits)
    return _convert_validated(value, max_bits=checked_bits)


# Descriptive aliases keep the boundary easy to discover without changing its
# one-argument callback shape.
convert_interval = convert
convert_endpoints = convert
scalar_converter = make_converter
DEFAULT_RECEIPT = conversion_receipt()
RECEIPT = DEFAULT_RECEIPT


# The default conversion function also has a receipt for callers that choose
# the direct ``convert(value)`` API rather than constructing a closure.
convert.receipt = DEFAULT_RECEIPT  # type: ignore[attr-defined]
convert.metadata_receipt = DEFAULT_RECEIPT  # type: ignore[attr-defined]


__all__ = [
    "SCHEMA",
    "RECEIPT_SCHEMA",
    "MPMATH_VERSION_FAMILY",
    "SUPPORTED_VERSION",
    "SUPPORTED_VERSIONS",
    "DEFAULT_MAX_BITS",
    "DEFAULT_RECEIPT",
    "RECEIPT",
    "conversion_receipt",
    "make_converter",
    "scalar_converter",
    "convert",
    "convert_interval",
    "convert_endpoints",
]
