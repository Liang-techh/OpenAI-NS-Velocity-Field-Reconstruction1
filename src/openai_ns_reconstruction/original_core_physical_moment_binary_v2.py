"""Bounded signed-dyadic physical moment scale ledger.

Input endpoints are exact mpmath tuples ``[sign, mantissa, exponent,
bitcount]``.  Their value is ``(-1)**sign * mantissa * 2**exponent``.  The
arithmetic below never expands a power of two according to its exponent and
rounds every result outward to at most 256 significant binary bits.

Only ordinary Z0/Z1 finite intervals, explicit independent error caps, and a
caller-supplied interval ``L=F0_Z/F0`` are consumed.  Epsilon and F0 are
represented only by their integer powers in the output; no scale value,
midpoint, source identity, Q4, closure, or recursion is selected here.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any


SCHEMA = "original_core_physical_moment_binary_v2"
MOMENTS = ("H", "M", "K", "A", "B", "C")
PHYSICAL_MOMENTS = ("I_V", "I_V2", "I_RF", "I_RFV", "I_RF2", "I_F2")
MAX_BITS = 256
GUARD_BITS = 8
EXACT_SHIFT_LIMIT = 2 * MAX_BITS + GUARD_BITS

_ROWS = {
    "I_V": ("M", 4, (1, 0)),
    "I_V2": ("A", 4, (1, 0)),
    "I_RF": ("H", 8, (2, 1)),
    "I_RFV": ("K", 8, (2, 1)),
    "I_RF2": ("B", 16, (2, 2)),
    "I_F2": ("C", 4, (1, 2)),
}


class BinaryMomentInputError(ValueError):
    """Raised when a signed-dyadic packet is malformed or unsafe."""


@dataclass(frozen=True)
class _D:
    sign: int
    mantissa: int
    exponent: int

    def __post_init__(self) -> None:
        if self.sign not in (0, 1):
            raise BinaryMomentInputError("sign must be 0 or 1")
        if self.mantissa < 0:
            raise BinaryMomentInputError("mantissa must be nonnegative")
        if self.mantissa == 0 and self.sign:
            raise BinaryMomentInputError("zero must have sign 0")

    @property
    def is_zero(self) -> bool:
        return self.mantissa == 0

    @property
    def bits(self) -> int:
        return self.mantissa.bit_length()

    def neg(self) -> "_D":
        return _D(0 if self.is_zero else 1 - self.sign, self.mantissa, self.exponent)

    def endpoint(self) -> list[int]:
        signed = -self.mantissa if self.sign else self.mantissa
        return [signed, self.exponent]


_ZERO = _D(0, 0, 0)


def _int(value: Any, label: str, *, nonnegative: bool = False) -> int:
    if isinstance(value, bool) or type(value) is not int:
        raise BinaryMomentInputError(f"{label} must be an integer")
    if nonnegative and value < 0:
        raise BinaryMomentInputError(f"{label} must be nonnegative")
    return value


def _endpoint(value: Any, label: str) -> _D:
    if isinstance(value, _D):
        return value
    if isinstance(value, float) or not isinstance(value, (tuple, list)):
        raise BinaryMomentInputError(f"{label} must be an exact dyadic endpoint")
    if len(value) == 2:
        signed_mantissa = _int(value[0], f"{label}.signed_mantissa")
        exponent = _int(value[1], f"{label}.exponent")
        if signed_mantissa == 0:
            return _ZERO
        return _D(1 if signed_mantissa < 0 else 0, abs(signed_mantissa), exponent)
    if len(value) != 4:
        raise BinaryMomentInputError(f"{label} must be [sign,mantissa,exponent,bitcount]")
    sign = _int(value[0], f"{label}.sign")
    mantissa = _int(value[1], f"{label}.mantissa", nonnegative=True)
    exponent = _int(value[2], f"{label}.exponent")
    bitcount = _int(value[3], f"{label}.bitcount", nonnegative=True)
    if sign not in (0, 1):
        raise BinaryMomentInputError(f"{label}.sign must be 0 or 1")
    if mantissa == 0:
        if sign != 0 or bitcount != 0:
            raise BinaryMomentInputError(f"{label}: malformed zero")
        return _ZERO
    if bitcount == 0 or bitcount != mantissa.bit_length():
        raise BinaryMomentInputError(f"{label}.bitcount must match mantissa")
    trailing = (mantissa & -mantissa).bit_length() - 1
    if trailing:
        mantissa >>= trailing
        exponent += trailing
    return _D(sign, mantissa, exponent)


def _cmp_abs(left: _D, right: _D) -> int:
    if left.is_zero:
        return 0 if right.is_zero else -1
    if right.is_zero:
        return 1
    left_top = left.bits + left.exponent
    right_top = right.bits + right.exponent
    if left_top != right_top:
        return 1 if left_top > right_top else -1
    shift = left.exponent - right.exponent
    if shift >= 0:
        a, b = left.mantissa << shift, right.mantissa
    else:
        a, b = left.mantissa, right.mantissa << (-shift)
    return (a > b) - (a < b)


def _cmp(left: _D, right: _D) -> int:
    if left.sign != right.sign:
        if left.is_zero and right.is_zero:
            return 0
        return -1 if left.sign else 1
    result = _cmp_abs(left, right)
    return -result if left.sign else result


def _round_mag(mantissa: int, exponent: int, *, upward: bool) -> _D:
    """Round a nonnegative exact dyadic magnitude outward to MAX_BITS."""

    if mantissa == 0:
        return _ZERO
    bits = mantissa.bit_length()
    if bits > MAX_BITS:
        drop = bits - MAX_BITS
        quotient, remainder = divmod(mantissa, 1 << drop)
        if upward and remainder:
            quotient += 1
        exponent += drop
        mantissa = quotient
        if mantissa.bit_length() > MAX_BITS:
            mantissa >>= 1
            exponent += 1
    trailing = (mantissa & -mantissa).bit_length() - 1
    if trailing:
        mantissa >>= trailing
        exponent += trailing
    return _D(0, mantissa, exponent)


def _signed_bounds(sign: int, mantissa: int, exponent: int) -> tuple[_D, _D]:
    lower_mag = _round_mag(mantissa, exponent, upward=False)
    upper_mag = _round_mag(mantissa, exponent, upward=True)
    if sign == 0:
        return lower_mag, upper_mag
    return upper_mag.neg(), lower_mag.neg()


def _mag_plus_unit(value: _D) -> _D:
    """Increase a nonnegative rounded magnitude by one representation unit."""

    if value.is_zero:
        return _D(0, 1, 0)
    return _round_mag(value.mantissa + 1, value.exponent, upward=True)


def _mag_half(value: _D) -> _D:
    if value.is_zero:
        return _ZERO
    return _D(0, value.mantissa, value.exponent - 1)


def _grid_quotient_bounds(value: _D, grid_exponent: int) -> tuple[int, int]:
    """Return floor/ceil of ``value / 2**grid_exponent`` as integers."""

    if value.is_zero:
        return 0, 0
    shift = value.exponent - grid_exponent
    signed_mantissa = -value.mantissa if value.sign else value.mantissa
    if shift >= 0:
        exact = signed_mantissa << shift
        return exact, exact
    divisor_shift = -shift
    # If the denominator is much larger than the mantissa, the quotient is
    # either 0 (positive) or -1 (negative) after directed rounding.  This
    # avoids constructing a 2**|exponent| denominator.
    if divisor_shift > value.bits + 1:
        return (-1, 0) if value.sign else (0, 1)
    divisor = 1 << divisor_shift
    quotient, remainder = divmod(value.mantissa, divisor)
    if value.sign:
        return -quotient - (1 if remainder else 0), -quotient
    return quotient, quotient + (1 if remainder else 0)


def _floor_from_integer(value: int, exponent: int) -> _D:
    if value == 0:
        return _ZERO
    if value > 0:
        return _round_mag(value, exponent, upward=False)
    return _round_mag(-value, exponent, upward=True).neg()


def _ceil_from_integer(value: int, exponent: int) -> _D:
    if value == 0:
        return _ZERO
    if value > 0:
        return _round_mag(value, exponent, upward=True)
    return _round_mag(-value, exponent, upward=False).neg()


def _add_bounds(left: _D, right: _D) -> tuple[_D, _D]:
    """Directed bounds using a fixed grid, never an exponent-sized shift."""

    if left.is_zero:
        return _signed_bounds(right.sign, right.mantissa, right.exponent)
    if right.is_zero:
        return _signed_bounds(left.sign, left.mantissa, left.exponent)
    top = max(left.bits + left.exponent, right.bits + right.exponent)
    grid_exponent = top - MAX_BITS - 4
    left_floor, left_ceil = _grid_quotient_bounds(left, grid_exponent)
    right_floor, right_ceil = _grid_quotient_bounds(right, grid_exponent)
    lower_integer = left_floor + right_floor
    upper_integer = left_ceil + right_ceil
    return (
        _floor_from_integer(lower_integer, grid_exponent),
        _ceil_from_integer(upper_integer, grid_exponent),
    )


def _mul_bounds(left: _D, right: _D) -> tuple[_D, _D]:
    if left.is_zero or right.is_zero:
        return _ZERO, _ZERO
    sign = left.sign ^ right.sign
    return _signed_bounds(sign, left.mantissa * right.mantissa, left.exponent + right.exponent)


def _iv_add(left: tuple[_D, _D], right: tuple[_D, _D]) -> tuple[_D, _D]:
    return _add_bounds(left[0], right[0])[0], _add_bounds(left[1], right[1])[1]


def _iv_mul(left: tuple[_D, _D], right: tuple[_D, _D]) -> tuple[_D, _D]:
    products = [_mul_bounds(a, b) for a in left for b in right]
    lower = products[0][0]
    upper = products[0][1]
    for candidate in (item[0] for item in products):
        if _cmp(candidate, lower) < 0:
            lower = candidate
    for candidate in (item[1] for item in products):
        if _cmp(candidate, upper) > 0:
            upper = candidate
    return lower, upper


def _iv_scale(value: tuple[_D, _D], factor: int) -> tuple[_D, _D]:
    if factor < 0:
        return _neg_interval(_iv_scale(value, -factor))
    factor_d = _D(0, factor, 0) if factor else _ZERO
    return _iv_mul(value, (factor_d, factor_d))


def _neg_interval(value: tuple[_D, _D]) -> tuple[_D, _D]:
    return value[1].neg(), value[0].neg()


def _point(value: _D) -> tuple[_D, _D]:
    return value, value


def _max_abs(interval: tuple[_D, _D]) -> _D:
    left, right = interval
    return left if _cmp_abs(left, right) >= 0 else right


def _interval(value: Any, label: str) -> tuple[_D, _D]:
    if isinstance(value, Mapping):
        lower_raw = value.get("lower", value.get("lo"))
        upper_raw = value.get("upper", value.get("hi"))
        lower = _endpoint(lower_raw, f"{label}.lower")
        upper = _endpoint(upper_raw, f"{label}.upper")
    elif isinstance(value, (tuple, list)) and len(value) == 2 and all(type(x) is int for x in value):
        lower = upper = _endpoint(value, label)
    elif isinstance(value, (tuple, list)) and len(value) == 2:
        lower = _endpoint(value[0], f"{label}.lower")
        upper = _endpoint(value[1], f"{label}.upper")
    else:
        point = _endpoint(value, label)
        lower = point
        upper = point
    if _cmp(lower, upper) > 0:
        raise BinaryMomentInputError(f"{label} endpoints are reversed")
    return lower, upper


def _error(value: Any, label: str) -> _D:
    if isinstance(value, Mapping):
        nested = value.get("error_upper", value.get("error", value.get("upper")))
        if nested is None:
            raise BinaryMomentInputError(f"{label} requires error_upper")
        value = nested
    if isinstance(value, (tuple, list)) and len(value) == 2 and all(
        isinstance(item, (tuple, list)) for item in value
    ):
        lower, upper = _interval(value, label)
        if lower.sign or upper.sign:
            raise BinaryMomentInputError(f"{label} must be nonnegative")
        return upper
    result = _endpoint(value, label)
    if result.sign:
        raise BinaryMomentInputError(f"{label} must be nonnegative")
    return result


def _moment_name(value: Any) -> str | None:
    if not isinstance(value, str):
        return None
    return {"h": "H", "m": "M", "k": "K", "a": "A", "b": "B", "c": "C"}.get(value.strip().lower())


def _raw_zone(row: Any, zone: str) -> Any:
    if not isinstance(row, Mapping):
        return None
    names = (
        ("z0", "Z0", "finite_Z0", "finite_z0")
        if zone == "Z0"
        else ("z1", "Z1", "finite_Z1", "finite_z1")
    )
    raw = next((row[name] for name in names if name in row), None)
    if raw is None and isinstance(row.get("finite"), Mapping):
        raw = next((row["finite"][name] for name in names if name in row["finite"]), None)
    return raw


def _zone_finite(row: Any, moment: str, zone: str) -> tuple[_D, _D]:
    raw = _raw_zone(row, zone)
    if raw is None:
        raise BinaryMomentInputError(f"{moment}.{zone}.finite is required")
    if isinstance(raw, Mapping):
        raw = raw.get("finite", raw.get("interval", raw.get("bounds")))
    if raw is None:
        raise BinaryMomentInputError(f"{moment}.{zone}.finite is required")
    return _interval(raw, f"{moment}.{zone}.finite")


def _zone_error(row: Any, moment: str, zone: str) -> _D:
    raw = _raw_zone(row, zone)
    if isinstance(raw, Mapping):
        error = raw.get("error_upper", raw.get("error"))
        if error is not None:
            return _error(error, f"{moment}.{zone}.error_upper")
    names = (
        ("error_Z0", "error_z0", "error0", "z0_error")
        if zone == "Z0"
        else ("error_Z1", "error_z1", "error1", "z1_error")
    )
    for name in names:
        if name in row:
            return _error(row[name], f"{moment}.{zone}.error_upper")
    raise BinaryMomentInputError(f"{moment}.{zone}.error_upper is required")


def _zone_family(row: Any, zone: str) -> Any:
    raw = _raw_zone(row, zone)
    if isinstance(raw, Mapping):
        return raw.get("source_family", raw.get("source_id"))
    return None


def _row_family(row: Any) -> Any:
    if isinstance(row, Mapping):
        return row.get("source_family", row.get("source_id"))
    return None


def _source_value(source_families: Any, moment: str) -> Any:
    if not isinstance(source_families, Mapping):
        return None
    value = source_families.get(moment, source_families.get(moment.lower()))
    if isinstance(value, Mapping):
        return value.get("source_family", value.get("source_id"))
    return value


def _same_family(values: list[Any]) -> Any:
    present = [item for item in values if item is not None]
    if not present:
        return None
    first = present[0]
    if any(item != first for item in present[1:]):
        raise BinaryMomentInputError("source family labels must all match")
    return first


def _log_interval(value: Any) -> tuple[tuple[_D, _D], Any]:
    family = None
    if isinstance(value, Mapping):
        family = value.get("source_family", value.get("source_id"))
        value = value.get("finite", value.get("L", value.get("interval", value)))
    return _interval(value, "L"), family


def _public_interval(value: tuple[_D, _D]) -> dict[str, list[int]]:
    return {"lower": value[0].endpoint(), "upper": value[1].endpoint()}


def _public_error(value: _D) -> list[int]:
    return value.endpoint()


def _widen(value: tuple[_D, _D], error: _D) -> tuple[_D, _D]:
    lower = _add_bounds(value[0], error.neg())[0]
    upper = _add_bounds(value[1], error)[1]
    return lower, upper


def interval_add(left: Any, right: Any) -> dict[str, list[int]]:
    """Public directed interval addition for exact dyadic endpoints."""

    return _public_interval(_iv_add(_interval(left, "left"), _interval(right, "right")))


def interval_mul(left: Any, right: Any) -> dict[str, list[int]]:
    """Public directed interval multiplication for exact dyadic endpoints."""

    return _public_interval(_iv_mul(_interval(left, "left"), _interval(right, "right")))


def interval_scale(value: Any, factor: int) -> dict[str, list[int]]:
    """Public directed integer scaling, including negative integers."""

    factor = _int(factor, "factor")
    return _public_interval(_iv_scale(_interval(value, "value"), factor))


def interval_maxabs(value: Any) -> list[int]:
    """Return the nonnegative endpoint with greatest absolute magnitude."""

    interval = _interval(value, "value")
    selected = interval[0] if _cmp_abs(interval[0], interval[1]) >= 0 else interval[1]
    return _D(0, selected.mantissa, selected.exponent).endpoint()


def scale(
    moments: Mapping[str, Any],
    log_derivative: Any,
    *,
    source_family: Any = None,
    source_families: Any = None,
    errors: Any = None,
) -> dict[str, Any]:
    """Return JSON-ready finite/error/enclosure coefficient rows.

    Each ``moments[name]['z0'|'z1']`` row must contain a ``finite`` interval
    ``[[sign,mantissa,exponent,bitcount], [sign,mantissa,exponent,bitcount]]``
    and an explicit ``error_upper`` endpoint tuple.  The second argument is
    the same two-endpoint interval for ``L=F0_Z/F0``.  Every output endpoint
    is ``[signed_mantissa, exponent]``; no exponent-sized integer is created.
    """

    if not isinstance(moments, Mapping):
        raise BinaryMomentInputError("moments must be a mapping")
    log_interval, log_family = _log_interval(log_derivative)
    sup_abs_raw = log_interval[0] if _cmp_abs(log_interval[0], log_interval[1]) >= 0 else log_interval[1]
    sup_abs = _D(0, sup_abs_raw.mantissa, sup_abs_raw.exponent)
    families: list[Any] = [source_family, log_family]
    ordinary: dict[str, tuple[tuple[_D, _D], tuple[_D, _D], _D, _D]] = {}
    for raw_name in MOMENTS:
        row = moments.get(raw_name, moments.get(raw_name.lower()))
        if row is None:
            raise BinaryMomentInputError(f"missing moment {raw_name}")
        z0 = _zone_finite(row, raw_name, "Z0")
        z1 = _zone_finite(row, raw_name, "Z1")
        e0 = _zone_error(row, raw_name, "Z0")
        e1 = _zone_error(row, raw_name, "Z1")
        ordinary[raw_name] = (z0, z1, e0, e1)
        families.extend([_row_family(row), _zone_family(row, "Z0"), _zone_family(row, "Z1"), _source_value(source_families, raw_name)])
        if isinstance(errors, Mapping):
            external = errors.get(raw_name, errors.get(raw_name.lower()))
            if isinstance(external, Mapping):
                for zone, current in (("Z0", e0), ("Z1", e1)):
                    extra = external.get(zone, external.get(zone.lower()))
                    if extra is not None:
                        addition = _error(extra, f"{raw_name}.{zone}.external_error")
                        if zone == "Z0":
                            e0 = _add_bounds(current, addition)[1]
                        else:
                            e1 = _add_bounds(current, addition)[1]
                ordinary[raw_name] = (z0, z1, e0, e1)
            elif external is not None:
                raise BinaryMomentInputError("external errors must specify Z0 or Z1")
    family = _same_family(families)

    finite: dict[str, Any] = {}
    finite_z1: dict[str, Any] = {}
    error_upper: dict[str, Any] = {}
    error_upper_z1: dict[str, Any] = {}
    enclosure: dict[str, Any] = {}
    enclosure_z1: dict[str, Any] = {}
    powers: dict[str, list[int]] = {}
    records: dict[str, Any] = {}

    for physical_name in PHYSICAL_MOMENTS:
        moment, numerical, power_pair = _ROWS[physical_name]
        z0, z1, e0, e1 = ordinary[moment]
        factor = _iv_scale(z0, numerical)
        derivative_inside = _iv_add(z1, _iv_scale(_iv_mul(log_interval, z0), power_pair[1]))
        derivative = _iv_scale(derivative_inside, numerical)
        term = _iv_mul(_point(sup_abs), _point(e0))
        if power_pair[1]:
            term = _iv_scale(term, power_pair[1])
        derivative_error = _add_bounds(e1, term[1] if power_pair[1] else _ZERO)[1]
        derivative_error = _iv_scale(_point(derivative_error), numerical)[1]
        value_error = _iv_scale(_point(e0), numerical)[1]
        finite[physical_name] = _public_interval(factor)
        finite_z1[physical_name] = _public_interval(derivative)
        error_upper[physical_name] = _public_error(value_error)
        error_upper_z1[physical_name] = _public_error(derivative_error)
        enclosure[physical_name] = _public_interval(_widen(factor, value_error))
        enclosure_z1[physical_name] = _public_interval(_widen(derivative, derivative_error))
        powers[physical_name] = [power_pair[0], power_pair[1]]
        records[physical_name] = {
            "ordinary_moment": moment,
            "numerical_coefficient": numerical,
            "epsilon_power": power_pair[0],
            "F0_power": power_pair[1],
            "scale": powers[physical_name],
            "finite": finite[physical_name],
            "error_upper": error_upper[physical_name],
            "enclosure": enclosure[physical_name],
            "finite_Z1": finite_z1[physical_name],
            "error_upper_Z1": error_upper_z1[physical_name],
            "enclosure_Z1": enclosure_z1[physical_name],
        }

    finite["Z1"] = finite_z1
    error_upper["Z1"] = error_upper_z1
    enclosure["Z1"] = enclosure_z1
    for name in PHYSICAL_MOMENTS:
        finite[f"{name}_Z1"] = finite_z1[name]
        error_upper[f"{name}_Z1"] = error_upper_z1[name]
        enclosure[f"{name}_Z1"] = enclosure_z1[name]

    return {
        "schema": SCHEMA,
        "L": _public_interval(log_interval),
        "sup_abs_L": _public_error(sup_abs),
        "source_family": family,
        "source_family_consistent": True,
        "finite": finite,
        "error_upper": error_upper,
        "enclosure": enclosure,
        "total": {**enclosure},
        "finite_Z1": finite_z1,
        "error_upper_Z1": error_upper_z1,
        "enclosure_Z1": enclosure_z1,
        "scale": powers,
        "scale_powers": powers,
        "records": records,
        "precision_bits": MAX_BITS,
        "assumptions": {
            "directed_dyadic_rounding": True,
            "epsilon_and_F0_values_materialized": False,
            "epsilon_Z_constant_is_caller_premise": True,
            "F0_positive_is_caller_premise": True,
            "L_is_caller_supplied_interval": True,
            "input_errors_explicit_and_independent": True,
            "old_errors_cancelled": False,
            "source_identity_certified": False,
            "q4_certified": False,
            "closure_proved": False,
            "time_recursion_proved": False,
            "midpoint_selected": False,
        },
    }


evaluate = scale
scale.interval_add = interval_add
scale.interval_mul = interval_mul
scale.interval_scale = interval_scale
scale.interval_maxabs = interval_maxabs


__all__ = [
    "SCHEMA",
    "MOMENTS",
    "PHYSICAL_MOMENTS",
    "MAX_BITS",
    "BinaryMomentInputError",
    "interval_add",
    "interval_mul",
    "interval_scale",
    "interval_maxabs",
    "scale",
    "evaluate",
]



