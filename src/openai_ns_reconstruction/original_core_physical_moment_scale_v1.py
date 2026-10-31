"""Exact Fraction conversion from six ordinary moments to physical rows.

This module is an algebraic scale ledger, not an evaluator of a physical
field.  It consumes caller-supplied ordinary ``Z0``/``Z1`` finite intervals,
independent nonnegative error caps, and the interval
``L = F0_Z/F0``.  It returns coefficient-space finite/error/enclosure rows and
the positive scale powers ``(epsilon_power, F0_power)`` separately.  No
epsilon or F0 value is materialized, no midpoint is selected, and no
integration, recursion, or ancestor helper is called.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from fractions import Fraction
from typing import Any


SCHEMA = "original_core_physical_moment_scale_v1"
MOMENTS = ("H", "M", "K", "A", "B", "C")
PHYSICAL_MOMENTS = ("I_V", "I_V2", "I_RF", "I_RFV", "I_RF2", "I_F2")

# physical row: ordinary moment, numerical coefficient, (epsilon power, F0 power)
_ROWS = {
    "I_V": ("M", 4, (1, 0)),
    "I_V2": ("A", 4, (1, 0)),
    "I_RF": ("H", 8, (2, 1)),
    "I_RFV": ("K", 8, (2, 1)),
    "I_RF2": ("B", 16, (2, 2)),
    "I_F2": ("C", 4, (1, 2)),
}


class PhysicalMomentScaleInputError(ValueError):
    """Raised when a scale packet is not exact and admissible."""


@dataclass(frozen=True)
class RationalInterval:
    """An ordered interval with exact Fraction endpoints."""

    lower: Fraction
    upper: Fraction

    def __post_init__(self) -> None:
        if type(self.lower) is not Fraction or type(self.upper) is not Fraction:
            raise TypeError("interval endpoints must be Fraction values")
        if self.lower > self.upper:
            raise ValueError("interval lower endpoint exceeds upper endpoint")

    @property
    def lo(self) -> Fraction:
        return self.lower

    @property
    def hi(self) -> Fraction:
        return self.upper

    def __iter__(self):
        yield self.lower
        yield self.upper

    def __len__(self) -> int:
        return 2

    def __getitem__(self, index: int) -> Fraction:
        if index == 0:
            return self.lower
        if index == 1:
            return self.upper
        raise IndexError(index)

    def __eq__(self, other: object) -> bool:
        if isinstance(other, RationalInterval):
            return self.lower == other.lower and self.upper == other.upper
        if isinstance(other, (tuple, list)) and len(other) == 2:
            return self.lower == other[0] and self.upper == other[1]
        if isinstance(other, Mapping):
            return self.lower == other.get("lower", other.get("lo")) and self.upper == other.get("upper", other.get("hi"))
        return NotImplemented

    def as_tuple(self) -> tuple[Fraction, Fraction]:
        return self.lower, self.upper


_MISSING = object()


def _lookup(mapping: Mapping[Any, Any], *names: str) -> Any:
    for name in names:
        if name in mapping:
            return mapping[name]
    return _MISSING


def _fraction(value: Any, label: str) -> Fraction:
    if isinstance(value, bool):
        raise PhysicalMomentScaleInputError(f"{label} must not be bool")
    if isinstance(value, float):
        raise PhysicalMomentScaleInputError(f"{label} must not be float")
    if type(value) is Fraction:
        return value
    if type(value) is int:
        return Fraction(value)
    if isinstance(value, str):
        text = value.strip()
        if not text:
            raise PhysicalMomentScaleInputError(f"{label} is empty")
        try:
            return Fraction(text)
        except (TypeError, ValueError, ZeroDivisionError) as exc:
            raise PhysicalMomentScaleInputError(f"{label} must be an exact rational") from exc
    raise PhysicalMomentScaleInputError(f"{label} must be an exact rational")


def _interval(value: Any, label: str) -> RationalInterval:
    if isinstance(value, RationalInterval):
        return value
    if isinstance(value, Mapping):
        lower = _lookup(value, "lower", "lo", "left")
        upper = _lookup(value, "upper", "hi", "right")
        if lower is _MISSING or upper is _MISSING:
            raise PhysicalMomentScaleInputError(f"{label} needs lower and upper")
        lo = _fraction(lower, f"{label}.lower")
        hi = _fraction(upper, f"{label}.upper")
    elif isinstance(value, (tuple, list)) and len(value) == 2:
        lo = _fraction(value[0], f"{label}.lower")
        hi = _fraction(value[1], f"{label}.upper")
    else:
        lo = _fraction(value, label)
        hi = lo
    if lo > hi:
        raise PhysicalMomentScaleInputError(f"{label} is reversed")
    return RationalInterval(lo, hi)


def _error(value: Any, label: str) -> Fraction:
    """Decode a nonnegative error cap without cancelling any old error."""

    if isinstance(value, Mapping):
        nested = _lookup(value, "error", "error_upper", "upper", "bound", "tail")
        if nested is _MISSING:
            if not value:
                return Fraction(0)
            return sum((_error(item, f"{label}[{key!r}]") for key, item in value.items()), Fraction(0))
        value = nested
    elif isinstance(value, (tuple, list)) and len(value) == 2:
        # An endpoint error interval denotes its upper cap; both endpoints must
        # themselves be nonnegative.
        interval = _interval(value, label)
        if interval.lower < 0:
            raise PhysicalMomentScaleInputError(f"{label} must be nonnegative")
        return interval.upper
    elif isinstance(value, (tuple, list)):
        return sum((_error(item, f"{label}[{index}]") for index, item in enumerate(value)), Fraction(0))
    result = _fraction(value, label)
    if result < 0:
        raise PhysicalMomentScaleInputError(f"{label} must be nonnegative")
    return result


def _iv_scale(value: RationalInterval, factor: Fraction) -> RationalInterval:
    if factor >= 0:
        return RationalInterval(value.lower * factor, value.upper * factor)
    return RationalInterval(value.upper * factor, value.lower * factor)


def _iv_add(left: RationalInterval, right: RationalInterval) -> RationalInterval:
    return RationalInterval(left.lower + right.lower, left.upper + right.upper)


def _iv_mul(left: RationalInterval, right: RationalInterval) -> RationalInterval:
    values = (
        left.lower * right.lower,
        left.lower * right.upper,
        left.upper * right.lower,
        left.upper * right.upper,
    )
    return RationalInterval(min(values), max(values))


def _moment_name(value: Any) -> str | None:
    if not isinstance(value, str):
        return None
    return {"h": "H", "m": "M", "k": "K", "a": "A", "b": "B", "c": "C"}.get(value.strip().lower())


def _zone_name(value: Any) -> str | None:
    if not isinstance(value, str):
        return None
    key = value.strip().lower().replace("_", "")
    if key in ("z0", "zero", "finite0", "j", "value"):
        return "Z0"
    if key in ("z1", "one", "finite1", "jz", "derivative", "ordinaryz1"):
        return "Z1"
    return None


def _raw_zone(row: Any, zone: str) -> Any:
    if not isinstance(row, Mapping):
        return _MISSING
    names = (
        ("z0", "Z0", "finite_Z0", "finite_z0", "finite0", "J", "value")
        if zone == "Z0"
        else ("z1", "Z1", "finite_Z1", "finite_z1", "finite1", "J_Z", "derivative", "ordinary_Z1")
    )
    direct = _lookup(row, *names)
    if direct is _MISSING:
        finite = _lookup(row, "finite", "finite_terms", "ordinary")
        if isinstance(finite, Mapping):
            direct = _lookup(finite, *names)
    return direct


def _zone_interval(row: Any, moment: str, zone: str) -> RationalInterval | Any:
    raw = _raw_zone(row, zone)
    if raw is _MISSING:
        return _MISSING
    if isinstance(raw, Mapping):
        nested = _lookup(raw, "finite", "interval", "bounds", "value")
        if nested is not _MISSING:
            raw = nested
    return _interval(raw, f"{moment}.{zone}")


def _zone_error(row: Any, moment: str, zone: str) -> Fraction | Any:
    if not isinstance(row, Mapping):
        return _MISSING
    raw = _raw_zone(row, zone)
    if isinstance(raw, Mapping):
        nested = _lookup(raw, "error", "error_upper", "tail", "bound")
        if nested is not _MISSING:
            return _error(nested, f"{moment}.{zone}.error")
    names = (
        ("error_Z0", "error_z0", "error0", "z0_error", "e0")
        if zone == "Z0"
        else ("error_Z1", "error_z1", "error1", "z1_error", "e1")
    )
    direct = _lookup(row, *names)
    if direct is not _MISSING:
        return _error(direct, f"{moment}.{zone}.error")
    errors = _lookup(row, "errors", "error_upper", "tail_errors")
    if isinstance(errors, Mapping):
        nested = _lookup(errors, *("Z0", "z0") if zone == "Z0" else ("Z1", "z1"))
        if nested is not _MISSING:
            return _error(nested, f"{moment}.{zone}.error")
    return _MISSING


def _zone_family(row: Any, zone: str) -> Any:
    raw = _raw_zone(row, zone)
    if isinstance(raw, Mapping):
        return _lookup(raw, "source_family", "source_id", "family")
    return _MISSING


def _row_family(row: Any) -> Any:
    if not isinstance(row, Mapping):
        return _MISSING
    return _lookup(row, "source_family", "source_id", "family")


def _map_value(mapping: Any, moment: str, zone: str) -> Any:
    if not isinstance(mapping, Mapping):
        return _MISSING
    direct = _lookup(mapping, moment, moment.lower())
    if direct is _MISSING:
        return _MISSING
    if isinstance(direct, Mapping):
        direct_zone = _lookup(
            direct,
            *("Z0", "z0", "finite_Z0", "finite_z0")
            if zone == "Z0"
            else ("Z1", "z1", "finite_Z1", "finite_z1")
        )
        if direct_zone is not _MISSING:
            return direct_zone
    return direct


def _validate_families(values: list[Any]) -> Any:
    present = [value for value in values if value is not _MISSING and value is not None]
    if not present:
        return None
    first = present[0]
    try:
        consistent = all(bool(first == value) for value in present[1:])
    except Exception:
        consistent = all(first is value for value in present[1:])
    if not consistent:
        raise PhysicalMomentScaleInputError("source_family labels must all match")
    return first


def _log_interval(log_derivative: Any) -> tuple[RationalInterval, Any]:
    family = _MISSING
    value = log_derivative
    if isinstance(value, Mapping):
        family = _lookup(value, "source_family", "source_id", "family")
        nested = _lookup(value, "finite", "interval", "bounds", "L", "l", "value")
        if nested is not _MISSING:
            value = nested
    return _interval(value, "L"), family


def _normalise_moments(moments: Any) -> dict[str, Any]:
    if not isinstance(moments, Mapping):
        raise PhysicalMomentScaleInputError("moments must be a mapping")
    rows: dict[str, Any] = {}
    for raw_name, row in moments.items():
        name = _moment_name(raw_name)
        if name is not None:
            rows[name] = row
    missing = [name for name in MOMENTS if name not in rows]
    if missing:
        raise PhysicalMomentScaleInputError("missing moment(s): " + ", ".join(missing))
    return rows


def scale(
    moments: Any,
    log_derivative: Any,
    *,
    source_family: Any = None,
    source_families: Any = None,
    errors: Any = None,
) -> dict[str, Any]:
    """Return coefficient rows for the six physical moment scales.

    ``moments[name]`` contains ``z0`` and ``z1`` rows.  Each zone row may be
    ``{"finite": [lower, upper], "error_upper": error}``; ``Z0``/``Z1`` and
    ``finite_Z0``/``finite_Z1`` spellings are accepted as well.  The second
    argument is the exact interval for ``L=F0_Z/F0`` (or a mapping containing
    that interval under ``finite``/``L``).  ``epsilon`` and ``F0`` are not
    materialized: their powers are carried in ``scale`` as exact integer
    pairs.  ``finite``, ``error_upper``, and ``enclosure`` are coefficient
    space; ``total`` is an independent alias of the widened coefficient
    enclosure.
    """

    rows = _normalise_moments(moments)
    log_interval, log_family = _log_interval(log_derivative)
    sup_abs_L = max(abs(log_interval.lower), abs(log_interval.upper))
    family_values: list[Any] = [source_family, log_family]

    z0: dict[str, RationalInterval] = {}
    z1: dict[str, RationalInterval] = {}
    e0: dict[str, Fraction] = {}
    e1: dict[str, Fraction] = {}
    for moment in MOMENTS:
        row = rows[moment]
        raw_z0 = _zone_interval(row, moment, "Z0")
        raw_z1 = _zone_interval(row, moment, "Z1")
        if raw_z0 is _MISSING or raw_z1 is _MISSING:
            raise PhysicalMomentScaleInputError(
                f"{moment} requires finite ordinary Z0 and Z1 intervals"
            )
        z0[moment] = raw_z0
        z1[moment] = raw_z1

        row_e0 = _zone_error(row, moment, "Z0")
        row_e1 = _zone_error(row, moment, "Z1")
        external_e0 = _map_value(errors, moment, "Z0")
        external_e1 = _map_value(errors, moment, "Z1")
        if external_e0 is not _MISSING:
            extra = _error(external_e0, f"{moment}.Z0.external_error")
            row_e0 = extra if row_e0 is _MISSING else row_e0 + extra
        if external_e1 is not _MISSING:
            extra = _error(external_e1, f"{moment}.Z1.external_error")
            row_e1 = extra if row_e1 is _MISSING else row_e1 + extra
        if row_e0 is _MISSING:
            raise PhysicalMomentScaleInputError(
                f"{moment}.Z0 requires an explicit error_upper (zero must be supplied explicitly)"
            )
        if row_e1 is _MISSING:
            raise PhysicalMomentScaleInputError(
                f"{moment}.Z1 requires an explicit error_upper (zero must be supplied explicitly)"
            )
        e0[moment] = row_e0
        e1[moment] = row_e1

        family_values.extend(
            [
                _row_family(row),
                _zone_family(row, "Z0"),
                _zone_family(row, "Z1"),
                _map_value(source_families, moment, "Z0"),
            ]
        )
    family = _validate_families(family_values)

    finite: dict[str, Any] = {}
    finite_z1: dict[str, RationalInterval] = {}
    error_upper: dict[str, Any] = {}
    error_upper_z1: dict[str, Fraction] = {}
    enclosure: dict[str, Any] = {}
    enclosure_z1: dict[str, RationalInterval] = {}
    scale_powers: dict[str, tuple[int, int]] = {}
    records: dict[str, Any] = {}

    for physical_name in PHYSICAL_MOMENTS:
        moment, numerical, powers = _ROWS[physical_name]
        coefficient = Fraction(numerical)
        base = _iv_scale(z0[moment], coefficient)
        derivative_inside = _iv_add(
            z1[moment],
            _iv_scale(_iv_mul(log_interval, z0[moment]), Fraction(powers[1])),
        )
        derivative = _iv_scale(derivative_inside, coefficient)
        base_error = coefficient * e0[moment]
        derivative_error = coefficient * (
            e1[moment] + Fraction(powers[1]) * sup_abs_L * e0[moment]
        )
        widened = RationalInterval(base.lower - base_error, base.upper + base_error)
        widened_z1 = RationalInterval(
            derivative.lower - derivative_error,
            derivative.upper + derivative_error,
        )
        finite[physical_name] = base
        finite_z1[physical_name] = derivative
        error_upper[physical_name] = base_error
        error_upper_z1[physical_name] = derivative_error
        enclosure[physical_name] = widened
        enclosure_z1[physical_name] = widened_z1
        scale_powers[physical_name] = powers
        records[physical_name] = {
            "ordinary_moment": moment,
            "numerical_coefficient": numerical,
            "epsilon_power": powers[0],
            "F0_power": powers[1],
            "scale": powers,
            "finite_coefficient": base,
            "error_upper_coefficient": base_error,
            "enclosure_coefficient": widened,
            "finite_coefficient_Z1": derivative,
            "error_upper_coefficient_Z1": derivative_error,
            "enclosure_coefficient_Z1": widened_z1,
        }

    finite["Z1"] = finite_z1
    error_upper["Z1"] = error_upper_z1
    enclosure["Z1"] = enclosure_z1
    for name in PHYSICAL_MOMENTS:
        finite[f"{name}_Z1"] = finite_z1[name]
        error_upper[f"{name}_Z1"] = error_upper_z1[name]
        enclosure[f"{name}_Z1"] = enclosure_z1[name]

    # ``total`` deliberately remains coefficient-space.  Applying epsilon or
    # F0 requires caller-owned positive values and is outside this leaf.
    total = {name: enclosure[name] for name in PHYSICAL_MOMENTS}
    total_z1 = dict(enclosure_z1)
    total["Z1"] = total_z1
    for name in PHYSICAL_MOMENTS:
        total[f"{name}_Z1"] = total_z1[name]

    return {
        "schema": SCHEMA,
        "L": log_interval,
        "sup_abs_L": sup_abs_L,
        "source_family": family,
        "source_family_consistent": True,
        "finite": finite,
        "error_upper": error_upper,
        "enclosure": enclosure,
        "total": total,
        "finite_Z1": finite_z1,
        "error_upper_Z1": error_upper_z1,
        "enclosure_Z1": enclosure_z1,
        "scale": scale_powers,
        "scale_powers": scale_powers,
        "records": records,
        "assumptions": {
            "epsilon_power_and_F0_power_only": True,
            "epsilon_Z_constant_is_caller_premise": True,
            "F0_positive_is_caller_premise": True,
            "L_is_caller_supplied_interval": True,
            "input_finite_rows_are_caller_enclosures": True,
            "input_errors_are_independent_caps": True,
            "old_errors_cancelled": False,
            "source_family_is_caller_premise": True,
            "source_identity_certified": False,
            "q4_certified": False,
            "closure_proved": False,
            "time_recursion_proved": False,
            "midpoint_selected": False,
            "physical_scale_materialized": False,
        },
    }


evaluate = scale
scale_physical_moments = scale


__all__ = [
    "SCHEMA",
    "MOMENTS",
    "PHYSICAL_MOMENTS",
    "PhysicalMomentScaleInputError",
    "RationalInterval",
    "scale",
    "evaluate",
    "scale_physical_moments",
]


