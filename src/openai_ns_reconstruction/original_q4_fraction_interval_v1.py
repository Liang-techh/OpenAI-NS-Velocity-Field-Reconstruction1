"""Exact rational interval arithmetic for a conditional Q4 integral layer.

This module is deliberately a small numerical boundary.  ``Interval`` keeps
both directed endpoints as :class:`fractions.Fraction` values and never
turns an enclosure, tail bound, or component result into a point value.  The
range integrator only calls an injected function on complete rational cells;
it does not sample a midpoint and it has no derivative-cap or source-provider
knowledge.

The resulting :class:`IntegralEnclosure` is therefore conditional: its
meaning depends on the caller having supplied a callable that returns a valid
enclosure for every complete cell.  No source graph, binder, provider,
rounding procedure, quadrature proof, or physical error is authenticated here.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from fractions import Fraction
from typing import Any, Callable


SCHEMA = "original_q4_fraction_interval_v1"

# A finite limit prevents an accidental request from allocating an
# unbounded coverage-cell list.  It is a resource guard, not an error bound.
MAX_PARTITIONS = 100_000
_MAX_PARTITIONS = MAX_PARTITIONS
_MISSING = object()


def _fraction(value: Any, label: str) -> Fraction:
    """Coerce only an exact integer or ``Fraction``.

    In particular, ``bool`` is rejected even though it is an ``int``
    subclass, and no ``float``/``Decimal``/string conversion is attempted.
    """

    if type(value) is int:
        return Fraction(value, 1)
    if isinstance(value, Fraction):
        return Fraction(value.numerator, value.denominator)
    raise TypeError(f"{label} must be an exact int or Fraction")


@dataclass(frozen=True, slots=True)
class Interval:
    """A closed exact-rational interval tagged as an enclosure.

    The constructor accepts only integer or :class:`Fraction` endpoints.
    Integers are normalized to ``Fraction`` so every subsequent operation is
    exact.  The tag is fixed and cannot be changed to a selected point value.
    """

    lo: Fraction | int
    hi: Fraction | int
    kind: str = field(default="enclosure", init=False)

    def __post_init__(self) -> None:
        lower = _fraction(self.lo, "lo")
        upper = _fraction(self.hi, "hi")
        if lower > upper:
            raise ValueError("interval lower endpoint exceeds upper endpoint")
        object.__setattr__(self, "lo", lower)
        object.__setattr__(self, "hi", upper)

    @property
    def lower(self) -> Fraction:
        """Alias for the lower directed endpoint."""

        return self.lo

    @property
    def upper(self) -> Fraction:
        """Alias for the upper directed endpoint."""

        return self.hi

    @property
    def width(self) -> Fraction:
        """Exact closed-interval width."""

        return self.hi - self.lo

    def _operand(self, value: Any) -> "Interval":
        if isinstance(value, Interval):
            if value.kind != "enclosure":
                # This is defensive for a foreign subclass/property.  A
                # normal Interval can never carry another tag.
                raise TypeError("interval operand is not an enclosure")
            return value
        scalar = _fraction(value, "scalar operand")
        return Interval(scalar, scalar)

    def __add__(self, other: Any) -> "Interval":
        rhs = self._operand(other)
        return Interval(self.lo + rhs.lo, self.hi + rhs.hi)

    def __radd__(self, other: Any) -> "Interval":
        return self + other

    def __sub__(self, other: Any) -> "Interval":
        rhs = self._operand(other)
        return Interval(self.lo - rhs.hi, self.hi - rhs.lo)

    def __rsub__(self, other: Any) -> "Interval":
        lhs = self._operand(other)
        return lhs - self

    def __mul__(self, other: Any) -> "Interval":
        rhs = self._operand(other)
        products = (
            self.lo * rhs.lo,
            self.lo * rhs.hi,
            self.hi * rhs.lo,
            self.hi * rhs.hi,
        )
        return Interval(min(products), max(products))

    def __rmul__(self, other: Any) -> "Interval":
        return self * other

    @staticmethod
    def _positive_denominator(value: "Interval | Fraction | int") -> "Interval":
        if isinstance(value, Interval):
            denominator = value
        else:
            scalar = _fraction(value, "denominator")
            denominator = Interval(scalar, scalar)
        if denominator.lo == 0:
            raise ZeroDivisionError("interval denominator contains zero")
        if denominator.lo < 0:
            raise ValueError("interval denominator must be strictly positive")
        return denominator

    def __truediv__(self, other: Any) -> "Interval":
        denominator = self._positive_denominator(other)
        # The reciprocal reverses the positive endpoint order.  Multiplying
        # by it then performs the full four-corner enclosure arithmetic.
        reciprocal = Interval(Fraction(1, denominator.hi), Fraction(1, denominator.lo))
        return self * reciprocal

    def __rtruediv__(self, other: Any) -> "Interval":
        numerator = self._operand(other)
        return numerator / self

    def __neg__(self) -> "Interval":
        return Interval(-self.hi, -self.lo)

    def __abs__(self) -> "Interval":
        if self.lo >= 0:
            return self
        if self.hi <= 0:
            return Interval(-self.hi, -self.lo)
        return Interval(Fraction(0), max(-self.lo, self.hi))

    def abs(self) -> "Interval":
        """Return the exact enclosure of the pointwise absolute value."""

        return abs(self)

    def is_positive(self) -> bool:
        """Return whether every point in the interval is strictly positive."""

        return self.lo > 0

    @property
    def is_strictly_positive(self) -> bool:
        """Property spelling useful to callers that inspect an interval."""

        return self.is_positive()

    def __float__(self) -> float:
        raise TypeError("an interval enclosure cannot be converted to float")

    def __int__(self) -> int:
        raise TypeError("an interval enclosure cannot be converted to int")

    def __index__(self) -> int:
        raise TypeError("an interval enclosure cannot be used as an integer")

    def __bool__(self) -> bool:
        raise TypeError("an interval enclosure has no scalar truth value")


@dataclass(frozen=True, slots=True)
class IntegralEnclosure:
    """Conditional enclosure returned by :func:`range_integrate`.

    ``coverage_cells`` records the complete rational cells passed to the
    callback.  The fixed flags intentionally make no claim that the callback
    came from an authenticated Q4 provider or that a point value was chosen.
    """

    interval: Interval
    coverage_cells: tuple[Interval, ...] = field(default_factory=tuple)
    kind: str = field(default="conditional_enclosure", init=False)
    provider_certified: bool = field(default=False, init=False)
    actual_value_selected: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        if not isinstance(self.interval, Interval) or self.interval.kind != "enclosure":
            raise TypeError("IntegralEnclosure.interval must be an Interval enclosure")
        cells = tuple(self.coverage_cells)
        if any(not isinstance(cell, Interval) or cell.kind != "enclosure" for cell in cells):
            raise TypeError("coverage_cells must contain only Interval enclosures")
        object.__setattr__(self, "coverage_cells", cells)

    @property
    def coverage_count(self) -> int:
        return len(self.coverage_cells)

    @property
    def count(self) -> int:
        """Short alias for the number of covered cells."""

        return self.coverage_count

    @property
    def cells(self) -> tuple[Interval, ...]:
        return self.coverage_cells

    @property
    def conditional(self) -> bool:
        return True

    def receipt(self) -> dict[str, Any]:
        """Return a non-certifying receipt for callers that need to retain it."""

        return {
            "schema": SCHEMA,
            "kind": self.kind,
            "coverage_count": self.coverage_count,
            "provider_certified": False,
            "actual_value_selected": False,
            "condition": (
                "caller supplied a callable returning a valid enclosure on "
                "each complete rational cell; source and numerical claims "
                "remain unverified"
            ),
        }

    def __float__(self) -> float:
        raise TypeError("an integral enclosure cannot be converted to float")

    def __bool__(self) -> bool:
        raise TypeError("an integral enclosure has no scalar truth value")


class Backend:
    """Backend adapter whose numeric values are always ``Interval`` tags."""

    @staticmethod
    def _value(value: Any) -> Interval:
        if isinstance(value, Interval):
            if value.kind != "enclosure":
                raise TypeError("backend received a non-enclosure interval")
            return value
        scalar = _fraction(value, "backend value")
        return Interval(scalar, scalar)

    def rational(self, numerator: Any, denominator: Any = 1) -> Interval:
        num = _fraction(numerator, "rational numerator")
        den = _fraction(denominator, "rational denominator")
        if den == 0:
            raise ZeroDivisionError("rational denominator is zero")
        if den < 0:
            raise ValueError("rational denominator must be strictly positive")
        value = num / den
        return Interval(value, value)

    def add(self, left: Any, right: Any) -> Interval:
        return self._value(left) + self._value(right)

    def subtract(self, left: Any, right: Any) -> Interval:
        return self._value(left) - self._value(right)

    def multiply(self, left: Any, right: Any) -> Interval:
        return self._value(left) * self._value(right)

    def negative(self, value: Any) -> Interval:
        return -self._value(value)

    def divide(self, numerator: Any, denominator: Any) -> Interval:
        return self._value(numerator) / self._value(denominator)

    def unary(self, name: Any, value: Any) -> Interval:
        if name != "abs":
            raise ValueError("only unary('abs', value) is supported")
        return abs(self._value(value))

    def is_positive(self, value: Any) -> bool:
        result = self._value(value).is_positive()
        # Keep the public result an actual bool even if a foreign subclass
        # supplies an unusual comparison implementation.
        if type(result) is not bool:
            raise TypeError("interval positivity did not return bool")
        return result


FractionIntervalBackend = Backend
ExactIntervalBackend = Backend
IntervalBackend = Backend


def _integration_endpoint(value: Any, label: str) -> Fraction:
    if isinstance(value, Interval):
        if value.kind != "enclosure":
            raise TypeError(f"{label} must be an enclosure interval")
        if value.lo != value.hi:
            raise ValueError(f"{label} Interval must be degenerate")
        return value.lo
    return _fraction(value, label)


def _partition_count(value: Any) -> int:
    if type(value) is not int or value <= 0:
        raise ValueError("partitions must be a positive integer")
    if value > MAX_PARTITIONS:
        raise ValueError(f"partitions exceeds the limit {MAX_PARTITIONS}")
    return value


def range_integrate(
    integrand: Callable[[Interval], Interval],
    lo: Any,
    hi: Any,
    *,
    partitions: int,
) -> IntegralEnclosure:
    """Enclose an injected integrand over a rational closed range.

    Every callback call receives the entire cell ``[a, b]``.  A callback that
    returns a scalar, a midpoint, a tail cap, or another untagged object is
    rejected.  The caller is responsible for the callback's source,
    integrability, and enclosure validity; this function records that
    condition but does not authenticate it.
    """

    if not callable(integrand):
        raise TypeError("integrand must be callable")
    left = _integration_endpoint(lo, "lo")
    right = _integration_endpoint(hi, "hi")
    count = _partition_count(partitions)
    if left > right:
        raise ValueError("integration lower endpoint exceeds upper endpoint")

    width = (right - left) / count
    total = Interval(Fraction(0), Fraction(0))
    cells: list[Interval] = []
    for index in range(count):
        a = left + width * index
        b = left + width * (index + 1)
        cell = Interval(a, b)
        value = integrand(cell)
        if not isinstance(value, Interval) or value.kind != "enclosure":
            raise TypeError("integrand must return an Interval enclosure for every cell")
        cells.append(cell)
        total = total + value * width
    return IntegralEnclosure(total, tuple(cells))


def range_integrate_interval(
    integrand: Callable[[Interval], Interval],
    lo: Any,
    hi: Any,
    *,
    partitions: int,
    receipt_sink: list[IntegralEnclosure],
) -> Interval:
    """Compatibility view returning only the tagged integral interval.

    The receipt sink is mandatory so the conditional coverage record is not
    silently discarded.  No source, provider or quadrature is certified.
    """

    if not isinstance(receipt_sink, list):
        raise TypeError("a list receipt sink is required")
    result = range_integrate(integrand, lo, hi, partitions=partitions)
    receipt_sink.append(result)
    return result.interval


class RangeIntegrator:
    """Oracle-compatible integration with retained conditional receipts."""

    def __init__(self, partitions: int):
        self.partitions = _partition_count(partitions)
        self.receipts: list[IntegralEnclosure] = []

    def integrate(self, integrand: Callable[[Interval], Interval], lo: Any, hi: Any) -> Interval:
        return range_integrate_interval(integrand, lo, hi, partitions=self.partitions,
                                        receipt_sink=self.receipts)


def _endpoint_pair(value: Any) -> tuple[Any, Any]:
    if isinstance(value, Interval):
        return value.lo, value.hi
    if isinstance(value, (tuple, list)) and len(value) == 2:
        return value[0], value[1]
    raise TypeError(
        "to_fraction_endpoints must return an Interval or a two-endpoint tuple"
    )


@dataclass(frozen=True, slots=True)
class ComponentIntervalExtraction:
    """Interval conversion whose source and directed rounding are unverified."""
    interval: Interval
    kind: str = field(default="conditional_enclosure", init=False)
    endpoint_conversion_verified: bool = field(default=False, init=False)
    source_identity_verified: bool = field(default=False, init=False)
    actual_value_selected: bool = field(default=False, init=False)


def extract_component_interval(
    component: Any,
    to_fraction_endpoints: Callable[[Any], Any],
) -> ComponentIntervalExtraction:
    """Extract exact endpoints from a caller-supplied component enclosure.

    The conversion callback is intentionally explicit and remains an
    unverified caller claim.  It must return exact endpoints (or an existing
    ``Interval``); no float conversion, midpoint, provider lookup, or source
    authentication is attempted here.
    """

    if not callable(to_fraction_endpoints):
        raise TypeError("an explicit endpoint conversion callback is required")
    if component is None:
        raise TypeError("component enclosure is required")
    if getattr(component, "kind", _MISSING) != "enclosure":
        raise TypeError("component must be tagged kind='enclosure'")
    if getattr(component, "actual_value_selected", _MISSING) is not False:
        raise ValueError("component must retain enclosure semantics")
    raw = getattr(component, "interval", _MISSING)
    if raw is _MISSING:
        raise TypeError("component has no complete interval member")
    converted = to_fraction_endpoints(raw)
    lower, upper = _endpoint_pair(converted)
    return ComponentIntervalExtraction(Interval(lower, upper))


# Explicit aliases make the caller-directed boundary discoverable without
# coupling this module to the existing adapter's implementation module.
extract_component_enclosure = extract_component_interval
extract_component_enclosure_interval = extract_component_interval
extract_interval = extract_component_interval
extract = extract_component_interval


__all__ = [
    "SCHEMA",
    "MAX_PARTITIONS",
    "Interval",
    "IntegralEnclosure",
    "Backend",
    "FractionIntervalBackend",
    "ExactIntervalBackend",
    "IntervalBackend",
    "range_integrate",
    "range_integrate_interval",
    "RangeIntegrator",
    "ComponentIntervalExtraction",
    "extract_component_interval",
    "extract_component_enclosure",
    "extract_component_enclosure_interval",
    "extract_interval",
    "extract",
]
