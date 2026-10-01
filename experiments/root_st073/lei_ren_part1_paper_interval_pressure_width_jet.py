"""Directed interval rectangular pressure/width atom arithmetic.

The legacy :class:`PressureWidthJet` stores ordinary ``mp.mpf``
coefficients.  This companion ring keeps the same sparse normalized atom
layout while storing every coefficient as an :class:`IntervalDifference` in
a caller-supplied ``MPIntervalContext``.  The pressure and width atoms remain
separate throughout all operations; evaluating at ``(1, 1)`` is only an
explicit output operation and is never used internally.

The nonlinear functions are finite nilpotent Taylor series in the retained
rectangle.  They enclose the nominal and common perturbation components, but
do not add a remainder for atoms outside the declared rectangle.
"""

from __future__ import annotations

from collections.abc import Mapping
import operator
from typing import Any

import mpmath as mp

from lei_ren_part1_paper_interval_difference import IntervalDifference


DEFAULT_PRESSURE_ORDER = 3
DEFAULT_WIDTH_ORDER = 2
PRESSURE_WIDTH_RECTANGULAR_TRUNCATED = True
PRESSURE_PARAMETER_ORDER_TRUNCATED = True
EXIT_WIDTH_ORDER_TRUNCATED = True
NORMALIZED_BIVARIATE_CONVENTION = (
    "coefficient of pressure**p * width**w (derivatives divided by p! * w!)"
)


def _endpoints(value: Any) -> tuple[mp.mpf, mp.mpf]:
    """Return exact MP endpoints for an interval value."""

    return tuple(mp.make_mpf(item) for item in value._mpi_)


def _as_interval(ctx: Any, value: Any) -> Any:
    """Convert a scalar or interval into the supplied interval context."""

    if hasattr(value, "_mpi_"):
        # Rebuild through the supplied context so operations never silently
        # use an interval object belonging to another context.
        lower, upper = _endpoints(value)
        return ctx.mpf([lower, upper])
    if isinstance(value, float):
        # Avoid importing binary-float roundoff into exact stored inputs.
        return ctx.mpf(str(value))
    return ctx.mpf(value)


def _validate_order(value: Any, name: str) -> int:
    try:
        order = operator.index(value)
    except TypeError as exc:
        raise TypeError(f"{name} must be an integer") from exc
    if order < 0:
        raise ValueError(f"{name} must be nonnegative")
    return int(order)


def _index(value: Any, name: str) -> int:
    try:
        result = operator.index(value)
    except TypeError:
        try:
            candidate = mp.mpf(value)
        except (TypeError, ValueError) as exc:
            raise TypeError(f"{name} must be an integer") from exc
        if not mp.isint(candidate):
            raise TypeError(f"{name} must be an integer")
        result = int(candidate)
    return int(result)


def _key(raw_key: Any) -> tuple[int, int]:
    """Normalize an atom key; an integer key means ``(pressure, 0)``."""

    if isinstance(raw_key, tuple):
        if len(raw_key) != 2:
            raise TypeError("Bivariate Taylor atom keys must be (pressure, width)")
        pressure = _index(raw_key[0], "pressure power")
        width = _index(raw_key[1], "width power")
    else:
        pressure = _index(raw_key, "pressure power")
        width = 0
    if pressure < 0 or width < 0:
        raise ValueError("Bivariate Taylor powers must be nonnegative")
    return pressure, width


def _exact_zero(value: Any) -> bool:
    lower, upper = _endpoints(value)
    return lower == 0 and upper == 0


def _pair_zero(ctx: Any) -> IntervalDifference:
    zero = ctx.mpf(0)
    return IntervalDifference(ctx, zero, zero)


def _pair_one(ctx: Any) -> IntervalDifference:
    return IntervalDifference(ctx, ctx.mpf(1), ctx.mpf(0))


def _pair_rational(ctx: Any, numerator: int, denominator: int) -> IntervalDifference:
    """Build a directed pair for an exact rational coefficient."""

    if denominator == 0:
        raise ZeroDivisionError("rational denominator is zero")
    value = ctx.mpf(numerator) / ctx.mpf(denominator)
    return IntervalDifference(ctx, value, ctx.mpf(0))


def _pair(ctx: Any, value: Any) -> IntervalDifference:
    if isinstance(value, IntervalDifference):
        if value.ctx is not ctx:
            raise ValueError("IntervalDifference contexts must match")
        return value
    return IntervalDifference(ctx, _as_interval(ctx, value), ctx.mpf(0))


class IntervalPressureWidthJet:
    """Sparse rectangular pressure/width jet with directed pair coefficients."""

    rectangular_truncation = PRESSURE_WIDTH_RECTANGULAR_TRUNCATED
    normalized_coefficient_convention = NORMALIZED_BIVARIATE_CONVENTION
    pressure_parameter_order_truncated = PRESSURE_PARAMETER_ORDER_TRUNCATED
    exit_width_order_truncated = EXIT_WIDTH_ORDER_TRUNCATED

    def __init__(
        self,
        ctx: Any,
        value: Any = 0,
        *,
        pressure_order: int = DEFAULT_PRESSURE_ORDER,
        width_order: int = DEFAULT_WIDTH_ORDER,
    ) -> None:
        self.ctx = ctx
        self.pressure_order = _validate_order(pressure_order, "pressure_order")
        self.width_order = _validate_order(width_order, "width_order")
        self._truncated_input_atoms: tuple[tuple[int, int], ...] = ()

        if isinstance(value, IntervalPressureWidthJet):
            if value.ctx is not ctx:
                raise ValueError("Interval pressure/width contexts must match")
            if value.orders != (self.pressure_order, self.width_order):
                raise ValueError(
                    "Bivariate Taylor orders do not match: "
                    f"{(self.pressure_order, self.width_order)} != {value.orders}"
                )
            source: Any = value.atoms
        elif hasattr(value, "atoms"):
            # This accepts the legacy PressureWidthJet and PressurePolynomial
            # atom maps without importing either implementation.
            source_orders = getattr(value, "orders", None)
            if source_orders is not None and tuple(source_orders) != (
                self.pressure_order,
                self.width_order,
            ):
                raise ValueError(
                    "Bivariate Taylor orders do not match: "
                    f"{(self.pressure_order, self.width_order)} != {source_orders}"
                )
            source = value.atoms
        elif isinstance(value, Mapping):
            source = value
        else:
            source = {(0, 0): value}

        if not isinstance(source, Mapping):
            raise TypeError("Taylor atoms must be supplied by a mapping")

        atoms: dict[tuple[int, int], IntervalDifference] = {}
        discarded: list[tuple[int, int]] = []
        for raw_key, raw_coefficient in source.items():
            key = _key(raw_key)
            pressure, width = key
            if pressure > self.pressure_order or width > self.width_order:
                discarded.append(key)
                continue
            coefficient = _pair(ctx, raw_coefficient)
            # A coefficient is pruned only when both its nominal and its
            # perturbation are exact zero intervals.
            if not (_exact_zero(coefficient.nominal) and _exact_zero(coefficient.difference)):
                atoms[key] = coefficient
        self.atoms = atoms
        self._truncated_input_atoms = tuple(sorted(set(discarded)))

    @classmethod
    def _from_atoms(
        cls,
        ctx: Any,
        atoms: Mapping[tuple[int, int], Any],
        pressure_order: int,
        width_order: int,
    ) -> "IntervalPressureWidthJet":
        result = cls.__new__(cls)
        result.ctx = ctx
        result.pressure_order = pressure_order
        result.width_order = width_order
        result._truncated_input_atoms = ()
        result.atoms = {}
        for key, raw_value in atoms.items():
            coefficient = _pair(ctx, raw_value)
            if not (_exact_zero(coefficient.nominal) and _exact_zero(coefficient.difference)):
                result.atoms[key] = coefficient
        return result

    @property
    def orders(self) -> tuple[int, int]:
        return self.pressure_order, self.width_order

    @classmethod
    def converter(
        cls,
        ctx: Any,
        *,
        pressure_order: int = DEFAULT_PRESSURE_ORDER,
        width_order: int = DEFAULT_WIDTH_ORDER,
    ):
        """Return a context/order-preserving one-argument constructor."""

        def convert(value: Any) -> "IntervalPressureWidthJet":
            if isinstance(value, cls):
                if value.ctx is not ctx:
                    raise ValueError("Interval pressure/width contexts must match")
                if value.orders != (pressure_order, width_order):
                    raise ValueError("Bivariate Taylor orders do not match converter")
                return value
            return cls(
                ctx,
                value,
                pressure_order=pressure_order,
                width_order=width_order,
            )

        return convert

    factory = converter
    scalar_converter = converter

    @property
    def constant(self) -> IntervalDifference:
        return self.component(0, 0)

    @property
    def truncated_input_atoms(self) -> tuple[tuple[int, int], ...]:
        return self._truncated_input_atoms

    def component(self, pressure_power: Any, width_power: Any = 0) -> IntervalDifference:
        pressure = _index(pressure_power, "pressure component index")
        width = _index(width_power, "width component index")
        if pressure < 0 or width < 0:
            return _pair_zero(self.ctx)
        if pressure > self.pressure_order or width > self.width_order:
            return _pair_zero(self.ctx)
        return self.atoms.get((pressure, width), _pair_zero(self.ctx))

    coefficient = component

    def evaluate(self, pressure: Any = 1, width: Any = 1) -> IntervalDifference:
        """Evaluate retained atoms at explicit scalar/pair parameters."""

        pressure_value = _pair(self.ctx, pressure)
        width_value = _pair(self.ctx, width)
        result = _pair_zero(self.ctx)
        for (pressure_power, width_power), coefficient in self.atoms.items():
            result = result + coefficient * pressure_value**pressure_power * width_value**width_power
        return result

    __call__ = evaluate

    def as_dict(self) -> dict[str, Any]:
        def describe(value: IntervalDifference) -> dict[str, str]:
            nlo, nhi = _endpoints(value.nominal)
            dlo, dhi = _endpoints(value.difference)
            return {
                "nominal_lower": mp.nstr(nlo, 30),
                "nominal_upper": mp.nstr(nhi, 30),
                "difference_lower": mp.nstr(dlo, 30),
                "difference_upper": mp.nstr(dhi, 30),
            }

        return {
            "pressure_order": self.pressure_order,
            "width_order": self.width_order,
            "atoms": {
                f"{pressure},{width}": describe(coefficient)
                for (pressure, width), coefficient in sorted(self.atoms.items())
            },
            "rectangular_truncation": True,
            "normalized_coefficient_convention": NORMALIZED_BIVARIATE_CONVENTION,
            "truncated_input_atoms": [list(key) for key in self.truncated_input_atoms],
        }

    def __repr__(self) -> str:
        return (
            f"IntervalPressureWidthJet({self.as_dict()['atoms']!r}, "
            f"pressure_order={self.pressure_order}, width_order={self.width_order})"
        )

    def _coerce_operand(self, other: Any) -> "IntervalPressureWidthJet":
        if isinstance(other, IntervalPressureWidthJet):
            if other.ctx is not self.ctx:
                raise ValueError("Interval pressure/width contexts must match")
            if other.orders != self.orders:
                raise ValueError(f"Bivariate Taylor orders do not match: {self.orders} != {other.orders}")
            return other
        return IntervalPressureWidthJet(
            self.ctx,
            other,
            pressure_order=self.pressure_order,
            width_order=self.width_order,
        )

    def _scale(self, scalar: Any) -> "IntervalPressureWidthJet":
        coefficient = _pair(self.ctx, scalar)
        return self._from_atoms(
            self.ctx,
            {key: value * coefficient for key, value in self.atoms.items()},
            self.pressure_order,
            self.width_order,
        )

    def _nonconstant(self, factor: Any = None) -> "IntervalPressureWidthJet":
        """Return the structural nilpotent part, omitting ``(0, 0)``.

        Directed interval subtraction cannot reliably cancel a constant
        interval against itself.  Nonlinear finite-series operations must
        remove that atom by key before multiplying powers.
        """

        if factor is None:
            coefficient_factor = None
        else:
            coefficient_factor = _pair(self.ctx, factor)
        atoms = {
            key: value if coefficient_factor is None else value * coefficient_factor
            for key, value in self.atoms.items()
            if key != (0, 0)
        }
        return self._from_atoms(
            self.ctx,
            atoms,
            self.pressure_order,
            self.width_order,
        )

    def __add__(self, other: Any) -> "IntervalPressureWidthJet":
        right = self._coerce_operand(other)
        terms: dict[tuple[int, int], list[IntervalDifference]] = {}
        for key, value in self.atoms.items():
            terms.setdefault(key, []).append(value)
        for key, value in right.atoms.items():
            terms.setdefault(key, []).append(value)
        return self._from_atoms(
            self.ctx,
            {key: sum(values[1:], values[0]) for key, values in terms.items()},
            self.pressure_order,
            self.width_order,
        )

    __radd__ = __add__

    def __neg__(self) -> "IntervalPressureWidthJet":
        return self._from_atoms(
            self.ctx,
            {key: -value for key, value in self.atoms.items()},
            self.pressure_order,
            self.width_order,
        )

    def __sub__(self, other: Any) -> "IntervalPressureWidthJet":
        return self + (-self._coerce_operand(other))

    def __rsub__(self, other: Any) -> "IntervalPressureWidthJet":
        return self._coerce_operand(other) - self

    def __mul__(self, other: Any) -> "IntervalPressureWidthJet":
        right = self._coerce_operand(other)
        terms: dict[tuple[int, int], list[IntervalDifference]] = {}
        for (pressure_a, width_a), value_a in self.atoms.items():
            for (pressure_b, width_b), value_b in right.atoms.items():
                key = pressure_a + pressure_b, width_a + width_b
                if key[0] <= self.pressure_order and key[1] <= self.width_order:
                    terms.setdefault(key, []).append(value_a * value_b)
        return self._from_atoms(
            self.ctx,
            {key: sum(values[1:], values[0]) for key, values in terms.items()},
            self.pressure_order,
            self.width_order,
        )

    __rmul__ = __mul__

    def reciprocal(self) -> "IntervalPressureWidthJet":
        """Return the rectangular finite-series reciprocal."""

        constant = self.constant
        # IntervalDifference.reciprocal performs both nominal and total
        # denominator checks before any nilpotent expansion is attempted.
        constant_inverse = constant.reciprocal()
        nonconstant = self._nonconstant(constant_inverse)
        result = self._from_atoms(self.ctx, {(0, 0): _pair_one(self.ctx)}, *self.orders)
        term = result
        for power in range(1, self.pressure_order + self.width_order + 1):
            term = term * nonconstant
            if not term.atoms:
                break
            result = result + term._scale(_pair_rational(self.ctx, (-1) ** power, 1))
        return result._scale(constant_inverse)

    inverse = reciprocal
    inv = reciprocal

    def __truediv__(self, other: Any) -> "IntervalPressureWidthJet":
        return self * self._coerce_operand(other).reciprocal()

    def __rtruediv__(self, other: Any) -> "IntervalPressureWidthJet":
        return self._coerce_operand(other) * self.reciprocal()

    def __pow__(self, exponent: Any) -> "IntervalPressureWidthJet":
        try:
            integer = operator.index(exponent)
        except TypeError:
            try:
                candidate = mp.mpf(exponent)
            except (TypeError, ValueError) as exc:
                raise ValueError("Taylor exponent must be an integer") from exc
            if not mp.isint(candidate):
                raise ValueError("Taylor exponent must be an integer")
            integer = int(candidate)
        integer = int(integer)
        if integer < 0:
            return self.reciprocal() ** (-integer)
        result = self._from_atoms(self.ctx, {(0, 0): _pair_one(self.ctx)}, *self.orders)
        base = self
        while integer:
            if integer & 1:
                result = result * base
            integer //= 2
            if integer:
                base = base * base
        return result

    def exp(self) -> "IntervalPressureWidthJet":
        """Return ``exp(self)`` through the retained rectangle."""

        constant = self.constant
        scale = constant.exp()
        nonconstant = self._nonconstant()
        result = self._from_atoms(self.ctx, {(0, 0): scale}, *self.orders)
        term = self._from_atoms(self.ctx, {(0, 0): _pair_one(self.ctx)}, *self.orders)
        factorial = 1
        for power in range(1, self.pressure_order + self.width_order + 1):
            term = term * nonconstant
            if not term.atoms:
                break
            factorial *= power
            result = result + term._scale(
                scale * _pair_rational(self.ctx, 1, factorial)
            )
        return result

    def log(self) -> "IntervalPressureWidthJet":
        """Return the real analytic logarithm with strict positive constant."""

        constant = self.constant
        scale = constant.log()
        nonconstant = self._nonconstant(constant.reciprocal())
        result = self._from_atoms(self.ctx, {(0, 0): scale}, *self.orders)
        term = self._from_atoms(self.ctx, {(0, 0): _pair_one(self.ctx)}, *self.orders)
        for power in range(1, self.pressure_order + self.width_order + 1):
            term = term * nonconstant
            if not term.atoms:
                break
            result = result + term._scale(
                _pair_rational(self.ctx, (-1) ** (power + 1), power)
            )
        return result

    def sqrt(self) -> "IntervalPressureWidthJet":
        """Return the real analytic square root with strict positive constant."""

        constant = self.constant
        root = constant.sqrt()
        nonconstant = self._nonconstant(constant.reciprocal())
        result = self._from_atoms(self.ctx, {(0, 0): root}, *self.orders)
        term = self._from_atoms(self.ctx, {(0, 0): _pair_one(self.ctx)}, *self.orders)
        binomial = _pair_one(self.ctx)
        half = _pair_rational(self.ctx, 1, 2)
        for power in range(1, self.pressure_order + self.width_order + 1):
            term = term * nonconstant
            if not term.atoms:
                break
            binomial = binomial * (half - (power - 1)) / power
            result = result + term._scale(root * binomial)
        return result


# Descriptive aliases for callers that put the interval qualifier first or
# last.  The implementation remains one class so order/context checks stay
# identical.
PressureWidthIntervalJet = IntervalPressureWidthJet
IntervalPressureWidthAtomRing = IntervalPressureWidthJet


__all__ = [
    "IntervalPressureWidthJet",
    "PressureWidthIntervalJet",
    "IntervalPressureWidthAtomRing",
    "DEFAULT_PRESSURE_ORDER",
    "DEFAULT_WIDTH_ORDER",
]
