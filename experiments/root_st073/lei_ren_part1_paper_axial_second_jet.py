"""Second-order axial automatic differentiation over ``PressureWidthJet``.

``AxialSecondJet`` stores a value, its actual first derivative with respect to
the axial coordinate ``Z``, and its actual second derivative with respect to
``Z``.  Every slot is a rectangular :class:`PressureWidthJet`; the pressure
tail and exit-width atoms therefore remain separate while the axial chain
rules are propagated through the three slots.

The second slot is explicit.  A first-order :class:`AxialDual` is not silently
promoted to this ring because it has no information from which a second
derivative can be recovered.  A caller may promote one only by supplying the
missing ``second`` argument explicitly.
"""

from __future__ import annotations

import json
import operator
from typing import Any

import mpmath as mp

try:  # package import
    from .lei_ren_part1_paper_axial_dual import AxialDual
    from .lei_ren_part1_paper_pressure_width_jet import PressureWidthJet
except (ImportError, ValueError):  # direct experiment execution
    from lei_ren_part1_paper_axial_dual import AxialDual
    from lei_ren_part1_paper_pressure_width_jet import PressureWidthJet


DEFAULT_PRESSURE_ORDER = 9
DEFAULT_WIDTH_ORDER = 2
AUTOMATIC_FIRST_Z_DERIVATIVE = True
AUTOMATIC_SECOND_Z_DERIVATIVE = True

_MISSING = object()


def _integer_exponent(exponent: Any) -> int:
    """Convert an exponent accepted by the finite Taylor rings to ``int``."""

    try:
        result = operator.index(exponent)
    except TypeError:
        try:
            candidate = mp.mpf(exponent)
        except (TypeError, ValueError) as exc:
            raise ValueError("AxialSecondJet exponent must be an integer") from exc
        if not mp.isint(candidate):
            raise ValueError("AxialSecondJet exponent must be an integer")
        result = int(candidate)
    return int(result)


def _validate_order(value: Any, name: str) -> int:
    try:
        result = operator.index(value)
    except TypeError as exc:
        raise TypeError(f"{name} must be an integer") from exc
    if result < 0:
        raise ValueError(f"{name} must be nonnegative")
    return int(result)


def _jet_orders(value: Any) -> tuple[int, int] | None:
    if isinstance(value, PressureWidthJet):
        return value.orders
    if isinstance(value, AxialSecondJet):
        return value.orders
    return None


class AxialSecondJet:
    """A value, first ``Z`` derivative, and second ``Z`` derivative.

    The ``tangent`` and ``second`` slots contain actual derivatives, rather
    than Taylor coefficients divided by factorials.  In particular, a
    quadratic scalar ``f(Z)`` is represented by ``second=f''(Z)``.

    Scalar-only construction uses the same effective pressure/width defaults
    as :class:`AxialDual`.  If any ``PressureWidthJet`` slot is supplied, its
    rectangular orders are inferred and all slots must match them.
    """

    automatic_first_Z_derivative = AUTOMATIC_FIRST_Z_DERIVATIVE
    automatic_second_Z_derivative = AUTOMATIC_SECOND_Z_DERIVATIVE
    actual_derivative_slots = True

    def __init__(
        self,
        value: Any = 0,
        tangent: Any = _MISSING,
        second: Any = _MISSING,
        *,
        pressure_order: int | None = None,
        width_order: int | None = None,
    ):
        # Reusing an already complete second jet is lossless.  Explicit slot
        # arguments replace only the slots supplied by the caller.
        if isinstance(value, AxialSecondJet):
            source = value
            value = source.value
            if tangent is _MISSING:
                tangent = source.tangent
            if second is _MISSING:
                second = source.second

        # A first-order dual may be used as the value only when its unknown
        # second derivative is provided explicitly.  Its value and tangent
        # slots remain authoritative unless the caller replaces tangent.
        elif isinstance(value, AxialDual):
            if second is _MISSING:
                raise TypeError(
                    "AxialDual cannot be promoted to AxialSecondJet without "
                    "an explicit second derivative"
                )
            source = value
            value = source.value
            if tangent is _MISSING:
                tangent = source.tangent

        # A first-order dual in either derivative slot is also an unknown
        # second-order quantity.  Do not unwrap it or invent a curvature.
        if isinstance(tangent, AxialDual) or isinstance(second, AxialDual):
            raise TypeError(
                "AxialDual cannot supply an AxialSecondJet derivative slot; "
                "provide a PressureWidthJet or scalar second derivative"
            )
        if isinstance(tangent, AxialSecondJet) or isinstance(second, AxialSecondJet):
            raise TypeError("AxialSecondJet derivative slots must be PressureWidthJet values")

        if tangent is _MISSING:
            tangent = 0
        if second is _MISSING:
            second = 0

        value_orders = _jet_orders(value)
        tangent_orders = _jet_orders(tangent)
        second_orders = _jet_orders(second)
        inferred_orders = value_orders or tangent_orders or second_orders
        supplied_orders = tuple(
            orders
            for orders in (value_orders, tangent_orders, second_orders)
            if orders is not None
        )
        if any(orders != supplied_orders[0] for orders in supplied_orders[1:]):
            raise ValueError(
                "AxialSecondJet value, tangent, and second orders must match: "
                + ", ".join(str(orders) for orders in supplied_orders)
            )

        if pressure_order is None:
            pressure = (
                inferred_orders[0]
                if inferred_orders is not None
                else DEFAULT_PRESSURE_ORDER
            )
        else:
            pressure = _validate_order(pressure_order, "pressure_order")
        if width_order is None:
            width = (
                inferred_orders[1]
                if inferred_orders is not None
                else DEFAULT_WIDTH_ORDER
            )
        else:
            width = _validate_order(width_order, "width_order")

        orders = dict(pressure_order=pressure, width_order=width)
        self.pressure_order = pressure
        self.width_order = width
        self.value = (
            value if isinstance(value, PressureWidthJet) else PressureWidthJet(value, **orders)
        )
        self.tangent = (
            tangent
            if isinstance(tangent, PressureWidthJet)
            else PressureWidthJet(tangent, **orders)
        )
        self.second = (
            second
            if isinstance(second, PressureWidthJet)
            else PressureWidthJet(second, **orders)
        )
        if any(slot.orders != (pressure, width) for slot in (self.value, self.tangent, self.second)):
            raise ValueError("AxialSecondJet slots must use the declared orders")

    @classmethod
    def _from_slots(
        cls,
        value: PressureWidthJet,
        tangent: PressureWidthJet,
        second: PressureWidthJet,
    ) -> "AxialSecondJet":
        """Build a result without coercing or evaluating any slot."""

        if not all(isinstance(slot, PressureWidthJet) for slot in (value, tangent, second)):
            raise TypeError("AxialSecondJet slots must be PressureWidthJet instances")
        if not (value.orders == tangent.orders == second.orders):
            raise ValueError("AxialSecondJet slots must have matching orders")
        result = cls.__new__(cls)
        result.pressure_order, result.width_order = value.orders
        result.value = value
        result.tangent = tangent
        result.second = second
        return result

    @property
    def orders(self) -> tuple[int, int]:
        return self.pressure_order, self.width_order

    @property
    def constant(self) -> mp.mpf:
        """Return the value slot's constant pressure/width coefficient."""

        return self.value.constant

    @property
    def derivative(self) -> PressureWidthJet:
        """Return the actual first axial derivative slot."""

        return self.tangent

    @property
    def second_derivative(self) -> PressureWidthJet:
        """Return the actual second axial derivative slot."""

        return self.second

    @property
    def curvature(self) -> PressureWidthJet:
        """Descriptive alias for :attr:`second_derivative`."""

        return self.second

    def component(self, pressure_power: Any, width_power: Any = 0) -> mp.mpf:
        """Return a value coefficient in the pressure/width jet."""

        return self.value.component(pressure_power, width_power)

    def tangent_component(self, pressure_power: Any, width_power: Any = 0) -> mp.mpf:
        """Return a first-``Z`` derivative coefficient."""

        return self.tangent.component(pressure_power, width_power)

    def second_component(self, pressure_power: Any, width_power: Any = 0) -> mp.mpf:
        """Return a second-``Z`` derivative coefficient."""

        return self.second.component(pressure_power, width_power)

    curvature_component = second_component

    def evaluate(self, pressure: Any = 1, width: Any = 1) -> "AxialSecondJet":
        """Evaluate all three pressure/width slots at one resolved scale."""

        orders = self._order_kwargs()
        return self._from_slots(
            PressureWidthJet(self.value.evaluate(pressure, width), **orders),
            PressureWidthJet(self.tangent.evaluate(pressure, width), **orders),
            PressureWidthJet(self.second.evaluate(pressure, width), **orders),
        )

    __call__ = evaluate

    def as_dict(self) -> dict[str, Any]:
        return {
            "pressure_order": self.pressure_order,
            "width_order": self.width_order,
            "value": self.value.as_dict(),
            "tangent": self.tangent.as_dict(),
            "second": self.second.as_dict(),
            "automatic_first_Z_derivative": True,
            "automatic_second_Z_derivative": True,
            "actual_derivative_slots": True,
        }

    def __repr__(self) -> str:
        return (
            f"AxialSecondJet(value={self.value!r}, tangent={self.tangent!r}, "
            f"second={self.second!r})"
        )

    def _order_kwargs(self) -> dict[str, int]:
        return {
            "pressure_order": self.pressure_order,
            "width_order": self.width_order,
        }

    def _coerce_operand(self, other: Any) -> "AxialSecondJet":
        if isinstance(other, AxialSecondJet):
            if other.orders != self.orders:
                raise ValueError(
                    f"AxialSecondJet orders must match: {self.orders} != {other.orders}"
                )
            return other
        if isinstance(other, AxialDual):
            raise TypeError(
                "AxialDual cannot be promoted to AxialSecondJet without an "
                "explicit second derivative"
            )
        if isinstance(other, PressureWidthJet):
            if other.orders != self.orders:
                raise ValueError(
                    "AxialSecondJet base orders must match: "
                    f"{self.orders} != {other.orders}"
                )
            zero = PressureWidthJet(0, **self._order_kwargs())
            return self._from_slots(other, zero, zero)
        return AxialSecondJet(other, **self._order_kwargs())

    def __add__(self, other: Any) -> "AxialSecondJet":
        right = self._coerce_operand(other)
        return self._from_slots(
            self.value + right.value,
            self.tangent + right.tangent,
            self.second + right.second,
        )

    __radd__ = __add__

    def __neg__(self) -> "AxialSecondJet":
        return self._from_slots(-self.value, -self.tangent, -self.second)

    def __sub__(self, other: Any) -> "AxialSecondJet":
        return self + (-self._coerce_operand(other))

    def __rsub__(self, other: Any) -> "AxialSecondJet":
        return self._coerce_operand(other) - self

    def __mul__(self, other: Any) -> "AxialSecondJet":
        right = self._coerce_operand(other)
        return self._from_slots(
            self.value * right.value,
            self.tangent * right.value + self.value * right.tangent,
            self.second * right.value
            + 2 * self.tangent * right.tangent
            + self.value * right.second,
        )

    __rmul__ = __mul__

    def reciprocal(self) -> "AxialSecondJet":
        """Return ``1 / self`` using the first and second chain rules."""

        value_reciprocal = self.value.reciprocal()
        value_reciprocal_squared = value_reciprocal**2
        value_reciprocal_cubed = value_reciprocal**3
        return self._from_slots(
            value_reciprocal,
            -self.tangent * value_reciprocal_squared,
            2 * self.tangent**2 * value_reciprocal_cubed
            - self.second * value_reciprocal_squared,
        )

    inverse = reciprocal
    inv = reciprocal

    def __truediv__(self, other: Any) -> "AxialSecondJet":
        return self * self._coerce_operand(other).reciprocal()

    def __rtruediv__(self, other: Any) -> "AxialSecondJet":
        return self._coerce_operand(other) * self.reciprocal()

    def __pow__(self, exponent: Any) -> "AxialSecondJet":
        integer = _integer_exponent(exponent)
        if integer == 0:
            return self._from_slots(
                PressureWidthJet(1, **self._order_kwargs()),
                PressureWidthJet(0, **self._order_kwargs()),
                PressureWidthJet(0, **self._order_kwargs()),
            )
        if integer < 0:
            return self.reciprocal() ** (-integer)
        if integer == 1:
            # Avoid evaluating an unnecessary reciprocal in the general
            # second-chain formula; a first power is valid even for a zero
            # constant jet.
            return self._from_slots(self.value, self.tangent, self.second)

        value_power = self.value**integer
        tangent_power = (
            integer * self.value ** (integer - 1) * self.tangent
        )
        second_power = (
            integer * (integer - 1) * self.value ** (integer - 2) * self.tangent**2
            + integer * self.value ** (integer - 1) * self.second
        )
        return self._from_slots(value_power, tangent_power, second_power)

    def exp(self) -> "AxialSecondJet":
        value_exp = self.value.exp()
        return self._from_slots(
            value_exp,
            value_exp * self.tangent,
            value_exp * (self.second + self.tangent**2),
        )

    def log(self) -> "AxialSecondJet":
        value_log = self.value.log()
        return self._from_slots(
            value_log,
            self.tangent / self.value,
            self.second / self.value - (self.tangent / self.value) ** 2,
        )

    def sqrt(self) -> "AxialSecondJet":
        value_sqrt = self.value.sqrt()
        return self._from_slots(
            value_sqrt,
            self.tangent / (2 * value_sqrt),
            self.second / (2 * value_sqrt)
            - self.tangent**2 / (4 * value_sqrt**3),
        )


def main() -> None:
    print(
        json.dumps(
            {
                "class": "AxialSecondJet",
                "pressure_order": DEFAULT_PRESSURE_ORDER,
                "width_order": DEFAULT_WIDTH_ORDER,
                "actual_derivative_slots": True,
            },
            indent=2,
            sort_keys=True,
        )
    )


__all__ = ["AxialSecondJet"]


if __name__ == "__main__":
    main()
