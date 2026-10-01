"""Second-order axial jet over the directed pressure/width interval ring.

``IntervalAxialSecondJet`` stores three *actual* derivatives with respect to
the axial variable ``Z``: ``value``, ``tangent`` and ``second``.  Each slot is
an :class:`IntervalPressureWidthJet`, so the pressure and exit-width atoms
remain available independently while ``IntervalDifference`` encloses nominal
coefficients and their perturbations.

The constructor deliberately requires an interval context.  Pressure and
width orders may be inferred only from an interval-ring slot; scalar-only
construction must name both orders.  A first-order ``AxialDual`` cannot be
promoted because it does not contain a second derivative.
"""

from __future__ import annotations

import json
import operator
from typing import Any

import mpmath as mp

try:  # package import
    from .lei_ren_part1_paper_axial_dual import AxialDual
    from .lei_ren_part1_paper_interval_pressure_width_jet import (
        DEFAULT_PRESSURE_ORDER,
        DEFAULT_WIDTH_ORDER,
        IntervalPressureWidthJet,
    )
except (ImportError, ValueError):  # direct experiment execution
    from lei_ren_part1_paper_axial_dual import AxialDual
    from lei_ren_part1_paper_interval_pressure_width_jet import (
        DEFAULT_PRESSURE_ORDER,
        DEFAULT_WIDTH_ORDER,
        IntervalPressureWidthJet,
    )


AUTOMATIC_FIRST_Z_DERIVATIVE = True
AUTOMATIC_SECOND_Z_DERIVATIVE = True
ACTUAL_DERIVATIVE_SLOTS = True
_MISSING = object()


def _integer_exponent(exponent: Any) -> int:
    """Convert an integer exponent accepted by the interval ring."""

    try:
        result = operator.index(exponent)
    except TypeError:
        try:
            candidate = mp.mpf(exponent)
        except (TypeError, ValueError) as exc:
            raise ValueError("IntervalAxialSecondJet exponent must be an integer") from exc
        if not mp.isint(candidate):
            raise ValueError("IntervalAxialSecondJet exponent must be an integer")
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


def _slot_orders(value: Any) -> tuple[int, int] | None:
    if isinstance(value, IntervalPressureWidthJet):
        return value.orders
    return None


class IntervalAxialSecondJet:
    """A value, first ``Z`` derivative, and second ``Z`` derivative.

    The derivative slots are actual derivatives rather than factorial-scaled
    Taylor coefficients.  ``ctx`` is always the first argument.  If any slot
    is an interval pressure/width jet, omitted orders are inferred from that
    ring; scalar-only construction requires explicit ``pressure_order`` and
    ``width_order``.
    """

    automatic_first_Z_derivative = AUTOMATIC_FIRST_Z_DERIVATIVE
    automatic_second_Z_derivative = AUTOMATIC_SECOND_Z_DERIVATIVE
    actual_derivative_slots = ACTUAL_DERIVATIVE_SLOTS

    def __init__(
        self,
        ctx: Any,
        value: Any = 0,
        tangent: Any = _MISSING,
        second: Any = _MISSING,
        *,
        pressure_order: int | None = None,
        width_order: int | None = None,
    ) -> None:
        if ctx is None:
            raise TypeError("IntervalAxialSecondJet requires an interval context")

        # Reusing a complete interval second jet is lossless.  Its context and
        # orders remain authoritative unless the caller supplies matching ones.
        if isinstance(value, IntervalAxialSecondJet):
            if value.ctx is not ctx:
                raise ValueError("Interval axial jet contexts must match")
            source = value
            value = source.value
            if tangent is _MISSING:
                tangent = source.tangent
            if second is _MISSING:
                second = source.second

        # A first-order dual has no information from which curvature can be
        # recovered.  Refuse it even if a caller happens to supply a second
        # argument; an explicit interval-ring value is the supported bridge.
        if isinstance(value, AxialDual) or isinstance(tangent, AxialDual) or isinstance(second, AxialDual):
            raise TypeError(
                "AxialDual cannot be promoted to IntervalAxialSecondJet; "
                "provide an IntervalPressureWidthJet value and explicit slots"
            )

        if tangent is _MISSING:
            tangent = 0
        if second is _MISSING:
            second = 0

        ring_slots = tuple(
            slot
            for slot in (value, tangent, second)
            if isinstance(slot, IntervalPressureWidthJet)
        )
        inferred_orders = ring_slots[0].orders if ring_slots else None
        if any(slot.ctx is not ctx for slot in ring_slots):
            raise ValueError("Interval pressure/width contexts must match")
        if any(slot.orders != inferred_orders for slot in ring_slots):
            raise ValueError(
                "IntervalAxialSecondJet value, tangent, and second orders must match"
            )

        if pressure_order is None:
            if inferred_orders is None:
                raise TypeError(
                    "pressure_order is required when no interval-ring slot is supplied"
                )
            pressure = inferred_orders[0]
        else:
            pressure = _validate_order(pressure_order, "pressure_order")
        if width_order is None:
            if inferred_orders is None:
                raise TypeError(
                    "width_order is required when no interval-ring slot is supplied"
                )
            width = inferred_orders[1]
        else:
            width = _validate_order(width_order, "width_order")
        orders = (pressure, width)
        if inferred_orders is not None and inferred_orders != orders:
            raise ValueError(
                "IntervalAxialSecondJet slot orders do not match declared orders: "
                f"{inferred_orders} != {orders}"
            )

        self.ctx = ctx
        self.pressure_order, self.width_order = orders
        self.value = self._coerce_slot(value, "value")
        self.tangent = self._coerce_slot(tangent, "tangent")
        self.second = self._coerce_slot(second, "second")

    def _coerce_slot(self, slot: Any, name: str) -> IntervalPressureWidthJet:
        if isinstance(slot, AxialDual):
            raise TypeError(
                f"AxialDual cannot supply the IntervalAxialSecondJet {name} slot"
            )
        if isinstance(slot, IntervalAxialSecondJet):
            raise TypeError(
                f"IntervalAxialSecondJet cannot supply the {name} slot; use its ring slot"
            )
        if isinstance(slot, IntervalPressureWidthJet):
            if slot.ctx is not self.ctx:
                raise ValueError("Interval pressure/width contexts must match")
            if slot.orders != self.orders:
                raise ValueError(
                    "IntervalAxialSecondJet slots must use the declared orders"
                )
            return slot
        return IntervalPressureWidthJet(
            self.ctx,
            slot,
            pressure_order=self.pressure_order,
            width_order=self.width_order,
        )

    @classmethod
    def _from_slots(
        cls,
        ctx: Any,
        value: IntervalPressureWidthJet,
        tangent: IntervalPressureWidthJet,
        second: IntervalPressureWidthJet,
    ) -> "IntervalAxialSecondJet":
        """Build a result without coercing any slot."""

        slots = (value, tangent, second)
        if not all(isinstance(slot, IntervalPressureWidthJet) for slot in slots):
            raise TypeError("IntervalAxialSecondJet slots must be interval-ring values")
        if any(slot.ctx is not ctx for slot in slots):
            raise ValueError("Interval pressure/width contexts must match")
        if not (value.orders == tangent.orders == second.orders):
            raise ValueError("IntervalAxialSecondJet slots must have matching orders")
        result = cls.__new__(cls)
        result.ctx = ctx
        result.pressure_order, result.width_order = value.orders
        result.value = value
        result.tangent = tangent
        result.second = second
        return result

    @property
    def orders(self) -> tuple[int, int]:
        return self.pressure_order, self.width_order

    @property
    def constant(self):
        return self.value.constant

    @property
    def derivative(self) -> IntervalPressureWidthJet:
        return self.tangent

    @property
    def second_derivative(self) -> IntervalPressureWidthJet:
        return self.second

    @property
    def curvature(self) -> IntervalPressureWidthJet:
        return self.second

    def component(self, pressure_power: Any, width_power: Any = 0):
        return self.value.component(pressure_power, width_power)

    def tangent_component(self, pressure_power: Any, width_power: Any = 0):
        return self.tangent.component(pressure_power, width_power)

    def second_component(self, pressure_power: Any, width_power: Any = 0):
        return self.second.component(pressure_power, width_power)

    curvature_component = second_component

    def evaluate(self, pressure: Any = 1, width: Any = 1) -> "IntervalAxialSecondJet":
        """Evaluate all retained pressure/width atoms at explicit parameters."""

        orders = self._order_kwargs()
        return self._from_slots(
            self.ctx,
            IntervalPressureWidthJet(self.ctx, self.value.evaluate(pressure, width), **orders),
            IntervalPressureWidthJet(self.ctx, self.tangent.evaluate(pressure, width), **orders),
            IntervalPressureWidthJet(self.ctx, self.second.evaluate(pressure, width), **orders),
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
            "directed_interval_arithmetic": True,
        }

    def __repr__(self) -> str:
        return (
            f"IntervalAxialSecondJet(value={self.value!r}, "
            f"tangent={self.tangent!r}, second={self.second!r})"
        )

    def _order_kwargs(self) -> dict[str, int]:
        return {
            "pressure_order": self.pressure_order,
            "width_order": self.width_order,
        }

    def _coerce_operand(self, other: Any) -> "IntervalAxialSecondJet":
        if isinstance(other, IntervalAxialSecondJet):
            if other.ctx is not self.ctx:
                raise ValueError("Interval axial jet contexts must match")
            if other.orders != self.orders:
                raise ValueError(
                    f"IntervalAxialSecondJet orders must match: {self.orders} != {other.orders}"
                )
            return other
        if isinstance(other, AxialDual):
            raise TypeError(
                "AxialDual cannot be promoted to IntervalAxialSecondJet without explicit slots"
            )
        if isinstance(other, IntervalPressureWidthJet):
            if other.ctx is not self.ctx:
                raise ValueError("Interval pressure/width contexts must match")
            if other.orders != self.orders:
                raise ValueError(
                    "IntervalAxialSecondJet base orders must match: "
                    f"{self.orders} != {other.orders}"
                )
            zero = IntervalPressureWidthJet(self.ctx, 0, **self._order_kwargs())
            return self._from_slots(self.ctx, other, zero, zero)
        return IntervalAxialSecondJet(
            self.ctx,
            other,
            pressure_order=self.pressure_order,
            width_order=self.width_order,
        )

    def __add__(self, other: Any) -> "IntervalAxialSecondJet":
        right = self._coerce_operand(other)
        return self._from_slots(
            self.ctx,
            self.value + right.value,
            self.tangent + right.tangent,
            self.second + right.second,
        )

    __radd__ = __add__

    def __neg__(self) -> "IntervalAxialSecondJet":
        return self._from_slots(self.ctx, -self.value, -self.tangent, -self.second)

    def __sub__(self, other: Any) -> "IntervalAxialSecondJet":
        return self + (-self._coerce_operand(other))

    def __rsub__(self, other: Any) -> "IntervalAxialSecondJet":
        return self._coerce_operand(other) - self

    def __mul__(self, other: Any) -> "IntervalAxialSecondJet":
        right = self._coerce_operand(other)
        return self._from_slots(
            self.ctx,
            self.value * right.value,
            self.tangent * right.value + self.value * right.tangent,
            self.second * right.value
            + 2 * self.tangent * right.tangent
            + self.value * right.second,
        )

    __rmul__ = __mul__

    def reciprocal(self) -> "IntervalAxialSecondJet":
        """Return ``1 / self`` using actual first and second chain rules."""

        value_reciprocal = self.value.reciprocal()
        value_reciprocal_squared = value_reciprocal**2
        value_reciprocal_cubed = value_reciprocal**3
        return self._from_slots(
            self.ctx,
            value_reciprocal,
            -self.tangent * value_reciprocal_squared,
            2 * self.tangent**2 * value_reciprocal_cubed
            - self.second * value_reciprocal_squared,
        )

    inverse = reciprocal
    inv = reciprocal

    def __truediv__(self, other: Any) -> "IntervalAxialSecondJet":
        return self * self._coerce_operand(other).reciprocal()

    def __rtruediv__(self, other: Any) -> "IntervalAxialSecondJet":
        return self._coerce_operand(other) * self.reciprocal()

    def __pow__(self, exponent: Any) -> "IntervalAxialSecondJet":
        integer = _integer_exponent(exponent)
        if integer == 0:
            return self._from_slots(
                self.ctx,
                IntervalPressureWidthJet(self.ctx, 1, **self._order_kwargs()),
                IntervalPressureWidthJet(self.ctx, 0, **self._order_kwargs()),
                IntervalPressureWidthJet(self.ctx, 0, **self._order_kwargs()),
            )
        if integer < 0:
            return self.reciprocal() ** (-integer)
        if integer == 1:
            return self._from_slots(self.ctx, self.value, self.tangent, self.second)

        value_power = self.value**integer
        tangent_power = integer * self.value ** (integer - 1) * self.tangent
        second_power = (
            integer * (integer - 1) * self.value ** (integer - 2) * self.tangent**2
            + integer * self.value ** (integer - 1) * self.second
        )
        return self._from_slots(self.ctx, value_power, tangent_power, second_power)

    def exp(self) -> "IntervalAxialSecondJet":
        value_exp = self.value.exp()
        return self._from_slots(
            self.ctx,
            value_exp,
            value_exp * self.tangent,
            value_exp * (self.second + self.tangent**2),
        )

    def log(self) -> "IntervalAxialSecondJet":
        value_log = self.value.log()
        return self._from_slots(
            self.ctx,
            value_log,
            self.tangent / self.value,
            self.second / self.value - (self.tangent / self.value) ** 2,
        )

    def sqrt(self) -> "IntervalAxialSecondJet":
        value_sqrt = self.value.sqrt()
        return self._from_slots(
            self.ctx,
            value_sqrt,
            self.tangent / (2 * value_sqrt),
            self.second / (2 * value_sqrt)
            - self.tangent**2 / (4 * value_sqrt**3),
        )


# Alternate names make the interval qualifier discoverable without changing
# the one implementation or its strict context/order checks.
AxialIntervalSecondJet = IntervalAxialSecondJet
IntervalSecondAxialJet = IntervalAxialSecondJet


def main() -> None:
    print(
        json.dumps(
            {
                "class": "IntervalAxialSecondJet",
                "pressure_order": DEFAULT_PRESSURE_ORDER,
                "width_order": DEFAULT_WIDTH_ORDER,
                "actual_derivative_slots": True,
                "directed_interval_arithmetic": True,
            },
            indent=2,
            sort_keys=True,
        )
    )


__all__ = [
    "IntervalAxialSecondJet",
    "AxialIntervalSecondJet",
    "IntervalSecondAxialJet",
]


if __name__ == "__main__":
    main()
