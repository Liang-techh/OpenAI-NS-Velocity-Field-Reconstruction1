"""First-order axial automatic differentiation over ``PressureWidthJet``.

``AxialDual`` stores a value and its first derivative with respect to the
axial variable ``Z``.  Both slots are fixed-order
``PressureWidthJet`` instances, so the formal pressure-tail and actual
exit-width atoms remain separate while the axial derivative follows the
product, quotient, and chain rules.

The wrapper deliberately does not expose a top-level ``atoms`` mapping.  A
dual value is distinct from its base bivariate jet; callers that need a
coefficient use ``component(p, w)`` for the value or
``tangent_component(p, w)`` for its axial derivative.
"""

from __future__ import annotations

import json
import operator
from typing import Any

import mpmath as mp

try:  # package import
    from .lei_ren_part1_paper_pressure_width_jet import PressureWidthJet
except (ImportError, ValueError):  # direct experiment execution
    from lei_ren_part1_paper_pressure_width_jet import PressureWidthJet


DEFAULT_PRESSURE_ORDER = 9
DEFAULT_WIDTH_ORDER = 2
AUTOMATIC_FIRST_Z_DERIVATIVE = True


def _is_zero_scalar(value: Any) -> bool:
    try:
        return value == 0
    except Exception:
        return False


def _integer_exponent(exponent: Any) -> int:
    try:
        result = operator.index(exponent)
    except TypeError:
        try:
            candidate = mp.mpf(exponent)
        except (TypeError, ValueError) as exc:
            raise ValueError("AxialDual exponent must be an integer") from exc
        if not mp.isint(candidate):
            raise ValueError("AxialDual exponent must be an integer")
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
    if isinstance(value, AxialDual):
        return value.orders
    return None


class AxialDual:
    """A value and first ``Z`` derivative over a fixed bivariate jet ring.

    If either input is already a ``PressureWidthJet``, its orders are inferred
    unless explicit orders are supplied.  Scalar-only construction uses the
    effective defaults ``pressure_order=9`` and ``width_order=2``.
    """

    automatic_first_Z_derivative = AUTOMATIC_FIRST_Z_DERIVATIVE

    def __init__(
        self,
        value: Any = 0,
        tangent: Any = 0,
        *,
        pressure_order: int | None = None,
        width_order: int | None = None,
    ):
        # Preserve a dual passed as the value when no replacement tangent is
        # supplied.  This keeps construction useful in generic coercion code.
        if isinstance(value, AxialDual) and _is_zero_scalar(tangent):
            inferred = value.orders
            pressure = (
                inferred[0]
                if pressure_order is None
                else _validate_order(pressure_order, "pressure_order")
            )
            width = (
                inferred[1]
                if width_order is None
                else _validate_order(width_order, "width_order")
            )
            if (pressure, width) != inferred:
                raise ValueError(
                    f"AxialDual orders must match: {(pressure, width)} != {inferred}"
                )
            self.pressure_order = pressure
            self.width_order = width
            self.value = value.value
            self.tangent = value.tangent
            return

        value_orders = _jet_orders(value)
        tangent_orders = _jet_orders(tangent)
        inferred_orders = value_orders or tangent_orders
        if value_orders is not None and tangent_orders is not None:
            if value_orders != tangent_orders:
                raise ValueError(
                    "AxialDual value and tangent orders must match: "
                    f"{value_orders} != {tangent_orders}"
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
                inferred_orders[1] if inferred_orders is not None else DEFAULT_WIDTH_ORDER
            )
        else:
            width = _validate_order(width_order, "width_order")

        orders = dict(pressure_order=pressure, width_order=width)
        self.pressure_order = pressure
        self.width_order = width
        self.value = value if isinstance(value, PressureWidthJet) else PressureWidthJet(value, **orders)
        self.tangent = (
            tangent
            if isinstance(tangent, PressureWidthJet)
            else PressureWidthJet(tangent, **orders)
        )
        if self.value.orders != (pressure, width) or self.tangent.orders != (pressure, width):
            raise ValueError("AxialDual base values must use the declared orders")

    @classmethod
    def _from_slots(
        cls,
        value: PressureWidthJet,
        tangent: PressureWidthJet,
    ) -> "AxialDual":
        result = cls.__new__(cls)
        result.pressure_order, result.width_order = value.orders
        if tangent.orders != value.orders:
            raise ValueError("AxialDual slots must have matching orders")
        result.value = value
        result.tangent = tangent
        return result

    @property
    def orders(self) -> tuple[int, int]:
        return self.pressure_order, self.width_order

    @property
    def derivative(self) -> PressureWidthJet:
        return self.tangent

    def component(self, pressure_power: Any, width_power: Any = 0) -> mp.mpf:
        """Return a value coefficient; use ``tangent_component`` for ``Z``."""

        return self.value.component(pressure_power, width_power)

    def tangent_component(self, pressure_power: Any, width_power: Any = 0) -> mp.mpf:
        return self.tangent.component(pressure_power, width_power)

    def evaluate(self, pressure: Any = 1, width: Any = 1) -> "AxialDual":
        """Evaluate the two base jets while retaining the axial dual slots."""

        return self._from_slots(
            PressureWidthJet(
                self.value.evaluate(pressure, width),
                pressure_order=self.pressure_order,
                width_order=self.width_order,
            ),
            PressureWidthJet(
                self.tangent.evaluate(pressure, width),
                pressure_order=self.pressure_order,
                width_order=self.width_order,
            ),
        )

    __call__ = evaluate

    def as_dict(self) -> dict[str, Any]:
        return {
            "pressure_order": self.pressure_order,
            "width_order": self.width_order,
            "value": self.value.as_dict(),
            "tangent": self.tangent.as_dict(),
            "automatic_first_Z_derivative": True,
        }

    def __repr__(self) -> str:
        return (
            f"AxialDual(value={self.value!r}, tangent={self.tangent!r})"
        )

    def _coerce_operand(self, other: Any) -> "AxialDual":
        if isinstance(other, AxialDual):
            if other.orders != self.orders:
                raise ValueError(
                    f"AxialDual orders must match: {self.orders} != {other.orders}"
                )
            return other
        if isinstance(other, PressureWidthJet):
            if other.orders != self.orders:
                raise ValueError(
                    f"AxialDual base orders must match: {self.orders} != {other.orders}"
                )
            return self._from_slots(other, PressureWidthJet(0, **self._order_kwargs()))
        return AxialDual(other, pressure_order=self.pressure_order, width_order=self.width_order)

    def _order_kwargs(self) -> dict[str, int]:
        return {
            "pressure_order": self.pressure_order,
            "width_order": self.width_order,
        }

    def __add__(self, other: Any) -> "AxialDual":
        right = self._coerce_operand(other)
        return self._from_slots(
            self.value + right.value,
            self.tangent + right.tangent,
        )

    __radd__ = __add__

    def __neg__(self) -> "AxialDual":
        return self._from_slots(-self.value, -self.tangent)

    def __sub__(self, other: Any) -> "AxialDual":
        return self + (-self._coerce_operand(other))

    def __rsub__(self, other: Any) -> "AxialDual":
        return self._coerce_operand(other) - self

    def __mul__(self, other: Any) -> "AxialDual":
        right = self._coerce_operand(other)
        return self._from_slots(
            self.value * right.value,
            self.tangent * right.value + self.value * right.tangent,
        )

    __rmul__ = __mul__

    def reciprocal(self) -> "AxialDual":
        value_reciprocal = self.value.reciprocal()
        return self._from_slots(
            value_reciprocal,
            -self.tangent * value_reciprocal**2,
        )

    inverse = reciprocal
    inv = reciprocal

    def __truediv__(self, other: Any) -> "AxialDual":
        return self * self._coerce_operand(other).reciprocal()

    def __rtruediv__(self, other: Any) -> "AxialDual":
        return self._coerce_operand(other) * self.reciprocal()

    def __pow__(self, exponent: Any) -> "AxialDual":
        integer = _integer_exponent(exponent)
        if integer == 0:
            return self._from_slots(
                PressureWidthJet(1, **self._order_kwargs()),
                PressureWidthJet(0, **self._order_kwargs()),
            )
        if integer < 0:
            return (self.reciprocal()) ** (-integer)
        value_power = self.value**integer
        tangent_power = (
            self.tangent
            * self.value ** (integer - 1)
            * integer
        )
        return self._from_slots(value_power, tangent_power)

    def exp(self) -> "AxialDual":
        value_exp = self.value.exp()
        return self._from_slots(value_exp, value_exp * self.tangent)

    def log(self) -> "AxialDual":
        return self._from_slots(self.value.log(), self.tangent / self.value)

    def sqrt(self) -> "AxialDual":
        value_sqrt = self.value.sqrt()
        return self._from_slots(value_sqrt, self.tangent / (2 * value_sqrt))


def _max_scaled_error(actual: PressureWidthJet, expected: PressureWidthJet) -> mp.mpf:
    errors = []
    for pressure in range(actual.pressure_order + 1):
        for width in range(actual.width_order + 1):
            difference = abs(actual.component(pressure, width) - expected.component(pressure, width))
            errors.append(difference / max(mp.mpf(1), abs(expected.component(pressure, width))))
    return max(errors, default=mp.mpf(0))


def run_fixtures(*, precision: int = 110) -> dict[str, Any]:
    """Check dual product/chain rules at resolved and separated scales."""

    with mp.workdps(max(60, int(precision))):
        orders = dict(pressure_order=3, width_order=2)
        value = PressureWidthJet(
            {
                (0, 0): mp.mpf("1.7"),
                (1, 0): mp.mpf("-0.11"),
                (0, 1): mp.mpf("0.08"),
                (1, 1): mp.mpf("0.015"),
                (2, 0): mp.mpf("0.022"),
                (0, 2): mp.mpf("0.01"),
            },
            **orders,
        )
        tangent = PressureWidthJet(
            {
                (0, 0): mp.mpf("0.23"),
                (1, 0): mp.mpf("0.04"),
                (0, 1): mp.mpf("-0.03"),
                (1, 1): mp.mpf("0.006"),
                (2, 0): mp.mpf("-0.008"),
                (0, 2): mp.mpf("0.004"),
            },
            **orders,
        )
        dual = AxialDual(value, tangent)

        quotient = (dual * dual + 2) / (dual + 3)
        quotient_value = (value * value + 2) / (value + 3)
        quotient_tangent = (
            (2 * value * tangent) * (value + 3)
            - (value * value + 2) * tangent
        ) / (value + 3) ** 2
        product_ok = _max_scaled_error(quotient.value, quotient_value)
        product_tangent_ok = _max_scaled_error(quotient.tangent, quotient_tangent)

        exp_dual = dual.exp()
        log_dual = dual.log()
        sqrt_dual = dual.sqrt()
        exp_value_error = _max_scaled_error(exp_dual.value, value.exp())
        exp_tangent_error = _max_scaled_error(exp_dual.tangent, value.exp() * tangent)
        log_value_error = _max_scaled_error(log_dual.value, value.log())
        log_tangent_error = _max_scaled_error(log_dual.tangent, tangent / value)
        sqrt_value_error = _max_scaled_error(sqrt_dual.value, value.sqrt())
        sqrt_tangent_error = _max_scaled_error(
            sqrt_dual.tangent, tangent / (2 * value.sqrt())
        )

        tiny_pressure = mp.exp(-mp.mpf("1e28"))
        tiny_width = mp.exp(-mp.mpf("1e154"))
        mixed = mp.mpf(3) * tiny_pressure * tiny_width
        tiny_value = PressureWidthJet(
            {
                (0, 0): mp.mpf(2),
                (1, 0): tiny_pressure,
                (0, 1): tiny_width,
                (1, 1): mixed,
            },
            pressure_order=3,
            width_order=2,
        )
        tiny_tangent = PressureWidthJet(
            {
                (0, 0): mp.mpf(5),
                (1, 0): mp.mpf(2) * tiny_pressure,
                (0, 1): mp.mpf(3) * tiny_width,
                (1, 1): mp.mpf(7) * mixed,
            },
            pressure_order=3,
            width_order=2,
        )
        tiny_dual = AxialDual(tiny_value, tiny_tangent)
        tiny_exp = tiny_dual.exp()
        tiny_checks = {
            "value_pressure_atom": tiny_dual.component(1, 0) == tiny_pressure,
            "value_width_atom": tiny_dual.component(0, 1) == tiny_width,
            "value_mixed_atom": tiny_dual.component(1, 1) == mixed,
            "tangent_pressure_atom": tiny_dual.tangent_component(1, 0)
            == 2 * tiny_pressure,
            "tangent_width_atom": tiny_dual.tangent_component(0, 1)
            == 3 * tiny_width,
            "tangent_mixed_atom": tiny_dual.tangent_component(1, 1)
            == 7 * mixed,
            "exp_tangent_pressure_atom": tiny_exp.tangent_component(1, 0)
            == mp.exp(2) * 2 * tiny_pressure
            + mp.exp(2) * tiny_pressure * 5,
            "exp_tangent_width_atom": tiny_exp.tangent_component(0, 1)
            == mp.exp(2) * 3 * tiny_width
            + mp.exp(2) * tiny_width * 5,
        }
        tolerance = mp.mpf(10) ** (-(max(60, int(precision)) // 2))
        checks = {
            "quotient_value": product_ok <= tolerance,
            "quotient_tangent": product_tangent_ok <= tolerance,
            "exp_value": exp_value_error <= tolerance,
            "exp_tangent": exp_tangent_error <= tolerance,
            "log_value": log_value_error <= tolerance,
            "log_tangent": log_tangent_error <= tolerance,
            "sqrt_value": sqrt_value_error <= tolerance,
            "sqrt_tangent": sqrt_tangent_error <= tolerance,
            "tiny_atoms": all(tiny_checks.values()),
        }
        return {
            "passed": all(checks.values()),
            "checks": checks,
            "automatic_first_Z_derivative": True,
            "pressure_order": 3,
            "width_order": 2,
            "quotient_value_max_scaled_error": mp.nstr(product_ok, 12),
            "quotient_tangent_max_scaled_error": mp.nstr(product_tangent_ok, 12),
            "exp_value_max_scaled_error": mp.nstr(exp_value_error, 12),
            "exp_tangent_max_scaled_error": mp.nstr(exp_tangent_error, 12),
            "log_value_max_scaled_error": mp.nstr(log_value_error, 12),
            "log_tangent_max_scaled_error": mp.nstr(log_tangent_error, 12),
            "sqrt_value_max_scaled_error": mp.nstr(sqrt_value_error, 12),
            "sqrt_tangent_max_scaled_error": mp.nstr(sqrt_tangent_error, 12),
            "tiny_checks": tiny_checks,
            "limitations": [
                "The wrapped pressure-width ring remains rectangular and finite-order.",
                "AxialDual carries first Z derivatives only; it is not a global ODE error enclosure.",
                "No field installation, prescribed bridge, or Navier--Stokes closure is claimed.",
            ],
        }


def run_fixture(*, precision: int = 110) -> dict[str, Any]:
    return run_fixtures(precision=precision)


def fixture(*, precision: int = 110) -> dict[str, Any]:
    return run_fixture(precision=precision)


def main() -> None:
    print(json.dumps(run_fixtures(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
