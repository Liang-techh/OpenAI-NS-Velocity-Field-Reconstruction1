"""Finite Taylor algebra for the pressure-tail parameter ``lambda``.

The component radial core uses a formal pressure-tail polynomial.  The exit
ODE needs a small, bounded Taylor ring when that polynomial is carried through
additional nonlinear operations.  :class:`PressureParameterJet` stores the
ordinary (derivative divided by ``k!``) Taylor coefficients around
``lambda = 0`` and truncates every operation at a declared order.

The coefficients remain independent ``mpmath`` atoms.  In particular, a term
such as ``mp.exp(-1e28)`` is not converted to a Python float and is not hidden
by evaluating the whole jet at ``lambda = 1``.  The nominal evaluation is
provided for callers that need it, while ``component(k)`` is the diagnostic
interface for separated scales.

This is a finite Taylor ring, rather than exact rational-function arithmetic.
The module therefore reports ``pressure_parameter_order_truncated=True`` and
does not claim exact propagation of pressure-tail rational functions or any
global Navier--Stokes closure.
"""

from __future__ import annotations

from collections.abc import Mapping
import json
import operator
from typing import Any

import mpmath as mp


DEFAULT_ORDER = 9
PRESSURE_PARAMETER_ORDER_TRUNCATED = True
NORMALIZED_COEFFICIENT_CONVENTION = "coefficient of lambda**k (derivative divided by k!)"


def _as_mp_scalar(value: Any) -> mp.mpf | mp.mpc:
    """Convert a scalar without routing it through binary floating point."""

    if isinstance(value, (mp.mpf, mp.mpc)):
        return value
    try:
        return mp.mpf(value)
    except (TypeError, ValueError):
        # mpmath's complex scalar is useful for the algebraic log branch.  Do
        # not catch arbitrary exceptions here: malformed atom values should
        # fail at the call site instead of being silently stringified.
        return mp.mpc(value)


def _validate_order(order: Any) -> int:
    try:
        n = operator.index(order)
    except TypeError as exc:
        raise TypeError("Taylor order must be an integer") from exc
    if n < 0:
        raise ValueError("Taylor order must be nonnegative")
    return int(n)


def _power_key(key: Any) -> int:
    try:
        power = operator.index(key)
    except TypeError:
        try:
            candidate = mp.mpf(key)
        except (TypeError, ValueError) as exc:
            raise TypeError(f"Taylor atom power must be an integer: {key!r}") from exc
        if not mp.isint(candidate):
            raise TypeError(f"Taylor atom power must be an integer: {key!r}")
        power = int(candidate)
    power = int(power)
    if power < 0:
        raise ValueError("Taylor atom powers must be nonnegative")
    return power


def _is_zero(value: mp.mpf | mp.mpc) -> bool:
    # Equality is intentional.  An MP atom with an exponent such as -1e28 is
    # nonzero and must remain present in the sparse atom dictionary.
    return value == 0


class PressureParameterJet:
    """A finite normalized Taylor jet in the pressure-tail parameter.

    ``value`` may be an MP scalar, a mapping from nonnegative powers to MP
    coefficients, or an existing pressure polynomial exposing ``.atoms``.
    Every operation has the same order as its left operand.  A scalar on
    either side is promoted to that order; two jets with different orders are
    rejected so that accidental precision/order mixing is visible.
    """

    pressure_parameter_order_truncated = PRESSURE_PARAMETER_ORDER_TRUNCATED
    normalized_coefficient_convention = NORMALIZED_COEFFICIENT_CONVENTION

    def __init__(self, value: Any = 0, *, order: int = DEFAULT_ORDER):
        self.order = _validate_order(order)
        self._truncated_input_powers: tuple[int, ...] = ()

        if isinstance(value, PressureParameterJet):
            if value.order != self.order:
                raise ValueError(
                    f"Taylor orders must match: {self.order} != {value.order}"
                )
            source: Any = value.atoms
        elif hasattr(value, "atoms") and isinstance(getattr(value, "atoms"), Mapping):
            # This is deliberately duck typed.  It accepts the existing
            # PressurePolynomial without importing or changing that module.
            source = value.atoms
        elif isinstance(value, Mapping):
            source = value
        else:
            source = {0: value}

        atoms: dict[int, mp.mpf | mp.mpc] = {}
        discarded: list[int] = []
        for raw_power, raw_coefficient in source.items():
            power = _power_key(raw_power)
            if power > self.order:
                discarded.append(power)
                continue
            coefficient = _as_mp_scalar(raw_coefficient)
            if not _is_zero(coefficient):
                atoms[power] = coefficient
        self.atoms = atoms
        self._truncated_input_powers = tuple(sorted(set(discarded)))

    @classmethod
    def _from_coefficients(
        cls,
        coefficients: list[mp.mpf | mp.mpc] | tuple[mp.mpf | mp.mpc, ...],
        order: int,
    ) -> "PressureParameterJet":
        """Construct an operation result without reinterpreting its atoms."""

        result = cls.__new__(cls)
        result.order = order
        result._truncated_input_powers = ()
        result.atoms = {
            power: coefficient
            for power, coefficient in enumerate(coefficients[: order + 1])
            if not _is_zero(coefficient)
        }
        return result

    def _coefficients(self) -> list[mp.mpf | mp.mpc]:
        return [self.component(power) for power in range(self.order + 1)]

    def _coerce_operand(self, other: Any) -> "PressureParameterJet":
        if isinstance(other, PressureParameterJet):
            if other.order != self.order:
                raise ValueError(
                    f"Taylor orders must match: {self.order} != {other.order}"
                )
            return other
        return PressureParameterJet(other, order=self.order)

    @property
    def constant(self) -> mp.mpf | mp.mpc:
        return self.component(0)

    @property
    def truncated_input_powers(self) -> tuple[int, ...]:
        return self._truncated_input_powers

    def component(self, power: Any) -> mp.mpf | mp.mpc:
        """Return the normalized coefficient of ``lambda**power``."""

        try:
            index = operator.index(power)
        except TypeError:
            candidate = mp.mpf(power)
            if not mp.isint(candidate):
                raise TypeError("Taylor component index must be an integer")
            index = int(candidate)
        index = int(index)
        if index < 0 or index > self.order:
            return mp.mpf(0)
        return self.atoms.get(index, mp.mpf(0))

    coefficient = component

    def evaluate(self, parameter: Any = 1) -> mp.mpf | mp.mpc:
        """Evaluate the retained polynomial at a supplied parameter.

        Component atoms should be used for diagnostics when scales differ by
        more than the working precision.  A scalar evaluation necessarily
        obeys the rounding of the active ``mpmath`` context.
        """

        p = _as_mp_scalar(parameter)
        return mp.fsum(
            self.component(power) * p**power for power in range(self.order + 1)
        )

    __call__ = evaluate

    def as_dict(self) -> dict[str, Any]:
        """Return a JSON-friendly receipt without replacing atom values."""

        return {
            "order": self.order,
            "atoms": {
                str(power): mp.nstr(coefficient, 30)
                for power, coefficient in sorted(self.atoms.items())
            },
            "pressure_parameter_order_truncated": True,
            "normalized_coefficient_convention": NORMALIZED_COEFFICIENT_CONVENTION,
            "truncated_input_powers": list(self.truncated_input_powers),
        }

    def __repr__(self) -> str:
        return (
            f"PressureParameterJet({self.as_dict()['atoms']!r}, "
            f"order={self.order})"
        )

    def __add__(self, other: Any) -> "PressureParameterJet":
        right = self._coerce_operand(other)
        coefficients = [
            self.component(power) + right.component(power)
            for power in range(self.order + 1)
        ]
        return self._from_coefficients(coefficients, self.order)

    __radd__ = __add__

    def __neg__(self) -> "PressureParameterJet":
        return self._from_coefficients(
            [-self.component(power) for power in range(self.order + 1)], self.order
        )

    def __sub__(self, other: Any) -> "PressureParameterJet":
        return self + (-self._coerce_operand(other))

    def __rsub__(self, other: Any) -> "PressureParameterJet":
        return self._coerce_operand(other) - self

    def __mul__(self, other: Any) -> "PressureParameterJet":
        right = self._coerce_operand(other)
        coefficients: list[mp.mpf | mp.mpc] = []
        for power in range(self.order + 1):
            coefficients.append(
                mp.fsum(
                    self.component(left)
                    * right.component(power - left)
                    for left in range(power + 1)
                )
            )
        return self._from_coefficients(coefficients, self.order)

    __rmul__ = __mul__

    def reciprocal(self) -> "PressureParameterJet":
        """Return the finite series reciprocal, requiring a nonzero constant."""

        constant = self.constant
        if _is_zero(constant):
            raise ZeroDivisionError(
                "PressureParameterJet reciprocal requires a nonzero constant term"
            )
        coefficients: list[mp.mpf | mp.mpc] = [1 / constant]
        for power in range(1, self.order + 1):
            convolution = mp.fsum(
                self.component(index) * coefficients[power - index]
                for index in range(1, power + 1)
            )
            coefficients.append(-convolution / constant)
        return self._from_coefficients(coefficients, self.order)

    inverse = reciprocal
    inv = reciprocal

    def __truediv__(self, other: Any) -> "PressureParameterJet":
        return self * self._coerce_operand(other).reciprocal()

    def __rtruediv__(self, other: Any) -> "PressureParameterJet":
        return self._coerce_operand(other) * self.reciprocal()

    def __pow__(self, exponent: Any) -> "PressureParameterJet":
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
            return (self.reciprocal()) ** (-integer)
        result = PressureParameterJet(1, order=self.order)
        base = self
        while integer:
            if integer & 1:
                result = result * base
            integer //= 2
            if integer:
                base = base * base
        return result

    def exp(self) -> "PressureParameterJet":
        """Return ``exp(self)`` through the retained Taylor order."""

        coefficients: list[mp.mpf | mp.mpc] = [mp.exp(self.constant)]
        for power in range(1, self.order + 1):
            numerator = mp.fsum(
                index
                * self.component(index)
                * coefficients[power - index]
                for index in range(1, power + 1)
            )
            coefficients.append(numerator / power)
        return self._from_coefficients(coefficients, self.order)

    def log(self) -> "PressureParameterJet":
        """Return ``log(self)`` through the retained Taylor order."""

        constant = self.constant
        if _is_zero(constant):
            raise ZeroDivisionError(
                "PressureParameterJet log requires a nonzero constant term"
            )
        coefficients: list[mp.mpf | mp.mpc] = [mp.log(constant)]
        for power in range(1, self.order + 1):
            convolution = mp.fsum(
                self.component(index)
                * (power - index)
                * coefficients[power - index]
                for index in range(1, power)
            )
            coefficients.append(
                (power * self.component(power) - convolution) / (power * constant)
            )
        return self._from_coefficients(coefficients, self.order)


def _max_scaled_error(
    actual: PressureParameterJet,
    expected: list[mp.mpf | mp.mpc],
) -> tuple[mp.mpf, mp.mpf]:
    absolute = [
        abs(actual.component(power) - expected[power])
        for power in range(actual.order + 1)
    ]
    scaled = [
        abs(actual.component(power) - expected[power])
        / max(mp.mpf(1), abs(expected[power]))
        for power in range(actual.order + 1)
    ]
    return max(absolute, default=mp.mpf(0)), max(scaled, default=mp.mpf(0))


def run_fixtures(*, order: int = DEFAULT_ORDER, precision: int = 110) -> dict[str, Any]:
    """Run independent resolved-scale checks for the finite Taylor ring.

    ``mp.taylor`` is used as the reference for the three nonlinear
    functions.  The final check deliberately inspects a component rather
    than a nominal scalar sum, because the latter cannot display an atom
    below the active precision once it is added to a unit constant.
    """

    order = _validate_order(order)
    precision = max(40, int(precision))
    with mp.workdps(precision):
        inverse_source = PressureParameterJet(
            {
                0: mp.mpf("1.7"),
                1: mp.mpf("-0.31"),
                2: mp.mpf("0.047"),
                3: mp.mpf("-0.009"),
            },
            order=order,
        )
        log_source = PressureParameterJet(
            {
                0: mp.mpf("1.3"),
                1: mp.mpf("0.22"),
                2: mp.mpf("-0.04"),
                3: mp.mpf("0.006"),
            },
            order=order,
        )
        exp_source = PressureParameterJet(
            {
                0: mp.mpf("-0.6"),
                1: mp.mpf("0.41"),
                2: mp.mpf("-0.08"),
                3: mp.mpf("0.015"),
            },
            order=order,
        )

        def polynomial(source: PressureParameterJet, parameter: mp.mpf) -> mp.mpf:
            return mp.fsum(
                source.component(power) * parameter**power
                for power in range(source.order + 1)
            )

        expected_inverse = mp.taylor(
            lambda parameter: 1 / polynomial(inverse_source, parameter), 0, order
        )
        expected_log = mp.taylor(
            lambda parameter: mp.log(polynomial(log_source, parameter)), 0, order
        )
        expected_exp = mp.taylor(
            lambda parameter: mp.exp(polynomial(exp_source, parameter)), 0, order
        )

        inverse_abs, inverse_scaled = _max_scaled_error(
            inverse_source.reciprocal(), expected_inverse
        )
        log_abs, log_scaled = _max_scaled_error(log_source.log(), expected_log)
        exp_abs, exp_scaled = _max_scaled_error(exp_source.exp(), expected_exp)

        tiny = mp.exp(-mp.mpf("1e28"))
        tiny_source = PressureParameterJet({0: mp.mpf(1), 1: tiny}, order=order)
        tiny_exp = tiny_source.exp()
        tiny_preserved = tiny_source.component(1) == tiny
        # ``exp(1 + tiny*lambda)`` has first atom ``e*tiny``.  Compare that
        # atom directly; comparing with ``tiny`` would miss the analytic
        # constant multiplier while still exercising the same extreme scale.
        tiny_exp_preserved = tiny_exp.component(1) == mp.exp(1) * tiny

        tolerance = mp.mpf(10) ** (-(precision // 2))
        checks = {
            "inverse": inverse_scaled <= tolerance,
            "log": log_scaled <= tolerance,
            "exp": exp_scaled <= tolerance,
            "tiny_atom": tiny_preserved and tiny_exp_preserved,
        }
        return {
            "passed": all(checks.values()),
            "checks": checks,
            "order": order,
            "precision": precision,
            "pressure_parameter_order_truncated": True,
            "exact_rational_propagation": False,
            "normalized_coefficient_convention": NORMALIZED_COEFFICIENT_CONVENTION,
            "inverse_max_absolute_error": mp.nstr(inverse_abs, 12),
            "inverse_max_scaled_error": mp.nstr(inverse_scaled, 12),
            "log_max_absolute_error": mp.nstr(log_abs, 12),
            "log_max_scaled_error": mp.nstr(log_scaled, 12),
            "exp_max_absolute_error": mp.nstr(exp_abs, 12),
            "exp_max_scaled_error": mp.nstr(exp_scaled, 12),
            "tiny_atom": mp.nstr(tiny, 12),
            "tiny_atom_preserved": tiny_preserved,
            "tiny_exp_atom_preserved": tiny_exp_preserved,
            "limitations": [
                "Operations retain powers through the declared finite order only.",
                "The nominal lambda=1 evaluation is a rounded scalar sum and can hide separated atoms.",
                "The module does not certify exact rational propagation, global field closure, or Navier--Stokes compatibility.",
            ],
        }


def run_fixture(*, order: int = DEFAULT_ORDER, precision: int = 110) -> dict[str, Any]:
    """Singular alias retained for small callers and command-line checks."""

    return run_fixtures(order=order, precision=precision)


def fixture(*, order: int = DEFAULT_ORDER, precision: int = 110) -> dict[str, Any]:
    """Compatibility alias for experiment modules that expose ``fixture``."""

    return run_fixtures(order=order, precision=precision)


def main() -> None:
    print(json.dumps(run_fixtures(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
