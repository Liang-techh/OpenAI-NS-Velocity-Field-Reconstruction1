"""Sparse finite Taylor algebra in pressure-tail and exit-width parameters.

``PressureWidthJet`` keeps a rectangular collection of normalized
coefficients

``sum(atoms[(p, w)] * pressure**p * width**w)``.

The first parameter is the formal pressure-tail parameter used by the
component radial core.  The second is the actual small exit-width parameter.
The two parameters are kept as separate atoms, which matters when one is on
the scale of ``exp(-1e28)`` and the other is on the scale of
``exp(-1e154)``.  Every arithmetic result is truncated to the declared
rectangle ``0 <= p <= pressure_order`` and ``0 <= w <= width_order``.

The nonlinear functions use the nilpotent nonconstant part of a jet.  This
keeps the operations coefficient based and avoids first evaluating the whole
jet at ``(pressure, width) = (1, 1)``.  The result is a finite Taylor ring,
not exact rational-function propagation or a global field certificate.
"""

from __future__ import annotations

from collections.abc import Mapping
import json
import operator
from typing import Any

import mpmath as mp


DEFAULT_PRESSURE_ORDER = 3
DEFAULT_WIDTH_ORDER = 2
PRESSURE_WIDTH_RECTANGULAR_TRUNCATED = True
PRESSURE_PARAMETER_ORDER_TRUNCATED = True
EXIT_WIDTH_ORDER_TRUNCATED = True
NORMALIZED_BIVARIATE_CONVENTION = (
    "coefficient of pressure**p * width**w (derivatives divided by p! * w!)"
)


def _as_real_mp(value: Any) -> mp.mpf:
    """Convert a real scalar without routing it through binary floating point."""

    if isinstance(value, mp.mpc):
        if value.imag != 0:
            raise TypeError("PressureWidthJet accepts real MP atoms only")
        return mp.mpf(value.real)
    return mp.mpf(value)


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


def _is_zero(value: mp.mpf) -> bool:
    # Exact equality is intentional: arbitrary MP exponents are unbounded, so
    # exp(-1e154) is a nonzero atom even though it is far below unit scale.
    return value == 0


def _key(raw_key: Any) -> tuple[int, int]:
    """Normalize a sparse key; integer keys mean pressure power, width 0."""

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


class PressureWidthJet:
    """Sparse real bivariate Taylor jet with rectangular truncation."""

    rectangular_truncation = PRESSURE_WIDTH_RECTANGULAR_TRUNCATED
    normalized_coefficient_convention = NORMALIZED_BIVARIATE_CONVENTION
    pressure_parameter_order_truncated = PRESSURE_PARAMETER_ORDER_TRUNCATED
    exit_width_order_truncated = EXIT_WIDTH_ORDER_TRUNCATED

    def __init__(
        self,
        value: Any = 0,
        *,
        pressure_order: int = DEFAULT_PRESSURE_ORDER,
        width_order: int = DEFAULT_WIDTH_ORDER,
    ):
        self.pressure_order = _validate_order(pressure_order, "pressure_order")
        self.width_order = _validate_order(width_order, "width_order")
        self._truncated_input_atoms: tuple[tuple[int, int], ...] = ()

        if isinstance(value, PressureWidthJet):
            if (
                value.pressure_order != self.pressure_order
                or value.width_order != self.width_order
            ):
                raise ValueError(
                    "Bivariate Taylor orders must match: "
                    f"({self.pressure_order}, {self.width_order}) != "
                    f"({value.pressure_order}, {value.width_order})"
                )
            source: Any = value.atoms
        elif hasattr(value, "atoms"):
            # Existing PressurePolynomial and PressureParameterJet objects
            # expose integer pressure powers.  _key maps those to (p, 0).
            source = value.atoms
        elif isinstance(value, Mapping):
            source = value
        else:
            source = {(0, 0): value}

        if not isinstance(source, Mapping):
            raise TypeError("Taylor atoms must be supplied by a mapping")

        atoms: dict[tuple[int, int], mp.mpf] = {}
        discarded: list[tuple[int, int]] = []
        for raw_key, raw_coefficient in source.items():
            key = _key(raw_key)
            pressure, width = key
            if pressure > self.pressure_order or width > self.width_order:
                discarded.append(key)
                continue
            coefficient = _as_real_mp(raw_coefficient)
            if not _is_zero(coefficient):
                atoms[key] = coefficient
        self.atoms = atoms
        self._truncated_input_atoms = tuple(sorted(set(discarded)))

    @classmethod
    def _from_atoms(
        cls,
        atoms: Mapping[tuple[int, int], mp.mpf],
        pressure_order: int,
        width_order: int,
    ) -> "PressureWidthJet":
        result = cls.__new__(cls)
        result.pressure_order = pressure_order
        result.width_order = width_order
        result._truncated_input_atoms = ()
        result.atoms = {
            key: coefficient
            for key, coefficient in atoms.items()
            if not _is_zero(coefficient)
        }
        return result

    @property
    def orders(self) -> tuple[int, int]:
        return self.pressure_order, self.width_order

    @property
    def constant(self) -> mp.mpf:
        return self.component(0, 0)

    @property
    def truncated_input_atoms(self) -> tuple[tuple[int, int], ...]:
        return self._truncated_input_atoms

    def component(self, pressure_power: Any, width_power: Any = 0) -> mp.mpf:
        """Return the normalized coefficient of pressure**p * width**w."""

        pressure = _index(pressure_power, "pressure component index")
        width = _index(width_power, "width component index")
        if pressure < 0 or width < 0:
            return mp.mpf(0)
        if pressure > self.pressure_order or width > self.width_order:
            return mp.mpf(0)
        return self.atoms.get((pressure, width), mp.mpf(0))

    coefficient = component

    def evaluate(self, pressure: Any = 1, width: Any = 1) -> mp.mpf:
        """Evaluate the retained rectangle at the requested parameters."""

        pressure_value = _as_real_mp(pressure)
        width_value = _as_real_mp(width)
        return mp.fsum(
            coefficient
            * pressure_value**key[0]
            * width_value**key[1]
            for key, coefficient in self.atoms.items()
        )

    __call__ = evaluate

    def as_dict(self) -> dict[str, Any]:
        return {
            "pressure_order": self.pressure_order,
            "width_order": self.width_order,
            "atoms": {
                f"{pressure},{width}": mp.nstr(coefficient, 30)
                for (pressure, width), coefficient in sorted(self.atoms.items())
            },
            "rectangular_truncation": True,
            "normalized_coefficient_convention": NORMALIZED_BIVARIATE_CONVENTION,
            "truncated_input_atoms": [list(key) for key in self.truncated_input_atoms],
        }

    def __repr__(self) -> str:
        return (
            f"PressureWidthJet({self.as_dict()['atoms']!r}, "
            f"pressure_order={self.pressure_order}, width_order={self.width_order})"
        )

    def _coerce_operand(self, other: Any) -> "PressureWidthJet":
        if isinstance(other, PressureWidthJet):
            if other.orders != self.orders:
                raise ValueError(
                    "Bivariate Taylor orders must match: "
                    f"{self.orders} != {other.orders}"
                )
            return other
        return PressureWidthJet(
            other,
            pressure_order=self.pressure_order,
            width_order=self.width_order,
        )

    def _scale(self, scalar: Any) -> "PressureWidthJet":
        coefficient = _as_real_mp(scalar)
        return self._from_atoms(
            {key: value * coefficient for key, value in self.atoms.items()},
            self.pressure_order,
            self.width_order,
        )

    def __add__(self, other: Any) -> "PressureWidthJet":
        right = self._coerce_operand(other)
        terms: dict[tuple[int, int], list[mp.mpf]] = {}
        for key, value in self.atoms.items():
            terms.setdefault(key, []).append(value)
        for key, value in right.atoms.items():
            terms.setdefault(key, []).append(value)
        return self._from_atoms(
            {key: mp.fsum(values) for key, values in terms.items()},
            self.pressure_order,
            self.width_order,
        )

    __radd__ = __add__

    def __neg__(self) -> "PressureWidthJet":
        return self._from_atoms(
            {key: -value for key, value in self.atoms.items()},
            self.pressure_order,
            self.width_order,
        )

    def __sub__(self, other: Any) -> "PressureWidthJet":
        return self + (-self._coerce_operand(other))

    def __rsub__(self, other: Any) -> "PressureWidthJet":
        return self._coerce_operand(other) - self

    def __mul__(self, other: Any) -> "PressureWidthJet":
        right = self._coerce_operand(other)
        terms: dict[tuple[int, int], list[mp.mpf]] = {}
        for (pressure_a, width_a), value_a in self.atoms.items():
            for (pressure_b, width_b), value_b in right.atoms.items():
                key = (pressure_a + pressure_b, width_a + width_b)
                if (
                    key[0] <= self.pressure_order
                    and key[1] <= self.width_order
                ):
                    terms.setdefault(key, []).append(value_a * value_b)
        return self._from_atoms(
            {key: mp.fsum(values) for key, values in terms.items()},
            self.pressure_order,
            self.width_order,
        )

    __rmul__ = __mul__

    def reciprocal(self) -> "PressureWidthJet":
        """Return the rectangular finite-series reciprocal."""

        constant = self.constant
        if _is_zero(constant):
            raise ZeroDivisionError(
                "PressureWidthJet reciprocal requires a nonzero constant term"
            )
        nonconstant = self._scale(1 / constant) - 1
        result = PressureWidthJet(
            1,
            pressure_order=self.pressure_order,
            width_order=self.width_order,
        )
        term = result
        for power in range(1, self.pressure_order + self.width_order + 1):
            term = term * nonconstant
            if not term.atoms:
                break
            result = result + term._scale((-1) ** power)
        return result._scale(1 / constant)

    inverse = reciprocal
    inv = reciprocal

    def __truediv__(self, other: Any) -> "PressureWidthJet":
        return self * self._coerce_operand(other).reciprocal()

    def __rtruediv__(self, other: Any) -> "PressureWidthJet":
        return self._coerce_operand(other) * self.reciprocal()

    def __pow__(self, exponent: Any) -> "PressureWidthJet":
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
        result = PressureWidthJet(
            1,
            pressure_order=self.pressure_order,
            width_order=self.width_order,
        )
        base = self
        while integer:
            if integer & 1:
                result = result * base
            integer //= 2
            if integer:
                base = base * base
        return result

    def exp(self) -> "PressureWidthJet":
        """Return ``exp(self)`` through the retained rectangle."""

        constant = self.constant
        scale = mp.exp(constant)
        nonconstant = self - constant
        result = PressureWidthJet(
            scale,
            pressure_order=self.pressure_order,
            width_order=self.width_order,
        )
        term = PressureWidthJet(
            1,
            pressure_order=self.pressure_order,
            width_order=self.width_order,
        )
        factorial = mp.mpf(1)
        for power in range(1, self.pressure_order + self.width_order + 1):
            term = term * nonconstant
            if not term.atoms:
                break
            factorial *= power
            result = result + term._scale(scale / factorial)
        return result

    def log(self) -> "PressureWidthJet":
        """Return the real analytic logarithm, requiring a positive constant."""

        constant = self.constant
        if constant <= 0:
            raise ValueError(
                "PressureWidthJet log requires a positive constant term"
            )
        nonconstant = (self / constant) - 1
        result = PressureWidthJet(
            mp.log(constant),
            pressure_order=self.pressure_order,
            width_order=self.width_order,
        )
        term = PressureWidthJet(
            1,
            pressure_order=self.pressure_order,
            width_order=self.width_order,
        )
        for power in range(1, self.pressure_order + self.width_order + 1):
            term = term * nonconstant
            if not term.atoms:
                break
            result = result + term._scale(
                mp.mpf((-1) ** (power + 1)) / power
            )
        return result

    def sqrt(self) -> "PressureWidthJet":
        """Return the real analytic square root, requiring a positive constant."""

        constant = self.constant
        if constant <= 0:
            raise ValueError(
                "PressureWidthJet sqrt requires a positive constant term"
            )
        nonconstant = self / constant - 1
        root = mp.sqrt(constant)
        result = PressureWidthJet(
            root,
            pressure_order=self.pressure_order,
            width_order=self.width_order,
        )
        term = PressureWidthJet(
            1,
            pressure_order=self.pressure_order,
            width_order=self.width_order,
        )
        binomial = mp.mpf(1)
        half = mp.mpf("0.5")
        for power in range(1, self.pressure_order + self.width_order + 1):
            term = term * nonconstant
            if not term.atoms:
                break
            binomial *= (half - (power - 1)) / power
            result = result + term._scale(root * binomial)
        return result


def _nested_taylor(
    source: PressureWidthJet,
    function: Any,
) -> dict[tuple[int, int], mp.mpf]:
    """Reference bivariate coefficients using nested ``mp.taylor`` calls."""

    def resolved(pressure: mp.mpf, width: mp.mpf) -> mp.mpf:
        return function(source.evaluate(pressure, width))

    result: dict[tuple[int, int], mp.mpf] = {}
    for width_power in range(source.width_order + 1):
        outer = mp.taylor(
            lambda pressure: mp.taylor(
                lambda width: resolved(pressure, width), 0, source.width_order
            )[width_power],
            0,
            source.pressure_order,
        )
        for pressure_power in range(source.pressure_order + 1):
            result[(pressure_power, width_power)] = outer[pressure_power]
    return result


def _max_scaled_error(
    actual: PressureWidthJet,
    expected: Mapping[tuple[int, int], mp.mpf],
) -> tuple[mp.mpf, mp.mpf]:
    absolute: list[mp.mpf] = []
    scaled: list[mp.mpf] = []
    for pressure in range(actual.pressure_order + 1):
        for width in range(actual.width_order + 1):
            error = abs(actual.component(pressure, width) - expected[(pressure, width)])
            absolute.append(error)
            scaled.append(error / max(mp.mpf(1), abs(expected[(pressure, width)])))
    return max(absolute, default=mp.mpf(0)), max(scaled, default=mp.mpf(0))


def run_fixtures(
    *,
    pressure_order: int = DEFAULT_PRESSURE_ORDER,
    width_order: int = DEFAULT_WIDTH_ORDER,
    precision: int = 110,
) -> dict[str, Any]:
    """Run resolved-scale nonlinear checks and separated tiny-atom checks."""

    pressure_order = _validate_order(pressure_order, "pressure_order")
    width_order = _validate_order(width_order, "width_order")
    precision = max(50, int(precision))
    with mp.workdps(precision):
        source = PressureWidthJet(
            {
                (0, 0): mp.mpf("1.45"),
                (1, 0): mp.mpf("-0.21"),
                (0, 1): mp.mpf("0.17"),
                (2, 0): mp.mpf("0.035"),
                (1, 1): mp.mpf("-0.026"),
                (0, 2): mp.mpf("0.019"),
                (3, 1): mp.mpf("0.004"),
                (2, 2): mp.mpf("-0.002"),
            },
            pressure_order=pressure_order,
            width_order=width_order,
        )
        expected_inverse = _nested_taylor(source, lambda value: 1 / value)
        expected_log = _nested_taylor(source, mp.log)
        expected_exp = _nested_taylor(source, mp.exp)
        expected_sqrt = _nested_taylor(source, mp.sqrt)
        inverse_abs, inverse_scaled = _max_scaled_error(
            source.reciprocal(), expected_inverse
        )
        log_abs, log_scaled = _max_scaled_error(source.log(), expected_log)
        exp_abs, exp_scaled = _max_scaled_error(source.exp(), expected_exp)
        sqrt_abs, sqrt_scaled = _max_scaled_error(source.sqrt(), expected_sqrt)

        # Integer pressure keys model both legacy PressurePolynomial and the
        # one-variable PressureParameterJet atom interface.
        legacy = PressureWidthJet(
            {0: mp.mpf("1.25"), 1: mp.mpf("-0.4")},
            pressure_order=pressure_order,
            width_order=width_order,
        )
        legacy_mapping_ok = (
            legacy.component(0, 0) == mp.mpf("1.25")
            and legacy.component(1, 0) == mp.mpf("-0.4")
            and legacy.component(0, 1) == 0
        )

        tiny_pressure = mp.exp(-mp.mpf("1e28"))
        tiny_width = mp.exp(-mp.mpf("1e154"))
        mixed = mp.mpf(3) * tiny_pressure * tiny_width
        extreme = PressureWidthJet(
            {
                (0, 0): mp.mpf(1),
                (1, 0): tiny_pressure,
                (0, 1): tiny_width,
                (1, 1): mixed,
            },
            pressure_order=max(1, pressure_order),
            width_order=max(1, width_order),
        )
        extreme_exp = extreme.exp()
        expected_extreme_mixed_exp = mp.exp(1) * (
            mixed + tiny_pressure * tiny_width
        )
        tiny_checks = {
            "pressure_atom_preserved": extreme.component(1, 0) == tiny_pressure,
            "width_atom_preserved": extreme.component(0, 1) == tiny_width,
            "mixed_atom_preserved": extreme.component(1, 1) == mixed,
            "exp_pressure_atom_preserved": extreme_exp.component(1, 0)
            == mp.exp(1) * tiny_pressure,
            "exp_width_atom_preserved": extreme_exp.component(0, 1)
            == mp.exp(1) * tiny_width,
            "exp_mixed_atom_preserved": extreme_exp.component(1, 1)
            == expected_extreme_mixed_exp,
        }
        tolerance = mp.mpf(10) ** (-(precision // 2))
        checks = {
            "inverse": inverse_scaled <= tolerance,
            "log": log_scaled <= tolerance,
            "exp": exp_scaled <= tolerance,
            "sqrt": sqrt_scaled <= tolerance,
            "legacy_integer_key_mapping": legacy_mapping_ok,
            "tiny_atoms": all(tiny_checks.values()),
        }
        return {
            "passed": all(checks.values()),
            "checks": checks,
            "pressure_order": pressure_order,
            "width_order": width_order,
            "precision": precision,
            "rectangular_truncation": PRESSURE_WIDTH_RECTANGULAR_TRUNCATED,
            "pressure_parameter_order_truncated": PRESSURE_PARAMETER_ORDER_TRUNCATED,
            "exit_width_order_truncated": EXIT_WIDTH_ORDER_TRUNCATED,
            "normalized_coefficient_convention": NORMALIZED_BIVARIATE_CONVENTION,
            "inverse_max_absolute_error": mp.nstr(inverse_abs, 12),
            "inverse_max_scaled_error": mp.nstr(inverse_scaled, 12),
            "log_max_absolute_error": mp.nstr(log_abs, 12),
            "log_max_scaled_error": mp.nstr(log_scaled, 12),
            "exp_max_absolute_error": mp.nstr(exp_abs, 12),
            "exp_max_scaled_error": mp.nstr(exp_scaled, 12),
            "sqrt_max_absolute_error": mp.nstr(sqrt_abs, 12),
            "sqrt_max_scaled_error": mp.nstr(sqrt_scaled, 12),
            "tiny_pressure_atom": mp.nstr(tiny_pressure, 12),
            "tiny_width_atom": mp.nstr(tiny_width, 12),
            "tiny_checks": tiny_checks,
            "limitations": [
                "Products and nonlinear functions retain the declared rectangular powers only.",
                "Nominal evaluation at pressure=1 and width=1 can round away separated atoms.",
                "This module does not enclose the rectangular remainder or certify exact rational propagation.",
                "It does not install a field or claim prescribed bridge, global, or Navier--Stokes closure.",
            ],
        }


def run_fixture(
    *,
    pressure_order: int = DEFAULT_PRESSURE_ORDER,
    width_order: int = DEFAULT_WIDTH_ORDER,
    precision: int = 110,
) -> dict[str, Any]:
    return run_fixtures(
        pressure_order=pressure_order,
        width_order=width_order,
        precision=precision,
    )


def fixture(
    *,
    pressure_order: int = DEFAULT_PRESSURE_ORDER,
    width_order: int = DEFAULT_WIDTH_ORDER,
    precision: int = 110,
) -> dict[str, Any]:
    return run_fixture(
        pressure_order=pressure_order,
        width_order=width_order,
        precision=precision,
    )


def main() -> None:
    print(json.dumps(run_fixtures(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
