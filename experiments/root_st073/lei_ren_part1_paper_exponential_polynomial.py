"""Coefficient-valued exponential polynomials in the axial coordinate.

An :class:`ExponentialPolynomial` stores terms

``c[a, b] * exp(a*y) * y**b``

with integer exponential exponent ``a`` and nonnegative polynomial power
``b``.  Coefficients stay in their supplied ring: they may be MP scalars,
``PressureWidthJet`` instances, or ``AxialDual`` instances.  In particular,
``evaluate`` and ``integral_zero_to`` can use a pressure-width jet or an
axial dual as ``y`` without first materialising its atoms at one nominal
point.

Integration is exact term by term.  For nonzero integer ``a`` it uses the
finite integration-by-parts recurrence; for ``a = 0`` it raises the ordinary
monomial degree.  ``antiderivative_zero_at_0`` returns the same algebra type,
with the lower-endpoint constant included explicitly.
"""

from __future__ import annotations

from collections.abc import Mapping
import json
import operator
from typing import Any

import mpmath as mp

try:
    from .lei_ren_part1_paper_pressure_width_jet import PressureWidthJet
    from .lei_ren_part1_paper_axial_dual import AxialDual
except (ImportError, ValueError):
    from lei_ren_part1_paper_pressure_width_jet import PressureWidthJet
    from lei_ren_part1_paper_axial_dual import AxialDual


def _integer(value: Any, name: str) -> int:
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
    if isinstance(raw_key, tuple):
        if len(raw_key) != 2:
            raise TypeError("Exponential-polynomial keys must be (a, b)")
        exponent = _integer(raw_key[0], "exponential exponent")
        power = _integer(raw_key[1], "polynomial power")
    else:
        # A bare integer is a useful shorthand for exp(a*y) * y**0.
        exponent = _integer(raw_key, "exponential exponent")
        power = 0
    if power < 0:
        raise ValueError("polynomial powers must be nonnegative")
    return exponent, power


def _looks_like_coefficient_ring(value: Any) -> bool:
    if isinstance(value, (PressureWidthJet, AxialDual, mp.mpf, mp.mpc)):
        return True
    # Keep compatible ring objects opaque if this module is imported through
    # a different module path than the experiment's canonical one.
    return hasattr(value, "component") and (
        hasattr(value, "evaluate") or hasattr(value, "value")
    )


def _coefficient(value: Any) -> Any:
    """Convert only scalar constants; preserve coefficient-valued rings."""

    if _looks_like_coefficient_ring(value):
        return value
    try:
        return mp.mpf(value)
    except (TypeError, ValueError):
        # A caller-supplied coefficient ring can remain opaque when it
        # supports the arithmetic protocol used below.
        return value


def _add(left: Any, right: Any) -> Any:
    try:
        return left + right
    except TypeError as first_error:
        try:
            return right + left
        except TypeError:
            raise first_error


def _mul(left: Any, right: Any) -> Any:
    try:
        return left * right
    except TypeError as first_error:
        try:
            return right * left
        except TypeError:
            raise first_error


def _scale(value: Any, scalar: Any) -> Any:
    return _mul(value, mp.mpf(scalar))


def _typed_zero(value: Any) -> Any:
    """Derive a zero from a coefficient rather than generic ``0`` coercion."""

    return _mul(value, mp.mpf(0))


def _prepare_y(value: Any) -> Any:
    if isinstance(value, (PressureWidthJet, AxialDual, mp.mpf, mp.mpc)):
        return value
    if hasattr(value, "exp") and hasattr(value, "__pow__"):
        return value
    return mp.mpf(value)


def _exp_at(exponent: int, y: Any) -> Any:
    if exponent == 0:
        return mp.mpf(1)
    if isinstance(y, (mp.mpf, mp.mpc)):
        return mp.exp(mp.mpf(exponent) * y)
    argument = y if exponent == 1 else y * exponent
    method = getattr(argument, "exp", None)
    if callable(method):
        return method()
    return mp.exp(argument)


def _power_at(y: Any, power: int) -> Any:
    if power == 0:
        return mp.mpf(1)
    if hasattr(y, "__pow__"):
        return y**power
    return mp.mpf(y) ** power


class ExponentialPolynomial:
    """Sparse coefficient-valued sum of ``exp(a*y) * y**b`` terms."""

    def __init__(self, value: Any = 0):
        if isinstance(value, ExponentialPolynomial):
            source: Any = value.atoms
        elif hasattr(value, "atoms") and isinstance(value.atoms, Mapping):
            source = value.atoms
        elif isinstance(value, Mapping):
            source = value
        else:
            source = {(0, 0): value}

        atoms: dict[tuple[int, int], Any] = {}
        for raw_key, raw_value in source.items():
            key = _key(raw_key)
            coefficient = _coefficient(raw_value)
            if key in atoms:
                atoms[key] = _add(atoms[key], coefficient)
            else:
                atoms[key] = coefficient
        self.atoms = atoms

    @property
    def terms(self) -> dict[tuple[int, int], Any]:
        return self.atoms

    def coefficient(self, exponent: Any, power: Any = 0) -> Any:
        return self.atoms.get(
            (_integer(exponent, "exponential exponent"), _integer(power, "polynomial power")),
            mp.mpf(0),
        )

    def __repr__(self) -> str:
        return f"ExponentialPolynomial({self.atoms!r})"

    def _coerce(self, other: Any) -> "ExponentialPolynomial":
        if isinstance(other, ExponentialPolynomial):
            return other
        return ExponentialPolynomial({(0, 0): _coefficient(other)})

    def __add__(self, other: Any) -> "ExponentialPolynomial":
        right = self._coerce(other)
        atoms = dict(self.atoms)
        for key, value in right.atoms.items():
            if key in atoms:
                atoms[key] = _add(atoms[key], value)
            else:
                atoms[key] = value
        return ExponentialPolynomial(atoms)

    __radd__ = __add__

    def __neg__(self) -> "ExponentialPolynomial":
        return ExponentialPolynomial({key: -value for key, value in self.atoms.items()})

    def __sub__(self, other: Any) -> "ExponentialPolynomial":
        return self + (-self._coerce(other))

    def __rsub__(self, other: Any) -> "ExponentialPolynomial":
        return self._coerce(other) - self

    def __mul__(self, other: Any) -> "ExponentialPolynomial":
        right = self._coerce(other)
        atoms: dict[tuple[int, int], Any] = {}
        for (a_left, b_left), left in self.atoms.items():
            for (a_right, b_right), right_value in right.atoms.items():
                key = (a_left + a_right, b_left + b_right)
                value = _mul(left, right_value)
                if key in atoms:
                    atoms[key] = _add(atoms[key], value)
                else:
                    atoms[key] = value
        if not atoms and (self.atoms or right.atoms):
            sample = next(iter(self.atoms.values()), None)
            if sample is None:
                sample = next(iter(right.atoms.values()))
            atoms[(0, 0)] = _typed_zero(sample)
        return ExponentialPolynomial(atoms)

    __rmul__ = __mul__

    def __truediv__(self, other: Any) -> "ExponentialPolynomial":
        # Division is intentionally limited to a scalar coefficient.  A
        # general exponential-polynomial quotient is outside this algebra.
        if isinstance(other, ExponentialPolynomial):
            raise TypeError("exponential-polynomial division is not supported")
        denominator = _coefficient(other)
        return ExponentialPolynomial({key: value / denominator for key, value in self.atoms.items()})

    def evaluate(self, y: Any) -> Any:
        """Evaluate without materialising coefficient or y atoms."""

        argument = _prepare_y(y)
        total: Any = None
        for (exponent, power), coefficient in self.atoms.items():
            term = _mul(coefficient, _exp_at(exponent, argument))
            term = _mul(term, _power_at(argument, power))
            total = term if total is None else _add(total, term)
        if total is None:
            return mp.mpf(0)
        return total

    __call__ = evaluate

    def derivative(self) -> "ExponentialPolynomial":
        """Differentiate term by term with respect to the polynomial argument."""

        atoms: dict[tuple[int, int], Any] = {}
        for (exponent, power), coefficient in self.atoms.items():
            if exponent != 0:
                key = (exponent, power)
                value = _scale(coefficient, exponent)
                atoms[key] = _add(atoms[key], value) if key in atoms else value
            if power > 0:
                key = (exponent, power - 1)
                value = _scale(coefficient, power)
                atoms[key] = _add(atoms[key], value) if key in atoms else value
        if not atoms and self.atoms:
            sample = next(iter(self.atoms.values()))
            atoms[(0, 0)] = _typed_zero(sample)
        return ExponentialPolynomial(atoms)

    differentiate = derivative

    def antiderivative_zero_at_0(self) -> "ExponentialPolynomial":
        """Return the exact antiderivative whose value at ``y=0`` is zero."""

        atoms: dict[tuple[int, int], Any] = {}
        for (exponent, power), coefficient in self.atoms.items():
            if exponent == 0:
                key = (0, power + 1)
                value = _scale(coefficient, mp.mpf(1) / (power + 1))
                atoms[key] = _add(atoms[key], value) if key in atoms else value
                continue

            denominator = mp.mpf(exponent)
            factorial = mp.factorial(power)
            for index in range(power + 1):
                lower_power = power - index
                scalar = (
                    (-1) ** index
                    * factorial
                    / mp.factorial(lower_power)
                    / denominator ** (index + 1)
                )
                key = (exponent, lower_power)
                value = _scale(coefficient, scalar)
                atoms[key] = _add(atoms[key], value) if key in atoms else value

            # The exponential polynomial at y=0 contributes only its
            # zero-power term, including the index=power endpoint above.
            endpoint_scalar = (
                (-1) ** power * factorial / denominator ** (power + 1)
            )
            key = (0, 0)
            value = _scale(coefficient, -endpoint_scalar)
            atoms[key] = _add(atoms[key], value) if key in atoms else value
        return ExponentialPolynomial(atoms)

    antiderivative_zero_at0 = antiderivative_zero_at_0
    integral_polynomial = antiderivative_zero_at_0

    def integral_zero_to(self, y: Any) -> Any:
        """Evaluate the exact zero-based antiderivative at ``y``."""

        return self.antiderivative_zero_at_0().evaluate(y)


ExponentialPolynomialAlgebra = ExponentialPolynomial
CoefficientExponentialPolynomial = ExponentialPolynomial


def _scalar_fixture() -> dict[str, Any]:
    modes = {-2: 1, -1: 2, 0: 3, 1: 4, 2: 5}
    terms = {
        (exponent, power): mp.mpf(magnitude) / (power + 2)
        for exponent, magnitude in modes.items()
        for power in range(4)
    }
    polynomial = ExponentialPolynomial(terms)
    upper = mp.mpf(".37")
    lower = mp.mpf("-.41")
    analytic_upper = polynomial.integral_zero_to(upper)
    analytic_lower = polynomial.integral_zero_to(lower)
    numeric_upper = mp.quad(lambda value: polynomial.evaluate(value), [0, upper])
    numeric_lower = mp.quad(lambda value: polynomial.evaluate(value), [0, lower])
    integral_error = max(
        abs(analytic_upper - numeric_upper), abs(analytic_lower - numeric_lower)
    )
    endpoint = polynomial.antiderivative_zero_at_0().evaluate(mp.mpf(0))

    other = ExponentialPolynomial(
        {
            (exponent, power): mp.mpf(power + 1) / (abs(exponent) + 3)
            for exponent in (-2, -1, 0, 1, 2)
            for power in range(4)
        }
    )
    product = polynomial * other
    point = mp.mpf(".23")
    multiplication_error = abs(
        product.evaluate(point)
        - polynomial.evaluate(point) * other.evaluate(point)
    )
    derivative_error = abs(
        product.derivative().evaluate(point)
        - mp.diff(lambda value: product.evaluate(value), point)
    )
    return {
        "integral_error": integral_error,
        "endpoint_zero": endpoint,
        "multiplication_error": multiplication_error,
        "derivative_error": derivative_error,
        "mode_exponents": sorted(modes),
        "polynomial_powers": list(range(4)),
    }


def _tiny_atom_fixture() -> dict[str, Any]:
    tiny_pressure = mp.exp(-mp.mpf("1e28"))
    tiny_width = mp.exp(-mp.mpf("1e154"))
    mixed = mp.mpf(3) * tiny_pressure * tiny_width
    base = PressureWidthJet(
        {
            (0, 0): mp.mpf(1),
            (1, 0): tiny_pressure,
            (0, 1): tiny_width,
            (1, 1): mixed,
        },
        pressure_order=3,
        width_order=2,
    )
    tangent = PressureWidthJet(
        {(1, 0): 2 * tiny_pressure, (0, 1): 3 * tiny_width, (1, 1): 5 * mixed},
        pressure_order=3,
        width_order=2,
    )
    dual = AxialDual(base, tangent)
    polynomial = ExponentialPolynomial({(0, 0): dual})
    value = polynomial.evaluate(mp.mpf(1))
    integral = polynomial.integral_zero_to(mp.mpf(1))
    checks = {
        "pressure_value_atom": value.value.component(1, 0) == tiny_pressure,
        "width_value_atom": value.value.component(0, 1) == tiny_width,
        "mixed_value_atom": value.value.component(1, 1) == mixed,
        "pressure_tangent_atom": value.tangent.component(1, 0) == 2 * tiny_pressure,
        "width_tangent_atom": value.tangent.component(0, 1) == 3 * tiny_width,
        "mixed_tangent_atom": value.tangent.component(1, 1) == 5 * mixed,
        "integral_pressure_atom": integral.value.component(1, 0) == tiny_pressure,
        "integral_width_atom": integral.value.component(0, 1) == tiny_width,
        "integral_mixed_atom": integral.value.component(1, 1) == mixed,
    }
    return {
        "checks": checks,
        "passed": all(checks.values()),
        "tiny_pressure": mp.nstr(tiny_pressure, 12),
        "tiny_width": mp.nstr(tiny_width, 12),
    }


def run_fixture(*, precision: int = 100) -> dict[str, Any]:
    with mp.workdps(max(70, int(precision))):
        scalar = _scalar_fixture()
        tiny = _tiny_atom_fixture()
        tolerance = mp.mpf(10) ** (-(max(70, int(precision)) // 2))
        checks = {
            "quadrature": scalar["integral_error"] <= tolerance,
            "endpoint_zero": abs(scalar["endpoint_zero"]) <= tolerance,
            "multiplication": scalar["multiplication_error"] <= tolerance,
            "derivative": scalar["derivative_error"] <= tolerance,
            "tiny_atoms": tiny["passed"],
        }
        return {
            "passed": all(checks.values()),
            "checks": checks,
            "quadrature_max_absolute_error": mp.nstr(scalar["integral_error"], 18),
            "endpoint_zero_error": mp.nstr(abs(scalar["endpoint_zero"]), 18),
            "multiplication_error": mp.nstr(scalar["multiplication_error"], 18),
            "derivative_error": mp.nstr(scalar["derivative_error"], 18),
            "mode_exponents": scalar["mode_exponents"],
            "polynomial_powers": scalar["polynomial_powers"],
            "tiny": tiny,
            "coefficient_rings": ["mp.mpf", "PressureWidthJet", "AxialDual"],
            "lower_endpoint_constant_included": True,
            "limitations": [
                "The algebra is finite and coefficient-wise; it does not provide a generic exponential of a polynomial exponent.",
                "Quadrature checks are resolved-scale fixtures and are not error enclosures for production continuation.",
            ],
        }


def fixture(**kwargs: Any) -> dict[str, Any]:
    return run_fixture(**kwargs)


def main() -> None:
    print(json.dumps(run_fixture(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
