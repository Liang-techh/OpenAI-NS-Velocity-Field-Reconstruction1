"""Synthetic exact checks for the correlated pulse density delta kernel.

This fixture is source-free.  It uses explicit Fraction-valued Z polynomials,
an independently expanded reference, and a separate mpmath path with a
1e-100-scale supplied difference jet.  It does not identify any physical
source, select a physical value, or establish equality of two source fields.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
import json

import mpmath as mp


from openai_ns_reconstruction.pulse_correlated_density_delta import (
    FieldJets,
    source_density_delta_and_Z1,
)


def _fraction(value: int | Fraction) -> Fraction:
    return value if isinstance(value, Fraction) else Fraction(value)


@dataclass(frozen=True)
class Polynomial:
    """Small exact Fraction polynomial helper for the independent reference."""

    coefficients: tuple[Fraction, ...]

    def __post_init__(self) -> None:
        values = tuple(_fraction(value) for value in self.coefficients)
        while len(values) > 1 and values[-1] == 0:
            values = values[:-1]
        object.__setattr__(self, "coefficients", values or (Fraction(0),))

    @classmethod
    def constant(cls, value: int | Fraction) -> "Polynomial":
        return cls((_fraction(value),))

    @classmethod
    def variable(cls) -> "Polynomial":
        return cls((Fraction(0), Fraction(1)))

    def __add__(self, other: "Polynomial | int | Fraction") -> "Polynomial":
        other_poly = _as_polynomial(other)
        width = max(len(self.coefficients), len(other_poly.coefficients))
        values = tuple(
            (self.coefficients[index] if index < len(self.coefficients) else 0)
            + (other_poly.coefficients[index] if index < len(other_poly.coefficients) else 0)
            for index in range(width)
        )
        return type(self)(values)

    def __radd__(self, other: "Polynomial | int | Fraction") -> "Polynomial":
        return self + other

    def __neg__(self) -> "Polynomial":
        return type(self)(tuple(-value for value in self.coefficients))

    def __sub__(self, other: "Polynomial | int | Fraction") -> "Polynomial":
        return self + (-_as_polynomial(other))

    def __rsub__(self, other: "Polynomial | int | Fraction") -> "Polynomial":
        return _as_polynomial(other) - self

    def __mul__(self, other: "Polynomial | int | Fraction") -> "Polynomial":
        other_poly = _as_polynomial(other)
        values = [Fraction(0) for _ in range(len(self.coefficients) + len(other_poly.coefficients) - 1)]
        for left_index, left in enumerate(self.coefficients):
            for right_index, right in enumerate(other_poly.coefficients):
                values[left_index + right_index] += left * right
        return type(self)(tuple(values))

    def __rmul__(self, other: "Polynomial | int | Fraction") -> "Polynomial":
        return self * other

    def derivative(self) -> "Polynomial":
        if len(self.coefficients) == 1:
            return type(self).constant(0)
        return type(self)(tuple(index * value for index, value in enumerate(self.coefficients[1:], 1)))

    def evaluate(self, value: Fraction) -> Fraction:
        result = Fraction(0)
        for coefficient in reversed(self.coefficients):
            result = result * value + coefficient
        return result


def _as_polynomial(value: Polynomial | int | Fraction) -> Polynomial:
    if isinstance(value, Polynomial):
        return value
    return Polynomial.constant(value)


_JET_NAMES = ("U", "U_Z", "Q", "Q_Z", "Q_ZZ", "Q_y", "Q_yZ")


def _add_jets(left: FieldJets, right: FieldJets) -> FieldJets:
    return FieldJets(**{name: getattr(left, name) + getattr(right, name) for name in _JET_NAMES})


def _negate_jets(value: FieldJets) -> FieldJets:
    return FieldJets(**{name: -getattr(value, name) for name in _JET_NAMES})


def _evaluate_polynomial_jets(polynomials: dict[str, Polynomial], Z: Fraction) -> FieldJets:
    return FieldJets(**{name: polynomials[name].evaluate(Z) for name in _JET_NAMES})


def _assert_coherent_polynomial_jets(polynomials: dict[str, Polynomial]) -> None:
    if polynomials["U"].derivative() != polynomials["U_Z"]:
        raise AssertionError("U_Z polynomial is not the Z derivative of U")
    if polynomials["Q"].derivative() != polynomials["Q_Z"]:
        raise AssertionError("Q_Z polynomial is not the Z derivative of Q")
    if polynomials["Q_Z"].derivative() != polynomials["Q_ZZ"]:
        raise AssertionError("Q_ZZ polynomial is not the Z derivative of Q_Z")
    if polynomials["Q_y"].derivative() != polynomials["Q_yZ"]:
        raise AssertionError("Q_yZ polynomial is not the Z derivative of Q_y")


def _direct_source_density_and_Z1(
    jets: FieldJets,
    *,
    Z: Fraction,
    delta: Fraction,
) -> tuple[Fraction, Fraction]:
    """Independent direct value/derivative expansion for one source side."""

    h = 1 - delta
    L = 1 - delta * Z * Z
    L_Z = -2 * delta * Z
    A = jets.Q + h * Z * jets.Q_Z / 2 + jets.Q_y
    A_Z = h * (jets.Q_Z + Z * jets.Q_ZZ) / 2 + jets.Q_Z + jets.Q_yZ
    B = jets.Q * (jets.Q / 2 + jets.Q_y)
    B_Z = jets.Q_Z * (jets.Q / 2 + jets.Q_y) + jets.Q * (
        jets.Q_Z / 2 + jets.Q_yZ
    )
    C = -2 * Z * jets.Q + (1 - Z * Z) * jets.Q_Z - 2 * Z * jets.Q_y
    C_Z = (
        -2 * jets.Q
        - 4 * Z * jets.Q_Z
        + (1 - Z * Z) * jets.Q_ZZ
        - 2 * jets.Q_y
        - 2 * Z * jets.Q_yZ
    )
    N = A + jets.U * C
    N_Z = A_Z + jets.U_Z * C + jets.U * C_Z
    S = N / L + B
    S_Z = N_Z / L + B_Z - N * L_Z / (L * L)
    return S, S_Z


def _direct_components(
    jets: FieldJets,
    *,
    Z: Fraction,
) -> dict[str, Fraction]:
    """Independent direct B and U*C component values and Z derivatives."""

    B = jets.Q * (jets.Q / 2 + jets.Q_y)
    B_Z = jets.Q_Z * (jets.Q / 2 + jets.Q_y) + jets.Q * (
        jets.Q_Z / 2 + jets.Q_yZ
    )
    C = -2 * Z * jets.Q + (1 - Z * Z) * jets.Q_Z - 2 * Z * jets.Q_y
    C_Z = (
        -2 * jets.Q
        - 4 * Z * jets.Q_Z
        + (1 - Z * Z) * jets.Q_ZZ
        - 2 * jets.Q_y
        - 2 * Z * jets.Q_yZ
    )
    return {
        "B": B,
        "B_Z": B_Z,
        "UC": jets.U * C,
        "UC_Z": jets.U_Z * C + jets.U * C_Z,
    }


def _symbolic_source_density_and_Z1(
    polynomials: dict[str, Polynomial],
    *,
    Z: Fraction,
    delta: Fraction,
) -> tuple[Fraction, Fraction]:
    """Exact polynomial/rational-function reference with quotient derivative."""

    half = (delta * 0 + 1) / (delta * 0 + 2)
    z = Polynomial.variable()
    h = 1 - delta
    L = 1 - delta * z * z
    Q = polynomials["Q"]
    Q_Z = polynomials["Q_Z"]
    Q_y = polynomials["Q_y"]
    U = polynomials["U"]
    A = Q + half * h * z * Q_Z + Q_y
    B = Q * (half * Q + Q_y)
    C = -2 * z * Q + (1 - z * z) * Q_Z - 2 * z * Q_y
    numerator = A + U * C
    numerator_Z = numerator.derivative()
    B_Z = B.derivative()
    L_Z = L.derivative()
    L_value = L.evaluate(Z)
    numerator_value = numerator.evaluate(Z)
    S = numerator_value / L_value + B.evaluate(Z)
    S_Z = (
        numerator_Z.evaluate(Z) / L_value
        - numerator_value * L_Z.evaluate(Z) / (L_value * L_value)
        + B_Z.evaluate(Z)
    )
    return S, S_Z


def _mpf_from_fraction(value: Fraction) -> mp.mpf:
    return mp.mpf(value.numerator) / value.denominator


def _mpf_jets(value: FieldJets) -> FieldJets:
    return FieldJets(**{name: _mpf_from_fraction(getattr(value, name)) for name in _JET_NAMES})


def _mpf_scale_jets(value: FieldJets, scale: mp.mpf) -> FieldJets:
    return FieldJets(**{name: getattr(value, name) * scale for name in _JET_NAMES})


def _make_frame_polynomials() -> dict[str, Polynomial]:
    return {
        "U": Polynomial((Fraction(1, 3), Fraction(-2, 5), Fraction(1, 7))),
        "U_Z": Polynomial((Fraction(-2, 5), Fraction(2, 7))),
        "Q": Polynomial((Fraction(2, 5), Fraction(1, 6), Fraction(-3, 10), Fraction(1, 11))),
        "Q_Z": Polynomial((Fraction(1, 6), Fraction(-3, 5), Fraction(3, 11))),
        "Q_ZZ": Polynomial((Fraction(-3, 5), Fraction(6, 11))),
        "Q_y": Polynomial((Fraction(-1, 7), Fraction(2, 9), Fraction(1, 13))),
        "Q_yZ": Polynomial((Fraction(2, 9), Fraction(2, 13))),
    }


def _make_delta_polynomials() -> dict[str, Polynomial]:
    epsilon = Fraction(1, 10**100)
    return {
        "U": epsilon * Polynomial((Fraction(1), Fraction(2), Fraction(-1))),
        "U_Z": epsilon * Polynomial((Fraction(2), Fraction(-2))),
        "Q": epsilon * Polynomial((Fraction(2), Fraction(-1), Fraction(3))),
        "Q_Z": epsilon * Polynomial((Fraction(-1), Fraction(6))),
        "Q_ZZ": epsilon * Polynomial((Fraction(6),)),
        "Q_y": epsilon * Polynomial((Fraction(1, 3), Fraction(4))),
        "Q_yZ": epsilon * Polynomial((Fraction(4),)),
    }


def run_fixture() -> dict[str, object]:
    """Return one JSON-compatible synthetic result without source claims."""

    frame_polynomials = _make_frame_polynomials()
    delta_polynomials = _make_delta_polynomials()
    _assert_coherent_polynomial_jets(frame_polynomials)
    _assert_coherent_polynomial_jets(delta_polynomials)

    Z = Fraction(2, 7)
    delta_parameter = Fraction(3, 10)
    frame = _evaluate_polynomial_jets(frame_polynomials, Z)
    difference = _evaluate_polynomial_jets(delta_polynomials, Z)
    physical = _add_jets(frame, difference)

    result = source_density_delta_and_Z1(
        frame,
        difference,
        Z=Z,
        delta=delta_parameter,
        delta_z_independent=True,
    )

    frame_direct = _direct_source_density_and_Z1(frame, Z=Z, delta=delta_parameter)
    physical_direct = _direct_source_density_and_Z1(physical, Z=Z, delta=delta_parameter)
    direct_difference = (
        physical_direct[0] - frame_direct[0],
        physical_direct[1] - frame_direct[1],
    )
    frame_symbolic = _symbolic_source_density_and_Z1(
        frame_polynomials,
        Z=Z,
        delta=delta_parameter,
    )
    physical_symbolic = _symbolic_source_density_and_Z1(
        {name: frame_polynomials[name] + delta_polynomials[name] for name in _JET_NAMES},
        Z=Z,
        delta=delta_parameter,
    )
    symbolic_difference = (
        physical_symbolic[0] - frame_symbolic[0],
        physical_symbolic[1] - frame_symbolic[1],
    )

    if (result.delta_S, result.delta_S_Z) != direct_difference:
        raise AssertionError("correlated delta disagrees with independent direct expansion")
    if direct_difference != symbolic_difference:
        raise AssertionError("direct expansion disagrees with exact polynomial derivative")

    frame_components = _direct_components(frame, Z=Z)
    physical_components = _direct_components(physical, Z=Z)
    if result.delta_B != physical_components["B"] - frame_components["B"]:
        raise AssertionError("returned Delta B disagrees with direct components")
    if result.delta_B_Z != physical_components["B_Z"] - frame_components["B_Z"]:
        raise AssertionError("returned Delta B_Z disagrees with direct components")
    if result.delta_UC != physical_components["UC"] - frame_components["UC"]:
        raise AssertionError("returned Delta(U*C) disagrees with direct components")
    if result.delta_UC_Z != physical_components["UC_Z"] - frame_components["UC_Z"]:
        raise AssertionError("returned Delta(U*C)_Z disagrees with direct components")

    swapped = source_density_delta_and_Z1(
        physical,
        _negate_jets(difference),
        Z=Z,
        delta=delta_parameter,
        delta_z_independent=True,
    )
    if (swapped.delta_S, swapped.delta_S_Z) != (-result.delta_S, -result.delta_S_Z):
        raise AssertionError("swapped physical/frame difference did not reverse sign")

    zero = FieldJets(*(Fraction(0) for _ in _JET_NAMES))
    zero_result = source_density_delta_and_Z1(
        frame,
        zero,
        Z=Z,
        delta=delta_parameter,
        delta_z_independent=True,
    )
    if zero_result.delta_S != 0 or zero_result.delta_S_Z != 0:
        raise AssertionError("zero supplied difference was not retained as zero")

    try:
        source_density_delta_and_Z1(
            frame,
            difference,
            Z=Z,
            delta=delta_parameter,
            delta_z_independent=False,
        )
    except ValueError:
        delta_independence_rejected = True
    else:
        raise AssertionError("delta_z_independent=False was accepted")

    try:
        source_density_delta_and_Z1(
            frame,
            zero,
            Z=Fraction(1),
            delta=Fraction(1),
            delta_z_independent=True,
        )
    except ValueError:
        singular_L_rejected = True
    else:
        raise AssertionError("L == 0 was accepted")

    with mp.workdps(240):
        frame_mpf = _mpf_jets(frame)
        difference_mpf = _mpf_scale_jets(
            _mpf_jets(_evaluate_polynomial_jets(
                {name: delta_polynomials[name] * 10**100 for name in _JET_NAMES},
                Z,
            )),
            mp.mpf("1e-100"),
        )
        Z_mpf = mp.mpf(2) / 7
        delta_parameter_mpf = mp.mpf(3) / 10
        physical_mpf = _add_jets(frame_mpf, difference_mpf)
        frame_reference = _direct_source_density_and_Z1(
            frame_mpf,
            Z=Z_mpf,
            delta=delta_parameter_mpf,
        )
        physical_reference = _direct_source_density_and_Z1(
            physical_mpf,
            Z=Z_mpf,
            delta=delta_parameter_mpf,
        )
        mp_reference_difference = (
            physical_reference[0] - frame_reference[0],
            physical_reference[1] - frame_reference[1],
        )
        with mp.workdps(180):
            mp_result = source_density_delta_and_Z1(
                frame_mpf,
                difference_mpf,
                Z=Z_mpf,
                delta=delta_parameter_mpf,
                delta_z_independent=True,
            )
        if mp_result.delta_S == 0 or mp_result.delta_S_Z == 0:
            raise AssertionError("mpmath 1e-100 supplied difference was lost")
        if not isinstance(mp_result.delta_S, type(mp.mpf("0"))):
            raise AssertionError("mpmath numeric type was not preserved")
        mp_reference_error = max(
            abs(mp_result.delta_S - mp_reference_difference[0]),
            abs(mp_result.delta_S_Z - mp_reference_difference[1]),
        )
        if mp_reference_error > mp.mpf("1e-170"):
            raise AssertionError("mpmath correlated result disagrees with high-precision direct reference")
        mp_small_difference = str(mp_result.delta_S)
        mp_reference_error_text = str(mp_reference_error)

    controls = {
        "fraction_exact_direct_difference": True,
        "fraction_exact_symbolic_Z1": True,
        "direct_delta_B_components": True,
        "direct_delta_B_Z_components": True,
        "direct_delta_UC_components": True,
        "direct_delta_UC_Z_components": True,
        "tiny_10^-100_difference_preserved": result.delta_S != 0 and result.delta_S_Z != 0,
        "mpmath_numeric_type_preserved": True,
        "mpmath_high_precision_direct_reference": True,
        "swapped_pair_reverses_sign": True,
        "zero_difference_is_zero": True,
        "delta_Z_independence_required": delta_independence_rejected,
        "singular_L_rejected": singular_L_rejected,
    }
    return {
        "schema": "pulse_correlated_density_delta_synthetic_fixture_v1",
        "synthetic_only": True,
        "source_values_selected": False,
        "physical_frame_equality_proved": False,
        "accepted_as_source_binding": False,
        "Z": str(Z),
        "delta_parameter": str(delta_parameter),
        "tiny_difference_scale": "1/10^100",
        "delta_S": str(result.delta_S),
        "delta_S_Z1": str(result.delta_S_Z),
        "mpmath_delta_S": mp_small_difference,
        "mpmath_reference_error": mp_reference_error_text,
        "controls": controls,
        "passed": all(controls.values()),
    }


if __name__ == "__main__":
    print(json.dumps(run_fixture(), indent=2, sort_keys=True))
