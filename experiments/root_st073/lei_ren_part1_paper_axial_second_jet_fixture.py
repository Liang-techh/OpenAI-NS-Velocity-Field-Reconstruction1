"""Independent checks for the second-order axial pressure/width jet.

The resolved checks compare the three stored slots with ``mp.diff`` applied to
the corresponding scalar rational, exponential, logarithmic, and square-root
functions.  A separate atom check keeps the pressure-tail and exit-width
coefficients distinct, including the very small mixed atom used by the
first-order jet fixture.  These are local algebra checks; they do not enclose
the rectangular jet remainder or install a global field.
"""

from __future__ import annotations

from pathlib import Path
import json
from typing import Any, Callable

import mpmath as mp

try:  # package import
    from .lei_ren_part1_paper_axial_dual import AxialDual
    from .lei_ren_part1_paper_axial_second_jet import AxialSecondJet
    from .lei_ren_part1_paper_pressure_width_jet import PressureWidthJet
except (ImportError, ValueError):  # direct experiment execution
    from lei_ren_part1_paper_axial_dual import AxialDual
    from lei_ren_part1_paper_axial_second_jet import AxialSecondJet
    from lei_ren_part1_paper_pressure_width_jet import PressureWidthJet


PRECISION = 160
PRESSURE_ORDER = 4
WIDTH_ORDER = 3
RESOLVED_Z = mp.mpf("0.37")
ERROR_THRESHOLD = mp.mpf("1e-120")


def _scaled(actual: Any, expected: Any) -> mp.mpf:
    actual_value = mp.mpf(actual)
    expected_value = mp.mpf(expected)
    return abs(actual_value - expected_value) / max(
        mp.mpf(1), abs(actual_value), abs(expected_value)
    )


def _oracle(
    function: Callable[[mp.mpf], mp.mpf], z: mp.mpf
) -> tuple[mp.mpf, mp.mpf, mp.mpf]:
    """Return ``f(z)``, ``f'(z)``, and ``f''(z)`` from independent ``mp.diff``."""

    return (
        function(z),
        mp.diff(function, z, 1),
        mp.diff(function, z, 2),
    )


def _base_scalar(z: mp.mpf) -> mp.mpf:
    return (
        mp.mpf("1.75")
        + mp.mpf("0.23") * z
        + mp.mpf("0.047") * z**2
        + mp.mpf("0.009") * z**3
    )


def _expression_functions() -> dict[str, Callable[[mp.mpf], mp.mpf]]:
    def rational(z: mp.mpf) -> mp.mpf:
        x = _base_scalar(z)
        return (x * x + 2) / (x + 3)

    def reciprocal(z: mp.mpf) -> mp.mpf:
        return 2 / _base_scalar(z)

    def power(z: mp.mpf) -> mp.mpf:
        return _base_scalar(z) ** -3

    def exponential(z: mp.mpf) -> mp.mpf:
        return mp.exp(_base_scalar(z))

    def logarithm(z: mp.mpf) -> mp.mpf:
        return mp.log(_base_scalar(z))

    def square_root(z: mp.mpf) -> mp.mpf:
        return mp.sqrt(_base_scalar(z))

    def composite(z: mp.mpf) -> mp.mpf:
        x = _base_scalar(z)
        rational_value = (x * x + 2) / (x + 3)
        return (rational_value + mp.exp(x) * mp.log(x)) / (mp.sqrt(x) + 3)

    return {
        "rational": rational,
        "reciprocal": reciprocal,
        "integer_power": power,
        "exp": exponential,
        "log": logarithm,
        "sqrt": square_root,
        "composite": composite,
    }


def _jet_expression(base: AxialSecondJet) -> dict[str, AxialSecondJet]:
    rational = (base * base + 2) / (base + 3)
    return {
        "rational": rational,
        "reciprocal": 2 / base,
        "integer_power": base**-3,
        "exp": base.exp(),
        "log": base.log(),
        "sqrt": base.sqrt(),
        "composite": (rational + base.exp() * base.log()) / (base.sqrt() + 3),
    }


def _atom_slot(
    constant: Any,
    pressure: Any,
    width: Any,
    mixed: Any,
) -> PressureWidthJet:
    return PressureWidthJet(
        {
            (0, 0): constant,
            (1, 0): pressure,
            (0, 1): width,
            (1, 1): mixed,
        },
        pressure_order=PRESSURE_ORDER,
        width_order=WIDTH_ORDER,
    )


def run_fixture(*, precision: int = PRECISION) -> dict[str, Any]:
    """Run the resolved analytic and separated pressure/width atom checks."""

    precision = max(100, int(precision))
    with mp.workdps(precision):
        # Parse the decimal at the fixture precision rather than retaining
        # the module-import precision of the convenience constant.
        z0 = mp.mpf("0.37")
        value = _base_scalar(z0)
        tangent = mp.diff(_base_scalar, z0, 1)
        second = mp.diff(_base_scalar, z0, 2)
        base = AxialSecondJet(
            value,
            tangent,
            second,
            pressure_order=PRESSURE_ORDER,
            width_order=WIDTH_ORDER,
        )

        functions = _expression_functions()
        jets = _jet_expression(base)
        oracle_errors: dict[str, dict[str, mp.mpf]] = {}
        for name, function in functions.items():
            expected = _oracle(function, z0)
            actual_jet = jets[name]
            actual = (
                actual_jet.value.evaluate(pressure=1, width=1),
                actual_jet.tangent.evaluate(pressure=1, width=1),
                actual_jet.second.evaluate(pressure=1, width=1),
            )
            oracle_errors[name] = {
                "value": _scaled(actual[0], expected[0]),
                "tangent": _scaled(actual[1], expected[1]),
                "second": _scaled(actual[2], expected[2]),
            }

        # Use separated arbitrary-exponent atoms.  Exact equality here is
        # intentional: the PW ring must retain nonzero atoms even when they
        # are far below the resolved scale.
        tiny_pressure = mp.exp(-mp.mpf("1e28"))
        tiny_width = mp.exp(-mp.mpf("1e154"))
        mixed = mp.mpf(3) * tiny_pressure * tiny_width
        atom_value = _atom_slot(2, tiny_pressure, tiny_width, mixed)
        atom_tangent = _atom_slot(
            5,
            2 * tiny_pressure,
            3 * tiny_width,
            7 * mixed,
        )
        atom_second = _atom_slot(
            11,
            5 * tiny_pressure,
            2 * tiny_width,
            13 * mixed,
        )
        atom_jet = AxialSecondJet(atom_value, atom_tangent, atom_second)
        atom_checks = {
            "value_pressure_atom": atom_jet.component(1, 0) == tiny_pressure,
            "value_width_atom": atom_jet.component(0, 1) == tiny_width,
            "value_mixed_atom": atom_jet.component(1, 1) == mixed,
            "tangent_pressure_atom": atom_jet.tangent_component(1, 0)
            == 2 * tiny_pressure,
            "tangent_width_atom": atom_jet.tangent_component(0, 1)
            == 3 * tiny_width,
            "tangent_mixed_atom": atom_jet.tangent_component(1, 1)
            == 7 * mixed,
            "second_pressure_atom": atom_jet.second_component(1, 0)
            == 5 * tiny_pressure,
            "second_width_atom": atom_jet.second_component(0, 1)
            == 2 * tiny_width,
            "second_mixed_atom": atom_jet.second_component(1, 1)
            == 13 * mixed,
        }

        # Check a non-linear atom coefficient independently from the class's
        # second-slot assembly.  For E=exp(value), E_pw is
        # exp(value_00) * (value_pw + value_p * value_w).
        atom_exp = atom_jet.exp()
        atom_exp_errors = {
            "value_pressure_atom": _scaled(
                atom_exp.component(1, 0), mp.exp(2) * tiny_pressure
            ),
            "value_width_atom": _scaled(
                atom_exp.component(0, 1), mp.exp(2) * tiny_width
            ),
            "value_mixed_atom": _scaled(
                atom_exp.component(1, 1),
                mp.exp(2) * (mixed + tiny_pressure * tiny_width),
            ),
            "tangent_pressure_atom": _scaled(
                atom_exp.tangent_component(1, 0),
                mp.exp(2) * (2 * tiny_pressure + 5 * tiny_pressure),
            ),
            "second_pressure_atom": _scaled(
                atom_exp.second_component(1, 0),
                mp.exp(2)
                * (5 * tiny_pressure + (2 * tiny_pressure) * (2 * tiny_pressure)),
            ),
        }

        # A first-order dual is deliberately rejected unless the missing
        # second slot is named explicitly.  The explicit promotion is checked
        # separately so a caller can bridge a known analytic second derivative.
        dual = AxialDual(atom_value, atom_tangent)
        unknown_second_rejected = False
        try:
            AxialSecondJet(dual)
        except TypeError:
            unknown_second_rejected = True
        explicit_promotion = AxialSecondJet(dual, second=atom_second)
        explicit_promotion_ok = (
            explicit_promotion.value.atoms == atom_value.atoms
            and explicit_promotion.tangent.atoms == atom_tangent.atoms
            and explicit_promotion.second.atoms == atom_second.atoms
        )
        dual_operand_rejected = False
        try:
            atom_jet + dual
        except TypeError:
            dual_operand_rejected = True

        checks = {
            "oracle": all(
                error <= ERROR_THRESHOLD
                for errors in oracle_errors.values()
                for error in errors.values()
            ),
            "atom_slots": all(atom_checks.values()),
            "atom_exp": all(error <= ERROR_THRESHOLD for error in atom_exp_errors.values()),
            "unknown_second_rejected": unknown_second_rejected,
            "explicit_dual_promotion": explicit_promotion_ok,
            "dual_operand_rejected": dual_operand_rejected,
            "orders": base.orders == (PRESSURE_ORDER, WIDTH_ORDER),
            "constant": base.constant == value,
        }
        maximum_oracle_error = max(
            (error for errors in oracle_errors.values() for error in errors.values()),
            default=mp.mpf(0),
        )
        maximum_atom_exp_error = max(atom_exp_errors.values(), default=mp.mpf(0))
        if not all(checks.values()):
            failed = [name for name, passed in checks.items() if not passed]
            raise AssertionError(f"AxialSecondJet fixture checks failed: {failed}")

        report = {
            "passed": True,
            "class": "AxialSecondJet",
            "precision": precision,
            "resolved_Z": mp.nstr(z0, 30),
            "pressure_order": PRESSURE_ORDER,
            "width_order": WIDTH_ORDER,
            "actual_derivative_slots": True,
            "automatic_first_Z_derivative": True,
            "automatic_second_Z_derivative": True,
            "checks": checks,
            "oracle": {
                name: {
                    key: mp.nstr(error, 18) for key, error in errors.items()
                }
                for name, errors in oracle_errors.items()
            },
            "maximum_scaled_oracle_error": mp.nstr(maximum_oracle_error, 30),
            "atom_checks": atom_checks,
            "atom_exp_errors": {
                key: mp.nstr(error, 18)
                for key, error in atom_exp_errors.items()
            },
            "maximum_scaled_atom_exp_error": mp.nstr(maximum_atom_exp_error, 30),
            "tiny_pressure_atom": mp.nstr(tiny_pressure, 18),
            "tiny_width_atom": mp.nstr(tiny_width, 18),
            "unsupported_calls": [
                "AxialDual cannot be promoted without an explicit second slot.",
                "AxialDual cannot be used as an arithmetic operand without an explicit second slot.",
                "Noninteger powers are unsupported; integer powers only.",
                "The rectangular pressure/width remainder is not enclosed.",
            ],
            "limitations": [
                "The three axial slots are local actual derivatives, not a global ODE enclosure.",
                "Pressure and width products are truncated to the declared rectangle.",
                "No field installation, bridge, stress-cone, or Navier--Stokes closure is claimed.",
            ],
        }
        Path(__file__).with_suffix(".json").write_text(
            json.dumps(report, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        return report


def fixture(**kwargs: Any) -> dict[str, Any]:
    return run_fixture(**kwargs)


def run_fixtures(**kwargs: Any) -> dict[str, Any]:
    """Plural compatibility alias used by the first-order jet fixture."""

    return run_fixture(**kwargs)


def main() -> None:
    print(json.dumps(run_fixture(), indent=2, sort_keys=True))


__all__ = ["fixture", "run_fixture", "run_fixtures"]


if __name__ == "__main__":
    main()
