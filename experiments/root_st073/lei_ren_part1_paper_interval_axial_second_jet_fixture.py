"""Bounded checks for the directed interval axial second jet.

The fixture checks every retained pressure/width atom in each of the three
actual ``Z`` derivative slots.  Nominal coefficients are compared with the
legacy :class:`AxialSecondJet`, while independent coefficient perturbations
are propagated through both implementations and checked against the
``IntervalDifference.difference`` intervals.  A second pass widens nominal
coefficient intervals so the checks exercise enclosure rather than only
point nominal arithmetic.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Callable

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext

try:  # package import
    from .lei_ren_part1_paper_axial_dual import AxialDual
    from .lei_ren_part1_paper_axial_second_jet import AxialSecondJet
    from .lei_ren_part1_paper_interval_axial_second_jet import (
        IntervalAxialSecondJet,
    )
    from .lei_ren_part1_paper_interval_difference import IntervalDifference
    from .lei_ren_part1_paper_interval_pressure_width_jet import (
        IntervalPressureWidthJet,
    )
    from .lei_ren_part1_paper_pressure_width_jet import PressureWidthJet
except (ImportError, ValueError):  # direct experiment execution
    from lei_ren_part1_paper_axial_dual import AxialDual
    from lei_ren_part1_paper_axial_second_jet import AxialSecondJet
    from lei_ren_part1_paper_interval_axial_second_jet import (
        IntervalAxialSecondJet,
    )
    from lei_ren_part1_paper_interval_difference import IntervalDifference
    from lei_ren_part1_paper_interval_pressure_width_jet import (
        IntervalPressureWidthJet,
    )
    from lei_ren_part1_paper_pressure_width_jet import PressureWidthJet


PRECISION = 140
PRESSURE_ORDER = 2
WIDTH_ORDER = 2
ERROR_THRESHOLD = mp.mpf("1e-100")


def _bounds(value: Any) -> tuple[mp.mpf, mp.mpf]:
    return tuple(mp.make_mpf(item) for item in value._mpi_)


def _inside(value: Any, interval: Any) -> bool:
    lower, upper = _bounds(interval)
    return lower <= value <= upper


def _exact_zero(value: Any) -> bool:
    lower, upper = _bounds(value)
    return lower == 0 and upper == 0


def _scaled(actual: Any, expected: Any) -> mp.mpf:
    actual_value = mp.mpf(actual)
    expected_value = mp.mpf(expected)
    return abs(actual_value - expected_value) / max(
        mp.mpf(1), abs(actual_value), abs(expected_value)
    )


def _str(value: Any) -> str:
    return mp.nstr(mp.mpf(value), 100)


def _base_scalar(z: mp.mpf) -> mp.mpf:
    return (
        mp.mpf("2.15")
        + mp.mpf("0.27") * z
        + mp.mpf("0.041") * z**2
        + mp.mpf("0.008") * z**3
    )


def _scalar_functions() -> dict[str, Callable[[mp.mpf], mp.mpf]]:
    def rational(z: mp.mpf) -> mp.mpf:
        x = _base_scalar(z)
        return (x * x + 2) / (x + 3)

    def subtract(z: mp.mpf) -> mp.mpf:
        return 3 - _base_scalar(z)

    def reciprocal(z: mp.mpf) -> mp.mpf:
        return 2 / _base_scalar(z)

    def positive_power(z: mp.mpf) -> mp.mpf:
        return _base_scalar(z) ** 3

    def negative_power(z: mp.mpf) -> mp.mpf:
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
        "identity": _base_scalar,
        "add": lambda z: _base_scalar(z) + 2,
        "subtract": subtract,
        "multiply": lambda z: _base_scalar(z) * (_base_scalar(z) + 2),
        "rational": rational,
        "reciprocal": reciprocal,
        "positive_power": positive_power,
        "negative_power": negative_power,
        "exp": exponential,
        "log": logarithm,
        "sqrt": square_root,
        "composite": composite,
    }


def _jet_functions(base: Any) -> dict[str, Any]:
    rational = (base * base + 2) / (base + 3)
    return {
        "identity": base,
        "add": base + 2,
        "subtract": 3 - base,
        "multiply": base * (base + 2),
        "rational": rational,
        "reciprocal": 2 / base,
        "positive_power": base**3,
        "negative_power": base**-3,
        "exp": base.exp(),
        "log": base.log(),
        "sqrt": base.sqrt(),
        "composite": (rational + base.exp() * base.log()) / (base.sqrt() + 3),
    }


def _legacy_jet(
    value: dict[tuple[int, int], Any],
    tangent: dict[tuple[int, int], Any],
    second: dict[tuple[int, int], Any],
) -> AxialSecondJet:
    orders = {"pressure_order": PRESSURE_ORDER, "width_order": WIDTH_ORDER}
    return AxialSecondJet(
        PressureWidthJet(value, **orders),
        PressureWidthJet(tangent, **orders),
        PressureWidthJet(second, **orders),
    )


def _interval_slot(
    ctx: MPIntervalContext,
    mapping: dict[tuple[int, int], Any],
    errors: dict[tuple[int, int], mp.mpf],
    *,
    wide: bool,
    zero_difference: bool = False,
) -> IntervalPressureWidthJet:
    atoms: dict[tuple[int, int], IntervalDifference] = {}
    for key, value in mapping.items():
        center = mp.mpf(value)
        if wide:
            half_width = mp.mpf("2e-9") * max(mp.mpf(1), abs(center))
            nominal = ctx.mpf([_str(center - half_width), _str(center + half_width)])
        else:
            nominal = ctx.mpf(_str(center))
        if zero_difference:
            difference = ctx.mpf(0)
        else:
            error = errors[key]
            difference = ctx.mpf([_str(-error), _str(error)])
        atoms[key] = IntervalDifference(ctx, nominal, difference)
    return IntervalPressureWidthJet(
        ctx,
        atoms,
        pressure_order=PRESSURE_ORDER,
        width_order=WIDTH_ORDER,
    )


def _interval_jet(
    ctx: MPIntervalContext,
    value: dict[tuple[int, int], Any],
    tangent: dict[tuple[int, int], Any],
    second: dict[tuple[int, int], Any],
    errors: dict[tuple[int, int], mp.mpf],
    *,
    wide: bool,
    zero_difference: bool = False,
) -> IntervalAxialSecondJet:
    return IntervalAxialSecondJet(
        ctx,
        _interval_slot(ctx, value, errors, wide=wide, zero_difference=zero_difference),
        _interval_slot(ctx, tangent, errors, wide=wide, zero_difference=zero_difference),
        _interval_slot(ctx, second, errors, wide=wide, zero_difference=zero_difference),
        pressure_order=PRESSURE_ORDER,
        width_order=WIDTH_ORDER,
    )


def _check_expression_containment(
    interval_expressions: dict[str, IntervalAxialSecondJet],
    legacy_expressions: dict[str, AxialSecondJet],
    perturbed_expressions: list[dict[str, AxialSecondJet]],
    *,
    compare_difference: bool,
) -> tuple[int, mp.mpf]:
    checks = 0
    maximum_scaled_error = mp.mpf(0)
    slot_names = ("value", "tangent", "second")
    for name, interval_result in interval_expressions.items():
        legacy_result = legacy_expressions[name]
        for slot_name in slot_names:
            interval_slot = getattr(interval_result, slot_name)
            legacy_slot = getattr(legacy_result, slot_name)
            for pressure in range(PRESSURE_ORDER + 1):
                for width in range(WIDTH_ORDER + 1):
                    coefficient = interval_slot.component(pressure, width)
                    expected = legacy_slot.component(pressure, width)
                    if not _inside(expected, coefficient.nominal):
                        raise AssertionError((name, slot_name, pressure, width, "nominal"))
                    checks += 1
                    maximum_scaled_error = max(
                        maximum_scaled_error,
                        _scaled(
                            (sum(_bounds(coefficient.nominal)) / 2),
                            expected,
                        ),
                    )
                    if compare_difference:
                        for perturbed in perturbed_expressions:
                            target = (
                                getattr(perturbed[name], slot_name).component(pressure, width)
                                - expected
                            )
                            if not _inside(target, coefficient.difference):
                                raise AssertionError(
                                    (name, slot_name, pressure, width, "difference", target)
                                )
                            checks += 1
    return checks, maximum_scaled_error


def _maps(z0: mp.mpf) -> tuple[
    dict[tuple[int, int], mp.mpf],
    dict[tuple[int, int], mp.mpf],
    dict[tuple[int, int], mp.mpf],
]:
    value = _base_scalar(z0)
    tangent = mp.diff(_base_scalar, z0, 1)
    second = mp.diff(_base_scalar, z0, 2)
    value_map = {
        (0, 0): value,
        (1, 0): mp.mpf("0.13"),
        (2, 0): mp.mpf("0.009"),
        (0, 1): mp.mpf("0.071"),
        (1, 1): mp.mpf("0.006"),
        (2, 1): mp.mpf("0.0012"),
        (0, 2): mp.mpf("0.004"),
        (1, 2): mp.mpf("0.0007"),
        (2, 2): mp.mpf("0.0001"),
    }
    tangent_map = {
        (0, 0): tangent,
        (1, 0): mp.mpf("0.021"),
        (2, 0): mp.mpf("0.002"),
        (0, 1): mp.mpf("0.017"),
        (1, 1): mp.mpf("0.0011"),
        (2, 1): mp.mpf("0.0004"),
        (0, 2): mp.mpf("0.0017"),
        (1, 2): mp.mpf("0.0002"),
        (2, 2): mp.mpf("0.00005"),
    }
    second_map = {
        (0, 0): second,
        (1, 0): mp.mpf("0.011"),
        (2, 0): mp.mpf("0.001"),
        (0, 1): mp.mpf("0.008"),
        (1, 1): mp.mpf("0.0008"),
        (2, 1): mp.mpf("0.0002"),
        (0, 2): mp.mpf("0.0009"),
        (1, 2): mp.mpf("0.0001"),
        (2, 2): mp.mpf("0.00002"),
    }
    return value_map, tangent_map, second_map


def run_fixture(*, precision: int = PRECISION) -> dict[str, Any]:
    precision = max(100, int(precision))
    with mp.workdps(precision):
        ctx = MPIntervalContext()
        ctx.dps = precision
        z0 = mp.mpf("0.37")
        value_map, tangent_map, second_map = _maps(z0)
        keys = tuple(value_map)
        errors = {key: mp.mpf("2e-10") * (1 + sum(key)) for key in keys}

        legacy = _legacy_jet(value_map, tangent_map, second_map)
        interval = _interval_jet(
            ctx,
            value_map,
            tangent_map,
            second_map,
            errors,
            wide=False,
        )
        legacy_expressions = _jet_functions(legacy)
        interval_expressions = _jet_functions(interval)

        scales = (
            {key: mp.mpf(value) for key, value in zip(keys, ("-.7", ".2", ".8", ".4", "-.3", ".6", "-.5", ".1", ".9"))},
            {key: mp.mpf(value) for key, value in zip(keys, (".35", "-.6", ".15", "-.8", ".7", "-.2", ".25", "-.9", ".45"))},
            {key: mp.mpf(value) for key, value in zip(keys, (".9", ".4", "-.5", ".1", "-.75", ".3", ".6", "-.25", ".55"))},
        )
        perturbed_legacy = []
        for scale in scales:
            resolved_value = {key: value_map[key] + errors[key] * scale[key] for key in keys}
            resolved_tangent = {key: tangent_map[key] + errors[key] * scale[key] for key in keys}
            resolved_second = {key: second_map[key] + errors[key] * scale[key] for key in keys}
            perturbed_legacy.append(
                _jet_functions(_legacy_jet(resolved_value, resolved_tangent, resolved_second))
            )

        point_checks, maximum_point_error = _check_expression_containment(
            interval_expressions,
            legacy_expressions,
            perturbed_legacy,
            compare_difference=True,
        )

        # Widen nominal intervals independently of their perturbation ranges.
        # This verifies that chain-rule outputs retain a genuine enclosure when
        # the nominal coefficient itself is not a point.
        wide_interval = _interval_jet(
            ctx,
            value_map,
            tangent_map,
            second_map,
            errors,
            wide=True,
        )
        wide_expressions = _jet_functions(wide_interval)
        wide_checks, maximum_wide_error = _check_expression_containment(
            wide_expressions,
            legacy_expressions,
            perturbed_legacy,
            compare_difference=True,
        )

        # Exact zero perturbations must remain structurally exact in every
        # retained atom and all three slots.
        zero_interval = _interval_jet(
            ctx,
            value_map,
            tangent_map,
            second_map,
            errors,
            wide=False,
            zero_difference=True,
        )
        zero_expressions = _jet_functions(zero_interval)
        zero_difference_checks = 0
        for result in zero_expressions.values():
            for slot_name in ("value", "tangent", "second"):
                slot = getattr(result, slot_name)
                for pressure in range(PRESSURE_ORDER + 1):
                    for width in range(WIDTH_ORDER + 1):
                        assert _exact_zero(slot.component(pressure, width).difference)
                        zero_difference_checks += 1

        # Independent scalar derivative checks use the constant atom only and
        # ``mp.diff``; these do not share the jet chain-rule implementation.
        oracle_errors: dict[str, dict[str, mp.mpf]] = {}
        scalar_functions = _scalar_functions()
        for name, function in scalar_functions.items():
            expected = (
                function(z0),
                mp.diff(function, z0, 1),
                mp.diff(function, z0, 2),
            )
            actual_result = interval_expressions[name]
            actual = tuple(
                getattr(actual_result, slot_name).component(0, 0).nominal
                for slot_name in ("value", "tangent", "second")
            )
            oracle_errors[name] = {
                slot_name: _scaled((sum(_bounds(actual_value)) / 2), expected[index])
                for index, slot_name in enumerate(("value", "tangent", "second"))
                for actual_value in (actual[index],)
            }
            for index, slot_name in enumerate(("value", "tangent", "second")):
                assert _inside(expected[index], actual[index])

        # Strict domain failures are delegated from the interval ring.
        invalid = IntervalPressureWidthJet(
            ctx,
            {(0, 0): IntervalDifference(ctx, ctx.mpf(["-.1", ".1"]))},
            pressure_order=PRESSURE_ORDER,
            width_order=WIDTH_ORDER,
        )
        invalid_jet = IntervalAxialSecondJet(
            ctx,
            invalid,
            pressure_order=PRESSURE_ORDER,
            width_order=WIDTH_ORDER,
        )
        domain_failures = {}
        for name, method in (
            ("reciprocal", invalid_jet.reciprocal),
            ("log", invalid_jet.log),
            ("sqrt", invalid_jet.sqrt),
        ):
            try:
                method()
            except (ValueError, ZeroDivisionError):
                domain_failures[name] = True
            else:
                domain_failures[name] = False
        assert all(domain_failures.values())

        # Context/order checks and first-order dual rejection are intentional
        # API guards, not interval arithmetic tests.
        ctx_mismatch = MPIntervalContext()
        ctx_mismatch.dps = precision
        mismatch_ring = IntervalPressureWidthJet(
            ctx_mismatch,
            2,
            pressure_order=PRESSURE_ORDER,
            width_order=WIDTH_ORDER,
        )
        context_mismatch_rejected = False
        try:
            IntervalAxialSecondJet(
                ctx,
                mismatch_ring,
                pressure_order=PRESSURE_ORDER,
                width_order=WIDTH_ORDER,
            )
        except ValueError:
            context_mismatch_rejected = True
        order_mismatch_rejected = False
        mismatch_order_ring = IntervalPressureWidthJet(
            ctx,
            2,
            pressure_order=PRESSURE_ORDER + 1,
            width_order=WIDTH_ORDER,
        )
        try:
            interval + mismatch_order_ring
        except ValueError:
            order_mismatch_rejected = True
        dual = AxialDual(2, mp.mpf(".3"))
        dual_promotion_rejected = False
        try:
            IntervalAxialSecondJet(
                ctx,
                dual,
                second=IntervalPressureWidthJet(
                    ctx,
                    0,
                    pressure_order=PRESSURE_ORDER,
                    width_order=WIDTH_ORDER,
                ),
                pressure_order=PRESSURE_ORDER,
                width_order=WIDTH_ORDER,
            )
        except TypeError:
            dual_promotion_rejected = True
        dual_operand_rejected = False
        try:
            interval + dual
        except TypeError:
            dual_operand_rejected = True
        scalar_orders_required = False
        try:
            IntervalAxialSecondJet(ctx, 2, mp.mpf(".3"), mp.mpf(".1"))
        except TypeError:
            scalar_orders_required = True

        checks = {
            "point_nominal_and_perturbation_containment": point_checks > 0,
            "wide_nominal_and_perturbation_containment": wide_checks > 0,
            "zero_difference_exact": zero_difference_checks > 0,
            "independent_mp_diff_oracle": all(
                error <= ERROR_THRESHOLD
                for errors_for_function in oracle_errors.values()
                for error in errors_for_function.values()
            ),
            "strict_domain_failures": all(domain_failures.values()),
            "context_mismatch_rejected": context_mismatch_rejected,
            "order_mismatch_rejected": order_mismatch_rejected,
            "first_order_dual_promotion_rejected": dual_promotion_rejected,
            "first_order_dual_operand_rejected": dual_operand_rejected,
            "scalar_orders_required": scalar_orders_required,
        }
        if not all(checks.values()):
            failed = [name for name, passed in checks.items() if not passed]
            raise AssertionError(f"IntervalAxialSecondJet fixture checks failed: {failed}")

        report = {
            "passed": True,
            "class": "IntervalAxialSecondJet",
            "precision": precision,
            "pressure_order": PRESSURE_ORDER,
            "width_order": WIDTH_ORDER,
            "actual_derivative_slots": True,
            "automatic_first_Z_derivative": True,
            "automatic_second_Z_derivative": True,
            "directed_interval_arithmetic": True,
            "checks": checks,
            "operations": list(legacy_expressions),
            "point_containment_checks": point_checks,
            "wide_containment_checks": wide_checks,
            "zero_difference_checks": zero_difference_checks,
            "maximum_scaled_point_nominal_rounding": mp.nstr(maximum_point_error, 30),
            "maximum_scaled_wide_nominal_rounding": mp.nstr(maximum_wide_error, 30),
            "oracle": {
                name: {
                    key: mp.nstr(error, 18) for key, error in errors_for_function.items()
                }
                for name, errors_for_function in oracle_errors.items()
            },
            "domain_failures": domain_failures,
            "limitations": [
                "Pressure/width products are truncated to the declared rectangle.",
                "The three axial slots are local actual derivatives; no global ODE remainder is claimed.",
                "No full core, collar, field, or paper-parameter rerun is performed.",
            ],
        }
        output = Path(__file__).with_suffix(".json")
        output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        return report


def fixture(**kwargs: Any) -> dict[str, Any]:
    return run_fixture(**kwargs)


def run_fixtures(**kwargs: Any) -> dict[str, Any]:
    return run_fixture(**kwargs)


def main() -> None:
    print(json.dumps(run_fixture(), indent=2, sort_keys=True))


__all__ = ["fixture", "run_fixture", "run_fixtures"]


if __name__ == "__main__":
    main()
