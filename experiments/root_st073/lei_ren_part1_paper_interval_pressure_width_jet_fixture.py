"""Independent checks for the interval pressure/width atom ring."""

from __future__ import annotations

import json
from pathlib import Path

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext

from lei_ren_part1_paper_interval_difference import IntervalDifference
from lei_ren_part1_paper_interval_pressure_width_jet import IntervalPressureWidthJet
from lei_ren_part1_paper_pressure_width_jet import PressureWidthJet


def _bounds(value):
    return tuple(mp.make_mpf(item) for item in value._mpi_)


def _inside(value, interval):
    lower, upper = _bounds(interval)
    return lower <= value <= upper


def _exact_zero(value):
    lower, upper = _bounds(value)
    return lower == 0 and upper == 0


def _resolved(mapping, errors, scale):
    return {
        key: value + errors[key] * scale[key]
        for key, value in mapping.items()
    }


def main():
    # The decimal strings are retained as the source of truth for both the
    # old resolved ring and the directed interval ring.
    with mp.workdps(120):
        ctx = MPIntervalContext()
        ctx.dps = 120
        pressure_order = 2
        width_order = 2
        mapping = {
            (0, 0): "2.0",
            (1, 0): "0.1",
            (2, 0): "1e-82",
            (0, 1): "1e-80",
            (1, 1): "0.003",
            (2, 1): "7e-5",
            (0, 2): "2e-90",
            (1, 2): "-2e-6",
        }
        errors = {
            (0, 0): "1e-7",
            (1, 0): "2e-7",
            (2, 0): "3e-87",
            (0, 1): "1e-85",
            (1, 1): "2e-8",
            (2, 1): "4e-9",
            (0, 2): "1e-95",
            (1, 2): "5e-10",
        }
        interval_atoms = {
            key: IntervalDifference(
                ctx,
                ctx.mpf(value),
                ctx.mpf([-ctx.mpf(errors[key]), ctx.mpf(errors[key])]),
            )
            for key, value in mapping.items()
        }
        interval_jet = IntervalPressureWidthJet(
            ctx,
            interval_atoms,
            pressure_order=pressure_order,
            width_order=width_order,
        )
        assert (2, 0) in interval_jet.atoms
        assert (0, 1) in interval_jet.atoms and (0, 2) in interval_jet.atoms
        assert (2, 1) in interval_jet.atoms and (1, 2) in interval_jet.atoms

        def operation(jet, name):
            if name == "add_sub":
                return (jet + jet * ".25") - jet * ".25"
            if name == "multiply":
                return jet * jet
            if name == "divide":
                return jet / (jet + 2)
            if name == "reciprocal":
                return jet.reciprocal()
            if name == "power":
                return jet**3
            if name == "exp":
                return jet.exp()
            if name == "log":
                return jet.log()
            if name == "sqrt":
                return jet.sqrt()
            if name == "mixed":
                return ((jet * jet + 1) / (jet + 2)) ** 3
            raise AssertionError(name)

        names = (
            "add_sub",
            "multiply",
            "divide",
            "reciprocal",
            "power",
            "exp",
            "log",
            "sqrt",
            "mixed",
        )
        # Distinct deterministic perturbation vectors exercise independent
        # coefficient errors, including the tiny width atoms.
        scales = (
            {
                key: mp.mpf(value)
                for key, value in zip(mapping, ("-.7", ".2", ".8", ".4", "-.3", ".6", "-.5", ".1"))
            },
            {
                key: mp.mpf(value)
                for key, value in zip(mapping, (".35", "-.6", ".15", "-.8", ".7", "-.2", ".25", "-.9"))
            },
            {
                key: mp.mpf(value)
                for key, value in zip(mapping, (".9", ".4", "-.5", ".1", "-.75", ".3", ".6", "-.25"))
            },
        )

        coefficient_checks = 0
        perturbation_checks = 0
        for name in names:
            interval_result = operation(interval_jet, name)
            nominal_result = operation(
                PressureWidthJet(
                    mapping,
                    pressure_order=pressure_order,
                    width_order=width_order,
                ),
                name,
            )
            for pressure in range(pressure_order + 1):
                for width in range(width_order + 1):
                    interval_coefficient = interval_result.component(pressure, width)
                    assert _inside(
                        nominal_result.component(pressure, width),
                        interval_coefficient.nominal,
                    ), (name, pressure, width)
                    coefficient_checks += 1
            for scale in scales:
                resolved = PressureWidthJet(
                    _resolved(
                        {key: mp.mpf(value) for key, value in mapping.items()},
                        {key: mp.mpf(value) for key, value in errors.items()},
                        scale,
                    ),
                    pressure_order=pressure_order,
                    width_order=width_order,
                )
                resolved_result = operation(resolved, name)
                for pressure in range(pressure_order + 1):
                    for width in range(width_order + 1):
                        target_difference = (
                            resolved_result.component(pressure, width)
                            - nominal_result.component(pressure, width)
                        )
                        assert _inside(
                            target_difference,
                            interval_result.component(pressure, width).difference,
                        ), (name, pressure, width, target_difference)
                        perturbation_checks += 1

        # A wide nominal constant exposes accidental interval cancellation:
        # subtracting an interval from itself can leave a nonzero enclosure.
        # The nonlinear ring must remove the constant atom structurally.
        shared_difference = ctx.mpf(["-.08", ".08"])
        wide_atoms = {
            (0, 0): IntervalDifference(
                ctx,
                ctx.mpf(["1.8", "2.2"]),
                shared_difference,
            ),
            (1, 0): IntervalDifference(
                ctx,
                ctx.mpf(".1"),
                shared_difference / 10,
            ),
            (2, 0): IntervalDifference(ctx, ctx.mpf("1e-82"), ctx.mpf(["-1e-87", "1e-87"])),
            (0, 1): IntervalDifference(ctx, ctx.mpf("1e-80"), ctx.mpf(["-1e-85", "1e-85"])),
            (1, 1): IntervalDifference(ctx, ctx.mpf(".003"), shared_difference / 100),
        }
        wide_jet = IntervalPressureWidthJet(
            ctx,
            wide_atoms,
            pressure_order=pressure_order,
            width_order=width_order,
        )
        assert (0, 0) not in wide_jet._nonconstant().atoms
        assert (0, 0) not in wide_jet._nonconstant(wide_jet.constant.reciprocal()).atoms

        wide_points = (
            {
                (0, 0): mp.mpf("1.74"),
                (1, 0): mp.mpf(".094"),
                (2, 0): mp.mpf("1e-82"),
                (0, 1): mp.mpf("1e-80"),
                (1, 1): mp.mpf(".0024"),
            },
            {
                (0, 0): mp.mpf("2.0"),
                (1, 0): mp.mpf(".1"),
                (2, 0): mp.mpf("1.0005e-82"),
                (0, 1): mp.mpf("1.0005e-80"),
                (1, 1): mp.mpf(".003"),
            },
            {
                (0, 0): mp.mpf("2.25"),
                (1, 0): mp.mpf(".106"),
                (2, 0): mp.mpf(".9995e-82"),
                (0, 1): mp.mpf(".9995e-80"),
                (1, 1): mp.mpf(".0036"),
            },
        )
        wide_operation_names = ("reciprocal", "exp", "log", "sqrt", "mixed")
        wide_point_checks = 0
        for name in wide_operation_names:
            interval_result = operation(wide_jet, name)
            for point in wide_points:
                resolved_result = operation(
                    PressureWidthJet(
                        point,
                        pressure_order=pressure_order,
                        width_order=width_order,
                    ),
                    name,
                )
                for pressure in range(pressure_order + 1):
                    for width in range(width_order + 1):
                        assert _inside(
                            resolved_result.component(pressure, width),
                            interval_result.component(pressure, width).value,
                        ), (name, pressure, width)
                        wide_point_checks += 1

        zero_atoms = {
            key: IntervalDifference(ctx, ctx.mpf(value), ctx.mpf(0))
            for key, value in mapping.items()
        }
        zero_jet = IntervalPressureWidthJet(
            ctx,
            zero_atoms,
            pressure_order=pressure_order,
            width_order=width_order,
        )
        zero_difference_exact = True
        for name in names:
            result = operation(zero_jet, name)
            for pressure in range(pressure_order + 1):
                for width in range(width_order + 1):
                    zero_difference_exact &= _exact_zero(
                        result.component(pressure, width).difference
                    )
        assert zero_difference_exact

        # Strict domains must reject a nominal or total constant interval that
        # reaches zero.  Reciprocal retains its separate nonzero check.
        invalid = IntervalPressureWidthJet(
            ctx,
            {(0, 0): IntervalDifference(ctx, ctx.mpf(["-.1", ".1"]))},
            pressure_order=pressure_order,
            width_order=width_order,
        )
        strict_domain_rejections = 0
        for method in (invalid.log, invalid.sqrt):
            try:
                method()
            except ValueError:
                strict_domain_rejections += 1
        assert strict_domain_rejections == 2

        report = {
            "all_checks_passed": True,
            "directed_interval_arithmetic": True,
            "coefficient_checks": coefficient_checks,
            "perturbation_checks": perturbation_checks,
            "wide_nominal_point_checks": wide_point_checks,
            "operation_names": list(names),
            "mixed_atoms_preserved": True,
            "tiny_pressure_width_atoms_preserved": True,
            "nonconstant_constant_atom_removed_structurally": True,
            "zero_difference_exact": True,
            "strict_positive_domain_rejections": strict_domain_rejections,
            "rectangular_truncation": True,
            "source_certificate": False,
        }
    output = Path(__file__).with_suffix(".json")
    output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"all_checks_passed": True, "output": str(output)}))


if __name__ == "__main__":
    main()
