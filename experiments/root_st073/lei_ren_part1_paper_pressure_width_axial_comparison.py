"""Axial automatic differentiation for the pressure/width comparison.

The pressure-tail and exit-width atoms stay in ``PressureWidthJet`` while an
``AxialDual`` carries the first derivative with respect to the physical
coordinate ``Z``.  The finite Section 9.23 state equations are inherited from
``PressureWidthComparison``; this adapter only supplies the differentiated
coefficient rows, core evaluation, and the stress evaluator's axial
coordinate.

The component core must retain at least three Z Taylor coefficients per
radial row.  The third coefficient is needed to make the tangent of the
first-Z row available while the radial core equations are evaluated.
"""

from __future__ import annotations

from collections.abc import Mapping
from math import comb
from types import SimpleNamespace
from typing import Any

import mpmath as mp

from lei_ren_part1_paper_axial_dual import AxialDual
from lei_ren_part1_paper_component_pressure_core import (
    PressurePolynomial,
    evaluate_component_core_coefficients,
)
from lei_ren_part1_paper_core_recursion import core_coefficients
from lei_ren_part1_paper_pressure_width_comparison import PressureWidthComparison
from lei_ren_part1_paper_pressure_width_jet import PressureWidthJet


class AxialPressureWidthComparison(PressureWidthComparison):
    """Pressure/width comparison with an automatic first Z derivative."""

    automatic_Z_tangent = True
    required_component_Z_depth = 2
    required_component_Z_row_length = 3

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        super().__init__(*args, **kwargs)
        self._axial_coeff_cache: dict[str, dict[str, Any]] = {}

    def _base_jet(self, value: Any = 0) -> PressureWidthJet:
        """Coerce a value through the parent bivariate ring only."""

        return PressureWidthComparison.jet(self, value)

    def _dual(self, value: Any = 0, tangent: Any = 0) -> AxialDual:
        """Construct a dual with this comparison's fixed bivariate orders."""

        if isinstance(value, AxialDual):
            value = value.value
        if isinstance(tangent, AxialDual):
            tangent = tangent.tangent
        return AxialDual(self._base_jet(value), self._base_jet(tangent))

    def jet(self, value: Any = 0) -> AxialDual:
        """Construct an axial dual over this comparison's pressure/width ring."""

        if isinstance(value, AxialDual):
            if value.orders == (self.pressure_order, self.width_order):
                return value
            return self._dual(value.value, value.tangent)
        return self._dual(value, 0)

    def _stress_axial_coordinate(self, Z: Any) -> AxialDual:
        """Return the physical axial coordinate with derivative one."""

        return self._dual(Z, 1)

    def _differentiate_rows(
        self, rows: list[list[Any] | tuple[Any, ...]]
    ) -> list[list[AxialDual]]:
        """Turn each Z Taylor row into value/first-derivative dual entries."""

        differentiated: list[list[AxialDual]] = []
        for row_number, row in enumerate(rows):
            if len(row) < self.required_component_Z_row_length:
                raise ValueError(
                    "component pressure core requires at least three Z Taylor "
                    f"coefficients in radial row {row_number}; got {len(row)}"
                )
            output: list[AxialDual] = []
            for k, coefficient in enumerate(row):
                coefficient_value = (
                    coefficient.value
                    if isinstance(coefficient, AxialDual)
                    else coefficient
                )
                next_coefficient = (
                    row[k + 1]
                    if k + 1 < len(row)
                    else self._base_jet(0)
                )
                next_value = (
                    next_coefficient.value
                    if isinstance(next_coefficient, AxialDual)
                    else next_coefficient
                )
                output.append(
                    self._dual(coefficient_value, (k + 1) * next_value)
                )
            differentiated.append(output)
        return differentiated

    def core_coefficients(self, Z: Any) -> dict[str, Any]:
        """Return component rows whose entries carry automatic Z tangents."""

        with mp.workdps(self.precision):
            key = mp.nstr(mp.mpf(str(Z)), self.precision)
            cached = self._axial_coeff_cache.get(key)
            if cached is not None:
                return cached

            # The parent converts only the pressure polynomial into the
            # PressureWidthJet ring and leaves scalar metadata untouched.
            base = super().core_coefficients(Z)
            result = dict(base)
            for name in ("F", "Uz", "P"):
                rows = base.get(name)
                if not isinstance(rows, list):
                    raise TypeError(f"component core {name} rows must be a list")
                result[name] = self._differentiate_rows(rows)
            result["automatic_Z_tangent"] = True
            result["required_component_Z_depth"] = self.required_component_Z_depth
            result["required_component_Z_row_length"] = self.required_component_Z_row_length
            self._axial_coeff_cache[key] = result
            return result

    def _mixed_rows(
        self,
        coefficients: Mapping[str, Any],
        radius: AxialDual,
    ) -> tuple[AxialDual, AxialDual]:
        """Compute F_RZ and Uz_RZ as duals, including their Z tangents."""

        f_rz = self.jet(0)
        uz_rz = self.jet(0)
        for n, row in enumerate(coefficients["F"]):
            if n:
                f_rz += n * row[1] * radius ** (n - 1)
        for n, row in enumerate(coefficients["Uz"]):
            if n:
                uz_rz += n * row[1] * radius ** (n - 1)
        return f_rz, uz_rz

    def _evaluate_core(
        self,
        coefficients: dict[str, Any],
        radius: AxialDual,
        Z: mp.mpf,
    ) -> dict[str, Any]:
        """Evaluate the shared core equations in the axial dual ring."""

        axial_coordinate = self._stress_axial_coordinate(Z)
        value = evaluate_component_core_coefficients(
            coefficients,
            radius,
            axial_coordinate,
            self.delta,
        )
        f_rz, uz_rz = self._mixed_rows(coefficients, radius)
        value["F_RZ"] = f_rz
        value["Uz_RZ"] = uz_rz
        value["R"] = radius
        value["Z"] = axial_coordinate
        return value

    def metadata(self) -> dict[str, Any]:
        result = super().metadata()
        result.update(
            automatic_Z_tangent=True,
            required_component_Z_depth=self.required_component_Z_depth,
            required_component_Z_row_length=self.required_component_Z_row_length,
            axial_dual_ring=True,
            axial_tangent_remainder_enclosed=False,
        )
        return result


# A descriptive alias for callers that name the adapter by its two jets.
PressureWidthAxialComparison = AxialPressureWidthComparison


def fixture() -> dict[str, Any]:
    """Compare dual driver tangents with an independent resolved replay."""

    def shifted(formal: Mapping[str, Any], offset: mp.mpf) -> dict[str, Any]:
        result = dict(formal)
        for name in ("F", "Uz", "P"):
            rows = []
            for row in formal[name]:
                rows.append(
                    [
                        sum(
                            (
                                mp.mpf(comb(k, j))
                                * offset ** (k - j)
                                * row[k]
                                for k in range(j, len(row))
                            ),
                            PressurePolynomial(0),
                        )
                        for j in range(len(row))
                    ]
                )
            result[name] = rows
        return result

    with mp.workdps(100):
        z0 = mp.mpf(".3")
        finite_difference_step = mp.mpf("1e-7")
        h_b = mp.mpf("1e-4")
        f0 = [mp.mpf(v) for v in ("2", ".3", ".02", ".001", ".0004", ".0001")]
        u0 = [mp.mpf(v) for v in (".4", ".2", ".01", ".001", ".0004", ".0001")]
        p0 = [mp.mpf(v) for v in ("-1", ".2", ".03", ".002", ".0003", ".0001")]
        formal = core_coefficients(
            z0,
            ".01",
            F0_Z_taylor=f0,
            U0_Z_taylor=u0,
            P0_Z_taylor=[PressurePolynomial({0: value}) for value in p0],
            radial_degree=3,
            precision=100,
            scalar_converter=PressurePolynomial,
        )
        axis = SimpleNamespace(
            precision=100,
            Lambda=mp.mpf(10),
            delta=mp.mpf(".01"),
        )
        component = SimpleNamespace(
            axis=axis,
            precision=100,
            coefficients=lambda _Z: formal,
        )
        comparison = AxialPressureWidthComparison(
            {
                "axis": axis,
                "Lambda": axis.Lambda,
                "precision": 100,
                "component_pressure_core": component,
            },
            h_b=h_b,
            pressure_order=2,
            width_order=1,
            transition_steps=4,
        )
        s = mp.mpf("1.5")
        dual_output = comparison.evaluate(s, z0)

        def scalar_driver(offset: mp.mpf, name: str) -> mp.mpf:
            scalar_component = SimpleNamespace(
                axis=axis,
                precision=100,
                coefficients=lambda _Z, data=shifted(formal, offset): data,
            )
            scalar_comparison = PressureWidthComparison(
                {
                    "axis": axis,
                    "Lambda": axis.Lambda,
                    "precision": 100,
                    "component_pressure_core": scalar_component,
                },
                h_b=h_b,
                pressure_order=2,
                width_order=1,
                transition_steps=4,
            )
            value = scalar_comparison.evaluate(s, z0 + offset)[name]
            return value.evaluate(pressure=0, width=1)

        errors: dict[str, mp.mpf] = {}
        for name in ("F", "Uz", "P", "I_theta", "I_z", "D", "E"):
            plus = scalar_driver(finite_difference_step, name)
            minus = scalar_driver(-finite_difference_step, name)
            finite_difference = (plus - minus) / (2 * finite_difference_step)
            tangent = dual_output[name].tangent.evaluate(pressure=0, width=1)
            errors[name] = abs(tangent - finite_difference) / max(
                mp.mpf(1), abs(finite_difference)
            )
        maximum = max(errors.values())
        if maximum >= mp.mpf("1e-10"):
            raise AssertionError(maximum)
        return {
            "available": True,
            "maximum_scaled_difference": mp.nstr(maximum, 30),
            "drivers_checked": list(errors),
            "automatic_Z_tangent": True,
            "required_component_Z_depth": 2,
            "required_component_Z_row_length": 3,
            "pressure_width_orders_remain_separate": True,
            "axial_tangent_remainder_enclosed": False,
            "global_field_installed": False,
        }


__all__ = [
    "AxialPressureWidthComparison",
    "PressureWidthAxialComparison",
    "fixture",
]


if __name__ == "__main__":
    print(fixture())
