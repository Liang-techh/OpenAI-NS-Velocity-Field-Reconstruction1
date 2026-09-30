"""Independent finite-difference fixture for the axial pressure-width bridge.

The production bridge carries first ``Z`` tangents through the finite
comparison and exit RK construction.  This module builds a resolved local
polynomial core with at least three ``Z`` coefficients in every radial row,
then compares those tangents with fourth-order centered finite differences of
fresh scalar ``PressureWidthComparison`` plus ``PressureWidthExitBridge``
replays at shifted ``Z`` values.

The radial velocity ``Ur`` is reconstructed independently at the resolved
point from the prescribed bridge value ``Uz``, moment ``M_z``, and a
fourth-order finite difference of the prescribed ``M_{z,Z}`` from the same
scalar replays.  The production bridge exposes ``Ur`` as a value; its
second-prescribed-Z input is kept out of this derivative fixture.

This is a resolved derivative check of one finite algorithm.  It is not an
ODE remainder bound, a formal radial derivative certificate, or a global
field or Navier--Stokes closure claim.
"""

from __future__ import annotations

from collections.abc import Mapping
import json
from math import comb
from types import SimpleNamespace
from pathlib import Path
from typing import Any

import mpmath as mp

try:
    from .lei_ren_part1_paper_component_pressure_core import PressurePolynomial
    from .lei_ren_part1_paper_core_recursion import core_coefficients
    from .lei_ren_part1_paper_pressure_width_axial_bridge import (
        AxialPressureWidthExitBridge,
    )
    from .lei_ren_part1_paper_pressure_width_axial_comparison import (
        AxialPressureWidthComparison,
    )
    from .lei_ren_part1_paper_pressure_width_bridge import PressureWidthExitBridge
    from .lei_ren_part1_paper_pressure_width_comparison import PressureWidthComparison
except (ImportError, ValueError):
    from lei_ren_part1_paper_component_pressure_core import PressurePolynomial
    from lei_ren_part1_paper_core_recursion import core_coefficients
    from lei_ren_part1_paper_pressure_width_axial_bridge import (
        AxialPressureWidthExitBridge,
    )
    from lei_ren_part1_paper_pressure_width_axial_comparison import (
        AxialPressureWidthComparison,
    )
    from lei_ren_part1_paper_pressure_width_bridge import PressureWidthExitBridge
    from lei_ren_part1_paper_pressure_width_comparison import PressureWidthComparison


PRESSURE_ORDER = 3
WIDTH_ORDER = 2
REQUIRED_Z_DEPTH = 2
REQUIRED_ROW_LENGTH = REQUIRED_Z_DEPTH + 1
FINITE_DIFFERENCE_STEP = mp.mpf("1e-4")
FINITE_DIFFERENCE_THRESHOLD = mp.mpf("1e-8")


def _shifted_rows(formal: Mapping[str, Any], offset: mp.mpf) -> dict[str, Any]:
    """Re-center the formal coefficient rows at ``Z + offset``."""

    result = dict(formal)
    for name in ("F", "Uz", "P"):
        rows: list[list[PressurePolynomial]] = []
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


def _make_formal_core(
    *,
    z0: mp.mpf,
    precision: int,
) -> tuple[SimpleNamespace, dict[str, Any]]:
    """Make a resolved pressure-tail family with depth-two Z rows."""

    f0 = [mp.mpf(value) for value in ("2", ".3", ".02", ".001", ".0004", ".0001")]
    u0 = [mp.mpf(value) for value in (".4", ".2", ".01", ".001", ".0004", ".0001")]
    p0 = [mp.mpf(value) for value in ("-1", ".2", ".03", ".002", ".0003", ".0001")]
    pressure_tail = [
        mp.mpf(value)
        for value in (".01", ".03", "-.02", ".001", ".0003", ".0001")
    ]
    formal = core_coefficients(
        z0,
        ".01",
        F0_Z_taylor=f0,
        U0_Z_taylor=u0,
        P0_Z_taylor=[
            PressurePolynomial({0: value, 1: tail})
            for value, tail in zip(p0, pressure_tail)
        ],
        radial_degree=3,
        precision=precision,
        scalar_converter=PressurePolynomial,
    )
    row_lengths = {
        name: [len(row) for row in formal[name]] for name in ("F", "Uz", "P")
    }
    if min(length for lengths in row_lengths.values() for length in lengths) < REQUIRED_ROW_LENGTH:
        raise AssertionError(f"Insufficient Z row depth: {row_lengths}")
    axis = SimpleNamespace(
        precision=precision,
        Lambda=mp.mpf(10),
        delta=mp.mpf(".01"),
    )
    component = SimpleNamespace(
        axis=axis,
        precision=precision,
        Z_jet_depth=REQUIRED_Z_DEPTH,
        radial_degree=3,
        coefficients=lambda _Z, data=formal: data,
    )
    return component, formal


def _bundle(axis: Any, component: Any) -> dict[str, Any]:
    return {
        "axis": axis,
        "Lambda": axis.Lambda,
        "precision": component.precision,
        "radial_degree": 3,
        "component_pressure_core": component,
    }


def _plain_bridge_for_offset(
    axis: Any,
    formal: Mapping[str, Any],
    offset: mp.mpf,
    *,
    precision: int,
    h_b: mp.mpf,
) -> PressureWidthExitBridge:
    component = SimpleNamespace(
        axis=axis,
        precision=precision,
        radial_degree=3,
        Z_jet_depth=REQUIRED_Z_DEPTH,
        coefficients=lambda _Z, data=_shifted_rows(formal, offset): data,
    )
    comparison = PressureWidthComparison(
        _bundle(axis, component),
        h_b=h_b,
        pressure_order=PRESSURE_ORDER,
        width_order=WIDTH_ORDER,
        transition_steps=4,
    )
    return PressureWidthExitBridge(comparison, steps=8)


def _fourth_difference(values: Mapping[int, mp.mpf], step: mp.mpf) -> mp.mpf:
    return (
        -values[2]
        + 8 * values[1]
        - 8 * values[-1]
        + values[-2]
    ) / (12 * step)


def run_fixture(
    *,
    precision: int = 100,
    finite_difference_step: Any = FINITE_DIFFERENCE_STEP,
) -> dict[str, Any]:
    """Compare axial bridge tangents against independent fourth-order replays."""

    precision = max(70, int(precision))
    step = mp.mpf(str(finite_difference_step))
    if step <= 0:
        raise ValueError("finite_difference_step must be positive")
    with mp.workdps(precision):
        z0 = mp.mpf(".3")
        h_b = mp.mpf("1e-4")
        component, formal = _make_formal_core(z0=z0, precision=precision)
        axis = component.axis
        base = _bundle(axis, component)
        axial_comparison = AxialPressureWidthComparison(
            base,
            h_b=h_b,
            pressure_order=PRESSURE_ORDER,
            width_order=WIDTH_ORDER,
            transition_steps=4,
        )
        axial_bridge = AxialPressureWidthExitBridge(axial_comparison, steps=8)

        offsets = {-2: -2 * step, -1: -step, 1: step, 2: 2 * step}
        plain_bridges = {
            index: _plain_bridge_for_offset(
                axis,
                formal,
                offset,
                precision=precision,
                h_b=h_b,
            )
            for index, offset in offsets.items()
        }

        fields = ("F", "Uz", "P")
        moment_names = ("theta", "z", "theta_z", "z_theta", "p")
        sample_s = (mp.mpf(".5"), mp.mpf("1"), mp.mpf("1.5"), mp.mpf("2"))
        errors: dict[str, mp.mpf] = {}
        ur_value_errors: dict[str, mp.mpf] = {}
        source_pressure_atoms: list[mp.mpf] = []
        source_width_atoms: list[mp.mpf] = []

        for s in sample_s:
            axial_output = axial_bridge.evaluate(s, z0)
            field_tangents = {
                name: axial_output[name + "_Z"].evaluate(
                    pressure=1,
                    width=1,
                )
                for name in fields
            }
            source_pressure_atoms.append(
                axial_output["P"].component(1, 0)
            )
            source_width_atoms.append(axial_output["R"].component(0, 1))

            for name in fields:
                values = {
                    index: plain_bridges[index]
                    .evaluate(s, z0 + offsets[index])[name]
                    .evaluate(pressure=1, width=1)
                    for index in offsets
                }
                finite_difference = _fourth_difference(values, step)
                error = abs(field_tangents[name] - finite_difference) / max(
                    mp.mpf(1), abs(finite_difference), abs(field_tangents[name])
                )
                errors[f"{name}@s={mp.nstr(s, 3)}"] = error

            for name in moment_names:
                tangent = axial_output["moments_Z"][name].evaluate(
                    pressure=1,
                    width=1,
                )
                values = {
                    index: plain_bridges[index]
                    .evaluate(s, z0 + offsets[index])["moments"][name]
                    .evaluate(pressure=1, width=1)
                    for index in offsets
                }
                finite_difference = _fourth_difference(values, step)
                error = abs(tangent - finite_difference) / max(
                    mp.mpf(1), abs(finite_difference), abs(tangent)
                )
                errors[f"moment:{name}@s={mp.nstr(s, 3)}"] = error

            # The axial bridge exposes Ur as a value.  Reconstruct it from
            # the prescribed endpoint fields and M_z,Z obtained by the same
            # scalar replays.  This avoids substituting auxiliary comparison
            # moments_Z for the prescribed bridge moment derivative.
            mz_values = {
                index: plain_bridges[index]
                .evaluate(s, z0 + offsets[index])["moments"]["z"]
                .evaluate(pressure=1, width=1)
                for index in offsets
            }
            prescribed_mz_z = _fourth_difference(mz_values, step)
            radius = axial_output["R"].evaluate(pressure=1, width=1)
            uz = axial_output["Uz"].evaluate(pressure=1, width=1)
            mz = axial_output["moments"]["z"].evaluate(
                pressure=1,
                width=1,
            )
            z_value = z0
            delta = axis.delta
            radial_root = mp.sqrt(2 * radius)
            L = 1 - delta * z_value * z_value
            d = 1 - z_value * z_value
            reconstructed_ur = (
                2 * z_value * radius * uz
                - (1 - delta) * z_value * mz
                - d * prescribed_mz_z
            ) / (L * radial_root)
            ur_value = axial_output["Ur"].evaluate(pressure=1, width=1)
            ur_value_errors[f"s={mp.nstr(s, 3)}"] = abs(
                ur_value - reconstructed_ur
            ) / max(mp.mpf(1), abs(ur_value), abs(reconstructed_ur))

        chart = axial_bridge.physical_chart(
            mp.mpf("1"),
            z0,
            mp.mpf("-4"),
            nu=mp.mpf(".01"),
            phi=mp.mpf(".2"),
        )
        chart_values = tuple(
            value.evaluate(pressure=1, width=1)
            if hasattr(value, "evaluate")
            else mp.mpf(value)
            for value in chart["uvw"]
        )
        chart_smoke = len(chart["xyz"]) == 3 and len(chart_values) == 3 and all(
            mp.isfinite(value) for value in chart_values
        )
        maximum = max(errors.values(), default=mp.mpf(0))
        maximum_ur_value = max(ur_value_errors.values(), default=mp.mpf(0))
        meaningful = maximum < FINITE_DIFFERENCE_THRESHOLD
        if not meaningful:
            raise AssertionError(
                f"axial bridge tangent finite-difference error {mp.nstr(maximum, 20)}"
            )
        if maximum_ur_value >= mp.mpf("1e-12"):
            raise AssertionError(
                "reconstructed Ur value disagrees with axial bridge: "
                f"{mp.nstr(maximum_ur_value, 20)}"
            )
        if not all(value != 0 for value in source_pressure_atoms):
            raise AssertionError("pressure-tail P atom vanished in resolved fixture")
        if not all(value != 0 for value in source_width_atoms):
            raise AssertionError("width R atom vanished in resolved fixture")
        if not chart_smoke:
            raise AssertionError("physical chart smoke failed")

        report = {
            "passed": True,
            "maximum_scaled_tangent_error": mp.nstr(maximum, 30),
            "maximum_reconstructed_Ur_value_error": mp.nstr(maximum_ur_value, 30),
            "errors": {
                key: mp.nstr(value, 18) for key, value in sorted(errors.items())
            },
            "ur_value_errors": {
                key: mp.nstr(value, 18)
                for key, value in sorted(ur_value_errors.items())
            },
            "pressure_order": PRESSURE_ORDER,
            "width_order": WIDTH_ORDER,
            "pressure_parameter": 1,
            "width_parameter": 1,
            "finite_difference_step": mp.nstr(step, 18),
            "finite_difference_formula": "(-f(2h)+8f(h)-8f(-h)+f(-2h))/(12h)",
            "required_component_Z_depth": REQUIRED_Z_DEPTH,
            "required_component_Z_row_length": REQUIRED_ROW_LENGTH,
            "all_radial_rows_have_depth_two": True,
            "pressure_tail_atoms_nonzero": True,
            "width_atoms_nonzero": True,
            "physical_chart_smoke": chart_smoke,
            "physical_chart_uvw": [mp.nstr(value, 18) for value in chart_values],
            "derivative_scope": "finite comparison plus finite exit-RK algorithm replay",
            "axial_tangent_remainder_enclosed": False,
            "ode_remainder_enclosed": False,
            "global_field_installed": False,
            "limitations": [
                "The finite differences validate the local finite algorithm at resolved scales.",
                "No ODE remainder, radial truncation remainder, global matching, or stress-cone bound is enclosed.",
            ],
        }
        Path(__file__).with_suffix(".json").write_text(
            json.dumps(report, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        return report


def fixture(**kwargs: Any) -> dict[str, Any]:
    return run_fixture(**kwargs)


def main() -> None:
    print(json.dumps(run_fixture(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
