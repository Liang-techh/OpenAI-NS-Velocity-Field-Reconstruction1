"""Independent resolved check for the finite pressure/width continuation.

The production continuation uses coefficient-valued exponential polynomials.
This fixture materializes the same retained pressure/width atoms at
``(pressure, width) = (1, 1)`` and integrates the prescribed equations with
an independent scalar RK4 replay in ``y = log(R / R_b)``.  It checks the two
fields, pressure, and all five cumulative moments without comparing against
the legacy frozen prescribed branch.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import mpmath as mp

try:
    from .lei_ren_part1_paper_axial_dual import AxialDual
    from .lei_ren_part1_paper_pressure_width_axial_bridge import (
        AxialPressureWidthExitBridge,
    )
    from .lei_ren_part1_paper_pressure_width_axial_bridge_fixture import (
        _bundle,
        _make_formal_core,
    )
    from .lei_ren_part1_paper_pressure_width_axial_comparison import (
        AxialPressureWidthComparison,
    )
    from .lei_ren_part1_paper_pressure_width_continuation import (
        PressureWidthExitContinuation,
    )
    from .lei_ren_part1_paper_pressure_width_jet import PressureWidthJet
except (ImportError, ValueError):
    from lei_ren_part1_paper_axial_dual import AxialDual
    from lei_ren_part1_paper_pressure_width_axial_bridge import (
        AxialPressureWidthExitBridge,
    )
    from lei_ren_part1_paper_pressure_width_axial_bridge_fixture import (
        _bundle,
        _make_formal_core,
    )
    from lei_ren_part1_paper_pressure_width_axial_comparison import (
        AxialPressureWidthComparison,
    )
    from lei_ren_part1_paper_pressure_width_continuation import (
        PressureWidthExitContinuation,
    )
    from lei_ren_part1_paper_pressure_width_jet import PressureWidthJet


PRECISION = 100
PRESSURE_ORDER = 3
WIDTH_ORDER = 2
TRANSITION_STEPS = 4
BRIDGE_STEPS = 8
RK_STEPS = 512
TARGET_RADIUS = ".8"
ERROR_THRESHOLD = mp.mpf("1e-9")


def _scalar(value: Any, *, pressure: Any = 1, width: Any = 1) -> mp.mpf:
    """Materialize an axial dual or pressure/width jet at one retained point."""

    if isinstance(value, AxialDual):
        value = value.value
    if isinstance(value, PressureWidthJet):
        return value.evaluate(pressure=pressure, width=width)
    return mp.mpf(str(value))


def _ring_component(value: Any, pressure: int, width: int) -> mp.mpf:
    if isinstance(value, AxialDual):
        value = value.value
    if isinstance(value, PressureWidthJet):
        return value.component(pressure, width)
    return mp.mpf(str(value)) if (pressure, width) == (0, 0) else mp.mpf(0)


def _add_state(left: list[mp.mpf], right: list[mp.mpf], scale: mp.mpf) -> list[mp.mpf]:
    return [a + scale * b for a, b in zip(left, right)]


def _rk4(
    functions: dict[str, Any],
    *,
    rb: mp.mpf,
    epsilon: mp.mpf,
    f_start: mp.mpf,
    f_bar: mp.mpf,
    initial_pressure: mp.mpf,
    initial_moments: dict[str, mp.mpf],
    initial_uz: mp.mpf,
    target_radius: mp.mpf,
    steps: int,
) -> dict[str, mp.mpf | dict[str, mp.mpf]]:
    """Independent scalar RK4 replay of the prescribed continuation."""

    moment_names = ("theta", "z", "theta_z", "z_theta", "p")
    y_target = mp.log(target_radius / rb)
    state = [mp.mpf(0), initial_uz] + [initial_moments[name] for name in moment_names]
    step = y_target / steps

    def rhs(y: mp.mpf, current: list[mp.mpf]) -> list[mp.mpf]:
        g, uz = current[:2]
        radius = rb * mp.exp(y)
        field = f_start * mp.exp(g)
        D = _scalar(functions["D"].evaluate(y))
        J = _scalar(functions["J"].evaluate(y))
        return [
            -epsilon * D / 2,
            -epsilon * (f_start / f_bar) * mp.exp(g) * J,
            2 * radius * radius * field,
            radius * uz,
            2 * radius * radius * field * uz,
            radius * uz * uz - radius * radius * field * field,
            radius * field * field,
        ]

    y = mp.mpf(0)
    for _ in range(steps):
        k1 = rhs(y, state)
        k2 = rhs(y + step / 2, _add_state(state, k1, step / 2))
        k3 = rhs(y + step / 2, _add_state(state, k2, step / 2))
        k4 = rhs(y + step, _add_state(state, k3, step))
        state = [
            value + step * (a + 2 * b + 2 * c + d) / 6
            for value, a, b, c, d in zip(state, k1, k2, k3, k4)
        ]
        y += step

    moments = dict(zip(moment_names, state[2:]))
    return {
        "g": state[0],
        "Uz": state[1],
        "F": f_start * mp.exp(state[0]),
        "moments": moments,
        "P": initial_pressure + state[-1] - initial_moments["p"],
        "y": y_target,
    }


def run_fixture(
    *,
    precision: int = PRECISION,
    target_radius: Any = TARGET_RADIUS,
    steps: int = RK_STEPS,
) -> dict[str, Any]:
    """Run the resolved finite-ring continuation check and save JSON evidence."""

    precision = max(70, int(precision))
    if int(steps) != steps or steps < 8:
        raise ValueError("steps must be an integer at least 8")

    with mp.workdps(precision):
        target = mp.mpf(str(target_radius))
        z0 = mp.mpf(".3")
        component, _formal = _make_formal_core(z0=z0, precision=precision)
        axis = component.axis
        comparison = AxialPressureWidthComparison(
            _bundle(axis, component),
            h_b=mp.mpf("1e-5"),
            pressure_order=PRESSURE_ORDER,
            width_order=WIDTH_ORDER,
            transition_steps=TRANSITION_STEPS,
        )
        bridge = AxialPressureWidthExitBridge(comparison, steps=BRIDGE_STEPS)
        continuation = PressureWidthExitContinuation(bridge)
        functions = continuation.functions(z0)

        rb = _scalar(functions["Rb"])
        epsilon = _scalar(comparison.W)
        start = functions["start"]
        f_start = _scalar(start["F"])
        f_bar = _scalar(continuation.auxiliary.endpoint(z0)["F"])
        initial_pressure = _scalar(start["P"])
        initial_uz = _scalar(start["Uz"])
        initial_moments = {
            name: _scalar(start["moments"][name])
            for name in ("theta", "z", "theta_z", "z_theta", "p")
        }
        independent = _rk4(
            functions,
            rb=rb,
            epsilon=epsilon,
            f_start=f_start,
            f_bar=f_bar,
            initial_pressure=initial_pressure,
            initial_moments=initial_moments,
            initial_uz=initial_uz,
            target_radius=target,
            steps=steps,
        )
        produced = continuation.evaluate_R(target, z0)

        errors: dict[str, mp.mpf] = {}
        for name in ("F", "Uz", "P"):
            actual = _scalar(produced[name])
            expected = independent[name]
            errors[name] = abs(actual - expected) / max(
                mp.mpf(1), abs(actual), abs(expected)
            )
        for name in ("theta", "z", "theta_z", "z_theta", "p"):
            actual = _scalar(produced["moments"][name])
            expected = independent["moments"][name]
            errors[f"moment:{name}"] = abs(actual - expected) / max(
                mp.mpf(1), abs(actual), abs(expected)
            )
        maximum = max(errors.values())
        analytic_defects={name:max((abs(v) for v in produced[name].atoms.values()),default=mp.mpf(0))
            for name in ('analytic_divergence_numerator','moment_radial_identity_defect','Uz_y_driver_defect')}
        if max(analytic_defects.values())>=mp.mpf('1e-80'):
            raise AssertionError(analytic_defects)

        start_f_ring = start["F"]
        start_uz_ring = start["Uz"]
        produced_f_width_change = _ring_component(produced["F"], 0, 1) - _ring_component(
            start_f_ring, 0, 1
        )
        produced_uz_width_change = _ring_component(produced["Uz"], 0, 1) - _ring_component(
            start_uz_ring, 0, 1
        )
        if maximum >= ERROR_THRESHOLD:
            raise AssertionError(
                f"continuation scalar RK error {mp.nstr(maximum, 20)}"
            )
        if produced_f_width_change == 0 or produced_uz_width_change == 0:
            raise AssertionError("positive-width continuation field change vanished")

        report = {
            "passed": True,
            "target_radius": mp.nstr(target, 30),
            "endpoint_radius": mp.nstr(rb, 30),
            "epsilon": mp.nstr(epsilon, 30),
            "rk_steps": steps,
            "pressure_order": PRESSURE_ORDER,
            "width_order": WIDTH_ORDER,
            "maximum_scaled_difference": mp.nstr(maximum, 30),
            "analytic_primitive_defects":{k:mp.nstr(v,30) for k,v in analytic_defects.items()},
            "independent_Cartesian_divergence_verified":False,
            "errors": {key: mp.nstr(value, 20) for key, value in sorted(errors.items())},
            "positive_width_F_change": mp.nstr(produced_f_width_change, 30),
            "positive_width_Uz_change": mp.nstr(produced_uz_width_change, 30),
            "continuation_metadata": continuation.metadata(),
            "independent_replay": "scalar RK4 in y=log(R/Rb), materialized at pressure=width=1",
            "pressure_width_truncation_remainder_enclosed": False,
            "rk_remainder_enclosed": False,
            "global_field_installed": False,
            "legacy_frozen_prescribed_comparison_used": False,
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
