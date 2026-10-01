"""Small protocol check for the optional analytic exit-driver hooks.

The check uses scalar fake drivers only.  It verifies the actual ``rhs`` and
switch ``_rhs``/``_drivers_at`` assemblies, while the provider raises if it is
ever called at an off-center axial coordinate.  It also exercises the old
finite-difference branches through locally replaced scalar driver functions;
no production source or core data are generated here.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import mpmath as mp


HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from lei_ren_part1_paper_exit_continuation import ExitContinuation
from lei_ren_part1_paper_exit_switches import ExitSwitches
from lei_ren_part1_paper_exit_tangents import ExitTangents


class _Core:
    Lambda = mp.mpf("4")
    delta = mp.mpf("0")
    precision = 100


class _Comparison:
    precision = 100
    h_b = mp.mpf(".01")
    core = _Core()


class _ExitProvider:
    def __init__(self, center: mp.mpf) -> None:
        self.center = center
        self.calls: list[tuple[mp.mpf, mp.mpf]] = []

    def exit_driver(self, y: mp.mpf, Z: mp.mpf) -> dict[str, mp.mpf]:
        y = mp.mpf(y)
        Z = mp.mpf(Z)
        self.calls.append((y, Z))
        if Z != self.center:
            raise AssertionError("analytic exit hook received an off-center Z")
        return {"A": mp.mpf(".2"), "B": mp.mpf(".3"),
                "A_Z": mp.mpf(".4"), "B_Z": mp.mpf(".5")}


class _SwitchProvider:
    def __init__(self, center: mp.mpf) -> None:
        self.center = center
        self.calls: list[tuple[mp.mpf, mp.mpf, mp.mpf]] = []

    def switch_driver(self, x: mp.mpf, Z: mp.mpf, blend: mp.mpf) -> dict[str, mp.mpf]:
        x = mp.mpf(x)
        Z = mp.mpf(Z)
        blend = mp.mpf(blend)
        self.calls.append((x, Z, blend))
        if Z != self.center:
            raise AssertionError("analytic switch hook received an off-center Z")
        return {
            "D": mp.mpf("2"),
            "D_Z": mp.mpf(".6"),
            "I_z": mp.mpf("3"),
            "I_z_Z": mp.mpf(".7"),
            "Fbar": mp.mpf("5"),
            "Fbar_Z": mp.mpf(".8"),
            "F0": mp.mpf("7"),
            "F0_Z": mp.mpf(".9"),
            "barE": mp.mpf("1.1"),
        }


def _assert_close(actual: mp.mpf, expected: mp.mpf, label: str) -> None:
    if abs(actual - expected) > mp.mpf("1e-80"):
        raise AssertionError(f"{label}: {actual} != {expected}")


def run() -> dict[str, object]:
    with mp.workdps(100):
        center = mp.mpf(".3")
        comparison = _Comparison()
        exit_provider = _ExitProvider(center)
        tangent = ExitTangents(
            comparison,
            epsilon=".01",
            steps=8,
            analytic_driver_provider=exit_provider,
        )
        exit_rhs = tangent.rhs(0, [mp.mpf(0)] * 16, center)
        _assert_close(exit_rhs[0], mp.mpf(".2"), "exit A")
        _assert_close(exit_rhs[1], mp.mpf(".3"), "exit B")
        _assert_close(exit_rhs[8], mp.mpf(".4"), "exit A_Z")
        _assert_close(exit_rhs[9], mp.mpf(".5"), "exit B_Z")

        continuation = ExitContinuation(tangent)
        switch_provider = _SwitchProvider(center)
        switches = ExitSwitches(
            continuation,
            steps=8,
            analytic_driver_provider=switch_provider,
        )
        switch_rhs = switches._rhs(0, [mp.mpf(0)] * 16, center, 1)
        _assert_close(switch_rhs[0], -mp.mpf(".01"), "switch a")
        _assert_close(switch_rhs[1], mp.mpf(".01") * 0 +
                      (-switches.epsilon * mp.sqrt(mp.mpf("50")) *
                       (mp.mpf("7") / 5) * 3), "switch B")
        _assert_close(switch_rhs[8], -switches.epsilon * mp.mpf(".6") / 2,
                      "switch a_Z")
        switch_driver = switches._drivers_at(0, center, "switch_1")
        _assert_close(switch_driver["a_Z"], switches.epsilon * mp.mpf(".6"),
                      "drivers_at a_Z")
        if not switch_driver["analytic_driver_used"]:
            raise AssertionError("switch metadata did not mark analytic hook")

        # Exercise fallback branch selection without invoking a real core.
        fallback_tangent = ExitTangents(comparison, epsilon=".01", steps=8)
        fallback_tangent._normalized_driver = lambda y, z: {
            "A": z, "B": 2 * z,
        }
        fallback = fallback_tangent._driver_tangent(0, center)
        if fallback["analytic_driver_used"]:
            raise AssertionError("tangent fallback unexpectedly used analytic hook")
        if "centered" not in fallback["derivative_method"]:
            raise AssertionError("tangent fallback method marker changed")

        fallback_switches = ExitSwitches(continuation, steps=8)
        fallback_switches._bar_driver = lambda x, z, blend: {
            "R": mp.mpf(100), "D": z, "B": z, "blend": blend,
            "F0": mp.mpf(7), "F0_Z": mp.mpf(".9"),
            "bar": {"D": z, "E": mp.mpf(1), "I_z": z,
                    "F": mp.mpf(5), "FZ": mp.mpf(".8")},
        }
        fallback_switches.comparison.evaluate = lambda y, z: {
            "D": z, "I_z": z, "F": mp.mpf(5), "FZ": mp.mpf(".8"),
        }
        fallback_switch = fallback_switches._driver_tangent(0, center, 1)
        if fallback_switch["analytic_driver_used"]:
            raise AssertionError("switch fallback unexpectedly used analytic hook")
        if "centered" not in fallback_switch["derivative_method"]:
            raise AssertionError("switch fallback method marker changed")

        result = {
            "protocol": {
                "exit": "provider.exit_driver(y, Z) -> A,B,A_Z,B_Z",
                "switch": "provider.switch_driver(x, Z, blend) -> D,D_Z,I_z,I_z_Z,Fbar,Fbar_Z,F0,F0_Z plus optional bar/barE",
            },
            "center": "0.3",
            "exit_provider_calls": len(exit_provider.calls),
            "switch_provider_calls": len(switch_provider.calls),
            "exit_rhs": [mp.nstr(v, 30) for v in exit_rhs[:10]],
            "switch_rhs": [mp.nstr(v, 30) for v in switch_rhs[:10]],
            "analytic_exit_used": True,
            "analytic_switch_used": True,
            "fallback_exit_used": True,
            "fallback_switch_used": True,
            "off_center_provider_calls": 0,
            "source_scope": "Hook plumbing only; no core, pressure, transition, or cone certification.",
        }
    Path(__file__).with_suffix(".json").write_text(
        json.dumps(result, indent=2) + "\n", encoding="utf-8"
    )
    return result


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))
