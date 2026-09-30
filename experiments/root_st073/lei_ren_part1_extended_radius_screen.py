"""Read-only extended-radius necessary-condition screen.

The saved finite core is loaded through
``lei_ren_part1_pressure_core.load_core()``.  For the full grid
``Z = linspace(-1, 1, 33)``, this module evaluates the exact supplied heat
angular target at four candidate joins and compares it with the conditional
bound

``Mtheta_heat_target(Rb,Z) <= F_core(.005,Z) * Rb^2``.

It also checks the endpoint condition ``F_heat(Rb,Z) < F_core(.005,Z)`` and
records the Cauchy minimum pressure moment
``3 Mtheta_heat_target^2/(4 Rb^3)``.  ``F_axis = F_core(0,Z)`` is retained in
every row so that the small endpoint value at ``Z = +/-1`` is explicit.

This is a finite parameter screen for the saved independent core.  Passing
these necessary inequalities is not sufficient for angular closure,
admissibility, a PDE, global energy, the paper exterior/collar, or a chosen
final parameter.  No profile is modified.

Source pinned to Z. Lei and X. Ren, arXiv:2609.35406v1, Section 2.5.  Heat
angular targets use the local finite-tail implementation tied to Sections
4.23 and 5.1.
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
LOCAL = ROOT / "experiments" / "root_st073"
for path in (SRC, LOCAL):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from lei_ren_part1_heat_moments import heat_tail_moments  # noqa: E402
from lei_ren_part1_pressure_core import load_core  # noqa: E402


SOURCE = "https://arxiv.org/html/2609.35406v1"
SOURCE_VERSION = "2609.35406v1"
SOURCE_SECTION = "Section 2.5"
RADIUS_SCREEN = (64.0, 256.0, 1024.0, 2048.0)
CORE_EXIT_RADIUS = 0.005
TARGET_SERIES_ORDER = 12


def _finite(value: Any, name: str) -> float:
    """Convert one value to a finite real scalar."""

    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be finite") from exc
    if not math.isfinite(result):
        raise ValueError(f"{name} must be finite")
    return result


def _row(core: Any, joined: Any, radius: float, z: float) -> dict[str, Any]:
    """Evaluate one ``(Rjoin,Z)`` screen row."""

    f_axis = _finite(core.F(0.0, z), "F_axis")
    f_core_exit = _finite(core.F(CORE_EXIT_RADIUS, z), "F_core_exit")
    target = _finite(
        heat_tail_moments(
            radius, z, joined.heat, series_order=TARGET_SERIES_ORDER
        )["angular_target"],
        "Mtheta_heat_target",
    )
    f_heat_endpoint = _finite(joined.F_heat(radius, z), "F_heat_endpoint")
    angular_bound = f_core_exit * radius * radius
    cauchy_lower = 3.0 * target * target / (4.0 * radius**3)
    return {
        "Z": float(z),
        "F_axis": f_axis,
        "F_core_exit_at_R=.005": f_core_exit,
        "Mtheta_heat_target": target,
        "monotone_angular_bound_Fcore_exit_R2": float(angular_bound),
        "angular_bound_minus_target": float(angular_bound - target),
        "angular_bound_satisfied": bool(target <= angular_bound),
        "F_heat_endpoint": f_heat_endpoint,
        "endpoint_minus_Fcore_exit": float(f_core_exit - f_heat_endpoint),
        "endpoint_condition_satisfied": bool(f_heat_endpoint < f_core_exit),
        "cauchy_minimum_Mp_for_heat_target": float(cauchy_lower),
    }


def _summary(radius: float, rows: list[dict[str, Any]]) -> dict[str, Any]:
    """Summarize pass/fail and worst axial locations for one radius."""

    angular_worst = min(rows, key=lambda row: row["angular_bound_minus_target"])
    endpoint_worst = min(rows, key=lambda row: row["endpoint_minus_Fcore_exit"])
    cauchy_worst = max(
        rows, key=lambda row: row["cauchy_minimum_Mp_for_heat_target"]
    )
    angular_ok = all(bool(row["angular_bound_satisfied"]) for row in rows)
    endpoint_ok = all(bool(row["endpoint_condition_satisfied"]) for row in rows)
    return {
        "R_join": float(radius),
        "all_angular_bounds_satisfied": angular_ok,
        "all_endpoint_conditions_satisfied": endpoint_ok,
        "all_screen_conditions_satisfied": bool(angular_ok and endpoint_ok),
        "worst_angular_Z": float(angular_worst["Z"]),
        "worst_angular_bound_minus_target": float(
            angular_worst["angular_bound_minus_target"]
        ),
        "worst_endpoint_Z": float(endpoint_worst["Z"]),
        "worst_endpoint_minus_Fcore_exit": float(
            endpoint_worst["endpoint_minus_Fcore_exit"]
        ),
        "max_cauchy_minimum_Mp": float(
            cauchy_worst["cauchy_minimum_Mp_for_heat_target"]
        ),
        "max_cauchy_minimum_Mp_Z": float(cauchy_worst["Z"]),
    }


def run() -> dict[str, Any]:
    """Run the fixed-radius, full-axial-grid screen and write its receipt."""

    core, joined = load_core()
    z_values = np.linspace(-1.0, 1.0, 33)
    screens: list[dict[str, Any]] = []
    summaries: list[dict[str, Any]] = []
    for radius in RADIUS_SCREEN:
        rows = [_row(core, joined, radius, float(z)) for z in z_values]
        screens.append({"R_join": radius, "rows": rows})
        summaries.append(_summary(radius, rows))

    passing = [
        row["R_join"] for row in summaries if row["all_screen_conditions_satisfied"]
    ]
    smallest = min(passing) if passing else None
    report: dict[str, Any] = {
        "status": "completed",
        "source": SOURCE,
        "source_version": SOURCE_VERSION,
        "source_section": SOURCE_SECTION,
        "profile_source": "lei_ren_part1_pressure_core.load_core()",
        "profile_json": str(
            Path(__file__).with_name("lei_ren_part1_pressure_core.json")
        ),
        "core_parameters": core.reference.metadata()["parameters"],
        "heat": {"h": float(joined.heat.h), "c_inf": float(joined.heat.c_inf)},
        "Z_grid": z_values.tolist(),
        "core_exit_radius": CORE_EXIT_RADIUS,
        "target_series_order": TARGET_SERIES_ORDER,
        "tested_R_join": list(RADIUS_SCREEN),
        "screen_summaries": summaries,
        "smallest_tested_radius_satisfying_all_necessary_screen_conditions": smallest,
        "screens": screens,
        "condition_definitions": {
            "angular": "Mtheta_heat_target(Rjoin,Z) <= F_core(.005,Z)*Rjoin^2",
            "endpoint": "F_heat(Rjoin,Z) < F_core(.005,Z)",
            "pressure": "reported lower bound 3*Mtheta_heat_target^2/(4*Rjoin^3); no pressure budget pass/fail is asserted",
        },
        "worst_Z_interpretation": (
            "Worst axial locations are reported separately for angular gap, "
            "endpoint gap, and the maximum Cauchy lower bound."
        ),
        "necessary_conditions_only": True,
        "sufficient_angular_closure": False,
        "admissibility_validated": False,
        "pde_validated": False,
        "chosen_final_parameter": False,
        "global_impossibility_claim": False,
        "scope": (
            "Read-only full-Z screen for four extended joins using the saved "
            "core axis and core-exit values. Passing the listed inequalities "
            "does not establish angular closure, admissibility, PDE, global "
            "energy, or a final outer/collar parameter."
        ),
    }
    output = Path(__file__).with_suffix(".json")
    output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report))
    return report


if __name__ == "__main__":
    run()


__all__ = ["run", "SOURCE", "SOURCE_SECTION", "SOURCE_VERSION"]
