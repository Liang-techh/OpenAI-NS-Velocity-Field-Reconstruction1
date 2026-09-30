"""Threshold radii for the saved-core necessary inequalities.

Using the saved ``load_core()`` profile, this module solves on the 33-point
grid ``Z = linspace(-1, 1, 33)`` the conditional angular equation

``F_core(.005,Z) R^2 = Mtheta_heat_target(R,Z)``

and the heat endpoint equation
``F_heat(R,Z) = F_core(.005,Z)``.  Roots are searched on ``1 <= R <= 4096``;
when an inequality already holds at ``R=1``, the reported threshold is the
lower search endpoint.  The receipt also records the absolute safety gaps at
``R=1024`` and ``R=2048``.

These are necessary-condition thresholds for the saved independent core.
They do not establish sufficient angular closure, admissibility, PDE
validity, global energy, or a chosen final outer parameter.

Source pinned to Z. Lei and X. Ren, arXiv:2609.35406v1, Section 2.5.  Heat
angular targets use the local finite-tail implementation tied to Sections
4.23 and 5.1.
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path
from typing import Any, Callable

import numpy as np
from scipy.optimize import brentq


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
CORE_EXIT_RADIUS = 0.005
SEARCH_LOWER = 1.0
SEARCH_UPPER = 4096.0
TARGET_SERIES_ORDER = 12
MARGIN_RADII = (1024.0, 2048.0)


def _finite(value: Any, name: str) -> float:
    """Convert one value to a finite scalar."""

    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be finite") from exc
    if not math.isfinite(result):
        raise ValueError(f"{name} must be finite")
    return result


def _root_threshold(
    gap: Callable[[float], float],
    *,
    lower: float = SEARCH_LOWER,
    upper: float = SEARCH_UPPER,
) -> dict[str, Any]:
    """Find the first bracketed sign change in the prescribed interval."""

    lower_gap = _finite(gap(lower), "gap_at_lower")
    upper_gap = _finite(gap(upper), "gap_at_upper")
    if lower_gap >= 0.0:
        return {
            "status": "already_satisfied_at_lower",
            "critical_radius_in_search_interval": float(lower),
            "R_lower": float(lower),
            "R_upper": float(lower),
            "gap_at_lower": lower_gap,
            "gap_at_upper": lower_gap,
            "lower_clipped": True,
        }
    if upper_gap < 0.0:
        return {
            "status": "not_bracketed_in_search_interval",
            "critical_radius_in_search_interval": None,
            "R_lower": float(lower),
            "R_upper": float(upper),
            "gap_at_lower": lower_gap,
            "gap_at_upper": upper_gap,
            "lower_clipped": False,
        }
    root = brentq(gap, lower, upper, xtol=1.0e-9, rtol=1.0e-12)
    root_gap = _finite(gap(root), "gap_at_root")
    return {
        "status": "bracketed",
        "critical_radius_in_search_interval": float(root),
        "R_lower": float(lower),
        "R_upper": float(upper),
        "gap_at_lower": lower_gap,
        "gap_at_upper": upper_gap,
        "gap_at_root": root_gap,
        "lower_clipped": False,
    }


def _target(heat: Any, radius: float, z: float) -> float:
    """Evaluate the stable heat angular target at one finite radius."""

    return _finite(
        heat_tail_moments(
            radius, z, heat, series_order=TARGET_SERIES_ORDER
        )["angular_target"],
        "Mtheta_heat_target",
    )


def _row(core: Any, joined: Any, z: float) -> dict[str, Any]:
    """Solve both threshold equations for one axial coordinate."""

    f_axis = _finite(core.F(0.0, z), "F_axis")
    f_core_exit = _finite(core.F(CORE_EXIT_RADIUS, z), "F_core_exit")

    def angular_gap(radius: float) -> float:
        return f_core_exit * radius * radius - _target(joined.heat, radius, z)

    def endpoint_gap(radius: float) -> float:
        return f_core_exit - _finite(
            joined.F_heat(radius, z), "F_heat_endpoint"
        )

    angular = _root_threshold(angular_gap)
    endpoint = _root_threshold(endpoint_gap)
    return {
        "Z": float(z),
        "F_axis": f_axis,
        "F_core_exit_at_R=.005": f_core_exit,
        "angular_threshold": angular,
        "endpoint_threshold": endpoint,
    }


def _margin_summary(core: Any, joined: Any, radius: float, z_values: np.ndarray) -> dict[str, Any]:
    """Report minimum angular/endpoint gaps at one requested join."""

    rows = []
    for z_value in z_values:
        z = float(z_value)
        f_core_exit = _finite(core.F(CORE_EXIT_RADIUS, z), "F_core_exit")
        target = _target(joined.heat, radius, z)
        f_heat_endpoint = _finite(joined.F_heat(radius, z), "F_heat_endpoint")
        angular_gap = f_core_exit * radius * radius - target
        endpoint_gap = f_core_exit - f_heat_endpoint
        cauchy = 3.0 * target * target / (4.0 * radius**3)
        rows.append(
            {
                "Z": z,
                "angular_gap": float(angular_gap),
                "endpoint_gap": float(endpoint_gap),
                "cauchy_minimum_Mp": float(cauchy),
            }
        )
    angular_worst = min(rows, key=lambda row: row["angular_gap"])
    endpoint_worst = min(rows, key=lambda row: row["endpoint_gap"])
    cauchy_worst = max(rows, key=lambda row: row["cauchy_minimum_Mp"])
    return {
        "R_join": float(radius),
        "minimum_angular_safety_gap": float(angular_worst["angular_gap"]),
        "minimum_angular_safety_gap_Z": float(angular_worst["Z"]),
        "minimum_endpoint_safety_gap": float(endpoint_worst["endpoint_gap"]),
        "minimum_endpoint_safety_gap_Z": float(endpoint_worst["Z"]),
        "maximum_cauchy_minimum_Mp": float(cauchy_worst["cauchy_minimum_Mp"]),
        "maximum_cauchy_minimum_Mp_Z": float(cauchy_worst["Z"]),
        "angular_condition_satisfied_on_grid": bool(
            all(row["angular_gap"] >= 0.0 for row in rows)
        ),
        "endpoint_condition_satisfied_on_grid": bool(
            all(row["endpoint_gap"] > 0.0 for row in rows)
        ),
    }


def run() -> dict[str, Any]:
    """Compute the 33-point threshold receipt and write its JSON artifact."""

    core, joined = load_core()
    z_values = np.linspace(-1.0, 1.0, 33)
    rows = [_row(core, joined, float(z)) for z in z_values]
    angular_values = [
        row["angular_threshold"]["critical_radius_in_search_interval"]
        for row in rows
        if row["angular_threshold"]["critical_radius_in_search_interval"] is not None
    ]
    endpoint_values = [
        row["endpoint_threshold"]["critical_radius_in_search_interval"]
        for row in rows
        if row["endpoint_threshold"]["critical_radius_in_search_interval"] is not None
    ]
    angular_max = max(angular_values) if angular_values else None
    endpoint_max = max(endpoint_values) if endpoint_values else None
    angular_worst = max(
        rows,
        key=lambda row: row["angular_threshold"][
            "critical_radius_in_search_interval"
        ]
        if row["angular_threshold"]["critical_radius_in_search_interval"] is not None
        else float("-inf"),
    )
    endpoint_worst = max(
        rows,
        key=lambda row: row["endpoint_threshold"][
            "critical_radius_in_search_interval"
        ]
        if row["endpoint_threshold"]["critical_radius_in_search_interval"] is not None
        else float("-inf"),
    )
    global_threshold_candidates = [
        (float(angular_max), "angular", float(angular_worst["Z"]))
        for _ in [0]
        if angular_max is not None
    ] + [
        (float(endpoint_max), "endpoint", float(endpoint_worst["Z"]))
        for _ in [0]
        if endpoint_max is not None
    ]
    global_threshold, global_metric, global_z = max(
        global_threshold_candidates, key=lambda item: item[0]
    )
    margins = [
        _margin_summary(core, joined, radius, z_values)
        for radius in (1024.0, 2048.0)
    ]
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
        "search_interval": [SEARCH_LOWER, SEARCH_UPPER],
        "target_series_order": TARGET_SERIES_ORDER,
        "rows": rows,
        "angular_threshold_max_radius": angular_max,
        "angular_threshold_worst_Z": float(angular_worst["Z"]),
        "endpoint_threshold_max_radius": endpoint_max,
        "endpoint_threshold_worst_Z": float(endpoint_worst["Z"]),
        "global_max_threshold_radius": float(global_threshold),
        "global_max_threshold_metric": global_metric,
        "global_max_threshold_Z": global_z,
        "safety_margins_at_1024_2048": margins,
        "necessary_conditions_only": True,
        "sufficient_angular_closure": False,
        "admissibility_validated": False,
        "pde_validated": False,
        "chosen_final_parameter": False,
        "global_impossibility_claim": False,
        "scope": (
            "Read-only threshold quantification for the saved core over the "
            "full 33-point axial grid. Roots and safety gaps only measure the "
            "listed necessary inequalities; they do not establish sufficient "
            "closure, admissibility, PDE, global energy, or a final outer "
            "parameter."
        ),
    }
    output = Path(__file__).with_suffix(".json")
    output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report))
    return report


if __name__ == "__main__":
    run()


__all__ = ["run", "SOURCE", "SOURCE_SECTION", "SOURCE_VERSION"]
