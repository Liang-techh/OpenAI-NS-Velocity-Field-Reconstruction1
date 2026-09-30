"""Finite angular and pressure budget diagnostics for the saved joined core.

This diagnostic reads the saved ``load_core()`` result from
``lei_ren_part1_pressure_core``.  It computes the finite inner moments

``Mtheta_inner = 2 integral_0^Rjoin rho F(rho,Z) d rho`` and
``Mp_inner = integral_0^Rjoin F(rho,Z)^2 d rho``

with Gauss--Legendre rules split at the saved core/join interface.  The
finite inner angular value is compared with the exact supplied heat-tail
target at the same ``Rjoin``.  It also records two necessary-condition
diagnostics: the monotone-cone estimate
``Mtheta <= F_axis Rjoin^2`` when ``F_R <= 0`` and ``F_axis > 0``, and the
Cauchy lower bound ``Mp >= 3 Mtheta^2/(4 Rjoin^3)``.

The profile and heat target are independent finite numerical objects.  These
comparisons do not prove impossibility, admissibility, PDE validity, global
energy, or the Lei--Ren exterior/collar construction.  The rough radius
roots only locate where the fixed saved ``F_axis`` could first meet the
monotone-cone inequality; the inequality is necessary and not sufficient.

Source pinned to Z. Lei and X. Ren, arXiv:2609.35406v1, Section 2.5.  The
heat target is evaluated by the local Section 4.23/5.1 tail implementation.
"""

from __future__ import annotations

import json
import math
import operator
import sys
from pathlib import Path
from typing import Any, Callable

import numpy as np
from numpy.polynomial.legendre import leggauss
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


def _finite(value: Any, name: str) -> float:
    """Return one finite scalar or raise a diagnostic input error."""

    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be finite") from exc
    if not math.isfinite(result):
        raise ValueError(f"{name} must be finite")
    return result


def _positive_order(value: Any) -> int:
    """Validate a positive integer Gauss order."""

    if isinstance(value, (bool, np.bool_)):
        raise TypeError("n must be a positive integer")
    try:
        result = operator.index(value)
    except TypeError as exc:
        raise TypeError("n must be a positive integer") from exc
    if result <= 0:
        raise ValueError("n must be a positive integer")
    return int(result)


def _gauss_split(
    function: Callable[[float], float],
    left: float,
    right: float,
    split: float,
    *,
    n: Any,
) -> float:
    """Integrate a scalar function by Gauss rules on core and blend pieces."""

    left = _finite(left, "left")
    right = _finite(right, "right")
    split = _finite(split, "split")
    order = _positive_order(n)
    if left < 0.0 or not left < split < right:
        raise ValueError("require 0 <= left < split < right")
    nodes, weights = leggauss(order)
    total = 0.0
    for a, b in ((left, split), (split, right)):
        points = (b - a) * (nodes + 1.0) / 2.0 + a
        values = np.empty(order, dtype=float)
        for index, point in enumerate(points):
            try:
                value = float(function(float(point)))
            except (TypeError, ValueError, IndexError) as exc:
                raise ValueError("quadrature integrand returned an invalid value") from exc
            if not math.isfinite(value):
                raise ValueError("quadrature integrand returned a non-finite value")
            values[index] = value
        total += float((b - a) * np.dot(weights, values) / 2.0)
    if not math.isfinite(total):
        raise ArithmeticError("split Gauss quadrature returned a non-finite value")
    return float(total)


def _inner_moments(joined: Any, z: float, *, n: int) -> tuple[float, float]:
    """Return finite ``(Mtheta, Mp)`` for one saved joined profile."""

    r_join = _finite(joined.R_join, "R_join")
    r_core = _finite(joined.R_core, "R_core")
    mtheta = _gauss_split(
        lambda rho: 2.0 * rho * float(joined.F(rho, z)),
        0.0,
        r_join,
        r_core,
        n=n,
    )
    mp = _gauss_split(
        lambda rho: float(joined.F(rho, z)) ** 2,
        0.0,
        r_join,
        r_core,
        n=n,
    )
    return mtheta, mp


def _sample_monotone(joined: Any, z: float, *, count: int = 129) -> dict[str, Any]:
    """Record a finite-grid monotonicity diagnostic for the saved profile."""

    if count < 2:
        raise ValueError("count must be at least 2")
    r_join = float(joined.R_join)
    radii = np.linspace(0.0, r_join, count)
    values = np.asarray([float(joined.F(float(radius), z)) for radius in radii])
    increments = np.diff(values)
    max_increase = float(np.max(increments)) if increments.size else 0.0
    return {
        "sample_count": int(count),
        "max_positive_increment": max(0.0, max_increase),
        "sampled_nonincreasing": bool(max_increase <= 1.0e-10),
        "F_axis": float(values[0]),
        "F_at_Rjoin": float(values[-1]),
    }


def _rough_minimum_join(
    heat: Any,
    z: float,
    f_axis: float,
    r_start: float,
    *,
    series_order: int = 8,
    max_doublings: int = 12,
) -> dict[str, Any]:
    """Find a first bracketed radius satisfying ``target <= F_axis R^2``."""

    def gap(radius: float) -> float:
        target = heat_tail_moments(
            radius, z, heat, series_order=series_order
        )["angular_target"]
        return f_axis * radius * radius - target

    if f_axis <= 0.0:
        return {
            "status": "not_applicable_Faxis_nonpositive",
            "R_lower": float(r_start),
            "R_upper": float(r_start),
            "minimum_R_rough": None,
            "gap_at_lower": None,
            "gap_at_upper": None,
            "series_order": int(series_order),
        }

    lower = float(r_start)
    lower_gap = float(gap(lower))
    if lower_gap >= 0.0:
        return {
            "status": "already_satisfied_at_Rjoin",
            "R_lower": lower,
            "R_upper": lower,
            "minimum_R_rough": lower,
            "gap_at_lower": lower_gap,
            "gap_at_upper": lower_gap,
            "series_order": int(series_order),
        }

    upper = lower
    upper_gap = lower_gap
    for _ in range(max_doublings):
        upper *= 2.0
        upper_gap = float(gap(upper))
        if upper_gap >= 0.0:
            root = brentq(gap, lower, upper, xtol=1.0e-10, rtol=1.0e-12)
            return {
                "status": "bracketed",
                "R_lower": lower,
                "R_upper": upper,
                "minimum_R_rough": float(root),
                "gap_at_lower": lower_gap,
                "gap_at_upper": upper_gap,
                "series_order": int(series_order),
            }
        lower = upper
        lower_gap = upper_gap
    return {
        "status": "not_bracketed",
        "R_lower": float(r_start),
        "R_upper": upper,
        "minimum_R_rough": None,
        "gap_at_lower": float(gap(r_start)),
        "gap_at_upper": upper_gap,
        "series_order": int(series_order),
    }


def run(*, n: int = 64) -> dict[str, Any]:
    """Compute the nine-point angular/pressure compatibility receipt."""

    core, joined = load_core()
    z_values = np.linspace(-0.5, 0.5, 9)
    r_core = float(joined.R_core)
    r_join = float(joined.R_join)
    heat = joined.heat

    rows: list[dict[str, Any]] = []
    for z_value in z_values:
        z = float(z_value)
        mtheta, mp = _inner_moments(joined, z, n=n)
        mtheta_32, mp_32 = _inner_moments(joined, z, n=32)
        target = heat_tail_moments(
            r_join, z, heat, series_order=12
        )["angular_target"]
        target_order8 = heat_tail_moments(
            r_join, z, heat, series_order=8
        )["angular_target"]
        f_axis = float(core.F(0.0, z))
        cone_bound = f_axis * r_join * r_join
        cauchy_bound = 3.0 * target * target / (4.0 * r_join**3)
        monotone = _sample_monotone(joined, z)
        monotone_condition_applicable = bool(
            monotone["sampled_nonincreasing"] and f_axis > 0.0
        )
        row = {
            "Z": z,
            "F_axis": f_axis,
            "F_at_Rjoin": float(joined.F(r_join, z)),
            "Mtheta_inner": mtheta,
            "Mp_inner": mp,
            "Mtheta_heat_target": float(target),
            "Mtheta_target_order8_minus_order12": float(target_order8 - target),
            "Mtheta_inner_minus_heat_target": float(mtheta - target),
            "Mp_inner_minus_cauchy_minimum": float(mp - cauchy_bound),
            "monotone_cone_bound": cone_bound,
            "heat_target_minus_monotone_cone_bound": float(target - cone_bound),
            "monotone_cone_bound_satisfied_for_target": bool(target <= cone_bound),
            "monotone_cone_condition_applicable_on_sample": monotone_condition_applicable,
            "monotone_cone_failure_under_sampled_assumption": bool(
                monotone_condition_applicable and target > cone_bound
            ),
            "cauchy_minimum_Mp_for_heat_target": cauchy_bound,
            "cauchy_pressure_bound_satisfied_by_inner_Mp": bool(mp + 1.0e-12 >= cauchy_bound),
            "quadrature_n": int(n),
            "quadrature_n32_Mtheta": mtheta_32,
            "quadrature_n32_Mp": mp_32,
            "quadrature_n_minus_n32_Mtheta": float(mtheta - mtheta_32),
            "quadrature_n_minus_n32_Mp": float(mp - mp_32),
            "monotonicity_sample": monotone,
            "rough_minimum_Rjoin_for_fixed_Faxis": _rough_minimum_join(
                heat, z, f_axis, r_join
            ),
        }
        rows.append(row)

    report: dict[str, Any] = {
        "status": "completed",
        "source": SOURCE,
        "source_version": SOURCE_VERSION,
        "source_section": SOURCE_SECTION,
        "profile_source": "lei_ren_part1_pressure_core.load_core()",
        "profile_json": str(Path(__file__).with_name("lei_ren_part1_pressure_core.json")),
        "R_core": r_core,
        "R_join": r_join,
        "heat": {"h": float(heat.h), "c_inf": float(heat.c_inf)},
        "Z": z_values.tolist(),
        "quadrature": {
            "primary_order": int(n),
            "comparison_order": 32,
            "split_edges": [0.0, r_core, r_join],
        },
        "rows": rows,
        "necessary_condition_scope": {
            "monotone_cone": (
                "If F_R<=0 on [0,Rjoin] and F_axis>0, then "
                "Mtheta_inner<=F_axis*Rjoin^2."
            ),
            "cauchy_pressure": (
                "For a requested Mtheta target on [0,Rjoin], "
                "Mp>=3*Mtheta^2/(4*Rjoin^3)."
            ),
            "rough_radius_roots": (
                "The reported roots use fixed saved F_axis and heat angular "
                "target; the inequality is necessary, not sufficient."
            ),
        },
        "pde_validated": False,
        "global_impossibility_claim": False,
        "scale_recursion_established": False,
        "full_five_moment_matching": False,
        "scope": (
            "Finite compatibility diagnostic for the saved independent short "
            "core/blend at Rjoin=.2. It does not alter the profile or certify "
            "the paper exterior/collar, admissible stress, PDE, or global "
            "energy."
        ),
    }
    output = Path(__file__).with_suffix(".json")
    output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report))
    return report


if __name__ == "__main__":
    run()


__all__ = ["run", "SOURCE", "SOURCE_SECTION", "SOURCE_VERSION"]
