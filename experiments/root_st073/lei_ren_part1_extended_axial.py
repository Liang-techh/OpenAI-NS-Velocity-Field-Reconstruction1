"""Finite axial quadratic repair on the saved ExtendedSwirl candidate.

This module binds the reusable finite axial repair primitive to the saved
shared-pressure candidate returned by
``lei_ren_part1_extended_pressure_core.load_extended_profile()``.  It keeps
the candidate's actual ``F(R,Z)`` and defines a finite axial seed from the
new core's ``U``: core ``U`` is retained through ``R=.005``, multiplied by a
flat C-infinity cutoff to zero by ``R=.02``, and set to zero afterwards.

For each requested axial slice, the quadratic target is the actual
``exterior_targets`` value at the collar inner radius ``R_a``.  Integration
is split at the core/cutoff/anchor/cross geometry and on logarithmic points
through the long radial interval.  The repair remains a finite three-moment
candidate for the shared-pressure ExtendedSwirl profile.  It does not claim a
final field, Z-derivative closure, the whole cone, a PDE, or complete
five-moment matching.

Source pinned to Z. Lei and X. Ren, arXiv:2609.35406v1, Section 2.5.
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

from lei_ren_part1_axial_quadratic_repair import (  # noqa: E402
    SOURCE,
    evaluate_repair,
    solve_axial_quadratic_repair,
)
from lei_ren_part1_exterior_targets import exterior_targets  # noqa: E402
from lei_ren_part1_extended_pressure_core import load_extended_profile  # noqa: E402


R_CORE = 0.005
R_ANCHOR = 0.01
R_CUTOFF = 0.02
Z_VALUES = (-1.0, -0.5, 0.0, 0.5, 1.0)
SOLVE_ORDER = 64
INDEPENDENT_ORDER = 128
LOG_BREAK_COUNT = 24
RESIDUAL_TOLERANCE = 1.0e-8


def _flat_step(value: float) -> float:
    """C-infinity step from zero to one on the unit interval."""

    if value <= 0.0:
        return 0.0
    if value >= 1.0:
        return 1.0
    left = math.exp(-1.0 / value)
    right = math.exp(-1.0 / (1.0 - value))
    return left / (left + right)


def _base_U(core: Any, R: float, Z: float) -> float:
    """Retain core ``U`` then taper it flatly to zero by ``R=.02``."""

    radius = float(R)
    z = float(Z)
    if radius <= R_CORE:
        return float(core.U(radius, z))
    if radius >= R_CUTOFF:
        return 0.0
    cutoff = 1.0 - _flat_step((radius - R_CORE) / (R_CUTOFF - R_CORE))
    return float(cutoff * core.U(radius, z))


def _quadrature_breakpoints(profile: Any, Z: float, Rmax: float) -> tuple[float, ...]:
    """Return geometry and logarithmic breaks required by the long profile."""

    cross = float(profile._cross(float(Z)))
    log_points = np.geomspace(R_CUTOFF, Rmax, LOG_BREAK_COUNT)
    candidates = [R_CORE, R_ANCHOR, R_CUTOFF, cross]
    candidates.extend(float(value) for value in log_points[1:-1])
    return tuple(sorted({value for value in candidates if 0.0 < value < Rmax}))


def _target(profile: Any, Rmax: float, Z: float) -> dict[str, Any]:
    """Get the actual collar-to-heat target at one axial slice."""

    return exterior_targets(
        profile.collar.F,
        profile.collar.heat,
        Rmax,
        profile.collar.R_b,
        float(Z),
        n=64,
    )


def run() -> dict[str, Any]:
    """Solve and independently reintegrate the five actual axial slices."""

    profile = load_extended_profile()
    core = profile.core
    Rmax = float(profile.collar.collar_inner_radius)
    rows: list[dict[str, Any]] = []
    for z in Z_VALUES:
        z = float(z)
        target = _target(profile, Rmax, z)
        breakpoints = _quadrature_breakpoints(profile, z, Rmax)
        repair = solve_axial_quadratic_repair(
            profile.F,
            lambda R, Z, c=core: _base_U(c, R, Z),
            Rmax,
            z,
            target["quadratic_target"],
            quadrature_breakpoints=breakpoints,
            quadrature_order=SOLVE_ORDER,
        )
        if not repair.solvable:
            raise AssertionError(
                f"actual ExtendedSwirl quadratic repair is unsolved at Z={z}: {repair.status}"
            )
        evaluated_64 = evaluate_repair(repair, quadrature_order=SOLVE_ORDER)
        evaluated_128 = evaluate_repair(repair, quadrature_order=INDEPENDENT_ORDER)
        rows.append(
            {
                "Z": z,
                "Rmax_Ra": Rmax,
                "cross_radius": float(profile._cross(z)),
                "breakpoint_count": len(breakpoints),
                "quadrature_breakpoints": list(breakpoints),
                "exterior_targets": {
                    key: value
                    for key, value in target.items()
                    if key
                    in {
                        "R_inner",
                        "R_outer",
                        "Z",
                        "angular_target",
                        "quadratic_target",
                        "finite_exterior_angular",
                        "finite_exterior_quadratic",
                        "pure_swirl_exterior_assumed",
                    }
                },
                "repair": repair.metadata(),
                "evaluation_order64": evaluated_64,
                "independent_evaluation_order128": evaluated_128,
                "base_U_support": [R_CORE, R_CUTOFF],
                "actual_profile_source": "ExtendedSwirl(load_extended_profile())",
            }
        )

    max_residual_64 = max(
        max(
            abs(row["evaluation_order64"][name])
            for name in ("mass_moment", "mixed_moment", "quadratic_residual")
        )
        for row in rows
    )
    max_residual_128 = max(
        max(
            abs(row["independent_evaluation_order128"][name])
            for name in ("mass_moment", "mixed_moment", "quadratic_residual")
        )
        for row in rows
    )
    checks_passed = bool(
        np.isfinite(max_residual_64)
        and np.isfinite(max_residual_128)
        and max_residual_64 <= RESIDUAL_TOLERANCE
        and max_residual_128 <= RESIDUAL_TOLERANCE
    )
    report = {
        "status": "completed" if checks_passed else "failed",
        "checks_passed": checks_passed,
        "checks": {
            "residual_tolerance": RESIDUAL_TOLERANCE,
            "max_abs_residual_order64": float(max_residual_64),
            "max_abs_residual_order128": float(max_residual_128),
        },
        "source": SOURCE,
        "source_version": "2609.35406v1",
        "profile_source": "lei_ren_part1_extended_pressure_core.load_extended_profile()",
        "profile_candidate_scope": "shared-pressure ExtendedSwirl intermediate candidate",
        "R_core": float(profile.R_core),
        "R_anchor": float(profile.R_anchor),
        "R_a_collar_inner": Rmax,
        "R_b_collar_outer": float(profile.collar.R_b),
        "Z_values": list(Z_VALUES),
        "solve_quadrature_order": SOLVE_ORDER,
        "independent_quadrature_order": INDEPENDENT_ORDER,
        "log_break_count": LOG_BREAK_COUNT,
        "base_U_definition": "core.U through R=.005; C-infinity flat cutoff to zero by R=.02; zero outside",
        "rows": rows,
        "max_abs_residual_order64": float(max_residual_64),
        "max_abs_residual_order128": float(max_residual_128),
        "actual_three_moment_repair_only": True,
        "final_field_claim": False,
        "z_derivatives_closed": False,
        "whole_cone_validated": False,
        "pde_validated": False,
        "scope": (
            "Finite axial repair of the saved ExtendedSwirl intermediate using "
            "actual collar quadratic targets. This does not certify the final "
            "field, Z derivatives, whole stress cone, PDE, or full five-moment "
            "outer construction."
        ),
    }
    output = Path(__file__).with_suffix(".json")
    output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report))
    if not checks_passed:
        raise AssertionError(
            "actual ExtendedSwirl axial repair residual exceeded "
            f"tolerance {RESIDUAL_TOLERANCE:g}: "
            f"order64={max_residual_64:.6g}, order128={max_residual_128:.6g}"
        )
    return report


if __name__ == "__main__":
    run()
