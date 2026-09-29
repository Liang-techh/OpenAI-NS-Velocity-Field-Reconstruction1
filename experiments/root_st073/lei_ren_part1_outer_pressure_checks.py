"""Focused checks for the finite core-to-heat pressure binding.

Checks use the requested h=0.001, c_inf=0.1 PaperCoreReference seed on
|Z|<=0.5.  They compare n=32/64 finite Gauss quadrature, verify terminal
pressure matching against the existing checked HeatExterior quadrature, and
differentiate the joined pressure independently away from the axis and joins.
No Uz, moment closure, stress cone, PDE, global-energy, or recursion claim is
made here.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
LOCAL = Path(__file__).resolve().parent
for path in (SRC, LOCAL):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from lei_ren_part1_outer_pressure import (
    SOURCE,
    SOURCE_SECTION,
    SOURCE_VERSION,
    build_reference_joined_profile,
    finite_pressure_moment,
)


def _values(function, z_values):
    return np.asarray([function(float(z)) for z in z_values], dtype=float)


def run() -> dict[str, object]:
    profile, reference = build_reference_joined_profile(h=0.001, c_inf=0.1)
    rjoin = profile.R_join
    z_values = np.linspace(-0.5, 0.5, 9)

    heat_pressure_join = _values(
        lambda z: profile.heat_pressure(rjoin, z), z_values
    )
    finite_increment = {
        str(order): _values(
            lambda z: finite_pressure_moment(profile.F, rjoin, z, n=order),
            z_values,
        )
        for order in (32, 64)
    }
    p0 = {
        str(order): _values(lambda z: profile.P0(z, n=order), z_values)
        for order in (32, 64)
    }
    p0_refinement = p0["64"] - p0["32"]
    increment_refinement = finite_increment["64"] - finite_increment["32"]
    p0_match = p0["64"] + finite_increment["64"]
    terminal_match_error = p0_match - heat_pressure_join

    old_autonomous_p0 = _values(lambda z: reference.Pi(0.0, z), z_values)
    autonomous_variation = p0["64"] - old_autonomous_p0

    # Central differences use pressure evaluation independent of the direct
    # F^2 quadrature identity and avoid R=0, R_core, and R_join.
    derivative_radii = np.asarray([0.025, 0.075, 0.15, 0.19, 0.25, 0.40])
    step = 1.0e-5
    derivative_rows = []
    derivative_errors = []
    for radius in derivative_radii:
        row_errors = []
        for z in z_values:
            derivative = (
                profile.pressure(radius + step, z, n=64)
                - profile.pressure(radius - step, z, n=64)
            ) / (2.0 * step)
            target = profile.F(radius, z) ** 2
            row_errors.append(float(derivative - target))
            derivative_errors.append(abs(float(derivative - target)))
        derivative_rows.append(
            {
                "R": float(radius),
                "max_abs_error": float(np.max(np.abs(row_errors))),
                "errors_by_Z": row_errors,
            }
        )

    # The axis uses only the core branch; the join uses the exact heat branch
    # value, so these checks also guard the intended branch protections.
    axis_values = _values(lambda z: profile.F(0.0, z), z_values)
    join_profile_error = _values(
        lambda z: profile.F(rjoin, z) - profile.F_heat(rjoin, z), z_values
    )
    exterior_profile_error = _values(
        lambda z: profile.F(0.4, z) - profile.F_heat(0.4, z), z_values
    )

    assert np.max(np.abs(terminal_match_error)) < 2.0e-12
    # The flat blend is smooth but the heat-factor calls are themselves
    # adaptive scalar quadratures; n=32 -> 64 therefore gives a finite,
    # measurable refinement rather than machine-zero agreement.
    assert np.max(np.abs(p0_refinement)) < 2.0e-8
    assert np.max(np.abs(join_profile_error)) < 2.0e-14
    assert np.max(np.abs(exterior_profile_error)) < 2.0e-14
    assert np.all(np.isfinite(axis_values))
    assert max(derivative_errors) < 2.0e-7

    reference_metadata = reference.metadata()
    report = {
        "status": "completed",
        "source": SOURCE,
        "source_version": SOURCE_VERSION,
        "source_section": SOURCE_SECTION,
        "profile": "PaperCoreReference F_core joined to HeatExterior F_heat",
        "parameters": {
            "h": 0.001,
            "c_inf": 0.1,
            "R_core": float(profile.R_core),
            "R_join": float(profile.R_join),
            "Z_interval": [-0.5, 0.5],
            "Z_count": int(len(z_values)),
        },
        "Z": z_values.tolist(),
        "quadrature_orders": [32, 64],
        "heat_pressure_join_checkedquad": heat_pressure_join.tolist(),
        "finite_pressure_increment_n32": finite_increment["32"].tolist(),
        "finite_pressure_increment_n64": finite_increment["64"].tolist(),
        "finite_pressure_increment_refinement_n64_minus_n32": (
            increment_refinement.tolist()
        ),
        "P0_n32": p0["32"].tolist(),
        "P0_n64": p0["64"].tolist(),
        "P0_refinement_n64_minus_n32": p0_refinement.tolist(),
        "pressure_match_at_Rjoin_n64_minus_heat": terminal_match_error.tolist(),
        "pressure_match_max_abs": float(np.max(np.abs(terminal_match_error))),
        "old_autonomous_P0": old_autonomous_p0.tolist(),
        "variation_from_old_autonomous_P0_n64": autonomous_variation.tolist(),
        "variation_from_old_autonomous_P0_max_abs": float(
            np.max(np.abs(autonomous_variation))
        ),
        "pressure_derivative_check": {
            "formula": "d_R P = F(R,Z)^2",
            "central_difference_step": step,
            "radii_away_from_axis_and_joins": derivative_radii.tolist(),
            "rows": derivative_rows,
            "max_abs_error": float(max(derivative_errors)),
        },
        "axis_F_values": axis_values.tolist(),
        "joined_profile_max_abs_error": float(np.max(np.abs(join_profile_error))),
        "exterior_profile_max_abs_error": float(
            np.max(np.abs(exterior_profile_error))
        ),
        "profile_metadata": reference_metadata,
        "uz_owned": False,
        "five_moment_closure": False,
        "admissible_stress": False,
        "finite_energy_certified": False,
        "global_profile_certified": False,
        "pde_validated": False,
        "scale_recursion_established": False,
        "scope": (
            "Finite local joined swirl-pressure seed on |Z|<=0.5: core F "
            "through R=0.05, flat blend to R=0.2, exact HeatExterior pressure "
            "outside. The pressure datum is assembled from this same F and an "
            "autonomous core reference; this is not the paper O(1)--O(8) "
            "exterior/collar, moment closure, admissibility, PDE, global-energy, "
            "or scale-recursion certification."
        ),
    }
    output = Path(__file__).with_suffix(".json")
    output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report))
    return report


if __name__ == "__main__":
    run()
