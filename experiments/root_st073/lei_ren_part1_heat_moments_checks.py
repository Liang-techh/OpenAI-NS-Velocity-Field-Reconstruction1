"""Focused checks for the Lei--Ren exact heat-tail moment evaluator."""

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

from lei_ren_part1_heat_moments import (
    SOURCE,
    SOURCE_SECTIONS,
    SOURCE_VERSION,
    direct_heat_tail_integrals,
    heat_tail_moments,
)
from openai_ns_reconstruction.heat_exterior import HeatExterior


def _moments_on_grid(R, z_values, heat, order):
    return [heat_tail_moments(R, float(z), heat, series_order=order)
            for z in z_values]


def _array(rows, key):
    return np.asarray([row[key] for row in rows], dtype=float)


def _json_rows(rows):
    output = []
    for row in rows:
        output.append({
            key: (int(value) if key == "series_order" else float(value))
            for key, value in row.items()
            if key in {
                "angular", "angular_target", "axial", "mixed", "quadratic",
                "quadratic_tail", "pressure", "pressure_tail",
                "angular_tail_difference", "angular_tail_difference_bracket",
                "quadratic_integral", "series_order",
                "R", "Z",
            }
        })
    return output


def run() -> dict[str, object]:
    radius = 0.2
    z_values = np.linspace(-0.5, 0.5, 9)
    primary_heat = HeatExterior(h=0.001, c_inf=0.1)
    order8 = _moments_on_grid(radius, z_values, primary_heat, 8)
    order12 = _moments_on_grid(radius, z_values, primary_heat, 12)
    refinement = {
        key: _array(order12, key) - _array(order8, key)
        for key in ("angular", "quadratic", "pressure", "angular_tail_difference",
                    "angular_tail_difference_bracket", "quadratic_integral")
    }
    refinement_max = {
        key: float(np.max(np.abs(value))) for key, value in refinement.items()
    }

    # The direct comparison deliberately uses a larger h where direct H-1
    # quadrature is numerically less delicate.  It is geometry-only evidence.
    comparison_z = np.asarray([-0.4, 0.0, 0.4])
    larger_heat = HeatExterior(h=0.05, c_inf=0.1)
    larger_primary = _moments_on_grid(radius, comparison_z, larger_heat, 8)
    larger_direct = [
        direct_heat_tail_integrals(radius, float(z), larger_heat)
        for z in comparison_z
    ]
    direct_errors = {
        key: np.asarray(
            [primary[key] - direct[key] for primary, direct in
             zip(larger_primary, larger_direct)],
            dtype=float,
        )
        for key in ("angular", "quadratic", "angular_tail_difference",
                    "angular_tail_difference_bracket", "quadratic_integral")
    }
    direct_error_max = {
        key: float(np.max(np.abs(value))) for key, value in direct_errors.items()
    }

    # Independent radial derivative identities at the requested R and Z.
    derivative_rows = []
    derivative_errors = []
    step = 1.0e-5
    for z in (-0.4, 0.0, 0.4):
        left = heat_tail_moments(radius - step, z, primary_heat, series_order=12)
        right = heat_tail_moments(radius + step, z, primary_heat, series_order=12)
        F = primary_heat.profile_E(radius, z) / np.sqrt(2.0 * radius)
        derivative_angular = (right["angular"] - left["angular"]) / (2.0 * step)
        derivative_quadratic = (
            right["quadratic_tail"] - left["quadratic_tail"]
        ) / (2.0 * step)
        derivative_pressure = (
            right["pressure_tail"] - left["pressure_tail"]
        ) / (2.0 * step)
        errors = {
            "angular": float(derivative_angular - 2.0 * radius * F),
            "quadratic": float(derivative_quadratic - radius * F * F),
            "pressure": float(derivative_pressure + F * F),
        }
        derivative_errors.extend(abs(value) for value in errors.values())
        derivative_rows.append({
            "Z": float(z),
            "dR_angular_minus_2RF": errors["angular"],
            "dR_quadratic_minus_RF2": errors["quadratic"],
            "dR_pressure_plus_F2": errors["pressure"],
        })

    validation_checks = {}
    for label, thunk in {
        "R_zero": lambda: heat_tail_moments(0.0, 0.0, primary_heat),
        "R_nonfinite": lambda: heat_tail_moments(np.inf, 0.0, primary_heat),
        "Z_outside": lambda: heat_tail_moments(radius, 1.1, primary_heat),
        "series_order_zero": lambda: heat_tail_moments(
            radius, 0.0, primary_heat, series_order=0
        ),
        "wrong_heat_type": lambda: heat_tail_moments(radius, 0.0, object()),
    }.items():
        try:
            thunk()
        except (TypeError, ValueError):
            validation_checks[label] = "passed"
        else:
            raise AssertionError(f"tail input validation did not reject {label}")

    assert max(refinement_max.values()) < 2.0e-10
    assert max(direct_error_max.values()) < 2.0e-9
    assert max(derivative_errors) < 2.0e-7
    assert all(row["axial"] == 0.0 and row["mixed"] == 0.0 for row in order12)

    report = {
        "status": "completed",
        "source": SOURCE,
        "source_version": SOURCE_VERSION,
        "source_sections": SOURCE_SECTIONS,
        "primary_heat": {"h": 0.001, "c_inf": 0.1},
        "R": radius,
        "Z": z_values.tolist(),
        "moments_series_order_8": _json_rows(order8),
        "moments_series_order_12": _json_rows(order12),
        "series_refinement_order12_minus_order8": {
            key: value.tolist() for key, value in refinement.items()
        },
        "series_refinement_max_abs": refinement_max,
        "larger_h_direct_comparison": {
            "h": 0.05,
            "c_inf": 0.1,
            "Z": comparison_z.tolist(),
            "primary_order": 8,
            "primary": _json_rows(larger_primary),
            "direct": larger_direct,
            "errors_primary_minus_direct": {
                key: value.tolist() for key, value in direct_errors.items()
            },
            "max_abs_error": direct_error_max,
            "scope": "Heat geometry diagnostic only; h=.05 is excluded from Part I certification.",
        },
        "radial_derivative_identities": {
            "R": radius,
            "central_difference_step": step,
            "rows": derivative_rows,
            "max_abs_error": float(max(derivative_errors)),
            "identities": {
                "angular": "d_R Mtheta_target = 2 R F_heat",
                "quadratic": "d_R Mztheta_tail = R F_heat^2",
                "pressure": "d_R Mp_tail = -F_heat^2",
            },
        },
        "finite_input_validation": validation_checks,
        "numerical_errors_are_not_rigorous_global_bounds": True,
        "uz_tail": 0.0,
        "mixed_tail": 0.0,
        "pde_validated": False,
        "scale_recursion_established": False,
        "global_certificate": False,
        "scope": (
            "Exact supplied heat-tail moment formulas at finite R=.2 on "
            "|Z|<=.5, with stable small-xi remainder quadrature and an h=.05 "
            "direct diagnostic. This is tail accounting only: no core match, "
            "full five-moment closure, admissible stress, PDE, global energy, "
            "or scale-recursion certificate is claimed."
        ),
    }
    output = Path(__file__).with_suffix(".json")
    output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report))
    return report


if __name__ == "__main__":
    run()
