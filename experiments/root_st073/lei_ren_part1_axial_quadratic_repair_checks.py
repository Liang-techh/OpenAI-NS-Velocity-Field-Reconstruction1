"""Focused numerical checks for the finite axial quadratic repair primitive."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "experiments" / "root_st073"
if str(LOCAL) not in sys.path:
    sys.path.insert(0, str(LOCAL))

from lei_ren_part1_axial_quadratic_repair import (  # noqa: E402
    SOURCE,
    SOURCE_SECTION,
    SOURCE_VERSION,
    evaluate_repair,
    solve_axial_quadratic_repair,
)


RMAX = 1.0
Z_VALUES = (-0.7, -0.3, 0.0, 0.35, 0.8)
TARGET_SCALE = 2.46463


def _F(R: float, Z: float) -> float:
    """Synthetic nonconstant swirl used only for this repair check."""

    return 0.7 + 0.35 * R + 0.12 * Z + 0.08 * np.sin(2.5 * R + 0.4 * Z)


def _base_U(R: float, Z: float) -> float:
    """Synthetic nonzero axial seed used only for this repair check."""

    return 0.3 + 0.4 * R - 0.2 * Z + 0.12 * np.cos(np.pi * R + 0.3 * Z)


def _quadratic_target(Z: float) -> float:
    """Nonzero test target near the reported collar scale, synthetic only."""

    return TARGET_SCALE * (1.0 + 0.04 * Z)


def _finite_ranges() -> dict[str, float]:
    radii = np.linspace(0.0, RMAX, 65)
    z_values = np.asarray(Z_VALUES, dtype=float)
    f_values = np.asarray([_F(float(r), float(z)) for r in radii for z in z_values])
    u_values = np.asarray(
        [_base_U(float(r), float(z)) for r in radii for z in z_values]
    )
    return {
        "F_min": float(np.min(f_values)),
        "F_max": float(np.max(f_values)),
        "F_range": float(np.ptp(f_values)),
        "base_U_min": float(np.min(u_values)),
        "base_U_max": float(np.max(u_values)),
        "base_U_range": float(np.ptp(u_values)),
    }


def run() -> dict[str, object]:
    rows = []
    all_n64_residuals = []
    all_n32_residuals = []
    coefficient_differences = []
    l2_differences = []
    for z in Z_VALUES:
        target = _quadratic_target(z)
        repair64 = solve_axial_quadratic_repair(
            _F,
            _base_U,
            RMAX,
            z,
            _quadratic_target,
            quadrature_order=64,
        )
        repair32 = solve_axial_quadratic_repair(
            _F,
            _base_U,
            RMAX,
            z,
            _quadratic_target,
            quadrature_order=32,
        )
        if not repair64.solvable or not repair32.solvable:
            raise AssertionError(f"synthetic target was not solvable at Z={z}")
        eval32 = evaluate_repair(repair64, quadrature_order=32)
        eval64 = evaluate_repair(repair64, quadrature_order=64)
        coefficient_difference = float(
            np.max(np.abs(repair64.coefficients - repair32.coefficients))
        )
        l2_difference = float(abs(eval64["l2_change"] - eval32["l2_change"]))
        all_n64_residuals.extend(
            [abs(eval64["mass_moment"]), abs(eval64["mixed_moment"]), abs(eval64["quadratic_residual"])]
        )
        all_n32_residuals.extend(
            [abs(eval32["mass_moment"]), abs(eval32["mixed_moment"]), abs(eval32["quadratic_residual"])]
        )
        coefficient_differences.append(coefficient_difference)
        l2_differences.append(l2_difference)
        rows.append(
            {
                "Z": float(z),
                "quadratic_target": float(target),
                "repair_order64": repair64.metadata(),
                "repair_order32": repair32.metadata(),
                "evaluation_order32": eval32,
                "evaluation_order64": eval64,
                "max_abs_coefficient_difference_order64_minus_order32": coefficient_difference,
                "abs_l2_change_difference_order64_minus_order32": l2_difference,
            }
        )

    unsolvable = solve_axial_quadratic_repair(
        _F, _base_U, RMAX, 0.0, -100.0, quadrature_order=64
    )
    validation = {}
    for label, thunk in {
        "nonpositive_Rmax": lambda: solve_axial_quadratic_repair(
            _F, _base_U, 0.0, 0.0, 1.0
        ),
        "noncallable_F": lambda: solve_axial_quadratic_repair(
            object(), _base_U, RMAX, 0.0, 1.0
        ),
        "too_few_bumps": lambda: solve_axial_quadratic_repair(
            _F, _base_U, RMAX, 0.0, 1.0,
            bump_intervals=((0.1, 0.2), (0.4, 0.5)),
        ),
    }.items():
        try:
            thunk()
        except (TypeError, ValueError):
            validation[label] = "passed"
        else:
            raise AssertionError(f"validation did not reject {label}")

    max_n64 = float(max(all_n64_residuals))
    max_n32 = float(max(all_n32_residuals))
    max_coeff_diff = float(max(coefficient_differences))
    max_l2_diff = float(max(l2_differences))
    assert max_n64 < 2.0e-9
    assert max_n32 < 2.0e-8
    assert max_coeff_diff < 2.0e-8
    assert max_l2_diff < 2.0e-8
    assert unsolvable.status == "quadratic_target_unreachable"
    assert all(row["repair_order64"]["linear_rank"] == 2 for row in rows)
    assert all(row["repair_order64"]["nullspace_dimension"] == 2 for row in rows)

    report = {
        "status": "completed",
        "source": SOURCE,
        "source_version": SOURCE_VERSION,
        "source_section": SOURCE_SECTION,
        "Rmax": RMAX,
        "Z_values": list(Z_VALUES),
        "bump_count": 4,
        "bump_supports": rows[0]["repair_order64"]["bump_intervals"],
        "synthetic_profiles": {
            "F": "0.7 + 0.35 R + 0.12 Z + 0.08 sin(2.5 R + 0.4 Z)",
            "base_U": "0.3 + 0.4 R - 0.2 Z + 0.12 cos(pi R + 0.3 Z)",
            "ranges": _finite_ranges(),
        },
        "synthetic_target": {
            "formula": "2.46463 * (1 + 0.04 Z)",
            "reference_scale": TARGET_SCALE,
            "actual_shared_background": False,
            "nonzero": True,
        },
        "rows": rows,
        "independent_quadrature_summary": {
            "max_abs_residual_order64": max_n64,
            "max_abs_residual_order32": max_n32,
            "max_abs_coefficient_difference_order64_minus_order32": max_coeff_diff,
            "max_abs_l2_change_difference_order64_minus_order32": max_l2_diff,
            "residual_components": [
                "integral U",
                "integral 2 R F U",
                "integral (U^2 - R F^2) - quadratic_target",
            ],
        },
        "unsolvable_target_diagnostic": unsolvable.metadata(),
        "finite_input_validation": validation,
        "minimal_L2_root_rule": (
            "The selected nullspace root is nearest the correction-L2 center "
            "on the quadratic level ellipse."
        ),
        "solves_generic_finite_subproblem_only": True,
        "full_outer_completed": False,
        "pde_validated": False,
        "admissibility_validated": False,
        "scope": (
            "Generic finite axial repair building block with four compact C-infinity "
            "bumps. The synthetic target is near the reported actual-collar "
            "quadratic scale but is not the actual shared background. No outer "
            "connection, five-moment closure, PDE, stress, or global certificate "
            "is claimed."
        ),
    }
    output = Path(__file__).with_suffix(".json")
    output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report))
    return report


if __name__ == "__main__":
    run()
