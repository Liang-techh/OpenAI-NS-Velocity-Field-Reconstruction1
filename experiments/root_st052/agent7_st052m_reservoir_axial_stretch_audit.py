"""Audit axial-stretch/core-radius capacity of the frozen outer-reservoir basis.

Preregistered in issue #784 and stacked on exact Agent-7 PR #775 head.
This increment changes no candidate velocity.  It analyzes the axial component of
the same one-channel vector-potential correction already used by #775.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

import agent7_st052m_outer_axial_return_reservoir_preflight as reservoir

TASK_ID = "CR003-ST052M-RESERVOIR-AXIAL-STRETCH-AUDIT-111"
PREREG_ISSUE = 784
SOURCE_PARENT_PR = 775
SOURCE_PARENT_HEAD = "93887a59729d22113badf2ae4dac2f7d868e2703"
SOURCE_RESERVOIR_PR = 767

Q_OUTER = 64.0 / 25.0
EXPECTED_Q_ROOT = (13.0 + 5.0 * math.sqrt(17.0)) / 30.0
EXPECTED_R_ROOT = math.sqrt(EXPECTED_Q_ROOT)
PROBE_ABS_Z = 1.55
PROBE_RADII = (0.6, 0.9, 1.2, 1.5)
RADIAL_GRIDS = (5001, 20001, 80001)
FLUX_TOL = 1.0e-9
REFINEMENT_TOL = 1.0e-9

TRUTH = {
    "candidate_velocity_changed": False,
    "canonical_velocity_changed": False,
    "saved_velocity_changed": False,
    "basis_dimension_changed": False,
    "candidate_coefficient_changed": False,
    "time_law_changed": False,
    "pressure_or_force_changed": False,
    "trajectory_replay_performed": False,
    "fresh_714_path_data_used": False,
    "held_out_pde_residual_evaluated": False,
    "public_image_numeric_target_used": False,
    "core_thickness_source_threshold_defined": False,
    "visualization_ready": False,
    "visual_correspondence_verified": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


def _assert_source_lock() -> None:
    reservoir._assert_source_lock()
    if reservoir.RADIAL_Q_A != -1.0 or reservoir.RADIAL_Q_B != Q_OUTER:
        raise RuntimeError("frozen reservoir radial bump drifted")
    if reservoir.R_SUPPORT_MAX != 1.6:
        raise RuntimeError("frozen reservoir radial support drifted")
    if reservoir.Z_VISIBLE_OUTER != 1.8 or reservoir.Z_GLOBAL_OUTER != 2.0:
        raise RuntimeError("frozen reservoir axial support drifted")


def axial_shape_factor(radius: np.ndarray | float) -> np.ndarray:
    """Return C_z/Z = 2(R+q R_q) for the live frozen radial bump."""
    rr = np.asarray(radius, dtype=float)
    q = rr * rr
    radial, radial_q = reservoir._radial_profile_and_derivative(q)
    return 2.0 * (radial + q * radial_q)


def analytic_root_geometry() -> dict[str, Any]:
    # 25 * [(q+1)(64/25-q) + 4q(39/25-2q)]
    # = -225 q^2 + 195 q + 64.
    roots = np.roots(np.asarray([-225.0, 195.0, 64.0], dtype=float))
    roots = np.sort(np.asarray(roots, dtype=float))
    physical = [float(q) for q in roots if 0.0 <= q <= Q_OUTER]
    if len(physical) != 1:
        raise RuntimeError("expected exactly one physical axial sign-change root")
    q_root = physical[0]
    r_root = math.sqrt(q_root)
    if abs(q_root - EXPECTED_Q_ROOT) > 5.0e-14:
        raise RuntimeError("analytic q-root identity drifted")
    if abs(r_root - EXPECTED_R_ROOT) > 5.0e-14:
        raise RuntimeError("analytic r-root identity drifted")
    return {
        "polynomial_coefficients": [-225, 195, 64],
        "all_q_roots": roots.tolist(),
        "physical_q_root": q_root,
        "exact_q_root_formula": "(13+5*sqrt(17))/30",
        "physical_r_root": r_root,
        "exact_r_root_formula": "sqrt((13+5*sqrt(17))/30)",
        "support_radius": reservoir.R_SUPPORT_MAX,
        "core_cross_section_area_fraction": float(q_root / Q_OUTER),
        "sign_rule_positive_z": "axial outward/stretch for r<r*, axial return for r>r*",
        "sign_rule_negative_z": "axial outward/stretch away from midplane for r<r*, opposite annular return",
    }


def live_probe_signs() -> dict[str, Any]:
    rows: dict[str, Any] = {}
    for radius in PROBE_RADII:
        plus = np.asarray(reservoir.reservoir_correction(np.asarray([[radius, 0.0, PROBE_ABS_Z]])), float)[0]
        minus = np.asarray(reservoir.reservoir_correction(np.asarray([[radius, 0.0, -PROBE_ABS_Z]])), float)[0]
        expected_core = radius < EXPECTED_R_ROOT
        gate = bool(
            (plus[2] > 0.0 and minus[2] < 0.0)
            if expected_core
            else (plus[2] < 0.0 and minus[2] > 0.0)
        )
        if not gate:
            raise RuntimeError(f"live axial sign gate failed at r={radius}")
        if not (plus[0] < 0.0 and minus[0] < 0.0):
            raise RuntimeError(f"visible radial inwardness drifted at r={radius}")
        rows[f"r={radius:.2f}"] = {
            "positive_z_correction": plus.tolist(),
            "negative_z_correction": minus.tolist(),
            "inside_axial_stretch_core": expected_core,
            "axial_sign_gate": gate,
            "radial_inward_both_z_signs": True,
        }
    return {
        "abs_z": PROBE_ABS_Z,
        "r_star": EXPECTED_R_ROOT,
        "by_radius": rows,
        "all_gates_pass": True,
    }


def radial_peak_geometry() -> dict[str, float]:
    radii = np.linspace(0.0, reservoir.R_SUPPORT_MAX, 200001)
    response = np.asarray(axial_shape_factor(radii), dtype=float)
    i_pos = int(np.argmax(response))
    i_neg = int(np.argmin(response))
    if not (radii[i_pos] < EXPECTED_R_ROOT < radii[i_neg]):
        raise RuntimeError("axial response peak ordering drifted")
    return {
        "positive_core_peak_radius": float(radii[i_pos]),
        "positive_core_peak_Cz_over_Z": float(response[i_pos]),
        "negative_annular_peak_radius": float(radii[i_neg]),
        "negative_annular_peak_Cz_over_Z": float(response[i_neg]),
    }


def _flux_level(n: int) -> dict[str, float]:
    radii = np.linspace(0.0, reservoir.R_SUPPORT_MAX, int(n))
    response = np.asarray(axial_shape_factor(radii), dtype=float)
    integrand = 2.0 * math.pi * radii * response
    flux = float(np.trapezoid(integrand, radii))
    positive_flux = float(np.trapezoid(2.0 * math.pi * radii * np.maximum(response, 0.0), radii))
    negative_flux = float(np.trapezoid(2.0 * math.pi * radii * np.minimum(response, 0.0), radii))
    return {
        "signed_flux_per_unit_Z": flux,
        "positive_core_flux_per_unit_Z": positive_flux,
        "negative_annular_flux_per_unit_Z": negative_flux,
        "positive_plus_negative": positive_flux + negative_flux,
    }


def axial_flux_refinement() -> dict[str, Any]:
    levels = {str(n): _flux_level(n) for n in RADIAL_GRIDS}
    medium = levels[str(RADIAL_GRIDS[-2])]
    fine = levels[str(RADIAL_GRIDS[-1])]
    max_fine_balance = max(
        abs(fine["signed_flux_per_unit_Z"]),
        abs(fine["positive_plus_negative"]),
    )
    max_refinement_change = max(
        abs(float(fine[key]) - float(medium[key]))
        for key in (
            "signed_flux_per_unit_Z",
            "positive_core_flux_per_unit_Z",
            "negative_annular_flux_per_unit_Z",
        )
    )
    passes = bool(max_fine_balance <= FLUX_TOL and max_refinement_change <= REFINEMENT_TOL)
    if not passes:
        raise RuntimeError("axial flux refinement gate failed")
    return {
        "analytic_identity": "integral C_z 2*pi*r dr = 2*pi*Z*[q*R(q)]_0^(64/25) = 0",
        "levels": levels,
        "max_abs_finest_balance": float(max_fine_balance),
        "max_medium_to_fine_change": float(max_refinement_change),
        "flux_tolerance": FLUX_TOL,
        "refinement_tolerance": REFINEMENT_TOL,
        "passes": passes,
    }


def run(out: Path) -> dict[str, Any]:
    _assert_source_lock()
    root = analytic_root_geometry()
    probes = live_probe_signs()
    peaks = radial_peak_geometry()
    flux = axial_flux_refinement()
    report: dict[str, Any] = {
        "task_id": TASK_ID,
        "prereg_issue": PREREG_ISSUE,
        "source_parent": {"pr": SOURCE_PARENT_PR, "head": SOURCE_PARENT_HEAD},
        "source_reservoir_pr": SOURCE_RESERVOIR_PR,
        "analytic_root_geometry": root,
        "live_visible_probe_signs": probes,
        "radial_peak_geometry": peaks,
        "axial_flux_refinement": flux,
        "capacity_conclusion": {
            "reservoir_channel_has_visible_inward_radial_and_core_axial_stretch_capacity": True,
            "axial_stretch_core_radius_is_fixed_by_current_radial_profile": True,
            "core_radius": root["physical_r_root"],
            "core_area_fraction": root["core_cross_section_area_fraction"],
            "same_dimension_radial_envelope_reshape_is_next_if_core_thickness_is_wrong": True,
            "second_poloidal_basis_justified_by_this_audit": False,
            "retune_775_from_this_audit_allowed": False,
        },
        **TRUTH,
    }
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    return report


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    report = run(args.out)
    print(json.dumps(report["capacity_conclusion"], indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
