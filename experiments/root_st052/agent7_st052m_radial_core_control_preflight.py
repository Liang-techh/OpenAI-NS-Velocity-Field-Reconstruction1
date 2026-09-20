"""Audit same-dimension radial control of the ST052-M reservoir stretch core.

Frozen in issue #800 and stacked on exact Agent-7 PR #794 head.  This is an
expression-capacity diagnostic only.  It changes no candidate velocity.

For the existing outer-reservoir vector-potential channel

    A=(-y*f,x*f,0),  f=R_k(r^2) Z(z),

keep the axial reservoir Z, support q<64/25, and quartic endpoint flatness, but
vary only the nonphysical inner root of the radial bump:

    R_k(q) propto ((q+k)(b-q))^4,  b=64/25.

The live #767 shape is k=1.  The axial response per unit Z is

    C_z/Z = 2 (R_k + q R_k'),

so its sign-change radius is a direct same-dimension core-thickness coordinate.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
from scipy.integrate import simpson
from scipy.optimize import brentq

import agent7_st052m_outer_axial_return_reservoir_preflight as reservoir

TASK_ID = "CR003-ST052M-RADIAL-CORE-CONTROL-PREFLIGHT-114"
PREREG_ISSUE = 800
SOURCE_PARENT_PR = 794
SOURCE_PARENT_HEAD = "5a3992e15f56aafa81e934c9ec306a2a10fa4e60"
SOURCE_RESERVOIR_PR = 767
SOURCE_UNIT_STRETCH_AUDIT_PR = 785

Q_OUTER = 64.0 / 25.0
R_SUPPORT = 1.6
EXPONENT = 4
K_VALUES = (0.50, 0.75, 1.00, 1.50, 2.00)
LIVE_K = 1.0
ROOT_TOL = 1.0e-10
SENSITIVITY_REL_TOL = 1.0e-6
FLUX_TOL = 1.0e-9
REFINEMENT_TOL = 1.0e-9
CONDITION_MAX = 25.0
RADIAL_GRIDS = (5001, 20001, 80001)
RESPONSE_RADIAL_N = 1601
PEAK_RADIAL_N = 50001
SENSITIVITY_H = 1.0e-3

TRUTH = {
    "candidate_velocity_changed": False,
    "canonical_velocity_changed": False,
    "saved_velocity_changed": False,
    "basis_dimension_changed": False,
    "candidate_coefficient_changed": False,
    "time_law_changed": False,
    "pressure_or_force_changed": False,
    "parameter_fit_performed": False,
    "candidate_selection_performed": False,
    "fresh_714_path_data_used": False,
    "held_out_pde_residual_evaluated": False,
    "public_image_numeric_target_used": False,
    "public_numeric_core_radius_target_defined": False,
    "visualization_ready": False,
    "visual_correspondence_verified": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


def _assert_source_lock() -> None:
    reservoir._assert_source_lock()
    if reservoir.RADIAL_Q_A != -1.0:
        raise RuntimeError("live reservoir pseudo-inner root drifted")
    if reservoir.RADIAL_Q_B != Q_OUTER:
        raise RuntimeError("live reservoir radial support q drifted")
    if reservoir.R_SUPPORT_MAX != R_SUPPORT:
        raise RuntimeError("live reservoir radial support radius drifted")
    if reservoir.Z_VISIBLE_OUTER != 1.8 or reservoir.Z_GLOBAL_OUTER != 2.0:
        raise RuntimeError("live reservoir axial semantics drifted")


def radial_profile_k(q: np.ndarray | float, k: float) -> tuple[np.ndarray, np.ndarray]:
    k = float(k)
    if k <= 0.0:
        raise ValueError("k must be positive")
    return reservoir.vp._poly_bump(np.asarray(q, dtype=float), -k, Q_OUTER)


def axial_shape_q(q: np.ndarray | float, k: float) -> np.ndarray:
    qq = np.asarray(q, dtype=float)
    radial, radial_q = radial_profile_k(qq, k)
    return 2.0 * (radial + qq * radial_q)


def axial_shape_radius(radius: np.ndarray | float, k: float) -> np.ndarray:
    rr = np.asarray(radius, dtype=float)
    return axial_shape_q(rr * rr, k)


def analytic_q_root(k: float) -> float:
    k = float(k)
    disc = 25.0 * (Q_OUTER - k) ** 2 + 36.0 * k * Q_OUTER
    return float((5.0 * (Q_OUTER - k) + math.sqrt(disc)) / 18.0)


def numerical_q_root(k: float) -> float:
    grid = np.linspace(0.0, Q_OUTER * 0.999, 20001)
    response = np.asarray(axial_shape_q(grid, k), dtype=float)
    idx = np.flatnonzero((response[:-1] > 0.0) & (response[1:] <= 0.0))
    if len(idx) != 1:
        raise RuntimeError(f"expected one physical live-response crossing for k={k}, got {len(idx)}")
    i = int(idx[0])
    return float(brentq(lambda q: float(axial_shape_q(q, k)), float(grid[i]), float(grid[i + 1])))


def root_geometry(k: float) -> dict[str, float]:
    q_analytic = analytic_q_root(k)
    q_numeric = numerical_q_root(k)
    if not (0.0 < q_analytic < Q_OUTER):
        raise RuntimeError(f"analytic root outside physical interval for k={k}")
    mismatch = abs(q_analytic - q_numeric)
    if mismatch > ROOT_TOL:
        raise RuntimeError(f"analytic/live numerical root mismatch for k={k}: {mismatch}")
    return {
        "k": float(k),
        "analytic_q_root": q_analytic,
        "numerical_q_root": q_numeric,
        "analytic_numeric_q_mismatch": mismatch,
        "r_star": math.sqrt(q_analytic),
        "core_area_fraction": q_analytic / Q_OUTER,
    }


def peak_geometry(k: float) -> dict[str, float]:
    radii = np.linspace(0.0, R_SUPPORT, PEAK_RADIAL_N)
    response = np.asarray(axial_shape_radius(radii, k), dtype=float)
    i_pos = int(np.argmax(response))
    i_neg = int(np.argmin(response))
    r_star = math.sqrt(analytic_q_root(k))
    if not (radii[i_pos] < r_star < radii[i_neg]):
        raise RuntimeError(f"peak ordering failed for k={k}")
    return {
        "positive_core_peak_radius": float(radii[i_pos]),
        "positive_core_peak_Cz_over_Z": float(response[i_pos]),
        "negative_annular_peak_radius": float(radii[i_neg]),
        "negative_annular_peak_Cz_over_Z": float(response[i_neg]),
    }


def response_geometry(k: float) -> dict[str, Any]:
    radii = np.linspace(0.0, R_SUPPORT * 0.999, RESPONSE_RADIAL_N)
    live = np.asarray(axial_shape_radius(radii, LIVE_K), dtype=float)
    trial = np.asarray(axial_shape_radius(radii, k), dtype=float)
    n_live = float(np.linalg.norm(live))
    n_trial = float(np.linalg.norm(trial))
    if min(n_live, n_trial) <= np.finfo(float).tiny:
        raise RuntimeError("degenerate response column")
    cosine = float(np.dot(live, trial) / (n_live * n_trial))
    if float(k) == LIVE_K:
        return {
            "response_cosine_vs_live_k1": 1.0,
            "two_column_rank_vs_live_k1": 1,
            "two_column_condition_vs_live_k1": None,
            "identical_to_live_shape": True,
        }
    matrix = np.column_stack((live / n_live, trial / n_trial))
    singular = np.linalg.svd(matrix, compute_uv=False)
    rank = int(np.linalg.matrix_rank(matrix, tol=1.0e-10))
    condition = float(singular[0] / singular[-1])
    if rank != 2 or condition > CONDITION_MAX:
        raise RuntimeError(f"shape response ill-conditioned for k={k}: rank={rank}, cond={condition}")
    return {
        "response_cosine_vs_live_k1": cosine,
        "two_column_rank_vs_live_k1": rank,
        "two_column_condition_vs_live_k1": condition,
        "two_column_singular_values": singular.tolist(),
        "identical_to_live_shape": False,
    }


def _flux_level(k: float, n: int) -> float:
    radii = np.linspace(0.0, R_SUPPORT, int(n))
    response = np.asarray(axial_shape_radius(radii, k), dtype=float)
    return float(simpson(2.0 * math.pi * radii * response, x=radii))


def flux_refinement(k: float) -> dict[str, Any]:
    levels = {str(n): _flux_level(k, n) for n in RADIAL_GRIDS}
    fine = float(levels[str(RADIAL_GRIDS[-1])])
    medium = float(levels[str(RADIAL_GRIDS[-2])])
    change = abs(fine - medium)
    if abs(fine) > FLUX_TOL or change > REFINEMENT_TOL:
        raise RuntimeError(f"signed axial-flux identity failed for k={k}")
    return {
        "analytic_identity": "integral C_z 2*pi*r dr = 2*pi*Z*[q*R_k(q)]_0^(64/25) = 0",
        "levels": levels,
        "finest_abs_signed_flux": abs(fine),
        "medium_to_fine_change": change,
        "passes": True,
    }


def support_check(k: float) -> dict[str, float]:
    radii = np.asarray([R_SUPPORT, 1.61, 1.75], dtype=float)
    response = np.asarray(axial_shape_radius(radii, k), dtype=float)
    max_abs = float(np.max(np.abs(response)))
    if max_abs != 0.0:
        raise RuntimeError(f"radial support leakage for k={k}")
    return {"outside_or_boundary_max_abs_Cz_over_Z": max_abs}


def sensitivity_at_live() -> dict[str, float]:
    k = LIVE_K
    q = analytic_q_root(k)
    f_q = -18.0 * q + 5.0 * (Q_OUTER - k)
    f_k = -5.0 * q + Q_OUTER
    dq_dk = -f_k / f_q
    dr_dk = dq_dk / (2.0 * math.sqrt(q))
    h = SENSITIVITY_H
    fd = (math.sqrt(analytic_q_root(k + h)) - math.sqrt(analytic_q_root(k - h))) / (2.0 * h)
    rel = abs(fd - dr_dk) / max(abs(dr_dk), np.finfo(float).tiny)
    if rel > SENSITIVITY_REL_TOL:
        raise RuntimeError("analytic-vs-FD core-radius sensitivity mismatch")
    return {
        "k": k,
        "analytic_dq_star_dk": dq_dk,
        "analytic_dr_star_dk": dr_dk,
        "fd_step": h,
        "fd_dr_star_dk": fd,
        "relative_mismatch": rel,
        "monotone_direction": "increasing k narrows the axial-stretch core" if dr_dk < 0.0 else "unexpected",
    }


def run(out: Path) -> dict[str, Any]:
    _assert_source_lock()
    rows: dict[str, Any] = {}
    r_stars = []
    for k in K_VALUES:
        root = root_geometry(k)
        r_stars.append(root["r_star"])
        rows[f"k={k:.2f}"] = {
            "root": root,
            "peaks": peak_geometry(k),
            "response_vs_live_k1": response_geometry(k),
            "support": support_check(k),
            "axial_flux": flux_refinement(k),
        }
    if not all(a > b for a, b in zip(r_stars, r_stars[1:])):
        raise RuntimeError("frozen k grid does not give strictly decreasing r_star")

    sensitivity = sensitivity_at_live()
    report: dict[str, Any] = {
        "task_id": TASK_ID,
        "prereg_issue": PREREG_ISSUE,
        "source_parent": {"pr": SOURCE_PARENT_PR, "head": SOURCE_PARENT_HEAD},
        "source_reservoir_pr": SOURCE_RESERVOIR_PR,
        "source_unit_stretch_audit_pr": SOURCE_UNIT_STRETCH_AUDIT_PR,
        "family": {
            "formula": "R_k(q) proportional to ((q+k)(64/25-q))^4 on 0<=q<64/25",
            "k_values": list(K_VALUES),
            "live_k": LIVE_K,
            "physical_support_radius": R_SUPPORT,
            "exponent": EXPONENT,
            "candidate_selection": None,
            "public_openai_core_radius_target": None,
        },
        "by_k": rows,
        "r_star_sequence": r_stars,
        "r_star_range": {"minimum": min(r_stars), "maximum": max(r_stars), "span": max(r_stars) - min(r_stars)},
        "live_sensitivity": sensitivity,
        "capacity_conclusion": {
            "same_dimension_core_thickness_control_demonstrated": True,
            "strict_monotone_core_radius_control_on_frozen_grid": True,
            "second_poloidal_basis_justified_by_this_audit": False,
            "candidate_k_selected": False,
            "total_child_visual_improvement_established": False,
            "next_if_fixed_render_shows_core_width_error": "freeze one radial-envelope reshape direction before adding a second poloidal basis",
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
    print(json.dumps({"r_star_range": report["r_star_range"], "live_sensitivity": report["live_sensitivity"], "capacity_conclusion": report["capacity_conclusion"]}, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
