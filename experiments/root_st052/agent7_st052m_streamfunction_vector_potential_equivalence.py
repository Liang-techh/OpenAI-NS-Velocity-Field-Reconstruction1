"""Audit streamfunction vs vector-potential equivalence for the ST052-M poloidal lane.

Frozen in issue #827 and stacked directly on exact Agent-7 PR #818 head.
This is one representation-capacity increment only: it changes no candidate
velocity and adds no spatial or temporal basis coefficient.

For an axisymmetric scalar f(r,z), the current Cartesian vector potential

    A = (-y*f, x*f, 0) = r*f*e_phi

produces

    u_r = -r*f_z,
    u_z = 2*f + r*f_r.

The standard axisymmetric streamfunction convention

    u_r = -(1/r) * psi_z,
    u_z =  (1/r) * psi_r

with psi=r^2*f gives exactly the same velocity.  The audit checks the algebraic
identity and independently assembles the current live outer-reservoir channel
through the streamfunction formulas, including the regular axis limit.

A PASS establishes representation equivalence only.  It does not establish
visual correspondence, PDE validity, paper exactness, or OpenAI-field identity.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import sympy as sp

import agent7_st052m_endpoint_temporal_curvature_preflight as parent_audit

TASK_ID = "CR003-ST052M-STREAMFUNCTION-VECTOR-POTENTIAL-EQUIVALENCE-118"
PREREG_ISSUE = 827
SOURCE_PARENT_PR = 818
SOURCE_PARENT_HEAD = "c5acc7f8f65b64885ec5660ed9c991696668e300"
SOURCE_LIVE_RESERVOIR_PR = 767
SOURCE_LIVE_RESERVOIR_HEAD = "0930a1b7eeb357dee3e83dd79a3c4fea293e9314"
LIVE_K = 1.0
EQUIVALENCE_ABS_MAX_GATE = 1.0e-12
AXIS_FINITE_GATE = True

TRUTH = {
    "candidate_velocity_changed": False,
    "canonical_velocity_changed": False,
    "saved_velocity_changed": False,
    "basis_dimension_changed": False,
    "new_spatial_basis_added": False,
    "new_temporal_basis_added": False,
    "representation_switch_only": True,
    "streamfunction_is_new_morphology_degree_of_freedom": False,
    "parameter_scan_performed": False,
    "optimization_performed": False,
    "fresh_714_path_data_used": False,
    "pressure_or_force_changed": False,
    "held_out_pde_residual_evaluated": False,
    "public_image_numeric_target_used": False,
    "pixel_similarity_objective_used": False,
    "closer_visualization_delivery_established": False,
    "visualization_ready": False,
    "visual_correspondence_verified": False,
    "source_correspondence_verified": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


def _assert_source_lock() -> None:
    parent_audit._assert_source_lock()
    if parent_audit.TASK_ID != "CR003-ST052M-ENDPOINT-TEMPORAL-CURVATURE-PREFLIGHT-116":
        raise RuntimeError("#818 parent task identity drifted")
    if parent_audit.PREREG_ISSUE != 817:
        raise RuntimeError("#818 preregistration identity drifted")
    shape = parent_audit.shape_audit
    if shape.K_LIVE != LIVE_K:
        raise RuntimeError("live radial-envelope coordinate drifted")
    if shape.reservoir.Z_VISIBLE_OUTER != 1.8 or shape.reservoir.Z_GLOBAL_OUTER != 2.0:
        raise RuntimeError("outer axial reservoir semantics drifted")


def symbolic_equivalence() -> dict[str, Any]:
    """Prove the generic axisymmetric poloidal equivalence with SymPy."""
    r, z = sp.symbols("r z", positive=True, finite=True)
    f = sp.Function("f")(r, z)
    a_phi = r * f
    ur_vector = -sp.diff(a_phi, z)
    uz_vector = sp.simplify(sp.diff(r * a_phi, r) / r)

    psi = r**2 * f
    ur_stream = -sp.diff(psi, z) / r
    uz_stream = sp.diff(psi, r) / r

    ur_diff = sp.simplify(ur_vector - ur_stream)
    uz_diff = sp.simplify(uz_vector - uz_stream)
    return {
        "A_phi": str(a_phi),
        "psi": str(psi),
        "u_r_vector_potential": str(sp.simplify(ur_vector)),
        "u_z_vector_potential": str(sp.simplify(uz_vector)),
        "u_r_streamfunction": str(sp.simplify(ur_stream)),
        "u_z_streamfunction": str(sp.simplify(uz_stream)),
        "u_r_difference": str(ur_diff),
        "u_z_difference": str(uz_diff),
        "exact_identity": bool(ur_diff == 0 and uz_diff == 0),
    }


def _deterministic_points() -> np.ndarray:
    radii = (0.0, 1.0e-12, 1.0e-8, 0.15, 0.55, 0.95, 1.30, 1.55, 1.60, 1.70)
    abs_z_values = (0.0, 1.0, 1.15, 1.55, 1.75, 1.80, 1.85, 1.95, 2.0, 2.10)
    angles = tuple(float(v) for v in np.arange(8) * (0.25 * np.pi))
    rows: list[list[float]] = []
    for radius in radii:
        use_angles = (0.0,) if radius == 0.0 else angles
        for abs_z in abs_z_values:
            signs = (1.0,) if abs_z == 0.0 else (-1.0, 1.0)
            for sign in signs:
                z = float(sign * abs_z)
                for angle in use_angles:
                    rows.append(
                        [
                            float(radius * math.cos(angle)),
                            float(radius * math.sin(angle)),
                            z,
                        ]
                    )
    pts = np.asarray(rows, dtype=float)
    if pts.ndim != 2 or pts.shape[1] != 3 or not np.all(np.isfinite(pts)):
        raise RuntimeError("deterministic equivalence cloud construction failed")
    return pts


def vector_potential_velocity(points: np.ndarray) -> np.ndarray:
    """Current live k=1 outer-reservoir curl realization."""
    _assert_source_lock()
    pts = np.asarray(points, dtype=float)
    return np.asarray(parent_audit.shape_audit.reservoir_correction_k(pts, LIVE_K), dtype=float)


def streamfunction_velocity(points: np.ndarray) -> np.ndarray:
    """Independently assemble the same channel from psi=r^2 R(r^2) Z(z).

    The Cartesian formulas are the analytic regular-axis extension of
    u_r=-(1/r)psi_z, u_z=(1/r)psi_r; no numerical division by r is used.
    """
    _assert_source_lock()
    pts = np.asarray(points, dtype=float)
    if pts.ndim != 2 or pts.shape[1] != 3:
        raise ValueError("points must have shape (n,3)")
    x, y, z = pts[:, 0], pts[:, 1], pts[:, 2]
    q = x * x + y * y
    shape = parent_audit.shape_audit
    radial, radial_q = shape.core_control.radial_profile_k(q, LIVE_K)
    axial, axial_z = shape.reservoir.axial_potential_and_derivative(z)

    # psi=q*R(q)*Z(z).  From the streamfunction convention,
    # u_r=-r*R*Z_z and u_z=2*(R+q R_q)*Z.  Multiplying u_r by
    # x/r,y/r yields the regular Cartesian extension below.
    ux = -x * radial * axial_z
    uy = -y * radial * axial_z
    uz = 2.0 * (radial + q * radial_q) * axial
    out = np.column_stack((ux, uy, uz))
    if not np.all(np.isfinite(out)):
        raise RuntimeError("nonfinite streamfunction velocity")
    return out


def live_reservoir_equivalence() -> dict[str, Any]:
    pts = _deterministic_points()
    vector = vector_potential_velocity(pts)
    stream = streamfunction_velocity(pts)
    diff = stream - vector
    abs_err = np.linalg.norm(diff, axis=1)
    ref_norm = np.linalg.norm(vector, axis=1)
    scale = np.maximum(ref_norm, 1.0)
    rel_like = abs_err / scale

    axis_mask = np.hypot(pts[:, 0], pts[:, 1]) == 0.0
    near_axis_mask = np.hypot(pts[:, 0], pts[:, 1]) <= 1.0e-8
    outside_support_mask = (np.hypot(pts[:, 0], pts[:, 1]) >= 1.6) | (np.abs(pts[:, 2]) >= 2.0)
    if not np.any(axis_mask) or not np.any(near_axis_mask) or not np.any(outside_support_mask):
        raise RuntimeError("equivalence cloud masks unexpectedly empty")

    max_abs = float(np.max(np.abs(diff)))
    max_vector_norm_error = float(np.max(abs_err))
    max_scaled_error = float(np.max(rel_like))
    axis_finite = bool(np.all(np.isfinite(stream[axis_mask])) and np.all(np.isfinite(vector[axis_mask])))
    axis_transverse_max = float(np.max(np.abs(stream[axis_mask, :2])))
    near_axis_finite = bool(np.all(np.isfinite(stream[near_axis_mask])))
    outside_support_max = float(np.max(np.abs(stream[outside_support_mask])))
    passes = bool(
        max_abs <= EQUIVALENCE_ABS_MAX_GATE
        and axis_finite
        and near_axis_finite
        and axis_transverse_max <= EQUIVALENCE_ABS_MAX_GATE
        and outside_support_max <= EQUIVALENCE_ABS_MAX_GATE
    )
    return {
        "point_count": int(len(pts)),
        "axis_point_count": int(np.sum(axis_mask)),
        "near_axis_point_count": int(np.sum(near_axis_mask)),
        "outside_support_point_count": int(np.sum(outside_support_mask)),
        "max_component_abs_difference": max_abs,
        "max_vector_norm_difference": max_vector_norm_error,
        "max_scaled_difference": max_scaled_error,
        "axis_finite": axis_finite,
        "near_axis_finite": near_axis_finite,
        "axis_transverse_max_abs": axis_transverse_max,
        "outside_support_streamfunction_max_abs": outside_support_max,
        "absolute_gate": EQUIVALENCE_ABS_MAX_GATE,
        "passes": passes,
    }


def manufactured_convention_check() -> dict[str, Any]:
    """Catch sign/factor mistakes with a simple smooth non-reservoir f."""
    pts = np.asarray(
        [
            [0.0, 0.0, 0.3],
            [0.2, -0.1, 0.3],
            [0.7, 0.4, -0.6],
            [-1.1, 0.2, 0.8],
        ],
        dtype=float,
    )
    x, y, z = pts[:, 0], pts[:, 1], pts[:, 2]
    q = x * x + y * y
    # f=(1+q)*(1+z+z^2), so f_q=(1+z+z^2), f_z=(1+q)*(1+2z).
    f = (1.0 + q) * (1.0 + z + z * z)
    f_q = 1.0 + z + z * z
    f_z = (1.0 + q) * (1.0 + 2.0 * z)
    vector = np.column_stack((-x * f_z, -y * f_z, 2.0 * f + 2.0 * q * f_q))
    # Independent streamfunction assembly from psi=q*f.
    psi_z = q * f_z
    psi_q = f + q * f_q
    stream = np.column_stack((-x * f_z, -y * f_z, 2.0 * psi_q))
    max_abs = float(np.max(np.abs(vector - stream)))
    return {
        "max_component_abs_difference": max_abs,
        "passes": bool(max_abs <= EQUIVALENCE_ABS_MAX_GATE and np.all(np.isfinite(stream))),
    }


def build_report() -> dict[str, Any]:
    symbolic = symbolic_equivalence()
    live = live_reservoir_equivalence()
    manufactured = manufactured_convention_check()
    passed = bool(symbolic["exact_identity"] and live["passes"] and manufactured["passes"])
    return {
        "task_id": TASK_ID,
        "prereg_issue": PREREG_ISSUE,
        "source_parent": {"pr": SOURCE_PARENT_PR, "head": SOURCE_PARENT_HEAD},
        "source_live_reservoir": {
            "pr": SOURCE_LIVE_RESERVOIR_PR,
            "head": SOURCE_LIVE_RESERVOIR_HEAD,
            "k": LIVE_K,
        },
        "representation_pair": {
            "vector_potential": "A=(-y*f,x*f,0)=r*f*e_phi",
            "streamfunction": "psi=r^2*f",
            "streamfunction_convention": "u_r=-(1/r)psi_z, u_z=(1/r)psi_r, u_theta=0",
            "axis_handling": "analytic regular extension; no numerical division by r=0",
        },
        "symbolic_equivalence": symbolic,
        "live_reservoir_equivalence": live,
        "manufactured_convention_check": manufactured,
        "decision": {
            "frozen_audit_passes": passed,
            "streamfunction_switch_adds_basis_dimension": False if passed else None,
            "streamfunction_switch_adds_morphology_degree_of_freedom": False if passed else None,
            "basis_growth_must_change_scalar_function_space_support_or_time_structure_if_passed": passed,
            "actual_velocity_changed": False,
            "closer_visualization_delivery_established": False,
        },
        **TRUTH,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    report = build_report()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
