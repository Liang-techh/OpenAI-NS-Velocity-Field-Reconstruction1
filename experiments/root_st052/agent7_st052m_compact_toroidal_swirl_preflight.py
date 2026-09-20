"""Preflight one compact toroidal/swirl basis against the ST052-M poloidal lane.

Frozen in issue #836 and stacked directly on exact Agent-7 PR #828 head.
This is one expression-capacity increment only: it changes no candidate velocity
and selects no new coefficient.

The existing outer-reservoir channel is poloidal.  This audit asks whether one
minimal axisymmetric toroidal direction can change swirl/circulation strength
without instantaneous radial/axial leakage:

    C_theta = (-y H(r^2,z), x H(r^2,z), 0),
    H(r^2,z) = R_1(r^2) W(z).

W is the even magnitude of the already-frozen C3 outer-reservoir axial
potential.  Therefore the new channel is compact on the same r<1.6, |z|<2
support, is axis regular, and has cylindrical components

    u_r = 0,  u_z = 0,  u_theta = r H.

For any axisymmetric H(r^2,z), div C_theta is identically zero.  A PASS only
establishes a well-conditioned morphology capacity direction.  It does not show
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

import agent7_st052m_streamfunction_vector_potential_equivalence as parent_equiv

TASK_ID = "CR003-ST052M-COMPACT-TOROIDAL-SWIRL-PREFLIGHT-119"
PREREG_ISSUE = 836
SOURCE_PARENT_PR = 828
SOURCE_PARENT_HEAD = "034a5e04df8900d8def05099fdbd1a287b43b4a0"
SOURCE_POLAR_RESERVOIR_PR = 767
SOURCE_POLAR_RESERVOIR_HEAD = "0930a1b7eeb357dee3e83dd79a3c4fea293e9314"
RELATED_FRESH_HOLDOUT_PR = 714

LIVE_K = 1.0
TIME_START = 0.25
RESPONSE_RADII = (0.20, 0.55, 0.90, 1.20, 1.50)
RESPONSE_ABS_Z = (1.15, 1.55, 1.75, 1.85, 1.95)
RESPONSE_ANGLES = tuple(float(v) for v in np.arange(8) * (0.25 * np.pi))
RESPONSE_TIMES = (0.375, 0.50, 0.625, 0.75)

ZERO_GATE = 1.0e-12
NONTRIVIAL_RMS_MIN = 1.0e-8
COSINE_MAX_ABS = 1.0e-12
CONDITION_MAX = 1.01

TRUTH = {
    "candidate_velocity_changed": False,
    "canonical_velocity_changed": False,
    "saved_velocity_changed": False,
    "basis_dimension_changed": False,
    "hypothetical_spatial_basis_capacity_tested": True,
    "toroidal_coefficient_selected": False,
    "new_temporal_basis_added": False,
    "parameter_scan_performed": False,
    "optimization_performed": False,
    "post_result_retuning_performed": False,
    "fresh_714_path_data_used": False,
    "pressure_or_force_changed": False,
    "held_out_pde_residual_evaluated": False,
    "public_image_numeric_target_used": False,
    "pixel_similarity_objective_used": False,
    "source_numeric_swirl_target_used": False,
    "closer_visualization_delivery_established": False,
    "visualization_ready": False,
    "visual_correspondence_verified": False,
    "source_correspondence_verified": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


def _shape_audit() -> Any:
    return parent_equiv.parent_audit.shape_audit


def _assert_source_lock() -> None:
    parent_equiv._assert_source_lock()
    if parent_equiv.TASK_ID != "CR003-ST052M-STREAMFUNCTION-VECTOR-POTENTIAL-EQUIVALENCE-118":
        raise RuntimeError("#828 parent task identity drifted")
    if parent_equiv.PREREG_ISSUE != 827:
        raise RuntimeError("#828 parent preregistration drifted")
    if parent_equiv.SOURCE_LIVE_RESERVOIR_PR != SOURCE_POLAR_RESERVOIR_PR:
        raise RuntimeError("live reservoir PR identity drifted")
    if parent_equiv.SOURCE_LIVE_RESERVOIR_HEAD != SOURCE_POLAR_RESERVOIR_HEAD:
        raise RuntimeError("live reservoir head drifted")
    shape = _shape_audit()
    if shape.K_LIVE != LIVE_K:
        raise RuntimeError("live radial-envelope coordinate drifted")
    if tuple(shape.RESPONSE_TIMES) != RESPONSE_TIMES:
        raise RuntimeError("registered response-time semantics drifted")
    if shape.reservoir.R_SUPPORT_MAX != 1.6:
        raise RuntimeError("radial support drifted")
    if shape.reservoir.Z_INNER != 1.0 or shape.reservoir.Z_VISIBLE_OUTER != 1.8:
        raise RuntimeError("visible axial support semantics drifted")
    if shape.reservoir.Z_GLOBAL_OUTER != 2.0:
        raise RuntimeError("global axial support semantics drifted")
    if RELATED_FRESH_HOLDOUT_PR != 714:
        raise RuntimeError("fresh holdout identity drifted")


def activation_g1(time: float) -> float:
    """Existing linear activation used only to make the spacetime response cloud."""
    return float(2.0 * (float(time) - TIME_START))


def even_axial_window(z: np.ndarray | float) -> np.ndarray:
    """Frozen even C3 window W=|Z_reservoir| on 1<|z|<2."""
    zz = np.asarray(z, dtype=float)
    axial, _ = _shape_audit().reservoir.axial_potential_and_derivative(zz)
    return np.abs(np.asarray(axial, dtype=float))


def toroidal_unit_velocity(points: np.ndarray) -> np.ndarray:
    """Compact axisymmetric toroidal unit channel in Cartesian coordinates."""
    _assert_source_lock()
    pts = np.asarray(points, dtype=float)
    if pts.ndim != 2 or pts.shape[1] != 3:
        raise ValueError("points must have shape (n,3)")
    if not np.all(np.isfinite(pts)):
        raise ValueError("points must be finite")
    x, y, z = pts[:, 0], pts[:, 1], pts[:, 2]
    q = x * x + y * y
    radial, _ = _shape_audit().core_control.radial_profile_k(q, LIVE_K)
    h = np.asarray(radial, dtype=float) * even_axial_window(z)
    out = np.column_stack((-y * h, x * h, np.zeros_like(h)))
    if not np.all(np.isfinite(out)):
        raise RuntimeError("nonfinite toroidal unit velocity")
    return out


def poloidal_unit_velocity(points: np.ndarray) -> np.ndarray:
    """Exact live k=1 outer-reservoir poloidal unit channel from #828/#767."""
    _assert_source_lock()
    return np.asarray(parent_equiv.vector_potential_velocity(np.asarray(points, dtype=float)), dtype=float)


def cylindrical_components(points: np.ndarray, vectors: np.ndarray) -> np.ndarray:
    pts = np.asarray(points, dtype=float)
    vec = np.asarray(vectors, dtype=float)
    if pts.ndim != 2 or pts.shape[1] != 3 or vec.shape != pts.shape:
        raise ValueError("points/vectors must both have shape (n,3)")
    x, y = pts[:, 0], pts[:, 1]
    r = np.hypot(x, y)
    ur = np.zeros_like(r)
    utheta = np.zeros_like(r)
    nz = r > 0.0
    ur[nz] = (x[nz] * vec[nz, 0] + y[nz] * vec[nz, 1]) / r[nz]
    utheta[nz] = (-y[nz] * vec[nz, 0] + x[nz] * vec[nz, 1]) / r[nz]
    return np.column_stack((ur, utheta, vec[:, 2]))


def symbolic_toroidal_divergence() -> dict[str, Any]:
    """Exact divergence identity for (-y H(q,z), x H(q,z), 0)."""
    x, y, z = sp.symbols("x y z", real=True, finite=True)
    q = x**2 + y**2
    h = sp.Function("H")(q, z)
    vx = -y * h
    vy = x * h
    divergence = sp.simplify(sp.diff(vx, x) + sp.diff(vy, y))
    return {
        "velocity": "(-y*H(x^2+y^2,z), x*H(x^2+y^2,z), 0)",
        "divergence": str(divergence),
        "exact_identity": bool(divergence == 0),
    }


def _response_points() -> np.ndarray:
    rows: list[list[float]] = []
    for radius in RESPONSE_RADII:
        for abs_z in RESPONSE_ABS_Z:
            for sign in (-1.0, 1.0):
                for angle in RESPONSE_ANGLES:
                    rows.append(
                        [
                            float(radius * math.cos(angle)),
                            float(radius * math.sin(angle)),
                            float(sign * abs_z),
                        ]
                    )
    points = np.asarray(rows, dtype=float)
    if points.shape != (len(RESPONSE_RADII) * len(RESPONSE_ABS_Z) * 2 * len(RESPONSE_ANGLES), 3):
        raise RuntimeError("response cloud construction failed")
    return points


def _normalize_column(values: np.ndarray) -> tuple[np.ndarray, float]:
    arr = np.asarray(values, dtype=float).ravel()
    norm = float(np.linalg.norm(arr))
    if not np.isfinite(norm) or norm <= 0.0:
        raise RuntimeError("zero/nonfinite response column")
    return arr / norm, norm


def response_identifiability() -> dict[str, Any]:
    """Check that toroidal and live poloidal spacetime response directions differ."""
    _assert_source_lock()
    points = _response_points()
    polar_spatial = poloidal_unit_velocity(points)
    toroidal_spatial = toroidal_unit_velocity(points)

    polar_blocks: list[np.ndarray] = []
    toroidal_blocks: list[np.ndarray] = []
    for time in RESPONSE_TIMES:
        gain = activation_g1(time)
        polar_blocks.append((gain * polar_spatial).ravel())
        toroidal_blocks.append((gain * toroidal_spatial).ravel())

    polar_col, polar_norm = _normalize_column(np.concatenate(polar_blocks))
    toroidal_col, toroidal_norm = _normalize_column(np.concatenate(toroidal_blocks))
    matrix = np.column_stack((polar_col, toroidal_col))
    singular = np.linalg.svd(matrix, compute_uv=False)
    rank = int(np.linalg.matrix_rank(matrix, tol=1.0e-12))
    condition = float(singular[0] / singular[-1])
    cosine = float(np.dot(polar_col, toroidal_col))
    passes = bool(rank == 2 and abs(cosine) <= COSINE_MAX_ABS and condition <= CONDITION_MAX)
    return {
        "point_count_per_time": int(len(points)),
        "time_count": int(len(RESPONSE_TIMES)),
        "flattened_component_count": int(matrix.shape[0]),
        "poloidal_response_norm": polar_norm,
        "toroidal_response_norm": toroidal_norm,
        "normalized_column_rank": rank,
        "normalized_column_cosine": cosine,
        "normalized_column_condition": condition,
        "singular_values": [float(v) for v in singular],
        "cosine_abs_gate": COSINE_MAX_ABS,
        "condition_gate": CONDITION_MAX,
        "passes": passes,
    }


def channel_selectivity() -> dict[str, Any]:
    """Measure cylindrical leakage and nontriviality on the frozen response cloud."""
    points = _response_points()
    polar = cylindrical_components(points, poloidal_unit_velocity(points))
    toroidal = cylindrical_components(points, toroidal_unit_velocity(points))

    toroidal_ur_max = float(np.max(np.abs(toroidal[:, 0])))
    toroidal_uz_max = float(np.max(np.abs(toroidal[:, 2])))
    toroidal_utheta_rms = float(np.sqrt(np.mean(toroidal[:, 1] ** 2)))
    polar_utheta_max = float(np.max(np.abs(polar[:, 1])))
    polar_poloidal_rms = float(np.sqrt(np.mean(polar[:, 0] ** 2 + polar[:, 2] ** 2)))
    passes = bool(
        toroidal_ur_max <= ZERO_GATE
        and toroidal_uz_max <= ZERO_GATE
        and toroidal_utheta_rms > NONTRIVIAL_RMS_MIN
        and polar_utheta_max <= ZERO_GATE
        and polar_poloidal_rms > NONTRIVIAL_RMS_MIN
    )
    return {
        "toroidal_u_r_max_abs": toroidal_ur_max,
        "toroidal_u_z_max_abs": toroidal_uz_max,
        "toroidal_u_theta_rms": toroidal_utheta_rms,
        "poloidal_u_theta_max_abs": polar_utheta_max,
        "poloidal_rz_rms": polar_poloidal_rms,
        "zero_gate": ZERO_GATE,
        "nontrivial_rms_min": NONTRIVIAL_RMS_MIN,
        "passes": passes,
    }


def support_and_axis_check() -> dict[str, Any]:
    """Check exact support boundary/outside and regular-axis zeros for both channels."""
    rows: list[list[float]] = []
    for radius in (0.0, 1.60, 1.70):
        angles = (0.0,) if radius == 0.0 else (0.0, 0.7)
        for z in (0.0, 1.55, 2.0, 2.10, -1.55, -2.0, -2.10):
            for angle in angles:
                rows.append([radius * math.cos(angle), radius * math.sin(angle), z])
    points = np.asarray(rows, dtype=float)
    radii = np.hypot(points[:, 0], points[:, 1])
    outside = (radii >= 1.6) | (np.abs(points[:, 2]) >= 2.0)
    axis = radii == 0.0
    toroidal = toroidal_unit_velocity(points)
    polar = poloidal_unit_velocity(points)
    toroidal_outside_max = float(np.max(np.abs(toroidal[outside])))
    polar_outside_max = float(np.max(np.abs(polar[outside])))
    axis_transverse_max = float(np.max(np.abs(toroidal[axis, :2])))
    passes = bool(
        toroidal_outside_max <= ZERO_GATE
        and polar_outside_max <= ZERO_GATE
        and axis_transverse_max <= ZERO_GATE
        and np.all(np.isfinite(toroidal))
        and np.all(np.isfinite(polar))
    )
    return {
        "probe_count": int(len(points)),
        "outside_probe_count": int(np.sum(outside)),
        "axis_probe_count": int(np.sum(axis)),
        "toroidal_outside_max_abs": toroidal_outside_max,
        "poloidal_outside_max_abs": polar_outside_max,
        "toroidal_axis_transverse_max_abs": axis_transverse_max,
        "passes": passes,
    }


def build_report() -> dict[str, Any]:
    symbolic = symbolic_toroidal_divergence()
    selectivity = channel_selectivity()
    support = support_and_axis_check()
    identifiability = response_identifiability()
    passed = bool(
        symbolic["exact_identity"]
        and selectivity["passes"]
        and support["passes"]
        and identifiability["passes"]
    )
    return {
        "task_id": TASK_ID,
        "prereg_issue": PREREG_ISSUE,
        "source_parent": {"pr": SOURCE_PARENT_PR, "head": SOURCE_PARENT_HEAD},
        "source_poloidal_reservoir": {
            "pr": SOURCE_POLAR_RESERVOIR_PR,
            "head": SOURCE_POLAR_RESERVOIR_HEAD,
            "k": LIVE_K,
        },
        "frozen_basis": {
            "cartesian": "(-y*H(r^2,z), x*H(r^2,z), 0)",
            "H": "R_k=1(r^2) * abs(Z_outer_reservoir(z))",
            "cylindrical": "u_r=0, u_theta=r*H, u_z=0",
            "radial_support": "r<1.6",
            "axial_support": "1<|z|<2",
            "time_response_gain": "g1(t)=2*(t-.25); no new temporal basis",
        },
        "symbolic_divergence": symbolic,
        "channel_selectivity": selectivity,
        "support_and_axis": support,
        "response_identifiability": identifiability,
        "decision": {
            "frozen_preflight_passes": passed,
            "toroidal_direction_is_independent_of_live_poloidal_direction": passed,
            "instantaneous_swirl_control_without_poloidal_leakage_established": passed,
            "coefficient_selected": False,
            "actual_velocity_changed": False,
            "closer_visualization_delivery_established": False,
            "next_use_requires_preregistered_swirl_or_circulation_mismatch": True,
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
