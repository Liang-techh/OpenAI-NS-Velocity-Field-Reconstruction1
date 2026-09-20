"""Audit the vorticity-morphology response of the frozen #846 toroidal child.

Preregistered in issue #854 before execution and stacked directly on exact
Agent-7 PR #846 head.  This is one minimal expression-capacity increment:
it changes no velocity, coefficient, pressure, forcing, renderer, or source
target.  It asks whether the already-frozen compact toroidal velocity channel
adds a complementary vorticity direction relative to the live outer-reservoir
poloidal channel, and whether the actual #846 child realizes that response.

For a generic axisymmetric H(q,z), q=x^2+y^2,

    C_theta = (-y H, x H, 0)

has exact curl

    curl C_theta = (-x H_z, -y H_z, 2 H + 2 q H_q).

Thus its vorticity is poloidal (omega_r,0,omega_z).  Conversely an
axisymmetric poloidal velocity (u_r,0,u_z) has purely azimuthal vorticity.
The numerical checks below verify those representation claims on the live
implementations and quantify the total-child response without inventing an
OpenAI numerical vorticity target.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any, Callable

import numpy as np
import sympy as sp

import agent7_st052m_toroidal_child_morphology_replay as child

TASK_ID = "CR003-ST052M-TOROIDAL-VORTICITY-RESPONSE-AUDIT-121"
PREREG_ISSUE = 854
SOURCE_PARENT_PR = 846
SOURCE_PARENT_HEAD = "55ea4f1ca34d81b386b1f9d647cfaae20b1863e0"
SOURCE_TOROIDAL_PREFLIGHT_PR = 837
SOURCE_TOROIDAL_PREFLIGHT_HEAD = "4c821c302cf89a7cc7a41730ce7029bdfbe6a7ef"
SOURCE_RESERVOIR_CHILD_PR = 775
SOURCE_RESERVOIR_CHILD_HEAD = "93887a59729d22113badf2ae4dac2f7d868e2703"
RELATED_FRESH_HOLDOUT_PR = 714

PROBE_RADII = (0.35, 0.70, 1.05, 1.35)
PROBE_ABS_Z = (1.20, 1.55, 1.75, 1.90)
PROBE_ANGLES = tuple(float(k) * math.pi / 4.0 for k in range(8))
AUDIT_TIME = 0.50
CURL_STEP = 2.0e-5
LEAKAGE_GATE = 2.0e-8
NONTRIVIAL_RMS_MIN = 1.0e-7
COSINE_MAX_ABS = 2.0e-6
CONDITION_MAX = 1.01
DELTA_REL_RMS_MAX = 2.0e-5
DELTA_REL_SAMPLED_MAX = 5.0e-5
RELATIVE_DENOMINATOR_FLOOR = 1.0e-10

TRUTH = {
    "candidate_velocity_changed": False,
    "canonical_velocity_changed": False,
    "saved_velocity_changed": False,
    "new_basis_added": False,
    "new_coefficient_selected": False,
    "existing_846_child_audited": True,
    "coefficient_scan_performed": False,
    "optimization_performed": False,
    "post_result_retuning_performed": False,
    "fresh_714_path_data_used": False,
    "pressure_or_force_changed": False,
    "held_out_pde_residual_evaluated": False,
    "st006_comparison_performed": False,
    "public_image_numeric_target_used": False,
    "source_numeric_vorticity_target_used": False,
    "renderer_or_camera_used": False,
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

SpatialVelocityFn = Callable[[np.ndarray], np.ndarray]


def _assert_source_lock() -> None:
    child._assert_source_lock()
    if child.TASK_ID != "CR003-ST052M-TOROIDAL-CHILD-MORPHOLOGY-REPLAY-120":
        raise RuntimeError("#846 parent task identity drifted")
    if child.PREREG_ISSUE != 845:
        raise RuntimeError("#846 parent preregistration drifted")
    if child.SOURCE_PARENT_PR != SOURCE_TOROIDAL_PREFLIGHT_PR:
        raise RuntimeError("#837 source PR drifted")
    if child.SOURCE_PARENT_HEAD != SOURCE_TOROIDAL_PREFLIGHT_HEAD:
        raise RuntimeError("#837 source head drifted")
    if child.SOURCE_RESERVOIR_CHILD_PR != SOURCE_RESERVOIR_CHILD_PR:
        raise RuntimeError("#775 source PR drifted")
    if child.SOURCE_RESERVOIR_CHILD_HEAD != SOURCE_RESERVOIR_CHILD_HEAD:
        raise RuntimeError("#775 source head drifted")
    if child.RELATED_FRESH_HOLDOUT_PR != RELATED_FRESH_HOLDOUT_PR:
        raise RuntimeError("fresh holdout identity drifted")


def _probe_points() -> np.ndarray:
    rows: list[list[float]] = []
    for radius in PROBE_RADII:
        for abs_z in PROBE_ABS_Z:
            for sign in (-1.0, 1.0):
                for angle in PROBE_ANGLES:
                    rows.append(
                        [
                            float(radius * math.cos(angle)),
                            float(radius * math.sin(angle)),
                            float(sign * abs_z),
                        ]
                    )
    points = np.asarray(rows, dtype=float)
    expected = len(PROBE_RADII) * len(PROBE_ABS_Z) * 2 * len(PROBE_ANGLES)
    if points.shape != (expected, 3):
        raise RuntimeError("frozen vorticity probe cloud drifted")
    radii = np.hypot(points[:, 0], points[:, 1])
    if np.any(radii + CURL_STEP >= 1.6):
        raise RuntimeError("radial curl stencil escaped compact support")
    if np.any(np.abs(points[:, 2]) + CURL_STEP >= 2.0):
        raise RuntimeError("axial curl stencil escaped compact support")
    if np.any(np.abs(points[:, 2]) - CURL_STEP <= 1.0):
        raise RuntimeError("axial curl stencil crossed inner support boundary")
    return points


def symbolic_toroidal_curl_identity() -> dict[str, Any]:
    x, y, z, q_symbol = sp.symbols("x y z q_symbol", real=True, finite=True)
    q = x**2 + y**2
    h_function = sp.Function("H")
    h = h_function(q, z)
    velocity = sp.Matrix([-y * h, x * h, 0])
    omega = sp.Matrix(
        [
            sp.diff(velocity[2], y) - sp.diff(velocity[1], z),
            sp.diff(velocity[0], z) - sp.diff(velocity[2], x),
            sp.diff(velocity[1], x) - sp.diff(velocity[0], y),
        ]
    )
    h_q = sp.Subs(sp.Derivative(h_function(q_symbol, z), q_symbol), q_symbol, q)
    target = sp.Matrix(
        [
            -x * sp.diff(h, z),
            -y * sp.diff(h, z),
            2 * h + 2 * q * h_q,
        ]
    )
    residual = sp.Matrix([sp.simplify(omega[i] - target[i]) for i in range(3)])
    exact = all(value == 0 for value in residual)
    return {
        "velocity": "(-y*H(q,z), x*H(q,z), 0), q=x^2+y^2",
        "curl": "(-x*H_z, -y*H_z, 2*H+2*q*H_q)",
        "exact_identity": bool(exact),
        "residual": [str(v) for v in residual],
    }


def curl_fd(velocity_fn: SpatialVelocityFn, points: np.ndarray, step: float = CURL_STEP) -> np.ndarray:
    pts = np.asarray(points, dtype=float)
    if pts.ndim != 2 or pts.shape[1] != 3:
        raise ValueError("points must have shape (n,3)")
    if not np.all(np.isfinite(pts)):
        raise ValueError("points must be finite")
    h = float(step)
    if not np.isfinite(h) or h <= 0.0:
        raise ValueError("step must be positive and finite")

    derivs: list[np.ndarray] = []
    for axis in range(3):
        shift = np.zeros_like(pts)
        shift[:, axis] = h
        plus = np.asarray(velocity_fn(pts + shift), dtype=float)
        minus = np.asarray(velocity_fn(pts - shift), dtype=float)
        if plus.shape != pts.shape or minus.shape != pts.shape:
            raise RuntimeError("velocity function returned wrong shape")
        deriv = (plus - minus) / (2.0 * h)
        if not np.all(np.isfinite(deriv)):
            raise RuntimeError("nonfinite curl derivative")
        derivs.append(deriv)

    dx, dy, dz = derivs
    omega = np.column_stack(
        (
            dy[:, 2] - dz[:, 1],
            dz[:, 0] - dx[:, 2],
            dx[:, 1] - dy[:, 0],
        )
    )
    if not np.all(np.isfinite(omega)):
        raise RuntimeError("nonfinite vorticity")
    return omega


def _rms(values: np.ndarray) -> float:
    arr = np.asarray(values, dtype=float)
    return float(np.sqrt(np.mean(arr * arr)))


def _vector_rms(vectors: np.ndarray) -> float:
    vec = np.asarray(vectors, dtype=float)
    return float(np.sqrt(np.mean(np.sum(vec * vec, axis=1))))


def _cylindrical(points: np.ndarray, vectors: np.ndarray) -> np.ndarray:
    return np.asarray(child.toroidal.cylindrical_components(points, vectors), dtype=float)


def unit_vorticity_selectivity() -> dict[str, Any]:
    _assert_source_lock()
    points = _probe_points()
    toroidal_omega = curl_fd(child.toroidal.toroidal_unit_velocity, points)
    poloidal_omega = curl_fd(child.toroidal.poloidal_unit_velocity, points)
    toroidal_cyl = _cylindrical(points, toroidal_omega)
    poloidal_cyl = _cylindrical(points, poloidal_omega)

    toroidal_theta_leak = float(np.max(np.abs(toroidal_cyl[:, 1])))
    toroidal_poloidal_rms = float(
        np.sqrt(np.mean(toroidal_cyl[:, 0] ** 2 + toroidal_cyl[:, 2] ** 2))
    )
    poloidal_rz_leak = float(
        np.max(np.sqrt(poloidal_cyl[:, 0] ** 2 + poloidal_cyl[:, 2] ** 2))
    )
    poloidal_theta_rms = _rms(poloidal_cyl[:, 1])

    tor_col = toroidal_omega.ravel()
    pol_col = poloidal_omega.ravel()
    tor_norm = float(np.linalg.norm(tor_col))
    pol_norm = float(np.linalg.norm(pol_col))
    if tor_norm <= 0.0 or pol_norm <= 0.0:
        raise RuntimeError("degenerate unit-vorticity response")
    tor_normalized = tor_col / tor_norm
    pol_normalized = pol_col / pol_norm
    matrix = np.column_stack((pol_normalized, tor_normalized))
    singular = np.linalg.svd(matrix, compute_uv=False)
    rank = int(np.linalg.matrix_rank(matrix, tol=1.0e-12))
    cosine = float(np.dot(pol_normalized, tor_normalized))
    condition = float(singular[0] / singular[-1])

    selectivity_pass = bool(
        toroidal_theta_leak <= LEAKAGE_GATE
        and toroidal_poloidal_rms > NONTRIVIAL_RMS_MIN
        and poloidal_rz_leak <= LEAKAGE_GATE
        and poloidal_theta_rms > NONTRIVIAL_RMS_MIN
    )
    identifiability_pass = bool(
        rank == 2 and abs(cosine) <= COSINE_MAX_ABS and condition <= CONDITION_MAX
    )
    return {
        "probe_count": int(len(points)),
        "curl_step": CURL_STEP,
        "toroidal": {
            "omega_theta_max_abs": toroidal_theta_leak,
            "omega_rz_rms": toroidal_poloidal_rms,
            "cartesian_vorticity_rms": _vector_rms(toroidal_omega),
        },
        "poloidal": {
            "omega_rz_max_abs": poloidal_rz_leak,
            "omega_theta_rms": poloidal_theta_rms,
            "cartesian_vorticity_rms": _vector_rms(poloidal_omega),
        },
        "normalized_response": {
            "rank": rank,
            "cosine": cosine,
            "condition": condition,
            "singular_values": [float(v) for v in singular],
            "cosine_abs_gate": COSINE_MAX_ABS,
            "condition_gate": CONDITION_MAX,
            "passes": identifiability_pass,
        },
        "leakage_gate": LEAKAGE_GATE,
        "nontrivial_rms_min": NONTRIVIAL_RMS_MIN,
        "selectivity_passes": selectivity_pass,
        "passes": bool(selectivity_pass and identifiability_pass),
    }


def _component_rms(cylindrical: np.ndarray) -> dict[str, float]:
    cyl = np.asarray(cylindrical, dtype=float)
    return {
        "omega_r_rms": _rms(cyl[:, 0]),
        "omega_theta_rms": _rms(cyl[:, 1]),
        "omega_z_rms": _rms(cyl[:, 2]),
        "omega_total_rms": _vector_rms(cyl),
    }


def actual_child_vorticity_response() -> dict[str, Any]:
    _assert_source_lock()
    points = _probe_points()
    parent_fn, alpha_res, alpha_calibration = child._build_reservoir_parent()
    beta, beta_calibration = child.derive_beta(parent_fn)

    def parent_spatial(pts: np.ndarray) -> np.ndarray:
        return np.asarray(parent_fn(pts, AUDIT_TIME), dtype=float)

    def child_spatial(pts: np.ndarray) -> np.ndarray:
        return np.asarray(child.child_velocity(parent_fn, pts, AUDIT_TIME, beta), dtype=float)

    parent_omega = curl_fd(parent_spatial, points)
    child_omega = curl_fd(child_spatial, points)
    unit_toroidal_omega = curl_fd(child.toroidal.toroidal_unit_velocity, points)
    gain = float(child.toroidal.activation_g1(AUDIT_TIME))
    predicted_delta = float(beta) * gain * unit_toroidal_omega
    observed_delta = child_omega - parent_omega
    error = observed_delta - predicted_delta

    pred_rms = max(_vector_rms(predicted_delta), RELATIVE_DENOMINATOR_FLOOR)
    pred_sample_norm = np.linalg.norm(predicted_delta, axis=1)
    err_sample_norm = np.linalg.norm(error, axis=1)
    pred_max = max(float(np.max(pred_sample_norm)), RELATIVE_DENOMINATOR_FLOOR)
    relative_rms = float(_vector_rms(error) / pred_rms)
    relative_sampled_max = float(np.max(err_sample_norm) / pred_max)

    parent_cyl = _cylindrical(points, parent_omega)
    child_cyl = _cylindrical(points, child_omega)
    delta_cyl = _cylindrical(points, observed_delta)
    predicted_cyl = _cylindrical(points, predicted_delta)

    passes = bool(
        beta_calibration["nontrivial"]
        and relative_rms <= DELTA_REL_RMS_MAX
        and relative_sampled_max <= DELTA_REL_SAMPLED_MAX
    )
    return {
        "time": AUDIT_TIME,
        "probe_count": int(len(points)),
        "alpha_res": float(alpha_res),
        "alpha_res_calibration_fresh_714_data_used": bool(
            alpha_calibration.get(
                "fresh_714_data_used",
                alpha_calibration.get("fresh_714_path_data_used", False),
            )
        ),
        "beta_theta": float(beta),
        "beta_calibration": beta_calibration,
        "time_gain_g1": gain,
        "relative_delta_consistency": {
            "relative_rms": relative_rms,
            "relative_sampled_max": relative_sampled_max,
            "relative_rms_gate": DELTA_REL_RMS_MAX,
            "relative_sampled_max_gate": DELTA_REL_SAMPLED_MAX,
            "denominator_floor": RELATIVE_DENOMINATOR_FLOOR,
            "passes": bool(
                relative_rms <= DELTA_REL_RMS_MAX
                and relative_sampled_max <= DELTA_REL_SAMPLED_MAX
            ),
        },
        "parent_total_vorticity": _component_rms(parent_cyl),
        "child_total_vorticity": _component_rms(child_cyl),
        "observed_child_minus_parent": _component_rms(delta_cyl),
        "predicted_toroidal_increment": _component_rms(predicted_cyl),
        "observed_delta_omega_theta_max_abs": float(np.max(np.abs(delta_cyl[:, 1]))),
        "source_numeric_vorticity_target_used": False,
        "fresh_714_data_used": False,
        "held_out_pde_residual_evaluated": False,
        "passes": passes,
    }


def build_report() -> dict[str, Any]:
    symbolic = symbolic_toroidal_curl_identity()
    selectivity = unit_vorticity_selectivity()
    actual = actual_child_vorticity_response()
    passed = bool(symbolic["exact_identity"] and selectivity["passes"] and actual["passes"])
    return {
        "task_id": TASK_ID,
        "prereg_issue": PREREG_ISSUE,
        "source_parent": {"pr": SOURCE_PARENT_PR, "head": SOURCE_PARENT_HEAD},
        "source_toroidal_preflight": {
            "pr": SOURCE_TOROIDAL_PREFLIGHT_PR,
            "head": SOURCE_TOROIDAL_PREFLIGHT_HEAD,
        },
        "source_reservoir_child": {
            "pr": SOURCE_RESERVOIR_CHILD_PR,
            "head": SOURCE_RESERVOIR_CHILD_HEAD,
        },
        "frozen_protocol": {
            "probe_radii": list(PROBE_RADII),
            "probe_abs_z": list(PROBE_ABS_Z),
            "probe_azimuth_count": len(PROBE_ANGLES),
            "audit_time": AUDIT_TIME,
            "curl_step": CURL_STEP,
            "leakage_gate": LEAKAGE_GATE,
            "nontrivial_rms_min": NONTRIVIAL_RMS_MIN,
            "cosine_abs_gate": COSINE_MAX_ABS,
            "condition_gate": CONDITION_MAX,
            "delta_relative_rms_gate": DELTA_REL_RMS_MAX,
            "delta_relative_sampled_max_gate": DELTA_REL_SAMPLED_MAX,
            "post_result_retuning_allowed": False,
        },
        "symbolic_toroidal_curl": symbolic,
        "unit_vorticity_selectivity": selectivity,
        "actual_846_child_vorticity_response": actual,
        "held_out_pde_residual": {
            "evaluated": False,
            "reason": (
                "#846 changes velocity without rebuilding matched pressure/forcing; "
                "this increment does not manufacture an unmatched pseudo-residual"
            ),
            "st006_comparison_performed": False,
        },
        "truth": dict(TRUTH),
        "pressure_or_force_changed": False,
        "visual_correspondence_verified": False,
        "source_correspondence_verified": False,
        "pde_validated": False,
        "paper_exact": False,
        "openai_field_identified": False,
        "toroidal_vorticity_response_audit_passed": passed,
        "visualization_delivery_effect": (
            "No new velocity is introduced in this increment. PASS would show that the already-frozen "
            "#846 toroidal child controls a vorticity direction complementary to the live poloidal "
            "reservoir channel, improving morphology controllability but not proving closer OpenAI correspondence."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    report = build_report()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    if not report["toroidal_vorticity_response_audit_passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
