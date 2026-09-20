"""Preflight one same-dimension axial-envelope reshape for the ST052-M tip channel.

Preregistered in issue #731 before execution and stacked directly on Agent-7
PR #723.  PR #723 proved that the frozen #701 single compact tip-local
poloidal correction must contain some radial return flow.  This increment does
not add a second basis direction.  It asks whether the unavoidable return flow
can be moved into a thinner outer collar by reshaping only the axial envelope.

With q=z^2, a=1, b=81/25 and

    B_{p,m}(q) propto (q-a)^p (b-q)^m,
    f=R(r^2) z B_{p,m}(z^2),
    A=(-y f, x f, 0),
    C=curl(A),

the sign of f_z inside the support is the sign of

    (q-a)(b-q) + 2q[p(b-q)-m(q-a)].

The outer exponent is frozen at m=4.  The inner exponent is selected by the
preregistered deterministic rule: the smallest integer p>=4 whose unique
physical sign change satisfies |z|_* >= 1.60.  No continuous scan, child
coefficient fit, trajectory replay, pressure/forcing change or visual/PDE
threshold is used here.
"""
from __future__ import annotations

import argparse
import json
import math
from fractions import Fraction
from pathlib import Path
from typing import Any

import numpy as np

import agent7_st052m_tip_odd_poloidal_screen as screen
import agent7_st052m_tip_return_flow_obstruction as obstruction
import agent7_st052m_vector_potential_poloidal_screen as vp

TASK_ID = "CR003-ST052M-TIP-RETURN-FLOW-RELOCATION-101"
PREREG_ISSUE = 731
SOURCE_PARENT_PR = 723
SOURCE_PARENT_HEAD = "1b16147b93e7b5279c833e8629e10aff1ef4f5a6"
SOURCE_BASIS_PR = 701
RELATED_FRESH_HOLDOUT_PR = 714

AXIAL_Q_A = Fraction(1, 1)
AXIAL_Q_B = Fraction(81, 25)
AXIAL_Z_INNER = 1.0
AXIAL_Z_OUTER = 1.8
OUTER_EXPONENT = 4
MIN_INNER_EXPONENT = 4
TARGET_SIGN_CHANGE_ABS_Z = 1.60
MAX_INNER_EXPONENT = 64

SPOT_RADIUS = 1.0
SHOULDER_Z_ABS = 0.85
INWARD_SPOT_Z_ABS = (1.25, 1.55)
RETURN_SPOT_Z_ABS = 1.70

TRUTH = {
    "canonical_velocity_changed": False,
    "saved_velocity_changed": False,
    "production_candidate_selected": False,
    "basis_dimension_changed": False,
    "new_spatial_basis_direction_added": False,
    "new_temporal_basis_added": False,
    "same_dimension_axial_envelope_preflighted": True,
    "nonlinear_child_evaluated": False,
    "coefficient_alpha_fitted": False,
    "trajectory_result_used": False,
    "held_out_pde_residual_evaluated": False,
    "pressure_or_force_changed": False,
    "parameter_continuous_scan_performed": False,
    "public_image_numeric_target_used": False,
    "pixel_similarity_objective_used": False,
    "visualization_ready": False,
    "visual_correspondence_verified": False,
    "pde_validated": False,
    "source_correspondence_verified": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


def _assert_source_lock() -> None:
    obstruction._assert_source_lock()
    expected = {
        "Z2_BUMP_A": 1.0,
        "Z2_BUMP_B": 3.24,
        "Z_SUPPORT_INNER": 1.0,
        "Z_SUPPORT_OUTER": 1.8,
        "R_SUPPORT_MAX": 1.6,
    }
    for name, value in expected.items():
        actual = float(getattr(screen, name))
        if actual != value:
            raise RuntimeError(f"frozen #701 source constant drifted: {name}={actual!r}")


def sign_polynomial_coefficients(p: int, m: int = OUTER_EXPONENT) -> tuple[Fraction, Fraction, Fraction]:
    """Return coefficients of the exact f_z sign polynomial in q=z^2."""
    p = int(p)
    m = int(m)
    if p < 1 or m < 1:
        raise ValueError("axial bump exponents must be positive integers")
    a, b = AXIAL_Q_A, AXIAL_Q_B
    # (q-a)(b-q) + 2q[p(b-q)-m(q-a)]
    c2 = Fraction(-(1 + 2 * (p + m)), 1)
    c1 = (a + b) + 2 * (p * b + m * a)
    c0 = -(a * b)
    return c2, c1, c0


def sign_roots(p: int, m: int = OUTER_EXPONENT) -> tuple[float, float]:
    c2, c1, c0 = sign_polynomial_coefficients(p, m)
    aa, bb, cc = float(c2), float(c1), float(c0)
    disc = bb * bb - 4.0 * aa * cc
    if disc <= 0.0:
        raise RuntimeError("asymmetric sign polynomial unexpectedly lacks two real roots")
    roots = sorted(((-bb - math.sqrt(disc)) / (2.0 * aa), (-bb + math.sqrt(disc)) / (2.0 * aa)))
    return float(roots[0]), float(roots[1])


def physical_sign_change(p: int, m: int = OUTER_EXPONENT) -> tuple[float, float]:
    roots = sign_roots(p, m)
    inside = [q for q in roots if float(AXIAL_Q_A) < q < float(AXIAL_Q_B)]
    if len(inside) != 1:
        raise RuntimeError(f"expected one physical sign root for p={p}, m={m}; got {inside}")
    q_star = float(inside[0])
    return q_star, float(math.sqrt(q_star))


def select_inner_exponent() -> tuple[int, float, float]:
    """Apply the preregistered discrete minimal-p rule; this is not a fit."""
    for p in range(MIN_INNER_EXPONENT, MAX_INNER_EXPONENT + 1):
        q_star, z_star = physical_sign_change(p, OUTER_EXPONENT)
        if z_star >= TARGET_SIGN_CHANGE_ABS_Z:
            if p > MIN_INNER_EXPONENT:
                _, previous_z = physical_sign_change(p - 1, OUTER_EXPONENT)
                if previous_z >= TARGET_SIGN_CHANGE_ABS_Z:
                    raise RuntimeError("selected p is not the minimal preregistered integer")
            return int(p), q_star, z_star
    raise RuntimeError("no preregistered integer exponent reaches the target sign location")


def _asymmetric_axial_bump(q: np.ndarray, p: int, m: int) -> tuple[np.ndarray, np.ndarray]:
    """Return normalized B(q) and dB/dq with exact compact support."""
    q = np.asarray(q, dtype=float)
    out = np.zeros_like(q)
    deriv = np.zeros_like(q)
    a, b = float(AXIAL_Q_A), float(AXIAL_Q_B)
    inside = (q > a) & (q < b)
    if not np.any(inside):
        return out, deriv

    peak = (float(p) * b + float(m) * a) / float(p + m)
    peak_raw = (peak - a) ** p * (b - peak) ** m
    if not (np.isfinite(peak_raw) and peak_raw > 0.0):
        raise RuntimeError("invalid asymmetric axial-bump normalization")

    qi = q[inside]
    left = qi - a
    right = b - qi
    raw = left**p * right**m
    raw_q = (
        float(p) * left ** (p - 1) * right**m
        - float(m) * left**p * right ** (m - 1)
    )
    out[inside] = raw / peak_raw
    deriv[inside] = raw_q / peak_raw
    return out, deriv


def reshaped_tip_correction(points: np.ndarray, p: int | None = None) -> np.ndarray:
    """Analytic curl(A) for the preregistered same-dimension reshaped channel."""
    pts = np.asarray(points, dtype=float)
    if pts.ndim != 2 or pts.shape[1] != 3:
        raise ValueError("points must have shape (n,3)")
    if p is None:
        p = select_inner_exponent()[0]
    p = int(p)

    x, y, z = pts[:, 0], pts[:, 1], pts[:, 2]
    r2 = x * x + y * y
    z2 = z * z
    radial, radial_s = vp._poly_bump(r2, screen.R2_BUMP_A, screen.R2_BUMP_B)
    axial_bump, axial_q = _asymmetric_axial_bump(z2, p, OUTER_EXPONENT)

    axial = z * axial_bump
    axial_z = axial_bump + 2.0 * z2 * axial_q
    f = radial * axial
    f_s = radial_s * axial
    f_z = radial * axial_z
    correction = np.column_stack((-x * f_z, -y * f_z, 2.0 * f + 2.0 * r2 * f_s))
    if not np.isfinite(correction).all():
        raise RuntimeError("nonfinite reshaped tip correction")
    return correction


def _radial_component_at(z: float, p: int) -> tuple[float, float]:
    point = np.asarray([[SPOT_RADIUS, 0.0, float(z)]], dtype=float)
    correction = np.asarray(reshaped_tip_correction(point, p), dtype=float)[0]
    return float(correction[0]), float(np.linalg.norm(correction))


def evaluate() -> dict[str, Any]:
    _assert_source_lock()
    baseline_q, baseline_z = obstruction.physical_sign_change()
    p, q_star, z_star = select_inner_exponent()
    roots = sign_roots(p, OUTER_EXPONENT)
    coeffs = sign_polynomial_coefficients(p, OUTER_EXPONENT)

    baseline_collar = AXIAL_Z_OUTER - baseline_z
    reshaped_collar = AXIAL_Z_OUTER - z_star
    if not (0.0 < reshaped_collar < baseline_collar):
        raise RuntimeError("preregistered reshape did not thin the outer return collar")
    collar_reduction_fraction = 1.0 - reshaped_collar / baseline_collar

    if not (AXIAL_Z_INNER < z_star < AXIAL_Z_OUTER and z_star >= TARGET_SIGN_CHANGE_ABS_Z):
        raise RuntimeError("selected sign change violates preregistered location gate")

    inward_spots: dict[str, float] = {}
    return_spots: dict[str, float] = {}
    for sign in (-1.0, 1.0):
        for abs_z in INWARD_SPOT_Z_ABS:
            z = sign * abs_z
            radial, norm = _radial_component_at(z, p)
            if not (radial < 0.0 and norm > 0.0):
                raise RuntimeError(f"reshaped channel is not inward at preregistered z={z}")
            inward_spots[f"z={z:+.2f}"] = radial
        z = sign * RETURN_SPOT_Z_ABS
        radial, norm = _radial_component_at(z, p)
        if not (radial > 0.0 and norm > 0.0):
            raise RuntimeError(f"reshaped channel lacks outer return flow at z={z}")
        return_spots[f"z={z:+.2f}"] = radial

    shoulder_checks = {}
    for sign in (-1.0, 1.0):
        z = sign * SHOULDER_Z_ABS
        radial, norm = _radial_component_at(z, p)
        if not (radial == 0.0 and norm == 0.0):
            raise RuntimeError("frozen shoulder probe is no longer exactly outside tip support")
        shoulder_checks[f"z={z:+.2f}"] = norm

    outside = np.asarray(
        [
            [1.61, 0.0, 1.25],
            [-1.61, 0.0, -1.25],
            [1.0, 0.0, 0.99],
            [1.0, 0.0, -0.99],
            [1.0, 0.0, 1.81],
            [1.0, 0.0, -1.81],
        ],
        dtype=float,
    )
    outside_max = float(np.max(np.abs(reshaped_tip_correction(outside, p))))
    if outside_max != 0.0:
        raise RuntimeError("reshaped correction leaked outside frozen support")

    # Compact-lobe return flow remains mathematically unavoidable.  The new
    # result only relocates its sign change closer to the outer support edge.
    compact_lobe_integral_identity = 0.0

    return {
        "task_id": TASK_ID,
        "prereg_issue": PREREG_ISSUE,
        "source_parent": {"pr": SOURCE_PARENT_PR, "head": SOURCE_PARENT_HEAD},
        "source_basis_pr": SOURCE_BASIS_PR,
        "related_fresh_holdout_pr": RELATED_FRESH_HOLDOUT_PR,
        "frozen_rule": {
            "family": "B_p_m(q) proportional to (q-a)^p*(b-q)^m",
            "a": float(AXIAL_Q_A),
            "b": float(AXIAL_Q_B),
            "outer_exponent_m": OUTER_EXPONENT,
            "minimum_inner_exponent_p": MIN_INNER_EXPONENT,
            "target_sign_change_abs_z": TARGET_SIGN_CHANGE_ABS_Z,
            "selection": "smallest integer p>=4 with physical |z| sign change >=1.60",
            "continuous_parameter_scan": False,
        },
        "selected_envelope": {
            "inner_exponent_p": p,
            "outer_exponent_m": OUTER_EXPONENT,
            "sign_polynomial_coefficients": [str(v) for v in coeffs],
            "q_roots": [float(v) for v in roots],
            "physical_q_sign_change": q_star,
            "physical_abs_z_sign_change": z_star,
            "inward_correction_subband_abs_z": [AXIAL_Z_INNER, z_star],
            "outward_return_subband_abs_z": [z_star, AXIAL_Z_OUTER],
            "outer_zero_extension_exponent_not_weakened": bool(OUTER_EXPONENT >= 4),
        },
        "baseline_comparison": {
            "baseline_symmetric_quartic_q_sign_change": baseline_q,
            "baseline_symmetric_quartic_abs_z_sign_change": baseline_z,
            "baseline_return_collar_thickness": baseline_collar,
            "reshaped_return_collar_thickness": reshaped_collar,
            "return_collar_reduction_fraction": collar_reduction_fraction,
        },
        "implementation_spot_checks": {
            "inward_radial_components": inward_spots,
            "outward_return_radial_components": return_spots,
            "shoulder_correction_norms": shoulder_checks,
            "outside_support_correction_max_abs": outside_max,
        },
        "structural_scope": {
            "compact_lobe_integral_Cr_exact": compact_lobe_integral_identity,
            "return_flow_removed": False,
            "return_flow_relocated_outward": True,
            "full_total_child_tip_inwardness_established": False,
        },
        "decision": {
            "classification": "same_dimension_return_flow_relocation_preflight_pass",
            "same_dimension_envelope_reshape_justified_for_one_nonlinear_replay": True,
            "second_poloidal_channel_justified_now": False,
            "next_minimal_increment": (
                "freeze one exact nonlinear #701-derived child using this p=9,m=4 envelope, "
                "then replay the existing radial/turns/axial-stretch guards before adding "
                "any second poloidal channel"
            ),
            "velocity_changed_this_increment": False,
            "closer_visualization_delivery_this_increment": False,
        },
        **TRUTH,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    report = evaluate()
    text = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.out is None:
        print(text, end="")
        return
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(text, encoding="utf-8")
    print(args.out)


if __name__ == "__main__":
    main()
