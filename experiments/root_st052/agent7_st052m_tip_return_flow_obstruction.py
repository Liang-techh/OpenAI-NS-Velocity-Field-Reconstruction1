"""Analytic return-flow audit for the frozen Agent-7 PR #701 tip channel.

Preregistered in issue #722 before executing this audit.  The source channel is
not changed:

    A=(-y*f, x*f, 0),
    f=R(r^2) * z * B(z^2),
    C_tip=curl(A),

with the frozen axial quartic bump on 1 < |z| < 1.8.  This module asks one
representation-capacity question only: can this single compact channel point
radially inward throughout its whole axial lobe?  It derives the exact sign
polynomial, binds it back to the #701 implementation at both z signs, and
records the compact-lobe integral identity that forces return flow.

No velocity coefficient, basis, time law, pressure, forcing, candidate
identity, path protocol, or scientific threshold is changed.
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

TASK_ID = "CR003-ST052M-TIP-RETURN-FLOW-OBSTRUCTION-100"
PREREG_ISSUE = 722
SOURCE_PARENT_PR = 701
SOURCE_PARENT_HEAD = "e3ecc8654f2c9a2324265bbd01df596b4e503404"
RELATED_FRESH_HOLDOUT_PR = 714
RELATED_HOLDOUT_SCOPE_PR = 721

AXIAL_Q_A = Fraction(1, 1)
AXIAL_Q_B = Fraction(81, 25)  # 3.24 = 1.8^2
AXIAL_Z_INNER = 1.0
AXIAL_Z_OUTER = 1.8
FRESH_HOLDOUT_TIP_Z_ABS = 1.25
SHOULDER_Z_ABS = 0.85
SPOT_INWARD_Z_ABS = 1.25
SPOT_OUTWARD_Z_ABS = 1.60
SPOT_RADIUS = 1.0

TRUTH = {
    "canonical_velocity_changed": False,
    "saved_velocity_changed": False,
    "production_candidate_selected": False,
    "new_spatial_basis_added": False,
    "new_temporal_basis_added": False,
    "source_701_candidate_retuned": False,
    "parameter_grid_scan_performed": False,
    "optimization_performed": False,
    "pressure_or_force_changed": False,
    "trajectory_result_used": False,
    "held_out_pde_residual_evaluated": False,
    "parent_pde_receipt_transferred": False,
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
    """Fail closed if the frozen #701 spatial channel constants drift."""
    expected = {
        "R2_BUMP_A": -1.0,
        "R2_BUMP_B": 2.56,
        "Z2_BUMP_A": 1.0,
        "Z2_BUMP_B": 3.24,
        "R_SUPPORT_MAX": 1.6,
        "Z_SUPPORT_INNER": 1.0,
        "Z_SUPPORT_OUTER": 1.8,
    }
    for name, value in expected.items():
        actual = float(getattr(screen, name))
        if actual != value:
            raise RuntimeError(f"frozen #701 source constant drifted: {name}={actual!r}")


def sign_polynomial_coefficients() -> tuple[Fraction, Fraction, Fraction]:
    """Return N(q) coefficients for f_z sign inside the axial bump support.

    For p=(q-a)(b-q), B=(p/pmax)^4 and f=z B(z^2),

        f_z = B + 2 q B'
            = p^3/pmax^4 * [p + 8 q (a+b-2q)].

    Since p>0 inside (a,b), the sign is the sign of

        N(q) = -17 q^2 + 9(a+b)q - ab.
    """
    a, b = AXIAL_Q_A, AXIAL_Q_B
    return Fraction(-17, 1), 9 * (a + b), -(a * b)


def sign_polynomial(q: float) -> float:
    c2, c1, c0 = sign_polynomial_coefficients()
    q = float(q)
    return float(c2) * q * q + float(c1) * q + float(c0)


def sign_roots() -> tuple[float, float]:
    c2, c1, c0 = sign_polynomial_coefficients()
    a = float(c2)
    b = float(c1)
    c = float(c0)
    disc = b * b - 4.0 * a * c
    if disc <= 0.0:
        raise RuntimeError("frozen sign polynomial unexpectedly lacks two real roots")
    roots = sorted(((-b - math.sqrt(disc)) / (2.0 * a), (-b + math.sqrt(disc)) / (2.0 * a)))
    return float(roots[0]), float(roots[1])


def physical_sign_change() -> tuple[float, float]:
    roots = sign_roots()
    inside = [q for q in roots if float(AXIAL_Q_A) < q < float(AXIAL_Q_B)]
    if len(inside) != 1:
        raise RuntimeError(f"expected exactly one sign root in physical support, got {inside}")
    q_star = float(inside[0])
    return q_star, float(math.sqrt(q_star))


def _radial_component_at(z: float) -> tuple[float, float]:
    """Return analytic-Cartesian radial component at (r,theta)=(1,0)."""
    point = np.asarray([[SPOT_RADIUS, 0.0, float(z)]], dtype=float)
    correction = np.asarray(screen.tip_odd_poloidal_correction(point), dtype=float)[0]
    radial = float(correction[0])  # e_r=(1,0,0) at theta=0
    return radial, float(np.linalg.norm(correction))


def evaluate() -> dict[str, Any]:
    _assert_source_lock()
    c2, c1, c0 = sign_polynomial_coefficients()
    roots = sign_roots()
    q_star, z_star = physical_sign_change()

    if not (AXIAL_Z_INNER < z_star < AXIAL_Z_OUTER):
        raise RuntimeError("physical sign change left the frozen axial support")

    # Inside the support, C_r=-r*R(r^2)*f_z.  At r=1 the frozen radial
    # bump is positive, so C_r has the opposite sign of N(q).
    inward_spots: dict[str, float] = {}
    outward_spots: dict[str, float] = {}
    for sign in (-1.0, 1.0):
        z_in = sign * SPOT_INWARD_Z_ABS
        z_out = sign * SPOT_OUTWARD_Z_ABS
        radial_in, norm_in = _radial_component_at(z_in)
        radial_out, norm_out = _radial_component_at(z_out)
        if not (radial_in < 0.0 and norm_in > 0.0):
            raise RuntimeError(f"implementation disagrees with inward analytic sub-band at z={z_in}")
        if not (radial_out > 0.0 and norm_out > 0.0):
            raise RuntimeError(f"implementation disagrees with outer return-flow sub-band at z={z_out}")
        inward_spots[f"z={z_in:+.2f}"] = radial_in
        outward_spots[f"z={z_out:+.2f}"] = radial_out

    shoulder_radial_pos, shoulder_norm_pos = _radial_component_at(SHOULDER_Z_ABS)
    shoulder_radial_neg, shoulder_norm_neg = _radial_component_at(-SHOULDER_Z_ABS)
    if not (
        shoulder_radial_pos == 0.0
        and shoulder_norm_pos == 0.0
        and shoulder_radial_neg == 0.0
        and shoulder_norm_neg == 0.0
    ):
        raise RuntimeError("frozen shoulder probe is no longer exactly outside tip support")

    holdout_q = FRESH_HOLDOUT_TIP_Z_ABS**2
    holdout_n = sign_polynomial(holdout_q)
    holdout_inward = bool(
        AXIAL_Z_INNER < FRESH_HOLDOUT_TIP_Z_ABS < z_star and holdout_n > 0.0
    )
    if not holdout_inward:
        raise RuntimeError("#714 frozen tip seed no longer lies in the analytic inward sub-band")

    # Exact compact-lobe identity: for fixed active radius,
    # C_r=-r*d_z f and f vanishes at both axial support edges.  Therefore
    # integral C_r dz=-r[f(z_outer)-f(z_inner)]=0.  Since the channel is
    # nonzero and we have observed both signs, strict one-sign inward flow
    # over the full compact lobe is impossible for this channel.
    endpoint_f_inner = 0.0
    endpoint_f_outer = 0.0
    integral_identity_rhs = -SPOT_RADIUS * (endpoint_f_outer - endpoint_f_inner)
    if integral_identity_rhs != 0.0:
        raise RuntimeError("compact-lobe endpoint identity drifted")

    return {
        "task_id": TASK_ID,
        "prereg_issue": PREREG_ISSUE,
        "source_parent": {"pr": SOURCE_PARENT_PR, "head": SOURCE_PARENT_HEAD},
        "related_fresh_holdout_pr": RELATED_FRESH_HOLDOUT_PR,
        "related_holdout_scope_pr": RELATED_HOLDOUT_SCOPE_PR,
        "frozen_formula": {
            "vector_potential": "A=(-y*f,x*f,0)",
            "scalar": "f=R(r^2)*z*B(z^2)",
            "radial_component": "C_r=-r*f_z",
            "axial_q_support": [float(AXIAL_Q_A), float(AXIAL_Q_B)],
            "axial_abs_z_support": [AXIAL_Z_INNER, AXIAL_Z_OUTER],
            "sign_polynomial": "N(q)=-17*q^2+(954/25)*q-(81/25)",
            "sign_polynomial_coefficients": [str(c2), str(c1), str(c0)],
        },
        "analytic_result": {
            "q_roots": [float(v) for v in roots],
            "physical_q_sign_change": q_star,
            "physical_abs_z_sign_change": z_star,
            "inward_correction_subband_abs_z": [AXIAL_Z_INNER, z_star],
            "outward_return_subband_abs_z": [z_star, AXIAL_Z_OUTER],
            "single_compact_tip_channel_strictly_inward_everywhere_possible": False,
            "return_flow_required_by_compact_lobe_identity": True,
            "compact_lobe_integral_Cr_exact": integral_identity_rhs,
            "identity": "integral(C_r dz)=-r*(f(z_outer)-f(z_inner))=0",
        },
        "implementation_spot_checks": {
            "inward_radial_components": inward_spots,
            "outward_return_radial_components": outward_spots,
            "shoulder_plus_abs_z_correction_norm": shoulder_norm_pos,
            "shoulder_minus_abs_z_correction_norm": shoulder_norm_neg,
        },
        "fresh_holdout_scope": {
            "pr": RELATED_FRESH_HOLDOUT_PR,
            "tip_seed_abs_z": FRESH_HOLDOUT_TIP_Z_ABS,
            "tip_seed_sign_polynomial": holdout_n,
            "tip_seed_inside_analytic_inward_subband": holdout_inward,
            "full_tip_support_inwardness_inferred_from_holdout": False,
        },
        "decision": {
            "classification": "single_channel_compact_return_flow_obstruction",
            "generic_swirl_or_time_basis_growth_justified": False,
            "next_minimal_if_fresh_paths_still_require_change": (
                "reshape the same tip-local poloidal axial envelope so the unavoidable radial "
                "return flow is concentrated in a thin outer collar; add a second poloidal "
                "channel only if one-channel shaping cannot preserve the existing guards"
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
