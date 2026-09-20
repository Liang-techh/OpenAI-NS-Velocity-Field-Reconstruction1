"""Audit whether radial localization alone can move compact tip return-flow debt.

Preregistered in issue #758 and stacked on Agent-7 PR #750.  This is an
expression-capacity diagnostic only.  It changes no candidate velocity.

For the axisymmetric vector potential A=(-y F, x F, 0), the cylindrical radial
correction is C_r=-r*dF/dz.  If every axial channel vanishes at the same compact
lobe endpoints, then at every fixed radius

    integral C_r dz = -r [F(z_out)-F(z_in)] = 0.

The live numerical witness deliberately uses two linearly independent channels:
a core-connected p=9 channel and an annular p=4 channel.  The audit checks that
independence does not remove the fixed-radius signed-balance identity.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

import agent7_st052m_tip_odd_poloidal_screen as screen
import agent7_st052m_tip_return_flow_relocation as relocation
import agent7_st052m_vector_potential_poloidal_screen as vp

TASK_ID = "CR003-ST052M-FIXED-RADIUS-RETURN-BALANCE-AUDIT-106"
ISSUE = 758
SOURCE_PARENT_PR = 750
SOURCE_PARENT_HEAD = "376223a3bb2dae0d53756ea48da26937a95b6e58"
SOURCE_RESHAPE_PR = 732
SOURCE_OBSTRUCTION_PR = 723
RELATED_NONLINEAR_PR = 740
RELATED_FRESH_HOLDOUT_PR = 714

Z_INNER = 1.0
Z_OUTER = 1.8
CORE_RADIAL_A = -1.0
CORE_RADIAL_B = 2.56
ANNULAR_INNER_RADIUS = 1.05
ANNULAR_RADIAL_A = ANNULAR_INNER_RADIUS**2
ANNULAR_RADIAL_B = 2.56
CORE_AXIAL_P = 9
ANNULAR_AXIAL_P = 4
AXIAL_M = 4
RADII = (0.55, 0.95, 1.30, 1.50)
GRID_SIZES = (5001, 20001, 80001)
COEFFICIENT_WITNESSES = ((1.0, 0.0), (0.0, 1.0), (1.0, 1.0), (1.0, -0.7), (2.3, -1.1))
SIGNED_INTEGRAL_TOL = 1.0e-10
REFINEMENT_CHANGE_TOL = 1.0e-10
RANK_RADIAL_N = 48
RANK_AXIAL_N = 257
RANK_COSINE_MAX_ABS = 0.999

TRUTH = {
    "canonical_velocity_changed": False,
    "saved_velocity_changed": False,
    "candidate_velocity_changed": False,
    "basis_dimension_changed": False,
    "second_poloidal_basis_added_to_candidate": False,
    "new_temporal_basis_added": False,
    "coefficient_fit_performed": False,
    "trajectory_replay_performed": False,
    "fresh_714_path_data_used": False,
    "held_out_pde_residual_evaluated": False,
    "pressure_or_force_changed": False,
    "public_image_numeric_target_used": False,
    "visualization_ready": False,
    "visual_correspondence_verified": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


def _assert_source_lock() -> None:
    relocation._assert_source_lock()
    expected = {
        "R2_BUMP_A": CORE_RADIAL_A,
        "R2_BUMP_B": CORE_RADIAL_B,
        "Z2_BUMP_A": Z_INNER**2,
        "Z2_BUMP_B": Z_OUTER**2,
    }
    for name, value in expected.items():
        actual = float(getattr(screen, name))
        if actual != float(value):
            raise RuntimeError(f"frozen source constant drifted: {name}={actual!r}")
    if relocation.select_inner_exponent()[0] != CORE_AXIAL_P:
        raise RuntimeError("frozen p=9 reshape identity drifted")


def _radial_profile(r: np.ndarray, channel: int) -> np.ndarray:
    r = np.asarray(r, dtype=float)
    s = r * r
    if int(channel) == 0:
        return vp._poly_bump(s, CORE_RADIAL_A, CORE_RADIAL_B)[0]
    if int(channel) == 1:
        return vp._poly_bump(s, ANNULAR_RADIAL_A, ANNULAR_RADIAL_B)[0]
    raise ValueError("channel must be 0 or 1")


def _axial_bump_and_q(z: np.ndarray, p: int) -> tuple[np.ndarray, np.ndarray]:
    z = np.asarray(z, dtype=float)
    return relocation._asymmetric_axial_bump(z * z, int(p), AXIAL_M)


def _axial_potential(z: np.ndarray, p: int) -> np.ndarray:
    z = np.asarray(z, dtype=float)
    bump, _ = _axial_bump_and_q(z, int(p))
    return z * bump


def _axial_response(z: np.ndarray, p: int) -> np.ndarray:
    z = np.asarray(z, dtype=float)
    bump, bump_q = _axial_bump_and_q(z, int(p))
    return bump + 2.0 * z * z * bump_q


def channel_radial_response(r: float | np.ndarray, z: np.ndarray, channel: int) -> np.ndarray:
    """Return C_r for a unit coefficient at fixed/array radius."""
    rr = np.asarray(r, dtype=float)
    zz = np.asarray(z, dtype=float)
    if int(channel) == 0:
        p = CORE_AXIAL_P
    elif int(channel) == 1:
        p = ANNULAR_AXIAL_P
    else:
        raise ValueError("channel must be 0 or 1")
    return -rr * _radial_profile(rr, int(channel)) * _axial_response(zz, p)


def endpoint_potential(r: float, z: float, channel: int) -> float:
    if int(channel) == 0:
        p = CORE_AXIAL_P
    elif int(channel) == 1:
        p = ANNULAR_AXIAL_P
    else:
        raise ValueError("channel must be 0 or 1")
    return float(_radial_profile(np.asarray([r]), int(channel))[0] * _axial_potential(np.asarray([z]), p)[0])


def _fixed_radius_integrals(n: int) -> dict[str, Any]:
    z = np.linspace(Z_INNER, Z_OUTER, int(n))
    rows: dict[str, Any] = {}
    for radius in RADII:
        c0 = np.asarray(channel_radial_response(radius, z, 0), dtype=float)
        c1 = np.asarray(channel_radial_response(radius, z, 1), dtype=float)
        column_integrals = [float(np.trapezoid(c0, z)), float(np.trapezoid(c1, z))]
        combos = {}
        for beta0, beta1 in COEFFICIENT_WITNESSES:
            key = f"beta=({beta0:+.3g},{beta1:+.3g})"
            response = float(beta0) * c0 + float(beta1) * c1
            combos[key] = float(np.trapezoid(response, z))
        rows[f"r={radius:.2f}"] = {
            "core_radial_profile": float(_radial_profile(np.asarray([radius]), 0)[0]),
            "annular_radial_profile": float(_radial_profile(np.asarray([radius]), 1)[0]),
            "column_signed_integrals": column_integrals,
            "combination_signed_integrals": combos,
            "max_abs_signed_integral": float(
                max(abs(v) for v in [*column_integrals, *combos.values()])
            ),
        }
    return rows


def _response_geometry() -> dict[str, Any]:
    radii = np.linspace(0.1, 1.55, RANK_RADIAL_N)
    z = np.linspace(1.0001, 1.7999, RANK_AXIAL_N)
    rr, zz = np.meshgrid(radii, z, indexing="ij")
    c0 = np.asarray(channel_radial_response(rr, zz, 0), dtype=float).ravel()
    c1 = np.asarray(channel_radial_response(rr, zz, 1), dtype=float).ravel()
    n0, n1 = float(np.linalg.norm(c0)), float(np.linalg.norm(c1))
    if n0 <= np.finfo(float).tiny or n1 <= np.finfo(float).tiny:
        raise RuntimeError("degenerate frozen response column")
    matrix = np.column_stack((c0 / n0, c1 / n1))
    singular = np.linalg.svd(matrix, compute_uv=False)
    rank = int(np.linalg.matrix_rank(matrix, tol=1.0e-10))
    cosine = float(np.dot(matrix[:, 0], matrix[:, 1]))
    return {
        "rank": rank,
        "column_cosine": cosine,
        "condition_number": float(singular[0] / singular[-1]),
        "singular_values": singular.tolist(),
        "core_response_norm": n0,
        "annular_response_norm": n1,
        "cloud": {
            "radial_nodes": RANK_RADIAL_N,
            "axial_nodes": RANK_AXIAL_N,
            "r_interval": [0.1, 1.55],
            "z_interval": [1.0001, 1.7999],
        },
    }


def evaluate() -> dict[str, Any]:
    _assert_source_lock()

    endpoint_checks = {}
    for radius in RADII:
        per_channel = {}
        for channel in (0, 1):
            inner = endpoint_potential(radius, Z_INNER, channel)
            outer = endpoint_potential(radius, Z_OUTER, channel)
            exact_integral_from_endpoints = float(-radius * (outer - inner))
            if inner != 0.0 or outer != 0.0 or exact_integral_from_endpoints != 0.0:
                raise RuntimeError("compact endpoint identity drifted")
            per_channel[f"channel_{channel}"] = {
                "F_inner": inner,
                "F_outer": outer,
                "exact_integral_from_endpoints": exact_integral_from_endpoints,
            }
        endpoint_checks[f"r={radius:.2f}"] = per_channel

    levels = {str(n): _fixed_radius_integrals(n) for n in GRID_SIZES}
    medium = levels[str(GRID_SIZES[-2])]
    fine = levels[str(GRID_SIZES[-1])]
    max_fine = max(row["max_abs_signed_integral"] for row in fine.values())
    max_refinement_change = 0.0
    for key in fine:
        frow = fine[key]
        mrow = medium[key]
        max_refinement_change = max(
            max_refinement_change,
            abs(frow["max_abs_signed_integral"] - mrow["max_abs_signed_integral"]),
        )
    numerical_identity_pass = bool(
        max_fine <= SIGNED_INTEGRAL_TOL
        and max_refinement_change <= REFINEMENT_CHANGE_TOL
    )
    if not numerical_identity_pass:
        raise RuntimeError("fixed-radius signed-balance refinement gate failed")

    geometry = _response_geometry()
    independent = bool(
        geometry["rank"] == 2
        and abs(float(geometry["column_cosine"])) < RANK_COSINE_MAX_ABS
    )
    if not independent:
        raise RuntimeError("frozen two-channel witness is not independently resolved")

    inner_annular_zero = all(
        fine[f"r={r:.2f}"]["annular_radial_profile"] == 0.0 for r in RADII[:2]
    )
    outer_annular_active = all(
        fine[f"r={r:.2f}"]["annular_radial_profile"] > 0.0 for r in RADII[2:]
    )
    if not (inner_annular_zero and outer_annular_active):
        raise RuntimeError("annular radial-localization witness drifted")

    return {
        "task_id": TASK_ID,
        "issue": ISSUE,
        "source_parent": {"pr": SOURCE_PARENT_PR, "head": SOURCE_PARENT_HEAD},
        "source_context": {
            "reshape_pr": SOURCE_RESHAPE_PR,
            "obstruction_pr": SOURCE_OBSTRUCTION_PR,
            "nonlinear_pr": RELATED_NONLINEAR_PR,
            "fresh_holdout_pr": RELATED_FRESH_HOLDOUT_PR,
        },
        "frozen_family": {
            "vector_potential": "A=(-y*F,x*F,0)",
            "radial_component": "C_r=-r*partial_z(F)",
            "positive_tip_lobe": [Z_INNER, Z_OUTER],
            "core_radial_q_support": [CORE_RADIAL_A, CORE_RADIAL_B],
            "annular_radial_q_support": [ANNULAR_RADIAL_A, ANNULAR_RADIAL_B],
            "core_axial_exponents": [CORE_AXIAL_P, AXIAL_M],
            "annular_axial_exponents": [ANNULAR_AXIAL_P, AXIAL_M],
            "radii": list(RADII),
            "grid_sizes": list(GRID_SIZES),
            "coefficient_witnesses": [list(v) for v in COEFFICIENT_WITNESSES],
        },
        "exact_fixed_radius_identity": {
            "formula": "integral C_r dz = -r*(F(z_out)-F(z_in)) = 0",
            "endpoint_checks": endpoint_checks,
            "applies_to_arbitrary_finite_linear_combination": True,
        },
        "numerical_refinement": levels,
        "numerical_gate": {
            "max_abs_finest_signed_integral": float(max_fine),
            "max_medium_to_fine_change": float(max_refinement_change),
            "signed_integral_tolerance": SIGNED_INTEGRAL_TOL,
            "refinement_change_tolerance": REFINEMENT_CHANGE_TOL,
            "passes": numerical_identity_pass,
        },
        "response_geometry": geometry,
        "radial_localization_witness": {
            "annular_zero_at_inner_probe_radii": inner_annular_zero,
            "annular_active_at_outer_probe_radii": outer_annular_active,
        },
        "decision": {
            "independent_two_channel_witness": independent,
            "fixed_radius_signed_balance_survives_basis_growth": True,
            "radial_only_second_channel_can_move_return_debt_between_radii": False,
            "radial_localization_only_is_structural_escape": False,
            "add_second_radial_localized_channel_now": False,
            "next_if_740_fails": (
                "test one axial return-reservoir extension that keeps the visible 1<|z|<1.8 "
                "tip region more inward while placing compensating correction nearer the global "
                "support collar 1.8<|z|<2.0; preregister it before any child replay"
            ),
            "closer_visualization_delivery_this_increment": False,
        },
        **TRUTH,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("outputs/agent7-st052m-fixed-radius-return-balance-audit/report.json"),
    )
    args = parser.parse_args()
    report = evaluate()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
