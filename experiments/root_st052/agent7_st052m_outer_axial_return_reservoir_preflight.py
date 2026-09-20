"""Preflight one outer axial return reservoir for the ST052-M tip channel.

Preregistered in issue #766 and stacked on exact Agent-7 PR #759 head.
This is an expression-capacity diagnostic only: it changes no candidate velocity.

PR #759 establishes that radial localization alone cannot escape the fixed-radius
signed-balance identity while all compact axial channels share endpoints 1 and
1.8.  The global support, however, already extends to |z|<2.  This preflight
changes only the axial potential semantics and asks whether the mandatory return
flow can be placed in 1.8<|z|<2 while the visible 1<|z|<1.8 tip lobe has an
inward radial correction.

For A=(-y*f,x*f,0), f=R(r^2) Z(z), C=curl(A),

    C_r = -r R(r^2) Z'(z).

Z is an odd C3 piecewise-septic potential. On each axial side its magnitude
rises smoothly from 0 to 1 over 1<|z|<1.8 and falls smoothly back to 0 over
1.8<|z|<2. Thus Z' is even, positive in the visible lobe and negative in the
outer reservoir. The full [1,2] signed radial integral remains exactly zero;
the representation relocates, rather than removes, the required return flow.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np

import agent7_st052m_fixed_radius_return_balance_audit as fixed
import agent7_st052m_tip_odd_poloidal_screen as screen
import agent7_st052m_tip_return_flow_relocation as relocation
import agent7_st052m_vector_potential_poloidal_screen as vp

TASK_ID = "CR003-ST052M-OUTER-AXIAL-RETURN-RESERVOIR-PREFLIGHT-108"
ISSUE = 766
SOURCE_PARENT_PR = 759
SOURCE_PARENT_HEAD = "40c7c1cd0ed611788b6766a0e5b3a5782a0418d7"
SOURCE_P9_PR = 740
SOURCE_RESHAPE_PR = 732
RELATED_FRESH_HOLDOUT_PR = 714

Z_INNER = 1.0
Z_VISIBLE_OUTER = 1.8
Z_GLOBAL_OUTER = 2.0
VISIBLE_WIDTH = Z_VISIBLE_OUTER - Z_INNER
RESERVOIR_WIDTH = Z_GLOBAL_OUTER - Z_VISIBLE_OUTER
RADIAL_Q_A = -1.0
RADIAL_Q_B = 2.56
R_SUPPORT_MAX = 1.6

VISIBLE_PROBES = (1.15, 1.55, 1.75)
RESERVOIR_PROBES = (1.85, 1.95)
SHOULDER_PROBE = 0.85
RADII = (0.55, 0.95, 1.30, 1.50)
GRID_SIZES = (5001, 20001, 80001)
RESPONSE_RADIAL_N = 48
RESPONSE_AXIAL_N = 321
INTEGRAL_TOL = 1.0e-9
REFINEMENT_TOL = 1.0e-9
COSINE_MAX_ABS = 0.999
CONDITION_MAX = 25.0

TRUTH = {
    "canonical_velocity_changed": False,
    "saved_velocity_changed": False,
    "candidate_velocity_changed": False,
    "basis_dimension_changed": False,
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
    fixed._assert_source_lock()
    expected = {
        "R2_BUMP_A": RADIAL_Q_A,
        "R2_BUMP_B": RADIAL_Q_B,
        "R_SUPPORT_MAX": R_SUPPORT_MAX,
        "Z_SUPPORT_INNER": Z_INNER,
        "Z_SUPPORT_OUTER": Z_VISIBLE_OUTER,
    }
    for name, value in expected.items():
        actual = float(getattr(screen, name))
        if actual != float(value):
            raise RuntimeError(f"frozen source constant drifted: {name}={actual!r}")
    if relocation.select_inner_exponent()[0] != 9:
        raise RuntimeError("frozen p=9 reshape identity drifted")


def septic_smoothstep(x: np.ndarray) -> np.ndarray:
    """C3 endpoint-flat step S:[0,1]->[0,1]."""
    x = np.asarray(x, dtype=float)
    return 35.0 * x**4 - 84.0 * x**5 + 70.0 * x**6 - 20.0 * x**7


def septic_smoothstep_prime(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    return 140.0 * x**3 * (1.0 - x) ** 3


def axial_potential_and_derivative(z: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return odd Z(z) and even dZ/dz for the frozen reservoir construction."""
    z = np.asarray(z, dtype=float)
    s = np.abs(z)
    h = np.zeros_like(s)
    hp = np.zeros_like(s)

    visible = (s > Z_INNER) & (s < Z_VISIBLE_OUTER)
    if np.any(visible):
        x = (s[visible] - Z_INNER) / VISIBLE_WIDTH
        h[visible] = septic_smoothstep(x)
        hp[visible] = septic_smoothstep_prime(x) / VISIBLE_WIDTH

    reservoir = (s > Z_VISIBLE_OUTER) & (s < Z_GLOBAL_OUTER)
    if np.any(reservoir):
        x = (s[reservoir] - Z_VISIBLE_OUTER) / RESERVOIR_WIDTH
        h[reservoir] = 1.0 - septic_smoothstep(x)
        hp[reservoir] = -septic_smoothstep_prime(x) / RESERVOIR_WIDTH

    join = s == Z_VISIBLE_OUTER
    h[join] = 1.0
    hp[join] = 0.0

    potential = np.sign(z) * h
    derivative = hp
    return potential, derivative


def _radial_profile_and_derivative(r2: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    return vp._poly_bump(np.asarray(r2, dtype=float), RADIAL_Q_A, RADIAL_Q_B)


def reservoir_correction(points: np.ndarray) -> np.ndarray:
    """Analytic curl(A) for the frozen one-channel reservoir preflight."""
    pts = np.asarray(points, dtype=float)
    if pts.ndim != 2 or pts.shape[1] != 3:
        raise ValueError("points must have shape (n,3)")
    x, y, z = pts[:, 0], pts[:, 1], pts[:, 2]
    r2 = x * x + y * y
    radial, radial_s = _radial_profile_and_derivative(r2)
    axial, axial_z = axial_potential_and_derivative(z)
    f = radial * axial
    f_s = radial_s * axial
    f_z = radial * axial_z
    correction = np.column_stack((-x * f_z, -y * f_z, 2.0 * f + 2.0 * r2 * f_s))
    if not np.isfinite(correction).all():
        raise RuntimeError("nonfinite reservoir correction")
    return correction


def radial_response(radius: float | np.ndarray, z: np.ndarray) -> np.ndarray:
    rr = np.asarray(radius, dtype=float)
    zz = np.asarray(z, dtype=float)
    radial, _ = _radial_profile_and_derivative(rr * rr)
    _, axial_z = axial_potential_and_derivative(zz)
    return -rr * radial * axial_z


def _trapz(values: np.ndarray, grid: np.ndarray) -> float:
    if hasattr(np, "trapezoid"):
        return float(np.trapezoid(values, grid))
    return float(np.trapz(values, grid))


def _probe_radial(correction_fn: Any, radius: float, z: float) -> float:
    point = np.asarray([[float(radius), 0.0, float(z)]], dtype=float)
    return float(np.asarray(correction_fn(point), dtype=float)[0, 0])


def _directional_spot_checks() -> dict[str, Any]:
    radius = 1.0
    visible: dict[str, float] = {}
    reservoir: dict[str, float] = {}
    symmetry_error = 0.0
    for abs_z in VISIBLE_PROBES:
        plus = _probe_radial(reservoir_correction, radius, abs_z)
        minus = _probe_radial(reservoir_correction, radius, -abs_z)
        if not (plus < 0.0 and minus < 0.0):
            raise RuntimeError(f"reservoir channel not inward at frozen visible |z|={abs_z}")
        symmetry_error = max(symmetry_error, abs(plus - minus))
        visible[f"|z|={abs_z:.2f}"] = plus
    for abs_z in RESERVOIR_PROBES:
        plus = _probe_radial(reservoir_correction, radius, abs_z)
        minus = _probe_radial(reservoir_correction, radius, -abs_z)
        if not (plus > 0.0 and minus > 0.0):
            raise RuntimeError(f"reservoir channel lacks outward return at |z|={abs_z}")
        symmetry_error = max(symmetry_error, abs(plus - minus))
        reservoir[f"|z|={abs_z:.2f}"] = plus

    p9_at_175 = _probe_radial(relocation.reshaped_tip_correction, radius, 1.75)
    reservoir_at_175 = _probe_radial(reservoir_correction, radius, 1.75)
    if not (p9_at_175 > 0.0 and reservoir_at_175 < 0.0):
        raise RuntimeError("primary p=9-vs-reservoir |z|=1.75 directional witness failed")

    return {
        "visible_inward_radial_components": visible,
        "reservoir_outward_radial_components": reservoir,
        "positive_negative_z_max_radial_mismatch": float(symmetry_error),
        "p9_radial_at_abs_z_1p75": p9_at_175,
        "reservoir_radial_at_abs_z_1p75": reservoir_at_175,
        "p9_outward_reservoir_inward_at_1p75": True,
    }


def _support_checks() -> dict[str, float]:
    shoulder = np.asarray([[1.0, 0.0, SHOULDER_PROBE], [1.0, 0.0, -SHOULDER_PROBE]], dtype=float)
    axial_outside = np.asarray(
        [[1.0, 0.0, 2.0], [1.0, 0.0, -2.0], [1.0, 0.0, 2.05], [1.0, 0.0, -2.05]],
        dtype=float,
    )
    radial_outside = np.asarray([[1.61, 0.0, 1.55], [-1.61, 0.0, -1.55]], dtype=float)
    shoulder_max = float(np.max(np.abs(reservoir_correction(shoulder))))
    axial_outside_max = float(np.max(np.abs(reservoir_correction(axial_outside))))
    radial_outside_max = float(np.max(np.abs(reservoir_correction(radial_outside))))
    if max(shoulder_max, axial_outside_max, radial_outside_max) != 0.0:
        raise RuntimeError("reservoir correction leaked outside frozen support semantics")
    return {
        "shoulder_max_abs": shoulder_max,
        "axial_outside_max_abs": axial_outside_max,
        "radial_outside_max_abs": radial_outside_max,
    }


def _integral_level(n: int) -> dict[str, Any]:
    z = np.linspace(Z_INNER, Z_GLOBAL_OUTER, int(n))
    split = int(round((Z_VISIBLE_OUTER - Z_INNER) * (int(n) - 1)))
    if not np.isclose(z[split], Z_VISIBLE_OUTER, rtol=0.0, atol=1.0e-15):
        raise RuntimeError("frozen refinement grid no longer resolves z=1.8")
    rows: dict[str, Any] = {}
    for radius in RADII:
        response = np.asarray(radial_response(radius, z), dtype=float)
        full = _trapz(response, z)
        visible = _trapz(response[: split + 1], z[: split + 1])
        reservoir = _trapz(response[split:], z[split:])
        radial, _ = _radial_profile_and_derivative(np.asarray([radius * radius]))
        exact_visible = -float(radius) * float(radial[0])
        exact_reservoir = -exact_visible
        if not (visible < 0.0 and reservoir > 0.0):
            raise RuntimeError("visible/reservoir signed-integral directions failed")
        rows[f"r={radius:.2f}"] = {
            "full_signed_integral": full,
            "visible_signed_integral": visible,
            "reservoir_signed_integral": reservoir,
            "visible_plus_reservoir": visible + reservoir,
            "exact_visible_from_potential_endpoints": exact_visible,
            "exact_reservoir_from_potential_endpoints": exact_reservoir,
            "visible_endpoint_error": visible - exact_visible,
            "reservoir_endpoint_error": reservoir - exact_reservoir,
        }
    return rows


def _integral_refinement() -> dict[str, Any]:
    levels = {str(n): _integral_level(n) for n in GRID_SIZES}
    medium = levels[str(GRID_SIZES[-2])]
    fine = levels[str(GRID_SIZES[-1])]
    max_full = max(abs(row["full_signed_integral"]) for row in fine.values())
    max_balance = max(abs(row["visible_plus_reservoir"]) for row in fine.values())
    max_endpoint_error = max(
        max(abs(row["visible_endpoint_error"]), abs(row["reservoir_endpoint_error"]))
        for row in fine.values()
    )
    max_refinement_change = 0.0
    for key in fine:
        for field in ("full_signed_integral", "visible_signed_integral", "reservoir_signed_integral"):
            max_refinement_change = max(
                max_refinement_change,
                abs(float(fine[key][field]) - float(medium[key][field])),
            )
    passes = bool(
        max_full <= INTEGRAL_TOL
        and max_balance <= INTEGRAL_TOL
        and max_endpoint_error <= INTEGRAL_TOL
        and max_refinement_change <= REFINEMENT_TOL
    )
    if not passes:
        raise RuntimeError("outer-reservoir signed-integral refinement gate failed")
    return {
        "levels": levels,
        "gate": {
            "max_abs_finest_full_integral": float(max_full),
            "max_abs_finest_visible_plus_reservoir": float(max_balance),
            "max_abs_finest_endpoint_error": float(max_endpoint_error),
            "max_medium_to_fine_change": float(max_refinement_change),
            "integral_tolerance": INTEGRAL_TOL,
            "refinement_tolerance": REFINEMENT_TOL,
            "passes": passes,
        },
    }


def _response_geometry() -> dict[str, Any]:
    radii = np.linspace(0.1, 1.55, RESPONSE_RADIAL_N)
    z = np.linspace(1.0001, 1.9999, RESPONSE_AXIAL_N)
    rr, zz = np.meshgrid(radii, z, indexing="ij")
    points = np.column_stack((rr.ravel(), np.zeros(rr.size), zz.ravel()))
    new = np.asarray(reservoir_correction(points), dtype=float)[:, 0]
    old = np.asarray(relocation.reshaped_tip_correction(points), dtype=float)[:, 0]
    n_new = float(np.linalg.norm(new))
    n_old = float(np.linalg.norm(old))
    if n_new <= np.finfo(float).tiny or n_old <= np.finfo(float).tiny:
        raise RuntimeError("degenerate frozen response column")
    matrix = np.column_stack((old / n_old, new / n_new))
    singular = np.linalg.svd(matrix, compute_uv=False)
    rank = int(np.linalg.matrix_rank(matrix, tol=1.0e-10))
    cosine = float(np.dot(matrix[:, 0], matrix[:, 1]))
    condition = float(singular[0] / singular[-1])
    passes = bool(rank == 2 and abs(cosine) < COSINE_MAX_ABS and condition <= CONDITION_MAX)
    if not passes:
        raise RuntimeError("outer-reservoir response is not independently resolved")
    return {
        "rank": rank,
        "column_cosine": cosine,
        "condition_number": condition,
        "singular_values": singular.tolist(),
        "p9_response_norm": n_old,
        "reservoir_response_norm": n_new,
        "cloud": {
            "radial_nodes": RESPONSE_RADIAL_N,
            "axial_nodes": RESPONSE_AXIAL_N,
            "r_interval": [0.1, 1.55],
            "abs_z_interval": [1.0001, 1.9999],
        },
        "passes": passes,
    }


def evaluate() -> dict[str, Any]:
    _assert_source_lock()

    z_checks = np.asarray([1.0, 1.8, 2.0, -1.0, -1.8, -2.0], dtype=float)
    potential, derivative = axial_potential_and_derivative(z_checks)
    endpoint_join = {
        "z_values": z_checks.tolist(),
        "potential": potential.tolist(),
        "derivative": derivative.tolist(),
    }
    if not (
        potential[0] == 0.0
        and potential[1] == 1.0
        and potential[2] == 0.0
        and potential[3] == 0.0
        and potential[4] == -1.0
        and potential[5] == 0.0
        and np.max(np.abs(derivative)) == 0.0
    ):
        raise RuntimeError("frozen septic endpoint/join semantics drifted")

    directional = _directional_spot_checks()
    support = _support_checks()
    refinement = _integral_refinement()
    geometry = _response_geometry()

    decision_pass = bool(
        directional["p9_outward_reservoir_inward_at_1p75"]
        and directional["positive_negative_z_max_radial_mismatch"] == 0.0
        and refinement["gate"]["passes"]
        and geometry["passes"]
        and max(support.values()) == 0.0
    )
    if not decision_pass:
        raise RuntimeError("preregistered outer axial reservoir preflight failed")

    return {
        "task_id": TASK_ID,
        "issue": ISSUE,
        "source_parent": {"pr": SOURCE_PARENT_PR, "head": SOURCE_PARENT_HEAD},
        "source_context": {
            "nonlinear_p9_pr": SOURCE_P9_PR,
            "reshape_pr": SOURCE_RESHAPE_PR,
            "fresh_holdout_pr": RELATED_FRESH_HOLDOUT_PR,
        },
        "frozen_representation": {
            "vector_potential": "A=(-y*f,x*f,0)",
            "f": "R(r^2)*Z(z)",
            "radial_component": "C_r=-r*R(r^2)*Z_prime(z)",
            "smoothstep": "35*x^4-84*x^5+70*x^6-20*x^7",
            "visible_abs_z_interval": [Z_INNER, Z_VISIBLE_OUTER],
            "return_reservoir_abs_z_interval": [Z_VISIBLE_OUTER, Z_GLOBAL_OUTER],
            "radial_support": "r<1.6",
            "global_axial_support": "|z|<2",
        },
        "exact_endpoint_join_checks": endpoint_join,
        "directional_spot_checks": directional,
        "support_checks": support,
        "signed_integral_refinement": refinement,
        "response_geometry_vs_p9": geometry,
        "structural_scope": {
            "full_positive_lobe_integral_identity": "integral_1^2 C_r dz = -r*R*(Z(2)-Z(1)) = 0",
            "visible_integral_is_inward": True,
            "outer_reservoir_integral_is_outward": True,
            "return_flow_removed": False,
            "return_flow_relocated_outside_visible_1_to_1p8_lobe": True,
            "total_velocity_tip_inwardness_established": False,
        },
        "decision": {
            "classification": "outer_axial_return_reservoir_preflight_pass",
            "outer_axial_return_reservoir_is_structural_escape_from_1p8_endpoint_obstruction": True,
            "one_later_frozen_nonlinear_child_replay_justified": True,
            "add_second_basis_dimension_now": False,
            "closer_visualization_delivery_this_increment": False,
            "next_minimal_increment": (
                "only after this exact-head preflight is verified, freeze one nonlinear child using "
                "the same reservoir shape and a preregistered coefficient rule; then test total-field "
                "tip radial motion together with winding, axial-stretch, support/divergence and fresh-data boundaries"
            ),
        },
        **TRUTH,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("outputs/agent7-st052m-outer-axial-return-reservoir-preflight/report.json"),
    )
    args = parser.parse_args()
    report = evaluate()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
