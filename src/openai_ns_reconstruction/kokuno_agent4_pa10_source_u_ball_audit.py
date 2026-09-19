"""Independent Agent-4 audit of the source-compatible PA.10 u-ball.

This validator consumes only the public Agent-1 #662 u-center / u-ball API and
public source-compatible parameters.  It does not call Agent-1's private
Fraction/Cauchy helpers and does not use its receipt as an oracle.

The independent path reconstructs the displayed rational axis formulas,
rebuilds the complex-neighborhood Cauchy majorant at high precision with an
independent infinite-sum evaluation, and stress-tests the public center on fresh
real and complex off-grid points.  It audits only the ``u_0`` / radius-one
``u``-ball seam.  It is not a Navier--Stokes residual validation.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import mpmath as mp
import numpy as np

from .kokuno_pa10_source_u_ball import KokunoPA10SourceUBallBounds

SCHEMA = "kokuno-agent4-pa10-source-u-ball-independent-audit-v1"
AGENT1_PR = 662
AGENT1_HEAD = "5aa5eaa60644d5303fe5ab14c94881231b160075"
SOURCE_AXIS_A4_PR = 671
SOURCE_AXIS_A4_HEAD = "6d469766f7cf54730fcd0af3ccc597a2a300d4d4"
SOURCE_AXIS_A4_RUN = 35456226819
SEED = 9173351
MP_DPS = 90
REAL_SAMPLES = 4096
COMPLEX_SAMPLES = 20000
COMPLEX_SAMPLE_LADDER = (5000, 10000, 20000)
MUTATION_FACTOR = 0.999999


def _canonical_json(payload: dict[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _mp(value: float) -> mp.mpf:
    return mp.mpf(repr(float(value)))


def _independent_complex_majorant(domain) -> dict[str, mp.mpf]:
    """Rebuild a guaranteed Cauchy majorant without Agent-1 private helpers."""
    h = _mp(domain.h)
    A = mp.mpf("0.5") + h
    D = mp.mpf("0.5") - h
    j0 = _mp(domain.j0)
    pressure_square = _mp(domain.pressure_square)
    margin = _mp(domain.enlarged_real_margin)
    radius = _mp(domain.cauchy_radius)
    rho = _mp(domain.coefficient_rho)
    if not (0 < rho < radius < 1):
        raise ValueError("expected 0 < rho < cauchy_radius < 1")

    z_abs = 1 + margin + radius
    d_abs = 1 + z_abs * z_abs
    L_lower = 1 - 2 * h * z_abs * z_abs
    if L_lower <= 0:
        raise ValueError("independent Cauchy neighborhood does not keep L nonzero")

    U_abs = 4 * z_abs + abs(j0)
    one_minus_2zU_abs = 1 + 2 * z_abs * U_abs
    H_abs = abs(D) * z_abs + d_abs * U_abs

    # 1+z^2=(z-i)(z+i); each real center is distance >=1 from either pole.
    pole_factor_lower = 1 - radius
    one_plus_z2_lower = pole_factor_lower * pole_factor_lower
    Pi_abs = pressure_square / (one_plus_z2_lower**2)
    Pi_eta_abs = 4 * pressure_square * z_abs / (one_plus_z2_lower**3)

    Z_abs = (
        A * one_minus_2zU_abs * U_abs
        + 4 * H_abs
        + d_abs * Pi_eta_abs
        + 4 * A * z_abs * Pi_abs
    )
    slope_abs = Z_abs / (2 * L_lower)
    x = rho / radius
    # Independent high-precision infinite summation rather than the upstream
    # exact-Fraction closed form (1-x)^-2.
    alpha1_weight_sum = mp.nsum(lambda k: (k + 1) * x**k, [0, mp.inf])
    u0_norm = 80 * slope_abs * alpha1_weight_sum
    return {
        "z_abs_upper": z_abs,
        "d_abs_upper": d_abs,
        "L_abs_lower": L_lower,
        "U_star_abs_upper": U_abs,
        "H_star_abs_upper": H_abs,
        "Pi0_abs_upper": Pi_abs,
        "Pi0_eta_abs_upper": Pi_eta_abs,
        "Z_star_abs_upper": Z_abs,
        "u0_Y_abs_upper": slope_abs,
        "rho_over_cauchy_radius": x,
        "alpha1_weight_sum": alpha1_weight_sum,
        "u0_coefficient_norm_upper": u0_norm,
        "radius_one_u_ball_norm_upper": u0_norm + 1,
    }


def _independent_axis_values(domain, Y: np.ndarray, eta: np.ndarray) -> dict[str, np.ndarray]:
    """Direct public-formula replay independent of ``domain.axis_state``."""
    h = float(domain.h)
    A = 0.5 + h
    D = 0.5 - h
    j0 = float(domain.j0)
    P2 = float(domain.pressure_square)

    d = 1.0 - eta * eta
    L = 1.0 - 2.0 * h * eta * eta
    U = 4.0 * eta + j0
    H = D * eta + d * U
    Pi = -P2 / (1.0 + eta * eta) ** 2
    Pi_eta = 4.0 * P2 * eta / (1.0 + eta * eta) ** 3
    Z = -A * (1.0 - 2.0 * eta * U) * U - 4.0 * H - d * Pi_eta + 4.0 * A * eta * Pi
    slope = -Z / (2.0 * L)
    return {
        "u_0": Y * slope,
        "u_0_Y": slope,
        "Y_u_0_Y": Y * slope,
        "Z_star": Z,
        "L": L,
    }


def _complex_slope_values(domain, z: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    h = float(domain.h)
    A = 0.5 + h
    D = 0.5 - h
    j0 = float(domain.j0)
    P2 = float(domain.pressure_square)
    d = 1.0 - z * z
    L = 1.0 - 2.0 * h * z * z
    U = 4.0 * z + j0
    H = D * z + d * U
    Pi = -P2 / (1.0 + z * z) ** 2
    Pi_eta = 4.0 * P2 * z / (1.0 + z * z) ** 3
    Z = -A * (1.0 - 2.0 * z * U) * U - 4.0 * H - d * Pi_eta + 4.0 * A * z * Pi
    return -Z / (2.0 * L), L


def _relative_rms(error: np.ndarray, reference: np.ndarray) -> float:
    num = float(np.sqrt(np.mean(np.square(np.asarray(error, dtype=float)))))
    den = float(np.sqrt(np.mean(np.square(np.asarray(reference, dtype=float)))))
    return num / max(den, 1.0e-300)


def run_audit(*, pr_head: str = "unbound", checkout_head: str = "unbound") -> dict[str, Any]:
    mp.mp.dps = MP_DPS
    datum = KokunoPA10SourceUBallBounds()
    domain = datum.domain
    public_cert = datum.u0_certificate()
    public_ball = datum.u_ball()
    independent = _independent_complex_majorant(domain)

    rng = np.random.default_rng(SEED)
    E = 1.0 + float(domain.enlarged_real_margin)
    eta = rng.uniform(-E, E, size=REAL_SAMPLES)
    Y = rng.uniform(0.0, 4.1, size=REAL_SAMPLES)
    # Force an explicit near-axis/center subset without changing the frozen seed.
    Y[:8] = np.asarray([0.0, 1e-14, 1e-12, 1e-10, 1e-8, 1e-6, 1e-4, 1e-2])
    eta[:8] = np.asarray([-1.0, -0.5, -1e-12, 0.0, 1e-12, 0.5, 0.999999999999, 1.0])

    public_values = datum.center_values(Y, eta)
    reference_values = _independent_axis_values(domain, Y, eta)
    value_errors: dict[str, dict[str, float]] = {}
    for name in ("u_0", "u_0_Y", "Y_u_0_Y", "Z_star", "L"):
        err = np.asarray(public_values[name], dtype=float) - reference_values[name]
        ref = reference_values[name]
        value_errors[name] = {
            "max_abs": float(np.max(np.abs(err))),
            "relative_rms": _relative_rms(err, ref),
        }

    # Fresh complex off-grid stress cloud over the full certified Cauchy union.
    centers = rng.uniform(-E, E, size=COMPLEX_SAMPLES)
    radii = float(domain.cauchy_radius) * np.sqrt(rng.random(COMPLEX_SAMPLES))
    angles = rng.uniform(0.0, 2.0 * math.pi, size=COMPLEX_SAMPLES)
    z = centers.astype(complex) + radii * np.exp(1j * angles)
    slope_complex, L_complex = _complex_slope_values(domain, z)
    slope_abs = np.abs(slope_complex)
    ladder = {
        str(n): float(np.max(slope_abs[:n])) for n in COMPLEX_SAMPLE_LADDER
    }
    complex_max_slope = float(np.max(slope_abs))
    complex_min_abs_L = float(np.min(np.abs(L_complex)))

    public_slope = float(public_cert["cauchy_u0_Y_abs_upper"])
    public_weight = float(public_cert["alpha1_cauchy_weight_sum_upper"])
    public_u0_norm = float(public_cert["u0_coefficient_norm_upper"])
    public_ball_norm = float(public_ball.norm)
    public_ball_lip = float(public_ball.lipschitz)

    ind_slope = float(independent["u0_Y_abs_upper"])
    ind_weight = float(independent["alpha1_weight_sum"])
    ind_u0_norm = float(independent["u0_coefficient_norm_upper"])
    ind_ball_norm = float(independent["radius_one_u_ball_norm_upper"])

    guards = {
        "public_center_u0_matches_independent_formula": bool(
            value_errors["u_0"]["relative_rms"] <= 2e-14
            and value_errors["u_0"]["max_abs"] <= 2e-10 * max(1.0, float(np.max(np.abs(reference_values["u_0"]))))
        ),
        "public_center_u0Y_matches_independent_formula": bool(
            value_errors["u_0_Y"]["relative_rms"] <= 2e-14
        ),
        "public_center_Z_L_match_independent_formula": bool(
            value_errors["Z_star"]["relative_rms"] <= 2e-14
            and value_errors["L"]["relative_rms"] <= 2e-14
        ),
        "axis_origin_u0_is_exact_zero": bool(float(np.asarray(public_values["u_0"])[0]) == 0.0),
        "center_is_nontrivial": bool(float(np.sqrt(np.mean(reference_values["u_0_Y"] ** 2))) > 0.0),
        "complex_L_stays_nonzero_on_fresh_cloud": bool(complex_min_abs_L > 0.0),
        "complex_offgrid_slope_below_independent_majorant": bool(
            complex_max_slope <= ind_slope * (1.0 + 2.0e-12)
        ),
        "public_slope_upper_dominates_independent_majorant": bool(
            public_slope + 5.0e-13 * max(1.0, public_slope) >= ind_slope
        ),
        "public_weight_sum_dominates_independent_nsum": bool(
            public_weight + 5.0e-14 >= ind_weight
        ),
        "public_u0_norm_dominates_independent_bound": bool(
            public_u0_norm + 5.0e-13 * max(1.0, public_u0_norm) >= ind_u0_norm
        ),
        "public_radius_one_ball_dominates_independent_center_plus_one": bool(
            public_ball_norm + 5.0e-13 * max(1.0, public_ball_norm) >= ind_ball_norm
        ),
        "public_u_projection_lipschitz_is_one": bool(public_ball_lip == 1.0),
    }

    mutation = {
        "slope_0p999999_detected": bool(MUTATION_FACTOR * public_slope < ind_slope),
        "weight_0p999999_detected": bool(MUTATION_FACTOR * public_weight < ind_weight),
        "u0_norm_0p999999_detected": bool(MUTATION_FACTOR * public_u0_norm < ind_u0_norm),
        "radius_increment_drop_to_0p9_detected": bool(public_u0_norm + 0.9 < ind_ball_norm),
        "center_sign_flip_detected": bool(
            _relative_rms(-np.asarray(public_values["u_0_Y"]) - reference_values["u_0_Y"], reference_values["u_0_Y"]) > 1.9
        ),
    }

    failed = [name for name, passed in guards.items() if not passed]
    if not all(mutation.values()):
        failed.append("mutation_detection")

    report: dict[str, Any] = {
        "schema": SCHEMA,
        "pr_head": pr_head,
        "checkout_head": checkout_head,
        "upstream": {
            "agent1_pr": AGENT1_PR,
            "agent1_head": AGENT1_HEAD,
            "source_axis_agent4_pr": SOURCE_AXIS_A4_PR,
            "source_axis_agent4_head": SOURCE_AXIS_A4_HEAD,
            "source_axis_agent4_dedicated_run": SOURCE_AXIS_A4_RUN,
            "source_repository": "KokunoYumeto/yang-mills-interacting-workbench",
            "source_commit": "143f6773feb424ad9ed3a8d116653200f20346b7",
            "source_path": "navier-stokes/navier_stokes_workbench.tex",
        },
        "frozen_protocol": {
            "seed": SEED,
            "mpmath_dps": MP_DPS,
            "real_offgrid_samples": REAL_SAMPLES,
            "complex_offgrid_samples": COMPLEX_SAMPLES,
            "complex_sample_ladder": list(COMPLEX_SAMPLE_LADDER),
            "independent_bound_path": "90-digit mpmath + infinite nsum + direct rational formulas",
            "mutation_factor": MUTATION_FACTOR,
            "final_project_gates_unchanged": {
                "normalized_momentum_max": 1.0e-3,
                "normalized_momentum_L2": 1.0e-3,
                "divergence_max": 1.0e-5,
                "divergence_L2": 1.0e-5,
            },
        },
        "independent_complex_majorant": {name: float(value) for name, value in independent.items()},
        "public_certificate": {
            "cauchy_u0_Y_abs_upper": public_slope,
            "alpha1_cauchy_weight_sum_upper": public_weight,
            "u0_coefficient_norm_upper": public_u0_norm,
            "radius_one_u_ball_norm": public_ball_norm,
            "radius_one_u_ball_lipschitz": public_ball_lip,
        },
        "public_to_independent_ratios": {
            "slope": public_slope / ind_slope,
            "weight_sum": public_weight / ind_weight,
            "u0_norm": public_u0_norm / ind_u0_norm,
            "radius_one_ball_norm": public_ball_norm / ind_ball_norm,
        },
        "fresh_real_value_path": {
            "errors": value_errors,
            "u0_Y_rms": float(np.sqrt(np.mean(reference_values["u_0_Y"] ** 2))),
            "max_abs_u0": float(np.max(np.abs(reference_values["u_0"]))),
        },
        "fresh_complex_offgrid_stress": {
            "max_abs_u0_Y_by_nested_sample_count": ladder,
            "max_abs_u0_Y": complex_max_slope,
            "min_abs_L": complex_min_abs_L,
            "max_to_independent_majorant_ratio": complex_max_slope / ind_slope,
        },
        "guards": guards,
        "mutation": mutation,
        "failed_guards": failed,
        "source_u_ball_independent_preflight_passed": not failed,
        "truth_boundary": {
            "source_axis_domain_independently_audited_by_A4_671": True,
            "source_u0_formula_independently_audited": not failed,
            "source_u0_coefficient_norm_independently_audited": not failed,
            "source_u_radius_one_ball_norm_independently_audited": not failed,
            "source_u_radius_one_ball_lipschitz_independently_audited": not failed,
            "agent1_684_algebraic_R2_slots_independently_audited": False,
            "all_R2_ordinary_slots_source_bound": False,
            "all_R1_ordinary_slots_source_bound": False,
            "source_operator_constant_M_machine_bound": False,
            "source_operator_constant_K_machine_bound": False,
            "global_pressure_matched": False,
            "global_leading_profile_reconstructed": False,
            "leading_only_ns_residual_assessed": False,
            "leading_plus_oscillatory_ns_residual_assessed": False,
            "after_correction_ns_residual_assessed": False,
            "heldout_ns_momentum_residual_assessed": False,
            "formal_full_domain_pde_gate_assessed": False,
            "pde_validated": False,
        },
    }
    identity = dict(report)
    report["receipt_sha256"] = hashlib.sha256(
        _canonical_json(identity).encode("utf-8")
    ).hexdigest()
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    parser.add_argument("--pr-head", default="unbound")
    parser.add_argument("--checkout-head", default="unbound")
    args = parser.parse_args()
    path = Path(args.output)
    path.parent.mkdir(parents=True, exist_ok=True)
    report = run_audit(pr_head=args.pr_head, checkout_head=args.checkout_head)
    path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
