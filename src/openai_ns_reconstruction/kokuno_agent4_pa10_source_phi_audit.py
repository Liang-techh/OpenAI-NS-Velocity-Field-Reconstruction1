"""Independent audit of the source-compatible PA.10 Phi-ball envelope.

This validator consumes the public scalar identity and public bound values from
Agent-1 PR #653, but it does not call the upstream private chi/inverse-series
helpers or use the upstream receipt as an oracle.  It reconstructs the source
H_* polynomial variation bound independently, stress-tests the complex chi
ratio on fresh off-grid complex points, evaluates the Cauchy coefficient sum
with high-precision infinite summation, and evaluates the positive inverse
series through the modified-Bessel closed form rather than Agent 1's rational
recurrence/tail-majorant implementation.

The result is only a prerequisite audit for the Phi component of the PA.10
radius-one contraction ball.  It is not a Navier--Stokes residual and cannot
set ``pde_validated``.
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

from .kokuno_pa10_source_phi_ball import KokunoPA10SourcePhiBallBounds

SCHEMA = "kokuno-agent4-pa10-source-phi-independent-audit-v1"
AGENT1_PR = 653
AGENT1_HEAD = "cc36aec054f1c3de744acb58929e9991ad823cbe"
SEED = 9173321
MP_DPS = 80
COMPLEX_SAMPLES = 20000
MUTATION_FACTOR = 0.999
INVERSE_MUTATION_FACTOR = 0.999999


def _canonical_json(payload: dict[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _mp(value: float) -> mp.mpf:
    return mp.mpf(repr(float(value)))


def _H(z: np.ndarray | complex, *, D: float, j0: float):
    return j0 + (D + 4.0) * z - j0 * z * z - 4.0 * z * z * z


def _independent_variation_bound(domain) -> mp.mpf:
    """Bound |H(z)-H(eta)| on the public tube via an independent derivative sup."""
    radius = (
        mp.mpf("1")
        + _mp(domain.enlarged_real_margin)
        + _mp(domain.complex_tube_radius)
    )
    derivative_sup = (
        abs(_mp(domain.D) + 4)
        + 2 * abs(_mp(domain.j0)) * radius
        + 12 * radius * radius
    )
    return _mp(domain.complex_tube_radius) * derivative_sup


def _factor_exact_max(e: mp.mpf) -> tuple[mp.mpf, mp.mpf]:
    if not 0 < e < 1:
        raise ValueError("relative H variation e must lie in (0,1)")
    # Independent calculus replay: solve r'(t)=0 for the unique t>=0 and
    # evaluate r there; do not consume the upstream claimed closed-form value.
    t_star = (1 - e * e) / (2 * e)
    ratio = (t_star + e) / (mp.sqrt(t_star * t_star + 1) - e)
    return t_star, ratio


def _infinite_weight_sum(x: mp.mpf) -> mp.mpf:
    return mp.nsum(lambda k: (k + 1) ** 2 * x**k, [0, mp.inf])


def _inverse_series_bessel(x: mp.mpf) -> mp.mpf:
    # sum_{k>=0} x^k/(k!(k+1)!) = I_1(2 sqrt(x))/sqrt(x)
    if x == 0:
        return mp.mpf(1)
    root = mp.sqrt(x)
    return mp.besseli(1, 2 * root) / root


def _relative_error(public: float, reference: mp.mpf) -> float:
    ref = float(reference)
    return abs(float(public) - ref) / max(abs(ref), 1.0e-300)


def run_audit() -> dict[str, Any]:
    mp.mp.dps = MP_DPS
    datum = KokunoPA10SourcePhiBallBounds()
    domain = datum.domain
    public_chi = datum.chi_certificate()
    public_phi = datum.phi_ball_certificate()

    sigma = _mp(domain.sigma_star)
    variation = _independent_variation_bound(domain)
    e = variation / sigma
    t_star, factor_max = _factor_exact_max(e)
    chi_sup = factor_max * factor_max

    x = _mp(domain.coefficient_rho) / _mp(domain.cauchy_radius)
    weight_sum = _infinite_weight_sum(x)
    multiplier = chi_sup * weight_sum
    inverse_x = 40 * multiplier
    inverse_exact = _inverse_series_bessel(inverse_x)
    phi0_exact = inverse_exact
    phi_ball_exact = phi0_exact + 1

    rng = np.random.default_rng(SEED)
    lo = -1.0 - float(domain.enlarged_real_margin)
    hi = 1.0 + float(domain.enlarged_real_margin)
    eta = rng.uniform(lo, hi, size=COMPLEX_SAMPLES)
    radii = float(domain.complex_tube_radius) * np.sqrt(rng.random(COMPLEX_SAMPLES))
    angles = rng.uniform(0.0, 2.0 * math.pi, size=COMPLEX_SAMPLES)
    dz = radii * np.exp(1j * angles)
    z = eta.astype(complex) + dz
    H_eta = _H(eta, D=float(domain.D), j0=float(domain.j0))
    H_z = _H(z, D=float(domain.D), j0=float(domain.j0))
    actual_variation = np.abs(H_z - H_eta)
    sigma_f = float(domain.sigma_star)
    chi_complex = np.abs(H_z * H_z / (H_z * H_z + sigma_f * sigma_f))
    direct_max_variation = float(np.max(actual_variation))
    direct_max_chi = float(np.max(chi_complex))

    public_factor = float(public_chi["single_factor_ratio_upper"])
    public_chi_sup = float(public_chi["chi_complex_sup_upper"])
    public_weight = float(public_chi["cauchy_weight_sum_upper"])
    public_multiplier = float(public_chi["chi_multiplier_norm_upper"])
    public_inverse = float(public_phi["inverse_one_plus_T_absolute_series_upper"])
    public_phi0 = float(public_phi["Phi0_coefficient_norm_upper"])
    public_phi_ball = float(public_phi["Phi_radius_one_ball_norm_upper"])
    public_lipschitz = float(public_phi["Phi_radius_one_ball_lipschitz_upper"])

    guards = {
        "independent_variation_is_subsigma": bool(0 < variation < sigma),
        "direct_complex_variation_below_independent_bound": bool(
            direct_max_variation <= float(variation) * (1.0 + 5.0e-13)
        ),
        "public_factor_dominates_independent_exact_max": bool(
            public_factor + 5.0e-14 >= float(factor_max)
        ),
        "public_chi_sup_dominates_independent_exact_max": bool(
            public_chi_sup + 5.0e-14 >= float(chi_sup)
        ),
        "direct_complex_chi_below_public_bound": bool(
            direct_max_chi <= public_chi_sup * (1.0 + 5.0e-13)
        ),
        "public_weight_sum_dominates_high_precision_nsum": bool(
            public_weight + 5.0e-14 >= float(weight_sum)
        ),
        "public_multiplier_dominates_independent_product": bool(
            public_multiplier + 5.0e-13 >= float(multiplier)
        ),
        "public_inverse_dominates_bessel_closed_form": bool(
            public_inverse + 5.0e-12 * max(1.0, public_inverse) >= float(inverse_exact)
        ),
        "public_phi0_dominates_independent_center": bool(
            public_phi0 + 5.0e-12 * max(1.0, public_phi0) >= float(phi0_exact)
        ),
        "public_radius_one_ball_dominates_center_plus_one": bool(
            public_phi_ball + 5.0e-12 * max(1.0, public_phi_ball) >= float(phi_ball_exact)
        ),
        "public_phi_projection_lipschitz_is_one": bool(public_lipschitz == 1.0),
    }

    mutation = {
        "factor_0p999_detected": bool(MUTATION_FACTOR * public_factor < float(factor_max)),
        "chi_sup_0p999_detected": bool(MUTATION_FACTOR * public_chi_sup < float(chi_sup)),
        "weight_sum_0p999_detected": bool(MUTATION_FACTOR * public_weight < float(weight_sum)),
        "multiplier_0p999_detected": bool(
            MUTATION_FACTOR * public_multiplier < float(multiplier)
        ),
        "inverse_0p999999_detected": bool(
            INVERSE_MUTATION_FACTOR * public_inverse < float(inverse_exact)
        ),
        "phi_ball_radius_drop_detected": bool(public_phi0 + 0.9 < float(phi_ball_exact)),
    }

    failed = [name for name, passed in guards.items() if not passed]
    if not all(mutation.values()):
        failed.append("mutation_detection")

    report: dict[str, Any] = {
        "schema": SCHEMA,
        "upstream": {
            "agent1_pr": AGENT1_PR,
            "agent1_head": AGENT1_HEAD,
            "source_repository": "KokunoYumeto/yang-mills-interacting-workbench",
            "source_commit": "143f6773feb424ad9ed3a8d116653200f20346b7",
            "source_path": "navier-stokes/navier_stokes_workbench.tex",
        },
        "frozen_protocol": {
            "seed": SEED,
            "mpmath_dps": MP_DPS,
            "complex_offgrid_samples": COMPLEX_SAMPLES,
            "independent_inverse_method": "mpmath besseli closed form",
            "independent_weight_method": "mpmath infinite nsum",
            "mutation_factor": MUTATION_FACTOR,
            "inverse_mutation_factor": INVERSE_MUTATION_FACTOR,
            "final_project_gates_unchanged": {
                "normalized_momentum_max": 1.0e-3,
                "normalized_momentum_L2": 1.0e-3,
                "divergence_max": 1.0e-5,
                "divergence_L2": 1.0e-5,
            },
        },
        "independent_chi": {
            "variation_bound": float(variation),
            "relative_variation_e": float(e),
            "stationary_t": float(t_star),
            "single_factor_exact_max": float(factor_max),
            "chi_complex_exact_sup": float(chi_sup),
            "public_single_factor_upper": public_factor,
            "public_chi_complex_sup_upper": public_chi_sup,
            "factor_relative_error": _relative_error(public_factor, factor_max),
            "chi_sup_relative_error": _relative_error(public_chi_sup, chi_sup),
            "direct_offgrid_max_H_variation": direct_max_variation,
            "direct_offgrid_max_abs_chi": direct_max_chi,
        },
        "independent_coefficient_multiplier": {
            "rho_over_cauchy": float(x),
            "high_precision_infinite_weight_sum": float(weight_sum),
            "public_weight_sum_upper": public_weight,
            "weight_relative_error": _relative_error(public_weight, weight_sum),
            "independent_chi_multiplier": float(multiplier),
            "public_chi_multiplier_upper": public_multiplier,
            "multiplier_relative_error": _relative_error(public_multiplier, multiplier),
        },
        "independent_inverse_and_ball": {
            "series_argument_40_Mchi": float(inverse_x),
            "bessel_inverse_exact": float(inverse_exact),
            "public_inverse_upper": public_inverse,
            "inverse_relative_overage": (public_inverse - float(inverse_exact))
            / max(float(inverse_exact), 1.0e-300),
            "independent_Phi0_norm": float(phi0_exact),
            "public_Phi0_upper": public_phi0,
            "independent_radius_one_ball_norm": float(phi_ball_exact),
            "public_radius_one_ball_upper": public_phi_ball,
            "public_Phi_lipschitz_upper": public_lipschitz,
        },
        "guards": guards,
        "mutation": mutation,
        "failed_guards": failed,
        "source_phi_ball_independent_preflight_passed": not failed,
        "truth_boundary": {
            "source_axis_domain_independent_A4_admission_required_separately": True,
            "source_Phi_ball_independently_audited": not failed,
            "source_Phi_radius_one_ball_norm_machine_bound_independently_audited": not failed,
            "source_Phi_radius_one_ball_lipschitz_machine_bound_independently_audited": not failed,
            "source_mixed_Y_Phi_Y_radius_one_ball_bound_machine_bound": False,
            "source_mixed_Y_u_Y_radius_one_ball_bound_machine_bound": False,
            "source_R1_R2_machine_bound": False,
            "source_operator_M_K_machine_bound": False,
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
    args = parser.parse_args()
    path = Path(args.output)
    path.parent.mkdir(parents=True, exist_ok=True)
    report = run_audit()
    path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
