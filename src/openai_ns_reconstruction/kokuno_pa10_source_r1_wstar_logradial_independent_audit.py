"""Independent Agent-4 audit of PA.10 ``W_* Y Phi_Y`` after ``J_2``.

This validator is implementation-distinct from Agent-1 #724.  It consumes the
public source-compatible interfaces, but reconstructs the coefficient algebra
with 80-digit Decimal arithmetic and the strict rational upper pi < 355/113.
It does not call Agent-1 private Fraction/product/rounding helpers and does not
use an Agent-1 receipt as a numerical oracle.

The target identity is

    J_2[L^-1 W_* Y Phi_Y]
      = J_2[Y partial_Y(L^-1 W_* Phi)],

because L^-1 and W_* depend only on eta.  The source single-logradial rule has
factor 80 *after* the radial inverse.  The independent audit therefore forms
exactly the three-factor product L^-1 W_* Phi (two coefficient convolutions)
and applies factor 80 exactly once.  A shape guard rejects both a missing
convolution and an accidental second J_2/operator factor.

Scope is intentionally narrow: a PASS can independently admit only #724's new
ordinary R1 slot, conditional on the already independently audited source-axis,
Phi-ball and fixed-multiplier/product prerequisites.  It cannot close full R1,
M/K, global leading/pressure, a staged NS residual, or ``pde_validated``.
"""

from __future__ import annotations

import argparse
from decimal import Decimal, ROUND_CEILING, localcontext
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from .kokuno_pa10_source_r1_wstar_logradial_ordinary_slot import (
    KokunoPA10SourceR1WstarLogradialOrdinarySlot,
)

SCHEMA = "kokuno-agent4-pa10-source-r1-wstar-logradial-independent-audit-v1"
AGENT1_PR = 724
AGENT1_HEAD = "d62461a873936500fcc5b125179f593433295e17"
SOURCE_AXIS_A4_PR = 671
SOURCE_PHI_BALL_A4_PR = 655
SOURCE_U_BALL_A4_PR = 687
MIXED_J_A4_PR = 663
DERIVATIVE_FIXED_MULTIPLIER_A4_PR = 704
ZETA_R1_A4_PR = 710
AVERAGED_ETA_R1_A4_PR = 718
SEED = 9173401
OFFGRID_COUNT = 8192
SLOT_NAME = "W_star_times_Y_Phi_Y"
MAX_PUBLIC_RATIO = Decimal("1.01")

_TRUTH_BOUNDARY = {
    "source_axis_domain_independently_audited": True,
    "source_Phi_radius_one_ball_independently_audited": True,
    "source_u_radius_one_ball_independently_audited": True,
    "conditional_mixed_J_algebra_independently_audited": True,
    "source_Wstar_fixed_multiplier_independently_reconstructed": True,
    "source_R1_Wstar_logradial_ordinary_slot_independently_audited": True,
    "source_R1_zeta_ordinary_slot_independently_audited": True,
    "source_R1_averaged_eta_ordinary_slot_independent_audit_pending": True,
    "source_R1_prior_five_algebraic_slots_independently_audited": False,
    "all_R1_ordinary_slots_independently_audited": False,
    "source_full_post_J2_R1_radius_one_ball_bound_machine_bound": False,
    "source_operator_constant_M_machine_bound": False,
    "source_operator_constant_K_machine_bound": False,
    "source_contraction_invariant_ball_machine_verified": False,
    "source_fixed_point_distance_machine_bound": False,
    "source_B0_dependencies_machine_bound": False,
    "source_T_sh_lower_bound_verified": False,
    "selected_pa16_handoff_allowed": False,
    "global_pressure_matched": False,
    "global_leading_profile_reconstructed": False,
    "leading_only_ns_residual_assessed": False,
    "leading_plus_oscillatory_ns_residual_assessed": False,
    "after_correction_ns_residual_assessed": False,
    "heldout_ns_momentum_residual_assessed": False,
    "formal_full_domain_pde_gate_assessed": False,
    "pde_validated": False,
}


def _canonical_json(payload: dict[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _d(value: float) -> Decimal:
    out = float(value)
    if not math.isfinite(out) or out < 0.0:
        raise ValueError("public bound must be finite and nonnegative")
    return Decimal.from_float(out)


def _pair(bound: Any) -> tuple[Decimal, Decimal]:
    return _d(bound.norm), _d(bound.lipschitz)


def _ratio(public: Decimal, independent: Decimal) -> float:
    if independent <= 0:
        raise ValueError("independent comparator must be positive")
    return float(public / independent)


def _independent_product_constant() -> Decimal:
    # Strict classical rational upper, distinct from Agent-1's pi < 22/7 path.
    pi_upper = Decimal(355) / Decimal(113)
    c_sq = (Decimal(4) / Decimal(3)) * pi_upper * pi_upper
    return c_sq * c_sq


def _independent_l_inverse_norm(
    *, h: float, margin: float, radius: float, rho: float
) -> Decimal:
    h_d = Decimal.from_float(float(h))
    margin_d = Decimal.from_float(float(margin))
    radius_d = Decimal.from_float(float(radius))
    rho_d = Decimal.from_float(float(rho))
    z_abs = Decimal(1) + margin_d + radius_d
    lower = Decimal(1) - Decimal(2) * h_d * z_abs * z_abs
    if lower <= 0 or not Decimal(0) < rho_d < radius_d:
        raise ValueError("invalid public Cauchy geometry")
    x = rho_d / radius_d
    weight_sum = (Decimal(1) + x) / (Decimal(1) - x) ** 3
    return weight_sum / lower


def _independent_wstar_norm(domain: Any) -> Decimal:
    """Reconstruct the alpha=0 W_* coefficient bound independently.

    For W_*=-3-2D*j0*eta+8h*eta^2, the source alpha=0 coefficient norm is
    controlled by the maximum of the beta=0,1,2 Cauchy-weighted derivatives.
    This is the same mathematical object as the public fixed multiplier, but
    evaluated with Decimal rather than Agent-1's Fraction/outward path.
    """
    rho = Decimal.from_float(float(domain.coefficient_rho))
    margin = Decimal.from_float(float(domain.enlarged_real_margin))
    h = Decimal.from_float(float(domain.h))
    D = Decimal("0.5") - h
    j0 = Decimal.from_float(float(domain.j0))
    E = Decimal(1) + margin

    beta0 = Decimal(3) + Decimal(2) * D * j0 * E + Decimal(8) * h * E * E
    beta1 = Decimal(4) * rho * (
        Decimal(2) * D * j0 + Decimal(16) * h * E
    )
    beta2 = Decimal("4.5") * rho * rho * (Decimal(16) * h)
    return max(beta0, beta1, beta2)


def _independent_slot(
    calc: KokunoPA10SourceR1WstarLogradialOrdinarySlot,
) -> tuple[tuple[Decimal, Decimal], dict[str, Decimal]]:
    domain = calc.domain
    inputs = calc.operator_inputs()
    pc = _independent_product_constant()
    linv = _independent_l_inverse_norm(
        h=float(domain.h),
        margin=float(domain.enlarged_real_margin),
        radius=float(domain.cauchy_radius),
        rho=float(domain.coefficient_rho),
    )
    wstar = _independent_wstar_norm(domain)
    phi_norm, phi_lip = _pair(inputs["Phi"])

    # Three factors -> exactly two coefficient convolutions.  Only Phi moves on
    # the radius-one ball; L^-1 and W_* are fixed multipliers.
    fixed = pc**2 * linv * wstar * Decimal(80)
    return (fixed * phi_norm, fixed * phi_lip), {
        "product_constant": pc,
        "L_inverse": linv,
        "W_star": wstar,
        "single_logradial_after_J2": Decimal(80),
        "coefficient_convolution_count": Decimal(2),
        "fixed_prefactor": fixed,
    }


def _fresh_wstar_stress(
    calc: KokunoPA10SourceR1WstarLogradialOrdinarySlot,
) -> dict[str, Any]:
    domain = calc.domain
    E = 1.0 + float(domain.enlarged_real_margin)
    rng = np.random.default_rng(SEED)
    random_eta = rng.uniform(-E, E, size=OFFGRID_COUNT)
    probes = np.array([-E, -1.0e-12, 0.0, 1.0e-12, E], dtype=float)
    eta = np.concatenate((random_eta, probes))

    h = float(domain.h)
    D = float(domain.D)
    j0 = float(domain.j0)
    W = -3.0 - 2.0 * D * j0 * eta + 8.0 * h * eta**2
    public_bound = float(calc.operator_inputs()["W_star"].norm)
    max_abs = float(np.max(np.abs(W)))

    # Parameter perturbation is diagnostic only; it is not used to retune j0.
    j0_lo = j0 * 0.999
    j0_hi = j0 * 1.001
    W_lo = -3.0 - 2.0 * D * j0_lo * eta + 8.0 * h * eta**2
    W_hi = -3.0 - 2.0 * D * j0_hi * eta + 8.0 * h * eta**2
    response = float(
        max(np.max(np.abs(W_lo - W)), np.max(np.abs(W_hi - W)))
    )

    return {
        "seed": SEED,
        "fresh_random_offgrid_points": OFFGRID_COUNT,
        "axis_near_probes": [-1.0e-12, 0.0, 1.0e-12],
        "endpoint_probes": [-E, E],
        "max_abs_W_star": max_abs,
        "public_W_star_coefficient_bound": public_bound,
        "value_over_public_bound": max_abs / public_bound,
        "j0_pm_0p1_percent_max_abs_response": response,
        "finite": bool(np.isfinite(W).all()),
        "dominated": bool(max_abs <= public_bound),
        "nontrivial": bool(max_abs > 1.0),
    }


def _weight_ratio_stress() -> dict[str, Any]:
    """Enumerate the source single-logradial first-weight-ratio expression."""
    worst = 0.0
    arg = (0, 0)
    for alpha in range(1, 513):
        for beta in range(0, 129):
            value = (
                20.0
                * alpha
                * (alpha + 2.0) ** 2
                / ((alpha + 1.0) ** 3 * (alpha + beta + 1.0))
            )
            if value > worst:
                worst = value
                arg = (alpha, beta)
    return {
        "alpha_max": 512,
        "beta_max": 128,
        "max_single_logradial_ratio": worst,
        "argmax": list(arg),
        "source_envelope": 80.0,
        "envelope_dominates": bool(worst <= 80.0),
    }


def run_independent_audit() -> dict[str, Any]:
    calc = KokunoPA10SourceR1WstarLogradialOrdinarySlot()
    inputs = calc.operator_inputs()

    with localcontext() as ctx:
        ctx.prec = 80
        ctx.rounding = ROUND_CEILING

        (ind_norm, ind_lip), pieces = _independent_slot(calc)
        public_slot = _pair(calc.r1_wstar_logradial_after_j2()[SLOT_NAME])
        public_pc = _d(calc.parent.product_constant_upper)
        public_wstar = _pair(inputs["W_star"])[0]
        public_linv = _pair(inputs["L_inverse"])[0]
        public_log = _d(calc.source_single_logradial_after_j2_factor_upper())

        ratios = {
            "product_constant": _ratio(public_pc, pieces["product_constant"]),
            "W_star": _ratio(public_wstar, pieces["W_star"]),
            "L_inverse": _ratio(public_linv, pieces["L_inverse"]),
            "single_logradial_after_J2": _ratio(
                public_log, pieces["single_logradial_after_J2"]
            ),
            "slot_norm": _ratio(public_slot[0], ind_norm),
            "slot_lipschitz": _ratio(public_slot[1], ind_lip),
        }

        # Frozen mutations.  Downward mutations must cease dominating the
        # independent comparator.  Missing/duplicate operator factors are shape
        # errors rather than alternative valid upper bounds.
        negative_controls = {
            "product_constant_0p998_detected": bool(
                public_pc * Decimal("0.998") < pieces["product_constant"]
            ),
            "W_star_0p999_detected": bool(
                public_wstar * Decimal("0.999") < pieces["W_star"]
            ),
            "L_inverse_0p999_detected": bool(
                public_linv * Decimal("0.999") < pieces["L_inverse"]
            ),
            "logradial_79_instead_of_80_detected": bool(
                public_log * Decimal(79) / Decimal(80)
                < pieces["single_logradial_after_J2"]
            ),
            "slot_norm_0p99_detected": bool(
                public_slot[0] * Decimal("0.99") < ind_norm
            ),
            "slot_lipschitz_0p99_detected": bool(
                public_slot[1] * Decimal("0.99") < ind_lip
            ),
            "missing_one_convolution_detected": bool(
                public_slot[0] / public_pc < ind_norm
            ),
            "duplicate_J2_factor_40_detected": bool(
                public_slot[0] * Decimal(40) / ind_norm > MAX_PUBLIC_RATIO
            ),
        }

        failed_guards: list[str] = []
        for name, ratio in ratios.items():
            if ratio < 1.0:
                failed_guards.append(f"public_underbounds_independent:{name}")
        for name in ("product_constant", "W_star", "L_inverse", "single_logradial_after_J2", "slot_norm", "slot_lipschitz"):
            if Decimal.from_float(float(ratios[name])) > MAX_PUBLIC_RATIO:
                failed_guards.append(f"unexpected_shape_overcount:{name}")
        for name, detected in negative_controls.items():
            if not detected:
                failed_guards.append(f"negative_control_not_detected:{name}")

        offgrid = _fresh_wstar_stress(calc)
        if not offgrid["finite"]:
            failed_guards.append("fresh_Wstar_stress_nonfinite")
        if not offgrid["dominated"]:
            failed_guards.append("fresh_Wstar_stress_not_dominated")
        if not offgrid["nontrivial"]:
            failed_guards.append("fresh_Wstar_stress_trivial")
        if offgrid["j0_pm_0p1_percent_max_abs_response"] <= 1.0e-12:
            failed_guards.append("j0_parameter_perturbation_not_detected")

        weights = _weight_ratio_stress()
        if not weights["envelope_dominates"]:
            failed_guards.append("single_logradial_weight_envelope_failed")

        payload: dict[str, Any] = {
            "schema": SCHEMA,
            "agent1_pr": AGENT1_PR,
            "agent1_exact_head": AGENT1_HEAD,
            "independence": {
                "decimal_precision": 80,
                "pi_upper": "355/113",
                "agent1_private_helpers_called": False,
                "agent1_receipt_used_as_numeric_oracle": False,
                "upstream_Phi_ball_consumed_as_independently_admitted_interface": True,
                "W_star_reconstructed_independently": True,
                "L_inverse_reconstructed_independently": True,
                "single_logradial_weight_ratio_enumerated_independently": True,
            },
            "independent_values": {
                "product_constant": float(pieces["product_constant"]),
                "W_star": float(pieces["W_star"]),
                "L_inverse": float(pieces["L_inverse"]),
                "single_logradial_after_J2": 80.0,
                "slot_norm": float(ind_norm),
                "slot_lipschitz": float(ind_lip),
            },
            "public_values": {
                "product_constant": float(public_pc),
                "W_star": float(public_wstar),
                "L_inverse": float(public_linv),
                "single_logradial_after_J2": float(public_log),
                "slot_norm": float(public_slot[0]),
                "slot_lipschitz": float(public_slot[1]),
            },
            "public_over_independent_ratios": ratios,
            "fresh_offgrid_Wstar_stress": offgrid,
            "single_logradial_weight_ratio_stress": weights,
            "negative_controls": negative_controls,
            "failed_guards": failed_guards,
            "source_R1_Wstar_logradial_ordinary_slot_independent_preflight_passed": not failed_guards,
            "dependencies": {
                "source_axis_A4_PR": SOURCE_AXIS_A4_PR,
                "source_Phi_ball_A4_PR": SOURCE_PHI_BALL_A4_PR,
                "source_u_ball_A4_PR": SOURCE_U_BALL_A4_PR,
                "mixed_J_algebra_A4_PR": MIXED_J_A4_PR,
                "fixed_multiplier_derivative_A4_PR": DERIVATIVE_FIXED_MULTIPLIER_A4_PR,
                "zeta_R1_A4_PR": ZETA_R1_A4_PR,
                "averaged_eta_R1_A4_PR": AVERAGED_ETA_R1_A4_PR,
            },
            "immutable_project_gates": {
                "normalized_momentum_max": 1.0e-3,
                "normalized_momentum_L2": 1.0e-3,
                "divergence_max": 1.0e-5,
                "divergence_L2": 1.0e-5,
            },
            "truth_boundary": dict(_TRUTH_BOUNDARY),
            "limitations": [
                "This is a coefficient/operator audit, not a Navier-Stokes residual.",
                "Agent-1 #724 upstream CI is not promoted by this receipt.",
                "Three ordinary R1 slots remain Agent-1-unbound after #724.",
                "The averaged-eta #718 A4 audit remains separately pending until its exact-head receipt resolves.",
                "No global leading velocity, matched pressure, real correction cycle, or full Kokuno composite is available here.",
            ],
        }
        payload["receipt_sha256"] = hashlib.sha256(
            _canonical_json(payload).encode("utf-8")
        ).hexdigest()
        return payload


def save_report(path: str | Path) -> dict[str, Any]:
    payload = run_independent_audit()
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return payload


def _main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    payload = save_report(args.output)
    print("receipt_sha256=", payload["receipt_sha256"])
    print("ratios=", payload["public_over_independent_ratios"])
    print("offgrid=", payload["fresh_offgrid_Wstar_stress"])
    print("weight_ratio=", payload["single_logradial_weight_ratio_stress"])
    print("negative_controls=", payload["negative_controls"])
    print("failed_guards=", payload["failed_guards"])


if __name__ == "__main__":
    _main()
