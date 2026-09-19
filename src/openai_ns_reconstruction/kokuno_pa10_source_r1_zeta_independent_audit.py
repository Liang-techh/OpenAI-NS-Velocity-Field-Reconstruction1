"""Independent Agent-4 audit of the PA.10 ``d u zeta_* Phi`` R1 slot.

This validator is deliberately implementation-distinct from Agent-1 #707.  It
consumes only public value/bound interfaces, reconstructs the fixed multiplier
``zeta_*=-L H_*/(H_*^2+sigma_*^2)`` with 80-digit Decimal arithmetic, stress
samples the actual complex Cauchy neighborhood, and rebuilds the numerator and
post-J2 bounds with the stricter rational upper ``pi < 355/113`` rather than
Agent-1's ``22/7`` product path.  It never calls Agent-1 private Fraction or
rounding helpers and never uses an Agent-1 receipt as a numerical oracle.

Scope is intentionally narrow: a PASS can independently admit only #707's one
R1 zeta ordinary slot, conditional on the already admitted source-axis, Phi and
u-ball prerequisites.  The other R1 ordinary slots, full R1/R2, M/K, global
leading/pressure and PDE gates remain closed.
"""

from __future__ import annotations

import argparse
from decimal import Decimal, ROUND_CEILING, ROUND_FLOOR, localcontext
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from .kokuno_pa10_source_r1_zeta_ordinary_slot import (
    KokunoPA10SourceR1ZetaOrdinarySlot,
)

SCHEMA = "kokuno-agent4-pa10-source-r1-zeta-independent-audit-v1"
AGENT1_PR = 707
AGENT1_HEAD = "3068062f92dcf9c9eb1276a014bed1f36bb8aa86"
SOURCE_AXIS_A4_PR = 671
SOURCE_PHI_A4_PR = 655
SOURCE_U_A4_PR = 687
R2_ALGEBRAIC_A4_PR = 695
R2_DERIVATIVE_A4_PR = 704
SEED = 9173381
REAL_OFFGRID_COUNT = 8192
COMPLEX_OFFGRID_COUNT = 20000
MUTATION_FIXED = Decimal("0.999")
MUTATION_SLOT = Decimal("0.99")

_TRUTH_BOUNDARY = {
    "source_axis_domain_independently_audited": True,
    "source_Phi_radius_one_ball_independently_audited": True,
    "source_u_radius_one_ball_independently_audited": True,
    "source_R1_zeta_ordinary_slot_independently_audited": True,
    "source_R1_other_ordinary_slots_independently_audited": False,
    "source_full_post_J2_R1_radius_one_ball_bound_machine_bound": False,
    "full_R1_independent_agent4_admission": False,
    "full_R2_independent_agent4_admission": False,
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


def _scale(pair: tuple[Decimal, Decimal], factor: Decimal) -> tuple[Decimal, Decimal]:
    factor = abs(factor)
    return pair[0] * factor, pair[1] * factor


def _product(
    *pairs: tuple[Decimal, Decimal], product_constant: Decimal
) -> tuple[Decimal, Decimal]:
    if not pairs:
        return Decimal(1), Decimal(0)
    algebra = product_constant ** max(len(pairs) - 1, 0)
    norm_product = Decimal(1)
    for norm, _ in pairs:
        norm_product *= norm
    lip_sum = Decimal(0)
    for i, (_, lip) in enumerate(pairs):
        if lip == 0:
            continue
        others = Decimal(1)
        for j, (norm, _) in enumerate(pairs):
            if i != j:
                others *= norm
        lip_sum += lip * others
    return algebra * norm_product, algebra * lip_sum


def _ratio(public: Decimal, independent: Decimal) -> float:
    if independent <= 0:
        raise ValueError("independent comparator must be positive")
    return float(public / independent)


def _independent_product_constant() -> Decimal:
    pi_upper = Decimal(355) / Decimal(113)
    c_sq = (Decimal(4) / Decimal(3)) * pi_upper * pi_upper
    return c_sq * c_sq


def _independent_l_inverse_norm(*, h: float, margin: float, radius: float, rho: float) -> Decimal:
    with localcontext() as ctx:
        ctx.prec = 80
        ctx.rounding = ROUND_CEILING
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


def _independent_zeta_majorant(domain: Any) -> dict[str, Decimal]:
    rho = Decimal.from_float(float(domain.coefficient_rho))
    radius = Decimal.from_float(float(domain.cauchy_radius))
    margin = Decimal.from_float(float(domain.enlarged_real_margin))
    sigma = Decimal.from_float(float(domain.sigma_star))
    h = Decimal.from_float(float(domain.h))
    j0 = Decimal.from_float(float(domain.j0))
    if not Decimal(0) < rho < radius or sigma <= 0:
        raise ValueError("invalid source-compatible zeta geometry")

    with localcontext() as up:
        up.prec = 80
        up.rounding = ROUND_CEILING
        D = Decimal("0.5") - h
        z_abs = Decimal(1) + margin + radius
        j0_abs = abs(j0)
        H_upper = (
            j0_abs
            + (D + Decimal(4)) * z_abs
            + j0_abs * z_abs * z_abs
            + Decimal(4) * z_abs**3
        )
        H_prime_upper = (
            D
            + Decimal(4)
            + Decimal(2) * j0_abs * z_abs
            + Decimal(12) * z_abs**2
        )
        H_variation_upper = radius * H_prime_upper
        L_upper = Decimal(1) + Decimal(2) * h * z_abs**2

    with localcontext() as down:
        down.prec = 80
        down.rounding = ROUND_FLOOR
        H_factor_lower = sigma - H_variation_upper
        if H_factor_lower <= 0:
            raise ValueError("Cauchy H variation reaches sigma_star")
        denominator_lower = H_factor_lower * H_factor_lower

    with localcontext() as up:
        up.prec = 80
        up.rounding = ROUND_CEILING
        zeta_complex_sup = L_upper * H_upper / denominator_lower
        x = rho / radius
        weight_sum = (Decimal(1) + x) / (Decimal(1) - x) ** 3
        coefficient_norm = zeta_complex_sup * weight_sum

    return {
        "cauchy_z_abs_upper": z_abs,
        "cauchy_H_abs_upper": H_upper,
        "cauchy_H_prime_abs_upper": H_prime_upper,
        "cauchy_H_variation_abs_upper": H_variation_upper,
        "H_factor_abs_lower_on_cauchy": H_factor_lower,
        "H_square_plus_sigma_square_abs_lower_on_cauchy": denominator_lower,
        "cauchy_L_abs_upper": L_upper,
        "rho_over_cauchy_radius": x,
        "cauchy_weight_sum_upper": weight_sum,
        "zeta_star_complex_sup_upper": zeta_complex_sup,
        "zeta_star_coefficient_norm_upper": coefficient_norm,
    }


def _fresh_real_replay(calc: KokunoPA10SourceR1ZetaOrdinarySlot) -> dict[str, Any]:
    domain = calc.domain
    E = 1.0 + float(domain.enlarged_real_margin)
    rng = np.random.default_rng(SEED)
    random_eta = rng.uniform(-E, E, size=REAL_OFFGRID_COUNT)
    probes = np.array([-E, -1.0e-12, 0.0, 1.0e-12, E], dtype=float)
    eta = np.concatenate((random_eta, probes))

    h = float(domain.h)
    D = 0.5 - h
    j0 = float(domain.j0)
    sigma = float(domain.sigma_star)
    L = 1.0 - 2.0 * h * eta**2
    H = j0 + (D + 4.0) * eta - j0 * eta**2 - 4.0 * eta**3
    independent = -L * H / (H * H + sigma * sigma)
    public = np.asarray(calc.zeta_star_values(eta), dtype=float)
    scale = np.maximum(1.0, np.abs(independent))
    rel = np.abs(public - independent) / scale
    sign_mutation_rel = np.abs(-public - independent) / scale

    sigma_low = sigma * 0.999
    sigma_high = sigma * 1.001
    low = -L * H / (H * H + sigma_low * sigma_low)
    high = -L * H / (H * H + sigma_high * sigma_high)
    sensitivity = np.abs(high - low) / scale

    return {
        "seed": SEED,
        "fresh_random_offgrid_points": REAL_OFFGRID_COUNT,
        "axis_near_probes": [-1.0e-12, 0.0, 1.0e-12],
        "endpoint_probes": [-E, E],
        "relative_rms_mismatch": float(np.sqrt(np.mean(rel * rel))),
        "relative_max_mismatch": float(np.max(rel)),
        "zeta_rms": float(np.sqrt(np.mean(independent * independent))),
        "zeta_max_abs": float(np.max(np.abs(independent))),
        "sign_flip_relative_rms": float(np.sqrt(np.mean(sign_mutation_rel**2))),
        "sigma_pm_0p1_percent_max_relative_response": float(np.max(sensitivity)),
        "finite": bool(
            np.isfinite(independent).all()
            and np.isfinite(public).all()
            and np.isfinite(low).all()
            and np.isfinite(high).all()
        ),
    }


def _fresh_complex_stress(
    calc: KokunoPA10SourceR1ZetaOrdinarySlot, independent_complex_sup: Decimal
) -> dict[str, Any]:
    domain = calc.domain
    E = 1.0 + float(domain.enlarged_real_margin)
    radius = float(domain.cauchy_radius)
    rng = np.random.default_rng(SEED + 1)
    centers = rng.uniform(-E, E, size=COMPLEX_OFFGRID_COUNT)
    radii = radius * np.sqrt(rng.random(COMPLEX_OFFGRID_COUNT))
    angles = rng.uniform(0.0, 2.0 * np.pi, size=COMPLEX_OFFGRID_COUNT)
    z = centers + radii * np.exp(1j * angles)

    h = float(domain.h)
    D = 0.5 - h
    j0 = float(domain.j0)
    sigma = float(domain.sigma_star)
    L = 1.0 - 2.0 * h * z**2
    H = j0 + (D + 4.0) * z - j0 * z**2 - 4.0 * z**3
    zeta = -L * H / (H * H + sigma * sigma)
    max_abs = float(np.max(np.abs(zeta)))
    upper = float(independent_complex_sup)
    return {
        "seed": SEED + 1,
        "fresh_complex_offgrid_points": COMPLEX_OFFGRID_COUNT,
        "max_abs_zeta": max_abs,
        "independent_complex_sup_upper": upper,
        "stress_ratio": max_abs / upper,
        "finite": bool(np.isfinite(zeta.real).all() and np.isfinite(zeta.imag).all()),
        "dominated": bool(max_abs <= upper),
    }


def run_independent_audit() -> dict[str, Any]:
    calc = KokunoPA10SourceR1ZetaOrdinarySlot()
    domain = calc.domain

    with localcontext() as ctx:
        ctx.prec = 80
        ctx.rounding = ROUND_CEILING
        pc = _independent_product_constant()
        majorant = _independent_zeta_majorant(domain)
        linv = _independent_l_inverse_norm(
            h=float(domain.h),
            margin=float(domain.enlarged_real_margin),
            radius=float(domain.cauchy_radius),
            rho=float(domain.coefficient_rho),
        )

        d = (Decimal(1), Decimal(0))
        u = _pair(calc.algebraic_r1.u_ball())
        phi = _pair(calc.algebraic_r1.phi_ball())
        zeta = (majorant["zeta_star_coefficient_norm_upper"], Decimal(0))
        numerator_ind = _product(d, u, zeta, phi, product_constant=pc)
        after_j2_ind = _scale(
            _product((linv, Decimal(0)), numerator_ind, product_constant=pc),
            Decimal(40),
        )

        public_cert = calc.zeta_star_certificate()
        public_zeta = _pair(calc.zeta_star_coefficient_ball())
        public_num = _pair(calc.r1_zeta_numerator_slot()["d_u_zeta_star_times_Phi"])
        public_after = _pair(calc.r1_zeta_after_j2()["d_u_zeta_star_times_Phi"])

        fixed_ratios = {
            "zeta_complex_sup": _ratio(
                _d(public_cert["zeta_star_complex_sup_upper"]),
                majorant["zeta_star_complex_sup_upper"],
            ),
            "zeta_coefficient_norm": _ratio(
                public_zeta[0], majorant["zeta_star_coefficient_norm_upper"]
            ),
        }
        slot_ratios = {
            "numerator_norm": _ratio(public_num[0], numerator_ind[0]),
            "numerator_lipschitz": _ratio(public_num[1], numerator_ind[1]),
            "post_J2_norm": _ratio(public_after[0], after_j2_ind[0]),
            "post_J2_lipschitz": _ratio(public_after[1], after_j2_ind[1]),
        }

        real_replay = _fresh_real_replay(calc)
        complex_stress = _fresh_complex_stress(
            calc, majorant["zeta_star_complex_sup_upper"]
        )

        negative_controls = {
            "fixed_0p999_detected": {
                "zeta_complex_sup": bool(
                    _d(public_cert["zeta_star_complex_sup_upper"]) * MUTATION_FIXED
                    < majorant["zeta_star_complex_sup_upper"]
                ),
                "zeta_coefficient_norm": bool(
                    public_zeta[0] * MUTATION_FIXED
                    < majorant["zeta_star_coefficient_norm_upper"]
                ),
            },
            "slot_0p99_detected": {
                "numerator": bool(
                    public_num[0] * MUTATION_SLOT < numerator_ind[0]
                    and public_num[1] * MUTATION_SLOT < numerator_ind[1]
                ),
                "post_J2": bool(
                    public_after[0] * MUTATION_SLOT < after_j2_ind[0]
                    and public_after[1] * MUTATION_SLOT < after_j2_ind[1]
                ),
            },
            "real_value_sign_flip_detected": bool(
                real_replay["sign_flip_relative_rms"] > 1.0e-3
            ),
        }

        failed_guards: list[str] = []
        if min(fixed_ratios.values()) < 1.0:
            failed_guards.append("zeta fixed bound under independent reconstruction")
        if min(slot_ratios.values()) < 1.0:
            failed_guards.append("R1 zeta slot under independent reconstruction")
        if not real_replay["finite"]:
            failed_guards.append("fresh real off-grid zeta replay nonfinite")
        if real_replay["relative_max_mismatch"] > 1.0e-12:
            failed_guards.append("fresh real off-grid zeta replay mismatch")
        if real_replay["zeta_rms"] <= 1.0e-8:
            failed_guards.append("fresh real off-grid zeta replay trivial")
        if real_replay["sigma_pm_0p1_percent_max_relative_response"] <= 1.0e-8:
            failed_guards.append("zeta parameter perturbation response trivial")
        if not complex_stress["finite"] or not complex_stress["dominated"]:
            failed_guards.append("fresh complex Cauchy-neighborhood zeta stress failed")
        if not all(negative_controls["fixed_0p999_detected"].values()):
            failed_guards.append("fixed-bound negative control not detected")
        if not all(negative_controls["slot_0p99_detected"].values()):
            failed_guards.append("slot negative control not detected")
        if not negative_controls["real_value_sign_flip_detected"]:
            failed_guards.append("real-value sign mutation not detected")

        independent_values = {
            "product_constant_pi355_over_113": float(pc),
            "L_inverse_coefficient_norm_upper": float(linv),
            "zeta_majorant": {name: float(value) for name, value in majorant.items()},
            "upstream_public_A4_admitted_inputs": {
                "u_ball": {"norm": float(u[0]), "lipschitz": float(u[1])},
                "Phi_ball": {"norm": float(phi[0]), "lipschitz": float(phi[1])},
            },
            "R1_zeta_numerator_independent": {
                "norm": float(numerator_ind[0]),
                "lipschitz": float(numerator_ind[1]),
            },
            "R1_zeta_post_J2_independent": {
                "norm": float(after_j2_ind[0]),
                "lipschitz": float(after_j2_ind[1]),
            },
        }

    payload: dict[str, Any] = {
        "schema": SCHEMA,
        "agent1_pr": AGENT1_PR,
        "agent1_exact_head": AGENT1_HEAD,
        "seed": SEED,
        "independent_path": {
            "arithmetic": "Decimal precision=80 with directed rounding",
            "product_constant": "pi<355/113, not Agent-1 pi<22/7",
            "fresh_real_offgrid_points": REAL_OFFGRID_COUNT,
            "fresh_complex_offgrid_points": COMPLEX_OFFGRID_COUNT,
            "private_agent1_helpers_called": False,
            "agent1_receipt_used_as_oracle": False,
            "training_internal_tensors_read": False,
        },
        "independent_values": independent_values,
        "public_over_independent_ratios": {
            "fixed": fixed_ratios,
            "R1_zeta_slot": slot_ratios,
        },
        "fresh_real_offgrid_replay": real_replay,
        "fresh_complex_cauchy_stress": complex_stress,
        "negative_controls": negative_controls,
        "failed_guards": failed_guards,
        "source_R1_zeta_ordinary_slot_independent_preflight_passed": not failed_guards,
        "dependency_state": {
            "source_axis_A4_PR": SOURCE_AXIS_A4_PR,
            "source_Phi_ball_A4_PR": SOURCE_PHI_A4_PR,
            "source_u_ball_A4_PR": SOURCE_U_A4_PR,
            "R2_algebraic_A4_PR": R2_ALGEBRAIC_A4_PR,
            "R2_derivative_A4_PR": R2_DERIVATIVE_A4_PR,
            "full_R2_independent_admission_still_pending": True,
            "other_R1_ordinary_slots_independent_audit_still_required": True,
        },
        "immutable_project_gates": {
            "normalized_momentum_max": 1.0e-3,
            "normalized_momentum_L2": 1.0e-3,
            "divergence_max": 1.0e-5,
            "divergence_L2": 1.0e-5,
        },
        "truth_boundary": dict(_TRUTH_BOUNDARY),
    }
    basis = dict(payload)
    payload["receipt_sha256"] = hashlib.sha256(
        _canonical_json(basis).encode("utf-8")
    ).hexdigest()
    return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    payload = run_independent_audit()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    if payload["failed_guards"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
