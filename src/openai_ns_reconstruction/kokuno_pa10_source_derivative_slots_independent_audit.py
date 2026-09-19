"""Independent Agent-4 audit of PA.10 derivative-bearing ordinary R2 slots.

This validator is deliberately implementation-distinct from Agent-1 #693.  It
consumes only public coefficient-ball/value interfaces, reconstructs the four
source derivative-bearing post-J1 slots with 80-digit Decimal arithmetic, and
uses the stricter rational upper pi < 355/113 rather than Agent-1's 22/7 path.
It never calls Agent-1 private Fraction/product/rounding helpers and never uses
an Agent-1 receipt as a numerical oracle.

Scope is intentionally narrow: a PASS can admit only #693's four derivative
ordinary R2 slots.  Full R2 remains conditional on the separate Agent-4 audit
of #684's algebraic four-slot subset; R1/M/K/global leading/pressure/PDE gates
remain closed.
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

from .kokuno_pa10_source_derivative_ordinary_slots import (
    KokunoPA10SourceDerivativeOrdinarySlots,
)

SCHEMA = "kokuno-agent4-pa10-source-derivative-slots-independent-audit-v1"
AGENT1_PR = 693
AGENT1_HEAD = "f2f79dc27259d35f0ca341d0a5fa03fe97d6c51e"
ALGEBRAIC_A4_PR = 695
SOURCE_AXIS_A4_PR = 671
PRESSURE_SLOTS_A4_PR = 678
SOURCE_U_BALL_A4_PR = 687
MIXED_J_A4_PR = 663
SEED = 9173371
OFFGRID_COUNT = 8192
MUTATION_FACTOR_FIXED = Decimal("0.999")
MUTATION_FACTOR_SLOT = Decimal("0.99")

_TRUTH_BOUNDARY = {
    "source_axis_domain_independently_audited": True,
    "source_R2_pressure_ordinary_slots_independently_audited": True,
    "source_u_radius_one_ball_independently_audited": True,
    "conditional_mixed_J_algebra_independently_audited": True,
    "source_R2_derivative_ordinary_slots_independently_audited": True,
    "source_R2_algebraic_ordinary_slots_independent_audit_pending": True,
    "all_R2_ordinary_slots_independently_audited": False,
    "all_R1_ordinary_slots_source_bound": False,
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


def _scale(pair: tuple[Decimal, Decimal], factor: Decimal) -> tuple[Decimal, Decimal]:
    if factor < 0:
        factor = -factor
    return pair[0] * factor, pair[1] * factor


def _add(*pairs: tuple[Decimal, Decimal]) -> tuple[Decimal, Decimal]:
    return sum((p[0] for p in pairs), Decimal(0)), sum(
        (p[1] for p in pairs), Decimal(0)
    )


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
    # Strict classical rational upper, independent from Agent-1's pi < 22/7.
    pi_upper = Decimal(355) / Decimal(113)
    c_sq = (Decimal(4) / Decimal(3)) * pi_upper * pi_upper
    return c_sq * c_sq


def _independent_l_inverse_norm(*, h: float, margin: float, radius: float, rho: float) -> Decimal:
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


def _independent_w_h_bounds(domain: Any) -> tuple[Decimal, Decimal]:
    rho = Decimal.from_float(float(domain.coefficient_rho))
    margin = Decimal.from_float(float(domain.enlarged_real_margin))
    h = Decimal.from_float(float(domain.h))
    D = Decimal("0.5") - h
    j0 = Decimal.from_float(float(domain.j0))
    E = Decimal(1) + margin

    w0 = Decimal(3) + Decimal(2) * D * j0 * E + Decimal(8) * h * E * E
    w1 = Decimal(4) * rho * (Decimal(2) * D * j0 + Decimal(16) * h * E)
    w2 = Decimal("4.5") * rho * rho * (Decimal(16) * h)
    w_bound = max(w0, w1, w2)

    h0 = j0 + (D + Decimal(4)) * E + j0 * E * E + Decimal(4) * E**3
    h1 = Decimal(4) * rho * (
        D + Decimal(4) + Decimal(2) * j0 * E + Decimal(12) * E * E
    )
    h2 = Decimal("4.5") * rho * rho * (Decimal(2) * j0 + Decimal(24) * E)
    h3 = (Decimal(8) / Decimal(3)) * rho**3 * Decimal(24)
    h_bound = max(h0, h1, h2, h3)
    return w_bound, h_bound


def _fresh_polynomial_stress(calc: KokunoPA10SourceDerivativeOrdinarySlots) -> dict[str, Any]:
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
    H = j0 + (D + 4.0) * eta - j0 * eta**2 - 4.0 * eta**3
    public_w = float(calc.wstar_coefficient_ball().norm)
    public_h = float(calc.hstar_coefficient_ball().norm)
    max_w = float(np.max(np.abs(W)))
    max_h = float(np.max(np.abs(H)))

    return {
        "seed": SEED,
        "fresh_random_offgrid_points": OFFGRID_COUNT,
        "axis_near_probes": [-1.0e-12, 0.0, 1.0e-12],
        "endpoint_probes": [-E, E],
        "max_abs_W_star": max_w,
        "max_abs_H_star": max_h,
        "public_W_star_bound": public_w,
        "public_H_star_bound": public_h,
        "W_star_stress_ratio": max_w / public_w,
        "H_star_stress_ratio": max_h / public_h,
        "finite": bool(np.isfinite(W).all() and np.isfinite(H).all()),
        "dominated": bool(max_w <= public_w and max_h <= public_h),
        "nontrivial": bool(max_w > 1.0 and max_h > 1.0),
    }


def _weight_ratio_stress() -> dict[str, Any]:
    # Directly enumerate the exact single-logradial J1 weight-ratio expression
    # documented by the public source-facing implementation.  The shared source
    # envelope remains 80; enumeration is an implementation-distinct sanity check.
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
        "shared_source_envelope": 80.0,
        "envelope_dominates": bool(worst <= 80.0),
    }


def run_independent_audit() -> dict[str, Any]:
    calc = KokunoPA10SourceDerivativeOrdinarySlots()
    domain = calc.domain

    with localcontext() as ctx:
        ctx.prec = 80
        ctx.rounding = ROUND_CEILING

        pc = _independent_product_constant()
        rho = Decimal.from_float(float(domain.coefficient_rho))
        D = Decimal.from_float(float(domain.D))
        linv = _independent_l_inverse_norm(
            h=float(domain.h),
            margin=float(domain.enlarged_real_margin),
            radius=float(domain.cauchy_radius),
            rho=float(domain.coefficient_rho),
        )
        W_ind, H_ind = _independent_w_h_bounds(domain)

        u = _pair(calc.bridge.u_source.u_ball())
        eta = _pair(calc.algebraic.fixed_multiplier_balls()["eta"])
        d = (Decimal(1), Decimal(0))
        Linv = (linv, Decimal(0))
        W = (W_ind, Decimal(0))
        H = (H_ind, Decimal(0))

        j1_log_u = _scale(u, Decimal(80))
        j1_eta_u = _scale(u, Decimal(80) / rho)

        Au = u
        Au_u = _product(Au, u, product_constant=pc)
        u2 = _product(u, u, product_constant=pc)
        j1_Au_YuY = _add(
            _scale(Au_u, Decimal(80)),
            _scale(u2, Decimal(80)),
            _scale(Au_u, Decimal(80)),
        )
        j1_u_ueta = _scale(u2, Decimal(40) / rho)

        independent_slots = {
            "W_star_times_Y_u_Y": _product(Linv, W, j1_log_u, product_constant=pc),
            "lambda_inv_2D_eta_Au_times_Y_u_Y": _scale(
                _product(Linv, eta, j1_Au_YuY, product_constant=pc),
                Decimal(2) * D,
            ),
            "H_star_times_u_eta": _product(Linv, H, j1_eta_u, product_constant=pc),
            "lambda_inv_d_u_u_eta": _product(Linv, d, j1_u_ueta, product_constant=pc),
        }

        public_slots = calc.r2_derivative_after_j1()
        slot_ratios: dict[str, dict[str, float]] = {}
        for name, ind in independent_slots.items():
            pub = _pair(public_slots[name])
            slot_ratios[name] = {
                "norm": _ratio(pub[0], ind[0]),
                "lipschitz": _ratio(pub[1], ind[1]),
            }

        public_w = _pair(calc.wstar_coefficient_ball())[0]
        public_h = _pair(calc.hstar_coefficient_ball())[0]
        public_linv = _pair(calc.bridge.l_inverse_coefficient_ball())[0]
        fixed_ratios = {
            "W_star": _ratio(public_w, W_ind),
            "H_star": _ratio(public_h, H_ind),
            "L_inverse": _ratio(public_linv, linv),
        }

        stress = _fresh_polynomial_stress(calc)
        weight_ratio = _weight_ratio_stress()

        mutation = {
            "fixed_0p999_detected": {
                "W_star": bool(public_w * MUTATION_FACTOR_FIXED < W_ind),
                "H_star": bool(public_h * MUTATION_FACTOR_FIXED < H_ind),
                "L_inverse": bool(public_linv * MUTATION_FACTOR_FIXED < linv),
            },
            "slot_0p99_detected": {},
        }
        for name, ind in independent_slots.items():
            pub = _pair(public_slots[name])
            mutation["slot_0p99_detected"][name] = bool(
                pub[0] * MUTATION_FACTOR_SLOT < ind[0]
                and pub[1] * MUTATION_FACTOR_SLOT < ind[1]
            )

        failed_guards: list[str] = []
        if min(fixed_ratios.values()) < 1.0:
            failed_guards.append("fixed coefficient/Cauchy bound under independent reconstruction")
        for name, ratios in slot_ratios.items():
            if min(ratios.values()) < 1.0:
                failed_guards.append(f"{name} under independent reconstruction")
        if not stress["finite"] or not stress["dominated"] or not stress["nontrivial"]:
            failed_guards.append("fresh W_star/H_star off-grid stress failed")
        if not weight_ratio["envelope_dominates"]:
            failed_guards.append("single-logradial source envelope failed")
        if not all(mutation["fixed_0p999_detected"].values()):
            failed_guards.append("fixed-bound negative control not detected")
        if not all(mutation["slot_0p99_detected"].values()):
            failed_guards.append("derivative-slot negative control not detected")

        independent_values = {
            "product_constant_pi355_over_113": float(pc),
            "rho": float(rho),
            "L_inverse": float(linv),
            "W_star_coefficient_bound": float(W_ind),
            "H_star_coefficient_bound": float(H_ind),
            "source_u_ball_public_upstream": {
                "norm": float(u[0]),
                "lipschitz": float(u[1]),
                "independent_admission_PR": SOURCE_U_BALL_A4_PR,
            },
            "R2_derivative_ordinary_after_J1": {
                name: {"norm": float(pair[0]), "lipschitz": float(pair[1])}
                for name, pair in independent_slots.items()
            },
        }

    payload: dict[str, Any] = {
        "schema": SCHEMA,
        "agent1_pr": AGENT1_PR,
        "agent1_exact_head": AGENT1_HEAD,
        "seed": SEED,
        "independent_path": {
            "arithmetic": "Decimal precision=80 ROUND_CEILING",
            "product_constant": "pi<355/113, not Agent-1 pi<22/7",
            "private_agent1_helpers_called": False,
            "agent1_receipt_used_as_oracle": False,
            "training_internal_tensors_read": False,
        },
        "independent_values": independent_values,
        "public_over_independent_ratios": {
            "fixed": fixed_ratios,
            "derivative_slots": slot_ratios,
        },
        "fresh_offgrid_polynomial_stress": stress,
        "weight_ratio_stress": weight_ratio,
        "negative_controls": mutation,
        "failed_guards": failed_guards,
        "source_R2_derivative_ordinary_slots_independent_preflight_passed": not failed_guards,
        "dependency_state": {
            "source_axis_A4_PR": SOURCE_AXIS_A4_PR,
            "pressure_slots_A4_PR": PRESSURE_SLOTS_A4_PR,
            "source_u_ball_A4_PR": SOURCE_U_BALL_A4_PR,
            "mixed_J_A4_PR": MIXED_J_A4_PR,
            "algebraic_slots_A4_PR": ALGEBRAIC_A4_PR,
            "full_R2_independent_admission_requires_algebraic_PR695_PASS": True,
        },
        "immutable_project_gates": {
            "normalized_momentum_max": 1.0e-3,
            "normalized_momentum_L2": 1.0e-3,
            "divergence_max": 1.0e-5,
            "divergence_L2": 1.0e-5,
        },
        "truth_boundary": dict(_TRUTH_BOUNDARY),
    }
    receipt_basis = dict(payload)
    payload["receipt_sha256"] = hashlib.sha256(
        _canonical_json(receipt_basis).encode()
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
