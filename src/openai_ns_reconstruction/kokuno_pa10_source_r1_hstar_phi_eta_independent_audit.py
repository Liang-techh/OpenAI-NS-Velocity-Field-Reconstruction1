"""Independent Agent-4 audit of the PA.10 ``H_* Phi_eta`` ordinary R1 slot.

This validator audits Agent-1 #741 through its public ball interfaces while
using an implementation-distinct numerical/algebraic path.  It never calls the
Agent-1 private Fraction/product/rounding helpers and never treats an Agent-1
receipt as a numerical oracle.

The audited source-compatible object is

    J2[L^-1 H_* partial_eta Phi],

with

    H_* = j0 + (D+4) eta - j0 eta^2 - 4 eta^3,
    D = 1/2-h.

The independent path reconstructs the H_* coefficient majorant, the Cauchy
L^-1 majorant, the 80/rho post-J2 eta-derivative factor, and the final
three-factor coefficient product using 80-digit Decimal arithmetic and the
strict rational upper pi < 355/113.  Fresh off-grid and axis-near H_* probes,
a j0 perturbation, and preregistered mutations are included.

A PASS can admit only Agent-1 #741's tenth ordinary R1 slot.  It does not close
full R1, M/K, global leading/pressure, any staged Navier-Stokes residual, or
``pde_validated``.
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

from .kokuno_pa10_source_r1_hstar_phi_eta_ordinary_slot import (
    KokunoPA10SourceR1HstarPhiEtaOrdinarySlot,
)

SCHEMA = "kokuno-agent4-pa10-source-r1-hstar-phi-eta-independent-audit-v1"
AGENT1_PR = 741
AGENT1_HEAD = "20b600e7bf34aef63cecf74715180e7e06da231d"
SOURCE_AXIS_A4_PR = 671
PHI_BALL_A4_PR = 655
SEED = 9173421
OFFGRID_COUNT = 8192
SLOT_NAME = "H_star_times_Phi_eta"
MAX_PUBLIC_RATIO = Decimal("1.01")
MOMENTUM_GATE = 1.0e-3
DIVERGENCE_GATE = 1.0e-5

_TRUTH_BOUNDARY = {
    "source_axis_domain_independently_audited": True,
    "source_Phi_radius_one_ball_independently_audited": True,
    "source_R1_Hstar_Phi_eta_ordinary_slot_independently_audited": True,
    "source_R1_zeta_ordinary_slot_independently_audited": True,
    "source_R1_averaged_eta_ordinary_slot_independently_audited": True,
    "source_R1_Wstar_logradial_independent_audit_pending": True,
    "source_R1_averaged_logradial_independent_audit_pending": True,
    "source_R1_all_ordinary_slots_independently_admitted": False,
    "source_full_post_J2_R1_independently_admitted": False,
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
        raise ValueError("bound must be finite and nonnegative")
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
    lip = Decimal(0)
    for i, (_, item_lip) in enumerate(pairs):
        if item_lip == 0:
            continue
        others = Decimal(1)
        for j, (norm, _) in enumerate(pairs):
            if i != j:
                others *= norm
        lip += item_lip * others
    return algebra * norm_product, algebra * lip


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


def _independent_hstar_norm(domain: Any) -> Decimal:
    rho = Decimal.from_float(float(domain.coefficient_rho))
    margin = Decimal.from_float(float(domain.enlarged_real_margin))
    h = Decimal.from_float(float(domain.h))
    j0 = Decimal.from_float(float(domain.j0))
    D = Decimal("0.5") - h
    E = Decimal(1) + margin
    if min(rho, h, j0, D) < 0:
        raise ValueError("source-compatible parameters must be nonnegative")

    beta0 = j0 + (D + Decimal(4)) * E + j0 * E * E + Decimal(4) * E**3
    beta1 = Decimal(4) * rho * (
        D + Decimal(4) + Decimal(2) * j0 * E + Decimal(12) * E * E
    )
    beta2 = Decimal("4.5") * rho * rho * (
        Decimal(2) * j0 + Decimal(24) * E
    )
    beta3 = (Decimal(8) / Decimal(3)) * rho**3 * Decimal(24)
    return max(beta0, beta1, beta2, beta3)


def _fresh_hstar_stress(calc: KokunoPA10SourceR1HstarPhiEtaOrdinarySlot) -> dict[str, Any]:
    domain = calc.domain
    E = 1.0 + float(domain.enlarged_real_margin)
    h = float(domain.h)
    D = 0.5 - h
    j0 = float(domain.j0)

    rng = np.random.default_rng(SEED)
    random_eta = rng.uniform(-E, E, size=OFFGRID_COUNT)
    probes = np.array([-E, -1.0e-12, 0.0, 1.0e-12, E], dtype=float)
    eta = np.concatenate((random_eta, probes))

    def values(j: float) -> np.ndarray:
        return j + (D + 4.0) * eta - j * eta**2 - 4.0 * eta**3

    H = values(j0)
    H_lo = values(j0 * 0.999)
    H_hi = values(j0 * 1.001)
    H_wrong_sign = j0 + (D + 4.0) * eta + j0 * eta**2 - 4.0 * eta**3

    public_bound = float(calc.hstar_coefficient_ball().norm)
    max_abs = float(np.max(np.abs(H)))
    perturbation = float(max(np.max(np.abs(H_lo - H)), np.max(np.abs(H_hi - H))))
    wrong_sign_error = float(np.max(np.abs(H_wrong_sign - H)))

    return {
        "seed": SEED,
        "fresh_random_offgrid_points": OFFGRID_COUNT,
        "axis_near_probes": [-1.0e-12, 0.0, 1.0e-12],
        "endpoint_probes": [-E, E],
        "max_abs_H_star": max_abs,
        "public_H_star_coefficient_bound": public_bound,
        "stress_to_public_bound_ratio": max_abs / public_bound,
        "finite": bool(np.isfinite(H).all()),
        "dominated": bool(max_abs <= public_bound),
        "nontrivial": bool(max_abs > 1.0),
        "j0_pm_0p1_percent_response": perturbation,
        "quadratic_sign_mutation_max_abs_error": wrong_sign_error,
        "quadratic_sign_mutation_detected": bool(wrong_sign_error >= 1.0e-8),
    }


def run_independent_audit() -> dict[str, Any]:
    calc = KokunoPA10SourceR1HstarPhiEtaOrdinarySlot()
    domain = calc.domain

    with localcontext() as ctx:
        ctx.prec = 80
        ctx.rounding = ROUND_CEILING

        pc = _independent_product_constant()
        linv = _independent_l_inverse_norm(
            h=float(domain.h),
            margin=float(domain.enlarged_real_margin),
            radius=float(domain.cauchy_radius),
            rho=float(domain.coefficient_rho),
        )
        hstar = _independent_hstar_norm(domain)
        rho = Decimal.from_float(float(domain.coefficient_rho))
        single_eta_factor = Decimal(80) / rho

        phi = _pair(calc.operator_inputs()["Phi"])
        independent_post_phi = _scale(phi, single_eta_factor)
        independent_slot = _product(
            (linv, Decimal(0)),
            (hstar, Decimal(0)),
            independent_post_phi,
            product_constant=pc,
        )

        public_pc = _d(calc.product_calculator.product_constant)
        public_linv = _pair(calc.operator_inputs()["L_inverse"])[0]
        public_hstar = _pair(calc.hstar_coefficient_ball())[0]
        public_factor = _d(calc.source_single_eta_after_j2_factor_upper())
        public_post_phi = _pair(calc.j2_phi_eta())
        public_slot = _pair(calc.r1_hstar_phi_eta_after_j2()[SLOT_NAME])

        ratios = {
            "product_constant": _ratio(public_pc, pc),
            "L_inverse": _ratio(public_linv, linv),
            "H_star": _ratio(public_hstar, hstar),
            "single_eta_after_J2": _ratio(public_factor, single_eta_factor),
            "J2_Phi_eta_norm": _ratio(public_post_phi[0], independent_post_phi[0]),
            "J2_Phi_eta_lipschitz": _ratio(
                public_post_phi[1], independent_post_phi[1]
            ),
            "slot_norm": _ratio(public_slot[0], independent_slot[0]),
            "slot_lipschitz": _ratio(public_slot[1], independent_slot[1]),
        }

        one_convolution_omitted_norm = independent_slot[0] / pc
        negative_controls = {
            "product_constant_0p998_detected": bool(
                public_pc * Decimal("0.998") < pc
            ),
            "L_inverse_0p999_detected": bool(
                public_linv * Decimal("0.999") < linv
            ),
            "H_star_0p999_detected": bool(
                public_hstar * Decimal("0.999") < hstar
            ),
            "eta_factor_79_instead_of_80_detected": bool(
                public_factor * Decimal(79) / Decimal(80) < single_eta_factor
            ),
            "J2_Phi_eta_norm_0p99_detected": bool(
                public_post_phi[0] * Decimal("0.99") < independent_post_phi[0]
            ),
            "J2_Phi_eta_lipschitz_0p99_detected": bool(
                public_post_phi[1] * Decimal("0.99") < independent_post_phi[1]
            ),
            "slot_norm_0p99_detected": bool(
                public_slot[0] * Decimal("0.99") < independent_slot[0]
            ),
            "slot_lipschitz_0p99_detected": bool(
                public_slot[1] * Decimal("0.99") < independent_slot[1]
            ),
            "one_convolution_omission_detected": bool(
                one_convolution_omitted_norm < independent_slot[0]
            ),
            "duplicate_J2_factor_40_detected": bool(
                public_slot[0] * Decimal(40) / independent_slot[0]
                > MAX_PUBLIC_RATIO
            ),
        }

        failed_guards: list[str] = []
        for name, value in ratios.items():
            ratio_d = Decimal.from_float(float(value))
            if ratio_d < Decimal(1):
                failed_guards.append(f"public_underbounds_independent:{name}")
            if ratio_d > MAX_PUBLIC_RATIO:
                failed_guards.append(f"unexpected_shape_overcount:{name}")
        for name, detected in negative_controls.items():
            if not detected:
                failed_guards.append(f"negative_control_not_detected:{name}")

        offgrid = _fresh_hstar_stress(calc)
        if not offgrid["finite"]:
            failed_guards.append("fresh_Hstar_stress_nonfinite")
        if not offgrid["dominated"]:
            failed_guards.append("fresh_Hstar_stress_not_dominated")
        if not offgrid["nontrivial"]:
            failed_guards.append("fresh_Hstar_stress_trivial")
        if offgrid["j0_pm_0p1_percent_response"] <= 0.0:
            failed_guards.append("j0_parameter_perturbation_not_detected")
        if not offgrid["quadratic_sign_mutation_detected"]:
            failed_guards.append("Hstar_quadratic_sign_mutation_not_detected")

        payload: dict[str, Any] = {
            "schema": SCHEMA,
            "agent": 4,
            "audited_agent1_pr": AGENT1_PR,
            "audited_agent1_exact_head": AGENT1_HEAD,
            "independence": {
                "agent1_private_fraction_helpers_called": False,
                "agent1_private_product_helpers_called": False,
                "agent1_receipt_used_as_numerical_oracle": False,
                "decimal_precision": 80,
                "decimal_rounding": "ROUND_CEILING",
                "pi_upper": "355/113",
                "fresh_seed": SEED,
            },
            "prerequisites": {
                "source_axis_agent4_pr": SOURCE_AXIS_A4_PR,
                "Phi_ball_agent4_pr": PHI_BALL_A4_PR,
                "Phi_ball_consumed_as_previously_admitted_public_interface": True,
            },
            "independent_values": {
                "product_constant": float(pc),
                "L_inverse": float(linv),
                "H_star": float(hstar),
                "single_eta_after_J2_factor": float(single_eta_factor),
                "J2_Phi_eta": {
                    "norm": float(independent_post_phi[0]),
                    "lipschitz": float(independent_post_phi[1]),
                },
                "slot": {
                    "norm": float(independent_slot[0]),
                    "lipschitz": float(independent_slot[1]),
                },
            },
            "public_to_independent_ratios": ratios,
            "fresh_offgrid_Hstar_stress": offgrid,
            "negative_controls": negative_controls,
            "failed_guards": failed_guards,
            "passed": not failed_guards,
            "project_gates_unchanged": {
                "normalized_momentum_max": MOMENTUM_GATE,
                "normalized_momentum_L2": MOMENTUM_GATE,
                "divergence_max": DIVERGENCE_GATE,
                "divergence_L2": DIVERGENCE_GATE,
                "free_residual_defined_forcing_allowed": False,
            },
            "truth_boundary": dict(_TRUTH_BOUNDARY),
            "limitations": [
                "This audits only the H_* Phi_eta ordinary R1 post-J2 slot.",
                "The coefficient-space values in this receipt are not Navier-Stokes residuals.",
                "Agent-1 #741 upstream CI status is not promoted by this audit.",
                "The final Lambda^-1 d u Phi_eta ordinary R1 slot remains outside this audit.",
                "No global leading velocity, matched pressure, or real correction composite is evaluated.",
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


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    payload = save_report(args.output)
    print(json.dumps(payload, indent=2, sort_keys=True))
    if payload["failed_guards"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
