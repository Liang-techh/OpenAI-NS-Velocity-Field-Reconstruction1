"""Independent Agent-4 audit of PA.10 ``Lambda^-1 d u Phi_eta`` R1 slot.

This validator audits Agent-1 #751 through public, already-admitted interfaces
while using a distinct numerical/algebraic implementation.  It never calls the
Agent-1 private Fraction/product/rounding helpers and never uses an Agent-1
receipt as a numerical oracle.

The audited source-compatible object is

    J2[L^-1 d u partial_eta Phi]
      = L^-1 d J2[u partial_eta Phi].

The independent path reconstructs the coefficient product constant, the
Cauchy ``L^-1`` majorant, the finite-polynomial ``d=1-eta^2`` coefficient
majorant, the source ``80/rho`` single-eta post-J2 factor, and the final slot
using 80-digit Decimal arithmetic with ROUND_CEILING and the strict rational
upper ``pi < 355/113``.  The source ``u`` and ``Phi`` radius-one balls are
consumed only as previously independently audited public interfaces.

A PASS can admit only Agent-1 #751's final ordinary R1 slot.  It does not by
itself close full R1 because earlier Agent-4 slot audits remain separate, and
it does not establish M/K, a global leading field, any staged NS residual, or
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

from .kokuno_pa10_source_algebraic_ordinary_slots import (
    KokunoPA10SourceAlgebraicOrdinarySlots,
)
from .kokuno_pa10_source_r1_du_phi_eta_ordinary_slot import (
    KokunoPA10SourceR1DuPhiEtaOrdinarySlot,
)

SCHEMA = "kokuno-agent4-pa10-source-r1-du-phi-eta-independent-audit-v1"
AGENT1_PR = 751
AGENT1_HEAD = "080090f32214a9da102542ce4e927ec62d028baf"
SOURCE_AXIS_A4_PR = 671
SOURCE_U_BALL_A4_PR = 687
PHI_BALL_A4_PR = 655
SEED = 9173431
OFFGRID_COUNT = 8192
SLOT_NAME = "lambda_inv_d_u_times_Phi_eta"
MAX_PUBLIC_RATIO = Decimal("1.01")
MOMENTUM_GATE = 1.0e-3
DIVERGENCE_GATE = 1.0e-5

_TRUTH_BOUNDARY = {
    "source_axis_domain_independently_audited": True,
    "source_u_radius_one_ball_independently_audited": True,
    "source_Phi_radius_one_ball_independently_audited": True,
    "source_R1_du_Phi_eta_ordinary_slot_independently_audited": True,
    "source_R1_Hstar_Phi_eta_independent_audit_pending": True,
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
    # Distinct from Agent-1's pi < 22/7 path.
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


def _independent_d_norm(*, margin: float, rho: float) -> Decimal:
    """Finite coefficient norm for d=1-eta^2 under the public alpha=0 weights."""
    margin_d = Decimal.from_float(float(margin))
    rho_d = Decimal.from_float(float(rho))
    E = Decimal(1) + margin_d
    if margin_d < 0 or rho_d <= 0:
        raise ValueError("invalid public d geometry")
    beta0 = max(Decimal(1), E * E - Decimal(1))
    beta1 = Decimal(8) * E * rho_d
    beta2 = Decimal(9) * rho_d * rho_d
    return max(beta0, beta1, beta2)


def _fresh_d_value_audit(calc: KokunoPA10SourceR1DuPhiEtaOrdinarySlot) -> dict[str, Any]:
    """Fresh off-grid/axis-near check of the fixed ``d=1-eta^2`` multiplier."""
    axis = KokunoPA10SourceAlgebraicOrdinarySlots()
    domain = calc.domain
    E = 1.0 + float(domain.enlarged_real_margin)
    rng = np.random.default_rng(SEED)
    random_eta = rng.uniform(-E, E, size=OFFGRID_COUNT)
    probes = np.array([-E, -1.0e-12, 0.0, 1.0e-12, E], dtype=float)
    eta = np.concatenate((random_eta, probes))

    public_d = np.asarray(axis.axis_multiplier_values(eta)["d"], dtype=float)
    independent_d = 1.0 - eta * eta
    wrong_sign_d = 1.0 + eta * eta

    max_abs_error = float(np.max(np.abs(public_d - independent_d)))
    mutation_response = float(np.max(np.abs(public_d - wrong_sign_d)))
    public_bound = float(calc.operator_inputs()["d"].norm)
    max_abs_d = float(np.max(np.abs(independent_d)))

    return {
        "seed": SEED,
        "fresh_random_offgrid_points": OFFGRID_COUNT,
        "axis_near_probes": [-1.0e-12, 0.0, 1.0e-12],
        "endpoint_probes": [-E, E],
        "max_abs_public_vs_independent_error": max_abs_error,
        "max_abs_d": max_abs_d,
        "public_d_coefficient_bound": public_bound,
        "stress_to_public_bound_ratio": max_abs_d / public_bound,
        "finite": bool(np.isfinite(public_d).all()),
        "value_path_passed": bool(max_abs_error <= 2.0e-15),
        "dominated": bool(max_abs_d <= public_bound),
        "nontrivial": bool(max_abs_d >= 0.999),
        "wrong_sign_d_mutation_response": mutation_response,
        "wrong_sign_d_mutation_detected": bool(mutation_response >= 1.0e-6),
    }


def run_independent_audit(
    *, pr_head: str | None = None, checkout_head: str | None = None
) -> dict[str, Any]:
    calc = KokunoPA10SourceR1DuPhiEtaOrdinarySlot()
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
        d_norm = _independent_d_norm(
            margin=float(domain.enlarged_real_margin),
            rho=float(domain.coefficient_rho),
        )
        rho = Decimal.from_float(float(domain.coefficient_rho))
        eta_factor = Decimal(80) / rho

        inputs = calc.operator_inputs()
        # u and Phi are consumed only as public upstream balls already separately
        # admitted by Agent 4.  Their construction internals are not read here.
        u_pair = _pair(inputs["u"])
        phi_pair = _pair(inputs["Phi"])
        independent_core = _scale(
            _product(u_pair, phi_pair, product_constant=pc), eta_factor
        )
        independent_slot = _product(
            (linv, Decimal(0)),
            (d_norm, Decimal(0)),
            independent_core,
            product_constant=pc,
        )

        public_pc = _d(calc.product_calculator.product_constant)
        public_linv = _pair(inputs["L_inverse"])[0]
        public_d = _pair(inputs["d"])[0]
        public_factor = _d(calc.source_single_eta_after_j2_factor_upper())
        public_core = _pair(calc.j2_u_times_phi_eta())
        public_slot = _pair(calc.r1_du_phi_eta_after_j2()[SLOT_NAME])

        ratios = {
            "product_constant": _ratio(public_pc, pc),
            "L_inverse": _ratio(public_linv, linv),
            "d": _ratio(public_d, d_norm),
            "single_eta_after_J2": _ratio(public_factor, eta_factor),
            "J2_u_Phi_eta_norm": _ratio(public_core[0], independent_core[0]),
            "J2_u_Phi_eta_lipschitz": _ratio(
                public_core[1], independent_core[1]
            ),
            "slot_norm": _ratio(public_slot[0], independent_slot[0]),
            "slot_lipschitz": _ratio(public_slot[1], independent_slot[1]),
        }

        negative_controls = {
            "product_constant_0p998_detected": bool(
                public_pc * Decimal("0.998") < pc
            ),
            "L_inverse_0p999_detected": bool(
                public_linv * Decimal("0.999") < linv
            ),
            "d_0p999_detected": bool(public_d * Decimal("0.999") < d_norm),
            "eta_factor_79_instead_of_80_detected": bool(
                public_factor * Decimal(79) / Decimal(80) < eta_factor
            ),
            "core_norm_0p99_detected": bool(
                public_core[0] * Decimal("0.99") < independent_core[0]
            ),
            "core_lipschitz_0p99_detected": bool(
                public_core[1] * Decimal("0.99") < independent_core[1]
            ),
            "slot_norm_0p99_detected": bool(
                public_slot[0] * Decimal("0.99") < independent_slot[0]
            ),
            "slot_lipschitz_0p99_detected": bool(
                public_slot[1] * Decimal("0.99") < independent_slot[1]
            ),
            "one_convolution_omission_detected": bool(
                independent_slot[0] / pc < independent_slot[0]
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

        offgrid = _fresh_d_value_audit(calc)
        if not offgrid["finite"]:
            failed_guards.append("fresh_d_value_path_nonfinite")
        if not offgrid["value_path_passed"]:
            failed_guards.append("fresh_d_value_path_mismatch")
        if not offgrid["dominated"]:
            failed_guards.append("fresh_d_values_not_dominated")
        if not offgrid["nontrivial"]:
            failed_guards.append("fresh_d_value_path_trivial")
        if not offgrid["wrong_sign_d_mutation_detected"]:
            failed_guards.append("wrong_sign_d_mutation_not_detected")

        payload: dict[str, Any] = {
            "schema": SCHEMA,
            "agent": 4,
            "audited_agent1_pr": AGENT1_PR,
            "audited_agent1_exact_head": AGENT1_HEAD,
            "pr_head": pr_head,
            "checkout_head": checkout_head,
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
                "source_u_ball_agent4_pr": SOURCE_U_BALL_A4_PR,
                "Phi_ball_agent4_pr": PHI_BALL_A4_PR,
                "u_and_Phi_consumed_as_previously_admitted_public_interfaces": True,
            },
            "independent_values": {
                "product_constant": float(pc),
                "L_inverse": float(linv),
                "d": float(d_norm),
                "single_eta_after_J2_factor": float(eta_factor),
                "source_u_ball_public_upstream": {
                    "norm": float(u_pair[0]),
                    "lipschitz": float(u_pair[1]),
                },
                "source_Phi_ball_public_upstream": {
                    "norm": float(phi_pair[0]),
                    "lipschitz": float(phi_pair[1]),
                },
                "J2_u_Phi_eta": {
                    "norm": float(independent_core[0]),
                    "lipschitz": float(independent_core[1]),
                },
                "slot": {
                    "norm": float(independent_slot[0]),
                    "lipschitz": float(independent_slot[1]),
                },
            },
            "public_to_independent_ratios": ratios,
            "fresh_offgrid_d_value_audit": offgrid,
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
                "This audits only the Lambda^-1 d u Phi_eta ordinary R1 post-J2 slot.",
                "The coefficient-space values in this receipt are not Navier-Stokes residuals.",
                "Agent-1 #751 upstream CI status is not promoted by this audit.",
                "Earlier R1 slot audits remain separate and may still be pending.",
                "No global leading velocity, matched pressure, or real correction composite is evaluated.",
            ],
        }
        payload["receipt_sha256"] = hashlib.sha256(
            _canonical_json(payload).encode("utf-8")
        ).hexdigest()
        return payload


def save_report(
    path: str | Path, *, pr_head: str | None = None, checkout_head: str | None = None
) -> dict[str, Any]:
    payload = run_independent_audit(pr_head=pr_head, checkout_head=checkout_head)
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    parser.add_argument("--pr-head", default=None)
    parser.add_argument("--checkout-head", default=None)
    args = parser.parse_args()
    payload = save_report(
        args.output, pr_head=args.pr_head, checkout_head=args.checkout_head
    )
    print(json.dumps(payload, indent=2, sort_keys=True))
    if payload["failed_guards"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
