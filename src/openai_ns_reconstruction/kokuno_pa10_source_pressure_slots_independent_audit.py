"""Independent Agent-4 audit of PA.10 source pressure ordinary R2 slots.

This validator is intentionally implementation-distinct from Agent-1 #676.
It consumes only the public source-pressure slot API and its public upstream
bounds, then reconstructs the three pressure monomials and post-J1 factors with
80-digit Decimal arithmetic and the independent rational upper pi < 355/113.
It does not call Agent-1 private product/rounding helpers and does not use the
Agent-1 receipt as an oracle.

Scientific scope is narrow: this can admit only the three pressure ordinary R2
slots.  It is not a full R2/R1, leading field, pressure match, correction, or
Navier-Stokes residual validation.
"""

from __future__ import annotations

import argparse
import copy
from decimal import Decimal, ROUND_CEILING, localcontext
import hashlib
import json
import math
from pathlib import Path
from typing import Any

from .kokuno_pa10_source_pressure_ordinary_slots import (
    KokunoPA10SourcePressureOrdinarySlots,
)

SCHEMA = "kokuno-agent4-pa10-source-pressure-slots-independent-audit-v1"
AGENT1_PR = 676
AGENT1_HEAD = "6637d6f6ed4b80c6531d4db1f0b415e2c6c77edc"
SOURCE_AXIS_A4_PR = 671
SOURCE_PHI_A4_PR = 655
MIXED_J_A4_PR = 663
MUTATION_FACTOR_MULTIPLIER = Decimal("0.999")
MUTATION_FACTOR_SLOT = Decimal("0.99")

_TRUTH_BOUNDARY = {
    "source_axis_domain_independently_audited": True,
    "source_Phi_ball_conditionally_independently_audited": True,
    "conditional_mixed_J_algebra_independently_audited": True,
    "source_R2_pressure_ordinary_slots_independently_audited": True,
    "all_R2_ordinary_slots_source_bound": False,
    "all_R1_ordinary_slots_source_bound": False,
    "source_full_post_J1_R2_radius_one_ball_bound_machine_bound": False,
    "source_full_post_J2_R1_radius_one_ball_bound_machine_bound": False,
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
}


def _canonical_json(payload: dict[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _d(value: float) -> Decimal:
    out = float(value)
    if not math.isfinite(out) or out < 0.0:
        raise ValueError("all public bounds must be finite and nonnegative")
    return Decimal.from_float(out)


def _pair(bound: Any) -> tuple[Decimal, Decimal]:
    return _d(bound.norm), _d(bound.lipschitz)


def _ratio(public: Decimal, independent: Decimal) -> float:
    if independent <= 0:
        raise ValueError("independent comparator must be positive")
    return float(public / independent)


def _scale(pair: tuple[Decimal, Decimal], factor: Decimal) -> tuple[Decimal, Decimal]:
    if factor < 0:
        factor = -factor
    return pair[0] * factor, pair[1] * factor


def _product2(
    left: tuple[Decimal, Decimal],
    right: tuple[Decimal, Decimal],
    product_constant: Decimal,
) -> tuple[Decimal, Decimal]:
    return (
        product_constant * left[0] * right[0],
        product_constant * (left[1] * right[0] + right[1] * left[0]),
    )


def _independent_product_constant() -> Decimal:
    # Classical strict rational upper pi < 355/113, distinct from Agent-1 22/7.
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
        raise ValueError("public source-compatible Cauchy geometry is invalid")
    x = rho_d / radius_d
    weight_sum = (Decimal(1) + x) / (Decimal(1) - x) ** 3
    return weight_sum / lower


def run_independent_audit(*, pr_head: str | None = None, checkout_head: str | None = None) -> dict[str, Any]:
    upstream = KokunoPA10SourcePressureOrdinarySlots()
    domain = upstream.domain

    with localcontext() as ctx:
        ctx.prec = 80
        ctx.rounding = ROUND_CEILING

        product_constant = _independent_product_constant()
        rho = Decimal.from_float(float(domain.coefficient_rho))
        margin = Decimal.from_float(float(domain.enlarged_real_margin))
        h = Decimal.from_float(float(domain.h))
        eta_ind = max(Decimal(1) + margin, Decimal(4) * rho)
        linv_ind = _independent_l_inverse_norm(
            h=float(domain.h),
            margin=float(domain.enlarged_real_margin),
            radius=float(domain.cauchy_radius),
            rho=float(domain.coefficient_rho),
        )

        public_eta = _pair(upstream.eta_coefficient_ball())
        public_linv = _pair(upstream.bridge.l_inverse_coefficient_ball())
        pressure_public = upstream.source_pressure_ball()
        numerator_public = upstream.r2_pressure_numerator_slots()
        after_public = upstream.r2_pressure_after_j1()

        eta_pair = (eta_ind, Decimal(0))
        # The fixed d=1-eta^2 coefficient-ball norm 1.0 was already separately
        # cross-audited in A4 #663; do not reuse Agent-1's private helper here.
        d_pair = (Decimal(1), Decimal(0))
        p_pair = _pair(pressure_public["p"])
        p_eta_pair = _pair(pressure_public["p_eta"])
        yp_pair = _pair(pressure_public["Y_p_Y"])

        A = Decimal("0.5") + h
        num_ind: dict[str, tuple[Decimal, Decimal]] = {}
        num_ind["4A_eta_p"] = _scale(
            _product2(eta_pair, p_pair, product_constant), Decimal(4) * A
        )
        num_ind["d_p_eta"] = _product2(d_pair, p_eta_pair, product_constant)
        num_ind["2eta_Y_p_Y"] = _scale(
            _product2(eta_pair, yp_pair, product_constant), Decimal(2)
        )

        # For ordinary vanishing order b=0 and nu=1, J1 factor is exactly 80.
        post_factor = product_constant * Decimal(80) * linv_ind
        after_ind = {
            name: _scale(pair, post_factor) for name, pair in num_ind.items()
        }

        multiplier_ratios = {
            "eta_norm": _ratio(public_eta[0], eta_ind),
            "L_inverse_norm": _ratio(public_linv[0], linv_ind),
        }
        numerator_ratios: dict[str, dict[str, float]] = {}
        after_ratios: dict[str, dict[str, float]] = {}
        for name in num_ind:
            pub_num = _pair(numerator_public[name])
            pub_after = _pair(after_public[name])
            numerator_ratios[name] = {
                "norm": _ratio(pub_num[0], num_ind[name][0]),
                "lipschitz": _ratio(pub_num[1], num_ind[name][1]),
            }
            after_ratios[name] = {
                "norm": _ratio(pub_after[0], after_ind[name][0]),
                "lipschitz": _ratio(pub_after[1], after_ind[name][1]),
            }

        failed_guards: list[str] = []
        if any(value < 1.0 for value in multiplier_ratios.values()):
            failed_guards.append("public multiplier underbounds independent oracle")
        for name, ratios in numerator_ratios.items():
            if min(ratios.values()) < 1.0:
                failed_guards.append(f"{name} numerator underbounds independent oracle")
        for name, ratios in after_ratios.items():
            if min(ratios.values()) < 1.0:
                failed_guards.append(f"{name} post-J1 underbounds independent oracle")

        mutation = {
            "eta_0p999_detected": public_eta[0] * MUTATION_FACTOR_MULTIPLIER < eta_ind,
            "L_inverse_0p999_detected": public_linv[0] * MUTATION_FACTOR_MULTIPLIER < linv_ind,
            "numerator_0p99_detected": {},
            "post_J1_0p99_detected": {},
        }
        for name in num_ind:
            pub_num = _pair(numerator_public[name])
            pub_after = _pair(after_public[name])
            mutation["numerator_0p99_detected"][name] = bool(
                pub_num[0] * MUTATION_FACTOR_SLOT < num_ind[name][0]
                and pub_num[1] * MUTATION_FACTOR_SLOT < num_ind[name][1]
            )
            mutation["post_J1_0p99_detected"][name] = bool(
                pub_after[0] * MUTATION_FACTOR_SLOT < after_ind[name][0]
                and pub_after[1] * MUTATION_FACTOR_SLOT < after_ind[name][1]
            )
        if not mutation["eta_0p999_detected"]:
            failed_guards.append("eta mutation not detected")
        if not mutation["L_inverse_0p999_detected"]:
            failed_guards.append("L inverse mutation not detected")
        if not all(mutation["numerator_0p99_detected"].values()):
            failed_guards.append("one or more numerator mutations not detected")
        if not all(mutation["post_J1_0p99_detected"].values()):
            failed_guards.append("one or more post-J1 mutations not detected")

        independent_values = {
            "product_constant_pi355_over_113": float(product_constant),
            "eta_norm": float(eta_ind),
            "L_inverse_norm": float(linv_ind),
            "R2_pressure_ordinary_numerator_slots": {
                name: {"norm": float(pair[0]), "lipschitz": float(pair[1])}
                for name, pair in num_ind.items()
            },
            "R2_pressure_ordinary_after_J1": {
                name: {"norm": float(pair[0]), "lipschitz": float(pair[1])}
                for name, pair in after_ind.items()
            },
        }

    payload: dict[str, Any] = {
        "schema": SCHEMA,
        "agent1_pr": AGENT1_PR,
        "agent1_head": AGENT1_HEAD,
        "pr_head": pr_head,
        "checkout_head": checkout_head,
        "independent_method": {
            "arithmetic": "80-digit Decimal with ROUND_CEILING",
            "pi_upper": "355/113",
            "agent1_private_helpers_called": False,
            "agent1_receipt_used_as_oracle": False,
            "pressure_ball_treated_as_previously_audited_public_upstream_input": True,
            "d_ball_source": "independent A4 #663 admitted fixed d=1-eta^2 coefficient norm 1",
        },
        "upstream_independent_admissions": {
            "source_axis_A4_PR": SOURCE_AXIS_A4_PR,
            "source_Phi_A4_PR": SOURCE_PHI_A4_PR,
            "mixed_J_A4_PR": MIXED_J_A4_PR,
        },
        "independent_values": independent_values,
        "public_to_independent_ratios": {
            "multipliers": multiplier_ratios,
            "numerator_slots": numerator_ratios,
            "post_J1_slots": after_ratios,
        },
        "mutation": mutation,
        "failed_guards": failed_guards,
        "source_R2_pressure_ordinary_slots_independent_preflight_passed": not failed_guards,
        "final_project_gates_unchanged": {
            "normalized_momentum_max": 1e-3,
            "normalized_momentum_L2": 1e-3,
            "divergence_max": 1e-5,
            "divergence_L2": 1e-5,
        },
        "truth_boundary": copy.deepcopy(_TRUTH_BOUNDARY),
    }
    payload["receipt_sha256"] = hashlib.sha256(
        _canonical_json(payload).encode("utf-8")
    ).hexdigest()
    return payload


def save_report(path: str | Path, *, pr_head: str | None = None, checkout_head: str | None = None) -> dict[str, Any]:
    payload = run_independent_audit(pr_head=pr_head, checkout_head=checkout_head)
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return payload


def _main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--pr-head", type=str, default=None)
    parser.add_argument("--checkout-head", type=str, default=None)
    args = parser.parse_args()
    print(
        json.dumps(
            save_report(
                args.output,
                pr_head=args.pr_head,
                checkout_head=args.checkout_head,
            ),
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    _main()
