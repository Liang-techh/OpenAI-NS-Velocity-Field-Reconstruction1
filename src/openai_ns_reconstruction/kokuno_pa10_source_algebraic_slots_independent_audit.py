"""Independent Agent-4 audit of PA.10 algebraic ordinary R2 slots.

This validator is intentionally implementation-distinct from Agent-1 #684.
It consumes only the public source-algebraic slot API plus already-admitted
upstream public bounds, then reconstructs the four algebraic monomials and
post-J1 factors with 80-digit Decimal arithmetic and the independent rational
upper pi < 355/113. It does not call Agent-1 private Fraction/product/rounding
helpers and does not use the Agent-1 receipt as an oracle.

Scientific scope is narrow: this can admit only #684's four algebraic ordinary
R2 slots, conditional on the already independently admitted source-axis and
source-u-ball prerequisites. It is not a full R2/R1, leading-field, pressure,
correction, or Navier-Stokes residual validation.
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

import numpy as np

from .kokuno_pa10_source_algebraic_ordinary_slots import (
    KokunoPA10SourceAlgebraicOrdinarySlots,
)

SCHEMA = "kokuno-agent4-pa10-source-algebraic-slots-independent-audit-v1"
AGENT1_PR = 684
AGENT1_HEAD = "e3db3d68363ddc5b666a5111a86634ea360b3bdd"
SOURCE_AXIS_A4_PR = 671
PRESSURE_SLOTS_A4_PR = 678
SOURCE_U_BALL_A4_PR = 687
MIXED_J_A4_PR = 663
SEED = 9173361
OFFGRID_COUNT = 4096
MUTATION_FACTOR_MULTIPLIER = Decimal("0.999")
MUTATION_FACTOR_SLOT = Decimal("0.99")

_TRUTH_BOUNDARY = {
    "source_axis_domain_independently_audited": True,
    "source_R2_pressure_ordinary_slots_independently_audited": True,
    "source_u_radius_one_ball_independently_audited": True,
    "conditional_mixed_J_algebra_independently_audited": True,
    "source_R2_algebraic_ordinary_slots_independently_audited": True,
    "all_R2_ordinary_slots_source_bound": False,
    "all_R2_ordinary_slots_independently_audited": False,
    "all_R1_ordinary_slots_source_bound": False,
    "source_full_post_J1_R2_radius_one_ball_bound_machine_bound": False,
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


def _fresh_axis_value_audit(upstream: KokunoPA10SourceAlgebraicOrdinarySlots) -> dict[str, Any]:
    """Fresh off-grid/axis-near replay of the public fixed multiplier API."""
    domain = upstream.domain
    E = 1.0 + float(domain.enlarged_real_margin)
    rng = np.random.default_rng(SEED)
    random_eta = rng.uniform(-E, E, size=OFFGRID_COUNT)
    probes = np.array([-E, -1.0e-12, 0.0, 1.0e-12, E], dtype=float)
    eta = np.concatenate((random_eta, probes))
    public = upstream.axis_multiplier_values(eta)
    h = float(domain.h)
    j0 = float(domain.j0)
    independent = {
        "eta": eta,
        "A": np.full_like(eta, 0.5 + h),
        "U_star": 4.0 * eta + j0,
        "U_star_eta": np.full_like(eta, 4.0),
        "d": 1.0 - eta * eta,
    }
    max_abs_error = {
        name: float(np.max(np.abs(np.asarray(public[name]) - values)))
        for name, values in independent.items()
    }
    mutated_ustar = 3.99 * eta + j0
    mutation_response = float(np.max(np.abs(np.asarray(public["U_star"]) - mutated_ustar)))
    return {
        "seed": SEED,
        "fresh_random_offgrid_points": OFFGRID_COUNT,
        "axis_near_probes": [-1.0e-12, 0.0, 1.0e-12],
        "endpoint_probes": [-E, E],
        "max_abs_error": max_abs_error,
        "U_star_slope_4_to_3p99_mutation_response": mutation_response,
        "value_path_passed": max(max_abs_error.values()) <= 2.0e-15,
        "value_path_mutation_detected": mutation_response >= 5.0e-3,
    }


def run_independent_audit(
    *, pr_head: str | None = None, checkout_head: str | None = None
) -> dict[str, Any]:
    upstream = KokunoPA10SourceAlgebraicOrdinarySlots()
    domain = upstream.domain

    with localcontext() as ctx:
        ctx.prec = 80
        ctx.rounding = ROUND_CEILING

        product_constant = _independent_product_constant()
        rho = Decimal.from_float(float(domain.coefficient_rho))
        margin = Decimal.from_float(float(domain.enlarged_real_margin))
        h = Decimal.from_float(float(domain.h))
        j0 = Decimal.from_float(float(domain.j0))
        E = Decimal(1) + margin
        A_ind = Decimal("0.5") + h
        eta_ind = max(E, Decimal(4) * rho)
        ustar_ind = max(Decimal(4) * E + j0, Decimal(16) * rho)
        linv_ind = _independent_l_inverse_norm(
            h=float(domain.h),
            margin=float(domain.enlarged_real_margin),
            radius=float(domain.cauchy_radius),
            rho=float(domain.coefficient_rho),
        )

        public_fixed = upstream.fixed_multiplier_balls()
        public_linv = _pair(upstream.bridge.l_inverse_coefficient_ball())
        public_num = upstream.r2_algebraic_numerator_slots()
        public_after = upstream.r2_algebraic_after_j1()

        # #687 independently admitted this source-compatible u ball. Treat it as
        # a frozen public upstream input rather than recomputing Agent-1 internals.
        u_pair = _pair(upstream.bridge.u_source.u_ball())
        A_pair = (A_ind, Decimal(0))
        eta_pair = (eta_ind, Decimal(0))
        ustar_pair = (ustar_ind, Decimal(0))
        d_pair = (Decimal(1), Decimal(0))

        num_ind: dict[str, tuple[Decimal, Decimal]] = {
            "A_times_u": _scale(u_pair, A_ind),
            "4A_eta_Ustar_times_u": _scale(
                _product(
                    eta_pair,
                    ustar_pair,
                    u_pair,
                    product_constant=product_constant,
                ),
                Decimal(4) * A_ind,
            ),
            "d_Ustar_eta_times_u": _scale(
                _product(d_pair, u_pair, product_constant=product_constant),
                Decimal(4),
            ),
            "2A_lambda_inv_eta_u_squared": _scale(
                _product(
                    eta_pair,
                    u_pair,
                    u_pair,
                    product_constant=product_constant,
                ),
                Decimal(2) * A_ind,
            ),
        }

        # Ordinary b=0, nu=1: J1 factor is exactly 80. Multiplying by L^-1
        # costs one additional product-algebra factor. Lambda^-1<=1 is already
        # absorbed in the fourth numerator's fail-safe source-regime bound.
        post_factor = product_constant * Decimal(80) * linv_ind
        after_ind = {name: _scale(pair, post_factor) for name, pair in num_ind.items()}

        multiplier_public = {
            "A": _pair(public_fixed["A"])[0],
            "eta": _pair(public_fixed["eta"])[0],
            "U_star": _pair(public_fixed["U_star"])[0],
            "U_star_eta": _pair(public_fixed["U_star_eta"])[0],
            "d": _pair(public_fixed["d"])[0],
            "L_inverse": public_linv[0],
        }
        multiplier_ind = {
            "A": A_ind,
            "eta": eta_ind,
            "U_star": ustar_ind,
            "U_star_eta": Decimal(4),
            "d": Decimal(1),
            "L_inverse": linv_ind,
        }
        multiplier_ratios = {
            name: _ratio(multiplier_public[name], multiplier_ind[name])
            for name in multiplier_ind
        }

        numerator_ratios: dict[str, dict[str, float]] = {}
        after_ratios: dict[str, dict[str, float]] = {}
        for name in num_ind:
            pub_num = _pair(public_num[name])
            pub_after = _pair(public_after[name])
            numerator_ratios[name] = {
                "norm": _ratio(pub_num[0], num_ind[name][0]),
                "lipschitz": _ratio(pub_num[1], num_ind[name][1]),
            }
            after_ratios[name] = {
                "norm": _ratio(pub_after[0], after_ind[name][0]),
                "lipschitz": _ratio(pub_after[1], after_ind[name][1]),
            }

        value_audit = _fresh_axis_value_audit(upstream)
        failed_guards: list[str] = []
        if any(value < 1.0 for value in multiplier_ratios.values()):
            failed_guards.append("public fixed multiplier underbounds independent oracle")
        for name, ratios in numerator_ratios.items():
            if min(ratios.values()) < 1.0:
                failed_guards.append(f"{name} numerator underbounds independent oracle")
        for name, ratios in after_ratios.items():
            if min(ratios.values()) < 1.0:
                failed_guards.append(f"{name} post-J1 underbounds independent oracle")
        if not value_audit["value_path_passed"]:
            failed_guards.append("fresh off-grid fixed-axis public value replay failed")
        if not value_audit["value_path_mutation_detected"]:
            failed_guards.append("fresh off-grid U_star slope mutation not detected")

        mutation = {
            "multiplier_0p999_detected": {
                name: bool(
                    multiplier_public[name] * MUTATION_FACTOR_MULTIPLIER
                    < multiplier_ind[name]
                )
                for name in ("A", "eta", "U_star", "L_inverse")
            },
            "numerator_0p99_detected": {},
            "post_J1_0p99_detected": {},
            "axis_value_slope_mutation_detected": bool(
                value_audit["value_path_mutation_detected"]
            ),
        }
        for name in num_ind:
            pub_num = _pair(public_num[name])
            pub_after = _pair(public_after[name])
            mutation["numerator_0p99_detected"][name] = bool(
                pub_num[0] * MUTATION_FACTOR_SLOT < num_ind[name][0]
                and pub_num[1] * MUTATION_FACTOR_SLOT < num_ind[name][1]
            )
            mutation["post_J1_0p99_detected"][name] = bool(
                pub_after[0] * MUTATION_FACTOR_SLOT < after_ind[name][0]
                and pub_after[1] * MUTATION_FACTOR_SLOT < after_ind[name][1]
            )
        if not all(mutation["multiplier_0p999_detected"].values()):
            failed_guards.append("one or more fixed-multiplier mutations not detected")
        if not all(mutation["numerator_0p99_detected"].values()):
            failed_guards.append("one or more numerator mutations not detected")
        if not all(mutation["post_J1_0p99_detected"].values()):
            failed_guards.append("one or more post-J1 mutations not detected")

        independent_values = {
            "product_constant_pi355_over_113": float(product_constant),
            "fixed_multiplier_norms": {
                name: float(value) for name, value in multiplier_ind.items()
            },
            "source_u_ball_public_upstream": {
                "norm": float(u_pair[0]),
                "lipschitz": float(u_pair[1]),
                "independent_admission_PR": SOURCE_U_BALL_A4_PR,
            },
            "R2_algebraic_ordinary_numerator_slots": {
                name: {"norm": float(pair[0]), "lipschitz": float(pair[1])}
                for name, pair in num_ind.items()
            },
            "R2_algebraic_ordinary_after_J1": {
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
            "source_u_ball_treated_as_independently_admitted_public_input": True,
            "fresh_offgrid_axis_seed": SEED,
            "fresh_offgrid_axis_random_count": OFFGRID_COUNT,
        },
        "upstream_independent_admissions": {
            "source_axis_A4_PR": SOURCE_AXIS_A4_PR,
            "pressure_slots_A4_PR": PRESSURE_SLOTS_A4_PR,
            "source_u_ball_A4_PR": SOURCE_U_BALL_A4_PR,
            "mixed_J_A4_PR": MIXED_J_A4_PR,
        },
        "independent_values": independent_values,
        "fresh_axis_value_audit": value_audit,
        "public_to_independent_ratios": {
            "multipliers": multiplier_ratios,
            "numerator_slots": numerator_ratios,
            "post_J1_slots": after_ratios,
        },
        "mutation": mutation,
        "failed_guards": failed_guards,
        "source_R2_algebraic_ordinary_slots_independent_preflight_passed": not failed_guards,
        "latest_agent1_derivative_R2_self_certificate_PR": 693,
        "latest_agent1_derivative_R2_subset_independently_audited_here": False,
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


def save_report(
    path: str | Path, *, pr_head: str | None = None, checkout_head: str | None = None
) -> dict[str, Any]:
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
