"""Independent Agent-4 audit of the PA.10 averaged-eta ordinary R1 slot.

This validator is deliberately implementation-distinct from Agent-1 #715.  It
consumes only the public source-compatible ball/value interfaces, reconstructs

    J_2[L^-1 d (partial_eta A(u)) Phi]

with 80-digit Decimal arithmetic, and uses the stricter rational upper
``pi < 355/113`` rather than Agent-1's ``22/7`` product-algebra path.  It never
calls Agent-1 private Fraction/product/rounding helpers and never uses an
Agent-1 receipt as a numerical oracle.

The source rule already includes the radial inverse and gives the averaged,
no-logradial specialization with factor ``80/rho``.  This audit therefore also
checks that the public bound has exactly three undifferentiated convolutions and
has not silently re-applied ``J_2``.  Scope is intentionally narrow: a PASS can
admit only Agent-1 #715's new averaged-eta post-J2 R1 slot.  It cannot promote
prior R1 slots, full R1/R2, M/K, global leading/pressure, any staged NS residual,
or ``pde_validated``.
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

from .kokuno_pa10_source_r1_averaged_eta_ordinary_slot import (
    KokunoPA10SourceR1AveragedEtaOrdinarySlot,
)

SCHEMA = "kokuno-agent4-pa10-source-r1-averaged-eta-independent-audit-v1"
AGENT1_PR = 715
AGENT1_HEAD = "e35f60335dcbb53a0227e5e37acb76ad21ae0dd0"
SOURCE_AXIS_A4_PR = 671
SOURCE_PHI_BALL_A4_PR = 655
SOURCE_U_BALL_A4_PR = 687
MIXED_J_A4_PR = 663
ZETA_R1_A4_PR = 710
SEED = 9173391
OFFGRID_COUNT = 8192
SLOT_NAME = "lambda_inv_d_detaAu_times_Phi"
MAX_EXPECTED_PUBLIC_RATIO = Decimal("1.01")

_TRUTH_BOUNDARY = {
    "source_axis_domain_independently_audited": True,
    "source_Phi_radius_one_ball_independently_audited": True,
    "source_u_radius_one_ball_independently_audited": True,
    "conditional_mixed_J_algebra_independently_audited": True,
    "source_R1_averaged_eta_ordinary_slot_independently_audited": True,
    "source_R1_zeta_ordinary_slot_independent_audit_pending": True,
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
    # Strict classical rational upper, independent from Agent-1's pi < 22/7.
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


def _independent_slot(
    calc: KokunoPA10SourceR1AveragedEtaOrdinarySlot,
) -> tuple[tuple[Decimal, Decimal], dict[str, Decimal]]:
    domain = calc.domain
    inputs = calc.operator_inputs()
    rho = Decimal.from_float(float(domain.coefficient_rho))
    product_constant = _independent_product_constant()
    derivative_after_j2 = Decimal(80) / rho
    linv = _independent_l_inverse_norm(
        h=float(domain.h),
        margin=float(domain.enlarged_real_margin),
        radius=float(domain.cauchy_radius),
        rho=float(domain.coefficient_rho),
    )

    d_norm = Decimal(1)
    u_norm, u_lip = _pair(inputs["u"])
    phi_norm, phi_lip = _pair(inputs["Phi"])
    fixed = product_constant**3 * derivative_after_j2 * d_norm * linv
    norm = fixed * u_norm * phi_norm
    lip = fixed * (u_lip * phi_norm + u_norm * phi_lip)
    return (norm, lip), {
        "rho": rho,
        "product_constant": product_constant,
        "derivative_after_j2": derivative_after_j2,
        "d_norm": d_norm,
        "L_inverse": linv,
        "fixed_prefactor": fixed,
    }


def _fresh_fixed_multiplier_stress(
    calc: KokunoPA10SourceR1AveragedEtaOrdinarySlot,
) -> dict[str, Any]:
    domain = calc.domain
    inputs = calc.operator_inputs()
    extent = 1.0 + float(domain.enlarged_real_margin)
    rng = np.random.default_rng(SEED)
    random_eta = rng.uniform(-extent, extent, size=OFFGRID_COUNT)
    probes = np.array([-extent, -1.0e-12, 0.0, 1.0e-12, extent], dtype=float)
    eta = np.concatenate((random_eta, probes))

    h = float(domain.h)
    d_value = 1.0 - eta * eta
    L_value = 1.0 - 2.0 * h * eta * eta
    Linv_value = 1.0 / L_value
    d_public = float(inputs["d"].norm)
    Linv_public = float(inputs["L_inverse"].norm)
    max_d = float(np.max(np.abs(d_value)))
    max_linv = float(np.max(np.abs(Linv_value)))

    return {
        "seed": SEED,
        "fresh_random_offgrid_points": OFFGRID_COUNT,
        "axis_near_probes": [-1.0e-12, 0.0, 1.0e-12],
        "endpoint_probes": [-extent, extent],
        "max_abs_d": max_d,
        "max_abs_L_inverse": max_linv,
        "public_d_coefficient_bound": d_public,
        "public_L_inverse_coefficient_bound": Linv_public,
        "d_value_over_public_bound": max_d / d_public,
        "L_inverse_value_over_public_bound": max_linv / Linv_public,
        "finite": bool(np.isfinite(d_value).all() and np.isfinite(Linv_value).all()),
        "dominated": bool(max_d <= d_public and max_linv <= Linv_public),
        "axis_nontrivial": bool(abs(float(Linv_value[-3])) > 0.0),
    }


def _rho_local_operator_sensitivity(
    calc: KokunoPA10SourceR1AveragedEtaOrdinarySlot,
) -> dict[str, float]:
    """Record a local operator-only rho perturbation; this is not a source re-fit."""
    domain = calc.domain
    h = float(domain.h)
    margin = float(domain.enlarged_real_margin)
    radius = float(domain.cauchy_radius)
    rho0 = float(domain.coefficient_rho)
    pc = _independent_product_constant()

    def prefactor(rho: float) -> float:
        linv = _independent_l_inverse_norm(
            h=h, margin=margin, radius=radius, rho=rho
        )
        return float(pc**3 * (Decimal(80) / Decimal.from_float(rho)) * linv)

    lo = rho0 * 0.999
    hi = rho0 * 1.001
    return {
        "rho_minus_0p1pct": lo,
        "prefactor_minus_0p1pct": prefactor(lo),
        "rho_nominal": rho0,
        "prefactor_nominal": prefactor(rho0),
        "rho_plus_0p1pct": hi,
        "prefactor_plus_0p1pct": prefactor(hi),
    }


def run_independent_audit() -> dict[str, Any]:
    calc = KokunoPA10SourceR1AveragedEtaOrdinarySlot()
    inputs = calc.operator_inputs()

    with localcontext() as ctx:
        ctx.prec = 80
        ctx.rounding = ROUND_CEILING

        (ind_norm, ind_lip), pieces = _independent_slot(calc)
        public_slot = _pair(calc.r1_averaged_eta_after_j2()[SLOT_NAME])
        public_product_constant = _d(calc.product_constant_upper)
        public_derivative = _d(calc.source_single_eta_after_j2_factor_upper())
        public_d = _pair(inputs["d"])[0]
        public_linv = _pair(inputs["L_inverse"])[0]

        ratios = {
            "product_constant": _ratio(
                public_product_constant, pieces["product_constant"]
            ),
            "derivative_after_j2": _ratio(
                public_derivative, pieces["derivative_after_j2"]
            ),
            "d_norm": _ratio(public_d, pieces["d_norm"]),
            "L_inverse": _ratio(public_linv, pieces["L_inverse"]),
            "slot_norm": _ratio(public_slot[0], ind_norm),
            "slot_lipschitz": _ratio(public_slot[1], ind_lip),
        }

        negative_controls = {
            "product_constant_0p998_detected": bool(
                public_product_constant * Decimal("0.998")
                < pieces["product_constant"]
            ),
            "derivative_79_over_80_detected": bool(
                public_derivative * Decimal(79) / Decimal(80)
                < pieces["derivative_after_j2"]
            ),
            "L_inverse_0p999_detected": bool(
                public_linv * Decimal("0.999") < pieces["L_inverse"]
            ),
            "slot_norm_0p99_detected": bool(
                public_slot[0] * Decimal("0.99") < ind_norm
            ),
            "slot_lipschitz_0p99_detected": bool(
                public_slot[1] * Decimal("0.99") < ind_lip
            ),
            "missing_one_convolution_detected": bool(
                public_slot[0] / public_product_constant < ind_norm
            ),
            "duplicate_J2_detected_by_shape_guard": bool(
                public_slot[0] * Decimal(40) / ind_norm
                > MAX_EXPECTED_PUBLIC_RATIO
            ),
        }

        failed_guards: list[str] = []
        if min(ratios.values()) < 1.0:
            failed_guards.append("public bound under independent reconstruction")
        if ratios["slot_norm"] > float(MAX_EXPECTED_PUBLIC_RATIO):
            failed_guards.append("slot norm is too loose; possible duplicated operator factor")
        if ratios["slot_lipschitz"] > float(MAX_EXPECTED_PUBLIC_RATIO):
            failed_guards.append(
                "slot Lipschitz is too loose; possible duplicated operator factor"
            )
        if ratios["product_constant"] > float(MAX_EXPECTED_PUBLIC_RATIO):
            failed_guards.append("product constant drift exceeds preregistered shape guard")
        if ratios["derivative_after_j2"] > 1.000000000001:
            failed_guards.append("80/rho derivative factor drifted")
        if ratios["d_norm"] > 1.000000000001:
            failed_guards.append("fixed d coefficient norm drifted")
        if ratios["L_inverse"] > 1.000000000001:
            failed_guards.append("L-inverse coefficient norm drifted")

        stress = _fresh_fixed_multiplier_stress(calc)
        if not stress["finite"] or not stress["dominated"] or not stress["axis_nontrivial"]:
            failed_guards.append("fresh off-grid fixed-multiplier stress failed")
        if not all(negative_controls.values()):
            failed_guards.append("one or more preregistered negative controls were not detected")

        independent_values = {
            "rho": float(pieces["rho"]),
            "product_constant_pi355_over_113": float(pieces["product_constant"]),
            "post_J2_single_eta_factor_80_over_rho": float(
                pieces["derivative_after_j2"]
            ),
            "d_coefficient_norm": float(pieces["d_norm"]),
            "L_inverse_coefficient_norm": float(pieces["L_inverse"]),
            "fixed_operator_prefactor": float(pieces["fixed_prefactor"]),
            "slot_after_J2": {
                "norm": float(ind_norm),
                "lipschitz": float(ind_lip),
            },
            "upstream_public_independently_admitted_inputs": {
                "u": {
                    "norm": float(_pair(inputs["u"])[0]),
                    "lipschitz": float(_pair(inputs["u"])[1]),
                    "A4_PR": SOURCE_U_BALL_A4_PR,
                },
                "Phi": {
                    "norm": float(_pair(inputs["Phi"])[0]),
                    "lipschitz": float(_pair(inputs["Phi"])[1]),
                    "A4_PR": SOURCE_PHI_BALL_A4_PR,
                },
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
            "J2_reapplied": False,
            "convolution_count": 3,
        },
        "independent_values": independent_values,
        "public_over_independent_ratios": ratios,
        "preregistered_max_expected_public_ratio": float(MAX_EXPECTED_PUBLIC_RATIO),
        "fresh_offgrid_fixed_multiplier_stress": stress,
        "rho_local_operator_sensitivity": _rho_local_operator_sensitivity(calc),
        "negative_controls": negative_controls,
        "failed_guards": failed_guards,
        "source_R1_averaged_eta_ordinary_slot_independent_preflight_passed": not failed_guards,
        "dependency_state": {
            "source_axis_A4_PR": SOURCE_AXIS_A4_PR,
            "source_Phi_ball_A4_PR": SOURCE_PHI_BALL_A4_PR,
            "source_u_ball_A4_PR": SOURCE_U_BALL_A4_PR,
            "mixed_J_algebra_A4_PR": MIXED_J_A4_PR,
            "zeta_R1_slot_A4_PR": ZETA_R1_A4_PR,
            "zeta_R1_slot_A4_audit_pending_at_freeze": True,
            "prior_five_R1_algebraic_slots_have_independent_A4_receipt": False,
            "full_R1_independent_admission": False,
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
        _canonical_json(receipt_basis).encode("utf-8")
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
