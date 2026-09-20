"""Independent Agent-4 audit of PA.10 averaged-logradial ordinary R1 slot.

Audits Agent-1 #733:
    J2[Lambda^-1 2D eta A(u) Y Phi_Y]
using the public identity
    A(u) Y Phi_Y = Y d_Y(A(u) Phi) - u Phi + A(u) Phi.

The implementation deliberately avoids Agent-1 private Fraction/product/rounding
helpers and receipt values as numerical oracles.  It reconstructs the
coefficient algebra with 80-digit Decimal arithmetic and pi < 355/113, adds
fresh off-grid / axis-near stress and a manufactured averaging-identity
mutation test, and preserves the fixed project PDE gates.
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

from .kokuno_pa10_source_r1_averaged_logradial_ordinary_slot import (
    KokunoPA10SourceR1AveragedLogradialOrdinarySlot,
)

SCHEMA = "kokuno-agent4-pa10-source-r1-averaged-logradial-independent-audit-v1"
AGENT1_PR = 733
AGENT1_HEAD = "e2333474234e0f509d6cfab04b94b01f5a952b04"
SEED = 9173411
IDENTITY_SEED = 9173412
OFFGRID_COUNT = 8192
IDENTITY_CASES = 512
SLOT_NAME = "lambda_inv_2D_eta_Au_times_Y_Phi_Y"
MAX_PUBLIC_RATIO = Decimal("1.01")


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


def _add(*pairs: tuple[Decimal, Decimal]) -> tuple[Decimal, Decimal]:
    return (
        sum((item[0] for item in pairs), Decimal(0)),
        sum((item[1] for item in pairs), Decimal(0)),
    )


def _product(
    *pairs: tuple[Decimal, Decimal], product_constant: Decimal
) -> tuple[Decimal, Decimal]:
    if not pairs:
        return Decimal(1), Decimal(0)
    algebra = product_constant ** max(len(pairs) - 1, 0)
    norm_product = math.prod(item[0] for item in pairs)
    lip = Decimal(0)
    for index, item in enumerate(pairs):
        if item[1] == 0:
            continue
        others = math.prod(
            other[0] for j, other in enumerate(pairs) if j != index
        )
        lip += item[1] * others
    return algebra * norm_product, algebra * lip


def _ratio(public: Decimal, independent: Decimal) -> float:
    if independent <= 0:
        raise ValueError("independent comparator must be positive")
    return float(public / independent)


def _independent_product_constant() -> Decimal:
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
        raise ValueError("invalid Cauchy geometry")
    x = rho_d / radius_d
    weight_sum = (Decimal(1) + x) / (Decimal(1) - x) ** 3
    return weight_sum / lower


def _independent_values(
    calc: KokunoPA10SourceR1AveragedLogradialOrdinarySlot,
) -> tuple[dict[str, tuple[Decimal, Decimal]], dict[str, Decimal]]:
    values = calc.operator_inputs()
    pc = _independent_product_constant()
    linv = _independent_l_inverse_norm(
        h=float(calc.domain.h),
        margin=float(calc.domain.enlarged_real_margin),
        radius=float(calc.domain.cauchy_radius),
        rho=float(calc.domain.coefficient_rho),
    )
    Au = _pair(values["A_u"])
    u = _pair(values["u"])
    phi = _pair(values["Phi"])
    eta = _pair(values["eta"])

    Au_phi = _product(Au, phi, product_constant=pc)
    u_phi = _product(u, phi, product_constant=pc)
    core = _add(
        _scale(Au_phi, Decimal(80)),
        _scale(u_phi, Decimal(40)),
        _scale(Au_phi, Decimal(40)),
    )
    linv_pair = (linv, Decimal(0))
    with_eta_linv = _product(linv_pair, eta, core, product_constant=pc)

    h = Decimal.from_float(float(calc.domain.h))
    D = Decimal("0.5") - h
    if D <= 0:
        raise ValueError("D must be positive")
    slot = _scale(with_eta_linv, Decimal(2) * D)

    return {"core": core, "slot": slot}, {
        "product_constant": pc,
        "L_inverse": linv,
        "single_logradial_after_J2": Decimal(80),
        "plain_after_J2": Decimal(40),
        "D": D,
        "effective_convolution_count_per_decomposition_term": Decimal(3),
    }


def _fresh_eta_stress(
    calc: KokunoPA10SourceR1AveragedLogradialOrdinarySlot,
) -> dict[str, Any]:
    domain = calc.domain
    E = 1.0 + float(domain.enlarged_real_margin)
    rng = np.random.default_rng(SEED)
    random_eta = rng.uniform(-E, E, size=OFFGRID_COUNT)
    probes = np.array([-E, -1.0e-12, 0.0, 1.0e-12, E], dtype=float)
    eta = np.concatenate((random_eta, probes))
    public_eta = float(calc.operator_inputs()["eta"].norm)
    max_abs = float(np.max(np.abs(eta)))

    h = float(domain.h)
    D = 0.5 - h
    coeff = 2.0 * D
    coeff_lo = 2.0 * (0.5 - h * 0.999)
    coeff_hi = 2.0 * (0.5 - h * 1.001)
    response = max(abs(coeff_lo - coeff), abs(coeff_hi - coeff))

    return {
        "seed": SEED,
        "fresh_random_offgrid_points": OFFGRID_COUNT,
        "axis_near_probes": [-1.0e-12, 0.0, 1.0e-12],
        "endpoint_probes": [-E, E],
        "max_abs_eta": max_abs,
        "public_eta_bound": public_eta,
        "finite": bool(np.isfinite(eta).all()),
        "dominated": bool(max_abs <= public_eta),
        "axis_nontrivial": bool(np.any(np.abs(probes[1:4]) < 1.0e-11)),
        "h_pm_0p1_percent_2D_response": float(response),
    }


def _manufactured_averaging_identity_stress() -> dict[str, Any]:
    """Check the source averaging identity without Agent-1 coefficient code."""
    rng = np.random.default_rng(IDENTITY_SEED)
    max_error = 0.0
    wrong_sign_error = 0.0
    max_signal = 0.0
    for _ in range(IDENTITY_CASES):
        degree_u = int(rng.integers(1, 7))
        degree_phi = int(rng.integers(1, 7))
        cu = rng.normal(size=degree_u + 1)
        cp = rng.normal(size=degree_phi + 1)
        Y = float(rng.uniform(0.05, 1.75))

        powers_u = np.array([Y**k for k in range(degree_u + 1)])
        powers_p = np.array([Y**k for k in range(degree_phi + 1)])
        u = float(np.dot(cu, powers_u))
        phi = float(np.dot(cp, powers_p))
        phi_y = float(
            sum(k * cp[k] * Y ** (k - 1) for k in range(1, degree_phi + 1))
        )
        Au = float(
            sum(cu[k] * Y**k / (k + 1) for k in range(degree_u + 1))
        )
        Au_y = float(
            sum(k * cu[k] * Y ** (k - 1) / (k + 1) for k in range(1, degree_u + 1))
        )
        lhs = Au * Y * phi_y
        rhs = Y * (Au_y * phi + Au * phi_y) - u * phi + Au * phi
        wrong = Y * (Au_y * phi + Au * phi_y) + u * phi + Au * phi

        max_error = max(max_error, abs(lhs - rhs))
        wrong_sign_error = max(wrong_sign_error, abs(lhs - wrong))
        max_signal = max(max_signal, abs(lhs), abs(rhs))

    return {
        "seed": IDENTITY_SEED,
        "cases": IDENTITY_CASES,
        "max_abs_identity_error": max_error,
        "wrong_sign_mutation_max_abs_error": wrong_sign_error,
        "identity_tolerance": 5.0e-12,
        "identity_passed": bool(max_error <= 5.0e-12),
        "wrong_sign_mutation_detected": bool(wrong_sign_error >= 1.0e-6),
        "nontrivial": bool(max_signal > 1.0e-3),
    }


def _weight_ratio_stress() -> dict[str, Any]:
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
    calc = KokunoPA10SourceR1AveragedLogradialOrdinarySlot()

    with localcontext() as ctx:
        ctx.prec = 80
        ctx.rounding = ROUND_CEILING

        independent, pieces = _independent_values(calc)
        public_core = _pair(calc.j2_average_times_logradial_phi())
        public_slot = _pair(calc.r1_averaged_logradial_after_j2()[SLOT_NAME])
        public_pc = _d(calc.product_calculator.product_constant)
        public_linv = _pair(calc.operator_inputs()["L_inverse"])[0]
        public_log = _d(calc.source_single_logradial_after_j2_factor_upper())
        public_plain = _d(calc.source_plain_after_j2_factor_upper())

        ratios = {
            "product_constant": _ratio(public_pc, pieces["product_constant"]),
            "L_inverse": _ratio(public_linv, pieces["L_inverse"]),
            "single_logradial_after_J2": _ratio(
                public_log, pieces["single_logradial_after_J2"]
            ),
            "plain_after_J2": _ratio(public_plain, pieces["plain_after_J2"]),
            "core_norm": _ratio(public_core[0], independent["core"][0]),
            "core_lipschitz": _ratio(public_core[1], independent["core"][1]),
            "slot_norm": _ratio(public_slot[0], independent["slot"][0]),
            "slot_lipschitz": _ratio(public_slot[1], independent["slot"][1]),
        }

        negative_controls = {
            "product_constant_0p998_detected": bool(
                public_pc * Decimal("0.998") < pieces["product_constant"]
            ),
            "L_inverse_0p999_detected": bool(
                public_linv * Decimal("0.999") < pieces["L_inverse"]
            ),
            "logradial_79_instead_of_80_detected": bool(
                public_log * Decimal(79) / Decimal(80)
                < pieces["single_logradial_after_J2"]
            ),
            "plain_39_instead_of_40_detected": bool(
                public_plain * Decimal(39) / Decimal(40) < pieces["plain_after_J2"]
            ),
            "core_norm_0p99_detected": bool(
                public_core[0] * Decimal("0.99") < independent["core"][0]
            ),
            "core_lipschitz_0p99_detected": bool(
                public_core[1] * Decimal("0.99") < independent["core"][1]
            ),
            "slot_norm_0p99_detected": bool(
                public_slot[0] * Decimal("0.99") < independent["slot"][0]
            ),
            "slot_lipschitz_0p99_detected": bool(
                public_slot[1] * Decimal("0.99") < independent["slot"][1]
            ),
            "duplicate_J2_factor_40_detected": bool(
                public_slot[0] * Decimal(40) / independent["slot"][0]
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

        offgrid = _fresh_eta_stress(calc)
        if not offgrid["finite"]:
            failed_guards.append("fresh_eta_stress_nonfinite")
        if not offgrid["dominated"]:
            failed_guards.append("fresh_eta_stress_not_dominated")
        if not offgrid["axis_nontrivial"]:
            failed_guards.append("axis_near_probe_missing")
        if offgrid["h_pm_0p1_percent_2D_response"] <= 0.0:
            failed_guards.append("h_parameter_perturbation_not_detected")

        identity = _manufactured_averaging_identity_stress()
        if not identity["identity_passed"]:
            failed_guards.append("averaging_identity_stress_failed")
        if not identity["wrong_sign_mutation_detected"]:
            failed_guards.append("averaging_identity_wrong_sign_mutation_not_detected")
        if not identity["nontrivial"]:
            failed_guards.append("averaging_identity_stress_trivial")

        weights = _weight_ratio_stress()
        if not weights["envelope_dominates"]:
            failed_guards.append("single_logradial_weight_envelope_failed")

        truth = {
            "source_axis_domain_independently_audited": True,
            "source_Phi_radius_one_ball_independently_audited": True,
            "source_u_radius_one_ball_independently_audited": True,
            "conditional_mixed_J_algebra_independently_audited": True,
            "source_R1_zeta_ordinary_slot_independently_audited": True,
            "source_R1_averaged_eta_ordinary_slot_independently_audited": True,
            "source_R1_Wstar_logradial_ordinary_slot_independent_audit_pending": True,
            "source_R1_averaged_logradial_ordinary_slot_independently_audited": not failed_guards,
            "source_R1_prior_five_algebraic_slots_independently_audited": False,
            "all_R1_ordinary_slots_independently_audited": False,
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

        payload: dict[str, Any] = {
            "schema": SCHEMA,
            "agent1_pr": AGENT1_PR,
            "agent1_exact_head": AGENT1_HEAD,
            "independence": {
                "decimal_precision": 80,
                "pi_upper": "355/113",
                "agent1_private_helpers_called": False,
                "agent1_receipt_used_as_numeric_oracle": False,
                "averaging_identity_checked_with_manufactured_polynomials": True,
                "fresh_offgrid_eta_sampling": True,
            },
            "independent_values": {
                "core": [float(v) for v in independent["core"]],
                "slot": [float(v) for v in independent["slot"]],
                "product_constant": float(pieces["product_constant"]),
                "L_inverse": float(pieces["L_inverse"]),
                "single_logradial_after_J2": 80.0,
                "plain_after_J2": 40.0,
                "effective_convolution_count_per_decomposition_term": 3,
            },
            "public_over_independent_ratios": ratios,
            "negative_controls": negative_controls,
            "fresh_offgrid_eta_stress": offgrid,
            "manufactured_averaging_identity_stress": identity,
            "weight_ratio_stress": weights,
            "failed_guards": failed_guards,
            "source_R1_averaged_logradial_ordinary_slot_independent_preflight_passed": (
                not failed_guards
            ),
            "truth_boundary": truth,
            "immutable_project_gates": {
                "normalized_momentum_max": 1.0e-3,
                "normalized_momentum_L2": 1.0e-3,
                "divergence_max": 1.0e-5,
                "divergence_L2": 1.0e-5,
                "free_residual_defined_forcing_allowed": False,
            },
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
    print("negative_controls=", payload["negative_controls"])
    print("identity_stress=", payload["manufactured_averaging_identity_stress"])
    print("offgrid_stress=", payload["fresh_offgrid_eta_stress"])
    if payload["failed_guards"]:
        raise SystemExit(1)


if __name__ == "__main__":
    _main()
