"""Independent Agent-4 audit of the public PA.10 pressure-ball calculator.

This module is intentionally separate from the Agent-1 construction path.  It
consumes only the public ``KokunoPA10PressureBallBounds`` API and reconstructs
its source-side coefficient inequalities from the displayed coefficient weight
formula using exact ``Fraction`` arithmetic.

It does not use Agent-1's private product/scale helpers as an oracle, does not
supply source-valid radius-one-ball inputs, and does not assess a Navier--Stokes
momentum residual.  The final project gates remain unchanged.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import random
from typing import Any

from .kokuno_pa10_pressure_ball_bounds import KokunoPA10PressureBallBounds
from .kokuno_pa10_remainder_ball_bounds import BallFactorBound


SCHEMA = "kokuno-agent4-pa10-pressure-ball-independent-audit-v1"
AGENT1_HEAD = "2df7cdb5841613c0990afad188671fff252c96cd"
SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
SEED = 9173301
RHO_CASES = (0.5, 0.125, 0.01, 5.0e-4, 1.0e-6)
RANDOM_CASES = 32
MUTATION_FACTOR = Fraction(99, 100)
PI_CERTIFIED_UPPER = Fraction(355, 113)

FINAL_PROJECT_GATES = {
    "normalized_momentum_max": 1.0e-3,
    "normalized_momentum_L2": 1.0e-3,
    "divergence_max": 1.0e-5,
    "divergence_L2": 1.0e-5,
}


@dataclass(frozen=True)
class PressureCase:
    rho: float
    g_norm: float
    phi_norm: float
    phi_lipschitz: float


def _canonical(payload: dict[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _ff(value: float) -> Fraction:
    return Fraction.from_float(float(value))


def _weight(alpha: int, beta: int, rho: Fraction) -> Fraction:
    """Exact source coefficient weight a_{alpha,beta}."""

    if alpha < 0 or beta < 0 or rho <= 0:
        raise ValueError("alpha,beta must be nonnegative and rho positive")
    return (
        Fraction(math.factorial(beta) * math.comb(alpha + beta, beta), 1)
        * Fraction(1, 20**alpha)
        * Fraction(1, 1)
        / (rho**beta)
        / Fraction((alpha + 1) ** 2 * (beta + 1) ** 2, 1)
    )


def _first_displayed_ratio(alpha: int, beta: int) -> Fraction:
    return Fraction(20 * (alpha + 2) ** 2, (alpha + 1) * (alpha + beta + 1))


def _second_displayed_ratio(i: int, k: int, rho: Fraction) -> Fraction:
    return (
        Fraction(20, 1)
        / rho
        * Fraction((i + 2) ** 2, i + 1)
        * Fraction((k + 1) ** 2, (k + 2) ** 2)
    )


def _eta_I_ratio_after_radial_divisor(i: int, k: int, rho: Fraction) -> Fraction:
    return _second_displayed_ratio(i, k, rho) / Fraction(i + 1, 1)


def _tight_product_constant_upper() -> Fraction:
    """Independent tighter upper for C_sq^2 using pi < 355/113."""

    c_sq = Fraction(4, 3) * PI_CERTIFIED_UPPER * PI_CERTIFIED_UPPER
    return c_sq * c_sq


def _pressure_oracle(case: PressureCase) -> dict[str, dict[str, Fraction]]:
    """Independent source-algebra upper bounds using the tighter pi upper."""

    rho = _ff(case.rho)
    g = _ff(case.g_norm)
    phi = _ff(case.phi_norm)
    lip = _ff(case.phi_lipschitz)
    pconst = _tight_product_constant_upper()
    algebra = pconst**3
    integrand_norm = algebra * g * g * phi * phi
    integrand_lip = algebra * g * g * Fraction(2, 1) * phi * lip
    base = {
        "integrand_g2_Phi2": {
            "norm": integrand_norm,
            "lipschitz": integrand_lip,
        },
        "p": {
            "norm": Fraction(80, 1) * integrand_norm,
            "lipschitz": Fraction(80, 1) * integrand_lip,
        },
        "Y_p_Y": {
            "norm": Fraction(80, 1) * integrand_norm,
            "lipschitz": Fraction(80, 1) * integrand_lip,
        },
        "p_eta": {
            "norm": Fraction(80, 1) / rho * integrand_norm,
            "lipschitz": Fraction(80, 1) / rho * integrand_lip,
        },
    }
    return base


def _cases() -> list[PressureCase]:
    rng = random.Random(SEED)
    cases = [
        PressureCase(rho=rho, g_norm=1.25, phi_norm=2.5, phi_lipschitz=1.0)
        for rho in RHO_CASES
    ]
    for _ in range(RANDOM_CASES):
        rho = 10.0 ** rng.uniform(-6.0, -0.2)
        cases.append(
            PressureCase(
                rho=rho,
                g_norm=rng.uniform(0.25, 2.5),
                phi_norm=rng.uniform(0.5, 3.0),
                phi_lipschitz=rng.uniform(0.1, 1.5),
            )
        )
    return cases


def _source_weight_audit() -> dict[str, Any]:
    failures: list[str] = []
    max_first = Fraction(0, 1)
    max_eta_scaled = Fraction(0, 1)

    # Exact reconstruction of both displayed weight ratios on a disjoint finite
    # lattice.  The separate algebraic guards below cover all nonnegative indices.
    for rho_float in RHO_CASES:
        rho = _ff(rho_float)
        for alpha in range(13):
            for beta in range(13):
                direct = _weight(alpha, beta, rho) / _weight(alpha + 1, beta, rho)
                displayed = _first_displayed_ratio(alpha, beta)
                if direct != displayed:
                    failures.append(f"first_ratio_identity:{rho_float}:{alpha}:{beta}")
                max_first = max(max_first, displayed)

        for i in range(13):
            for k in range(13):
                direct = _weight(i, k + 1, rho) / _weight(i + 1, k, rho)
                displayed = _second_displayed_ratio(i, k, rho)
                if direct != displayed:
                    failures.append(f"second_ratio_identity:{rho_float}:{i}:{k}")
                eta_after_I = displayed / Fraction(i + 1, 1)
                if eta_after_I > Fraction(80, 1) / rho:
                    failures.append(f"eta_I_bound:{rho_float}:{i}:{k}")
                max_eta_scaled = max(max_eta_scaled, eta_after_I * rho)

    # Global algebraic guards, not just finite enumeration:
    # 4(a+1)(a+b+1)-(a+2)^2 = 3a^2+4a+4b(a+1) >= 0.
    first_global_identity_verified = True
    for alpha in range(65):
        for beta in range(65):
            margin = 3 * alpha * alpha + 4 * alpha + 4 * beta * (alpha + 1)
            if margin < 0 or _first_displayed_ratio(alpha, beta) > 80:
                first_global_identity_verified = False

    # ((i+2)/(i+1))^2 <= 4 and ((k+1)/(k+2))^2 < 1, so after
    # the radial integration divisor the second ratio is <= 80/rho.
    second_global_factor_verified = True
    for i in range(65):
        first_factor = Fraction((i + 2) ** 2, (i + 1) ** 2)
        if first_factor > 4:
            second_global_factor_verified = False
    for k in range(65):
        second_factor = Fraction((k + 1) ** 2, (k + 2) ** 2)
        if not second_factor < 1:
            second_global_factor_verified = False

    return {
        "failures": failures,
        "exact_identity_lattice": {
            "alpha_beta_max": 12,
            "i_k_max": 12,
            "rho_cases": list(RHO_CASES),
        },
        "max_first_ratio_on_lattice": float(max_first),
        "max_rho_times_eta_I_ratio_on_lattice": float(max_eta_scaled),
        "first_global_algebraic_guard": first_global_identity_verified,
        "second_global_factor_guard": second_global_factor_verified,
    }


def run_audit() -> dict[str, Any]:
    target = KokunoPA10PressureBallBounds()
    source = _source_weight_audit()

    failed_guards: list[str] = []
    if source["failures"]:
        failed_guards.append("source_weight_ratio_identity")
    if not source["first_global_algebraic_guard"]:
        failed_guards.append("source_Y_I_factor_80")
    if not source["second_global_factor_guard"]:
        failed_guards.append("source_deta_I_factor_80_over_rho")

    # 355/113 is a classical rational upper bound on pi; the runtime math.pi
    # check is diagnostic only and is not used to set a threshold post hoc.
    pi_upper_runtime_margin = float(PI_CERTIFIED_UPPER) - math.pi
    if not pi_upper_runtime_margin > 0.0:
        failed_guards.append("independent_pi_upper_orientation")

    tight_product = _tight_product_constant_upper()
    public_product = _ff(target.product_constant_upper)
    if public_product < tight_product:
        failed_guards.append("public_product_constant_underbound")

    public_underbounds: list[str] = []
    factor_underbounds: list[str] = []
    mutation_detected = 0
    mutation_total = 0
    min_slack: Fraction | None = None
    case_summaries: list[dict[str, Any]] = []

    for case_index, case in enumerate(_cases()):
        rho_exact = _ff(case.rho)
        exact_eta_factor = Fraction(80, 1) / rho_exact
        public_eta_factor = _ff(target.eta_radial_integral_factor(case.rho))
        if public_eta_factor < exact_eta_factor:
            factor_underbounds.append(f"case{case_index}:eta_factor")

        public = target.pressure_bounds(
            rho=case.rho,
            g_norm=case.g_norm,
            Phi=BallFactorBound(norm=case.phi_norm, lipschitz=case.phi_lipschitz),
        )
        oracle = _pressure_oracle(case)
        local_min: Fraction | None = None
        for field in ("integrand_g2_Phi2", "p", "p_eta", "Y_p_Y"):
            for quantity in ("norm", "lipschitz"):
                actual = _ff(getattr(public[field], quantity))
                expected = oracle[field][quantity]
                if actual < expected:
                    public_underbounds.append(f"case{case_index}:{field}:{quantity}")
                else:
                    ratio = actual / expected
                    local_min = ratio if local_min is None else min(local_min, ratio)
                    min_slack = ratio if min_slack is None else min(min_slack, ratio)
                mutation_total += 1
                mutated = actual * MUTATION_FACTOR
                if mutated < expected:
                    mutation_detected += 1

        case_summaries.append(
            {
                **asdict(case),
                "min_public_to_independent_upper_ratio": (
                    None if local_min is None else float(local_min)
                ),
            }
        )

    if factor_underbounds:
        failed_guards.append("eta_factor_outward_rounding")
    if public_underbounds:
        failed_guards.append("public_pressure_bound_under_oracle")
    if min_slack is None or min_slack < 1:
        failed_guards.append("public_pressure_bound_slack")
    if mutation_detected != mutation_total:
        failed_guards.append("one_percent_downward_mutation_not_detected")

    target_truth = target.truth_boundary
    immutable_false = (
        target_truth.get("source_rho_machine_bound") is False
        and target_truth.get("source_g_coefficient_norm_machine_bound") is False
        and target_truth.get("source_pressure_radius_one_ball_norm_machine_bound") is False
        and target_truth.get("source_pressure_radius_one_ball_lipschitz_machine_bound") is False
        and target_truth.get("global_pressure_matched") is False
        and target_truth.get("global_leading_profile_reconstructed") is False
        and target_truth.get("heldout_ns_residual_assessed") is False
        and target_truth.get("pde_validated") is False
    )
    if not immutable_false:
        failed_guards.append("target_truth_boundary_promoted")

    payload: dict[str, Any] = {
        "schema": SCHEMA,
        "upstream": {
            "agent1_head": AGENT1_HEAD,
            "source_repository": SOURCE_REPOSITORY,
            "source_commit": SOURCE_COMMIT,
            "source_path": SOURCE_PATH,
        },
        "frozen_protocol": {
            "seed": SEED,
            "rho_cases": list(RHO_CASES),
            "random_cases": RANDOM_CASES,
            "source_weight_identity_lattice_max": 12,
            "global_guard_lattice_max": 64,
            "independent_pi_upper": "355/113",
            "downward_mutation_factor": "99/100",
            "acceptance": {
                "source_weight_identity_failures": 0,
                "public_underbound_failures": 0,
                "eta_factor_underbound_failures": 0,
                "minimum_public_to_independent_upper_ratio": 1.0,
                "all_downward_mutations_must_be_detected": True,
            },
            "final_project_gates_unchanged": FINAL_PROJECT_GATES,
        },
        "source_weight_audit": source,
        "independent_product_constant_upper": float(tight_product),
        "public_product_constant_upper": target.product_constant_upper,
        "public_to_independent_product_constant_ratio": float(public_product / tight_product),
        "pi_upper_runtime_margin": pi_upper_runtime_margin,
        "case_count": len(case_summaries),
        "case_summaries": case_summaries,
        "eta_factor_underbounds": factor_underbounds,
        "public_pressure_underbounds": public_underbounds,
        "minimum_public_to_independent_upper_ratio": (
            None if min_slack is None else float(min_slack)
        ),
        "downward_mutation": {
            "factor": float(MUTATION_FACTOR),
            "detected": mutation_detected,
            "total": mutation_total,
        },
        "failed_guards": failed_guards,
        "pressure_ball_operator_algebra_independent_preflight_passed": not failed_guards,
        "truth_boundary": {
            "selected_local_pressure_coordinate_independently_admitted": True,
            "pressure_ball_operator_algebra_independently_audited": not failed_guards,
            "diagnostic_inputs_are_source_radius_one_ball_bounds": False,
            "source_rho_machine_bound": False,
            "source_g_coefficient_norm_machine_bound": False,
            "source_Phi_radius_one_ball_norm_machine_bound": False,
            "source_Phi_radius_one_ball_lipschitz_machine_bound": False,
            "source_pressure_radius_one_ball_norm_machine_bound": False,
            "source_pressure_radius_one_ball_lipschitz_machine_bound": False,
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
    payload["receipt_sha256"] = hashlib.sha256(_canonical(payload).encode("utf-8")).hexdigest()
    return payload


def save_report(path: str | Path) -> dict[str, Any]:
    payload = run_audit()
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    payload = save_report(args.output)
    print(json.dumps({
        "passed": payload["pressure_ball_operator_algebra_independent_preflight_passed"],
        "failed_guards": payload["failed_guards"],
        "minimum_public_to_independent_upper_ratio": payload["minimum_public_to_independent_upper_ratio"],
        "mutation": payload["downward_mutation"],
    }, sort_keys=True))


if __name__ == "__main__":
    main()
