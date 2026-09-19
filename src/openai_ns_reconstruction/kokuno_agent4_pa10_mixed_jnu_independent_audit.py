"""Independent Agent-4 audit of the PA.10 post-J_nu mixed seams.

This validator intentionally does not call Agent 1 private exact-Fraction helpers.
It treats the public mixed-seam API as the object under test and rebuilds a
strictly tighter admissible product constant with Decimal arithmetic and the
classical rational upper bound pi < 355/113.  Random positive ball envelopes
then exercise norm/Lipschitz propagation, the fixed d=1-eta^2 factor, and the
conditional R1/R2 mixed shapes.

Passing this audit clears only the conditional mixed-operator algebra.  It does
not manufacture the still-missing source u ball, full R1/R2, M/K, global
leading pressure/velocity, correction velocity, or any Navier-Stokes residual.
"""

from __future__ import annotations

import argparse
import copy
from decimal import Decimal, getcontext
import hashlib
import json
import math
from pathlib import Path
import random
from typing import Any

from .kokuno_pa10_mixed_jnu_seams import KokunoPA10MixedJnuSeamBounds
from .kokuno_pa10_remainder_ball_bounds import BallFactorBound

SCHEMA = "kokuno-agent4-pa10-mixed-jnu-independent-audit-v1"
SEED = 9173331
RANDOM_CASES = 64
MUTATION_FACTOR = Decimal("0.99")
MOMENTUM_GATE = Decimal("1e-3")
DIVERGENCE_GATE = Decimal("1e-5")


def _d(x: float | int | str | Decimal) -> Decimal:
    if isinstance(x, Decimal):
        return x
    return Decimal(str(x))


def _canonical_json(payload: dict[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _tight_product_constant() -> Decimal:
    # Independent from Agent 1's exact Fraction + pi<22/7 path.
    pi_upper = Decimal(355) / Decimal(113)
    c_sq = Decimal(4) * pi_upper * pi_upper / Decimal(3)
    return c_sq * c_sq


def _tight_mixed_factor(rho: Decimal, extras: int) -> Decimal:
    return (Decimal(80) / rho) * _tight_product_constant() ** (1 + extras)


def _independent_ball(
    *, rho: Decimal, factors: tuple[BallFactorBound, ...], extras: int
) -> tuple[Decimal, Decimal]:
    constant = _tight_mixed_factor(rho, extras)
    norms = tuple(_d(x.norm) for x in factors)
    lips = tuple(_d(x.lipschitz) for x in factors)
    norm = constant
    for value in norms:
        norm *= value
    lip = Decimal(0)
    for i, Li in enumerate(lips):
        if Li == 0:
            continue
        term = constant * Li
        for j, Nj in enumerate(norms):
            if i != j:
                term *= Nj
        lip += term
    return norm, lip


def _ratio(public: float, independent_upper: Decimal) -> float:
    if independent_upper == 0:
        return math.inf if public > 0 else 1.0
    return float(_d(public) / independent_upper)


def run_audit() -> dict[str, Any]:
    getcontext().prec = 80
    seam = KokunoPA10MixedJnuSeamBounds()
    rho = _d(seam.rho)
    rng = random.Random(SEED)

    tight_product = _tight_product_constant()
    public_product = _d(seam.product_constant_upper)
    product_ratio = public_product / tight_product

    mixed_ratios: list[float] = []
    mixed_mutations_detected: list[bool] = []
    for extras in range(4):
        independent = _tight_mixed_factor(rho, extras)
        public = _d(seam.mixed_operator_factor(extras))
        mixed_ratios.append(float(public / independent))
        mixed_mutations_detected.append(public * MUTATION_FACTOR < independent)

    domain = seam.phi_source.domain
    E = Decimal(1) + _d(domain.enlarged_real_margin)
    d_independent = max(
        Decimal(1),
        E * E - Decimal(1),
        Decimal(8) * E * rho,
        Decimal(9) * rho * rho,
    )
    d_public = seam.d_coefficient_ball()
    d_ratio = _ratio(d_public.norm, d_independent)

    min_norm_ratio = math.inf
    min_lip_ratio = math.inf
    random_records: list[dict[str, Any]] = []
    downward_mutations_detected = 0
    downward_mutations_total = 0

    for case in range(RANDOM_CASES):
        nu = 1 + (case % 2)
        extras_count = case % 4
        u = BallFactorBound(
            norm=10 ** rng.uniform(-1.0, 1.0),
            lipschitz=10 ** rng.uniform(-2.0, 0.7),
        )
        G = BallFactorBound(
            norm=10 ** rng.uniform(-1.0, 1.0),
            lipschitz=10 ** rng.uniform(-2.0, 0.7),
        )
        extras = tuple(
            BallFactorBound(
                norm=10 ** rng.uniform(-1.0, 1.0),
                lipschitz=(0.0 if j % 2 == 0 else 10 ** rng.uniform(-2.0, 0.3)),
            )
            for j in range(extras_count)
        )
        public = seam.averaged_eta_logradial_after_jnu(
            u=u, G=G, nu=nu, extra_factors=extras
        )
        independent_norm, independent_lip = _independent_ball(
            rho=rho, factors=(u, G, *extras), extras=extras_count
        )
        norm_ratio = _ratio(public.norm, independent_norm)
        lip_ratio = _ratio(public.lipschitz, independent_lip)
        min_norm_ratio = min(min_norm_ratio, norm_ratio)
        min_lip_ratio = min(min_lip_ratio, lip_ratio)

        norm_mutation = _d(public.norm) * MUTATION_FACTOR < independent_norm
        downward_mutations_total += 1
        downward_mutations_detected += int(norm_mutation)
        lip_mutation: bool | None = None
        if independent_lip > 0:
            lip_mutation = _d(public.lipschitz) * MUTATION_FACTOR < independent_lip
            downward_mutations_total += 1
            downward_mutations_detected += int(lip_mutation)

        if case < 8:
            random_records.append(
                {
                    "case": case,
                    "nu": nu,
                    "extra_factor_count": extras_count,
                    "norm_ratio_public_over_independent": norm_ratio,
                    "lipschitz_ratio_public_over_independent": lip_ratio,
                    "norm_mutation_detected": norm_mutation,
                    "lipschitz_mutation_detected": lip_mutation,
                }
            )

    phi = seam.source_phi_ball()
    r_shape_ratios: dict[str, dict[str, float]] = {}
    for label, u in {
        "unit": BallFactorBound(1.0, 1.0),
        "offgrid_a": BallFactorBound(0.731, 0.219),
        "offgrid_b": BallFactorBound(2.417, 0.083),
    }.items():
        r1_public = seam.r1_mixed_after_j2(u=u)
        r1_norm, r1_lip = _independent_ball(
            rho=rho, factors=(u, phi, d_public), extras=1
        )
        r2_public = seam.r2_mixed_after_j1(u=u)
        r2_norm, r2_lip = _independent_ball(
            rho=rho, factors=(u, u, d_public), extras=1
        )
        r_shape_ratios[label] = {
            "r1_norm_ratio": _ratio(r1_public.norm, r1_norm),
            "r1_lipschitz_ratio": _ratio(r1_public.lipschitz, r1_lip),
            "r2_norm_ratio": _ratio(r2_public.norm, r2_norm),
            "r2_lipschitz_ratio": _ratio(r2_public.lipschitz, r2_lip),
        }

    truth = seam.truth_boundary
    required_false = (
        "source_u_radius_one_ball_norm_machine_bound",
        "source_u_radius_one_ball_lipschitz_machine_bound",
        "source_R1_radius_one_ball_norm_machine_bound",
        "source_R2_radius_one_ball_norm_machine_bound",
        "source_operator_constant_M_machine_bound",
        "source_operator_constant_K_machine_bound",
        "global_pressure_matched",
        "global_leading_profile_reconstructed",
        "heldout_ns_residual_assessed",
        "pde_validated",
    )
    truth_boundary_fail_closed = all(truth[key] is False for key in required_false)

    all_r_ratios = [v for block in r_shape_ratios.values() for v in block.values()]
    guards = {
        "public_product_dominates_independent_tighter_upper": public_product >= tight_product,
        "all_public_mixed_factors_dominate_independent_tighter_upper": min(mixed_ratios) >= 1.0,
        "all_preregistered_mixed_factor_mutations_detected": all(mixed_mutations_detected),
        "public_d_ball_dominates_independent_coefficient_check": d_ratio >= 1.0,
        "random_public_norm_bounds_dominate_independent_bounds": min_norm_ratio >= 1.0,
        "random_public_lipschitz_bounds_dominate_independent_bounds": min_lip_ratio >= 1.0,
        "all_random_downward_mutations_detected": downward_mutations_detected == downward_mutations_total,
        "conditional_r1_r2_shapes_dominate_independent_bounds": min(all_r_ratios) >= 1.0,
        "truth_boundary_remains_fail_closed": truth_boundary_fail_closed,
    }
    failed = sorted(key for key, passed in guards.items() if not passed)

    payload: dict[str, Any] = {
        "schema": SCHEMA,
        "protocol": {
            "seed": SEED,
            "random_cases": RANDOM_CASES,
            "decimal_precision_digits": 80,
            "independent_pi_upper": "355/113",
            "agent1_pi_upper_not_reused": True,
            "mutation_factor": float(MUTATION_FACTOR),
            "final_momentum_gate_unchanged": float(MOMENTUM_GATE),
            "final_divergence_gate_unchanged": float(DIVERGENCE_GATE),
        },
        "operator_constant": {
            "rho": float(rho),
            "independent_tight_product_constant": float(tight_product),
            "public_product_constant": float(public_product),
            "product_ratio_public_over_independent": float(product_ratio),
            "mixed_factor_ratios_extras_0_to_3": mixed_ratios,
            "mixed_factor_downward_mutations_detected": mixed_mutations_detected,
        },
        "d_ball": {
            "independent_norm_upper": float(d_independent),
            "public_norm_upper": d_public.norm,
            "ratio_public_over_independent": d_ratio,
            "public_lipschitz": d_public.lipschitz,
        },
        "random_ball_stress": {
            "min_norm_ratio_public_over_independent": min_norm_ratio,
            "min_lipschitz_ratio_public_over_independent": min_lip_ratio,
            "downward_mutations_detected": downward_mutations_detected,
            "downward_mutations_total": downward_mutations_total,
            "sample_records": random_records,
        },
        "conditional_r1_r2_shapes": r_shape_ratios,
        "guards": guards,
        "failed_guards": failed,
        "mixed_jnu_independent_preflight_passed": not failed,
        "truth_boundary": {
            "conditional_mixed_operator_algebra_independently_audited": not failed,
            "source_u_ball_independently_audited": False,
            "full_R1_R2_independently_audited": False,
            "M_K_independently_audited": False,
            "global_pressure_independently_audited": False,
            "global_leading_velocity_independently_audited": False,
            "leading_only_ns_residual_assessed": False,
            "leading_plus_oscillatory_ns_residual_assessed": False,
            "after_correction_ns_residual_assessed": False,
            "heldout_ns_momentum_residual_assessed": False,
            "formal_full_domain_pde_gate_assessed": False,
            "pde_validated": False,
        },
    }
    identity = copy.deepcopy(payload)
    payload["receipt_sha256"] = hashlib.sha256(
        _canonical_json(identity).encode("utf-8")
    ).hexdigest()
    return payload


def save_report(path: str | Path) -> dict[str, Any]:
    report = run_audit()
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = run_audit()
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, indent=2, sort_keys=True))
    raise SystemExit(0 if report["mixed_jnu_independent_preflight_passed"] else 1)


if __name__ == "__main__":
    main()
