"""Independent high-precision audit of the source band-covering schedule.

The audit treats Agent 2's ``KokunoSourceBandCovering`` object as a value-only
public contract. Reference values are rebuilt with 80-digit Decimal arithmetic,
including adversarial cases whose covering-level logarithm lies just to either
side of an integer. Agent-2 bound/interaction helper methods are not used.

This is a local structural audit, not a Navier--Stokes residual calculation.
"""
from __future__ import annotations

import argparse
from decimal import Decimal, localcontext
import json
import math
from pathlib import Path
from typing import Any

from .kokuno_source_band_covering import KokunoSourceBandCovering

SCHEMA = "kokuno-agent4-source-band-covering-independent-audit-v1"
TASK_ID = "KOKUNO-A4-SOURCE-BAND-COVERING-INDEPENDENT-AUDIT-025"
BASE_PR = 430
BASE_HEAD = "167d8e750d7ba5e4abaeaec7fc93ccf70145b7fd"
DECIMAL_PRECISION = 80
H_CASES = (0.0003, 0.004, 0.0097)
ELL_CASES = (5, 6, 7, 8, 9, 10, 12, 16, 20, 32, 48, 58, 64, 83, 96, 128, 160, 192)
TRANSITION_BASES = ((12, 2), (18, 4), (45, 14), (83, 29))
TRANSITION_DELTAS = (-1.0e-9, -1.0e-12, 1.0e-12, 1.0e-9)
INTERACTION_H_CASES = (0.0003, 0.004, 0.0097)
INTERACTION_ELL_RANGE = (12, 96)

# Frozen before execution; local guards only.
VALUE_RELATIVE_GUARD = 5.0e-12
TRANSITION_LEVEL_MISMATCH_GUARD = 0
INTERACTION_BOUND_RELATIVE_GUARD = 5.0e-12
INTERACTION_VIOLATION_GUARD = 0
OFF_BY_ONE_MUTATION_RELATIVE_FLOOR = 1.0
MESH_EXPONENT_MUTATION_RELATIVE_FLOOR = 4.0
FORMAL_MOMENTUM_GATE = 1.0e-3
FORMAL_DIVERGENCE_GATE = 1.0e-5


def _d(value: float | int | str) -> Decimal:
    return Decimal(str(value))


def _decimal_constants() -> dict[str, Decimal]:
    with localcontext() as ctx:
        ctx.prec = DECIMAL_PRECISION
        root2 = Decimal(2).sqrt()
        t_g = Decimal(4) + root2
        lambda_g = Decimal(4) - root2
        rho_g = lambda_g.ln() / t_g.ln()
        return {
            "T_g": +t_g,
            "Lambda_g": +lambda_g,
            "rho_g": +rho_g,
            "kappa_s": Decimal("1e-5"),
        }


def _pow_decimal(base: Decimal, exponent: Decimal) -> Decimal:
    return (exponent * base.ln()).exp()


def _reference(ell: int, h: float) -> dict[str, Any]:
    constants = _decimal_constants()
    with localcontext() as ctx:
        ctx.prec = DECIMAL_PRECISION
        hh = _d(h)
        q = Decimal(2) ** (-int(ell))
        s_star = Decimal(int(ell) * int(ell))
        epsilon = _pow_decimal(q, hh)
        target = _pow_decimal(q, -(Decimal(1) + hh)) / s_star
        t_g = constants["T_g"]
        level = 0
        power = Decimal(1)
        while power * t_g <= target:
            power *= t_g
            level += 1
        c_i = (t_g ** level) * _pow_decimal(q, Decimal(1) + hh)
        d_r = Decimal(2) * (
            (Decimal(1) + hh) * constants["rho_g"] - hh * constants["kappa_s"]
        )
        m_i = (constants["Lambda_g"] ** level) * _pow_decimal(q, d_r / Decimal(2))
        mesh = s_star ** Decimal(-3)
        c_scaled = c_i * s_star
        m_scaled = m_i * _pow_decimal(epsilon, constants["kappa_s"]) * _pow_decimal(s_star, constants["rho_g"])
        return {
            "Q": +q,
            "epsilon": +epsilon,
            "S_star": +s_star,
            "mesh": +mesh,
            "covering_level": int(level),
            "c_i": +c_i,
            "M_i": +m_i,
            "rho_g": +constants["rho_g"],
            "d_r": +d_r,
            "c_scaled": +c_scaled,
            "M_scaled": +m_scaled,
            "target": +target,
        }


def _relative_error(observed: float, reference: Decimal) -> float:
    ref = float(reference)
    scale = max(abs(ref), 1.0e-300)
    return abs(float(observed) - ref) / scale


def _transition_h(ell: int, integer_level: int, delta: float) -> float:
    root2 = math.sqrt(2.0)
    t_g = 4.0 + root2
    return float(
        (((integer_level + delta) * math.log(t_g) + 2.0 * math.log(float(ell)))
         / (float(ell) * math.log(2.0)))
        - 1.0
    )


def _sweep_report() -> dict[str, Any]:
    rows: list[dict[str, Any]] = []
    worst_value_relative = 0.0
    level_mismatches = 0
    for h in H_CASES:
        for ell in ELL_CASES:
            ref = _reference(ell, h)
            public = KokunoSourceBandCovering(ell=ell, h=h)
            errors = {
                "Q": _relative_error(public.Q, ref["Q"]),
                "epsilon": _relative_error(public.epsilon, ref["epsilon"]),
                "product_grid_mesh": _relative_error(public.product_grid_mesh, ref["mesh"]),
                "c_i": _relative_error(public.c_i, ref["c_i"]),
                "M_i": _relative_error(public.M_i, ref["M_i"]),
            }
            worst_value_relative = max(worst_value_relative, *errors.values())
            mismatch = int(public.covering_level != ref["covering_level"])
            level_mismatches += mismatch
            rows.append({
                "ell": ell,
                "h": h,
                "public_level": public.covering_level,
                "reference_level": ref["covering_level"],
                "level_mismatch": bool(mismatch),
                "relative_errors": errors,
                "reference_c_scaled": float(ref["c_scaled"]),
                "reference_M_scaled": float(ref["M_scaled"]),
            })
    return {
        "case_count": len(rows),
        "worst_value_relative_error": worst_value_relative,
        "covering_level_mismatch_count": level_mismatches,
        "rows": rows,
    }


def _transition_report() -> dict[str, Any]:
    rows: list[dict[str, Any]] = []
    mismatches = 0
    minimum_distance = math.inf
    for ell, level in TRANSITION_BASES:
        for delta in TRANSITION_DELTAS:
            h = _transition_h(ell, level, delta)
            ref = _reference(ell, h)
            public = KokunoSourceBandCovering(ell=ell, h=h)
            constants = _decimal_constants()
            with localcontext() as ctx:
                ctx.prec = DECIMAL_PRECISION
                real_level = ref["target"].ln() / constants["T_g"].ln()
                nearest = real_level.to_integral_value()
                distance = abs(real_level - nearest)
            minimum_distance = min(minimum_distance, float(distance))
            mismatch = int(public.covering_level != ref["covering_level"])
            mismatches += mismatch
            rows.append({
                "ell": ell,
                "target_integer_level": level,
                "requested_delta": delta,
                "h": h,
                "reference_real_level": float(real_level),
                "reference_floor": ref["covering_level"],
                "public_floor": public.covering_level,
                "mismatch": bool(mismatch),
                "distance_to_nearest_integer": float(distance),
            })
    return {
        "case_count": len(rows),
        "mismatch_count": mismatches,
        "minimum_distance_to_integer": minimum_distance,
        "rows": rows,
    }


def _interaction_bound_decimal(h: float, ell0: int) -> Decimal:
    constants = _decimal_constants()
    with localcontext() as ctx:
        ctx.prec = DECIMAL_PRECISION
        hh = _d(h)
        numerator = Decimal(4) * (Decimal(1) + hh) * Decimal(2).ln() + Decimal(8) / Decimal(ell0)
        return +(Decimal(1) + numerator / constants["T_g"].ln())


def _interaction_report() -> dict[str, Any]:
    lo, hi = INTERACTION_ELL_RANGE
    checked = 0
    violations = 0
    worst_bound_relative = 0.0
    maximum_level_difference = 0
    tightest_slack = math.inf
    for h in INTERACTION_H_CASES:
        bands = {ell: KokunoSourceBandCovering(ell=ell, h=h) for ell in range(lo, hi + 1)}
        bound_ref = _interaction_bound_decimal(h, lo)
        bound_float = float(bound_ref)
        public_formula = 1.0 + (
            4.0 * (1.0 + h) * math.log(2.0) + 8.0 / lo
        ) / math.log(4.0 + math.sqrt(2.0))
        worst_bound_relative = max(
            worst_bound_relative,
            abs(public_formula - bound_float) / abs(bound_float),
        )
        for ell in range(lo, hi + 1):
            for delta in range(0, 5):
                other = ell + delta
                if other > hi:
                    continue
                actual = abs(bands[other].covering_level - bands[ell].covering_level)
                maximum_level_difference = max(maximum_level_difference, actual)
                slack = bound_float - actual
                tightest_slack = min(tightest_slack, slack)
                checked += 1
                if actual > bound_float:
                    violations += 1
    return {
        "pair_count": checked,
        "violation_count": violations,
        "maximum_level_difference": maximum_level_difference,
        "tightest_bound_slack": tightest_slack,
        "worst_bound_formula_relative_error": worst_bound_relative,
    }


def _mutation_report() -> dict[str, float]:
    off_by_one_changes: list[float] = []
    mesh_changes: list[float] = []
    constants = _decimal_constants()
    for h in H_CASES:
        for ell in (5, 12, 32, 83, 192):
            ref = _reference(ell, h)
            with localcontext() as ctx:
                ctx.prec = DECIMAL_PRECISION
                wrong_c = ref["c_i"] * constants["T_g"]
                wrong_m = ref["M_i"] * constants["Lambda_g"]
            off_by_one_changes.append(
                min(
                    _relative_error(float(wrong_c), ref["c_i"]),
                    _relative_error(float(wrong_m), ref["M_i"]),
                )
            )
            wrong_mesh = float(ell ** -5)
            mesh_changes.append(abs(wrong_mesh - float(ref["mesh"])) / float(ref["mesh"]))
    return {
        "weakest_off_by_one_level_relative_change": min(off_by_one_changes),
        "weakest_wrong_mesh_exponent_relative_change": min(mesh_changes),
    }


def generate_report() -> dict[str, Any]:
    sweep = _sweep_report()
    transition = _transition_report()
    interaction = _interaction_report()
    mutation = _mutation_report()
    checks = {
        "decimal_value_contract": sweep["worst_value_relative_error"] <= VALUE_RELATIVE_GUARD,
        "decimal_covering_levels": sweep["covering_level_mismatch_count"] <= TRANSITION_LEVEL_MISMATCH_GUARD,
        "near_transition_floor_cases": transition["mismatch_count"] <= TRANSITION_LEVEL_MISMATCH_GUARD,
        "interacting_band_bound_formula": interaction["worst_bound_formula_relative_error"] <= INTERACTION_BOUND_RELATIVE_GUARD,
        "interacting_band_bound_exhaustive": interaction["violation_count"] <= INTERACTION_VIOLATION_GUARD,
        "off_by_one_level_mutation_detected": mutation["weakest_off_by_one_level_relative_change"] >= OFF_BY_ONE_MUTATION_RELATIVE_FLOOR,
        "mesh_exponent_mutation_detected": mutation["weakest_wrong_mesh_exponent_relative_change"] >= MESH_EXPONENT_MUTATION_RELATIVE_FLOOR,
    }
    return {
        "schema": SCHEMA,
        "task_id": TASK_ID,
        "base_pr": BASE_PR,
        "base_head": BASE_HEAD,
        "independent_operator": "80-digit Decimal schedule reconstruction + transition-near floor tests; Agent-2 bound/interaction helpers unused",
        "sample_contract": {
            "h_cases": list(H_CASES),
            "ell_cases": list(ELL_CASES),
            "transition_bases": [list(x) for x in TRANSITION_BASES],
            "transition_deltas": list(TRANSITION_DELTAS),
            "interaction_ell_range": list(INTERACTION_ELL_RANGE),
            "uses_agent2_bound_diagnostics": False,
            "uses_agent2_interaction_receipt": False,
            "uses_training_tensor_or_loss": False,
            "uses_pressure_or_forcing_fit": False,
        },
        "frozen_guards": {
            "value_relative": VALUE_RELATIVE_GUARD,
            "transition_level_mismatches": TRANSITION_LEVEL_MISMATCH_GUARD,
            "interaction_bound_relative": INTERACTION_BOUND_RELATIVE_GUARD,
            "interaction_violations": INTERACTION_VIOLATION_GUARD,
            "off_by_one_mutation_relative_floor": OFF_BY_ONE_MUTATION_RELATIVE_FLOOR,
            "mesh_exponent_mutation_relative_floor": MESH_EXPONENT_MUTATION_RELATIVE_FLOOR,
        },
        "sweep": sweep,
        "transition": transition,
        "interaction": interaction,
        "mutation": mutation,
        "checks": checks,
        "local_structural_preflight_passed": bool(all(checks.values())),
        "formal_gates": {
            "normalized_momentum_max_and_L2": FORMAL_MOMENTUM_GATE,
            "divergence_max_and_L2": FORMAL_DIVERGENCE_GATE,
            "formal_full_domain_pde_gate_assessed": False,
            "pde_validated": False,
        },
        "truth_boundary": {
            "source_band_covering_schedule_independently_audited": bool(all(checks.values())),
            "actual_positive_order_background_bound": False,
            "actual_auxiliary_torus_mode_family_bound": False,
            "public_xyz_t_velocity_correction_materialized": False,
            "full_leading_oscillatory_correction_composite_available": False,
            "formal_full_domain_pde_gate_assessed": False,
            "pde_validated": False,
            "paper_exact": False,
            "openai_field_identified": False,
            "blowup_proved": False,
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = generate_report()
    text = json.dumps(report, indent=2, sort_keys=True)
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text + "\n", encoding="utf-8")
    print(text)


if __name__ == "__main__":
    main()
