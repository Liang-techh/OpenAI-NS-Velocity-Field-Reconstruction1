"""Representation-only repair for the Agent-4 PA.10 ``H_* Phi_eta`` R1 audit.

The original Agent-4 audit on exact head
``a45b9964f8bb36355c0a3bca999aeb2bee163c40`` did not produce a scientific
receipt.  Its focused regression stopped before receipt generation because the
validator tried to read ``calc.product_calculator.product_constant``.  The
stable public product calculator exposes ``product(...)`` but intentionally has
no ``product_constant`` attribute.

This repair changes only how Agent 4 observes the public two-factor product
multiplier.  It exercises the public ``product(...)`` API on the already frozen
``L^-1`` and ``H_*`` coefficient bounds and divides the returned norm by the two
input norms.  The independent oracle remains the failed audit's 80-digit
Decimal implementation with strict ``pi < 355/113``.  The original seed,
off-grid/axis-near probes, parameter perturbation, mutation factors, ratio gate,
and final PDE gates are unchanged.

A PASS here can repair only the #744 ``H_* Phi_eta`` ordinary post-J2 R1 slot.
It is not a Navier-Stokes residual validation and does not promote full R1,
M/K, global leading velocity, pressure, correction, or ``pde_validated``.
"""

from __future__ import annotations

import argparse
from decimal import Decimal, ROUND_CEILING, localcontext
import hashlib
import json
from pathlib import Path
import subprocess
from typing import Any

from .kokuno_pa10_source_r1_hstar_phi_eta_independent_audit import (
    AGENT1_HEAD,
    AGENT1_PR,
    DIVERGENCE_GATE,
    MAX_PUBLIC_RATIO,
    MOMENTUM_GATE,
    OFFGRID_COUNT,
    PHI_BALL_A4_PR,
    SEED,
    SLOT_NAME,
    SOURCE_AXIS_A4_PR,
    _canonical_json,
    _d,
    _fresh_hstar_stress,
    _independent_hstar_norm,
    _independent_l_inverse_norm,
    _independent_product_constant,
    _pair,
    _product,
    _ratio,
    _scale,
)
from .kokuno_pa10_source_r1_hstar_phi_eta_ordinary_slot import (
    KokunoPA10SourceR1HstarPhiEtaOrdinarySlot,
)

SCHEMA = "kokuno-agent4-pa10-source-r1-hstar-phi-eta-api-repair-v1"
FAILED_A4_HEAD = "a45b9964f8bb36355c0a3bca999aeb2bee163c40"
FAILED_A4_RUN = 35482164500


def _observable_public_product_constant(
    calc: KokunoPA10SourceR1HstarPhiEtaOrdinarySlot,
) -> Decimal:
    """Observe the public two-factor norm multiplier via the stable API only."""

    values = calc.operator_inputs()
    linv = _pair(values["L_inverse"])[0]
    hstar = _pair(values["H_star"])[0]
    if linv <= 0 or hstar <= 0:
        raise ValueError("public L^-1/H_* norms must be positive")
    public_product = calc.product_calculator.product(
        values["L_inverse"], values["H_star"]
    )
    observed = _pair(public_product)[0] / (linv * hstar)
    if observed <= 0:
        raise ValueError("observed public product multiplier must be positive")
    return observed


def run_repaired_independent_audit() -> dict[str, Any]:
    """Re-run the frozen #744 protocol with only its invalid public read repaired."""

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

        public_pc = _observable_public_product_constant(calc)
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

        try:
            exact_head = subprocess.check_output(
                ["git", "rev-parse", "HEAD"], text=True
            ).strip()
        except Exception:
            exact_head = None

        truth = {
            "failed_a4_head_remains_failed_without_receipt": True,
            "repair_is_representation_only": True,
            "source_R1_Hstar_Phi_eta_ordinary_slot_independently_audited": (
                not failed_guards
            ),
            "source_R1_all_ordinary_slots_independently_admitted": False,
            "source_full_post_J2_R1_independently_admitted": False,
            "source_operator_constant_M_machine_bound": False,
            "source_operator_constant_K_machine_bound": False,
            "global_pressure_matched": False,
            "global_leading_profile_reconstructed": False,
            "complete_kokuno_composite_velocity": False,
            "leading_only_ns_residual_assessed": False,
            "leading_plus_oscillatory_ns_residual_assessed": False,
            "after_correction_ns_residual_assessed": False,
            "heldout_ns_momentum_residual_assessed": False,
            "formal_full_domain_pde_gate_assessed": False,
            "pde_validated": False,
        }

        payload: dict[str, Any] = {
            "schema": SCHEMA,
            "exact_head": exact_head,
            "agent": 4,
            "audited_agent1_pr": AGENT1_PR,
            "audited_agent1_exact_head": AGENT1_HEAD,
            "superseded_failed_a4": {
                "exact_head": FAILED_A4_HEAD,
                "workflow_run": FAILED_A4_RUN,
                "receipt_generated": False,
                "failure_class": "validator_public_api_read_bug",
                "scientific_reject": False,
            },
            "repair": {
                "changed_scientific_formula": False,
                "changed_seed": False,
                "changed_mutation_factor": False,
                "changed_ratio_gate": False,
                "changed_final_pde_gate": False,
                "old_invalid_read": "product_calculator.product_constant",
                "new_public_observation": (
                    "product(L_inverse,H_star).norm/"
                    "(L_inverse.norm*H_star.norm)"
                ),
            },
            "independence": {
                "agent1_private_fraction_helpers_called": False,
                "agent1_private_product_helpers_called": False,
                "agent1_receipt_used_as_numerical_oracle": False,
                "public_product_api_used_only_as_public_observation": True,
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
            "frozen_protocol": {
                "offgrid_seed": SEED,
                "offgrid_count": OFFGRID_COUNT,
                "max_public_ratio": float(MAX_PUBLIC_RATIO),
                "Hstar_parameter_perturbation": "j0 +/-0.1%",
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
            "observed_public_product_constant": float(public_pc),
            "public_to_independent_ratios": ratios,
            "fresh_offgrid_Hstar_stress": offgrid,
            "negative_controls": negative_controls,
            "failed_guards": failed_guards,
            "source_R1_Hstar_Phi_eta_ordinary_slot_independent_preflight_passed": (
                not failed_guards
            ),
            "project_gates_unchanged": {
                "normalized_momentum_max": MOMENTUM_GATE,
                "normalized_momentum_L2": MOMENTUM_GATE,
                "divergence_max": DIVERGENCE_GATE,
                "divergence_L2": DIVERGENCE_GATE,
                "free_residual_defined_forcing_allowed": False,
            },
            "truth_boundary": truth,
            "limitations": [
                "This repair audits only the H_* Phi_eta ordinary post-J2 R1 slot.",
                "The coefficient-space values in this receipt are not Navier-Stokes residuals.",
                "The failed #744 head remains failed and has no scientific receipt.",
                "No global leading velocity, matched pressure, or real correction composite is evaluated.",
            ],
        }
        payload["receipt_sha256"] = hashlib.sha256(
            _canonical_json(payload).encode("utf-8")
        ).hexdigest()
        return payload


def save_report(path: str | Path) -> dict[str, Any]:
    payload = run_repaired_independent_audit()
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    payload = save_report(args.output)
    print("receipt_sha256=", payload["receipt_sha256"])
    print("ratios=", payload["public_to_independent_ratios"])
    print("negative_controls=", payload["negative_controls"])
    print("failed_guards=", payload["failed_guards"])
    if payload["failed_guards"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
