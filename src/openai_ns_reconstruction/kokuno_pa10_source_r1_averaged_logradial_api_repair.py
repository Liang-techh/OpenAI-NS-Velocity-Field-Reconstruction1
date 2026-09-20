"""Representation-only repair for the Agent-4 PA.10 averaged-logradial R1 audit.

The original independent audit on exact head
``1e5b66dba633a868fc0b898178d712f3d640f022`` never produced a scientific
receipt.  Its focused regression stopped before receipt generation because the
validator tried to read ``calc.product_calculator.product_constant``.  The
public ``KokunoPA10PressureBallBounds`` product calculator intentionally exposes
``product(...)`` but no ``product_constant`` attribute.

This repair changes only how Agent 4 *observes the public product multiplier*.
It infers the effective multiplier from one ordinary public two-factor product
on the already frozen ``A(u)`` and ``Phi`` bounds, while the independent oracle
continues to use 80-digit Decimal arithmetic and the stricter rational
``pi < 355/113`` path from the failed audit.  All seeds, mutation factors,
shape guards, and final PDE thresholds remain unchanged.

No Agent-1 private Fraction/product/rounding helper or Agent-1 receipt is used as
a numerical oracle.  A PASS here can repair only the #736 averaged-logradial
slot admission; it is not a Navier-Stokes residual validation.
"""

from __future__ import annotations

import argparse
from decimal import Decimal, ROUND_CEILING, localcontext
import hashlib
import json
from pathlib import Path
import subprocess
from typing import Any

from .kokuno_pa10_source_r1_averaged_logradial_independent_audit import (
    AGENT1_HEAD,
    AGENT1_PR,
    IDENTITY_CASES,
    IDENTITY_SEED,
    MAX_PUBLIC_RATIO,
    OFFGRID_COUNT,
    SEED,
    SLOT_NAME,
    _canonical_json,
    _d,
    _fresh_eta_stress,
    _independent_values,
    _manufactured_averaging_identity_stress,
    _pair,
    _ratio,
    _weight_ratio_stress,
)
from .kokuno_pa10_source_r1_averaged_logradial_ordinary_slot import (
    KokunoPA10SourceR1AveragedLogradialOrdinarySlot,
)

SCHEMA = "kokuno-agent4-pa10-source-r1-averaged-logradial-api-repair-v1"
FAILED_A4_HEAD = "1e5b66dba633a868fc0b898178d712f3d640f022"
FAILED_A4_RUN = 35479780213


def _observable_public_product_constant(
    calc: KokunoPA10SourceR1AveragedLogradialOrdinarySlot,
) -> Decimal:
    """Infer the public two-factor norm multiplier without private attributes.

    Agent 1's stable public product API is exercised on the exact frozen
    ``A(u)`` and ``Phi`` bounds already used by the audited slot.  Dividing the
    returned public norm by the two public input norms observes the effective
    upward-rounded multiplier.  This number is used only on the *public* side
    of the public/independent ratio; the independent side remains the separate
    ``355/113`` Decimal reconstruction.
    """

    values = calc.operator_inputs()
    au = _pair(values["A_u"])[0]
    phi = _pair(values["Phi"])[0]
    if au <= 0 or phi <= 0:
        raise ValueError("public A(u)/Phi norms must be positive")
    public_product = calc.product_calculator.product(values["A_u"], values["Phi"])
    observed = _pair(public_product)[0] / (au * phi)
    if observed <= 0:
        raise ValueError("observed public product multiplier must be positive")
    return observed


def run_repaired_independent_audit() -> dict[str, Any]:
    """Re-run the frozen #736 protocol with only the public-API read repaired."""

    calc = KokunoPA10SourceR1AveragedLogradialOrdinarySlot()

    with localcontext() as ctx:
        ctx.prec = 80
        ctx.rounding = ROUND_CEILING

        independent, pieces = _independent_values(calc)
        public_core = _pair(calc.j2_average_times_logradial_phi())
        public_slot = _pair(calc.r1_averaged_logradial_after_j2()[SLOT_NAME])
        public_pc = _observable_public_product_constant(calc)
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
                public_plain * Decimal(39) / Decimal(40)
                < pieces["plain_after_J2"]
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

        try:
            exact_head = subprocess.check_output(
                ["git", "rev-parse", "HEAD"], text=True
            ).strip()
        except Exception:
            exact_head = None

        truth = {
            "failed_a4_head_remains_failed_without_receipt": True,
            "repair_is_representation_only": True,
            "source_R1_averaged_logradial_ordinary_slot_independently_audited": (
                not failed_guards
            ),
            "all_R1_ordinary_slots_independently_audited": False,
            "source_full_post_J2_R1_radius_one_ball_bound_machine_bound": False,
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
            "agent1_pr": AGENT1_PR,
            "agent1_exact_head": AGENT1_HEAD,
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
                    "product(A_u,Phi).norm/(A_u.norm*Phi.norm)"
                ),
            },
            "independence": {
                "decimal_precision": 80,
                "pi_upper": "355/113",
                "agent1_private_helpers_called": False,
                "agent1_receipt_used_as_numeric_oracle": False,
                "public_product_api_used_only_as_public_observation": True,
                "averaging_identity_checked_with_manufactured_polynomials": True,
                "fresh_offgrid_eta_sampling": True,
            },
            "frozen_protocol": {
                "offgrid_seed": SEED,
                "offgrid_count": OFFGRID_COUNT,
                "identity_seed": IDENTITY_SEED,
                "identity_cases": IDENTITY_CASES,
                "max_public_ratio": float(MAX_PUBLIC_RATIO),
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
            "observed_public_product_constant": float(public_pc),
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
    payload = run_repaired_independent_audit()
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
    print("failed_guards=", payload["failed_guards"])
    if payload["failed_guards"]:
        raise SystemExit(1)


if __name__ == "__main__":
    _main()
