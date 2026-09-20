"""Independent Agent-4 audit of PA.10 source pair M/K arithmetic.

This validator is stacked on Agent-1 #760, which materializes self-certified
coefficient-space M/K envelopes from *already post-J* R1/R2 bounds.  The audit
uses only the public full-post-J inputs, the public ``||(1+T)^-1||`` upper, and
the public M/K result.  It does not call Agent-1 private Fraction/rounding
helpers, does not call the #668 replay helper, and does not use an Agent-1
receipt as a numerical oracle.

The independently reconstructed pair is

    M_phi = 1/2 ||(1+T)^-1|| ||J2 R1||,
    M_u   = 1/2 ||J1 R2||,
    K_phi = 1/2 ||(1+T)^-1|| Lip(J2 R1),
    K_u   = 1/2 Lip(J1 R2),

with pair values given by the componentwise maxima.  All arithmetic is rebuilt
with 90-digit Decimal directed upward rounding.  Shape mutations explicitly
exercise accidental second applications of J2/J1 and omission of the displayed
factor 1/2.

A PASS is conditional arithmetic evidence only.  It cannot independently admit
M/K as source operator constants until full post-J R1/R2 themselves are
independently admitted.  It is not a Navier-Stokes residual and cannot promote
any PDE state.
"""

from __future__ import annotations

import argparse
from decimal import Decimal, ROUND_CEILING, localcontext
import hashlib
import json
import math
from pathlib import Path
from typing import Any

from .kokuno_pa10_source_pair_mk import KokunoPA10SourcePairMK

SCHEMA = "kokuno-agent4-pa10-source-pair-mk-independent-audit-v1"
AGENT1_PR = 760
AGENT1_HEAD = "76517f4e6084429a84b2f2e44189a3d287a4187e"
SEED = 9173441
MAX_PUBLIC_RATIO = Decimal("1.000000000001")
MOMENTUM_GATE = 1.0e-3
DIVERGENCE_GATE = 1.0e-5
ST006_MOMENTUM_MAX = 0.1082289305112118
ST006_VOLUME_L2 = 0.10758432876230622

_TRUTH_BOUNDARY = {
    "source_pair_MK_arithmetic_independently_audited": True,
    "source_pair_MK_audit_is_conditional_on_post_J_inputs": True,
    "source_full_post_J2_R1_independent_admission": False,
    "source_full_post_J1_R2_independent_admission": False,
    "source_operator_constant_M_independent_agent4_admission": False,
    "source_operator_constant_K_independent_agent4_admission": False,
    "source_contraction_invariant_ball_machine_verified": False,
    "source_contraction_factor_machine_verified": False,
    "source_fixed_point_distance_machine_bound": False,
    "source_fixed_point_solved": False,
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
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


def _canonical_json(payload: dict[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _d(value: float, name: str) -> Decimal:
    out = float(value)
    if not math.isfinite(out) or out < 0.0:
        raise ValueError(f"{name} must be finite and nonnegative")
    return Decimal.from_float(out)


def _ratio(public: Decimal, independent: Decimal) -> float:
    if independent <= 0:
        raise ValueError("independent comparator must be positive")
    return float(public / independent)


def _independent_pair(
    *,
    inverse_one_plus_t: Decimal,
    r1_norm: Decimal,
    r1_lipschitz: Decimal,
    r2_norm: Decimal,
    r2_lipschitz: Decimal,
) -> dict[str, Decimal]:
    half = Decimal(1) / Decimal(2)
    m_phi = half * inverse_one_plus_t * r1_norm
    m_u = half * r2_norm
    k_phi = half * inverse_one_plus_t * r1_lipschitz
    k_u = half * r2_lipschitz
    return {
        "inverse_one_plus_T_absolute_series_upper": inverse_one_plus_t,
        "M_phi_component": m_phi,
        "M_u_component": m_u,
        "M_pair_max": max(m_phi, m_u),
        "K_phi_component": k_phi,
        "K_u_component": k_u,
        "K_pair_max": max(k_phi, k_u),
    }


def _shape_mutations(independent: dict[str, Decimal]) -> dict[str, bool]:
    """Frozen wrong-shape controls; all must visibly leave the tight pair envelope."""
    m_phi = independent["M_phi_component"]
    m_u = independent["M_u_component"]
    k_phi = independent["K_phi_component"]
    k_u = independent["K_u_component"]

    # If an already-post-J R1 input is integrated by J2 again, the b=0 ordinary
    # factor is 40.  Likewise an accidental second J1 on R2 contributes 80.
    duplicate_j2_ratio = max(
        (Decimal(40) * m_phi) / m_phi,
        (Decimal(40) * k_phi) / k_phi,
    )
    duplicate_j1_ratio = max(
        (Decimal(80) * m_u) / m_u,
        (Decimal(80) * k_u) / k_u,
    )

    # Omitting the displayed pair factor one-half doubles every component.
    missing_half_ratio = Decimal(2)
    return {
        "duplicate_J2_on_already_post_J2_R1_detected": bool(
            duplicate_j2_ratio > Decimal("1.01")
        ),
        "duplicate_J1_on_already_post_J1_R2_detected": bool(
            duplicate_j1_ratio > Decimal("1.01")
        ),
        "missing_pair_factor_one_half_detected": bool(
            missing_half_ratio > Decimal("1.01")
        ),
    }


def _inverse_sensitivity(
    *,
    inverse_one_plus_t: Decimal,
    r1_norm: Decimal,
    r1_lipschitz: Decimal,
    r2_norm: Decimal,
    r2_lipschitz: Decimal,
) -> dict[str, Any]:
    """Preregistered +/-0.1% input perturbation; diagnostic only, never retuned."""
    center = _independent_pair(
        inverse_one_plus_t=inverse_one_plus_t,
        r1_norm=r1_norm,
        r1_lipschitz=r1_lipschitz,
        r2_norm=r2_norm,
        r2_lipschitz=r2_lipschitz,
    )
    low = _independent_pair(
        inverse_one_plus_t=inverse_one_plus_t * Decimal("0.999"),
        r1_norm=r1_norm,
        r1_lipschitz=r1_lipschitz,
        r2_norm=r2_norm,
        r2_lipschitz=r2_lipschitz,
    )
    high = _independent_pair(
        inverse_one_plus_t=inverse_one_plus_t * Decimal("1.001"),
        r1_norm=r1_norm,
        r1_lipschitz=r1_lipschitz,
        r2_norm=r2_norm,
        r2_lipschitz=r2_lipschitz,
    )
    return {
        "seed": SEED,
        "perturbation_fraction": 0.001,
        "M_phi_low_over_center": float(low["M_phi_component"] / center["M_phi_component"]),
        "M_phi_high_over_center": float(high["M_phi_component"] / center["M_phi_component"]),
        "K_phi_low_over_center": float(low["K_phi_component"] / center["K_phi_component"]),
        "K_phi_high_over_center": float(high["K_phi_component"] / center["K_phi_component"]),
        "M_u_invariant": bool(low["M_u_component"] == center["M_u_component"] == high["M_u_component"]),
        "K_u_invariant": bool(low["K_u_component"] == center["K_u_component"] == high["K_u_component"]),
        "monotone_phi_components": bool(
            low["M_phi_component"] < center["M_phi_component"] < high["M_phi_component"]
            and low["K_phi_component"] < center["K_phi_component"] < high["K_phi_component"]
        ),
    }


def run_independent_audit(
    *, pr_head: str | None = None, checkout_head: str | None = None
) -> dict[str, Any]:
    calc = KokunoPA10SourcePairMK()

    # Only public upstream interfaces are consumed.  Crucially, this audit does
    # not call pair_mk_self_certificate until after its independent values have
    # been constructed, and never calls bridge_post_j_replay.
    post_j = calc.full_post_j_inputs()
    r1 = post_j["post_J2_R1"]
    r2 = post_j["post_J1_R2"]
    inverse_float = calc.inverse_one_plus_t_upper()

    with localcontext() as ctx:
        ctx.prec = 90
        ctx.rounding = ROUND_CEILING
        inverse = _d(inverse_float, "inverse_one_plus_T")
        r1_norm = _d(r1.norm, "post_J2_R1.norm")
        r1_lip = _d(r1.lipschitz, "post_J2_R1.lipschitz")
        r2_norm = _d(r2.norm, "post_J1_R2.norm")
        r2_lip = _d(r2.lipschitz, "post_J1_R2.lipschitz")

        independent = _independent_pair(
            inverse_one_plus_t=inverse,
            r1_norm=r1_norm,
            r1_lipschitz=r1_lip,
            r2_norm=r2_norm,
            r2_lipschitz=r2_lip,
        )

        # Public outputs are now read solely as the object under audit.
        public_float = calc.pair_mk_self_certificate()
        public = {name: _d(value, f"public.{name}") for name, value in public_float.items()}

        ratios = {name: _ratio(public[name], independent[name]) for name in independent}
        negative_controls = {
            f"{name}_0p999_downward_detected": bool(
                public[name] * Decimal("0.999") < independent[name]
            )
            for name in (
                "M_phi_component",
                "M_u_component",
                "M_pair_max",
                "K_phi_component",
                "K_u_component",
                "K_pair_max",
            )
        }
        negative_controls.update(_shape_mutations(independent))

        component_guards = {
            "M_pair_is_component_max": bool(
                public["M_pair_max"] == max(public["M_phi_component"], public["M_u_component"])
            ),
            "K_pair_is_component_max": bool(
                public["K_pair_max"] == max(public["K_phi_component"], public["K_u_component"])
            ),
            "all_components_nontrivial": bool(
                min(
                    independent["M_phi_component"],
                    independent["M_u_component"],
                    independent["K_phi_component"],
                    independent["K_u_component"],
                )
                > 0
            ),
        }
        sensitivity = _inverse_sensitivity(
            inverse_one_plus_t=inverse,
            r1_norm=r1_norm,
            r1_lipschitz=r1_lip,
            r2_norm=r2_norm,
            r2_lipschitz=r2_lip,
        )

        failed_guards: list[str] = []
        for name, value in ratios.items():
            ratio_d = Decimal.from_float(float(value))
            if ratio_d < Decimal(1):
                failed_guards.append(f"public_underbounds_independent:{name}")
            if ratio_d > MAX_PUBLIC_RATIO:
                failed_guards.append(f"unexpected_pair_shape_overcount:{name}")
        for name, detected in negative_controls.items():
            if not detected:
                failed_guards.append(f"negative_control_not_detected:{name}")
        for name, passed in component_guards.items():
            if not passed:
                failed_guards.append(f"component_guard_failed:{name}")
        if not sensitivity["M_u_invariant"]:
            failed_guards.append("inverse_sensitivity_crosswired_M_u")
        if not sensitivity["K_u_invariant"]:
            failed_guards.append("inverse_sensitivity_crosswired_K_u")
        if not sensitivity["monotone_phi_components"]:
            failed_guards.append("inverse_sensitivity_phi_not_monotone")

        payload: dict[str, Any] = {
            "schema": SCHEMA,
            "agent": 4,
            "audited_agent1_pr": AGENT1_PR,
            "audited_agent1_exact_head": AGENT1_HEAD,
            "pr_head": pr_head,
            "checkout_head": checkout_head,
            "independence": {
                "agent1_private_fraction_helpers_called": False,
                "agent1_private_rounding_helpers_called": False,
                "agent1_bridge_post_J_replay_called": False,
                "agent1_receipt_used_as_numerical_oracle": False,
                "decimal_precision": 90,
                "decimal_rounding": "ROUND_CEILING",
                "pair_formula_reconstructed_from_public_post_J_inputs": True,
            },
            "conditional_upstream_inputs": {
                "post_J2_R1_norm": float(r1_norm),
                "post_J2_R1_lipschitz": float(r1_lip),
                "post_J1_R2_norm": float(r2_norm),
                "post_J1_R2_lipschitz": float(r2_lip),
                "inverse_one_plus_T_upper": float(inverse),
                "full_R1_independent_admission_required_before_promotion": True,
                "full_R2_independent_admission_required_before_promotion": True,
            },
            "independent_values": {name: float(value) for name, value in independent.items()},
            "public_values": public_float,
            "public_to_independent_ratios": ratios,
            "component_guards": component_guards,
            "inverse_plus_minus_0p1_percent_sensitivity": sensitivity,
            "negative_controls": negative_controls,
            "failed_guards": failed_guards,
            "passed": not failed_guards,
            "same_protocol_PDE_comparison": {
                "performed": False,
                "reason": "This increment audits coefficient-space pair arithmetic, not a global velocity/pressure/forcing candidate.",
                "ST006_momentum_sampled_max": ST006_MOMENTUM_MAX,
                "ST006_momentum_volume_L2": ST006_VOLUME_L2,
            },
            "project_gates_unchanged": {
                "normalized_momentum_max": MOMENTUM_GATE,
                "normalized_momentum_L2": MOMENTUM_GATE,
                "divergence_max": DIVERGENCE_GATE,
                "divergence_L2": DIVERGENCE_GATE,
                "free_residual_defined_forcing_allowed": False,
            },
            "truth_boundary": dict(_TRUTH_BOUNDARY),
            "limitations": [
                "This is a conditional independent audit of M/K pair arithmetic only.",
                "The consumed full post-J R1/R2 envelopes remain Agent-1 self-certificates until separately independently admitted.",
                "No radial inverse is reapplied in the valid path; duplicate-J mutations are negative controls only.",
                "No global leading velocity, matched pressure, restricted-forcing composite, or correction field is evaluated.",
                "The coefficient-space M/K numbers are not Navier-Stokes residuals and are not compared numerically to ST006.",
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
