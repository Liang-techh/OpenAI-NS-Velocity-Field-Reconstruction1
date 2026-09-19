"""Source-compatible PA.10 ``d partial_eta A(u) Phi`` R1 slot after ``J_2``.

Pinned public provenance is KokunoYumeto/yang-mills-interacting-workbench at
``143f6773feb424ad9ed3a8d116653200f20346b7``, file
``navier-stokes/navier_stokes_workbench.tex`` (corrected reader dated
2026-09-09, Zenodo 22678406).

The corrected reconstruction expands

    W = W_* + Lambda^-1[-2D eta A(u) - d partial_eta A(u)]

inside

    L R1 = [W + h(1-2 eta U) + d u zeta_*] Phi
           + W Y Phi_Y + H_c Phi_eta.

Hence one ordinary R1 monomial is

    Lambda^-1 d (partial_eta A(u)) Phi.

The source coefficient proof gives, after the radial inverse,

    ||J_nu[(partial_eta F)(Y partial_Y G)]||_rho
        <= (80/rho) C_sq^2 ||F||_rho ||G||_rho.

It then states explicitly that if the logarithmic radial derivative is absent
its degree factor is replaced by one, and if the differentiated factor is
averaged the ``(i+1)^-1`` averaging divisor cancels the added-degree factor.
Further undifferentiated factors only add the usual convolution/product
constants.  Therefore the present slot is bounded *after* J_2 without
inventing a standalone coefficient norm for ``partial_eta A(u)``.

This module consumes the existing source-compatible radius-one ``u`` and
``Phi`` balls, fixed ``d`` and ``L^-1`` coefficient balls, and only the source
regime inequality ``Lambda^-1<=1``.  It stacks on Agent-1 #707, advances the
post-J2 R1 ordinary ledger from 6/11 to 7/11, and deliberately leaves the four
``Y Phi_Y`` / ``Phi_eta`` ordinary slots open.  It is an Agent-1 source-bound
self-certificate, not independent validation and not an NS residual.
"""

from __future__ import annotations

import argparse
import copy
from dataclasses import asdict, dataclass, field
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
from typing import Any

from .kokuno_pa10_post_j_remainder_bridge import R1_ORDINARY_TERMS
from .kokuno_pa10_remainder_ball_bounds import BallFactorBound
from .kokuno_pa10_source_r1_zeta_ordinary_slot import (
    R1_SOURCE_BOUND_THROUGH_ZETA_TERMS,
    KokunoPA10SourceR1ZetaOrdinarySlot,
)

SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
CORRECTED_RELEASE = "zenodo:22678406"
CORRECTED_RELEASE_DATE = "2026-09-09"
SCHEMA = "kokuno-pa10-source-r1-averaged-eta-ordinary-slot-v1"

R1_AVERAGED_ETA_ORDINARY_TERMS = ("lambda_inv_d_detaAu_times_Phi",)
_BOUND_NAME_SET = set(R1_SOURCE_BOUND_THROUGH_ZETA_TERMS) | set(
    R1_AVERAGED_ETA_ORDINARY_TERMS
)
R1_POST_J2_SOURCE_BOUND_THROUGH_AVERAGED_ETA_TERMS = tuple(
    name for name in R1_ORDINARY_TERMS if name in _BOUND_NAME_SET
)
R1_REMAINING_AFTER_AVERAGED_ETA_TERMS = tuple(
    name for name in R1_ORDINARY_TERMS if name not in _BOUND_NAME_SET
)

_SOURCE_FORMULAS = {
    "R1": "L R1=[W+h(1-2eta U)+d u zeta_*]Phi+W YPhi_Y+H_c Phi_eta",
    "W": "W=W_*+Lambda^-1[-2D eta A(u)-d partial_eta A(u)]",
    "ordinary_slot": "Lambda^-1 d (partial_eta A(u)) Phi",
    "coefficient_weight": (
        "a_{alpha,beta}=20^{-alpha}rho^{-beta}beta!binom(alpha+beta,beta)/"
        "((alpha+1)^2(beta+1)^2)"
    ),
    "source_post_J_derivative_bound": (
        "||J_nu[(partial_eta F)(Y partial_Y G)]||_rho <= "
        "(80/rho) C_sq^2 ||F||_rho ||G||_rho"
    ),
    "no_logradial_specialization": (
        "replace the logarithmic-radial degree factor by one"
    ),
    "averaged_specialization": (
        "A(u) contributes (i+1)^-1, cancelling the added-degree factor"
    ),
    "undifferentiated_factors": (
        "d, Phi and L^-1 enter through ordinary convolution/product bounds"
    ),
    "Lambda": "Lambda>=1, hence Lambda^-1<=1",
}

_TRUTH_BOUNDARY = {
    "public_reconstruction_source": True,
    "source_single_eta_derivative_post_J2_rule_executable": True,
    "source_averaging_divisor_cancellation_used": True,
    "source_compatible_R1_averaged_eta_ordinary_post_J2_slot_machine_bound": True,
    "source_compatible_R1_post_J2_ordinary_source_bound_count_is_seven": True,
    "standalone_deta_Au_coefficient_norm_invented": False,
    "standalone_Y_Phi_Y_coefficient_norm_invented": False,
    "standalone_Phi_eta_coefficient_norm_invented": False,
    "uniform_Lambda_inverse_upper_for_Lambda_ge_1_used": True,
    "this_R1_increment_independent_agent4_audit_required": True,
    "source_compatible_all_R1_ordinary_post_J2_slots_machine_bound": False,
    "source_full_post_J2_R1_radius_one_ball_bound_machine_bound": False,
    "full_R1_independent_agent4_admission": False,
    "source_full_post_J1_R2_independent_admission": False,
    "source_operator_constant_M_machine_bound": False,
    "source_operator_constant_K_machine_bound": False,
    "source_contraction_invariant_ball_machine_verified": False,
    "source_fixed_point_distance_machine_bound": False,
    "source_B0_dependencies_machine_bound": False,
    "source_T_sh_lower_bound_verified": False,
    "selected_pa16_handoff_allowed": False,
    "global_pressure_matched": False,
    "global_leading_profile_reconstructed": False,
    "heldout_ns_residual_assessed": False,
    "formal_full_domain_pde_gate_assessed": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


def _canonical_json(payload: dict[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _q_nonnegative(value: float, name: str) -> Fraction:
    out = float(value)
    if not math.isfinite(out) or out < 0.0:
        raise ValueError(f"{name} must be finite and nonnegative")
    return Fraction.from_float(out)


def _upper(value: Fraction) -> float:
    if value < 0:
        raise ValueError("upper-bound helper requires a nonnegative value")
    out = float(value)
    if not math.isfinite(out):
        raise OverflowError("upper bound is outside binary64 range")
    if Fraction.from_float(out) < value:
        out = math.nextafter(out, math.inf)
    return float(out)


def _product_constant_exact() -> Fraction:
    # Elementary pi < 22/7.  C_sq=4*pi^2/3 and the convolution constant is C_sq^2.
    pi_upper = Fraction(22, 7)
    c_sq_upper = Fraction(4, 3) * pi_upper * pi_upper
    return c_sq_upper * c_sq_upper


@dataclass(frozen=True)
class KokunoPA10SourceR1AveragedEtaOrdinarySlot:
    """Bind ``J2[L^-1 d (partial_eta A(u)) Phi]`` without a derivative norm."""

    parent: KokunoPA10SourceR1ZetaOrdinarySlot = field(
        default_factory=KokunoPA10SourceR1ZetaOrdinarySlot
    )

    def __post_init__(self) -> None:
        if not isinstance(self.parent, KokunoPA10SourceR1ZetaOrdinarySlot):
            raise TypeError("parent must be KokunoPA10SourceR1ZetaOrdinarySlot")
        if not 0.0 < self.rho < self.parent.domain.cauchy_radius:
            raise ValueError("coefficient rho must lie inside the Cauchy radius")

    @property
    def bridge(self):
        return self.parent.bridge

    @property
    def domain(self):
        return self.parent.domain

    @property
    def rho(self) -> float:
        return float(self.domain.coefficient_rho)

    @property
    def product_constant_upper(self) -> float:
        return _upper(_product_constant_exact())

    def source_single_eta_after_j2_factor_upper(self) -> float:
        """Return the source ``80/rho`` derivative/radial-inverse factor."""
        return _upper(Fraction(80, 1) / Fraction.from_float(self.rho))

    def operator_inputs(self) -> dict[str, BallFactorBound]:
        """Return only undifferentiated/source-ball inputs to the post-J2 rule."""
        d = self.parent.algebraic_r1.algebraic.fixed_multiplier_balls()["d"]
        u = self.parent.algebraic_r1.u_ball()
        phi = self.parent.algebraic_r1.phi_ball()
        Linv = self.bridge.l_inverse_coefficient_ball()
        return {"d": d, "u": u, "Phi": phi, "L_inverse": Linv}

    def r1_averaged_eta_after_j2(self) -> dict[str, BallFactorBound]:
        """Bound the seventh ordinary R1 slot directly after J2.

        The source post-J derivative proof already includes the radial inverse.
        Consequently this method MUST NOT call ``ordinary_term_after_jnu``:
        doing so would apply J2 twice.  Four coefficient factors are present
        after specializing the differentiated factor to averaged ``u``:
        ``u, Phi, d, L^-1``.  Their three convolutions contribute
        ``(C_sq^2)^3``.  ``d`` and ``L^-1`` are fixed; the radius-one pair
        Lipschitz bound varies only ``u`` and ``Phi`` one factor at a time.
        """
        values = self.operator_inputs()
        u = values["u"]
        phi = values["Phi"]
        d = values["d"]
        Linv = values["L_inverse"]

        rho = Fraction.from_float(self.rho)
        derivative_j2 = Fraction(80, 1) / rho
        convolution = _product_constant_exact() ** 3
        fixed = (
            convolution
            * derivative_j2
            * _q_nonnegative(d.norm, "d norm")
            * _q_nonnegative(Linv.norm, "L inverse norm")
        )
        norm = (
            fixed
            * _q_nonnegative(u.norm, "u norm")
            * _q_nonnegative(phi.norm, "Phi norm")
        )
        lip = fixed * (
            _q_nonnegative(u.lipschitz, "u lipschitz")
            * _q_nonnegative(phi.norm, "Phi norm")
            + _q_nonnegative(u.norm, "u norm")
            * _q_nonnegative(phi.lipschitz, "Phi lipschitz")
        )
        return {
            "lambda_inv_d_detaAu_times_Phi": BallFactorBound(
                _upper(norm), _upper(lip)
            )
        }

    def combined_r1_after_j2_subset(self) -> dict[str, BallFactorBound]:
        """Merge the prior six post-J2 slots with this direct post-J2 slot."""
        previous = self.parent.combined_r1_after_j2_subset()
        current = self.r1_averaged_eta_after_j2()
        overlap = set(previous).intersection(current)
        if overlap:
            raise RuntimeError(f"R1 post-J2 source-bound subset overlap: {sorted(overlap)}")
        merged = {**previous, **current}
        out = {
            name: merged[name]
            for name in R1_POST_J2_SOURCE_BOUND_THROUGH_AVERAGED_ETA_TERMS
        }
        if tuple(out) != R1_POST_J2_SOURCE_BOUND_THROUGH_AVERAGED_ETA_TERMS:
            raise RuntimeError("R1 post-J2 source-bound slot order drifted")
        if set(out) | set(R1_REMAINING_AFTER_AVERAGED_ETA_TERMS) != set(
            R1_ORDINARY_TERMS
        ):
            raise RuntimeError("R1 post-J2 ledger partition drifted")
        return out

    @property
    def truth_boundary(self) -> dict[str, bool]:
        return copy.deepcopy(_TRUTH_BOUNDARY)

    def report(self) -> dict[str, Any]:
        inputs = self.operator_inputs()
        new_slot = self.r1_averaged_eta_after_j2()
        combined = self.combined_r1_after_j2_subset()
        payload: dict[str, Any] = {
            "schema": SCHEMA,
            "source": {
                "repository": SOURCE_REPOSITORY,
                "commit": SOURCE_COMMIT,
                "path": SOURCE_PATH,
                "corrected_release": CORRECTED_RELEASE,
                "corrected_release_date": CORRECTED_RELEASE_DATE,
            },
            "source_formulas": copy.deepcopy(_SOURCE_FORMULAS),
            "upstream_identity": {
                "R1_six_slot_parent_receipt_sha256": self.parent.sha256,
                "post_J_bridge_receipt_sha256": self.bridge.sha256,
                "source_axis_receipt_sha256": self.domain.sha256,
            },
            "selected_source_compatible_inputs": {
                "rho": self.rho,
                "Lambda_inverse_uniform_upper": 1.0,
                "choice_status": (
                    "repository-autonomous source-compatible existence choice; "
                    "not recovered hidden source data"
                ),
            },
            "operator_certificate": {
                "post_J2_derivative_factor_upper": self.source_single_eta_after_j2_factor_upper(),
                "product_constant_Csq_squared_upper": self.product_constant_upper,
                "undifferentiated_convolution_count": 3,
                "averaged_factor_divisor_cancellation_used": True,
                "logradial_derivative_absent": True,
                "J2_already_in_source_derivative_factor": True,
                "ordinary_bridge_J2_reapplied": False,
            },
            "operator_inputs": {name: asdict(bound) for name, bound in inputs.items()},
            "R1_averaged_eta_ordinary_after_J2": {
                name: asdict(bound) for name, bound in new_slot.items()
            },
            "R1_combined_source_bound_after_J2": {
                name: asdict(bound) for name, bound in combined.items()
            },
            "ledger_progress": {
                "R1_post_J2_source_bound_ordinary_slots_count": len(combined),
                "R1_total_ordinary_slots": len(R1_ORDINARY_TERMS),
                "R1_ordinary_slots_still_required": len(
                    R1_REMAINING_AFTER_AVERAGED_ETA_TERMS
                ),
                "R1_remaining_ordinary_slots": list(
                    R1_REMAINING_AFTER_AVERAGED_ETA_TERMS
                ),
                "R1_numerator_space_source_bound_count_remains_six": True,
                "reason": (
                    "the new eta-derivative slot is certified only after J2; "
                    "no standalone partial_eta A(u) norm is invented"
                ),
                "R1_mixed_slot_handled_upstream_by_A1_668": True,
                "R2_all_ordinary_slots_agent1_self_certificate_exists": True,
                "full_R1_still_open": True,
                "M_K_not_promoted": True,
            },
            "independent_audit_boundary": {
                "this_averaged_eta_R1_slot_is_agent1_self_certificate": True,
                "this_averaged_eta_R1_slot_requires_independent_A4_audit": True,
                "R1_prior_six_slots_not_promoted_by_this_receipt": True,
                "full_R1_independent_A4_admission": False,
                "full_R2_independent_A4_admission_not_inferred_here": True,
            },
            "truth_boundary": self.truth_boundary,
            "limitations": [
                "This is a coefficient/operator upper bound, not an NS residual.",
                "The differentiated averaged factor is bounded only after J2; no standalone derivative norm exists.",
                "Four ordinary R1 YPhi_Y/Phi_eta slots remain open.",
                "Agent-4 independent audit is required before admission.",
                "No M/K, fixed point, B0/T_sh, PA.16, global pressure, or global leading field is promoted.",
            ],
        }
        receipt_payload = copy.deepcopy(payload)
        payload["receipt_sha256"] = hashlib.sha256(
            _canonical_json(receipt_payload).encode("utf-8")
        ).hexdigest()
        return payload

    @property
    def sha256(self) -> str:
        return str(self.report()["receipt_sha256"])

    def save(self, path: str | Path) -> Path:
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(json.dumps(self.report(), indent=2, sort_keys=True) + "\n")
        return target


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Bind the Kokuno PA.10 averaged-eta ordinary R1 slot after J2."
    )
    parser.add_argument("--output", type=Path, default=None)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _build_parser().parse_args(argv)
    calc = KokunoPA10SourceR1AveragedEtaOrdinarySlot()
    payload = calc.report()
    if args.output is not None:
        calc.save(args.output)
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
