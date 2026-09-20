"""Source-compatible PA.10 ``A(u) Y Phi_Y`` R1 slot after ``J_2``.

Pinned public provenance is KokunoYumeto/yang-mills-interacting-workbench at
``143f6773feb424ad9ed3a8d116653200f20346b7``, file
``navier-stokes/navier_stokes_workbench.tex`` (corrected reader dated
2026-09-09, Zenodo 22678406).

The corrected reconstruction writes

    L R1 = [W+h(1-2 eta U)+d u zeta_*] Phi
           + W Y Phi_Y + H_c Phi_eta,

with

    W = W_* + Lambda^-1[-2D eta A(u)-d partial_eta A(u)].

After Agent-1 #724, one of the three remaining ordinary R1 monomials is

    Lambda^-1 2D eta A(u) Y Phi_Y.

The source averaging operator satisfies

    Y partial_Y A(u) = u - A(u).

Hence

    A(u) Y Phi_Y
      = Y partial_Y(A(u) Phi) - u Phi + A(u) Phi.

This identity is useful because the source coefficient proof controls a single
logarithmic-radial derivative only *after* the radial inverse.  It lets us
bound the target without inventing a standalone ``Y Phi_Y`` coefficient norm:
apply the source post-J2 logradial factor to ``A(u) Phi`` and the ordinary J2
factor to the two algebraic correction terms.  The eta-only multipliers
``L^-1`` and ``eta`` commute with the radial integration/logarithmic derivative
and are multiplied afterward.  ``Lambda^-1<=1`` is used uniformly for the
public source regime ``Lambda>=1``; no hidden numerical Lambda is selected.

This module stacks on Agent-1 #724 and advances the source-compatible Agent-1
post-J2 ordinary R1 ledger from 8/11 to 9/11.  It remains a self-certificate:
independent Agent-4 admission, full R1, M/K, the global leading velocity, and
held-out Navier-Stokes validation all remain open.
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
from .kokuno_pa10_source_r1_wstar_logradial_ordinary_slot import (
    R1_POST_J2_SOURCE_BOUND_THROUGH_WSTAR_LOGRADIAL_TERMS,
    KokunoPA10SourceR1WstarLogradialOrdinarySlot,
)

SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
CORRECTED_RELEASE = "zenodo:22678406"
CORRECTED_RELEASE_DATE = "2026-09-09"
SCHEMA = "kokuno-pa10-source-r1-averaged-logradial-ordinary-slot-v1"

R1_AVERAGED_LOGRADIAL_ORDINARY_TERMS = (
    "lambda_inv_2D_eta_Au_times_Y_Phi_Y",
)
_BOUND_NAME_SET = set(R1_POST_J2_SOURCE_BOUND_THROUGH_WSTAR_LOGRADIAL_TERMS) | set(
    R1_AVERAGED_LOGRADIAL_ORDINARY_TERMS
)
R1_POST_J2_SOURCE_BOUND_THROUGH_AVERAGED_LOGRADIAL_TERMS = tuple(
    name for name in R1_ORDINARY_TERMS if name in _BOUND_NAME_SET
)
R1_REMAINING_AFTER_AVERAGED_LOGRADIAL_TERMS = tuple(
    name for name in R1_ORDINARY_TERMS if name not in _BOUND_NAME_SET
)

_SOURCE_FORMULAS = {
    "R1": "L R1=[W+h(1-2eta U)+d u zeta_*]Phi+W YPhi_Y+H_c Phi_eta",
    "W": "W=W_*+Lambda^-1[-2D eta A(u)-d partial_eta A(u)]",
    "ordinary_slot": "Lambda^-1 2D eta A(u) Y Phi_Y",
    "average": "A(u)=Y^-1 Iu; Y partial_Y A(u)=u-A(u)",
    "decomposition": (
        "A(u)Y Phi_Y=Y partial_Y(A(u)Phi)-u Phi+A(u)Phi"
    ),
    "single_logradial_rule": (
        "source first-weight-ratio specialization: "
        "||J_2[Y partial_Y F]||_rho <= 80 ||F||_rho"
    ),
    "plain_J2_rule": "source vanishing-order b=0: ||J_2 F||_rho <= 40 ||F||_rho",
    "product": "||FG||_rho <= C_sq^2 ||F||_rho ||G||_rho",
    "Lambda": "Lambda>=1, hence Lambda^-1<=1",
}

_TRUTH_BOUNDARY = {
    "public_reconstruction_source": True,
    "source_averaging_identity_executable": True,
    "source_single_logradial_post_J2_rule_executable": True,
    "source_plain_post_J2_rule_executable": True,
    "source_compatible_R1_averaged_logradial_ordinary_post_J2_slot_machine_bound": True,
    "source_compatible_R1_post_J2_ordinary_source_bound_count_is_nine": True,
    "uniform_Lambda_inverse_upper_for_Lambda_ge_1_used": True,
    "standalone_Y_Phi_Y_coefficient_norm_invented": False,
    "standalone_Phi_eta_coefficient_norm_invented": False,
    "ordinary_bridge_J2_reapplied": False,
    "this_R1_increment_independent_agent4_audit_required": True,
    "source_compatible_all_R1_ordinary_post_J2_slots_machine_bound": False,
    "source_full_post_J2_R1_radius_one_ball_bound_machine_bound": False,
    "full_R1_independent_agent4_admission": False,
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


def _scale_exact(bound: BallFactorBound, scalar: Fraction) -> BallFactorBound:
    if not isinstance(bound, BallFactorBound):
        raise TypeError("bound must be a BallFactorBound")
    if scalar < 0:
        scalar = -scalar
    return BallFactorBound(
        _upper(scalar * _q_nonnegative(bound.norm, "bound norm")),
        _upper(scalar * _q_nonnegative(bound.lipschitz, "bound lipschitz")),
    )


def _add_exact(*bounds: BallFactorBound) -> BallFactorBound:
    if not bounds:
        return BallFactorBound(0.0, 0.0)
    if not all(isinstance(item, BallFactorBound) for item in bounds):
        raise TypeError("all bounds must be BallFactorBound values")
    norm = sum((_q_nonnegative(item.norm, "norm") for item in bounds), Fraction(0))
    lip = sum(
        (_q_nonnegative(item.lipschitz, "lipschitz") for item in bounds),
        Fraction(0),
    )
    return BallFactorBound(_upper(norm), _upper(lip))


@dataclass(frozen=True)
class KokunoPA10SourceR1AveragedLogradialOrdinarySlot:
    """Bind ``J2[L^-1 2D eta A(u) Y Phi_Y]`` after radial inversion."""

    parent: KokunoPA10SourceR1WstarLogradialOrdinarySlot = field(
        default_factory=KokunoPA10SourceR1WstarLogradialOrdinarySlot
    )

    def __post_init__(self) -> None:
        if not isinstance(self.parent, KokunoPA10SourceR1WstarLogradialOrdinarySlot):
            raise TypeError("parent must be KokunoPA10SourceR1WstarLogradialOrdinarySlot")

    @property
    def bridge(self):
        return self.parent.bridge

    @property
    def domain(self):
        return self.parent.domain

    @property
    def algebraic_r1(self):
        return self.parent.algebraic_r1

    @property
    def product_calculator(self):
        return self.algebraic_r1.product_calculator

    def source_single_logradial_after_j2_factor_upper(self) -> float:
        return _upper(Fraction(80, 1))

    def source_plain_after_j2_factor_upper(self) -> float:
        return _upper(Fraction(40, 1))

    def operator_inputs(self) -> dict[str, BallFactorBound]:
        fixed = self.algebraic_r1.fixed_multiplier_balls()
        return {
            "eta": fixed["eta"],
            "A_u": self.algebraic_r1.radial_average_u_ball(),
            "u": self.algebraic_r1.u_ball(),
            "Phi": self.algebraic_r1.phi_ball(),
            "L_inverse": self.bridge.l_inverse_coefficient_ball(),
        }

    def j2_average_times_logradial_phi(self) -> BallFactorBound:
        """Bound ``J2[A(u) Y Phi_Y]`` through the source averaging identity.

        ``A(u)Y Phi_Y = Y d_Y(A(u)Phi) - u Phi + A(u)Phi``.  Signs are
        irrelevant for the norm/Lipschitz envelope, so the triangle bound uses
        one post-J2 logradial term plus two ordinary post-J2 terms.  No
        standalone derivative norm is introduced.
        """
        values = self.operator_inputs()
        Au_phi = self.product_calculator.product(values["A_u"], values["Phi"])
        u_phi = self.product_calculator.product(values["u"], values["Phi"])
        return _add_exact(
            _scale_exact(Au_phi, Fraction(80, 1)),
            _scale_exact(u_phi, Fraction(40, 1)),
            _scale_exact(Au_phi, Fraction(40, 1)),
        )

    def r1_averaged_logradial_after_j2(self) -> dict[str, BallFactorBound]:
        """Return exactly the ninth ordinary R1 post-J2 source-compatible slot."""
        values = self.operator_inputs()
        core_after_j2 = self.j2_average_times_logradial_phi()
        with_eta_linv = self.product_calculator.product(
            values["L_inverse"], values["eta"], core_after_j2
        )
        h = Fraction.from_float(float(self.domain.h))
        D = Fraction(1, 2) - h
        if h < 0 or D <= 0:
            raise ValueError("source h/D parameters are outside their stated regime")
        return {
            "lambda_inv_2D_eta_Au_times_Y_Phi_Y": _scale_exact(
                with_eta_linv, 2 * D
            )
        }

    def combined_r1_after_j2_subset(self) -> dict[str, BallFactorBound]:
        previous = self.parent.combined_r1_after_j2_subset()
        current = self.r1_averaged_logradial_after_j2()
        overlap = set(previous).intersection(current)
        if overlap:
            raise RuntimeError(f"R1 post-J2 source-bound subset overlap: {sorted(overlap)}")
        merged = {**previous, **current}
        out = {
            name: merged[name]
            for name in R1_POST_J2_SOURCE_BOUND_THROUGH_AVERAGED_LOGRADIAL_TERMS
        }
        if tuple(out) != R1_POST_J2_SOURCE_BOUND_THROUGH_AVERAGED_LOGRADIAL_TERMS:
            raise RuntimeError("R1 post-J2 source-bound slot order drifted")
        if set(out) | set(R1_REMAINING_AFTER_AVERAGED_LOGRADIAL_TERMS) != set(
            R1_ORDINARY_TERMS
        ):
            raise RuntimeError("R1 post-J2 ledger partition drifted")
        return out

    @property
    def truth_boundary(self) -> dict[str, bool]:
        return copy.deepcopy(_TRUTH_BOUNDARY)

    def report(self) -> dict[str, Any]:
        inputs = self.operator_inputs()
        core = self.j2_average_times_logradial_phi()
        new_slot = self.r1_averaged_logradial_after_j2()
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
                "R1_eight_slot_parent_receipt_sha256": self.parent.sha256,
                "post_J_bridge_receipt_sha256": self.bridge.sha256,
                "source_u_ball_receipt_sha256": self.bridge.u_source.sha256,
                "source_Phi_ball_receipt_sha256": self.bridge.u_source.mixed.phi_source.sha256,
                "source_axis_receipt_sha256": self.domain.sha256,
            },
            "selected_source_compatible_inputs": {
                "h": float(self.domain.h),
                "D": float(self.domain.D),
                "rho": float(self.domain.coefficient_rho),
                "Lambda_inverse_upper": 1.0,
                "choice_status": (
                    "repository-autonomous source-compatible existence choice; "
                    "not recovered hidden source data"
                ),
            },
            "operator_certificate": {
                "post_J2_single_logradial_factor_upper": (
                    self.source_single_logradial_after_j2_factor_upper()
                ),
                "post_J2_plain_factor_upper": self.source_plain_after_j2_factor_upper(),
                "averaging_identity_used": (
                    "A(u)Y Phi_Y=Y partial_Y(A(u)Phi)-u Phi+A(u)Phi"
                ),
                "eta_and_L_inverse_commute_with_radial_operator": True,
                "standalone_Y_Phi_Y_norm_used": False,
                "J2_already_applied_in_decomposition_terms": True,
                "ordinary_bridge_J2_reapplied": False,
                "decomposition_terms_after_J2": 3,
            },
            "operator_inputs": {name: asdict(bound) for name, bound in inputs.items()},
            "J2_Au_times_Y_Phi_Y": asdict(core),
            "R1_averaged_logradial_ordinary_after_J2": {
                name: asdict(bound) for name, bound in new_slot.items()
            },
            "R1_combined_source_bound_after_J2": {
                name: asdict(bound) for name, bound in combined.items()
            },
            "ledger_progress": {
                "R1_post_J2_source_bound_ordinary_slots_count": len(combined),
                "R1_total_ordinary_slots": len(R1_ORDINARY_TERMS),
                "R1_ordinary_slots_still_required": len(
                    R1_REMAINING_AFTER_AVERAGED_LOGRADIAL_TERMS
                ),
                "R1_remaining_ordinary_slots": list(
                    R1_REMAINING_AFTER_AVERAGED_LOGRADIAL_TERMS
                ),
                "R1_numerator_space_source_bound_count_remains_six": True,
                "reason": (
                    "this slot is certified only after J2; no standalone "
                    "Y Phi_Y coefficient norm is invented"
                ),
                "R1_mixed_slot_handled_upstream_by_A1_668": True,
                "R2_all_ordinary_slots_agent1_self_certificate_exists": True,
                "full_R1_still_open": True,
                "M_K_not_promoted": True,
            },
            "independent_audit_boundary": {
                "this_averaged_logradial_R1_slot_is_agent1_self_certificate": True,
                "this_averaged_logradial_R1_slot_requires_independent_A4_audit": True,
                "R1_prior_eight_slots_not_promoted_by_this_receipt": True,
                "full_R1_independent_A4_admission": False,
            },
            "truth_boundary": self.truth_boundary,
            "limitations": [
                "Two ordinary R1 Phi_eta slots remain source-unbound after this increment.",
                "This coefficient/operator envelope is not a Navier-Stokes residual.",
                "No hidden Kokuno/OpenAI parameter is recovered or claimed paper-exact.",
            ],
        }
        payload["receipt_sha256"] = hashlib.sha256(
            _canonical_json(payload).encode("utf-8")
        ).hexdigest()
        return payload

    @property
    def sha256(self) -> str:
        return str(self.report()["receipt_sha256"])

    def save_report(self, path: str | Path) -> dict[str, Any]:
        payload = self.report()
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        return payload


def _main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    payload = KokunoPA10SourceR1AveragedLogradialOrdinarySlot().save_report(args.output)
    print("receipt_sha256=", payload["receipt_sha256"])
    print(
        "R1_averaged_logradial_after_J2=",
        payload["R1_averaged_logradial_ordinary_after_J2"],
    )
    print("progress=", payload["ledger_progress"])


if __name__ == "__main__":
    _main()
