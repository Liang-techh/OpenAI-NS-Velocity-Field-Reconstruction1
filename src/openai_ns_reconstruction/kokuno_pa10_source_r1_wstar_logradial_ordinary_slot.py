"""Source-compatible PA.10 ``W_* Y Phi_Y`` R1 slot after ``J_2``.

Pinned public provenance is KokunoYumeto/yang-mills-interacting-workbench at
``143f6773feb424ad9ed3a8d116653200f20346b7``, file
``navier-stokes/navier_stokes_workbench.tex`` (corrected reader dated
2026-09-09, Zenodo 22678406).

The corrected reconstruction contains

    L R1 = [W+h(1-2 eta U)+d u zeta_*] Phi
           + W Y Phi_Y + H_c Phi_eta,

with ``W=W_*+Lambda^-1[-2D eta A(u)-d partial_eta A(u)]``.  One remaining
ordinary R1 monomial is therefore

    W_* Y Phi_Y.

The source coefficient proof controls logarithmic-radial derivatives only
after the radial inverse.  In the no-eta-derivative specialization it uses the
first weight ratio; the common source constant 80 is a valid upper bound for
``J_nu[Y partial_Y F]``.  Since ``L^-1`` and ``W_*`` depend only on eta,

    J_2[L^-1 W_* Y Phi_Y]
      = J_2[Y partial_Y(L^-1 W_* Phi)].

This module uses that identity directly.  It first forms the ordinary
coefficient-space product ``L^-1 W_* Phi`` and then applies the single
logarithmic-radial post-``J_2`` factor exactly once.  It never invents a
standalone ``Y Phi_Y`` coefficient norm and never applies the ordinary J2
bridge a second time.

This stacks on Agent-1 #715 and advances the Agent-1 post-J2 ordinary R1 ledger
from 7/11 to 8/11.  It remains a source-compatible self-certificate only: no
full R1, fixed point, global leading velocity, or PDE validation is inferred.
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
from .kokuno_pa10_source_r1_averaged_eta_ordinary_slot import (
    R1_POST_J2_SOURCE_BOUND_THROUGH_AVERAGED_ETA_TERMS,
    KokunoPA10SourceR1AveragedEtaOrdinarySlot,
)

SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
CORRECTED_RELEASE = "zenodo:22678406"
CORRECTED_RELEASE_DATE = "2026-09-09"
SCHEMA = "kokuno-pa10-source-r1-wstar-logradial-ordinary-slot-v1"

R1_WSTAR_LOGRADIAL_ORDINARY_TERMS = ("W_star_times_Y_Phi_Y",)
_BOUND_NAME_SET = set(R1_POST_J2_SOURCE_BOUND_THROUGH_AVERAGED_ETA_TERMS) | set(
    R1_WSTAR_LOGRADIAL_ORDINARY_TERMS
)
R1_POST_J2_SOURCE_BOUND_THROUGH_WSTAR_LOGRADIAL_TERMS = tuple(
    name for name in R1_ORDINARY_TERMS if name in _BOUND_NAME_SET
)
R1_REMAINING_AFTER_WSTAR_LOGRADIAL_TERMS = tuple(
    name for name in R1_ORDINARY_TERMS if name not in _BOUND_NAME_SET
)

_SOURCE_FORMULAS = {
    "R1": "L R1=[W+h(1-2eta U)+d u zeta_*]Phi+W YPhi_Y+H_c Phi_eta",
    "W": "W=W_*+Lambda^-1[-2D eta A(u)-d partial_eta A(u)]",
    "ordinary_slot": "W_* Y Phi_Y",
    "commuted_post_J2_identity": (
        "J_2[L^-1 W_* YPhi_Y]=J_2[Y partial_Y(L^-1 W_* Phi)]"
    ),
    "coefficient_weight": (
        "a_{alpha,beta}=20^{-alpha}rho^{-beta}beta!binom(alpha+beta,beta)/"
        "((alpha+1)^2(beta+1)^2)"
    ),
    "single_logradial_rule": (
        "source first-weight-ratio specialization: "
        "||J_nu[Y partial_Y F]||_rho <= 80 ||F||_rho"
    ),
    "product": "||FG||_rho <= C_sq^2 ||F||_rho ||G||_rho",
}

_TRUTH_BOUNDARY = {
    "public_reconstruction_source": True,
    "source_single_logradial_post_J2_rule_executable": True,
    "source_compatible_R1_Wstar_logradial_ordinary_post_J2_slot_machine_bound": True,
    "source_compatible_R1_post_J2_ordinary_source_bound_count_is_eight": True,
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


@dataclass(frozen=True)
class KokunoPA10SourceR1WstarLogradialOrdinarySlot:
    """Bind ``J2[L^-1 W_* Y Phi_Y]`` without a standalone derivative norm."""

    parent: KokunoPA10SourceR1AveragedEtaOrdinarySlot = field(
        default_factory=KokunoPA10SourceR1AveragedEtaOrdinarySlot
    )

    def __post_init__(self) -> None:
        if not isinstance(self.parent, KokunoPA10SourceR1AveragedEtaOrdinarySlot):
            raise TypeError("parent must be KokunoPA10SourceR1AveragedEtaOrdinarySlot")

    @property
    def bridge(self):
        return self.parent.bridge

    @property
    def domain(self):
        return self.parent.domain

    @property
    def algebraic_r1(self):
        return self.parent.parent.algebraic_r1

    @property
    def product_calculator(self):
        return self.algebraic_r1.product_calculator

    def source_single_logradial_after_j2_factor_upper(self) -> float:
        """Return the source-common upper bound for ``J2[Y partial_Y F]``."""
        return _upper(Fraction(80, 1))

    def operator_inputs(self) -> dict[str, BallFactorBound]:
        """Return the eta-only fixed factors and the radius-one Phi ball."""
        fixed = self.algebraic_r1.fixed_multiplier_balls()
        return {
            "W_star": fixed["W_star"],
            "Phi": self.algebraic_r1.phi_ball(),
            "L_inverse": self.bridge.l_inverse_coefficient_ball(),
        }

    def r1_wstar_logradial_after_j2(self) -> dict[str, BallFactorBound]:
        """Bound the eighth ordinary R1 slot directly after J2.

        Because ``L^-1`` and ``W_*`` are eta-only multipliers, the target is
        exactly ``J2[Y partial_Y(L^-1 W_* Phi)]``.  We therefore apply the
        ordinary product algebra to ``L^-1 W_* Phi`` and then the source
        single-logradial post-J factor once.  No standalone ``Y Phi_Y`` norm
        and no second ordinary J2 bridge are used.
        """
        values = self.operator_inputs()
        base = self.product_calculator.product(
            values["L_inverse"], values["W_star"], values["Phi"]
        )
        return {
            "W_star_times_Y_Phi_Y": _scale_exact(
                base, Fraction(80, 1)
            )
        }

    def combined_r1_after_j2_subset(self) -> dict[str, BallFactorBound]:
        """Merge the prior seven post-J2 slots with this new logradial slot."""
        previous = self.parent.combined_r1_after_j2_subset()
        current = self.r1_wstar_logradial_after_j2()
        overlap = set(previous).intersection(current)
        if overlap:
            raise RuntimeError(f"R1 post-J2 source-bound subset overlap: {sorted(overlap)}")
        merged = {**previous, **current}
        out = {
            name: merged[name]
            for name in R1_POST_J2_SOURCE_BOUND_THROUGH_WSTAR_LOGRADIAL_TERMS
        }
        if tuple(out) != R1_POST_J2_SOURCE_BOUND_THROUGH_WSTAR_LOGRADIAL_TERMS:
            raise RuntimeError("R1 post-J2 source-bound slot order drifted")
        if set(out) | set(R1_REMAINING_AFTER_WSTAR_LOGRADIAL_TERMS) != set(
            R1_ORDINARY_TERMS
        ):
            raise RuntimeError("R1 post-J2 ledger partition drifted")
        return out

    @property
    def truth_boundary(self) -> dict[str, bool]:
        return copy.deepcopy(_TRUTH_BOUNDARY)

    def report(self) -> dict[str, Any]:
        inputs = self.operator_inputs()
        new_slot = self.r1_wstar_logradial_after_j2()
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
                "R1_seven_slot_parent_receipt_sha256": self.parent.sha256,
                "post_J_bridge_receipt_sha256": self.bridge.sha256,
                "source_axis_receipt_sha256": self.domain.sha256,
            },
            "selected_source_compatible_inputs": {
                "rho": float(self.domain.coefficient_rho),
                "choice_status": (
                    "repository-autonomous source-compatible existence choice; "
                    "not recovered hidden source data"
                ),
            },
            "operator_certificate": {
                "post_J2_single_logradial_factor_upper": (
                    self.source_single_logradial_after_j2_factor_upper()
                ),
                "eta_only_multipliers_commuted_through_Y_partial_Y": True,
                "standalone_Y_Phi_Y_norm_used": False,
                "J2_already_in_source_logradial_factor": True,
                "ordinary_bridge_J2_reapplied": False,
                "coefficient_product_factor_count_before_post_J2_rule": 3,
            },
            "operator_inputs": {name: asdict(bound) for name, bound in inputs.items()},
            "R1_Wstar_logradial_ordinary_after_J2": {
                name: asdict(bound) for name, bound in new_slot.items()
            },
            "R1_combined_source_bound_after_J2": {
                name: asdict(bound) for name, bound in combined.items()
            },
            "ledger_progress": {
                "R1_post_J2_source_bound_ordinary_slots_count": len(combined),
                "R1_total_ordinary_slots": len(R1_ORDINARY_TERMS),
                "R1_ordinary_slots_still_required": len(
                    R1_REMAINING_AFTER_WSTAR_LOGRADIAL_TERMS
                ),
                "R1_remaining_ordinary_slots": list(
                    R1_REMAINING_AFTER_WSTAR_LOGRADIAL_TERMS
                ),
                "R1_numerator_space_source_bound_count_remains_six": True,
                "reason": (
                    "the W_* YPhi_Y slot is certified only after J2; no "
                    "standalone YPhi_Y coefficient norm is invented"
                ),
                "R1_mixed_slot_handled_upstream_by_A1_668": True,
                "R2_all_ordinary_slots_agent1_self_certificate_exists": True,
                "full_R1_still_open": True,
                "M_K_not_promoted": True,
            },
            "independent_audit_boundary": {
                "this_Wstar_logradial_R1_slot_is_agent1_self_certificate": True,
                "this_Wstar_logradial_R1_slot_requires_independent_A4_audit": True,
                "R1_prior_seven_slots_not_promoted_by_this_receipt": True,
                "full_R1_independent_A4_admission": False,
            },
            "truth_boundary": self.truth_boundary,
            "limitations": [
                "Three ordinary R1 slots remain source-unbound after this increment.",
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
    payload = KokunoPA10SourceR1WstarLogradialOrdinarySlot().save_report(args.output)
    print("receipt_sha256=", payload["receipt_sha256"])
    print("R1_Wstar_logradial_after_J2=", payload["R1_Wstar_logradial_ordinary_after_J2"])
    print("progress=", payload["ledger_progress"])


if __name__ == "__main__":
    _main()
