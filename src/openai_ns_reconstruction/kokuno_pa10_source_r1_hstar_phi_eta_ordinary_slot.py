"""Source-compatible PA.10 ``H_* Phi_eta`` R1 slot after ``J_2``.

Pinned public provenance is KokunoYumeto/yang-mills-interacting-workbench at
``143f6773feb424ad9ed3a8d116653200f20346b7``, file
``navier-stokes/navier_stokes_workbench.tex`` (corrected reader dated
2026-09-09, Zenodo 22678406).

The corrected reconstruction writes

    L R1 = [W+h(1-2 eta U)+d u zeta_*] Phi
           + W Y Phi_Y + H_c Phi_eta,

with

    H_c = H_* + Lambda^-1 d u.

After Agent-1 #733, the first of the two remaining ordinary R1 monomials is

    H_* Phi_eta.

The public coefficient proof controls a single eta derivative only after the
radial inverse.  We use the same conservative source envelope already executed
on the R2 derivative route,

    ||J_2[partial_eta F]||_rho <= (80/rho) ||F||_rho.

Because ``L^-1`` and ``H_*`` depend only on eta while ``J_2`` integrates only
the radial variable, they can be multiplied around the already-post-J2
``Phi_eta`` envelope.  This never invents a standalone ``Phi_eta`` coefficient
norm and never reapplies the ordinary J2 bridge.

The fixed polynomial multiplier is the public axis profile

    H_* = D eta + (1-eta^2)(4 eta+j_0)
        = j_0 + (D+4)eta - j_0 eta^2 - 4 eta^3.

Its executable coefficient ball is reused from the already-present Agent-1 R2
derivative source module rather than duplicated here.

This module stacks on Agent-1 #733 and advances the source-compatible Agent-1
post-J2 ordinary R1 ledger from 9/11 to 10/11.  It remains a self-certificate:
independent Agent-4 admission, the last ``Lambda^-1 d u Phi_eta`` ordinary slot,
full R1, M/K, global leading velocity, and held-out Navier-Stokes validation all
remain open.
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
from .kokuno_pa10_source_r1_averaged_logradial_ordinary_slot import (
    R1_POST_J2_SOURCE_BOUND_THROUGH_AVERAGED_LOGRADIAL_TERMS,
    KokunoPA10SourceR1AveragedLogradialOrdinarySlot,
)

SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
CORRECTED_RELEASE = "zenodo:22678406"
CORRECTED_RELEASE_DATE = "2026-09-09"
SCHEMA = "kokuno-pa10-source-r1-hstar-phi-eta-ordinary-slot-v1"

R1_HSTAR_PHI_ETA_ORDINARY_TERMS = ("H_star_times_Phi_eta",)
_BOUND_NAME_SET = set(
    R1_POST_J2_SOURCE_BOUND_THROUGH_AVERAGED_LOGRADIAL_TERMS
) | set(R1_HSTAR_PHI_ETA_ORDINARY_TERMS)
R1_POST_J2_SOURCE_BOUND_THROUGH_HSTAR_PHI_ETA_TERMS = tuple(
    name for name in R1_ORDINARY_TERMS if name in _BOUND_NAME_SET
)
R1_REMAINING_AFTER_HSTAR_PHI_ETA_TERMS = tuple(
    name for name in R1_ORDINARY_TERMS if name not in _BOUND_NAME_SET
)

_SOURCE_FORMULAS = {
    "R1": "L R1=[W+h(1-2eta U)+d u zeta_*]Phi+W YPhi_Y+H_c Phi_eta",
    "H_c": "H_c=H_*+Lambda^-1 d u",
    "ordinary_slot": "H_* Phi_eta",
    "H_star": "H_*=Deta+(1-eta^2)(4eta+j_0)",
    "single_eta_rule": (
        "source single-eta specialization: "
        "||J_2[partial_eta F]||_rho <= (80/rho)||F||_rho"
    ),
    "radial_commutation": (
        "J_2[L^-1 H_* partial_eta Phi]"
        "=L^-1 H_* J_2[partial_eta Phi]"
    ),
    "product": "||FG||_rho <= C_sq^2 ||F||_rho ||G||_rho",
}

_TRUTH_BOUNDARY = {
    "public_reconstruction_source": True,
    "source_Hstar_coefficient_norm_machine_bound_reused": True,
    "source_single_eta_post_J2_rule_executable": True,
    "source_compatible_R1_Hstar_Phi_eta_ordinary_post_J2_slot_machine_bound": True,
    "source_compatible_R1_post_J2_ordinary_source_bound_count_is_ten": True,
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
class KokunoPA10SourceR1HstarPhiEtaOrdinarySlot:
    """Bind ``J2[L^-1 H_* Phi_eta]`` without a standalone eta derivative norm."""

    parent: KokunoPA10SourceR1AveragedLogradialOrdinarySlot = field(
        default_factory=KokunoPA10SourceR1AveragedLogradialOrdinarySlot
    )

    def __post_init__(self) -> None:
        if not isinstance(
            self.parent, KokunoPA10SourceR1AveragedLogradialOrdinarySlot
        ):
            raise TypeError(
                "parent must be KokunoPA10SourceR1AveragedLogradialOrdinarySlot"
            )
        if not 0.0 < self.rho < float(self.domain.cauchy_radius):
            raise ValueError("coefficient rho must lie inside the Cauchy radius")

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

    @property
    def rho(self) -> float:
        return float(self.domain.coefficient_rho)

    def source_single_eta_after_j2_factor_upper(self) -> float:
        """Return the source-common ``80/rho`` post-J2 single-eta envelope."""
        rho = Fraction.from_float(self.rho)
        return _upper(Fraction(80, 1) / rho)

    def hstar_coefficient_ball(self) -> BallFactorBound:
        """Reuse the executable public-axis ``H_*`` coefficient certificate."""
        return self.algebraic_r1.r2_source.hstar_coefficient_ball()

    def operator_inputs(self) -> dict[str, BallFactorBound]:
        return {
            "H_star": self.hstar_coefficient_ball(),
            "Phi": self.algebraic_r1.phi_ball(),
            "L_inverse": self.bridge.l_inverse_coefficient_ball(),
        }

    def j2_phi_eta(self) -> BallFactorBound:
        """Bound ``J2[Phi_eta]`` directly; no standalone ``Phi_eta`` norm."""
        return _scale_exact(
            self.operator_inputs()["Phi"],
            Fraction(80, 1) / Fraction.from_float(self.rho),
        )

    def r1_hstar_phi_eta_after_j2(self) -> dict[str, BallFactorBound]:
        """Return exactly the tenth ordinary R1 source-compatible post-J2 slot."""
        values = self.operator_inputs()
        post_j2_phi_eta = self.j2_phi_eta()
        return {
            "H_star_times_Phi_eta": self.product_calculator.product(
                values["L_inverse"],
                values["H_star"],
                post_j2_phi_eta,
            )
        }

    def combined_r1_after_j2_subset(self) -> dict[str, BallFactorBound]:
        """Merge the prior nine post-J2 slots with ``H_* Phi_eta``."""
        previous = self.parent.combined_r1_after_j2_subset()
        current = self.r1_hstar_phi_eta_after_j2()
        overlap = set(previous).intersection(current)
        if overlap:
            raise RuntimeError(
                f"R1 post-J2 source-bound subset overlap: {sorted(overlap)}"
            )
        merged = {**previous, **current}
        out = {
            name: merged[name]
            for name in R1_POST_J2_SOURCE_BOUND_THROUGH_HSTAR_PHI_ETA_TERMS
        }
        if tuple(out) != R1_POST_J2_SOURCE_BOUND_THROUGH_HSTAR_PHI_ETA_TERMS:
            raise RuntimeError("R1 post-J2 source-bound slot order drifted")
        if set(out) | set(R1_REMAINING_AFTER_HSTAR_PHI_ETA_TERMS) != set(
            R1_ORDINARY_TERMS
        ):
            raise RuntimeError("R1 post-J2 ledger partition drifted")
        return out

    @property
    def truth_boundary(self) -> dict[str, bool]:
        return copy.deepcopy(_TRUTH_BOUNDARY)

    def report(self) -> dict[str, Any]:
        inputs = self.operator_inputs()
        post_j2_phi_eta = self.j2_phi_eta()
        new_slot = self.r1_hstar_phi_eta_after_j2()
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
                "R1_nine_slot_parent_receipt_sha256": self.parent.sha256,
                "R2_derivative_source_receipt_sha256": self.algebraic_r1.r2_source.sha256,
                "post_J_bridge_receipt_sha256": self.bridge.sha256,
                "source_Phi_ball_receipt_sha256": (
                    self.bridge.u_source.mixed.phi_source.sha256
                ),
                "source_axis_receipt_sha256": self.domain.sha256,
            },
            "selected_source_compatible_inputs": {
                "rho": self.rho,
                "choice_status": (
                    "repository-autonomous source-compatible existence choice; "
                    "not recovered hidden source data"
                ),
            },
            "operator_certificate": {
                "post_J2_single_eta_factor_upper": (
                    self.source_single_eta_after_j2_factor_upper()
                ),
                "radial_inverse_commutes_with_eta_only_multipliers": True,
                "standalone_Phi_eta_norm_used": False,
                "J2_already_in_source_single_eta_factor": True,
                "ordinary_bridge_J2_reapplied": False,
                "coefficient_product_factor_count_after_J2_rule": 3,
                "Hstar_coefficient_bound_reused_not_duplicated": True,
            },
            "operator_inputs": {
                name: asdict(bound) for name, bound in inputs.items()
            },
            "J2_Phi_eta": asdict(post_j2_phi_eta),
            "R1_Hstar_Phi_eta_ordinary_after_J2": {
                name: asdict(bound) for name, bound in new_slot.items()
            },
            "R1_combined_source_bound_after_J2": {
                name: asdict(bound) for name, bound in combined.items()
            },
            "ledger_progress": {
                "R1_post_J2_source_bound_ordinary_slots_count": len(combined),
                "R1_total_ordinary_slots": len(R1_ORDINARY_TERMS),
                "R1_ordinary_slots_still_required": len(
                    R1_REMAINING_AFTER_HSTAR_PHI_ETA_TERMS
                ),
                "R1_remaining_ordinary_slots": list(
                    R1_REMAINING_AFTER_HSTAR_PHI_ETA_TERMS
                ),
                "R1_numerator_space_source_bound_count_remains_six": True,
                "reason": (
                    "this slot is certified only after J2; no standalone "
                    "Phi_eta coefficient norm is invented"
                ),
                "R1_mixed_slot_handled_upstream_by_A1_668": True,
                "R2_all_ordinary_slots_agent1_self_certificate_exists": True,
                "full_R1_still_open": True,
                "M_K_not_promoted": True,
            },
            "independent_audit_boundary": {
                "this_Hstar_Phi_eta_R1_slot_is_agent1_self_certificate": True,
                "this_Hstar_Phi_eta_R1_slot_requires_independent_A4_audit": True,
                "R1_prior_nine_slots_not_promoted_by_this_receipt": True,
                "full_R1_independent_A4_admission": False,
            },
            "truth_boundary": self.truth_boundary,
            "limitations": [
                (
                    "The ordinary R1 slot Lambda^-1 d u Phi_eta remains "
                    "source-unbound after this increment."
                ),
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
    payload = KokunoPA10SourceR1HstarPhiEtaOrdinarySlot().save_report(args.output)
    print("receipt_sha256=", payload["receipt_sha256"])
    print(
        "R1_Hstar_Phi_eta_after_J2=",
        payload["R1_Hstar_Phi_eta_ordinary_after_J2"],
    )
    print("progress=", payload["ledger_progress"])


if __name__ == "__main__":
    _main()
