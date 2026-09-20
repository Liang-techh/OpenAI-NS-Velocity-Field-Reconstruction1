"""Source-compatible PA.10 ``Lambda^-1 d u Phi_eta`` R1 slot after ``J_2``.

Pinned public provenance is KokunoYumeto/yang-mills-interacting-workbench at
``143f6773feb424ad9ed3a8d116653200f20346b7``, file
``navier-stokes/navier_stokes_workbench.tex`` (corrected reader dated
2026-09-09, Zenodo 22678406).

The corrected reconstruction writes

    L R1 = [W+h(1-2 eta U)+d u zeta_*] Phi
           + W Y Phi_Y + H_c Phi_eta,

with

    H_c = H_* + Lambda^-1 d u.

After Agent-1 #741 the only ordinary R1 monomial still open is therefore

    Lambda^-1 d u Phi_eta.

The public coefficient proof controls derivative-bearing monomials after the
radial inverse.  Specializing its single-eta-derivative rule to an
undifferentiated factor ``u`` gives the conservative source envelope

    ||J_2[u partial_eta Phi]||_rho
        <= (80/rho) C_sq^2 ||u||_rho ||Phi||_rho.

The eta-only factors ``L^-1`` and ``d`` commute with the radial inverse and are
then added by the ordinary coefficient-product bound.  No standalone
``Phi_eta`` coefficient norm is invented and the ordinary J2 bridge is not
reapplied.  ``Lambda^-1`` is bounded only by the public regime
``Lambda>=1 => Lambda^-1<=1``.

This closes the Agent-1 source-compatible ordinary R1 post-J2 ledger at 11/11
and materializes the corresponding Agent-1 full-R1 self-envelope by adding the
already-existing mixed seam exactly once.  Independent Agent-4 admission,
M/K, fixed-point closure, global leading velocity, and held-out NS validation
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
from .kokuno_pa10_source_r1_hstar_phi_eta_ordinary_slot import (
    R1_POST_J2_SOURCE_BOUND_THROUGH_HSTAR_PHI_ETA_TERMS,
    KokunoPA10SourceR1HstarPhiEtaOrdinarySlot,
)

SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
CORRECTED_RELEASE = "zenodo:22678406"
CORRECTED_RELEASE_DATE = "2026-09-09"
SCHEMA = "kokuno-pa10-source-r1-du-phi-eta-ordinary-slot-v1"

R1_DU_PHI_ETA_ORDINARY_TERMS = ("lambda_inv_d_u_times_Phi_eta",)
_BOUND_NAME_SET = set(
    R1_POST_J2_SOURCE_BOUND_THROUGH_HSTAR_PHI_ETA_TERMS
) | set(R1_DU_PHI_ETA_ORDINARY_TERMS)
R1_ALL_ORDINARY_POST_J2_SOURCE_BOUND_TERMS = tuple(
    name for name in R1_ORDINARY_TERMS if name in _BOUND_NAME_SET
)
R1_REMAINING_AFTER_DU_PHI_ETA_TERMS = tuple(
    name for name in R1_ORDINARY_TERMS if name not in _BOUND_NAME_SET
)

_SOURCE_FORMULAS = {
    "R1": "L R1=[W+h(1-2eta U)+d u zeta_*]Phi+W YPhi_Y+H_c Phi_eta",
    "H_c": "H_c=H_*+Lambda^-1 d u",
    "ordinary_slot": "Lambda^-1 d u Phi_eta",
    "single_eta_product_rule": (
        "single-eta specialization after J_2: "
        "||J_2[u partial_eta Phi]||_rho <= "
        "(80/rho) C_sq^2 ||u||_rho ||Phi||_rho"
    ),
    "radial_commutation": (
        "J_2[L^-1 d u partial_eta Phi]"
        "=L^-1 d J_2[u partial_eta Phi]"
    ),
    "product": "||FG||_rho <= C_sq^2 ||F||_rho ||G||_rho",
    "Lambda": "Lambda>=1, hence Lambda^-1<=1",
}

_TRUTH_BOUNDARY = {
    "public_reconstruction_source": True,
    "source_single_eta_product_post_J2_rule_executable": True,
    "source_compatible_R1_du_Phi_eta_ordinary_post_J2_slot_machine_bound": True,
    "source_compatible_all_R1_ordinary_post_J2_slots_machine_bound": True,
    "source_compatible_full_post_J2_R1_radius_one_ball_bound_agent1_self_certificate": True,
    "source_compatible_R1_post_J2_ordinary_source_bound_count_is_eleven": True,
    "standalone_Phi_eta_coefficient_norm_invented": False,
    "ordinary_bridge_J2_reapplied": False,
    "uniform_Lambda_inverse_upper_for_Lambda_ge_1_used": True,
    "mixed_R1_slot_added_exactly_once": True,
    "this_R1_increment_independent_agent4_audit_required": True,
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
    norm = sum(
        (_q_nonnegative(item.norm, "bound norm") for item in bounds),
        Fraction(0),
    )
    lip = sum(
        (_q_nonnegative(item.lipschitz, "bound lipschitz") for item in bounds),
        Fraction(0),
    )
    return BallFactorBound(_upper(norm), _upper(lip))


@dataclass(frozen=True)
class KokunoPA10SourceR1DuPhiEtaOrdinarySlot:
    """Bind the final ordinary R1 slot directly after J2."""

    parent: KokunoPA10SourceR1HstarPhiEtaOrdinarySlot = field(
        default_factory=KokunoPA10SourceR1HstarPhiEtaOrdinarySlot
    )

    def __post_init__(self) -> None:
        if not isinstance(self.parent, KokunoPA10SourceR1HstarPhiEtaOrdinarySlot):
            raise TypeError("parent must be KokunoPA10SourceR1HstarPhiEtaOrdinarySlot")
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
        return self.parent.product_calculator

    @property
    def rho(self) -> float:
        return float(self.domain.coefficient_rho)

    def source_single_eta_after_j2_factor_upper(self) -> float:
        """Return the source-common ``80/rho`` post-J2 single-eta factor."""
        return _upper(Fraction(80, 1) / Fraction.from_float(self.rho))

    def operator_inputs(self) -> dict[str, BallFactorBound]:
        fixed = self.algebraic_r1.algebraic.fixed_multiplier_balls()
        return {
            "d": fixed["d"],
            "u": self.algebraic_r1.u_ball(),
            "Phi": self.algebraic_r1.phi_ball(),
            "L_inverse": self.bridge.l_inverse_coefficient_ball(),
        }

    def j2_u_times_phi_eta(self) -> BallFactorBound:
        """Bound ``J2[u Phi_eta]`` without a standalone ``Phi_eta`` norm."""
        values = self.operator_inputs()
        product = self.product_calculator.product(values["u"], values["Phi"])
        return _scale_exact(
            product,
            Fraction(80, 1) / Fraction.from_float(self.rho),
        )

    def r1_du_phi_eta_after_j2(self) -> dict[str, BallFactorBound]:
        """Return exactly the eleventh ordinary R1 source-compatible slot."""
        values = self.operator_inputs()
        core = self.j2_u_times_phi_eta()
        return {
            "lambda_inv_d_u_times_Phi_eta": self.product_calculator.product(
                values["L_inverse"],
                values["d"],
                core,
            )
        }

    def all_r1_ordinary_after_j2(self) -> dict[str, BallFactorBound]:
        """Merge the prior ten slots with the final ordinary post-J2 slot."""
        previous = self.parent.combined_r1_after_j2_subset()
        current = self.r1_du_phi_eta_after_j2()
        overlap = set(previous).intersection(current)
        if overlap:
            raise RuntimeError(
                f"R1 post-J2 source-bound subset overlap: {sorted(overlap)}"
            )
        merged = {**previous, **current}
        if set(merged) != set(R1_ORDINARY_TERMS):
            missing = sorted(set(R1_ORDINARY_TERMS) - set(merged))
            extra = sorted(set(merged) - set(R1_ORDINARY_TERMS))
            raise RuntimeError(
                f"R1 post-J2 ordinary ledger mismatch; missing={missing}, extra={extra}"
            )
        out = {name: merged[name] for name in R1_ORDINARY_TERMS}
        if tuple(out) != R1_ALL_ORDINARY_POST_J2_SOURCE_BOUND_TERMS:
            raise RuntimeError("R1 post-J2 ordinary slot order drifted")
        if R1_REMAINING_AFTER_DU_PHI_ETA_TERMS:
            raise RuntimeError("R1 post-J2 ordinary ledger did not close")
        return out

    def full_r1_after_j2(self) -> BallFactorBound:
        """Add all eleven ordinary slots and the existing mixed seam once."""
        ordinary = self.all_r1_ordinary_after_j2()
        mixed = self.bridge.full_mixed_terms()[
            "lambda_inv_Linv_d_detaAu_Y_Phi_Y_after_J2"
        ]
        return _add_exact(*ordinary.values(), mixed)

    @property
    def truth_boundary(self) -> dict[str, bool]:
        return copy.deepcopy(_TRUTH_BOUNDARY)

    def report(self) -> dict[str, Any]:
        inputs = self.operator_inputs()
        core = self.j2_u_times_phi_eta()
        new_slot = self.r1_du_phi_eta_after_j2()
        all_ordinary = self.all_r1_ordinary_after_j2()
        mixed = self.bridge.full_mixed_terms()[
            "lambda_inv_Linv_d_detaAu_Y_Phi_Y_after_J2"
        ]
        full = self.full_r1_after_j2()
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
                "R1_ten_slot_parent_receipt_sha256": self.parent.sha256,
                "post_J_bridge_receipt_sha256": self.bridge.sha256,
                "source_u_ball_receipt_sha256": self.bridge.u_source.sha256,
                "source_Phi_ball_receipt_sha256": (
                    self.bridge.u_source.mixed.phi_source.sha256
                ),
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
                "post_J2_single_eta_factor_upper": (
                    self.source_single_eta_after_j2_factor_upper()
                ),
                "single_undifferentiated_radial_factor_inside_J2_rule": "u",
                "eta_only_multipliers_outside_J2_rule": ["L_inverse", "d"],
                "standalone_Phi_eta_norm_used": False,
                "J2_already_in_source_single_eta_product_factor": True,
                "ordinary_bridge_J2_reapplied": False,
                "Lambda_inverse_uniform_upper": 1.0,
                "mixed_R1_slot_added_exactly_once_in_full_envelope": True,
            },
            "operator_inputs": {
                name: asdict(bound) for name, bound in inputs.items()
            },
            "J2_u_times_Phi_eta": asdict(core),
            "R1_du_Phi_eta_ordinary_after_J2": {
                name: asdict(bound) for name, bound in new_slot.items()
            },
            "R1_all_ordinary_after_J2": {
                name: asdict(bound) for name, bound in all_ordinary.items()
            },
            "R1_mixed_after_J2": asdict(mixed),
            "R1_full_after_J2_agent1_self_certificate": asdict(full),
            "ledger_progress": {
                "R1_post_J2_source_bound_ordinary_slots_count": len(all_ordinary),
                "R1_total_ordinary_slots": len(R1_ORDINARY_TERMS),
                "R1_ordinary_slots_still_required": 0,
                "R1_remaining_ordinary_slots": [],
                "R1_numerator_space_source_bound_count_remains_six": True,
                "R1_mixed_slot_handled_upstream_by_A1_668": True,
                "full_R1_agent1_self_certificate_materialized": True,
                "full_R1_independent_admission_still_open": True,
                "R2_all_ordinary_slots_agent1_self_certificate_exists": True,
                "M_K_not_promoted": True,
            },
            "independent_audit_boundary": {
                "this_du_Phi_eta_R1_slot_is_agent1_self_certificate": True,
                "this_du_Phi_eta_R1_slot_requires_independent_A4_audit": True,
                "R1_prior_ten_slots_not_promoted_by_this_receipt": True,
                "full_R1_agent1_self_certificate_is_not_independent_admission": True,
                "full_R1_independent_A4_admission": False,
            },
            "truth_boundary": self.truth_boundary,
            "limitations": [
                (
                    "All eleven ordinary R1 slots now have Agent-1 source-compatible "
                    "post-J2 bounds, but independent Agent-4 admission is still open."
                ),
                (
                    "The full R1 value here is only the sum of those Agent-1 "
                    "self-bounds plus the already-existing mixed seam."
                ),
                "This coefficient/operator envelope is not a Navier-Stokes residual.",
                "M/K and fixed-point closure are not computed or promoted here.",
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
    payload = KokunoPA10SourceR1DuPhiEtaOrdinarySlot().save_report(args.output)
    print("receipt_sha256=", payload["receipt_sha256"])
    print(
        "R1_du_Phi_eta_after_J2=",
        payload["R1_du_Phi_eta_ordinary_after_J2"],
    )
    print(
        "R1_full_after_J2_agent1_self_certificate=",
        payload["R1_full_after_J2_agent1_self_certificate"],
    )
    print("progress=", payload["ledger_progress"])


if __name__ == "__main__":
    _main()
