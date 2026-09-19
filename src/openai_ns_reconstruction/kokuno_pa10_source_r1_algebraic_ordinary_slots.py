"""Source-compatible algebraic ordinary PA.10 R1 slots after ``J_2``.

Pinned public provenance is KokunoYumeto/yang-mills-interacting-workbench at
``143f6773feb424ad9ed3a8d116653200f20346b7``, file
``navier-stokes/navier_stokes_workbench.tex`` (corrected reader dated
2026-09-09, Zenodo 22678406).

The corrected reconstruction writes

    L R1=[W+h(1-2 eta U)+d u zeta_*] Phi
         + W Y Phi_Y + H_c Phi_eta,

with

    U=U_*+Lambda^-1 u,
    W=W_*+Lambda^-1[-2D eta A(u)-d partial_eta A(u)],
    H_c=H_*+Lambda^-1 d u.

After expansion, five ordinary R1 monomials require only already-materialized
coefficient balls and no unintegrated eta/log-radial derivative:

    W_* Phi,
    Lambda^-1 2D eta A(u) Phi,
    h Phi,
    2h eta U_* Phi,
    2h Lambda^-1 eta u Phi.

This module binds exactly those five numerator slots and applies the existing
``L^-1`` + ``J_2`` bridge exactly once.  It consumes the source-compatible Phi
radius-one ball, the independently admitted source-compatible u ball, fixed
axis multipliers, and only the uniform source-regime inequality
``Lambda^-1<=1``.  It does not invent standalone ``Y Phi_Y`` or ``Phi_eta``
coefficient norms and deliberately leaves the remaining six ordinary R1 slots
open.

The selected sigma/rho/pressure datum remains a repository-autonomous
source-compatible existence choice, not recovered hidden Kokuno/OpenAI data.
This is an Agent-1 source-bound subset, not independent validation and not an
NS residual.
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
from .kokuno_pa10_source_derivative_ordinary_slots import (
    KokunoPA10SourceDerivativeOrdinarySlots,
)

SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
CORRECTED_RELEASE = "zenodo:22678406"
CORRECTED_RELEASE_DATE = "2026-09-09"
SCHEMA = "kokuno-pa10-source-r1-algebraic-ordinary-slots-v1"

R1_ALGEBRAIC_ORDINARY_TERMS = (
    "W_star_times_Phi",
    "lambda_inv_2D_eta_Au_times_Phi",
    "h_times_Phi",
    "2h_eta_Ustar_times_Phi",
    "2h_lambda_inv_eta_u_times_Phi",
)

R1_REMAINING_ORDINARY_TERMS = tuple(
    name for name in R1_ORDINARY_TERMS if name not in R1_ALGEBRAIC_ORDINARY_TERMS
)

_SOURCE_FORMULAS = {
    "R1": (
        "L R1=[W+h(1-2eta U)+d u zeta_*]Phi+W YPhi_Y+H_c Phi_eta"
    ),
    "rescaling": "U=U_*+Lambda^-1u",
    "W": "W=W_*+Lambda^-1[-2D eta A(u)-d partial_eta A(u)]",
    "H_c": "H_c=H_*+Lambda^-1 d u",
    "algebraic_subset": (
        "W_*Phi - Lambda^-1 2D eta A(u)Phi + hPhi "
        "-2h eta U_*Phi -2h Lambda^-1 eta u Phi"
    ),
    "average": "A(u)=Y^-1 Iu; ||A(u)||_rho<=||u||_rho",
    "Lambda": "Lambda>=1, hence Lambda^-1<=1",
    "post_J2": "J_2[L^-1 R1_numerator]",
}

_TRUTH_BOUNDARY = {
    "public_reconstruction_source": True,
    "source_compatible_R1_algebraic_ordinary_numerator_slots_machine_bound": True,
    "source_compatible_R1_algebraic_ordinary_post_J2_slots_machine_bound": True,
    "source_compatible_Phi_radius_one_ball_consumed": True,
    "source_compatible_u_radius_one_ball_consumed": True,
    "source_radial_average_contraction_consumed": True,
    "uniform_Lambda_inverse_upper_for_Lambda_ge_1_used": True,
    "standalone_Y_Phi_Y_coefficient_norm_invented": False,
    "standalone_Phi_eta_coefficient_norm_invented": False,
    "this_R1_subset_independent_agent4_audit_required": True,
    "source_compatible_all_R1_ordinary_post_J2_slots_machine_bound": False,
    "source_full_post_J2_R1_radius_one_ball_bound_machine_bound": False,
    "full_R1_independent_agent4_admission": False,
    "full_R2_independent_agent4_admission": False,
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
class KokunoPA10SourceR1AlgebraicOrdinarySlots:
    """Bind the five derivative-free/algebraic ordinary PA.10 R1 slots."""

    r2_source: KokunoPA10SourceDerivativeOrdinarySlots = field(
        default_factory=KokunoPA10SourceDerivativeOrdinarySlots
    )

    def __post_init__(self) -> None:
        if not isinstance(self.r2_source, KokunoPA10SourceDerivativeOrdinarySlots):
            raise TypeError("r2_source must be KokunoPA10SourceDerivativeOrdinarySlots")
        if not 0.0 < self.domain.coefficient_rho < self.domain.cauchy_radius:
            raise ValueError("coefficient rho must lie inside the Cauchy radius")

    @property
    def bridge(self):
        return self.r2_source.bridge

    @property
    def algebraic(self):
        return self.r2_source.algebraic

    @property
    def domain(self):
        return self.r2_source.domain

    @property
    def product_calculator(self):
        return self.algebraic.product_calculator

    def phi_ball(self) -> BallFactorBound:
        """Return the existing source-compatible radius-one Phi envelope."""
        return self.bridge.u_source.mixed.source_phi_ball()

    def u_ball(self) -> BallFactorBound:
        return self.bridge.u_source.u_ball()

    def radial_average_u_ball(self) -> BallFactorBound:
        """Consume the source averaging contraction already used by the R2 route."""
        return self.r2_source.radial_average_u_ball()

    def fixed_multiplier_balls(self) -> dict[str, BallFactorBound]:
        fixed = self.algebraic.fixed_multiplier_balls()
        return {
            "eta": fixed["eta"],
            "U_star": fixed["U_star"],
            "W_star": self.r2_source.wstar_coefficient_ball(),
            "L_inverse": self.bridge.l_inverse_coefficient_ball(),
        }

    def r1_algebraic_numerator_slots(self) -> dict[str, BallFactorBound]:
        """Return exactly five source-compatible ordinary R1 numerator slots."""
        fixed = self.fixed_multiplier_balls()
        phi = self.phi_ball()
        u = self.u_ball()
        Au = self.radial_average_u_ball()
        calc = self.product_calculator

        h = Fraction.from_float(float(self.domain.h))
        D = Fraction(1, 2) - h
        if h < 0 or D <= 0:
            raise ValueError("source h/D parameters are outside their stated regime")

        eta_Au_phi = calc.product(fixed["eta"], Au, phi)
        eta_Ustar_phi = calc.product(fixed["eta"], fixed["U_star"], phi)
        eta_u_phi = calc.product(fixed["eta"], u, phi)

        return {
            "W_star_times_Phi": calc.product(fixed["W_star"], phi),
            # Source regime Lambda>=1 is used only as Lambda^-1<=1.
            "lambda_inv_2D_eta_Au_times_Phi": _scale_exact(
                eta_Au_phi, 2 * D
            ),
            "h_times_Phi": _scale_exact(phi, h),
            "2h_eta_Ustar_times_Phi": _scale_exact(
                eta_Ustar_phi, 2 * h
            ),
            "2h_lambda_inv_eta_u_times_Phi": _scale_exact(
                eta_u_phi, 2 * h
            ),
        }

    def r1_algebraic_after_j2(self) -> dict[str, BallFactorBound]:
        """Apply the fixed L^-1 multiplier and J2 exactly once."""
        slots = self.r1_algebraic_numerator_slots()
        return {
            name: self.bridge.ordinary_term_after_jnu(bound, nu=2)
            for name, bound in slots.items()
        }

    def combined_r1_source_bound_subset(self) -> dict[str, BallFactorBound]:
        out = self.r1_algebraic_numerator_slots()
        if tuple(out) != R1_ALGEBRAIC_ORDINARY_TERMS:
            raise RuntimeError("R1 algebraic source-bound slot order drifted")
        if not set(out).issubset(set(R1_ORDINARY_TERMS)):
            raise RuntimeError("R1 source-bound subset contains an unregistered slot")
        return out

    @property
    def truth_boundary(self) -> dict[str, bool]:
        return copy.deepcopy(_TRUTH_BOUNDARY)

    def report(self) -> dict[str, Any]:
        fixed = self.fixed_multiplier_balls()
        phi = self.phi_ball()
        u = self.u_ball()
        Au = self.radial_average_u_ball()
        numerator = self.r1_algebraic_numerator_slots()
        after_j2 = self.r1_algebraic_after_j2()
        subset = self.combined_r1_source_bound_subset()
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
                "R2_agent1_self_certificate_receipt_sha256": self.r2_source.sha256,
                "R2_algebraic_receipt_sha256": self.algebraic.sha256,
                "post_J_bridge_receipt_sha256": self.bridge.sha256,
                "source_u_ball_receipt_sha256": self.bridge.u_source.sha256,
                "source_Phi_ball_receipt_sha256": self.bridge.u_source.mixed.phi_source.sha256,
                "source_axis_receipt_sha256": self.domain.sha256,
            },
            "selected_source_compatible_inputs": {
                "h": float(self.domain.h),
                "D": float(self.domain.D),
                "j0": float(self.domain.j0),
                "rho": float(self.domain.coefficient_rho),
                "Lambda_inverse_upper": 1.0,
                "choice_status": (
                    "repository-autonomous source-compatible existence choice; "
                    "not recovered hidden source data"
                ),
            },
            "fixed_multiplier_bounds": {
                name: asdict(bound) for name, bound in fixed.items()
            },
            "source_compatible_radius_one_Phi_ball": asdict(phi),
            "source_compatible_radius_one_u_ball": asdict(u),
            "source_compatible_radial_average_u_ball": asdict(Au),
            "R1_algebraic_ordinary_numerator_slots": {
                name: asdict(bound) for name, bound in numerator.items()
            },
            "R1_algebraic_ordinary_after_J2": {
                name: asdict(bound) for name, bound in after_j2.items()
            },
            "ledger_progress": {
                "R1_algebraic_slots_source_bound": list(R1_ALGEBRAIC_ORDINARY_TERMS),
                "R1_source_bound_ordinary_slots_count": len(subset),
                "R1_total_ordinary_slots": len(R1_ORDINARY_TERMS),
                "R1_ordinary_slots_still_required": len(R1_REMAINING_ORDINARY_TERMS),
                "R1_remaining_ordinary_slots": list(R1_REMAINING_ORDINARY_TERMS),
                "R1_mixed_slot_handled_upstream_by_A1_668": True,
                "R2_agent1_self_certificate_materialized_upstream": True,
                "full_R1_still_open": True,
                "M_K_not_promoted": True,
            },
            "independent_audit_boundary": {
                "A4_source_axis_independent_PASS_exists": True,
                "A4_source_Phi_ball_conditional_PASS_exists": True,
                "A4_source_u_ball_independent_PASS_exists": True,
                "R2_full_independent_A4_admission": False,
                "this_R1_subset_is_agent1_self_certificate": True,
                "this_R1_subset_requires_independent_A4_audit": True,
            },
            "truth_boundary": self.truth_boundary,
            "heldout_ns_residual_assessed": False,
            "pde_validated": False,
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
    payload = KokunoPA10SourceR1AlgebraicOrdinarySlots().save_report(args.output)
    print("receipt_sha256=", payload["receipt_sha256"])
    print("R1_after_J2=", payload["R1_algebraic_ordinary_after_J2"])
    print("progress=", payload["ledger_progress"])


if __name__ == "__main__":
    _main()
