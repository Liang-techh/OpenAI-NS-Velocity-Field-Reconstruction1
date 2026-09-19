"""Source-compatible ``d u zeta_* Phi`` PA.10 R1 slot after ``J_2``.

Pinned public provenance is KokunoYumeto/yang-mills-interacting-workbench at
``143f6773feb424ad9ed3a8d116653200f20346b7``, file
``navier-stokes/navier_stokes_workbench.tex`` (corrected reader dated
2026-09-09, Zenodo 22678406).

The corrected reconstruction defines

    zeta_* = -L H_* / (H_*^2 + sigma_*^2)

and uses the exact identity ``H_* zeta_*/L=-chi``.  In PA.10 the R1 numerator
contains the ordinary monomial ``d u zeta_* Phi``.  Agent-1 #702 already
source-binds five derivative-free R1 ordinary monomials.  This module binds
exactly this sixth ordinary monomial and applies the existing ``L^-1`` +
``J_2`` bridge once.

The coefficient bound for the fixed analytic multiplier ``zeta_*`` is obtained
on the already-certified source-compatible Cauchy neighborhood.  It uses the
same factorization

    H_*^2 + sigma_*^2 = (H_*-i sigma_*)(H_*+i sigma_*)

as the source-domain certificate.  Around the real enlarged interval each
factor has modulus at least ``sigma_*-r sup|H_*'|`` on the Cauchy radius ``r``.
A direct polynomial majorant supplies ``sup|L H_*|``.  The resulting complex
supremum is converted to the source alpha=0 coefficient norm with

    sum_{beta>=0} (beta+1)^2 x^beta = (1+x)/(1-x)^3,
    x=rho/r.

All newly introduced bound arithmetic uses exact ``Fraction`` values of the
configured binary64 inputs and outward rounding.  The selected sigma/rho datum
remains a repository-autonomous source-compatible existence choice, not a
recovery of hidden Kokuno/OpenAI parameters.  This is an Agent-1 source-bound
increment, not independent validation and not an NS residual.
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

import numpy as np

from .kokuno_pa10_post_j_remainder_bridge import R1_ORDINARY_TERMS
from .kokuno_pa10_remainder_ball_bounds import BallFactorBound
from .kokuno_pa10_source_r1_algebraic_ordinary_slots import (
    R1_ALGEBRAIC_ORDINARY_TERMS,
    R1_REMAINING_ORDINARY_TERMS,
    KokunoPA10SourceR1AlgebraicOrdinarySlots,
)

SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
CORRECTED_RELEASE = "zenodo:22678406"
CORRECTED_RELEASE_DATE = "2026-09-09"
SCHEMA = "kokuno-pa10-source-r1-zeta-ordinary-slot-v1"

R1_ZETA_ORDINARY_TERMS = ("d_u_zeta_star_times_Phi",)
R1_SOURCE_BOUND_THROUGH_ZETA_TERMS = (
    *R1_ALGEBRAIC_ORDINARY_TERMS,
    *R1_ZETA_ORDINARY_TERMS,
)
R1_REMAINING_AFTER_ZETA_TERMS = tuple(
    name for name in R1_REMAINING_ORDINARY_TERMS if name not in R1_ZETA_ORDINARY_TERMS
)

_SOURCE_FORMULAS = {
    "R1": "L R1=[W+h(1-2eta U)+d u zeta_*]Phi+W YPhi_Y+H_c Phi_eta",
    "zeta_star": "zeta_*=-L H_*/(H_*^2+sigma_*^2)",
    "zeta_identity": "H_* zeta_*/L=-chi",
    "denominator_factorization": (
        "H_*^2+sigma_*^2=(H_*-i sigma_*)(H_*+i sigma_*)"
    ),
    "ordinary_slot": "d u zeta_* Phi",
    "coefficient_weight_alpha0": (
        "a_{0,beta}=rho^{-beta} beta!/((beta+1)^2)"
    ),
    "cauchy_sum": "sum_(beta>=0) (beta+1)^2 x^beta=(1+x)/(1-x)^3",
    "post_J2": "J_2[L^-1(d u zeta_* Phi)]",
}

_TRUTH_BOUNDARY = {
    "public_reconstruction_source": True,
    "source_compatible_zeta_star_value_api_executable": True,
    "source_compatible_zeta_star_coefficient_norm_machine_bound": True,
    "source_compatible_R1_zeta_ordinary_numerator_slot_machine_bound": True,
    "source_compatible_R1_zeta_ordinary_post_J2_slot_machine_bound": True,
    "source_compatible_R1_ordinary_source_bound_count_is_six": True,
    "source_compatible_all_R1_ordinary_post_J2_slots_machine_bound": False,
    "standalone_Y_Phi_Y_coefficient_norm_invented": False,
    "standalone_Phi_eta_coefficient_norm_invented": False,
    "this_R1_increment_independent_agent4_audit_required": True,
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


def _q(value: float, name: str) -> Fraction:
    out = float(value)
    if not math.isfinite(out):
        raise ValueError(f"{name} must be finite")
    return Fraction.from_float(out)


def _q_nonnegative(value: float, name: str) -> Fraction:
    out = _q(value, name)
    if out < 0:
        raise ValueError(f"{name} must be nonnegative")
    return out


def _upper(value: Fraction) -> float:
    if value < 0:
        raise ValueError("upper-bound helper requires a nonnegative value")
    out = float(value)
    if not math.isfinite(out):
        raise OverflowError("upper bound is outside binary64 range")
    if Fraction.from_float(out) < value:
        out = math.nextafter(out, math.inf)
    return float(out)


def _lower(value: Fraction) -> float:
    if value < 0:
        raise ValueError("lower-bound helper requires a nonnegative value")
    out = float(value)
    if not math.isfinite(out):
        raise OverflowError("lower bound is outside binary64 range")
    if Fraction.from_float(out) > value:
        out = math.nextafter(out, -math.inf)
    return float(out)


@dataclass(frozen=True)
class KokunoPA10SourceR1ZetaOrdinarySlot:
    """Bind ``d u zeta_* Phi`` and merge it with the five #702 R1 slots."""

    algebraic_r1: KokunoPA10SourceR1AlgebraicOrdinarySlots = field(
        default_factory=KokunoPA10SourceR1AlgebraicOrdinarySlots
    )

    def __post_init__(self) -> None:
        if not isinstance(self.algebraic_r1, KokunoPA10SourceR1AlgebraicOrdinarySlots):
            raise TypeError("algebraic_r1 must be KokunoPA10SourceR1AlgebraicOrdinarySlots")
        if not 0.0 < self.domain.coefficient_rho < self.domain.cauchy_radius:
            raise ValueError("coefficient rho must lie inside the Cauchy radius")
        cert = self.zeta_star_certificate()
        if cert["H_factor_abs_lower_on_cauchy"] <= 0.0:
            raise ValueError("Cauchy neighborhood does not keep H_* +/- i sigma away from zero")

    @property
    def bridge(self):
        return self.algebraic_r1.bridge

    @property
    def domain(self):
        return self.algebraic_r1.domain

    @property
    def product_calculator(self):
        return self.algebraic_r1.product_calculator

    def zeta_star_values(self, eta: Any) -> np.ndarray:
        """Vectorized public-source ``zeta_*`` on the certified real interval."""
        return np.asarray(self.domain.axis_state(eta)["zeta_star"], dtype=float)

    def _zeta_star_exact_certificate(self) -> dict[str, Fraction]:
        rho = _q_nonnegative(self.domain.coefficient_rho, "coefficient_rho")
        radius = _q_nonnegative(self.domain.cauchy_radius, "cauchy_radius")
        margin = _q_nonnegative(self.domain.enlarged_real_margin, "real margin")
        sigma = _q_nonnegative(self.domain.sigma_star, "sigma_star")
        h = _q_nonnegative(self.domain.h, "h")
        D = Fraction(1, 2) - h
        j0 = _q(self.domain.j0, "j0")
        if not 0 < rho < radius or sigma <= 0 or D <= 0:
            raise ValueError("invalid source-compatible Cauchy parameters")

        z_abs = Fraction(1) + margin + radius
        j0_abs = abs(j0)

        # H_*=j0+(D+4)z-j0 z^2-4z^3 and
        # H_*'=(D+4)-2j0 z-12z^2.
        H_upper = (
            j0_abs
            + (D + 4) * z_abs
            + j0_abs * z_abs * z_abs
            + 4 * z_abs**3
        )
        H_prime_upper = (D + 4) + 2 * j0_abs * z_abs + 12 * z_abs**2
        H_variation = radius * H_prime_upper
        H_factor_lower = sigma - H_variation
        if H_factor_lower <= 0:
            raise ValueError("Cauchy H variation reaches sigma_star")
        denominator_lower = H_factor_lower**2

        # L=1-2h z^2, hence |L|<=1+2h|z|^2.
        L_upper = Fraction(1) + 2 * h * z_abs**2
        zeta_complex_sup = L_upper * H_upper / denominator_lower

        x = rho / radius
        cauchy_weight_sum = (Fraction(1) + x) / (Fraction(1) - x) ** 3
        coefficient_norm = zeta_complex_sup * cauchy_weight_sum
        return {
            "rho_over_cauchy_radius": x,
            "cauchy_z_abs_upper": z_abs,
            "cauchy_H_abs_upper": H_upper,
            "cauchy_H_prime_abs_upper": H_prime_upper,
            "cauchy_H_variation_abs_upper": H_variation,
            "H_factor_abs_lower_on_cauchy": H_factor_lower,
            "H_square_plus_sigma_square_abs_lower_on_cauchy": denominator_lower,
            "cauchy_L_abs_upper": L_upper,
            "zeta_star_complex_sup_upper": zeta_complex_sup,
            "cauchy_weight_sum_upper": cauchy_weight_sum,
            "zeta_star_coefficient_norm_upper": coefficient_norm,
        }

    def zeta_star_certificate(self) -> dict[str, float]:
        exact = self._zeta_star_exact_certificate()
        out = {name: _upper(value) for name, value in exact.items()}
        # These two are genuine lower bounds, not upper bounds.
        out["H_factor_abs_lower_on_cauchy"] = _lower(
            exact["H_factor_abs_lower_on_cauchy"]
        )
        out["H_square_plus_sigma_square_abs_lower_on_cauchy"] = _lower(
            exact["H_square_plus_sigma_square_abs_lower_on_cauchy"]
        )
        return out

    def zeta_star_coefficient_ball(self) -> BallFactorBound:
        cert = self._zeta_star_exact_certificate()
        return BallFactorBound(_upper(cert["zeta_star_coefficient_norm_upper"]), 0.0)

    def r1_zeta_numerator_slot(self) -> dict[str, BallFactorBound]:
        fixed = self.algebraic_r1.algebraic.fixed_multiplier_balls()
        d = fixed["d"]
        u = self.algebraic_r1.u_ball()
        phi = self.algebraic_r1.phi_ball()
        zeta = self.zeta_star_coefficient_ball()
        return {
            "d_u_zeta_star_times_Phi": self.product_calculator.product(d, u, zeta, phi)
        }

    def r1_zeta_after_j2(self) -> dict[str, BallFactorBound]:
        slots = self.r1_zeta_numerator_slot()
        return {
            name: self.bridge.ordinary_term_after_jnu(bound, nu=2)
            for name, bound in slots.items()
        }

    def combined_r1_source_bound_subset(self) -> dict[str, BallFactorBound]:
        base = self.algebraic_r1.combined_r1_source_bound_subset()
        zeta = self.r1_zeta_numerator_slot()
        overlap = set(base).intersection(zeta)
        if overlap:
            raise RuntimeError(f"R1 source-bound subset overlap: {sorted(overlap)}")
        out = {**base, **zeta}
        if tuple(out) != R1_SOURCE_BOUND_THROUGH_ZETA_TERMS:
            raise RuntimeError("R1 source-bound slot order drifted")
        if not set(out).issubset(set(R1_ORDINARY_TERMS)):
            raise RuntimeError("R1 source-bound subset contains an unregistered slot")
        return out

    def combined_r1_after_j2_subset(self) -> dict[str, BallFactorBound]:
        base = self.algebraic_r1.r1_algebraic_after_j2()
        zeta = self.r1_zeta_after_j2()
        out = {**base, **zeta}
        if tuple(out) != R1_SOURCE_BOUND_THROUGH_ZETA_TERMS:
            raise RuntimeError("R1 post-J2 source-bound slot order drifted")
        return out

    @property
    def truth_boundary(self) -> dict[str, bool]:
        return copy.deepcopy(_TRUTH_BOUNDARY)

    def report(self) -> dict[str, Any]:
        zeta_cert = self.zeta_star_certificate()
        zeta_ball = self.zeta_star_coefficient_ball()
        numerator = self.r1_zeta_numerator_slot()
        after_j2 = self.r1_zeta_after_j2()
        combined = self.combined_r1_source_bound_subset()
        combined_after = self.combined_r1_after_j2_subset()
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
                "R1_algebraic_five_slot_receipt_sha256": self.algebraic_r1.sha256,
                "post_J_bridge_receipt_sha256": self.bridge.sha256,
                "source_axis_receipt_sha256": self.domain.sha256,
            },
            "selected_source_compatible_inputs": {
                "sigma_star": float(self.domain.sigma_star),
                "rho": float(self.domain.coefficient_rho),
                "cauchy_radius": float(self.domain.cauchy_radius),
                "enlarged_real_margin": float(self.domain.enlarged_real_margin),
                "choice_status": (
                    "repository-autonomous source-compatible existence choice; "
                    "not recovered hidden source data"
                ),
            },
            "zeta_star_certificate": zeta_cert,
            "zeta_star_coefficient_ball": asdict(zeta_ball),
            "R1_zeta_ordinary_numerator_slot": {
                name: asdict(bound) for name, bound in numerator.items()
            },
            "R1_zeta_ordinary_after_J2": {
                name: asdict(bound) for name, bound in after_j2.items()
            },
            "R1_source_bound_after_J2_subset": {
                name: asdict(bound) for name, bound in combined_after.items()
            },
            "ledger_progress": {
                "R1_source_bound_ordinary_slots": list(combined),
                "R1_source_bound_ordinary_slots_count": len(combined),
                "R1_total_ordinary_slots": len(R1_ORDINARY_TERMS),
                "R1_ordinary_slots_still_required": len(R1_REMAINING_AFTER_ZETA_TERMS),
                "R1_remaining_ordinary_slots": list(R1_REMAINING_AFTER_ZETA_TERMS),
                "R1_mixed_slot_handled_upstream_by_A1_668": True,
                "R2_agent1_self_certificate_materialized_upstream": True,
                "full_R1_still_open": True,
                "M_K_not_promoted": True,
            },
            "independent_audit_boundary": {
                "A4_source_axis_independent_PASS_exists": True,
                "A4_source_Phi_ball_conditional_PASS_exists": True,
                "A4_source_u_ball_independent_PASS_exists": True,
                "this_zeta_R1_slot_is_agent1_self_certificate": True,
                "this_zeta_R1_slot_requires_independent_A4_audit": True,
                "R2_full_independent_A4_admission": False,
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
    payload = KokunoPA10SourceR1ZetaOrdinarySlot().save_report(args.output)
    print("receipt_sha256=", payload["receipt_sha256"])
    print("zeta_star_certificate=", payload["zeta_star_certificate"])
    print("R1_zeta_after_J2=", payload["R1_zeta_ordinary_after_J2"])
    print("progress=", payload["ledger_progress"])


if __name__ == "__main__":
    _main()
