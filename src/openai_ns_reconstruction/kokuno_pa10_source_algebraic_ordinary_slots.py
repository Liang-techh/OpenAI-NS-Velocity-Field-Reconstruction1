"""Source-compatible PA.10 algebraic ordinary R2 slots after ``J_1``.

Pinned public provenance is KokunoYumeto/yang-mills-interacting-workbench at
``143f6773feb424ad9ed3a8d116653200f20346b7``, file
``navier-stokes/navier_stokes_workbench.tex`` (corrected reader dated
2026-09-09, Zenodo 22678406).

After Agent-1 #676 bound the three pressure monomials in the PA.10 R2 ledger,
the next source-displayed ordinary terms that require no logarithmic-radial or
unintegrated eta derivative are

    A u,
    4 A eta U_* u,
    d U_*' u,
    2 A Lambda^-1 eta u^2,

where ``A=1/2+h``, ``U_*=4 eta+j_0``, ``U_*'=4``, ``d=1-eta^2`` and the source
regime gives ``Lambda^-1<=1``.  This module machine-bounds exactly those four
terms from the already materialized source-compatible u-ball and fixed axis
multipliers, then applies the #668 ``L^-1`` + ``J_1`` bridge exactly once.

Together with #676 this fills 7/11 ordinary R2 slots.  The remaining four R2
ordinary slots carry ``Y u_Y`` or ``u_eta`` and still require their own
post-radial-inverse degree/derivative cancellation.  Full R2, R1, M/K, the
fixed point, global matched pressure/velocity and PDE validation remain false.

The selected sigma/rho/pressure datum is a repository-autonomous
source-compatible existence choice, not recovered hidden Kokuno/OpenAI data.
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

from .kokuno_pa10_post_j_remainder_bridge import (
    R2_ORDINARY_TERMS,
    KokunoPA10PostJRemainderBridge,
)
from .kokuno_pa10_remainder_ball_bounds import BallFactorBound
from .kokuno_pa10_source_pressure_ordinary_slots import (
    R2_PRESSURE_ORDINARY_TERMS,
    KokunoPA10SourcePressureOrdinarySlots,
)

SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
CORRECTED_RELEASE = "zenodo:22678406"
CORRECTED_RELEASE_DATE = "2026-09-09"
SCHEMA = "kokuno-pa10-source-algebraic-ordinary-slots-v1"

R2_ALGEBRAIC_ORDINARY_TERMS = (
    "A_times_u",
    "4A_eta_Ustar_times_u",
    "d_Ustar_eta_times_u",
    "2A_lambda_inv_eta_u_squared",
)

R2_REMAINING_DERIVATIVE_ORDINARY_TERMS = (
    "W_star_times_Y_u_Y",
    "lambda_inv_2D_eta_Au_times_Y_u_Y",
    "H_star_times_u_eta",
    "lambda_inv_d_u_u_eta",
)

_SOURCE_FORMULAS = {
    "R2_algebraic_subset": (
        "[A(1-4eta U_*)+d U_*']u-2Aeta Lambda^-1u^2"
    ),
    "A": "A=1/2+h",
    "U_star": "U_*=4eta+j_0",
    "U_star_eta": "U_*'=4",
    "d": "d=1-eta^2",
    "Lambda": "Lambda>=1, hence Lambda^-1<=1",
    "post_J1": "J_1[L^-1 R2_numerator]",
    "coefficient_weight_alpha0": (
        "a_{0,beta}=rho^{-beta} beta!/((beta+1)^2)"
    ),
}

_TRUTH_BOUNDARY = {
    "public_reconstruction_source": True,
    "source_compatible_R2_algebraic_ordinary_numerator_slots_machine_bound": True,
    "source_compatible_R2_algebraic_ordinary_post_J1_slots_machine_bound": True,
    "source_compatible_Ustar_coefficient_norm_machine_bound": True,
    "source_compatible_u_radius_one_ball_consumed": True,
    "uniform_Lambda_inverse_upper_for_Lambda_ge_1_used": True,
    "source_axis_domain_independent_agent4_pass_exists": True,
    "source_pressure_ordinary_slots_independent_agent4_pass_exists": True,
    "this_algebraic_slot_increment_independent_agent4_audit_required": True,
    "all_R2_ordinary_slots_source_bound": False,
    "all_R1_ordinary_slots_source_bound": False,
    "source_full_post_J1_R2_radius_one_ball_bound_machine_bound": False,
    "source_full_post_J2_R1_radius_one_ball_bound_machine_bound": False,
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
class KokunoPA10SourceAlgebraicOrdinarySlots:
    """Bind the four non-derivative/non-logradial PA.10 R2 ordinary slots."""

    pressure_slots: KokunoPA10SourcePressureOrdinarySlots = field(
        default_factory=KokunoPA10SourcePressureOrdinarySlots
    )

    def __post_init__(self) -> None:
        if not isinstance(self.pressure_slots, KokunoPA10SourcePressureOrdinarySlots):
            raise TypeError("pressure_slots must be KokunoPA10SourcePressureOrdinarySlots")
        if not 0.0 < self.domain.coefficient_rho < self.domain.cauchy_radius:
            raise ValueError("coefficient rho must lie inside the Cauchy radius")

    @property
    def bridge(self) -> KokunoPA10PostJRemainderBridge:
        return self.pressure_slots.bridge

    @property
    def domain(self):
        return self.bridge.domain

    @property
    def product_calculator(self):
        return self.pressure_slots.pressure_calculator

    def axis_multiplier_values(self, eta: Any) -> dict[str, np.ndarray]:
        """Vectorized executable values of the fixed PA.10 axis multipliers."""
        eta_arr = np.asarray(eta, dtype=float)
        if np.any(~np.isfinite(eta_arr)):
            raise ValueError("eta must be finite")
        E = 1.0 + float(self.domain.enlarged_real_margin)
        if np.any(np.abs(eta_arr) > E):
            raise ValueError("eta lies outside the certified enlarged real interval")
        h = float(self.domain.h)
        j0 = float(self.domain.j0)
        return {
            "eta": eta_arr,
            "A": np.full_like(eta_arr, 0.5 + h, dtype=float),
            "U_star": 4.0 * eta_arr + j0,
            "U_star_eta": np.full_like(eta_arr, 4.0, dtype=float),
            "d": 1.0 - eta_arr * eta_arr,
        }

    def ustar_coefficient_ball(self) -> BallFactorBound:
        """Exact finite-derivative coefficient bound for ``U_*=4eta+j0``.

        For alpha=0, ``a_{0,0}=1`` and ``a_{0,1}=1/(4rho)``.  ``U_*`` has
        only beta=0,1 derivatives, hence

            ||U_*||_rho <= max(4 E + j0, 16 rho),

        where ``E=1+margin`` on the certified enlarged real interval.
        """
        rho = Fraction.from_float(float(self.domain.coefficient_rho))
        margin = Fraction.from_float(float(self.domain.enlarged_real_margin))
        j0 = Fraction.from_float(float(self.domain.j0))
        if rho <= 0 or margin < 0 or j0 < 0:
            raise ValueError("source-compatible rho/margin/j0 must be valid")
        E = Fraction(1) + margin
        return BallFactorBound(_upper(max(4 * E + j0, 16 * rho)), 0.0)

    def fixed_multiplier_balls(self) -> dict[str, BallFactorBound]:
        h = Fraction.from_float(float(self.domain.h))
        if h < 0:
            raise ValueError("h must be nonnegative")
        A = Fraction(1, 2) + h
        return {
            "A": BallFactorBound(_upper(A), 0.0),
            "eta": self.pressure_slots.eta_coefficient_ball(),
            "U_star": self.ustar_coefficient_ball(),
            "U_star_eta": BallFactorBound(4.0, 0.0),
            "d": self.bridge.u_source.mixed.d_coefficient_ball(),
        }

    def r2_algebraic_numerator_slots(self) -> dict[str, BallFactorBound]:
        """Return source-compatible envelopes for exactly four R2 monomials."""
        fixed = self.fixed_multiplier_balls()
        u = self.bridge.u_source.u_ball()
        calc = self.product_calculator
        h = Fraction.from_float(float(self.domain.h))
        A = Fraction(1, 2) + h

        eta_U_u = calc.product(fixed["eta"], fixed["U_star"], u)
        d_u = calc.product(fixed["d"], u)
        eta_u2 = calc.product(fixed["eta"], u, u)

        return {
            "A_times_u": _scale_exact(u, A),
            "4A_eta_Ustar_times_u": _scale_exact(eta_U_u, 4 * A),
            "d_Ustar_eta_times_u": _scale_exact(d_u, Fraction(4)),
            # The source regime Lambda>=1 is used only as Lambda^-1<=1.
            "2A_lambda_inv_eta_u_squared": _scale_exact(eta_u2, 2 * A),
        }

    def r2_algebraic_after_j1(self) -> dict[str, BallFactorBound]:
        slots = self.r2_algebraic_numerator_slots()
        return {
            name: self.bridge.ordinary_term_after_jnu(bound, nu=1)
            for name, bound in slots.items()
        }

    def combined_source_bound_r2_subset(self) -> dict[str, BallFactorBound]:
        """Merge #676 pressure slots and this increment without placeholders."""
        pressure = self.pressure_slots.r2_pressure_numerator_slots()
        algebraic = self.r2_algebraic_numerator_slots()
        overlap = set(pressure).intersection(algebraic)
        if overlap:
            raise RuntimeError(f"R2 source-bound subset overlap: {sorted(overlap)}")
        out = {**algebraic, **pressure}
        if not set(out).issubset(set(R2_ORDINARY_TERMS)):
            raise RuntimeError("source-bound subset contains an unregistered R2 slot")
        return out

    @property
    def truth_boundary(self) -> dict[str, bool]:
        return copy.deepcopy(_TRUTH_BOUNDARY)

    def report(self) -> dict[str, Any]:
        fixed = self.fixed_multiplier_balls()
        u = self.bridge.u_source.u_ball()
        numerator = self.r2_algebraic_numerator_slots()
        after_j1 = self.r2_algebraic_after_j1()
        combined = self.combined_source_bound_r2_subset()
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
                "pressure_slots_receipt_sha256": self.pressure_slots.sha256,
                "post_J_bridge_receipt_sha256": self.bridge.sha256,
                "source_u_ball_receipt_sha256": self.bridge.u_source.sha256,
                "source_axis_receipt_sha256": self.domain.sha256,
            },
            "selected_source_compatible_inputs": {
                "h": float(self.domain.h),
                "j0": float(self.domain.j0),
                "rho": float(self.domain.coefficient_rho),
                "enlarged_real_margin": float(self.domain.enlarged_real_margin),
                "Lambda_inverse_upper": 1.0,
                "choice_status": (
                    "repository-autonomous source-compatible existence choice; "
                    "not recovered hidden source data"
                ),
            },
            "fixed_multiplier_bounds": {
                name: asdict(bound) for name, bound in fixed.items()
            },
            "source_compatible_radius_one_u_ball": asdict(u),
            "R2_algebraic_ordinary_numerator_slots": {
                name: asdict(bound) for name, bound in numerator.items()
            },
            "R2_algebraic_ordinary_after_J1": {
                name: asdict(bound) for name, bound in after_j1.items()
            },
            "ledger_progress": {
                "R2_pressure_slots_source_bound": list(R2_PRESSURE_ORDINARY_TERMS),
                "R2_algebraic_slots_source_bound": list(R2_ALGEBRAIC_ORDINARY_TERMS),
                "R2_source_bound_ordinary_slots_count": len(combined),
                "R2_total_ordinary_slots": len(R2_ORDINARY_TERMS),
                "R2_ordinary_slots_still_required": len(R2_ORDINARY_TERMS) - len(combined),
                "R2_remaining_derivative_slots": list(
                    R2_REMAINING_DERIVATIVE_ORDINARY_TERMS
                ),
                "R1_ordinary_slots_still_required": 11,
                "mixed_slots_handled_upstream_by_A1_668": True,
                "full_R1_R2_not_promoted": True,
            },
            "independent_audit_boundary": {
                "A4_source_axis_repair_PR": 671,
                "A4_source_axis_dedicated_run": 35456226819,
                "A4_source_axis_independent_PASS_exists": True,
                "A4_pressure_slots_audit_PR": 678,
                "A4_pressure_slots_dedicated_run": 35459240658,
                "A4_pressure_slots_independent_PASS_exists": True,
                "this_new_algebraic_subset_is_agent1_self_certificate": True,
                "this_new_algebraic_subset_requires_independent_A4_audit": True,
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
    print(
        json.dumps(
            KokunoPA10SourceAlgebraicOrdinarySlots().save_report(args.output),
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    _main()
