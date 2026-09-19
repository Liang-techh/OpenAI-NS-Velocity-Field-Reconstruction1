"""Source-compatible PA.10 derivative-bearing ordinary R2 slots after ``J_1``.

Pinned public provenance is KokunoYumeto/yang-mills-interacting-workbench at
``143f6773feb424ad9ed3a8d116653200f20346b7``, file
``navier-stokes/navier_stokes_workbench.tex`` (corrected reader dated
2026-09-09, Zenodo 22678406).

Agent-1 #684 source-binds seven of the eleven ordinary PA.10 R2 numerator
slots.  The remaining four are

    W_* Y u_Y,
    Lambda^-1 2D eta A(u) Y u_Y,
    H_* u_eta,
    Lambda^-1 d u u_eta.

The source coefficient-space proof states that every non-mixed monomial is
controlled by multiplication plus a radial inverse, with at most one of the
two derivative types.  In the displayed weight-ratio proof, removing the
parameter derivative uses the first weight ratio, removing the logarithmic
radial derivative replaces its degree factor by one, and averaging contributes
an ``(i+1)^-1`` divisor that cancels the added-degree factor.

This module makes exactly those four post-``J_1`` operators executable without
inventing standalone ``Y u_Y`` or ``u_eta`` coefficient norms.  It also binds
fixed polynomial coefficient balls for ``W_*`` and ``H_*`` from

    U_*=4 eta+j_0,
    H_*=D eta+d U_*,
    W_*=1-d U_*'-2D eta U_*.

The selected sigma/rho datum remains a repository-autonomous source-compatible
existence choice.  This closes an Agent-1 self-certificate for the complete
post-``J_1`` R2 radius-one-ball envelope; independent Agent-4 admission is
still required, and all R1 / M,K / global velocity / PDE gates remain closed.
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

from .kokuno_pa10_post_j_remainder_bridge import R2_ORDINARY_TERMS
from .kokuno_pa10_remainder_ball_bounds import BallFactorBound
from .kokuno_pa10_source_algebraic_ordinary_slots import (
    R2_ALGEBRAIC_ORDINARY_TERMS,
    R2_REMAINING_DERIVATIVE_ORDINARY_TERMS,
    KokunoPA10SourceAlgebraicOrdinarySlots,
)
from .kokuno_pa10_source_pressure_ordinary_slots import R2_PRESSURE_ORDINARY_TERMS

SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
CORRECTED_RELEASE = "zenodo:22678406"
CORRECTED_RELEASE_DATE = "2026-09-09"
SCHEMA = "kokuno-pa10-source-derivative-ordinary-slots-v1"

R2_DERIVATIVE_ORDINARY_TERMS = R2_REMAINING_DERIVATIVE_ORDINARY_TERMS

_SOURCE_FORMULAS = {
    "axis": (
        "U_*=4eta+j_0; H_*=Deta+dU_*; "
        "W_*=1-dU_*'-2DetaU_*; d=1-eta^2"
    ),
    "R2_derivative_subset": (
        "W_*Yu_Y + Lambda^-1[-2Deta A(u)]Yu_Y + "
        "H_*u_eta + Lambda^-1 d u u_eta"
    ),
    "coefficient_weight": (
        "a_{alpha,beta}=20^{-alpha}rho^{-beta}beta!binom(alpha+beta,beta)/"
        "((alpha+1)^2(beta+1)^2)"
    ),
    "product": "||FG||_rho <= C_sq^2 ||F||_rho ||G||_rho",
    "radial_inverse": (
        "(J_nu F)_{alpha+1}=F_alpha/((alpha+1)(alpha+nu))"
    ),
    "single_derivative_rule": (
        "source weight-ratio proof: if logarithmic derivative is absent replace "
        "its degree factor by one; if parameter derivative is absent use the "
        "first weight ratio; averaged F contributes (i+1)^-1 cancellation"
    ),
    "average": "A F=Y^-1 I F; ||A F||_rho<=||F||_rho",
    "u_u_eta_identity": "u u_eta=(1/2) partial_eta(u^2)",
    "Au_YuY_identity": (
        "A(u) Yu_Y = Y partial_Y(A(u)u) - u^2 + A(u)u"
    ),
    "Lambda": "Lambda>=1, hence Lambda^-1<=1",
}

_TRUTH_BOUNDARY = {
    "public_reconstruction_source": True,
    "source_compatible_Wstar_coefficient_norm_machine_bound": True,
    "source_compatible_Hstar_coefficient_norm_machine_bound": True,
    "source_single_derivative_post_J1_rules_executable": True,
    "source_averaging_contraction_executable": True,
    "source_compatible_R2_derivative_ordinary_post_J1_slots_machine_bound": True,
    "source_compatible_all_R2_ordinary_post_J1_slots_machine_bound": True,
    "source_compatible_full_post_J1_R2_radius_one_ball_bound_agent1_self_certificate": True,
    "standalone_Y_u_Y_coefficient_norm_invented": False,
    "standalone_u_eta_coefficient_norm_invented": False,
    "uniform_Lambda_inverse_upper_for_Lambda_ge_1_used": True,
    "source_axis_domain_independent_agent4_pass_exists": True,
    "source_pressure_ordinary_slots_independent_agent4_pass_exists": True,
    "source_u_ball_independent_agent4_pass_exists": True,
    "algebraic_slot_increment_independent_agent4_audit_required": True,
    "this_derivative_slot_increment_independent_agent4_audit_required": True,
    "full_R2_independent_agent4_admission": False,
    "all_R1_ordinary_slots_source_bound": False,
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


def _product_constant_exact() -> Fraction:
    # Elementary pi < 22/7.  C_sq=4*pi^2/3 and product constant=C_sq^2.
    pi_upper = Fraction(22, 7)
    c_sq_upper = Fraction(4, 3) * pi_upper * pi_upper
    return c_sq_upper * c_sq_upper


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


def _product_exact(*bounds: BallFactorBound) -> BallFactorBound:
    if not bounds:
        return BallFactorBound(1.0, 0.0)
    if not all(isinstance(item, BallFactorBound) for item in bounds):
        raise TypeError("all factors must be BallFactorBound values")
    constant = _product_constant_exact() ** max(len(bounds) - 1, 0)
    norm = constant
    for item in bounds:
        norm *= _q_nonnegative(item.norm, "factor norm")
    lip = Fraction(0)
    for i, item in enumerate(bounds):
        if item.lipschitz == 0.0:
            continue
        term = constant * _q_nonnegative(item.lipschitz, "factor lipschitz")
        for j, other in enumerate(bounds):
            if i != j:
                term *= _q_nonnegative(other.norm, "factor norm")
        lip += term
    return BallFactorBound(_upper(norm), _upper(lip))


@dataclass(frozen=True)
class KokunoPA10SourceDerivativeOrdinarySlots:
    """Bind the four derivative-bearing ordinary PA.10 R2 slots after J1."""

    algebraic: KokunoPA10SourceAlgebraicOrdinarySlots = field(
        default_factory=KokunoPA10SourceAlgebraicOrdinarySlots
    )

    def __post_init__(self) -> None:
        if not isinstance(self.algebraic, KokunoPA10SourceAlgebraicOrdinarySlots):
            raise TypeError("algebraic must be KokunoPA10SourceAlgebraicOrdinarySlots")
        if not 0.0 < self.rho < self.domain.cauchy_radius:
            raise ValueError("coefficient rho must lie inside the Cauchy radius")

    @property
    def bridge(self):
        return self.algebraic.bridge

    @property
    def domain(self):
        return self.algebraic.domain

    @property
    def rho(self) -> float:
        return float(self.domain.coefficient_rho)

    @property
    def product_constant_upper(self) -> float:
        return _upper(_product_constant_exact())

    def radial_average_u_ball(self) -> BallFactorBound:
        """Averaging divides each radial coefficient by alpha+1, so is contractive."""
        u = self.bridge.u_source.u_ball()
        return BallFactorBound(u.norm, u.lipschitz)

    def wstar_coefficient_ball(self) -> BallFactorBound:
        """Finite-derivative coefficient bound for W_*=1-dU_*'-2DetaU_*.

        With U_*=4eta+j0 and D=1/2-h,
        W_*=-3-2D*j0*eta+8h*eta^2.
        """
        rho = Fraction.from_float(self.rho)
        margin = Fraction.from_float(float(self.domain.enlarged_real_margin))
        h = Fraction.from_float(float(self.domain.h))
        D = Fraction(1, 2) - h
        j0 = Fraction.from_float(float(self.domain.j0))
        if min(rho, h, D, j0) < 0:
            raise ValueError("source-compatible parameters must be nonnegative")
        E = Fraction(1) + margin
        beta0 = 3 + 2 * D * j0 * E + 8 * h * E * E
        beta1 = 4 * rho * (2 * D * j0 + 16 * h * E)
        beta2 = Fraction(9, 2) * rho * rho * (16 * h)
        return BallFactorBound(_upper(max(beta0, beta1, beta2)), 0.0)

    def hstar_coefficient_ball(self) -> BallFactorBound:
        """Finite-derivative coefficient bound for H_*=Deta+dU_*.

        H_*=j0+(D+4)eta-j0*eta^2-4eta^3.  The alpha=0 source
        weights are a_{0,beta}=rho^-beta beta!/(beta+1)^2.
        """
        rho = Fraction.from_float(self.rho)
        margin = Fraction.from_float(float(self.domain.enlarged_real_margin))
        h = Fraction.from_float(float(self.domain.h))
        D = Fraction(1, 2) - h
        j0 = Fraction.from_float(float(self.domain.j0))
        if min(rho, h, D, j0) < 0:
            raise ValueError("source-compatible parameters must be nonnegative")
        E = Fraction(1) + margin
        beta0 = j0 + (D + 4) * E + j0 * E * E + 4 * E**3
        beta1 = 4 * rho * ((D + 4) + 2 * j0 * E + 12 * E * E)
        beta2 = Fraction(9, 2) * rho * rho * (2 * j0 + 24 * E)
        beta3 = Fraction(8, 3) * rho**3 * 24
        return BallFactorBound(_upper(max(beta0, beta1, beta2, beta3)), 0.0)

    def _j1_plain(self, bound: BallFactorBound) -> BallFactorBound:
        # Source J_1 vanishing-order b=0 factor: 80/((1)(1))=80.
        return _scale_exact(bound, Fraction(80))

    def _j1_logradial(self, bound: BallFactorBound) -> BallFactorBound:
        """Bound J1[Y d_Y F] directly from the source weight ratio.

        For alpha>=1 the exact output/input weight factor is at most 40 and
        therefore at most the shared source constant 80.  No standalone
        ``Y d_Y F`` norm is introduced.
        """
        return _scale_exact(bound, Fraction(80))

    def _j1_eta_derivative(self, bound: BallFactorBound) -> BallFactorBound:
        """Bound J1[partial_eta F] by the source 80/rho weight loss."""
        rho = Fraction.from_float(self.rho)
        return _scale_exact(bound, Fraction(80, 1) / rho)

    def _j1_average_times_logradial_u(self) -> BallFactorBound:
        """Bound J1[A(u) * Y u_Y] without a standalone Y u_Y norm.

        Use Y d_Y A(u)=u-A(u), hence

          A(u) Y u_Y = Y d_Y(A(u)u) - u^2 + A(u)u.

        Each product is bounded in coefficient space, then the source J1 or
        J1-logradial bound is applied exactly once.
        """
        u = self.bridge.u_source.u_ball()
        Au = self.radial_average_u_ball()
        Au_u = _product_exact(Au, u)
        u2 = _product_exact(u, u)
        return _add_exact(
            self._j1_logradial(Au_u),
            self._j1_plain(u2),
            self._j1_plain(Au_u),
        )

    def _j1_u_times_u_eta(self) -> BallFactorBound:
        """Bound J1[u u_eta]=(1/2) J1[partial_eta(u^2)]."""
        u = self.bridge.u_source.u_ball()
        u2 = _product_exact(u, u)
        return _scale_exact(self._j1_eta_derivative(u2), Fraction(1, 2))

    def r2_derivative_after_j1(self) -> dict[str, BallFactorBound]:
        """Return exactly the four missing source-compatible post-J1 R2 slots."""
        u = self.bridge.u_source.u_ball()
        Linv = self.bridge.l_inverse_coefficient_ball()
        fixed = self.algebraic.fixed_multiplier_balls()
        Wstar = self.wstar_coefficient_ball()
        Hstar = self.hstar_coefficient_ball()
        D = Fraction(1, 2) - Fraction.from_float(float(self.domain.h))

        w_yuy = _product_exact(Linv, Wstar, self._j1_logradial(u))
        avg_yuy = _product_exact(
            Linv,
            fixed["eta"],
            self._j1_average_times_logradial_u(),
        )
        avg_yuy = _scale_exact(avg_yuy, 2 * D)
        h_ueta = _product_exact(Linv, Hstar, self._j1_eta_derivative(u))
        u_ueta = _product_exact(
            Linv,
            fixed["d"],
            self._j1_u_times_u_eta(),
        )

        return {
            "W_star_times_Y_u_Y": w_yuy,
            "lambda_inv_2D_eta_Au_times_Y_u_Y": avg_yuy,
            "H_star_times_u_eta": h_ueta,
            "lambda_inv_d_u_u_eta": u_ueta,
        }

    def all_r2_ordinary_after_j1(self) -> dict[str, BallFactorBound]:
        """Merge the seven upstream slots with these four, all already post-J1."""
        pressure_after = self.algebraic.pressure_slots.r2_pressure_after_j1()
        algebraic_after = self.algebraic.r2_algebraic_after_j1()
        derivative_after = self.r2_derivative_after_j1()
        out = {**algebraic_after, **derivative_after, **pressure_after}
        if set(out) != set(R2_ORDINARY_TERMS):
            missing = sorted(set(R2_ORDINARY_TERMS) - set(out))
            extra = sorted(set(out) - set(R2_ORDINARY_TERMS))
            raise RuntimeError(f"R2 post-J1 ledger mismatch; missing={missing}, extra={extra}")
        return {name: out[name] for name in R2_ORDINARY_TERMS}

    def full_r2_after_j1(self) -> BallFactorBound:
        """Add all 11 ordinary post-J1 slots and the upstream mixed post-J1 slot."""
        ordinary = self.all_r2_ordinary_after_j1()
        mixed = self.bridge.full_mixed_terms()[
            "lambda_inv_Linv_d_detaAu_Y_u_Y_after_J1"
        ]
        return _add_exact(*ordinary.values(), mixed)

    @property
    def truth_boundary(self) -> dict[str, bool]:
        return copy.deepcopy(_TRUTH_BOUNDARY)

    def report(self) -> dict[str, Any]:
        u = self.bridge.u_source.u_ball()
        Au = self.radial_average_u_ball()
        Wstar = self.wstar_coefficient_ball()
        Hstar = self.hstar_coefficient_ball()
        derivative = self.r2_derivative_after_j1()
        all_ordinary = self.all_r2_ordinary_after_j1()
        mixed = self.bridge.full_mixed_terms()[
            "lambda_inv_Linv_d_detaAu_Y_u_Y_after_J1"
        ]
        full = self.full_r2_after_j1()
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
                "algebraic_slots_receipt_sha256": self.algebraic.sha256,
                "pressure_slots_receipt_sha256": self.algebraic.pressure_slots.sha256,
                "post_J_bridge_receipt_sha256": self.bridge.sha256,
                "source_u_ball_receipt_sha256": self.bridge.u_source.sha256,
                "source_axis_receipt_sha256": self.domain.sha256,
            },
            "selected_source_compatible_inputs": {
                "rho": self.rho,
                "h": float(self.domain.h),
                "D": float(self.domain.D),
                "j0": float(self.domain.j0),
                "Lambda_inverse_upper": 1.0,
                "product_constant_upper": self.product_constant_upper,
                "choice_status": (
                    "repository-autonomous source-compatible existence choice; "
                    "not recovered hidden source data"
                ),
            },
            "fixed_multiplier_bounds": {
                "L_inverse": asdict(self.bridge.l_inverse_coefficient_ball()),
                "W_star": asdict(Wstar),
                "H_star": asdict(Hstar),
                "eta": asdict(self.algebraic.fixed_multiplier_balls()["eta"]),
                "d": asdict(self.algebraic.fixed_multiplier_balls()["d"]),
            },
            "source_compatible_radius_one_u_ball": asdict(u),
            "source_compatible_radial_average_u_ball": asdict(Au),
            "R2_derivative_ordinary_after_J1": {
                name: asdict(bound) for name, bound in derivative.items()
            },
            "R2_all_ordinary_after_J1": {
                name: asdict(bound) for name, bound in all_ordinary.items()
            },
            "R2_mixed_after_J1": asdict(mixed),
            "R2_full_after_J1_agent1_self_certificate": asdict(full),
            "ledger_progress": {
                "R2_pressure_slots_source_bound": list(R2_PRESSURE_ORDINARY_TERMS),
                "R2_algebraic_slots_source_bound": list(R2_ALGEBRAIC_ORDINARY_TERMS),
                "R2_derivative_slots_source_bound": list(R2_DERIVATIVE_ORDINARY_TERMS),
                "R2_source_bound_ordinary_slots_count": len(all_ordinary),
                "R2_total_ordinary_slots": len(R2_ORDINARY_TERMS),
                "R2_ordinary_slots_still_required": 0,
                "R1_ordinary_slots_still_required": 11,
                "mixed_R2_slot_handled_upstream_by_A1_668": True,
                "full_R2_agent1_self_certificate_materialized": True,
                "full_R2_independent_admission": False,
                "full_R1_still_open": True,
            },
            "independent_audit_boundary": {
                "A4_source_axis_repair_PR": 671,
                "A4_source_axis_independent_PASS_exists": True,
                "A4_pressure_slots_audit_PR": 678,
                "A4_pressure_slots_independent_PASS_exists": True,
                "A4_source_u_ball_audit_PR": 687,
                "A4_source_u_ball_dedicated_run": 35462565728,
                "A4_source_u_ball_independent_PASS_exists": True,
                "A1_algebraic_slots_PR": 684,
                "A1_algebraic_slots_still_requires_independent_A4_audit": True,
                "this_derivative_subset_is_agent1_self_certificate": True,
                "this_derivative_subset_requires_independent_A4_audit": True,
                "full_R2_independent_A4_admission": False,
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
    print(json.dumps(KokunoPA10SourceDerivativeOrdinarySlots().save_report(args.output), indent=2, sort_keys=True))


if __name__ == "__main__":
    _main()
