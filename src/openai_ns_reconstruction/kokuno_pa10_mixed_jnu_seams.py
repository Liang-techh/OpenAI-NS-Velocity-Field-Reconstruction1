"""Executable post-J_nu bounds for Kokuno PA.10 mixed derivative seams.

Pinned provenance: KokunoYumeto/yang-mills-interacting-workbench
@143f6773feb424ad9ed3a8d116653200f20346b7,
navier-stokes/navier_stokes_workbench.tex, corrected 2026-09-09 reader.

The source isolates (d_eta A u) Y Phi_Y and (d_eta A u) Y u_Y and proves,
for nu=1,2,

 ||J_nu[(d_eta F)(Y d_Y G)]||_rho
 <= (80/rho) C_sq^2 ||F||_rho ||G||_rho.

If F is averaged, its (i+1)^-1 divisor cancels the added-degree factor.
Extra undifferentiated factors add ordinary product-convolution constants.
This module exposes that post-radial-inverse estimate; it never invents
standalone Y Phi_Y or Y u_Y bounds.  The u-ball remains caller evidence, so
R1/R2, M/K and PDE validation stay fail-closed.
"""

from __future__ import annotations

import argparse
import copy
from dataclasses import dataclass, field
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Iterable

from .kokuno_pa10_remainder_ball_bounds import BallFactorBound
from .kokuno_pa10_source_phi_ball import KokunoPA10SourcePhiBallBounds

SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
CORRECTED_RELEASE = "zenodo:22678406"
CORRECTED_RELEASE_DATE = "2026-09-09"
SCHEMA = "kokuno-pa10-mixed-jnu-seams-v1"

_SOURCE_FORMULAS = {
    "mixed": (
        "||J_nu[(partial_eta F)(Y partial_Y G)]||_rho <= "
        "(80/rho) C_sq^2 ||F||_rho ||G||_rho"
    ),
    "averaged_F": "A F=Y^-1 I F; divisor i+1 cancels added-degree factor",
    "extra_factors": "each undifferentiated factor adds product convolution",
    "R1_seam": "(partial_eta A u) Y Phi_Y",
    "R2_seam": "(partial_eta A u) Y u_Y",
    "d": "d=1-eta^2",
}

_TRUTH_BOUNDARY = {
    "public_reconstruction_source": True,
    "source_post_Jnu_mixed_derivative_formula_recorded": True,
    "source_post_Jnu_mixed_derivative_bound_executable": True,
    "source_averaged_factor_cancellation_executable": True,
    "source_compatible_d_coefficient_norm_machine_bound": True,
    "source_compatible_Phi_ball_consumed": True,
    "standalone_Y_Phi_Y_coefficient_norm_invented": False,
    "standalone_Y_u_Y_coefficient_norm_invented": False,
    "source_u_radius_one_ball_norm_machine_bound": False,
    "source_u_radius_one_ball_lipschitz_machine_bound": False,
    "source_post_J2_R1_mixed_numeric_ball_bound_complete": False,
    "source_post_J1_R2_mixed_numeric_ball_bound_complete": False,
    "source_R1_radius_one_ball_norm_machine_bound": False,
    "source_R2_radius_one_ball_norm_machine_bound": False,
    "source_operator_constant_M_machine_bound": False,
    "source_operator_constant_K_machine_bound": False,
    "source_fixed_point_distance_machine_bound": False,
    "source_B0_dependencies_machine_bound": False,
    "source_T_sh_lower_bound_verified": False,
    "selected_pa16_handoff_allowed": False,
    "global_pressure_matched": False,
    "global_leading_profile_reconstructed": False,
    "heldout_ns_residual_assessed": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


def _canonical_json(payload: dict[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _q(value: float, name: str) -> Fraction:
    out = float(value)
    if not math.isfinite(out) or out < 0.0:
        raise ValueError(f"{name} must be finite and nonnegative")
    return Fraction.from_float(out)


def _upper(value: Fraction) -> float:
    if value < 0:
        raise ValueError("upper bound must be nonnegative")
    out = float(value)
    if not math.isfinite(out):
        raise OverflowError("upper bound left binary64 range")
    if Fraction.from_float(out) < value:
        out = math.nextafter(out, math.inf)
    return float(out)


@dataclass(frozen=True)
class KokunoPA10MixedJnuSeamBounds:
    """Post-J_nu mixed estimate with explicit caller-owned u ball."""

    phi_source: KokunoPA10SourcePhiBallBounds = field(
        default_factory=KokunoPA10SourcePhiBallBounds
    )

    def __post_init__(self) -> None:
        if not isinstance(self.phi_source, KokunoPA10SourcePhiBallBounds):
            raise TypeError("phi_source must be KokunoPA10SourcePhiBallBounds")
        if self.rho <= 0.0:
            raise ValueError("coefficient rho must be positive")

    @property
    def rho(self) -> float:
        return float(self.phi_source.domain.coefficient_rho)

    @staticmethod
    def _product_constant_exact() -> Fraction:
        pi_upper = Fraction(22, 7)  # elementary pi < 22/7
        c_sq_upper = Fraction(4, 3) * pi_upper * pi_upper
        return c_sq_upper * c_sq_upper

    @property
    def product_constant_upper(self) -> float:
        return _upper(self._product_constant_exact())

    def _mixed_factor_exact(self, extra_count: int) -> Fraction:
        if isinstance(extra_count, bool) or not isinstance(extra_count, int):
            raise TypeError("extra_count must be an integer")
        if extra_count < 0:
            raise ValueError("extra_count must be nonnegative")
        return (
            Fraction(80) / Fraction.from_float(self.rho)
            * self._product_constant_exact() ** (1 + extra_count)
        )

    def mixed_operator_factor(self, extra_factor_count: int = 0) -> float:
        return _upper(self._mixed_factor_exact(extra_factor_count))

    @staticmethod
    def _check_nu(nu: int) -> None:
        if isinstance(nu, bool) or not isinstance(nu, int):
            raise TypeError("nu must be an integer")
        if nu not in (1, 2):
            raise ValueError("nu must be 1 or 2")

    def averaged_eta_logradial_after_jnu(
        self,
        *,
        u: BallFactorBound,
        G: BallFactorBound,
        nu: int,
        extra_factors: Iterable[BallFactorBound] = (),
    ) -> BallFactorBound:
        """Bound J_nu[(d_eta A u)(Y d_Y G)*extras] after radial inversion."""
        self._check_nu(nu)
        if not isinstance(u, BallFactorBound) or not isinstance(G, BallFactorBound):
            raise TypeError("u and G must be BallFactorBound values")
        extras = tuple(extra_factors)
        if not all(isinstance(x, BallFactorBound) for x in extras):
            raise TypeError("extra_factors must be BallFactorBound values")
        factors = (u, G, *extras)
        constant = self._mixed_factor_exact(len(extras))
        norm = constant
        for factor in factors:
            norm *= _q(factor.norm, "factor norm")
        lip = Fraction(0)
        for i, factor in enumerate(factors):
            if factor.lipschitz == 0.0:
                continue
            term = constant * _q(factor.lipschitz, "factor lipschitz")
            for j, other in enumerate(factors):
                if i != j:
                    term *= _q(other.norm, "factor norm")
            lip += term
        return BallFactorBound(_upper(norm), _upper(lip))

    def d_coefficient_ball(self) -> BallFactorBound:
        """Coefficient-norm certificate for fixed d=1-eta^2 on #644 domain."""
        domain = self.phi_source.domain
        rho = Fraction.from_float(float(domain.coefficient_rho))
        E = Fraction(1) + Fraction.from_float(float(domain.enlarged_real_margin))
        # alpha=0. Only beta=0,1,2 are nonzero.
        beta0 = max(Fraction(1), E * E - 1)
        beta1 = 8 * E * rho            # |d'| / a_{0,1}
        beta2 = 9 * rho * rho          # |d''| / a_{0,2}
        return BallFactorBound(_upper(max(beta0, beta1, beta2)), 0.0)

    def source_phi_ball(self) -> BallFactorBound:
        cert = self.phi_source.phi_ball_certificate()
        return BallFactorBound(
            float(cert["Phi_radius_one_ball_norm_upper"]),
            float(cert["Phi_radius_one_ball_lipschitz_upper"]),
        )

    def r1_mixed_after_j2(self, *, u: BallFactorBound) -> BallFactorBound:
        """Conditional bound for J_2[d(d_eta A u)Y Phi_Y]."""
        return self.averaged_eta_logradial_after_jnu(
            u=u,
            G=self.source_phi_ball(),
            nu=2,
            extra_factors=(self.d_coefficient_ball(),),
        )

    def r2_mixed_after_j1(self, *, u: BallFactorBound) -> BallFactorBound:
        """Conditional bound for J_1[d(d_eta A u)Y u_Y]."""
        return self.averaged_eta_logradial_after_jnu(
            u=u,
            G=u,
            nu=1,
            extra_factors=(self.d_coefficient_ball(),),
        )

    @property
    def truth_boundary(self) -> dict[str, bool]:
        return copy.deepcopy(_TRUTH_BOUNDARY)

    def report(self) -> dict[str, Any]:
        unit_u = BallFactorBound(1.0, 1.0)
        d = self.d_coefficient_ball()
        phi = self.source_phi_ball()
        r1 = self.r1_mixed_after_j2(u=unit_u)
        r2 = self.r2_mixed_after_j1(u=unit_u)
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
            "source_compatible_inputs": {
                "source_phi_receipt_sha256": self.phi_source.sha256,
                "source_axis_receipt_sha256": self.phi_source.domain.sha256,
                "rho": self.rho,
                "product_constant_upper": self.product_constant_upper,
                "mixed_operator_factor_no_extra": self.mixed_operator_factor(),
                "mixed_operator_factor_one_extra": self.mixed_operator_factor(1),
                "d_ball": {"norm": d.norm, "lipschitz": d.lipschitz},
                "Phi_ball": {"norm": phi.norm, "lipschitz": phi.lipschitz},
            },
            "diagnostic_unit_u_fixture": {
                "status": "arithmetic-only; not a source u-ball certificate",
                "post_J2_d_detaAu_YPhiY": {"norm": r1.norm, "lipschitz": r1.lipschitz},
                "post_J1_d_detaAu_YuY": {"norm": r2.norm, "lipschitz": r2.lipschitz},
            },
            "integration_boundary": {
                "mixed_terms_bounded_after_Jnu": True,
                "standalone_Y_Phi_Y_bound_required": False,
                "standalone_Y_u_Y_bound_required": False,
                "source_u_ball_still_required": True,
                "pre_J_R1_R2_calculator_not_relabelled_as_closed": True,
                "future_termwise_post_J_R1_R2_bridge_required": True,
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
    print(json.dumps(KokunoPA10MixedJnuSeamBounds().save_report(args.output), indent=2, sort_keys=True))


if __name__ == "__main__":
    _main()
