"""Termwise post-radial-inverse bridge for Kokuno PA.10 R1/R2.

Pinned public provenance: KokunoYumeto/yang-mills-interacting-workbench
@143f6773feb424ad9ed3a8d116653200f20346b7,
``navier-stokes/navier_stokes_workbench.tex``, corrected 2026-09-09 reader
(Zenodo 22678406).

The corrected reconstruction defines

    L R1=[W+h(1-2 eta U)+d u zeta_*] Phi
         + W Y Phi_Y + H_c Phi_eta,

    L R2=[A(1-4 eta U_*)+d U_*']u -2A eta Lambda^-1 u^2
         + W Y u_Y + H_* u_eta + d Lambda^-1 u u_eta
         -4A eta p + d p_eta -2 eta Y p_Y,

with

    W=W_*+Lambda^-1[-2D eta A(u)-d partial_eta A(u)].

The source explicitly isolates only two monomials carrying *both* an
unintegrated eta derivative and a logarithmic radial derivative:

    (partial_eta A u) Y Phi_Y,   (partial_eta A u) Y u_Y.

Agent-1 #660 made their post-J_nu operator algebra executable and #662 supplied
a source-compatible radius-one u ball.  Those receipts intentionally stopped
short of a full R1/R2 claim: they did not include the outer L^-1 factor and did
not assemble every ordinary monomial after J_2/J_1.

This module supplies that missing structural bridge without inventing the still
missing ordinary-term coefficient bounds.  It:

* records an exact exhaustive ledger of the ordinary numerator monomials after
  expanding W and U;
* requires a caller bound for every ordinary monomial, so omission fails closed;
* multiplies each ordinary numerator bound by a machine-bounded coefficient
  norm for L^-1 and then applies the source J_nu vanishing-order estimate;
* recomputes the two mixed terms with BOTH d and L^-1 as undifferentiated
  factors, using #660's derivative-loss estimate and #662's actual u ball;
* uses only Lambda^-1<=1, valid uniformly for the source regime Lambda>=1, so
  no circular or hidden numerical Lambda is introduced;
* optionally forms the fixed-point pair envelopes directly from the already
  post-J bounds, avoiding the double-J mistake that would result from feeding
  them to the older raw-R1/R2 helper.

The ordinary-term inputs in this increment are caller evidence.  The built-in
unit fixture exists only to execute every ledger slot in CI.  Therefore full
source R1/R2, M/K, the fixed point, B0/T_sh, PA.16, the global leading field,
and PDE validation remain false.
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
from typing import Any, Mapping

from .kokuno_pa10_remainder_ball_bounds import BallFactorBound
from .kokuno_pa10_source_u_ball import KokunoPA10SourceUBallBounds

SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
CORRECTED_RELEASE = "zenodo:22678406"
CORRECTED_RELEASE_DATE = "2026-09-09"
SCHEMA = "kokuno-pa10-post-j-remainder-bridge-v1"

R1_ORDINARY_TERMS = (
    "W_star_times_Phi",
    "lambda_inv_2D_eta_Au_times_Phi",
    "lambda_inv_d_detaAu_times_Phi",
    "h_times_Phi",
    "2h_eta_Ustar_times_Phi",
    "2h_lambda_inv_eta_u_times_Phi",
    "d_u_zeta_star_times_Phi",
    "W_star_times_Y_Phi_Y",
    "lambda_inv_2D_eta_Au_times_Y_Phi_Y",
    "H_star_times_Phi_eta",
    "lambda_inv_d_u_times_Phi_eta",
)

R2_ORDINARY_TERMS = (
    "A_times_u",
    "4A_eta_Ustar_times_u",
    "d_Ustar_eta_times_u",
    "2A_lambda_inv_eta_u_squared",
    "W_star_times_Y_u_Y",
    "lambda_inv_2D_eta_Au_times_Y_u_Y",
    "H_star_times_u_eta",
    "lambda_inv_d_u_u_eta",
    "4A_eta_p",
    "d_p_eta",
    "2eta_Y_p_Y",
)

_SOURCE_FORMULAS = {
    "R1": (
        "L R1=[W+h(1-2eta U)+d u zeta_*]Phi+W YPhi_Y+H_c Phi_eta"
    ),
    "R2": (
        "L R2=[A(1-4eta U_*)+dU_*']u-2Aeta Lambda^-1u^2+W Yu_Y+"
        "H_*u_eta+d Lambda^-1u u_eta-4Aeta p+d p_eta-2eta Yp_Y"
    ),
    "W": "W=W_*+Lambda^-1[-2D eta A(u)-d partial_eta A(u)]",
    "mixed_R1": "Lambda^-1 L^-1 d(partial_eta A u)Y Phi_Y",
    "mixed_R2": "Lambda^-1 L^-1 d(partial_eta A u)Y u_Y",
    "J_nu": "||J_nu F||_rho <= 80/((b+1)(b+nu)) ||F||_rho",
    "mixed_bound": (
        "||J_nu[(partial_eta F)(Y partial_Y G)]||_rho <= "
        "(80/rho) C_sq^2 ||F||_rho ||G||_rho before extra factors"
    ),
    "fixed_point_pair": "((1+T)^-1 J_2 R1/2, J_1 R2/2)",
}

_TRUTH_BOUNDARY = {
    "public_reconstruction_source": True,
    "exact_expanded_R1_R2_term_ledger_recorded": True,
    "termwise_post_J_remainder_bridge_executable": True,
    "source_compatible_L_inverse_coefficient_norm_machine_bound": True,
    "source_compatible_full_R1_mixed_post_J2_term_executable": True,
    "source_compatible_full_R2_mixed_post_J1_term_executable": True,
    "uniform_Lambda_inverse_upper_for_Lambda_ge_1_used": True,
    "standalone_Y_Phi_Y_coefficient_norm_invented": False,
    "standalone_Y_u_Y_coefficient_norm_invented": False,
    "ordinary_term_inputs_are_source_certified": False,
    "source_axis_domain_independent_agent4_admission_required": True,
    "source_u_ball_independent_agent4_admission_required": True,
    "source_full_post_J2_R1_radius_one_ball_bound_machine_bound": False,
    "source_full_post_J1_R2_radius_one_ball_bound_machine_bound": False,
    "source_R1_radius_one_ball_norm_machine_bound": False,
    "source_R2_radius_one_ball_norm_machine_bound": False,
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


def _add(*bounds: BallFactorBound) -> BallFactorBound:
    norm = sum((_q_nonnegative(item.norm, "norm") for item in bounds), Fraction(0))
    lip = sum(
        (_q_nonnegative(item.lipschitz, "lipschitz") for item in bounds), Fraction(0)
    )
    return BallFactorBound(_upper(norm), _upper(lip))


@dataclass(frozen=True)
class KokunoPA10PostJRemainderBridge:
    """Assemble exact PA.10 term slots after the required radial inverses."""

    u_source: KokunoPA10SourceUBallBounds = field(
        default_factory=KokunoPA10SourceUBallBounds
    )

    def __post_init__(self) -> None:
        if not isinstance(self.u_source, KokunoPA10SourceUBallBounds):
            raise TypeError("u_source must be KokunoPA10SourceUBallBounds")
        if not 0.0 < self.rho < self.cauchy_radius:
            raise ValueError("coefficient rho must lie inside the Cauchy radius")

    @property
    def domain(self):
        return self.u_source.domain

    @property
    def rho(self) -> float:
        return float(self.domain.coefficient_rho)

    @property
    def cauchy_radius(self) -> float:
        return float(self.domain.cauchy_radius)

    @staticmethod
    def _product_constant_exact() -> Fraction:
        # Elementary pi < 22/7; C_sq=4*pi^2/3 and product constant=C_sq^2.
        pi_upper = Fraction(22, 7)
        c_sq_upper = Fraction(4, 3) * pi_upper * pi_upper
        return c_sq_upper * c_sq_upper

    @staticmethod
    def _j_nu_factor_exact(nu: int, vanishing_order: int = 0) -> Fraction:
        if isinstance(nu, bool) or nu not in (1, 2):
            raise ValueError("nu must be 1 or 2")
        if isinstance(vanishing_order, bool) or not isinstance(vanishing_order, int):
            raise TypeError("vanishing_order must be an integer")
        if vanishing_order < 0:
            raise ValueError("vanishing_order must be nonnegative")
        b = int(vanishing_order)
        return Fraction(80, (b + 1) * (b + nu))

    def l_inverse_coefficient_ball(self) -> BallFactorBound:
        """Cauchy coefficient-norm bound for the fixed multiplier L^-1."""
        h = Fraction.from_float(float(self.domain.h))
        margin = Fraction.from_float(float(self.domain.enlarged_real_margin))
        radius = Fraction.from_float(self.cauchy_radius)
        rho = Fraction.from_float(self.rho)
        z_abs = Fraction(1) + margin + radius
        lower = Fraction(1) - 2 * h * z_abs * z_abs
        if lower <= 0:
            raise ValueError("Cauchy neighborhood does not keep L away from zero")
        x = rho / radius
        # sum_{beta>=0} (beta+1)^2 x^beta=(1+x)/(1-x)^3,
        # a conservative majorant for the degree-zero coefficient sup norm.
        cauchy_weight_sum = (1 + x) / (1 - x) ** 3
        norm = cauchy_weight_sum / lower
        return BallFactorBound(_upper(norm), 0.0)

    @staticmethod
    def _validate_term_map(
        values: Mapping[str, BallFactorBound], expected: tuple[str, ...], name: str
    ) -> dict[str, BallFactorBound]:
        if not isinstance(values, Mapping):
            raise TypeError(f"{name} must be a mapping")
        actual = set(values)
        wanted = set(expected)
        if actual != wanted:
            missing = sorted(wanted - actual)
            extra = sorted(actual - wanted)
            raise ValueError(f"{name} term ledger mismatch; missing={missing}, extra={extra}")
        out = dict(values)
        if not all(isinstance(out[key], BallFactorBound) for key in expected):
            raise TypeError(f"every {name} term must be a BallFactorBound")
        return out

    def ordinary_term_after_jnu(
        self,
        numerator_term: BallFactorBound,
        *,
        nu: int,
        vanishing_order: int = 0,
    ) -> BallFactorBound:
        """Apply fixed L^-1 and then the source J_nu norm factor.

        ``numerator_term`` is the complete coefficient-space bound for one
        expanded *ordinary* numerator monomial, including any displayed scalar
        such as h, A, D or Lambda^-1<=1.  It must not be one of the two mixed
        eta/log-radial seams.
        """
        if not isinstance(numerator_term, BallFactorBound):
            raise TypeError("numerator_term must be a BallFactorBound")
        Linv = self.l_inverse_coefficient_ball()
        product = self._product_constant_exact()
        j = self._j_nu_factor_exact(nu, vanishing_order)
        factor = product * j * _q_nonnegative(Linv.norm, "L inverse norm")
        return BallFactorBound(
            _upper(factor * _q_nonnegative(numerator_term.norm, "term norm")),
            _upper(
                factor * _q_nonnegative(numerator_term.lipschitz, "term lipschitz")
            ),
        )

    def full_mixed_terms(self) -> dict[str, BallFactorBound]:
        """Return the two complete mixed post-J terms, including L^-1.

        The source estimates are uniform for Lambda>=1.  We therefore use the
        exact inequality Lambda^-1<=1 rather than selecting a circular Lambda.
        """
        mixed = self.u_source.mixed
        u = self.u_source.u_ball()
        phi = mixed.source_phi_ball()
        d = mixed.d_coefficient_ball()
        Linv = self.l_inverse_coefficient_ball()
        r1 = mixed.averaged_eta_logradial_after_jnu(
            u=u, G=phi, nu=2, extra_factors=(d, Linv)
        )
        r2 = mixed.averaged_eta_logradial_after_jnu(
            u=u, G=u, nu=1, extra_factors=(d, Linv)
        )
        return {
            "lambda_inv_Linv_d_detaAu_Y_Phi_Y_after_J2": r1,
            "lambda_inv_Linv_d_detaAu_Y_u_Y_after_J1": r2,
        }

    def assemble(
        self,
        *,
        r1_ordinary_numerator_terms: Mapping[str, BallFactorBound],
        r2_ordinary_numerator_terms: Mapping[str, BallFactorBound],
        r1_vanishing_order: int = 0,
        r2_vanishing_order: int = 0,
    ) -> dict[str, Any]:
        """Assemble conditional full post-J2 R1 and post-J1 R2 envelopes."""
        r1_in = self._validate_term_map(
            r1_ordinary_numerator_terms, R1_ORDINARY_TERMS, "R1 ordinary"
        )
        r2_in = self._validate_term_map(
            r2_ordinary_numerator_terms, R2_ORDINARY_TERMS, "R2 ordinary"
        )
        r1_after = {
            name: self.ordinary_term_after_jnu(
                r1_in[name], nu=2, vanishing_order=r1_vanishing_order
            )
            for name in R1_ORDINARY_TERMS
        }
        r2_after = {
            name: self.ordinary_term_after_jnu(
                r2_in[name], nu=1, vanishing_order=r2_vanishing_order
            )
            for name in R2_ORDINARY_TERMS
        }
        mixed = self.full_mixed_terms()
        r1_total = _add(
            *r1_after.values(),
            mixed["lambda_inv_Linv_d_detaAu_Y_Phi_Y_after_J2"],
        )
        r2_total = _add(
            *r2_after.values(),
            mixed["lambda_inv_Linv_d_detaAu_Y_u_Y_after_J1"],
        )
        return {
            "R1_ordinary_after_J2": {name: asdict(value) for name, value in r1_after.items()},
            "R2_ordinary_after_J1": {name: asdict(value) for name, value in r2_after.items()},
            "mixed_full_after_Jnu": {name: asdict(value) for name, value in mixed.items()},
            "post_J2_R1": asdict(r1_total),
            "post_J1_R2": asdict(r2_total),
            "ordinary_inputs_source_certified": False,
        }

    def diagnostic_pair_envelope_from_post_j(self, assembly: Mapping[str, Any]) -> dict[str, float]:
        """Form the fixed-point pair directly from already-post-J envelopes.

        This is diagnostic until every ordinary term input is independently
        source-certified.  It intentionally does NOT call the older helper that
        would apply J_2/J_1 a second time.
        """
        r1 = assembly["post_J2_R1"]
        r2 = assembly["post_J1_R2"]
        phi_cert = self.u_source.mixed.phi_source.phi_ball_certificate()
        inverse = _q_nonnegative(
            phi_cert["inverse_one_plus_T_absolute_series_upper"],
            "inverse_one_plus_T_absolute_series_upper",
        )
        m_phi = Fraction(1, 2) * inverse * _q_nonnegative(r1["norm"], "post J2 R1 norm")
        m_u = Fraction(1, 2) * _q_nonnegative(r2["norm"], "post J1 R2 norm")
        k_phi = Fraction(1, 2) * inverse * _q_nonnegative(
            r1["lipschitz"], "post J2 R1 lipschitz"
        )
        k_u = Fraction(1, 2) * _q_nonnegative(r2["lipschitz"], "post J1 R2 lipschitz")
        return {
            "inverse_one_plus_T_absolute_series_upper": _upper(inverse),
            "diagnostic_M_phi_component": _upper(m_phi),
            "diagnostic_M_u_component": _upper(m_u),
            "diagnostic_M_pair_max": _upper(max(m_phi, m_u)),
            "diagnostic_K_phi_component": _upper(k_phi),
            "diagnostic_K_u_component": _upper(k_u),
            "diagnostic_K_pair_max": _upper(max(k_phi, k_u)),
        }

    @staticmethod
    def diagnostic_unit_ordinary_terms() -> tuple[dict[str, BallFactorBound], dict[str, BallFactorBound]]:
        """Arithmetic-only fixture; every ordinary slot receives (1,1)."""
        r1 = {name: BallFactorBound(1.0, 1.0) for name in R1_ORDINARY_TERMS}
        r2 = {name: BallFactorBound(1.0, 1.0) for name in R2_ORDINARY_TERMS}
        return r1, r2

    @property
    def truth_boundary(self) -> dict[str, bool]:
        return copy.deepcopy(_TRUTH_BOUNDARY)

    def report(self) -> dict[str, Any]:
        r1, r2 = self.diagnostic_unit_ordinary_terms()
        assembly = self.assemble(
            r1_ordinary_numerator_terms=r1,
            r2_ordinary_numerator_terms=r2,
        )
        pair = self.diagnostic_pair_envelope_from_post_j(assembly)
        Linv = self.l_inverse_coefficient_ball()
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
                "source_u_ball_receipt_sha256": self.u_source.sha256,
                "mixed_jnu_receipt_sha256": self.u_source.mixed.sha256,
                "source_phi_receipt_sha256": self.u_source.mixed.phi_source.sha256,
                "source_axis_receipt_sha256": self.domain.sha256,
            },
            "term_ledger": {
                "R1_ordinary_numerator_terms": list(R1_ORDINARY_TERMS),
                "R2_ordinary_numerator_terms": list(R2_ORDINARY_TERMS),
                "R1_mixed_term": "lambda_inv_Linv_d_detaAu_Y_Phi_Y_after_J2",
                "R2_mixed_term": "lambda_inv_Linv_d_detaAu_Y_u_Y_after_J1",
            },
            "source_compatible_fixed_inputs": {
                "rho": self.rho,
                "cauchy_radius": self.cauchy_radius,
                "L_inverse_coefficient_ball": asdict(Linv),
                "Lambda_inverse_uniform_upper_for_Lambda_ge_1": 1.0,
            },
            "source_compatible_full_mixed_terms": assembly["mixed_full_after_Jnu"],
            "diagnostic_unit_ordinary_fixture": {
                "status": "arithmetic-only; NOT source ordinary-term bounds",
                "assembly": assembly,
                "pair_envelope": pair,
            },
            "integration_boundary": {
                "full_mixed_terms_now_include_L_inverse": True,
                "Lambda_inverse_is_uniform_source_bound_not_selected_parameter": True,
                "ordinary_term_ledger_is_exhaustive_and_fail_closed": True,
                "ordinary_term_source_bounds_still_required": True,
                "post_J_pair_helper_does_not_apply_J_twice": True,
                "full_R1_R2_M_K_not_promoted": True,
                "independent_agent4_admission_required": True,
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
    print(json.dumps(KokunoPA10PostJRemainderBridge().save_report(args.output), indent=2, sort_keys=True))


if __name__ == "__main__":
    _main()
