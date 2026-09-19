"""Source-compatible PA.10 chi/Phi radius-one-ball bounds.

Pinned public provenance is KokunoYumeto/yang-mills-interacting-workbench at
``143f6773feb424ad9ed3a8d116653200f20346b7``, file
``navier-stokes/navier_stokes_workbench.tex`` (corrected reader dated
2026-09-09, Zenodo 22678406).

The corrected reconstruction uses the coefficient norm

    a_{alpha,beta} = 20^{-alpha} rho^{-beta} beta!
                     binom(alpha+beta,beta)
                     / ((alpha+1)^2 (beta+1)^2),

sets ``T=J_2 chi/2``, and proves

    ||T^k|| <= (40 M_chi)^k / (k! (k+1)!),
    Phi_0 = (1+T)^(-1) 1.

It then takes the closed radius-one ball about ``(Phi_0,u_0)``.  Agent-1
#644 supplies one explicit *source-compatible existence choice* for PA.8 and a
complex tube/Cauchy radius/rho; it is not recovered hidden source data.  This
module consumes only that certified domain and turns it into:

* a complex-supremum and coefficient multiplier bound for ``chi``;
* an exact-positive-series upper bound for ``||(1+T)^(-1)||``;
* a coefficient norm bound for ``Phi_0``;
* the induced norm/Lipschitz envelope for the Phi component on the source
  radius-one contraction ball;
* completed inputs for the already-existing conditional PA.10 pressure-ball
  calculator.

The complex ``chi`` estimate avoids the very loose numerator/denominator
quotient.  If ``a=H(eta)`` is real, ``w=H(z)`` and ``|w-a|<=epsilon<sigma``, then
for either sign

    |w|/|w +/- i sigma|
      <= (t+e)/(sqrt(t^2+1)-e),
      t=|a|/sigma, e=epsilon/sigma.

The right side has exact maximum ``(1+e^2)/(1-e^2)``.  Therefore

    sup |chi| <= ((1+e^2)/(1-e^2))^2.

All bound arithmetic introduced here is performed with ``Fraction`` values of
the configured binary64 inputs and converted back only with outward rounding.
This remains a reconstruction/source-compatible certificate, not a unique
Kokuno parameter recovery, paper-exact/OpenAI-field claim, or PDE validation.
The mixed ``Y Phi_Y`` / ``Y u_Y`` seams and hence R1/R2 -> M,K remain open.
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
from typing import Any

from .kokuno_pa10_pressure_ball_bounds import KokunoPA10PressureBallBounds
from .kokuno_pa10_remainder_ball_bounds import BallFactorBound
from .kokuno_pa10_source_axis_domain import KokunoPA10SourceCompatibleAxisDomain


SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
CORRECTED_RELEASE = "zenodo:22678406"
CORRECTED_RELEASE_DATE = "2026-09-09"
SCHEMA = "kokuno-pa10-source-compatible-phi-ball-v1"

_SOURCE_FORMULAS = {
    "coefficient_weight": (
        "a_{alpha,beta}=20^{-alpha}rho^{-beta}beta!binom(alpha+beta,beta)/"
        "((alpha+1)^2(beta+1)^2)"
    ),
    "chi": "chi=H_*^2/(H_*^2+sigma_*^2)",
    "T": "T=J_2 chi/2, multiplication by chi preceding J_2",
    "T_power": "||T^k|| <= (40 M_chi)^k/(k!(k+1)!)",
    "Phi_center": "Phi_0=f_0(Y chi)=(1+T)^(-1)1",
    "contraction_ball": "closed radius-one ball about (Phi_0,u_0)",
    "pressure": "p=I(g^2 Phi^2)",
}

_TRUTH_BOUNDARY = {
    "public_reconstruction_source": True,
    "source_compatible_axis_domain_consumed": True,
    "source_compatible_chi_complex_sup_machine_bound": True,
    "source_compatible_chi_multiplier_norm_machine_bound": True,
    "source_compatible_inverse_one_plus_T_machine_bound": True,
    "source_compatible_Phi0_coefficient_norm_machine_bound": True,
    "source_compatible_Phi_radius_one_ball_norm_machine_bound": True,
    "source_compatible_Phi_radius_one_ball_lipschitz_machine_bound": True,
    "source_compatible_pressure_operator_inputs_complete": True,
    "source_compatible_pressure_ball_agent1_envelope_executable": True,
    "source_axis_domain_independent_agent4_admission_required": True,
    "source_Phi_ball_independent_agent4_admission_required": True,
    "old_selected_sigma_star_equals_source_compatible_choice": False,
    "old_selected_local_pressure_audit_transfers_to_new_sigma": False,
    "source_hidden_sigma_star_recovered": False,
    "source_unique_rho_recovered": False,
    "source_mixed_Y_Phi_Y_radius_one_ball_bound_machine_bound": False,
    "source_mixed_Y_u_Y_radius_one_ball_bound_machine_bound": False,
    "source_R1_radius_one_ball_norm_machine_bound": False,
    "source_R2_radius_one_ball_norm_machine_bound": False,
    "source_R1_radius_one_ball_lipschitz_machine_bound": False,
    "source_R2_radius_one_ball_lipschitz_machine_bound": False,
    "source_operator_constant_M_machine_bound": False,
    "source_operator_constant_K_machine_bound": False,
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


@dataclass(frozen=True)
class KokunoPA10SourcePhiBallBounds:
    """Bind the Phi component of the source-compatible PA.10 unit ball."""

    domain: KokunoPA10SourceCompatibleAxisDomain = field(
        default_factory=KokunoPA10SourceCompatibleAxisDomain
    )

    def __post_init__(self) -> None:
        if not isinstance(self.domain, KokunoPA10SourceCompatibleAxisDomain):
            raise TypeError("domain must be KokunoPA10SourceCompatibleAxisDomain")
        complex_cert = self.domain.complex_domain_certificate()
        if not complex_cert["source_complex_neighborhood_machine_certified"]:
            raise ValueError("source-compatible complex-domain certificate is not valid")
        if not self.domain.real_interval_certificate()["source_PA8_implication_machine_certified"]:
            raise ValueError("source-compatible PA.8 certificate is not valid")
        variation = _q_nonnegative(
            complex_cert["tube_H_variation_abs_upper"], "tube_H_variation_abs_upper"
        )
        sigma = _q_nonnegative(self.domain.sigma_star, "sigma_star")
        if sigma <= 0 or variation >= sigma:
            raise ValueError("chi complex-sup bound requires tube H variation < sigma_star")

    def _chi_exact_bounds(self) -> dict[str, Fraction]:
        cert = self.domain.complex_domain_certificate()
        sigma = _q_nonnegative(self.domain.sigma_star, "sigma_star")
        variation = _q_nonnegative(
            cert["tube_H_variation_abs_upper"], "tube_H_variation_abs_upper"
        )
        e = variation / sigma
        e2 = e * e
        factor_ratio = (1 + e2) / (1 - e2)
        chi_complex_sup = factor_ratio * factor_ratio

        rho = _q_nonnegative(self.domain.coefficient_rho, "coefficient_rho")
        cauchy = _q_nonnegative(self.domain.cauchy_radius, "cauchy_radius")
        if not 0 < rho < cauchy:
            raise ValueError("coefficient rho must lie strictly inside Cauchy radius")
        x = rho / cauchy
        cauchy_weight_sum = (1 + x) / (1 - x) ** 3
        multiplier_norm = chi_complex_sup * cauchy_weight_sum
        return {
            "relative_H_variation_e": e,
            "single_factor_ratio_upper": factor_ratio,
            "chi_complex_sup_upper": chi_complex_sup,
            "rho_over_cauchy_radius": x,
            "cauchy_weight_sum_upper": cauchy_weight_sum,
            "chi_multiplier_norm_upper": multiplier_norm,
        }

    def chi_certificate(self) -> dict[str, float]:
        exact = self._chi_exact_bounds()
        return {name: _upper(value) for name, value in exact.items()}

    @staticmethod
    def _inverse_series_upper_exact(
        multiplier_norm_upper: Fraction,
        *,
        relative_tail: Fraction = Fraction(1, 10**18),
        max_terms: int = 4096,
    ) -> tuple[Fraction, int]:
        """Exact-rational positive-series upper for ``||(1+T)^-1||``."""
        if multiplier_norm_upper < 0:
            raise ValueError("multiplier norm must be nonnegative")
        if relative_tail <= 0:
            raise ValueError("relative_tail must be positive")
        if max_terms < 2:
            raise ValueError("max_terms must be >=2")
        if multiplier_norm_upper == 0:
            return Fraction(1), 0

        x = 40 * multiplier_norm_upper
        total = Fraction(1)
        term = Fraction(1)
        for k in range(1, max_terms + 1):
            term *= x / (k * (k + 1))
            total += term
            next_ratio = x / ((k + 1) * (k + 2))
            if next_ratio < 1:
                next_term = term * next_ratio
                tail = next_term / (1 - next_ratio)
                if tail <= relative_tail * total:
                    return total + tail, k
        raise RuntimeError("positive inverse series failed to obtain a finite tail bound")

    def phi_ball_certificate(self) -> dict[str, Any]:
        chi = self._chi_exact_bounds()
        inverse, terms = self._inverse_series_upper_exact(
            chi["chi_multiplier_norm_upper"]
        )
        # The source coefficient weight has a_{0,0}=1, hence ||1||_rho=1.
        phi0 = inverse
        phi_ball = phi0 + 1
        # Coordinate projection from the product contraction ball to Phi is
        # 1-Lipschitz for the repository/source pair norm used by the M/K bridge.
        return {
            "inverse_one_plus_T_absolute_series_upper": _upper(inverse),
            "inverse_series_terms_before_tail_majorant": int(terms),
            "Phi0_coefficient_norm_upper": _upper(phi0),
            "radius_one_ball_increment": 1.0,
            "Phi_radius_one_ball_norm_upper": _upper(phi_ball),
            "Phi_radius_one_ball_lipschitz_upper": 1.0,
        }

    def completed_pressure_operator_inputs(self) -> dict[str, Any]:
        axis_inputs = self.domain.pressure_operator_inputs()
        phi = self.phi_ball_certificate()
        return {
            "rho": float(axis_inputs["rho"]),
            "g_norm_upper": float(axis_inputs["g_norm_upper"]),
            "g_normalization_contract": str(axis_inputs["g_normalization_contract"]),
            "Phi_radius_one_ball_norm": float(phi["Phi_radius_one_ball_norm_upper"]),
            "Phi_radius_one_ball_lipschitz": float(
                phi["Phi_radius_one_ball_lipschitz_upper"]
            ),
            "ready_for_source_compatible_pressure_operator": True,
            "mixed_Y_Phi_Y_still_required_for_R1": True,
            "mixed_Y_u_Y_still_required_for_R2": True,
        }

    def pressure_ball_envelope(self) -> dict[str, dict[str, float]]:
        """Execute #635 pressure algebra with the now-complete Phi inputs."""
        inputs = self.completed_pressure_operator_inputs()
        calculator = KokunoPA10PressureBallBounds()
        pressure = calculator.pressure_bounds(
            rho=inputs["rho"],
            g_norm=inputs["g_norm_upper"],
            Phi=BallFactorBound(
                norm=inputs["Phi_radius_one_ball_norm"],
                lipschitz=inputs["Phi_radius_one_ball_lipschitz"],
            ),
        )
        return {
            name: {"norm": float(bound.norm), "lipschitz": float(bound.lipschitz)}
            for name, bound in pressure.items()
        }

    @property
    def truth_boundary(self) -> dict[str, bool]:
        return copy.deepcopy(_TRUTH_BOUNDARY)

    def report(self) -> dict[str, Any]:
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
            "source_compatible_axis_identity": {
                "axis_receipt_sha256": self.domain.sha256,
                "sigma_star": self.domain.sigma_star,
                "coefficient_rho": self.domain.coefficient_rho,
                "cauchy_radius": self.domain.cauchy_radius,
                "complex_tube_radius": self.domain.complex_tube_radius,
                "choice_status": (
                    "repository-autonomous source-compatible existence choice; "
                    "not recovered hidden source data"
                ),
            },
            "chi_certificate": self.chi_certificate(),
            "phi_ball_certificate": self.phi_ball_certificate(),
            "completed_pressure_operator_inputs": self.completed_pressure_operator_inputs(),
            "source_compatible_pressure_ball_envelope": self.pressure_ball_envelope(),
            "integration_boundary": {
                "conditional_pressure_algebra_reused_from_A1_635": True,
                "source_axis_independent_A4_audit_pending": True,
                "source_Phi_ball_independent_A4_audit_required": True,
                "mixed_Y_Phi_Y_bound_still_required": True,
                "mixed_Y_u_Y_bound_still_required": True,
                "R1_R2_M_K_not_promoted": True,
                "old_sigma_0p5_local_pressure_interface_reused": False,
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
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    payload = KokunoPA10SourcePhiBallBounds().save_report(args.output)
    chi = payload["chi_certificate"]
    phi = payload["phi_ball_certificate"]
    pressure = payload["source_compatible_pressure_ball_envelope"]
    print("receipt_sha256=", payload["receipt_sha256"])
    print("chi_complex_sup_upper=", chi["chi_complex_sup_upper"])
    print("chi_multiplier_norm_upper=", chi["chi_multiplier_norm_upper"])
    print("Phi_ball_norm_upper=", phi["Phi_radius_one_ball_norm_upper"])
    print("pressure_p_eta=", pressure["p_eta"])


if __name__ == "__main__":
    _main()
