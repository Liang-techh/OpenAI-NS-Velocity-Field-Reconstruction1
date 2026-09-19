"""Source-compatible PA.10 u_0 and radius-one u-ball certificate.

Pinned public provenance is KokunoYumeto/yang-mills-interacting-workbench at
``143f6773feb424ad9ed3a8d116653200f20346b7``, file
``navier-stokes/navier_stokes_workbench.tex`` (corrected reader dated
2026-09-09, Zenodo 22678406).

The corrected reconstruction centers the PA.10 contraction at

    Phi_0 = f_0(Y chi),
    u_0   = -Y Z_*/(2 L),

and takes the closed radius-one ball about ``(Phi_0,u_0)`` in the displayed
coefficient norm

    a_{alpha,beta}=20^{-alpha} rho^{-beta} beta!
                    binom(alpha+beta,beta)/((alpha+1)^2(beta+1)^2).

For ``u_0=Y f(eta)`` only alpha=1 is present, so

    a_{1,beta}=beta! rho^{-beta}/(80(beta+1)).

On the #644 certified complex tube, Cauchy's estimate at radius ``r`` gives

    ||u_0||_rho <= 80 sup|f| sup_beta (beta+1)(rho/r)^beta
                 <= 80 sup|f|/(1-rho/r)^2.

This module bounds ``f=-Z_*/(2L)`` directly from the public rational axis
formulas on the certified Cauchy neighborhood, adds the source radius-one ball
increment, and feeds the resulting u envelope into Agent-1 #660's post-J_nu
mixed-seam operators.  It therefore removes the diagnostic unit-u placeholder
from those two mixed terms without inventing standalone ``Y Phi_Y`` or
``Y u_Y`` norms.

The selected sigma/rho/pressure datum remain repository-autonomous
*source-compatible existence choices*, not recovered hidden Kokuno/OpenAI
parameters.  A termwise post-J_nu bridge for all R1/R2 terms, independent
Agent-4 admission, M/K closure, global matched pressure/velocity and held-out
PDE validation remain open.
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

import numpy as np

from .kokuno_pa10_mixed_jnu_seams import KokunoPA10MixedJnuSeamBounds
from .kokuno_pa10_remainder_ball_bounds import BallFactorBound


SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
CORRECTED_RELEASE = "zenodo:22678406"
CORRECTED_RELEASE_DATE = "2026-09-09"
SCHEMA = "kokuno-pa10-source-compatible-u-ball-v1"

_SOURCE_FORMULAS = {
    "axis_center": "u_0(Y,eta)=-Y Z_*(eta)/(2L(eta))",
    "L": "L=1-2h eta^2",
    "U_star": "U_*=4eta+j_0",
    "H_star": "H_*=Deta+(1-eta^2)U_*",
    "Z_star": (
        "Z_*=-A(1-2etaU_*)U_*-H_*U_*'-(1-eta^2)Pi_0'"
        "+4AetaPi_0; U_*'=4"
    ),
    "Pi_0": "Pi_0=-P^2/(1+eta^2)^2 on the selected source-compatible datum",
    "coefficient_weight_alpha1": (
        "a_{1,beta}=beta! rho^{-beta}/(80(beta+1))"
    ),
    "cauchy_u0": (
        "||u_0||_rho <= 80 sup|Z_*/(2L)|/(1-rho/r)^2"
    ),
    "contraction_ball": "closed radius-one ball about (Phi_0,u_0)",
    "mixed_R1": "J_2[d(partial_eta A u)Y Phi_Y]",
    "mixed_R2": "J_1[d(partial_eta A u)Y u_Y]",
}

_TRUTH_BOUNDARY = {
    "public_reconstruction_source": True,
    "source_compatible_axis_domain_consumed": True,
    "source_compatible_u0_formula_executable": True,
    "source_compatible_u0_coefficient_norm_machine_bound": True,
    "source_compatible_u_radius_one_ball_norm_machine_bound": True,
    "source_compatible_u_radius_one_ball_lipschitz_machine_bound": True,
    "source_compatible_post_J2_R1_mixed_numeric_ball_bound_executable": True,
    "source_compatible_post_J1_R2_mixed_numeric_ball_bound_executable": True,
    "standalone_Y_Phi_Y_coefficient_norm_invented": False,
    "standalone_Y_u_Y_coefficient_norm_invented": False,
    "source_axis_domain_independent_agent4_admission_required": True,
    "source_Phi_ball_independent_agent4_admission_required": True,
    "source_u_ball_independent_agent4_admission_required": True,
    "source_hidden_u0_parameters_recovered": False,
    "source_unique_sigma_rho_pressure_datum_recovered": False,
    "source_full_post_J2_R1_radius_one_ball_bound_machine_bound": False,
    "source_full_post_J1_R2_radius_one_ball_bound_machine_bound": False,
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


def _q_positive(value: float, name: str) -> Fraction:
    out = _q_nonnegative(value, name)
    if out <= 0:
        raise ValueError(f"{name} must be positive")
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


@dataclass(frozen=True)
class KokunoPA10SourceUBallBounds:
    """Bind the u component of the selected source-compatible PA.10 unit ball."""

    mixed: KokunoPA10MixedJnuSeamBounds = field(
        default_factory=KokunoPA10MixedJnuSeamBounds
    )

    def __post_init__(self) -> None:
        if not isinstance(self.mixed, KokunoPA10MixedJnuSeamBounds):
            raise TypeError("mixed must be KokunoPA10MixedJnuSeamBounds")
        domain = self.domain
        if not domain.real_interval_certificate()["source_PA8_implication_machine_certified"]:
            raise ValueError("source-compatible PA.8 certificate is not valid")
        if not domain.complex_domain_certificate()["source_complex_neighborhood_machine_certified"]:
            raise ValueError("source-compatible complex-domain certificate is not valid")
        if not 0.0 < domain.coefficient_rho < domain.cauchy_radius:
            raise ValueError("coefficient rho must lie inside the Cauchy radius")

    @property
    def domain(self):
        return self.mixed.phi_source.domain

    def center_values(self, Y: Any, eta: Any) -> dict[str, np.ndarray]:
        """Vectorized source-center ``u_0`` and its exact radial derivative.

        This is the displayed contraction center on the physical inner rectangle
        ``0<=Y<=4.1`` and the #644 certified enlarged real eta interval.
        """
        Y_arr, eta_arr = np.broadcast_arrays(
            np.asarray(Y, dtype=float), np.asarray(eta, dtype=float)
        )
        if np.any(~np.isfinite(Y_arr)) or np.any(~np.isfinite(eta_arr)):
            raise ValueError("Y and eta must be finite")
        if np.any((Y_arr < 0.0) | (Y_arr > 4.1)):
            raise ValueError("Y must lie in [0,4.1]")
        state = self.domain.axis_state(eta_arr)
        slope = -state["Z_star"] / (2.0 * state["L"])
        u0 = Y_arr * slope
        return {
            "u_0": u0,
            "u_0_Y": slope,
            "Y_u_0_Y": u0,
            "Z_star": state["Z_star"],
            "L": state["L"],
        }

    def _complex_u0_exact_bounds(self) -> dict[str, Fraction]:
        domain = self.domain
        h = _q_positive(domain.h, "h")
        A = Fraction(1, 2) + h
        D = Fraction(1, 2) - h
        j0 = _q_positive(domain.j0, "j0")
        pressure_square = _q_positive(domain.pressure_square, "pressure_square")
        margin = _q_positive(domain.enlarged_real_margin, "enlarged_real_margin")
        radius = _q_positive(domain.cauchy_radius, "cauchy_radius")
        rho = _q_positive(domain.coefficient_rho, "coefficient_rho")
        if radius >= 1:
            raise ValueError("Cauchy radius must be below the Pi_0 pole distance")
        if rho >= radius:
            raise ValueError("coefficient rho must be smaller than Cauchy radius")

        z_abs = Fraction(1) + margin + radius
        d_abs = Fraction(1) + z_abs * z_abs
        L_lower = Fraction(1) - 2 * h * z_abs * z_abs
        if L_lower <= 0:
            raise ValueError("Cauchy neighborhood does not keep L away from zero")

        U_abs = 4 * z_abs + j0
        B_abs = Fraction(1) + 2 * z_abs * U_abs
        H_abs = abs(D) * z_abs + d_abs * U_abs

        # 1+z^2=(z-i)(z+i).  Every real center has distance >=1 from
        # either pole, so a radius-r perturbation keeps each factor >=1-r.
        pole_factor_lower = Fraction(1) - radius
        one_plus_z2_lower = pole_factor_lower * pole_factor_lower
        Pi_abs = pressure_square / (one_plus_z2_lower**2)
        Pi_eta_abs = 4 * pressure_square * z_abs / (one_plus_z2_lower**3)

        Z_abs = (
            A * B_abs * U_abs
            + 4 * H_abs
            + d_abs * Pi_eta_abs
            + 4 * A * z_abs * Pi_abs
        )
        slope_abs = Z_abs / (2 * L_lower)

        x = rho / radius
        cauchy_alpha1_sum = Fraction(1) / (Fraction(1) - x) ** 2
        u0_norm = 80 * slope_abs * cauchy_alpha1_sum
        return {
            "cauchy_z_abs_upper": z_abs,
            "cauchy_d_abs_upper": d_abs,
            "cauchy_L_abs_lower": L_lower,
            "cauchy_U_star_abs_upper": U_abs,
            "cauchy_B_abs_upper": B_abs,
            "cauchy_H_star_abs_upper": H_abs,
            "cauchy_one_plus_z2_abs_lower": one_plus_z2_lower,
            "cauchy_Pi0_abs_upper": Pi_abs,
            "cauchy_Pi0_eta_abs_upper": Pi_eta_abs,
            "cauchy_Z_star_abs_upper": Z_abs,
            "cauchy_u0_Y_abs_upper": slope_abs,
            "rho_over_cauchy_radius": x,
            "alpha1_cauchy_weight_sum_upper": cauchy_alpha1_sum,
            "u0_coefficient_norm_upper": u0_norm,
        }

    def u0_certificate(self) -> dict[str, float]:
        return {
            name: _upper(value)
            for name, value in self._complex_u0_exact_bounds().items()
        }

    def u_ball(self) -> BallFactorBound:
        exact = self._complex_u0_exact_bounds()
        norm = exact["u0_coefficient_norm_upper"] + 1
        # Projection from the source pair ball to its u coordinate is
        # 1-Lipschitz in the pair norm used by the PA.10 contraction map.
        return BallFactorBound(norm=_upper(norm), lipschitz=1.0)

    def mixed_handoff(self) -> dict[str, dict[str, float]]:
        u = self.u_ball()
        r1 = self.mixed.r1_mixed_after_j2(u=u)
        r2 = self.mixed.r2_mixed_after_j1(u=u)
        return {
            "post_J2_d_detaAu_YPhiY": {
                "norm": float(r1.norm), "lipschitz": float(r1.lipschitz)
            },
            "post_J1_d_detaAu_YuY": {
                "norm": float(r2.norm), "lipschitz": float(r2.lipschitz)
            },
        }

    @property
    def truth_boundary(self) -> dict[str, bool]:
        return copy.deepcopy(_TRUTH_BOUNDARY)

    def report(self) -> dict[str, Any]:
        u0 = self.u0_certificate()
        u_ball = self.u_ball()
        mixed = self.mixed_handoff()
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
                "mixed_jnu_receipt_sha256": self.mixed.sha256,
                "source_phi_receipt_sha256": self.mixed.phi_source.sha256,
                "source_axis_receipt_sha256": self.domain.sha256,
            },
            "selected_source_compatible_inputs": {
                "h": self.domain.h,
                "j0": self.domain.j0,
                "pressure_square": self.domain.pressure_square,
                "sigma_star": self.domain.sigma_star,
                "rho": self.domain.coefficient_rho,
                "cauchy_radius": self.domain.cauchy_radius,
                "choice_status": (
                    "repository-autonomous source-compatible existence choice; "
                    "not recovered hidden data"
                ),
            },
            "u0_coefficient_certificate": u0,
            "source_compatible_radius_one_u_ball": {
                "norm": float(u_ball.norm),
                "lipschitz": float(u_ball.lipschitz),
                "radius_increment": 1.0,
            },
            "source_compatible_mixed_post_Jnu_handoff": mixed,
            "integration_boundary": {
                "diagnostic_unit_u_fixture_replaced_for_these_two_mixed_terms": True,
                "standalone_Y_Phi_Y_bound_required": False,
                "standalone_Y_u_Y_bound_required": False,
                "termwise_post_J_R1_R2_bridge_still_required": True,
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
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    payload = KokunoPA10SourceUBallBounds().save_report(args.output)
    print("receipt_sha256=", payload["receipt_sha256"])
    print("u_ball=", payload["source_compatible_radius_one_u_ball"])
    print("mixed_post_Jnu=", payload["source_compatible_mixed_post_Jnu_handoff"])


if __name__ == "__main__":
    _main()
