"""Selected-data analytic certificate for Kokuno's PA.10 multiplier ``chi``.

Pinned public provenance is KokunoYumeto/yang-mills-interacting-workbench at
``143f6773feb424ad9ed3a8d116653200f20346b7``, corrected reader dated
2026-09-09 (Zenodo 22678406).  The source axis construction defines

    H_* = D eta + (1-eta^2) U_*,    U_*=4 eta+j_0,
    chi = H_*^2 / (H_*^2 + sigma_*^2),

and the coefficient-space contraction later lets ``M_chi`` denote the bounded
multiplier norm of this same ``chi``.  The source does not publish one numeric
coefficient radius ``rho`` or one numeric ``M_chi``.

This module closes one narrower executable seam for the repository's *selected*
axis datum.  It proves by elementary complex-distance bounds that the selected
rational ``chi`` is holomorphic on a concrete complex tube around [-1,1], then
uses Cauchy's estimate and the source coefficient weights to obtain a finite
multiplier-norm upper bound for any explicitly selected ``rho`` below that tube
radius.  The derivation is

    |chi^(j)| <= j! M_R / R^j,

and, after Leibniz plus the exact source weight ratios,

    ||chi F||_rho <= M_R sum_{j>=0}(j+1)^2 (rho/R)^j ||F||_rho
                  = M_R (1+x)/(1-x)^3 ||F||_rho,   x=rho/R<1.

The certificate is intentionally labelled selected/autonomous.  In particular,
the repository's ``sigma_*`` has not been proved to be Kokuno's existential
compact-set choice and the source's slightly enlarged real interval ``I`` is not
numerically recovered.  Therefore this does not promote ``source_rho`` or
``source_chi_multiplier_norm`` to certified source data, and it does not open
PA.16 or claim PDE validation.
"""

from __future__ import annotations

from dataclasses import dataclass
import copy
import hashlib
import json
import math
from pathlib import Path
from typing import Any

from .kokuno_pa10_operator_primitives import KokunoPA10OperatorPrimitiveBounds
from .kokuno_rescaled_core_seed import KokunoSourceRescaledCoreSeed


SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
CORRECTED_RELEASE = "zenodo:22678406"
CORRECTED_RELEASE_DATE = "2026-09-09"
SCHEMA = "kokuno-pa10-selected-chi-multiplier-v1"

_SOURCE_FORMULAS = {
    "axis_H": "H_*=D eta+(1-eta^2)U_*; U_*=4 eta+j_0; D=1/2-h",
    "chi": "chi=H_*^2/(H_*^2+sigma_*^2)",
    "coefficient_weight": (
        "a_{alpha,beta}=20^{-alpha} rho^{-beta} beta! "
        "binom(alpha+beta,beta)/((alpha+1)^2(beta+1)^2)"
    ),
    "source_multiplier_role": "M_chi is the bounded multiplier norm of chi",
    "T_definition": "T=J_2 chi/2, multiplication by chi before integration",
}

_DERIVED_FORMULAS = {
    "tube_H_derivative_bound": (
        "sup_tube |H_*'| <= |D+4|+2|j_0|(1+R)+12(1+R)^2"
    ),
    "tube_H_variation": "|H_*(z)-H_*(eta)| <= R sup_tube |H_*'|",
    "chi_denominator_lower": (
        "|H_*^2+sigma_*^2|=|H_*-i sigma_*||H_*+i sigma_*| "
        ">= (sigma_*-delta_H)^2"
    ),
    "cauchy": "|chi^(j)| <= j! M_R/R^j",
    "multiplier_bound": (
        "M_chi <= M_R sum_{j>=0}(j+1)^2(rho/R)^j "
        "= M_R(1+rho/R)/(1-rho/R)^3"
    ),
}

_TRUTH_BOUNDARY = {
    "public_reconstruction_source": True,
    "source_chi_formula_recorded": True,
    "selected_axis_chi_complex_tube_machine_bound": True,
    "selected_coefficient_rho_executable": True,
    "selected_chi_multiplier_norm_machine_bound": True,
    "selected_chi_bound_feeds_operator_primitives": True,
    "selected_sigma_star_is_autonomous_repository_choice": True,
    "source_sigma_star_admissibility_verified": False,
    "source_enlarged_real_interval_I_recovered": False,
    "source_rho_machine_bound": False,
    "source_chi_multiplier_norm_machine_bound": False,
    "source_R1_radius_one_ball_norm_machine_bound": False,
    "source_R2_radius_one_ball_norm_machine_bound": False,
    "source_R1_radius_one_ball_lipschitz_machine_bound": False,
    "source_R2_radius_one_ball_lipschitz_machine_bound": False,
    "source_operator_constant_M_machine_bound": False,
    "source_operator_constant_K_machine_bound": False,
    "source_fixed_point_distance_machine_bound": False,
    "source_B0_dependencies_machine_bound": False,
    "source_B0_analytic_bound_proved": False,
    "source_T_sh_lower_bound_verified": False,
    "selected_pa16_handoff_allowed": False,
    "inner_to_outer_join_completed": False,
    "global_pressure_matched": False,
    "global_leading_profile_reconstructed": False,
    "formal_full_domain_pde_gate_assessed": False,
    "heldout_ns_residual_assessed": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


def _canonical_json(payload: dict[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _positive_finite(value: float, name: str) -> float:
    out = float(value)
    if not math.isfinite(out) or out <= 0.0:
        raise ValueError(f"{name} must be finite and >0")
    return out


@dataclass(frozen=True)
class KokunoPA10SelectedChiMultiplier:
    """Complex-tube and multiplier certificate for the selected core datum.

    ``analytic_tube_radius`` is a repository-selected complex distance from the
    closed real interval [-1,1]. ``coefficient_rho_fraction`` sets
    ``rho=fraction*R``.  Both remain autonomous numerical choices; the source
    only requires sufficiently small positive radii.
    """

    analytic_tube_radius: float = 5.0e-3
    coefficient_rho_fraction: float = 0.1

    def __post_init__(self) -> None:
        R = _positive_finite(self.analytic_tube_radius, "analytic_tube_radius")
        fraction = _positive_finite(
            self.coefficient_rho_fraction, "coefficient_rho_fraction"
        )
        if R >= 0.25:
            raise ValueError("analytic_tube_radius must be <0.25")
        if fraction >= 1.0:
            raise ValueError("coefficient_rho_fraction must be <1")
        object.__setattr__(self, "analytic_tube_radius", R)
        object.__setattr__(self, "coefficient_rho_fraction", fraction)
        # Fail closed at construction time if the elementary tube proof does
        # not actually certify all explicit selected axis denominators.
        certificate = self._certificate_unchecked()
        if not certificate["selected_axis_complex_domain_certified"]:
            raise ValueError("selected complex tube is not certified by the explicit bounds")

    @property
    def core(self) -> KokunoSourceRescaledCoreSeed:
        return KokunoSourceRescaledCoreSeed()

    @property
    def h(self) -> float:
        return float(self.core.h)

    @property
    def D(self) -> float:
        return float(self.core.D)

    @property
    def j0(self) -> float:
        return float(self.core.j0)

    @property
    def sigma_star(self) -> float:
        return float(self.core.sigma_star)

    @property
    def coefficient_rho(self) -> float:
        return float(self.analytic_tube_radius * self.coefficient_rho_fraction)

    def H_star_complex(self, z: complex) -> complex:
        zz = complex(z)
        return self.D * zz + (1.0 - zz * zz) * (4.0 * zz + self.j0)

    def chi_complex(self, z: complex) -> complex:
        H = self.H_star_complex(z)
        denominator = H * H + self.sigma_star * self.sigma_star
        if denominator == 0.0:
            raise ZeroDivisionError("selected chi denominator vanished")
        return H * H / denominator

    def _certificate_unchecked(self) -> dict[str, float | bool]:
        R = self.analytic_tube_radius
        rho = self.coefficient_rho
        z_abs_upper = 1.0 + R

        # H'(z)=D+4-2*j0*z-12*z^2.  Any point on a straight segment
        # joining eta in [-1,1] to its tube point has |z|<=1+R.
        H_prime_abs_upper = (
            abs(self.D + 4.0)
            + 2.0 * abs(self.j0) * z_abs_upper
            + 12.0 * z_abs_upper * z_abs_upper
        )
        H_variation_upper = R * H_prime_abs_upper

        # On the real interval, |D eta|<=|D| and
        # |(1-eta^2)(4eta+j0)|<=4+|j0|.
        H_real_abs_upper = abs(self.D) + 4.0 + abs(self.j0)
        H_tube_abs_upper = H_real_abs_upper + H_variation_upper

        factor_lower = self.sigma_star - H_variation_upper
        chi_denominator_abs_lower = factor_lower * factor_lower
        chi_sup_abs_upper = (
            H_tube_abs_upper * H_tube_abs_upper / chi_denominator_abs_lower
            if factor_lower > 0.0
            else math.inf
        )

        # Other explicit axis denominators used by the same selected core.
        # L=1-2h eta^2 and Pi0 uses (1+eta^2)^-2.
        L_abs_lower = 1.0 - 2.0 * self.h * z_abs_upper * z_abs_upper
        one_plus_eta_squared_abs_lower = (1.0 - R) * (1.0 - R)

        x = rho / R
        multiplier_series_factor = (1.0 + x) / ((1.0 - x) ** 3)
        multiplier_norm_upper = chi_sup_abs_upper * multiplier_series_factor

        certified = bool(
            R < 1.0
            and rho < R
            and factor_lower > 0.0
            and L_abs_lower > 0.0
            and one_plus_eta_squared_abs_lower > 0.0
            and math.isfinite(multiplier_norm_upper)
        )
        return {
            "analytic_tube_radius": float(R),
            "coefficient_rho": float(rho),
            "rho_over_tube_radius": float(x),
            "z_abs_upper": float(z_abs_upper),
            "H_real_abs_upper": float(H_real_abs_upper),
            "H_prime_abs_upper_on_tube": float(H_prime_abs_upper),
            "H_variation_abs_upper": float(H_variation_upper),
            "H_plusminus_i_sigma_factor_abs_lower": float(factor_lower),
            "chi_denominator_abs_lower": float(chi_denominator_abs_lower),
            "H_tube_abs_upper": float(H_tube_abs_upper),
            "chi_sup_abs_upper": float(chi_sup_abs_upper),
            "L_abs_lower": float(L_abs_lower),
            "one_plus_eta_squared_abs_lower": float(one_plus_eta_squared_abs_lower),
            "coefficient_multiplier_series_factor": float(multiplier_series_factor),
            "selected_multiplier_norm_chi_upper": float(multiplier_norm_upper),
            "selected_axis_complex_domain_certified": certified,
        }

    def certificate(self) -> dict[str, float | bool]:
        return copy.deepcopy(self._certificate_unchecked())

    @property
    def multiplier_norm_chi_upper(self) -> float:
        return float(self._certificate_unchecked()["selected_multiplier_norm_chi_upper"])

    def operator_input_receipt(self) -> dict[str, float]:
        """Feed the selected rho/M_chi bound into #569 universal arithmetic.

        This is an executable selected-data receipt only.  It deliberately does
        not mutate #569's source truth flags.
        """
        operators = KokunoPA10OperatorPrimitiveBounds()
        mchi = self.multiplier_norm_chi_upper
        return {
            "coefficient_rho": self.coefficient_rho,
            "selected_multiplier_norm_chi_upper": mchi,
            "mixed_derivative_template_constant": operators.mixed_derivative_template_constant(
                self.coefficient_rho
            ),
            "inverse_one_plus_T_absolute_bound": operators.inverse_one_plus_t_bound(mchi),
        }

    @property
    def truth_boundary(self) -> dict[str, bool]:
        return copy.deepcopy(_TRUTH_BOUNDARY)

    def report(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "source": {
                "repository": SOURCE_REPOSITORY,
                "commit": SOURCE_COMMIT,
                "path": SOURCE_PATH,
                "corrected_release": CORRECTED_RELEASE,
                "corrected_release_date": CORRECTED_RELEASE_DATE,
            },
            "source_formulas": copy.deepcopy(_SOURCE_FORMULAS),
            "derived_formulas": copy.deepcopy(_DERIVED_FORMULAS),
            "selected_inputs": {
                "eta_interval": [-1.0, 1.0],
                "analytic_tube_radius": self.analytic_tube_radius,
                "coefficient_rho_fraction": self.coefficient_rho_fraction,
                "h": self.h,
                "D": self.D,
                "j0": self.j0,
                "sigma_star": self.sigma_star,
                "sigma_star_status": (
                    "selected repository B.13 value; source compact-set admissibility not certified"
                ),
            },
            "certificate": self.certificate(),
            "operator_input_receipt": self.operator_input_receipt(),
            "interpretation": (
                "rigorous elementary bound for the repository-selected explicit chi on [-1,1]; "
                "not recovery of Kokuno's hidden/enlarged-I rho or existential sigma_* choice"
            ),
            "truth_boundary": self.truth_boundary,
        }

    def to_payload(self) -> dict[str, Any]:
        payload = self.report()
        payload["sha256"] = hashlib.sha256(_canonical_json(payload).encode()).hexdigest()
        return payload

    @property
    def sha256(self) -> str:
        return str(self.to_payload()["sha256"])

    def save_json(self, path: str | Path) -> None:
        Path(path).write_text(json.dumps(self.to_payload(), indent=2, sort_keys=True) + "\n")

    @classmethod
    def from_payload(cls, payload: dict[str, Any]) -> "KokunoPA10SelectedChiMultiplier":
        raw = copy.deepcopy(payload)
        supplied_sha = raw.pop("sha256", None)
        expected_sha = hashlib.sha256(_canonical_json(raw).encode()).hexdigest()
        if supplied_sha != expected_sha:
            raise ValueError("payload sha256 mismatch")
        if raw.get("schema") != SCHEMA:
            raise ValueError("schema mismatch")
        selected = raw.get("selected_inputs")
        if not isinstance(selected, dict):
            raise ValueError("selected input metadata mismatch")
        candidate = cls(
            analytic_tube_radius=float(selected.get("analytic_tube_radius")),
            coefficient_rho_fraction=float(selected.get("coefficient_rho_fraction")),
        )
        baseline = candidate.report()
        for key in (
            "source",
            "source_formulas",
            "derived_formulas",
            "selected_inputs",
            "certificate",
            "operator_input_receipt",
            "interpretation",
            "truth_boundary",
        ):
            if raw.get(key) != baseline[key]:
                raise ValueError(f"{key.replace('_', ' ')} mismatch")
        return candidate

    @classmethod
    def load_json(cls, path: str | Path) -> "KokunoPA10SelectedChiMultiplier":
        return cls.from_payload(json.loads(Path(path).read_text()))
