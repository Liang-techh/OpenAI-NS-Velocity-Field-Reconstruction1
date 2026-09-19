"""Executable source-compatible PA.10 axis-domain certificate.

Public provenance: KokunoYumeto/yang-mills-interacting-workbench at
143f6773feb424ad9ed3a8d116653200f20346b7, corrected 2026-09-09 reader
(Zenodo 22678406).  The source chooses delta_* and then sufficiently small
sigma_* so that chi>0.99 on K_delta={|Z_*|<=delta_*}, enlarges [-1,1]
slightly to I, chooses a complex neighborhood Omega where L and
H_*^2+sigma_*^2 do not vanish and Pi_0 is holomorphic, then normalizes
C>=sup_Omega|phi_*| so g=phi_*/C has complex supremum <=1.

The repository's older B.13 selected center used autonomous sigma_*=0.5 and
never certified it as the source PA.8 choice.  This module makes a new explicit
source-compatible choice sigma_*=1e-3, certifies PA.8 by an exact-rational
finite cover, certifies a concrete complex tube/Cauchy radius, and exposes a
source-compatible coefficient rho and normalized-g coefficient norm.  These
are admissible reconstruction choices, not recovered hidden OpenAI/Kokuno
numbers.  Phi ball bounds, R1/R2, M/K, global pressure/velocity and PDE
validation remain open.
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

from .kokuno_pressure_datum_binding import KokunoPressureDatumBinding

SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
CORRECTED_RELEASE = "zenodo:22678406"
CORRECTED_RELEASE_DATE = "2026-09-09"
SCHEMA = "kokuno-pa10-source-compatible-axis-domain-v1"

_SOURCE_FORMULAS = {
    "axis": (
        "U_*=4eta+j_0; H_*=Deta+(1-eta^2)U_*; "
        "Z_*=-A(1-2etaU_*)U_*-H_*U_*'-dPi_0'+4AetaPi_0"
    ),
    "PA8": (
        "K_delta={|Z_*|<=delta_*}; if |H_*|>=m on K_delta and "
        "sigma_*^2<m^2/99 then chi=H_*^2/(H_*^2+sigma_*^2)>0.99"
    ),
    "complex_domain": (
        "strict margins allow I superset [-1,1] and bounded simply connected "
        "Omega with L,H_*^2+sigma_*^2 nonzero and Pi_0 holomorphic"
    ),
    "g_normalization": "C>=sup_Omega|phi_*|; g=phi_*/C; sup_Omega|g|<=1",
    "coefficient_weight": (
        "a_{alpha,beta}=20^{-alpha}rho^{-beta}beta!binom(alpha+beta,beta)/"
        "((alpha+1)^2(beta+1)^2)"
    ),
}

_TRUTH_BOUNDARY = {
    "public_reconstruction_source": True,
    "selected_source_compatible_sigma_PA8_machine_certified": True,
    "selected_enlarged_real_interval_machine_certified": True,
    "selected_complex_neighborhood_machine_certified": True,
    "selected_source_compatible_coefficient_rho_machine_bound": True,
    "source_normalized_g_coefficient_norm_upper_executable": True,
    "source_normalized_g_requires_later_C_choice": True,
    "old_selected_sigma_star_equals_source_compatible_choice": False,
    "old_selected_local_pressure_audit_transfers_to_new_sigma": False,
    "source_hidden_sigma_star_recovered": False,
    "source_unique_rho_recovered": False,
    "source_hidden_C_recovered": False,
    "source_Phi_radius_one_ball_norm_machine_bound": False,
    "source_Phi_radius_one_ball_lipschitz_machine_bound": False,
    "source_mixed_Y_Phi_Y_radius_one_ball_bound_machine_bound": False,
    "source_mixed_Y_u_Y_radius_one_ball_bound_machine_bound": False,
    "source_pressure_radius_one_ball_norm_machine_bound": False,
    "source_pressure_radius_one_ball_lipschitz_machine_bound": False,
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
    "formal_full_domain_pde_gate_assessed": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


def _canonical_json(payload: dict[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _q(value: float) -> Fraction:
    out = float(value)
    if not math.isfinite(out):
        raise ValueError("bound input must be finite")
    return Fraction.from_float(out)


def _positive(value: float, name: str) -> float:
    out = float(value)
    if not math.isfinite(out) or out <= 0.0:
        raise ValueError(f"{name} must be finite and positive")
    return out


def _upper(value: Fraction) -> float:
    out = float(value)
    if not math.isfinite(out):
        raise OverflowError("upper bound is outside binary64 range")
    if Fraction.from_float(out) < value:
        out = math.nextafter(out, math.inf)
    return out


def _lower(value: Fraction) -> float:
    out = float(value)
    if not math.isfinite(out):
        raise OverflowError("lower bound is outside binary64 range")
    if Fraction.from_float(out) > value:
        out = math.nextafter(out, -math.inf)
    return out


@dataclass(frozen=True)
class KokunoPA10SourceCompatibleAxisDomain:
    """Concrete PA.8 / complex-neighborhood realization."""

    binding: KokunoPressureDatumBinding = field(default_factory=KokunoPressureDatumBinding)
    sigma_star: float = 1.0e-3
    delta_star: float = 0.5
    enlarged_real_margin: float = 1.0e-4
    partition_intervals: int = 8192
    complex_tube_radius: float = 2.0e-5
    cauchy_radius: float = 1.0e-5
    coefficient_rho: float = 1.0e-6

    def __post_init__(self) -> None:
        if not isinstance(self.binding, KokunoPressureDatumBinding):
            raise TypeError("binding must be a KokunoPressureDatumBinding")
        for name in (
            "sigma_star", "delta_star", "enlarged_real_margin",
            "complex_tube_radius", "cauchy_radius", "coefficient_rho",
        ):
            object.__setattr__(self, name, _positive(getattr(self, name), name))
        if self.sigma_star >= 0.1:
            raise ValueError("sigma_star must be <0.1")
        if self.enlarged_real_margin >= 0.05:
            raise ValueError("enlarged_real_margin must be <0.05")
        if not self.cauchy_radius < self.complex_tube_radius:
            raise ValueError("cauchy_radius must be < complex_tube_radius")
        if not self.coefficient_rho < self.cauchy_radius:
            raise ValueError("coefficient_rho must be < cauchy_radius")
        if isinstance(self.partition_intervals, bool) or not isinstance(
            self.partition_intervals, (int, np.integer)
        ):
            raise TypeError("partition_intervals must be an integer")
        count = int(self.partition_intervals)
        if not 1024 <= count <= 65536:
            raise ValueError("partition_intervals must lie in [1024,65536]")
        object.__setattr__(self, "partition_intervals", count)
        if not self._real_certificate()["source_PA8_implication_machine_certified"]:
            raise ValueError("configured PA.8 finite cover failed")
        if not self._complex_certificate()["source_complex_neighborhood_machine_certified"]:
            raise ValueError("configured complex-domain certificate failed")

    @property
    def h(self) -> float:
        return float(self.binding.outer_schedule.h)

    @property
    def A(self) -> float:
        return 0.5 + self.h

    @property
    def D(self) -> float:
        return 0.5 - self.h

    @property
    def j0(self) -> float:
        return float(self.binding.j0)

    @property
    def pressure_square(self) -> float:
        value = float(self.binding.selected_pressure_scale**2)
        if not math.isfinite(value) or value <= 0.0:
            raise OverflowError("selected pressure square is invalid")
        return value

    @property
    def old_selected_sigma_star(self) -> float:
        return float(self.binding.sigma)

    def axis_state(self, eta: Any) -> dict[str, np.ndarray]:
        """Vectorized new-sigma axis datum on the certified enlarged interval."""
        e = np.asarray(eta, dtype=float)
        if np.any(~np.isfinite(e)):
            raise ValueError("eta must be finite")
        if np.any(np.abs(e) > 1.0 + self.enlarged_real_margin):
            raise ValueError("eta lies outside the certified enlarged interval")
        d = 1.0 - e * e
        L = 1.0 - 2.0 * self.h * e * e
        U = 4.0 * e + self.j0
        H = self.D * e + d * U
        B = 1.0 - 2.0 * e * U
        den = 1.0 + e * e
        Pi = -self.pressure_square / den**2
        Pi_eta = 4.0 * self.pressure_square * e / den**3
        Z = -self.A * B * U - 4.0 * H - d * Pi_eta + 4.0 * self.A * e * Pi
        sigma2 = self.sigma_star**2
        return {
            "d": d, "L": L, "U_star": U, "H_star": H, "Z_star": Z,
            "Pi_0": Pi, "Pi_0_eta": Pi_eta,
            "chi": H * H / (H * H + sigma2),
            "zeta_star": -L * H / (H * H + sigma2),
        }

    def _constants(self) -> dict[str, Fraction]:
        h = _q(self.h)
        return {
            "one": Fraction(1), "h": h, "A": Fraction(1, 2) + h,
            "D": Fraction(1, 2) - h, "j0": _q(self.j0),
            "p2": _q(self.pressure_square),
        }

    def _H(self, eta: Fraction) -> Fraction:
        c = self._constants()
        return c["D"] * eta + (c["one"] - eta * eta) * (4 * eta + c["j0"])

    def _Z(self, eta: Fraction) -> Fraction:
        c = self._constants()
        d = c["one"] - eta * eta
        U = 4 * eta + c["j0"]
        H = c["D"] * eta + d * U
        B = c["one"] - 2 * eta * U
        den = c["one"] + eta * eta
        Pi = -c["p2"] / den**2
        Pi_eta = 4 * c["p2"] * eta / den**3
        return -c["A"] * B * U - 4 * H - d * Pi_eta + 4 * c["A"] * eta * Pi

    def _real_derivative_bounds(self) -> tuple[Fraction, Fraction]:
        c = self._constants()
        E = c["one"] + _q(self.enlarged_real_margin)
        H_prime = abs(c["D"] + 4) + 2 * abs(c["j0"]) * E + 12 * E * E
        U_abs = 4 * E + abs(c["j0"])
        B_abs = c["one"] + 2 * E * U_abs
        B_eta_abs = 2 * U_abs + 8 * E
        d_abs = c["one"] + E * E
        Pi_eta_abs = 4 * c["p2"] * E
        Pi_etaeta_abs = 4 * c["p2"] * (c["one"] + 5 * E * E)
        Z_prime = (
            c["A"] * (B_eta_abs * U_abs + 4 * B_abs) + 4 * H_prime
            + 2 * E * Pi_eta_abs + d_abs * Pi_etaeta_abs
            + 4 * c["A"] * (c["p2"] + E * Pi_eta_abs)
        )
        return H_prime, Z_prime

    def _real_certificate(self) -> dict[str, Any]:
        one = Fraction(1)
        margin, sigma, delta = map(_q, (
            self.enlarged_real_margin, self.sigma_star, self.delta_star
        ))
        lower, upper = -one - margin, one + margin
        step = (upper - lower) / self.partition_intervals
        half = step / 2
        H_prime, Z_prime = self._real_derivative_bounds()
        required_H = 10 * sigma
        h_count = z_count = failed = 0
        min_H: Fraction | None = None
        min_Z: Fraction | None = None
        for index in range(self.partition_intervals):
            mid = lower + Fraction(2 * index + 1, 2) * step
            H_lower = abs(self._H(mid)) - H_prime * half
            if H_lower > required_H:
                h_count += 1
                min_H = H_lower if min_H is None or H_lower < min_H else min_H
                continue
            Z_lower = self._Z(mid) - Z_prime * half
            if Z_lower > 2 * delta:
                z_count += 1
                min_Z = Z_lower if min_Z is None or Z_lower < min_Z else min_Z
            else:
                failed += 1
        chi_lower = required_H**2 / (required_H**2 + sigma**2)
        certified = bool(
            failed == 0 and h_count > 0 and z_count > 0
            and 99 * sigma**2 < required_H**2
            and chi_lower > Fraction(99, 100)
        )
        return {
            "enlarged_real_interval": [_lower(lower), _upper(upper)],
            "partition_intervals": self.partition_intervals,
            "partition_step_upper": _upper(step),
            "global_H_prime_abs_upper": _upper(H_prime),
            "global_Z_prime_abs_upper": _upper(Z_prime),
            "delta_star": self.delta_star,
            "sigma_star": self.sigma_star,
            "required_H_abs_lower_on_K": _upper(required_H),
            "PA8_chi_lower_from_dichotomy": _lower(chi_lower),
            "H_large_intervals": h_count,
            "Z_positive_exclusion_intervals": z_count,
            "failed_intervals": failed,
            "minimum_certified_H_abs_lower_on_H_branch": _lower(min_H),
            "minimum_certified_Z_lower_on_H_small_branch": _lower(min_Z),
            "source_PA8_implication_machine_certified": certified,
        }

    def real_interval_certificate(self) -> dict[str, Any]:
        return copy.deepcopy(self._real_certificate())

    def _complex_certificate(self) -> dict[str, Any]:
        c = self._constants()
        one = c["one"]
        margin, tube, cauchy, rho, sigma = map(_q, (
            self.enlarged_real_margin, self.complex_tube_radius,
            self.cauchy_radius, self.coefficient_rho, self.sigma_star,
        ))
        z_abs = one + margin + tube
        H_prime = abs(c["D"] + 4) + 2 * abs(c["j0"]) * z_abs + 12 * z_abs**2
        H_variation = tube * H_prime
        factor_lower = sigma - H_variation
        denominator_lower = factor_lower**2
        L_lower = one - 2 * c["h"] * z_abs**2
        Pi_factor_lower = (one - tube)**2
        x = rho / cauchy
        g_norm = (one + x) / (one - x)**3
        certified = bool(
            tube < one and cauchy < tube and rho < cauchy
            and factor_lower > 0 and L_lower > 0 and Pi_factor_lower > 0
        )
        return {
            "complex_tube_radius": self.complex_tube_radius,
            "cauchy_radius": self.cauchy_radius,
            "coefficient_rho": self.coefficient_rho,
            "rho_over_cauchy_radius": _upper(x),
            "tube_z_abs_upper": _upper(z_abs),
            "tube_H_prime_abs_upper": _upper(H_prime),
            "tube_H_variation_abs_upper": _upper(H_variation),
            "H_plusminus_i_sigma_abs_lower": _lower(factor_lower),
            "H_square_plus_sigma_square_abs_lower": _lower(denominator_lower),
            "L_abs_lower": _lower(L_lower),
            "Pi_0_one_plus_z_squared_abs_lower": _lower(Pi_factor_lower),
            "source_normalized_g_complex_sup_upper": 1.0,
            "source_normalized_g_coefficient_norm_upper": _upper(g_norm),
            "source_C_condition": "C >= sup_{eta in Omega}|phi_*(eta)|",
            "source_C_numeric_value_materialized": False,
            "source_complex_neighborhood_machine_certified": certified,
        }

    def complex_domain_certificate(self) -> dict[str, Any]:
        return copy.deepcopy(self._complex_certificate())

    def pressure_operator_inputs(self) -> dict[str, Any]:
        cert = self._complex_certificate()
        return {
            "rho": self.coefficient_rho,
            "g_norm_upper": cert["source_normalized_g_coefficient_norm_upper"],
            "g_normalization_contract": cert["source_C_condition"],
            "Phi_radius_one_ball_norm": None,
            "Phi_radius_one_ball_lipschitz": None,
            "ready_for_pressure_operator_after_Phi_bound": True,
            "ready_for_source_pressure_ball_now": False,
        }

    @property
    def truth_boundary(self) -> dict[str, bool]:
        return copy.deepcopy(_TRUTH_BOUNDARY)

    def report(self) -> dict[str, Any]:
        payload: dict[str, Any] = {
            "schema": SCHEMA,
            "source": {
                "repository": SOURCE_REPOSITORY, "commit": SOURCE_COMMIT,
                "path": SOURCE_PATH, "corrected_release": CORRECTED_RELEASE,
                "corrected_release_date": CORRECTED_RELEASE_DATE,
            },
            "source_formulas": copy.deepcopy(_SOURCE_FORMULAS),
            "selected_source_compatible_inputs": {
                "h": self.h, "j0": self.j0, "pressure_square": self.pressure_square,
                "sigma_star": self.sigma_star, "delta_star": self.delta_star,
                "old_repository_selected_sigma_star": self.old_selected_sigma_star,
                "old_sigma_reused": False,
                "enlarged_real_margin": self.enlarged_real_margin,
                "complex_tube_radius": self.complex_tube_radius,
                "cauchy_radius": self.cauchy_radius,
                "coefficient_rho": self.coefficient_rho,
                "choice_status": (
                    "repository-autonomous source-compatible existence choice; "
                    "not recovered hidden data"
                ),
            },
            "real_interval_certificate": self.real_interval_certificate(),
            "complex_domain_certificate": self.complex_domain_certificate(),
            "pressure_operator_inputs": self.pressure_operator_inputs(),
            "integration_boundary": {
                "old_sigma_0p5_local_pressure_interface_reused": False,
                "new_sigma_requires_new_center_pressure_preflight": True,
                "Phi_ball_bound_still_required": True,
                "mixed_Y_Phi_Y_bound_still_required": True,
                "mixed_Y_u_Y_bound_still_required": True,
                "R1_R2_M_K_not_promoted": True,
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
    payload = KokunoPA10SourceCompatibleAxisDomain().save_report(args.output)
    real = payload["real_interval_certificate"]
    domain = payload["complex_domain_certificate"]
    print("receipt_sha256=", payload["receipt_sha256"])
    print("PA8=", real["source_PA8_implication_machine_certified"])
    print("H/Z intervals=", real["H_large_intervals"], real["Z_positive_exclusion_intervals"])
    print("rho/g_norm=", domain["coefficient_rho"], domain["source_normalized_g_coefficient_norm_upper"])


if __name__ == "__main__":
    _main()
