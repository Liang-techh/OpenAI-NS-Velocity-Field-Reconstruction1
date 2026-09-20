"""Map the public PA.10 contraction center into physical leading profiles.

Pinned provenance: KokunoYumeto/yang-mills-interacting-workbench at
143f6773feb424ad9ed3a8d116653200f20346b7, corrected 2026-09-09 reader.

The public reconstruction states

    Y = Lambda X,
    phi = phi_* Phi,
    U = U_* + Lambda^{-1} u,
    F = phi/C,
    phi_*(eta) = exp(Lambda * integral_0^eta zeta_*(w) dw),

with the physical leading map

    E = sqrt(2X) F,
    u_theta = q^{-A} E,
    u_z = q^{-A} U,
    r u_r = V_0 = X v_0,

and incompressibility determining V_0 from U.  This module applies those exact
mappings to the already executable PA.10 contraction center (Phi_0,u_0).  It
therefore exposes a source-native physical center profile (F_0,U_0,v_0), not
the final corrected fixed point and not yet a Cartesian spacetime velocity.

The real-axis phi_* primitive is evaluated by a deterministic partial-fraction
antiderivative of the displayed rational zeta_* formula.  This is repository
numerical realization, not hidden source data.  The configured C is inherited
unchanged from the upstream explicit autonomous configuration.  Whether that C
satisfies the source complex-domain normalization C>=sup_Omega|phi_*| remains a
separate fail-closed admission condition.
"""

from __future__ import annotations

import copy
from dataclasses import dataclass, field
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from .kokuno_pa10_phi0_profile_contract import (
    CORRECTED_RELEASE,
    CORRECTED_RELEASE_DATE,
    SOURCE_COMMIT,
    SOURCE_PATH,
    SOURCE_REPOSITORY,
    KokunoPA10Phi0ProfileContract,
)


SCHEMA = "kokuno-pa10-physical-center-profile-contract-v1"

_SOURCE_FORMULAS = {
    "similarity_radial_rescaling": "Y=Lambda X",
    "profile_rescaling": "phi=phi_* Phi; U=U_*+Lambda^-1 u",
    "phi_star": "phi_*=exp(Lambda int_0^eta zeta_*(w)dw)",
    "physical_swirl_profile": "F=phi/C=(phi_*/C)Phi",
    "physical_azimuthal_profile": "E=sqrt(2X)F",
    "physical_velocity_scaling": (
        "u_theta=q^-A E; u_z=q^-A U; r u_r=V_0=X v_0"
    ),
    "incompressibility": (
        "V0_X=L^-1(2A eta U-d U_eta+2 eta X U_X); "
        "V0=(2eta XU-2Deta M-d M_eta)/L, M=int_0^X U dx"
    ),
    "center": "Phi=Phi_0=f_0(Y chi); u=u_0=-Y Z_*/(2L)",
}

_TRUTH_BOUNDARY = {
    "public_reconstruction_source": True,
    "source_Y_equals_Lambda_X_mapping_resolved": True,
    "physical_F_from_Phi_mapping_resolved": True,
    "physical_U_from_u_mapping_resolved": True,
    "source_native_physical_center_profiles_executable": True,
    "source_native_physical_center_profiles_vectorized": True,
    "source_native_physical_center_first_derivatives_executable": True,
    "source_native_physical_center_configuration_serializable": True,
    "center_incompressibility_assembly_executable": True,
    "center_axis_regular": True,
    "real_phi_star_primitive_numerically_materialized": True,
    "source_complex_C_normalization_certified": False,
    "source_center_is_final_corrected_fixed_point": False,
    "fixed_point_correction_materialized": False,
    "global_leading_profile_reconstructed": False,
    "cartesian_spacetime_velocity_materialized": False,
    "matched_global_pressure_materialized": False,
    "complete_kokuno_composite_velocity": False,
    "heldout_ns_residual_assessed": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


def _canonical_json(payload: Mapping[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _finite_array(value: Any, name: str) -> np.ndarray:
    out = np.asarray(value, dtype=float)
    if np.any(~np.isfinite(out)):
        raise ValueError(f"{name} must contain only finite values")
    return out


@dataclass(frozen=True)
class KokunoPA10PhysicalCenterProfileContract:
    """Executable physical (F_0,U_0,v_0) map of the PA.10 center."""

    center_profiles: KokunoPA10Phi0ProfileContract = field(
        default_factory=KokunoPA10Phi0ProfileContract
    )

    def __post_init__(self) -> None:
        if not isinstance(self.center_profiles, KokunoPA10Phi0ProfileContract):
            raise TypeError("center_profiles must be KokunoPA10Phi0ProfileContract")

    @property
    def axis_profiles(self):
        return self.center_profiles.axis_profiles

    @property
    def Lambda(self) -> float:
        return float(self.axis_profiles.domain.binding.Lambda)

    @property
    def C(self) -> float:
        return float(self.axis_profiles.domain.binding.outer_schedule.C)

    @property
    def source_X_interval(self) -> tuple[float, float]:
        y0, y1 = self.center_profiles.source_Y_interval
        return (y0 / self.Lambda, y1 / self.Lambda)

    @property
    def eta_interval(self) -> tuple[float, float]:
        return self.center_profiles.eta_interval

    def _broadcast(self, X: Any, eta: Any) -> tuple[np.ndarray, np.ndarray]:
        X_arr, eta_arr = np.broadcast_arrays(
            _finite_array(X, "X"), _finite_array(eta, "eta")
        )
        x0, x1 = self.source_X_interval
        if np.any((X_arr < x0) | (X_arr > x1)):
            raise ValueError(f"X must lie in [{x0},{x1}]")
        e0, e1 = self.eta_interval
        if np.any((eta_arr < e0) | (eta_arr > e1)):
            raise ValueError("eta lies outside the certified enlarged real interval")
        return X_arr, eta_arr

    def _zeta_rational_data(self) -> tuple[np.ndarray, np.ndarray]:
        """Return poles and residues for the public real-axis zeta_* rational form."""
        domain = self.axis_profiles.domain
        h = float(domain.h)
        D = float(domain.D)
        j0 = float(domain.j0)
        sigma = float(domain.sigma_star)

        # Ascending polynomial coefficients.
        H = np.array([j0, D + 4.0, -j0, -4.0], dtype=float)
        L = np.array([1.0, 0.0, -2.0 * h], dtype=float)
        denominator = np.polynomial.polynomial.polymul(H, H)
        denominator[0] += sigma * sigma
        numerator = -np.polynomial.polynomial.polymul(L, H)
        derivative = np.arange(1, denominator.size, dtype=float) * denominator[1:]
        roots = np.polynomial.polynomial.polyroots(denominator)
        residues = np.array(
            [
                np.polynomial.polynomial.polyval(root, numerator)
                / np.polynomial.polynomial.polyval(root, derivative)
                for root in roots
            ],
            dtype=complex,
        )
        if np.any(np.abs(roots.imag) < 1.0e-14):
            raise RuntimeError("zeta_* denominator unexpectedly has a real pole")
        return roots, residues

    def zeta_primitive(self, eta: Any) -> np.ndarray:
        """Evaluate integral_0^eta zeta_*(w)dw on the certified real interval."""
        eta_arr = _finite_array(eta, "eta")
        e0, e1 = self.eta_interval
        if np.any((eta_arr < e0) | (eta_arr > e1)):
            raise ValueError("eta lies outside the certified enlarged real interval")
        roots, residues = self._zeta_rational_data()
        result = np.zeros(eta_arr.shape, dtype=complex)
        for root, residue in zip(roots, residues, strict=True):
            result += residue * np.log((eta_arr - root) / (-root))
        imag_max = float(np.max(np.abs(result.imag))) if result.size else 0.0
        if imag_max > 2.0e-10:
            raise RuntimeError("real-axis zeta_* primitive lost conjugate-pair cancellation")
        return result.real

    def phi_star(self, eta: Any) -> np.ndarray:
        exponent = self.Lambda * self.zeta_primitive(eta)
        if np.any(exponent > math.log(np.finfo(float).max)):
            raise OverflowError("configured phi_* exceeds binary64 range")
        return np.exp(exponent)

    def values(self, X: Any, eta: Any) -> dict[str, np.ndarray]:
        """Return vectorized physical center profiles on X=Y/Lambda."""
        X_arr, eta_arr = self._broadcast(X, eta)
        Y = self.Lambda * X_arr
        center = self.center_profiles.values(Y, eta_arr)
        axis = self.axis_profiles.values(Y, eta_arr)
        axis_d = self.axis_profiles.derivatives(Y, eta_arr)

        phi_star = self.phi_star(eta_arr)
        g = phi_star / self.C
        F = g * center["Phi_0"]
        U = axis["U_star"] + center["u_0"] / self.Lambda

        # At the contraction center u_0=Y*s with s=-Z_*/(2L), hence
        # U=U_*+X*s and the radial primitive is analytic.
        slope = axis["u_0_Y"]
        slope_eta = axis_d["u_0_Y_eta"]
        M_over_X = axis["U_star"] + 0.5 * X_arr * slope
        M_eta_over_X = axis["U_star_eta"] + 0.5 * X_arr * slope_eta
        numerator_over_X = (
            2.0 * eta_arr * U
            - 2.0 * float(self.axis_profiles.domain.D) * eta_arr * M_over_X
            - axis["d"] * M_eta_over_X
        )
        v0 = numerator_over_X / axis["L"]
        V0 = X_arr * v0
        E = np.sqrt(2.0 * X_arr) * F
        M = X_arr * M_over_X

        return {
            "X": X_arr,
            "Y": Y,
            "eta": eta_arr,
            "Lambda": np.full_like(X_arr, self.Lambda),
            "C": np.full_like(X_arr, self.C),
            "phi_star": phi_star,
            "g": g,
            "Phi_0": center["Phi_0"],
            "u_0": center["u_0"],
            "F_0": F,
            "E_0": E,
            "U_0": U,
            "M_0": M,
            "M_0_over_X": M_over_X,
            "v_0": v0,
            "V_0": V0,
        }

    def derivatives(self, X: Any, eta: Any) -> dict[str, np.ndarray]:
        """Return analytic first derivatives needed by later velocity assembly."""
        X_arr, eta_arr = self._broadcast(X, eta)
        Y = self.Lambda * X_arr
        center = self.center_profiles.values(Y, eta_arr)
        center_d = self.center_profiles.derivatives(Y, eta_arr)
        axis = self.axis_profiles.values(Y, eta_arr)
        axis_d = self.axis_profiles.derivatives(Y, eta_arr)
        phi_star = self.phi_star(eta_arr)
        g = phi_star / self.C
        g_eta = self.Lambda * axis["zeta_star"] * g

        slope = axis["u_0_Y"]
        slope_eta = axis_d["u_0_Y_eta"]
        v0_X = (
            (2.0 - float(self.axis_profiles.domain.D)) * eta_arr * slope
            - 0.5 * axis["d"] * slope_eta
        ) / axis["L"]

        return {
            "phi_star_eta": self.Lambda * axis["zeta_star"] * phi_star,
            "g_eta": g_eta,
            "Phi_0_X": self.Lambda * center_d["Phi_0_Y"],
            "Phi_0_eta": center_d["Phi_0_eta"],
            "F_0_X": g * self.Lambda * center_d["Phi_0_Y"],
            "F_0_eta": g_eta * center["Phi_0"] + g * center_d["Phi_0_eta"],
            "U_0_X": slope,
            "U_0_eta": axis["U_star_eta"] + center_d["u_0_eta"] / self.Lambda,
            "v_0_X": v0_X,
        }

    def incompressibility_defect(self, X: Any, eta: Any) -> np.ndarray:
        """Replay the scalar V0_X identity using the mapped center profiles."""
        values = self.values(X, eta)
        derivatives = self.derivatives(X, eta)
        X_arr = values["X"]
        eta_arr = values["eta"]
        axis = self.axis_profiles.values(values["Y"], eta_arr)
        V0_X = values["v_0"] + X_arr * derivatives["v_0_X"]
        rhs = (
            2.0 * float(self.axis_profiles.domain.A) * eta_arr * values["U_0"]
            - axis["d"] * derivatives["U_0_eta"]
            + 2.0 * eta_arr * X_arr * derivatives["U_0_X"]
        ) / axis["L"]
        return V0_X - rhs

    def configuration(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "center_profile_configuration": self.center_profiles.configuration(),
        }

    @classmethod
    def from_configuration(
        cls, payload: Mapping[str, Any]
    ) -> "KokunoPA10PhysicalCenterProfileContract":
        if not isinstance(payload, Mapping):
            raise TypeError("configuration must be a mapping")
        if payload.get("schema") != SCHEMA:
            raise ValueError("configuration schema mismatch")
        return cls(
            center_profiles=KokunoPA10Phi0ProfileContract.from_configuration(
                payload["center_profile_configuration"]
            )
        )

    def save_configuration(self, path: str | Path) -> dict[str, Any]:
        payload = self.configuration()
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        return payload

    @classmethod
    def load_configuration(
        cls, path: str | Path
    ) -> "KokunoPA10PhysicalCenterProfileContract":
        return cls.from_configuration(json.loads(Path(path).read_text()))

    @property
    def truth_boundary(self) -> dict[str, bool]:
        return copy.deepcopy(_TRUTH_BOUNDARY)

    def report(self) -> dict[str, Any]:
        x0, x1 = self.source_X_interval
        eta_probe = np.linspace(-1.0, 1.0, 8193)
        phi = self.phi_star(eta_probe)
        g = phi / self.C
        X_grid = np.linspace(x0, x1, 17)
        eta_grid = np.linspace(-1.0, 1.0, 19)
        XX, EE = np.meshgrid(X_grid, eta_grid, indexing="ij")
        values = self.values(XX, EE)
        div_defect = self.incompressibility_defect(XX, EE)
        real_phi_max = float(np.max(phi))
        real_g_max = float(np.max(np.abs(g)))
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
            "mapping_contract": {
                "source_Y_equals_Lambda_X": True,
                "Lambda": self.Lambda,
                "source_X_interval": list(self.source_X_interval),
                "source_Y_interval": list(self.center_profiles.source_Y_interval),
                "physical_F_equals_phi_star_over_C_times_Phi": True,
                "physical_U_equals_Ustar_plus_u_over_Lambda": True,
                "configured_C": self.C,
                "configured_C_source_complex_normalization_certified": False,
            },
            "configuration": self.configuration(),
            "machine_checks": {
                "Y_equals_Lambda_X_exact_on_probe": bool(
                    np.all(values["Y"] == self.Lambda * values["X"])
                ),
                "axis_V0_exact_zero_on_probe": bool(
                    np.all(self.values(np.zeros_like(eta_grid), eta_grid)["V_0"] == 0.0)
                ),
                "profiles_nontrivial": bool(
                    np.any(np.abs(values["F_0"]) > 0.0)
                    and np.any(np.abs(values["U_0"]) > 0.0)
                ),
                "all_profile_values_finite": bool(
                    all(np.all(np.isfinite(v)) for v in values.values())
                ),
                "incompressibility_scalar_identity_max_abs_defect": float(
                    np.max(np.abs(div_defect))
                ),
                "real_phi_star_probe_max": real_phi_max,
                "real_g_probe_max": real_g_max,
                "configured_C_exceeds_real_probe_phi_star": bool(self.C >= real_phi_max),
            },
            "normalization_limitation": {
                "source_condition": "C >= sup_{eta in Omega}|phi_*(eta)|",
                "complex_domain_condition_certified_here": False,
                "real_probe_is_not_complex_domain_proof": True,
            },
            "truth_boundary": self.truth_boundary,
        }
        payload["receipt_sha256"] = hashlib.sha256(
            _canonical_json(payload).encode("utf-8")
        ).hexdigest()
        return payload

    def save_report(self, path: str | Path) -> dict[str, Any]:
        payload = self.report()
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        return payload


def _main() -> None:
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--config-output", type=Path)
    args = parser.parse_args()
    profiles = KokunoPA10PhysicalCenterProfileContract()
    payload = profiles.save_report(args.output)
    if args.config_output is not None:
        profiles.save_configuration(args.config_output)
    print("receipt_sha256=", payload["receipt_sha256"])
    print("machine_checks=", payload["machine_checks"])


if __name__ == "__main__":
    _main()
