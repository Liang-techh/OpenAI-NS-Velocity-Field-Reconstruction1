"""Source-native PA.10 Phi_0 contraction-center profile contract.

Pinned public provenance is KokunoYumeto/yang-mills-interacting-workbench at
``143f6773feb424ad9ed3a8d116653200f20346b7``,
``navier-stokes/navier_stokes_workbench.tex`` (corrected reader dated
2026-09-09, Zenodo 22678406).

The corrected public reconstruction gives the PA.10 contraction center

    Phi_0(Y,eta) = f_0(Y chi(eta)),
    f_0(z) = sum_{alpha>=0} (-z/2)^alpha/(alpha! (alpha+1)!),
    u_0(Y,eta) = -Y Z_*(eta)/(2 L(eta)),

with ``2(z f_0'' + 2 f_0') = -f_0`` and ``f_0(0)=1``.  This module
makes the Phi_0 half of that pair directly callable on the source-native
``(Y,eta)`` rectangle, while reusing the already executable Agent-1 axis/base
contract for ``chi`` and ``u_0``.

The finite polynomial used to evaluate the entire function ``f_0`` is a
repository numerical realization of the displayed infinite series.  The
source-native radial coordinate is still ``Y``: this module does not identify
repository ``X`` with source ``Y``, does not materialize the fixed-point
correction ``(Phi-Phi_0,u-u_0)``, and does not claim a Cartesian/global
velocity or PDE validation.
"""

from __future__ import annotations

import copy
from dataclasses import dataclass, field
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from .kokuno_pa10_leading_axis_profile_contract import (
    CORRECTED_RELEASE,
    CORRECTED_RELEASE_DATE,
    SOURCE_COMMIT,
    SOURCE_PATH,
    SOURCE_REPOSITORY,
    KokunoPA10LeadingAxisProfileContract,
)


SCHEMA = "kokuno-pa10-phi0-profile-contract-v1"
DEFAULT_SERIES_TERMS = 64

_SOURCE_FORMULAS = {
    "Phi_0": "Phi_0(Y,eta)=f_0(Y chi(eta))",
    "f_0": (
        "f_0(z)=sum_{alpha>=0}(-z/2)^alpha/"
        "(alpha!(alpha+1)!)"
    ),
    "f_0_ode": "2(z f_0''+2 f_0')=-f_0; f_0(0)=1",
    "u_0": "u_0(Y,eta)=-Y Z_*(eta)/(2L(eta))",
    "fixed_point_pair": (
        "(Phi,u)=(Phi_0,u_0)+Lambda^-1((1+T)^-1 J_2 R_1/2,"
        "J_1 R_2/2)"
    ),
}

_TRUTH_BOUNDARY = {
    "public_reconstruction_source": True,
    "source_native_phi0_contraction_center_executable": True,
    "source_native_phi0_vectorized": True,
    "source_native_phi0_analytic_derivatives_executable": True,
    "source_native_phi0_configuration_serializable": True,
    "source_phi0_axis_regular": True,
    "source_phi0_is_contraction_center_not_final_corrected_Phi": True,
    "source_u0_is_contraction_center_not_final_corrected_u": True,
    "finite_f0_series_is_repository_numerical_realization": True,
    "repository_X_identified_with_source_Y": False,
    "physical_F_from_Phi_mapping_resolved": False,
    "fixed_point_correction_materialized": False,
    "global_leading_profile_reconstructed": False,
    "complete_kokuno_composite_velocity": False,
    "unified_cartesian_velocity_export_ready": False,
    "heldout_ns_residual_assessed": False,
    "formal_full_domain_pde_gate_assessed": False,
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
class KokunoPA10Phi0ProfileContract:
    """Executable source-native PA.10 contraction-center pair ``(Phi_0,u_0)``."""

    axis_profiles: KokunoPA10LeadingAxisProfileContract = field(
        default_factory=KokunoPA10LeadingAxisProfileContract
    )
    series_terms: int = DEFAULT_SERIES_TERMS

    def __post_init__(self) -> None:
        if not isinstance(self.axis_profiles, KokunoPA10LeadingAxisProfileContract):
            raise TypeError("axis_profiles must be KokunoPA10LeadingAxisProfileContract")
        if isinstance(self.series_terms, bool) or int(self.series_terms) != self.series_terms:
            raise TypeError("series_terms must be an integer")
        if int(self.series_terms) < 16:
            raise ValueError("series_terms must be at least 16")

    @property
    def source_Y_interval(self) -> tuple[float, float]:
        return self.axis_profiles.source_Y_interval

    @property
    def eta_interval(self) -> tuple[float, float]:
        return self.axis_profiles.eta_interval

    def _coefficients(self) -> np.ndarray:
        """Return the first ``series_terms`` source-series coefficients."""
        coeff = np.empty(int(self.series_terms), dtype=float)
        coeff[0] = 1.0
        for alpha in range(int(self.series_terms) - 1):
            coeff[alpha + 1] = (
                coeff[alpha] * (-0.5) / ((alpha + 1) * (alpha + 2))
            )
        return coeff

    def f0_triplet(self, z: Any) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Evaluate ``f_0``, ``f_0'`` and ``f_0''`` by one frozen series.

        Derivative coefficient arrays are derived from the same finite
        polynomial, so the public derivative path is exactly the derivative of
        the numerical value path rather than a separately fitted approximation.
        """
        z_arr = _finite_array(z, "z")
        if np.any(z_arr < 0.0):
            raise ValueError("source PA.10 f_0 evaluator expects z=Y*chi >= 0")
        coeff = self._coefficients()
        first = np.polynomial.polynomial.polyval(z_arr, coeff)
        first_coeff = np.arange(1, coeff.size, dtype=float) * coeff[1:]
        second_coeff = (
            np.arange(2, coeff.size, dtype=float)
            * np.arange(1, coeff.size - 1, dtype=float)
            * coeff[2:]
        )
        derivative = np.polynomial.polynomial.polyval(z_arr, first_coeff)
        second = np.polynomial.polynomial.polyval(z_arr, second_coeff)
        return first, derivative, second

    def _tail_upper(self, z: np.ndarray) -> np.ndarray:
        """Bound the omitted absolute f_0 tail on the real source interval."""
        coeff = self._coefficients()
        n = int(self.series_terms)
        first_omitted_coeff = coeff[-1] * (-0.5) / (n * (n + 1))
        first_omitted = abs(first_omitted_coeff) * np.power(z, n)
        ratio = (0.5 * z) / ((n + 1) * (n + 2))
        if np.any(ratio >= 1.0):
            raise RuntimeError("series tail ratio is not contractive")
        return first_omitted / (1.0 - ratio)

    def values(self, Y: Any, eta: Any) -> dict[str, np.ndarray]:
        """Return vectorized source-native contraction-center profile values."""
        axis = self.axis_profiles.values(Y, eta)
        chi = axis["chi"]
        z = axis["Y"] * chi
        Phi0, f0_prime, f0_second = self.f0_triplet(z)
        return {
            "Y": axis["Y"],
            "eta": axis["eta"],
            "chi": chi,
            "f0_argument": z,
            "Phi_0": Phi0,
            "u_0": axis["u_0"],
            "f0_prime": f0_prime,
            "f0_second": f0_second,
            "f0_tail_upper": self._tail_upper(z),
        }

    def derivatives(self, Y: Any, eta: Any) -> dict[str, np.ndarray]:
        """Return analytic first/mixed derivatives of the numerical profile path."""
        axis = self.axis_profiles.values(Y, eta)
        axis_derivatives = self.axis_profiles.derivatives(Y, eta)
        chi = axis["chi"]
        chi_eta = axis_derivatives["chi_eta"]
        z = axis["Y"] * chi
        _, f0_prime, f0_second = self.f0_triplet(z)
        z_eta = axis["Y"] * chi_eta
        return {
            "Phi_0_Y": f0_prime * chi,
            "Phi_0_eta": f0_prime * z_eta,
            "Phi_0_YY": f0_second * chi * chi,
            "Phi_0_Y_eta": chi_eta * (f0_prime + z * f0_second),
            "u_0_Y": axis_derivatives["u_0_Y"],
            "u_0_eta": axis_derivatives["u_0_eta"],
            "u_0_Y_eta": axis_derivatives["u_0_Y_eta"],
        }

    def configuration(self) -> dict[str, Any]:
        """Return a canonical reloadable configuration for this profile contract."""
        return {
            "schema": SCHEMA,
            "series_terms": int(self.series_terms),
            "axis_profile_configuration": self.axis_profiles.configuration(),
        }

    @classmethod
    def from_configuration(
        cls, payload: Mapping[str, Any]
    ) -> "KokunoPA10Phi0ProfileContract":
        if not isinstance(payload, Mapping):
            raise TypeError("configuration must be a mapping")
        if payload.get("schema") != SCHEMA:
            raise ValueError("configuration schema mismatch")
        axis_profiles = KokunoPA10LeadingAxisProfileContract.from_configuration(
            payload["axis_profile_configuration"]
        )
        return cls(
            axis_profiles=axis_profiles,
            series_terms=int(payload["series_terms"]),
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
    ) -> "KokunoPA10Phi0ProfileContract":
        return cls.from_configuration(json.loads(Path(path).read_text()))

    def report(self) -> dict[str, Any]:
        probe_eta = np.linspace(-1.0, 1.0, 33)
        axis = self.values(np.zeros_like(probe_eta), probe_eta)
        off_axis = self.values(np.full_like(probe_eta, 4.1), probe_eta)
        z = off_axis["f0_argument"]
        f0 = off_axis["Phi_0"]
        fp = off_axis["f0_prime"]
        fpp = off_axis["f0_second"]
        recurrence_error = np.max(np.abs(2.0 * (z * fpp + 2.0 * fp) + f0))
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
            "coordinate_contract": {
                "native_similarity_coordinates": ["Y", "eta"],
                "Y_interval": list(self.source_Y_interval),
                "eta_interval": list(self.eta_interval),
                "repository_X_identification_asserted": False,
                "physical_F_from_Phi_mapping_resolved": False,
            },
            "numerical_realization": {
                "kind": "fixed finite polynomial from displayed entire f_0 series",
                "series_terms": int(self.series_terms),
                "maximum_probe_tail_upper": float(np.max(off_axis["f0_tail_upper"])),
            },
            "configuration": self.configuration(),
            "machine_checks": {
                "Phi0_axis_value_exact_one_on_probe": bool(
                    np.all(axis["Phi_0"] == 1.0)
                ),
                "Phi0_nontrivial_off_axis_on_probe": bool(
                    np.any(np.abs(off_axis["Phi_0"] - 1.0) > 0.0)
                ),
                "f0_source_ode_max_abs_error": float(recurrence_error),
                "all_profile_values_finite": bool(
                    all(np.all(np.isfinite(v)) for v in off_axis.values())
                ),
            },
            "truth_boundary": copy.deepcopy(_TRUTH_BOUNDARY),
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
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--config-output", type=Path)
    args = parser.parse_args()
    profiles = KokunoPA10Phi0ProfileContract()
    payload = profiles.save_report(args.output)
    if args.config_output is not None:
        profiles.save_configuration(args.config_output)
    print("receipt_sha256=", payload["receipt_sha256"])
    print("machine_checks=", payload["machine_checks"])


if __name__ == "__main__":
    _main()
