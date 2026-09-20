"""Candidate-facing source-native PA.10 leading-axis profile contract.

Pinned public provenance is KokunoYumeto/yang-mills-interacting-workbench at
``143f6773feb424ad9ed3a8d116653200f20346b7``,
``navier-stokes/navier_stokes_workbench.tex`` (corrected reader dated
2026-09-09, Zenodo 22678406).

This increment collects the already-used public axis/base formulas behind one
small executable profile API on the source-native similarity rectangle
``0 <= Y <= 4.1`` and the certified enlarged real eta interval.  It exposes

    U_*=4 eta+j_0,
    H_*=D eta+(1-eta^2) U_*,
    W_*=1-(1-eta^2) U_*' - 2D eta U_*,
    chi=H_*^2/(H_*^2+sigma_*^2),
    zeta_*=-L H_*/(H_*^2+sigma_*^2),
    u_0(Y,eta)=-Y Z_*/(2L),

together with analytic eta/radial derivatives needed by a later candidate
evaluator.  The existing source-compatible axis-domain object remains the
single value source for d,L,U_*,H_*,Z_*,Pi_0,chi,zeta_*; this module adds only
the candidate-facing assembly, derivatives, and deterministic configuration
serialization.

``Y`` is the reconstruction's native inner radial similarity coordinate.
No claim ``Y == repository X`` is made here; the fixed axis profiles are
Y-independent and ``u_0`` is source-native.  The selected sigma/rho/pressure
datum and outer free parameters are explicit repository-autonomous
source-compatible choices, not recovered hidden Kokuno/OpenAI values.
``u_0`` is the contraction center, not the final corrected velocity.  Global
velocity assembly and held-out PDE validation remain fail-closed.
"""

from __future__ import annotations

import copy
from dataclasses import dataclass, field
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from .kokuno_outer_reserved_patch_schedule import KokunoOuterReservedPatchSchedule
from .kokuno_pressure_datum_binding import KokunoPressureDatumBinding
from .kokuno_pa10_source_axis_domain import KokunoPA10SourceCompatibleAxisDomain


SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
CORRECTED_RELEASE = "zenodo:22678406"
CORRECTED_RELEASE_DATE = "2026-09-09"
SCHEMA = "kokuno-pa10-leading-axis-profile-contract-v1"

_SOURCE_FORMULAS = {
    "d": "d=1-eta^2",
    "L": "L=1-2h eta^2",
    "U_star": "U_*=4eta+j_0",
    "H_star": "H_*=Deta+dU_*",
    "W_star": "W_*=1-dU_*'-2DetaU_*; U_*'=4",
    "chi": "chi=H_*^2/(H_*^2+sigma_*^2)",
    "zeta_star": "zeta_*=-L H_*/(H_*^2+sigma_*^2)",
    "u0": "u_0(Y,eta)=-Y Z_*(eta)/(2L(eta))",
    "Z_star": (
        "Z_*=-A(1-2etaU_*)U_*-H_*U_*'-(1-eta^2)Pi_0'"
        "+4AetaPi_0; U_*'=4"
    ),
    "Pi_0": "Pi_0=-P^2/(1+eta^2)^2 on the selected source-compatible datum",
}

_TRUTH_BOUNDARY = {
    "public_reconstruction_source": True,
    "source_native_leading_axis_profile_api_executable": True,
    "source_native_profile_vectorized": True,
    "source_native_profile_analytic_derivatives_executable": True,
    "source_native_profile_configuration_serializable": True,
    "axis_regular_u0_center_executable": True,
    "source_u0_is_contraction_center_not_final_corrected_u": True,
    "repository_X_identified_with_source_Y": False,
    "autonomous_source_compatible_parameter_values_explicit": True,
    "source_hidden_numeric_choices_recovered": False,
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
class KokunoPA10LeadingAxisProfileContract:
    """Executable source-native fixed leading profiles plus the PA.10 u_0 center."""

    domain: KokunoPA10SourceCompatibleAxisDomain = field(
        default_factory=KokunoPA10SourceCompatibleAxisDomain
    )

    def __post_init__(self) -> None:
        if not isinstance(self.domain, KokunoPA10SourceCompatibleAxisDomain):
            raise TypeError("domain must be KokunoPA10SourceCompatibleAxisDomain")

    @property
    def source_Y_interval(self) -> tuple[float, float]:
        return (0.0, 4.1)

    @property
    def eta_interval(self) -> tuple[float, float]:
        E = 1.0 + float(self.domain.enlarged_real_margin)
        return (-E, E)

    def _broadcast(self, Y: Any, eta: Any) -> tuple[np.ndarray, np.ndarray]:
        Y_arr, eta_arr = np.broadcast_arrays(
            _finite_array(Y, "Y"), _finite_array(eta, "eta")
        )
        y0, y1 = self.source_Y_interval
        if np.any((Y_arr < y0) | (Y_arr > y1)):
            raise ValueError(f"Y must lie in [{y0},{y1}]")
        e0, e1 = self.eta_interval
        if np.any((eta_arr < e0) | (eta_arr > e1)):
            raise ValueError("eta lies outside the certified enlarged real interval")
        return Y_arr, eta_arr

    def values(self, Y: Any, eta: Any) -> dict[str, np.ndarray]:
        """Return vectorized source-native leading-axis/base values.

        The fixed multipliers are independent of Y.  ``u_0`` is the source
        PA.10 contraction center and vanishes linearly at the radial axis.
        """
        Y_arr, eta_arr = self._broadcast(Y, eta)
        state = self.domain.axis_state(eta_arr)
        h = float(self.domain.h)
        D = float(self.domain.D)
        j0 = float(self.domain.j0)
        U = state["U_star"]
        d = state["d"]
        W = 1.0 - 4.0 * d - 2.0 * D * eta_arr * U
        slope = -state["Z_star"] / (2.0 * state["L"])
        return {
            "Y": Y_arr,
            "eta": eta_arr,
            "d": state["d"],
            "L": state["L"],
            "U_star": state["U_star"],
            "H_star": state["H_star"],
            "W_star": W,
            "Z_star": state["Z_star"],
            "Pi_0": state["Pi_0"],
            "Pi_0_eta": state["Pi_0_eta"],
            "chi": state["chi"],
            "zeta_star": state["zeta_star"],
            "u_0": Y_arr * slope,
            "u_0_Y": slope,
            "U_star_eta": np.full_like(eta_arr, 4.0, dtype=float),
            "j0": np.full_like(eta_arr, j0, dtype=float),
            "h": np.full_like(eta_arr, h, dtype=float),
        }

    def derivatives(self, Y: Any, eta: Any) -> dict[str, np.ndarray]:
        """Return analytic eta/radial derivatives of the exported profiles."""
        Y_arr, eta_arr = self._broadcast(Y, eta)
        state = self.domain.axis_state(eta_arr)
        h = float(self.domain.h)
        A = float(self.domain.A)
        D = float(self.domain.D)
        j0 = float(self.domain.j0)
        p2 = float(self.domain.pressure_square)
        sigma2 = float(self.domain.sigma_star) ** 2

        U = state["U_star"]
        H = state["H_star"]
        L = state["L"]
        Z = state["Z_star"]
        d = state["d"]
        Pi = state["Pi_0"]
        Pi_eta = state["Pi_0_eta"]

        d_eta = -2.0 * eta_arr
        L_eta = -4.0 * h * eta_arr
        U_eta = np.full_like(eta_arr, 4.0, dtype=float)
        H_eta = D - 2.0 * eta_arr * U + 4.0 * d
        W_eta = -2.0 * D * j0 + 16.0 * h * eta_arr

        den = 1.0 + eta_arr * eta_arr
        Pi_eta_eta = 4.0 * p2 * (1.0 - 5.0 * eta_arr * eta_arr) / den**4
        B = 1.0 - 2.0 * eta_arr * U
        B_eta = -2.0 * U - 8.0 * eta_arr
        Z_eta = (
            -A * (B_eta * U + 4.0 * B)
            - 4.0 * H_eta
            - d_eta * Pi_eta
            - d * Pi_eta_eta
            + 4.0 * A * (Pi + eta_arr * Pi_eta)
        )

        Q = H * H + sigma2
        chi_eta = 2.0 * H * H_eta * sigma2 / Q**2
        zeta_eta = (
            -(L_eta * H + L * H_eta) * Q
            + 2.0 * L * H * H * H_eta
        ) / Q**2

        slope = -Z / (2.0 * L)
        slope_eta = -(Z_eta * L - Z * L_eta) / (2.0 * L * L)
        return {
            "d_eta": d_eta,
            "L_eta": L_eta,
            "U_star_eta": U_eta,
            "H_star_eta": H_eta,
            "W_star_eta": W_eta,
            "Z_star_eta": Z_eta,
            "Pi_0_eta_eta": Pi_eta_eta,
            "chi_eta": chi_eta,
            "zeta_star_eta": zeta_eta,
            "u_0_Y": slope,
            "u_0_eta": Y_arr * slope_eta,
            "u_0_Y_eta": slope_eta,
        }

    def configuration(self) -> dict[str, Any]:
        """Return a deterministic, complete scalar configuration for reload."""
        schedule = self.domain.binding.outer_schedule
        binding = self.domain.binding
        return {
            "schema": SCHEMA,
            "outer_schedule": {
                "M_d": float(schedule.M_d),
                "log_p_star_margin": float(schedule.log_p_star_margin),
                "C": float(schedule.C),
                "lambda_outer": float(schedule.lambda_outer),
                "h": float(schedule.h),
                "repair_log_offset": float(schedule.repair_log_offset),
            },
            "pressure_binding": {
                "pressure_margin": float(binding.pressure_margin),
                "j0": float(binding.j0),
                "sigma": float(binding.sigma),
                "Lambda": float(binding.Lambda),
                "log_transition_width": float(binding.log_transition_width),
                "quadrature_points": int(binding.quadrature_points),
            },
            "source_axis_domain": {
                "sigma_star": float(self.domain.sigma_star),
                "delta_star": float(self.domain.delta_star),
                "enlarged_real_margin": float(self.domain.enlarged_real_margin),
                "partition_intervals": int(self.domain.partition_intervals),
                "complex_tube_radius": float(self.domain.complex_tube_radius),
                "cauchy_radius": float(self.domain.cauchy_radius),
                "coefficient_rho": float(self.domain.coefficient_rho),
            },
        }

    @classmethod
    def from_configuration(
        cls, payload: Mapping[str, Any]
    ) -> "KokunoPA10LeadingAxisProfileContract":
        if not isinstance(payload, Mapping):
            raise TypeError("configuration must be a mapping")
        if payload.get("schema") != SCHEMA:
            raise ValueError("configuration schema mismatch")
        outer = dict(payload["outer_schedule"])
        pressure = dict(payload["pressure_binding"])
        domain_values = dict(payload["source_axis_domain"])
        schedule = KokunoOuterReservedPatchSchedule(**outer)
        binding = KokunoPressureDatumBinding(outer_schedule=schedule, **pressure)
        domain = KokunoPA10SourceCompatibleAxisDomain(
            binding=binding, **domain_values
        )
        return cls(domain=domain)

    def save_configuration(self, path: str | Path) -> dict[str, Any]:
        payload = self.configuration()
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        return payload

    @classmethod
    def load_configuration(
        cls, path: str | Path
    ) -> "KokunoPA10LeadingAxisProfileContract":
        payload = json.loads(Path(path).read_text())
        return cls.from_configuration(payload)

    def report(self) -> dict[str, Any]:
        probe_eta = np.asarray([-1.0, 0.0, 1.0], dtype=float)
        axis = self.values(np.zeros_like(probe_eta), probe_eta)
        off_axis = self.values(np.ones_like(probe_eta), probe_eta)
        derivatives = self.derivatives(np.ones_like(probe_eta), probe_eta)
        identity_error = np.max(
            np.abs(
                axis["H_star"] * axis["zeta_star"] / axis["L"] + axis["chi"]
            )
        )
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
                "fixed_axis_profiles_are_Y_independent": True,
                "u_0_is_source_PA10_contraction_center": True,
            },
            "configuration": self.configuration(),
            "selected_choice_status": (
                "repository-autonomous source-compatible existence choice; "
                "not recovered hidden source data"
            ),
            "machine_checks": {
                "axis_regular_u0_exact_zero_on_probe": bool(
                    np.all(axis["u_0"] == 0.0)
                ),
                "u0_nontrivial_off_axis_on_probe": bool(
                    np.any(np.abs(off_axis["u_0"]) > 0.0)
                ),
                "H_zeta_over_L_equals_minus_chi_abs_error": float(identity_error),
                "all_profile_values_finite": bool(
                    all(np.all(np.isfinite(v)) for v in off_axis.values())
                ),
                "all_analytic_derivatives_finite": bool(
                    all(np.all(np.isfinite(v)) for v in derivatives.values())
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
    profiles = KokunoPA10LeadingAxisProfileContract()
    payload = profiles.save_report(args.output)
    if args.config_output is not None:
        profiles.save_configuration(args.config_output)
    print("receipt_sha256=", payload["receipt_sha256"])
    print("machine_checks=", payload["machine_checks"])


if __name__ == "__main__":
    _main()
