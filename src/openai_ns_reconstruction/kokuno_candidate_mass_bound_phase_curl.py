"""Bridge candidate h_sigma masses into Kokuno's phase-bound complete curls.

The h_sigma formula and the phase/curl formulas are source-displayed. Numerical
pulse samples, h_sigma derivatives, background values and mode prototypes remain
candidate inputs. This bridge therefore materializes candidate Q-scaled velocity
arrays without claiming source-exact data or PDE validation.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any
import math
import numpy as np

from .kokuno_source_band_covering import KokunoSourceBandCovering
from .kokuno_source_phase_bound_signed_curl import KokunoSourcePhaseBoundSignedCurlFamily

_ALLOWED_DERIVATIVES = {
    "candidate_derived_from_same_pulse_model",
    "explicit_frozen_local_zero",
}
_RESERVED = {
    "h_plus", "h_minus", "D_r_h_plus", "D_r_h_minus",
    "D_z_h_plus", "D_z_h_minus",
}


def _labels(value: Any) -> tuple[tuple[Any, ...], ...]:
    try:
        out = tuple(tuple(item) for item in value)
    except TypeError as exc:
        raise ValueError("beta labels must be a sequence") from exc
    if not out:
        raise ValueError("beta labels must not be empty")
    return out


def _matrix(value: Any, name: str, n_beta: int) -> np.ndarray:
    out = np.asarray(value, dtype=float)
    if out.shape != (n_beta, 2) or not np.all(np.isfinite(out)):
        raise ValueError(f"{name} must be finite with shape (n_beta,2)")
    return out


@dataclass(frozen=True)
class KokunoCandidateMassBoundPhaseCurlFamily:
    h: float = 0.005
    partition_atol: float = 2e-12
    derivative_rtol: float = 2e-10
    reality_atol: float = 2e-12
    covariance_consistency_atol: float = 2e-12
    mass_consistency_atol: float = 2e-12

    def __post_init__(self) -> None:
        family = self.phase_family
        atol = float(self.mass_consistency_atol)
        if not math.isfinite(atol) or not 0.0 < atol <= 1e-8:
            raise ValueError("mass_consistency_atol must lie in (0,1e-8]")
        object.__setattr__(self, "h", family.h)
        object.__setattr__(self, "mass_consistency_atol", atol)

    @property
    def phase_family(self) -> KokunoSourcePhaseBoundSignedCurlFamily:
        return KokunoSourcePhaseBoundSignedCurlFamily(
            h=self.h,
            partition_atol=self.partition_atol,
            derivative_rtol=self.derivative_rtol,
            reality_atol=self.reality_atol,
            covariance_consistency_atol=self.covariance_consistency_atol,
        )

    def _mass_values(self, mass: dict[str, Any], beta_labels: Any, L_s: Any) -> np.ndarray:
        required = (
            "beta_labels", "sign_labels", "Q_by_beta", "epsilon_by_beta",
            "L_s_by_beta", "h_sigma_by_beta_sign", "pulse_binding",
            "quadrature_binding", "source_formula_bound",
            "source_actual_h_sigma_numeric_recovered", "paper_exact",
        )
        if not isinstance(mass, dict) or any(k not in mass for k in required):
            raise ValueError("materialized_mass must be a complete signed covariance-mass result")
        if mass["source_formula_bound"] is not True:
            raise ValueError("mass artifact lost its displayed source-formula binding")
        if mass["source_actual_h_sigma_numeric_recovered"] is not False or mass["paper_exact"] is not False:
            raise ValueError("candidate mass artifact may not claim recovered source numerics or paper_exact")
        if tuple(mass["sign_labels"]) != ("sigma_plus", "sigma_minus"):
            raise ValueError("mass sign labels must be sigma_plus/sigma_minus")

        labels = _labels(beta_labels)
        if labels != _labels(mass["beta_labels"]):
            raise ValueError("phase/curl beta_labels must exactly match mass beta_labels")
        n = len(labels)
        h_sigma = _matrix(mass["h_sigma_by_beta_sign"], "h_sigma_by_beta_sign", n)
        if np.any(h_sigma <= 0.0):
            raise ValueError("h_sigma masses must be strictly positive")

        Q = np.asarray(mass["Q_by_beta"], dtype=float)
        eps = np.asarray(mass["epsilon_by_beta"], dtype=float)
        mass_L = np.asarray(mass["L_s_by_beta"], dtype=float)
        call_L = np.asarray(L_s, dtype=float)
        for name, arr in (("Q", Q), ("epsilon", eps), ("L_s", mass_L), ("call L_s", call_L)):
            if arr.shape != (n,) or not np.all(np.isfinite(arr)) or np.any(arr <= 0.0):
                raise ValueError(f"{name} must provide one positive finite value per beta")
        for j, label in enumerate(labels):
            schedule = KokunoSourceBandCovering(int(label[0]), self.h)
            for name, actual, expected in (("Q", Q[j], schedule.Q), ("epsilon", eps[j], schedule.epsilon)):
                scale = max(1.0, abs(float(expected)))
                if abs(float(actual) - float(expected)) > self.mass_consistency_atol * scale:
                    raise ValueError(f"mass {name} disagrees with the source band schedule")
            scale = max(1.0, abs(float(mass_L[j])))
            if abs(float(call_L[j]) - float(mass_L[j])) > self.mass_consistency_atol * scale:
                raise ValueError("phase/curl L_s must match mass L_s_by_beta")
        return h_sigma

    def physical_family_from_materialized_mass(
        self,
        materialized_mass: dict[str, Any],
        *,
        D_r_h_by_beta_sign: Any,
        D_z_h_by_beta_sign: Any,
        h_derivative_binding: str,
        **phase_curl_inputs: Any,
    ) -> dict[str, Any]:
        overlap = _RESERVED.intersection(phase_curl_inputs)
        if overlap:
            raise ValueError("h_sigma values and derivatives are owned by this bridge")
        if "beta_labels" not in phase_curl_inputs or "L_s" not in phase_curl_inputs:
            raise ValueError("phase_curl_inputs must include beta_labels and L_s")
        if h_derivative_binding not in _ALLOWED_DERIVATIVES:
            raise ValueError("h_derivative_binding must identify the candidate derivative route")

        labels = _labels(phase_curl_inputs["beta_labels"])
        masses = self._mass_values(materialized_mass, labels, phase_curl_inputs["L_s"])
        D_r_h = _matrix(D_r_h_by_beta_sign, "D_r_h_by_beta_sign", len(labels))
        D_z_h = _matrix(D_z_h_by_beta_sign, "D_z_h_by_beta_sign", len(labels))
        if h_derivative_binding == "explicit_frozen_local_zero" and (
            np.any(D_r_h != 0.0) or np.any(D_z_h != 0.0)
        ):
            raise ValueError("explicit_frozen_local_zero requires exactly zero h_sigma derivatives")

        out = self.phase_family.physical_family(
            **phase_curl_inputs,
            h_plus=masses[:, 0], h_minus=masses[:, 1],
            D_r_h_plus=D_r_h[:, 0], D_r_h_minus=D_r_h[:, 1],
            D_z_h_plus=D_z_h[:, 0], D_z_h_minus=D_z_h[:, 1],
        )
        return {
            **out,
            "candidate_h_sigma_by_beta_sign": masses.copy(),
            "candidate_D_r_h_by_beta_sign": D_r_h.copy(),
            "candidate_D_z_h_by_beta_sign": D_z_h.copy(),
            "h_sigma_pulse_binding": materialized_mass["pulse_binding"],
            "h_sigma_quadrature_binding": materialized_mass["quadrature_binding"],
            "h_sigma_derivative_binding": h_derivative_binding,
            "candidate_h_sigma_bound_into_complete_curl": True,
            "candidate_velocity_osc_cartesian_by_beta_sign": out["velocity_physical_cartesian_by_beta_sign"],
            "candidate_velocity_osc_cartesian_by_beta": out["velocity_physical_cartesian_by_beta"],
            "candidate_velocity_osc_cartesian_total": out["velocity_physical_cartesian_total"],
            "actual_source_h_sigma_pulse_integrals_bound": False,
            "actual_positive_order_background_bound": False,
            "actual_auxiliary_torus_mode_family_bound": False,
            "public_xyz_t_velocity_correction_materialized": False,
            "pde_validated": False,
            "paper_exact": False,
        }

    def receipt(self) -> dict[str, Any]:
        return {
            "schema": "kokuno-candidate-mass-bound-phase-curl-v1",
            "source_formula_chain": (
                "source h_sigma law -> candidate quadrature values -> explicit slow h_sigma derivatives -> "
                "source phase/frame -> localized complete curls -> Q_ell^(-A) physical velocity arrays"
            ),
            "numerical_boundary": {
                "h_sigma_values": "candidate numerical quadrature from the signed mass adapter",
                "h_sigma_derivatives": "explicit candidate inputs; never silently assumed zero",
                "background_and_modes": "caller supplied",
            },
            "truth_boundary": {
                "candidate_h_sigma_bound_into_complete_curl": True,
                "candidate_velocity_osc_arrays_materialized": True,
                "actual_source_h_sigma_pulse_integrals_bound": False,
                "actual_positive_order_background_bound": False,
                "actual_auxiliary_torus_mode_family_bound": False,
                "public_xyz_t_velocity_correction_materialized": False,
                "formal_full_domain_pde_gate_assessed": False,
                "pde_validated": False,
                "paper_exact": False,
            },
        }
