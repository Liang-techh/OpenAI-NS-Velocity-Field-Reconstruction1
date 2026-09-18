"""Agent-3 bridge from physical mean-covariance targets to Kokuno's signed inverse.

Agent 2's source-facing signed pair exposes the corrected reconstruction's
reference map

    H_ref y = T_0,                  y_sigma = a_sigma**2,
    C(W_0) = epsilon H_ref y.

For Agent 3 the important unit seam is therefore not ``H_ref y = Delta C``
but

    H_ref y = Delta C / epsilon.

This module makes that epsilon conversion explicit, keeps target provenance and
structural-vs-real-defect status visible, and refuses surrogate-defect promotion.
It remains a *reference-level* bridge: the actual h_sigma data, signed auxiliary
rectangles, complete-curl physical velocity and finite correction cycle are not
materialized here.
"""
from __future__ import annotations

from dataclasses import dataclass, field
import json
from typing import Any

import numpy as np

from .kokuno_source_signed_covariance_pair import KokunoSourceSignedCovariancePair

SCHEMA = "kokuno-agent3-source-signed-covariance-bridge-v1"
TARGET_UNITS = "physical_phase_mean_covariance_correction"
_TARGET_KINDS = {"structural_calibration", "real_candidate_phase_mean_defect"}


def _finite(value: Any, name: str) -> np.ndarray:
    out = np.asarray(value, dtype=float)
    if not np.all(np.isfinite(out)):
        raise ValueError(f"{name} must contain only finite real values")
    return out


@dataclass(frozen=True)
class KokunoSourceSignedCovarianceBridge:
    """Fail-closed Agent-3 reference inverse with explicit epsilon units."""

    pair: KokunoSourceSignedCovariancePair = field(default_factory=KokunoSourceSignedCovariancePair)
    reconstruction_atol: float = 2.0e-12

    def __post_init__(self) -> None:
        atol = float(self.reconstruction_atol)
        if not np.isfinite(atol) or not (0.0 < atol <= 1.0e-8):
            raise ValueError("reconstruction_atol must lie in (0,1e-8]")
        object.__setattr__(self, "reconstruction_atol", atol)

    def screen_target(
        self,
        *,
        A_c: Any,
        u_star: Any,
        h_plus: Any,
        h_minus: Any,
        epsilon: Any,
        target_N: Any,
        target_K: Any,
        target_kind: str,
        target_provenance: str,
        target_units: str = TARGET_UNITS,
        surrogate_defect_used: bool = False,
        direction_gap_eta: float | None = None,
    ) -> dict[str, Any]:
        """Map one physical covariance correction target into source y_sigma.

        ``target_N,target_K`` are the *physical* phase-mean covariance/stress
        correction requested from the primary wave.  Since the source identity
        is ``C(W_0)=epsilon H_ref y``, Agent 3 must solve against
        ``(target_N,target_K)/epsilon``.  This conversion is the main contract of
        this bridge.
        """
        if surrogate_defect_used:
            raise ValueError("surrogate defects cannot be promoted through the signed source bridge")
        if target_kind not in _TARGET_KINDS:
            raise ValueError(f"target_kind must be one of {sorted(_TARGET_KINDS)}")
        if not isinstance(target_provenance, str) or not target_provenance.strip():
            raise ValueError("target_provenance must be a nonempty explicit source/measurement label")
        if target_units != TARGET_UNITS:
            raise ValueError(f"target_units must equal {TARGET_UNITS!r}")

        eps = _finite(epsilon, "epsilon")
        if np.any(eps <= 0.0):
            raise ValueError("epsilon must be strictly positive")
        tN = _finite(target_N, "target_N")
        tK = _finite(target_K, "target_K")
        try:
            desired_shape = np.broadcast_shapes(eps.shape, tN.shape, tK.shape)
        except ValueError as exc:
            raise ValueError("epsilon and target components must be broadcast-compatible") from exc
        eps0 = np.broadcast_to(eps, desired_shape)
        desired = np.stack((np.broadcast_to(tN, desired_shape), np.broadcast_to(tK, desired_shape)), axis=-1)
        source_target = desired / eps0[..., None]

        inverse = self.pair.solve_reference_amplitudes(
            A_c=A_c,
            u_star=u_star,
            h_plus=h_plus,
            h_minus=h_minus,
            T_N=source_target[..., 0],
            T_K=source_target[..., 1],
            direction_gap_eta=direction_gap_eta,
        )

        full_shape = np.asarray(inverse["A_c"]).shape
        eps_full = np.broadcast_to(eps0, full_shape)
        desired_full = np.broadcast_to(desired, full_shape + (2,))
        reconstructed = eps_full[..., None] * np.asarray(inverse["reconstructed_target"])
        error = reconstructed - desired_full
        scale = np.maximum(1.0, np.max(np.abs(desired_full), axis=-1, keepdims=True))
        if not np.all(np.abs(error) <= self.reconstruction_atol * scale):
            raise RuntimeError("epsilon-scaled source inverse failed physical covariance reconstruction")

        p = np.asarray(inverse["p"])
        q = np.asarray(inverse["q"])
        cone_margin_fraction = (p - np.abs(q)) / p
        if np.any(~np.isfinite(cone_margin_fraction)) or np.any(cone_margin_fraction <= 0.0):
            raise RuntimeError("signed source inverse returned a nonpositive covariance-cone margin")

        singular_values = np.asarray(inverse["reference_column_singular_values"])
        condition_number = singular_values[..., 0] / singular_values[..., 1]
        amplitudes = np.asarray(inverse["amplitudes"])
        squared_amplitudes = np.asarray(inverse["squared_amplitudes"])

        real_target = target_kind == "real_candidate_phase_mean_defect"
        return {
            "schema": SCHEMA,
            "target_kind": target_kind,
            "target_provenance": target_provenance.strip(),
            "target_units": target_units,
            "surrogate_defect_used": False,
            "epsilon": eps_full,
            "physical_target": desired_full,
            "source_reference_target": np.asarray(inverse["target"]),
            "epsilon_scaling_identity": "physical_target = epsilon * H_ref * y",
            "reference_column_matrix": np.asarray(inverse["reference_column_matrix"]),
            "squared_amplitudes": squared_amplitudes,
            "amplitudes": amplitudes,
            "physical_reconstructed_target": reconstructed,
            "physical_reconstruction_error": error,
            "max_abs_physical_reconstruction_error": float(np.max(np.abs(error))),
            "cone_margin_fraction": cone_margin_fraction,
            "min_cone_margin_fraction": float(np.min(cone_margin_fraction)),
            "reference_column_singular_values": singular_values,
            "reference_column_condition_number": condition_number,
            "max_reference_column_condition_number": float(np.max(condition_number)),
            "reference_covariance_rank_two": bool(inverse["reference_covariance_rank_two"]),
            "positive_reference_coefficients": bool(inverse["positive_reference_coefficients"]),
            "reference_inverse_preflight_passed": True,
            "real_candidate_defect_consumed": real_target,
            "actual_source_h_sigma_pulse_integrals_bound": False,
            "actual_signed_auxiliary_rectangles_bound": False,
            "physical_complete_curl_signed_family_bound": False,
            "public_velocity_correction_materialized": False,
            "finite_correction_cycle_rerun_allowed": False,
            "heldout_ns_residual_assessed": False,
            "residual_reduction_claimed": False,
            "pde_validated": False,
        }


def structural_receipt() -> dict[str, Any]:
    """Deterministic source-structure receipt; deliberately not a real NS defect."""
    bridge = KokunoSourceSignedCovarianceBridge()
    epsilon = np.array([4.0e-2, 1.0e-2, 2.5e-3])
    out = bridge.screen_target(
        A_c=2.0,
        u_star=3.0,
        h_plus=1.0,
        h_minus=1.0,
        epsilon=epsilon,
        target_N=-8.0e-2,
        target_K=2.4e-2,
        target_kind="structural_calibration",
        target_provenance="repository_predeclared_signed_reference_unit_calibration",
        direction_gap_eta=0.2,
    )
    return {
        "schema": out["schema"],
        "epsilon": np.asarray(out["epsilon"]).tolist(),
        "squared_amplitudes": np.asarray(out["squared_amplitudes"]).tolist(),
        "amplitudes": np.asarray(out["amplitudes"]).tolist(),
        "min_cone_margin_fraction": out["min_cone_margin_fraction"],
        "max_reference_column_condition_number": out["max_reference_column_condition_number"],
        "max_abs_physical_reconstruction_error": out["max_abs_physical_reconstruction_error"],
        "reference_covariance_rank_two": out["reference_covariance_rank_two"],
        "reference_inverse_preflight_passed": out["reference_inverse_preflight_passed"],
        "real_candidate_defect_consumed": out["real_candidate_defect_consumed"],
        "physical_complete_curl_signed_family_bound": out["physical_complete_curl_signed_family_bound"],
        "finite_correction_cycle_rerun_allowed": out["finite_correction_cycle_rerun_allowed"],
        "residual_reduction_claimed": out["residual_reduction_claimed"],
        "pde_validated": out["pde_validated"],
    }


def main() -> None:
    print(json.dumps(structural_receipt(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
