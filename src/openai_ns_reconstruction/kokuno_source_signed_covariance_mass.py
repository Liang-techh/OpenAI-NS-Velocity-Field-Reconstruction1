"""Numerically materialize Kokuno's signed pulse covariance masses.

The corrected 2026-09-09 reconstruction displays, for each auxiliary-rectangle
sign sigma inside one slow beta box,

    I0_sigma = int_0^Ls psi(v)^2 x_sigma(v)^2 dv,
    D_g = (1+b_g^2)/2 * int_R chi_g(xi)^2 dxi,
    h_sigma = D_g c_i I0_sigma,
    L_s = 2 r0 / c_i.

The formula is source-displayed.  The actual numerical pulse x_sigma, cutoff
samples, and hidden source rectangle data are not recovered by this repository.
This module therefore evaluates the displayed mass law only on an explicitly
provenance-labelled candidate realization, and keeps that separate from source
truth.  The resulting positive h_+/h_- can be handed to the already existing
signed reference covariance inverse; doing so does not promote the candidate
pulse to paper-exact data.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
import math
from typing import Any, Sequence

import numpy as np

from .kokuno_physical_evaluation import (
    CORRECTED_RELEASE,
    CORRECTED_RELEASE_DATE,
    SOURCE_B_G,
    SOURCE_COMMIT,
    SOURCE_PATH,
    SOURCE_REPOSITORY,
)
from .kokuno_source_band_covering import KokunoSourceBandCovering
from .kokuno_source_signed_covariance_pair import KokunoSourceSignedCovariancePair

SCHEMA = "kokuno-source-signed-covariance-mass-v1"

_SOURCE_FORMULAS = {
    "pulse_mass": "I0_sigma=int_0^L_s psi(v)^2 x_sigma(v)^2 dv",
    "transverse_factor": "D_g=(1+b_g^2)/2 * int_R chi_g(xi)^2 dxi",
    "signed_mass": "h_sigma=D_g*c_i*I0_sigma",
    "pulse_length": "L_s=2*r0/c_i",
    "reference_columns": "H_ref=[[-A_c h_+,-A_c h_-],[-u_* h_+,+u_* h_-]]",
}

_TRUTH_BOUNDARY = {
    "source_h_sigma_integral_formula_executable": True,
    "source_displayed_D_g_formula_executable": True,
    "autonomous_geometry_can_feed_h_sigma_evaluation": True,
    "candidate_h_sigma_numerically_materializable": True,
    "candidate_h_sigma_can_feed_reference_covariance_pair": True,
    "actual_source_h_sigma_pulse_integrals_bound": False,
    "actual_source_pulse_samples_recovered": False,
    "actual_positive_order_background_bound": False,
    "actual_auxiliary_torus_mode_family_bound": False,
    "public_xyz_t_velocity_correction_materialized": False,
    "actual_source_physical_covariance_rank_two_assessed": False,
    "formal_full_domain_pde_gate_assessed": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}

_ALLOWED_PULSE_BINDINGS = {
    "repository_autonomous_source_compatible",
    "candidate_derived_from_executable_background",
}


def _canonical_json(payload: dict[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _finite_1d(value: Any, name: str, *, minimum_size: int = 2) -> np.ndarray:
    out = np.asarray(value, dtype=float)
    if out.ndim != 1 or out.size < minimum_size:
        raise ValueError(f"{name} must be a one-dimensional array with at least {minimum_size} samples")
    if not np.all(np.isfinite(out)):
        raise ValueError(f"{name} must contain only finite real values")
    return out


def _strict_grid(value: Any, name: str) -> np.ndarray:
    out = _finite_1d(value, name)
    if np.any(np.diff(out) <= 0.0):
        raise ValueError(f"{name} must be strictly increasing")
    return out


def _trapz_strict(grid: np.ndarray, values: np.ndarray) -> float:
    widths = np.diff(grid)
    return float(np.sum(0.5 * widths * (values[:-1] + values[1:])))


def _as_sequence(value: Sequence[Any], name: str, expected: int) -> tuple[Any, ...]:
    try:
        seq = tuple(value)
    except TypeError as exc:
        raise ValueError(f"{name} must be a sequence over beta labels") from exc
    if len(seq) != expected:
        raise ValueError(f"{name} length must match beta label count")
    return seq


@dataclass(frozen=True)
class KokunoSourceSignedCovarianceMass:
    """Evaluate source-displayed h_sigma on labelled numerical pulse samples."""

    h: float = 0.005
    consistency_atol: float = 2.0e-12

    def __post_init__(self) -> None:
        h = float(self.h)
        atol = float(self.consistency_atol)
        if not math.isfinite(h) or not (0.0 < h < 0.01):
            raise ValueError("h must satisfy the corrected-reader range 0<h<1/100")
        if not math.isfinite(atol) or not (0.0 < atol <= 1.0e-8):
            raise ValueError("consistency_atol must lie in (0,1e-8]")
        object.__setattr__(self, "h", h)
        object.__setattr__(self, "consistency_atol", atol)

    @staticmethod
    def sign_labels() -> tuple[str, str]:
        return ("sigma_plus", "sigma_minus")

    def transverse_factor_from_samples(self, xi: Any, chi_g: Any) -> dict[str, float]:
        xi_grid = _strict_grid(xi, "xi")
        chi = _finite_1d(chi_g, "chi_g")
        if chi.shape != xi_grid.shape:
            raise ValueError("chi_g samples must have the same shape as xi")
        chi_sq_integral = _trapz_strict(xi_grid, chi * chi)
        if not math.isfinite(chi_sq_integral) or chi_sq_integral <= 0.0:
            raise ValueError("integral chi_g^2 dxi must be strictly positive")
        D_g = 0.5 * (1.0 + SOURCE_B_G * SOURCE_B_G) * chi_sq_integral
        return {
            "chi_g_squared_integral": chi_sq_integral,
            "D_g": D_g,
            "quadrature_binding": "repository_numerical_trapezoid_not_source_exact",
        }

    def materialize_from_autonomous_geometry(
        self,
        geometry: dict[str, Any],
        *,
        v_by_beta: Sequence[Any],
        psi_by_beta: Sequence[Any],
        x_by_beta_sign: Sequence[Any],
        xi: Any,
        chi_g: Any,
        pulse_binding: str,
    ) -> dict[str, Any]:
        """Evaluate h_+/h_- for every beta in an Agent-2 geometry result.

        ``geometry`` is the public result of
        ``KokunoAutonomousSignedRectangleGeometry.instantiate``.  Pulse and
        cutoff arrays are numerical candidate inputs.  Their provenance must be
        named explicitly and is never upgraded to recovered source data.
        """
        if pulse_binding not in _ALLOWED_PULSE_BINDINGS:
            raise ValueError("pulse_binding must explicitly identify a permitted candidate realization")
        if not isinstance(geometry, dict):
            raise TypeError("geometry must be the public Agent-2 signed-rectangle geometry result")
        required = (
            "beta_labels",
            "Q_by_beta",
            "epsilon_by_beta",
            "L_s_by_beta",
            "covering_level_by_beta",
            "r0",
            "r0_binding",
            "rectangle_geometry_is_autonomous",
            "source_actual_rectangle_labels_instantiated",
        )
        if any(key not in geometry for key in required):
            raise ValueError("geometry is missing required signed-rectangle schedule fields")
        if geometry["rectangle_geometry_is_autonomous"] is not True:
            raise ValueError("this adapter expects the audited autonomous/source-compatible geometry route")
        if geometry["source_actual_rectangle_labels_instantiated"] is not False:
            raise ValueError("autonomous geometry may not be relabelled as actual source rectangles")

        labels = tuple(geometry["beta_labels"])
        n_beta = len(labels)
        if n_beta == 0:
            raise ValueError("geometry must contain at least one beta label")
        v_seq = _as_sequence(v_by_beta, "v_by_beta", n_beta)
        psi_seq = _as_sequence(psi_by_beta, "psi_by_beta", n_beta)
        x_seq = _as_sequence(x_by_beta_sign, "x_by_beta_sign", n_beta)

        Q = np.asarray(geometry["Q_by_beta"], dtype=float)
        epsilon = np.asarray(geometry["epsilon_by_beta"], dtype=float)
        L_s = np.asarray(geometry["L_s_by_beta"], dtype=float)
        levels = np.asarray(geometry["covering_level_by_beta"])
        if Q.shape != (n_beta,) or epsilon.shape != (n_beta,) or L_s.shape != (n_beta,) or levels.shape != (n_beta,):
            raise ValueError("geometry schedule arrays must have one entry per beta")
        if not np.all(np.isfinite(Q)) or not np.all(np.isfinite(epsilon)) or not np.all(np.isfinite(L_s)):
            raise ValueError("geometry schedule values must be finite")
        if np.any(Q <= 0.0) or np.any(epsilon <= 0.0) or np.any(L_s <= 0.0):
            raise ValueError("geometry Q, epsilon and L_s must be strictly positive")

        r0 = float(geometry["r0"])
        if not math.isfinite(r0) or r0 <= 0.0:
            raise ValueError("geometry r0 must be finite and positive")

        transverse = self.transverse_factor_from_samples(xi, chi_g)
        D_g = transverse["D_g"]
        I0 = np.empty((n_beta, 2), dtype=float)
        h_sigma = np.empty((n_beta, 2), dtype=float)
        c_i = np.empty(n_beta, dtype=float)

        for j, label in enumerate(labels):
            try:
                ell = int(label[0])
            except (TypeError, ValueError, IndexError) as exc:
                raise ValueError("each beta label must expose integer ell in position zero") from exc
            band = KokunoSourceBandCovering(ell=ell, h=self.h)
            c_i[j] = band.c_i
            expected_L = 2.0 * r0 / band.c_i
            scale = max(1.0, abs(expected_L))
            if abs(float(L_s[j]) - expected_L) > self.consistency_atol * scale:
                raise ValueError("geometry L_s is inconsistent with source L_s=2*r0/c_i")
            if abs(float(Q[j]) - band.Q) > self.consistency_atol * max(1.0, abs(band.Q)):
                raise ValueError("geometry Q is inconsistent with the source band schedule")
            if abs(float(epsilon[j]) - band.epsilon) > self.consistency_atol * max(1.0, abs(band.epsilon)):
                raise ValueError("geometry epsilon is inconsistent with the source band schedule")
            if int(levels[j]) != band.covering_level:
                raise ValueError("geometry covering level is inconsistent with the source band schedule")

            v = _strict_grid(v_seq[j], f"v_by_beta[{j}]")
            psi = _finite_1d(psi_seq[j], f"psi_by_beta[{j}]")
            x = np.asarray(x_seq[j], dtype=float)
            if psi.shape != v.shape:
                raise ValueError("psi samples must align with the beta v grid")
            if x.shape != (v.size, 2) or not np.all(np.isfinite(x)):
                raise ValueError("x_by_beta_sign entries must have shape (n_v,2) and be finite")
            if np.any(psi < -self.consistency_atol) or np.any(psi > 1.0 + self.consistency_atol):
                raise ValueError("psi must stay in the source cutoff range [0,1]")
            if np.any(x <= 0.0):
                raise ValueError("homogeneous radial pulse samples x_sigma must be strictly positive")
            endpoint_scale = max(1.0, float(L_s[j]))
            if abs(float(v[0])) > self.consistency_atol * endpoint_scale:
                raise ValueError("pulse v grid must start at zero")
            if abs(float(v[-1]) - float(L_s[j])) > self.consistency_atol * endpoint_scale:
                raise ValueError("pulse v grid must end at geometry L_s")

            for sign_index in range(2):
                integrand = psi * psi * x[:, sign_index] * x[:, sign_index]
                I0[j, sign_index] = _trapz_strict(v, integrand)
                h_sigma[j, sign_index] = D_g * band.c_i * I0[j, sign_index]
            if np.any(I0[j] <= 0.0) or np.any(h_sigma[j] <= 0.0):
                raise ValueError("source covariance masses must be strictly positive")

        return {
            "beta_labels": labels,
            "sign_labels": self.sign_labels(),
            "Q_by_beta": Q.copy(),
            "epsilon_by_beta": epsilon.copy(),
            "L_s_by_beta": L_s.copy(),
            "c_i_by_beta": c_i,
            "I0_by_beta_sign": I0,
            "h_sigma_by_beta_sign": h_sigma,
            "D_g": D_g,
            "chi_g_squared_integral": transverse["chi_g_squared_integral"],
            "pulse_binding": pulse_binding,
            "quadrature_binding": transverse["quadrature_binding"],
            "geometry_r0_binding": geometry["r0_binding"],
            "source_formula_bound": True,
            "source_actual_h_sigma_numeric_recovered": False,
            "paper_exact": False,
        }

    def reference_columns(
        self,
        materialized: dict[str, Any],
        *,
        A_c: Any,
        u_star: Any,
    ) -> dict[str, Any]:
        """Feed candidate h_+/h_- into the existing source reference-column map."""
        if not isinstance(materialized, dict) or "h_sigma_by_beta_sign" not in materialized:
            raise ValueError("materialized signed covariance masses are required")
        masses = np.asarray(materialized["h_sigma_by_beta_sign"], dtype=float)
        if masses.ndim != 2 or masses.shape[1] != 2 or np.any(masses <= 0.0):
            raise ValueError("h_sigma_by_beta_sign must have shape (n_beta,2) and be positive")
        n_beta = masses.shape[0]
        A = np.broadcast_to(np.asarray(A_c, dtype=float), (n_beta,))
        u = np.broadcast_to(np.asarray(u_star, dtype=float), (n_beta,))
        if not np.all(np.isfinite(A)) or not np.all(np.isfinite(u)) or np.any(A <= 0.0) or np.any(u <= 0.0):
            raise ValueError("A_c and u_star must broadcast to positive finite beta arrays")
        pair = KokunoSourceSignedCovariancePair(consistency_atol=self.consistency_atol)
        H = pair.reference_column_matrix(A, u, masses[:, 0], masses[:, 1])
        determinant = np.linalg.det(H)
        expected = -2.0 * A * u * masses[:, 0] * masses[:, 1]
        if not np.allclose(determinant, expected, rtol=0.0, atol=self.consistency_atol * np.maximum(1.0, np.abs(expected))):
            raise RuntimeError("reference determinant failed the displayed source identity")
        return {
            "reference_column_matrix_by_beta": H,
            "reference_column_determinant_by_beta": determinant,
            "reference_formula_only": True,
            "actual_complete_curl_covariance_rank_assessed": False,
            "actual_source_mode_family_bound": False,
        }

    def receipt(self) -> dict[str, Any]:
        payload = {
            "schema": SCHEMA,
            "source": {
                "repository": SOURCE_REPOSITORY,
                "commit": SOURCE_COMMIT,
                "path": SOURCE_PATH,
                "release": CORRECTED_RELEASE,
                "release_date": CORRECTED_RELEASE_DATE,
                "reader_location": "expanded corrected reader around lines 9503-9554",
            },
            "source_formulas": dict(_SOURCE_FORMULAS),
            "numerical_contract": {
                "quadrature": "composite trapezoid on caller/candidate supplied samples",
                "pulse_bindings": sorted(_ALLOWED_PULSE_BINDINGS),
                "geometry": "consumes audited repository-autonomous/source-compatible Agent-2 geometry",
                "source_hidden_numeric_values_are_not_inferred": True,
            },
            "truth_boundary": dict(_TRUTH_BOUNDARY),
        }
        payload["sha256"] = hashlib.sha256(_canonical_json(payload).encode("utf-8")).hexdigest()
        return payload
