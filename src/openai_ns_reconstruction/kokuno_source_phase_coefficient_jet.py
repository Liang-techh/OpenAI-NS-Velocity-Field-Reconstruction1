"""Source-shaped phase and complete-curl coefficient directional jets.

This module is a narrow continuation of the executable corrected-source complete
curl reference.  It removes one caller-supplied layer: ``D_r C_m`` and
``D_z C_m`` are now derived analytically from the corrected-source phase-vector
jet and a transverse amplitude directional jet.

The corrected source uses, on an auxiliary band rectangle,

    Phi = p theta + p_z Z/epsilon + x_0 R - v H_Phi,
    H_Phi = p F + p_z G,
    n_Phi = (x_0-v(H_Phi)_R, p/R,
             p_z-epsilon v(H_Phi)_Z),

with ``D_r v = D_z v = 0`` and auxiliary-independent base fields.  Consequently
``D_r`` acts as ``partial_R`` on the base and ``D_z = epsilon partial_Z``.
This module executes those phase/vector first directional jets from declared
background derivatives, then differentiates

    C_m = i (n_Phi x t_m) / (k m |n_Phi|^2)

by the exact product/quotient rule.

It does *not* solve the source amplitude ODE for ``t_m`` and therefore still
requires ``t_m, D_r t_m, D_z t_m`` from an upstream source-amplitude
realization.  It also does not map the repository's autonomous oscillatory
runtime to the source chart/support, create a Cartesian project-domain velocity,
materialize pressure/forcing, or assess a Navier--Stokes residual.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np

from openai_ns_reconstruction.kokuno_corrected_oscillation_source_ledger import (
    SOURCE_COMMIT,
    SOURCE_REPOSITORY,
    SOURCE_WORKBENCH_BLOB,
    default_kokuno_corrected_oscillation_source_ledger,
)
from openai_ns_reconstruction.kokuno_source_normalized_complete_curl_reference import (
    CompleteHarmonicAmplitude,
    complete_harmonic_amplitude_from_coefficient_jet,
    harmonic_vector_potential_coefficient,
)


PARENT_A2_HEAD = "ecc2be96eee441308cc0216bf187b43bbd7fd3a0"
PARENT_COMPLETE_CURL_REFERENCE_BLOB = "42809fc1b10937b4246c65121ac1dbfdc2986822"
DIRECTIONAL_TRANSVERSALITY_ATOL = 2.0e-11
DIRECTIONAL_TRANSVERSALITY_RTOL = 2.0e-11


@dataclass(frozen=True)
class SourceBackgroundPhaseJet:
    """Auxiliary-independent base values needed for ``Phi`` and ``n_Phi`` jets.

    Derivative suffixes ``_r`` and ``_z`` denote ordinary ``partial_R`` and
    ``partial_Z`` derivatives of the imported base fields.  The conversion to
    the normalized source operator ``D_z=epsilon partial_Z`` is performed by
    :func:`source_phase_vector_jet`.
    """

    f: Any
    g: Any
    f_r: Any
    g_r: Any
    f_z: Any
    g_z: Any
    f_rr: Any
    g_rr: Any
    f_rz: Any
    g_rz: Any
    f_zz: Any
    g_zz: Any


@dataclass(frozen=True)
class SourcePhaseVectorJet:
    """Corrected-source phase, phase vector and normalized directional jets."""

    phase: np.ndarray
    h_phi: np.ndarray
    n_phi: np.ndarray
    dr_n_phi: np.ndarray
    dz_n_phi: np.ndarray


@dataclass(frozen=True)
class SourceAmplitudeDirectionalJet:
    """One transverse harmonic amplitude and its normalized directional jets."""

    t_m: Any
    dr_t_m: Any
    dz_t_m: Any


@dataclass(frozen=True)
class SourceCoefficientDirectionalJet:
    """``C_m`` together with analytically derived ``D_r C_m`` / ``D_z C_m``."""

    coefficient: np.ndarray
    dr_coefficient: np.ndarray
    dz_coefficient: np.ndarray


@dataclass(frozen=True)
class SourceCompleteHarmonicFromPhaseJet:
    """Phase/coefficient directional jets plus the executable complete amplitude."""

    phase_jet: SourcePhaseVectorJet
    coefficient_jet: SourceCoefficientDirectionalJet
    complete_amplitude: CompleteHarmonicAmplitude


def _finite(value: Any, *, name: str, dtype: Any = float) -> np.ndarray:
    array = np.asarray(value, dtype=dtype)
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must be finite")
    return array


def _positive_radius(radius: Any) -> np.ndarray:
    array = _finite(radius, name="radius")
    if np.any(array <= 0.0):
        raise ValueError("corrected-source phase chart requires radius > 0")
    return array


def _positive_epsilon(epsilon: float) -> float:
    if not np.isfinite(epsilon) or epsilon <= 0.0:
        raise ValueError("epsilon must be finite and positive")
    return float(epsilon)


def _finite_scalar(value: float, *, name: str) -> float:
    if not np.isfinite(value):
        raise ValueError(f"{name} must be finite")
    return float(value)


def _vector(value: Any, *, name: str, dtype: Any) -> np.ndarray:
    array = _finite(value, name=name, dtype=dtype)
    if array.ndim == 0 or array.shape[-1] != 3:
        raise ValueError(f"{name} must have final dimension 3")
    return array


def source_phase_vector_jet(
    radius: Any,
    theta: Any,
    z_normalized: Any,
    pulse_v: Any,
    background: SourceBackgroundPhaseJet,
    *,
    epsilon: float,
    p: float,
    p_z: float,
    x_0: float,
) -> SourcePhaseVectorJet:
    """Evaluate the corrected-source ``Phi``, ``n_Phi`` and first ``D_r/D_z`` jets.

    This is the explicit auxiliary-independent specialization used in the source
    band geometry: ``D_r v=D_z v=0``, ``D_r=partial_R`` on the imported base and
    ``D_z=epsilon partial_Z``.  No finite-difference step is introduced.
    """

    radius = _positive_radius(radius)
    epsilon = _positive_epsilon(epsilon)
    p = _finite_scalar(p, name="p")
    p_z = _finite_scalar(p_z, name="p_z")
    x_0 = _finite_scalar(x_0, name="x_0")

    values = [
        _finite(theta, name="theta"),
        _finite(z_normalized, name="z_normalized"),
        _finite(pulse_v, name="pulse_v"),
        _finite(background.f, name="background.f"),
        _finite(background.g, name="background.g"),
        _finite(background.f_r, name="background.f_r"),
        _finite(background.g_r, name="background.g_r"),
        _finite(background.f_z, name="background.f_z"),
        _finite(background.g_z, name="background.g_z"),
        _finite(background.f_rr, name="background.f_rr"),
        _finite(background.g_rr, name="background.g_rr"),
        _finite(background.f_rz, name="background.f_rz"),
        _finite(background.g_rz, name="background.g_rz"),
        _finite(background.f_zz, name="background.f_zz"),
        _finite(background.g_zz, name="background.g_zz"),
    ]
    (
        radius,
        theta,
        z_normalized,
        pulse_v,
        f,
        g,
        f_r,
        g_r,
        f_z,
        g_z,
        f_rr,
        g_rr,
        f_rz,
        g_rz,
        f_zz,
        g_zz,
    ) = np.broadcast_arrays(radius, *values)

    h_phi = p * f + p_z * g
    h_r = p * f_r + p_z * g_r
    h_z = p * f_z + p_z * g_z
    h_rr = p * f_rr + p_z * g_rr
    h_rz = p * f_rz + p_z * g_rz
    h_zz = p * f_zz + p_z * g_zz

    phase = (
        p * theta
        + p_z * z_normalized / epsilon
        + x_0 * radius
        - pulse_v * h_phi
    )
    n_phi = np.stack(
        (
            x_0 - pulse_v * h_r,
            p / radius,
            p_z - epsilon * pulse_v * h_z,
        ),
        axis=-1,
    )
    dr_n_phi = np.stack(
        (
            -pulse_v * h_rr,
            -p / radius**2,
            -epsilon * pulse_v * h_rz,
        ),
        axis=-1,
    )
    dz_n_phi = np.stack(
        (
            -epsilon * pulse_v * h_rz,
            np.zeros_like(radius),
            -(epsilon**2) * pulse_v * h_zz,
        ),
        axis=-1,
    )

    norm_sq = np.sum(n_phi * n_phi, axis=-1)
    if np.any(norm_sq <= 0.0) or not np.all(np.isfinite(norm_sq)):
        raise ValueError("source phase vector n_phi must be nonzero")

    return SourcePhaseVectorJet(
        phase=phase,
        h_phi=h_phi,
        n_phi=n_phi,
        dr_n_phi=dr_n_phi,
        dz_n_phi=dz_n_phi,
    )


def _check_directional_transversality(
    n_phi: np.ndarray,
    dn_phi: np.ndarray,
    t_m: np.ndarray,
    dt_m: np.ndarray,
    *,
    label: str,
) -> None:
    derivative = np.sum(dn_phi * t_m + n_phi * dt_m, axis=-1)
    scale = np.maximum(
        1.0,
        np.linalg.norm(dn_phi, axis=-1) * np.linalg.norm(t_m, axis=-1)
        + np.linalg.norm(n_phi, axis=-1) * np.linalg.norm(dt_m, axis=-1),
    )
    allowed = DIRECTIONAL_TRANSVERSALITY_ATOL + DIRECTIONAL_TRANSVERSALITY_RTOL * scale
    if np.any(np.abs(derivative) > allowed):
        raise ValueError(f"{label} must preserve D(n_phi dot t_m)=0")


def _differentiate_coefficient(
    n_phi: np.ndarray,
    dn_phi: np.ndarray,
    t_m: np.ndarray,
    dt_m: np.ndarray,
    *,
    k: float,
    m: int,
) -> np.ndarray:
    norm_sq = np.sum(n_phi * n_phi, axis=-1)
    cross = np.cross(n_phi, t_m)
    d_cross = np.cross(dn_phi, t_m) + np.cross(n_phi, dt_m)
    d_norm_sq = 2.0 * np.sum(n_phi * dn_phi, axis=-1)
    prefactor = 1j / float(k * m)
    return prefactor * (
        d_cross / norm_sq[..., None]
        - cross * d_norm_sq[..., None] / (norm_sq**2)[..., None]
    )


def harmonic_coefficient_directional_jet(
    phase_jet: SourcePhaseVectorJet,
    amplitude_jet: SourceAmplitudeDirectionalJet,
    *,
    k: float,
    m: int,
) -> SourceCoefficientDirectionalJet:
    """Derive ``C_m, D_r C_m, D_z C_m`` by the exact source chain rule."""

    n = _vector(phase_jet.n_phi, name="phase_jet.n_phi", dtype=float)
    dr_n = _vector(phase_jet.dr_n_phi, name="phase_jet.dr_n_phi", dtype=float)
    dz_n = _vector(phase_jet.dz_n_phi, name="phase_jet.dz_n_phi", dtype=float)
    t = _vector(amplitude_jet.t_m, name="amplitude_jet.t_m", dtype=np.complex128)
    dr_t = _vector(
        amplitude_jet.dr_t_m, name="amplitude_jet.dr_t_m", dtype=np.complex128
    )
    dz_t = _vector(
        amplitude_jet.dz_t_m, name="amplitude_jet.dz_t_m", dtype=np.complex128
    )
    n, dr_n, dz_n, t, dr_t, dz_t = np.broadcast_arrays(
        n, dr_n, dz_n, t, dr_t, dz_t
    )

    coefficient = harmonic_vector_potential_coefficient(n, t, k=k, m=m)
    _check_directional_transversality(n, dr_n, t, dr_t, label="D_r amplitude jet")
    _check_directional_transversality(n, dz_n, t, dz_t, label="D_z amplitude jet")

    dr_coefficient = _differentiate_coefficient(n, dr_n, t, dr_t, k=k, m=m)
    dz_coefficient = _differentiate_coefficient(n, dz_n, t, dz_t, k=k, m=m)
    return SourceCoefficientDirectionalJet(
        coefficient=coefficient,
        dr_coefficient=dr_coefficient,
        dz_coefficient=dz_coefficient,
    )


def complete_harmonic_from_source_phase_jet(
    radius: Any,
    phase_jet: SourcePhaseVectorJet,
    amplitude_jet: SourceAmplitudeDirectionalJet,
    *,
    k: float,
    m: int,
) -> SourceCompleteHarmonicFromPhaseJet:
    """Execute the #1018 complete-curl amplitude without caller-supplied C jets."""

    coefficient_jet = harmonic_coefficient_directional_jet(
        phase_jet,
        amplitude_jet,
        k=k,
        m=m,
    )
    complete = complete_harmonic_amplitude_from_coefficient_jet(
        radius,
        phase_jet.n_phi,
        amplitude_jet.t_m,
        coefficient_jet.dr_coefficient,
        coefficient_jet.dz_coefficient,
        k=k,
        m=m,
    )
    if not np.allclose(
        complete.coefficient,
        coefficient_jet.coefficient,
        rtol=0.0,
        atol=2.0e-14,
    ):
        raise RuntimeError("complete-curl coefficient drifted from source chain-rule jet")
    return SourceCompleteHarmonicFromPhaseJet(
        phase_jet=phase_jet,
        coefficient_jet=coefficient_jet,
        complete_amplitude=complete,
    )


def source_phase_coefficient_jet_contract() -> dict[str, Any]:
    """Fail-closed provenance and delivery boundary for this source adapter."""

    ledger = default_kokuno_corrected_oscillation_source_ledger()
    if ledger.phase != "Phi = p theta + p_z Z/epsilon + x_0 R - v H_Phi":
        raise RuntimeError("corrected-source phase formula drifted")
    if ledger.phase_vector != (
        "(n_Phi)_r = x_0 - v (H_Phi)_R",
        "(n_Phi)_theta = p/R",
        "(n_Phi)_z = p_z - epsilon v (H_Phi)_Z",
    ):
        raise RuntimeError("corrected-source phase-vector formula drifted")

    return {
        "schema": "kokuno-source-phase-coefficient-jet-v1",
        "parent_a2_head": PARENT_A2_HEAD,
        "parent_complete_curl_reference_blob": PARENT_COMPLETE_CURL_REFERENCE_BLOB,
        "source_repository": SOURCE_REPOSITORY,
        "source_commit": SOURCE_COMMIT,
        "source_workbench_blob": SOURCE_WORKBENCH_BLOB,
        "source_phase_vector_jet_executable": True,
        "source_coefficient_directional_jets_derived_analytically": True,
        "source_complete_curl_consumes_derived_coefficient_jets": True,
        "caller_supplies_coefficient_directional_jets": False,
        "caller_supplies_background_phase_derivative_jet": True,
        "caller_supplies_amplitude_directional_jet": True,
        "source_amplitude_ode_materialized": False,
        "source_support_cutoff_runtime_materialized": False,
        "self_contained_velocity_xyzt_provider": False,
        "current_runtime_phase_source_exact": False,
        "current_runtime_support_source_exact": False,
        "current_runtime_complete_curl_source_equivalence_verified": False,
        "source_to_runtime_parameter_map_complete": False,
        "paper_exact": False,
        "matched_pressure_materialized": False,
        "restricted_forcing_materialized": False,
        "complete_ns_residual_assessed": False,
        "residual_reduction_claimed": False,
        "pde_validated": False,
    }
