"""Executable corrected-source projected harmonic-amplitude ODE algebra.

This module advances the Kokuno Agent-2 source route by turning the corrected
2026-09-09 pulse equation into a batch-safe executable right-hand side.  It
consumes an already-materialized source phase vector and background jet and
implements, for each nonzero harmonic ``m``, the corrected formulas

    t_m' + K t_m + m^2 d t_m + i k m n_Phi pi_m = -f_m,
    n_Phi . t_m = 0,

with

    d = epsilon k^2 |n_Phi|^2,
    pi_m = i/(k m) *
        (n_Phi . K t_m - n_Phi' . t_m + n_Phi . f_m) / |n_Phi|^2.

Eliminating ``pi_m`` gives the projected finite-dimensional ODE and the exact
constraint evolution

    (n_Phi . t_m)' = -m^2 d (n_Phi . t_m).

Hence zero transversality persists.  The pulse derivative ``n_Phi'`` is derived
analytically from the same corrected-source phase law at fixed background
labels; it is not a caller-supplied vector.

This is deliberately only an ODE *RHS/reference surface*.  It does not yet
integrate the zero-data Duhamel/IVP problem, generate ``D_r t_m`` or
``D_z t_m``, materialize the source support/cutoff, or provide a self-contained
Cartesian velocity.  It therefore does not make the repository oscillation
paper-exact and does not assess a Navier--Stokes residual.
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
from openai_ns_reconstruction.kokuno_source_phase_coefficient_jet import (
    SourceBackgroundPhaseJet,
    SourcePhaseVectorJet,
)


PARENT_A2_HEAD = "f3a8ac969832f73aa13d9af260593ea0bd730e46"
PARENT_BACKGROUND_PHASE_JET_BLOB = "c3b455923bf697b15a4c20d5bb6fd678b526d475"


@dataclass(frozen=True)
class SourceProjectedAmplitudeODEResult:
    """One evaluation of the corrected projected harmonic-amplitude ODE."""

    k_matrix: np.ndarray
    damping: np.ndarray
    n_phi_prime: np.ndarray
    pressure: np.ndarray
    derivative: np.ndarray
    full_equation_residual: np.ndarray
    constraint: np.ndarray
    constraint_derivative: np.ndarray
    expected_constraint_derivative: np.ndarray


def _finite(value: Any, *, name: str, dtype: Any = float) -> np.ndarray:
    array = np.asarray(value, dtype=dtype)
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must be finite")
    return array


def _vector(value: Any, *, name: str, dtype: Any) -> np.ndarray:
    array = _finite(value, name=name, dtype=dtype)
    if array.ndim == 0 or array.shape[-1] != 3:
        raise ValueError(f"{name} must have final dimension 3")
    return array


def _positive_radius(radius: Any) -> np.ndarray:
    array = _finite(radius, name="radius")
    if np.any(array <= 0.0):
        raise ValueError("corrected-source amplitude chart requires radius > 0")
    return array


def _positive_scalar(value: float, *, name: str) -> float:
    if not np.isfinite(value) or value <= 0.0:
        raise ValueError(f"{name} must be finite and positive")
    return float(value)


def _finite_scalar(value: float, *, name: str) -> float:
    if not np.isfinite(value):
        raise ValueError(f"{name} must be finite")
    return float(value)


def _nonzero_integer(value: int, *, name: str) -> int:
    if isinstance(value, (bool, np.bool_)):
        raise ValueError(f"{name} must be a nonzero integer")
    integer = int(value)
    if integer != value or integer == 0:
        raise ValueError(f"{name} must be a nonzero integer")
    return integer


def _dot(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    return np.sum(a * b, axis=-1)


def source_amplitude_matrix(
    radius: Any,
    background: SourceBackgroundPhaseJet,
) -> np.ndarray:
    """Build the corrected source matrix ``K`` from ``F, F_R, G_R``.

    ``K = [[0,-2F,0], [2F+R F_R,0,0], [G_R,0,0]]``.
    """

    radius = _positive_radius(radius)
    f = _finite(background.f, name="background.f")
    f_r = _finite(background.f_r, name="background.f_r")
    g_r = _finite(background.g_r, name="background.g_r")
    radius, f, f_r, g_r = np.broadcast_arrays(radius, f, f_r, g_r)
    zero = np.zeros_like(f)
    row0 = np.stack((zero, -2.0 * f, zero), axis=-1)
    row1 = np.stack((2.0 * f + radius * f_r, zero, zero), axis=-1)
    row2 = np.stack((g_r, zero, zero), axis=-1)
    return np.stack((row0, row1, row2), axis=-2)


def source_phase_vector_pulse_derivative(
    background: SourceBackgroundPhaseJet,
    *,
    epsilon: float,
    p: float,
    p_z: float,
) -> np.ndarray:
    """Derive ``n_Phi' = partial_v n_Phi`` at fixed source labels.

    From

    ``n_Phi=(x_0-v H_R, p/R, p_z-epsilon v H_Z)``

    and the source pulse geometry where the background is fixed along ``v``, the
    pulse derivative is ``(-H_R, 0, -epsilon H_Z)``.
    """

    epsilon = _positive_scalar(epsilon, name="epsilon")
    p = _finite_scalar(p, name="p")
    p_z = _finite_scalar(p_z, name="p_z")
    f_r = _finite(background.f_r, name="background.f_r")
    g_r = _finite(background.g_r, name="background.g_r")
    f_z = _finite(background.f_z, name="background.f_z")
    g_z = _finite(background.g_z, name="background.g_z")
    f_r, g_r, f_z, g_z = np.broadcast_arrays(f_r, g_r, f_z, g_z)
    h_r = p * f_r + p_z * g_r
    h_z = p * f_z + p_z * g_z
    return np.stack(
        (-h_r, np.zeros_like(h_r), -epsilon * h_z),
        axis=-1,
    )


def source_projected_amplitude_ode_rhs(
    radius: Any,
    background: SourceBackgroundPhaseJet,
    phase_jet: SourcePhaseVectorJet,
    t_m: Any,
    forcing_m: Any,
    *,
    epsilon: float,
    p: float,
    p_z: float,
    k: float,
    m: int,
) -> SourceProjectedAmplitudeODEResult:
    """Evaluate the corrected projected amplitude equation and source pressure.

    ``t_m`` and ``forcing_m`` may be complex and support NumPy broadcasting.
    The source matrix, damping, pulse derivative of the phase normal, pressure,
    projected ODE derivative, full-equation closure residual, and exact
    constraint-evolution terms are returned together.

    This function intentionally does not treat ``t_m`` as a fitted velocity
    component.  It is the state of the source finite-dimensional pulse ODE; a
    zero-data IVP/Duhamel solver remains a separate unresolved increment.
    """

    radius = _positive_radius(radius)
    epsilon = _positive_scalar(epsilon, name="epsilon")
    p = _finite_scalar(p, name="p")
    p_z = _finite_scalar(p_z, name="p_z")
    k = _positive_scalar(k, name="k")
    m = _nonzero_integer(m, name="m")

    n_phi = _vector(phase_jet.n_phi, name="phase_jet.n_phi", dtype=float)
    t_m = _vector(t_m, name="t_m", dtype=np.complex128)
    forcing_m = _vector(forcing_m, name="forcing_m", dtype=np.complex128)
    k_matrix = source_amplitude_matrix(radius, background)
    n_phi_prime = source_phase_vector_pulse_derivative(
        background,
        epsilon=epsilon,
        p=p,
        p_z=p_z,
    )

    n_phi, t_m, forcing_m, n_phi_prime = np.broadcast_arrays(
        n_phi, t_m, forcing_m, n_phi_prime
    )
    batch_shape = n_phi.shape[:-1]
    k_matrix = np.broadcast_to(k_matrix, batch_shape + (3, 3))

    norm_sq = _dot(n_phi, n_phi)
    if np.any(norm_sq <= 0.0) or not np.all(np.isfinite(norm_sq)):
        raise ValueError("source phase vector n_phi must be nonzero")

    damping = epsilon * (k**2) * norm_sq
    k_t = np.einsum("...ij,...j->...i", k_matrix, t_m)
    n_dot_k_t = _dot(n_phi, k_t)
    n_prime_dot_t = _dot(n_phi_prime, t_m)
    n_dot_forcing = _dot(n_phi, forcing_m)

    pressure = (
        1j
        / float(k * m)
        * (n_dot_k_t - n_prime_dot_t + n_dot_forcing)
        / norm_sq
    )

    derivative = (
        -k_t
        - (float(m * m) * damping)[..., None] * t_m
        - forcing_m
        - 1j * float(k * m) * n_phi * pressure[..., None]
    )

    full_equation_residual = (
        derivative
        + k_t
        + (float(m * m) * damping)[..., None] * t_m
        + 1j * float(k * m) * n_phi * pressure[..., None]
        + forcing_m
    )

    constraint = _dot(n_phi, t_m)
    constraint_derivative = _dot(n_phi_prime, t_m) + _dot(n_phi, derivative)
    expected_constraint_derivative = -float(m * m) * damping * constraint

    return SourceProjectedAmplitudeODEResult(
        k_matrix=k_matrix,
        damping=damping,
        n_phi_prime=n_phi_prime,
        pressure=pressure,
        derivative=derivative,
        full_equation_residual=full_equation_residual,
        constraint=constraint,
        constraint_derivative=constraint_derivative,
        expected_constraint_derivative=expected_constraint_derivative,
    )


def source_projected_amplitude_ode_contract() -> dict[str, Any]:
    """Fail-closed provenance and scientific boundary for this ODE surface."""

    ledger = default_kokuno_corrected_oscillation_source_ledger()
    if ledger.constrained_amplitude_equation != (
        "t_m' + K t_m + m^2 d t_m + i k m n_Phi pi_m = -f_m"
    ):
        raise RuntimeError("corrected-source amplitude equation drifted")
    if ledger.source_pressure_relation != (
        "pi_m = (i/(k m)) (n_Phi dot K t_m - n_Phi' dot t_m + "
        "n_Phi dot f_m) / |n_Phi|^2"
    ):
        raise RuntimeError("corrected-source pressure relation drifted")

    return {
        "schema": "kokuno-source-projected-amplitude-ode-v1",
        "parent_a2_head": PARENT_A2_HEAD,
        "parent_background_phase_jet_blob": PARENT_BACKGROUND_PHASE_JET_BLOB,
        "source_repository": SOURCE_REPOSITORY,
        "source_commit": SOURCE_COMMIT,
        "source_workbench_blob": SOURCE_WORKBENCH_BLOB,
        "source_amplitude_matrix_executable": True,
        "source_phase_vector_pulse_derivative_analytic": True,
        "source_pressure_elimination_executable": True,
        "source_projected_amplitude_rhs_materialized": True,
        "source_constraint_evolution_identity_executable": True,
        "caller_supplies_n_phi_prime": False,
        "caller_supplies_mode_pressure_pi_m": False,
        "caller_supplies_mode_state_t_m": True,
        "caller_supplies_mode_forcing_f_m": True,
        "source_zero_data_ivp_solver_materialized": False,
        "source_amplitude_directional_jet_materialized": False,
        "source_amplitude_ode_materialized": False,
        "actual_corrected_background_provider_materialized": False,
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
