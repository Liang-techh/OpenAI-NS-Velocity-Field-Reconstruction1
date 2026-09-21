"""Directional sensitivities for the corrected Kokuno zero-data amplitude IVP.

This module differentiates the already-executable corrected projected harmonic
amplitude equation used by :mod:`kokuno_source_zero_data_amplitude_ivp` with
respect to the normalized source directions ``D_r`` and ``D_z``.  It therefore
generates

    t_m, D_r t_m, D_z t_m

from one zero-data IVP realization on the same fixed RK4 panel ladder.

The corrected source has
``D_r = partial_R + M d_r R^(d_r-1) L`` and ``D_z = epsilon partial_Z``.
The background/phase coefficients consumed here are the already-materialized
auxiliary-independent specialization, so ``D_r`` on those coefficients reduces
to ``partial_R``.  The caller-supplied ``D_r f_m`` is nevertheless the *full*
normalized source derivative and may contain the auxiliary ``L`` contribution.
Thus this solver does not assume the unresolved source forcing is auxiliary
independent.

The sensitivity equations are a repository-side analytic differentiation of
the corrected projected ODE.  They are not quoted as an additional source
formula, and they do not make the numerical RK4 realization paper-exact.
Callers still provide the corrected background jet and the upstream harmonic
forcing together with its normalized directional derivatives.  No source
forcing provider, source background provider, cutoff/support runtime, global
Cartesian velocity, matched pressure/forcing for the complete NS candidate, or
same-protocol Navier--Stokes residual is materialized here.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable

import numpy as np

from openai_ns_reconstruction.kokuno_corrected_oscillation_source_ledger import (
    SOURCE_COMMIT,
    SOURCE_REPOSITORY,
    SOURCE_WORKBENCH_BLOB,
)
from openai_ns_reconstruction.kokuno_source_phase_coefficient_jet import (
    SourceAmplitudeDirectionalJet,
    SourceBackgroundPhaseJet,
    SourcePhaseVectorJet,
    source_phase_vector_jet,
)
from openai_ns_reconstruction.kokuno_source_projected_amplitude_ode import (
    SourceProjectedAmplitudeODEResult,
    source_amplitude_matrix,
    source_projected_amplitude_ode_rhs,
)
from openai_ns_reconstruction.kokuno_source_zero_data_amplitude_ivp import (
    IVP_PANEL_LADDER,
)


PARENT_A2_HEAD = "b3ac194b410073820d8160a470824034e08e9b5a"
PARENT_ZERO_DATA_IVP_BLOB = "1e13e81bea6c586938c91c610e0812f7283f131d"
PARENT_PROJECTED_ODE_BLOB = "ee55f4cb00eee67a11bcf0361aadd1c0c3012441"


@dataclass(frozen=True)
class SourceModeForcingDirectionalJet:
    """One harmonic source and its normalized ``D_r`` / ``D_z`` derivatives.

    ``dr_value`` is the full normalized source ``D_r f_m``.  If the unresolved
    source has auxiliary dependence, its ``L`` contribution belongs there.
    """

    value: Any
    dr_value: Any
    dz_value: Any


ModeForcingDirectionalJet = Callable[[float], SourceModeForcingDirectionalJet]


@dataclass(frozen=True)
class SourceAmplitudeSensitivityEndpoint:
    """Endpoint mechanics for one zero-data sensitivity realization."""

    base_ode: SourceProjectedAmplitudeODEResult
    dr_pressure: np.ndarray
    dz_pressure: np.ndarray
    dr_constraint: np.ndarray
    dz_constraint: np.ndarray


@dataclass(frozen=True)
class SourceZeroDataAmplitudeSensitivityResult:
    """Fixed-ladder solution for ``t_m, D_r t_m, D_z t_m``."""

    pulse_v: float
    panels: tuple[int, int, int]
    amplitude_jets: tuple[
        SourceAmplitudeDirectionalJet,
        SourceAmplitudeDirectionalJet,
        SourceAmplitudeDirectionalJet,
    ]
    amplitude_jet: SourceAmplitudeDirectionalJet
    t_coarse_medium_abs_max: float
    t_medium_fine_abs_max: float
    dr_coarse_medium_abs_max: float
    dr_medium_fine_abs_max: float
    dz_coarse_medium_abs_max: float
    dz_medium_fine_abs_max: float
    endpoint: SourceAmplitudeSensitivityEndpoint


def _finite_scalar(value: float, *, name: str) -> float:
    value = float(value)
    if not np.isfinite(value):
        raise ValueError(f"{name} must be finite")
    return value


def _finite_positive(value: float, *, name: str) -> float:
    value = _finite_scalar(value, name=name)
    if value <= 0.0:
        raise ValueError(f"{name} must be positive")
    return value


def _nonzero_integer(value: int, *, name: str) -> int:
    if isinstance(value, (bool, np.bool_)):
        raise ValueError(f"{name} must be a nonzero integer")
    integer = int(value)
    if integer != value or integer == 0:
        raise ValueError(f"{name} must be a nonzero integer")
    return integer


def _complex_vector(value: Any, *, name: str) -> np.ndarray:
    array = np.asarray(value, dtype=np.complex128)
    if array.ndim == 0 or array.shape[-1] != 3:
        raise ValueError(f"{name} must have final dimension 3")
    if not np.all(np.isfinite(array.real)) or not np.all(np.isfinite(array.imag)):
        raise ValueError(f"{name} must be finite")
    return array


def _real_vector(value: Any, *, name: str) -> np.ndarray:
    array = np.asarray(value, dtype=float)
    if array.ndim == 0 or array.shape[-1] != 3:
        raise ValueError(f"{name} must have final dimension 3")
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must be finite")
    return array


def _abs_max(value: np.ndarray) -> float:
    return float(np.max(np.abs(np.asarray(value))))


def _dot(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    return np.sum(a * b, axis=-1)


def _directional_matrix_derivatives(
    radius: Any,
    background: SourceBackgroundPhaseJet,
    *,
    epsilon: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Return analytic normalized ``D_r K`` and ``D_z K``.

    The corrected source matrix is
    ``K=[[0,-2F,0],[2F+R F_R,0,0],[G_R,0,0]]``.  Its imported background
    fields are auxiliary independent, hence the auxiliary part of ``D_r``
    vanishes on these coefficients while ``D_z=epsilon partial_Z``.
    """

    radius = np.asarray(radius, dtype=float)
    if not np.all(np.isfinite(radius)) or np.any(radius <= 0.0):
        raise ValueError("corrected-source amplitude chart requires radius > 0")

    f_r = np.asarray(background.f_r, dtype=float)
    f_z = np.asarray(background.f_z, dtype=float)
    f_rr = np.asarray(background.f_rr, dtype=float)
    f_rz = np.asarray(background.f_rz, dtype=float)
    g_rr = np.asarray(background.g_rr, dtype=float)
    g_rz = np.asarray(background.g_rz, dtype=float)
    values = (f_r, f_z, f_rr, f_rz, g_rr, g_rz)
    if not all(np.all(np.isfinite(value)) for value in values):
        raise ValueError("background directional derivatives must be finite")
    radius, f_r, f_z, f_rr, f_rz, g_rr, g_rz = np.broadcast_arrays(
        radius, f_r, f_z, f_rr, f_rz, g_rr, g_rz
    )
    zero = np.zeros_like(radius)
    dr_row0 = np.stack((zero, -2.0 * f_r, zero), axis=-1)
    dr_row1 = np.stack((3.0 * f_r + radius * f_rr, zero, zero), axis=-1)
    dr_row2 = np.stack((g_rr, zero, zero), axis=-1)
    dr_k = np.stack((dr_row0, dr_row1, dr_row2), axis=-2)

    dz_row0 = np.stack((zero, -2.0 * epsilon * f_z, zero), axis=-1)
    dz_row1 = np.stack(
        (epsilon * (2.0 * f_z + radius * f_rz), zero, zero),
        axis=-1,
    )
    dz_row2 = np.stack((epsilon * g_rz, zero, zero), axis=-1)
    dz_k = np.stack((dz_row0, dz_row1, dz_row2), axis=-2)
    return dr_k, dz_k


def _directional_pulse_normal_derivatives(
    background: SourceBackgroundPhaseJet,
    *,
    epsilon: float,
    p: float,
    p_z: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Return ``D_r n_Phi'`` and ``D_z n_Phi'`` analytically."""

    f_rr = np.asarray(background.f_rr, dtype=float)
    g_rr = np.asarray(background.g_rr, dtype=float)
    f_rz = np.asarray(background.f_rz, dtype=float)
    g_rz = np.asarray(background.g_rz, dtype=float)
    f_zz = np.asarray(background.f_zz, dtype=float)
    g_zz = np.asarray(background.g_zz, dtype=float)
    values = (f_rr, g_rr, f_rz, g_rz, f_zz, g_zz)
    if not all(np.all(np.isfinite(value)) for value in values):
        raise ValueError("background second derivatives must be finite")
    f_rr, g_rr, f_rz, g_rz, f_zz, g_zz = np.broadcast_arrays(*values)
    h_rr = p * f_rr + p_z * g_rr
    h_rz = p * f_rz + p_z * g_rz
    h_zz = p * f_zz + p_z * g_zz
    zero = np.zeros_like(h_rr)
    dr_n_prime = np.stack((-h_rr, zero, -epsilon * h_rz), axis=-1)
    dz_n_prime = np.stack(
        (-epsilon * h_rz, zero, -(epsilon**2) * h_zz),
        axis=-1,
    )
    return dr_n_prime, dz_n_prime


def _directional_rhs(
    *,
    radius: Any,
    background: SourceBackgroundPhaseJet,
    phase: SourcePhaseVectorJet,
    base: SourceProjectedAmplitudeODEResult,
    state: np.ndarray,
    sensitivity: np.ndarray,
    forcing: np.ndarray,
    forcing_directional: np.ndarray,
    dn: np.ndarray,
    dk_matrix: np.ndarray,
    dn_prime: np.ndarray,
    epsilon: float,
    k: float,
    m: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Differentiate the corrected projected ODE in one normalized direction."""

    n = _real_vector(phase.n_phi, name="phase.n_phi")
    n_prime = _real_vector(base.n_phi_prime, name="base.n_phi_prime")
    dn = _real_vector(dn, name="dn")
    dn_prime = _real_vector(dn_prime, name="dn_prime")
    state = _complex_vector(state, name="state")
    sensitivity = _complex_vector(sensitivity, name="sensitivity")
    forcing = _complex_vector(forcing, name="forcing")
    forcing_directional = _complex_vector(
        forcing_directional, name="forcing_directional"
    )

    n, n_prime, dn, dn_prime, state, sensitivity, forcing, forcing_directional = (
        np.broadcast_arrays(
            n,
            n_prime,
            dn,
            dn_prime,
            state,
            sensitivity,
            forcing,
            forcing_directional,
        )
    )
    batch_shape = n.shape[:-1]
    k_matrix = np.broadcast_to(
        source_amplitude_matrix(radius, background), batch_shape + (3, 3)
    )
    dk_matrix = np.broadcast_to(dk_matrix, batch_shape + (3, 3))

    norm_sq = _dot(n, n)
    if np.any(norm_sq <= 0.0):
        raise ValueError("source phase vector n_phi must be nonzero")
    d_norm_sq = 2.0 * _dot(n, dn)

    kt = np.einsum("...ij,...j->...i", k_matrix, state)
    d_kt = (
        np.einsum("...ij,...j->...i", dk_matrix, state)
        + np.einsum("...ij,...j->...i", k_matrix, sensitivity)
    )

    numerator = _dot(n, kt) - _dot(n_prime, state) + _dot(n, forcing)
    d_numerator = (
        _dot(dn, kt)
        + _dot(n, d_kt)
        - _dot(dn_prime, state)
        - _dot(n_prime, sensitivity)
        + _dot(dn, forcing)
        + _dot(n, forcing_directional)
    )
    q = numerator / norm_sq
    dq = d_numerator / norm_sq - numerator * d_norm_sq / (norm_sq**2)

    damping = epsilon * (k**2) * norm_sq
    d_damping = epsilon * (k**2) * d_norm_sq
    mu = float(m * m) * damping
    d_mu = float(m * m) * d_damping

    derivative = (
        -d_kt
        - d_mu[..., None] * state
        - mu[..., None] * sensitivity
        - forcing_directional
        + dn * q[..., None]
        + n * dq[..., None]
    )
    pressure_directional = 1j / float(k * m) * dq
    return derivative, pressure_directional


def solve_source_zero_data_amplitude_sensitivity(
    pulse_v: float,
    radius: Any,
    theta: Any,
    z_normalized: Any,
    background: SourceBackgroundPhaseJet,
    forcing_jet_m: ModeForcingDirectionalJet,
    *,
    epsilon: float,
    p: float,
    p_z: float,
    x_0: float,
    k: float,
    m: int,
) -> SourceZeroDataAmplitudeSensitivityResult:
    """Integrate ``t_m, D_r t_m, D_z t_m`` from common zero data.

    The coupled state is advanced by classical RK4 on the same immutable
    ``(32,64,128)`` panel ladder used by the parent zero-data IVP.  The upstream
    callback supplies the source harmonic forcing and its full normalized
    directional derivatives; it does *not* supply amplitude directional jets.
    """

    pulse_v = _finite_scalar(pulse_v, name="pulse_v")
    epsilon = _finite_positive(epsilon, name="epsilon")
    p = _finite_scalar(p, name="p")
    p_z = _finite_scalar(p_z, name="p_z")
    x_0 = _finite_scalar(x_0, name="x_0")
    k = _finite_positive(k, name="k")
    m = _nonzero_integer(m, name="m")
    if not callable(forcing_jet_m):
        raise ValueError("forcing_jet_m must be callable")

    phase_zero = source_phase_vector_jet(
        radius,
        theta,
        z_normalized,
        0.0,
        background,
        epsilon=epsilon,
        p=p,
        p_z=p_z,
        x_0=x_0,
    )
    zero = np.zeros_like(phase_zero.n_phi, dtype=np.complex128)
    dr_k, dz_k = _directional_matrix_derivatives(
        radius, background, epsilon=epsilon
    )
    dr_n_prime, dz_n_prime = _directional_pulse_normal_derivatives(
        background, epsilon=epsilon, p=p, p_z=p_z
    )

    def coerce_forcing(v_value: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        supplied = forcing_jet_m(float(v_value))
        if not isinstance(supplied, SourceModeForcingDirectionalJet):
            raise ValueError(
                "forcing_jet_m(v) must return SourceModeForcingDirectionalJet"
            )
        return (
            _complex_vector(supplied.value, name="forcing_jet_m(v).value"),
            _complex_vector(supplied.dr_value, name="forcing_jet_m(v).dr_value"),
            _complex_vector(supplied.dz_value, name="forcing_jet_m(v).dz_value"),
        )

    def rhs(
        v_value: float,
        state: tuple[np.ndarray, np.ndarray, np.ndarray],
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        t_m, dr_t_m, dz_t_m = state
        phase = source_phase_vector_jet(
            radius,
            theta,
            z_normalized,
            v_value,
            background,
            epsilon=epsilon,
            p=p,
            p_z=p_z,
            x_0=x_0,
        )
        forcing, dr_forcing, dz_forcing = coerce_forcing(v_value)
        base = source_projected_amplitude_ode_rhs(
            radius,
            background,
            phase,
            t_m,
            forcing,
            epsilon=epsilon,
            p=p,
            p_z=p_z,
            k=k,
            m=m,
        )
        dr_rhs, _ = _directional_rhs(
            radius=radius,
            background=background,
            phase=phase,
            base=base,
            state=t_m,
            sensitivity=dr_t_m,
            forcing=forcing,
            forcing_directional=dr_forcing,
            dn=phase.dr_n_phi,
            dk_matrix=dr_k,
            dn_prime=dr_n_prime,
            epsilon=epsilon,
            k=k,
            m=m,
        )
        dz_rhs, _ = _directional_rhs(
            radius=radius,
            background=background,
            phase=phase,
            base=base,
            state=t_m,
            sensitivity=dz_t_m,
            forcing=forcing,
            forcing_directional=dz_forcing,
            dn=phase.dz_n_phi,
            dk_matrix=dz_k,
            dn_prime=dz_n_prime,
            epsilon=epsilon,
            k=k,
            m=m,
        )
        return base.derivative, dr_rhs, dz_rhs

    def add_scaled(
        state: tuple[np.ndarray, np.ndarray, np.ndarray],
        derivative: tuple[np.ndarray, np.ndarray, np.ndarray],
        scale: float,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        return tuple(
            value + scale * increment
            for value, increment in zip(state, derivative, strict=True)
        )  # type: ignore[return-value]

    def integrate(panels: int) -> SourceAmplitudeDirectionalJet:
        state = (zero.copy(), zero.copy(), zero.copy())
        if pulse_v == 0.0:
            return SourceAmplitudeDirectionalJet(
                t_m=state[0], dr_t_m=state[1], dz_t_m=state[2]
            )
        step = pulse_v / float(panels)
        v_value = 0.0
        for _ in range(panels):
            k1 = rhs(v_value, state)
            k2 = rhs(v_value + 0.5 * step, add_scaled(state, k1, 0.5 * step))
            k3 = rhs(v_value + 0.5 * step, add_scaled(state, k2, 0.5 * step))
            k4 = rhs(v_value + step, add_scaled(state, k3, step))
            state = tuple(
                value
                + (step / 6.0)
                * (d1 + 2.0 * d2 + 2.0 * d3 + d4)
                for value, d1, d2, d3, d4 in zip(
                    state, k1, k2, k3, k4, strict=True
                )
            )  # type: ignore[assignment]
            if any(
                not np.all(np.isfinite(value.real))
                or not np.all(np.isfinite(value.imag))
                for value in state
            ):
                raise FloatingPointError(
                    "zero-data source amplitude sensitivity became non-finite"
                )
            v_value += step
        return SourceAmplitudeDirectionalJet(
            t_m=state[0], dr_t_m=state[1], dz_t_m=state[2]
        )

    jets = tuple(integrate(panels) for panels in IVP_PANEL_LADDER)
    finest = jets[-1]
    endpoint_phase = source_phase_vector_jet(
        radius,
        theta,
        z_normalized,
        pulse_v,
        background,
        epsilon=epsilon,
        p=p,
        p_z=p_z,
        x_0=x_0,
    )
    endpoint_forcing, endpoint_dr_forcing, endpoint_dz_forcing = coerce_forcing(
        pulse_v
    )
    endpoint_base = source_projected_amplitude_ode_rhs(
        radius,
        background,
        endpoint_phase,
        finest.t_m,
        endpoint_forcing,
        epsilon=epsilon,
        p=p,
        p_z=p_z,
        k=k,
        m=m,
    )
    _, dr_pressure = _directional_rhs(
        radius=radius,
        background=background,
        phase=endpoint_phase,
        base=endpoint_base,
        state=finest.t_m,
        sensitivity=finest.dr_t_m,
        forcing=endpoint_forcing,
        forcing_directional=endpoint_dr_forcing,
        dn=endpoint_phase.dr_n_phi,
        dk_matrix=dr_k,
        dn_prime=dr_n_prime,
        epsilon=epsilon,
        k=k,
        m=m,
    )
    _, dz_pressure = _directional_rhs(
        radius=radius,
        background=background,
        phase=endpoint_phase,
        base=endpoint_base,
        state=finest.t_m,
        sensitivity=finest.dz_t_m,
        forcing=endpoint_forcing,
        forcing_directional=endpoint_dz_forcing,
        dn=endpoint_phase.dz_n_phi,
        dk_matrix=dz_k,
        dn_prime=dz_n_prime,
        epsilon=epsilon,
        k=k,
        m=m,
    )
    dr_constraint = _dot(endpoint_phase.dr_n_phi, finest.t_m) + _dot(
        endpoint_phase.n_phi, finest.dr_t_m
    )
    dz_constraint = _dot(endpoint_phase.dz_n_phi, finest.t_m) + _dot(
        endpoint_phase.n_phi, finest.dz_t_m
    )

    return SourceZeroDataAmplitudeSensitivityResult(
        pulse_v=pulse_v,
        panels=IVP_PANEL_LADDER,
        amplitude_jets=(jets[0], jets[1], jets[2]),
        amplitude_jet=finest,
        t_coarse_medium_abs_max=_abs_max(jets[1].t_m - jets[0].t_m),
        t_medium_fine_abs_max=_abs_max(jets[2].t_m - jets[1].t_m),
        dr_coarse_medium_abs_max=_abs_max(jets[1].dr_t_m - jets[0].dr_t_m),
        dr_medium_fine_abs_max=_abs_max(jets[2].dr_t_m - jets[1].dr_t_m),
        dz_coarse_medium_abs_max=_abs_max(jets[1].dz_t_m - jets[0].dz_t_m),
        dz_medium_fine_abs_max=_abs_max(jets[2].dz_t_m - jets[1].dz_t_m),
        endpoint=SourceAmplitudeSensitivityEndpoint(
            base_ode=endpoint_base,
            dr_pressure=dr_pressure,
            dz_pressure=dz_pressure,
            dr_constraint=dr_constraint,
            dz_constraint=dz_constraint,
        ),
    )


def source_zero_data_amplitude_sensitivity_contract() -> dict[str, Any]:
    """Provenance and hard truth boundary for the sensitivity increment."""

    return {
        "schema": "kokuno-source-zero-data-amplitude-sensitivity-v1",
        "parent_a2_head": PARENT_A2_HEAD,
        "parent_zero_data_ivp_blob": PARENT_ZERO_DATA_IVP_BLOB,
        "parent_projected_ode_blob": PARENT_PROJECTED_ODE_BLOB,
        "source_repository": SOURCE_REPOSITORY,
        "source_commit": SOURCE_COMMIT,
        "source_workbench_blob": SOURCE_WORKBENCH_BLOB,
        "sensitivity_derivation": (
            "repository analytic differentiation of corrected projected amplitude ODE"
        ),
        "source_normalized_dr_operator_has_auxiliary_L_term": True,
        "background_phase_coefficients_auxiliary_independent": True,
        "forcing_dr_input_is_full_normalized_source_Dr": True,
        "panel_ladder": list(IVP_PANEL_LADDER),
        "caller_can_tune_integrator_resolution": False,
        "caller_supplies_mode_state_t_m": False,
        "caller_supplies_amplitude_directional_jet": False,
        "caller_supplies_mode_forcing_directional_jet": True,
        "caller_supplies_corrected_background": True,
        "source_zero_data_ivp_numerical_solver_materialized": True,
        "source_amplitude_directional_jet_materialized": True,
        "source_complete_curl_directional_input_compatible": True,
        "source_forcing_provider_materialized": False,
        "full_source_forcing_auxiliary_provider_materialized": False,
        "actual_corrected_background_provider_materialized": False,
        "source_exact_duhamel_frame_B_materialized": False,
        "source_exact_duhamel_propagator_Vm_materialized": False,
        "source_amplitude_mode_provider_complete": False,
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
