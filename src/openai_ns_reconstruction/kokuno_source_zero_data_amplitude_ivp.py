"""Deterministic zero-data numerical IVP for the corrected Kokuno mode equation.

The corrected 2026-09-09 reconstruction records a zero-data Duhamel solution in
its moving transverse frame,

    z_m(v) = integral_0^v V_m(v,w) g_m(w) dw,

with ``t_m = B z_m``.  The repository does not yet materialize the exact source
frame ``B`` or propagator ``V_m``.  This module therefore makes a narrower,
truthful increment: it integrates the already-executable corrected projected
``t_m`` equation from :mod:`kokuno_source_projected_amplitude_ode` with
``t_m(0)=0``.  The integration is a deterministic numerical approximation to
the source IVP, not a paper-exact implementation of the source Duhamel frame.

The source phase vector is regenerated analytically at every Runge--Kutta stage
from the fixed background jet and the corrected phase formula.  Callers still
supply the corrected background data and the mode source ``f_m(v)``.  The panel
ladder is fixed in this module and is intentionally not caller-tunable, so this
reference cannot acquire hidden accuracy/fit knobs.

No ``D_r t_m``/``D_z t_m`` sensitivity equation, source cutoff, global
Cartesian velocity, matched full-candidate pressure/forcing, or Navier--Stokes
residual is claimed here.
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
    SourceBackgroundPhaseJet,
    source_phase_vector_jet,
)
from openai_ns_reconstruction.kokuno_source_projected_amplitude_ode import (
    SourceProjectedAmplitudeODEResult,
    source_projected_amplitude_ode_rhs,
)


PARENT_A2_HEAD = "eee26916e5103ac18ac57422be7f937ec5db6b56"
PARENT_PROJECTED_ODE_BLOB = "ee55f4cb00eee67a11bcf0361aadd1c0c3012441"
IVP_PANEL_LADDER = (32, 64, 128)
SOURCE_ZERO_DATA_DUHAMEL_FORMULA = (
    "z_m(v) = integral_0^v V_m(v,w) g_m(w) dw; t_m = B z_m"
)


ModeForcing = Callable[[float], Any]


@dataclass(frozen=True)
class SourceZeroDataAmplitudeIVPResult:
    """Fixed-resolution-ladder result for one zero-data harmonic IVP."""

    pulse_v: float
    panels: tuple[int, int, int]
    states: tuple[np.ndarray, np.ndarray, np.ndarray]
    state: np.ndarray
    coarse_medium_abs_max: float
    medium_fine_abs_max: float
    endpoint: SourceProjectedAmplitudeODEResult


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


def _abs_max(value: np.ndarray) -> float:
    return float(np.max(np.abs(np.asarray(value))))


def solve_source_zero_data_amplitude_ivp(
    pulse_v: float,
    radius: Any,
    theta: Any,
    z_normalized: Any,
    background: SourceBackgroundPhaseJet,
    forcing_m: ModeForcing,
    *,
    epsilon: float,
    p: float,
    p_z: float,
    x_0: float,
    k: float,
    m: int,
) -> SourceZeroDataAmplitudeIVPResult:
    """Integrate the corrected projected mode equation from ``t_m(0)=0``.

    A classical RK4 solve is repeated on the frozen panel ladder
    ``(32, 64, 128)``.  The ladder is part of the implementation contract and
    cannot be changed by the caller.  ``radius``, the corrected background jet,
    and the other slow/transverse labels are held fixed along the pulse path, as
    in the source Duhamel construction.  The phase vector itself is *not* held
    fixed: it is regenerated at each stage from the corrected ``Phi`` law.

    ``forcing_m(v)`` is the unresolved upstream corrected-source harmonic
    source.  Supplying that source does not make the result a self-contained
    velocity provider.
    """

    pulse_v = _finite_scalar(pulse_v, name="pulse_v")
    epsilon = _finite_positive(epsilon, name="epsilon")
    p = _finite_scalar(p, name="p")
    p_z = _finite_scalar(p_z, name="p_z")
    x_0 = _finite_scalar(x_0, name="x_0")
    k = _finite_positive(k, name="k")
    m = _nonzero_integer(m, name="m")
    if not callable(forcing_m):
        raise ValueError("forcing_m must be callable")

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
    zero_state = np.zeros_like(phase_zero.n_phi, dtype=np.complex128)

    def rhs(v_value: float, state: np.ndarray) -> np.ndarray:
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
        forcing = _complex_vector(forcing_m(float(v_value)), name="forcing_m(v)")
        result = source_projected_amplitude_ode_rhs(
            radius,
            background,
            phase,
            state,
            forcing,
            epsilon=epsilon,
            p=p,
            p_z=p_z,
            k=k,
            m=m,
        )
        return result.derivative

    def integrate(panels: int) -> np.ndarray:
        state = zero_state.copy()
        if pulse_v == 0.0:
            return state
        step = pulse_v / float(panels)
        v_value = 0.0
        for _ in range(panels):
            k1 = rhs(v_value, state)
            k2 = rhs(v_value + 0.5 * step, state + 0.5 * step * k1)
            k3 = rhs(v_value + 0.5 * step, state + 0.5 * step * k2)
            k4 = rhs(v_value + step, state + step * k3)
            state = state + (step / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
            if not np.all(np.isfinite(state.real)) or not np.all(np.isfinite(state.imag)):
                raise FloatingPointError("zero-data source amplitude IVP became non-finite")
            v_value += step
        return state

    states = tuple(integrate(panels) for panels in IVP_PANEL_LADDER)
    coarse_medium = _abs_max(states[1] - states[0])
    medium_fine = _abs_max(states[2] - states[1])
    finest = states[-1]

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
    endpoint_forcing = _complex_vector(forcing_m(pulse_v), name="forcing_m(pulse_v)")
    endpoint = source_projected_amplitude_ode_rhs(
        radius,
        background,
        endpoint_phase,
        finest,
        endpoint_forcing,
        epsilon=epsilon,
        p=p,
        p_z=p_z,
        k=k,
        m=m,
    )

    return SourceZeroDataAmplitudeIVPResult(
        pulse_v=pulse_v,
        panels=IVP_PANEL_LADDER,
        states=(states[0], states[1], states[2]),
        state=finest,
        coarse_medium_abs_max=coarse_medium,
        medium_fine_abs_max=medium_fine,
        endpoint=endpoint,
    )


def source_zero_data_amplitude_ivp_contract() -> dict[str, Any]:
    """Provenance and hard truth boundary for the numerical IVP increment."""

    return {
        "schema": "kokuno-source-zero-data-amplitude-ivp-v1",
        "parent_a2_head": PARENT_A2_HEAD,
        "parent_projected_ode_blob": PARENT_PROJECTED_ODE_BLOB,
        "source_repository": SOURCE_REPOSITORY,
        "source_commit": SOURCE_COMMIT,
        "source_workbench_blob": SOURCE_WORKBENCH_BLOB,
        "source_zero_data_duhamel_formula": SOURCE_ZERO_DATA_DUHAMEL_FORMULA,
        "integration_method": "classical-rk4-fixed-panel-ladder",
        "panel_ladder": list(IVP_PANEL_LADDER),
        "caller_can_tune_integrator_resolution": False,
        "source_phase_jet_regenerated_analytically_along_v": True,
        "caller_supplies_mode_state_t_m": False,
        "caller_supplies_mode_forcing_f_m": True,
        "caller_supplies_corrected_background": True,
        "source_projected_amplitude_rhs_materialized": True,
        "source_zero_data_ivp_numerical_solver_materialized": True,
        "source_exact_duhamel_frame_B_materialized": False,
        "source_exact_duhamel_propagator_Vm_materialized": False,
        "source_amplitude_directional_jet_materialized": False,
        "source_amplitude_mode_provider_complete": False,
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
