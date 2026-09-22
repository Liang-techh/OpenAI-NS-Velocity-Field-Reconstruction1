"""Assemble one corrected-source localized complete-curl harmonic from zero data.

This is a narrow Kokuno Agent-2 adapter.  It connects three already executable
source-facing pieces without inventing the missing numerical source data:

* the fixed-ladder zero-data amplitude/sensitivity solve, which generates
  ``t_m, D_r t_m, D_z t_m`` from caller-supplied corrected background and
  harmonic forcing jets;
* the corrected phase/coefficient chain
  ``(Phi,n_Phi,t_m) -> (C_m,D_r C_m,D_z C_m)``;
* potential-level localization before the complete cylindrical curl.

For one positive harmonic representative, let ``eta`` be the source slow
partition factor and ``psi(v)`` the source pulse cutoff.  The corrected source
has ``D_r v = D_z v = 0`` on the band chart, hence ``D_r psi = D_z psi = 0``
for this spatial curl.  With ``w = psi eta`` we therefore form

    C_hat = w C_m,
    D_r C_hat = w D_r C_m + psi (D_r eta) C_m,
    D_z C_hat = w D_z C_m + psi (D_z eta) C_m,

and feed those product-rule jets into the *complete* normalized curl.  This is
not the invalid shortcut ``u_hat = w u``: cutoff-gradient curl terms are
retained whenever ``D eta`` is nonzero.

The returned harmonic is still a corrected-source *chart* object.  The
background, forcing and concrete partition/cutoff values are upstream inputs;
the zero-data solve is a repository RK4 realization rather than the source's
exact B/V_m Duhamel propagator; physical Q-scaling and the project-domain
Cartesian map are not applied.  Accordingly this module does not claim a
self-contained ``velocity(x,y,z,t)`` provider, source/runtime equivalence,
paper exactness, a Navier--Stokes residual, or PDE validation.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np

from .kokuno_corrected_oscillation_source_ledger import (
    SOURCE_COMMIT,
    SOURCE_REPOSITORY,
    SOURCE_WORKBENCH_BLOB,
    default_kokuno_corrected_oscillation_source_ledger,
)
from .kokuno_source_normalized_complete_curl_reference import (
    CompleteHarmonicAmplitude,
    complete_harmonic_amplitude_from_coefficient_jet,
)
from .kokuno_source_phase_coefficient_jet import (
    SourceBackgroundPhaseJet,
    SourceCoefficientDirectionalJet,
    SourcePhaseVectorJet,
    harmonic_coefficient_directional_jet,
    source_phase_vector_jet,
)
from .kokuno_source_zero_data_amplitude_sensitivity import (
    ModeForcingDirectionalJet,
    SourceZeroDataAmplitudeSensitivityResult,
    solve_source_zero_data_amplitude_sensitivity,
    source_zero_data_amplitude_sensitivity_contract,
)


TASK = "K2-OSC-098"
PARENT_A2_HEAD = "c40d8ddecd2971544a6e07dab093436b423cf326"
AMPLITUDE_SENSITIVITY_BLOB = "eba3c00703bd5763b12e6f48a8c34eb918dec82e"
PHASE_COEFFICIENT_BLOB = "d3c9684d47cec5dfabf713db01ffab45ea9461c4"
COMPLETE_CURL_REFERENCE_BLOB = "42809fc1b10937b4246c65121ac1dbfdc2986822"
SOURCE_LOCALIZER_BLOB = "831d498c957493ea8fca197ed45738385935e655"
PARTITION_ATOL = 2.0e-12
COEFFICIENT_CLOSURE_ATOL = 2.0e-12


@dataclass(frozen=True)
class SourceLocalizationJet:
    """Caller-supplied source slow support data for one band/harmonic.

    ``eta``, ``dr_eta`` and ``dz_eta`` use the source normalized directional
    convention.  ``pulse_cutoff`` is ``psi(v)`` at the already-fixed pulse
    coordinate.  The source band geometry has ``D_r v=D_z v=0``, so no spatial
    derivatives of ``pulse_cutoff`` are requested here.
    """

    eta: Any
    dr_eta: Any
    dz_eta: Any
    pulse_cutoff: float


@dataclass(frozen=True)
class SourceLocalizedZeroDataHarmonic:
    """Executable localized harmonic and its identity-preserving ingredients."""

    amplitude_sensitivity: SourceZeroDataAmplitudeSensitivityResult
    phase_jet: SourcePhaseVectorJet
    coefficient_jet: SourceCoefficientDirectionalJet
    complete_amplitude: CompleteHarmonicAmplitude
    localization_weight: np.ndarray
    dr_localization_weight: np.ndarray
    dz_localization_weight: np.ndarray
    localized_coefficient: np.ndarray
    localized_dr_coefficient: np.ndarray
    localized_dz_coefficient: np.ndarray
    complex_vector_potential: np.ndarray
    complex_velocity_cylindrical: np.ndarray
    real_pair_velocity_cylindrical: np.ndarray
    real_pair_velocity_chart_cartesian: np.ndarray


def _finite_real(value: Any, *, name: str) -> np.ndarray:
    array = np.asarray(value, dtype=float)
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must be finite")
    return array


def _positive_harmonic(value: int) -> int:
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)):
        raise ValueError("m must be a positive integer representative of a conjugate pair")
    m = int(value)
    if m <= 0:
        raise ValueError("m must be positive when materializing the real conjugate pair")
    return m


def _positive_scalar(value: float, *, name: str) -> float:
    value = float(value)
    if not np.isfinite(value) or value <= 0.0:
        raise ValueError(f"{name} must be finite and positive")
    return value


def _localization_arrays(
    localization: SourceLocalizationJet,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    if not isinstance(localization, SourceLocalizationJet):
        raise ValueError("localization must be SourceLocalizationJet")
    eta = _finite_real(localization.eta, name="localization.eta")
    dr_eta = _finite_real(localization.dr_eta, name="localization.dr_eta")
    dz_eta = _finite_real(localization.dz_eta, name="localization.dz_eta")
    eta, dr_eta, dz_eta = np.broadcast_arrays(eta, dr_eta, dz_eta)
    if np.any(np.abs(eta) > 1.0 + PARTITION_ATOL):
        raise ValueError("eta is incompatible with a real squared partition: |eta| must be <=1")
    psi = float(localization.pulse_cutoff)
    if not np.isfinite(psi) or psi < 0.0 or psi > 1.0:
        raise ValueError("pulse_cutoff must be finite and lie in [0,1]")
    return eta, dr_eta, dz_eta, psi


def _chart_cartesian_from_cylindrical(theta: Any, vector: np.ndarray) -> np.ndarray:
    theta = _finite_real(theta, name="theta")
    vector = np.asarray(vector, dtype=float)
    if vector.ndim == 0 or vector.shape[-1] != 3 or not np.all(np.isfinite(vector)):
        raise ValueError("cylindrical chart velocity must be a finite 3-vector field")
    shape = np.broadcast_shapes(theta.shape, vector.shape[:-1])
    theta = np.broadcast_to(theta, shape)
    vector = np.broadcast_to(vector, shape + (3,))
    c = np.cos(theta)
    s = np.sin(theta)
    u_r = vector[..., 0]
    u_theta = vector[..., 1]
    u_z = vector[..., 2]
    return np.stack(
        (u_r * c - u_theta * s, u_r * s + u_theta * c, u_z),
        axis=-1,
    )


def materialize_source_zero_data_localized_harmonic(
    pulse_v: float,
    radius: Any,
    theta: Any,
    z_normalized: Any,
    background: SourceBackgroundPhaseJet,
    forcing_jet_m: ModeForcingDirectionalJet,
    localization: SourceLocalizationJet,
    *,
    epsilon: float,
    p: float,
    p_z: float,
    x_0: float,
    k: float,
    m: int,
) -> SourceLocalizedZeroDataHarmonic:
    """Generate, localize, completely curl and real-pair one source harmonic.

    No caller controls the parent RK4 panel ladder or a numerical tolerance.
    ``m`` is required to be positive because this surface emits the physical
    real contribution of the ``+/-m`` conjugate pair as ``2*Re(u_m)``.
    """

    m = _positive_harmonic(m)
    k = _positive_scalar(k, name="k")
    eta, dr_eta, dz_eta, psi = _localization_arrays(localization)

    sensitivity = solve_source_zero_data_amplitude_sensitivity(
        pulse_v,
        radius,
        theta,
        z_normalized,
        background,
        forcing_jet_m,
        epsilon=epsilon,
        p=p,
        p_z=p_z,
        x_0=x_0,
        k=k,
        m=m,
    )
    phase = source_phase_vector_jet(
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
    coefficient = harmonic_coefficient_directional_jet(
        phase,
        sensitivity.amplitude_jet,
        k=k,
        m=m,
    )

    scalar_shape = np.broadcast_shapes(
        np.shape(phase.phase),
        np.shape(coefficient.coefficient)[:-1],
        eta.shape,
    )
    eta = np.broadcast_to(eta, scalar_shape)
    dr_eta = np.broadcast_to(dr_eta, scalar_shape)
    dz_eta = np.broadcast_to(dz_eta, scalar_shape)
    c = np.broadcast_to(coefficient.coefficient, scalar_shape + (3,))
    dr_c = np.broadcast_to(coefficient.dr_coefficient, scalar_shape + (3,))
    dz_c = np.broadcast_to(coefficient.dz_coefficient, scalar_shape + (3,))
    t_m = np.broadcast_to(
        np.asarray(sensitivity.amplitude_jet.t_m, dtype=np.complex128),
        scalar_shape + (3,),
    )
    n_phi = np.broadcast_to(
        np.asarray(phase.n_phi, dtype=float),
        scalar_shape + (3,),
    )
    radius_b = np.broadcast_to(_finite_real(radius, name="radius"), scalar_shape)
    phase_b = np.broadcast_to(_finite_real(phase.phase, name="phase"), scalar_shape)

    weight = psi * eta
    dr_weight = psi * dr_eta
    dz_weight = psi * dz_eta
    c_local = weight[..., None] * c
    dr_c_local = weight[..., None] * dr_c + dr_weight[..., None] * c
    dz_c_local = weight[..., None] * dz_c + dz_weight[..., None] * c
    t_local = weight[..., None] * t_m

    complete = complete_harmonic_amplitude_from_coefficient_jet(
        radius_b,
        n_phi,
        t_local,
        dr_c_local,
        dz_c_local,
        k=k,
        m=m,
    )
    if not np.allclose(
        complete.coefficient,
        c_local,
        rtol=0.0,
        atol=COEFFICIENT_CLOSURE_ATOL,
    ):
        raise RuntimeError("localized coefficient drifted from potential-level product rule")

    exponential = np.exp(1j * (k * float(m)) * phase_b)
    complex_potential = c_local * exponential[..., None]
    complex_velocity = complete.amplitude * exponential[..., None]
    real_cylindrical = 2.0 * np.real(complex_velocity)
    real_chart_cartesian = _chart_cartesian_from_cylindrical(theta, real_cylindrical)

    return SourceLocalizedZeroDataHarmonic(
        amplitude_sensitivity=sensitivity,
        phase_jet=phase,
        coefficient_jet=coefficient,
        complete_amplitude=complete,
        localization_weight=weight,
        dr_localization_weight=dr_weight,
        dz_localization_weight=dz_weight,
        localized_coefficient=c_local,
        localized_dr_coefficient=dr_c_local,
        localized_dz_coefficient=dz_c_local,
        complex_vector_potential=complex_potential,
        complex_velocity_cylindrical=complex_velocity,
        real_pair_velocity_cylindrical=real_cylindrical,
        real_pair_velocity_chart_cartesian=real_chart_cartesian,
    )


def source_zero_data_localized_harmonic_contract() -> dict[str, Any]:
    """Return the source/realization boundary for this adapter."""

    ledger = default_kokuno_corrected_oscillation_source_ledger()
    sensitivity = source_zero_data_amplitude_sensitivity_contract()
    if ledger.harmonic_vector_potential != "A_m = C_m exp(i k m Phi)":
        raise RuntimeError("corrected-source harmonic vector-potential formula drifted")
    if "conjugate harmonic pairs" not in ledger.real_wave_pairing:
        raise RuntimeError("corrected-source real-wave pairing statement drifted")
    if sensitivity["source_amplitude_directional_jet_materialized"] is not True:
        raise RuntimeError("required zero-data amplitude directional jet is unavailable")

    return {
        "schema": "kokuno-source-zero-data-localized-harmonic-v1",
        "task": TASK,
        "parent_a2_head": PARENT_A2_HEAD,
        "amplitude_sensitivity_blob": AMPLITUDE_SENSITIVITY_BLOB,
        "phase_coefficient_blob": PHASE_COEFFICIENT_BLOB,
        "complete_curl_reference_blob": COMPLETE_CURL_REFERENCE_BLOB,
        "source_localizer_blob": SOURCE_LOCALIZER_BLOB,
        "source_repository": SOURCE_REPOSITORY,
        "source_commit": SOURCE_COMMIT,
        "source_workbench_blob": SOURCE_WORKBENCH_BLOB,
        "source_amplitude_directional_jet_consumed": True,
        "potential_level_product_rule_localization_executable": True,
        "pulse_cutoff_spatial_derivatives_zero_from_Dr_v_Dz_v_zero": True,
        "cutoff_gradient_complete_curl_terms_retained": True,
        "real_conjugate_pair_velocity_materialized_in_source_chart": True,
        "source_chart_cylindrical_to_cartesian_basis_rotation_materialized": True,
        "caller_supplies_mode_forcing_directional_jet": True,
        "caller_supplies_corrected_background": True,
        "caller_supplies_source_partition_jet": True,
        "caller_supplies_pulse_cutoff_value": True,
        "concrete_source_partition_bumps_reconstructed": False,
        "source_forcing_provider_materialized": False,
        "actual_corrected_background_provider_materialized": False,
        "source_exact_duhamel_frame_B_materialized": False,
        "source_exact_duhamel_propagator_Vm_materialized": False,
        "source_physical_Q_scaling_applied": False,
        "project_domain_coordinate_map_applied": False,
        "global_axis_safe_source_velocity_materialized": False,
        "self_contained_velocity_xyzt_provider": False,
        "current_autonomous_runtime_source_equivalence_verified": False,
        "complete_ns_residual_assessed": False,
        "residual_reduction_claimed": False,
        "paper_exact": False,
        "pde_validated": False,
    }
