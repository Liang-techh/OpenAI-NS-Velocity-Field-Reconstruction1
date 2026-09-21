from __future__ import annotations

import inspect

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_source_normalized_complete_curl_reference import (
    complete_harmonic_amplitude_from_coefficient_jet,
    harmonic_vector_potential_coefficient,
)
from openai_ns_reconstruction.kokuno_source_phase_coefficient_jet import (
    PARENT_A2_HEAD,
    PARENT_COMPLETE_CURL_REFERENCE_BLOB,
    SourceAmplitudeDirectionalJet,
    SourceBackgroundPhaseJet,
    complete_harmonic_from_source_phase_jet,
    harmonic_coefficient_directional_jet,
    source_phase_coefficient_jet_contract,
    source_phase_vector_jet,
)


EPSILON = 0.2
P = 1.25
P_Z = -0.4
X_0 = 0.8
PULSE_V = 0.35
K = 7.0
M = 1
B = np.array([0.3, -0.7, 1.1])


def _background(radius, z_normalized) -> SourceBackgroundPhaseJet:
    radius = np.asarray(radius, dtype=float)
    z_normalized = np.asarray(z_normalized, dtype=float)
    radius, z_normalized = np.broadcast_arrays(radius, z_normalized)
    return SourceBackgroundPhaseJet(
        f=radius**2 + z_normalized**2,
        g=radius * z_normalized,
        f_r=2.0 * radius,
        g_r=z_normalized,
        f_z=2.0 * z_normalized,
        g_z=radius,
        f_rr=2.0 * np.ones_like(radius),
        g_rr=np.zeros_like(radius),
        f_rz=np.zeros_like(radius),
        g_rz=np.ones_like(radius),
        f_zz=2.0 * np.ones_like(radius),
        g_zz=np.zeros_like(radius),
    )


def _phase(radius, z_normalized, theta=0.37):
    return source_phase_vector_jet(
        radius,
        theta,
        z_normalized,
        PULSE_V,
        _background(radius, z_normalized),
        epsilon=EPSILON,
        p=P,
        p_z=P_Z,
        x_0=X_0,
    )


def _amplitude(phase_jet) -> SourceAmplitudeDirectionalJet:
    t_m = np.cross(phase_jet.n_phi, B)
    dr_t_m = np.cross(phase_jet.dr_n_phi, B)
    dz_t_m = np.cross(phase_jet.dz_n_phi, B)
    return SourceAmplitudeDirectionalJet(t_m=t_m, dr_t_m=dr_t_m, dz_t_m=dz_t_m)


def _coefficient_at(radius: float, z_normalized: float) -> np.ndarray:
    phase_jet = _phase(radius, z_normalized)
    t_m = np.cross(phase_jet.n_phi, B)
    return harmonic_vector_potential_coefficient(phase_jet.n_phi, t_m, k=K, m=M)


def test_source_phase_vector_jet_replays_corrected_source_formulas() -> None:
    radius = np.array([0.73, 1.17, 1.91])
    zed = np.array([-0.31, 0.14, 0.42])
    theta = np.array([-0.7, 0.2, 1.1])
    actual = source_phase_vector_jet(
        radius,
        theta,
        zed,
        PULSE_V,
        _background(radius, zed),
        epsilon=EPSILON,
        p=P,
        p_z=P_Z,
        x_0=X_0,
    )

    f = radius**2 + zed**2
    g = radius * zed
    h = P * f + P_Z * g
    h_r = P * (2.0 * radius) + P_Z * zed
    h_z = P * (2.0 * zed) + P_Z * radius
    h_rr = 2.0 * P
    h_rz = P_Z
    h_zz = 2.0 * P

    expected_phase = P * theta + P_Z * zed / EPSILON + X_0 * radius - PULSE_V * h
    expected_n = np.stack(
        (
            X_0 - PULSE_V * h_r,
            P / radius,
            P_Z - EPSILON * PULSE_V * h_z,
        ),
        axis=-1,
    )
    expected_dr_n = np.stack(
        (
            np.full_like(radius, -PULSE_V * h_rr),
            -P / radius**2,
            np.full_like(radius, -EPSILON * PULSE_V * h_rz),
        ),
        axis=-1,
    )
    expected_dz_n = np.stack(
        (
            np.full_like(radius, -EPSILON * PULSE_V * h_rz),
            np.zeros_like(radius),
            np.full_like(radius, -(EPSILON**2) * PULSE_V * h_zz),
        ),
        axis=-1,
    )

    np.testing.assert_allclose(actual.h_phi, h, rtol=0.0, atol=2.0e-14)
    np.testing.assert_allclose(actual.phase, expected_phase, rtol=0.0, atol=2.0e-14)
    np.testing.assert_allclose(actual.n_phi, expected_n, rtol=0.0, atol=2.0e-14)
    np.testing.assert_allclose(actual.dr_n_phi, expected_dr_n, rtol=0.0, atol=2.0e-14)
    np.testing.assert_allclose(actual.dz_n_phi, expected_dz_n, rtol=0.0, atol=2.0e-14)


def test_coefficient_directional_jets_match_independent_finite_difference() -> None:
    radius = 1.23
    zed = -0.27
    phase_jet = _phase(radius, zed)
    amplitude_jet = _amplitude(phase_jet)
    actual = harmonic_coefficient_directional_jet(phase_jet, amplitude_jet, k=K, m=M)

    h = 1.0e-6
    fd_r = (_coefficient_at(radius + h, zed) - _coefficient_at(radius - h, zed)) / (2.0 * h)
    fd_z = EPSILON * (
        _coefficient_at(radius, zed + h) - _coefficient_at(radius, zed - h)
    ) / (2.0 * h)

    np.testing.assert_allclose(actual.dr_coefficient, fd_r, rtol=2.0e-8, atol=2.0e-9)
    np.testing.assert_allclose(actual.dz_coefficient, fd_z, rtol=2.0e-8, atol=2.0e-9)


def test_complete_harmonic_consumes_derived_coefficient_jets() -> None:
    radius = np.array([0.82, 1.06, 1.44])
    zed = np.array([-0.19, 0.08, 0.33])
    phase_jet = _phase(radius, zed)
    amplitude_jet = _amplitude(phase_jet)

    result = complete_harmonic_from_source_phase_jet(
        radius, phase_jet, amplitude_jet, k=K, m=M
    )
    direct = complete_harmonic_amplitude_from_coefficient_jet(
        radius,
        phase_jet.n_phi,
        amplitude_jet.t_m,
        result.coefficient_jet.dr_coefficient,
        result.coefficient_jet.dz_coefficient,
        k=K,
        m=M,
    )

    np.testing.assert_allclose(
        result.complete_amplitude.coefficient,
        result.coefficient_jet.coefficient,
        rtol=0.0,
        atol=2.0e-14,
    )
    np.testing.assert_allclose(
        result.complete_amplitude.remainder, direct.remainder, rtol=0.0, atol=2.0e-14
    )
    np.testing.assert_allclose(
        result.complete_amplitude.amplitude, direct.amplitude, rtol=0.0, atol=2.0e-14
    )
    np.testing.assert_allclose(
        result.complete_amplitude.longitudinal_amplitude,
        result.complete_amplitude.longitudinal_remainder,
        rtol=0.0,
        atol=2.0e-12,
    )


def test_batch_broadcast_and_source_operator_scaling_are_preserved() -> None:
    radius = np.array([[0.71], [1.21]])
    zed = np.array([[-0.3, 0.0, 0.4]])
    phase_jet = _phase(radius, zed, theta=np.array([[0.1, 0.5, 0.9]]))
    amplitude_jet = _amplitude(phase_jet)
    coefficient_jet = harmonic_coefficient_directional_jet(
        phase_jet, amplitude_jet, k=K, m=M
    )

    assert phase_jet.n_phi.shape == (2, 3, 3)
    assert coefficient_jet.coefficient.shape == (2, 3, 3)
    assert coefficient_jet.dr_coefficient.shape == (2, 3, 3)
    assert coefficient_jet.dz_coefficient.shape == (2, 3, 3)


def test_directional_transversality_drift_fails_closed() -> None:
    phase_jet = _phase(1.1, 0.2)
    good = _amplitude(phase_jet)
    bad = SourceAmplitudeDirectionalJet(
        t_m=good.t_m,
        dr_t_m=np.asarray(good.dr_t_m) + np.asarray(phase_jet.n_phi),
        dz_t_m=good.dz_t_m,
    )
    with pytest.raises(ValueError, match="D_r amplitude jet"):
        harmonic_coefficient_directional_jet(phase_jet, bad, k=K, m=M)


def test_phase_adapter_fails_closed_on_axis_invalid_epsilon_and_nonfinite_inputs() -> None:
    with pytest.raises(ValueError, match="radius > 0"):
        _phase(0.0, 0.1)
    with pytest.raises(ValueError, match="epsilon"):
        source_phase_vector_jet(
            1.0,
            0.2,
            0.1,
            PULSE_V,
            _background(1.0, 0.1),
            epsilon=0.0,
            p=P,
            p_z=P_Z,
            x_0=X_0,
        )
    with pytest.raises(ValueError, match="finite"):
        _phase(1.0, np.nan)


def test_contract_advances_only_coefficient_jet_input_completeness() -> None:
    contract = source_phase_coefficient_jet_contract()
    assert PARENT_A2_HEAD == "ecc2be96eee441308cc0216bf187b43bbd7fd3a0"
    assert PARENT_COMPLETE_CURL_REFERENCE_BLOB == "42809fc1b10937b4246c65121ac1dbfdc2986822"
    assert contract["source_phase_vector_jet_executable"] is True
    assert contract["source_coefficient_directional_jets_derived_analytically"] is True
    assert contract["source_complete_curl_consumes_derived_coefficient_jets"] is True
    assert contract["caller_supplies_coefficient_directional_jets"] is False
    assert contract["caller_supplies_background_phase_derivative_jet"] is True
    assert contract["caller_supplies_amplitude_directional_jet"] is True

    for key in (
        "source_amplitude_ode_materialized",
        "source_support_cutoff_runtime_materialized",
        "self_contained_velocity_xyzt_provider",
        "current_runtime_phase_source_exact",
        "current_runtime_support_source_exact",
        "current_runtime_complete_curl_source_equivalence_verified",
        "source_to_runtime_parameter_map_complete",
        "paper_exact",
        "matched_pressure_materialized",
        "restricted_forcing_materialized",
        "complete_ns_residual_assessed",
        "residual_reduction_claimed",
        "pde_validated",
    ):
        assert contract[key] is False

    for function in (
        source_phase_vector_jet,
        harmonic_coefficient_directional_jet,
        complete_harmonic_from_source_phase_jet,
    ):
        parameters = inspect.signature(function).parameters
        assert "tolerance" not in parameters
        assert "derivative_step" not in parameters
        assert "residual" not in parameters
