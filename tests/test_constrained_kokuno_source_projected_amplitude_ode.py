from __future__ import annotations

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_source_phase_coefficient_jet import (
    SourceBackgroundPhaseJet,
    source_phase_vector_jet,
)
from openai_ns_reconstruction.kokuno_source_projected_amplitude_ode import (
    source_amplitude_matrix,
    source_phase_vector_pulse_derivative,
    source_projected_amplitude_ode_contract,
    source_projected_amplitude_ode_rhs,
)


def _background(radius: np.ndarray, zed: np.ndarray) -> SourceBackgroundPhaseJet:
    # F=R^2+Z^2, G=R Z gives a nontrivial source K and phase normal.
    return SourceBackgroundPhaseJet(
        f=radius**2 + zed**2,
        g=radius * zed,
        f_r=2.0 * radius,
        g_r=zed,
        f_z=2.0 * zed,
        g_z=radius,
        f_rr=2.0 * np.ones_like(radius),
        g_rr=np.zeros_like(radius),
        f_rz=np.zeros_like(radius),
        g_rz=np.ones_like(radius),
        f_zz=2.0 * np.ones_like(radius),
        g_zz=np.zeros_like(radius),
    )


def _phase(
    radius: np.ndarray,
    zed: np.ndarray,
    pulse_v: np.ndarray,
    *,
    epsilon: float,
    p: float,
    p_z: float,
):
    return source_phase_vector_jet(
        radius,
        np.array([0.2, -0.6, 0.4]),
        zed,
        pulse_v,
        _background(radius, zed),
        epsilon=epsilon,
        p=p,
        p_z=p_z,
        x_0=0.41,
    )


def _transverse_state(n_phi: np.ndarray) -> np.ndarray:
    seed = np.broadcast_to(
        np.array([0.31 + 0.17j, -0.22 + 0.09j, 0.43 - 0.11j]),
        n_phi.shape,
    )
    return np.cross(n_phi, seed)


def test_source_k_matrix_matches_corrected_formula() -> None:
    radius = np.array([0.83, 1.27, 1.61])
    zed = np.array([-0.31, 0.44, 0.12])
    bg = _background(radius, zed)
    got = source_amplitude_matrix(radius, bg)

    expected = np.zeros((3, 3, 3))
    expected[:, 0, 1] = -2.0 * bg.f
    expected[:, 1, 0] = 2.0 * bg.f + radius * bg.f_r
    expected[:, 2, 0] = bg.g_r
    np.testing.assert_allclose(got, expected, rtol=0.0, atol=2e-14)


def test_phase_pulse_derivative_is_analytic_and_matches_independent_fd() -> None:
    radius = np.array([0.83, 1.27, 1.61])
    zed = np.array([-0.31, 0.44, 0.12])
    pulse_v = np.array([0.35, 0.58, 0.41])
    epsilon = 0.17
    p = 2.0
    p_z = -0.75
    h = 2.0e-6
    bg = _background(radius, zed)

    got = source_phase_vector_pulse_derivative(
        bg, epsilon=epsilon, p=p, p_z=p_z
    )
    plus = _phase(
        radius, zed, pulse_v + h, epsilon=epsilon, p=p, p_z=p_z
    ).n_phi
    minus = _phase(
        radius, zed, pulse_v - h, epsilon=epsilon, p=p, p_z=p_z
    ).n_phi
    fd = (plus - minus) / (2.0 * h)
    np.testing.assert_allclose(got, fd, rtol=0.0, atol=4e-10)


def test_projected_rhs_matches_independent_projection_and_full_equation() -> None:
    radius = np.array([0.83, 1.27, 1.61])
    zed = np.array([-0.31, 0.44, 0.12])
    pulse_v = np.array([0.35, 0.58, 0.41])
    epsilon = 0.17
    p = 2.0
    p_z = -0.75
    k = 5.0
    m = -2
    bg = _background(radius, zed)
    phase = _phase(radius, zed, pulse_v, epsilon=epsilon, p=p, p_z=p_z)
    n = phase.n_phi
    t_m = _transverse_state(n)
    forcing = np.broadcast_to(
        np.array([0.07 - 0.03j, -0.11 + 0.05j, 0.09 + 0.02j]),
        n.shape,
    ).copy()

    got = source_projected_amplitude_ode_rhs(
        radius,
        bg,
        phase,
        t_m,
        forcing,
        epsilon=epsilon,
        p=p,
        p_z=p_z,
        k=k,
        m=m,
    )

    norm_sq = np.sum(n * n, axis=-1)
    k_t = np.einsum("...ij,...j->...i", got.k_matrix, t_m)
    n_dot_k_t = np.sum(n * k_t, axis=-1)
    nprime_dot_t = np.sum(got.n_phi_prime * t_m, axis=-1)
    n_dot_f = np.sum(n * forcing, axis=-1)
    projection_f = forcing - n * (n_dot_f / norm_sq)[..., None]
    a_t = -k_t + n * ((n_dot_k_t - nprime_dot_t) / norm_sq)[..., None]
    expected_rhs = (
        a_t - (m * m * got.damping)[..., None] * t_m - projection_f
    )
    expected_pressure = (
        1j / (k * m) * (n_dot_k_t - nprime_dot_t + n_dot_f) / norm_sq
    )

    np.testing.assert_allclose(got.derivative, expected_rhs, rtol=0.0, atol=3e-13)
    np.testing.assert_allclose(got.pressure, expected_pressure, rtol=0.0, atol=3e-13)
    assert float(np.max(np.abs(got.full_equation_residual))) <= 3e-13


def test_constraint_evolution_keeps_zero_transversality_zero() -> None:
    radius = np.array([0.91, 1.33])
    zed = np.array([0.21, -0.37])
    pulse_v = np.array([0.28, 0.63])
    epsilon = 0.13
    p = -3.0
    p_z = 0.62
    bg = _background(radius, zed)
    phase = source_phase_vector_jet(
        radius,
        np.array([0.1, -0.2]),
        zed,
        pulse_v,
        bg,
        epsilon=epsilon,
        p=p,
        p_z=p_z,
        x_0=-0.27,
    )
    t_m = _transverse_state(phase.n_phi)
    forcing = np.broadcast_to(
        np.array([0.04 + 0.01j, -0.03 + 0.02j, 0.02 - 0.05j]),
        phase.n_phi.shape,
    ).copy()

    got = source_projected_amplitude_ode_rhs(
        radius,
        bg,
        phase,
        t_m,
        forcing,
        epsilon=epsilon,
        p=p,
        p_z=p_z,
        k=7.0,
        m=1,
    )

    assert float(np.max(np.abs(got.constraint))) <= 2e-15
    np.testing.assert_allclose(
        got.constraint_derivative,
        got.expected_constraint_derivative,
        rtol=0.0,
        atol=4e-13,
    )
    assert float(np.max(np.abs(got.constraint_derivative))) <= 4e-13


def test_constraint_decay_identity_also_holds_off_constraint_as_algebra_check() -> None:
    radius = np.asarray(1.07)
    zed = np.asarray(-0.24)
    pulse_v = np.asarray(0.46)
    epsilon = 0.19
    p = 2.0
    p_z = 0.3
    bg = _background(radius, zed)
    phase = source_phase_vector_jet(
        radius,
        0.2,
        zed,
        pulse_v,
        bg,
        epsilon=epsilon,
        p=p,
        p_z=p_z,
        x_0=0.18,
    )
    t_m = np.array([0.2 + 0.1j, -0.4 + 0.05j, 0.3 - 0.2j])
    forcing = np.array([0.01, -0.02j, 0.03 + 0.01j])

    got = source_projected_amplitude_ode_rhs(
        radius,
        bg,
        phase,
        t_m,
        forcing,
        epsilon=epsilon,
        p=p,
        p_z=p_z,
        k=4.0,
        m=3,
    )
    assert abs(complex(got.constraint)) > 1e-6
    np.testing.assert_allclose(
        got.constraint_derivative,
        got.expected_constraint_derivative,
        rtol=0.0,
        atol=3e-13,
    )


@pytest.mark.parametrize(
    ("radius", "epsilon", "k", "m"),
    [
        (0.0, 0.17, 5.0, 1),
        (-1.0, 0.17, 5.0, 1),
        (1.0, 0.0, 5.0, 1),
        (1.0, 0.17, 0.0, 1),
        (1.0, 0.17, 5.0, 0),
    ],
)
def test_projected_rhs_fails_closed_on_invalid_source_domain(
    radius: float, epsilon: float, k: float, m: int
) -> None:
    rr = np.asarray(1.0)
    zz = np.asarray(0.2)
    bg = _background(rr, zz)
    phase = source_phase_vector_jet(
        rr,
        0.1,
        zz,
        0.3,
        bg,
        epsilon=0.17,
        p=2.0,
        p_z=-0.5,
        x_0=0.4,
    )
    with pytest.raises(ValueError):
        source_projected_amplitude_ode_rhs(
            radius,
            bg,
            phase,
            np.array([0.1, 0.2, 0.3]),
            np.array([0.01, 0.02, 0.03]),
            epsilon=epsilon,
            p=2.0,
            p_z=-0.5,
            k=k,
            m=m,
        )


def test_contract_keeps_ivp_directional_jet_and_velocity_claims_false() -> None:
    contract = source_projected_amplitude_ode_contract()
    assert contract["source_amplitude_matrix_executable"] is True
    assert contract["source_phase_vector_pulse_derivative_analytic"] is True
    assert contract["source_pressure_elimination_executable"] is True
    assert contract["source_projected_amplitude_rhs_materialized"] is True
    assert contract["source_constraint_evolution_identity_executable"] is True
    assert contract["caller_supplies_n_phi_prime"] is False
    assert contract["caller_supplies_mode_pressure_pi_m"] is False
    assert contract["caller_supplies_mode_state_t_m"] is True
    assert contract["caller_supplies_mode_forcing_f_m"] is True
    assert contract["source_zero_data_ivp_solver_materialized"] is False
    assert contract["source_amplitude_directional_jet_materialized"] is False
    assert contract["source_amplitude_ode_materialized"] is False
    assert contract["self_contained_velocity_xyzt_provider"] is False
    assert contract["paper_exact"] is False
    assert contract["complete_ns_residual_assessed"] is False
    assert contract["residual_reduction_claimed"] is False
    assert contract["pde_validated"] is False
