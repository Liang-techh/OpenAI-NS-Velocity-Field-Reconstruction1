from __future__ import annotations

import inspect

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_source_phase_coefficient_jet import (
    SourceBackgroundPhaseJet,
    complete_harmonic_from_source_phase_jet,
    source_phase_vector_jet,
)
from openai_ns_reconstruction.kokuno_source_zero_data_amplitude_ivp import (
    solve_source_zero_data_amplitude_ivp,
)
from openai_ns_reconstruction.kokuno_source_zero_data_amplitude_sensitivity import (
    SourceModeForcingDirectionalJet,
    solve_source_zero_data_amplitude_sensitivity,
    source_zero_data_amplitude_sensitivity_contract,
)


def _background(
    radius: np.ndarray | float,
    z_value: np.ndarray | float,
) -> SourceBackgroundPhaseJet:
    r = np.asarray(radius, dtype=float)
    z = np.asarray(z_value, dtype=float)
    r, z = np.broadcast_arrays(r, z)
    f = 0.12 * r**2 + 0.07 * r * z + 0.04 * z**2 + 0.03 * r + 0.02 * z + 0.1
    g = -0.05 * r**2 + 0.09 * r * z + 0.03 * z**2 - 0.02 * r + 0.04 * z
    return SourceBackgroundPhaseJet(
        f=f,
        g=g,
        f_r=0.24 * r + 0.07 * z + 0.03,
        g_r=-0.10 * r + 0.09 * z - 0.02,
        f_z=0.07 * r + 0.08 * z + 0.02,
        g_z=0.09 * r + 0.06 * z + 0.04,
        f_rr=np.full_like(r, 0.24),
        g_rr=np.full_like(r, -0.10),
        f_rz=np.full_like(r, 0.07),
        g_rz=np.full_like(r, 0.09),
        f_zz=np.full_like(r, 0.08),
        g_zz=np.full_like(r, 0.06),
    )


_FORCING_COEFFICIENT = np.array(
    [0.20 + 0.03j, -0.11 + 0.07j, 0.05 - 0.02j],
    dtype=np.complex128,
)


def _forcing_value(
    radius: np.ndarray | float,
    z_value: np.ndarray | float,
    pulse_v: float,
) -> np.ndarray:
    r = np.asarray(radius, dtype=float)
    z = np.asarray(z_value, dtype=float)
    r, z = np.broadcast_arrays(r, z)
    shape = np.stack(
        (
            1.0 + 0.10 * r + 0.05 * z + 0.20 * pulse_v,
            1.0 - 0.08 * r + 0.03 * z - 0.10 * pulse_v,
            1.0 + 0.02 * r - 0.04 * z + 0.15 * pulse_v,
        ),
        axis=-1,
    )
    return shape * _FORCING_COEFFICIENT


def _forcing_jet(
    radius: np.ndarray | float,
    z_value: np.ndarray | float,
    *,
    epsilon: float,
):
    r = np.asarray(radius, dtype=float)
    z = np.asarray(z_value, dtype=float)
    r, z = np.broadcast_arrays(r, z)
    ones = np.ones(r.shape + (1,), dtype=float)
    dr = ones * np.array([0.10, -0.08, 0.02]) * _FORCING_COEFFICIENT
    dz = (
        epsilon
        * ones
        * np.array([0.05, 0.03, -0.04])
        * _FORCING_COEFFICIENT
    )

    def evaluate(pulse_v: float) -> SourceModeForcingDirectionalJet:
        return SourceModeForcingDirectionalJet(
            value=_forcing_value(r, z, pulse_v),
            dr_value=dr,
            dz_value=dz,
        )

    return evaluate


def _solve_parent(
    radius: float,
    z_value: float,
    *,
    pulse_v: float,
    epsilon: float,
    p: float,
    p_z: float,
    x_0: float,
    k: float,
    m: int,
) -> np.ndarray:
    return solve_source_zero_data_amplitude_ivp(
        pulse_v,
        radius,
        0.27,
        z_value,
        _background(radius, z_value),
        lambda v: _forcing_value(radius, z_value, v),
        epsilon=epsilon,
        p=p,
        p_z=p_z,
        x_0=x_0,
        k=k,
        m=m,
    ).state


def test_sensitivity_matches_independent_parent_ivp_finite_differences() -> None:
    radius = 1.2
    z_value = -0.25
    pulse_v = 0.4
    epsilon = 0.15
    p = 1.1
    p_z = -0.4
    x_0 = 0.22
    k = 2.2
    m = 1

    got = solve_source_zero_data_amplitude_sensitivity(
        pulse_v,
        radius,
        0.27,
        z_value,
        _background(radius, z_value),
        _forcing_jet(radius, z_value, epsilon=epsilon),
        epsilon=epsilon,
        p=p,
        p_z=p_z,
        x_0=x_0,
        k=k,
        m=m,
    )
    center = _solve_parent(
        radius,
        z_value,
        pulse_v=pulse_v,
        epsilon=epsilon,
        p=p,
        p_z=p_z,
        x_0=x_0,
        k=k,
        m=m,
    )
    np.testing.assert_allclose(got.amplitude_jet.t_m, center, rtol=0.0, atol=2.0e-13)

    h = 3.0e-5
    plus_r = _solve_parent(
        radius + h,
        z_value,
        pulse_v=pulse_v,
        epsilon=epsilon,
        p=p,
        p_z=p_z,
        x_0=x_0,
        k=k,
        m=m,
    )
    minus_r = _solve_parent(
        radius - h,
        z_value,
        pulse_v=pulse_v,
        epsilon=epsilon,
        p=p,
        p_z=p_z,
        x_0=x_0,
        k=k,
        m=m,
    )
    dr_reference = (plus_r - minus_r) / (2.0 * h)

    plus_z = _solve_parent(
        radius,
        z_value + h,
        pulse_v=pulse_v,
        epsilon=epsilon,
        p=p,
        p_z=p_z,
        x_0=x_0,
        k=k,
        m=m,
    )
    minus_z = _solve_parent(
        radius,
        z_value - h,
        pulse_v=pulse_v,
        epsilon=epsilon,
        p=p,
        p_z=p_z,
        x_0=x_0,
        k=k,
        m=m,
    )
    dz_reference = epsilon * (plus_z - minus_z) / (2.0 * h)

    np.testing.assert_allclose(
        got.amplitude_jet.dr_t_m, dr_reference, rtol=0.0, atol=5.0e-8
    )
    np.testing.assert_allclose(
        got.amplitude_jet.dz_t_m, dz_reference, rtol=0.0, atol=5.0e-8
    )
    assert float(np.max(np.abs(got.endpoint.dr_constraint))) <= 5.0e-10
    assert float(np.max(np.abs(got.endpoint.dz_constraint))) <= 5.0e-10
    assert float(np.max(np.abs(got.endpoint.base_ode.full_equation_residual))) <= 5.0e-12


def test_generated_directional_jet_feeds_existing_complete_curl_path() -> None:
    radius = 1.18
    z_value = -0.19
    epsilon = 0.14
    p = 1.0
    p_z = -0.35
    x_0 = 0.2
    k = 2.1
    m = -1
    pulse_v = 0.33

    result = solve_source_zero_data_amplitude_sensitivity(
        pulse_v,
        radius,
        -0.31,
        z_value,
        _background(radius, z_value),
        _forcing_jet(radius, z_value, epsilon=epsilon),
        epsilon=epsilon,
        p=p,
        p_z=p_z,
        x_0=x_0,
        k=k,
        m=m,
    )
    phase = source_phase_vector_jet(
        radius,
        -0.31,
        z_value,
        pulse_v,
        _background(radius, z_value),
        epsilon=epsilon,
        p=p,
        p_z=p_z,
        x_0=x_0,
    )
    complete = complete_harmonic_from_source_phase_jet(
        radius,
        phase,
        result.amplitude_jet,
        k=k,
        m=m,
    )
    assert np.all(np.isfinite(complete.complete_amplitude.coefficient))
    assert np.all(np.isfinite(complete.complete_amplitude.remainder))
    assert np.all(np.isfinite(complete.complete_amplitude.amplitude))


def test_sensitivity_solver_is_batch_safe() -> None:
    radius = np.array([0.95, 1.25, 1.55])
    z_value = np.array([-0.21, 0.08, 0.27])
    epsilon = 0.13
    got = solve_source_zero_data_amplitude_sensitivity(
        0.28,
        radius,
        np.array([0.1, -0.2, 0.3]),
        z_value,
        _background(radius, z_value),
        _forcing_jet(radius, z_value, epsilon=epsilon),
        epsilon=epsilon,
        p=1.2,
        p_z=-0.3,
        x_0=0.18,
        k=2.3,
        m=2,
    )
    assert got.amplitude_jet.t_m.shape == (3, 3)
    assert got.amplitude_jet.dr_t_m.shape == (3, 3)
    assert got.amplitude_jet.dz_t_m.shape == (3, 3)
    assert np.all(np.isfinite(got.amplitude_jet.t_m))
    assert np.all(np.isfinite(got.amplitude_jet.dr_t_m))
    assert np.all(np.isfinite(got.amplitude_jet.dz_t_m))
    assert float(np.max(np.abs(got.endpoint.dr_constraint))) <= 2.0e-9
    assert float(np.max(np.abs(got.endpoint.dz_constraint))) <= 2.0e-9


def test_zero_pulse_returns_exact_zero_amplitude_and_directional_jets() -> None:
    epsilon = 0.12
    got = solve_source_zero_data_amplitude_sensitivity(
        0.0,
        1.1,
        0.2,
        -0.1,
        _background(1.1, -0.1),
        _forcing_jet(1.1, -0.1, epsilon=epsilon),
        epsilon=epsilon,
        p=1.0,
        p_z=0.2,
        x_0=0.3,
        k=2.0,
        m=1,
    )
    np.testing.assert_array_equal(got.amplitude_jet.t_m, np.zeros(3, dtype=complex))
    np.testing.assert_array_equal(got.amplitude_jet.dr_t_m, np.zeros(3, dtype=complex))
    np.testing.assert_array_equal(got.amplitude_jet.dz_t_m, np.zeros(3, dtype=complex))


def test_integrator_accuracy_is_not_a_caller_knob() -> None:
    parameters = inspect.signature(
        solve_source_zero_data_amplitude_sensitivity
    ).parameters
    for forbidden in ("panels", "steps", "rtol", "atol", "tol", "max_step"):
        assert forbidden not in parameters


def test_forcing_directional_jet_is_typed_and_fails_closed() -> None:
    common = dict(
        pulse_v=0.2,
        radius=1.0,
        theta=0.1,
        z_normalized=-0.2,
        background=_background(1.0, -0.2),
        epsilon=0.15,
        p=1.0,
        p_z=0.3,
        x_0=0.2,
        k=2.0,
        m=1,
    )
    with pytest.raises(ValueError):
        solve_source_zero_data_amplitude_sensitivity(
            forcing_jet_m=lambda _v: np.ones(3), **common  # type: ignore[arg-type]
        )
    with pytest.raises(ValueError):
        solve_source_zero_data_amplitude_sensitivity(
            forcing_jet_m=lambda _v: SourceModeForcingDirectionalJet(
                value=np.ones(3),
                dr_value=np.array([1.0, np.nan, 0.0]),
                dz_value=np.zeros(3),
            ),
            **common,
        )


def test_contract_closes_amplitude_jet_only_not_upstream_source_providers() -> None:
    contract = source_zero_data_amplitude_sensitivity_contract()
    assert contract["source_commit"] == "143f6773feb424ad9ed3a8d116653200f20346b7"
    assert contract["panel_ladder"] == [32, 64, 128]
    assert contract["source_normalized_dr_operator_has_auxiliary_L_term"] is True
    assert contract["background_phase_coefficients_auxiliary_independent"] is True
    assert contract["forcing_dr_input_is_full_normalized_source_Dr"] is True
    assert contract["caller_can_tune_integrator_resolution"] is False
    assert contract["caller_supplies_mode_state_t_m"] is False
    assert contract["caller_supplies_amplitude_directional_jet"] is False
    assert contract["caller_supplies_mode_forcing_directional_jet"] is True
    assert contract["caller_supplies_corrected_background"] is True
    assert contract["source_amplitude_directional_jet_materialized"] is True
    assert contract["source_complete_curl_directional_input_compatible"] is True
    assert contract["source_forcing_provider_materialized"] is False
    assert contract["full_source_forcing_auxiliary_provider_materialized"] is False
    assert contract["actual_corrected_background_provider_materialized"] is False
    assert contract["source_amplitude_mode_provider_complete"] is False
    assert contract["self_contained_velocity_xyzt_provider"] is False
    assert contract["paper_exact"] is False
    assert contract["complete_ns_residual_assessed"] is False
    assert contract["residual_reduction_claimed"] is False
    assert contract["pde_validated"] is False
