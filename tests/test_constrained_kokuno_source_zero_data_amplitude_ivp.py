from __future__ import annotations

import inspect

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_source_phase_coefficient_jet import (
    SourceBackgroundPhaseJet,
)
from openai_ns_reconstruction.kokuno_source_zero_data_amplitude_ivp import (
    IVP_PANEL_LADDER,
    solve_source_zero_data_amplitude_ivp,
    source_zero_data_amplitude_ivp_contract,
)


def _zero_background(shape: tuple[int, ...] = ()) -> SourceBackgroundPhaseJet:
    zero = np.zeros(shape, dtype=float)
    return SourceBackgroundPhaseJet(
        f=zero,
        g=zero,
        f_r=zero,
        g_r=zero,
        f_z=zero,
        g_z=zero,
        f_rr=zero,
        g_rr=zero,
        f_rz=zero,
        g_rz=zero,
        f_zz=zero,
        g_zz=zero,
    )


def _constant_transverse_forcing(
    radius: np.ndarray,
    *,
    p: float,
    p_z: float,
    x_0: float,
) -> np.ndarray:
    radius = np.asarray(radius, dtype=float)
    n = np.stack(
        (
            np.broadcast_to(x_0, radius.shape),
            p / radius,
            np.broadcast_to(p_z, radius.shape),
        ),
        axis=-1,
    )
    seed = np.broadcast_to(
        np.array([0.31 + 0.07j, -0.22 + 0.11j, 0.43 - 0.09j]),
        n.shape,
    )
    return np.cross(n, seed)


def test_zero_data_ivp_matches_independent_constant_coefficient_solution() -> None:
    radius = np.asarray(1.1)
    epsilon = 0.18
    p = 1.3
    p_z = -0.4
    x_0 = 0.2
    k = 2.5
    m = 1
    pulse_v = 0.8
    forcing = _constant_transverse_forcing(radius, p=p, p_z=p_z, x_0=x_0)

    got = solve_source_zero_data_amplitude_ivp(
        pulse_v,
        radius,
        0.25,
        -0.2,
        _zero_background(),
        lambda _v: forcing,
        epsilon=epsilon,
        p=p,
        p_z=p_z,
        x_0=x_0,
        k=k,
        m=m,
    )

    n = np.array([x_0, p / float(radius), p_z])
    gamma = (m * m) * epsilon * (k**2) * float(np.dot(n, n))
    exact = -(1.0 - np.exp(-gamma * pulse_v)) / gamma * forcing
    errors = [float(np.max(np.abs(state - exact))) for state in got.states]

    assert got.panels == IVP_PANEL_LADDER == (32, 64, 128)
    assert errors[2] < errors[1] < errors[0]
    assert got.medium_fine_abs_max < got.coarse_medium_abs_max
    assert errors[2] <= 2.0e-10
    np.testing.assert_allclose(got.state, exact, rtol=0.0, atol=2.0e-10)
    assert abs(complex(np.dot(n, got.state))) <= 2.0e-13
    assert abs(complex(got.endpoint.constraint)) <= 2.0e-13
    assert float(np.max(np.abs(got.endpoint.full_equation_residual))) <= 3.0e-13


def test_zero_data_ivp_is_batch_safe_and_preserves_transversality() -> None:
    radius = np.array([0.9, 1.4, 1.8])
    p = -1.7
    p_z = 0.35
    x_0 = -0.28
    forcing = _constant_transverse_forcing(radius, p=p, p_z=p_z, x_0=x_0)

    got = solve_source_zero_data_amplitude_ivp(
        0.55,
        radius,
        np.array([0.1, -0.3, 0.5]),
        np.array([-0.2, 0.4, 0.1]),
        _zero_background((3,)),
        lambda _v: forcing,
        epsilon=0.14,
        p=p,
        p_z=p_z,
        x_0=x_0,
        k=3.0,
        m=-2,
    )

    n = np.stack(
        (
            np.full_like(radius, x_0),
            p / radius,
            np.full_like(radius, p_z),
        ),
        axis=-1,
    )
    assert got.state.shape == (3, 3)
    assert np.all(np.isfinite(got.state.real))
    assert np.all(np.isfinite(got.state.imag))
    assert float(np.max(np.abs(np.sum(n * got.state, axis=-1)))) <= 3.0e-12
    assert float(np.max(np.abs(got.endpoint.constraint))) <= 3.0e-12


def test_zero_target_returns_exact_zero_state() -> None:
    radius = np.asarray(1.2)
    forcing = np.array([0.03 + 0.01j, -0.02j, 0.04 - 0.01j])
    got = solve_source_zero_data_amplitude_ivp(
        0.0,
        radius,
        0.1,
        0.2,
        _zero_background(),
        lambda _v: forcing,
        epsilon=0.16,
        p=1.0,
        p_z=0.2,
        x_0=0.3,
        k=2.0,
        m=1,
    )
    np.testing.assert_array_equal(got.state, np.zeros(3, dtype=np.complex128))
    assert got.coarse_medium_abs_max == 0.0
    assert got.medium_fine_abs_max == 0.0


def test_signed_negative_pulse_path_is_supported_by_same_fixed_integrator() -> None:
    radius = np.asarray(1.05)
    epsilon = 0.12
    p = 1.1
    p_z = -0.25
    x_0 = 0.15
    k = 2.0
    m = 1
    pulse_v = -0.35
    forcing = _constant_transverse_forcing(radius, p=p, p_z=p_z, x_0=x_0)
    got = solve_source_zero_data_amplitude_ivp(
        pulse_v,
        radius,
        -0.2,
        0.3,
        _zero_background(),
        lambda _v: forcing,
        epsilon=epsilon,
        p=p,
        p_z=p_z,
        x_0=x_0,
        k=k,
        m=m,
    )
    n = np.array([x_0, p / float(radius), p_z])
    gamma = epsilon * (k**2) * float(np.dot(n, n))
    exact = -(1.0 - np.exp(-gamma * pulse_v)) / gamma * forcing
    np.testing.assert_allclose(got.state, exact, rtol=0.0, atol=2.0e-10)


def test_integrator_resolution_and_tolerances_are_not_caller_knobs() -> None:
    parameters = inspect.signature(solve_source_zero_data_amplitude_ivp).parameters
    for forbidden in ("panels", "steps", "rtol", "atol", "tol", "max_step"):
        assert forbidden not in parameters


def test_zero_data_ivp_fails_closed_on_nonfinite_or_invalid_inputs() -> None:
    bg = _zero_background()
    common = dict(
        radius=1.0,
        theta=0.1,
        z_normalized=0.2,
        background=bg,
        epsilon=0.15,
        p=1.0,
        p_z=0.3,
        x_0=0.2,
        k=2.0,
        m=1,
    )
    with pytest.raises(ValueError):
        solve_source_zero_data_amplitude_ivp(
            np.nan, forcing_m=lambda _v: np.ones(3), **common
        )
    with pytest.raises(ValueError):
        solve_source_zero_data_amplitude_ivp(
            0.2, forcing_m=None, **common  # type: ignore[arg-type]
        )
    with pytest.raises(ValueError):
        solve_source_zero_data_amplitude_ivp(
            0.2, forcing_m=lambda _v: np.array([1.0, np.nan, 0.0]), **common
        )


@pytest.mark.parametrize(("epsilon", "k", "m"), [(0.0, 2.0, 1), (0.2, 0.0, 1), (0.2, 2.0, 0)])
def test_zero_data_ivp_rejects_invalid_source_scalars(
    epsilon: float, k: float, m: int
) -> None:
    with pytest.raises(ValueError):
        solve_source_zero_data_amplitude_ivp(
            0.2,
            1.0,
            0.1,
            0.2,
            _zero_background(),
            lambda _v: np.ones(3),
            epsilon=epsilon,
            p=1.0,
            p_z=0.3,
            x_0=0.2,
            k=k,
            m=m,
        )


def test_contract_labels_numerical_ivp_without_promoting_source_exactness() -> None:
    contract = source_zero_data_amplitude_ivp_contract()
    assert contract["source_commit"] == "143f6773feb424ad9ed3a8d116653200f20346b7"
    assert contract["panel_ladder"] == [32, 64, 128]
    assert contract["caller_can_tune_integrator_resolution"] is False
    assert contract["source_phase_jet_regenerated_analytically_along_v"] is True
    assert contract["caller_supplies_mode_state_t_m"] is False
    assert contract["caller_supplies_mode_forcing_f_m"] is True
    assert contract["caller_supplies_corrected_background"] is True
    assert contract["source_zero_data_ivp_numerical_solver_materialized"] is True
    assert contract["source_exact_duhamel_frame_B_materialized"] is False
    assert contract["source_exact_duhamel_propagator_Vm_materialized"] is False
    assert contract["source_amplitude_directional_jet_materialized"] is False
    assert contract["source_amplitude_mode_provider_complete"] is False
    assert contract["self_contained_velocity_xyzt_provider"] is False
    assert contract["paper_exact"] is False
    assert contract["complete_ns_residual_assessed"] is False
    assert contract["residual_reduction_claimed"] is False
    assert contract["pde_validated"] is False
