from __future__ import annotations

from dataclasses import replace

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_source_background_phase_jet import (
    SourceActualBackgroundVGJet,
    source_background_phase_jet_contract,
    source_background_phase_jet_from_vg,
    source_phase_vector_jet_from_vg,
)
from openai_ns_reconstruction.kokuno_source_phase_coefficient_jet import (
    SourceBackgroundPhaseJet,
    source_phase_vector_jet,
)


def _manufactured_vg(radius: np.ndarray, z: np.ndarray) -> SourceActualBackgroundVGJet:
    # V = R (R^2 + Z^2), hence the source identity F=V/R gives F=R^2+Z^2.
    v = radius * (radius**2 + z**2)
    return SourceActualBackgroundVGJet(
        v=v,
        g=radius * z,
        v_r=3.0 * radius**2 + z**2,
        g_r=z,
        v_z=2.0 * radius * z,
        g_z=radius,
        v_rr=6.0 * radius,
        g_rr=np.zeros_like(radius),
        v_rz=2.0 * z,
        g_rz=np.ones_like(radius),
        v_zz=2.0 * radius,
        g_zz=np.zeros_like(radius),
    )


def _expected_fg(radius: np.ndarray, z: np.ndarray) -> SourceBackgroundPhaseJet:
    return SourceBackgroundPhaseJet(
        f=radius**2 + z**2,
        g=radius * z,
        f_r=2.0 * radius,
        g_r=z,
        f_z=2.0 * z,
        g_z=radius,
        f_rr=2.0 * np.ones_like(radius),
        g_rr=np.zeros_like(radius),
        f_rz=np.zeros_like(radius),
        g_rz=np.ones_like(radius),
        f_zz=2.0 * np.ones_like(radius),
        g_zz=np.zeros_like(radius),
    )


def _assert_background_close(
    got: SourceBackgroundPhaseJet, expected: SourceBackgroundPhaseJet
) -> None:
    for name in (
        "f",
        "g",
        "f_r",
        "g_r",
        "f_z",
        "g_z",
        "f_rr",
        "g_rr",
        "f_rz",
        "g_rz",
        "f_zz",
        "g_zz",
    ):
        np.testing.assert_allclose(
            getattr(got, name), getattr(expected, name), rtol=0.0, atol=2e-13
        )


def test_f_equals_v_over_r_and_all_required_derivatives_are_analytic() -> None:
    radius = np.array([0.71, 1.13, 1.82])
    z = np.array([-0.37, 0.22, 0.49])

    got = source_background_phase_jet_from_vg(radius, _manufactured_vg(radius, z))
    expected = _expected_fg(radius, z)
    _assert_background_close(got, expected)


def test_background_adapter_replays_parent_phase_vector_jet_without_free_h_derivatives() -> None:
    radius = np.array([0.83, 1.27])
    z = np.array([-0.31, 0.44])
    theta = np.array([0.2, -0.6])
    pulse_v = np.array([0.35, 0.58])
    epsilon = 0.17
    p = 2.0
    p_z = -0.75
    x_0 = 0.41

    actual = _manufactured_vg(radius, z)
    adapted = source_phase_vector_jet_from_vg(
        radius,
        theta,
        z,
        pulse_v,
        actual,
        epsilon=epsilon,
        p=p,
        p_z=p_z,
        x_0=x_0,
    )
    direct = source_phase_vector_jet(
        radius,
        theta,
        z,
        pulse_v,
        _expected_fg(radius, z),
        epsilon=epsilon,
        p=p,
        p_z=p_z,
        x_0=x_0,
    )

    _assert_background_close(adapted.background_phase_jet, _expected_fg(radius, z))
    for name in ("phase", "h_phi", "n_phi", "dr_n_phi", "dz_n_phi"):
        np.testing.assert_allclose(
            getattr(adapted.phase_vector_jet, name),
            getattr(direct, name),
            rtol=0.0,
            atol=2e-13,
        )


def test_f_derivatives_agree_with_independent_centered_differences() -> None:
    r = 1.19
    z = -0.28
    h = 2.0e-5

    def f_value(rr: float, zz: float) -> float:
        vv = rr * (rr**2 + zz**2)
        return vv / rr

    actual = _manufactured_vg(np.asarray(r), np.asarray(z))
    got = source_background_phase_jet_from_vg(r, actual)

    fd_r = (f_value(r + h, z) - f_value(r - h, z)) / (2.0 * h)
    fd_z = (f_value(r, z + h) - f_value(r, z - h)) / (2.0 * h)
    fd_rr = (
        f_value(r + h, z) - 2.0 * f_value(r, z) + f_value(r - h, z)
    ) / h**2
    fd_rz = (
        f_value(r + h, z + h)
        - f_value(r + h, z - h)
        - f_value(r - h, z + h)
        + f_value(r - h, z - h)
    ) / (4.0 * h**2)
    fd_zz = (
        f_value(r, z + h) - 2.0 * f_value(r, z) + f_value(r, z - h)
    ) / h**2

    assert abs(float(got.f_r) - fd_r) < 2e-9
    assert abs(float(got.f_z) - fd_z) < 2e-9
    assert abs(float(got.f_rr) - fd_rr) < 2e-6
    assert abs(float(got.f_rz) - fd_rz) < 2e-6
    assert abs(float(got.f_zz) - fd_zz) < 2e-6


@pytest.mark.parametrize("radius", [0.0, -0.4, np.nan, np.inf])
def test_background_adapter_fails_closed_off_annular_positive_radius_chart(
    radius: float,
) -> None:
    actual = _manufactured_vg(np.asarray(1.0), np.asarray(0.2))
    with pytest.raises(ValueError):
        source_background_phase_jet_from_vg(radius, actual)


def test_background_adapter_rejects_nonfinite_actual_background_input() -> None:
    actual = _manufactured_vg(np.asarray(1.0), np.asarray(0.2))
    with pytest.raises(ValueError):
        source_background_phase_jet_from_vg(
            1.0,
            replace(actual, v_rr=np.nan),
        )


def test_contract_keeps_source_realization_boundary_hard_false() -> None:
    contract = source_background_phase_jet_contract()

    assert contract["source_f_equals_v_over_r_executed"] is True
    assert contract["source_f_first_second_derivatives_derived_analytically"] is True
    assert contract["caller_supplies_background_phase_derivative_jet"] is False
    assert contract["caller_supplies_actual_background_vg_derivative_jet"] is True
    assert contract["actual_corrected_background_provider_materialized"] is False
    assert contract["agent1_leading_substituted_for_corrected_background"] is False
    assert contract["source_amplitude_ode_materialized"] is False
    assert contract["self_contained_velocity_xyzt_provider"] is False
    assert contract["paper_exact"] is False
    assert contract["complete_ns_residual_assessed"] is False
    assert contract["residual_reduction_claimed"] is False
    assert contract["pde_validated"] is False
