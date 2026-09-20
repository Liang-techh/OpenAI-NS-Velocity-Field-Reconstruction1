from __future__ import annotations

import copy
import math

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_pa10_cartesian_center_velocity_dt import (
    KokunoPA10CartesianCenterVelocityTimeDerivative,
)


def _physical_probe(derivative: KokunoPA10CartesianCenterVelocityTimeDerivative):
    _, xmax = derivative.source_X_interval
    X = np.asarray([0.04, 0.13, 0.29, 0.47], dtype=float) * xmax
    eta = np.asarray([0.07, -0.19, 0.31, -0.43], dtype=float)
    t = np.asarray([0.31, 0.43, 0.57, 0.69], dtype=float)
    theta = np.asarray([0.2, 0.8, 1.4, 2.2], dtype=float)
    return X, eta, derivative.field.cartesian_from_similarity(X, eta, t, theta)


def _fd4_time(derivative, physical, step):
    f = derivative.field.velocity
    x, y, z, t = (
        physical["x"],
        physical["y"],
        physical["z"],
        physical["t"],
    )
    return (
        f(x, y, z, t - 2.0 * step)
        - 8.0 * f(x, y, z, t - step)
        + 8.0 * f(x, y, z, t + step)
        - f(x, y, z, t + 2.0 * step)
    ) / (12.0 * step)


def _relative_max(actual, reference):
    actual = np.asarray(actual, dtype=float)
    reference = np.asarray(reference, dtype=float)
    scale = np.maximum(np.maximum(np.abs(actual), np.abs(reference)), 1.0)
    return float(np.max(np.abs(actual - reference) / scale))


def test_analytic_coordinate_time_derivatives_against_independent_fd():
    derivative = KokunoPA10CartesianCenterVelocityTimeDerivative()
    _, _, physical = _physical_probe(derivative)
    analytic = derivative.coordinate_time_derivatives(
        physical["x"], physical["y"], physical["z"], physical["t"]
    )
    step = 1.0e-5
    plus = derivative.field.similarity_coordinates(
        physical["x"], physical["y"], physical["z"], physical["t"] + step
    )
    minus = derivative.field.similarity_coordinates(
        physical["x"], physical["y"], physical["z"], physical["t"] - step
    )
    q_t_fd = (plus["q"] - minus["q"]) / (2.0 * step)
    X_t_fd = (plus["X"] - minus["X"]) / (2.0 * step)
    eta_t_fd = (plus["eta"] - minus["eta"]) / (2.0 * step)
    assert _relative_max(analytic["q_t"], q_t_fd) < 2e-8
    assert _relative_max(analytic["X_t"], X_t_fd) < 2e-8
    assert _relative_max(analytic["eta_t"], eta_t_fd) < 2e-8

    source_identity = 1.0 - 2.0 * derivative.h * analytic["eta"] ** 2
    np.testing.assert_array_equal(analytic["L"], source_identity)
    np.testing.assert_allclose(
        analytic["q_t"], -1.0 / source_identity, rtol=0.0, atol=0.0
    )


def test_analytic_v0_eta_against_independent_profile_difference():
    derivative = KokunoPA10CartesianCenterVelocityTimeDerivative()
    X, eta, _ = _physical_probe(derivative)
    analytic = derivative.v0_eta(X, eta)
    step = 2.0e-5
    plus = derivative.physical_profiles.values(X, eta + step)["v_0"]
    minus = derivative.physical_profiles.values(X, eta - step)["v_0"]
    centered = (plus - minus) / (2.0 * step)
    assert _relative_max(analytic, centered) < 5e-7
    assert np.all(np.isfinite(analytic))


def test_velocity_dt_matches_independent_fixed_cartesian_fd4():
    derivative = KokunoPA10CartesianCenterVelocityTimeDerivative()
    _, _, physical = _physical_probe(derivative)
    analytic = derivative.velocity_dt(
        physical["x"], physical["y"], physical["z"], physical["t"]
    )
    coarse = _fd4_time(derivative, physical, 2.0e-4)
    fine = _fd4_time(derivative, physical, 1.0e-4)
    assert analytic.shape == (4, 3)
    assert np.all(np.isfinite(analytic))
    assert np.any(np.abs(analytic) > 0.0)
    assert _relative_max(analytic, fine) < 2e-5
    assert _relative_max(fine, coarse) < 2e-5

    values = derivative.values(
        physical["x"], physical["y"], physical["z"], physical["t"]
    )
    np.testing.assert_array_equal(
        analytic, np.stack((values["u_t"], values["v_t"], values["w_t"]), axis=-1)
    )
    np.testing.assert_array_equal(
        derivative.velocity(
            physical["x"], physical["y"], physical["z"], physical["t"]
        ),
        derivative.field.velocity(
            physical["x"], physical["y"], physical["z"], physical["t"]
        ),
    )


def test_velocity_dt_axis_regular_and_rotation_equivariant():
    derivative = KokunoPA10CartesianCenterVelocityTimeDerivative()
    axis = derivative.field.cartesian_from_similarity(
        np.zeros(3),
        np.asarray([-0.25, 0.0, 0.25]),
        np.asarray([0.32, 0.5, 0.68]),
    )
    axis_dt = derivative.velocity_dt(axis["x"], axis["y"], axis["z"], axis["t"])
    assert np.all(axis_dt[..., :2] == 0.0)
    assert np.all(np.isfinite(axis_dt[..., 2]))

    _, xmax = derivative.source_X_interval
    X = 0.25 * xmax
    eta = 0.21
    t = 0.51
    angle = 0.71
    p0 = derivative.field.cartesian_from_similarity(X, eta, t, 0.0)
    p1 = derivative.field.cartesian_from_similarity(X, eta, t, angle)
    d0 = derivative.velocity_dt(p0["x"], p0["y"], p0["z"], p0["t"])
    d1 = derivative.velocity_dt(p1["x"], p1["y"], p1["z"], p1["t"])
    c, s = math.cos(angle), math.sin(angle)
    expected = np.asarray([c * d0[0] - s * d0[1], s * d0[0] + c * d0[1], d0[2]])
    np.testing.assert_allclose(d1, expected, rtol=4e-13, atol=1e-8)


def test_configuration_roundtrip_and_fail_closed_truth_boundary():
    derivative = KokunoPA10CartesianCenterVelocityTimeDerivative()
    payload = derivative.configuration()
    replay = KokunoPA10CartesianCenterVelocityTimeDerivative.from_configuration(
        copy.deepcopy(payload)
    )
    assert replay.configuration() == payload
    assert replay.field_sha256 == derivative.field_sha256
    assert replay.derivative_sha256 == derivative.derivative_sha256

    _, _, physical = _physical_probe(derivative)
    np.testing.assert_array_equal(
        replay.velocity_dt(
            physical["x"], physical["y"], physical["z"], physical["t"]
        ),
        derivative.velocity_dt(
            physical["x"], physical["y"], physical["z"], physical["t"]
        ),
    )

    with pytest.raises(ValueError, match="registered interval"):
        derivative.velocity_dt(0.0, 0.0, 0.0, 0.9)

    truth = derivative.truth_boundary
    assert truth["source_coordinate_time_derivatives_analytic"] is True
    assert truth["center_v0_eta_analytic"] is True
    assert truth["inner_cartesian_center_velocity_dt_executable"] is True
    assert truth["velocity_dt_uses_finite_difference_in_production"] is False
    assert truth["source_center_is_final_corrected_fixed_point"] is False
    assert truth["global_cartesian_spacetime_leading_velocity_materialized"] is False
    assert truth["outer_join_localization_materialized"] is False
    assert truth["matched_global_pressure_materialized"] is False
    assert truth["restricted_forcing_materialized"] is False
    assert truth["complete_candidate_api_ready"] is False
    assert truth["heldout_ns_residual_assessed"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False
    assert truth["openai_field_identified"] is False


def test_report_keeps_pde_and_global_gates_closed():
    derivative = KokunoPA10CartesianCenterVelocityTimeDerivative()
    report = derivative.report()
    checks = report["machine_checks"]
    assert checks["velocity_dt_vs_fd4_fine_relative_max"] < 2e-5
    assert checks["fd4_fine_vs_coarse_relative_max"] < 2e-5
    assert checks["v0_eta_vs_centered_fd_relative_max"] < 5e-7
    assert checks["q_t_vs_centered_fd_relative_max"] < 2e-8
    assert checks["X_t_vs_centered_fd_relative_max"] < 2e-8
    assert checks["eta_t_vs_centered_fd_relative_max"] < 2e-8
    assert checks["axis_transverse_velocity_dt_exact_zero"] is True
    assert checks["velocity_dt_nontrivial_on_probe"] is True
    assert checks["all_velocity_dt_probe_values_finite"] is True

    gates = report["scientific_gates"]
    assert gates["momentum_max_l2"] == 1e-3
    assert gates["divergence_max_l2"] == 1e-5
    assert gates["free_residual_defined_forcing_forbidden"] is True
    truth = report["truth_boundary"]
    assert truth["complete_candidate_api_ready"] is False
    assert truth["heldout_ns_residual_assessed"] is False
    assert truth["pde_validated"] is False
