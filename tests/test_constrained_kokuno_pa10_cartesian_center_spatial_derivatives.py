from __future__ import annotations

import copy
import math

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_pa10_cartesian_center_spatial_derivatives import (
    KokunoPA10CartesianCenterSpatialDerivatives,
)


def _physical_probe(spatial: KokunoPA10CartesianCenterSpatialDerivatives):
    _, xmax = spatial.source_X_interval
    X = np.asarray([0.06, 0.14, 0.23, 0.31], dtype=float) * xmax
    eta = np.asarray([0.08, -0.17, 0.27, -0.36], dtype=float)
    t = np.asarray([0.34, 0.46, 0.54, 0.66], dtype=float)
    theta = np.asarray([0.25, 0.85, 1.45, 2.05], dtype=float)
    return X, eta, spatial.field.cartesian_from_similarity(X, eta, t, theta)


def _fd4_axis(spatial, physical, axis, step):
    coords = [
        np.asarray(physical["x"], dtype=float),
        np.asarray(physical["y"], dtype=float),
        np.asarray(physical["z"], dtype=float),
    ]

    def shifted(multiplier):
        moved = [item.copy() for item in coords]
        moved[axis] = moved[axis] + multiplier * step
        return spatial.field.velocity(moved[0], moved[1], moved[2], physical["t"])

    return (
        shifted(-2.0)
        - 8.0 * shifted(-1.0)
        + 8.0 * shifted(1.0)
        - shifted(2.0)
    ) / (12.0 * step)


def _fd4_jacobian(spatial, physical, step):
    return np.stack(
        [_fd4_axis(spatial, physical, axis, step) for axis in range(3)],
        axis=-1,
    )


def _relative_max(actual, reference):
    actual = np.asarray(actual, dtype=float)
    reference = np.asarray(reference, dtype=float)
    scale = np.maximum(np.maximum(np.abs(actual), np.abs(reference)), 1.0)
    return float(np.max(np.abs(actual - reference) / scale))


def test_analytic_coordinate_spatial_derivatives_against_independent_fd():
    spatial = KokunoPA10CartesianCenterSpatialDerivatives()
    _, _, physical = _physical_probe(spatial)
    analytic = spatial.coordinate_spatial_derivatives(
        physical["x"], physical["y"], physical["z"], physical["t"]
    )
    step = 1.0e-5
    plus = spatial.field.similarity_coordinates(
        physical["x"], physical["y"], physical["z"] + step, physical["t"]
    )
    minus = spatial.field.similarity_coordinates(
        physical["x"], physical["y"], physical["z"] - step, physical["t"]
    )
    q_z_fd = (plus["q"] - minus["q"]) / (2.0 * step)
    X_z_fd = (plus["X"] - minus["X"]) / (2.0 * step)
    eta_z_fd = (plus["eta"] - minus["eta"]) / (2.0 * step)
    assert _relative_max(analytic["q_z"], q_z_fd) < 2e-8
    assert _relative_max(analytic["X_z"], X_z_fd) < 2e-8
    assert _relative_max(analytic["eta_z"], eta_z_fd) < 2e-8

    np.testing.assert_array_equal(analytic["q_x"], np.zeros_like(analytic["q"]))
    np.testing.assert_array_equal(analytic["q_y"], np.zeros_like(analytic["q"]))
    np.testing.assert_array_equal(analytic["eta_x"], np.zeros_like(analytic["eta"]))
    np.testing.assert_array_equal(analytic["eta_y"], np.zeros_like(analytic["eta"]))
    np.testing.assert_allclose(
        analytic["X_x"], physical["x"] / analytic["q"], rtol=0.0, atol=0.0
    )
    np.testing.assert_allclose(
        analytic["X_y"], physical["y"] / analytic["q"], rtol=0.0, atol=0.0
    )


def test_velocity_jacobian_matches_independent_fixed_cartesian_fd4():
    spatial = KokunoPA10CartesianCenterSpatialDerivatives()
    _, _, physical = _physical_probe(spatial)
    analytic = spatial.velocity_jacobian(
        physical["x"], physical["y"], physical["z"], physical["t"]
    )
    coarse = _fd4_jacobian(spatial, physical, 2.0e-4)
    fine = _fd4_jacobian(spatial, physical, 1.0e-4)
    assert analytic.shape == (4, 3, 3)
    assert np.all(np.isfinite(analytic))
    assert np.any(np.abs(analytic) > 0.0)
    assert _relative_max(analytic, fine) < 5e-5
    assert _relative_max(fine, coarse) < 5e-5

    values = spatial.values(
        physical["x"], physical["y"], physical["z"], physical["t"]
    )
    np.testing.assert_array_equal(analytic, values["velocity_jacobian"])


def test_divergence_and_vorticity_are_analytic_and_consistent():
    spatial = KokunoPA10CartesianCenterSpatialDerivatives()
    _, _, physical = _physical_probe(spatial)
    jacobian = spatial.velocity_jacobian(
        physical["x"], physical["y"], physical["z"], physical["t"]
    )
    divergence = spatial.divergence(
        physical["x"], physical["y"], physical["z"], physical["t"]
    )
    vorticity = spatial.vorticity(
        physical["x"], physical["y"], physical["z"], physical["t"]
    )
    np.testing.assert_array_equal(
        divergence, np.trace(jacobian, axis1=-2, axis2=-1)
    )
    expected_vorticity = np.stack(
        (
            jacobian[..., 2, 1] - jacobian[..., 1, 2],
            jacobian[..., 0, 2] - jacobian[..., 2, 0],
            jacobian[..., 1, 0] - jacobian[..., 0, 1],
        ),
        axis=-1,
    )
    np.testing.assert_array_equal(vorticity, expected_vorticity)

    fd = _fd4_jacobian(spatial, physical, 1.0e-4)
    fd_curl = np.stack(
        (
            fd[..., 2, 1] - fd[..., 1, 2],
            fd[..., 0, 2] - fd[..., 2, 0],
            fd[..., 1, 0] - fd[..., 0, 1],
        ),
        axis=-1,
    )
    assert _relative_max(vorticity, fd_curl) < 5e-5
    assert np.any(np.abs(vorticity) > 0.0)

    scale = max(1.0, float(np.max(np.abs(jacobian))))
    assert float(np.max(np.abs(divergence))) / scale < 1e-9


def test_spatial_jacobian_axis_regular_and_rotation_covariant():
    spatial = KokunoPA10CartesianCenterSpatialDerivatives()
    axis = spatial.field.cartesian_from_similarity(
        np.zeros(3),
        np.asarray([-0.25, 0.0, 0.25]),
        np.asarray([0.36, 0.50, 0.64]),
    )
    axis_jacobian = spatial.velocity_jacobian(
        axis["x"], axis["y"], axis["z"], axis["t"]
    )
    assert np.all(np.isfinite(axis_jacobian))
    assert np.all(
        np.isfinite(spatial.vorticity(axis["x"], axis["y"], axis["z"], axis["t"]))
    )

    _, xmax = spatial.source_X_interval
    X = 0.20 * xmax
    eta = 0.19
    t = 0.51
    angle = 0.63
    p0 = spatial.field.cartesian_from_similarity(X, eta, t, 0.0)
    p1 = spatial.field.cartesian_from_similarity(X, eta, t, angle)
    j0 = spatial.velocity_jacobian(p0["x"], p0["y"], p0["z"], p0["t"])
    j1 = spatial.velocity_jacobian(p1["x"], p1["y"], p1["z"], p1["t"])
    w0 = spatial.vorticity(p0["x"], p0["y"], p0["z"], p0["t"])
    w1 = spatial.vorticity(p1["x"], p1["y"], p1["z"], p1["t"])

    c, s = math.cos(angle), math.sin(angle)
    rotation = np.asarray([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])
    np.testing.assert_allclose(j1, rotation @ j0 @ rotation.T, rtol=5e-10, atol=1e-7)
    np.testing.assert_allclose(w1, rotation @ w0, rtol=5e-10, atol=1e-7)


def test_configuration_roundtrip_and_fail_closed_truth_boundary():
    spatial = KokunoPA10CartesianCenterSpatialDerivatives()
    payload = spatial.configuration()
    replay = KokunoPA10CartesianCenterSpatialDerivatives.from_configuration(
        copy.deepcopy(payload)
    )
    assert replay.configuration() == payload
    assert replay.field_sha256 == spatial.field_sha256
    assert replay.temporal_derivative_sha256 == spatial.temporal_derivative_sha256
    assert replay.spatial_derivative_sha256 == spatial.spatial_derivative_sha256

    _, _, physical = _physical_probe(spatial)
    np.testing.assert_array_equal(
        replay.velocity_jacobian(
            physical["x"], physical["y"], physical["z"], physical["t"]
        ),
        spatial.velocity_jacobian(
            physical["x"], physical["y"], physical["z"], physical["t"]
        ),
    )

    with pytest.raises(ValueError, match="registered interval"):
        spatial.velocity_jacobian(0.0, 0.0, 0.0, 0.9)

    truth = spatial.truth_boundary
    assert truth["source_coordinate_spatial_derivatives_analytic"] is True
    assert truth["inner_cartesian_center_velocity_spatial_jacobian_executable"] is True
    assert truth["inner_cartesian_center_divergence_executable"] is True
    assert truth["inner_cartesian_center_vorticity_executable"] is True
    assert truth["spatial_derivatives_use_finite_difference_in_production"] is False
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
    spatial = KokunoPA10CartesianCenterSpatialDerivatives()
    report = spatial.report()
    checks = report["machine_checks"]
    assert checks["velocity_jacobian_vs_fd4_fine_relative_max"] < 5e-5
    assert checks["fd4_fine_vs_coarse_relative_max"] < 5e-5
    assert checks["vorticity_vs_fd4_curl_relative_max"] < 5e-5
    assert checks["analytic_divergence_relative_to_jacobian_scale"] < 1e-9
    assert checks["q_z_vs_centered_fd_relative_max"] < 2e-8
    assert checks["X_z_vs_centered_fd_relative_max"] < 2e-8
    assert checks["eta_z_vs_centered_fd_relative_max"] < 2e-8
    assert checks["rotation_covariance_relative_max"] < 5e-10
    assert checks["axis_jacobian_all_finite"] is True
    assert checks["vorticity_nontrivial_on_probe"] is True
    assert checks["all_spatial_derivative_probe_values_finite"] is True

    gates = report["scientific_gates"]
    assert gates["momentum_max_l2"] == 1e-3
    assert gates["divergence_max_l2"] == 1e-5
    assert gates["free_residual_defined_forcing_forbidden"] is True
    truth = report["truth_boundary"]
    assert truth["complete_candidate_api_ready"] is False
    assert truth["heldout_ns_residual_assessed"] is False
    assert truth["pde_validated"] is False
