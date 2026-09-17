import numpy as np
import pytest

from openai_ns_reconstruction.eq45_supported_delivery import Eq45SupportedDeliveryField
from openai_ns_reconstruction.kokuno_complete_curl import (
    KokunoCompleteCurlCorrection,
    compose_velocity,
)


def _correction():
    return KokunoCompleteCurlCorrection(
        amplitude=0.35,
        wave_vector=(5.0, -3.0, 4.0),
        polarization=(1.0, 2.0, 1.0),
        center=(0.1, -0.2, 0.15),
        half_widths=(0.9, 0.8, 0.7),
        omega=2.3,
        phase=0.4,
    )


def _fd_curl_of_potential(field, point, step):
    x, y, z, time = point

    def potential(xx, yy, zz):
        return field.vector_potential(xx, yy, zz, time)

    dx = (potential(x + step, y, z) - potential(x - step, y, z)) / (2.0 * step)
    dy = (potential(x, y + step, z) - potential(x, y - step, z)) / (2.0 * step)
    dz = (potential(x, y, z + step) - potential(x, y, z - step)) / (2.0 * step)
    return np.asarray(
        [
            dy[2] - dz[1],
            dz[0] - dx[2],
            dx[1] - dy[0],
        ]
    )


def _fd_divergence(field, point, step):
    x, y, z, time = point

    def velocity(xx, yy, zz):
        return field.velocity(xx, yy, zz, time)

    dx = (velocity(x + step, y, z) - velocity(x - step, y, z)) / (2.0 * step)
    dy = (velocity(x, y + step, z) - velocity(x, y - step, z)) / (2.0 * step)
    dz = (velocity(x, y, z + step) - velocity(x, y, z - step)) / (2.0 * step)
    return float(dx[0] + dy[1] + dz[2])


def _fd_spatial_jacobian(field, point, step):
    x, y, z, time = point
    columns = []
    for axis in range(3):
        plus = [x, y, z]
        minus = [x, y, z]
        plus[axis] += step
        minus[axis] -= step
        columns.append(
            (
                field.velocity(*plus, time) - field.velocity(*minus, time)
            )
            / (2.0 * step)
        )
    return np.stack(columns, axis=-1)


def _fd_time_derivative(field, point, step):
    x, y, z, time = point
    return (
        field.velocity(x, y, z, time + step)
        - field.velocity(x, y, z, time - step)
    ) / (2.0 * step)


def _fd_laplacian(field, point, step):
    x, y, z, time = point
    center = field.velocity(x, y, z, time)
    total = np.zeros(3)
    for axis in range(3):
        plus = [x, y, z]
        minus = [x, y, z]
        plus[axis] += step
        minus[axis] -= step
        total += field.velocity(*plus, time) - 2.0 * center + field.velocity(*minus, time)
    return total / (step * step)


def test_complete_curl_matches_independent_finite_difference_curl_at_three_resolutions():
    field = _correction()
    point = (0.23, -0.07, 0.31, 0.37)
    exact = field.velocity(*point)

    errors = []
    for step in (0.02, 0.01, 0.005):
        numerical = _fd_curl_of_potential(field, point, step)
        errors.append(float(np.linalg.norm(numerical - exact)))

    assert errors[1] < 0.27 * errors[0]
    assert errors[2] < 0.27 * errors[1]
    assert errors[2] < 7.0e-5


def test_complete_curl_divergence_converges_to_zero_at_three_resolutions():
    field = _correction()
    point = (0.23, -0.07, 0.31, 0.37)
    errors = [abs(_fd_divergence(field, point, step)) for step in (0.02, 0.01, 0.005)]

    assert errors[1] < 0.27 * errors[0]
    assert errors[2] < 0.27 * errors[1]
    assert errors[2] < 5.0e-6


def test_analytic_spatial_jacobian_matches_independent_three_level_fd():
    field = _correction()
    point = (0.23, -0.07, 0.31, 0.37)
    exact = field.spatial_jacobian(*point)
    errors = [
        float(np.linalg.norm(_fd_spatial_jacobian(field, point, step) - exact))
        for step in (0.02, 0.01, 0.005)
    ]

    assert exact.shape == (3, 3)
    assert errors[1] < 0.30 * errors[0]
    assert errors[2] < 0.30 * errors[1]
    assert errors[2] < 5.0e-3


def test_analytic_time_derivative_matches_independent_three_level_fd():
    field = _correction()
    point = (0.23, -0.07, 0.31, 0.37)
    exact = field.time_derivative(*point)
    errors = [
        float(np.linalg.norm(_fd_time_derivative(field, point, step) - exact))
        for step in (0.02, 0.01, 0.005)
    ]

    assert errors[1] < 0.30 * errors[0]
    assert errors[2] < 0.30 * errors[1]
    assert errors[2] < 2.0e-4


def test_analytic_laplacian_matches_independent_three_level_fd():
    field = _correction()
    point = (0.23, -0.07, 0.31, 0.37)
    exact = field.laplacian(*point)
    errors = [
        float(np.linalg.norm(_fd_laplacian(field, point, step) - exact))
        for step in (0.02, 0.01, 0.005)
    ]

    assert errors[1] < 0.30 * errors[0]
    assert errors[2] < 0.30 * errors[1]
    assert errors[2] < 2.0e-2


def test_analytic_vorticity_and_divergence_are_consistent_with_public_velocity():
    field = _correction()
    point = (0.23, -0.07, 0.31, 0.37)
    jacobian = field.spatial_jacobian(*point)
    numerical = _fd_spatial_jacobian(field, point, 0.005)
    numerical_curl = np.asarray(
        [
            numerical[2, 1] - numerical[1, 2],
            numerical[0, 2] - numerical[2, 0],
            numerical[1, 0] - numerical[0, 1],
        ]
    )

    assert np.linalg.norm(field.vorticity(*point) - numerical_curl) < 8.0e-3
    assert abs(field.divergence(*point)) < 2.0e-12
    assert abs(np.trace(jacobian)) < 2.0e-12


def test_support_is_exact_and_batch_evaluation_is_finite():
    field = _correction()
    points = np.asarray(
        [
            [0.23, -0.07, 0.31],
            [0.10, -0.20, 0.15],
            [1.01, -0.20, 0.15],
            [0.10, 0.61, 0.15],
            [0.10, -0.20, 0.86],
        ]
    )
    times = np.asarray([0.37, 0.42, 0.37, 0.37, 0.37])
    values = field.at_points(points, time=times)
    jacobian = field.spatial_jacobian(points[:, 0], points[:, 1], points[:, 2], times)
    time_derivative = field.time_derivative(points[:, 0], points[:, 1], points[:, 2], times)
    vorticity = field.vorticity(points[:, 0], points[:, 1], points[:, 2], times)
    laplacian = field.laplacian(points[:, 0], points[:, 1], points[:, 2], times)

    assert values.shape == (5, 3)
    assert jacobian.shape == (5, 3, 3)
    assert time_derivative.shape == (5, 3)
    assert vorticity.shape == (5, 3)
    assert laplacian.shape == (5, 3)
    for array in (values, jacobian, time_derivative, vorticity, laplacian):
        assert np.all(np.isfinite(array))
    assert np.linalg.norm(values[0]) > 0.0
    assert np.linalg.norm(values[1]) > 0.0
    assert np.array_equal(values[2:], np.zeros((3, 3)))
    assert np.array_equal(jacobian[2:], np.zeros((3, 3, 3)))
    assert np.array_equal(time_derivative[2:], np.zeros((3, 3)))
    assert np.array_equal(vorticity[2:], np.zeros((3, 3)))
    assert np.array_equal(laplacian[2:], np.zeros((3, 3)))


def test_public_supported_velocity_can_be_composed_without_mutating_base_candidate():
    base = Eq45SupportedDeliveryField()
    correction = _correction()
    combined = compose_velocity(base.velocity, correction)

    inside = (0.23, -0.07, 0.31, 0.37)
    outside = (1.2, -0.2, 0.15, 0.37)

    base_inside = base.velocity(*inside)
    combined_inside = combined.velocity(*inside)
    base_outside = base.velocity(*outside)
    combined_outside = combined.velocity(*outside)

    assert base_inside.shape == (3,)
    assert np.all(np.isfinite(combined_inside))
    assert np.linalg.norm(combined_inside - base_inside) > 0.0
    assert np.array_equal(combined_outside, base_outside)
    assert combined.velocity(*inside).shape == base.velocity(*inside).shape


def test_parameter_bounds_and_provenance_truth_boundary():
    field = _correction()
    metadata = field.metadata()

    assert metadata["source"]["edition"] == "2026.09.09-consolidated"
    assert metadata["source"]["doi"] == "10.5281/zenodo.22678406"
    assert metadata["paper_exact"] is False
    assert metadata["openai_field_identified"] is False
    assert "Cartesian affine phase" in metadata["autonomous_choices"]
    assert metadata["analytic_surrogate_derivatives"] == [
        "partial_t velocity",
        "Cartesian spatial Jacobian",
        "vorticity",
        "componentwise Laplacian",
    ]

    with pytest.raises(ValueError):
        KokunoCompleteCurlCorrection(amplitude=2.01)
    with pytest.raises(ValueError):
        KokunoCompleteCurlCorrection(half_widths=(0.01, 0.8, 0.7))
    with pytest.raises(ValueError):
        KokunoCompleteCurlCorrection(wave_vector=(0.0, 0.0, 0.0))
    with pytest.raises(ValueError):
        KokunoCompleteCurlCorrection(wave_vector=(1.0, 0.0, 0.0), polarization=(2.0, 0.0, 0.0))
