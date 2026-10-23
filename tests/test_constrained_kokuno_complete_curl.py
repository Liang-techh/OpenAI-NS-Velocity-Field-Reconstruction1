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
    values = field.at_points(points, time=np.asarray([0.37, 0.42, 0.37, 0.37, 0.37]))

    assert values.shape == (5, 3)
    assert np.all(np.isfinite(values))
    assert np.linalg.norm(values[0]) > 0.0
    assert np.linalg.norm(values[1]) > 0.0
    assert np.array_equal(values[2:], np.zeros((3, 3)))


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

    with pytest.raises(ValueError):
        KokunoCompleteCurlCorrection(amplitude=2.01)
    with pytest.raises(ValueError):
        KokunoCompleteCurlCorrection(half_widths=(0.01, 0.8, 0.7))
    with pytest.raises(ValueError):
        KokunoCompleteCurlCorrection(wave_vector=(0.0, 0.0, 0.0))
    with pytest.raises(ValueError):
        KokunoCompleteCurlCorrection(wave_vector=(1.0, 0.0, 0.0), polarization=(2.0, 0.0, 0.0))
