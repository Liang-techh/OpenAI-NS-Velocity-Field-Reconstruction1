import math

import numpy as np

from openai_ns_reconstruction.kokuno_native_x_curl import (
    KokunoNativeXCompleteCurlCorrection,
)


def _fd_curl_of_potential(correction, point, time, step):
    point = np.asarray(point, dtype=float)
    jacobian = np.empty((3, 3), dtype=float)
    for axis in range(3):
        shift = np.zeros(3, dtype=float)
        shift[axis] = step
        plus = correction.vector_potential(*(point + shift), time)
        minus = correction.vector_potential(*(point - shift), time)
        jacobian[:, axis] = (plus - minus) / (2.0 * step)
    return np.array(
        [
            jacobian[2, 1] - jacobian[1, 2],
            jacobian[0, 2] - jacobian[2, 0],
            jacobian[1, 0] - jacobian[0, 1],
        ]
    )


def _fd_divergence(correction, point, time, step):
    point = np.asarray(point, dtype=float)
    divergence = 0.0
    for axis in range(3):
        shift = np.zeros(3, dtype=float)
        shift[axis] = step
        plus = correction.velocity(*(point + shift), time)
        minus = correction.velocity(*(point - shift), time)
        divergence += (plus[axis] - minus[axis]) / (2.0 * step)
    return float(divergence)


def test_native_x_support_is_axis_safe_and_compact():
    correction = KokunoNativeXCompleteCurlCorrection()
    time = 0.5
    # At z=0, q=1-t=.5 and therefore X=r^2.  The predeclared support is (.4,.8).
    radii = [0.0, math.sqrt(0.2), math.sqrt(0.6), math.sqrt(0.9)]
    velocity = correction.velocity(np.asarray(radii), 0.0, 0.0, time)
    assert velocity.shape == (4, 3)
    assert np.all(np.isfinite(velocity))
    assert np.array_equal(velocity[0], np.zeros(3))
    assert np.array_equal(velocity[1], np.zeros(3))
    assert np.linalg.norm(velocity[2]) > 1.0e-5
    assert np.array_equal(velocity[3], np.zeros(3))


def test_native_x_window_moves_with_similarity_q():
    correction = KokunoNativeXCompleteCurlCorrection()
    target_X = 0.6
    values = []
    recovered_X = []
    for time in (0.375, 0.5, 0.625):
        q = 1.0 - time  # exact at z=0
        radius = math.sqrt(2.0 * q * target_X)
        X, envelope, gradient = correction.envelope_data(radius, 0.0, 0.0, time)
        recovered_X.append(float(X))
        values.append(float(envelope))
        assert gradient.shape == (3,)
        assert np.all(np.isfinite(gradient))
    assert np.allclose(recovered_X, target_X, rtol=0.0, atol=2.0e-13)
    assert np.allclose(values, 1.0, rtol=0.0, atol=2.0e-13)


def test_analytic_velocity_is_complete_curl_with_fd_convergence():
    correction = KokunoNativeXCompleteCurlCorrection()
    point = np.array([0.75, 0.1, 0.01])
    time = 0.5
    exact = correction.velocity(*point, time)
    errors = []
    for step in (0.02, 0.01, 0.005):
        numerical = _fd_curl_of_potential(correction, point, time, step)
        errors.append(float(np.linalg.norm(numerical - exact)))
    assert errors[2] < errors[1] < errors[0]
    assert errors[0] / errors[1] > 2.5
    assert errors[1] / errors[2] > 2.5
    assert errors[2] < 5.0e-4


def test_divergence_of_complete_curl_converges_to_zero():
    correction = KokunoNativeXCompleteCurlCorrection()
    point = np.array([0.75, 0.1, 0.01])
    time = 0.5
    errors = [abs(_fd_divergence(correction, point, time, h)) for h in (0.02, 0.01, 0.005)]
    assert errors[2] < errors[1] < errors[0]
    assert errors[0] / errors[1] > 2.0
    assert errors[1] / errors[2] > 2.0
    assert errors[2] < 5.0e-4


def test_batch_compose_and_truth_boundary():
    correction = KokunoNativeXCompleteCurlCorrection()
    points = np.array(
        [
            [0.75, 0.10, 0.01],
            [-0.72, 0.18, -0.02],
            [0.00, 0.00, 0.00],
        ]
    )
    direct = correction.at_points(points, 0.5)
    assert direct.shape == (3, 3)
    assert np.all(np.isfinite(direct))

    def base(x, y, z, time):
        x, y, z, time = np.broadcast_arrays(x, y, z, time)
        return np.stack((x * 0.0 + 1.0, y * 0.0 - 2.0, z * 0.0 + time), axis=-1)

    composed = correction.compose_velocity(base)
    result = composed(points[:, 0], points[:, 1], points[:, 2], 0.5)
    expected = base(points[:, 0], points[:, 1], points[:, 2], 0.5) + direct
    assert np.allclose(result, expected, rtol=0.0, atol=1.0e-14)

    payload = correction.to_payload()
    assert payload["source"]["corrected_release"] == "zenodo:22678406"
    assert payload["autonomous_choices"]["X_support"] == [0.4, 0.8]
    assert payload["truth_boundary"]["source_native_localization_coordinate"] is True
    assert payload["truth_boundary"]["source_exact_support"] is False
    assert payload["truth_boundary"]["source_exact_phase_frame"] is False
    assert payload["truth_boundary"]["parameter_selection_performed"] is False
    assert payload["truth_boundary"]["pde_validated"] is False
    assert payload["truth_boundary"]["paper_exact"] is False
