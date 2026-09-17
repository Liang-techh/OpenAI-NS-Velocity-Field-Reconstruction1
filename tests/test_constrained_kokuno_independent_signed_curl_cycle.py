import numpy as np

from openai_ns_reconstruction.kokuno_independent_signed_curl_cycle import (
    _mutation_calibration,
    _sample_annulus,
    evaluate_fd2,
)


def _zero_pressure(points, time):
    del time
    return np.zeros(len(points), dtype=float)


def _zero_velocity(points, time):
    del time
    return np.zeros_like(points, dtype=float)


def test_fd2_recovers_manufactured_pressure_gradient():
    points = np.asarray([
        [0.13, -0.08, 0.05],
        [-0.17, 0.09, -0.04],
        [0.11, 0.12, 0.03],
    ])

    def pressure(x, time):
        del time
        return 0.02 * x[:, 0]

    result = evaluate_fd2(_zero_velocity, pressure, points, 0.5, nu=0.01, step=0.001)
    expected = np.column_stack((np.full(len(points), 0.02), np.zeros(len(points)), np.zeros(len(points))))
    assert np.allclose(result["residual"], expected, rtol=0.0, atol=2.0e-14)
    assert np.allclose(result["divergence"], 0.0, rtol=0.0, atol=1.0e-14)


def test_fd2_recovers_linear_velocity_divergence_and_convection():
    points = np.asarray([
        [0.12, 0.02, -0.03],
        [-0.19, 0.06, 0.04],
        [0.15, -0.05, 0.02],
    ])

    def velocity(x, time):
        del time
        values = np.zeros_like(x, dtype=float)
        values[:, 0] = 0.02 * x[:, 0]
        return values

    result = evaluate_fd2(velocity, _zero_pressure, points, 0.5, nu=0.0, step=0.001)
    expected_x = 0.0004 * points[:, 0]
    assert np.allclose(result["divergence"], 0.02, rtol=0.0, atol=2.0e-14)
    assert np.allclose(result["residual"][:, 0], expected_x, rtol=0.0, atol=2.0e-14)
    assert np.allclose(result["residual"][:, 1:], 0.0, rtol=0.0, atol=2.0e-14)


def test_sampling_and_mutation_calibration_are_deterministic_and_sensitive():
    first = _sample_annulus(9173041, 16)
    second = _sample_annulus(9173041, 16)
    assert np.array_equal(first, second)
    radius = np.hypot(first[:, 0], first[:, 1])
    assert np.all((radius >= 0.10) & (radius <= 0.20))
    assert np.all(np.abs(first[:, 2]) <= 0.12)

    calibration = _mutation_calibration(0.001)
    assert calibration["max_divergence_shift_error"] < 1.0e-10
    assert calibration["max_pressure_x_momentum_shift_error"] < 1.0e-10
