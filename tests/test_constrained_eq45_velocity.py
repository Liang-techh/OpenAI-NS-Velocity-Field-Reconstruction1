import numpy as np
import pytest

from openai_ns_reconstruction.constrained_eq45_q import solve_eq45_q
from openai_ns_reconstruction.constrained_eq45_velocity import Eq45VelocityBackbone


def _constant_profiles(X, eta):
    shape = np.shape(X)
    return np.stack(
        (
            np.ones(shape),
            2.0 * np.ones(shape),
            3.0 * np.ones(shape),
        ),
        axis=-1,
    )


def test_vectorized_formula_matches_direct_cartesian_mixing():
    field = Eq45VelocityBackbone(_constant_profiles, h=0.25)
    points = np.array(
        [
            [0.2, -0.3, -0.4],
            [0.5, 0.1, 0.0],
            [-0.4, 0.6, 0.7],
        ]
    )
    times = np.array([0.2, 0.5, 0.7])

    got = field.velocity(points, times)
    q = solve_eq45_q(points[:, 2], times, 0.25)
    x = points[:, 0]
    y = points[:, 1]
    expected = np.stack(
        (
            x / (2.0 * q) - y * q ** (-1.25) * 2.0,
            y / (2.0 * q) + x * q ** (-1.25) * 2.0,
            q ** (-0.75) * 3.0,
        ),
        axis=-1,
    )
    np.testing.assert_allclose(got, expected, rtol=2e-13, atol=2e-13)


def test_coordinates_satisfy_source_identities_and_axis_is_regular():
    field = Eq45VelocityBackbone(_constant_profiles, h=0.005)
    points = np.array(
        [
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.8],
            [0.3, -0.2, -0.7],
        ]
    )
    time = 0.4

    coordinates = field.coordinates(points, time)
    D = 0.5 - field.h
    np.testing.assert_allclose(
        points[:, 2],
        coordinates.q**D * coordinates.eta,
        rtol=2e-12,
        atol=2e-12,
    )
    np.testing.assert_allclose(
        1.0 - time,
        coordinates.q * (1.0 - coordinates.eta**2),
        rtol=2e-11,
        atol=2e-11,
    )

    velocity = field.velocity(points, time)
    np.testing.assert_allclose(velocity[:2, :2], 0.0, rtol=0.0, atol=0.0)
    assert np.all(np.isfinite(velocity))


def test_fail_closed_on_bad_profiles_inputs_and_time_broadcast():
    bad_shape = Eq45VelocityBackbone(
        lambda X, eta: np.zeros(np.shape(X) + (2,)), h=0.01
    )
    with pytest.raises(ValueError, match="profile_provider"):
        bad_shape.velocity(np.zeros((4, 3)), 0.5)

    nonfinite = Eq45VelocityBackbone(
        lambda X, eta: np.full(np.shape(X) + (3,), np.nan), h=0.01
    )
    with pytest.raises(ValueError, match="nonfinite"):
        nonfinite.velocity(np.zeros((2, 3)), 0.5)

    good = Eq45VelocityBackbone(_constant_profiles, h=0.01)
    with pytest.raises(ValueError, match="broadcast"):
        good.velocity(np.zeros((2, 3)), np.zeros(3))
    with pytest.raises(ValueError, match="0 < h < 1/2"):
        Eq45VelocityBackbone(_constant_profiles, h=0.5)
