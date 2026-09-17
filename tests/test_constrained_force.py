import numpy as np
import pytest

from openai_ns_reconstruction.constrained_force import RestrictedForce


def test_force_agrees_with_independent_potential_curl():
    force = RestrictedForce(a=1.3, c=0.7)
    points = np.random.default_rng(78).uniform(-1.7, 1.7, (40, 3))
    points = np.vstack((points, [0, 0, 0.4], [0, 0, 0]))
    h = 1e-5
    jac = np.stack([
        (force.potential(points + np.eye(3)[j]*h, 0.4)
         - force.potential(points - np.eye(3)[j]*h, 0.4))/(2*h)
        for j in range(3)
    ], axis=-1)
    curl = np.stack((jac[:, 2, 1]-jac[:, 1, 2],
                     jac[:, 0, 2]-jac[:, 2, 0],
                     jac[:, 1, 0]-jac[:, 0, 1]), axis=-1)
    np.testing.assert_allclose(force(points, 0.4), curl, atol=2e-8, rtol=2e-7)


def test_support_axis_and_coefficient_restrictions():
    f = RestrictedForce()
    np.testing.assert_array_equal(f([[2, 0, 0], [0, 0, 2], [3, 3, 3]], 0.5), 0)
    np.testing.assert_array_equal(f([0.1, 0.1, 0.1], [0, 1, -1, 2]), 0)
    assert np.isfinite(f([0, 0, 0.2], 0.5)).all()
    x, y, z = f([0.1, 0, 0.1], 0.5)
    assert x < 0 < y and z > 0
    for value in (-1, 11, np.nan, np.inf):
        with pytest.raises(ValueError):
            RestrictedForce(a=value)
