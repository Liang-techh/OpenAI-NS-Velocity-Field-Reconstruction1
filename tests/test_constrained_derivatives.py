import numpy as np
import pytest

from openai_ns_reconstruction.constrained_derivatives import (
    SupportsAnalyticVelocityDerivatives,
    VelocityDerivativeBundle,
    reference_velocity_derivatives,
)


class ManufacturedDivergenceFree:
    def velocity(self, points, time):
        p = np.asarray(points, dtype=float)
        t = np.asarray(time, dtype=float)
        x, y, z, t = np.broadcast_arrays(*np.moveaxis(p, -1, 0), t)
        return np.stack(
            (
                2.0 * x * y + t * z,
                -x * x + t * x,
                -2.0 * y * z + t * y,
            ),
            axis=-1,
        )


class WrongShape:
    def velocity(self, points, time):
        p = np.asarray(points, dtype=float)
        return np.zeros(p.shape[:-1] + (2,))


def test_reference_bundle_matches_manufactured_first_second_and_time_derivatives():
    field = ManufacturedDivergenceFree()
    points = np.array([[0.2, -0.3, 0.4], [-0.7, 0.25, -0.1]])
    t = np.array([0.4, 0.6])

    d = reference_velocity_derivatives(field, points, t, spatial_step=2e-4, time_step=3e-5)

    x, y, z = points.T
    expected_grad = np.zeros((2, 3, 3))
    expected_grad[:, 0, 0] = 2 * y
    expected_grad[:, 0, 1] = 2 * x
    expected_grad[:, 0, 2] = t
    expected_grad[:, 1, 0] = -2 * x + t
    expected_grad[:, 2, 1] = -2 * z + t
    expected_grad[:, 2, 2] = -2 * y
    np.testing.assert_allclose(d.gradient, expected_grad, atol=2e-10, rtol=0)

    expected_dt = np.stack((z, x, y), axis=-1)
    np.testing.assert_allclose(d.time, expected_dt, atol=2e-10, rtol=0)
    np.testing.assert_allclose(d.divergence, 0.0, atol=3e-10, rtol=0)
    np.testing.assert_allclose(d.laplacian, [[0, -2, 0], [0, -2, 0]], atol=2e-8, rtol=0)

    assert np.allclose(d.hessian[:, 0, 0, 1], 2.0, atol=2e-8, rtol=0)
    assert np.allclose(d.hessian[:, 1, 0, 0], -2.0, atol=2e-8, rtol=0)
    assert np.allclose(d.hessian[:, 2, 1, 2], -2.0, atol=2e-8, rtol=0)
    np.testing.assert_allclose(
        d.hessian, np.swapaxes(d.hessian, -1, -2), atol=0, rtol=0
    )


def test_scalar_point_shape_and_bundle_protocol_surface():
    d = reference_velocity_derivatives(
        ManufacturedDivergenceFree(), [0.1, 0.2, -0.3], 0.5
    )
    assert isinstance(d, VelocityDerivativeBundle)
    assert d.value.shape == (3,)
    assert d.time.shape == (3,)
    assert d.gradient.shape == (3, 3)
    assert d.hessian.shape == (3, 3, 3)
    assert d.laplacian.shape == (3,)
    assert not isinstance(ManufacturedDivergenceFree(), SupportsAnalyticVelocityDerivatives)


def test_invalid_inputs_fail_closed():
    field = ManufacturedDivergenceFree()
    with pytest.raises(ValueError):
        reference_velocity_derivatives(field, [1.0, 2.0], 0.5)
    with pytest.raises(ValueError):
        reference_velocity_derivatives(field, [1.0, 2.0, np.nan], 0.5)
    with pytest.raises(TypeError):
        reference_velocity_derivatives(field, [1.0, 2.0, 3.0], 0.5, spatial_step=True)
    with pytest.raises(ValueError):
        reference_velocity_derivatives(field, [1.0, 2.0, 3.0], 0.5, time_step=0)
    with pytest.raises(ValueError):
        reference_velocity_derivatives(WrongShape(), [1.0, 2.0, 3.0], 0.5)
    with pytest.raises(TypeError):
        reference_velocity_derivatives(object(), [1.0, 2.0, 3.0], 0.5)
