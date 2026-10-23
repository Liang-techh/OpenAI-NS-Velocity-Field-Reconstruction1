import numpy as np
import pytest

from openai_ns_reconstruction.constrained_poloidal_block import CompactPoloidalBlock


def _block():
    return CompactPoloidalBlock(
        coefficients=(0.20, -0.15, 0.10, 0.05),
        coefficient_bound=0.25,
    )


def test_axis_regular_divergence_free_and_exact_support_zero():
    block = _block()
    rng = np.random.default_rng(7301)
    points = rng.uniform((-1.5, -1.5, -1.5), (1.5, 1.5, 1.5), size=(96, 3))
    points = np.vstack((points, [[0.0, 0.0, 0.0], [0.0, 0.0, 1.1]]))
    out = block.kinematics(points, 0.47)

    assert np.all(np.isfinite(out.velocity))
    assert np.max(np.abs(out.divergence)) < 2.0e-15
    assert np.all(out.velocity[-2:, :2] == 0.0)

    exterior = np.array(
        [[2.0, 0.0, 0.0], [2.1, 0.0, 0.2], [0.3, 0.2, 2.0], [0.3, 0.2, -2.1]]
    )
    assert np.array_equal(block.velocity(exterior, 0.47), np.zeros((4, 3)))


def test_analytic_first_and_time_derivatives_match_independent_differences():
    block = _block()
    points = np.array([[0.35, -0.42, 0.31], [0.71, 0.18, -0.57], [-0.41, 0.63, 0.82]])
    time = 0.51
    out = block.kinematics(points, time)
    h = 2.0e-6

    for axis in range(3):
        delta = np.zeros(3)
        delta[axis] = h
        fd = (block.velocity(points + delta, time) - block.velocity(points - delta, time)) / (2.0 * h)
        assert np.max(np.abs(fd - out.spatial_jacobian[:, :, axis])) < 3.0e-9

    fd_t = (block.velocity(points, time + h) - block.velocity(points, time - h)) / (2.0 * h)
    assert np.max(np.abs(fd_t - out.time_derivative)) < 3.0e-9


def test_parameter_bounds_serialization_and_fail_closed_inputs():
    block = _block()
    rebuilt = CompactPoloidalBlock.from_dict(block.to_dict())
    assert rebuilt == block
    assert rebuilt.parameter_count == 4
    assert rebuilt.is_nontrivial

    with pytest.raises(ValueError):
        CompactPoloidalBlock(coefficients=(0.3, 0.0, 0.0, 0.0), coefficient_bound=0.25)
    with pytest.raises(ValueError):
        CompactPoloidalBlock(coefficients=(np.nan, 0.0, 0.0, 0.0))
    with pytest.raises(ValueError):
        block.velocity(np.zeros((4, 2)), 0.5)
    with pytest.raises(ValueError):
        block.velocity(np.zeros((4, 3)), 0.9)
    with pytest.raises(ValueError):
        CompactPoloidalBlock.from_dict({**block.to_dict(), "schema": "wrong"})
