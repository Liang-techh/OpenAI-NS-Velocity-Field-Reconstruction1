import numpy as np
import pytest

from openai_ns_reconstruction.constrained_eq45_channel_capacity import (
    analytic_eq45_channel_norms,
    diagnose_eq45_channel_capacity,
    eq45_profile_mixing_matrix,
)


def test_mixing_matrix_reproduces_eq45_velocity_and_orthogonal_channels():
    x, y, q, h = 0.3, -0.4, 0.6, 0.005
    v0, F, U = 1.2, -0.7, 0.9
    matrix = eq45_profile_mixing_matrix(x, y, q, h)
    got = matrix @ np.array([v0, F, U])
    expected = np.array([
        x / (2*q) * v0 - y * q**(-1-h) * F,
        y / (2*q) * v0 + x * q**(-1-h) * F,
        q**(-0.5-h) * U,
    ])
    np.testing.assert_allclose(got, expected, rtol=2e-15, atol=2e-15)

    report = diagnose_eq45_channel_capacity(x, y, q, h)
    np.testing.assert_allclose(
        report.channel_norms,
        analytic_eq45_channel_norms(x, y, q, h),
        rtol=2e-15,
        atol=2e-15,
    )
    assert report.numerical_rank == 3
    assert report.max_normalized_cross_dot < 2e-15


def test_axis_rank_loss_is_structural_and_near_axis_conditioning_is_quantified():
    q, h = 0.5, 0.005
    axis = diagnose_eq45_channel_capacity(0.0, 0.0, q, h)
    assert axis.numerical_rank == 1
    assert np.isinf(axis.raw_condition_number)
    np.testing.assert_allclose(axis.channel_norms[:2], 0.0, atol=0.0)
    assert axis.channel_norms[2] > 0.0

    radii = np.array([0.01, 0.05, 0.25, 0.5])
    off_axis = diagnose_eq45_channel_capacity(radii, 0.0, q, h)
    assert np.all(off_axis.numerical_rank == 3)
    assert np.all(np.diff(off_axis.raw_condition_number) < 0.0)
    np.testing.assert_allclose(
        off_axis.channel_norms[:, 1] / off_axis.channel_norms[:, 0],
        2.0 * q**(-h),
        rtol=2e-15,
    )
    assert off_axis.raw_condition_number[0] > 100.0
    assert off_axis.raw_condition_number[-1] < 3.0


def test_vectorization_and_fail_closed_inputs():
    report = diagnose_eq45_channel_capacity(
        np.array([0.2, 0.4]),
        np.array([0.1, -0.2]),
        np.array([0.4, 0.7]),
        0.005,
    )
    assert report.singular_values.shape == (2, 3)
    assert report.channel_norms.shape == (2, 3)
    assert np.all(report.numerical_rank == 3)

    with pytest.raises(ValueError, match="strictly positive"):
        diagnose_eq45_channel_capacity(0.1, 0.2, 0.0, 0.005)
    with pytest.raises(ValueError, match="finite"):
        diagnose_eq45_channel_capacity(np.nan, 0.2, 0.5, 0.005)
    with pytest.raises(ValueError, match="rank_rtol"):
        diagnose_eq45_channel_capacity(0.1, 0.2, 0.5, 0.005, rank_rtol=0.0)
