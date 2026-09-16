import numpy as np
import pytest

from openai_ns_reconstruction.constrained_eq45_q import (
    eq45_q_residual,
    solve_eq45_q,
)


def test_exact_h_zero_and_axis_broadcast():
    z = np.array([0.0, 0.2, -0.7, 1.5])
    t = np.array([[0.25], [0.5], [0.75]])
    q = solve_eq45_q(z, t, h=0.0)
    expected = (1.0 - t) + z * z
    np.testing.assert_allclose(q, expected, rtol=0.0, atol=0.0)

    times = np.array([0.1, 0.4, 0.9])
    q_axis = solve_eq45_q(np.zeros(3), times, h=0.005)
    np.testing.assert_allclose(q_axis, 1.0 - times, rtol=0.0, atol=0.0)


def test_h_half_matches_closed_form_and_residual():
    z = np.linspace(-2.0, 2.0, 41)
    t = np.linspace(0.05, 0.95, 41)
    q = solve_eq45_q(z, t, h=0.5)
    tau = 1.0 - t
    expected = 0.5 * (tau + np.sqrt(tau * tau + 4.0 * z * z))
    np.testing.assert_allclose(q, expected, rtol=2e-12, atol=2e-14)
    assert np.max(np.abs(eq45_q_residual(q, z, t, 0.5))) < 2e-12
    assert np.all(q >= tau)


def test_small_h_grid_stability_monotonicity_and_fail_closed():
    z = np.linspace(-3.0, 3.0, 121)[:, None]
    t = np.array([0.25, 0.5, 0.75])[None, :]
    h = 0.005
    q = solve_eq45_q(z, t, h)
    residual = eq45_q_residual(q, z, t, h)
    assert np.max(np.abs(residual)) < 5e-12
    assert np.all(q > 0.0)
    assert np.all(q >= 1.0 - t)

    np.testing.assert_allclose(q, q[::-1], rtol=5e-13, atol=5e-14)
    half = q[60:, :]
    assert np.all(np.diff(half, axis=0) >= -2e-13)

    q_loose = solve_eq45_q(z, t, h, rtol=1e-9, atol=1e-12, max_iter=80)
    np.testing.assert_allclose(q_loose, q, rtol=2e-9, atol=2e-11)

    with pytest.raises(ValueError):
        solve_eq45_q(0.0, 1.0, h)
    with pytest.raises(ValueError):
        solve_eq45_q(0.0, 0.5, -0.1)
    with pytest.raises(ValueError):
        solve_eq45_q(np.nan, 0.5, h)
    with pytest.raises(ValueError):
        solve_eq45_q(0.0, 0.5, h, max_iter=0)
