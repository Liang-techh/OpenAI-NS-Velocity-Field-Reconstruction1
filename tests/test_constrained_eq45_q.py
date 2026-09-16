import numpy as np
import pytest

from openai_ns_reconstruction.constrained_eq45_q import (
    eq45_q_residual,
    solve_eq45_q,
)


def test_quarter_h_matches_closed_form_and_axis_broadcast():
    # For h=1/4, y=sqrt(q) solves y^2-z^2*y-tau=0.
    z = np.linspace(-2.0, 2.0, 41)
    t = np.linspace(0.05, 0.95, 41)
    q = solve_eq45_q(z, t, h=0.25)
    tau = 1.0 - t
    y = 0.5 * (z * z + np.sqrt(z**4 + 4.0 * tau))
    expected = y * y
    np.testing.assert_allclose(q, expected, rtol=3e-12, atol=3e-14)
    assert np.max(np.abs(eq45_q_residual(q, z, t, 0.25))) < 3e-11

    times = np.array([0.1, 0.4, 0.9])
    q_axis = solve_eq45_q(np.zeros(3), times, h=0.005)
    np.testing.assert_allclose(q_axis, 1.0 - times, rtol=0.0, atol=0.0)


def test_small_h_grid_physical_branch_symmetry_and_tolerance_stability():
    z = np.linspace(-3.0, 3.0, 121)[:, None]
    t = np.array([0.25, 0.5, 0.75])[None, :]
    h = 0.005
    q = solve_eq45_q(z, t, h)
    residual = eq45_q_residual(q, z, t, h)
    assert np.max(np.abs(residual)) < 5e-11
    assert np.all(q >= 1.0 - t)
    np.testing.assert_allclose(q, q[::-1], rtol=5e-13, atol=5e-14)

    D = 0.5 - h
    eta = z / np.power(q, D)
    assert np.max(np.abs(eta)) < 1.0

    loose = solve_eq45_q(z, t, h, rtol=1e-9, atol=1e-12)
    assert np.max(np.abs(loose - q)) < 2e-8


def test_fail_closed_source_domain_and_inputs():
    with pytest.raises(ValueError, match="0 < h < 1/2"):
        solve_eq45_q(0.2, 0.5, 0.0)
    with pytest.raises(ValueError, match="0 < h < 1/2"):
        solve_eq45_q(0.2, 0.5, 0.5)
    with pytest.raises(ValueError):
        solve_eq45_q(np.inf, 0.5, 0.005)
    with pytest.raises(ValueError):
        solve_eq45_q(0.2, 1.0, 0.005)
    with pytest.raises(ValueError):
        solve_eq45_q(0.2, 0.5, 0.005, max_iter=0)
