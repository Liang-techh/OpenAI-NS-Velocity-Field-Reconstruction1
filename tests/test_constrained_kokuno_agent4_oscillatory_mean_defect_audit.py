import math

import numpy as np

from openai_ns_reconstruction.kokuno_agent4_oscillatory_mean_defect_audit import (
    AGENT3_RECEIPT_TARGET,
    ANGULAR_ORDER,
    FD6_STEPS,
    GUARDS,
    RADIAL_COUNT,
    _fd6_spatial_operator,
    _ring_cloud,
)


def _manufactured_velocity(x, y, z, t):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    z = np.asarray(z, dtype=float)
    t = np.broadcast_to(np.asarray(t, dtype=float), x.shape)
    return np.stack((x**4 + t, y**4 + 2.0 * t, z**4 - t), axis=-1)


def _manufactured_dt(x, y, z, t):
    x = np.asarray(x, dtype=float)
    return np.broadcast_to(np.asarray((1.0, 2.0, -1.0), dtype=float), x.shape + (3,)).copy()


def test_independent_fd6_operator_matches_polynomial_manufactured_solution():
    x = np.asarray((0.21, -0.37, 0.48), dtype=float)
    y = np.asarray((-0.31, 0.42, 0.18), dtype=float)
    z = np.asarray((0.27, -0.16, 0.39), dtype=float)
    t = np.asarray((0.35, 0.51, 0.63), dtype=float)
    nu = 0.01
    result = _fd6_spatial_operator(
        _manufactured_velocity,
        _manufactured_dt,
        x,
        y,
        z,
        t,
        step=0.013,
        nu=nu,
    )
    u = _manufactured_velocity(x, y, z, t)
    dt = _manufactured_dt(x, y, z, t)
    exact_gradient = np.zeros((x.size, 3, 3), dtype=float)
    exact_gradient[:, 0, 0] = 4.0 * x**3
    exact_gradient[:, 1, 1] = 4.0 * y**3
    exact_gradient[:, 2, 2] = 4.0 * z**3
    exact_laplacian = np.stack((12.0 * x**2, 12.0 * y**2, 12.0 * z**2), axis=-1)
    exact_quadratic = np.einsum("...j,...jk->...k", u, exact_gradient)
    exact_raw = dt + exact_quadratic - nu * exact_laplacian

    np.testing.assert_allclose(result["gradient"], exact_gradient, rtol=2e-11, atol=2e-11)
    np.testing.assert_allclose(result["laplacian"], exact_laplacian, rtol=2e-9, atol=2e-9)
    np.testing.assert_allclose(result["raw_self_residual"], exact_raw, rtol=2e-10, atol=2e-10)


def test_independent_ring_quadrature_and_fd6_stencils_are_inside_support():
    radii, theta, weights, x, y, z, t = _ring_cloud()
    assert radii.shape == (RADIAL_COUNT - 2,)
    assert theta.shape == weights.shape == (ANGULAR_ORDER,)
    assert abs(float(np.sum(weights)) - 1.0) <= 64.0 * math.ulp(1.0)
    assert x.shape == y.shape == z.shape == t.shape == ((RADIAL_COUNT - 2) * ANGULAR_ORDER,)
    assert abs(float(radii[0]) - 0.2) < 1e-14
    assert abs(float(radii[-1]) - 1.3) < 1e-14
    assert 3.0 * max(FD6_STEPS) < min(radii[0] - 0.15, 1.35 - radii[-1])


def test_receipt_targets_and_scientific_guards_remain_frozen():
    assert AGENT3_RECEIPT_TARGET["raw_mean_vector_rms"] == 96194924.96271932
    assert AGENT3_RECEIPT_TARGET["quadratic_mean_vector_rms"] == 96194924.96273328
    assert tuple(FD6_STEPS) == (0.01, 0.005, 0.0025)
    assert GUARDS["finest_receipt_raw_mean_rms_relative"] == 5e-3
    assert GUARDS["finest_fd6_profile_change"] == 2e-4
    assert GUARDS["minimum_fd6_profile_refinement_ratio"] == 8.0
    assert GUARDS["minimum_quadratic_sign_flip_response"] == 1.0
