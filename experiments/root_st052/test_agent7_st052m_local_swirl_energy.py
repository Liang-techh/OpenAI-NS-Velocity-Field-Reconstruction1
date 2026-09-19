import numpy as np
import pytest

import agent7_st052m_local_swirl_energy as m


def test_shoulder_bump_is_disjoint_and_central_identity():
    z = np.linspace(-2.0, 2.0, 4001)
    qt, _ = m.failed._bump_and_derivative(z, m.prior.WINDOW)
    qs, dqs = m.shoulder_profile(z)
    assert not np.any((qt > 0.0) & (qs > 0.0))
    central = np.abs(z) / m.failed.AXIAL_SUPPORT <= 0.35
    assert np.all(qs[central] == 0.0)
    assert np.all(dqs[central] == 0.0)
    assert qs[0] == qs[-1] == 0.0


def test_swirl_compensation_changes_only_theta_component():
    pts = np.array(
        [
            [0.70, 0.20, 0.84],
            [0.40, -0.50, -0.90],
            [0.50, 0.10, 0.20],
        ],
        float,
    )
    vel = np.array(
        [
            [-0.20, 0.80, 0.30],
            [-0.40, 0.60, -0.20],
            [-0.10, 0.50, 0.25],
        ],
        float,
    )
    beta = 0.4
    out = m.apply_swirl_compensation(pts, vel, beta)
    q, _ = m.shoulder_profile(pts[:, 2])
    r = np.hypot(pts[:, 0], pts[:, 1])
    er = np.column_stack((pts[:, 0] / r, pts[:, 1] / r, np.zeros(len(pts))))
    et = np.column_stack((-pts[:, 1] / r, pts[:, 0] / r, np.zeros(len(pts))))
    radial0 = np.sum(vel * er, axis=1)
    radial1 = np.sum(out * er, axis=1)
    theta0 = np.sum(vel * et, axis=1)
    theta1 = np.sum(out * et, axis=1)
    np.testing.assert_allclose(radial1, radial0, rtol=0.0, atol=2e-16)
    np.testing.assert_allclose(out[:, 2], vel[:, 2], rtol=0.0, atol=0.0)
    np.testing.assert_allclose(theta1, (1.0 - beta * q) * theta0, rtol=2e-15, atol=2e-15)
    np.testing.assert_array_equal(out[2], vel[2])


def test_frozen_bracket_keeps_swirl_factor_positive_and_truth_boundary():
    z = np.linspace(-2.0, 2.0, 8001)
    q, _ = m.shoulder_profile(z)
    factor = 1.0 - m.BETA_BRACKET[1] * q
    assert np.min(factor) >= 0.1 - 1e-12
    assert np.max(factor) <= 1.0 + 1e-12
    assert m.PREREG_ISSUE == 557
    assert m.TAPER_TAU == 0.05
    assert m.SHOULDER_WINDOW == (0.36, 0.49)
    assert m.BETA_BRACKET == (0.0, 0.9)
    assert m.TRUTH["pde_validated"] is False
    assert m.TRUTH["visual_correspondence_verified"] is False
    assert m.TRUTH["material_paths_integrated"] is False


def test_invalid_beta_is_rejected():
    pts = np.array([[0.5, 0.0, 0.85]])
    vel = np.array([[0.0, 1.0, 0.0]])
    with pytest.raises(ValueError, match="beta outside frozen bracket"):
        m.apply_swirl_compensation(pts, vel, -1e-6)
    with pytest.raises(ValueError, match="beta outside frozen bracket"):
        m.apply_swirl_compensation(pts, vel, 0.900001)
