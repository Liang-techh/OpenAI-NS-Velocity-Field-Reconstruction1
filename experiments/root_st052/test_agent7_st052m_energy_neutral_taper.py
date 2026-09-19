import numpy as np

import agent7_st052m_energy_neutral_taper as m


def test_two_bumps_are_disjoint_and_central_identity():
    z = np.linspace(-2.0, 2.0, 4001)
    qt, dqt = m._bump_and_derivative(z, m.TIP_WINDOW)
    qs, dqs = m._bump_and_derivative(z, m.SHOULDER_WINDOW)
    assert not np.any((qt > 0.0) & (qs > 0.0))
    central = np.abs(z) / m.AXIAL_SUPPORT <= 0.35
    assert np.all(qt[central] == 0.0)
    assert np.all(qs[central] == 0.0)
    assert np.all(dqt[central] == 0.0)
    assert np.all(dqs[central] == 0.0)
    assert qt[0] == qt[-1] == 0.0
    assert qs[0] == qs[-1] == 0.0


def test_bump_derivatives_match_centered_difference():
    h = 1.0e-6
    for window, z in [
        (m.TIP_WINDOW, np.array([-1.45, -1.30, 1.22, 1.47])),
        (m.SHOULDER_WINDOW, np.array([-0.95, -0.82, 0.80, 0.94])),
    ]:
        q, dq = m._bump_and_derivative(z, window)
        assert np.all(q > 0.0)
        qp, _ = m._bump_and_derivative(z + h, window)
        qm, _ = m._bump_and_derivative(z - h, window)
        fd = (qp - qm) / (2.0 * h)
        np.testing.assert_allclose(dq, fd, rtol=3e-7, atol=3e-9)


def test_frozen_bracket_preserves_orientation_and_contract():
    z = np.linspace(-2.0, 2.0, 8001)
    for alpha in m.ENERGY_ALPHA_BRACKET:
        q, _ = m.neutral_profile_and_derivative(z, alpha)
        a = 1.0 + m.TAPER_TAU * q
        assert np.min(a) >= 0.8 - 1e-12
        assert np.max(a) <= 1.05 + 1e-12
    assert m.PREREG_ISSUE == 548
    assert m.TAPER_TAU == 0.05
    assert m.TIP_WINDOW == (0.50, 0.82)
    assert m.SHOULDER_WINDOW == (0.36, 0.49)
    assert m.ENERGY_ALPHA_BRACKET == (0.0, 4.0)
    assert m.TRUTH["pde_validated"] is False
    assert m.TRUTH["visual_correspondence_verified"] is False
    assert m.TRUTH["material_paths_integrated"] is False


def test_neutral_map_is_exact_identity_in_registered_center():
    pts = np.array([
        [0.7, 0.2, 0.0],
        [0.4, -0.3, 0.3],
        [0.5, 0.1, -0.6],
    ])
    mapped, a, ap = m.neutral_map(pts, alpha_e=2.0)
    np.testing.assert_array_equal(mapped, pts)
    np.testing.assert_array_equal(a, np.ones(len(pts)))
    np.testing.assert_array_equal(ap, np.zeros(len(pts)))
