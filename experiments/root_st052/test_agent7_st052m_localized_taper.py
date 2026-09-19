import numpy as np

import agent7_st052m_localized_taper as m


def test_localized_profile_is_exactly_central_and_endpoint_identity():
    z = np.array([-2.0, -1.8, -0.9, -0.4, 0.0, 0.4, 0.9, 1.8, 2.0])
    q, dq = m.localized_profile_and_derivative(z)
    central = np.abs(z) / m.AXIAL_SUPPORT <= 0.35
    assert np.all(q[central] == 0.0)
    assert np.all(dq[central] == 0.0)
    assert q[0] == 0.0 and q[-1] == 0.0
    assert dq[0] == 0.0 and dq[-1] == 0.0
    assert np.max(q) <= 1.0 + 1e-15


def test_profile_derivative_matches_centered_difference():
    z = np.array([-1.45, -1.30, 1.22, 1.47])
    q, dq = m.localized_profile_and_derivative(z)
    assert np.all(q > 0.0)
    h = 1e-6
    qp, _ = m.localized_profile_and_derivative(z + h)
    qm, _ = m.localized_profile_and_derivative(z - h)
    fd = (qp - qm) / (2.0 * h)
    np.testing.assert_allclose(dq, fd, rtol=2e-7, atol=2e-9)


def test_map_orientation_and_frozen_contract():
    pts = np.array([[0.7, 0.2, 0.0], [0.7, 0.2, 1.30], [0.7, 0.2, 1.95]])
    mapped, a, ap = m.taper_map(pts)
    assert m.PREREG_ISSUE == 542
    assert m.TAPER_TAU == 0.05
    assert m.WINDOW == (0.50, 0.82)
    assert np.all(a > 0.0)
    assert a[0] == 1.0 and ap[0] == 0.0
    np.testing.assert_allclose(mapped[0], pts[0], rtol=0, atol=0)
    assert m.TRUTH["pde_validated"] is False
    assert m.TRUTH["visual_correspondence_verified"] is False
    assert m.TRUTH["material_paths_integrated"] is False
