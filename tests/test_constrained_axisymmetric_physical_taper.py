import numpy as np
import pytest

from openai_ns_reconstruction.constrained_axisymmetric_physical_taper import (
    AxisymmetricPhysicalTaper,
)


def _base_cylindrical(points):
    points = np.asarray(points, dtype=float)
    x = points[..., 0]
    y = points[..., 1]
    z = points[..., 2]
    r = np.hypot(x, y)

    shape = 1.0 + 0.2 * z + 0.1 * z * z
    psi = r * r * shape
    u_r = -r * (0.2 + 0.2 * z)
    u_z = 2.0 * shape
    u_theta = 0.7 * r * (1.0 + 0.1 * z)
    return psi, u_r, u_theta, u_z


def _field(points, taper):
    psi, u_r, u_theta, u_z = _base_cylindrical(points)
    return taper.apply_cylindrical(
        points, psi=psi, u_r=u_r, u_theta=u_theta, u_z=u_z
    )


def _base_cartesian(points):
    points = np.asarray(points, dtype=float)
    x = points[..., 0]
    y = points[..., 1]
    r = np.hypot(x, y)
    _, u_r, u_theta, u_z = _base_cylindrical(points)
    cos_theta = np.divide(x, r, out=np.ones_like(r), where=r > 0.0)
    sin_theta = np.divide(y, r, out=np.zeros_like(r), where=r > 0.0)
    u = u_r * cos_theta - u_theta * sin_theta
    v = u_r * sin_theta + u_theta * cos_theta
    u = np.where(r == 0.0, 0.0, u)
    v = np.where(r == 0.0, 0.0, v)
    return np.stack((u, v, u_z), axis=-1)


def test_inner_plateau_is_exact_and_axis_is_regular():
    taper = AxisymmetricPhysicalTaper()
    points = np.array([
        [0.5, 0.4, 0.5],
        [-0.8, 0.3, -0.7],
        [0.0, 0.0, 0.5],
    ])
    got = _field(points, taper)
    np.testing.assert_allclose(got, _base_cartesian(points), rtol=0.0, atol=0.0)
    assert np.all(np.isfinite(got))
    assert got[-1, 0] == 0.0
    assert got[-1, 1] == 0.0


def test_taper_is_exactly_zero_at_and_outside_physical_support():
    taper = AxisymmetricPhysicalTaper(radial_support=2.0, axial_half_height=2.0)
    points = np.array([
        [2.0, 0.0, 0.0],
        [2.1, 0.0, 0.0],
        [1.5, 1.5, 0.0],
        [0.0, 0.0, 2.0],
        [0.0, 0.0, -2.2],
        [2.0, 0.0, 2.0],
    ])
    np.testing.assert_allclose(_field(points, taper), 0.0, rtol=0.0, atol=0.0)


def test_streamfunction_taper_keeps_divergence_small_through_transition():
    taper = AxisymmetricPhysicalTaper()
    rng = np.random.default_rng(20260916)
    points = rng.uniform(-1.85, 1.85, size=(96, 3))

    step = 1e-5
    divergence = np.zeros(points.shape[0])
    for axis in range(3):
        delta = np.zeros(3)
        delta[axis] = step
        plus = _field(points + delta, taper)
        minus = _field(points - delta, taper)
        divergence += (plus[:, axis] - minus[:, axis]) / (2.0 * step)

    assert np.max(np.abs(divergence)) < 1e-6
    assert np.sqrt(np.mean(divergence * divergence)) < 2e-7


def test_serialization_and_fail_closed_axis_regular_contract():
    taper = AxisymmetricPhysicalTaper(
        radial_plateau_q=0.7,
        axial_plateau_q=0.75,
    )
    restored = AxisymmetricPhysicalTaper.from_dict(taper.to_dict())
    assert restored == taper
    payload = taper.to_dict()
    assert payload["classification"] == "autonomous_design"
    assert payload["pde_validated"] is False
    assert payload["paper_exact"] is False

    with pytest.raises(ValueError, match="streamfunction"):
        taper.apply_cylindrical(
            np.array([[0.0, 0.0, 0.2]]),
            psi=np.array([1e-3]),
            u_r=np.array([0.0]),
            u_theta=np.array([0.0]),
            u_z=np.array([0.0]),
        )

    with pytest.raises(ValueError, match="plateau"):
        AxisymmetricPhysicalTaper(radial_plateau_q=1.0)
    with pytest.raises(ValueError, match="shape"):
        taper.factors(np.zeros((4, 2)))
