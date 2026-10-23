import json

import numpy as np

from openai_ns_reconstruction.constrained_eq45_structure_receipt import (
    build_eq45_structure_receipt,
    write_eq45_structure_receipt,
)


def test_symbolic_receipt_verifies_four_local_identities(tmp_path):
    receipt = build_eq45_structure_receipt()
    assert receipt["status"] == "four_strict_local_identities_verified"
    assert len(receipt["identities"]) == 4
    assert receipt["states"]["pde_validated"] is False
    assert receipt["states"]["paper_exact"] is False
    path = tmp_path / "receipt.json"
    write_eq45_structure_receipt(path)
    assert json.loads(path.read_text())["schema"] == "eq45_structure_receipt_v1"


def _mix(points, q, h, profile):
    points = np.asarray(points, dtype=float)
    q = np.asarray(q, dtype=float)
    v0, F, U = (np.asarray(profile[..., i], dtype=float) for i in range(3))
    x, y = points[..., 0], points[..., 1]
    ss = q ** (-1.0 - h)
    aa = q ** (-0.5 - h)
    return np.stack((x*v0/(2*q)-y*ss*F, y*v0/(2*q)+x*ss*F, aa*U), axis=-1)


def test_cartesian_mixing_matches_cylindrical_projection_without_axis_division():
    rng = np.random.default_rng(914071)
    points = rng.normal(size=(128, 3))
    points[:4, :2] = 0.0
    q = rng.uniform(0.3, 1.4, size=128)
    profile = rng.normal(size=(128, 3))
    h = 0.005
    vel = _mix(points, q, h, profile)
    x, y = points[:, 0], points[:, 1]
    r2 = x*x + y*y
    np.testing.assert_allclose(x*vel[:,0] + y*vel[:,1], r2*profile[:,0]/(2*q), rtol=2e-14, atol=2e-14)
    np.testing.assert_allclose(-y*vel[:,0] + x*vel[:,1], r2*q**(-1-h)*profile[:,1], rtol=2e-14, atol=2e-14)
    np.testing.assert_allclose(vel[:4,:2], 0.0, atol=0.0)


def test_rotation_equivariance_for_rotation_invariant_profiles():
    rng = np.random.default_rng(914073)
    points = rng.normal(size=(64, 3))
    q = rng.uniform(0.4, 1.2, size=64)
    profile = rng.normal(size=(64, 3))
    h = 0.005
    angle = 0.731
    c, s = np.cos(angle), np.sin(angle)
    rot = np.array([[c,-s,0.0],[s,c,0.0],[0.0,0.0,1.0]])
    rotated_points = points @ rot.T
    base = _mix(points, q, h, profile)
    rotated = _mix(rotated_points, q, h, profile)
    np.testing.assert_allclose(rotated, base @ rot.T, rtol=2e-14, atol=2e-14)
