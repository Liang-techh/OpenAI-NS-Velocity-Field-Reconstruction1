from pathlib import Path

import numpy as np
import pytest

from openai_ns_reconstruction.constrained_eq45_candidate import Eq45VelocityCandidate
from openai_ns_reconstruction.constrained_eq45_cylindrical_state import (
    Eq45CylindricalState,
    eq45_candidate_cylindrical_state,
)


ROOT = Path(__file__).resolve().parents[1]
SEED = ROOT / "artifacts" / "constrained" / "eq45_velocity_candidate_seed.json"


def _seed() -> Eq45VelocityCandidate:
    return Eq45VelocityCandidate.load_json(SEED)


def test_cylindrical_state_roundtrips_public_velocity_and_is_axis_regular():
    candidate = _seed()
    points = np.array(
        [
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.25],
            [0.35, -0.2, 0.15],
            [-0.45, 0.3, -0.2],
            [0.7, 0.1, 0.35],
        ]
    )
    times = np.array([0.25, 0.5, 0.75, 0.5, 0.25])

    state = eq45_candidate_cylindrical_state(candidate, points, times)
    expected = candidate.at_points(points, times)
    reconstructed = state.cartesian(points)

    np.testing.assert_allclose(reconstructed, expected, rtol=2e-14, atol=2e-14)
    assert state.q.shape == points.shape[:-1]
    assert np.all(state.q > 0.0)
    assert np.all(state.X >= 0.0)

    axis = np.hypot(points[:, 0], points[:, 1]) == 0.0
    np.testing.assert_allclose(state.psi[axis], 0.0, rtol=0.0, atol=0.0)
    np.testing.assert_allclose(state.u_r[axis], 0.0, rtol=0.0, atol=0.0)
    np.testing.assert_allclose(state.u_theta[axis], 0.0, rtol=0.0, atol=0.0)
    assert np.all(np.isfinite(state.u_z[axis]))

    taper_inputs = state.taper_inputs()
    assert set(taper_inputs) == {"psi", "u_r", "u_theta", "u_z"}
    for name in taper_inputs:
        np.testing.assert_array_equal(taper_inputs[name], getattr(state, name))


def test_cylindrical_streamfunction_matches_candidate_poloidal_derivatives():
    candidate = _seed()
    point = np.array([0.55, 0.0, 0.2])
    time = 0.5
    step = 2e-6

    center = eq45_candidate_cylindrical_state(candidate, point, time)

    radial_plus = point.copy()
    radial_minus = point.copy()
    radial_plus[0] += step
    radial_minus[0] -= step
    psi_r_plus = eq45_candidate_cylindrical_state(candidate, radial_plus, time).psi
    psi_r_minus = eq45_candidate_cylindrical_state(candidate, radial_minus, time).psi
    dpsi_dr = (psi_r_plus - psi_r_minus) / (2.0 * step)

    axial_plus = point.copy()
    axial_minus = point.copy()
    axial_plus[2] += step
    axial_minus[2] -= step
    psi_z_plus = eq45_candidate_cylindrical_state(candidate, axial_plus, time).psi
    psi_z_minus = eq45_candidate_cylindrical_state(candidate, axial_minus, time).psi
    dpsi_dz = (psi_z_plus - psi_z_minus) / (2.0 * step)

    radius = point[0]
    assert center.u_z == pytest.approx(dpsi_dr / radius, rel=3e-6, abs=3e-8)
    assert center.u_r == pytest.approx(-dpsi_dz / radius, rel=3e-6, abs=3e-8)


def test_cylindrical_adapter_preserves_candidate_identity_and_truth_boundary():
    candidate = _seed()
    before_sha = candidate.sha256
    before_payload = candidate.to_dict()

    points = np.array([[0.25, 0.15, -0.2], [0.4, -0.3, 0.25]])
    state = eq45_candidate_cylindrical_state(candidate, points, 0.5)
    assert isinstance(state, Eq45CylindricalState)

    assert candidate.sha256 == before_sha
    assert candidate.to_dict() == before_payload
    truth = candidate.to_dict()["truth_boundary"]
    assert truth["velocity_export_ready"] is True
    assert truth["physical_support_validated"] is False
    assert truth["visualization_ready"] is False
    assert truth["visual_correspondence_verified"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False
    assert truth["openai_field_identified"] is False


def test_cylindrical_adapter_fails_closed_on_invalid_inputs():
    candidate = _seed()

    with pytest.raises(TypeError, match="Eq45VelocityCandidate"):
        eq45_candidate_cylindrical_state(object(), np.zeros(3), 0.5)
    with pytest.raises(ValueError, match="declared interval"):
        eq45_candidate_cylindrical_state(candidate, np.zeros(3), 0.9)

    state = eq45_candidate_cylindrical_state(candidate, np.zeros(3), 0.5)
    with pytest.raises(ValueError, match="shape"):
        state.cartesian(np.zeros((2, 3)))
    with pytest.raises(ValueError, match="axis_tolerance"):
        state.cartesian(np.zeros(3), axis_tolerance=-1.0)
