from pathlib import Path

import numpy as np
import pytest

from openai_ns_reconstruction.constrained_axisymmetric_physical_taper import (
    AxisymmetricPhysicalTaper,
)
from openai_ns_reconstruction.constrained_eq45_candidate import Eq45VelocityCandidate
from openai_ns_reconstruction.constrained_eq45_supported_candidate import (
    Eq45SupportedVelocityCandidate,
)


ROOT = Path(__file__).resolve().parents[1]
SEED = ROOT / "artifacts" / "constrained" / "eq45_velocity_candidate_seed.json"


def _parent() -> Eq45VelocityCandidate:
    return Eq45VelocityCandidate.load_json(SEED)


def _supported() -> Eq45SupportedVelocityCandidate:
    return Eq45SupportedVelocityCandidate(parent=_parent())


def test_supported_candidate_preserves_plateau_and_zeroes_physical_exterior():
    parent = _parent()
    child = Eq45SupportedVelocityCandidate(parent=parent)

    # Default q=0.64 plateau means r<=1.6 and |z|<=1.6 are exactly untouched.
    plateau = np.array(
        [
            [0.0, 0.0, 0.0],
            [0.25, -0.1, 0.3],
            [0.8, 0.6, -0.4],
            [-1.0, 0.5, 0.75],
            [1.45, 0.2, -1.2],
        ]
    )
    times = np.array([0.25, 0.5, 0.75, 0.5, 0.25])
    np.testing.assert_allclose(
        child.at_points(plateau, times),
        parent.at_points(plateau, times),
        rtol=3e-14,
        atol=3e-14,
    )

    exterior = np.array(
        [
            [2.0, 0.0, 0.0],
            [2.15, 0.0, 0.2],
            [0.0, -2.0, -0.1],
            [0.0, 0.0, 2.0],
            [0.25, 0.2, -2.25],
            [1.6, 1.6, 0.5],
        ]
    )
    zero = child.at_points(exterior, 0.5)
    np.testing.assert_allclose(zero, 0.0, rtol=0.0, atol=0.0)

    # The exterior connection is a genuine non-identity child in the collar.
    collar = np.array(
        [
            [1.75, 0.0, 0.25],
            [1.65, 0.25, -0.5],
            [0.6, 0.2, 1.75],
        ]
    )
    difference = child.at_points(collar, 0.5) - parent.at_points(collar, 0.5)
    assert np.max(np.linalg.norm(difference, axis=-1)) > 1e-6


def test_supported_candidate_public_grid_serialization_and_child_identity(tmp_path):
    parent = _parent()
    child = Eq45SupportedVelocityCandidate(parent=parent)

    assert child.parent_sha256 == parent.sha256
    assert child.sha256 != parent.sha256
    payload = child.to_dict()
    truth = payload["truth_boundary"]
    assert truth["velocity_export_ready"] is True
    assert truth["physical_support_connection_implemented"] is True
    assert truth["requires_physical_support_connection"] is False
    assert truth["physical_support_validated"] is False
    assert truth["visualization_ready"] is False
    assert truth["visual_correspondence_verified"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False
    assert truth["openai_field_identified"] is False

    target = tmp_path / "supported.json"
    child.save_json(target)
    loaded = Eq45SupportedVelocityCandidate.load_json(target)
    assert loaded.sha256 == child.sha256
    assert loaded.parent_sha256 == parent.sha256

    points = np.array([[0.3, 0.2, -0.4], [1.75, 0.0, 0.3], [2.0, 0.0, 0.0]])
    np.testing.assert_allclose(
        loaded.at_points(points, np.array([0.25, 0.5, 0.75])),
        child.at_points(points, np.array([0.25, 0.5, 0.75])),
        rtol=0.0,
        atol=0.0,
    )

    axis = np.linspace(-2.0, 2.0, 5)
    grid = loaded.grid(axis, axis, axis, np.array([0.25, 0.75]))
    assert grid.shape == (2, 5, 5, 5, 3)
    assert np.all(np.isfinite(grid))
    # Every r=2 or |z|=2 face entry is exactly zero under the physical connection.
    np.testing.assert_allclose(grid[:, 0, :, :, :], 0.0, rtol=0.0, atol=0.0)
    np.testing.assert_allclose(grid[:, -1, :, :, :], 0.0, rtol=0.0, atol=0.0)
    np.testing.assert_allclose(grid[:, :, 0, :, :], 0.0, rtol=0.0, atol=0.0)
    np.testing.assert_allclose(grid[:, :, -1, :, :], 0.0, rtol=0.0, atol=0.0)
    np.testing.assert_allclose(grid[:, :, :, 0, :], 0.0, rtol=0.0, atol=0.0)
    np.testing.assert_allclose(grid[:, :, :, -1, :], 0.0, rtol=0.0, atol=0.0)


def _centered_divergence(candidate, points, time, step=2e-5):
    points = np.asarray(points, dtype=float)
    divergence = np.zeros(points.shape[0], dtype=float)
    for axis in range(3):
        plus = points.copy()
        minus = points.copy()
        plus[:, axis] += step
        minus[:, axis] -= step
        values_plus = candidate.at_points(plus, time)
        values_minus = candidate.at_points(minus, time)
        divergence += (values_plus[:, axis] - values_minus[:, axis]) / (2.0 * step)
    return divergence


def test_supported_candidate_production_collar_retains_numerical_incompressibility():
    child = _supported()
    # Off-axis probes exercise radial and axial taper collars without touching a
    # plateau/support interface, where centered finite differences are cleanest.
    points = np.array(
        [
            [1.70, 0.20, 0.20],
            [1.82, -0.15, 0.45],
            [1.65, 0.35, -0.65],
            [0.80, 0.25, 1.70],
            [0.55, -0.35, 1.82],
            [0.95, 0.20, -1.72],
            [1.68, 0.18, 1.68],
            [1.76, -0.12, -1.70],
        ]
    )
    divergence = _centered_divergence(child, points, 0.5)
    assert float(np.sqrt(np.mean(divergence * divergence))) < 2e-5
    assert float(np.max(np.abs(divergence))) < 5e-5


def test_supported_candidate_fails_closed_on_identity_and_metadata_tampering():
    parent = _parent()
    child = Eq45SupportedVelocityCandidate(parent=parent)

    with pytest.raises(TypeError, match="Eq45VelocityCandidate"):
        Eq45SupportedVelocityCandidate(parent=object())
    with pytest.raises(TypeError, match="AxisymmetricPhysicalTaper"):
        Eq45SupportedVelocityCandidate(parent=parent, taper=object())
    with pytest.raises(ValueError, match="declared interval"):
        child.at_points(np.zeros(3), 0.9)

    payload = child.to_dict()
    payload["parent_sha256"] = "0" * 64
    with pytest.raises(ValueError, match="parent_sha256"):
        Eq45SupportedVelocityCandidate.from_dict(payload)

    payload = child.to_dict()
    payload["truth_boundary"]["pde_validated"] = True
    with pytest.raises(ValueError, match="truth-boundary"):
        Eq45SupportedVelocityCandidate.from_dict(payload)
