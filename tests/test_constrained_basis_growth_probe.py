import importlib.util
from pathlib import Path
import sys

import numpy as np
import pytest

from openai_ns_reconstruction.constrained_candidate import CompactCandidate

PROBE_PATH = Path(__file__).parents[1] / "experiments" / "constrained_basis_growth_probe.py"
spec = importlib.util.spec_from_file_location("constrained_basis_growth_probe", PROBE_PATH)
probe = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = probe
spec.loader.exec_module(probe)


def finite_difference_divergence(field, points, time, step=2e-5):
    total = np.zeros(len(points))
    for j in range(3):
        delta = np.eye(3)[j] * step
        plus = field(points + delta, time)[:, j]
        minus = field(points - delta, time)[:, j]
        total += (plus - minus) / (2 * step)
    return total


def test_all_modes_are_numerically_divergence_free_and_compact():
    candidate = CompactCandidate(radial_width=1.05, axial_width=1.02)
    points = np.array([
        [0.23, -0.31, 0.17],
        [0.61, 0.42, -0.29],
        [0.91, -0.18, 0.37],
    ])
    for mode in range(len(probe.MODE_NAMES)):
        field = lambda x, t, mode=mode: probe.correction_mode(candidate, mode, x, t)
        divergence = finite_difference_divergence(field, points, 0.53)
        assert np.max(np.abs(divergence)) < 2e-7
        exterior = field(np.array([[5.0, 0.0, 0.0], [0.0, 0.0, 5.0]]), 0.53)
        assert np.array_equal(exterior, np.zeros((2, 3)))


def test_mode_validation_fails_closed():
    candidate = CompactCandidate()
    with pytest.raises(TypeError):
        probe.correction_mode(candidate, True, [[0.1, 0.0, 0.1]], 0.5)
    with pytest.raises(ValueError):
        probe.correction_mode(candidate, 8, [[0.1, 0.0, 0.1]], 0.5)
    with pytest.raises(ValueError):
        probe.correction_mode(candidate, 0, [[0.1, 0.0, 0.1]], 0.9)


def test_cylindrical_projection_separates_swirl():
    points = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    momentum = np.array([[
        [0.0, 2.0, 3.0],
        [-2.0, 0.0, 3.0],
    ]])
    result = probe._cylindrical_rms(momentum, points)
    assert result["radial"] == 0.0
    assert result["azimuthal"] == 2.0
    assert result["axial"] == 3.0
