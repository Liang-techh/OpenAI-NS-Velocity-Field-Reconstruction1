import json
from pathlib import Path

import numpy as np

from openai_ns_reconstruction.constrained_channel_capacity import (
    PEAK_RADIUS,
    PEAK_Z,
    perturbed_velocity,
    poloidal_curl_probe,
    pure_swirl_probe,
)
from openai_ns_reconstruction.constrained_force import RestrictedForce
from openai_ns_reconstruction.constrained_inner_swirl import InnerSwirlCandidate
from openai_ns_reconstruction.constrained_validation import residual


CANDIDATE = Path("artifacts/constrained/inner_swirl_pressure/candidate.json")
VALIDATION = Path("artifacts/constrained/inner_swirl_pressure/validation.json")


def _candidate_force():
    candidate = InnerSwirlCandidate.load(CANDIDATE)
    force_data = json.loads(VALIDATION.read_text())["force"]
    return candidate, RestrictedForce(**force_data)


def test_pure_swirl_is_numerically_null_in_axial_momentum_channel():
    candidate, force = _candidate_force()
    points = np.array(
        [
            [PEAK_RADIUS, 0.0, PEAK_Z],
            [0.95, 0.0, 0.62],
            [1.08, 0.0, 0.74],
            [1.15, 0.0, 0.55],
        ]
    )
    kwargs = dict(points=points, time=0.748, nu=0.01, step=0.005, time_bounds=(0.25, 0.75))
    base = residual(candidate.velocity, candidate.pressure, force, **kwargs)["momentum"]
    changed = residual(
        perturbed_velocity(candidate.velocity, pure_swirl_probe, 0.35),
        candidate.pressure,
        force,
        **kwargs,
    )["momentum"]
    assert np.max(np.abs(changed[:, 2] - base[:, 2])) < 1.0e-10
    assert np.max(np.abs(changed[:, :2] - base[:, :2])) > 1.0e-4


def test_one_poloidal_curl_direction_has_nonzero_axial_sensitivity():
    candidate, force = _candidate_force()
    points = np.array(
        [
            [PEAK_RADIUS, 0.0, PEAK_Z],
            [0.98, 0.0, 0.66],
            [1.05, 0.0, 0.72],
        ]
    )
    kwargs = dict(points=points, time=0.748, nu=0.01, step=0.005, time_bounds=(0.25, 0.75))
    base = residual(candidate.velocity, candidate.pressure, force, **kwargs)["momentum"]
    plus = residual(
        perturbed_velocity(candidate.velocity, poloidal_curl_probe, 1.0e-3),
        candidate.pressure,
        force,
        **kwargs,
    )["momentum"]
    minus = residual(
        perturbed_velocity(candidate.velocity, poloidal_curl_probe, -1.0e-3),
        candidate.pressure,
        force,
        **kwargs,
    )["momentum"]
    axial_jacobian = (plus[:, 2] - minus[:, 2]) / (2.0e-3)
    assert np.max(np.abs(axial_jacobian)) > 1.0e-2
    assert np.max(np.abs(base[:, 2])) > 0.1


def test_poloidal_probe_is_divergence_free_and_compact_on_interior_sample():
    rng = np.random.default_rng(7703)
    radius = rng.uniform(0.7, 1.3, 24)
    theta = rng.uniform(0.0, 2.0 * np.pi, 24)
    z = rng.uniform(0.3, 1.1, 24)
    points = np.column_stack((radius * np.cos(theta), radius * np.sin(theta), z))
    step = 2.0e-6
    divergence = np.zeros(len(points))
    for axis in range(3):
        direction = np.eye(3)[axis] * step
        plus = poloidal_curl_probe(points + direction, 0.748)
        minus = poloidal_curl_probe(points - direction, 0.748)
        divergence += (plus[:, axis] - minus[:, axis]) / (2.0 * step)
    assert np.max(np.abs(divergence)) < 2.0e-6

    exterior = np.array([[2.0, 0.0, 0.0], [0.9, 0.0, 2.0], [2.2, 0.0, 0.7]])
    assert np.max(np.abs(poloidal_curl_probe(exterior, 0.748))) == 0.0
    assert np.max(np.abs(pure_swirl_probe(exterior, 0.748))) == 0.0


def test_invalid_probe_inputs_fail_closed():
    with np.testing.assert_raises(ValueError):
        pure_swirl_probe(np.zeros((2, 2)), 0.5)
    with np.testing.assert_raises(ValueError):
        poloidal_curl_probe(np.array([[np.nan, 0.0, 0.0]]), 0.5)
    with np.testing.assert_raises(ValueError):
        pure_swirl_probe(np.zeros((2, 3)), np.zeros(3))
    with np.testing.assert_raises(ValueError):
        perturbed_velocity(lambda p, t: np.zeros_like(p), pure_swirl_probe, np.inf)
