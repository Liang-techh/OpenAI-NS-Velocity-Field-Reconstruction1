import numpy as np
import pytest

from openai_ns_reconstruction.constrained_pressure_compatibility import (
    pressure_compatibility_report,
    pressure_poisson_terms,
)


def _zero_force(points, times):
    return np.zeros_like(points, dtype=float)


def test_solid_rotation_has_compatible_pressure():
    omega = 1.7

    def velocity(points, times):
        x, y, _ = points.T
        return np.column_stack((-omega * y, omega * x, np.zeros_like(x)))

    def pressure(points, times):
        x, y, _ = points.T
        return 0.5 * omega**2 * (x * x + y * y)

    points = np.array([[0.2, -0.3, 0.1], [-0.5, 0.4, -0.2], [0.7, 0.1, 0.3]])
    report = pressure_compatibility_report(
        velocity, pressure, _zero_force, points, np.array([0.3, 0.5, 0.7]), step=1e-3
    )
    assert report.poisson_defect_max_abs < 2e-8
    assert report.velocity_divergence_max_abs < 1e-12


def test_wrong_pressure_sign_is_detected():
    omega = 1.25

    def velocity(points, times):
        x, y, _ = points.T
        return np.column_stack((-omega * y, omega * x, np.zeros_like(x)))

    def wrong_pressure(points, times):
        x, y, _ = points.T
        return -0.5 * omega**2 * (x * x + y * y)

    points = np.array([[0.2, 0.3, 0.1], [0.6, -0.4, 0.2]])
    terms = pressure_poisson_terms(velocity, wrong_pressure, _zero_force, points, 0.5)
    assert np.allclose(terms["defect"], -4.0 * omega**2, atol=2e-8, rtol=0.0)


def test_force_divergence_enters_with_fixed_sign():
    def zero_velocity(points, times):
        return np.zeros_like(points, dtype=float)

    def zero_pressure(points, times):
        return np.zeros(len(points), dtype=float)

    def force(points, times):
        out = np.zeros_like(points, dtype=float)
        out[:, 0] = points[:, 0]
        return out

    points = np.array([[0.1, 0.2, 0.3], [-0.4, 0.5, -0.1]])
    terms = pressure_poisson_terms(zero_velocity, zero_pressure, force, points, 0.4)
    assert np.allclose(terms["force_divergence"], 1.0, atol=1e-12, rtol=0.0)
    assert np.allclose(terms["defect"], -1.0, atol=1e-12, rtol=0.0)


def test_exterior_gauge_is_reported_separately():
    def zero_velocity(points, times):
        return np.zeros_like(points, dtype=float)

    def compact_pressure(points, times):
        inside = np.max(np.abs(points), axis=1) < 1.0
        return np.where(inside, 0.25, 0.0)

    points = np.array([[0.1, 0.2, 0.3], [-0.2, 0.1, 0.2]])
    exterior = np.array([[1.5, 0.0, 0.0], [0.0, -1.2, 0.0]])
    report = pressure_compatibility_report(
        zero_velocity,
        compact_pressure,
        _zero_force,
        points,
        0.5,
        exterior_points=exterior,
    )
    assert report.exterior_pressure_max_abs == 0.0


def test_invalid_samples_and_field_shapes_fail_closed():
    def zero_velocity(points, times):
        return np.zeros_like(points, dtype=float)

    def zero_pressure(points, times):
        return np.zeros(len(points), dtype=float)

    points = np.array([[0.1, 0.2, 0.3]])
    with pytest.raises(ValueError):
        pressure_poisson_terms(zero_velocity, zero_pressure, _zero_force, points, 0.5, step=0.0)
    with pytest.raises(ValueError):
        pressure_poisson_terms(zero_velocity, zero_pressure, _zero_force, np.zeros((2, 2)), 0.5)
    with pytest.raises(ValueError):
        pressure_poisson_terms(zero_velocity, lambda p, t: np.zeros((len(p), 1)), _zero_force, points, 0.5)
