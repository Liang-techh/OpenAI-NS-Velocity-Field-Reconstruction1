import numpy as np
import pytest

from openai_ns_reconstruction.constrained_pressure_gauge import (
    canonicalize_exterior_pressure_gauge,
)


def _compact_pressure_with_time_gauge(points, time):
    points = np.asarray(points, dtype=float)
    time = np.asarray(time, dtype=float)
    x, y, z = np.moveaxis(points, -1, 0)
    r2 = x*x + y*y
    inside = (r2 < 4.0) & (np.abs(z) < 2.0)
    value = 2.75 + 0.4*time
    core = np.zeros_like(value)
    core[inside] = ((1.0-r2[inside]/4.0)**4
                    *(1.0-z[inside]**2/4.0)**4
                    *(x[inside]+0.25*z[inside]+0.1*time[inside]))
    return value + core


def _probes():
    return np.array([
        [2.0, 0.0, 0.1],
        [0.0, -2.1, -0.2],
        [0.2, 0.3, 2.0],
        [-0.3, 0.4, -2.2],
    ])


def test_time_dependent_exterior_gauge_is_zero_and_preserves_pressure_gradient():
    gauged, receipt = canonicalize_exterior_pressure_gauge(
        _compact_pressure_with_time_gauge, _probes(), [0.25, 0.5, 0.75]
    )
    for time in (0.25, 0.5, 0.75):
        np.testing.assert_allclose(gauged(_probes(), time), 0.0, atol=1e-14)

    points = np.array([[0.2, 0.3, 0.1], [-0.4, 0.2, -0.3]])
    time = np.array([0.4, 0.6])
    h = 1e-6
    for axis in range(3):
        delta = np.zeros_like(points)
        delta[:, axis] = h
        base_gradient = (
            _compact_pressure_with_time_gauge(points+delta, time)
            - _compact_pressure_with_time_gauge(points-delta, time)
        )/(2*h)
        gauged_gradient = (
            gauged(points+delta, time)-gauged(points-delta, time)
        )/(2*h)
        np.testing.assert_allclose(gauged_gradient, base_gradient, rtol=0, atol=1e-9)

    assert receipt.max_abs_gauged_exterior <= 1e-10
    assert receipt.reference_offsets == (2.85, 2.95, 3.05)
    assert receipt.gradient_invariant_by_construction
    assert not receipt.velocity_changed
    assert not receipt.forcing_changed
    assert not receipt.pde_validated


def test_nonconstant_exterior_tail_fails_closed():
    def bad_pressure(points, time):
        points = np.asarray(points, dtype=float)
        time = np.asarray(time, dtype=float)
        return 1.0 + 0.2*time + 1e-5*points[:, 0]

    with pytest.raises(ValueError, match="spatially nonconstant"):
        canonicalize_exterior_pressure_gauge(
            bad_pressure, _probes(), [0.25, 0.5, 0.75], tolerance=1e-10
        )


def test_probe_inside_support_is_rejected():
    probes = _probes()
    probes[1] = [0.2, 0.2, 0.2]
    with pytest.raises(ValueError, match="outside compact support"):
        canonicalize_exterior_pressure_gauge(
            _compact_pressure_with_time_gauge, probes, [0.5]
        )


def test_nonfinite_or_wrong_shape_pressure_output_is_rejected():
    with pytest.raises(ValueError, match="shape"):
        canonicalize_exterior_pressure_gauge(
            lambda points, time: np.zeros((len(points), 1)), _probes(), [0.5]
        )
    with pytest.raises(ValueError, match="nonfinite"):
        canonicalize_exterior_pressure_gauge(
            lambda points, time: np.full(len(points), np.nan), _probes(), [0.5]
        )
