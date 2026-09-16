import numpy as np
import pytest

from openai_ns_reconstruction.constrained_boundary_support import boundary_support_report


def compact_velocity(points, time):
    points = np.asarray(points, dtype=float)
    bump = np.maximum(0.0, 1.0 - np.sum(points * points, axis=-1)) ** 3
    swirl = np.stack(
        (-points[..., 1], points[..., 0], 0 * points[..., 2]), axis=-1
    )
    return bump[..., None] * swirl * (1 + time)


def leaky_velocity(points, time):
    points = np.asarray(points, dtype=float)
    return compact_velocity(points, time) + np.array([1e-5, 0.0, 0.0])


def test_compact_field_has_zero_box_boundary_and_exterior_shell():
    report = boundary_support_report(
        compact_velocity,
        time=0.5,
        support_half_width=1.0,
        probe_half_width=1.25,
        points_per_axis=15,
        tolerance=0.0,
    )
    assert report.face_max_speed == 0.0
    assert report.exterior_shell_max_speed == 0.0
    assert report.sampled_support_pass


def test_leak_mutation_is_detected():
    report = boundary_support_report(
        leaky_velocity,
        time=0.5,
        support_half_width=1.0,
        probe_half_width=1.25,
        points_per_axis=11,
        tolerance=1e-8,
    )
    assert report.face_max_speed >= 1e-5
    assert report.exterior_shell_max_speed >= 1e-5
    assert not report.sampled_support_pass


def test_fail_closed_on_invalid_geometry_and_velocity_shape():
    with pytest.raises(ValueError):
        boundary_support_report(
            compact_velocity,
            time=0.5,
            support_half_width=1,
            probe_half_width=1,
        )
    with pytest.raises(ValueError):
        boundary_support_report(
            compact_velocity,
            time=0.5,
            support_half_width=1,
            probe_half_width=2,
            points_per_axis=True,
        )
    with pytest.raises(ValueError):
        boundary_support_report(
            compact_velocity,
            time=0.5,
            support_half_width=1,
            probe_half_width=2,
            tolerance=-1,
        )
    with pytest.raises(ValueError):
        boundary_support_report(
            lambda p, t: np.zeros((len(p), 2)),
            time=0.5,
            support_half_width=1,
            probe_half_width=2,
        )
