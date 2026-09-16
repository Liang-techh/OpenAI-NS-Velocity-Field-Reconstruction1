import numpy as np
import pytest

from openai_ns_reconstruction.constrained_temporal_snapshot_consistency import (
    audit_temporal_snapshot_consistency,
)


class _Field:
    def __init__(self, mode: str = "quadratic"):
        self.mode = mode

    def at_points(self, points, time):
        points = np.asarray(points, dtype=float)
        t = float(time)
        x, y, z = points.T
        if self.mode == "affine":
            return np.column_stack(
                [(1.0 + x) * (1.0 + 0.4 * t), y - 0.7 * t, z + 0.2 * t]
            )
        return np.column_stack(
            [
                (1.0 + x) * t**2,
                (0.5 + y * y) * t**3,
                (1.0 + z) * np.sin(2.0 * t),
            ]
        )


def _points():
    rng = np.random.default_rng(914083)
    return rng.uniform(-0.8, 0.8, size=(64, 3))


def _times():
    return np.linspace(0.271, 0.729, 23)


def test_affine_in_time_is_reproduced_to_roundoff():
    report = audit_temporal_snapshot_consistency(
        _Field("affine"), _points(), _times(), snapshot_counts=(3, 5, 9)
    )
    assert max(level.vector_error_max for level in report.levels) < 2.0e-15
    assert report.claim_scope == "visualization_time_sampling_only"
    assert report.pde_validated is False
    assert report.visual_correspondence_verified is False


def test_nonlinear_time_error_contracts_under_snapshot_refinement():
    report = audit_temporal_snapshot_consistency(
        _Field(), _points(), _times(), snapshot_counts=(3, 5, 9)
    )
    errors = [level.vector_error_rms for level in report.levels]
    assert errors[0] > errors[1] > errors[2] > 0.0
    assert errors[1] / errors[0] < 0.35
    assert errors[2] / errors[1] < 0.35
    assert [level.max_time_spacing for level in report.levels] == pytest.approx(
        [0.25, 0.125, 0.0625]
    )


def test_fail_closed_inputs_and_zero_reference():
    with pytest.raises(ValueError):
        audit_temporal_snapshot_consistency(
            _Field(), _points(), _times(), snapshot_counts=(3, 5)
        )
    with pytest.raises(ValueError):
        audit_temporal_snapshot_consistency(
            _Field(), _points(), [0.3, 0.3, 0.4], snapshot_counts=(3, 5, 9)
        )

    class Zero:
        def at_points(self, points, time):
            return np.zeros((len(points), 3))

    with pytest.raises(ValueError):
        audit_temporal_snapshot_consistency(
            Zero(), _points(), _times(), snapshot_counts=(3, 5, 9)
        )

    class Bad:
        def at_points(self, points, time):
            return np.zeros((len(points), 2))

    with pytest.raises(ValueError):
        audit_temporal_snapshot_consistency(
            Bad(), _points(), _times(), snapshot_counts=(3, 5, 9)
        )


def test_packaged_velocity_public_api_smoke():
    from openai_ns_reconstruction.velocity_components import VelocityField

    rng = np.random.default_rng(914089)
    points = rng.uniform(-0.7, 0.7, size=(24, 3))
    probe_times = np.linspace(0.271, 0.729, 17)
    report = audit_temporal_snapshot_consistency(
        VelocityField(),
        points,
        probe_times,
        snapshot_counts=(3, 5, 9),
        time_interval=(0.25, 0.75),
    )
    assert report.point_count == 24
    assert report.probe_time_count == 17
    assert all(np.isfinite(level.relative_rms) for level in report.levels)
    assert report.pde_validated is False
    assert report.visual_correspondence_verified is False
