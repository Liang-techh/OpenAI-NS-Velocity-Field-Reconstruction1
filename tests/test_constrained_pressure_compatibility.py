import numpy as np
import pytest

from openai_ns_reconstruction.constrained_eq45_candidate import Eq45VelocityCandidate
from openai_ns_reconstruction.constrained_pressure_compatibility import (
    audit_pressure_compatibility,
)


def _rigid_rotation(points, times):
    points = np.asarray(points, dtype=float)
    return np.stack(
        (-points[:, 1], points[:, 0], np.zeros(points.shape[0])), axis=-1
    )


def _viscous_incompatible(points, times):
    points = np.asarray(points, dtype=float)
    return np.stack(
        (
            np.sin(points[:, 1]),
            np.zeros(points.shape[0]),
            np.zeros(points.shape[0]),
        ),
        axis=-1,
    )


def test_pressure_compatible_rigid_rotation_has_zero_curl_obstruction():
    points = np.array(
        [
            [0.10, -0.40, 0.10],
            [0.20, -0.10, -0.20],
            [-0.30, 0.20, 0.30],
            [0.40, 0.35, -0.10],
        ]
    )
    report = audit_pressure_compatibility(_rigid_rotation, points, 0.5)

    assert [level.spatial_step for level in report.levels] == [0.02, 0.01, 0.005]
    assert all(level.momentum_rms > 0.0 for level in report.levels)
    assert all(level.curl_momentum_max < 1e-12 for level in report.levels)
    assert all(level.normalized_curl_rms < 1e-12 for level in report.levels)


def test_incompatible_viscous_field_recovers_nonzero_pressure_obstruction():
    points = np.array(
        [
            [0.10, -0.40, 0.10],
            [0.20, -0.10, -0.20],
            [-0.30, 0.20, 0.30],
            [0.40, 0.35, -0.10],
        ]
    )
    nu = 0.01
    report = audit_pressure_compatibility(_viscous_incompatible, points, 0.5, nu=nu)
    expected_rms = nu * np.sqrt(np.mean(np.cos(points[:, 1]) ** 2))
    errors = [abs(level.curl_momentum_rms - expected_rms) for level in report.levels]

    assert errors[2] < errors[1] < errors[0]
    assert report.levels[-1].curl_momentum_rms == pytest.approx(expected_rms, rel=3e-5)
    assert report.levels[-1].curl_momentum_max > 9e-3
    assert report.levels[-1].normalized_curl_rms > 1.0


def test_report_is_fail_closed_about_scope_and_stencil_configuration():
    points = np.array([[0.10, 0.15, -0.10], [-0.20, 0.25, 0.20]])
    report = audit_pressure_compatibility(
        _rigid_rotation,
        points,
        np.array([0.4, 0.6]),
        nu=0.01,
        derivative_steps=(0.02, 0.01),
    )
    payload = report.to_dict()

    assert payload["forcing_mode"] == "zero_forcing_subcase"
    assert payload["pressure_model"] == "scalar_pressure_gradient_only"
    assert payload["velocity_changed"] is False
    assert payload["forcing_fitted"] is False
    assert payload["pressure_fitted"] is False
    assert payload["pde_validated"] is False
    assert "necessary" in payload["interpretation"]

    with pytest.raises(ValueError, match="strictly decreasing"):
        audit_pressure_compatibility(
            _rigid_rotation, points, 0.5, derivative_steps=(0.01, 0.02)
        )
    with pytest.raises(ValueError, match="at least two"):
        audit_pressure_compatibility(
            _rigid_rotation, points, 0.5, derivative_steps=(0.01,)
        )


def test_eq45_candidate_public_velocity_path_pins_obstruction_calibration():
    candidate = Eq45VelocityCandidate.seed()
    points = np.array(
        [
            [0.15, 0.20, -0.30],
            [-0.20, 0.25, 0.35],
        ]
    )
    times = np.array([0.40, 0.60])
    report = audit_pressure_compatibility(candidate, points, times, nu=0.01)

    assert report.point_count == 2
    assert report.nu == pytest.approx(0.01)
    for level in report.levels:
        values = np.array(
            [
                level.momentum_rms,
                level.momentum_max,
                level.curl_momentum_rms,
                level.curl_momentum_max,
                level.normalized_curl_rms,
            ]
        )
        assert np.all(np.isfinite(values))
        assert np.all(values >= 0.0)
        assert level.momentum_rms > 0.0

    expected = (
        (1.3849572125, 3.8671981168, 2.7922870700),
        (1.3860084226, 3.8181367081, 2.7547716491),
        (1.3862710976, 3.8060632118, 2.7455403336),
    )
    for level, (momentum_rms, curl_rms, normalized) in zip(report.levels, expected):
        assert level.momentum_rms == pytest.approx(momentum_rms, rel=2e-6)
        assert level.curl_momentum_rms == pytest.approx(curl_rms, rel=2e-6)
        assert level.normalized_curl_rms == pytest.approx(normalized, rel=2e-6)

    assert report.levels[-1].curl_momentum_max == pytest.approx(5.1334835987, rel=2e-6)
