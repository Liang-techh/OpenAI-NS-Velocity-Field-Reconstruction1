import numpy as np
import pytest

from openai_ns_reconstruction.constrained_streamline_resample import (
    resample_streamlines_for_rendering,
)


def _irregular_helix(phase=0.0):
    u = np.linspace(0.0, 1.0, 91) ** 1.7
    theta = phase + 5.0 * np.pi * u
    z = -1.2 + 2.4 * u
    return np.column_stack((0.7 * np.cos(theta), 0.7 * np.sin(theta), z))


def test_pchip_arclength_resample_preserves_dense_helix_geometry():
    line = _irregular_helix()
    result = resample_streamlines_for_rendering(
        [line], points_per_line=301, color_mode="height"
    )

    out = result.points[0]
    radius = np.linalg.norm(out[:, :2], axis=1)
    np.testing.assert_allclose(out[0], line[0], rtol=0.0, atol=0.0)
    np.testing.assert_allclose(out[-1], line[-1], rtol=0.0, atol=0.0)
    assert np.max(np.abs(radius - 0.7)) < 6.0e-3
    assert np.all(np.diff(out[:, 2]) > 0.0)
    assert result.segment_length_cv[0] < 0.03
    np.testing.assert_allclose(result.scalars[0], out[:, 2])
    assert result.points.shape == (1, 301, 3)
    assert not result.points.flags.writeable
    assert result.claim_scope == "visualization_polyline_resampling_only"
    assert result.visualization_ready is False
    assert result.pde_validated is False


def test_speed_coloring_uses_only_resampled_points_and_declared_time():
    lines = [_irregular_helix(0.0), _irregular_helix(0.4)]

    def velocity(points, time):
        p = np.asarray(points)
        return np.column_stack((-p[:, 1], p[:, 0], np.full(len(p), 0.5 + time)))

    result = resample_streamlines_for_rendering(
        lines,
        points_per_line=129,
        color_mode="speed",
        velocity=velocity,
        time=0.5,
    )
    flat = result.points.reshape(-1, 3)
    expected = np.linalg.norm(velocity(flat, 0.5), axis=1).reshape(2, 129)
    np.testing.assert_allclose(result.scalars, expected, rtol=1e-13, atol=1e-13)
    assert result.time == 0.5
    assert result.points.shape == (2, 129, 3)


def test_consecutive_duplicate_points_are_collapsed_without_endpoint_drift():
    line = _irregular_helix()
    duplicated = np.insert(line, [10, 40], [line[9], line[39]], axis=0)
    result = resample_streamlines_for_rendering([duplicated], points_per_line=80)
    np.testing.assert_allclose(result.points[0, 0], line[0], rtol=0.0, atol=0.0)
    np.testing.assert_allclose(result.points[0, -1], line[-1], rtol=0.0, atol=0.0)


def test_fail_closed_on_bad_lines_and_speed_contract():
    with pytest.raises(ValueError):
        resample_streamlines_for_rendering([])
    with pytest.raises(ValueError):
        resample_streamlines_for_rendering([np.zeros((3, 2))])
    with pytest.raises(ValueError):
        resample_streamlines_for_rendering([np.zeros((3, 3))])
    bad = _irregular_helix()
    bad[4, 1] = np.nan
    with pytest.raises(ValueError):
        resample_streamlines_for_rendering([bad])
    with pytest.raises(ValueError):
        resample_streamlines_for_rendering([_irregular_helix()], color_mode="speed")

    def wrong_velocity(points, time):
        return np.zeros((len(points), 2))

    with pytest.raises(ValueError):
        resample_streamlines_for_rendering(
            [_irregular_helix()], color_mode="speed", velocity=wrong_velocity, time=0.5
        )
