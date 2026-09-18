import numpy as np
import pytest

import agent7_st048s_piola_temporal_screen as temporal
import agent7_st048s_piola_warp_screen as fixed


class _SyntheticParent:
    def at_points(self, points, time):
        points = np.asarray(points, dtype=float)
        x, y, z = points.T
        # Smooth nonzero field used only for representation-contract checks.
        return np.column_stack((1.0 + 0.1 * y, -0.2 + 0.05 * x, 0.3 + 0.02 * z + 0.01 * float(time)))


def test_beta_schedule_and_bounds():
    assert temporal.beta_at_time(0.50, 0.125) == pytest.approx(0.075)
    assert temporal.beta_at_time(0.25, 0.125) == pytest.approx(0.20)
    assert temporal.beta_at_time(0.75, 0.125) == pytest.approx(0.20)
    assert temporal.beta_at_time(0.375, 0.05) == pytest.approx(0.0875)
    with pytest.raises(ValueError):
        temporal.beta_at_time(0.50, 0.126)
    with pytest.raises(ValueError):
        temporal.beta_at_time(0.20, 0.05)


def test_gamma_zero_replays_fixed_beta_piola_field():
    parent = _SyntheticParent()
    temporal_field = temporal.TemporalPiolaWarpField(parent, 0.0)
    fixed_field = fixed._PiolaWarpField(parent, temporal.BASE_BETA)
    points = np.array([
        [0.1, 0.2, 0.3],
        [0.7, -0.4, -0.8],
        [1.1, 0.1, 1.4],
        [0.2, -0.3, -1.7],
    ])
    for time in temporal.DEFAULT_TIMES:
        np.testing.assert_allclose(
            temporal_field.at_points(points, time),
            fixed_field.at_points(points, time),
            rtol=0.0,
            atol=0.0,
        )


def test_temporal_schedule_preserves_piola_spatial_divergence_for_simple_field():
    class HorizontalConstant:
        def at_points(self, points, time):
            points = np.asarray(points, dtype=float)
            return np.column_stack((np.ones(len(points)), np.full(len(points), 0.25), np.zeros(len(points))))

    field = temporal.TemporalPiolaWarpField(HorizontalConstant(), 0.125)
    points = np.array([[0.31, 0.22, 0.17], [0.62, -0.27, -0.41], [0.88, 0.31, 0.67]])
    h = 1.0e-6
    for time in temporal.DEFAULT_TIMES:
        divergence = np.zeros(len(points))
        for axis in range(3):
            offset = np.zeros_like(points)
            offset[:, axis] = h
            plus = field.at_points(points + offset, time)
            minus = field.at_points(points - offset, time)
            divergence += (plus[:, axis] - minus[:, axis]) / (2.0 * h)
        assert np.max(np.abs(divergence)) < 1.0e-9
