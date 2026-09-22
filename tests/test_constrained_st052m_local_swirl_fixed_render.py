import numpy as np

from openai_ns_reconstruction import constrained_st052m_local_swirl_fixed_render as m


def _solid_rotation(points, time):
    del time
    p = np.asarray(points, dtype=float)
    return np.column_stack((-p[:, 1], p[:, 0], np.zeros(len(p))))


def _helical_inflow(points, time):
    del time
    p = np.asarray(points, dtype=float)
    return np.column_stack(
        (
            -p[:, 1] - 0.10 * p[:, 0],
            p[:, 0] - 0.10 * p[:, 1],
            np.full(len(p), 0.20),
        )
    )


def test_frozen_seed_source_and_truth_contract():
    seeds = m.seed_points()
    assert seeds.shape == (48, 3)
    assert len(np.unique(np.round(seeds, 12), axis=0)) == 48
    assert m.REFERENCE_TIMES == (0.25, 0.50, 0.75)
    assert m.GRID_RESOLUTION == 33
    assert m.LOCAL_SWIRL_HEAD == "39b106ad8cb8df2064cabead3a12682089575e74"
    assert m.MORPHOLOGY_CONVERGENCE_HEAD == "0ddf5f6b321ba78e27aa3ca3e073ef593e55869a"
    assert abs(m.EXPECTED_BETA - 0.08837490297155456) < 1e-16
    assert not any(m.TRUTH.values())


def test_cartesian_vorticity_on_solid_rotation():
    axis, fields = m.sample_velocity_grid(_solid_rotation, resolution=17)
    spacing = float(axis[1] - axis[0])
    omega, mag = m.vorticity(fields[0], spacing)
    np.testing.assert_allclose(omega[..., 0], 0.0, atol=2e-14, rtol=0)
    np.testing.assert_allclose(omega[..., 1], 0.0, atol=2e-14, rtol=0)
    np.testing.assert_allclose(omega[..., 2], 2.0, atol=2e-14, rtol=0)
    np.testing.assert_allclose(mag, 2.0, atol=2e-14, rtol=0)


def test_identical_pair_has_zero_relative_diagnostics():
    axis, fields = m.sample_velocity_grid(_helical_inflow, resolution=17)
    report = m.analyze_pair(axis, fields, fields.copy())
    assert report["contract"]["streamline_count"] == 48
    assert report["contract"]["image_or_openai_numeric_target_used"] is False
    assert report["contract"]["visual_acceptance_threshold_used"] is False
    assert not any(report["truth"].values())
    for row in report["per_time"]:
        assert row["control"]["streamlines"]["line_count"] == 48
        assert row["child"]["vorticity"]["selected_points"] > 0
        for value in row["relative"].values():
            assert abs(value) < 5e-13


def test_threshold_free_vorticity_metrics_are_finite():
    axis, fields = m.sample_velocity_grid(_solid_rotation, resolution=17)
    spacing = float(axis[1] - axis[0])
    _, mag = m.vorticity(fields[0], spacing)
    metrics = m.vorticity_metrics(mag, axis)
    assert metrics["enstrophy_axial_rms"] > 0.0
    assert metrics["enstrophy_radial_rms"] > 0.0
    assert metrics["enstrophy_aspect"] > 0.0
    assert np.isfinite(list(metrics.values())).all()


def test_grid_interpolator_returns_finite_values():
    axis, fields = m.sample_velocity_grid(_solid_rotation, resolution=17)
    evaluate = m.grid_interpolator(axis, fields[1])
    points = np.array([[0.25, -0.5, 0.1], [-0.7, 0.3, -0.2]])
    out = evaluate(points)
    assert out.shape == (2, 3)
    assert np.isfinite(out).all()
    np.testing.assert_allclose(out[:, :2], _solid_rotation(points, 0.5)[:, :2], atol=1e-12)
