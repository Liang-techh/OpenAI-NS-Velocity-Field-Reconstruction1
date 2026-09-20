from __future__ import annotations

import numpy as np

import agent7_st052m_outer_reservoir_nonlinear_replay as replay


def _radial_component(points, velocity):
    points = np.asarray(points, dtype=float)
    velocity = np.asarray(velocity, dtype=float)
    r = np.hypot(points[:, 0], points[:, 1])
    return (points[:, 0] * velocity[:, 0] + points[:, 1] * velocity[:, 1]) / r


def test_source_lock_binds_exact_767_reservoir_parent():
    replay._assert_source_lock()
    assert replay.PREREG_ISSUE == 774
    assert replay.SOURCE_PARENT_PR == 767
    assert replay.SOURCE_PARENT_HEAD == "0930a1b7eeb357dee3e83dd79a3c4fea293e9314"
    assert replay.reservoir.Z_VISIBLE_OUTER == 1.8
    assert replay.reservoir.Z_GLOBAL_OUTER == 2.0
    assert replay.reservoir.R_SUPPORT_MAX == 1.6


def test_reservoir_channel_moves_return_flow_beyond_visible_tip_lobe():
    pts = np.asarray(
        [
            [1.0, 0.0, 1.75],
            [1.0, 0.0, -1.75],
            [1.0, 0.0, 1.85],
            [1.0, 0.0, -1.85],
        ],
        dtype=float,
    )
    new_ur = _radial_component(pts, replay.reservoir.reservoir_correction(pts))
    assert np.all(new_ur[:2] < 0.0)
    assert np.all(new_ur[2:] > 0.0)

    old_pts = pts[:2]
    p9_ur = _radial_component(
        old_pts, replay.p9.relocation.reshaped_tip_correction(old_pts, replay.p9.SELECTED_P)
    )
    assert np.all(p9_ur > 0.0)


def test_reservoir_scalar_coordinate_is_derived_only_from_historical_unit_response():
    times = np.linspace(0.25, 0.75, 33)
    positions = np.zeros((33, 48, 3), dtype=float)
    positions[:, :, 0] = 1.0
    positions[:, :, 2] = 1.25
    metadata = [{"band": "tip" if i < 24 else "shoulder"} for i in range(48)]

    def source_velocity(points, time):
        points = np.asarray(points, dtype=float)
        out = np.zeros_like(points)
        out[:, 0] = 0.2
        return out

    alpha, calibration = replay.derive_reservoir_alpha(
        source_velocity, times, positions, metadata
    )
    pts = positions[replay.screen.MID_INDEX, :24]
    base_ur = replay.screen.loc._radial_velocity(source_velocity, pts, 0.50)
    unit_ur = replay.screen.loc._radial_velocity(
        replay.reservoir_unit_time_correction, pts, 0.50
    )
    assert np.isclose(alpha, -np.mean(base_ur) / np.mean(unit_ur))
    assert abs(calibration["fixed_position_mean_after_linear_cancellation"]) < 1.0e-12
    assert calibration["fresh_714_data_used"] is False
    assert calibration["coefficient_grid_scan_performed"] is False
    assert calibration["optimizer_used"] is False


def test_visible_and_outer_probe_bands_are_frozen_before_replay():
    assert replay.VISIBLE_ABS_Z == (1.55, 1.75)
    assert replay.RESERVOIR_ABS_Z == (1.85, 1.95)
    assert replay.PROBE_TIMES == (0.375, 0.50, 0.625, 0.75)
    for abs_z in (*replay.VISIBLE_ABS_Z, *replay.RESERVOIR_ABS_Z):
        pts = replay._probe_points(abs_z)
        assert pts.shape == (24, 3)
        assert set(np.round(np.hypot(pts[:, 0], pts[:, 1]), 12)) == {0.6, 0.9, 1.2}
        assert set(np.round(np.abs(pts[:, 2]), 12)) == {abs_z}


def test_reservoir_spatial_channel_is_zero_at_shoulder_and_global_support_edge():
    pts = np.asarray(
        [
            [0.9, 0.0, 0.85],
            [0.9, 0.0, -0.85],
            [0.9, 0.0, 2.0],
            [0.9, 0.0, -2.0],
            [1.61, 0.0, 1.55],
            [-1.61, 0.0, -1.55],
        ],
        dtype=float,
    )
    assert np.max(np.abs(replay.reservoir.reservoir_correction(pts))) == 0.0


def test_truth_boundary_forbids_promotion_retuning_and_fresh_data_leakage():
    assert replay.TRUTH["basis_dimension_changed"] is False
    assert replay.TRUTH["single_tip_poloidal_channel_replayed"] is True
    assert replay.TRUTH["second_poloidal_basis_added"] is False
    assert replay.TRUTH["new_temporal_basis_added"] is False
    assert replay.TRUTH["historical_development_paths_used_for_coefficient"] is True
    assert replay.TRUTH["fresh_714_path_data_used"] is False
    assert replay.TRUTH["parameter_grid_scan_performed"] is False
    assert replay.TRUTH["optimization_performed"] is False
    assert replay.TRUTH["post_result_damping_or_retuning_performed"] is False
    assert replay.TRUTH["pressure_or_force_changed"] is False
    assert replay.TRUTH["held_out_pde_residual_evaluated"] is False
    assert replay.TRUTH["visual_correspondence_verified"] is False
    assert replay.TRUTH["pde_validated"] is False
    assert replay.TRUTH["paper_exact"] is False
    assert replay.TRUTH["openai_field_identified"] is False
