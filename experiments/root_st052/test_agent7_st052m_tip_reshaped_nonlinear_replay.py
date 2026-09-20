from __future__ import annotations

import numpy as np

import agent7_st052m_tip_reshaped_nonlinear_replay as replay


def _radial_component(points, velocity):
    points = np.asarray(points, dtype=float)
    velocity = np.asarray(velocity, dtype=float)
    r = np.hypot(points[:, 0], points[:, 1])
    return (points[:, 0] * velocity[:, 0] + points[:, 1] * velocity[:, 1]) / r


def test_source_lock_uses_preregistered_p9_m4_coordinate():
    replay._assert_source_lock()
    p, _, z_star = replay.relocation.select_inner_exponent()
    assert p == 9
    assert replay.SELECTED_M == 4
    assert z_star >= 1.60
    assert replay.SOURCE_PARENT_PR == 738
    assert replay.SOURCE_PARENT_HEAD == "8072b099b8886e7c8716347ef2a60066d5985206"


def test_reshape_relocates_radial_sign_at_preregistered_z155_probe():
    pts = np.array(
        [
            [0.6, 0.0, 1.55],
            [0.9, 0.0, 1.55],
            [1.2, 0.0, 1.55],
            [0.6, 0.0, -1.55],
            [0.9, 0.0, -1.55],
            [1.2, 0.0, -1.55],
        ],
        dtype=float,
    )
    old_ur = _radial_component(pts, replay.screen.tip_odd_poloidal_correction(pts))
    new_ur = _radial_component(pts, replay.relocation.reshaped_tip_correction(pts, 9))
    assert np.all(old_ur > 0.0)
    assert np.all(new_ur < 0.0)


def test_new_scalar_coordinate_is_derived_from_new_unit_response():
    times = np.linspace(0.25, 0.75, 33)
    positions = np.zeros((33, 48, 3), dtype=float)
    positions[:, :, 0] = 1.0
    positions[:, :, 2] = 1.25
    metadata = [
        {"band": "tip" if i < 24 else "shoulder"}
        for i in range(48)
    ]

    def source_velocity(points, time):
        points = np.asarray(points, dtype=float)
        out = np.zeros_like(points)
        out[:, 0] = 0.2
        return out

    alpha, calibration = replay.derive_reshaped_alpha(
        source_velocity, times, positions, metadata
    )
    pts = positions[replay.screen.MID_INDEX, :24]
    base_ur = replay.screen.loc._radial_velocity(source_velocity, pts, 0.50)
    unit_ur = replay.screen.loc._radial_velocity(replay.reshaped_unit_time_correction, pts, 0.50)
    assert np.isclose(alpha, -np.mean(base_ur) / np.mean(unit_ur))
    assert abs(calibration["fixed_position_mean_after_linear_cancellation"]) < 1.0e-12
    assert calibration["fresh_714_data_used"] is False
    assert calibration["coefficient_grid_scan_performed"] is False


def test_targeted_outer_tip_cloud_is_frozen_and_disjoint_from_alpha_inputs_by_contract():
    pts = replay._target_probe_points()
    assert pts.shape == (24, 3)
    assert set(np.round(np.hypot(pts[:, 0], pts[:, 1]), 12)) == {0.6, 0.9, 1.2}
    assert set(np.round(np.abs(pts[:, 2]), 12)) == {1.55}
    assert replay.TARGET_PROBE_TIMES == (0.375, 0.50, 0.625, 0.75)
    assert replay.RELATED_FRESH_HOLDOUT_PR == 714


def test_truth_boundary_forbids_promotion_and_retuning():
    assert replay.PREREG_ISSUE == 739
    assert replay.TRUTH["basis_dimension_changed"] is False
    assert replay.TRUTH["same_dimension_envelope_replayed"] is True
    assert replay.TRUTH["second_poloidal_basis_added"] is False
    assert replay.TRUTH["new_temporal_basis_added"] is False
    assert replay.TRUTH["historical_development_paths_used_for_coefficient"] is True
    assert replay.TRUTH["fresh_714_path_data_used"] is False
    assert replay.TRUTH["parameter_grid_scan_performed"] is False
    assert replay.TRUTH["optimization_performed"] is False
    assert replay.TRUTH["pressure_or_force_changed"] is False
    assert replay.TRUTH["held_out_pde_residual_evaluated"] is False
    assert replay.TRUTH["visual_correspondence_verified"] is False
    assert replay.TRUTH["pde_validated"] is False
    assert replay.TRUTH["paper_exact"] is False
    assert replay.TRUTH["openai_field_identified"] is False
