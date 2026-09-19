import numpy as np

from openai_ns_reconstruction import constrained_st052m_offmidplane_material_path as m


def test_offmidplane_seed_contract_hits_frozen_active_windows():
    seeds, metadata = m._offmidplane_seed_table()
    assert seeds.shape == (48, 3)
    assert len(metadata) == 48
    assert sum(row["band"] == "shoulder" for row in metadata) == 24
    assert sum(row["band"] == "tip" for row in metadata) == 24
    assert set(row["radius"] for row in metadata) == set(m.SEED_RADII)
    assert set(row["angle_index"] for row in metadata) == set(range(m.SEED_ANGLES))
    assert set(row["sign"] for row in metadata) == {-1, 1}
    for row in metadata:
        scaled = abs(row["z"]) / 2.0
        if row["band"] == "shoulder":
            assert m.SHOULDER_WINDOW[0] < scaled < m.SHOULDER_WINDOW[1]
            assert not (m.TIP_WINDOW[0] < scaled < m.TIP_WINDOW[1])
        else:
            assert m.TIP_WINDOW[0] < scaled < m.TIP_WINDOW[1]
            assert not (m.SHOULDER_WINDOW[0] < scaled < m.SHOULDER_WINDOW[1])


def _manufactured_inward_swirl(points, time):
    del time
    points = np.asarray(points, dtype=float)
    x, y, z = points[:, 0], points[:, 1], points[:, 2]
    a = 0.03
    omega = 0.12
    b = 0.02
    return np.column_stack((-a * x - omega * y, omega * x - a * y, b * z))


def test_measure_offmidplane_paths_has_expected_population_and_active_exposure():
    result = m.measure_offmidplane_paths(_manufactured_inward_swirl)
    assert result["path_count"] == 48
    assert result["paired_material_line_count"] == 24
    assert result["inward_path_count"] == 48
    assert result["pair_axial_separation_growth_count"] == 24
    assert result["pair_axial_separation_shrink_count"] == 0
    assert result["mean_absolute_turns"] > 0.0
    assert result["by_seed_band"]["shoulder"]["path_count"] == 24
    assert result["by_seed_band"]["tip"]["path_count"] == 24
    assert result["by_seed_band"]["shoulder"]["mean_shoulder_window_sample_fraction"] > 0.95
    assert result["by_seed_band"]["tip"]["mean_tip_window_sample_fraction"] > 0.95


def test_identity_comparison_is_zero_without_acceptance_semantics():
    result = m.measure_offmidplane_paths(_manufactured_inward_swirl)
    comparison = m.compare_measurements(result, result)
    assert comparison["mean_absolute_turns_relative"] == 0.0
    assert comparison["maximum_absolute_turns_relative"] == 0.0
    assert comparison["radial_contraction_magnitude_relative"] == 0.0
    assert comparison["mean_pair_axial_separation_change_relative"] == 0.0
    assert comparison["mean_abs_z_change_delta"] == 0.0
    assert comparison["inward_path_count_delta"] == 0
    assert comparison["pair_growth_count_delta"] == 0
    assert comparison["pair_shrink_count_delta"] == 0


def test_truth_boundary_remains_fail_closed():
    false_keys = (
        "canonical_velocity_changed",
        "production_candidate_selected",
        "production_taper_selected",
        "production_compensation_selected",
        "pressure_or_force_refit_in_this_increment",
        "held_out_pde_residual_recomputed_in_this_increment",
        "parent_pde_receipt_transferred",
        "upstream_eulerian_receipt_promoted_to_pde_truth",
        "public_image_used_as_numeric_target",
        "hidden_openai_time_camera_seed_or_velocity_used",
        "visual_acceptance_threshold_defined",
        "trajectory_acceptance_threshold_defined",
        "visualization_ready",
        "visual_correspondence_verified",
        "pde_validated",
        "source_correspondence_verified",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    )
    for key in false_keys:
        assert m.TRUTH_BOUNDARY[key] is False
    assert m.TRUTH_BOUNDARY["comparison_is_descriptive_not_acceptance"] is True
