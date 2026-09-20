from __future__ import annotations

import numpy as np

import agent7_st052m_total_shape_coordinate_sensitivity as audit


def test_local_family_is_frozen_and_live_k_replays_reservoir_channel() -> None:
    audit._assert_source_lock()
    assert audit.K_VALUES == (0.90, 1.00, 1.10)
    assert audit.K_LIVE == 1.0
    pts = np.array(
        [
            [0.55, 0.0, 1.25],
            [0.90, 0.10, -1.55],
            [1.20, -0.20, 1.75],
            [1.45, 0.0, 1.90],
            [1.61, 0.0, 1.55],
        ],
        dtype=float,
    )
    live = audit.reservoir_correction_k(pts, 1.0)
    source = audit.reservoir.reservoir_correction(pts)
    assert np.array_equal(live, source)
    assert np.max(np.abs(audit.reservoir_correction_k(pts[-1:], 0.9))) == 0.0
    assert np.max(np.abs(audit.reservoir_correction_k(pts[-1:], 1.1))) == 0.0


def test_positive_weighted_effective_radius_is_well_defined() -> None:
    profile = np.linspace(1.0, -1.0, len(audit.total_axial.RADII))
    value = audit._positive_weighted_effective_radius(profile)
    assert value is not None
    assert float(audit.total_axial.RADII[0]) <= value <= float(audit.total_axial.RADII[-1])
    assert audit._positive_weighted_effective_radius(-np.ones_like(profile)) is None


def test_recalibrated_total_shape_tangent_is_independent() -> None:
    parent_fn, alphas, calibrations, children = audit.build_local_family()
    assert set(alphas) == {0.9, 1.0, 1.1}
    assert all(np.isfinite(v) and v != 0.0 for v in alphas.values())
    assert calibrations[1.0]["relative_mismatch_vs_exact_775"] <= audit.ALPHA_REPLAY_TOL
    assert all(row["tip_path_count"] == 24 for row in calibrations.values())
    assert all(row["fresh_714_data_used"] is False for row in calibrations.values())
    diagnostic = audit.total_family_identifiability(parent_fn, alphas, children)
    assert diagnostic["normalized_column_rank"] == 2
    assert abs(diagnostic["normalized_column_cosine"]) < audit.IDENTIFIABILITY_COSINE_MAX_ABS
    assert diagnostic["normalized_two_column_condition"] <= audit.IDENTIFIABILITY_CONDITION_MAX
    assert diagnostic["passes"] is True


def test_truth_boundary_constants() -> None:
    assert audit.PREREG_ISSUE == 809
    assert audit.SOURCE_PARENT_PR == 801
    assert audit.RELATED_FRESH_HOLDOUT_PR == 714
    assert audit.TRUTH["candidate_velocity_changed"] is False
    assert audit.TRUTH["basis_dimension_changed"] is False
    assert audit.TRUTH["candidate_k_selected"] is False
    assert audit.TRUTH["parameter_scan_performed"] is False
    assert audit.TRUTH["fresh_714_path_data_used"] is False
    assert audit.TRUTH["held_out_pde_residual_evaluated"] is False
    assert audit.TRUTH["public_image_numeric_target_used"] is False
    assert audit.TRUTH["closer_visualization_delivery_established"] is False
    assert audit.TRUTH["visual_correspondence_verified"] is False
    assert audit.TRUTH["pde_validated"] is False
    assert audit.TRUTH["openai_field_identified"] is False
