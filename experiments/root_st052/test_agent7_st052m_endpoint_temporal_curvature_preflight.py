from __future__ import annotations

import numpy as np

import agent7_st052m_endpoint_temporal_curvature_preflight as audit


def test_time_shape_is_frozen_and_endpoint_preserving() -> None:
    audit._assert_source_lock()
    assert audit.temporal_curvature(0.25) == 0.0
    assert audit.temporal_curvature(0.50) == 1.0
    assert audit.temporal_curvature(0.75) == 0.0
    assert audit.temporal_curvature(0.375) == audit.temporal_curvature(0.625)
    try:
        audit.temporal_curvature(0.80)
    except ValueError:
        pass
    else:
        raise AssertionError("out-of-window temporal basis evaluation must fail")


def test_temporal_tangent_is_exactly_zero_at_registered_endpoints() -> None:
    _parent, alphas, calibrations, _children = audit.shape_audit.build_local_family()
    assert calibrations[audit.shape_audit.K_LIVE]["fresh_714_data_used"] is False
    alpha = float(alphas[audit.shape_audit.K_LIVE])
    pts = audit.shape_audit._response_points()
    assert np.max(np.abs(audit.temporal_curvature_tangent(pts, 0.25, alpha))) == 0.0
    assert np.max(np.abs(audit.temporal_curvature_tangent(pts, 0.75, alpha))) == 0.0
    assert np.linalg.norm(audit.temporal_curvature_tangent(pts, 0.50, alpha)) > 0.0


def test_three_column_response_is_locally_identifiable() -> None:
    diagnostic = audit.three_column_identifiability()
    assert diagnostic["live_alpha_replay_relative_mismatch"] <= audit.shape_audit.ALPHA_REPLAY_TOL
    assert diagnostic["normalized_column_rank"] == 3
    assert diagnostic["normalized_three_column_condition"] <= audit.IDENTIFIABILITY_CONDITION_MAX
    assert diagnostic["max_abs_pairwise_cosine"] < audit.IDENTIFIABILITY_COSINE_MAX_ABS
    assert diagnostic["passes"] is True


def test_truth_boundary_constants() -> None:
    assert audit.PREREG_ISSUE == 817
    assert audit.SOURCE_PARENT_PR == 810
    assert audit.RELATED_FRESH_HOLDOUT_PR == 714
    assert audit.TRUTH["candidate_velocity_changed"] is False
    assert audit.TRUTH["basis_dimension_changed"] is False
    assert audit.TRUTH["second_spatial_poloidal_basis_added"] is False
    assert audit.TRUTH["temporal_basis_dimension_changed"] is False
    assert audit.TRUTH["temporal_coefficient_selected"] is False
    assert audit.TRUTH["parameter_scan_performed"] is False
    assert audit.TRUTH["fresh_714_path_data_used"] is False
    assert audit.TRUTH["held_out_pde_residual_evaluated"] is False
    assert audit.TRUTH["public_image_numeric_target_used"] is False
    assert audit.TRUTH["closer_visualization_delivery_established"] is False
    assert audit.TRUTH["visual_correspondence_verified"] is False
    assert audit.TRUTH["pde_validated"] is False
    assert audit.TRUTH["openai_field_identified"] is False
