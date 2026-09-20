from __future__ import annotations

import numpy as np

import agent7_st052m_streamfunction_vector_potential_equivalence as audit


def test_symbolic_axisymmetric_equivalence_is_exact() -> None:
    audit._assert_source_lock()
    result = audit.symbolic_equivalence()
    assert result["exact_identity"] is True
    assert result["u_r_difference"] == "0"
    assert result["u_z_difference"] == "0"


def test_live_outer_reservoir_matches_streamfunction_realization() -> None:
    result = audit.live_reservoir_equivalence()
    assert result["point_count"] > 100
    assert result["axis_point_count"] > 0
    assert result["near_axis_point_count"] > result["axis_point_count"]
    assert result["outside_support_point_count"] > 0
    assert result["max_component_abs_difference"] <= audit.EQUIVALENCE_ABS_MAX_GATE
    assert result["axis_finite"] is True
    assert result["near_axis_finite"] is True
    assert result["axis_transverse_max_abs"] <= audit.EQUIVALENCE_ABS_MAX_GATE
    assert result["outside_support_streamfunction_max_abs"] <= audit.EQUIVALENCE_ABS_MAX_GATE
    assert result["passes"] is True


def test_manufactured_field_locks_sign_and_factor_convention() -> None:
    result = audit.manufactured_convention_check()
    assert result["max_component_abs_difference"] <= audit.EQUIVALENCE_ABS_MAX_GATE
    assert result["passes"] is True


def test_vector_and_streamfunction_paths_are_nontrivial_and_finite() -> None:
    pts = audit._deterministic_points()
    vector = audit.vector_potential_velocity(pts)
    stream = audit.streamfunction_velocity(pts)
    assert vector.shape == stream.shape == pts.shape
    assert np.all(np.isfinite(vector))
    assert np.all(np.isfinite(stream))
    assert np.linalg.norm(vector) > 0.0
    assert np.linalg.norm(stream) > 0.0


def test_truth_boundary_and_routing_are_fail_closed() -> None:
    report = audit.build_report()
    assert report["prereg_issue"] == 827
    assert report["source_parent"] == {
        "pr": 818,
        "head": "c5acc7f8f65b64885ec5660ed9c991696668e300",
    }
    assert report["decision"]["frozen_audit_passes"] is True
    assert report["decision"]["streamfunction_switch_adds_basis_dimension"] is False
    assert report["decision"]["streamfunction_switch_adds_morphology_degree_of_freedom"] is False
    assert report["decision"]["actual_velocity_changed"] is False
    assert report["decision"]["closer_visualization_delivery_established"] is False
    assert audit.TRUTH["candidate_velocity_changed"] is False
    assert audit.TRUTH["basis_dimension_changed"] is False
    assert audit.TRUTH["new_spatial_basis_added"] is False
    assert audit.TRUTH["new_temporal_basis_added"] is False
    assert audit.TRUTH["parameter_scan_performed"] is False
    assert audit.TRUTH["fresh_714_path_data_used"] is False
    assert audit.TRUTH["held_out_pde_residual_evaluated"] is False
    assert audit.TRUTH["public_image_numeric_target_used"] is False
    assert audit.TRUTH["closer_visualization_delivery_established"] is False
    assert audit.TRUTH["visual_correspondence_verified"] is False
    assert audit.TRUTH["pde_validated"] is False
    assert audit.TRUTH["openai_field_identified"] is False
