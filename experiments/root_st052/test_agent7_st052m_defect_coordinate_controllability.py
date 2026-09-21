from __future__ import annotations

import math

import numpy as np

import agent7_st052m_defect_coordinate_controllability as audit


def test_source_lock_and_preregistered_protocol_are_frozen() -> None:
    audit._assert_source_lock()
    assert audit.TASK_ID == "CR003-ST052M-DEFECT-COORDINATE-CONTROLLABILITY-128"
    assert audit.PREREG_ISSUE == 1077
    assert audit.SOURCE_PARENT_PR == 1069
    assert audit.SOURCE_PARENT_HEAD == "a564fbc8bfa6a5c630b5839fa208e41cb7bc72d7"
    assert audit.SUPPORT_RADIUS == 2.0
    assert audit.EPSILONS == (1.0e-3, 5.0e-4)
    assert audit.PRIMARY_EPSILON == 5.0e-4
    assert audit.RANK_TARGET == 5
    assert audit.CONDITION_MAX == 25.0
    assert audit.COSINE_MAX_ABS == 0.995
    assert audit.DERIVATIVE_DRIFT_MAX == 0.05
    assert audit.DERIVATIVE_COSINE_MIN == 0.999
    assert audit.ROW_RESPONSE_FLOOR == 1.0e-12
    assert len(audit.DEFECT_COORDINATES) == 8


def test_defect_coordinates_are_finite_and_have_frozen_sign_semantics() -> None:
    base_fn, _, _ = audit.parent._build_field_and_tangents()
    speed_scale = audit.parent._baseline_speed_scale(base_fn)
    coordinates = audit.defect_coordinate_vector(base_fn, speed_scale)
    assert tuple(coordinates) == audit.DEFECT_COORDINATES
    assert all(math.isfinite(float(value)) for value in coordinates.values())
    assert set(audit.COORDINATE_SEMANTICS) == set(audit.DEFECT_COORDINATES)
    assert "negative means contraction" in audit.COORDINATE_SEMANTICS[
        "core_width_change_025_to_075"
    ]["positive_direction"]
    assert "speed increases" in audit.COORDINATE_SEMANTICS[
        "core_speed_change_025_to_075"
    ]["public_observable_relation"]


def test_coordinate_helpers_use_only_frozen_renderer_independent_observables() -> None:
    base_fn, _, _ = audit.parent._build_field_and_tangents()
    speed_scale = audit.parent._baseline_speed_scale(base_fn)
    obs = audit.parent.observable_vector(base_fn)
    radial = audit._ring_values(obs, "radial_ring", 0.50)
    swirl = audit._ring_values(obs, "swirl_ring", 0.50)
    axial = audit._ring_values(obs, "axial_ring", 0.50)
    assert radial.shape == swirl.shape == axial.shape == (16,)
    assert np.all(np.isfinite(radial))
    assert np.all(np.isfinite(swirl))
    assert np.all(np.isfinite(axial))
    variation = audit._angular_rotation_variation(obs, speed_scale)
    assert math.isfinite(variation) and variation >= 0.0
    assert math.isfinite(audit._core_width(obs, 0.25))
    assert math.isfinite(audit._core_width(obs, 0.50))
    assert math.isfinite(audit._core_width(obs, 0.75))


def test_controllability_receipts_scientific_failure_instead_of_hiding_it() -> None:
    receipt = audit.controllability_audit()
    assert receipt["defect_coordinate_count"] == 8
    assert receipt["defect_coordinate_order"] == list(audit.DEFECT_COORDINATES)
    assert receipt["rank_target"] == 5
    assert 0 <= receipt["normalized_column_rank"] <= 5
    assert isinstance(receipt["compressed_routing_geometry_gate_passes"], bool)
    assert isinstance(receipt["all_derivative_stability_gates_pass"], bool)
    assert isinstance(receipt["all_coordinates_have_numeric_response"], bool)
    assert math.isfinite(float(receipt["baseline_rms_speed_scale"]))
    assert receipt["baseline_rms_speed_scale"] > 0.0
    if receipt["normalized_condition"] is not None:
        assert math.isfinite(float(receipt["normalized_condition"]))
    if receipt["max_abs_pairwise_cosine"] is not None:
        assert math.isfinite(float(receipt["max_abs_pairwise_cosine"]))

    channels = tuple(audit.parent._joint().CHANNELS)
    for channel in channels:
        stability = receipt["derivative_stability"][channel]
        assert stability["coarse_epsilon"] == 1.0e-3
        assert stability["fine_epsilon"] == 5.0e-4
        assert isinstance(stability["passes"], bool)
        if stability["relative_column_drift"] is not None:
            assert math.isfinite(float(stability["relative_column_drift"]))
        if stability["coarse_fine_cosine"] is not None:
            assert math.isfinite(float(stability["coarse_fine_cosine"]))

    for name in audit.DEFECT_COORDINATES:
        row = receipt["coordinates"][name]
        assert math.isfinite(float(row["baseline_value"]))
        assert math.isfinite(float(row["fine_raw_row_l2"]))
        assert isinstance(row["structurally_unresponsive_at_numeric_floor"], bool)
        assert set(row["fine_raw_derivative_by_channel"]) == set(channels)
        assert set(row["normalized_column_alignment_by_channel"]) == set(channels)
        fractions = row["normalized_alignment_energy_fraction_by_channel"]
        assert set(fractions) == set(channels)
        assert all(math.isfinite(float(value)) and value >= 0.0 for value in fractions.values())
        if row["dominant_aligned_existing_channel"] is not None:
            assert row["dominant_aligned_existing_channel"] in channels
            assert abs(sum(fractions.values()) - 1.0) <= 2.0e-12


def test_report_keeps_candidate_basis_visual_and_pde_truth_fail_closed() -> None:
    report = audit.build_report()
    protocol = report["frozen_protocol"]
    assert protocol["channels"] == [
        "amplitude",
        "radial_shape",
        "axial_turnover",
        "temporal_curvature",
        "toroidal_swirl",
    ]
    assert protocol["defect_coordinates"] == list(audit.DEFECT_COORDINATES)
    assert protocol["new_basis_dimension"] == 0
    assert protocol["coefficient_selected"] is None
    assert protocol["public_openai_numeric_target"] is None
    assert report["decision"]["additional_basis_dimension_justified_by_this_audit"] is False
    assert report["decision"]["candidate_mutation_authorized_by_this_audit"] is False
    assert report["decision"]["actual_velocity_changed"] is False
    assert report["decision"]["direct_visualization_fingerprint_improvement"] == 0.0
    assert report["decision"]["closer_visualization_delivery_established"] is False
    assert report["held_out_pde_residual"]["evaluated"] is False
    assert report["held_out_pde_residual"]["st006_comparison_performed"] is False
    assert report["candidate_velocity_changed"] is False
    assert report["basis_dimension_changed"] is False
    assert report["pressure_or_force_changed"] is False
    assert report["public_image_numeric_target_used"] is False
    assert report["visual_correspondence_verified"] is False
    assert report["pde_validated"] is False
    assert report["paper_exact"] is False
    assert report["openai_field_identified"] is False
