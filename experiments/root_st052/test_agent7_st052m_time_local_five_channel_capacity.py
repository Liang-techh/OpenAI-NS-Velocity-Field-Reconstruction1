from __future__ import annotations

import math

import agent7_st052m_time_local_five_channel_capacity as audit


def test_source_lock_and_preregistered_protocol_are_frozen() -> None:
    audit._assert_source_lock()
    assert audit.TASK_ID == "CR003-ST052M-TIME-LOCAL-FIVE-CHANNEL-CAPACITY-125"
    assert audit.PREREG_ISSUE == 1048
    assert audit.SOURCE_PARENT_PR == 1042
    assert audit.SOURCE_PARENT_HEAD == "63a114c466041f8e3281e0f80e8268611505828d"
    assert audit.TIME_SLICES == (0.25, 0.375, 0.50, 0.625, 0.75)
    assert audit.ACTIVE_RELATIVE_NORM_MIN == 1.0e-12
    assert audit.CONDITION_MAX == 25.0
    assert audit.COSINE_MAX_ABS == 0.995
    assert audit.EXPECTED_INACTIVE["0.250"] == tuple(audit.regional.joint.CHANNELS)
    assert audit.EXPECTED_INACTIVE["0.750"] == ("temporal_curvature",)


def test_endpoint_inactivity_semantics_are_exposed_not_hidden() -> None:
    receipt = audit.time_local_capacity()
    start = receipt["by_time"]["0.250"]
    end = receipt["by_time"]["0.750"]

    assert start["active_channels"] == []
    assert start["inactive_channels"] == list(audit.regional.joint.CHANNELS)
    assert start["active_subspace_rank"] == 0
    assert start["active_subspace_condition"] is None
    assert start["inactive_set_matches_preregistered_expectation"] is True
    assert receipt["start_endpoint_all_current_corrections_inactive"] is True

    assert end["inactive_channels"] == ["temporal_curvature"]
    assert end["active_channels"] == [
        "amplitude",
        "radial_shape",
        "axial_turnover",
        "toroidal_swirl",
    ]
    assert end["inactive_set_matches_preregistered_expectation"] is True
    assert receipt["end_endpoint_temporal_curvature_inactive"] is True


def test_interior_slices_keep_all_five_channels_active() -> None:
    receipt = audit.time_local_capacity()
    for key in ("0.375", "0.500", "0.625"):
        row = receipt["by_time"][key]
        assert row["active_channels"] == list(audit.regional.joint.CHANNELS)
        assert row["inactive_channels"] == []
        assert row["inactive_set_matches_preregistered_expectation"] is True
        assert row["active_channel_count"] == 5
        assert 0 <= row["active_subspace_rank"] <= 5
        assert len(row["singular_values_active_subspace"]) == 5
        assert len(row["pairwise_cosines_active_subspace"]) == 10
        if row["active_subspace_condition"] is not None:
            assert math.isfinite(float(row["active_subspace_condition"]))
        if row["max_abs_pairwise_cosine_active_subspace"] is not None:
            assert math.isfinite(float(row["max_abs_pairwise_cosine_active_subspace"]))


def test_receipt_is_complete_even_if_scientific_rank_gate_fails() -> None:
    receipt = audit.time_local_capacity()
    assert receipt["time_order"] == list(audit.TIME_SLICES)
    assert set(receipt["full_spacetime_channel_norms"]) == set(audit.regional.joint.CHANNELS)
    assert isinstance(receipt["all_time_slices_match_preregistered_capacity_semantics"], bool)
    assert isinstance(receipt["failing_time_slices"], list)
    assert isinstance(receipt["unexpected_inactive_time_slices"], list)

    for key, row in receipt["by_time"].items():
        assert set(row["raw_column_norms"]) == set(audit.regional.joint.CHANNELS)
        assert set(row["relative_norm_of_full_spacetime_channel"]) == set(
            audit.regional.joint.CHANNELS
        )
        assert set(row["component_fingerprints"]) == set(audit.regional.joint.CHANNELS)
        assert row["selectivity_guard_pass"] is True
        assert isinstance(row["active_subspace_pass"], bool)
        assert isinstance(row["passes_preregistered_time_slice_gate"], bool)
        for value in row["relative_norm_of_full_spacetime_channel"].values():
            assert math.isfinite(float(value)) and float(value) >= 0.0


def test_report_keeps_basis_growth_and_scientific_truth_fail_closed() -> None:
    report = audit.build_report()
    assert report["frozen_protocol"]["new_basis_dimension"] == 0
    assert report["frozen_protocol"]["coefficient_selected"] is None
    assert report["frozen_protocol"]["public_openai_numeric_target"] is None
    assert report["decision"]["start_endpoint_repair_available_in_current_correction_family"] is False
    assert report["decision"]["endpoint_preserving_temporal_channel_can_change_t075"] is False
    assert report["decision"]["additional_basis_dimension_justified_by_this_audit"] is False
    assert report["decision"]["actual_velocity_changed"] is False
    assert report["decision"]["direct_visualization_fingerprint_improvement"] == 0.0
    assert report["decision"]["closer_visualization_delivery_established"] is False
    assert report["held_out_pde_residual"]["evaluated"] is False
    assert report["held_out_pde_residual"]["st006_comparison_performed"] is False
    assert report["candidate_velocity_changed"] is False
    assert report["basis_dimension_changed"] is False
    assert report["pressure_or_force_changed"] is False
    assert report["visual_correspondence_verified"] is False
    assert report["pde_validated"] is False
    assert report["paper_exact"] is False
    assert report["openai_field_identified"] is False
