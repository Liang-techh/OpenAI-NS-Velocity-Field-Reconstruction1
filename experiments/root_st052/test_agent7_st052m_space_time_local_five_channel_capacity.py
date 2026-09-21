from __future__ import annotations

import math

import agent7_st052m_space_time_local_five_channel_capacity as audit


def test_source_lock_and_preregistered_protocol_are_frozen() -> None:
    audit._assert_source_lock()
    assert audit.TASK_ID == "CR003-ST052M-SPACE-TIME-LOCAL-FIVE-CHANNEL-CAPACITY-126"
    assert audit.PREREG_ISSUE == 1057
    assert audit.SOURCE_PARENT_PR == 1050
    assert audit.SOURCE_PARENT_HEAD == "3b94159939ddc48811069a34bc84e5fe18e2fbb2"
    assert audit.REGIONS == ("inner_core", "outer_radial", "midplane", "axial_tip")
    assert audit.TIME_SLICES == (0.25, 0.375, 0.50, 0.625, 0.75)
    assert audit.ACTIVE_RELATIVE_NORM_MIN == 1.0e-12
    assert audit.CONDITION_MAX == 25.0
    assert audit.COSINE_MAX_ABS == 0.995
    assert audit.LEAKAGE_MAX == 1.0e-12


def test_all_twenty_cells_are_materialized_with_expected_start_obstruction() -> None:
    receipt = audit.space_time_local_capacity()
    assert receipt["cell_count"] == 20
    assert receipt["region_order"] == list(audit.REGIONS)
    assert receipt["time_order"] == list(audit.TIME_SLICES)
    assert receipt["unexpected_global_inactive_time_slices"] == []
    assert receipt["all_global_time_activity_semantics_match"] is True

    expected_start = [f"0.250/{region}" for region in audit.REGIONS]
    assert receipt["expected_start_endpoint_obstruction_cells"] == expected_start
    start = receipt["by_time"]["0.250"]
    assert start["global_active_channels"] == []
    assert start["global_inactive_channels"] == list(audit.time_local.regional.joint.CHANNELS)

    for region in audit.REGIONS:
        row = start["regions"][region]
        assert row["active_subspace_rank"] == 0
        assert row["active_subspace_condition"] is None
        assert row["expected_start_endpoint_obstruction"] is True
        assert row["local_representation_obstruction"] is True
        assert row["passes_preregistered_cell_gate"] is True


def test_nonendpoint_cells_expose_rank_conditioning_and_component_fingerprints() -> None:
    receipt = audit.space_time_local_capacity()
    channels = list(audit.time_local.regional.joint.CHANNELS)
    expected_counts = {"0.375": 5, "0.500": 5, "0.625": 5, "0.750": 4}

    for key, expected_count in expected_counts.items():
        time_row = receipt["by_time"][key]
        assert len(time_row["global_active_channels"]) == expected_count
        for region in audit.REGIONS:
            row = time_row["regions"][region]
            assert row["point_count"] > 0
            assert 0.0 < row["point_fraction"] <= 1.0
            assert 0 <= row["active_subspace_rank"] <= expected_count
            assert set(row["raw_regional_column_norms"]) == set(channels)
            assert set(row["regional_norm_over_full_spacetime_channel_norm"]) == set(channels)
            assert set(row["component_fingerprints"]) == set(channels)
            assert set(row["locally_negligible_globally_active_channels"]).issubset(
                set(time_row["global_active_channels"])
            )
            assert len(row["top_two_channels_by_raw_regional_response"]) == 2
            assert isinstance(row["active_subspace_pass"], bool)
            assert isinstance(row["selectivity_guard_pass"], bool)
            assert isinstance(row["local_representation_obstruction"], bool)
            assert isinstance(row["passes_preregistered_cell_gate"], bool)
            if row["active_subspace_condition"] is not None:
                assert math.isfinite(float(row["active_subspace_condition"]))
            if row["max_abs_pairwise_cosine_active_subspace"] is not None:
                assert math.isfinite(float(row["max_abs_pairwise_cosine_active_subspace"]))
            for value in row["regional_norm_over_full_spacetime_channel_norm"].values():
                assert math.isfinite(float(value)) and float(value) >= 0.0


def test_scientific_failure_is_receipted_not_hidden() -> None:
    receipt = audit.space_time_local_capacity()
    assert isinstance(receipt["local_representation_obstruction_cells"], list)
    assert isinstance(receipt["unexpected_capacity_failure_cells"], list)
    assert isinstance(receipt["locally_negligible_channels_by_cell"], dict)
    assert isinstance(receipt["all_nonendpoint_cells_pass_preregistered_gates"], bool)
    expected_start = set(receipt["expected_start_endpoint_obstruction_cells"])
    assert not expected_start.intersection(receipt["unexpected_capacity_failure_cells"])


def test_report_keeps_basis_growth_and_scientific_truth_fail_closed() -> None:
    report = audit.build_report()
    protocol = report["frozen_protocol"]
    assert protocol["cell_count"] == 20
    assert protocol["new_basis_dimension"] == 0
    assert protocol["coefficient_selected"] is None
    assert protocol["public_openai_numeric_target"] is None
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
