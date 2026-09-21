from __future__ import annotations

import math

import numpy as np

import agent7_st052m_regional_five_channel_capacity_localization as audit


def test_source_lock_and_preregistered_protocol_are_frozen() -> None:
    audit._assert_source_lock()
    assert audit.TASK_ID == "CR003-ST052M-REGIONAL-FIVE-CHANNEL-CAPACITY-124"
    assert audit.PREREG_ISSUE == 1041
    assert audit.SOURCE_PARENT_PR == 1033
    assert audit.SOURCE_PARENT_HEAD == "c692a445e4875bcc8965d649804a3b1a33a0ecec"
    assert audit.REGIONS == ("inner_core", "outer_radial", "midplane", "axial_tip")
    assert audit.MIN_REGION_FRACTION == 0.25
    assert audit.RANK_TARGET == 5
    assert audit.CONDITION_MAX == 25.0
    assert audit.COSINE_MAX_ABS == 0.995


def test_median_regions_are_deterministic_and_nonvacuous() -> None:
    points, _, _ = audit.joint._build_tangent_blocks()
    masks, thresholds = audit._region_masks(points)
    assert set(masks) == set(audit.REGIONS)
    assert math.isfinite(thresholds["r50"]) and thresholds["r50"] >= 0.0
    assert math.isfinite(thresholds["abs_z50"]) and thresholds["abs_z50"] >= 0.0
    for name in audit.REGIONS:
        mask = masks[name]
        assert mask.dtype == np.bool_
        assert mask.shape == (points.shape[0],)
        assert float(np.mean(mask)) >= audit.MIN_REGION_FRACTION


def test_regional_receipt_is_complete_even_if_scientific_gate_fails() -> None:
    receipt = audit.regional_capacity_localization()
    assert receipt["region_order"] == list(audit.REGIONS)
    assert set(receipt["regions"]) == set(audit.REGIONS)
    assert set(receipt["full_cloud_channel_norms"]) == set(audit.joint.CHANNELS)
    assert all(
        math.isfinite(float(value)) and float(value) > 0.0
        for value in receipt["full_cloud_channel_norms"].values()
    )
    assert isinstance(receipt["all_regions_pass"], bool)
    assert isinstance(receipt["all_region_selectivity_guards_pass"], bool)

    for name in audit.REGIONS:
        region = receipt["regions"][name]
        assert region["point_fraction"] >= audit.MIN_REGION_FRACTION
        assert 0 <= region["normalized_column_rank"] <= 5
        assert len(region["singular_values"]) == 5
        assert len(region["pairwise_cosines"]) == 10
        assert len(region["top_two_channels_by_raw_norm"]) == 2
        assert set(region["raw_column_norms"]) == set(audit.joint.CHANNELS)
        assert set(region["channel_nonzero"]) == set(audit.joint.CHANNELS)
        assert set(region["regional_energy_fraction_of_full_channel"]) == set(audit.joint.CHANNELS)
        assert set(region["component_fingerprints"]) == set(audit.joint.CHANNELS)
        assert isinstance(region["passes"], bool)
        if region["normalized_condition"] is not None:
            assert math.isfinite(float(region["normalized_condition"]))
        if region["max_abs_pairwise_cosine"] is not None:
            assert math.isfinite(float(region["max_abs_pairwise_cosine"]))
        for value in region["regional_energy_fraction_of_full_channel"].values():
            assert -1.0e-14 <= float(value) <= 1.0 + 1.0e-12


def test_regional_component_selectivity_remains_representation_consistent() -> None:
    receipt = audit.regional_capacity_localization()
    for region in receipt["regions"].values():
        fingerprints = region["component_fingerprints"]
        tor = fingerprints["toroidal_swirl"]["cylindrical_component_rms"]
        assert abs(tor["u_r"]) <= audit.LEAKAGE_MAX
        assert abs(tor["u_z"]) <= audit.LEAKAGE_MAX
        for name in ("amplitude", "radial_shape", "axial_turnover", "temporal_curvature"):
            cyl = fingerprints[name]["cylindrical_component_rms"]
            assert abs(cyl["u_theta"]) <= audit.LEAKAGE_MAX
        assert region["selectivity_guard_pass"] is True


def test_report_keeps_basis_growth_and_scientific_truth_fail_closed() -> None:
    report = audit.build_report()
    assert report["frozen_protocol"]["new_basis_dimension"] == 0
    assert report["frozen_protocol"]["coefficient_selected"] is None
    assert report["frozen_protocol"]["public_openai_numeric_target"] is None
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
