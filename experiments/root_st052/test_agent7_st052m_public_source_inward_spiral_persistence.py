from __future__ import annotations

import numpy as np

import agent7_st052m_public_source_inward_spiral_persistence as audit


def _records(*, child_turns_gain=0.01, child_radius=-0.02, child_count=48):
    keys = ("0.250-0.375", "0.375-0.500", "0.500-0.625", "0.625-0.750")
    control_segments = {}
    linear_segments = {}
    child_segments = {}
    for i, key in enumerate(keys):
        cturn = 0.10 + 0.01 * i
        lturn = 0.95 * cturn
        hturn = lturn + child_turns_gain * cturn
        base_row = {
            "mean_radius_change": -0.03,
            "contraction_magnitude": 0.03,
            "inward_path_count": 47,
            "mean_absolute_turns": cturn,
            "maximum_absolute_turns": 0.2,
        }
        control_segments[key] = dict(base_row)
        linear_segments[key] = {
            **base_row,
            "mean_radius_change": -0.025,
            "contraction_magnitude": 0.025,
            "inward_path_count": 47,
            "mean_absolute_turns": lturn,
        }
        child_segments[key] = {
            **base_row,
            "mean_radius_change": child_radius,
            "contraction_magnitude": abs(child_radius),
            "inward_path_count": child_count,
            "mean_absolute_turns": hturn,
        }

    # Pick whole-interval values that replay the already-frozen #667 turns
    # increment exactly while preserving the expected inward-count change.
    whole_control_turns = 1.0
    whole_linear_turns = 0.9
    whole_child_turns = 0.9 + audit.WHOLE_TURNS_INCREMENT_EXPECTED
    whole_common = {
        "mean_radius_change": -0.1,
        "contraction_magnitude": 0.1,
        "inward_path_count": 48,
        "maximum_absolute_turns": 1.0,
    }
    return (
        {
            "segments": control_segments,
            "whole_interval": {**whole_common, "mean_absolute_turns": whole_control_turns},
        },
        {
            "segments": linear_segments,
            "whole_interval": {**whole_common, "mean_absolute_turns": whole_linear_turns},
        },
        {
            "segments": child_segments,
            "whole_interval": {**whole_common, "mean_absolute_turns": whole_child_turns},
        },
    )


def test_checkpoint_indices_align_with_frozen_33_samples():
    times = np.linspace(audit.base.TIME_INTERVAL[0], audit.base.TIME_INTERVAL[1], audit.base.OUTPUT_SAMPLES)
    assert audit.base.OUTPUT_SAMPLES == 33
    assert audit.CHECKPOINT_INDICES == (0, 8, 16, 24, 32)
    np.testing.assert_allclose(
        times[list(audit.CHECKPOINT_INDICES)],
        np.asarray(audit.CHECKPOINT_TIMES),
        rtol=0.0,
        atol=1.0e-15,
    )


def test_preregistered_rule_accepts_all_segment_improvement():
    control, linear, child = _records()
    result = audit.persistence_decision(control, linear, child)
    assert result["all_segments_preferred"] is True
    assert result["whole_interval_consistency"]["guard"] is True
    assert result["public_source_temporal_inward_spiral_preferred"] is True


def test_rule_rejects_one_segment_turns_regression():
    control, linear, child = _records()
    child["segments"]["0.500-0.625"]["mean_absolute_turns"] = (
        linear["segments"]["0.500-0.625"]["mean_absolute_turns"] - 1.0e-4
    )
    result = audit.persistence_decision(control, linear, child)
    assert result["segments"]["0.500-0.625"]["segment_preferred"] is False
    assert result["public_source_temporal_inward_spiral_preferred"] is False


def test_rule_rejects_outward_segment_or_inward_count_loss():
    control, linear, child = _records()
    child["segments"]["0.375-0.500"]["mean_radius_change"] = 1.0e-6
    child["segments"]["0.250-0.375"]["inward_path_count"] = 46
    result = audit.persistence_decision(control, linear, child)
    assert result["segments"]["0.375-0.500"]["child_mean_radial_change_negative_guard"] is False
    assert result["segments"]["0.250-0.375"]["inward_path_count_not_lower_guard"] is False
    assert result["public_source_temporal_inward_spiral_preferred"] is False


def test_rule_fail_closes_whole_interval_identity_drift():
    control, linear, child = _records()
    child["whole_interval"]["mean_absolute_turns"] += 1.0e-4
    result = audit.persistence_decision(control, linear, child)
    assert result["whole_interval_consistency"]["guard"] is False
    assert result["public_source_temporal_inward_spiral_preferred"] is False


def test_truth_boundary_and_frozen_witness_are_nonpromotional():
    assert audit.PREREG_ISSUE == 682
    assert audit.SOURCE_WITNESS_PR == 652
    assert audit.SOURCE_TEMPORAL_PR == 587
    assert audit.witness.SWIRL_A == 0.065899695471146
    assert abs(audit.witness.SHOULDER_LAMBDA + 0.02) <= 1.0e-15
    assert audit.TRUTH["new_spatial_basis_added"] is False
    assert audit.TRUTH["new_temporal_basis_added"] is False
    assert audit.TRUTH["witness_retuned"] is False
    assert audit.TRUTH["visual_correspondence_verified"] is False
    assert audit.TRUTH["pde_validated"] is False
    assert audit.TRUTH["openai_field_identified"] is False
