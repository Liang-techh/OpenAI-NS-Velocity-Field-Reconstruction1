from __future__ import annotations

import copy

import agent7_st052m_public_source_temporal_persistence as audit


def _records():
    times = audit.TIMES
    linear = [0.0, 0.005, 0.010, 0.015, 0.020]
    child = [0.0, 0.00502, 0.01003, 0.01504, 0.02005]
    out = {}
    for t, l, c in zip(times, linear, child):
        out[f"{t:.3f}"] = {
            "time": t,
            "relative_to_control": {
                "linear": {"full_aspect_ratio": l},
                "child": {"full_aspect_ratio": c},
            },
            "child_minus_linear": {
                "full_aspect_ratio": c - l,
            },
        }
    return out


def test_temporal_preference_passes_only_with_persistent_advantage():
    d = audit.temporal_decision(_records())
    assert d["start_identity_guard"] is True
    assert d["challenger_aspect_advantage_all_post_start_times"] is True
    assert d["linear_increasing_elongation_guard"] is True
    assert d["child_increasing_elongation_guard"] is True
    assert d["public_source_temporal_axial_stretch_preferred"] is True


def test_one_nonpositive_post_start_increment_fails_closed():
    r = _records()
    r["0.625"]["relative_to_control"]["child"]["full_aspect_ratio"] = 0.015
    r["0.625"]["child_minus_linear"]["full_aspect_ratio"] = 0.0
    d = audit.temporal_decision(r)
    assert d["challenger_aspect_advantage_all_post_start_times"] is False
    assert d["public_source_temporal_axial_stretch_preferred"] is False


def test_nonmonotone_linear_or_child_sequence_fails_closed():
    r = _records()
    r["0.625"]["relative_to_control"]["linear"]["full_aspect_ratio"] = 0.009
    r["0.625"]["child_minus_linear"]["full_aspect_ratio"] = (
        r["0.625"]["relative_to_control"]["child"]["full_aspect_ratio"] - 0.009
    )
    d = audit.temporal_decision(r)
    assert d["linear_increasing_elongation_guard"] is False
    assert d["public_source_temporal_axial_stretch_preferred"] is False

    r = _records()
    r["0.625"]["relative_to_control"]["child"]["full_aspect_ratio"] = 0.0099
    r["0.625"]["child_minus_linear"]["full_aspect_ratio"] = -0.0051
    d = audit.temporal_decision(r)
    assert d["child_increasing_elongation_guard"] is False
    assert d["public_source_temporal_axial_stretch_preferred"] is False


def test_start_identity_guard_uses_frozen_tolerance():
    r = _records()
    r["0.250"]["relative_to_control"]["child"]["full_aspect_ratio"] = 2.0e-12
    d = audit.temporal_decision(r)
    assert d["start_identity_guard"] is False
    assert d["public_source_temporal_axial_stretch_preferred"] is False


def test_nonfinite_value_fails_closed_and_truth_boundary_stays_false():
    r = _records()
    r["0.500"]["relative_to_control"]["child"]["extra"] = float("nan")
    d = audit.temporal_decision(r)
    assert d["all_morphology_values_finite"] is False
    assert d["public_source_temporal_axial_stretch_preferred"] is False

    for key in (
        "canonical_velocity_changed",
        "saved_velocity_changed",
        "production_candidate_selected",
        "new_spatial_basis_added",
        "new_temporal_basis_added",
        "witness_retuned",
        "optimization_performed",
        "public_image_numeric_target_used",
        "pixel_similarity_objective_used",
        "visual_correspondence_verified",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    ):
        assert audit.TRUTH[key] is False


def test_secondary_fields_do_not_change_primary_decision():
    r = _records()
    baseline = audit.temporal_decision(r)
    r2 = copy.deepcopy(r)
    for rec in r2.values():
        rec["secondary_tip_metric"] = -999.0
    assert audit.temporal_decision(r2) == baseline
