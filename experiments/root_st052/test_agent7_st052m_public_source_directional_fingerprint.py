from __future__ import annotations

import copy

import agent7_st052m_public_source_directional_fingerprint as audit


def _grid_report():
    return {
        "grid_robustness_decision": {
            "aspect_increment_by_resolution": {
                "33": 2.0e-5,
                "41": 2.1e-5,
                "49": 2.2e-5,
            }
        }
    }


def _witness_report():
    return {
        "path_relative_to_control": {
            "linear": {
                "mean_absolute_turns_relative": -0.009,
                "inward_path_count_delta": 0,
                "contraction_magnitude_relative": 0.167,
            },
            "child": {
                "mean_absolute_turns_relative": -0.008,
                "inward_path_count_delta": 0,
                "contraction_magnitude_relative": 0.167 + 1.0e-12,
            },
        },
        "structure_preflight": {
            "support_max_abs": 0.0,
            "divergence_fd_max": 1.0e-8,
        },
    }


def test_public_source_directional_pass_requires_all_primary_guards():
    d = audit.directional_decision(_grid_report(), _witness_report())
    assert d["axial_stretching_guard"] is True
    assert d["inward_spiraling_guard"] is True
    assert d["central_contraction_retained_guard"] is True
    assert d["structural_guard"] is True
    assert d["public_source_directionally_preferred"] is True


def test_one_bad_grid_rejects_axial_stretching():
    g = _grid_report()
    g["grid_robustness_decision"]["aspect_increment_by_resolution"]["49"] = -1.0e-12
    d = audit.directional_decision(g, _witness_report())
    assert d["axial_stretching_guard"] is False
    assert d["public_source_directionally_preferred"] is False


def test_spiral_and_contraction_guards_fail_closed():
    w = _witness_report()
    w["path_relative_to_control"]["child"]["mean_absolute_turns_relative"] = -0.010
    d = audit.directional_decision(_grid_report(), w)
    assert d["inward_spiraling_guard"] is False
    assert d["public_source_directionally_preferred"] is False

    w = _witness_report()
    w["path_relative_to_control"]["child"]["contraction_magnitude_relative"] += 2.0e-10
    d = audit.directional_decision(_grid_report(), w)
    assert d["central_contraction_retained_guard"] is False
    assert d["public_source_directionally_preferred"] is False


def test_structure_and_truth_boundary_are_fail_closed():
    w = _witness_report()
    w["structure_preflight"]["support_max_abs"] = 1.0e-15
    d = audit.directional_decision(_grid_report(), w)
    assert d["support_exact_zero_guard"] is False
    assert d["public_source_directionally_preferred"] is False

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


def test_tip_thinning_is_not_a_primary_decision_input():
    # The frozen public-source decision consumes only aspect, path contraction/
    # spiraling and structure; autonomous tip-thinning diagnostics are excluded.
    g = _grid_report()
    w = _witness_report()
    baseline = audit.directional_decision(g, w)
    w2 = copy.deepcopy(w)
    w2["arbitrary_tip_thinning_diagnostic"] = -999.0
    mutated = audit.directional_decision(g, w2)
    assert mutated == baseline
