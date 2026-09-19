from __future__ import annotations

import numpy as np

import agent7_st052m_combined_pareto_witness as mod


def test_preregistered_witness_is_frozen_and_not_a_scan():
    assert mod.PREREG_ISSUE == 651
    assert mod.SOURCE_CONE_PR == 642
    assert mod.SOURCE_CONE_HEAD == "f523c221847e1a587ec9af4ec053f6f864ce90c4"
    assert mod.SWIRL_A == 0.0658997
    assert mod.SHOULDER_B == -1.0
    assert mod.SHOULDER_LAMBDA == -0.02
    assert mod.PARETO_TOL == 5.0e-5


def test_fractional_swirl_addition_is_zero_at_early_endpoint():
    pts = np.array(
        [
            [0.8, 0.0, 1.0],
            [0.0, 0.8, -1.0],
            [0.4, 0.4, 1.2],
        ],
        dtype=float,
    )
    phased = np.arange(9, dtype=float).reshape(3, 3)
    out = mod.add_compact_swirl(phased, pts, 0.25, epsilon=0.017)
    assert np.array_equal(out, phased)


def test_fractional_swirl_addition_respects_frozen_compact_support():
    pts = np.array(
        [
            [1.36, 0.0, 1.0],
            [0.8, 0.0, 0.69],
            [0.8, 0.0, 1.46],
        ],
        dtype=float,
    )
    phased = np.zeros((3, 3), dtype=float)
    out = mod.add_compact_swirl(phased, pts, 0.50, epsilon=0.017)
    assert np.max(np.abs(out)) == 0.0


def test_nonlinear_pareto_rule_accepts_tolerated_tradeoff_with_real_gain():
    increments = {
        "aspect_gain_increment": 2.0e-5,
        "tip_thinning_increment": -5.0e-5,
        "turns_fidelity_increment": 7.8e-4,
        "axial_pair_fidelity_increment": 0.0,
    }
    result = mod.nonlinear_pareto_rule(increments, structural_pass=True)
    assert result["all_desirabilities_nonworsening"] is True
    assert result["strict_improvement_any"] is True
    assert result["passes"] is True


def test_nonlinear_pareto_rule_rejects_posthoc_tradeoff_beyond_tolerance():
    increments = {
        "aspect_gain_increment": 2.0e-5,
        "tip_thinning_increment": -5.1e-5,
        "turns_fidelity_increment": 7.8e-4,
        "axial_pair_fidelity_increment": 0.0,
    }
    result = mod.nonlinear_pareto_rule(increments, structural_pass=True)
    assert result["all_desirabilities_nonworsening"] is False
    assert result["passes"] is False


def test_truth_boundary_stays_fail_closed():
    assert mod.TRUTH["canonical_velocity_changed"] is False
    assert mod.TRUTH["saved_velocity_changed"] is False
    assert mod.TRUTH["new_spatial_basis_added"] is False
    assert mod.TRUTH["new_temporal_basis_added"] is False
    assert mod.TRUTH["parameter_grid_scan_performed"] is False
    assert mod.TRUTH["optimization_performed"] is False
    assert mod.TRUTH["held_out_pde_residual_evaluated"] is False
    assert mod.TRUTH["visual_correspondence_verified"] is False
    assert mod.TRUTH["pde_validated"] is False
    assert mod.TRUTH["openai_field_identified"] is False
