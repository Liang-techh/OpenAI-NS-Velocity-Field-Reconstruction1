import numpy as np

from openai_ns_reconstruction.kokuno_physical_evaluation_independent import (
    AXIS_NEAR_RADII,
    BASE_HEAD,
    FINE_RELATIVE_RMS_GUARD,
    FORMAL_DIVERGENCE_GATE,
    FORMAL_MOMENTUM_GATE,
    MIN_REFINEMENT_RATIO_GUARD,
    OMIT_TORUS_MUTATION_RELATIVE_RMS_FLOOR,
    SEAM_RADII,
    STEPS,
    generate_report,
)


def test_independent_physical_evaluation_preflight_passes_frozen_guards():
    report = generate_report()
    summary = report["summary"]
    assert report["base_head"] == BASE_HEAD
    assert tuple(report["sample_contract"]["resolution_ladder"]) == STEPS
    assert report["sample_contract"]["uses_agent2_fd4"] is False
    assert summary["worst_finest_primary_relative_rms"] <= FINE_RELATIVE_RMS_GUARD
    assert summary["worst_refinement_ratio"] >= MIN_REFINEMENT_RATIO_GUARD
    assert (
        summary["weakest_omit_torus_mutation_relative_rms"]
        >= OMIT_TORUS_MUTATION_RELATIVE_RMS_FLOOR
    )
    assert summary["local_structural_preflight_passed"] is True


def test_fresh_off_grid_axis_near_and_mod1_seam_samples_are_exercised():
    report = generate_report()
    for case in report["cases"]:
        assert case["random_off_grid_count"] == 18
        assert case["axis_near_count"] == len(AXIS_NEAR_RADII)
        assert case["designed_torus_seam_count"] == len(SEAM_RADII)
        assert case["observed_torus_seam_count"] >= len(SEAM_RADII)
        assert case["minimum_radius"] > 3.0 * max(STEPS)
        assert len(case["primary_relative_rms_ladder"]) == 3
        assert len(case["refinement_ratios"]) == 2
        assert np.all(np.isfinite(case["primary_relative_rms_ladder"]))


def test_omitting_torus_chain_terms_is_strongly_detected():
    report = generate_report()
    for mutation in report["mutation"]:
        assert mutation["omit_torus_grad_r_relative_rms"] >= OMIT_TORUS_MUTATION_RELATIVE_RMS_FLOOR
        assert mutation["omit_torus_grad_t_relative_rms"] >= OMIT_TORUS_MUTATION_RELATIVE_RMS_FLOOR


def test_truth_boundary_keeps_formal_pde_gate_false_and_thresholds_fixed():
    report = generate_report()
    assert report["frozen_guards"]["formal_momentum_gate_unchanged"] == FORMAL_MOMENTUM_GATE
    assert report["frozen_guards"]["formal_divergence_gate_unchanged"] == FORMAL_DIVERGENCE_GATE
    truth = report["truth_boundary"]
    assert truth["auxiliary_torus_pullback_independently_preflighted"] is True
    assert truth["actual_positive_order_background_bound"] is False
    assert truth["actual_auxiliary_torus_mode_family_bound"] is False
    assert truth["public_xyz_t_velocity_correction_materialized"] is False
    assert truth["formal_full_domain_pde_gate_assessed"] is False
    assert truth["pde_validated"] is False
