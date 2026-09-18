import numpy as np

from openai_ns_reconstruction.kokuno_source_torus_binding_independent import (
    EXPONENT_MUTATION_RELATIVE_RMS_FLOOR,
    FINE_RELATIVE_RMS_GUARD,
    MIN_REFINEMENT_RATIO_GUARD,
    PHASE_EMBEDDING_MAX_GUARD,
    SIGN_MUTATION_RELATIVE_RMS_FLOOR,
    SOURCE_BINDING_RELATIVE_GUARD,
    generate_report,
)


def test_source_torus_binding_independent_audit_passes_frozen_guards():
    report = generate_report()
    assert report["local_structural_preflight_passed"] is True
    summary = report["summary"]
    assert summary["worst_phase_embedding_max_abs_error"] <= PHASE_EMBEDDING_MAX_GUARD
    assert summary["worst_source_binding_relative_error"] <= SOURCE_BINDING_RELATIVE_GUARD
    assert summary["worst_finest_primary_relative_rms"] <= FINE_RELATIVE_RMS_GUARD
    assert summary["worst_refinement_ratio"] >= MIN_REFINEMENT_RATIO_GUARD
    assert summary["weakest_wrong_sign_mutation_relative_rms"] >= SIGN_MUTATION_RELATIVE_RMS_FLOOR
    assert summary["weakest_exponent_mutation_relative_rms"] >= EXPONENT_MUTATION_RELATIVE_RMS_FLOOR


def test_source_torus_binding_audit_exercises_seams_axis_and_h_perturbations():
    report = generate_report()
    assert len(report["cases"]) == 3
    assert len({row["h"] for row in report["cases"]}) == 3
    for row in report["cases"]:
        assert row["observed_seam_count"] >= row["designed_seam_count"]
        assert row["minimum_radius"] > 3.0 * max(report["sample_contract"]["resolution_ladder"])
        assert row["binding"] == "corrected_source_displayed_constants"
        assert len(row["ladder"]) == 3
        assert np.all(np.isfinite(row["refinement_ratios"]))


def test_source_torus_binding_audit_keeps_pde_truth_boundary_false():
    report = generate_report()
    gates = report["formal_gates"]
    truth = report["truth_boundary"]
    assert gates["normalized_momentum_max_and_L2"] == 1.0e-3
    assert gates["divergence_max_and_L2"] == 1.0e-5
    assert gates["formal_full_domain_pde_gate_assessed"] is False
    assert gates["pde_validated"] is False
    assert truth["actual_positive_order_background_bound"] is False
    assert truth["actual_auxiliary_torus_mode_family_bound"] is False
    assert truth["public_xyz_t_velocity_correction_materialized"] is False
    assert truth["full_leading_oscillatory_correction_composite_available"] is False
