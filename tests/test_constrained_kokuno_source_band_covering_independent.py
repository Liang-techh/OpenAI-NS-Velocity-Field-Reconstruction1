from openai_ns_reconstruction.kokuno_source_band_covering_independent import (
    INTERACTION_BOUND_RELATIVE_GUARD,
    INTERACTION_VIOLATION_GUARD,
    MESH_EXPONENT_MUTATION_RELATIVE_FLOOR,
    OFF_BY_ONE_MUTATION_RELATIVE_FLOOR,
    TRANSITION_LEVEL_MISMATCH_GUARD,
    VALUE_RELATIVE_GUARD,
    generate_report,
)


def test_source_band_covering_independent_audit_passes_frozen_guards():
    report = generate_report()
    assert report["local_structural_preflight_passed"] is True
    assert report["sweep"]["worst_value_relative_error"] <= VALUE_RELATIVE_GUARD
    assert report["sweep"]["covering_level_mismatch_count"] <= TRANSITION_LEVEL_MISMATCH_GUARD
    assert report["transition"]["mismatch_count"] <= TRANSITION_LEVEL_MISMATCH_GUARD
    assert report["interaction"]["worst_bound_formula_relative_error"] <= INTERACTION_BOUND_RELATIVE_GUARD
    assert report["interaction"]["violation_count"] <= INTERACTION_VIOLATION_GUARD
    assert report["mutation"]["weakest_off_by_one_level_relative_change"] >= OFF_BY_ONE_MUTATION_RELATIVE_FLOOR
    assert report["mutation"]["weakest_wrong_mesh_exponent_relative_change"] >= MESH_EXPONENT_MUTATION_RELATIVE_FLOOR


def test_source_band_covering_independent_audit_stresses_floor_transitions_and_interactions():
    report = generate_report()
    assert report["transition"]["case_count"] == 16
    assert report["transition"]["minimum_distance_to_integer"] < 2.0e-12
    assert report["interaction"]["pair_count"] > 1000
    assert report["interaction"]["maximum_level_difference"] >= 1
    assert report["sample_contract"]["uses_agent2_bound_diagnostics"] is False
    assert report["sample_contract"]["uses_agent2_interaction_receipt"] is False


def test_source_band_covering_independent_audit_keeps_pde_truth_boundary_false():
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
