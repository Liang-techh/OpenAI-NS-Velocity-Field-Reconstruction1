import inspect

from openai_ns_reconstruction.kokuno_agent4_pa10_source_axis_audit import (
    AGENT1_HEAD,
    OFFGRID_CELL_DISTANCE_FLOOR,
    PRIOR_A4_HEAD,
    PRIOR_FAILED_WORKFLOW,
    REAL_LEVELS,
    TARGETED_OFFGRID_COUNT,
    _point_satisfies_pa8_separation,
    _targeted_offgrid_region,
    run_audit,
)


def test_repaired_independent_source_axis_audit_passes_frozen_guards():
    report = run_audit()
    assert report["failed_guards"] == []
    assert report["source_axis_domain_independent_preflight_passed"] is True
    assert report["upstream"]["agent1_head"] == AGENT1_HEAD
    real = report["real_axis"]
    assert real["nonvacuous_K_and_Hsmall"] is True
    assert real["PA8_independent_sampling_passed"] is True
    assert [level["count"] for level in real["levels"]] == list(REAL_LEVELS)
    assert real["levels"][-1]["min_abs_H_on_K"] > 0.01
    assert real["levels"][-1]["min_chi_on_K"] > 0.99
    assert real["levels"][-1]["min_Z_on_Hsmall"] > 1.0


def test_targeted_fresh_offgrid_cloud_is_nonvacuous_and_preserves_pa8():
    report = run_audit()
    targeted = report["targeted_offgrid"]
    K = targeted["K_delta"]
    Hsmall = targeted["H_small"]
    assert targeted["targeted_offgrid_PA8_passed"] is True
    assert K["selected_count"] == TARGETED_OFFGRID_COUNT
    assert Hsmall["selected_count"] == TARGETED_OFFGRID_COUNT
    assert K["all_points_strictly_off_dense_grid"] is True
    assert Hsmall["all_points_strictly_off_dense_grid"] is True
    assert K["minimum_distance_from_dense_grid_in_cells"] > OFFGRID_CELL_DISTANCE_FLOOR
    assert Hsmall["minimum_distance_from_dense_grid_in_cells"] > OFFGRID_CELL_DISTANCE_FLOOR
    assert K["max_abs_Z"] <= report["public_parameters"]["delta_star"]
    assert K["min_abs_H"] > 10.0 * report["public_parameters"]["sigma_star"]
    assert K["min_chi"] > 0.99
    assert Hsmall["max_abs_H"] <= 10.0 * report["public_parameters"]["sigma_star"]
    assert Hsmall["min_Z"] > 2.0 * report["public_parameters"]["delta_star"]


def test_protocol_repair_records_old_sampling_failure_without_threshold_retuning():
    report = run_audit()
    repair = report["protocol_repair"]
    assert repair["prior_agent4_head"] == PRIOR_A4_HEAD
    assert repair["prior_failed_workflow"] == PRIOR_FAILED_WORKFLOW
    assert "zero K_delta samples" in repair["prior_failure"]
    assert repair["source_parameters_changed"] is False
    assert repair["scientific_thresholds_changed"] is False
    assert repair["final_project_gates_changed"] is False
    gates = report["frozen_protocol"]["final_project_gates_unchanged"]
    assert gates == {
        "normalized_momentum_max": 1.0e-3,
        "normalized_momentum_L2": 1.0e-3,
        "divergence_max": 1.0e-5,
        "divergence_L2": 1.0e-5,
    }


def test_complex_roots_and_g_series_are_independent_and_mutation_sensitive():
    report = run_audit()
    complex_report = report["complex_zero_exclusion"]
    assert complex_report["independent_zero_exclusion_passed"] is True
    assert complex_report["min_H_root_distance_over_tube"] > 1.0
    g = report["g_coefficient_norm"]
    assert g["independent_norm_domination_passed"] is True
    assert g["public_to_independent_ratio"] >= 1.0
    assert g["partial_sums"]["64"] <= g["independent_closed_form"]
    assert report["mutation"]["complex_tube_x100_detected"] is True
    assert report["mutation"]["g_bound_0p99_detected"] is True


def test_synthetic_pa8_bad_point_is_actually_rejected_by_shared_predicate():
    report = run_audit()
    sigma = report["public_parameters"]["sigma_star"]
    delta = report["public_parameters"]["delta_star"]
    assert _point_satisfies_pa8_separation(0.0, 0.0, sigma=sigma, delta=delta) is False
    assert _point_satisfies_pa8_separation(0.02, 0.0, sigma=sigma, delta=delta) is True
    assert report["mutation"]["synthetic_PA8_bad_point_detected"] is True


def test_validator_does_not_call_upstream_certificate_or_axis_state_paths():
    src = inspect.getsource(run_audit) + inspect.getsource(_targeted_offgrid_region)
    for forbidden in (
        "real_interval_certificate(",
        "complex_domain_certificate(",
        "axis_state(",
        "_real_certificate(",
        "_complex_certificate(",
    ):
        assert forbidden not in src


def test_truth_boundary_keeps_ns_and_missing_leading_objects_closed():
    truth = run_audit()["truth_boundary"]
    assert truth["source_axis_domain_independently_audited"] is True
    assert truth["source_Phi_radius_one_ball_norm_machine_bound"] is False
    assert truth["source_Phi_radius_one_ball_lipschitz_machine_bound"] is False
    assert truth["source_R1_R2_machine_bound"] is False
    assert truth["source_operator_M_K_machine_bound"] is False
    assert truth["global_pressure_matched"] is False
    assert truth["global_leading_profile_reconstructed"] is False
    assert truth["leading_only_ns_residual_assessed"] is False
    assert truth["leading_plus_oscillatory_ns_residual_assessed"] is False
    assert truth["after_correction_ns_residual_assessed"] is False
    assert truth["heldout_ns_momentum_residual_assessed"] is False
    assert truth["formal_full_domain_pde_gate_assessed"] is False
    assert truth["pde_validated"] is False
