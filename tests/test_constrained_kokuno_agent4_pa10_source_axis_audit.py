import inspect

from openai_ns_reconstruction.kokuno_agent4_pa10_source_axis_audit import (
    run_audit,
)


def test_independent_source_axis_audit_passes_frozen_guards():
    report = run_audit()
    assert report["failed_guards"] == []
    assert report["source_axis_domain_independent_preflight_passed"] is True
    real = report["real_axis"]
    assert real["nonvacuous_K_and_Hsmall"] is True
    assert real["PA8_independent_sampling_passed"] is True
    assert len(real["levels"]) == 3
    assert real["levels"][-1]["min_abs_H_on_K"] > 0.01
    assert real["levels"][-1]["min_chi_on_K"] > 0.99
    assert real["levels"][-1]["min_Z_on_Hsmall"] > 1.0


def test_complex_roots_and_g_series_are_independent_and_mutation_sensitive():
    report = run_audit()
    complex_report = report["complex_zero_exclusion"]
    assert complex_report["independent_zero_exclusion_passed"] is True
    assert complex_report["min_H_root_distance_over_tube"] > 1.0
    g = report["g_coefficient_norm"]
    assert g["independent_norm_domination_passed"] is True
    assert g["public_to_independent_ratio"] >= 1.0
    assert g["partial_sums"]["64"] <= g["independent_closed_form"]
    assert all(report["mutation"].values())


def test_validator_does_not_call_upstream_certificate_or_axis_state_paths():
    src = inspect.getsource(run_audit)
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
    assert truth["source_Phi_radius_one_ball_norm_machine_bound"] is False
    assert truth["source_Phi_radius_one_ball_lipschitz_machine_bound"] is False
    assert truth["source_R1_R2_machine_bound"] is False
    assert truth["source_operator_M_K_machine_bound"] is False
    assert truth["global_pressure_matched"] is False
    assert truth["global_leading_profile_reconstructed"] is False
    assert truth["heldout_ns_momentum_residual_assessed"] is False
    assert truth["formal_full_domain_pde_gate_assessed"] is False
    assert truth["pde_validated"] is False
