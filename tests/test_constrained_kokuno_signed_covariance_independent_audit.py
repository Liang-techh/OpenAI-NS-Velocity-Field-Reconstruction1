from openai_ns_reconstruction.kokuno_signed_covariance_independent_audit import (
    FD5_STEPS,
    FROZEN_GUARDS,
    run_audit,
)


def test_signed_covariance_independent_audit_passes_frozen_guards():
    report = run_audit()
    metrics = report["metrics"]
    assert report["local_structural_preflight_passed"] is True
    assert metrics["sample_count"] == 48
    assert metrics["value_relative_error_max"] <= FROZEN_GUARDS["value_relative_error_max"]
    assert metrics["directional_derivative_relative_rms_by_step"][str(FD5_STEPS[-1])] <= FROZEN_GUARDS[
        "finest_directional_derivative_relative_rms_max"
    ]
    assert min(metrics["directional_derivative_refinement_ratios"]) >= FROZEN_GUARDS[
        "directional_derivative_refinement_ratio_min"
    ]
    assert metrics["semantic_sigma_swap_relative_max"] <= FROZEN_GUARDS[
        "semantic_sigma_swap_relative_max"
    ]
    assert metrics["wrong_q_sign_mutation_relative_rms_min"] >= FROZEN_GUARDS[
        "wrong_q_sign_mutation_relative_rms_min"
    ]
    assert metrics["drop_pulse_derivative_mutation_relative_rms_min"] >= FROZEN_GUARDS[
        "drop_pulse_derivative_mutation_relative_rms_min"
    ]
    assert metrics["cone_negative_rejection_fraction"] == 1.0
    assert metrics["direction_gap_negative_rejection_fraction"] == 1.0


def test_signed_covariance_audit_keeps_full_pde_gate_fail_closed():
    report = run_audit()
    truth = report["truth_boundary"]
    gates = report["formal_project_gates"]
    assert truth["reference_signed_covariance_inverse_independently_audited"] is True
    assert truth["actual_positive_order_background_bound"] is False
    assert truth["actual_source_h_sigma_pulse_integrals_bound"] is False
    assert truth["actual_auxiliary_torus_signed_mode_family_bound"] is False
    assert truth["public_source_bound_xyz_t_oscillatory_velocity_ready"] is False
    assert truth["genuinely_independent_second_covariance_column_ready"] is False
    assert truth["correction_ready"] is False
    assert gates["normalized_momentum_max_l2"] == 1.0e-3
    assert gates["normalized_divergence_max_l2"] == 1.0e-5
    assert gates["formal_full_domain_pde_gate_assessed"] is False
    assert gates["pde_validated"] is False
