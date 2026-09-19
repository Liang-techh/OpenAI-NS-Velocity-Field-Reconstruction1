from __future__ import annotations

from openai_ns_reconstruction.kokuno_pa10_source_r1_averaged_eta_independent_audit import (
    AGENT1_HEAD,
    AGENT1_PR,
    MAX_EXPECTED_PUBLIC_RATIO,
    SLOT_NAME,
    run_independent_audit,
)


def test_independent_audit_passes_without_reusing_agent1_private_path() -> None:
    report = run_independent_audit()
    assert report["agent1_pr"] == AGENT1_PR == 715
    assert report["agent1_exact_head"] == AGENT1_HEAD
    assert report["failed_guards"] == []
    assert report["source_R1_averaged_eta_ordinary_slot_independent_preflight_passed"] is True

    path = report["independent_path"]
    assert path["private_agent1_helpers_called"] is False
    assert path["agent1_receipt_used_as_oracle"] is False
    assert path["training_internal_tensors_read"] is False
    assert path["J2_reapplied"] is False
    assert path["convolution_count"] == 3


def test_public_slot_dominates_tighter_independent_reconstruction_without_shape_drift() -> None:
    report = run_independent_audit()
    ratios = report["public_over_independent_ratios"]
    assert min(ratios.values()) >= 1.0
    assert ratios["slot_norm"] <= float(MAX_EXPECTED_PUBLIC_RATIO)
    assert ratios["slot_lipschitz"] <= float(MAX_EXPECTED_PUBLIC_RATIO)
    assert ratios["product_constant"] <= float(MAX_EXPECTED_PUBLIC_RATIO)
    assert ratios["derivative_after_j2"] <= 1.000000000001
    assert ratios["d_norm"] <= 1.000000000001
    assert ratios["L_inverse"] <= 1.000000000001

    slot = report["independent_values"]["slot_after_J2"]
    assert slot["norm"] > 0.0
    assert slot["lipschitz"] > 0.0
    assert SLOT_NAME == "lambda_inv_d_detaAu_times_Phi"


def test_fresh_offgrid_and_axis_near_fixed_multiplier_stress_is_nonvacuous() -> None:
    report = run_independent_audit()
    stress = report["fresh_offgrid_fixed_multiplier_stress"]
    assert stress["fresh_random_offgrid_points"] == 8192
    assert stress["axis_near_probes"] == [-1e-12, 0.0, 1e-12]
    assert stress["finite"] is True
    assert stress["dominated"] is True
    assert stress["axis_nontrivial"] is True
    assert 0.0 < stress["d_value_over_public_bound"] <= 1.0
    assert 0.0 < stress["L_inverse_value_over_public_bound"] <= 1.0


def test_preregistered_negative_controls_detect_underbounds_and_double_j2() -> None:
    report = run_independent_audit()
    negative = report["negative_controls"]
    assert negative
    assert all(negative.values())
    assert negative["product_constant_0p998_detected"] is True
    assert negative["derivative_79_over_80_detected"] is True
    assert negative["slot_norm_0p99_detected"] is True
    assert negative["slot_lipschitz_0p99_detected"] is True
    assert negative["missing_one_convolution_detected"] is True
    assert negative["duplicate_J2_detected_by_shape_guard"] is True


def test_local_rho_perturbation_is_recorded_without_changing_project_gates() -> None:
    report = run_independent_audit()
    sensitivity = report["rho_local_operator_sensitivity"]
    assert sensitivity["rho_minus_0p1pct"] < sensitivity["rho_nominal"]
    assert sensitivity["rho_nominal"] < sensitivity["rho_plus_0p1pct"]
    assert sensitivity["prefactor_minus_0p1pct"] > 0.0
    assert sensitivity["prefactor_nominal"] > 0.0
    assert sensitivity["prefactor_plus_0p1pct"] > 0.0

    gates = report["immutable_project_gates"]
    assert gates["normalized_momentum_max"] == 1e-3
    assert gates["normalized_momentum_L2"] == 1e-3
    assert gates["divergence_max"] == 1e-5
    assert gates["divergence_L2"] == 1e-5


def test_truth_boundary_keeps_full_r1_and_pde_fail_closed() -> None:
    report = run_independent_audit()
    truth = report["truth_boundary"]
    assert truth["source_R1_averaged_eta_ordinary_slot_independently_audited"] is True
    assert truth["source_R1_zeta_ordinary_slot_independent_audit_pending"] is True
    assert truth["source_R1_prior_five_algebraic_slots_independently_audited"] is False
    assert truth["all_R1_ordinary_slots_independently_audited"] is False
    assert truth["source_operator_constant_M_machine_bound"] is False
    assert truth["source_operator_constant_K_machine_bound"] is False
    assert truth["global_pressure_matched"] is False
    assert truth["global_leading_profile_reconstructed"] is False
    assert truth["leading_only_ns_residual_assessed"] is False
    assert truth["leading_plus_oscillatory_ns_residual_assessed"] is False
    assert truth["after_correction_ns_residual_assessed"] is False
    assert truth["heldout_ns_momentum_residual_assessed"] is False
    assert truth["formal_full_domain_pde_gate_assessed"] is False
    assert truth["pde_validated"] is False


def test_receipt_is_deterministic() -> None:
    first = run_independent_audit()
    second = run_independent_audit()
    assert first == second
    assert len(first["receipt_sha256"]) == 64
