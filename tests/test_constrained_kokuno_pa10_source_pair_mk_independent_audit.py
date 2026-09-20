from __future__ import annotations

from openai_ns_reconstruction.kokuno_pa10_source_pair_mk_independent_audit import (
    AGENT1_HEAD,
    AGENT1_PR,
    DIVERGENCE_GATE,
    MOMENTUM_GATE,
    run_independent_audit,
)


def test_pair_mk_independent_audit_passes_preregistered_guards() -> None:
    report = run_independent_audit()
    assert report["audited_agent1_pr"] == AGENT1_PR == 760
    assert report["audited_agent1_exact_head"] == AGENT1_HEAD
    assert report["passed"] is True
    assert report["failed_guards"] == []


def test_pair_mk_public_values_are_tight_uppers_of_decimal_rebuild() -> None:
    report = run_independent_audit()
    ratios = report["public_to_independent_ratios"]
    assert set(ratios) == {
        "inverse_one_plus_T_absolute_series_upper",
        "M_phi_component",
        "M_u_component",
        "M_pair_max",
        "K_phi_component",
        "K_u_component",
        "K_pair_max",
    }
    assert all(1.0 <= value <= 1.000000000001 for value in ratios.values())


def test_no_double_j_and_half_factor_mutations_are_detected() -> None:
    report = run_independent_audit()
    controls = report["negative_controls"]
    assert controls
    assert all(controls.values())
    assert controls["duplicate_J2_on_already_post_J2_R1_detected"] is True
    assert controls["duplicate_J1_on_already_post_J1_R2_detected"] is True
    assert controls["missing_pair_factor_one_half_detected"] is True


def test_component_max_and_inverse_sensitivity_shape_are_correct() -> None:
    report = run_independent_audit()
    guards = report["component_guards"]
    assert guards["M_pair_is_component_max"] is True
    assert guards["K_pair_is_component_max"] is True
    assert guards["all_components_nontrivial"] is True

    sensitivity = report["inverse_plus_minus_0p1_percent_sensitivity"]
    assert sensitivity["perturbation_fraction"] == 0.001
    assert sensitivity["M_phi_low_over_center"] == 0.999
    assert sensitivity["M_phi_high_over_center"] == 1.001
    assert sensitivity["K_phi_low_over_center"] == 0.999
    assert sensitivity["K_phi_high_over_center"] == 1.001
    assert sensitivity["M_u_invariant"] is True
    assert sensitivity["K_u_invariant"] is True
    assert sensitivity["monotone_phi_components"] is True


def test_truth_boundary_keeps_upstream_admission_and_pde_fail_closed() -> None:
    report = run_independent_audit()
    truth = report["truth_boundary"]
    assert truth["source_pair_MK_arithmetic_independently_audited"] is True
    assert truth["source_pair_MK_audit_is_conditional_on_post_J_inputs"] is True
    assert truth["source_full_post_J2_R1_independent_admission"] is False
    assert truth["source_full_post_J1_R2_independent_admission"] is False
    assert truth["source_operator_constant_M_independent_agent4_admission"] is False
    assert truth["source_operator_constant_K_independent_agent4_admission"] is False
    assert truth["global_pressure_matched"] is False
    assert truth["global_leading_profile_reconstructed"] is False
    assert truth["leading_only_ns_residual_assessed"] is False
    assert truth["leading_plus_oscillatory_ns_residual_assessed"] is False
    assert truth["after_correction_ns_residual_assessed"] is False
    assert truth["heldout_ns_momentum_residual_assessed"] is False
    assert truth["pde_validated"] is False

    gates = report["project_gates_unchanged"]
    assert gates["normalized_momentum_max"] == MOMENTUM_GATE == 1.0e-3
    assert gates["normalized_momentum_L2"] == MOMENTUM_GATE
    assert gates["divergence_max"] == DIVERGENCE_GATE == 1.0e-5
    assert gates["divergence_L2"] == DIVERGENCE_GATE
    assert gates["free_residual_defined_forcing_allowed"] is False

    comparison = report["same_protocol_PDE_comparison"]
    assert comparison["performed"] is False


def test_receipt_is_deterministic_for_fixed_head_binding() -> None:
    first = run_independent_audit(pr_head="head", checkout_head="head")
    second = run_independent_audit(pr_head="head", checkout_head="head")
    assert first == second
    assert first["receipt_sha256"] == second["receipt_sha256"]
    assert first["pr_head"] == "head"
    assert first["checkout_head"] == "head"
