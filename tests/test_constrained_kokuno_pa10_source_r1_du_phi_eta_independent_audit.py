from __future__ import annotations

from openai_ns_reconstruction.kokuno_pa10_source_r1_du_phi_eta_independent_audit import (
    AGENT1_HEAD,
    AGENT1_PR,
    DIVERGENCE_GATE,
    MOMENTUM_GATE,
    run_independent_audit,
)


def test_independent_du_phi_eta_audit_passes_all_preregistered_guards() -> None:
    report = run_independent_audit()
    assert report["audited_agent1_pr"] == AGENT1_PR == 751
    assert report["audited_agent1_exact_head"] == AGENT1_HEAD
    assert report["passed"] is True
    assert report["failed_guards"] == []


def test_public_bounds_are_tight_upper_bounds_of_independent_reconstruction() -> None:
    report = run_independent_audit()
    ratios = report["public_to_independent_ratios"]
    assert set(ratios) == {
        "product_constant",
        "L_inverse",
        "d",
        "single_eta_after_J2",
        "J2_u_Phi_eta_norm",
        "J2_u_Phi_eta_lipschitz",
        "slot_norm",
        "slot_lipschitz",
    }
    assert all(1.0 <= value <= 1.01 for value in ratios.values())


def test_fresh_offgrid_d_path_and_axis_near_probes_are_nonvacuous() -> None:
    audit = run_independent_audit()["fresh_offgrid_d_value_audit"]
    assert audit["fresh_random_offgrid_points"] == 8192
    assert audit["axis_near_probes"] == [-1.0e-12, 0.0, 1.0e-12]
    assert audit["finite"] is True
    assert audit["value_path_passed"] is True
    assert audit["dominated"] is True
    assert audit["nontrivial"] is True
    assert audit["max_abs_d"] >= 0.999
    assert audit["wrong_sign_d_mutation_detected"] is True
    assert audit["wrong_sign_d_mutation_response"] > 0.0


def test_all_frozen_negative_controls_are_detected() -> None:
    controls = run_independent_audit()["negative_controls"]
    assert controls
    assert all(controls.values())


def test_truth_boundary_keeps_full_r1_and_pde_fail_closed() -> None:
    report = run_independent_audit()
    truth = report["truth_boundary"]
    assert truth["source_R1_du_Phi_eta_ordinary_slot_independently_audited"] is True
    assert truth["source_R1_Hstar_Phi_eta_independent_audit_pending"] is True
    assert truth["source_R1_Wstar_logradial_independent_audit_pending"] is True
    assert truth["source_R1_averaged_logradial_independent_audit_pending"] is True
    assert truth["source_R1_all_ordinary_slots_independently_admitted"] is False
    assert truth["source_full_post_J2_R1_independently_admitted"] is False
    assert truth["source_operator_constant_M_machine_bound"] is False
    assert truth["source_operator_constant_K_machine_bound"] is False
    assert truth["global_pressure_matched"] is False
    assert truth["global_leading_profile_reconstructed"] is False
    assert truth["heldout_ns_momentum_residual_assessed"] is False
    assert truth["pde_validated"] is False

    gates = report["project_gates_unchanged"]
    assert gates["normalized_momentum_max"] == MOMENTUM_GATE == 1.0e-3
    assert gates["normalized_momentum_L2"] == MOMENTUM_GATE
    assert gates["divergence_max"] == DIVERGENCE_GATE == 1.0e-5
    assert gates["divergence_L2"] == DIVERGENCE_GATE
    assert gates["free_residual_defined_forcing_allowed"] is False


def test_receipt_is_deterministic_for_fixed_head_binding() -> None:
    first = run_independent_audit(pr_head="head", checkout_head="head")
    second = run_independent_audit(pr_head="head", checkout_head="head")
    assert first == second
    assert first["receipt_sha256"] == second["receipt_sha256"]
    assert first["pr_head"] == "head"
    assert first["checkout_head"] == "head"
