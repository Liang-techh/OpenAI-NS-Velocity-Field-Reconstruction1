from __future__ import annotations

from openai_ns_reconstruction.kokuno_pa10_source_r1_averaged_logradial_independent_audit import (
    OFFGRID_COUNT,
    _manufactured_averaging_identity_stress,
    run_independent_audit,
)


def test_averaged_logradial_independent_audit_passes_frozen_guards() -> None:
    payload = run_independent_audit()
    assert payload["failed_guards"] == []
    assert (
        payload[
            "source_R1_averaged_logradial_ordinary_slot_independent_preflight_passed"
        ]
        is True
    )

    ratios = payload["public_over_independent_ratios"]
    assert min(ratios.values()) >= 1.0
    assert max(ratios.values()) <= 1.01
    assert all(payload["negative_controls"].values())

    offgrid = payload["fresh_offgrid_eta_stress"]
    assert offgrid["fresh_random_offgrid_points"] == OFFGRID_COUNT
    assert offgrid["finite"] is True
    assert offgrid["dominated"] is True
    assert offgrid["axis_nontrivial"] is True
    assert offgrid["h_pm_0p1_percent_2D_response"] > 0.0

    identity = payload["manufactured_averaging_identity_stress"]
    assert identity["identity_passed"] is True
    assert identity["wrong_sign_mutation_detected"] is True
    assert identity["nontrivial"] is True

    gates = payload["immutable_project_gates"]
    assert gates["normalized_momentum_max"] == 1e-3
    assert gates["normalized_momentum_L2"] == 1e-3
    assert gates["divergence_max"] == 1e-5
    assert gates["divergence_L2"] == 1e-5
    assert gates["free_residual_defined_forcing_allowed"] is False


def test_manufactured_averaging_identity_detects_sign_mutation() -> None:
    stress = _manufactured_averaging_identity_stress()
    assert stress["cases"] == 512
    assert stress["max_abs_identity_error"] <= stress["identity_tolerance"]
    assert stress["wrong_sign_mutation_max_abs_error"] >= 1e-6
    assert stress["wrong_sign_mutation_detected"] is True


def test_truth_boundary_remains_fail_closed() -> None:
    truth = run_independent_audit()["truth_boundary"]
    assert truth["source_R1_averaged_logradial_ordinary_slot_independently_audited"] is True
    assert truth["source_R1_Wstar_logradial_ordinary_slot_independent_audit_pending"] is True
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
