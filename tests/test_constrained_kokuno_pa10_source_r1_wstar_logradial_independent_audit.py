from __future__ import annotations

import pytest

from openai_ns_reconstruction.kokuno_pa10_source_r1_wstar_logradial_independent_audit import (
    AGENT1_HEAD,
    AGENT1_PR,
    OFFGRID_COUNT,
    SEED,
    run_independent_audit,
)


def test_independent_wstar_logradial_audit_passes_frozen_guards() -> None:
    report = run_independent_audit()
    assert report["agent1_pr"] == AGENT1_PR == 724
    assert report["agent1_exact_head"] == AGENT1_HEAD
    assert report["failed_guards"] == []
    assert report[
        "source_R1_Wstar_logradial_ordinary_slot_independent_preflight_passed"
    ] is True

    ratios = report["public_over_independent_ratios"]
    assert min(ratios.values()) >= 1.0
    assert max(ratios.values()) <= 1.01
    assert ratios["single_logradial_after_J2"] == pytest.approx(1.0, abs=1e-15)


def test_fresh_offgrid_and_weight_ratio_stresses_are_nonvacuous() -> None:
    report = run_independent_audit()
    stress = report["fresh_offgrid_Wstar_stress"]
    assert stress["seed"] == SEED
    assert stress["fresh_random_offgrid_points"] == OFFGRID_COUNT == 8192
    assert stress["finite"] is True
    assert stress["dominated"] is True
    assert stress["nontrivial"] is True
    assert 0.0 < stress["value_over_public_bound"] <= 1.0
    assert stress["j0_pm_0p1_percent_max_abs_response"] > 1e-12
    assert stress["axis_near_probes"] == [-1e-12, 0.0, 1e-12]

    weights = report["single_logradial_weight_ratio_stress"]
    assert weights["alpha_max"] == 512
    assert weights["beta_max"] == 128
    assert weights["max_single_logradial_ratio"] > 0.0
    assert weights["source_envelope"] == 80.0
    assert weights["envelope_dominates"] is True


def test_preregistered_mutations_are_detected() -> None:
    report = run_independent_audit()
    neg = report["negative_controls"]
    assert neg
    assert all(neg.values())


def test_truth_boundary_keeps_full_pde_gate_closed() -> None:
    report = run_independent_audit()
    truth = report["truth_boundary"]
    assert truth["source_R1_Wstar_logradial_ordinary_slot_independently_audited"] is True
    assert truth["source_R1_zeta_ordinary_slot_independently_audited"] is True
    assert truth["source_R1_averaged_eta_ordinary_slot_independent_audit_pending"] is True
    assert truth["source_R1_prior_five_algebraic_slots_independently_audited"] is False
    assert truth["all_R1_ordinary_slots_independently_audited"] is False
    assert truth["source_operator_constant_M_machine_bound"] is False
    assert truth["source_operator_constant_K_machine_bound"] is False
    assert truth["global_pressure_matched"] is False
    assert truth["global_leading_profile_reconstructed"] is False
    assert truth["heldout_ns_momentum_residual_assessed"] is False
    assert truth["pde_validated"] is False

    gates = report["immutable_project_gates"]
    assert gates == {
        "normalized_momentum_max": 1e-3,
        "normalized_momentum_L2": 1e-3,
        "divergence_max": 1e-5,
        "divergence_L2": 1e-5,
    }
