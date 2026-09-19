from __future__ import annotations

from openai_ns_reconstruction.kokuno_pa10_source_r1_zeta_independent_audit import (
    AGENT1_HEAD,
    AGENT1_PR,
    run_independent_audit,
)


def test_independent_r1_zeta_audit_passes_all_preregistered_guards() -> None:
    report = run_independent_audit()
    assert report["agent1_pr"] == AGENT1_PR == 707
    assert report["agent1_exact_head"] == AGENT1_HEAD
    assert report["failed_guards"] == []
    assert report["source_R1_zeta_ordinary_slot_independent_preflight_passed"] is True

    path = report["independent_path"]
    assert path["private_agent1_helpers_called"] is False
    assert path["agent1_receipt_used_as_oracle"] is False
    assert path["training_internal_tensors_read"] is False

    ratios = report["public_over_independent_ratios"]
    assert min(ratios["fixed"].values()) >= 1.0
    assert min(ratios["R1_zeta_slot"].values()) >= 1.0


def test_fresh_real_and_complex_offgrid_checks_are_nontrivial() -> None:
    report = run_independent_audit()
    real = report["fresh_real_offgrid_replay"]
    assert real["finite"] is True
    assert real["relative_max_mismatch"] <= 1e-12
    assert real["zeta_rms"] > 1e-8
    assert real["sigma_pm_0p1_percent_max_relative_response"] > 1e-8

    complex_stress = report["fresh_complex_cauchy_stress"]
    assert complex_stress["finite"] is True
    assert complex_stress["dominated"] is True
    assert 0.0 < complex_stress["stress_ratio"] <= 1.0


def test_negative_controls_are_detected() -> None:
    report = run_independent_audit()
    negative = report["negative_controls"]
    assert all(negative["fixed_0p999_detected"].values())
    assert all(negative["slot_0p99_detected"].values())
    assert negative["real_value_sign_flip_detected"] is True


def test_truth_boundary_keeps_full_r1_r2_and_pde_fail_closed() -> None:
    report = run_independent_audit()
    truth = report["truth_boundary"]
    assert truth["source_R1_zeta_ordinary_slot_independently_audited"] is True
    assert truth["source_R1_other_ordinary_slots_independently_audited"] is False
    assert truth["full_R1_independent_agent4_admission"] is False
    assert truth["full_R2_independent_agent4_admission"] is False
    assert truth["global_pressure_matched"] is False
    assert truth["global_leading_profile_reconstructed"] is False
    assert truth["heldout_ns_momentum_residual_assessed"] is False
    assert truth["pde_validated"] is False

    gates = report["immutable_project_gates"]
    assert gates["normalized_momentum_max"] == 1e-3
    assert gates["normalized_momentum_L2"] == 1e-3
    assert gates["divergence_max"] == 1e-5
    assert gates["divergence_L2"] == 1e-5
