from __future__ import annotations

import math

from openai_ns_reconstruction.kokuno_pa10_source_derivative_slots_independent_audit import (
    AGENT1_HEAD,
    AGENT1_PR,
    run_independent_audit,
)


def test_independent_derivative_slot_audit_passes_all_preregistered_guards() -> None:
    report = run_independent_audit()
    assert report["agent1_pr"] == AGENT1_PR == 693
    assert report["agent1_exact_head"] == AGENT1_HEAD
    assert report["failed_guards"] == []
    assert report["source_R2_derivative_ordinary_slots_independent_preflight_passed"] is True

    path = report["independent_path"]
    assert path["private_agent1_helpers_called"] is False
    assert path["agent1_receipt_used_as_oracle"] is False
    assert path["training_internal_tensors_read"] is False

    ratios = report["public_over_independent_ratios"]
    assert min(ratios["fixed"].values()) >= 1.0
    assert len(ratios["derivative_slots"]) == 4
    for slot in ratios["derivative_slots"].values():
        assert slot["norm"] >= 1.0
        assert slot["lipschitz"] >= 1.0


def test_negative_controls_are_nonvacuous_and_detected() -> None:
    report = run_independent_audit()
    negative = report["negative_controls"]
    assert all(negative["fixed_0p999_detected"].values())
    assert all(negative["slot_0p99_detected"].values())

    stress = report["fresh_offgrid_polynomial_stress"]
    assert stress["finite"] is True
    assert stress["dominated"] is True
    assert stress["nontrivial"] is True
    assert 0.0 < stress["W_star_stress_ratio"] <= 1.0
    assert 0.0 < stress["H_star_stress_ratio"] <= 1.0


def test_weight_ratio_stress_independently_recovers_safe_logradial_envelope() -> None:
    report = run_independent_audit()
    stress = report["weight_ratio_stress"]
    assert stress["envelope_dominates"] is True
    assert stress["shared_source_envelope"] == 80.0
    assert math.isclose(
        stress["max_single_logradial_ratio"], 11.25, rel_tol=0.0, abs_tol=1e-12
    )


def test_truth_boundary_keeps_full_r2_and_pde_fail_closed() -> None:
    report = run_independent_audit()
    truth = report["truth_boundary"]
    assert truth["source_R2_derivative_ordinary_slots_independently_audited"] is True
    assert truth["source_R2_algebraic_ordinary_slots_independent_audit_pending"] is True
    assert truth["all_R2_ordinary_slots_independently_audited"] is False
    assert truth["all_R1_ordinary_slots_source_bound"] is False
    assert truth["global_pressure_matched"] is False
    assert truth["global_leading_profile_reconstructed"] is False
    assert truth["heldout_ns_momentum_residual_assessed"] is False
    assert truth["pde_validated"] is False

    gates = report["immutable_project_gates"]
    assert gates["normalized_momentum_max"] == 1e-3
    assert gates["normalized_momentum_L2"] == 1e-3
    assert gates["divergence_max"] == 1e-5
    assert gates["divergence_L2"] == 1e-5
