from __future__ import annotations

import json

from openai_ns_reconstruction.kokuno_agent4_pa10_source_u_ball_audit import run_audit


def test_independent_source_u_ball_audit_passes_without_promoting_pde() -> None:
    report = run_audit(pr_head="test-head", checkout_head="test-head")
    diagnostic = json.dumps(
        {
            "failed_guards": report["failed_guards"],
            "guards": report["guards"],
            "mutation": report["mutation"],
            "ratios": report["public_to_independent_ratios"],
            "public_certificate": report["public_certificate"],
            "fresh_real_value_path": report["fresh_real_value_path"],
            "fresh_complex_offgrid_stress": report["fresh_complex_offgrid_stress"],
        },
        sort_keys=True,
    )
    assert report["source_u_ball_independent_preflight_passed"] is True, diagnostic
    assert report["failed_guards"] == [], diagnostic
    assert all(report["guards"].values()), diagnostic
    assert all(report["mutation"].values()), diagnostic

    ratios = report["public_to_independent_ratios"]
    assert ratios["slope"] >= 1.0
    assert ratios["weight_sum"] >= 1.0
    assert ratios["u0_norm"] >= 1.0
    assert ratios["radius_one_ball_norm"] >= 1.0

    stress = report["fresh_complex_offgrid_stress"]
    assert stress["max_to_independent_majorant_ratio"] < 1.0
    assert stress["min_abs_L"] > 0.0

    truth = report["truth_boundary"]
    assert truth["source_axis_domain_independently_audited_by_A4_671"] is True
    assert truth["source_u_radius_one_ball_norm_independently_audited"] is True
    assert truth["agent1_684_algebraic_R2_slots_independently_audited"] is False
    assert truth["global_pressure_matched"] is False
    assert truth["global_leading_profile_reconstructed"] is False
    assert truth["heldout_ns_momentum_residual_assessed"] is False
    assert truth["pde_validated"] is False

    gates = report["frozen_protocol"]["final_project_gates_unchanged"]
    assert gates == {
        "normalized_momentum_max": 1.0e-3,
        "normalized_momentum_L2": 1.0e-3,
        "divergence_max": 1.0e-5,
        "divergence_L2": 1.0e-5,
    }
