from __future__ import annotations

from openai_ns_reconstruction.kokuno_agent4_source_signed_covariance_mass_independent_audit import (
    FROZEN_GUARDS,
    H_CASES,
    LADDER,
    run_audit,
)


def test_independent_signed_covariance_mass_audit_passes_frozen_guards() -> None:
    report = run_audit()
    metrics = report["metrics"]

    assert report["local_structural_preflight_passed"] is True
    assert tuple(report["h_cases"]) == H_CASES
    assert tuple(report["resolution_ladder_samples"]) == LADDER
    assert metrics["pulse_relative_rms_ladder"][-1] <= FROZEN_GUARDS[
        "pulse_finest_relative_rms"
    ]
    assert min(metrics["pulse_refinement_ratios"]) >= FROZEN_GUARDS[
        "pulse_refinement_ratio_min"
    ]
    assert metrics["transverse_relative_ladder"][-1] <= FROZEN_GUARDS[
        "transverse_finest_relative"
    ]
    assert min(metrics["transverse_refinement_ratios"]) >= FROZEN_GUARDS[
        "transverse_refinement_ratio_min"
    ]
    assert metrics["combined_h_sigma_relative_rms_ladder"][-1] <= FROZEN_GUARDS[
        "combined_finest_relative_rms"
    ]
    assert metrics["combined_h_sigma_relative_max_ladder"][-1] <= FROZEN_GUARDS[
        "combined_finest_relative_max"
    ]


def test_independent_audit_mutations_remain_detectable() -> None:
    report = run_audit()
    metrics = report["metrics"]

    assert metrics["wrong_psi_power_mutation_relative_rms"] >= FROZEN_GUARDS[
        "wrong_psi_power_mutation_relative_rms_min"
    ]
    assert metrics["wrong_transverse_prefactor_mutation_relative_rms"] >= FROZEN_GUARDS[
        "wrong_transverse_prefactor_mutation_relative_min"
    ]


def test_independent_audit_preserves_pde_truth_boundary() -> None:
    report = run_audit()
    oracle = report["oracle"]
    truth = report["truth_boundary"]

    assert oracle["agent2_trapezoid_helper_reused"] is False
    assert oracle["agent2_band_covering_helper_reused"] is False
    assert oracle["agent2_signed_inverse_reused"] is False
    assert oracle["agent2_complete_curl_reused"] is False
    assert oracle["training_tensor_or_loss_read"] is False
    assert oracle["pressure_or_forcing_fit_used"] is False
    assert oracle["pulse_and_transverse_resolution_varied_separately"] is True

    assert report["formal_full_domain_pde_gate_assessed"] is False
    assert report["heldout_ns_residual_assessed"] is False
    assert report["residual_reduction_claimed"] is False
    assert report["pde_validated"] is False
    assert truth["candidate_quadrature_only"] is True
    assert truth["actual_source_pulse_samples_recovered"] is False
    assert truth["actual_source_h_sigma_bound"] is False
    assert truth["actual_positive_order_background_bound"] is False
    assert truth["actual_auxiliary_torus_mode_family_bound"] is False
    assert truth["public_source_bound_velocity_osc_materialized"] is False
    assert truth["public_velocity_correction_materialized"] is False
    assert truth["full_ns_momentum_gate_value"] is None
    assert truth["full_ns_divergence_gate_value"] is None
