from __future__ import annotations

from openai_ns_reconstruction.kokuno_independent_terminal_tail_audit import run_audit


def test_terminal_tail_independent_audit_passes_frozen_local_guards() -> None:
    report = run_audit()
    summary = report["summary"]

    assert summary["structural_preflight_passed"] is True
    assert all(summary["guard_results"].values())
    assert summary["max_release_relative_error"] <= 5.0e-9
    assert summary["max_terminal_q_relative_error"] <= 5.0e-9
    assert summary["max_log_scale_absolute_error"] <= 5.0e-8
    assert summary["max_independent_medium_to_fine_relative_change"] <= 1.0e-8
    assert summary["minimum_transition_observed_order"] >= 1.8
    assert summary["minimum_mutation_half_q_p_log_X_tail_shift"] >= 0.5
    assert summary["materialization_classification_matches"] is True

    cases = {case["name"]: case for case in report["cases"]}
    assert cases["default_extreme"]["public_materialized"] is False
    assert cases["moderate_materializable"]["public_materialized"] is True
    assert cases["perturbed_outer"]["materialization_classification_match"] is True


def test_terminal_tail_audit_keeps_pde_and_cross_route_truth_fail_closed() -> None:
    report = run_audit()
    truth = report["truth_boundary"]
    baseline = report["cross_route_baseline"]

    assert report["operator_independence"]["training_loss_or_tensor_read"] is False
    assert report["operator_independence"]["production_internal_ode_state_read"] is False
    assert report["operator_independence"]["free_residual_canceling_forcing_used"] is False
    assert report["formal_gates_unchanged"]["normalized_momentum"] == 1.0e-3
    assert report["formal_gates_unchanged"]["divergence_max"] == 1.0e-5
    assert baseline["directly_comparable"] is False
    assert baseline["held_out_momentum_sampled_max"] == 0.1082289305112118
    assert baseline["held_out_momentum_volume_L2"] == 0.10758432876230622

    assert truth["terminal_tail_schedule_independently_preflighted"] is True
    assert truth["global_leading_profile_reconstructed"] is False
    assert truth["complete_kokuno_composite_velocity"] is False
    assert truth["formal_full_domain_pde_gate_assessed"] is False
    assert truth["normalized_ns_residual_le_1e-3_claimed"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False
    assert truth["openai_field_identified"] is False
    assert truth["blowup_proved"] is False
