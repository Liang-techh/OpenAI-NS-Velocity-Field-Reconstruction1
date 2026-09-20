from __future__ import annotations

from openai_ns_reconstruction.kokuno_pa10_source_r1_hstar_phi_eta_independent_audit import (
    AGENT1_HEAD,
    DIVERGENCE_GATE,
    MAX_PUBLIC_RATIO,
    MOMENTUM_GATE,
    run_independent_audit,
)
from openai_ns_reconstruction.kokuno_pa10_source_r1_hstar_phi_eta_ordinary_slot import (
    KokunoPA10SourceR1HstarPhiEtaOrdinarySlot,
)


def test_hstar_phi_eta_independent_audit_passes_frozen_guards() -> None:
    payload = run_independent_audit()

    assert payload["audited_agent1_exact_head"] == AGENT1_HEAD
    assert payload["passed"] is True
    assert payload["failed_guards"] == []

    for ratio in payload["public_to_independent_ratios"].values():
        assert 1.0 <= ratio <= float(MAX_PUBLIC_RATIO)

    stress = payload["fresh_offgrid_Hstar_stress"]
    assert stress["fresh_random_offgrid_points"] == 8192
    assert stress["finite"] is True
    assert stress["dominated"] is True
    assert stress["nontrivial"] is True
    assert stress["j0_pm_0p1_percent_response"] > 0.0
    assert stress["quadratic_sign_mutation_detected"] is True

    assert payload["negative_controls"]
    assert all(payload["negative_controls"].values())


def test_hstar_phi_eta_audit_preserves_agent4_truth_boundary() -> None:
    payload = run_independent_audit()
    truth = payload["truth_boundary"]

    assert truth["source_R1_Hstar_Phi_eta_ordinary_slot_independently_audited"] is True
    assert truth["source_R1_all_ordinary_slots_independently_admitted"] is False
    assert truth["source_full_post_J2_R1_independently_admitted"] is False
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

    gates = payload["project_gates_unchanged"]
    assert gates["normalized_momentum_max"] == MOMENTUM_GATE == 1.0e-3
    assert gates["normalized_momentum_L2"] == MOMENTUM_GATE
    assert gates["divergence_max"] == DIVERGENCE_GATE == 1.0e-5
    assert gates["divergence_L2"] == DIVERGENCE_GATE
    assert gates["free_residual_defined_forcing_allowed"] is False


def test_audited_agent1_contract_does_not_invent_phi_eta_norm_or_double_j2() -> None:
    calc = KokunoPA10SourceR1HstarPhiEtaOrdinarySlot()
    report = calc.report()
    certificate = report["operator_certificate"]
    progress = report["ledger_progress"]

    assert certificate["standalone_Phi_eta_norm_used"] is False
    assert certificate["J2_already_in_source_single_eta_factor"] is True
    assert certificate["ordinary_bridge_J2_reapplied"] is False
    assert certificate["coefficient_product_factor_count_after_J2_rule"] == 3
    assert progress["R1_post_J2_source_bound_ordinary_slots_count"] == 10
    assert progress["R1_total_ordinary_slots"] == 11
    assert progress["R1_ordinary_slots_still_required"] == 1
    assert progress["R1_remaining_ordinary_slots"] == [
        "lambda_inv_d_u_times_Phi_eta"
    ]
    assert progress["full_R1_still_open"] is True
