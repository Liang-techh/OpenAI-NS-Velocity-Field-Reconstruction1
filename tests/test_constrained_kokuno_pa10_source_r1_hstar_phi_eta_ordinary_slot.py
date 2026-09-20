from __future__ import annotations

import math

from openai_ns_reconstruction.kokuno_pa10_post_j_remainder_bridge import (
    R1_ORDINARY_TERMS,
)
from openai_ns_reconstruction.kokuno_pa10_source_r1_hstar_phi_eta_ordinary_slot import (
    R1_HSTAR_PHI_ETA_ORDINARY_TERMS,
    R1_POST_J2_SOURCE_BOUND_THROUGH_HSTAR_PHI_ETA_TERMS,
    R1_REMAINING_AFTER_HSTAR_PHI_ETA_TERMS,
    KokunoPA10SourceR1HstarPhiEtaOrdinarySlot,
)


def test_hstar_phi_eta_increment_is_exactly_one_registered_r1_slot() -> None:
    calc = KokunoPA10SourceR1HstarPhiEtaOrdinarySlot()
    current = calc.r1_hstar_phi_eta_after_j2()
    assert tuple(current) == R1_HSTAR_PHI_ETA_ORDINARY_TERMS
    assert set(current).issubset(set(R1_ORDINARY_TERMS))

    combined = calc.combined_r1_after_j2_subset()
    assert tuple(combined) == R1_POST_J2_SOURCE_BOUND_THROUGH_HSTAR_PHI_ETA_TERMS
    assert len(combined) == 10
    assert R1_REMAINING_AFTER_HSTAR_PHI_ETA_TERMS == (
        "lambda_inv_d_u_times_Phi_eta",
    )
    assert set(combined).isdisjoint(set(R1_REMAINING_AFTER_HSTAR_PHI_ETA_TERMS))
    assert set(combined) | set(R1_REMAINING_AFTER_HSTAR_PHI_ETA_TERMS) == set(
        R1_ORDINARY_TERMS
    )


def test_source_single_eta_post_j2_factor_is_executable_and_frozen() -> None:
    calc = KokunoPA10SourceR1HstarPhiEtaOrdinarySlot()
    factor = calc.source_single_eta_after_j2_factor_upper()
    expected = 80.0 / calc.rho
    assert math.isfinite(factor)
    assert factor >= expected
    assert factor <= math.nextafter(expected, math.inf)


def test_hstar_bound_is_reused_from_existing_source_derivative_route() -> None:
    calc = KokunoPA10SourceR1HstarPhiEtaOrdinarySlot()
    reused = calc.hstar_coefficient_ball()
    upstream = calc.algebraic_r1.r2_source.hstar_coefficient_ball()
    assert reused == upstream
    assert reused.norm > 0.0
    assert reused.lipschitz == 0.0

    inputs = calc.operator_inputs()
    assert tuple(inputs) == ("H_star", "Phi", "L_inverse")
    assert inputs["H_star"] == reused
    assert inputs["Phi"].norm > 1.0 and inputs["Phi"].lipschitz > 0.0
    assert inputs["L_inverse"].norm > 0.0
    assert inputs["L_inverse"].lipschitz == 0.0


def test_hstar_phi_eta_produces_finite_nontrivial_post_j2_bound() -> None:
    calc = KokunoPA10SourceR1HstarPhiEtaOrdinarySlot()
    phi_eta = calc.j2_phi_eta()
    slot = calc.r1_hstar_phi_eta_after_j2()["H_star_times_Phi_eta"]

    assert math.isfinite(phi_eta.norm) and phi_eta.norm > 0.0
    assert math.isfinite(phi_eta.lipschitz) and phi_eta.lipschitz > 0.0
    assert math.isfinite(slot.norm) and slot.norm > 0.0
    assert math.isfinite(slot.lipschitz) and slot.lipschitz > 0.0

    report = calc.report()
    certificate = report["operator_certificate"]
    assert certificate["post_J2_single_eta_factor_upper"] >= 80.0 / calc.rho
    assert certificate["radial_inverse_commutes_with_eta_only_multipliers"] is True
    assert certificate["standalone_Phi_eta_norm_used"] is False
    assert certificate["J2_already_in_source_single_eta_factor"] is True
    assert certificate["ordinary_bridge_J2_reapplied"] is False
    assert certificate["coefficient_product_factor_count_after_J2_rule"] == 3
    assert certificate["Hstar_coefficient_bound_reused_not_duplicated"] is True


def test_combined_post_j2_subset_preserves_prior_nine_slots_exactly() -> None:
    calc = KokunoPA10SourceR1HstarPhiEtaOrdinarySlot()
    prior = calc.parent.combined_r1_after_j2_subset()
    combined = calc.combined_r1_after_j2_subset()
    for name, bound in prior.items():
        assert combined[name] == bound
    assert "H_star_times_Phi_eta" not in prior
    assert "H_star_times_Phi_eta" in combined


def test_report_is_deterministic_and_keeps_only_du_phi_eta_open() -> None:
    calc = KokunoPA10SourceR1HstarPhiEtaOrdinarySlot()
    first = calc.report()
    second = calc.report()
    assert first == second
    assert first["receipt_sha256"] == calc.sha256

    progress = first["ledger_progress"]
    assert progress["R1_post_J2_source_bound_ordinary_slots_count"] == 10
    assert progress["R1_total_ordinary_slots"] == 11
    assert progress["R1_ordinary_slots_still_required"] == 1
    assert tuple(progress["R1_remaining_ordinary_slots"]) == (
        "lambda_inv_d_u_times_Phi_eta",
    )
    assert progress["R1_numerator_space_source_bound_count_remains_six"] is True
    assert progress["R1_mixed_slot_handled_upstream_by_A1_668"] is True
    assert progress["full_R1_still_open"] is True
    assert progress["M_K_not_promoted"] is True

    truth = first["truth_boundary"]
    assert truth["source_Hstar_coefficient_norm_machine_bound_reused"] is True
    assert truth["source_single_eta_post_J2_rule_executable"] is True
    assert truth[
        "source_compatible_R1_Hstar_Phi_eta_ordinary_post_J2_slot_machine_bound"
    ] is True
    assert truth["standalone_Phi_eta_coefficient_norm_invented"] is False
    assert truth["ordinary_bridge_J2_reapplied"] is False
    assert truth["source_compatible_all_R1_ordinary_post_J2_slots_machine_bound"] is False
    assert truth["source_operator_constant_M_machine_bound"] is False
    assert truth["source_operator_constant_K_machine_bound"] is False
    assert truth["global_leading_profile_reconstructed"] is False
    assert truth["heldout_ns_residual_assessed"] is False
    assert truth["pde_validated"] is False
