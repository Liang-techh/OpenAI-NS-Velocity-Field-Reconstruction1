from __future__ import annotations

import math

from openai_ns_reconstruction.kokuno_pa10_post_j_remainder_bridge import (
    R1_ORDINARY_TERMS,
)
from openai_ns_reconstruction.kokuno_pa10_source_r1_wstar_logradial_ordinary_slot import (
    R1_POST_J2_SOURCE_BOUND_THROUGH_WSTAR_LOGRADIAL_TERMS,
    R1_REMAINING_AFTER_WSTAR_LOGRADIAL_TERMS,
    R1_WSTAR_LOGRADIAL_ORDINARY_TERMS,
    KokunoPA10SourceR1WstarLogradialOrdinarySlot,
)


def test_wstar_logradial_increment_is_exactly_one_registered_r1_slot() -> None:
    calc = KokunoPA10SourceR1WstarLogradialOrdinarySlot()
    current = calc.r1_wstar_logradial_after_j2()
    assert tuple(current) == R1_WSTAR_LOGRADIAL_ORDINARY_TERMS
    assert set(current).issubset(set(R1_ORDINARY_TERMS))

    combined = calc.combined_r1_after_j2_subset()
    assert tuple(combined) == R1_POST_J2_SOURCE_BOUND_THROUGH_WSTAR_LOGRADIAL_TERMS
    assert len(combined) == 8
    assert len(R1_REMAINING_AFTER_WSTAR_LOGRADIAL_TERMS) == 3
    assert set(combined).isdisjoint(set(R1_REMAINING_AFTER_WSTAR_LOGRADIAL_TERMS))
    assert set(combined) | set(R1_REMAINING_AFTER_WSTAR_LOGRADIAL_TERMS) == set(
        R1_ORDINARY_TERMS
    )


def test_source_single_logradial_post_j2_factor_is_executable() -> None:
    calc = KokunoPA10SourceR1WstarLogradialOrdinarySlot()
    observed = calc.source_single_logradial_after_j2_factor_upper()
    assert math.isfinite(observed)
    assert observed >= 80.0
    assert observed <= math.nextafter(80.0, math.inf)


def test_operator_inputs_are_fixed_eta_multipliers_and_phi_ball() -> None:
    calc = KokunoPA10SourceR1WstarLogradialOrdinarySlot()
    inputs = calc.operator_inputs()
    assert tuple(inputs) == ("W_star", "Phi", "L_inverse")
    assert inputs["W_star"].norm > 0.0
    assert inputs["W_star"].lipschitz == 0.0
    assert inputs["L_inverse"].norm > 0.0
    assert inputs["L_inverse"].lipschitz == 0.0
    assert inputs["Phi"].norm > 1.0
    assert inputs["Phi"].lipschitz > 0.0


def test_new_post_j2_slot_is_finite_nontrivial_and_not_double_integrated() -> None:
    calc = KokunoPA10SourceR1WstarLogradialOrdinarySlot()
    slot = calc.r1_wstar_logradial_after_j2()["W_star_times_Y_Phi_Y"]
    assert math.isfinite(slot.norm) and slot.norm > 0.0
    assert math.isfinite(slot.lipschitz) and slot.lipschitz > 0.0

    report = calc.report()
    certificate = report["operator_certificate"]
    assert certificate["post_J2_single_logradial_factor_upper"] >= 80.0
    assert certificate["eta_only_multipliers_commuted_through_Y_partial_Y"] is True
    assert certificate["standalone_Y_Phi_Y_norm_used"] is False
    assert certificate["J2_already_in_source_logradial_factor"] is True
    assert certificate["ordinary_bridge_J2_reapplied"] is False
    assert certificate["coefficient_product_factor_count_before_post_J2_rule"] == 3


def test_combined_post_j2_subset_preserves_prior_seven_slots_exactly() -> None:
    calc = KokunoPA10SourceR1WstarLogradialOrdinarySlot()
    prior = calc.parent.combined_r1_after_j2_subset()
    combined = calc.combined_r1_after_j2_subset()
    for name, bound in prior.items():
        assert combined[name] == bound
    assert "W_star_times_Y_Phi_Y" not in prior
    assert "W_star_times_Y_Phi_Y" in combined


def test_report_is_deterministic_and_keeps_three_r1_slots_open() -> None:
    calc = KokunoPA10SourceR1WstarLogradialOrdinarySlot()
    first = calc.report()
    second = calc.report()
    assert first == second
    assert first["receipt_sha256"] == calc.sha256

    progress = first["ledger_progress"]
    assert progress["R1_post_J2_source_bound_ordinary_slots_count"] == 8
    assert progress["R1_total_ordinary_slots"] == 11
    assert progress["R1_ordinary_slots_still_required"] == 3
    assert tuple(progress["R1_remaining_ordinary_slots"]) == (
        "lambda_inv_2D_eta_Au_times_Y_Phi_Y",
        "H_star_times_Phi_eta",
        "lambda_inv_d_u_times_Phi_eta",
    )
    assert progress["R1_numerator_space_source_bound_count_remains_six"] is True
    assert progress["R1_mixed_slot_handled_upstream_by_A1_668"] is True
    assert progress["full_R1_still_open"] is True
    assert progress["M_K_not_promoted"] is True

    truth = first["truth_boundary"]
    assert truth["source_single_logradial_post_J2_rule_executable"] is True
    assert truth[
        "source_compatible_R1_Wstar_logradial_ordinary_post_J2_slot_machine_bound"
    ] is True
    assert truth["standalone_Y_Phi_Y_coefficient_norm_invented"] is False
    assert truth["standalone_Phi_eta_coefficient_norm_invented"] is False
    assert truth["ordinary_bridge_J2_reapplied"] is False
    assert truth["source_compatible_all_R1_ordinary_post_J2_slots_machine_bound"] is False
    assert truth["source_operator_constant_M_machine_bound"] is False
    assert truth["source_operator_constant_K_machine_bound"] is False
    assert truth["global_leading_profile_reconstructed"] is False
    assert truth["heldout_ns_residual_assessed"] is False
    assert truth["pde_validated"] is False
