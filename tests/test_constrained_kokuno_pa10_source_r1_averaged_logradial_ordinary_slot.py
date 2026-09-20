from __future__ import annotations

import math

from openai_ns_reconstruction.kokuno_pa10_post_j_remainder_bridge import (
    R1_ORDINARY_TERMS,
)
from openai_ns_reconstruction.kokuno_pa10_source_r1_averaged_logradial_ordinary_slot import (
    R1_AVERAGED_LOGRADIAL_ORDINARY_TERMS,
    R1_POST_J2_SOURCE_BOUND_THROUGH_AVERAGED_LOGRADIAL_TERMS,
    R1_REMAINING_AFTER_AVERAGED_LOGRADIAL_TERMS,
    KokunoPA10SourceR1AveragedLogradialOrdinarySlot,
)


def test_averaged_logradial_increment_is_exactly_one_registered_r1_slot() -> None:
    calc = KokunoPA10SourceR1AveragedLogradialOrdinarySlot()
    current = calc.r1_averaged_logradial_after_j2()
    assert tuple(current) == R1_AVERAGED_LOGRADIAL_ORDINARY_TERMS
    assert set(current).issubset(set(R1_ORDINARY_TERMS))

    combined = calc.combined_r1_after_j2_subset()
    assert tuple(combined) == R1_POST_J2_SOURCE_BOUND_THROUGH_AVERAGED_LOGRADIAL_TERMS
    assert len(combined) == 9
    assert len(R1_REMAINING_AFTER_AVERAGED_LOGRADIAL_TERMS) == 2
    assert set(combined).isdisjoint(set(R1_REMAINING_AFTER_AVERAGED_LOGRADIAL_TERMS))
    assert set(combined) | set(R1_REMAINING_AFTER_AVERAGED_LOGRADIAL_TERMS) == set(
        R1_ORDINARY_TERMS
    )


def test_source_post_j2_factors_are_executable_and_frozen() -> None:
    calc = KokunoPA10SourceR1AveragedLogradialOrdinarySlot()
    logradial = calc.source_single_logradial_after_j2_factor_upper()
    plain = calc.source_plain_after_j2_factor_upper()
    assert math.isfinite(logradial) and math.isfinite(plain)
    assert 80.0 <= logradial <= math.nextafter(80.0, math.inf)
    assert 40.0 <= plain <= math.nextafter(40.0, math.inf)


def test_operator_inputs_use_existing_source_balls_and_fixed_multipliers() -> None:
    calc = KokunoPA10SourceR1AveragedLogradialOrdinarySlot()
    inputs = calc.operator_inputs()
    assert tuple(inputs) == ("eta", "A_u", "u", "Phi", "L_inverse")
    assert inputs["eta"].norm > 0.0 and inputs["eta"].lipschitz == 0.0
    assert inputs["L_inverse"].norm > 0.0
    assert inputs["L_inverse"].lipschitz == 0.0
    assert inputs["A_u"].norm == inputs["u"].norm
    assert inputs["A_u"].lipschitz == inputs["u"].lipschitz
    assert inputs["u"].norm > 1.0 and inputs["u"].lipschitz > 0.0
    assert inputs["Phi"].norm > 1.0 and inputs["Phi"].lipschitz > 0.0


def test_averaging_decomposition_produces_finite_nontrivial_post_j2_bound() -> None:
    calc = KokunoPA10SourceR1AveragedLogradialOrdinarySlot()
    core = calc.j2_average_times_logradial_phi()
    slot = calc.r1_averaged_logradial_after_j2()[
        "lambda_inv_2D_eta_Au_times_Y_Phi_Y"
    ]
    assert math.isfinite(core.norm) and core.norm > 0.0
    assert math.isfinite(core.lipschitz) and core.lipschitz > 0.0
    assert math.isfinite(slot.norm) and slot.norm > 0.0
    assert math.isfinite(slot.lipschitz) and slot.lipschitz > 0.0

    report = calc.report()
    certificate = report["operator_certificate"]
    assert certificate["post_J2_single_logradial_factor_upper"] >= 80.0
    assert certificate["post_J2_plain_factor_upper"] >= 40.0
    assert certificate["averaging_identity_used"] == (
        "A(u)Y Phi_Y=Y partial_Y(A(u)Phi)-u Phi+A(u)Phi"
    )
    assert certificate["eta_and_L_inverse_commute_with_radial_operator"] is True
    assert certificate["standalone_Y_Phi_Y_norm_used"] is False
    assert certificate["J2_already_applied_in_decomposition_terms"] is True
    assert certificate["ordinary_bridge_J2_reapplied"] is False
    assert certificate["decomposition_terms_after_J2"] == 3


def test_combined_post_j2_subset_preserves_prior_eight_slots_exactly() -> None:
    calc = KokunoPA10SourceR1AveragedLogradialOrdinarySlot()
    prior = calc.parent.combined_r1_after_j2_subset()
    combined = calc.combined_r1_after_j2_subset()
    for name, bound in prior.items():
        assert combined[name] == bound
    assert "lambda_inv_2D_eta_Au_times_Y_Phi_Y" not in prior
    assert "lambda_inv_2D_eta_Au_times_Y_Phi_Y" in combined


def test_report_is_deterministic_and_keeps_only_phi_eta_slots_open() -> None:
    calc = KokunoPA10SourceR1AveragedLogradialOrdinarySlot()
    first = calc.report()
    second = calc.report()
    assert first == second
    assert first["receipt_sha256"] == calc.sha256

    progress = first["ledger_progress"]
    assert progress["R1_post_J2_source_bound_ordinary_slots_count"] == 9
    assert progress["R1_total_ordinary_slots"] == 11
    assert progress["R1_ordinary_slots_still_required"] == 2
    assert tuple(progress["R1_remaining_ordinary_slots"]) == (
        "H_star_times_Phi_eta",
        "lambda_inv_d_u_times_Phi_eta",
    )
    assert progress["R1_numerator_space_source_bound_count_remains_six"] is True
    assert progress["R1_mixed_slot_handled_upstream_by_A1_668"] is True
    assert progress["full_R1_still_open"] is True
    assert progress["M_K_not_promoted"] is True

    truth = first["truth_boundary"]
    assert truth["source_averaging_identity_executable"] is True
    assert truth["source_single_logradial_post_J2_rule_executable"] is True
    assert truth["source_plain_post_J2_rule_executable"] is True
    assert truth[
        "source_compatible_R1_averaged_logradial_ordinary_post_J2_slot_machine_bound"
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
