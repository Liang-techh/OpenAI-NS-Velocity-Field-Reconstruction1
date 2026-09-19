from __future__ import annotations

import math

from openai_ns_reconstruction.kokuno_pa10_post_j_remainder_bridge import (
    R1_ORDINARY_TERMS,
)
from openai_ns_reconstruction.kokuno_pa10_source_r1_algebraic_ordinary_slots import (
    R1_ALGEBRAIC_ORDINARY_TERMS,
    R1_REMAINING_ORDINARY_TERMS,
    KokunoPA10SourceR1AlgebraicOrdinarySlots,
)


def test_r1_algebraic_subset_is_exactly_five_registered_slots() -> None:
    calc = KokunoPA10SourceR1AlgebraicOrdinarySlots()
    slots = calc.r1_algebraic_numerator_slots()
    assert tuple(slots) == R1_ALGEBRAIC_ORDINARY_TERMS
    assert len(slots) == 5
    assert set(slots).issubset(set(R1_ORDINARY_TERMS))
    assert set(slots).isdisjoint(set(R1_REMAINING_ORDINARY_TERMS))
    assert set(slots) | set(R1_REMAINING_ORDINARY_TERMS) == set(R1_ORDINARY_TERMS)


def test_r1_algebraic_slots_consume_nontrivial_phi_and_u_balls() -> None:
    calc = KokunoPA10SourceR1AlgebraicOrdinarySlots()
    phi = calc.phi_ball()
    u = calc.u_ball()
    Au = calc.radial_average_u_ball()
    assert phi.norm > 1.0
    assert phi.lipschitz == 1.0
    assert u.norm > 1.0
    assert u.lipschitz == 1.0
    assert Au == u


def test_fixed_multiplier_inputs_are_source_compatible_and_finite() -> None:
    calc = KokunoPA10SourceR1AlgebraicOrdinarySlots()
    fixed = calc.fixed_multiplier_balls()
    assert set(fixed) == {"eta", "U_star", "W_star", "L_inverse"}
    for bound in fixed.values():
        assert math.isfinite(bound.norm) and bound.norm > 0.0
        assert bound.lipschitz == 0.0
    assert calc.domain.coefficient_rho < calc.domain.cauchy_radius
    assert 0.0 < calc.domain.h < 0.01
    assert calc.domain.D > 0.0


def test_five_r1_numerator_bounds_are_finite_nonzero_and_lipschitz() -> None:
    calc = KokunoPA10SourceR1AlgebraicOrdinarySlots()
    slots = calc.r1_algebraic_numerator_slots()
    for bound in slots.values():
        assert math.isfinite(bound.norm) and bound.norm > 0.0
        assert math.isfinite(bound.lipschitz) and bound.lipschitz > 0.0


def test_post_j2_bridge_applies_once_to_exact_same_five_slots() -> None:
    calc = KokunoPA10SourceR1AlgebraicOrdinarySlots()
    numerator = calc.r1_algebraic_numerator_slots()
    after = calc.r1_algebraic_after_j2()
    assert tuple(after) == R1_ALGEBRAIC_ORDINARY_TERMS
    for name in R1_ALGEBRAIC_ORDINARY_TERMS:
        direct = calc.bridge.ordinary_term_after_jnu(numerator[name], nu=2)
        assert after[name] == direct
        assert after[name].norm > numerator[name].norm
        assert after[name].lipschitz > numerator[name].lipschitz


def test_uniform_lambda_inverse_is_only_upper_not_hidden_lambda_choice() -> None:
    calc = KokunoPA10SourceR1AlgebraicOrdinarySlots()
    report = calc.report()
    selected = report["selected_source_compatible_inputs"]
    assert selected["Lambda_inverse_upper"] == 1.0
    assert "Lambda" not in selected
    truth = report["truth_boundary"]
    assert truth["uniform_Lambda_inverse_upper_for_Lambda_ge_1_used"] is True


def test_report_is_deterministic_and_keeps_six_r1_slots_open() -> None:
    calc = KokunoPA10SourceR1AlgebraicOrdinarySlots()
    first = calc.report()
    second = calc.report()
    assert first == second
    assert first["receipt_sha256"] == calc.sha256

    progress = first["ledger_progress"]
    assert progress["R1_source_bound_ordinary_slots_count"] == 5
    assert progress["R1_total_ordinary_slots"] == 11
    assert progress["R1_ordinary_slots_still_required"] == 6
    assert tuple(progress["R1_remaining_ordinary_slots"]) == R1_REMAINING_ORDINARY_TERMS
    assert progress["R1_mixed_slot_handled_upstream_by_A1_668"] is True
    assert progress["R2_agent1_self_certificate_materialized_upstream"] is True
    assert progress["full_R1_still_open"] is True
    assert progress["M_K_not_promoted"] is True

    audit = first["independent_audit_boundary"]
    assert audit["A4_source_axis_independent_PASS_exists"] is True
    assert audit["A4_source_Phi_ball_conditional_PASS_exists"] is True
    assert audit["A4_source_u_ball_independent_PASS_exists"] is True
    assert audit["R2_full_independent_A4_admission"] is False
    assert audit["this_R1_subset_requires_independent_A4_audit"] is True

    truth = first["truth_boundary"]
    assert truth["source_compatible_R1_algebraic_ordinary_post_J2_slots_machine_bound"] is True
    assert truth["standalone_Y_Phi_Y_coefficient_norm_invented"] is False
    assert truth["standalone_Phi_eta_coefficient_norm_invented"] is False
    assert truth["source_compatible_all_R1_ordinary_post_J2_slots_machine_bound"] is False
    assert truth["source_full_post_J2_R1_radius_one_ball_bound_machine_bound"] is False
    assert truth["source_operator_constant_M_machine_bound"] is False
    assert truth["source_operator_constant_K_machine_bound"] is False
    assert truth["global_pressure_matched"] is False
    assert truth["global_leading_profile_reconstructed"] is False
    assert truth["heldout_ns_residual_assessed"] is False
    assert truth["pde_validated"] is False
