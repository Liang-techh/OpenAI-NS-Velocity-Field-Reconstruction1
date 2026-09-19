from __future__ import annotations

import math

import numpy as np

from openai_ns_reconstruction.kokuno_pa10_post_j_remainder_bridge import (
    R1_ORDINARY_TERMS,
)
from openai_ns_reconstruction.kokuno_pa10_source_r1_zeta_ordinary_slot import (
    R1_REMAINING_AFTER_ZETA_TERMS,
    R1_SOURCE_BOUND_THROUGH_ZETA_TERMS,
    R1_ZETA_ORDINARY_TERMS,
    KokunoPA10SourceR1ZetaOrdinarySlot,
)


def test_zeta_increment_is_exactly_one_registered_r1_slot() -> None:
    calc = KokunoPA10SourceR1ZetaOrdinarySlot()
    slots = calc.r1_zeta_numerator_slot()
    assert tuple(slots) == R1_ZETA_ORDINARY_TERMS
    assert set(slots).issubset(set(R1_ORDINARY_TERMS))
    combined = calc.combined_r1_source_bound_subset()
    assert tuple(combined) == R1_SOURCE_BOUND_THROUGH_ZETA_TERMS
    assert len(combined) == 6
    assert len(R1_REMAINING_AFTER_ZETA_TERMS) == 5
    assert set(combined).isdisjoint(set(R1_REMAINING_AFTER_ZETA_TERMS))
    assert set(combined) | set(R1_REMAINING_AFTER_ZETA_TERMS) == set(R1_ORDINARY_TERMS)


def test_zeta_value_api_replays_public_axis_formula() -> None:
    calc = KokunoPA10SourceR1ZetaOrdinarySlot()
    eta = np.array([-1.0, -0.731, -0.125, 0.0, 0.417, 0.999])
    public = calc.zeta_star_values(eta)
    axis = np.asarray(calc.domain.axis_state(eta)["zeta_star"], dtype=float)
    np.testing.assert_array_equal(public, axis)
    assert np.all(np.isfinite(public))
    assert float(np.max(np.abs(public))) > 0.0


def test_zeta_cauchy_certificate_is_finite_positive_and_nonvacuous() -> None:
    calc = KokunoPA10SourceR1ZetaOrdinarySlot()
    cert = calc.zeta_star_certificate()
    for value in cert.values():
        assert math.isfinite(value) and value > 0.0
    assert cert["rho_over_cauchy_radius"] <= 0.10000000000000002
    assert cert["H_factor_abs_lower_on_cauchy"] > 0.0
    assert cert["H_square_plus_sigma_square_abs_lower_on_cauchy"] > 0.0
    assert cert["zeta_star_coefficient_norm_upper"] >= cert["zeta_star_complex_sup_upper"]


def test_zeta_coefficient_bound_dominates_fresh_real_value_stress() -> None:
    calc = KokunoPA10SourceR1ZetaOrdinarySlot()
    edge = 1.0 + float(calc.domain.enlarged_real_margin)
    eta = np.linspace(-edge, edge, 4097)
    observed = float(np.max(np.abs(calc.zeta_star_values(eta))))
    bound = calc.zeta_star_coefficient_ball()
    assert observed > 0.0
    assert observed <= bound.norm
    assert bound.lipschitz == 0.0


def test_zeta_r1_numerator_consumes_d_u_phi_and_is_nontrivial() -> None:
    calc = KokunoPA10SourceR1ZetaOrdinarySlot()
    slot = calc.r1_zeta_numerator_slot()["d_u_zeta_star_times_Phi"]
    assert math.isfinite(slot.norm) and slot.norm > 0.0
    assert math.isfinite(slot.lipschitz) and slot.lipschitz > 0.0
    assert calc.algebraic_r1.u_ball().norm > 1.0
    assert calc.algebraic_r1.phi_ball().norm > 1.0
    assert calc.algebraic_r1.algebraic.fixed_multiplier_balls()["d"].norm > 0.0


def test_post_j2_is_applied_exactly_once_to_zeta_slot() -> None:
    calc = KokunoPA10SourceR1ZetaOrdinarySlot()
    numerator = calc.r1_zeta_numerator_slot()["d_u_zeta_star_times_Phi"]
    after = calc.r1_zeta_after_j2()["d_u_zeta_star_times_Phi"]
    direct = calc.bridge.ordinary_term_after_jnu(numerator, nu=2)
    assert after == direct
    assert after.norm > numerator.norm
    assert after.lipschitz > numerator.lipschitz
    combined_after = calc.combined_r1_after_j2_subset()
    assert tuple(combined_after) == R1_SOURCE_BOUND_THROUGH_ZETA_TERMS
    assert combined_after["d_u_zeta_star_times_Phi"] == after


def test_report_is_deterministic_and_keeps_five_r1_slots_open() -> None:
    calc = KokunoPA10SourceR1ZetaOrdinarySlot()
    first = calc.report()
    second = calc.report()
    assert first == second
    assert first["receipt_sha256"] == calc.sha256

    progress = first["ledger_progress"]
    assert progress["R1_source_bound_ordinary_slots_count"] == 6
    assert progress["R1_total_ordinary_slots"] == 11
    assert progress["R1_ordinary_slots_still_required"] == 5
    assert tuple(progress["R1_remaining_ordinary_slots"]) == R1_REMAINING_AFTER_ZETA_TERMS
    assert progress["R1_mixed_slot_handled_upstream_by_A1_668"] is True
    assert progress["R2_agent1_self_certificate_materialized_upstream"] is True
    assert progress["full_R1_still_open"] is True
    assert progress["M_K_not_promoted"] is True

    audit = first["independent_audit_boundary"]
    assert audit["this_zeta_R1_slot_is_agent1_self_certificate"] is True
    assert audit["this_zeta_R1_slot_requires_independent_A4_audit"] is True
    assert audit["R2_full_independent_A4_admission"] is False

    truth = first["truth_boundary"]
    assert truth["source_compatible_zeta_star_coefficient_norm_machine_bound"] is True
    assert truth["source_compatible_R1_zeta_ordinary_post_J2_slot_machine_bound"] is True
    assert truth["source_compatible_R1_ordinary_source_bound_count_is_six"] is True
    assert truth["source_compatible_all_R1_ordinary_post_J2_slots_machine_bound"] is False
    assert truth["standalone_Y_Phi_Y_coefficient_norm_invented"] is False
    assert truth["standalone_Phi_eta_coefficient_norm_invented"] is False
    assert truth["source_operator_constant_M_machine_bound"] is False
    assert truth["source_operator_constant_K_machine_bound"] is False
    assert truth["global_leading_profile_reconstructed"] is False
    assert truth["heldout_ns_residual_assessed"] is False
    assert truth["pde_validated"] is False
