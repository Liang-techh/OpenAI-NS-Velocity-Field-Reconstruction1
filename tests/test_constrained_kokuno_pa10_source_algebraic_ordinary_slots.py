from __future__ import annotations

from dataclasses import asdict
import math

import numpy as np

from openai_ns_reconstruction.kokuno_pa10_post_j_remainder_bridge import (
    R2_ORDINARY_TERMS,
)
from openai_ns_reconstruction.kokuno_pa10_source_algebraic_ordinary_slots import (
    R2_ALGEBRAIC_ORDINARY_TERMS,
    R2_REMAINING_DERIVATIVE_ORDINARY_TERMS,
    KokunoPA10SourceAlgebraicOrdinarySlots,
)
from openai_ns_reconstruction.kokuno_pa10_source_pressure_ordinary_slots import (
    R2_PRESSURE_ORDINARY_TERMS,
)


def test_axis_multiplier_values_are_vectorized_and_match_source_formulas() -> None:
    calc = KokunoPA10SourceAlgebraicOrdinarySlots()
    eta = np.array([-0.7, 0.0, 0.4])
    values = calc.axis_multiplier_values(eta)
    assert np.allclose(values["eta"], eta, rtol=0.0, atol=0.0)
    assert np.allclose(values["A"], 0.5 + calc.domain.h, rtol=0.0, atol=0.0)
    assert np.allclose(values["U_star"], 4.0 * eta + calc.domain.j0, rtol=0.0, atol=0.0)
    assert np.allclose(values["U_star_eta"], 4.0, rtol=0.0, atol=0.0)
    assert np.allclose(values["d"], 1.0 - eta * eta, rtol=0.0, atol=0.0)


def test_axis_multiplier_values_reject_outside_certified_interval() -> None:
    calc = KokunoPA10SourceAlgebraicOrdinarySlots()
    E = 1.0 + calc.domain.enlarged_real_margin
    try:
        calc.axis_multiplier_values(E + 1e-3)
    except ValueError as exc:
        assert "certified enlarged real interval" in str(exc)
    else:
        raise AssertionError("outside-interval eta must fail closed")


def test_ustar_coefficient_ball_matches_alpha0_source_weight() -> None:
    calc = KokunoPA10SourceAlgebraicOrdinarySlots()
    bound = calc.ustar_coefficient_ball()
    expected = max(
        4.0 * (1.0 + calc.domain.enlarged_real_margin) + calc.domain.j0,
        16.0 * calc.domain.coefficient_rho,
    )
    assert bound.norm >= expected
    assert bound.norm <= math.nextafter(expected, math.inf)
    assert bound.lipschitz == 0.0


def test_algebraic_slots_fill_only_registered_r2_subset() -> None:
    calc = KokunoPA10SourceAlgebraicOrdinarySlots()
    slots = calc.r2_algebraic_numerator_slots()
    assert tuple(slots) == R2_ALGEBRAIC_ORDINARY_TERMS
    assert set(slots).issubset(set(R2_ORDINARY_TERMS))
    assert set(slots).isdisjoint(set(R2_PRESSURE_ORDINARY_TERMS))
    assert len(slots) == 4
    for bound in slots.values():
        assert math.isfinite(bound.norm) and bound.norm > 0.0
        assert math.isfinite(bound.lipschitz) and bound.lipschitz > 0.0


def test_four_algebraic_slots_replay_public_product_algebra() -> None:
    calc = KokunoPA10SourceAlgebraicOrdinarySlots()
    slots = calc.r2_algebraic_numerator_slots()
    fixed = calc.fixed_multiplier_balls()
    u = calc.bridge.u_source.u_ball()
    product = calc.product_calculator
    A = 0.5 + calc.domain.h

    expected_Au = product.scale(u, A)
    expected_eta_U_u = product.scale(
        product.product(fixed["eta"], fixed["U_star"], u), 4.0 * A
    )
    expected_d_u = product.scale(product.product(fixed["d"], u), 4.0)
    expected_eta_u2 = product.scale(
        product.product(fixed["eta"], u, u), 2.0 * A
    )

    # The implementation uses exact-Fraction scaling while the public helper
    # treats the same binary64 scalar as exact.  They must therefore replay.
    assert asdict(slots["A_times_u"]) == asdict(expected_Au)
    assert asdict(slots["4A_eta_Ustar_times_u"]) == asdict(expected_eta_U_u)
    assert asdict(slots["d_Ustar_eta_times_u"]) == asdict(expected_d_u)
    assert asdict(slots["2A_lambda_inv_eta_u_squared"]) == asdict(expected_eta_u2)


def test_post_j1_applies_bridge_exactly_once() -> None:
    calc = KokunoPA10SourceAlgebraicOrdinarySlots()
    numerator = calc.r2_algebraic_numerator_slots()
    after = calc.r2_algebraic_after_j1()
    assert tuple(after) == R2_ALGEBRAIC_ORDINARY_TERMS
    for name in R2_ALGEBRAIC_ORDINARY_TERMS:
        direct = calc.bridge.ordinary_term_after_jnu(numerator[name], nu=1)
        assert asdict(after[name]) == asdict(direct)
        assert after[name].norm > numerator[name].norm
        assert after[name].lipschitz > numerator[name].lipschitz


def test_combined_pressure_and_algebraic_subset_is_exactly_seven_of_eleven() -> None:
    calc = KokunoPA10SourceAlgebraicOrdinarySlots()
    combined = calc.combined_source_bound_r2_subset()
    assert len(combined) == 7
    assert set(combined) == set(R2_PRESSURE_ORDINARY_TERMS).union(
        R2_ALGEBRAIC_ORDINARY_TERMS
    )
    assert set(R2_ORDINARY_TERMS) - set(combined) == set(
        R2_REMAINING_DERIVATIVE_ORDINARY_TERMS
    )


def test_report_is_deterministic_and_fail_closed_beyond_seven_slots() -> None:
    calc = KokunoPA10SourceAlgebraicOrdinarySlots()
    first = calc.report()
    second = calc.report()
    assert first == second
    assert first["receipt_sha256"] == calc.sha256
    progress = first["ledger_progress"]
    assert progress["R2_source_bound_ordinary_slots_count"] == 7
    assert progress["R2_total_ordinary_slots"] == 11
    assert progress["R2_ordinary_slots_still_required"] == 4
    assert progress["R1_ordinary_slots_still_required"] == 11
    assert progress["R2_remaining_derivative_slots"] == list(
        R2_REMAINING_DERIVATIVE_ORDINARY_TERMS
    )

    independent = first["independent_audit_boundary"]
    assert independent["A4_source_axis_independent_PASS_exists"] is True
    assert independent["A4_pressure_slots_independent_PASS_exists"] is True
    assert independent["this_new_algebraic_subset_requires_independent_A4_audit"] is True

    truth = first["truth_boundary"]
    assert truth["source_compatible_R2_algebraic_ordinary_post_J1_slots_machine_bound"] is True
    assert truth["all_R2_ordinary_slots_source_bound"] is False
    assert truth["all_R1_ordinary_slots_source_bound"] is False
    assert truth["source_full_post_J1_R2_radius_one_ball_bound_machine_bound"] is False
    assert truth["source_operator_constant_M_machine_bound"] is False
    assert truth["global_pressure_matched"] is False
    assert truth["global_leading_profile_reconstructed"] is False
    assert truth["heldout_ns_residual_assessed"] is False
    assert truth["pde_validated"] is False
