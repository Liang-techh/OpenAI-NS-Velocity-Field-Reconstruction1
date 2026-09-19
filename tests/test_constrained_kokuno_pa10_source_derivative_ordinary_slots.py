from __future__ import annotations

import math

import numpy as np

from openai_ns_reconstruction.kokuno_pa10_post_j_remainder_bridge import (
    R2_ORDINARY_TERMS,
)
from openai_ns_reconstruction.kokuno_pa10_source_derivative_ordinary_slots import (
    R2_DERIVATIVE_ORDINARY_TERMS,
    KokunoPA10SourceDerivativeOrdinarySlots,
)


def test_fixed_Wstar_Hstar_polynomial_bounds_dominate_source_derivatives() -> None:
    calc = KokunoPA10SourceDerivativeOrdinarySlots()
    E = 1.0 + calc.domain.enlarged_real_margin
    eta = np.linspace(-E, E, 4001)
    h = calc.domain.h
    D = calc.domain.D
    j0 = calc.domain.j0
    rho = calc.rho

    W = -3.0 - 2.0 * D * j0 * eta + 8.0 * h * eta**2
    W1 = -2.0 * D * j0 + 16.0 * h * eta
    W2 = np.full_like(eta, 16.0 * h)
    W_ratios = [
        np.max(np.abs(W)),
        4.0 * rho * np.max(np.abs(W1)),
        4.5 * rho**2 * np.max(np.abs(W2)),
    ]
    assert calc.wstar_coefficient_ball().norm >= max(W_ratios)

    H = j0 + (D + 4.0) * eta - j0 * eta**2 - 4.0 * eta**3
    H1 = (D + 4.0) - 2.0 * j0 * eta - 12.0 * eta**2
    H2 = -2.0 * j0 - 24.0 * eta
    H3 = np.full_like(eta, -24.0)
    H_ratios = [
        np.max(np.abs(H)),
        4.0 * rho * np.max(np.abs(H1)),
        4.5 * rho**2 * np.max(np.abs(H2)),
        (8.0 / 3.0) * rho**3 * np.max(np.abs(H3)),
    ]
    assert calc.hstar_coefficient_ball().norm >= max(H_ratios)
    assert calc.wstar_coefficient_ball().lipschitz == 0.0
    assert calc.hstar_coefficient_ball().lipschitz == 0.0


def test_radial_average_is_contractively_bound_by_same_u_ball() -> None:
    calc = KokunoPA10SourceDerivativeOrdinarySlots()
    u = calc.bridge.u_source.u_ball()
    Au = calc.radial_average_u_ball()
    assert Au == u


def test_source_weight_ratio_supports_single_logradial_80_factor() -> None:
    # For J1[Y d_Y F], the exact weight ratio is
    # 20*a*(a+2)^2/((a+1)^3*(a+b+1)); it is safely <=80.
    worst = 0.0
    for alpha in range(1, 256):
        for beta in range(0, 64):
            value = (
                20.0
                * alpha
                * (alpha + 2.0) ** 2
                / ((alpha + 1.0) ** 3 * (alpha + beta + 1.0))
            )
            worst = max(worst, value)
    assert worst < 80.0
    assert math.isclose(worst, 11.25, rel_tol=0.0, abs_tol=1e-12)


def test_derivative_slots_fill_exactly_the_four_remaining_r2_slots() -> None:
    calc = KokunoPA10SourceDerivativeOrdinarySlots()
    slots = calc.r2_derivative_after_j1()
    assert tuple(slots) == R2_DERIVATIVE_ORDINARY_TERMS
    assert set(slots).issubset(set(R2_ORDINARY_TERMS))
    assert len(slots) == 4
    for bound in slots.values():
        assert math.isfinite(bound.norm) and bound.norm > 0.0
        assert math.isfinite(bound.lipschitz) and bound.lipschitz > 0.0


def test_eta_derivative_slots_are_rho_loss_but_no_standalone_u_eta_norm() -> None:
    calc = KokunoPA10SourceDerivativeOrdinarySlots()
    u = calc.bridge.u_source.u_ball()
    direct = calc._j1_eta_derivative(u)
    assert math.isclose(direct.norm / u.norm, 80.0 / calc.rho, rel_tol=2e-15)
    assert math.isclose(
        direct.lipschitz / u.lipschitz, 80.0 / calc.rho, rel_tol=2e-15
    )
    truth = calc.truth_boundary
    assert truth["standalone_u_eta_coefficient_norm_invented"] is False
    assert truth["standalone_Y_u_Y_coefficient_norm_invented"] is False


def test_quadratic_derivative_identities_produce_finite_post_j1_bounds() -> None:
    calc = KokunoPA10SourceDerivativeOrdinarySlots()
    avg_log = calc._j1_average_times_logradial_u()
    u_ueta = calc._j1_u_times_u_eta()
    for bound in (avg_log, u_ueta):
        assert math.isfinite(bound.norm) and bound.norm > 0.0
        assert math.isfinite(bound.lipschitz) and bound.lipschitz > 0.0
    assert u_ueta.norm > avg_log.norm


def test_full_r2_ordinary_post_j1_ledger_is_exactly_eleven_of_eleven() -> None:
    calc = KokunoPA10SourceDerivativeOrdinarySlots()
    all_slots = calc.all_r2_ordinary_after_j1()
    assert tuple(all_slots) == R2_ORDINARY_TERMS
    assert len(all_slots) == 11
    assert set(all_slots) == set(R2_ORDINARY_TERMS)
    for bound in all_slots.values():
        assert math.isfinite(bound.norm) and bound.norm > 0.0
        assert math.isfinite(bound.lipschitz) and bound.lipschitz > 0.0


def test_full_r2_self_certificate_adds_mixed_slot_without_double_j1() -> None:
    calc = KokunoPA10SourceDerivativeOrdinarySlots()
    all_slots = calc.all_r2_ordinary_after_j1()
    mixed = calc.bridge.full_mixed_terms()[
        "lambda_inv_Linv_d_detaAu_Y_u_Y_after_J1"
    ]
    total = calc.full_r2_after_j1()
    assert total.norm >= mixed.norm
    assert total.lipschitz >= mixed.lipschitz
    assert total.norm >= max(item.norm for item in all_slots.values())
    assert total.lipschitz >= max(item.lipschitz for item in all_slots.values())


def test_report_is_deterministic_closes_r2_only_and_keeps_science_gates_closed() -> None:
    calc = KokunoPA10SourceDerivativeOrdinarySlots()
    first = calc.report()
    second = calc.report()
    assert first == second
    assert first["receipt_sha256"] == calc.sha256

    progress = first["ledger_progress"]
    assert progress["R2_source_bound_ordinary_slots_count"] == 11
    assert progress["R2_total_ordinary_slots"] == 11
    assert progress["R2_ordinary_slots_still_required"] == 0
    assert progress["R1_ordinary_slots_still_required"] == 11
    assert progress["full_R2_agent1_self_certificate_materialized"] is True
    assert progress["full_R2_independent_admission"] is False
    assert progress["full_R1_still_open"] is True

    audit = first["independent_audit_boundary"]
    assert audit["A4_source_axis_independent_PASS_exists"] is True
    assert audit["A4_pressure_slots_independent_PASS_exists"] is True
    assert audit["A4_source_u_ball_independent_PASS_exists"] is True
    assert audit["A1_algebraic_slots_still_requires_independent_A4_audit"] is True
    assert audit["this_derivative_subset_requires_independent_A4_audit"] is True
    assert audit["full_R2_independent_A4_admission"] is False

    truth = first["truth_boundary"]
    assert truth["source_compatible_all_R2_ordinary_post_J1_slots_machine_bound"] is True
    assert truth[
        "source_compatible_full_post_J1_R2_radius_one_ball_bound_agent1_self_certificate"
    ] is True
    assert truth["full_R2_independent_agent4_admission"] is False
    assert truth["all_R1_ordinary_slots_source_bound"] is False
    assert truth["source_full_post_J2_R1_radius_one_ball_bound_machine_bound"] is False
    assert truth["source_operator_constant_M_machine_bound"] is False
    assert truth["source_operator_constant_K_machine_bound"] is False
    assert truth["global_pressure_matched"] is False
    assert truth["global_leading_profile_reconstructed"] is False
    assert truth["heldout_ns_residual_assessed"] is False
    assert truth["pde_validated"] is False
