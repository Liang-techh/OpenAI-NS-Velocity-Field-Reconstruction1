from __future__ import annotations

import copy
import math

import pytest

from openai_ns_reconstruction.kokuno_coupled_c_appendix_b_normalization import (
    KokunoCoupledCAppendixBNormalization,
)
from openai_ns_reconstruction.kokuno_pa10_source_choice_order_guard import (
    KokunoPA10SourceChoiceOrderGuard,
)


def test_source_ell_i_cancels_later_normalization_C_exactly():
    guard = KokunoPA10SourceChoiceOrderGuard()
    log_phi_i = -7.25

    small = guard.verify_C_cancellation(log_phi_i, math.log(2.0))
    huge = guard.verify_C_cancellation(log_phi_i, 4.5e23, atol=1.0e-7)

    assert small["cancellation_passed"] is True
    assert huge["cancellation_passed"] is True
    assert small["ell_i_direct_C_independent"] == pytest.approx(
        guard.source_ell_i(log_phi_i), abs=1e-15
    )
    assert huge["ell_i_direct_C_independent"] == pytest.approx(
        guard.source_ell_i(log_phi_i), abs=1e-15
    )


def test_displayed_B0_and_Tsh_formulas_are_executable_but_not_promoted():
    guard = KokunoPA10SourceChoiceOrderGuard()
    B0 = guard.displayed_B0_upper_bound(
        combined_profile_C0_norm=1.25,
        M0=2.0,
        L0=3.0,
    )
    expected = 1.25 + 0.5 * (2.0 + 0.8) * 3.0 + 0.5 * math.log(220.0)
    assert B0 == pytest.approx(expected, rel=2e-16)
    assert guard.displayed_T_sh_lower_bound_from_B0(B0) == pytest.approx(
        20.0 * guard.sigma_prime_sup * (B0 + guard.log_f_sup), rel=2e-16
    )

    truth = guard.report()["truth_boundary"]
    assert truth["displayed_B0_bound_formula_executable"] is True
    assert truth["source_B0_dependencies_machine_bound"] is False
    assert truth["source_B0_analytic_bound_proved"] is False
    assert truth["source_T_sh_lower_bound_verified"] is False
    assert truth["selected_pa16_handoff_allowed"] is False

    with pytest.raises(ValueError, match="M0"):
        guard.displayed_B0_upper_bound(1.0, -1.0, 2.0)
    with pytest.raises(ValueError, match="L0"):
        guard.displayed_B0_upper_bound(1.0, 2.0, float("inf"))


@pytest.fixture(scope="module")
def coupled_candidate() -> KokunoCoupledCAppendixBNormalization:
    return KokunoCoupledCAppendixBNormalization(log_C_factor=32.0)


def test_factor32_pointwise_pass_does_not_promote_source_pa10_or_pa16(
    coupled_candidate,
):
    guard = KokunoPA10SourceChoiceOrderGuard()
    classification = guard.classify_coupled_candidate(coupled_candidate)

    # Preserve the actual #521 autonomous result while enforcing the newly
    # reread public-source dependency order.
    assert classification["autonomous_coupled_C_pointwise_screen_passed"] is True
    assert classification["autonomous_pointwise_strict_margin"] > 0.0
    assert classification[
        "source_choice_order_allows_C_to_reparameterize_upstream_phi"
    ] is False
    assert classification[
        "source_choice_order_allows_C_to_discharge_B0_or_T_sh"
    ] is False
    assert classification["source_B0_dependencies_machine_bound"] is False
    assert classification["source_B0_analytic_bound_proved"] is False
    assert classification["source_T_sh_lower_bound_verified"] is False
    assert classification["selected_pa16_handoff_allowed"] is False


def test_payload_round_trip_and_source_order_truth_boundary_fail_closed():
    guard = KokunoPA10SourceChoiceOrderGuard()
    payload = guard.to_payload()
    replay = KokunoPA10SourceChoiceOrderGuard.from_payload(payload)
    assert replay.to_payload() == payload
    assert len(guard.sha256) == 64

    tampered = copy.deepcopy(payload)
    tampered["source_choice_order"][3] = "C_and_X_R"
    with pytest.raises(ValueError, match="choice order"):
        KokunoPA10SourceChoiceOrderGuard.from_payload(tampered)

    tampered = copy.deepcopy(payload)
    tampered["truth_boundary"]["selected_pa16_handoff_allowed"] = True
    with pytest.raises(ValueError, match="truth-boundary"):
        KokunoPA10SourceChoiceOrderGuard.from_payload(tampered)
