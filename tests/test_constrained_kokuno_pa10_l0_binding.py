from __future__ import annotations

import copy
import math

import pytest

from openai_ns_reconstruction.kokuno_pa10_l0_binding import KokunoPA10L0Binding
from openai_ns_reconstruction.kokuno_pa10_source_choice_order_guard import (
    KokunoPA10SourceChoiceOrderGuard,
)


def test_selected_lambda_binds_source_X0_and_L0_exactly():
    binding = KokunoPA10L0Binding()
    replay = binding.verify_identity()

    assert binding.selected_Lambda > 0.0
    assert binding.X_0 == pytest.approx(4.0 / binding.selected_Lambda, rel=0.0, abs=0.0)
    assert binding.L_0 == pytest.approx(
        math.log(binding.X_i / 4.0) + math.log(binding.selected_Lambda),
        rel=2e-16,
    )
    assert replay["identity_within_8ulp"] is True
    assert replay["identity_abs_error"] <= 8.0 * math.ulp(binding.L_0)


def test_l0_is_the_only_newly_resolved_B0_dependency_and_pa16_stays_closed():
    binding = KokunoPA10L0Binding()
    report = binding.report()
    truth = report["truth_boundary"]

    assert report["resolved_source_B0_dependencies"] == ["L_0"]
    assert report["unresolved_source_B0_dependencies"] == [
        "||log phi_*+log Phi(4,.)||_{C^0_eta}",
        "M_0",
    ]
    assert truth["source_L0_machine_bound"] is True
    assert truth["selected_Lambda_is_autonomous_conditioning_choice"] is True
    assert truth["source_hidden_Lambda_recovered"] is False
    assert truth["combined_profile_C0_norm_machine_bound"] is False
    assert truth["source_M0_machine_bound"] is False
    assert truth["source_B0_dependencies_machine_bound"] is False
    assert truth["source_B0_analytic_bound_proved"] is False
    assert truth["source_T_sh_lower_bound_verified"] is False
    assert truth["selected_pa16_handoff_allowed"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False


def test_bound_l0_feeds_displayed_B0_formula_without_promoting_it():
    binding = KokunoPA10L0Binding()
    guard = KokunoPA10SourceChoiceOrderGuard(X_i=binding.X_i)

    diagnostic = guard.displayed_B0_upper_bound(
        combined_profile_C0_norm=1.0,
        M0=2.0,
        L0=binding.L_0,
    )
    expected = 1.0 + 0.5 * (2.0 + 0.8) * binding.L_0 + guard.half_log_2X_i
    assert diagnostic == pytest.approx(expected, rel=2e-16)

    # The executable substitution is only a diagnostic until the other two
    # source dependencies are genuinely bounded.
    assert binding.report()["truth_boundary"]["source_B0_analytic_bound_proved"] is False
    assert guard.report()["truth_boundary"]["source_B0_analytic_bound_proved"] is False


def test_payload_round_trip_and_truth_boundary_fail_closed():
    binding = KokunoPA10L0Binding()
    payload = binding.to_payload()
    replay = KokunoPA10L0Binding.from_payload(payload)
    assert replay.to_payload() == payload
    assert len(binding.sha256) == 64

    tampered = copy.deepcopy(payload)
    tampered["truth_boundary"]["selected_pa16_handoff_allowed"] = True
    with pytest.raises(ValueError, match="truth-boundary"):
        KokunoPA10L0Binding.from_payload(tampered)

    tampered = copy.deepcopy(payload)
    tampered["source_formulas"]["inner_radius"] = "X_0=1/Lambda"
    with pytest.raises(ValueError, match="source formulas"):
        KokunoPA10L0Binding.from_payload(tampered)
