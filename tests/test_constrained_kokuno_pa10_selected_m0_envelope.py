from __future__ import annotations

import copy
import math

import pytest

from openai_ns_reconstruction.kokuno_pa10_selected_m0_envelope import (
    KokunoPA10SelectedM0Envelope,
)


@pytest.fixture(scope="module")
def selected_m0() -> KokunoPA10SelectedM0Envelope:
    return KokunoPA10SelectedM0Envelope()


def test_selected_reference_m0_envelope_is_finite_and_refinement_guarded(selected_m0):
    report = selected_m0.report()
    numerical = report["numerical_envelope"]
    coarse = numerical["coarse"]
    fine = numerical["fine"]

    assert math.isfinite(fine["sample_max_abs_p1_reference"])
    assert fine["sample_max_abs_p1_reference"] > 3.0
    assert fine["sample_min_p1_reference"] > 0.0
    assert 0.0 <= fine["sample_max_X"] <= 110.0
    assert numerical["M0_candidate_upper_envelope"] > fine["sample_max_abs_p1_reference"]
    assert numerical["M0_candidate_upper_envelope"] > coarse["sample_max_abs_p1_reference"]
    assert numerical["nested_grid_relative_drift"] <= numerical[
        "max_allowed_refinement_drift_fraction"
    ]
    assert numerical["numerical_envelope_guard_passed"] is True
    assert numerical["continuum_supremum_proved"] is False


def test_candidate_envelope_feeds_displayed_B0_formula_without_source_promotion(selected_m0):
    value = selected_m0.candidate_B0_diagnostic(combined_profile_C0_norm=1.0)
    assert math.isfinite(value)
    assert value > selected_m0.envelope["M0_candidate_upper_envelope"]

    truth = selected_m0.report()["truth_boundary"]
    assert truth["selected_reference_M0_numerical_envelope_materialized"] is True
    assert truth["selected_reference_is_source_fixed_point"] is False
    assert truth["selected_M0_envelope_is_continuum_interval_proof"] is False
    assert truth["source_M0_machine_bound"] is False
    assert truth["source_B0_dependencies_machine_bound"] is False
    assert truth["source_B0_analytic_bound_proved"] is False
    assert truth["source_T_sh_lower_bound_verified"] is False
    assert truth["selected_pa16_handoff_allowed"] is False
    assert truth["heldout_ns_residual_assessed"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False


def test_source_definition_and_full_log_radial_interval_are_recorded(selected_m0):
    report = selected_m0.report()
    formulas = report["source_formulas"]
    geometry = report["geometry"]

    assert formulas["reference_primitive"] == "D_X p1_r=X*S_q,r/L-l_r*p1_r"
    assert geometry["X_0"] > 0.0
    assert geometry["X_i"] == 110.0
    assert geometry["y_i"] > geometry["reference_transition_y_end"] > 0.0
    assert report["resolved_candidate_side_B0_inputs"] == [
        "L_0",
        "selected M_0 numerical envelope",
    ]
    assert report["unresolved_source_B0_dependencies"] == [
        "||log phi_*+log Phi(4,.)||_{C^0_eta}",
        "M_0 analytic/source-fixed-point bound",
    ]


def test_payload_round_trip_and_truth_boundary_fail_closed(selected_m0):
    payload = selected_m0.to_payload()
    replay = KokunoPA10SelectedM0Envelope.from_payload(payload)
    assert replay.to_payload() == payload
    assert len(selected_m0.sha256) == 64

    tampered = copy.deepcopy(payload)
    tampered["truth_boundary"]["source_M0_machine_bound"] = True
    with pytest.raises(ValueError, match="truth-boundary"):
        KokunoPA10SelectedM0Envelope.from_payload(tampered)

    tampered = copy.deepcopy(payload)
    tampered["source_formulas"]["reference_primitive"] = "p1_r=0"
    with pytest.raises(ValueError, match="source formulas"):
        KokunoPA10SelectedM0Envelope.from_payload(tampered)
