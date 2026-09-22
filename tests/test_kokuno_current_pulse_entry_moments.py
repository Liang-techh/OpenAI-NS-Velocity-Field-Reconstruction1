from __future__ import annotations

import copy
import math

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_current_pulse_entry_moments import (
    MOMENT_ORDER,
    PARENT_EXACT_HEAD,
    KokunoCurrentPulseEntryMoments,
)


def test_public_J_definition_and_current_truth_boundary() -> None:
    candidate = KokunoCurrentPulseEntryMoments()
    formulas = candidate.source_formulas
    truth = candidate.truth_boundary

    assert MOMENT_ORDER == ("M", "J", "I", "S", "C_p")
    assert "J=int_0^X U H dx" in formulas["moments"]
    assert formulas["row_slopes"] == "s1=1/2-lambda; s2=1/2-2lambda"
    assert candidate.i1_gate.closure_tolerance == 5.0e-7
    assert truth["current_lineage_J_entry_materialized"] is True
    assert truth["current_J_from_frozen_I1_closure_residual"] is True
    assert truth["current_J_assumed_zero"] is False
    assert truth["source_exact_amplitude_root_materialized"] is False
    assert truth["current_cartesian_end_compensation_composed"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False
    assert truth["openai_field_identified"] is False


def test_entry_ratios_are_exact_conversion_of_frozen_I1_residual() -> None:
    candidate = KokunoCurrentPulseEntryMoments()
    eta = np.asarray([-0.71, -0.37, 0.19, 0.63], dtype=float)
    values = candidate.entry_ratios(eta)

    residual = np.asarray(values["post_I1_residual_normalized"], dtype=float)
    f = candidate.source_f(eta)
    delta_log = candidate.repair.log_X_1 - candidate.leading.log_X_p
    scale_m = math.exp(candidate.compensator.s1 * delta_log)
    scale_j = math.exp(candidate.compensator.s2 * delta_log) / math.sqrt(2.0)

    expected_m = scale_m * residual[..., 0] / f
    expected_j = scale_j * residual[..., 1] / (f * f)
    assert np.array_equal(values["m_entry_ratio"], expected_m)
    assert np.array_equal(values["j_entry_ratio"], expected_j)
    assert np.all(np.isfinite(values["M_physical_at_pulse_entry"]))
    assert np.all(np.isfinite(values["J_physical_at_pulse_entry"]))

    # The current M primitive supplies an independent implementation path.  Do
    # not force exact bit identity because the closure and primitive routes use
    # different quadrature/bookkeeping; this is a diagnostic rather than PDE evidence.
    direct = np.asarray(values["m_entry_ratio_direct_current_primitive"], dtype=float)
    assert np.all(np.isfinite(direct))
    assert np.all(np.isfinite(values["m_entry_crosscheck_relative_difference"]))


def test_analytic_eta_jets_match_centered_replay_and_symmetry() -> None:
    candidate = KokunoCurrentPulseEntryMoments()
    eta = np.asarray([0.23, 0.47, 0.76], dtype=float)
    h = 2.0e-6

    center = candidate.entry_ratios(eta)
    plus = candidate.entry_ratios(eta + h)
    minus = candidate.entry_ratios(eta - h)
    fd_m = (plus["m_entry_ratio"] - minus["m_entry_ratio"]) / (2.0 * h)
    fd_j = (plus["j_entry_ratio"] - minus["j_entry_ratio"]) / (2.0 * h)

    np.testing.assert_allclose(
        center["m_entry_ratio_eta"], fd_m, rtol=3.0e-4, atol=2.0e-13
    )
    np.testing.assert_allclose(
        center["j_entry_ratio_eta"], fd_j, rtol=5.0e-4, atol=2.0e-13
    )

    positive = candidate.entry_ratios(eta)
    negative = candidate.entry_ratios(-eta)
    np.testing.assert_allclose(
        positive["m_entry_ratio"], negative["m_entry_ratio"], rtol=2.0e-11, atol=2.0e-14
    )
    np.testing.assert_allclose(
        positive["j_entry_ratio"], negative["j_entry_ratio"], rtol=2.0e-11, atol=2.0e-14
    )
    np.testing.assert_allclose(
        positive["m_entry_ratio_eta"], -negative["m_entry_ratio_eta"], rtol=3.0e-9, atol=2.0e-13
    )
    np.testing.assert_allclose(
        positive["j_entry_ratio_eta"], -negative["j_entry_ratio_eta"], rtol=3.0e-9, atol=2.0e-13
    )


def test_current_entry_data_feed_public_end_algebra_without_promoting_velocity() -> None:
    candidate = KokunoCurrentPulseEntryMoments()
    eta = np.asarray([-0.63, -0.21, 0.0, 0.41, 0.79], dtype=float)
    inputs = candidate.compensator_inputs(eta)
    solution, c1_eta, c2_eta = candidate.compensator.solve_eta_jet(
        inputs["amplitude"],
        inputs["m_entry_ratio"],
        inputs["j_entry_ratio"],
        inputs["amplitude_eta"],
        inputs["m_entry_ratio_eta"],
        inputs["j_entry_ratio_eta"],
    )

    assert np.max(np.abs(solution.residual_row1)) < 2.0e-15
    assert np.max(np.abs(solution.residual_row2)) < 2.0e-15
    assert np.all(np.isfinite(solution.c1))
    assert np.all(np.isfinite(solution.c2))
    assert np.all(np.isfinite(c1_eta))
    assert np.all(np.isfinite(c2_eta))
    assert candidate.truth_boundary["current_cartesian_end_compensation_composed"] is False
    assert candidate.truth_boundary["source_exact_amplitude_root_materialized"] is False


def test_configuration_roundtrip_and_truth_mutations_fail_closed(tmp_path) -> None:
    candidate = KokunoCurrentPulseEntryMoments()
    path = tmp_path / "current-pulse-entry-moments.json"
    config = candidate.save_configuration(path)
    loaded = KokunoCurrentPulseEntryMoments.load_configuration(path)

    assert config["parent_exact_head"] == PARENT_EXACT_HEAD
    assert loaded.semantic_sha256 == candidate.semantic_sha256
    np.testing.assert_allclose(
        loaded.entry_ratios(np.asarray([0.17, 0.58]))["j_entry_ratio"],
        candidate.entry_ratios(np.asarray([0.17, 0.58]))["j_entry_ratio"],
        rtol=0.0,
        atol=0.0,
    )

    mutated = copy.deepcopy(config)
    mutated["binding"]["J_policy"] = "assume_zero"
    with pytest.raises(ValueError):
        KokunoCurrentPulseEntryMoments.from_configuration(mutated)

    mutated = copy.deepcopy(config)
    mutated["truth_boundary"]["current_cartesian_end_compensation_composed"] = True
    with pytest.raises(ValueError):
        KokunoCurrentPulseEntryMoments.from_configuration(mutated)

    mutated = copy.deepcopy(config)
    mutated["source"]["commit"] = "paper-exact-hidden-field"
    with pytest.raises(ValueError):
        KokunoCurrentPulseEntryMoments.from_configuration(mutated)
