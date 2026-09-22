from __future__ import annotations

from dataclasses import replace
import inspect
import math

import pytest

import openai_ns_reconstruction.kokuno_rf49_source_exponent_gain_firewall as rf49
from openai_ns_reconstruction.kokuno_rf44_rf49_postupdate_recompute import (
    materialize_rf44_rf49_postupdate_recompute as parent_recompute,
)
from openai_ns_reconstruction.kokuno_rf49_source_exponent_gain_firewall import (
    AUTONOMOUS_DEFECT_NORM_NONEXPANSION_LIMIT,
    AUTONOMOUS_RATIO_ROUNDOFF,
    FINAL_NORMALIZED_DIVERGENCE_GATE,
    FINAL_NORMALIZED_MOMENTUM_GATE,
    RF49SourceGainError,
    SOURCE_KAPPA_S,
    SOURCE_REQUIRED_EXPONENT_GAIN,
    _RF49MechanicsBackend,
    assess_rf49_source_exponent_gain_firewall,
    truth_boundary,
)


def test_source_exponent_arithmetic_is_exactly_separated_from_norm_ratio() -> None:
    r = assess_rf49_source_exponent_gain_firewall(_RF49MechanicsBackend(), None)
    a = r.source_arithmetic
    assert a.stage_index == 0
    assert a.kappa_s == pytest.approx(SOURCE_KAPPA_S, abs=0.0)
    assert a.sigma_j == pytest.approx(0.2)
    assert a.B_j == pytest.approx(0.7)
    assert a.C_j == pytest.approx(1.2)
    assert a.wave_residual_exponent_gain == pytest.approx(0.39999)
    assert a.tangential_mean_raw_margin_018_minus_2k == pytest.approx(0.17998)
    assert a.tangential_mean_raw_margin_sigma_minus_3k == pytest.approx(0.19997)
    assert a.tangential_mean_raw_margin_1_minus_4k == pytest.approx(0.99996)
    assert a.retained_tangential_mean_exponent_gain == pytest.approx(0.17)
    assert a.defect_exponent_gain == pytest.approx(0.89996)
    assert a.source_required_exponent_gain == pytest.approx(0.1)
    assert a.source_exponent_gain_certified is True
    assert min(
        a.wave_residual_exponent_gain,
        a.retained_tangential_mean_exponent_gain,
        a.defect_exponent_gain,
    ) > SOURCE_REQUIRED_EXPONENT_GAIN
    # The observed norm ratio is a separate numerical diagnostic, not one of
    # the source exponent gains.
    assert math.isfinite(r.observed_rf44_defect_norm_ratio_after_over_before)


def test_mechanics_handoff_records_nonexpansion_result_without_scientific_promotion() -> None:
    r = assess_rf49_source_exponent_gain_firewall(_RF49MechanicsBackend(), None)
    expected = r.observed_rf44_defect_norm_ratio_after_over_before <= (
        AUTONOMOUS_DEFECT_NORM_NONEXPANSION_LIMIT + AUTONOMOUS_RATIO_ROUNDOFF
    )
    assert r.observed_defect_nonexpansion is expected
    assert r.autonomous_divergence_guard_passed is expected
    assert r.source_specific_mean_handoff_accepted is expected
    assert r.repository_candidate_gain_evidence is False
    assert r.current_i4_candidate_gain_evidence is False
    assert r.complete_ns_cycle_gain is False
    assert r.heldout_ns_residual_assessed is False
    assert r.pde_validated is False


def test_observed_defect_growth_is_retained_as_rejected_receipt(monkeypatch) -> None:
    backend = _RF49MechanicsBackend()
    baseline = parent_recompute(backend, None)
    inflated = replace(baseline, defect_gain_ratio_2=1.25)
    monkeypatch.setattr(
        rf49,
        "materialize_rf44_rf49_postupdate_recompute",
        lambda backend, candidate: inflated,
    )
    r = rf49.assess_rf49_source_exponent_gain_firewall(backend, None)
    assert r.source_arithmetic.source_exponent_gain_certified is True
    assert r.observed_rf44_defect_norm_ratio_after_over_before == pytest.approx(1.25)
    assert r.autonomous_divergence_guard_passed is False
    assert r.source_specific_mean_handoff_accepted is False
    assert r.rejection_reasons == ("observed_rf44_defect_norm_increased",)


def test_heldout_selected_gain_context_fails_closed() -> None:
    class Bad(_RF49MechanicsBackend):
        def rf49_source_gain_context(self, candidate, recompute):
            return replace(
                super().rf49_source_gain_context(candidate, recompute),
                heldout_samples_used_to_select_stage_or_gain=True,
            )

    with pytest.raises(RF49SourceGainError, match="held-out"):
        assess_rf49_source_exponent_gain_firewall(Bad(), None)


def test_caller_supplied_gain_context_fails_closed() -> None:
    class Bad(_RF49MechanicsBackend):
        def rf49_source_gain_context(self, candidate, recompute):
            return replace(
                super().rf49_source_gain_context(candidate, recompute),
                caller_supplied_gain_used=True,
            )

    with pytest.raises(RF49SourceGainError, match="caller supplied"):
        assess_rf49_source_exponent_gain_firewall(Bad(), None)


def test_identity_drift_fails_closed() -> None:
    class Bad(_RF49MechanicsBackend):
        def rf49_source_gain_context(self, candidate, recompute):
            return replace(
                super().rf49_source_gain_context(candidate, recompute),
                source_chart_sha256="0" * 64,
            )

    with pytest.raises(RF49SourceGainError, match="identity mismatch"):
        assess_rf49_source_exponent_gain_firewall(Bad(), None)


def test_stage_must_be_frozen_and_nonnegative() -> None:
    class NotFrozen(_RF49MechanicsBackend):
        def rf49_source_gain_context(self, candidate, recompute):
            return replace(
                super().rf49_source_gain_context(candidate, recompute),
                stage_index_frozen_before_defect_evaluation=False,
            )

    with pytest.raises(RF49SourceGainError, match="applicability provenance"):
        assess_rf49_source_exponent_gain_firewall(NotFrozen(), None)

    class Negative(_RF49MechanicsBackend):
        def rf49_source_gain_context(self, candidate, recompute):
            return replace(
                super().rf49_source_gain_context(candidate, recompute),
                stage_index=-1,
            )

    with pytest.raises(RF49SourceGainError, match="nonnegative integer"):
        assess_rf49_source_exponent_gain_firewall(Negative(), None)


def test_public_entry_point_has_no_gain_residual_or_threshold_knobs() -> None:
    assert tuple(inspect.signature(assess_rf49_source_exponent_gain_firewall).parameters) == (
        "backend",
        "candidate",
    )


def test_truth_boundary_keeps_source_gain_and_scientific_acceptance_distinct() -> None:
    t = truth_boundary()
    assert t["source_exponent_gain_firewall_materialized"] is True
    assert t["source_exponent_gain_is_multiplicative_norm_contraction"] is False
    assert t["observed_rf44_defect_norm_nonexpansion_guard_materialized"] is True
    assert t["observed_nonexpansion_guard_role"] == "repository-autonomous engineering rejection only"
    assert t["new_project_scientific_acceptance_threshold_introduced"] is False
    assert t["source_required_exponent_gain"] == pytest.approx(0.1)
    for key in (
        "current_i4_source_chart_backend_materialized",
        "current_i4_candidate_gain_evidence",
        "cartesian_correction_velocity_materialized",
        "complete_ns_defect",
        "complete_ns_cycle_gain",
        "real_candidate_finite_correction_cycle_run",
        "heldout_normalized_ns_residual_assessed",
        "residual_reduction_claimed",
        "same_protocol_comparable_to_st006",
        "pde_validated",
    ):
        assert t[key] is False, key
    assert t["final_normalized_momentum_gate"] == FINAL_NORMALIZED_MOMENTUM_GATE == 1e-3
    assert t["final_normalized_divergence_gate"] == FINAL_NORMALIZED_DIVERGENCE_GATE == 1e-5
    assert t["residual_defined_free_forcing_allowed"] is False
    assert t["finite_stage_small_residual_is_blowup_proof"] is False
