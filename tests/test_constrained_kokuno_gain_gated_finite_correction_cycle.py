from __future__ import annotations

from dataclasses import replace
import inspect

import pytest

from openai_ns_reconstruction.kokuno_finite_correction_stage import ValidationPartition
from openai_ns_reconstruction.kokuno_gain_gated_finite_correction_stage import (
    GainGatedFiniteCorrectionError,
)
from openai_ns_reconstruction.kokuno_gain_gated_finite_correction_cycle import (
    FINAL_NORMALIZED_DIVERGENCE_GATE,
    FINAL_NORMALIZED_MOMENTUM_GATE,
    MAX_FINITE_STAGES,
    _CycleMechanicsBackend,
    build_mechanics_report,
    run_gain_gated_finite_correction_cycle,
    truth_boundary,
)


def partitions() -> tuple[ValidationPartition, ValidationPartition]:
    return (
        ValidationPartition("held-in", ("train-0", "train-1", "train-2")),
        ValidationPartition("held-out", ("test-0", "test-1", "test-2")),
    )


def test_cycle_refits_and_accepts_four_stagewise_half_corrections() -> None:
    held_in, held_out = partitions()
    result = run_gain_gated_finite_correction_cycle(
        _CycleMechanicsBackend(),
        {"amplitude": 1.0},
        held_in,
        held_out,
    )
    assert result.candidate_after_last_accepted_stage == {"amplitude": 0.0625}
    receipt = result.receipt
    assert receipt.fixed_max_stages == MAX_FINITE_STAGES == 4
    assert receipt.stages_attempted == 4
    assert receipt.stages_accepted == 4
    assert receipt.stop_reason == "fixed_stage_budget_exhausted"
    assert receipt.candidate_lineage_preserved_across_accepted_stages
    assert receipt.physical_contract_preserved_across_accepted_stages
    assert receipt.residual_protocol_preserved_across_accepted_stages
    assert receipt.held_out_excluded_from_every_correction_fit
    assert receipt.mechanics_only
    assert not receipt.candidate_residual_evidence
    assert not receipt.real_candidate_finite_correction_cycle_run
    assert not receipt.repository_numeric_gate_met
    assert not receipt.finite_stage_numeric_gate_is_blowup_proof
    assert not receipt.pde_validated

    expected_l2 = [0.5, 0.25, 0.125, 0.0625]
    for index, (round_receipt, expected) in enumerate(zip(receipt.rounds, expected_l2, strict=True)):
        assert round_receipt.stage_index == index
        assert round_receipt.stage_accepted
        assert round_receipt.local_five_moment_gain_gate_passed
        assert round_receipt.correction_held_in_l2 == pytest.approx(expected)
        assert round_receipt.correction_held_out_l2 == pytest.approx(expected)
        assert round_receipt.correction_held_in_nontriviality == pytest.approx(expected)
        assert round_receipt.correction_held_out_nontriviality == pytest.approx(expected)
        assert round_receipt.held_in_momentum_contraction_factor == pytest.approx(0.5)
        assert round_receipt.held_out_momentum_contraction_factor == pytest.approx(0.5)
        assert round_receipt.held_in_divergence_contraction_factor == pytest.approx(0.5)
        assert round_receipt.held_out_divergence_contraction_factor == pytest.approx(0.5)
        assert not round_receipt.heldout_repository_numeric_gate_met

    for left, right in zip(receipt.rounds[:-1], receipt.rounds[1:], strict=True):
        assert left.candidate_after_id == right.candidate_before_id
        assert left.candidate_after_sha256 == right.candidate_before_sha256
        assert left.physical_contract_sha256 == right.physical_contract_sha256
        assert left.residual_protocol_sha256 == right.residual_protocol_sha256


def test_public_cycle_has_no_stage_budget_residual_gain_or_physical_knobs() -> None:
    signature = inspect.signature(run_gain_gated_finite_correction_cycle)
    assert tuple(signature.parameters) == ("backend", "candidate", "held_in", "held_out")
    forbidden = {
        "correction",
        "residual",
        "defect",
        "max_stages",
        "stage_budget",
        "B",
        "Q",
        "d",
        "nu",
        "q",
        "gain",
        "damping",
        "alpha",
        "pressure",
        "forcing",
        "viscosity",
        "threshold",
        "scientific_threshold",
    }
    assert forbidden.isdisjoint(signature.parameters)


class RejectSecondStageBackend(_CycleMechanicsBackend):
    def five_moment_correction_system(self, candidate):
        system = super().five_moment_correction_system(candidate)
        if float(candidate["amplitude"]) <= 0.5:
            d = list(system.discrepancy_d)
            d[0] = 0.1
            return replace(system, discrepancy_d=tuple(d))
        return system


def test_rejected_second_stage_is_recorded_but_not_applied() -> None:
    held_in, held_out = partitions()
    result = run_gain_gated_finite_correction_cycle(
        RejectSecondStageBackend(),
        {"amplitude": 1.0},
        held_in,
        held_out,
    )
    receipt = result.receipt
    assert receipt.stages_attempted == 2
    assert receipt.stages_accepted == 1
    assert receipt.stop_reason == "stage_rejected"
    assert receipt.rounds[0].stage_accepted
    assert not receipt.rounds[1].stage_accepted
    assert "parent:held_in_momentum_increased" in receipt.rounds[1].rejection_reasons
    assert "parent:held_out_momentum_increased" in receipt.rounds[1].rejection_reasons
    assert result.candidate_after_last_accepted_stage == {"amplitude": 0.5}


class LeakOnSecondStageBackend(_CycleMechanicsBackend):
    def five_moment_correction_system(self, candidate):
        system = super().five_moment_correction_system(candidate)
        if float(candidate["amplitude"]) <= 0.5:
            return replace(system, fit_sample_ids=("train-0", "test-0"))
        return system


def test_later_stage_cannot_bypass_held_out_fit_isolation() -> None:
    held_in, held_out = partitions()
    with pytest.raises(GainGatedFiniteCorrectionError, match="subset of held-in|leaked"):
        run_gain_gated_finite_correction_cycle(
            LeakOnSecondStageBackend(),
            {"amplitude": 1.0},
            held_in,
            held_out,
        )


def test_truth_boundary_keeps_real_cycle_and_pde_claims_closed() -> None:
    truth = truth_boundary()
    assert truth["bounded_multi_stage_cycle_implemented"]
    assert truth["fixed_max_stages"] == 4
    assert not truth["stage_budget_is_public_input"]
    assert truth["each_stage_recomputes_typed_five_moment_system"]
    assert truth["candidate_lineage_checked_across_stages"]
    assert truth["physical_contract_checked_across_stages"]
    assert truth["residual_protocol_checked_across_stages"]
    assert truth["held_out_excluded_from_every_accepted_correction_fit"]
    assert truth["per_round_residual_and_correction_metrics_recorded"]
    assert truth["forbidden_public_parameters_absent"]
    assert not truth["analytic_full_ns_cycle_gain_available"]
    assert not truth["current_real_five_moment_discrepancy_materialized"]
    assert not truth["current_real_ns_correction_velocity_materialized"]
    assert not truth["real_candidate_finite_correction_cycle_run"]
    assert not truth["heldout_normalized_ns_residual_assessed"]
    assert not truth["residual_reduction_claimed"]
    assert not truth["same_protocol_comparable_to_st006"]
    assert not truth["finite_stage_numeric_gate_is_blowup_proof"]
    assert not truth["pde_validated"]
    assert truth["final_normalized_momentum_gate"] == FINAL_NORMALIZED_MOMENTUM_GATE
    assert truth["final_normalized_divergence_gate"] == FINAL_NORMALIZED_DIVERGENCE_GATE


def test_mechanics_report_is_not_candidate_residual_evidence() -> None:
    report = build_mechanics_report()
    assert report["mechanics_only"]
    assert not report["candidate_residual_evidence"]
    assert report["initial_amplitude"] == 1.0
    assert report["final_amplitude"] == pytest.approx(0.0625)
    assert report["cycle_receipt"]["stages_accepted"] == 4
    assert not report["cycle_receipt"]["real_candidate_finite_correction_cycle_run"]
    assert not report["truth_boundary"]["pde_validated"]
