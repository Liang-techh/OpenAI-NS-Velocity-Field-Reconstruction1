from __future__ import annotations

import inspect

import pytest

from openai_ns_reconstruction.kokuno_finite_correction_stage import (
    FINAL_NORMALIZED_DIVERGENCE_GATE,
    FINAL_NORMALIZED_MOMENTUM_GATE,
    KOKUNO_ANALYTIC_GAIN_AVAILABLE,
    KOKUNO_ANALYTIC_GAIN_BINDINGS_COMPLETE,
    OBSERVED_CONTRACTION_REJECT_AT,
    PARENT_AGENT3_HEAD,
    CandidateProvenance,
    CorrectionProvenance,
    FiniteCorrectionContractError,
    FullNSResidualMetrics,
    ValidationPartition,
    _MechanicsBackend,
    build_mechanics_report,
    current_agent3_scoped_source_cycle_admission,
    run_finite_correction_stage,
    truth_boundary,
)


def _partitions() -> tuple[ValidationPartition, ValidationPartition]:
    return (
        ValidationPartition("held-in", ("train-0", "train-1", "train-2")),
        ValidationPartition("held-out", ("test-0", "test-1", "test-2")),
    )


def test_mechanics_stage_recomputes_disjoint_before_after_and_contracts():
    held_in, held_out = _partitions()
    result = run_finite_correction_stage(
        _MechanicsBackend(),
        {"amplitude": 1.0},
        {"delta": -0.5},
        held_in,
        held_out,
    )
    receipt = result.receipt
    assert receipt.stage_accepted is True
    assert receipt.rejection_reasons == ()
    assert result.candidate_after_if_accepted == {"amplitude": 0.5}
    assert receipt.held_in_momentum_contraction_factor == pytest.approx(0.5)
    assert receipt.held_out_momentum_contraction_factor == pytest.approx(0.5)
    assert receipt.held_in_divergence_contraction_factor == pytest.approx(0.5)
    assert receipt.held_out_divergence_contraction_factor == pytest.approx(0.5)
    assert receipt.correction_held_in.l2_norm == pytest.approx(0.5)
    assert receipt.correction_held_out.nontriviality == pytest.approx(0.5)
    assert receipt.mechanics_only is True
    assert receipt.analytic_kokuno_gain_available is False
    assert receipt.observed_gain_is_analytic_kokuno_bound is False
    assert receipt.heldout_repository_numeric_gate_met is False
    assert receipt.pde_validated is False


def test_heldout_fit_leakage_is_rejected_before_any_success_claim():
    class LeakyBackend(_MechanicsBackend):
        def correction_provenance(self, candidate, correction):
            source = self.candidate_provenance(candidate)
            return CorrectionProvenance(
                correction_id="leaky",
                source_candidate_id=source.candidate_id,
                source_candidate_sha256=source.candidate_sha256,
                correction_operator="mechanics-leaky",
                fit_partition_id="held-in",
                fit_sample_ids=("train-0", "test-0"),
                derived_from_complete_ns_defect=True,
                residual_as_forcing_shortcut_used=False,
            )

    held_in, held_out = _partitions()
    with pytest.raises(FiniteCorrectionContractError, match="held-in samples"):
        run_finite_correction_stage(
            LeakyBackend(),
            {"amplitude": 1.0},
            {"delta": -0.5},
            held_in,
            held_out,
        )


def test_incomplete_candidate_and_residual_as_forcing_shortcut_fail_closed():
    class IncompleteBackend(_MechanicsBackend):
        def candidate_provenance(self, candidate):
            base = super().candidate_provenance(candidate)
            return CandidateProvenance(
                candidate_id=base.candidate_id,
                candidate_sha256=base.candidate_sha256,
                evidence_kind="real-candidate",
                complete_ns_defect=False,
                corrected_global_leading_join_complete=False,
                matched_pressure_gradient_included=False,
                restricted_forcing_included=False,
                restricted_forcing_preregistered=False,
                residual_as_forcing_shortcut_used=True,
                repository_acceptance_protocol=False,
            )

    held_in, held_out = _partitions()
    with pytest.raises(FiniteCorrectionContractError, match="not admitted"):
        run_finite_correction_stage(
            IncompleteBackend(),
            {"amplitude": 1.0},
            {"delta": -0.5},
            held_in,
            held_out,
        )


def test_observed_heldout_growth_rejects_correction_without_returning_candidate():
    class GrowingBackend(_MechanicsBackend):
        def apply_correction(self, candidate, correction):
            return {"amplitude": 1.5}

    held_in, held_out = _partitions()
    result = run_finite_correction_stage(
        GrowingBackend(),
        {"amplitude": 1.0},
        {"delta": 0.5},
        held_in,
        held_out,
    )
    assert result.receipt.stage_accepted is False
    assert result.candidate_after_if_accepted is None
    assert "held_in_momentum_increased" in result.receipt.rejection_reasons
    assert "held_out_momentum_increased" in result.receipt.rejection_reasons
    assert result.receipt.held_out_momentum_contraction_factor > OBSERVED_CONTRACTION_REJECT_AT


def test_residual_protocol_change_is_rejected():
    class ProtocolDriftBackend(_MechanicsBackend):
        calls = 0

        def evaluate_full_ns_residual(self, candidate, partition):
            metrics = super().evaluate_full_ns_residual(candidate, partition)
            self.calls += 1
            if self.calls < 3:
                return metrics
            return FullNSResidualMetrics(
                protocol_sha256="changed-after-correction",
                partition_id=metrics.partition_id,
                sample_count=metrics.sample_count,
                normalized_momentum_sample_max=metrics.normalized_momentum_sample_max,
                normalized_momentum_grid_l2=metrics.normalized_momentum_grid_l2,
                normalized_momentum_volume_l2=metrics.normalized_momentum_volume_l2,
                normalized_divergence_sample_max=metrics.normalized_divergence_sample_max,
                normalized_divergence_grid_l2=metrics.normalized_divergence_grid_l2,
                normalized_divergence_volume_l2=metrics.normalized_divergence_volume_l2,
                complete_ns_residual=True,
            )

    held_in, held_out = _partitions()
    with pytest.raises(FiniteCorrectionContractError, match="protocol changed"):
        run_finite_correction_stage(
            ProtocolDriftBackend(),
            {"amplitude": 1.0},
            {"delta": -0.5},
            held_in,
            held_out,
        )


def test_current_scoped_radial_stress_is_explicitly_not_admitted_to_cycle():
    admission = current_agent3_scoped_source_cycle_admission()
    assert admission["parent_agent3_head"] == PARENT_AGENT3_HEAD
    assert admission["admitted_to_real_finite_correction_cycle"] is False
    assert admission["surrogate_defect_substitution_allowed"] is False
    assert admission["residual_as_forcing_shortcut_allowed"] is False
    missing = set(admission["missing_requirements"])
    assert "complete_ns_defect" in missing
    assert "pressure_gradient_included" in missing
    assert "restricted_forcing_included" in missing
    assert "scoped_transport_stress_authorized_as_correction_target" in missing


def test_public_api_has_no_residual_gain_forcing_or_threshold_inputs():
    parameters = inspect.signature(run_finite_correction_stage).parameters
    assert list(parameters) == ["backend", "candidate", "correction", "held_in", "held_out"]
    forbidden = {
        "residual",
        "defect",
        "mean",
        "stress",
        "pressure",
        "forcing",
        "target",
        "gain",
        "rho0",
        "rho_star",
        "alpha",
        "damping",
        "viscosity",
        "nu",
        "delta_a",
        "normalized_score",
        "momentum_gate",
        "divergence_gate",
        "scientific_threshold",
    }
    assert forbidden.isdisjoint(parameters)


def test_truth_boundary_does_not_promote_mechanics_to_real_cycle_or_gain_bound():
    report = build_mechanics_report()
    boundary = truth_boundary()
    assert report["mechanics_only"] is True
    assert report["candidate_residual_evidence"] is False
    assert report["mechanics_stage_receipt"]["stage_accepted"] is True
    assert boundary["typed_full_ns_residual_contract_defined"] is True
    assert boundary["held_in_held_out_identity_disjointness_enforced"] is True
    assert boundary["correction_fit_must_be_subset_of_held_in"] is True
    assert boundary["held_out_fit_leakage_rejected"] is True
    assert boundary["residuals_recomputed_by_backend_before_and_after"] is True
    assert boundary["caller_supplied_residual_or_gain_allowed"] is False
    assert boundary["residual_as_forcing_shortcut_allowed"] is False
    assert boundary["analytic_kokuno_gain_bindings_complete"] is KOKUNO_ANALYTIC_GAIN_BINDINGS_COMPLETE is False
    assert boundary["analytic_kokuno_gain_available"] is KOKUNO_ANALYTIC_GAIN_AVAILABLE is False
    assert boundary["observed_gain_is_analytic_kokuno_bound"] is False
    assert boundary["current_scoped_transport_source_admitted"] is False
    assert boundary["real_candidate_finite_correction_cycle_run"] is False
    assert boundary["heldout_normalized_ns_residual_assessed"] is False
    assert boundary["residual_reduction_claimed"] is False
    assert boundary["same_protocol_comparable_to_st006"] is False
    assert boundary["pde_validated"] is False
    assert boundary["final_normalized_momentum_gate"] == FINAL_NORMALIZED_MOMENTUM_GATE == 1.0e-3
    assert boundary["final_normalized_divergence_gate"] == FINAL_NORMALIZED_DIVERGENCE_GATE == 1.0e-5
