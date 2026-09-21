from __future__ import annotations

from dataclasses import replace
import inspect

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_five_moment_correction_gain import (
    MOMENT_DIMENSION,
    FiveMomentCorrectionSystem,
    MomentCorrectionGainError,
)
from openai_ns_reconstruction.kokuno_finite_correction_stage import (
    ValidationPartition,
)
from openai_ns_reconstruction.kokuno_gain_gated_finite_correction_stage import (
    FINAL_NORMALIZED_DIVERGENCE_GATE,
    FINAL_NORMALIZED_MOMENTUM_GATE,
    GainGatedCorrectionBinding,
    GainGatedFiniteCorrectionError,
    _GainGatedMechanicsBackend,
    _coefficient_sha256,
    build_mechanics_report,
    run_gain_gated_finite_correction_stage,
    truth_boundary,
)


def partitions() -> tuple[ValidationPartition, ValidationPartition]:
    return (
        ValidationPartition("held-in", ("train-0", "train-1", "train-2")),
        ValidationPartition("held-out", ("test-0", "test-1", "test-2")),
    )


def test_gain_gate_materializes_internal_correction_then_observes_half_contraction() -> None:
    held_in, held_out = partitions()
    result = run_gain_gated_finite_correction_stage(
        _GainGatedMechanicsBackend(),
        {"amplitude": 1.0},
        held_in,
        held_out,
    )
    assert result.candidate_after_if_accepted == {"amplitude": 0.5}
    receipt = result.receipt
    assert receipt.local_five_moment_gain_gate_passed
    assert receipt.local_moment_solve.local_moment_fixed_point_converged
    assert receipt.local_moment_solve.coefficient_vector[0] == pytest.approx(-0.5)
    assert receipt.correction_materialized_only_after_gain_gate
    assert receipt.correction_bound_to_exact_solve
    assert receipt.held_out_excluded_from_correction_fit
    assert receipt.fixed_physical_stage.physical_contract_identity_preserved
    assert receipt.fixed_physical_stage.stage_accepted
    assert receipt.stage_accepted
    assert receipt.mechanics_only
    assert not receipt.candidate_residual_evidence
    assert (
        receipt.fixed_physical_stage.parent_stage.held_out_momentum_contraction_factor
        == pytest.approx(0.5)
    )
    assert not receipt.analytic_full_ns_cycle_gain_available
    assert not receipt.observed_gain_is_analytic_full_ns_bound
    assert not receipt.pde_validated


def test_public_stage_exposes_no_caller_correction_gain_or_scientific_knobs() -> None:
    signature = inspect.signature(run_gain_gated_finite_correction_stage)
    assert tuple(signature.parameters) == (
        "backend",
        "candidate",
        "held_in",
        "held_out",
    )
    forbidden = {
        "correction",
        "residual",
        "defect",
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


class LeakageBackend(_GainGatedMechanicsBackend):
    materialize_called = False

    def five_moment_correction_system(
        self, candidate: dict[str, float]
    ) -> FiveMomentCorrectionSystem:
        system = super().five_moment_correction_system(candidate)
        return replace(system, fit_sample_ids=("train-0", "test-0"))

    def materialize_velocity_correction(self, candidate, solve_receipt):
        self.materialize_called = True
        return super().materialize_velocity_correction(candidate, solve_receipt)


def test_held_out_leakage_is_rejected_before_solve_materialization() -> None:
    held_in, held_out = partitions()
    backend = LeakageBackend()
    with pytest.raises(GainGatedFiniteCorrectionError, match="subset of held-in|leaked"):
        run_gain_gated_finite_correction_stage(
            backend, {"amplitude": 1.0}, held_in, held_out
        )
    assert not backend.materialize_called


class IncompleteDefectBackend(_GainGatedMechanicsBackend):
    materialize_called = False

    def five_moment_correction_system(
        self, candidate: dict[str, float]
    ) -> FiveMomentCorrectionSystem:
        return replace(
            super().five_moment_correction_system(candidate),
            discrepancy_from_complete_ns_defect=False,
        )

    def materialize_velocity_correction(self, candidate, solve_receipt):
        self.materialize_called = True
        return super().materialize_velocity_correction(candidate, solve_receipt)


def test_incomplete_defect_moment_system_is_rejected_before_materialization() -> None:
    held_in, held_out = partitions()
    backend = IncompleteDefectBackend()
    with pytest.raises(
        GainGatedFiniteCorrectionError, match="not derived from a complete NS defect"
    ):
        run_gain_gated_finite_correction_stage(
            backend, {"amplitude": 1.0}, held_in, held_out
        )
    assert not backend.materialize_called


class NonContractingBackend(_GainGatedMechanicsBackend):
    materialize_called = False

    def five_moment_correction_system(
        self, candidate: dict[str, float]
    ) -> FiveMomentCorrectionSystem:
        system = super().five_moment_correction_system(candidate)
        Q = np.zeros(
            (MOMENT_DIMENSION, MOMENT_DIMENSION, MOMENT_DIMENSION), dtype=float
        )
        for i in range(MOMENT_DIMENSION):
            Q[i, i, i] = 20.0
        return replace(
            system,
            bilinear_Q=tuple(
                tuple(tuple(float(v) for v in row) for row in plane)
                for plane in Q
            ),
        )

    def materialize_velocity_correction(self, candidate, solve_receipt):
        self.materialize_called = True
        return super().materialize_velocity_correction(candidate, solve_receipt)


def test_noncontracting_local_map_is_rejected_before_correction_materialization() -> None:
    held_in, held_out = partitions()
    backend = NonContractingBackend()
    with pytest.raises(MomentCorrectionGainError, match="rejected before iteration"):
        run_gain_gated_finite_correction_stage(
            backend, {"amplitude": 1.0}, held_in, held_out
        )
    assert not backend.materialize_called


class TamperedBindingBackend(_GainGatedMechanicsBackend):
    def gain_gated_correction_provenance(
        self, candidate: dict[str, float], correction: dict[str, object]
    ) -> GainGatedCorrectionBinding:
        good = super().gain_gated_correction_provenance(candidate, correction)
        bad_vector = tuple(
            value + (0.01 if index == 0 else 0.0)
            for index, value in enumerate(good.coefficient_vector)
        )
        return replace(
            good,
            coefficient_vector=bad_vector,
            coefficient_sha256=_coefficient_sha256(bad_vector),
        )


def test_materialized_correction_must_bind_exact_solved_coefficients() -> None:
    held_in, held_out = partitions()
    with pytest.raises(
        GainGatedFiniteCorrectionError, match="differs from solved coefficients"
    ):
        run_gain_gated_finite_correction_stage(
            TamperedBindingBackend(),
            {"amplitude": 1.0},
            held_in,
            held_out,
        )


def test_binding_dataclass_rejects_forged_coefficient_hash() -> None:
    coefficients = (-0.5, 0.0, 0.0, 0.0, 0.0)
    with pytest.raises(ValueError, match="coefficient_sha256"):
        GainGatedCorrectionBinding(
            correction_id="x",
            source_candidate_id="c",
            source_candidate_sha256="h",
            moment_system_sha256="m",
            coefficient_sha256="forged",
            coefficient_vector=coefficients,
            fit_partition_id="held-in",
            fit_sample_ids=("train-0",),
            derived_from_complete_ns_defect=True,
            residual_as_forcing_shortcut_used=False,
        )


def test_truth_boundary_keeps_local_gain_separate_from_full_ns_cycle() -> None:
    truth = truth_boundary()
    assert truth["local_five_moment_analytic_contraction_gate_integrated"]
    assert truth["correction_materialized_only_from_authenticated_solve"]
    assert truth["five_moment_fit_must_be_held_in_only"]
    assert truth["held_out_fit_leakage_rejected_before_solve"]
    assert truth["complete_ns_defect_required_before_solve"]
    assert truth["fixed_physical_contract_stage_required_after_materialization"]
    assert truth["forbidden_public_parameters_absent"]
    assert not truth["analytic_full_ns_cycle_gain_available"]
    assert not truth["current_real_five_moment_discrepancy_materialized"]
    assert not truth["current_real_ns_correction_velocity_materialized"]
    assert not truth["real_candidate_finite_correction_cycle_run"]
    assert not truth["heldout_normalized_ns_residual_assessed"]
    assert not truth["residual_reduction_claimed"]
    assert not truth["same_protocol_comparable_to_st006"]
    assert not truth["pde_validated"]
    assert truth["final_normalized_momentum_gate"] == FINAL_NORMALIZED_MOMENTUM_GATE
    assert truth["final_normalized_divergence_gate"] == FINAL_NORMALIZED_DIVERGENCE_GATE


def test_mechanics_report_is_explicitly_not_candidate_residual_evidence() -> None:
    report = build_mechanics_report()
    assert report["mechanics_only"]
    assert not report["candidate_residual_evidence"]
    assert report["mechanics_stage_receipt"]["stage_accepted"]
    assert report["observed_held_out_momentum_contraction_factor"] == pytest.approx(0.5)
    assert not report["truth_boundary"]["pde_validated"]
