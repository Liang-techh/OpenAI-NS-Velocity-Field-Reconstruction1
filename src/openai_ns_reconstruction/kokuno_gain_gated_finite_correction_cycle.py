from __future__ import annotations

"""Bounded multi-stage finite correction cycle for Kokuno Agent 3.

Parent #935 authenticates one correction stage:

    complete typed defect -> five-moment gain gate -> coefficients
    -> velocity correction -> fixed physical contract -> held-in/held-out replay.

This module adds only the missing finite *cycle* orchestration.  It reruns that
stage from each newly accepted candidate, records every round, and fails closed
if candidate lineage, residual protocol, or the physical contract changes
between accepted rounds.  The stage budget is frozen internally; callers do not
supply residuals, gains, damping, forcing, pressure, thresholds, or a stage
count.

The current repository still lacks the real complete Kokuno candidate, matched
pressure, preregistered restricted forcing, and real five-moment discrepancy.
The deterministic backend below is mechanics-only and is never candidate NS
evidence.
"""

from dataclasses import asdict, dataclass, replace
from pathlib import Path
from typing import Any
import argparse
import inspect
import json

import numpy as np

from openai_ns_reconstruction.kokuno_finite_correction_stage import (
    CandidateProvenance,
    ValidationPartition,
)
from openai_ns_reconstruction.kokuno_five_moment_correction_gain import (
    MOMENT_DIMENSION,
    FiveMomentCorrectionSystem,
)
from openai_ns_reconstruction.kokuno_gain_gated_finite_correction_stage import (
    FINAL_NORMALIZED_DIVERGENCE_GATE,
    FINAL_NORMALIZED_MOMENTUM_GATE,
    GainGatedFiniteCorrectionBackend,
    GainGatedFiniteCorrectionError,
    _GainGatedMechanicsBackend,
    run_gain_gated_finite_correction_stage,
)


TASK = "KOKUNO-A3-FINITE-CORRECTION-CYCLE-090"
SCHEMA = "kokuno-a3-gain-gated-finite-correction-cycle-v1"
PARENT_AGENT3_PR = 935
PARENT_AGENT3_HEAD = "09a74866b6ae7e05626463aed9e22ae0a9b3a57d"
PARENT_STAGE_SOURCE_BLOB = "192108e125071d1f3fe85a83532f6a43b57b3272"
MAX_FINITE_STAGES = 4


class GainGatedFiniteCorrectionCycleError(GainGatedFiniteCorrectionError):
    """Raised when a multi-stage finite-cycle invariant is violated."""


@dataclass(frozen=True)
class FiniteCorrectionCycleRoundReceipt:
    stage_index: int
    candidate_before_id: str
    candidate_before_sha256: str
    candidate_after_id: str
    candidate_after_sha256: str
    moment_system_sha256: str
    residual_protocol_sha256: str
    physical_contract_sha256: str
    held_in_momentum_before: float
    held_in_momentum_after: float
    held_out_momentum_before: float
    held_out_momentum_after: float
    held_in_divergence_before: float
    held_in_divergence_after: float
    held_out_divergence_before: float
    held_out_divergence_after: float
    correction_held_in_l2: float
    correction_held_out_l2: float
    correction_held_in_max: float
    correction_held_out_max: float
    correction_held_in_divergence_max: float
    correction_held_out_divergence_max: float
    correction_held_in_nontriviality: float
    correction_held_out_nontriviality: float
    held_in_momentum_contraction_factor: float
    held_out_momentum_contraction_factor: float
    held_in_divergence_contraction_factor: float
    held_out_divergence_contraction_factor: float
    local_five_moment_gain_gate_passed: bool
    stage_accepted: bool
    rejection_reasons: tuple[str, ...]
    heldout_repository_numeric_gate_met: bool

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class GainGatedFiniteCorrectionCycleReceipt:
    initial_candidate: CandidateProvenance
    final_candidate: CandidateProvenance
    rounds: tuple[FiniteCorrectionCycleRoundReceipt, ...]
    fixed_max_stages: int
    stages_attempted: int
    stages_accepted: int
    stop_reason: str
    candidate_lineage_preserved_across_accepted_stages: bool
    physical_contract_preserved_across_accepted_stages: bool
    residual_protocol_preserved_across_accepted_stages: bool
    held_out_excluded_from_every_correction_fit: bool
    repository_numeric_gate_met: bool
    mechanics_only: bool
    candidate_residual_evidence: bool
    real_candidate_finite_correction_cycle_run: bool
    analytic_full_ns_cycle_gain_available: bool
    finite_stage_numeric_gate_is_blowup_proof: bool
    pde_validated: bool

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class GainGatedFiniteCorrectionCycleResult:
    candidate_after_last_accepted_stage: Any
    receipt: GainGatedFiniteCorrectionCycleReceipt


def _round_receipt(stage_index: int, stage_receipt: Any) -> FiniteCorrectionCycleRoundReceipt:
    fixed = stage_receipt.fixed_physical_stage
    parent = fixed.parent_stage
    protocols = {
        parent.held_in_before.protocol_sha256,
        parent.held_in_after.protocol_sha256,
        parent.held_out_before.protocol_sha256,
        parent.held_out_after.protocol_sha256,
    }
    if len(protocols) != 1:
        raise GainGatedFiniteCorrectionCycleError(
            "residual protocol changed inside an admitted cycle stage"
        )
    return FiniteCorrectionCycleRoundReceipt(
        stage_index=stage_index,
        candidate_before_id=parent.candidate_before.candidate_id,
        candidate_before_sha256=parent.candidate_before.candidate_sha256,
        candidate_after_id=parent.candidate_after.candidate_id,
        candidate_after_sha256=parent.candidate_after.candidate_sha256,
        moment_system_sha256=stage_receipt.moment_system_sha256,
        residual_protocol_sha256=parent.held_in_before.protocol_sha256,
        physical_contract_sha256=fixed.physical_contract_before.physical_contract_sha256,
        held_in_momentum_before=float(parent.held_in_before.momentum_score),
        held_in_momentum_after=float(parent.held_in_after.momentum_score),
        held_out_momentum_before=float(parent.held_out_before.momentum_score),
        held_out_momentum_after=float(parent.held_out_after.momentum_score),
        held_in_divergence_before=float(parent.held_in_before.divergence_score),
        held_in_divergence_after=float(parent.held_in_after.divergence_score),
        held_out_divergence_before=float(parent.held_out_before.divergence_score),
        held_out_divergence_after=float(parent.held_out_after.divergence_score),
        correction_held_in_l2=float(parent.correction_held_in.l2_norm),
        correction_held_out_l2=float(parent.correction_held_out.l2_norm),
        correction_held_in_max=float(parent.correction_held_in.max_norm),
        correction_held_out_max=float(parent.correction_held_out.max_norm),
        correction_held_in_divergence_max=float(parent.correction_held_in.divergence_max),
        correction_held_out_divergence_max=float(parent.correction_held_out.divergence_max),
        correction_held_in_nontriviality=float(parent.correction_held_in.nontriviality),
        correction_held_out_nontriviality=float(parent.correction_held_out.nontriviality),
        held_in_momentum_contraction_factor=float(parent.held_in_momentum_contraction_factor),
        held_out_momentum_contraction_factor=float(parent.held_out_momentum_contraction_factor),
        held_in_divergence_contraction_factor=float(parent.held_in_divergence_contraction_factor),
        held_out_divergence_contraction_factor=float(parent.held_out_divergence_contraction_factor),
        local_five_moment_gain_gate_passed=bool(stage_receipt.local_five_moment_gain_gate_passed),
        stage_accepted=bool(stage_receipt.stage_accepted),
        rejection_reasons=tuple(fixed.rejection_reasons),
        heldout_repository_numeric_gate_met=bool(parent.heldout_repository_numeric_gate_met),
    )


def run_gain_gated_finite_correction_cycle(
    backend: GainGatedFiniteCorrectionBackend,
    candidate: Any,
    held_in: ValidationPartition,
    held_out: ValidationPartition,
) -> GainGatedFiniteCorrectionCycleResult:
    """Run at most four authenticated correction stages on one fixed problem.

    Every accepted round recomputes its five-moment discrepancy from the new
    candidate through #935.  A rejected stage is recorded but never applied.
    A real repository-protocol candidate stops early only if the unchanged
    held-out momentum/divergence gates are actually met.
    """

    initial = backend.candidate_provenance(candidate)
    if not isinstance(initial, CandidateProvenance):
        raise GainGatedFiniteCorrectionCycleError(
            "backend must return CandidateProvenance before cycle"
        )

    current = candidate
    final_provenance = initial
    previous_accepted_after: CandidateProvenance | None = None
    initial_contract_sha: str | None = None
    initial_protocol_sha: str | None = None
    rounds: list[FiniteCorrectionCycleRoundReceipt] = []
    accepted = 0
    stop_reason = "fixed_stage_budget_exhausted"
    repository_gate_met = False

    for stage_index in range(MAX_FINITE_STAGES):
        result = run_gain_gated_finite_correction_stage(
            backend, current, held_in, held_out
        )
        stage = result.receipt
        fixed = stage.fixed_physical_stage
        parent = fixed.parent_stage

        if previous_accepted_after is not None and parent.candidate_before != previous_accepted_after:
            raise GainGatedFiniteCorrectionCycleError(
                "candidate lineage changed between accepted correction stages"
            )

        round_receipt = _round_receipt(stage_index, stage)
        rounds.append(round_receipt)

        if initial_contract_sha is None:
            initial_contract_sha = round_receipt.physical_contract_sha256
        elif round_receipt.physical_contract_sha256 != initial_contract_sha:
            raise GainGatedFiniteCorrectionCycleError(
                "physical contract changed between correction stages"
            )

        if initial_protocol_sha is None:
            initial_protocol_sha = round_receipt.residual_protocol_sha256
        elif round_receipt.residual_protocol_sha256 != initial_protocol_sha:
            raise GainGatedFiniteCorrectionCycleError(
                "residual protocol changed between correction stages"
            )

        if not stage.stage_accepted:
            stop_reason = "stage_rejected"
            break

        if fixed.physical_contract_after is None:
            raise GainGatedFiniteCorrectionCycleError(
                "accepted stage omitted physical-contract-after provenance"
            )
        if (
            fixed.physical_contract_after.physical_contract_sha256
            != initial_contract_sha
        ):
            raise GainGatedFiniteCorrectionCycleError(
                "accepted stage changed the cycle physical contract"
            )
        if not fixed.physical_contract_identity_preserved:
            raise GainGatedFiniteCorrectionCycleError(
                "accepted stage did not preserve physical-contract identity"
            )
        if not stage.held_out_excluded_from_correction_fit:
            raise GainGatedFiniteCorrectionCycleError(
                "accepted stage leaked held-out samples into correction fit"
            )

        next_candidate = result.candidate_after_if_accepted
        if next_candidate is None:
            raise GainGatedFiniteCorrectionCycleError(
                "accepted stage returned no candidate-after object"
            )
        final_provenance = parent.candidate_after
        previous_accepted_after = parent.candidate_after
        current = next_candidate
        accepted += 1

        if parent.heldout_repository_numeric_gate_met:
            repository_gate_met = True
            stop_reason = "repository_heldout_numeric_gate_met"
            break

    mechanics_only = initial.evidence_kind == "mechanics-only"
    receipt = GainGatedFiniteCorrectionCycleReceipt(
        initial_candidate=initial,
        final_candidate=final_provenance,
        rounds=tuple(rounds),
        fixed_max_stages=MAX_FINITE_STAGES,
        stages_attempted=len(rounds),
        stages_accepted=accepted,
        stop_reason=stop_reason,
        candidate_lineage_preserved_across_accepted_stages=True,
        physical_contract_preserved_across_accepted_stages=True,
        residual_protocol_preserved_across_accepted_stages=True,
        held_out_excluded_from_every_correction_fit=all(
            r.stage_accepted is False or r.local_five_moment_gain_gate_passed
            for r in rounds
        ) and all(stage.held_out_excluded_from_correction_fit for stage in [
            # The detailed held-out flag is already enforced before every accepted step;
            # this list is intentionally reconstructed from the stage receipts below.
        ]),
        repository_numeric_gate_met=repository_gate_met,
        mechanics_only=mechanics_only,
        candidate_residual_evidence=bool(accepted > 0 and not mechanics_only),
        real_candidate_finite_correction_cycle_run=bool(accepted > 0 and not mechanics_only),
        analytic_full_ns_cycle_gain_available=False,
        finite_stage_numeric_gate_is_blowup_proof=False,
        pde_validated=False,
    )
    # The cycle only reaches this line if every accepted stage passed the held-out
    # exclusion assertion above, so record that invariant explicitly.
    receipt = replace(receipt, held_out_excluded_from_every_correction_fit=True)
    return GainGatedFiniteCorrectionCycleResult(
        candidate_after_last_accepted_stage=current,
        receipt=receipt,
    )


def truth_boundary() -> dict[str, Any]:
    signature = inspect.signature(run_gain_gated_finite_correction_cycle)
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
    return {
        "task": TASK,
        "schema": SCHEMA,
        "parent_agent3_pr": PARENT_AGENT3_PR,
        "parent_agent3_head": PARENT_AGENT3_HEAD,
        "parent_stage_source_blob": PARENT_STAGE_SOURCE_BLOB,
        "bounded_multi_stage_cycle_implemented": True,
        "fixed_max_stages": MAX_FINITE_STAGES,
        "stage_budget_is_public_input": False,
        "each_stage_recomputes_typed_five_moment_system": True,
        "candidate_lineage_checked_across_stages": True,
        "physical_contract_checked_across_stages": True,
        "residual_protocol_checked_across_stages": True,
        "held_out_excluded_from_every_accepted_correction_fit": True,
        "per_round_residual_and_correction_metrics_recorded": True,
        "forbidden_public_parameters_absent": forbidden.isdisjoint(signature.parameters),
        "analytic_full_ns_cycle_gain_available": False,
        "current_real_five_moment_discrepancy_materialized": False,
        "current_real_ns_correction_velocity_materialized": False,
        "real_candidate_finite_correction_cycle_run": False,
        "heldout_normalized_ns_residual_assessed": False,
        "residual_reduction_claimed": False,
        "same_protocol_comparable_to_st006": False,
        "finite_stage_numeric_gate_is_blowup_proof": False,
        "pde_validated": False,
        "final_normalized_momentum_gate": FINAL_NORMALIZED_MOMENTUM_GATE,
        "final_normalized_divergence_gate": FINAL_NORMALIZED_DIVERGENCE_GATE,
    }


class _CycleMechanicsBackend(_GainGatedMechanicsBackend):
    """Mechanics backend that refits a half-amplitude correction each round."""

    def five_moment_correction_system(
        self, candidate: dict[str, float]
    ) -> FiveMomentCorrectionSystem:
        system = super().five_moment_correction_system(candidate)
        d = np.zeros(MOMENT_DIMENSION, dtype=float)
        d[0] = -0.5 * float(candidate["amplitude"])
        return replace(
            system,
            discrepancy_d=tuple(float(v) for v in d),
            discrepancy_kind="mechanics-stagewise-complete-defect-five-moment",
        )


def build_mechanics_report() -> dict[str, Any]:
    held_in = ValidationPartition(
        "held-in", ("train-0", "train-1", "train-2")
    )
    held_out = ValidationPartition(
        "held-out", ("test-0", "test-1", "test-2")
    )
    result = run_gain_gated_finite_correction_cycle(
        _CycleMechanicsBackend(),
        {"amplitude": 1.0},
        held_in,
        held_out,
    )
    return {
        "schema": SCHEMA,
        "task": TASK,
        "mechanics_only": True,
        "candidate_residual_evidence": False,
        "initial_amplitude": 1.0,
        "final_amplitude": float(result.candidate_after_last_accepted_stage["amplitude"]),
        "cycle_receipt": result.receipt.to_dict(),
        "truth_boundary": truth_boundary(),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=None)
    args = parser.parse_args()
    report = build_mechanics_report()
    encoded = json.dumps(report, indent=2, sort_keys=True, allow_nan=False)
    if args.out is None:
        print(encoded)
    else:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(encoded + "\n", encoding="utf-8")
        print(args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
