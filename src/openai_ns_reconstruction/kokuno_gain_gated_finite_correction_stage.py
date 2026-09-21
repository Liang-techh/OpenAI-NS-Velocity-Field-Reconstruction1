from __future__ import annotations

"""Gain-gated finite-correction integration for Kokuno Agent 3.

This module connects the local five-moment contraction firewall from Agent-3
#927 to the fixed-physical-contract finite-correction stage from #895.

The public stage is intentionally narrow:

    run_gain_gated_finite_correction_stage(backend, candidate, held_in, held_out)

The caller cannot supply a correction, residual, defect, gain, damping factor,
pressure, forcing, viscosity, or scientific threshold.  The backend must expose
a typed five-moment system tied to the candidate and held-in samples.  This
module first solves that system through #927's source-derived contraction
certificate, then asks the backend to materialize a velocity correction from
the authenticated solve receipt, verifies that the correction is bound to the
exact system and coefficient vector, and only then evaluates it through the
existing fixed-physical-contract held-in/held-out stage.

This integration does not turn the local five-moment contraction bound into an
analytic full Navier--Stokes correction-cycle gain.  The current repository
still lacks a complete Kokuno global candidate / matched pressure / registered
forcing / real five-moment discrepancy handoff, so the real-candidate truth
flags remain false.
"""

from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Protocol
import argparse
import hashlib
import inspect
import json

import numpy as np

from openai_ns_reconstruction.kokuno_five_moment_correction_gain import (
    MOMENT_DIMENSION,
    CandidateIdentity,
    FiveMomentCorrectionSolveReceipt,
    FiveMomentCorrectionSystem,
    _canonical_system_sha256,
    solve_five_moment_correction,
)
from openai_ns_reconstruction.kokuno_finite_correction_stage import (
    CandidateProvenance,
    FiniteCorrectionContractError,
    ValidationPartition,
)
from openai_ns_reconstruction.kokuno_finite_correction_fixed_physical_contract import (
    FixedPhysicalContractStageReceipt,
    _FixedContractMechanicsBackend,
    run_fixed_physical_contract_finite_correction_stage,
)


TASK = "KOKUNO-A3-GAIN-GATED-FINITE-CORRECTION-STAGE-089"
SCHEMA = "kokuno-a3-gain-gated-finite-correction-stage-v1"
PARENT_AGENT3_PR = 927
PARENT_AGENT3_HEAD = "5029b6f88bfbbc4a498c3a9f9522f0ad2d79ecc1"
PARENT_GAIN_SOURCE_BLOB = "6625db71ec8478db5d6ec09cabb7033a03e658a4"
FIXED_CONTRACT_SOURCE_BLOB = "0cd41dfcc998064c37075b6086a246e93831ebfd"
SOURCE_READER_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_READER_HEAD = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_READER_TEX_BLOB = "205a99807302e21a51c5eaf223390c0dfc42bcd0"
FINAL_NORMALIZED_MOMENTUM_GATE = 1.0e-3
FINAL_NORMALIZED_DIVERGENCE_GATE = 1.0e-5


class GainGatedFiniteCorrectionError(FiniteCorrectionContractError):
    """Raised when local-gain provenance cannot authorize the finite stage."""


def _coefficient_sha256(values: tuple[float, ...]) -> str:
    payload = json.dumps(
        [float(v) for v in values],
        separators=(",", ":"),
        allow_nan=False,
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


@dataclass(frozen=True)
class GainGatedCorrectionBinding:
    """Identity binding from one exact five-moment solve to one correction."""

    correction_id: str
    source_candidate_id: str
    source_candidate_sha256: str
    moment_system_sha256: str
    coefficient_sha256: str
    coefficient_vector: tuple[float, ...]
    fit_partition_id: str
    fit_sample_ids: tuple[str, ...]
    derived_from_complete_ns_defect: bool
    residual_as_forcing_shortcut_used: bool

    def __post_init__(self) -> None:
        for value, name in (
            (self.correction_id, "correction_id"),
            (self.source_candidate_id, "source_candidate_id"),
            (self.source_candidate_sha256, "source_candidate_sha256"),
            (self.moment_system_sha256, "moment_system_sha256"),
            (self.coefficient_sha256, "coefficient_sha256"),
            (self.fit_partition_id, "fit_partition_id"),
        ):
            if not value:
                raise ValueError(f"{name} must be nonempty")
        if len(self.coefficient_vector) != MOMENT_DIMENSION:
            raise ValueError("coefficient_vector must contain exactly five values")
        if not np.all(np.isfinite(np.asarray(self.coefficient_vector, dtype=float))):
            raise ValueError("coefficient_vector must contain only finite values")
        if not self.fit_sample_ids or any(not item for item in self.fit_sample_ids):
            raise ValueError("fit_sample_ids must be nonempty")
        if len(set(self.fit_sample_ids)) != len(self.fit_sample_ids):
            raise ValueError("fit_sample_ids must be unique")
        if self.coefficient_sha256 != _coefficient_sha256(self.coefficient_vector):
            raise ValueError("coefficient_sha256 does not match coefficient_vector")


class GainGatedFiniteCorrectionBackend(Protocol):
    def materialize_velocity_correction(
        self, candidate: Any, solve_receipt: FiveMomentCorrectionSolveReceipt
    ) -> Any: ...

    def gain_gated_correction_provenance(
        self, candidate: Any, correction: Any
    ) -> GainGatedCorrectionBinding: ...


@dataclass(frozen=True)
class GainGatedFiniteCorrectionStageReceipt:
    candidate_before: CandidateProvenance
    moment_system_sha256: str
    local_moment_solve: FiveMomentCorrectionSolveReceipt
    correction_binding: GainGatedCorrectionBinding
    fixed_physical_stage: FixedPhysicalContractStageReceipt
    local_five_moment_gain_gate_passed: bool
    correction_materialized_only_after_gain_gate: bool
    correction_bound_to_exact_solve: bool
    held_out_excluded_from_correction_fit: bool
    stage_accepted: bool
    mechanics_only: bool
    candidate_residual_evidence: bool
    analytic_full_ns_cycle_gain_available: bool
    observed_gain_is_analytic_full_ns_bound: bool
    pde_validated: bool

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class GainGatedFiniteCorrectionStageResult:
    candidate_after_if_accepted: Any | None
    receipt: GainGatedFiniteCorrectionStageReceipt


def _require_backend(backend: Any) -> None:
    required = (
        "candidate_provenance",
        "candidate_identity",
        "five_moment_correction_system",
        "materialize_velocity_correction",
        "gain_gated_correction_provenance",
        "correction_provenance",
        "evaluate_full_ns_residual",
        "evaluate_correction",
        "apply_correction",
        "physical_contract_provenance",
        "correction_update_provenance",
        "residual_physical_contract_sha256",
    )
    missing = [name for name in required if not callable(getattr(backend, name, None))]
    if missing:
        raise GainGatedFiniteCorrectionError(
            "backend is missing gain-gated finite-correction methods: "
            + ", ".join(missing)
        )


def _require_partition_isolation(
    held_in: ValidationPartition, held_out: ValidationPartition
) -> None:
    if held_in.partition_id == held_out.partition_id:
        raise GainGatedFiniteCorrectionError(
            "held-in and held-out partition ids must differ"
        )
    overlap = set(held_in.sample_ids) & set(held_out.sample_ids)
    if overlap:
        raise GainGatedFiniteCorrectionError(
            "held-in and held-out sample ids overlap"
        )


def _require_complete_candidate(provenance: CandidateProvenance) -> None:
    if not isinstance(provenance, CandidateProvenance):
        raise GainGatedFiniteCorrectionError(
            "backend must return CandidateProvenance"
        )
    missing: list[str] = []
    if not provenance.complete_ns_defect:
        missing.append("complete_ns_defect")
    if not provenance.corrected_global_leading_join_complete:
        missing.append("corrected_global_leading_join_complete")
    if not provenance.matched_pressure_gradient_included:
        missing.append("matched_pressure_gradient_included")
    if not provenance.restricted_forcing_included:
        missing.append("restricted_forcing_included")
    if not provenance.restricted_forcing_preregistered:
        missing.append("restricted_forcing_preregistered")
    if provenance.residual_as_forcing_shortcut_used:
        missing.append("residual_as_forcing_shortcut_used=false")
    if missing:
        raise GainGatedFiniteCorrectionError(
            "candidate is not admitted before five-moment solve: "
            + ", ".join(missing)
        )


def _typed_moment_system(
    backend: GainGatedFiniteCorrectionBackend,
    candidate: Any,
    candidate_before: CandidateProvenance,
    held_in: ValidationPartition,
    held_out: ValidationPartition,
) -> tuple[CandidateIdentity, FiveMomentCorrectionSystem, str]:
    identity = backend.candidate_identity(candidate)
    if not isinstance(identity, CandidateIdentity):
        raise GainGatedFiniteCorrectionError(
            "backend must return CandidateIdentity for five-moment gain"
        )
    if identity.candidate_id != candidate_before.candidate_id:
        raise GainGatedFiniteCorrectionError(
            "five-moment candidate id differs from finite-cycle candidate"
        )
    if identity.candidate_sha256 != candidate_before.candidate_sha256:
        raise GainGatedFiniteCorrectionError(
            "five-moment candidate hash differs from finite-cycle candidate"
        )

    system = backend.five_moment_correction_system(candidate)
    if not isinstance(system, FiveMomentCorrectionSystem):
        raise GainGatedFiniteCorrectionError(
            "backend must return FiveMomentCorrectionSystem"
        )
    if system.source_candidate_id != candidate_before.candidate_id:
        raise GainGatedFiniteCorrectionError(
            "five-moment system source candidate id mismatch"
        )
    if system.source_candidate_sha256 != candidate_before.candidate_sha256:
        raise GainGatedFiniteCorrectionError(
            "five-moment system source candidate hash mismatch"
        )
    if system.fit_partition_id != held_in.partition_id:
        raise GainGatedFiniteCorrectionError(
            "five-moment discrepancy was not fit on declared held-in partition"
        )

    fit_ids = set(system.fit_sample_ids)
    held_in_ids = set(held_in.sample_ids)
    held_out_ids = set(held_out.sample_ids)
    if not fit_ids.issubset(held_in_ids):
        raise GainGatedFiniteCorrectionError(
            "five-moment fit samples are not a subset of held-in samples"
        )
    if fit_ids & held_out_ids:
        raise GainGatedFiniteCorrectionError(
            "held-out samples leaked into five-moment discrepancy fitting"
        )
    if not system.discrepancy_from_complete_ns_defect:
        raise GainGatedFiniteCorrectionError(
            "five-moment discrepancy is not derived from a complete NS defect"
        )
    if system.residual_as_forcing_shortcut_used:
        raise GainGatedFiniteCorrectionError(
            "five-moment system used forbidden residual-as-forcing shortcut"
        )

    return identity, system, _canonical_system_sha256(identity, system)


def _verify_correction_binding(
    binding: GainGatedCorrectionBinding,
    solve_receipt: FiveMomentCorrectionSolveReceipt,
    system: FiveMomentCorrectionSystem,
    candidate_before: CandidateProvenance,
) -> None:
    if not isinstance(binding, GainGatedCorrectionBinding):
        raise GainGatedFiniteCorrectionError(
            "backend must return GainGatedCorrectionBinding"
        )
    if binding.source_candidate_id != candidate_before.candidate_id:
        raise GainGatedFiniteCorrectionError(
            "materialized correction source candidate id mismatch"
        )
    if binding.source_candidate_sha256 != candidate_before.candidate_sha256:
        raise GainGatedFiniteCorrectionError(
            "materialized correction source candidate hash mismatch"
        )
    if binding.moment_system_sha256 != solve_receipt.diagnostic.system_sha256:
        raise GainGatedFiniteCorrectionError(
            "materialized correction is not bound to solved five-moment system"
        )
    solved = tuple(float(v) for v in solve_receipt.coefficient_vector)
    if tuple(float(v) for v in binding.coefficient_vector) != solved:
        raise GainGatedFiniteCorrectionError(
            "materialized correction coefficient vector differs from solved coefficients"
        )
    if binding.coefficient_sha256 != _coefficient_sha256(solved):
        raise GainGatedFiniteCorrectionError(
            "materialized correction coefficient hash mismatch"
        )
    if binding.fit_partition_id != system.fit_partition_id:
        raise GainGatedFiniteCorrectionError(
            "materialized correction fit partition differs from five-moment system"
        )
    if tuple(binding.fit_sample_ids) != tuple(system.fit_sample_ids):
        raise GainGatedFiniteCorrectionError(
            "materialized correction fit samples differ from five-moment system"
        )
    if not binding.derived_from_complete_ns_defect:
        raise GainGatedFiniteCorrectionError(
            "materialized correction lost complete-defect provenance"
        )
    if binding.residual_as_forcing_shortcut_used:
        raise GainGatedFiniteCorrectionError(
            "materialized correction used forbidden residual-as-forcing shortcut"
        )


def run_gain_gated_finite_correction_stage(
    backend: GainGatedFiniteCorrectionBackend,
    candidate: Any,
    held_in: ValidationPartition,
    held_out: ValidationPartition,
) -> GainGatedFiniteCorrectionStageResult:
    """Run one correction stage only after the local five-moment gain gate passes.

    The correction object is *not* a public argument.  It is materialized only
    from the authenticated solve receipt after the source-derived local
    contraction certificate succeeds.
    """

    _require_backend(backend)
    _require_partition_isolation(held_in, held_out)

    candidate_before = backend.candidate_provenance(candidate)
    _require_complete_candidate(candidate_before)

    _identity, system, system_sha256 = _typed_moment_system(
        backend, candidate, candidate_before, held_in, held_out
    )

    # #927 fails before iteration if the public-source smallness/self-map/
    # Lipschitz certificate does not hold.
    solve_receipt = solve_five_moment_correction(backend, candidate)
    if solve_receipt.diagnostic.system_sha256 != system_sha256:
        raise GainGatedFiniteCorrectionError(
            "five-moment payload changed between admission and solve"
        )
    if not solve_receipt.local_moment_fixed_point_converged:
        raise GainGatedFiniteCorrectionError(
            "five-moment solve returned without local fixed-point convergence"
        )
    if not solve_receipt.local_moment_correction_accepted:
        raise GainGatedFiniteCorrectionError(
            "five-moment solve did not accept the local correction"
        )

    correction = backend.materialize_velocity_correction(candidate, solve_receipt)
    binding = backend.gain_gated_correction_provenance(candidate, correction)
    _verify_correction_binding(binding, solve_receipt, system, candidate_before)

    fixed_result = run_fixed_physical_contract_finite_correction_stage(
        backend, candidate, correction, held_in, held_out
    )
    fixed_receipt = fixed_result.receipt

    mechanics_only = candidate_before.evidence_kind == "mechanics-only"
    stage_accepted = bool(fixed_receipt.stage_accepted)
    held_out_excluded = not (
        set(system.fit_sample_ids) & set(held_out.sample_ids)
    )
    receipt = GainGatedFiniteCorrectionStageReceipt(
        candidate_before=candidate_before,
        moment_system_sha256=system_sha256,
        local_moment_solve=solve_receipt,
        correction_binding=binding,
        fixed_physical_stage=fixed_receipt,
        local_five_moment_gain_gate_passed=True,
        correction_materialized_only_after_gain_gate=True,
        correction_bound_to_exact_solve=True,
        held_out_excluded_from_correction_fit=held_out_excluded,
        stage_accepted=stage_accepted,
        mechanics_only=mechanics_only,
        candidate_residual_evidence=bool(stage_accepted and not mechanics_only),
        analytic_full_ns_cycle_gain_available=False,
        observed_gain_is_analytic_full_ns_bound=False,
        pde_validated=False,
    )
    return GainGatedFiniteCorrectionStageResult(
        candidate_after_if_accepted=(
            fixed_result.candidate_after_if_accepted if stage_accepted else None
        ),
        receipt=receipt,
    )


def truth_boundary() -> dict[str, Any]:
    signature = inspect.signature(run_gain_gated_finite_correction_stage)
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
    return {
        "task": TASK,
        "schema": SCHEMA,
        "parent_agent3_pr": PARENT_AGENT3_PR,
        "parent_agent3_head": PARENT_AGENT3_HEAD,
        "parent_gain_source_blob": PARENT_GAIN_SOURCE_BLOB,
        "fixed_contract_source_blob": FIXED_CONTRACT_SOURCE_BLOB,
        "source_reader_repository": SOURCE_READER_REPOSITORY,
        "source_reader_head": SOURCE_READER_HEAD,
        "source_reader_tex_blob": SOURCE_READER_TEX_BLOB,
        "source_local_gain_equation": "B c + Q(c,c) = d",
        "source_local_gain_smallness_rule": "8*nu^2*q*||d||<=1",
        "correction_is_public_input": False,
        "correction_materialized_only_from_authenticated_solve": True,
        "five_moment_fit_must_be_held_in_only": True,
        "held_out_fit_leakage_rejected_before_solve": True,
        "complete_ns_defect_required_before_solve": True,
        "fixed_physical_contract_stage_required_after_materialization": True,
        "forbidden_public_parameters_absent": forbidden.isdisjoint(
            signature.parameters
        ),
        "local_five_moment_analytic_contraction_gate_integrated": True,
        "analytic_full_ns_cycle_gain_available": False,
        "observed_gain_is_analytic_full_ns_bound": False,
        "current_real_five_moment_discrepancy_materialized": False,
        "current_real_ns_correction_velocity_materialized": False,
        "real_candidate_finite_correction_cycle_run": False,
        "heldout_normalized_ns_residual_assessed": False,
        "residual_reduction_claimed": False,
        "same_protocol_comparable_to_st006": False,
        "pde_validated": False,
        "paper_exact": False,
        "openai_field_identified": False,
        "blowup_proved": False,
        "final_normalized_momentum_gate": FINAL_NORMALIZED_MOMENTUM_GATE,
        "final_normalized_divergence_gate": FINAL_NORMALIZED_DIVERGENCE_GATE,
    }


class _GainGatedMechanicsBackend(_FixedContractMechanicsBackend):
    """Mechanics-only backend proving stage plumbing, never NS evidence."""

    def candidate_identity(self, candidate: dict[str, float]) -> CandidateIdentity:
        finite = self.candidate_provenance(candidate)
        return CandidateIdentity(
            candidate_id=finite.candidate_id,
            candidate_sha256=finite.candidate_sha256,
            evidence_kind="mechanics-only",
        )

    def five_moment_correction_system(
        self, candidate: dict[str, float]
    ) -> FiveMomentCorrectionSystem:
        finite = self.candidate_provenance(candidate)
        B = np.eye(MOMENT_DIMENSION, dtype=float)
        Q = np.zeros(
            (MOMENT_DIMENSION, MOMENT_DIMENSION, MOMENT_DIMENSION), dtype=float
        )
        d = np.zeros(MOMENT_DIMENSION, dtype=float)
        d[0] = -0.5
        return FiveMomentCorrectionSystem(
            source_candidate_id=finite.candidate_id,
            source_candidate_sha256=finite.candidate_sha256,
            fit_partition_id="held-in",
            fit_sample_ids=("train-0", "train-1"),
            matrix_B=tuple(tuple(float(v) for v in row) for row in B),
            bilinear_Q=tuple(
                tuple(tuple(float(v) for v in row) for row in plane)
                for plane in Q
            ),
            discrepancy_d=tuple(float(v) for v in d),
            discrepancy_kind="mechanics-complete-defect-five-moment",
            discrepancy_from_complete_ns_defect=True,
            residual_as_forcing_shortcut_used=False,
        )

    def materialize_velocity_correction(
        self,
        candidate: dict[str, float],
        solve_receipt: FiveMomentCorrectionSolveReceipt,
    ) -> dict[str, Any]:
        coefficients = tuple(float(v) for v in solve_receipt.coefficient_vector)
        return {
            "delta": coefficients[0],
            "moment_system_sha256": solve_receipt.diagnostic.system_sha256,
            "coefficient_vector": coefficients,
            "coefficient_sha256": _coefficient_sha256(coefficients),
        }

    def gain_gated_correction_provenance(
        self, candidate: dict[str, float], correction: dict[str, Any]
    ) -> GainGatedCorrectionBinding:
        finite = self.candidate_provenance(candidate)
        system = self.five_moment_correction_system(candidate)
        coefficients = tuple(float(v) for v in correction["coefficient_vector"])
        return GainGatedCorrectionBinding(
            correction_id="mechanics-five-moment-to-half-amplitude",
            source_candidate_id=finite.candidate_id,
            source_candidate_sha256=finite.candidate_sha256,
            moment_system_sha256=str(correction["moment_system_sha256"]),
            coefficient_sha256=str(correction["coefficient_sha256"]),
            coefficient_vector=coefficients,
            fit_partition_id=system.fit_partition_id,
            fit_sample_ids=system.fit_sample_ids,
            derived_from_complete_ns_defect=True,
            residual_as_forcing_shortcut_used=False,
        )


def build_mechanics_report() -> dict[str, Any]:
    held_in = ValidationPartition(
        "held-in", ("train-0", "train-1", "train-2")
    )
    held_out = ValidationPartition(
        "held-out", ("test-0", "test-1", "test-2")
    )
    result = run_gain_gated_finite_correction_stage(
        _GainGatedMechanicsBackend(),
        {"amplitude": 1.0},
        held_in,
        held_out,
    )
    receipt = result.receipt
    return {
        "schema": SCHEMA,
        "task": TASK,
        "mechanics_only": True,
        "candidate_residual_evidence": False,
        "mechanics_stage_receipt": receipt.to_dict(),
        "observed_held_out_momentum_contraction_factor": (
            receipt.fixed_physical_stage.parent_stage.held_out_momentum_contraction_factor
        ),
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
