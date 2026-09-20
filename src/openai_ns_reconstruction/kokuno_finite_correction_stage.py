from __future__ import annotations

"""Typed one-step finite-correction machinery for Kokuno Agent 3.

This module deliberately does *not* invent a full Navier--Stokes defect for the
current strict-inner transport source.  Instead it defines the smallest
executable contract needed once such a defect exists:

    u^(k+1) = u^(k) + delta_u_k.

Residuals are evaluated by the backend on disjoint held-in / held-out
partitions before and after applying the correction.  Callers cannot pass
precomputed residuals, gains, damping factors, forcing values, pressure values,
or scientific thresholds into :func:`run_finite_correction_stage`.

The current Agent-3 radial-stress artifact is explicitly checked below and is
not admitted as a correction target because it is still based only on
pressure/forcing-free strict-inner transport.
"""

from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Protocol
import argparse
import json
import math


PARENT_AGENT3_HEAD = "29b385398d0c7636ed8d4049821a534d80728b33"
PARENT_AGENT3_PR = 882
FINAL_NORMALIZED_MOMENTUM_GATE = 1.0e-3
FINAL_NORMALIZED_DIVERGENCE_GATE = 1.0e-5
OBSERVED_CONTRACTION_REJECT_AT = 1.0

# The live constrained checkpoint reports that the replay-side finite-cycle
# closure scan has zero evaluable closures and unresolved bindings including
# rho0, rho_star, Delta E, v_q, L0 and nu.  Until those bindings are executable,
# an analytic Kokuno gain is not fabricated here.
KOKUNO_ANALYTIC_GAIN_BINDINGS_COMPLETE = False
KOKUNO_ANALYTIC_GAIN_AVAILABLE = False


class FiniteCorrectionContractError(RuntimeError):
    """Raised when a candidate/correction violates the typed cycle contract."""


@dataclass(frozen=True)
class ValidationPartition:
    partition_id: str
    sample_ids: tuple[str, ...]

    def __post_init__(self) -> None:
        if not self.partition_id:
            raise ValueError("partition_id must be nonempty")
        if not self.sample_ids:
            raise ValueError("validation partition must contain at least one sample")
        if any(not value for value in self.sample_ids):
            raise ValueError("sample ids must be nonempty")
        if len(set(self.sample_ids)) != len(self.sample_ids):
            raise ValueError("sample ids must be unique within a partition")


@dataclass(frozen=True)
class CandidateProvenance:
    candidate_id: str
    candidate_sha256: str
    evidence_kind: str
    complete_ns_defect: bool
    corrected_global_leading_join_complete: bool
    matched_pressure_gradient_included: bool
    restricted_forcing_included: bool
    restricted_forcing_preregistered: bool
    residual_as_forcing_shortcut_used: bool
    repository_acceptance_protocol: bool

    def __post_init__(self) -> None:
        if not self.candidate_id or not self.candidate_sha256:
            raise ValueError("candidate identity must be nonempty")
        if self.evidence_kind not in {"mechanics-only", "real-candidate"}:
            raise ValueError("evidence_kind must be mechanics-only or real-candidate")


@dataclass(frozen=True)
class CorrectionProvenance:
    correction_id: str
    source_candidate_id: str
    source_candidate_sha256: str
    correction_operator: str
    fit_partition_id: str
    fit_sample_ids: tuple[str, ...]
    derived_from_complete_ns_defect: bool
    residual_as_forcing_shortcut_used: bool

    def __post_init__(self) -> None:
        for value, name in (
            (self.correction_id, "correction_id"),
            (self.source_candidate_id, "source_candidate_id"),
            (self.source_candidate_sha256, "source_candidate_sha256"),
            (self.correction_operator, "correction_operator"),
            (self.fit_partition_id, "fit_partition_id"),
        ):
            if not value:
                raise ValueError(f"{name} must be nonempty")
        if not self.fit_sample_ids:
            raise ValueError("correction provenance must identify held-in fit samples")
        if len(set(self.fit_sample_ids)) != len(self.fit_sample_ids):
            raise ValueError("fit_sample_ids must be unique")


@dataclass(frozen=True)
class FullNSResidualMetrics:
    protocol_sha256: str
    partition_id: str
    sample_count: int
    normalized_momentum_sample_max: float
    normalized_momentum_grid_l2: float
    normalized_momentum_volume_l2: float
    normalized_divergence_sample_max: float
    normalized_divergence_grid_l2: float
    normalized_divergence_volume_l2: float
    complete_ns_residual: bool

    def __post_init__(self) -> None:
        if not self.protocol_sha256 or not self.partition_id:
            raise ValueError("residual protocol and partition identities must be nonempty")
        if self.sample_count <= 0:
            raise ValueError("sample_count must be positive")
        for name in (
            "normalized_momentum_sample_max",
            "normalized_momentum_grid_l2",
            "normalized_momentum_volume_l2",
            "normalized_divergence_sample_max",
            "normalized_divergence_grid_l2",
            "normalized_divergence_volume_l2",
        ):
            value = float(getattr(self, name))
            if not math.isfinite(value) or value < 0.0:
                raise ValueError(f"{name} must be finite and nonnegative")

    @property
    def momentum_score(self) -> float:
        return max(
            self.normalized_momentum_sample_max,
            self.normalized_momentum_grid_l2,
            self.normalized_momentum_volume_l2,
        )

    @property
    def divergence_score(self) -> float:
        return max(
            self.normalized_divergence_sample_max,
            self.normalized_divergence_grid_l2,
            self.normalized_divergence_volume_l2,
        )


@dataclass(frozen=True)
class CorrectionMetrics:
    partition_id: str
    l2_norm: float
    max_norm: float
    divergence_max: float
    nontriviality: float

    def __post_init__(self) -> None:
        if not self.partition_id:
            raise ValueError("correction metric partition must be nonempty")
        for name in ("l2_norm", "max_norm", "divergence_max", "nontriviality"):
            value = float(getattr(self, name))
            if not math.isfinite(value) or value < 0.0:
                raise ValueError(f"{name} must be finite and nonnegative")


class FiniteCorrectionBackend(Protocol):
    def candidate_provenance(self, candidate: Any) -> CandidateProvenance: ...

    def correction_provenance(
        self, candidate: Any, correction: Any
    ) -> CorrectionProvenance: ...

    def evaluate_full_ns_residual(
        self, candidate: Any, partition: ValidationPartition
    ) -> FullNSResidualMetrics: ...

    def evaluate_correction(
        self, correction: Any, partition: ValidationPartition
    ) -> CorrectionMetrics: ...

    def apply_correction(self, candidate: Any, correction: Any) -> Any: ...


@dataclass(frozen=True)
class FiniteCorrectionStageReceipt:
    candidate_before: CandidateProvenance
    candidate_after: CandidateProvenance
    correction: CorrectionProvenance
    held_in_partition: ValidationPartition
    held_out_partition: ValidationPartition
    held_in_before: FullNSResidualMetrics
    held_in_after: FullNSResidualMetrics
    held_out_before: FullNSResidualMetrics
    held_out_after: FullNSResidualMetrics
    correction_held_in: CorrectionMetrics
    correction_held_out: CorrectionMetrics
    held_in_momentum_contraction_factor: float
    held_out_momentum_contraction_factor: float
    held_in_divergence_contraction_factor: float
    held_out_divergence_contraction_factor: float
    stage_accepted: bool
    rejection_reasons: tuple[str, ...]
    mechanics_only: bool
    analytic_kokuno_gain_available: bool
    observed_gain_is_analytic_kokuno_bound: bool
    heldout_repository_numeric_gate_met: bool
    pde_validated: bool

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class FiniteCorrectionStageResult:
    candidate_after_if_accepted: Any | None
    receipt: FiniteCorrectionStageReceipt


def _ratio(after: float, before: float) -> float:
    if before == 0.0:
        return 0.0 if after == 0.0 else math.inf
    return after / before


def _require_backend_contract(backend: Any) -> None:
    required = (
        "candidate_provenance",
        "correction_provenance",
        "evaluate_full_ns_residual",
        "evaluate_correction",
        "apply_correction",
    )
    missing = [name for name in required if not callable(getattr(backend, name, None))]
    if missing:
        raise FiniteCorrectionContractError(
            "backend is missing typed finite-cycle methods: " + ", ".join(missing)
        )


def _validate_candidate_admission(provenance: CandidateProvenance) -> None:
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
        raise FiniteCorrectionContractError(
            "candidate is not admitted to a finite correction cycle: " + ", ".join(missing)
        )


def _validate_correction_provenance(
    provenance: CorrectionProvenance,
    candidate: CandidateProvenance,
    held_in: ValidationPartition,
    held_out: ValidationPartition,
) -> None:
    if provenance.source_candidate_id != candidate.candidate_id:
        raise FiniteCorrectionContractError("correction source candidate id mismatch")
    if provenance.source_candidate_sha256 != candidate.candidate_sha256:
        raise FiniteCorrectionContractError("correction source candidate hash mismatch")
    if provenance.fit_partition_id != held_in.partition_id:
        raise FiniteCorrectionContractError("correction was not fit on declared held-in partition")
    fit_ids = set(provenance.fit_sample_ids)
    held_in_ids = set(held_in.sample_ids)
    held_out_ids = set(held_out.sample_ids)
    if not fit_ids.issubset(held_in_ids):
        raise FiniteCorrectionContractError("correction fit samples are not a subset of held-in samples")
    if fit_ids & held_out_ids:
        raise FiniteCorrectionContractError("held-out samples leaked into correction fitting")
    if not provenance.derived_from_complete_ns_defect:
        raise FiniteCorrectionContractError("correction was not derived from a complete NS defect")
    if provenance.residual_as_forcing_shortcut_used:
        raise FiniteCorrectionContractError("correction provenance used forbidden residual-as-forcing shortcut")


def _validate_partition_pair(
    held_in: ValidationPartition, held_out: ValidationPartition
) -> None:
    if held_in.partition_id == held_out.partition_id:
        raise FiniteCorrectionContractError("held-in and held-out partition ids must differ")
    overlap = set(held_in.sample_ids) & set(held_out.sample_ids)
    if overlap:
        raise FiniteCorrectionContractError(
            "held-in and held-out sample ids overlap: " + ", ".join(sorted(overlap))
        )


def _validate_residual(
    metrics: FullNSResidualMetrics,
    partition: ValidationPartition,
    protocol_sha256: str | None,
) -> str:
    if not isinstance(metrics, FullNSResidualMetrics):
        raise FiniteCorrectionContractError("backend must return FullNSResidualMetrics")
    if metrics.partition_id != partition.partition_id:
        raise FiniteCorrectionContractError("residual partition identity mismatch")
    if metrics.sample_count != len(partition.sample_ids):
        raise FiniteCorrectionContractError("residual sample_count does not match partition")
    if not metrics.complete_ns_residual:
        raise FiniteCorrectionContractError("backend returned an incomplete NS residual")
    if protocol_sha256 is not None and metrics.protocol_sha256 != protocol_sha256:
        raise FiniteCorrectionContractError("residual protocol changed within the correction stage")
    return metrics.protocol_sha256


def _validate_correction_metrics(
    metrics: CorrectionMetrics, partition: ValidationPartition
) -> None:
    if not isinstance(metrics, CorrectionMetrics):
        raise FiniteCorrectionContractError("backend must return CorrectionMetrics")
    if metrics.partition_id != partition.partition_id:
        raise FiniteCorrectionContractError("correction metric partition identity mismatch")


def run_finite_correction_stage(
    backend: FiniteCorrectionBackend,
    candidate: Any,
    correction: Any,
    held_in: ValidationPartition,
    held_out: ValidationPartition,
) -> FiniteCorrectionStageResult:
    """Evaluate and conditionally admit one finite correction step.

    All residual and correction metrics are recomputed by ``backend``.  The
    caller supplies candidate/correction objects and disjoint partition
    identities only; there is no route for caller-provided residual scalars or
    gain values.
    """

    _require_backend_contract(backend)
    _validate_partition_pair(held_in, held_out)

    candidate_before = backend.candidate_provenance(candidate)
    if not isinstance(candidate_before, CandidateProvenance):
        raise FiniteCorrectionContractError("backend must return CandidateProvenance")
    _validate_candidate_admission(candidate_before)

    correction_provenance = backend.correction_provenance(candidate, correction)
    if not isinstance(correction_provenance, CorrectionProvenance):
        raise FiniteCorrectionContractError("backend must return CorrectionProvenance")
    _validate_correction_provenance(
        correction_provenance, candidate_before, held_in, held_out
    )

    protocol: str | None = None
    held_in_before = backend.evaluate_full_ns_residual(candidate, held_in)
    protocol = _validate_residual(held_in_before, held_in, protocol)
    held_out_before = backend.evaluate_full_ns_residual(candidate, held_out)
    protocol = _validate_residual(held_out_before, held_out, protocol)

    correction_held_in = backend.evaluate_correction(correction, held_in)
    correction_held_out = backend.evaluate_correction(correction, held_out)
    _validate_correction_metrics(correction_held_in, held_in)
    _validate_correction_metrics(correction_held_out, held_out)

    candidate_after_object = backend.apply_correction(candidate, correction)
    candidate_after = backend.candidate_provenance(candidate_after_object)
    if not isinstance(candidate_after, CandidateProvenance):
        raise FiniteCorrectionContractError("backend must return CandidateProvenance after correction")
    _validate_candidate_admission(candidate_after)
    if candidate_after.candidate_sha256 == candidate_before.candidate_sha256:
        raise FiniteCorrectionContractError("correction did not change the candidate identity")

    held_in_after = backend.evaluate_full_ns_residual(candidate_after_object, held_in)
    _validate_residual(held_in_after, held_in, protocol)
    held_out_after = backend.evaluate_full_ns_residual(candidate_after_object, held_out)
    _validate_residual(held_out_after, held_out, protocol)

    in_momentum_factor = _ratio(
        held_in_after.momentum_score, held_in_before.momentum_score
    )
    out_momentum_factor = _ratio(
        held_out_after.momentum_score, held_out_before.momentum_score
    )
    in_divergence_factor = _ratio(
        held_in_after.divergence_score, held_in_before.divergence_score
    )
    out_divergence_factor = _ratio(
        held_out_after.divergence_score, held_out_before.divergence_score
    )

    rejection_reasons: list[str] = []
    if correction_held_in.nontriviality <= 0.0 or correction_held_out.nontriviality <= 0.0:
        rejection_reasons.append("trivial_correction")
    if correction_held_in.l2_norm <= 0.0 or correction_held_out.l2_norm <= 0.0:
        rejection_reasons.append("zero_correction_l2_norm")
    if (
        correction_held_in.divergence_max > FINAL_NORMALIZED_DIVERGENCE_GATE
        or correction_held_out.divergence_max > FINAL_NORMALIZED_DIVERGENCE_GATE
    ):
        rejection_reasons.append("correction_divergence_exceeds_repository_gate")
    if in_momentum_factor > OBSERVED_CONTRACTION_REJECT_AT:
        rejection_reasons.append("held_in_momentum_increased")
    if out_momentum_factor > OBSERVED_CONTRACTION_REJECT_AT:
        rejection_reasons.append("held_out_momentum_increased")
    if in_divergence_factor > OBSERVED_CONTRACTION_REJECT_AT:
        rejection_reasons.append("held_in_divergence_increased")
    if out_divergence_factor > OBSERVED_CONTRACTION_REJECT_AT:
        rejection_reasons.append("held_out_divergence_increased")

    stage_accepted = not rejection_reasons
    mechanics_only = candidate_before.evidence_kind == "mechanics-only"
    heldout_repository_numeric_gate_met = bool(
        stage_accepted
        and not mechanics_only
        and candidate_after.repository_acceptance_protocol
        and held_out_after.momentum_score <= FINAL_NORMALIZED_MOMENTUM_GATE
        and held_out_after.divergence_score <= FINAL_NORMALIZED_DIVERGENCE_GATE
    )

    receipt = FiniteCorrectionStageReceipt(
        candidate_before=candidate_before,
        candidate_after=candidate_after,
        correction=correction_provenance,
        held_in_partition=held_in,
        held_out_partition=held_out,
        held_in_before=held_in_before,
        held_in_after=held_in_after,
        held_out_before=held_out_before,
        held_out_after=held_out_after,
        correction_held_in=correction_held_in,
        correction_held_out=correction_held_out,
        held_in_momentum_contraction_factor=float(in_momentum_factor),
        held_out_momentum_contraction_factor=float(out_momentum_factor),
        held_in_divergence_contraction_factor=float(in_divergence_factor),
        held_out_divergence_contraction_factor=float(out_divergence_factor),
        stage_accepted=stage_accepted,
        rejection_reasons=tuple(rejection_reasons),
        mechanics_only=mechanics_only,
        analytic_kokuno_gain_available=KOKUNO_ANALYTIC_GAIN_AVAILABLE,
        observed_gain_is_analytic_kokuno_bound=False,
        heldout_repository_numeric_gate_met=heldout_repository_numeric_gate_met,
        pde_validated=False,
    )
    return FiniteCorrectionStageResult(
        candidate_after_if_accepted=candidate_after_object if stage_accepted else None,
        receipt=receipt,
    )


def current_agent3_scoped_source_cycle_admission() -> dict[str, Any]:
    """Evaluate whether the exact parent Agent-3 scoped source may drive a cycle."""

    from openai_ns_reconstruction.kokuno_strict_inner_transport_first_cell_nested_refinement import (
        truth_boundary as parent_truth_boundary,
    )

    parent = parent_truth_boundary()
    required_parent_flags = {
        "complete_ns_defect": bool(parent.get("complete_ns_defect", False)),
        "pressure_gradient_included": bool(parent.get("pressure_gradient_included", False)),
        "restricted_forcing_included": bool(parent.get("restricted_forcing_included", False)),
        "scoped_transport_stress_authorized_as_correction_target": bool(
            parent.get("scoped_transport_stress_authorized_as_correction_target", False)
        ),
    }
    missing = tuple(name for name, value in required_parent_flags.items() if not value)
    return {
        "parent_agent3_pr": PARENT_AGENT3_PR,
        "parent_agent3_head": PARENT_AGENT3_HEAD,
        "source_kind": "strict-inner pressure/forcing-free transport radial stress",
        "required_parent_flags": required_parent_flags,
        "admitted_to_real_finite_correction_cycle": not missing,
        "missing_requirements": list(missing),
        "surrogate_defect_substitution_allowed": False,
        "residual_as_forcing_shortcut_allowed": False,
    }


def truth_boundary() -> dict[str, Any]:
    current = current_agent3_scoped_source_cycle_admission()
    return {
        "parent_agent3_pr": PARENT_AGENT3_PR,
        "parent_agent3_head": PARENT_AGENT3_HEAD,
        "typed_full_ns_residual_contract_defined": True,
        "held_in_held_out_identity_disjointness_enforced": True,
        "correction_fit_must_be_subset_of_held_in": True,
        "held_out_fit_leakage_rejected": True,
        "residuals_recomputed_by_backend_before_and_after": True,
        "caller_supplied_residual_or_gain_allowed": False,
        "residual_as_forcing_shortcut_allowed": False,
        "observed_heldout_contraction_is_engineering_diagnostic": True,
        "analytic_kokuno_gain_bindings_complete": KOKUNO_ANALYTIC_GAIN_BINDINGS_COMPLETE,
        "analytic_kokuno_gain_available": KOKUNO_ANALYTIC_GAIN_AVAILABLE,
        "observed_gain_is_analytic_kokuno_bound": False,
        "current_scoped_transport_source_admitted": current[
            "admitted_to_real_finite_correction_cycle"
        ],
        "current_scoped_transport_missing_requirements": current["missing_requirements"],
        "real_candidate_finite_correction_cycle_run": False,
        "heldout_normalized_ns_residual_assessed": False,
        "residual_reduction_claimed": False,
        "same_protocol_comparable_to_st006": False,
        "pde_validated": False,
        "final_normalized_momentum_gate": FINAL_NORMALIZED_MOMENTUM_GATE,
        "final_normalized_divergence_gate": FINAL_NORMALIZED_DIVERGENCE_GATE,
    }


class _MechanicsBackend:
    """Deterministic synthetic backend used only to exercise stage mechanics."""

    protocol_sha256 = "mechanics-protocol-v1"

    def candidate_provenance(self, candidate: dict[str, float]) -> CandidateProvenance:
        amplitude = float(candidate["amplitude"])
        return CandidateProvenance(
            candidate_id="mechanics-candidate",
            candidate_sha256=f"mechanics-amplitude-{amplitude:.16g}",
            evidence_kind="mechanics-only",
            complete_ns_defect=True,
            corrected_global_leading_join_complete=True,
            matched_pressure_gradient_included=True,
            restricted_forcing_included=True,
            restricted_forcing_preregistered=True,
            residual_as_forcing_shortcut_used=False,
            repository_acceptance_protocol=False,
        )

    def correction_provenance(
        self, candidate: dict[str, float], correction: dict[str, float]
    ) -> CorrectionProvenance:
        source = self.candidate_provenance(candidate)
        return CorrectionProvenance(
            correction_id="mechanics-half-amplitude",
            source_candidate_id=source.candidate_id,
            source_candidate_sha256=source.candidate_sha256,
            correction_operator="mechanics-linear-contraction",
            fit_partition_id="held-in",
            fit_sample_ids=("train-0", "train-1"),
            derived_from_complete_ns_defect=True,
            residual_as_forcing_shortcut_used=False,
        )

    def evaluate_full_ns_residual(
        self, candidate: dict[str, float], partition: ValidationPartition
    ) -> FullNSResidualMetrics:
        amplitude = abs(float(candidate["amplitude"]))
        scale = 1.0 if partition.partition_id == "held-in" else 1.2
        momentum = scale * amplitude
        divergence = scale * amplitude * 1.0e-8
        return FullNSResidualMetrics(
            protocol_sha256=self.protocol_sha256,
            partition_id=partition.partition_id,
            sample_count=len(partition.sample_ids),
            normalized_momentum_sample_max=1.10 * momentum,
            normalized_momentum_grid_l2=momentum,
            normalized_momentum_volume_l2=1.05 * momentum,
            normalized_divergence_sample_max=1.10 * divergence,
            normalized_divergence_grid_l2=divergence,
            normalized_divergence_volume_l2=1.05 * divergence,
            complete_ns_residual=True,
        )

    def evaluate_correction(
        self, correction: dict[str, float], partition: ValidationPartition
    ) -> CorrectionMetrics:
        magnitude = abs(float(correction["delta"]))
        return CorrectionMetrics(
            partition_id=partition.partition_id,
            l2_norm=magnitude,
            max_norm=magnitude,
            divergence_max=0.0,
            nontriviality=magnitude,
        )

    def apply_correction(
        self, candidate: dict[str, float], correction: dict[str, float]
    ) -> dict[str, float]:
        return {"amplitude": float(candidate["amplitude"]) + float(correction["delta"])}


def build_mechanics_report() -> dict[str, Any]:
    held_in = ValidationPartition("held-in", ("train-0", "train-1", "train-2"))
    held_out = ValidationPartition("held-out", ("test-0", "test-1", "test-2"))
    result = run_finite_correction_stage(
        _MechanicsBackend(),
        {"amplitude": 1.0},
        {"delta": -0.5},
        held_in,
        held_out,
    )
    return {
        "task": "KOKUNO-A3-FINITE-CORRECTION-STAGE-CONTRACT-077",
        "mechanics_only": True,
        "candidate_residual_evidence": False,
        "current_scoped_source_admission": current_agent3_scoped_source_cycle_admission(),
        "mechanics_stage_receipt": result.receipt.to_dict(),
        "truth_boundary": truth_boundary(),
    }


def _main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = build_mechanics_report()
    encoded = json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n"
    if args.output is None:
        print(encoded, end="")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(encoded, encoding="utf-8")


if __name__ == "__main__":
    _main()
