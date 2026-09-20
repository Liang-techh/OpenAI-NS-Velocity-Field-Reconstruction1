"""Evidence-bound canonical quadrature gate for the Kokuno full-NS validator.

This module is the Agent-4 response to CR002 PR #898.  Parent A4 PR #896
correctly separates held-out/sample RMS from canonical volume quadrature, but its
low-level API accepts a bare ``quadrature_ladder_assessed`` boolean.  A boolean
cannot carry the identity or numerical evidence of the registered CR001
24/48/96 volume-quadrature ladder.

The successor admission path below therefore requires a checksum-bound typed
quadrature receipt tied to the exact candidate, fixed physical contract,
implementation-distinct A4 operator identity, held-out protocol identity,
derivative ladder, validation times, per-time/per-order momentum/divergence
volume-L2 and kinetic-energy values, and frozen thresholds.  The old bare
boolean has no input on this public successor path.

No current Kokuno candidate is promoted by this module.  The strict-inner stack
still lacks the global/outer leading join, matched pressure, preregistered
restricted forcing, and a real correction velocity.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, replace
import hashlib
import inspect
import json
from typing import Any, Mapping, Sequence

import numpy as np

from . import kokuno_a4_blackbox_full_ns_validator as parent


TASK = "KOKUNO-A4-CANONICAL-QUADRATURE-EVIDENCE-082"
SCHEMA = "kokuno-a4-canonical-quadrature-evidence-v1"
RECEIPT_SCHEMA = "canonical-cr001-quadrature-receipt-v1"
PARENT_A4_PR = 896
PARENT_A4_HEAD = "e8a7712515f75bb8a0a5d86b9a6177ab9447f2b1"
PARENT_A4_SOURCE_BLOB = "a9e46daf60d909a63261dd10821e1a08101dd872"
CR002_AUDIT_PR = 898
CR002_AUDIT_HEAD = "3b34375303c7237eb8a300a6053aab72e7826bf6"
CR002_SCOPE_BLOB = "70665fbfd56fa582c2c448069f724d946fefcbd8"

QUADRATURE_ORDERS = (24, 48, 96)
# Additional fail-closed engineering stability gate, frozen before any real
# Kokuno canonical-quadrature output.  It only makes admission stricter; it does
# not replace or relax the registered 1e-3 / 1e-5 scientific gates.
QUADRATURE_STABILITY_GATE = 5.0e-2


def _canonical_json(payload: Mapping[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _sha256(payload: Mapping[str, Any]) -> str:
    return hashlib.sha256(_canonical_json(payload).encode("utf-8")).hexdigest()


VALIDATOR_OPERATOR_IDENTITY = {
    "parent_a4_pr": PARENT_A4_PR,
    "parent_a4_head": PARENT_A4_HEAD,
    "parent_a4_source_blob": PARENT_A4_SOURCE_BLOB,
    "equation": "u_t + (u dot grad)u + grad(p) - 0.01 Laplacian(u) - f",
    "derivative_operator": "Cartesian FD4 five-point spatial plus endpoint-aware FD4 time",
    "derivative_steps": list(parent.DERIVATIVE_STEPS),
    "canonical_volume_quadrature": "tensor-product registered per-axis ladder",
    "quadrature_orders_per_axis": list(QUADRATURE_ORDERS),
}
VALIDATOR_OPERATOR_SHA256 = _sha256(VALIDATOR_OPERATOR_IDENTITY)

HELDOUT_PROTOCOL_IDENTITY = {
    "seed": parent.HELDOUT_SEED,
    "heldout_points": parent.HELDOUT_POINTS,
    "validation_times": list(parent.VALIDATION_TIMES),
    "derivative_steps": list(parent.DERIVATIVE_STEPS),
    "momentum_gate": parent.MOMENTUM_GATE,
    "divergence_gate": parent.DIVERGENCE_GATE,
    "boundary_velocity_gate": parent.BOUNDARY_VELOCITY_GATE,
    "boundary_pressure_gate": parent.BOUNDARY_PRESSURE_GATE,
    "energy_reference_time": parent.ENERGY_REFERENCE_TIME,
    "energy_reference": parent.ENERGY_REFERENCE,
    "energy_reference_tolerance": parent.ENERGY_REFERENCE_TOLERANCE,
}
HELDOUT_PROTOCOL_SHA256 = _sha256(HELDOUT_PROTOCOL_IDENTITY)


@dataclass(frozen=True)
class QuadratureTimeMetrics:
    time: float
    momentum_volume_l2: float
    divergence_volume_l2: float
    kinetic_energy: float


@dataclass(frozen=True)
class QuadratureOrderMetrics:
    order_per_axis: int
    by_time: tuple[QuadratureTimeMetrics, ...]


@dataclass(frozen=True)
class FrozenQuadratureThresholds:
    momentum_volume_l2: float = parent.MOMENTUM_GATE
    divergence_volume_l2: float = parent.DIVERGENCE_GATE
    reference_energy: float = parent.ENERGY_REFERENCE
    reference_energy_abs_tolerance: float = parent.ENERGY_REFERENCE_TOLERANCE
    energy_min: float = parent.ENERGY_MIN
    energy_max: float = parent.ENERGY_MAX
    stability_relative_change: float = QUADRATURE_STABILITY_GATE


@dataclass(frozen=True)
class CanonicalQuadratureEvidence:
    schema: str
    task: str
    candidate_sha256: str
    stage: str
    physical_contract_sha256: str
    validator_operator_sha256: str
    heldout_protocol_sha256: str
    derivative_steps: tuple[float, float, float]
    validation_times: tuple[float, ...]
    quadrature_orders_per_axis: tuple[int, int, int]
    levels: tuple[QuadratureOrderMetrics, ...]
    thresholds: FrozenQuadratureThresholds
    max_scaled_change_48_to_96: float
    convergence_assessed: bool
    independent_of_training_and_held_in_selection: bool
    residual_defined_forcing: bool
    receipt_sha256: str

    def payload_without_digest(self) -> dict[str, Any]:
        payload = asdict(self)
        payload.pop("receipt_sha256", None)
        return payload

    def recomputed_sha256(self) -> str:
        return _sha256(self.payload_without_digest())


@dataclass(frozen=True)
class EvidenceBoundValidationReceipt:
    schema: str
    task: str
    candidate_id: str
    candidate_sha256: str
    stage: str
    physical_contract_sha256: str
    parent_validator_receipt: Mapping[str, Any]
    canonical_quadrature_evidence_sha256: str | None
    canonical_quadrature_assessed: bool
    quadrature_gate_failures: tuple[str, ...]
    combined_gate_failures: tuple[str, ...]
    heldout_normalized_ns_residual_assessed: bool
    pde_validated: bool

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def _normalize_levels(levels: Sequence[QuadratureOrderMetrics]) -> tuple[QuadratureOrderMetrics, ...]:
    ordered = tuple(sorted(levels, key=lambda item: int(item.order_per_axis)))
    if tuple(level.order_per_axis for level in ordered) != QUADRATURE_ORDERS:
        raise ValueError("canonical quadrature levels must be exactly 24/48/96 per axis")
    expected_times = tuple(float(v) for v in parent.VALIDATION_TIMES)
    normalized: list[QuadratureOrderMetrics] = []
    for level in ordered:
        rows = tuple(sorted(level.by_time, key=lambda row: float(row.time)))
        if tuple(float(row.time) for row in rows) != expected_times:
            raise ValueError("every quadrature order must contain all registered validation times")
        for row in rows:
            values = (
                float(row.time),
                float(row.momentum_volume_l2),
                float(row.divergence_volume_l2),
                float(row.kinetic_energy),
            )
            if not all(np.isfinite(v) for v in values):
                raise ValueError("canonical quadrature metrics must be finite")
            if row.momentum_volume_l2 < 0.0 or row.divergence_volume_l2 < 0.0 or row.kinetic_energy < 0.0:
                raise ValueError("canonical quadrature norms/energy must be nonnegative")
        normalized.append(QuadratureOrderMetrics(level.order_per_axis, rows))
    return tuple(normalized)


def _scaled_change(medium: float, fine: float, scale: float) -> float:
    return abs(float(medium) - float(fine)) / max(abs(float(fine)), float(scale), 1.0e-300)


def _max_scaled_change(levels: tuple[QuadratureOrderMetrics, ...]) -> float:
    medium = levels[1].by_time
    fine = levels[2].by_time
    changes: list[float] = []
    for m, f in zip(medium, fine):
        changes.append(_scaled_change(m.momentum_volume_l2, f.momentum_volume_l2, parent.MOMENTUM_GATE))
        changes.append(_scaled_change(m.divergence_volume_l2, f.divergence_volume_l2, parent.DIVERGENCE_GATE))
        changes.append(_scaled_change(m.kinetic_energy, f.kinetic_energy, parent.ENERGY_REFERENCE))
    return float(max(changes, default=0.0))


def build_canonical_quadrature_evidence(
    *,
    candidate_sha256: str,
    stage: str,
    physical_contract_sha256: str,
    levels: Sequence[QuadratureOrderMetrics],
) -> CanonicalQuadratureEvidence:
    """Build a checksum-bound receipt from an already-run canonical quadrature.

    This function does not claim to perform the expensive 24/48/96 run itself.
    It binds its per-time/per-order outputs so downstream PDE admission cannot be
    authorized by a bare truth value or by the held-out Monte-Carlo weights.
    """
    if not str(candidate_sha256):
        raise ValueError("candidate_sha256 must be nonempty")
    if stage not in ("leading_only", "leading_plus_oscillatory", "after_correction"):
        raise ValueError("unsupported staged candidate identity")
    if not str(physical_contract_sha256):
        raise ValueError("physical_contract_sha256 must be nonempty")
    normalized = _normalize_levels(levels)
    change = _max_scaled_change(normalized)
    draft = CanonicalQuadratureEvidence(
        schema=RECEIPT_SCHEMA,
        task=TASK,
        candidate_sha256=str(candidate_sha256),
        stage=stage,
        physical_contract_sha256=str(physical_contract_sha256),
        validator_operator_sha256=VALIDATOR_OPERATOR_SHA256,
        heldout_protocol_sha256=HELDOUT_PROTOCOL_SHA256,
        derivative_steps=parent.DERIVATIVE_STEPS,
        validation_times=parent.VALIDATION_TIMES,
        quadrature_orders_per_axis=QUADRATURE_ORDERS,
        levels=normalized,
        thresholds=FrozenQuadratureThresholds(),
        max_scaled_change_48_to_96=change,
        convergence_assessed=True,
        independent_of_training_and_held_in_selection=True,
        residual_defined_forcing=False,
        receipt_sha256="",
    )
    return replace(draft, receipt_sha256=draft.recomputed_sha256())


def canonical_quadrature_failures(
    evidence: CanonicalQuadratureEvidence,
    *,
    candidate_sha256: str,
    stage: str,
    physical_contract_sha256: str,
) -> list[str]:
    failures: list[str] = []
    if evidence.schema != RECEIPT_SCHEMA:
        failures.append("quadrature_receipt_schema")
    if evidence.task != TASK:
        failures.append("quadrature_receipt_task")
    if evidence.receipt_sha256 != evidence.recomputed_sha256():
        failures.append("quadrature_receipt_digest")
    if evidence.candidate_sha256 != candidate_sha256:
        failures.append("quadrature_candidate_identity")
    if evidence.stage != stage:
        failures.append("quadrature_stage_identity")
    if evidence.physical_contract_sha256 != physical_contract_sha256:
        failures.append("quadrature_physical_contract_identity")
    if evidence.validator_operator_sha256 != VALIDATOR_OPERATOR_SHA256:
        failures.append("quadrature_operator_identity")
    if evidence.heldout_protocol_sha256 != HELDOUT_PROTOCOL_SHA256:
        failures.append("quadrature_heldout_protocol_identity")
    if tuple(evidence.derivative_steps) != parent.DERIVATIVE_STEPS:
        failures.append("quadrature_derivative_ladder")
    if tuple(evidence.validation_times) != parent.VALIDATION_TIMES:
        failures.append("quadrature_validation_times")
    if tuple(evidence.quadrature_orders_per_axis) != QUADRATURE_ORDERS:
        failures.append("quadrature_order_ladder")
    if evidence.thresholds != FrozenQuadratureThresholds():
        failures.append("quadrature_threshold_identity")
    if not evidence.independent_of_training_and_held_in_selection:
        failures.append("quadrature_not_independent")
    if evidence.residual_defined_forcing:
        failures.append("quadrature_residual_defined_forcing")
    if not evidence.convergence_assessed:
        failures.append("quadrature_convergence_unassessed")

    try:
        levels = _normalize_levels(evidence.levels)
    except ValueError:
        failures.append("quadrature_metrics_shape")
        return failures

    recomputed_change = _max_scaled_change(levels)
    if not np.isclose(
        evidence.max_scaled_change_48_to_96,
        recomputed_change,
        rtol=0.0,
        atol=1.0e-15,
    ):
        failures.append("quadrature_convergence_metric_identity")
    if recomputed_change > QUADRATURE_STABILITY_GATE:
        failures.append("quadrature_resolution_stability")

    fine = levels[-1]
    for row in fine.by_time:
        if row.momentum_volume_l2 > parent.MOMENTUM_GATE:
            failures.append(f"quadrature_momentum_volume_l2@{row.time:.10g}")
        if row.divergence_volume_l2 > parent.DIVERGENCE_GATE:
            failures.append(f"quadrature_divergence_volume_l2@{row.time:.10g}")
        if row.kinetic_energy < parent.ENERGY_MIN or row.kinetic_energy > parent.ENERGY_MAX:
            failures.append(f"quadrature_energy_range@{row.time:.10g}")
        if np.isclose(row.time, parent.ENERGY_REFERENCE_TIME, rtol=0.0, atol=1.0e-14):
            if abs(row.kinetic_energy - parent.ENERGY_REFERENCE) > parent.ENERGY_REFERENCE_TOLERANCE:
                failures.append("quadrature_reference_energy")
    return failures


def validate_complete_candidate_with_canonical_quadrature(
    candidate: parent.BlackBoxFullNSCandidate,
    contract: parent.FixedPhysicalContract,
    x: Any,
    y: Any,
    z: Any,
    t: Any,
    *,
    expected_stage: str,
    volume_weights: Any | None = None,
    canonical_quadrature_evidence: CanonicalQuadratureEvidence | None = None,
) -> EvidenceBoundValidationReceipt:
    """Authoritative successor admission path for eventual PDE promotion.

    Parent #896 is always invoked with its quadrature boolean forced ``False``.
    Thus caller-supplied sample weights can support held-out diagnostics but can
    never stand in for canonical quadrature evidence on this path.
    """
    parent_receipt = parent.validate_complete_candidate(
        candidate,
        contract,
        x,
        y,
        z,
        t,
        expected_stage=expected_stage,
        volume_weights=volume_weights,
        quadrature_ladder_assessed=False,
    )
    metadata = dict(candidate.validation_metadata())
    parent_failures = tuple(
        failure
        for failure in parent_receipt.gate_failures
        if failure != "quadrature_ladder_unassessed"
    )
    if canonical_quadrature_evidence is None:
        quadrature_failures = ("canonical_quadrature_evidence_missing",)
        evidence_sha = None
    else:
        quadrature_failures = tuple(
            canonical_quadrature_failures(
                canonical_quadrature_evidence,
                candidate_sha256=str(metadata["candidate_sha256"]),
                stage=expected_stage,
                physical_contract_sha256=contract.physical_contract_sha256,
            )
        )
        evidence_sha = canonical_quadrature_evidence.receipt_sha256
    combined = parent_failures + quadrature_failures
    quadrature_assessed = canonical_quadrature_evidence is not None and not quadrature_failures
    pde_validated = bool(
        parent_receipt.heldout_normalized_ns_residual_assessed
        and quadrature_assessed
        and not combined
    )
    return EvidenceBoundValidationReceipt(
        schema=SCHEMA,
        task=TASK,
        candidate_id=str(metadata["candidate_id"]),
        candidate_sha256=str(metadata["candidate_sha256"]),
        stage=expected_stage,
        physical_contract_sha256=contract.physical_contract_sha256,
        parent_validator_receipt=parent_receipt.to_dict(),
        canonical_quadrature_evidence_sha256=evidence_sha,
        canonical_quadrature_assessed=quadrature_assessed,
        quadrature_gate_failures=quadrature_failures,
        combined_gate_failures=combined,
        heldout_normalized_ns_residual_assessed=parent_receipt.heldout_normalized_ns_residual_assessed,
        pde_validated=pde_validated,
    )


def current_ineligibility() -> dict[str, Any]:
    state = dict(parent.current_strict_inner_ineligibility())
    state.update(
        {
            "canonical_quadrature_evidence_available": False,
            "canonical_quadrature_receipt_sha256": None,
            "bare_quadrature_boolean_authorizes_pde_promotion": False,
            "pde_validated": False,
        }
    )
    return state


def public_contract() -> dict[str, Any]:
    signature = inspect.signature(validate_complete_candidate_with_canonical_quadrature)
    return {
        "schema": SCHEMA,
        "task": TASK,
        "parent_a4_pr": PARENT_A4_PR,
        "parent_a4_head": PARENT_A4_HEAD,
        "parent_a4_source_blob": PARENT_A4_SOURCE_BLOB,
        "cr002_audit_pr": CR002_AUDIT_PR,
        "cr002_audit_head": CR002_AUDIT_HEAD,
        "cr002_scope_blob": CR002_SCOPE_BLOB,
        "validator_operator_sha256": VALIDATOR_OPERATOR_SHA256,
        "heldout_protocol_sha256": HELDOUT_PROTOCOL_SHA256,
        "derivative_steps": list(parent.DERIVATIVE_STEPS),
        "validation_times": list(parent.VALIDATION_TIMES),
        "quadrature_orders_per_axis": list(QUADRATURE_ORDERS),
        "quadrature_stability_gate": QUADRATURE_STABILITY_GATE,
        "frozen_gates": {
            "momentum_max": parent.MOMENTUM_GATE,
            "momentum_volume_l2": parent.MOMENTUM_GATE,
            "divergence_max": parent.DIVERGENCE_GATE,
            "divergence_volume_l2": parent.DIVERGENCE_GATE,
            "reference_energy_abs_tolerance": parent.ENERGY_REFERENCE_TOLERANCE,
        },
        "successor_validator_parameters": list(signature.parameters),
        "bare_quadrature_ladder_assessed_parameter_present": (
            "quadrature_ladder_assessed" in signature.parameters
        ),
        "caller_sample_weights_are_canonical_quadrature": False,
        "typed_canonical_quadrature_evidence_required_for_pde_promotion": True,
        "current_ineligibility": current_ineligibility(),
        "truth_boundary": {
            "real_canonical_quadrature_run_completed": False,
            "real_complete_kokuno_candidate_available": False,
            "training_loss_is_validation": False,
            "residual_defined_free_forcing_allowed": False,
            "threshold_relaxation_allowed": False,
            "pde_validated": False,
        },
    }
