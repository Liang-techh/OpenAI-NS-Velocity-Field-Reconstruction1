"""Three-resolution slow-derivative audit for the exact current-I4 A2 provider.

This Kokuno Agent-2 diagnostic is intentionally stacked on PR #1198.  It does
not create a new oscillatory field and does not modify Agent 3's admission pin.
Instead it consumes the raw auxiliary-``T^2`` wave already materialized by
``KokunoCurrentI4AuxiliaryT2Provider`` and recomputes its slow ``R``/``Z``
derivatives with an implementation-distinct centered-FD4 helper on a frozen
three-resolution ladder.

The public #1198 derivative realization uses relative step ``2^-12``.  The
present audit freezes ``2^-10, 2^-11, 2^-12`` before evaluation, reports the
successive numerical changes, and independently replays the production
``2^-12`` derivative from the raw lifted complete-curl wave.  Auxiliary torus
samples are shifted away from the uniform-grid nodes used by the parent tests.

This is numerical/operator-consistency evidence only.  The corrected Kokuno
source supplies the auxiliary-T2/Haar organization, but neither the repository
sign-channel phase lift nor these finite-difference derivatives are claimed to
be source-exact.  No covariance, RF30/RF49 defect, pressure, forcing, residual,
or correction coefficient is accepted by this diagnostic.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import inspect
import json
import math
from typing import Any

import numpy as np

from .kokuno_current_i4_auxiliary_t2_provider import (
    FD4_RELATIVE_STEP,
    KokunoCurrentI4AuxiliaryT2Provider,
)

TASK = "K2-OSC-112"
SCHEMA = "kokuno-a2-current-i4-auxiliary-t2-derivative-diagnostics-v1"
PARENT_AGENT2_PR = 1198
PARENT_AGENT2_HEAD = "5765c2b3bab7482df514efa87f9b6bba48e04b7f"
PARENT_PROVIDER_SOURCE_BLOB = "96168ac6583ad5e71aa044bd558c6feeadf051f4"

# Frozen before execution.  These are numerical audit resolutions, not NS
# acceptance thresholds and are not caller-configurable.
AUDIT_RELATIVE_STEPS = (2.0 ** -10, 2.0 ** -11, 2.0 ** -12)
AUDIT_TORUS_ORDER = 12
AUDIT_Y1_CELL_SHIFT = 0.37
AUDIT_Y2_CELL_SHIFT = 0.13

# Mechanical replay gate only: the independently written finest FD4 operator
# must reproduce the parent provider's frozen FD4 output.  This is not a PDE or
# source-correspondence acceptance threshold.
OPERATOR_REPLAY_RTOL = 5.0e-12
OPERATOR_REPLAY_ATOL = 5.0e-13


class AuxiliaryT2DerivativeDiagnosticError(RuntimeError):
    """Raised when the exact parent-provider audit contract is violated."""


def _canonical_json(payload: Any) -> str:
    return json.dumps(
        payload,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
        allow_nan=False,
    )


def _sha256(payload: Any) -> str:
    return hashlib.sha256(_canonical_json(payload).encode("utf-8")).hexdigest()


def _rms(value: Any) -> float:
    arr = np.asarray(value, dtype=float)
    if not np.all(np.isfinite(arr)):
        raise AuxiliaryT2DerivativeDiagnosticError("diagnostic array became non-finite")
    if arr.size == 0:
        raise AuxiliaryT2DerivativeDiagnosticError("diagnostic array is empty")
    return float(np.sqrt(np.mean(arr * arr)))


def _relative_rms_difference(new: Any, old: Any) -> float:
    a = np.asarray(new, dtype=float)
    b = np.asarray(old, dtype=float)
    if a.shape != b.shape or not np.all(np.isfinite(a)) or not np.all(np.isfinite(b)):
        raise AuxiliaryT2DerivativeDiagnosticError("relative-difference arrays are invalid")
    numerator = _rms(a - b)
    denominator = max(_rms(a), np.finfo(float).tiny)
    return float(numerator / denominator)


def _fixed_offgrid_torus() -> tuple[np.ndarray, np.ndarray]:
    n = AUDIT_TORUS_ORDER
    y1 = (np.arange(n, dtype=float) + AUDIT_Y1_CELL_SHIFT) / float(n)
    y2 = (np.arange(n, dtype=float) + AUDIT_Y2_CELL_SHIFT) / float(n)
    return np.meshgrid(y1, y2, indexing="ij")


def _stack_production_derivatives(sample: Any) -> tuple[np.ndarray, np.ndarray]:
    dR = np.stack(
        (
            np.asarray(sample.dR_w_r, dtype=float),
            np.asarray(sample.dR_w_theta, dtype=float),
            np.asarray(sample.dR_w_z, dtype=float),
        ),
        axis=-1,
    )
    dZ = np.stack(
        (
            np.asarray(sample.dZ_w_r, dtype=float),
            np.asarray(sample.dZ_w_theta, dtype=float),
            np.asarray(sample.dZ_w_z, dtype=float),
        ),
        axis=-1,
    )
    if not np.all(np.isfinite(dR)) or not np.all(np.isfinite(dZ)):
        raise AuxiliaryT2DerivativeDiagnosticError("production derivatives are non-finite")
    return dR, dZ


def _audit_fd4_from_raw_wave(
    provider: KokunoCurrentI4AuxiliaryT2Provider,
    *,
    relative_step: float,
    y1: np.ndarray,
    y2: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Independent FD4 realization using only the parent's raw lifted wave."""
    step = float(relative_step)
    if not math.isfinite(step) or step <= 0.0:
        raise AuxiliaryT2DerivativeDiagnosticError("relative_step must be positive and finite")

    chart = provider.chart
    R = np.asarray(chart.radius_R, dtype=float)
    Z = float(chart.Z)
    T = float(chart.T)

    hR = step * np.maximum(1.0, np.abs(R))
    if np.any(R - 2.0 * hR <= 0.0):
        raise AuxiliaryT2DerivativeDiagnosticError("audit radial FD4 stencil would cross R=0")

    # Deliberately do not call provider._fd4_derivatives here.  All derivative
    # values below are rebuilt from the raw lifted complete-curl wave.
    fp2 = provider._raw_wave(R + 2.0 * hR, Z, T, y1, y2)
    fp1 = provider._raw_wave(R + hR, Z, T, y1, y2)
    fm1 = provider._raw_wave(R - hR, Z, T, y1, y2)
    fm2 = provider._raw_wave(R - 2.0 * hR, Z, T, y1, y2)
    dR = (-fp2 + 8.0 * fp1 - 8.0 * fm1 + fm2) / (
        12.0 * hR[:, None, None, None]
    )

    hZ = step * max(1.0, abs(Z))
    zp2 = provider._raw_wave(R, Z + 2.0 * hZ, T, y1, y2)
    zp1 = provider._raw_wave(R, Z + hZ, T, y1, y2)
    zm1 = provider._raw_wave(R, Z - hZ, T, y1, y2)
    zm2 = provider._raw_wave(R, Z - 2.0 * hZ, T, y1, y2)
    dZ = (-zp2 + 8.0 * zp1 - 8.0 * zm1 + zm2) / (12.0 * hZ)

    if not np.all(np.isfinite(dR)) or not np.all(np.isfinite(dZ)):
        raise AuxiliaryT2DerivativeDiagnosticError("independent FD4 derivative became non-finite")
    return dR, dZ


@dataclass(frozen=True)
class AuxiliaryT2DerivativeDiagnosticReceipt:
    parent_agent2_pr: int
    parent_agent2_head: str
    parent_provider_source_blob_sha1: str
    candidate_semantic_sha256: str
    oscillatory_runtime_sha256: str
    leading_semantic_sha256: str
    provider_semantic_sha256: str
    provider_source_blob_sha1: str
    relative_steps: tuple[float, ...]
    torus_order: int
    y1_cell_shift: float
    y2_cell_shift: float
    raw_wave_rms: float
    dR_rms_by_step: tuple[float, ...]
    dZ_rms_by_step: tuple[float, ...]
    dR_successive_relative_rms_difference: tuple[float, ...]
    dZ_successive_relative_rms_difference: tuple[float, ...]
    production_dR_replay_max_abs: float
    production_dZ_replay_max_abs: float
    production_dR_replay_relative_rms: float
    production_dZ_replay_relative_rms: float
    production_derivative_operator_consistent: bool
    receipt_sha256: str

    def to_dict(self) -> dict[str, Any]:
        return {
            **asdict(self),
            "schema": SCHEMA,
            "task": TASK,
            "truth_boundary": truth_boundary(),
        }


def materialize_current_i4_auxiliary_t2_derivative_diagnostics(
    provider: KokunoCurrentI4AuxiliaryT2Provider,
) -> AuxiliaryT2DerivativeDiagnosticReceipt:
    """Audit exact #1198 slow derivatives without exposing tuning controls."""
    if not isinstance(provider, KokunoCurrentI4AuxiliaryT2Provider):
        raise TypeError("provider must be KokunoCurrentI4AuxiliaryT2Provider")

    metadata = provider.metadata
    if metadata.source_blob_sha1 != PARENT_PROVIDER_SOURCE_BLOB:
        raise AuxiliaryT2DerivativeDiagnosticError("exact #1198 provider source blob drifted")
    if not metadata.recomputed_from_actual_candidate:
        raise AuxiliaryT2DerivativeDiagnosticError("provider is not recomputed from the actual candidate")
    if metadata.surrogate_or_preaveraged_covariance_used:
        raise AuxiliaryT2DerivativeDiagnosticError("surrogate/pre-averaged covariance provider is forbidden")
    if metadata.heldout_data_used or metadata.residual_as_forcing_used:
        raise AuxiliaryT2DerivativeDiagnosticError("held-out or residual-as-forcing provider is forbidden")

    y1, y2 = _fixed_offgrid_torus()
    chart = provider.chart
    R = np.asarray(chart.radius_R, dtype=float)
    production = provider.sample_auxiliary_wave(
        R=R, Z=chart.Z, T=chart.T, y1=y1, y2=y2
    )
    production_dR, production_dZ = _stack_production_derivatives(production)
    raw_wave = np.stack(
        (
            np.asarray(production.w_r, dtype=float),
            np.asarray(production.w_theta, dtype=float),
            np.asarray(production.w_z, dtype=float),
        ),
        axis=-1,
    )

    audit_dR: list[np.ndarray] = []
    audit_dZ: list[np.ndarray] = []
    for step in AUDIT_RELATIVE_STEPS:
        dR, dZ = _audit_fd4_from_raw_wave(
            provider, relative_step=step, y1=y1, y2=y2
        )
        audit_dR.append(dR)
        audit_dZ.append(dZ)

    finest_dR = audit_dR[-1]
    finest_dZ = audit_dZ[-1]
    dR_abs = float(np.max(np.abs(finest_dR - production_dR)))
    dZ_abs = float(np.max(np.abs(finest_dZ - production_dZ)))
    dR_rel = _relative_rms_difference(finest_dR, production_dR)
    dZ_rel = _relative_rms_difference(finest_dZ, production_dZ)
    operator_consistent = bool(
        np.allclose(
            finest_dR,
            production_dR,
            rtol=OPERATOR_REPLAY_RTOL,
            atol=OPERATOR_REPLAY_ATOL,
        )
        and np.allclose(
            finest_dZ,
            production_dZ,
            rtol=OPERATOR_REPLAY_RTOL,
            atol=OPERATOR_REPLAY_ATOL,
        )
        and math.isclose(AUDIT_RELATIVE_STEPS[-1], FD4_RELATIVE_STEP, rel_tol=0.0, abs_tol=0.0)
    )

    core: dict[str, Any] = {
        "parent_agent2_pr": PARENT_AGENT2_PR,
        "parent_agent2_head": PARENT_AGENT2_HEAD,
        "parent_provider_source_blob_sha1": PARENT_PROVIDER_SOURCE_BLOB,
        "candidate_semantic_sha256": metadata.candidate_semantic_sha256,
        "oscillatory_runtime_sha256": metadata.oscillatory_runtime_sha256,
        "leading_semantic_sha256": metadata.leading_semantic_sha256,
        "provider_semantic_sha256": metadata.provider_semantic_sha256,
        "provider_source_blob_sha1": metadata.source_blob_sha1,
        "relative_steps": tuple(float(v) for v in AUDIT_RELATIVE_STEPS),
        "torus_order": AUDIT_TORUS_ORDER,
        "y1_cell_shift": AUDIT_Y1_CELL_SHIFT,
        "y2_cell_shift": AUDIT_Y2_CELL_SHIFT,
        "raw_wave_rms": _rms(raw_wave),
        "dR_rms_by_step": tuple(_rms(v) for v in audit_dR),
        "dZ_rms_by_step": tuple(_rms(v) for v in audit_dZ),
        "dR_successive_relative_rms_difference": tuple(
            _relative_rms_difference(audit_dR[i + 1], audit_dR[i])
            for i in range(len(audit_dR) - 1)
        ),
        "dZ_successive_relative_rms_difference": tuple(
            _relative_rms_difference(audit_dZ[i + 1], audit_dZ[i])
            for i in range(len(audit_dZ) - 1)
        ),
        "production_dR_replay_max_abs": dR_abs,
        "production_dZ_replay_max_abs": dZ_abs,
        "production_dR_replay_relative_rms": dR_rel,
        "production_dZ_replay_relative_rms": dZ_rel,
        "production_derivative_operator_consistent": operator_consistent,
    }
    return AuxiliaryT2DerivativeDiagnosticReceipt(
        **core,
        receipt_sha256=_sha256(core),
    )


def truth_boundary() -> dict[str, Any]:
    params = tuple(
        inspect.signature(materialize_current_i4_auxiliary_t2_derivative_diagnostics).parameters
    )
    forbidden = {
        "residual", "defect", "forcing", "pressure", "target", "gain", "damping",
        "threshold", "viscosity", "nu", "covariance", "correction", "fd_step",
        "relative_step", "phase", "amplitude", "heldout",
    }
    return {
        "schema": SCHEMA,
        "task": TASK,
        "parent_agent2_pr": PARENT_AGENT2_PR,
        "parent_agent2_head": PARENT_AGENT2_HEAD,
        "parent_provider_source_blob_sha1": PARENT_PROVIDER_SOURCE_BLOB,
        "public_parameters": params,
        "forbidden_scientific_controls_exposed": bool(forbidden.intersection(params)),
        "three_resolution_slow_derivative_diagnostic_materialized": True,
        "audit_derivatives_recomputed_from_raw_wave": True,
        "production_private_fd4_helper_used_by_audit_operator": False,
        "offgrid_auxiliary_phase_sampling_used": True,
        "production_fd4_finest_step_replayed": True,
        "numerical_stability_observations_are_pde_acceptance_gates": False,
        "source_exact_slow_derivatives_claimed": False,
        "source_exact_auxiliary_phase_assignment_recovered": False,
        "agent3_provider_pin_modified_here": False,
        "agent3_rf30_defect_materialized_here": False,
        "agent3_mean_correction_materialized_here": False,
        "forcing_or_pressure_added": False,
        "heldout_ns_residual_assessed": False,
        "residual_reduction_claimed": False,
        "paper_exact": False,
        "pde_validated": False,
    }


__all__ = [
    "AUDIT_RELATIVE_STEPS",
    "AUDIT_TORUS_ORDER",
    "AuxiliaryT2DerivativeDiagnosticError",
    "AuxiliaryT2DerivativeDiagnosticReceipt",
    "materialize_current_i4_auxiliary_t2_derivative_diagnostics",
    "truth_boundary",
]
