"""Current-I4 adapter from typed RF30 defect to the source RF34--RF39 compact correction.

This Kokuno Agent-3 increment deliberately reuses the two existing scientific
operators rather than inventing another correction:

* #1181: current-I4 fixed-Q + normalized auxiliary-T2 Haar -> typed RF30 state;
* #1134: typed RF30/RF31 backend -> RF34--RF39 compact mean correction.

The only new job is the identity-preserving adapter between them.  The
reserved-I4 patch geometry and the five RF31 compact basis columns are frozen
from the authenticated A1 outer schedule *before* the RF30 defect is evaluated.
No caller may provide P/J targets, a defect, basis columns, correction
coefficients, pressure, forcing, gain, held-out samples, or a scientific
threshold.

The corrected source relation used to map the reserved I4 interval from the
fixed-Q similarity coordinate to the RF34 profile radius is

    X = R^2/(2 s_Q),  x = R/sqrt(s_Q),  hence X = x^2/2.

The concrete C-infinity bump realization remains the autonomous repository
choice already frozen by #1134; this adapter does not claim source-exact bump
recovery.

At the current stack #1172 still has no checksum-pinned same-identity
repository auxiliary-T2 provider.  Consequently the executable correction is
mechanics/provenance evidence only.  It is not applied to the Cartesian
candidate here and it is not an NS residual improvement.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import inspect
import json
import math
from typing import Any

import numpy as np

from .kokuno_current_i4_nonlinear_mean_attribution import (
    AGENT1_LEADING_PR,
    AGENT2_COMPOSITE_PR,
    ExactCurrentI4NonlinearBackend,
)
from .kokuno_current_i4_rf30_fixedq_preflight import (
    CurrentI4FixedQObservablePreflight,
    materialize_current_i4_rf30_fixedq_preflight,
)
from .kokuno_current_i4_rf30_typed_defect import (
    CurrentI4RF30TypedDefectReceipt,
    materialize_current_i4_rf30_typed_defect,
)
from .kokuno_rf30_auxiliary_t2_haar_bridge import (
    PINNED_REPOSITORY_PROVIDER_BLOB,
    SourceAuxiliaryT2WaveProvider,
)
from .kokuno_rf30_rf31_typed_mean_correction import (
    CandidateIdentity,
    RF30PreMeanState,
    RF31CorrectionBasis,
)
from .kokuno_rf34_rf39_compact_mean_correction import (
    RF34RF39CorrectionReceipt,
    RF34RF39PatchContext,
    _canonical_chart_basis,
    materialize_rf34_rf39_compact_mean_correction,
)

TASK = "KOKUNO-A3-CURRENT-I4-RF34-RF39-CORRECTION-129"
SCHEMA = "kokuno-a3-current-i4-rf34-rf39-correction-v1"

PARENT_AGENT3_PR = 1181
PARENT_AGENT3_HEAD = "d157d2c718d7cb27935abb86229f7a2f0e8fbda3"
PARENT_TYPED_RF30_BLOB = "9b9704bcd188fe2e03a77e3d39425ea34668b1ed"
RF34_RF39_PARENT_PR = 1134
RF34_RF39_PARENT_BLOB = "392975eafd3787c91e2dc82ad9b735567513271d"

SOURCE_READER_REPO = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_READER_HEAD = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_READER_DATE = "2026-09-09"
SOURCE_READER_PATH = "navier-stokes/navier_stokes_workbench.tex"
SOURCE_READER_BLOB = "205a99807302e21a51c5eaf223390c0dfc42bcd0"
SOURCE_FORMULAS = ("RF30", "RF31", "RF34", "RF35", "RF36", "RF37", "RF38", "RF39", "RF43")

FINAL_NORMALIZED_MOMENTUM_GATE = 1.0e-3
FINAL_NORMALIZED_DIVERGENCE_GATE = 1.0e-5
SCHEDULE_ABS_TOL = 5.0e-12


class CurrentI4RF34RF39CorrectionError(RuntimeError):
    """Raised when current-I4 correction identity/provenance fails closed."""


@dataclass(frozen=True)
class CurrentI4RF34RF39CorrectionReceipt:
    candidate: CandidateIdentity
    typed_rf30_receipt_sha256: str
    source_chart_id: str
    source_chart_sha256: str
    patch_context: RF34RF39PatchContext
    correction: RF34RF39CorrectionReceipt
    basis_frozen_before_defect_evaluation: bool
    patch_context_frozen_before_defect_evaluation: bool
    repository_provider_blob_pinned: bool
    repository_candidate_correction_evidence: bool
    current_i4_rf31_system_materialized: bool
    current_i4_rf34_rf39_correction_materialized: bool
    correction_applied_to_candidate: bool
    rf44_rf49_nonlinear_remainder_recomputed: bool
    cartesian_correction_velocity_materialized: bool
    heldout_ns_residual_assessed: bool
    pde_validated: bool
    receipt_sha256: str

    def to_dict(self) -> dict[str, object]:
        return {
            **asdict(self),
            "schema": SCHEMA,
            "task": TASK,
            "source_provenance": source_provenance(),
            "truth_boundary": truth_boundary(),
        }


def _canonical_sha(payload: dict[str, object]) -> str:
    raw = json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)
    return hashlib.sha256(raw.encode("utf-8")).hexdigest()


def source_provenance() -> dict[str, object]:
    return {
        "repository": SOURCE_READER_REPO,
        "commit": SOURCE_READER_HEAD,
        "date": SOURCE_READER_DATE,
        "path": SOURCE_READER_PATH,
        "blob_sha1": SOURCE_READER_BLOB,
        "source_formulas": SOURCE_FORMULAS,
        "fixed_q_to_profile_radius": "X=R^2/(2s_Q), x=R/sqrt(s_Q), therefore X=x^2/2",
        "concrete_bump_realization": "autonomous_repository_choice_reused_from_PR_1134",
        "classification": "public_structure_provenance_not_independent_validation",
    }


def _chart_identity(
    preflight: CurrentI4FixedQObservablePreflight,
) -> tuple[str, str]:
    """Reproduce the #1181 current-I4 RF30 chart identity before defect evaluation."""
    payload: dict[str, object] = {
        "schema": "kokuno-a3-current-i4-fixed-q-rf30-chart-v1",
        "candidate_semantic_sha256": preflight.candidate_semantic_sha256,
        "fixed_q_preflight_sha256": preflight.preflight_sha256,
        "Q": preflight.Q,
        "s_Q": preflight.s_Q,
        "A": preflight.A,
        "D": preflight.D,
        "eta": preflight.eta,
        "Z": preflight.Z,
        "T": preflight.T,
        "radius_R": list(preflight.radius_R),
    }
    sha = _canonical_sha(payload)
    return f"current-I4-fixed-Q-{sha[:16]}", sha


def _freeze_patch_context_and_basis(
    backend: ExactCurrentI4NonlinearBackend,
    preflight: CurrentI4FixedQObservablePreflight,
) -> tuple[CandidateIdentity, RF34RF39PatchContext, RF31CorrectionBasis]:
    """Freeze all source-patch/basis choices without observing the RF30 defect."""
    leading = getattr(backend.composite_field, "leading_backend", None)
    if leading is None:
        raise CurrentI4RF34RF39CorrectionError("current-I4 leading backend is missing")
    schedule = getattr(leading, "outer_schedule", None)
    if schedule is None:
        raise CurrentI4RF34RF39CorrectionError("A1 leading lacks source-bound outer schedule")
    for name in ("lambda_outer", "h", "log_c_patch", "reserved_log_intervals"):
        if not hasattr(schedule, name):
            raise CurrentI4RF34RF39CorrectionError(
                f"A1 outer schedule lacks correction datum {name}"
            )

    lam = float(schedule.lambda_outer)
    h = float(schedule.h)
    log_c_patch = float(schedule.log_c_patch)
    if not all(math.isfinite(v) for v in (lam, h, log_c_patch)):
        raise CurrentI4RF34RF39CorrectionError("source patch parameters are non-finite")
    if not (0.0 < lam < 1.0):
        raise CurrentI4RF34RF39CorrectionError("source lambda must lie in (0,1)")
    if not (0.0 < h < 1.0e-2):
        raise CurrentI4RF34RF39CorrectionError("source similarity h must lie in (0,1e-2)")
    if not math.isclose(preflight.A - 0.5, h, rel_tol=0.0, abs_tol=SCHEDULE_ABS_TOL):
        raise CurrentI4RF34RF39CorrectionError("preflight A detached from source h")
    if not math.isclose(preflight.D, 0.5 - h, rel_tol=0.0, abs_tol=SCHEDULE_ABS_TOL):
        raise CurrentI4RF34RF39CorrectionError("preflight D detached from source h")
    if not preflight.all_strict_I4:
        raise CurrentI4RF34RF39CorrectionError("every correction probe must lie strictly in I4")

    intervals = schedule.reserved_log_intervals()
    if "I4" not in intervals:
        raise CurrentI4RF34RF39CorrectionError("source schedule lacks reserved mean interval I4")
    log_X_lo, log_X_hi = (float(v) for v in intervals["I4"])
    if not (math.isfinite(log_X_lo) and math.isfinite(log_X_hi) and log_X_lo < log_X_hi):
        raise CurrentI4RF34RF39CorrectionError("reserved I4 log-X interval is invalid")

    X = np.asarray(preflight.X, dtype=float)
    X_lo = math.exp(log_X_lo)
    X_hi = math.exp(log_X_hi)
    if np.any(~np.isfinite(X)) or np.any(X <= X_lo) or np.any(X >= X_hi):
        raise CurrentI4RF34RF39CorrectionError(
            "fixed-Q radial probes are not strictly inside source I4"
        )

    # X=x^2/2 on the fixed-Q chart. The patch is therefore fixed from source
    # I4 endpoints alone, without consulting P/J or any covariance.
    x_lo = math.sqrt(2.0 * X_lo)
    x_hi = math.sqrt(2.0 * X_hi)
    c_patch = math.exp(log_c_patch)
    if not (math.isfinite(c_patch) and c_patch > 0.0):
        raise CurrentI4RF34RF39CorrectionError("source c_patch must be positive")

    chart_id, chart_sha = _chart_identity(preflight)
    identity = CandidateIdentity(
        candidate_id=f"current-I4-A2-{AGENT2_COMPOSITE_PR}-A1-{AGENT1_LEADING_PR}",
        candidate_sha256=preflight.candidate_semantic_sha256,
        evidence_kind="mechanics-only",
    )
    context = RF34RF39PatchContext(
        source_candidate_id=identity.candidate_id,
        source_candidate_sha256=identity.candidate_sha256,
        source_chart_id=chart_id,
        source_chart_sha256=chart_sha,
        mean_patch_x_min=x_lo,
        mean_patch_x_max=x_hi,
        scale_s_q_over_Q=float(preflight.s_Q),
        similarity_A=float(preflight.A),
        eta=float(preflight.eta),
        lambda_value=lam,
        c_patch=c_patch,
        reserved_mean_patch_verified=True,
        context_frozen_before_defect_evaluation=True,
        heldout_samples_used_to_choose_context=False,
        source_exact_bump_claimed=False,
    )

    R = np.asarray(preflight.radius_R, dtype=float)
    dv, gd, _, _ = _canonical_chart_basis(context, R)
    basis = RF31CorrectionBasis(
        source_candidate_id=identity.candidate_id,
        source_candidate_sha256=identity.candidate_sha256,
        source_chart_id=chart_id,
        source_chart_sha256=chart_sha,
        radius_R=tuple(float(v) for v in R),
        delta_v_columns=tuple(tuple(float(v) for v in row) for row in dv),
        gamma_d_columns=tuple(tuple(float(v) for v in row) for row in gd),
        basis_selected_before_defect_evaluation=True,
        compact_support_verified=True,
        torus_independent_verified=True,
        surrogate_target_used_to_choose_basis=False,
        heldout_samples_used_to_choose_basis=False,
    )
    return identity, context, basis


class _FrozenCurrentI4RF34RF39Backend:
    """Typed adapter consumed by the already-existing #1134 source operator."""

    def __init__(
        self,
        *,
        identity: CandidateIdentity,
        state: RF30PreMeanState,
        basis: RF31CorrectionBasis,
        context: RF34RF39PatchContext,
    ) -> None:
        self._identity = identity
        self._state = state
        self._basis = basis
        self._context = context

    def candidate_identity(self, candidate: Any) -> CandidateIdentity:
        return self._identity

    def rf30_pre_mean_state(self, candidate: Any) -> RF30PreMeanState:
        return self._state

    def rf31_correction_basis(
        self, candidate: Any, state: RF30PreMeanState
    ) -> RF31CorrectionBasis:
        if state != self._state:
            raise CurrentI4RF34RF39CorrectionError("RF30 state changed after basis freeze")
        return self._basis

    def rf34_rf39_patch_context(
        self, candidate: Any, state: RF30PreMeanState
    ) -> RF34RF39PatchContext:
        if state != self._state:
            raise CurrentI4RF34RF39CorrectionError("RF30 state changed after context freeze")
        return self._context


def _bind_typed_state(
    provisional: CandidateIdentity,
    basis: RF31CorrectionBasis,
    context: RF34RF39PatchContext,
    typed: CurrentI4RF30TypedDefectReceipt,
    preflight: CurrentI4FixedQObservablePreflight,
) -> tuple[CandidateIdentity, RF31CorrectionBasis, RF34RF39PatchContext]:
    if typed.fixed_q_preflight_sha256 != preflight.preflight_sha256:
        raise CurrentI4RF34RF39CorrectionError("typed RF30 receipt detached from preflight")
    if typed.candidate.candidate_id != provisional.candidate_id:
        raise CurrentI4RF34RF39CorrectionError("typed RF30 candidate id drifted")
    if typed.candidate.candidate_sha256 != provisional.candidate_sha256:
        raise CurrentI4RF34RF39CorrectionError("typed RF30 candidate hash drifted")
    if typed.source_chart_id != context.source_chart_id:
        raise CurrentI4RF34RF39CorrectionError("typed RF30 source chart id drifted")
    if typed.source_chart_sha256 != context.source_chart_sha256:
        raise CurrentI4RF34RF39CorrectionError("typed RF30 source chart hash drifted")
    if typed.state.source_chart_sha256 != context.source_chart_sha256:
        raise CurrentI4RF34RF39CorrectionError("typed RF30 state/chart mismatch")
    if tuple(typed.state.radius_R) != tuple(basis.radius_R):
        raise CurrentI4RF34RF39CorrectionError("typed RF30 radial grid changed after basis freeze")

    parent_pinned = PINNED_REPOSITORY_PROVIDER_BLOB is not None
    if typed.repository_candidate_defect_evidence and not parent_pinned:
        raise CurrentI4RF34RF39CorrectionError(
            "repository RF30 evidence cannot promote while parent provider blob is unpinned"
        )

    # Frozen geometry is unchanged; only the evidence kind is inherited from
    # the authenticated #1181 receipt after its defect has been computed.
    identity = CandidateIdentity(
        candidate_id=provisional.candidate_id,
        candidate_sha256=provisional.candidate_sha256,
        evidence_kind=typed.candidate.evidence_kind,
    )
    return identity, basis, context


def materialize_current_i4_rf34_rf39_correction(
    backend: ExactCurrentI4NonlinearBackend,
    provider: SourceAuxiliaryT2WaveProvider,
    radius: Any,
    z: Any,
    t: Any,
) -> CurrentI4RF34RF39CorrectionReceipt:
    """Freeze the current-I4 compact basis, then consume actual typed RF30 mechanics."""
    if not isinstance(backend, ExactCurrentI4NonlinearBackend):
        raise TypeError("backend must be ExactCurrentI4NonlinearBackend")

    # Freeze source patch + basis *before* any RF30 defect materialization.
    preflight = materialize_current_i4_rf30_fixedq_preflight(backend, radius, z, t)
    provisional, context, basis = _freeze_patch_context_and_basis(backend, preflight)

    # Only now evaluate the raw-wave Haar covariance and RF30 P/J defect.
    typed = materialize_current_i4_rf30_typed_defect(
        backend, provider, radius, z, t
    )
    identity, basis, context = _bind_typed_state(
        provisional, basis, context, typed, preflight
    )

    adapter = _FrozenCurrentI4RF34RF39Backend(
        identity=identity,
        state=typed.state,
        basis=basis,
        context=context,
    )
    correction = materialize_rf34_rf39_compact_mean_correction(
        adapter, typed.receipt_sha256
    )

    if correction.candidate != identity:
        raise CurrentI4RF34RF39CorrectionError("RF34-RF39 candidate identity drifted")
    if correction.source_chart_sha256 != typed.source_chart_sha256:
        raise CurrentI4RF34RF39CorrectionError("RF34-RF39 source chart drifted")
    if (
        correction.repository_candidate_correction_evidence
        != typed.repository_candidate_defect_evidence
    ):
        raise CurrentI4RF34RF39CorrectionError("correction evidence kind drifted from RF30")

    repository_evidence = bool(
        correction.repository_candidate_correction_evidence
        and typed.repository_candidate_defect_evidence
        and PINNED_REPOSITORY_PROVIDER_BLOB is not None
    )
    payload: dict[str, object] = {
        "schema": SCHEMA,
        "task": TASK,
        "candidate": asdict(identity),
        "typed_rf30_receipt_sha256": typed.receipt_sha256,
        "source_chart_id": typed.source_chart_id,
        "source_chart_sha256": typed.source_chart_sha256,
        "patch_context": asdict(context),
        "correction": asdict(correction),
        "basis_frozen_before_defect_evaluation": True,
        "repository_provider_blob_pinned": PINNED_REPOSITORY_PROVIDER_BLOB is not None,
        "repository_candidate_correction_evidence": repository_evidence,
    }
    receipt_sha = _canonical_sha(payload)

    return CurrentI4RF34RF39CorrectionReceipt(
        candidate=identity,
        typed_rf30_receipt_sha256=typed.receipt_sha256,
        source_chart_id=typed.source_chart_id,
        source_chart_sha256=typed.source_chart_sha256,
        patch_context=context,
        correction=correction,
        basis_frozen_before_defect_evaluation=True,
        patch_context_frozen_before_defect_evaluation=True,
        repository_provider_blob_pinned=PINNED_REPOSITORY_PROVIDER_BLOB is not None,
        repository_candidate_correction_evidence=repository_evidence,
        current_i4_rf31_system_materialized=True,
        current_i4_rf34_rf39_correction_materialized=True,
        correction_applied_to_candidate=False,
        rf44_rf49_nonlinear_remainder_recomputed=False,
        cartesian_correction_velocity_materialized=False,
        heldout_ns_residual_assessed=False,
        pde_validated=False,
        receipt_sha256=receipt_sha,
    )


def truth_boundary() -> dict[str, object]:
    signature = inspect.signature(materialize_current_i4_rf34_rf39_correction)
    forbidden = {
        "P", "J_theta", "J_z", "defect", "residual", "covariance", "basis",
        "delta_v", "gamma_d", "coefficient", "bump", "forcing", "pressure",
        "gain", "damping", "threshold", "viscosity", "nu", "correction",
        "heldout", "lambda_value", "c_patch",
    }
    return {
        "schema": SCHEMA,
        "task": TASK,
        "parent_agent3_pr": PARENT_AGENT3_PR,
        "parent_agent3_head": PARENT_AGENT3_HEAD,
        "parent_typed_rf30_blob": PARENT_TYPED_RF30_BLOB,
        "rf34_rf39_parent_pr": RF34_RF39_PARENT_PR,
        "rf34_rf39_parent_blob": RF34_RF39_PARENT_BLOB,
        "source_reader_head": SOURCE_READER_HEAD,
        "source_formulas": SOURCE_FORMULAS,
        "public_parameters": tuple(signature.parameters),
        "forbidden_scientific_controls_exposed": bool(
            forbidden.intersection(signature.parameters)
        ),
        "current_i4_rf30_typed_defect_consumed": True,
        "rf31_basis_derived_from_source_i4_schedule_before_defect": True,
        "rf34_rf39_source_specific_compact_realization_reused": True,
        "caller_supplied_rf31_basis_allowed": False,
        "caller_supplied_defect_allowed": False,
        "agent2_curl_or_jacobian_reimplemented": False,
        "repository_provider_blob_pinned_in_parent": (
            PINNED_REPOSITORY_PROVIDER_BLOB is not None
        ),
        "current_i4_repository_candidate_rf34_rf39_correction_materialized": False,
        "correction_applied_to_candidate": False,
        "rf44_rf49_nonlinear_remainder_recomputed": False,
        "cartesian_correction_velocity_materialized": False,
        "finite_correction_cycle_run": False,
        "heldout_ns_residual_assessed": False,
        "residual_reduction_claimed": False,
        "same_protocol_comparable_to_st006": False,
        "final_normalized_momentum_gate": FINAL_NORMALIZED_MOMENTUM_GATE,
        "final_normalized_divergence_gate": FINAL_NORMALIZED_DIVERGENCE_GATE,
        "residual_defined_free_forcing_allowed": False,
        "finite_stage_small_residual_is_blowup_proof": False,
        "paper_exact": False,
        "openai_field_identified": False,
        "pde_validated": False,
    }


__all__ = [
    "CurrentI4RF34RF39CorrectionError",
    "CurrentI4RF34RF39CorrectionReceipt",
    "materialize_current_i4_rf34_rf39_correction",
    "source_provenance",
    "truth_boundary",
]
