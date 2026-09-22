"""Current-I4 typed RF30 defect from the normalized auxiliary-T2 Haar bridge.

This is one narrow Kokuno Agent-3 adapter stacked on PR #1172.  The parent
bridge already enforces the important semantic distinction between a physical
azimuthal average and the corrected-reader RF30 normalized Haar mean on the
source auxiliary two-torus.  This module advances exactly one seam:

    exact current-I4 fixed-Q candidate
      -> raw auxiliary-T2 wave provider
      -> #1172 normalized-Haar covariance/slow derivatives
      -> #1126 RF30PreMeanState
      -> #1126 RF30 (P, J_theta, J_z) defect tuple.

The reserved-I4 source background is not caller supplied.  It is reconstructed
from the exact Agent-1 #1079 source-bound outer schedule:

    G = 0,
    V = a(eta) s_Q^(lambda-h) R^(-1-2 lambda),
    a(eta) = 2^(1/2+lambda) c_patch / (1+eta^2),
    h = A - 1/2.

No RF31 basis is chosen here, no correction coefficient is solved, and no
correction is applied.  The public entry point does not accept P/J defects,
pre-averaged covariance, residual, forcing, pressure, gain, held-out samples,
or scientific thresholds.  At this exact stack, #1172 has no checksum-pinned
repository auxiliary-T2 provider, so repository-candidate RF30 evidence remains
fail-closed false even though the mechanics path is executable.
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
    AGENT1_LEADING_HEAD,
    AGENT1_LEADING_PR,
    AGENT2_COMPOSITE_HEAD,
    AGENT2_COMPOSITE_PR,
    ExactCurrentI4NonlinearBackend,
)
from .kokuno_current_i4_rf30_fixedq_preflight import (
    CurrentI4FixedQObservablePreflight,
    materialize_current_i4_rf30_fixedq_preflight,
)
from .kokuno_rf30_auxiliary_t2_haar_bridge import (
    PINNED_REPOSITORY_PROVIDER_BLOB,
    AuxiliaryT2HaarCovarianceReceipt,
    SourceAuxiliaryT2WaveProvider,
    materialize_current_i4_rf30_auxiliary_haar_bridge,
)
from .kokuno_rf30_rf31_typed_mean_correction import (
    CandidateIdentity,
    RF30DefectTuple,
    RF30PreMeanState,
    _rf30_defect,
)

TASK = "KOKUNO-A3-CURRENT-I4-RF30-TYPED-DEFECT-128"
SCHEMA = "kokuno-a3-current-i4-rf30-typed-defect-v1"
PARENT_AGENT3_PR = 1172
PARENT_AGENT3_HEAD = "9b0cecdf7eb9195bd468c3ddcd026c9d28d8e8d7"
PARENT_HAAR_BRIDGE_BLOB = "5b20d4945f3aeb177e327e8d6b5709a6726c885a"
RF30_TYPED_PARENT_PR = 1126
RF30_TYPED_PARENT_HEAD = "edcf3a3a2e87f335c268460ca0d4b3b23f4a9847"
RF30_TYPED_PARENT_BLOB = "fa31fdf12adceaa5b295db502cf0ccaaad778307"
AGENT1_I4_SOURCE_BLOB = "6f04ce0a856b44430402576dad88438da90d1ebb"

SOURCE_READER_REPO = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_READER_HEAD = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_READER_DATE = "2026-09-09"
SOURCE_READER_PATH = "navier-stokes/navier_stokes_workbench.tex"
SOURCE_READER_BLOB = "205a99807302e21a51c5eaf223390c0dfc42bcd0"
SOURCE_BACKGROUND_FORMULA = (
    "RF43 on I_mean=I4: G=0; "
    "V=a(eta)*s_Q^(lambda-h)*R^(-1-2lambda); "
    "a(eta)=2^(1/2+lambda)*c_patch/(1+eta^2)"
)

FINAL_NORMALIZED_MOMENTUM_GATE = 1.0e-3
FINAL_NORMALIZED_DIVERGENCE_GATE = 1.0e-5
BACKGROUND_REL_TOL = 5.0e-12


class CurrentI4RF30TypedDefectError(RuntimeError):
    """Raised when current-I4 RF30 identity/background provenance fails closed."""


def _canonical_sha(payload: dict[str, object]) -> str:
    raw = json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)
    return hashlib.sha256(raw.encode("utf-8")).hexdigest()


def _finite_scalar(value: Any, label: str) -> float:
    out = float(value)
    if not math.isfinite(out):
        raise CurrentI4RF30TypedDefectError(f"{label} must be finite")
    return out


def _relative_error(a: float, b: float) -> float:
    return abs(a - b) / max(abs(a), abs(b), np.finfo(float).tiny)


@dataclass(frozen=True)
class CurrentI4RF30TypedDefectReceipt:
    candidate: CandidateIdentity
    source_chart_id: str
    source_chart_sha256: str
    fixed_q_preflight_sha256: str
    auxiliary_haar_receipt_sha256: str
    provider_semantic_sha256: str
    provider_source_blob_sha1: str
    provider_kind: str
    state: RF30PreMeanState
    defect: RF30DefectTuple
    lambda_value: float
    h_value: float
    log_c_patch: float
    base_V_min_abs: float
    base_V_max_abs: float
    base_G_max_abs: float
    source_background_from_exact_a1_i4_schedule: bool
    normalized_source_auxiliary_t2_haar_mean_used: bool
    typed_rf30_state_materialized: bool
    rf30_defect_materialized: bool
    repository_candidate_defect_evidence: bool
    rf31_five_row_system_materialized: bool
    correction_applied: bool
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


def source_provenance() -> dict[str, object]:
    return {
        "repository": SOURCE_READER_REPO,
        "commit": SOURCE_READER_HEAD,
        "date": SOURCE_READER_DATE,
        "path": SOURCE_READER_PATH,
        "blob_sha1": SOURCE_READER_BLOB,
        "background_formula": SOURCE_BACKGROUND_FORMULA,
        "classification": "public_structure_provenance_not_independent_validation",
    }


def _exact_i4_background(
    backend: ExactCurrentI4NonlinearBackend,
    preflight: CurrentI4FixedQObservablePreflight,
) -> tuple[np.ndarray, np.ndarray, float, float, float]:
    leading = getattr(backend.composite_field, "leading_backend", None)
    if leading is None:
        raise CurrentI4RF30TypedDefectError("exact current-I4 leading backend is missing")
    if str(getattr(leading, "semantic_sha256", "")) != preflight.leading_semantic_sha256:
        raise CurrentI4RF30TypedDefectError("A1 leading semantic identity drifted")

    schedule = getattr(leading, "outer_schedule", None)
    if schedule is None:
        raise CurrentI4RF30TypedDefectError("A1 leading lacks the source-bound outer schedule")
    for name in ("lambda_outer", "h", "log_c_patch", "reserved_log_intervals"):
        if not hasattr(schedule, name):
            raise CurrentI4RF30TypedDefectError(
                f"A1 outer schedule lacks source background datum {name}"
            )

    lam = _finite_scalar(schedule.lambda_outer, "lambda_outer")
    h = _finite_scalar(schedule.h, "h")
    log_c_patch = _finite_scalar(schedule.log_c_patch, "log_c_patch")
    if not (0.0 < lam < 1.0):
        raise CurrentI4RF30TypedDefectError("source lambda must lie in (0,1)")
    if not (0.0 < h < 1.0e-2):
        raise CurrentI4RF30TypedDefectError("source similarity h must lie in (0,1e-2)")
    if _relative_error(preflight.A - 0.5, h) > BACKGROUND_REL_TOL:
        raise CurrentI4RF30TypedDefectError("A1 similarity A is detached from outer-schedule h")
    if _relative_error(preflight.D, 0.5 - h) > BACKGROUND_REL_TOL:
        raise CurrentI4RF30TypedDefectError("A1 similarity D is detached from outer-schedule h")

    intervals = schedule.reserved_log_intervals()
    if "I4" not in intervals:
        raise CurrentI4RF30TypedDefectError("A1 outer schedule lacks reserved interval I4")
    i4_lo, i4_hi = (float(v) for v in intervals["I4"])
    leading_lo = _finite_scalar(getattr(leading, "X_I4_start", np.nan), "X_I4_start")
    leading_hi = _finite_scalar(getattr(leading, "X_I4_end", np.nan), "X_I4_end")
    if not (
        math.isclose(math.log(leading_lo), i4_lo, rel_tol=0.0, abs_tol=5.0e-12)
        and math.isclose(math.log(leading_hi), i4_hi, rel_tol=0.0, abs_tol=5.0e-12)
    ):
        raise CurrentI4RF30TypedDefectError("A1 I4 interval is detached from the outer schedule")
    if not preflight.all_strict_I4:
        raise CurrentI4RF30TypedDefectError("RF30 typed state requires every probe to lie strictly in I4")

    R = np.asarray(preflight.radius_R, dtype=float)
    if R.ndim != 1 or np.any(~np.isfinite(R)) or np.any(R <= 0.0):
        raise CurrentI4RF30TypedDefectError("fixed-Q R grid is invalid")
    eta = _finite_scalar(preflight.eta, "eta")
    s_q = _finite_scalar(preflight.s_Q, "s_Q")
    c_patch = math.exp(log_c_patch)
    a_eta = (2.0 ** (0.5 + lam)) * c_patch / (1.0 + eta * eta)
    V = a_eta * (s_q ** (lam - h)) * np.power(R, -1.0 - 2.0 * lam)
    G = np.zeros_like(R)
    if np.any(~np.isfinite(V)) or np.any(V <= 0.0):
        raise CurrentI4RF30TypedDefectError("source RF43 base V became non-finite/non-positive")
    return V, G, lam, h, log_c_patch


def _typed_state(
    backend: ExactCurrentI4NonlinearBackend,
    preflight: CurrentI4FixedQObservablePreflight,
    haar: AuxiliaryT2HaarCovarianceReceipt,
    provider: SourceAuxiliaryT2WaveProvider,
) -> tuple[CandidateIdentity, RF30PreMeanState, str, str, float, float, float]:
    if haar.fixed_q_preflight_sha256 != preflight.preflight_sha256:
        raise CurrentI4RF30TypedDefectError("#1172 Haar receipt detached from fixed-Q preflight")
    if haar.candidate_semantic_sha256 != preflight.candidate_semantic_sha256:
        raise CurrentI4RF30TypedDefectError("#1172 Haar candidate identity drifted")
    if haar.leading_semantic_sha256 != preflight.leading_semantic_sha256:
        raise CurrentI4RF30TypedDefectError("#1172 Haar leading identity drifted")
    if tuple(haar.radius_R) != tuple(preflight.radius_R):
        raise CurrentI4RF30TypedDefectError("#1172 Haar R grid drifted")

    V, G, lam, h, log_c_patch = _exact_i4_background(backend, preflight)
    chart_payload: dict[str, object] = {
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
    chart_sha = _canonical_sha(chart_payload)
    chart_id = f"current-I4-fixed-Q-{chart_sha[:16]}"

    authorized = bool(haar.rf30_repository_candidate_state_authorized)
    identity = CandidateIdentity(
        candidate_id=f"current-I4-A2-{AGENT2_COMPOSITE_PR}-A1-{AGENT1_LEADING_PR}",
        candidate_sha256=preflight.candidate_semantic_sha256,
        evidence_kind="repository-candidate" if authorized else "mechanics-only",
    )
    metadata = getattr(provider, "metadata", None)
    actual_candidate_recomputed = bool(
        getattr(metadata, "provider_kind", None) == "repository_candidate"
        and getattr(metadata, "recomputed_from_actual_candidate", False)
    )

    state = RF30PreMeanState(
        source_candidate_id=identity.candidate_id,
        source_candidate_sha256=identity.candidate_sha256,
        source_chart_id=chart_id,
        source_chart_sha256=chart_sha,
        radius_R=tuple(float(v) for v in preflight.radius_R),
        base_V=tuple(float(v) for v in V),
        base_G=tuple(float(v) for v in G),
        mean_W_rr=haar.mean_W_rr_auxiliary_haar,
        mean_W_zr=haar.mean_W_zr_auxiliary_haar,
        mean_W_thetatheta=haar.mean_W_thetatheta_auxiliary_haar,
        mean_W_ztheta=haar.mean_W_ztheta_auxiliary_haar,
        mean_W_zz=haar.mean_W_zz_auxiliary_haar,
        dR_mean_W_rr=haar.dR_mean_W_rr_auxiliary_haar,
        dZ_mean_W_zr=haar.dZ_mean_W_zr_auxiliary_haar,
        actual_candidate_recomputed=actual_candidate_recomputed,
        oscillatory_covariance_recomputed=True,
        normalized_haar_mean_used=haar.normalized_source_auxiliary_t2_haar_mean_used,
        source_fixed_q_chart_used=True,
        surrogate_defect_used=False,
        residual_as_forcing_shortcut_used=False,
        heldout_samples_used_to_construct_state=False,
    )
    return identity, state, chart_id, chart_sha, lam, h, log_c_patch


def materialize_current_i4_rf30_typed_defect(
    backend: ExactCurrentI4NonlinearBackend,
    provider: SourceAuxiliaryT2WaveProvider,
    radius: Any,
    z: Any,
    t: Any,
) -> CurrentI4RF30TypedDefectReceipt:
    """Materialize current-I4 RF30 state/defect from raw auxiliary-T2 waves only."""
    if not isinstance(backend, ExactCurrentI4NonlinearBackend):
        raise TypeError("backend must be ExactCurrentI4NonlinearBackend")

    preflight = materialize_current_i4_rf30_fixedq_preflight(backend, radius, z, t)
    haar = materialize_current_i4_rf30_auxiliary_haar_bridge(
        backend, provider, radius, z, t
    )
    identity, state, chart_id, chart_sha, lam, h, log_c_patch = _typed_state(
        backend, preflight, haar, provider
    )
    defect, _, V, G = _rf30_defect(state)

    repository_evidence = bool(
        identity.evidence_kind == "repository-candidate"
        and haar.rf30_repository_candidate_state_authorized
        and state.actual_candidate_recomputed
        and state.normalized_haar_mean_used
        and state.source_fixed_q_chart_used
    )
    payload: dict[str, object] = {
        "schema": SCHEMA,
        "task": TASK,
        "candidate": asdict(identity),
        "source_chart_id": chart_id,
        "source_chart_sha256": chart_sha,
        "fixed_q_preflight_sha256": preflight.preflight_sha256,
        "auxiliary_haar_receipt_sha256": haar.receipt_sha256,
        "provider_semantic_sha256": haar.provider_semantic_sha256,
        "provider_source_blob_sha1": haar.provider_source_blob_sha1,
        "provider_kind": haar.provider_kind,
        "state": asdict(state),
        "defect": asdict(defect),
        "lambda_value": lam,
        "h_value": h,
        "log_c_patch": log_c_patch,
        "repository_candidate_defect_evidence": repository_evidence,
    }
    receipt_sha = _canonical_sha(payload)
    return CurrentI4RF30TypedDefectReceipt(
        candidate=identity,
        source_chart_id=chart_id,
        source_chart_sha256=chart_sha,
        fixed_q_preflight_sha256=preflight.preflight_sha256,
        auxiliary_haar_receipt_sha256=haar.receipt_sha256,
        provider_semantic_sha256=haar.provider_semantic_sha256,
        provider_source_blob_sha1=haar.provider_source_blob_sha1,
        provider_kind=haar.provider_kind,
        state=state,
        defect=defect,
        lambda_value=lam,
        h_value=h,
        log_c_patch=log_c_patch,
        base_V_min_abs=float(np.min(np.abs(V))),
        base_V_max_abs=float(np.max(np.abs(V))),
        base_G_max_abs=float(np.max(np.abs(G))),
        source_background_from_exact_a1_i4_schedule=True,
        normalized_source_auxiliary_t2_haar_mean_used=True,
        typed_rf30_state_materialized=True,
        rf30_defect_materialized=True,
        repository_candidate_defect_evidence=repository_evidence,
        rf31_five_row_system_materialized=False,
        correction_applied=False,
        heldout_ns_residual_assessed=False,
        pde_validated=False,
        receipt_sha256=receipt_sha,
    )


def truth_boundary() -> dict[str, object]:
    signature = inspect.signature(materialize_current_i4_rf30_typed_defect)
    forbidden = {
        "P", "J_theta", "J_z", "defect", "residual", "covariance", "forcing",
        "pressure", "gain", "damping", "threshold", "viscosity", "nu",
        "correction", "heldout", "base_V", "base_G", "lambda_value", "c_patch",
    }
    return {
        "schema": SCHEMA,
        "task": TASK,
        "parent_agent3_pr": PARENT_AGENT3_PR,
        "parent_agent3_head": PARENT_AGENT3_HEAD,
        "parent_haar_bridge_blob": PARENT_HAAR_BRIDGE_BLOB,
        "rf30_typed_parent_pr": RF30_TYPED_PARENT_PR,
        "rf30_typed_parent_head": RF30_TYPED_PARENT_HEAD,
        "rf30_typed_parent_blob": RF30_TYPED_PARENT_BLOB,
        "agent1_i4_pr": AGENT1_LEADING_PR,
        "agent1_i4_head": AGENT1_LEADING_HEAD,
        "agent1_i4_source_blob": AGENT1_I4_SOURCE_BLOB,
        "source_reader_head": SOURCE_READER_HEAD,
        "public_parameters": tuple(signature.parameters),
        "forbidden_scientific_controls_exposed": bool(forbidden.intersection(signature.parameters)),
        "normalized_auxiliary_t2_haar_bridge_consumed": True,
        "physical_theta_mean_promoted_to_source_auxiliary_haar": False,
        "source_i4_background_derived_from_exact_a1_schedule": True,
        "caller_supplied_base_background_allowed": False,
        "caller_supplied_rf30_defect_allowed": False,
        "typed_rf30_state_adapter_executable": True,
        "rf30_defect_operator_reused_from_1126": True,
        "agent2_curl_or_jacobian_reimplemented": False,
        "repository_provider_blob_pinned_in_parent": PINNED_REPOSITORY_PROVIDER_BLOB is not None,
        "current_i4_repository_candidate_rf30_defect_materialized": False,
        "rf31_five_row_system_materialized": False,
        "correction_applied": False,
        "finite_correction_cycle_run": False,
        "heldout_ns_residual_assessed": False,
        "residual_reduction_claimed": False,
        "final_normalized_momentum_gate": FINAL_NORMALIZED_MOMENTUM_GATE,
        "final_normalized_divergence_gate": FINAL_NORMALIZED_DIVERGENCE_GATE,
        "residual_defined_free_forcing_allowed": False,
        "finite_stage_small_residual_is_blowup_proof": False,
        "paper_exact": False,
        "openai_field_identified": False,
        "pde_validated": False,
    }


__all__ = [
    "CurrentI4RF30TypedDefectError",
    "CurrentI4RF30TypedDefectReceipt",
    "materialize_current_i4_rf30_typed_defect",
    "source_provenance",
    "truth_boundary",
]
