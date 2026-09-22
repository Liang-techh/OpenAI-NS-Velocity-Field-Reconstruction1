"""Current-I4 handoff from the compact mean correction to RF44--RF49 recomputation.

This module is a thin Kokuno Agent-3 adapter.  It does not introduce a new
mean-defect formula, correction solve, curl/Jacobian, or gain criterion.
Instead it composes two already-delivered operators on one exact current-I4
identity:

* #1190: current-I4 RF30 -> RF31 -> RF34--RF39 compact chart correction;
* #1142: typed RF44--RF49 post-update nonlinear mean-defect recomputation.

The same provider object must supply the raw auxiliary-T2 wave interface used by
#1190 and candidate-bound RF44 pre/post states plus the compact correction
increment required by #1142.  #1142 independently verifies the pre/post state,
normalized Haar/source fixed-Q provenance, exact RF34--RF39 delta-v/delta-gamma
match, solenoidal delta-beta provenance, and RF47 post-update identities.

At the current stack the Agent-2 same-identity repository auxiliary-T2 provider
is still not checksum-pinned.  Therefore this handoff can exercise mechanics
but cannot promote repository-candidate correction evidence or an NS residual
claim.  No held-out sample, pressure, forcing, residual, gain, damping, or
scientific threshold is accepted by the public API.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import inspect
import json
from typing import Any, Protocol

from .kokuno_current_i4_nonlinear_mean_attribution import ExactCurrentI4NonlinearBackend
from .kokuno_current_i4_rf30_fixedq_preflight import materialize_current_i4_rf30_fixedq_preflight
from .kokuno_current_i4_rf30_typed_defect import materialize_current_i4_rf30_typed_defect
from .kokuno_current_i4_rf34_rf39_correction import (
    CurrentI4RF34RF39CorrectionReceipt,
    _FrozenCurrentI4RF34RF39Backend,
    _bind_typed_state,
    _freeze_patch_context_and_basis,
    materialize_current_i4_rf34_rf39_correction,
)
from .kokuno_rf30_auxiliary_t2_haar_bridge import (
    PINNED_REPOSITORY_PROVIDER_BLOB,
    SourceAuxiliaryT2WaveProvider,
)
from .kokuno_rf34_rf39_compact_mean_correction import RF34RF39CorrectionReceipt
from .kokuno_rf44_rf49_postupdate_recompute import (
    RF44CorrectionIncrement,
    RF44MeanState,
    RF44RF49RecomputeReceipt,
    _sha as _rf44_sha,
    materialize_rf44_rf49_postupdate_recompute,
)
from .kokuno_rf30_rf31_typed_mean_correction import CandidateIdentity

TASK = "KOKUNO-A3-CURRENT-I4-RF44-RF49-HANDOFF-130"
SCHEMA = "kokuno-a3-current-i4-rf44-rf49-handoff-v1"

PARENT_AGENT3_PR = 1190
PARENT_AGENT3_HEAD = "21b0ade966fda3b7344273c7bff2efe7e7b9a93e"
PARENT_CURRENT_I4_SOURCE_BLOB = "fe08f93f48aeb2d4055b12889325a6ae5e05a5a6"
RF44_RF49_PR = 1142
RF44_RF49_SOURCE_BLOB = "f6dac0be5cf88f39ca94a6066cdbf1f9270293fb"

SOURCE_READER_REPO = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_READER_HEAD = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_READER_DATE = "2026-09-09"
SOURCE_READER_PATH = "navier-stokes/navier_stokes_workbench.tex"
SOURCE_READER_BLOB = "205a99807302e21a51c5eaf223390c0dfc42bcd0"
SOURCE_FORMULAS = ("RF30", "RF31", "RF34", "RF39", "RF44", "RF45", "RF46", "RF47", "RF48", "RF49")

FINAL_NORMALIZED_MOMENTUM_GATE = 1.0e-3
FINAL_NORMALIZED_DIVERGENCE_GATE = 1.0e-5


class CurrentI4RF44RF49HandoffError(RuntimeError):
    """Raised when the current-I4 RF44/RF49 identity handoff fails closed."""


class CurrentI4RF44StateProvider(SourceAuxiliaryT2WaveProvider, Protocol):
    """Same provider identity used for raw Haar waves and RF44 pre/post states."""

    def rf44_pre_update_state(
        self,
        *,
        current: CurrentI4RF34RF39CorrectionReceipt,
        correction: RF34RF39CorrectionReceipt,
    ) -> RF44MeanState: ...

    def rf44_post_update_state(
        self,
        *,
        current: CurrentI4RF34RF39CorrectionReceipt,
        correction: RF34RF39CorrectionReceipt,
    ) -> RF44MeanState: ...

    def rf44_correction_increment(
        self,
        *,
        current: CurrentI4RF34RF39CorrectionReceipt,
        correction: RF34RF39CorrectionReceipt,
        pre_state: RF44MeanState,
        post_state: RF44MeanState,
    ) -> RF44CorrectionIncrement: ...


@dataclass(frozen=True)
class CurrentI4RF44RF49HandoffReceipt:
    candidate: CandidateIdentity
    current_i4_rf34_rf39_receipt_sha256: str
    source_chart_id: str
    source_chart_sha256: str
    provider_semantic_sha256: str
    rf44_rf49: RF44RF49RecomputeReceipt
    exact_rf34_rf39_replay_verified: bool
    rf44_state_provider_same_object_as_auxiliary_wave_provider: bool
    repository_provider_blob_pinned: bool
    rf44_rf49_operator_recompute_executed: bool
    current_i4_repository_candidate_rf44_rf49_remainder_recomputed: bool
    chart_correction_applied_in_recompute: bool
    cartesian_correction_velocity_materialized: bool
    finite_correction_cycle_run: bool
    heldout_ns_residual_assessed: bool
    residual_reduction_claimed: bool
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
        "rf44_rf49_formula_implementation": "reused_exactly_from_A3_PR_1142",
        "current_i4_correction_implementation": "reused_exactly_from_A3_PR_1190",
        "classification": "public_structure_provenance_not_independent_validation",
    }


def _require_state_provider(provider: Any) -> None:
    for name in (
        "rf44_pre_update_state",
        "rf44_post_update_state",
        "rf44_correction_increment",
    ):
        if not callable(getattr(provider, name, None)):
            raise CurrentI4RF44RF49HandoffError(
                f"provider lacks typed current-I4 RF44 method {name}"
            )


class _CurrentI4RF44Backend(_FrozenCurrentI4RF34RF39Backend):
    """Adapter consumed by the existing #1142 recomputation operator."""

    def __init__(self, *, identity, state, basis, context, provider, current) -> None:
        super().__init__(identity=identity, state=state, basis=basis, context=context)
        self._provider = provider
        self._current = current

    def rf44_pre_update_state(
        self, candidate: Any, correction: RF34RF39CorrectionReceipt
    ) -> RF44MeanState:
        return self._provider.rf44_pre_update_state(
            current=self._current, correction=correction
        )

    def rf44_post_update_state(
        self, candidate: Any, correction: RF34RF39CorrectionReceipt
    ) -> RF44MeanState:
        return self._provider.rf44_post_update_state(
            current=self._current, correction=correction
        )

    def rf44_correction_increment(
        self,
        candidate: Any,
        correction: RF34RF39CorrectionReceipt,
        pre_state: RF44MeanState,
        post_state: RF44MeanState,
    ) -> RF44CorrectionIncrement:
        return self._provider.rf44_correction_increment(
            current=self._current,
            correction=correction,
            pre_state=pre_state,
            post_state=post_state,
        )


def materialize_current_i4_rf44_rf49_handoff(
    backend: ExactCurrentI4NonlinearBackend,
    provider: CurrentI4RF44StateProvider,
    radius: Any,
    z: Any,
    t: Any,
) -> CurrentI4RF44RF49HandoffReceipt:
    """Route one exact current-I4 compact correction through existing RF44--RF49."""
    if not isinstance(backend, ExactCurrentI4NonlinearBackend):
        raise TypeError("backend must be ExactCurrentI4NonlinearBackend")
    _require_state_provider(provider)

    # First materialize the public #1190 current-I4 correction receipt.
    current = materialize_current_i4_rf34_rf39_correction(
        backend, provider, radius, z, t
    )

    # Rebuild only the typed inputs required by #1142 so that #1142 itself
    # independently rematerializes the exact same #1134 RF34--RF39 correction.
    # No RF44/RF47 formula is reproduced in this adapter.
    preflight = materialize_current_i4_rf30_fixedq_preflight(backend, radius, z, t)
    provisional, context, basis = _freeze_patch_context_and_basis(backend, preflight)
    typed = materialize_current_i4_rf30_typed_defect(
        backend, provider, radius, z, t
    )
    identity, basis, context = _bind_typed_state(
        provisional, basis, context, typed, preflight
    )

    if identity != current.candidate:
        raise CurrentI4RF44RF49HandoffError(
            "current-I4 identity changed between correction and RF44 handoff"
        )
    if typed.source_chart_id != current.source_chart_id:
        raise CurrentI4RF44RF49HandoffError("current-I4 source chart id changed")
    if typed.source_chart_sha256 != current.source_chart_sha256:
        raise CurrentI4RF44RF49HandoffError("current-I4 source chart hash changed")

    adapter = _CurrentI4RF44Backend(
        identity=identity,
        state=typed.state,
        basis=basis,
        context=context,
        provider=provider,
        current=current,
    )
    recomputed = materialize_rf44_rf49_postupdate_recompute(
        adapter, typed.receipt_sha256
    )

    expected_correction_sha = _rf44_sha(asdict(current.correction))
    if recomputed.parent_correction_sha256 != expected_correction_sha:
        raise CurrentI4RF44RF49HandoffError(
            "RF44/RF49 replay did not use the exact #1190 RF34--RF39 correction"
        )
    if recomputed.candidate != current.candidate:
        raise CurrentI4RF44RF49HandoffError("RF44/RF49 candidate identity drifted")
    if recomputed.source_chart_id != current.source_chart_id:
        raise CurrentI4RF44RF49HandoffError("RF44/RF49 source chart id drifted")
    if recomputed.source_chart_sha256 != current.source_chart_sha256:
        raise CurrentI4RF44RF49HandoffError("RF44/RF49 source chart hash drifted")

    pinned = PINNED_REPOSITORY_PROVIDER_BLOB is not None
    repository_evidence = bool(
        pinned
        and current.repository_candidate_correction_evidence
        and recomputed.repository_candidate_recompute_evidence
        and recomputed.candidate.evidence_kind == "repository-candidate"
    )
    provider_sha = str(provider.metadata.provider_semantic_sha256)
    payload: dict[str, object] = {
        "schema": SCHEMA,
        "task": TASK,
        "candidate": asdict(current.candidate),
        "current_i4_rf34_rf39_receipt_sha256": current.receipt_sha256,
        "source_chart_id": current.source_chart_id,
        "source_chart_sha256": current.source_chart_sha256,
        "provider_semantic_sha256": provider_sha,
        "rf44_parent_correction_sha256": recomputed.parent_correction_sha256,
        "rf44_pre_state_sha256": recomputed.pre_defect.state_sha256,
        "rf44_post_state_sha256": recomputed.post_defect.state_sha256,
        "rf44_defect_gain_ratio_2": recomputed.defect_gain_ratio_2,
        "repository_provider_blob_pinned": pinned,
        "repository_candidate_recompute_evidence": repository_evidence,
    }
    receipt_sha = _canonical_sha(payload)

    return CurrentI4RF44RF49HandoffReceipt(
        candidate=current.candidate,
        current_i4_rf34_rf39_receipt_sha256=current.receipt_sha256,
        source_chart_id=current.source_chart_id,
        source_chart_sha256=current.source_chart_sha256,
        provider_semantic_sha256=provider_sha,
        rf44_rf49=recomputed,
        exact_rf34_rf39_replay_verified=True,
        rf44_state_provider_same_object_as_auxiliary_wave_provider=True,
        repository_provider_blob_pinned=pinned,
        rf44_rf49_operator_recompute_executed=True,
        current_i4_repository_candidate_rf44_rf49_remainder_recomputed=repository_evidence,
        chart_correction_applied_in_recompute=True,
        cartesian_correction_velocity_materialized=False,
        finite_correction_cycle_run=False,
        heldout_ns_residual_assessed=False,
        residual_reduction_claimed=False,
        pde_validated=False,
        receipt_sha256=receipt_sha,
    )


def truth_boundary() -> dict[str, object]:
    signature = inspect.signature(materialize_current_i4_rf44_rf49_handoff)
    forbidden = {
        "P", "J_theta", "J_z", "defect", "residual", "covariance", "basis",
        "delta_beta", "delta_v", "delta_gamma", "coefficient", "correction",
        "pressure", "forcing", "gain", "damping", "threshold", "viscosity",
        "nu", "heldout", "lambda_value", "c_patch",
    }
    return {
        "schema": SCHEMA,
        "task": TASK,
        "parent_agent3_pr": PARENT_AGENT3_PR,
        "parent_agent3_head": PARENT_AGENT3_HEAD,
        "parent_current_i4_source_blob": PARENT_CURRENT_I4_SOURCE_BLOB,
        "rf44_rf49_pr": RF44_RF49_PR,
        "rf44_rf49_source_blob": RF44_RF49_SOURCE_BLOB,
        "source_reader_head": SOURCE_READER_HEAD,
        "source_formulas": SOURCE_FORMULAS,
        "public_parameters": tuple(signature.parameters),
        "forbidden_scientific_controls_exposed": bool(
            forbidden.intersection(signature.parameters)
        ),
        "current_i4_rf34_rf39_correction_consumed": True,
        "existing_rf44_rf49_operator_reused_without_formula_duplication": True,
        "same_provider_object_required_for_auxiliary_t2_and_rf44_states": True,
        "caller_supplied_defect_allowed": False,
        "caller_supplied_correction_allowed": False,
        "agent2_curl_or_jacobian_reimplemented": False,
        "repository_provider_blob_pinned_in_parent": (
            PINNED_REPOSITORY_PROVIDER_BLOB is not None
        ),
        "rf44_rf49_operator_recompute_executable": True,
        "current_i4_repository_candidate_rf44_rf49_remainder_recomputed": False,
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
    "CurrentI4RF44RF49HandoffError",
    "CurrentI4RF44StateProvider",
    "CurrentI4RF44RF49HandoffReceipt",
    "materialize_current_i4_rf44_rf49_handoff",
    "source_provenance",
    "truth_boundary",
]
