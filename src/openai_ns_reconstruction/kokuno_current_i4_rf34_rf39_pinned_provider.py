"""Promote the exact pinned current-I4 RF30 state into RF31/RF34--RF39 evidence.

This Kokuno Agent-3 adapter closes one provenance seam only.  PR #1207 already
pins the exact same-identity Agent-2 #1198 auxiliary-T2 provider and proves that
its normalized-Haar RF30 defect is repository-candidate evidence.  PR #1190
already freezes the source-I4 correction patch and RF31 basis before observing
that defect, then reuses the source-specific #1134 RF34--RF39 compact correction.

The present module executes both parents on the same backend/provider/probes and
requires them to share the *same deterministic typed-RF30 receipt SHA*, candidate
identity and fixed-Q chart identity.  Only then is the already-computed RF31 /
RF34--RF39 correction promoted from mechanics provenance to repository-candidate
correction evidence.

No RF30/RF31/RF34--RF39 formula is reimplemented here.  No caller may supply a
defect, covariance, P/J target, basis, coefficient, correction array, pressure,
forcing, gain, held-out sample or scientific threshold.  The correction is not
applied to Cartesian velocity in this increment; RF44--RF49 recomputation,
Cartesian delta-u, a finite correction cycle and held-out complete NS residual
remain separate downstream obligations.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import inspect
import json
import math
from typing import Any

from .kokuno_current_i4_nonlinear_mean_attribution import ExactCurrentI4NonlinearBackend
from .kokuno_current_i4_rf30_pinned_provider import (
    PINNED_AGENT2_HEAD,
    PINNED_AGENT2_PR,
    PINNED_AGENT2_PROVIDER_BLOB,
    CurrentI4RF30PinnedProviderReceipt,
    materialize_current_i4_rf30_pinned_provider_evidence,
)
from .kokuno_current_i4_rf34_rf39_correction import (
    CurrentI4RF34RF39CorrectionReceipt,
    materialize_current_i4_rf34_rf39_correction,
)
from .kokuno_rf30_auxiliary_t2_haar_bridge import SourceAuxiliaryT2WaveProvider

TASK = "KOKUNO-A3-CURRENT-I4-PINNED-RF34-RF39-147"
SCHEMA = "kokuno-a3-current-i4-pinned-rf34-rf39-v1"

PARENT_PIN_PR = 1207
PARENT_PIN_HEAD = "5fc495d322edd0c683834f888d8185ae9bcae368"
PARENT_PIN_BLOB = "82009711a2e80dfb5afbd2b262f3258e0b2f4ac6"
PARENT_RF34_RF39_PR = 1190
PARENT_RF34_RF39_HEAD = "21b0ade966fda3b7344273c7bff2efe7e7b9a93e"
PARENT_RF34_RF39_BLOB = "fe08f93f48aeb2d4055b12889325a6ae5e05a5a6"
SOURCE_RF34_RF39_PR = 1134
SOURCE_RF34_RF39_BLOB = "392975eafd3787c91e2dc82ad9b735567513271d"

SOURCE_READER_REPO = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_READER_HEAD = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_READER_DATE = "2026-09-09"
SOURCE_READER_PATH = "navier-stokes/navier_stokes_workbench.tex"
SOURCE_READER_BLOB = "205a99807302e21a51c5eaf223390c0dfc42bcd0"
SOURCE_FORMULAS = ("RF30", "RF31", "RF34", "RF35", "RF36", "RF37", "RF38", "RF39", "RF43")

FINAL_NORMALIZED_MOMENTUM_GATE = 1.0e-3
FINAL_NORMALIZED_DIVERGENCE_GATE = 1.0e-5


class CurrentI4PinnedRF34RF39Error(RuntimeError):
    """Raised when the pinned-RF30 to compact-correction handoff drifts."""


def _canonical_sha(payload: dict[str, object]) -> str:
    raw = json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)
    return hashlib.sha256(raw.encode("utf-8")).hexdigest()


@dataclass(frozen=True)
class CurrentI4PinnedRF34RF39Receipt:
    candidate_id: str
    candidate_sha256: str
    source_chart_id: str
    source_chart_sha256: str
    provider_semantic_sha256: str
    provider_source_blob_sha1: str
    pinned_rf30_receipt_sha256: str
    shared_typed_rf30_receipt_sha256: str
    parent_rf34_rf39_receipt_sha256: str
    rf31_system_sha256: str
    correction_payload_sha256: str
    coefficients_u0_u1_u2_s0_s1: tuple[float, float, float, float, float]
    correction_chart_l2: float
    correction_chart_max: float
    five_row_closure_max_relative: float
    coefficient_direct_solve_max_relative: float
    exact_provider_blob_pinned: bool
    same_typed_rf30_receipt_preserved: bool
    same_candidate_identity_preserved: bool
    same_source_chart_identity_preserved: bool
    basis_frozen_before_defect_evaluation: bool
    patch_context_frozen_before_defect_evaluation: bool
    compact_support_preserved: bool
    two_zero_moments_preserved: bool
    repository_candidate_rf31_system_evidence: bool
    repository_candidate_rf34_rf39_correction_evidence: bool
    correction_applied_to_candidate: bool
    rf44_rf49_nonlinear_remainder_recomputed: bool
    cartesian_correction_velocity_materialized: bool
    finite_correction_cycle_run: bool
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
        "source_formulas": SOURCE_FORMULAS,
        "agent2_provider": {
            "pr": PINNED_AGENT2_PR,
            "head": PINNED_AGENT2_HEAD,
            "blob_sha1": PINNED_AGENT2_PROVIDER_BLOB,
            "classification": "repository_autonomous_same_identity_auxiliary_T2_candidate_lift_not_source_exact_mode_recovery",
        },
        "compact_bump_realization": "autonomous_repository_choice_reused_from_PR_1134",
        "classification": "public_structure_provenance_not_independent_validation",
    }


def _require_parent_alignment(
    pinned: CurrentI4RF30PinnedProviderReceipt,
    parent: CurrentI4RF34RF39CorrectionReceipt,
) -> None:
    if not pinned.repository_candidate_rf30_defect_evidence:
        raise CurrentI4PinnedRF34RF39Error("pinned parent did not authorize repository RF30 evidence")
    if not pinned.exact_provider_blob_pinned:
        raise CurrentI4PinnedRF34RF39Error("exact Agent-2 provider blob is not pinned")
    if pinned.provider_source_blob_sha1 != PINNED_AGENT2_PROVIDER_BLOB:
        raise CurrentI4PinnedRF34RF39Error("pinned RF30 provider blob drifted")

    if parent.typed_rf30_receipt_sha256 != pinned.parent_typed_rf30_receipt_sha256:
        raise CurrentI4PinnedRF34RF39Error(
            "RF34-RF39 path did not consume the exact typed RF30 receipt authorized by #1207"
        )
    if parent.candidate.candidate_id != pinned.candidate_id:
        raise CurrentI4PinnedRF34RF39Error("candidate id drifted between pinned RF30 and correction")
    if parent.candidate.candidate_sha256 != pinned.candidate_sha256:
        raise CurrentI4PinnedRF34RF39Error("candidate hash drifted between pinned RF30 and correction")
    if parent.source_chart_id != pinned.source_chart_id:
        raise CurrentI4PinnedRF34RF39Error("source chart id drifted between RF30 and correction")
    if parent.source_chart_sha256 != pinned.source_chart_sha256:
        raise CurrentI4PinnedRF34RF39Error("source chart hash drifted between RF30 and correction")

    if not parent.basis_frozen_before_defect_evaluation:
        raise CurrentI4PinnedRF34RF39Error("RF31 basis was not frozen before defect evaluation")
    if not parent.patch_context_frozen_before_defect_evaluation:
        raise CurrentI4PinnedRF34RF39Error("RF34-RF39 patch context was not frozen before defect evaluation")
    if not parent.current_i4_rf31_system_materialized:
        raise CurrentI4PinnedRF34RF39Error("current-I4 RF31 system is not materialized")
    if not parent.current_i4_rf34_rf39_correction_materialized:
        raise CurrentI4PinnedRF34RF39Error("current-I4 RF34-RF39 correction is not materialized")

    correction = parent.correction
    if not correction.source_rf34_rf39_correction_materialized:
        raise CurrentI4PinnedRF34RF39Error("source RF34-RF39 correction receipt is incomplete")
    if not correction.compact_support_preserved:
        raise CurrentI4PinnedRF34RF39Error("compact-support condition failed")
    if not correction.two_zero_moments_preserved:
        raise CurrentI4PinnedRF34RF39Error("RF31 zero-moment conditions failed")

    values = (
        *correction.coefficients_u0_u1_u2_s0_s1,
        correction.correction_chart_l2,
        correction.correction_chart_max,
        correction.five_row_closure_max_relative,
        correction.coefficient_direct_solve_max_relative,
    )
    if not all(math.isfinite(float(v)) for v in values):
        raise CurrentI4PinnedRF34RF39Error("compact correction diagnostics became non-finite")

    # These downstream states must remain false in this one-increment adapter.
    if parent.correction_applied_to_candidate or correction.correction_applied_to_candidate:
        raise CurrentI4PinnedRF34RF39Error("parent unexpectedly applies correction to candidate")
    if parent.rf44_rf49_nonlinear_remainder_recomputed or correction.rf44_rf49_nonlinear_remainder_recomputed:
        raise CurrentI4PinnedRF34RF39Error("parent unexpectedly promotes RF44-RF49 recomputation")
    if parent.cartesian_correction_velocity_materialized or correction.cartesian_correction_velocity_materialized:
        raise CurrentI4PinnedRF34RF39Error("parent unexpectedly materializes Cartesian correction")
    if parent.heldout_ns_residual_assessed or correction.heldout_ns_residual_assessed:
        raise CurrentI4PinnedRF34RF39Error("parent unexpectedly assesses held-out NS residual")
    if parent.pde_validated or correction.pde_validated:
        raise CurrentI4PinnedRF34RF39Error("parent unexpectedly promotes PDE validation")


def materialize_current_i4_rf34_rf39_pinned_provider_correction(
    backend: ExactCurrentI4NonlinearBackend,
    provider: SourceAuxiliaryT2WaveProvider,
    radius: Any,
    z: Any,
    t: Any,
) -> CurrentI4PinnedRF34RF39Receipt:
    """Promote one exact pinned RF30 receipt through the existing compact correction."""
    pinned = materialize_current_i4_rf30_pinned_provider_evidence(
        backend, provider, radius, z, t
    )
    parent = materialize_current_i4_rf34_rf39_correction(
        backend, provider, radius, z, t
    )
    _require_parent_alignment(pinned, parent)

    correction = parent.correction
    correction_payload = {
        "parent_system_sha256": correction.parent_system_sha256,
        "source_chart_sha256": correction.source_chart_sha256,
        "patch_context_sha256": correction.patch_context_sha256,
        "autonomous_bump_family_sha256": correction.autonomous_bump_family_sha256,
        "coefficients": correction.coefficients_u0_u1_u2_s0_s1,
        "delta_v_chart": correction.delta_v_chart,
        "gamma_d_chart": correction.gamma_d_chart,
    }
    correction_payload_sha = _canonical_sha(correction_payload)

    payload: dict[str, object] = {
        "schema": SCHEMA,
        "task": TASK,
        "candidate_id": pinned.candidate_id,
        "candidate_sha256": pinned.candidate_sha256,
        "source_chart_id": pinned.source_chart_id,
        "source_chart_sha256": pinned.source_chart_sha256,
        "provider_semantic_sha256": pinned.provider_semantic_sha256,
        "provider_source_blob_sha1": pinned.provider_source_blob_sha1,
        "pinned_rf30_receipt_sha256": pinned.receipt_sha256,
        "shared_typed_rf30_receipt_sha256": pinned.parent_typed_rf30_receipt_sha256,
        "parent_rf34_rf39_receipt_sha256": parent.receipt_sha256,
        "rf31_system_sha256": correction.parent_system_sha256,
        "correction_payload_sha256": correction_payload_sha,
        "repository_candidate_rf31_system_evidence": True,
        "repository_candidate_rf34_rf39_correction_evidence": True,
    }
    receipt_sha = _canonical_sha(payload)

    return CurrentI4PinnedRF34RF39Receipt(
        candidate_id=pinned.candidate_id,
        candidate_sha256=pinned.candidate_sha256,
        source_chart_id=pinned.source_chart_id,
        source_chart_sha256=pinned.source_chart_sha256,
        provider_semantic_sha256=pinned.provider_semantic_sha256,
        provider_source_blob_sha1=pinned.provider_source_blob_sha1,
        pinned_rf30_receipt_sha256=pinned.receipt_sha256,
        shared_typed_rf30_receipt_sha256=pinned.parent_typed_rf30_receipt_sha256,
        parent_rf34_rf39_receipt_sha256=parent.receipt_sha256,
        rf31_system_sha256=correction.parent_system_sha256,
        correction_payload_sha256=correction_payload_sha,
        coefficients_u0_u1_u2_s0_s1=tuple(
            float(v) for v in correction.coefficients_u0_u1_u2_s0_s1
        ),
        correction_chart_l2=float(correction.correction_chart_l2),
        correction_chart_max=float(correction.correction_chart_max),
        five_row_closure_max_relative=float(correction.five_row_closure_max_relative),
        coefficient_direct_solve_max_relative=float(
            correction.coefficient_direct_solve_max_relative
        ),
        exact_provider_blob_pinned=True,
        same_typed_rf30_receipt_preserved=True,
        same_candidate_identity_preserved=True,
        same_source_chart_identity_preserved=True,
        basis_frozen_before_defect_evaluation=True,
        patch_context_frozen_before_defect_evaluation=True,
        compact_support_preserved=True,
        two_zero_moments_preserved=True,
        repository_candidate_rf31_system_evidence=True,
        repository_candidate_rf34_rf39_correction_evidence=True,
        correction_applied_to_candidate=False,
        rf44_rf49_nonlinear_remainder_recomputed=False,
        cartesian_correction_velocity_materialized=False,
        finite_correction_cycle_run=False,
        heldout_ns_residual_assessed=False,
        pde_validated=False,
        receipt_sha256=receipt_sha,
    )


def truth_boundary() -> dict[str, object]:
    signature = inspect.signature(
        materialize_current_i4_rf34_rf39_pinned_provider_correction
    )
    forbidden = {
        "P", "J_theta", "J_z", "defect", "residual", "covariance", "basis",
        "delta_v", "gamma_d", "coefficient", "bump", "forcing", "pressure",
        "gain", "damping", "threshold", "viscosity", "nu", "correction",
        "heldout", "lambda_value", "c_patch",
    }
    return {
        "schema": SCHEMA,
        "task": TASK,
        "parent_pin_pr": PARENT_PIN_PR,
        "parent_pin_head": PARENT_PIN_HEAD,
        "parent_pin_blob": PARENT_PIN_BLOB,
        "parent_rf34_rf39_pr": PARENT_RF34_RF39_PR,
        "parent_rf34_rf39_head": PARENT_RF34_RF39_HEAD,
        "parent_rf34_rf39_blob": PARENT_RF34_RF39_BLOB,
        "source_rf34_rf39_pr": SOURCE_RF34_RF39_PR,
        "source_rf34_rf39_blob": SOURCE_RF34_RF39_BLOB,
        "public_parameters": tuple(signature.parameters),
        "forbidden_scientific_controls_exposed": bool(
            forbidden.intersection(signature.parameters)
        ),
        "exact_agent2_provider_blob_pin_reused": True,
        "same_typed_rf30_receipt_required_for_promotion": True,
        "rf30_rf31_rf34_rf39_formulas_reused_without_duplication": True,
        "agent2_curl_or_jacobian_reimplemented": False,
        "repository_candidate_rf31_system_evidence_materializable": True,
        "repository_candidate_rf34_rf39_correction_evidence_materializable": True,
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
    "CurrentI4PinnedRF34RF39Error",
    "CurrentI4PinnedRF34RF39Receipt",
    "materialize_current_i4_rf34_rf39_pinned_provider_correction",
    "source_provenance",
    "truth_boundary",
]
