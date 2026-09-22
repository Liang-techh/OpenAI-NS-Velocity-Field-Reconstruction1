"""Checksum-pinned current-I4 RF30 evidence from the exact Agent-2 T2 provider.

This Kokuno Agent-3 adapter is deliberately stacked above the immutable #1172
Haar bridge and #1181 typed RF30 defect adapter.  Those parents remain unchanged.
A2 PR #1198 now supplies an implementation-distinct, same-identity raw
auxiliary-T2 provider.  This module pins its exact Git blob and promotes only
that provider's already-recomputed RF30 state/defect from mechanics provenance
to repository-candidate RF30 evidence.

No covariance, P/J defect, forcing, pressure, correction coefficient, held-out
sample, gain, damping, viscosity, or scientific threshold is caller supplied.
The A2 auxiliary phase lift remains the repository-autonomous candidate lift
stated by A2; this module does not relabel it as source-exact Kokuno mode data.
It also does not execute RF31/RF34-RF39, RF44-RF49, a Cartesian correction, or a
finite correction cycle.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import inspect
import json
import math
from typing import Any

from .kokuno_current_i4_nonlinear_mean_attribution import ExactCurrentI4NonlinearBackend
from .kokuno_current_i4_rf30_typed_defect import (
    CurrentI4RF30TypedDefectReceipt,
    materialize_current_i4_rf30_typed_defect,
)
from .kokuno_rf30_auxiliary_t2_haar_bridge import SourceAuxiliaryT2WaveProvider

TASK = "KOKUNO-A3-CURRENT-I4-AUXT2-PIN-146"
SCHEMA = "kokuno-a3-current-i4-rf30-pinned-provider-v1"
PARENT_AGENT3_PR = 1197
PARENT_AGENT3_HEAD = "aa6b1a77bee404e1e852835fe048a7365b58f379"
PARENT_HAAR_BRIDGE_BLOB = "5b20d4945f3aeb177e327e8d6b5709a6726c885a"
PARENT_TYPED_RF30_BLOB = "9b9704bcd188fe2e03a77e3d39425ea34668b1ed"

PINNED_AGENT2_PR = 1198
PINNED_AGENT2_HEAD = "5765c2b3bab7482df514efa87f9b6bba48e04b7f"
PINNED_AGENT2_PROVIDER_BLOB = "96168ac6583ad5e71aa044bd558c6feeadf051f4"

SOURCE_READER_REPO = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_READER_HEAD = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_READER_DATE = "2026-09-09"
SOURCE_READER_PATH = "navier-stokes/navier_stokes_workbench.tex"
SOURCE_READER_BLOB = "205a99807302e21a51c5eaf223390c0dfc42bcd0"

FINAL_NORMALIZED_MOMENTUM_GATE = 1.0e-3
FINAL_NORMALIZED_DIVERGENCE_GATE = 1.0e-5


class CurrentI4RF30PinnedProviderError(RuntimeError):
    """Raised when exact provider or parent RF30 provenance fails closed."""


def _canonical_sha(payload: dict[str, object]) -> str:
    raw = json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)
    return hashlib.sha256(raw.encode("utf-8")).hexdigest()


@dataclass(frozen=True)
class CurrentI4RF30PinnedProviderReceipt:
    candidate_id: str
    candidate_sha256: str
    source_chart_id: str
    source_chart_sha256: str
    provider_semantic_sha256: str
    provider_source_blob_sha1: str
    pinned_agent2_pr: int
    pinned_agent2_head: str
    pinned_agent2_provider_blob: str
    parent_typed_rf30_receipt_sha256: str
    defect_P: float
    defect_J_theta: float
    defect_J_z: float
    defect_tuple_l2: float
    exact_provider_blob_pinned: bool
    same_candidate_identity_preserved: bool
    actual_candidate_recomputed: bool
    normalized_source_auxiliary_t2_haar_mean_used: bool
    source_fixed_q_chart_used: bool
    surrogate_defect_used: bool
    residual_as_forcing_shortcut_used: bool
    heldout_samples_used_to_construct_state: bool
    repository_candidate_rf30_defect_evidence: bool
    rf31_five_row_system_materialized: bool
    rf34_rf39_correction_materialized: bool
    rf44_rf49_remainder_recomputed: bool
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
        "scope": "RF30 normalized auxiliary-T2 Haar mean and mean-defect moments",
        "agent2_provider": {
            "pr": PINNED_AGENT2_PR,
            "head": PINNED_AGENT2_HEAD,
            "blob_sha1": PINNED_AGENT2_PROVIDER_BLOB,
            "classification": "repository_autonomous_same_identity_auxiliary_T2_candidate_lift_not_source_exact_mode_recovery",
        },
        "classification": "public_structure_provenance_not_independent_validation",
    }


def _validate_exact_provider(provider: SourceAuxiliaryT2WaveProvider) -> Any:
    metadata = getattr(provider, "metadata", None)
    required = (
        "candidate_semantic_sha256",
        "provider_semantic_sha256",
        "source_blob_sha1",
        "provider_kind",
        "source_auxiliary_t2_field_materialized",
        "normalized_haar_measure_total_mass_one",
        "recomputed_from_actual_candidate",
        "surrogate_or_preaveraged_covariance_used",
        "heldout_data_used",
        "residual_as_forcing_used",
    )
    if metadata is None or any(not hasattr(metadata, name) for name in required):
        raise CurrentI4RF30PinnedProviderError("provider metadata contract is incomplete")
    if metadata.provider_kind != "repository_candidate":
        raise CurrentI4RF30PinnedProviderError("provider must be repository_candidate")
    if metadata.source_blob_sha1 != PINNED_AGENT2_PROVIDER_BLOB:
        raise CurrentI4RF30PinnedProviderError("Agent-2 auxiliary-T2 provider blob is not the exact pinned implementation")
    if not metadata.source_auxiliary_t2_field_materialized:
        raise CurrentI4RF30PinnedProviderError("source auxiliary-T2 field is not materialized")
    if not metadata.normalized_haar_measure_total_mass_one:
        raise CurrentI4RF30PinnedProviderError("normalized Haar convention drifted")
    if not metadata.recomputed_from_actual_candidate:
        raise CurrentI4RF30PinnedProviderError("provider is not recomputed from the actual candidate")
    if metadata.surrogate_or_preaveraged_covariance_used:
        raise CurrentI4RF30PinnedProviderError("surrogate/pre-averaged covariance is forbidden")
    if metadata.heldout_data_used:
        raise CurrentI4RF30PinnedProviderError("held-out data cannot construct RF30 evidence")
    if metadata.residual_as_forcing_used:
        raise CurrentI4RF30PinnedProviderError("residual-as-forcing shortcut is forbidden")
    return metadata


def materialize_current_i4_rf30_pinned_provider_evidence(
    backend: ExactCurrentI4NonlinearBackend,
    provider: SourceAuxiliaryT2WaveProvider,
    radius: Any,
    z: Any,
    t: Any,
) -> CurrentI4RF30PinnedProviderReceipt:
    """Promote only the exact A2 #1198 provider's parent RF30 recomputation."""
    if not isinstance(backend, ExactCurrentI4NonlinearBackend):
        raise TypeError("backend must be ExactCurrentI4NonlinearBackend")
    metadata = _validate_exact_provider(provider)
    typed: CurrentI4RF30TypedDefectReceipt = materialize_current_i4_rf30_typed_defect(
        backend, provider, radius, z, t
    )

    if typed.provider_source_blob_sha1 != PINNED_AGENT2_PROVIDER_BLOB:
        raise CurrentI4RF30PinnedProviderError("typed RF30 receipt provider blob drifted")
    if typed.provider_semantic_sha256 != metadata.provider_semantic_sha256:
        raise CurrentI4RF30PinnedProviderError("typed RF30 receipt provider semantic identity drifted")
    if typed.candidate.candidate_sha256 != metadata.candidate_semantic_sha256:
        raise CurrentI4RF30PinnedProviderError("typed RF30 candidate identity drifted from provider")
    state = typed.state
    if not state.actual_candidate_recomputed:
        raise CurrentI4RF30PinnedProviderError("RF30 state is not recomputed from actual candidate")
    if not state.oscillatory_covariance_recomputed:
        raise CurrentI4RF30PinnedProviderError("RF30 covariance was not recomputed")
    if not state.normalized_haar_mean_used or not state.source_fixed_q_chart_used:
        raise CurrentI4RF30PinnedProviderError("RF30 source-chart/Haar contract drifted")
    if state.surrogate_defect_used:
        raise CurrentI4RF30PinnedProviderError("surrogate RF30 defect is forbidden")
    if state.residual_as_forcing_shortcut_used:
        raise CurrentI4RF30PinnedProviderError("residual-as-forcing shortcut is forbidden")
    if state.heldout_samples_used_to_construct_state:
        raise CurrentI4RF30PinnedProviderError("held-out data leaked into RF30 state")
    if not typed.rf30_defect_materialized:
        raise CurrentI4RF30PinnedProviderError("parent RF30 defect was not materialized")
    if typed.repository_candidate_defect_evidence:
        raise CurrentI4RF30PinnedProviderError(
            "parent #1181 unexpectedly self-promoted; exact-stack provenance changed"
        )

    values = (float(typed.defect.P), float(typed.defect.J_theta), float(typed.defect.J_z))
    if not all(math.isfinite(v) for v in values):
        raise CurrentI4RF30PinnedProviderError("RF30 defect tuple became non-finite")
    defect_l2 = math.sqrt(sum(v * v for v in values))

    payload: dict[str, object] = {
        "schema": SCHEMA,
        "task": TASK,
        "candidate_id": typed.candidate.candidate_id,
        "candidate_sha256": typed.candidate.candidate_sha256,
        "source_chart_id": typed.source_chart_id,
        "source_chart_sha256": typed.source_chart_sha256,
        "provider_semantic_sha256": typed.provider_semantic_sha256,
        "provider_source_blob_sha1": typed.provider_source_blob_sha1,
        "pinned_agent2_pr": PINNED_AGENT2_PR,
        "pinned_agent2_head": PINNED_AGENT2_HEAD,
        "pinned_agent2_provider_blob": PINNED_AGENT2_PROVIDER_BLOB,
        "parent_typed_rf30_receipt_sha256": typed.receipt_sha256,
        "defect": values,
        "defect_tuple_l2": defect_l2,
        "repository_candidate_rf30_defect_evidence": True,
    }
    receipt_sha = _canonical_sha(payload)
    return CurrentI4RF30PinnedProviderReceipt(
        candidate_id=typed.candidate.candidate_id,
        candidate_sha256=typed.candidate.candidate_sha256,
        source_chart_id=typed.source_chart_id,
        source_chart_sha256=typed.source_chart_sha256,
        provider_semantic_sha256=typed.provider_semantic_sha256,
        provider_source_blob_sha1=typed.provider_source_blob_sha1,
        pinned_agent2_pr=PINNED_AGENT2_PR,
        pinned_agent2_head=PINNED_AGENT2_HEAD,
        pinned_agent2_provider_blob=PINNED_AGENT2_PROVIDER_BLOB,
        parent_typed_rf30_receipt_sha256=typed.receipt_sha256,
        defect_P=values[0],
        defect_J_theta=values[1],
        defect_J_z=values[2],
        defect_tuple_l2=defect_l2,
        exact_provider_blob_pinned=True,
        same_candidate_identity_preserved=True,
        actual_candidate_recomputed=True,
        normalized_source_auxiliary_t2_haar_mean_used=True,
        source_fixed_q_chart_used=True,
        surrogate_defect_used=False,
        residual_as_forcing_shortcut_used=False,
        heldout_samples_used_to_construct_state=False,
        repository_candidate_rf30_defect_evidence=True,
        rf31_five_row_system_materialized=False,
        rf34_rf39_correction_materialized=False,
        rf44_rf49_remainder_recomputed=False,
        cartesian_correction_velocity_materialized=False,
        finite_correction_cycle_run=False,
        heldout_ns_residual_assessed=False,
        pde_validated=False,
        receipt_sha256=receipt_sha,
    )


def truth_boundary() -> dict[str, object]:
    signature = inspect.signature(materialize_current_i4_rf30_pinned_provider_evidence)
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
        "parent_typed_rf30_blob": PARENT_TYPED_RF30_BLOB,
        "pinned_agent2_pr": PINNED_AGENT2_PR,
        "pinned_agent2_head": PINNED_AGENT2_HEAD,
        "pinned_agent2_provider_blob": PINNED_AGENT2_PROVIDER_BLOB,
        "public_parameters": tuple(signature.parameters),
        "forbidden_scientific_controls_exposed": bool(forbidden.intersection(signature.parameters)),
        "exact_provider_blob_pin_enforced": True,
        "parent_haar_and_rf30_operators_reused_without_formula_duplication": True,
        "repository_candidate_rf30_promotion_adapter_executable": True,
        "source_exact_auxiliary_mode_recovery_claimed": False,
        "agent2_curl_or_jacobian_reimplemented": False,
        "rf31_five_row_system_materialized": False,
        "rf34_rf39_correction_materialized": False,
        "rf44_rf49_remainder_recomputed": False,
        "cartesian_correction_velocity_materialized": False,
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
    "CurrentI4RF30PinnedProviderError",
    "CurrentI4RF30PinnedProviderReceipt",
    "materialize_current_i4_rf30_pinned_provider_evidence",
    "source_provenance",
    "truth_boundary",
]
