"""Require an exact runtime rebind before current-I4 correction evidence is promoted.

Kokuno Agent 3 previously pinned the A2 auxiliary-T2 provider and routed its
RF30 defect through RF31/RF34--RF39.  That path still accepted an already-built
``ExactCurrentI4NonlinearBackend`` dataclass.  A dataclass instance can be
constructed directly, so a mechanics fixture can carry plausible hashes while
using a synthetic leading field or differential function.  Such a fixture is
useful for interface tests but is not evidence that the real A2 #1080 / #960
candidate was recomputed.

This successor is therefore fail-closed: it first reconstructs the backend via
``ExactCurrentI4NonlinearBackend.bind`` from the supplied composite object and
differential function.  ``bind`` authenticates the concrete module/class/function
identities and Git blob SHAs.  The correction is materialized only on that
rebound backend, never on the caller's unchecked dataclass instance.

No RF30/RF31/RF34--RF39 formula is changed here.  No defect, covariance,
coefficient, correction array, forcing, pressure, held-out sample or scientific
threshold is accepted from the caller.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import inspect
import json
import math
from typing import Any

from .kokuno_current_i4_nonlinear_mean_attribution import ExactCurrentI4NonlinearBackend
from .kokuno_current_i4_rf34_rf39_pinned_provider import (
    CurrentI4PinnedRF34RF39Receipt,
    materialize_current_i4_rf34_rf39_pinned_provider_correction,
)
from .kokuno_rf30_auxiliary_t2_haar_bridge import SourceAuxiliaryT2WaveProvider

TASK = "KOKUNO-A3-CURRENT-I4-EXACT-BACKEND-FIREWALL-148"
SCHEMA = "kokuno-a3-current-i4-exact-backend-firewall-v1"

PARENT_PR = 1210
PARENT_HEAD = "e75940c33127b0725ddc703eaae97590158f2a56"
PARENT_BLOB = "b99b3bee2b2dad6003c11acae7946c9be53347fb"

FINAL_NORMALIZED_MOMENTUM_GATE = 1.0e-3
FINAL_NORMALIZED_DIVERGENCE_GATE = 1.0e-5


class CurrentI4ExactBackendFirewallError(RuntimeError):
    """Raised when an unchecked/synthetic backend reaches the correction seam."""


def _canonical_sha(payload: dict[str, object]) -> str:
    raw = json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)
    return hashlib.sha256(raw.encode("utf-8")).hexdigest()


def _backend_identity(backend: ExactCurrentI4NonlinearBackend) -> dict[str, object]:
    return {
        "composite_semantic_sha256": backend.composite_semantic_sha256,
        "oscillatory_runtime_sha256": backend.oscillatory_runtime_sha256,
        "differential_semantic_sha256": backend.differential_semantic_sha256,
        "composite_source_blob": backend.composite_source_blob,
        "differential_source_blob": backend.differential_source_blob,
        "typed_receipt": backend.to_receipt(),
    }


def _require_exact_rebind(
    backend: ExactCurrentI4NonlinearBackend,
) -> ExactCurrentI4NonlinearBackend:
    if not isinstance(backend, ExactCurrentI4NonlinearBackend):
        raise TypeError("backend must be ExactCurrentI4NonlinearBackend")

    try:
        rebound = ExactCurrentI4NonlinearBackend.bind(
            backend.composite_field,
            backend.differential_function,
        )
    except Exception as exc:  # exact bind failures are scientific fail-closed events
        raise CurrentI4ExactBackendFirewallError(
            "exact A2 #1080/#960 backend rebind failed; mechanics fixtures are not correction evidence"
        ) from exc

    if _backend_identity(rebound) != _backend_identity(backend):
        raise CurrentI4ExactBackendFirewallError(
            "caller backend identity disagrees with exact rebound backend"
        )
    return rebound


@dataclass(frozen=True)
class CurrentI4ExactBackendCorrectionReceipt:
    candidate_id: str
    candidate_sha256: str
    source_chart_id: str
    source_chart_sha256: str
    exact_backend_receipt_sha256: str
    parent_correction_receipt_sha256: str
    correction_payload_sha256: str
    coefficients_u0_u1_u2_s0_s1: tuple[float, float, float, float, float]
    correction_chart_l2: float
    correction_chart_max: float
    five_row_closure_max_relative: float
    coefficient_direct_solve_max_relative: float
    exact_backend_rebind_executed: bool
    concrete_a2_composite_identity_authenticated: bool
    concrete_a2_differential_identity_authenticated: bool
    caller_dataclass_identity_replayed_exactly: bool
    pinned_auxiliary_t2_provider_consumed: bool
    repository_candidate_rf30_defect_evidence: bool
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
            "truth_boundary": truth_boundary(),
        }


def materialize_current_i4_rf34_rf39_exact_backend_correction(
    backend: ExactCurrentI4NonlinearBackend,
    provider: SourceAuxiliaryT2WaveProvider,
    radius: Any,
    z: Any,
    t: Any,
) -> CurrentI4ExactBackendCorrectionReceipt:
    """Materialize #1210 only after exact concrete backend authentication."""
    rebound = _require_exact_rebind(backend)
    exact_backend_payload = _backend_identity(rebound)
    exact_backend_sha = _canonical_sha(exact_backend_payload)

    parent: CurrentI4PinnedRF34RF39Receipt = (
        materialize_current_i4_rf34_rf39_pinned_provider_correction(
            rebound,
            provider,
            radius,
            z,
            t,
        )
    )

    required_true = (
        parent.exact_provider_blob_pinned,
        parent.same_typed_rf30_receipt_preserved,
        parent.same_candidate_identity_preserved,
        parent.same_source_chart_identity_preserved,
        parent.repository_candidate_rf31_system_evidence,
        parent.repository_candidate_rf34_rf39_correction_evidence,
    )
    if not all(required_true):
        raise CurrentI4ExactBackendFirewallError(
            "parent correction did not preserve the pinned repository-candidate chain"
        )

    required_false = (
        parent.correction_applied_to_candidate,
        parent.rf44_rf49_nonlinear_remainder_recomputed,
        parent.cartesian_correction_velocity_materialized,
        parent.finite_correction_cycle_run,
        parent.heldout_ns_residual_assessed,
        parent.pde_validated,
    )
    if any(required_false):
        raise CurrentI4ExactBackendFirewallError(
            "parent unexpectedly promoted a downstream correction/PDE state"
        )

    diagnostics = (
        *parent.coefficients_u0_u1_u2_s0_s1,
        parent.correction_chart_l2,
        parent.correction_chart_max,
        parent.five_row_closure_max_relative,
        parent.coefficient_direct_solve_max_relative,
    )
    if not all(math.isfinite(float(v)) for v in diagnostics):
        raise CurrentI4ExactBackendFirewallError("parent correction diagnostics are non-finite")

    payload: dict[str, object] = {
        "schema": SCHEMA,
        "task": TASK,
        "candidate_id": parent.candidate_id,
        "candidate_sha256": parent.candidate_sha256,
        "source_chart_id": parent.source_chart_id,
        "source_chart_sha256": parent.source_chart_sha256,
        "exact_backend_receipt_sha256": exact_backend_sha,
        "parent_correction_receipt_sha256": parent.receipt_sha256,
        "correction_payload_sha256": parent.correction_payload_sha256,
        "exact_backend_rebind_executed": True,
        "repository_candidate_rf34_rf39_correction_evidence": True,
    }
    receipt_sha = _canonical_sha(payload)

    return CurrentI4ExactBackendCorrectionReceipt(
        candidate_id=parent.candidate_id,
        candidate_sha256=parent.candidate_sha256,
        source_chart_id=parent.source_chart_id,
        source_chart_sha256=parent.source_chart_sha256,
        exact_backend_receipt_sha256=exact_backend_sha,
        parent_correction_receipt_sha256=parent.receipt_sha256,
        correction_payload_sha256=parent.correction_payload_sha256,
        coefficients_u0_u1_u2_s0_s1=tuple(
            float(v) for v in parent.coefficients_u0_u1_u2_s0_s1
        ),
        correction_chart_l2=float(parent.correction_chart_l2),
        correction_chart_max=float(parent.correction_chart_max),
        five_row_closure_max_relative=float(parent.five_row_closure_max_relative),
        coefficient_direct_solve_max_relative=float(
            parent.coefficient_direct_solve_max_relative
        ),
        exact_backend_rebind_executed=True,
        concrete_a2_composite_identity_authenticated=True,
        concrete_a2_differential_identity_authenticated=True,
        caller_dataclass_identity_replayed_exactly=True,
        pinned_auxiliary_t2_provider_consumed=True,
        repository_candidate_rf30_defect_evidence=True,
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
        materialize_current_i4_rf34_rf39_exact_backend_correction
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
        "parent_pr": PARENT_PR,
        "parent_head": PARENT_HEAD,
        "parent_blob": PARENT_BLOB,
        "public_parameters": tuple(signature.parameters),
        "forbidden_scientific_controls_exposed": bool(
            forbidden.intersection(signature.parameters)
        ),
        "exact_backend_rebind_required": True,
        "exact_backend_bind_authenticates_concrete_module_class_function_and_blobs": True,
        "manual_ExactCurrentI4NonlinearBackend_dataclass_construction_is_sufficient": False,
        "zero_or_synthetic_backend_fixture_is_scientific_correction_evidence": False,
        "parent_1210_promotion_without_exact_rebind_is_sufficient": False,
        "parent_rf30_rf31_rf34_rf39_formulas_changed": False,
        "agent2_curl_or_jacobian_reimplemented": False,
        "repository_candidate_rf34_rf39_correction_materializable_only_after_exact_rebind": True,
        "correction_applied_to_candidate": False,
        "rf44_rf49_nonlinear_remainder_recomputed": False,
        "cartesian_correction_velocity_materialized": False,
        "finite_correction_cycle_run": False,
        "heldout_ns_residual_assessed": False,
        "residual_reduction_claimed": False,
        "residual_defined_free_forcing_allowed": False,
        "final_normalized_momentum_gate": FINAL_NORMALIZED_MOMENTUM_GATE,
        "final_normalized_divergence_gate": FINAL_NORMALIZED_DIVERGENCE_GATE,
        "pde_validated": False,
    }
