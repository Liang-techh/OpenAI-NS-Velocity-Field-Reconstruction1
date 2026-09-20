"""Fail-closed binding of an Agent-2 typed complete-curl backend into Agent 3.

Agent 3 owns actual-defect -> mean/radial -> signed-amplitude correction mechanics.
Agent 2 owns the complete-curl realization.  This module intentionally contains
no curl formula.  It accepts only the narrow payload emitted by the external
Agent-2 backend and turns that payload into Agent-3's existing
``CompleteCurlCorrectionAdapter``.

The pinned sibling interface is Agent-2 PR #734 exact head
``c535eeec3267869f630c1327506b3683d0417969``.  The binding itself does not
promote that backend to source-certified status; the certification bit is copied
verbatim and remains false for #734 pending independent audit.
"""
from __future__ import annotations

from typing import Any

from .kokuno_same_cycle_defect_contract import CycleIdentity
from .kokuno_typed_mean_velocity_correction_handoff import CompleteCurlCorrectionAdapter

TASK = "KOKUNO-A3-AGENT2-COMPLETE-CURL-BINDING-059"
SCHEMA = "kokuno-a3-agent2-complete-curl-binding-v1"
PARENT_AGENT3_PR = 726
PARENT_AGENT3_HEAD = "87c39f48113908c0755497de973c05ed66f1ba21"
AGENT2_BACKEND_PR = 734
AGENT2_BACKEND_HEAD = "c535eeec3267869f630c1327506b3683d0417969"

_EXPECTED_PAYLOAD_KEYS = frozenset(
    {
        "identity",
        "radii",
        "producer_kind",
        "provenance",
        "velocity_evaluator",
        "velocity_dt_evaluator",
        "source_agent2_complete_curl_certified",
    }
)


def bind_agent2_complete_curl_backend(
    backend: Any,
    identity: CycleIdentity,
) -> CompleteCurlCorrectionAdapter:
    """Bind one external A2 backend without importing/reimplementing its physics.

    The backend must expose ``agent3_adapter_kwargs(identity)`` and return exactly
    the seven fields in Agent-3's adapter contract.  Extra fields are rejected so
    a future backend cannot silently smuggle residuals, targets, gains, forcing,
    or scientific thresholds through this seam.
    """

    if not isinstance(identity, CycleIdentity):
        raise TypeError("identity must be a CycleIdentity")
    factory = getattr(backend, "agent3_adapter_kwargs", None)
    if not callable(factory):
        raise TypeError("Agent-2 backend must expose callable agent3_adapter_kwargs(identity)")
    payload = factory(identity)
    if not isinstance(payload, dict):
        raise TypeError("Agent-2 adapter payload must be a dict")
    keys = frozenset(payload)
    if keys != _EXPECTED_PAYLOAD_KEYS:
        missing = sorted(_EXPECTED_PAYLOAD_KEYS - keys)
        extra = sorted(keys - _EXPECTED_PAYLOAD_KEYS)
        raise ValueError(f"Agent-2 adapter payload schema mismatch; missing={missing}, extra={extra}")
    if payload["identity"] != identity:
        raise ValueError("Agent-2 adapter payload changed the requested cycle identity")
    return CompleteCurlCorrectionAdapter(**payload)


def truth_boundary() -> dict[str, object]:
    return {
        "schema": SCHEMA,
        "task": TASK,
        "parent_agent3_pr": PARENT_AGENT3_PR,
        "parent_agent3_head": PARENT_AGENT3_HEAD,
        "agent2_backend_pr": AGENT2_BACKEND_PR,
        "agent2_backend_head": AGENT2_BACKEND_HEAD,
        "agent2_complete_curl_reimplemented_by_agent3": False,
        "payload_schema_fail_closed": True,
        "caller_supplied_residual_allowed": False,
        "caller_supplied_defect_allowed": False,
        "caller_supplied_stress_allowed": False,
        "caller_supplied_target_allowed": False,
        "caller_supplied_gain_allowed": False,
        "caller_supplied_scientific_threshold_allowed": False,
        "source_agent2_complete_curl_certified_by_this_binding": False,
        "real_full_candidate_bound": False,
        "real_candidate_finite_correction_cycle_run": False,
        "heldout_normalized_ns_residual_assessed": False,
        "residual_reduction_claimed": False,
        "pde_validated": False,
    }
