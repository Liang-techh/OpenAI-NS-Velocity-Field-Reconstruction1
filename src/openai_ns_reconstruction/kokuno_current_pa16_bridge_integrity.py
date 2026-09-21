"""Checksum guard for the current-Xi bridge consumed by the PA.16 repair lane.

This is a CR002 representation/provenance guard.  It does not change the
five-moment data, PA.16 solver, velocity, pressure, forcing, or any scientific
threshold.  The guard only replays the exact checksum payload emitted by the
Agent-3 #949 current-Xi bridge before allowing the Agent-3 #957 PA.16 repair
binder to consume that receipt.
"""
from __future__ import annotations

from dataclasses import asdict
import hashlib
import json
from typing import Any

from .kokuno_current_pa16_repair_receipt import (
    CurrentPA16RepairBackend,
    CurrentPA16RepairReceipt,
    CurrentPA16RepairReceiptError,
    materialize_current_pa16_repair_receipt as _materialize_parent,
)
from .kokuno_current_xi_moment_discrepancy_bridge import (
    CurrentXiMomentBridgeReceipt,
    SCHEMA as XI_BRIDGE_SCHEMA,
)


TASK = "CR002-CURRENT-PA16-UPSTREAM-BRIDGE-INTEGRITY-093"
SCHEMA = "cr002-current-pa16-upstream-bridge-integrity-v1"
PARENT_PR = 957
PARENT_HEAD = "bb8adedbb4fd32519a00b78765c0fa62f16bfe52"
SOURCE_BRIDGE_PR = 949
SOURCE_BRIDGE_HEAD = "87c1ef765d6bf3613c7c6c301669ffccf5260117"


def _canonical_json(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _sha256_json(value: Any) -> str:
    return hashlib.sha256(_canonical_json(value).encode("utf-8")).hexdigest()


def replay_current_xi_handoff_sha256(receipt: CurrentXiMomentBridgeReceipt) -> str:
    """Replay the exact #949 handoff checksum payload from typed receipt fields."""
    if not isinstance(receipt, CurrentXiMomentBridgeReceipt):
        raise CurrentPA16RepairReceiptError(
            "current-Xi bridge integrity guard requires CurrentXiMomentBridgeReceipt"
        )
    payload = {
        "schema": XI_BRIDGE_SCHEMA,
        "upstream": asdict(receipt.upstream),
        "upstream_semantic_sha256": receipt.upstream_semantic_sha256,
        "upstream_report_sha256": receipt.upstream_report_sha256,
        "outer_schedule_sha256": receipt.outer_schedule_sha256,
        "samples": [asdict(sample) for sample in receipt.samples],
    }
    return _sha256_json(payload)


def validate_current_xi_bridge_handoff(
    receipt: CurrentXiMomentBridgeReceipt,
) -> CurrentXiMomentBridgeReceipt:
    """Fail closed if the stored #949 handoff SHA does not replay from its payload."""
    expected = replay_current_xi_handoff_sha256(receipt)
    if receipt.handoff_sha256 != expected:
        raise CurrentPA16RepairReceiptError(
            "current-Xi bridge handoff SHA does not replay from checksum-bound payload"
        )

    if receipt.samples:
        max_raw = max(sample.raw_l2_norm for sample in receipt.samples)
        max_transformed = max(
            sample.pa16_transformed_l2_norm for sample in receipt.samples
        )
    else:
        raise CurrentPA16RepairReceiptError(
            "current-Xi bridge integrity guard requires at least one sample"
        )
    if receipt.max_raw_l2_norm != max_raw:
        raise CurrentPA16RepairReceiptError(
            "current-Xi bridge max_raw_l2_norm does not replay from samples"
        )
    if receipt.max_pa16_transformed_l2_norm != max_transformed:
        raise CurrentPA16RepairReceiptError(
            "current-Xi bridge max_pa16_transformed_l2_norm does not replay from samples"
        )
    if receipt.any_correction_needed != any(
        sample.correction_needed for sample in receipt.samples
    ):
        raise CurrentPA16RepairReceiptError(
            "current-Xi bridge correction-needed summary does not replay from samples"
        )
    return receipt


class _PinnedBridgeBackend:
    """Delegate all #957 backend methods except the already-validated bridge."""

    def __init__(self, backend: CurrentPA16RepairBackend, bridge: CurrentXiMomentBridgeReceipt):
        self._backend = backend
        self._bridge = bridge

    def current_xi_bridge_receipt(self) -> CurrentXiMomentBridgeReceipt:
        return self._bridge

    def upstream_identity(self):
        return self._backend.upstream_identity()

    def tsh_report(self):
        return self._backend.tsh_report()

    def pa16_solve_at_eta(self, eta: float):
        return self._backend.pa16_solve_at_eta(eta)


def materialize_checksum_verified_current_pa16_repair_receipt(
    backend: CurrentPA16RepairBackend,
) -> CurrentPA16RepairReceipt:
    """Validate and pin one #949 bridge receipt before delegating to #957."""
    method = getattr(backend, "current_xi_bridge_receipt", None)
    if not callable(method):
        raise CurrentPA16RepairReceiptError(
            "backend is missing typed PA.16 method: current_xi_bridge_receipt"
        )
    bridge = validate_current_xi_bridge_handoff(method())
    return _materialize_parent(_PinnedBridgeBackend(backend, bridge))


def truth_boundary() -> dict[str, bool]:
    return {
        "source_bridge_checksum_payload_replayed": True,
        "source_bridge_checksum_failure_fails_closed": True,
        "validated_bridge_pinned_across_parent_delegation": True,
        "parent_957_entrypoint_itself_modified": False,
        "source_moment_receipt_promoted_to_complete_ns_defect": False,
        "authorized_for_gain_gated_ns_stage": False,
        "velocity_modified": False,
        "pressure_modified": False,
        "forcing_modified": False,
        "scientific_threshold_modified": False,
        "pde_validated": False,
        "visual_correspondence_verified": False,
        "paper_exact": False,
        "openai_field_identified": False,
    }
