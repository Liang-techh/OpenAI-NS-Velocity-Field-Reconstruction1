"""Bind Agent-3 correction admission to one actual Agent-4 audit execution.

The inherited Agent-4 validator intentionally stayed byte-identical across the
axial-support rerun.  Its JSON receipt therefore still carries the historical
``parent_agent2_head`` of the protocol-origin provider, not the actual Agent-2
candidate checkout exercised by the newer workflow.  That field is useful for
protocol provenance, but it is not sufficient evidence for opening a correction
cycle on a later candidate.

This module adds a separate fail-closed execution-identity pin around the
existing Agent-3 scientific admission.  It binds the exact Agent-4 receipt bytes
and canonical payload to the actual audited Agent-2 head, Agent-4 head, workflow
run, and artifact identity.  It does not change any scientific guard.  The
currently pinned audit is a scientific REJECT, so correction ingestion remains
closed.  Any future Agent-4 PASS must receive a new Agent-3 identity pin instead
of reusing this receipt under a different candidate ancestry.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

from .kokuno_independent_oscillatory_admission import (
    EXPECTED_AGENT2_HEAD as PROTOCOL_ORIGIN_AGENT2_HEAD,
    EXPECTED_AUDIT_SCHEMA,
    evaluate_independent_oscillatory_admission,
)

TASK = "KOKUNO-A3-OSCILLATORY-AUDIT-IDENTITY-036"
SCHEMA = "kokuno-agent3-oscillatory-audit-identity-admission-v1"

# Actual execution independently recorded by Agent 4 PR #543.  These are not
# inferred from the legacy parent_agent2_head embedded in the byte-identical
# validator receipt.
ACTUAL_AUDITED_AGENT2_PR = 539
ACTUAL_AUDITED_AGENT2_HEAD = "985d1c3fc43081d2d87feb1f0c27643cc4ebe75b"
ACTUAL_AUDITOR_AGENT4_PR = 543
ACTUAL_AUDITOR_AGENT4_HEAD = "9498209376fd9dc5abf228248967958759e9f64a"
AUDIT_WORKFLOW_RUN_ID = 35416733758
AUDIT_ARTIFACT_ID = 10576031714
AUDIT_ARTIFACT_ZIP_DIGEST = (
    "sha256:8d7dbdc7654fc74820fb6916c5726984b42116055df146936b76607e6d1a0981"
)
AUDIT_RECEIPT_FILENAME = "kokuno_agent4_public_oscillatory_independent_audit.json"
AUDIT_RECEIPT_RAW_SHA256 = "c73f3aafe101923314cee0e9b6f95ef443409ffde2dba95ad0fa84a254683706"
AUDIT_RECEIPT_CANONICAL_SHA256 = "388582082b93c9476b8490fb8f4e483d7161ca22823f18f8fc9fed55c846aa52"


def _sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _canonical_payload_bytes(payload: Mapping[str, Any]) -> bytes:
    return json.dumps(
        payload,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
    ).encode("utf-8")


def _parse_and_bind_receipt(raw_receipt: bytes) -> Mapping[str, Any]:
    if not isinstance(raw_receipt, bytes):
        raise TypeError("raw_receipt must be bytes")

    raw_sha = _sha256_bytes(raw_receipt)
    if raw_sha != AUDIT_RECEIPT_RAW_SHA256:
        raise ValueError("audit receipt bytes do not match the pinned Agent-4 artifact member")

    try:
        payload = json.loads(raw_receipt.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ValueError("audit receipt is not valid UTF-8 JSON") from exc
    if not isinstance(payload, Mapping):
        raise ValueError("audit receipt JSON must be an object")

    canonical_sha = _sha256_bytes(_canonical_payload_bytes(payload))
    if canonical_sha != AUDIT_RECEIPT_CANONICAL_SHA256:
        raise ValueError("audit receipt canonical payload does not match the pinned execution")
    if payload.get("schema") != EXPECTED_AUDIT_SCHEMA:
        raise ValueError("pinned audit schema changed")
    if payload.get("parent_agent2_head") != PROTOCOL_ORIGIN_AGENT2_HEAD:
        raise ValueError("legacy protocol-origin Agent-2 head changed")
    return payload


def evaluate_pinned_oscillatory_audit_admission(raw_receipt: bytes) -> dict[str, Any]:
    """Recompute the v1 science gate and bind it to the exact #543 execution.

    The exact receipt hash is part of the admission contract.  Consequently a
    future scientific PASS cannot be substituted into this function and inherit
    #543's workflow/artifact identity.  Agent 3 must explicitly pin the future
    execution in a later increment.
    """

    audit = _parse_and_bind_receipt(raw_receipt)
    scientific = evaluate_independent_oscillatory_admission(audit)

    if scientific.get("consumed_audit_schema") != EXPECTED_AUDIT_SCHEMA:
        raise ValueError("parent admission consumed an unexpected audit schema")
    if scientific.get("consumed_agent2_head") != PROTOCOL_ORIGIN_AGENT2_HEAD:
        raise ValueError("parent admission lost the protocol-origin identity")

    # The exact trusted #543 receipt is a scientific REJECT.  If inherited code
    # ever interprets these exact bytes as PASS, fail closed rather than opening
    # correction machinery through an accidental threshold/semantic drift.
    if scientific.get("independent_public_oscillatory_preflight_passed") is not False:
        raise ValueError("pinned Agent-4 #543 receipt unexpectedly changed scientific verdict")
    if scientific.get("correction_ingest_allowed") is not False:
        raise ValueError("pinned rejected execution cannot open correction ingestion")

    return {
        "task": TASK,
        "schema": SCHEMA,
        "execution_identity_bound": True,
        "protocol_origin_agent2_head": PROTOCOL_ORIGIN_AGENT2_HEAD,
        "protocol_origin_semantics": "validator_protocol_origin_metadata_only",
        "actual_audited_agent2_pr": ACTUAL_AUDITED_AGENT2_PR,
        "actual_audited_agent2_head": ACTUAL_AUDITED_AGENT2_HEAD,
        "actual_auditor_agent4_pr": ACTUAL_AUDITOR_AGENT4_PR,
        "actual_auditor_agent4_head": ACTUAL_AUDITOR_AGENT4_HEAD,
        "audit_workflow_run_id": AUDIT_WORKFLOW_RUN_ID,
        "audit_artifact_id": AUDIT_ARTIFACT_ID,
        "audit_artifact_zip_digest": AUDIT_ARTIFACT_ZIP_DIGEST,
        "audit_receipt_filename": AUDIT_RECEIPT_FILENAME,
        "audit_receipt_raw_sha256": AUDIT_RECEIPT_RAW_SHA256,
        "audit_receipt_canonical_sha256": AUDIT_RECEIPT_CANONICAL_SHA256,
        "scientific_failed_guards": tuple(scientific["failed_guards"]),
        "independent_public_oscillatory_preflight_passed": False,
        "correction_ingest_allowed": False,
        "same_cycle_requested_stress_materialization_allowed": False,
        "candidate_finite_head_mean_debt_materialization_allowed": False,
        "requested_stress_actual_state_values_materialized": False,
        "finite_head_mean_debt_materialized": False,
        "real_candidate_defect_consumed": False,
        "signed_mean_inverse_input_ready": False,
        "public_velocity_correction_materialized": False,
        "finite_correction_cycle_rerun_allowed": False,
        "finite_correction_cycle_run": False,
        "heldout_ns_residual_assessed": False,
        "residual_reduction_claimed": False,
        "pde_validated": False,
        "surrogate_defect_used": False,
        "future_pass_requires_new_agent3_identity_pin": True,
        "scope": (
            "execution-provenance hardening only; exact Agent-4 #543 receipt remains a "
            "scientific reject and no correction state is materialized"
        ),
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--audit", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)

    report = evaluate_pinned_oscillatory_audit_admission(args.audit.read_bytes())
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(report, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
