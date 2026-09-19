"""Agent-3 admission for the independently accepted public-z oscillatory field.

Agent 4 PR #563 reran the unchanged public black-box protocol on Agent 2 PR #561
and independently accepted the oscillatory component.  This module pins that
exact execution identity, rechecks the frozen scientific guards from the exact
receipt bytes, and opens only the next Agent-3 ingest stage.

A PASS here authorizes materializing same-cycle theta/axial ``requestedStress``
and candidate-specific finite-head mean debt.  It does not materialize either
quantity, does not construct a correction velocity, and does not run or claim a
finite correction cycle or full Navier--Stokes validation.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

from .kokuno_carrier_resolved_audit_admission import (
    FD4_STEPS,
    FROZEN_GUARDS,
    PHASE_RESOLUTIONS,
    SEED,
    _canonical_payload_bytes,
    _recompute_scientific_verdict,
    _require_frozen_protocol,
)

TASK = "KOKUNO-A3-PUBLIC-Z-PASS-IDENTITY-038"
SCHEMA = "kokuno-agent3-public-z-pass-admission-v1"
EXPECTED_AUDIT_SCHEMA = "kokuno-agent4-public-z-pullback-independent-audit-v1"

ACTUAL_AUDITED_AGENT2_PR = 561
ACTUAL_AUDITED_AGENT2_HEAD = "732800ce4990464b49c8aa32d0dff4580f6684d4"
ACTUAL_AUDITOR_AGENT4_PR = 563
ACTUAL_AUDITOR_AGENT4_HEAD = "9b0f86012c53fa8e32a19f766dbc150931870425"
PROTOCOL_ORIGIN_AGENT4_HEAD = "9498209376fd9dc5abf228248967958759e9f64a"
PREVIOUS_AGENT4_EXECUTION_HEAD = "2a18a87db0ef8c6682f10fc5a6355c90b4168c7c"
AUDIT_WORKFLOW_RUN_ID = 35422203621
AUDIT_ARTIFACT_ID = 10577644750
AUDIT_ARTIFACT_ZIP_DIGEST = (
    "sha256:7205e041c31d2a68b345235bbd7b0f3316a50ec6c41fa3012c107789c6a07029"
)
AUDIT_RECEIPT_FILENAME = "agent4_563_public_z_pullback_independent_audit.json"
AUDIT_RECEIPT_RAW_SHA256 = "5bc3512cc6b185e510680fd97873352dc0170bd917cdf86100483b608372a6a8"
AUDIT_RECEIPT_CANONICAL_SHA256 = "93ceb24263a63a0c310aa6f19d28565acaf7c8e42baf577a6f10a3d883f34ded"


def _sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _parse_and_bind_receipt(raw_receipt: bytes) -> Mapping[str, Any]:
    if not isinstance(raw_receipt, bytes):
        raise TypeError("raw_receipt must be bytes")
    if _sha256_bytes(raw_receipt) != AUDIT_RECEIPT_RAW_SHA256:
        raise ValueError("audit receipt bytes do not match pinned Agent-4 #563 artifact member")
    try:
        payload = json.loads(raw_receipt.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ValueError("audit receipt is not valid UTF-8 JSON") from exc
    if not isinstance(payload, Mapping):
        raise ValueError("audit receipt JSON must be an object")
    if _sha256_bytes(_canonical_payload_bytes(payload)) != AUDIT_RECEIPT_CANONICAL_SHA256:
        raise ValueError("audit receipt canonical payload does not match pinned execution")
    if payload.get("schema") != EXPECTED_AUDIT_SCHEMA:
        raise ValueError("pinned audit schema changed")
    if payload.get("audited_agent2_head") != ACTUAL_AUDITED_AGENT2_HEAD:
        raise ValueError("pinned audit Agent-2 candidate identity changed")
    if payload.get("frozen_protocol_origin_agent4_head") != PROTOCOL_ORIGIN_AGENT4_HEAD:
        raise ValueError("pinned audit protocol-origin Agent-4 head changed")
    if payload.get("previous_agent4_execution_head") != PREVIOUS_AGENT4_EXECUTION_HEAD:
        raise ValueError("pinned audit previous Agent-4 execution identity changed")
    return payload


def evaluate_public_z_pass_admission(raw_receipt: bytes) -> dict[str, Any]:
    """Bind Agent-4 #563 and open only the next same-cycle defect ingest stage."""
    audit = _parse_and_bind_receipt(raw_receipt)
    _require_frozen_protocol(audit)
    scientific = _recompute_scientific_verdict(audit)

    declared_pairs = (
        ("local_divergence_passed", audit["divergence"]["local_divergence_passed"]),
        (
            "covariance_rank_preflight_passed",
            audit["phase_mean_covariance"]["covariance_rank_preflight_passed"],
        ),
        ("project_support_preflight_passed", audit["project_support_preflight_passed"]),
        ("public_oscillatory_preflight_passed", audit["public_oscillatory_preflight_passed"]),
    )
    for key, declared in declared_pairs:
        if bool(declared) is not bool(scientific[key]):
            raise ValueError(f"Agent-4 declared {key} disagrees with recomputed frozen guards")

    if tuple(float(row["step"]) for row in audit["divergence"]["rows"]) != FD4_STEPS:
        raise ValueError("pinned FD4 ladder changed")
    if tuple(
        int(row["phase_resolution"]) for row in audit["phase_mean_covariance"]["rows"]
    ) != PHASE_RESOLUTIONS:
        raise ValueError("pinned phase-resolution ladder changed")
    if audit.get("seed") != SEED:
        raise ValueError("pinned independent seed changed")
    if audit.get("guards_frozen_before_actions") != FROZEN_GUARDS:
        raise ValueError("pinned frozen guard dictionary changed")

    truth = audit["truth_boundary"]
    if truth.get("oscillatory_component_only") is not True:
        raise ValueError("pinned receipt must remain oscillatory-component only")
    if float(truth["formal_momentum_normalized_max_l2_gate"]) != 1.0e-3:
        raise ValueError("formal momentum gate changed")
    if float(truth["formal_divergence_max_l2_gate"]) != 1.0e-5:
        raise ValueError("formal divergence gate changed")
    for key in (
        "formal_full_domain_pde_gate_assessed",
        "global_leading_velocity_pressure_available",
        "heldout_ns_residual_assessed",
        "materialized_correction_available",
        "pde_validated",
    ):
        if truth.get(key) is not False:
            raise ValueError(f"component audit cannot promote {key}")

    if scientific["failed_guards"] != ():
        raise ValueError("pinned #563 receipt no longer passes every frozen local oscillatory guard")
    if scientific["local_divergence_passed"] is not True:
        raise ValueError("pinned #563 receipt lost independent divergence PASS")
    if scientific["covariance_rank_preflight_passed"] is not True:
        raise ValueError("pinned #563 receipt lost independent covariance PASS")
    if scientific["project_support_preflight_passed"] is not True:
        raise ValueError("pinned #563 receipt lost independent support PASS")
    if scientific["public_oscillatory_preflight_passed"] is not True:
        raise ValueError("pinned #563 receipt lost overall local oscillatory PASS")

    return {
        "task": TASK,
        "schema": SCHEMA,
        "execution_identity_bound": True,
        "actual_audited_agent2_pr": ACTUAL_AUDITED_AGENT2_PR,
        "actual_audited_agent2_head": ACTUAL_AUDITED_AGENT2_HEAD,
        "actual_auditor_agent4_pr": ACTUAL_AUDITOR_AGENT4_PR,
        "actual_auditor_agent4_head": ACTUAL_AUDITOR_AGENT4_HEAD,
        "protocol_origin_agent4_head": PROTOCOL_ORIGIN_AGENT4_HEAD,
        "previous_agent4_execution_head": PREVIOUS_AGENT4_EXECUTION_HEAD,
        "audit_workflow_run_id": AUDIT_WORKFLOW_RUN_ID,
        "audit_artifact_id": AUDIT_ARTIFACT_ID,
        "audit_artifact_zip_digest": AUDIT_ARTIFACT_ZIP_DIGEST,
        "audit_receipt_filename": AUDIT_RECEIPT_FILENAME,
        "audit_receipt_raw_sha256": AUDIT_RECEIPT_RAW_SHA256,
        "audit_receipt_canonical_sha256": AUDIT_RECEIPT_CANONICAL_SHA256,
        **scientific,
        "independent_public_z_oscillatory_preflight_passed": True,
        "correction_receipt_identity_pinned_to_current_pass": True,
        "correction_ingest_allowed": True,
        "same_cycle_requested_stress_materialization_allowed": True,
        "candidate_finite_head_mean_debt_materialization_allowed": True,
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
        "scope": (
            "exact Agent-2 #561 / Agent-4 #563 public-z execution is independently admitted; "
            "only same-cycle requestedStress and candidate finite-head mean-debt materialization "
            "are newly authorized, while correction construction and PDE claims remain closed"
        ),
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--audit", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    result = evaluate_public_z_pass_admission(args.audit.read_bytes())
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(result, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
