"""Freshness gate for Agent-7 ST052 morphology/basis routing.

PR #1115 correctly requires a fresh axial measurement after A9's stable semantic
identity split, but it pins an older exact A9 #1112 head/run. A9 #1112 has since
advanced by two identity-governance commits without changing the frozen ST052
velocity. This module makes that handoff fail closed against stale upstream
pins while preserving the existing rule: stable identity success still requires
one fresh axial-aspect receipt before any existing control or new basis can be
selected.

This is governance/routing only. It does not modify velocity, basis dimension,
pressure, forcing, residuals, thresholds, or visual targets.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path
from typing import Any

TASK_ID = "CR003-ST052M-STABLE-IDENTITY-FRESHNESS-133"
PREREG_ISSUE = 1122
PARENT_PR = 1115
PARENT_HEAD = "78c269db3daa90a2fe561bbcf953f111a765d80b"

STABLE_PR = 1112
OLD_STABLE_HEAD = "e0bf8eabafa34553cbd6b96460bfda9c003abe50"
OLD_STABLE_RUN = 35681756248
LATEST_STABLE_HEAD = "97630f3af72d3b3690968132a6a16a3592d3a8ba"
LATEST_STABLE_RUN = 35685567737
LATEST_STABLE_TESTS_RUN = 35685567705
EXPECTED_AHEAD_BY = 2
EXPECTED_CHANGED_PATHS = {
    "src/openai_ns_reconstruction/st052_stable_semantic_identity.py",
    "tests/test_st052_stable_semantic_identity.py",
}

STABLE_SCHEMA = "st052-stable-semantic-identity-split/v1"
STABLE_TASK_ID = "CR-A9-105"
CANDIDATE_ID = "ST052-M-linear-temporal-child-v1"
CONSTRAINED_RUNTIME_REPOSITORY = "Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1"
CONSTRAINED_RUNTIME_IMPLEMENTATION_COMMIT = "6a293b3c870d70b8d8aece1b12d68685c9b8b010"
FAILED_AXIAL_PR = 1095
FAILED_AXIAL_RUN = 35674391631
FAILED_AXIAL_ASPECT = 1.0813103938217585

_SHA256 = re.compile(r"^[0-9a-f]{64}$")
_SHA1 = re.compile(r"^[0-9a-f]{40}$")

TRUTH = {
    "candidate_velocity_changed": False,
    "basis_dimension_changed": False,
    "coefficient_selected": False,
    "existing_control_selected": False,
    "targeted_basis_preflight_authorized": False,
    "candidate_mutation_authorized": False,
    "public_openai_numeric_target_used": False,
    "historical_1095_measurement_admitted": False,
    "held_out_pde_residual_evaluated": False,
    "direct_visualization_fingerprint_improvement": 0.0,
    "closer_visualization_delivery_established": False,
    "visualization_ready": False,
    "visual_correspondence_verified": False,
    "source_correspondence_verified": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


def _canonical_sha256(value: Any) -> str:
    data = json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=False, allow_nan=False
    ).encode("utf-8")
    return hashlib.sha256(data).hexdigest()


def _require_dict(value: Any, label: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise ValueError(f"{label} must be an object")
    return value


def _require_sha256(value: Any, label: str) -> str:
    if not isinstance(value, str) or _SHA256.fullmatch(value) is None:
        raise ValueError(f"{label} must be lowercase 64-hex SHA-256")
    return value


def _require_sha1(value: Any, label: str) -> str:
    if not isinstance(value, str) or _SHA1.fullmatch(value) is None:
        raise ValueError(f"{label} must be lowercase 40-hex SHA-1")
    return value


def validate_upstream_pr(pr: dict[str, Any]) -> None:
    if int(pr.get("number")) != STABLE_PR:
        raise ValueError("unexpected stable-identity PR number")
    head = _require_dict(pr.get("head"), "stable PR head")
    if head.get("sha") != LATEST_STABLE_HEAD:
        raise ValueError("stable-identity PR head rolled again; preregistration is stale")
    if pr.get("state") != "open":
        raise ValueError("stable-identity PR state changed")


def validate_upstream_compare(compare: dict[str, Any]) -> dict[str, Any]:
    base = _require_dict(compare.get("base_commit"), "base_commit")
    head = compare.get("head_commit")
    if head is None:
        head_sha = LATEST_STABLE_HEAD
    else:
        head_sha = _require_dict(head, "head_commit").get("sha")
    if base.get("sha") != OLD_STABLE_HEAD:
        raise ValueError("unexpected compare base")
    if head_sha != LATEST_STABLE_HEAD:
        raise ValueError("unexpected compare head")
    if compare.get("status") != "ahead":
        raise ValueError("latest stable-identity head is not a strict descendant")
    if int(compare.get("ahead_by", -1)) != EXPECTED_AHEAD_BY:
        raise ValueError("unexpected stable-identity upstream commit count")
    if int(compare.get("behind_by", -1)) != 0:
        raise ValueError("stable-identity upstream history diverged")
    files = compare.get("files")
    if not isinstance(files, list):
        raise ValueError("compare files must be a list")
    changed = {f.get("filename") for f in files if isinstance(f, dict)}
    if changed != EXPECTED_CHANGED_PATHS:
        raise ValueError("unexpected upstream changed path; candidate/velocity semantics require fresh audit")
    for f in files:
        if f.get("status") not in {"modified"}:
            raise ValueError("upstream freshness compare must only modify identity module/tests")
    return {
        "ahead_by": EXPECTED_AHEAD_BY,
        "behind_by": 0,
        "changed_paths": sorted(changed),
        "candidate_or_velocity_implementation_changed": False,
    }


def validate_latest_run(run: dict[str, Any]) -> tuple[str, Any]:
    if int(run.get("id")) != LATEST_STABLE_RUN:
        raise ValueError("unexpected latest stable-identity run id")
    if run.get("head_sha") != LATEST_STABLE_HEAD:
        raise ValueError("latest stable-identity run head mismatch")
    status = run.get("status")
    conclusion = run.get("conclusion")
    if status not in {"queued", "in_progress", "completed"}:
        raise ValueError("unexpected latest stable-identity run status")
    if status != "completed" and conclusion is not None:
        raise ValueError("non-completed run cannot carry a conclusion")
    return str(status), conclusion


def validate_stable_receipt(receipt: dict[str, Any]) -> dict[str, str]:
    if receipt.get("schema") != STABLE_SCHEMA or receipt.get("task_id") != STABLE_TASK_ID:
        raise ValueError("unexpected stable-identity receipt schema/task")
    if receipt.get("candidate_id") != CANDIDATE_ID:
        raise ValueError("stable-identity candidate id drift")
    stable = _require_dict(receipt.get("stable_identity"), "stable_identity")
    truth = _require_dict(receipt.get("truth_boundary"), "truth_boundary")
    candidate_id = _require_sha256(
        stable.get("candidate_semantic_identity_sha256"), "candidate semantic identity"
    )
    velocity_id = _require_sha256(
        stable.get("velocity_semantic_identity_sha256"), "velocity semantic identity"
    )
    velocity_payload = _require_dict(
        stable.get("velocity_semantic_payload"), "velocity_semantic_payload"
    )
    if velocity_payload.get("candidate_semantic_identity_sha256") != candidate_id:
        raise ValueError("stable velocity/candidate linkage mismatch")
    runtime = _require_dict(
        velocity_payload.get("constrained_callable_runtime"), "constrained_callable_runtime"
    )
    if runtime.get("repository") != CONSTRAINED_RUNTIME_REPOSITORY:
        raise ValueError("stable callable repository drift")
    implementation_commit = _require_sha1(runtime.get("commit"), "stable callable commit")
    if implementation_commit != CONSTRAINED_RUNTIME_IMPLEMENTATION_COMMIT:
        raise ValueError("stable callable implementation commit drift")
    if truth.get("candidate_changed") is not False:
        raise ValueError("stable identity receipt claims candidate change")
    if truth.get("velocity_coefficients_changed") is not False:
        raise ValueError("stable identity receipt claims velocity coefficient change")
    if truth.get("constrained_callable_implementation_identity_bound") is not True:
        raise ValueError("latest stable callable identity is not implementation-bound")
    for key in (
        "velocity_export_ready",
        "visualization_ready",
        "visual_correspondence_verified",
        "source_correspondence_verified",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    ):
        if truth.get(key) is not False:
            raise ValueError(f"stable identity truth promotion forbidden: {key}")
    return {
        "candidate_semantic_identity_sha256": candidate_id,
        "velocity_semantic_identity_sha256": velocity_id,
        "constrained_callable_implementation_commit": implementation_commit,
    }


def build_receipt(
    *,
    upstream_pr: dict[str, Any],
    upstream_compare: dict[str, Any],
    latest_run: dict[str, Any],
    stable_receipt: dict[str, Any] | None = None,
) -> dict[str, Any]:
    validate_upstream_pr(upstream_pr)
    compare_summary = validate_upstream_compare(upstream_compare)
    status, conclusion = validate_latest_run(latest_run)

    decision: dict[str, Any]
    stable_ids: dict[str, str] | None = None
    if status == "completed" and conclusion == "success":
        if stable_receipt is None:
            raise ValueError("successful latest stable-identity run requires its receipt")
        stable_ids = validate_stable_receipt(stable_receipt)
        decision = {
            "status": "latest_stable_identity_ready_fresh_axial_rerun_required",
            "fresh_axial_rerun_required": True,
            "fresh_axial_measurement_admitted": False,
            "handoff_to_existing_axial_router_authorized": False,
            "existing_control_selection_authorized": False,
            "targeted_basis_preflight_authorized": False,
            "candidate_mutation_authorized": False,
        }
    else:
        if stable_receipt is not None:
            raise ValueError("unresolved/failed latest run cannot admit a stable receipt")
        decision = {
            "status": "blocked_latest_stable_identity_ci",
            "fresh_axial_rerun_required": True,
            "fresh_axial_measurement_admitted": False,
            "handoff_to_existing_axial_router_authorized": False,
            "existing_control_selection_authorized": False,
            "targeted_basis_preflight_authorized": False,
            "candidate_mutation_authorized": False,
        }

    receipt: dict[str, Any] = {
        "schema": "agent7-st052m-stable-identity-freshness/v1",
        "task_id": TASK_ID,
        "prereg_issue": PREREG_ISSUE,
        "parent": {"pr": PARENT_PR, "head": PARENT_HEAD},
        "upstream_rollover": {
            "pr": STABLE_PR,
            "old_head": OLD_STABLE_HEAD,
            "old_run": OLD_STABLE_RUN,
            "latest_head": LATEST_STABLE_HEAD,
            "latest_run": LATEST_STABLE_RUN,
            "latest_tests_run": LATEST_STABLE_TESTS_RUN,
            "old_pin_superseded_for_future_agent7_routing": True,
            "compare": compare_summary,
        },
        "latest_run_state": {"status": status, "conclusion": conclusion},
        "latest_stable_identity": stable_ids,
        "historical_failed_axial_evidence": {
            "pr": FAILED_AXIAL_PR,
            "run": FAILED_AXIAL_RUN,
            "descriptive_aspect": FAILED_AXIAL_ASPECT,
            "admitted_for_basis_routing": False,
            "retroactive_relabelling_forbidden": True,
        },
        "decision": decision,
        "next_legal_step": (
            "wait_for_latest_stable_identity_ci_then_rerun_axial_aspect_on_exact_latest_velocity_semantic_identity"
        ),
        "truth": dict(TRUTH),
    }
    unsigned = dict(receipt)
    receipt["receipt_sha256"] = _canonical_sha256(unsigned)
    validate_receipt(receipt)
    return receipt


def validate_receipt(receipt: dict[str, Any]) -> None:
    if receipt.get("schema") != "agent7-st052m-stable-identity-freshness/v1":
        raise ValueError("unexpected Agent-7 freshness schema")
    if receipt.get("task_id") != TASK_ID or receipt.get("prereg_issue") != PREREG_ISSUE:
        raise ValueError("task/prereg drift")
    if receipt.get("truth") != TRUTH:
        raise ValueError("truth boundary drift")
    rollover = _require_dict(receipt.get("upstream_rollover"), "upstream_rollover")
    if rollover.get("latest_head") != LATEST_STABLE_HEAD or rollover.get("latest_run") != LATEST_STABLE_RUN:
        raise ValueError("latest upstream pin drift")
    if rollover.get("old_pin_superseded_for_future_agent7_routing") is not True:
        raise ValueError("stale upstream pin was not superseded")
    historical = _require_dict(
        receipt.get("historical_failed_axial_evidence"), "historical_failed_axial_evidence"
    )
    if historical.get("admitted_for_basis_routing") is not False:
        raise ValueError("historical failed axial evidence cannot be admitted")
    decision = _require_dict(receipt.get("decision"), "decision")
    if decision.get("status") not in {
        "blocked_latest_stable_identity_ci",
        "latest_stable_identity_ready_fresh_axial_rerun_required",
    }:
        raise ValueError("unexpected freshness decision")
    if decision.get("fresh_axial_rerun_required") is not True:
        raise ValueError("fresh axial rerun must remain required")
    for key in (
        "fresh_axial_measurement_admitted",
        "handoff_to_existing_axial_router_authorized",
        "existing_control_selection_authorized",
        "targeted_basis_preflight_authorized",
        "candidate_mutation_authorized",
    ):
        if decision.get(key) is not False:
            raise ValueError(f"premature routing promotion: {key}")
    recorded = _require_sha256(receipt.get("receipt_sha256"), "receipt_sha256")
    unsigned = dict(receipt)
    unsigned.pop("receipt_sha256")
    if recorded != _canonical_sha256(unsigned):
        raise ValueError("receipt checksum mismatch")


def _load(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"{path} root must be an object")
    return value


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--upstream-pr", type=Path, required=True)
    parser.add_argument("--upstream-compare", type=Path, required=True)
    parser.add_argument("--latest-run", type=Path, required=True)
    parser.add_argument("--stable-receipt", type=Path)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    receipt = build_receipt(
        upstream_pr=_load(args.upstream_pr),
        upstream_compare=_load(args.upstream_compare),
        latest_run=_load(args.latest_run),
        stable_receipt=_load(args.stable_receipt) if args.stable_receipt else None,
    )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"status": receipt["decision"]["status"], "out": str(args.out)}))


if __name__ == "__main__":
    main()
