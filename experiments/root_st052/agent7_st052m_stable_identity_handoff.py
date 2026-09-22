"""Gate ST052 axial basis routing on the repaired stable semantic identity.

Preregistered in Issue #1114 and stacked exactly on Agent-7 PR #1106.
A9 PR #1112 introduces a stable semantic/callable identity that deliberately
excludes regenerated materialization hashes.  That repair is necessary, but it
does not retroactively admit the failed #1095 axial-aspect measurement.

This module therefore permits only the following progression:
  upstream stable-identity CI -> fresh stable-ID-bound axial measurement ->
  existing Agent-7 axial router/Jacobian.

No candidate, coefficient, basis, pressure, forcing, or scientific threshold is
changed here.
"""
from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path
from typing import Any

TASK_ID = "CR003-ST052M-STABLE-IDENTITY-HANDOFF-132"
PREREG_ISSUE = 1114
SOURCE_PARENT_PR = 1106
SOURCE_PARENT_HEAD = "f90de7bb4126934823370ee64d1f1c9625d6ece2"

STABLE_PR = 1112
STABLE_HEAD = "e0bf8eabafa34553cbd6b96460bfda9c003abe50"
STABLE_RUN = 35681756248
STABLE_SCHEMA = "st052-stable-semantic-identity-split/v1"
STABLE_TASK_ID = "CR-A9-105"
CANDIDATE_ID = "ST052-M-linear-temporal-child-v1"

FAILED_AXIAL_PR = 1095
FAILED_AXIAL_RUN = 35674391631
FAILED_AXIAL_ASPECT = 1.0813103938217585
LEGACY_ADMITTED_VELOCITY_SHA256 = "b0654bd8955644845fcb3789628299318b43bdb72bb467c284822235ffbeb902"
LEGACY_REPLAY_VELOCITY_SHA256 = "de28f374386e78679f79006e2dbdbdf238179267d13f8470796c44d7f36ebb42"
_SHA256 = re.compile(r"^[0-9a-f]{64}$")

TRUTH = {
    "candidate_velocity_changed": False,
    "basis_dimension_changed": False,
    "coefficient_selected": False,
    "existing_control_selected": False,
    "public_openai_numeric_target_used": False,
    "failed_1095_measurement_retroactively_admitted": False,
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


def _load(path: str | Path) -> dict[str, Any]:
    value = json.loads(Path(path).read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError("JSON root must be an object")
    return value


def _sha256(value: Any, label: str) -> str:
    if not isinstance(value, str) or _SHA256.fullmatch(value) is None:
        raise ValueError(f"{label} must be lowercase 64-hex SHA-256")
    return value


def validate_stable_run(run: dict[str, Any]) -> tuple[str, Any]:
    if int(run.get("id")) != STABLE_RUN:
        raise ValueError("unexpected stable-identity run id")
    if run.get("head_sha") != STABLE_HEAD:
        raise ValueError("unexpected stable-identity run head")
    status = run.get("status")
    conclusion = run.get("conclusion")
    if status not in {"queued", "in_progress", "completed"}:
        raise ValueError("unexpected stable-identity run status")
    if status != "completed" and conclusion is not None:
        raise ValueError("non-completed run cannot carry a conclusion")
    return str(status), conclusion


def validate_stable_receipt(receipt: dict[str, Any]) -> dict[str, str]:
    if receipt.get("schema") != STABLE_SCHEMA or receipt.get("task_id") != STABLE_TASK_ID:
        raise ValueError("unexpected stable-identity receipt schema/task")
    if receipt.get("candidate_id") != CANDIDATE_ID:
        raise ValueError("stable-identity candidate id drift")
    stable = receipt.get("stable_identity")
    truth = receipt.get("truth_boundary")
    evidence = receipt.get("materialization_evidence")
    migration = receipt.get("migration_boundary")
    if not all(isinstance(v, dict) for v in (stable, truth, evidence, migration)):
        raise ValueError("malformed stable-identity receipt")

    candidate_id = _sha256(
        stable.get("candidate_semantic_identity_sha256"),
        "candidate_semantic_identity_sha256",
    )
    velocity_id = _sha256(
        stable.get("velocity_semantic_identity_sha256"),
        "velocity_semantic_identity_sha256",
    )
    if velocity_id in {LEGACY_ADMITTED_VELOCITY_SHA256, LEGACY_REPLAY_VELOCITY_SHA256}:
        raise ValueError("legacy velocity hash was relabelled as stable semantic identity")
    if truth.get("stable_semantic_identity_ready") is not True:
        raise ValueError("stable semantic identity not marked ready")
    for key in (
        "velocity_export_ready",
        "visualization_ready",
        "visual_correspondence_verified",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    ):
        if truth.get(key) is not False:
            raise ValueError(f"stable receipt promoted forbidden state: {key}")
    if evidence.get("included_in_stable_candidate_identity") is not False:
        raise ValueError("materialization evidence leaked into stable candidate identity")
    if evidence.get("included_in_stable_velocity_identity") is not False:
        raise ValueError("materialization evidence leaked into stable velocity identity")
    if migration.get("legacy_identity_reused_as_stable_identity") is not False:
        raise ValueError("legacy identity was reused as stable identity")
    return {"candidate_semantic_identity_sha256": candidate_id, "velocity_semantic_identity_sha256": velocity_id}


def validate_fresh_axial_receipt(
    receipt: dict[str, Any], *, expected_velocity_semantic_identity_sha256: str
) -> dict[str, Any]:
    """Validate the minimum contract for a future stable-ID-bound A9 rerun.

    The future producer may add fields, but it must state the stable callable
    identity explicitly; legacy whole-candidate hashes are insufficient.
    """
    bound = receipt.get("stable_velocity_semantic_identity_sha256")
    if bound is None and isinstance(receipt.get("candidate_identity"), dict):
        bound = receipt["candidate_identity"].get("stable_velocity_semantic_identity_sha256")
    bound = _sha256(bound, "stable_velocity_semantic_identity_sha256")

    measurement = receipt.get("measurement")
    truth = receipt.get("truth_boundary")
    protocol = receipt.get("protocol")
    if not all(isinstance(v, dict) for v in (measurement, truth, protocol)):
        raise ValueError("malformed fresh axial receipt")
    if protocol.get("source_numeric_targets_used") is not False:
        raise ValueError("fresh axial receipt used a public numerical target")
    if truth.get("visualization_ready") is not False or truth.get("pde_validated") is not False:
        raise ValueError("fresh axial receipt promoted visual/PDE state")
    metrics = measurement.get("metrics")
    if not isinstance(metrics, dict):
        raise ValueError("fresh axial receipt lacks metrics")
    aspect = float(metrics.get("full_aspect_ratio"))
    if not math.isfinite(aspect) or aspect <= 0.0:
        raise ValueError("fresh axial aspect must be finite and positive")
    return {
        "stable_velocity_semantic_identity_sha256": bound,
        "identity_matches": bound == expected_velocity_semantic_identity_sha256,
        "full_aspect_ratio": aspect,
    }


def classify(
    stable_run: dict[str, Any],
    stable_receipt: dict[str, Any] | None = None,
    fresh_axial_receipt: dict[str, Any] | None = None,
) -> dict[str, Any]:
    status, conclusion = validate_stable_run(stable_run)
    report: dict[str, Any] = {
        "task_id": TASK_ID,
        "prereg_issue": PREREG_ISSUE,
        "source_parent": {"pr": SOURCE_PARENT_PR, "head": SOURCE_PARENT_HEAD},
        "stable_identity_upstream": {
            "pr": STABLE_PR,
            "head": STABLE_HEAD,
            "run": STABLE_RUN,
            "status": status,
            "conclusion": conclusion,
        },
        "historical_failed_axial_evidence": {
            "pr": FAILED_AXIAL_PR,
            "run": FAILED_AXIAL_RUN,
            "descriptive_aspect": FAILED_AXIAL_ASPECT,
            "admitted_for_basis_routing": False,
        },
        "truth": dict(TRUTH),
    }

    if status != "completed" or conclusion != "success":
        report["decision"] = {
            "status": "blocked_stable_identity_ci",
            "stable_identity_admitted": False,
            "fresh_axial_rerun_required": True,
            "fresh_axial_measurement_admitted": False,
            "handoff_to_existing_axial_router_authorized": False,
            "existing_control_selection_authorized": False,
            "targeted_basis_preflight_authorized": False,
            "candidate_mutation_authorized": False,
        }
        return report

    if stable_receipt is None:
        raise ValueError("successful stable-identity run requires its exact receipt")
    stable_ids = validate_stable_receipt(stable_receipt)
    report["stable_semantic_identity"] = stable_ids

    if fresh_axial_receipt is None:
        report["decision"] = {
            "status": "stable_identity_ready_fresh_axial_rerun_required",
            "stable_identity_admitted": True,
            "fresh_axial_rerun_required": True,
            "fresh_axial_measurement_admitted": False,
            "handoff_to_existing_axial_router_authorized": False,
            "existing_control_selection_authorized": False,
            "targeted_basis_preflight_authorized": False,
            "candidate_mutation_authorized": False,
        }
        return report

    fresh = validate_fresh_axial_receipt(
        fresh_axial_receipt,
        expected_velocity_semantic_identity_sha256=stable_ids[
            "velocity_semantic_identity_sha256"
        ],
    )
    report["fresh_axial_measurement"] = fresh
    if not fresh["identity_matches"]:
        report["decision"] = {
            "status": "stable_identity_measurement_mismatch",
            "stable_identity_admitted": True,
            "fresh_axial_rerun_required": True,
            "fresh_axial_measurement_admitted": False,
            "handoff_to_existing_axial_router_authorized": False,
            "existing_control_selection_authorized": False,
            "targeted_basis_preflight_authorized": False,
            "candidate_mutation_authorized": False,
        }
        return report

    report["decision"] = {
        "status": "stable_identity_axial_measurement_ready_for_existing_router",
        "stable_identity_admitted": True,
        "fresh_axial_rerun_required": False,
        "fresh_axial_measurement_admitted": True,
        "handoff_to_existing_axial_router_authorized": True,
        "existing_control_selection_authorized": False,
        "targeted_basis_preflight_authorized": False,
        "candidate_mutation_authorized": False,
        "next": "feed the fresh stable-ID-bound axial coordinate into #1098/#1078; do not add a sixth basis unless that existing-control route demonstrates an obstruction",
    }
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stable-run-meta", type=Path, required=True)
    parser.add_argument("--stable-receipt", type=Path)
    parser.add_argument("--fresh-axial-receipt", type=Path)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    report = classify(
        _load(args.stable_run_meta),
        _load(args.stable_receipt) if args.stable_receipt else None,
        _load(args.fresh_axial_receipt) if args.fresh_axial_receipt else None,
    )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(report["decision"], indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
