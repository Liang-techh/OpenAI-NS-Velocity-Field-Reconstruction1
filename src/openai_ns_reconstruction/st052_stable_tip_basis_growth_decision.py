"""Fail-closed Agent-7 routing from the exact successful ST052 tip-taper receipt.

Preregistered in issue #1239.  This module does not mutate the candidate.  It
only decides whether the coarse candidate-side defect "tips are not narrower
than the core" is present, and if present whether the already-preregistered
existing-control receipt is available.  Visual resemblance is not PDE validity
or source correspondence.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Any

TASK_ID = "CR003-ST052M-STABLE-TIP-DECISION-146"
PREREG_ISSUE = 1239
SOURCE_PR = 1195
SOURCE_HEAD = "cb7bd1e4374ddb84cd6273307024c940cf6319a5"
SOURCE_RUN = 35726365645
SOURCE_ARTIFACT_ID = 10702308370
SOURCE_ARTIFACT_DIGEST = "sha256:85868e2c371c54ee1c805e62bb1d9e6e5efaa4711c6377aed1f2a6a093e45d68"
SOURCE_TASK_ID = "CR003-ST052M-STABLE-TIP-TAPER-PROXY-141"
SOURCE_PREREG_ISSUE = 1194
EXPECTED_CANDIDATE_IDENTITY = "2b2e743e72e31f69ad4232340f52d48b5759ac8aad1dfb735708112276c0b9c4"
EXPECTED_VELOCITY_IDENTITY = "7738013a4fe4cacbfebfb4e8822bb4ee815b1d3bab744b90f9473ad9424f4f41"
CONTROL_PR = 1078
CONTROL_HEAD = "849713ba95cde8002fc1002d8b6856ea3f1dc264"
CONTROL_RUN = 35666562373


def canonical_sha256(value: Any) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode()
    return hashlib.sha256(payload).hexdigest()


def validate_source_receipt(receipt: dict[str, Any]) -> dict[str, Any]:
    if receipt.get("task_id") != SOURCE_TASK_ID or receipt.get("prereg_issue") != SOURCE_PREREG_ISSUE:
        raise ValueError("unexpected source task/preregistration")
    unsigned = dict(receipt)
    recorded = unsigned.pop("receipt_sha256", None)
    if not isinstance(recorded, str) or recorded != canonical_sha256(unsigned):
        raise ValueError("source receipt checksum mismatch")

    ident = receipt.get("stable_candidate_identity")
    if not isinstance(ident, dict):
        raise ValueError("missing stable identity")
    if ident.get("candidate_semantic_identity_sha256") != EXPECTED_CANDIDATE_IDENTITY:
        raise ValueError("candidate semantic identity drift")
    if ident.get("velocity_semantic_identity_sha256") != EXPECTED_VELOCITY_IDENTITY:
        raise ValueError("velocity semantic identity drift")

    public = receipt.get("public_context")
    if not isinstance(public, dict) or public.get("public_numerical_target") is not None:
        raise ValueError("public numerical target boundary drift")
    truth = receipt.get("truth_boundary")
    if not isinstance(truth, dict):
        raise ValueError("missing source truth boundary")
    required_false = (
        "candidate_velocity_changed", "basis_dimension_changed", "coefficient_selected",
        "visualization_ready", "source_correspondence_verified", "pde_validated",
        "paper_exact", "openai_field_identified", "blowup_proved",
    )
    if any(truth.get(k) is not False for k in required_false):
        raise ValueError("source truth boundary promoted")
    if float(truth.get("direct_visualization_fingerprint_improvement", math.nan)) != 0.0:
        raise ValueError("source direct-visualization delta drift")

    measurement = receipt.get("measurement")
    if not isinstance(measurement, dict) or receipt.get("measurement_sha256") != canonical_sha256(measurement):
        raise ValueError("measurement checksum mismatch")
    metrics = measurement.get("metrics")
    if not isinstance(metrics, dict):
        raise ValueError("missing tip metrics")
    keys = ("tip_radial_rms_upper", "tip_radial_rms_lower", "central_radial_rms",
            "worst_tip_radial_rms", "worst_tip_to_central_ratio")
    vals = {k: float(metrics.get(k, math.nan)) for k in keys}
    if not all(math.isfinite(v) and v > 0.0 for v in vals.values()):
        raise ValueError("nonfinite/degenerate tip metrics")
    worst = max(vals["tip_radial_rms_upper"], vals["tip_radial_rms_lower"])
    if not math.isclose(vals["worst_tip_radial_rms"], worst, rel_tol=0.0, abs_tol=2e-15):
        raise ValueError("worst-tip arithmetic drift")
    ratio = worst / vals["central_radial_rms"]
    if not math.isclose(vals["worst_tip_to_central_ratio"], ratio, rel_tol=0.0, abs_tol=2e-15):
        raise ValueError("tip/core ratio arithmetic drift")
    tapered = worst < vals["central_radial_rms"]
    if metrics.get("both_tips_radially_tapered_proxy") is not tapered:
        raise ValueError("tip-taper proxy drift")
    return metrics


def make_decision(receipt: dict[str, Any], control_run: dict[str, Any]) -> dict[str, Any]:
    metrics = validate_source_receipt(receipt)
    if control_run.get("id") != CONTROL_RUN or control_run.get("head_sha") != CONTROL_HEAD:
        raise ValueError("existing-control run identity drift")
    status = control_run.get("status")
    conclusion = control_run.get("conclusion")
    tapered = bool(metrics["both_tips_radially_tapered_proxy"])

    if tapered:
        classification = "coarse_blunt_tip_basis_trigger_closed"
        next_action = "require_distinct_stable_identity_bound_trajectory_or_return_flow_discrepancy"
    elif status == "completed" and conclusion == "success":
        classification = "tip_taper_deficit_existing_control_receipt_ready"
        next_action = "consume_pr1078_tip_radial_thickness_t050_before_any_new_basis"
    else:
        classification = "tip_taper_deficit_existing_control_receipt_blocked"
        next_action = "wait_for_exact_pr1078_tip_control_receipt_before_any_new_basis"

    decision = {
        "schema": "st052-stable-tip-basis-growth-decision/v1",
        "task_id": TASK_ID,
        "prereg_issue": PREREG_ISSUE,
        "source": {
            "pr": SOURCE_PR, "head": SOURCE_HEAD, "run": SOURCE_RUN,
            "artifact_id": SOURCE_ARTIFACT_ID, "artifact_digest": SOURCE_ARTIFACT_DIGEST,
            "receipt_sha256": receipt["receipt_sha256"],
        },
        "stable_candidate_identity": receipt["stable_candidate_identity"],
        "tip_metrics": metrics,
        "existing_control_dependency": {
            "pr": CONTROL_PR, "head": CONTROL_HEAD, "run": CONTROL_RUN,
            "status": status, "conclusion": conclusion,
        },
        "routing": {
            "classification": classification,
            "new_tip_basis_authorized": False,
            "axial_turnover_change_authorized": False,
            "coefficient_selection_authorized": False,
            "next_minimal_action": next_action,
        },
        "truth_boundary": {
            "candidate_velocity_changed": False,
            "basis_dimension_changed": False,
            "coefficient_selected": False,
            "held_out_pde_residual_evaluated": False,
            "direct_visualization_fingerprint_improvement": 0.0,
            "closer_visualization_delivery_established": False,
            "visualization_ready": False,
            "source_correspondence_verified": False,
            "pde_validated": False,
            "paper_exact": False,
            "openai_field_identified": False,
        },
    }
    decision["decision_sha256"] = canonical_sha256(decision)
    return decision


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--receipt", required=True)
    p.add_argument("--control-run", required=True)
    p.add_argument("--out", required=True)
    args = p.parse_args()
    receipt = json.loads(Path(args.receipt).read_text())
    control_run = json.loads(Path(args.control_run).read_text())
    decision = make_decision(receipt, control_run)
    out = Path(args.out); out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(decision, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"classification": decision["routing"]["classification"],
                      "ratio": decision["tip_metrics"]["worst_tip_to_central_ratio"],
                      "control_status": decision["existing_control_dependency"]["status"],
                      "control_conclusion": decision["existing_control_dependency"]["conclusion"]}, sort_keys=True))


if __name__ == "__main__":
    main()
