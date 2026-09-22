"""Quarantine identity-drifted ST052 axial morphology evidence from basis routing.

Preregistered in Issue #1105 and stacked exactly on Agent-7 PR #1098.
A9 #1095 produced a valid numerical axial/vorticity-aspect measurement but its
exact-head workflow failed the same-identity assertion.  This module records
that result as descriptive, unadmitted evidence and blocks any control/basis
selection until semantic identity is repaired upstream.

It does not modify the candidate, repair identity construction, introduce an
OpenAI numerical target, select a coefficient, or add a basis.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

TASK_ID = "CR003-ST052M-AXIAL-EVIDENCE-ADMISSIBILITY-131"
PREREG_ISSUE = 1105
SOURCE_PARENT_PR = 1098
SOURCE_PARENT_HEAD = "cd7ef4531712bd369674aaf391e53759d97ec559"

A9_PR = 1095
A9_HEAD = "fe646b09b6f310068ca85f32b237891a71965f91"
A9_RUN = 35674391631
A9_ARTIFACT_ID = 10673062099
A9_ARTIFACT_NAME = "ST052-M-identity-bound-axial-vorticity-aspect"
A9_ARTIFACT_DIGEST = "sha256:d360434a480a4bd735d71c4256304dcb9145db29e22dace5d6ed5b37d0b7de85"
A9_SCHEMA = "st052-identity-bound-axial-vorticity-aspect/v1"
A9_TASK_ID = "CR-A9-104"

ADMITTED_CANDIDATE_SHA256 = "be6eeb6d15aa2e146e87a6feeaae3627ac9ce3c4be6b3fc19b4cf8393bf5973c"
ADMITTED_VELOCITY_SHA256 = "b0654bd8955644845fcb3789628299318b43bdb72bb467c284822235ffbeb902"
REPLAY_CANDIDATE_SHA256 = "5918fb93aca3f82ca63cf508510c17c02ddd151dd6eb5c9aaa2d2f253939e21e"
REPLAY_VELOCITY_SHA256 = "de28f374386e78679f79006e2dbdbdf238179267d13f8470796c44d7f36ebb42"
MEASUREMENT_SHA256 = "92802a6cd92184869756a673cf52f22648b74df8b95c8b5c8d555023df8839a4"
EXPECTED_AXIAL_RMS = 1.0237804347017856
EXPECTED_RADIAL_RMS = 0.9467960731269397
EXPECTED_ASPECT = 1.0813103938217585

TRUTH = {
    "candidate_velocity_changed": False,
    "basis_dimension_changed": False,
    "coefficient_selected": False,
    "public_openai_numeric_target_used": False,
    "held_out_pde_residual_evaluated": False,
    "measurement_observed": True,
    "measurement_admitted_for_basis_routing": False,
    "same_identity_gate_passed": False,
    "axial_proxy_no_deficit_conclusion_admitted": False,
    "existing_control_selection_authorized": False,
    "basis_preflight_authorized": False,
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


def _finite_positive(value: Any, name: str) -> float:
    x = float(value)
    if not math.isfinite(x) or x <= 0.0:
        raise ValueError(f"{name} must be finite and positive")
    return x


def validate_failed_run(run: dict[str, Any]) -> None:
    if int(run.get("id")) != A9_RUN:
        raise ValueError("unexpected A9 run id")
    if run.get("head_sha") != A9_HEAD:
        raise ValueError("unexpected A9 run head")
    if run.get("status") != "completed" or run.get("conclusion") != "failure":
        raise ValueError("this quarantine applies only to the exact completed/failed A9 run")


def validate_failed_receipt(receipt: dict[str, Any]) -> dict[str, Any]:
    if receipt.get("schema") != A9_SCHEMA or receipt.get("task_id") != A9_TASK_ID:
        raise ValueError("unexpected A9 receipt schema/task")

    identity = receipt.get("candidate_identity")
    protocol = receipt.get("protocol")
    measurement = receipt.get("measurement")
    public = receipt.get("public_observable")
    truth = receipt.get("truth_boundary")
    if not all(isinstance(v, dict) for v in (identity, protocol, measurement, public, truth)):
        raise ValueError("malformed A9 receipt")

    if identity.get("candidate_sha256") != REPLAY_CANDIDATE_SHA256:
        raise ValueError("unexpected replay candidate identity")
    if identity.get("velocity_identity_sha256") != REPLAY_VELOCITY_SHA256:
        raise ValueError("unexpected replay velocity identity")
    if REPLAY_CANDIDATE_SHA256 == ADMITTED_CANDIDATE_SHA256:
        raise ValueError("candidate identity unexpectedly did not drift")
    if REPLAY_VELOCITY_SHA256 == ADMITTED_VELOCITY_SHA256:
        raise ValueError("velocity identity unexpectedly did not drift")

    if float(protocol.get("time")) != 0.5:
        raise ValueError("time protocol drift")
    if int(protocol.get("grid_resolution")) != 25:
        raise ValueError("grid protocol drift")
    if [float(v) for v in protocol.get("box", [])] != [-2.0, 2.0]:
        raise ValueError("box protocol drift")
    if protocol.get("source_numeric_targets_used") is not False:
        raise ValueError("public numerical target laundering")
    if public.get("numerical_target") is not None:
        raise ValueError("public numerical target must remain null")
    if truth.get("axial_stretching_correspondence_verified") is not False:
        raise ValueError("receipt promoted axial correspondence")
    if truth.get("visualization_ready") is not False or truth.get("pde_validated") is not False:
        raise ValueError("receipt promoted visual/PDE state")

    if receipt.get("measurement_sha256") != MEASUREMENT_SHA256:
        raise ValueError("unexpected measurement digest")
    metrics = measurement.get("metrics")
    if not isinstance(metrics, dict):
        raise ValueError("missing aspect metrics")
    axial = _finite_positive(metrics.get("full_axial_rms"), "full_axial_rms")
    radial = _finite_positive(metrics.get("full_radial_rms"), "full_radial_rms")
    aspect = _finite_positive(metrics.get("full_aspect_ratio"), "full_aspect_ratio")
    elongated = metrics.get("axially_elongated_proxy")
    if not math.isclose(axial, EXPECTED_AXIAL_RMS, rel_tol=0.0, abs_tol=2e-15):
        raise ValueError("observed axial RMS drift")
    if not math.isclose(radial, EXPECTED_RADIAL_RMS, rel_tol=0.0, abs_tol=2e-15):
        raise ValueError("observed radial RMS drift")
    if not math.isclose(aspect, EXPECTED_ASPECT, rel_tol=0.0, abs_tol=2e-15):
        raise ValueError("observed aspect drift")
    if not math.isclose(aspect, axial / radial, rel_tol=0.0, abs_tol=2e-12):
        raise ValueError("aspect arithmetic drift")
    if elongated is not True or elongated is not (axial > radial):
        raise ValueError("observed axial proxy drift")

    return {
        "full_axial_rms": axial,
        "full_radial_rms": radial,
        "full_aspect_ratio": aspect,
        "axially_elongated_proxy": True,
        "measurement_sha256": MEASUREMENT_SHA256,
    }


def classify(run: dict[str, Any], receipt: dict[str, Any]) -> dict[str, Any]:
    validate_failed_run(run)
    measurement = validate_failed_receipt(receipt)
    return {
        "task_id": TASK_ID,
        "prereg_issue": PREREG_ISSUE,
        "source_parent": {"pr": SOURCE_PARENT_PR, "head": SOURCE_PARENT_HEAD},
        "exact_upstream": {
            "pr": A9_PR,
            "head": A9_HEAD,
            "run": A9_RUN,
            "run_status": "completed",
            "run_conclusion": "failure",
            "artifact_id": A9_ARTIFACT_ID,
            "artifact_name": A9_ARTIFACT_NAME,
            "artifact_digest": A9_ARTIFACT_DIGEST,
        },
        "identity_admissibility": {
            "previously_admitted_candidate_sha256": ADMITTED_CANDIDATE_SHA256,
            "replay_candidate_sha256": REPLAY_CANDIDATE_SHA256,
            "previously_admitted_velocity_identity_sha256": ADMITTED_VELOCITY_SHA256,
            "replay_velocity_identity_sha256": REPLAY_VELOCITY_SHA256,
            "candidate_identity_matches": False,
            "velocity_identity_matches": False,
            "same_identity_gate_passed": False,
            "classification": "identity_drifted_measurement_observed_not_admitted",
        },
        "descriptive_unadmitted_measurement": measurement,
        "decision": {
            "status": "blocked_identity_admissibility",
            "action": "repair_semantic_identity_before_basis_decision",
            "measurement_observed": True,
            "measurement_admitted_for_basis_routing": False,
            "axial_proxy_no_deficit_conclusion_admitted": False,
            "recommended_existing_control": None,
            "coefficient_magnitude_selected": None,
            "existing_control_selection_authorized": False,
            "targeted_basis_preflight_authorized": False,
            "candidate_mutation_authorized": False,
            "minimum_upstream_repair": (
                "Separate stable semantic candidate/velocity identity (source replay identity + temporal transform semantics + exact runtime/acceptance contract) "
                "from regenerated materialization/evidence checksums; continue recording/checking those checksums and fail closed on true semantic drift."
            ),
            "hash_replacement_alone_is_sufficient": False,
        },
        "truth": dict(TRUTH),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-meta", type=Path, required=True)
    parser.add_argument("--axial-receipt", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    report = classify(_load(args.run_meta), _load(args.axial_receipt))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(report["decision"], indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
