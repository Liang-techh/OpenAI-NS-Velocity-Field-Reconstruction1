"""Route frozen-ST052 axial proxy evidence to the minimum existing morphology control.

Preregistered in Issue #1097. This module combines two already-defined receipts:
A9's same-identity frozen-ST052 axial/vorticity aspect measurement and Agent-7's
five-control defect-coordinate Jacobian. It does not alter the candidate,
introduce an OpenAI numerical target, choose a coefficient magnitude, or claim
visual/PDE correspondence.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

TASK_ID = "CR003-ST052M-AXIAL-PROXY-CONTROL-ROUTER-130"
PREREG_ISSUE = 1097
SOURCE_PARENT_PR = 1087
SOURCE_PARENT_HEAD = "88db8fde63a637f184d8cd7544e5997eef16a49d"

A9_SOURCE_PR = 1095
A9_SOURCE_HEAD = "fe646b09b6f310068ca85f32b237891a71965f91"
A9_SOURCE_RUN = 35674391631
A9_ARTIFACT_NAME = "ST052-M-identity-bound-axial-vorticity-aspect"
A9_SCHEMA = "st052-identity-bound-axial-vorticity-aspect/v1"
A9_TASK_ID = "CR-A9-104"

A7_SOURCE_PR = 1078
A7_SOURCE_HEAD = "849713ba95cde8002fc1002d8b6856ea3f1dc264"
A7_SOURCE_RUN = 35666562373
A7_ARTIFACT_NAME = "agent7-st052m-defect-coordinate-controllability"
A7_TASK_ID = "CR003-ST052M-DEFECT-COORDINATE-CONTROLLABILITY-128"
A7_PREREG_ISSUE = 1077
A7_PARENT_PR = 1069
A7_PARENT_HEAD = "a564fbc8bfa6a5c630b5839fa208e41cb7bc72d7"

EXPECTED_CANDIDATE_ID = "ST052-M-linear-temporal-child-v1"
EXPECTED_CANDIDATE_SHA256 = "be6eeb6d15aa2e146e87a6feeaae3627ac9ce3c4be6b3fc19b4cf8393bf5973c"
EXPECTED_VELOCITY_IDENTITY_SHA256 = "b0654bd8955644845fcb3789628299318b43bdb72bb467c284822235ffbeb902"

CHANNELS = (
    "amplitude",
    "radial_shape",
    "axial_turnover",
    "temporal_curvature",
    "toroidal_swirl",
)
COORDINATE = "axial_vorticity_aspect_t050"

TRUTH = {
    "candidate_velocity_changed": False,
    "basis_dimension_changed": False,
    "coefficient_magnitude_selected": False,
    "public_openai_numeric_target_used": False,
    "renderer_or_camera_fit_used": False,
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


def _finite_positive(value: Any, name: str) -> float:
    x = float(value)
    if not math.isfinite(x) or x <= 0.0:
        raise ValueError(f"{name} must be finite and positive")
    return x


def _finite(value: Any, name: str) -> float:
    x = float(value)
    if not math.isfinite(x):
        raise ValueError(f"{name} must be finite")
    return x


def _load_json(path: str | Path) -> dict[str, Any]:
    obj = json.loads(Path(path).read_text(encoding="utf-8"))
    if not isinstance(obj, dict):
        raise ValueError("receipt root must be an object")
    return obj


def validate_axial_receipt(receipt: dict[str, Any]) -> dict[str, Any]:
    if receipt.get("schema") != A9_SCHEMA or receipt.get("task_id") != A9_TASK_ID:
        raise ValueError("unexpected A9 axial receipt schema/task")
    identity = receipt.get("candidate_identity")
    protocol = receipt.get("protocol")
    public = receipt.get("public_observable")
    measurement = receipt.get("measurement")
    truth = receipt.get("truth_boundary")
    if not all(isinstance(v, dict) for v in (identity, protocol, public, measurement, truth)):
        raise ValueError("malformed A9 axial receipt")
    expected_identity = {
        "candidate_id": EXPECTED_CANDIDATE_ID,
        "candidate_sha256": EXPECTED_CANDIDATE_SHA256,
        "velocity_identity_sha256": EXPECTED_VELOCITY_IDENTITY_SHA256,
    }
    if {k: identity.get(k) for k in expected_identity} != expected_identity:
        raise ValueError("frozen ST052 identity drift")
    if float(protocol.get("time")) != 0.5:
        raise ValueError("A9 time protocol drift")
    if int(protocol.get("grid_resolution")) != 25:
        raise ValueError("A9 grid protocol drift")
    if [float(v) for v in protocol.get("box", [])] != [-2.0, 2.0]:
        raise ValueError("A9 box protocol drift")
    if protocol.get("source_numeric_targets_used") is not False:
        raise ValueError("public numerical target laundering")
    if protocol.get("renderer_or_camera_used") is not False or protocol.get("pixel_loss_used") is not False:
        raise ValueError("renderer/pixel target laundering")
    if public.get("numerical_target") is not None:
        raise ValueError("A9 public numerical target must remain null")
    if truth.get("axial_stretching_correspondence_verified") is not False:
        raise ValueError("A9 receipt promoted axial correspondence")
    if truth.get("visualization_ready") is not False or truth.get("pde_validated") is not False:
        raise ValueError("A9 receipt promoted visual/PDE state")

    metrics = measurement.get("metrics")
    if not isinstance(metrics, dict):
        raise ValueError("missing A9 axial metrics")
    axial = _finite_positive(metrics.get("full_axial_rms"), "full_axial_rms")
    radial = _finite_positive(metrics.get("full_radial_rms"), "full_radial_rms")
    aspect = _finite_positive(metrics.get("full_aspect_ratio"), "full_aspect_ratio")
    if not math.isclose(aspect, axial / radial, rel_tol=0.0, abs_tol=2.0e-12):
        raise ValueError("A9 aspect arithmetic drift")
    elongated = metrics.get("axially_elongated_proxy")
    if not isinstance(elongated, bool) or elongated is not (axial > radial):
        raise ValueError("A9 axial proxy drift")
    return {
        "candidate_id": EXPECTED_CANDIDATE_ID,
        "candidate_sha256": EXPECTED_CANDIDATE_SHA256,
        "velocity_identity_sha256": EXPECTED_VELOCITY_IDENTITY_SHA256,
        "full_axial_rms": axial,
        "full_radial_rms": radial,
        "full_aspect_ratio": aspect,
        "axially_elongated_proxy": elongated,
        "measurement_sha256": receipt.get("measurement_sha256"),
    }


def validate_controllability_receipt(receipt: dict[str, Any]) -> dict[str, Any]:
    if receipt.get("task_id") != A7_TASK_ID or receipt.get("prereg_issue") != A7_PREREG_ISSUE:
        raise ValueError("unexpected A7 controllability receipt identity")
    if receipt.get("source_parent") != {"pr": A7_PARENT_PR, "head": A7_PARENT_HEAD}:
        raise ValueError("A7 controllability parent drift")
    protocol = receipt.get("frozen_protocol")
    audit = receipt.get("defect_coordinate_controllability")
    decision = receipt.get("decision")
    if not all(isinstance(v, dict) for v in (protocol, audit, decision)):
        raise ValueError("malformed A7 controllability receipt")
    if tuple(protocol.get("channels", ())) != CHANNELS:
        raise ValueError("A7 five-channel order drift")
    if COORDINATE not in tuple(protocol.get("defect_coordinates", ())):
        raise ValueError("A7 axial coordinate missing")
    if protocol.get("new_basis_dimension") != 0 or protocol.get("coefficient_selected") is not None:
        raise ValueError("A7 parent unexpectedly changed basis/coefficient")
    if protocol.get("public_openai_numeric_target") is not None:
        raise ValueError("A7 public numerical target laundering")
    if decision.get("candidate_mutation_authorized_by_this_audit") is not False:
        raise ValueError("A7 parent unexpectedly authorized candidate mutation")
    if float(decision.get("direct_visualization_fingerprint_improvement")) != 0.0:
        raise ValueError("A7 parent unexpectedly claims visual improvement")

    rows = audit.get("coordinates")
    if not isinstance(rows, dict) or not isinstance(rows.get(COORDINATE), dict):
        raise ValueError("A7 axial coordinate row missing")
    row = rows[COORDINATE]
    derivatives = row.get("fine_raw_derivative_by_channel")
    if not isinstance(derivatives, dict) or tuple(derivatives.keys()) != CHANNELS:
        raise ValueError("A7 axial derivative channel order drift")
    derivs = {channel: _finite(derivatives[channel], f"derivative[{channel}]") for channel in CHANNELS}
    row_l2 = _finite(row.get("fine_raw_row_l2"), "fine_raw_row_l2")
    if row_l2 < 0.0:
        raise ValueError("negative A7 row norm")
    unresponsive = row.get("structurally_unresponsive_at_numeric_floor")
    if not isinstance(unresponsive, bool):
        raise ValueError("invalid A7 unresponsive flag")
    dominant = row.get("dominant_aligned_existing_channel")
    if dominant is not None and dominant not in CHANNELS:
        raise ValueError("invalid A7 dominant channel")
    geometry_pass = audit.get("compressed_routing_geometry_gate_passes")
    stability_pass = audit.get("all_derivative_stability_gates_pass")
    if not isinstance(geometry_pass, bool) or not isinstance(stability_pass, bool):
        raise ValueError("invalid A7 routing/stability gate")
    condition = audit.get("normalized_condition")
    if condition is not None:
        condition = _finite(condition, "normalized_condition")
    rank = int(audit.get("normalized_column_rank"))
    if not 0 <= rank <= 5:
        raise ValueError("invalid A7 normalized rank")
    return {
        "rank": rank,
        "condition": condition,
        "routing_geometry_gate_passes": geometry_pass,
        "derivative_stability_gates_pass": stability_pass,
        "row_l2": row_l2,
        "structurally_unresponsive_at_numeric_floor": unresponsive,
        "dominant_aligned_existing_channel": dominant,
        "fine_raw_derivative_by_channel": derivs,
    }


def route(axial_receipt: dict[str, Any], controllability_receipt: dict[str, Any]) -> dict[str, Any]:
    axial = validate_axial_receipt(axial_receipt)
    control = validate_controllability_receipt(controllability_receipt)
    deficit = not axial["axially_elongated_proxy"]

    action: str
    recommended_control: str | None = None
    recommended_direction: str | None = None
    existing_control_test_justified = False
    targeted_basis_preflight_justified = False

    if not deficit:
        action = "no_axial_proxy_deficit_identified"
        explanation = (
            "The frozen candidate is axially elongated under the preregistered candidate-side "
            "RMS proxy. Without a public numerical aspect target, this does not prove source "
            "matching and does not justify changing an existing control or adding a basis."
        )
    elif control["structurally_unresponsive_at_numeric_floor"]:
        action = "targeted_axial_tip_basis_preflight_only"
        targeted_basis_preflight_justified = True
        explanation = (
            "The candidate-side axial proxy is not elongated and the existing five-control "
            "axial-aspect row is numerically unresponsive. One narrowly targeted axial/tip "
            "basis preflight is justified for testing only; candidate promotion is not."
        )
    elif not (control["routing_geometry_gate_passes"] and control["derivative_stability_gates_pass"]):
        action = "reparameterize_existing_controls_before_basis_growth"
        explanation = (
            "The candidate-side axial proxy is not elongated, but the existing-control Jacobian "
            "is not sufficiently conditioned/stable for a reliable minimum-control choice. "
            "Repair or reparameterize the current control geometry before adding a basis."
        )
    else:
        dominant = control["dominant_aligned_existing_channel"]
        if dominant is None:
            raise ValueError("stable responsive A7 row has no dominant existing channel")
        derivative = control["fine_raw_derivative_by_channel"][dominant]
        if derivative == 0.0:
            raise ValueError("dominant existing channel has zero axial-aspect derivative")
        recommended_control = dominant
        recommended_direction = "increase" if derivative > 0.0 else "decrease"
        existing_control_test_justified = True
        action = "test_existing_control_first"
        explanation = (
            "The candidate-side axial proxy is not elongated and the existing five-control "
            "Jacobian is responsive, conditioned and locally stable. Test the dominant existing "
            "control in the signed direction that increases axial aspect before any basis growth."
        )

    return {
        "task_id": TASK_ID,
        "prereg_issue": PREREG_ISSUE,
        "source_parent": {"pr": SOURCE_PARENT_PR, "head": SOURCE_PARENT_HEAD},
        "upstream_evidence": {
            "a9_axial": {"pr": A9_SOURCE_PR, "head": A9_SOURCE_HEAD, "run": A9_SOURCE_RUN, "artifact": A9_ARTIFACT_NAME},
            "a7_controllability": {"pr": A7_SOURCE_PR, "head": A7_SOURCE_HEAD, "run": A7_SOURCE_RUN, "artifact": A7_ARTIFACT_NAME},
        },
        "frozen_protocol": {
            "candidate_id": EXPECTED_CANDIDATE_ID,
            "candidate_sha256": EXPECTED_CANDIDATE_SHA256,
            "velocity_identity_sha256": EXPECTED_VELOCITY_IDENTITY_SHA256,
            "coordinate": COORDINATE,
            "controls": list(CHANNELS),
            "geometric_aspect_equality_boundary": 1.0,
            "public_openai_numeric_aspect_target": None,
        },
        "axial_proxy": axial,
        "existing_control_response": control,
        "decision": {
            "status": "ready",
            "candidate_side_axial_proxy_deficit": deficit,
            "action": action,
            "recommended_existing_control": recommended_control,
            "recommended_control_direction_to_increase_axial_aspect": recommended_direction,
            "coefficient_magnitude_selected": None,
            "bounded_existing_control_child_test_justified": existing_control_test_justified,
            "targeted_axial_tip_basis_preflight_justified": targeted_basis_preflight_justified,
            "candidate_mutation_authorized_by_this_router": False,
            "basis_growth_authorized_for_candidate_promotion": False,
            "explanation": explanation,
        },
        "truth": dict(TRUTH),
    }


def blocked_report(upstream_status: dict[str, Any]) -> dict[str, Any]:
    return {
        "task_id": TASK_ID,
        "prereg_issue": PREREG_ISSUE,
        "source_parent": {"pr": SOURCE_PARENT_PR, "head": SOURCE_PARENT_HEAD},
        "upstream_status": upstream_status,
        "decision": {
            "status": "blocked_upstream",
            "action": "wait_for_exact_upstream_receipts",
            "recommended_existing_control": None,
            "recommended_control_direction_to_increase_axial_aspect": None,
            "coefficient_magnitude_selected": None,
            "bounded_existing_control_child_test_justified": False,
            "targeted_axial_tip_basis_preflight_justified": False,
            "candidate_mutation_authorized_by_this_router": False,
            "basis_growth_authorized_for_candidate_promotion": False,
        },
        "truth": dict(TRUTH),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--axial-receipt", type=Path, required=True)
    parser.add_argument("--controllability-receipt", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    report = route(_load_json(args.axial_receipt), _load_json(args.controllability_receipt))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(report["decision"], indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
