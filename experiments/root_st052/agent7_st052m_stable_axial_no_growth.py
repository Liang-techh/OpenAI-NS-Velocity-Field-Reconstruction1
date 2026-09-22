"""Close axial-only basis growth against fresh stable-identity ST052 evidence.

Preregistered in Issue #1186 and stacked exactly on Agent-7 PR #1178 head.
This module consumes only the merged A9 #1163 axial/vorticity receipt.  It does
not modify the candidate, choose a coefficient, add a basis, or infer public
visual correspondence from the autonomous candidate-side aspect proxy.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Any

TASK_ID = "CR003-ST052M-STABLE-AXIAL-NO-GROWTH-140"
PREREG_ISSUE = 1186
STACK_BASE_PR = 1178
STACK_BASE_HEAD = "16a52f0733b50c9a2d0ef05222479b410f71f694"

A9_SOURCE_PR = 1163
A9_SOURCE_HEAD = "8fec2517677adb37754df8c0164330d6f346db57"
A9_MAIN_MERGE = "b0093e04039946923de10e421e001cd363213dcc"
A9_SOURCE_RUN = 35708397160
A9_ARTIFACT = "ST052-M-stable-identity-bound-axial-vorticity-aspect"
A9_TASK_ID = "CR-A9-109"
A9_SCHEMA = "st052-stable-identity-bound-axial-vorticity-aspect/v1"
EXPECTED_CANDIDATE_ID = "ST052-M-linear-temporal-child-v1"
EXPECTED_CANDIDATE_SEMANTIC_ID = (
    "2b2e743e72e31f69ad4232340f52d48b5759ac8aad1dfb735708112276c0b9c4"
)
EXPECTED_VELOCITY_SEMANTIC_ID = (
    "7738013a4fe4cacbfebfb4e8822bb4ee815b1d3bab744b90f9473ad9424f4f41"
)
EXPECTED_MEASUREMENT_SHA256 = (
    "ee0439ef09569aec9d0288b41663e73c1dc4d3473adf2debc140f121560e2a3b"
)


def _canonical_sha256(value: Any) -> str:
    payload = json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
        allow_nan=False,
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _finite_positive(value: Any, label: str) -> float:
    _require(isinstance(value, (int, float)) and not isinstance(value, bool), f"{label} must be numeric")
    x = float(value)
    _require(math.isfinite(x) and x > 0.0, f"{label} must be finite and positive")
    return x


def validate_upstream_receipt(receipt: dict[str, Any]) -> dict[str, float | bool]:
    """Authenticate exactly the preregistered A9 #1163 scientific boundary."""

    _require(receipt.get("task_id") == A9_TASK_ID, "unexpected upstream task_id")
    _require(receipt.get("schema") == A9_SCHEMA, "unexpected upstream schema")

    identity = receipt.get("stable_candidate_identity")
    _require(isinstance(identity, dict), "missing stable_candidate_identity")
    _require(identity.get("candidate_id") == EXPECTED_CANDIDATE_ID, "candidate id drift")
    _require(
        identity.get("candidate_semantic_identity_sha256") == EXPECTED_CANDIDATE_SEMANTIC_ID,
        "candidate semantic identity drift",
    )
    _require(
        identity.get("velocity_semantic_identity_sha256") == EXPECTED_VELOCITY_SEMANTIC_ID,
        "velocity semantic identity drift",
    )

    measurement = receipt.get("measurement")
    _require(isinstance(measurement, dict), "missing measurement")
    measurement_sha = receipt.get("measurement_sha256")
    _require(measurement_sha == EXPECTED_MEASUREMENT_SHA256, "unexpected measurement SHA")
    _require(_canonical_sha256(measurement) == measurement_sha, "measurement content/hash mismatch")

    protocol = receipt.get("protocol")
    _require(isinstance(protocol, dict), "missing protocol")
    _require(protocol.get("time") == 0.5, "diagnostic time drift")
    _require(protocol.get("grid_resolution") == 25, "grid resolution drift")
    _require(protocol.get("box") == [-2.0, 2.0], "box drift")
    _require(protocol.get("cartesian_curl_order") == 2, "curl order drift")
    _require(protocol.get("integration_rule") == "tensor_product_trapezoid", "integration drift")
    _require(protocol.get("source_numeric_targets_used") is False, "public numerical target introduced")
    _require(protocol.get("renderer_or_camera_used") is False, "renderer/camera unexpectedly used")
    _require(protocol.get("pixel_loss_used") is False, "pixel loss unexpectedly used")

    public = receipt.get("public_observable")
    _require(isinstance(public, dict), "missing public observable boundary")
    _require(public.get("id") == "axial_stretching_trajectories", "public observable drift")
    _require(public.get("numerical_target") is None, "public axial numerical target must stay null")
    _require(public.get("source_contract_admitted_by_this_increment") is False, "source contract promotion")

    truth = receipt.get("truth_boundary")
    _require(isinstance(truth, dict), "missing truth boundary")
    _require(truth.get("stable_semantic_identity_bound") is True, "stable identity not admitted")
    _require(truth.get("axial_vorticity_aspect_measured") is True, "axial aspect not measured")
    _require(truth.get("historical_cr_a9_104_receipt_reused") is False, "historical stale receipt reused")
    for key in (
        "candidate_changed",
        "velocity_coefficients_changed",
        "pressure_changed",
        "forcing_changed",
        "scientific_threshold_changed",
        "public_source_numeric_targets_used",
        "renderer_or_camera_used",
        "pixel_loss_used",
        "axial_stretching_correspondence_verified",
        "visual_correspondence_verified",
        "source_correspondence_verified",
        "visualization_ready",
        "velocity_export_ready",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    ):
        _require(truth.get(key) is False, f"forbidden truth promotion: {key}")

    metrics = measurement.get("metrics")
    _require(isinstance(metrics, dict), "missing measurement metrics")
    axial = _finite_positive(metrics.get("full_axial_rms"), "full_axial_rms")
    radial = _finite_positive(metrics.get("full_radial_rms"), "full_radial_rms")
    aspect = _finite_positive(metrics.get("full_aspect_ratio"), "full_aspect_ratio")
    enstrophy = _finite_positive(metrics.get("full_enstrophy_trapezoid"), "full_enstrophy_trapezoid")
    proxy = metrics.get("axially_elongated_proxy")
    _require(isinstance(proxy, bool), "axially_elongated_proxy must be bool")
    _require(math.isclose(aspect, axial / radial, rel_tol=2.0e-14, abs_tol=2.0e-14), "aspect arithmetic drift")
    _require(proxy == (aspect > 1.0), "proxy/aspect inconsistency")

    return {
        "full_axial_rms": axial,
        "full_radial_rms": radial,
        "full_aspect_ratio": aspect,
        "full_enstrophy_trapezoid": enstrophy,
        "axially_elongated_proxy": proxy,
    }


def build_report(receipt: dict[str, Any]) -> dict[str, Any]:
    metrics = validate_upstream_receipt(receipt)
    elongated = bool(metrics["axially_elongated_proxy"])

    if elongated:
        status = "no_candidate_side_axial_extent_deficit_established"
        action = "hold_axial_specific_growth_until_distinct_stable_identity_bound_defect"
        existing_control_router_required = False
        rationale = (
            "The admitted candidate-side autonomous aspect proxy is already axial>radial. "
            "That removes axial extent alone as a present basis-growth trigger; a distinct "
            "stable-identity-bound tip/return-flow or trajectory defect is required before "
            "an axial-specific mutation is reconsidered."
        )
    else:
        status = "candidate_side_axial_proxy_deficit_requires_existing_control_routing"
        action = "route_to_existing_five_control_controllability_before_basis_preflight"
        existing_control_router_required = True
        rationale = (
            "A candidate-side proxy deficit would first be routed through the existing five "
            "controls. This increment never promotes a sixth basis or chooses a coefficient."
        )

    return {
        "task_id": TASK_ID,
        "prereg_issue": PREREG_ISSUE,
        "stack_base_pr": STACK_BASE_PR,
        "stack_base_head": STACK_BASE_HEAD,
        "upstream": {
            "pr": A9_SOURCE_PR,
            "head": A9_SOURCE_HEAD,
            "main_merge": A9_MAIN_MERGE,
            "run": A9_SOURCE_RUN,
            "artifact": A9_ARTIFACT,
            "candidate_id": EXPECTED_CANDIDATE_ID,
            "candidate_semantic_identity_sha256": EXPECTED_CANDIDATE_SEMANTIC_ID,
            "velocity_semantic_identity_sha256": EXPECTED_VELOCITY_SEMANTIC_ID,
            "measurement_sha256": EXPECTED_MEASUREMENT_SHA256,
        },
        "measurement": metrics,
        "decision": {
            "status": status,
            "action": action,
            "rationale": rationale,
            "existing_control_router_required": existing_control_router_required,
            "axial_specific_basis_growth_authorized": False,
            "axial_turnover_adjustment_authorized": False,
            "candidate_mutation_authorized": False,
            "coefficient_selected": None,
            "distinct_stable_identity_bound_defect_required_for_future_axial_mutation": elongated,
        },
        "truth": {
            "candidate_velocity_changed": False,
            "basis_dimension_changed": False,
            "coefficient_selected": False,
            "pressure_or_force_changed": False,
            "free_residual_force_used": False,
            "amplitude_collapse_used": False,
            "held_out_pde_residual_evaluated": False,
            "public_openai_numeric_target_used": False,
            "visual_correspondence_verified": False,
            "source_correspondence_verified": False,
            "visualization_ready": False,
            "pde_validated": False,
            "paper_exact": False,
            "openai_field_identified": False,
            "blowup_proved": False,
            "direct_visualization_fingerprint_improvement": 0.0,
        },
        "interpretation": (
            "This closes only an axial-extent basis-growth trigger under one autonomous "
            "candidate-side proxy. It does not assert that the public OpenAI visualization "
            "has aspect>1 numerically, nor that the frozen field visually corresponds to it."
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--axial-receipt", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    receipt = json.loads(args.axial_receipt.read_text(encoding="utf-8"))
    report = build_report(receipt)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(report["decision"], sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
