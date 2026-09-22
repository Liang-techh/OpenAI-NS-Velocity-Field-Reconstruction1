"""Audit public-observable evidence coverage for frozen ST052-M morphology.

Preregistered in issue #1086 and stacked exactly on Agent-7 PR #1078 head
849713ba95cde8002fc1002d8b6856ea3f1dc264.

This is an evidence/routing diagnostic only. It consumes the exact successful
Agent-9 #1056 identity-bound morphology receipt and asks which of the six
qualitative public observables allowed by the live constrained routing document
are directly represented by candidate-bound measurements. It does not fit a
public image, introduce a numerical target, select a coefficient, or alter the
candidate velocity/basis/PDE state.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any, Iterable

TASK_ID = "CR003-ST052M-PUBLIC-OBSERVABLE-COVERAGE-129"
PREREG_ISSUE = 1086
SOURCE_PARENT_PR = 1078
SOURCE_PARENT_HEAD = "849713ba95cde8002fc1002d8b6856ea3f1dc264"

A9_SOURCE_PR = 1056
A9_SOURCE_HEAD = "751d954a33d420fb3c9b3381aef7d4ba274e422a"
A9_SOURCE_RUN = 35655112496
A9_ARTIFACT_ID = 10664106855
A9_ARTIFACT_SHA256 = "754d3109e1bfc3ff9e5855f73e40fb499ab6e8394522c4afc3190efb17b491a1"
EXPECTED_CANDIDATE_ID = "ST052-M-linear-temporal-child-v1"
EXPECTED_CANDIDATE_SHA256 = "be6eeb6d15aa2e146e87a6feeaae3627ac9ce3c4be6b3fc19b4cf8393bf5973c"
EXPECTED_VELOCITY_IDENTITY_SHA256 = "b0654bd8955644845fcb3789628299318b43bdb72bb467c284822235ffbeb902"
EXPECTED_MEASUREMENT_SHA256 = "7eb2ae35e12cd375eeb24f5e424855c888dd2a573ac8b85653440d6d38d4ed22"

PUBLIC_OBSERVABLES = (
    "vortex_swirl_presence",
    "inward_spiraling_trajectories",
    "axial_stretching_elongation",
    "shrinking_central_region_while_speed_increases",
    "spatial_variation_in_angular_rotation",
    "radius_dependent_circulation_speed",
)

PUBLIC_LABEL_TEXT = {
    "vortex_swirl_presence": "vortex / swirl presence",
    "inward_spiraling_trajectories": "inward-spiraling trajectories",
    "axial_stretching_elongation": "axial stretching / elongation",
    "shrinking_central_region_while_speed_increases": "shrinking central region while speed increases",
    "spatial_variation_in_angular_rotation": "spatial variation in angular rotation",
    "radius_dependent_circulation_speed": "radius-dependent circulation speed",
}

AXIAL_ASPECT_KEYS = {
    "axial_vorticity_aspect",
    "vorticity_aspect",
    "enstrophy_aspect",
    "full_aspect_ratio",
    "axial_to_transverse_aspect",
}

TRUTH = {
    "candidate_velocity_changed": False,
    "basis_dimension_changed": False,
    "coefficient_selected": False,
    "optimization_performed": False,
    "public_image_numeric_target_used": False,
    "pixel_similarity_objective_used": False,
    "renderer_or_camera_fit_used": False,
    "held_out_pde_residual_evaluated": False,
    "st006_comparison_performed": False,
    "candidate_mutation_authorized": False,
    "basis_growth_authorized": False,
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


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def _finite_positive(values: Iterable[Any]) -> bool:
    vals = [float(v) for v in values]
    return bool(vals) and all(math.isfinite(v) and v > 0.0 for v in vals)


def _walk_keys(value: Any) -> Iterable[tuple[str, Any]]:
    if isinstance(value, dict):
        for key, child in value.items():
            yield str(key), child
            yield from _walk_keys(child)
    elif isinstance(value, list):
        for child in value:
            yield from _walk_keys(child)


def _validate_identity(receipt: dict[str, Any]) -> dict[str, Any]:
    top = receipt.get("candidate_identity")
    morph = receipt.get("morphology_receipt")
    _require(isinstance(top, dict), "missing top-level candidate identity")
    _require(isinstance(morph, dict), "missing morphology_receipt")
    inner = morph.get("candidate_identity")
    _require(isinstance(inner, dict), "missing nested morphology candidate identity")

    for where, ident in (("top", top), ("morphology", inner)):
        _require(ident.get("candidate_id") == EXPECTED_CANDIDATE_ID, f"{where} candidate_id drift")
        _require(
            ident.get("candidate_sha256") == EXPECTED_CANDIDATE_SHA256,
            f"{where} candidate_sha256 drift",
        )
        _require(
            ident.get("velocity_identity_sha256") == EXPECTED_VELOCITY_IDENTITY_SHA256,
            f"{where} velocity identity drift",
        )

    _require(morph.get("candidate_identity_bound") is True, "morphology receipt is not identity-bound")
    _require(
        morph.get("measurement_sha256") == EXPECTED_MEASUREMENT_SHA256,
        "morphology measurement digest drift",
    )
    _require(
        morph.get("cylindrical_morphology_diagnostic_ready") is True,
        "morphology diagnostic not ready",
    )
    return morph


def audit_receipt(receipt: dict[str, Any]) -> dict[str, Any]:
    morph = _validate_identity(receipt)
    measurements = morph.get("measurements")
    _require(isinstance(measurements, dict), "missing morphology measurements")

    radial_variation = measurements.get("radial_variation")
    ring_rows = measurements.get("ring_rows")
    time_proxy = measurements.get("time_evolution_proxy")
    _require(isinstance(radial_variation, list) and radial_variation, "missing radial_variation rows")
    _require(isinstance(ring_rows, list) and ring_rows, "missing ring_rows")
    _require(isinstance(time_proxy, dict), "missing time_evolution_proxy")

    all_swirl = all(row.get("all_rings_swirl_nonzero") is True for row in radial_variation)
    positive_swirl_fraction = _finite_positive(row.get("swirl_nonzero_fraction", 0.0) for row in ring_rows)
    swirl_covered = bool(all_swirl and positive_swirl_fraction)

    all_inward = all(row.get("all_rings_inward") is True for row in radial_variation)
    positive_inward_swirl_fraction = _finite_positive(
        row.get("inward_swirl_fraction", 0.0) for row in ring_rows
    )
    inward_spiral_covered = bool(all_inward and all_swirl and positive_inward_swirl_fraction)

    shrink_speed_covered = bool(
        time_proxy.get("speed_weighted_radius_decreased") is True
        and time_proxy.get("peak_mean_speed_increased") is True
    )

    angular_variation_covered = _finite_positive(
        row.get("angular_rate_span_across_radii", 0.0) for row in radial_variation
    )
    circulation_variation_covered = _finite_positive(
        row.get("circulation_speed_span_across_radii", 0.0) for row in radial_variation
    )

    explicit_axial_aspect = [
        (key, value)
        for key, value in _walk_keys(measurements)
        if key in AXIAL_ASPECT_KEYS and isinstance(value, (int, float)) and math.isfinite(float(value))
    ]
    axial_status = (
        "measured_but_no_public_magnitude_target"
        if explicit_axial_aspect
        else "not_measured_by_current_receipt"
    )

    coverage = {
        "vortex_swirl_presence": {
            "public_label": PUBLIC_LABEL_TEXT["vortex_swirl_presence"],
            "status": "covered_by_candidate_bound_proxy" if swirl_covered else "not_covered_by_current_receipt",
            "proxy": "all sampled rings report nonzero swirl",
            "scope_limit": "presence only; no public magnitude target",
        },
        "inward_spiraling_trajectories": {
            "public_label": PUBLIC_LABEL_TEXT["inward_spiraling_trajectories"],
            "status": "covered_by_candidate_bound_proxy" if inward_spiral_covered else "not_covered_by_current_receipt",
            "proxy": "all sampled rings are inward while swirl is nonzero",
            "scope_limit": "renderer-independent ring proxy only; not a streamline proof",
        },
        "axial_stretching_elongation": {
            "public_label": PUBLIC_LABEL_TEXT["axial_stretching_elongation"],
            "status": axial_status,
            "proxy": None if not explicit_axial_aspect else "explicit candidate-bound axial/vorticity aspect observable present",
            "explicit_aspect_keys_found": [key for key, _ in explicit_axial_aspect],
            "scope_limit": "do not infer elongation from axial velocity signs or ring data",
        },
        "shrinking_central_region_while_speed_increases": {
            "public_label": PUBLIC_LABEL_TEXT["shrinking_central_region_while_speed_increases"],
            "status": "covered_by_candidate_bound_proxy" if shrink_speed_covered else "not_covered_by_current_receipt",
            "proxy": "speed-weighted radius decreases while peak mean speed increases from first to last registered time",
            "scope_limit": "trend sign only; no public rate/magnitude target or frame-to-time correspondence",
        },
        "spatial_variation_in_angular_rotation": {
            "public_label": PUBLIC_LABEL_TEXT["spatial_variation_in_angular_rotation"],
            "status": "covered_by_candidate_bound_proxy" if angular_variation_covered else "not_covered_by_current_receipt",
            "proxy": "strictly positive angular-rate span across sampled radii in every recorded row",
            "scope_limit": "variation presence only; no public magnitude target",
        },
        "radius_dependent_circulation_speed": {
            "public_label": PUBLIC_LABEL_TEXT["radius_dependent_circulation_speed"],
            "status": "covered_by_candidate_bound_proxy" if circulation_variation_covered else "not_covered_by_current_receipt",
            "proxy": "strictly positive circulation-speed span across sampled radii in every recorded row",
            "scope_limit": "radius dependence presence only; no public magnitude target",
        },
    }

    _require(tuple(coverage) == PUBLIC_OBSERVABLES, "public-observable ordering drifted")
    covered = [name for name, row in coverage.items() if row["status"] == "covered_by_candidate_bound_proxy"]
    unresolved = [name for name, row in coverage.items() if row["status"] != "covered_by_candidate_bound_proxy"]
    only_axial_gap = unresolved == ["axial_stretching_elongation"]

    routing = {
        "only_uncovered_public_label_is_axial_stretching": only_axial_gap,
        "next_existing_agent7_coordinate": "axial_vorticity_aspect_t050" if only_axial_gap else None,
        "next_existing_agent7_parent_pr": 1078 if only_axial_gap else None,
        "candidate_mutation_authorized": False,
        "basis_growth_authorized": False,
        "required_before_candidate_change": [
            "candidate-bound axial/vorticity aspect evidence producing a concrete qualitative discrepancy",
            "numerically resolved exact-head #1078 defect-coordinate Jacobian identifying an existing control or genuine unresponsive obstruction",
        ],
        "interpretation": (
            "If axial stretching remains the sole evidence gap, measure the already-defined axial_vorticity_aspect_t050 coordinate before adding any basis."
            if only_axial_gap
            else "Resolve every uncovered qualitative label before selecting a morphology coefficient or adding a basis."
        ),
    }

    return {
        "task_id": TASK_ID,
        "prereg_issue": PREREG_ISSUE,
        "source_parent": {"pr": SOURCE_PARENT_PR, "head": SOURCE_PARENT_HEAD},
        "a9_source": {
            "pr": A9_SOURCE_PR,
            "head": A9_SOURCE_HEAD,
            "successful_run": A9_SOURCE_RUN,
            "artifact_id": A9_ARTIFACT_ID,
            "artifact_sha256": A9_ARTIFACT_SHA256,
        },
        "candidate_identity": {
            "candidate_id": EXPECTED_CANDIDATE_ID,
            "candidate_sha256": EXPECTED_CANDIDATE_SHA256,
            "velocity_identity_sha256": EXPECTED_VELOCITY_IDENTITY_SHA256,
            "measurement_sha256": EXPECTED_MEASUREMENT_SHA256,
        },
        "coverage": coverage,
        "coverage_summary": {
            "public_observable_count": len(PUBLIC_OBSERVABLES),
            "covered_by_candidate_bound_proxy_count": len(covered),
            "covered": covered,
            "unresolved_or_unmeasured": unresolved,
            "visual_correspondence_established": False,
        },
        "routing": routing,
        "truth": dict(TRUTH),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--receipt", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    receipt = json.loads(args.receipt.read_text())
    report = audit_receipt(receipt)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report["coverage_summary"], sort_keys=True))
    print(json.dumps(report["routing"], sort_keys=True))


if __name__ == "__main__":
    main()
