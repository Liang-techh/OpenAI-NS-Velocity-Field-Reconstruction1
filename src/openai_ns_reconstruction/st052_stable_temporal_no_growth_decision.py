"""Fail-closed Agent-7 basis-growth routing from the admitted ST052 core-trend receipt.

This module does not measure or modify the velocity field.  It authenticates the exact
successful #1222 morphology receipt and turns only that already-admitted evidence into
one bounded decision: the coarse public qualitative ingredient "central region shrinks
while speed increases" does not, by itself, justify adding a new temporal basis.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

SCHEMA = "st052-stable-temporal-no-growth-decision/v1"
TASK_ID = "CR003-ST052M-STABLE-TEMPORAL-NO-GROWTH-145"
ISSUE = 1231
BASE_MAIN = "08fb5f27adfb6811474f3f82e2ec4766196ac07c"
UPSTREAM_PR = 1222
UPSTREAM_RUN = 35745851137
UPSTREAM_TASK = "CR003-ST052M-STABLE-CORE-CONTRACTION-SPEED-144"
UPSTREAM_MEASUREMENT_SHA256 = "f8595de192fe1e54c56d92269a35ebfb874ee2658726a98071e50ad954b32c7c"
CANDIDATE_SEMANTIC_SHA256 = "2b2e743e72e31f69ad4232340f52d48b5759ac8aad1dfb735708112276c0b9c4"
VELOCITY_SEMANTIC_SHA256 = "7738013a4fe4cacbfebfb4e8822bb4ee815b1d3bab744b90f9473ad9424f4f41"
EXPECTED_WIDTHS = (0.9805978797661573, 0.9298981687005765, 0.8816987778759415)
EXPECTED_SPEEDS = (0.2750591930576447, 0.304524609725179, 0.33361562589768623)
EXPECTED_WIDTH_FRACTION = -0.10085592058775428
EXPECTED_SPEED_FRACTION = 0.21288665973716334
ATOL = 2e-15


def _close(a: float, b: float) -> bool:
    return math.isfinite(a) and math.isfinite(b) and abs(a - b) <= ATOL


def _require_close_sequence(actual: Any, expected: tuple[float, ...], name: str) -> list[float]:
    if not isinstance(actual, list) or len(actual) != len(expected):
        raise ValueError(f"{name} shape drift")
    values = [float(v) for v in actual]
    if not all(_close(v, e) for v, e in zip(values, expected, strict=True)):
        raise ValueError(f"{name} numerical drift")
    return values


def validate_upstream_receipt(receipt: dict[str, Any]) -> dict[str, Any]:
    """Authenticate the exact admitted #1222 receipt needed for this decision."""
    if receipt.get("task_id") != UPSTREAM_TASK:
        raise ValueError("unexpected upstream task")
    if receipt.get("issue") != 1220:
        raise ValueError("unexpected upstream issue")
    if receipt.get("measurement_sha256") != UPSTREAM_MEASUREMENT_SHA256:
        raise ValueError("upstream measurement identity drift")

    identity = receipt.get("stable_candidate_identity")
    if not isinstance(identity, dict):
        raise ValueError("missing stable candidate identity")
    if identity.get("candidate_semantic_identity_sha256") != CANDIDATE_SEMANTIC_SHA256:
        raise ValueError("candidate semantic identity drift")
    if identity.get("velocity_semantic_identity_sha256") != VELOCITY_SEMANTIC_SHA256:
        raise ValueError("velocity semantic identity drift")

    observable = receipt.get("public_observable")
    if not isinstance(observable, dict) or observable.get("id") != "shrinking_central_region_while_speed_increases":
        raise ValueError("public observable drift")
    if observable.get("numerical_target") is not None:
        raise ValueError("invented public numerical target")

    measurement = receipt.get("measurement")
    if not isinstance(measurement, dict):
        raise ValueError("missing upstream measurement")
    derived = measurement.get("derived")
    routing = measurement.get("routing")
    if not isinstance(derived, dict) or not isinstance(routing, dict):
        raise ValueError("malformed upstream measurement")

    widths = _require_close_sequence(
        derived.get("speed_weighted_rms_radius_by_time"), EXPECTED_WIDTHS, "core widths"
    )
    speeds = _require_close_sequence(
        derived.get("peak_ring_mean_speed_by_time"), EXPECTED_SPEEDS, "peak speeds"
    )
    if not _close(float(derived.get("endpoint_width_fraction")), EXPECTED_WIDTH_FRACTION):
        raise ValueError("endpoint width fraction drift")
    if not _close(float(derived.get("endpoint_peak_speed_fraction")), EXPECTED_SPEED_FRACTION):
        raise ValueError("endpoint peak-speed fraction drift")

    if not (widths[0] > widths[1] > widths[2]):
        raise ValueError("admitted core width is no longer monotonically decreasing")
    if not (speeds[0] < speeds[1] < speeds[2]):
        raise ValueError("admitted peak speed is no longer monotonically increasing")
    for key in (
        "endpoint_shrink_and_speedup_proxy",
        "three_time_monotone_proxy",
        "three_time_width_monotone_decrease",
        "three_time_peak_speed_monotone_increase",
    ):
        if derived.get(key) is not True:
            raise ValueError(f"upstream qualitative trend flag drift: {key}")
    if routing.get("sixth_basis_authorized") is not False or routing.get("candidate_mutation_authorized") is not False:
        raise ValueError("upstream routing improperly authorizes mutation/basis growth")

    truth = receipt.get("truth_boundary")
    if not isinstance(truth, dict):
        raise ValueError("missing upstream truth boundary")
    required_false = (
        "candidate_changed",
        "basis_dimension_changed",
        "coefficient_selected",
        "visualization_ready",
        "visual_correspondence_verified",
        "source_correspondence_verified",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    )
    if any(truth.get(key) is not False for key in required_false):
        raise ValueError("upstream truth promotion detected")
    if float(truth.get("visualization_fingerprint_direct_improvement")) != 0.0:
        raise ValueError("upstream direct visualization-improvement boundary drift")

    return {
        "widths": widths,
        "peak_speeds": speeds,
        "endpoint_width_fraction": EXPECTED_WIDTH_FRACTION,
        "endpoint_peak_speed_fraction": EXPECTED_SPEED_FRACTION,
    }


def decide(receipt: dict[str, Any]) -> dict[str, Any]:
    admitted = validate_upstream_receipt(receipt)
    return {
        "schema": SCHEMA,
        "task_id": TASK_ID,
        "issue": ISSUE,
        "base_main": BASE_MAIN,
        "upstream": {
            "pr": UPSTREAM_PR,
            "successful_exact_head_run": UPSTREAM_RUN,
            "measurement_sha256": UPSTREAM_MEASUREMENT_SHA256,
            "candidate_semantic_identity_sha256": CANDIDATE_SEMANTIC_SHA256,
            "velocity_semantic_identity_sha256": VELOCITY_SEMANTIC_SHA256,
        },
        "authenticated_evidence": admitted,
        "decision": {
            "coarse_missing_shrink_speedup_trigger_closed": True,
            "new_temporal_basis_authorized_from_this_observable": False,
            "sixth_basis_authorized": False,
            "coefficient_selection_authorized": False,
            "candidate_mutation_authorized": False,
            "existing_offgrid_temporal_capacity_pr": 1178,
            "offgrid_capacity_promotion_authorized": False,
            "next_step": (
                "retain the existing temporal basis for this observable; reconsider #1178 only after a separately bound temporal magnitude/profile discrepancy"
            ),
        },
        "interpretation": (
            "The stable ST052 callable already exhibits monotone candidate-side core shrink and peak-speed increase at t=.25/.50/.75. "
            "This closes only the coarse absence trigger; no OpenAI numerical time-law or magnitude target is claimed."
        ),
        "truth_boundary": {
            "velocity_changed": False,
            "basis_dimension_changed": False,
            "coefficient_selected": False,
            "pressure_changed": False,
            "forcing_changed": False,
            "scientific_threshold_changed": False,
            "visualization_fingerprint_direct_improvement": 0.0,
            "visualization_ready": False,
            "visual_correspondence_verified": False,
            "source_correspondence_verified": False,
            "pde_validated": False,
            "paper_exact": False,
            "openai_field_identified": False,
            "blowup_proved": False,
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--receipt", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    receipt = json.loads(args.receipt.read_text())
    report = decide(receipt)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report["decision"], sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
