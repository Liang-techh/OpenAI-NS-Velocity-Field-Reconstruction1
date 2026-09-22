"""Stable-identity-bound core contraction / speed-up diagnostic for frozen ST052-M.

This Agent-7 increment measures only the unchanged callable velocity.  It reuses the
already-merged renderer-independent midplane core-speed proxy and binds the three-time
trend to the admitted stable semantic candidate/callable identities.  It does not add
a basis, choose a coefficient, fit pixels, infer a public numerical target, or perform
Navier--Stokes acceptance.
"""
from __future__ import annotations

import argparse
import importlib
import json
from pathlib import Path
from typing import Any, Protocol

import numpy as np

from . import st052_identity_bound_morphology_execution as identity_bridge
from . import st052_stable_semantic_identity as stable_identity
from . import st054_cylindrical_morphology as morphology

SCHEMA = "st052-stable-identity-bound-core-contraction-speed/v1"
TASK_ID = "CR003-ST052M-STABLE-CORE-CONTRACTION-SPEED-144"
CANDIDATE_ID = stable_identity.CANDIDATE_ID
ISSUE = 1220
BASE_MAIN = "c340901736c3508e4e8dc25ffa01ba90d9145882"
STABLE_IDENTITY_MAIN_MERGE = "0713513575178f11b60a59fc7bcd1d1ecd27b1da"
EXPECTED_STABLE_CANDIDATE_IDENTITY = "2b2e743e72e31f69ad4232340f52d48b5759ac8aad1dfb735708112276c0b9c4"
EXPECTED_STABLE_VELOCITY_IDENTITY = "7738013a4fe4cacbfebfb4e8822bb4ee815b1d3bab744b90f9473ad9424f4f41"
CONSTRAINED_CALLABLE_IMPLEMENTATION = "6a293b3c870d70b8d8aece1b12d68685c9b8b010"
EXACT_SOURCE_RUNTIME = "b3b8bfdbe1077f9ec967d158602951997d81e17d"
MORPHOLOGY_SOURCE_BLOB = "0ab62bdbabd753ec3d180d5c2bab00e5802da2a2"
TIMES = (0.25, 0.50, 0.75)
CORE_RADII = tuple(float(v) for v in np.linspace(0.2, 1.6, 15))
AZIMUTH_COUNT = 16

PUBLIC_OBSERVABLE = {
    "id": "shrinking_central_region_while_speed_increases",
    "kind": "qualitative_time_evolution",
    "publisher": "OpenAI",
    "public_observation": "The displayed central region shrinks while speed increases.",
    "numerical_target": None,
    "source_contract_pr": 1219,
    "source_contract_admitted_by_this_increment": False,
}

METRIC_PROVENANCE = {
    "implementation": "openai_ns_reconstruction.st054_cylindrical_morphology:_core_speed_proxy",
    "implementation_blob_sha": MORPHOLOGY_SOURCE_BLOB,
    "classification": "reuse_existing_renderer_independent_candidate_metric",
    "migration_scope": "midplane ring-mean speed, speed-weighted RMS radius, and peak ring-mean speed only",
    "retuned": False,
}

TRUTH_BOUNDARY = {
    "candidate_changed": False,
    "basis_dimension_changed": False,
    "coefficient_selected": False,
    "velocity_coefficients_changed": False,
    "pressure_changed": False,
    "forcing_changed": False,
    "scientific_threshold_changed": False,
    "public_source_numeric_targets_used": False,
    "renderer_or_camera_used": False,
    "pixel_loss_used": False,
    "stable_semantic_identity_bound": True,
    "core_contraction_speed_trend_measured": True,
    "visualization_fingerprint_direct_improvement": 0.0,
    "visualization_ready": False,
    "visual_correspondence_verified": False,
    "source_correspondence_verified": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


class VelocityField(Protocol):
    def velocity(self, x: Any, y: Any, z: Any, t: float) -> np.ndarray:
        """Return Cartesian velocity with final component axis of length three."""


def _canonical_sha256(value: Any) -> str:
    return stable_identity.canonical_sha256(value)


def _is_sha256(value: Any) -> bool:
    return isinstance(value, str) and len(value) == 64 and all(c in "0123456789abcdef" for c in value)


def protocol() -> dict[str, Any]:
    return {
        "times": list(TIMES),
        "z": 0.0,
        "core_radii": list(CORE_RADII),
        "azimuth_count": AZIMUTH_COUNT,
        "ring_statistic": "mean Cartesian speed magnitude over equally spaced azimuths",
        "core_width_proxy": "sqrt(sum(mean_speed_i*r_i^2)/sum(mean_speed_i))",
        "speed_proxy": "max ring mean speed",
        "endpoint_width_fraction": "(R_w(t=.75)-R_w(t=.25))/R_w(t=.25)",
        "endpoint_peak_speed_fraction": "(S_peak(t=.75)-S_peak(t=.25))/S_peak(t=.25)",
        "source_numeric_targets_used": False,
        "renderer_or_camera_used": False,
        "pixel_loss_used": False,
    }


def _routing(endpoint_proxy: bool, monotone_proxy: bool, width_decreased: bool, speed_increased: bool) -> dict[str, Any]:
    if endpoint_proxy and monotone_proxy:
        decision = "close_coarse_absence_trigger_no_temporal_basis_growth_from_this_presence_proxy"
        next_step = "retain existing basis unless a separately bound magnitude/profile discrepancy appears"
    elif endpoint_proxy:
        decision = "route_temporal_shape_defect_to_existing_temporal_curvature_and_offgrid_time_skewed_capacity"
        next_step = "use existing Agent-7 temporal controls/capacity before any new spatial basis"
    else:
        decision = "route_endpoint_sign_defect_to_existing_five_control_temporal_sensitivity"
        next_step = "a sixth basis remains unauthorized until unresponsive or ill-conditioned existing-control evidence exists"
    failures: list[str] = []
    if not width_decreased:
        failures.append("speed_weighted_core_radius_did_not_decrease")
    if not speed_increased:
        failures.append("peak_ring_mean_speed_did_not_increase")
    if endpoint_proxy and not monotone_proxy:
        failures.append("three_time_midpoint_monotonicity_failed")
    return {
        "decision": decision,
        "next_step": next_step,
        "candidate_mutation_authorized": False,
        "sixth_basis_authorized": False,
        "coefficient_selection_authorized": False,
        "failures": failures,
    }


def measure_field(field: VelocityField) -> dict[str, Any]:
    """Measure the preregistered three-time core-width / speed trend."""
    core = [
        morphology._core_speed_proxy(  # noqa: SLF001 - deliberate exact metric reuse
            field,
            time=t,
            radii=CORE_RADII,
            azimuth_count=AZIMUTH_COUNT,
        )
        for t in TIMES
    ]
    widths = [item["speed_weighted_rms_radius"] for item in core]
    speeds = [float(item["peak_mean_speed"]) for item in core]
    if any(value is None for value in widths):
        raise ValueError("core speed proxy produced an undefined weighted radius")
    widths_f = [float(value) for value in widths]
    if not all(np.isfinite(v) and v > 0.0 for v in widths_f + speeds):
        raise ValueError("core trend requires finite positive width and peak-speed proxies")

    width_decreased = bool(widths_f[-1] < widths_f[0])
    speed_increased = bool(speeds[-1] > speeds[0])
    endpoint_proxy = bool(width_decreased and speed_increased)
    monotone_width = bool(widths_f[0] > widths_f[1] > widths_f[2])
    monotone_speed = bool(speeds[0] < speeds[1] < speeds[2])
    monotone_proxy = bool(monotone_width and monotone_speed)

    derived = {
        "speed_weighted_rms_radius_by_time": widths_f,
        "peak_ring_mean_speed_by_time": speeds,
        "endpoint_width_fraction": float((widths_f[-1] - widths_f[0]) / widths_f[0]),
        "endpoint_peak_speed_fraction": float((speeds[-1] - speeds[0]) / speeds[0]),
        "endpoint_width_decreased": width_decreased,
        "endpoint_peak_speed_increased": speed_increased,
        "endpoint_shrink_and_speedup_proxy": endpoint_proxy,
        "three_time_width_monotone_decrease": monotone_width,
        "three_time_peak_speed_monotone_increase": monotone_speed,
        "three_time_monotone_proxy": monotone_proxy,
    }
    return {
        "core_speed_proxy": core,
        "derived": derived,
        "routing": _routing(endpoint_proxy, monotone_proxy, width_decreased, speed_increased),
        "interpretation": (
            "Candidate-side renderer-independent temporal morphology only. Strict ordering has no OpenAI numerical magnitude target and is not a source-correspondence or PDE-validity verdict."
        ),
    }


def _measurement_binding_sha256(identity: dict[str, Any], measurement_sha256: str) -> str:
    return _canonical_sha256(
        {
            "schema": SCHEMA,
            "stable_candidate_identity": identity,
            "protocol": protocol(),
            "metric_provenance": METRIC_PROVENANCE,
            "measurement_sha256": measurement_sha256,
        }
    )


def validate_measurement(measurement: dict[str, Any]) -> None:
    core = measurement.get("core_speed_proxy")
    derived = measurement.get("derived")
    routing = measurement.get("routing")
    if not isinstance(core, list) or len(core) != 3 or not isinstance(derived, dict) or not isinstance(routing, dict):
        raise ValueError("malformed core contraction/speed measurement")
    if [float(item.get("time")) for item in core] != list(TIMES):
        raise ValueError("core trend time grid drift")
    widths = [float(v) for v in derived.get("speed_weighted_rms_radius_by_time", [])]
    speeds = [float(v) for v in derived.get("peak_ring_mean_speed_by_time", [])]
    if len(widths) != 3 or len(speeds) != 3 or not all(np.isfinite(v) and v > 0.0 for v in widths + speeds):
        raise ValueError("invalid core trend values")
    expected_width_frac = (widths[-1] - widths[0]) / widths[0]
    expected_speed_frac = (speeds[-1] - speeds[0]) / speeds[0]
    if not np.isclose(float(derived.get("endpoint_width_fraction")), expected_width_frac, rtol=0.0, atol=2e-15):
        raise ValueError("endpoint width-fraction arithmetic drift")
    if not np.isclose(float(derived.get("endpoint_peak_speed_fraction")), expected_speed_frac, rtol=0.0, atol=2e-15):
        raise ValueError("endpoint speed-fraction arithmetic drift")
    width_decreased = widths[-1] < widths[0]
    speed_increased = speeds[-1] > speeds[0]
    endpoint = width_decreased and speed_increased
    monotone_width = widths[0] > widths[1] > widths[2]
    monotone_speed = speeds[0] < speeds[1] < speeds[2]
    monotone = monotone_width and monotone_speed
    expected_flags = {
        "endpoint_width_decreased": width_decreased,
        "endpoint_peak_speed_increased": speed_increased,
        "endpoint_shrink_and_speedup_proxy": endpoint,
        "three_time_width_monotone_decrease": monotone_width,
        "three_time_peak_speed_monotone_increase": monotone_speed,
        "three_time_monotone_proxy": monotone,
    }
    for key, value in expected_flags.items():
        if derived.get(key) is not bool(value):
            raise ValueError(f"core trend flag drift: {key}")
    expected_routing = _routing(endpoint, monotone, width_decreased, speed_increased)
    if routing != expected_routing:
        raise ValueError("core trend routing drift")


def validate_receipt(receipt: dict[str, Any]) -> None:
    if receipt.get("schema") != SCHEMA or receipt.get("task_id") != TASK_ID or receipt.get("issue") != ISSUE:
        raise ValueError("unexpected core contraction/speed receipt schema/task")
    identity = receipt.get("stable_candidate_identity")
    materialization = receipt.get("materialization_join")
    measurement = receipt.get("measurement")
    truth = receipt.get("truth_boundary")
    if not all(isinstance(v, dict) for v in (identity, materialization, measurement, truth)):
        raise ValueError("malformed core contraction/speed receipt")
    expected_identity = {
        "candidate_id": CANDIDATE_ID,
        "candidate_semantic_identity_sha256": EXPECTED_STABLE_CANDIDATE_IDENTITY,
        "velocity_semantic_identity_sha256": EXPECTED_STABLE_VELOCITY_IDENTITY,
    }
    if identity != expected_identity:
        raise ValueError("stable candidate/callable identity drift")
    if receipt.get("protocol") != protocol():
        raise ValueError("frozen core trend protocol drift")
    if receipt.get("metric_provenance") != METRIC_PROVENANCE:
        raise ValueError("core metric provenance drift")
    if receipt.get("public_observable") != PUBLIC_OBSERVABLE:
        raise ValueError("public qualitative observable drift")
    if truth != TRUTH_BOUNDARY:
        raise ValueError("core trend truth boundary drift")
    if materialization.get("stable_identity_main_merge") != STABLE_IDENTITY_MAIN_MERGE:
        raise ValueError("stable identity ancestry drift")
    for key in ("stable_identity_receipt_sha256", "legacy_whole_candidate_identity_sha256", "materialization_evidence_sha256"):
        if not _is_sha256(materialization.get(key)):
            raise ValueError(f"invalid materialization evidence: {key}")
    if materialization.get("included_in_measurement_binding") is not False:
        raise ValueError("materialization evidence leaked into stable semantic measurement binding")
    validate_measurement(measurement)
    measurement_sha = _canonical_sha256(measurement)
    if receipt.get("measurement_sha256") != measurement_sha:
        raise ValueError("core trend measurement checksum mismatch")
    if receipt.get("candidate_measurement_binding_sha256") != _measurement_binding_sha256(identity, measurement_sha):
        raise ValueError("stable candidate/core-trend binding mismatch")
    if receipt.get("receipt_sha256") is not None:
        unsigned = dict(receipt)
        recorded = unsigned.pop("receipt_sha256")
        if recorded != _canonical_sha256(unsigned):
            raise ValueError("core trend receipt checksum mismatch")


def execute(*, constrained_root: str | Path, source_root: str | Path, bundle_dir: str | Path) -> dict[str, Any]:
    constrained_root = Path(constrained_root).resolve()
    source_root = Path(source_root).resolve()
    bundle_dir = Path(bundle_dir).resolve()

    stable_receipt = stable_identity.execute(constrained_root=constrained_root, bundle_dir=bundle_dir)
    stable_identity.validate_receipt(stable_receipt)
    stable = stable_receipt["stable_identity"]
    candidate_semantic_id = stable["candidate_semantic_identity_sha256"]
    velocity_semantic_id = stable["velocity_semantic_identity_sha256"]
    if candidate_semantic_id != EXPECTED_STABLE_CANDIDATE_IDENTITY:
        raise ValueError("fresh stable candidate semantic identity does not match admitted identity")
    if velocity_semantic_id != EXPECTED_STABLE_VELOCITY_IDENTITY:
        raise ValueError("fresh stable callable semantic identity does not match admitted identity")

    identity_bridge._install_constrained_package_path(constrained_root)
    whole = importlib.import_module("openai_ns_reconstruction.st052_linear_temporal_capsule")
    candidate = whole.load_bundle_runtime(bundle_dir, exact_source_root=source_root)
    legacy_id = str(candidate.identity_sha256)
    stable_legacy_id = stable_receipt["materialization_evidence"]["legacy_whole_candidate_identity_sha256"]
    if legacy_id != stable_legacy_id:
        raise ValueError("loaded callable materialization disagrees with stable-identity-verified bundle")

    identity = {
        "candidate_id": CANDIDATE_ID,
        "candidate_semantic_identity_sha256": candidate_semantic_id,
        "velocity_semantic_identity_sha256": velocity_semantic_id,
    }
    measurement = measure_field(candidate)
    validate_measurement(measurement)
    measurement_sha = _canonical_sha256(measurement)
    receipt: dict[str, Any] = {
        "schema": SCHEMA,
        "task_id": TASK_ID,
        "issue": ISSUE,
        "base_main": BASE_MAIN,
        "stable_candidate_identity": identity,
        "materialization_join": {
            "stable_identity_main_merge": STABLE_IDENTITY_MAIN_MERGE,
            "stable_identity_receipt_sha256": stable_receipt["receipt_sha256"],
            "legacy_whole_candidate_identity_sha256": legacy_id,
            "materialization_evidence_sha256": _canonical_sha256(stable_receipt["materialization_evidence"]),
            "included_in_measurement_binding": False,
        },
        "protocol": protocol(),
        "metric_provenance": METRIC_PROVENANCE,
        "public_observable": PUBLIC_OBSERVABLE,
        "measurement": measurement,
        "measurement_sha256": measurement_sha,
        "candidate_measurement_binding_sha256": _measurement_binding_sha256(identity, measurement_sha),
        "direct_contribution": (
            "measures whether the unchanged stable ST052 callable already has endpoint and three-time core shrink+speed-up structure before any temporal basis growth"
        ),
        "truth_boundary": dict(TRUTH_BOUNDARY),
    }
    receipt["receipt_sha256"] = _canonical_sha256(receipt)
    validate_receipt(receipt)
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--constrained-root", type=Path, required=True)
    parser.add_argument("--source-root", type=Path, required=True)
    parser.add_argument("--bundle-dir", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    receipt = execute(
        constrained_root=args.constrained_root,
        source_root=args.source_root,
        bundle_dir=args.bundle_dir,
    )
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(receipt, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    d = receipt["measurement"]["derived"]
    print(json.dumps({
        "candidate_semantic_identity_sha256": receipt["stable_candidate_identity"]["candidate_semantic_identity_sha256"],
        "velocity_semantic_identity_sha256": receipt["stable_candidate_identity"]["velocity_semantic_identity_sha256"],
        "endpoint_width_fraction": d["endpoint_width_fraction"],
        "endpoint_peak_speed_fraction": d["endpoint_peak_speed_fraction"],
        "endpoint_shrink_and_speedup_proxy": d["endpoint_shrink_and_speedup_proxy"],
        "three_time_monotone_proxy": d["three_time_monotone_proxy"],
        "routing": receipt["measurement"]["routing"]["decision"],
        "measurement_sha256": receipt["measurement_sha256"],
        "pde_validated": receipt["truth_boundary"]["pde_validated"],
    }, sort_keys=True))


if __name__ == "__main__":
    main()
