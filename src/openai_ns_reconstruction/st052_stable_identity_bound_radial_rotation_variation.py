"""Stable-identity-bound radial angular-rotation variation diagnostic for ST052-M.

Preregistered in issue #1202.  This is a candidate-side morphology diagnostic only:
it evaluates the unchanged frozen ST052-M callable through the same exact bundle/load
path used by the admitted stable-identity axial diagnostic.  It does not add a basis,
select a coefficient, fit pixels, infer an OpenAI numerical target, or perform PDE
acceptance.
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

SCHEMA = "st052-stable-identity-bound-radial-rotation-variation/v1"
TASK_ID = "CR003-ST052M-STABLE-RADIAL-ROTATION-VARIATION-142"
PREREG_ISSUE = 1202
BASE_MAIN = "c340901736c3508e4e8dc25ffa01ba90d9145882"
CANDIDATE_ID = stable_identity.CANDIDATE_ID
EXPECTED_STABLE_CANDIDATE_IDENTITY = "2b2e743e72e31f69ad4232340f52d48b5759ac8aad1dfb735708112276c0b9c4"
EXPECTED_STABLE_VELOCITY_IDENTITY = "7738013a4fe4cacbfebfb4e8822bb4ee815b1d3bab744b90f9473ad9424f4f41"
STABLE_IDENTITY_MAIN_MERGE = "0713513575178f11b60a59fc7bcd1d1ecd27b1da"
STABLE_AXIAL_MAIN_MERGE = "b0093e04039946923de10e421e001cd363213dcc"

TIME = 0.50
TIMES = (0.25, 0.50, 0.75)
RING_RADII = (0.40, 0.70, 1.00, 1.30)
RING_ABS_Z = (0.30, 0.80)
AZIMUTH_COUNT = 16
CORE_RADII = tuple(float(v) for v in np.linspace(0.20, 1.60, 15))
SUPPORT_RADIUS = 2.0
NONZERO_VARIATION_FLOOR = 1.0e-10

PUBLIC_OBSERVABLE = {
    "id": "radius_dependent_circulation_speed",
    "kind": "qualitative_velocity_structure",
    "publisher": "OpenAI",
    "source_title": "On the Navier-Stokes Millennium Prize Problem",
    "source_url": "https://openai.com/index/navier-stokes-solution/",
    "public_observation": "The public visualization describes spatial/radial variation in angular rotation or circulation speed.",
    "numerical_target": None,
}

METRIC_PROVENANCE = {
    "coordinate_source_pr": 1078,
    "coordinate_source_head": "849713ba95cde8002fc1002d8b6856ea3f1dc264",
    "coordinate_source_path": "experiments/root_st052/agent7_st052m_defect_coordinate_controllability.py",
    "observable_source_pr": 1069,
    "observable_source_head": "a564fbc8bfa6a5c630b5839fa208e41cb7bc72d7",
    "observable_source_path": "experiments/root_st052/agent7_st052m_observable_space_five_channel_jacobian.py",
    "classification": "internal_protocol_reimplementation",
    "migration_scope": "angular_rotation_radial_variation_t050 and its frozen baseline-speed/ring sampling family only",
}

TRUTH_BOUNDARY = {
    "candidate_changed": False,
    "velocity_coefficients_changed": False,
    "basis_dimension_changed": False,
    "new_spatial_basis_added": False,
    "new_temporal_basis_added": False,
    "coefficient_selected": False,
    "pressure_changed": False,
    "forcing_changed": False,
    "held_out_pde_residual_evaluated": False,
    "public_source_numeric_targets_used": False,
    "renderer_or_camera_used": False,
    "pixel_loss_used": False,
    "stable_semantic_identity_bound": True,
    "radial_rotation_variation_measured": True,
    "direct_visualization_fingerprint_improvement": 0.0,
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


def _azimuths() -> np.ndarray:
    return np.arange(AZIMUTH_COUNT, dtype=float) * (2.0 * np.pi / AZIMUTH_COUNT)


def _ring_points(radius: float, z: float) -> np.ndarray:
    radius = float(radius)
    angles = _azimuths()
    return np.column_stack(
        (
            radius * np.cos(angles),
            radius * np.sin(angles),
            np.full(AZIMUTH_COUNT, float(z), dtype=float),
        )
    )


def _velocity(field: VelocityField, points: np.ndarray, time: float) -> np.ndarray:
    pts = np.asarray(points, dtype=float)
    values = np.asarray(
        field.velocity(pts[:, 0], pts[:, 1], pts[:, 2], float(time)),
        dtype=float,
    )
    if values.shape != (len(pts), 3) or not np.isfinite(values).all():
        raise ValueError("velocity must return finite shape-(n,3) data")
    return values


def _cylindrical_components(points: np.ndarray, velocity: np.ndarray) -> np.ndarray:
    pts = np.asarray(points, dtype=float)
    vel = np.asarray(velocity, dtype=float)
    if pts.ndim != 2 or pts.shape[1] != 3 or vel.shape != pts.shape:
        raise ValueError("points/velocity must both have shape (n,3)")
    radius = np.hypot(pts[:, 0], pts[:, 1])
    if np.any(radius <= np.finfo(float).tiny):
        raise ValueError("ring diagnostic requires positive radius")
    radial = (pts[:, 0] * vel[:, 0] + pts[:, 1] * vel[:, 1]) / radius
    swirl = (-pts[:, 1] * vel[:, 0] + pts[:, 0] * vel[:, 1]) / radius
    axial = vel[:, 2]
    out = np.column_stack((radial, swirl, axial))
    if not np.isfinite(out).all():
        raise ValueError("nonfinite cylindrical velocity")
    return out


def _baseline_speed_scale(field: VelocityField) -> float:
    """Reproduce the #1069 baseline RMS-speed sample family exactly."""
    samples: list[np.ndarray] = []
    for time in TIMES:
        for abs_z in RING_ABS_Z:
            for sign in (-1.0, 1.0):
                for radius in RING_RADII:
                    pts = _ring_points(radius, sign * abs_z)
                    samples.append(_velocity(field, pts, time))
        for radius in CORE_RADII:
            pts = _ring_points(radius, 0.0)
            samples.append(_velocity(field, pts, time))
    values = np.vstack(samples)
    scale = float(np.sqrt(np.mean(np.sum(values * values, axis=1))))
    if not np.isfinite(scale) or scale <= np.finfo(float).tiny:
        raise ValueError("invalid baseline RMS-speed scale")
    return scale


def rotation_variation_metrics(field: VelocityField) -> dict[str, Any]:
    """Evaluate the frozen #1078 angular-rotation-variation coordinate."""
    baseline_speed = _baseline_speed_scale(field)
    rows: list[dict[str, float]] = []
    omega_star: list[float] = []
    for abs_z in RING_ABS_Z:
        for sign in (-1.0, 1.0):
            z = float(sign * abs_z)
            for radius in RING_RADII:
                pts = _ring_points(radius, z)
                cyl = _cylindrical_components(pts, _velocity(field, pts, TIME))
                mean_swirl = float(np.mean(cyl[:, 1]))
                omega = float(SUPPORT_RADIUS * mean_swirl / (float(radius) * baseline_speed))
                rows.append(
                    {
                        "radius": float(radius),
                        "z": z,
                        "mean_signed_u_theta": mean_swirl,
                        "omega_star": omega,
                    }
                )
                omega_star.append(omega)

    values = np.asarray(omega_star, dtype=float)
    if values.size != 16 or not np.isfinite(values).all():
        raise ValueError("expected sixteen finite angular-speed samples")
    variation = float(np.std(values, ddof=0))
    mean_abs = float(np.mean(np.abs(values)))
    minimum = float(np.min(values))
    maximum = float(np.max(values))
    sample_range = float(maximum - minimum)
    return {
        "baseline_rms_speed": baseline_speed,
        "angular_rotation_radial_variation_t050": variation,
        "omega_star_mean_abs": mean_abs,
        "omega_star_min": minimum,
        "omega_star_max": maximum,
        "omega_star_range": sample_range,
        "nonzero_variation_floor": NONZERO_VARIATION_FLOOR,
        "nonzero_radius_dependent_rotation_proxy": bool(variation > NONZERO_VARIATION_FLOOR),
        "samples": rows,
    }


def _protocol() -> dict[str, Any]:
    return {
        "time": TIME,
        "times_for_baseline_speed": list(TIMES),
        "ring_radii": list(RING_RADII),
        "ring_abs_z": list(RING_ABS_Z),
        "azimuth_count": AZIMUTH_COUNT,
        "core_radii": list(CORE_RADII),
        "support_radius_scale": SUPPORT_RADIUS,
        "nonzero_variation_floor": NONZERO_VARIATION_FLOOR,
        "source_numeric_targets_used": False,
        "renderer_or_camera_used": False,
        "pixel_loss_used": False,
    }


def _decision(metrics: dict[str, Any]) -> dict[str, Any]:
    nonzero = bool(metrics["nonzero_radius_dependent_rotation_proxy"])
    return {
        "radius_invariant_rotation_deficit_established": bool(not nonzero),
        "new_basis_authorized": False,
        "coefficient_change_authorized": False,
        "candidate_mutation_authorized": False,
        "actual_velocity_changed": False,
        "basis_dimension_change": 0,
        "direct_visualization_fingerprint_improvement": 0.0,
        "next_route": (
            "flat-rotation trigger closed; any future change requires a separately bound magnitude/trajectory discrepancy"
            if nonzero
            else "route the flat-rotation deficit first through existing #1078 radial_shape/toroidal_swirl controls before any new basis"
        ),
    }


def _measurement_binding_sha256(identity: dict[str, Any], measurement_sha256: str) -> str:
    return _canonical_sha256(
        {
            "schema": SCHEMA,
            "stable_candidate_identity": identity,
            "protocol": _protocol(),
            "measurement_sha256": measurement_sha256,
            "metric_provenance": METRIC_PROVENANCE,
        }
    )


def validate_receipt(receipt: dict[str, Any]) -> None:
    if receipt.get("schema") != SCHEMA or receipt.get("task_id") != TASK_ID:
        raise ValueError("unexpected radial-rotation receipt schema/task")
    if receipt.get("preregister_issue") != PREREG_ISSUE or receipt.get("base_main") != BASE_MAIN:
        raise ValueError("preregistration/base drift")
    identity = receipt.get("stable_candidate_identity")
    materialization = receipt.get("materialization_join")
    measurement = receipt.get("measurement")
    truth = receipt.get("truth_boundary")
    if not all(isinstance(v, dict) for v in (identity, materialization, measurement, truth)):
        raise ValueError("malformed radial-rotation receipt")
    expected_identity = {
        "candidate_id": CANDIDATE_ID,
        "candidate_semantic_identity_sha256": EXPECTED_STABLE_CANDIDATE_IDENTITY,
        "velocity_semantic_identity_sha256": EXPECTED_STABLE_VELOCITY_IDENTITY,
    }
    if identity != expected_identity:
        raise ValueError("stable candidate/callable identity drift")
    if receipt.get("protocol") != _protocol():
        raise ValueError("frozen radial-rotation protocol drift")
    if receipt.get("public_observable") != PUBLIC_OBSERVABLE:
        raise ValueError("public qualitative observable drift")
    if receipt.get("metric_provenance") != METRIC_PROVENANCE:
        raise ValueError("internal metric provenance drift")
    if truth != TRUTH_BOUNDARY:
        raise ValueError("truth boundary drift")
    for key in ("stable_identity_receipt_sha256", "legacy_whole_candidate_identity_sha256", "materialization_evidence_sha256"):
        if not _is_sha256(materialization.get(key)):
            raise ValueError(f"invalid materialization join {key}")
    if materialization.get("included_in_measurement_binding") is not False:
        raise ValueError("materialization evidence leaked into semantic measurement binding")
    if materialization.get("stable_identity_main_merge") != STABLE_IDENTITY_MAIN_MERGE:
        raise ValueError("stable identity main ancestry drift")
    if materialization.get("stable_axial_main_merge") != STABLE_AXIAL_MAIN_MERGE:
        raise ValueError("stable axial main ancestry drift")

    metrics = measurement.get("metrics")
    if not isinstance(metrics, dict):
        raise ValueError("missing radial-rotation metrics")
    variation = float(metrics.get("angular_rotation_radial_variation_t050"))
    min_v = float(metrics.get("omega_star_min"))
    max_v = float(metrics.get("omega_star_max"))
    range_v = float(metrics.get("omega_star_range"))
    if not all(np.isfinite(v) for v in (variation, min_v, max_v, range_v)) or variation < 0.0:
        raise ValueError("nonfinite/negative radial-rotation metric")
    if not np.isclose(range_v, max_v - min_v, rtol=0.0, atol=2.0e-14):
        raise ValueError("angular-speed range arithmetic drift")
    if metrics.get("nonzero_radius_dependent_rotation_proxy") is not bool(variation > NONZERO_VARIATION_FLOOR):
        raise ValueError("nonzero-variation proxy drift")
    samples = metrics.get("samples")
    if not isinstance(samples, list) or len(samples) != 16:
        raise ValueError("radial-rotation sample count drift")

    if receipt.get("decision") != _decision(metrics):
        raise ValueError("routing decision drift")
    measurement_sha = _canonical_sha256(measurement)
    if receipt.get("measurement_sha256") != measurement_sha:
        raise ValueError("measurement checksum mismatch")
    if receipt.get("candidate_measurement_binding_sha256") != _measurement_binding_sha256(identity, measurement_sha):
        raise ValueError("stable candidate/measurement binding mismatch")
    if receipt.get("receipt_sha256") is not None:
        unsigned = dict(receipt)
        recorded = unsigned.pop("receipt_sha256")
        if recorded != _canonical_sha256(unsigned):
            raise ValueError("receipt checksum mismatch")


def execute(*, constrained_root: str | Path, source_root: str | Path, bundle_dir: str | Path) -> dict[str, Any]:
    """Measure candidate-side radial rotation variation on the admitted stable ST052 identity."""
    constrained_root = Path(constrained_root).resolve()
    source_root = Path(source_root).resolve()
    bundle_dir = Path(bundle_dir).resolve()

    stable_receipt = stable_identity.execute(constrained_root=constrained_root, bundle_dir=bundle_dir)
    stable_identity.validate_receipt(stable_receipt)
    stable = stable_receipt["stable_identity"]
    candidate_semantic_id = stable["candidate_semantic_identity_sha256"]
    velocity_semantic_id = stable["velocity_semantic_identity_sha256"]
    if candidate_semantic_id != EXPECTED_STABLE_CANDIDATE_IDENTITY:
        raise ValueError("fresh stable candidate identity does not match admitted ST052 identity")
    if velocity_semantic_id != EXPECTED_STABLE_VELOCITY_IDENTITY:
        raise ValueError("fresh stable callable identity does not match admitted ST052 identity")

    identity_bridge._install_constrained_package_path(constrained_root)
    whole = importlib.import_module("openai_ns_reconstruction.st052_linear_temporal_capsule")
    candidate = whole.load_bundle_runtime(bundle_dir, exact_source_root=source_root)
    legacy_id = str(candidate.identity_sha256)
    materialization_evidence = stable_receipt["materialization_evidence"]
    if legacy_id != materialization_evidence["legacy_whole_candidate_identity_sha256"]:
        raise ValueError("loaded callable materialization disagrees with stable-identity-verified bundle")

    identity = {
        "candidate_id": CANDIDATE_ID,
        "candidate_semantic_identity_sha256": candidate_semantic_id,
        "velocity_semantic_identity_sha256": velocity_semantic_id,
    }
    measurement = {
        "time": TIME,
        "metrics": rotation_variation_metrics(candidate),
        "proxy_interpretation": (
            "nonzero variation means only that the frozen candidate is not radius/spatially rigid under this autonomous ring diagnostic; it is not a numerical OpenAI match or adequacy threshold"
        ),
    }
    measurement_sha = _canonical_sha256(measurement)
    receipt: dict[str, Any] = {
        "schema": SCHEMA,
        "task_id": TASK_ID,
        "preregister_issue": PREREG_ISSUE,
        "base_main": BASE_MAIN,
        "stable_candidate_identity": identity,
        "materialization_join": {
            "stable_identity_receipt_sha256": stable_receipt["receipt_sha256"],
            "legacy_whole_candidate_identity_sha256": legacy_id,
            "materialization_evidence_sha256": _canonical_sha256(materialization_evidence),
            "stable_identity_main_merge": STABLE_IDENTITY_MAIN_MERGE,
            "stable_axial_main_merge": STABLE_AXIAL_MAIN_MERGE,
            "included_in_measurement_binding": False,
        },
        "protocol": _protocol(),
        "public_observable": PUBLIC_OBSERVABLE,
        "metric_provenance": METRIC_PROVENANCE,
        "measurement": measurement,
        "measurement_sha256": measurement_sha,
        "candidate_measurement_binding_sha256": _measurement_binding_sha256(identity, measurement_sha),
        "decision": _decision(measurement["metrics"]),
        "truth_boundary": TRUTH_BOUNDARY,
    }
    receipt["receipt_sha256"] = _canonical_sha256(receipt)
    validate_receipt(receipt)
    return receipt


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--constrained-root", required=True)
    parser.add_argument("--source-root", required=True)
    parser.add_argument("--bundle-dir", required=True)
    parser.add_argument("--report", required=True)
    args = parser.parse_args(argv)
    receipt = execute(
        constrained_root=args.constrained_root,
        source_root=args.source_root,
        bundle_dir=args.bundle_dir,
    )
    path = Path(args.report)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    metrics = receipt["measurement"]["metrics"]
    print("angular_rotation_radial_variation_t050=", metrics["angular_rotation_radial_variation_t050"])
    print("nonzero_radius_dependent_rotation_proxy=", metrics["nonzero_radius_dependent_rotation_proxy"])
    print("next_route=", receipt["decision"]["next_route"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
