"""Stable-identity-bound fixed-seed streamline morphology for frozen ST052-M.

Preregistered in issue #1211.  This is a candidate-side morphology diagnostic only:
it evaluates the unchanged admitted ST052-M callable through the same authenticated
bundle/load path as the merged stable-identity diagnostics.  It does not add a basis,
select a coefficient, infer an OpenAI numerical target, or perform PDE acceptance.
"""
from __future__ import annotations

import argparse
import importlib
import json
from pathlib import Path
from typing import Any, Protocol

import numpy as np
from scipy.interpolate import RegularGridInterpolator

from . import st052_identity_bound_morphology_execution as identity_bridge
from . import st052_stable_semantic_identity as stable_identity

SCHEMA = "st052-stable-identity-bound-streamline-trajectory/v1"
TASK_ID = "CR003-ST052M-STABLE-STREAMLINE-TRAJECTORY-143"
PREREG_ISSUE = 1211
BASE_MAIN = "c340901736c3508e4e8dc25ffa01ba90d9145882"
CANDIDATE_ID = stable_identity.CANDIDATE_ID
EXPECTED_STABLE_CANDIDATE_IDENTITY = "2b2e743e72e31f69ad4232340f52d48b5759ac8aad1dfb735708112276c0b9c4"
EXPECTED_STABLE_VELOCITY_IDENTITY = "7738013a4fe4cacbfebfb4e8822bb4ee815b1d3bab744b90f9473ad9424f4f41"
STABLE_IDENTITY_MAIN_MERGE = "0713513575178f11b60a59fc7bcd1d1ecd27b1da"
STABLE_AXIAL_MAIN_MERGE = "b0093e04039946923de10e421e001cd363213dcc"

TIME = 0.50
GRID_RESOLUTION = 33
BOX = (-2.0, 2.0)
SEED_RADII = (0.6, 0.9, 1.2)
SEED_Z = (-0.3, 0.3)
SEED_AZIMUTHS = 8
STREAM_ARCLENGTH = 2.25
STREAM_OUTPUTS_PER_SIDE = 91
SPEED_FLOOR = 1.0e-8
ACTIVE_BOX_MAX = 1.995
NONDEGENERACY_FLOOR = 1.0e-8

PUBLIC_OBSERVABLES = (
    {
        "id": "inward_spiraling_trajectories",
        "kind": "qualitative_trajectory",
        "publisher": "OpenAI",
        "public_observation": "Displayed trajectories spiral while moving inward toward the central vortex region.",
        "numerical_target": None,
    },
    {
        "id": "axial_stretching_trajectories",
        "kind": "qualitative_trajectory",
        "publisher": "OpenAI",
        "public_observation": "Displayed trajectories have non-planar axial extent along the vortex axis.",
        "numerical_target": None,
    },
)

METRIC_PROVENANCE = {
    "streamline_protocol_prs": [583, 601],
    "streamline_protocol": "48 fixed seeds; normalized-velocity bidirectional RK4; 2.25 arclength/side; 91 states/side",
    "public_proxy_audit_pr": 1087,
    "public_proxy_gap": "the admitted inward+swirl evidence was a ring proxy, not a streamline proof",
    "classification": "internal_protocol_reimplementation",
    "new_external_method_migrated": False,
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
    "streamline_trajectory_measured": True,
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


def seed_points() -> np.ndarray:
    points: list[tuple[float, float, float]] = []
    for radius in SEED_RADII:
        for z in SEED_Z:
            for index in range(SEED_AZIMUTHS):
                theta = 2.0 * np.pi * index / SEED_AZIMUTHS
                points.append((radius * np.cos(theta), radius * np.sin(theta), z))
    out = np.asarray(points, dtype=np.float64)
    if out.shape != (48, 3):
        raise AssertionError("frozen seed contract drift")
    return out


def sample_velocity_grid(field: VelocityField) -> tuple[np.ndarray, np.ndarray]:
    axis = np.linspace(BOX[0], BOX[1], GRID_RESOLUTION, dtype=np.float64)
    xx, yy, zz = np.meshgrid(axis, axis, axis, indexing="ij")
    values = np.asarray(
        field.velocity(xx.ravel(), yy.ravel(), zz.ravel(), TIME), dtype=np.float64
    )
    if values.shape != (GRID_RESOLUTION**3, 3) or not np.isfinite(values).all():
        raise ValueError("velocity must return finite shape-(n^3,3) data")
    if not np.any(values):
        raise ValueError("zero sampled field is not a valid trajectory diagnostic")
    return axis, values.reshape(GRID_RESOLUTION, GRID_RESOLUTION, GRID_RESOLUTION, 3)


def grid_interpolator(axis: np.ndarray, field: np.ndarray):
    coords = np.asarray(axis, dtype=np.float64)
    values = np.asarray(field, dtype=np.float64)
    expected = (len(coords), len(coords), len(coords), 3)
    if values.shape != expected or not np.isfinite(values).all():
        raise ValueError("invalid grid field")
    components = [
        RegularGridInterpolator(
            (coords, coords, coords), values[..., component], method="linear",
            bounds_error=False, fill_value=0.0
        )
        for component in range(3)
    ]

    def evaluate(points: np.ndarray) -> np.ndarray:
        pts = np.asarray(points, dtype=np.float64)
        one = pts.ndim == 1
        pts = np.atleast_2d(pts)
        out = np.column_stack([fn(pts) for fn in components])
        return out[0] if one else out

    return evaluate


def _normalized_direction(evaluate, points: np.ndarray, sign: float) -> np.ndarray:
    values = np.asarray(evaluate(points), dtype=np.float64)
    values = np.atleast_2d(values)
    speeds = np.linalg.norm(values, axis=1)
    out = np.zeros_like(values)
    active = speeds >= SPEED_FLOOR
    out[active] = float(sign) * values[active] / speeds[active, None]
    return out


def _rk4_side(evaluate, seeds: np.ndarray, sign: float) -> np.ndarray:
    states = np.asarray(seeds, dtype=np.float64).copy()
    step = STREAM_ARCLENGTH / float(STREAM_OUTPUTS_PER_SIDE - 1)
    history = [states.copy()]
    active = np.max(np.abs(states), axis=1) < ACTIVE_BOX_MAX
    for _ in range(STREAM_OUTPUTS_PER_SIDE - 1):
        k1 = _normalized_direction(evaluate, states, sign)
        k2 = _normalized_direction(evaluate, states + 0.5 * step * k1, sign)
        k3 = _normalized_direction(evaluate, states + 0.5 * step * k2, sign)
        k4 = _normalized_direction(evaluate, states + step * k3, sign)
        proposal = states + (step / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
        proposal_active = np.max(np.abs(proposal), axis=1) < ACTIVE_BOX_MAX
        keep = active & proposal_active
        states[keep] = proposal[keep]
        active = keep
        history.append(states.copy())
    return np.stack(history, axis=1)


def integrate_streamlines(evaluate, seeds: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    backward = _rk4_side(evaluate, seeds, -1.0)
    forward = _rk4_side(evaluate, seeds, +1.0)
    full = np.concatenate((backward[:, :0:-1, :], forward), axis=1)
    return full, forward


def _turn_count(line: np.ndarray) -> float:
    angle = np.unwrap(np.arctan2(line[:, 1], line[:, 0]))
    return float(np.sum(np.abs(np.diff(angle))) / (2.0 * np.pi))


def full_line_metrics(lines: np.ndarray) -> dict[str, float | int]:
    turns: list[float] = []
    axial_spans: list[float] = []
    radial_spans: list[float] = []
    lengths: list[float] = []
    for line in np.asarray(lines, dtype=np.float64):
        radius = np.hypot(line[:, 0], line[:, 1])
        turns.append(_turn_count(line))
        axial_spans.append(float(np.ptp(line[:, 2])))
        radial_spans.append(float(np.ptp(radius)))
        lengths.append(float(np.sum(np.linalg.norm(np.diff(line, axis=0), axis=1))))
    return {
        "line_count": int(len(lines)),
        "mean_absolute_turns": float(np.mean(turns)),
        "max_absolute_turns": float(np.max(turns)),
        "mean_axial_span": float(np.mean(axial_spans)),
        "mean_radial_span": float(np.mean(radial_spans)),
        "mean_rendered_length": float(np.mean(lengths)),
    }


def _forward_subset_metrics(forward: np.ndarray, seeds: np.ndarray) -> dict[str, float | int]:
    fwd = np.asarray(forward, dtype=np.float64)
    src = np.asarray(seeds, dtype=np.float64)
    if fwd.ndim != 3 or fwd.shape[0] != src.shape[0] or fwd.shape[2] != 3:
        raise ValueError("forward streamline/seed shape mismatch")
    seed_radius = np.hypot(src[:, 0], src[:, 1])
    end_radius = np.hypot(fwd[:, -1, 0], fwd[:, -1, 1])
    displacement = end_radius - seed_radius
    turns = np.asarray([_turn_count(line) for line in fwd], dtype=np.float64)
    axial_excursion = np.asarray(
        [np.max(np.abs(line[:, 2] - seed[2])) for line, seed in zip(fwd, src)], dtype=np.float64
    )
    return {
        "line_count": int(len(fwd)),
        "mean_endpoint_radial_displacement": float(np.mean(displacement)),
        "min_endpoint_radial_displacement": float(np.min(displacement)),
        "max_endpoint_radial_displacement": float(np.max(displacement)),
        "inward_endpoint_fraction": float(np.mean(displacement < 0.0)),
        "mean_absolute_forward_turns": float(np.mean(turns)),
        "mean_axial_excursion": float(np.mean(axial_excursion)),
    }


def trajectory_metrics(field: VelocityField) -> dict[str, Any]:
    axis, grid = sample_velocity_grid(field)
    evaluate = grid_interpolator(axis, grid)
    seeds = seed_points()
    full, forward = integrate_streamlines(evaluate, seeds)
    overall = _forward_subset_metrics(forward, seeds)
    upper_mask = seeds[:, 2] > 0.0
    lower_mask = seeds[:, 2] < 0.0
    upper = _forward_subset_metrics(forward[upper_mask], seeds[upper_mask])
    lower = _forward_subset_metrics(forward[lower_mask], seeds[lower_mask])
    complete = full_line_metrics(full)

    side_rows = (upper, lower)
    inward = bool(
        overall["mean_endpoint_radial_displacement"] < 0.0
        and all(row["mean_endpoint_radial_displacement"] < 0.0 for row in side_rows)
    )
    winding = bool(
        overall["mean_absolute_forward_turns"] > NONDEGENERACY_FLOOR
        and all(row["mean_absolute_forward_turns"] > NONDEGENERACY_FLOOR for row in side_rows)
    )
    axial = bool(
        overall["mean_axial_excursion"] > NONDEGENERACY_FLOOR
        and all(row["mean_axial_excursion"] > NONDEGENERACY_FLOOR for row in side_rows)
    )
    return {
        "complete_lines": complete,
        "forward_overall": overall,
        "forward_upper_seeds": upper,
        "forward_lower_seeds": lower,
        "proxies": {
            "inward_forward_trajectory_proxy": inward,
            "nonzero_spiral_winding_proxy": winding,
            "nonplanar_axial_trajectory_proxy": axial,
            "coarse_inward_spiraling_axial_structure_proxy": bool(inward and winding and axial),
        },
    }


def _protocol() -> dict[str, Any]:
    return {
        "time": TIME,
        "grid_resolution": GRID_RESOLUTION,
        "box": list(BOX),
        "seed_radii": list(SEED_RADII),
        "seed_z": list(SEED_Z),
        "seed_azimuths": SEED_AZIMUTHS,
        "stream_arclength_per_side": STREAM_ARCLENGTH,
        "stream_outputs_per_side": STREAM_OUTPUTS_PER_SIDE,
        "speed_floor": SPEED_FLOOR,
        "active_box_max": ACTIVE_BOX_MAX,
        "nondegeneracy_floor": NONDEGENERACY_FLOOR,
        "interpolation": "scipy RegularGridInterpolator linear; fill_value=0 outside box",
        "integrator": "normalized-velocity deterministic RK4; bidirectional complete lines plus physical-velocity forward branch",
        "source_numeric_targets_used": False,
        "renderer_or_camera_used": False,
        "pixel_loss_used": False,
    }


def _decision(metrics: dict[str, Any]) -> dict[str, Any]:
    proxies = metrics["proxies"]
    failures = [name for name, passed in proxies.items() if name != "coarse_inward_spiraling_axial_structure_proxy" and not passed]
    return {
        "candidate_side_trajectory_structure_defect_established": bool(failures),
        "failed_structural_proxies": failures,
        "new_basis_authorized": False,
        "coefficient_change_authorized": False,
        "candidate_mutation_authorized": False,
        "actual_velocity_changed": False,
        "basis_dimension_change": 0,
        "direct_visualization_fingerprint_improvement": 0.0,
        "next_route": (
            "coarse absence-of-trajectory-structure triggers closed; any future mutation requires a separately bound magnitude/profile/trajectory discrepancy"
            if not failures
            else "route the concrete streamline-structure deficit first through existing #1078 five-control sensitivity evidence before any basis growth"
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
        raise ValueError("unexpected stable streamline receipt schema/task")
    if receipt.get("preregister_issue") != PREREG_ISSUE or receipt.get("base_main") != BASE_MAIN:
        raise ValueError("preregistration/base drift")
    identity = receipt.get("stable_candidate_identity")
    materialization = receipt.get("materialization_join")
    measurement = receipt.get("measurement")
    truth = receipt.get("truth_boundary")
    if not all(isinstance(v, dict) for v in (identity, materialization, measurement, truth)):
        raise ValueError("malformed stable streamline receipt")
    expected_identity = {
        "candidate_id": CANDIDATE_ID,
        "candidate_semantic_identity_sha256": EXPECTED_STABLE_CANDIDATE_IDENTITY,
        "velocity_semantic_identity_sha256": EXPECTED_STABLE_VELOCITY_IDENTITY,
    }
    if identity != expected_identity:
        raise ValueError("stable candidate/callable identity drift")
    if receipt.get("protocol") != _protocol():
        raise ValueError("frozen streamline protocol drift")
    if tuple(receipt.get("public_observables", ())) != PUBLIC_OBSERVABLES:
        raise ValueError("public qualitative observable drift")
    if receipt.get("metric_provenance") != METRIC_PROVENANCE:
        raise ValueError("metric provenance drift")
    if truth != TRUTH_BOUNDARY:
        raise ValueError("truth boundary drift")
    for key in (
        "stable_identity_receipt_sha256",
        "legacy_whole_candidate_identity_sha256",
        "materialization_evidence_sha256",
    ):
        if not _is_sha256(materialization.get(key)):
            raise ValueError(f"invalid materialization join {key}")
    if materialization.get("included_in_measurement_binding") is not False:
        raise ValueError("materialization evidence leaked into semantic measurement binding")
    if materialization.get("stable_identity_main_merge") != STABLE_IDENTITY_MAIN_MERGE:
        raise ValueError("stable identity ancestry drift")
    if materialization.get("stable_axial_main_merge") != STABLE_AXIAL_MAIN_MERGE:
        raise ValueError("stable axial ancestry drift")

    metrics = measurement.get("metrics")
    if not isinstance(metrics, dict):
        raise ValueError("missing trajectory metrics")
    for section in ("complete_lines", "forward_overall", "forward_upper_seeds", "forward_lower_seeds"):
        if not isinstance(metrics.get(section), dict):
            raise ValueError(f"missing trajectory section {section}")
    if metrics["complete_lines"].get("line_count") != 48:
        raise ValueError("complete-line count drift")
    if metrics["forward_overall"].get("line_count") != 48:
        raise ValueError("forward-line count drift")
    if metrics["forward_upper_seeds"].get("line_count") != 24 or metrics["forward_lower_seeds"].get("line_count") != 24:
        raise ValueError("upper/lower seed split drift")
    for section in ("complete_lines", "forward_overall", "forward_upper_seeds", "forward_lower_seeds"):
        for key, value in metrics[section].items():
            if key != "line_count" and (not np.isfinite(float(value))):
                raise ValueError(f"nonfinite trajectory metric {section}.{key}")
    expected_decision = _decision(metrics)
    if receipt.get("decision") != expected_decision:
        raise ValueError("trajectory routing decision drift")
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
    """Measure fixed-seed trajectories and bind them to the admitted ST052 semantic identity."""
    constrained_root = Path(constrained_root).resolve()
    source_root = Path(source_root).resolve()
    bundle_dir = Path(bundle_dir).resolve()

    stable_receipt = stable_identity.execute(constrained_root=constrained_root, bundle_dir=bundle_dir)
    stable_identity.validate_receipt(stable_receipt)
    stable = stable_receipt["stable_identity"]
    candidate_semantic_id = stable["candidate_semantic_identity_sha256"]
    velocity_semantic_id = stable["velocity_semantic_identity_sha256"]
    if candidate_semantic_id != EXPECTED_STABLE_CANDIDATE_IDENTITY:
        raise ValueError("fresh stable candidate semantic identity drift")
    if velocity_semantic_id != EXPECTED_STABLE_VELOCITY_IDENTITY:
        raise ValueError("fresh stable callable identity drift")

    identity_bridge._install_constrained_package_path(constrained_root)
    whole = importlib.import_module("openai_ns_reconstruction.st052_linear_temporal_capsule")
    candidate = whole.load_bundle_runtime(bundle_dir, exact_source_root=source_root)
    legacy_id = str(candidate.identity_sha256)
    evidence = stable_receipt["materialization_evidence"]
    if legacy_id != evidence["legacy_whole_candidate_identity_sha256"]:
        raise ValueError("loaded callable materialization disagrees with stable-identity bundle")

    identity = {
        "candidate_id": CANDIDATE_ID,
        "candidate_semantic_identity_sha256": candidate_semantic_id,
        "velocity_semantic_identity_sha256": velocity_semantic_id,
    }
    measurement = {
        "time": TIME,
        "metrics": trajectory_metrics(candidate),
        "interpretation": "target-free candidate-side fixed-seed trajectory structure; not an OpenAI magnitude/profile match",
    }
    measurement_sha = _canonical_sha256(measurement)
    materialization_join = {
        "stable_identity_receipt_sha256": stable_receipt["receipt_sha256"],
        "legacy_whole_candidate_identity_sha256": legacy_id,
        "materialization_evidence_sha256": _canonical_sha256(evidence),
        "included_in_measurement_binding": False,
        "stable_identity_main_merge": STABLE_IDENTITY_MAIN_MERGE,
        "stable_axial_main_merge": STABLE_AXIAL_MAIN_MERGE,
    }
    receipt: dict[str, Any] = {
        "schema": SCHEMA,
        "task_id": TASK_ID,
        "preregister_issue": PREREG_ISSUE,
        "base_main": BASE_MAIN,
        "stable_candidate_identity": identity,
        "materialization_join": materialization_join,
        "protocol": _protocol(),
        "public_observables": list(PUBLIC_OBSERVABLES),
        "metric_provenance": METRIC_PROVENANCE,
        "measurement": measurement,
        "measurement_sha256": measurement_sha,
        "candidate_measurement_binding_sha256": _measurement_binding_sha256(identity, measurement_sha),
        "decision": _decision(measurement["metrics"]),
        "truth_boundary": TRUTH_BOUNDARY,
        "receipt_sha256": None,
    }
    unsigned = dict(receipt)
    unsigned.pop("receipt_sha256")
    receipt["receipt_sha256"] = _canonical_sha256(unsigned)
    validate_receipt(receipt)
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--constrained-root", required=True)
    parser.add_argument("--source-root", required=True)
    parser.add_argument("--bundle-dir", required=True)
    parser.add_argument("--report", required=True)
    args = parser.parse_args()
    receipt = execute(
        constrained_root=args.constrained_root,
        source_root=args.source_root,
        bundle_dir=args.bundle_dir,
    )
    path = Path(args.report)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "task_id": TASK_ID,
        "stable_candidate_identity": receipt["stable_candidate_identity"],
        "metrics": receipt["measurement"]["metrics"],
        "decision": receipt["decision"],
        "measurement_sha256": receipt["measurement_sha256"],
        "receipt_sha256": receipt["receipt_sha256"],
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
