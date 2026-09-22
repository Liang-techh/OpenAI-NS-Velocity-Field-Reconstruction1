"""Map the admitted ST052 inward-streamline defect through five existing controls.

Preregistered in issue #1230.  This increment is diagnostic/routing only: it
binds the successful stable-identity streamline receipt from PR #1212 to the
already-defined five Agent-7 control tangents from exact PR #1078.  It does not
add a basis, select a coefficient magnitude, mutate the candidate, or perform
PDE acceptance.
"""
from __future__ import annotations

import argparse
import importlib
import json
from pathlib import Path
from typing import Any

import numpy as np

from . import st052_identity_bound_morphology_execution as identity_bridge
from . import st052_stable_identity_bound_streamline_trajectory as trajectory
from . import st052_stable_semantic_identity as stable_identity

SCHEMA = "st052-stable-streamline-inward-control-sensitivity/v1"
TASK_ID = "CR003-ST052M-STABLE-STREAMLINE-CONTROL-SENSITIVITY-145"
PREREG_ISSUE = 1230
SOURCE_STREAMLINE_PR = 1212
SOURCE_STREAMLINE_HEAD = "4193677861bc2a322ecb3128e359dd528022624a"
SOURCE_STREAMLINE_RUN = 35739391629
SOURCE_STREAMLINE_ARTIFACT_ID = 10704265735
SOURCE_STREAMLINE_ARTIFACT_DIGEST = "sha256:c4207944219c54e16908a24b3723eedc479f1ff7ce591d9393e8711e32ea34d1"
SOURCE_STREAMLINE_MEASUREMENT_SHA256 = "e9324a890621efa2b2569d3caeb3e3e03356e1be0b10b15f8fbcca8436e7daef"
SOURCE_CONTROL_PR = 1078
SOURCE_CONTROL_HEAD = "849713ba95cde8002fc1002d8b6856ea3f1dc264"

CHANNELS = (
    "amplitude",
    "radial_shape",
    "axial_turnover",
    "temporal_curvature",
    "toroidal_swirl",
)
TIMES_FOR_TRANSFER = (0.25, 0.50, 0.75)
EPSILONS = (1.0e-3, 5.0e-4)
PRIMARY_EPSILON = 5.0e-4
SUPPORT_RADIUS = 2.0
TRANSFER_ABS_MAX = 1.0e-9
TRANSFER_REL_MAX = 1.0e-9
DERIVATIVE_DRIFT_MAX = 0.05
DERIVATIVE_COSINE_MIN = 0.999
RESPONSE_FLOOR = 1.0e-10

EXPECTED_BASELINE_DISPLACEMENT = {
    "overall": 0.02211242744995166,
    "upper": 0.022238527286895324,
    "lower": 0.021986327613007978,
}

TRUTH_BOUNDARY = {
    "candidate_changed": False,
    "velocity_coefficients_changed": False,
    "basis_dimension_changed": False,
    "new_spatial_basis_added": False,
    "new_temporal_basis_added": False,
    "coefficient_selected": False,
    "coefficient_magnitude_selected": False,
    "pressure_changed": False,
    "forcing_changed": False,
    "held_out_pde_residual_evaluated": False,
    "public_source_numeric_targets_used": False,
    "renderer_or_camera_used": False,
    "pixel_loss_used": False,
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


def _canonical_sha256(value: Any) -> str:
    return stable_identity.canonical_sha256(value)


def _protocol() -> dict[str, Any]:
    return {
        "source_streamline_pr": SOURCE_STREAMLINE_PR,
        "source_streamline_head": SOURCE_STREAMLINE_HEAD,
        "source_streamline_run": SOURCE_STREAMLINE_RUN,
        "source_streamline_artifact_id": SOURCE_STREAMLINE_ARTIFACT_ID,
        "source_streamline_artifact_digest": SOURCE_STREAMLINE_ARTIFACT_DIGEST,
        "source_streamline_measurement_sha256": SOURCE_STREAMLINE_MEASUREMENT_SHA256,
        "source_control_pr": SOURCE_CONTROL_PR,
        "source_control_head": SOURCE_CONTROL_HEAD,
        "channels": list(CHANNELS),
        "transfer_times": list(TIMES_FOR_TRANSFER),
        "transfer_abs_max": TRANSFER_ABS_MAX,
        "transfer_relative_max": TRANSFER_REL_MAX,
        "trajectory_protocol": trajectory._protocol(),
        "epsilons": list(EPSILONS),
        "primary_epsilon": PRIMARY_EPSILON,
        "support_radius_normalization": SUPPORT_RADIUS,
        "derivative_drift_max": DERIVATIVE_DRIFT_MAX,
        "derivative_cosine_min": DERIVATIVE_COSINE_MIN,
        "response_floor": RESPONSE_FLOOR,
        "public_openai_numeric_target": None,
    }


def _load_control_bundle(path: str | Path) -> dict[str, np.ndarray]:
    with np.load(Path(path), allow_pickle=False) as data:
        bundle = {name: np.asarray(data[name], dtype=np.float64) for name in data.files}
    required = {"seed_points", "baseline_seed_velocity", "axis", *[f"tangent_{c}" for c in CHANNELS]}
    if set(bundle) != required:
        raise ValueError(f"control bundle keys drifted: {sorted(bundle)}")
    seeds = trajectory.seed_points()
    if bundle["seed_points"].shape != seeds.shape or not np.array_equal(bundle["seed_points"], seeds):
        raise ValueError("control-bundle seed contract drift")
    if bundle["baseline_seed_velocity"].shape != (len(TIMES_FOR_TRANSFER), len(seeds), 3):
        raise ValueError("control baseline seed-velocity shape drift")
    axis = bundle["axis"]
    expected_axis = np.linspace(trajectory.BOX[0], trajectory.BOX[1], trajectory.GRID_RESOLUTION)
    if axis.shape != expected_axis.shape or not np.array_equal(axis, expected_axis):
        raise ValueError("control tangent grid axis drift")
    shape = (trajectory.GRID_RESOLUTION, trajectory.GRID_RESOLUTION, trajectory.GRID_RESOLUTION, 3)
    for channel in CHANNELS:
        grid = bundle[f"tangent_{channel}"]
        if grid.shape != shape or not np.isfinite(grid).all():
            raise ValueError(f"invalid tangent grid for {channel}")
    if not np.isfinite(bundle["baseline_seed_velocity"]).all():
        raise ValueError("nonfinite control baseline seed velocities")
    return bundle


def _load_control_metadata(path: str | Path) -> dict[str, Any]:
    metadata = json.loads(Path(path).read_text())
    expected = {
        "source_control_pr": SOURCE_CONTROL_PR,
        "source_control_head": SOURCE_CONTROL_HEAD,
        "channels": list(CHANNELS),
        "times": list(TIMES_FOR_TRANSFER),
        "time_for_tangent_grids": trajectory.TIME,
        "grid_resolution": trajectory.GRID_RESOLUTION,
    }
    for key, value in expected.items():
        if metadata.get(key) != value:
            raise ValueError(f"control metadata drift for {key}")
    return metadata


def _load_source_receipt(path: str | Path) -> dict[str, Any]:
    receipt = json.loads(Path(path).read_text())
    trajectory.validate_receipt(receipt)
    if receipt.get("measurement_sha256") != SOURCE_STREAMLINE_MEASUREMENT_SHA256:
        raise ValueError("source streamline measurement SHA drift")
    metrics = receipt["measurement"]["metrics"]
    actual = {
        "overall": float(metrics["forward_overall"]["mean_endpoint_radial_displacement"]),
        "upper": float(metrics["forward_upper_seeds"]["mean_endpoint_radial_displacement"]),
        "lower": float(metrics["forward_lower_seeds"]["mean_endpoint_radial_displacement"]),
    }
    for key, expected in EXPECTED_BASELINE_DISPLACEMENT.items():
        if not np.isclose(actual[key], expected, rtol=0.0, atol=1.0e-15):
            raise ValueError(f"source streamline displacement drift for {key}")
    if metrics["proxies"]["inward_forward_trajectory_proxy"] is not False:
        raise ValueError("source receipt no longer establishes inward-trajectory defect")
    if metrics["proxies"]["nonzero_spiral_winding_proxy"] is not True:
        raise ValueError("source receipt winding truth drift")
    if metrics["proxies"]["nonplanar_axial_trajectory_proxy"] is not True:
        raise ValueError("source receipt axial-trajectory truth drift")
    return receipt


def _load_stable_candidate(*, constrained_root: str | Path, source_root: str | Path, bundle_dir: str | Path):
    constrained_root = Path(constrained_root).resolve()
    source_root = Path(source_root).resolve()
    bundle_dir = Path(bundle_dir).resolve()
    stable_receipt = stable_identity.execute(constrained_root=constrained_root, bundle_dir=bundle_dir)
    stable_identity.validate_receipt(stable_receipt)
    stable = stable_receipt["stable_identity"]
    if stable["candidate_semantic_identity_sha256"] != trajectory.EXPECTED_STABLE_CANDIDATE_IDENTITY:
        raise ValueError("stable candidate semantic identity drift")
    if stable["velocity_semantic_identity_sha256"] != trajectory.EXPECTED_STABLE_VELOCITY_IDENTITY:
        raise ValueError("stable velocity semantic identity drift")
    identity_bridge._install_constrained_package_path(constrained_root)
    whole = importlib.import_module("openai_ns_reconstruction.st052_linear_temporal_capsule")
    candidate = whole.load_bundle_runtime(bundle_dir, exact_source_root=source_root)
    evidence = stable_receipt["materialization_evidence"]
    if str(candidate.identity_sha256) != evidence["legacy_whole_candidate_identity_sha256"]:
        raise ValueError("stable callable materialization mismatch")
    return candidate, stable_receipt


def _stable_seed_velocities(candidate) -> np.ndarray:
    seeds = trajectory.seed_points()
    values = []
    for time in TIMES_FOR_TRANSFER:
        values.append(np.asarray(candidate.velocity(seeds[:, 0], seeds[:, 1], seeds[:, 2], time), dtype=np.float64))
    out = np.stack(values, axis=0)
    if out.shape != (len(TIMES_FOR_TRANSFER), len(seeds), 3) or not np.isfinite(out).all():
        raise ValueError("invalid stable seed velocities")
    return out


def _compatibility(stable_values: np.ndarray, control_values: np.ndarray) -> dict[str, Any]:
    delta = np.asarray(control_values, dtype=np.float64) - np.asarray(stable_values, dtype=np.float64)
    point_abs = np.max(np.abs(delta), axis=-1)
    speed = np.linalg.norm(stable_values, axis=-1)
    point_rel = point_abs / np.maximum(1.0, speed)
    max_abs = float(np.max(point_abs))
    max_rel = float(np.max(point_rel))
    return {
        "max_absolute_velocity_mismatch": max_abs,
        "max_relative_velocity_mismatch": max_rel,
        "absolute_gate": TRANSFER_ABS_MAX,
        "relative_gate": TRANSFER_REL_MAX,
        "passes": bool(max_abs <= TRANSFER_ABS_MAX and max_rel <= TRANSFER_REL_MAX),
    }


def _forward_metrics_from_grid(axis: np.ndarray, grid: np.ndarray) -> dict[str, Any]:
    evaluate = trajectory.grid_interpolator(axis, grid)
    seeds = trajectory.seed_points()
    _full, forward = trajectory.integrate_streamlines(evaluate, seeds)
    upper = seeds[:, 2] > 0.0
    lower = seeds[:, 2] < 0.0
    return {
        "overall": trajectory._forward_subset_metrics(forward, seeds),
        "upper": trajectory._forward_subset_metrics(forward[upper], seeds[upper]),
        "lower": trajectory._forward_subset_metrics(forward[lower], seeds[lower]),
    }


def _coordinates(metrics: dict[str, Any]) -> np.ndarray:
    return np.asarray(
        [
            float(metrics["overall"]["mean_endpoint_radial_displacement"]) / SUPPORT_RADIUS,
            float(metrics["upper"]["mean_endpoint_radial_displacement"]) / SUPPORT_RADIUS,
            float(metrics["lower"]["mean_endpoint_radial_displacement"]) / SUPPORT_RADIUS,
        ],
        dtype=np.float64,
    )


def _collateral(metrics: dict[str, Any]) -> np.ndarray:
    return np.asarray(
        [
            float(metrics["overall"]["mean_absolute_forward_turns"]),
            float(metrics["overall"]["mean_axial_excursion"]) / SUPPORT_RADIUS,
        ],
        dtype=np.float64,
    )


def _derivative(plus: np.ndarray, minus: np.ndarray, epsilon: float) -> np.ndarray:
    out = (np.asarray(plus, dtype=np.float64) - np.asarray(minus, dtype=np.float64)) / (2.0 * float(epsilon))
    if not np.isfinite(out).all():
        raise ValueError("nonfinite local derivative")
    return out


def _channel_stability(fine: np.ndarray, coarse: np.ndarray) -> dict[str, Any]:
    fine = np.asarray(fine, dtype=np.float64)
    coarse = np.asarray(coarse, dtype=np.float64)
    nf = float(np.linalg.norm(fine))
    nc = float(np.linalg.norm(coarse))
    if nf <= RESPONSE_FLOOR or nc <= RESPONSE_FLOOR:
        drift = None
        cosine = None
        passes = False
    else:
        drift = float(np.linalg.norm(fine - coarse) / nf)
        cosine = float(np.dot(fine, coarse) / (nf * nc))
        passes = bool(drift <= DERIVATIVE_DRIFT_MAX and cosine >= DERIVATIVE_COSINE_MIN)
    return {
        "fine_l2": nf,
        "coarse_l2": nc,
        "relative_drift": drift,
        "coarse_fine_cosine": cosine,
        "drift_gate": DERIVATIVE_DRIFT_MAX,
        "cosine_gate": DERIVATIVE_COSINE_MIN,
        "passes": passes,
    }


def _coherent_leverage(fine: np.ndarray, stability: dict[str, Any]) -> tuple[bool, str | None]:
    vec = np.asarray(fine, dtype=np.float64)
    nonzero = bool(np.all(np.abs(vec) > RESPONSE_FLOOR))
    same_positive = bool(np.all(vec > RESPONSE_FLOOR))
    same_negative = bool(np.all(vec < -RESPONSE_FLOOR))
    coherent = bool(nonzero and (same_positive or same_negative) and stability["passes"])
    if not coherent:
        return False, None
    return True, ("negative" if same_positive else "positive")


def _sensitivity_audit(axis: np.ndarray, stable_grid: np.ndarray, control: dict[str, np.ndarray]) -> dict[str, Any]:
    baseline_metrics = _forward_metrics_from_grid(axis, stable_grid)
    baseline_coords = _coordinates(baseline_metrics)
    baseline_collateral = _collateral(baseline_metrics)
    derivatives: dict[float, dict[str, dict[str, np.ndarray]]] = {}
    for epsilon in EPSILONS:
        per_channel: dict[str, dict[str, np.ndarray]] = {}
        for channel in CHANNELS:
            tangent = control[f"tangent_{channel}"]
            plus_metrics = _forward_metrics_from_grid(axis, stable_grid + float(epsilon) * tangent)
            minus_metrics = _forward_metrics_from_grid(axis, stable_grid - float(epsilon) * tangent)
            per_channel[channel] = {
                "displacement": _derivative(_coordinates(plus_metrics), _coordinates(minus_metrics), epsilon),
                "collateral": _derivative(_collateral(plus_metrics), _collateral(minus_metrics), epsilon),
            }
        derivatives[float(epsilon)] = per_channel

    fine = derivatives[PRIMARY_EPSILON]
    coarse = derivatives[EPSILONS[0]]
    channels: dict[str, Any] = {}
    coherent_channels: list[str] = []
    for channel in CHANNELS:
        stability = _channel_stability(fine[channel]["displacement"], coarse[channel]["displacement"])
        coherent, reducing_sign = _coherent_leverage(fine[channel]["displacement"], stability)
        if coherent:
            coherent_channels.append(channel)
        channels[channel] = {
            "coarse_displacement_derivative": coarse[channel]["displacement"].tolist(),
            "fine_displacement_derivative": fine[channel]["displacement"].tolist(),
            "fine_collateral_derivative": {
                "mean_absolute_forward_turns": float(fine[channel]["collateral"][0]),
                "mean_axial_excursion_over_support_radius": float(fine[channel]["collateral"][1]),
            },
            "stability": stability,
            "coherent_existing_control_leverage": coherent,
            "local_coefficient_sign_that_reduces_positive_displacement": reducing_sign,
        }

    matrix = np.column_stack([fine[channel]["displacement"] for channel in CHANNELS])
    side_matrix = matrix[1:, :]
    singular = np.linalg.svd(side_matrix, compute_uv=False)
    rank = int(np.linalg.matrix_rank(side_matrix, tol=1.0e-10))
    condition = None if singular[-1] <= np.finfo(float).tiny else float(singular[0] / singular[-1])
    return {
        "baseline": {
            "displacement_coordinates_over_support_radius": baseline_coords.tolist(),
            "mean_absolute_forward_turns": float(baseline_collateral[0]),
            "mean_axial_excursion_over_support_radius": float(baseline_collateral[1]),
        },
        "fine_signed_displacement_jacobian_3x5": matrix.tolist(),
        "upper_lower_jacobian_2x5": side_matrix.tolist(),
        "upper_lower_singular_values": singular.tolist(),
        "upper_lower_rank": rank,
        "upper_lower_condition": condition,
        "channels": channels,
        "coherent_existing_control_channels": coherent_channels,
        "existing_control_leverage_established": bool(coherent_channels),
    }


def _decision(compatibility: dict[str, Any], audit: dict[str, Any] | None) -> dict[str, Any]:
    if not compatibility["passes"]:
        route = "baseline transfer failed; rebuild the five existing controls directly on the stable callable identity before any basis conclusion"
        obstruction = False
        existing = False
    elif audit is not None and audit["existing_control_leverage_established"]:
        route = "at least one existing control has stable coherent local inward-trajectory leverage; test an existing aligned control before any basis growth"
        obstruction = False
        existing = True
    else:
        route = "all five existing controls lack stable coherent local leverage on the admitted inward-trajectory defect; a later targeted inward-trajectory basis screen is permitted, but no basis is added here"
        obstruction = True
        existing = False
    return {
        "baseline_transfer_passed": bool(compatibility["passes"]),
        "existing_control_leverage_established": existing,
        "local_five_control_trajectory_obstruction_established": obstruction,
        "new_basis_added": False,
        "new_basis_authorized_in_this_increment": False,
        "coefficient_selected": False,
        "candidate_mutation_authorized": False,
        "actual_velocity_changed": False,
        "basis_dimension_change": 0,
        "direct_visualization_fingerprint_improvement": 0.0,
        "next_route": route,
    }


def validate_receipt(receipt: dict[str, Any]) -> None:
    if receipt.get("schema") != SCHEMA or receipt.get("task_id") != TASK_ID:
        raise ValueError("unexpected control-sensitivity receipt schema/task")
    if receipt.get("preregister_issue") != PREREG_ISSUE:
        raise ValueError("preregistration identity drift")
    if receipt.get("protocol") != _protocol():
        raise ValueError("frozen control-sensitivity protocol drift")
    if receipt.get("truth_boundary") != TRUTH_BOUNDARY:
        raise ValueError("truth boundary drift")
    source = receipt.get("source_streamline_receipt")
    compatibility = receipt.get("baseline_transfer")
    if not isinstance(source, dict) or not isinstance(compatibility, dict):
        raise ValueError("malformed source/transfer section")
    if source.get("measurement_sha256") != SOURCE_STREAMLINE_MEASUREMENT_SHA256:
        raise ValueError("source measurement binding drift")
    audit = receipt.get("trajectory_control_sensitivity")
    expected_decision = _decision(compatibility, audit)
    if receipt.get("decision") != expected_decision:
        raise ValueError("routing decision drift")
    if compatibility.get("passes"):
        if not isinstance(audit, dict):
            raise ValueError("missing sensitivity audit after compatible transfer")
        if len(audit.get("coherent_existing_control_channels", [])) > len(CHANNELS):
            raise ValueError("invalid coherent-control list")
        for channel in CHANNELS:
            row = audit["channels"].get(channel)
            if not isinstance(row, dict):
                raise ValueError(f"missing channel result {channel}")
    else:
        if audit is not None:
            raise ValueError("sensitivity result emitted despite failed baseline transfer")
    if receipt.get("receipt_sha256") is not None:
        unsigned = dict(receipt)
        recorded = unsigned.pop("receipt_sha256")
        if recorded != _canonical_sha256(unsigned):
            raise ValueError("receipt checksum mismatch")


def execute(
    *,
    constrained_root: str | Path,
    source_root: str | Path,
    bundle_dir: str | Path,
    source_receipt: str | Path,
    control_bundle: str | Path,
    control_metadata: str | Path,
) -> dict[str, Any]:
    source = _load_source_receipt(source_receipt)
    controls = _load_control_bundle(control_bundle)
    control_meta = _load_control_metadata(control_metadata)
    candidate, stable_receipt = _load_stable_candidate(
        constrained_root=constrained_root, source_root=source_root, bundle_dir=bundle_dir
    )
    stable_seed = _stable_seed_velocities(candidate)
    compatibility = _compatibility(stable_seed, controls["baseline_seed_velocity"])

    audit = None
    if compatibility["passes"]:
        axis, stable_grid = trajectory.sample_velocity_grid(candidate)
        if not np.array_equal(axis, controls["axis"]):
            raise ValueError("stable/control grid axis mismatch")
        audit = _sensitivity_audit(axis, stable_grid, controls)

    receipt: dict[str, Any] = {
        "schema": SCHEMA,
        "task_id": TASK_ID,
        "preregister_issue": PREREG_ISSUE,
        "protocol": _protocol(),
        "stable_candidate_identity": {
            "candidate_id": trajectory.CANDIDATE_ID,
            "candidate_semantic_identity_sha256": trajectory.EXPECTED_STABLE_CANDIDATE_IDENTITY,
            "velocity_semantic_identity_sha256": trajectory.EXPECTED_STABLE_VELOCITY_IDENTITY,
            "stable_identity_receipt_sha256": stable_receipt["receipt_sha256"],
        },
        "source_streamline_receipt": {
            "run_id": SOURCE_STREAMLINE_RUN,
            "artifact_id": SOURCE_STREAMLINE_ARTIFACT_ID,
            "artifact_digest": SOURCE_STREAMLINE_ARTIFACT_DIGEST,
            "measurement_sha256": source["measurement_sha256"],
            "receipt_sha256": source["receipt_sha256"],
            "baseline_mean_endpoint_radial_displacement": dict(EXPECTED_BASELINE_DISPLACEMENT),
        },
        "source_control_lineage": control_meta,
        "baseline_transfer": compatibility,
        "trajectory_control_sensitivity": audit,
        "decision": _decision(compatibility, audit),
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
    parser.add_argument("--source-receipt", required=True)
    parser.add_argument("--control-bundle", required=True)
    parser.add_argument("--control-metadata", required=True)
    parser.add_argument("--report", required=True)
    args = parser.parse_args()
    receipt = execute(
        constrained_root=args.constrained_root,
        source_root=args.source_root,
        bundle_dir=args.bundle_dir,
        source_receipt=args.source_receipt,
        control_bundle=args.control_bundle,
        control_metadata=args.control_metadata,
    )
    path = Path(args.report)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "task_id": TASK_ID,
        "baseline_transfer": receipt["baseline_transfer"],
        "trajectory_control_sensitivity": receipt["trajectory_control_sensitivity"],
        "decision": receipt["decision"],
        "receipt_sha256": receipt["receipt_sha256"],
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
