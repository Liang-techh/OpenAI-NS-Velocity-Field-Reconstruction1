"""Actual-field replay on a disjoint composite-midpoint spatial grid.

The selected dual-grid tangent candidate and the prior feasible candidate are
evaluated with the same full finite-difference field replay.  The grid uses
positive cylindrical cell weights, composite midpoints, and a new angular
offset, so this is a separate sampled spatial check rather than another fit
grid.  No PDE, trajectory, or scale-recursion acceptance is inferred.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree


ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from affine_momentum import jets, momentum  # noqa: E402
from enriched_shape_replay import build_field  # noqa: E402


DUAL_CANDIDATE = ROOT / "dual_grid_peak_tangent.json"
OLD_CANDIDATE = ROOT / "enriched_mean_endpoint_tangent.json"
FROZEN_GEOMETRY = ROOT / "full_wave_frozen_cache.json"
DENSE_REPORT = ROOT / "full_wave_dense_tangent.json"
REFINED_CACHE = ROOT / "refined_wave_momentum_cache.npz"
OUTPUT_PATH = ROOT / "dual_grid_disjoint_replay.json"

AXIAL_ORDER = 11
RADIAL_ORDER = 9
ANGLES = 12
ANGLE_SHIFT = 0.137


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _metric(residual, weights):
    residual = np.asarray(residual, dtype=float)
    weights = np.asarray(weights, dtype=float)
    magnitudes = np.linalg.norm(residual, axis=1)
    index = int(np.argmax(magnitudes))
    return {
        "point_count": int(len(residual)),
        "momentum_max": float(magnitudes[index]),
        "momentum_max_index": index,
        "momentum_max_point": [],
        "momentum_max_vector": residual[index].tolist(),
        "momentum_volume_L2": float(
            np.linalg.norm(residual.reshape(-1) * np.repeat(np.sqrt(weights), 3))
        ),
        "momentum_volume_RMS": float(
            np.sqrt(np.sum(weights * np.sum(residual * residual, axis=1)) / np.sum(weights))
        ),
        "physical_volume": float(np.sum(weights)),
    }


def _panel_cuts(mean, center, widths, tau, breaks):
    """Replicate frozen radial support panels without reusing quadrature nodes."""

    from joined_field import coordinates

    rlo, rhi = np.asarray(center, dtype=float)[0] + np.array([-1.0, 1.0]) * float(widths[0])
    q = float(coordinates(0.0, float(center[1]) / np.sqrt(mean.nu), tau, mean.inner.h)["q"])
    ri = np.sqrt(2.0 * mean.nu * q * mean.join_X)
    cuts = sorted(
        set(
            [float(rlo), float(rhi)]
            + [
                float(r)
                for r in ri * (1.0 + (mean.ratio - 1.0) * np.asarray(breaks, dtype=float))
                if rlo < r < rhi
            ]
        )
    )
    return cuts


def _composite_midpoint_grid(mean, center, widths, tau, breaks):
    """Return a positive-weight grid with no Gauss nodes.

    Axial cells are uniform composite midpoints.  Each frozen radial support
    panel is subdivided uniformly and sampled at radial midpoints.  We assign
    each cell its exact cylindrical annular volume divided by the angular
    count; this keeps every weight positive and the total domain volume exact.
    """

    center = np.asarray(center, dtype=float)
    widths = np.asarray(widths, dtype=float)
    zlo, zhi = center[1] - widths[1], center[1] + widths[1]
    z_edges = np.linspace(zlo, zhi, AXIAL_ORDER + 1)
    z_nodes = 0.5 * (z_edges[:-1] + z_edges[1:])
    dz = float(z_edges[1] - z_edges[0])
    cuts = _panel_cuts(mean, center, widths, tau, breaks)
    theta = ANGLE_SHIFT + np.arange(ANGLES, dtype=float) * 2.0 * np.pi / ANGLES
    points = []
    weights = []
    for z in z_nodes:
        for lo, hi in zip(cuts[:-1], cuts[1:]):
            radial_edges = np.linspace(lo, hi, RADIAL_ORDER + 1)
            for rlo, rhi in zip(radial_edges[:-1], radial_edges[1:]):
                r = 0.5 * (rlo + rhi)
                cell_weight = np.pi * (rhi * rhi - rlo * rlo) * dz / ANGLES
                points.extend(
                    np.column_stack(
                        (
                            np.full(ANGLES, r * 0.0 + r) * np.cos(theta),
                            np.full(ANGLES, r * 0.0 + r) * np.sin(theta),
                            np.full(ANGLES, z),
                        )
                    )
                )
                weights.extend([cell_weight] * ANGLES)
    points = np.asarray(points, dtype=float)
    weights = np.asarray(weights, dtype=float)
    expected_volume = np.pi * ((center[0] + widths[0]) ** 2 - (center[0] - widths[0]) ** 2) * (2.0 * widths[1])
    if points.shape != (AXIAL_ORDER * (len(cuts) - 1) * RADIAL_ORDER * ANGLES, 3):
        raise ValueError(f"Unexpected disjoint grid shape {points.shape}")
    if np.any(weights <= 0.0):
        raise ValueError("Composite midpoint grid has nonpositive weights")
    if not np.isclose(np.sum(weights), expected_volume, rtol=1e-12, atol=0.0):
        raise ValueError("Composite midpoint weights do not sum to the physical volume")
    return points, weights, cuts, expected_volume


def _actual_replay(candidate, points, frozen):
    field, snapshot = build_field(candidate)
    if snapshot["inputs"]["mean"]["coefficients"] != frozen["inputs"]["mean"]["coefficients"]:
        raise ValueError("Replay mean differs from frozen mean")
    tau = float(snapshot["inputs"]["mean"]["tau"])
    hspace = float(snapshot["timesteps"]["hspace"])
    htime = float(snapshot["timesteps"]["htime"])
    residual = momentum(jets(field, points, tau, hspace, htime))
    return residual, {"tau": tau, "hspace": hspace, "htime": htime}


def _set_peak_location(metric, points):
    index = int(metric["momentum_max_index"])
    metric["momentum_max_point"] = np.asarray(points[index], dtype=float).tolist()
    return metric


def run(output_path=OUTPUT_PATH):
    started = time.perf_counter()
    output_path = Path(output_path)
    dual = json.loads(DUAL_CANDIDATE.read_text(encoding="utf-8"))
    old = json.loads(OLD_CANDIDATE.read_text(encoding="utf-8"))
    frozen = json.loads(FROZEN_GEOMETRY.read_text(encoding="utf-8"))
    dense = json.loads(DENSE_REPORT.read_text(encoding="utf-8"))
    if dual.get("status") != "completed" or not dual.get("assembled_feasible"):
        raise ValueError("Dual-grid candidate is not assembled feasible")
    if old.get("status") != "completed" or not old.get("assembled_feasible"):
        raise ValueError("Old candidate is not assembled feasible")
    if not dense.get("new_frozen_cache"):
        raise ValueError("Dense report has no prior order-13 grid")

    center = np.asarray(frozen["inputs"]["wave"]["center"], dtype=float)
    widths = np.asarray(frozen["inputs"]["wave"]["widths"], dtype=float)
    tau = float(frozen["inputs"]["mean"]["tau"])
    breaks = list(frozen["inputs"]["geometry"]["radial_breaks"])
    mean, _ = __import__("broad_meridional_constrained", fromlist=["load_saved_field"]).load_saved_field()
    points, weights, cuts, expected_volume = _composite_midpoint_grid(
        mean, center, widths, tau, breaks
    )
    report = {
        "status": "running",
        "accepted": False,
        "pde_validated": False,
        "scale_recursion_established": False,
        "constraints_maintained": False,
        "scope": (
            "Actual full finite-difference replay of the dual-grid selected "
            "candidate and old feasible reference on a new composite-midpoint "
            "cylindrical grid. No optimization, PDE, trajectory, or recursion "
            "acceptance is claimed."
        ),
        "sources": {
            "dual_candidate": {"path": DUAL_CANDIDATE.name, "sha256": _sha(DUAL_CANDIDATE)},
            "old_candidate": {"path": OLD_CANDIDATE.name, "sha256": _sha(OLD_CANDIDATE)},
            "frozen_geometry": {"path": FROZEN_GEOMETRY.name, "sha256": _sha(FROZEN_GEOMETRY)},
            "dense_prior_grid": {"path": DENSE_REPORT.name, "sha256": _sha(DENSE_REPORT)},
            "refined_prior_grid": {"path": REFINED_CACHE.name, "sha256": _sha(REFINED_CACHE)},
        },
        "grid": {
            "construction": "uniform composite midpoints in axial cells and every frozen radial support panel",
            "axial_order": AXIAL_ORDER,
            "radial_order_per_panel": RADIAL_ORDER,
            "angles": ANGLES,
            "angle_shift": ANGLE_SHIFT,
            "radial_panel_count": len(cuts) - 1,
            "point_count": int(len(points)),
            "weight_sum": float(np.sum(weights)),
            "expected_physical_volume": float(expected_volume),
            "minimum_weight": float(np.min(weights)),
            "all_weights_positive": bool(np.all(weights > 0.0)),
            "radial_cuts": [float(v) for v in cuts],
            "source_definition": "new disjoint grid; not a Gauss-Legendre node set",
        },
        "replay_method": {
            "field": "enriched_shape_replay.build_field",
            "jets": "affine_momentum.jets five-point spatial/temporal FD",
            "same_grid_for_candidates": True,
            "base_mean_reused": False,
            "full_actual_replay_selected_and_reference": True,
        },
    }
    output_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"stage": "grid_built", "points": len(points), "weight_sum": float(np.sum(weights))}), flush=True)

    # Verify separation from the two prior sampled grids before replay.  The
    # nearest-neighbour distances are recorded rather than relying only on the
    # different quadrature labels.
    prior_dense = np.asarray(dense["new_frozen_cache"]["points"], dtype=float)
    with np.load(REFINED_CACHE, allow_pickle=False) as cache:
        prior_refined = np.asarray(cache["points"], dtype=float)
    dense_distance = cKDTree(prior_dense).query(points, k=1)[0]
    refined_distance = cKDTree(prior_refined).query(points, k=1)[0]
    exact_dense = int(np.sum(dense_distance <= 1.0e-14))
    exact_refined = int(np.sum(refined_distance <= 1.0e-14))
    report["disjoint_check"] = {
        "prior_dense_point_count": int(len(prior_dense)),
        "prior_refined_point_count": int(len(prior_refined)),
        "nearest_distance_to_dense_min": float(np.min(dense_distance)),
        "nearest_distance_to_refined_min": float(np.min(refined_distance)),
        "exact_coordinate_matches_at_1e-14_dense": exact_dense,
        "exact_coordinate_matches_at_1e-14_refined": exact_refined,
        "passes_exact_match_check": bool(exact_dense == 0 and exact_refined == 0),
    }
    if exact_dense or exact_refined:
        raise ValueError("Disjoint grid unexpectedly shares prior coordinates")
    report["status"] = "grid_disjoint_check_complete"
    output_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"stage": "disjoint_check_complete",
                      "min_dense_distance": float(np.min(dense_distance)),
                      "min_refined_distance": float(np.min(refined_distance))}), flush=True)

    print(json.dumps({"stage": "actual_replays_started", "points": len(points)}), flush=True)
    dual_residual, steps = _actual_replay(dual, points, frozen)
    old_residual, _ = _actual_replay(old, points, frozen)
    dual_metric = _set_peak_location(_metric(dual_residual, weights), points)
    old_metric = _set_peak_location(_metric(old_residual, weights), points)
    report["replays"] = {
        "dual_selected": dual_metric,
        "old_reference": old_metric,
        "finite_difference_steps": steps,
    }
    report["comparison"] = {
        "selected_minus_old_L2": float(dual_metric["momentum_volume_L2"] - old_metric["momentum_volume_L2"]),
        "selected_minus_old_peak": float(dual_metric["momentum_max"] - old_metric["momentum_max"]),
        "selected_L2_improves": bool(dual_metric["momentum_volume_L2"] < old_metric["momentum_volume_L2"]),
        "selected_peak_regresses": bool(dual_metric["momentum_max"] > old_metric["momentum_max"]),
        "adoption_decision": "do_not_adopt_if_peak_regresses",
    }
    report["status"] = "completed"
    report["elapsed_seconds"] = float(time.perf_counter() - started)
    output_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": report["status"],
        "dual_l2": dual_metric["momentum_volume_L2"],
        "old_l2": old_metric["momentum_volume_L2"],
        "dual_max": dual_metric["momentum_max"],
        "old_max": old_metric["momentum_max"],
        "selected_peak_regresses": report["comparison"]["selected_peak_regresses"],
        "elapsed_seconds": report["elapsed_seconds"],
    }), flush=True)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT_PATH)
    args = parser.parse_args()
    run(args.output)
