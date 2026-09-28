"""Independent finite-difference replay of the constrained acceleration candidate.

The 18,720-point grid is copied from ``full_wave_dense_tangent.json``'s
``new_frozen_cache``.  It is used only as a holdout grid: the balanced parent
and the selected constrained acceleration candidate are evaluated directly
with the Cartesian five-point finite-difference momentum implementation.
No cached momentum/design response is used to form either reported residual.
"""

from __future__ import annotations

import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from affine_momentum import jets, momentum  # noqa: E402
from broad_meridional_constrained import load_saved_field  # noqa: E402
from endpoint_acceleration_projection import (  # noqa: E402
    AccelerationCorrectionField,
    MixedAffineField,
    _metric,
    _unpack_acceleration_control,
    _unpack_mixed_affine_control,
    _unpack_wave,
)
from grouped_joined_field import install_in_field  # noqa: E402


OUTPUT_PATH = ROOT / "constrained_acceleration_holdout.json"
GRID_PATH = ROOT / "full_wave_dense_tangent.json"
BALANCED_PATH = ROOT / "balanced_refined_tangent.json"
FROZEN_PATH = ROOT / "full_wave_frozen_cache.json"
ACCELERATION_PATH = ROOT / "shape_constrained_acceleration.json"
PROJECTION_PATH = ROOT / "endpoint_acceleration_projection.py"

DELTA_K = 1.0e-6


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _save(path: Path, report: dict) -> None:
    path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")


def _array_sha256(*arrays: np.ndarray) -> str:
    digest = hashlib.sha256()
    for array in arrays:
        value = np.ascontiguousarray(np.asarray(array))
        digest.update(str(value.dtype).encode("ascii"))
        digest.update(np.asarray(value.shape, dtype=np.int64).tobytes())
        digest.update(value.tobytes())
    return digest.hexdigest()


def _max_detail(value, points, weights):
    magnitude = np.linalg.norm(np.asarray(value), axis=1)
    index = int(np.argmax(magnitude))
    return {
        "index": index,
        "point": np.asarray(points[index], dtype=float).tolist(),
        "weight": float(weights[index]),
        "vector": np.asarray(value[index], dtype=float).tolist(),
        "norm": float(magnitude[index]),
    }


def _load_holdout():
    source = json.loads(GRID_PATH.read_text(encoding="utf-8"))
    cache = source.get("new_frozen_cache")
    if not isinstance(cache, dict):
        raise ValueError("full_wave_dense_tangent.json has no new_frozen_cache")
    points = np.asarray(cache["points"], dtype=float)
    weights = np.asarray(cache["weights"], dtype=float)
    frozen_residual = np.asarray(cache["residual"], dtype=float)
    if points.shape != (18720, 3):
        raise ValueError(f"Expected holdout points (18720,3), got {points.shape}")
    if weights.shape != (18720,):
        raise ValueError(f"Expected holdout weights (18720,), got {weights.shape}")
    if frozen_residual.shape != (18720, 3):
        raise ValueError(
            f"Expected cached reference residual (18720,3), got {frozen_residual.shape}"
        )
    if not np.all(np.isfinite(points)) or not np.all(np.isfinite(weights)):
        raise ValueError("Holdout grid contains nonfinite points or weights")
    if np.any(weights <= 0.0):
        raise ValueError("Holdout weights must be positive")
    return source, cache, points, weights, frozen_residual


def _load_inputs():
    balanced = json.loads(BALANCED_PATH.read_text(encoding="utf-8"))
    frozen = json.loads(FROZEN_PATH.read_text(encoding="utf-8"))
    constrained = json.loads(ACCELERATION_PATH.read_text(encoding="utf-8"))
    if balanced.get("status") != "completed":
        raise ValueError("balanced_refined_tangent.json is not complete")
    if constrained.get("status") != "completed":
        raise ValueError("shape_constrained_acceleration.json is not complete")
    selected = constrained.get("selected")
    if not isinstance(selected, dict):
        raise ValueError("shape constrained report has no selected candidate")
    acceleration_coefficients = np.asarray(
        selected.get("acceleration_coefficients"), dtype=float
    )
    if acceleration_coefficients.shape != (324,):
        raise ValueError(
            "Expected selected acceleration coefficients with shape (324,), "
            f"got {acceleration_coefficients.shape}"
        )
    balanced_selected = balanced.get("selected")
    if not isinstance(balanced_selected, dict):
        raise ValueError("balanced report has no selected parent")
    base_control = np.asarray(
        balanced_selected.get("tangent_coefficients"), dtype=float
    )
    if base_control.shape != (264,):
        raise ValueError(
            "Expected balanced tangent coefficients with shape (264,), "
            f"got {base_control.shape}"
        )
    wave = _unpack_wave(balanced_selected.get("coefficients_original"))
    geometry = frozen["inputs"]["wave"]
    center = tuple(float(value) for value in geometry["center"])
    widths = tuple(float(value) for value in geometry["widths"])
    carrier = np.asarray(geometry["carrier"], dtype=float)
    tau0 = float(frozen["inputs"]["mean"]["tau"])
    k0 = float(frozen["inputs"]["mean"]["k"])
    tau1 = 0.5 * 2.0 ** (-(k0 + DELTA_K))
    hspace = float(frozen["timesteps"]["hspace"])
    htime = float(frozen["timesteps"]["htime"])
    return (
        balanced,
        frozen,
        constrained,
        base_control,
        acceleration_coefficients,
        wave,
        center,
        widths,
        carrier,
        tau0,
        tau1,
        k0,
        hspace,
        htime,
    )


def run(output_path=OUTPUT_PATH):
    started = time.perf_counter()
    output_path = Path(output_path)
    (
        grid_source,
        grid_cache,
        points,
        weights,
        cached_residual,
    ) = _load_holdout()
    (
        balanced,
        frozen,
        constrained,
        base_control,
        acceleration_coefficients,
        wave,
        center,
        widths,
        carrier,
        tau0,
        tau1,
        k0,
        hspace,
        htime,
    ) = _load_inputs()
    dt = tau0 - tau1
    grid_hash = _array_sha256(points, weights)
    report = {
        "status": "initializing",
        "accepted": False,
        "pde_validated": False,
        "trajectory_validated": False,
        "scale_recursion_established": False,
        "constraints_maintained": False,
        "scope": (
            "Independent actual Cartesian finite-difference momentum replay on the "
            "18,720-point new_frozen_cache holdout grid. Both the balanced parent "
            "and the selected constrained acceleration candidate are evaluated "
            "directly; no cached design or predicted residual supplies the metrics. "
            "This is a spatial holdout check only, with no PDE, trajectory, or "
            "recursion acceptance."
        ),
        "sources": {
            "holdout_grid_report": {
                "path": GRID_PATH.name,
                "sha256": _sha256(GRID_PATH),
            },
            "balanced_parent": {
                "path": BALANCED_PATH.name,
                "sha256": _sha256(BALANCED_PATH),
            },
            "frozen_geometry": {
                "path": FROZEN_PATH.name,
                "sha256": _sha256(FROZEN_PATH),
            },
            "selected_acceleration": {
                "path": ACCELERATION_PATH.name,
                "sha256": _sha256(ACCELERATION_PATH),
            },
            "momentum_helper": {
                "path": PROJECTION_PATH.name,
                "sha256": _sha256(PROJECTION_PATH),
            },
        },
        "grid": {
            "origin": "full_wave_dense_tangent.json:new_frozen_cache",
            "point_count": int(len(points)),
            "point_shape": list(points.shape),
            "weight_shape": list(weights.shape),
            "weighting": grid_cache.get("weighting"),
            "array_sha256": grid_hash,
            "cached_residual_is_reference_only": True,
            "cached_reference_metric": _metric(cached_residual, weights),
        },
        "inputs": {
            "delta_k": DELTA_K,
            "k0": k0,
            "tau0": tau0,
            "tau1": tau1,
            "physical_time_increment": dt,
            "hspace": hspace,
            "htime": htime,
            "viscosity": float(frozen["inputs"]["mean"]["nu"]),
            "base_control_count": int(base_control.size),
            "acceleration_control_count": int(acceleration_coefficients.size),
            "center": list(center),
            "widths": list(widths),
            "carrier": carrier.tolist(),
            "replay_method": "affine_momentum.jets five-point Cartesian FD then affine_momentum.momentum",
            "cached_design_used_for_metrics": False,
        },
        "candidate_provenance": {
            "balanced_parent_selected_reason": balanced["selected"].get("reason"),
            "acceleration_selected_reason": constrained["selected"].get("reason"),
            "acceleration_reference_state_unchanged": constrained["selected"].get(
                "reference_state_unchanged"
            ),
            "acceleration_linearized_gates_only": True,
        },
    }
    _save(output_path, report)
    print(json.dumps({"stage": "sources_loaded", "point_count": len(points), "grid_sha256": grid_hash}), flush=True)

    mean, mean_report = load_saved_field()
    replacements = install_in_field(mean)
    base_field = MixedAffineField(mean, center, widths, carrier, wave, base_control, tau0)
    corrected_data = _unpack_acceleration_control(acceleration_coefficients)
    corrected_field = AccelerationCorrectionField(
        base_field,
        center,
        widths,
        carrier,
        {mode: value[0] for mode, value in corrected_data.items()},
        {mode: value[1] for mode, value in corrected_data.items()},
        tau0,
    )
    report["mean_loader"] = {
        "status": mean_report.get("status", "loaded"),
        "grouped_backend_replacements": int(replacements),
    }
    report["status"] = "fields_ready"
    _save(output_path, report)
    print(json.dumps({"stage": "fields_ready", "grouped_backend_replacements": replacements}), flush=True)

    print(json.dumps({"stage": "balanced_parent_replay_started", "point_count": len(points)}), flush=True)
    base_jets = jets(base_field, points, tau1, hspace, htime)
    base_residual = momentum(base_jets)
    report["balanced_parent"] = {
        "method": "actual MixedAffineField plus affine_momentum.jets five-point Cartesian FD",
        "metric": _metric(base_residual, weights),
        "max_residual": _max_detail(base_residual, points, weights),
    }
    report["status"] = "balanced_parent_replayed"
    _save(output_path, report)
    print(json.dumps({"stage": "balanced_parent_replayed", "l2": report["balanced_parent"]["metric"]["momentum_volume_L2"], "max": report["balanced_parent"]["metric"]["momentum_max"]}), flush=True)

    print(json.dumps({"stage": "corrected_candidate_replay_started", "point_count": len(points)}), flush=True)
    corrected_jets = jets(corrected_field, points, tau1, hspace, htime)
    corrected_residual = momentum(corrected_jets)
    difference = corrected_residual - base_residual
    report["corrected_acceleration"] = {
        "method": "actual AccelerationCorrectionField plus affine_momentum.jets five-point Cartesian FD",
        "metric": _metric(corrected_residual, weights),
        "max_residual": _max_detail(corrected_residual, points, weights),
        "difference_from_balanced_parent": _metric(difference, weights),
        "difference_max_detail": _max_detail(difference, points, weights),
        "quadratic_remainder_not_inferred": True,
    }
    report["comparison"] = {
        "l2_improvement": float(
            report["balanced_parent"]["metric"]["momentum_volume_L2"]
            - report["corrected_acceleration"]["metric"]["momentum_volume_L2"]
        ),
        "l2_relative_change": float(
            (
                report["corrected_acceleration"]["metric"]["momentum_volume_L2"]
                - report["balanced_parent"]["metric"]["momentum_volume_L2"]
            )
            / max(report["balanced_parent"]["metric"]["momentum_volume_L2"], 1.0e-300)
        ),
        "peak_change": float(
            report["corrected_acceleration"]["metric"]["momentum_max"]
            - report["balanced_parent"]["metric"]["momentum_max"]
        ),
    }
    report["status"] = "completed"
    report["elapsed_seconds"] = float(time.perf_counter() - started)
    _save(output_path, report)
    print(json.dumps({
        "stage": "completed",
        "balanced_l2": report["balanced_parent"]["metric"]["momentum_volume_L2"],
        "corrected_l2": report["corrected_acceleration"]["metric"]["momentum_volume_L2"],
        "corrected_max": report["corrected_acceleration"]["metric"]["momentum_max"],
        "elapsed_seconds": report["elapsed_seconds"],
    }), flush=True)
    return report


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT_PATH)
    args = parser.parse_args()
    run(args.output)
