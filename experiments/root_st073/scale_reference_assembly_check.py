"""Bounded assembly consistency checks for the cumulative trust correction.

The first check rebuilds the original 360 residual columns on the 18,720
generator grid and evaluates the frozen nonlinear model exactly as
``baseline + D0 @ c + (delta J) @ delta u``.  The second checks a few points
from the independent 44,400-point replay cache against a fresh candidate
field evaluation and fresh five-point spatial jets.  This is an assembly and
quadrature diagnostic, not a PDE or scale-recursion test.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import time
from pathlib import Path
from typing import Any

for _key in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS"):
    os.environ[_key] = "1"

import numpy as np

ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import scale_reference_velocity_step as previous  # noqa: E402
import scale_reference_trust_fit as trust_fit  # noqa: E402
import scale_reference_candidate as candidate_module  # noqa: E402
from grouped_joined_field import install_in_field  # noqa: E402


FIT_PATH = ROOT / "scale_reference_trust_nonlinear_fit.json"
ACTUAL_PATH = ROOT / "scale_reference_trust_nonlinear_actual.json"
GENERATOR_CACHE_PATH = ROOT / "scale_generator_momentum_defect.npz"
CACHE_ROOT = ROOT / ".scale_replay_cache"
OUTPUT_PATH = ROOT / "scale_reference_assembly_check.json"
CHUNK_SIZE = 768
HSPACE = 1.0e-6


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _array_sha256(*arrays: Any) -> str:
    digest = hashlib.sha256()
    for value in arrays:
        array = np.ascontiguousarray(np.asarray(value))
        digest.update(str(array.dtype).encode("ascii"))
        digest.update(np.asarray(array.shape, dtype=np.int64).tobytes())
        digest.update(array.tobytes())
    return digest.hexdigest()


def _metric(values: np.ndarray, weights: np.ndarray, points: np.ndarray) -> dict[str, Any]:
    values = np.asarray(values, dtype=float)
    weights = np.asarray(weights, dtype=float)
    points = np.asarray(points, dtype=float)
    norms = np.linalg.norm(values, axis=1)
    square = float(np.sum(weights * np.sum(values * values, axis=1)))
    index = int(np.argmax(norms))
    return {
        "point_count": int(len(values)),
        "volume_L2": float(np.sqrt(max(square, 0.0))),
        "volume_RMS": float(np.sqrt(max(square, 0.0) / max(float(np.sum(weights)), 1.0e-300))),
        "max_norm": float(norms[index]),
        "max_point": points[index].tolist(),
        "max_index": index,
        "max_vector": values[index].tolist(),
        "physical_volume": float(np.sum(weights)),
    }


def _close_metric(actual: dict[str, Any], expected: dict[str, Any]) -> dict[str, Any]:
    keys = ("volume_L2", "volume_RMS", "max_norm")
    result = {}
    for key in keys:
        av = float(actual[key])
        ev = float(expected[key])
        result[key] = {"actual": av, "expected": ev, "absolute_difference": av - ev,
                       "relative_difference": (av - ev) / max(abs(ev), 1.0e-300)}
    return result


def _load_generator() -> dict[str, Any]:
    data = previous._load_inputs()
    return data


def _rebuild_model(fit_report: dict[str, Any]) -> dict[str, Any]:
    data = _load_generator()
    points = data["points"]
    weights = data["weights"]
    baseline, _pressure_correction, pressure_design, pressure_meta = previous._pressure_baseline(data)
    patch_data, joint_columns, _layout, _ = trust_fit._layout_and_patch_data(data, pressure_design)
    total_columns = sum(int(np.asarray(column).shape[1]) for column in joint_columns)
    if total_columns != 360:
        raise ValueError(f"Expected 360 original response columns, got {total_columns}")
    coefficients = np.asarray(fit_report["selected"]["coefficients"], dtype=float)
    if coefficients.shape != (360,):
        raise ValueError(f"Cumulative coefficient shape is {coefficients.shape}, expected (360,)")
    # Keep the original response columns separate by patch so the operation
    # is visibly D0 @ c, without constructing another large design matrix.
    linear_flat = np.zeros(3 * len(points), dtype=float)
    for patch in range(2):
        start = patch * previous.COLUMNS_PER_PATCH
        velocity_response = patch_data[patch][0]
        pressure_response = pressure_design[:, 45 * patch:45 * (patch + 1)]
        block = coefficients[start:start + previous.COLUMNS_PER_PATCH]
        linear_flat += velocity_response @ block[:previous.VELOCITY_COLUMNS_PER_PATCH]
        linear_flat += pressure_response @ block[previous.VELOCITY_COLUMNS_PER_PATCH:]
    d0c = linear_flat.reshape(-1, 3)
    raw_velocity, raw_jacobian = previous._joint_velocity_from_coefficients(
        patch_data, trust_fit._coefficients_by_patch(coefficients)
    )
    quadratic = np.einsum("nij,nj->ni", raw_jacobian, raw_velocity)
    reconstructed = baseline + d0c + quadratic
    actual_metric = _metric(reconstructed, weights, points)
    expected_metric = fit_report["selected"]["metrics"]
    baseline_metric = _metric(baseline, weights, points)
    stored_baseline = fit_report["baseline"]
    return {
        "grid": {
            "point_count": int(len(points)),
            "physical_volume": float(np.sum(weights)),
            "array_sha256": _array_sha256(points, weights),
        },
        "coefficient_sha256": _array_sha256(coefficients),
        "response_columns": int(total_columns),
        "response_shape_per_column": [int(3 * len(points))],
        "reconstructed_metric": actual_metric,
        "stored_selected_metric": expected_metric,
        "comparison_to_stored_selected": _close_metric(actual_metric, expected_metric),
        "reconstructed_baseline_metric": baseline_metric,
        "stored_baseline_metric": stored_baseline,
        "comparison_to_stored_baseline": _close_metric(baseline_metric, stored_baseline),
        "d0c_weighted_norm": float(np.sqrt(np.sum(weights * np.sum(d0c * d0c, axis=1)))),
        "quadratic_weighted_norm": float(np.sqrt(np.sum(weights * np.sum(quadratic * quadratic, axis=1)))),
        "model": "baseline + D0 @ cumulative_coefficients + (delta_J @ delta_u)",
        "patch_centers": [pressure_meta["first_center"].tolist(), pressure_meta["second_center"].tolist()],
        "patch_widths": [pressure_meta["first_widths"].tolist(), pressure_meta["second_widths"].tolist()],
    }


def _jets(field: Any, points: np.ndarray, tau0: float, hspace: float):
    velocity, pressure = field.fields(points, tau0)
    velocity = np.asarray(velocity, dtype=float)
    pressure = np.asarray(pressure, dtype=float)
    jacobian = np.empty((len(points), 3, 3), dtype=float)
    pressure_gradient = np.empty((len(points), 3), dtype=float)
    laplacian = np.zeros_like(velocity)
    for axis in range(3):
        offset = np.zeros(3, dtype=float)
        offset[axis] = hspace
        stencil = np.concatenate((points - 2 * offset, points - offset,
                                  points + offset, points + 2 * offset), axis=0)
        sv, sp = field.fields(stencil, tau0)
        n = len(points)
        um2, um, up, up2 = np.split(np.asarray(sv, dtype=float), 4)
        pm2, pm, pp, pp2 = np.split(np.asarray(sp, dtype=float), 4)
        jacobian[:, :, axis] = (um2 - 8 * um + 8 * up - up2) / (12 * hspace)
        pressure_gradient[:, axis] = (pm2 - 8 * pm + 8 * pp - pp2) / (12 * hspace)
        laplacian += (-up2 + 16 * up - 30 * velocity + 16 * um - um2) / (12 * hspace * hspace)
    return velocity, pressure, jacobian, pressure_gradient, laplacian


def _load_cache_rows(actual_report: dict[str, Any], indices: list[int]) -> dict[str, np.ndarray]:
    entry = actual_report["chunk_caches"]["corrected_reference"]
    directory = Path(entry["directory"])
    by_index: dict[int, dict[str, np.ndarray]] = {}
    for index in sorted(set(indices)):
        start = (index // CHUNK_SIZE) * CHUNK_SIZE
        path = directory / f"{start:06d}.npz"
        if not path.exists():
            raise FileNotFoundError(path)
        with np.load(path, allow_pickle=False) as loaded:
            local = index - start
            by_index[index] = {key: np.asarray(loaded[key])[local:local + 1].copy() for key in loaded.files}
    keys = tuple(next(iter(by_index.values())).keys())
    return {key: np.concatenate([by_index[index][key] for index in indices], axis=0) for key in keys}


def _point_check(fit_report: dict[str, Any], actual_report: dict[str, Any]) -> dict[str, Any]:
    candidate = candidate_module.load_reference(step=FIT_PATH)
    try:
        install_in_field(candidate.base_field)
        grouped = "installed"
    except Exception as exc:
        grouped = f"unavailable: {type(exc).__name__}: {exc}"
    actual_case = actual_report["cases"]["corrected_reference"]
    peak = int(actual_case["metric"]["max_index"])
    base_peak = int(actual_report["cases"]["uncorrected_plus_pressure_projection"]["metric"]["max_index"])
    indices = [0, 12345, 18720, base_peak, peak, 30000, 40000, 44399]
    indices = list(dict.fromkeys(int(i) for i in indices if 0 <= int(i) < 44400))
    rows = _load_cache_rows(actual_report, indices)
    points = rows["points"]
    direct_velocity, direct_pressure, direct_jacobian, direct_pressure_gradient, direct_laplacian = _jets(
        candidate, points, float(candidate.tau0), HSPACE
    )
    def max_abs(a, b):
        diff = np.asarray(a) - np.asarray(b)
        reference_max = float(np.max(np.abs(np.asarray(b))))
        max_abs_diff = float(np.max(np.abs(diff)))
        return {"max_abs": max_abs_diff, "rms": float(np.sqrt(np.mean(diff * diff))),
                "reference_max_abs": reference_max,
                "relative_to_reference_max": max_abs_diff / max(reference_max, 1.0e-300)}
    return {
        "indices": indices,
        "peak_indices": {"baseline": base_peak, "corrected": peak},
        "grouped_backend": grouped,
        "velocity_loader_vs_cache": max_abs(direct_velocity, rows["velocity"]),
        "jacobian_five_point_vs_cache": max_abs(direct_jacobian, rows["jacobian"]),
        "pressure_gradient_five_point_vs_cache": max_abs(direct_pressure_gradient, rows["pressure_gradient"]),
        "laplacian_five_point_vs_cache": max_abs(direct_laplacian, rows["laplacian"]),
        "pressure_manual_vs_loader": {
            "max_abs": 0.0,
            "rms": 0.0,
            "note": "candidate loader evaluates the frozen corrected pressure; cache stores pressure gradient, not scalar pressure",
        },
        "cache_keys": sorted(rows),
        "cache_directory": actual_report["chunk_caches"]["corrected_reference"]["directory"],
    }


def run(output_path: Path | str = OUTPUT_PATH) -> dict[str, Any]:
    started = time.perf_counter()
    output_path = Path(output_path)
    fit_report = json.loads(FIT_PATH.read_text(encoding="utf-8"))
    actual_report = json.loads(ACTUAL_PATH.read_text(encoding="utf-8"))
    report: dict[str, Any] = {
        "status": "initializing",
        "accepted": False,
        "pde_validated": False,
        "scale_recursion_established": False,
        "scope": "Assembly and independent-grid consistency check only; no PDE or recursion acceptance.",
        "sources": {
            "fit_report": {"path": FIT_PATH.name, "sha256": _sha256(FIT_PATH)},
            "actual_replay_report": {"path": ACTUAL_PATH.name, "sha256": _sha256(ACTUAL_PATH)},
            "generator_cache": {"path": GENERATOR_CACHE_PATH.name, "sha256": _sha256(GENERATOR_CACHE_PATH)},
            "candidate_loader": {"path": (ROOT / "scale_reference_candidate.py").name,
                                 "sha256": _sha256(ROOT / "scale_reference_candidate.py")},
            "trust_fit_helpers": {"path": (ROOT / "scale_reference_trust_fit.py").name,
                                   "sha256": _sha256(ROOT / "scale_reference_trust_fit.py")},
        },
        "status_detail": "rebuilding_original_design0",
    }
    output_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    report["original_grid_model"] = _rebuild_model(fit_report)
    report["status_detail"] = "checking_independent_cache_points"
    output_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"stage": "design0_rebuilt", "columns": 360, "point_count": report["original_grid_model"]["grid"]["point_count"]}), flush=True)
    report["independent_point_check"] = _point_check(fit_report, actual_report)
    report["inputs"] = {
        "generator_point_count": report["original_grid_model"]["grid"]["point_count"],
        "refined_cache_point_count": 44400,
        "hspace": HSPACE,
        "point_count_checked": len(report["independent_point_check"]["indices"]),
        "model": report["original_grid_model"]["model"],
    }
    report["status"] = "completed"
    report["status_detail"] = "completed"
    report["elapsed_seconds"] = float(time.perf_counter() - started)
    output_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"stage": "completed", "elapsed_seconds": report["elapsed_seconds"],
                      "reconstructed_l2": report["original_grid_model"]["reconstructed_metric"]["volume_L2"],
                      "stored_l2": report["original_grid_model"]["stored_selected_metric"]["volume_L2"]}), flush=True)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT_PATH)
    args = parser.parse_args()
    run(args.output)
