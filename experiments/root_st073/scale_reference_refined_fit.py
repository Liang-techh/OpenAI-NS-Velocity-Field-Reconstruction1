"""One peak-aware trust correction on the independent 44,400-point replay grid.

The cached corrected replay supplies the current velocity, spatial jets, and
momentum residual.  Existing analytic degree-2 compact-patch columns are
rebuilt around that state.  A single weighted trust solve proposes an
increment, and the exact quadratic convection remainder is used during a
short line search.  A trial is retained only when both true volume L2 and
true pointwise maximum decrease.
"""

from __future__ import annotations

import argparse
import gc
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
from scale_generator_pressure_projection import _pressure_patch_design  # noqa: E402
from supported_fourier_analytic_jets import basis_jets  # noqa: E402


FIT_PATH = ROOT / "scale_reference_trust_nonlinear_fit.json"
ACTUAL_PATH = ROOT / "scale_reference_trust_nonlinear_actual.json"
GENERATOR_CACHE_PATH = ROOT / "scale_generator_momentum_defect.npz"
OUTPUT_PATH = ROOT / "scale_reference_refined_fit.json"
CHUNK_SIZE = 768
DEGREE = 2
MODES = (0, 1, 2)
VELOCITY_COLUMNS_PER_PATCH = 135
PRESSURE_COLUMNS_PER_PATCH = 45
COLUMNS_PER_PATCH = 180
TOTAL_COLUMNS = 360
PEAK_COUNT = 64
PEAK_WEIGHT_MULTIPLIER = 8.0
RIDGE_RELATIVE = 1.0e-10


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


def _save(path: Path, report: dict[str, Any]) -> None:
    path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")


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


def _load_replay_cache() -> tuple[np.ndarray, ...]:
    actual = json.loads(ACTUAL_PATH.read_text(encoding="utf-8"))
    entry = actual["chunk_caches"]["corrected_reference"]
    directory = Path(entry["directory"])
    starts = []
    for path in directory.glob("*.npz"):
        try:
            starts.append((int(path.stem), path))
        except ValueError:
            pass
    starts.sort()
    if not starts:
        raise FileNotFoundError(f"No replay chunks in {directory}")
    pieces = {key: [] for key in ("points", "weights", "velocity", "jacobian", "pressure_gradient", "laplacian", "mean_swirl", "residual")}
    for start, path in starts:
        with np.load(path, allow_pickle=False) as loaded:
            missing = sorted(set(pieces).difference(loaded.files))
            if missing:
                raise ValueError(f"Replay chunk {path.name} missing {missing}")
            for key in pieces:
                pieces[key].append(np.asarray(loaded[key], dtype=float))
    arrays = tuple(np.concatenate(pieces[key], axis=0) for key in pieces)
    points, weights, velocity, jacobian, pressure_gradient, laplacian, mean_swirl, residual = arrays
    count = len(points)
    if points.shape != (count, 3) or weights.shape != (count,) or velocity.shape != (count, 3):
        raise ValueError("Replay cache has invalid point/velocity shapes")
    if jacobian.shape != (count, 3, 3) or pressure_gradient.shape != (count, 3) or laplacian.shape != (count, 3):
        raise ValueError("Replay cache has invalid spatial-jet shapes")
    if mean_swirl.shape != (count,) or residual.shape != (count, 3):
        raise ValueError("Replay cache has invalid generator/residual shapes")
    if count != 44400 or np.any(weights <= 0.0) or not all(np.all(np.isfinite(a)) for a in arrays):
        raise ValueError("Replay cache is not the frozen 44,400-point finite grid")
    return arrays


def _base_inputs() -> dict[str, Any]:
    with np.load(GENERATOR_CACHE_PATH, allow_pickle=False) as loaded:
        required = {"points", "weights", "velocity"}
        if not required.issubset(loaded.files):
            raise ValueError("Generator cache lacks original reference velocity")
        points = np.asarray(loaded["points"], dtype=float)
        weights = np.asarray(loaded["weights"], dtype=float)
        velocity = np.asarray(loaded["velocity"], dtype=float)
    if points.shape != (len(points), 3) or weights.shape != (len(points),) or velocity.shape != (len(points), 3):
        raise ValueError("Invalid original generator cache shapes")
    return {"points": points, "weights": weights, "velocity": velocity}


def _build_columns(points: np.ndarray, velocity: np.ndarray, jacobian: np.ndarray,
                   laplacian: np.ndarray, carrier: np.ndarray, centers: list[np.ndarray],
                   widths: list[np.ndarray], nu: float, h: float, tau0: float) -> tuple[np.ndarray, np.ndarray]:
    """Build D0 and the physical velocity-value matrix in cache-sized chunks."""

    count = len(points)
    design = np.empty((3 * count, TOTAL_COLUMNS), dtype=float)
    values_all = np.empty((3 * count, 2 * VELOCITY_COLUMNS_PER_PATCH), dtype=float)
    for begin in range(0, count, CHUNK_SIZE):
        end = min(begin + CHUNK_SIZE, count)
        row = slice(3 * begin, 3 * end)
        for patch, (center, width) in enumerate(zip(centers, widths)):
            response, values, _jac, _lap, _swirl, _layout = previous._mode_velocity_columns(
                points[begin:end], center, width, carrier, nu, h,
                velocity[begin:end], jacobian[begin:end], laplacian[begin:end], tau0,
            )
            pressure = _pressure_patch_design(points[begin:end], center, width, carrier)
            c0 = patch * COLUMNS_PER_PATCH
            v0 = patch * VELOCITY_COLUMNS_PER_PATCH
            design[row, c0:c0 + VELOCITY_COLUMNS_PER_PATCH] = response
            design[row, c0 + VELOCITY_COLUMNS_PER_PATCH:c0 + COLUMNS_PER_PATCH] = pressure
            values_all[row, v0:v0 + VELOCITY_COLUMNS_PER_PATCH] = values.reshape(3 * (end - begin), VELOCITY_COLUMNS_PER_PATCH)
        if begin == 0 or end == count or (begin // CHUNK_SIZE) % 10 == 0:
            print(json.dumps({"stage": "basis_chunk", "processed": int(end), "total": int(count)}), flush=True)
    return design, values_all


def _solve_trust(design: np.ndarray, values: np.ndarray, residual: np.ndarray,
                 weights: np.ndarray, peak_indices: np.ndarray, trust_bound: float) -> dict[str, Any]:
    row_weight = np.repeat(np.sqrt(weights), 3)
    peak_weight = PEAK_WEIGHT_MULTIPLIER * np.sqrt(float(np.sum(weights)))
    active_rows = np.concatenate([3 * peak_indices + component for component in range(3)])
    row_weight[active_rows] = np.maximum(row_weight[active_rows], peak_weight)
    weighted_design = design * row_weight[:, None]
    rhs = -residual.reshape(-1) * row_weight
    velocity_indices = np.r_[np.arange(VELOCITY_COLUMNS_PER_PATCH),
                             np.arange(COLUMNS_PER_PATCH, COLUMNS_PER_PATCH + VELOCITY_COLUMNS_PER_PATCH)]
    pressure_indices = np.r_[np.arange(VELOCITY_COLUMNS_PER_PATCH, COLUMNS_PER_PATCH),
                             np.arange(COLUMNS_PER_PATCH + VELOCITY_COLUMNS_PER_PATCH, TOTAL_COLUMNS)]
    av = weighted_design[:, velocity_indices]
    ap = weighted_design[:, pressure_indices]
    pressure_u, pressure_s, pressure_vh = np.linalg.svd(ap, full_matrices=False)
    pressure_cutoff = RIDGE_RELATIVE * max(float(pressure_s[0]), 1.0e-300)
    pressure_rank = int(np.sum(pressure_s > pressure_cutoff))
    pressure_u = pressure_u[:, :pressure_rank]
    pressure_s = pressure_s[:pressure_rank]
    pressure_vh = pressure_vh[:pressure_rank, :]
    pressure_projection = pressure_u.T @ rhs
    projected_rhs = rhs - pressure_u @ pressure_projection
    av_projected = av - pressure_u @ (pressure_u.T @ av)
    metric_row_weight = np.repeat(np.sqrt(weights), 3)
    weighted_values = values * metric_row_weight[:, None]
    velocity_scales = np.maximum(np.linalg.norm(weighted_values, axis=0), 1.0e-30)
    normalized_metric = weighted_values / velocity_scales[None, :]
    trust_r = np.linalg.qr(normalized_metric, mode="r")
    trust_u, trust_s, trust_vh = np.linalg.svd(trust_r, full_matrices=False)
    trust_cutoff = RIDGE_RELATIVE * max(float(trust_s[0]), 1.0e-300)
    trust_rank = int(np.sum(trust_s > trust_cutoff))
    trust_v = trust_vh[:trust_rank, :].T
    trust_s_kept = trust_s[:trust_rank]
    av_scaled = av_projected / velocity_scales[None, :]
    a_tilde = (av_scaled @ trust_v) / trust_s_kept[None, :]
    solver = trust_fit._solve_multiplier(a_tilde, projected_rhs, trust_bound)
    y = solver["coordinates"]
    z = trust_v @ (y / trust_s_kept)
    velocity_coefficients = z / velocity_scales
    pressure_rhs = pressure_projection - (pressure_u.T @ av) @ velocity_coefficients
    pressure_coefficients = pressure_vh.T @ (pressure_rhs / pressure_s)
    coefficients = np.zeros(TOTAL_COLUMNS, dtype=float)
    coefficients[velocity_indices] = velocity_coefficients
    coefficients[pressure_indices] = pressure_coefficients
    return {
        "coefficients": coefficients,
        "velocity_coefficients": velocity_coefficients,
        "pressure_coefficients": pressure_coefficients,
        "trust_bound": float(trust_bound),
        "pressure_rank": pressure_rank,
        "trust_rank": trust_rank,
        "velocity_scales_min": float(np.min(velocity_scales)),
        "velocity_scales_max": float(np.max(velocity_scales)),
        "peak_count": int(len(peak_indices)),
        "peak_weight": float(peak_weight),
        "active_peak_indices": peak_indices.tolist(),
        "solver": {key: float(value) if isinstance(value, (float, np.floating)) else int(value)
                   for key, value in solver.items() if key in ("multiplier", "unconstrained_norm", "constrained_norm", "iterations")},
    }


def _basis_delta(points: np.ndarray, coefficients: np.ndarray, centers: list[np.ndarray],
                 widths: list[np.ndarray], carrier: np.ndarray, nu: float) -> tuple[np.ndarray, np.ndarray]:
    delta_u = np.zeros((len(points), 3), dtype=float)
    delta_j = np.zeros((len(points), 3, 3), dtype=float)
    for patch, (center, width) in enumerate(zip(centers, widths)):
        block = coefficients[patch * COLUMNS_PER_PATCH:patch * COLUMNS_PER_PATCH + VELOCITY_COLUMNS_PER_PATCH]
        cursor = 0
        for mode in MODES:
            mode_carrier = np.zeros(2, dtype=float) if mode == 0 else mode * carrier
            values, gradients, _viscous, _pressure, _pressure_gradient = basis_jets(
                points, center, width, mode, DEGREE, mode_carrier, nu
            )
            for component in range(3 * 9):
                if mode == 0:
                    coeff = block[cursor]
                    value = values[:, :, component].real
                    jacobian = gradients[:, :, :, component].real
                    cursor += 1
                    delta_u += coeff * value
                    delta_j += coeff * jacobian
                else:
                    coeff_real = block[cursor]
                    coeff_imag = block[cursor + 1]
                    value = values[:, :, component].real * coeff_real + values[:, :, component].imag * coeff_imag
                    jacobian = gradients[:, :, :, component].real * coeff_real + gradients[:, :, :, component].imag * coeff_imag
                    cursor += 2
                    delta_u += value
                    delta_j += jacobian
        if cursor != VELOCITY_COLUMNS_PER_PATCH:
            raise ValueError(f"velocity basis decode consumed {cursor}, expected {VELOCITY_COLUMNS_PER_PATCH}")
    return delta_u, delta_j


def _metrics_for_trial(base_residual: np.ndarray, linear: np.ndarray, quadratic: np.ndarray,
                       alpha: float, points: np.ndarray, weights: np.ndarray) -> dict[str, Any]:
    residual = base_residual + alpha * linear + alpha * alpha * quadratic
    return _metric(residual, weights, points)


def _candidate_blocks(coefficients: np.ndarray) -> dict[str, Any]:
    return {
        "velocity_coefficients_patch1": coefficients[:135].tolist(),
        "pressure_coefficients_patch1": coefficients[135:180].tolist(),
        "velocity_coefficients_patch2": coefficients[180:315].tolist(),
        "pressure_coefficients_patch2": coefficients[315:360].tolist(),
    }


def run(output_path: Path | str = OUTPUT_PATH) -> dict[str, Any]:
    started = time.perf_counter()
    output_path = Path(output_path)
    fit_report = json.loads(FIT_PATH.read_text(encoding="utf-8"))
    actual_report = json.loads(ACTUAL_PATH.read_text(encoding="utf-8"))
    points, weights, base_velocity, base_jacobian, base_pressure_gradient, base_laplacian, base_mean_swirl, base_residual = _load_replay_cache()
    frozen = fit_report["inputs"]
    tau0 = float(frozen["tau0"])
    h = float(frozen["similarity_exponent_h"])
    nu = float(frozen["viscosity"])
    # The historical nonlinear-fit report predates the explicit carrier field;
    # its helper report carries the same mode-1 carrier used to assemble the
    # selected coefficients.  Keep this diagnostic reproducible across both
    # report schemas.
    carrier_value = frozen.get("carrier")
    if carrier_value is None:
        helper_report = json.loads((ROOT / "scale_reference_velocity_step.json").read_text(encoding="utf-8"))
        carrier_value = helper_report["inputs"]["carrier"]
    carrier = np.asarray(carrier_value, dtype=float)
    centers = [np.asarray(fit_report["patches"][name]["center"], dtype=float) for name in ("first", "second")]
    widths = [np.asarray(fit_report["patches"][name]["widths"], dtype=float) for name in ("first", "second")]
    original = _base_inputs()
    original_velocity_norm = float(np.sqrt(np.sum(original["weights"] * np.sum(original["velocity"] ** 2, axis=1))))
    trust_bound = 0.01 * original_velocity_norm
    base_metric = _metric(base_residual, weights, points)
    residual_norms = np.linalg.norm(base_residual, axis=1)
    peak_indices = np.argsort(residual_norms)[-PEAK_COUNT:][::-1]
    report: dict[str, Any] = {
        "status": "initializing",
        "accepted": False,
        "pde_validated": False,
        "scale_recursion_established": False,
        "constraints_maintained": False,
        "scope": "One refined-grid peak-aware reference correction diagnostic; no PDE, trajectory, or recursion acceptance.",
        "sources": {
            "fit_report": {"path": FIT_PATH.name, "sha256": _sha256(FIT_PATH)},
            "actual_replay_report": {"path": ACTUAL_PATH.name, "sha256": _sha256(ACTUAL_PATH)},
            "generator_cache": {"path": GENERATOR_CACHE_PATH.name, "sha256": _sha256(GENERATOR_CACHE_PATH)},
            "velocity_step_helpers": {"path": (ROOT / "scale_reference_velocity_step.py").name, "sha256": _sha256(ROOT / "scale_reference_velocity_step.py")},
            "trust_helpers": {"path": (ROOT / "scale_reference_trust_fit.py").name, "sha256": _sha256(ROOT / "scale_reference_trust_fit.py")},
            "pressure_design": {"path": (ROOT / "scale_generator_pressure_projection.py").name, "sha256": _sha256(ROOT / "scale_generator_pressure_projection.py")},
        },
        "inputs": {
            "point_count": int(len(points)),
            "physical_volume": float(np.sum(weights)),
            "grid_array_sha256": _array_sha256(points, weights),
            "tau0": tau0,
            "similarity_exponent_h": h,
            "viscosity": nu,
            "carrier": carrier.tolist(),
            "degree": DEGREE,
            "modes": list(MODES),
            "patch_count": 2,
            "velocity_columns_per_patch": VELOCITY_COLUMNS_PER_PATCH,
            "pressure_columns_per_patch": PRESSURE_COLUMNS_PER_PATCH,
            "control_count_per_patch": COLUMNS_PER_PATCH,
            "joint_column_count": TOTAL_COLUMNS,
            "original_reference_velocity_norm": original_velocity_norm,
            "additional_velocity_trust_fraction": 0.01,
            "additional_velocity_trust_bound": trust_bound,
            "peak_count": PEAK_COUNT,
            "peak_weight_multiplier": PEAK_WEIGHT_MULTIPLIER,
            "cache_chunk_count": len(list(Path(actual_report["chunk_caches"]["corrected_reference"]["directory"]).glob("*.npz"))),
            "correction_kind": "initial reference curl-potential velocity increment around refined corrected state",
        },
        "baseline": base_metric,
        "status_detail": "building_refined_columns",
    }
    _save(output_path, report)
    print(json.dumps({"stage": "cache_loaded", "point_count": int(len(points)), "baseline_l2": base_metric["volume_L2"], "baseline_max": base_metric["max_norm"]}), flush=True)
    design, values = _build_columns(points, base_velocity, base_jacobian, base_laplacian, carrier, centers, widths, nu, h, tau0)
    report["status_detail"] = "solving_peak_weighted_trust"
    report["design"] = {"shape": [int(v) for v in design.shape], "values_shape": [int(v) for v in values.shape]}
    _save(output_path, report)
    solved = _solve_trust(design, values, base_residual, weights, peak_indices, trust_bound)
    increment = solved["coefficients"]
    delta_u = values @ solved["velocity_coefficients"]
    linear = design @ increment
    # Keep the checkpoint JSON scalar-only.  The coefficient arrays remain in
    # local NumPy variables and are written to the selected record only if a
    # line-search trial passes both refined-grid tests.
    solver_summary = {key: value for key, value in solved.items()
                      if key not in ("coefficients", "velocity_coefficients", "pressure_coefficients")}
    report["solver"] = solver_summary | {
        "increment_velocity_metric": _metric(delta_u.reshape(-1, 3), weights, points),
        "increment_pressure_norm": float(np.linalg.norm(solved["pressure_coefficients"])),
    }
    report["status_detail"] = "evaluating_exact_quadratic_trials"
    _save(output_path, report)
    # Rebuild only delta_u and delta_J once for the solved direction.  The
    # exact convection remainder then scales as alpha^2 in the line search.
    quadratic = np.empty_like(base_residual)
    for begin in range(0, len(points), CHUNK_SIZE):
        end = min(begin + CHUNK_SIZE, len(points))
        du, dj = _basis_delta(points[begin:end], increment, centers, widths, carrier, nu)
        quadratic[begin:end] = np.einsum("nij,nj->ni", dj, du)
    trials = []
    selected = None
    alphas = [1.0 / (2.0 ** i) for i in range(13)]
    for alpha in alphas:
        metric = _metrics_for_trial(base_residual, linear.reshape(-1, 3), quadratic, alpha, points, weights)
        velocity_metric = _metric(alpha * delta_u.reshape(-1, 3), weights, points)
        row = {
            "alpha": alpha,
            "metric": metric,
            "velocity_metric": velocity_metric,
            "trust_region_pass": bool(velocity_metric["volume_L2"] <= trust_bound * (1.0 + 1.0e-10)),
            "improves_volume_L2": bool(metric["volume_L2"] < base_metric["volume_L2"]),
            "improves_max": bool(metric["max_norm"] < base_metric["max_norm"]),
        }
        trials.append(row)
        if selected is None and row["trust_region_pass"] and row["improves_volume_L2"] and row["improves_max"]:
            selected = row
    report["trials"] = trials
    report["linear_direction"] = {
        "weighted_L2": float(np.sqrt(np.sum(weights * np.sum(linear.reshape(-1, 3) ** 2, axis=1)))),
        "quadratic_remainder": _metric(quadratic, weights, points),
    }
    old_coefficients = np.asarray(fit_report["selected"]["coefficients"], dtype=float)
    if old_coefficients.shape != (360,):
        raise ValueError("Frozen cumulative coefficient vector is not 360-dimensional")
    report["cumulative_velocity_before"] = fit_report["selected"]["velocity_correction"]
    if selected is None:
        report["selected"] = None
        report["best_l2_trial"] = min(trials, key=lambda row: row["metric"]["volume_L2"])
        report["best_max_trial"] = min(trials, key=lambda row: row["metric"]["max_norm"])
        report["reason"] = "No line-search trial improved both true refined-grid L2 and true refined-grid maximum within the additional one-percent trust bound."
    else:
        alpha = float(selected["alpha"])
        cumulative = old_coefficients + alpha * increment
        report["selected"] = {
            "feasible": True,
            "alpha": alpha,
            "increment_coefficients": (alpha * increment).tolist(),
            "summed_candidate_coefficients": cumulative.tolist(),
            **_candidate_blocks(cumulative),
            "increment_blocks": _candidate_blocks(alpha * increment),
            "metric": selected["metric"],
            "velocity_metric": selected["velocity_metric"],
            "cumulative_velocity_norm_estimate": float(np.sqrt(np.sum(weights * np.sum((base_velocity + alpha * delta_u.reshape(-1, 3)) ** 2, axis=1)))),
            "cumulative_trust_fraction_of_original": float(np.sqrt(np.sum(weights * np.sum((base_velocity + alpha * delta_u.reshape(-1, 3)) ** 2, axis=1))) / max(original_velocity_norm, 1.0e-300)),
        }
        report["reason"] = "First exact-quadratic trial improved both refined-grid L2 and maximum within the additional one-percent velocity trust bound."
    # These fields make a selected report directly consumable by the strict
    # scale_reference_candidate loader, while preserving the increment too.
    report["control_semantics"] = "initial velocity correction"
    report["sources"].update({
        "two_patch_candidate": {"path": "localized_drift_1800_fit.json", "sha256": _sha256(ROOT / "localized_drift_1800_fit.json")},
        "geometry": {"path": "full_wave_frozen_cache.json", "sha256": _sha256(ROOT / "full_wave_frozen_cache.json")},
        "pressure_report": {"path": "scale_generator_pressure_projection.json", "sha256": _sha256(ROOT / "scale_generator_pressure_projection.json")},
    })
    report["status"] = "completed"
    report["status_detail"] = "completed"
    report["elapsed_seconds"] = float(time.perf_counter() - started)
    _save(output_path, report)
    print(json.dumps({"stage": "completed", "selected": selected is not None, "baseline_l2": base_metric["volume_L2"], "selected_l2": None if selected is None else selected["metric"]["volume_L2"], "baseline_max": base_metric["max_norm"], "selected_max": None if selected is None else selected["metric"]["max_norm"], "elapsed_seconds": report["elapsed_seconds"]}), flush=True)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT_PATH)
    args = parser.parse_args()
    run(args.output)
