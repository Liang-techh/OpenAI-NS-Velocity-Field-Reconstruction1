"""Balanced refined-grid reference correction diagnostic.

This companion to ``scale_reference_refined_fit`` assembles the same 360
initial-correction columns once, then solves four trust-region directions with
different peak weights.  Each direction is tested with the exact quadratic
convection remainder on the frozen 44,400-point replay cache.  A candidate is
reported only when both the volume L2 and the pointwise maximum decrease.
"""

from __future__ import annotations

import argparse
import gc
import json
import time
from pathlib import Path
from typing import Any

import numpy as np

import scale_reference_refined_fit as refined
import scale_reference_trust_fit as trust_fit


ROOT = Path(__file__).resolve().parent
OUTPUT_PATH = ROOT / "scale_reference_refined_balanced.json"
PEAK_MULTIPLIERS = (0.0, 0.002, 0.01, 0.05)
PEAK_COUNT = refined.PEAK_COUNT
RIDGE_RELATIVE = refined.RIDGE_RELATIVE
CHUNK_SIZE = refined.CHUNK_SIZE
TOTAL_COLUMNS = refined.TOTAL_COLUMNS
COLUMNS_PER_PATCH = refined.COLUMNS_PER_PATCH
VELOCITY_COLUMNS_PER_PATCH = refined.VELOCITY_COLUMNS_PER_PATCH
PRESSURE_COLUMNS_PER_PATCH = refined.PRESSURE_COLUMNS_PER_PATCH


def _save(path: Path, report: dict[str, Any]) -> None:
    path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")


def _candidate_blocks(coefficients: np.ndarray) -> dict[str, Any]:
    return {
        "velocity_coefficients_patch1": coefficients[:135].tolist(),
        "pressure_coefficients_patch1": coefficients[135:180].tolist(),
        "velocity_coefficients_patch2": coefficients[180:315].tolist(),
        "pressure_coefficients_patch2": coefficients[315:360].tolist(),
    }


def _solve_direction(design: np.ndarray, values: np.ndarray, residual: np.ndarray,
                     weights: np.ndarray, peak_indices: np.ndarray,
                     trust_bound: float, peak_multiplier: float) -> dict[str, Any]:
    """Solve one weighted least-squares direction under a velocity trust ball."""

    row_weight = np.repeat(np.sqrt(weights), 3)
    peak_weight = float(peak_multiplier) * np.sqrt(float(np.sum(weights)))
    if peak_multiplier > 0.0:
        active_rows = np.concatenate([3 * peak_indices + component for component in range(3)])
        row_weight[active_rows] = np.maximum(row_weight[active_rows], peak_weight)
    weighted_design = design * row_weight[:, None]
    rhs = -residual.reshape(-1) * row_weight
    velocity_indices = np.r_[
        np.arange(VELOCITY_COLUMNS_PER_PATCH),
        np.arange(COLUMNS_PER_PATCH, COLUMNS_PER_PATCH + VELOCITY_COLUMNS_PER_PATCH),
    ]
    pressure_indices = np.r_[
        np.arange(VELOCITY_COLUMNS_PER_PATCH, COLUMNS_PER_PATCH),
        np.arange(COLUMNS_PER_PATCH + VELOCITY_COLUMNS_PER_PATCH, TOTAL_COLUMNS),
    ]
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
    _trust_u, trust_s, trust_vh = np.linalg.svd(trust_r, full_matrices=False)
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
    scalar_solver = {
        key: float(value) if isinstance(value, (float, np.floating)) else int(value)
        for key, value in solver.items()
        if key in ("multiplier", "unconstrained_norm", "constrained_norm", "iterations")
    }
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
        "peak_multiplier": float(peak_multiplier),
        "peak_weight": float(peak_weight),
        "active_peak_indices": peak_indices.tolist(),
        "solver": scalar_solver,
    }


def _full_basis_delta(points: np.ndarray, coefficients: np.ndarray,
                      centers: list[np.ndarray], widths: list[np.ndarray],
                      carrier: np.ndarray, nu: float) -> tuple[np.ndarray, np.ndarray]:
    delta_u = np.empty((len(points), 3), dtype=float)
    delta_j = np.empty((len(points), 3, 3), dtype=float)
    for begin in range(0, len(points), CHUNK_SIZE):
        end = min(begin + CHUNK_SIZE, len(points))
        delta_u[begin:end], delta_j[begin:end] = refined._basis_delta(
            points[begin:end], coefficients, centers, widths, carrier, nu
        )
    return delta_u, delta_j


def _directional_derivatives(residual: np.ndarray, linear: np.ndarray,
                             weights: np.ndarray) -> dict[str, Any]:
    norms = np.linalg.norm(residual, axis=1)
    l2 = float(np.sqrt(np.sum(weights * np.sum(residual * residual, axis=1))))
    l2_derivative = float(np.sum(weights * np.sum(residual * linear, axis=1)) / max(l2, 1.0e-300))
    max_norm = float(np.max(norms))
    active = np.flatnonzero(norms >= max_norm * (1.0 - 1.0e-12))
    point_derivative = np.sum(residual[active] * linear[active], axis=1) / np.maximum(norms[active], 1.0e-300)
    return {
        "volume_L2_at_zero": l2,
        "volume_L2_directional_derivative": l2_derivative,
        "active_max_count": int(len(active)),
        "active_max_indices": active.tolist(),
        "active_max_directional_derivative_min": float(np.min(point_derivative)),
        "active_max_directional_derivative_max": float(np.max(point_derivative)),
        "active_max_directional_derivatives": point_derivative.tolist(),
    }


def run(output_path: Path | str = OUTPUT_PATH) -> dict[str, Any]:
    started = time.perf_counter()
    output_path = Path(output_path)
    fit_report = json.loads(refined.FIT_PATH.read_text(encoding="utf-8"))
    actual_report = json.loads(refined.ACTUAL_PATH.read_text(encoding="utf-8"))
    points, weights, base_velocity, base_jacobian, _base_pressure_gradient, base_laplacian, _base_mean_swirl, base_residual = refined._load_replay_cache()
    frozen = fit_report["inputs"]
    tau0 = float(frozen["tau0"])
    h = float(frozen["similarity_exponent_h"])
    nu = float(frozen["viscosity"])
    carrier_value = frozen.get("carrier")
    if carrier_value is None:
        helper = json.loads((ROOT / "scale_reference_velocity_step.json").read_text(encoding="utf-8"))
        carrier_value = helper["inputs"]["carrier"]
    carrier = np.asarray(carrier_value, dtype=float)
    centers = [np.asarray(fit_report["patches"][name]["center"], dtype=float) for name in ("first", "second")]
    widths = [np.asarray(fit_report["patches"][name]["widths"], dtype=float) for name in ("first", "second")]
    old_coefficients = np.asarray(fit_report["selected"]["coefficients"], dtype=float)
    if old_coefficients.shape != (TOTAL_COLUMNS,):
        raise ValueError("Frozen cumulative coefficient vector is not 360-dimensional")

    # Compute norms of the actual correction, rather than accidentally treating
    # the total field as the correction.  This is also the refined-grid norm
    # used for the additional one-percent trust ball.
    old_delta_u, _old_delta_j = _full_basis_delta(points, old_coefficients, centers, widths, carrier, nu)
    raw_reference = base_velocity - old_delta_u
    raw_reference_norm = float(np.sqrt(np.sum(weights * np.sum(raw_reference * raw_reference, axis=1))))
    old_correction_norm = float(np.sqrt(np.sum(weights * np.sum(old_delta_u * old_delta_u, axis=1))))
    trust_bound = 0.01 * raw_reference_norm
    coarse = refined._base_inputs()
    coarse_reference_norm = float(np.sqrt(np.sum(coarse["weights"] * np.sum(coarse["velocity"] * coarse["velocity"], axis=1))))
    base_metric = refined._metric(base_residual, weights, points)
    residual_norms = np.linalg.norm(base_residual, axis=1)
    peak_indices = np.argsort(residual_norms)[-PEAK_COUNT:][::-1]

    report: dict[str, Any] = {
        "status": "initializing",
        "accepted": False,
        "pde_validated": False,
        "scale_recursion_established": False,
        "constraints_maintained": False,
        "scope": "Balanced peak/L2 refined-grid reference correction diagnostic; no PDE, trajectory, or recursion acceptance.",
        "sources": {
            "fit_report": {"path": refined.FIT_PATH.name, "sha256": refined._sha256(refined.FIT_PATH)},
            "actual_replay_report": {"path": refined.ACTUAL_PATH.name, "sha256": refined._sha256(refined.ACTUAL_PATH)},
            "generator_cache": {"path": refined.GENERATOR_CACHE_PATH.name, "sha256": refined._sha256(refined.GENERATOR_CACHE_PATH)},
            "refined_fit_helpers": {"path": (ROOT / "scale_reference_refined_fit.py").name, "sha256": refined._sha256(ROOT / "scale_reference_refined_fit.py")},
            "velocity_step_helpers": {"path": (ROOT / "scale_reference_velocity_step.py").name, "sha256": refined._sha256(ROOT / "scale_reference_velocity_step.py")},
            "trust_helpers": {"path": (ROOT / "scale_reference_trust_fit.py").name, "sha256": refined._sha256(ROOT / "scale_reference_trust_fit.py")},
            "pressure_design": {"path": (ROOT / "scale_generator_pressure_projection.py").name, "sha256": refined._sha256(ROOT / "scale_generator_pressure_projection.py")},
        },
        "inputs": {
            "point_count": int(len(points)),
            "physical_volume": float(np.sum(weights)),
            "grid_array_sha256": refined._array_sha256(points, weights),
            "tau0": tau0,
            "similarity_exponent_h": h,
            "viscosity": nu,
            "carrier": carrier.tolist(),
            "degree": refined.DEGREE,
            "modes": list(refined.MODES),
            "patch_count": 2,
            "velocity_columns_per_patch": VELOCITY_COLUMNS_PER_PATCH,
            "pressure_columns_per_patch": PRESSURE_COLUMNS_PER_PATCH,
            "joint_column_count": TOTAL_COLUMNS,
            "peak_multipliers": list(PEAK_MULTIPLIERS),
            "peak_count": PEAK_COUNT,
            "additional_velocity_trust_fraction": 0.01,
            "original_reference_velocity_norm_refined": raw_reference_norm,
            "original_reference_velocity_norm_coarse": coarse_reference_norm,
            "additional_velocity_trust_bound_refined": trust_bound,
            "old_correction_norm_refined": old_correction_norm,
            "cache_chunk_count": len(list(Path(actual_report["chunk_caches"]["corrected_reference"]["directory"]).glob("*.npz"))),
            "correction_kind": "initial reference curl-potential velocity increment around refined corrected state",
        },
        "patches": fit_report.get("patches"),
        "baseline": base_metric,
        "status_detail": "building_balanced_refined_columns",
    }
    # The shared helper's decoder check is copied into this companion report
    # so a selected balanced direction carries the same sign-convention
    # evidence as the primary refined fit.
    if refined.OUTPUT_PATH.exists():
        try:
            helper_report = json.loads(refined.OUTPUT_PATH.read_text(encoding="utf-8"))
            report["consistency_check"] = helper_report.get("consistency_check")
        except Exception as exc:  # pragma: no cover - diagnostic provenance only
            report["consistency_check_error"] = repr(exc)
    _save(output_path, report)
    print(json.dumps({"stage": "cache_loaded", "point_count": int(len(points)), "baseline_l2": base_metric["volume_L2"], "baseline_max": base_metric["max_norm"], "trust_bound": trust_bound}), flush=True)

    design, values = refined._build_columns(points, base_velocity, base_jacobian, base_laplacian, carrier, centers, widths, nu, h, tau0)
    report["design"] = {"shape": [int(v) for v in design.shape], "values_shape": [int(v) for v in values.shape]}
    report["status_detail"] = "solving_balanced_directions"
    _save(output_path, report)

    directions: list[dict[str, Any]] = []
    selected: dict[str, Any] | None = None
    alphas = [1.0 / (2.0 ** i) for i in range(13)]
    for direction_index, peak_multiplier in enumerate(PEAK_MULTIPLIERS):
        print(json.dumps({"stage": "solve_direction", "direction_index": direction_index, "peak_multiplier": peak_multiplier}), flush=True)
        solved = _solve_direction(design, values, base_residual, weights, peak_indices, trust_bound, peak_multiplier)
        increment = solved["coefficients"]
        delta_u = (values @ solved["velocity_coefficients"]).reshape(-1, 3)
        linear = (design @ increment).reshape(-1, 3)
        quadratic = np.empty_like(base_residual)
        for begin in range(0, len(points), CHUNK_SIZE):
            end = min(begin + CHUNK_SIZE, len(points))
            du, dj = refined._basis_delta(points[begin:end], increment, centers, widths, carrier, nu)
            quadratic[begin:end] = np.einsum("nij,nj->ni", dj, du)
        derivatives = _directional_derivatives(base_residual, linear, weights)
        trials = []
        direction_selected: dict[str, Any] | None = None
        for alpha in alphas:
            metric = refined._metrics_for_trial(base_residual, linear, quadratic, alpha, points, weights)
            velocity_metric = refined._metric(alpha * delta_u, weights, points)
            row = {
                "alpha": alpha,
                "metric": metric,
                "velocity_metric": velocity_metric,
                "trust_region_pass": bool(velocity_metric["volume_L2"] <= trust_bound * (1.0 + 1.0e-10)),
                "improves_volume_L2": bool(metric["volume_L2"] < base_metric["volume_L2"]),
                "improves_max": bool(metric["max_norm"] < base_metric["max_norm"]),
            }
            trials.append(row)
            if direction_selected is None and row["trust_region_pass"] and row["improves_volume_L2"] and row["improves_max"]:
                direction_selected = row
        old_correction_plus = old_delta_u + (float(direction_selected["alpha"]) if direction_selected else 0.0) * delta_u
        best_l2_trial = min(trials, key=lambda row: row["metric"]["volume_L2"])
        best_l2_alpha = float(best_l2_trial["alpha"])
        best_l2_correction = old_delta_u + best_l2_alpha * delta_u
        cumulative_best_coefficients = old_coefficients + best_l2_alpha * increment
        entry: dict[str, Any] = {
            "peak_multiplier": float(peak_multiplier),
            "solver": {key: value for key, value in solved.items() if key not in ("coefficients", "velocity_coefficients", "pressure_coefficients")},
            "increment_coefficients": increment.tolist(),
            "increment_blocks": _candidate_blocks(increment),
            "trials": trials,
            "directional_derivatives": derivatives,
            "best_l2_trial": best_l2_trial,
            "best_l2_summed_candidate_coefficients": cumulative_best_coefficients.tolist(),
            "best_l2_summed_candidate_blocks": _candidate_blocks(cumulative_best_coefficients),
            "best_l2_cumulative_correction_norm_refined": float(np.sqrt(np.sum(weights * np.sum(best_l2_correction * best_l2_correction, axis=1)))),
            "direction_selected": direction_selected,
        }
        directions.append(entry)
        print(json.dumps({"stage": "direction_done", "peak_multiplier": peak_multiplier, "directional_L2": derivatives["volume_L2_directional_derivative"], "directional_max_min": derivatives["active_max_directional_derivative_min"], "best_l2": best_l2_trial["metric"]["volume_L2"], "best_max": min(trials, key=lambda row: row["metric"]["max_norm"])["metric"]["max_norm"], "joint": direction_selected is not None}), flush=True)
        if direction_selected is not None:
            alpha = float(direction_selected["alpha"])
            cumulative = old_coefficients + alpha * increment
            cumulative_correction = old_delta_u + alpha * delta_u
            selected = {
                "feasible": True,
                "peak_multiplier": float(peak_multiplier),
                "alpha": alpha,
                "coefficients": cumulative.tolist(),
                "increment_coefficients": (alpha * increment).tolist(),
                "summed_candidate_coefficients": cumulative.tolist(),
                **_candidate_blocks(cumulative),
                "increment_blocks": _candidate_blocks(alpha * increment),
                "metric": direction_selected["metric"],
                "velocity_metric": direction_selected["velocity_metric"],
                "cumulative_correction_norm_refined": float(np.sqrt(np.sum(weights * np.sum(cumulative_correction * cumulative_correction, axis=1)))),
                "cumulative_correction_fraction_of_original_refined": float(np.sqrt(np.sum(weights * np.sum(cumulative_correction * cumulative_correction, axis=1))) / max(raw_reference_norm, 1.0e-300)),
            }
            break
        del solved, delta_u, linear, quadratic
        gc.collect()

    report["directions"] = directions
    report["selected"] = selected
    report["control_semantics"] = "initial velocity correction"
    report["sources"].update({
        "two_patch_candidate": {"path": "localized_drift_1800_fit.json", "sha256": refined._sha256(ROOT / "localized_drift_1800_fit.json")},
        "geometry": {"path": "full_wave_frozen_cache.json", "sha256": refined._sha256(ROOT / "full_wave_frozen_cache.json")},
        "pressure_report": {"path": "scale_generator_pressure_projection.json", "sha256": refined._sha256(ROOT / "scale_generator_pressure_projection.json")},
    })
    report["status"] = "completed"
    report["status_detail"] = "completed"
    report["accepted"] = False
    report["elapsed_seconds"] = float(time.perf_counter() - started)
    report["reason"] = (
        "First balanced direction with a trial reducing both refined-grid L2 and maximum is recorded; "
        "this remains a diagnostic candidate without PDE, trajectory, or recursion validation."
        if selected is not None else
        "No tested balanced direction had a trust-feasible line-search trial reducing both refined-grid L2 and maximum; directional derivatives and rejected coefficients are retained."
    )
    _save(output_path, report)
    print(json.dumps({"stage": "completed", "selected": selected is not None, "directions": len(directions), "baseline_l2": base_metric["volume_L2"], "selected_l2": None if selected is None else selected["metric"]["volume_L2"], "baseline_max": base_metric["max_norm"], "selected_max": None if selected is None else selected["metric"]["max_norm"], "elapsed_seconds": report["elapsed_seconds"]}), flush=True)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT_PATH)
    args = parser.parse_args()
    run(args.output)
