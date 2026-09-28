"""Bounded Gauss--Newton refinement of the frozen scale reference velocity.

The one-step reference correction is a useful warm start, but its response
columns were formed around the uncorrected velocity.  This driver rebuilds
those analytic columns around each accepted velocity/Jacobian state and uses
the exact quadratic convection remainder supplied by
``scale_generator_linearization.velocity_correction_response``.  Each update
is limited to one percent of the original velocity norm and the cumulative
correction is limited to five percent.  A trial is accepted only when the
recomputed nonlinear residual decreases.

This is a bounded scale-generator diagnostic.  It does not establish support,
moment, cone, trajectory, PDE, or recursive-scale acceptance.
"""

from __future__ import annotations

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

import scale_reference_velocity_step as one_step  # noqa: E402


GENERATOR_REPORT_PATH = ROOT / "scale_generator_momentum_defect.json"
GENERATOR_CACHE_PATH = ROOT / "scale_generator_momentum_defect.npz"
PRESSURE_REPORT_PATH = ROOT / "scale_generator_pressure_projection.json"
DRIFT_REPORT_PATH = ROOT / "localized_drift_1800_fit.json"
FROZEN_PATH = ROOT / "full_wave_frozen_cache.json"
LINEARIZATION_PATH = ROOT / "scale_generator_linearization.py"
BASIS_PATH = ROOT / "supported_fourier_analytic_jets.py"
STEP_PATH = ROOT / "scale_reference_velocity_step.json"
TRUST_PATH = ROOT / "scale_reference_trust_fit.json"
OUTPUT_PATH = ROOT / "scale_reference_nonlinear_fit.json"

DEGREE = one_step.DEGREE
MODES = one_step.MODES
VELOCITY_COLUMNS_PER_PATCH = one_step.VELOCITY_COLUMNS_PER_PATCH
PRESSURE_COLUMNS_PER_PATCH = one_step.PRESSURE_COLUMNS_PER_PATCH
COLUMNS_PER_PATCH = one_step.COLUMNS_PER_PATCH
JOINT_COLUMNS = 2 * COLUMNS_PER_PATCH
MAX_STEPS = 5
STEP_TRUST_FRACTION = 1.0e-2
CUMULATIVE_TRUST_FRACTION = 5.0e-2
MAX_BACKTRACK = 12


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _save(report: dict[str, Any]) -> None:
    OUTPUT_PATH.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")


def _metric(values: np.ndarray, weights: np.ndarray, points: np.ndarray) -> dict[str, Any]:
    return one_step._metric(values, weights, points)


def _finite_coefficients(value: Any, count: int, name: str) -> np.ndarray:
    result = np.asarray(value, dtype=float)
    if result.shape != (count,) or not np.all(np.isfinite(result)):
        raise ValueError(f"{name} must contain {count} finite real values")
    return result.copy()


def _split_controls(values: Any, name: str) -> tuple[tuple[np.ndarray, np.ndarray], tuple[np.ndarray, np.ndarray]]:
    """Decode patch 1/2 blocks in [135 velocity, 45 pressure] order."""

    array = None if isinstance(values, dict) else np.asarray(values, dtype=float)
    if array is not None and array.shape == (2, COLUMNS_PER_PATCH):
        blocks = [array[0], array[1]]
    elif array is not None and array.shape == (2 * COLUMNS_PER_PATCH,):
        blocks = [array[:COLUMNS_PER_PATCH], array[COLUMNS_PER_PATCH:]]
    else:
        if not isinstance(values, dict):
            raise ValueError(f"{name} must contain a flat 360-vector or patch fields")
        blocks = []
        for patch in (1, 2):
            blocks.append(
                np.concatenate((
                    _finite_coefficients(
                        value[f"velocity_coefficients_patch{patch}"],
                        VELOCITY_COLUMNS_PER_PATCH,
                        f"{name} patch {patch} velocity",
                    ),
                    _finite_coefficients(
                        value[f"pressure_coefficients_patch{patch}"],
                        PRESSURE_COLUMNS_PER_PATCH,
                        f"{name} patch {patch} pressure",
                    ),
                ))
            )
    if len(blocks) != 2:
        raise ValueError(f"{name} must contain two patch blocks")
    decoded = []
    for patch, block in enumerate(blocks, 1):
        block = _finite_coefficients(block, COLUMNS_PER_PATCH, f"{name} patch {patch}")
        decoded.append((block[:VELOCITY_COLUMNS_PER_PATCH].copy(),
                        block[VELOCITY_COLUMNS_PER_PATCH:].copy()))
    return (decoded[0], decoded[1])


def _coefficients_from_report(report: dict[str, Any], name: str) -> tuple[tuple[np.ndarray, np.ndarray], tuple[np.ndarray, np.ndarray]]:
    section = report.get("selected")
    if not isinstance(section, dict):
        section = report.get("fit")
    if not isinstance(section, dict):
        raise ValueError(f"{name} has no selected or fit coefficient object")
    if "coefficients" in section:
        return _split_controls(section["coefficients"], name)
    return _split_controls(section, name)


def _load_seed() -> tuple[Path, dict[str, Any], tuple[tuple[np.ndarray, np.ndarray], tuple[np.ndarray, np.ndarray]]]:
    # The grouped trust worker may provide a better one-percent warm start.
    # The existing frozen step remains a deterministic fallback until that
    # report is available; no source is silently mixed with another schema.
    path = TRUST_PATH if TRUST_PATH.exists() else STEP_PATH
    report = json.loads(path.read_text(encoding="utf-8"))
    if report.get("status") != "completed":
        raise ValueError(f"Warm-start report is not completed: {path.name}")
    return path, report, _coefficients_from_report(report, path.name)


def _patch_geometry(data: dict[str, Any]) -> tuple[tuple[np.ndarray, np.ndarray], tuple[np.ndarray, np.ndarray]]:
    inputs = data["drift_report"]["inputs"]
    return (
        (np.asarray(inputs["first_patch_center"], dtype=float), np.asarray(inputs["first_patch_widths"], dtype=float)),
        (np.asarray(inputs["second_patch_center"], dtype=float), np.asarray(inputs["second_patch_widths"], dtype=float)),
    )


def _build_design(data: dict[str, Any], velocity: np.ndarray, jacobian: np.ndarray,
                  patch_geometry: tuple[tuple[np.ndarray, np.ndarray], tuple[np.ndarray, np.ndarray]],
                  pressure_design: np.ndarray) -> tuple[np.ndarray, list[tuple[Any, ...]], list[dict[str, Any]]]:
    points = data["points"]
    patch_data: list[tuple[Any, ...]] = []
    blocks: list[np.ndarray] = []
    layout: list[dict[str, Any]] = []
    for patch_index, (center, widths) in enumerate(patch_geometry, 1):
        velocity_response, values, jacobians, laplacians, mean_swirl, velocity_layout = one_step._mode_velocity_columns(
            points,
            center,
            widths,
            data["carrier"],
            data["nu"],
            data["h"],
            velocity,
            jacobian,
            data["laplacian"],
            data["tau0"],
        )
        p0 = (patch_index - 1) * PRESSURE_COLUMNS_PER_PATCH
        p1 = p0 + PRESSURE_COLUMNS_PER_PATCH
        blocks.extend((velocity_response, pressure_design[:, p0:p1]))
        patch_data.append((velocity_response, values, jacobians, laplacians, mean_swirl))
        for entry in velocity_layout:
            layout.append({"patch": patch_index, **entry})
        for index in range(PRESSURE_COLUMNS_PER_PATCH):
            layout.append({"patch": patch_index, "kind": "pressure_gradient", "index": index})
    design = np.column_stack(blocks)
    if design.shape != (3 * len(points), JOINT_COLUMNS):
        raise ValueError(f"Unexpected nonlinear response shape {design.shape}")
    return design, patch_data, layout


def _weighted_lstsq(design: np.ndarray, residual: np.ndarray,
                    weights: np.ndarray) -> tuple[np.ndarray, dict[str, Any]]:
    row_weight = np.repeat(np.sqrt(weights), 3)
    weighted_design = design * row_weight[:, None]
    weighted_rhs = -np.asarray(residual, dtype=float).reshape(-1) * row_weight
    scales = np.maximum(np.linalg.norm(weighted_design, axis=0), 1.0e-30)
    normalized = weighted_design / scales[None, :]
    singular = np.linalg.svd(normalized, compute_uv=False)
    cutoff = 1.0e-10 * max(float(singular[0]), 1.0e-300)
    rank = int(np.sum(singular > cutoff))
    solution, _, lstsq_rank, _ = np.linalg.lstsq(normalized, weighted_rhs, rcond=1.0e-10)
    return solution / scales, {
        "rank": rank,
        "lstsq_rank": int(lstsq_rank),
        "singular_values": singular.tolist(),
        "column_scale_min": float(np.min(scales)),
        "column_scale_max": float(np.max(scales)),
    }


def _velocity_and_jacobian(patch_data: list[tuple[Any, ...]], controls: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    velocity_blocks = []
    jacobian_blocks = []
    for patch, data in enumerate(patch_data):
        start = patch * COLUMNS_PER_PATCH
        coefficients = controls[start:start + VELOCITY_COLUMNS_PER_PATCH]
        velocity_blocks.append(np.einsum("niq,q->ni", data[1], coefficients))
        jacobian_blocks.append(np.einsum("nijq,q->nij", data[2], coefficients))
    return sum(velocity_blocks), sum(jacobian_blocks)


def _quadratic(patch_data: list[tuple[Any, ...]], controls: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    velocity, jacobian = _velocity_and_jacobian(patch_data, controls)
    return velocity, jacobian


def _cap_alpha(current_change: np.ndarray, proposed: np.ndarray, alpha: float,
               bound: float, weights: np.ndarray) -> float:
    def norm(value: np.ndarray) -> float:
        return float(np.sqrt(np.sum(weights * np.sum(value * value, axis=1))))
    if norm(current_change + alpha * proposed) <= bound * (1.0 + 1.0e-12):
        return alpha
    lo, hi = 0.0, alpha
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if norm(current_change + mid * proposed) <= bound:
            lo = mid
        else:
            hi = mid
    return lo * (1.0 - 1.0e-10)


def _selected_fields(coefficients: np.ndarray) -> dict[str, Any]:
    blocks = [coefficients[:COLUMNS_PER_PATCH], coefficients[COLUMNS_PER_PATCH:]]
    return {
        "coefficients": coefficients.tolist(),
        "velocity_coefficients_patch1": blocks[0][:VELOCITY_COLUMNS_PER_PATCH].tolist(),
        "pressure_coefficients_patch1": blocks[0][VELOCITY_COLUMNS_PER_PATCH:].tolist(),
        "velocity_coefficients_patch2": blocks[1][:VELOCITY_COLUMNS_PER_PATCH].tolist(),
        "pressure_coefficients_patch2": blocks[1][VELOCITY_COLUMNS_PER_PATCH:].tolist(),
    }


def run(output_path: Path = OUTPUT_PATH) -> dict[str, Any]:
    started = time.perf_counter()
    data = one_step._load_inputs()
    source_path, seed_report, seed_blocks = _load_seed()
    points = data["points"]
    weights = data["weights"]
    baseline, pressure_correction, pressure_design, pressure_meta = one_step._pressure_baseline(data)
    patch_geometry = _patch_geometry(data)
    design0, patch_data0, layout = _build_design(
        data, data["velocity"], data["jacobian"], patch_geometry, pressure_design
    )
    seed_controls = np.concatenate((
        np.concatenate(seed_blocks[0]), np.concatenate(seed_blocks[1])
    ))
    seed_velocity, seed_jacobian = _velocity_and_jacobian(patch_data0, seed_controls)
    seed_linear = design0 @ seed_controls
    seed_quadratic = np.einsum("nij,nj->ni", seed_jacobian, seed_velocity)
    current_velocity = data["velocity"] + seed_velocity
    current_jacobian = data["jacobian"] + seed_jacobian
    current_residual = baseline + seed_linear.reshape(-1, 3) + seed_quadratic
    cumulative_velocity = seed_velocity.copy()
    cumulative_controls = seed_controls.copy()
    velocity_norm = float(np.sqrt(np.sum(weights * np.sum(data["velocity"] ** 2, axis=1))))
    step_bound = STEP_TRUST_FRACTION * velocity_norm
    cumulative_bound = CUMULATIVE_TRUST_FRACTION * velocity_norm
    baseline_metric = _metric(baseline, weights, points)
    report: dict[str, Any] = {
        "status": "running",
        "accepted": False,
        "pde_validated": False,
        "scale_recursion_established": False,
        "constraints_maintained": False,
        "scope": "At most five analytic Gauss-Newton reference-velocity updates on the frozen scale-generator cache; no PDE or recursive-scale acceptance.",
        "sources": {
            "generator_report": {"path": GENERATOR_REPORT_PATH.name, "sha256": _sha256(GENERATOR_REPORT_PATH)},
            "generator_cache": {"path": GENERATOR_CACHE_PATH.name, "sha256": _sha256(GENERATOR_CACHE_PATH)},
            "pressure_report": {"path": PRESSURE_REPORT_PATH.name, "sha256": _sha256(PRESSURE_REPORT_PATH)},
            "drift_report": {"path": DRIFT_REPORT_PATH.name, "sha256": _sha256(DRIFT_REPORT_PATH)},
            "frozen_geometry": {"path": FROZEN_PATH.name, "sha256": _sha256(FROZEN_PATH)},
            "linearization": {"path": LINEARIZATION_PATH.name, "sha256": _sha256(LINEARIZATION_PATH)},
            "analytic_basis_jets": {"path": BASIS_PATH.name, "sha256": _sha256(BASIS_PATH)},
            "warm_start": {"path": source_path.name, "sha256": _sha256(source_path)},
        },
        "inputs": {
            "patch_count": 2,
            "correction_kind": "initial reference curl-potential velocity; time-tangent coefficients are not reused as velocity increments",
            "point_count": int(len(points)),
            "tau0": data["tau0"],
            "similarity_exponent_h": data["h"],
            "viscosity": data["nu"],
            "degree": DEGREE,
            "modes": list(MODES),
            "velocity_columns_per_patch": VELOCITY_COLUMNS_PER_PATCH,
            "pressure_columns_per_patch": PRESSURE_COLUMNS_PER_PATCH,
            "joint_column_count": JOINT_COLUMNS,
            "step_trust_fraction": STEP_TRUST_FRACTION,
            "cumulative_trust_fraction": CUMULATIVE_TRUST_FRACTION,
            "warm_start": source_path.name,
            "warm_start_accepted": bool(seed_report.get("selected", {}).get("feasible", False)),
        },
        "pressure_baseline": {
            "coefficients_patch1": pressure_meta["coefficients"][:PRESSURE_COLUMNS_PER_PATCH].tolist(),
            "coefficients_patch2": pressure_meta["coefficients"][PRESSURE_COLUMNS_PER_PATCH:].tolist(),
            "source": PRESSURE_REPORT_PATH.name,
        },
        "patches": {
            "first": {"center": patch_geometry[0][0].tolist(), "widths": patch_geometry[0][1].tolist()},
            "second": {"center": patch_geometry[1][0].tolist(), "widths": patch_geometry[1][1].tolist()},
            "joint_order": "patch 1 [135 velocity, 45 pressure], then patch 2 [135 velocity, 45 pressure]",
        },
        "baseline": baseline_metric,
        "warm_start": {"metrics": _metric(current_residual, weights, points), **_selected_fields(seed_controls)},
        "iterations": [],
    }
    _save(report)

    for iteration in range(1, MAX_STEPS + 1):
        design, patch_data, _ = _build_design(
            data, current_velocity, current_jacobian, patch_geometry, pressure_design
        )
        raw_controls, matrix = _weighted_lstsq(design, current_residual, weights)
        raw_velocity, raw_jacobian = _velocity_and_jacobian(patch_data, raw_controls)
        raw_quadratic = np.einsum("nij,nj->ni", raw_jacobian, raw_velocity)
        raw_norm = float(np.sqrt(np.sum(weights * np.sum(raw_velocity * raw_velocity, axis=1))))
        alpha_cap = min(1.0, step_bound / max(raw_norm, 1.0e-300))
        alpha_cap = _cap_alpha(cumulative_velocity, raw_velocity, alpha_cap, cumulative_bound, weights)
        current_metric = _metric(current_residual, weights, points)
        trials = []
        selected = None
        for backtrack in range(MAX_BACKTRACK):
            alpha = alpha_cap * (0.5 ** backtrack)
            trial_residual = current_residual + alpha * (design @ raw_controls).reshape(-1, 3) + alpha * alpha * raw_quadratic
            trial_velocity = alpha * raw_velocity
            trial_metric = _metric(trial_residual, weights, points)
            cumulative_trial = cumulative_velocity + trial_velocity
            step_metric = _metric(trial_velocity, weights, points)
            cumulative_metric = _metric(cumulative_trial, weights, points)
            row = {
                "backtrack": backtrack,
                "alpha": float(alpha),
                "residual": trial_metric,
                "step_velocity": step_metric,
                "cumulative_velocity": cumulative_metric,
                "quadratic_remainder": _metric(alpha * alpha * raw_quadratic, weights, points),
                "decreased": bool(trial_metric["volume_L2"] < current_metric["volume_L2"]),
                "step_trust_pass": bool(step_metric["volume_L2"] <= step_bound * (1.0 + 1.0e-10)),
                "cumulative_trust_pass": bool(cumulative_metric["volume_L2"] <= cumulative_bound * (1.0 + 1.0e-10)),
            }
            trials.append(row)
            if selected is None and row["decreased"] and row["step_trust_pass"] and row["cumulative_trust_pass"]:
                selected = (alpha, trial_residual, trial_velocity, alpha * raw_jacobian, trial_metric, step_metric, cumulative_metric, backtrack)
        if selected is None:
            report["iterations"].append({
                "iteration": iteration,
                "accepted": False,
                "current": current_metric,
                "matrix": matrix,
                "raw_velocity": _metric(raw_velocity, weights, points),
                "raw_alpha_cap": float(alpha_cap),
                "trials": trials,
            })
            _save(report)
            break
        alpha, trial_residual, trial_velocity, trial_jacobian, trial_metric, step_metric, cumulative_metric, backtrack = selected
        current_velocity += trial_velocity
        current_jacobian += trial_jacobian
        current_residual = trial_residual
        cumulative_velocity += trial_velocity
        cumulative_controls += alpha * raw_controls
        report["iterations"].append({
            "iteration": iteration,
            "accepted": True,
            "backtrack": int(backtrack),
            "alpha": float(alpha),
            "before": current_metric,
            "after": trial_metric,
            "step_velocity": step_metric,
            "cumulative_velocity": cumulative_metric,
            "raw_velocity": _metric(raw_velocity, weights, points),
            "raw_alpha_cap": float(alpha_cap),
            "quadratic_remainder": _metric(alpha * alpha * raw_quadratic, weights, points),
            "matrix": matrix,
            "trials": trials,
        })
        report["warm_start"]["metrics"] = _metric(current_residual, weights, points)
        _save(report)

    report["selected"] = {
        "feasible": bool(any(row.get("accepted", False) for row in report["iterations"])),
        "iterations_accepted": int(sum(bool(row.get("accepted", False)) for row in report["iterations"])),
        "metrics": _metric(current_residual, weights, points),
        "velocity_correction": _metric(cumulative_velocity, weights, points),
        **_selected_fields(cumulative_controls),
    }
    report["status"] = "completed"
    report["elapsed_seconds"] = float(time.perf_counter() - started)
    _save(report)
    print(json.dumps({
        "status": report["status"],
        "accepted_steps": report["selected"]["iterations_accepted"],
        "initial_l2": report["baseline"]["volume_L2"],
        "final_l2": report["selected"]["metrics"]["volume_L2"],
        "final_max": report["selected"]["metrics"]["max_norm"],
        "elapsed_seconds": report["elapsed_seconds"],
    }), flush=True)
    return report


if __name__ == "__main__":
    run()
