"""Analytic transfer of the refined wave candidate to the original coarse grid.

The 18,720-point generator cache already stores the base velocity jets and
residual.  This module rebuilds the current trust-nonlinear parent residual by
adding its saved compact local velocity and pressure controls analytically,
then adds the selected refined wave amplitude response and pressure increment.
No finite differences or new field replay are used.  The result is a
cross-grid comparison, not a fresh independent holdout.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scale_generator_linearization import velocity_correction_response  # noqa: E402
from scale_reference_candidate import _unpack_velocity_only  # noqa: E402
from scale_reference_refined_wave_adapter import load_reference  # noqa: E402
from scale_reference_velocity_step import _metric, _pressure_patch_design  # noqa: E402
from supported_fourier_analytic_jets import basis_jets  # noqa: E402


ROOT_CACHE_PATH = ROOT / "scale_generator_momentum_defect.npz"
ROOT_CACHE_REPORT_PATH = ROOT / "scale_generator_momentum_defect.json"
TRUST_STEP_PATH = ROOT / "scale_reference_trust_nonlinear_fit.json"
PRESSURE_PATH = ROOT / "scale_generator_pressure_projection.json"
DRIFT_PATH = ROOT / "localized_drift_1800_fit.json"
GEOMETRY_PATH = ROOT / "full_wave_frozen_cache.json"
REFINED_REPORT_PATH = ROOT / "scale_reference_refined_wave.json"
ADAPTER_PATH = ROOT / "scale_reference_refined_wave_adapter.py"
OUTPUT_PATH = ROOT / "scale_reference_refined_wave_transfer.json"
MEAN_ANGLES = 12


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _save(path: Path, report: dict) -> None:
    path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")


def _mean_swirl(points: np.ndarray, velocity: np.ndarray) -> np.ndarray:
    count = len(points)
    if count % MEAN_ANGLES:
        raise ValueError("Coarse grid is not ordered in 12-angle rings")
    radius = np.linalg.norm(points[:, :2], axis=1)
    if np.any(radius <= 0.0):
        raise ValueError("Coarse grid contains an axis point")
    e_theta = np.column_stack((-points[:, 1] / radius, points[:, 0] / radius, np.zeros(count)))
    swirl = np.einsum("ni,ni->n", e_theta, velocity)
    return np.repeat(np.mean(swirl.reshape(-1, MEAN_ANGLES), axis=1), MEAN_ANGLES)


def _local_jets(points: np.ndarray, controls: np.ndarray,
                center: np.ndarray, widths: np.ndarray,
                carrier: np.ndarray, degree: int, nu: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    modes = _unpack_velocity_only(controls, "local velocity controls")
    value_sum = np.zeros((len(points), 3), dtype=float)
    gradient_sum = np.zeros((len(points), 3, 3), dtype=float)
    laplacian_sum = np.zeros((len(points), 3), dtype=float)
    for mode in (0, 1, 2):
        mode_carrier = np.zeros(2, dtype=float) if mode == 0 else mode * carrier
        values, gradients, viscous, _, _ = basis_jets(
            points, center, widths, mode, degree, mode_carrier, nu
        )
        coefficient = modes[mode]
        value_sum += np.einsum("niq,q->ni", values, coefficient).real
        gradient_sum += np.einsum("nijq,q->nij", gradients, coefficient).real
        laplacian_sum += -np.einsum("niq,q->ni", viscous, coefficient).real / nu
    return value_sum, gradient_sum, laplacian_sum


def _load_cache() -> tuple[np.ndarray, ...]:
    with np.load(ROOT_CACHE_PATH, allow_pickle=False) as loaded:
        required = {"points", "weights", "velocity", "gradient", "pressure_gradient", "laplacian", "generator_residual"}
        missing = sorted(required.difference(loaded.files))
        if missing:
            raise ValueError(f"Coarse generator cache missing {missing}")
        return tuple(np.asarray(loaded[key], dtype=float) for key in (
            "points", "weights", "velocity", "gradient", "pressure_gradient", "laplacian", "generator_residual"
        ))


def run(output_path: Path | str = OUTPUT_PATH) -> dict:
    started = time.perf_counter()
    output_path = Path(output_path)
    refined = json.loads(REFINED_REPORT_PATH.read_text(encoding="utf-8"))
    if refined.get("status") != "completed" or not refined["selected"].get("accepted", False):
        raise ValueError("Refined wave report is not a completed accepted candidate")
    points, weights, base_velocity, base_jacobian, _, base_laplacian, generator_residual = _load_cache()
    root_report = json.loads(ROOT_CACHE_REPORT_PATH.read_text(encoding="utf-8"))
    if root_report.get("status") not in ("completed", "fields_ready"):
        raise ValueError("Coarse generator cache report is not completed")
    step = json.loads(TRUST_STEP_PATH.read_text(encoding="utf-8"))
    geometry = json.loads(GEOMETRY_PATH.read_text(encoding="utf-8"))["inputs"]["wave"]
    drift = json.loads(DRIFT_PATH.read_text(encoding="utf-8"))
    pressure_report = json.loads(PRESSURE_PATH.read_text(encoding="utf-8"))
    tau0 = float(root_report["inputs"]["tau0"])
    h = float(root_report["inputs"]["similarity_exponent_h"])
    nu = float(root_report["inputs"]["viscosity"])
    carrier = np.asarray(geometry["carrier"], dtype=float)
    degree = int(geometry["degree"])
    selected_step = step["selected"]
    patch_centers = (
        np.asarray(drift["inputs"]["first_patch_center"], dtype=float),
        np.asarray(drift["inputs"]["second_patch_center"], dtype=float),
    )
    patch_widths = (
        np.asarray(drift["inputs"]["first_patch_widths"], dtype=float),
        np.asarray(drift["inputs"]["second_patch_widths"], dtype=float),
    )
    pressure_design = np.column_stack((
        _pressure_patch_design(points, patch_centers[0], patch_widths[0], carrier),
        _pressure_patch_design(points, patch_centers[1], patch_widths[1], carrier),
    ))
    if pressure_design.shape != (3 * len(points), 90):
        raise ValueError(f"Unexpected coarse pressure design shape {pressure_design.shape}")
    pressure_base = np.asarray(pressure_report["fit"]["coefficients"], dtype=float)
    pressure_step = np.concatenate((
        np.asarray(selected_step["pressure_coefficients_patch1"], dtype=float),
        np.asarray(selected_step["pressure_coefficients_patch2"], dtype=float),
    ))
    if pressure_base.shape != (90,) or pressure_step.shape != (90,):
        raise ValueError("Expected 90 baseline and 90 selected trust pressure coefficients")
    local_velocity = np.zeros_like(base_velocity)
    local_jacobian = np.zeros_like(base_jacobian)
    local_laplacian = np.zeros_like(base_laplacian)
    local_blocks = []
    for number, (center, widths) in enumerate(zip(patch_centers, patch_widths), 1):
        controls = np.asarray(selected_step[f"velocity_coefficients_patch{number}"], dtype=float)
        value, gradient, laplacian = _local_jets(points, controls, center, widths, carrier, degree, nu)
        local_velocity += value
        local_jacobian += gradient
        local_laplacian += laplacian
        local_blocks.append({"patch": number, "control_count": int(len(controls)), "coefficient_norm": float(np.linalg.norm(controls))})
    local_mean_swirl = _mean_swirl(points, local_velocity)
    local_linear, local_quadratic = velocity_correction_response(
        points,
        base_velocity,
        base_jacobian,
        local_velocity,
        local_jacobian,
        local_laplacian,
        tau0,
        h,
        local_mean_swirl,
        nu,
    )
    parent_residual = (
        generator_residual
        + local_linear
        + local_quadratic
        + (pressure_design @ (pressure_base + pressure_step)).reshape(-1, 3)
    )
    parent_metric = _metric(parent_residual, weights, points)
    expected_parent_metric = step["selected"]["metrics"]
    adapter = load_reference()
    wave_velocity, wave_jacobian, wave_laplacian = adapter.wave_jets(points)
    wave_linear, wave_quadratic = velocity_correction_response(
        points,
        base_velocity + local_velocity,
        base_jacobian + local_jacobian,
        wave_velocity,
        wave_jacobian,
        wave_laplacian,
        tau0,
        h,
        np.zeros(len(points), dtype=float),
        nu,
    )
    amplitude_shift = float(refined["selected"]["amplitude_shift"])
    refined_pressure = np.asarray(refined["selected"]["pressure_coefficients"], dtype=float)
    selected_residual = (
        parent_residual
        + amplitude_shift * wave_linear
        + amplitude_shift * amplitude_shift * wave_quadratic
        + (pressure_design @ refined_pressure).reshape(-1, 3)
    )
    selected_metric = _metric(selected_residual, weights, points)
    report = {
        "status": "completed",
        "accepted": False,
        "pde_validated": False,
        "scale_recursion_established": False,
        "scope": (
            "Analytic transfer of the accepted refined original-wave amplitude candidate "
            "to the original 18,720-point generator cache. The current trust-nonlinear "
            "parent is reconstructed from cached base jets plus its saved local and pressure "
            "controls; the selected wave and pressure increments are then added analytically. "
            "This is a cross-grid comparison, not a fresh holdout or PDE acceptance test."
        ),
        "sources": {
            "coarse_cache": {"path": ROOT_CACHE_PATH.name, "sha256": _sha256(ROOT_CACHE_PATH)},
            "coarse_report": {"path": ROOT_CACHE_REPORT_PATH.name, "sha256": _sha256(ROOT_CACHE_REPORT_PATH)},
            "trust_step": {"path": TRUST_STEP_PATH.name, "sha256": _sha256(TRUST_STEP_PATH)},
            "pressure_projection": {"path": PRESSURE_PATH.name, "sha256": _sha256(PRESSURE_PATH)},
            "drift_report": {"path": DRIFT_PATH.name, "sha256": _sha256(DRIFT_PATH)},
            "geometry": {"path": GEOMETRY_PATH.name, "sha256": _sha256(GEOMETRY_PATH)},
            "refined_wave_report": {"path": REFINED_REPORT_PATH.name, "sha256": _sha256(REFINED_REPORT_PATH)},
            "adapter": {"path": ADAPTER_PATH.name, "sha256": _sha256(ADAPTER_PATH)},
        },
        "inputs": {
            "point_count": int(len(points)),
            "physical_volume": float(np.sum(weights)),
            "tau0": tau0,
            "similarity_exponent_h": h,
            "viscosity": nu,
            "grid_role": "original 18,720-point generator grid; cross-grid transfer, not fresh holdout",
            "pressure_semantics": "coarse parent includes pressure projection plus trust-step pressure; refined pressure is an increment to that parent",
            "amplitude_shift": amplitude_shift,
        },
        "parent_reconstruction": {
            "metric": parent_metric,
            "expected_trust_step_metric": expected_parent_metric,
            "volume_L2_abs_difference": float(parent_metric["volume_L2"] - expected_parent_metric["volume_L2"]),
            "max_abs_difference": float(parent_metric["max_norm"] - expected_parent_metric["max_norm"]),
            "local_blocks": local_blocks,
            "pressure_base_coefficient_count": 90,
            "pressure_step_coefficient_count": 90,
        },
        "selected_wave_transfer": {
            "metric": selected_metric,
            "parent_metric": parent_metric,
            "volume_L2_change": float(selected_metric["volume_L2"] - parent_metric["volume_L2"]),
            "max_change": float(selected_metric["max_norm"] - parent_metric["max_norm"]),
            "volume_L2_relative_change": float(selected_metric["volume_L2"] / max(parent_metric["volume_L2"], 1.0e-300) - 1.0),
            "max_relative_change": float(selected_metric["max_norm"] / max(parent_metric["max_norm"], 1.0e-300) - 1.0),
            "both_metrics_improve_on_coarse_transfer": bool(selected_metric["volume_L2"] < parent_metric["volume_L2"] and selected_metric["max_norm"] < parent_metric["max_norm"]),
            "wave_linear_metric": _metric(wave_linear, weights, points),
            "wave_quadratic_metric": _metric(wave_quadratic, weights, points),
            "pressure_increment_coefficient_count": int(len(refined_pressure)),
        },
        "limitations": [
            "The comparison reuses the original generator grid and is not a new spatial holdout.",
            "No finite-difference replay, moments, cones, shape, trajectory, PDE, or scale-recursion gate is evaluated here.",
            "The parent reconstruction relies on saved analytic compact-control coefficients and cached base jets; its agreement with the frozen trust-step metrics is reported explicitly.",
        ],
        "elapsed_seconds": float(time.perf_counter() - started),
    }
    _save(output_path, report)
    print(json.dumps({
        "stage": "completed",
        "parent_l2": parent_metric["volume_L2"],
        "selected_l2": selected_metric["volume_L2"],
        "parent_max": parent_metric["max_norm"],
        "selected_max": selected_metric["max_norm"],
        "both_improve": report["selected_wave_transfer"]["both_metrics_improve_on_coarse_transfer"],
        "parent_rebuild_l2_abs": report["parent_reconstruction"]["volume_L2_abs_difference"],
        "elapsed_seconds": report["elapsed_seconds"],
    }), flush=True)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT_PATH)
    args = parser.parse_args()
    run(args.output)
