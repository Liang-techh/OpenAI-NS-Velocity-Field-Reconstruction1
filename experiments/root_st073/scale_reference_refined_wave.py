"""Bounded original mode-1 wave-amplitude fit on the refined replay grid.

The frozen corrected reference and its Cartesian jets are loaded from the
completed 44,400-point replay cache.  Only the amplitude of the original
degree-2, mode-1 wave is varied; the mean and two compact local velocity
corrections remain fixed.  A 90-column compact pressure-gradient increment is
eliminated by an orthogonal SVD projection.  The projected scalar objective is
quartic in the amplitude increment because the momentum response is

    R(a) = R_0 + a L + a^2 Q.

All stationary roots and trust-region endpoints are replayed with the exact
quadratic residual.  A candidate is accepted only when both refined-grid
volume L2 and maximum decrease strictly against the current corrected field.
This is a spatial diagnostic and does not establish PDE validity or scale
recursion.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import time
from pathlib import Path

for _key in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS"):
    os.environ[_key] = "1"

import numpy as np

ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scale_generator_linearization import velocity_correction_response  # noqa: E402
from scale_reference_velocity_step import _metric, _pressure_patch_design  # noqa: E402
from supported_fourier_analytic_jets import basis_jets  # noqa: E402


ACTUAL_REPLAY_PATH = ROOT / "scale_reference_trust_nonlinear_actual.json"
GEOMETRY_PATH = ROOT / "full_wave_frozen_cache.json"
DRIFT_PATH = ROOT / "localized_drift_1800_fit.json"
COARSE_PROBE_PATH = ROOT / "scale_wave_amplitude_probe.json"
OUTPUT_PATH = ROOT / "scale_reference_refined_wave.json"
RIDGE_RELATIVE = 1.0e-12


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _array_sha256(*arrays) -> str:
    digest = hashlib.sha256()
    for array in arrays:
        value = np.ascontiguousarray(np.asarray(array))
        digest.update(str(value.dtype).encode("ascii"))
        digest.update(np.asarray(value.shape, dtype=np.int64).tobytes())
        digest.update(value.tobytes())
    return digest.hexdigest()


def _save(path: Path, report: dict) -> None:
    path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")


def _load_cache(report: dict) -> tuple[np.ndarray, ...]:
    cache = report["chunk_caches"]["corrected_reference"]
    cache_dir = Path(cache["directory"])
    paths = sorted(cache_dir.glob("*.npz"), key=lambda value: int(value.stem))
    if len(paths) != int(report["cases"]["corrected_reference"]["chunks_completed"]):
        raise ValueError("Refined cache chunk count does not match completed replay")
    arrays = {key: [] for key in ("points", "weights", "velocity", "jacobian", "laplacian", "pressure_gradient", "residual")}
    expected = 0
    for path in paths:
        begin = int(path.stem)
        if begin != expected:
            raise ValueError(f"Unexpected chunk offset {begin}; expected {expected}")
        with np.load(path, allow_pickle=False) as loaded:
            missing = sorted(set(arrays).difference(loaded.files))
            if missing:
                raise ValueError(f"Chunk {path.name} missing {missing}")
            for key in arrays:
                arrays[key].append(np.asarray(loaded[key], dtype=float))
            expected += len(loaded["points"])
    result = tuple(np.concatenate(arrays[key], axis=0) for key in arrays)
    if result[0].shape != (int(report["grid"]["point_count"]), 3):
        raise ValueError(f"Unexpected refined grid shape {result[0].shape}")
    return result


def _metric_from_flat(values: np.ndarray, weights: np.ndarray, points: np.ndarray) -> dict:
    return _metric(np.asarray(values, dtype=float).reshape(-1, 3), weights, points)


def _pressure_columns(points: np.ndarray, carrier: np.ndarray) -> tuple[np.ndarray, dict]:
    drift = json.loads(DRIFT_PATH.read_text(encoding="utf-8"))
    inputs = drift["inputs"]
    first_center = np.asarray(inputs["first_patch_center"], dtype=float)
    first_widths = np.asarray(inputs["first_patch_widths"], dtype=float)
    second_center = np.asarray(inputs["second_patch_center"], dtype=float)
    second_widths = np.asarray(inputs["second_patch_widths"], dtype=float)
    first = _pressure_patch_design(points, first_center, first_widths, carrier)
    second = _pressure_patch_design(points, second_center, second_widths, carrier)
    design = np.column_stack((first, second))
    if design.shape != (3 * len(points), 90):
        raise ValueError(f"Unexpected refined pressure design shape {design.shape}")
    return design, {
        "first_center": first_center.tolist(),
        "first_widths": first_widths.tolist(),
        "second_center": second_center.tolist(),
        "second_widths": second_widths.tolist(),
    }


def _fit_pressure_projection(raw: np.ndarray, pressure_design: np.ndarray,
                             row_weight: np.ndarray, scales: np.ndarray,
                             left: np.ndarray, singular: np.ndarray,
                             right: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    weighted_rhs = -np.asarray(raw, dtype=float).reshape(-1) * row_weight
    projected = left.T @ weighted_rhs
    normalized_coefficients = right.T @ (projected / singular)
    pressure_coefficients = normalized_coefficients / scales
    corrected = np.asarray(raw, dtype=float).reshape(-1) + pressure_design @ pressure_coefficients
    return corrected.reshape(-1, 3), pressure_coefficients


def run(output_path: Path | str = OUTPUT_PATH,
        replay_path: Path | str = ACTUAL_REPLAY_PATH) -> dict:
    started = time.perf_counter()
    output_path = Path(output_path)
    replay_path = Path(replay_path)
    replay = json.loads(replay_path.read_text(encoding="utf-8"))
    if replay.get("status") != "completed":
        raise ValueError("Actual refined replay is not completed")
    if replay["cases"]["corrected_reference"].get("status") != "completed":
        raise ValueError("Corrected refined replay case is not completed")
    points, weights, base_velocity, base_jacobian, base_laplacian, base_pressure_gradient, baseline = _load_cache(replay)
    geometry = json.loads(GEOMETRY_PATH.read_text(encoding="utf-8"))["inputs"]["wave"]
    degree = int(geometry["degree"])
    carrier = np.asarray(geometry["carrier"], dtype=float)
    wave_packed = np.asarray(geometry["coefficients_original"], dtype=float)
    if wave_packed.shape != (27, 2):
        raise ValueError(f"Unexpected original wave coefficient shape {wave_packed.shape}")
    wave_coefficients = wave_packed[:, 0] + 1j * wave_packed[:, 1]
    tau0 = float(replay["inputs"]["tau0"])
    h = float(replay["inputs"]["similarity_exponent_h"])
    nu = float(replay["inputs"]["viscosity"])
    report = {
        "status": "initializing",
        "accepted": False,
        "pde_validated": False,
        "scale_recursion_established": False,
        "scope": (
            "Refined-grid scalar amplitude diagnostic. The original degree-2 mode-1 wave "
            "in full_wave_frozen_cache is varied around the frozen corrected reference; "
            "the mean and two compact local velocity corrections are held fixed. A compact "
            "90-column pressure-gradient increment is fit relative to the current corrected "
            "pressure. The exact quadratic convection response is replayed on all 44,400 "
            "points. No moments, cones, shape, trajectory, PDE or scale-recursion gate is "
            "claimed."
        ),
        "sources": {
            "actual_replay": {"path": replay_path.name, "sha256": _sha256(replay_path)},
            "geometry": {"path": GEOMETRY_PATH.name, "sha256": _sha256(GEOMETRY_PATH)},
            "drift_geometry": {"path": DRIFT_PATH.name, "sha256": _sha256(DRIFT_PATH)},
            "coarse_normalization_probe": {"path": COARSE_PROBE_PATH.name, "sha256": _sha256(COARSE_PROBE_PATH)},
            "linearization": {"path": (ROOT / "scale_generator_linearization.py").name, "sha256": _sha256(ROOT / "scale_generator_linearization.py")},
            "analytic_basis_jets": {"path": (ROOT / "supported_fourier_analytic_jets.py").name, "sha256": _sha256(ROOT / "supported_fourier_analytic_jets.py")},
            "pressure_design": {"path": (ROOT / "scale_reference_velocity_step.py").name, "sha256": _sha256(ROOT / "scale_reference_velocity_step.py")},
        },
        "inputs": {
            "point_count": int(len(points)),
            "physical_volume": float(np.sum(weights)),
            "grid_array_sha256": _array_sha256(points, weights),
            "tau0": tau0,
            "similarity_exponent_h": h,
            "viscosity": nu,
            "degree": degree,
            "mode": 1,
            "carrier": carrier.tolist(),
            "wave_center": np.asarray(geometry["center"], dtype=float).tolist(),
            "wave_widths": np.asarray(geometry["widths"], dtype=float).tolist(),
            "wave_coefficient_count": int(len(wave_coefficients)),
            "extra_velocity_norm_reference": "scale_wave_amplitude_probe.json reference_velocity_L2 on the coarse 18,720-point normalization",
            "trust_fraction": 1.0e-2,
            "pressure_columns": 90,
            "pressure_semantics": "pressure coefficients are increments added to the current corrected reference pressure; the existing corrected pressure remains the baseline",
        },
        "wave": {
            "coefficients_original": wave_packed.tolist(),
            "coefficient_source": "full_wave_frozen_cache.inputs.wave.coefficients_original",
            "velocity_column_semantics": "one real scalar amplitude increment times the full original mode-1 wave",
        },
    }
    _save(output_path, report)
    print(json.dumps({"stage": "sources_loaded", "point_count": len(points)}), flush=True)

    values, gradients, viscous, _, _ = basis_jets(
        points, geometry["center"], geometry["widths"], 1, degree, carrier, nu
    )
    wave_velocity = np.einsum("niq,q->ni", values, wave_coefficients).real
    wave_jacobian = np.einsum("nijq,q->nij", gradients, wave_coefficients).real
    wave_laplacian = np.einsum("niq,q->ni", viscous, wave_coefficients).real
    wave_laplacian = -wave_laplacian / nu
    wave_linear, wave_quadratic = velocity_correction_response(
        points,
        base_velocity,
        base_jacobian,
        wave_velocity,
        wave_jacobian,
        wave_laplacian,
        tau0,
        h,
        np.zeros(len(points), dtype=float),
        nu,
    )
    row_weight = np.repeat(np.sqrt(weights), 3)
    pressure_design, pressure_geometry = _pressure_columns(points, carrier)
    weighted_pressure = pressure_design * row_weight[:, None]
    pressure_scales = np.maximum(np.linalg.norm(weighted_pressure, axis=0), 1.0e-30)
    pressure_normalized = weighted_pressure / pressure_scales[None, :]
    pressure_left, pressure_singular, pressure_right = np.linalg.svd(
        pressure_normalized, full_matrices=False
    )
    pressure_cutoff = RIDGE_RELATIVE * max(float(pressure_singular[0]), 1.0e-300)
    keep = pressure_singular > pressure_cutoff
    pressure_rank = int(np.sum(keep))
    pressure_left = pressure_left[:, keep]
    pressure_singular = pressure_singular[keep]
    pressure_right = pressure_right[keep, :]
    baseline_weighted = baseline.reshape(-1) * row_weight
    linear_weighted = wave_linear.reshape(-1) * row_weight
    quadratic_weighted = wave_quadratic.reshape(-1) * row_weight
    projected_vectors = []
    for vector in (baseline_weighted, linear_weighted, quadratic_weighted):
        projected_vectors.append(vector - pressure_left @ (pressure_left.T @ vector))
    projected_baseline, projected_linear, projected_quadratic = projected_vectors
    polynomial = np.array([
        projected_baseline @ projected_baseline,
        2.0 * (projected_baseline @ projected_linear),
        projected_linear @ projected_linear + 2.0 * (projected_baseline @ projected_quadratic),
        2.0 * (projected_linear @ projected_quadratic),
        projected_quadratic @ projected_quadratic,
    ], dtype=float)
    coarse_probe = json.loads(COARSE_PROBE_PATH.read_text(encoding="utf-8"))
    coarse_reference_norm = float(coarse_probe["reference_velocity_L2"])
    wave_velocity_norm = float(np.sqrt(np.sum(weights * np.sum(wave_velocity * wave_velocity, axis=1))))
    trust_bound = 1.0e-2 * coarse_reference_norm
    amplitude_bound = trust_bound / max(wave_velocity_norm, 1.0e-300)
    derivative_polynomial = polynomial[1:] * np.arange(1, 5, dtype=float)
    roots = np.polynomial.polynomial.polyroots(derivative_polynomial)
    candidate_shifts = [-amplitude_bound, 0.0, amplitude_bound]
    for root in roots:
        if abs(float(root.imag)) <= 1.0e-9 * max(1.0, abs(float(root.real))):
            value = float(root.real)
            if -amplitude_bound <= value <= amplitude_bound:
                candidate_shifts.append(value)
    candidate_shifts = sorted(set(float(value) for value in candidate_shifts))
    baseline_metric = _metric(baseline, weights, points)
    evaluations = []
    for shift in candidate_shifts:
        raw = baseline + shift * wave_linear + shift * shift * wave_quadratic
        corrected, pressure_coefficients = _fit_pressure_projection(
            raw, pressure_design, row_weight, pressure_scales,
            pressure_left, pressure_singular, pressure_right,
        )
        metric = _metric(corrected, weights, points)
        evaluations.append({
            "amplitude_shift": float(shift),
            "wave_amplitude": float(1.0 + shift),
            "pressure_coefficients": pressure_coefficients.tolist(),
            "velocity_change_volume_L2": float(abs(shift) * wave_velocity_norm),
            "trust_region_pass": bool(abs(shift) * wave_velocity_norm <= trust_bound * (1.0 + 1.0e-12)),
            "residual": metric,
            "l2_change": float(metric["volume_L2"] - baseline_metric["volume_L2"]),
            "max_change": float(metric["max_norm"] - baseline_metric["max_norm"]),
            "l2_relative_change": float(metric["volume_L2"] / max(baseline_metric["volume_L2"], 1.0e-300) - 1.0),
            "max_relative_change": float(metric["max_norm"] / max(baseline_metric["max_norm"], 1.0e-300) - 1.0),
            "both_metrics_improve": bool(metric["volume_L2"] < baseline_metric["volume_L2"] and metric["max_norm"] < baseline_metric["max_norm"]),
        })
    accepted = [row for row in evaluations if row["trust_region_pass"] and row["both_metrics_improve"]]
    if accepted:
        selected_row = min(accepted, key=lambda row: (row["residual"]["volume_L2"], row["residual"]["max_norm"]))
        accepted_flag = True
        reason = "selected minimum refined-grid volume L2 among candidates improving both L2 and maximum"
    else:
        selected_row = {
            "amplitude_shift": 0.0,
            "wave_amplitude": 1.0,
            "pressure_coefficients": np.zeros(90, dtype=float).tolist(),
            "velocity_change_volume_L2": 0.0,
            "trust_region_pass": True,
            "residual": baseline_metric,
            "l2_change": 0.0,
            "max_change": 0.0,
            "l2_relative_change": 0.0,
            "max_relative_change": 0.0,
            "both_metrics_improve": False,
        }
        accepted_flag = False
        reason = "rejected: no stationary root or endpoint improved both refined-grid L2 and maximum"
    report["pressure_projection"] = {
        "rank": pressure_rank,
        "column_count": 90,
        "cutoff": float(pressure_cutoff),
        "singular_values": pressure_singular.tolist(),
        "column_scale_min": float(np.min(pressure_scales)),
        "column_scale_max": float(np.max(pressure_scales)),
        "geometry": pressure_geometry,
    }
    report["trust_region"] = {
        "coarse_reference_velocity_norm": coarse_reference_norm,
        "refined_wave_velocity_norm": wave_velocity_norm,
        "velocity_trust_bound": trust_bound,
        "amplitude_shift_bound": amplitude_bound,
        "wave_velocity_metric": _metric(wave_velocity, weights, points),
    }
    report["quartic"] = {
        "objective_definition": "||Proj_pressure(R0 + a L + a^2 Q)||_weighted^2",
        "coefficients_ascending": polynomial.tolist(),
        "derivative_coefficients_ascending": derivative_polynomial.tolist(),
        "stationary_roots": [complex(root).real if abs(complex(root).imag) < 1.0e-9 else [complex(root).real, complex(root).imag] for root in roots],
        "candidate_shifts": candidate_shifts,
        "pressure_projection_is_orthogonal": True,
    }
    report["baseline"] = baseline_metric
    report["candidate_evaluations"] = evaluations
    report["selected"] = {
        **selected_row,
        "accepted": accepted_flag,
        "reason": reason,
        "baseline_pressure_semantics": "current corrected pressure from frozen replay; selected pressure coefficients are increments",
    }
    report["status"] = "completed"
    report["elapsed_seconds"] = float(time.perf_counter() - started)
    _save(output_path, report)
    print(json.dumps({
        "stage": "completed",
        "accepted": accepted_flag,
        "amplitude_shift": selected_row["amplitude_shift"],
        "baseline_l2": baseline_metric["volume_L2"],
        "selected_l2": selected_row["residual"]["volume_L2"],
        "baseline_max": baseline_metric["max_norm"],
        "selected_max": selected_row["residual"]["max_norm"],
        "candidate_count": len(evaluations),
        "elapsed_seconds": report["elapsed_seconds"],
    }), flush=True)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT_PATH)
    parser.add_argument("--replay", type=Path, default=ACTUAL_REPLAY_PATH)
    args = parser.parse_args()
    run(args.output, args.replay)
