"""Pressure-only projection of the frozen scale-generator momentum defect.

The grouped generator diagnostic supplies an independent grid and its frozen
generator residual.  This module builds only the pressure-gradient columns of
the two compact degree-2 local patches (45 columns per patch, 90 total) and
solves a weighted least-squares projection.  Velocity and generator values are
held fixed.  The result is a pressure-capacity diagnostic, with no cone,
full-domain, PDE, forcing, or recursive-scale acceptance claim.
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

from localized_residual_enrichment import _mode_block_columns  # noqa: E402


DRIFT_PATH = ROOT / "localized_drift_1800_fit.json"
DRIFT_SOURCE_PATH = ROOT / "localized_drift_1800_fit.py"
LOCAL_BASIS_SOURCE_PATH = ROOT / "localized_residual_enrichment.py"
GENERATOR_SCRIPT_PATH = ROOT / "scale_generator_momentum_defect.py"
OUTPUT_PATH = ROOT / "scale_generator_pressure_projection.json"
GENERATOR_REPORT_CANDIDATES = (
    ROOT / "scale_generator_momentum_defect.json",
    ROOT / "scale_generator_momentum_defect_report.json",
)
GENERATOR_CACHE_CANDIDATES = (
    ROOT / "scale_generator_momentum_defect.npz",
    ROOT / "scale_generator_momentum_defect_cache.npz",
    ROOT / "scale_generator_momentum_defect_grid.npz",
)
DEGREE = 2
MODES = (0, 1, 2)
PRESSURE_COLUMN_COUNT_PER_PATCH = 45
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


def _metric(values, weights, points):
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


def _pressure_layout():
    layout = []
    q = (DEGREE + 1) ** 2
    for mode, components in ((0, ("real",)), (1, ("real", "imag")), (2, ("real", "imag"))):
        for index in range(q):
            for component in components:
                layout.append({"patch_local_mode": mode, "kind": "pressure_gradient", "component": component, "index": index})
    if len(layout) != PRESSURE_COLUMN_COUNT_PER_PATCH:
        raise ValueError(f"Unexpected pressure layout length {len(layout)}")
    return layout


def _pressure_patch_design(points, center, widths, carrier):
    """Return only pressure-gradient columns in the known 180-column layout."""

    mode0_full = _mode_block_columns(points, center, widths, 0, DEGREE, np.zeros(2, dtype=float))
    mode1_full = _mode_block_columns(points, center, widths, 1, DEGREE, carrier)
    mode2_full = _mode_block_columns(points, center, widths, 2, DEGREE, 2.0 * carrier)
    # _mode_block_columns orders 27 velocity columns then 9 pressure-gradient
    # columns, with real/-imag interleaving for nonzero modes.  Mode 0's local
    # layout keeps only the real member of each pair.
    mode0_pressure = mode0_full[:, ::2][:, 27:36]
    mode1_pressure = mode1_full[:, 54:72]
    mode2_pressure = mode2_full[:, 54:72]
    design = np.column_stack((mode0_pressure, mode1_pressure, mode2_pressure))
    if design.shape[1] != PRESSURE_COLUMN_COUNT_PER_PATCH:
        raise ValueError(f"Unexpected pressure design shape {design.shape}")
    return design


def _find_report_and_cache():
    report_path = next((path for path in GENERATOR_REPORT_CANDIDATES if path.exists()), None)
    cache_path = next((path for path in GENERATOR_CACHE_CANDIDATES if path.exists()), None)
    return report_path, cache_path


def run(output_path=OUTPUT_PATH):
    started = time.perf_counter()
    output_path = Path(output_path)
    report_path, cache_path = _find_report_and_cache()
    if report_path is None or cache_path is None:
        raise FileNotFoundError(
            "Waiting for grouped generator report/NPZ with points, weights, and generator_residual"
        )
    generator_report = json.loads(report_path.read_text(encoding="utf-8"))
    if generator_report.get("status") not in ("completed", "fields_ready"):
        raise ValueError(f"Generator report is not frozen: {generator_report.get('status')}")
    with np.load(cache_path, allow_pickle=False) as loaded:
        required = {"points", "weights", "generator_residual"}
        missing = sorted(required.difference(loaded.files))
        if missing:
            raise ValueError(f"Generator cache missing required arrays {missing}; keys={loaded.files}")
        points = np.asarray(loaded["points"], dtype=float)
        weights = np.asarray(loaded["weights"], dtype=float)
        generator_residual = np.asarray(loaded["generator_residual"], dtype=float)
    if points.ndim != 2 or points.shape[1] != 3 or weights.shape != (len(points),):
        raise ValueError(f"Unexpected generator grid shapes {points.shape} {weights.shape}")
    if generator_residual.shape != (len(points), 3):
        raise ValueError(f"Unexpected generator residual shape {generator_residual.shape}")
    if not np.all(np.isfinite(points)) or not np.all(np.isfinite(weights)) or np.any(weights <= 0.0):
        raise ValueError("Generator grid contains invalid points or weights")
    if not np.all(np.isfinite(generator_residual)):
        raise ValueError("Generator residual contains nonfinite values")

    drift_report = json.loads(DRIFT_PATH.read_text(encoding="utf-8"))
    if drift_report.get("status") != "completed" or not drift_report.get("selected", {}).get("feasible", False):
        raise ValueError("Frozen 1800 drift report is not completed and feasible")
    inputs = drift_report["inputs"]
    first_center = np.asarray(inputs["first_patch_center"], dtype=float)
    first_widths = np.asarray(inputs["first_patch_widths"], dtype=float)
    second_center = np.asarray(inputs["second_patch_center"], dtype=float)
    second_widths = np.asarray(inputs["second_patch_widths"], dtype=float)
    carrier = np.asarray(json.loads((ROOT / "full_wave_frozen_cache.json").read_text(encoding="utf-8"))["inputs"]["wave"]["carrier"], dtype=float)

    pressure_1 = _pressure_patch_design(points, first_center, first_widths, carrier)
    pressure_2 = _pressure_patch_design(points, second_center, second_widths, carrier)
    design = np.column_stack((pressure_1, pressure_2))
    if design.shape != (3 * len(points), 90):
        raise ValueError(f"Unexpected joint pressure design shape {design.shape}")
    row_weight = np.repeat(np.sqrt(weights), 3)
    weighted_design = design * row_weight[:, None]
    weighted_rhs = -generator_residual.reshape(-1) * row_weight
    column_scales = np.maximum(np.linalg.norm(weighted_design, axis=0), 1.0e-30)
    normalized_design = weighted_design / column_scales[None, :]
    singular_values = np.linalg.svd(normalized_design, compute_uv=False)
    cutoff = RIDGE_RELATIVE * max(float(singular_values[0]), 1.0e-300)
    rank = int(np.sum(singular_values > cutoff))
    normalized_coefficients, _, _, _ = np.linalg.lstsq(normalized_design, weighted_rhs, rcond=RIDGE_RELATIVE)
    coefficients = normalized_coefficients / column_scales
    correction = (design @ coefficients).reshape(-1, 3)
    corrected = generator_residual + correction
    before = _metric(generator_residual, weights, points)
    after = _metric(corrected, weights, points)
    weighted_after = corrected.reshape(-1) * row_weight
    normal_residual = normalized_design.T @ weighted_after
    rhs_norm = max(float(np.linalg.norm(weighted_rhs)), 1.0e-300)
    design_norm = max(float(np.linalg.norm(normalized_design)), 1.0e-300)
    projection_orthogonality = {
        "normalized_design_transpose_weighted_residual_l2": float(np.linalg.norm(normal_residual)),
        "relative_to_weighted_rhs": float(np.linalg.norm(normal_residual) / rhs_norm),
        "relative_to_design_norm_times_residual": float(np.linalg.norm(normal_residual) / (design_norm * max(float(np.linalg.norm(weighted_after)), 1.0e-300))),
    }
    source_paths = {
        "script": Path(__file__),
        "generator_report": report_path,
        "generator_cache": cache_path,
        "generator_script": GENERATOR_SCRIPT_PATH,
        "drift_report": DRIFT_PATH,
        "drift_source": DRIFT_SOURCE_PATH,
        "local_basis_source": LOCAL_BASIS_SOURCE_PATH,
    }
    report = {
        "status": "completed",
        "accepted": False,
        "pde_validated": False,
        "scale_recursion_established": False,
        "constraints_maintained": False,
        "scope": "Weighted pressure-only projection of the frozen scale-generator momentum defect. Velocity, generator residual, cones, full-domain constraints, PDE, forcing, and recursive-scale acceptance are outside this diagnostic.",
        "source": Path(__file__).name,
        "source_sha256": _sha256(Path(__file__)),
        "sources": {name: {"path": path.name, "sha256": _sha256(path)} for name, path in source_paths.items()},
        "inputs": {
            "point_count": int(len(points)),
            "physical_volume": float(np.sum(weights)),
            "grid_array_sha256": _array_sha256(points, weights),
            "generator_residual_array_sha256": _array_sha256(generator_residual),
            "patch_count": 2,
            "pressure_columns_per_patch": PRESSURE_COLUMN_COUNT_PER_PATCH,
            "pressure_column_count": 90,
            "degree": DEGREE,
            "modes": list(MODES),
            "first_patch_center": first_center.tolist(),
            "first_patch_widths": first_widths.tolist(),
            "second_patch_center": second_center.tolist(),
            "second_patch_widths": second_widths.tolist(),
            "carrier": carrier.tolist(),
            "weighted_least_squares": True,
            "column_normalization": "weighted L2 norm",
            "svd_relative_cutoff": RIDGE_RELATIVE,
            "generator_report_status": generator_report.get("status"),
        },
        "pressure_layout": {
            "per_patch": _pressure_layout(),
            "joint_order": "first patch 45 pressure columns, then second patch 45 pressure columns",
            "velocity_columns_changed": False,
            "generator_values_changed": False,
        },
        "before": before,
        "after": after,
        "change": {
            "volume_L2": float(after["volume_L2"] - before["volume_L2"]),
            "volume_L2_relative": float((after["volume_L2"] - before["volume_L2"]) / max(before["volume_L2"], 1.0e-300)),
            "max_norm": float(after["max_norm"] - before["max_norm"]),
            "max_norm_relative": float((after["max_norm"] - before["max_norm"]) / max(before["max_norm"], 1.0e-300)),
        },
        "fit": {
            "design_shape": list(design.shape),
            "weighted_design_shape": list(weighted_design.shape),
            "normalized_singular_values": singular_values.tolist(),
            "normalized_rank": rank,
            "column_scale_min": float(np.min(column_scales)),
            "column_scale_max": float(np.max(column_scales)),
            "coefficient_norm": float(np.linalg.norm(coefficients)),
            "coefficients": coefficients.tolist(),
            "first_patch_coefficients": coefficients[:45].tolist(),
            "second_patch_coefficients": coefficients[45:].tolist(),
            "correction_metric": _metric(correction, weights, points),
            "weighted_residual_orthogonality": projection_orthogonality,
        },
        "elapsed_seconds": float(time.perf_counter() - started),
    }
    output_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": report["status"],
        "before_l2": before["volume_L2"],
        "after_l2": after["volume_L2"],
        "before_max": before["max_norm"],
        "after_max": after["max_norm"],
        "rank": rank,
        "orthogonality_relative": projection_orthogonality["relative_to_weighted_rhs"],
        "elapsed_seconds": report["elapsed_seconds"],
    }), flush=True)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT_PATH)
    args = parser.parse_args()
    run(args.output)
