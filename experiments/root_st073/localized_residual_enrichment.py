"""Single localized exact-curl residual enrichment on the frozen wave cache.

The refined 44,400-point cache already contains the frozen base residual and
264-column tangent response.  This diagnostic keeps those old controls fixed
and fits one compact degree-2 patch at the largest current residual point.
Modes 0, 1, and 2 contribute 36, 72, and 72 real columns respectively.  The
patch is an exact-curl velocity-derivative block plus pressure-gradient block,
so it changes only the reference momentum residual in this linear diagnostic.

No mean replay, moment/cone check, endpoint-shape check, PDE check, time step,
or recursive-scale acceptance is performed here.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
if str(ROOT) not in __import__("sys").path:
    __import__("sys").path.insert(0, str(ROOT))

from supported_fourier_basis import basis_data  # noqa: E402
from wave_residual_harmonics import budget as harmonic_budget  # noqa: E402


CACHE_PATH = ROOT / "refined_wave_momentum_cache.npz"
CACHE_REPORT_PATH = ROOT / "refined_wave_momentum_cache.json"
OUTPUT_PATH = ROOT / "localized_residual_enrichment.json"

PATCH_CENTER = np.asarray((9.4e-4, -5.49e-5), dtype=float)
PATCH_WIDTHS = np.asarray((2.5e-4, 2.0e-4), dtype=float)
DEGREE = 2
MODES = (0, 1, 2)
RIDGE_LAMBDA = 1.0e-8
SVD_RANK_RELATIVE_CUTOFF = 1.0e-12
ANGLES = 12


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _metric(residual: np.ndarray, weights: np.ndarray) -> dict:
    residual = np.asarray(residual, dtype=float)
    weights = np.asarray(weights, dtype=float)
    magnitudes = np.linalg.norm(residual, axis=1)
    weighted_square = float(np.sum(weights * np.sum(residual * residual, axis=1)))
    return {
        "point_count": int(len(residual)),
        "momentum_max": float(np.max(magnitudes)),
        "momentum_volume_L2": float(np.sqrt(max(weighted_square, 0.0))),
        "momentum_volume_RMS": float(np.sqrt(max(weighted_square, 0.0) / np.sum(weights))),
        "physical_volume": float(np.sum(weights)),
    }


def _max_residual_detail(residual: np.ndarray, points: np.ndarray, weights: np.ndarray) -> dict:
    magnitudes = np.linalg.norm(np.asarray(residual), axis=1)
    index = int(np.argmax(magnitudes))
    point = np.asarray(points[index], dtype=float)
    return {
        "index": index,
        "point": point.tolist(),
        "radius": float(np.hypot(point[0], point[1])),
        "z": float(point[2]),
        "weight": float(weights[index]),
        "vector": np.asarray(residual[index], dtype=float).tolist(),
        "norm": float(magnitudes[index]),
    }


def _mode_block_columns(
    points: np.ndarray,
    center: np.ndarray,
    widths: np.ndarray,
    mode: int,
    degree: int,
    carrier: np.ndarray,
) -> np.ndarray:
    """Return interleaved [Re, -Im] velocity and pressure-gradient columns."""

    velocity, _, pressure_gradient = basis_data(
        points, center, widths, mode, degree, carrier
    )
    blocks = []
    for tensor in (velocity, pressure_gradient):
        for index in range(tensor.shape[-1]):
            blocks.append(tensor[:, :, index].real.reshape(-1))
            blocks.append(-tensor[:, :, index].imag.reshape(-1))
    return np.stack(blocks, axis=1)


def _patch_design(points: np.ndarray, center: np.ndarray, widths: np.ndarray, carrier: np.ndarray):
    mode0_full = _mode_block_columns(points, center, widths, 0, DEGREE, np.zeros(2))
    mode0 = mode0_full[:, ::2]
    mode1 = _mode_block_columns(points, center, widths, 1, DEGREE, carrier)
    mode2 = _mode_block_columns(points, center, widths, 2, DEGREE, 2.0 * carrier)
    design = np.column_stack((mode0, mode1, mode2))
    if design.shape[1] != 180:
        raise ValueError(f"Expected 180 local columns, got {design.shape}")
    layout = []
    q = (DEGREE + 1) ** 2
    for mode, components in ((0, ("real",)), (1, ("real", "imag")), (2, ("real", "imag"))):
        for kind, count in (("velocity", 3 * q), ("pressure_gradient", q)):
            for index in range(count):
                for component in components:
                    layout.append({
                        "mode": mode,
                        "kind": kind,
                        "component": component,
                        "index": index,
                    })
    if len(layout) != 180:
        raise ValueError(f"Unexpected local layout length {len(layout)}")
    return design, layout


def _support_report(cache_center: np.ndarray, cache_widths: np.ndarray) -> dict:
    original_lower = cache_center - cache_widths
    original_upper = cache_center + cache_widths
    patch_lower = PATCH_CENTER - PATCH_WIDTHS
    patch_upper = PATCH_CENTER + PATCH_WIDTHS
    lower_margin = patch_lower - original_lower
    upper_margin = original_upper - patch_upper
    strict_inside = bool(np.all(lower_margin > 0.0) and np.all(upper_margin > 0.0))
    return {
        "coordinates": ["radius", "z"],
        "patch_center": PATCH_CENTER.tolist(),
        "patch_widths": PATCH_WIDTHS.tolist(),
        "patch_bounds": {"lower": patch_lower.tolist(), "upper": patch_upper.tolist()},
        "original_wave_center": np.asarray(cache_center, dtype=float).tolist(),
        "original_wave_widths": np.asarray(cache_widths, dtype=float).tolist(),
        "original_wave_bounds": {"lower": original_lower.tolist(), "upper": original_upper.tolist()},
        "strictly_inside_original_support": strict_inside,
        "lower_containment_margin": lower_margin.tolist(),
        "upper_containment_margin": upper_margin.tolist(),
        "patch_outside_original_support": bool(not strict_inside),
        "added_support_outside_original_support": bool(not strict_inside),
    }


def _source_entry(path: Path) -> dict:
    return {"path": path.name, "sha256": _sha256(path)}


def run(output_path=OUTPUT_PATH):
    started = time.perf_counter()
    output_path = Path(output_path)
    cache_report = json.loads(CACHE_REPORT_PATH.read_text(encoding="utf-8"))
    if cache_report.get("status") != "completed":
        raise ValueError("Refined momentum cache report is not completed")
    with np.load(CACHE_PATH, allow_pickle=False) as loaded:
        points = np.asarray(loaded["points"], dtype=float)
        weights = np.asarray(loaded["weights"], dtype=float)
        base_residual = np.asarray(loaded["residual"], dtype=float)
        old_design = np.asarray(loaded["tangent_design"], dtype=float)
        old_control = np.asarray(loaded["tangent_coefficients"], dtype=float)
        cache_center = np.asarray(loaded["center"], dtype=float)
        cache_widths = np.asarray(loaded["widths"], dtype=float)
        carrier = np.asarray(loaded["carrier"], dtype=float)

    if points.shape != (44400, 3) or weights.shape != (44400,):
        raise ValueError(f"Unexpected refined grid shapes {points.shape} {weights.shape}")
    if base_residual.shape != (44400, 3) or old_design.shape != (133200, 264):
        raise ValueError(f"Unexpected refined cache arrays {base_residual.shape} {old_design.shape}")
    if old_control.shape != (264,):
        raise ValueError(f"Expected 264 frozen controls, got {old_control.shape}")
    reported_control = np.asarray(cache_report["inputs"]["tangent_coefficients"], dtype=float)
    if reported_control.shape != (264,) or not np.array_equal(reported_control, old_control):
        raise ValueError("Refined cache control does not match its report")

    # The old candidate is held fixed.  This is the exact frozen residual plus
    # its stored 264-column affine tangent response.
    current_residual = base_residual + (old_design @ old_control).reshape(-1, 3)
    del old_design, base_residual
    before_metric = _metric(current_residual, weights)
    before_peak = _max_residual_detail(current_residual, points, weights)
    before_budget = harmonic_budget(points, weights, current_residual, ANGLES)

    design, layout = _patch_design(points, PATCH_CENTER, PATCH_WIDTHS, carrier)
    sqrt_weights = np.sqrt(weights)
    row_weight = np.repeat(sqrt_weights, 3)
    weighted_design = design * row_weight[:, None]
    del design
    column_scales = np.linalg.norm(weighted_design, axis=0)
    if np.any(column_scales <= 1.0e-30):
        raise ValueError("Local patch contains a zero weighted column")
    weighted_design /= column_scales[None, :]
    weighted_rhs = -current_residual.reshape(-1) * row_weight

    # Singular values document the normalized design rank.  The fit itself is
    # the requested normalized ridge solve at lambda = 1e-8.
    singular_values = np.linalg.svd(weighted_design, compute_uv=False)
    cutoff = SVD_RANK_RELATIVE_CUTOFF * max(float(singular_values[0]), 1.0e-300)
    rank = int(np.sum(singular_values > cutoff))
    gram = weighted_design.T @ weighted_design
    rhs = weighted_design.T @ weighted_rhs
    normalized_coefficients = np.linalg.solve(
        gram + RIDGE_LAMBDA * np.eye(weighted_design.shape[1]), rhs
    )
    patch_coefficients = normalized_coefficients / column_scales
    weighted_correction = weighted_design @ normalized_coefficients
    correction = (weighted_correction / row_weight).reshape(-1, 3)
    corrected_residual = current_residual + correction
    after_metric = _metric(corrected_residual, weights)
    after_peak = _max_residual_detail(corrected_residual, points, weights)
    after_budget = harmonic_budget(points, weights, corrected_residual, ANGLES)

    linked_sources = {
        "refined_momentum_cache": _source_entry(CACHE_PATH),
        "refined_momentum_cache_report": _source_entry(CACHE_REPORT_PATH),
    }
    for key, entry in cache_report.get("sources", {}).items():
        source_path = ROOT / entry["path"]
        if source_path.exists():
            linked_sources[key] = _source_entry(source_path)

    report = {
        "status": "completed",
        "accepted": False,
        "pde_validated": False,
        "scale_recursion_established": False,
        "constraints_maintained": False,
        "scope": (
            "One localized degree-2 exact-curl velocity-derivative and pressure-"
            "gradient patch fitted against the frozen 44400-point residual. The "
            "264 previous controls are held fixed. No mean finite-difference replay, "
            "moments, cones, endpoint shape, PDE, time step, or recursion claim."
        ),
        "sources": linked_sources,
        "inputs": {
            "point_count": int(len(points)),
            "old_design_shape": [133200, 264],
            "old_control_count": 264,
            "new_control_count": 180,
            "degree": DEGREE,
            "modes": list(MODES),
            "carrier_mode1": carrier.tolist(),
            "carrier_mode2": (2.0 * carrier).tolist(),
            "angles": ANGLES,
            "old_controls_held_fixed": True,
            "initial_velocity_unchanged_at_reference": True,
            "new_layout": layout,
        },
        "support_containment": _support_report(cache_center, cache_widths),
        "target_peak_before_fit": before_peak,
        "fit": {
            "method": "weighted least squares of current residual onto new columns only",
            "weighted_row_definition": "sqrt(volume weight) repeated for 3 Cartesian components",
            "column_normalization": "weighted L2 norm per column",
            "ridge_lambda_normalized": RIDGE_LAMBDA,
            "svd_relative_rank_cutoff": SVD_RANK_RELATIVE_CUTOFF,
            "rank": rank,
            "column_count": int(weighted_design.shape[1]),
            "singular_values": singular_values.tolist(),
            "singular_value_min": float(np.min(singular_values)),
            "singular_value_max": float(np.max(singular_values)),
            "condition_number": float(np.max(singular_values) / np.min(singular_values)),
            "column_scales": column_scales.tolist(),
            "coefficients": patch_coefficients.tolist(),
            "coefficient_l2": float(np.linalg.norm(patch_coefficients)),
            "coefficient_max_abs": float(np.max(np.abs(patch_coefficients))),
        },
        "momentum_before": before_metric,
        "momentum_after": after_metric,
        "momentum_change": {
            "max_absolute": float(after_metric["momentum_max"] - before_metric["momentum_max"]),
            "max_fractional": float(after_metric["momentum_max"] / before_metric["momentum_max"] - 1.0),
            "L2_absolute": float(after_metric["momentum_volume_L2"] - before_metric["momentum_volume_L2"]),
            "L2_fractional": float(after_metric["momentum_volume_L2"] / before_metric["momentum_volume_L2"] - 1.0),
        },
        "peak_after_fit": after_peak,
        "harmonic_budget_before": before_budget,
        "harmonic_budget_after": after_budget,
        "harmonic_budget_scope": "Sampled cylindrical Fourier decomposition on the frozen 12-angle grid; no continuum certificate.",
        "elapsed_seconds": float(time.perf_counter() - started),
    }
    output_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": report["status"],
        "rank": rank,
        "before_max": before_metric["momentum_max"],
        "after_max": after_metric["momentum_max"],
        "before_L2": before_metric["momentum_volume_L2"],
        "after_L2": after_metric["momentum_volume_L2"],
        "max_fractional": report["momentum_change"]["max_fractional"],
        "L2_fractional": report["momentum_change"]["L2_fractional"],
        "strictly_inside_original_support": report["support_containment"]["strictly_inside_original_support"],
        "elapsed_seconds": report["elapsed_seconds"],
    }), flush=True)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT_PATH)
    args = parser.parse_args()
    run(args.output)
