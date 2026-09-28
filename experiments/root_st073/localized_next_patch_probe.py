"""One additional localized degree-2 residual enrichment probe.

The fixed 264-control enriched candidate and the assembled-feasible 180-control
first localized patch are reconstructed from frozen arrays.  A second compact
degree-2 mode-0/1/2 patch is then fit by normalized weighted ridge least
squares only.  This is a basis-enrichment diagnostic; no moments, cones,
endpoint shape, finite-difference replay, PDE, or adoption decision is made.
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

from localized_residual_enrichment import (  # noqa: E402
    _max_residual_detail,
    _metric,
    _patch_design,
)
from wave_residual_harmonics import budget as harmonic_budget  # noqa: E402


CACHE_PATH = ROOT / "refined_wave_momentum_cache.npz"
CACHE_REPORT_PATH = ROOT / "refined_wave_momentum_cache.json"
FIRST_RESIDUAL_REPORT_PATH = ROOT / "localized_residual_enrichment.json"
FIRST_CONSTRAINT_REPORT_PATH = ROOT / "localized_constraint_rows.json"
FIRST_FIT_REPORT_PATH = ROOT / "localized_constrained_tangent.json"
CANDIDATE_PATH = ROOT / "enriched_mean_endpoint_tangent.json"
GEOMETRY_PATH = ROOT / "full_wave_frozen_cache.json"
OUTPUT_PATH = ROOT / "localized_next_patch_probe.json"

NEXT_CENTER = np.asarray((0.001180209079, 1.4192077977e-5), dtype=float)
NEXT_WIDTHS = np.asarray((3.0e-4, 2.5e-4), dtype=float)
DEGREE = 2
MODES = (0, 1, 2)
RIDGE_LAMBDA = 1.0e-8
ANGLES = 12


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _source_entry(path: Path) -> dict:
    return {"path": path.name, "sha256": _sha256(path)}


def _support_report(cache_center, cache_widths):
    cache_center = np.asarray(cache_center, dtype=float)
    cache_widths = np.asarray(cache_widths, dtype=float)
    original_lower = cache_center - cache_widths
    original_upper = cache_center + cache_widths
    patch_lower = NEXT_CENTER - NEXT_WIDTHS
    patch_upper = NEXT_CENTER + NEXT_WIDTHS
    lower_margin = patch_lower - original_lower
    upper_margin = original_upper - patch_upper
    strict_inside = bool(np.all(lower_margin > 0.0) and np.all(upper_margin > 0.0))
    return {
        "coordinates": ["radius", "z"],
        "patch_center": NEXT_CENTER.tolist(),
        "patch_widths": NEXT_WIDTHS.tolist(),
        "patch_bounds": {"lower": patch_lower.tolist(), "upper": patch_upper.tolist()},
        "original_wave_center": cache_center.tolist(),
        "original_wave_widths": cache_widths.tolist(),
        "original_wave_bounds": {"lower": original_lower.tolist(), "upper": original_upper.tolist()},
        "strictly_inside_original_support": strict_inside,
        "patch_outside_original_support": bool(not strict_inside),
        "lower_containment_margin": lower_margin.tolist(),
        "upper_containment_margin": upper_margin.tolist(),
    }


def _normalized_patch_location(point, center, widths):
    point = np.asarray(point, dtype=float)
    radius = float(np.hypot(point[0], point[1]))
    return {
        "radius": radius,
        "z": float(point[2]),
        "normalized_radius": float((radius - center[0]) / widths[0]),
        "normalized_z": float((point[2] - center[1]) / widths[1]),
        "distance_to_nearest_edge": float(
            min(1.0 - abs((radius - center[0]) / widths[0]),
                1.0 - abs((point[2] - center[1]) / widths[1]))
        ),
    }


def run(output_path=OUTPUT_PATH):
    started = time.perf_counter()
    output_path = Path(output_path)
    cache_report = json.loads(CACHE_REPORT_PATH.read_text(encoding="utf-8"))
    first_residual_report = json.loads(FIRST_RESIDUAL_REPORT_PATH.read_text(encoding="utf-8"))
    first_constraint_report = json.loads(FIRST_CONSTRAINT_REPORT_PATH.read_text(encoding="utf-8"))
    first_fit_report = json.loads(FIRST_FIT_REPORT_PATH.read_text(encoding="utf-8"))
    candidate = json.loads(CANDIDATE_PATH.read_text(encoding="utf-8"))
    geometry = json.loads(GEOMETRY_PATH.read_text(encoding="utf-8"))["inputs"]["wave"]
    if cache_report.get("status") != "completed":
        raise ValueError("Refined cache is not completed")
    if first_fit_report.get("status") != "completed" or not first_fit_report.get("assembled_feasible"):
        raise ValueError("First localized constrained fit is not assembled-feasible")

    with np.load(CACHE_PATH, allow_pickle=False) as loaded:
        points = np.asarray(loaded["points"], dtype=float)
        weights = np.asarray(loaded["weights"], dtype=float)
        base_residual = np.asarray(loaded["residual"], dtype=float)
        old_design = np.asarray(loaded["tangent_design"], dtype=float)
        old_control = np.asarray(loaded["tangent_coefficients"], dtype=float)
        cache_wave = np.asarray(loaded["coefficients_original"], dtype=float)
        cache_center = np.asarray(loaded["center"], dtype=float)
        cache_widths = np.asarray(loaded["widths"], dtype=float)
        carrier = np.asarray(loaded["carrier"], dtype=float)
    wave_packed = np.asarray(candidate["selected"]["coefficients_original"], dtype=float)
    if not np.array_equal(cache_wave, wave_packed):
        raise ValueError("Refined cache and parent candidate initial waves differ")
    local_coefficients = np.asarray(first_fit_report["selected"]["local_patch_coefficients"], dtype=float)
    if local_coefficients.shape != (180,):
        raise ValueError(f"Expected first localized 180-vector, got {local_coefficients.shape}")
    if points.shape != (44400, 3) or weights.shape != (44400,):
        raise ValueError("Unexpected refined grid dimensions")
    if old_design.shape != (133200, 264) or old_control.shape != (264,):
        raise ValueError("Unexpected frozen parent design/control dimensions")

    first_center = np.asarray(first_constraint_report["support_containment"]["patch_center"], dtype=float)
    first_widths = np.asarray(first_constraint_report["support_containment"]["patch_widths"], dtype=float)
    first_design, first_layout = _patch_design(points, first_center, first_widths, carrier)
    if first_design.shape != (133200, 180) or len(first_layout) != 180:
        raise ValueError("Unexpected first localized response dimensions")
    parent_residual = base_residual + (old_design @ old_control).reshape(-1, 3)
    current_residual = parent_residual + (first_design @ local_coefficients).reshape(-1, 3)
    del base_residual, old_design, first_design
    parent_metric = _metric(parent_residual, weights)
    current_metric = _metric(current_residual, weights)
    parent_peak = _max_residual_detail(parent_residual, points, weights)
    current_peak = _max_residual_detail(current_residual, points, weights)
    parent_budget = harmonic_budget(points, weights, parent_residual, ANGLES)
    current_budget = harmonic_budget(points, weights, current_residual, ANGLES)

    next_design, next_layout = _patch_design(points, NEXT_CENTER, NEXT_WIDTHS, carrier)
    sqrt_weights = np.sqrt(weights)
    row_weight = np.repeat(sqrt_weights, 3)
    weighted_design = next_design * row_weight[:, None]
    del next_design
    column_scales = np.maximum(np.linalg.norm(weighted_design, axis=0), 1.0e-30)
    normalized_design = weighted_design / column_scales[None, :]
    weighted_rhs = -current_residual.reshape(-1) * row_weight
    singular_values = np.linalg.svd(normalized_design, compute_uv=False)
    gram = normalized_design.T @ normalized_design
    rhs = normalized_design.T @ weighted_rhs
    normalized_coefficients = np.linalg.solve(
        gram + RIDGE_LAMBDA * np.eye(normalized_design.shape[1]), rhs
    )
    next_coefficients = normalized_coefficients / column_scales
    weighted_correction = normalized_design @ normalized_coefficients
    correction = (weighted_correction / row_weight).reshape(-1, 3)
    enriched_residual = current_residual + correction
    enriched_metric = _metric(enriched_residual, weights)
    enriched_peak = _max_residual_detail(enriched_residual, points, weights)
    enriched_budget = harmonic_budget(points, weights, enriched_residual, ANGLES)

    source_paths = {
        "refined_cache": CACHE_PATH,
        "refined_cache_report": CACHE_REPORT_PATH,
        "parent_candidate": CANDIDATE_PATH,
        "first_residual_enrichment": FIRST_RESIDUAL_REPORT_PATH,
        "first_constraint_rows": FIRST_CONSTRAINT_REPORT_PATH,
        "first_constrained_fit": FIRST_FIT_REPORT_PATH,
        "geometry": GEOMETRY_PATH,
    }
    original_support = _support_report(cache_center, cache_widths)
    report = {
        "status": "completed",
        "accepted": False,
        "pde_validated": False,
        "scale_recursion_established": False,
        "constraints_maintained": False,
        "source": FIRST_FIT_REPORT_PATH.name,
        "source_sha256": _sha256(FIRST_FIT_REPORT_PATH),
        "sources": {name: _source_entry(path) for name, path in source_paths.items()},
        "scope": (
            "One additional compact degree-2 mode-0/1/2 weighted ridge probe. "
            "The fixed 264 parent and first assembled-feasible local 180 patch are "
            "reconstructed from frozen cache responses. No mean finite-difference "
            "replay, moment/cone/shape enforcement, PDE, or adoption claim."
        ),
        "inputs": {
            "point_count": int(len(points)),
            "parent_control_count": 264,
            "first_local_control_count": 180,
            "additional_control_count": 180,
            "degree": DEGREE,
            "modes": list(MODES),
            "angles": ANGLES,
            "carrier_mode1": carrier.tolist(),
            "carrier_mode2": (2.0 * carrier).tolist(),
            "ridge_lambda_normalized": RIDGE_LAMBDA,
            "old_parent_controls_fixed": True,
            "first_local_controls_fixed": True,
        },
        "support_containment": original_support,
        "additional_patch": {
            "center": NEXT_CENTER.tolist(),
            "widths": NEXT_WIDTHS.tolist(),
            "bounds": original_support["patch_bounds"],
            "layout": next_layout,
            "velocity_columns": 81,
            "pressure_gradient_columns": 27,
            "total_columns": 180,
        },
        "parent_residual": {
            "fixed_264_metric": parent_metric,
            "fixed_264_peak": parent_peak,
            "fixed_264_harmonic_budget": parent_budget,
        },
        "current_reference_after_first_patch": {
            "metric": current_metric,
            "peak": current_peak,
            "harmonic_budget": current_budget,
            "first_patch_coefficient_norm": float(np.linalg.norm(local_coefficients)),
            "peak_location_relative_to_additional_patch": _normalized_patch_location(
                current_peak["point"], NEXT_CENTER, NEXT_WIDTHS
            ),
        },
        "fit": {
            "method": "weighted least squares of current post-first-patch residual onto additional patch only",
            "weighted_row_definition": "sqrt(volume weight) repeated for 3 Cartesian components",
            "column_normalization": "weighted L2 norm per additional column",
            "ridge_lambda_normalized": RIDGE_LAMBDA,
            "column_count": 180,
            "rank_relative_cutoff": 1.0e-12,
            "rank": int(np.sum(singular_values > 1.0e-12 * max(float(singular_values[0]), 1.0e-300))),
            "singular_value_min": float(np.min(singular_values)),
            "singular_value_max": float(np.max(singular_values)),
            "condition_number": float(np.max(singular_values) / np.min(singular_values)),
            "column_scale_min": float(np.min(column_scales)),
            "column_scale_max": float(np.max(column_scales)),
            "coefficients": next_coefficients.tolist(),
            "coefficient_l2": float(np.linalg.norm(next_coefficients)),
            "coefficient_max_abs": float(np.max(np.abs(next_coefficients))),
        },
        "after_additional_patch": {
            "metric": enriched_metric,
            "peak": enriched_peak,
            "harmonic_budget": enriched_budget,
            "peak_location_relative_to_additional_patch": _normalized_patch_location(
                enriched_peak["point"], NEXT_CENTER, NEXT_WIDTHS
            ),
        },
        "change_from_current_reference": {
            "max_absolute": float(enriched_metric["momentum_max"] - current_metric["momentum_max"]),
            "max_fractional": float(enriched_metric["momentum_max"] / current_metric["momentum_max"] - 1.0),
            "L2_absolute": float(enriched_metric["momentum_volume_L2"] - current_metric["momentum_volume_L2"]),
            "L2_fractional": float(enriched_metric["momentum_volume_L2"] / current_metric["momentum_volume_L2"] - 1.0),
        },
        "angular_diagnosis": {
            "current_mode_L2": [row["volume_L2"] for row in current_budget["modes"][:3]],
            "current_mode_squared_fractions": [row["squared_L2_fraction"] for row in current_budget["modes"][:3]],
            "after_mode_L2": [row["volume_L2"] for row in enriched_budget["modes"][:3]],
            "after_mode_squared_fractions": [row["squared_L2_fraction"] for row in enriched_budget["modes"][:3]],
            "interpretation": (
                "The current peak is inside the additional patch support; compare its "
                "normalized edge distances with the post-fit peak. The Fourier budget "
                "is sampled angular evidence only: mode-0/mode-2 dominance is already "
                "within the represented basis and does not establish missing higher "
                "angular modes. A large local reduction with the peak moving outside "
                "the patch indicates residual spatial structure in those represented "
                "modes. No continuum conclusion is made."
            ),
        },
        "elapsed_seconds": float(time.perf_counter() - started),
    }
    output_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": report["status"],
        "current_l2": current_metric["momentum_volume_L2"],
        "after_l2": enriched_metric["momentum_volume_L2"],
        "current_max": current_metric["momentum_max"],
        "after_max": enriched_metric["momentum_max"],
        "rank": report["fit"]["rank"],
        "current_peak": current_peak["point"],
        "after_peak": enriched_peak["point"],
        "elapsed_seconds": report["elapsed_seconds"],
    }), flush=True)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT_PATH)
    args = parser.parse_args()
    run(args.output)
