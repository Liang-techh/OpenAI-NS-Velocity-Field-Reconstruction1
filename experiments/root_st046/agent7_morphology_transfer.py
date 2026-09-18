"""Target-free 3-D vorticity morphology transfer audit for ST046-A.

This is Constrained Agent 7 representation/visualization evidence only.  It
reconstructs the frozen ST045-H parent and frozen ST046-A child, then measures
them under one preregistered Cartesian vorticity contract.  It does not fit a
basis, select a visual target, alter pressure/forcing, rerun held-out momentum,
or promote a candidate.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Iterable

import numpy as np

import replay_st046

TASK_ID = "CR003-ST046-MORPHOLOGY-TRANSFER-SCREEN-058"
DEFAULT_TIMES = (0.25, 0.50, 0.75)
DEFAULT_GRID_SIZES = (25, 33)
BOX_HALF_WIDTH = 1.95
RADIAL_SUPPORT = 2.0
AXIAL_SUPPORT = 2.0
OUTER_COLLAR_FRACTION = 0.75
TRUTH_BOUNDARY = {
    "velocity_changed": False,
    "candidate_artifact_changed": False,
    "new_basis_added": False,
    "pressure_or_force_changed": False,
    "held_out_pde_residual_recomputed": False,
    "public_image_used": False,
    "visual_score_defined": False,
    "visualization_ready": False,
    "visual_correspondence_verified": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
}


class _FrozenField:
    def __init__(self, family, raw, ident: str):
        self.family = family
        self.raw = np.asarray(raw, dtype=float)
        self.ident = str(ident)

    def at_points(self, points, time: float) -> np.ndarray:
        points = np.asarray(points, dtype=float)
        if points.ndim != 2 or points.shape[1] != 3:
            raise ValueError("points must have shape (n,3)")
        velocity, _pressure = self.family.fields(self.raw, points, float(time))
        velocity = np.asarray(velocity, dtype=float)
        if velocity.shape != points.shape or not np.all(np.isfinite(velocity)):
            raise RuntimeError("unexpected/nonfinite reconstructed velocity")
        return velocity


def _weighted_quantile(values: np.ndarray, weights: np.ndarray, quantile: float) -> float:
    values = np.asarray(values, dtype=float).ravel()
    weights = np.asarray(weights, dtype=float).ravel()
    quantile = float(quantile)
    if values.shape != weights.shape or not 0.0 <= quantile <= 1.0:
        raise ValueError("invalid weighted quantile inputs")
    if not np.all(np.isfinite(values)) or not np.all(np.isfinite(weights)) or np.any(weights < 0.0):
        raise ValueError("weighted quantile inputs must be finite with nonnegative weights")
    order = np.argsort(values)
    values = values[order]
    weights = weights[order]
    total = float(np.sum(weights))
    if total <= np.finfo(float).tiny:
        raise ValueError("weighted quantile requires positive total weight")
    cumulative = np.cumsum(weights)
    index = int(np.searchsorted(cumulative, quantile * total, side="left"))
    return float(values[min(index, len(values) - 1)])


def _morphology(field: _FrozenField, *, time: float, grid_size: int) -> dict[str, float]:
    grid_size = int(grid_size)
    if grid_size < 17 or grid_size % 2 == 0:
        raise ValueError("grid_size must be odd and >=17")
    axis = np.linspace(-BOX_HALF_WIDTH, BOX_HALF_WIDTH, grid_size)
    spacing = float(axis[1] - axis[0])
    x, y, z = np.meshgrid(axis, axis, axis, indexing="ij")
    points = np.column_stack((x.ravel(), y.ravel(), z.ravel()))
    velocity_flat = field.at_points(points, float(time))
    velocity = velocity_flat.reshape(grid_size, grid_size, grid_size, 3)

    u, v, w = velocity[..., 0], velocity[..., 1], velocity[..., 2]
    two_h = 2.0 * spacing
    dw_dy = (w[1:-1, 2:, 1:-1] - w[1:-1, :-2, 1:-1]) / two_h
    dv_dz = (v[1:-1, 1:-1, 2:] - v[1:-1, 1:-1, :-2]) / two_h
    du_dz = (u[1:-1, 1:-1, 2:] - u[1:-1, 1:-1, :-2]) / two_h
    dw_dx = (w[2:, 1:-1, 1:-1] - w[:-2, 1:-1, 1:-1]) / two_h
    dv_dx = (v[2:, 1:-1, 1:-1] - v[:-2, 1:-1, 1:-1]) / two_h
    du_dy = (u[1:-1, 2:, 1:-1] - u[1:-1, :-2, 1:-1]) / two_h
    omega_x = dw_dy - dv_dz
    omega_y = du_dz - dw_dx
    omega_z = dv_dx - du_dy
    omega_sq = omega_x * omega_x + omega_y * omega_y + omega_z * omega_z
    if not np.all(np.isfinite(omega_sq)):
        raise RuntimeError("vorticity reconstruction became nonfinite")
    total = float(np.sum(omega_sq))
    if total <= np.finfo(float).tiny:
        raise RuntimeError("vorticity is numerically inactive")

    xi = x[1:-1, 1:-1, 1:-1]
    yi = y[1:-1, 1:-1, 1:-1]
    zi = z[1:-1, 1:-1, 1:-1]
    radius_sq = xi * xi + yi * yi
    abs_z = np.abs(zi)
    radial_rms = float(np.sqrt(np.sum(radius_sq * omega_sq) / total))
    axial_rms = float(np.sqrt(np.sum(zi * zi * omega_sq) / total))
    q90 = _weighted_quantile(abs_z, omega_sq, 0.90)
    q99 = _weighted_quantile(abs_z, omega_sq, 0.99)
    collar = (
        (np.sqrt(radius_sq) >= OUTER_COLLAR_FRACTION * RADIAL_SUPPORT)
        | (abs_z >= OUTER_COLLAR_FRACTION * AXIAL_SUPPORT)
    )
    speed_sq = np.sum(velocity_flat * velocity_flat, axis=1)
    return {
        "radial_rms_over_support": radial_rms / RADIAL_SUPPORT,
        "axial_rms_over_support": axial_rms / AXIAL_SUPPORT,
        "aspect_z_over_r": axial_rms / max(radial_rms, 1.0e-300),
        "axial_q90_over_support": q90 / AXIAL_SUPPORT,
        "axial_q99_over_support": q99 / AXIAL_SUPPORT,
        "outer_half_enstrophy_fraction": float(np.sum(omega_sq[abs_z >= 0.50 * AXIAL_SUPPORT]) / total),
        "outer_065_enstrophy_fraction": float(np.sum(omega_sq[abs_z >= 0.65 * AXIAL_SUPPORT]) / total),
        "outer_075_enstrophy_fraction": float(np.sum(omega_sq[abs_z >= 0.75 * AXIAL_SUPPORT]) / total),
        "outer_support_collar_enstrophy_fraction": float(np.sum(omega_sq[collar]) / total),
        "vorticity_rms": float(np.sqrt(np.mean(omega_sq))),
        "vorticity_max": float(np.sqrt(np.max(omega_sq))),
        "velocity_rms": float(np.sqrt(np.mean(speed_sq))),
    }


def _structure_checks(field: _FrozenField, times: Iterable[float]) -> dict[str, object]:
    support_points = np.array(
        [[2.01, 0.0, 0.0], [-2.01, 0.0, 0.0], [0.4, 0.2, 2.01], [0.4, -0.2, -2.01]],
        dtype=float,
    )
    support_max = max(float(np.max(np.abs(field.at_points(support_points, time)))) for time in times)
    rows = []
    for time in times:
        for radius in (0.05, 0.10, 0.20):
            for z in (-0.10, 0.10):
                point = np.array([[radius, 0.0, z]], dtype=float)
                ux, uy, uz = field.at_points(point, float(time))[0]
                rows.append(bool(ux < 0.0 and uy > 0.0 and z * uz > 0.0))
    return {
        "outside_support_max_abs_velocity": support_max,
        "representative_inward_swirl_bipolar_signs_all_pass": bool(all(rows)),
        "representative_sign_checks": len(rows),
    }


def _delta(child: dict[str, float], parent: dict[str, float]) -> dict[str, float]:
    out = {}
    for key in parent:
        out[key + "_delta"] = float(child[key] - parent[key])
        out[key + "_relative"] = float(child[key] / max(abs(parent[key]), 1.0e-300) - 1.0)
    return out


def audit_morphology_transfer(
    *,
    times: Iterable[float] = DEFAULT_TIMES,
    grid_sizes: Iterable[int] = DEFAULT_GRID_SIZES,
) -> dict[str, object]:
    times = tuple(float(value) for value in times)
    grid_sizes = tuple(int(value) for value in grid_sizes)
    if not times or any(not np.isfinite(value) or value < 0.25 or value > 0.75 for value in times):
        raise ValueError("times must lie in [0.25,0.75]")
    if len(grid_sizes) != 2 or not 17 <= grid_sizes[0] < grid_sizes[1] or any(value % 2 == 0 for value in grid_sizes):
        raise ValueError("grid_sizes must contain two increasing odd values >=17")

    parent_family, parent_raw = replay_st046.previous_load("ST045-H")
    child_family, child_raw = replay_st046.reconstruct("ST046-A")
    parent = _FrozenField(parent_family, parent_raw, "ST045-H")
    child = _FrozenField(child_family, child_raw, "ST046-A")

    rows = []
    by_grid = {}
    for grid_size in grid_sizes:
        grid_rows = []
        for time in times:
            parent_metrics = _morphology(parent, time=time, grid_size=grid_size)
            child_metrics = _morphology(child, time=time, grid_size=grid_size)
            grid_rows.append(
                {
                    "time": time,
                    "parent": parent_metrics,
                    "child": child_metrics,
                    "child_minus_parent": _delta(child_metrics, parent_metrics),
                }
            )
        by_grid[grid_size] = grid_rows
        rows.extend(grid_rows)

    coarse = by_grid[grid_sizes[0]]
    fine = by_grid[grid_sizes[1]]
    refinement = {}
    for ident in ("parent", "child"):
        keys = tuple(coarse[0][ident].keys())
        refinement[ident] = {
            key: max(
                abs(coarse_row[ident][key] - fine_row[ident][key])
                / max(abs(fine_row[ident][key]), 1.0e-12)
                for coarse_row, fine_row in zip(coarse, fine)
            )
            for key in keys
        }

    fine_deltas = [row["child_minus_parent"] for row in fine]
    axial_rel = [row["axial_rms_over_support_relative"] for row in fine_deltas]
    radial_rel = [row["radial_rms_over_support_relative"] for row in fine_deltas]
    q90_delta = [row["axial_q90_over_support_delta"] for row in fine_deltas]
    q99_delta = [row["axial_q99_over_support_delta"] for row in fine_deltas]
    collar_ratio = [
        row["child"]["outer_support_collar_enstrophy_fraction"]
        / max(row["parent"]["outer_support_collar_enstrophy_fraction"], 1.0e-300)
        for row in fine
    ]

    morphology_extension_signal = bool(
        all(value >= 0.0 for value in axial_rel)
        and any(value > 1.0e-12 for value in q90_delta + q99_delta)
    )
    low_radial_collar_cost = bool(max(radial_rel) <= 0.05 and max(collar_ratio) <= 2.0)
    clean_transfer = bool(morphology_extension_signal and low_radial_collar_cost)

    if clean_transfer:
        routing = (
            "ST046-A already carries a clean axial morphology gain under this frozen target-free contract. "
            "Use ST046-A as the next visualization backbone before adding the Piola warp from PR #367."
        )
    else:
        routing = (
            "ST046-A's residual improvement does not by itself satisfy the preregistered clean axial-morphology "
            "transfer criterion. Keep ST046-A as the PDE/residual backbone and retain only the single Piola axial-"
            "warp degree from PR #367 as the next minimal visualization-geometry handoff; do not resume basis growth."
        )

    recipes = json.loads((Path(__file__).resolve().parent / "recipes.json").read_text())
    return {
        "task_id": TASK_ID,
        "parent_id": "ST045-H",
        "child_id": "ST046-A",
        "child_original_raw_sha256": recipes["ST046-A"].get("original_candidate_sha256"),
        "contract": {
            "times": list(times),
            "grid_sizes": list(grid_sizes),
            "box_half_width": BOX_HALF_WIDTH,
            "radial_support": RADIAL_SUPPORT,
            "axial_support": AXIAL_SUPPORT,
            "outer_collar_fraction_of_support": OUTER_COLLAR_FRACTION,
            "classification": "autonomous target-free visualization diagnostic",
            "fit_or_parameter_selection": False,
        },
        "rows_by_grid": {str(key): value for key, value in by_grid.items()},
        "refinement_max_relative_change": refinement,
        "parent_structure": _structure_checks(parent, times),
        "child_structure": _structure_checks(child, times),
        "fine_grid_summary": {
            "axial_rms_relative_changes": axial_rel,
            "radial_rms_relative_changes": radial_rel,
            "q90_deltas_over_support": q90_delta,
            "q99_deltas_over_support": q99_delta,
            "collar_enstrophy_ratios": collar_ratio,
            "morphology_extension_signal": morphology_extension_signal,
            "low_radial_collar_cost": low_radial_collar_cost,
            "clean_axial_morphology_transfer": clean_transfer,
        },
        "routing": routing,
        "truth_boundary": dict(TRUTH_BOUNDARY),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", default="artifacts/st046_agent7_morphology_transfer/report.json")
    args = parser.parse_args()
    report = audit_morphology_transfer()
    path = Path(args.output)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
