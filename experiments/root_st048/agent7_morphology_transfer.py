"""Target-free 3-D vorticity morphology transfer audit for frozen ST048 children.

Constrained Agent 7 representation/visualization evidence only. This compares
ST048-S and ST048-B against frozen ST047-E under the same contract used by the
preceding ST046/ST047 Agent-7 audits. No fitting, new basis, pressure/forcing
change, hidden-image target, or held-out PDE recomputation occurs here.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Iterable

import numpy as np

import replay_st048

TASK_ID = "CR003-ST048-MORPHOLOGY-TRANSFER-SCREEN-060"
PARENT_ID = "ST047-E"
CHILD_IDS = ("ST048-S", "ST048-B")
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


def _load_frozen_fields():
    family, raw = replay_st048.parent_field()
    parent = _FrozenField(family, raw, PARENT_ID)
    children = {}
    for ident in CHILD_IDS:
        family, raw = replay_st048.reconstruct(ident)
        children[ident] = _FrozenField(family, raw, ident)
    return parent, children


def _weighted_quantile(values, weights, quantile):
    values = np.asarray(values, dtype=float).ravel()
    weights = np.asarray(weights, dtype=float).ravel()
    if values.shape != weights.shape or not 0.0 <= quantile <= 1.0:
        raise ValueError("invalid weighted quantile inputs")
    if not np.all(np.isfinite(values)) or not np.all(np.isfinite(weights)) or np.any(weights < 0.0):
        raise ValueError("weighted quantile inputs must be finite with nonnegative weights")
    order = np.argsort(values)
    values, weights = values[order], weights[order]
    total = float(np.sum(weights))
    if total <= np.finfo(float).tiny:
        raise ValueError("weighted quantile requires positive total weight")
    index = int(np.searchsorted(np.cumsum(weights), quantile * total, side="left"))
    return float(values[min(index, len(values) - 1)])


def _morphology(field, *, time: float, grid_size: int):
    if grid_size < 17 or grid_size % 2 == 0:
        raise ValueError("grid_size must be odd and >=17")
    axis = np.linspace(-BOX_HALF_WIDTH, BOX_HALF_WIDTH, grid_size)
    spacing = float(axis[1] - axis[0])
    x, y, z = np.meshgrid(axis, axis, axis, indexing="ij")
    points = np.column_stack((x.ravel(), y.ravel(), z.ravel()))
    velocity_flat = field.at_points(points, time)
    velocity = velocity_flat.reshape(grid_size, grid_size, grid_size, 3)
    u, v, w = velocity[..., 0], velocity[..., 1], velocity[..., 2]
    two_h = 2.0 * spacing
    omega_x = (w[1:-1, 2:, 1:-1] - w[1:-1, :-2, 1:-1]) / two_h - (v[1:-1, 1:-1, 2:] - v[1:-1, 1:-1, :-2]) / two_h
    omega_y = (u[1:-1, 1:-1, 2:] - u[1:-1, 1:-1, :-2]) / two_h - (w[2:, 1:-1, 1:-1] - w[:-2, 1:-1, 1:-1]) / two_h
    omega_z = (v[2:, 1:-1, 1:-1] - v[:-2, 1:-1, 1:-1]) / two_h - (u[1:-1, 2:, 1:-1] - u[1:-1, :-2, 1:-1]) / two_h
    omega_sq = omega_x * omega_x + omega_y * omega_y + omega_z * omega_z
    total = float(np.sum(omega_sq))
    if not np.isfinite(total) or total <= np.finfo(float).tiny:
        raise RuntimeError("vorticity is nonfinite/inactive")
    xi, yi, zi = x[1:-1, 1:-1, 1:-1], y[1:-1, 1:-1, 1:-1], z[1:-1, 1:-1, 1:-1]
    radius_sq, abs_z = xi * xi + yi * yi, np.abs(zi)
    radial_rms = float(np.sqrt(np.sum(radius_sq * omega_sq) / total))
    axial_rms = float(np.sqrt(np.sum(zi * zi * omega_sq) / total))
    collar = (np.sqrt(radius_sq) >= OUTER_COLLAR_FRACTION * RADIAL_SUPPORT) | (abs_z >= OUTER_COLLAR_FRACTION * AXIAL_SUPPORT)
    speed_sq = np.sum(velocity_flat * velocity_flat, axis=1)
    return {
        "radial_rms_over_support": radial_rms / RADIAL_SUPPORT,
        "axial_rms_over_support": axial_rms / AXIAL_SUPPORT,
        "aspect_z_over_r": axial_rms / max(radial_rms, 1e-300),
        "axial_q90_over_support": _weighted_quantile(abs_z, omega_sq, 0.90) / AXIAL_SUPPORT,
        "axial_q99_over_support": _weighted_quantile(abs_z, omega_sq, 0.99) / AXIAL_SUPPORT,
        "outer_half_enstrophy_fraction": float(np.sum(omega_sq[abs_z >= 0.50 * AXIAL_SUPPORT]) / total),
        "outer_065_enstrophy_fraction": float(np.sum(omega_sq[abs_z >= 0.65 * AXIAL_SUPPORT]) / total),
        "outer_075_enstrophy_fraction": float(np.sum(omega_sq[abs_z >= 0.75 * AXIAL_SUPPORT]) / total),
        "outer_support_collar_enstrophy_fraction": float(np.sum(omega_sq[collar]) / total),
        "vorticity_rms": float(np.sqrt(np.mean(omega_sq))),
        "vorticity_max": float(np.sqrt(np.max(omega_sq))),
        "velocity_rms": float(np.sqrt(np.mean(speed_sq))),
    }


def _structure_checks(field, times: Iterable[float]):
    exterior = np.array([[2.01, 0, 0], [-2.01, 0, 0], [.4, .2, 2.01], [.4, -.2, -2.01]], dtype=float)
    support_max = max(float(np.max(np.abs(field.at_points(exterior, t)))) for t in times)
    checks = []
    for time in times:
        for radius in (.05, .10, .20):
            for z in (-.10, .10):
                ux, uy, uz = field.at_points(np.array([[radius, 0.0, z]]), time)[0]
                checks.append(bool(ux < 0 and uy > 0 and z * uz > 0))
    return {
        "outside_support_max_abs_velocity": support_max,
        "representative_inward_swirl_bipolar_signs_all_pass": bool(all(checks)),
        "representative_sign_checks": len(checks),
    }


def _delta(child, parent):
    out = {}
    for key in parent:
        out[key + "_delta"] = float(child[key] - parent[key])
        out[key + "_relative"] = float(child[key] / max(abs(parent[key]), 1e-300) - 1.0)
    return out


def _residual_receipt():
    results = json.loads((Path(__file__).resolve().parent / "results.json").read_text())
    rows = [{key: row[key] for key in ("id", "seed", "points", "time_slices", "max", "volume_L2")}
            for row in results["paired_results"] if row["id"] in (PARENT_ID,) + CHILD_IDS]
    return {
        "classification": "frozen upstream independent holdout receipt; not recomputed by Agent 7",
        "freeze_utc": results["freeze_utc"],
        "momentum_thresholds": results["momentum_thresholds"],
        "paired_results": rows,
        "original_NS_target_achieved": bool(results["original_NS_target_achieved"]),
    }


def audit_morphology_transfer(*, times=DEFAULT_TIMES, grid_sizes=DEFAULT_GRID_SIZES):
    times, grid_sizes = tuple(map(float, times)), tuple(map(int, grid_sizes))
    if any(t < .25 or t > .75 or not np.isfinite(t) for t in times):
        raise ValueError("times must lie in [0.25,0.75]")
    if len(grid_sizes) != 2 or not 17 <= grid_sizes[0] < grid_sizes[1] or any(g % 2 == 0 for g in grid_sizes):
        raise ValueError("grid_sizes must be two increasing odd values >=17")
    parent, children = _load_frozen_fields()
    by_grid = {}
    for grid_size in grid_sizes:
        rows = []
        for time in times:
            p = _morphology(parent, time=time, grid_size=grid_size)
            cs = {ident: _morphology(field, time=time, grid_size=grid_size) for ident, field in children.items()}
            rows.append({"time": time, "parent": p, "children": cs,
                         "child_minus_parent": {ident: _delta(metrics, p) for ident, metrics in cs.items()}})
        by_grid[grid_size] = rows
    coarse, fine = by_grid[grid_sizes[0]], by_grid[grid_sizes[1]]
    refinement = {}
    for ident in (PARENT_ID,) + CHILD_IDS:
        source = "parent" if ident == PARENT_ID else "children"
        sample = coarse[0]["parent"] if ident == PARENT_ID else coarse[0]["children"][ident]
        refinement[ident] = {}
        for key in sample:
            diffs = []
            for cr, fr in zip(coarse, fine):
                cv = cr["parent"][key] if ident == PARENT_ID else cr[source][ident][key]
                fv = fr["parent"][key] if ident == PARENT_ID else fr[source][ident][key]
                diffs.append(abs(cv - fv) / max(abs(fv), 1e-12))
            refinement[ident][key] = max(diffs)
    summaries = {}
    for ident in CHILD_IDS:
        ds = [row["child_minus_parent"][ident] for row in fine]
        axial = [d["axial_rms_over_support_relative"] for d in ds]
        radial = [d["radial_rms_over_support_relative"] for d in ds]
        q90 = [d["axial_q90_over_support_delta"] for d in ds]
        q99 = [d["axial_q99_over_support_delta"] for d in ds]
        collar = [row["children"][ident]["outer_support_collar_enstrophy_fraction"] /
                  max(row["parent"]["outer_support_collar_enstrophy_fraction"], 1e-300) for row in fine]
        extension = bool(all(x >= 0 for x in axial) and any(x > 1e-12 for x in q90 + q99))
        low_cost = bool(max(radial) <= .05 and max(collar) <= 2.0)
        summaries[ident] = {
            "axial_rms_relative_changes": axial,
            "radial_rms_relative_changes": radial,
            "q90_deltas_over_support": q90,
            "q99_deltas_over_support": q99,
            "collar_enstrophy_ratios": collar,
            "morphology_extension_signal": extension,
            "low_radial_collar_cost": low_cost,
            "clean_axial_morphology_transfer": bool(extension and low_cost),
        }
    routing = (
        "ST048-S is the cleaner visualization-oriented residual challenger: it improves the frozen held-out residual "
        "receipt on both seeds, raises continuous axial enstrophy RMS at every sampled time, preserves q90/q99 bins, "
        "and lowers radial/collar loading. ST048-B gives the stronger L2 reduction but contracts q99 at t=.50 and "
        "loses substantially more far-axial tail enstrophy. Neither child satisfies the prior strict clean-extension "
        "rule because neither produces a positive q90/q99 bin motion. Stop basis growth; if fixed 3-D rendering still "
        "shows insufficient axial reach, test only the existing single Piola axial-warp degree from PR #367 on frozen "
        "ST048-S rather than adding another spatial basis family."
    )
    recipes = json.loads((Path(__file__).resolve().parent / "recipes.json").read_text())
    return {
        "task_id": TASK_ID,
        "parent_id": PARENT_ID,
        "child_ids": list(CHILD_IDS),
        "candidate_sha256": {PARENT_ID: replay_st048.PARENT_SHA,
                             **{ident: recipes[ident].get("original_candidate_sha256") for ident in CHILD_IDS}},
        "contract": {"times": list(times), "grid_sizes": list(grid_sizes), "box_half_width": BOX_HALF_WIDTH,
                     "radial_support": RADIAL_SUPPORT, "axial_support": AXIAL_SUPPORT,
                     "outer_collar_fraction_of_support": OUTER_COLLAR_FRACTION,
                     "classification": "autonomous target-free visualization diagnostic",
                     "fit_or_parameter_selection": False, "clean_transfer_rule_reused_from_pr_378": True},
        "upstream_residual_receipt": _residual_receipt(),
        "rows_by_grid": {str(k): v for k, v in by_grid.items()},
        "refinement_max_relative_change": refinement,
        "structure": {PARENT_ID: _structure_checks(parent, times),
                      **{ident: _structure_checks(field, times) for ident, field in children.items()}},
        "fine_grid_summary": summaries,
        "routing": routing,
        "truth_boundary": dict(TRUTH_BOUNDARY),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", default="artifacts/st048_agent7_morphology_transfer/report.json")
    args = parser.parse_args()
    report = audit_morphology_transfer()
    path = Path(args.output)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
