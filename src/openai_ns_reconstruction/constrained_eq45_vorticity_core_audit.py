"""Resolution audit for Eq45 vorticity concentration and core geometry.

The production path deliberately consumes only the public
``Eq45VelocityCandidate.at_points(points, time) -> [u,v,w]`` interface.  Curl is
then reconstructed with independent centered finite differences on a Cartesian
cell-centred grid.  No profile derivatives or optimizer internals are reused.

All grid, box, isovalue and summary choices here are autonomous visualization
 diagnostics.  The report is a finite-window morphology fingerprint only: it is
not Navier--Stokes validation, a blow-up claim, or identification of OpenAI's
hidden field.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Callable, Iterable

import numpy as np

from .constrained_eq45_candidate import Eq45VelocityCandidate


DEFAULT_BOX_HALF_WIDTH = 2.0
DEFAULT_GRID_SIZES = (48, 64, 80)
DEFAULT_TIMES = (0.25, 0.50, 0.75)
DEFAULT_HALF_MAX_FRACTION = 0.5

_TRUTH_BOUNDARY = {
    "velocity_changed": False,
    "physical_support_validated": False,
    "visualization_ready": False,
    "visual_correspondence_verified": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


def _cell_centres(grid_size: int, box_half_width: float) -> tuple[np.ndarray, float]:
    if not isinstance(grid_size, (int, np.integer)) or grid_size < 5:
        raise ValueError("grid_size must be an integer >= 5")
    if not np.isfinite(box_half_width) or box_half_width <= 0.0:
        raise ValueError("box_half_width must be positive and finite")
    dx = 2.0 * float(box_half_width) / int(grid_size)
    axis = np.linspace(
        -box_half_width + 0.5 * dx,
        box_half_width - 0.5 * dx,
        int(grid_size),
    )
    return axis, dx


def centered_vorticity(velocity: np.ndarray, spacing: float) -> np.ndarray:
    """Return centered-difference curl on the strict interior of a cubic grid."""
    values = np.asarray(velocity, dtype=float)
    if values.ndim != 4 or values.shape[-1] != 3:
        raise ValueError("velocity must have shape (n,n,n,3)")
    if len(set(values.shape[:3])) != 1 or values.shape[0] < 5:
        raise ValueError("velocity grid must be cubic with n >= 5")
    if not np.all(np.isfinite(values)):
        raise ValueError("velocity must be finite")
    if not np.isfinite(spacing) or spacing <= 0.0:
        raise ValueError("spacing must be positive and finite")

    u = values[..., 0]
    v = values[..., 1]
    w = values[..., 2]
    scale = 1.0 / (2.0 * float(spacing))

    du_dy = (u[:, 2:, :] - u[:, :-2, :]) * scale
    du_dz = (u[:, :, 2:] - u[:, :, :-2]) * scale
    dv_dx = (v[2:, :, :] - v[:-2, :, :]) * scale
    dv_dz = (v[:, :, 2:] - v[:, :, :-2]) * scale
    dw_dx = (w[2:, :, :] - w[:-2, :, :]) * scale
    dw_dy = (w[:, 2:, :] - w[:, :-2, :]) * scale

    omega_x = dw_dy[1:-1, :, 1:-1] - dv_dz[1:-1, 1:-1, :]
    omega_y = du_dz[1:-1, 1:-1, :] - dw_dx[:, 1:-1, 1:-1]
    omega_z = dv_dx[:, 1:-1, 1:-1] - du_dy[1:-1, :, 1:-1]
    out = np.stack((omega_x, omega_y, omega_z), axis=-1)
    if not np.all(np.isfinite(out)):
        raise RuntimeError("vorticity reconstruction became nonfinite")
    return out


def _weighted_quantile(values: np.ndarray, weights: np.ndarray, fraction: float) -> float:
    if not 0.0 <= fraction <= 1.0:
        raise ValueError("fraction must lie in [0,1]")
    values = np.asarray(values, dtype=float).ravel()
    weights = np.asarray(weights, dtype=float).ravel()
    if values.shape != weights.shape:
        raise ValueError("values and weights must have the same size")
    mask = np.isfinite(values) & np.isfinite(weights) & (weights > 0.0)
    if not np.any(mask):
        raise ValueError("weighted quantile requires positive finite weight")
    values = values[mask]
    weights = weights[mask]
    order = np.argsort(values)
    values = values[order]
    cumulative = np.cumsum(weights[order])
    index = np.searchsorted(cumulative, fraction * cumulative[-1], side="left")
    return float(values[min(index, values.size - 1)])


def vorticity_core_row(
    velocity_fn: Callable[[np.ndarray, float], np.ndarray],
    *,
    grid_size: int,
    time: float,
    box_half_width: float = DEFAULT_BOX_HALF_WIDTH,
    half_max_fraction: float = DEFAULT_HALF_MAX_FRACTION,
) -> dict[str, float | int]:
    """Measure one finite-window vorticity/concentration fingerprint."""
    if not callable(velocity_fn):
        raise TypeError("velocity_fn must be callable")
    if not np.isfinite(time):
        raise ValueError("time must be finite")
    if not 0.0 < half_max_fraction < 1.0:
        raise ValueError("half_max_fraction must lie strictly in (0,1)")

    axis, dx = _cell_centres(grid_size, box_half_width)
    x, y, z = np.meshgrid(axis, axis, axis, indexing="ij")
    points = np.stack((x, y, z), axis=-1)
    velocity = np.asarray(velocity_fn(points, float(time)), dtype=float)
    if velocity.shape != points.shape or not np.all(np.isfinite(velocity)):
        raise ValueError("velocity_fn must return finite values with shape points.shape")

    omega = centered_vorticity(velocity, dx)
    magnitude = np.linalg.norm(omega, axis=-1)
    enstrophy_density = magnitude * magnitude
    total_weight = float(np.sum(enstrophy_density))
    if total_weight <= 0.0:
        raise ValueError("vorticity audit requires nonzero enstrophy")

    xi = x[1:-1, 1:-1, 1:-1]
    yi = y[1:-1, 1:-1, 1:-1]
    zi = z[1:-1, 1:-1, 1:-1]
    radius = np.sqrt(xi * xi + yi * yi)
    abs_z = np.abs(zi)

    radial_rms = float(np.sqrt(np.sum(radius * radius * enstrophy_density) / total_weight))
    axial_rms = float(np.sqrt(np.sum(zi * zi * enstrophy_density) / total_weight))
    omega_max = float(np.max(magnitude))
    high = magnitude >= half_max_fraction * omega_max

    return {
        "grid_size": int(grid_size),
        "time": float(time),
        "box_half_width": float(box_half_width),
        "spacing": float(dx),
        "omega_rms": float(np.sqrt(np.mean(enstrophy_density))),
        "omega_max": omega_max,
        "enstrophy_integral": float(np.sum(enstrophy_density) * dx**3),
        "enstrophy_radial_q50": _weighted_quantile(radius, enstrophy_density, 0.50),
        "enstrophy_radial_q90": _weighted_quantile(radius, enstrophy_density, 0.90),
        "enstrophy_abs_z_q50": _weighted_quantile(abs_z, enstrophy_density, 0.50),
        "enstrophy_abs_z_q90": _weighted_quantile(abs_z, enstrophy_density, 0.90),
        "enstrophy_radial_rms": radial_rms,
        "enstrophy_axial_rms": axial_rms,
        "enstrophy_rms_aspect_ratio": axial_rms / radial_rms,
        "halfmax_radial_extent": float(np.max(radius[high])),
        "halfmax_abs_z_extent": float(np.max(abs_z[high])),
        "halfmax_volume_fraction": float(np.mean(high)),
    }


def _relative_delta(value: float, reference: float) -> float:
    denominator = max(abs(float(reference)), np.finfo(float).tiny)
    return abs(float(value) - float(reference)) / denominator


def audit_eq45_vorticity_core(
    candidate: Eq45VelocityCandidate,
    *,
    grid_sizes: Iterable[int] = DEFAULT_GRID_SIZES,
    times: Iterable[float] = DEFAULT_TIMES,
    box_half_width: float = DEFAULT_BOX_HALF_WIDTH,
) -> dict[str, object]:
    """Run a three-or-more-resolution audit through the public velocity API."""
    if not isinstance(candidate, Eq45VelocityCandidate):
        raise TypeError("candidate must be Eq45VelocityCandidate")
    grids = tuple(int(value) for value in grid_sizes)
    if len(grids) < 3 or any(b <= a for a, b in zip(grids, grids[1:])):
        raise ValueError("grid_sizes must contain at least three strictly increasing levels")
    audit_times = tuple(float(value) for value in times)
    if not audit_times or not all(np.isfinite(audit_times)):
        raise ValueError("times must be nonempty and finite")
    if min(audit_times) < candidate.time_start or max(audit_times) > candidate.time_end:
        raise ValueError("times must lie inside the candidate delivery interval")

    rows = [
        vorticity_core_row(
            candidate.at_points,
            grid_size=grid_size,
            time=time,
            box_half_width=box_half_width,
        )
        for grid_size in grids
        for time in audit_times
    ]

    summaries = []
    for time in audit_times:
        time_rows = [row for row in rows if row["time"] == time]
        middle = time_rows[-2]
        finest = time_rows[-1]
        summaries.append(
            {
                "time": float(time),
                "finest_grid_size": int(finest["grid_size"]),
                "finest_omega_max": float(finest["omega_max"]),
                "finest_enstrophy_integral": float(finest["enstrophy_integral"]),
                "finest_enstrophy_radial_rms": float(finest["enstrophy_radial_rms"]),
                "finest_enstrophy_axial_rms": float(finest["enstrophy_axial_rms"]),
                "finest_enstrophy_rms_aspect_ratio": float(
                    finest["enstrophy_rms_aspect_ratio"]
                ),
                "finest_enstrophy_radial_q90": float(finest["enstrophy_radial_q90"]),
                "finest_enstrophy_abs_z_q90": float(finest["enstrophy_abs_z_q90"]),
                "mid_to_finest_omega_max_relative_delta": _relative_delta(
                    middle["omega_max"], finest["omega_max"]
                ),
                "mid_to_finest_enstrophy_relative_delta": _relative_delta(
                    middle["enstrophy_integral"], finest["enstrophy_integral"]
                ),
                "mid_to_finest_radial_rms_relative_delta": _relative_delta(
                    middle["enstrophy_radial_rms"], finest["enstrophy_radial_rms"]
                ),
                "mid_to_finest_axial_rms_relative_delta": _relative_delta(
                    middle["enstrophy_axial_rms"], finest["enstrophy_axial_rms"]
                ),
                "mid_to_finest_aspect_relative_delta": _relative_delta(
                    middle["enstrophy_rms_aspect_ratio"],
                    finest["enstrophy_rms_aspect_ratio"],
                ),
                "mid_to_finest_radial_q90_relative_delta": _relative_delta(
                    middle["enstrophy_radial_q90"], finest["enstrophy_radial_q90"]
                ),
                "mid_to_finest_abs_z_q90_relative_delta": _relative_delta(
                    middle["enstrophy_abs_z_q90"], finest["enstrophy_abs_z_q90"]
                ),
            }
        )

    return {
        "schema": "eq45_public_vorticity_core_audit_v1",
        "claim_scope": "finite_window_visualization_vorticity_fingerprint_only",
        "candidate_sha256": candidate.sha256,
        "box_half_width": float(box_half_width),
        "grid_sizes": list(grids),
        "times": list(audit_times),
        "curl_method": "second_order_centered_on_cell_centred_cartesian_grid",
        "halfmax_fraction": DEFAULT_HALF_MAX_FRACTION,
        "rows": rows,
        "time_summaries": summaries,
        "truth_boundary": dict(_TRUTH_BOUNDARY),
    }


def run_seed_audit(candidate_path: str | Path) -> dict[str, object]:
    candidate = Eq45VelocityCandidate.load_json(candidate_path)
    report = audit_eq45_vorticity_core(candidate)
    report["candidate_path"] = str(candidate_path)
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("candidate", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = run_seed_audit(args.candidate)
    text = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output is None:
        print(text, end="")
    else:
        args.output.write_text(text, encoding="utf-8")


if __name__ == "__main__":
    main()
