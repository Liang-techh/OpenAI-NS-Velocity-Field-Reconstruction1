"""Whole-domain vorticity-envelope audit for the support-connected Eq45 field.

This diagnostic is downstream of the public velocity interface.  It samples the
parent and support-connected child only through ``at_points(...)->[u,v,w]``,
reconstructs Cartesian curl on structured grids, and compares a fixed
vorticity-superlevel envelope.  It complements the pointwise collar audit by
asking whether taper-generated vorticity can alter a visualization-facing
isosurface envelope.

The result is morphology/resolution evidence only.  It is not Navier--Stokes
validation and it does not identify OpenAI's hidden velocity field.
"""
from __future__ import annotations

from collections.abc import Iterable

import numpy as np

from .constrained_eq45_supported_candidate import Eq45SupportedVelocityCandidate
from .constrained_eq45_supported_collar_vorticity import governed_supported_seed


DEFAULT_RESOLUTIONS = (49, 65, 81)
DEFAULT_TIMES = (0.25, 0.50, 0.75)
DEFAULT_HALF_WIDTH = 2.0
DEFAULT_SUPERLEVEL_FRACTION = 0.25
DEFAULT_PLATEAU_RADIUS = 1.5
DEFAULT_PLATEAU_HALF_HEIGHT = 1.5
DEFAULT_COLLAR_START = 1.6

_TRUTH_BOUNDARY = {
    "velocity_changed": False,
    "physical_support_connection_implemented": True,
    "physical_support_validated": False,
    "visualization_ready": False,
    "visual_correspondence_verified": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


def structured_vorticity_magnitude(candidate, axis: np.ndarray, time: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return ``|curl u|``, cylindrical radius and ``|z|`` on one Cartesian grid."""
    axis = np.asarray(axis, dtype=float)
    if axis.ndim != 1 or axis.size < 5 or not np.all(np.isfinite(axis)):
        raise ValueError("axis must be a finite 1D grid with at least five points")
    spacing = np.diff(axis)
    if not np.all(spacing > 0.0) or not np.allclose(spacing, spacing[0], rtol=1e-12, atol=1e-14):
        raise ValueError("axis must be strictly increasing and uniformly spaced")
    if not callable(getattr(candidate, "at_points", None)):
        raise TypeError("candidate must expose at_points(points,time)")

    x, y, z = np.meshgrid(axis, axis, axis, indexing="ij")
    points = np.stack((x, y, z), axis=-1)
    velocity = np.asarray(candidate.at_points(points, float(time)), dtype=float)
    if velocity.shape != points.shape or not np.all(np.isfinite(velocity)):
        raise ValueError("candidate.at_points must return finite [...,3] velocity")

    dx = float(spacing[0])
    du_dx, du_dy, du_dz = np.gradient(velocity[..., 0], dx, dx, dx, edge_order=2)
    dv_dx, dv_dy, dv_dz = np.gradient(velocity[..., 1], dx, dx, dx, edge_order=2)
    dw_dx, dw_dy, dw_dz = np.gradient(velocity[..., 2], dx, dx, dx, edge_order=2)
    omega_x = dw_dy - dv_dz
    omega_y = du_dz - dw_dx
    omega_z = dv_dx - du_dy
    magnitude = np.sqrt(omega_x * omega_x + omega_y * omega_y + omega_z * omega_z)
    radius = np.hypot(x, y)
    abs_z = np.abs(z)
    if not np.all(np.isfinite(magnitude)):
        raise RuntimeError("structured curl produced nonfinite vorticity")
    return magnitude, radius, abs_z


def _metrics(magnitude, radius, abs_z, threshold: float, collar_start: float) -> dict[str, float | int]:
    active = magnitude >= float(threshold)
    if not np.any(active):
        raise RuntimeError("vorticity superlevel set is empty")
    r = radius[active]
    z = abs_z[active]
    collar = (radius > collar_start) | (abs_z > collar_start)
    total_vorticity2 = float(np.sum(magnitude * magnitude))
    collar_vorticity2 = float(np.sum((magnitude[collar]) ** 2))
    return {
        "active_voxels": int(np.count_nonzero(active)),
        "radial_q99": float(np.quantile(r, 0.99)),
        "axial_q99": float(np.quantile(z, 0.99)),
        "radial_max": float(np.max(r)),
        "axial_max": float(np.max(z)),
        "aspect_q99": float(np.quantile(z, 0.99) / max(np.quantile(r, 0.99), 1e-15)),
        "superlevel_collar_voxel_fraction": float(np.mean(collar[active])),
        "whole_grid_collar_vorticity2_fraction": collar_vorticity2 / max(total_vorticity2, 1e-300),
    }


def audit_supported_vorticity_envelope(
    child: Eq45SupportedVelocityCandidate | None = None,
    *,
    resolutions: Iterable[int] = DEFAULT_RESOLUTIONS,
    times: Iterable[float] = DEFAULT_TIMES,
    half_width: float = DEFAULT_HALF_WIDTH,
    superlevel_fraction: float = DEFAULT_SUPERLEVEL_FRACTION,
    plateau_radius: float = DEFAULT_PLATEAU_RADIUS,
    plateau_half_height: float = DEFAULT_PLATEAU_HALF_HEIGHT,
    collar_start: float = DEFAULT_COLLAR_START,
) -> dict[str, object]:
    """Compare parent/child 25%-core vorticity envelopes over three grid levels."""
    child = governed_supported_seed() if child is None else child
    if not isinstance(child, Eq45SupportedVelocityCandidate):
        raise TypeError("child must be Eq45SupportedVelocityCandidate")
    resolution_values = tuple(int(value) for value in resolutions)
    if len(resolution_values) < 3 or any(value < 5 for value in resolution_values):
        raise ValueError("resolutions must contain at least three levels >= 5")
    if not all(b > a for a, b in zip(resolution_values, resolution_values[1:])):
        raise ValueError("resolutions must be strictly increasing")
    time_values = tuple(float(value) for value in times)
    if not time_values or not all(np.isfinite(time_values)):
        raise ValueError("times must be nonempty and finite")
    if min(time_values) < child.time_start or max(time_values) > child.time_end:
        raise ValueError("times must lie inside the candidate delivery interval")
    for value, name in ((half_width, "half_width"), (plateau_radius, "plateau_radius"), (plateau_half_height, "plateau_half_height"), (collar_start, "collar_start")):
        if not np.isfinite(value) or value <= 0.0:
            raise ValueError(f"{name} must be positive and finite")
    if not np.isfinite(superlevel_fraction) or not (0.0 < superlevel_fraction < 1.0):
        raise ValueError("superlevel_fraction must lie in (0,1)")

    rows: list[dict[str, object]] = []
    for resolution in resolution_values:
        axis = np.linspace(-half_width, half_width, resolution)
        dx = float(axis[1] - axis[0])
        for time in time_values:
            parent_mag, radius, abs_z = structured_vorticity_magnitude(child.parent, axis, time)
            child_mag, child_radius, child_abs_z = structured_vorticity_magnitude(child, axis, time)
            plateau = (radius <= plateau_radius) & (abs_z <= plateau_half_height)
            core_peak = float(np.max(parent_mag[plateau]))
            threshold = float(superlevel_fraction * core_peak)
            rows.append({
                "resolution": resolution,
                "dx": dx,
                "time": time,
                "core_peak": core_peak,
                "threshold": threshold,
                "parent": _metrics(parent_mag, radius, abs_z, threshold, collar_start),
                "child": _metrics(child_mag, child_radius, child_abs_z, threshold, collar_start),
            })

    summaries: list[dict[str, float]] = []
    for time in time_values:
        fine = next(row for row in rows if row["resolution"] == resolution_values[-1] and row["time"] == time)
        middle = next(row for row in rows if row["resolution"] == resolution_values[-2] and row["time"] == time)
        p = fine["parent"]
        c = fine["child"]
        c_mid = middle["child"]
        summaries.append({
            "time": time,
            "child_parent_radial_q99_ratio": float(c["radial_q99"] / max(p["radial_q99"], 1e-15)),
            "child_parent_axial_q99_ratio": float(c["axial_q99"] / max(p["axial_q99"], 1e-15)),
            "child_parent_aspect_ratio": float(c["aspect_q99"] / max(p["aspect_q99"], 1e-15)),
            "child_radial_q99_mid_to_fine_relative_change": float(abs(c["radial_q99"] - c_mid["radial_q99"]) / max(abs(c["radial_q99"]), 1e-15)),
            "parent_whole_grid_collar_vorticity2_fraction": float(p["whole_grid_collar_vorticity2_fraction"]),
            "child_whole_grid_collar_vorticity2_fraction": float(c["whole_grid_collar_vorticity2_fraction"]),
            "child_minus_parent_collar_vorticity2_fraction": float(c["whole_grid_collar_vorticity2_fraction"] - p["whole_grid_collar_vorticity2_fraction"]),
        })

    return {
        "schema": "eq45_supported_vorticity_envelope_audit_v1",
        "claim_scope": "whole_domain_visualization_morphology_only",
        "parent_sha256": child.parent_sha256,
        "supported_child_sha256": child.sha256,
        "resolutions": list(resolution_values),
        "times": list(time_values),
        "half_width": float(half_width),
        "superlevel_fraction": float(superlevel_fraction),
        "threshold_reference": "fraction_of_parent_vorticity_peak_inside_r<=1.5_absz<=1.5",
        "curl_method": "numpy_second_order_structured_cartesian_gradient",
        "rows": rows,
        "time_summaries": summaries,
        "truth_boundary": dict(_TRUTH_BOUNDARY),
    }
