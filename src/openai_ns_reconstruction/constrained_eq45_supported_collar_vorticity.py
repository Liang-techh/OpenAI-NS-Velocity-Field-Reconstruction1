"""Independent collar-vorticity audit for the support-connected Eq45 candidate.

This diagnostic is intentionally downstream of the public velocity interfaces.
It samples only ``candidate.at_points(points, time) -> [u,v,w]`` and reconstructs
Cartesian curl with centered finite differences.  It does not reuse profile
jets, streamfunction derivatives, taper derivatives, training tensors, or PDE
operators.

The purpose is narrow: after the physical-space support transform changes the
exterior velocity, check whether the visualization-facing vorticity field
acquires a resolution-dependent or disproportionately large collar artifact.
The result is morphology evidence only, not Navier--Stokes validation and not
an identification of OpenAI's hidden field.
"""
from __future__ import annotations

from collections.abc import Iterable

import numpy as np

from .constrained_eq45_candidate import Eq45VelocityCandidate
from .constrained_eq45_supported_candidate import Eq45SupportedVelocityCandidate


DEFAULT_SPATIAL_STEPS = (0.02, 0.01, 0.005)
DEFAULT_TIMES = (0.25, 0.50, 0.75)

PARENT_SHA256 = "48f1845fcfd71ec95495c98bb7aac3fca4653e748a8856a5421234e9601525a7"
SUPPORTED_SHA256 = "2fdff812c22131d56eac1b7d6e455187ff3207500c6385aa9abf282a8e3d7b1d"

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


def _polar_point(radius: float, azimuth: float, z: float) -> tuple[float, float, float]:
    return (
        float(radius * np.cos(azimuth)),
        float(radius * np.sin(azimuth)),
        float(z),
    )


def default_region_probes() -> dict[str, np.ndarray]:
    """Return deterministic probes safely separated from taper break surfaces."""
    quarter = 0.25 * np.pi
    probes = {
        "plateau": np.asarray(
            [
                _polar_point(0.45, 0.0, 0.00),
                _polar_point(0.80, quarter, 0.25),
                _polar_point(1.10, 0.0, -0.35),
                _polar_point(1.25, quarter, 0.55),
            ],
            dtype=float,
        ),
        "radial_collar": np.asarray(
            [
                _polar_point(1.70, 0.0, 0.00),
                _polar_point(1.80, quarter, 0.25),
                _polar_point(1.88, 0.0, -0.40),
                _polar_point(1.92, quarter, 0.55),
            ],
            dtype=float,
        ),
        "axial_collar": np.asarray(
            [
                _polar_point(0.45, 0.0, 1.70),
                _polar_point(0.80, quarter, -1.80),
                _polar_point(1.10, 0.0, 1.88),
                _polar_point(1.25, quarter, -1.92),
            ],
            dtype=float,
        ),
        "corner_collar": np.asarray(
            [
                _polar_point(1.70, 0.0, 1.70),
                _polar_point(1.78, quarter, -1.78),
                _polar_point(1.86, 0.0, 1.86),
                _polar_point(1.90, quarter, -1.90),
            ],
            dtype=float,
        ),
    }
    for name, values in probes.items():
        if values.shape != (4, 3) or not np.all(np.isfinite(values)):
            raise RuntimeError(f"invalid default probe set {name!r}")
    return probes


def centered_point_vorticity(velocity_fn, points, time: float, step: float) -> np.ndarray:
    """Reconstruct curl at arbitrary points using independent centered stencils."""
    if not callable(velocity_fn):
        raise TypeError("velocity_fn must be callable")
    points_arr = np.asarray(points, dtype=float)
    if points_arr.ndim != 2 or points_arr.shape[1] != 3:
        raise ValueError("points must have shape (n,3)")
    if not np.all(np.isfinite(points_arr)):
        raise ValueError("points must be finite")
    if not np.isfinite(time):
        raise ValueError("time must be finite")
    if not np.isfinite(step) or step <= 0.0:
        raise ValueError("step must be positive and finite")

    derivatives = np.empty((points_arr.shape[0], 3, 3), dtype=float)
    for axis in range(3):
        offset = np.zeros(3, dtype=float)
        offset[axis] = float(step)
        plus = np.asarray(velocity_fn(points_arr + offset, float(time)), dtype=float)
        minus = np.asarray(velocity_fn(points_arr - offset, float(time)), dtype=float)
        if plus.shape != points_arr.shape or minus.shape != points_arr.shape:
            raise ValueError("velocity_fn must return shape (n,3)")
        if not np.all(np.isfinite(plus)) or not np.all(np.isfinite(minus)):
            raise ValueError("velocity_fn returned nonfinite values")
        derivatives[:, :, axis] = (plus - minus) / (2.0 * float(step))

    omega = np.empty_like(points_arr)
    omega[:, 0] = derivatives[:, 2, 1] - derivatives[:, 1, 2]
    omega[:, 1] = derivatives[:, 0, 2] - derivatives[:, 2, 0]
    omega[:, 2] = derivatives[:, 1, 0] - derivatives[:, 0, 1]
    return omega


def _rms(values: np.ndarray) -> float:
    arr = np.asarray(values, dtype=float)
    return float(np.sqrt(np.mean(arr * arr)))


def _vector_rms(vectors: np.ndarray) -> float:
    arr = np.asarray(vectors, dtype=float)
    if arr.ndim != 2 or arr.shape[1] != 3:
        raise ValueError("vectors must have shape (n,3)")
    return float(np.sqrt(np.mean(np.sum(arr * arr, axis=1))))


def _relative_change(value: float, reference: float, scale_floor: float = 1e-14) -> float:
    denominator = max(abs(float(reference)), float(scale_floor))
    return abs(float(value) - float(reference)) / denominator


def audit_supported_collar_vorticity(
    child: Eq45SupportedVelocityCandidate,
    *,
    steps: Iterable[float] = DEFAULT_SPATIAL_STEPS,
    times: Iterable[float] = DEFAULT_TIMES,
    regions: dict[str, np.ndarray] | None = None,
) -> dict[str, object]:
    """Compare parent/child vorticity across plateau and support-collar regions."""
    if not isinstance(child, Eq45SupportedVelocityCandidate):
        raise TypeError("child must be Eq45SupportedVelocityCandidate")
    if child.parent_sha256 != PARENT_SHA256 or child.sha256 != SUPPORTED_SHA256:
        raise ValueError("audit is pinned to the governed frozen parent/child identities")

    step_values = tuple(float(value) for value in steps)
    if len(step_values) < 3 or not all(np.isfinite(step_values)):
        raise ValueError("steps must contain at least three finite levels")
    if not all(value > 0.0 for value in step_values):
        raise ValueError("steps must be positive")
    if not all(b < a for a, b in zip(step_values, step_values[1:])):
        raise ValueError("steps must be strictly decreasing")

    time_values = tuple(float(value) for value in times)
    if not time_values or not all(np.isfinite(time_values)):
        raise ValueError("times must be nonempty and finite")
    if min(time_values) < child.time_start or max(time_values) > child.time_end:
        raise ValueError("times must lie inside the candidate delivery interval")

    region_points = default_region_probes() if regions is None else regions
    required_regions = {"plateau", "radial_collar", "axial_collar", "corner_collar"}
    if set(region_points) != required_regions:
        raise ValueError(f"regions must contain exactly {sorted(required_regions)}")

    rows: list[dict[str, object]] = []
    for time in time_values:
        for step in step_values:
            for region in ("plateau", "radial_collar", "axial_collar", "corner_collar"):
                points = np.asarray(region_points[region], dtype=float)
                parent_omega = centered_point_vorticity(
                    child.parent.at_points, points, time, step
                )
                child_omega = centered_point_vorticity(child.at_points, points, time, step)
                delta = child_omega - parent_omega
                rows.append(
                    {
                        "time": float(time),
                        "step": float(step),
                        "region": region,
                        "parent_vorticity_rms": _vector_rms(parent_omega),
                        "child_vorticity_rms": _vector_rms(child_omega),
                        "taper_induced_vorticity_rms": _vector_rms(delta),
                        "parent_vorticity_max": float(
                            np.max(np.linalg.norm(parent_omega, axis=1))
                        ),
                        "child_vorticity_max": float(
                            np.max(np.linalg.norm(child_omega, axis=1))
                        ),
                    }
                )

    summaries: list[dict[str, object]] = []
    for time in time_values:
        finest_rows = {
            row["region"]: row
            for row in rows
            if row["time"] == time and row["step"] == step_values[-1]
        }
        middle_rows = {
            row["region"]: row
            for row in rows
            if row["time"] == time and row["step"] == step_values[-2]
        }
        plateau_scale = max(
            float(finest_rows["plateau"]["child_vorticity_rms"]), 1e-14
        )
        region_summary: dict[str, object] = {
            "time": float(time),
            "plateau_parent_child_relative_mismatch": _relative_change(
                float(finest_rows["plateau"]["child_vorticity_rms"]),
                float(finest_rows["plateau"]["parent_vorticity_rms"]),
            ),
        }
        for region in ("radial_collar", "axial_collar", "corner_collar"):
            finest = finest_rows[region]
            middle = middle_rows[region]
            region_summary[f"{region}_child_to_plateau_rms_ratio"] = (
                float(finest["child_vorticity_rms"]) / plateau_scale
            )
            region_summary[f"{region}_taper_delta_to_plateau_rms_ratio"] = (
                float(finest["taper_induced_vorticity_rms"]) / plateau_scale
            )
            region_summary[f"{region}_mid_to_fine_relative_delta"] = _relative_change(
                float(middle["child_vorticity_rms"]),
                float(finest["child_vorticity_rms"]),
                scale_floor=1e-12 * plateau_scale,
            )
        summaries.append(region_summary)

    return {
        "schema": "eq45_supported_collar_vorticity_audit_v1",
        "claim_scope": "post_support_transform_visualization_morphology_only",
        "parent_sha256": child.parent_sha256,
        "supported_child_sha256": child.sha256,
        "steps": list(step_values),
        "times": list(time_values),
        "regions": {name: values.tolist() for name, values in region_points.items()},
        "curl_method": "second_order_centered_cartesian_point_stencil",
        "rows": rows,
        "time_summaries": summaries,
        "truth_boundary": dict(_TRUTH_BOUNDARY),
    }


def governed_supported_seed() -> Eq45SupportedVelocityCandidate:
    """Return the currently governed support-connected frozen Eq45 seed."""
    parent = Eq45VelocityCandidate.seed()
    if parent.sha256 != PARENT_SHA256:
        raise RuntimeError("frozen Eq45 seed identity drifted")
    child = Eq45SupportedVelocityCandidate(parent=parent)
    if child.sha256 != SUPPORTED_SHA256:
        raise RuntimeError("support-connected Eq45 child identity drifted")
    return child
