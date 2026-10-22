"""Truth-bounded 3-D vorticity morphology audit for retained ST006.

This module samples the published ``velocity(x,y,z,t)->[u,v,w]`` field on
fixed Cartesian grids and describes the spatial distribution of vorticity.
It is a visualization diagnostic only: it does not fit the candidate, a
camera, a hidden frame time, pressure, forcing, or any acceptance threshold.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np

from . import CANDIDATE_SHA256, load_best

SCHEMA = "st006_vorticity_morphology_timeseries_v1"
REFERENCE_TIMES = (0.25, 0.50, 0.75)
RESOLUTION_LADDER = (17, 25, 33)
BOX = (-2.0, 2.0)
PROVENANCE = "retained_ST006_public_velocity_api_fullbox_vorticity_morphology"

_FALSE_TRUTH_FLAGS = {
    "visualization_ready": False,
    "visual_correspondence_verified": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
    "hidden_frame_time_inferred": False,
    "camera_fitted": False,
    "visual_acceptance_threshold_selected": False,
    "used_for_pde_acceptance": False,
}


def _digest(payload: dict[str, Any]) -> str:
    raw = json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
    return hashlib.sha256(raw).hexdigest()


def _valid_sha256(value: str) -> bool:
    if not isinstance(value, str) or len(value) != 64:
        return False
    return all(c in "0123456789abcdef" for c in value)


def _canonical_unsigned_axis(vector: np.ndarray) -> np.ndarray:
    axis = np.asarray(vector, dtype=float).copy()
    index = int(np.argmax(np.abs(axis)))
    if axis[index] < 0.0:
        axis *= -1.0
    return axis


def _diagnose_one_grid(
    velocity: Any,
    *,
    time: float,
    resolution: int,
    candidate_sha256: str,
    provenance: str,
) -> dict[str, Any]:
    """Measure one fixed Cartesian grid without making a visual pass/fail claim."""
    if not _valid_sha256(candidate_sha256):
        raise ValueError("candidate_sha256 must be 64 lowercase hexadecimal characters")
    if not isinstance(provenance, str) or not provenance.strip():
        raise ValueError("provenance must be a nonempty string")
    if not isinstance(resolution, (int, np.integer)) or resolution < 9 or resolution > 65:
        raise ValueError("resolution must be an integer in [9,65]")
    if resolution % 2 == 0:
        raise ValueError("resolution must be odd")
    time = float(time)
    if not np.isfinite(time) or not 0.25 <= time <= 0.75:
        raise ValueError("time must lie in the registered [0.25,0.75] interval")

    axis = np.linspace(BOX[0], BOX[1], int(resolution), dtype=float)
    spacing = float(axis[1] - axis[0])
    x, y, z = np.meshgrid(axis, axis, axis, indexing="ij")
    points = np.column_stack((x.ravel(), y.ravel(), z.ravel()))
    uvw = np.asarray(velocity(points, time), dtype=float)
    if uvw.shape != points.shape:
        raise ValueError("velocity output must have shape (n,3)")
    if not np.isfinite(uvw).all():
        raise ValueError("velocity output must be finite")
    velocity_grid = uvw.reshape(resolution, resolution, resolution, 3)
    u = velocity_grid[..., 0]
    v = velocity_grid[..., 1]
    w = velocity_grid[..., 2]

    inv_2h = 0.5 / spacing
    ux = (u[2:, 1:-1, 1:-1] - u[:-2, 1:-1, 1:-1]) * inv_2h
    uy = (u[1:-1, 2:, 1:-1] - u[1:-1, :-2, 1:-1]) * inv_2h
    uz = (u[1:-1, 1:-1, 2:] - u[1:-1, 1:-1, :-2]) * inv_2h
    vx = (v[2:, 1:-1, 1:-1] - v[:-2, 1:-1, 1:-1]) * inv_2h
    vy = (v[1:-1, 2:, 1:-1] - v[1:-1, :-2, 1:-1]) * inv_2h
    vz = (v[1:-1, 1:-1, 2:] - v[1:-1, 1:-1, :-2]) * inv_2h
    wx = (w[2:, 1:-1, 1:-1] - w[:-2, 1:-1, 1:-1]) * inv_2h
    wy = (w[1:-1, 2:, 1:-1] - w[1:-1, :-2, 1:-1]) * inv_2h
    wz = (w[1:-1, 1:-1, 2:] - w[1:-1, 1:-1, :-2]) * inv_2h

    omega = np.stack((wy - vz, uz - wx, vx - uy), axis=-1)
    omega_mag = np.linalg.norm(omega, axis=-1)
    weights = omega_mag * omega_mag
    total_weight = float(np.sum(weights))
    velocity_rms = float(np.sqrt(np.mean(uvw * uvw)))
    omega_rms = float(np.sqrt(np.mean(weights)))
    omega_max = float(np.max(omega_mag))
    if not np.isfinite(velocity_rms) or velocity_rms <= 0.0:
        raise ValueError("sampled velocity must be nonzero and finite")
    weight_floor = float(
        256.0
        * np.finfo(float).eps
        * max(1.0, omega_max * omega_max)
        * weights.size
    )
    if not np.isfinite(total_weight) or total_weight <= weight_floor:
        raise ValueError("sampled vorticity is zero/nonfinite at floating-point scale")

    interior_axis = axis[1:-1]
    xi, yi, zi = np.meshgrid(interior_axis, interior_axis, interior_axis, indexing="ij")
    positions = np.stack((xi, yi, zi), axis=-1).reshape(-1, 3)
    flat_weights = weights.reshape(-1)
    centroid = np.sum(flat_weights[:, None] * positions, axis=0) / total_weight
    centered = positions - centroid
    covariance = np.einsum("n,ni,nj->ij", flat_weights, centered, centered) / total_weight
    covariance = 0.5 * (covariance + covariance.T)

    eigenvalues, eigenvectors = np.linalg.eigh(covariance)
    scale = max(1.0, float(np.max(np.abs(eigenvalues))))
    negative_tolerance = 256.0 * np.finfo(float).eps * scale
    if float(np.min(eigenvalues)) < -negative_tolerance:
        raise ValueError("enstrophy covariance is not positive semidefinite")
    eigenvalues = np.maximum(eigenvalues, 0.0)
    transverse_variance = float(np.mean(eigenvalues[:2]))
    if transverse_variance <= 256.0 * np.finfo(float).eps * scale:
        raise ValueError("transverse enstrophy extent is numerically degenerate")

    principal_axis = _canonical_unsigned_axis(eigenvectors[:, -1])
    principal_rms_extent = float(np.sqrt(eigenvalues[-1]))
    transverse_rms_extent = float(np.sqrt(transverse_variance))
    aspect_ratio = float(principal_rms_extent / transverse_rms_extent)
    covariance_volume_proxy = float(np.sqrt(np.prod(eigenvalues)))
    divergence_rms = float(np.sqrt(np.mean((ux + vy + wz) ** 2)))
    enstrophy_integral = float(np.sum(weights) * spacing**3)

    return {
        "time": time,
        "resolution": int(resolution),
        "grid_spacing": spacing,
        "velocity_rms": velocity_rms,
        "vorticity_rms": omega_rms,
        "vorticity_max": omega_max,
        "integral_omega_squared": enstrophy_integral,
        "same_operator_divergence_rms": divergence_rms,
        "enstrophy_weighted_centroid": [float(x) for x in centroid],
        "enstrophy_covariance": [[float(x) for x in row] for row in covariance],
        "enstrophy_covariance_eigenvalues": [float(x) for x in eigenvalues],
        "principal_axis_registered_coordinates": [float(x) for x in principal_axis],
        "principal_rms_extent": principal_rms_extent,
        "transverse_rms_extent": transverse_rms_extent,
        "principal_to_transverse_aspect_ratio": aspect_ratio,
        "covariance_rms_volume_proxy": covariance_volume_proxy,
    }


def _safe_ratio(final: float, initial: float, name: str) -> float:
    if not np.isfinite(initial) or not np.isfinite(final) or initial <= 0.0:
        raise ValueError(f"cannot form endpoint ratio for {name}")
    return float(final / initial)


def _measure_field(field: Any) -> dict[str, Any]:
    if getattr(field, "sha256", None) != CANDIDATE_SHA256:
        raise ValueError("vorticity-morphology audit is bound to retained ST006 identity")

    by_resolution: list[dict[str, Any]] = []
    finest_by_time: dict[float, dict[str, Any]] = {}
    for resolution in RESOLUTION_LADDER:
        rows = []
        for time in REFERENCE_TIMES:
            row = _diagnose_one_grid(
                field.at_points,
                time=time,
                resolution=resolution,
                candidate_sha256=CANDIDATE_SHA256,
                provenance=PROVENANCE,
            )
            rows.append(row)
            if resolution == RESOLUTION_LADDER[-1]:
                finest_by_time[time] = row
        by_resolution.append({"resolution": resolution, "per_time": rows})

    finest = [finest_by_time[t] for t in REFERENCE_TIMES]
    first = finest[0]
    last = finest[-1]
    metric_names = (
        "velocity_rms",
        "vorticity_rms",
        "vorticity_max",
        "integral_omega_squared",
        "principal_rms_extent",
        "transverse_rms_extent",
        "principal_to_transverse_aspect_ratio",
        "covariance_rms_volume_proxy",
    )
    endpoint_trend = {
        name: {
            "initial": float(first[name]),
            "final": float(last[name]),
            "delta": float(last[name] - first[name]),
            "ratio_final_over_initial": _safe_ratio(float(last[name]), float(first[name]), name),
        }
        for name in metric_names
    }

    resolution_drift: list[dict[str, Any]] = []
    finest_resolution = RESOLUTION_LADDER[-1]
    for t in REFERENCE_TIMES:
        finest_row = finest_by_time[t]
        comparisons = []
        for block in by_resolution[:-1]:
            coarse = next(row for row in block["per_time"] if row["time"] == t)
            comparisons.append(
                {
                    "resolution": int(block["resolution"]),
                    "relative_to_finest": {
                        name: float(
                            abs(coarse[name] - finest_row[name])
                            / max(abs(finest_row[name]), np.finfo(float).tiny)
                        )
                        for name in (
                            "vorticity_rms",
                            "principal_rms_extent",
                            "transverse_rms_extent",
                            "principal_to_transverse_aspect_ratio",
                            "covariance_rms_volume_proxy",
                        )
                    },
                }
            )
        resolution_drift.append(
            {
                "time": float(t),
                "finest_resolution": finest_resolution,
                "coarser_comparisons": comparisons,
            }
        )

    payload: dict[str, Any] = {
        "schema": SCHEMA,
        "candidate": "ST006",
        "candidate_sha256": CANDIDATE_SHA256,
        "provenance": PROVENANCE,
        "measurement_contract": {
            "physical_box": [[BOX[0], BOX[1]]] * 3,
            "times": [float(t) for t in REFERENCE_TIMES],
            "resolutions": [int(n) for n in RESOLUTION_LADDER],
            "derivative": "centered second-order Cartesian finite differences on interior nodes",
            "vorticity_weight": "|curl(u)|^2",
            "geometry": "weighted centroid and symmetric 3x3 covariance eigensystem",
            "principal_axis_semantics": "unsigned axis in the registered physical coordinates; no camera registration",
            "target_values_or_visual_thresholds": None,
        },
        "by_resolution": by_resolution,
        "finest_resolution_endpoint_trend": endpoint_trend,
        "resolution_drift": resolution_drift,
        "interpretation_boundary": {
            "candidate_side_only": True,
            "descriptive_uses": [
                "3D vortex-core elongation/contraction trend review",
                "vorticity concentration and principal-axis visualization planning",
                "same-contract comparison with a later production candidate",
            ],
            "does_not_infer": [
                "OpenAI numerical velocity",
                "OpenAI physical coordinate scale",
                "OpenAI hidden frame times",
                "OpenAI camera pose or projection",
                "OpenAI streamline seed locations",
                "OpenAI hidden parameters",
            ],
            "exported_grid_derivatives_are_not_pde_acceptance_evidence": True,
            **_FALSE_TRUTH_FLAGS,
        },
    }
    payload["report_sha256"] = _digest(payload)
    return payload


def measure_retained_st006_vorticity_morphology() -> dict[str, Any]:
    """Measure the retained ST006 field on the frozen 3-D morphology contract."""
    return _measure_field(load_best())


def write_report(path: str | Path) -> dict[str, Any]:
    path = Path(path)
    if path.exists():
        raise FileExistsError(f"refusing to overwrite existing report: {path}")
    report = measure_retained_st006_vorticity_morphology()
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n")
    return report


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args(argv)
    report = write_report(args.output) if args.output is not None else measure_retained_st006_vorticity_morphology()
    print(json.dumps(report, indent=2, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
