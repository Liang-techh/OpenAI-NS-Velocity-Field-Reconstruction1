"""Energy-concentration fingerprint for the bounded optimized Eq45 candidate.

This module is a production-side visualization/research diagnostic.  It compares
exactly the frozen canonical Eq45 candidate and the bounded profile update already
recorded by CR005/CR007.  Every velocity sample is obtained through
``Eq45VelocityCandidate.at_points(...)->[u,v,w]``.

The quantities here are finite-window kinetic-energy concentration metrics.  They
are not the CR001 physical-energy normalization, do not validate physical compact
support, and cannot establish visual correspondence, Navier--Stokes validity,
paper exactness, or recovery of an OpenAI hidden field.
"""
from __future__ import annotations

from typing import Iterable

import numpy as np

from .constrained_eq45_candidate import Eq45VelocityCandidate
from .constrained_eq45_optimized_energy import (
    DEFAULT_HALF_WIDTH,
    DEFAULT_RESOLUTIONS,
    DEFAULT_TIMES,
    _candidate_pair,
    _midpoint_axis,
    _validated_resolutions,
    _validated_times,
)


def _weighted_quantile(values: np.ndarray, weights: np.ndarray, quantile: float) -> float:
    values = np.asarray(values, dtype=float).reshape(-1)
    weights = np.asarray(weights, dtype=float).reshape(-1)
    quantile = float(quantile)
    if values.shape != weights.shape or values.size == 0:
        raise ValueError("values and weights must be nonempty arrays with equal shape")
    if not np.all(np.isfinite(values)) or not np.all(np.isfinite(weights)):
        raise ValueError("values and weights must be finite")
    if np.any(weights < 0.0) or not (0.0 < quantile < 1.0):
        raise ValueError("weights must be nonnegative and quantile must lie in (0,1)")
    total = float(np.sum(weights))
    if not np.isfinite(total) or total <= 0.0:
        raise ValueError("weights must have positive finite total")
    order = np.argsort(values, kind="mergesort")
    sorted_values = values[order]
    cumulative = np.cumsum(weights[order])
    index = int(np.searchsorted(cumulative, quantile * total, side="left"))
    return float(sorted_values[min(index, sorted_values.size - 1)])


def finite_window_energy_concentration(
    candidate: Eq45VelocityCandidate,
    time: float,
    *,
    half_width: float = DEFAULT_HALF_WIDTH,
    resolution: int = 32,
) -> dict[str, float]:
    """Return energy-weighted radial/axial extents on one Cartesian midpoint grid."""
    if not isinstance(candidate, Eq45VelocityCandidate):
        raise TypeError("candidate must be Eq45VelocityCandidate")
    time = float(time)
    half_width = float(half_width)
    if not np.isfinite(time):
        raise ValueError("time must be finite")
    if not np.isfinite(half_width) or half_width <= 0.0:
        raise ValueError("half_width must be positive and finite")
    if not isinstance(resolution, (int, np.integer)) or int(resolution) < 4:
        raise ValueError("resolution must be an integer >=4")
    resolution = int(resolution)

    axis = _midpoint_axis(half_width, resolution)
    dx = 2.0 * half_width / resolution
    xx, yy = np.meshgrid(axis, axis, indexing="ij")
    radius = np.sqrt(xx * xx + yy * yy).reshape(-1)

    radial_samples: list[np.ndarray] = []
    axial_samples: list[np.ndarray] = []
    energy_weights: list[np.ndarray] = []
    weighted_r2 = 0.0
    weighted_z2 = 0.0
    weight_sum = 0.0

    for z_value in axis:
        zz = np.full_like(xx, z_value)
        points = np.stack((xx, yy, zz), axis=-1).reshape(-1, 3)
        velocity = np.asarray(candidate.at_points(points, time), dtype=float)
        if velocity.shape != points.shape or not np.all(np.isfinite(velocity)):
            raise RuntimeError("candidate returned malformed/nonfinite velocity")
        weight = np.sum(velocity * velocity, axis=-1)
        if np.any(weight < 0.0) or not np.all(np.isfinite(weight)):
            raise RuntimeError("kinetic-energy weights became invalid")
        radial_samples.append(radius)
        axial_samples.append(np.full(radius.shape, abs(float(z_value)), dtype=float))
        energy_weights.append(weight)
        weight_sum += float(np.sum(weight))
        weighted_r2 += float(np.dot(weight, radius * radius))
        weighted_z2 += float(np.sum(weight)) * float(z_value * z_value)

    if not np.isfinite(weight_sum) or weight_sum <= 0.0:
        raise RuntimeError("candidate must have positive finite sampled kinetic energy")

    radial = np.concatenate(radial_samples)
    axial = np.concatenate(axial_samples)
    weights = np.concatenate(energy_weights)
    radial_rms = float(np.sqrt(weighted_r2 / weight_sum))
    axial_rms = float(np.sqrt(weighted_z2 / weight_sum))
    aspect = axial_rms / max(radial_rms, np.finfo(float).tiny)

    return {
        "half_width": half_width,
        "resolution": resolution,
        "cell_width": dx,
        "total_energy": 0.5 * weight_sum * dx**3,
        "radial_q50": _weighted_quantile(radial, weights, 0.50),
        "radial_q90": _weighted_quantile(radial, weights, 0.90),
        "axial_abs_q50": _weighted_quantile(axial, weights, 0.50),
        "axial_abs_q90": _weighted_quantile(axial, weights, 0.90),
        "radial_rms": radial_rms,
        "axial_rms": axial_rms,
        "rms_aspect_z_over_r": aspect,
    }


def _relative_change(after: float, before: float) -> float:
    return float((after - before) / max(abs(before), np.finfo(float).tiny))


def _relative_delta(a: float, b: float) -> float:
    return float(abs(a - b) / max(abs(b), np.finfo(float).tiny))


def audit_optimized_eq45_energy_concentration(
    *,
    times: Iterable[float] = DEFAULT_TIMES,
    half_width: float = DEFAULT_HALF_WIDTH,
    resolutions: Iterable[int] = DEFAULT_RESOLUTIONS,
    candidate_path=None,
    optimization_path=None,
) -> dict:
    """Compare canonical/optimized energy-support geometry over >=3 resolutions."""
    checked_times = _validated_times(times)
    checked_resolutions = _validated_resolutions(resolutions)
    half_width = float(half_width)
    if not np.isfinite(half_width) or half_width <= 0.0:
        raise ValueError("half_width must be positive and finite")

    canonical, optimized, optimization = _candidate_pair(candidate_path, optimization_path)
    reports: list[dict] = []
    for time in checked_times:
        levels: list[dict] = []
        for resolution in checked_resolutions:
            levels.append(
                {
                    "resolution": resolution,
                    "canonical": finite_window_energy_concentration(
                        canonical, time, half_width=half_width, resolution=resolution
                    ),
                    "optimized": finite_window_energy_concentration(
                        optimized, time, half_width=half_width, resolution=resolution
                    ),
                }
            )
        reports.append({"time": time, "levels": levels})

    finest_comparison: list[dict] = []
    resolution_stability: list[dict] = []
    bulk_metrics = ("radial_rms", "axial_rms", "rms_aspect_z_over_r")
    quantile_metrics = ("radial_q50", "radial_q90", "axial_abs_q50", "axial_abs_q90")
    for report in reports:
        previous = report["levels"][-2]
        finest = report["levels"][-1]
        c = finest["canonical"]
        o = finest["optimized"]
        finest_comparison.append(
            {
                "time": report["time"],
                "resolution": finest["resolution"],
                "canonical": c,
                "optimized": o,
                "radial_rms_relative_change": _relative_change(o["radial_rms"], c["radial_rms"]),
                "axial_rms_relative_change": _relative_change(o["axial_rms"], c["axial_rms"]),
                "rms_aspect_relative_change": _relative_change(
                    o["rms_aspect_z_over_r"], c["rms_aspect_z_over_r"]
                ),
            }
        )
        stability = {"time": report["time"]}
        for label in ("canonical", "optimized"):
            stability[label] = {
                "bulk_relative_delta_to_finest": {
                    metric: _relative_delta(previous[label][metric], finest[label][metric])
                    for metric in bulk_metrics
                },
                "quantile_relative_delta_to_finest": {
                    metric: _relative_delta(previous[label][metric], finest[label][metric])
                    for metric in quantile_metrics
                },
            }
        resolution_stability.append(stability)

    max_bulk_delta = max(
        delta
        for row in resolution_stability
        for label in ("canonical", "optimized")
        for delta in row[label]["bulk_relative_delta_to_finest"].values()
    )
    max_quantile_delta = max(
        delta
        for row in resolution_stability
        for label in ("canonical", "optimized")
        for delta in row[label]["quantile_relative_delta_to_finest"].values()
    )

    return {
        "schema": "eq45_optimized_energy_concentration_audit_v1",
        "task_id": "CR007-EQ45-OPTIMIZED-ENERGY-CONCENTRATION-016",
        "canonical_candidate_sha256": canonical.sha256,
        "optimized_candidate_sha256": optimized.sha256,
        "source_optimization_task_id": optimization["task_id"],
        "times": list(checked_times),
        "half_width": half_width,
        "resolutions": list(checked_resolutions),
        "weight_definition": "kinetic-energy density proportional to |u|^2 on a Cartesian midpoint grid",
        "reports": reports,
        "finest_comparison": finest_comparison,
        "resolution_stability": resolution_stability,
        "max_bulk_rms_relative_delta_penultimate_to_finest": max_bulk_delta,
        "max_quantile_relative_delta_penultimate_to_finest": max_quantile_delta,
        "interpretation": (
            "finite-window visualization/research fingerprint; RMS spatial moments are "
            "the primary smooth concentration metric, while grid-quantized q50/q90 are "
            "reported with their explicit resolution sensitivity"
        ),
        "velocity_changed_by_source_optimization": True,
        "velocity_changed_by_this_audit": False,
        "physical_support_validated": False,
        "cr001_energy_normalization_status": "pending_unknown",
        "visualization_ready": False,
        "visual_correspondence_verified": False,
        "pde_validated": False,
        "paper_exact": False,
        "openai_field_identified": False,
        "blowup_proved": False,
    }
