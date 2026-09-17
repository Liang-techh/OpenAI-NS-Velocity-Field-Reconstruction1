"""Resolution audit for the selected supported Eq45 Phi(1,0) temporal trial.

Agent 1's nonlinear replay selects slope 1.4 because it strongly reduces one
resolved early radial-outer vorticity feature on a 41^3 screen.  This module
asks the independent Agent-4 question before that trial is treated as a useful
visualization candidate: does the changed public ``[u,v,w]`` produce a stable
whole-domain vorticity envelope across finer grids, and what happens to radial
versus axial extent?

The comparison is target-free.  It does not use an OpenAI image and it does not
change or promote the velocity field.  At each time and resolution a threshold
is fixed from the *baseline supported child* inner-core peak, then applied
unchanged to both baseline and temporal trial.  All curl quantities are rebuilt
from public ``at_points(...)->[u,v,w]`` samples.
"""
from __future__ import annotations

from collections.abc import Iterable

import numpy as np

from .constrained_eq45_supported_candidate import Eq45SupportedVelocityCandidate
from .constrained_eq45_supported_collar_vorticity import governed_supported_seed
from .constrained_eq45_supported_phi10_temporal_replay import build_phi10_temporal_trial
from .constrained_eq45_supported_vorticity_envelope import structured_vorticity_magnitude


DEFAULT_SLOPE = 1.4
DEFAULT_RESOLUTIONS = (49, 65, 81)
DEFAULT_TIMES = (0.25, 0.50, 0.75)
DEFAULT_HALF_WIDTH = 2.0
DEFAULT_SUPERLEVEL_FRACTION = 0.25
DEFAULT_PLATEAU_RADIUS = 1.5
DEFAULT_PLATEAU_HALF_HEIGHT = 1.5
DEFAULT_COLLAR_START = 1.6

_TRUTH_BOUNDARY = {
    "canonical_velocity_changed": False,
    "trial_velocity_changed": True,
    "production_slope_promoted": False,
    "openai_image_fitted": False,
    "physical_support_validated": False,
    "visualization_ready": False,
    "visual_correspondence_verified": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


def _metrics(
    magnitude: np.ndarray,
    radius: np.ndarray,
    abs_z: np.ndarray,
    *,
    threshold: float,
    collar_start: float,
) -> dict[str, float | int]:
    active = magnitude >= float(threshold)
    if not np.any(active):
        raise RuntimeError("vorticity superlevel set is empty")

    active_r = radius[active]
    active_z = abs_z[active]
    radial_q99 = float(np.quantile(active_r, 0.99))
    axial_q99 = float(np.quantile(active_z, 0.99))

    weight = magnitude * magnitude
    total = float(np.sum(weight))
    if not np.isfinite(total) or total <= 1.0e-300:
        raise RuntimeError("vorticity-squared weight is zero or nonfinite")
    radial_rms = float(np.sqrt(np.sum(weight * radius * radius) / total))
    axial_rms = float(np.sqrt(np.sum(weight * abs_z * abs_z) / total))

    radial_outer = radius > float(collar_start)
    collar = radial_outer | (abs_z > float(collar_start))
    return {
        "active_voxels": int(np.count_nonzero(active)),
        "radial_q99": radial_q99,
        "axial_q99": axial_q99,
        "aspect_q99": axial_q99 / max(radial_q99, 1.0e-300),
        "radial_max": float(np.max(active_r)),
        "axial_max": float(np.max(active_z)),
        "vorticity2_weighted_radial_rms": radial_rms,
        "vorticity2_weighted_axial_rms": axial_rms,
        "vorticity2_weighted_aspect": axial_rms / max(radial_rms, 1.0e-300),
        "radial_outer_vorticity2_fraction": float(np.sum(weight[radial_outer]) / total),
        "whole_grid_collar_vorticity2_fraction": float(np.sum(weight[collar]) / total),
        "superlevel_collar_voxel_fraction": float(np.mean(collar[active])),
    }


def _relative_change(fine: float, middle: float) -> float:
    return float(abs(float(fine) - float(middle)) / max(abs(float(fine)), 1.0e-300))


def audit_supported_phi10_temporal_envelope(
    base: Eq45SupportedVelocityCandidate | None = None,
    *,
    slope: float = DEFAULT_SLOPE,
    resolutions: Iterable[int] = DEFAULT_RESOLUTIONS,
    times: Iterable[float] = DEFAULT_TIMES,
    half_width: float = DEFAULT_HALF_WIDTH,
    superlevel_fraction: float = DEFAULT_SUPERLEVEL_FRACTION,
    plateau_radius: float = DEFAULT_PLATEAU_RADIUS,
    plateau_half_height: float = DEFAULT_PLATEAU_HALF_HEIGHT,
    collar_start: float = DEFAULT_COLLAR_START,
) -> dict[str, object]:
    """Compare the selected temporal trial with its static supported baseline."""
    child = governed_supported_seed() if base is None else base
    if not isinstance(child, Eq45SupportedVelocityCandidate):
        raise TypeError("base must be Eq45SupportedVelocityCandidate")
    if not np.isfinite(slope) or float(slope) == 0.0:
        raise ValueError("slope must be finite and nonzero")

    resolution_values = tuple(int(value) for value in resolutions)
    if len(resolution_values) < 3 or any(value < 9 for value in resolution_values):
        raise ValueError("resolutions must contain at least three levels >= 9")
    if not all(b > a for a, b in zip(resolution_values, resolution_values[1:])):
        raise ValueError("resolutions must be strictly increasing")

    time_values = tuple(float(value) for value in times)
    if len(time_values) != 3 or not all(np.isfinite(value) for value in time_values):
        raise ValueError("times must contain exactly three finite values")
    if time_values != tuple(sorted(time_values)):
        raise ValueError("times must be increasing")
    if time_values[0] < child.time_start or time_values[-1] > child.time_end:
        raise ValueError("times must lie inside the candidate delivery interval")
    midpoint = 0.5 * (child.time_start + child.time_end)
    if not np.isclose(time_values[1], midpoint, rtol=0.0, atol=1.0e-15):
        raise ValueError("middle time must equal the candidate midpoint")

    for value, name in (
        (half_width, "half_width"),
        (plateau_radius, "plateau_radius"),
        (plateau_half_height, "plateau_half_height"),
        (collar_start, "collar_start"),
    ):
        if not np.isfinite(value) or value <= 0.0:
            raise ValueError(f"{name} must be positive and finite")
    if not np.isfinite(superlevel_fraction) or not (0.0 < superlevel_fraction < 1.0):
        raise ValueError("superlevel_fraction must lie in (0,1)")

    trial = build_phi10_temporal_trial(child, slope=float(slope))
    rows: list[dict[str, object]] = []
    for resolution in resolution_values:
        axis = np.linspace(-float(half_width), float(half_width), resolution)
        for time in time_values:
            baseline_mag, radius, abs_z = structured_vorticity_magnitude(child, axis, time)
            trial_mag, trial_radius, trial_abs_z = structured_vorticity_magnitude(trial, axis, time)
            plateau = (radius <= float(plateau_radius)) & (abs_z <= float(plateau_half_height))
            core_peak = float(np.max(baseline_mag[plateau]))
            threshold = float(superlevel_fraction * core_peak)
            rows.append(
                {
                    "resolution": resolution,
                    "dx": float(axis[1] - axis[0]),
                    "time": time,
                    "baseline_core_peak": core_peak,
                    "shared_threshold": threshold,
                    "baseline": _metrics(
                        baseline_mag,
                        radius,
                        abs_z,
                        threshold=threshold,
                        collar_start=collar_start,
                    ),
                    "trial": _metrics(
                        trial_mag,
                        trial_radius,
                        trial_abs_z,
                        threshold=threshold,
                        collar_start=collar_start,
                    ),
                }
            )

    summaries: list[dict[str, float]] = []
    for time in time_values:
        fine = next(
            row for row in rows
            if row["resolution"] == resolution_values[-1] and row["time"] == time
        )
        middle = next(
            row for row in rows
            if row["resolution"] == resolution_values[-2] and row["time"] == time
        )
        baseline = fine["baseline"]
        trial_metrics = fine["trial"]
        trial_middle = middle["trial"]
        summaries.append(
            {
                "time": time,
                "trial_baseline_radial_q99_ratio": float(
                    trial_metrics["radial_q99"] / max(baseline["radial_q99"], 1.0e-300)
                ),
                "trial_baseline_axial_q99_ratio": float(
                    trial_metrics["axial_q99"] / max(baseline["axial_q99"], 1.0e-300)
                ),
                "trial_baseline_aspect_q99_ratio": float(
                    trial_metrics["aspect_q99"] / max(baseline["aspect_q99"], 1.0e-300)
                ),
                "trial_baseline_radial_rms_ratio": float(
                    trial_metrics["vorticity2_weighted_radial_rms"]
                    / max(baseline["vorticity2_weighted_radial_rms"], 1.0e-300)
                ),
                "trial_baseline_axial_rms_ratio": float(
                    trial_metrics["vorticity2_weighted_axial_rms"]
                    / max(baseline["vorticity2_weighted_axial_rms"], 1.0e-300)
                ),
                "trial_baseline_weighted_aspect_ratio": float(
                    trial_metrics["vorticity2_weighted_aspect"]
                    / max(baseline["vorticity2_weighted_aspect"], 1.0e-300)
                ),
                "trial_baseline_radial_outer_vorticity2_ratio": float(
                    trial_metrics["radial_outer_vorticity2_fraction"]
                    / max(baseline["radial_outer_vorticity2_fraction"], 1.0e-300)
                ),
                "trial_baseline_collar_vorticity2_ratio": float(
                    trial_metrics["whole_grid_collar_vorticity2_fraction"]
                    / max(baseline["whole_grid_collar_vorticity2_fraction"], 1.0e-300)
                ),
                "trial_radial_q99_mid_to_fine_relative_change": _relative_change(
                    trial_metrics["radial_q99"], trial_middle["radial_q99"]
                ),
                "trial_axial_q99_mid_to_fine_relative_change": _relative_change(
                    trial_metrics["axial_q99"], trial_middle["axial_q99"]
                ),
                "trial_radial_rms_mid_to_fine_relative_change": _relative_change(
                    trial_metrics["vorticity2_weighted_radial_rms"],
                    trial_middle["vorticity2_weighted_radial_rms"],
                ),
                "trial_axial_rms_mid_to_fine_relative_change": _relative_change(
                    trial_metrics["vorticity2_weighted_axial_rms"],
                    trial_middle["vorticity2_weighted_axial_rms"],
                ),
            }
        )

    return {
        "schema": "eq45_supported_phi10_temporal_envelope_audit_v1",
        "claim_scope": "target_free_resolution_stable_visualization_morphology_diagnostic",
        "base_supported_sha256": child.sha256,
        "trial_candidate_sha256": trial.sha256,
        "mode": {"family": "phi", "index": [1, 0]},
        "slope": float(slope),
        "resolutions": list(resolution_values),
        "times": list(time_values),
        "half_width": float(half_width),
        "superlevel_fraction": float(superlevel_fraction),
        "threshold_reference": "same baseline-supported inner-core peak per time/resolution",
        "curl_method": "numpy_second_order_structured_cartesian_gradient_from_public_velocity",
        "collar_start": float(collar_start),
        "rows": rows,
        "time_summaries": summaries,
        "truth_boundary": dict(_TRUTH_BOUNDARY),
    }
