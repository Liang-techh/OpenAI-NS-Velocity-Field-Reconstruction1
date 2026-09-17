"""Three-resolution vorticity morphology audit for the bounded axial taper family.

Agent 1 materialized a single support-localized control that changes only the
axial identity plateau of the existing divergence-preserving physical taper.
Agent 2 found that this control has substantial local public-velocity leverage
with essentially flat fixed-probe pressure-free vorticity obstruction, while
Agent 3 found no non-convergent divergence defect.  This module asks the next
visualization-facing question without selecting a taper value: does moving the
identity plateau actually change the resolved whole-domain vortex envelope or
mainly alter low-vorticity velocity in the outer axial collar?

All vorticity is reconstructed independently from public
``at_points(...)->[u,v,w]`` samples.  At each grid resolution one threshold is
frozen from the baseline q=0.64 inner-core vorticity peak and then applied
unchanged to every taper value.  No public image, force, pressure, PDE residual,
or hidden OpenAI parameter enters this diagnostic.
"""
from __future__ import annotations

from collections.abc import Iterable

import numpy as np

from .constrained_eq45_supported_collar_vorticity import governed_supported_seed
from .constrained_eq45_supported_phi10_blend_axial_taper_mode import (
    BASE_AXIAL_PLATEAU_Q,
    MAX_AXIAL_PLATEAU_Q,
    Eq45SupportedPhi10BlendAxialTaperCandidate,
)
from .constrained_eq45_supported_phi10_compact_quartic_blend_temporal_mode import (
    Eq45SupportedPhi10CompactQuarticBlendTemporalCandidate,
)
from .constrained_eq45_supported_vorticity_envelope import (
    structured_vorticity_magnitude,
)


TASK_ID = "CR008-EQ45-SUPPORTED-BLEND-AXIAL-TAPER-VORTICITY-MORPHOLOGY-024"
DEFAULT_RESOLUTIONS = (49, 65, 81)
DEFAULT_AXIAL_PLATEAU_Q_VALUES = (0.64, 0.7225, 0.81)
DEFAULT_BLEND_WEIGHT = 0.5
DEFAULT_EARLY_DELTA = -1.4
DEFAULT_TIME = 0.3125
DEFAULT_HALF_WIDTH = 2.0
DEFAULT_SUPERLEVEL_FRACTION = 0.25
DEFAULT_PLATEAU_RADIUS = 1.5
DEFAULT_PLATEAU_HALF_HEIGHT = 1.5
DEFAULT_FIXED_COLLAR_START = 1.6

_TRUTH_BOUNDARY = {
    "velocity_changed_by_diagnostic": False,
    "visualization_candidate_only": True,
    "axial_taper_value_selected": False,
    "canonical_velocity_changed": False,
    "force_or_pressure_fitted": False,
    "pde_objective_used_to_choose_axial_taper": False,
    "public_image_fitted": False,
    "physical_support_validated": False,
    "visualization_ready": False,
    "visual_correspondence_verified": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


def _relative_change(fine: float, middle: float) -> float:
    return float(abs(float(fine) - float(middle)) / max(abs(float(fine)), 1.0e-300))


def _ratio(value: float, baseline: float) -> float:
    return float(float(value) / max(abs(float(baseline)), 1.0e-300))


def _morphology_metrics(
    magnitude: np.ndarray,
    radius: np.ndarray,
    abs_z: np.ndarray,
    *,
    threshold: float,
    fixed_collar_start: float,
) -> dict[str, float | int]:
    active = magnitude >= float(threshold)
    if not np.any(active):
        raise RuntimeError("vorticity superlevel set is empty")

    active_r = radius[active]
    active_z = abs_z[active]
    weight = magnitude * magnitude
    total_weight = float(np.sum(weight))
    if not np.isfinite(total_weight) or total_weight <= 1.0e-300:
        raise RuntimeError("vorticity-squared weight is zero or nonfinite")

    radial_rms = float(np.sqrt(np.sum(weight * radius * radius) / total_weight))
    axial_rms = float(np.sqrt(np.sum(weight * abs_z * abs_z) / total_weight))
    radial_outer = radius > float(fixed_collar_start)
    axial_outer = abs_z > float(fixed_collar_start)
    whole_collar = radial_outer | axial_outer

    return {
        "active_voxels": int(np.count_nonzero(active)),
        "radial_q99": float(np.quantile(active_r, 0.99)),
        "axial_q99": float(np.quantile(active_z, 0.99)),
        "radial_max": float(np.max(active_r)),
        "axial_max": float(np.max(active_z)),
        "aspect_q99": float(
            np.quantile(active_z, 0.99) / max(np.quantile(active_r, 0.99), 1.0e-300)
        ),
        "vorticity2_weighted_radial_rms": radial_rms,
        "vorticity2_weighted_axial_rms": axial_rms,
        "vorticity2_weighted_aspect": axial_rms / max(radial_rms, 1.0e-300),
        "radial_outer_vorticity2_fraction": float(np.sum(weight[radial_outer]) / total_weight),
        "axial_outer_vorticity2_fraction": float(np.sum(weight[axial_outer]) / total_weight),
        "whole_grid_collar_vorticity2_fraction": float(
            np.sum(weight[whole_collar]) / total_weight
        ),
        "superlevel_axial_collar_voxel_fraction": float(np.mean(axial_outer[active])),
        "superlevel_whole_collar_voxel_fraction": float(np.mean(whole_collar[active])),
    }


def _build_family(
    *,
    blend_weight: float,
    early_delta: float,
    axial_plateau_q: float,
) -> Eq45SupportedPhi10BlendAxialTaperCandidate:
    base = Eq45SupportedPhi10CompactQuarticBlendTemporalCandidate(
        base=governed_supported_seed(),
        blend_weight=float(blend_weight),
        early_delta=float(early_delta),
    )
    return Eq45SupportedPhi10BlendAxialTaperCandidate(
        base=base,
        axial_plateau_q=float(axial_plateau_q),
    )


def audit_axial_taper_vorticity_morphology(
    *,
    resolutions: Iterable[int] = DEFAULT_RESOLUTIONS,
    axial_plateau_q_values: Iterable[float] = DEFAULT_AXIAL_PLATEAU_Q_VALUES,
    blend_weight: float = DEFAULT_BLEND_WEIGHT,
    early_delta: float = DEFAULT_EARLY_DELTA,
    time: float = DEFAULT_TIME,
    half_width: float = DEFAULT_HALF_WIDTH,
    superlevel_fraction: float = DEFAULT_SUPERLEVEL_FRACTION,
    plateau_radius: float = DEFAULT_PLATEAU_RADIUS,
    plateau_half_height: float = DEFAULT_PLATEAU_HALF_HEIGHT,
    fixed_collar_start: float = DEFAULT_FIXED_COLLAR_START,
) -> dict[str, object]:
    """Audit resolved vorticity geometry across a frozen axial-taper grid."""
    resolution_values = tuple(int(value) for value in resolutions)
    if len(resolution_values) < 3 or any(value < 9 for value in resolution_values):
        raise ValueError("resolutions must contain at least three levels >= 9")
    if not all(b > a for a, b in zip(resolution_values, resolution_values[1:])):
        raise ValueError("resolutions must be strictly increasing")

    q_values = tuple(float(value) for value in axial_plateau_q_values)
    if len(q_values) < 2 or not all(np.isfinite(value) for value in q_values):
        raise ValueError("axial_plateau_q_values must contain at least two finite values")
    if q_values != tuple(sorted(q_values)) or len(set(q_values)) != len(q_values):
        raise ValueError("axial_plateau_q_values must be unique and strictly increasing")
    if q_values[0] != BASE_AXIAL_PLATEAU_Q:
        raise ValueError("axial taper grid must begin at the governed q=0.64 baseline")
    if q_values[-1] > MAX_AXIAL_PLATEAU_Q:
        raise ValueError("axial taper grid exceeds the autonomous q bound")

    weight = float(blend_weight)
    early = float(early_delta)
    scalar_time = float(time)
    if not np.isfinite(weight) or not (0.0 <= weight <= 1.0):
        raise ValueError("blend_weight must be finite and lie in [0,1]")
    if not np.isfinite(early) or early == 0.0:
        raise ValueError("early_delta must be finite and nonzero")
    if not np.isfinite(scalar_time):
        raise ValueError("time must be finite")

    for value, name in (
        (half_width, "half_width"),
        (plateau_radius, "plateau_radius"),
        (plateau_half_height, "plateau_half_height"),
        (fixed_collar_start, "fixed_collar_start"),
    ):
        if not np.isfinite(value) or float(value) <= 0.0:
            raise ValueError(f"{name} must be positive and finite")
    if not np.isfinite(superlevel_fraction) or not (0.0 < superlevel_fraction < 1.0):
        raise ValueError("superlevel_fraction must lie in (0,1)")

    candidates = {
        q: _build_family(
            blend_weight=weight,
            early_delta=early,
            axial_plateau_q=q,
        )
        for q in q_values
    }
    baseline_candidate = candidates[q_values[0]]
    if not (baseline_candidate.time_start <= scalar_time <= baseline_candidate.time_end):
        raise ValueError("time must lie inside the candidate delivery interval")

    rows: list[dict[str, object]] = []
    for resolution in resolution_values:
        axis = np.linspace(-float(half_width), float(half_width), resolution)
        baseline_mag, radius, abs_z = structured_vorticity_magnitude(
            baseline_candidate, axis, scalar_time
        )
        plateau = (radius <= float(plateau_radius)) & (
            abs_z <= float(plateau_half_height)
        )
        core_peak = float(np.max(baseline_mag[plateau]))
        threshold = float(superlevel_fraction * core_peak)

        for q in q_values:
            candidate = candidates[q]
            if q == q_values[0]:
                magnitude = baseline_mag
                candidate_radius = radius
                candidate_abs_z = abs_z
            else:
                magnitude, candidate_radius, candidate_abs_z = structured_vorticity_magnitude(
                    candidate, axis, scalar_time
                )
                np.testing.assert_allclose(candidate_radius, radius, rtol=0.0, atol=0.0)
                np.testing.assert_allclose(candidate_abs_z, abs_z, rtol=0.0, atol=0.0)

            rows.append(
                {
                    "resolution": resolution,
                    "dx": float(axis[1] - axis[0]),
                    "time": scalar_time,
                    "axial_plateau_q": q,
                    "axial_identity_half_height": float(candidate.axial_identity_half_height),
                    "baseline_core_peak": core_peak,
                    "shared_threshold": threshold,
                    "metrics": _morphology_metrics(
                        magnitude,
                        radius,
                        abs_z,
                        threshold=threshold,
                        fixed_collar_start=fixed_collar_start,
                    ),
                }
            )

    fine_resolution = resolution_values[-1]
    middle_resolution = resolution_values[-2]
    baseline_fine = next(
        row
        for row in rows
        if row["resolution"] == fine_resolution
        and row["axial_plateau_q"] == q_values[0]
    )
    baseline_metrics = baseline_fine["metrics"]

    ratio_keys = (
        "radial_q99",
        "axial_q99",
        "aspect_q99",
        "vorticity2_weighted_radial_rms",
        "vorticity2_weighted_axial_rms",
        "vorticity2_weighted_aspect",
        "radial_outer_vorticity2_fraction",
        "axial_outer_vorticity2_fraction",
        "whole_grid_collar_vorticity2_fraction",
    )
    sensitivity_keys = (
        "radial_q99",
        "axial_q99",
        "vorticity2_weighted_radial_rms",
        "vorticity2_weighted_axial_rms",
        "vorticity2_weighted_aspect",
        "axial_outer_vorticity2_fraction",
        "whole_grid_collar_vorticity2_fraction",
    )

    summaries: list[dict[str, float]] = []
    for q in q_values:
        fine = next(
            row
            for row in rows
            if row["resolution"] == fine_resolution and row["axial_plateau_q"] == q
        )
        middle = next(
            row
            for row in rows
            if row["resolution"] == middle_resolution and row["axial_plateau_q"] == q
        )
        fine_metrics = fine["metrics"]
        middle_metrics = middle["metrics"]
        summary: dict[str, float] = {
            "axial_plateau_q": float(q),
            "axial_identity_half_height": float(fine["axial_identity_half_height"]),
            "superlevel_axial_collar_voxel_fraction": float(
                fine_metrics["superlevel_axial_collar_voxel_fraction"]
            ),
            "superlevel_whole_collar_voxel_fraction": float(
                fine_metrics["superlevel_whole_collar_voxel_fraction"]
            ),
        }
        for key in ratio_keys:
            summary[f"baseline_{key}_ratio"] = _ratio(
                fine_metrics[key], baseline_metrics[key]
            )
        for key in sensitivity_keys:
            summary[f"{key}_mid_to_fine_relative_change"] = _relative_change(
                fine_metrics[key], middle_metrics[key]
            )
        summaries.append(summary)

    return {
        "schema": "eq45_supported_phi10_blend_axial_taper_vorticity_morphology_v1",
        "task_id": TASK_ID,
        "claim_scope": "target_free_three_resolution_axial_taper_visualization_morphology_diagnostic",
        "base_supported_sha256": baseline_candidate.supported_base_sha256,
        "base_blend_sha256": baseline_candidate.base_sha256,
        "blend_weight": weight,
        "early_delta": early,
        "time": scalar_time,
        "resolutions": list(resolution_values),
        "axial_plateau_q_values": list(q_values),
        "half_width": float(half_width),
        "superlevel_fraction": float(superlevel_fraction),
        "threshold_reference": "same q=0.64 inner-core peak per resolution",
        "curl_method": "numpy_second_order_structured_cartesian_gradient_from_public_velocity",
        "fixed_collar_start": float(fixed_collar_start),
        "fixed_collar_semantics": "baseline physical identity plateau boundary; held fixed for every q",
        "rows": rows,
        "finest_summaries": summaries,
        "truth_boundary": dict(_TRUTH_BOUNDARY),
    }
