"""Three-resolution whole-domain morphology audit for the early-localized Phi10 trial.

The early-localized quadratic schedule is already identical to the selected affine
Phi(1,0) trial at the delivery start and identical to the static supported field
at the midpoint and delivery end.  Its only new visual question is therefore the
intermediate excursion.  This module measures that excursion at ``t=0.625`` on
three structured grids without changing the candidate or using any public image
as an optimization target.

All vorticity is reconstructed independently from public
``at_points(...)->[u,v,w]`` samples.  A single threshold derived from the static
supported baseline is applied unchanged to the static, quadratic and affine
fields at each resolution.
"""
from __future__ import annotations

from collections.abc import Iterable

import numpy as np

from .constrained_eq45_supported_candidate import Eq45SupportedVelocityCandidate
from .constrained_eq45_supported_collar_vorticity import governed_supported_seed
from .constrained_eq45_supported_phi10_early_localized_temporal_mode import (
    Eq45SupportedPhi10EarlyLocalizedTemporalCandidate,
)
from .constrained_eq45_supported_phi10_temporal_envelope import _metrics
from .constrained_eq45_supported_phi10_temporal_replay import build_phi10_temporal_trial
from .constrained_eq45_supported_vorticity_envelope import structured_vorticity_magnitude


DEFAULT_RESOLUTIONS = (49, 65, 81)
DEFAULT_TIME = 0.625
DEFAULT_AFFINE_SLOPE = 1.4
DEFAULT_EARLY_DELTA = -1.4
DEFAULT_HALF_WIDTH = 2.0
DEFAULT_SUPERLEVEL_FRACTION = 0.25
DEFAULT_PLATEAU_RADIUS = 1.5
DEFAULT_PLATEAU_HALF_HEIGHT = 1.5
DEFAULT_COLLAR_START = 1.6

_TRUTH_BOUNDARY = {
    "canonical_velocity_changed": False,
    "diagnostic_velocity_changed": False,
    "visualization_candidate_only": True,
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


def _collateral_fraction(quadratic: float, baseline: float, affine: float) -> float:
    denominator = abs(float(affine) - float(baseline))
    if denominator <= 1.0e-300:
        return 0.0 if abs(float(quadratic) - float(baseline)) <= 1.0e-300 else float("inf")
    return float(abs(float(quadratic) - float(baseline)) / denominator)


def audit_early_localized_intermediate_morphology(
    base: Eq45SupportedVelocityCandidate | None = None,
    *,
    resolutions: Iterable[int] = DEFAULT_RESOLUTIONS,
    time: float = DEFAULT_TIME,
    affine_slope: float = DEFAULT_AFFINE_SLOPE,
    early_delta: float = DEFAULT_EARLY_DELTA,
    half_width: float = DEFAULT_HALF_WIDTH,
    superlevel_fraction: float = DEFAULT_SUPERLEVEL_FRACTION,
    plateau_radius: float = DEFAULT_PLATEAU_RADIUS,
    plateau_half_height: float = DEFAULT_PLATEAU_HALF_HEIGHT,
    collar_start: float = DEFAULT_COLLAR_START,
) -> dict[str, object]:
    """Compare static, early-localized quadratic and affine fields at one time."""
    child = governed_supported_seed() if base is None else base
    if not isinstance(child, Eq45SupportedVelocityCandidate):
        raise TypeError("base must be Eq45SupportedVelocityCandidate")

    resolution_values = tuple(int(value) for value in resolutions)
    if len(resolution_values) < 3 or any(value < 9 for value in resolution_values):
        raise ValueError("resolutions must contain at least three levels >= 9")
    if not all(b > a for a, b in zip(resolution_values, resolution_values[1:])):
        raise ValueError("resolutions must be strictly increasing")

    scalar_time = float(time)
    if not np.isfinite(scalar_time) or not (child.time_start < scalar_time < child.time_end):
        raise ValueError("time must be finite and strictly inside the delivery interval")
    if np.isclose(scalar_time, child.time_midpoint if hasattr(child, "time_midpoint") else 0.5 * (child.time_start + child.time_end), rtol=0.0, atol=1.0e-15):
        raise ValueError("time must not be the static midpoint")
    if not np.isfinite(affine_slope) or float(affine_slope) == 0.0:
        raise ValueError("affine_slope must be finite and nonzero")
    if not np.isfinite(early_delta) or float(early_delta) == 0.0:
        raise ValueError("early_delta must be finite and nonzero")

    for value, name in (
        (half_width, "half_width"),
        (plateau_radius, "plateau_radius"),
        (plateau_half_height, "plateau_half_height"),
        (collar_start, "collar_start"),
    ):
        if not np.isfinite(value) or float(value) <= 0.0:
            raise ValueError(f"{name} must be positive and finite")
    if not np.isfinite(superlevel_fraction) or not (0.0 < superlevel_fraction < 1.0):
        raise ValueError("superlevel_fraction must lie in (0,1)")

    quadratic = Eq45SupportedPhi10EarlyLocalizedTemporalCandidate(
        base=child, early_delta=float(early_delta)
    )
    affine = build_phi10_temporal_trial(child, slope=float(affine_slope))

    rows: list[dict[str, object]] = []
    for resolution in resolution_values:
        axis = np.linspace(-float(half_width), float(half_width), resolution)
        baseline_mag, radius, abs_z = structured_vorticity_magnitude(child, axis, scalar_time)
        quadratic_mag, quadratic_radius, quadratic_abs_z = structured_vorticity_magnitude(
            quadratic, axis, scalar_time
        )
        affine_mag, affine_radius, affine_abs_z = structured_vorticity_magnitude(
            affine, axis, scalar_time
        )
        np.testing.assert_allclose(quadratic_radius, radius, rtol=0.0, atol=0.0)
        np.testing.assert_allclose(quadratic_abs_z, abs_z, rtol=0.0, atol=0.0)
        np.testing.assert_allclose(affine_radius, radius, rtol=0.0, atol=0.0)
        np.testing.assert_allclose(affine_abs_z, abs_z, rtol=0.0, atol=0.0)

        plateau = (radius <= float(plateau_radius)) & (abs_z <= float(plateau_half_height))
        core_peak = float(np.max(baseline_mag[plateau]))
        threshold = float(superlevel_fraction * core_peak)
        rows.append(
            {
                "resolution": resolution,
                "dx": float(axis[1] - axis[0]),
                "time": scalar_time,
                "baseline_core_peak": core_peak,
                "shared_threshold": threshold,
                "baseline": _metrics(
                    baseline_mag,
                    radius,
                    abs_z,
                    threshold=threshold,
                    collar_start=collar_start,
                ),
                "quadratic": _metrics(
                    quadratic_mag,
                    radius,
                    abs_z,
                    threshold=threshold,
                    collar_start=collar_start,
                ),
                "affine": _metrics(
                    affine_mag,
                    radius,
                    abs_z,
                    threshold=threshold,
                    collar_start=collar_start,
                ),
            }
        )

    fine = rows[-1]
    middle = rows[-2]
    baseline = fine["baseline"]
    quadratic_metrics = fine["quadratic"]
    affine_metrics = fine["affine"]
    quadratic_middle = middle["quadratic"]

    compare_keys = (
        "radial_q99",
        "axial_q99",
        "aspect_q99",
        "vorticity2_weighted_radial_rms",
        "vorticity2_weighted_axial_rms",
        "vorticity2_weighted_aspect",
        "radial_outer_vorticity2_fraction",
        "whole_grid_collar_vorticity2_fraction",
    )
    summary: dict[str, float] = {}
    for key in compare_keys:
        summary[f"quadratic_baseline_{key}_ratio"] = _ratio(
            quadratic_metrics[key], baseline[key]
        )
        summary[f"affine_baseline_{key}_ratio"] = _ratio(affine_metrics[key], baseline[key])
        summary[f"quadratic_affine_{key}_collateral_fraction"] = _collateral_fraction(
            quadratic_metrics[key], baseline[key], affine_metrics[key]
        )

    for key in (
        "radial_q99",
        "axial_q99",
        "vorticity2_weighted_radial_rms",
        "vorticity2_weighted_axial_rms",
        "vorticity2_weighted_aspect",
        "radial_outer_vorticity2_fraction",
        "whole_grid_collar_vorticity2_fraction",
    ):
        summary[f"quadratic_{key}_mid_to_fine_relative_change"] = _relative_change(
            quadratic_metrics[key], quadratic_middle[key]
        )

    return {
        "schema": "eq45_supported_phi10_early_localized_intermediate_morphology_v1",
        "claim_scope": "target_free_three_resolution_intermediate_visualization_morphology_diagnostic",
        "base_supported_sha256": child.sha256,
        "quadratic_candidate_sha256": quadratic.sha256,
        "affine_candidate_sha256": affine.sha256,
        "mode": {"family": "phi", "index": [1, 0]},
        "time": scalar_time,
        "tau": float(quadratic.tau(scalar_time)),
        "coefficients": {
            "static": float(quadratic.midpoint_coefficient),
            "quadratic": float(quadratic.coefficient_at(scalar_time)),
            "affine": float(affine.coefficient_at(scalar_time)),
        },
        "quadratic_early_delta": float(early_delta),
        "affine_slope": float(affine_slope),
        "resolutions": list(resolution_values),
        "half_width": float(half_width),
        "superlevel_fraction": float(superlevel_fraction),
        "threshold_reference": "same static-supported inner-core peak per resolution",
        "curl_method": "numpy_second_order_structured_cartesian_gradient_from_public_velocity",
        "collar_start": float(collar_start),
        "rows": rows,
        "finest_summary": summary,
        "truth_boundary": dict(_TRUTH_BOUNDARY),
    }
