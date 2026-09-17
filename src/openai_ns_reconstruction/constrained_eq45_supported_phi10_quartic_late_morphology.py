"""Three-resolution late off-keyframe morphology audit for quartic Phi(1,0).

The derivative-balanced quartic candidate preserves the selected early snapshot
and the static supported snapshots at t=.50/.625/.75.  Its remaining visual
question is whether the stronger temporal localization also reduces the actual
whole-domain 3-D vorticity-geometry collateral between those keyframes.

At t=.6875 this diagnostic compares the frozen static supported field, the cubic
localized field, and the quartic derivative-balanced field on three Cartesian
resolutions.  Vorticity is reconstructed independently from public
``at_points(...)->[u,v,w]`` values and all three fields use one threshold derived
from the static inner-core peak at each resolution.  No public image, PDE
residual, force, pressure, or candidate coefficient is fitted here.
"""
from __future__ import annotations

from collections.abc import Iterable

import numpy as np

from .constrained_eq45_supported_candidate import Eq45SupportedVelocityCandidate
from .constrained_eq45_supported_collar_vorticity import governed_supported_seed
from .constrained_eq45_supported_phi10_cubic_temporal_mode import (
    Eq45SupportedPhi10CubicLocalizedTemporalCandidate,
)
from .constrained_eq45_supported_phi10_quartic_temporal_mode import (
    Eq45SupportedPhi10QuarticDerivativeBalancedTemporalCandidate,
)
from .constrained_eq45_supported_phi10_temporal_envelope import _metrics
from .constrained_eq45_supported_vorticity_envelope import structured_vorticity_magnitude


DEFAULT_RESOLUTIONS = (49, 65, 81)
DEFAULT_TIME = 0.6875
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


def _excursion_fraction(quartic: float, baseline: float, cubic: float) -> float:
    denominator = abs(float(cubic) - float(baseline))
    if denominator <= 1.0e-300:
        return 0.0 if abs(float(quartic) - float(baseline)) <= 1.0e-300 else float("inf")
    return float(abs(float(quartic) - float(baseline)) / denominator)


def audit_quartic_late_offkeyframe_morphology(
    base: Eq45SupportedVelocityCandidate | None = None,
    *,
    resolutions: Iterable[int] = DEFAULT_RESOLUTIONS,
    time: float = DEFAULT_TIME,
    early_delta: float = DEFAULT_EARLY_DELTA,
    half_width: float = DEFAULT_HALF_WIDTH,
    superlevel_fraction: float = DEFAULT_SUPERLEVEL_FRACTION,
    plateau_radius: float = DEFAULT_PLATEAU_RADIUS,
    plateau_half_height: float = DEFAULT_PLATEAU_HALF_HEIGHT,
    collar_start: float = DEFAULT_COLLAR_START,
) -> dict[str, object]:
    """Compare static, cubic and quartic supported fields at one late off-keyframe."""
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

    cubic = Eq45SupportedPhi10CubicLocalizedTemporalCandidate(
        base=child, early_delta=float(early_delta)
    )
    quartic = Eq45SupportedPhi10QuarticDerivativeBalancedTemporalCandidate(
        base=child, early_delta=float(early_delta)
    )

    rows: list[dict[str, object]] = []
    for resolution in resolution_values:
        axis = np.linspace(-float(half_width), float(half_width), resolution)
        baseline_mag, radius, abs_z = structured_vorticity_magnitude(child, axis, scalar_time)
        cubic_mag, cubic_radius, cubic_abs_z = structured_vorticity_magnitude(
            cubic, axis, scalar_time
        )
        quartic_mag, quartic_radius, quartic_abs_z = structured_vorticity_magnitude(
            quartic, axis, scalar_time
        )
        np.testing.assert_allclose(cubic_radius, radius, rtol=0.0, atol=0.0)
        np.testing.assert_allclose(cubic_abs_z, abs_z, rtol=0.0, atol=0.0)
        np.testing.assert_allclose(quartic_radius, radius, rtol=0.0, atol=0.0)
        np.testing.assert_allclose(quartic_abs_z, abs_z, rtol=0.0, atol=0.0)

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
                "cubic": _metrics(
                    cubic_mag,
                    radius,
                    abs_z,
                    threshold=threshold,
                    collar_start=collar_start,
                ),
                "quartic": _metrics(
                    quartic_mag,
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
    cubic_metrics = fine["cubic"]
    quartic_metrics = fine["quartic"]
    quartic_middle = middle["quartic"]

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
        summary[f"cubic_baseline_{key}_ratio"] = _ratio(cubic_metrics[key], baseline[key])
        summary[f"quartic_baseline_{key}_ratio"] = _ratio(quartic_metrics[key], baseline[key])
        summary[f"quartic_cubic_{key}_excursion_fraction"] = _excursion_fraction(
            quartic_metrics[key], baseline[key], cubic_metrics[key]
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
        summary[f"quartic_{key}_mid_to_fine_relative_change"] = _relative_change(
            quartic_metrics[key], quartic_middle[key]
        )

    return {
        "schema": "eq45_supported_phi10_quartic_late_offkeyframe_morphology_v1",
        "task_id": "CR008-EQ45-SUPPORTED-PHI10-QUARTIC-LATE-OFFKEYFRAME-MORPHOLOGY-021",
        "claim_scope": "target_free_three_resolution_late_offkeyframe_visualization_morphology_diagnostic",
        "base_supported_sha256": child.sha256,
        "cubic_candidate_sha256": cubic.sha256,
        "quartic_candidate_sha256": quartic.sha256,
        "mode": {"family": "phi", "index": [1, 0]},
        "time": scalar_time,
        "tau": float(quartic.tau(scalar_time)),
        "coefficients": {
            "static": float(quartic.midpoint_coefficient),
            "cubic": float(cubic.coefficient_at(scalar_time)),
            "quartic": float(quartic.coefficient_at(scalar_time)),
        },
        "quartic_nullspace_coefficient": float(quartic.nullspace_coefficient),
        "early_delta": float(early_delta),
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
