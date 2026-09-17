"""Three-resolution whole-domain morphology audit for the compact C2 Phi10 screen.

Agent 7's compact C2 temporal screen preserves the selected early Phi(1,0)
snapshot, returns exactly to the static supported field before the temporal
midpoint, and is deliberately more intrusive than the derivative-balanced
quartic around an early off-keyframe. Agent 2 independently found a measurable
pressure-free vorticity-equation cost for that shorter temporal window. This
module asks the remaining visualization-facing question at one frozen early
off-keyframe: does the extra temporal excursion buy a material, resolved change
in actual 3-D vorticity morphology?

Every vorticity value is reconstructed from public ``at_points(...)->[u,v,w]``
samples. A single threshold derived from the static supported inner-core peak is
applied unchanged to static, quartic and compact fields at each resolution. No
public image, PDE residual, force, pressure or hidden OpenAI parameter enters the
diagnostic.
"""
from __future__ import annotations

from collections.abc import Iterable

import numpy as np

from .constrained_eq45_supported_candidate import Eq45SupportedVelocityCandidate
from .constrained_eq45_supported_collar_vorticity import governed_supported_seed
from .constrained_eq45_supported_phi10_c2_compact_temporal_capacity import (
    DEFAULT_EARLY_DELTA,
    compact_phi10_snapshot,
    compact_window_parameters,
)
from .constrained_eq45_supported_phi10_quartic_temporal_mode import (
    Eq45SupportedPhi10QuarticDerivativeBalancedTemporalCandidate,
)
from .constrained_eq45_supported_phi10_temporal_envelope import _metrics
from .constrained_eq45_supported_vorticity_envelope import structured_vorticity_magnitude


DEFAULT_RESOLUTIONS = (49, 65, 81)
DEFAULT_TIME = 0.3125
DEFAULT_HALF_WIDTH = 2.0
DEFAULT_SUPERLEVEL_FRACTION = 0.25
DEFAULT_PLATEAU_RADIUS = 1.5
DEFAULT_PLATEAU_HALF_HEIGHT = 1.5
DEFAULT_COLLAR_START = 1.6

_TRUTH_BOUNDARY = {
    "canonical_velocity_changed": False,
    "diagnostic_velocity_changed": False,
    "visualization_candidate_only": True,
    "compact_temporal_shape_promoted": False,
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


def _excursion_fraction(value: float, baseline: float, reference: float) -> float:
    denominator = abs(float(reference) - float(baseline))
    numerator = abs(float(value) - float(baseline))
    if denominator <= 1.0e-300:
        return 0.0 if numerator <= 1.0e-300 else float("inf")
    return float(numerator / denominator)


def _phi10_coefficient(candidate: Eq45SupportedVelocityCandidate) -> float:
    basis = candidate.parent.profile_basis
    mode = (1, 0)
    if mode not in basis.mode_indices:
        raise RuntimeError("supported candidate does not contain Phi(1,0)")
    return float(basis.phi_coefficients[basis.mode_indices.index(mode)])


def audit_c2_compact_early_morphology(
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
    """Compare static, derivative-balanced quartic and compact C2 morphology."""
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
    early = float(early_delta)
    if not np.isfinite(early) or early == 0.0:
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

    compact_parameters = compact_window_parameters(early_delta=early)
    midpoint = 0.5 * (float(child.time_start) + float(child.time_end))
    half_time_width = 0.5 * (float(child.time_end) - float(child.time_start))
    tau = float((scalar_time - midpoint) / half_time_width)
    if not (-1.0 < tau < float(compact_parameters["return_tau"])):
        raise ValueError("time must lie strictly inside the active compact C2 transition")

    quartic = Eq45SupportedPhi10QuarticDerivativeBalancedTemporalCandidate(
        base=child, early_delta=early
    )
    compact_snapshot = compact_phi10_snapshot(child, scalar_time, early_delta=early)

    rows: list[dict[str, object]] = []
    for resolution in resolution_values:
        axis = np.linspace(-float(half_width), float(half_width), resolution)
        baseline_mag, radius, abs_z = structured_vorticity_magnitude(child, axis, scalar_time)
        quartic_mag, quartic_radius, quartic_abs_z = structured_vorticity_magnitude(
            quartic, axis, scalar_time
        )
        compact_mag, compact_radius, compact_abs_z = structured_vorticity_magnitude(
            compact_snapshot, axis, scalar_time
        )
        np.testing.assert_allclose(quartic_radius, radius, rtol=0.0, atol=0.0)
        np.testing.assert_allclose(quartic_abs_z, abs_z, rtol=0.0, atol=0.0)
        np.testing.assert_allclose(compact_radius, radius, rtol=0.0, atol=0.0)
        np.testing.assert_allclose(compact_abs_z, abs_z, rtol=0.0, atol=0.0)

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
                "quartic": _metrics(
                    quartic_mag,
                    radius,
                    abs_z,
                    threshold=threshold,
                    collar_start=collar_start,
                ),
                "compact": _metrics(
                    compact_mag,
                    radius,
                    abs_z,
                    threshold=threshold,
                    collar_start=collar_start,
                ),
            }
        )

    fine = rows[-1]
    middle_row = rows[-2]
    baseline = fine["baseline"]
    quartic_metrics = fine["quartic"]
    compact_metrics = fine["compact"]
    compact_middle = middle_row["compact"]

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
        summary[f"compact_baseline_{key}_ratio"] = _ratio(compact_metrics[key], baseline[key])
        summary[f"quartic_baseline_{key}_ratio"] = _ratio(quartic_metrics[key], baseline[key])
        summary[f"compact_quartic_{key}_ratio"] = _ratio(compact_metrics[key], quartic_metrics[key])
        summary[f"compact_vs_quartic_{key}_excursion_fraction"] = _excursion_fraction(
            compact_metrics[key], baseline[key], quartic_metrics[key]
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
        summary[f"compact_{key}_mid_to_fine_relative_change"] = _relative_change(
            compact_metrics[key], compact_middle[key]
        )

    return {
        "schema": "eq45_supported_phi10_c2_compact_early_morphology_v1",
        "claim_scope": "target_free_three_resolution_early_off_keyframe_visualization_morphology_diagnostic",
        "base_supported_sha256": child.sha256,
        "quartic_candidate_sha256": quartic.sha256,
        "compact_snapshot_sha256": compact_snapshot.sha256,
        "mode": {"family": "phi", "index": [1, 0]},
        "time": scalar_time,
        "tau": tau,
        "compact_return_tau": float(compact_parameters["return_tau"]),
        "compact_return_time": float(midpoint + half_time_width * float(compact_parameters["return_tau"])),
        "coefficients": {
            "static": _phi10_coefficient(child),
            "quartic": float(quartic.coefficient_at(scalar_time)),
            "compact": _phi10_coefficient(compact_snapshot),
        },
        "early_delta": early,
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
