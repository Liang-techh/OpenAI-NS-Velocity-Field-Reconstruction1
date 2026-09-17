"""Nonlinear replay of the existing supported Eq45 Phi(1,0) temporal direction.

This module consumes two already-owned ingredients instead of reimplementing them:
Agent 1's bounded support-connected affine temporal wrapper and the whole-domain
structured-curl diagnostic used by the supported vorticity-envelope lane.  It
asks one narrow representation question: does a bounded affine slope on the
*existing* ``Phi(1,0)`` mode actually suppress the resolved early-time radial
outer-vorticity branch in the nonlinear public ``[u,v,w]`` field?

The replay is target-free: it does not fit OpenAI imagery and it does not promote
any slope to the canonical candidate.  A selected row is only a visualization
trial for downstream render/fingerprint review.  Visual improvement is not PDE
validation and is not identification of an exact OpenAI field.
"""
from __future__ import annotations

from collections.abc import Iterable

import numpy as np

from .constrained_eq45_supported_candidate import Eq45SupportedVelocityCandidate
from .constrained_eq45_supported_collar_vorticity import governed_supported_seed
from .constrained_eq45_supported_temporal_mode import (
    Eq45SupportedAffineTemporalModeCandidate,
)
from .constrained_eq45_supported_vorticity_envelope import (
    structured_vorticity_magnitude,
)


DEFAULT_SLOPES = (0.8, 1.1, 1.4, 1.7, 2.0)
DEFAULT_TIMES = (0.25, 0.50, 0.75)
DEFAULT_RESOLUTION = 41
DEFAULT_HALF_WIDTH = 2.0
DEFAULT_COLLAR_START = 1.6
DEFAULT_LATE_ABSOLUTE_GUARD = 1.0e-4
DEFAULT_LATE_RELATIVE_GUARD = 10.0

_TRUTH_BOUNDARY = {
    "canonical_velocity_changed": False,
    "trial_velocity_changed": True,
    "production_slope_promoted": False,
    "new_spatial_basis_added": False,
    "openai_image_fitted": False,
    "physical_support_validated": False,
    "visualization_ready": False,
    "visual_correspondence_verified": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


def build_phi10_temporal_trial(
    base: Eq45SupportedVelocityCandidate | None = None,
    *,
    slope: float,
) -> Eq45SupportedAffineTemporalModeCandidate:
    """Build one bounded support-connected affine ``Phi(1,0)`` trial."""
    child = governed_supported_seed() if base is None else base
    if not isinstance(child, Eq45SupportedVelocityCandidate):
        raise TypeError("base must be Eq45SupportedVelocityCandidate")
    return Eq45SupportedAffineTemporalModeCandidate(
        base=child,
        family="phi",
        mode_i=1,
        mode_j=0,
        slope=float(slope),
    )


def _fingerprint(candidate, axis: np.ndarray, time: float, collar_start: float) -> dict[str, float]:
    magnitude, radius, abs_z = structured_vorticity_magnitude(candidate, axis, time)
    weight = magnitude * magnitude
    total = float(np.sum(weight))
    if not np.isfinite(total) or total <= 1.0e-300:
        raise RuntimeError("vorticity fingerprint has zero/nonfinite weight")
    radial_outer = radius > float(collar_start)
    radial_rms = float(np.sqrt(np.sum(weight * radius * radius) / total))
    axial_rms = float(np.sqrt(np.sum(weight * abs_z * abs_z) / total))
    return {
        "radial_outer_vorticity2_fraction": float(np.sum(weight[radial_outer]) / total),
        "vorticity_weighted_radial_rms": radial_rms,
        "vorticity_weighted_axial_rms": axial_rms,
        "vorticity_weighted_aspect": axial_rms / max(radial_rms, 1.0e-300),
    }


def replay_supported_phi10_temporal_slopes(
    base: Eq45SupportedVelocityCandidate | None = None,
    *,
    slopes: Iterable[float] = DEFAULT_SLOPES,
    times: Iterable[float] = DEFAULT_TIMES,
    resolution: int = DEFAULT_RESOLUTION,
    half_width: float = DEFAULT_HALF_WIDTH,
    collar_start: float = DEFAULT_COLLAR_START,
    late_absolute_guard: float = DEFAULT_LATE_ABSOLUTE_GUARD,
    late_relative_guard: float = DEFAULT_LATE_RELATIVE_GUARD,
) -> dict[str, object]:
    """Replay a small bounded slope grid and rank target-free morphology trials.

    The objective is deliberately narrow: reduce the early radial-outer
    ``|curl u|^2`` fraction while preserving the exact temporal midpoint and
    keeping the late radial-outer fraction below a preregistered guard.  This is
    a support-artifact suppression screen, not an OpenAI-image fit.
    """
    child = governed_supported_seed() if base is None else base
    if not isinstance(child, Eq45SupportedVelocityCandidate):
        raise TypeError("base must be Eq45SupportedVelocityCandidate")

    slope_values = tuple(float(value) for value in slopes)
    if not slope_values or any((not np.isfinite(value)) or value <= 0.0 for value in slope_values):
        raise ValueError("slopes must be a nonempty collection of positive finite values")
    if len(set(slope_values)) != len(slope_values):
        raise ValueError("slopes must be unique")

    time_values = tuple(float(value) for value in times)
    if len(time_values) != 3 or not all(np.isfinite(time_values)):
        raise ValueError("times must contain exactly three finite delivery times")
    if time_values != tuple(sorted(time_values)):
        raise ValueError("times must be increasing")
    midpoint = 0.5 * (child.time_start + child.time_end)
    if not np.isclose(time_values[1], midpoint, rtol=0.0, atol=1.0e-15):
        raise ValueError("middle replay time must be the candidate midpoint")
    if time_values[0] < child.time_start or time_values[-1] > child.time_end:
        raise ValueError("times must lie in the candidate delivery interval")

    resolution = int(resolution)
    if resolution < 9:
        raise ValueError("resolution must be at least 9")
    for value, name in (
        (half_width, "half_width"),
        (collar_start, "collar_start"),
        (late_absolute_guard, "late_absolute_guard"),
        (late_relative_guard, "late_relative_guard"),
    ):
        if not np.isfinite(value) or value <= 0.0:
            raise ValueError(f"{name} must be positive and finite")

    axis = np.linspace(-float(half_width), float(half_width), resolution)
    baseline_rows = [_fingerprint(child, axis, time, collar_start) for time in time_values]
    baseline_early = float(baseline_rows[0]["radial_outer_vorticity2_fraction"])
    baseline_late = float(baseline_rows[-1]["radial_outer_vorticity2_fraction"])
    late_limit = max(
        float(late_absolute_guard),
        float(late_relative_guard) * baseline_late,
    )

    trials: list[dict[str, object]] = []
    for slope in slope_values:
        candidate = build_phi10_temporal_trial(child, slope=slope)
        rows = [
            _fingerprint(candidate.snapshot(time), axis, time, collar_start)
            for time in time_values
        ]
        early = float(rows[0]["radial_outer_vorticity2_fraction"])
        late = float(rows[-1]["radial_outer_vorticity2_fraction"])
        midpoint_velocity = candidate.at_points(
            np.array([[0.37, 0.11, -0.22], [1.25, 0.0, 0.35]], dtype=float),
            midpoint,
        )
        base_midpoint_velocity = child.at_points(
            np.array([[0.37, 0.11, -0.22], [1.25, 0.0, 0.35]], dtype=float),
            midpoint,
        )
        midpoint_exact = bool(np.array_equal(midpoint_velocity, base_midpoint_velocity))
        feasible = bool(midpoint_exact and late <= late_limit)
        trials.append(
            {
                "slope": slope,
                "candidate_sha256": candidate.sha256,
                "coefficient_start": float(candidate.coefficient_at(time_values[0])),
                "coefficient_midpoint": float(candidate.coefficient_at(time_values[1])),
                "coefficient_end": float(candidate.coefficient_at(time_values[-1])),
                "fingerprints": rows,
                "early_radial_outer_fraction": early,
                "early_absolute_change": early - baseline_early,
                "early_relative_change": (early - baseline_early) / max(baseline_early, 1.0e-300),
                "late_radial_outer_fraction": late,
                "late_absolute_change": late - baseline_late,
                "midpoint_exact": midpoint_exact,
                "late_guard_limit": late_limit,
                "feasible": feasible,
            }
        )

    feasible_trials = [row for row in trials if bool(row["feasible"])]
    selected = min(
        feasible_trials,
        key=lambda row: float(row["early_radial_outer_fraction"]),
    ) if feasible_trials else None

    return {
        "schema": "eq45_supported_phi10_temporal_replay_v1",
        "claim_scope": "target_free_nonlinear_visualization_morphology_trial",
        "base_supported_sha256": child.sha256,
        "mode": {"family": "phi", "index": [1, 0]},
        "times": list(time_values),
        "resolution": resolution,
        "half_width": float(half_width),
        "collar_start": float(collar_start),
        "baseline_fingerprints": baseline_rows,
        "baseline_early_radial_outer_fraction": baseline_early,
        "baseline_late_radial_outer_fraction": baseline_late,
        "late_guard_limit": late_limit,
        "trials": trials,
        "selected_trial": selected,
        "selection_rule": "minimum_early_radial_outer_fraction_subject_to_exact_midpoint_and_late_guard",
        "truth_boundary": dict(_TRUTH_BOUNDARY),
    }
