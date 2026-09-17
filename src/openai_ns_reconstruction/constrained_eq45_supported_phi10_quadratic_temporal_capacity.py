"""Capacity screen for an early-localized temporal shape on supported Eq45 Phi(1,0).

The selected affine Phi(1,0) trial removes a resolved early support-induced
radial vorticity branch, but a later whole-domain audit also finds nontrivial
central morphology change at the late endpoint.  This module asks the smallest
representation question created by those two results: can one additional
*temporal* basis direction keep the successful early endpoint while returning
to the static supported field at both the midpoint and the late endpoint?

No spatial basis term is added.  The screen uses the existing Phi(1,0) mode and
compares the selected affine schedule against the unique quadratic schedule

    delta c(tau) = 0.5 * d_early * tau * (tau - 1),

where tau=-1,0,+1 at the delivery start/mid/end.  Hence the quadratic trial has
(delta_start, delta_mid, delta_end)=(d_early,0,0).  With d_early=-1.4 it keeps
exactly the early coefficient used by the selected affine slope-1.4 trial while
removing that trial's +1.4 late excursion.

All morphology is rebuilt from public support-connected ``[u,v,w]`` samples.
This is target-free representation-capacity evidence.  It does not fit an
OpenAI image, promote a production candidate, validate the PDE, or identify an
exact OpenAI field.
"""
from __future__ import annotations

from dataclasses import replace
from collections.abc import Iterable

import numpy as np

from .constrained_eq45_supported_candidate import Eq45SupportedVelocityCandidate
from .constrained_eq45_supported_collar_vorticity import governed_supported_seed
from .constrained_eq45_supported_phi10_temporal_envelope import _metrics
from .constrained_eq45_supported_phi10_temporal_replay import build_phi10_temporal_trial
from .constrained_eq45_supported_vorticity_envelope import structured_vorticity_magnitude


DEFAULT_EARLY_DELTA = -1.4
DEFAULT_AFFINE_SLOPE = 1.4
DEFAULT_RESOLUTIONS = (49, 65, 81)
DEFAULT_REFERENCE_TIMES = (0.25, 0.50, 0.75)
DEFAULT_INTERMEDIATE_TIME = 0.625
DEFAULT_HALF_WIDTH = 2.0
DEFAULT_SUPERLEVEL_FRACTION = 0.25
DEFAULT_PLATEAU_RADIUS = 1.5
DEFAULT_PLATEAU_HALF_HEIGHT = 1.5
DEFAULT_COLLAR_START = 1.6

_TRUTH_BOUNDARY = {
    "canonical_velocity_changed": False,
    "screen_trial_velocity_changed": True,
    "production_temporal_shape_promoted": False,
    "new_spatial_basis_added": False,
    "temporal_basis_dimension_increment": 1,
    "openai_image_fitted": False,
    "physical_support_validated": False,
    "visualization_ready": False,
    "visual_correspondence_verified": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


def _tau(base: Eq45SupportedVelocityCandidate, time: float) -> float:
    value = float(time)
    if not np.isfinite(value):
        raise ValueError("time must be finite")
    if value < base.time_start or value > base.time_end:
        raise ValueError("time must lie inside the candidate delivery interval")
    midpoint = 0.5 * (base.time_start + base.time_end)
    half_width = 0.5 * (base.time_end - base.time_start)
    return (value - midpoint) / half_width


def quadratic_phi10_delta(tau: float | np.ndarray, *, early_delta: float = DEFAULT_EARLY_DELTA):
    """Return the early-localized quadratic coefficient displacement."""
    early = float(early_delta)
    if not np.isfinite(early) or early == 0.0:
        raise ValueError("early_delta must be finite and nonzero")
    tau_arr = np.asarray(tau, dtype=float)
    if not np.all(np.isfinite(tau_arr)) or np.any((tau_arr < -1.0) | (tau_arr > 1.0)):
        raise ValueError("tau must be finite and lie in [-1,1]")
    return 0.5 * early * tau_arr * (tau_arr - 1.0)


def _phi10_midpoint_coefficient(base: Eq45SupportedVelocityCandidate) -> tuple[int, float]:
    basis = base.parent.profile_basis
    mode = (1, 0)
    if mode not in basis.mode_indices:
        raise RuntimeError("base profile basis does not contain Phi(1,0)")
    index = basis.mode_indices.index(mode)
    return index, float(basis.phi_coefficients[index])


def quadratic_phi10_coefficient(
    base: Eq45SupportedVelocityCandidate,
    time: float,
    *,
    early_delta: float = DEFAULT_EARLY_DELTA,
) -> float:
    """Return the bounded Phi(1,0) coefficient for the quadratic screen."""
    if not isinstance(base, Eq45SupportedVelocityCandidate):
        raise TypeError("base must be Eq45SupportedVelocityCandidate")
    _, midpoint = _phi10_midpoint_coefficient(base)
    return midpoint + float(quadratic_phi10_delta(_tau(base, time), early_delta=early_delta))


def _schedule_extrema(base: Eq45SupportedVelocityCandidate, early_delta: float) -> dict[str, float]:
    _, midpoint = _phi10_midpoint_coefficient(base)
    tau_values = np.array([-1.0, 0.0, 0.5, 1.0], dtype=float)
    coefficients = midpoint + quadratic_phi10_delta(tau_values, early_delta=early_delta)
    limit = float(base.parent.profile_basis.coefficient_limit)
    lo = float(np.min(coefficients))
    hi = float(np.max(coefficients))
    if lo < -limit or hi > limit:
        raise ValueError(
            "quadratic Phi(1,0) schedule exceeds the existing profile coefficient bound"
        )
    return {
        "minimum": lo,
        "maximum": hi,
        "coefficient_limit": limit,
        "vertex_tau": 0.5,
        "vertex_coefficient": float(
            midpoint + quadratic_phi10_delta(0.5, early_delta=early_delta)
        ),
    }


def quadratic_phi10_snapshot(
    base: Eq45SupportedVelocityCandidate,
    time: float,
    *,
    early_delta: float = DEFAULT_EARLY_DELTA,
) -> Eq45SupportedVelocityCandidate:
    """Materialize one static supported snapshot of the quadratic screen."""
    if not isinstance(base, Eq45SupportedVelocityCandidate):
        raise TypeError("base must be Eq45SupportedVelocityCandidate")
    _schedule_extrema(base, float(early_delta))
    index, _ = _phi10_midpoint_coefficient(base)
    coefficient = quadratic_phi10_coefficient(base, float(time), early_delta=early_delta)
    phi = list(base.parent.profile_basis.phi_coefficients)
    phi[index] = coefficient
    basis = replace(base.parent.profile_basis, phi_coefficients=tuple(phi))
    parent = replace(base.parent, profile_basis=basis)
    return Eq45SupportedVelocityCandidate(parent=parent, taper=base.taper)


def _velocity_delta_rms(
    candidate: Eq45SupportedVelocityCandidate,
    baseline: Eq45SupportedVelocityCandidate,
    points: np.ndarray,
    time: float,
) -> float:
    delta = candidate.at_points(points, time) - baseline.at_points(points, time)
    return float(np.sqrt(np.mean(delta * delta)))


def _time_design() -> dict[str, object]:
    tau = np.array([-1.0, -0.5, 0.0, 0.5, 1.0], dtype=float)
    design = np.column_stack((tau, tau * tau))
    singular_values = np.linalg.svd(design, compute_uv=False)
    return {
        "sample_tau": tau.tolist(),
        "basis": ["tau", "tau_squared"],
        "rank": int(np.linalg.matrix_rank(design)),
        "singular_values": singular_values.tolist(),
        "condition_number": float(singular_values[0] / singular_values[-1]),
        "affine_temporal_basis_dimension": 1,
        "screen_temporal_basis_dimension": 2,
        "basis_dimension_increment": 1,
    }


def audit_supported_phi10_early_localized_temporal_capacity(
    base: Eq45SupportedVelocityCandidate | None = None,
    *,
    early_delta: float = DEFAULT_EARLY_DELTA,
    affine_slope: float = DEFAULT_AFFINE_SLOPE,
    resolutions: Iterable[int] = DEFAULT_RESOLUTIONS,
    reference_times: Iterable[float] = DEFAULT_REFERENCE_TIMES,
    intermediate_time: float = DEFAULT_INTERMEDIATE_TIME,
    half_width: float = DEFAULT_HALF_WIDTH,
    superlevel_fraction: float = DEFAULT_SUPERLEVEL_FRACTION,
    plateau_radius: float = DEFAULT_PLATEAU_RADIUS,
    plateau_half_height: float = DEFAULT_PLATEAU_HALF_HEIGHT,
    collar_start: float = DEFAULT_COLLAR_START,
) -> dict[str, object]:
    """Compare affine Phi10 against one early-localized quadratic time shape."""
    child = governed_supported_seed() if base is None else base
    if not isinstance(child, Eq45SupportedVelocityCandidate):
        raise TypeError("base must be Eq45SupportedVelocityCandidate")

    early_delta = float(early_delta)
    affine_slope = float(affine_slope)
    if not np.isfinite(affine_slope) or affine_slope == 0.0:
        raise ValueError("affine_slope must be finite and nonzero")
    if not np.isclose(early_delta, -affine_slope, rtol=0.0, atol=1.0e-15):
        raise ValueError("early_delta must equal -affine_slope for an endpoint-matched screen")
    schedule_extrema = _schedule_extrema(child, early_delta)

    resolution_values = tuple(int(value) for value in resolutions)
    if len(resolution_values) < 3 or any(value < 9 for value in resolution_values):
        raise ValueError("resolutions must contain at least three levels >= 9")
    if not all(b > a for a, b in zip(resolution_values, resolution_values[1:])):
        raise ValueError("resolutions must be strictly increasing")

    time_values = tuple(float(value) for value in reference_times)
    if len(time_values) != 3 or not all(np.isfinite(value) for value in time_values):
        raise ValueError("reference_times must contain exactly three finite values")
    if time_values != tuple(sorted(time_values)):
        raise ValueError("reference_times must be increasing")
    midpoint = 0.5 * (child.time_start + child.time_end)
    if not np.isclose(time_values[0], child.time_start, rtol=0.0, atol=1.0e-15):
        raise ValueError("first reference time must equal the delivery start")
    if not np.isclose(time_values[1], midpoint, rtol=0.0, atol=1.0e-15):
        raise ValueError("middle reference time must equal the delivery midpoint")
    if not np.isclose(time_values[2], child.time_end, rtol=0.0, atol=1.0e-15):
        raise ValueError("last reference time must equal the delivery end")
    intermediate_time = float(intermediate_time)
    if not (midpoint < intermediate_time < child.time_end):
        raise ValueError("intermediate_time must lie strictly between midpoint and end")

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

    affine = build_phi10_temporal_trial(child, slope=affine_slope)
    rows: list[dict[str, object]] = []
    for resolution in resolution_values:
        axis = np.linspace(-float(half_width), float(half_width), resolution)
        for time in time_values:
            baseline_mag, radius, abs_z = structured_vorticity_magnitude(child, axis, time)
            affine_mag, affine_radius, affine_abs_z = structured_vorticity_magnitude(
                affine.snapshot(time), axis, time
            )
            quadratic = quadratic_phi10_snapshot(child, time, early_delta=early_delta)
            quadratic_mag, quadratic_radius, quadratic_abs_z = structured_vorticity_magnitude(
                quadratic, axis, time
            )
            plateau = (radius <= float(plateau_radius)) & (abs_z <= float(plateau_half_height))
            threshold = float(superlevel_fraction * np.max(baseline_mag[plateau]))
            rows.append(
                {
                    "resolution": resolution,
                    "dx": float(axis[1] - axis[0]),
                    "time": time,
                    "shared_threshold": threshold,
                    "baseline": _metrics(
                        baseline_mag, radius, abs_z, threshold=threshold, collar_start=collar_start
                    ),
                    "affine": _metrics(
                        affine_mag,
                        affine_radius,
                        affine_abs_z,
                        threshold=threshold,
                        collar_start=collar_start,
                    ),
                    "quadratic": _metrics(
                        quadratic_mag,
                        quadratic_radius,
                        quadratic_abs_z,
                        threshold=threshold,
                        collar_start=collar_start,
                    ),
                }
            )

    fine_resolution = resolution_values[-1]
    fine_rows = {
        float(row["time"]): row for row in rows if int(row["resolution"]) == fine_resolution
    }
    summaries: list[dict[str, float | bool]] = []
    for time in time_values:
        row = fine_rows[time]
        baseline = row["baseline"]
        affine_metrics = row["affine"]
        quadratic_metrics = row["quadratic"]
        summaries.append(
            {
                "time": time,
                "quadratic_baseline_radial_q99_ratio": float(
                    quadratic_metrics["radial_q99"] / max(baseline["radial_q99"], 1.0e-300)
                ),
                "quadratic_baseline_axial_q99_ratio": float(
                    quadratic_metrics["axial_q99"] / max(baseline["axial_q99"], 1.0e-300)
                ),
                "quadratic_baseline_weighted_aspect_ratio": float(
                    quadratic_metrics["vorticity2_weighted_aspect"]
                    / max(baseline["vorticity2_weighted_aspect"], 1.0e-300)
                ),
                "quadratic_baseline_radial_outer_vorticity2_ratio": float(
                    quadratic_metrics["radial_outer_vorticity2_fraction"]
                    / max(baseline["radial_outer_vorticity2_fraction"], 1.0e-300)
                ),
                "quadratic_affine_radial_q99_ratio": float(
                    quadratic_metrics["radial_q99"] / max(affine_metrics["radial_q99"], 1.0e-300)
                ),
                "quadratic_affine_axial_q99_ratio": float(
                    quadratic_metrics["axial_q99"] / max(affine_metrics["axial_q99"], 1.0e-300)
                ),
                "quadratic_affine_radial_outer_vorticity2_ratio": float(
                    quadratic_metrics["radial_outer_vorticity2_fraction"]
                    / max(affine_metrics["radial_outer_vorticity2_fraction"], 1.0e-300)
                ),
            }
        )

    probe_points = np.array(
        [
            [0.37, 0.11, -0.22],
            [0.95, -0.35, 0.48],
            [1.45, 0.15, -0.62],
            [1.72, 0.0, 0.30],
        ],
        dtype=float,
    )
    endpoint_equivalence = {}
    for time, expected in (
        (time_values[0], "affine"),
        (time_values[1], "baseline"),
        (time_values[2], "baseline"),
    ):
        quadratic_velocity = quadratic_phi10_snapshot(
            child, time, early_delta=early_delta
        ).at_points(probe_points, time)
        expected_velocity = (
            affine.snapshot(time).at_points(probe_points, time)
            if expected == "affine"
            else child.at_points(probe_points, time)
        )
        endpoint_equivalence[str(time)] = {
            "expected": expected,
            "exact_public_velocity_match": bool(np.array_equal(quadratic_velocity, expected_velocity)),
            "max_abs_difference": float(np.max(np.abs(quadratic_velocity - expected_velocity))),
        }

    quadratic_intermediate = quadratic_phi10_snapshot(
        child, intermediate_time, early_delta=early_delta
    )
    affine_intermediate = affine.snapshot(intermediate_time)
    quadratic_delta_rms = _velocity_delta_rms(
        quadratic_intermediate, child, probe_points, intermediate_time
    )
    affine_delta_rms = _velocity_delta_rms(
        affine_intermediate, child, probe_points, intermediate_time
    )
    velocity_delta_ratio = quadratic_delta_rms / max(affine_delta_rms, 1.0e-300)

    intermediate_axis = np.linspace(-float(half_width), float(half_width), fine_resolution)
    baseline_mag, radius, abs_z = structured_vorticity_magnitude(
        child, intermediate_axis, intermediate_time
    )
    affine_mag, affine_radius, affine_abs_z = structured_vorticity_magnitude(
        affine_intermediate, intermediate_axis, intermediate_time
    )
    quadratic_mag, quadratic_radius, quadratic_abs_z = structured_vorticity_magnitude(
        quadratic_intermediate, intermediate_axis, intermediate_time
    )
    plateau = (radius <= float(plateau_radius)) & (abs_z <= float(plateau_half_height))
    intermediate_threshold = float(superlevel_fraction * np.max(baseline_mag[plateau]))
    intermediate = {
        "time": intermediate_time,
        "tau": _tau(child, intermediate_time),
        "baseline_coefficient": _phi10_midpoint_coefficient(child)[1],
        "affine_coefficient": float(affine.coefficient_at(intermediate_time)),
        "quadratic_coefficient": quadratic_phi10_coefficient(
            child, intermediate_time, early_delta=early_delta
        ),
        "quadratic_to_affine_public_velocity_delta_rms_ratio": float(velocity_delta_ratio),
        "baseline": _metrics(
            baseline_mag, radius, abs_z, threshold=intermediate_threshold, collar_start=collar_start
        ),
        "affine": _metrics(
            affine_mag,
            affine_radius,
            affine_abs_z,
            threshold=intermediate_threshold,
            collar_start=collar_start,
        ),
        "quadratic": _metrics(
            quadratic_mag,
            quadratic_radius,
            quadratic_abs_z,
            threshold=intermediate_threshold,
            collar_start=collar_start,
        ),
    }

    tau_coeff = -0.5 * early_delta
    tau2_coeff = 0.5 * early_delta
    return {
        "schema": "eq45_supported_phi10_early_localized_temporal_capacity_v1",
        "claim_scope": "target_free_temporal_basis_capacity_screen",
        "base_supported_sha256": child.sha256,
        "mode": {"family": "phi", "index": [1, 0]},
        "reference_times": list(time_values),
        "resolutions": list(resolution_values),
        "early_delta": early_delta,
        "affine_slope": affine_slope,
        "quadratic_schedule": {
            "formula": "delta=0.5*early_delta*tau*(tau-1)",
            "tau_coefficient": tau_coeff,
            "tau_squared_coefficient": tau2_coeff,
            "start_delta": float(quadratic_phi10_delta(-1.0, early_delta=early_delta)),
            "midpoint_delta": float(quadratic_phi10_delta(0.0, early_delta=early_delta)),
            "end_delta": float(quadratic_phi10_delta(1.0, early_delta=early_delta)),
            "extrema": schedule_extrema,
        },
        "time_design": _time_design(),
        "rows": rows,
        "fine_reference_summaries": summaries,
        "endpoint_equivalence": endpoint_equivalence,
        "intermediate_collateral": intermediate,
        "pde_validation_rerun": False,
        "pde_note": "PDE tradeoff is owned by sibling CR005/CR009 lanes; this screen changes no production candidate and does not inherit a PDE pass.",
        "truth_boundary": dict(_TRUTH_BOUNDARY),
    }
