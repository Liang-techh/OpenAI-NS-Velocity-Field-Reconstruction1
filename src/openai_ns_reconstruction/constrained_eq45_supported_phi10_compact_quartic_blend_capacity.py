"""Screen one scalar blend between existing quartic and compact Phi10 schedules.

The derivative-balanced quartic and slope-capped compact C2 schedules already
span a useful early-time morphology tradeoff on the *same* spatial mode
``Phi(1,0)``.  The compact endpoint tightens the early radial vorticity envelope
and support collar more strongly, while independent work reports a larger
pressure-free vorticity-equation cost.  Before growing the spatial basis, this
module asks the smaller representation question: is the morphology gap between
those two existing schedules smoothly controllable by one scalar temporal blend?

For blend weight ``lambda in [0,1]`` we use

    delta_blend = (1-lambda) delta_quartic + lambda delta_compact.

No image, PDE residual, force, pressure, or hidden OpenAI parameter selects
``lambda`` here.  Every morphology value is reconstructed from public
``at_points(...)->[u,v,w]`` samples.  The screen is capacity evidence only.
"""
from __future__ import annotations

from dataclasses import replace
from collections.abc import Iterable

import numpy as np

from .constrained_eq45_supported_candidate import Eq45SupportedVelocityCandidate
from .constrained_eq45_supported_collar_vorticity import governed_supported_seed
from .constrained_eq45_supported_phi10_c2_compact_temporal_capacity import (
    DEFAULT_EARLY_DELTA,
    compact_phi10_delta,
    compact_phi10_delta_dtau,
)
from .constrained_eq45_supported_phi10_cubic_temporal_capacity import DEFAULT_PROBES
from .constrained_eq45_supported_phi10_quartic_derivative_capacity import (
    quartic_phi10_delta,
    quartic_phi10_delta_dtau,
)
from .constrained_eq45_supported_phi10_temporal_envelope import _metrics
from .constrained_eq45_supported_vorticity_envelope import structured_vorticity_magnitude


_MODE = (1, 0)
DEFAULT_TIME = 0.3125
DEFAULT_BLEND_WEIGHTS = (0.0, 0.25, 0.5, 0.75, 1.0)
DEFAULT_RESOLUTION = 49
DEFAULT_HALF_WIDTH = 2.0
DEFAULT_SUPERLEVEL_FRACTION = 0.25
DEFAULT_PLATEAU_RADIUS = 1.5
DEFAULT_PLATEAU_HALF_HEIGHT = 1.5
DEFAULT_COLLAR_START = 1.6

_TRUTH_BOUNDARY = {
    "canonical_velocity_changed": False,
    "diagnostic_velocity_changed": False,
    "new_spatial_basis_added": False,
    "screen_temporal_blend_degree_added": True,
    "blend_weight_fitted": False,
    "force_or_pressure_fitted": False,
    "pde_objective_used_to_choose_blend": False,
    "public_image_fitted": False,
    "production_temporal_shape_promoted": False,
    "visualization_candidate_only": True,
    "visualization_ready": False,
    "visual_correspondence_verified": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


def _validate_weight(weight: float) -> float:
    value = float(weight)
    if not np.isfinite(value) or not (0.0 <= value <= 1.0):
        raise ValueError("blend_weight must be finite and lie in [0,1]")
    return value


def _tau(base: Eq45SupportedVelocityCandidate, time: float) -> float:
    if not isinstance(base, Eq45SupportedVelocityCandidate):
        raise TypeError("base must be Eq45SupportedVelocityCandidate")
    value = float(time)
    if not np.isfinite(value) or not (base.time_start <= value <= base.time_end):
        raise ValueError("time must be finite and lie inside the delivery interval")
    midpoint = 0.5 * (float(base.time_start) + float(base.time_end))
    half_width = 0.5 * (float(base.time_end) - float(base.time_start))
    return float((value - midpoint) / half_width)


def blended_phi10_delta(
    tau: float | np.ndarray,
    *,
    blend_weight: float,
    early_delta: float = DEFAULT_EARLY_DELTA,
):
    weight = _validate_weight(blend_weight)
    quartic = np.asarray(quartic_phi10_delta(tau, early_delta=early_delta), dtype=float)
    compact = np.asarray(compact_phi10_delta(tau, early_delta=early_delta), dtype=float)
    return (1.0 - weight) * quartic + weight * compact


def blended_phi10_delta_dtau(
    tau: float | np.ndarray,
    *,
    blend_weight: float,
    early_delta: float = DEFAULT_EARLY_DELTA,
):
    weight = _validate_weight(blend_weight)
    quartic = np.asarray(quartic_phi10_delta_dtau(tau, early_delta=early_delta), dtype=float)
    compact = np.asarray(compact_phi10_delta_dtau(tau, early_delta=early_delta), dtype=float)
    return (1.0 - weight) * quartic + weight * compact


def _phi10_index_and_midpoint(base: Eq45SupportedVelocityCandidate) -> tuple[int, float]:
    basis = base.parent.profile_basis
    if _MODE not in basis.mode_indices:
        raise RuntimeError("supported candidate does not contain Phi(1,0)")
    index = basis.mode_indices.index(_MODE)
    return index, float(basis.phi_coefficients[index])


def blended_phi10_snapshot(
    base: Eq45SupportedVelocityCandidate,
    time: float,
    *,
    blend_weight: float,
    early_delta: float = DEFAULT_EARLY_DELTA,
) -> Eq45SupportedVelocityCandidate:
    """Return one support-connected snapshot on the quartic/compact blend path."""
    weight = _validate_weight(blend_weight)
    tau = _tau(base, time)
    index, midpoint = _phi10_index_and_midpoint(base)
    delta = float(
        blended_phi10_delta(tau, blend_weight=weight, early_delta=early_delta)
    )
    coefficient = midpoint + delta
    limit = float(base.parent.profile_basis.coefficient_limit)
    # Pointwise convexity also implies the whole blended schedule stays inside
    # the common endpoint-family bound whenever both endpoints do.
    if not np.isfinite(coefficient) or not (-limit <= coefficient <= limit):
        raise ValueError("blended Phi(1,0) coefficient exceeds the existing bound")
    phi = list(base.parent.profile_basis.phi_coefficients)
    phi[index] = coefficient
    basis = replace(base.parent.profile_basis, phi_coefficients=tuple(phi))
    parent = replace(base.parent, profile_basis=basis)
    return Eq45SupportedVelocityCandidate(parent=parent, taper=base.taper)


def _rms(values: np.ndarray) -> float:
    array = np.asarray(values, dtype=float)
    return float(np.sqrt(np.mean(array * array)))


def _progress(value: float, quartic: float, compact: float) -> float:
    denominator = float(compact) - float(quartic)
    if abs(denominator) <= np.finfo(float).tiny:
        return 0.0
    return float((float(value) - float(quartic)) / denominator)


def audit_supported_phi10_compact_quartic_blend_capacity(
    base: Eq45SupportedVelocityCandidate | None = None,
    *,
    time: float = DEFAULT_TIME,
    blend_weights: Iterable[float] = DEFAULT_BLEND_WEIGHTS,
    early_delta: float = DEFAULT_EARLY_DELTA,
    resolution: int = DEFAULT_RESOLUTION,
    probes: np.ndarray = DEFAULT_PROBES,
    half_width: float = DEFAULT_HALF_WIDTH,
    superlevel_fraction: float = DEFAULT_SUPERLEVEL_FRACTION,
    plateau_radius: float = DEFAULT_PLATEAU_RADIUS,
    plateau_half_height: float = DEFAULT_PLATEAU_HALF_HEIGHT,
    collar_start: float = DEFAULT_COLLAR_START,
) -> dict[str, object]:
    """Quantify one-degree morphology control between quartic and compact endpoints."""
    child = governed_supported_seed() if base is None else base
    if not isinstance(child, Eq45SupportedVelocityCandidate):
        raise TypeError("base must be Eq45SupportedVelocityCandidate")
    scalar_time = float(time)
    tau = _tau(child, scalar_time)
    if not (-1.0 < tau < 0.0):
        raise ValueError("time must lie in the early transition tau in (-1,0)")

    weights = tuple(_validate_weight(value) for value in blend_weights)
    if len(weights) < 3:
        raise ValueError("blend_weights must contain at least three values")
    if weights[0] != 0.0 or weights[-1] != 1.0:
        raise ValueError("blend_weights must include quartic=0 and compact=1 endpoints")
    if not all(b > a for a, b in zip(weights, weights[1:])):
        raise ValueError("blend_weights must be strictly increasing")

    grid_n = int(resolution)
    if grid_n < 17 or grid_n % 2 == 0:
        raise ValueError("resolution must be an odd integer >= 17")
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

    point_array = np.asarray(probes, dtype=float)
    if point_array.ndim != 2 or point_array.shape[1] != 3 or point_array.shape[0] < 3:
        raise ValueError("probes must have shape (n,3) with n>=3")
    if not np.all(np.isfinite(point_array)):
        raise ValueError("probes must be finite")

    axis = np.linspace(-float(half_width), float(half_width), grid_n)
    baseline_mag, radius, abs_z = structured_vorticity_magnitude(child, axis, scalar_time)
    plateau = (radius <= float(plateau_radius)) & (abs_z <= float(plateau_half_height))
    core_peak = float(np.max(baseline_mag[plateau]))
    threshold = float(superlevel_fraction * core_peak)
    baseline_metrics = _metrics(
        baseline_mag,
        radius,
        abs_z,
        threshold=threshold,
        collar_start=collar_start,
    )

    quartic_snapshot = blended_phi10_snapshot(
        child, scalar_time, blend_weight=0.0, early_delta=early_delta
    )
    compact_snapshot = blended_phi10_snapshot(
        child, scalar_time, blend_weight=1.0, early_delta=early_delta
    )
    quartic_probe = quartic_snapshot.at_points(point_array, scalar_time)
    compact_probe = compact_snapshot.at_points(point_array, scalar_time)

    rows: list[dict[str, object]] = []
    for weight in weights:
        snapshot = blended_phi10_snapshot(
            child, scalar_time, blend_weight=weight, early_delta=early_delta
        )
        magnitude, row_radius, row_abs_z = structured_vorticity_magnitude(
            snapshot, axis, scalar_time
        )
        np.testing.assert_allclose(row_radius, radius, rtol=0.0, atol=0.0)
        np.testing.assert_allclose(row_abs_z, abs_z, rtol=0.0, atol=0.0)
        metrics = _metrics(
            magnitude,
            radius,
            abs_z,
            threshold=threshold,
            collar_start=collar_start,
        )
        public_velocity = snapshot.at_points(point_array, scalar_time)
        affine_reference = (1.0 - weight) * quartic_probe + weight * compact_probe
        difference = public_velocity - affine_reference
        rows.append(
            {
                "blend_weight": weight,
                "phi10_coefficient": float(
                    _phi10_index_and_midpoint(snapshot)[1]
                ),
                "metrics": metrics,
                "public_velocity_affine_linearity_rms": _rms(difference),
                "public_velocity_affine_linearity_max_abs": float(np.max(np.abs(difference))),
            }
        )

    endpoint_quartic = rows[0]["metrics"]
    endpoint_compact = rows[-1]["metrics"]
    progress_keys = (
        "vorticity2_weighted_radial_rms",
        "vorticity2_weighted_axial_rms",
        "vorticity2_weighted_aspect",
        "whole_grid_collar_vorticity2_fraction",
    )
    max_progress_nonlinearity: dict[str, float] = {}
    for key in progress_keys:
        deviations = []
        for row in rows:
            progress = _progress(
                row["metrics"][key], endpoint_quartic[key], endpoint_compact[key]
            )
            row[f"{key}_endpoint_progress"] = progress
            deviations.append(abs(progress - float(row["blend_weight"])))
        max_progress_nonlinearity[key] = float(max(deviations))

    dense_tau = np.linspace(-1.0, 1.0, 4001)
    quartic_peak = float(
        np.max(np.abs(quartic_phi10_delta_dtau(dense_tau, early_delta=early_delta)))
    )
    compact_peak = float(
        np.max(np.abs(compact_phi10_delta_dtau(dense_tau, early_delta=early_delta)))
    )
    derivative_rows = []
    for weight in weights:
        derivative = np.asarray(
            blended_phi10_delta_dtau(
                dense_tau, blend_weight=weight, early_delta=early_delta
            ),
            dtype=float,
        )
        derivative_rows.append(
            {
                "blend_weight": weight,
                "peak_abs_d_delta_dtau": float(np.max(np.abs(derivative))),
                "integral_d_delta_dtau_squared": float(np.trapz(derivative * derivative, dense_tau)),
            }
        )

    return {
        "schema": "eq45_supported_phi10_compact_quartic_blend_capacity_v1",
        "task_id": "CR003-EQ45-SUPPORTED-PHI10-COMPACT-QUARTIC-BLEND-CAPACITY-031",
        "claim_scope": "target_free_single_temporal_blend_degree_capacity_screen",
        "base_supported_sha256": child.sha256,
        "mode": {"family": "phi", "index": [1, 0]},
        "time": scalar_time,
        "tau": tau,
        "early_delta": float(early_delta),
        "blend_weights": list(weights),
        "resolution": grid_n,
        "shared_threshold": threshold,
        "baseline_metrics": baseline_metrics,
        "rows": rows,
        "max_endpoint_progress_minus_weight": max_progress_nonlinearity,
        "temporal_derivative": {
            "quartic_peak_abs_d_delta_dtau": quartic_peak,
            "compact_peak_abs_d_delta_dtau": compact_peak,
            "rows": derivative_rows,
            "peak_cap_preserved_by_convexity": True,
        },
        "representation": {
            "spatial_basis_increment": 0,
            "temporal_scalar_degree_increment_if_materialized": 1,
            "blend_weight_selected": False,
            "endpoint_families": ["derivative_balanced_quartic", "slope_capped_c2_compact"],
        },
        "truth_boundary": dict(_TRUTH_BOUNDARY),
    }
