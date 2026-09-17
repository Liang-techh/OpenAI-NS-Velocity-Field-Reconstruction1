"""Local capacity audit for the supported Eq45 early radial vorticity branch.

The whole-domain envelope audit shows that the physical support transform can
produce a larger early-time radial high-vorticity envelope than the untapered
parent.  This module asks a narrower representation question: can the two
*already present* bounded poloidal modes ``Phi(1,0)`` and ``Phi(1,2)`` move a
smooth radial-collar vorticity fingerprint mainly at the early delivery time,
or is a time-shaped coefficient direction more efficient?

No OpenAI image is fitted here and no saved candidate is changed.  All
vorticity is reconstructed from public ``at_points(...)->[u,v,w]`` samples via
the independent structured-grid curl helper used by the whole-domain audit.
"""
from __future__ import annotations

from dataclasses import replace

import numpy as np

from .constrained_eq45_candidate import Eq45VelocityCandidate
from .constrained_eq45_supported_candidate import Eq45SupportedVelocityCandidate
from .constrained_eq45_supported_vorticity_envelope import (
    structured_vorticity_magnitude,
)


MODES = ((1, 0), (1, 2))
DEFAULT_TIMES = (0.25, 0.50, 0.75)
DEFAULT_STEPS = (0.04, 0.02)
DEFAULT_RESOLUTION = 41
DEFAULT_HALF_WIDTH = 2.0
DEFAULT_COLLAR_START = 1.6

_TRUTH_BOUNDARY = {
    "velocity_changed": False,
    "production_coefficients_changed": False,
    "new_basis_added": False,
    "temporal_mode_implemented": False,
    "forcing_or_pressure_refit": False,
    "visualization_ready": False,
    "visual_correspondence_verified": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


def _replace_phi_mode(
    candidate: Eq45VelocityCandidate,
    mode: tuple[int, int],
    delta: float,
) -> Eq45VelocityCandidate:
    basis = candidate.profile_basis
    try:
        index = basis.mode_indices.index(mode)
    except ValueError as exc:
        raise ValueError(f"candidate does not contain Phi mode {mode}") from exc
    coefficients = list(basis.phi_coefficients)
    coefficients[index] += float(delta)
    new_basis = replace(basis, phi_coefficients=tuple(coefficients))
    return replace(candidate, profile_basis=new_basis)


def _perturbed_child(
    child: Eq45SupportedVelocityCandidate,
    deltas: dict[tuple[int, int], float],
) -> Eq45SupportedVelocityCandidate:
    parent = child.parent
    for mode, delta in deltas.items():
        parent = _replace_phi_mode(parent, mode, delta)
    return Eq45SupportedVelocityCandidate(parent=parent, taper=child.taper)


def _radial_fingerprint(
    candidate: Eq45SupportedVelocityCandidate,
    axis: np.ndarray,
    time: float,
    collar_start: float,
) -> dict[str, float]:
    magnitude, radius, abs_z = structured_vorticity_magnitude(candidate, axis, time)
    weight = magnitude * magnitude
    total = float(np.sum(weight))
    if not np.isfinite(total) or total <= 1e-300:
        raise RuntimeError("vorticity fingerprint has zero/nonfinite weight")

    radial_outer = radius > float(collar_start)
    radial_outer_fraction = float(np.sum(weight[radial_outer]) / total)
    radial_rms = float(np.sqrt(np.sum(weight * radius * radius) / total))
    axial_rms = float(np.sqrt(np.sum(weight * abs_z * abs_z) / total))
    return {
        "radial_outer_vorticity2_fraction": radial_outer_fraction,
        "vorticity_weighted_radial_rms": radial_rms,
        "vorticity_weighted_axial_rms": axial_rms,
        "vorticity_weighted_aspect": axial_rms / max(radial_rms, 1e-300),
    }


def _feature_vector(
    candidate: Eq45SupportedVelocityCandidate,
    axis: np.ndarray,
    times: tuple[float, ...],
    collar_start: float,
) -> tuple[np.ndarray, list[dict[str, float]]]:
    rows = [
        _radial_fingerprint(candidate, axis, time, collar_start)
        for time in times
    ]
    return (
        np.asarray([row["radial_outer_vorticity2_fraction"] for row in rows], dtype=float),
        rows,
    )


def _response_matrix(
    child: Eq45SupportedVelocityCandidate,
    axis: np.ndarray,
    times: tuple[float, ...],
    collar_start: float,
    step: float,
) -> np.ndarray:
    columns = []
    for mode in MODES:
        plus = _perturbed_child(child, {mode: step})
        minus = _perturbed_child(child, {mode: -step})
        fp, _ = _feature_vector(plus, axis, times, collar_start)
        fm, _ = _feature_vector(minus, axis, times, collar_start)
        columns.append((fp - fm) / (2.0 * step))
    return np.column_stack(columns)


def _affine_temporal_screen(
    child: Eq45SupportedVelocityCandidate,
    times: tuple[float, ...],
    fine_matrix: np.ndarray,
    baseline: np.ndarray,
    local_budget: float,
    static_predicted_delta: np.ndarray,
) -> dict[str, object]:
    """Screen one affine coefficient slope without implementing a new candidate.

    For the standard affine coordinate ``tau=-1,0,+1`` at delivery start,
    midpoint and end, the local response of mode ``i`` is simply
    ``tau * d(feature)/d(coefficient_i)``.  This consumes the measured static
    Jacobian rather than duplicating Agent 1's production temporal wrapper.
    """
    midpoint = 0.5 * (child.time_start + child.time_end)
    tau = 2.0 * (np.asarray(times, dtype=float) - midpoint) / (
        child.time_end - child.time_start
    )
    temporal_matrix = fine_matrix * tau[:, None]
    static_early_suppression = max(-float(static_predicted_delta[0]), np.finfo(float).tiny)

    rows: dict[str, dict[str, float | list[float]]] = {}
    basis = child.parent.profile_basis
    best_mode = None
    best_same_budget_gain = -np.inf
    for index, mode in enumerate(MODES):
        signature = temporal_matrix[:, index]
        slope_sign = -1.0 if signature[0] > 0.0 else 1.0
        signed_signature = slope_sign * signature
        early_suppression = max(-float(signed_signature[0]), 0.0)
        late_leak = float(np.linalg.norm(signed_signature[1:]))
        selectivity = early_suppression / max(late_leak, np.finfo(float).tiny)

        coefficient_index = basis.mode_indices.index(mode)
        coefficient = float(basis.phi_coefficients[coefficient_index])
        slope_headroom = float(basis.coefficient_limit - abs(coefficient))
        slope_to_zero = (
            float(baseline[0] / early_suppression)
            if early_suppression > np.finfo(float).tiny
            else float("inf")
        )
        same_budget_gain = float(local_budget * early_suppression)
        efficiency_vs_static = same_budget_gain / static_early_suppression
        rows[str(mode)] = {
            "slope_sign_to_suppress_early": slope_sign,
            "signed_response_per_unit_abs_slope": signed_signature.tolist(),
            "early_suppression_per_unit_abs_slope": early_suppression,
            "late_leak_l2_per_unit_abs_slope": late_leak,
            "early_to_late_selectivity": selectivity,
            "max_abs_slope_under_existing_bound": slope_headroom,
            "linearized_abs_slope_to_zero_early_feature": slope_to_zero,
            "linearized_zero_early_within_bound": bool(slope_to_zero <= slope_headroom),
            "same_local_budget_early_suppression": same_budget_gain,
            "same_local_budget_efficiency_vs_static_cancellation": efficiency_vs_static,
        }
        if same_budget_gain > best_same_budget_gain:
            best_mode = mode
            best_same_budget_gain = same_budget_gain

    return {
        "tau": tau.tolist(),
        "local_abs_slope_budget": float(local_budget),
        "modes": rows,
        "recommended_existing_mode_for_temporal_followup": str(best_mode),
        "production_temporal_candidate_implemented_here": False,
    }


def audit_supported_early_radial_capacity(
    child: Eq45SupportedVelocityCandidate,
    *,
    times: tuple[float, ...] = DEFAULT_TIMES,
    coefficient_steps: tuple[float, float] = DEFAULT_STEPS,
    resolution: int = DEFAULT_RESOLUTION,
    half_width: float = DEFAULT_HALF_WIDTH,
    collar_start: float = DEFAULT_COLLAR_START,
) -> dict[str, object]:
    """Audit whether existing static Phi modes can isolate an early radial correction.

    The smooth morphology feature is the fraction of whole-grid ``|curl u|^2``
    lying at physical radius ``r > collar_start``.  The region itself is fixed,
    avoiding coefficient-dependent voxel jumps of a hard vorticity superlevel
    extent.  We compare two ways to obtain an early-localized change: static
    cancellation between the existing Phi modes, and the local signature that
    one affine-in-time coefficient slope would provide.  The latter is only a
    Jacobian screen; no temporal production field is constructed here.
    """
    if not isinstance(child, Eq45SupportedVelocityCandidate):
        raise TypeError("child must be Eq45SupportedVelocityCandidate")
    times = tuple(float(value) for value in times)
    if len(times) != 3 or not all(np.isfinite(times)):
        raise ValueError("times must contain exactly three finite delivery times")
    if not (times[0] < times[1] < times[2]):
        raise ValueError("times must be strictly increasing")
    if times[0] < child.time_start or times[-1] > child.time_end:
        raise ValueError("times must lie inside the child delivery interval")
    if len(coefficient_steps) != 2:
        raise ValueError("coefficient_steps must contain coarse and fine steps")
    coarse, fine = (float(value) for value in coefficient_steps)
    if not np.isfinite(coarse) or not np.isfinite(fine) or not (coarse > fine > 0.0):
        raise ValueError("coefficient steps must satisfy finite coarse > fine > 0")
    if not isinstance(resolution, int) or resolution < 21:
        raise ValueError("resolution must be an integer >= 21")
    if not np.isfinite(half_width) or half_width <= 0.0:
        raise ValueError("half_width must be positive and finite")
    if not np.isfinite(collar_start) or not (0.0 < collar_start < half_width):
        raise ValueError("collar_start must lie in (0, half_width)")

    basis = child.parent.profile_basis
    for mode in MODES:
        if mode not in basis.mode_indices:
            raise ValueError(f"parent is missing required existing Phi mode {mode}")
        index = basis.mode_indices.index(mode)
        if abs(basis.phi_coefficients[index]) + coarse > basis.coefficient_limit:
            raise ValueError(f"coefficient step leaves registered bound for Phi mode {mode}")

    axis = np.linspace(-half_width, half_width, resolution)
    baseline, baseline_rows = _feature_vector(child, axis, times, collar_start)
    coarse_matrix = _response_matrix(child, axis, times, collar_start, coarse)
    fine_matrix = _response_matrix(child, axis, times, collar_start, fine)

    singular_values = np.linalg.svd(fine_matrix, compute_uv=False)
    tolerance = np.finfo(float).eps * max(fine_matrix.shape) * singular_values[0]
    rank = int(np.count_nonzero(singular_values > tolerance))
    condition = float(singular_values[0] / max(singular_values[-1], np.finfo(float).tiny))
    refinement = float(
        np.linalg.norm(fine_matrix - coarse_matrix)
        / max(np.linalg.norm(fine_matrix), np.finfo(float).tiny)
    )

    # Project an early-only suppression signature onto the static two-mode span.
    target = np.asarray((-1.0, 0.0, 0.0), dtype=float)
    direction, _, _, _ = np.linalg.lstsq(fine_matrix, target, rcond=None)
    projected = fine_matrix @ direction
    projection_fraction = float(np.linalg.norm(projected) / np.linalg.norm(target))
    target_residual_fraction = float(np.linalg.norm(target - projected) / np.linalg.norm(target))
    projected_early_gain = float(-projected[0])
    projected_late_leak = float(np.linalg.norm(projected[1:]))

    # Turn the least-squares direction into one bounded 0.02-sized local trial,
    # then evaluate the actual nonlinear fingerprint rather than only J*delta.
    max_abs = float(np.max(np.abs(direction)))
    if max_abs <= np.finfo(float).tiny:
        scaled_direction = np.zeros_like(direction)
    else:
        scaled_direction = direction * (fine / max_abs)
    combined = _perturbed_child(
        child,
        {mode: float(scaled_direction[i]) for i, mode in enumerate(MODES)},
    )
    trial_features, trial_rows = _feature_vector(combined, axis, times, collar_start)
    actual_delta = trial_features - baseline
    predicted_delta = fine_matrix @ scaled_direction
    linearization_error = float(
        np.linalg.norm(actual_delta - predicted_delta)
        / max(np.linalg.norm(actual_delta), np.finfo(float).tiny)
    )
    static_relative_early_reduction = float(
        -actual_delta[0] / max(baseline[0], np.finfo(float).tiny)
    )

    temporal_screen = _affine_temporal_screen(
        child,
        times,
        fine_matrix,
        baseline,
        fine,
        predicted_delta,
    )

    return {
        "schema": "eq45_supported_early_radial_capacity_v1",
        "claim_scope": "local_visualization_morphology_capacity_only",
        "parent_sha256": child.parent_sha256,
        "supported_child_sha256": child.sha256,
        "modes": [list(mode) for mode in MODES],
        "parameter_count_added": 0,
        "existing_parameters_audited": 2,
        "times": list(times),
        "coefficient_steps": [coarse, fine],
        "resolution": resolution,
        "half_width": float(half_width),
        "collar_start": float(collar_start),
        "feature": "whole_grid_radial_outer_vorticity2_fraction",
        "baseline_feature": baseline.tolist(),
        "baseline_morphology": baseline_rows,
        "response_matrix": fine_matrix.tolist(),
        "singular_values": [float(value) for value in singular_values],
        "numerical_rank": rank,
        "condition_number": condition,
        "response_refinement_relative_change": refinement,
        "early_only_projection": {
            "target": target.tolist(),
            "least_squares_direction": direction.tolist(),
            "projected_signature": projected.tolist(),
            "projection_fraction": projection_fraction,
            "target_residual_fraction": target_residual_fraction,
            "projected_early_gain": projected_early_gain,
            "projected_late_leak_l2": projected_late_leak,
        },
        "bounded_static_local_trial": {
            "coefficient_deltas": {
                str(mode): float(scaled_direction[i]) for i, mode in enumerate(MODES)
            },
            "feature": trial_features.tolist(),
            "actual_delta": actual_delta.tolist(),
            "predicted_delta": predicted_delta.tolist(),
            "relative_early_reduction": static_relative_early_reduction,
            "linearization_relative_error": linearization_error,
            "morphology": trial_rows,
        },
        "affine_temporal_screen": temporal_screen,
        "interpretation": {
            "new_spatial_basis_added": False,
            "public_visual_target_fitted": False,
            "recommended_minimal_followup": "affine_existing_Phi(1,0)_screen_on_production_temporal_wrapper",
            "reason": "static_modes_span_early_signature_but_require_ill_conditioned_cancellation;_Phi(1,0)_affine_signature_is_more_efficient",
        },
        "truth_boundary": dict(_TRUTH_BOUNDARY),
    }
