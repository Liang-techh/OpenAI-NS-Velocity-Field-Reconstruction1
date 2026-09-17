"""Local capacity audit for the supported Eq45 early radial vorticity branch.

The whole-domain envelope audit shows that the physical support transform can
produce a larger early-time radial high-vorticity envelope than the untapered
parent.  This module asks a narrower representation question: can the two
*already present* bounded poloidal modes ``Phi(1,0)`` and ``Phi(1,2)`` move a
smooth radial-collar vorticity fingerprint mainly at the early delivery time,
or is a time-shaped coefficient direction needed?

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

    The differentiable-ish morphology feature is the fraction of whole-grid
    ``|curl u|^2`` lying at physical radius ``r > collar_start``.  The region
    itself is fixed, avoiding the coefficient-dependent voxel jumps of a hard
    vorticity superlevel extent.  A least-squares projection of the idealized
    time signature ``[-1,0,0]`` onto the two static response columns measures
    how selectively these existing modes can suppress the early radial branch
    while leaving later delivery times unchanged to first order.
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
        "bounded_local_trial": {
            "coefficient_deltas": {
                str(mode): float(scaled_direction[i]) for i, mode in enumerate(MODES)
            },
            "feature": trial_features.tolist(),
            "actual_delta": actual_delta.tolist(),
            "predicted_delta": predicted_delta.tolist(),
            "linearization_relative_error": linearization_error,
            "morphology": trial_rows,
        },
        "interpretation": {
            "new_spatial_basis_added": False,
            "public_visual_target_fitted": False,
            "purpose": "screen_existing_static_Phi_modes_before_temporal_or_basis_growth",
        },
        "truth_boundary": dict(_TRUTH_BOUNDARY),
    }
