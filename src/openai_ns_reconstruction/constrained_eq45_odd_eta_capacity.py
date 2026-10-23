"""Target-free capacity audit for dormant odd-eta Eq. (4.5) profile modes.

The frozen Eq45 seed uses zero coefficients for the existing ``(0,1)`` Phi/F
modes, which enforces a strict z-reflection parity.  This module asks the narrow
representation question needed by the visualization lane: can thawing exactly
one already-registered odd-eta coefficient break the relevant morphology
constraint without broad basis growth?

All fingerprints are reconstructed from ``candidate.at_points(...)->[u,v,w]``.
No target image is fitted, the frozen candidate is not modified, and no PDE,
visual-correspondence, paper-exact, or hidden-profile claim is made.
"""
from __future__ import annotations

from dataclasses import replace
from typing import Any, Iterable

import numpy as np

from .constrained_eq45_candidate import Eq45VelocityCandidate


CLAIM_SCOPE = "eq45_existing_odd_eta_parity_capacity_only"
PARAMETER_LABELS = ("Phi(0,1)", "F(0,1)")
METRIC_LABELS = (
    "midplane_radial_velocity",
    "odd_z_axial_velocity",
    "odd_z_swirl_velocity",
)
TRUTH_BOUNDARY = {
    "velocity_changed": False,
    "odd_eta_mode_activated_in_production": False,
    "visualization_ready": False,
    "visual_correspondence_verified": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


def _validated_sequence(name: str, values: Iterable[float]) -> tuple[float, ...]:
    out = tuple(float(value) for value in values)
    if not out or not all(np.isfinite(value) for value in out):
        raise ValueError(f"{name} must contain finite values")
    return out


def _with_odd_eta_delta(
    candidate: Eq45VelocityCandidate, *, channel: str, delta: float
) -> Eq45VelocityCandidate:
    if not np.isfinite(delta):
        raise ValueError("coefficient delta must be finite")
    basis = candidate.profile_basis
    try:
        index = basis.mode_indices.index((0, 1))
    except ValueError as exc:
        raise ValueError("candidate basis does not contain the required (0,1) mode") from exc

    if channel == "phi":
        coefficients = list(basis.phi_coefficients)
        coefficients[index] += float(delta)
        updated = replace(basis, phi_coefficients=tuple(coefficients))
    elif channel == "swirl":
        coefficients = list(basis.swirl_coefficients)
        coefficients[index] += float(delta)
        updated = replace(basis, swirl_coefficients=tuple(coefficients))
    else:
        raise ValueError("channel must be 'phi' or 'swirl'")
    return replace(candidate, profile_basis=updated)


def _cylindrical_components(points: np.ndarray, velocity: np.ndarray):
    x = points[..., 0]
    y = points[..., 1]
    radius = np.hypot(x, y)
    if np.any(radius <= 0.0):
        raise ValueError("parity probes must stay away from the symmetry axis")
    radial = (x * velocity[..., 0] + y * velocity[..., 1]) / radius
    swirl = (-y * velocity[..., 0] + x * velocity[..., 1]) / radius
    return radial, swirl, velocity[..., 2]


def _fingerprint_raw(
    candidate: Eq45VelocityCandidate,
    *,
    times: tuple[float, ...],
    radii: tuple[float, ...],
    z_offsets: tuple[float, ...],
    azimuth_count: int,
) -> tuple[np.ndarray, float]:
    angles = 2.0 * np.pi * np.arange(azimuth_count, dtype=float) / azimuth_count
    midplane_radial: list[float] = []
    odd_axial: list[float] = []
    odd_swirl: list[float] = []
    speed_samples: list[float] = []

    for time in times:
        for radius in radii:
            mid = np.column_stack(
                (
                    radius * np.cos(angles),
                    radius * np.sin(angles),
                    np.zeros_like(angles),
                )
            )
            mid_velocity = candidate.at_points(mid, time)
            radial, _, _ = _cylindrical_components(mid, mid_velocity)
            midplane_radial.extend(radial.tolist())
            speed_samples.extend(np.linalg.norm(mid_velocity, axis=-1).tolist())

            for z_offset in z_offsets:
                plus = np.column_stack(
                    (
                        radius * np.cos(angles),
                        radius * np.sin(angles),
                        np.full_like(angles, z_offset),
                    )
                )
                minus = plus.copy()
                minus[:, 2] *= -1.0
                velocity_plus = candidate.at_points(plus, time)
                velocity_minus = candidate.at_points(minus, time)
                _, swirl_plus, axial_plus = _cylindrical_components(plus, velocity_plus)
                _, swirl_minus, axial_minus = _cylindrical_components(minus, velocity_minus)
                odd_axial.extend((0.5 * (axial_plus - axial_minus)).tolist())
                odd_swirl.extend((0.5 * (swirl_plus - swirl_minus)).tolist())
                speed_samples.extend(np.linalg.norm(velocity_plus, axis=-1).tolist())
                speed_samples.extend(np.linalg.norm(velocity_minus, axis=-1).tolist())

    metrics = np.asarray(
        [np.mean(midplane_radial), np.mean(odd_axial), np.mean(odd_swirl)],
        dtype=float,
    )
    reference_speed_rms = float(np.sqrt(np.mean(np.square(speed_samples))))
    if not np.all(np.isfinite(metrics)) or not np.isfinite(reference_speed_rms):
        raise RuntimeError("odd-eta parity fingerprint became nonfinite")
    if reference_speed_rms <= np.finfo(float).tiny:
        raise ValueError("baseline public velocity is numerically trivial")
    return metrics, reference_speed_rms


def audit_eq45_odd_eta_capacity(
    candidate: Eq45VelocityCandidate,
    *,
    times: Iterable[float] = (0.25, 0.5, 0.75),
    radii: Iterable[float] = (0.35, 0.70),
    z_offsets: Iterable[float] = (0.25, 0.50),
    azimuth_count: int = 8,
    coefficient_steps: Iterable[float] = (0.04, 0.02),
) -> dict[str, Any]:
    """Measure local parity-breaking capacity of existing ``(0,1)`` Phi/F modes."""
    if not isinstance(candidate, Eq45VelocityCandidate):
        raise TypeError("candidate must be Eq45VelocityCandidate")
    times = _validated_sequence("times", times)
    radii = _validated_sequence("radii", radii)
    z_offsets = _validated_sequence("z_offsets", z_offsets)
    steps = _validated_sequence("coefficient_steps", coefficient_steps)
    if any(b >= a for a, b in zip(steps, steps[1:])) or any(step <= 0.0 for step in steps):
        raise ValueError("coefficient_steps must be positive and strictly decreasing")
    if times[0] < candidate.time_start or times[-1] > candidate.time_end:
        raise ValueError("times must stay inside the candidate delivery interval")
    if any(radius <= 0.0 for radius in radii) or any(offset <= 0.0 for offset in z_offsets):
        raise ValueError("radii and z_offsets must be positive")
    if not isinstance(azimuth_count, (int, np.integer)) or azimuth_count < 4:
        raise ValueError("azimuth_count must be an integer >= 4")

    basis = candidate.profile_basis
    limit = float(basis.coefficient_limit)
    if steps[0] > limit:
        raise ValueError("coefficient step exceeds registered coefficient bound")
    odd_index = basis.mode_indices.index((0, 1))
    if abs(basis.phi_coefficients[odd_index]) > 1e-14 or abs(basis.swirl_coefficients[odd_index]) > 1e-14:
        raise ValueError("audit requires a frozen candidate with dormant (0,1) Phi/F modes")

    baseline_raw, reference_speed = _fingerprint_raw(
        candidate,
        times=times,
        radii=radii,
        z_offsets=z_offsets,
        azimuth_count=int(azimuth_count),
    )
    baseline = baseline_raw / reference_speed

    matrices: list[np.ndarray] = []
    for step in steps:
        columns = []
        for channel in ("phi", "swirl"):
            plus = _with_odd_eta_delta(candidate, channel=channel, delta=step)
            minus = _with_odd_eta_delta(candidate, channel=channel, delta=-step)
            plus_raw, _ = _fingerprint_raw(
                plus, times=times, radii=radii, z_offsets=z_offsets, azimuth_count=int(azimuth_count)
            )
            minus_raw, _ = _fingerprint_raw(
                minus, times=times, radii=radii, z_offsets=z_offsets, azimuth_count=int(azimuth_count)
            )
            columns.append((plus_raw - minus_raw) / (2.0 * step * reference_speed))
        matrices.append(np.column_stack(columns))

    response = matrices[-1]
    singular_values = np.linalg.svd(response, compute_uv=False)
    tolerance = max(response.shape) * np.finfo(float).eps * singular_values[0]
    rank = int(np.sum(singular_values > tolerance))
    condition_number = float(singular_values[0] / singular_values[-1]) if rank == 2 else float("inf")
    response_norm = float(np.linalg.norm(response))
    refinement_change = float(np.linalg.norm(matrices[0] - response) / response_norm)

    phi_cross_talk = float(abs(response[2, 0]) / max(np.linalg.norm(response[:, 0]), np.finfo(float).tiny))
    swirl_cross_talk = float(np.linalg.norm(response[:2, 1]) / max(np.linalg.norm(response[:, 1]), np.finfo(float).tiny))

    return {
        "claim_scope": CLAIM_SCOPE,
        "candidate_sha256": candidate.sha256,
        "times": list(times),
        "radii": list(radii),
        "z_offsets": list(z_offsets),
        "azimuth_count": int(azimuth_count),
        "coefficient_limit": limit,
        "coefficient_steps": list(steps),
        "parameter_labels": list(PARAMETER_LABELS),
        "metric_labels": list(METRIC_LABELS),
        "baseline_reference_speed_rms": reference_speed,
        "baseline_normalized_parity_fingerprint": baseline.tolist(),
        "finest_normalized_response_per_unit_coefficient": response.tolist(),
        "singular_values": singular_values.tolist(),
        "numerical_rank": rank,
        "condition_number": condition_number,
        "step_refinement_relative_change": refinement_change,
        "phi_to_odd_swirl_cross_talk_fraction": phi_cross_talk,
        "swirl_to_poloidal_cross_talk_fraction": swirl_cross_talk,
        "recommended_minimal_extension": (
            "If a public target trace requires nonzero midplane poloidal/radial flow or z-odd axial motion, thaw only bounded Phi(0,1) first. "
            "Thaw bounded F(0,1) only if the target independently requires z-odd swirl. Do not thaw a broad odd-eta block without new evidence."
        ),
        "truth_boundary": dict(TRUTH_BOUNDARY),
    }
