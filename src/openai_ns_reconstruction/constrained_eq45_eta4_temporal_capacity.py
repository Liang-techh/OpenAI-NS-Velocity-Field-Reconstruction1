"""Candidate-specific temporal capacity audit for the Eq. (4.5) eta^4 lift.

The static production lift exposes two bounded spatial correction coefficients,
``Phi(0,4)`` and ``F(0,4)``.  This module asks one narrower follow-up question:
would adding exactly one affine-in-time slope to ``Phi(0,4)`` provide a genuinely
new control direction over the public-velocity vorticity-core aspect evolution?

This is a local expression-capacity diagnostic.  It does not implement a new
production velocity artifact, fit a target image, change the frozen candidate,
or establish visualization correspondence or Navier--Stokes validity.
"""
from __future__ import annotations

from typing import Any, Iterable

import numpy as np

from .constrained_eq45_axial_correction import lift_eq45_axial_eta4
from .constrained_eq45_candidate import Eq45VelocityCandidate
from .constrained_eq45_eta4_vorticity_capacity import vorticity_core_rms


CLAIM_SCOPE = "eq45_eta4_affine_temporal_expression_capacity_only"
PARAMETER_LABELS = ("Phi(0,4):static", "F(0,4):static", "Phi(0,4):time_slope")
TRUTH_BOUNDARY = {
    "velocity_changed": False,
    "time_dependent_lift_implemented": False,
    "visualization_ready": False,
    "visual_correspondence_verified": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


def _validate_grid(grid_size: int, box_half_width: float) -> tuple[int, float]:
    if not isinstance(grid_size, (int, np.integer)) or grid_size < 17:
        raise ValueError("grid_size must be an integer >= 17")
    if grid_size % 2 == 0:
        raise ValueError("grid_size must be odd so the symmetry axis is sampled")
    box_half_width = float(box_half_width)
    if not np.isfinite(box_half_width) or box_half_width <= 0.0:
        raise ValueError("box_half_width must be positive and finite")
    return int(grid_size), box_half_width


def _validate_times(
    candidate: Eq45VelocityCandidate, times: Iterable[float]
) -> tuple[tuple[float, ...], np.ndarray]:
    values = tuple(float(value) for value in times)
    if len(values) < 3 or not all(np.isfinite(value) for value in values):
        raise ValueError("times must contain at least three finite values")
    if any(b <= a for a, b in zip(values, values[1:])):
        raise ValueError("times must be strictly increasing")
    if values[0] < candidate.time_start or values[-1] > candidate.time_end:
        raise ValueError("times must stay inside the candidate delivery interval")

    midpoint = 0.5 * (values[0] + values[-1])
    half_span = 0.5 * (values[-1] - values[0])
    if half_span <= 0.0:
        raise ValueError("time span must be positive")
    tau = (np.asarray(values, dtype=float) - midpoint) / half_span
    return values, tau


def _validate_steps(limit: float, coefficient_steps: Iterable[float]) -> tuple[float, ...]:
    steps = tuple(float(value) for value in coefficient_steps)
    if len(steps) < 2 or not all(np.isfinite(value) and value > 0.0 for value in steps):
        raise ValueError("coefficient_steps must contain at least two positive finite values")
    if any(b >= a for a, b in zip(steps, steps[1:])):
        raise ValueError("coefficient_steps must be strictly decreasing")
    if steps[0] > limit:
        raise ValueError("coefficient step exceeds the registered coefficient bound")
    return steps


def _aspect_vector(
    candidate: Eq45VelocityCandidate,
    times: tuple[float, ...],
    *,
    grid_size: int,
    box_half_width: float,
) -> np.ndarray:
    return np.asarray(
        [
            vorticity_core_rms(
                candidate,
                time,
                grid_size=grid_size,
                box_half_width=box_half_width,
            ).aspect_ratio_z_over_r
            for time in times
        ],
        dtype=float,
    )


def _static_response_column(
    candidate: Eq45VelocityCandidate,
    times: tuple[float, ...],
    baseline: np.ndarray,
    *,
    step: float,
    channel: str,
    grid_size: int,
    box_half_width: float,
) -> np.ndarray:
    if channel == "phi":
        plus = lift_eq45_axial_eta4(candidate, phi_eta4=step)
        minus = lift_eq45_axial_eta4(candidate, phi_eta4=-step)
    elif channel == "swirl":
        plus = lift_eq45_axial_eta4(candidate, swirl_eta4=step)
        minus = lift_eq45_axial_eta4(candidate, swirl_eta4=-step)
    else:  # pragma: no cover - internal programming guard
        raise ValueError("unsupported static response channel")

    plus_aspect = _aspect_vector(
        plus, times, grid_size=grid_size, box_half_width=box_half_width
    )
    minus_aspect = _aspect_vector(
        minus, times, grid_size=grid_size, box_half_width=box_half_width
    )
    return (plus_aspect - minus_aspect) / (2.0 * step * baseline)


def _temporal_phi_response_column(
    candidate: Eq45VelocityCandidate,
    times: tuple[float, ...],
    tau: np.ndarray,
    baseline: np.ndarray,
    *,
    step: float,
    grid_size: int,
    box_half_width: float,
) -> np.ndarray:
    plus_values: list[float] = []
    minus_values: list[float] = []
    for time, tau_value in zip(times, tau, strict=True):
        plus = lift_eq45_axial_eta4(candidate, phi_eta4=float(tau_value * step))
        minus = lift_eq45_axial_eta4(candidate, phi_eta4=float(-tau_value * step))
        plus_values.append(
            vorticity_core_rms(
                plus,
                time,
                grid_size=grid_size,
                box_half_width=box_half_width,
            ).aspect_ratio_z_over_r
        )
        minus_values.append(
            vorticity_core_rms(
                minus,
                time,
                grid_size=grid_size,
                box_half_width=box_half_width,
            ).aspect_ratio_z_over_r
        )
    return (
        np.asarray(plus_values, dtype=float) - np.asarray(minus_values, dtype=float)
    ) / (2.0 * step * baseline)


def audit_eq45_eta4_temporal_capacity(
    candidate: Eq45VelocityCandidate,
    *,
    times: Iterable[float] = (0.25, 0.5, 0.75),
    grid_size: int = 17,
    box_half_width: float = 2.0,
    coefficient_steps: Iterable[float] = (0.04, 0.02),
) -> dict[str, Any]:
    """Measure the local novelty of one affine ``Phi(0,4)`` time slope.

    The normalized time coordinate is ``tau=(t-midpoint)/half_span``.  The two
    existing static eta^4 coefficients form the first two response columns; the
    proposed third parameter uses ``Phi(0,4)=tau*slope``.  Every morphology value
    is reconstructed from ``candidate.at_points(...)->[u,v,w]`` through the
    existing independent Cartesian vorticity-core metric.

    The report is target-free: rank, condition number and the fraction of the
    slope response outside the static two-column span quantify expression
    capacity only.  No hidden visual target or PDE residual is fitted.
    """
    if not isinstance(candidate, Eq45VelocityCandidate):
        raise TypeError("candidate must be Eq45VelocityCandidate")
    grid_size, box_half_width = _validate_grid(grid_size, box_half_width)
    times, tau = _validate_times(candidate, times)
    limit = float(candidate.profile_basis.coefficient_limit)
    steps = _validate_steps(limit, coefficient_steps)
    if float(np.max(np.abs(tau))) > 1.0 + 32.0 * np.finfo(float).eps:
        raise RuntimeError("normalized time coordinate escaped [-1,1]")

    baseline = _aspect_vector(
        candidate, times, grid_size=grid_size, box_half_width=box_half_width
    )
    if not np.all(np.isfinite(baseline)) or np.any(baseline <= 0.0):
        raise ValueError("baseline vorticity-core aspect must be finite and positive")

    response_matrices: list[np.ndarray] = []
    for step in steps:
        phi = _static_response_column(
            candidate,
            times,
            baseline,
            step=step,
            channel="phi",
            grid_size=grid_size,
            box_half_width=box_half_width,
        )
        swirl = _static_response_column(
            candidate,
            times,
            baseline,
            step=step,
            channel="swirl",
            grid_size=grid_size,
            box_half_width=box_half_width,
        )
        temporal_phi = _temporal_phi_response_column(
            candidate,
            times,
            tau,
            baseline,
            step=step,
            grid_size=grid_size,
            box_half_width=box_half_width,
        )
        response_matrices.append(np.column_stack((phi, swirl, temporal_phi)))

    response = response_matrices[-1]
    singular_values = np.linalg.svd(response, compute_uv=False)
    tolerance = max(response.shape) * np.finfo(float).eps * singular_values[0]
    rank = int(np.sum(singular_values > tolerance))
    condition = (
        float(singular_values[0] / singular_values[-1])
        if rank == len(PARAMETER_LABELS)
        else float("inf")
    )

    static = response[:, :2]
    static_singular_values = np.linalg.svd(static, compute_uv=False)
    static_tolerance = max(static.shape) * np.finfo(float).eps * static_singular_values[0]
    static_rank = int(np.sum(static_singular_values > static_tolerance))
    static_condition = (
        float(static_singular_values[0] / static_singular_values[-1])
        if static_rank == 2
        else float("inf")
    )

    slope = response[:, 2]
    static_fit_coefficients, *_ = np.linalg.lstsq(static, slope, rcond=None)
    orthogonal = slope - static @ static_fit_coefficients
    slope_norm = float(np.linalg.norm(slope))
    if slope_norm <= np.finfo(float).tiny:
        raise ValueError("affine temporal response is numerically inactive")
    novelty_fraction = float(np.linalg.norm(orthogonal) / slope_norm)

    response_norm = float(np.linalg.norm(response))
    if response_norm <= np.finfo(float).tiny:
        raise ValueError("eta4 temporal response is numerically inactive")
    refinement_change = float(
        np.linalg.norm(response_matrices[0] - response) / response_norm
    )

    return {
        "claim_scope": CLAIM_SCOPE,
        "candidate_sha256": candidate.sha256,
        "times": list(times),
        "normalized_time_tau": tau.tolist(),
        "grid_size": grid_size,
        "box_half_width": box_half_width,
        "coefficient_limit": limit,
        "coefficient_steps": list(steps),
        "parameter_labels": list(PARAMETER_LABELS),
        "existing_static_parameter_count": 2,
        "proposed_new_parameter_count": 1,
        "baseline_aspect_ratio_z_over_r": baseline.tolist(),
        "finest_relative_aspect_response": response.tolist(),
        "singular_values": singular_values.tolist(),
        "numerical_rank": rank,
        "condition_number": condition,
        "static_singular_values": static_singular_values.tolist(),
        "static_numerical_rank": static_rank,
        "static_condition_number": static_condition,
        "temporal_slope_novelty_fraction": novelty_fraction,
        "temporal_slope_static_projection_coefficients": static_fit_coefficients.tolist(),
        "coarse_to_finest_response_frobenius_change": refinement_change,
        "truth_boundary": dict(TRUTH_BOUNDARY),
        "interpretation": (
            "local target-free expression-capacity evidence only; one affine Phi(0,4) "
            "time slope is screened against the two existing static eta4 controls, but "
            "no time-dependent production candidate or scientific readiness state is created"
        ),
    }
