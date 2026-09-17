"""Bounded vorticity-core capacity audit for the Eq. (4.5) eta^4 lift.

This module asks one narrow representation question: after the production
``lift_eq45_axial_eta4`` exposes only ``Phi(0,4)`` and ``F(0,4)``, how much
independent control do those two bounded coefficients provide over a robust
vorticity-core morphology measured from the public ``[u,v,w]`` interface?

It is an expression-capacity diagnostic only.  It does not fit OpenAI's hidden
profiles, change the frozen candidate, validate Navier--Stokes, or establish
visual correspondence.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from itertools import product
from typing import Any, Iterable

import numpy as np

from .constrained_eq45_axial_correction import lift_eq45_axial_eta4
from .constrained_eq45_candidate import Eq45VelocityCandidate


CLAIM_SCOPE = "eq45_eta4_vorticity_expression_capacity_only"
TRUTH_BOUNDARY = {
    "velocity_changed": False,
    "visualization_ready": False,
    "visual_correspondence_verified": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


@dataclass(frozen=True)
class VorticityCoreRms:
    time: float
    grid_size: int
    box_half_width: float
    radial_rms: float
    axial_rms: float
    aspect_ratio_z_over_r: float
    omega_max: float
    enstrophy_integral: float


def _validated_grid(grid_size: int, box_half_width: float) -> tuple[int, float]:
    if not isinstance(grid_size, (int, np.integer)) or grid_size < 17:
        raise ValueError("grid_size must be an integer >= 17")
    if grid_size % 2 == 0:
        raise ValueError("grid_size must be odd so the symmetry axis is sampled")
    box_half_width = float(box_half_width)
    if not np.isfinite(box_half_width) or box_half_width <= 0.0:
        raise ValueError("box_half_width must be positive and finite")
    return int(grid_size), box_half_width


def vorticity_core_rms(
    candidate: Eq45VelocityCandidate,
    time: float,
    *,
    grid_size: int = 33,
    box_half_width: float = 2.0,
) -> VorticityCoreRms:
    """Measure enstrophy-weighted core RMS scales from public velocity samples.

    A uniform Cartesian grid is sampled through ``candidate.at_points`` only.
    Curl is reconstructed independently with second-order centered differences
    on the one-cell interior.  The RMS scales are deliberately bulk quantities;
    they are less grid-sensitive than peak-vorticity or one chosen isosurface.
    """
    if not isinstance(candidate, Eq45VelocityCandidate):
        raise TypeError("candidate must be Eq45VelocityCandidate")
    grid_size, box_half_width = _validated_grid(grid_size, box_half_width)
    time = float(time)
    if not np.isfinite(time):
        raise ValueError("time must be finite")

    axis = np.linspace(-box_half_width, box_half_width, grid_size, dtype=float)
    spacing = float(axis[1] - axis[0])
    x, y, z = np.meshgrid(axis, axis, axis, indexing="ij")
    points = np.stack((x, y, z), axis=-1)
    velocity = np.asarray(candidate.at_points(points, time), dtype=float)
    if velocity.shape != points.shape or not np.all(np.isfinite(velocity)):
        raise ValueError("public candidate velocity must be finite with shape (...,3)")

    u = velocity[..., 0]
    v = velocity[..., 1]
    w = velocity[..., 2]
    two_h = 2.0 * spacing

    dw_dy = (w[1:-1, 2:, 1:-1] - w[1:-1, :-2, 1:-1]) / two_h
    dv_dz = (v[1:-1, 1:-1, 2:] - v[1:-1, 1:-1, :-2]) / two_h
    du_dz = (u[1:-1, 1:-1, 2:] - u[1:-1, 1:-1, :-2]) / two_h
    dw_dx = (w[2:, 1:-1, 1:-1] - w[:-2, 1:-1, 1:-1]) / two_h
    dv_dx = (v[2:, 1:-1, 1:-1] - v[:-2, 1:-1, 1:-1]) / two_h
    du_dy = (u[1:-1, 2:, 1:-1] - u[1:-1, :-2, 1:-1]) / two_h

    omega_x = dw_dy - dv_dz
    omega_y = du_dz - dw_dx
    omega_z = dv_dx - du_dy
    omega_sq = omega_x * omega_x + omega_y * omega_y + omega_z * omega_z
    if not np.all(np.isfinite(omega_sq)):
        raise RuntimeError("vorticity reconstruction became nonfinite")

    weight_sum = float(np.sum(omega_sq))
    if not np.isfinite(weight_sum) or weight_sum <= np.finfo(float).tiny:
        raise ValueError("candidate has numerically inactive vorticity on this grid")

    xi = x[1:-1, 1:-1, 1:-1]
    yi = y[1:-1, 1:-1, 1:-1]
    zi = z[1:-1, 1:-1, 1:-1]
    radial_sq = xi * xi + yi * yi
    radial_rms = float(np.sqrt(np.sum(omega_sq * radial_sq) / weight_sum))
    axial_rms = float(np.sqrt(np.sum(omega_sq * zi * zi) / weight_sum))
    if radial_rms <= np.finfo(float).tiny:
        raise ValueError("radial vorticity-core RMS is numerically zero")

    return VorticityCoreRms(
        time=time,
        grid_size=grid_size,
        box_half_width=box_half_width,
        radial_rms=radial_rms,
        axial_rms=axial_rms,
        aspect_ratio_z_over_r=axial_rms / radial_rms,
        omega_max=float(np.sqrt(np.max(omega_sq))),
        enstrophy_integral=float(weight_sum * spacing**3),
    )


def _validated_times(candidate: Eq45VelocityCandidate, times: Iterable[float]) -> tuple[float, ...]:
    values = tuple(float(value) for value in times)
    if len(values) < 3 or not all(np.isfinite(value) for value in values):
        raise ValueError("times must contain at least three finite values")
    if any(b <= a for a, b in zip(values, values[1:])):
        raise ValueError("times must be strictly increasing")
    if values[0] < candidate.time_start or values[-1] > candidate.time_end:
        raise ValueError("times must stay inside the candidate delivery interval")
    return values


def _validated_steps(limit: float, steps: Iterable[float]) -> tuple[float, ...]:
    values = tuple(float(value) for value in steps)
    if len(values) < 2 or not all(np.isfinite(value) and value > 0.0 for value in values):
        raise ValueError("coefficient_steps must contain at least two positive finite values")
    if any(b >= a for a, b in zip(values, values[1:])):
        raise ValueError("coefficient_steps must be strictly decreasing")
    if values[0] > limit:
        raise ValueError("coefficient step exceeds the registered coefficient bound")
    return values


def _metrics(candidate: Eq45VelocityCandidate, times: tuple[float, ...], grid_size: int, box: float):
    return tuple(
        vorticity_core_rms(candidate, time, grid_size=grid_size, box_half_width=box)
        for time in times
    )


def _relative_response(
    plus: tuple[VorticityCoreRms, ...],
    minus: tuple[VorticityCoreRms, ...],
    baseline: tuple[VorticityCoreRms, ...],
    step: float,
) -> tuple[np.ndarray, np.ndarray]:
    bulk: list[float] = []
    aspect: list[float] = []
    for p, m, b in zip(plus, minus, baseline, strict=True):
        bulk.extend(
            [
                (p.radial_rms - m.radial_rms) / (2.0 * step * b.radial_rms),
                (p.axial_rms - m.axial_rms) / (2.0 * step * b.axial_rms),
            ]
        )
        aspect.append(
            (p.aspect_ratio_z_over_r - m.aspect_ratio_z_over_r)
            / (2.0 * step * b.aspect_ratio_z_over_r)
        )
    return np.asarray(bulk, dtype=float), np.asarray(aspect, dtype=float)


def audit_eq45_eta4_vorticity_capacity(
    candidate: Eq45VelocityCandidate,
    *,
    times: Iterable[float] = (0.25, 0.5, 0.75),
    grid_size: int = 33,
    box_half_width: float = 2.0,
    coefficient_steps: Iterable[float] = (0.02, 0.01),
) -> dict[str, Any]:
    """Audit local rank and full-bound morphology reach of the two eta^4 modes.

    The finest central-difference response uses radial and axial RMS core scales
    as independent rows.  Aspect ratio is reported separately because it is a
    derived quotient of those two rows.  Full-bound corners answer a different
    question: how far can this *specific* two-mode lift move core anisotropy
    without relaxing the existing coefficient bounds?
    """
    if not isinstance(candidate, Eq45VelocityCandidate):
        raise TypeError("candidate must be Eq45VelocityCandidate")
    grid_size, box_half_width = _validated_grid(grid_size, box_half_width)
    times = _validated_times(candidate, times)
    limit = float(candidate.profile_basis.coefficient_limit)
    steps = _validated_steps(limit, coefficient_steps)

    baseline = _metrics(candidate, times, grid_size, box_half_width)
    response_matrices: list[np.ndarray] = []
    aspect_responses: list[np.ndarray] = []

    for step in steps:
        columns: list[np.ndarray] = []
        aspect_columns: list[np.ndarray] = []
        for channel in ("phi_eta4", "swirl_eta4"):
            kwargs_plus = {channel: step}
            kwargs_minus = {channel: -step}
            plus = _metrics(lift_eq45_axial_eta4(candidate, **kwargs_plus), times, grid_size, box_half_width)
            minus = _metrics(lift_eq45_axial_eta4(candidate, **kwargs_minus), times, grid_size, box_half_width)
            bulk, aspect = _relative_response(plus, minus, baseline, step)
            columns.append(bulk)
            aspect_columns.append(aspect)
        response_matrices.append(np.column_stack(columns))
        aspect_responses.append(np.column_stack(aspect_columns))

    response = response_matrices[-1]
    singular_values = np.linalg.svd(response, compute_uv=False)
    tolerance = max(response.shape) * np.finfo(float).eps * singular_values[0]
    rank = int(np.sum(singular_values > tolerance))
    condition = float(singular_values[0] / singular_values[-1]) if rank == 2 else float("inf")
    denom = float(np.linalg.norm(response))
    refinement_change = float(np.linalg.norm(response_matrices[0] - response) / denom)

    endpoints: dict[str, Any] = {}
    for channel, phi_value, swirl_value in (
        ("phi_minus", -limit, 0.0),
        ("phi_plus", limit, 0.0),
        ("swirl_minus", 0.0, -limit),
        ("swirl_plus", 0.0, limit),
    ):
        current = _metrics(
            lift_eq45_axial_eta4(candidate, phi_eta4=phi_value, swirl_eta4=swirl_value),
            times,
            grid_size,
            box_half_width,
        )
        endpoints[channel] = [asdict(item) for item in current]

    corners: list[dict[str, Any]] = []
    baseline_aspect = np.asarray([item.aspect_ratio_z_over_r for item in baseline])
    for phi_value, swirl_value in product((-limit, limit), repeat=2):
        current = _metrics(
            lift_eq45_axial_eta4(candidate, phi_eta4=phi_value, swirl_eta4=swirl_value),
            times,
            grid_size,
            box_half_width,
        )
        aspect = np.asarray([item.aspect_ratio_z_over_r for item in current])
        corners.append(
            {
                "phi_eta4": float(phi_value),
                "swirl_eta4": float(swirl_value),
                "aspect_ratio_z_over_r": aspect.tolist(),
                "relative_aspect_change": (aspect / baseline_aspect - 1.0).tolist(),
                "mean_aspect_ratio": float(np.mean(aspect)),
            }
        )
    max_corner = max(corners, key=lambda item: item["mean_aspect_ratio"])
    min_corner = min(corners, key=lambda item: item["mean_aspect_ratio"])

    return {
        "claim_scope": CLAIM_SCOPE,
        "candidate_sha256": candidate.sha256,
        "times": list(times),
        "grid_size": grid_size,
        "box_half_width": box_half_width,
        "coefficient_limit": limit,
        "coefficient_steps": list(steps),
        "channels": ["Phi(0,4)", "F(0,4)"],
        "baseline": [asdict(item) for item in baseline],
        "finest_relative_bulk_response": response.tolist(),
        "finest_relative_aspect_response": aspect_responses[-1].tolist(),
        "singular_values": singular_values.tolist(),
        "numerical_rank": rank,
        "condition_number": condition,
        "coarse_to_finest_response_frobenius_change": refinement_change,
        "single_channel_bound_endpoints": endpoints,
        "joint_bound_corners": corners,
        "max_mean_aspect_corner": max_corner,
        "min_mean_aspect_corner": min_corner,
        "truth_boundary": dict(TRUTH_BOUNDARY),
        "interpretation": (
            "bounded local expression-capacity evidence only; no target image was fitted and "
            "no scientific readiness state is promoted"
        ),
    }
