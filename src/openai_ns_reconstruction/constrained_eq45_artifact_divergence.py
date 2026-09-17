"""Independent divergence audit for a frozen public Eq. (4.5) velocity artifact.

The audit deliberately treats the candidate as a black-box public field. Every
velocity sample is obtained through ``field.at_points(points, time)``; no
optimizer tensors, profile derivatives, or streamfunction identities are reused
by the finite-difference operator.

This module checks only sampled incompressibility and spatial-difference
convergence. It does not evaluate the momentum equation and must never promote a
candidate to PDE-valid, paper-exact, or OpenAI-field-identified status.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any, Iterable

import numpy as np

from .constrained_eq45_candidate import Eq45VelocityCandidate


CLAIM_SCOPE = "independent_public_velocity_divergence_only"
MIN_ACTIVITY_RMS = 1e-12


def _validated_points_times(points, times) -> tuple[np.ndarray, np.ndarray]:
    pts = np.asarray(points, dtype=float)
    if pts.ndim != 2 or pts.shape[1] != 3 or pts.shape[0] == 0:
        raise ValueError("points must have shape (n, 3) with n > 0")
    if not np.all(np.isfinite(pts)):
        raise ValueError("points must be finite")

    t = np.asarray(times, dtype=float)
    if not np.all(np.isfinite(t)):
        raise ValueError("times must be finite")
    try:
        t = np.broadcast_to(t, (pts.shape[0],)).copy()
    except ValueError as exc:
        raise ValueError("times must be scalar or broadcast to one value per point") from exc
    return pts, t


def _validated_steps(steps: Iterable[float]) -> tuple[float, ...]:
    values = tuple(float(step) for step in steps)
    if len(values) < 3:
        raise ValueError("at least three spatial derivative steps are required")
    if not all(np.isfinite(step) and step > 0.0 for step in values):
        raise ValueError("spatial derivative steps must be finite and positive")
    if not all(values[i + 1] < values[i] for i in range(len(values) - 1)):
        raise ValueError("spatial derivative steps must be strictly decreasing")
    return values


def _public_velocity(field, points: np.ndarray, times: np.ndarray) -> np.ndarray:
    at_points = getattr(field, "at_points", None)
    if not callable(at_points):
        raise TypeError("field must expose a callable at_points(points, time) interface")
    values = np.asarray(at_points(points, times), dtype=float)
    if values.shape != points.shape:
        raise ValueError("public velocity output must have shape (n, 3)")
    if not np.all(np.isfinite(values)):
        raise ValueError("public velocity output must be finite")
    return values


def _finite_difference_gradient(
    field,
    points: np.ndarray,
    times: np.ndarray,
    step: float,
) -> np.ndarray:
    gradient = np.empty((points.shape[0], 3, 3), dtype=float)
    for axis in range(3):
        offset = np.zeros(3, dtype=float)
        offset[axis] = step
        plus = _public_velocity(field, points + offset, times)
        minus = _public_velocity(field, points - offset, times)
        gradient[:, :, axis] = (plus - minus) / (2.0 * step)
    return gradient


def audit_public_velocity_divergence(
    field,
    points,
    times,
    *,
    steps=(0.08, 0.04, 0.02),
) -> dict[str, Any]:
    """Audit sampled divergence with an operator independent of candidate internals.

    ``normalized_rms`` is the divergence RMS divided by the RMS Frobenius norm
    of the independently sampled velocity gradient at the same level. This is a
    dimensionless diagnostic normalization, not an acceptance threshold.
    """

    pts, t = _validated_points_times(points, times)
    derivative_steps = _validated_steps(steps)

    base_velocity = _public_velocity(field, pts, t)
    velocity_rms = float(np.sqrt(np.mean(np.sum(base_velocity * base_velocity, axis=1))))
    if not np.isfinite(velocity_rms) or velocity_rms <= MIN_ACTIVITY_RMS:
        raise ValueError("reference public velocity is numerically inactive")

    levels: list[dict[str, float]] = []
    for step in derivative_steps:
        gradient = _finite_difference_gradient(field, pts, t, step)
        divergence = gradient[:, 0, 0] + gradient[:, 1, 1] + gradient[:, 2, 2]
        max_abs = float(np.max(np.abs(divergence)))
        rms = float(np.sqrt(np.mean(divergence * divergence)))
        gradient_rms = float(
            np.sqrt(np.mean(np.sum(gradient * gradient, axis=(1, 2))))
        )
        if not np.isfinite(gradient_rms) or gradient_rms <= 0.0:
            raise ValueError("sampled velocity-gradient norm is not usable")
        levels.append(
            {
                "step": step,
                "max_abs": max_abs,
                "rms": rms,
                "gradient_frobenius_rms": gradient_rms,
                "normalized_rms": rms / gradient_rms,
            }
        )

    observed_orders: list[float | None] = []
    for coarse, fine in zip(levels[:-1], levels[1:]):
        if coarse["rms"] <= 0.0 or fine["rms"] <= 0.0:
            observed_orders.append(None)
            continue
        observed_orders.append(
            float(
                np.log(coarse["rms"] / fine["rms"])
                / np.log(coarse["step"] / fine["step"])
            )
        )

    return {
        "claim_scope": CLAIM_SCOPE,
        "operator": "centered_second_order_cartesian_space",
        "velocity_access": "at_points_only",
        "point_count": int(pts.shape[0]),
        "velocity_rms": velocity_rms,
        "levels": levels,
        "observed_rms_orders": observed_orders,
        "pde_validated": False,
        "visual_correspondence_verified": False,
        "paper_exact": False,
        "openai_field_identified": False,
        "blowup_proved": False,
    }


def audit_eq45_candidate_artifact_divergence(
    candidate_path: str | Path,
    points,
    times,
    *,
    steps=(0.08, 0.04, 0.02),
) -> dict[str, Any]:
    """Reload one serialized Eq. (4.5) candidate and audit its public velocity."""

    path = Path(candidate_path)
    field = Eq45VelocityCandidate.load_json(path)
    report = audit_public_velocity_divergence(field, points, times, steps=steps)
    report["candidate_sha256"] = field.sha256
    report["artifact_reloaded"] = True
    report["physical_support_validated"] = False
    report["full_momentum_residual_assessed"] = False
    return report
