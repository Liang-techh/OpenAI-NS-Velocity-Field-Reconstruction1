"""Independent finite-window boundary/support diagnostics for constrained candidates.

This module deliberately does not know a candidate formula. It accepts a velocity
callable so it can be used as an independent CR007 post-processing path rather
than reusing a candidate's internal support or energy logic.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np

VelocityFn = Callable[[np.ndarray, float], np.ndarray]


@dataclass(frozen=True)
class BoundarySupportReport:
    time: float
    support_half_width: float
    probe_half_width: float
    points_per_axis: int
    face_max_speed: float
    exterior_shell_max_speed: float
    face_rms_speed: float
    exterior_shell_rms_speed: float
    sampled_support_pass: bool
    tolerance: float


def _validate_scalar(name: str, value: float, *, positive: bool = False) -> float:
    value = float(value)
    if not np.isfinite(value):
        raise ValueError(f"{name} must be finite")
    if positive and value <= 0:
        raise ValueError(f"{name} must be positive")
    return value


def _eval_velocity(velocity: VelocityFn, points: np.ndarray, time: float) -> np.ndarray:
    values = np.asarray(velocity(points, time), dtype=float)
    if values.shape != points.shape:
        raise ValueError("velocity must return shape (N,3) for input shape (N,3)")
    if not np.all(np.isfinite(values)):
        raise ValueError("velocity returned nonfinite values")
    return values


def _face_points(half_width: float, n: int) -> np.ndarray:
    axis = np.linspace(-half_width, half_width, n)
    a, b = np.meshgrid(axis, axis, indexing="ij")
    faces = []
    for dim in range(3):
        for sign in (-1.0, 1.0):
            p = np.zeros((n, n, 3), dtype=float)
            p[..., dim] = sign * half_width
            p[..., (dim + 1) % 3] = a
            p[..., (dim + 2) % 3] = b
            faces.append(p.reshape(-1, 3))
    return np.concatenate(faces, axis=0)


def _exterior_shell_points(inner: float, outer: float, n: int) -> np.ndarray:
    axis = np.linspace(-outer, outer, n)
    x, y, z = np.meshgrid(axis, axis, axis, indexing="ij")
    points = np.stack((x, y, z), axis=-1).reshape(-1, 3)
    mask = np.max(np.abs(points), axis=1) > inner
    return points[mask]


def boundary_support_report(
    velocity: VelocityFn,
    *,
    time: float,
    support_half_width: float,
    probe_half_width: float,
    points_per_axis: int = 17,
    tolerance: float = 1e-12,
) -> BoundarySupportReport:
    """Probe the declared support boundary and a strictly larger exterior shell.

    ``sampled_support_pass`` is only a sampling statement. Even exact zeros on
    the sampled shell are not promoted to an analytic compact-support proof.
    """
    time = _validate_scalar("time", time)
    support_half_width = _validate_scalar(
        "support_half_width", support_half_width, positive=True
    )
    probe_half_width = _validate_scalar(
        "probe_half_width", probe_half_width, positive=True
    )
    tolerance = _validate_scalar("tolerance", tolerance)
    if tolerance < 0:
        raise ValueError("tolerance must be nonnegative")
    if (
        isinstance(points_per_axis, bool)
        or not isinstance(points_per_axis, int)
        or points_per_axis < 3
    ):
        raise ValueError("points_per_axis must be an integer >= 3")
    if probe_half_width <= support_half_width:
        raise ValueError("probe_half_width must exceed support_half_width")

    faces = _face_points(support_half_width, points_per_axis)
    shell = _exterior_shell_points(
        support_half_width, probe_half_width, points_per_axis
    )
    face_speed = np.linalg.norm(_eval_velocity(velocity, faces, time), axis=1)
    shell_speed = np.linalg.norm(_eval_velocity(velocity, shell, time), axis=1)

    face_max = float(np.max(face_speed))
    shell_max = float(np.max(shell_speed))
    return BoundarySupportReport(
        time=time,
        support_half_width=support_half_width,
        probe_half_width=probe_half_width,
        points_per_axis=points_per_axis,
        face_max_speed=face_max,
        exterior_shell_max_speed=shell_max,
        face_rms_speed=float(np.sqrt(np.mean(face_speed**2))),
        exterior_shell_rms_speed=float(np.sqrt(np.mean(shell_speed**2))),
        sampled_support_pass=bool(face_max <= tolerance and shell_max <= tolerance),
        tolerance=tolerance,
    )
