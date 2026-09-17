"""Truth-bounded arc-length resampling for dense streamline visualization.

This module changes only the polyline representation used for rendering. It does
not modify the underlying velocity candidate, pressure, forcing, PDE residual,
validation thresholds, or scientific claim state.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Iterable

import numpy as np
from scipy.interpolate import PchipInterpolator

VelocityCallable = Callable[[np.ndarray, float], np.ndarray]


@dataclass(frozen=True)
class ResampledStreamlineBundle:
    points: np.ndarray
    scalars: np.ndarray
    lengths: np.ndarray
    segment_length_cv: np.ndarray
    points_per_line: int
    color_mode: str
    time: float | None
    claim_scope: str = "visualization_polyline_resampling_only"
    visualization_ready: bool = False
    visual_correspondence_verified: bool = False
    pde_validated: bool = False
    paper_exact: bool = False
    openai_field_identified: bool = False
    blowup_proved: bool = False


def _clean_polyline(polyline: np.ndarray, *, min_length: float) -> tuple[np.ndarray, np.ndarray, float]:
    pts = np.asarray(polyline, dtype=float)
    if pts.ndim != 2 or pts.shape[1] != 3 or pts.shape[0] < 2:
        raise ValueError("each streamline must have shape (n, 3) with n >= 2")
    if not np.all(np.isfinite(pts)):
        raise ValueError("streamline points must be finite")

    step = np.linalg.norm(np.diff(pts, axis=0), axis=1)
    keep = np.concatenate(([True], step > 0.0))
    pts = pts[keep]
    if pts.shape[0] < 2:
        raise ValueError("streamline collapses to one unique point")

    step = np.linalg.norm(np.diff(pts, axis=0), axis=1)
    s = np.concatenate(([0.0], np.cumsum(step)))
    length = float(s[-1])
    if not np.isfinite(length) or length <= min_length:
        raise ValueError("streamline arclength is too small for visualization resampling")
    return pts, s, length


def _evaluate_velocity(
    velocity: VelocityCallable,
    points: np.ndarray,
    time: float,
) -> np.ndarray:
    values = np.asarray(velocity(points, float(time)), dtype=float)
    if values.shape != points.shape:
        raise ValueError("velocity must return shape (n, 3)")
    if not np.all(np.isfinite(values)):
        raise ValueError("velocity returned non-finite values")
    return values


def resample_streamlines_for_rendering(
    streamlines: Iterable[np.ndarray],
    *,
    points_per_line: int = 256,
    color_mode: str = "height",
    velocity: VelocityCallable | None = None,
    time: float | None = None,
    min_length: float = 1.0e-10,
) -> ResampledStreamlineBundle:
    """Resample streamlines on uniform original-polyline arclength coordinates.

    ``PchipInterpolator`` is used independently on x/y/z as a smooth,
    shape-preserving cubic interpolator. Every output line has the same point
    count, which keeps dense 3D rendering bounded to one graphics object per
    streamline instead of one object per small segment.

    ``color_mode='height'`` uses z directly. ``color_mode='speed'`` evaluates
    the caller-supplied frozen velocity only at the resampled points.
    Resampling is a visualization operation and is never PDE validation.
    """

    if not isinstance(points_per_line, (int, np.integer)) or points_per_line < 2:
        raise ValueError("points_per_line must be an integer >= 2")
    if color_mode not in {"height", "speed"}:
        raise ValueError("color_mode must be 'height' or 'speed'")
    if not np.isfinite(min_length) or min_length <= 0:
        raise ValueError("min_length must be positive and finite")
    if color_mode == "speed":
        if velocity is None:
            raise ValueError("velocity is required for speed coloring")
        if time is None or not np.isfinite(time):
            raise ValueError("finite time is required for speed coloring")
    elif time is not None and not np.isfinite(time):
        raise ValueError("time must be finite when provided")

    raw = list(streamlines)
    if not raw:
        raise ValueError("at least one streamline is required")

    lines: list[np.ndarray] = []
    lengths: list[float] = []
    cvs: list[float] = []

    for polyline in raw:
        pts, s, length = _clean_polyline(polyline, min_length=float(min_length))
        target_s = np.linspace(0.0, length, int(points_per_line), dtype=float)
        interp = PchipInterpolator(s, pts, axis=0, extrapolate=False)
        out = np.asarray(interp(target_s), dtype=float)
        if out.shape != (int(points_per_line), 3) or not np.all(np.isfinite(out)):
            raise ValueError("interpolator returned malformed/non-finite streamline")

        # Preserve the supplied geometry endpoints exactly.
        out[0] = pts[0]
        out[-1] = pts[-1]

        seg = np.linalg.norm(np.diff(out, axis=0), axis=1)
        mean_seg = float(np.mean(seg))
        if not np.isfinite(mean_seg) or mean_seg <= 0:
            raise ValueError("resampled streamline has degenerate segments")
        cv = float(np.std(seg) / mean_seg)

        lines.append(out)
        lengths.append(length)
        cvs.append(cv)

    points = np.stack(lines, axis=0)
    if color_mode == "height":
        scalars = np.array(points[..., 2], copy=True)
    else:
        flat = points.reshape(-1, 3)
        values = _evaluate_velocity(velocity, flat, float(time))  # type: ignore[arg-type]
        scalars = np.linalg.norm(values, axis=1).reshape(points.shape[:2])

    lengths_arr = np.asarray(lengths, dtype=float)
    cvs_arr = np.asarray(cvs, dtype=float)
    if not np.all(np.isfinite(scalars)):
        raise ValueError("visualization scalar contains non-finite values")

    points.setflags(write=False)
    scalars.setflags(write=False)
    lengths_arr.setflags(write=False)
    cvs_arr.setflags(write=False)

    return ResampledStreamlineBundle(
        points=points,
        scalars=scalars,
        lengths=lengths_arr,
        segment_length_cv=cvs_arr,
        points_per_line=int(points_per_line),
        color_mode=color_mode,
        time=None if time is None else float(time),
    )
