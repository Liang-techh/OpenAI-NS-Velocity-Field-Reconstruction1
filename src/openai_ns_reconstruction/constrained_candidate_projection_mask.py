"""Truth-bounded fixed-frame rasterization for selected candidate observables.

This module does not select velocity/vorticity thresholds and does not fit a
camera or registration. It only converts an already-selected 3-D observable
point cloud into a fixed orthographic x-z or y-z binary occupancy mask.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping

import numpy as np
from scipy.ndimage import binary_dilation


_TRUTH_BOUNDARY = {
    "segmentation_performed": False,
    "threshold_selected": False,
    "camera_fitted": False,
    "registration_fitted": False,
    "velocity_changed": False,
    "visualization_ready": False,
    "visual_correspondence_verified": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


@dataclass(frozen=True)
class CandidateProjectionMask:
    mask: np.ndarray
    raw_occupancy: np.ndarray
    projection: str
    raster_shape: tuple[int, int]
    transverse_bounds: tuple[float, float]
    axial_bounds: tuple[float, float]
    point_radius_pixels: int
    point_count: int
    occupied_pixels_before_dilation: int
    occupied_pixels_after_dilation: int
    touches_image_frame: bool
    candidate_id: str
    point_selection_provenance: str
    frame_provenance: str
    truth_boundary: Mapping[str, bool]


def _nonempty_text(value: str, name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{name} must be a non-empty string")
    return value.strip()


def _bounds(value: tuple[float, float], name: str) -> tuple[float, float]:
    try:
        lo, hi = value
    except Exception as exc:  # pragma: no cover - defensive shape guard
        raise ValueError(f"{name} must contain exactly two values") from exc
    lo = float(lo)
    hi = float(hi)
    if not np.isfinite([lo, hi]).all() or not lo < hi:
        raise ValueError(f"{name} must be finite and strictly increasing")
    return lo, hi


def _shape(value: tuple[int, int]) -> tuple[int, int]:
    if len(value) != 2:
        raise ValueError("shape must be (rows, columns)")
    rows, cols = value
    if isinstance(rows, bool) or isinstance(cols, bool):
        raise ValueError("shape entries must be integers")
    if int(rows) != rows or int(cols) != cols:
        raise ValueError("shape entries must be integers")
    rows = int(rows)
    cols = int(cols)
    if rows < 3 or cols < 3:
        raise ValueError("shape entries must both be at least 3")
    return rows, cols


def _disk(radius: int) -> np.ndarray:
    yy, xx = np.ogrid[-radius : radius + 1, -radius : radius + 1]
    return (xx * xx + yy * yy) <= radius * radius


def rasterize_candidate_projection_points(
    points_xyz: np.ndarray,
    *,
    projection: str,
    shape: tuple[int, int],
    transverse_bounds: tuple[float, float],
    axial_bounds: tuple[float, float],
    point_radius_pixels: int,
    candidate_id: str,
    point_selection_provenance: str,
    frame_provenance: str,
) -> CandidateProjectionMask:
    """Rasterize selected 3-D candidate points in a declared fixed frame.

    ``points_xyz`` must already represent a selected public-observable candidate
    structure (for example a frozen vorticity-core point cloud or sampled
    streamline cloud). This function does not choose that selection threshold.

    The returned image uses ordinary array coordinates: row 0 is the largest
    axial coordinate, and columns increase with the declared transverse axis.
    All points must lie inside the declared frame; silent cropping is rejected.
    ``point_radius_pixels`` is caller-declared and is never optimized here.
    """

    points = np.asarray(points_xyz, dtype=float)
    if points.ndim != 2 or points.shape[1] != 3 or points.shape[0] == 0:
        raise ValueError("points_xyz must have shape (N, 3) with N >= 1")
    if not np.isfinite(points).all():
        raise ValueError("points_xyz must contain only finite values")

    if projection not in {"xz", "yz"}:
        raise ValueError("projection must be exactly 'xz' or 'yz'")
    rows, cols = _shape(shape)
    t_lo, t_hi = _bounds(transverse_bounds, "transverse_bounds")
    a_lo, a_hi = _bounds(axial_bounds, "axial_bounds")

    if isinstance(point_radius_pixels, bool) or int(point_radius_pixels) != point_radius_pixels:
        raise ValueError("point_radius_pixels must be a non-negative integer")
    radius = int(point_radius_pixels)
    if radius < 0:
        raise ValueError("point_radius_pixels must be a non-negative integer")

    candidate_id = _nonempty_text(candidate_id, "candidate_id")
    point_selection_provenance = _nonempty_text(
        point_selection_provenance, "point_selection_provenance"
    )
    frame_provenance = _nonempty_text(frame_provenance, "frame_provenance")

    transverse = points[:, 0] if projection == "xz" else points[:, 1]
    axial = points[:, 2]
    in_frame = (
        (transverse >= t_lo)
        & (transverse <= t_hi)
        & (axial >= a_lo)
        & (axial <= a_hi)
    )
    if not bool(np.all(in_frame)):
        bad = int(np.count_nonzero(~in_frame))
        raise ValueError(
            f"{bad} selected point(s) fall outside the declared fixed projection frame; "
            "silent cropping is forbidden"
        )

    # histogram2d gives axial bins from low to high. Flip the first axis so row
    # zero is visually 'up' (largest axial coordinate), matching image arrays.
    counts, _, _ = np.histogram2d(
        axial,
        transverse,
        bins=(rows, cols),
        range=((a_lo, a_hi), (t_lo, t_hi)),
    )
    raw = counts[::-1, :] > 0

    if radius:
        mask = binary_dilation(raw, structure=_disk(radius), iterations=1)
    else:
        mask = raw.copy()

    if not bool(np.any(mask)):
        raise ValueError("rasterization produced an empty mask")
    if bool(np.all(mask)):
        raise ValueError("rasterization filled the entire frame; comparison would be degenerate")

    frame_touch = bool(
        np.any(mask[0, :])
        or np.any(mask[-1, :])
        or np.any(mask[:, 0])
        or np.any(mask[:, -1])
    )

    raw = np.asarray(raw, dtype=bool)
    mask = np.asarray(mask, dtype=bool)
    raw.setflags(write=False)
    mask.setflags(write=False)

    return CandidateProjectionMask(
        mask=mask,
        raw_occupancy=raw,
        projection=projection,
        raster_shape=(rows, cols),
        transverse_bounds=(t_lo, t_hi),
        axial_bounds=(a_lo, a_hi),
        point_radius_pixels=radius,
        point_count=int(points.shape[0]),
        occupied_pixels_before_dilation=int(np.count_nonzero(raw)),
        occupied_pixels_after_dilation=int(np.count_nonzero(mask)),
        touches_image_frame=frame_touch,
        candidate_id=candidate_id,
        point_selection_provenance=point_selection_provenance,
        frame_provenance=frame_provenance,
        truth_boundary=dict(_TRUTH_BOUNDARY),
    )
