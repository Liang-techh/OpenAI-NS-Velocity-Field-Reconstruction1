"""Truth-bounded meridional projection envelope for candidate observables.

This module does not inspect or modify velocity values.  It takes 3-D points
that another diagnostic has already identified as a public-facing candidate
observable (for example, a vorticity-core point cloud or a streamline outer
cloud), projects them onto one declared meridional plane, and returns the
outer convex envelope.

The projection is deliberately restricted to x-z or y-z.  There is no camera
fit, rotation fit, reflection, anisotropic scale fit, threshold selection, or
point correspondence fit here.  The result is visualization geometry only;
it is not PDE validation and it does not identify OpenAI's hidden field.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.spatial import ConvexHull, QhullError


_PROJECTION_AXES = {
    "xz": (0, 2),
    "yz": (1, 2),
}


@dataclass(frozen=True)
class CandidateProjectionEnvelope:
    """Deterministic outer envelope of a declared candidate point cloud."""

    candidate_id: str
    observable: str
    projection: str
    source_point_count: int
    unique_projected_point_count: int
    boundary_points: np.ndarray
    projected_area: float
    projected_perimeter: float
    claim_scope: str = "candidate_visualization_geometry_only"
    camera_fitted: bool = False
    rotation_fitted: bool = False
    reflection_fitted: bool = False
    anisotropic_scale_fitted: bool = False
    velocity_changed: bool = False
    visualization_ready: bool = False
    visual_correspondence_verified: bool = False
    pde_validated: bool = False
    paper_exact: bool = False
    openai_field_identified: bool = False
    blowup_proved: bool = False


def _nonempty_text(value: str, name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{name} must be a non-empty string")
    return value.strip()


def _canonicalize_boundary(points: np.ndarray) -> np.ndarray:
    """Return a deterministic CCW polygon starting at the lexicographic minimum."""

    x = points[:, 0]
    y = points[:, 1]
    signed_twice_area = float(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))
    ordered = points if signed_twice_area > 0.0 else points[::-1]

    start = int(np.lexsort((ordered[:, 1], ordered[:, 0]))[0])
    ordered = np.roll(ordered, -start, axis=0).copy()
    ordered.setflags(write=False)
    return ordered


def project_candidate_meridional_envelope(
    points_xyz: np.ndarray,
    *,
    candidate_id: str,
    observable: str,
    projection: str = "xz",
) -> CandidateProjectionEnvelope:
    """Project a candidate 3-D observable and return its outer 2-D envelope.

    Parameters
    ----------
    points_xyz:
        ``(N, 3)`` Cartesian points already selected by a separate candidate
        diagnostic.  This function does not choose a vorticity/speed threshold.
    candidate_id:
        Frozen candidate identity carried into the result.
    observable:
        Human-readable provenance for the supplied point cloud, such as
        ``"vorticity_core_superlevel"`` or ``"streamline_outer_cloud"``.
    projection:
        Exactly ``"xz"`` or ``"yz"``.  The second output coordinate is always
        physical ``z``; no camera fitting is performed.

    Notes
    -----
    The convex hull is intentionally an *outer-envelope* diagnostic.  It will
    fill concavities and therefore must not be presented as a detailed vortex
    boundary.  Use it only for coarse public silhouette/aspect comparisons.
    """

    candidate_id = _nonempty_text(candidate_id, "candidate_id")
    observable = _nonempty_text(observable, "observable")
    if projection not in _PROJECTION_AXES:
        raise ValueError("projection must be exactly 'xz' or 'yz'")

    points = np.asarray(points_xyz, dtype=float)
    if points.ndim != 2 or points.shape[1] != 3:
        raise ValueError("points_xyz must have shape (N, 3)")
    if points.shape[0] < 3:
        raise ValueError("at least three 3-D points are required")
    if not np.all(np.isfinite(points)):
        raise ValueError("points_xyz must contain only finite values")

    projected = points[:, _PROJECTION_AXES[projection]]
    unique = np.unique(projected, axis=0)
    if unique.shape[0] < 3:
        raise ValueError("projection must contain at least three unique points")

    centered = unique - np.mean(unique, axis=0, keepdims=True)
    if np.linalg.matrix_rank(centered) < 2:
        raise ValueError("projected points are collinear or otherwise degenerate")

    try:
        hull = ConvexHull(unique)
    except QhullError as exc:
        raise ValueError("projected points do not define a stable 2-D convex hull") from exc

    boundary = _canonicalize_boundary(unique[hull.vertices])
    next_boundary = np.roll(boundary, -1, axis=0)
    projected_perimeter = float(np.sum(np.linalg.norm(next_boundary - boundary, axis=1)))
    x = boundary[:, 0]
    z = boundary[:, 1]
    projected_area = 0.5 * abs(float(np.dot(x, np.roll(z, -1)) - np.dot(z, np.roll(x, -1))))

    if not np.isfinite(projected_area) or projected_area <= 0.0:
        raise ValueError("projected hull area must be finite and positive")
    if not np.isfinite(projected_perimeter) or projected_perimeter <= 0.0:
        raise ValueError("projected hull perimeter must be finite and positive")

    return CandidateProjectionEnvelope(
        candidate_id=candidate_id,
        observable=observable,
        projection=projection,
        source_point_count=int(points.shape[0]),
        unique_projected_point_count=int(unique.shape[0]),
        boundary_points=boundary,
        projected_area=projected_area,
        projected_perimeter=projected_perimeter,
    )
