"""Truth-bounded comparison of public 2D projection geometry.

This module compares two already-extracted point clouds in one declared image
frame.  It deliberately does not infer a camera, rotate, mirror, or fit the
velocity field.  The only nuisance normalization is independent translation
and isotropic scale, so axial/transverse aspect and orientation remain visible.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np
from scipy.spatial import cKDTree


CLAIM_SCOPE = "public_observable_projection_geometry_only"


@dataclass(frozen=True)
class ProjectionShapeComparison:
    candidate_count: int
    reference_count: int
    candidate_center: tuple[float, float]
    reference_center: tuple[float, float]
    candidate_scale: float
    reference_scale: float
    candidate_axial_to_transverse_aspect: float
    reference_axial_to_transverse_aspect: float
    log_aspect_error: float
    symmetric_chamfer_mean: float
    symmetric_chamfer_rms: float
    symmetric_hausdorff: float
    candidate_to_reference_q90: float
    reference_to_candidate_q90: float
    alignment: str = "bbox_center_plus_isotropic_bbox_diagonal_scale_only"
    point_convention: str = "columns=(transverse_image_coordinate, axial_image_coordinate)"
    claim_scope: str = CLAIM_SCOPE
    velocity_changed: bool = False
    visualization_ready: bool = False
    visual_correspondence_verified: bool = False
    pde_validated: bool = False
    paper_exact: bool = False
    openai_field_identified: bool = False
    blowup_proved: bool = False


def _validated_points(name: str, points: object, min_points: int) -> np.ndarray:
    array = np.asarray(points, dtype=float)
    if array.ndim != 2 or array.shape[1] != 2:
        raise ValueError(f"{name} must have shape (n, 2)")
    if array.shape[0] < min_points:
        raise ValueError(f"{name} must contain at least {min_points} points")
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must contain only finite values")

    unique = np.unique(array, axis=0)
    if unique.shape[0] < min_points:
        raise ValueError(f"{name} must contain at least {min_points} unique points")

    spans = np.ptp(array, axis=0)
    if np.any(spans <= np.finfo(float).eps):
        raise ValueError(f"{name} must have nonzero transverse and axial extent")
    return array


def _normalize(points: np.ndarray) -> tuple[np.ndarray, tuple[float, float], float, float]:
    lo = np.min(points, axis=0)
    hi = np.max(points, axis=0)
    center = 0.5 * (lo + hi)
    spans = hi - lo
    scale = float(np.linalg.norm(spans))
    if not math.isfinite(scale) or scale <= np.finfo(float).eps:
        raise ValueError("point-cloud normalization scale must be positive and finite")
    normalized = (points - center) / scale
    aspect = float(spans[1] / spans[0])
    return normalized, (float(center[0]), float(center[1])), scale, aspect


def compare_public_projection_shapes(
    candidate_points: object,
    reference_points: object,
    *,
    min_points: int = 8,
) -> ProjectionShapeComparison:
    """Compare two visible 2D point clouds without fitting hidden geometry.

    Parameters
    ----------
    candidate_points, reference_points:
        Arrays of shape ``(n, 2)`` in a common image-axis convention.  Column
        zero is transverse and column one is axial.  The caller is responsible
        for extracting these points from a candidate rendering and a public
        reference image, respectively.
    min_points:
        Minimum number of unique points required in each cloud.

    Notes
    -----
    Each cloud is centered by its own bounding-box midpoint and divided by its
    own bounding-box diagonal.  No rotation, reflection, anisotropic scaling,
    camera optimization, or point correspondence is fitted.  Nearest-neighbor
    distances are computed with SciPy's public ``cKDTree`` API.  The returned
    numbers are scale-free visualization diagnostics only; no pass threshold is
    defined here.
    """

    if isinstance(min_points, bool) or not isinstance(min_points, (int, np.integer)):
        raise TypeError("min_points must be an integer")
    if int(min_points) < 4:
        raise ValueError("min_points must be at least 4")
    min_points = int(min_points)

    candidate = _validated_points("candidate_points", candidate_points, min_points)
    reference = _validated_points("reference_points", reference_points, min_points)

    cand_norm, cand_center, cand_scale, cand_aspect = _normalize(candidate)
    ref_norm, ref_center, ref_scale, ref_aspect = _normalize(reference)

    reference_tree = cKDTree(ref_norm, copy_data=True)
    candidate_tree = cKDTree(cand_norm, copy_data=True)
    cand_to_ref = np.asarray(reference_tree.query(cand_norm, k=1, workers=1)[0], dtype=float)
    ref_to_cand = np.asarray(candidate_tree.query(ref_norm, k=1, workers=1)[0], dtype=float)

    if not (np.all(np.isfinite(cand_to_ref)) and np.all(np.isfinite(ref_to_cand))):
        raise ValueError("nearest-neighbor distances must be finite")

    chamfer_mean = 0.5 * (float(np.mean(cand_to_ref)) + float(np.mean(ref_to_cand)))
    chamfer_rms = math.sqrt(
        0.5 * (float(np.mean(cand_to_ref**2)) + float(np.mean(ref_to_cand**2)))
    )
    hausdorff = max(float(np.max(cand_to_ref)), float(np.max(ref_to_cand)))
    log_aspect_error = abs(math.log(cand_aspect / ref_aspect))

    return ProjectionShapeComparison(
        candidate_count=int(candidate.shape[0]),
        reference_count=int(reference.shape[0]),
        candidate_center=cand_center,
        reference_center=ref_center,
        candidate_scale=cand_scale,
        reference_scale=ref_scale,
        candidate_axial_to_transverse_aspect=cand_aspect,
        reference_axial_to_transverse_aspect=ref_aspect,
        log_aspect_error=float(log_aspect_error),
        symmetric_chamfer_mean=float(chamfer_mean),
        symmetric_chamfer_rms=float(chamfer_rms),
        symmetric_hausdorff=float(hausdorff),
        candidate_to_reference_q90=float(np.quantile(cand_to_ref, 0.9)),
        reference_to_candidate_q90=float(np.quantile(ref_to_cand, 0.9)),
    )
