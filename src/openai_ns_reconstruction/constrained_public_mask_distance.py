"""Truth-bounded same-frame comparison of presegmented public/candidate masks.

This module compares two already-segmented 2-D visible-region masks in one
fixed image frame.  It deliberately performs no segmentation, registration,
camera fit, rotation, reflection, or scale fit.  Surface distances use
SciPy's public Euclidean distance-transform API.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Iterable

import numpy as np
from scipy.ndimage import distance_transform_edt


CLAIM_SCOPE = "public_observable_same_frame_mask_geometry_only"


@dataclass(frozen=True)
class PublicMaskDistance:
    shape: tuple[int, int]
    pixel_spacing_axial_transverse: tuple[float, float]
    candidate_foreground_count: int
    reference_foreground_count: int
    candidate_area_fraction: float
    reference_area_fraction: float
    area_ratio_candidate_to_reference: float
    intersection_over_union: float
    dice_coefficient: float
    symmetric_surface_mean: float
    symmetric_surface_rms: float
    symmetric_surface_q90: float
    symmetric_hausdorff: float
    normalized_surface_mean: float
    normalized_surface_rms: float
    normalized_surface_q90: float
    normalized_hausdorff: float
    candidate_touches_image_frame: bool
    reference_touches_image_frame: bool
    alignment: str = "none_same_pixel_frame_required"
    segmentation_performed: bool = False
    threshold_selected: bool = False
    camera_fitted: bool = False
    velocity_changed: bool = False
    visualization_ready: bool = False
    visual_correspondence_verified: bool = False
    pde_validated: bool = False
    paper_exact: bool = False
    openai_field_identified: bool = False
    blowup_proved: bool = False
    claim_scope: str = CLAIM_SCOPE


def _validated_mask(name: str, mask: object) -> np.ndarray:
    array = np.asarray(mask)
    if array.ndim != 2:
        raise ValueError(f"{name} must be a 2-D mask")
    if min(array.shape) < 3:
        raise ValueError(f"{name} must be at least 3 pixels along each axis")

    if array.dtype == np.bool_:
        result = np.array(array, dtype=bool, copy=True)
    elif np.issubdtype(array.dtype, np.number):
        numeric = np.asarray(array, dtype=float)
        if not np.all(np.isfinite(numeric)):
            raise ValueError(f"{name} must contain only finite values")
        if not np.all((numeric == 0.0) | (numeric == 1.0)):
            raise ValueError(f"{name} must be Boolean or exact 0/1; grayscale thresholding is not allowed")
        result = numeric.astype(bool)
    else:
        raise TypeError(f"{name} must be Boolean or numeric exact 0/1")

    foreground = int(np.count_nonzero(result))
    if foreground == 0:
        raise ValueError(f"{name} must contain a nonempty foreground")
    if foreground == result.size:
        raise ValueError(f"{name} cannot be full-frame foreground")
    return result


def _validated_spacing(pixel_spacing: Iterable[float]) -> tuple[float, float]:
    values = tuple(pixel_spacing)
    if len(values) != 2:
        raise ValueError("pixel_spacing must contain exactly (axial, transverse) spacing")
    try:
        axial, transverse = (float(values[0]), float(values[1]))
    except (TypeError, ValueError) as exc:
        raise TypeError("pixel_spacing entries must be real numbers") from exc
    if not (math.isfinite(axial) and math.isfinite(transverse)):
        raise ValueError("pixel_spacing entries must be finite")
    if axial <= 0.0 or transverse <= 0.0:
        raise ValueError("pixel_spacing entries must be strictly positive")
    return axial, transverse


def _boundary4(mask: np.ndarray) -> np.ndarray:
    padded = np.pad(mask, 1, mode="constant", constant_values=False)
    center = padded[1:-1, 1:-1]
    interior = (
        center
        & padded[:-2, 1:-1]
        & padded[2:, 1:-1]
        & padded[1:-1, :-2]
        & padded[1:-1, 2:]
    )
    boundary = center & ~interior
    if not np.any(boundary):
        raise ValueError("mask must have a nonempty 4-connected boundary")
    return boundary


def _touches_frame(mask: np.ndarray) -> bool:
    return bool(
        np.any(mask[0, :])
        or np.any(mask[-1, :])
        or np.any(mask[:, 0])
        or np.any(mask[:, -1])
    )


def compare_presegmented_masks(
    candidate_mask: object,
    reference_mask: object,
    *,
    pixel_spacing: Iterable[float] = (1.0, 1.0),
) -> PublicMaskDistance:
    """Compare two presegmented visible-region masks in exactly one image frame.

    Parameters
    ----------
    candidate_mask, reference_mask:
        Boolean or exact 0/1 arrays with identical shape.  The caller is
        responsible for segmentation and for placing both masks in the same
        declared projection/pixel frame.  This function performs no alignment.
    pixel_spacing:
        Positive ``(axial, transverse)`` pixel spacing.  The default reports
        distances in pixels.  No spacing is inferred from the images.

    Notes
    -----
    IoU/Dice compare occupied image area.  Surface distances are bidirectional
    distances between 4-connected foreground boundaries, evaluated through
    SciPy's ``distance_transform_edt``.  Distances are also divided by the
    physical image diagonal for a dimensionless report.  There is deliberately
    no pass threshold and no promotion of visualization/PDE readiness.
    """

    candidate = _validated_mask("candidate_mask", candidate_mask)
    reference = _validated_mask("reference_mask", reference_mask)
    if candidate.shape != reference.shape:
        raise ValueError("candidate_mask and reference_mask must have identical shape")
    spacing = _validated_spacing(pixel_spacing)

    candidate_boundary = _boundary4(candidate)
    reference_boundary = _boundary4(reference)

    candidate_to_reference_field = distance_transform_edt(
        ~reference_boundary, sampling=spacing
    )
    reference_to_candidate_field = distance_transform_edt(
        ~candidate_boundary, sampling=spacing
    )
    candidate_to_reference = np.asarray(
        candidate_to_reference_field[candidate_boundary], dtype=float
    )
    reference_to_candidate = np.asarray(
        reference_to_candidate_field[reference_boundary], dtype=float
    )
    if not (
        np.all(np.isfinite(candidate_to_reference))
        and np.all(np.isfinite(reference_to_candidate))
    ):
        raise ValueError("surface distances must be finite")

    bidirectional = np.concatenate((candidate_to_reference, reference_to_candidate))
    surface_mean = 0.5 * (
        float(np.mean(candidate_to_reference))
        + float(np.mean(reference_to_candidate))
    )
    surface_rms = math.sqrt(
        0.5
        * (
            float(np.mean(candidate_to_reference**2))
            + float(np.mean(reference_to_candidate**2))
        )
    )
    surface_q90 = float(np.quantile(bidirectional, 0.9))
    hausdorff = max(
        float(np.max(candidate_to_reference)),
        float(np.max(reference_to_candidate)),
    )

    intersection = int(np.count_nonzero(candidate & reference))
    union = int(np.count_nonzero(candidate | reference))
    candidate_count = int(np.count_nonzero(candidate))
    reference_count = int(np.count_nonzero(reference))
    if union == 0 or reference_count == 0:
        raise ValueError("mask overlap statistics are undefined")
    iou = float(intersection / union)
    dice = float(2.0 * intersection / (candidate_count + reference_count))

    rows, cols = candidate.shape
    image_diagonal = math.hypot((rows - 1) * spacing[0], (cols - 1) * spacing[1])
    if not math.isfinite(image_diagonal) or image_diagonal <= 0.0:
        raise ValueError("image diagonal must be positive and finite")

    return PublicMaskDistance(
        shape=(int(rows), int(cols)),
        pixel_spacing_axial_transverse=spacing,
        candidate_foreground_count=candidate_count,
        reference_foreground_count=reference_count,
        candidate_area_fraction=float(candidate_count / candidate.size),
        reference_area_fraction=float(reference_count / reference.size),
        area_ratio_candidate_to_reference=float(candidate_count / reference_count),
        intersection_over_union=iou,
        dice_coefficient=dice,
        symmetric_surface_mean=float(surface_mean),
        symmetric_surface_rms=float(surface_rms),
        symmetric_surface_q90=surface_q90,
        symmetric_hausdorff=float(hausdorff),
        normalized_surface_mean=float(surface_mean / image_diagonal),
        normalized_surface_rms=float(surface_rms / image_diagonal),
        normalized_surface_q90=float(surface_q90 / image_diagonal),
        normalized_hausdorff=float(hausdorff / image_diagonal),
        candidate_touches_image_frame=_touches_frame(candidate),
        reference_touches_image_frame=_touches_frame(reference),
    )
