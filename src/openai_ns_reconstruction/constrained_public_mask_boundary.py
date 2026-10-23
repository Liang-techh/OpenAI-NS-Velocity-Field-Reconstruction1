"""Truth-bounded boundary extraction from an already-segmented public image mask.

This module does not segment, threshold, infer camera geometry, or recover hidden
velocity values.  The caller supplies a binary mask representing a visible
public-image region.  We expose its one-pixel 4-connected boundary as a 2D point
cloud that can be passed to the public projection-shape comparator.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import ndimage


CLAIM_SCOPE = "public_observable_presegmented_mask_boundary_only"


@dataclass(frozen=True)
class PublicMaskBoundary:
    points: np.ndarray
    image_shape: tuple[int, int]
    foreground_pixels: int
    boundary_pixels: int
    foreground_fraction: float
    touches_image_frame: bool
    connectivity: int = 4
    point_convention: str = "columns=(transverse_pixel, axial_pixel_up); pixel_centers"
    claim_scope: str = CLAIM_SCOPE
    segmentation_performed: bool = False
    camera_fitted: bool = False
    velocity_changed: bool = False
    visualization_ready: bool = False
    visual_correspondence_verified: bool = False
    pde_validated: bool = False
    paper_exact: bool = False
    openai_field_identified: bool = False
    blowup_proved: bool = False


def _validated_binary_mask(mask: object) -> np.ndarray:
    array = np.asarray(mask)
    if array.ndim != 2:
        raise ValueError("mask must be a two-dimensional array")
    if min(array.shape) < 3:
        raise ValueError("mask dimensions must each be at least 3 pixels")
    if np.issubdtype(array.dtype, np.complexfloating):
        raise TypeError("mask must not be complex")

    if array.dtype == np.bool_:
        binary = array.copy()
    else:
        if not np.issubdtype(array.dtype, np.number):
            raise TypeError("mask must be boolean or numeric 0/1")
        values = np.asarray(array, dtype=float)
        if not np.all(np.isfinite(values)):
            raise ValueError("mask must contain only finite values")
        if not np.all((values == 0.0) | (values == 1.0)):
            raise ValueError(
                "mask must already be segmented as exact 0/1 values; "
                "grayscale thresholding is intentionally out of scope"
            )
        binary = values.astype(bool)

    if not np.any(binary):
        raise ValueError("mask must contain foreground pixels")
    if np.all(binary):
        raise ValueError("mask must contain background pixels as well as foreground")
    return binary


def extract_public_mask_boundary(
    mask: object,
    *,
    min_boundary_points: int = 8,
) -> PublicMaskBoundary:
    """Extract a 4-connected one-pixel boundary from a supplied binary mask.

    The mask is treated as an already-made public-observable segmentation.
    No image threshold, contour smoothing, rotation, reflection, camera fit, or
    velocity-field fit is performed here.  Image columns become the transverse
    coordinate.  Image rows are flipped so larger returned axial coordinate is
    upward in the displayed image.

    SciPy's public ``ndimage.binary_erosion`` API supplies the morphology
    primitive.  The boundary is defined locally as ``mask & ~erode(mask)`` with
    the rank-2, connectivity-1 structuring element.
    """
    if isinstance(min_boundary_points, bool) or not isinstance(
        min_boundary_points, (int, np.integer)
    ):
        raise TypeError("min_boundary_points must be an integer")
    min_boundary_points = int(min_boundary_points)
    if min_boundary_points < 4:
        raise ValueError("min_boundary_points must be at least 4")

    binary = _validated_binary_mask(mask)
    structure = ndimage.generate_binary_structure(2, 1)
    eroded = ndimage.binary_erosion(
        binary,
        structure=structure,
        iterations=1,
        border_value=0,
    )
    boundary = binary & ~np.asarray(eroded, dtype=bool)
    rows, cols = np.nonzero(boundary)
    count = int(rows.size)
    if count < min_boundary_points:
        raise ValueError(
            f"mask boundary has {count} pixels, fewer than required {min_boundary_points}"
        )

    height, width = binary.shape
    axial_up = (height - 1) - rows
    points = np.column_stack((cols, axial_up)).astype(float, copy=False)
    points.setflags(write=False)

    touches = bool(
        np.any(binary[0, :])
        or np.any(binary[-1, :])
        or np.any(binary[:, 0])
        or np.any(binary[:, -1])
    )
    foreground_pixels = int(np.count_nonzero(binary))
    return PublicMaskBoundary(
        points=points,
        image_shape=(int(height), int(width)),
        foreground_pixels=foreground_pixels,
        boundary_pixels=count,
        foreground_fraction=float(foreground_pixels / binary.size),
        touches_image_frame=touches,
    )
