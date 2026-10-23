"""Truth-bounded morphology features for presegmented public/candidate masks.

This module measures an already-segmented binary mask in a fixed image frame.
It does not segment imagery, choose a threshold, fit a camera, register frames,
or modify a velocity field.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping

import numpy as np
from scipy.ndimage import binary_fill_holes


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
class MaskMorphologyFeatures:
    foreground_pixels: int
    area_fraction: float
    physical_area: float
    transverse_centroid: float
    axial_centroid: float
    transverse_extent: float
    axial_extent: float
    transverse_rms: float
    axial_rms: float
    bbox_aspect_axial_over_transverse: float
    rms_aspect_axial_over_transverse: float
    hole_pixels: int
    hole_area_fraction_of_filled: float
    touches_image_frame: bool
    pixel_spacing: tuple[float, float]
    mask_id: str
    segmentation_provenance: str
    frame_provenance: str
    truth_boundary: Mapping[str, bool]


def _nonempty_text(value: str, name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{name} must be a non-empty string")
    return value.strip()


def _binary_mask(mask: np.ndarray) -> np.ndarray:
    arr = np.asarray(mask)
    if arr.ndim != 2 or min(arr.shape) < 3:
        raise ValueError("mask must be a 2-D array with both dimensions at least 3")
    if arr.dtype == np.bool_:
        out = arr.astype(bool, copy=False)
    else:
        if not np.issubdtype(arr.dtype, np.number):
            raise ValueError("mask must be boolean or exact numeric 0/1")
        if not np.isfinite(arr).all():
            raise ValueError("mask must contain only finite values")
        if not bool(np.all((arr == 0) | (arr == 1))):
            raise ValueError("mask must be boolean or exact numeric 0/1")
        out = arr.astype(bool)
    if not bool(np.any(out)):
        raise ValueError("mask must contain foreground pixels")
    if bool(np.all(out)):
        raise ValueError("mask may not fill the entire frame")
    return out


def _spacing(pixel_spacing: tuple[float, float]) -> tuple[float, float]:
    try:
        axial, transverse = pixel_spacing
    except Exception as exc:  # pragma: no cover - defensive shape guard
        raise ValueError(
            "pixel_spacing must contain exactly (axial_spacing, transverse_spacing)"
        ) from exc
    axial = float(axial)
    transverse = float(transverse)
    if not np.isfinite([axial, transverse]).all() or axial <= 0.0 or transverse <= 0.0:
        raise ValueError("pixel_spacing values must be finite and strictly positive")
    return axial, transverse


def measure_presegmented_mask_morphology(
    mask: np.ndarray,
    *,
    pixel_spacing: tuple[float, float] = (1.0, 1.0),
    mask_id: str,
    segmentation_provenance: str,
    frame_provenance: str,
) -> MaskMorphologyFeatures:
    """Measure axis-fixed morphology of an explicitly presegmented mask.

    Array columns are the transverse axis. Array rows are the axial axis, with
    row zero treated as visually upward. Coordinates are centered on the image
    frame; no translation, rotation, reflection, scaling, registration, or
    camera fit is performed.

    ``pixel_spacing`` is ``(axial_spacing, transverse_spacing)`` and must be
    supplied from the declared comparison frame when physical units are known.
    """
    binary = _binary_mask(mask)
    axial_spacing, transverse_spacing = _spacing(pixel_spacing)
    mask_id = _nonempty_text(mask_id, "mask_id")
    segmentation_provenance = _nonempty_text(
        segmentation_provenance, "segmentation_provenance"
    )
    frame_provenance = _nonempty_text(frame_provenance, "frame_provenance")

    rows, cols = binary.shape
    rr, cc = np.nonzero(binary)

    transverse_coords = (cc - 0.5 * (cols - 1)) * transverse_spacing
    axial_coords = (0.5 * (rows - 1) - rr) * axial_spacing

    transverse_centroid = float(np.mean(transverse_coords))
    axial_centroid = float(np.mean(axial_coords))
    transverse_rms = float(
        np.sqrt(np.mean((transverse_coords - transverse_centroid) ** 2))
    )
    axial_rms = float(np.sqrt(np.mean((axial_coords - axial_centroid) ** 2)))

    transverse_extent = float(
        (int(np.max(cc)) - int(np.min(cc)) + 1) * transverse_spacing
    )
    axial_extent = float(
        (int(np.max(rr)) - int(np.min(rr)) + 1) * axial_spacing
    )
    if transverse_extent <= 0.0 or axial_extent <= 0.0:
        raise ValueError("mask extents must be strictly positive")
    if transverse_rms <= 0.0 or axial_rms <= 0.0:
        raise ValueError(
            "mask must span at least two distinct pixel centers along both axes"
        )

    filled = binary_fill_holes(binary)
    holes = np.asarray(filled, dtype=bool) & ~binary
    hole_pixels = int(np.count_nonzero(holes))
    filled_pixels = int(np.count_nonzero(filled))
    hole_fraction = float(hole_pixels / filled_pixels)

    touches = bool(
        np.any(binary[0, :])
        or np.any(binary[-1, :])
        or np.any(binary[:, 0])
        or np.any(binary[:, -1])
    )

    foreground_pixels = int(rr.size)
    pixel_area = axial_spacing * transverse_spacing
    return MaskMorphologyFeatures(
        foreground_pixels=foreground_pixels,
        area_fraction=float(foreground_pixels / binary.size),
        physical_area=float(foreground_pixels * pixel_area),
        transverse_centroid=transverse_centroid,
        axial_centroid=axial_centroid,
        transverse_extent=transverse_extent,
        axial_extent=axial_extent,
        transverse_rms=transverse_rms,
        axial_rms=axial_rms,
        bbox_aspect_axial_over_transverse=float(axial_extent / transverse_extent),
        rms_aspect_axial_over_transverse=float(axial_rms / transverse_rms),
        hole_pixels=hole_pixels,
        hole_area_fraction_of_filled=hole_fraction,
        touches_image_frame=touches,
        pixel_spacing=(axial_spacing, transverse_spacing),
        mask_id=mask_id,
        segmentation_provenance=segmentation_provenance,
        frame_provenance=frame_provenance,
        truth_boundary=dict(_TRUTH_BOUNDARY),
    )
