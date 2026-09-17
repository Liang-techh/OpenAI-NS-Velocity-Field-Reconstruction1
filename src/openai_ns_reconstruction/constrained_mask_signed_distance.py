"""Truth-bounded signed-distance residuals for same-frame morphology fitting.

This module converts two explicitly presegmented 2-D masks in one fixed
projection/pixel frame into a dimensionless signed-distance residual vector.
The residual is a visual morphology fitting signal only. It does not segment
images, fit a camera, align shapes, change velocity values, or validate the
Navier--Stokes equations.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Any

import numpy as np
from scipy.ndimage import distance_transform_edt


@dataclass(frozen=True)
class SignedDistanceVisualResidual:
    """Same-frame visual residual and governance metadata."""

    signed_distance_candidate: np.ndarray
    signed_distance_reference: np.ndarray
    normalized_delta: np.ndarray
    weighted_residual_vector: np.ndarray
    foreground_rms: float
    background_rms: float
    balanced_rms: float
    max_abs: float
    pixel_spacing: tuple[float, float]
    frame_diagonal: float
    candidate_provenance: str
    reference_provenance: str
    truth_boundary: dict[str, Any]


def _strict_binary_mask(mask: Any, name: str) -> np.ndarray:
    array = np.asarray(mask)
    if array.ndim != 2 or min(array.shape) < 3:
        raise ValueError(f"{name} must be a 2-D mask with each dimension >= 3")
    try:
        finite = np.isfinite(array)
    except TypeError as exc:
        raise ValueError(f"{name} must contain numeric or boolean values") from exc
    if not np.all(finite):
        raise ValueError(f"{name} must contain only finite values")
    if array.dtype == np.bool_:
        binary = array.copy()
    else:
        values = np.unique(array)
        if not np.all(np.isin(values, (0, 1))):
            raise ValueError(f"{name} must be boolean or contain exact 0/1 values")
        binary = array.astype(bool, copy=True)
    if not np.any(binary):
        raise ValueError(f"{name} must contain foreground pixels")
    if np.all(binary):
        raise ValueError(f"{name} must contain background pixels")
    return binary


def _spacing2(pixel_spacing: Any) -> tuple[float, float]:
    spacing = np.asarray(pixel_spacing, dtype=float)
    if spacing.shape != (2,) or not np.all(np.isfinite(spacing)) or np.any(spacing <= 0):
        raise ValueError("pixel_spacing must contain two finite positive values")
    return float(spacing[0]), float(spacing[1])


def _provenance(value: Any, name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{name} must be a non-empty string")
    return value.strip()


def _signed_distance(mask: np.ndarray, spacing: tuple[float, float]) -> np.ndarray:
    inside = distance_transform_edt(mask, sampling=spacing)
    outside = distance_transform_edt(~mask, sampling=spacing)
    return inside - outside


def build_signed_distance_visual_residual(
    candidate_mask: Any,
    reference_mask: Any,
    *,
    pixel_spacing: tuple[float, float] = (1.0, 1.0),
    candidate_provenance: str,
    reference_provenance: str,
) -> SignedDistanceVisualResidual:
    """Build a fixed-frame signed-distance morphology residual.

    Both masks must already be explicitly segmented and must use the same pixel
    frame. No translation, rotation, reflection, scaling, camera fitting,
    threshold selection, or time alignment is performed.

    The returned residual vector uses reference foreground/background balancing:
    each class contributes one half of the squared objective. This prevents a
    large background from numerically erasing a smaller visible structure while
    keeping the residual deterministic and free of a fitted tolerance.
    """
    candidate = _strict_binary_mask(candidate_mask, "candidate_mask")
    reference = _strict_binary_mask(reference_mask, "reference_mask")
    if candidate.shape != reference.shape:
        raise ValueError("candidate_mask and reference_mask must have the same shape")

    spacing = _spacing2(pixel_spacing)
    candidate_source = _provenance(candidate_provenance, "candidate_provenance")
    reference_source = _provenance(reference_provenance, "reference_provenance")

    height, width = candidate.shape
    frame_diagonal = math.hypot((height - 1) * spacing[0], (width - 1) * spacing[1])
    if not math.isfinite(frame_diagonal) or frame_diagonal <= 0:
        raise ValueError("pixel frame diagonal must be finite and positive")

    candidate_sdf = _signed_distance(candidate, spacing)
    reference_sdf = _signed_distance(reference, spacing)
    delta = (candidate_sdf - reference_sdf) / frame_diagonal

    foreground = reference
    background = ~reference
    fg_sq = np.square(delta[foreground])
    bg_sq = np.square(delta[background])
    foreground_rms = float(np.sqrt(np.mean(fg_sq)))
    background_rms = float(np.sqrt(np.mean(bg_sq)))
    balanced_mse = 0.5 * (float(np.mean(fg_sq)) + float(np.mean(bg_sq)))
    balanced_rms = math.sqrt(balanced_mse)
    max_abs = float(np.max(np.abs(delta)))

    weights = np.empty(delta.shape, dtype=float)
    weights[foreground] = 0.5 / int(np.count_nonzero(foreground))
    weights[background] = 0.5 / int(np.count_nonzero(background))
    weighted = (delta * np.sqrt(weights)).ravel()

    for array in (candidate_sdf, reference_sdf, delta, weighted):
        array.setflags(write=False)

    truth_boundary = {
        "purpose": "same_frame_public_observable_morphology_fit_only",
        "segmentation_performed": False,
        "alignment_performed": False,
        "camera_fitted": False,
        "time_alignment_performed": False,
        "velocity_changed": False,
        "pde_threshold_changed": False,
        "defines_visual_acceptance_threshold": False,
        "visualization_ready": False,
        "visual_correspondence_verified": False,
        "pde_validated": False,
        "paper_exact": False,
        "openai_field_identified": False,
        "blowup_proved": False,
    }

    return SignedDistanceVisualResidual(
        signed_distance_candidate=candidate_sdf,
        signed_distance_reference=reference_sdf,
        normalized_delta=delta,
        weighted_residual_vector=weighted,
        foreground_rms=foreground_rms,
        background_rms=background_rms,
        balanced_rms=balanced_rms,
        max_abs=max_abs,
        pixel_spacing=spacing,
        frame_diagonal=frame_diagonal,
        candidate_provenance=candidate_source,
        reference_provenance=reference_source,
        truth_boundary=truth_boundary,
    )
