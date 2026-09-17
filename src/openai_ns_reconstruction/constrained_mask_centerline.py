"""Truth-bounded centroidal centerline observables for presegmented masks.

This module measures an already-presegmented 2-D visual mask in a fixed image
frame.  It does not segment imagery, choose thresholds, rotate/translate/scale
frames, fit a camera, infer hidden velocity data, or define a visual acceptance
threshold.  The returned line is a rowwise foreground center-of-mass descriptor,
not a fluid trajectory, vortex line, or recovered hidden centerline.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.ndimage import center_of_mass


@dataclass(frozen=True)
class MaskCenterlineProfile:
    """Fixed-frame row-centroid centerline observables for one binary mask."""

    normalized_axial: NDArray[np.float64]
    axial_position_physical: NDArray[np.float64]
    transverse_center_physical: NDArray[np.float64]
    chord_residual_physical: NDArray[np.float64]
    active_row_count: int
    axial_extent_physical: float
    endpoint_axial_separation_physical: float
    frame_transverse_extent_physical: float
    global_transverse_centroid_physical: float
    bottom_center_physical: float
    midplane_center_physical: float
    top_center_physical: float
    end_to_end_tilt_slope: float
    lateral_drift_physical: float
    lateral_drift_frame_fraction: float
    rms_chord_deviation_physical: float
    max_chord_deviation_physical: float
    rms_chord_deviation_frame_fraction: float
    centerline_arclength_physical: float
    endpoint_chord_length_physical: float
    straightness_ratio: float
    touches_image_frame: bool
    mask_provenance: str
    frame_provenance: str
    segmentation_performed: bool = False
    threshold_selected: bool = False
    camera_fitted: bool = False
    registration_fitted: bool = False
    rotation_fitted: bool = False
    velocity_changed: bool = False
    visualization_ready: bool = False
    visual_correspondence_verified: bool = False
    pde_validated: bool = False
    paper_exact: bool = False
    openai_field_identified: bool = False
    blowup_proved: bool = False


def _binary_mask(mask: ArrayLike) -> NDArray[np.bool_]:
    array = np.asarray(mask)
    if array.ndim != 2 or min(array.shape) < 3:
        raise ValueError("mask must be a 2-D array with at least 3 pixels per axis")
    if array.dtype == np.bool_:
        binary = array.astype(bool, copy=True)
    else:
        try:
            finite = np.isfinite(array)
        except TypeError as exc:
            raise ValueError("mask must contain finite Boolean or exact 0/1 values") from exc
        if not bool(np.all(finite)):
            raise ValueError("mask must contain only finite values")
        if not bool(np.all((array == 0) | (array == 1))):
            raise ValueError("mask must be Boolean or exact 0/1; grayscale masks are not accepted")
        binary = array.astype(bool, copy=True)
    if not bool(np.any(binary)):
        raise ValueError("empty masks are rejected")
    if bool(np.all(binary)):
        raise ValueError("frame-filling masks are rejected")
    return binary


def _positive_spacing(pixel_spacing: tuple[float, float]) -> tuple[float, float]:
    if len(pixel_spacing) != 2:
        raise ValueError("pixel_spacing must be (axial_spacing, transverse_spacing)")
    axial, transverse = (float(pixel_spacing[0]), float(pixel_spacing[1]))
    if not np.isfinite(axial) or not np.isfinite(transverse) or axial <= 0.0 or transverse <= 0.0:
        raise ValueError("pixel spacing values must be finite and strictly positive")
    return axial, transverse


def _sample_count(value: int) -> int:
    if isinstance(value, bool) or not isinstance(value, (int, np.integer)):
        raise ValueError("sample_count must be an odd integer")
    count = int(value)
    if count < 5 or count > 4097 or count % 2 == 0:
        raise ValueError("sample_count must be an odd integer in [5, 4097]")
    return count


def _provenance(value: str, name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{name} must be a non-empty provenance string")
    return value.strip()


def _readonly(array: NDArray[np.float64]) -> NDArray[np.float64]:
    out = np.asarray(array, dtype=np.float64)
    out.setflags(write=False)
    return out


def _row_centers(binary: NDArray[np.bool_], active_rows: NDArray[np.int64]) -> NDArray[np.float64]:
    labels = np.zeros(binary.shape, dtype=np.int32)
    for label_value, row in enumerate(active_rows, start=1):
        labels[row, binary[row]] = label_value
    indices = np.arange(1, active_rows.size + 1, dtype=np.int32)
    centers = center_of_mass(binary.astype(np.float64), labels=labels, index=indices)
    columns = np.asarray([center[1] for center in centers], dtype=np.float64)
    if columns.shape != (active_rows.size,) or not bool(np.all(np.isfinite(columns))):
        raise ValueError("row center-of-mass calculation produced invalid coordinates")
    return columns


def measure_presegmented_mask_centerline(
    mask: ArrayLike,
    *,
    pixel_spacing: tuple[float, float] = (1.0, 1.0),
    sample_count: int = 129,
    mask_provenance: str,
    frame_provenance: str,
) -> MaskCenterlineProfile:
    """Measure a fixed-frame row-centroid centerline from a presegmented mask.

    Foreground-bearing rows must form one contiguous axial interval.  For each
    active row, the foreground column center of mass is measured in the original
    image frame.  Image row order (top -> bottom) is converted to physical
    axial-up order (bottom -> top), then linearly resampled on a fixed descriptor
    grid.  No translation, rotation, reflection, anisotropic scaling, camera fit,
    or registration is performed.

    ``chord_residual_physical`` subtracts only the straight line joining the
    measured bottom/top centers.  The corresponding end-to-end tilt is returned
    separately, so this intrinsic bending diagnostic cannot silently rotate a
    tilted candidate into agreement.  This is a visual morphology descriptor;
    it is not a fluid streamline or hidden-field reconstruction.
    """

    binary = _binary_mask(mask)
    axial_spacing, transverse_spacing = _positive_spacing(pixel_spacing)
    count = _sample_count(sample_count)
    mask_source = _provenance(mask_provenance, "mask_provenance")
    frame_source = _provenance(frame_provenance, "frame_provenance")

    active_rows = np.flatnonzero(np.any(binary, axis=1)).astype(np.int64, copy=False)
    if active_rows.size < 3:
        raise ValueError("at least three active axial rows are required")
    if bool(np.any(np.diff(active_rows) != 1)):
        raise ValueError("foreground rows must form one contiguous axial interval")

    row_center_columns = _row_centers(binary, active_rows)
    row_transverse = (row_center_columns + 0.5) * transverse_spacing
    row_axial = (binary.shape[0] - active_rows.astype(np.float64) - 0.5) * axial_spacing

    raw_transverse = row_transverse[::-1]
    raw_axial = row_axial[::-1]
    endpoint_separation = float(raw_axial[-1] - raw_axial[0])
    if not np.isfinite(endpoint_separation) or endpoint_separation <= 0.0:
        raise ValueError("active axial rows do not define a positive endpoint separation")

    raw_normalized = 2.0 * (raw_axial - raw_axial[0]) / endpoint_separation - 1.0
    sample_axis = np.linspace(-1.0, 1.0, count, dtype=np.float64)
    sample_axial = np.linspace(raw_axial[0], raw_axial[-1], count, dtype=np.float64)
    sampled_transverse = np.interp(sample_axis, raw_normalized, raw_transverse)
    if not bool(np.all(np.isfinite(sampled_transverse))):
        raise ValueError("centerline resampling produced nonfinite coordinates")

    chord = np.linspace(sampled_transverse[0], sampled_transverse[-1], count, dtype=np.float64)
    chord_residual = sampled_transverse - chord
    rms_chord = float(np.sqrt(np.mean(chord_residual * chord_residual)))
    max_chord = float(np.max(np.abs(chord_residual)))

    d_axial = np.diff(sample_axial)
    d_transverse = np.diff(sampled_transverse)
    centerline_arclength = float(np.sum(np.hypot(d_axial, d_transverse)))
    endpoint_chord = float(np.hypot(endpoint_separation, sampled_transverse[-1] - sampled_transverse[0]))
    if endpoint_chord <= 0.0 or not np.isfinite(endpoint_chord):
        raise ValueError("endpoint chord length must be finite and positive")
    straightness_ratio = centerline_arclength / endpoint_chord
    tolerance = 128.0 * np.finfo(float).eps
    if straightness_ratio < 1.0 - tolerance:
        raise ValueError("centerline arclength fell below endpoint chord length")
    straightness_ratio = max(1.0, float(straightness_ratio))

    global_center = center_of_mass(binary.astype(np.float64))
    global_col = float(global_center[1])
    if not np.isfinite(global_col):
        raise ValueError("global center-of-mass calculation produced a nonfinite coordinate")
    global_transverse = (global_col + 0.5) * transverse_spacing

    frame_transverse_extent = float(binary.shape[1] * transverse_spacing)
    lateral_drift = float(np.max(sampled_transverse) - np.min(sampled_transverse))
    end_to_end_tilt = float((sampled_transverse[-1] - sampled_transverse[0]) / endpoint_separation)
    touches_frame = bool(
        np.any(binary[0])
        or np.any(binary[-1])
        or np.any(binary[:, 0])
        or np.any(binary[:, -1])
    )

    return MaskCenterlineProfile(
        normalized_axial=_readonly(sample_axis),
        axial_position_physical=_readonly(sample_axial),
        transverse_center_physical=_readonly(sampled_transverse),
        chord_residual_physical=_readonly(chord_residual),
        active_row_count=int(active_rows.size),
        axial_extent_physical=float(active_rows.size * axial_spacing),
        endpoint_axial_separation_physical=endpoint_separation,
        frame_transverse_extent_physical=frame_transverse_extent,
        global_transverse_centroid_physical=float(global_transverse),
        bottom_center_physical=float(sampled_transverse[0]),
        midplane_center_physical=float(sampled_transverse[count // 2]),
        top_center_physical=float(sampled_transverse[-1]),
        end_to_end_tilt_slope=end_to_end_tilt,
        lateral_drift_physical=lateral_drift,
        lateral_drift_frame_fraction=float(lateral_drift / frame_transverse_extent),
        rms_chord_deviation_physical=rms_chord,
        max_chord_deviation_physical=max_chord,
        rms_chord_deviation_frame_fraction=float(rms_chord / frame_transverse_extent),
        centerline_arclength_physical=centerline_arclength,
        endpoint_chord_length_physical=endpoint_chord,
        straightness_ratio=straightness_ratio,
        touches_image_frame=touches_frame,
        mask_provenance=mask_source,
        frame_provenance=frame_source,
    )
