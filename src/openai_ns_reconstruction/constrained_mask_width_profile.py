"""Truth-bounded axial width profiles for explicitly presegmented visual masks.

This module measures an already-presegmented 2-D mask.  It does not segment
imagery, fit a camera, register frames, infer hidden velocity data, or define a
visual acceptance threshold.  The normalized axial coordinate is a descriptor
for taper/waist shape only; the physical axial extent is returned separately so
normalization cannot be mistaken for a fitted geometric rescaling.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.interpolate import PchipInterpolator


@dataclass(frozen=True)
class MaskWidthProfile:
    """Deterministic width-profile observables for one presegmented mask."""

    normalized_axial: NDArray[np.float64]
    span_width_physical: NDArray[np.float64]
    occupied_width_physical: NDArray[np.float64]
    span_width_fraction: NDArray[np.float64]
    occupancy_fraction: NDArray[np.float64]
    active_row_count: int
    axial_extent_physical: float
    maximum_span_width_physical: float
    minimum_span_width_physical: float
    bottom_tip_span_fraction: float
    midplane_span_fraction: float
    top_tip_span_fraction: float
    central_waist_span_fraction: float
    maximum_width_axial_location: float
    touches_image_frame: bool
    mask_provenance: str
    frame_provenance: str
    segmentation_performed: bool = False
    threshold_selected: bool = False
    camera_fitted: bool = False
    registration_fitted: bool = False
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


def _provenance(value: str, name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{name} must be a non-empty provenance string")
    return value.strip()


def _sample_count(value: int) -> int:
    if isinstance(value, bool) or not isinstance(value, (int, np.integer)):
        raise ValueError("sample_count must be an odd integer")
    value = int(value)
    if value < 5 or value > 4097 or value % 2 == 0:
        raise ValueError("sample_count must be an odd integer in [5, 4097]")
    return value


def _readonly(array: NDArray[np.float64]) -> NDArray[np.float64]:
    out = np.asarray(array, dtype=np.float64)
    out.setflags(write=False)
    return out


def measure_presegmented_mask_width_profile(
    mask: ArrayLike,
    *,
    pixel_spacing: tuple[float, float] = (1.0, 1.0),
    sample_count: int = 129,
    mask_provenance: str,
    frame_provenance: str,
) -> MaskWidthProfile:
    """Measure fixed-axis taper/waist observables from a presegmented mask.

    Rows containing foreground must form one contiguous axial interval.  For
    each active row, ``span_width`` is the distance from the leftmost through
    rightmost foreground pixel (inclusive), while ``occupied_width`` counts
    only foreground pixels.  The latter therefore retains evidence of visible
    holes/concavities inside the outer span.

    Rowwise profiles are sampled on a fixed axial-up coordinate from -1
    (bottom of the visible object) to +1 (top) with SciPy PCHIP.  This
    normalization is descriptor-only: ``axial_extent_physical`` and physical
    widths remain explicit outputs and no image registration is performed.
    """

    binary = _binary_mask(mask)
    axial_spacing, transverse_spacing = _positive_spacing(pixel_spacing)
    count = _sample_count(sample_count)
    mask_source = _provenance(mask_provenance, "mask_provenance")
    frame_source = _provenance(frame_provenance, "frame_provenance")

    active_rows = np.flatnonzero(np.any(binary, axis=1))
    if active_rows.size < 3:
        raise ValueError("at least three active axial rows are required")
    if np.any(np.diff(active_rows) != 1):
        raise ValueError("foreground rows must form one contiguous axial interval")

    span_pixels = np.empty(active_rows.size, dtype=np.float64)
    occupied_pixels = np.empty(active_rows.size, dtype=np.float64)
    for index, row in enumerate(active_rows):
        columns = np.flatnonzero(binary[row])
        if columns.size == 0:  # guarded by active_rows; kept fail-closed.
            raise ValueError("active-row extraction failed")
        span_pixels[index] = float(columns[-1] - columns[0] + 1)
        occupied_pixels[index] = float(columns.size)

    if float(np.max(span_pixels)) < 2.0:
        raise ValueError("mask is transversely degenerate")

    # Convert image row order (top -> bottom) into axial-up order (bottom -> top).
    span_width = span_pixels[::-1] * transverse_spacing
    occupied_width = occupied_pixels[::-1] * transverse_spacing
    raw_axis = np.linspace(-1.0, 1.0, active_rows.size, dtype=np.float64)
    sample_axis = np.linspace(-1.0, 1.0, count, dtype=np.float64)

    span_interp = PchipInterpolator(raw_axis, span_width, extrapolate=False)
    occupied_interp = PchipInterpolator(raw_axis, occupied_width, extrapolate=False)
    sampled_span = np.asarray(span_interp(sample_axis), dtype=np.float64)
    sampled_occupied = np.asarray(occupied_interp(sample_axis), dtype=np.float64)
    if not np.all(np.isfinite(sampled_span)) or not np.all(np.isfinite(sampled_occupied)):
        raise ValueError("profile interpolation produced nonfinite values")
    if np.any(sampled_span <= 0.0) or np.any(sampled_occupied < 0.0):
        raise ValueError("profile interpolation produced an invalid width")

    maximum_span = float(np.max(span_width))
    span_fraction = sampled_span / maximum_span
    occupancy_fraction = sampled_occupied / sampled_span
    tolerance = 64.0 * np.finfo(float).eps
    if np.any(occupancy_fraction < -tolerance) or np.any(occupancy_fraction > 1.0 + tolerance):
        raise ValueError("interpolated occupied width escaped the row span")
    occupancy_fraction = np.clip(occupancy_fraction, 0.0, 1.0)

    central = np.abs(sample_axis) <= 0.6
    if not bool(np.any(central)):
        raise ValueError("central profile window is empty")
    maximum_mask = np.isclose(sampled_span, np.max(sampled_span), rtol=0.0, atol=1e-12 * maximum_span)
    maximum_location = float(np.mean(sample_axis[maximum_mask]))

    touches_frame = bool(
        np.any(binary[0])
        or np.any(binary[-1])
        or np.any(binary[:, 0])
        or np.any(binary[:, -1])
    )

    return MaskWidthProfile(
        normalized_axial=_readonly(sample_axis),
        span_width_physical=_readonly(sampled_span),
        occupied_width_physical=_readonly(sampled_occupied),
        span_width_fraction=_readonly(span_fraction),
        occupancy_fraction=_readonly(occupancy_fraction),
        active_row_count=int(active_rows.size),
        axial_extent_physical=float(active_rows.size * axial_spacing),
        maximum_span_width_physical=maximum_span,
        minimum_span_width_physical=float(np.min(span_width)),
        bottom_tip_span_fraction=float(span_fraction[0]),
        midplane_span_fraction=float(span_fraction[count // 2]),
        top_tip_span_fraction=float(span_fraction[-1]),
        central_waist_span_fraction=float(np.min(span_fraction[central])),
        maximum_width_axial_location=maximum_location,
        touches_image_frame=touches_frame,
        mask_provenance=mask_source,
        frame_provenance=frame_source,
    )
