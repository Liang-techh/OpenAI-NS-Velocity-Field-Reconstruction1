"""Truth-bounded connectivity features for explicitly presegmented visual masks.

This module measures connectedness and axial/transverse fragmentation in an
already-presegmented 2-D mask.  It does not segment imagery, bridge gaps, clean
components, fit a camera, register frames, infer hidden velocity data, or
define a visual acceptance threshold.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.ndimage import label


_FOUR_CONNECTED = np.array(
    [[False, True, False], [True, True, True], [False, True, False]], dtype=bool
)


@dataclass(frozen=True)
class MaskConnectivityFeatures:
    """Deterministic fragmentation observables for one presegmented mask."""

    component_sizes_pixels: NDArray[np.int64]
    component_count: int
    foreground_pixels: int
    foreground_area_physical: float
    largest_component_pixels: int
    largest_component_fraction: float
    detached_foreground_fraction: float
    active_row_count: int
    axial_segment_count: int
    axial_gap_rows: int
    largest_axial_gap_rows: int
    largest_axial_gap_physical: float
    fragmented_active_row_fraction: float
    maximum_transverse_runs: int
    touches_image_frame: bool
    connectivity: str
    mask_provenance: str
    frame_provenance: str
    segmentation_performed: bool = False
    threshold_selected: bool = False
    gap_bridging_performed: bool = False
    component_cleanup_performed: bool = False
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


def _readonly_int(array: NDArray[np.int64]) -> NDArray[np.int64]:
    out = np.asarray(array, dtype=np.int64)
    out.setflags(write=False)
    return out


def _transverse_run_count(row: NDArray[np.bool_]) -> int:
    starts = int(row[0]) + int(np.count_nonzero((~row[:-1]) & row[1:]))
    return starts


def measure_presegmented_mask_connectivity(
    mask: ArrayLike,
    *,
    pixel_spacing: tuple[float, float] = (1.0, 1.0),
    mask_provenance: str,
    frame_provenance: str,
) -> MaskConnectivityFeatures:
    """Measure fixed 4-neighbor connectivity and fragmentation observables.

    Connected components use one fixed 4-neighbor structure.  Diagonal-only
    contacts therefore remain distinct.  This choice is intentionally not
    exposed as a tuning parameter so candidate/public comparisons cannot change
    topology post hoc.

    ``axial_segment_count`` counts contiguous groups of active image rows and
    ``axial_gap_rows`` counts blank rows between the first and last active row.
    ``fragmented_active_row_fraction`` records the fraction of active rows with
    more than one transverse foreground run, distinguishing lateral splitting
    or visible holes from a single uninterrupted row span.
    """

    binary = _binary_mask(mask)
    axial_spacing, transverse_spacing = _positive_spacing(pixel_spacing)
    mask_source = _provenance(mask_provenance, "mask_provenance")
    frame_source = _provenance(frame_provenance, "frame_provenance")

    labels, component_count = label(binary, structure=_FOUR_CONNECTED)
    component_count = int(component_count)
    if component_count < 1:
        raise ValueError("connected-component labeling returned no foreground component")

    sizes = np.bincount(labels.ravel(), minlength=component_count + 1)[1:].astype(np.int64)
    if sizes.size != component_count or bool(np.any(sizes <= 0)):
        raise ValueError("connected-component labeling produced an invalid component census")
    sizes = np.sort(sizes)[::-1]

    foreground_pixels = int(np.count_nonzero(binary))
    if int(np.sum(sizes)) != foreground_pixels:
        raise ValueError("component census does not conserve foreground pixels")

    active_rows = np.flatnonzero(np.any(binary, axis=1))
    if active_rows.size == 0:  # guarded by _binary_mask; kept fail-closed.
        raise ValueError("active-row extraction failed")

    row_diffs = np.diff(active_rows)
    gap_sizes = row_diffs[row_diffs > 1] - 1
    axial_segment_count = 1 + int(gap_sizes.size)
    axial_gap_rows = int(np.sum(gap_sizes)) if gap_sizes.size else 0
    largest_axial_gap_rows = int(np.max(gap_sizes)) if gap_sizes.size else 0

    run_counts = np.array(
        [_transverse_run_count(binary[row]) for row in active_rows], dtype=np.int64
    )
    if bool(np.any(run_counts < 1)):
        raise ValueError("active-row run counting produced an invalid result")
    fragmented_rows = int(np.count_nonzero(run_counts > 1))

    largest_component_pixels = int(sizes[0])
    largest_fraction = largest_component_pixels / foreground_pixels
    detached_fraction = 1.0 - largest_fraction

    touches_frame = bool(
        np.any(binary[0])
        or np.any(binary[-1])
        or np.any(binary[:, 0])
        or np.any(binary[:, -1])
    )

    return MaskConnectivityFeatures(
        component_sizes_pixels=_readonly_int(sizes),
        component_count=component_count,
        foreground_pixels=foreground_pixels,
        foreground_area_physical=float(foreground_pixels * axial_spacing * transverse_spacing),
        largest_component_pixels=largest_component_pixels,
        largest_component_fraction=float(largest_fraction),
        detached_foreground_fraction=float(detached_fraction),
        active_row_count=int(active_rows.size),
        axial_segment_count=axial_segment_count,
        axial_gap_rows=axial_gap_rows,
        largest_axial_gap_rows=largest_axial_gap_rows,
        largest_axial_gap_physical=float(largest_axial_gap_rows * axial_spacing),
        fragmented_active_row_fraction=float(fragmented_rows / active_rows.size),
        maximum_transverse_runs=int(np.max(run_counts)),
        touches_image_frame=touches_frame,
        connectivity="fixed_4_neighbor",
        mask_provenance=mask_source,
        frame_provenance=frame_source,
    )
