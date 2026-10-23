"""Truth-bounded packing for efficient colored 3D streamline rendering.

This module does not integrate streamlines and never changes a velocity field.
It packs already-computed polylines into collection-friendly arrays so dense
visualizations do not require one graphics object per tiny line segment.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np


_CLAIM_SCOPE = "visualization_render_packing_only"


def _readonly(array: np.ndarray) -> np.ndarray:
    out = np.asarray(array)
    out.setflags(write=False)
    return out


@dataclass(frozen=True)
class StreamlineCollectionPack:
    """Packed geometry and scalar data for collection-based rendering."""

    segments: np.ndarray
    segment_scalars: np.ndarray
    line_segment_offsets: np.ndarray
    matlab_x: np.ndarray
    matlab_y: np.ndarray
    matlab_z: np.ndarray
    matlab_c: np.ndarray
    line_count: int
    segment_count: int
    scalar_min: float
    scalar_max: float
    claim_scope: str = _CLAIM_SCOPE
    visualization_ready: bool = False
    visual_correspondence_verified: bool = False
    pde_validated: bool = False
    paper_exact: bool = False
    openai_field_identified: bool = False
    blowup_proved: bool = False


def _coerce_lines(lines: Sequence[np.ndarray]) -> tuple[np.ndarray, ...]:
    if isinstance(lines, np.ndarray):
        if lines.ndim != 3 or lines.shape[2] != 3:
            raise ValueError("array input must have shape (n_lines,n_points,3)")
        seq = tuple(np.asarray(line, dtype=float) for line in lines)
    else:
        try:
            seq = tuple(np.asarray(line, dtype=float) for line in lines)
        except TypeError as exc:
            raise ValueError("lines must be a non-empty sequence of polylines") from exc

    if not seq:
        raise ValueError("at least one streamline is required")

    checked: list[np.ndarray] = []
    for index, line in enumerate(seq):
        if line.ndim != 2 or line.shape[1] != 3:
            raise ValueError(f"line {index} must have shape (n_points,3)")
        if line.shape[0] < 2:
            raise ValueError(f"line {index} must contain at least two points")
        if not np.all(np.isfinite(line)):
            raise ValueError(f"line {index} contains non-finite coordinates")
        if not np.any(np.linalg.norm(np.diff(line, axis=0), axis=1) > 0.0):
            raise ValueError(f"line {index} is geometrically degenerate")
        checked.append(line)
    return tuple(checked)


def _coerce_point_scalars(
    lines: tuple[np.ndarray, ...], point_scalars: Sequence[np.ndarray] | np.ndarray | None
) -> tuple[np.ndarray, ...]:
    if point_scalars is None:
        return tuple(line[:, 2].astype(float, copy=True) for line in lines)

    if isinstance(point_scalars, np.ndarray) and point_scalars.ndim == 2:
        if point_scalars.shape[0] != len(lines):
            raise ValueError("point_scalars line count does not match lines")
        seq = tuple(np.asarray(row, dtype=float) for row in point_scalars)
    else:
        try:
            seq = tuple(np.asarray(values, dtype=float) for values in point_scalars)
        except TypeError as exc:
            raise ValueError("point_scalars must match the streamline sequence") from exc

    if len(seq) != len(lines):
        raise ValueError("point_scalars line count does not match lines")

    checked: list[np.ndarray] = []
    for index, (line, values) in enumerate(zip(lines, seq)):
        if values.ndim != 1 or values.shape[0] != line.shape[0]:
            raise ValueError(
                f"point_scalars[{index}] must have one finite scalar per streamline point"
            )
        if not np.all(np.isfinite(values)):
            raise ValueError(f"point_scalars[{index}] contains non-finite values")
        checked.append(values)
    return tuple(checked)


def _nan_join(values: tuple[np.ndarray, ...]) -> np.ndarray:
    total = sum(value.size for value in values) + len(values) - 1
    joined = np.empty(total, dtype=float)
    cursor = 0
    for index, value in enumerate(values):
        n = value.size
        joined[cursor : cursor + n] = value
        cursor += n
        if index + 1 < len(values):
            joined[cursor] = np.nan
            cursor += 1
    return joined


def pack_streamline_collection(
    lines: Sequence[np.ndarray] | np.ndarray,
    point_scalars: Sequence[np.ndarray] | np.ndarray | None = None,
) -> StreamlineCollectionPack:
    """Pack streamline polylines for one collection/surface graphics object.

    Parameters
    ----------
    lines:
        Sequence of ``(n_i,3)`` Cartesian polylines, or one rectangular
        ``(n_lines,n_points,3)`` array.
    point_scalars:
        Optional scalar value at every input point, for example frozen-time
        speed. If omitted, z/height is used. Scalars are only visualization
        data and never affect geometry.

    Returns
    -------
    StreamlineCollectionPack
        ``segments``/``segment_scalars`` can be passed to a single
        Matplotlib ``Line3DCollection``. ``matlab_x/y/z/c`` are ``(2,N)``
        NaN-separated arrays for a single MATLAB ``surface`` object with
        ``FaceColor='none'`` and ``EdgeColor='interp'``.
    """

    polylines = _coerce_lines(lines)
    scalars = _coerce_point_scalars(polylines, point_scalars)

    segment_blocks: list[np.ndarray] = []
    scalar_blocks: list[np.ndarray] = []
    offsets = [0]
    for line, values in zip(polylines, scalars):
        block = np.stack((line[:-1], line[1:]), axis=1)
        segment_blocks.append(block)
        scalar_blocks.append(0.5 * (values[:-1] + values[1:]))
        offsets.append(offsets[-1] + block.shape[0])

    segments = np.concatenate(segment_blocks, axis=0)
    segment_scalars = np.concatenate(scalar_blocks)
    if segments.shape != (segment_scalars.size, 2, 3):
        raise RuntimeError("internal segment packing mismatch")
    if not np.all(np.isfinite(segments)) or not np.all(np.isfinite(segment_scalars)):
        raise RuntimeError("internal packing produced non-finite collection data")

    joined_x = _nan_join(tuple(line[:, 0] for line in polylines))
    joined_y = _nan_join(tuple(line[:, 1] for line in polylines))
    joined_z = _nan_join(tuple(line[:, 2] for line in polylines))
    joined_c = _nan_join(scalars)

    matlab_x = np.vstack((joined_x, joined_x))
    matlab_y = np.vstack((joined_y, joined_y))
    matlab_z = np.vstack((joined_z, joined_z))
    matlab_c = np.vstack((joined_c, joined_c))

    return StreamlineCollectionPack(
        segments=_readonly(segments),
        segment_scalars=_readonly(segment_scalars),
        line_segment_offsets=_readonly(np.asarray(offsets, dtype=np.int64)),
        matlab_x=_readonly(matlab_x),
        matlab_y=_readonly(matlab_y),
        matlab_z=_readonly(matlab_z),
        matlab_c=_readonly(matlab_c),
        line_count=len(polylines),
        segment_count=int(segment_scalars.size),
        scalar_min=float(np.min(segment_scalars)),
        scalar_max=float(np.max(segment_scalars)),
    )


def make_matplotlib_line3d_collection(
    pack: StreamlineCollectionPack,
    *,
    linewidth: float = 0.35,
    cmap: str = "turbo",
):
    """Create one Matplotlib ``Line3DCollection`` from a packed bundle.

    Matplotlib is a development/visualization dependency in this repository;
    importing it is delayed so numerical users of the packing primitive do not
    require Matplotlib at module import time.
    """

    if not isinstance(pack, StreamlineCollectionPack):
        raise TypeError("pack must be a StreamlineCollectionPack")
    if not np.isfinite(linewidth) or linewidth <= 0.0:
        raise ValueError("linewidth must be finite and positive")
    if not isinstance(cmap, str) or not cmap.strip():
        raise ValueError("cmap must be a non-empty Matplotlib colormap name")

    from mpl_toolkits.mplot3d.art3d import Line3DCollection

    collection = Line3DCollection(pack.segments, linewidths=float(linewidth), cmap=cmap)
    collection.set_array(pack.segment_scalars)
    collection.set_clim(pack.scalar_min, pack.scalar_max)
    return collection
