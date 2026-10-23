from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Iterable, Sequence

import numpy as np


@dataclass(frozen=True)
class TemporalSnapshotLevel:
    snapshot_count: int
    max_time_spacing: float
    vector_error_max: float
    vector_error_rms: float
    relative_rms: float


@dataclass(frozen=True)
class TemporalSnapshotConsistencyReport:
    time_interval: tuple[float, float]
    point_count: int
    probe_time_count: int
    levels: tuple[TemporalSnapshotLevel, ...]
    direct_velocity_rms: float
    claim_scope: str = "visualization_time_sampling_only"
    pde_validated: bool = False
    visual_correspondence_verified: bool = False
    paper_exact: bool = False
    openai_field_identified: bool = False

    def to_dict(self) -> dict[str, Any]:
        return {
            "time_interval": list(self.time_interval),
            "point_count": self.point_count,
            "probe_time_count": self.probe_time_count,
            "direct_velocity_rms": self.direct_velocity_rms,
            "levels": [
                {
                    "snapshot_count": level.snapshot_count,
                    "max_time_spacing": level.max_time_spacing,
                    "vector_error_max": level.vector_error_max,
                    "vector_error_rms": level.vector_error_rms,
                    "relative_rms": level.relative_rms,
                }
                for level in self.levels
            ],
            "claim_scope": self.claim_scope,
            "pde_validated": self.pde_validated,
            "visual_correspondence_verified": self.visual_correspondence_verified,
            "paper_exact": self.paper_exact,
            "openai_field_identified": self.openai_field_identified,
        }


def _evaluate(field: Any, points: np.ndarray, time: float) -> np.ndarray:
    values = np.asarray(field.at_points(points, float(time)), dtype=float)
    if values.shape != (points.shape[0], 3):
        raise ValueError(
            "field.at_points(points, time) must return shape (point_count, 3)"
        )
    if not np.all(np.isfinite(values)):
        raise ValueError("field returned non-finite velocity")
    return values


def _validate_snapshot_counts(snapshot_counts: Iterable[int]) -> tuple[int, ...]:
    counts = tuple(int(value) for value in snapshot_counts)
    if len(counts) < 3:
        raise ValueError("at least three temporal snapshot levels are required")
    if any(value < 3 for value in counts):
        raise ValueError("each snapshot count must be at least 3")
    if any(float(raw) != int(raw) for raw in snapshot_counts):
        raise ValueError("snapshot counts must be integers")
    if any(b <= a for a, b in zip(counts, counts[1:])):
        raise ValueError("snapshot counts must be strictly increasing")
    return counts


def audit_temporal_snapshot_consistency(
    field: Any,
    points: Sequence[Sequence[float]] | np.ndarray,
    probe_times: Sequence[float] | np.ndarray,
    *,
    snapshot_counts: Sequence[int] = (3, 5, 9),
    time_interval: tuple[float, float] = (0.25, 0.75),
    min_reference_rms: float = 1.0e-12,
) -> TemporalSnapshotConsistencyReport:
    """Compare direct public velocity calls to linearly interpolated time snapshots.

    Spatial points are held fixed. Only temporal sampling density changes, so this
    audit isolates whether a visualization assembled from a finite number of saved
    time slices can hide or invent time evolution relative to direct field calls.
    It is not a PDE validator and does not define a visual-correspondence threshold.
    """

    pts = np.asarray(points, dtype=float)
    if pts.ndim != 2 or pts.shape[1] != 3 or pts.shape[0] == 0:
        raise ValueError("points must have shape (N, 3) with N > 0")
    if not np.all(np.isfinite(pts)):
        raise ValueError("points must be finite")

    times = np.asarray(probe_times, dtype=float)
    if times.ndim != 1 or times.size < 3:
        raise ValueError("probe_times must be a one-dimensional array with >= 3 values")
    if not np.all(np.isfinite(times)):
        raise ValueError("probe_times must be finite")
    if np.any(np.diff(times) <= 0.0):
        raise ValueError("probe_times must be strictly increasing")

    t0, t1 = (float(time_interval[0]), float(time_interval[1]))
    if not np.isfinite(t0) or not np.isfinite(t1) or not t0 < t1:
        raise ValueError("time_interval must be finite and strictly increasing")
    if times[0] <= t0 or times[-1] >= t1:
        raise ValueError("probe_times must lie strictly inside time_interval")
    if not np.isfinite(min_reference_rms) or min_reference_rms <= 0.0:
        raise ValueError("min_reference_rms must be positive and finite")

    counts = _validate_snapshot_counts(snapshot_counts)

    direct = np.stack([_evaluate(field, pts, t) for t in times], axis=0)
    direct_rms = float(np.sqrt(np.mean(direct * direct)))
    if direct_rms < min_reference_rms:
        raise ValueError("reference field is numerically inactive on the audit probes")

    levels: list[TemporalSnapshotLevel] = []
    for count in counts:
        snapshots = np.linspace(t0, t1, count, dtype=float)
        sampled = np.stack([_evaluate(field, pts, t) for t in snapshots], axis=0)

        right = np.searchsorted(snapshots, times, side="right")
        right = np.clip(right, 1, count - 1)
        left = right - 1
        t_left = snapshots[left]
        t_right = snapshots[right]
        alpha = ((times - t_left) / (t_right - t_left))[:, None, None]
        interpolated = (1.0 - alpha) * sampled[left] + alpha * sampled[right]

        error = interpolated - direct
        vector_error = np.linalg.norm(error, axis=-1)
        vector_error_max = float(np.max(vector_error))
        vector_error_rms = float(np.sqrt(np.mean(error * error)))
        levels.append(
            TemporalSnapshotLevel(
                snapshot_count=count,
                max_time_spacing=float(np.max(np.diff(snapshots))),
                vector_error_max=vector_error_max,
                vector_error_rms=vector_error_rms,
                relative_rms=vector_error_rms / direct_rms,
            )
        )

    return TemporalSnapshotConsistencyReport(
        time_interval=(t0, t1),
        point_count=int(pts.shape[0]),
        probe_time_count=int(times.size),
        levels=tuple(levels),
        direct_velocity_rms=direct_rms,
    )
