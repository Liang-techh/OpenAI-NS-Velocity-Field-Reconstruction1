from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping

import numpy as np
from scipy.spatial.distance import directed_hausdorff


@dataclass(frozen=True)
class StreamlineResolutionReport:
    resolutions: np.ndarray
    s_fraction: np.ndarray
    resampled_points: np.ndarray
    arclengths: np.ndarray
    arclength_relative_error_to_finest: np.ndarray
    symmetric_hausdorff: np.ndarray
    normalized_hausdorff_by_finest_arclength: np.ndarray
    pointwise_rms: np.ndarray
    normalized_pointwise_rms_by_finest_arclength: np.ndarray
    start_error: np.ndarray
    end_error: np.ndarray
    finest_resolution: int
    finest_arclength: float
    max_normalized_hausdorff: float
    max_normalized_pointwise_rms: float
    metadata: Mapping[str, object]


def _validate_path(points: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    p = np.asarray(points, dtype=float)
    if p.ndim != 2 or p.shape[1] != 3 or p.shape[0] < 4:
        raise ValueError("each streamline must have shape (n,3) with n>=4")
    if not np.all(np.isfinite(p)):
        raise ValueError("streamline points must be finite")
    seg = np.linalg.norm(np.diff(p, axis=0), axis=1)
    if np.any(~np.isfinite(seg)) or np.any(seg <= 0.0):
        raise ValueError("consecutive duplicate/zero-length streamline samples are not allowed")
    s = np.concatenate(([0.0], np.cumsum(seg)))
    length = float(s[-1])
    if not np.isfinite(length) or length <= 0.0:
        raise ValueError("streamline arclength must be positive")
    return p, s, length


def audit_streamline_grid_resolution(
    resolution_to_points: Mapping[int, np.ndarray],
    *,
    sample_count: int = 129,
    seed_xyz: tuple[float, float, float],
    time: float,
    field_identity: str,
    integration_contract: str,
    resolution_provenance: str,
    provenance: str,
) -> StreamlineResolutionReport:
    """Compare one same-seed 3-D streamline across frozen grid resolutions.

    Inputs are already-integrated streamlines. This routine performs no velocity
    interpolation, streamline integration, alignment, reversal, camera fit,
    candidate fitting, or visual pass/fail decision.
    """
    if not isinstance(resolution_to_points, Mapping) or len(resolution_to_points) < 3:
        raise ValueError("at least three resolution levels are required")
    if not isinstance(sample_count, int) or sample_count < 33 or sample_count % 2 == 0:
        raise ValueError("sample_count must be an odd integer >=33")
    if (
        not isinstance(field_identity, str)
        or not field_identity.strip()
        or not isinstance(integration_contract, str)
        or not integration_contract.strip()
        or not isinstance(resolution_provenance, str)
        or not resolution_provenance.strip()
        or not isinstance(provenance, str)
        or not provenance.strip()
    ):
        raise ValueError("field, integration, resolution, and comparison provenance are required")

    if isinstance(time, bool) or not np.isscalar(time) or not np.isfinite(float(time)):
        raise ValueError("time must be one finite scalar")
    time_value = float(time)

    seed = np.asarray(seed_xyz, dtype=float)
    if seed.shape != (3,) or not np.all(np.isfinite(seed)):
        raise ValueError("seed_xyz must contain exactly three finite coordinates")

    levels: list[int] = []
    for value in resolution_to_points:
        if isinstance(value, bool) or not isinstance(value, (int, np.integer)):
            raise ValueError("resolution keys must be positive integers")
        level = int(value)
        if level <= 1:
            raise ValueError("resolution keys must be integers greater than one")
        levels.append(level)
    levels.sort()
    if np.any(np.diff(np.asarray(levels, dtype=int)) <= 0):
        raise ValueError("resolution levels must be strictly increasing")

    s_fraction = np.linspace(0.0, 1.0, sample_count)
    resampled = np.empty((len(levels), sample_count, 3), dtype=float)
    arclengths = np.empty(len(levels), dtype=float)

    for i, level in enumerate(levels):
        p, s, length = _validate_path(resolution_to_points[level])
        arclengths[i] = length
        target = s_fraction * length
        for j in range(3):
            resampled[i, :, j] = np.interp(target, s, p[:, j])

    finest_points = resampled[-1]
    finest_length = float(arclengths[-1])
    if finest_length <= 0.0:
        raise ValueError("finest streamline arclength must be positive")

    arclength_relative_error = np.abs(arclengths - finest_length) / finest_length
    hausdorff = np.zeros(len(levels), dtype=float)
    pointwise_rms = np.zeros(len(levels), dtype=float)
    start_error = np.zeros(len(levels), dtype=float)
    end_error = np.zeros(len(levels), dtype=float)

    for i in range(len(levels) - 1):
        points = resampled[i]
        forward = float(directed_hausdorff(points, finest_points, 0)[0])
        backward = float(directed_hausdorff(finest_points, points, 0)[0])
        hausdorff[i] = max(forward, backward)
        delta = points - finest_points
        pointwise_rms[i] = float(np.sqrt(np.mean(np.sum(delta * delta, axis=1))))
        start_error[i] = float(np.linalg.norm(points[0] - finest_points[0]))
        end_error[i] = float(np.linalg.norm(points[-1] - finest_points[-1]))

    normalized_hausdorff = hausdorff / finest_length
    normalized_pointwise_rms = pointwise_rms / finest_length
    if np.any(~np.isfinite(normalized_hausdorff)) or np.any(~np.isfinite(normalized_pointwise_rms)):
        raise ValueError("nonfinite resolution discrepancy")

    for arr in (
        s_fraction,
        resampled,
        arclengths,
        arclength_relative_error,
        hausdorff,
        normalized_hausdorff,
        pointwise_rms,
        normalized_pointwise_rms,
        start_error,
        end_error,
    ):
        arr.setflags(write=False)
    resolutions = np.asarray(levels, dtype=int)
    resolutions.setflags(write=False)

    metadata = {
        "claim_scope": "visualization_streamline_grid_resolution_only",
        "same_seed_required": True,
        "same_time_required": True,
        "same_field_required": True,
        "same_integration_contract_required": True,
        "resampling": "piecewise-linear normalized-arclength resampling via numpy.interp",
        "hausdorff_metric": (
            "symmetric Euclidean Hausdorff = max of two "
            "scipy.spatial.distance.directed_hausdorff calls"
        ),
        "automatic_alignment": False,
        "automatic_reversal": False,
        "camera_or_registration_fit": False,
        "visual_pass_threshold": None,
        "velocity_changed": False,
        "grid_derivatives_used_for_pde": False,
        "visualization_ready": False,
        "visual_correspondence_verified": False,
        "pde_validated": False,
        "paper_exact": False,
        "openai_field_identified": False,
        "blowup_proved": False,
        "seed_xyz": tuple(float(v) for v in seed),
        "time": time_value,
        "field_identity": field_identity,
        "integration_contract": integration_contract,
        "resolution_provenance": resolution_provenance,
        "provenance": provenance,
    }

    return StreamlineResolutionReport(
        resolutions=resolutions,
        s_fraction=s_fraction,
        resampled_points=resampled,
        arclengths=arclengths,
        arclength_relative_error_to_finest=arclength_relative_error,
        symmetric_hausdorff=hausdorff,
        normalized_hausdorff_by_finest_arclength=normalized_hausdorff,
        pointwise_rms=pointwise_rms,
        normalized_pointwise_rms_by_finest_arclength=normalized_pointwise_rms,
        start_error=start_error,
        end_error=end_error,
        finest_resolution=int(levels[-1]),
        finest_arclength=finest_length,
        max_normalized_hausdorff=float(np.max(normalized_hausdorff[:-1])),
        max_normalized_pointwise_rms=float(np.max(normalized_pointwise_rms[:-1])),
        metadata=metadata,
    )
