"""Truth-bounded comparison of ordered public morphology trends.

The input is *already extracted* scalar morphology features over ordered frames.
This module does not segment images, infer hidden frame times, fit a camera,
register frames, or modify a velocity field.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np
from scipy.stats import spearmanr


@dataclass(frozen=True)
class MorphologyTrendFeature:
    name: str
    spearman_rho: float
    candidate_endpoint_direction: int
    reference_endpoint_direction: int
    endpoint_direction_agrees: bool
    normalized_increment_rms: float


@dataclass(frozen=True)
class MorphologyTrendReport:
    candidate_id: str
    reference_source: str
    pairing_provenance: str
    pairing_mode: str
    candidate_indices: tuple[int, ...]
    reference_indices: tuple[int, ...]
    candidate_times: tuple[float, ...] | None
    features: tuple[MorphologyTrendFeature, ...]
    mean_spearman_rho: float
    mean_normalized_increment_rms: float
    all_endpoint_directions_agree: bool
    time_registration_performed: bool = False
    hidden_frame_times_inferred: bool = False
    camera_fitted: bool = False
    segmentation_performed: bool = False
    velocity_changed: bool = False
    visualization_ready: bool = False
    visual_correspondence_verified: bool = False
    pde_validated: bool = False
    paper_exact: bool = False
    openai_field_identified: bool = False
    blowup_proved: bool = False


def _as_feature_matrix(values: object, label: str) -> np.ndarray:
    array = np.asarray(values, dtype=float)
    if array.ndim != 2:
        raise ValueError(f"{label} must have shape (frames, features)")
    if array.shape[0] < 3:
        raise ValueError(f"{label} must contain at least three frames")
    if array.shape[1] < 1:
        raise ValueError(f"{label} must contain at least one feature")
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{label} must contain only finite values")
    return array


def _validate_names(feature_names: Sequence[str], count: int) -> tuple[str, ...]:
    names = tuple(str(name).strip() for name in feature_names)
    if len(names) != count:
        raise ValueError("feature_names length must match the feature dimension")
    if any(not name for name in names):
        raise ValueError("feature_names must be non-empty")
    if len(set(names)) != len(names):
        raise ValueError("feature_names must be unique")
    return names


def _validate_times(candidate_times: object, count: int) -> tuple[float, ...] | None:
    if candidate_times is None:
        return None
    times = np.asarray(candidate_times, dtype=float)
    if times.shape != (count,):
        raise ValueError("candidate_times must have one value per candidate frame")
    if not np.all(np.isfinite(times)):
        raise ValueError("candidate_times must be finite")
    if not np.all(np.diff(times) > 0.0):
        raise ValueError("candidate_times must be strictly increasing")
    return tuple(float(value) for value in times)


def _resolve_pairs(
    candidate_count: int,
    reference_count: int,
    frame_pairs: object,
) -> tuple[str, np.ndarray, np.ndarray]:
    if frame_pairs is None:
        if candidate_count != reference_count:
            raise ValueError(
                "different frame counts require explicit frame_pairs; automatic time alignment is forbidden"
            )
        indices = np.arange(candidate_count, dtype=int)
        return "ordinal_equal_length", indices, indices.copy()

    raw = np.asarray(frame_pairs)
    if raw.ndim != 2 or raw.shape[1] != 2 or raw.shape[0] < 3:
        raise ValueError("frame_pairs must have shape (pairs, 2) with at least three pairs")
    if not np.issubdtype(raw.dtype, np.integer):
        if not np.all(np.isfinite(raw)) or not np.all(raw == np.floor(raw)):
            raise ValueError("frame_pairs must contain integer indices")
    pairs = raw.astype(int, copy=False)
    candidate_indices = pairs[:, 0]
    reference_indices = pairs[:, 1]
    if np.any(candidate_indices < 0) or np.any(candidate_indices >= candidate_count):
        raise ValueError("candidate frame index out of range")
    if np.any(reference_indices < 0) or np.any(reference_indices >= reference_count):
        raise ValueError("reference frame index out of range")
    if not np.all(np.diff(candidate_indices) > 0):
        raise ValueError("candidate frame indices must be strictly increasing")
    if not np.all(np.diff(reference_indices) > 0):
        raise ValueError("reference frame indices must be strictly increasing")
    return "explicit_monotone_pairs", candidate_indices, reference_indices


def _nonempty_text(value: object, label: str) -> str:
    text = str(value).strip()
    if not text:
        raise ValueError(f"{label} must be non-empty")
    return text


def compare_ordered_morphology_trends(
    candidate_features: object,
    reference_features: object,
    feature_names: Sequence[str],
    *,
    candidate_id: str,
    reference_source: str,
    pairing_provenance: str,
    candidate_times: object = None,
    frame_pairs: object = None,
) -> MorphologyTrendReport:
    """Compare ordered morphology evolution without fitting hidden frame times.

    Spearman rho tests whether each public-observable feature evolves in the same
    rank order. A second scale/offset-invariant quantity compares the signed
    increment profile on the *declared* frame pairing. Neither metric is a pass
    threshold, and neither establishes PDE validity or hidden-field identity.
    """

    candidate = _as_feature_matrix(candidate_features, "candidate_features")
    reference = _as_feature_matrix(reference_features, "reference_features")
    if candidate.shape[1] != reference.shape[1]:
        raise ValueError("candidate and reference feature dimensions must match")
    names = _validate_names(feature_names, candidate.shape[1])
    times = _validate_times(candidate_times, candidate.shape[0])
    pairing_mode, candidate_indices, reference_indices = _resolve_pairs(
        candidate.shape[0], reference.shape[0], frame_pairs
    )

    candidate_id_text = _nonempty_text(candidate_id, "candidate_id")
    reference_source_text = _nonempty_text(reference_source, "reference_source")
    provenance_text = _nonempty_text(pairing_provenance, "pairing_provenance")

    selected_candidate = candidate[candidate_indices]
    selected_reference = reference[reference_indices]
    if selected_candidate.shape[0] < 3:
        raise ValueError("at least three paired frames are required")

    feature_reports: list[MorphologyTrendFeature] = []
    for column, name in enumerate(names):
        candidate_series = selected_candidate[:, column]
        reference_series = selected_reference[:, column]
        candidate_diffs = np.diff(candidate_series)
        reference_diffs = np.diff(reference_series)
        candidate_total_variation = float(np.sum(np.abs(candidate_diffs)))
        reference_total_variation = float(np.sum(np.abs(reference_diffs)))
        if candidate_total_variation == 0.0:
            raise ValueError(f"candidate feature {name!r} is constant on paired frames")
        if reference_total_variation == 0.0:
            raise ValueError(f"reference feature {name!r} is constant on paired frames")

        rho = float(spearmanr(candidate_series, reference_series).statistic)
        if not np.isfinite(rho):
            raise ValueError(f"Spearman rho is undefined for feature {name!r}")

        candidate_increment_profile = candidate_diffs / candidate_total_variation
        reference_increment_profile = reference_diffs / reference_total_variation
        increment_rms = float(
            np.sqrt(np.mean((candidate_increment_profile - reference_increment_profile) ** 2))
        )
        candidate_direction = int(np.sign(candidate_series[-1] - candidate_series[0]))
        reference_direction = int(np.sign(reference_series[-1] - reference_series[0]))
        feature_reports.append(
            MorphologyTrendFeature(
                name=name,
                spearman_rho=rho,
                candidate_endpoint_direction=candidate_direction,
                reference_endpoint_direction=reference_direction,
                endpoint_direction_agrees=candidate_direction == reference_direction,
                normalized_increment_rms=increment_rms,
            )
        )

    mean_rho = float(np.mean([item.spearman_rho for item in feature_reports]))
    mean_increment_rms = float(
        np.mean([item.normalized_increment_rms for item in feature_reports])
    )
    return MorphologyTrendReport(
        candidate_id=candidate_id_text,
        reference_source=reference_source_text,
        pairing_provenance=provenance_text,
        pairing_mode=pairing_mode,
        candidate_indices=tuple(int(value) for value in candidate_indices),
        reference_indices=tuple(int(value) for value in reference_indices),
        candidate_times=times,
        features=tuple(feature_reports),
        mean_spearman_rho=mean_rho,
        mean_normalized_increment_rms=mean_increment_rms,
        all_endpoint_directions_agree=all(
            item.endpoint_direction_agrees for item in feature_reports
        ),
    )
