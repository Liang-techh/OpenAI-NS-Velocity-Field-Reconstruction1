"""Truth-bounded fidelity audit for frozen-grid velocity replay callables.

This audit compares already-frozen visualization/delivery replays against the
reference velocity callable at identical fixed interior probes. It measures
representation error only: no derivatives are taken and no result is PDE,
divergence, or OpenAI-correspondence evidence.
"""

from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
import math
import re
from typing import Callable, Mapping, Sequence

import numpy as np
from scipy.stats import qmc

VelocityCallable = Callable[[np.ndarray, np.ndarray, np.ndarray, float], np.ndarray]
_SHA256_RE = re.compile(r"^[0-9a-f]{64}$")


@dataclass(frozen=True)
class ReplayLevel:
    label: str
    resolution_xyz: tuple[int, int, int]
    velocity: VelocityCallable


@dataclass(frozen=True)
class ReplayFidelityMetrics:
    label: str
    resolution_xyz: tuple[int, int, int]
    vector_rms_error: float
    vector_max_error: float
    normalized_vector_rms_error: float
    normalized_vector_max_error: float
    component_rms_error: np.ndarray
    mean_speed_bias: float


@dataclass(frozen=True)
class VelocityReplayFidelityReport:
    bounds_xyz: np.ndarray
    times: np.ndarray
    probe_points_xyz: np.ndarray
    reference_rms_speed: float
    metrics: tuple[ReplayFidelityMetrics, ...]
    probe_sha256: str
    metadata: Mapping[str, object]


def _validate_velocity_output(name: str, values: object, count: int) -> np.ndarray:
    out = np.asarray(values, dtype=float)
    if out.shape != (count, 3):
        raise ValueError(f"{name} must return shape ({count},3), got {out.shape}")
    if not np.all(np.isfinite(out)):
        raise ValueError(f"{name} returned nonfinite velocity")
    return out


def _validate_levels(levels: Sequence[ReplayLevel]) -> tuple[ReplayLevel, ...]:
    levels = tuple(levels)
    if not levels:
        raise ValueError("at least one replay level is required")
    seen: set[str] = set()
    previous: tuple[int, int, int] | None = None
    for level in levels:
        if not isinstance(level, ReplayLevel):
            raise TypeError("replay_levels must contain ReplayLevel objects")
        if not isinstance(level.label, str) or not level.label.strip():
            raise ValueError("every replay level needs a nonempty label")
        if level.label in seen:
            raise ValueError("replay level labels must be unique")
        seen.add(level.label)
        resolution = level.resolution_xyz
        if (
            not isinstance(resolution, tuple)
            or len(resolution) != 3
            or any(not isinstance(v, int) or isinstance(v, bool) or v < 2 for v in resolution)
        ):
            raise ValueError("resolution_xyz must be a tuple of three integers >=2")
        if previous is not None and not all(a < b for a, b in zip(previous, resolution)):
            raise ValueError("replay resolutions must increase strictly in every axis")
        if not callable(level.velocity):
            raise TypeError("replay level velocity must be callable")
        previous = resolution
    return levels


def _probe_digest(
    bounds: np.ndarray,
    times: np.ndarray,
    points: np.ndarray,
    seed: int,
) -> str:
    h = sha256()
    h.update(b"velocity_replay_fidelity_probe_v1")
    for name, arr in (("bounds", bounds), ("times", times), ("points", points)):
        canonical = np.ascontiguousarray(arr, dtype="<f8")
        h.update(name.encode("ascii"))
        h.update(np.asarray(canonical.shape, dtype="<i8").tobytes())
        h.update(canonical.tobytes())
    h.update(np.asarray([seed], dtype="<i8").tobytes())
    return h.hexdigest()


def audit_velocity_replay_fidelity(
    reference_velocity: VelocityCallable,
    replay_levels: Sequence[ReplayLevel],
    *,
    bounds_xyz: np.ndarray,
    times: np.ndarray,
    points_per_time: int,
    seed: int,
    candidate_sha256: str,
    provenance: str,
) -> VelocityReplayFidelityReport:
    """Compare frozen-grid replay callables to one reference velocity callable.

    The same scrambled Sobol spatial probes and declared times are sent to the
    reference and every replay. No registration, rescaling, component fitting,
    extrapolation, derivative estimation, thresholding, or automatic replay
    selection is performed.
    """
    if not callable(reference_velocity):
        raise TypeError("reference_velocity must be callable")
    levels = _validate_levels(replay_levels)

    bounds = np.asarray(bounds_xyz, dtype=float)
    if bounds.shape != (3, 2) or not np.all(np.isfinite(bounds)):
        raise ValueError("bounds_xyz must be a finite array with shape (3,2)")
    if not np.all(bounds[:, 1] > bounds[:, 0]):
        raise ValueError("each bounds_xyz upper bound must exceed its lower bound")

    ts = np.asarray(times, dtype=float)
    if ts.ndim != 1 or ts.size == 0 or not np.all(np.isfinite(ts)):
        raise ValueError("times must be a nonempty finite 1-D array")
    if ts.size > 1 and not np.all(np.diff(ts) > 0.0):
        raise ValueError("times must be strictly increasing")

    if (
        not isinstance(points_per_time, int)
        or isinstance(points_per_time, bool)
        or points_per_time < 16
        or points_per_time & (points_per_time - 1)
    ):
        raise ValueError("points_per_time must be a power-of-two integer >=16")
    if not isinstance(seed, int) or isinstance(seed, bool) or seed < 0:
        raise ValueError("seed must be a nonnegative integer")
    if not isinstance(candidate_sha256, str) or _SHA256_RE.fullmatch(candidate_sha256) is None:
        raise ValueError("candidate_sha256 must be 64 lowercase hexadecimal characters")
    if not isinstance(provenance, str) or not provenance.strip():
        raise ValueError("nonempty provenance is required")

    exponent = int(math.log2(points_per_time))
    sampler = qmc.Sobol(d=3, scramble=True, seed=seed)
    unit = np.asarray(sampler.random_base2(exponent), dtype=float)
    points = qmc.scale(unit, bounds[:, 0], bounds[:, 1])
    if points.shape != (points_per_time, 3) or not np.all(np.isfinite(points)):
        raise RuntimeError("Sobol probe generation failed")

    reference_slices: list[np.ndarray] = []
    for t in ts:
        reference_slices.append(
            _validate_velocity_output(
                "reference_velocity",
                reference_velocity(points[:, 0], points[:, 1], points[:, 2], float(t)),
                points_per_time,
            )
        )
    reference = np.stack(reference_slices, axis=0)
    reference_rms_speed = float(np.sqrt(np.mean(np.sum(reference * reference, axis=-1))))
    if not np.isfinite(reference_rms_speed) or reference_rms_speed <= 0.0:
        raise ValueError("reference velocity is exact-zero on all declared probes/times")

    reports: list[ReplayFidelityMetrics] = []
    for level in levels:
        replay_slices: list[np.ndarray] = []
        for t in ts:
            replay_slices.append(
                _validate_velocity_output(
                    f"replay[{level.label}]",
                    level.velocity(points[:, 0], points[:, 1], points[:, 2], float(t)),
                    points_per_time,
                )
            )
        replay = np.stack(replay_slices, axis=0)
        error = replay - reference
        vector_error = np.linalg.norm(error, axis=-1)
        vector_rms = float(np.sqrt(np.mean(vector_error * vector_error)))
        vector_max = float(np.max(vector_error))
        component_rms = np.sqrt(np.mean(error * error, axis=(0, 1)))
        speed_bias = float(np.mean(np.linalg.norm(replay, axis=-1) - np.linalg.norm(reference, axis=-1)))
        component_rms = np.asarray(component_rms, dtype=float)
        component_rms.setflags(write=False)
        reports.append(
            ReplayFidelityMetrics(
                label=level.label,
                resolution_xyz=level.resolution_xyz,
                vector_rms_error=vector_rms,
                vector_max_error=vector_max,
                normalized_vector_rms_error=vector_rms / reference_rms_speed,
                normalized_vector_max_error=vector_max / reference_rms_speed,
                component_rms_error=component_rms,
                mean_speed_bias=speed_bias,
            )
        )

    bounds = np.array(bounds, copy=True)
    ts = np.array(ts, copy=True)
    points = np.array(points, copy=True)
    for arr in (bounds, ts, points):
        arr.setflags(write=False)

    metadata = {
        "claim_scope": "frozen_grid_velocity_representation_fidelity_only",
        "sampling": "scipy.stats.qmc.Sobol scrambled; random_base2; fixed common spatial probes",
        "points_per_time": points_per_time,
        "seed": seed,
        "candidate_sha256": candidate_sha256,
        "provenance": provenance.strip(),
        "same_probes_for_all_levels": True,
        "probe_role": "representation_audit_only_not_optimization_or_pde_validation",
        "optimization_or_validation_data_reused": False,
        "registration_or_camera_fit": False,
        "rescaling_or_component_fit": False,
        "derivatives_computed": False,
        "pde_acceptance_evidence": False,
        "automatic_resolution_selection": False,
        "visualization_ready": False,
        "visual_correspondence_verified": False,
        "pde_validated": False,
        "paper_exact": False,
        "openai_field_identified": False,
        "blowup_proved": False,
    }
    return VelocityReplayFidelityReport(
        bounds_xyz=bounds,
        times=ts,
        probe_points_xyz=points,
        reference_rms_speed=reference_rms_speed,
        metrics=tuple(reports),
        probe_sha256=_probe_digest(bounds, ts, points, seed),
        metadata=metadata,
    )
