"""Renderer-independent three-resolution morphology receipt for the strict-inner Kokuno candidate.

This Agent-2 increment consumes only the checksum-bound public ``velocity`` and
``vorticity`` surfaces from A2 #881.  It does not modify the candidate, introduce
a new derivative, or infer numerical targets from Kokuno/OpenAI figures.

The frozen protocol samples six source-safe PA.10 strict-inner rings through an
externally supplied coordinate mapper (pinned to exact Agent-1 #848 in dedicated
CI), evaluates angular resolutions 16/32/64, and records cylindrical velocity
and vorticity means/RMS plus low Fourier-mode magnitudes.  A phase-shifted fine
grid is an off-grid stability check.  This is a scoped morphology/quadrature
receipt, not whole-domain validation, source correspondence, or an NS residual.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import inspect
import json
import math
from typing import Any, Callable, Mapping, Protocol, runtime_checkable

import numpy as np

TASK = "KOKUNO-A2-STRICT-INNER-VELOCITY-VORTICITY-MORPHOLOGY-068"
SCHEMA = "kokuno-a2-strict-inner-velocity-vorticity-morphology-v1"
PARENT_AGENT2_PR = 881
PARENT_AGENT2_HEAD = "ded8c88dd03b784368cf1a30601faf76a3682c64"
PARENT_AGENT2_SOURCE_BLOB = "9dd3b7400fea8f2d4462cda86561103e29a1d053"
AGENT1_RUNTIME_PR = 848
AGENT1_RUNTIME_HEAD = "300919075c07396e4e466f553be78eaae0707df2"
SOURCE_CORRECTED_READER_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_CORRECTED_READER_DATE = "2026-09-09"

ANGULAR_RESOLUTIONS = (16, 32, 64)
OFFGRID_PHASE_FRACTION = 0.371
FOURIER_MODES = (1, 2)

# (source X, eta, t). These are repository validation choices, not source data.
SOURCE_RINGS = (
    (0.180, -0.180, 0.470),
    (0.180, +0.180, 0.530),
    (0.225, -0.120, 0.490),
    (0.225, +0.120, 0.510),
    (0.270, -0.060, 0.480),
    (0.270, +0.060, 0.520),
)

COARSE_MEDIUM_RELATIVE_RMS_MAX = 3.0e-2
MEDIUM_FINE_RELATIVE_RMS_MAX = 1.0e-2
OFFGRID_FINE_RELATIVE_RMS_MAX = 1.0e-2
VELOCITY_RMS_MIN = 1.0e-8
VORTICITY_RMS_MIN = 1.0e-10

_TRUTH_BOUNDARY = {
    "strict_inner_candidate_identity_preserved": True,
    "strict_inner_three_resolution_morphology_protocol_materialized": True,
    "renderer_or_camera_used": False,
    "source_numeric_targets_used": False,
    "pixel_loss_used": False,
    "new_velocity_or_vorticity_derivative_realization_introduced": False,
    "candidate_parameters_changed": False,
    "whole_domain_velocity_morphology_verified": False,
    "whole_domain_vorticity_morphology_verified": False,
    "visual_correspondence_verified": False,
    "global_leading_velocity_materialized": False,
    "outer_join_materialized": False,
    "matched_pressure_included": False,
    "restricted_forcing_included": False,
    "agent3_correction_velocity_included": False,
    "complete_kokuno_candidate_assembled": False,
    "complete_ns_momentum_residual_formed": False,
    "heldout_ns_residual_assessed": False,
    "st006_same_protocol_comparison_valid": False,
    "residual_reduction_claimed": False,
    "paper_exact": False,
    "pde_validated": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


@runtime_checkable
class VelocityVorticityCandidate(Protocol):
    @property
    def velocity_candidate_sha256(self) -> str: ...

    @property
    def vorticity_sha256(self) -> str: ...

    def velocity(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray: ...

    def vorticity(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray: ...


CoordinateMapper = Callable[
    [np.ndarray, np.ndarray, np.ndarray, np.ndarray],
    Mapping[str, Any],
]


def _canonical_json(payload: Mapping[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _sha256_payload(payload: Mapping[str, Any]) -> str:
    return hashlib.sha256(_canonical_json(payload).encode("utf-8")).hexdigest()


def _validate_sha256(value: Any, name: str) -> str:
    text = str(value)
    if len(text) != 64 or any(c not in "0123456789abcdef" for c in text):
        raise ValueError(f"{name} must be one lowercase hexadecimal sha256")
    return text


def _rms(values: Any) -> float:
    arr = np.asarray(values, dtype=float)
    return float(np.sqrt(np.mean(arr * arr))) if arr.size else 0.0


def _relative_rms(a: Any, b: Any) -> float:
    aa = np.asarray(a, dtype=float)
    bb = np.asarray(b, dtype=float)
    return _rms(aa - bb) / max(_rms(bb), 1.0e-30)


def frozen_protocol() -> dict[str, Any]:
    payload: dict[str, Any] = {
        "schema": SCHEMA,
        "parent_agent2_pr": PARENT_AGENT2_PR,
        "parent_agent2_head": PARENT_AGENT2_HEAD,
        "parent_agent2_source_blob": PARENT_AGENT2_SOURCE_BLOB,
        "agent1_runtime_pr": AGENT1_RUNTIME_PR,
        "agent1_runtime_head": AGENT1_RUNTIME_HEAD,
        "source_provenance": {
            "corrected_reader_commit": SOURCE_CORRECTED_READER_COMMIT,
            "corrected_reader_date": SOURCE_CORRECTED_READER_DATE,
        },
        "source_rings": [list(v) for v in SOURCE_RINGS],
        "angular_resolutions": list(ANGULAR_RESOLUTIONS),
        "offgrid_phase_fraction": OFFGRID_PHASE_FRACTION,
        "fourier_modes": list(FOURIER_MODES),
        "gates": {
            "coarse_medium_relative_rms_max": COARSE_MEDIUM_RELATIVE_RMS_MAX,
            "medium_fine_relative_rms_max": MEDIUM_FINE_RELATIVE_RMS_MAX,
            "offgrid_fine_relative_rms_max": OFFGRID_FINE_RELATIVE_RMS_MAX,
            "velocity_rms_min": VELOCITY_RMS_MIN,
            "vorticity_rms_min": VORTICITY_RMS_MIN,
        },
        "semantics": {
            "geometry": "strict-inner source-safe rings mapped by exact Agent-1 coordinate runtime",
            "resolution_axis": "azimuthal ring quadrature only",
            "morphology": "cylindrical means/RMS and |m|=1,2 Fourier magnitudes",
            "offgrid_check": "phase-shifted finest azimuthal grid",
            "comparison_scope": "same checksum-bound candidate identity",
        },
        "truth_boundary": dict(_TRUTH_BOUNDARY),
    }
    payload["protocol_sha256"] = _sha256_payload(payload)
    return payload


def _map_ring_points(
    coordinate_mapper: CoordinateMapper,
    n_theta: int,
    phase_fraction: float = 0.0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    if int(n_theta) != n_theta or n_theta < 8:
        raise ValueError("n_theta must be an integer >= 8")
    n_theta = int(n_theta)
    base = 2.0 * math.pi * np.arange(n_theta, dtype=float) / n_theta
    theta_1d = base + (2.0 * math.pi / n_theta) * float(phase_fraction)

    X = np.asarray([row[0] for row in SOURCE_RINGS], dtype=float)[:, None]
    eta = np.asarray([row[1] for row in SOURCE_RINGS], dtype=float)[:, None]
    tt = np.asarray([row[2] for row in SOURCE_RINGS], dtype=float)[:, None]
    theta = np.broadcast_to(theta_1d[None, :], (len(SOURCE_RINGS), n_theta))
    X = np.broadcast_to(X, theta.shape)
    eta = np.broadcast_to(eta, theta.shape)
    tt = np.broadcast_to(tt, theta.shape)

    raw = coordinate_mapper(X, eta, tt, theta)
    if not isinstance(raw, Mapping):
        raise TypeError("coordinate_mapper must return a mapping with x/y/z/t")
    values = []
    for name in ("x", "y", "z", "t"):
        arr = np.asarray(raw.get(name), dtype=float)
        if arr.shape != theta.shape:
            raise ValueError(
                f"coordinate_mapper {name} must have shape {theta.shape}, got {arr.shape}"
            )
        if not np.all(np.isfinite(arr)):
            raise ValueError(f"coordinate_mapper {name} must be finite")
        values.append(arr)
    return values[0], values[1], values[2], values[3], theta


def _cylindrical_components(vector: np.ndarray, x: np.ndarray, y: np.ndarray) -> np.ndarray:
    vec = np.asarray(vector, dtype=float)
    if vec.shape != x.shape + (3,):
        raise ValueError(f"vector must have shape {x.shape + (3,)}, got {vec.shape}")
    if not np.all(np.isfinite(vec)):
        raise ValueError("candidate vector values must be finite")
    phi = np.arctan2(y, x)
    c = np.cos(phi)
    s = np.sin(phi)
    return np.stack(
        (
            vec[..., 0] * c + vec[..., 1] * s,
            -vec[..., 0] * s + vec[..., 1] * c,
            vec[..., 2],
        ),
        axis=-1,
    )


def _field_descriptor(cyl: np.ndarray, phi: np.ndarray) -> tuple[dict[str, Any], np.ndarray]:
    mean = np.mean(cyl, axis=1)
    rms = np.sqrt(np.mean(cyl * cyl, axis=1))
    norm_rms = np.sqrt(np.mean(np.sum(cyl * cyl, axis=-1), axis=1))

    harmonics: dict[str, Any] = {}
    vectors = [mean.reshape(-1), rms.reshape(-1), norm_rms.reshape(-1)]
    for mode in FOURIER_MODES:
        phase = np.exp(-1j * mode * phi)[..., None]
        amp = np.abs(np.mean(cyl * phase, axis=1))
        harmonics[f"m{mode}_abs"] = amp.tolist()
        vectors.append(amp.reshape(-1))

    report = {
        "component_mean": mean.tolist(),
        "component_rms": rms.tolist(),
        "vector_norm_rms": norm_rms.tolist(),
        "harmonic_magnitudes": harmonics,
    }
    return report, np.concatenate(vectors)


def _evaluate_level(
    candidate: VelocityVorticityCandidate,
    coordinate_mapper: CoordinateMapper,
    n_theta: int,
    phase_fraction: float = 0.0,
) -> dict[str, Any]:
    x, y, z, t, _source_theta = _map_ring_points(
        coordinate_mapper, n_theta, phase_fraction
    )
    velocity = np.asarray(candidate.velocity(x, y, z, t), dtype=float)
    vorticity = np.asarray(candidate.vorticity(x, y, z, t), dtype=float)
    expected = x.shape + (3,)
    if velocity.shape != expected:
        raise ValueError(f"candidate velocity must have shape {expected}, got {velocity.shape}")
    if vorticity.shape != expected:
        raise ValueError(
            f"candidate vorticity must have shape {expected}, got {vorticity.shape}"
        )
    if not np.all(np.isfinite(velocity)) or not np.all(np.isfinite(vorticity)):
        raise ValueError("candidate velocity/vorticity must be finite")

    phi = np.arctan2(y, x)
    velocity_cyl = _cylindrical_components(velocity, x, y)
    vorticity_cyl = _cylindrical_components(vorticity, x, y)
    velocity_report, velocity_vector = _field_descriptor(velocity_cyl, phi)
    vorticity_report, vorticity_vector = _field_descriptor(vorticity_cyl, phi)

    radius = np.sqrt(x * x + y * y)
    geometry = {
        "radius_mean": np.mean(radius, axis=1).tolist(),
        "radius_span": (np.max(radius, axis=1) - np.min(radius, axis=1)).tolist(),
        "z_mean": np.mean(z, axis=1).tolist(),
        "t_mean": np.mean(t, axis=1).tolist(),
    }
    return {
        "n_theta": int(n_theta),
        "phase_fraction": float(phase_fraction),
        "velocity": velocity_report,
        "vorticity": vorticity_report,
        "velocity_vector": velocity_vector,
        "vorticity_vector": vorticity_vector,
        "velocity_rms": _rms(velocity),
        "vorticity_rms": _rms(vorticity),
        "geometry": geometry,
    }


def materialize_strict_inner_morphology(
    candidate: VelocityVorticityCandidate,
    coordinate_mapper: CoordinateMapper,
) -> dict[str, Any]:
    """Evaluate the frozen strict-inner three-resolution morphology protocol."""
    if not callable(getattr(candidate, "velocity", None)):
        raise TypeError("candidate must expose velocity(x,y,z,t)")
    if not callable(getattr(candidate, "vorticity", None)):
        raise TypeError("candidate must expose vorticity(x,y,z,t)")
    if not callable(coordinate_mapper):
        raise TypeError("coordinate_mapper must be callable")

    velocity_sha = _validate_sha256(
        getattr(candidate, "velocity_candidate_sha256", ""),
        "candidate.velocity_candidate_sha256",
    )
    vorticity_sha = _validate_sha256(
        getattr(candidate, "vorticity_sha256", ""),
        "candidate.vorticity_sha256",
    )

    levels: dict[str, dict[str, Any]] = {}
    raw_levels = []
    for n_theta in ANGULAR_RESOLUTIONS:
        level = _evaluate_level(candidate, coordinate_mapper, n_theta)
        raw_levels.append(level)
        levels[str(n_theta)] = {
            key: value
            for key, value in level.items()
            if key not in ("velocity_vector", "vorticity_vector")
        }

    coarse, medium, fine = raw_levels
    shifted = _evaluate_level(
        candidate,
        coordinate_mapper,
        ANGULAR_RESOLUTIONS[-1],
        OFFGRID_PHASE_FRACTION,
    )

    velocity_cm = _relative_rms(
        coarse["velocity_vector"], medium["velocity_vector"]
    )
    velocity_mf = _relative_rms(medium["velocity_vector"], fine["velocity_vector"])
    velocity_offgrid = _relative_rms(
        shifted["velocity_vector"], fine["velocity_vector"]
    )
    vorticity_cm = _relative_rms(
        coarse["vorticity_vector"], medium["vorticity_vector"]
    )
    vorticity_mf = _relative_rms(
        medium["vorticity_vector"], fine["vorticity_vector"]
    )
    vorticity_offgrid = _relative_rms(
        shifted["vorticity_vector"], fine["vorticity_vector"]
    )

    gate_values = {
        "velocity_coarse_medium_relative_rms": velocity_cm,
        "velocity_medium_fine_relative_rms": velocity_mf,
        "velocity_offgrid_fine_relative_rms": velocity_offgrid,
        "vorticity_coarse_medium_relative_rms": vorticity_cm,
        "vorticity_medium_fine_relative_rms": vorticity_mf,
        "vorticity_offgrid_fine_relative_rms": vorticity_offgrid,
        "fine_velocity_rms": fine["velocity_rms"],
        "fine_vorticity_rms": fine["vorticity_rms"],
    }
    gate_checks = {
        "velocity_coarse_medium": velocity_cm <= COARSE_MEDIUM_RELATIVE_RMS_MAX,
        "velocity_medium_fine": velocity_mf <= MEDIUM_FINE_RELATIVE_RMS_MAX,
        "velocity_offgrid": velocity_offgrid <= OFFGRID_FINE_RELATIVE_RMS_MAX,
        "vorticity_coarse_medium": vorticity_cm <= COARSE_MEDIUM_RELATIVE_RMS_MAX,
        "vorticity_medium_fine": vorticity_mf <= MEDIUM_FINE_RELATIVE_RMS_MAX,
        "vorticity_offgrid": vorticity_offgrid <= OFFGRID_FINE_RELATIVE_RMS_MAX,
        "velocity_nontrivial": fine["velocity_rms"] >= VELOCITY_RMS_MIN,
        "vorticity_nontrivial": fine["vorticity_rms"] >= VORTICITY_RMS_MIN,
    }
    passed = bool(all(gate_checks.values()))

    payload: dict[str, Any] = {
        "schema": SCHEMA,
        "task": TASK,
        "candidate_identity": {
            "velocity_candidate_sha256": velocity_sha,
            "vorticity_sha256": vorticity_sha,
        },
        "protocol": frozen_protocol(),
        "levels": levels,
        "shifted_fine": {
            key: value
            for key, value in shifted.items()
            if key not in ("velocity_vector", "vorticity_vector")
        },
        "gate_values": gate_values,
        "gate_checks": gate_checks,
        "strict_inner_three_resolution_morphology_gate_passed": passed,
        "truth_boundary": dict(_TRUTH_BOUNDARY),
    }
    payload["receipt_sha256"] = _sha256_payload(payload)
    return payload


def public_contract() -> dict[str, Any]:
    signature = inspect.signature(materialize_strict_inner_morphology)
    return {
        "task": TASK,
        "schema": SCHEMA,
        "materialize_signature": str(signature),
        "scientific_knobs_exposed": False,
        "protocol": frozen_protocol(),
        "truth_boundary": dict(_TRUTH_BOUNDARY),
    }
