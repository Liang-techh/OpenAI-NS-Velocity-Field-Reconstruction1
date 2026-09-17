"""Target-free interior axial-basis capacity audit for the bipolar Eq45 candidate.

Agent 4's physical-taper screen found that moving the exterior axial taper does
not move the main vorticity-core geometry. This module therefore re-tests one
old basis-growth idea, ``Phi(0,4)``, on the *current* source-aligned bipolar
supported candidate. The old screen predates the active ``Phi(0,1)`` central
poloidal seed, so its conclusion cannot simply be inherited.

The audit is diagnostic only. It embeds the current basis into an eta-degree-4
tensor basis without changing any existing coefficient, measures centered
public-velocity Jacobians for the existing ``Phi(1,2)`` direction and the new
``Phi(0,4)`` direction, checks their z-parity, and records one small whole-domain
vorticity morphology screen. It does not fit OpenAI pixels, force, pressure, or
a PDE residual, and it does not modify the packaged/canonical velocity field.
"""
from __future__ import annotations

from dataclasses import replace
from typing import Any, Iterable

import numpy as np

from .constrained_eq45_bipolar_seed import bipolar_seed
from .constrained_eq45_profile_basis import Eq45CompactProfileBasis
from .constrained_eq45_supported_candidate import Eq45SupportedVelocityCandidate
from .constrained_eq45_supported_vorticity_envelope import structured_vorticity_magnitude
from .eq45_supported_delivery import default_field


TASK_ID = "CR003-BIPOLAR-PHI04-INTERIOR-AXIAL-CAPACITY-034"
PHI12 = (1, 2)
PHI04 = (0, 4)
TRUTH_BOUNDARY = {
    "velocity_changed": False,
    "candidate_artifact_changed": False,
    "force_or_pressure_fitted": False,
    "held_out_pde_residual_evaluated": False,
    "visualization_ready": False,
    "visual_correspondence_verified": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


def _source_bipolar_field() -> Eq45SupportedVelocityCandidate:
    return bipolar_seed(default_field().candidate)


def _extend_eta4(basis: Eq45CompactProfileBasis) -> Eq45CompactProfileBasis:
    """Embed the current basis exactly into eta degree four with zero new modes."""
    if not isinstance(basis, Eq45CompactProfileBasis):
        raise TypeError("basis must be Eq45CompactProfileBasis")
    if basis.eta_degree > 4:
        raise ValueError("audit only supports source bases with eta_degree <= 4")

    modes = tuple((i, j) for i in range(basis.radial_degree + 1) for j in range(5))
    index = {mode: k for k, mode in enumerate(modes)}
    phi = np.zeros(len(modes), dtype=float)
    swirl = np.zeros(len(modes), dtype=float)
    for source_index, mode in enumerate(basis.mode_indices):
        destination = index[mode]
        phi[destination] = basis.phi_coefficients[source_index]
        swirl[destination] = basis.swirl_coefficients[source_index]

    return Eq45CompactProfileBasis(
        radial_degree=basis.radial_degree,
        eta_degree=4,
        phi_coefficients=tuple(float(v) for v in phi),
        swirl_coefficients=tuple(float(v) for v in swirl),
        x_cut=basis.x_cut,
        eta_cut=basis.eta_cut,
        cutoff_power=basis.cutoff_power,
        coefficient_limit=basis.coefficient_limit,
    )


def _embedded_bipolar_field() -> Eq45SupportedVelocityCandidate:
    base = _source_bipolar_field()
    extended = _extend_eta4(base.parent.profile_basis)
    return replace(base, parent=replace(base.parent, profile_basis=extended))


def _with_phi_delta(
    field: Eq45SupportedVelocityCandidate,
    mode: tuple[int, int],
    delta: float,
) -> Eq45SupportedVelocityCandidate:
    if not np.isfinite(delta):
        raise ValueError("delta must be finite")
    basis = field.parent.profile_basis
    try:
        index = basis.mode_indices.index(mode)
    except ValueError as exc:
        raise ValueError(f"basis does not contain mode {mode}") from exc
    values = list(basis.phi_coefficients)
    values[index] += float(delta)
    updated = replace(basis, phi_coefficients=tuple(values))
    return replace(field, parent=replace(field.parent, profile_basis=updated))


def _ring(radius: float, z: float, count: int = 12) -> np.ndarray:
    theta = 2.0 * np.pi * np.arange(count, dtype=float) / count
    return np.column_stack(
        (radius * np.cos(theta), radius * np.sin(theta), np.full(count, z))
    )


def _probe_sets() -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    midplane = np.concatenate([_ring(r, 0.0) for r in (0.40, 0.80, 1.20)], axis=0)
    shoulder = np.concatenate(
        [_ring(r, z) for r in (0.50, 0.90) for z in (-0.60, -0.40, 0.40, 0.60)],
        axis=0,
    )
    axial_tip = np.concatenate(
        [_ring(r, z) for r in (0.25, 0.50) for z in (-1.30, -1.00, 1.00, 1.30)],
        axis=0,
    )
    combined = np.concatenate((midplane, shoulder, axial_tip), axis=0)
    return midplane, shoulder, axial_tip, combined


def _response(
    field: Eq45SupportedVelocityCandidate,
    mode: tuple[int, int],
    step: float,
    points: np.ndarray,
    times: tuple[float, ...],
) -> np.ndarray:
    plus = _with_phi_delta(field, mode, step)
    minus = _with_phi_delta(field, mode, -step)
    chunks = []
    for time in times:
        up = np.asarray(plus.at_points(points, time), dtype=float)
        um = np.asarray(minus.at_points(points, time), dtype=float)
        chunks.append(((up - um) / (2.0 * step)).reshape(-1))
    out = np.concatenate(chunks)
    if not np.all(np.isfinite(out)):
        raise RuntimeError("public-velocity response became nonfinite")
    return out


def _rms(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    return float(np.sqrt(np.mean(values * values)))


def _w_parity_response(
    field: Eq45SupportedVelocityCandidate,
    mode: tuple[int, int],
    step: float,
    times: tuple[float, ...],
) -> dict[str, float]:
    plus_field = _with_phi_delta(field, mode, step)
    minus_field = _with_phi_delta(field, mode, -step)
    even_samples: list[np.ndarray] = []
    odd_samples: list[np.ndarray] = []
    for time in times:
        for radius in (0.25, 0.50):
            for z in (0.60, 1.00):
                p = _ring(radius, z)
                m = p.copy()
                m[:, 2] *= -1.0
                rp = (plus_field.at_points(p, time) - minus_field.at_points(p, time)) / (2.0 * step)
                rm = (plus_field.at_points(m, time) - minus_field.at_points(m, time)) / (2.0 * step)
                even_samples.append(0.5 * (rp[:, 2] + rm[:, 2]))
                odd_samples.append(0.5 * (rp[:, 2] - rm[:, 2]))
    even_rms = _rms(np.concatenate(even_samples))
    odd_rms = _rms(np.concatenate(odd_samples))
    return {
        "axial_velocity_even_z_response_rms": even_rms,
        "axial_velocity_odd_z_response_rms": odd_rms,
        "even_fraction_of_axial_parity_response": even_rms / max(even_rms + odd_rms, 1e-300),
    }


def _vorticity_metrics(
    field: Eq45SupportedVelocityCandidate,
    *,
    time: float,
    resolution: int,
) -> dict[str, float]:
    axis = np.linspace(-2.0, 2.0, int(resolution))
    magnitude, radius, abs_z = structured_vorticity_magnitude(field, axis, time)
    weights = magnitude * magnitude
    total = float(np.sum(weights))
    if not np.isfinite(total) or total <= np.finfo(float).tiny:
        raise RuntimeError("vorticity weight is numerically trivial")
    radial_rms = float(np.sqrt(np.sum(weights * radius * radius) / total))
    axial_rms = float(np.sqrt(np.sum(weights * abs_z * abs_z) / total))
    plateau = (radius <= 1.5) & (abs_z <= 1.5)
    threshold = 0.25 * float(np.max(magnitude[plateau]))
    active = magnitude >= threshold
    if not np.any(active):
        raise RuntimeError("25%-core vorticity set is empty")
    return {
        "weighted_radial_rms": radial_rms,
        "weighted_axial_rms": axial_rms,
        "weighted_aspect": axial_rms / max(radial_rms, 1e-300),
        "core_radial_q99": float(np.quantile(radius[active], 0.99)),
        "core_axial_q99": float(np.quantile(abs_z[active], 0.99)),
    }


def audit_bipolar_phi04_capacity(
    *,
    times: Iterable[float] = (0.25, 0.50, 0.75),
    coefficient_steps: Iterable[float] = (0.02, 0.01),
    morphology_step: float = 0.25,
    morphology_time: float = 0.50,
    morphology_resolution: int = 33,
) -> dict[str, Any]:
    """Measure whether one interior ``Phi(0,4)`` mode adds useful axial control."""
    time_values = tuple(float(v) for v in times)
    steps = tuple(float(v) for v in coefficient_steps)
    if not time_values or not all(np.isfinite(v) for v in time_values):
        raise ValueError("times must be nonempty and finite")
    if len(steps) != 2 or not (steps[0] > steps[1] > 0.0):
        raise ValueError("coefficient_steps must contain two positive decreasing steps")
    if not np.isfinite(morphology_step) or morphology_step <= 0.0:
        raise ValueError("morphology_step must be positive and finite")
    if not isinstance(morphology_resolution, (int, np.integer)) or morphology_resolution < 17:
        raise ValueError("morphology_resolution must be an integer >= 17")

    raw_base = _source_bipolar_field()
    field = _embedded_bipolar_field()
    if min(time_values) < field.time_start or max(time_values) > field.time_end:
        raise ValueError("times must lie inside the delivery interval")
    limit = float(field.parent.profile_basis.coefficient_limit)
    if max(steps[0], morphology_step) > limit:
        raise ValueError("screen perturbation exceeds coefficient_limit")

    midplane_points, shoulder_points, axial_tip_points, combined_points = _probe_sets()
    embedding_errors = []
    for time in time_values:
        a = np.asarray(raw_base.at_points(combined_points, time), dtype=float)
        b = np.asarray(field.at_points(combined_points, time), dtype=float)
        embedding_errors.append(float(np.max(np.abs(a - b))))

    matrices = []
    for step in steps:
        columns = [
            _response(field, PHI12, step, combined_points, time_values),
            _response(field, PHI04, step, combined_points, time_values),
        ]
        matrices.append(np.column_stack(columns))
    response = matrices[-1]
    singular_values = np.linalg.svd(response, compute_uv=False)
    tolerance = max(response.shape) * np.finfo(float).eps * singular_values[0]
    rank = int(np.sum(singular_values > tolerance))
    condition = (
        float(singular_values[0] / singular_values[-1])
        if rank == 2 and singular_values[-1] > 0.0
        else float("inf")
    )
    refinement = float(
        np.linalg.norm(matrices[0] - response) / max(np.linalg.norm(response), 1e-300)
    )

    phi12 = response[:, 0]
    phi04 = response[:, 1]
    projection = phi12 * (float(np.dot(phi12, phi04)) / max(float(np.dot(phi12, phi12)), 1e-300))
    novelty = float(np.linalg.norm(phi04 - projection) / max(np.linalg.norm(phi04), 1e-300))
    cosine = float(np.dot(phi12, phi04) / max(np.linalg.norm(phi12) * np.linalg.norm(phi04), 1e-300))

    finest = steps[-1]
    local: dict[str, Any] = {}
    for label, mode in (("Phi(1,2)", PHI12), ("Phi(0,4)", PHI04)):
        midplane_rms = _rms(_response(field, mode, finest, midplane_points, time_values))
        shoulder_rms = _rms(_response(field, mode, finest, shoulder_points, time_values))
        axial_rms = _rms(_response(field, mode, finest, axial_tip_points, time_values))
        local[label] = {
            "midplane_response_rms_per_unit_coefficient": midplane_rms,
            "shoulder_response_rms_per_unit_coefficient": shoulder_rms,
            "axial_tip_response_rms_per_unit_coefficient": axial_rms,
            "axial_tip_to_shoulder_response_ratio": axial_rms / max(shoulder_rms, 1e-300),
            "parity": _w_parity_response(field, mode, finest, time_values),
        }

    morphology_fields = {
        "minus": _with_phi_delta(field, PHI04, -morphology_step),
        "base": field,
        "plus": _with_phi_delta(field, PHI04, morphology_step),
    }
    morphology = {
        label: _vorticity_metrics(candidate, time=float(morphology_time), resolution=int(morphology_resolution))
        for label, candidate in morphology_fields.items()
    }
    centered_morphology_derivative = {
        key: float((morphology["plus"][key] - morphology["minus"][key]) / (2.0 * morphology_step))
        for key in morphology["base"]
    }

    base_parameter_count = raw_base.parent.profile_basis.parameter_count
    extended_parameter_count = field.parent.profile_basis.parameter_count
    return {
        "task_id": TASK_ID,
        "claim_scope": "local_visualization_expression_capacity_only",
        "active_candidate": "bipolar_seed(default_field().candidate)",
        "base_profile_parameter_count": int(base_parameter_count),
        "eta4_tensor_parameter_count": int(extended_parameter_count),
        "newly_active_parameters_in_proposed_minimal_extension": 1,
        "new_mode": "Phi(0,4)",
        "reference_mode": "Phi(1,2)",
        "times": list(time_values),
        "coefficient_steps": list(steps),
        "coefficient_limit": limit,
        "zero_padding_embedding_max_abs_velocity_error": max(embedding_errors),
        "singular_values": singular_values.tolist(),
        "numerical_rank": rank,
        "condition_number": condition,
        "phi04_novelty_fraction_vs_phi12": novelty,
        "phi12_phi04_response_cosine": cosine,
        "finite_difference_refinement_relative_change": refinement,
        "local_response": local,
        "morphology_step": float(morphology_step),
        "morphology_time": float(morphology_time),
        "morphology_resolution": int(morphology_resolution),
        "whole_domain_vorticity_morphology": morphology,
        "whole_domain_centered_derivative_per_unit_phi04": centered_morphology_derivative,
        "held_out_pde_residual": {
            "evaluated": False,
            "reason": "capacity-only lane; no force/pressure fit or residual retuning is permitted in this increment",
        },
        "truth_boundary": dict(TRUTH_BOUNDARY),
    }
