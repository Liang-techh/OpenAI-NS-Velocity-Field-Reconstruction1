"""Post-optimization capacity audit for existing even-radial Eq. (4.5) modes.

The bounded CR005 profile optimization moved only ``Phi(0,2)`` and ``F(0,2)``
close to their positive coefficient limits.  This module asks a narrower
representation question before growing the basis again: do the already-existing
symmetry-preserving ``Phi(1,0)`` and ``F(1,0)`` coefficients still provide
independent morphology control on the checked optimized candidate?

Every morphology quantity and validation residual is rebuilt from the public
``Eq45VelocityCandidate.at_points(...)->[u,v,w]`` path.  The registered
restricted force is frozen from the CR005 receipt and is never refit here.
Results are expression-capacity evidence only; this module does not modify the
saved candidate or promote visualization/PDE/OpenAI-field claims.
"""
from __future__ import annotations

from dataclasses import replace
import json
from typing import Any, Iterable

import numpy as np

from .constrained_eq45_candidate import Eq45VelocityCandidate
from .constrained_eq45_optimized_artifact import (
    OPTIMIZED_CANDIDATE_SHA256,
    default_optimization_path,
    load_checked_candidate,
)
from .constrained_force import RestrictedForce
from .constrained_pressure_compatibility import curl_pressure_free_momentum
from .constrained_restricted_force_curl_capacity import (
    REGISTERED_NU,
    _curl_force,
    _vector_rms,
    deterministic_probe_cloud,
)


TASK_ID = "CR003-EQ45-OPTIMIZED-EVEN-RADIAL-CAPACITY-015"
PARAMETER_LABELS = ("Phi(1,0)", "F(1,0)")
FEATURE_NAMES = (
    "radial_thickness_ratio",
    "axial_reach_ratio",
    "shoulder_swirl_energy_fraction",
    "shoulder_axial_poloidal_fraction",
)
REFERENCE_TIMES = (0.25, 0.50, 0.75)
OPTIMIZATION_SCHEMA = "eq45_profile_force_optimization_v1"
OPTIMIZATION_TASK_ID = "CR005-EQ45-PROFILE-FORCE-OPT-015"
TRUTH_BOUNDARY = {
    "velocity_changed": False,
    "production_coefficients_changed": False,
    "forcing_refit_on_validation": False,
    "new_basis_added": False,
    "visualization_ready": False,
    "visual_correspondence_verified": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


def _validated_sequence(name: str, values: Iterable[float]) -> tuple[float, ...]:
    out = tuple(float(value) for value in values)
    if not out or not all(np.isfinite(value) for value in out):
        raise ValueError(f"{name} must contain finite values")
    return out


def _with_existing_radial_delta(
    candidate: Eq45VelocityCandidate, *, channel: str, delta: float
) -> Eq45VelocityCandidate:
    """Return a temporary candidate with one existing ``(1,0)`` mode perturbed."""
    if not isinstance(candidate, Eq45VelocityCandidate):
        raise TypeError("candidate must be Eq45VelocityCandidate")
    if not np.isfinite(delta):
        raise ValueError("coefficient delta must be finite")

    basis = candidate.profile_basis
    try:
        index = basis.mode_indices.index((1, 0))
    except ValueError as exc:
        raise ValueError("candidate basis does not contain the required (1,0) mode") from exc

    if channel == "phi":
        values = list(basis.phi_coefficients)
        values[index] += float(delta)
        if abs(values[index]) > basis.coefficient_limit:
            raise ValueError("Phi(1,0) perturbation exceeds the registered coefficient bound")
        updated = replace(basis, phi_coefficients=tuple(values))
    elif channel == "swirl":
        values = list(basis.swirl_coefficients)
        values[index] += float(delta)
        if abs(values[index]) > basis.coefficient_limit:
            raise ValueError("F(1,0) perturbation exceeds the registered coefficient bound")
        updated = replace(basis, swirl_coefficients=tuple(values))
    else:
        raise ValueError("channel must be 'phi' or 'swirl'")
    return replace(candidate, profile_basis=updated)


def _ring_points(radius: float, z: float, n_angles: int) -> np.ndarray:
    angles = 2.0 * np.pi * np.arange(n_angles, dtype=float) / n_angles
    return np.column_stack(
        (
            radius * np.cos(angles),
            radius * np.sin(angles),
            np.full(n_angles, z, dtype=float),
        )
    )


def _rms(values: np.ndarray) -> float:
    return float(np.sqrt(np.mean(np.asarray(values, dtype=float) ** 2)))


def _cylindrical_components(points: np.ndarray, velocity: np.ndarray):
    radius = np.hypot(points[:, 0], points[:, 1])
    if np.any(radius <= 0.0):
        raise ValueError("morphology rings must stay away from the symmetry axis")
    er_x = points[:, 0] / radius
    er_y = points[:, 1] / radius
    radial = velocity[:, 0] * er_x + velocity[:, 1] * er_y
    swirl = -velocity[:, 0] * er_y + velocity[:, 1] * er_x
    return radial, swirl, velocity[:, 2]


def _morphology_features(
    candidate: Eq45VelocityCandidate,
    time: float,
    *,
    n_angles: int,
    inner_radius: float = 0.45,
    mid_radius: float = 0.65,
    outer_radius: float = 1.05,
    shoulder_z: float = 0.65,
) -> np.ndarray:
    """Return four dimensionless public-velocity morphology fingerprints."""
    inner = _ring_points(inner_radius, 0.0, n_angles)
    middle = _ring_points(mid_radius, 0.0, n_angles)
    outer = _ring_points(outer_radius, 0.0, n_angles)
    shoulder_plus = _ring_points(mid_radius, shoulder_z, n_angles)
    shoulder_minus = _ring_points(mid_radius, -shoulder_z, n_angles)
    shoulder = np.concatenate((shoulder_plus, shoulder_minus), axis=0)

    velocity_inner = np.asarray(candidate.at_points(inner, time), dtype=float)
    velocity_middle = np.asarray(candidate.at_points(middle, time), dtype=float)
    velocity_outer = np.asarray(candidate.at_points(outer, time), dtype=float)
    velocity_shoulder = np.concatenate(
        (
            np.asarray(candidate.at_points(shoulder_plus, time), dtype=float),
            np.asarray(candidate.at_points(shoulder_minus, time), dtype=float),
        ),
        axis=0,
    )
    for points, velocity in (
        (inner, velocity_inner),
        (middle, velocity_middle),
        (outer, velocity_outer),
        (shoulder, velocity_shoulder),
    ):
        if velocity.shape != points.shape or not np.all(np.isfinite(velocity)):
            raise RuntimeError("public velocity returned malformed morphology samples")

    speed_inner = np.linalg.norm(velocity_inner, axis=1)
    speed_middle = np.linalg.norm(velocity_middle, axis=1)
    speed_outer = np.linalg.norm(velocity_outer, axis=1)
    speed_shoulder = np.linalg.norm(velocity_shoulder, axis=1)
    inner_rms = _rms(speed_inner)
    middle_rms = _rms(speed_middle)
    if inner_rms <= 1e-14 or middle_rms <= 1e-14:
        raise ValueError("optimized candidate is numerically inactive on morphology rings")

    radial, swirl, axial = _cylindrical_components(shoulder, velocity_shoulder)
    total_energy = float(np.mean(speed_shoulder**2))
    poloidal_energy = float(np.mean(radial**2 + axial**2))
    if total_energy <= 1e-28 or poloidal_energy <= 1e-28:
        raise ValueError("optimized candidate is numerically inactive on shoulder rings")

    return np.asarray(
        (
            _rms(speed_outer) / inner_rms,
            _rms(speed_shoulder) / middle_rms,
            float(np.mean(swirl**2)) / total_energy,
            float(np.mean(axial**2)) / poloidal_energy,
        ),
        dtype=float,
    )


def _feature_vector(
    candidate: Eq45VelocityCandidate, times: tuple[float, ...], *, n_angles: int
) -> np.ndarray:
    return np.concatenate(
        [_morphology_features(candidate, time, n_angles=n_angles) for time in times]
    )


def _load_frozen_holdout_contract(candidate: Eq45VelocityCandidate) -> dict[str, Any]:
    receipt = json.loads(default_optimization_path().read_text(encoding="utf-8"))
    if receipt.get("schema") != OPTIMIZATION_SCHEMA or receipt.get("task_id") != OPTIMIZATION_TASK_ID:
        raise ValueError("unexpected CR005 optimization receipt")
    optimized_meta = receipt.get("optimized_candidate")
    contract = receipt.get("contract")
    fit = receipt.get("fit")
    holdout = receipt.get("holdout")
    if not isinstance(optimized_meta, dict) or not isinstance(contract, dict) or not isinstance(fit, dict):
        raise ValueError("optimization receipt is missing required sections")
    if not isinstance(holdout, list) or not holdout:
        raise ValueError("optimization receipt is missing holdout levels")
    if optimized_meta.get("sha256") != candidate.sha256:
        raise ValueError("candidate identity disagrees with the optimization receipt")
    if candidate.sha256 != OPTIMIZED_CANDIDATE_SHA256:
        raise ValueError("audit requires the checked optimized Eq45 candidate")
    if contract.get("training_holdout_separate") is not True:
        raise ValueError("optimization receipt lost train/holdout separation")

    steps = tuple(float(value) for value in contract.get("holdout_derivative_steps", ()))
    if len(steps) < 3 or any(step <= 0.0 for step in steps):
        raise ValueError("optimization receipt has an invalid derivative ladder")
    finest_step = steps[-1]
    matching = [row for row in holdout if float(row.get("step", -1.0)) == finest_step]
    if len(matching) != 1:
        raise ValueError("optimization receipt must contain exactly one finest holdout level")

    return {
        "holdout_seed": int(contract["holdout_seed"]),
        "sample_count": int(contract["sample_count"]),
        "finest_step": finest_step,
        "force": RestrictedForce(a=float(fit["force_a"]), c=float(fit["force_c"])),
        "receipt_finest_rms": float(matching[0]["optimized_curl_rms_after_force"]),
    }


def _frozen_force_holdout_rms(
    candidate: Eq45VelocityCandidate,
    points: np.ndarray,
    times: np.ndarray,
    *,
    step: float,
    force: RestrictedForce,
) -> float:
    target = curl_pressure_free_momentum(
        candidate,
        points,
        times,
        nu=REGISTERED_NU,
        spatial_step=step,
        time_step=step,
    )[1]
    after = target - _curl_force(force, points, times, step)
    return _vector_rms(after)


def audit_eq45_optimized_even_radial_capacity(
    candidate: Eq45VelocityCandidate | None = None,
    *,
    times: Iterable[float] = REFERENCE_TIMES,
    coefficient_steps: Iterable[float] = (0.04, 0.02),
    n_angles: int = 16,
) -> dict[str, Any]:
    """Audit two existing radial controls on the checked optimized candidate."""
    current = load_checked_candidate() if candidate is None else candidate
    if not isinstance(current, Eq45VelocityCandidate):
        raise TypeError("candidate must be Eq45VelocityCandidate")
    if current.sha256 != OPTIMIZED_CANDIDATE_SHA256:
        raise ValueError("audit requires the checked optimized Eq45 candidate identity")

    times_tuple = _validated_sequence("times", times)
    steps = _validated_sequence("coefficient_steps", coefficient_steps)
    if len(times_tuple) < 3 or any(a >= b for a, b in zip(times_tuple, times_tuple[1:])):
        raise ValueError("times must contain at least three strictly increasing values")
    if times_tuple[0] < current.time_start or times_tuple[-1] > current.time_end:
        raise ValueError("times must stay inside the candidate delivery interval")
    if len(steps) < 2 or any(step <= 0.0 for step in steps) or any(
        b >= a for a, b in zip(steps, steps[1:])
    ):
        raise ValueError("coefficient_steps must be positive and strictly decreasing")
    if not isinstance(n_angles, (int, np.integer)) or n_angles < 8:
        raise ValueError("n_angles must be an integer >= 8")

    basis = current.profile_basis
    mode_index = basis.mode_indices.index((1, 0))
    baseline_coefficients = np.asarray(
        (basis.phi_coefficients[mode_index], basis.swirl_coefficients[mode_index]), dtype=float
    )
    if np.any(np.abs(baseline_coefficients) + steps[0] > basis.coefficient_limit):
        raise ValueError("coefficient perturbation would exceed the registered bound")

    baseline_features = _feature_vector(current, times_tuple, n_angles=int(n_angles))
    matrices: list[np.ndarray] = []
    for step in steps:
        columns = []
        for channel in ("phi", "swirl"):
            plus = _with_existing_radial_delta(current, channel=channel, delta=step)
            minus = _with_existing_radial_delta(current, channel=channel, delta=-step)
            response = (
                _feature_vector(plus, times_tuple, n_angles=int(n_angles))
                - _feature_vector(minus, times_tuple, n_angles=int(n_angles))
            ) / (2.0 * step)
            columns.append(response)
        matrices.append(np.column_stack(columns))

    response = matrices[-1]
    singular_values = np.linalg.svd(response, compute_uv=False)
    tolerance = 1e-10 * singular_values[0]
    rank = int(np.count_nonzero(singular_values > tolerance))
    condition_number = (
        float(singular_values[0] / singular_values[-1]) if rank == response.shape[1] else float("inf")
    )
    response_norm = float(np.linalg.norm(response))
    if response_norm <= 0.0:
        raise RuntimeError("existing radial morphology response is numerically inactive")
    refinement_change = float(np.linalg.norm(matrices[0] - response) / response_norm)

    holdout_contract = _load_frozen_holdout_contract(current)
    holdout_points, holdout_times = deterministic_probe_cloud(
        holdout_contract["holdout_seed"], holdout_contract["sample_count"]
    )
    force = holdout_contract["force"]
    holdout_step = holdout_contract["finest_step"]
    baseline_holdout_rms = _frozen_force_holdout_rms(
        current,
        holdout_points,
        holdout_times,
        step=holdout_step,
        force=force,
    )
    if not np.isclose(
        baseline_holdout_rms,
        holdout_contract["receipt_finest_rms"],
        rtol=1e-10,
        atol=1e-10,
    ):
        raise RuntimeError("recomputed optimized holdout RMS disagrees with the checked CR005 receipt")

    baseline_velocity = np.asarray(current.at_points(holdout_points, holdout_times), dtype=float)
    baseline_velocity_rms = _vector_rms(baseline_velocity)
    finest_delta = steps[-1]
    residual_sensitivity: list[dict[str, float | str]] = []
    for column_index, channel in enumerate(("phi", "swirl")):
        plus = _with_existing_radial_delta(current, channel=channel, delta=finest_delta)
        minus = _with_existing_radial_delta(current, channel=channel, delta=-finest_delta)
        plus_rms = _frozen_force_holdout_rms(
            plus, holdout_points, holdout_times, step=holdout_step, force=force
        )
        minus_rms = _frozen_force_holdout_rms(
            minus, holdout_points, holdout_times, step=holdout_step, force=force
        )
        plus_velocity = np.asarray(plus.at_points(holdout_points, holdout_times), dtype=float)
        minus_velocity = np.asarray(minus.at_points(holdout_points, holdout_times), dtype=float)
        residual_sensitivity.append(
            {
                "parameter": PARAMETER_LABELS[column_index],
                "morphology_response_l2_per_unit": float(np.linalg.norm(response[:, column_index])),
                "plus_holdout_rms": float(plus_rms),
                "minus_holdout_rms": float(minus_rms),
                "best_fractional_holdout_change": float(
                    min(plus_rms, minus_rms) / baseline_holdout_rms - 1.0
                ),
                "worst_fractional_holdout_change": float(
                    max(plus_rms, minus_rms) / baseline_holdout_rms - 1.0
                ),
                "central_holdout_rms_derivative": float(
                    (plus_rms - minus_rms) / (2.0 * finest_delta)
                ),
                "plus_velocity_relative_rms_change": float(
                    _vector_rms(plus_velocity - baseline_velocity) / baseline_velocity_rms
                ),
                "minus_velocity_relative_rms_change": float(
                    _vector_rms(minus_velocity - baseline_velocity) / baseline_velocity_rms
                ),
            }
        )

    return {
        "task_id": TASK_ID,
        "claim_scope": "optimized_existing_even_radial_expression_capacity_only",
        "candidate_sha256": current.sha256,
        "times": list(times_tuple),
        "n_angles": int(n_angles),
        "coefficient_limit": float(basis.coefficient_limit),
        "coefficient_steps": list(steps),
        "parameter_labels": list(PARAMETER_LABELS),
        "baseline_parameter_values": baseline_coefficients.tolist(),
        "feature_names_per_time": list(FEATURE_NAMES),
        "baseline_feature_vector": baseline_features.tolist(),
        "finest_morphology_response_per_unit_coefficient": response.tolist(),
        "singular_values": singular_values.tolist(),
        "numerical_rank": rank,
        "condition_number": condition_number,
        "step_refinement_relative_change": refinement_change,
        "frozen_holdout": {
            "seed": holdout_contract["holdout_seed"],
            "sample_count": holdout_contract["sample_count"],
            "spatial_and_time_step": holdout_step,
            "force_a": float(force.a),
            "force_c": float(force.c),
            "baseline_curl_rms_after_force": float(baseline_holdout_rms),
            "force_refit": False,
        },
        "residual_sensitivity_at_finest_coefficient_step": residual_sensitivity,
        "recommended_use": (
            "Use these already-existing bounded radial controls before adding a new radial degree only if a public-observable morphology defect aligns with their response and independent validation tolerates the resulting residual change. This audit does not select a coefficient or promote the optimized artifact."
        ),
        "truth_boundary": dict(TRUTH_BOUNDARY),
    }
