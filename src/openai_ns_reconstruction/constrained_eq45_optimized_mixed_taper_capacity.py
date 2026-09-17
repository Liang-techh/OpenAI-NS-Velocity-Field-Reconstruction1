"""Post-optimization capacity audit for existing mixed radial/axial Eq. (4.5) modes.

The checked CR005 optimized candidate still contains the bounded ``Phi(1,2)``
and ``F(1,2)`` coefficients well inside their registered limits.  This module
asks one narrow representation question before adding any new basis channel:
do those already-present mixed radial×eta-even controls add a genuinely new
morphology direction beyond the existing ``(1,0)`` radial controls audited by
Constrained Agent 7?

All morphology and validation quantities are rebuilt from the public
``Eq45VelocityCandidate.at_points(...)->[u,v,w]`` path.  The frozen restricted
force from the CR005 receipt is never refit.  The saved optimized candidate is
not modified and no visualization/PDE/OpenAI-field claim is promoted.
"""
from __future__ import annotations

from dataclasses import replace
from typing import Any, Iterable

import numpy as np

from .constrained_eq45_candidate import Eq45VelocityCandidate
from .constrained_eq45_optimized_artifact import (
    OPTIMIZED_CANDIDATE_SHA256,
    load_checked_candidate,
)
from .constrained_eq45_optimized_even_radial_capacity import (
    FEATURE_NAMES,
    REFERENCE_TIMES,
    TRUTH_BOUNDARY as RADIAL_TRUTH_BOUNDARY,
    _feature_vector,
    _frozen_force_holdout_rms,
    _load_frozen_holdout_contract,
    _validated_sequence,
    audit_eq45_optimized_even_radial_capacity,
)
from .constrained_restricted_force_curl_capacity import (
    _vector_rms,
    deterministic_probe_cloud,
)


TASK_ID = "CR003-EQ45-OPTIMIZED-MIXED-TAPER-CAPACITY-016"
PARAMETER_LABELS = ("Phi(1,2)", "F(1,2)")
MODE = (1, 2)
TRUTH_BOUNDARY = dict(RADIAL_TRUTH_BOUNDARY)
TRUTH_BOUNDARY.update(
    {
        "mixed_mode_activated_in_production": False,
        "existing_candidate_artifact_modified": False,
    }
)


def _with_existing_mixed_delta(
    candidate: Eq45VelocityCandidate, *, channel: str, delta: float
) -> Eq45VelocityCandidate:
    """Return a temporary candidate with one existing ``(1,2)`` mode perturbed."""
    if not isinstance(candidate, Eq45VelocityCandidate):
        raise TypeError("candidate must be Eq45VelocityCandidate")
    if not np.isfinite(delta):
        raise ValueError("coefficient delta must be finite")

    basis = candidate.profile_basis
    try:
        index = basis.mode_indices.index(MODE)
    except ValueError as exc:
        raise ValueError("candidate basis does not contain the required (1,2) mode") from exc

    if channel == "phi":
        values = list(basis.phi_coefficients)
        values[index] += float(delta)
        if abs(values[index]) > basis.coefficient_limit:
            raise ValueError("Phi(1,2) perturbation exceeds the registered coefficient bound")
        updated = replace(basis, phi_coefficients=tuple(values))
    elif channel == "swirl":
        values = list(basis.swirl_coefficients)
        values[index] += float(delta)
        if abs(values[index]) > basis.coefficient_limit:
            raise ValueError("F(1,2) perturbation exceeds the registered coefficient bound")
        updated = replace(basis, swirl_coefficients=tuple(values))
    else:
        raise ValueError("channel must be 'phi' or 'swirl'")
    return replace(candidate, profile_basis=updated)


def _response_rank_condition(response: np.ndarray) -> tuple[list[float], int, float]:
    singular_values = np.linalg.svd(response, compute_uv=False)
    if singular_values.size == 0 or singular_values[0] <= 0.0:
        raise RuntimeError("morphology response is numerically inactive")
    tolerance = 1e-10 * singular_values[0]
    rank = int(np.count_nonzero(singular_values > tolerance))
    condition = (
        float(singular_values[0] / singular_values[-1])
        if rank == response.shape[1]
        else float("inf")
    )
    return singular_values.tolist(), rank, condition


def _novelty_fraction(column: np.ndarray, reference: np.ndarray) -> float:
    """Return the fraction of ``column`` orthogonal to the reference column span."""
    column = np.asarray(column, dtype=float)
    reference = np.asarray(reference, dtype=float)
    denominator = float(np.linalg.norm(column))
    if denominator <= 0.0:
        raise RuntimeError("mixed morphology response column is numerically inactive")
    q, _ = np.linalg.qr(reference, mode="reduced")
    residual = column - q @ (q.T @ column)
    return float(np.linalg.norm(residual) / denominator)


def audit_eq45_optimized_mixed_taper_capacity(
    candidate: Eq45VelocityCandidate | None = None,
    *,
    times: Iterable[float] = REFERENCE_TIMES,
    coefficient_steps: Iterable[float] = (0.04, 0.02),
    n_angles: int = 16,
) -> dict[str, Any]:
    """Audit existing ``(1,2)`` controls against the existing ``(1,0)`` span."""
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
    try:
        mode_index = basis.mode_indices.index(MODE)
    except ValueError as exc:
        raise ValueError("optimized candidate lacks the existing (1,2) mode") from exc
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
            plus = _with_existing_mixed_delta(current, channel=channel, delta=step)
            minus = _with_existing_mixed_delta(current, channel=channel, delta=-step)
            response = (
                _feature_vector(plus, times_tuple, n_angles=int(n_angles))
                - _feature_vector(minus, times_tuple, n_angles=int(n_angles))
            ) / (2.0 * step)
            columns.append(response)
        matrices.append(np.column_stack(columns))

    mixed_response = matrices[-1]
    mixed_singular_values, mixed_rank, mixed_condition = _response_rank_condition(mixed_response)
    response_norm = float(np.linalg.norm(mixed_response))
    refinement_change = float(np.linalg.norm(matrices[0] - mixed_response) / response_norm)

    radial_report = audit_eq45_optimized_even_radial_capacity(
        current,
        times=times_tuple,
        coefficient_steps=steps,
        n_angles=int(n_angles),
    )
    radial_response = np.asarray(
        radial_report["finest_morphology_response_per_unit_coefficient"], dtype=float
    )
    combined_response = np.column_stack((radial_response, mixed_response))
    combined_singular_values, combined_rank, combined_condition = _response_rank_condition(
        combined_response
    )
    novelty = [
        _novelty_fraction(mixed_response[:, column], radial_response)
        for column in range(mixed_response.shape[1])
    ]

    feature_count = len(FEATURE_NAMES)
    axial_reach_indices = np.arange(1, mixed_response.shape[0], feature_count)
    radial_thickness_indices = np.arange(0, mixed_response.shape[0], feature_count)
    channel_summary = []
    for column, label in enumerate(PARAMETER_LABELS):
        total = float(np.linalg.norm(mixed_response[:, column]))
        channel_summary.append(
            {
                "parameter": label,
                "morphology_response_l2_per_unit": total,
                "novelty_fraction_outside_existing_10_span": novelty[column],
                "axial_reach_response_fraction": float(
                    np.linalg.norm(mixed_response[axial_reach_indices, column]) / total
                ),
                "radial_thickness_response_fraction": float(
                    np.linalg.norm(mixed_response[radial_thickness_indices, column]) / total
                ),
            }
        )

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
    for column, channel in enumerate(("phi", "swirl")):
        plus = _with_existing_mixed_delta(current, channel=channel, delta=finest_delta)
        minus = _with_existing_mixed_delta(current, channel=channel, delta=-finest_delta)
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
                "parameter": PARAMETER_LABELS[column],
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
        "claim_scope": "optimized_existing_mixed_taper_expression_capacity_only",
        "candidate_sha256": current.sha256,
        "times": list(times_tuple),
        "n_angles": int(n_angles),
        "coefficient_limit": float(basis.coefficient_limit),
        "coefficient_steps": list(steps),
        "parameter_labels": list(PARAMETER_LABELS),
        "baseline_parameter_values": baseline_coefficients.tolist(),
        "feature_names_per_time": list(FEATURE_NAMES),
        "baseline_feature_vector": baseline_features.tolist(),
        "finest_morphology_response_per_unit_coefficient": mixed_response.tolist(),
        "mixed_singular_values": mixed_singular_values,
        "mixed_numerical_rank": mixed_rank,
        "mixed_condition_number": mixed_condition,
        "step_refinement_relative_change": refinement_change,
        "existing_radial_span_rank": int(radial_report["numerical_rank"]),
        "combined_10_12_singular_values": combined_singular_values,
        "combined_10_12_rank": combined_rank,
        "combined_10_12_condition_number": combined_condition,
        "channel_summary": channel_summary,
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
            "Treat Phi(1,2)/F(1,2) as already-available mixed radial×axial controls. "
            "Use them before adding another taper/radial-axial basis channel only when the public-observable defect aligns with their measured response; prefer channels with substantial novelty outside the existing (1,0) span and independently acceptable validation cost."
        ),
        "truth_boundary": dict(TRUTH_BOUNDARY),
    }
