"""Lift the frozen time-grid alias with one off-grid temporal morphology coordinate.

Preregistered in issue #1177 and stacked exactly on Agent-7 PR #1167 head.
This is an analysis-only capacity diagnostic. It changes no candidate velocity,
basis dimension, coefficient, pressure, forcing, renderer, or scientific state.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Callable

import numpy as np

import agent7_st052m_temporal_sampling_alias_audit as alias

TASK_ID = "CR003-ST052M-OFFGRID-TEMPORAL-MORPHOLOGY-139"
PREREG_ISSUE = 1177
STACK_BASE_PR = 1167
STACK_BASE_HEAD = "3fc6286c4b338126ce60eb167ed082a6af953335"

T_LEFT = 0.375
T_RIGHT = 0.625
SUPPORT_RADIUS_SCALE = 2.0
EPSILONS = (1.0e-3, 5.0e-4)
PRIMARY_EPSILON = 5.0e-4
PROBE_RESPONSE_FLOOR = 1.0e-10
DERIVATIVE_DRIFT_MAX = 0.05
RANK_TOL = 1.0e-10
COORDINATE_NAME = "core_width_skew_0375_to_0625"
PROBE_NAME = "time_skewed_poloidal_extension"

TRUTH = {
    "candidate_velocity_changed": False,
    "canonical_velocity_changed": False,
    "saved_velocity_changed": False,
    "basis_dimension_changed": False,
    "new_spatial_basis_added_to_candidate": False,
    "new_temporal_basis_added_to_candidate": False,
    "coefficient_selected": False,
    "optimization_performed": False,
    "post_result_retuning_performed": False,
    "pressure_or_force_changed": False,
    "free_residual_force_used": False,
    "amplitude_collapse_used": False,
    "held_out_pde_residual_evaluated": False,
    "public_image_numeric_target_used": False,
    "pixel_similarity_objective_used": False,
    "renderer_or_camera_fit_used": False,
    "direct_visualization_fingerprint_improvement": 0.0,
    "closer_visualization_delivery_established": False,
    "visualization_ready": False,
    "visual_correspondence_verified": False,
    "source_correspondence_verified": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


def _scientific_parent():
    return alias._scientific_parent()


def _observable_parent():
    return _scientific_parent().parent


def _assert_source_lock() -> None:
    alias._assert_source_lock()
    if alias.TASK_ID != "CR003-ST052M-TEMPORAL-SAMPLING-ALIAS-AUDIT-138":
        raise RuntimeError("#1167 alias parent identity drifted")
    if alias.PREREG_ISSUE != 1165 or alias.STACK_BASE_PR != 1156:
        raise RuntimeError("#1167 preregistration/lineage drifted")
    if tuple(alias.OBSERVATION_TIMES) != (0.25, 0.50, 0.75):
        raise RuntimeError("#1167 observation times drifted")
    if tuple(alias.OFFGRID_WITNESS_TIMES) != (T_LEFT, T_RIGHT):
        raise RuntimeError("#1167 off-grid witness times drifted")

    scientific = _scientific_parent()
    if scientific.TASK_ID != "CR003-ST052M-DEFECT-COORDINATE-CONTROLLABILITY-128":
        raise RuntimeError("#1078 scientific parent identity drifted")
    if tuple(scientific.EPSILONS) != EPSILONS:
        raise RuntimeError("#1078 finite-difference steps drifted")
    if scientific.PRIMARY_EPSILON != PRIMARY_EPSILON:
        raise RuntimeError("#1078 primary finite-difference step drifted")
    if scientific.SUPPORT_RADIUS != SUPPORT_RADIUS_SCALE:
        raise RuntimeError("#1078 support-radius normalization drifted")
    if PROBE_NAME not in alias.rep.PROBE_NAMES:
        raise RuntimeError("#1147 time-skewed probe disappeared")


def _core_width(field_fn: Callable, time: float) -> float:
    _profile, width = _observable_parent()._core_speed_profile(field_fn, float(time))
    value = float(width)
    if not np.isfinite(value) or value <= 0.0:
        raise RuntimeError("invalid off-grid core width")
    return value


def core_width_skew_coordinate(field_fn: Callable) -> float:
    """Frozen dimensionless off-grid core-width skew coordinate."""
    left = _core_width(field_fn, T_LEFT)
    right = _core_width(field_fn, T_RIGHT)
    value = (right - left) / SUPPORT_RADIUS_SCALE
    if not np.isfinite(value):
        raise RuntimeError("nonfinite off-grid core-width skew")
    return float(value)


def _perturbed_field(
    base_fn: Callable, tangent_fn: Callable, signed_epsilon: float
) -> Callable:
    def field(points: np.ndarray, time: float) -> np.ndarray:
        pts = np.asarray(points, dtype=float)
        return np.asarray(base_fn(pts, float(time)), dtype=float) + float(
            signed_epsilon
        ) * np.asarray(tangent_fn(pts, float(time)), dtype=float)

    return field


def coordinate_derivative(
    base_fn: Callable, tangent_fn: Callable, epsilon: float
) -> float:
    eps = float(epsilon)
    if not np.isfinite(eps) or eps <= 0.0:
        raise ValueError("epsilon must be finite and positive")
    plus = core_width_skew_coordinate(_perturbed_field(base_fn, tangent_fn, eps))
    minus = core_width_skew_coordinate(_perturbed_field(base_fn, tangent_fn, -eps))
    value = (plus - minus) / (2.0 * eps)
    if not np.isfinite(value):
        raise RuntimeError("nonfinite off-grid coordinate derivative")
    return float(value)


def _derivative_gate(coarse: float, fine: float) -> dict[str, Any]:
    coarse = float(coarse)
    fine = float(fine)
    if not np.isfinite(coarse) or not np.isfinite(fine):
        raise RuntimeError("nonfinite derivative gate input")
    response_nonzero = abs(fine) > PROBE_RESPONSE_FLOOR
    same_nonzero_sign = coarse * fine > 0.0
    relative_drift = (
        None
        if abs(fine) <= np.finfo(float).tiny
        else float(abs(fine - coarse) / abs(fine))
    )
    stable = bool(
        response_nonzero
        and same_nonzero_sign
        and relative_drift is not None
        and relative_drift <= DERIVATIVE_DRIFT_MAX
    )
    return {
        "coarse_derivative": coarse,
        "fine_derivative": fine,
        "absolute_fine_response": abs(fine),
        "response_floor": PROBE_RESPONSE_FLOOR,
        "response_nonzero": bool(response_nonzero),
        "coarse_fine_same_nonzero_sign": bool(same_nonzero_sign),
        "relative_derivative_drift": relative_drift,
        "relative_drift_gate": DERIVATIVE_DRIFT_MAX,
        "stable_nonzero_response": stable,
    }


def _column_normalized_rank(matrix: np.ndarray) -> tuple[int, float | None]:
    m = np.asarray(matrix, dtype=float)
    if m.ndim != 2 or not np.all(np.isfinite(m)):
        raise RuntimeError("invalid local sensitivity matrix")
    norms = np.linalg.norm(m, axis=0)
    normalized = np.zeros_like(m)
    good = norms > np.finfo(float).tiny
    normalized[:, good] = m[:, good] / norms[good]
    rank = int(np.linalg.matrix_rank(normalized, tol=RANK_TOL))
    singular = np.linalg.svd(normalized, compute_uv=False)
    condition = (
        None
        if singular.size == 0 or singular[-1] <= np.finfo(float).tiny
        else float(singular[0] / singular[-1])
    )
    return rank, condition


def build_report() -> dict[str, Any]:
    _assert_source_lock()
    scientific = _scientific_parent()
    observable = _observable_parent()
    base_fn, existing_tangents, metadata = observable._build_field_and_tangents()
    channels = tuple(observable._joint().CHANNELS)
    if channels != (
        "amplitude",
        "radial_shape",
        "axial_turnover",
        "temporal_curvature",
        "toroidal_swirl",
    ):
        raise RuntimeError("existing control order drifted")

    probe_tangent = alias.rep._probe_tangents()[PROBE_NAME]
    alias_receipt = alias.build_report()
    baseline_coordinate = core_width_skew_coordinate(base_fn)

    derivatives: dict[float, dict[str, float]] = {}
    for epsilon in EPSILONS:
        row = {
            channel: coordinate_derivative(
                base_fn, existing_tangents[channel], float(epsilon)
            )
            for channel in channels
        }
        row[PROBE_NAME] = coordinate_derivative(base_fn, probe_tangent, float(epsilon))
        derivatives[float(epsilon)] = row

    coarse = derivatives[EPSILONS[0]]
    fine = derivatives[PRIMARY_EPSILON]
    probe_gate = _derivative_gate(coarse[PROBE_NAME], fine[PROBE_NAME])
    existing_row = np.asarray([fine[channel] for channel in channels], dtype=float)
    existing_row_l2 = float(np.linalg.norm(existing_row))
    response_ratio = (
        None
        if existing_row_l2 <= np.finfo(float).tiny
        else float(abs(fine[PROBE_NAME]) / existing_row_l2)
    )

    baseline_speed = observable._baseline_speed_scale(base_fn)
    j8_existing = np.column_stack(
        [
            scientific._coordinate_derivative(
                base_fn,
                existing_tangents[channel],
                PRIMARY_EPSILON,
                baseline_speed,
            )
            for channel in channels
        ]
    )
    d8_probe = scientific._coordinate_derivative(
        base_fn, probe_tangent, PRIMARY_EPSILON, baseline_speed
    )
    if j8_existing.shape != (8, 5) or d8_probe.shape != (8,):
        raise RuntimeError("unexpected frozen defect-space shape")

    j9_existing = np.vstack((j8_existing, existing_row.reshape(1, -1)))
    d9_probe = np.concatenate((d8_probe, np.asarray([fine[PROBE_NAME]], dtype=float)))
    j9_augmented = np.column_stack((j9_existing, d9_probe))
    existing_rank, existing_condition = _column_normalized_rank(j9_existing)
    augmented_rank, augmented_condition = _column_normalized_rank(j9_augmented)

    alias_lifted = bool(
        alias_receipt["temporal_sampling_alias_detected"]
        and probe_gate["stable_nonzero_response"]
    )
    rank_gain = int(augmented_rank - existing_rank)

    return {
        "task_id": TASK_ID,
        "preregistered_issue": PREREG_ISSUE,
        "stack_base_pr": STACK_BASE_PR,
        "stack_base_head": STACK_BASE_HEAD,
        "frozen_protocol": {
            "coordinate_name": COORDINATE_NAME,
            "time_left": T_LEFT,
            "time_right": T_RIGHT,
            "support_radius_scale": SUPPORT_RADIUS_SCALE,
            "definition": "(speed_weighted_core_width(t=.625)-speed_weighted_core_width(t=.375))/2.0",
            "positive_direction": "wider core at the later symmetric off-grid time",
            "epsilons": list(EPSILONS),
            "primary_epsilon": PRIMARY_EPSILON,
            "probe_name": PROBE_NAME,
            "probe_response_floor": PROBE_RESPONSE_FLOOR,
            "derivative_drift_max": DERIVATIVE_DRIFT_MAX,
            "public_openai_numeric_target": None,
        },
        "base_metadata": metadata,
        "baseline_coordinate": float(baseline_coordinate),
        "coarse_derivative_by_direction": {
            key: float(value) for key, value in coarse.items()
        },
        "fine_derivative_by_direction": {
            key: float(value) for key, value in fine.items()
        },
        "existing_five_new_coordinate_row_l2": existing_row_l2,
        "time_skewed_probe_to_existing_row_response_ratio": response_ratio,
        "time_skewed_probe_gate": probe_gate,
        "frozen_eight_coordinate_probe_fine_l2": float(np.linalg.norm(d8_probe)),
        "nine_coordinate_existing_rank": existing_rank,
        "nine_coordinate_existing_condition": existing_condition,
        "nine_coordinate_augmented_rank": augmented_rank,
        "nine_coordinate_augmented_condition": augmented_condition,
        "nine_coordinate_rank_gain_from_time_skewed_probe": rank_gain,
        "parent_alias_detected_by_source_execution": bool(
            alias_receipt["temporal_sampling_alias_detected"]
        ),
        "temporal_sampling_alias_lifted_for_this_coordinate": alias_lifted,
        "decision": {
            "time_skewed_probe_may_reenter_capacity_comparison_for_this_offgrid_coordinate": alias_lifted,
            "new_basis_addition_authorized": False,
            "coefficient_selection_authorized": False,
            "candidate_mutation_authorized": False,
            "actual_velocity_changed": False,
            "direct_visualization_fingerprint_improvement": 0.0,
            "next_action": (
                "measure a real candidate-bound temporal morphology discrepancy on the same off-grid coordinate before any one-basis child"
                if alias_lifted
                else "retain the temporal probe as unresolved or negative evidence; do not retune the frozen diagnostic"
            ),
        },
        "held_out_pde_residual": {
            "evaluated": False,
            "reason": "no candidate velocity or coefficient changed in this diagnostic increment",
        },
        "truth": dict(TRUTH),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = build_report()
    text = json.dumps(report, indent=2, sort_keys=True)
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text + "\n", encoding="utf-8")
    print(text)


if __name__ == "__main__":
    main()
