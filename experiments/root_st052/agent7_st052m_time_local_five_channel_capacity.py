"""Time-local capacity audit for the existing five ST052-M morphology controls.

Preregistered in issue #1048 and stacked exactly on Agent-7 PR #1042 head
63a114c466041f8e3281e0f80e8268611505828d.

#1033 established one global five-channel response cloud and #1042 localized
that capacity in space.  This increment asks the orthogonal representation
question: at which registered times are the already-existing controls actually
available and jointly identifiable?

No basis or candidate is changed.  In particular, the audit makes the endpoint
semantics explicit: the current correction family preserves t=.25, while the
endpoint-preserving temporal-curvature tangent vanishes at both t=.25 and
 t=.75.  Those structural zeros are useful routing information for a future
candidate-bound morphology discrepancy; they are not visual correspondence.
"""
from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path
from typing import Any

import numpy as np

import agent7_st052m_regional_five_channel_capacity_localization as regional

TASK_ID = "CR003-ST052M-TIME-LOCAL-FIVE-CHANNEL-CAPACITY-125"
PREREG_ISSUE = 1048
SOURCE_PARENT_PR = 1042
SOURCE_PARENT_HEAD = "63a114c466041f8e3281e0f80e8268611505828d"

TIME_SLICES = (0.25, 0.375, 0.50, 0.625, 0.75)
ACTIVE_RELATIVE_NORM_MIN = 1.0e-12
CONDITION_MAX = regional.joint.CONDITION_MAX
COSINE_MAX_ABS = regional.joint.COSINE_MAX_ABS
LEAKAGE_MAX = regional.LEAKAGE_MAX

EXPECTED_INACTIVE = {
    "0.250": tuple(regional.joint.CHANNELS),
    "0.375": (),
    "0.500": (),
    "0.625": (),
    "0.750": ("temporal_curvature",),
}

TRUTH = {
    "candidate_velocity_changed": False,
    "canonical_velocity_changed": False,
    "saved_velocity_changed": False,
    "basis_dimension_changed": False,
    "new_spatial_basis_added": False,
    "new_temporal_basis_added": False,
    "coefficient_selected": False,
    "parameter_scan_performed": False,
    "optimization_performed": False,
    "post_result_retuning_performed": False,
    "pressure_or_force_changed": False,
    "held_out_pde_residual_evaluated": False,
    "st006_comparison_performed": False,
    "public_image_numeric_target_used": False,
    "pixel_similarity_objective_used": False,
    "time_local_visualization_sensitivity_audited": True,
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


def _time_key(time: float) -> str:
    return f"{float(time):.3f}"


def _assert_source_lock() -> None:
    regional._assert_source_lock()
    if regional.TASK_ID != "CR003-ST052M-REGIONAL-FIVE-CHANNEL-CAPACITY-124":
        raise RuntimeError("#1042 task identity drifted")
    if regional.PREREG_ISSUE != 1041:
        raise RuntimeError("#1042 preregistration drifted")
    if regional.SOURCE_PARENT_PR != 1033:
        raise RuntimeError("#1042 parent identity drifted")
    if tuple(regional.joint.CHANNELS) != (
        "amplitude",
        "radial_shape",
        "axial_turnover",
        "temporal_curvature",
        "toroidal_swirl",
    ):
        raise RuntimeError("five-channel set drifted")
    if tuple(regional.joint.axial.RESPONSE_TIMES) != (0.375, 0.50, 0.625, 0.75):
        raise RuntimeError("registered response times drifted")
    if regional.joint.temporal.TIME_START != 0.25 or regional.joint.temporal.TIME_END != 0.75:
        raise RuntimeError("registered temporal endpoints drifted")
    if CONDITION_MAX != 25.0 or COSINE_MAX_ABS != 0.995:
        raise RuntimeError("parent conditioning gates drifted")


def _build_blocks_at_time(
    time: float,
    points: np.ndarray,
    parent_fn: Any,
    s_alphas: dict[float, float],
    s_children: dict[float, Any],
    k_children: dict[float, Any],
    toroidal_spatial: np.ndarray,
) -> dict[str, np.ndarray]:
    """Evaluate the exact five #1033 tangents on one frozen time slice."""
    t = float(time)
    axial = regional.joint.axial
    radial = axial.radial_shape
    temporal = regional.joint.temporal
    toroidal = regional.joint.toroidal

    parent = np.asarray(parent_fn(points, t), dtype=float)
    live = np.asarray(s_children[axial.S_LIVE](points, t), dtype=float)
    k_minus = np.asarray(k_children[radial.K_MINUS](points, t), dtype=float)
    k_plus = np.asarray(k_children[radial.K_PLUS](points, t), dtype=float)
    s_minus = np.asarray(s_children[axial.S_MINUS](points, t), dtype=float)
    s_plus = np.asarray(s_children[axial.S_PLUS](points, t), dtype=float)
    live_alpha = float(s_alphas[axial.S_LIVE])

    blocks = {
        "amplitude": live - parent,
        "radial_shape": (k_plus - k_minus) / (radial.K_PLUS - radial.K_MINUS),
        "axial_turnover": (s_plus - s_minus) / (axial.S_PLUS - axial.S_MINUS),
        "temporal_curvature": temporal.temporal_curvature_tangent(points, t, live_alpha),
        "toroidal_swirl": toroidal.activation_g1(t) * toroidal_spatial,
    }
    for name, block in blocks.items():
        arr = np.asarray(block, dtype=float)
        if arr.shape != points.shape:
            raise RuntimeError(f"time-local response-shape drift: {name}")
        if not np.all(np.isfinite(arr)):
            raise RuntimeError(f"nonfinite time-local response: {name}")
    return blocks


def _component_fingerprint(points: np.ndarray, block: np.ndarray) -> dict[str, Any]:
    cyl = regional.joint.toroidal.cylindrical_components(points, np.asarray(block, dtype=float))
    rms = np.sqrt(np.mean(cyl * cyl, axis=0))
    total = float(np.sqrt(np.mean(np.sum(cyl * cyl, axis=1))))
    denom = max(total * total, np.finfo(float).tiny)
    labels = ("u_r", "u_theta", "u_z")
    return {
        "cylindrical_component_rms": {label: float(value) for label, value in zip(labels, rms)},
        "normalized_component_energy_fractions": {
            label: float(value * value / denom) for label, value in zip(labels, rms)
        },
        "total_vector_rms": total,
    }


def _audit_time_slice(
    time: float,
    points: np.ndarray,
    blocks: dict[str, np.ndarray],
    full_norms: dict[str, float],
) -> dict[str, Any]:
    channels = tuple(regional.joint.CHANNELS)
    raw_norms = {
        name: float(np.linalg.norm(np.asarray(blocks[name], dtype=float).ravel())) for name in channels
    }
    relative_norms = {
        name: float(raw_norms[name] / full_norms[name]) for name in channels
    }
    active = tuple(name for name in channels if relative_norms[name] > ACTIVE_RELATIVE_NORM_MIN)
    inactive = tuple(name for name in channels if name not in active)
    expected_inactive = EXPECTED_INACTIVE[_time_key(time)]
    inactive_matches_expected = bool(inactive == expected_inactive)

    if active:
        columns = []
        for name in active:
            vector = np.asarray(blocks[name], dtype=float).ravel()
            columns.append(vector / raw_norms[name])
        matrix = np.column_stack(columns)
        singular = np.linalg.svd(matrix, compute_uv=False)
        rank = int(np.linalg.matrix_rank(matrix, tol=1.0e-10))
        condition = float(singular[0] / singular[-1]) if singular[-1] > 0.0 else None
        gram = matrix.T @ matrix
        cosines: dict[str, float] = {}
        for i, j in itertools.combinations(range(len(active)), 2):
            cosines[f"{active[i]}__{active[j]}"] = float(gram[i, j])
        max_abs_cosine = float(max(abs(v) for v in cosines.values())) if cosines else 0.0
        active_subspace_pass = bool(
            rank == len(active)
            and condition is not None
            and condition <= CONDITION_MAX
            and max_abs_cosine < COSINE_MAX_ABS
        )
        singular_values = [float(v) for v in singular]
    else:
        rank = 0
        condition = None
        cosines = {}
        max_abs_cosine = None
        singular_values = []
        active_subspace_pass = bool(_time_key(time) == "0.250")

    fingerprints = {
        name: _component_fingerprint(points, blocks[name]) for name in channels
    }
    tor = fingerprints["toroidal_swirl"]["cylindrical_component_rms"]
    toroidal_rz_leakage = float(max(abs(tor["u_r"]), abs(tor["u_z"])))
    poloidal_theta_leakage = float(
        max(
            abs(fingerprints[name]["cylindrical_component_rms"]["u_theta"])
            for name in channels
            if name != "toroidal_swirl"
        )
    )
    selectivity_pass = bool(
        toroidal_rz_leakage <= LEAKAGE_MAX and poloidal_theta_leakage <= LEAKAGE_MAX
    )

    return {
        "time": float(time),
        "raw_column_norms": raw_norms,
        "relative_norm_of_full_spacetime_channel": relative_norms,
        "active_relative_norm_gate": ACTIVE_RELATIVE_NORM_MIN,
        "active_channels": list(active),
        "inactive_channels": list(inactive),
        "expected_inactive_channels": list(expected_inactive),
        "inactive_set_matches_preregistered_expectation": inactive_matches_expected,
        "active_channel_count": int(len(active)),
        "active_subspace_rank": rank,
        "active_subspace_condition": condition,
        "condition_gate": CONDITION_MAX,
        "pairwise_cosines_active_subspace": cosines,
        "max_abs_pairwise_cosine_active_subspace": max_abs_cosine,
        "cosine_abs_gate": COSINE_MAX_ABS,
        "singular_values_active_subspace": singular_values,
        "active_subspace_pass": active_subspace_pass,
        "component_fingerprints": fingerprints,
        "toroidal_rz_leakage_rms_max": toroidal_rz_leakage,
        "poloidal_theta_leakage_rms_max": poloidal_theta_leakage,
        "selectivity_guard_pass": selectivity_pass,
        "passes_preregistered_time_slice_gate": bool(
            inactive_matches_expected and active_subspace_pass and selectivity_pass
        ),
    }


def time_local_capacity() -> dict[str, Any]:
    """Audit active morphology capacity independently at each frozen time slice."""
    _assert_source_lock()
    joint = regional.joint
    points, full_blocks, metadata = joint._build_tangent_blocks()
    full_norms = {
        name: float(
            np.linalg.norm(
                np.concatenate([np.asarray(v, dtype=float).ravel() for v in full_blocks[name]])
            )
        )
        for name in joint.CHANNELS
    }
    if any((not np.isfinite(v)) or v <= np.finfo(float).tiny for v in full_norms.values()):
        raise RuntimeError("zero/nonfinite full-spacetime channel norm")

    parent_fn, s_alphas, _, s_children, _, k_children = joint.axial.build_local_families()
    toroidal_spatial = joint.toroidal.toroidal_unit_velocity(points)
    by_time: dict[str, Any] = {}
    for time in TIME_SLICES:
        blocks = _build_blocks_at_time(
            float(time), points, parent_fn, s_alphas, s_children, k_children, toroidal_spatial
        )
        by_time[_time_key(time)] = _audit_time_slice(
            float(time), points, blocks, full_norms
        )

    failing = [
        key for key, value in by_time.items() if not value["passes_preregistered_time_slice_gate"]
    ]
    unexpected_inactive = [
        key
        for key, value in by_time.items()
        if not value["inactive_set_matches_preregistered_expectation"]
    ]
    return {
        **metadata,
        "time_order": [float(v) for v in TIME_SLICES],
        "full_spacetime_channel_norms": full_norms,
        "active_relative_norm_gate": ACTIVE_RELATIVE_NORM_MIN,
        "expected_inactive_by_time": {
            key: list(value) for key, value in EXPECTED_INACTIVE.items()
        },
        "by_time": by_time,
        "all_time_slices_match_preregistered_capacity_semantics": bool(not failing),
        "failing_time_slices": failing,
        "unexpected_inactive_time_slices": unexpected_inactive,
        "start_endpoint_all_current_corrections_inactive": bool(
            not by_time["0.250"]["active_channels"]
        ),
        "end_endpoint_temporal_curvature_inactive": bool(
            "temporal_curvature" in by_time["0.750"]["inactive_channels"]
        ),
    }


def build_report() -> dict[str, Any]:
    capacity = time_local_capacity()
    start_locked = capacity["start_endpoint_all_current_corrections_inactive"]
    end_temporal_locked = capacity["end_endpoint_temporal_curvature_inactive"]
    return {
        "task_id": TASK_ID,
        "prereg_issue": PREREG_ISSUE,
        "source_parent": {"pr": SOURCE_PARENT_PR, "head": SOURCE_PARENT_HEAD},
        "frozen_protocol": {
            "channels": list(regional.joint.CHANNELS),
            "times": [float(v) for v in TIME_SLICES],
            "active_relative_norm_gate": ACTIVE_RELATIVE_NORM_MIN,
            "condition_gate": CONDITION_MAX,
            "cosine_abs_gate": COSINE_MAX_ABS,
            "expected_inactive_by_time": {
                key: list(value) for key, value in EXPECTED_INACTIVE.items()
            },
            "new_basis_dimension": 0,
            "coefficient_selected": None,
            "public_openai_numeric_target": None,
        },
        "time_local_capacity": capacity,
        "decision": {
            "start_endpoint_repair_available_in_current_correction_family": bool(not start_locked),
            "endpoint_preserving_temporal_channel_can_change_t075": bool(not end_temporal_locked),
            "additional_basis_dimension_justified_by_this_audit": False,
            "if_future_discrepancy_at_t025": (
                "inspect_parent_family_or_endpoint_nonpreserving_time_coordinate_before_basis_growth"
            ),
            "if_future_pure_time_discrepancy_at_t075": (
                "inspect_endpoint_nonzero_temporal_coordinate_or_parent_time_law_before_basis_growth"
            ),
            "if_future_interior_discrepancy": (
                "choose_smallest_existing_active_control_using_candidate_bound_discrepancy"
            ),
            "actual_velocity_changed": False,
            "direct_visualization_fingerprint_improvement": 0.0,
            "closer_visualization_delivery_established": False,
        },
        "held_out_pde_residual": {
            "evaluated": False,
            "st006_comparison_performed": False,
            "st006_momentum_max_reference": 0.1082289305112118,
            "st006_momentum_volume_l2_reference": 0.10758432876230622,
        },
        **TRUTH,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    report = build_report()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    cap = report["time_local_capacity"]
    print("all_time_slices_match=", cap["all_time_slices_match_preregistered_capacity_semantics"])
    print("failing_time_slices=", cap["failing_time_slices"])
    for key in (_time_key(v) for v in TIME_SLICES):
        row = cap["by_time"][key]
        print(
            key,
            "active=", row["active_channels"],
            "rank=", row["active_subspace_rank"],
            "condition=", row["active_subspace_condition"],
            "max_abs_cosine=", row["max_abs_pairwise_cosine_active_subspace"],
        )


if __name__ == "__main__":
    main()
