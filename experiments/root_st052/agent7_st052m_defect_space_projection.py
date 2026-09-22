"""Diagnose defect-space directions outside the existing five-control ST052-M tangent.

Preregistered in issue #1138 and stacked exactly on Agent-7 PR #1131.
This target-free basis-capacity audit reuses #1078's frozen 8x5 morphology
Jacobian.  It changes no velocity, coefficient, basis, renderer, pressure,
forcing, or PDE-validation state.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np

import agent7_st052m_minimum_existing_control_subset as stack_parent

TASK_ID = "CR003-ST052M-DEFECT-SPACE-PROJECTION-135"
PREREG_ISSUE = 1138
STACK_BASE_PR = 1131
STACK_BASE_HEAD = "ce93cda444a30add232251a895dd3953328fd9a1"
SCIENTIFIC_PARENT_PR = 1078
SCIENTIFIC_PARENT_TASK = "CR003-ST052M-DEFECT-COORDINATE-CONTROLLABILITY-128"

CHANNELS = (
    "amplitude",
    "radial_shape",
    "axial_turnover",
    "temporal_curvature",
    "toroidal_swirl",
)
DEFECT_COORDINATES = (
    "inward_radial_index_t050",
    "swirl_to_poloidal_ratio_t050",
    "axial_vorticity_aspect_t050",
    "core_width_t050",
    "tip_radial_thickness_t050",
    "angular_rotation_radial_variation_t050",
    "core_width_change_025_to_075",
    "core_speed_change_025_to_075",
)
ROW_RESPONSE_FLOOR = 1.0e-12
RANK_TOL = 1.0e-10

TRUTH = {
    "candidate_velocity_changed": False,
    "canonical_velocity_changed": False,
    "saved_velocity_changed": False,
    "basis_dimension_changed": False,
    "new_spatial_basis_added": False,
    "new_temporal_basis_added": False,
    "coefficient_selected": False,
    "optimization_performed": False,
    "pressure_or_force_changed": False,
    "held_out_pde_residual_evaluated": False,
    "st006_comparison_performed": False,
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
    return stack_parent.parent


def _assert_source_lock() -> None:
    stack_parent._assert_source_lock()
    if stack_parent.TASK_ID != "CR003-ST052M-MINIMUM-EXISTING-CONTROL-SUBSET-134":
        raise RuntimeError("#1131 stack parent task identity drifted")
    if stack_parent.PREREG_ISSUE != 1130:
        raise RuntimeError("#1131 preregistration identity drifted")
    scientific = _scientific_parent()
    if scientific.TASK_ID != SCIENTIFIC_PARENT_TASK:
        raise RuntimeError("#1078 scientific parent identity drifted")
    if scientific.PREREG_ISSUE != 1077:
        raise RuntimeError("#1078 preregistration identity drifted")
    if tuple(scientific.DEFECT_COORDINATES) != DEFECT_COORDINATES:
        raise RuntimeError("#1078 defect-coordinate order drifted")
    if tuple(scientific.parent._joint().CHANNELS) != CHANNELS:
        raise RuntimeError("#1078 control order drifted")
    if scientific.ROW_RESPONSE_FLOOR != ROW_RESPONSE_FLOOR:
        raise RuntimeError("#1078 row-response floor drifted")


def _raw_jacobian(audit: dict[str, Any]) -> np.ndarray:
    rows = []
    for coordinate in DEFECT_COORDINATES:
        try:
            mapping = audit["coordinates"][coordinate]["fine_raw_derivative_by_channel"]
            rows.append([float(mapping[channel]) for channel in CHANNELS])
        except KeyError as exc:
            raise RuntimeError(f"missing frozen #1078 Jacobian entry: {exc}") from exc
    matrix = np.asarray(rows, dtype=float)
    if matrix.shape != (len(DEFECT_COORDINATES), len(CHANNELS)):
        raise RuntimeError("unexpected frozen defect-Jacobian shape")
    if not np.all(np.isfinite(matrix)):
        raise RuntimeError("nonfinite frozen defect-Jacobian entry")
    return matrix


def _canonicalize_null_basis(basis: np.ndarray) -> np.ndarray:
    if basis.ndim != 2:
        raise ValueError("null basis must be a matrix")
    cols = []
    for j in range(basis.shape[1]):
        v = np.asarray(basis[:, j], dtype=float).copy()
        if v.size:
            pivot = int(np.argmax(np.abs(v)))
            if v[pivot] < 0.0:
                v *= -1.0
        cols.append(v)
    if not cols:
        return np.empty((basis.shape[0], 0), dtype=float)
    cols.sort(key=lambda v: tuple(np.round(v, 14).tolist()), reverse=True)
    return np.column_stack(cols)


def projection_audit_from_parent(audit: dict[str, Any]) -> dict[str, Any]:
    raw = _raw_jacobian(audit)
    row_norms = np.linalg.norm(raw, axis=1)
    if np.any(~np.isfinite(row_norms)):
        raise RuntimeError("nonfinite raw row norm")
    bad = [
        DEFECT_COORDINATES[i]
        for i, value in enumerate(row_norms)
        if float(value) <= ROW_RESPONSE_FLOOR
    ]
    if bad:
        raise RuntimeError(
            "cannot row-balance structurally unresponsive defect coordinate(s): "
            + ", ".join(bad)
        )

    balanced = raw / row_norms[:, None]
    u, singular, vh = np.linalg.svd(balanced, full_matrices=True)
    rank = int(np.linalg.matrix_rank(balanced, tol=RANK_TOL))
    pinv = np.linalg.pinv(balanced, rcond=RANK_TOL)
    projector = balanced @ pinv
    projector = 0.5 * (projector + projector.T)
    null_projector = np.eye(len(DEFECT_COORDINATES)) - projector
    null_projector = 0.5 * (null_projector + null_projector.T)

    projector_idempotence = float(np.linalg.norm(projector @ projector - projector, ord=2))
    null_idempotence = float(
        np.linalg.norm(null_projector @ null_projector - null_projector, ord=2)
    )
    orthogonality = float(np.linalg.norm(projector @ null_projector, ord=2))
    symmetry = float(np.linalg.norm(projector - projector.T, ord=2))

    null_basis = _canonicalize_null_basis(u[:, rank:])
    nullity = int(len(DEFECT_COORDINATES) - rank)

    coordinate_rows: dict[str, Any] = {}
    obstruction_diag = np.diag(null_projector)
    leverage_diag = np.diag(projector)
    for i, coordinate in enumerate(DEFECT_COORDINATES):
        target = np.zeros(len(DEFECT_COORDINATES), dtype=float)
        target[i] = 1.0
        coefficients = pinv @ target
        projected = projector @ target
        residual = target - projected

        off_target = np.abs(projected).copy()
        off_target[i] = 0.0
        j = int(np.argmax(off_target))
        spill_name = DEFECT_COORDINATES[j] if float(off_target[j]) > 0.0 else None
        spill_value = float(projected[j]) if spill_name is not None else 0.0

        coordinate_rows[coordinate] = {
            "raw_row_l2_before_balance": float(row_norms[i]),
            "leverage": float(leverage_diag[i]),
            "obstruction_exposure": float(obstruction_diag[i]),
            "best_attainable_projection_l2": float(np.linalg.norm(projected)),
            "unavoidable_orthogonal_residual_l2": float(np.linalg.norm(residual)),
            "minimum_norm_control_vector": {
                channel: float(coefficients[k]) for k, channel in enumerate(CHANNELS)
            },
            "minimum_norm_control_l2": float(np.linalg.norm(coefficients)),
            "minimum_norm_control_max_abs": float(np.max(np.abs(coefficients))),
            "projected_balanced_defect_vector": {
                name: float(projected[k])
                for k, name in enumerate(DEFECT_COORDINATES)
            },
            "orthogonal_residual_vector": {
                name: float(residual[k])
                for k, name in enumerate(DEFECT_COORDINATES)
            },
            "largest_off_target_projected_coordinate": spill_name,
            "largest_off_target_projected_value": spill_value,
            "largest_off_target_projected_abs": float(np.max(off_target)),
        }

    worst = max(
        DEFECT_COORDINATES,
        key=lambda name: coordinate_rows[name]["obstruction_exposure"],
    )
    best = min(
        DEFECT_COORDINATES,
        key=lambda name: coordinate_rows[name]["obstruction_exposure"],
    )

    return {
        "channel_order": list(CHANNELS),
        "defect_coordinate_order": list(DEFECT_COORDINATES),
        "raw_jacobian": raw.tolist(),
        "raw_row_l2": {
            name: float(row_norms[i]) for i, name in enumerate(DEFECT_COORDINATES)
        },
        "row_balanced_jacobian": balanced.tolist(),
        "row_balance_is_diagnostic_only": True,
        "rank": rank,
        "rank_tolerance": RANK_TOL,
        "singular_values": [float(x) for x in singular],
        "existing_control_tangent_dimension": rank,
        "balanced_defect_space_dimension": len(DEFECT_COORDINATES),
        "left_nullity": nullity,
        "expected_dimension_gap_if_rank_five": 3,
        "projector_trace": float(np.trace(projector)),
        "left_null_projector_trace": float(np.trace(null_projector)),
        "projector_idempotence_error_l2": projector_idempotence,
        "left_null_projector_idempotence_error_l2": null_idempotence,
        "projector_null_orthogonality_error_l2": orthogonality,
        "projector_symmetry_error_l2": symmetry,
        "existing_control_projector": projector.tolist(),
        "left_null_projector": null_projector.tolist(),
        "left_null_basis_svd_sign_canonicalized": null_basis.tolist(),
        "left_null_basis_orientation_authoritative": False,
        "projector_is_authoritative_under_null_basis_degeneracy": True,
        "coordinate_capacity": coordinate_rows,
        "most_obstruction_exposed_coordinate_for_isolated_change": worst,
        "least_obstruction_exposed_coordinate_for_isolated_change": best,
        "rank_below_five_existing_control_degeneracy": bool(rank < len(CHANNELS)),
        "independent_eight_coordinate_control_possible_with_five_local_directions": False,
    }


def build_report() -> dict[str, Any]:
    _assert_source_lock()
    parent_audit = _scientific_parent().controllability_audit()
    try:
        projection = projection_audit_from_parent(parent_audit)
        projection_status = "computed"
        projection_error = None
    except RuntimeError as exc:
        projection = None
        projection_status = "blocked_structural_or_numeric_input"
        projection_error = str(exc)
    return {
        "task_id": TASK_ID,
        "prereg_issue": PREREG_ISSUE,
        "stack_base": {"pr": STACK_BASE_PR, "head": STACK_BASE_HEAD},
        "scientific_parent": {
            "pr": SCIENTIFIC_PARENT_PR,
            "task": SCIENTIFIC_PARENT_TASK,
        },
        "frozen_protocol": {
            "source": "#1078 fine raw 8x5 defect-coordinate Jacobian only",
            "row_response_floor": ROW_RESPONSE_FLOOR,
            "rank_tolerance": RANK_TOL,
            "row_balance": "divide each raw defect-coordinate row by its own L2 norm",
            "scientific_pass_fail_threshold_on_projection": None,
            "public_openai_numeric_target": None,
            "new_basis_dimension": 0,
            "coefficient_selected": None,
        },
        "projection_status": projection_status,
        "projection_error": projection_error,
        "defect_space_projection_audit": projection,
        "decision": {
            "interpretation": (
                "this maps local defect combinations that the existing five-control tangent "
                "cannot independently realize; a future candidate-bound discrepancy must be "
                "projected onto this left-null space before one targeted basis is considered"
            ),
            "new_basis_dimension_justified_by_this_audit": False,
            "candidate_mutation_authorized_by_this_audit": False,
            "actual_velocity_changed": False,
            "direct_visualization_fingerprint_improvement": 0.0,
            "closer_visualization_delivery_established": False,
        },
        "held_out_pde_residual": {
            "evaluated": False,
            "st006_comparison_performed": False,
            "reason": "projection/capacity audit only; candidate velocity is unchanged",
        },
        **TRUTH,
    }


def run(out: Path) -> dict[str, Any]:
    report = build_report()
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return report


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    report = run(args.out)
    a = report["defect_space_projection_audit"]
    summary = {
        "projection_status": report["projection_status"],
        "projection_error": report["projection_error"],
        "decision": report["decision"],
    }
    if a is not None:
        summary.update({
            "rank": a["rank"],
            "left_nullity": a["left_nullity"],
            "most_obstruction_exposed_coordinate_for_isolated_change": a[
                "most_obstruction_exposed_coordinate_for_isolated_change"
            ],
            "least_obstruction_exposed_coordinate_for_isolated_change": a[
                "least_obstruction_exposed_coordinate_for_isolated_change"
            ],
            "rank_below_five_existing_control_degeneracy": a[
                "rank_below_five_existing_control_degeneracy"
            ],
        })
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
