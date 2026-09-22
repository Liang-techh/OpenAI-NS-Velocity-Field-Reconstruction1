"""Quantify how one already-screened sixth direction reduces morphology obstruction.

Preregistered in issue #1155 and stacked exactly on Agent-7 PR #1147 head.
This is a target-free routing audit over the already-frozen #1147 probes. It
adds no candidate basis, selects no coefficient, and changes no velocity.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np

import agent7_st052m_representation_archetype_left_null_preflight as parent

TASK_ID = "CR003-ST052M-ONE-DIRECTION-OBSTRUCTION-REDUCTION-137"
PREREG_ISSUE = 1155
STACK_BASE_PR = 1147
STACK_BASE_HEAD = "f20a883b0b2721d940e9012651b737d52583349b"

CHANNELS = tuple(parent.CHANNELS)
DEFECT_COORDINATES = tuple(parent.DEFECT_COORDINATES)
PROBE_NAMES = tuple(parent.PROBE_NAMES)
RANK_TOL = float(parent.RANK_TOL)
RESIDUAL_FLOOR = 1.0e-12

TRUTH = {
    "candidate_velocity_changed": False,
    "canonical_velocity_changed": False,
    "saved_velocity_changed": False,
    "basis_dimension_changed": False,
    "new_probe_formula_added": False,
    "new_spatial_basis_added_to_candidate": False,
    "new_temporal_basis_added_to_candidate": False,
    "coefficient_selected": False,
    "parameter_scan_performed": False,
    "optimization_performed": False,
    "post_result_retuning_performed": False,
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


def _assert_source_lock() -> None:
    parent._assert_source_lock()
    if parent.TASK_ID != "CR003-ST052M-REPRESENTATION-ARCHETYPE-LEFT-NULL-PREFLIGHT-136":
        raise RuntimeError("#1147 parent task identity drifted")
    if parent.PREREG_ISSUE != 1146:
        raise RuntimeError("#1147 preregistration identity drifted")
    if parent.STACK_BASE_PR != 1139:
        raise RuntimeError("#1147 stack lineage drifted")
    if tuple(parent.PROBE_NAMES) != PROBE_NAMES:
        raise RuntimeError("#1147 probe order drifted")
    if tuple(parent.DEFECT_COORDINATES) != DEFECT_COORDINATES:
        raise RuntimeError("defect-coordinate order drifted")
    if tuple(parent.CHANNELS) != CHANNELS:
        raise RuntimeError("existing-control order drifted")


def _projector(matrix: np.ndarray) -> tuple[np.ndarray, np.ndarray, int]:
    m = np.asarray(matrix, dtype=float)
    if m.ndim != 2 or m.shape[0] != len(DEFECT_COORDINATES):
        raise ValueError("unexpected defect-space tangent shape")
    if not np.all(np.isfinite(m)):
        raise ValueError("nonfinite tangent entry")
    rank = int(np.linalg.matrix_rank(m, tol=RANK_TOL))
    p = m @ np.linalg.pinv(m, rcond=RANK_TOL)
    p = 0.5 * (p + p.T)
    q = np.eye(m.shape[0]) - p
    q = 0.5 * (q + q.T)
    return p, q, rank


def obstruction_reduction_for_probe(
    existing_balanced: np.ndarray,
    probe_balanced: np.ndarray,
    parent_probe_result: dict[str, Any],
) -> dict[str, Any]:
    """Pure projector audit used by the real report and synthetic regressions."""
    j5 = np.asarray(existing_balanced, dtype=float)
    d = np.asarray(probe_balanced, dtype=float)
    if j5.shape != (len(DEFECT_COORDINATES), len(CHANNELS)):
        raise ValueError("existing tangent must be 8x5")
    if d.shape != (len(DEFECT_COORDINATES),):
        raise ValueError("probe tangent must have one value per defect coordinate")
    if not np.all(np.isfinite(j5)) or not np.all(np.isfinite(d)):
        raise ValueError("nonfinite tangent input")
    if float(np.linalg.norm(d)) <= np.finfo(float).tiny:
        raise ValueError("zero probe tangent")

    p5, q5, rank5 = _projector(j5)
    j6 = np.column_stack((j5, d))
    p6, q6, rank6 = _projector(j6)
    pinv6 = np.linalg.pinv(j6, rcond=RANK_TOL)

    coordinate_rows: dict[str, Any] = {}
    for i, coordinate in enumerate(DEFECT_COORDINATES):
        target = np.zeros(len(DEFECT_COORDINATES), dtype=float)
        target[i] = 1.0
        projected5 = p5 @ target
        residual5 = q5 @ target
        projected6 = p6 @ target
        residual6 = q6 @ target
        before = float(np.linalg.norm(residual5))
        after = float(np.linalg.norm(residual6))
        reduction = float(before - after)
        fraction = None if before <= RESIDUAL_FLOOR else float(reduction / before)

        coefficients = pinv6 @ target
        off_target = projected6.copy()
        off_target[i] = 0.0
        abs_off = np.abs(off_target)
        j = int(np.argmax(abs_off))
        spill_name = DEFECT_COORDINATES[j] if float(abs_off[j]) > 0.0 else None
        spill_value = float(projected6[j]) if spill_name is not None else 0.0

        coordinate_rows[coordinate] = {
            "unavoidable_residual_l2_before": before,
            "unavoidable_residual_l2_after": after,
            "absolute_obstruction_reduction_l2": reduction,
            "fractional_obstruction_reduction": fraction,
            "leverage_before": float(target @ projected5),
            "leverage_after": float(target @ projected6),
            "leverage_gain": float(target @ (projected6 - projected5)),
            "minimum_norm_augmented_control_vector": {
                **{
                    channel: float(coefficients[k])
                    for k, channel in enumerate(CHANNELS)
                },
                "appended_probe": float(coefficients[-1]),
            },
            "minimum_norm_augmented_control_l2": float(np.linalg.norm(coefficients)),
            "appended_probe_coefficient": float(coefficients[-1]),
            "projected_balanced_defect_vector_after": {
                name: float(projected6[k])
                for k, name in enumerate(DEFECT_COORDINATES)
            },
            "off_target_projected_l2_after": float(np.linalg.norm(off_target)),
            "largest_off_target_projected_coordinate_after": spill_name,
            "largest_off_target_projected_value_after": spill_value,
            "largest_off_target_projected_abs_after": float(np.max(abs_off)),
        }

    reductions = {
        name: float(row["absolute_obstruction_reduction_l2"])
        for name, row in coordinate_rows.items()
    }
    best_coordinate = max(DEFECT_COORDINATES, key=lambda name: reductions[name])
    trace_reduction = float(np.trace(q5) - np.trace(q6))
    rank_gain = int(rank6 - rank5)
    parent_eligible = bool(parent_probe_result.get("independent_preflight_eligible", False))

    return {
        "parent_independent_preflight_eligible": parent_eligible,
        "eligible_for_descriptive_obstruction_ranking": parent_eligible,
        "existing_rank": rank5,
        "augmented_rank": rank6,
        "rank_gain": rank_gain,
        "left_nullity_before": int(len(DEFECT_COORDINATES) - rank5),
        "left_nullity_after": int(len(DEFECT_COORDINATES) - rank6),
        "left_null_projector_trace_before": float(np.trace(q5)),
        "left_null_projector_trace_after": float(np.trace(q6)),
        "left_null_trace_reduction": trace_reduction,
        "one_added_direction_maximum_possible_nullity_reduction": 1,
        "at_least_two_uncontrolled_combinations_remain_if_rank5_to6": bool(
            rank5 == 5 and rank6 == 6 and len(DEFECT_COORDINATES) - rank6 >= 2
        ),
        "coordinate_obstruction_reduction": coordinate_rows,
        "best_isolated_coordinate_by_absolute_obstruction_reduction": best_coordinate,
        "maximum_absolute_obstruction_reduction_l2": reductions[best_coordinate],
        "mean_absolute_obstruction_reduction_l2": float(np.mean(list(reductions.values()))),
        "projector_monotonicity_min_reduction": float(min(reductions.values())),
        "candidate_promotion_authorized": False,
        "coefficient_selection_authorized": False,
    }


def build_report() -> dict[str, Any]:
    _assert_source_lock()
    parent_report = parent.build_report()
    scientific = parent._scientific_parent()
    scientific_audit = scientific.controllability_audit()
    projection = parent.parent.projection_audit_from_parent(scientific_audit)
    existing_balanced = np.asarray(projection["row_balanced_jacobian"], dtype=float)
    row_norms = np.asarray(
        [projection["raw_row_l2"][name] for name in DEFECT_COORDINATES], dtype=float
    )
    if existing_balanced.shape != (len(DEFECT_COORDINATES), len(CHANNELS)):
        raise RuntimeError("unexpected #1139 balanced tangent shape")
    if np.any(row_norms <= parent.ROW_RESPONSE_FLOOR):
        raise RuntimeError("invalid frozen #1139 row scale")

    results: dict[str, Any] = {}
    for name in PROBE_NAMES:
        parent_probe = parent_report["probe_results"][name]
        fine_raw = np.asarray(
            [
                float(parent_probe["fine_raw_derivative_by_coordinate"][coordinate])
                for coordinate in DEFECT_COORDINATES
            ],
            dtype=float,
        )
        balanced_probe = fine_raw / row_norms
        results[name] = obstruction_reduction_for_probe(
            existing_balanced, balanced_probe, parent_probe
        )
        results[name]["parent_augmented_normalized_column_condition"] = parent_probe[
            "augmented_normalized_column_condition"
        ]
        results[name]["parent_left_null_novelty_fraction_l2"] = parent_probe[
            "left_null_novelty_fraction_l2"
        ]

    eligible = [
        name
        for name in PROBE_NAMES
        if results[name]["eligible_for_descriptive_obstruction_ranking"]
    ]
    ranked = sorted(
        eligible,
        key=lambda name: (
            results[name]["maximum_absolute_obstruction_reduction_l2"],
            -float(results[name]["parent_augmented_normalized_column_condition"]),
        ),
        reverse=True,
    )

    return {
        "task_id": TASK_ID,
        "prereg_issue": PREREG_ISSUE,
        "stack_base": {"pr": STACK_BASE_PR, "head": STACK_BASE_HEAD},
        "parent_task": parent.TASK_ID,
        "frozen_protocol": {
            "existing_tangent_source": "#1139 row-balanced 8x5 Jacobian",
            "probe_source": "#1147 three already-frozen fine morphology tangents only",
            "probe_names": list(PROBE_NAMES),
            "defect_coordinate_order": list(DEFECT_COORDINATES),
            "existing_control_order": list(CHANNELS),
            "residual_floor_for_fraction": RESIDUAL_FLOOR,
            "ranking_requires_parent_independent_preflight_eligible": True,
            "ranking_metric": "maximum absolute isolated-coordinate obstruction reduction",
            "ranking_tiebreak": "lower #1147 augmented normalized-column condition",
            "scientific_pass_fail_threshold_on_obstruction_reduction": None,
            "new_probe_formula_count": 0,
            "candidate_basis_dimension_added": 0,
            "public_openai_numeric_target": None,
        },
        "probe_obstruction_reduction": results,
        "eligible_probe_archetypes": eligible,
        "eligible_probe_archetypes_ranked_by_obstruction_reduction": ranked,
        "decision": {
            "basis_candidate_promoted": None,
            "candidate_mutation_authorized": False,
            "coefficient_selection_authorized": False,
            "interpretation": (
                "this audit measures which isolated morphology coordinates each already-screened "
                "one-direction extension could make more reachable. A real candidate-bound "
                "public-observable discrepancy must still identify the needed defect direction "
                "before any one-basis child is justified."
            ),
            "one_direction_cannot_close_full_eight_coordinate_gap": True,
            "actual_velocity_changed": False,
            "direct_visualization_fingerprint_improvement": 0.0,
            "closer_visualization_delivery_established": False,
        },
        "held_out_pde_residual": {
            "evaluated": False,
            "st006_comparison_performed": False,
            "reason": "projector/routing audit only; candidate velocity is unchanged",
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
    print(
        json.dumps(
            {
                "eligible": report["eligible_probe_archetypes"],
                "ranked": report[
                    "eligible_probe_archetypes_ranked_by_obstruction_reduction"
                ],
                "best_coordinate_by_probe": {
                    name: result[
                        "best_isolated_coordinate_by_absolute_obstruction_reduction"
                    ]
                    for name, result in report["probe_obstruction_reduction"].items()
                },
                "decision": report["decision"],
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
