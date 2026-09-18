"""Bound a future second covariance column before another Kokuno correction cycle.

PR #282 measured the real two-component compact-stress remainder not representable
by Agent 2's current covariance column. A nonparallel second column spans that
remainder algebraically, but a nearly parallel or tiny response can require an
arbitrarily large signed amplitude correction. The corrected Kokuno reader uses
an invertible two-column covariance map H and a signed differential inverse.

For current response c1, future response c2, and measured orthogonal remainder m,
only c2_perp = c2 - proj_c1(c2) removes m. At a rank-two node,

    |delta a_2| = |m| / |c2_perp|,
    |det[c1,c2]| = |c1| |c2_perp|.

Thus a declared signed-coefficient budget B implies the necessary local floors
|c2_perp| >= |m|/B and |det H| >= |c1||m|/B. These are algebraic preflight
bounds, not Navier--Stokes residual contraction claims.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np

from .kokuno_correction_cycle_gain_guard import (
    ST006_FINEST_SPATIAL_STEP,
    ST006_MOMENTUM_MAX,
    ST006_SHA256,
    ST006_VALIDATION_POINTS,
    ST006_VALIDATION_SEED,
    ST006_VOLUME_L2,
)
from .kokuno_missing_covariance_column_target import (
    MissingCovarianceColumnTarget,
    build_missing_covariance_column_target,
    generate_actual_core_report as generate_missing_target_report,
)
from .kokuno_signed_covariance_inverse import ROUTED_OSCILLATORY_AMPLITUDE

TASK = "KOKUNO-A3-SECOND-COLUMN-BOUNDED-INVERSE-010"
SOURCE_READER = "KokunoYumeto corrected 208-page reconstruction"
SOURCE_EDITION = "2026.09.09-consolidated"
SOURCE_DOI = "10.5281/zenodo.22678406"
SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
SOURCE_SECTION = "Signed covariance and exact signed differential inverse"
RANK_RELATIVE_TOLERANCE = 1.0e-10
ALGEBRAIC_RELATIVE_RESIDUAL_TOLERANCE = 1.0e-10
BUDGET_ROUNDOFF_RELATIVE_TOLERANCE = 1.0e-12


def _vector_rms(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    if values.ndim != 2 or values.shape[1] != 2:
        raise ValueError("values must have shape (N,2)")
    return float(np.sqrt(np.mean(np.sum(values * values, axis=1))))


def _second_needed_mask(receipt: MissingCovarianceColumnTarget) -> np.ndarray:
    target_norm = np.linalg.norm(receipt.target_stress, axis=1)
    missing_norm = np.linalg.norm(receipt.missing_response, axis=1)
    target_max = max(float(np.max(target_norm)), np.finfo(float).tiny)
    return (
        receipt.active_target_mask
        & receipt.safe_response_mask
        & (missing_norm > receipt.target_floor_fraction * target_max)
    )


def current_signed_coefficient_budget(receipt: MissingCovarianceColumnTarget) -> float:
    """Use the largest existing rank-one signed coefficient as a data-bound scale."""
    if not isinstance(receipt, MissingCovarianceColumnTarget):
        raise TypeError("receipt must be MissingCovarianceColumnTarget")
    mask = receipt.safe_active_mask
    if not np.any(mask):
        raise ValueError("receipt has no safe active nodes")
    budget = float(np.max(np.abs(receipt.current_coefficient[mask])))
    if not np.isfinite(budget) or budget <= np.finfo(float).tiny:
        raise ValueError("current signed-coefficient scale is numerically zero")
    return budget


def required_second_column_envelope(
    receipt: MissingCovarianceColumnTarget,
    *,
    coefficient_budget: float,
) -> dict[str, Any]:
    """Return the transverse-response and determinant floors implied by a budget."""
    if not isinstance(receipt, MissingCovarianceColumnTarget):
        raise TypeError("receipt must be MissingCovarianceColumnTarget")
    budget = float(coefficient_budget)
    if not np.isfinite(budget) or budget <= 0.0:
        raise ValueError("coefficient_budget must be positive and finite")
    needed = _second_needed_mask(receipt)
    if not np.any(needed):
        raise ValueError("receipt has no nodes requiring a second covariance direction")
    current_norm = np.linalg.norm(receipt.current_response[needed], axis=1)
    missing_norm = np.linalg.norm(receipt.missing_response[needed], axis=1)
    if np.any(current_norm <= np.finfo(float).tiny):
        raise ValueError("second-needed node has a numerically zero current response")
    required_transverse = missing_norm / budget
    required_determinant = current_norm * required_transverse
    relative_to_current = required_transverse / current_norm
    return {
        "coefficient_budget": budget,
        "nodes_requiring_second_direction": int(np.count_nonzero(needed)),
        "required_transverse_response_rms": float(np.sqrt(np.mean(required_transverse**2))),
        "required_transverse_response_max": float(np.max(required_transverse)),
        "required_transverse_over_current_median": float(np.median(relative_to_current)),
        "required_transverse_over_current_max": float(np.max(relative_to_current)),
        "required_determinant_abs_min_median": float(np.median(required_determinant)),
        "required_determinant_abs_min_max": float(np.max(required_determinant)),
        "radii_requiring_second_direction": receipt.radii[needed].tolist(),
        "required_transverse_response": required_transverse.tolist(),
        "required_determinant_abs_min": required_determinant.tolist(),
        "interpretation": (
            "Necessary local scale for keeping the future second signed coefficient within the declared budget; "
            "not a sufficient condition for Navier-Stokes residual reduction."
        ),
    }


def evaluate_bounded_second_column(
    receipt: MissingCovarianceColumnTarget,
    second_response: np.ndarray,
    *,
    coefficient_budget: float,
    relative_rank_tolerance: float = RANK_RELATIVE_TOLERANCE,
    algebraic_relative_residual_tolerance: float = ALGEBRAIC_RELATIVE_RESIDUAL_TOLERANCE,
) -> dict[str, Any]:
    """Measure rank, stress coverage, conditioning and signed-update magnitude."""
    if not isinstance(receipt, MissingCovarianceColumnTarget):
        raise TypeError("receipt must be MissingCovarianceColumnTarget")
    second = np.asarray(second_response, dtype=float)
    if second.shape != receipt.target_stress.shape or not np.isfinite(second).all():
        raise ValueError("second_response must be finite and match target shape")
    budget = float(coefficient_budget)
    if not np.isfinite(budget) or budget <= 0.0:
        raise ValueError("coefficient_budget must be positive and finite")
    if not 0.0 < relative_rank_tolerance < 1.0e-2:
        raise ValueError("relative_rank_tolerance must lie in (0,1e-2)")
    if not 0.0 < algebraic_relative_residual_tolerance < 1.0e-4:
        raise ValueError("algebraic_relative_residual_tolerance must lie in (0,1e-4)")

    active = receipt.active_target_mask
    needed = _second_needed_mask(receipt)
    rank = np.zeros(len(receipt.radii), dtype=int)
    condition = np.full(len(receipt.radii), np.nan, dtype=float)
    determinant = np.zeros(len(receipt.radii), dtype=float)
    alpha = np.zeros(len(receipt.radii), dtype=float)
    beta = np.zeros(len(receipt.radii), dtype=float)
    transverse = np.zeros(len(receipt.radii), dtype=float)
    residuals: list[np.ndarray] = []

    for index in np.flatnonzero(active):
        c1 = receipt.current_response[index]
        c2 = second[index]
        target = receipt.target_stress[index]
        matrix = np.column_stack((c1, c2))
        singular = np.linalg.svd(matrix, compute_uv=False)
        if singular[0] > np.finfo(float).tiny:
            rank[index] = int(
                np.count_nonzero(singular > relative_rank_tolerance * singular[0])
            )
        determinant[index] = float(abs(np.linalg.det(matrix)))
        c1_norm_sq = float(np.dot(c1, c1))
        if c1_norm_sq > np.finfo(float).tiny:
            c2_perp = c2 - c1 * (float(np.dot(c1, c2)) / c1_norm_sq)
            transverse[index] = float(np.linalg.norm(c2_perp))
        coefficients, _, _, _ = np.linalg.lstsq(matrix, target, rcond=None)
        alpha[index] = float(coefficients[0])
        beta[index] = float(coefficients[1])
        residuals.append(matrix @ coefficients - target)
        if rank[index] == 2:
            condition[index] = float(singular[0] / singular[-1])

    target_active = receipt.target_stress[active]
    residual_active = np.asarray(residuals, dtype=float)
    relative_residual = _vector_rms(residual_active) / max(
        _vector_rms(target_active), np.finfo(float).tiny
    )
    rank2_needed = rank[needed] == 2
    budget_limit = budget * (1.0 + BUDGET_ROUNDOFF_RELATIVE_TOLERANCE) + np.finfo(float).eps
    second_bounded = np.abs(beta[needed]) <= budget_limit
    both_bounded = (
        (np.abs(alpha[needed]) <= budget_limit) & second_bounded & rank2_needed
    )
    finite_conditions = condition[needed & np.isfinite(condition)]

    return {
        "coefficient_budget": budget,
        "budget_roundoff_relative_tolerance": BUDGET_ROUNDOFF_RELATIVE_TOLERANCE,
        "active_target_nodes": int(np.count_nonzero(active)),
        "nodes_requiring_second_direction": int(np.count_nonzero(needed)),
        "rank2_required_nodes": int(np.count_nonzero(rank2_needed)),
        "rank_deficient_required_nodes": int(np.count_nonzero(~rank2_needed)),
        "second_coefficient_within_budget_nodes": int(np.count_nonzero(second_bounded & rank2_needed)),
        "both_coefficients_within_budget_nodes": int(np.count_nonzero(both_bounded)),
        "max_abs_first_signed_coefficient_on_required": float(np.max(np.abs(alpha[needed]))),
        "max_abs_second_signed_coefficient_on_required": float(np.max(np.abs(beta[needed]))),
        "min_transverse_response_on_required": float(np.min(transverse[needed])),
        "max_transverse_response_on_required": float(np.max(transverse[needed])),
        "min_determinant_abs_on_required": float(np.min(determinant[needed])),
        "condition_max_on_required": float(np.max(finite_conditions)) if len(finite_conditions) else None,
        "condition_median_on_required": float(np.median(finite_conditions)) if len(finite_conditions) else None,
        "two_column_relative_stress_residual_rms": float(relative_residual),
        "all_required_nodes_rank2": bool(np.all(rank2_needed)),
        "second_signed_update_within_budget_on_all_required": bool(np.all(second_bounded & rank2_needed)),
        "both_signed_updates_within_budget_on_all_required": bool(np.all(both_bounded)),
        "algebraic_stress_fit_within_tolerance": bool(
            relative_residual <= algebraic_relative_residual_tolerance
        ),
        "bounded_inverse_preflight_passed": bool(
            np.all(rank2_needed)
            and np.all(both_bounded)
            and relative_residual <= algebraic_relative_residual_tolerance
        ),
        "interpretation": (
            "Covariance/stress inverse preflight only. A passing map must still be materialized through public velocity "
            "and pass frozen held-in/held-out correction-cycle residual checks."
        ),
    }


def _receipt_from_target_report(report: dict[str, Any]) -> MissingCovarianceColumnTarget:
    target = report["missing_second_column_target"]
    inputs = report["inputs"]
    return build_missing_covariance_column_target(
        np.asarray(target["radii"], dtype=float),
        np.asarray(target["target_stress_theta_axial"], dtype=float),
        np.asarray(target["current_response_theta_axial"], dtype=float),
        target_floor_fraction=float(inputs["target_floor_fraction"]),
        response_floor_fraction=float(inputs["response_floor_fraction"]),
    )


def generate_actual_core_report(
    *,
    output: str | Path = "artifacts/kokuno_agent3/second_column_bounded_inverse_report.json",
    target_output: str | Path = "artifacts/kokuno_agent3/second_column_bounded_inverse_target.json",
) -> dict[str, Any]:
    """Derive bounded inverse requirements from the real #282 stress target."""
    target_report = generate_missing_target_report(output=target_output)
    receipt = _receipt_from_target_report(target_report)
    symmetric_budget = current_signed_coefficient_budget(receipt)
    routed_amplitude_budget = abs(float(ROUTED_OSCILLATORY_AMPLITUDE))
    symmetric_envelope = required_second_column_envelope(receipt, coefficient_budget=symmetric_budget)
    routed_envelope = required_second_column_envelope(receipt, coefficient_budget=routed_amplitude_budget)
    duplicate = evaluate_bounded_second_column(receipt, receipt.current_response, coefficient_budget=symmetric_budget)
    raw_missing = evaluate_bounded_second_column(receipt, receipt.missing_response, coefficient_budget=symmetric_budget)
    scaled_missing = evaluate_bounded_second_column(
        receipt,
        receipt.missing_response / symmetric_budget,
        coefficient_budget=symmetric_budget,
    )

    report: dict[str, Any] = {
        "task": TASK,
        "source": {
            "reader": SOURCE_READER,
            "edition": SOURCE_EDITION,
            "doi": SOURCE_DOI,
            "repository": SOURCE_REPOSITORY,
            "commit": SOURCE_COMMIT,
            "path": SOURCE_PATH,
            "section": SOURCE_SECTION,
            "source_structure": (
                "The source signed covariance step uses an invertible two-column response map H and signed differential inverse. "
                "The determinant/transverse-response floors below are engineering consequences of that inversion, not source numerical constants."
            ),
        },
        "inputs": {
            "leading_candidate_sha256": target_report["inputs"]["leading_candidate_sha256"],
            "oscillatory_amplitude": target_report["inputs"]["oscillatory_amplitude"],
            "oscillatory_phase": target_report["inputs"]["oscillatory_phase"],
            "profile_time": target_report["inputs"]["profile_time"],
            "profile_z": target_report["inputs"]["profile_z"],
            "profile_annulus": target_report["inputs"]["profile_annulus"],
            "surrogate_defect_used": False,
        },
        "measured_target_summary": target_report["missing_second_column_target"],
        "bounded_inverse": {
            "symmetric_measured_coefficient_budget": symmetric_budget,
            "symmetric_budget_definition": (
                "max |delta a_1| from the existing real rank-one target projection; data-bound engineering scale"
            ),
            "routed_primary_amplitude_budget": routed_amplitude_budget,
            "routed_primary_amplitude_budget_definition": (
                "absolute routed Agent-2 primary amplitude; broad engineering scale, not a source smallness theorem"
            ),
            "symmetric_budget_envelope": symmetric_envelope,
            "routed_amplitude_envelope": routed_envelope,
            "calibration": {
                "duplicate_existing_column": duplicate,
                "raw_orthogonal_remainder_as_second_response": raw_missing,
                "orthogonal_remainder_scaled_to_symmetric_budget": scaled_missing,
                "raw_or_scaled_remainder_used_as_candidate_velocity": False,
            },
        },
        "st006_cross_route_reference": {
            "candidate_sha256": ST006_SHA256,
            "momentum_max": ST006_MOMENTUM_MAX,
            "volume_l2": ST006_VOLUME_L2,
            "validation_seed": ST006_VALIDATION_SEED,
            "validation_points": ST006_VALIDATION_POINTS,
            "finest_spatial_step": ST006_FINEST_SPATIAL_STEP,
            "directly_comparable": False,
            "reason": "Local covariance/stress inverse scale is not ST006's full-domain normalized momentum protocol.",
        },
        "routing": {
            "second_public_covariance_column_available": False,
            "finite_correction_cycle_rerun_allowed": False,
            "next_required": (
                "When Agent 2 materializes a genuinely independent public covariance column, evaluate its measured response here. "
                "Only a bounded inverse may proceed to public-velocity materialization and frozen held-in/held-out residual tests."
            ),
        },
        "truth_boundary": {
            "real_candidate_defect_consumed": True,
            "surrogate_defect_used": False,
            "real_missing_covariance_target_consumed": True,
            "second_column_scale_requirement_measured": True,
            "second_public_covariance_column_available": False,
            "new_oscillatory_column_constructed": False,
            "agent2_complete_curl_reimplemented": False,
            "public_velocity_correction_materialized": False,
            "finite_correction_cycle_run": False,
            "residual_reduction_claimed": False,
            "formal_full_domain_normalized_pde_gate_assessed": False,
            "normalized_ns_residual_le_1e_minus_3_claimed": False,
            "pde_validated": False,
            "paper_exact": False,
            "openai_field_identified": False,
            "blowup_proved": False,
        },
    }
    path = Path(output)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    return report


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", default="artifacts/kokuno_agent3/second_column_bounded_inverse_report.json")
    parser.add_argument("--target-output", default="artifacts/kokuno_agent3/second_column_bounded_inverse_target.json")
    return parser


def main() -> None:
    args = _parser().parse_args()
    report = generate_actual_core_report(output=args.output, target_output=args.target_output)
    bounded = report["bounded_inverse"]
    summary = report["measured_target_summary"]
    print(
        json.dumps(
            {
                "task": report["task"],
                "missing_relative_vector_rms": summary["missing_relative_vector_rms"],
                "nodes_requiring_second_direction": summary["nodes_requiring_second_direction"],
                "symmetric_measured_coefficient_budget": bounded["symmetric_measured_coefficient_budget"],
                "required_transverse_over_current_max": bounded["symmetric_budget_envelope"]["required_transverse_over_current_max"],
                "raw_missing_preflight_passed": bounded["calibration"]["raw_orthogonal_remainder_as_second_response"]["bounded_inverse_preflight_passed"],
                "scaled_missing_preflight_passed": bounded["calibration"]["orthogonal_remainder_scaled_to_symmetric_budget"]["bounded_inverse_preflight_passed"],
                "finite_correction_cycle_rerun_allowed": report["routing"]["finite_correction_cycle_rerun_allowed"],
            },
            indent=2,
            allow_nan=False,
        )
    )


if __name__ == "__main__":
    main()
