"""Bound the second covariance column before another Kokuno correction cycle.

Agent-3 PR #282 exposed the real two-component compact-stress remainder that
cannot be represented by the current Agent-2 covariance column.  In two
components, any genuinely independent second column spans the missing direction
algebraically, but that fact alone is not sufficient: a nearly parallel or very
small column can require an arbitrarily large signed amplitude update.

The corrected Kokuno reader uses a genuinely invertible two-column covariance
map ``H`` and a signed differential inverse.  This module turns the scale part
of that requirement into a fail-closed engineering preflight.  For current
response ``c1``, future response ``c2`` and the measured orthogonal remainder
``m``, only

    c2_perp = c2 - proj_{c1}(c2)

can remove ``m``.  The exact second signed coefficient obeys

    |delta a_2| = |m| / |c2_perp|

at a rank-two node, while

    |det[c1,c2]| = |c1| |c2_perp|.

Consequently a declared coefficient budget ``B`` implies the necessary local
response floor

    |c2_perp| >= |m| / B

and determinant floor ``|det H| >= |c1| |m| / B``.  These are algebraic
preflight bounds, not Navier--Stokes residual contraction claims and not a
replacement for the later held-in/held-out finite-cycle test.
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
    generate_actual_core_report,
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
    """Return a data-bound symmetric budget from the existing rank-one solve.

    This is not a source constant.  It is the largest already-measured signed
    coefficient required by the current column on safe active nodes, and is
    used only as an engineering reference scale for a future second column.
    """
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
    """Compute the transverse-response floor required by a signed-amplitude budget."""
    if not isinstance(receipt, MissingCovarianceColumnTarget):
        raise TypeError("receipt must be MissingCovarianceColumnTarget")
    budget = float(coefficient_budget)
    if not np.isfinite(budget) or budget <= 0.0:
        raise ValueError("coefficient_budget must be positive and finite")

    needed = _second_needed_mask(receipt)
    indices = np.flatnonzero(needed)
    if len(indices) == 0:
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
        "nodes_requiring_second_direction": int(len(indices)),
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
            "Necessary local scale for keeping the future second signed coefficient within the declared budget. "
            "It does not guarantee a Navier-Stokes residual reduction."
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
    """Evaluate rank, stress coverage and signed-coefficient size for a second column."""
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
    active_indices = np.flatnonzero(active)
    target_active = receipt.target_stress[active]
    residual_active = np.zeros_like(target_active)

    rank_by_node = np.zeros(len(receipt.radii), dtype=int)
    condition_by_node = np.full(len(receipt.radii), np.nan, dtype=float)
    determinant_by_node = np.zeros(len(receipt.radii), dtype=float)
    alpha_by_node = np.zeros(len(receipt.radii), dtype=float)
    beta_by_node = np.zeros(len(receipt.radii), dtype=float)
    transverse_by_node = np.zeros(len(receipt.radii), dtype=float)

    for output_index, source_index in enumerate(active_indices):
        c1 = receipt.current_response[source_index]
        c2 = second[source_index]
        target = receipt.target_stress[source_index]
        matrix = np.column_stack((c1, c2))
        singular = np.linalg.svd(matrix, compute_uv=False)
        if singular[0] <= np.finfo(float).tiny:
            rank = 0
        else:
            rank = int(np.count_nonzero(singular > relative_rank_tolerance * singular[0]))
        rank_by_node[source_index] = rank
        determinant_by_node[source_index] = float(abs(np.linalg.det(matrix)))

        c1_norm_sq = float(np.dot(c1, c1))
        if c1_norm_sq > np.finfo(float).tiny:
            c2_perp = c2 - c1 * (float(np.dot(c1, c2)) / c1_norm_sq)
            transverse_by_node[source_index] = float(np.linalg.norm(c2_perp))

        coefficients, _, _, _ = np.linalg.lstsq(matrix, target, rcond=None)
        alpha_by_node[source_index] = float(coefficients[0])
        beta_by_node[source_index] = float(coefficients[1])
        residual_active[output_index] = matrix @ coefficients - target
        if rank == 2:
            condition_by_node[source_index] = float(singular[0] / singular[-1])

    target_scale = max(_vector_rms(target_active), np.finfo(float).tiny)
    residual_relative = _vector_rms(residual_active) / target_scale
    needed_indices = np.flatnonzero(needed)
    rank2_needed = rank_by_node[needed] == 2
    bounded_beta_needed = np.abs(beta_by_node[needed]) <= budget
    bounded_both_needed = (
        (np.abs(alpha_by_node[needed]) <= budget)
        & bounded_beta_needed
        & rank2_needed
    )
    finite_conditions = condition_by_node[needed & np.isfinite(condition_by_node)]

    return {
        "coefficient_budget": budget,
        "active_target_nodes": int(np.count_nonzero(active)),
        "nodes_requiring_second_direction": int(len(needed_indices)),
        "rank2_required_nodes": int(np.count_nonzero(rank2_needed)),
        "rank_deficient_required_nodes": int(len(needed_indices) - np.count_nonzero(rank2_needed)),
        "second_coefficient_within_budget_nodes": int(np.count_nonzero(bounded_beta_needed & rank2_needed)),
        "both_coefficients_within_budget_nodes": int(np.count_nonzero(bounded_both_needed)),
        "max_abs_first_signed_coefficient_on_required": float(np.max(np.abs(alpha_by_node[needed]))),
        "max_abs_second_signed_coefficient_on_required": float(np.max(np.abs(beta_by_node[needed]))),
        "min_transverse_response_on_required": float(np.min(transverse_by_node[needed])),
        "max_transverse_response_on_required": float(np.max(transverse_by_node[needed])),
        "min_determinant_abs_on_required": float(np.min(determinant_by_node[needed])),
        "condition_max_on_required": float(np.max(finite_conditions)) if len(finite_conditions) else None,
        "condition_median_on_required": float(np.median(finite_conditions)) if len(finite_conditions) else None,
        "two_column_relative_stress_residual_rms": float(residual_relative),
        "all_required_nodes_rank2": bool(np.all(rank2_needed)),
        "second_signed_update_within_budget_on_all_required": bool(np.all(bounded_beta_needed & rank2_needed)),
        "both_signed_updates_within_budget_on_all_required": bool(np.all(bounded_both_needed)),
        "algebraic_stress_fit_within_tolerance": bool(
            residual_relative <= algebraic_relative_residual_tolerance
        ),
        "bounded_inverse_preflight_passed": bool(
            np.all(rank2_needed)
            and np.all(bounded_both_needed)
            and residual_relative <= algebraic_relative_residual_tolerance
        ),
        "interpretation": (
            "This is only a covariance/stress inverse preflight. Passing it does not authorize a PDE-valid claim; "
            "the resulting public velocity must still pass the frozen held-in/held-out correction-cycle checks."
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
    """Measure the real target and derive scale requirements for Agent 2's future column."""
    target_report = generate_actual_core_report.__globals__["generate_actual_core_report"]
    # The globals lookup above would recurse after this function is defined.  Keep
    # the imported generator under an explicit local alias instead.
    raise RuntimeError("internal alias was not initialized")


def _generate_actual_core_report_impl(
    target_generator,
    *,
    output: str | Path,
    target_output: str | Path,
) -> dict[str, Any]:
    target_report = target_generator(output=target_output)
    receipt = _receipt_from_target_report(target_report)
    symmetric_budget = current_signed_coefficient_budget(receipt)
    routed_amplitude_budget = abs(float(ROUTED_OSCILLATORY_AMPLITUDE))

    symmetric_envelope = required_second_column_envelope(
        receipt, coefficient_budget=symmetric_budget
    )
    routed_envelope = required_second_column_envelope(
        receipt, coefficient_budget=routed_amplitude_budget
    )

    duplicate = evaluate_bounded_second_column(
        receipt,
        receipt.current_response,
        coefficient_budget=symmetric_budget,
    )
    raw_missing = evaluate_bounded_second_column(
        receipt,
        receipt.missing_response,
        coefficient_budget=symmetric_budget,
    )
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
                "Kokuno's signed covariance step uses an invertible two-column response map H and a signed differential inverse. "
                "The local determinant/transverse-response bounds reported here are engineering consequences of that inversion, "
                "not source numerical constants."
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
            "coefficient_budget_source": {
                "symmetric_measured_budget": (
                    "maximum absolute current-column signed coefficient on safe active nodes"
                ),
                "routed_amplitude_budget": (
                    "absolute routed Agent-2 primary amplitude; broad engineering scale only"
                ),
            },
        },
        "measured_target_summary": target_report["missing_second_column_target"],
        "bounded_inverse": {
            "symmetric_measured_coefficient_budget": symmetric_budget,
            "routed_primary_amplitude_budget": routed_amplitude_budget,
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
            "reason": (
                "This report measures a local covariance/stress inverse scale, not the ST006 full-domain normalized momentum protocol."
            ),
        },
        "routing": {
            "second_public_covariance_column_available": False,
            "finite_correction_cycle_rerun_allowed": False,
            "next_required": (
                "When Agent 2 materializes a genuinely independent public covariance column, evaluate its measured response here. "
                "Reject a rank-two direction if the required signed coefficients are unbounded at the preregistered scale; only a "
                "bounded covariance inverse may proceed to public-velocity materialization and frozen held-in/held-out residual tests."
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
            "normalized_ns_residual_le_1e3_claimed": False,
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


# Preserve an unambiguous alias to the #282 real-target generator, then replace
# the public report function above with the actual wrapper.  Keeping the alias
# explicit avoids copying the real defect/stress reconstruction machinery.
_target_report_generator = globals().pop("generate_actual_core_report")


def generate_actual_core_report(
    *,
    output: str | Path = "artifacts/kokuno_agent3/second_column_bounded_inverse_report.json",
    target_output: str | Path = "artifacts/kokuno_agent3/second_column_bounded_inverse_target.json",
) -> dict[str, Any]:
    return _generate_actual_core_report_impl(
        _target_report_generator,
        output=output,
        target_output=target_output,
    )


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        default="artifacts/kokuno_agent3/second_column_bounded_inverse_report.json",
    )
    parser.add_argument(
        "--target-output",
        default="artifacts/kokuno_agent3/second_column_bounded_inverse_target.json",
    )
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
