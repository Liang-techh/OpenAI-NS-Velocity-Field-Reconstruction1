"""Fail-closed bounded inverse for a multi-label Kokuno covariance family.

Agent 3 PR #302 bounds a future *single* second covariance column before another
finite correction cycle.  Agent 3 PR #382 then exposes the complete per-label
covariance Jacobian of an overlapping Agent-2 slow-label family.  This module
joins those interfaces without changing the real defect target or the existing
signed-update scale.

For one radial node let ``c0`` be the existing covariance response and let
``J=[c1,...,cm]`` contain additional family responses.  The repository solves

    [c0 J] da ~= target

with the minimum-l2 least-squares solution only as an algebraic preflight.  A
rank-two span and a tiny algebraic residual are necessary but not sufficient for
a useful correction.  In particular, many duplicate labels must not evade the
signed-update budget by splitting one large correction into many individually
small coefficients.  Therefore the additional family is conservatively guarded
by

    |da0| <= B0,        sum_j |daj| <= B_family.

For exactly one additional column and ``B0=B_family`` this is the same
per-coefficient budget rule used by PR #302.  For multiple labels the l1 guard is
an autonomous anti-budget-laundering engineering rule, not a Kokuno source
constant or a Navier--Stokes acceptance threshold.

The source only motivates the signed covariance inverse structure.  Multi-label
minimum-norm solving, the l1 family budget, floating tolerances and unit-mapping
guards are repository machinery.  No source coefficient normalization is
inferred here: callers must explicitly establish that the supplied coefficient
budgets use the same units as the covariance Jacobian before the preflight can
pass.
"""
from __future__ import annotations

import argparse
from dataclasses import replace
import json
from pathlib import Path
from typing import Any, Sequence

import numpy as np

from .kokuno_complete_curl import KokunoCompleteCurlCorrection
from .kokuno_covariance_tangent_preflight import complete_curl_phase_evaluator
from .kokuno_family_covariance_cross_terms import (
    measure_family_covariance_jacobian,
    single_column_family_evaluator,
)
from .kokuno_missing_covariance_column_target import (
    MissingCovarianceColumnTarget,
    build_missing_covariance_column_target,
    generate_actual_core_report as generate_missing_target_report,
)
from .kokuno_second_column_bounded_inverse import (
    ALGEBRAIC_RELATIVE_RESIDUAL_TOLERANCE,
    BUDGET_ROUNDOFF_RELATIVE_TOLERANCE,
    RANK_RELATIVE_TOLERANCE,
    current_signed_coefficient_budget,
)
from .kokuno_signed_covariance_inverse import (
    ROUTED_OSCILLATORY_AMPLITUDE,
    ROUTED_OSCILLATORY_PHASE,
)

TASK = "KOKUNO-A3-FAMILY-BOUNDED-INVERSE-019"
SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
SOURCE_EDITION = "2026.09.09-consolidated"
SOURCE_DOI = "10.5281/zenodo.22678406"
AGENT2_FAMILY_PR = 380
AGENT2_FAMILY_HEAD = "43f021e904dbb09c9bc2aba1d9ea2954083937e1"
AGENT3_FAMILY_PR = 382
AGENT3_FAMILY_HEAD = "c80d3d0993c0cc62e2f24cab5d427c5fcb12138d"


def _vector_rms(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    if values.ndim != 2 or values.shape[1] != 2 or not np.isfinite(values).all():
        raise ValueError("values must be finite with shape (N,2)")
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


def evaluate_family_bounded_inverse(
    receipt: MissingCovarianceColumnTarget,
    additional_responses: np.ndarray,
    *,
    coefficient_labels: Sequence[str],
    current_coefficient_budget: float,
    additional_family_l1_budget: float,
    budget_unit_matches_jacobian: bool,
    coefficient_unit: str,
    relative_rank_tolerance: float = RANK_RELATIVE_TOLERANCE,
    algebraic_relative_residual_tolerance: float = ALGEBRAIC_RELATIVE_RESIDUAL_TOLERANCE,
) -> dict[str, Any]:
    """Screen an arbitrary number of additional covariance response columns.

    ``additional_responses`` has shape ``(radial_nodes, labels, 2)``.  The
    existing response in ``receipt`` is always retained as column zero.  The
    least-squares coefficients are diagnostic only; a passing result additionally
    requires an explicit unit match and the aggregate l1 bound on all additional
    labels.
    """
    if not isinstance(receipt, MissingCovarianceColumnTarget):
        raise TypeError("receipt must be MissingCovarianceColumnTarget")
    responses = np.asarray(additional_responses, dtype=float)
    if responses.ndim != 3 or responses.shape[0] != len(receipt.radii) or responses.shape[2] != 2:
        raise ValueError("additional_responses must have shape (radial_nodes, labels, 2)")
    if responses.shape[1] < 1 or not np.isfinite(responses).all():
        raise ValueError("additional_responses must contain at least one finite label column")

    labels = tuple(str(label) for label in coefficient_labels)
    if len(labels) != responses.shape[1] or len(set(labels)) != len(labels) or any(not x for x in labels):
        raise ValueError("coefficient_labels must be unique, nonempty and match the label count")

    current_budget = float(current_coefficient_budget)
    family_budget = float(additional_family_l1_budget)
    if not np.isfinite([current_budget, family_budget]).all() or current_budget <= 0.0 or family_budget <= 0.0:
        raise ValueError("coefficient budgets must be positive and finite")
    if not isinstance(budget_unit_matches_jacobian, (bool, np.bool_)):
        raise TypeError("budget_unit_matches_jacobian must be boolean")
    coefficient_unit = str(coefficient_unit)
    if not coefficient_unit:
        raise ValueError("coefficient_unit must be nonempty")
    if not 0.0 < relative_rank_tolerance < 1.0e-2:
        raise ValueError("relative_rank_tolerance must lie in (0,1e-2)")
    if not 0.0 < algebraic_relative_residual_tolerance < 1.0e-4:
        raise ValueError("algebraic_relative_residual_tolerance must lie in (0,1e-4)")

    active = receipt.active_target_mask
    needed = _second_needed_mask(receipt)
    if not np.any(active) or not np.any(needed):
        raise ValueError("receipt must contain active nodes requiring a second direction")

    count = len(receipt.radii)
    label_count = responses.shape[1]
    rank = np.zeros(count, dtype=int)
    condition = np.full(count, np.nan, dtype=float)
    current_coeff = np.zeros(count, dtype=float)
    additional_coeff = np.zeros((count, label_count), dtype=float)
    additional_l1 = np.zeros(count, dtype=float)
    additional_linf = np.zeros(count, dtype=float)
    smallest_singular = np.zeros(count, dtype=float)
    residual_by_node = np.zeros((count, 2), dtype=float)

    for index in np.flatnonzero(active):
        matrix = np.concatenate(
            (receipt.current_response[index, :, None], responses[index].T),
            axis=1,
        )
        singular = np.linalg.svd(matrix, compute_uv=False)
        if singular[0] > np.finfo(float).tiny:
            rank[index] = int(np.count_nonzero(singular > relative_rank_tolerance * singular[0]))
        if len(singular) >= 2:
            smallest_singular[index] = float(singular[-1])
        if rank[index] == 2 and singular[-1] > np.finfo(float).tiny:
            condition[index] = float(singular[0] / singular[-1])
        coeffs, _, _, _ = np.linalg.lstsq(matrix, receipt.target_stress[index], rcond=None)
        current_coeff[index] = float(coeffs[0])
        additional_coeff[index] = coeffs[1:]
        additional_l1[index] = float(np.sum(np.abs(coeffs[1:])))
        additional_linf[index] = float(np.max(np.abs(coeffs[1:])))
        residual_by_node[index] = matrix @ coeffs - receipt.target_stress[index]

    target_active = receipt.target_stress[active]
    residual_active = residual_by_node[active]
    relative_residual = _vector_rms(residual_active) / max(
        _vector_rms(target_active), np.finfo(float).tiny
    )

    roundoff = BUDGET_ROUNDOFF_RELATIVE_TOLERANCE
    current_limit = current_budget * (1.0 + roundoff) + np.finfo(float).eps
    family_limit = family_budget * (1.0 + roundoff) + np.finfo(float).eps
    rank2_needed = rank[needed] == 2
    current_bounded = np.abs(current_coeff[needed]) <= current_limit
    family_bounded = additional_l1[needed] <= family_limit
    finite_condition = condition[needed & np.isfinite(condition)]

    algebraic_pass = bool(relative_residual <= algebraic_relative_residual_tolerance)
    rank_pass = bool(np.all(rank2_needed))
    current_budget_pass = bool(np.all(current_bounded & rank2_needed))
    family_budget_pass = bool(np.all(family_bounded & rank2_needed))
    units_pass = bool(budget_unit_matches_jacobian)
    passed = bool(
        units_pass
        and rank_pass
        and current_budget_pass
        and family_budget_pass
        and algebraic_pass
    )

    return {
        "coefficient_labels": list(labels),
        "coefficient_unit": coefficient_unit,
        "budget_unit_matches_jacobian": units_pass,
        "current_coefficient_budget": current_budget,
        "additional_family_l1_budget": family_budget,
        "budget_roundoff_relative_tolerance": roundoff,
        "active_target_nodes": int(np.count_nonzero(active)),
        "nodes_requiring_second_direction": int(np.count_nonzero(needed)),
        "rank2_required_nodes": int(np.count_nonzero(rank2_needed)),
        "rank_deficient_required_nodes": int(np.count_nonzero(~rank2_needed)),
        "max_abs_current_coefficient_on_required": float(np.max(np.abs(current_coeff[needed]))),
        "max_additional_coefficient_linf_on_required": float(np.max(additional_linf[needed])),
        "max_additional_coefficient_l1_on_required": float(np.max(additional_l1[needed])),
        "max_family_l1_over_budget_on_required": float(np.max(additional_l1[needed]) / family_budget),
        "min_smallest_singular_value_on_required": float(np.min(smallest_singular[needed])),
        "condition_max_on_required": float(np.max(finite_condition)) if len(finite_condition) else None,
        "condition_median_on_required": float(np.median(finite_condition)) if len(finite_condition) else None,
        "algebraic_relative_stress_residual_rms": float(relative_residual),
        "all_required_nodes_rank2": rank_pass,
        "current_update_within_budget_on_all_required": current_budget_pass,
        "additional_family_l1_within_budget_on_all_required": family_budget_pass,
        "algebraic_stress_fit_within_tolerance": algebraic_pass,
        "family_bounded_inverse_preflight_passed": passed,
        "anti_budget_laundering_rule": (
            "sum(abs(additional label coefficients)) <= additional_family_l1_budget; "
            "for one additional label this reduces to the PR #302 second-coefficient bound"
        ),
        "interpretation": (
            "Algebraic covariance/stress inverse preflight only. Passing still requires a materialized public "
            "velocity correction, retained radial d_z sigma_1, and frozen held-in/held-out NS residual checks."
        ),
    }


def generate_actual_core_report(
    *,
    output: str | Path = "artifacts/kokuno_agent3/family_bounded_inverse_report.json",
    target_output: str | Path = "artifacts/kokuno_agent3/family_bounded_inverse_target.json",
) -> dict[str, Any]:
    """Run a real-defect duplicate-column calibration through the family guard."""
    target_report = generate_missing_target_report(output=target_output)
    receipt = _receipt_from_target_report(target_report)

    unit = replace(
        KokunoCompleteCurlCorrection(),
        amplitude=1.0,
        phase=ROUTED_OSCILLATORY_PHASE,
    )
    family = single_column_family_evaluator(
        complete_curl_phase_evaluator(unit),
        physical_amplitude=ROUTED_OSCILLATORY_AMPLITUDE,
        label="existing_complete_curl",
    )
    measurement = measure_family_covariance_jacobian(
        family,
        receipt.radii,
        coefficient_tangent_scales=[1.0 / ROUTED_OSCILLATORY_AMPLITUDE],
        time=float(target_report["inputs"]["profile_time"]),
        z=float(target_report["inputs"]["profile_z"]),
    )
    jacobian = np.asarray(measurement["per_label_covariance_jacobian_theta_axial"], dtype=float)
    current = receipt.current_response
    calibration_relative_error = _vector_rms(jacobian[:, 0, :] - current) / max(
        _vector_rms(current), np.finfo(float).tiny
    )

    frozen_budget = current_signed_coefficient_budget(receipt)
    calibration = evaluate_family_bounded_inverse(
        receipt,
        jacobian,
        coefficient_labels=measurement["beta_labels"],
        current_coefficient_budget=frozen_budget,
        additional_family_l1_budget=frozen_budget,
        budget_unit_matches_jacobian=True,
        coefficient_unit="existing_complete_curl_physical_amplitude",
    )

    report: dict[str, Any] = {
        "task": TASK,
        "source": {
            "repository": SOURCE_REPOSITORY,
            "commit": SOURCE_COMMIT,
            "path": SOURCE_PATH,
            "edition": SOURCE_EDITION,
            "doi": SOURCE_DOI,
            "source_structure": (
                "The corrected reader supplies the signed covariance inverse structure. The N-column minimum-norm "
                "preflight, aggregate l1 budget and coefficient-unit guard in this report are repository engineering."
            ),
        },
        "handoff": {
            "agent2_family_pr": AGENT2_FAMILY_PR,
            "agent2_family_head": AGENT2_FAMILY_HEAD,
            "agent3_cross_term_pr": AGENT3_FAMILY_PR,
            "agent3_cross_term_head": AGENT3_FAMILY_HEAD,
        },
        "inputs": {
            "leading_candidate_sha256": target_report["inputs"]["leading_candidate_sha256"],
            "oscillatory_amplitude": target_report["inputs"]["oscillatory_amplitude"],
            "oscillatory_phase": target_report["inputs"]["oscillatory_phase"],
            "profile_time": target_report["inputs"]["profile_time"],
            "profile_z": target_report["inputs"]["profile_z"],
            "profile_annulus": target_report["inputs"]["profile_annulus"],
            "surrogate_defect_used": False,
            "frozen_signed_coefficient_budget": frozen_budget,
            "budget_changed": False,
        },
        "real_target_summary": target_report["missing_second_column_target"],
        "existing_one_label_family_calibration": {
            "coefficient_tangent_scale": 1.0 / ROUTED_OSCILLATORY_AMPLITUDE,
            "relative_error_to_existing_current_response": float(calibration_relative_error),
            "cross_covariance_relative_vector_rms": measurement["cross_covariance_relative_vector_rms"],
            "family_inverse": calibration,
            "expected_result": (
                "REJECT: the existing family Jacobian duplicates the current covariance direction and cannot "
                "supply the missing transverse stress response."
            ),
        },
        "routing": {
            "actual_source_multilabel_family_consumed": False,
            "source_coefficient_unit_mapping_available": False,
            "actual_source_family_bounded_inverse_assessed": False,
            "finite_correction_cycle_rerun_allowed": False,
            "next_required_handoff": (
                "Agent 2 must supply a concrete source-motivated by-beta public velocity family plus explicit "
                "coefficient normalization. Agent 3 must then measure the full cross-term-aware Jacobian, establish "
                "budget units, and rerun this guard before materializing a correction."
            ),
        },
        "truth_boundary": {
            "real_candidate_defect_target_consumed": True,
            "surrogate_defect_used": False,
            "family_bounded_inverse_contract_executable": True,
            "existing_one_label_duplicate_rejected": not calibration["family_bounded_inverse_preflight_passed"],
            "new_oscillatory_column_constructed": False,
            "public_velocity_correction_materialized": False,
            "finite_correction_cycle_run": False,
            "residual_reduction_claimed": False,
            "formal_full_domain_normalized_pde_gate_assessed": False,
            "pde_validated": False,
            "paper_exact": False,
            "openai_field_identified": False,
            "blowup_proved": False,
        },
    }

    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        default="artifacts/kokuno_agent3/family_bounded_inverse_report.json",
    )
    parser.add_argument(
        "--target-output",
        default="artifacts/kokuno_agent3/family_bounded_inverse_target.json",
    )
    args = parser.parse_args()
    report = generate_actual_core_report(output=args.output, target_output=args.target_output)
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
