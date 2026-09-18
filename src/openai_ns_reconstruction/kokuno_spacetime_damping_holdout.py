"""Validate one frozen Kokuno covariance damping step on disjoint spacetime cells.

Agent-3 PR #422 selected one shared scalar damping on a 3x3 real-defect
spacetime screen.  Those cells are therefore calibration data for that scalar,
not independent validation data.  This module freezes the selected update and
its damping, then evaluates them without retuning on a disjoint spacetime
screen.

The corrected Kokuno reader motivates retaining nonlinear correction
self-covariance after the signed covariance inverse.  The calibration/validation
split and fresh screen are repository validation engineering; they are not new
Kokuno formulas and do not alter any PDE or algebraic threshold.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Iterable, Mapping

import numpy as np

from .kokuno_bounded_coordinate_covariance import (
    EXPECTED_PARAMETER_NAMES,
    EXPECTED_PARAMETER_UNIT,
    _one_label_real_family,
    convert_amplitude_budget_to_fractional,
    measure_bounded_coordinate_covariance,
)
from .kokuno_bounded_family_coefficients import KokunoBoundedFamilyCoefficientCoordinates
from .kokuno_family_bounded_inverse import evaluate_family_bounded_inverse
from .kokuno_missing_covariance_column_target import (
    DEFAULT_ANGULAR_COUNT,
    DEFAULT_PHASE_COUNT,
    DEFAULT_RADIAL_COUNT,
    MissingCovarianceColumnTarget,
)
from .kokuno_multislice_covariance_preflight import (
    HELD_OUT_CYCLE_TIMES,
    REFERENCE_BUDGET_TIME,
    build_real_covariance_slice,
)
from .kokuno_quadratic_covariance_gain_guard import (
    QUADRATIC_IDENTITY_RELATIVE_TOLERANCE,
    measure_bounded_coordinate_quadratic_change,
)
from .kokuno_second_column_bounded_inverse import (
    ALGEBRAIC_RELATIVE_RESIDUAL_TOLERANCE,
    current_signed_coefficient_budget,
)
from .kokuno_signed_covariance_inverse import (
    PROFILE_ANNULUS,
    ROUTED_OSCILLATORY_AMPLITUDE,
    ROUTED_OSCILLATORY_PHASE,
)
from .kokuno_spacetime_covariance_preflight import REFERENCE_BUDGET_Z, SPATIAL_SCREEN_Z
from .kokuno_spacetime_damped_quadratic_gain import solve_shared_spacetime_damping

TASK = "KOKUNO-A3-DAMPING-HOLDOUT-TRANSFER-024"
SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
SOURCE_EDITION = "2026.09.09-consolidated"
SOURCE_DOI = "10.5281/zenodo.22678406"
SOURCE_SECTION = "signed covariance inverse; nonlinear self-covariance retained"
AGENT3_CALIBRATION_PR = 422
AGENT3_BOUNDED_INVERSE_PR = 393
AGENT3_RADIAL_SIDE_EFFECT_PR = 370

# Interleaved in z and outside the #422 calibration times.  These values are
# fixed repository validation choices, not source constants.
VALIDATION_TIMES = (0.3125, 0.6875)
VALIDATION_Z = (0.07, 0.09)
FROZEN_COMMON_UPDATE_SIGN = -1.0


def _vector_rms(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    if values.ndim != 2 or values.shape[1] != 2 or not np.isfinite(values).all():
        raise ValueError("values must be finite with shape (N,2)")
    return float(np.sqrt(np.mean(np.sum(values * values, axis=1))))


def _validate_axis(values: Iterable[float], *, name: str) -> tuple[float, ...]:
    result = tuple(float(value) for value in values)
    if len(result) < 1 or not np.isfinite(result).all():
        raise ValueError(f"{name} must contain at least one finite value")
    if any(right <= left for left, right in zip(result, result[1:])):
        raise ValueError(f"{name} must be strictly increasing")
    return result


def _fixed_relative_residual(
    target: np.ndarray,
    active: np.ndarray,
    linear: np.ndarray,
    quadratic: np.ndarray,
    damping: float,
) -> float:
    target_active = np.asarray(target, dtype=float)[np.asarray(active, dtype=bool)]
    linear_active = np.asarray(linear, dtype=float)[np.asarray(active, dtype=bool)]
    quadratic_active = np.asarray(quadratic, dtype=float)[np.asarray(active, dtype=bool)]
    denominator = max(_vector_rms(target_active), np.finfo(float).tiny)
    change = damping * linear_active + damping * damping * quadratic_active
    return float(_vector_rms(target_active - change) / denominator)


def evaluate_fixed_damping_on_disjoint_cells(
    cells: Iterable[Mapping[str, Any]],
    *,
    damping: float,
    calibration_pairs: Iterable[tuple[float, float]],
) -> dict[str, Any]:
    """Evaluate one fixed damping without fitting anything on validation cells."""
    lam = float(damping)
    if not np.isfinite(lam) or not (0.0 <= lam <= 1.0):
        raise ValueError("damping must lie in [0,1]")

    frozen_calibration = {(float(t), float(z)) for t, z in calibration_pairs}
    rows = tuple(dict(row) for row in cells)
    if not rows:
        raise ValueError("at least one validation cell is required")

    parsed: list[dict[str, Any]] = []
    seen: set[tuple[float, float]] = set()
    proposed_l1: float | None = None
    target_parts: list[np.ndarray] = []
    linear_parts: list[np.ndarray] = []
    quadratic_parts: list[np.ndarray] = []

    for row in rows:
        time = float(row["time"])
        z = float(row["z"])
        pair = (time, z)
        if not np.isfinite([time, z]).all() or pair in seen:
            raise ValueError("validation time/z pairs must be finite and unique")
        if pair in frozen_calibration:
            raise ValueError("validation cells must be disjoint from damping calibration cells")
        seen.add(pair)

        target = np.asarray(row["target_stress"], dtype=float)
        active = np.asarray(row["active_target_mask"], dtype=bool)
        measurement = dict(row["measurement"])
        linear = np.asarray(measurement["linear_covariance_change_theta_axial"], dtype=float)
        quadratic = np.asarray(measurement["self_covariance_theta_axial"], dtype=float)
        exact = np.asarray(measurement["exact_covariance_change_theta_axial"], dtype=float)
        if target.ndim != 2 or target.shape[1] != 2 or not np.isfinite(target).all():
            raise ValueError("each target_stress must be finite with shape (N,2)")
        if active.shape != (len(target),) or not np.any(active):
            raise ValueError("each active_target_mask must select at least one node")
        if linear.shape != target.shape or quadratic.shape != target.shape or exact.shape != target.shape:
            raise ValueError("measurement arrays must match each validation target")
        if not np.isfinite(linear).all() or not np.isfinite(quadratic).all() or not np.isfinite(exact).all():
            raise ValueError("measurement arrays must be finite")
        if not bool(measurement.get("quadratic_identity_passed", False)):
            raise ValueError("every validation cell must pass quadratic identity closure")

        closure_scale = max(_vector_rms((linear + quadratic)[active]), np.finfo(float).tiny)
        closure = _vector_rms((exact - linear - quadratic)[active]) / closure_scale
        if closure > QUADRATIC_IDENTITY_RELATIVE_TOLERANCE:
            raise ValueError("a validation cell is inconsistent with linear+self covariance")

        cell_l1 = measurement.get("aggregate_l1_update")
        if cell_l1 is None:
            raise ValueError("every validation cell must declare aggregate_l1_update")
        cell_l1 = float(cell_l1)
        if not np.isfinite(cell_l1) or cell_l1 < 0.0:
            raise ValueError("aggregate_l1_update must be finite and nonnegative")
        if proposed_l1 is None:
            proposed_l1 = cell_l1
        elif not np.isclose(cell_l1, proposed_l1, rtol=0.0, atol=64.0 * np.finfo(float).eps):
            raise ValueError("the same frozen update must be used in every validation cell")

        fixed_relative = _fixed_relative_residual(target, active, linear, quadratic, lam)
        full_relative = _fixed_relative_residual(target, active, linear, quadratic, 1.0)
        parsed.append(
            {
                "time": time,
                "z": z,
                "active_target_nodes": int(np.count_nonzero(active)),
                "target_vector_rms": _vector_rms(target[active]),
                "fixed_damping_relative_stress_residual_rms": fixed_relative,
                "full_step_relative_stress_residual_rms": full_relative,
                "fixed_damping_improves_over_zero_update": bool(fixed_relative < 1.0),
                "fixed_damping_not_worse_than_full_step": bool(
                    fixed_relative <= full_relative + 64.0 * np.finfo(float).eps
                ),
                "quadratic_identity_rechecked_relative_rms": float(closure),
            }
        )
        target_parts.append(target[active])
        linear_parts.append(linear[active])
        quadratic_parts.append(quadratic[active])

    assert proposed_l1 is not None
    stacked_target = np.concatenate(target_parts, axis=0)
    stacked_linear = np.concatenate(linear_parts, axis=0)
    stacked_quadratic = np.concatenate(quadratic_parts, axis=0)
    stacked_active = np.ones(len(stacked_target), dtype=bool)
    aggregate_fixed = _fixed_relative_residual(
        stacked_target, stacked_active, stacked_linear, stacked_quadratic, lam
    )
    aggregate_full = _fixed_relative_residual(
        stacked_target, stacked_active, stacked_linear, stacked_quadratic, 1.0
    )
    improves_every = all(row["fixed_damping_improves_over_zero_update"] for row in parsed)
    not_worse_every = all(row["fixed_damping_not_worse_than_full_step"] for row in parsed)

    return {
        "validation_cell_count": len(parsed),
        "fixed_damping": lam,
        "proposed_aggregate_l1_update": proposed_l1,
        "fixed_damped_aggregate_l1_update": float(lam * proposed_l1),
        "aggregate_fixed_damping_relative_stress_residual_rms": aggregate_fixed,
        "aggregate_full_step_relative_stress_residual_rms": aggregate_full,
        "aggregate_zero_update_relative_stress_residual_rms": 1.0,
        "fixed_damping_improves_every_validation_cell": bool(improves_every),
        "fixed_damping_not_worse_than_full_step_every_validation_cell": bool(not_worse_every),
        "validation_transfer_preflight_passed": bool(improves_every and not_worse_every),
        "heldout_damping_reoptimized": False,
        "validation_objective_used_for_selection": False,
        "cells": parsed,
    }


def _build_cells(
    *,
    times: tuple[float, ...],
    z_values: tuple[float, ...],
    family: Any,
    coordinates: KokunoBoundedFamilyCoefficientCoordinates,
    delta_common: float,
    radial_count: int,
    angular_count: int,
    phase_count: int,
) -> tuple[list[dict[str, Any]], list[tuple[float, float, MissingCovarianceColumnTarget, dict[str, Any]]]]:
    rows: list[dict[str, Any]] = []
    receipts: list[tuple[float, float, MissingCovarianceColumnTarget, dict[str, Any]]] = []
    for z in z_values:
        for time in times:
            receipt, metadata = build_real_covariance_slice(
                time=time,
                z=z,
                radial_count=radial_count,
                angular_count=angular_count,
                phase_count=phase_count,
            )
            measurement = measure_bounded_coordinate_quadratic_change(
                family,
                receipt.radii,
                coordinate_contract=coordinates,
                delta_common=delta_common,
                delta_band=0.0,
                time=time,
                z=z,
                angular_count=angular_count,
                phase_count=phase_count,
            )
            rows.append(
                {
                    "time": time,
                    "z": z,
                    "target_stress": receipt.target_stress,
                    "active_target_mask": receipt.active_target_mask,
                    "measurement": measurement,
                }
            )
            receipts.append((time, z, receipt, metadata))
    return rows, receipts


def generate_actual_holdout_report(
    *,
    output: str | Path = "artifacts/kokuno_agent3/spacetime_damping_holdout_report.json",
    validation_times: Iterable[float] = VALIDATION_TIMES,
    validation_z: Iterable[float] = VALIDATION_Z,
    radial_count: int = DEFAULT_RADIAL_COUNT,
    angular_count: int = DEFAULT_ANGULAR_COUNT,
    phase_count: int = DEFAULT_PHASE_COUNT,
) -> dict[str, Any]:
    """Calibrate on #422 cells, then transfer one frozen step to disjoint cells."""
    calibration_times = tuple(float(value) for value in HELD_OUT_CYCLE_TIMES)
    calibration_z = tuple(float(value) for value in SPATIAL_SCREEN_Z)
    fresh_times = _validate_axis(validation_times, name="validation_times")
    fresh_z = _validate_axis(validation_z, name="validation_z")
    calibration_pairs = {(t, z) for z in calibration_z for t in calibration_times}
    validation_pairs = {(t, z) for z in fresh_z for t in fresh_times}
    if calibration_pairs & validation_pairs:
        raise ValueError("validation spacetime grid must be disjoint from #422 calibration grid")

    reference_receipt, reference_metadata = build_real_covariance_slice(
        time=REFERENCE_BUDGET_TIME,
        z=REFERENCE_BUDGET_Z,
        radial_count=int(radial_count),
        angular_count=int(angular_count),
        phase_count=int(phase_count),
    )
    frozen_amplitude_budget = current_signed_coefficient_budget(reference_receipt)
    fractional_budget = convert_amplitude_budget_to_fractional(
        frozen_amplitude_budget, ROUTED_OSCILLATORY_AMPLITUDE
    )
    coordinates = KokunoBoundedFamilyCoefficientCoordinates(max_l1_update=fractional_budget)
    family = _one_label_real_family()

    reference_linear = measure_bounded_coordinate_covariance(
        family,
        reference_receipt.radii,
        coordinate_contract=coordinates,
        time=REFERENCE_BUDGET_TIME,
        z=REFERENCE_BUDGET_Z,
        angular_count=int(angular_count),
        phase_count=int(phase_count),
    )
    reference_jacobian = np.asarray(
        reference_linear["coordinate_covariance_jacobian_theta_axial"], dtype=float
    )
    reference_inverse = evaluate_family_bounded_inverse(
        reference_receipt,
        reference_jacobian,
        coefficient_labels=EXPECTED_PARAMETER_NAMES,
        current_coefficient_budget=frozen_amplitude_budget,
        additional_family_l1_budget=fractional_budget,
        budget_unit_matches_jacobian=True,
        coefficient_unit=EXPECTED_PARAMETER_UNIT,
    )

    delta_common = FROZEN_COMMON_UPDATE_SIGN * fractional_budget
    calibration_rows, calibration_receipts = _build_cells(
        times=calibration_times,
        z_values=calibration_z,
        family=family,
        coordinates=coordinates,
        delta_common=delta_common,
        radial_count=int(radial_count),
        angular_count=int(angular_count),
        phase_count=int(phase_count),
    )
    calibration = solve_shared_spacetime_damping(calibration_rows)
    frozen_damping = float(calibration["shared_selected_damping"])

    validation_rows, validation_receipts = _build_cells(
        times=fresh_times,
        z_values=fresh_z,
        family=family,
        coordinates=coordinates,
        delta_common=delta_common,
        radial_count=int(radial_count),
        angular_count=int(angular_count),
        phase_count=int(phase_count),
    )
    validation = evaluate_fixed_damping_on_disjoint_cells(
        validation_rows,
        damping=frozen_damping,
        calibration_pairs=calibration_pairs,
    )

    leading_sha = str(reference_metadata["leading_candidate_sha256"])
    for _, _, _, metadata in calibration_receipts + validation_receipts:
        if str(metadata["leading_candidate_sha256"]) != leading_sha:
            raise RuntimeError("leading candidate changed across calibration/validation reconstruction")

    rank_passed = bool(reference_inverse["family_bounded_inverse_preflight_passed"])
    transfer_passed = bool(validation["validation_transfer_preflight_passed"])
    report: dict[str, Any] = {
        "task": TASK,
        "source": {
            "repository": SOURCE_REPOSITORY,
            "commit": SOURCE_COMMIT,
            "path": SOURCE_PATH,
            "edition": SOURCE_EDITION,
            "doi": SOURCE_DOI,
            "section": SOURCE_SECTION,
            "source_scope": (
                "The corrected reader retains nonlinear self-covariance after the signed covariance inverse. "
                "The disjoint calibration/validation split is repository validation engineering."
            ),
        },
        "handoff": {
            "agent3_calibration_pr": AGENT3_CALIBRATION_PR,
            "agent3_bounded_inverse_pr": AGENT3_BOUNDED_INVERSE_PR,
            "agent3_radial_side_effect_pr": AGENT3_RADIAL_SIDE_EFFECT_PR,
            "agent2_latest_dependency": (
                "PR #430 executes the source band-covering schedule but still lacks the actual "
                "positive-order/background multi-band physical velocity family and a genuinely "
                "independent second covariance direction"
            ),
        },
        "inputs": {
            "leading_candidate_sha256": leading_sha,
            "oscillatory_amplitude": ROUTED_OSCILLATORY_AMPLITUDE,
            "oscillatory_phase": ROUTED_OSCILLATORY_PHASE,
            "profile_annulus": [float(PROFILE_ANNULUS[0]), float(PROFILE_ANNULUS[1])],
            "radial_count": int(radial_count),
            "angular_count": int(angular_count),
            "phase_count": int(phase_count),
            "surrogate_defect_used": False,
            "frozen_physical_amplitude_budget": float(frozen_amplitude_budget),
            "converted_fractional_budget": float(fractional_budget),
            "budget_changed": False,
            "algebraic_relative_residual_tolerance_reused_from_pr393": (
                ALGEBRAIC_RELATIVE_RESIDUAL_TOLERANCE
            ),
            "quadratic_identity_relative_tolerance_reused": QUADRATIC_IDENTITY_RELATIVE_TOLERANCE,
            "frozen_delta_common": float(delta_common),
            "frozen_delta_band": 0.0,
        },
        "damping_calibration": {
            "role": "calibration_only_for_damping_selection",
            "times": list(calibration_times),
            "z_values": list(calibration_z),
            "cell_count": len(calibration_rows),
            "shared_result": calibration,
        },
        "disjoint_validation": {
            "role": "validation_only_no_damping_selection",
            "times": list(fresh_times),
            "z_values": list(fresh_z),
            "cell_count": len(validation_rows),
            "result": validation,
        },
        "reference_one_band_bounded_inverse": reference_inverse,
        "routing": {
            "damping_validation_transfer_passed": transfer_passed,
            "genuinely_independent_second_covariance_column_ready": False,
            "upstream_bounded_inverse_passed": rank_passed,
            "public_velocity_correction_materialized": False,
            "radial_force_retained_in_materialized_correction": False,
            "finite_correction_cycle_rerun_allowed": False,
            "next_required": (
                "Agent 2 must expose an actual source-motivated multi-band physical velocity family with "
                "an independent covariance direction. Only after unit/rank/budget/spacetime/radial guards "
                "pass may the same disjoint damping-transfer protocol be applied to the actual correction "
                "before materialization and a fresh held-in/held-out finite NS correction cycle."
            ),
        },
        "truth_boundary": {
            "real_candidate_defect_targets_consumed": True,
            "surrogate_defect_used": False,
            "damping_calibration_validation_disjoint": True,
            "heldout_damping_reoptimized": False,
            "validation_objective_used_for_selection": False,
            "actual_source_multiband_family_consumed": False,
            "genuinely_independent_second_covariance_column_ready": False,
            "public_velocity_correction_materialized": False,
            "radial_force_retained_in_materialized_correction": False,
            "finite_correction_cycle_run": False,
            "residual_reduction_claimed": False,
            "formal_full_domain_normalized_pde_gate_assessed": False,
            "pde_validated": False,
            "paper_exact": False,
            "openai_field_identified": False,
            "blowup_proved": False,
        },
    }
    destination = Path(output)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    return report


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        default="artifacts/kokuno_agent3/spacetime_damping_holdout_report.json",
    )
    parser.add_argument("--radial-count", type=int, default=DEFAULT_RADIAL_COUNT)
    parser.add_argument("--angular-count", type=int, default=DEFAULT_ANGULAR_COUNT)
    parser.add_argument("--phase-count", type=int, default=DEFAULT_PHASE_COUNT)
    return parser


def main() -> None:
    args = _parser().parse_args()
    report = generate_actual_holdout_report(
        output=args.output,
        radial_count=args.radial_count,
        angular_count=args.angular_count,
        phase_count=args.phase_count,
    )
    result = report["disjoint_validation"]["result"]
    calibration = report["damping_calibration"]["shared_result"]
    summary = {
        "task": report["task"],
        "calibration_selected_damping": calibration["shared_selected_damping"],
        "calibration_relative_stress_residual_rms": calibration[
            "shared_best_relative_stress_residual_rms"
        ],
        "validation_fixed_damping_relative_stress_residual_rms": result[
            "aggregate_fixed_damping_relative_stress_residual_rms"
        ],
        "validation_full_step_relative_stress_residual_rms": result[
            "aggregate_full_step_relative_stress_residual_rms"
        ],
        "validation_transfer_preflight_passed": result["validation_transfer_preflight_passed"],
        "heldout_damping_reoptimized": result["heldout_damping_reoptimized"],
        "upstream_bounded_inverse_passed": report["routing"]["upstream_bounded_inverse_passed"],
        "finite_correction_cycle_rerun_allowed": report["routing"][
            "finite_correction_cycle_rerun_allowed"
        ],
    }
    print(json.dumps(summary, sort_keys=True, allow_nan=False))


if __name__ == "__main__":
    main()
