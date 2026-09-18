"""Require frozen damping to survive the source-required radial stress side effect.

Agent-3 PR #431 validates one frozen scalar damping on spacetime cells disjoint
from the cells used to select that damping, but its local objective contains only
the theta/z covariance-stress channels.  Agent-3 PRs #357/#370 separately show
that the corrected Kokuno stress has the additional source-required divergence

    (div T)_r = d_z sigma_1.

For one already-proposed physical update V, PRs #411/#415 retain the exact
quadratic covariance change

    Delta C(lambda) = lambda B(W,V) + lambda^2 C(V).

The axial component of that covariance change is the proposed sigma_1 change.
This module therefore differentiates the *same frozen damped update* in z on the
unchanged (.02,.01,.005) ladder and evaluates its signed radial divergence
against the real ring-averaged radial phase-mean defect on the disjoint #431
validation cells.

The routing rule is deliberately conservative and threshold-free: a step that
passes theta/z transfer is still rejected if its source-sign radial divergence
increases the gated radial-defect RMS on any fresh cell.  This is repository
finite-step engineering, not a new Kokuno formula and not a Navier--Stokes PDE
acceptance test.
"""
from __future__ import annotations

import argparse
from dataclasses import replace
import json
from pathlib import Path
from typing import Any, Iterable, Mapping

import numpy as np

from .kokuno_bounded_coordinate_covariance import _one_label_real_family
from .kokuno_bounded_family_coefficients import KokunoBoundedFamilyCoefficientCoordinates
from .kokuno_complete_curl import KokunoCompleteCurlCorrection
from .kokuno_leading_core_series import KokunoLeadingCoreSeriesCandidate
from .kokuno_mean_defect import PhaseMeanDefectContract
from .kokuno_missing_covariance_column_target import (
    DEFAULT_ANGULAR_COUNT,
    DEFAULT_PHASE_COUNT,
    DEFAULT_RADIAL_COUNT,
)
from .kokuno_multislice_covariance_preflight import build_real_covariance_slice
from .kokuno_quadratic_covariance_gain_guard import (
    measure_bounded_coordinate_quadratic_change,
)
from .kokuno_radial_force_completeness_audit import (
    Z_DERIVATIVE_STEP_LADDER,
    audit_radial_force_derivative_ladder,
    centered_radial_force,
    sample_real_radial_mean_defect,
)
from .kokuno_signed_covariance_inverse import (
    ROUTED_OSCILLATORY_AMPLITUDE,
    ROUTED_OSCILLATORY_PHASE,
)
from .kokuno_spacetime_damping_holdout import (
    VALIDATION_TIMES,
    VALIDATION_Z,
    generate_actual_holdout_report,
)

TASK = "KOKUNO-A3-RADIAL-AWARE-DAMPING-HOLDOUT-025"
SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
SOURCE_EDITION = "2026.09.09-consolidated"
SOURCE_DOI = "10.5281/zenodo.22678406"
SOURCE_FORMULAS = "R33-R34 and R41; signed covariance inverse with retained self-covariance"
AGENT3_DAMPING_HOLDOUT_PR = 431
AGENT3_RADIAL_FORCE_PR = 370
AGENT3_QUADRATIC_DAMPING_PR = 415


def _rms(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    if values.ndim != 1 or not np.isfinite(values).all():
        raise ValueError("values must be a finite one-dimensional array")
    return float(np.sqrt(np.mean(values * values)))


def damped_axial_covariance_stress(
    measurement: Mapping[str, Any],
    damping: float,
) -> np.ndarray:
    """Return axial ``lambda*B + lambda^2*C`` from one measured full update."""
    lam = float(damping)
    if not np.isfinite(lam) or not 0.0 <= lam <= 1.0:
        raise ValueError("damping must lie in [0,1]")
    linear = np.asarray(measurement["linear_covariance_change_theta_axial"], dtype=float)
    quadratic = np.asarray(measurement["self_covariance_theta_axial"], dtype=float)
    exact = np.asarray(measurement["exact_covariance_change_theta_axial"], dtype=float)
    if (
        linear.ndim != 2
        or linear.shape[1] != 2
        or quadratic.shape != linear.shape
        or exact.shape != linear.shape
    ):
        raise ValueError("covariance arrays must have common shape (N,2)")
    if not np.isfinite(linear).all() or not np.isfinite(quadratic).all() or not np.isfinite(exact).all():
        raise ValueError("covariance arrays must be finite")
    if not bool(measurement.get("quadratic_identity_passed", False)):
        raise ValueError("measurement must pass the retained quadratic covariance identity")
    return lam * linear[:, 1] + lam * lam * quadratic[:, 1]


def evaluate_fixed_damping_radial_cell(
    gated_radial_defect: np.ndarray,
    radial_force_derivatives: Iterable[np.ndarray],
    *,
    z_steps: Iterable[float] = Z_DERIVATIVE_STEP_LADDER,
    active_mask: np.ndarray | None = None,
) -> dict[str, Any]:
    """Evaluate the source-sign radial divergence of one frozen damped step.

    The compact stress contributes ``+ d_z sigma_1`` to the radial divergence,
    so the signed radial residual proxy is ``gated_defect + d_z sigma_1``.
    No sign is selected from validation data.
    """
    defect = np.asarray(gated_radial_defect, dtype=float)
    derivatives = tuple(np.asarray(value, dtype=float) for value in radial_force_derivatives)
    if defect.ndim != 1 or len(defect) < 9 or not np.isfinite(defect).all():
        raise ValueError("gated_radial_defect must be finite 1-D with at least 9 nodes")
    if not derivatives or any(value.shape != defect.shape for value in derivatives):
        raise ValueError("every radial-force derivative must match the defect grid")

    if active_mask is None:
        mask = np.ones(defect.shape, dtype=bool)
    else:
        mask = np.asarray(active_mask, dtype=bool)
        if mask.shape != defect.shape or not np.any(mask):
            raise ValueError("active_mask must select at least one radial node")

    derivative_audit = audit_radial_force_derivative_ladder(
        z_steps,
        derivatives,
        active_mask=mask,
    )
    finest = derivatives[-1]
    before = _rms(defect[mask])
    after_values = defect[mask] + finest[mask]
    after = _rms(after_values)
    scale = max(before, np.finfo(float).tiny)
    numerical_slack = 64.0 * np.finfo(float).eps * max(before, after, 1.0)
    nonworse = bool(after <= before + numerical_slack)
    improves = bool(after < before - numerical_slack)
    return {
        "source_radial_divergence_sign": "+d_z_sigma_1",
        "active_radial_nodes": int(np.count_nonzero(mask)),
        "gated_radial_defect_rms_before": before,
        "damped_radial_force_rms": _rms(finest[mask]),
        "gated_radial_defect_rms_after_source_sign": after,
        "radial_relative_residual_rms_after_source_sign": float(after / scale),
        "radial_source_sign_improves": improves,
        "radial_source_sign_not_worse": nonworse,
        "radial_force_derivative_audit": derivative_audit,
    }


def evaluate_radial_aware_holdout_transfer(
    tangential_validation: Mapping[str, Any],
    radial_cells: Iterable[Mapping[str, Any]],
) -> dict[str, Any]:
    """Combine #431 theta/z transfer with radial-channel fail-closed routing."""
    validation = dict(tangential_validation)
    cells = tuple(dict(row) for row in radial_cells)
    if not cells:
        raise ValueError("at least one radial validation cell is required")
    if validation.get("heldout_damping_reoptimized") is not False:
        raise ValueError("held-out damping must remain frozen")
    if validation.get("validation_objective_used_for_selection") is not False:
        raise ValueError("validation objective must not select damping")

    tangential_pairs = {
        (float(row["time"]), float(row["z"])) for row in validation.get("cells", ())
    }
    radial_pairs: set[tuple[float, float]] = set()
    for row in cells:
        pair = (float(row["time"]), float(row["z"]))
        if pair in radial_pairs:
            raise ValueError("radial validation cells must have unique time/z pairs")
        radial_pairs.add(pair)
        result = row.get("radial_result", {})
        if result.get("radial_force_derivative_audit", {}).get(
            "radial_force_derivative_stability_preflight_passed"
        ) is not True:
            raise ValueError("every radial validation cell must pass derivative stability")
    if tangential_pairs != radial_pairs:
        raise ValueError("radial and tangential validation cells must match exactly")

    tangential_passed = bool(validation.get("validation_transfer_preflight_passed", False))
    all_radial_nonworse = all(
        bool(row["radial_result"]["radial_source_sign_not_worse"]) for row in cells
    )
    all_radial_improve = all(
        bool(row["radial_result"]["radial_source_sign_improves"]) for row in cells
    )
    ratios = np.asarray(
        [
            float(row["radial_result"]["radial_relative_residual_rms_after_source_sign"])
            for row in cells
        ],
        dtype=float,
    )
    if not np.isfinite(ratios).all() or np.any(ratios < 0.0):
        raise ValueError("radial relative residuals must be finite and nonnegative")

    return {
        "validation_cell_count": len(cells),
        "tangential_validation_transfer_passed": tangential_passed,
        "radial_source_sign_not_worse_every_validation_cell": bool(all_radial_nonworse),
        "radial_source_sign_improves_every_validation_cell": bool(all_radial_improve),
        "radial_relative_residual_rms_min": float(np.min(ratios)),
        "radial_relative_residual_rms_max": float(np.max(ratios)),
        "radial_aware_damping_transfer_preflight_passed": bool(
            tangential_passed and all_radial_nonworse
        ),
        "heldout_damping_reoptimized": False,
        "validation_objective_used_for_selection": False,
        "cells": list(cells),
        "interpretation": (
            "This guard adds only the source-required signed radial divergence check to "
            "the already-disjoint compact-stress transfer. It is not a full Navier-Stokes "
            "residual and cannot authorize a finite correction cycle by itself."
        ),
    }


def _validate_axis(values: Iterable[float], *, name: str) -> tuple[float, ...]:
    result = tuple(float(value) for value in values)
    if not result or not np.isfinite(result).all():
        raise ValueError(f"{name} must contain finite values")
    if any(right <= left for left, right in zip(result, result[1:])):
        raise ValueError(f"{name} must be strictly increasing")
    return result


def generate_actual_core_report(
    *,
    output: str | Path = "artifacts/kokuno_agent3/radial_aware_damping_holdout_report.json",
    parent_output: str | Path = "artifacts/kokuno_agent3/radial_aware_damping_parent_431.json",
    validation_times: Iterable[float] = VALIDATION_TIMES,
    validation_z: Iterable[float] = VALIDATION_Z,
    radial_count: int = DEFAULT_RADIAL_COUNT,
    angular_count: int = DEFAULT_ANGULAR_COUNT,
    phase_count: int = DEFAULT_PHASE_COUNT,
) -> dict[str, Any]:
    """Replay #431, then evaluate the same frozen step's radial divergence."""
    fresh_times = _validate_axis(validation_times, name="validation_times")
    fresh_z = _validate_axis(validation_z, name="validation_z")
    parent = generate_actual_holdout_report(
        output=parent_output,
        validation_times=fresh_times,
        validation_z=fresh_z,
        radial_count=int(radial_count),
        angular_count=int(angular_count),
        phase_count=int(phase_count),
    )
    parent_validation = parent["disjoint_validation"]["result"]
    frozen_damping = float(parent["damping_calibration"]["shared_result"]["shared_selected_damping"])
    fractional_budget = float(parent["inputs"]["converted_fractional_budget"])
    delta_common = float(parent["inputs"]["frozen_delta_common"])

    coordinates = KokunoBoundedFamilyCoefficientCoordinates(max_l1_update=fractional_budget)
    family = _one_label_real_family()
    leading = KokunoLeadingCoreSeriesCandidate()
    oscillatory = replace(
        KokunoCompleteCurlCorrection(),
        amplitude=ROUTED_OSCILLATORY_AMPLITUDE,
        phase=ROUTED_OSCILLATORY_PHASE,
    )
    contract = PhaseMeanDefectContract(phase_count=int(phase_count))

    radial_rows: list[dict[str, Any]] = []
    for z in fresh_z:
        for time in fresh_times:
            receipt, metadata = build_real_covariance_slice(
                time=time,
                z=z,
                radial_count=int(radial_count),
                angular_count=int(angular_count),
                phase_count=int(phase_count),
            )
            radii = np.asarray(receipt.radii, dtype=float)
            active = np.ones(len(radii), dtype=bool)
            active[[0, -1]] = False
            _, gated_radial = sample_real_radial_mean_defect(
                leading.at_points,
                oscillatory,
                contract,
                radii,
                time=time,
                z=z,
                angular_count=int(angular_count),
            )

            derivatives: list[np.ndarray] = []
            level_receipts: list[dict[str, Any]] = []
            for step in Z_DERIVATIVE_STEP_LADDER:
                plus = measure_bounded_coordinate_quadratic_change(
                    family,
                    radii,
                    coordinate_contract=coordinates,
                    delta_common=delta_common,
                    delta_band=0.0,
                    time=time,
                    z=z + float(step),
                    angular_count=int(angular_count),
                    phase_count=int(phase_count),
                )
                minus = measure_bounded_coordinate_quadratic_change(
                    family,
                    radii,
                    coordinate_contract=coordinates,
                    delta_common=delta_common,
                    delta_band=0.0,
                    time=time,
                    z=z - float(step),
                    angular_count=int(angular_count),
                    phase_count=int(phase_count),
                )
                sigma_plus = damped_axial_covariance_stress(plus, frozen_damping)
                sigma_minus = damped_axial_covariance_stress(minus, frozen_damping)
                derivative = centered_radial_force(sigma_plus, sigma_minus, float(step))
                derivatives.append(derivative)
                level_receipts.append(
                    {
                        "z_step": float(step),
                        "z_minus": float(z - step),
                        "z_plus": float(z + step),
                        "plus_quadratic_identity_passed": bool(plus["quadratic_identity_passed"]),
                        "minus_quadratic_identity_passed": bool(minus["quadratic_identity_passed"]),
                    }
                )

            radial_result = evaluate_fixed_damping_radial_cell(
                gated_radial,
                derivatives,
                z_steps=Z_DERIVATIVE_STEP_LADDER,
                active_mask=active,
            )
            radial_rows.append(
                {
                    "time": float(time),
                    "z": float(z),
                    "leading_candidate_sha256": str(metadata["leading_candidate_sha256"]),
                    "surrogate_defect_used": False,
                    "z_derivative_levels": level_receipts,
                    "radial_result": radial_result,
                }
            )

    radial_transfer = evaluate_radial_aware_holdout_transfer(
        parent_validation,
        radial_rows,
    )
    parent_rank_passed = bool(parent["routing"]["upstream_bounded_inverse_passed"])
    report: dict[str, Any] = {
        "task": TASK,
        "source": {
            "repository": SOURCE_REPOSITORY,
            "commit": SOURCE_COMMIT,
            "path": SOURCE_PATH,
            "edition": SOURCE_EDITION,
            "doi": SOURCE_DOI,
            "formulas": SOURCE_FORMULAS,
            "scope": (
                "R33-R34/R41 require +d_z sigma_1 in the radial stress divergence. "
                "The fixed-damping quadratic covariance expansion is repository algebra "
                "already introduced in Agent-3 #415."
            ),
        },
        "handoff": {
            "agent3_damping_holdout_pr": AGENT3_DAMPING_HOLDOUT_PR,
            "agent3_radial_force_pr": AGENT3_RADIAL_FORCE_PR,
            "agent3_quadratic_damping_pr": AGENT3_QUADRATIC_DAMPING_PR,
            "agent2_latest_dependency": (
                "PR #439 certifies source-compatible rational rectangle separation but still "
                "does not materialize the actual positive-order/background multi-band physical "
                "mode family or public xyz,t oscillatory correction."
            ),
        },
        "inputs": {
            "validation_times": list(fresh_times),
            "validation_z": list(fresh_z),
            "z_derivative_steps": list(Z_DERIVATIVE_STEP_LADDER),
            "radial_count": int(radial_count),
            "angular_count": int(angular_count),
            "phase_count": int(phase_count),
            "frozen_damping_from_pr431": frozen_damping,
            "frozen_delta_common": delta_common,
            "frozen_delta_band": 0.0,
            "converted_fractional_budget": fractional_budget,
            "budget_changed": False,
            "heldout_damping_reoptimized": False,
            "validation_objective_used_for_selection": False,
            "surrogate_defect_used": False,
            "radial_active_mask": "all interior annulus nodes; endpoints excluded",
        },
        "parent_disjoint_tangential_transfer": parent_validation,
        "radial_aware_disjoint_transfer": radial_transfer,
        "routing": {
            "radial_aware_damping_transfer_preflight_passed": radial_transfer[
                "radial_aware_damping_transfer_preflight_passed"
            ],
            "genuinely_independent_second_covariance_column_ready": False,
            "upstream_bounded_inverse_passed": parent_rank_passed,
            "public_velocity_correction_materialized": False,
            "radial_force_retained_in_materialized_correction": False,
            "finite_correction_cycle_rerun_allowed": False,
            "next_required": (
                "Agent 2 must still supply the actual source-motivated multi-band physical "
                "velocity family with a genuine independent covariance direction. Only after "
                "the existing unit/rank/budget/spacetime/quadratic/radial guards pass may a "
                "public delta_u be materialized and a fresh held-in/held-out NS cycle run."
            ),
        },
        "truth_boundary": {
            "real_candidate_defect_targets_consumed": True,
            "surrogate_defect_used": False,
            "damping_calibration_validation_disjoint": True,
            "heldout_damping_reoptimized": False,
            "source_radial_force_sign_used": True,
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
    destination.write_text(
        json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    return report


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        default="artifacts/kokuno_agent3/radial_aware_damping_holdout_report.json",
    )
    parser.add_argument(
        "--parent-output",
        default="artifacts/kokuno_agent3/radial_aware_damping_parent_431.json",
    )
    parser.add_argument("--radial-count", type=int, default=DEFAULT_RADIAL_COUNT)
    parser.add_argument("--angular-count", type=int, default=DEFAULT_ANGULAR_COUNT)
    parser.add_argument("--phase-count", type=int, default=DEFAULT_PHASE_COUNT)
    return parser


def main() -> None:
    args = _parser().parse_args()
    report = generate_actual_core_report(
        output=args.output,
        parent_output=args.parent_output,
        radial_count=args.radial_count,
        angular_count=args.angular_count,
        phase_count=args.phase_count,
    )
    result = report["radial_aware_disjoint_transfer"]
    print(
        json.dumps(
            {
                "task": report["task"],
                "frozen_damping": report["inputs"]["frozen_damping_from_pr431"],
                "tangential_validation_transfer_passed": result[
                    "tangential_validation_transfer_passed"
                ],
                "radial_source_sign_not_worse_every_validation_cell": result[
                    "radial_source_sign_not_worse_every_validation_cell"
                ],
                "radial_relative_residual_rms_min": result[
                    "radial_relative_residual_rms_min"
                ],
                "radial_relative_residual_rms_max": result[
                    "radial_relative_residual_rms_max"
                ],
                "radial_aware_damping_transfer_preflight_passed": result[
                    "radial_aware_damping_transfer_preflight_passed"
                ],
                "finite_correction_cycle_rerun_allowed": report["routing"][
                    "finite_correction_cycle_rerun_allowed"
                ],
            },
            sort_keys=True,
            allow_nan=False,
        )
    )


if __name__ == "__main__":
    main()
