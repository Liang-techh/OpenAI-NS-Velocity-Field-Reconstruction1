"""Screen one shared damped Kokuno covariance step across frozen spacetime cells.

Agent 3 PR #415 chooses an exact scalar damping for one already-proposed bounded
physical update while retaining the nonlinear correction self-covariance.  That
receipt is local to the reference cell ``(t=.5,z=.08)``.  A finite correction
coefficient, however, is not allowed to be retuned independently at every
validation slice.

This module therefore reuses the existing frozen Agent-3 spacetime screen

    t = (.375, .5, .625),   z = (.06, .08, .10)

and chooses one *shared* ``lambda in [0,1]`` for the same predeclared physical
update at all nine cells.  The stacked squared mismatch is still a quartic in
``lambda`` because every cell has

    Delta C(lambda) = lambda B(W,V) + lambda^2 C(V).

The source provenance is deliberately narrower than this engineering rule: the
corrected Kokuno reader retains nonlinear self-covariance after the signed
covariance inverse.  Shared spacetime damping, the frozen 3x3 screen, and the
aggregate objective are repository engineering.  They do not alter any PDE or
algebraic acceptance threshold and they do not create Agent 2's still-missing
independent covariance direction.
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
from .kokuno_damped_quadratic_covariance_gain import (
    ROOT_INTERVAL_TOLERANCE,
    solve_damped_quadratic_covariance_gain,
)
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
from .kokuno_spacetime_covariance_preflight import (
    REFERENCE_BUDGET_Z,
    SPATIAL_SCREEN_Z,
)

TASK = "KOKUNO-A3-SPACETIME-DAMPED-QUADRATIC-GAIN-023"
SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
SOURCE_EDITION = "2026.09.09-consolidated"
SOURCE_DOI = "10.5281/zenodo.22678406"
SOURCE_SECTION = "signed covariance inverse; nonlinear self-covariance retained"
AGENT3_DAMPED_GAIN_PR = 415
AGENT3_SPACETIME_PREFLIGHT_PR = 338
AGENT3_RADIAL_SIDE_EFFECT_PR = 370


def _vector_rms(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    if values.ndim != 2 or values.shape[1] != 2 or not np.isfinite(values).all():
        raise ValueError("values must be finite with shape (N,2)")
    return float(np.sqrt(np.mean(np.sum(values * values, axis=1))))


def _validate_axis(values: Iterable[float], *, name: str) -> tuple[float, ...]:
    result = tuple(float(value) for value in values)
    if len(result) < 2 or not np.isfinite(result).all():
        raise ValueError(f"{name} must contain at least two finite values")
    if any(right <= left for left, right in zip(result, result[1:])):
        raise ValueError(f"{name} must be strictly increasing")
    return result


def _shared_relative_residual(
    target: np.ndarray,
    active: np.ndarray,
    linear: np.ndarray,
    quadratic: np.ndarray,
    damping: float,
) -> float:
    target = np.asarray(target, dtype=float)
    active = np.asarray(active, dtype=bool)
    linear = np.asarray(linear, dtype=float)
    quadratic = np.asarray(quadratic, dtype=float)
    if target.ndim != 2 or target.shape[1] != 2:
        raise ValueError("target must have shape (N,2)")
    if active.shape != (len(target),) or not np.any(active):
        raise ValueError("active mask must select at least one node")
    if linear.shape != target.shape or quadratic.shape != target.shape:
        raise ValueError("linear/quadratic arrays must match target")
    lam = float(damping)
    if not np.isfinite(lam) or lam < 0.0 or lam > 1.0:
        raise ValueError("damping must lie in [0,1]")
    target_active = target[active]
    change = lam * linear[active] + lam * lam * quadratic[active]
    denominator = max(_vector_rms(target_active), np.finfo(float).tiny)
    return float(_vector_rms(target_active - change) / denominator)


def solve_shared_spacetime_damping(
    cells: Iterable[Mapping[str, Any]],
) -> dict[str, Any]:
    """Choose one exact quartic damping parameter for every supplied cell.

    Each row must provide ``time``, ``z``, ``target_stress``,
    ``active_target_mask`` and a PR-#411-style full-step ``measurement``.  The
    proposed bounded step size must be identical in every cell; otherwise the
    caller has already retuned the correction and this function fails closed.
    """
    rows = tuple(dict(row) for row in cells)
    if len(rows) < 2:
        raise ValueError("at least two spacetime cells are required")

    targets: list[np.ndarray] = []
    active_masks: list[np.ndarray] = []
    linear_parts: list[np.ndarray] = []
    quadratic_parts: list[np.ndarray] = []
    exact_parts: list[np.ndarray] = []
    parsed: list[dict[str, Any]] = []
    pairs: set[tuple[float, float]] = set()
    proposed_l1: float | None = None

    for row in rows:
        time = float(row["time"])
        z = float(row["z"])
        if not np.isfinite([time, z]).all() or (time, z) in pairs:
            raise ValueError("cell time/z pairs must be finite and unique")
        pairs.add((time, z))
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
            raise ValueError("measurement arrays must match each cell target")
        if not np.isfinite(linear).all() or not np.isfinite(quadratic).all() or not np.isfinite(exact).all():
            raise ValueError("measurement arrays must be finite")
        if not bool(measurement.get("quadratic_identity_passed", False)):
            raise ValueError("every cell must pass quadratic identity closure")
        closure_scale = max(_vector_rms((linear + quadratic)[active]), np.finfo(float).tiny)
        closure = _vector_rms((exact - linear - quadratic)[active]) / closure_scale
        if closure > QUADRATIC_IDENTITY_RELATIVE_TOLERANCE:
            raise ValueError("a cell exact change is inconsistent with linear+self covariance")

        cell_l1 = measurement.get("aggregate_l1_update")
        if cell_l1 is None:
            raise ValueError("every cell must declare aggregate_l1_update")
        cell_l1 = float(cell_l1)
        if not np.isfinite(cell_l1) or cell_l1 < 0.0:
            raise ValueError("aggregate_l1_update must be finite and nonnegative")
        if proposed_l1 is None:
            proposed_l1 = cell_l1
        elif not np.isclose(cell_l1, proposed_l1, rtol=0.0, atol=64.0 * np.finfo(float).eps):
            raise ValueError("the same proposed bounded update must be used in every cell")

        targets.append(target)
        active_masks.append(active)
        linear_parts.append(linear)
        quadratic_parts.append(quadratic)
        exact_parts.append(exact)
        parsed.append(
            {
                "time": time,
                "z": z,
                "target": target,
                "active": active,
                "linear": linear,
                "quadratic": quadratic,
                "measurement": measurement,
                "quadratic_identity_rechecked_relative_rms": float(closure),
            }
        )

    assert proposed_l1 is not None
    stacked_measurement = {
        "linear_covariance_change_theta_axial": np.concatenate(linear_parts, axis=0),
        "self_covariance_theta_axial": np.concatenate(quadratic_parts, axis=0),
        "exact_covariance_change_theta_axial": np.concatenate(exact_parts, axis=0),
        "quadratic_identity_passed": True,
        "aggregate_l1_update": proposed_l1,
    }
    shared = solve_damped_quadratic_covariance_gain(
        np.concatenate(targets, axis=0),
        np.concatenate(active_masks, axis=0),
        stacked_measurement,
    )
    shared_lambda = float(shared["selected_damping"])

    cell_reports: list[dict[str, Any]] = []
    for row in parsed:
        local = solve_damped_quadratic_covariance_gain(
            row["target"], row["active"], row["measurement"]
        )
        shared_relative = _shared_relative_residual(
            row["target"], row["active"], row["linear"], row["quadratic"], shared_lambda
        )
        cell_reports.append(
            {
                "time": row["time"],
                "z": row["z"],
                "active_target_nodes": int(np.count_nonzero(row["active"])),
                "target_vector_rms": _vector_rms(row["target"][row["active"]]),
                "individually_optimal_damping": float(local["selected_damping"]),
                "individually_optimal_relative_stress_residual_rms": float(
                    local["best_damped_relative_stress_residual_rms"]
                ),
                "shared_damping_relative_stress_residual_rms": shared_relative,
                "full_step_relative_stress_residual_rms": float(
                    local["full_step_relative_stress_residual_rms"]
                ),
                "shared_damping_improves_over_zero_update": bool(shared_relative < 1.0),
                "shared_damping_not_worse_than_full_step": bool(
                    shared_relative
                    <= float(local["full_step_relative_stress_residual_rms"])
                    + 64.0 * np.finfo(float).eps
                ),
                "quadratic_identity_rechecked_relative_rms": row[
                    "quadratic_identity_rechecked_relative_rms"
                ],
            }
        )

    return {
        "cell_count": len(cell_reports),
        "shared_selected_damping": shared_lambda,
        "shared_best_relative_stress_residual_rms": float(
            shared["best_damped_relative_stress_residual_rms"]
        ),
        "shared_full_step_relative_stress_residual_rms": float(
            shared["full_step_relative_stress_residual_rms"]
        ),
        "shared_zero_update_relative_stress_residual_rms": 1.0,
        "shared_damping_improves_aggregate_over_zero_update": bool(
            shared["best_damping_improves_over_zero_update"]
        ),
        "shared_damping_improves_every_cell": bool(
            all(row["shared_damping_improves_over_zero_update"] for row in cell_reports)
        ),
        "shared_damping_not_worse_than_full_step_every_cell": bool(
            all(row["shared_damping_not_worse_than_full_step"] for row in cell_reports)
        ),
        "proposed_aggregate_l1_update": proposed_l1,
        "shared_damped_aggregate_l1_update": float(shared_lambda * proposed_l1),
        "stationary_and_endpoint_damping_candidates": shared[
            "stationary_and_endpoint_damping_candidates"
        ],
        "objective_weighting": "equal weight per active radial node after concatenating all cells",
        "cells": cell_reports,
    }


def generate_actual_core_report(
    *,
    output: str | Path = "artifacts/kokuno_agent3/spacetime_damped_quadratic_gain_report.json",
    times: Iterable[float] = HELD_OUT_CYCLE_TIMES,
    z_values: Iterable[float] = SPATIAL_SCREEN_Z,
    radial_count: int = DEFAULT_RADIAL_COUNT,
    angular_count: int = DEFAULT_ANGULAR_COUNT,
    phase_count: int = DEFAULT_PHASE_COUNT,
) -> dict[str, Any]:
    """Run shared damping on the current real one-band negative control."""
    frozen_times = _validate_axis(times, name="times")
    frozen_z = _validate_axis(z_values, name="z_values")
    if REFERENCE_BUDGET_TIME not in frozen_times or REFERENCE_BUDGET_Z not in frozen_z:
        raise ValueError("the frozen reference cell (t=.5,z=.08) must be present")

    cells: list[tuple[float, float, MissingCovarianceColumnTarget, dict[str, Any]]] = []
    reference_receipt: MissingCovarianceColumnTarget | None = None
    leading_sha: str | None = None
    for z in frozen_z:
        for time in frozen_times:
            receipt, metadata = build_real_covariance_slice(
                time=time,
                z=z,
                radial_count=int(radial_count),
                angular_count=int(angular_count),
                phase_count=int(phase_count),
            )
            if leading_sha is None:
                leading_sha = str(metadata["leading_candidate_sha256"])
            elif leading_sha != metadata["leading_candidate_sha256"]:
                raise RuntimeError("leading candidate changed across spacetime reconstruction")
            if time == REFERENCE_BUDGET_TIME and z == REFERENCE_BUDGET_Z:
                reference_receipt = receipt
            cells.append((time, z, receipt, metadata))
    if reference_receipt is None:
        raise RuntimeError("reference budget receipt was not constructed")

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

    boundary_reports: list[dict[str, Any]] = []
    for sign, name in ((1.0, "positive_common_boundary"), (-1.0, "negative_common_boundary")):
        rows: list[dict[str, Any]] = []
        for time, z, receipt, _ in cells:
            measurement = measure_bounded_coordinate_quadratic_change(
                family,
                receipt.radii,
                coordinate_contract=coordinates,
                delta_common=sign * fractional_budget,
                delta_band=0.0,
                time=time,
                z=z,
                angular_count=int(angular_count),
                phase_count=int(phase_count),
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
        boundary_reports.append(
            {
                "name": name,
                "delta_common": float(sign * fractional_budget),
                "delta_band": 0.0,
                "shared_spacetime_damping": solve_shared_spacetime_damping(rows),
            }
        )

    target_cells = [
        {
            "time": time,
            "z": z,
            "active_target_nodes": int(np.count_nonzero(receipt.active_target_mask)),
            "nodes_requiring_second_direction": int(
                metadata["metrics"]["nodes_requiring_second_direction"]
            ),
            "missing_relative_vector_rms": float(
                metadata["metrics"]["missing_relative_vector_rms"]
            ),
        }
        for time, z, receipt, metadata in cells
    ]
    upstream_passed = bool(reference_inverse["family_bounded_inverse_preflight_passed"])
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
                "The corrected reader retains nonlinear correction self-covariance after its signed "
                "covariance inverse. One shared damping across the frozen 3x3 spacetime screen is "
                "repository finite-step engineering, not an additional Kokuno source formula."
            ),
        },
        "handoff": {
            "agent3_damped_gain_pr": AGENT3_DAMPED_GAIN_PR,
            "agent3_spacetime_preflight_pr": AGENT3_SPACETIME_PREFLIGHT_PR,
            "agent3_radial_side_effect_pr": AGENT3_RADIAL_SIDE_EFFECT_PR,
            "agent2_latest_dependency": (
                "PR #420 binds displayed auxiliary-torus constants but still lacks the actual "
                "positive-order/background multi-band physical velocity family and a genuinely "
                "independent second covariance direction"
            ),
        },
        "inputs": {
            "leading_candidate_sha256": leading_sha,
            "oscillatory_amplitude": ROUTED_OSCILLATORY_AMPLITUDE,
            "oscillatory_phase": ROUTED_OSCILLATORY_PHASE,
            "times": list(frozen_times),
            "z_values": list(frozen_z),
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
            "damping_interval": [0.0, 1.0],
            "same_damping_required_across_all_cells": True,
        },
        "real_spacetime_targets": {
            "cell_count": len(target_cells),
            "cells": target_cells,
        },
        "reference_one_band_bounded_inverse": reference_inverse,
        "boundary_shared_damping_diagnostics": boundary_reports,
        "routing": {
            "genuinely_independent_second_covariance_column_ready": False,
            "upstream_bounded_inverse_passed": upstream_passed,
            "public_velocity_correction_materialized": False,
            "radial_force_retained_in_materialized_correction": False,
            "finite_correction_cycle_rerun_allowed": False,
            "next_required": (
                "Agent 2 must first expose a source-motivated multi-band physical velocity family with "
                "an independent covariance direction. After rank/unit/budget/spacetime/radial guards "
                "pass, rerun this shared damping on the actual proposed correction before materializing "
                "and measuring the frozen held-in/held-out finite NS correction cycle."
            ),
        },
        "truth_boundary": {
            "real_candidate_defect_targets_consumed": True,
            "surrogate_defect_used": False,
            "shared_spacetime_quadratic_damping_executable": True,
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
        default="artifacts/kokuno_agent3/spacetime_damped_quadratic_gain_report.json",
    )
    parser.add_argument("--radial-count", type=int, default=DEFAULT_RADIAL_COUNT)
    parser.add_argument("--angular-count", type=int, default=DEFAULT_ANGULAR_COUNT)
    parser.add_argument("--phase-count", type=int, default=DEFAULT_PHASE_COUNT)
    return parser


def main() -> None:
    args = _parser().parse_args()
    report = generate_actual_core_report(
        output=args.output,
        radial_count=args.radial_count,
        angular_count=args.angular_count,
        phase_count=args.phase_count,
    )
    summary = {
        "task": report["task"],
        "upstream_bounded_inverse_passed": report["routing"]["upstream_bounded_inverse_passed"],
        "finite_correction_cycle_rerun_allowed": report["routing"][
            "finite_correction_cycle_rerun_allowed"
        ],
        "boundary_shared_damping": [
            {
                "name": row["name"],
                "selected_damping": row["shared_spacetime_damping"]["shared_selected_damping"],
                "best_relative_stress_residual_rms": row["shared_spacetime_damping"][
                    "shared_best_relative_stress_residual_rms"
                ],
                "improves_every_cell": row["shared_spacetime_damping"][
                    "shared_damping_improves_every_cell"
                ],
            }
            for row in report["boundary_shared_damping_diagnostics"]
        ],
    }
    print(json.dumps(summary, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
