"""Choose a bounded damped covariance step without dropping Kokuno's quadratic term.

Agent 3 PR #411 makes the retained nonlinear correction self-covariance executable.
For one already-proposed physical bounded update ``V`` it measures

    Delta C(1) = B(W,V) + C(V)

against the real compact theta/axial stress target.  A full step can overshoot
because the self-covariance is quadratic.  The corrected Kokuno reader motivates
retaining that nonlinear term, but it does not prescribe this repository's
finite-dimensional damping rule.

This module adds only that engineering rule.  For a scalar damping
``lambda in [0,1]`` along a *fixed* proposed update direction,

    Delta C(lambda) = lambda B(W,V) + lambda^2 C(V).

Hence the squared target mismatch is a quartic polynomial in ``lambda``.  We
find its exact interval minimizer by evaluating the endpoints and every real
stationary point from the cubic derivative.  No stress/PDE threshold is changed,
no new oscillatory direction is synthesized, and a local covariance gain never
replaces the later held-in/held-out Navier--Stokes residual cycle.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from .kokuno_bounded_coordinate_covariance import (
    EXPECTED_PARAMETER_NAMES,
    EXPECTED_PARAMETER_UNIT,
    _one_label_real_family,
    convert_amplitude_budget_to_fractional,
    measure_bounded_coordinate_covariance,
)
from .kokuno_bounded_family_coefficients import KokunoBoundedFamilyCoefficientCoordinates
from .kokuno_family_bounded_inverse import (
    _receipt_from_target_report,
    evaluate_family_bounded_inverse,
)
from .kokuno_missing_covariance_column_target import (
    generate_actual_core_report as generate_missing_target_report,
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
    ROUTED_OSCILLATORY_AMPLITUDE,
    ROUTED_OSCILLATORY_PHASE,
)

TASK = "KOKUNO-A3-DAMPED-QUADRATIC-COVARIANCE-GAIN-022"
SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
SOURCE_EDITION = "2026.09.09-consolidated"
SOURCE_DOI = "10.5281/zenodo.22678406"
SOURCE_SECTION = "signed covariance inverse; nonlinear self-covariance retained"
AGENT3_QUADRATIC_GAIN_PR = 411
AGENT3_LINEAR_INVERSE_PR = 393
AGENT3_UNIT_BRIDGE_PR = 401
ROOT_IMAGINARY_TOLERANCE = 1.0e-10
ROOT_INTERVAL_TOLERANCE = 1.0e-12


def _vector_rms(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    if values.ndim != 2 or values.shape[1] != 2 or not np.isfinite(values).all():
        raise ValueError("values must be finite with shape (N,2)")
    return float(np.sqrt(np.mean(np.sum(values * values, axis=1))))


def _mean_inner(left: np.ndarray, right: np.ndarray) -> float:
    left = np.asarray(left, dtype=float)
    right = np.asarray(right, dtype=float)
    if left.shape != right.shape or left.ndim != 2 or left.shape[1] != 2:
        raise ValueError("inner-product inputs must share shape (N,2)")
    if not np.isfinite(left).all() or not np.isfinite(right).all():
        raise ValueError("inner-product inputs must be finite")
    return float(np.mean(np.sum(left * right, axis=1)))


def _quartic_coefficients(
    target: np.ndarray,
    linear: np.ndarray,
    quadratic: np.ndarray,
) -> np.ndarray:
    """Return low-to-high coefficients of ||T-lambda L-lambda^2 Q||_RMS^2."""
    c0 = _mean_inner(target, target)
    c1 = -2.0 * _mean_inner(target, linear)
    c2 = _mean_inner(linear, linear) - 2.0 * _mean_inner(target, quadratic)
    c3 = 2.0 * _mean_inner(linear, quadratic)
    c4 = _mean_inner(quadratic, quadratic)
    return np.asarray([c0, c1, c2, c3, c4], dtype=float)


def _polyval_low_to_high(coefficients: np.ndarray, value: float) -> float:
    coefficients = np.asarray(coefficients, dtype=float)
    x = float(value)
    total = 0.0
    for coefficient in coefficients[::-1]:
        total = total * x + float(coefficient)
    return float(total)


def _stationary_candidates(coefficients: np.ndarray) -> list[float]:
    coefficients = np.asarray(coefficients, dtype=float)
    if coefficients.shape != (5,) or not np.isfinite(coefficients).all():
        raise ValueError("quartic coefficients must be finite with shape (5,)")
    # Derivative coefficients in descending order for numpy.roots.
    derivative = np.asarray(
        [
            4.0 * coefficients[4],
            3.0 * coefficients[3],
            2.0 * coefficients[2],
            coefficients[1],
        ],
        dtype=float,
    )
    scale = max(float(np.max(np.abs(derivative))), 1.0)
    cutoff = 32.0 * np.finfo(float).eps * scale
    first = 0
    while first < len(derivative) and abs(float(derivative[first])) <= cutoff:
        first += 1

    candidates = [0.0, 1.0]
    if first < len(derivative):
        for root in np.roots(derivative[first:]):
            if abs(float(np.imag(root))) > ROOT_IMAGINARY_TOLERANCE:
                continue
            value = float(np.real(root))
            if -ROOT_INTERVAL_TOLERANCE <= value <= 1.0 + ROOT_INTERVAL_TOLERANCE:
                candidates.append(float(np.clip(value, 0.0, 1.0)))

    candidates.sort()
    unique: list[float] = []
    for candidate in candidates:
        if not unique or abs(candidate - unique[-1]) > ROOT_INTERVAL_TOLERANCE:
            unique.append(candidate)
    return unique


def solve_damped_quadratic_covariance_gain(
    target_stress: np.ndarray,
    active_target_mask: np.ndarray,
    measurement: Mapping[str, Any],
) -> dict[str, Any]:
    """Minimize the exact quadratic covariance mismatch along one fixed update.

    ``measurement`` must be the full-step output of PR #411's
    ``measure_bounded_coordinate_quadratic_change``.  Damping scales the already
    proposed update direction; it never creates a new covariance column.
    """
    target = np.asarray(target_stress, dtype=float)
    active = np.asarray(active_target_mask, dtype=bool)
    if target.ndim != 2 or target.shape[1] != 2 or not np.isfinite(target).all():
        raise ValueError("target_stress must be finite with shape (N,2)")
    if active.shape != (len(target),) or not np.any(active):
        raise ValueError("active_target_mask must select at least one node")

    linear = np.asarray(measurement["linear_covariance_change_theta_axial"], dtype=float)
    quadratic = np.asarray(measurement["self_covariance_theta_axial"], dtype=float)
    exact = np.asarray(measurement["exact_covariance_change_theta_axial"], dtype=float)
    if linear.shape != target.shape or quadratic.shape != target.shape or exact.shape != target.shape:
        raise ValueError("measurement covariance arrays must match target_stress")
    if not np.isfinite(linear).all() or not np.isfinite(quadratic).all() or not np.isfinite(exact).all():
        raise ValueError("measurement covariance arrays must be finite")
    if not bool(measurement.get("quadratic_identity_passed", False)):
        raise ValueError("measurement must pass PR #411 quadratic identity closure")

    expected_exact = linear + quadratic
    exact_denominator = max(_vector_rms(expected_exact[active]), np.finfo(float).tiny)
    exact_closure_relative = _vector_rms(exact[active] - expected_exact[active]) / exact_denominator
    if exact_closure_relative > QUADRATIC_IDENTITY_RELATIVE_TOLERANCE:
        raise ValueError("measurement exact change is inconsistent with linear+self covariance")

    target_active = target[active]
    linear_active = linear[active]
    quadratic_active = quadratic[active]
    coefficients = _quartic_coefficients(target_active, linear_active, quadratic_active)
    candidates = _stationary_candidates(coefficients)
    scored = [
        (max(_polyval_low_to_high(coefficients, value), 0.0), value)
        for value in candidates
    ]
    best_objective, best_damping = min(scored, key=lambda row: (row[0], row[1]))

    target_rms = max(_vector_rms(target_active), np.finfo(float).tiny)
    best_relative = float(np.sqrt(best_objective) / target_rms)
    full_objective = max(_polyval_low_to_high(coefficients, 1.0), 0.0)
    full_relative = float(np.sqrt(full_objective) / target_rms)
    zero_relative = 1.0
    best_change = best_damping * linear + best_damping * best_damping * quadratic
    best_residual = target - best_change

    proposed_l1 = measurement.get("aggregate_l1_update")
    damped_l1: float | None
    if proposed_l1 is None:
        damped_l1 = None
    else:
        proposed_l1 = float(proposed_l1)
        if not np.isfinite(proposed_l1) or proposed_l1 < 0.0:
            raise ValueError("aggregate_l1_update must be finite and nonnegative")
        damped_l1 = float(best_damping * proposed_l1)

    return {
        "quartic_objective_coefficients_low_to_high": coefficients.tolist(),
        "stationary_and_endpoint_damping_candidates": [float(x) for x in candidates],
        "selected_damping": float(best_damping),
        "zero_update_relative_stress_residual_rms": zero_relative,
        "full_step_relative_stress_residual_rms": float(full_relative),
        "best_damped_relative_stress_residual_rms": float(best_relative),
        "best_damped_covariance_change_theta_axial": best_change.tolist(),
        "best_damped_stress_residual_theta_axial": best_residual.tolist(),
        "best_damping_improves_over_zero_update": bool(best_relative < zero_relative),
        "best_damping_not_worse_than_full_step": bool(best_relative <= full_relative + 64.0 * np.finfo(float).eps),
        "damping_reduces_step": bool(best_damping < 1.0 - ROOT_INTERVAL_TOLERANCE),
        "proposed_aggregate_l1_update": proposed_l1,
        "damped_aggregate_l1_update": damped_l1,
        "quadratic_identity_rechecked_relative_rms": float(exact_closure_relative),
    }


def evaluate_damped_quadratic_covariance_gain(
    target_stress: np.ndarray,
    active_target_mask: np.ndarray,
    measurement: Mapping[str, Any],
    *,
    upstream_bounded_inverse_passed: bool,
    spacetime_preflight_passed: bool,
    radial_force_retained: bool,
    public_velocity_correction_materialized: bool,
    algebraic_relative_residual_tolerance: float = ALGEBRAIC_RELATIVE_RESIDUAL_TOLERANCE,
) -> dict[str, Any]:
    """Route only a genuinely improving damped exact-covariance fit.

    The algebraic tolerance is exactly the pre-existing PR #393 tolerance.  The
    damping minimization only chooses a smaller step along an already proposed
    bounded direction; it does not relax coefficient or PDE acceptance gates.
    """
    flags = (
        upstream_bounded_inverse_passed,
        spacetime_preflight_passed,
        radial_force_retained,
        public_velocity_correction_materialized,
    )
    if any(not isinstance(flag, (bool, np.bool_)) for flag in flags):
        raise TypeError("routing prerequisite flags must be boolean")
    tolerance = float(algebraic_relative_residual_tolerance)
    if not np.isfinite(tolerance) or not 0.0 < tolerance < 1.0e-4:
        raise ValueError("algebraic_relative_residual_tolerance must lie in (0,1e-4)")

    solved = solve_damped_quadratic_covariance_gain(
        target_stress,
        active_target_mask,
        measurement,
    )
    best_relative = float(solved["best_damped_relative_stress_residual_rms"])
    positive_step = bool(float(solved["selected_damping"]) > ROOT_INTERVAL_TOLERANCE)
    exact_fit_passed = bool(best_relative <= tolerance)
    gain_passed = bool(solved["best_damping_improves_over_zero_update"])
    local_passed = bool(
        upstream_bounded_inverse_passed
        and positive_step
        and gain_passed
        and exact_fit_passed
    )
    finite_allowed = bool(
        local_passed
        and spacetime_preflight_passed
        and radial_force_retained
        and public_velocity_correction_materialized
    )
    return {
        **solved,
        "algebraic_relative_residual_tolerance": tolerance,
        "positive_nontrivial_damped_step": positive_step,
        "best_damped_exact_stress_fit_within_existing_algebraic_tolerance": exact_fit_passed,
        "upstream_bounded_inverse_passed": bool(upstream_bounded_inverse_passed),
        "spacetime_preflight_passed": bool(spacetime_preflight_passed),
        "radial_force_retained": bool(radial_force_retained),
        "public_velocity_correction_materialized": bool(public_velocity_correction_materialized),
        "damped_quadratic_covariance_preflight_passed": local_passed,
        "finite_correction_cycle_rerun_allowed": finite_allowed,
        "interpretation": (
            "The selected damping is the exact quartic covariance-space minimizer along one fixed proposed "
            "bounded update. It does not create a second direction or replace the retained radial force and "
            "fresh held-in/held-out Navier-Stokes residual cycle."
        ),
    }


def generate_actual_core_report(
    *,
    output: str | Path = "artifacts/kokuno_agent3/damped_quadratic_covariance_gain_report.json",
    target_output: str | Path = "artifacts/kokuno_agent3/damped_quadratic_covariance_gain_target.json",
) -> dict[str, Any]:
    """Run the damping machinery on the current real one-band negative control."""
    target_report = generate_missing_target_report(output=target_output)
    receipt = _receipt_from_target_report(target_report)
    frozen_amplitude_budget = current_signed_coefficient_budget(receipt)
    fractional_budget = convert_amplitude_budget_to_fractional(
        frozen_amplitude_budget, ROUTED_OSCILLATORY_AMPLITUDE
    )
    coordinates = KokunoBoundedFamilyCoefficientCoordinates(max_l1_update=fractional_budget)
    family = _one_label_real_family()

    linear_measurement = measure_bounded_coordinate_covariance(
        family,
        receipt.radii,
        coordinate_contract=coordinates,
        time=float(target_report["inputs"]["profile_time"]),
        z=float(target_report["inputs"]["profile_z"]),
    )
    jacobian = np.asarray(
        linear_measurement["coordinate_covariance_jacobian_theta_axial"], dtype=float
    )
    upstream_inverse = evaluate_family_bounded_inverse(
        receipt,
        jacobian,
        coefficient_labels=EXPECTED_PARAMETER_NAMES,
        current_coefficient_budget=frozen_amplitude_budget,
        additional_family_l1_budget=fractional_budget,
        budget_unit_matches_jacobian=True,
        coefficient_unit=EXPECTED_PARAMETER_UNIT,
    )
    upstream_passed = bool(upstream_inverse["family_bounded_inverse_preflight_passed"])

    boundary_rows: list[dict[str, Any]] = []
    for sign, name in ((1.0, "positive_common_boundary"), (-1.0, "negative_common_boundary")):
        measurement = measure_bounded_coordinate_quadratic_change(
            family,
            receipt.radii,
            coordinate_contract=coordinates,
            delta_common=sign * fractional_budget,
            delta_band=0.0,
            time=float(target_report["inputs"]["profile_time"]),
            z=float(target_report["inputs"]["profile_z"]),
        )
        routing = evaluate_damped_quadratic_covariance_gain(
            receipt.target_stress,
            receipt.active_target_mask,
            measurement,
            upstream_bounded_inverse_passed=upstream_passed,
            spacetime_preflight_passed=False,
            radial_force_retained=False,
            public_velocity_correction_materialized=False,
        )
        boundary_rows.append({"name": name, "measurement": measurement, "damped_gain": routing})

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
                "The corrected reader retains nonlinear correction self-covariance after its signed linear "
                "covariance inverse. Scalar damping and exact quartic interval minimization are repository "
                "finite-step engineering, not an additional Kokuno formula."
            ),
        },
        "handoff": {
            "agent3_quadratic_gain_pr": AGENT3_QUADRATIC_GAIN_PR,
            "agent3_linear_inverse_pr": AGENT3_LINEAR_INVERSE_PR,
            "agent3_unit_bridge_pr": AGENT3_UNIT_BRIDGE_PR,
            "agent2_coordinate_contract": "PR #400 bounded physical-family coordinates",
        },
        "inputs": {
            "leading_candidate_sha256": target_report["inputs"]["leading_candidate_sha256"],
            "oscillatory_amplitude": ROUTED_OSCILLATORY_AMPLITUDE,
            "oscillatory_phase": ROUTED_OSCILLATORY_PHASE,
            "profile_time": target_report["inputs"]["profile_time"],
            "profile_z": target_report["inputs"]["profile_z"],
            "profile_annulus": target_report["inputs"]["profile_annulus"],
            "surrogate_defect_used": False,
            "frozen_physical_amplitude_budget": float(frozen_amplitude_budget),
            "converted_fractional_budget": float(fractional_budget),
            "budget_changed": False,
            "algebraic_relative_residual_tolerance_reused_from_pr393": ALGEBRAIC_RELATIVE_RESIDUAL_TOLERANCE,
            "damping_interval": [0.0, 1.0],
        },
        "upstream_linear_bounded_inverse": upstream_inverse,
        "boundary_damped_gain_diagnostics": boundary_rows,
        "truth_boundary": {
            "real_candidate_defect_target_consumed": True,
            "surrogate_defect_used": False,
            "damped_quadratic_covariance_gain_executable": True,
            "actual_source_multiband_family_consumed": False,
            "genuinely_independent_second_covariance_column_ready": False,
            "upstream_bounded_inverse_passed": upstream_passed,
            "public_velocity_correction_materialized": False,
            "radial_force_retained_in_materialized_correction": False,
            "finite_correction_cycle_run": False,
            "finite_correction_cycle_rerun_allowed": False,
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
    destination.write_text(json.dumps(report, indent=2, sort_keys=True), encoding="utf-8")
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        default="artifacts/kokuno_agent3/damped_quadratic_covariance_gain_report.json",
    )
    parser.add_argument(
        "--target-output",
        default="artifacts/kokuno_agent3/damped_quadratic_covariance_gain_target.json",
    )
    args = parser.parse_args()
    report = generate_actual_core_report(output=args.output, target_output=args.target_output)
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
