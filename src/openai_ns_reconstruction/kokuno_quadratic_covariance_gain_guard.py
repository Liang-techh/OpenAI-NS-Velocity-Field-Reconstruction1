"""Retain nonlinear self-covariance before routing a finite correction cycle.

The corrected Kokuno reconstruction uses a signed *linear* covariance inverse,
but it also retains the nonlinear self-covariance of the correction.  Agent 3's
PR #393 currently screens a future multi-label family through the linearized
map ``[c0 J] da ~= target``.  A tiny linear algebraic residual is therefore only
a necessary preflight: after a physical update is chosen, the exact covariance
change also contains the quadratic term of that update.

For an assembled physical oscillatory field ``W`` and a bounded-coordinate
velocity update ``V`` this module measures, directly from public velocity values,

    Delta C = C(W + V) - C(W)
            = B(W, V) + C(V),

where for the two Agent-3 mean-stress channels

    B(W,V) = (<V_r W_theta + W_r V_theta>,
              <V_r W_z     + W_r V_z>),
    C(V)   = (<V_r V_theta>, <V_r V_z>).

The product expansion is repository algebra.  The source provenance is narrower:
the corrected reader keeps the nonlinear self-covariance separately after its
signed differential inverse.  This module does not add a new source formula,
does not construct Agent-2 oscillatory curls, and does not replace the final
held-in/held-out Navier--Stokes residual gate.

The routing rule deliberately reuses PR #393's existing algebraic stress-fit
tolerance.  Thus a linear bounded inverse is not enough: the *exact* quadratic
covariance change must still fit the same real stress target before a finite
cycle may be considered.  No PDE threshold is changed.
"""
from __future__ import annotations

import argparse
import json
from math import pi
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from .kokuno_bounded_coordinate_covariance import (
    EXPECTED_PARAMETER_NAMES,
    EXPECTED_PARAMETER_UNIT,
    FamilyPhaseEvaluator,
    _one_label_real_family,
    convert_amplitude_budget_to_fractional,
    covariance_tangent_from_velocity_and_tangents,
    measure_bounded_coordinate_covariance,
)
from .kokuno_bounded_family_coefficients import KokunoBoundedFamilyCoefficientCoordinates
from .kokuno_family_bounded_inverse import (
    _receipt_from_target_report,
    evaluate_family_bounded_inverse,
)
from .kokuno_missing_covariance_column_target import (
    DEFAULT_ANGULAR_COUNT,
    DEFAULT_PHASE_COUNT,
    PROFILE_Z,
    generate_actual_core_report as generate_missing_target_report,
)
from .kokuno_second_column_bounded_inverse import (
    ALGEBRAIC_RELATIVE_RESIDUAL_TOLERANCE,
    current_signed_coefficient_budget,
)
from .kokuno_signed_covariance_inverse import (
    PROFILE_TIME,
    ROUTED_OSCILLATORY_AMPLITUDE,
    ROUTED_OSCILLATORY_PHASE,
)
from .kokuno_covariance_tangent_preflight import _ring_points

TASK = "KOKUNO-A3-QUADRATIC-COVARIANCE-GAIN-GUARD-021"
SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
SOURCE_EDITION = "2026.09.09-consolidated"
SOURCE_DOI = "10.5281/zenodo.22678406"
SOURCE_SECTION = "signed covariance inverse; nonlinear self-covariance retained"
AGENT3_LINEAR_INVERSE_PR = 393
AGENT3_UNIT_BRIDGE_PR = 401
QUADRATIC_IDENTITY_RELATIVE_TOLERANCE = 5.0e-12


def _vector_rms(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    if values.ndim != 2 or values.shape[1] != 2 or not np.isfinite(values).all():
        raise ValueError("values must be finite with shape (N,2)")
    return float(np.sqrt(np.mean(np.sum(values * values, axis=1))))


def _covariance(values: np.ndarray, *, angular_axis: int = 1) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    if values.ndim < 2 or values.shape[-1] != 3 or not np.isfinite(values).all():
        raise ValueError("cylindrical velocity must be finite with final dimension 3")
    axis = int(angular_axis)
    if axis < 0:
        axis += values.ndim - 1
    if not 0 <= axis < values.ndim - 1:
        raise ValueError("angular_axis must index a sample dimension")
    radial = values[..., 0]
    theta = values[..., 1]
    axial = values[..., 2]
    return np.stack(
        (
            np.mean(radial * theta, axis=axis),
            np.mean(radial * axial, axis=axis),
        ),
        axis=-1,
    )


def measure_bounded_coordinate_quadratic_change(
    family_velocity: FamilyPhaseEvaluator,
    radii: np.ndarray,
    *,
    coordinate_contract: KokunoBoundedFamilyCoefficientCoordinates,
    delta_common: float,
    delta_band: float,
    time: float = PROFILE_TIME,
    z: float = PROFILE_Z,
    angular_count: int = DEFAULT_ANGULAR_COUNT,
    phase_count: int = DEFAULT_PHASE_COUNT,
) -> dict[str, Any]:
    """Measure exact, linear and self-covariance changes for one bounded update.

    The physical family is consumed through Agent 2/PR #400's declared by-beta
    contract.  No covariance column is inferred from labels alone.  The exact
    modulated total is evaluated by the coordinate contract and compared against
    the product-rule decomposition on the same physical-angle/phase quadrature.
    """
    if not callable(family_velocity):
        raise TypeError("family_velocity must be callable")
    if not isinstance(coordinate_contract, KokunoBoundedFamilyCoefficientCoordinates):
        raise TypeError("coordinate_contract must be KokunoBoundedFamilyCoefficientCoordinates")
    radii = np.asarray(radii, dtype=float)
    if radii.ndim != 1 or len(radii) < 3 or not np.isfinite(radii).all():
        raise ValueError("radii must be a finite one-dimensional array")
    if angular_count < 8 or angular_count % 2:
        raise ValueError("angular_count must be an even integer >= 8")
    if phase_count < 4 or phase_count % 2:
        raise ValueError("phase_count must be an even integer >= 4")
    if not np.isfinite([time, z, delta_common, delta_band]).all():
        raise ValueError("time, z and coordinate updates must be finite")

    points = _ring_points(radii, angular_count=angular_count, z=float(z))
    flat = points.reshape(-1, 3)
    shifts = 2.0 * pi * np.arange(phase_count, dtype=float) / phase_count
    delta = np.asarray([float(delta_common), float(delta_band)], dtype=float)

    base_covariance = np.zeros((len(radii), 2), dtype=float)
    modulated_covariance = np.zeros((len(radii), 2), dtype=float)
    linear_change = np.zeros((len(radii), 2), dtype=float)
    self_covariance = np.zeros((len(radii), 2), dtype=float)
    labels: tuple[str, ...] | None = None
    band_active: bool | None = None

    for shift in shifts:
        sample = dict(family_velocity(flat, float(time), float(shift)))
        out = coordinate_contract.consume_physical_family(
            sample,
            delta_common=float(delta_common),
            delta_band=float(delta_band),
        )
        names = tuple(str(x) for x in out["parameter_names"])
        units = tuple(str(x) for x in out["parameter_units"])
        if names != EXPECTED_PARAMETER_NAMES:
            raise ValueError("unexpected bounded-coordinate parameter names")
        if units != (EXPECTED_PARAMETER_UNIT, EXPECTED_PARAMETER_UNIT):
            raise ValueError("unexpected bounded-coordinate parameter units")
        current_labels = tuple(str(x) for x in out["beta_labels"])
        current_band_active = bool(out["band_contrast_active"])
        if labels is None:
            labels = current_labels
            band_active = current_band_active
        elif current_labels != labels or current_band_active != band_active:
            raise ValueError("bounded family labels/active-band state changed across phase samples")

        base = np.asarray(out["velocity_physical_cylindrical_total_base"], dtype=float)
        modulated = np.asarray(out["velocity_physical_cylindrical_total_modulated"], dtype=float)
        tangents = np.asarray(out["velocity_parameter_tangent_cylindrical"], dtype=float)
        base = base.reshape(len(radii), angular_count, 3)
        modulated = modulated.reshape(len(radii), angular_count, 3)
        tangents = tangents.reshape(len(radii), angular_count, 2, 3)
        update = modulated - base

        tangent_covariance = covariance_tangent_from_velocity_and_tangents(
            base, tangents, sample_axis=1
        )
        linear_change += np.sum(tangent_covariance * delta[None, :, None], axis=1)
        self_covariance += _covariance(update, angular_axis=1)
        base_covariance += _covariance(base, angular_axis=1)
        modulated_covariance += _covariance(modulated, angular_axis=1)

    scale = float(phase_count)
    base_covariance /= scale
    modulated_covariance /= scale
    linear_change /= scale
    self_covariance /= scale
    exact_change = modulated_covariance - base_covariance
    closure = exact_change - linear_change - self_covariance
    denominator = max(
        _vector_rms(exact_change),
        _vector_rms(linear_change + self_covariance),
        np.finfo(float).tiny,
    )
    closure_relative = _vector_rms(closure) / denominator
    identity_passed = bool(closure_relative <= QUADRATIC_IDENTITY_RELATIVE_TOLERANCE)
    if not identity_passed:
        raise RuntimeError("exact covariance change failed the quadratic product identity")

    linear_rms = _vector_rms(linear_change)
    self_rms = _vector_rms(self_covariance)
    assert labels is not None and band_active is not None
    return {
        "radii": radii.tolist(),
        "beta_labels": list(labels),
        "parameter_names": list(EXPECTED_PARAMETER_NAMES),
        "parameter_unit": EXPECTED_PARAMETER_UNIT,
        "band_contrast_active": band_active,
        "delta_common": float(delta_common),
        "delta_band": float(delta_band),
        "aggregate_l1_update": float(abs(delta_common) + abs(delta_band)),
        "base_covariance_theta_axial": base_covariance.tolist(),
        "linear_covariance_change_theta_axial": linear_change.tolist(),
        "self_covariance_theta_axial": self_covariance.tolist(),
        "exact_covariance_change_theta_axial": exact_change.tolist(),
        "modulated_covariance_theta_axial": modulated_covariance.tolist(),
        "quadratic_identity_closure_theta_axial": closure.tolist(),
        "quadratic_identity_relative_rms": float(closure_relative),
        "quadratic_identity_passed": identity_passed,
        "linear_change_vector_rms": float(linear_rms),
        "self_covariance_vector_rms": float(self_rms),
        "self_to_linear_vector_rms_ratio": (
            float(self_rms / linear_rms) if linear_rms > np.finfo(float).tiny else None
        ),
    }


def evaluate_quadratic_covariance_gain(
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
    """Fail closed unless the exact quadratic covariance still fits the target.

    ``target_stress`` is the real compact theta/axial stress target already used
    by the linear bounded inverse.  The pre-correction residual is therefore the
    target itself.  The exact post-covariance residual is

        target - (linear_change + self_covariance).

    Passing this local guard is never sufficient for PDE acceptance; it only
    allows downstream finite-cycle routing when all separately declared upstream
    prerequisites are also true.
    """
    target = np.asarray(target_stress, dtype=float)
    active = np.asarray(active_target_mask, dtype=bool)
    if target.ndim != 2 or target.shape[1] != 2 or not np.isfinite(target).all():
        raise ValueError("target_stress must be finite with shape (N,2)")
    if active.shape != (len(target),) or not np.any(active):
        raise ValueError("active_target_mask must select at least one target node")
    linear = np.asarray(measurement["linear_covariance_change_theta_axial"], dtype=float)
    self_term = np.asarray(measurement["self_covariance_theta_axial"], dtype=float)
    exact = np.asarray(measurement["exact_covariance_change_theta_axial"], dtype=float)
    if linear.shape != target.shape or self_term.shape != target.shape or exact.shape != target.shape:
        raise ValueError("measurement covariance arrays must match target_stress")
    if not np.isfinite(linear).all() or not np.isfinite(self_term).all() or not np.isfinite(exact).all():
        raise ValueError("measurement covariance arrays must be finite")
    tolerance = float(algebraic_relative_residual_tolerance)
    if not np.isfinite(tolerance) or not 0.0 < tolerance < 1.0e-4:
        raise ValueError("algebraic_relative_residual_tolerance must lie in (0,1e-4)")

    flags = (
        upstream_bounded_inverse_passed,
        spacetime_preflight_passed,
        radial_force_retained,
        public_velocity_correction_materialized,
    )
    if any(not isinstance(flag, (bool, np.bool_)) for flag in flags):
        raise TypeError("routing prerequisite flags must be boolean")

    target_active = target[active]
    target_rms = max(_vector_rms(target_active), np.finfo(float).tiny)
    linear_residual = target_active - linear[active]
    exact_residual = target_active - exact[active]
    linear_relative = _vector_rms(linear_residual) / target_rms
    exact_relative = _vector_rms(exact_residual) / target_rms
    self_relative = _vector_rms(self_term[active]) / target_rms
    identity_passed = bool(measurement.get("quadratic_identity_passed", False))
    exact_fit_passed = bool(exact_relative <= tolerance)
    exact_improves = bool(exact_relative < 1.0)
    local_passed = bool(
        identity_passed and upstream_bounded_inverse_passed and exact_fit_passed
    )
    finite_allowed = bool(
        local_passed
        and spacetime_preflight_passed
        and radial_force_retained
        and public_velocity_correction_materialized
    )
    return {
        "algebraic_relative_residual_tolerance": tolerance,
        "target_vector_rms": float(target_rms),
        "linear_predicted_relative_stress_residual_rms": float(linear_relative),
        "exact_quadratic_relative_stress_residual_rms": float(exact_relative),
        "self_covariance_relative_to_target_rms": float(self_relative),
        "exact_covariance_improves_over_zero_update": exact_improves,
        "quadratic_identity_passed": identity_passed,
        "exact_quadratic_stress_fit_within_existing_algebraic_tolerance": exact_fit_passed,
        "upstream_bounded_inverse_passed": bool(upstream_bounded_inverse_passed),
        "spacetime_preflight_passed": bool(spacetime_preflight_passed),
        "radial_force_retained": bool(radial_force_retained),
        "public_velocity_correction_materialized": bool(public_velocity_correction_materialized),
        "quadratic_covariance_preflight_passed": local_passed,
        "finite_correction_cycle_rerun_allowed": finite_allowed,
        "interpretation": (
            "The exact quadratic covariance target check reuses the existing linear algebraic tolerance. "
            "Even a pass remains only a preflight; retained radial d_z sigma_1 and a fresh held-in/held-out "
            "Navier-Stokes residual cycle remain mandatory."
        ),
    }


def generate_actual_core_report(
    *,
    output: str | Path = "artifacts/kokuno_agent3/quadratic_covariance_gain_guard_report.json",
    target_output: str | Path = "artifacts/kokuno_agent3/quadratic_covariance_gain_guard_target.json",
) -> dict[str, Any]:
    """Calibrate the exact quadratic guard on the current real one-band family.

    The current family is intentionally a negative control: it has no active band
    direction, so PR #401/#393's upstream bounded inverse remains false.  We do
    not choose a favorable step.  Both signs of the already-frozen fractional
    budget boundary are measured only to quantify the nonlinear self-covariance
    scale in the declared coordinate units.
    """
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
    theory_ratio = 0.5 * fractional_budget
    for sign, name in ((1.0, "positive_common_boundary"), (-1.0, "negative_common_boundary")):
        delta_common = sign * fractional_budget
        measurement = measure_bounded_coordinate_quadratic_change(
            family,
            receipt.radii,
            coordinate_contract=coordinates,
            delta_common=delta_common,
            delta_band=0.0,
            time=float(target_report["inputs"]["profile_time"]),
            z=float(target_report["inputs"]["profile_z"]),
        )
        routing = evaluate_quadratic_covariance_gain(
            receipt.target_stress,
            receipt.active_target_mask,
            measurement,
            upstream_bounded_inverse_passed=upstream_passed,
            spacetime_preflight_passed=False,
            radial_force_retained=False,
            public_velocity_correction_materialized=False,
        )
        measured_ratio = measurement["self_to_linear_vector_rms_ratio"]
        if measured_ratio is None:
            raise RuntimeError("common-coordinate calibration unexpectedly has zero linear response")
        boundary_rows.append(
            {
                "name": name,
                "measurement": measurement,
                "routing": routing,
                "common_scaling_exact_self_to_linear_theory": float(theory_ratio),
                "common_scaling_self_to_linear_absolute_error": float(
                    abs(float(measured_ratio) - theory_ratio)
                ),
            }
        )

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
                "The corrected reader retains nonlinear correction self-covariance after the signed linear "
                "covariance inverse. The finite-dimensional bounded-coordinate implementation and exact product "
                "expansion diagnostics here are repository machinery."
            ),
        },
        "handoff": {
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
            "quadratic_identity_relative_tolerance": QUADRATIC_IDENTITY_RELATIVE_TOLERANCE,
            "algebraic_relative_residual_tolerance_reused_from_pr393": ALGEBRAIC_RELATIVE_RESIDUAL_TOLERANCE,
        },
        "upstream_linear_bounded_inverse": upstream_inverse,
        "boundary_diagnostics": boundary_rows,
        "truth_boundary": {
            "real_candidate_defect_target_consumed": True,
            "surrogate_defect_used": False,
            "nonlinear_self_covariance_retained_in_executable_guard": True,
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
        default="artifacts/kokuno_agent3/quadratic_covariance_gain_guard_report.json",
    )
    parser.add_argument(
        "--target-output",
        default="artifacts/kokuno_agent3/quadratic_covariance_gain_guard_target.json",
    )
    args = parser.parse_args()
    report = generate_actual_core_report(output=args.output, target_output=args.target_output)
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
