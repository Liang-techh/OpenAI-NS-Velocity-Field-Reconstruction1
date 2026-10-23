"""Target-free joint-energy-renormalization audit for the axial-cap poloidal mode.

PR #327 showed that the fixed support-adjacent axial-cap mode is the first small
basis direction in this lane that can move gross q90/q99 vorticity reach, but
that the useful +/-0.25 capacity probe lies outside the unchanged fixed-parent
reference-energy gate.  This module asks exactly one next representation
question without adding another basis: if the entire axial-cap trial velocity is
multiplied by one positive scalar so that the reference kinetic energy is
restored exactly, does the useful finite-amplitude cap geometry fit inside the
inherited Eq45 coefficient guard?

    lambda(a) = sqrt(E_target / E_raw(a))
    u_joint   = lambda(a) * (u_parent + a * delta_u_cap)

A common positive scale preserves support, parity, instantaneous streamline
directions, and normalized vorticity-location fingerprints.  It is not a
Navier--Stokes invariance, does not transfer any pressure/forcing result, and is
not materialized here as a new candidate.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Iterable

import numpy as np

from .constrained_eq45_bipolar_axial_cap_energy_envelope import (
    _constraints,
    _energy,
    _energy_quadratic,
    _refinement,
)
from .constrained_eq45_bipolar_axial_cap_poloidal_capacity import (
    AXIAL_CAP_POLOIDAL,
    _cap_velocity,
    _cap_vorticity_morphology,
    _source_field,
)

TASK_ID = "CR003-BIPOLAR-AXIAL-CAP-JOINT-RENORM-CAPACITY-047"
DEFAULT_TRIAL_MAGNITUDE = 0.25
TRUTH_BOUNDARY = {
    "velocity_changed": False,
    "candidate_artifact_changed": False,
    "new_basis_added": False,
    "joint_renormalized_child_materialized": False,
    "axial_cap_poloidal_coefficient_selected": False,
    "materialization_bound_selected": False,
    "candidate_sha_created": False,
    "pressure_or_force_changed": False,
    "held_out_pde_residual_evaluated": False,
    "public_image_used": False,
    "visualization_ready": False,
    "visual_correspondence_verified": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
}


def _normalization_scale(coefficients: dict[str, float], coefficient: float, target: float) -> float:
    raw = _energy(coefficients, float(coefficient))
    if not np.isfinite(raw) or raw <= 0.0 or not np.isfinite(target) or target <= 0.0:
        raise ValueError("energy normalization requires positive finite energies")
    return float(np.sqrt(float(target) / raw))


def _profile_bound_screen(field, scale: float, cap_coefficient: float) -> dict[str, Any]:
    """Screen the simplest coefficient-space realization of the common scale.

    The existing Eq45 velocity is linear in its stored Phi/F coefficients.  A
    common velocity scale therefore scales those stored coefficients too.  The
    cap coefficient is not yet a production parameter; its scaled value is only
    checked against the same inherited implementation guard used by the capacity
    module.  This is a representation preflight, not candidate acceptance.
    """
    basis = field.parent.profile_basis
    limit = float(basis.coefficient_limit)
    phi = np.asarray(basis.phi_coefficients, dtype=float)
    swirl = np.asarray(basis.swirl_coefficients, dtype=float)
    stored = np.concatenate((phi, swirl))
    scaled = float(scale) * stored
    effective_cap = float(scale) * float(cap_coefficient)
    max_before = float(np.max(np.abs(stored))) if stored.size else 0.0
    max_after = float(np.max(np.abs(scaled))) if scaled.size else 0.0
    tolerance = 64.0 * np.finfo(float).eps * max(1.0, limit)
    stored_ok = bool(max_after <= limit + tolerance)
    cap_guard_ok = bool(abs(effective_cap) <= limit + tolerance)
    return {
        "inherited_coefficient_limit": limit,
        "maximum_abs_stored_profile_coefficient_before_scale": max_before,
        "maximum_abs_stored_profile_coefficient_after_scale": max_after,
        "effective_axial_cap_coefficient_after_common_scale": effective_cap,
        "stored_profile_coefficients_within_inherited_limit": stored_ok,
        "effective_axial_cap_coefficient_within_inherited_implementation_guard": cap_guard_ok,
        "simple_common_scale_representation_preflight_passed": bool(stored_ok and cap_guard_ok),
        "scope": (
            "coefficient-space preflight only; a materialized jointly renormalized child needs a new identity, an explicit "
            "production cap bound/value, and fresh energy/sign/support/divergence/full-momentum validation"
        ),
    }


def _core_sign_screen(field, coefficient: float, scale: float, times: Iterable[float]) -> dict[str, Any]:
    rows = []
    for time in tuple(float(value) for value in times):
        tau = 1.0 - time
        point = np.array([[0.1 * np.sqrt(tau), 0.0, 0.1 * tau**0.495]], dtype=float)
        raw = np.asarray(_cap_velocity(field, point, time, coefficient)[0], dtype=float)
        velocity = float(scale) * raw
        rows.append({
            "time": time,
            "u_r": float(velocity[0]),
            "u_theta": float(velocity[1]),
            "u_z": float(velocity[2]),
            "signs_satisfied": bool(velocity[0] < 0.0 and velocity[1] > 0.0 and velocity[2] > 0.0),
        })
    return {
        "all_times_satisfied": bool(all(row["signs_satisfied"] for row in rows)),
        "positive_common_scale_preserves_signs_analytically": bool(scale > 0.0),
        "rows": rows,
    }


def _morphology_delta(base: dict[str, float], other: dict[str, float]) -> dict[str, float]:
    keys = (
        "axial_rms_over_Zp",
        "radial_rms_over_Rp",
        "axial_q90_over_Zp",
        "axial_q99_over_Zp",
        "outer_065_enstrophy_fraction",
        "outer_075_enstrophy_fraction",
    )
    return {f"{key}_delta": float(other[key] - base[key]) for key in keys}


def _trial_row(
    field,
    *,
    coefficient: float,
    reference_coefficients: dict[str, float],
    energy_by_time: dict[float, dict[str, float]],
    target: float,
    validation_energy_lower: float,
    validation_energy_upper: float,
    times: tuple[float, ...],
    morphology_time: float,
    morphology_grid_size: int,
    baseline_morphology: dict[str, float],
) -> dict[str, Any]:
    coefficient = float(coefficient)
    scale = _normalization_scale(reference_coefficients, coefficient, target)
    raw_reference_energy = _energy(reference_coefficients, coefficient)

    validation_rows = []
    all_validation_energy = True
    for time in times:
        raw = _energy(energy_by_time[time], coefficient)
        normalized = scale * scale * raw
        passed = bool(validation_energy_lower <= normalized <= validation_energy_upper)
        all_validation_energy &= passed
        validation_rows.append({
            "time": time,
            "raw_energy": float(raw),
            "jointly_renormalized_energy": float(normalized),
            "validation_energy_range_satisfied": passed,
        })

    morphology = _cap_vorticity_morphology(
        field,
        time=float(morphology_time),
        coefficient=coefficient,
        grid_size=int(morphology_grid_size),
    )
    delta = _morphology_delta(baseline_morphology, morphology)
    return {
        "coefficient_before_common_scale": coefficient,
        "raw_reference_energy": float(raw_reference_energy),
        "common_velocity_scale": scale,
        "common_scale_fractional_change": float(scale - 1.0),
        "jointly_renormalized_reference_energy": float(scale * scale * raw_reference_energy),
        "profile_bound_screen": _profile_bound_screen(field, scale, coefficient),
        "core_sign_screen": _core_sign_screen(field, coefficient, scale, times),
        "validation_time_energy": {
            "all_times_satisfied": bool(all_validation_energy),
            "required_range": [validation_energy_lower, validation_energy_upper],
            "rows": validation_rows,
        },
        "vorticity_morphology": morphology,
        "morphology_delta_from_baseline": delta,
        "gross_tip_reach_extended_on_this_grid": bool(
            delta["axial_q90_over_Zp_delta"] > 0.0 and delta["axial_q99_over_Zp_delta"] > 0.0
        ),
        "morphology_invariance_under_positive_common_scale": {
            "instantaneous_streamline_geometry_invariant": True,
            "normalized_vorticity_location_fingerprints_invariant": True,
            "reason": (
                "at fixed time the common positive scale multiplies velocity and vorticity uniformly, so streamline directions "
                "and normalized enstrophy-location metrics are unchanged from the raw cap trial"
            ),
            "navier_stokes_balance_not_invariant": True,
            "pathline_time_parameterization_not_claimed_invariant": True,
        },
    }


def audit_axial_cap_joint_renormalization(
    *,
    trial_magnitude: float = DEFAULT_TRIAL_MAGNITUDE,
    quadrature_orders: Iterable[int] = (96, 192),
    morphology_time: float = 0.5,
    morphology_grid_size: int = 33,
) -> dict[str, Any]:
    magnitude = float(trial_magnitude)
    if not np.isfinite(magnitude) or magnitude <= 0.0:
        raise ValueError("trial_magnitude must be positive and finite")
    if morphology_grid_size < 17 or morphology_grid_size % 2 == 0:
        raise ValueError("morphology_grid_size must be odd and >=17")

    orders = tuple(int(value) for value in quadrature_orders)
    if len(orders) < 2 or any(value < 16 for value in orders) or any(
        orders[index] >= orders[index + 1] for index in range(len(orders) - 1)
    ):
        raise ValueError("quadrature_orders must be strictly increasing and >=16")

    constraints = _constraints()
    nontriviality = constraints["nontriviality"]
    validation = constraints["validation"]
    reference_time = float(nontriviality["reference_time"])
    target = float(nontriviality["reference_energy"])
    validation_energy_lower = float(nontriviality["minimum_energy_each_validation_time"])
    validation_energy_upper = float(nontriviality["maximum_energy_each_validation_time"])
    times = tuple(float(value) for value in validation["times"])

    field = _source_field()
    implementation_guard = float(field.parent.profile_basis.coefficient_limit)
    if magnitude > implementation_guard:
        raise ValueError("trial_magnitude exceeds inherited implementation guard")

    rows_by_order = {
        order: {time: _energy_quadratic(field, time, order) for time in times}
        for order in orders
    }
    finest = rows_by_order[orders[-1]]
    reference = finest[reference_time]
    baseline_morphology = _cap_vorticity_morphology(
        field, time=float(morphology_time), coefficient=0.0, grid_size=int(morphology_grid_size)
    )

    negative = _trial_row(
        field,
        coefficient=-magnitude,
        reference_coefficients=reference,
        energy_by_time=finest,
        target=target,
        validation_energy_lower=validation_energy_lower,
        validation_energy_upper=validation_energy_upper,
        times=times,
        morphology_time=morphology_time,
        morphology_grid_size=morphology_grid_size,
        baseline_morphology=baseline_morphology,
    )
    positive = _trial_row(
        field,
        coefficient=magnitude,
        reference_coefficients=reference,
        energy_by_time=finest,
        target=target,
        validation_energy_lower=validation_energy_lower,
        validation_energy_upper=validation_energy_upper,
        times=times,
        morphology_time=morphology_time,
        morphology_grid_size=morphology_grid_size,
        baseline_morphology=baseline_morphology,
    )

    both_preflight = bool(
        negative["profile_bound_screen"]["simple_common_scale_representation_preflight_passed"]
        and positive["profile_bound_screen"]["simple_common_scale_representation_preflight_passed"]
    )
    both_tip = bool(
        negative["gross_tip_reach_extended_on_this_grid"] and positive["gross_tip_reach_extended_on_this_grid"]
    )

    return {
        "task_id": TASK_ID,
        "source_candidate_sha256": str(field.sha256),
        "mode": AXIAL_CAP_POLOIDAL,
        "basis_growth": {
            "new_basis_shapes_added": 0,
            "new_temporal_degrees_added": 0,
            "new_scalar_fit_parameters_added": 0,
            "coefficient_selected": False,
        },
        "energy_model": {
            "definition": "lambda(a)=sqrt(E_target/E_raw(a)); u_joint=lambda(a)*(u_parent+a*delta_u_cap)",
            "quadrature_orders": list(orders),
            "finest_order": int(orders[-1]),
            "reference_time": reference_time,
            "target_energy": target,
            "reference_energy_quadratic": dict(reference),
            "maximum_relative_quadratic_coefficient_change": _refinement(rows_by_order),
        },
        "diagnostic_trial_magnitude": magnitude,
        "baseline_vorticity_morphology": baseline_morphology,
        "negative_trial": negative,
        "positive_trial": positive,
        "both_signs_simple_representation_preflight_passed": both_preflight,
        "both_signs_gross_tip_reach_extended_on_this_grid": both_tip,
        "routing_contract": (
            "Do not add another axial-cap, higher odd-q, F30/F40, or swirl basis from this audit. If the useful +/-0.25 cap "
            "geometry survives exact common-energy renormalization and the inherited coefficient preflight, hand this single "
            "mode to Agent 1/2 for an explicitly bounded, newly identified child and fresh reference/validation energy, core-sign, "
            "support/divergence, pressure/force, and held-out full-momentum checks. If those candidate-level checks reject every "
            "useful nonzero coefficient, route the remaining tip defect to base support/taper geometry rather than basis growth."
        ),
        "truth_boundary": dict(TRUTH_BOUNDARY),
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--trial-magnitude", type=float, default=DEFAULT_TRIAL_MAGNITUDE)
    parser.add_argument("--quadrature-order", type=int, action="append", dest="orders")
    parser.add_argument("--morphology-grid-size", type=int, default=33)
    args = parser.parse_args(argv)
    orders = tuple(args.orders) if args.orders else (96, 192)
    report = audit_axial_cap_joint_renormalization(
        trial_magnitude=args.trial_magnitude,
        quadrature_orders=orders,
        morphology_grid_size=args.morphology_grid_size,
    )
    text = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output is None:
        print(text, end="")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
