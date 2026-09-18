"""Target-free joint-energy-renormalization audit for the compact poloidal mode.

This increment adds no basis and does not materialize a candidate.  PR #309
showed that the already-integrated ``COMPACT_C4_ODD_Z_POLOIDAL`` direction has
very little usable amplitude when the parent normalization is held fixed.  Here
we ask one narrower representation question: if the *entire* trial velocity is
multiplied by one positive scalar so that the reference kinetic energy is
restored exactly, does the previous finite-amplitude morphology become
representable without violating the inherited profile coefficient limit?

The common scale is

    lambda(a) = sqrt(E_target / E_raw(a)),

where ``E_raw(a)`` is the independently quadrature-evaluated quadratic from the
parent audit.  Positive common scaling preserves support, parity, instantaneous
streamline directions, and all normalized vorticity-location fingerprints; it
does not preserve Navier--Stokes momentum balance.  Pressure and forcing are not
changed or fitted here.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Iterable

import numpy as np

from .constrained_eq45_bipolar_compact_poloidal_capacity import (
    COMPACT_POLOIDAL,
    _compact_velocity,
    _source_field,
    _vorticity_morphology,
)
from .constrained_eq45_bipolar_compact_poloidal_energy_envelope import (
    _constraints,
    _energy_quadratic,
    _quadrature_refinement,
    _quadratic_value,
)

TASK_ID = "CR003-BIPOLAR-COMPACT-POLOIDAL-JOINT-RENORM-CAPACITY-045"
DEFAULT_TRIAL_MAGNITUDE = 0.25
TRUTH_BOUNDARY = {
    "velocity_changed": False,
    "candidate_artifact_changed": False,
    "new_basis_added": False,
    "joint_renormalized_child_materialized": False,
    "compact_poloidal_coefficient_selected": False,
    "materialization_bound_selected": False,
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
    raw = _quadratic_value(coefficients, float(coefficient))
    if not np.isfinite(raw) or raw <= 0.0 or not np.isfinite(target) or target <= 0.0:
        raise ValueError("energy normalization requires positive finite energies")
    return float(np.sqrt(float(target) / raw))


def _profile_bound_screen(field, scale: float, compact_coefficient: float) -> dict[str, Any]:
    """Check the simplest coefficient-space realization of a common velocity scale.

    The current Eq45 velocity is linear in the stored Phi/F coefficients.  This
    screen therefore asks whether multiplying those stored coefficients by the
    common positive scale would stay inside the inherited coefficient limit.
    It is only a representation preflight: a real child still needs a new
    identity and full validation.
    """
    basis = field.parent.profile_basis
    limit = float(basis.coefficient_limit)
    phi = np.asarray(basis.phi_coefficients, dtype=float)
    swirl = np.asarray(basis.swirl_coefficients, dtype=float)
    stored = np.concatenate((phi, swirl))
    scaled = float(scale) * stored
    effective_compact = float(scale) * float(compact_coefficient)
    max_before = float(np.max(np.abs(stored))) if stored.size else 0.0
    max_after = float(np.max(np.abs(scaled))) if scaled.size else 0.0
    tolerance = 64.0 * np.finfo(float).eps * max(1.0, limit)
    stored_ok = bool(max_after <= limit + tolerance)
    compact_guard_ok = bool(abs(effective_compact) <= limit + tolerance)
    return {
        "inherited_coefficient_limit": limit,
        "maximum_abs_stored_profile_coefficient_before_scale": max_before,
        "maximum_abs_stored_profile_coefficient_after_scale": max_after,
        "effective_compact_coefficient_after_common_scale": effective_compact,
        "stored_profile_coefficients_within_inherited_limit": stored_ok,
        "effective_compact_coefficient_within_inherited_implementation_guard": compact_guard_ok,
        "simple_common_scale_representation_preflight_passed": bool(stored_ok and compact_guard_ok),
        "scope": (
            "coefficient-space preflight only; it does not select a production bound or prove that a "
            "jointly renormalized child satisfies momentum/pressure/forcing constraints"
        ),
    }


def _core_sign_screen(field, coefficient: float, scale: float, times: Iterable[float]) -> dict[str, Any]:
    rows = []
    for time in tuple(float(v) for v in times):
        tau = 1.0 - time
        point = np.array([[0.1 * np.sqrt(tau), 0.0, 0.1 * tau**0.495]], dtype=float)
        raw = np.asarray(_compact_velocity(field, point, time, coefficient)[0], dtype=float)
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
        "rows": rows,
        "positive_common_scale_preserves_signs_analytically": bool(scale > 0.0),
    }


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
) -> dict[str, Any]:
    coefficient = float(coefficient)
    scale = _normalization_scale(reference_coefficients, coefficient, target)
    raw_reference_energy = _quadratic_value(reference_coefficients, coefficient)
    validation_rows = []
    all_validation_energy = True
    for time in times:
        raw = _quadratic_value(energy_by_time[time], coefficient)
        normalized = scale * scale * raw
        passed = bool(validation_energy_lower <= normalized <= validation_energy_upper)
        all_validation_energy &= passed
        validation_rows.append({
            "time": time,
            "raw_energy": float(raw),
            "jointly_renormalized_energy": float(normalized),
            "validation_energy_range_satisfied": passed,
        })

    morphology = _vorticity_morphology(
        field,
        time=float(morphology_time),
        mode=COMPACT_POLOIDAL,
        coefficient=coefficient,
        grid_size=int(morphology_grid_size),
    )
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
        "morphology_invariance_under_positive_common_scale": {
            "instantaneous_streamline_geometry_invariant": True,
            "normalized_vorticity_location_fingerprints_invariant": True,
            "reason": (
                "u is multiplied by one positive scalar at fixed time, so streamline directions are unchanged; "
                "curl(u) and enstrophy weights gain common scalar factors that cancel from normalized location metrics"
            ),
            "time_evolution_speed_not_claimed_invariant": True,
        },
    }


def audit_compact_poloidal_joint_renormalization(
    *,
    trial_magnitude: float = DEFAULT_TRIAL_MAGNITUDE,
    quadrature_orders: Iterable[int] = (96, 192),
    morphology_time: float = 0.5,
    morphology_grid_size: int = 25,
) -> dict[str, Any]:
    magnitude = float(trial_magnitude)
    if not np.isfinite(magnitude) or magnitude <= 0.0:
        raise ValueError("trial_magnitude must be positive and finite")

    constraints = _constraints()
    nontriviality = constraints["nontriviality"]
    validation = constraints["validation"]
    reference_time = float(nontriviality["reference_time"])
    target = float(nontriviality["reference_energy"])
    validation_energy_lower = float(nontriviality["minimum_energy_each_validation_time"])
    validation_energy_upper = float(nontriviality["maximum_energy_each_validation_time"])
    times = tuple(float(v) for v in validation["times"])
    orders = tuple(int(v) for v in quadrature_orders)
    if len(orders) < 2 or any(order < 16 for order in orders) or any(
        orders[index] >= orders[index + 1] for index in range(len(orders) - 1)
    ):
        raise ValueError("quadrature_orders must be strictly increasing and >=16")

    field = _source_field()
    rows_by_order = {
        order: {time: _energy_quadratic(field, time, order) for time in times}
        for order in orders
    }
    finest = rows_by_order[orders[-1]]
    reference = finest[reference_time]
    baseline_morphology = _vorticity_morphology(
        field,
        time=float(morphology_time),
        mode=None,
        coefficient=0.0,
        grid_size=int(morphology_grid_size),
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
    )

    def morphology_delta(row: dict[str, Any]) -> dict[str, float]:
        morph = row["vorticity_morphology"]
        return {
            "axial_rms_over_Zp_delta": float(morph["axial_rms_over_Zp"] - baseline_morphology["axial_rms_over_Zp"]),
            "radial_rms_over_Rp_delta": float(morph["radial_rms_over_Rp"] - baseline_morphology["radial_rms_over_Rp"]),
            "axial_q90_over_Zp_delta": float(morph["axial_q90_over_Zp"] - baseline_morphology["axial_q90_over_Zp"]),
            "axial_q99_over_Zp_delta": float(morph["axial_q99_over_Zp"] - baseline_morphology["axial_q99_over_Zp"]),
            "outer_axial_enstrophy_fraction_delta": float(
                morph["outer_axial_enstrophy_fraction"] - baseline_morphology["outer_axial_enstrophy_fraction"]
            ),
        }

    negative["morphology_delta_from_baseline"] = morphology_delta(negative)
    positive["morphology_delta_from_baseline"] = morphology_delta(positive)

    return {
        "task_id": TASK_ID,
        "source_candidate_sha256": field.sha256,
        "mode": COMPACT_POLOIDAL,
        "basis_growth": {
            "new_basis_shapes_added": 0,
            "new_temporal_degrees_added": 0,
            "new_scalar_fit_parameters_added": 0,
            "coefficient_selected": False,
        },
        "energy_model": {
            "definition": "lambda(a)=sqrt(E_target/E_raw(a)); renormalized velocity=lambda(a)*(u+a*delta_u)",
            "quadrature_orders": list(orders),
            "finest_order": orders[-1],
            "reference_time": reference_time,
            "target_energy": target,
            "reference_energy_quadratic": reference,
            "maximum_relative_quadratic_coefficient_change": _quadrature_refinement(rows_by_order),
        },
        "diagnostic_trial_magnitude": magnitude,
        "baseline_vorticity_morphology": baseline_morphology,
        "negative_trial": negative,
        "positive_trial": positive,
        "routing_contract": (
            "Do not add another center compact-poloidal basis from this audit.  If a finite-amplitude compact child is "
            "materialized, the simplest exact common-energy renormalization must respect the inherited coefficient bounds, "
            "receive a new representation identity, and undergo fresh pressure/force, support/divergence and held-out full-momentum "
            "validation.  A common positive scale preserves instantaneous geometry but is not a Navier--Stokes invariance."
        ),
        "truth_boundary": dict(TRUTH_BOUNDARY),
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--trial-magnitude", type=float, default=DEFAULT_TRIAL_MAGNITUDE)
    parser.add_argument("--quadrature-order", type=int, action="append", dest="orders")
    parser.add_argument("--morphology-grid-size", type=int, default=25)
    args = parser.parse_args(argv)
    orders = tuple(args.orders) if args.orders else (96, 192)
    report = audit_compact_poloidal_joint_renormalization(
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
