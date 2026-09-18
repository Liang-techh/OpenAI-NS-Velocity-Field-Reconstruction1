"""Joint-energy-renormalized screen for the axial-envelope geometry family.

PR #359 showed that every nonzero small positive ``alpha`` in the fixed-parent
axial-envelope family immediately violates the unchanged reference-energy gate,
while useful q90 leverage appears only near the upper end of that small-alpha
window and is not robust across the registered time slices.  This module asks
one deliberately narrow follow-up question without adding a basis:

    lambda(alpha) = sqrt(E_target / E_raw(alpha, t_ref))
    u_joint       = lambda(alpha) * u_alpha

For nonzero alpha the positive common scale restores the reference kinetic
energy exactly.  Alpha zero is the frozen parent control and is replayed without
rescaling when it already satisfies the registered reference-energy tolerance;
this prevents quadrature roundoff from spuriously pushing an already-bound
parent coefficient past its inherited limit.  Positive common scaling preserves
instantaneous streamline directions and normalized vorticity-location
fingerprints, but it is not a Navier--Stokes invariance.  The screen therefore
checks only representation/energy/core/support preflight plus target-free 3-D
vorticity morphology.  It does not create a production candidate, fit pressure
or forcing, or claim held-out PDE evidence.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Iterable

import numpy as np

from .constrained_eq45_axial_envelope_flattening_capacity import (
    _morphology,
    _source_field,
    _structure_checks,
    _validate_alpha,
)
from .constrained_eq45_axial_envelope_small_alpha_screen import _constraints, _energy

TASK_ID = "CR003-AXIAL-ENVELOPE-JOINT-RENORM-SCREEN-055"
DEFAULT_ALPHAS = (0.0, 0.25, 0.35, 0.50, 0.75, 1.0)
DEFAULT_MORPHOLOGY_TIMES = (0.25, 0.50, 0.75)
MORPHOLOGY_ZERO_TOLERANCE = 1.0e-12
TRUTH_BOUNDARY = {
    "velocity_changed": False,
    "canonical_velocity_changed": False,
    "candidate_artifact_changed": False,
    "new_spatial_basis_added": False,
    "joint_renormalized_child_materialized": False,
    "production_geometry_parameter_selected": False,
    "production_bound_selected": False,
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


def _validate_alpha_grid(alphas: Iterable[float]) -> tuple[float, ...]:
    values = tuple(_validate_alpha(float(value)) for value in alphas)
    if len(values) < 2 or values[0] != 0.0:
        raise ValueError("alpha grid must start at 0 and contain at least two values")
    if any(a >= b for a, b in zip(values, values[1:])):
        raise ValueError("alpha grid must be strictly increasing")
    if values[-1] > 1.0:
        raise ValueError("joint-renormalized screen is restricted to alpha<=1")
    return values


def _normalization_scale(raw_reference_energy: float, target: float) -> float:
    raw_reference_energy = float(raw_reference_energy)
    target = float(target)
    if not np.isfinite(raw_reference_energy) or raw_reference_energy <= 0.0:
        raise ValueError("raw reference energy must be positive and finite")
    if not np.isfinite(target) or target <= 0.0:
        raise ValueError("target energy must be positive and finite")
    return float(np.sqrt(target / raw_reference_energy))


def _profile_bound_screen(field, scale: float, alpha: float) -> dict[str, Any]:
    """Check the simplest common-scale realization against inherited Eq45 limits."""
    basis = field.parent.profile_basis
    limit = float(basis.coefficient_limit)
    stored = np.concatenate(
        (
            np.asarray(basis.phi_coefficients, dtype=float),
            np.asarray(basis.swirl_coefficients, dtype=float),
        )
    )
    scaled = float(scale) * stored
    max_before = float(np.max(np.abs(stored))) if stored.size else 0.0
    max_after = float(np.max(np.abs(scaled))) if stored.size else 0.0
    tolerance = 64.0 * np.finfo(float).eps * max(1.0, limit)
    stored_ok = bool(max_after <= limit + tolerance)
    alpha_ok = bool(0.0 <= float(alpha) <= 1.0)
    return {
        "inherited_coefficient_limit": limit,
        "maximum_abs_stored_profile_coefficient_before_scale": max_before,
        "maximum_abs_stored_profile_coefficient_after_scale": max_after,
        "stored_profile_coefficients_within_inherited_limit": stored_ok,
        "geometry_alpha_within_preregistered_screen_interval": alpha_ok,
        "simple_common_scale_representation_preflight_passed": bool(stored_ok and alpha_ok),
        "scope": (
            "representation preflight only: alpha changes the governed envelope geometry while the positive common scale "
            "multiplies the resulting velocity. A materialized child still needs a new identity and fresh pressure/force/PDE checks."
        ),
    }


def _morphology_delta(base: dict[str, float], trial: dict[str, float]) -> dict[str, float]:
    return {f"{key}_delta": float(trial[key] - base[key]) for key in base}


def audit_axial_envelope_joint_renorm_screen(
    *,
    alphas: Iterable[float] = DEFAULT_ALPHAS,
    quadrature_orders: Iterable[int] = (64, 96),
    morphology_times: Iterable[float] = DEFAULT_MORPHOLOGY_TIMES,
    morphology_grid_size: int = 33,
) -> dict[str, Any]:
    alphas = _validate_alpha_grid(alphas)
    orders = tuple(int(value) for value in quadrature_orders)
    if len(orders) != 2 or not 16 <= orders[0] < orders[1]:
        raise ValueError("quadrature_orders must contain two increasing orders >=16")
    morph_times = tuple(float(value) for value in morphology_times)
    if not morph_times or any(not np.isfinite(value) for value in morph_times):
        raise ValueError("morphology_times must be finite and nonempty")
    if morphology_grid_size < 17 or morphology_grid_size % 2 == 0:
        raise ValueError("morphology_grid_size must be odd and >=17")

    cfg = _constraints()
    nontriviality = cfg["nontriviality"]
    validation_times = tuple(float(value) for value in cfg["validation"]["times"])
    reference_time = float(nontriviality["reference_time"])
    target = float(nontriviality["reference_energy"])
    tolerance = float(nontriviality["reference_energy_abs_tolerance"])
    validation_lower = float(nontriviality["minimum_energy_each_validation_time"])
    validation_upper = float(nontriviality["maximum_energy_each_validation_time"])

    field = _source_field()
    energy_by_order: dict[int, dict[float, dict[float, float]]] = {}
    for order in orders:
        energy_by_order[order] = {
            alpha: {time: _energy(field, alpha=alpha, time=time, order=order) for time in validation_times}
            for alpha in alphas
        }
    coarse = energy_by_order[orders[0]]
    fine = energy_by_order[orders[1]]
    energy_refinement = max(
        abs(coarse[alpha][time] - fine[alpha][time]) / max(abs(fine[alpha][time]), 1.0e-12)
        for alpha in alphas
        for time in validation_times
    )

    baseline_by_time = {
        time: _morphology(field, time=time, alpha=0.0, grid_size=morphology_grid_size)
        for time in morph_times
    }

    rows = []
    robust_rows = []
    for alpha in alphas:
        raw_reference = float(fine[alpha][reference_time])
        baseline_parent_replay = bool(alpha == 0.0 and abs(raw_reference - target) <= tolerance)
        scale = 1.0 if baseline_parent_replay else _normalization_scale(raw_reference, target)
        normalized_reference = scale * scale * raw_reference
        validation_rows = []
        validation_ok = True
        for time in validation_times:
            raw = float(fine[alpha][time])
            normalized = scale * scale * raw
            passed = bool(validation_lower <= normalized <= validation_upper)
            validation_ok &= passed
            validation_rows.append(
                {
                    "time": time,
                    "raw_energy": raw,
                    "jointly_renormalized_energy": normalized,
                    "validation_energy_range_pass": passed,
                }
            )

        morphology_rows = []
        q90_gains = []
        q99_gains = []
        axial_rel = []
        radial_rel = []
        collar_ratio = []
        for time in morph_times:
            base = baseline_by_time[time]
            trial = base if alpha == 0.0 else _morphology(
                field, time=time, alpha=alpha, grid_size=morphology_grid_size
            )
            q90_gain = float(trial["axial_q90_over_Zp"] - base["axial_q90_over_Zp"])
            q99_gain = float(trial["axial_q99_over_Zp"] - base["axial_q99_over_Zp"])
            q90_gains.append(q90_gain)
            q99_gains.append(q99_gain)
            axial_rel.append(float(trial["axial_rms_over_Zp"] / base["axial_rms_over_Zp"] - 1.0))
            radial_rel.append(float(trial["radial_rms_over_Rp"] / base["radial_rms_over_Rp"] - 1.0))
            collar_ratio.append(
                float(
                    trial["physical_plateau_collar_enstrophy_fraction"]
                    / max(base["physical_plateau_collar_enstrophy_fraction"], 1.0e-300)
                )
            )
            morphology_rows.append(
                {
                    "time": time,
                    "metrics": trial,
                    "delta_from_alpha_zero": _morphology_delta(base, trial),
                }
            )

        profile = _profile_bound_screen(field, scale, alpha)
        reference_ok = bool(abs(normalized_reference - target) <= tolerance)
        robust_q90 = bool(all(value > MORPHOLOGY_ZERO_TOLERANCE for value in q90_gains))
        robust_q99 = bool(all(value > MORPHOLOGY_ZERO_TOLERANCE for value in q99_gains))
        preflight = bool(reference_ok and validation_ok and profile["simple_common_scale_representation_preflight_passed"])
        row = {
            "alpha": alpha,
            "baseline_parent_replay_without_common_scale": baseline_parent_replay,
            "raw_reference_energy": raw_reference,
            "common_velocity_scale": scale,
            "common_scale_fractional_change": float(scale - 1.0),
            "jointly_renormalized_reference_energy": float(normalized_reference),
            "reference_energy_gate_pass": reference_ok,
            "validation_energy_range_all_pass": bool(validation_ok),
            "validation_time_energy": validation_rows,
            "profile_bound_screen": profile,
            "representation_energy_preflight_passed": preflight,
            "morphology_by_time": morphology_rows,
            "minimum_q90_gain_over_Zp_across_morphology_times": float(min(q90_gains)),
            "minimum_q99_gain_over_Zp_across_morphology_times": float(min(q99_gains)),
            "q90_gain_over_Zp_by_time": q90_gains,
            "q99_gain_over_Zp_by_time": q99_gains,
            "q90_gain_positive_at_all_morphology_times": robust_q90,
            "q99_gain_positive_at_all_morphology_times": robust_q99,
            "minimum_axial_rms_relative_change_across_morphology_times": float(min(axial_rel)),
            "maximum_axial_rms_relative_change_across_morphology_times": float(max(axial_rel)),
            "maximum_radial_rms_relative_change_across_morphology_times": float(max(radial_rel)),
            "maximum_collar_enstrophy_fraction_ratio_across_morphology_times": float(max(collar_ratio)),
            "positive_common_scale_morphology_invariance": True,
            "navier_stokes_balance_invariant_under_common_scale": False,
        }
        if alpha > 0.0 and preflight and robust_q90:
            robust_rows.append(row)
        rows.append(row)

    smallest_robust = robust_rows[0]["alpha"] if robust_rows else None
    structure_at_max = _structure_checks(field, trial_alpha=alphas[-1])

    return {
        "task_id": TASK_ID,
        "source_candidate_sha256": str(field.sha256),
        "parent_screen_task": "CR003-AXIAL-ENVELOPE-SMALL-ALPHA-SCREEN-054",
        "screen_protocol": {
            "alphas": list(alphas),
            "quadrature_orders": list(orders),
            "morphology_times": list(morph_times),
            "morphology_grid_size": int(morphology_grid_size),
            "parameter_count_added_for_screen": 0,
            "new_basis_shapes_added": 0,
            "production_alpha_selected": False,
        },
        "energy_quadrature_max_relative_refinement_change": float(energy_refinement),
        "baseline_morphology_by_time": baseline_by_time,
        "rows": rows,
        "smallest_sampled_preflight_passing_multi_time_q90_mover": smallest_robust,
        "smallest_mover_is_production_selection": False,
        "max_alpha_structure_checks_before_common_scale": structure_at_max,
        "held_out_pde_residual": {
            "evaluated": False,
            "reason": (
                "The envelope proxy remains velocity-only and has no governed pressure/force child. Common velocity scaling "
                "is not an NS invariance, so candidate-level pressure/force and held-out momentum must be rebuilt after materialization."
            ),
        },
        "routing_contract": (
            "Do not add another swirl/cap/odd-q basis. If a moderate alpha passes the inherited representation/energy preflight "
            "and moves q90 at every sampled time without unacceptable radial/collar cost, hand at most this one envelope-geometry "
            "degree to Agent 1 for explicit bounded child materialization; Agent 2/3 must then rebuild compatible pressure/force "
            "and run fresh divergence/full-momentum checks plus actual 3-D streamline/vorticity comparison. If robust q90 motion "
            "appears only at a strongly distorted/collar-heavy alpha, stop this envelope route and escalate to larger support/coordinate geometry."
        ),
        "truth_boundary": dict(TRUTH_BOUNDARY),
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--morphology-grid-size", type=int, default=33)
    args = parser.parse_args(argv)
    report = audit_axial_envelope_joint_renorm_screen(
        morphology_grid_size=args.morphology_grid_size
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
