"""Bounded small-alpha screen for the fixed-support axial-envelope geometry.

The preceding capacity audit showed that replacing the current quadratic axial
compact envelope by a quartic-flattened interior can move q90/q99 vorticity
reach, but alpha=1 also loads the support collar strongly and was never a
materialization proposal.  This module makes the next deliberately small
increment: it scans a preregistered low-alpha line without adding another
basis, fitting a public image, or changing forcing/pressure.

The screen measures reference/validation-time kinetic energy and the same 3-D
vorticity morphology fingerprints used by the capacity audit.  It reports the
smallest *sampled* alpha that moves q90 on the fixed grid while satisfying the
existing energy gates.  That is a diagnostic crossing only, not a production
parameter selection.

A fresh held-out full Navier--Stokes residual is intentionally not fabricated:
the current supported Eq45 artifact is a velocity-only delivery object and this
geometry proxy has no governed pressure/force child.  Candidate-level PDE
screening must therefore follow materialization by the pressure/force lane.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Iterable

import numpy as np
from numpy.polynomial.legendre import leggauss

from .constrained_eq45_axial_envelope_flattening_capacity import (
    _EnvelopeFieldProxy,
    _morphology,
    _source_field,
    _structure_checks,
    _validate_alpha,
)

TASK_ID = "CR003-AXIAL-ENVELOPE-SMALL-ALPHA-SCREEN-054"
DEFAULT_ALPHAS = (0.0, 0.025, 0.05, 0.10, 0.15, 0.20, 0.25)
DEFAULT_MORPHOLOGY_TIMES = (0.25, 0.50, 0.75)
MORPHOLOGY_ZERO_TOLERANCE = 1.0e-12
TRUTH_BOUNDARY = {
    "velocity_changed": False,
    "canonical_velocity_changed": False,
    "candidate_artifact_changed": False,
    "new_spatial_basis_added": False,
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


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _constraints() -> dict[str, Any]:
    value = json.loads((_repo_root() / "configs" / "constraints.json").read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError("constraints.json must contain an object")
    return value


def _axisymmetric_grid(order: int) -> tuple[np.ndarray, np.ndarray]:
    """Gauss--Legendre grid for 2*pi*int_0^2 int_-2^2 (...) r dz dr."""
    order = int(order)
    if order < 16:
        raise ValueError("quadrature order must be >=16")
    nodes, weights = leggauss(order)
    r, z = np.meshgrid(nodes + 1.0, 2.0 * nodes, indexing="ij")
    wr, wz = np.meshgrid(weights, 2.0 * weights, indexing="ij")
    points = np.column_stack((r.ravel(), np.zeros(r.size), z.ravel()))
    volume = (2.0 * np.pi * r * wr * wz).ravel()
    return points, volume


def _energy(field, *, alpha: float, time: float, order: int) -> float:
    points, volume = _axisymmetric_grid(order)
    velocity = np.asarray(_EnvelopeFieldProxy(field, alpha).at_points(points, float(time)), dtype=float)
    if velocity.shape != (len(points), 3) or not np.all(np.isfinite(velocity)):
        raise RuntimeError("unexpected/nonfinite velocity in energy quadrature")
    return 0.5 * float(np.sum(volume * np.sum(velocity * velocity, axis=1)))


def _morphology_delta(base: dict[str, float], trial: dict[str, float]) -> dict[str, float]:
    return {f"{key}_delta": float(trial[key] - base[key]) for key in base}


def _validate_alpha_grid(alphas: Iterable[float]) -> tuple[float, ...]:
    values = tuple(_validate_alpha(float(value)) for value in alphas)
    if len(values) < 2 or values[0] != 0.0:
        raise ValueError("alpha grid must start at 0 and contain at least two values")
    if any(a >= b for a, b in zip(values, values[1:])):
        raise ValueError("alpha grid must be strictly increasing")
    if values[-1] > 0.25:
        raise ValueError("this bounded screen is restricted to alpha<=0.25")
    return values


def audit_axial_envelope_small_alpha_screen(
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
    energy_rows_by_order: dict[int, dict[float, dict[float, float]]] = {}
    for order in orders:
        energy_rows_by_order[order] = {
            alpha: {time: _energy(field, alpha=alpha, time=time, order=order) for time in validation_times}
            for alpha in alphas
        }
    coarse = energy_rows_by_order[orders[0]]
    fine = energy_rows_by_order[orders[1]]
    energy_refinement = max(
        abs(coarse[alpha][time] - fine[alpha][time]) / max(abs(fine[alpha][time]), 1.0e-12)
        for alpha in alphas for time in validation_times
    )

    energy_rows = []
    energy_pass_by_alpha: dict[float, bool] = {}
    for alpha in alphas:
        reference_energy = float(fine[alpha][reference_time])
        validation_values = [float(fine[alpha][time]) for time in validation_times]
        reference_pass = bool(abs(reference_energy - target) <= tolerance)
        validation_pass = bool(
            min(validation_values) >= validation_lower and max(validation_values) <= validation_upper
        )
        energy_pass_by_alpha[alpha] = reference_pass and validation_pass
        energy_rows.append(
            {
                "alpha": alpha,
                "reference_energy": reference_energy,
                "reference_energy_gate_pass": reference_pass,
                "validation_energy_min": min(validation_values),
                "validation_energy_max": max(validation_values),
                "validation_energy_range_pass": validation_pass,
                "combined_energy_gate_pass": bool(reference_pass and validation_pass),
                "validation_time_energies": [
                    {"time": time, "energy": float(fine[alpha][time])} for time in validation_times
                ],
            }
        )

    baseline_by_time = {
        time: _morphology(field, time=time, alpha=0.0, grid_size=morphology_grid_size)
        for time in morph_times
    }
    morphology_rows = []
    morphology_by_alpha_time: dict[float, dict[float, dict[str, float]]] = {}
    for alpha in alphas:
        morphology_by_alpha_time[alpha] = {}
        for time in morph_times:
            trial = (
                baseline_by_time[time]
                if alpha == 0.0
                else _morphology(field, time=time, alpha=alpha, grid_size=morphology_grid_size)
            )
            morphology_by_alpha_time[alpha][time] = trial
            baseline = baseline_by_time[time]
            morphology_rows.append(
                {
                    "alpha": alpha,
                    "time": time,
                    "metrics": trial,
                    "delta_from_alpha_zero": _morphology_delta(baseline, trial),
                }
            )

    diagnostic_time = min(morph_times, key=lambda value: abs(value - 0.5))
    baseline = baseline_by_time[diagnostic_time]
    q90_base = float(baseline["axial_q90_over_Zp"])
    q99_base = float(baseline["axial_q99_over_Zp"])
    axial_base = float(baseline["axial_rms_over_Zp"])
    radial_base = float(baseline["radial_rms_over_Rp"])
    collar_base = float(baseline["physical_plateau_collar_enstrophy_fraction"])

    diagnostic_rows = []
    admissible_q90_movers = []
    for alpha in alphas:
        metrics = morphology_by_alpha_time[alpha][diagnostic_time]
        q90_gain = float(metrics["axial_q90_over_Zp"] - q90_base)
        q99_gain = float(metrics["axial_q99_over_Zp"] - q99_base)
        row = {
            "alpha": alpha,
            "energy_admissible_on_registered_gates": bool(energy_pass_by_alpha[alpha]),
            "q90_gain_over_Zp": q90_gain,
            "q99_gain_over_Zp": q99_gain,
            "q90_moved_on_this_grid": bool(q90_gain > MORPHOLOGY_ZERO_TOLERANCE),
            "q99_moved_on_this_grid": bool(q99_gain > MORPHOLOGY_ZERO_TOLERANCE),
            "axial_rms_relative_change": float(metrics["axial_rms_over_Zp"] / axial_base - 1.0),
            "radial_rms_relative_change": float(metrics["radial_rms_over_Rp"] / radial_base - 1.0),
            "collar_enstrophy_fraction_absolute_change": float(
                metrics["physical_plateau_collar_enstrophy_fraction"] - collar_base
            ),
            "collar_enstrophy_fraction_ratio": float(
                metrics["physical_plateau_collar_enstrophy_fraction"] / max(collar_base, 1.0e-300)
            ),
        }
        if row["energy_admissible_on_registered_gates"] and row["q90_moved_on_this_grid"]:
            admissible_q90_movers.append(alpha)
        diagnostic_rows.append(row)

    smallest_mover = min(admissible_q90_movers) if admissible_q90_movers else None
    max_alpha = alphas[-1]
    structure_at_max = _structure_checks(field, trial_alpha=max_alpha)

    return {
        "task_id": TASK_ID,
        "source_candidate_sha256": str(field.sha256),
        "parent_capacity_task": "CR003-AXIAL-ENVELOPE-FLATTENING-CAPACITY-053",
        "screen_protocol": {
            "alphas": list(alphas),
            "quadrature_orders": list(orders),
            "morphology_times": list(morph_times),
            "morphology_grid_size": int(morphology_grid_size),
            "diagnostic_time_for_crossing": float(diagnostic_time),
            "parameter_count_added_for_screen": 0,
            "new_basis_shapes_added": 0,
        },
        "energy_quadrature_max_relative_refinement_change": float(energy_refinement),
        "energy_rows": energy_rows,
        "morphology_rows": morphology_rows,
        "diagnostic_time_rows": diagnostic_rows,
        "smallest_sampled_energy_admissible_q90_mover": smallest_mover,
        "smallest_mover_is_production_selection": False,
        "max_alpha_structure_checks": structure_at_max,
        "held_out_pde_residual": {
            "evaluated": False,
            "reason": (
                "The supported Eq45 artifact and envelope proxy provide velocity only; no governed pressure/force child "
                "exists for this geometry deformation. A residual with invented or freely fitted pressure/force would violate "
                "the project validation contract."
            ),
        },
        "routing_contract": (
            "Do not add another swirl/cap/odd-q basis from this screen. If a small sampled alpha moves q90 while keeping "
            "the inherited energy and structure gates, Agent 1 may materialize at most one explicitly bounded envelope-shape "
            "child near the smallest useful region, after which Agent 2/3 must construct/govern compatible pressure/force and "
            "run fresh held-out full-momentum validation plus actual 3-D streamline/vorticity comparison. If no bounded small "
            "alpha moves q90 economically, stop this envelope direction and escalate to larger-scale support/coordinate geometry."
        ),
        "truth_boundary": dict(TRUTH_BOUNDARY),
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--morphology-grid-size", type=int, default=33)
    args = parser.parse_args(argv)
    report = audit_axial_envelope_small_alpha_screen(morphology_grid_size=args.morphology_grid_size)
    text = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output is None:
        print(text, end="")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
