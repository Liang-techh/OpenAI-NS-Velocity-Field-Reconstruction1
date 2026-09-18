"""Fixed-parent energy-envelope audit for the Agent-7 axial-cap poloidal mode.

The preceding capacity screen showed that the fixed support-adjacent cap mode can
move gross q90/q99 vorticity reach at diagnostic coefficient +/-0.25.  That
amplitude was never declared materializable.  This module measures how much of
that morphology survives the unchanged preregistered E(0.25)=1+/-0.001 gate
when the parent normalization is held fixed.

No basis, coefficient, pressure, forcing, residual or public-image target is fit.
This is pre-materialization representation evidence, not PDE validation.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Iterable

import numpy as np
from numpy.polynomial.legendre import leggauss

from .constrained_eq45_bipolar_axial_cap_poloidal_capacity import (
    _axial_cap_basis_velocity,
    _cap_velocity,
    _cap_vorticity_morphology,
    _source_field,
)

TASK_ID = "CR003-BIPOLAR-AXIAL-CAP-ENERGY-ENVELOPE-046"
PREVIOUS_CAPACITY_TRIAL = 0.25
TRUTH_BOUNDARY = {
    "velocity_changed": False,
    "candidate_artifact_changed": False,
    "new_basis_added": False,
    "axial_cap_poloidal_coefficient_selected": False,
    "materialization_bound_selected": False,
    "candidate_sha_created": False,
    "held_out_pde_residual_evaluated": False,
    "force_or_pressure_changed": False,
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


def _grid(order: int) -> tuple[np.ndarray, np.ndarray]:
    order = int(order)
    if order < 16:
        raise ValueError("quadrature order must be >=16")
    nodes, weights = leggauss(order)
    r, z = np.meshgrid(nodes + 1.0, 2.0 * nodes, indexing="ij")
    wr, wz = np.meshgrid(weights, 2.0 * weights, indexing="ij")
    points = np.column_stack((r.ravel(), np.zeros(r.size), z.ravel()))
    volume = (2.0 * np.pi * r * wr * wz).ravel()
    return points, volume


def _energy_quadratic(field, time: float, order: int) -> dict[str, float]:
    points, volume = _grid(order)
    base = np.asarray(field.at_points(points, float(time)), dtype=float)
    basis = np.asarray(_axial_cap_basis_velocity(field, points), dtype=float)
    if base.shape != basis.shape or base.shape != (len(points), 3):
        raise RuntimeError("unexpected velocity shape")
    constant = 0.5 * float(np.sum(volume * np.sum(base * base, axis=1)))
    linear = float(np.sum(volume * np.sum(base * basis, axis=1)))
    quadratic = 0.5 * float(np.sum(volume * np.sum(basis * basis, axis=1)))
    if not np.all(np.isfinite((constant, linear, quadratic))) or quadratic <= 0.0:
        raise RuntimeError("invalid kinetic-energy quadratic")
    return {"constant": constant, "linear": linear, "quadratic": quadratic}


def _energy(coefficients: dict[str, float], coefficient: float) -> float:
    a = float(coefficient)
    return float(coefficients["constant"] + coefficients["linear"] * a + coefficients["quadratic"] * a * a)


def _direct_energy(field, time: float, coefficient: float, order: int) -> float:
    points, volume = _grid(order)
    velocity = np.asarray(_cap_velocity(field, points, float(time), float(coefficient)), dtype=float)
    return 0.5 * float(np.sum(volume * np.sum(velocity * velocity, axis=1)))


def _interval_extrema(coefficients: dict[str, float], width: float) -> tuple[float, float]:
    width = float(width)
    if not np.isfinite(width) or width < 0.0:
        raise ValueError("width must be finite and nonnegative")
    left, right = _energy(coefficients, -width), _energy(coefficients, width)
    minimum = min(left, right)
    vertex = -coefficients["linear"] / (2.0 * coefficients["quadratic"])
    if -width <= vertex <= width:
        minimum = min(minimum, _energy(coefficients, vertex))
    return float(minimum), float(max(left, right))


def _symmetric_envelope(coefficients: dict[str, float], lower: float, upper: float, guard: float) -> dict[str, Any]:
    lower, upper, guard = float(lower), float(upper), float(guard)
    if not lower < upper or not np.isfinite(guard) or guard <= 0.0:
        raise ValueError("invalid envelope inputs")
    baseline = _energy(coefficients, 0.0)
    if not lower <= baseline <= upper:
        return {"exists": False, "half_width": 0.0, "guard_limited": False,
                "minimum_energy": baseline, "maximum_energy": baseline}

    def ok(width: float) -> bool:
        minimum, maximum = _interval_extrema(coefficients, width)
        return minimum >= lower and maximum <= upper

    if ok(guard):
        minimum, maximum = _interval_extrema(coefficients, guard)
        return {"exists": True, "half_width": guard, "guard_limited": True,
                "minimum_energy": minimum, "maximum_energy": maximum}
    lo, hi = 0.0, guard
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if ok(mid):
            lo = mid
        else:
            hi = mid
    minimum, maximum = _interval_extrema(coefficients, lo)
    return {"exists": True, "half_width": float(lo), "guard_limited": False,
            "minimum_energy": minimum, "maximum_energy": maximum}


def _normalization_roots(coefficients: dict[str, float], target: float, guard: float) -> list[float]:
    q = float(coefficients["quadratic"])
    b = float(coefficients["linear"])
    c = float(coefficients["constant"] - target)
    disc = b * b - 4.0 * q * c
    if disc < 0.0:
        return []
    root = float(np.sqrt(max(disc, 0.0)))
    values = [(-b - root) / (2.0 * q), (-b + root) / (2.0 * q)]
    return sorted({float(value) for value in values if np.isfinite(value) and abs(value) <= guard})


def _refinement(rows: dict[int, dict[float, dict[str, float]]]) -> float:
    orders = sorted(rows)
    coarse, fine = rows[orders[-2]], rows[orders[-1]]
    return float(max(
        abs(coarse[time][key] - fine[time][key]) / max(abs(fine[time][key]), 1.0e-12)
        for time in fine for key in ("constant", "linear", "quadratic")
    ))


def _morph_delta(base: dict[str, float], other: dict[str, float]) -> dict[str, float]:
    keys = (
        "axial_rms_over_Zp", "radial_rms_over_Rp", "axial_q90_over_Zp", "axial_q99_over_Zp",
        "outer_065_enstrophy_fraction", "outer_075_enstrophy_fraction",
    )
    return {key: float(other[key] - base[key]) for key in keys}


def audit_axial_cap_energy_envelope(
    *, quadrature_orders: Iterable[int] = (96, 192), morphology_time: float = 0.5,
    morphology_grid_size: int = 33,
) -> dict[str, Any]:
    orders = tuple(int(value) for value in quadrature_orders)
    if len(orders) < 2 or any(value < 16 for value in orders) or any(a >= b for a, b in zip(orders, orders[1:])):
        raise ValueError("quadrature_orders must be strictly increasing and >=16")
    if morphology_grid_size < 17 or morphology_grid_size % 2 == 0:
        raise ValueError("morphology_grid_size must be odd and >=17")

    constraints = _constraints()
    nontriviality = constraints["nontriviality"]
    times = tuple(float(value) for value in constraints["validation"]["times"])
    reference_time = float(nontriviality["reference_time"])
    target = float(nontriviality["reference_energy"])
    tolerance = float(nontriviality["reference_energy_abs_tolerance"])
    validation_lower = float(nontriviality["minimum_energy_each_validation_time"])
    validation_upper = float(nontriviality["maximum_energy_each_validation_time"])

    field = _source_field()
    guard = float(field.parent.profile_basis.coefficient_limit)
    rows = {order: {time: _energy_quadratic(field, time, order) for time in times} for order in orders}
    finest = rows[orders[-1]]
    reference = finest[reference_time]
    envelope = _symmetric_envelope(reference, target - tolerance, target + tolerance, guard)
    half_width = float(envelope["half_width"])
    if not envelope["exists"] or half_width <= 0.0:
        raise RuntimeError("no nonzero fixed-parent energy envelope exists")

    validation_rows = []
    validation_ok = True
    for time in times:
        minimum, maximum = _interval_extrema(finest[time], half_width)
        passed = bool(minimum >= validation_lower and maximum <= validation_upper)
        validation_ok &= passed
        validation_rows.append({"time": time, "baseline_energy": _energy(finest[time], 0.0),
                                "minimum_energy": minimum, "maximum_energy": maximum,
                                "validation_energy_range_satisfied": passed})

    replay = []
    for coefficient in (-half_width, 0.0, half_width, -PREVIOUS_CAPACITY_TRIAL, PREVIOUS_CAPACITY_TRIAL):
        if abs(coefficient) <= guard:
            predicted = _energy(reference, coefficient)
            direct = _direct_energy(field, reference_time, coefficient, orders[-1])
            replay.append({"coefficient": float(coefficient), "quadratic_energy": predicted,
                           "direct_energy": direct, "absolute_difference": abs(predicted - direct)})

    morphology = {
        "baseline": _cap_vorticity_morphology(field, time=morphology_time, coefficient=0.0, grid_size=morphology_grid_size),
        "minus_envelope": _cap_vorticity_morphology(field, time=morphology_time, coefficient=-half_width, grid_size=morphology_grid_size),
        "plus_envelope": _cap_vorticity_morphology(field, time=morphology_time, coefficient=half_width, grid_size=morphology_grid_size),
        "minus_capacity_trial": _cap_vorticity_morphology(field, time=morphology_time, coefficient=-PREVIOUS_CAPACITY_TRIAL, grid_size=morphology_grid_size),
        "plus_capacity_trial": _cap_vorticity_morphology(field, time=morphology_time, coefficient=PREVIOUS_CAPACITY_TRIAL, grid_size=morphology_grid_size),
    }
    baseline = morphology["baseline"]
    morphology["minus_envelope_delta"] = _morph_delta(baseline, morphology["minus_envelope"])
    morphology["plus_envelope_delta"] = _morph_delta(baseline, morphology["plus_envelope"])
    for percentile in (90, 99):
        key = f"axial_q{percentile}_over_Zp"
        base = float(baseline[key])
        envelope_max = max(float(morphology["minus_envelope"][key]), float(morphology["plus_envelope"][key]))
        trial_max = max(float(morphology["minus_capacity_trial"][key]), float(morphology["plus_capacity_trial"][key]))
        morphology[f"baseline_q{percentile}_over_Zp"] = base
        morphology[f"energy_envelope_max_q{percentile}_over_Zp"] = envelope_max
        morphology[f"capacity_trial_max_q{percentile}_over_Zp"] = trial_max
        morphology[f"energy_envelope_q{percentile}_moved_on_this_grid"] = bool(envelope_max != base)

    return {
        "task_id": TASK_ID,
        "source_candidate_sha256": str(field.sha256),
        "basis": "AXIAL_CAP_BANDED_C4_ODD_Z_POLOIDAL",
        "basis_parameter_count": {"fixed_spatial_shapes": 1, "potential_scalar_coefficients": 1,
                                  "fitted_shape_parameters": 0, "temporal_degrees_added": 0},
        "quadrature_orders": list(orders),
        "quadrature_max_relative_coefficient_change": _refinement(rows),
        "reference_energy_quadratic": dict(reference),
        "reference_energy_target": target,
        "reference_energy_abs_tolerance": tolerance,
        "zero_centered_fixed_parent_energy_envelope": dict(envelope),
        "normalization_roots_within_diagnostic_guard": _normalization_roots(reference, target, guard),
        "previous_capacity_trial_abs_coefficient": PREVIOUS_CAPACITY_TRIAL,
        "capacity_trial_to_energy_envelope_amplitude_ratio": float(PREVIOUS_CAPACITY_TRIAL / half_width),
        "capacity_trial_inside_fixed_parent_energy_envelope": bool(PREVIOUS_CAPACITY_TRIAL <= half_width),
        "validation_time_energy_rows": validation_rows,
        "validation_time_energy_range_satisfied_over_envelope": bool(validation_ok),
        "direct_energy_replay_checks": replay,
        "morphology_time": float(morphology_time),
        "morphology_grid_size": int(morphology_grid_size),
        "morphology": morphology,
        "routing_contract": (
            "If gross q90/q99 movement survives inside this fixed-parent energy envelope, Agent 1 may materialize only this "
            "one bounded cap coefficient and Agent 2/3 must rerun fresh candidate-specific energy/sign/support/divergence and "
            "held-out full momentum checks. If gross tip movement disappears, do not add another cap, higher odd-q, or swirl "
            "basis: next test joint energy renormalization for this same mode; if that also fails, route the tip defect to base "
            "support/taper geometry."
        ),
        "truth_boundary": dict(TRUTH_BOUNDARY),
    }


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--morphology-grid-size", type=int, default=33)
    args = parser.parse_args(list(argv) if argv is not None else None)
    report = audit_axial_cap_energy_envelope(morphology_grid_size=args.morphology_grid_size)
    text = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output is None:
        print(text, end="")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding="utf-8")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
