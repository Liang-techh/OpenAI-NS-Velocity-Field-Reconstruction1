"""Energy-admissible morphology audit for the fixed axial-cap poloidal mode.

The support-adjacent axial-cap mode from the preceding Agent-7 capacity screen
can move gross q90/q99 vorticity reach at a diagnostic coefficient of +/-0.25.
That coefficient was deliberately *not* declared materializable.  This module
asks the next narrower representation question: how much of that morphology is
available while the parent field and the preregistered reference-energy gate
are held fixed?

No basis shape, pressure, forcing, residual, coefficient or acceptance threshold
is fitted here.  In particular this is not a candidate-selection or PDE screen.
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
    data = json.loads((_repo_root() / "configs" / "constraints.json").read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise ValueError("constraints.json must contain an object")
    return data


def _axisymmetric_grid(order: int) -> tuple[np.ndarray, np.ndarray]:
    order = int(order)
    if order < 16:
        raise ValueError("quadrature order must be >=16")
    nodes, weights = leggauss(order)
    r, z = np.meshgrid(nodes + 1.0, 2.0 * nodes, indexing="ij")
    wr, wz = np.meshgrid(weights, 2.0 * weights, indexing="ij")
    volume = 2.0 * np.pi * r * wr * wz
    points = np.column_stack((r.ravel(), np.zeros(r.size), z.ravel()))
    return points, volume.ravel()


def _energy_quadratic(field, time: float, order: int) -> dict[str, float]:
    """Coefficients of E(a)=constant+linear*a+quadratic*a^2."""
    points, volume = _axisymmetric_grid(order)
    base = np.asarray(field.at_points(points, float(time)), dtype=float)
    basis = np.asarray(_axial_cap_basis_velocity(field, points), dtype=float)
    if base.shape != basis.shape or base.shape != (len(points), 3):
        raise RuntimeError("unexpected velocity shape in axial-cap energy audit")
    constant = 0.5 * float(np.sum(volume * np.sum(base * base, axis=1)))
    linear = float(np.sum(volume * np.sum(base * basis, axis=1)))
    quadratic = 0.5 * float(np.sum(volume * np.sum(basis * basis, axis=1)))
    values = (constant, linear, quadratic)
    if not np.all(np.isfinite(values)) or quadratic <= 0.0:
        raise RuntimeError("invalid kinetic-energy quadratic")
    return {"constant": constant, "linear": linear, "quadratic": quadratic}


def _quadratic_value(coefficients: dict[str, float], coefficient: float) -> float:
    a = float(coefficient)
    return float(
        coefficients["constant"]
        + coefficients["linear"] * a
        + coefficients["quadratic"] * a * a
    )


def _direct_energy(field, time: float, coefficient: float, order: int) -> float:
    points, volume = _axisymmetric_grid(order)
    velocity = np.asarray(_cap_velocity(field, points, float(time), float(coefficient)), dtype=float)
    return 0.5 * float(np.sum(volume * np.sum(velocity * velocity, axis=1)))


def _interval_extrema(coefficients: dict[str, float], half_width: float) -> tuple[float, float]:
    width = float(half_width)
    if not np.isfinite(width) or width < 0.0:
        raise ValueError("half_width must be finite and nonnegative")
    q = float(coefficients["quadratic"])
    b = float(coefficients["linear"])
    left = _quadratic_value(coefficients, -width)
    right = _quadratic_value(coefficients, width)
    minimum = min(left, right)
    vertex = -b / (2.0 * q)
    if -width <= vertex <= width:
        minimum = min(minimum, _quadratic_value(coefficients, vertex))
    return float(minimum), float(max(left, right))


def _largest_symmetric_half_width(
    coefficients: dict[str, float], *, lower: float, upper: float, diagnostic_search_guard: float
) -> dict[str, Any]:
    lower, upper = float(lower), float(upper)
    guard = float(diagnostic_search_guard)
    if not lower < upper or not np.isfinite(guard) or guard <= 0.0:
        raise ValueError("invalid envelope inputs")
    zero = _quadratic_value(coefficients, 0.0)
    if not lower <= zero <= upper:
        return {
            "exists": False,
            "half_width": 0.0,
            "guard_limited": False,
            "minimum_energy": zero,
            "maximum_energy": zero,
        }

    def admissible(width: float) -> bool:
        minimum, maximum = _interval_extrema(coefficients, width)
        return minimum >= lower and maximum <= upper

    if admissible(guard):
        minimum, maximum = _interval_extrema(coefficients, guard)
        return {
            "exists": True,
            "half_width": guard,
            "guard_limited": True,
            "minimum_energy": minimum,
            "maximum_energy": maximum,
        }

    lo, hi = 0.0, guard
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if admissible(mid):
            lo = mid
        else:
            hi = mid
    minimum, maximum = _interval_extrema(coefficients, lo)
    return {
        "exists": True,
        "half_width": float(lo),
        "guard_limited": False,
        "minimum_energy": minimum,
        "maximum_energy": maximum,
    }


def _normalization_roots(
    coefficients: dict[str, float], *, target: float, diagnostic_search_guard: float
) -> list[float]:
    q = float(coefficients["quadratic"])
    b = float(coefficients["linear"])
    c = float(coefficients["constant"] - float(target))
    discriminant = b * b - 4.0 * q * c
    if discriminant < 0.0:
        return []
    root_disc = float(np.sqrt(max(discriminant, 0.0)))
    guard = float(diagnostic_search_guard)
    roots = [(-b - root_disc) / (2.0 * q), (-b + root_disc) / (2.0 * q)]
    return sorted({float(root) for root in roots if np.isfinite(root) and abs(root) <= guard})


def _quadrature_refinement(rows: dict[int, dict[float, dict[str, float]]]) -> float:
    orders = sorted(rows)
    coarse, fine = rows[orders[-2]], rows[orders[-1]]
    maximum = 0.0
    for time in fine:
        for key in ("constant", "linear", "quadratic"):
            denominator = max(abs(fine[time][key]), 1.0e-12)
            maximum = max(maximum, abs(coarse[time][key] - fine[time][key]) / denominator)
    return float(maximum)


def _morphology_delta(base: dict[str, float], other: dict[str, float]) -> dict[str, float]:
    keys = (
        "axial_rms_over_Zp",
        "radial_rms_over_Rp",
        "axial_q90_over_Zp",
        "axial_q99_over_Zp",
        "outer_065_enstrophy_fraction",
        "outer_075_enstrophy_fraction",
    )
    return {key: float(other[key] - base[key]) for key in keys}


def audit_axial_cap_energy_envelope(
    *,
    quadrature_orders: Iterable[int] = (96, 192),
    morphology_time: float = 0.5,
    morphology_grid_size: int = 33,
) -> dict[str, Any]:
    orders = tuple(int(value) for value in quadrature_orders)
    if len(orders) < 2 or any(order < 16 for order in orders):
        raise ValueError("quadrature_orders must contain at least two values >=16")
    if any(orders[index] >= orders[index + 1] for index in range(len(orders) - 1)):
        raise ValueError("quadrature_orders must be strictly increasing")
    if morphology_grid_size < 17 or morphology_grid_size % 2 == 0:
        raise ValueError("morphology_grid_size must be odd and >=17")

    constraints = _constraints()
    nontriviality = constraints["nontriviality"]
    validation = constraints["validation"]
    reference_time = float(nontriviality["reference_time"])
    target = float(nontriviality["reference_energy"])
    tolerance = float(nontriviality["reference_energy_abs_tolerance"])
    validation_energy_lower = float(nontriviality["minimum_energy_each_validation_time"])
    validation_energy_upper = float(nontriviality["maximum_energy_each_validation_time"])
    times = tuple(float(value) for value in validation["times"])

    field = _source_field()
    search_guard = float(field.parent.profile_basis.coefficient_limit)
    rows_by_order = {
        order: {time: _energy_quadratic(field, time, order) for time in times}
        for order in orders
    }
    finest = rows_by_order[orders[-1]]
    reference = finest[reference_time]
    envelope = _largest_symmetric_half_width(
        reference,
        lower=target - tolerance,
        upper=target + tolerance,
        diagnostic_search_guard=search_guard,
    )
    half_width = float(envelope["half_width"])
    if not envelope["exists"] or half_width <= 0.0:
        raise RuntimeError("no nonzero zero-centered axial-cap energy envelope exists")

    validation_rows = []
    validation_range_satisfied = True
    for time in times:
        minimum, maximum = _interval_extrema(finest[time], half_width)
        satisfied = bool(minimum >= validation_energy_lower and maximum <= validation_energy_upper)
        validation_range_satisfied &= satisfied
        validation_rows.append(
            {
                "time": time,
                "baseline_energy": _quadratic_value(finest[time], 0.0),
                "minimum_energy_over_zero_centered_envelope": minimum,
                "maximum_energy_over_zero_centered_envelope": maximum,
                "validation_energy_range_satisfied": satisfied,
            }
        )

    direct_checks = []
    for coefficient in (-half_width, 0.0, half_width, -PREVIOUS_CAPACITY_TRIAL, PREVIOUS_CAPACITY_TRIAL):
        if abs(coefficient) > search_guard:
            continue
        quadratic_energy = _quadratic_value(reference, coefficient)
        direct_energy = _direct_energy(field, reference_time, coefficient, orders[-1])
        direct_checks.append(
            {
                "coefficient": float(coefficient),
                "quadratic_energy": quadratic_energy,
                "direct_energy": direct_energy,
                "absolute_difference": abs(quadratic_energy - direct_energy),
            }
        )

    baseline_morphology = _cap_vorticity_morphology(
        field, time=float(morphology_time), coefficient=0.0, grid_size=int(morphology_grid_size)
    )
    minus_morphology = _cap_vorticity_morphology(
        field, time=float(morphology_time), coefficient=-half_width, grid_size=int(morphology_grid_size)
    )
    plus_morphology = _cap_vorticity_morphology(
        field, time=float(morphology_time), coefficient=half_width, grid_size=int(morphology_grid_size)
    )
    trial_minus = _cap_vorticity_morphology(
        field,
        time=float(morphology_time),
        coefficient=-PREVIOUS_CAPACITY_TRIAL,
        grid_size=int(morphology_grid_size),
    )
    trial_plus = _cap_vorticity_morphology(
        field,
        time=float(morphology_time),
        coefficient=PREVIOUS_CAPACITY_TRIAL,
        grid_size=int(morphology_grid_size),
    )

    q99_base = float(baseline_morphology["axial_q99_over_Zp"])
    envelope_q99 = max(float(minus_morphology["axial_q99_over_Zp"]), float(plus_morphology["axial_q99_over_Zp"]))
    trial_q99 = max(float(trial_minus["axial_q99_over_Zp"]), float(trial_plus["axial_q99_over_Zp"]))
    q90_base = float(baseline_morphology["axial_q90_over_Zp"])
    envelope_q90 = max(float(minus_morphology["axial_q90_over_Zp"]), float(plus_morphology["axial_q90_over_Zp"]))
    trial_q90 = max(float(trial_minus["axial_q90_over_Zp"]), float(trial_plus["axial_q90_over_Zp"]))

    capacity_amplitude_ratio = PREVIOUS_CAPACITY_TRIAL / half_width
    result = {
        "task_id": TASK_ID,
        "source_candidate_sha256": str(field.candidate_sha256),
        "basis": "AXIAL_CAP_BANDED_C4_ODD_Z_POLOIDAL",
        "basis_parameters_added": {
            "fixed_spatial_shapes": 1,
            "scalar_coefficients_if_materialized": 1,
            "fitted_shape_parameters": 0,
            "temporal_degrees_added": 0,
        },
        "quadrature_orders": list(orders),
        "quadrature_max_relative_coefficient_change": _quadrature_refinement(rows_by_order),
        "reference_time": reference_time,
        "reference_energy_target": target,
        "reference_energy_abs_tolerance": tolerance,
        "reference_energy_quadratic": dict(reference),
        "normalization_roots_within_diagnostic_guard": _normalization_roots(
            reference, target=target, diagnostic_search_guard=search_guard
        ),
        "zero_centered_energy_envelope": dict(envelope),
        "previous_capacity_trial_abs_coefficient": PREVIOUS_CAPACITY_TRIAL,
        "previous_capacity_trial_to_energy_envelope_amplitude_ratio": float(capacity_amplitude_ratio),
        "previous_capacity_trial_inside_fixed_parent_energy_envelope": bool(PREVIOUS_CAPACITY_TRIAL <= half_width),
        "validation_time_energy_rows": validation_rows,
        "validation_time_energy_range_satisfied_over_envelope": bool(validation_range_satisfied),
        "direct_energy_replay_checks": direct_checks,
        "morphology": {
            "time": float(morphology_time),
            "grid_size": int(morphology_grid_size),
            "baseline": baseline_morphology,
            "negative_energy_envelope_endpoint": minus_morphology,
            "positive_energy_envelope_endpoint": plus_morphology,
            "negative_endpoint_delta_from_baseline": _morphology_delta(baseline_morphology, minus_morphology),
            "positive_endpoint_delta_from_baseline": _morphology_delta(baseline_morphology, plus_morphology),
            "negative_previous_capacity_trial": trial_minus,
            "positive_previous_capacity_trial": trial_plus,
            "baseline_axial_q90_over_Zp": q90_base,
            "maximum_energy_envelope_axial_q90_over_Zp": envelope_q90,
            "maximum_previous_trial_axial_q90_over_Zp": trial_q90,
            "baseline_axial_q99_over_Zp": q99_base,
            "maximum_energy_envelope_axial_q99_over_Zp": envelope_q99,
            "maximum_previous_trial_axial_q99_over_Zp": trial_q99,
            "energy_envelope_q90_moved_on_this_grid": bool(abs(envelope_q90 - q90_base) > 0.0),
            "energy_envelope_q99_moved_on_this_grid": bool(abs(envelope_q99 - q99_base) > 0.0),
        },
        "routing_contract": (
            "This audit decides only whether the previously demonstrated cap-localized tip morphology survives the unchanged "
            "fixed-parent reference-energy gate. If useful q90/q99 movement survives inside the zero-centered envelope, "
            "Agent 1 may materialize only this one bounded coefficient and Agent 2/3 must run fresh candidate-specific "
            "energy/sign/support/divergence and held-out full momentum checks. If gross tip movement disappears inside the "
            "fixed-parent envelope, do not add another cap or odd-q basis: next test an explicitly governed joint energy "
            "renormalization for this same mode, and if that fails route the remaining tip defect to base support/taper geometry."
        ),
        "truth_boundary": dict(TRUTH_BOUNDARY),
    }
    return result


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
