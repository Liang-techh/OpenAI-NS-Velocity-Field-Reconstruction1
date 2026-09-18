"""Target-free energy-envelope audit for the integrated compact poloidal mode.

This Agent-7 increment adds no basis and materializes no candidate.  It asks a
narrow pre-materialization question for the already-integrated
``COMPACT_C4_ODD_Z_POLOIDAL`` response: how much zero-centered coefficient
freedom survives the unchanged reference-energy gate, and how much target-free
vorticity morphology can move inside that envelope?

The diagnostic uses the exact quadratic kinetic-energy dependence of
``u + a*delta_u`` under axisymmetric Gauss-Legendre quadrature.  The inherited
profile coefficient limit is used only as a finite diagnostic search guard; it
is *not* promoted to a materialization bound.  No pressure/force is fitted, no
held-out momentum residual is evaluated, no public image is used, and no
coefficient or candidate SHA is selected here.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Iterable

import numpy as np
from numpy.polynomial.legendre import leggauss

from .constrained_eq45_bipolar_compact_poloidal_capacity import (
    COMPACT_POLOIDAL,
    _compact_poloidal_basis_velocity,
    _compact_velocity,
    _source_field,
    _vorticity_morphology,
)

TASK_ID = "CR003-BIPOLAR-COMPACT-POLOIDAL-ENERGY-ENVELOPE-044"
PREVIOUS_CAPACITY_TRIAL = 0.25
TRUTH_BOUNDARY = {
    "velocity_changed": False,
    "candidate_artifact_changed": False,
    "new_basis_added": False,
    "compact_poloidal_coefficient_selected": False,
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
    """Return E(a)=constant + linear*a + quadratic*a^2."""
    points, volume = _axisymmetric_grid(order)
    base = np.asarray(field.at_points(points, float(time)), dtype=float)
    basis = np.asarray(_compact_poloidal_basis_velocity(field, points), dtype=float)
    if base.shape != basis.shape or base.shape != (len(points), 3):
        raise RuntimeError("unexpected velocity shape in energy audit")
    constant = 0.5 * float(np.sum(volume * np.sum(base * base, axis=1)))
    linear = float(np.sum(volume * np.sum(base * basis, axis=1)))
    quadratic = 0.5 * float(np.sum(volume * np.sum(basis * basis, axis=1)))
    if not np.all(np.isfinite((constant, linear, quadratic))) or quadratic <= 0.0:
        raise RuntimeError("invalid kinetic-energy quadratic")
    return {
        "constant": constant,
        "linear": linear,
        "quadratic": quadratic,
    }


def _quadratic_value(coefficients: dict[str, float], coefficient: float) -> float:
    a = float(coefficient)
    return float(
        coefficients["constant"]
        + coefficients["linear"] * a
        + coefficients["quadratic"] * a * a
    )


def _direct_energy(field, time: float, coefficient: float, order: int) -> float:
    points, volume = _axisymmetric_grid(order)
    velocity = np.asarray(_compact_velocity(field, points, float(time), float(coefficient)), dtype=float)
    return 0.5 * float(np.sum(volume * np.sum(velocity * velocity, axis=1)))


def _interval_extrema(coefficients: dict[str, float], half_width: float) -> tuple[float, float]:
    half_width = float(half_width)
    if not np.isfinite(half_width) or half_width < 0.0:
        raise ValueError("half_width must be finite and nonnegative")
    q = coefficients["quadratic"]
    b = coefficients["linear"]
    left = _quadratic_value(coefficients, -half_width)
    right = _quadratic_value(coefficients, half_width)
    minimum = min(left, right)
    if q > 0.0:
        vertex = -b / (2.0 * q)
        if -half_width <= vertex <= half_width:
            minimum = min(minimum, _quadratic_value(coefficients, vertex))
    return float(minimum), float(max(left, right))


def _largest_symmetric_half_width(
    coefficients: dict[str, float],
    *,
    lower: float,
    upper: float,
    diagnostic_search_guard: float,
) -> dict[str, Any]:
    lower = float(lower)
    upper = float(upper)
    guard = float(diagnostic_search_guard)
    if not np.isfinite(guard) or guard <= 0.0 or not lower < upper:
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
    coefficients: dict[str, float], target: float, diagnostic_search_guard: float
) -> list[float]:
    q = float(coefficients["quadratic"])
    b = float(coefficients["linear"])
    c = float(coefficients["constant"] - float(target))
    discriminant = b * b - 4.0 * q * c
    if discriminant < 0.0:
        return []
    root_disc = float(np.sqrt(max(discriminant, 0.0)))
    roots = [(-b - root_disc) / (2.0 * q), (-b + root_disc) / (2.0 * q)]
    guard = float(diagnostic_search_guard)
    return sorted({float(root) for root in roots if np.isfinite(root) and abs(root) <= guard})


def _core_signs(field, coefficient: float, times: Iterable[float]) -> dict[str, Any]:
    times = np.asarray(tuple(float(t) for t in times), dtype=float)
    tau = 1.0 - times
    points = np.column_stack(
        (0.1 * np.sqrt(tau), np.zeros_like(times), 0.1 * tau**0.495)
    )
    velocity = _compact_velocity(field, points, times, float(coefficient))
    # x>0,y=0, so Cartesian (u_x,u_y,u_z)=(u_r,u_theta,u_z).
    radial_ok = bool(np.all(velocity[:, 0] < 0.0))
    swirl_ok = bool(np.all(velocity[:, 1] > 0.0))
    axial_ok = bool(np.all(velocity[:, 2] > 0.0))
    return {
        "coefficient": float(coefficient),
        "radial_inward": radial_ok,
        "swirl_positive": swirl_ok,
        "axial_positive": axial_ok,
        "all_signs_satisfied": bool(radial_ok and swirl_ok and axial_ok),
        "minimum_sign_margins": {
            "minus_u_r": float(np.min(-velocity[:, 0])),
            "u_theta": float(np.min(velocity[:, 1])),
            "u_z": float(np.min(velocity[:, 2])),
        },
    }


def _quadrature_refinement(rows_by_order: dict[int, dict[float, dict[str, float]]]) -> dict[str, float]:
    orders = sorted(rows_by_order)
    if len(orders) < 2:
        return {"maximum_relative_coefficient_change": 0.0}
    coarse, fine = rows_by_order[orders[-2]], rows_by_order[orders[-1]]
    maximum = 0.0
    for time in fine:
        for key in ("constant", "linear", "quadratic"):
            denominator = max(abs(fine[time][key]), 1e-12)
            maximum = max(maximum, abs(coarse[time][key] - fine[time][key]) / denominator)
    return {"maximum_relative_coefficient_change": float(maximum)}


def audit_compact_poloidal_energy_envelope(
    *,
    quadrature_orders: Iterable[int] = (48, 96),
    morphology_time: float = 0.5,
    morphology_grid_size: int = 25,
) -> dict[str, Any]:
    constraints = _constraints()
    nontriviality = constraints["nontriviality"]
    validation = constraints["validation"]
    reference_time = float(nontriviality["reference_time"])
    target = float(nontriviality["reference_energy"])
    tolerance = float(nontriviality["reference_energy_abs_tolerance"])
    energy_lower = float(nontriviality["minimum_energy_each_validation_time"])
    energy_upper = float(nontriviality["maximum_energy_each_validation_time"])
    times = tuple(float(t) for t in validation["times"])
    orders = tuple(int(v) for v in quadrature_orders)
    if len(orders) < 2 or any(order < 16 for order in orders) or any(
        orders[i] >= orders[i + 1] for i in range(len(orders) - 1)
    ):
        raise ValueError("quadrature_orders must contain at least two strictly increasing orders >=16")

    field = _source_field()
    diagnostic_search_guard = float(field.parent.profile_basis.coefficient_limit)
    rows_by_order: dict[int, dict[float, dict[str, float]]] = {}
    for order in orders:
        rows_by_order[order] = {
            time: _energy_quadratic(field, time, order) for time in times
        }
    finest = rows_by_order[orders[-1]]
    if reference_time not in finest:
        finest[reference_time] = _energy_quadratic(field, reference_time, orders[-1])
    reference = finest[reference_time]
    envelope = _largest_symmetric_half_width(
        reference,
        lower=target - tolerance,
        upper=target + tolerance,
        diagnostic_search_guard=diagnostic_search_guard,
    )
    half_width = float(envelope["half_width"])
    roots = _normalization_roots(reference, target, diagnostic_search_guard)

    validation_rows = []
    validation_range_satisfied = True
    for time in times:
        minimum, maximum = _interval_extrema(finest[time], half_width)
        satisfied = bool(minimum >= energy_lower and maximum <= energy_upper)
        validation_range_satisfied = validation_range_satisfied and satisfied
        validation_rows.append(
            {
                "time": time,
                "baseline_energy": _quadratic_value(finest[time], 0.0),
                "minimum_energy_over_zero_centered_envelope": minimum,
                "maximum_energy_over_zero_centered_envelope": maximum,
                "validation_energy_range_satisfied": satisfied,
            }
        )

    baseline_morphology = _vorticity_morphology(
        field, time=float(morphology_time), grid_size=int(morphology_grid_size)
    )
    plus_morphology = _vorticity_morphology(
        field,
        time=float(morphology_time),
        mode=COMPACT_POLOIDAL,
        coefficient=half_width,
        grid_size=int(morphology_grid_size),
    )
    minus_morphology = _vorticity_morphology(
        field,
        time=float(morphology_time),
        mode=COMPACT_POLOIDAL,
        coefficient=-half_width,
        grid_size=int(morphology_grid_size),
    )

    core_sign_checks = [
        _core_signs(field, coefficient, times)
        for coefficient in (-half_width, 0.0, half_width)
    ]
    reference_baseline = _quadratic_value(reference, 0.0)
    zero_gate_satisfied = bool(target - tolerance <= reference_baseline <= target + tolerance)
    previous_trial_inside = bool(PREVIOUS_CAPACITY_TRIAL <= half_width + 1e-15)

    return {
        "task_id": TASK_ID,
        "source_candidate_sha256": field.sha256,
        "mode": COMPACT_POLOIDAL,
        "basis_growth": {
            "new_basis_shapes_added": 0,
            "new_temporal_degrees_added": 0,
            "coefficient_selected": False,
        },
        "energy_model": {
            "definition": "E(a)=E0+B*a+A*a^2 for u+a*delta_u",
            "axisymmetric_volume_quadrature": "2*pi*r dr dz on r in [0,2], z in [-2,2]",
            "quadrature_orders": list(orders),
            "finest_order": orders[-1],
            "refinement": _quadrature_refinement(rows_by_order),
        },
        "reference_energy_gate": {
            "time": reference_time,
            "target": target,
            "absolute_tolerance": tolerance,
            "allowed_interval": [target - tolerance, target + tolerance],
            "baseline_energy": reference_baseline,
            "zero_coefficient_gate_satisfied": zero_gate_satisfied,
            "quadratic_coefficients": reference,
            "diagnostic_search_guard_abs_coefficient": diagnostic_search_guard,
            "diagnostic_search_guard_is_materialization_bound": False,
            "largest_symmetric_zero_centered_envelope": envelope,
            "diagnostic_exact_normalization_roots_within_search_guard": roots,
            "normalization_roots_select_coefficient": False,
            "previous_capacity_trial_abs_coefficient": PREVIOUS_CAPACITY_TRIAL,
            "previous_capacity_trial_inside_zero_centered_energy_envelope": previous_trial_inside,
        },
        "validation_time_energy_envelope": {
            "required_range": [energy_lower, energy_upper],
            "all_times_satisfied": bool(validation_range_satisfied),
            "rows": validation_rows,
        },
        "core_sign_screen": {
            "all_envelope_endpoint_and_baseline_checks_satisfied": bool(
                all(row["all_signs_satisfied"] for row in core_sign_checks)
            ),
            "rows": core_sign_checks,
            "scope": "target-free core sign screen only; not full candidate acceptance",
        },
        "vorticity_morphology_at_energy_envelope": {
            "time": float(morphology_time),
            "grid_size": int(morphology_grid_size),
            "baseline": baseline_morphology,
            "negative_endpoint": minus_morphology,
            "positive_endpoint": plus_morphology,
            "axial_rms_span_over_Zp": float(
                abs(plus_morphology["axial_rms_over_Zp"] - minus_morphology["axial_rms_over_Zp"])
            ),
            "radial_rms_span_over_Rp": float(
                abs(plus_morphology["radial_rms_over_Rp"] - minus_morphology["radial_rms_over_Rp"])
            ),
            "axial_q90_span_over_Zp": float(
                abs(plus_morphology["axial_q90_over_Zp"] - minus_morphology["axial_q90_over_Zp"])
            ),
            "axial_q99_span_over_Zp": float(
                abs(plus_morphology["axial_q99_over_Zp"] - minus_morphology["axial_q99_over_Zp"])
            ),
            "outer_axial_enstrophy_fraction_span": float(
                abs(
                    plus_morphology["outer_axial_enstrophy_fraction"]
                    - minus_morphology["outer_axial_enstrophy_fraction"]
                )
            ),
            "scope": "target-free finite-amplitude fingerprint; not visual correspondence",
        },
        "routing_contract": (
            "Do not grow another basis before this already-integrated compact-poloidal direction is materialized and screened. "
            "The zero-centered energy envelope is diagnostic capacity evidence, not an admissible production bound. Agent 1 must "
            "declare any nonzero coefficient/bound and new candidate identity; Agent 2/3 must then rerun unchanged energy, sign, "
            "support/divergence and fresh full momentum validation. If useful morphology requires amplitudes outside this zero-centered "
            "reference-energy envelope, use an explicit joint renormalization/materialization experiment rather than silently widening "
            "the coefficient or relaxing E(0.25)=1+/-0.001."
        ),
        "truth_boundary": dict(TRUTH_BOUNDARY),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        default="artifacts/bipolar_compact_poloidal_energy_envelope/report.json",
    )
    parser.add_argument("--morphology-grid-size", type=int, default=25)
    args = parser.parse_args()
    report = audit_compact_poloidal_energy_envelope(
        morphology_grid_size=args.morphology_grid_size
    )
    path = Path(args.output)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
