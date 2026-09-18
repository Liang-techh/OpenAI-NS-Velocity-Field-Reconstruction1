"""Target-free time-ramp capacity audit for the existing axial-cap poloidal mode.

No new spatial basis is introduced.  The existing support-adjacent axial-cap
mode is driven by one fixed autonomous diagnostic schedule

    a_sigma(t) = sigma*A*(t-t0)/(t1-t0), sigma in {-1,+1}.

The schedule is exactly zero at the registered reference time t0, so it replays
the parent field there.  A, its sign, and this linear schedule are capacity
probes only: no production coefficient/bound or OpenAI hidden time is inferred.
A materialized time-dependent child would change u_t and therefore needs a new
identity and fresh pressure/force/divergence/full-momentum validation.
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

TASK_ID = "CR003-BIPOLAR-AXIAL-CAP-TEMPORAL-RAMP-CAPACITY-048"
DEFAULT_TERMINAL_MAGNITUDE = 0.25
MORPH_KEYS = (
    "axial_rms_over_Zp",
    "radial_rms_over_Rp",
    "axial_q90_over_Zp",
    "axial_q99_over_Zp",
    "outer_065_enstrophy_fraction",
    "outer_075_enstrophy_fraction",
)
TRUTH_BOUNDARY = {
    "velocity_changed": False,
    "canonical_velocity_changed": False,
    "candidate_artifact_changed": False,
    "new_spatial_basis_added": False,
    "temporal_schedule_materialized": False,
    "production_coefficient_selected": False,
    "production_sign_selected": False,
    "production_bound_selected": False,
    "candidate_sha_created": False,
    "pressure_or_force_changed": False,
    "held_out_pde_residual_evaluated": False,
    "public_image_used": False,
    "hidden_frame_time_inferred": False,
    "visualization_ready": False,
    "visual_correspondence_verified": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
}


def _linear_ramp_coefficient(
    time: float, *, start_time: float, end_time: float, terminal_coefficient: float
) -> float:
    values = tuple(float(v) for v in (time, start_time, end_time, terminal_coefficient))
    time, start_time, end_time, terminal_coefficient = values
    if not all(np.isfinite(v) for v in values):
        raise ValueError("ramp inputs must be finite")
    if end_time <= start_time:
        raise ValueError("end_time must be greater than start_time")
    tol = 64.0 * np.finfo(float).eps * max(1.0, abs(start_time), abs(end_time))
    if time < start_time - tol or time > end_time + tol:
        raise ValueError("time lies outside the registered ramp window")
    phase = (time - start_time) / (end_time - start_time)
    if abs(phase) <= tol:
        phase = 0.0
    elif abs(phase - 1.0) <= tol:
        phase = 1.0
    return float(terminal_coefficient * phase)


def _delta(base: dict[str, float], trial: dict[str, float]) -> dict[str, float]:
    return {f"{k}_delta": float(trial[k] - base[k]) for k in MORPH_KEYS}


def _trend(rows: list[dict[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key in ("axial_rms_over_Zp", "radial_rms_over_Rp", "axial_q90_over_Zp", "axial_q99_over_Zp"):
        trial = np.asarray([r["trial_morphology"][key] for r in rows], dtype=float)
        base = np.asarray([r["baseline_morphology"][key] for r in rows], dtype=float)
        tol = 256.0 * np.finfo(float).eps * max(1.0, float(np.max(np.abs(trial))))
        result[key] = {
            "trial_start_to_end_change": float(trial[-1] - trial[0]),
            "baseline_start_to_end_change": float(base[-1] - base[0]),
            "incremental_change_vs_baseline_trend": float(
                (trial[-1] - trial[0]) - (base[-1] - base[0])
            ),
            "nondecreasing_step_count": int(np.sum(np.diff(trial) >= -tol)),
            "step_count": int(trial.size - 1),
        }
    return result


def _core_row(field, time: float, coefficient: float) -> dict[str, Any]:
    tau = 1.0 - float(time)
    point = np.array([[0.1 * np.sqrt(tau), 0.0, 0.1 * tau**0.495]])
    u = np.asarray(_cap_velocity(field, point, float(time), float(coefficient))[0], dtype=float)
    return {
        "time": float(time),
        "coefficient": float(coefficient),
        "u_r": float(u[0]),
        "u_theta": float(u[1]),
        "u_z": float(u[2]),
        "signs_satisfied": bool(u[0] < 0.0 and u[1] > 0.0 and u[2] > 0.0),
    }


def _schedule_report(
    field,
    *,
    sign: int,
    magnitude: float,
    start_time: float,
    end_time: float,
    registered_times: tuple[float, ...],
    energy_by_time: dict[float, dict[str, float]],
    energy_lower: float,
    energy_upper: float,
    morphology_times: tuple[float, ...],
    morphology_grid_size: int,
    baseline_morphology: dict[float, dict[str, float]],
) -> dict[str, Any]:
    terminal = float(sign) * magnitude
    energy_rows = []
    core_rows = []
    for time in registered_times:
        coefficient = _linear_ramp_coefficient(
            time,
            start_time=start_time,
            end_time=end_time,
            terminal_coefficient=terminal,
        )
        energy = float(_energy(energy_by_time[time], coefficient))
        energy_rows.append(
            {
                "time": time,
                "coefficient": coefficient,
                "energy": energy,
                "validation_energy_range_satisfied": bool(energy_lower <= energy <= energy_upper),
            }
        )
        core_rows.append(_core_row(field, time, coefficient))

    morph_rows = []
    for time in morphology_times:
        coefficient = _linear_ramp_coefficient(
            time,
            start_time=start_time,
            end_time=end_time,
            terminal_coefficient=terminal,
        )
        trial = _cap_vorticity_morphology(
            field, time=time, coefficient=coefficient, grid_size=morphology_grid_size
        )
        base = baseline_morphology[time]
        morph_rows.append(
            {
                "time": time,
                "coefficient": coefficient,
                "baseline_morphology": base,
                "trial_morphology": trial,
                "delta_from_coefficient_zero_baseline": _delta(base, trial),
            }
        )

    first = morph_rows[0]
    final = morph_rows[-1]
    final_delta = final["delta_from_coefficient_zero_baseline"]
    replay_error = max(
        abs(float(first["trial_morphology"][k]) - float(first["baseline_morphology"][k]))
        for k in MORPH_KEYS
    )
    return {
        "sign": sign,
        "terminal_coefficient": terminal,
        "coefficient_time_derivative": float(terminal / (end_time - start_time)),
        "reference_time_coefficient": float(energy_rows[0]["coefficient"]),
        "reference_time_parent_replayed_exactly_by_schedule": bool(energy_rows[0]["coefficient"] == 0.0),
        "registered_time_energy": {
            "required_range": [energy_lower, energy_upper],
            "all_times_within_broad_registered_range": bool(
                all(r["validation_energy_range_satisfied"] for r in energy_rows)
            ),
            "rows": energy_rows,
        },
        "core_signs": {
            "all_times_satisfied": bool(all(r["signs_satisfied"] for r in core_rows)),
            "rows": core_rows,
        },
        "morphology_rows": morph_rows,
        "trend_summary": _trend(morph_rows),
        "reference_time_morphology_replay_max_abs_error": float(replay_error),
        "final_time_incremental_tip_reach": {
            "axial_q90_over_Zp_delta": float(final_delta["axial_q90_over_Zp_delta"]),
            "axial_q99_over_Zp_delta": float(final_delta["axial_q99_over_Zp_delta"]),
            "axial_rms_over_Zp_delta": float(final_delta["axial_rms_over_Zp_delta"]),
            "radial_rms_over_Rp_delta": float(final_delta["radial_rms_over_Rp_delta"]),
        },
        "gross_final_tip_reach_extended_on_this_grid": bool(
            final_delta["axial_q90_over_Zp_delta"] > 0.0
            and final_delta["axial_q99_over_Zp_delta"] > 0.0
        ),
    }


def audit_axial_cap_temporal_ramp_capacity(
    *,
    terminal_magnitude: float = DEFAULT_TERMINAL_MAGNITUDE,
    quadrature_orders: Iterable[int] = (64, 96),
    morphology_times: Iterable[float] | None = None,
    morphology_grid_size: int = 33,
) -> dict[str, Any]:
    magnitude = float(terminal_magnitude)
    if not np.isfinite(magnitude) or magnitude <= 0.0:
        raise ValueError("terminal_magnitude must be positive and finite")
    if morphology_grid_size < 17 or morphology_grid_size % 2 == 0:
        raise ValueError("morphology_grid_size must be odd and >=17")
    orders = tuple(int(v) for v in quadrature_orders)
    if len(orders) < 2 or any(v < 16 for v in orders) or any(a >= b for a, b in zip(orders, orders[1:])):
        raise ValueError("quadrature_orders must be strictly increasing and >=16")

    constraints = _constraints()
    nontriviality = constraints["nontriviality"]
    registered_times = tuple(float(v) for v in constraints["validation"]["times"])
    start_time, end_time = (float(v) for v in constraints["domain"]["time_interval"])
    reference_time = float(nontriviality["reference_time"])
    if not np.isclose(reference_time, start_time, rtol=0.0, atol=1e-15):
        raise ValueError("temporal-ramp audit requires reference_time == interval start")
    if registered_times[0] != start_time or registered_times[-1] != end_time:
        raise ValueError("validation times must span the registered interval endpoints")

    if morphology_times is None:
        morphology_times_tuple = (start_time, 0.5 * (start_time + end_time), end_time)
    else:
        morphology_times_tuple = tuple(float(v) for v in morphology_times)
    if len(morphology_times_tuple) < 2 or any(
        a >= b for a, b in zip(morphology_times_tuple, morphology_times_tuple[1:])
    ):
        raise ValueError("morphology_times must be strictly increasing")
    tol = 64.0 * np.finfo(float).eps
    if morphology_times_tuple[0] < start_time - tol or morphology_times_tuple[-1] > end_time + tol:
        raise ValueError("morphology_times must lie inside the registered interval")

    field = _source_field()
    if magnitude > float(field.parent.profile_basis.coefficient_limit):
        raise ValueError("terminal_magnitude exceeds inherited implementation guard")

    rows_by_order = {
        order: {time: _energy_quadratic(field, time, order) for time in registered_times}
        for order in orders
    }
    finest = rows_by_order[orders[-1]]
    baseline = {
        time: _cap_vorticity_morphology(
            field, time=time, coefficient=0.0, grid_size=morphology_grid_size
        )
        for time in morphology_times_tuple
    }
    common = dict(
        field=field,
        magnitude=magnitude,
        start_time=start_time,
        end_time=end_time,
        registered_times=registered_times,
        energy_by_time=finest,
        energy_lower=float(nontriviality["minimum_energy_each_validation_time"]),
        energy_upper=float(nontriviality["maximum_energy_each_validation_time"]),
        morphology_times=morphology_times_tuple,
        morphology_grid_size=morphology_grid_size,
        baseline_morphology=baseline,
    )
    negative = _schedule_report(sign=-1, **common)
    positive = _schedule_report(sign=1, **common)

    return {
        "task_id": TASK_ID,
        "source_candidate_sha256": str(field.sha256),
        "mode": AXIAL_CAP_POLOIDAL,
        "representation_increment": {
            "new_spatial_basis_shapes_added": 0,
            "diagnostic_temporal_shapes_added": 1,
            "temporal_shape": "linear registered-window ramp with exact zero at reference time",
            "new_production_fit_parameters_added": 0,
            "production_coefficient_selected": False,
            "production_sign_selected": False,
            "production_bound_selected": False,
        },
        "schedule": {
            "definition": "a_sigma(t)=sigma*A*(t-t0)/(t1-t0)",
            "start_time": start_time,
            "end_time": end_time,
            "reference_time": reference_time,
            "diagnostic_terminal_magnitude": magnitude,
            "morphology_times": list(morphology_times_tuple),
            "registered_validation_times": list(registered_times),
            "not_a_production_bound": True,
        },
        "energy_quadrature": {
            "orders": list(orders),
            "finest_order": orders[-1],
            "maximum_relative_quadratic_coefficient_change": _refinement(rows_by_order),
        },
        "baseline_vorticity_morphology": {str(t): baseline[t] for t in morphology_times_tuple},
        "negative_schedule": negative,
        "positive_schedule": positive,
        "both_signs_reference_time_parent_replay": bool(
            negative["reference_time_parent_replayed_exactly_by_schedule"]
            and positive["reference_time_parent_replayed_exactly_by_schedule"]
        ),
        "both_signs_broad_energy_range_passed": bool(
            negative["registered_time_energy"]["all_times_within_broad_registered_range"]
            and positive["registered_time_energy"]["all_times_within_broad_registered_range"]
        ),
        "both_signs_core_signs_passed": bool(
            negative["core_signs"]["all_times_satisfied"]
            and positive["core_signs"]["all_times_satisfied"]
        ),
        "routing_contract": (
            "Do not add another spatial cap/odd-q/swirl basis. If this one temporal ramp gives material finite-window "
            "elongation/tip-reach leverage while preserving registered energy/core checks, materialize one explicitly bounded "
            "time-dependent child and rerun fresh u_t-aware pressure/force/divergence/full-momentum validation plus actual 3D "
            "streamline/vorticity comparison. Static-cap PDE evidence cannot transfer because a'(t) changes u_t. If temporal "
            "leverage is weak, keep the cap coefficient time independent and route remaining mismatch to support/taper geometry."
        ),
        "truth_boundary": dict(TRUTH_BOUNDARY),
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--terminal-magnitude", type=float, default=DEFAULT_TERMINAL_MAGNITUDE)
    parser.add_argument("--quadrature-order", type=int, action="append", dest="orders")
    parser.add_argument("--morphology-grid-size", type=int, default=33)
    args = parser.parse_args(argv)
    report = audit_axial_cap_temporal_ramp_capacity(
        terminal_magnitude=args.terminal_magnitude,
        quadrature_orders=tuple(args.orders) if args.orders else (64, 96),
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
