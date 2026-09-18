"""Target-free temporal-capacity audit for the existing axial-cap poloidal mode.

The support-adjacent axial-cap direction has already demonstrated finite-amplitude
control of gross q90/q99 vorticity reach.  This module does not add another
spatial basis.  It asks one narrower representation question: can one fixed
linear-in-window coefficient schedule turn that already identified spatial
channel on over the registered time interval, while leaving the reference-time
field exactly unchanged?

For diagnostic terminal magnitude A and sign sigma in {-1,+1},

    a_sigma(t) = sigma * A * (t - t0) / (t1 - t0),

with t0=.25 and t1=.75 taken from the registered constraints.  Thus a(t0)=0
exactly.  The schedule, terminal magnitude, and signs are autonomous capacity
probes; no production coefficient, sign, bound, hidden frame time, or OpenAI
numerical target is selected here.

The audit measures registered-time energy/core-sign behavior and target-free 3D
vorticity morphology at t0, midpoint, and t1.  It is not a Navier--Stokes
validation.  A materialized time-dependent child would have a new u_t term and
must receive a new identity plus fresh pressure/force/divergence/full-momentum
validation.
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

_MORPH_KEYS = (
    "axial_rms_over_Zp",
    "radial_rms_over_Rp",
    "axial_q90_over_Zp",
    "axial_q99_over_Zp",
    "outer_065_enstrophy_fraction",
    "outer_075_enstrophy_fraction",
)


def _linear_ramp_coefficient(
    time: float,
    *,
    start_time: float,
    end_time: float,
    terminal_coefficient: float,
) -> float:
    time = float(time)
    start_time = float(start_time)
    end_time = float(end_time)
    terminal_coefficient = float(terminal_coefficient)
    values = (time, start_time, end_time, terminal_coefficient)
    if not all(np.isfinite(value) for value in values):
        raise ValueError("ramp inputs must be finite")
    if not end_time > start_time:
        raise ValueError("end_time must be greater than start_time")
    tolerance = 64.0 * np.finfo(float).eps * max(1.0, abs(start_time), abs(end_time))
    if time < start_time - tolerance or time > end_time + tolerance:
        raise ValueError("time lies outside the registered ramp window")
    phase = (time - start_time) / (end_time - start_time)
    if abs(phase) <= tolerance:
        phase = 0.0
    elif abs(phase - 1.0) <= tolerance:
        phase = 1.0
    return float(terminal_coefficient * phase)


def _morphology_delta(base: dict[str, float], trial: dict[str, float]) -> dict[str, float]:
    return {f"{key}_delta": float(trial[key] - base[key]) for key in _MORPH_KEYS}


def _trend_summary(rows: list[dict[str, Any]]) -> dict[str, Any]:
    if len(rows) < 2:
        raise ValueError("at least two morphology rows are required")
    metrics = {}
    for key in ("axial_rms_over_Zp", "radial_rms_over_Rp", "axial_q90_over_Zp", "axial_q99_over_Zp"):
        values = np.asarray([row["trial_morphology"][key] for row in rows], dtype=float)
        baseline = np.asarray([row["baseline_morphology"][key] for row in rows], dtype=float)
        tolerance = 256.0 * np.finfo(float).eps * max(1.0, float(np.max(np.abs(values))))
        metrics[key] = {
            "trial_start_to_end_change": float(values[-1] - values[0]),
            "baseline_start_to_end_change": float(baseline[-1] - baseline[0]),
            "incremental_change_vs_baseline_trend": float(
                (values[-1] - values[0]) - (baseline[-1] - baseline[0])
            ),
            "nondecreasing_step_count": int(np.sum(np.diff(values) >= -tolerance)),
            "step_count": int(values.size - 1),
        }
    return metrics


def _core_row(field, time: float, coefficient: float) -> dict[str, Any]:
    tau = 1.0 - float(time)
    point = np.array([[0.1 * np.sqrt(tau), 0.0, 0.1 * tau**0.495]], dtype=float)
    velocity = np.asarray(_cap_velocity(field, point, float(time), float(coefficient))[0], dtype=float)
    return {
        "time": float(time),
        "coefficient": float(coefficient),
        "u_r": float(velocity[0]),
        "u_theta": float(velocity[1]),
        "u_z": float(velocity[2]),
        "signs_satisfied": bool(velocity[0] < 0.0 and velocity[1] > 0.0 and velocity[2] > 0.0),
    }


def _signed_schedule_report(
    field,
    *,
    sign: float,
    magnitude: float,
    start_time: float,
    end_time: float,
    registered_times: tuple[float, ...],
    energy_by_time: dict[float, dict[str, float]],
    validation_energy_lower: float,
    validation_energy_upper: float,
    morphology_times: tuple[float, ...],
    morphology_grid_size: int,
    baseline_morphology: dict[float, dict[str, float]],
) -> dict[str, Any]:
    terminal = float(sign) * float(magnitude)
    slope = terminal / (end_time - start_time)

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
        energy_ok = bool(validation_energy_lower <= energy <= validation_energy_upper)
        energy_rows.append({
            "time": float(time),
            "coefficient": coefficient,
            "energy": energy,
            "validation_energy_range_satisfied": energy_ok,
        })
        core_rows.append(_core_row(field, time, coefficient))

    morphology_rows = []
    for time in morphology_times:
        coefficient = _linear_ramp_coefficient(
            time,
            start_time=start_time,
            end_time=end_time,
            terminal_coefficient=terminal,
        )
        trial = _cap_vorticity_morphology(
            field,
            time=float(time),
            coefficient=coefficient,
            grid_size=int(morphology_grid_size),
        )
        base = baseline_morphology[float(time)]
        morphology_rows.append({
            "time": float(time),
            "coefficient": coefficient,
            "baseline_morphology": base,
            "trial_morphology": trial,
            "delta_from_coefficient_zero_baseline": _morphology_delta(base, trial),
        })

    first = morphology_rows[0]
    final = morphology_rows[-1]
    initial_replay_error = max(
        abs(float(first["trial_morphology"][key]) - float(first["baseline_morphology"][key]))
        for key in _MORPH_KEYS
    )
    final_delta = final["delta_from_coefficient_zero_baseline"]
    return {
        "sign": int(np.sign(sign)),
        "terminal_coefficient": terminal,
        "coefficient_time_derivative": float(slope),
        "reference_time_coefficient": float(energy_rows[0]["coefficient"]),
        "reference_time_parent_replayed_exactly_by_schedule": bool(energy_rows[0]["coefficient"] == 0.0),
        "registered_time_energy": {
            "all_times_within_broad_registered_range": bool(
                all(row["validation_energy_range_satisfied"] for row in energy_rows)
            ),
            "required_range": [validation_energy_lower, validation_energy_upper],
            "rows": energy_rows,
        },
        "core_signs": {
            "all_times_satisfied": bool(all(row["signs_satisfied"] for row in core_rows)),
            "rows": core_rows,
        },
        "morphology_rows": morphology_rows,
        "trend_summary": _trend_summary(morphology_rows),
        "reference_time_morphology_replay_max_abs_error": float(initial_replay_error),
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

    orders = tuple(int(value) for value in quadrature_orders)
    if len(orders) < 2 or any(value < 16 for value in orders) or any(
        orders[index] >= orders[index + 1] for index in range(len(orders) - 1)
    ):
        raise ValueError("quadrature_orders must be strictly increasing and >=16")

    constraints = _constraints()
    nontriviality = constraints["nontriviality"]
    validation = constraints["validation"]
    start_time, end_time = (float(value) for value in constraints["time"]["interval"])
    registered_times = tuple(float(value) for value in validation["times"])
    reference_time = float(nontriviality["reference_time"])
    if not np.isclose(reference_time, start_time, rtol=0.0, atol=1e-15):
        raise ValueError("temporal-ramp audit requires reference_time == registered interval start")
    if registered_times[0] != start_time or registered_times[-1] != end_time:
        raise ValueError("validation times must span the registered interval endpoints")

    if morphology_times is None:
        morphology_times_tuple = (start_time, 0.5 * (start_time + end_time), end_time)
    else:
        morphology_times_tuple = tuple(float(value) for value in morphology_times)
    if len(morphology_times_tuple) < 2 or any(
        morphology_times_tuple[index] >= morphology_times_tuple[index + 1]
        for index in range(len(morphology_times_tuple) - 1)
    ):
        raise ValueError("morphology_times must be strictly increasing")
    tolerance = 64.0 * np.finfo(float).eps
    if morphology_times_tuple[0] < start_time - tolerance or morphology_times_tuple[-1] > end_time + tolerance:
        raise ValueError("morphology_times must lie inside the registered interval")

    field = _source_field()
    implementation_guard = float(field.parent.profile_basis.coefficient_limit)
    if magnitude > implementation_guard:
        raise ValueError("terminal_magnitude exceeds inherited implementation guard")

    rows_by_order = {
        order: {time: _energy_quadratic(field, time, order) for time in registered_times}
        for order in orders
    }
    finest = rows_by_order[orders[-1]]
    baseline_morphology = {
        float(time): _cap_vorticity_morphology(
            field,
            time=float(time),
            coefficient=0.0,
            grid_size=int(morphology_grid_size),
        )
        for time in morphology_times_tuple
    }

    energy_lower = float(nontriviality["minimum_energy_each_validation_time"])
    energy_upper = float(nontriviality["maximum_energy_each_validation_time"])
    negative = _signed_schedule_report(
        field,
        sign=-1.0,
        magnitude=magnitude,
        start_time=start_time,
        end_time=end_time,
        registered_times=registered_times,
        energy_by_time=finest,
        validation_energy_lower=energy_lower,
        validation_energy_upper=energy_upper,
        morphology_times=morphology_times_tuple,
        morphology_grid_size=morphology_grid_size,
        baseline_morphology=baseline_morphology,
    )
    positive = _signed_schedule_report(
        field,
        sign=1.0,
        magnitude=magnitude,
        start_time=start_time,
        end_time=end_time,
        registered_times=registered_times,
        energy_by_time=finest,
        validation_energy_lower=energy_lower,
        validation_energy_upper=energy_upper,
        morphology_times=morphology_times_tuple,
        morphology_grid_size=morphology_grid_size,
        baseline_morphology=baseline_morphology,
    )

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
            "finest_order": int(orders[-1]),
            "maximum_relative_quadratic_coefficient_change": _refinement(rows_by_order),
        },
        "baseline_vorticity_morphology": {
            str(time): baseline_morphology[float(time)] for time in morphology_times_tuple
        },
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
            "This audit adds no spatial basis. If one fixed temporal ramp makes the already identified axial-cap channel produce "
            "a materially stronger finite-window elongation/tip-reach trend without breaking registered energy/core checks, the "
            "next candidate-level step is to materialize one explicitly bounded time-dependent child and rerun fresh u_t-aware "
            "pressure/force/divergence/full-momentum validation plus actual 3D streamline/vorticity comparison. Do not transfer "
            "static-cap PDE evidence because a'(t) changes u_t. If the ramp adds little temporal morphology leverage, keep the cap "
            "coefficient time independent and route remaining mismatch to support/taper geometry rather than adding more basis."
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
    orders = tuple(args.orders) if args.orders else (64, 96)
    report = audit_axial_cap_temporal_ramp_capacity(
        terminal_magnitude=args.terminal_magnitude,
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
