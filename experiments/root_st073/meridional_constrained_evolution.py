"""Bounded state-dependent constrained evolution of the meridional controls.

This is a one-step local tangent experiment.  The velocity state starts at
zero at the saved ``k0`` and is advanced by ``delta_k * slope``.  At each
state a fresh :class:`StateCache` and :class:`StateConstraints` instance is
assembled; no initial cone normal or constraint matrix is reused.  The four
moment rows used by the fit come from the order-96 fixed-radius integrated
operator, while order-48 :class:`StateConstraints` supplies the current cone
geometry.  The order-48 moment rows are retained as a discrepancy diagnostic.

The complete vector residual is fitted with physical-volume quadrature and a
scaled equality/inequality solve.  Independent order-96 moment and order-64
cone replays are run on the actual endpoint ``StatefulMean`` fields.  This is
an instantaneous local trajectory diagnostic, not a PDE, continuum, global
smoothness, finite-time, or scale-recursion result.
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.optimize import linprog


ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from broad_meridional_constrained import constrained_volume_fit  # noqa: E402
from broad_meridional_momentum import (  # noqa: E402
    RADIAL_QUADRATURE_BREAKS,
    evaluate_field,
    quadrature_nodes,
)
from broad_shear_dynamic_control import load_saved_field  # noqa: E402
from grouped_joined_field import install_in_field  # noqa: E402
from integrated_state_moments import assemble as assemble_integrated_moments  # noqa: E402
from meridional_state_cache import (  # noqa: E402
    StateCache,
    StatefulMean,
    value_and_pressure_modes,
)
from meridional_state_constraints import StateConstraints  # noqa: E402
from midplane_integrated_moment_balance import evaluate as replay_moments  # noqa: E402
from midplane_outer_residual_source import outer_cones  # noqa: E402


OUTPUT_PATH = ROOT / "meridional_constrained_evolution.json"
MOMENT_ORDER = 96
CONSTRAINT_ORDER = 48
REPLAY_CONE_ORDER = 64
FIT_ETA_ORDER = 3
FIT_RADIAL_ORDER = 5
HOLDOUT_ETA_ORDER = 4
HOLDOUT_RADIAL_ORDER = 4
ANGLE_SHIFT = np.pi / 7.0
RADIAL_BOUNDS = (0.01, 0.99)
DELTA_K = 1.0e-4
OPTIONAL_DELTA_K = 1.0e-3


def _save(report):
    OUTPUT_PATH.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")


def _array(value):
    return np.asarray(value, dtype=float)


def _list(value):
    return _array(value).tolist()


def _metric_summary(metrics):
    keys = (
        "point_count",
        "momentum_max",
        "momentum_volume_L2",
        "momentum_volume_RMS",
        "divergence_max",
        "divergence_volume_L2",
    )
    return {key: int(metrics[key]) if key == "point_count"
            else float(metrics[key]) for key in keys}


def _locations(saved):
    return [(float(row["label"]["eta"]), float(row["label"]["y"]))
            for row in saved["cone_problem"]["diagnostics"]]


def _volume_grid(inner, tau, eta_order, radial_order, angle_shift=0.0):
    return quadrature_nodes(
        inner,
        tau,
        eta_order,
        radial_order,
        1,
        radial_bounds=RADIAL_BOUNDS,
        radial_breaks=RADIAL_QUADRATURE_BREAKS,
        angle_shift=angle_shift,
    )


def _constraints_payload(engine_result):
    return {
        "E": _list(engine_result["E"]),
        "m": _list(engine_result["m"]),
        "A": _list(engine_result["A"]),
        "b": _list(engine_result["b"]),
        "geometry_admissible": bool(engine_result["geometry_admissible"]),
        "diagnostics": engine_result["diagnostics"],
    }


def _check_constraints(E, m, A, b, coefficients):
    x = _array(coefficients)
    equality = _array(E) @ x + _array(m)
    cone = _array(A) @ x + _array(b)
    return {
        "equality_residual_max": float(np.max(np.abs(equality))),
        "minimum_cone_margin": float(np.min(cone)),
        "equality_residual": _list(equality),
        "cone_margin_min": float(np.min(cone)),
        "valid": bool(np.isfinite(x).all()
                      and np.max(np.abs(equality)) <= 1.0e-7
                      and np.min(cone) >= -1.0e-7),
    }


def _scaled_constraint_system(E, m, A, b, column_scales=None):
    """Scale the current finite constraints for a robust linprog phase one."""

    E = _array(E)
    m = _array(m)
    A = _array(A)
    b = _array(b)
    if column_scales is None:
        raw = np.vstack((E, A))
        column_scales = np.maximum(np.linalg.norm(raw, axis=0), 1.0e-12)
    else:
        column_scales = np.maximum(_array(column_scales), 1.0e-12)
    eq = E / column_scales[None, :]
    ub = A / column_scales[None, :]
    eq_rows = np.maximum(np.maximum(np.max(np.abs(eq), axis=1), np.abs(m)),
                         1.0e-12)
    ub_rows = np.maximum(np.maximum(np.max(np.abs(ub), axis=1), np.abs(b)),
                         1.0e-12)
    return (eq / eq_rows[:, None], -m / eq_rows,
            ub / ub_rows[:, None], b / ub_rows,
            column_scales, eq_rows, ub_rows)


def _phase_one(E, m, A, b, column_scales=None):
    """Find a current feasible point, with explicit equality/inequality slacks."""

    eq, eq_rhs, ub, ub_rhs, scales, eq_rows, ub_rows = \
        _scaled_constraint_system(E, m, A, b, column_scales)
    n = len(scales)
    neq = len(eq_rhs)
    nub = len(ub_rhs)
    # E z - e+ + e- = rhs; -A z - s <= b.  Positive phase-one slack
    # measures an actual constraint violation in the normalized rows.
    c = np.r_[np.zeros(n), np.ones(neq * 2 + nub)]
    A_eq = np.zeros((neq, n + 2 * neq + nub))
    A_eq[:, :n] = eq
    A_eq[:, n:n + neq] = -np.eye(neq)
    A_eq[:, n + neq:n + 2 * neq] = np.eye(neq)
    A_ub = np.zeros((nub, n + 2 * neq + nub))
    A_ub[:, :n] = -ub
    A_ub[:, n + 2 * neq:] = -np.eye(nub)
    bounds = [(None, None)] * n + [(0.0, None)] * (2 * neq + nub)
    phase = linprog(
        c,
        A_ub=A_ub,
        b_ub=ub_rhs,
        A_eq=A_eq,
        b_eq=eq_rhs,
        bounds=bounds,
        method="highs",
        options={"primal_feasibility_tolerance": 1.0e-9},
    )
    if phase.x is None or not np.isfinite(phase.x).all():
        return None, {
            "success": bool(phase.success),
            "status": int(phase.status),
            "message": str(phase.message),
            "phase_one_objective": float("inf"),
            "column_scales": _list(scales),
        }
    z = phase.x[:n]
    x = z / scales
    phase_objective = float(np.sum(phase.x[n:]))
    checked = _check_constraints(E, m, A, b, x)
    info = {
        "success": bool(phase.success),
        "status": int(phase.status),
        "message": str(phase.message),
        "phase_one_objective": phase_objective,
        "column_scales": _list(scales),
        "equality_row_scales": _list(eq_rows),
        "cone_row_scales": _list(ub_rows),
        "constraint_check": checked,
    }
    if not checked["valid"]:
        return None, info
    return x, info


def _fit_at_state(base, values, pressures, inner, k, state, locations,
                  report, label):
    """Assemble fresh state constraints and solve the volume fit."""

    tau = 0.5 * 2.0 ** (-float(k))
    points, weights = _volume_grid(inner, tau, FIT_ETA_ORDER, FIT_RADIAL_ORDER)
    cache_started = time.perf_counter()
    cache = StateCache(base, values, pressures, points, float(k))
    offset, columns = cache.derivative_problem(state)
    cache_seconds = time.perf_counter() - cache_started
    report["stages"].append({
        "stage": f"{label}_volume_cache",
        "point_count": int(len(points)),
        "cache_seconds": float(cache_seconds),
    })
    _save(report)

    constraint_started = time.perf_counter()
    state_engine = StateConstraints(
        base,
        values,
        pressures,
        float(k),
        locations,
        RADIAL_QUADRATURE_BREAKS,
        order=CONSTRAINT_ORDER,
    )
    builder = state_engine.assemble(state)
    report["stages"].append({
        "stage": f"{label}_state_constraints",
        "order": CONSTRAINT_ORDER,
        "point_count": int(len(state_engine.cache.points)),
        "assembly_seconds": float(time.perf_counter() - constraint_started),
        "geometry_admissible": bool(builder["geometry_admissible"]),
    })
    _save(report)
    if not builder["geometry_admissible"]:
        return None, {
            "status": "geometry_inadmissible",
            "k": float(k),
            "state": _list(state),
            "constraints": _constraints_payload(builder),
            "cache_point_count": int(len(points)),
        }

    # Use the high-order fixed-radius moment operator for the fit.  Retain
    # StateConstraints' pointwise-panel rows to expose the order discrepancy.
    moment_started = time.perf_counter()
    high_E, high_m = assemble_integrated_moments(
        base, values, pressures, _array(state), float(k),
        RADIAL_QUADRATURE_BREAKS, order=MOMENT_ORDER,
    )
    moment_seconds = time.perf_counter() - moment_started
    high_E, high_m = _array(high_E), _array(high_m)
    builder_E, builder_m = _array(builder["E"]), _array(builder["m"])
    moment_discrepancy = {
        "constraint_order": CONSTRAINT_ORDER,
        "replay_order": MOMENT_ORDER,
        "E_max_abs_difference": float(np.max(np.abs(builder_E - high_E))),
        "m_max_abs_difference": float(np.max(np.abs(builder_m - high_m))),
        "builder_moment_for_zero_control_max": float(np.max(np.abs(builder_m))),
        "high_order_moment_for_zero_control_max": float(np.max(np.abs(high_m))),
        "assembly_seconds": float(moment_seconds),
    }
    A, b = _array(builder["A"]), _array(builder["b"])
    feasible, feasible_info = _phase_one(high_E, high_m, A, b)
    stage = {
        "stage": f"{label}_phase_one",
        "phase_one": feasible_info,
        "moment_discrepancy": moment_discrepancy,
    }
    report["stages"].append(stage)
    _save(report)
    if feasible is None:
        return None, {
            "status": "constraints_infeasible_or_phase_one_failed",
            "k": float(k),
            "state": _list(state),
            "constraints": _constraints_payload(builder),
            "high_order_moment": {"E": _list(high_E), "m": _list(high_m)},
            "moment_discrepancy": moment_discrepancy,
            "phase_one": feasible_info,
            "cache_point_count": int(len(points)),
        }

    equality = {
        "E": high_E,
        "m": high_m,
        "known_point": _array(feasible),
    }
    cones = {"A": A, "b": b}
    volume = {
        "columns": _array(columns),
        "baseline_residual": _array(offset),
        "weights": _array(weights),
    }
    fit_started = time.perf_counter()
    coefficients, solve = constrained_volume_fit(volume, equality, cones)
    fit_seconds = time.perf_counter() - fit_started
    selected_check = _check_constraints(high_E, high_m, A, b, coefficients)
    if not selected_check["valid"]:
        # Keep the valid phase-one point as the state tangent if SLSQP moves
        # outside the current geometry.  The report records this guard.
        coefficients = _array(feasible)
        selected_check = _check_constraints(high_E, high_m, A, b, coefficients)
        solve["state_constraint_fallback"] = True
    else:
        solve["state_constraint_fallback"] = False
    solve["phase_one"] = feasible_info
    solve["selected_constraint_check"] = selected_check
    solve["fit_seconds"] = float(fit_seconds)
    record = {
        "status": "fit_complete",
        "k": float(k),
        "tau": float(tau),
        "state": _list(state),
        "slope": _list(coefficients[:len(values)]),
        "pressure": _list(coefficients[len(values):]),
        "coefficients": _list(coefficients),
        "fit_grid": {
            "eta_order": FIT_ETA_ORDER,
            "radial_order_per_split_interval": FIT_RADIAL_ORDER,
            "radial_split_breaks": list(RADIAL_QUADRATURE_BREAKS),
            "angle_order": 1,
            "point_count": int(len(points)),
            "physical_volume": float(np.sum(weights)),
            "jacobian": "2*pi*r*(dr/dy)*(dz/deta); mapped Gauss weights include interval factors",
        },
        "fit_cache": {
            "point_count": int(len(points)),
            "hspace": float(cache.hspace),
            "htime": float(cache.htime),
            "cache_seconds": float(cache_seconds),
        },
        "constraints": _constraints_payload(builder),
        "high_order_moment": {"E": _list(high_E), "m": _list(high_m)},
        "moment_discrepancy": moment_discrepancy,
        "phase_one": feasible_info,
        "solve": solve,
        "selected_constraint_check": selected_check,
    }
    report["stages"].append({
        "stage": f"{label}_fit_complete",
        "selected_constraint_check": selected_check,
        "solver_success": bool(solve["success"]),
        "used_solver_point": bool(solve["used_solver_point"]),
        "fit_seconds": float(fit_seconds),
    })
    _save(report)
    return record, None


def _stateful(base, values, pressures, state, slope, pressure, kref):
    return StatefulMean(
        base,
        values,
        pressures,
        _array(state),
        _array(slope),
        _array(pressure),
        float(kref),
    )


def _endpoint_metrics(field, inner, k, angle_shift=ANGLE_SHIFT):
    tau = 0.5 * 2.0 ** (-float(k))
    points, weights = _volume_grid(
        inner, tau, HOLDOUT_ETA_ORDER, HOLDOUT_RADIAL_ORDER, angle_shift
    )
    return _metric_summary(evaluate_field(field, points, tau, weights)), {
        "eta_order": HOLDOUT_ETA_ORDER,
        "radial_order_per_split_interval": HOLDOUT_RADIAL_ORDER,
        "radial_split_breaks": list(RADIAL_QUADRATURE_BREAKS),
        "angle_order": 1,
        "angle_shift": float(angle_shift),
        "point_count": int(len(points)),
        "physical_volume": float(np.sum(weights)),
        "jacobian": "2*pi*r*(dr/dy)*(dz/deta); theta integration is 2*pi",
    }


def _replays(field, inner, k, locations):
    """Run the independent high-order moment and cone replays."""

    moments = replay_moments(
        field,
        inner,
        float(k),
        MOMENT_ORDER,
        0.002,
        radial_breaks=RADIAL_QUADRATURE_BREAKS,
    )
    cones = outer_cones(
        field,
        float(k),
        order=REPLAY_CONE_ORDER,
        radial_breaks=RADIAL_QUADRATURE_BREAKS,
        locations=locations,
    )
    return {
        "moment_order": MOMENT_ORDER,
        "moments": _list(moments),
        "moment_max_abs": float(np.max(np.abs(moments))),
        "cone_order": REPLAY_CONE_ORDER,
        "cone_count": len(cones),
        "cone_pass_count": int(sum(row["cone_pass"] for row in cones)),
        "cone_rows": cones,
    }


def run():
    started = time.perf_counter()
    dynamic, source = load_saved_field()
    replacements = install_in_field(dynamic)
    k0 = float(source["k"])
    tau0 = 0.5 * 2.0 ** (-k0)
    base, values, pressures, value_names, pressure_names = \
        value_and_pressure_modes(dynamic, source)
    saved = json.loads((ROOT / "broad_meridional_constrained.json").read_text(
        encoding="utf-8"
    ))
    locations = _locations(saved)
    state0 = np.zeros(len(values), dtype=float)
    report = {
        "status": "initializing",
        "accepted": False,
        "pde_validated": False,
        "scale_recursion_established": False,
        "constraints_maintained": False,
        "k0": k0,
        "tau0": tau0,
        "delta_k": DELTA_K,
        "optional_delta_k": OPTIONAL_DELTA_K,
        "optional_delta_status": "not attempted; bounded one-step run only",
        "velocity_state_count": len(values),
        "slope_count": len(values),
        "pressure_count": len(pressures),
        "value_names": value_names,
        "pressure_names": pressure_names,
        "grouped_backend_replacements": int(replacements),
        "initial_state": _list(state0),
        "locations": [list(location) for location in locations],
        "constraint_builder": {
            "class": "StateConstraints",
            "order": CONSTRAINT_ORDER,
            "moment_fit_operator": "integrated_state_moments.assemble",
            "moment_order": MOMENT_ORDER,
            "cone_count": 27,
            "cone_row_count": 81,
            "radial_breaks": list(RADIAL_QUADRATURE_BREAKS),
            "geometry_rebuilt_each_k": True,
        },
        "fit_quadrature": {
            "eta_order": FIT_ETA_ORDER,
            "radial_order_per_split_interval": FIT_RADIAL_ORDER,
            "radial_breaks": list(RADIAL_QUADRATURE_BREAKS),
            "point_count": int(FIT_ETA_ORDER * FIT_RADIAL_ORDER
                                * (len(RADIAL_QUADRATURE_BREAKS) - 1)),
            "physical_weight": "2*pi*r*(dr/dy)*(dz/deta) dy deta dtheta",
        },
        "holdout_quadrature": {
            "eta_order": HOLDOUT_ETA_ORDER,
            "radial_order_per_split_interval": HOLDOUT_RADIAL_ORDER,
            "angle_shift": float(ANGLE_SHIFT),
            "radial_breaks": list(RADIAL_QUADRATURE_BREAKS),
            "point_count": int(HOLDOUT_ETA_ORDER * HOLDOUT_RADIAL_ORDER
                                * (len(RADIAL_QUADRATURE_BREAKS) - 1)),
            "physical_weight": "same mapped volume Jacobian with independent Gauss orders",
        },
        "stages": [],
        "scope": (
            "Bounded instantaneous local tangent step from zero velocity state at k0. "
            "At both k values StateCache and StateConstraints are rebuilt, with current "
            "cone geometry and high-order integrated moment rows. Independent physical "
            "holdout, moment, and cone replays are finite diagnostics. No global smooth "
            "trajectory, continuum, finite-time PDE, or scale-recursion acceptance claim."
        ),
    }
    _save(report)
    print(json.dumps({
        "stage": "initial_fit_started",
        "k0": k0,
        "delta_k": DELTA_K,
        "state_count": len(values),
        "pressure_count": len(pressures),
        "fit_point_count": report["fit_quadrature"]["point_count"],
    }), flush=True)

    initial, failure = _fit_at_state(
        base, values, pressures, dynamic.inner, k0, state0, locations, report,
        "initial",
    )
    if failure is not None:
        report["status"] = failure["status"]
        report["initial_failure"] = failure
        report["elapsed_seconds"] = float(time.perf_counter() - started)
        _save(report)
        print(json.dumps({"stage": "stopped", "status": report["status"]}), flush=True)
        return report
    report["initial_fit"] = initial
    report["status"] = "initial_fit_complete"
    _save(report)
    print(json.dumps({
        "stage": "initial_fit_complete",
        "phase_one_objective": initial["phase_one"]["phase_one_objective"],
        "moment_discrepancy": initial["moment_discrepancy"]["m_max_abs_difference"],
        "constraint_check": initial["selected_constraint_check"],
    }), flush=True)

    initial_field = _stateful(
        base,
        values,
        pressures,
        state0,
        initial["slope"],
        initial["pressure"],
        k0,
    )
    report["initial_metrics"], report["initial_metric_domain"] = \
        _endpoint_metrics(initial_field, dynamic.inner, k0)
    report["status"] = "initial_metrics_complete"
    _save(report)

    state1 = state0 + DELTA_K * _array(initial["slope"])
    k1 = k0 + DELTA_K
    report["endpoint_state"] = _list(state1)
    report["endpoint_k"] = float(k1)
    report["status"] = "endpoint_fit_started"
    _save(report)
    print(json.dumps({"stage": "endpoint_fit_started", "k": k1,
                      "state_norm": float(np.linalg.norm(state1))}), flush=True)

    endpoint, failure = _fit_at_state(
        base, values, pressures, dynamic.inner, k1, state1, locations, report,
        "endpoint",
    )
    if failure is not None:
        report["status"] = failure["status"]
        report["endpoint_failure"] = failure
        report["elapsed_seconds"] = float(time.perf_counter() - started)
        _save(report)
        print(json.dumps({"stage": "stopped", "status": report["status"]}), flush=True)
        return report
    report["endpoint_fit"] = endpoint
    report["status"] = "endpoint_fit_complete"
    _save(report)
    print(json.dumps({
        "stage": "endpoint_fit_complete",
        "phase_one_objective": endpoint["phase_one"]["phase_one_objective"],
        "moment_discrepancy": endpoint["moment_discrepancy"]["m_max_abs_difference"],
        "constraint_check": endpoint["selected_constraint_check"],
    }), flush=True)

    endpoint_field = _stateful(
        base,
        values,
        pressures,
        state1,
        endpoint["slope"],
        endpoint["pressure"],
        k1,
    )
    report["endpoint_metrics"], report["endpoint_metric_domain"] = \
        _endpoint_metrics(endpoint_field, dynamic.inner, k1)
    report["status"] = "endpoint_metrics_complete"
    _save(report)
    print(json.dumps({"stage": "endpoint_metrics_complete",
                      "momentum_rms": report["endpoint_metrics"]["momentum_volume_RMS"]}),
          flush=True)

    # Replay only after all fit data and endpoint metrics are checkpointed.
    report["status"] = "initial_replay_started"
    _save(report)
    report["initial_replay"] = _replays(initial_field, dynamic.inner, k0, locations)
    _save(report)
    print(json.dumps({"stage": "initial_replay_complete",
                      "moment_max": report["initial_replay"]["moment_max_abs"],
                      "cone_pass": report["initial_replay"]["cone_pass_count"]}),
          flush=True)

    report["status"] = "endpoint_replay_started"
    _save(report)
    report["endpoint_replay"] = _replays(endpoint_field, dynamic.inner, k1, locations)
    _save(report)
    report["constraints_maintained"] = bool(
        initial["selected_constraint_check"]["valid"]
        and endpoint["selected_constraint_check"]["valid"]
        and report["initial_replay"]["cone_pass_count"] == 27
        and report["endpoint_replay"]["cone_pass_count"] == 27
        and report["initial_replay"]["moment_max_abs"] < 1.0e-3
        and report["endpoint_replay"]["moment_max_abs"] < 1.0e-3
    )
    report["constraints_maintained_scope"] = (
        "Initial k0 and endpoint k0+delta_k checks only; no interior pointwise "
        "constraint trajectory or continuum certificate was evaluated."
    )
    report["status"] = "completed"
    report["elapsed_seconds"] = float(time.perf_counter() - started)
    _save(report)
    print(json.dumps({
        "stage": "completed",
        "status": report["status"],
        "constraints_maintained": report["constraints_maintained"],
        "initial_moment_max": report["initial_replay"]["moment_max_abs"],
        "endpoint_moment_max": report["endpoint_replay"]["moment_max_abs"],
        "initial_cone_pass": report["initial_replay"]["cone_pass_count"],
        "endpoint_cone_pass": report["endpoint_replay"]["cone_pass_count"],
        "elapsed_seconds": report["elapsed_seconds"],
    }), flush=True)
    return report


if __name__ == "__main__":
    run()
