"""Bounded state-dependent refit diagnostic for the meridional basis.

This experiment starts from the zero 25-component correction state at the
saved scale ``k0``.  At each short k step it rebuilds :class:`StateCache` at
the new scale and fits all 25 velocity slopes and 19 pressure directions to
the complete finite-difference momentum residual.  The state is then advanced
explicitly by ``state += dk * slope``.  A fixed-slope path is evaluated beside
the refreshed path at the same endpoint.

The fit is deliberately unconstrained: no old moment or cone row is reused at
the changed state, and no non-axisymmetric wave is included.  The endpoint
metrics are therefore a local tangent replay diagnostic rather than a PDE,
moment/cone, scale-recursion, or physical-acceptance result.
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from broad_meridional_momentum import (  # noqa: E402
    RADIAL_QUADRATURE_BREAKS,
    evaluate_field,
    quadrature_nodes,
)
from broad_shear_dynamic_control import load_saved_field  # noqa: E402
from grouped_joined_field import install_in_field  # noqa: E402
from meridional_state_cache import (  # noqa: E402
    StateCache,
    StatefulMean,
    value_and_pressure_modes,
)


OUTPUT_PATH = ROOT / "meridional_state_evolution.json"
# The first run used (2, 2) and is preserved in
# ``meridional_state_evolution_coarse.json``.  This denser same-domain fit
# resolves the eta-polynomial directions more faithfully without changing the
# endpoint steps or metric grid.
FIT_ETA_ORDER = 3
FIT_RADIAL_ORDER = 5
METRIC_ETA_ORDER = 4
METRIC_RADIAL_ORDER = 4
METRIC_ANGLE_SHIFT = np.pi / 7.0
RADIAL_BOUNDS = (0.01, 0.99)


def _save(report, output_path=OUTPUT_PATH):
    """Write a JSON-safe incremental report."""

    Path(output_path).write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")


def _as_float_list(value):
    return np.asarray(value, dtype=float).tolist()


def _metric_summary(metrics):
    """Drop per-point residual vectors while retaining physical diagnostics."""

    keys = (
        "point_count",
        "momentum_max",
        "momentum_volume_L2",
        "momentum_volume_RMS",
        "divergence_max",
        "divergence_volume_L2",
    )
    return {key: float(metrics[key]) if key != "point_count" else int(metrics[key])
            for key in keys}


def _fit_cache(cache, state, weights):
    """Fit the complete residual's slope and pressure columns.

    ``StateCache.slopes`` are physical-time derivatives of a local
    ``(k-kref)`` direction.  Therefore the fitted coefficient is ``da/dk``:
    its physical-time contribution is exactly ``(dk/dt) * da/dk``.  Weighted
    column scaling keeps the pressure and velocity directions numerically
    comparable without changing the physical residual objective.
    """

    state = np.asarray(state, dtype=float)
    weights = np.asarray(weights, dtype=float)
    offset, columns = cache.derivative_problem(state)
    if offset.shape[0] != weights.shape[0]:
        raise ValueError("StateCache points and weights have different lengths")
    count = columns.shape[-1]
    sqrt_weights = np.sqrt(weights)
    weighted_columns = columns * sqrt_weights[:, None, None]
    weighted_target = -offset * sqrt_weights[:, None]
    matrix = weighted_columns.reshape(-1, count)
    target = weighted_target.reshape(-1)
    scales = np.linalg.norm(matrix, axis=0)
    safe_scales = np.maximum(scales, 1.0e-30)
    scaled = matrix / safe_scales[None, :]
    # Dimensionless ridge only regularizes nearly invisible directions.  The
    # momentum residual remains the dominant objective.
    ridge = 1.0e-8
    augmented = np.vstack((scaled, np.sqrt(ridge) * np.eye(count)))
    augmented_target = np.concatenate((target, np.zeros(count)))
    unregularized_singular = np.linalg.svd(scaled, compute_uv=False)
    unregularized_rank = int(np.linalg.matrix_rank(scaled, tol=1.0e-12))
    scaled_coefficients, _, augmented_rank, augmented_singular = np.linalg.lstsq(
        augmented, augmented_target, rcond=1.0e-12
    )
    coefficients = scaled_coefficients / safe_scales
    predicted = offset + np.einsum("ncp,p->nc", columns, coefficients)
    before = float(np.sqrt(np.sum(weights[:, None] * offset * offset)))
    after = float(np.sqrt(np.sum(weights[:, None] * predicted * predicted)))
    positive = unregularized_singular[np.abs(unregularized_singular) > 1.0e-14]
    condition = float(positive[0] / positive[-1]) if len(positive) else float("inf")
    n_velocity = len(cache.slopes)
    return {
        "slope": coefficients[:n_velocity],
        "pressure": coefficients[n_velocity:],
        "fit_info": {
            "point_count": int(len(cache.points)),
            "weighted_l2_before": before,
            "weighted_l2_after": after,
            "residual_max_before": float(np.max(np.linalg.norm(offset, axis=1))),
            "residual_max_after": float(np.max(np.linalg.norm(predicted, axis=1))),
            "rank": unregularized_rank,
            "unregularized_rank": unregularized_rank,
            "augmented_rank": int(augmented_rank),
            "unregularized_singular_values": _as_float_list(unregularized_singular),
            "augmented_singular_values": _as_float_list(augmented_singular),
            "condition_scaled": condition,
            "ridge": ridge,
            "column_scales": _as_float_list(scales),
        },
    }


def _fit_record(k, state, result):
    return {
        "k": float(k),
        "tau": float(0.5 * 2.0 ** (-float(k))),
        "state": _as_float_list(state),
        "slope": _as_float_list(result["slope"]),
        "pressure": _as_float_list(result["pressure"]),
        "fit_info": result["fit_info"],
    }


def _fit_at(base, values, pressures, inner, k, state, fit_grid):
    """Rebuild a full cache and fit at one scale/state."""

    tau = 0.5 * 2.0 ** (-float(k))
    points, weights = quadrature_nodes(
        inner,
        tau,
        FIT_ETA_ORDER,
        FIT_RADIAL_ORDER,
        1,
        radial_bounds=RADIAL_BOUNDS,
        radial_breaks=RADIAL_QUADRATURE_BREAKS,
    )
    started = time.perf_counter()
    cache = StateCache(base, values, pressures, points, float(k))
    result = _fit_cache(cache, state, weights)
    result["assembly_seconds"] = float(time.perf_counter() - started)
    fit_grid.append({
        "k": float(k),
        "tau": float(tau),
        "point_count": int(len(points)),
        "assembly_seconds": result["assembly_seconds"],
    })
    return result


def _stateful(base, values, pressures, state, slope, pressure, kref):
    return StatefulMean(
        base,
        values,
        pressures,
        np.asarray(state, dtype=float),
        np.asarray(slope, dtype=float),
        np.asarray(pressure, dtype=float),
        float(kref),
    )


def _endpoint_metrics(field, inner, k):
    """Use a separate higher-order physical quadrature for endpoint checks."""

    tau = 0.5 * 2.0 ** (-float(k))
    points, weights = quadrature_nodes(
        inner,
        tau,
        METRIC_ETA_ORDER,
        METRIC_RADIAL_ORDER,
        1,
        radial_bounds=RADIAL_BOUNDS,
        radial_breaks=RADIAL_QUADRATURE_BREAKS,
        angle_shift=METRIC_ANGLE_SHIFT,
    )
    return _metric_summary(evaluate_field(field, points, tau, weights))


def _run_delta(report, delta, base, values, pressures, inner, k0, state0, fit0,
               *, nsteps=2, output_path=OUTPUT_PATH):
    """Run a fixed/refreshed comparison for one positive delta k."""

    delta = float(delta)
    nsteps = int(nsteps)
    if nsteps < 1:
        raise ValueError("nsteps must be positive")
    step = delta / nsteps
    key = f"{delta:.12g}"
    run = {
        "delta_k": delta,
        "step_count": nsteps,
        "step_delta_k": step,
        "status": "running",
        "refit_steps": [],
        "fit_grid": [],
    }
    report["delta_runs"][key] = run
    _save(report, output_path)

    state = np.asarray(state0, dtype=float).copy()
    k = float(k0)
    current = fit0
    for index in range(nsteps):
        # Explicit state update in k.  The new cache is assembled at the
        # changed state before the next slope/pressure refit.
        state = state + step * current["slope"]
        k = k + step
        current = _fit_at(base, values, pressures, inner, k, state, run["fit_grid"])
        record = _fit_record(k, state, current)
        record["step_index"] = int(index + 1)
        run["refit_steps"].append(record)
        run["state_norm"] = float(np.linalg.norm(state))
        run["slope_norm"] = float(np.linalg.norm(current["slope"]))
        _save(report, output_path)
        print(json.dumps({
            "stage": "refit_step",
            "delta_k": delta,
            "step": index + 1,
            "k": k,
            "state_norm": run["state_norm"],
            "slope_norm": run["slope_norm"],
            "assembly_seconds": current["assembly_seconds"],
        }), flush=True)

    fixed_state = np.asarray(state0, dtype=float) + delta * fit0["slope"]
    fixed_field = _stateful(
        base, values, pressures, fixed_state, fit0["slope"], fit0["pressure"], k0
    )
    refreshed_field = _stateful(
        base, values, pressures, state, current["slope"], current["pressure"], k
    )
    run["endpoint_k"] = float(k)
    run["fixed_state"] = _as_float_list(fixed_state)
    run["refreshed_state"] = _as_float_list(state)
    run["fixed_slope"] = _as_float_list(fit0["slope"])
    run["refreshed_slope"] = _as_float_list(current["slope"])
    run["fixed_pressure"] = _as_float_list(fit0["pressure"])
    run["refreshed_pressure"] = _as_float_list(current["pressure"])
    run["status"] = "endpoint_metrics_running"
    _save(report, output_path)
    print(json.dumps({"stage": "endpoint_metrics_started", "delta_k": delta}), flush=True)
    run["fixed_metrics"] = _endpoint_metrics(fixed_field, inner, k)
    _save(report, output_path)
    print(json.dumps({"stage": "fixed_metrics_done", "delta_k": delta}), flush=True)
    run["refreshed_metrics"] = _endpoint_metrics(refreshed_field, inner, k)
    run["status"] = "completed"
    run["refreshed_minus_fixed_rms"] = float(
        run["refreshed_metrics"]["momentum_volume_RMS"]
        - run["fixed_metrics"]["momentum_volume_RMS"]
    )
    _save(report, output_path)
    print(json.dumps({
        "stage": "delta_done",
        "delta_k": delta,
        "fixed_rms": run["fixed_metrics"]["momentum_volume_RMS"],
        "refreshed_rms": run["refreshed_metrics"]["momentum_volume_RMS"],
    }), flush=True)


def run(*, steps=2, deltas=(1.0e-3, 1.0e-2), output_path=OUTPUT_PATH):
    """Run the requested local refit deltas.

    The defaults reproduce the dense report.  ``steps`` and ``deltas`` are
    exposed for bounded sensitivity checks; an alternate output path is
    required when preserving an existing report.
    """

    steps = int(steps)
    if steps < 1:
        raise ValueError("steps must be positive")
    deltas = tuple(float(delta) for delta in deltas)
    if not deltas or any(delta <= 0.0 for delta in deltas):
        raise ValueError("deltas must contain positive values")
    output_path = Path(output_path)
    started = time.perf_counter()
    dynamic, source = load_saved_field()
    replacements = install_in_field(dynamic)
    k0 = float(source["k"])
    tau0 = 0.5 * 2.0 ** (-k0)
    base, values, pressures, value_names, pressure_names = value_and_pressure_modes(
        dynamic, source
    )
    state0 = np.zeros(len(values), dtype=float)
    report = {
        "status": "initializing",
        "accepted": False,
        "pde_validated": False,
        "scale_recursion_established": False,
        "constraints_maintained": False,
        "moment_cone_rows_reused": False,
        "nonaxisymmetric_waves_included": False,
        "k0": k0,
        "tau0": tau0,
        "velocity_state_count": int(len(values)),
        "slope_count": int(len(values)),
        "pressure_count": int(len(pressures)),
        "value_names": value_names,
        "pressure_names": pressure_names,
        "grouped_backend_replacements": int(replacements),
        "output_report": str(output_path),
        "requested_deltas": list(deltas),
        "refit_step_count": steps,
        "coarse_reference_report": "meridional_state_evolution_coarse.json",
        "fit_grid_revision": "dense_same_domain_eta3_radial5_165_points",
        "fit_quadrature": {
            "eta_order": FIT_ETA_ORDER,
            "radial_order_per_break": FIT_RADIAL_ORDER,
            "angle_order": 1,
            "radial_bounds": list(RADIAL_BOUNDS),
            "radial_breaks": list(RADIAL_QUADRATURE_BREAKS),
        },
        "metric_quadrature": {
            "eta_order": METRIC_ETA_ORDER,
            "radial_order_per_break": METRIC_RADIAL_ORDER,
            "angle_order": 1,
            "angle_shift": float(METRIC_ANGLE_SHIFT),
            "radial_bounds": list(RADIAL_BOUNDS),
            "radial_breaks": list(RADIAL_QUADRATURE_BREAKS),
        },
        "trajectory_interpretation": {
            "state_update": "explicit Euler in k: state_next = state + delta_k * local_slope",
            "fixed_path": "one initial slope and pressure fit held through the endpoint",
            "refreshed_path": f"a new StateCache, slope, and pressure fit at each endpoint of {steps} k steps",
            "endpoint_metrics_only": True,
            "globally_time_differentiable": False,
            "note": "Piecewise local tangent fields have slope jumps between refit knots; intermediate-time replay is not claimed.",
        },
        "delta_runs": {},
        "scope": (
            "Unconstrained state-dependent local tangent diagnostic. The 25 velocity "
            "state directions and 19 pressure directions are refit from the complete "
            "quadratic momentum residual at each changed state. Endpoint fields use "
            "independent physical quadrature. No old moment/cone constraints, waves, "
            "global trajectory smoothness, PDE acceptance, or scale recursion are claimed."
        ),
    }
    _save(report, output_path)
    print(json.dumps({
        "stage": "fit_started",
        "k0": k0,
        "tau0": tau0,
        "velocity_state_count": len(values),
        "pressure_count": len(pressures),
        "grouped_backend_replacements": replacements,
    }), flush=True)

    fit_grid = []
    fit0 = _fit_at(base, values, pressures, dynamic.inner, k0, state0, fit_grid)
    report["initial_fit"] = _fit_record(k0, state0, fit0)
    report["initial_fit"]["assembly_seconds"] = fit0["assembly_seconds"]
    report["initial_fit_grid"] = fit_grid
    report["status"] = "initial_fit_complete"
    _save(report, output_path)
    print(json.dumps({
        "stage": "initial_fit_complete",
        "assembly_seconds": fit0["assembly_seconds"],
        "weighted_before": fit0["fit_info"]["weighted_l2_before"],
        "weighted_after": fit0["fit_info"]["weighted_l2_after"],
        "rank": fit0["fit_info"]["rank"],
    }), flush=True)

    # One local endpoint check anchors the physical metric before advancing.
    initial_field = _stateful(
        base, values, pressures, state0, fit0["slope"], fit0["pressure"], k0
    )
    report["initial_metrics"] = _endpoint_metrics(initial_field, dynamic.inner, k0)
    report["status"] = "initial_metrics_complete"
    _save(report, output_path)
    print(json.dumps({"stage": "initial_metrics_complete",
                      "momentum_rms": report["initial_metrics"]["momentum_volume_RMS"]}),
          flush=True)

    # Keep requested scales bounded and use the explicit step count supplied
    # by the caller.
    for delta in deltas:
        _run_delta(report, delta, base, values, pressures, dynamic.inner,
                   k0, state0, fit0, nsteps=steps, output_path=output_path)

    coarse_path = ROOT / "meridional_state_evolution_coarse.json"
    if coarse_path.exists():
        coarse = json.loads(coarse_path.read_text(encoding="utf-8"))
        comparison = {
            "coarse_path": coarse_path.name,
            "coarse_fit_quadrature": coarse.get("fit_quadrature"),
            "dense_fit_quadrature": report["fit_quadrature"],
            "coarse_initial_metric_rms": coarse.get("initial_metrics", {}).get("momentum_volume_RMS"),
            "dense_initial_metric_rms": report.get("initial_metrics", {}).get("momentum_volume_RMS"),
            "endpoint_rms": {},
            "note": "The coarse report predates separate unregularized spectrum reporting; its saved rank/condition used the augmented ridge system.",
        }
        for delta_key, dense_run in report["delta_runs"].items():
            coarse_run = coarse.get("delta_runs", {}).get(delta_key, {})
            comparison["endpoint_rms"][delta_key] = {
                "coarse_fixed": coarse_run.get("fixed_metrics", {}).get("momentum_volume_RMS"),
                "coarse_refreshed": coarse_run.get("refreshed_metrics", {}).get("momentum_volume_RMS"),
                "dense_fixed": dense_run.get("fixed_metrics", {}).get("momentum_volume_RMS"),
                "dense_refreshed": dense_run.get("refreshed_metrics", {}).get("momentum_volume_RMS"),
            }
        report["coarse_comparison"] = comparison

    # Alternate-output sensitivity runs can compare directly against the
    # preserved dense default report without overwriting it.
    if output_path.resolve() != OUTPUT_PATH.resolve() and OUTPUT_PATH.exists():
        reference = json.loads(OUTPUT_PATH.read_text(encoding="utf-8"))
        reference_comparison = {
            "reference_path": OUTPUT_PATH.name,
            "reference_step_count": reference.get("refit_step_count", 2),
            "endpoint_state_and_metrics": {},
            "note": "This is a two-step versus requested-step sensitivity check; it does not establish convergence.",
        }
        for delta_key, refined_run in report["delta_runs"].items():
            reference_run = reference.get("delta_runs", {}).get(delta_key, {})
            reference_comparison["endpoint_state_and_metrics"][delta_key] = {
                "reference_refreshed_state": reference_run.get("refreshed_state"),
                "refined_refreshed_state": refined_run.get("refreshed_state"),
                "reference_refreshed_state_norm": (
                    float(np.linalg.norm(reference_run["refreshed_state"]))
                    if reference_run.get("refreshed_state") is not None else None
                ),
                "refined_refreshed_state_norm": (
                    float(np.linalg.norm(refined_run["refreshed_state"]))
                    if refined_run.get("refreshed_state") is not None else None
                ),
                "reference_refreshed_rms": reference_run.get("refreshed_metrics", {}).get("momentum_volume_RMS"),
                "refined_refreshed_rms": refined_run.get("refreshed_metrics", {}).get("momentum_volume_RMS"),
                "reference_fixed_rms": reference_run.get("fixed_metrics", {}).get("momentum_volume_RMS"),
                "refined_fixed_rms": refined_run.get("fixed_metrics", {}).get("momentum_volume_RMS"),
            }
        report["dense_reference_comparison"] = reference_comparison

    report["status"] = "completed"
    report["elapsed_seconds"] = float(time.perf_counter() - started)
    _save(report, output_path)
    print(json.dumps({"stage": "completed",
                      "elapsed_seconds": report["elapsed_seconds"]}), flush=True)
    return report


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--steps",
        type=int,
        default=2,
        help="number of explicit k refit steps per requested delta (default: 2)",
    )
    parser.add_argument(
        "--delta",
        type=float,
        action="append",
        dest="deltas",
        help="positive delta k to run; repeat for multiple values (default: 0.001 and 0.01)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=OUTPUT_PATH,
        help="JSON output path; use an alternate path to preserve an existing report",
    )
    args = parser.parse_args()
    run(
        steps=args.steps,
        deltas=tuple(args.deltas) if args.deltas else (1.0e-3, 1.0e-2),
        output_path=args.output,
    )
