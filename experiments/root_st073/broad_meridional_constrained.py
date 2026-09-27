"""Constrained full-momentum fit for the broad meridional experiment.

This follows ``broad_meridional_momentum.py`` but retains the four saved
integrated moment equations and the 81 sampled cone inequalities from the
nineteen-control broad-shear solve.  The saved raw dynamic field is a known
feasible point: the new compact-pressure baseline is returned to that field by
setting the compact primitive coefficient to ``-1``.

All directions have zero instantaneous velocity at ``k0``.  The complete
Cartesian momentum residual is therefore affine in the new coefficients at
the reference time.  This remains a finite-grid instantaneous diagnostic;
it does not establish a continuum cone, finite-time evolution, forcing, PDE
validity, or scale recursion.
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.optimize import LinearConstraint, minimize


ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from affine_momentum import jets, momentum  # noqa: E402
from adaptive_bridge_moment_fit import moment_slices, outer_moments  # noqa: E402
from broad_meridional_momentum import (  # noqa: E402
    CompactPressureDirection,
    DYNAMIC_PATH,
    OUTPUT_PATH as UNCONSTRAINED_PATH,
    build_compact_baseline,
    build_modes_with_parameters,
    evaluate_field,
    quadrature_nodes,
    RADIAL_QUADRATURE_BREAKS,
)
from broad_shear_dynamic_control import (  # noqa: E402
    load_saved_field as load_dynamic_field,
    old_outer_locations,
    wave_locations,
)
from grouped_outer_cache import outer_cache  # noqa: E402
from midplane_integrated_moment_balance import evaluate as replay_moments  # noqa: E402
from midplane_outer_residual_source import outer_cones  # noqa: E402
from midplane_resolved_feasibility import ZeroBackground  # noqa: E402


OUTPUT_PATH = ROOT / "broad_meridional_constrained.json"
MARGIN = np.array([0.02, 0.0, 0.0], dtype=float)
COMPACT_ORDER = 8
ASSEMBLY_ORDER = 48
REPLAY_ORDER = 96
REPLAY_CONE_ORDER = 64


def _save(report):
    OUTPUT_PATH.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")


def _locations(dynamic_report):
    wave_report = json.loads((ROOT / "broad_shear_wave_cone.json").read_text(encoding="utf-8"))
    old = old_outer_locations()
    wave = wave_locations(wave_report)
    locations = old + wave
    labels = ([{"eta": float(eta), "y": float(y), "set": "old_outer"}
               for eta, y in old]
              + [{"eta": float(eta), "y": float(y), "set": "wave_center"}
                 for eta, y in wave])
    diagnostics = dynamic_report["cone_problem"]["diagnostics"]
    if len(diagnostics) != len(locations):
        raise ValueError("Saved cone diagnostics do not match the 27-node layout")
    return locations, labels, diagnostics


def _raw_moment_from_saved(dynamic_report):
    problem = dynamic_report["moment_problem"]
    old_offset = np.asarray(problem["offset"], dtype=float)
    old_linear = np.asarray(problem["linear"], dtype=float)
    control = np.asarray(dynamic_report["control"], dtype=float)
    if old_linear.shape != (4, 19) or control.shape != (19,):
        raise ValueError("Saved dynamic moment problem has unexpected shape")
    return old_offset + old_linear @ control


def assemble_moment_rows(dynamic, dynamic_report, baseline, modes, k0,
                         radial_breaks):
    """Assemble new mode moment columns without re-solving old controls."""

    zero = ZeroBackground(dynamic.base)
    regular_modes = list(modes[:-1])
    data = moment_slices(
        dynamic.inner,
        dynamic.base,
        zero,
        orders=(k0,),
        n=ASSEMBLY_ORDER,
        unit_fields=regular_modes,
        unit_fields_are_deltas=True,
        radial_breaks=radial_breaks,
    )[0]
    mode_count = len(modes)
    regular_columns = np.column_stack([
        outer_moments(data, np.eye(mode_count - 1)[index])
        for index in range(mode_count - 1)
    ])
    # The compact primitive is pressure-only.  At k0, its complete residual
    # column equals compact-baseline residual minus raw dynamic residual.
    # Assemble that column with two direct jets rather than invoking the
    # expensive radial primitive for each finite-difference stencil.
    panel_points = []
    for panel in data["panels"]:
        eta = float(panel["eta"])
        z = float(dynamic.inner.from_similarity(
            [dynamic.inner.p.X_max], [eta], data["tau"]
        )[0, 2])
        panel_points.extend((float(radius), 0.0, z)
                            for radius in panel["radii"])
    panel_points = np.asarray(panel_points, dtype=float)
    hs = 5.0e-4 * np.sqrt(dynamic.nu * data["tau"])
    # The primitive contributes pressure only.  Its tangential moment column
    # needs only p_z, so evaluate the same fourth-order z stencil directly;
    # radial velocity/time jets would be redundant and very expensive.
    ez = np.array([0.0, 0.0, 1.0])
    pm2 = baseline.increment(panel_points - 2.0 * hs * ez, data["tau"])
    pm = baseline.increment(panel_points - hs * ez, data["tau"])
    pp = baseline.increment(panel_points + hs * ez, data["tau"])
    pp2 = baseline.increment(panel_points + 2.0 * hs * ez, data["tau"])
    compact_delta = np.zeros((len(panel_points), 3), dtype=float)
    compact_delta[:, 2] = (pm2 - 8.0 * pm + 8.0 * pp - pp2) / (12.0 * hs)
    compact_column = []
    for panel in data["panels"]:
        local = compact_delta[panel["slice"]]
        rr = panel["radii"]
        ww = panel["weights"]
        ro = panel["outer_radius"]
        compact_column.extend((
            -np.dot(ww * rr**2, local[:, 1]) / ro**2,
            -np.dot(ww * rr, local[:, 2]) / ro,
        ))
    columns = np.column_stack((regular_columns, np.asarray(compact_column)))
    raw_moment = _raw_moment_from_saved(dynamic_report)
    compact_index = mode_count - 1
    baseline_moment = raw_moment + columns[:, compact_index]
    known = np.zeros(mode_count, dtype=float)
    known[compact_index] = -1.0
    return {
        "E": columns,
        "m": baseline_moment,
        "raw_moment": raw_moment,
        "known_point": known,
        "known_residual": columns @ known + baseline_moment,
        "panel_count": len(data["panels"]),
        "quadrature_order": ASSEMBLY_ORDER,
        "compact_column_direct_fd": True,
    }


def assemble_cone_rows(dynamic, dynamic_report, modes, k0, radial_breaks):
    """Assemble H-normalized cone rows from actual new-direction jets."""

    locations, labels, diagnostics = _locations(dynamic_report)
    zero = ZeroBackground(dynamic.base)
    regular_modes = list(modes[:-1])
    cache = outer_cache(
        zero,
        regular_modes,
        k=k0,
        order=ASSEMBLY_ORDER,
        locations=locations,
        radial_breaks=radial_breaks,
    )
    mode_count = len(modes)
    columns = np.zeros((len(diagnostics) * 3, mode_count), dtype=float)
    offsets = np.zeros(len(diagnostics) * 3, dtype=float)
    # The final mode is the compact pressure primitive.  Only p_z enters the
    # cone's axial target; evaluate its fourth-order z derivative once over
    # the grouped cache and append it to the 43 regular mode columns.
    primitive = modes[-1].primitive
    hs = 5.0e-4 * np.sqrt(dynamic.nu * cache["tau"])
    ez = np.array([0.0, 0.0, 1.0])
    cache_points = cache["points"]
    primitive_pm2 = primitive.increment(cache_points - 2.0 * hs * ez, cache["tau"])
    primitive_pm = primitive.increment(cache_points - hs * ez, cache["tau"])
    primitive_pp = primitive.increment(cache_points + hs * ez, cache["tau"])
    primitive_pp2 = primitive.increment(cache_points + 2.0 * hs * ez, cache["tau"])
    primitive_linear = np.zeros((1, len(cache_points), 3), dtype=float)
    primitive_linear[0, :, 2] = (
        primitive_pm2 - 8.0 * primitive_pm + 8.0 * primitive_pp - primitive_pp2
    ) / (12.0 * hs)
    all_linear = np.concatenate((np.asarray(cache["modes"][2]), primitive_linear), axis=0)
    metadata = []
    control = np.asarray(dynamic_report["control"], dtype=float)
    for index, (panel, diagnostic) in enumerate(zip(
            cache["panels"], diagnostics)):
        sl, radii, weights, radius = panel
        # The saved source_T is for the raw dynamic field.  Its H matrix and
        # geometry are unchanged because every new velocity direction vanishes
        # at k0.  Recover the raw dynamic target from saved old-direction rows.
        source_T = np.asarray(diagnostic["source_T"], dtype=float)
        old_unit_target = np.asarray(diagnostic["unit_target"], dtype=float)
        raw_target = source_T + old_unit_target @ control
        mode_linear = np.asarray(all_linear[:, sl], dtype=float)
        mode_target = np.stack([
            -np.einsum("n,pn->p", weights * radii**2 / radius**2,
                       mode_linear[:, :, 1]),
            -np.einsum("n,pn->p", weights * radii / radius,
                       mode_linear[:, :, 2]),
        ], axis=1)
        # mode_target has shape (mode_count, 2).  The cone transform maps the
        # target's two tangential components to three inequalities.
        H = np.asarray(diagnostic["H"], dtype=float)
        columns[index * 3:(index + 1) * 3] = H @ mode_target.T
        compact_target = raw_target + mode_target[-1]
        offsets[index * 3:(index + 1) * 3] = H @ compact_target - MARGIN
        metadata.append({
            "index": int(index),
            "kind": diagnostic["kind"],
            "label": diagnostic["label"],
            "physical_center": diagnostic["physical_center"],
            "H": H.tolist(),
            "raw_target": raw_target.tolist(),
            "compact_baseline_target": compact_target.tolist(),
            "minimum_saved_margin": float(np.min(H @ raw_target - MARGIN)),
        })
    known = np.zeros(mode_count, dtype=float)
    known[-1] = -1.0
    return {
        "A": columns.reshape(len(diagnostics) * 3, mode_count),
        "b": offsets,
        "known_point": known,
        "known_residual": columns.reshape(len(diagnostics) * 3, mode_count) @ known + offsets,
        "diagnostics": metadata,
        "node_count": len(diagnostics),
        "order": ASSEMBLY_ORDER,
    }


def assemble_volume_rows(dynamic, baseline, modes, points, weights, tau):
    """Build the complete vector volume objective columns."""

    hspace = 5.0e-4 * np.sqrt(baseline.nu * tau)
    htime = 1.0e-4 * tau
    baseline_jet = jets(baseline, points, tau, hspace, htime)
    raw_jet = jets(dynamic, points, tau, hspace, htime)
    baseline_residual = momentum(baseline_jet)
    raw_residual = momentum(raw_jet)
    columns = []
    for mode in modes:
        if isinstance(mode, CompactPressureDirection):
            columns.append(baseline_residual - raw_residual)
        else:
            columns.append(momentum(jets(mode, points, tau, hspace, htime)))
    columns = np.stack(columns, axis=-1)
    return {
        "points": points,
        "weights": weights,
        "baseline_residual": baseline_residual,
        "raw_residual": raw_residual,
        "columns": columns,
        "hspace": hspace,
        "htime": htime,
    }


def constrained_volume_fit(volume, equality, cones):
    """Solve the scaled convex volume least-squares problem from known x0."""

    columns = volume["columns"]
    baseline_residual = volume["baseline_residual"]
    weights = volume["weights"]
    mode_count = columns.shape[-1]
    sqrt_weights = np.sqrt(weights)
    C = (columns * sqrt_weights[:, None, None]).reshape(-1, mode_count)
    r = (baseline_residual * sqrt_weights[:, None]).reshape(-1)
    scales = np.maximum(np.linalg.norm(C, axis=0), 1.0e-30)
    C_scaled = C / scales[None, :]
    E = np.asarray(equality["E"], dtype=float)
    m = np.asarray(equality["m"], dtype=float)
    A = np.asarray(cones["A"], dtype=float)
    b = np.asarray(cones["b"], dtype=float)
    x0 = np.asarray(equality["known_point"], dtype=float)
    z0 = scales * x0
    # Normalize both the residual target and the optimization coordinates.
    # The compact pressure coefficient has a very large physical column; raw
    # SLSQP coordinates otherwise mix O(1e5) residuals with O(1e5) variables.
    variable_scale = max(float(np.linalg.norm(z0)), 1.0)
    objective_scale = max(float(np.linalg.norm(r)),
                          float(np.linalg.norm(C_scaled @ z0)), 1.0)
    y0 = z0 / variable_scale
    C_objective = C_scaled * variable_scale / objective_scale
    r_objective = r / objective_scale
    eq_scale = np.maximum(np.max(np.abs(E / scales[None, :]), axis=1), 1.0e-12)
    E_scaled = (E / scales[None, :]) * variable_scale / eq_scale[:, None]
    m_scaled = m / eq_scale
    cone_scale = np.maximum(np.linalg.norm(A / scales[None, :], axis=1),
                            np.abs(b))
    cone_scale = np.maximum(cone_scale, 1.0e-12)
    A_scaled = (A / scales[None, :]) * variable_scale / cone_scale[:, None]
    b_scaled = b / cone_scale

    def objective(z):
        residual = r_objective + C_objective @ z
        # A tiny coefficient-centered ridge makes the selected solution
        # reproducible in null directions without changing the residual gate.
        dz = z - y0
        return 0.5 * float(residual @ residual) + 1.0e-10 * float(dz @ dz)

    def gradient(z):
        residual = r_objective + C_objective @ z
        return C_objective.T @ residual + 2.0e-10 * (z - y0)

    constraints = [
        LinearConstraint(E_scaled, -m_scaled, -m_scaled),
        LinearConstraint(A_scaled, -b_scaled, np.full_like(b_scaled, np.inf)),
    ]
    result = minimize(
        objective,
        y0,
        jac=gradient,
        method="SLSQP",
        constraints=constraints,
        options={"maxiter": 600, "ftol": 1.0e-12, "disp": False},
    )
    y = np.asarray(result.x if np.isfinite(result.x).all() else y0, dtype=float)
    z = variable_scale * y
    x = z / scales
    eq_residual = E @ x + m
    cone_margin = A @ x + b
    feasible = (np.isfinite(x).all() and np.max(np.abs(eq_residual)) <= 1.0e-7
                and np.min(cone_margin) >= -1.0e-7)
    solver_objective = float(objective(y))
    known_objective = float(objective(y0))
    objective_tolerance = 1.0e-10 * max(1.0, abs(known_objective))
    worse_than_known = bool(solver_objective > known_objective + objective_tolerance)
    selected_solver_point = bool(feasible and not worse_than_known)
    if not selected_solver_point:
        x = x0
        eq_residual = E @ x + m
        cone_margin = A @ x + b
        y_selected = y0
    else:
        y_selected = y
    return x, {
        "success": bool(result.success),
        "message": str(result.message),
        "status": int(result.status),
        "iterations": int(getattr(result, "nit", -1)),
        "feasible_solver_point": bool(feasible),
        "used_solver_point": selected_solver_point,
        "fallback_to_known_feasible": bool(not selected_solver_point),
        "selection_guard_triggered": bool(worse_than_known and feasible),
        "worse_than_known": worse_than_known,
        "equality_residual_max": float(np.max(np.abs(eq_residual))),
        "minimum_cone_margin": float(np.min(cone_margin)),
        "objective_at_solution": float(objective(y_selected)),
        "solver_objective": solver_objective,
        "known_feasible_objective": known_objective,
        "objective_tolerance": objective_tolerance,
        "column_scales": scales.tolist(),
        "equality_row_scales": eq_scale.tolist(),
        "cone_row_scales": cone_scale.tolist(),
        "variable_scale": variable_scale,
        "objective_scale": objective_scale,
    }


def _metrics_without_residual(field, points, weights, tau):
    result = evaluate_field(field, points, tau, weights)
    result.pop("residual_components", None)
    return result


def load_saved_field(path=OUTPUT_PATH):
    """Reconstruct the constrained candidate from saved coefficients."""

    report = json.loads(Path(path).read_text(encoding="utf-8"))
    if "coefficients" not in report:
        raise ValueError("Constrained experiment has no coefficients")
    dynamic, dynamic_report = load_dynamic_field(DYNAMIC_PATH)
    baseline = build_compact_baseline(dynamic, dynamic_report,
                                      order=int(report.get("compact_pressure_order", COMPACT_ORDER)))
    modes, _ = build_modes_with_parameters(dynamic, dynamic_report,
                                           float(report["k"]),
                                           compact_pressure=baseline)
    from broad_meridional_momentum import CorrectedField
    return CorrectedField(baseline, modes, report["coefficients"]), report


def run():
    started = time.perf_counter()
    dynamic, dynamic_report = load_dynamic_field(DYNAMIC_PATH)
    k0 = float(dynamic_report["k"])
    tau0 = float(dynamic_report["tau"])
    baseline = build_compact_baseline(dynamic, dynamic_report, order=COMPACT_ORDER)
    modes, mode_names = build_modes_with_parameters(
        dynamic, dynamic_report, k0, compact_pressure=baseline
    )
    mode_count = len(modes)
    radial_breaks = sorted(set(float(value) for value in dynamic_report["radial_breaks"])
                           | set(RADIAL_QUADRATURE_BREAKS))
    locations, labels, _ = _locations(dynamic_report)
    train_points, train_weights = quadrature_nodes(
        dynamic.inner, tau0, 3, 5, 1,
        radial_bounds=(0.01, 0.99), radial_breaks=radial_breaks,
    )
    holdout_points, holdout_weights = quadrature_nodes(
        dynamic.inner, tau0, 4, 4, 1,
        radial_bounds=(0.01, 0.99), radial_breaks=radial_breaks,
        angle_shift=np.pi / 7.0,
    )
    report = {
        "accepted": False,
        "pde_validated": False,
        "scale_recursion_established": False,
        "status": "assembling_constrained_matrices",
        "source_dynamic_report": DYNAMIC_PATH.name,
        "source_unconstrained_report": UNCONSTRAINED_PATH.name,
        "field_loader": "load_saved_field",
        "k": k0,
        "tau": tau0,
        "physical_time": "t=-tau",
        "dk_dt": float(1.0 / (tau0 * np.log(2.0))),
        "nu": float(dynamic.nu),
        "mode_count": mode_count,
        "mode_names": mode_names,
        "known_feasible_point": {"compact_pressure_coefficient": -1.0,
                                  "all_other_coefficients": 0.0},
        "compatibility_source": {
            "saved_control_count": 19,
            "saved_moment_equation_count": 4,
            "saved_cone_inequality_count": 81,
            "cone_node_count": 27,
            "cone_margin": MARGIN.tolist(),
            "reuse": "raw dynamic source_T and H/unit_target normalization from broad_shear_dynamic_control.json; new mode targets evaluated from actual jets",
        },
        "quadrature": {
            "eta_bounds": [-0.5, 0.5],
            "radial_fraction_bounds": [0.01, 0.99],
            "radial_split_breaks": list(radial_breaks),
            "training": {"eta_order": 3, "radial_order_per_interval": 5,
                         "angle_order": 1, "point_count": int(len(train_points))},
            "holdout": {"eta_order": 4, "radial_order_per_interval": 4,
                        "angle_order": 1, "angle_shift": float(np.pi / 7.0),
                        "point_count": int(len(holdout_points))},
            "jacobian": "2*pi*r*(dr/dy)*(dz/deta), theta order 1 carries exact 2*pi for axisymmetric fields",
        },
        "residual_definition": "R=partial_t u+(u dot grad)u+grad p-nu Delta u; fourth-order Cartesian FD jets; t=-tau",
    }
    _save(report)

    print(json.dumps({"stage": "moment_rows_started", "mode_count": mode_count}), flush=True)
    equality = assemble_moment_rows(dynamic, dynamic_report, baseline, modes, k0,
                                    radial_breaks)
    report["moment_problem"] = {
        "E": equality["E"].tolist(),
        "m": equality["m"].tolist(),
        "raw_moment": equality["raw_moment"].tolist(),
        "known_point": equality["known_point"].tolist(),
        "known_residual": equality["known_residual"].tolist(),
        "panel_count": equality["panel_count"],
        "quadrature_order": equality["quadrature_order"],
    }
    report["status"] = "moment_rows_assembled"
    _save(report)
    print(json.dumps({"stage": "cone_rows_started", "location_count": len(locations)}), flush=True)
    cones = assemble_cone_rows(dynamic, dynamic_report, modes, k0, radial_breaks)
    report["cone_problem"] = {
        "A": cones["A"].tolist(),
        "b": cones["b"].tolist(),
        "known_point": cones["known_point"].tolist(),
        "known_residual": cones["known_residual"].tolist(),
        "node_count": cones["node_count"],
        "quadrature_order": cones["order"],
        "diagnostics": cones["diagnostics"],
    }
    report["status"] = "moment_cone_matrices_assembled"
    _save(report)

    print(json.dumps({"stage": "volume_rows_started", "train_points": len(train_points)}), flush=True)
    volume = assemble_volume_rows(dynamic, baseline, modes,
                                  train_points, train_weights, tau0)
    # Save all equality/inequality matrices before any optimizer call.
    report["assembled_problem"] = {
        "mode_names": mode_names,
        "moment_shape": list(equality["E"].shape),
        "cone_shape": list(cones["A"].shape),
        "volume_column_shape": list(volume["columns"].shape),
        "volume_point_count": len(train_points),
        "volume_weights": train_weights.tolist(),
        "volume_baseline_residual": volume["baseline_residual"].tolist(),
        "volume_columns": volume["columns"].tolist(),
        "objective": "weighted physical-volume squared norm of complete vector residual",
    }
    report["status"] = "matrices_assembled"
    report["elapsed_before_solve_seconds"] = float(time.perf_counter() - started)
    _save(report)
    print(json.dumps({"stage": "matrices_saved", "moment_shape": list(equality["E"].shape),
                      "cone_shape": list(cones["A"].shape),
                      "volume_shape": list(volume["columns"].shape),
                      "known_moment_residual": float(np.max(np.abs(equality["known_residual"]))),
                      "known_cone_margin": float(np.min(cones["known_residual"]))}), flush=True)

    coefficients, solve = constrained_volume_fit(volume, equality, cones)
    report["coefficients"] = coefficients.tolist()
    report["solve"] = solve
    report["status"] = "constrained_fit_replayed"
    # Checkpoint the coefficients before the expensive independent replay.
    report["elapsed_after_solve_seconds"] = float(time.perf_counter() - started)
    _save(report)
    # Direct full-FD replay on train and independently sampled shifted nodes.
    from broad_meridional_momentum import CorrectedField
    candidate = CorrectedField(baseline, modes, coefficients)
    report["baseline_raw"] = {
        "training": _metrics_without_residual(dynamic, train_points, train_weights, tau0),
        "holdout": _metrics_without_residual(dynamic, holdout_points, holdout_weights, tau0),
    }
    report["baseline_compact"] = {
        "training": _metrics_without_residual(baseline, train_points, train_weights, tau0),
        "holdout": _metrics_without_residual(baseline, holdout_points, holdout_weights, tau0),
    }
    report["corrected"] = {
        "training": _metrics_without_residual(candidate, train_points, train_weights, tau0),
        "holdout": _metrics_without_residual(candidate, holdout_points, holdout_weights, tau0),
    }
    # Independent radial orders replay the saved compatibility diagnostics.
    replay = replay_moments(candidate, dynamic.inner, k0, REPLAY_ORDER,
                            0.002, radial_breaks=radial_breaks)
    cone_replay = outer_cones(
        candidate, k0, order=REPLAY_CONE_ORDER,
        radial_breaks=radial_breaks, locations=locations,
    )
    report["compatibility_replay"] = {
        "moment_order": REPLAY_ORDER,
        "moments": np.asarray(replay, dtype=float).tolist(),
        "moment_max_abs": float(np.max(np.abs(replay))),
        "cone_order": REPLAY_CONE_ORDER,
        "cone_rows": cone_replay,
        "cone_pass_count": int(sum(row["cone_pass"] for row in cone_replay)),
        "cone_count": len(cone_replay),
    }
    report["scope"] = (
        "Instantaneous k0 constrained complete-vector volume fit with four saved moment equations and 81 sampled cone inequalities. "
        "Independent split-quadrature FD holdout and higher-order moment/cone replay are reported. "
        "No continuum, finite-time, forcing, PDE or scale-recursion claim."
    )
    report["elapsed_seconds"] = float(time.perf_counter() - started)
    _save(report)
    print(json.dumps({
        "stage": "constrained_fit_complete",
        "status": report["status"],
        "solver_success": solve["success"],
        "used_solver_point": solve["used_solver_point"],
        "minimum_cone_margin": solve["minimum_cone_margin"],
        "moment_max_abs_replay": report["compatibility_replay"]["moment_max_abs"],
        "cone_pass_count": report["compatibility_replay"]["cone_pass_count"],
        "raw_holdout_max": report["baseline_raw"]["holdout"]["momentum_max"],
        "corrected_holdout_max": report["corrected"]["holdout"]["momentum_max"],
        "raw_holdout_l2": report["baseline_raw"]["holdout"]["momentum_volume_L2"],
        "corrected_holdout_l2": report["corrected"]["holdout"]["momentum_volume_L2"],
        "elapsed_seconds": report["elapsed_seconds"],
    }), flush=True)
    return report


if __name__ == "__main__":
    run()
