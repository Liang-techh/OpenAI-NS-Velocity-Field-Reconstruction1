"""Joint instantaneous mean-control and exact-curl wave-energy fit.

The saved 44-column constrained mean problem is reused at ``k0``.  One
additional column is the complete angular-mean force of the saved growing
mode-1 exact-curl wave, evaluated by ``wave_mean_flux.single_mode_force``.
Its coefficient ``e`` is constrained to be nonnegative and represents the
wave energy/amplitude squared relative to the saved potential coefficient
vector.  The four integrated moment rows and 81 cone rows receive matching
wave columns from the same radial grouped quadrature and saved cone ``H``.

The selected result corrects the angular mean residual only.  Oscillatory
wave momentum, continuum constraints, finite-time evolution, PDE validity,
and scale recursion remain outside scope.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import LinearConstraint, minimize


ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from affine_momentum import jets, momentum  # noqa: E402
from broad_meridional_constrained import load_saved_field as load_mean_candidate  # noqa: E402
from broad_meridional_momentum import (  # noqa: E402
    CorrectedField,
    quadrature_nodes,
)
from broad_shear_dynamic_control import load_saved_field as load_dynamic_field  # noqa: E402
from broad_wave_mean_fit import patch_nodes  # noqa: E402
from grouped_outer_cache import _group_indices_by_z, _group_quadrature  # noqa: E402
from meridional_state_cache import value_and_pressure_modes  # noqa: E402
from midplane_integrated_moment_balance import evaluate as replay_moments  # noqa: E402
from supported_fourier_basis import basis_data  # noqa: E402
from wave_mean_flux import single_mode_force  # noqa: E402


MEAN_PATH = ROOT / "broad_meridional_constrained.json"
GROWTH_PATH = ROOT / "broad_shear_growth.json"
DYNAMIC_PATH = ROOT / "broad_shear_dynamic_control.json"
OUTPUT_PATH = ROOT / "joint_wave_mean_fit.json"
MOMENT_ORDER_SOURCE = 48
MOMENT_ORDER_REPLAY = 96
MOMENT_ZFACTOR = 0.002
CONSTRAINT_ORDER_SOURCE = 48
HOLDOUT_ETA_ORDER = 4
HOLDOUT_RADIAL_ORDER = 4
HOLDOUT_ANGLE_SHIFT = np.pi / 7.0
MARGIN = np.array([0.02, 0.0, 0.0], dtype=float)


def _save(report):
    OUTPUT_PATH.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")


def _array(value):
    return np.asarray(value, dtype=float)


def _list(value):
    return _array(value).tolist()


def _decode_complex(value):
    data = np.asarray(value, dtype=float)
    if data.ndim != 2 or data.shape[1] != 2:
        raise ValueError("wave coefficients must be packed as [real, imaginary]")
    return data[:, 0] + 1j * data[:, 1]


def _wave_metadata(path):
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    if "selected_candidate" in data:
        selected = data["selected_candidate"]
        return {
            "mode": int(selected["mode"]),
            "center": data["center"],
            "widths": data["widths"],
            "degree": int(data["degree"]),
            "carrier": data["carriers"][str(selected["mode"])],
            "potential_coefficients": selected["potential_coefficients"],
            "growth_amplitude": selected.get("amplitude"),
            "source": Path(path).name,
        }
    # The co-designed wave report stores the selected complex potential in
    # ``selected.coefficients_original`` while keeping the geometry and
    # carrier at the report top level.  Accept that report as an optional
    # future wave input without changing the default saved-growth mode.
    if ("selected" in data
            and "coefficients_original" in data["selected"]):
        selected = data["selected"]
        return {
            "mode": int(data["mode"]),
            "center": data["center"],
            "widths": data["widths"],
            "degree": int(data["degree"]),
            "carrier": data["carrier"],
            "potential_coefficients": selected["coefficients_original"],
            "growth_amplitude": selected.get("growth_lambda"),
            "source": Path(path).name,
        }
    if "wave_field" in data:
        field = data["wave_field"]
        return {
            "mode": int(field["mode"]),
            "center": field["center"],
            "widths": field["widths"],
            "degree": int(field["degree"]),
            "carrier": field["carrier"],
            "potential_coefficients": field["potential_coefficients"],
            "growth_amplitude": field.get("growth_amplitude_metadata"),
            "source": Path(path).name,
        }
    raise ValueError("wave source has neither selected_candidate nor wave_field")


def _wave_force_cartesian(points, wave, nu, h):
    """Evaluate the axisymmetric angular-mean force in Cartesian coordinates."""

    points = np.asarray(points, dtype=float)
    radius = np.hypot(points[:, 0], points[:, 1])
    rz = np.column_stack((radius, points[:, 2]))
    force_cyl = single_mode_force(
        rz,
        wave["center"],
        wave["widths"],
        wave["mode"],
        wave["degree"],
        wave["carrier"],
        wave["coefficients"],
        nu,
        h,
    )
    cos_theta = np.divide(points[:, 0], radius,
                          out=np.ones_like(radius), where=radius > 0.0)
    sin_theta = np.divide(points[:, 1], radius,
                          out=np.zeros_like(radius), where=radius > 0.0)
    return np.column_stack((
        cos_theta * force_cyl[:, 0] - sin_theta * force_cyl[:, 1],
        sin_theta * force_cyl[:, 0] + cos_theta * force_cyl[:, 1],
        force_cyl[:, 2],
    ))


def _radial_panels(inner, k, eta_values, breaks, order):
    """Build the fixed-radius moment panels without reassembling mean jets."""

    nodes, weights = leggauss(int(order))
    tau = 0.5 * 2.0 ** (-float(k))
    panels = []
    for eta in eta_values:
        ends = inner.from_similarity(
            inner.p.X_max * np.array([1.0, 16.0 ** 2]),
            [eta, eta],
            tau,
        )
        ri, ro = map(float, ends[:, 0])
        z = float(ends[0, 2])
        edges = [0.0, ri, ro]
        edges.extend(ri + (ro - ri) * float(value) for value in breaks)
        edges = sorted(set(edges))
        radii = np.concatenate([
            0.5 * (lo + hi) + 0.5 * (hi - lo) * nodes
            for lo, hi in zip(edges[:-1], edges[1:])
        ])
        radial_weights = np.concatenate([
            0.5 * (hi - lo) * weights
            for lo, hi in zip(edges[:-1], edges[1:])
        ])
        panels.append({
            "eta": float(eta),
            "z": z,
            "outer_radius": ro,
            "radii": radii,
            "weights": radial_weights,
        })
    return panels


def _wave_moment_column(inner, k, breaks, wave, nu, h, order):
    panels = _radial_panels(inner, k, (-0.2, 0.2), breaks, order)
    values = []
    for panel in panels:
        points = np.column_stack((panel["radii"],
                                  np.zeros(len(panel["radii"])),
                                  np.full(len(panel["radii"]), panel["z"])))
        force = _wave_force_cartesian(points, wave, nu, h)
        rr = panel["radii"]
        ww = panel["weights"]
        ro = panel["outer_radius"]
        values.extend((
            -np.dot(ww * rr ** 2, force[:, 1]) / ro ** 2,
            -np.dot(ww * rr, force[:, 2]) / ro,
        ))
    return np.asarray(values), panels


def _grouped_cone_panels(saved, inner, k, breaks, order):
    diagnostics = saved["cone_problem"]["diagnostics"]
    centers = [np.asarray(row["physical_center"], dtype=float)
               for row in diagnostics]
    points = [center.tolist() for center in centers]
    panels_by_index = {}
    for group in _group_indices_by_z(centers):
        panels_by_index.update(_group_quadrature(
            centers, group, inner, 0.5 * 2.0 ** (-float(k)), int(order),
            points, breaks,
        ))
    panels = [panels_by_index[index] for index in range(len(centers))]
    return np.asarray(points, dtype=float), panels


def _wave_cone_column(saved, inner, k, breaks, wave, nu, h, order):
    points, panels = _grouped_cone_panels(saved, inner, k, breaks, order)
    force = _wave_force_cartesian(points, wave, nu, h)
    values = []
    diagnostics = saved["cone_problem"]["diagnostics"]
    for index, (sl, radii, weights, radius) in enumerate(panels):
        target = np.array([
            -np.dot(weights * radii ** 2, force[sl, 1]) / radius ** 2,
            -np.dot(weights * radii, force[sl, 2]) / radius,
        ])
        H = np.asarray(diagnostics[index]["H"], dtype=float)
        values.extend((H @ target).tolist())
    return np.asarray(values), {
        "order": int(order),
        "point_count": int(len(points)),
        "location_count": len(panels),
        "group_count": len(_group_indices_by_z([
            np.asarray(row["physical_center"], dtype=float)
            for row in diagnostics
        ])),
    }


def _constraint_check(E, m, A, b, x):
    equality = _array(E) @ _array(x) + _array(m)
    cone = _array(A) @ _array(x) + _array(b)
    return {
        "equality_residual_max": float(np.max(np.abs(equality))),
        "minimum_cone_margin": float(np.min(cone)),
        "cone_pass_count": int(np.sum(cone >= -1.0e-7)),
        "cone_count": int(len(cone)),
        "valid": bool(np.isfinite(x).all()
                      and np.max(np.abs(equality)) <= 1.0e-7
                      and np.min(cone) >= -1.0e-7),
        "equality_residual": _list(equality),
    }


def _joint_volume_fit(volume, equality, cones, known_point):
    """Scaled constrained LS with e >= 0 and a known feasible guard point."""

    columns = _array(volume["columns"])
    baseline = _array(volume["baseline_residual"])
    weights = _array(volume["weights"])
    E = _array(equality["E"])
    m = _array(equality["m"])
    A = _array(cones["A"])
    b = _array(cones["b"])
    x0 = _array(known_point)
    count = columns.shape[-1]
    sqrt_weights = np.sqrt(weights)
    C = (columns * sqrt_weights[:, None, None]).reshape(-1, count)
    r = (baseline * sqrt_weights[:, None]).reshape(-1)
    scales = np.maximum(np.linalg.norm(C, axis=0), 1.0e-30)
    C_scaled = C / scales[None, :]
    z0 = scales * x0
    variable_scale = max(float(np.linalg.norm(z0)), 1.0)
    objective_scale = max(float(np.linalg.norm(r)),
                          float(np.linalg.norm(C_scaled @ z0)), 1.0)
    y0 = z0 / variable_scale
    C_objective = C_scaled * variable_scale / objective_scale
    r_objective = r / objective_scale
    eq_scale = np.maximum(np.max(np.abs(E / scales[None, :]), axis=1), 1.0e-12)
    E_scaled = (E / scales[None, :]) * variable_scale / eq_scale[:, None]
    m_scaled = m / eq_scale
    cone_scale = np.maximum(np.linalg.norm(A / scales[None, :], axis=1), np.abs(b))
    cone_scale = np.maximum(cone_scale, 1.0e-12)
    A_scaled = (A / scales[None, :]) * variable_scale / cone_scale[:, None]
    b_scaled = b / cone_scale

    def objective(y):
        residual = r_objective + C_objective @ y
        dz = y - y0
        return 0.5 * float(residual @ residual) + 1.0e-10 * float(dz @ dz)

    def gradient(y):
        residual = r_objective + C_objective @ y
        return C_objective.T @ residual + 2.0e-10 * (y - y0)

    constraints = [
        LinearConstraint(E_scaled, -m_scaled, -m_scaled),
        LinearConstraint(A_scaled, -b_scaled, np.full_like(b_scaled, np.inf)),
    ]
    lower = np.full(count, -np.inf)
    lower[-1] = 0.0
    result = minimize(
        objective,
        y0,
        jac=gradient,
        method="SLSQP",
        bounds=list(zip(lower, np.full(count, np.inf))),
        constraints=constraints,
        options={"maxiter": 600, "ftol": 1.0e-12, "disp": False},
    )
    y_solver = np.asarray(result.x if np.isfinite(result.x).all() else y0)
    solver_x = (variable_scale * y_solver) / scales
    solver_check = _constraint_check(E, m, A, b, solver_x)
    solver_check["energy_nonnegative"] = bool(solver_x[-1] >= -1.0e-9)
    solver_objective = float(objective(y_solver))
    known_objective = float(objective(y0))
    tolerance = 1.0e-10 * max(1.0, abs(known_objective))
    worse = bool(solver_objective > known_objective + tolerance)
    use_solver = bool(solver_check["valid"] and solver_check["energy_nonnegative"]
                      and not worse)
    if use_solver:
        selected = solver_x
        selected_y = y_solver
    else:
        selected = x0.copy()
        selected_y = y0
    selected_check = _constraint_check(E, m, A, b, selected)
    selected_check["energy_nonnegative"] = bool(selected[-1] >= -1.0e-9)
    return selected, {
        "success": bool(result.success),
        "status": int(result.status),
        "message": str(result.message),
        "iterations": int(getattr(result, "nit", -1)),
        "feasible_solver_point": bool(solver_check["valid"]
                                       and solver_check["energy_nonnegative"]),
        "used_solver_point": use_solver,
        "fallback_to_known_feasible": bool(not use_solver),
        "worse_than_known": worse,
        "selection_guard_triggered": bool(worse),
        "solver_objective": solver_objective,
        "known_feasible_objective": known_objective,
        "objective_at_solution": float(objective(selected_y)),
        "objective_tolerance": tolerance,
        "variable_scale": variable_scale,
        "objective_scale": objective_scale,
        "column_scales": _list(scales),
        "solver_constraint_check": solver_check,
        "selected_constraint_check": selected_check,
    }


def _metrics_with_wave(field, points, weights, tau, wave, energy, nu, h):
    point_jet = jets(field, points, tau, h, 1.0e-4 * tau)
    residual = momentum(point_jet)
    force = _wave_force_cartesian(points, wave, nu, h)
    total = residual + float(energy) * force
    divergence = np.trace(point_jet[1], axis1=1, axis2=2)
    norms = np.linalg.norm(total, axis=1)
    weighted = np.sum(weights[:, None] * total * total)
    return {
        "point_count": int(len(points)),
        "momentum_max": float(np.max(norms)),
        "momentum_volume_L2": float(np.sqrt(weighted)),
        "momentum_volume_RMS": float(np.sqrt(weighted / np.sum(weights))),
        "divergence_max_mean_only": float(np.max(np.abs(divergence))),
        "divergence_volume_L2_mean_only": float(
            np.sqrt(np.sum(weights * divergence * divergence))
        ),
    }


def _mean_only_metrics(field, points, weights, tau):
    point_jet = jets(field, points, tau, 5.0e-4 * np.sqrt(field.nu * tau),
                     1.0e-4 * tau)
    residual = momentum(point_jet)
    divergence = np.trace(point_jet[1], axis1=1, axis2=2)
    norms = np.linalg.norm(residual, axis=1)
    weighted = np.sum(weights[:, None] * residual * residual)
    return {
        "point_count": int(len(points)),
        "momentum_max": float(np.max(norms)),
        "momentum_volume_L2": float(np.sqrt(weighted)),
        "momentum_volume_RMS": float(np.sqrt(weighted / np.sum(weights))),
        "divergence_max": float(np.max(np.abs(divergence))),
        "divergence_volume_L2": float(
            np.sqrt(np.sum(weights * divergence * divergence))
        ),
    }


def run(wave_path=GROWTH_PATH):
    started = time.perf_counter()
    saved = json.loads(MEAN_PATH.read_text(encoding="utf-8"))
    dynamic, dynamic_report = load_dynamic_field(DYNAMIC_PATH)
    k = float(saved["k"])
    tau = float(saved["tau"])
    nu = float(dynamic.nu)
    h = 5.0e-4 * np.sqrt(nu * tau)
    wave = _wave_metadata(wave_path)
    wave["coefficients"] = _decode_complex(wave["potential_coefficients"])
    if wave["mode"] != 1:
        raise ValueError("joint fit requires the saved growing mode 1")
    assembled = saved["assembled_problem"]
    baseline_residual = _array(assembled["volume_baseline_residual"])
    mean_columns = _array(assembled["volume_columns"])
    volume_weights = _array(assembled["volume_weights"])
    breaks = [float(value) for value in saved["quadrature"]["radial_split_breaks"]]
    train_points, train_weights = quadrature_nodes(
        dynamic.inner,
        tau,
        3,
        5,
        1,
        radial_bounds=(0.01, 0.99),
        radial_breaks=breaks,
    )
    if len(train_points) != len(volume_weights):
        raise ValueError("saved volume matrix and reconstructed train grid differ")
    weight_difference = float(np.max(np.abs(train_weights - volume_weights)))
    wave_volume_force = _wave_force_cartesian(train_points, wave, nu, h)
    volume_columns = np.concatenate((mean_columns, wave_volume_force[:, :, None]), axis=2)

    mean_E = _array(saved["moment_problem"]["E"])
    mean_m = _array(saved["moment_problem"]["m"])
    cone_A = _array(saved["cone_problem"]["A"])
    cone_b = _array(saved["cone_problem"]["b"])
    wave_moment, moment_panels = _wave_moment_column(
        dynamic.inner, k, breaks, wave, nu, h, MOMENT_ORDER_SOURCE
    )
    wave_moment_high, _ = _wave_moment_column(
        dynamic.inner, k, breaks, wave, nu, h, MOMENT_ORDER_REPLAY
    )
    wave_cone, cone_geometry = _wave_cone_column(
        saved, dynamic.inner, k, breaks, wave, nu, h, CONSTRAINT_ORDER_SOURCE
    )
    equality_E = np.column_stack((mean_E, wave_moment))
    cone_A_joint = np.column_stack((cone_A, wave_cone))
    known_point = np.r_[np.asarray(saved["coefficients"], dtype=float), 0.0]
    known_check = _constraint_check(equality_E, mean_m, cone_A_joint, cone_b,
                                    known_point)
    report = {
        "accepted": False,
        "pde_validated": False,
        "scale_recursion_established": False,
        "status": "matrices_assembled",
        "k": k,
        "tau": tau,
        "nu": nu,
        "mean_source": MEAN_PATH.name,
        "wave_source": wave["source"],
        "wave_field": {
            key: (value.tolist() if isinstance(value, np.ndarray) else value)
            for key, value in wave.items() if key != "coefficients"
        },
        "wave_coefficient_l2_norm": float(np.linalg.norm(wave["coefficients"])),
        "wave_energy_definition": "e >= 0 multiplies the angular-mean complete single_mode_force of the saved coefficient vector; amplitude = sqrt(e)",
        "instantaneous_geometry": "mean velocity is unchanged at k0 because all saved velocity directions vanish instantaneously; saved H is reused consistently for cone columns",
        "quadrature": {
            "volume_train": {
                "eta_order": 3,
                "radial_order_per_split_interval": 5,
                "point_count": int(len(train_points)),
                "physical_volume": float(np.sum(volume_weights)),
                "weight_reconstruction_max_abs_difference": weight_difference,
                "jacobian": "2*pi*r*(dr/dy)*(dz/deta); saved mapped physical weights reused",
            },
            "moment_source_order": MOMENT_ORDER_SOURCE,
            "moment_replay_order": MOMENT_ORDER_REPLAY,
            "cone_source_order": CONSTRAINT_ORDER_SOURCE,
            "holdout": {
                "eta_order": HOLDOUT_ETA_ORDER,
                "radial_order_per_split_interval": HOLDOUT_RADIAL_ORDER,
                "angle_shift": float(HOLDOUT_ANGLE_SHIFT),
            },
        },
        "assembled_problem": {
            "volume_column_shape": list(volume_columns.shape),
            "volume_point_count": int(len(train_points)),
            "volume_weights": volume_weights.tolist(),
            "volume_baseline_residual": baseline_residual.tolist(),
            "volume_columns": volume_columns.tolist(),
            "moment_E": equality_E.tolist(),
            "moment_m": mean_m.tolist(),
            "cone_A": cone_A_joint.tolist(),
            "cone_b": cone_b.tolist(),
            "wave_volume_force_column": wave_volume_force.tolist(),
            "wave_moment_column_source_order": wave_moment.tolist(),
            "wave_cone_column_source_order": wave_cone.tolist(),
        },
        "known_feasible_point": {
            "coefficients": known_point.tolist(),
            "mean_coefficients_source": "broad_meridional_constrained.json coefficients",
            "wave_energy": 0.0,
            "constraint_check": known_check,
        },
        "wave_moment_source_order_check": {
            "order_source": MOMENT_ORDER_SOURCE,
            "order_replay": MOMENT_ORDER_REPLAY,
            "source_column": wave_moment.tolist(),
            "replay_column": wave_moment_high.tolist(),
            "maximum_absolute_difference": float(np.max(np.abs(wave_moment_high - wave_moment))),
            "relative_l2_difference": float(
                np.linalg.norm(wave_moment_high - wave_moment)
                / max(np.linalg.norm(wave_moment_high), np.finfo(float).tiny)
            ),
            "moment_panel_count": len(moment_panels),
        },
        "cone_geometry": cone_geometry,
        "stages": [],
        "scope": "Instantaneous joint angular-mean fit: 44 saved mean controls plus nonnegative exact-curl mode-1 wave energy. Spatial metrics are mean residual plus e times single_mode_force; oscillatory wave momentum is unchecked. No PDE, finite-time, continuum, or scale-recursion acceptance.",
    }
    _save(report)
    print(json.dumps({"stage": "matrices_saved",
                      "volume_shape": list(volume_columns.shape),
                      "moment_shape": list(equality_E.shape),
                      "cone_shape": list(cone_A_joint.shape),
                      "known_feasible": known_check}), flush=True)

    coefficients, solve = _joint_volume_fit(
        {"columns": volume_columns,
         "baseline_residual": baseline_residual,
         "weights": volume_weights},
        {"E": equality_E, "m": mean_m},
        {"A": cone_A_joint, "b": cone_b},
        known_point,
    )
    report["coefficients"] = coefficients.tolist()
    report["mean_coefficients"] = coefficients[:-1].tolist()
    report["wave_energy"] = float(coefficients[-1])
    report["wave_amplitude"] = float(np.sqrt(max(coefficients[-1], 0.0)))
    report["solve"] = solve
    report["status"] = "fit_complete"
    report["elapsed_after_solve_seconds"] = float(time.perf_counter() - started)
    _save(report)
    print(json.dumps({"stage": "fit_complete",
                      "energy": report["wave_energy"],
                      "amplitude": report["wave_amplitude"],
                      "constraint_check": solve["selected_constraint_check"]}), flush=True)

    mean_field, _ = load_mean_candidate(MEAN_PATH)
    # Reconstruct arbitrary selected mean controls on the same compact baseline.
    dynamic_for_modes, source_for_modes = load_dynamic_field(DYNAMIC_PATH)
    base, values, pressures, _, _ = value_and_pressure_modes(
        dynamic_for_modes, source_for_modes
    )
    selected_mean = CorrectedField(
        base,
        __import__("broad_meridional_momentum").build_modes_with_parameters(
            dynamic_for_modes, source_for_modes, k, compact_pressure=base
        )[0],
        coefficients[:-1],
    )
    # The saved selected mean is retained as a direct baseline; selected_mean
    # is the reconstructed candidate for the joint coefficients.
    holdout_points, holdout_weights = quadrature_nodes(
        dynamic.inner,
        tau,
        HOLDOUT_ETA_ORDER,
        HOLDOUT_RADIAL_ORDER,
        1,
        radial_bounds=(0.01, 0.99),
        radial_breaks=breaks,
        angle_shift=HOLDOUT_ANGLE_SHIFT,
    )
    center = np.asarray(wave["center"], dtype=float)
    widths = np.asarray(wave["widths"], dtype=float)
    patch_breaks = sorted(set(breaks + [0.62]))
    patch_points, patch_weights, patch_geometry = patch_nodes(
        selected_mean,
        center,
        widths,
        tau,
        9,
        9,
        patch_breaks,
    )
    domains = {
        "annulus_holdout": (holdout_points, holdout_weights),
        "wave_patch_756": (patch_points, patch_weights),
    }
    metrics = {}
    for label, (points, weights) in domains.items():
        metrics[label] = {
            "mean_feasible_baseline": _metrics_with_wave(
                mean_field,
                points,
                weights,
                tau,
                wave,
                0.0,
                nu,
                h,
            ),
            "joint_corrected": _metrics_with_wave(
                selected_mean,
                points,
                weights,
                tau,
                wave,
                coefficients[-1],
                nu,
                h,
            ),
            "domain": {
                "point_count": int(len(points)),
                "physical_volume": float(np.sum(weights)),
                "weights": "mapped physical volume Jacobian",
            },
        }
    report["spatial_metrics"] = metrics
    report["patch_geometry"] = patch_geometry
    report["status"] = "spatial_holdout_complete"
    _save(report)

    # Independent high-order integrated moment replay adds the wave force
    # column to the actual mean residual replay; this exposes order mismatch.
    mean_replay = replay_moments(
        selected_mean,
        dynamic.inner,
        k,
        MOMENT_ORDER_REPLAY,
        MOMENT_ZFACTOR,
        radial_breaks=breaks,
    )
    joint_replay = np.asarray(mean_replay) + coefficients[-1] * wave_moment_high
    report["moment_replay"] = {
        "order": MOMENT_ORDER_REPLAY,
        "mean_only": np.asarray(mean_replay, dtype=float).tolist(),
        "joint_corrected": joint_replay.tolist(),
        "joint_max_abs": float(np.max(np.abs(joint_replay))),
        "source_order_joint_rows": (equality_E @ coefficients + mean_m).tolist(),
        "source_vs_replay_max_abs": float(np.max(
            np.abs(joint_replay - (equality_E @ coefficients + mean_m))
        )),
    }
    report["status"] = "completed"
    report["elapsed_seconds"] = float(time.perf_counter() - started)
    _save(report)
    print(json.dumps({
        "stage": "completed",
        "energy": report["wave_energy"],
        "annulus_mean_l2": metrics["annulus_holdout"]["mean_feasible_baseline"]["momentum_volume_L2"],
        "annulus_joint_l2": metrics["annulus_holdout"]["joint_corrected"]["momentum_volume_L2"],
        "patch_mean_l2": metrics["wave_patch_756"]["mean_feasible_baseline"]["momentum_volume_L2"],
        "patch_joint_l2": metrics["wave_patch_756"]["joint_corrected"]["momentum_volume_L2"],
        "joint_moment_max": report["moment_replay"]["joint_max_abs"],
        "elapsed_seconds": report["elapsed_seconds"],
    }), flush=True)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--wave-path", default=str(GROWTH_PATH))
    args = parser.parse_args()
    run(Path(args.wave_path))
