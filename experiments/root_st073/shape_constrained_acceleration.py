"""Constrained endpoint acceleration fit with fixed-cylinder shape gates.

This diagnostic uses the frozen 324-column endpoint momentum linearization and
the analytic acceleration shape oracle.  The acceleration is zero at the
reference time, so the parent state is preserved there by construction.  The
fit is only a frozen-grid tangent solve; it does not replay the nonlinear
endpoint or claim a trajectory/PDE result.
"""

from __future__ import annotations

import hashlib
import json
import os
import time
from pathlib import Path

for _name in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS"):
    os.environ[_name] = "1"

import numpy as np
from scipy.optimize import minimize

from acceleration_shape_oracle import AccelerationShapeOracle
from affine_momentum import jets, momentum
from broad_meridional_constrained import load_saved_field
from endpoint_acceleration_projection import (
    _unpack_wave,
    MixedAffineField,
)
from grouped_joined_field import install_in_field
from supported_fourier_analytic_jets import basis_jets


ROOT = Path(__file__).resolve().parent
REFINED_CACHE = ROOT / "refined_wave_momentum_cache.npz"
PROJECTION = ROOT / "endpoint_acceleration_projection.json"
FROZEN = ROOT / "full_wave_frozen_cache.json"
OUTPUT = ROOT / "shape_constrained_acceleration.json"
MOMENTUM_CACHE = ROOT / "acceleration_momentum_cache.npz"
PRESSURE_SEED = ROOT / "pressure_acceleration_seed.json"

CONTROL_COUNT = 324
RIDGE = 1.0e-6
PEAK_RELATIVE_CAP = 1.0e-6
INITIAL_POOL = 32
MAX_ROUNDS = 1
TRUST_RADIUS = 1.0e-2


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _metric(residual, weights):
    residual = np.asarray(residual, dtype=float)
    weights = np.asarray(weights, dtype=float)
    magnitudes = np.linalg.norm(residual, axis=1)
    return {
        "point_count": int(len(residual)),
        "momentum_max": float(np.max(magnitudes)),
        "momentum_volume_L2": float(np.linalg.norm(
            residual.reshape(-1) * np.repeat(np.sqrt(weights), 3)
        )),
        "momentum_volume_RMS": float(np.sqrt(
            np.sum(weights * np.sum(residual * residual, axis=1)) / np.sum(weights)
        )),
        "physical_volume": float(np.sum(weights)),
    }


def _signed_shape(values, reference, requested):
    values = np.asarray(values, dtype=float)
    reference = np.asarray(reference, dtype=float)
    signs = np.asarray((-1.0, 1.0, np.sign(reference[2])), dtype=float)
    fractional = signs * (values - reference) / np.maximum(np.abs(reference), 1.0e-30)
    return fractional - np.asarray(requested, dtype=float)


def _reconstruct_base(projection, frozen, points, hspace, htime):
    geometry = frozen["inputs"]["wave"]
    center = tuple(float(value) for value in geometry["center"])
    widths = tuple(float(value) for value in geometry["widths"])
    carrier = np.asarray(geometry["carrier"], dtype=float)
    tau0 = float(frozen["inputs"]["mean"]["tau"])
    tau1 = float(projection["inputs"]["tau1"])
    base_control = np.asarray(projection["inputs"]["base_tangent_coefficients"], dtype=float)
    if base_control.shape != (264,):
        raise ValueError(f"Expected 264 base controls, got {base_control.shape}")
    wave = _unpack_wave(projection["inputs"]["initial_wave_coefficients"])
    mean, mean_report = load_saved_field()
    install_in_field(mean)
    base = MixedAffineField(mean, center, widths, carrier, wave, base_control, tau0)
    endpoint_jets = jets(base, points, tau1, hspace, htime)
    residual = momentum(endpoint_jets)
    return {
        "base": base,
        "mean": mean,
        "mean_report": mean_report,
        "residual": np.asarray(residual, dtype=float),
        "velocity": np.asarray(endpoint_jets[0], dtype=float),
        "gradient": np.asarray(endpoint_jets[1], dtype=float),
        "tau0": tau0,
        "tau1": tau1,
        "base_control": base_control,
        "wave": wave,
        "center": center,
        "widths": widths,
        "carrier": carrier,
    }


def _pool_data(y, indices, base_residual, design, scales, coordinate_scale, peak_cap):
    indices = np.asarray(sorted(indices), dtype=int)
    point_design = design.reshape(len(base_residual), 3, CONTROL_COUNT)[indices]
    physical = coordinate_scale * np.asarray(y, dtype=float) / scales
    values = base_residual[indices] + np.einsum("nci,i->nc", point_design, physical)
    normalized_design = coordinate_scale * point_design / scales[None, None, :]
    margins = (peak_cap * peak_cap - np.sum(values * values, axis=1)) / (peak_cap * peak_cap)
    jacobian = np.empty((len(indices), CONTROL_COUNT), dtype=float)
    for row, (value, block) in enumerate(zip(values, normalized_design)):
        jacobian[row] = -2.0 * (value @ block) / (peak_cap * peak_cap)
    return margins, jacobian


def _load_momentum_cache(points, weights):
    if not MOMENTUM_CACHE.exists():
        return None
    try:
        with np.load(MOMENTUM_CACHE, allow_pickle=False) as loaded:
            required = {
                "points", "weights", "base_residual", "base_velocity", "base_gradient",
                "design", "design_convention", "projection_sha256", "frozen_sha256",
            }
            if not required.issubset(set(loaded.files)):
                return None
            cached_points = np.asarray(loaded["points"], dtype=float)
            cached_weights = np.asarray(loaded["weights"], dtype=float)
            if cached_points.shape != points.shape or cached_weights.shape != weights.shape:
                return None
            if not np.array_equal(cached_points, points) or not np.array_equal(cached_weights, weights):
                return None
            if str(np.asarray(loaded["design_convention"]).item()) != "real-minus-imag":
                return None
            if str(np.asarray(loaded["projection_sha256"]).item()) != _sha(PROJECTION):
                return None
            if str(np.asarray(loaded["frozen_sha256"]).item()) != _sha(FROZEN):
                return None
            arrays = {key: np.asarray(loaded[key]) for key in required - {"points", "weights"}}
            if arrays["base_residual"].shape != (len(points), 3):
                return None
            if arrays["base_velocity"].shape != (len(points), 3):
                return None
            if arrays["base_gradient"].shape != (len(points), 3, 3):
                return None
            if arrays["design"].shape != (3 * len(points), CONTROL_COUNT):
                return None
            return arrays
    except (OSError, ValueError, KeyError):
        return None


def _save_momentum_cache(points, weights, base_data, design):
    # Rebuildable local cache; the report records its path and convention.
    np.savez(
        MOMENTUM_CACHE,
        points=np.asarray(points, dtype=float),
        weights=np.asarray(weights, dtype=float),
        base_residual=np.asarray(base_data["residual"], dtype=float),
        base_velocity=np.asarray(base_data["velocity"], dtype=float),
        base_gradient=np.asarray(base_data["gradient"], dtype=float),
        design=np.asarray(design, dtype=float),
        design_convention=np.asarray("real-minus-imag"),
        projection_sha256=np.asarray(_sha(PROJECTION)),
        frozen_sha256=np.asarray(_sha(FROZEN)),
    )


def _replay_linear(control, base_residual, design, weights):
    residual = base_residual + (design @ np.asarray(control, dtype=float)).reshape(-1, 3)
    return residual, _metric(residual, weights)


def _corrected_acceleration_design(points, center, widths, carrier, dt,
                                   base_velocity, base_gradient, nu):
    """Build the endpoint momentum response with correct Re/Im columns.

    For a complex mode, the physical field is ``Re(V c)``.  Consequently the
    real coefficient column is ``Re(response)`` and the imaginary coefficient
    column is ``-Im(response)``.  The historical projection helper multiplied
    by ``-1`` and then always selected ``real``, producing ``-Re``; this local
    builder deliberately avoids that frozen bug.
    """

    points = np.asarray(points, dtype=float)
    base_velocity = np.asarray(base_velocity, dtype=float)
    base_gradient = np.asarray(base_gradient, dtype=float)
    modes = (0, 1, 2, 3, 4)
    degree = 2
    q = (degree + 1) ** 2
    velocity_count = 3 * q
    columns = []
    layout = []
    for mode in modes:
        values, gradients, diffusion, _, pressure_gradient = basis_jets(
            points, center, widths, mode, degree, mode * np.asarray(carrier), nu
        )
        components = (("real", 1.0),) if mode == 0 else (("real", 1.0), ("imag", -1.0))
        for kind, tensor, count in (
            ("velocity", values, velocity_count),
            ("pressure", pressure_gradient, q),
        ):
            for index in range(count):
                if kind == "velocity":
                    response = (
                        dt * tensor[:, :, index]
                        + 0.5 * dt**2 * diffusion[:, :, index]
                        + 0.5 * dt**2 * (
                            np.einsum("nij,nj->ni", gradients[:, :, :, index], base_velocity)
                            + np.einsum("nij,nj->ni", base_gradient, tensor[:, :, index])
                        )
                    )
                else:
                    response = dt * tensor[:, :, index]
                for component, sign in components:
                    columns.append((response.real if sign > 0.0 else -response.imag).reshape(-1))
                    layout.append({
                        "mode": int(mode),
                        "kind": kind,
                        "component": component,
                        "index": int(index),
                        "coefficient_units": "velocity acceleration" if kind == "velocity" else "pressure slope",
                    })
    design = np.stack(columns, axis=1)
    if design.shape != (3 * len(points), CONTROL_COUNT):
        raise ValueError(f"Corrected acceleration design has shape {design.shape}")
    return design, layout


def run(output_path=OUTPUT):
    started = time.perf_counter()
    output_path = Path(output_path)
    projection_raw = PROJECTION.read_bytes()
    frozen_raw = FROZEN.read_bytes()
    projection = json.loads(projection_raw)
    frozen = json.loads(frozen_raw)
    if projection.get("status") != "completed":
        raise ValueError("endpoint_acceleration_projection.json is not completed")

    with np.load(REFINED_CACHE, allow_pickle=False) as loaded:
        points = np.asarray(loaded["points"], dtype=float)
        weights = np.asarray(loaded["weights"], dtype=float)
    if points.shape != (44400, 3) or weights.shape != (44400,):
        raise ValueError(f"Unexpected refined cache grid {points.shape} {weights.shape}")
    hspace = float(projection["inputs"]["hspace"])
    htime = float(projection["inputs"]["htime"])
    nu = float(projection["inputs"]["viscosity"])
    report = {
        "status": "initializing",
        "accepted": False,
        "pde_validated": False,
        "constraints_maintained": False,
        "scope": (
            "Frozen-grid constrained endpoint acceleration diagnostic. The 324-real-column "
            "momentum tangent is fit with analytic radial-contraction, aspect-increase, "
            "weighted-spin, and sampled peak-norm gates. The acceleration vanishes at t0. "
            "No nonlinear endpoint replay, trajectory, PDE, or recursion acceptance."
        ),
        "sources": {
            "projection": {"path": PROJECTION.name, "sha256": _sha(PROJECTION)},
            "frozen_geometry": {"path": FROZEN.name, "sha256": _sha(FROZEN)},
            "refined_grid": {"path": REFINED_CACHE.name, "sha256": _sha(REFINED_CACHE)},
            "shape_oracle": {"path": "acceleration_shape_oracle.py", "sha256": _sha(ROOT / "acceleration_shape_oracle.py")},
        },
        "inputs": {
            "point_count": int(len(points)),
            "control_count": CONTROL_COUNT,
            "delta_k": float(projection["inputs"]["delta_k"]),
            "tau0": float(projection["inputs"]["tau0"]),
            "tau1": float(projection["inputs"]["tau1"]),
            "physical_time_increment": float(projection["inputs"]["physical_time_increment"]),
            "hspace": hspace,
            "htime": htime,
            "viscosity": nu,
            "ridge": RIDGE,
            "peak_relative_cap": PEAK_RELATIVE_CAP,
            "initial_pool_size": INITIAL_POOL,
            "max_rounds": MAX_ROUNDS,
            "trust_radius": TRUST_RADIUS,
        },
    }
    if PRESSURE_SEED.exists():
        report["sources"]["pressure_seed"] = {"path": PRESSURE_SEED.name, "sha256": _sha(PRESSURE_SEED)}
    output_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"stage": "sources_loaded", "point_count": len(points)}), flush=True)

    # Rebuild the actual endpoint parent once.  Subsequent constrained fits
    # reuse the local ignored cache, including the corrected Re/-Im design.
    cached = _load_momentum_cache(points, weights)
    cache_reused = cached is not None
    if cached is None:
        base_data = _reconstruct_base(projection, frozen, points, hspace, htime)
    else:
        geometry = frozen["inputs"]["wave"]
        base_data = {
            "base": None,
            "mean": None,
            "mean_report": {"status": "reused_momentum_cache"},
            "residual": cached["base_residual"],
            "velocity": cached["base_velocity"],
            "gradient": cached["base_gradient"],
            "tau0": float(projection["inputs"]["tau0"]),
            "tau1": float(projection["inputs"]["tau1"]),
            "base_control": np.asarray(projection["inputs"]["base_tangent_coefficients"], dtype=float),
            "wave": _unpack_wave(projection["inputs"]["initial_wave_coefficients"]),
            "center": tuple(float(value) for value in geometry["center"]),
            "widths": tuple(float(value) for value in geometry["widths"]),
            "carrier": np.asarray(geometry["carrier"], dtype=float),
        }
    base_residual = base_data["residual"]
    base_metric = _metric(base_residual, weights)
    report["base_endpoint"] = {
        "method": "actual MixedAffineField plus affine_momentum.jets five-point FD",
        "cache_reused": cache_reused,
        "cache_path": MOMENTUM_CACHE.name if cache_reused else None,
        "metric": base_metric,
        "reference_projection_metric": projection["base_endpoint"]["metric"],
        "metric_L2_difference_from_projection": float(
            base_metric["momentum_volume_L2"]
            - projection["base_endpoint"]["metric"]["momentum_volume_L2"]
        ),
        "max_residual": float(base_metric["momentum_max"]),
    }
    report["status"] = "base_endpoint_replayed"
    output_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"stage": "base_endpoint_replayed", "l2": base_metric["momentum_volume_L2"], "peak": base_metric["momentum_max"]}), flush=True)

    geometry = frozen["inputs"]["wave"]
    if cached is None:
        design, layout = _corrected_acceleration_design(
            points,
            base_data["center"],
            base_data["widths"],
            base_data["carrier"],
            base_data["tau0"] - base_data["tau1"],
            base_data["velocity"],
            base_data["gradient"],
            nu,
        )
        _save_momentum_cache(points, weights, base_data, design)
    else:
        design = np.asarray(cached["design"], dtype=float)
        layout = []
    if design.shape != (133200, CONTROL_COUNT):
        raise ValueError(f"Unexpected acceleration design shape {design.shape}")
    row_weight = np.repeat(np.sqrt(weights), 3)
    weighted_columns = design * row_weight[:, None]
    scales = np.maximum(np.linalg.norm(weighted_columns, axis=0), 1.0e-300)
    # Release the extra normalized copy before the optimizer; objective and
    # gradients below use the physical design with normalized coordinates.
    del weighted_columns
    report["acceleration_design"] = {
        "shape": list(design.shape),
        "builder": "shape_constrained_acceleration._corrected_acceleration_design",
        "imaginary_column_convention": "real column Re(response), imaginary column -Im(response)",
        "cache": {
            "path": MOMENTUM_CACHE.name,
            "reused": bool(cache_reused),
            "regenerable": True,
            "keys": ["points", "weights", "base_residual", "base_velocity", "base_gradient", "design"],
        },
        "layout": layout,
        "column_scales_min": float(np.min(scales)),
        "column_scales_max": float(np.max(scales)),
        "rank_rcond_1e-10": int(np.linalg.matrix_rank(design * row_weight[:, None] / scales[None, :], tol=1.0e-10)),
        "response_definition": (
            "dt*V + 0.5*dt^2*(-nu*lapV + (gradV)u_endpoint + "
            "(gradU_endpoint)V) for velocity acceleration; dt*gradP for pressure slope"
        ),
    }
    report["status"] = "acceleration_design_assembled"
    output_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"stage": "acceleration_design_assembled", "shape": list(design.shape)}), flush=True)

    oracle = AccelerationShapeOracle()
    reference = np.asarray(oracle.reference, dtype=float)
    physical_requested = np.asarray((1.0e-6, 1.0e-6, 1.0e-3), dtype=float)
    # The corrected oracle now agrees with the independent spatial FD check;
    # enforce the requested aspect margin directly.
    requested = np.asarray((1.0e-6, 1.0e-6, 1.0e-3), dtype=float)
    signs = np.asarray((-1.0, 1.0, np.sign(reference[2])), dtype=float)
    reference_scale = np.maximum(np.abs(reference), 1.0e-30)
    zero_control = np.zeros(CONTROL_COUNT, dtype=float)
    ridge_control = np.asarray(projection["ridge_fit"]["coefficients"], dtype=float)
    if ridge_control.shape != (CONTROL_COUNT,):
        raise ValueError("Frozen acceleration projection has unexpected coefficient count")
    zero_values, zero_jacobian = oracle.evaluate(zero_control)
    ridge_values, _ = oracle.evaluate(ridge_control)
    zero_shape_margin = _signed_shape(zero_values, reference, requested)
    ridge_shape_margin = _signed_shape(ridge_values, reference, requested)
    peak_cap = base_metric["momentum_max"] * (1.0 + PEAK_RELATIVE_CAP)
    zero_residual, zero_metric = _replay_linear(zero_control, base_residual, design, weights)
    ridge_residual, ridge_metric = _replay_linear(ridge_control, base_residual, design, weights)
    report["shape_constraints"] = {
        "reference_observables": reference.tolist(),
        "physical_requested_signed_fractional_change": physical_requested.tolist(),
        "oracle_enforced_signed_fractional_change": requested.tolist(),
        "sign_convention": "radial contraction = negative radial change; aspect and signed weighted spin increase",
        "zero_acceleration_observables": zero_values.tolist(),
        "zero_acceleration_margin": zero_shape_margin.tolist(),
        "ridge_observables": ridge_values.tolist(),
        "ridge_margin": ridge_shape_margin.tolist(),
        "zero_reference_state_unchanged": True,
        "oracle_point_count": int(oracle.point_count),
    }
    report["seed_metrics"] = {
        "zero_acceleration": zero_metric,
        "ridge_acceleration": ridge_metric,
        "ridge_l2_improvement_over_zero": float(zero_metric["momentum_volume_L2"] - ridge_metric["momentum_volume_L2"]),
        "ridge_peak_change_from_zero": float(ridge_metric["momentum_max"] - zero_metric["momentum_max"]),
        "peak_cap": float(peak_cap),
    }
    output_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"stage": "seeds_evaluated", "zero_l2": zero_metric["momentum_volume_L2"], "ridge_l2": ridge_metric["momentum_volume_L2"], "ridge_shape_margin": ridge_shape_margin.tolist()}), flush=True)

    coordinate_scale = max(zero_metric["momentum_volume_L2"], 1.0)
    objective_scale = coordinate_scale ** 2
    base_flat = base_residual.reshape(-1)

    def control_from_y(y):
        return coordinate_scale * np.asarray(y, dtype=float) / scales

    def objective(y):
        control = control_from_y(y)
        residual = base_flat + design @ control
        weighted = residual * row_weight
        return 0.5 * float(weighted @ weighted) / objective_scale + 0.5 * RIDGE * float(np.dot(y, y))

    def objective_gradient(y):
        control = control_from_y(y)
        residual = base_flat + design @ control
        weighted = residual * row_weight
        # The objective uses ||row_weight * residual||^2.  Apply the row
        # weights a second time when differentiating through that norm.
        return (coordinate_scale * (design.T @ (weighted * row_weight)) / scales / objective_scale
                + RIDGE * np.asarray(y, dtype=float))

    def shape_constraint(y):
        control = control_from_y(y)
        values, jacobian = oracle.evaluate(control)
        margin = _signed_shape(values, reference, requested)
        normalized_margin = margin / requested
        normalized_jacobian = (
            signs[:, None] * jacobian * coordinate_scale / scales[None, :]
            / reference_scale[:, None]
            / requested[:, None]
        )
        return normalized_margin, normalized_jacobian

    seed_magnitudes = np.linalg.norm(base_residual, axis=1)
    pool = set(np.argsort(seed_magnitudes)[-INITIAL_POOL:].tolist())

    def pool_constraint(y, indices):
        margins, jacobian = _pool_data(y, indices, base_residual, design, scales, coordinate_scale, peak_cap)
        return margins, jacobian

    def global_replay(control):
        residual = base_residual + (design @ np.asarray(control, dtype=float)).reshape(-1, 3)
        values, _ = oracle.evaluate(np.asarray(control, dtype=float))
        shape_margin = _signed_shape(values, reference, requested)
        return {
            "control": np.asarray(control, dtype=float),
            "residual": residual,
            "metrics": _metric(residual, weights),
            "shape_observables": values,
            "shape_margin": shape_margin,
            "peak_cap": peak_cap,
            "peak_margin": float(peak_cap - np.max(np.linalg.norm(residual, axis=1))),
            "feasible": bool(
                np.min(shape_margin) >= -1.0e-9
                and np.max(np.linalg.norm(residual, axis=1)) <= peak_cap * (1.0 + 1.0e-10)
            ),
        }

    pressure_seed = None
    pressure_seed_report = None
    if PRESSURE_SEED.exists():
        try:
            pressure_seed_report = json.loads(PRESSURE_SEED.read_text(encoding="utf-8"))
            pressure_control = np.asarray(pressure_seed_report["coefficients"], dtype=float)
            if (
                pressure_seed_report.get("status") == "completed"
                and pressure_control.shape == (CONTROL_COUNT,)
                and pressure_seed_report.get("design_convention") == "real-minus-imag"
                and pressure_seed_report.get("cache_sha256") == _sha(MOMENTUM_CACHE)
            ):
                pressure_seed = global_replay(pressure_control)
        except (OSError, ValueError, KeyError, json.JSONDecodeError):
            pressure_seed = None

    report["optimization"] = {
        "status": "running",
        "normalized_coordinate": "y = z / base_endpoint_L2; physical_control = y * base_endpoint_L2 / column_scale",
        "objective": "weighted endpoint momentum L2 plus 1e-6 normalized ridge",
        "peak_constraint": "normalized squared norm cap on the active point pool, expanded by global linear replay",
        "shape_constraint_tolerance": 1.0e-9,
        "trust_region": "componentwise y bounds around the pressure-only feasible seed",
        "rounds": [],
    }
    gradient_check_center = (
        pressure_seed["control"] * scales / coordinate_scale
        if pressure_seed is not None and pressure_seed["feasible"]
        else np.zeros(CONTROL_COUNT, dtype=float)
    )
    gradient_check_direction = np.linspace(-1.0, 1.0, CONTROL_COUNT, dtype=float)
    gradient_check_direction /= np.linalg.norm(gradient_check_direction)
    gradient_check_epsilon = 1.0e-6
    gradient_check_fd = (
        objective(gradient_check_center + gradient_check_epsilon * gradient_check_direction)
        - objective(gradient_check_center - gradient_check_epsilon * gradient_check_direction)
    ) / (2.0 * gradient_check_epsilon)
    gradient_check_analytic = float(
        objective_gradient(gradient_check_center) @ gradient_check_direction
    )
    report["optimization"]["objective_gradient_check"] = {
        "center": "pressure_only_feasible_seed" if pressure_seed is not None and pressure_seed["feasible"] else "zero",
        "epsilon": gradient_check_epsilon,
        "finite_difference_directional": float(gradient_check_fd),
        "analytic_directional": gradient_check_analytic,
        "absolute_error": float(abs(gradient_check_analytic - gradient_check_fd)),
        "relative_error": float(
            abs(gradient_check_analytic - gradient_check_fd)
            / max(abs(gradient_check_fd), 1.0e-12)
        ),
    }
    if pressure_seed is not None:
        report["optimization"]["pressure_seed"] = {
            "path": PRESSURE_SEED.name,
            "metric": pressure_seed["metrics"],
            "shape_margin": pressure_seed["shape_margin"].tolist(),
            "peak_margin": pressure_seed["peak_margin"],
            "feasible": pressure_seed["feasible"],
        }
    output_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"stage": "optimization_started", "pool_size": len(pool)}), flush=True)

    best = global_replay(zero_control)
    best_y = np.zeros(CONTROL_COUNT, dtype=float)
    best_reason = "zero_acceleration_seed"
    if pressure_seed is not None and pressure_seed["feasible"]:
        pressure_y = pressure_seed["control"] * scales / coordinate_scale
        if pressure_seed["metrics"]["momentum_volume_L2"] < best["metrics"]["momentum_volume_L2"]:
            best = pressure_seed
            best_y = pressure_y
            best_reason = "pressure_only_feasible_seed"
    best_peak = best
    best_peak_y = best_y.copy()
    for round_index in range(MAX_ROUNDS):
        indices = sorted(pool)
        constraints = [
            {"type": "ineq", "fun": lambda y, idx=indices: shape_constraint(y)[0],
             "jac": lambda y, idx=indices: shape_constraint(y)[1]},
            {"type": "ineq", "fun": lambda y, idx=indices: pool_constraint(y, idx)[0],
             "jac": lambda y, idx=indices: pool_constraint(y, idx)[1]},
        ]
        fit = minimize(
            objective,
            best_y,
            jac=objective_gradient,
            bounds=[(float(value - TRUST_RADIUS), float(value + TRUST_RADIUS)) for value in best_y],
            method="SLSQP",
            constraints=constraints,
            options={"maxiter": 100, "ftol": 1.0e-10, "disp": False},
        )
        candidate_y = np.asarray(fit.x if fit.x is not None else best_y, dtype=float)
        candidate_control = control_from_y(candidate_y)
        candidate = global_replay(candidate_control)
        magnitudes = np.linalg.norm(candidate["residual"], axis=1)
        violations = np.flatnonzero(magnitudes > peak_cap * (1.0 + 1.0e-10))
        if candidate["feasible"] and candidate["metrics"]["momentum_volume_L2"] < best["metrics"]["momentum_volume_L2"]:
            best = candidate
            best_y = candidate_y.copy()
            best_reason = "feasible_lower_L2"
        if candidate["feasible"] and candidate["metrics"]["momentum_max"] < best_peak["metrics"]["momentum_max"]:
            best_peak = candidate
            best_peak_y = candidate_y.copy()
        worst = np.argsort(magnitudes)[-10:][::-1]
        added = [int(index) for index in worst if int(index) not in pool]
        for index in added:
            pool.add(index)
        report["optimization"]["rounds"].append({
            "round": round_index + 1,
            "pool_size": len(pool),
            "optimizer_success": bool(fit.success),
            "optimizer_status": int(fit.status),
            "optimizer_message": str(fit.message),
            "optimizer_iterations": int(getattr(fit, "nit", -1)),
            "candidate_feasible": bool(candidate["feasible"]),
            "candidate_metric": candidate["metrics"],
            "candidate_shape_margin": candidate["shape_margin"].tolist(),
            "candidate_peak_margin": candidate["peak_margin"],
            "global_peak": float(np.max(magnitudes)),
            "newly_added_worst_points": added,
        })
        report["status"] = "optimization_round"
        output_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
        print(json.dumps({"stage": "optimization_round", "round": round_index + 1, "pool_size": len(pool), "success": bool(fit.success), "candidate_feasible": bool(candidate["feasible"]), "candidate_l2": candidate["metrics"]["momentum_volume_L2"], "candidate_shape_margin": candidate["shape_margin"].tolist(), "global_peak": float(np.max(magnitudes))}), flush=True)
        if candidate["feasible"] and not len(violations):
            break
        if not added and not len(violations):
            break

    zero_final = global_replay(zero_control)
    selected = best if best["feasible"] and best["metrics"]["momentum_volume_L2"] < zero_final["metrics"]["momentum_volume_L2"] else zero_final
    if selected is zero_final:
        best_reason = "zero_acceleration_seed_retained_no_feasible_improvement"
    selected_control = selected["control"]
    report["selected"] = {
        "reason": best_reason,
        "acceleration_coefficients": selected_control.tolist(),
        "momentum": selected["metrics"],
        "momentum_l2_improvement_over_zero": float(zero_final["metrics"]["momentum_volume_L2"] - selected["metrics"]["momentum_volume_L2"]),
        "momentum_max_change_from_zero": float(selected["metrics"]["momentum_max"] - zero_final["metrics"]["momentum_max"]),
        "momentum_peak_cap": float(peak_cap),
        "peak_margin": float(selected["peak_margin"]),
        "shape_observables": selected["shape_observables"].tolist(),
        "shape_margin": selected["shape_margin"].tolist(),
        "assembled_feasible": bool(selected["feasible"]),
        "reference_state_unchanged": True,
    }
    report["best_peak_feasible"] = {
        "momentum": best_peak["metrics"],
        "momentum_max_change_from_zero": float(best_peak["metrics"]["momentum_max"] - zero_final["metrics"]["momentum_max"]),
        "shape_margin": best_peak["shape_margin"].tolist(),
        "peak_margin": float(best_peak["peak_margin"]),
    }
    report["optimization"]["status"] = "completed"
    report["optimization"]["final_pool_size"] = len(pool)
    report["optimization"]["coordinate_scale_base_l2"] = float(coordinate_scale)
    report["optimization"]["selected_y_norm"] = float(np.linalg.norm(selected_control * scales / coordinate_scale))
    report["status"] = "completed"
    report["elapsed_seconds"] = float(time.perf_counter() - started)
    report["constraints_maintained"] = bool(selected["feasible"])
    report["constraints_scope"] = (
        "Frozen 44400-point linearized momentum replay, analytic endpoint shape oracle, "
        "and global linear peak check; no nonlinear endpoint replay or independent validation."
    )
    output_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "stage": "completed",
        "reason": best_reason,
        "selected_l2": selected["metrics"]["momentum_volume_L2"],
        "zero_l2": zero_final["metrics"]["momentum_volume_L2"],
        "selected_peak": selected["metrics"]["momentum_max"],
        "selected_shape_margin": selected["shape_margin"].tolist(),
        "elapsed_seconds": report["elapsed_seconds"],
    }), flush=True)
    return report


if __name__ == "__main__":
    run()
