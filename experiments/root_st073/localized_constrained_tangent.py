"""Bounded constrained fit for the localized degree-2 residual patch.

The enriched 264-control candidate remains fixed.  Only the 180 controls of
``localized_residual_enrichment`` are optimized against the frozen 44,400
point residual.  The fit enforces the four moment rows, 81 cone rows, the
cached endpoint shape inequalities, and a cutting-plane approximation to the
global refined-grid peak cap.  All responses are analytic or reconstructed
from frozen arrays; no mean finite-difference replay is performed.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.optimize import LinearConstraint, minimize

ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from enriched_endpoint_shape_cache import EnrichedEndpointShape  # noqa: E402
from localized_constraint_rows import _local_endpoint_response  # noqa: E402
from localized_residual_enrichment import _patch_design  # noqa: E402
from wave_endpoint_shape_cache import _observables_and_jacobian  # noqa: E402


CACHE_PATH = ROOT / "refined_wave_momentum_cache.npz"
CACHE_REPORT_PATH = ROOT / "refined_wave_momentum_cache.json"
LOCAL_RESIDUAL_REPORT_PATH = ROOT / "localized_residual_enrichment.json"
LOCAL_CONSTRAINT_REPORT_PATH = ROOT / "localized_constraint_rows.json"
MEAN_ROWS_PATH = ROOT / "enriched_mean_constraint_rows.json"
MOMENT_PATH = ROOT / "wave_dynamics_mean_compatibility.json"
CONE_PATH = ROOT / "wave_mean_cone_projection.json"
ENDPOINT_CACHE_PATH = ROOT / "enriched_endpoint_shape_cache.npz"
ENDPOINT_CACHE_REPORT_PATH = ROOT / "enriched_endpoint_shape_cache.json"
CANDIDATE_PATH = ROOT / "enriched_mean_endpoint_tangent.json"
GEOMETRY_PATH = ROOT / "full_wave_frozen_cache.json"
OUTPUT_PATH = ROOT / "localized_constrained_tangent.json"

MOMENT_TOLERANCE = 1.0e-5
CONE_TOLERANCE = 1.0e-7
CONE_OPTIMIZER_SAFETY = 1.0e-9
ENDPOINT_TOLERANCE = 1.0e-10
PEAK_RELATIVE_CAP = 1.0e-6
PEAK_FEASIBILITY_SLACK = 1.0e-10
MAX_ROUNDS = 5
POOL_INITIAL = 20
POOL_ADD = 10


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _metric(residual, weights):
    residual = np.asarray(residual, dtype=float)
    weights = np.asarray(weights, dtype=float)
    magnitudes = np.linalg.norm(residual, axis=1)
    square = float(np.sum(weights * np.sum(residual * residual, axis=1)))
    return {
        "point_count": int(len(residual)),
        "momentum_max": float(np.max(magnitudes)),
        "momentum_volume_L2": float(np.sqrt(max(square, 0.0))),
        "momentum_volume_RMS": float(np.sqrt(max(square, 0.0) / np.sum(weights))),
        "physical_volume": float(np.sum(weights)),
    }


def _max_detail(residual, points, weights):
    residual = np.asarray(residual, dtype=float)
    magnitudes = np.linalg.norm(residual, axis=1)
    index = int(np.argmax(magnitudes))
    point = np.asarray(points[index], dtype=float)
    return {
        "index": index,
        "point": point.tolist(),
        "radius": float(np.hypot(point[0], point[1])),
        "z": float(point[2]),
        "weight": float(weights[index]),
        "vector": residual[index].tolist(),
        "norm": float(magnitudes[index]),
    }


def _source_entry(path: Path) -> dict:
    return {"path": path.name, "sha256": _sha256(path)}


def run(output_path=OUTPUT_PATH, max_rounds=MAX_ROUNDS):
    started = time.perf_counter()
    output_path = Path(output_path)
    cache_report = json.loads(CACHE_REPORT_PATH.read_text(encoding="utf-8"))
    local_residual_report = json.loads(LOCAL_RESIDUAL_REPORT_PATH.read_text(encoding="utf-8"))
    local_constraint_report = json.loads(LOCAL_CONSTRAINT_REPORT_PATH.read_text(encoding="utf-8"))
    mean_rows_report = json.loads(MEAN_ROWS_PATH.read_text(encoding="utf-8"))
    moment_report = json.loads(MOMENT_PATH.read_text(encoding="utf-8"))
    cone_report = json.loads(CONE_PATH.read_text(encoding="utf-8"))
    candidate = json.loads(CANDIDATE_PATH.read_text(encoding="utf-8"))
    geometry = json.loads(GEOMETRY_PATH.read_text(encoding="utf-8"))["inputs"]["wave"]
    if cache_report.get("status") != "completed":
        raise ValueError("Refined momentum cache is not completed")
    if local_residual_report.get("status") != "completed":
        raise ValueError("Localized residual report is not completed")
    if local_constraint_report.get("status") != "completed":
        raise ValueError("Localized constraint rows are not completed")
    if candidate.get("status") != "completed" or candidate.get("control_count") != 264:
        raise ValueError("Expected completed enriched 264-control candidate")

    with np.load(CACHE_PATH, allow_pickle=False) as loaded:
        points = np.asarray(loaded["points"], dtype=float)
        weights = np.asarray(loaded["weights"], dtype=float)
        base_residual = np.asarray(loaded["residual"], dtype=float)
        old_design = np.asarray(loaded["tangent_design"], dtype=float)
        old_control_cache = np.asarray(loaded["tangent_coefficients"], dtype=float)
        cache_wave = np.asarray(loaded["coefficients_original"], dtype=float)
        center = np.asarray(loaded["center"], dtype=float)
        widths = np.asarray(loaded["widths"], dtype=float)
        carrier = np.asarray(loaded["carrier"], dtype=float)
    old_control = np.asarray(candidate["selected"]["tangent_coefficients"], dtype=float)
    wave_packed = np.asarray(candidate["selected"]["coefficients_original"], dtype=float)
    if points.shape != (44400, 3) or weights.shape != (44400,):
        raise ValueError(f"Unexpected refined grid shape {points.shape} {weights.shape}")
    if old_design.shape != (133200, 264) or base_residual.shape != (44400, 3):
        raise ValueError(f"Unexpected refined cache array shapes {old_design.shape} {base_residual.shape}")
    if old_control.shape != (264,) or not np.array_equal(old_control, old_control_cache):
        raise ValueError("Refined cache and enriched candidate use different fixed controls")
    if wave_packed.shape != cache_wave.shape or not np.array_equal(wave_packed, cache_wave):
        raise ValueError("Refined cache and enriched candidate use different initial waves")

    current_residual = base_residual + (old_design @ old_control).reshape(-1, 3)
    del old_design, base_residual
    before_metric = _metric(current_residual, weights)
    before_peak = _max_detail(current_residual, points, weights)
    peak_cap = before_metric["momentum_max"] * (1.0 + PEAK_RELATIVE_CAP)

    # Local rows and frozen baseline forms.
    local_moment_rows = np.asarray(local_constraint_report["moment_rows"], dtype=float)
    local_cone_rows = np.asarray(local_constraint_report["cone_rows"], dtype=float)
    old_moment_rows = np.asarray(mean_rows_report["moment_rows"], dtype=float)
    old_cone_rows = np.asarray(mean_rows_report["cone_control_rows"], dtype=float)
    if local_moment_rows.shape != (4, 180) or local_cone_rows.shape != (81, 180):
        raise ValueError("Localized row dimensions do not match the 180-control layout")
    if old_moment_rows.shape != (4, 64) or old_cone_rows.shape != (81, 64):
        raise ValueError("Enriched degree-3 row dimensions do not match fixed candidate")

    moment_linearization = moment_report["reusable_moment_linearization"]
    baseline_moments = np.asarray(moment_linearization["baseline_moments"], dtype=float)
    wave_moment_forms = np.asarray(moment_linearization["wave_moment_forms"], dtype=float)
    wave_x = np.r_[wave_packed[:, 0], wave_packed[:, 1]]
    moment_target = -baseline_moments - np.einsum("i,kij,j->k", wave_x, wave_moment_forms, wave_x)
    old_moment_error = old_moment_rows @ old_control[:64] - moment_target

    cone_baseline = np.asarray(cone_report["cone_baseline"], dtype=float)
    cone_wave_forms = np.asarray(cone_report["cone_wave_forms"], dtype=float)
    cone_lower = np.asarray(cone_report["cone_lower"], dtype=float)
    cone_wave_baseline = cone_baseline + np.einsum("i,kij,j->k", wave_x, cone_wave_forms, wave_x)
    old_cone_margin = cone_wave_baseline + old_cone_rows @ old_control[:64] - cone_lower

    # Build the analytic local momentum design and normalize both columns and
    # objective variables.  y = column_scale * local_control.
    local_design, local_layout = _patch_design(
        points,
        np.asarray(local_constraint_report["support_containment"]["patch_center"], dtype=float),
        np.asarray(local_constraint_report["support_containment"]["patch_widths"], dtype=float),
        carrier,
    )
    sqrt_weights = np.sqrt(weights)
    row_weight = np.repeat(sqrt_weights, 3)
    weighted_design = local_design * row_weight[:, None]
    del local_design
    column_scales = np.maximum(np.linalg.norm(weighted_design, axis=0), 1.0e-30)
    normalized_design = weighted_design / column_scales[None, :]
    del weighted_design
    weighted_base = current_residual.reshape(-1) * row_weight
    variable_scale = max(float(np.linalg.norm(weighted_base)), 1.0)
    objective_scale = variable_scale**2

    # Parent endpoint state and exact local analytic endpoint response.
    endpoint = EnrichedEndpointShape(ENDPOINT_CACHE_PATH)
    if endpoint.point_count != 7776:
        raise ValueError("Expected 7776 endpoint points")
    parent_velocity = endpoint.velocity_offset + np.einsum(
        "ncq,q->nc", endpoint.velocity_control, old_control
    )
    parent_gradient = endpoint.gradient_offset + np.einsum(
        "ncdq,q->ncd", endpoint.gradient_control, old_control
    )
    local_velocity, local_gradient, _, _ = _local_endpoint_response(
        endpoint.points,
        np.asarray(local_constraint_report["support_containment"]["patch_center"], dtype=float),
        np.asarray(local_constraint_report["support_containment"]["patch_widths"], dtype=float),
        np.asarray(geometry["carrier"], dtype=float),
        endpoint.dtau,
    )
    reference = np.asarray(endpoint.reference, dtype=float)
    signs = np.asarray((-1.0, 1.0, np.sign(reference[2])), dtype=float)
    reference_scale = np.maximum(np.abs(reference), 1.0e-30)
    endpoint_requested = np.asarray((1.0e-6, 1.0e-6, 1.0e-3), dtype=float)

    # Constraint rows act on y through local_control = y / column_scales.
    moment_matrix = variable_scale * local_moment_rows / column_scales[None, :]
    moment_row_norm = np.linalg.norm(moment_matrix, axis=1)
    moment_responsive = moment_row_norm > 1.0e-30
    if np.any(~moment_responsive & (np.abs(old_moment_error) > MOMENT_TOLERANCE)):
        raise ValueError("A fixed moment row is outside tolerance with zero local response")
    moment_matrix_n = moment_matrix[moment_responsive] / moment_row_norm[moment_responsive, None]
    moment_lower = (-MOMENT_TOLERANCE - old_moment_error[moment_responsive]) / moment_row_norm[moment_responsive]
    moment_upper = (MOMENT_TOLERANCE - old_moment_error[moment_responsive]) / moment_row_norm[moment_responsive]

    cone_matrix = variable_scale * local_cone_rows / column_scales[None, :]
    cone_row_norm = np.linalg.norm(cone_matrix, axis=1)
    cone_responsive = cone_row_norm > 1.0e-30
    fixed_cone_margin = old_cone_margin[~cone_responsive]
    if np.any(fixed_cone_margin < -CONE_TOLERANCE):
        raise ValueError("A fixed cone row is infeasible with zero local response")
    cone_matrix_n = cone_matrix[cone_responsive] / cone_row_norm[cone_responsive, None]
    cone_lower_n = (
        CONE_OPTIMIZER_SAFETY - old_cone_margin[cone_responsive]
    ) / cone_row_norm[cone_responsive]

    endpoint_velocity_response = local_velocity * (
        variable_scale / column_scales[None, None, :]
    )
    endpoint_gradient_response = local_gradient * (
        variable_scale / column_scales[None, None, None, :]
    )

    def local_control(y):
        return variable_scale * np.asarray(y, dtype=float) / column_scales

    def residual_at(y):
        weighted_correction = normalized_design @ (variable_scale * np.asarray(y, dtype=float))
        return current_residual + (weighted_correction / row_weight).reshape(-1, 3)

    def endpoint_eval(y):
        control = local_control(y)
        velocity = parent_velocity + np.einsum("ncq,q->nc", endpoint_velocity_response, y)
        gradient = parent_gradient + np.einsum("ncdq,q->ncd", endpoint_gradient_response, y)
        values, jacobian, _ = _observables_and_jacobian(
            velocity,
            gradient,
            endpoint.weights,
            endpoint.points,
            endpoint_velocity_response,
            endpoint_gradient_response,
        )
        margins = signs * (values - reference) / reference_scale - endpoint_requested
        margin_jacobian = signs[:, None] * jacobian / reference_scale[:, None]
        return margins, margin_jacobian, values

    def replay(y):
        control = local_control(y)
        residual = residual_at(y)
        endpoint_margin, _, endpoint_values = endpoint_eval(y)
        return {
            "y": np.asarray(y, dtype=float),
            "control": control,
            "residual": residual,
            "metrics": _metric(residual, weights),
            "moment_error": old_moment_error + local_moment_rows @ control,
            "cone_margin": old_cone_margin + local_cone_rows @ control,
            "endpoint_margin": endpoint_margin,
            "endpoint_values": endpoint_values,
        }

    def feasible(state):
        return bool(
            np.max(np.abs(state["moment_error"])) <= MOMENT_TOLERANCE
            and np.min(state["cone_margin"]) >= -CONE_TOLERANCE
            and np.min(state["endpoint_margin"]) >= -ENDPOINT_TOLERANCE
            and state["metrics"]["momentum_max"] <= peak_cap * (1.0 + PEAK_FEASIBILITY_SLACK)
        )

    zero_state = replay(np.zeros(180, dtype=float))
    if not feasible(zero_state):
        raise ValueError("Zero local controls are not feasible within documented tolerances")

    def objective(y):
        weighted_residual = weighted_base + normalized_design @ (variable_scale * np.asarray(y, dtype=float))
        return 0.5 * float(weighted_residual @ weighted_residual) / objective_scale

    def objective_jac(y):
        weighted_residual = weighted_base + normalized_design @ (variable_scale * np.asarray(y, dtype=float))
        return variable_scale * (normalized_design.T @ weighted_residual) / objective_scale

    objective_direction = np.linspace(-1.0, 1.0, 180, dtype=float)
    objective_direction /= np.linalg.norm(objective_direction)
    objective_fd_step = 1.0e-4
    objective_gradient_analytic = float(objective_jac(np.zeros(180, dtype=float)) @ objective_direction)
    objective_gradient_fd = float(
        (objective(objective_fd_step * objective_direction)
         - objective(-objective_fd_step * objective_direction))
        / (2.0 * objective_fd_step)
    )
    objective_gradient_check = {
        "step": objective_fd_step,
        "analytic_directional": objective_gradient_analytic,
        "finite_difference_directional": objective_gradient_fd,
        "absolute_error": abs(objective_gradient_fd - objective_gradient_analytic),
        "relative_error": abs(objective_gradient_fd - objective_gradient_analytic)
        / max(abs(objective_gradient_fd), 1.0e-30),
    }

    def endpoint_fun(y):
        return endpoint_eval(y)[0]

    def endpoint_jac(y):
        return endpoint_eval(y)[1]

    def pool_data(y, indices):
        indices = np.asarray(sorted(indices), dtype=int)
        rows = np.concatenate([np.arange(3 * index, 3 * index + 3) for index in indices])
        point_matrix = variable_scale * normalized_design[rows] / row_weight[rows, None]
        values = current_residual[indices] + (point_matrix @ np.asarray(y, dtype=float)).reshape(len(indices), 3)
        normalized_margin = 1.0 - np.sum(values * values, axis=1) / peak_cap**2
        point_blocks = point_matrix.reshape(len(indices), 3, -1)
        jacobian = np.empty((len(indices), normalized_design.shape[1]), dtype=float)
        for row, (value, block) in enumerate(zip(values, point_blocks)):
            jacobian[row] = -2.0 * (value @ block) / peak_cap**2
        return normalized_margin, jacobian

    def full_peak(y):
        residual = residual_at(y)
        return np.linalg.norm(residual, axis=1)

    pool = set(np.argsort(np.linalg.norm(current_residual, axis=1))[-POOL_INITIAL:].tolist())
    best = zero_state
    best_reason = "zero_local_controls_feasible"
    rounds = []

    for round_index in range(int(max_rounds)):
        indices = sorted(pool)

        def pool_fun(y):
            return pool_data(y, indices)[0]

        def pool_jac(y):
            return pool_data(y, indices)[1]

        constraints = []
        if np.any(moment_responsive):
            constraints.append(LinearConstraint(moment_matrix_n, moment_lower, moment_upper))
        if np.any(cone_responsive):
            constraints.append(LinearConstraint(cone_matrix_n, cone_lower_n, np.inf))
        constraints.extend([
            {"type": "ineq", "fun": endpoint_fun, "jac": endpoint_jac},
            {"type": "ineq", "fun": pool_fun, "jac": pool_jac},
        ])
        fit = minimize(
            objective,
            best["y"],
            jac=objective_jac,
            method="SLSQP",
            constraints=constraints,
            options={"maxiter": 100, "ftol": 1.0e-10, "disp": False},
        )
        candidate_y = np.asarray(fit.x if fit.x is not None else best["y"], dtype=float)
        candidate = replay(candidate_y)
        magnitudes = full_peak(candidate_y)
        violating = np.flatnonzero(magnitudes > peak_cap * (1.0 + PEAK_FEASIBILITY_SLACK))
        worst = np.argsort(magnitudes)[-POOL_ADD:][::-1]
        added = [int(index) for index in worst if int(index) not in pool]
        if feasible(candidate) and candidate["metrics"]["momentum_volume_L2"] < best["metrics"]["momentum_volume_L2"]:
            best = candidate
            best_reason = "feasible_lower_L2"
        for index in added:
            pool.add(index)
        rounds.append({
            "round": round_index + 1,
            "pool_size_before": len(indices),
            "pool_size_after": len(pool),
            "optimizer_success": bool(fit.success),
            "optimizer_status": int(fit.status),
            "optimizer_message": str(fit.message),
            "optimizer_iterations": int(getattr(fit, "nit", -1)),
            "candidate_feasible": bool(feasible(candidate)),
            "candidate_metrics": candidate["metrics"],
            "candidate_moment_max_abs": float(np.max(np.abs(candidate["moment_error"]))),
            "candidate_cone_min_margin": float(np.min(candidate["cone_margin"])),
            "candidate_endpoint_min_margin": float(np.min(candidate["endpoint_margin"])),
            "candidate_global_peak": float(np.max(magnitudes)),
            "peak_cap": float(peak_cap),
            "global_peak_violating_point_count": int(len(violating)),
            "newly_added_worst_points": added,
            "worst_point": int(np.argmax(magnitudes)),
        })
        if feasible(candidate) and len(violating) == 0:
            break
        if not added and len(violating) == 0:
            break

    selected = best
    selected_feasible = feasible(selected)
    source_paths = {
        "candidate": CANDIDATE_PATH,
        "refined_cache": CACHE_PATH,
        "refined_cache_report": CACHE_REPORT_PATH,
        "localized_residual": LOCAL_RESIDUAL_REPORT_PATH,
        "localized_constraints": LOCAL_CONSTRAINT_REPORT_PATH,
        "mean_rows": MEAN_ROWS_PATH,
        "moment": MOMENT_PATH,
        "cone": CONE_PATH,
        "endpoint_cache": ENDPOINT_CACHE_PATH,
        "endpoint_cache_report": ENDPOINT_CACHE_REPORT_PATH,
        "geometry": GEOMETRY_PATH,
    }
    report = {
        "status": "completed",
        "accepted": False,
        "pde_validated": False,
        "scale_recursion_established": False,
        "constraints_maintained": False,
        "assembled_feasible": bool(selected_feasible),
        "assembled_constraints_maintained": bool(selected_feasible),
        "source": CANDIDATE_PATH.name,
        "source_sha256": _sha256(CANDIDATE_PATH),
        "sources": {name: _source_entry(path) for name, path in source_paths.items()},
        "scope": (
            "Reference-time constrained local patch fit. The enriched 264-control "
            "candidate is fixed; only 180 localized controls vary. Four moments, "
            "81 cones, endpoint contraction/aspect/spin inequalities, and a refined "
            "grid peak cap are checked on frozen analytic/cache responses. No mean "
            "replay, finite-time, PDE, or recursive-scale acceptance."
        ),
        "inputs": {
            "point_count": int(len(points)),
            "candidate_control_count": 264,
            "local_control_count": 180,
            "patch_center": local_constraint_report["support_containment"]["patch_center"],
            "patch_widths": local_constraint_report["support_containment"]["patch_widths"],
            "patch_degree": 2,
            "patch_modes": [0, 1, 2],
            "old_controls_held_fixed": True,
            "column_normalization": "weighted L2 norm of each local momentum-response column",
            "variable_normalization": "local normalized variables divided by fixed-candidate weighted L2",
            "objective": "normalized weighted residual L2 using frozen base residual plus fixed 264 response",
            "peak_cap": float(peak_cap),
            "peak_cap_relative_to_fixed_candidate": PEAK_RELATIVE_CAP,
            "peak_cutting_plane_round_limit": int(max_rounds),
            "tolerances": {
                "moment_max_abs": MOMENT_TOLERANCE,
                "cone_min_margin": CONE_TOLERANCE,
                "cone_optimizer_safety": CONE_OPTIMIZER_SAFETY,
                "endpoint_min_margin": ENDPOINT_TOLERANCE,
                "peak_feasibility_slack": PEAK_FEASIBILITY_SLACK,
            },
            "endpoint_requested_signed_fractional_change": endpoint_requested.tolist(),
        },
        "zero_local_control_feasibility": {
            "feasible": bool(feasible(zero_state)),
            "moment_max_abs": float(np.max(np.abs(zero_state["moment_error"]))),
            "cone_min_margin": float(np.min(zero_state["cone_margin"])),
            "cone_pass_count": int(np.sum(zero_state["cone_margin"] >= -CONE_TOLERANCE)),
            "endpoint_margin": zero_state["endpoint_margin"].tolist(),
            "endpoint_min_margin": float(np.min(zero_state["endpoint_margin"])),
            "momentum": zero_state["metrics"],
            "peak_cap": float(peak_cap),
        },
        "selected": {
            "reason": best_reason,
            "feasible": bool(selected_feasible),
            "local_patch_coefficients": selected["control"].tolist(),
            "tangent_coefficients": selected["control"].tolist(),
            "local_patch_coefficient_norm": float(np.linalg.norm(selected["control"])),
            "momentum": selected["metrics"],
            "peak_detail": _max_detail(selected["residual"], points, weights),
            "momentum_l2_improvement": float(before_metric["momentum_volume_L2"] - selected["metrics"]["momentum_volume_L2"]),
            "momentum_max_change": float(selected["metrics"]["momentum_max"] - before_metric["momentum_max"]),
            "moment_max_abs": float(np.max(np.abs(selected["moment_error"]))),
            "cone_min_margin": float(np.min(selected["cone_margin"])),
            "cone_pass_count": int(np.sum(selected["cone_margin"] >= -CONE_TOLERANCE)),
            "endpoint_values": selected["endpoint_values"].tolist(),
            "endpoint_margin": selected["endpoint_margin"].tolist(),
            "endpoint_min_margin": float(np.min(selected["endpoint_margin"])),
            "old_controls_unchanged": True,
            "active_constraints": {
                "moment_rows_within_1e-8_of_tolerance": [
                    int(index)
                    for index, value in enumerate(selected["moment_error"])
                    if abs(value) >= MOMENT_TOLERANCE - 1.0e-8
                ],
                "cone_rows_within_1e-8_of_tolerance": [
                    int(index)
                    for index, value in enumerate(selected["cone_margin"])
                    if value <= -CONE_TOLERANCE + 1.0e-8
                ],
                "endpoint_rows_within_1e-8": [
                    int(index)
                    for index, value in enumerate(selected["endpoint_margin"])
                    if value <= 1.0e-8
                ],
                "peak_cap_fractional_slack": float(
                    peak_cap / selected["metrics"]["momentum_max"] - 1.0
                ),
            },
        },
        "fixed_candidate_before": {
            "momentum": before_metric,
            "peak_detail": before_peak,
            "moment_max_abs": float(np.max(np.abs(old_moment_error))),
            "cone_min_margin": float(np.min(old_cone_margin)),
            "cone_pass_count": int(np.sum(old_cone_margin >= -CONE_TOLERANCE)),
        },
        "cutting_plane_rounds": rounds,
        "matrix_diagnostics": {
            "local_design_shape": [133200, 180],
            "column_scale_min": float(np.min(column_scales)),
            "column_scale_max": float(np.max(column_scales)),
            "variable_scale_fixed_candidate_weighted_L2": float(variable_scale),
            "objective_gradient_check": objective_gradient_check,
            "moment_row_rank": int(np.linalg.matrix_rank(local_moment_rows)),
            "moment_responsive_rows": int(np.sum(moment_responsive)),
            "cone_responsive_rows": int(np.sum(cone_responsive)),
            "endpoint_response_shape": list(local_velocity.shape),
            "endpoint_gradient_response_shape": list(local_gradient.shape),
        },
        "elapsed_seconds": float(time.perf_counter() - started),
    }
    output_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": report["status"],
        "selected_feasible": selected_feasible,
        "best_reason": best_reason,
        "before_l2": before_metric["momentum_volume_L2"],
        "selected_l2": selected["metrics"]["momentum_volume_L2"],
        "before_max": before_metric["momentum_max"],
        "selected_max": selected["metrics"]["momentum_max"],
        "moment_max_abs": report["selected"]["moment_max_abs"],
        "cone_min_margin": report["selected"]["cone_min_margin"],
        "endpoint_min_margin": report["selected"]["endpoint_min_margin"],
        "rounds": len(rounds),
        "elapsed_seconds": report["elapsed_seconds"],
    }), flush=True)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT_PATH)
    parser.add_argument("--max-rounds", type=int, default=MAX_ROUNDS)
    args = parser.parse_args()
    run(args.output, args.max_rounds)
