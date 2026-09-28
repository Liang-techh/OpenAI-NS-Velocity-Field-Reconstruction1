"""Bounded two-patch fit with a target profile-drift cap of 1800.

The enriched 264-control candidate and the two-patch constrained candidate
are frozen inputs.  This diagnostic warm-starts from the latter and continues
the exact finite-interval drift cap down to 1800, while retaining
the four moments, 81 cones, endpoint shape inequalities, and first-patch peak
cap.  A target candidate is selected only when every gate and the 1800 drift
cap pass.  All responses are analytic or reconstructed from frozen arrays;
there is no mean replay, holdout fit, PDE, or recursive-scale acceptance.
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

from broad_shear_dynamic_control import load_saved_field  # noqa: E402
from enriched_endpoint_shape_cache import EnrichedEndpointShape  # noqa: E402
from grouped_joined_field import install_in_field  # noqa: E402
from localized_constraint_rows import _local_endpoint_response  # noqa: E402
from localized_residual_enrichment import _patch_design  # noqa: E402
from wave_dynamics_mean_compatibility import _mode0_moment_rows  # noqa: E402
from wave_endpoint_shape_cache import _observables_and_jacobian  # noqa: E402
from wave_mean_cone_projection import _cone_geometry, _mode0_rows  # noqa: E402


CACHE_PATH = ROOT / "refined_wave_momentum_cache.npz"
CACHE_REPORT_PATH = ROOT / "refined_wave_momentum_cache.json"
MEAN_ROWS_PATH = ROOT / "enriched_mean_constraint_rows.json"
MOMENT_PATH = ROOT / "wave_dynamics_mean_compatibility.json"
CONE_PATH = ROOT / "wave_mean_cone_projection.json"
ENDPOINT_CACHE_PATH = ROOT / "enriched_endpoint_shape_cache.npz"
ENDPOINT_CACHE_REPORT_PATH = ROOT / "enriched_endpoint_shape_cache.json"
CANDIDATE_PATH = ROOT / "enriched_mean_endpoint_tangent.json"
FIRST_FIT_PATH = ROOT / "localized_constrained_tangent.json"
FIRST_ROWS_PATH = ROOT / "localized_constraint_rows.json"
SECOND_PROBE_PATH = ROOT / "localized_next_patch_probe.json"
MEAN_PATH = ROOT / "broad_meridional_constrained.json"
DYNAMIC_PATH = ROOT / "broad_shear_dynamic_control.json"
GEOMETRY_PATH = ROOT / "full_wave_frozen_cache.json"
DRIFT_CACHE_PATH = ROOT / "localized_drift_constraint.npz"
DRIFT_REPORT_PATH = ROOT / "localized_drift_constraint.json"
TWO_PATCH_PATH = ROOT / "localized_parent_drift_fit.json"
OUTPUT_PATH = ROOT / "localized_drift_1800_fit.json"
TARGET_DRIFT_CAP = 1800.0

MOMENT_TOLERANCE = 1.0e-5
CONE_TOLERANCE = 1.0e-7
CONE_OPTIMIZER_SAFETY = 1.0e-9
ENDPOINT_TOLERANCE = 1.0e-10
PEAK_RELATIVE_CAP = 1.0e-6
PEAK_FEASIBILITY_SLACK = 1.0e-10
DRIFT_SOLVER_SAFETY = 1.0e-5
MAX_ROUNDS = 5
POOL_INITIAL = 20
POOL_ADD = 10


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _source_entry(path: Path) -> dict:
    return {"path": path.name, "sha256": _sha256(path)}


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


def _build_constraint_rows(mean_report, geometry, center, widths):
    """Build rows with the same frozen mode-0 definitions as the first patch."""

    dynamic, _ = load_saved_field(DYNAMIC_PATH)
    install_in_field(dynamic)
    tau = float(mean_report["tau"])
    k = float(mean_report["k"])
    breaks = [float(value) for value in mean_report["quadrature"]["radial_split_breaks"]]
    moment_rows, moment_geometry = _mode0_moment_rows(
        dynamic.inner, tau, breaks, np.asarray(center, dtype=float), np.asarray(widths, dtype=float)
    )
    cone_points, cone_panels, cone_H = _cone_geometry(
        mean_report, dynamic.inner, k, breaks
    )
    _, cone_rows = _mode0_rows(
        cone_points,
        cone_panels,
        cone_H,
        np.asarray(center, dtype=float),
        np.asarray(widths, dtype=float),
    )
    if moment_rows.shape != (4, 180) or cone_rows.shape != (81, 180):
        raise ValueError(f"Unexpected local row shapes {moment_rows.shape} {cone_rows.shape}")
    return np.asarray(moment_rows, dtype=float), np.asarray(cone_rows, dtype=float), moment_geometry, {
        "location_count": int(len(cone_panels)),
        "inequality_count": int(len(cone_rows)),
        "grouped_point_count": int(len(cone_points)),
        "cone_order": 64,
        "H_shape": list(np.asarray(cone_H).shape),
    }


def run(output_path=OUTPUT_PATH, max_rounds=MAX_ROUNDS):
    started = time.perf_counter()
    output_path = Path(output_path)
    prior_parent_report = None
    prior_parent_hash = None
    if output_path.exists():
        try:
            prior_parent_report = json.loads(output_path.read_text(encoding="utf-8"))
            prior_parent_hash = _sha256(output_path)
        except (OSError, json.JSONDecodeError):
            prior_parent_report = None
            prior_parent_hash = None
    cache_report = json.loads(CACHE_REPORT_PATH.read_text(encoding="utf-8"))
    mean_rows_report = json.loads(MEAN_ROWS_PATH.read_text(encoding="utf-8"))
    moment_report = json.loads(MOMENT_PATH.read_text(encoding="utf-8"))
    cone_report = json.loads(CONE_PATH.read_text(encoding="utf-8"))
    candidate = json.loads(CANDIDATE_PATH.read_text(encoding="utf-8"))
    first_fit = json.loads(FIRST_FIT_PATH.read_text(encoding="utf-8"))
    two_patch = json.loads(TWO_PATCH_PATH.read_text(encoding="utf-8"))
    first_rows = json.loads(FIRST_ROWS_PATH.read_text(encoding="utf-8"))
    second_probe = json.loads(SECOND_PROBE_PATH.read_text(encoding="utf-8"))
    mean_report = json.loads(MEAN_PATH.read_text(encoding="utf-8"))
    geometry = json.loads(GEOMETRY_PATH.read_text(encoding="utf-8"))["inputs"]["wave"]
    if cache_report.get("status") != "completed":
        raise ValueError("Refined momentum cache is not completed")
    if candidate.get("status") != "completed" or candidate.get("control_count") != 264:
        raise ValueError("Expected completed enriched 264-control candidate")
    if first_fit.get("status") != "completed":
        raise ValueError("First constrained local fit is not completed")
    if two_patch.get("status") != "completed" or not two_patch.get("selected", {}).get("feasible", False):
        raise ValueError("Frozen two-patch candidate is not completed and feasible")
    if first_rows.get("status") != "completed" or second_probe.get("status") != "completed":
        raise ValueError("Localized row/probe reports are not completed")

    with np.load(CACHE_PATH, allow_pickle=False) as loaded:
        points = np.asarray(loaded["points"], dtype=float)
        weights = np.asarray(loaded["weights"], dtype=float)
        base_residual = np.asarray(loaded["residual"], dtype=float)
        old_design = np.asarray(loaded["tangent_design"], dtype=float)
        old_control_cache = np.asarray(loaded["tangent_coefficients"], dtype=float)
        cache_wave = np.asarray(loaded["coefficients_original"], dtype=float)
        carrier = np.asarray(loaded["carrier"], dtype=float)

    old_control = np.asarray(candidate["selected"]["tangent_coefficients"], dtype=float)
    first_control = np.asarray(first_fit["selected"]["local_patch_coefficients"], dtype=float)
    warm_patch1 = np.asarray(two_patch["selected"]["patch1_coefficients"], dtype=float)
    warm_patch2 = np.asarray(two_patch["selected"]["patch2_coefficients"], dtype=float)
    warm_start_source = TWO_PATCH_PATH.name
    if prior_parent_report and prior_parent_report.get("status") == "completed":
        previous_controls = None
        previous_rounds = prior_parent_report.get("cutting_plane_rounds", [])
        if previous_rounds:
            previous_controls = previous_rounds[-1].get("trial_before_backtrack_controls")
        if previous_controls is None:
            previous_controls = prior_parent_report.get("selected", {}).get("joint_local_coefficients")
        if previous_controls is not None and len(previous_controls) == 360:
            warm_patch1 = np.asarray(previous_controls[:180], dtype=float)
            warm_patch2 = np.asarray(previous_controls[180:], dtype=float)
            warm_start_source = f"{output_path.name}:last_trial_before_backtrack"
    if old_control.shape != (264,) or not np.array_equal(old_control, old_control_cache):
        raise ValueError("Refined cache and enriched candidate use different fixed controls")
    if cache_wave.shape != np.asarray(candidate["selected"]["coefficients_original"]).shape:
        raise ValueError("Refined cache and enriched candidate wave shapes differ")
    if first_control.shape != (180,):
        raise ValueError("First constrained patch does not have 180 controls")
    if warm_patch1.shape != (180,) or warm_patch2.shape != (180,):
        raise ValueError("Frozen two-patch candidate does not have two 180-control blocks")
    if points.shape != (44400, 3) or weights.shape != (44400,):
        raise ValueError(f"Unexpected refined grid shape {points.shape} {weights.shape}")
    if old_design.shape != (133200, 264) or base_residual.shape != (44400, 3):
        raise ValueError(f"Unexpected refined cache array shapes {old_design.shape} {base_residual.shape}")

    fixed_residual = base_residual + (old_design @ old_control).reshape(-1, 3)
    del old_design, base_residual
    first_center = np.asarray(first_rows["support_containment"]["patch_center"], dtype=float)
    first_widths = np.asarray(first_rows["support_containment"]["patch_widths"], dtype=float)
    second_center = np.asarray(second_probe["support_containment"]["patch_center"], dtype=float)
    second_widths = np.asarray(second_probe["support_containment"]["patch_widths"], dtype=float)
    first_design, local_layout = _patch_design(points, first_center, first_widths, carrier)
    second_design, second_layout = _patch_design(points, second_center, second_widths, carrier)
    if first_design.shape != (133200, 180) or second_design.shape != (133200, 180):
        raise ValueError(f"Unexpected local design shapes {first_design.shape} {second_design.shape}")
    if local_layout != second_layout:
        raise ValueError("The two patches do not use the same 180-column layout")

    # The first constrained patch is the peak-cap baseline.  The warm start is
    # the separately frozen two-patch candidate; it is deliberately allowed to
    # be infeasible for the stricter parent drift target below.
    first_residual = fixed_residual + (first_design @ first_control).reshape(-1, 3)
    warm_residual = fixed_residual + (first_design @ warm_patch1).reshape(-1, 3) + (second_design @ warm_patch2).reshape(-1, 3)
    first_metric = _metric(first_residual, weights)
    warm_metric = _metric(warm_residual, weights)
    warm_peak = _max_detail(warm_residual, points, weights)
    first_peak = _max_detail(first_residual, points, weights)
    peak_cap = float(first_metric["momentum_max"] * (1.0 + PEAK_RELATIVE_CAP))
    if not DRIFT_CACHE_PATH.exists() or not DRIFT_REPORT_PATH.exists():
        raise ValueError("Localized drift cache/report is missing")
    with np.load(DRIFT_CACHE_PATH, allow_pickle=False) as drift_loaded:
        drift_gram = np.asarray(drift_loaded["gram"], dtype=float)
        drift_linear = np.asarray(drift_loaded["linear"], dtype=float)
        drift_constant = float(drift_loaded["constant"])
    drift_report = json.loads(DRIFT_REPORT_PATH.read_text(encoding="utf-8"))
    if drift_gram.shape != (360, 360) or drift_linear.shape != (360,):
        raise ValueError("Unexpected localized drift cache shapes")
    drift_cap_parent = float(np.sqrt(max(drift_constant, 0.0)))
    if abs(drift_cap_parent - float(drift_report["parent_drift"])) > 1.0e-8:
        raise ValueError("Parent drift cache/report values disagree")
    drift_cap = drift_cap_parent
    drift_solver_cap = drift_cap_parent

    sqrt_weights = np.sqrt(weights)
    row_weight = np.repeat(sqrt_weights, 3)
    weighted_design_1 = first_design * row_weight[:, None]
    weighted_design_2 = second_design * row_weight[:, None]
    del first_design, second_design
    column_scales_1 = np.maximum(np.linalg.norm(weighted_design_1, axis=0), 1.0e-30)
    column_scales_2 = np.maximum(np.linalg.norm(weighted_design_2, axis=0), 1.0e-30)
    normalized_design_1 = weighted_design_1 / column_scales_1[None, :]
    normalized_design_2 = weighted_design_2 / column_scales_2[None, :]
    del weighted_design_1, weighted_design_2
    # The normalized variable contains both physical patch controls.  Keep
    # the fixed 264 state as the affine objective base, so z_initial carries
    # the frozen two-patch warm start exactly once.
    weighted_base = fixed_residual.reshape(-1) * row_weight
    variable_scale = max(float(np.linalg.norm(warm_residual.reshape(-1) * row_weight)), 1.0)
    objective_scale = variable_scale**2
    z_initial = np.r_[
        warm_patch1 * column_scales_1 / variable_scale,
        warm_patch2 * column_scales_2 / variable_scale,
    ]

    # Build the second patch's exact moment/cone rows using the same reusable
    # machinery that produced localized_constraint_rows.json.
    second_moment_rows, second_cone_rows, second_moment_geometry, second_cone_geometry = (
        _build_constraint_rows(mean_report, geometry, second_center, second_widths)
    )
    first_moment_rows = np.asarray(first_rows["moment_rows"], dtype=float)
    first_cone_rows = np.asarray(first_rows["cone_rows"], dtype=float)
    old_moment_rows = np.asarray(mean_rows_report["moment_rows"], dtype=float)
    old_cone_rows = np.asarray(mean_rows_report["cone_control_rows"], dtype=float)
    if first_moment_rows.shape != (4, 180) or second_moment_rows.shape != (4, 180):
        raise ValueError("Unexpected four-moment row dimensions")
    if first_cone_rows.shape != (81, 180) or second_cone_rows.shape != (81, 180):
        raise ValueError("Unexpected cone row dimensions")
    if old_moment_rows.shape != (4, 64) or old_cone_rows.shape != (81, 64):
        raise ValueError("Unexpected fixed degree-3 row dimensions")

    wave_packed = np.asarray(candidate["selected"]["coefficients_original"], dtype=float)
    wave_x = np.r_[wave_packed[:, 0], wave_packed[:, 1]]
    moment_linearization = moment_report["reusable_moment_linearization"]
    baseline_moments = np.asarray(moment_linearization["baseline_moments"], dtype=float)
    wave_moment_forms = np.asarray(moment_linearization["wave_moment_forms"], dtype=float)
    moment_target = -baseline_moments - np.einsum("i,kij,j->k", wave_x, wave_moment_forms, wave_x)
    old_moment_error = old_moment_rows @ old_control[:64] - moment_target
    cone_baseline = np.asarray(cone_report["cone_baseline"], dtype=float)
    cone_wave_forms = np.asarray(cone_report["cone_wave_forms"], dtype=float)
    cone_lower = np.asarray(cone_report["cone_lower"], dtype=float)
    cone_wave_baseline = cone_baseline + np.einsum("i,kij,j->k", wave_x, cone_wave_forms, wave_x)
    old_cone_margin = cone_wave_baseline + old_cone_rows @ old_control[:64] - cone_lower
    local_moment_rows = np.concatenate((first_moment_rows, second_moment_rows), axis=1)
    local_cone_rows = np.concatenate((first_cone_rows, second_cone_rows), axis=1)

    # Endpoint parent and both analytic local responses.  The normalized
    # endpoint response is held in one 360-column array, so endpoint_eval is
    # exact in the same normalized variables as the momentum objective.
    endpoint = EnrichedEndpointShape(ENDPOINT_CACHE_PATH)
    if endpoint.point_count != 7776:
        raise ValueError("Expected 7776 endpoint points")
    old_endpoint_v = endpoint.velocity_offset + np.einsum(
        "ncq,q->nc", endpoint.velocity_control, old_control
    )
    old_endpoint_j = endpoint.gradient_offset + np.einsum(
        "ncdq,q->ncd", endpoint.gradient_control, old_control
    )
    endpoint_v_1, endpoint_j_1, _, _ = _local_endpoint_response(
        endpoint.points, first_center, first_widths, np.asarray(geometry["carrier"], dtype=float), endpoint.dtau
    )
    endpoint_v_2, endpoint_j_2, _, _ = _local_endpoint_response(
        endpoint.points, second_center, second_widths, np.asarray(geometry["carrier"], dtype=float), endpoint.dtau
    )
    endpoint_v_response = np.concatenate(
        (
            endpoint_v_1 * (variable_scale / column_scales_1)[None, None, :],
            endpoint_v_2 * (variable_scale / column_scales_2)[None, None, :],
        ),
        axis=-1,
    )
    endpoint_j_response = np.concatenate(
        (
            endpoint_j_1 * (variable_scale / column_scales_1)[None, None, None, :],
            endpoint_j_2 * (variable_scale / column_scales_2)[None, None, None, :],
        ),
        axis=-1,
    )
    del endpoint_v_1, endpoint_v_2, endpoint_j_1, endpoint_j_2
    reference = np.asarray(endpoint.reference, dtype=float)
    signs = np.asarray((-1.0, 1.0, np.sign(reference[2])), dtype=float)
    reference_scale = np.maximum(np.abs(reference), 1.0e-30)
    endpoint_requested = np.asarray((1.0e-6, 1.0e-6, 1.0e-3), dtype=float)

    # Row responses in normalized variables.  Rows with zero combined local
    # response remain fixed and are checked directly against the baseline.
    local_matrix_m = variable_scale * np.column_stack(
        (first_moment_rows / column_scales_1[None, :], second_moment_rows / column_scales_2[None, :])
    )
    local_matrix_c = variable_scale * np.column_stack(
        (first_cone_rows / column_scales_1[None, :], second_cone_rows / column_scales_2[None, :])
    )
    moment_row_norm = np.linalg.norm(local_matrix_m, axis=1)
    cone_row_norm = np.linalg.norm(local_matrix_c, axis=1)
    moment_responsive = moment_row_norm > 1.0e-30
    cone_responsive = cone_row_norm > 1.0e-30
    moment_initial_error = old_moment_error + first_moment_rows @ first_control
    cone_initial_margin = old_cone_margin + first_cone_rows @ first_control
    if np.any(~moment_responsive & (np.abs(moment_initial_error) > MOMENT_TOLERANCE)):
        raise ValueError("A fixed moment row is outside tolerance with zero two-patch response")
    if np.any(~cone_responsive & (cone_initial_margin < -CONE_TOLERANCE)):
        raise ValueError("A fixed cone row is infeasible with zero two-patch response")
    moment_matrix_n = local_matrix_m[moment_responsive] / moment_row_norm[moment_responsive, None]
    moment_lower = (-MOMENT_TOLERANCE - old_moment_error[moment_responsive]) / moment_row_norm[moment_responsive]
    moment_upper = (MOMENT_TOLERANCE - old_moment_error[moment_responsive]) / moment_row_norm[moment_responsive]
    cone_matrix_n = local_matrix_c[cone_responsive] / cone_row_norm[cone_responsive, None]
    cone_lower_n = (CONE_OPTIMIZER_SAFETY - old_cone_margin[cone_responsive]) / cone_row_norm[cone_responsive]

    def physical_controls(z):
        z = np.asarray(z, dtype=float)
        return np.r_[variable_scale * z[:180] / column_scales_1,
                     variable_scale * z[180:] / column_scales_2]

    def residual_at(z):
        z = np.asarray(z, dtype=float)
        weighted_correction = (
            normalized_design_1 @ (variable_scale * z[:180])
            + normalized_design_2 @ (variable_scale * z[180:])
        )
        return fixed_residual + (weighted_correction / row_weight).reshape(-1, 3)

    def endpoint_eval(z):
        z = np.asarray(z, dtype=float)
        velocity = old_endpoint_v + np.einsum("ncq,q->nc", endpoint_v_response, z)
        gradient = old_endpoint_j + np.einsum("ncdq,q->ncd", endpoint_j_response, z)
        values, jacobian, _ = _observables_and_jacobian(
            velocity,
            gradient,
            endpoint.weights,
            endpoint.points,
            endpoint_v_response,
            endpoint_j_response,
        )
        margins = signs * (values - reference) / reference_scale - endpoint_requested
        margin_jacobian = signs[:, None] * jacobian / reference_scale[:, None]
        return margins, margin_jacobian, values

    def replay(z):
        z = np.asarray(z, dtype=float)
        controls = physical_controls(z)
        residual = residual_at(z)
        endpoint_margin, _, endpoint_values = endpoint_eval(z)
        drift_square = None
        if drift_gram is not None:
            drift_square = float(
                controls @ drift_gram @ controls + 2.0 * drift_linear @ controls + drift_constant
            )
        return {
            "z": z.copy(),
            "controls": controls,
            "patch1_controls": controls[:180],
            "patch2_controls": controls[180:],
            "residual": residual,
            "metrics": _metric(residual, weights),
            "moment_error": old_moment_error + local_moment_rows @ controls,
            "cone_margin": old_cone_margin + local_cone_rows @ controls,
            "endpoint_margin": endpoint_margin,
            "endpoint_values": endpoint_values,
            "drift_square": drift_square,
        }

    def gate_feasible(state):
        return bool(
            np.max(np.abs(state["moment_error"])) <= MOMENT_TOLERANCE
            and np.min(state["cone_margin"]) >= -CONE_TOLERANCE
            and np.min(state["endpoint_margin"]) >= -ENDPOINT_TOLERANCE
            and state["metrics"]["momentum_max"] <= peak_cap * (1.0 + PEAK_FEASIBILITY_SLACK)
        )

    def feasible(state):
        return bool(
            gate_feasible(state)
            and state["drift_square"] <= drift_cap**2 * (1.0 + PEAK_FEASIBILITY_SLACK)
        )

    initial_state = replay(z_initial)
    if not gate_feasible(initial_state):
        raise ValueError(
            "Frozen two-patch warm start fails a non-drift gate: "
            f"moment={np.max(np.abs(initial_state['moment_error']))!r}, "
            f"cone={np.min(initial_state['cone_margin'])!r}, "
            f"endpoint={np.min(initial_state['endpoint_margin'])!r}, "
            f"peak={initial_state['metrics']['momentum_max']!r}, cap={peak_cap!r}"
        )

    def objective(z):
        z = np.asarray(z, dtype=float)
        weighted_residual = weighted_base + (
            normalized_design_1 @ (variable_scale * z[:180])
            + normalized_design_2 @ (variable_scale * z[180:])
        )
        return 0.5 * float(weighted_residual @ weighted_residual) / objective_scale

    def objective_jac(z):
        z = np.asarray(z, dtype=float)
        weighted_residual = weighted_base + (
            normalized_design_1 @ (variable_scale * z[:180])
            + normalized_design_2 @ (variable_scale * z[180:])
        )
        return variable_scale * np.r_[
            normalized_design_1.T @ weighted_residual,
            normalized_design_2.T @ weighted_residual,
        ] / objective_scale

    objective_direction = np.linspace(-1.0, 1.0, 360, dtype=float)
    objective_direction /= np.linalg.norm(objective_direction)
    objective_fd_step = 1.0e-4
    objective_gradient_analytic = float(objective_jac(z_initial) @ objective_direction)
    objective_gradient_fd = float(
        (objective(z_initial + objective_fd_step * objective_direction)
         - objective(z_initial - objective_fd_step * objective_direction))
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

    def endpoint_fun(z):
        return endpoint_eval(z)[0]

    def endpoint_jac(z):
        return endpoint_eval(z)[1]

    control_per_variable = np.r_[
        variable_scale / column_scales_1,
        variable_scale / column_scales_2,
    ]

    def drift_fun(z):
        if drift_gram is None:
            return np.ones(1, dtype=float)
        controls = physical_controls(z)
        square = float(controls @ drift_gram @ controls + 2.0 * drift_linear @ controls + drift_constant)
        return np.asarray(((drift_solver_cap**2 - square) / drift_solver_cap**2,), dtype=float)

    def drift_jac(z):
        if drift_gram is None:
            return np.zeros((1, 360), dtype=float)
        controls = physical_controls(z)
        gradient_controls = 2.0 * (drift_gram @ controls + drift_linear)
        return (-gradient_controls * control_per_variable / drift_solver_cap**2)[None, :]

    def retract_toward_zero(z):
        """Check the full quadratic line for a feasible target-cap interval."""

        controls = physical_controls(z)
        quadratic = float(controls @ drift_gram @ controls)
        linear = float(2.0 * drift_linear @ controls)
        constant = float(drift_constant - drift_cap**2)
        if quadratic <= 0.0:
            return None
        discriminant = linear * linear - 4.0 * quadratic * constant
        if discriminant <= 0.0:
            return None
        root = float(np.sqrt(discriminant))
        lower = (-linear - root) / (2.0 * quadratic)
        upper = (-linear + root) / (2.0 * quadratic)
        lower = max(0.0, lower)
        upper = min(1.0, upper)
        if upper <= lower:
            return None
        lam = float(upper * (1.0 - 1.0e-8))
        if lam <= lower:
            lam = float(lower * (1.0 + 1.0e-8))
        repaired_z = lam * np.asarray(z, dtype=float)
        repaired_state = replay(repaired_z)
        return lam, repaired_z, repaired_state, {
            "quadratic": quadratic,
            "linear": linear,
            "constant_relative_to_cap": constant,
            "discriminant": discriminant,
            "lower_root": float(lower),
            "upper_root": float(upper),
        }

    def pool_data(z, indices):
        indices = np.asarray(sorted(indices), dtype=int)
        rows = np.concatenate([np.arange(3 * index, 3 * index + 3) for index in indices])
        point_matrix = variable_scale * np.column_stack((
            normalized_design_1[rows], normalized_design_2[rows]
        )) / row_weight[rows, None]
        values = fixed_residual[indices] + (point_matrix @ np.asarray(z, dtype=float)).reshape(len(indices), 3)
        normalized_margin = 1.0 - np.sum(values * values, axis=1) / peak_cap**2
        point_blocks = point_matrix.reshape(len(indices), 3, -1)
        jacobian = np.empty((len(indices), 360), dtype=float)
        for row, (value, block) in enumerate(zip(values, point_blocks)):
            jacobian[row] = -2.0 * (value @ block) / peak_cap**2
        return normalized_margin, jacobian

    def full_peak(z):
        return np.linalg.norm(residual_at(z), axis=1)

    pool = set(np.argsort(np.linalg.norm(warm_residual, axis=1))[-POOL_INITIAL:].tolist())
    warm_drift = float(np.sqrt(max(initial_state["drift_square"], 0.0)))
    cap_schedule = []
    for value in (warm_drift * (1.0 + 1.0e-6), 1950.0, 1850.0, TARGET_DRIFT_CAP):
        value = float(max(value, TARGET_DRIFT_CAP))
        if not cap_schedule or value < cap_schedule[-1] - 1.0e-8:
            cap_schedule.append(value)
    if int(max_rounds) == 1:
        cap_schedule = [TARGET_DRIFT_CAP]
    else:
        cap_schedule = cap_schedule[:max(1, int(max_rounds))]
    rounds = []
    stage_results = []
    stage_start_z = z_initial.copy()
    target_best = None
    last_trial = initial_state
    best_reason = "no_target_cap_feasible_improving_candidate"

    def save_stage_progress():
        progress = {
            "status": "running",
            "target_drift_cap": TARGET_DRIFT_CAP,
            "parent_drift_cap_reference": drift_cap_parent,
            "warm_start_drift": warm_drift,
            "completed_stages": stage_results,
            "source": Path(__file__).name,
        }
        output_path.write_text(json.dumps(progress, indent=2) + "\n", encoding="utf-8")

    for stage_index, stage_cap in enumerate(cap_schedule):
        drift_cap = float(stage_cap)
        drift_solver_cap = max(0.0, float(stage_cap) - DRIFT_SOLVER_SAFETY)
        indices = sorted(pool)

        def pool_fun(z):
            return pool_data(z, indices)[0]

        def pool_jac(z):
            return pool_data(z, indices)[1]

        constraints = [
            {"type": "ineq", "fun": endpoint_fun, "jac": endpoint_jac},
            {"type": "ineq", "fun": pool_fun, "jac": pool_jac},
            {"type": "ineq", "fun": drift_fun, "jac": drift_jac},
        ]
        if np.any(moment_responsive):
            constraints.append(LinearConstraint(moment_matrix_n, moment_lower, moment_upper))
        if np.any(cone_responsive):
            constraints.append(LinearConstraint(cone_matrix_n, cone_lower_n, np.inf))
        fit = minimize(
            objective,
            stage_start_z,
            jac=objective_jac,
            method="SLSQP",
            constraints=constraints,
            options={"maxiter": 100, "ftol": 1.0e-10, "disp": False},
        )
        solver_candidate_z = np.asarray(fit.x if fit.x is not None else stage_start_z, dtype=float)
        solver_candidate_state = replay(solver_candidate_z)
        candidate_z = solver_candidate_z.copy()
        candidate_state = solver_candidate_state
        backtrack_alpha = None
        backtrack_steps = 0
        retraction_lambda = None
        retraction_diagnostic = None
        retraction_state = None
        if gate_feasible(candidate_state) and not feasible(candidate_state):
            retraction = retract_toward_zero(candidate_z)
            if retraction is not None:
                retraction_lambda, repaired_z, repaired_state, retraction_diagnostic = retraction
                if gate_feasible(repaired_state) and feasible(repaired_state):
                    candidate_z = repaired_z
                    candidate_state = repaired_state
                    retraction_state = repaired_state
        stage_start_state = replay(stage_start_z)
        if retraction_state is None and gate_feasible(candidate_state) and not feasible(candidate_state) and gate_feasible(stage_start_state) and feasible(stage_start_state):
            delta_z = candidate_z - stage_start_z
            lo = 0.0
            hi = 1.0
            candidate_state = stage_start_state
            for backtrack_steps in range(1, 61):
                alpha = 0.5 * (lo + hi)
                probe_z = stage_start_z + alpha * delta_z
                probe_state = replay(probe_z)
                if feasible(probe_state):
                    lo = alpha
                    candidate_z = probe_z
                    candidate_state = probe_state
                else:
                    hi = alpha
            if lo > 0.0:
                backtrack_alpha = float(lo)
        magnitudes = full_peak(candidate_z)
        violating = np.flatnonzero(magnitudes > peak_cap * (1.0 + PEAK_FEASIBILITY_SLACK))
        worst = np.argsort(magnitudes)[-POOL_ADD:][::-1]
        added = [int(index) for index in worst if int(index) not in pool]
        candidate_is_target = bool(stage_index == len(cap_schedule) - 1 and abs(stage_cap - TARGET_DRIFT_CAP) <= 1.0e-8)
        improving_vs_first = bool(candidate_state["metrics"]["momentum_volume_L2"] < first_metric["momentum_volume_L2"])
        if candidate_is_target and feasible(candidate_state) and improving_vs_first:
            target_best = candidate_state
            best_reason = "target_cap_feasible_lower_L2_than_first_patch"
        if feasible(candidate_state):
            stage_start_z = candidate_z.copy()
        else:
            stage_start_z = candidate_z.copy()
        last_trial = candidate_state
        for index in added:
            pool.add(index)
        rounds.append({
            "stage": stage_index + 1,
            "stage_cap": float(stage_cap),
            "solver_drift_cap": float(drift_solver_cap),
            "parent_target_stage": candidate_is_target,
            "pool_size_before": len(indices),
            "pool_size_after": len(pool),
            "optimizer_success": bool(fit.success),
            "optimizer_status": int(fit.status),
            "optimizer_message": str(fit.message),
            "optimizer_iterations": int(getattr(fit, "nit", -1)),
            "candidate_gate_feasible": bool(gate_feasible(candidate_state)),
            "candidate_feasible_at_stage_cap": bool(feasible(candidate_state)),
            "candidate_improving_vs_first": improving_vs_first,
            "candidate_metrics": candidate_state["metrics"],
            "candidate_moment_max_abs": float(np.max(np.abs(candidate_state["moment_error"]))),
            "candidate_cone_min_margin": float(np.min(candidate_state["cone_margin"])),
            "candidate_endpoint_min_margin": float(np.min(candidate_state["endpoint_margin"])),
            "candidate_global_peak": float(np.max(magnitudes)),
            "candidate_drift": float(np.sqrt(max(candidate_state["drift_square"], 0.0))),
            "trial_before_backtrack_feasible": bool(feasible(solver_candidate_state)),
            "backtrack_alpha_from_stage_start": backtrack_alpha,
            "backtrack_bisection_steps": int(backtrack_steps),
            "retraction_lambda_toward_zero": retraction_lambda,
            "retraction_quadratic_diagnostic": retraction_diagnostic,
            "trial_before_backtrack_metrics": solver_candidate_state["metrics"],
            "trial_before_backtrack_drift": float(np.sqrt(max(solver_candidate_state["drift_square"], 0.0))),
            "trial_before_backtrack_controls": solver_candidate_state["controls"].tolist(),
            "peak_cap": float(peak_cap),
            "global_peak_violating_point_count": int(len(violating)),
            "newly_added_worst_points": added,
            "worst_point": int(np.argmax(magnitudes)),
        })
        stage_results.append({
            "stage": stage_index + 1,
            "cap": float(stage_cap),
            "candidate_feasible": bool(feasible(candidate_state)),
            "candidate_gate_feasible": bool(gate_feasible(candidate_state)),
            "candidate_drift": float(np.sqrt(max(candidate_state["drift_square"], 0.0))),
            "candidate_l2": candidate_state["metrics"]["momentum_volume_L2"],
            "candidate_max": candidate_state["metrics"]["momentum_max"],
        })
        save_stage_progress()
        if candidate_is_target:
            break

    selected = target_best if target_best is not None else last_trial
    selected_feasible = bool(target_best is not None and feasible(target_best))
    source_paths = {
        "script": Path(__file__),
        "candidate": CANDIDATE_PATH,
        "first_fit": FIRST_FIT_PATH,
        "first_rows": FIRST_ROWS_PATH,
        "second_probe": SECOND_PROBE_PATH,
        "refined_cache": CACHE_PATH,
        "refined_cache_report": CACHE_REPORT_PATH,
        "mean_rows": MEAN_ROWS_PATH,
        "moment": MOMENT_PATH,
        "cone": CONE_PATH,
        "endpoint_cache": ENDPOINT_CACHE_PATH,
        "endpoint_cache_report": ENDPOINT_CACHE_REPORT_PATH,
        "mean": MEAN_PATH,
        "dynamic": DYNAMIC_PATH,
        "geometry": GEOMETRY_PATH,
        "two_patch_candidate": TWO_PATCH_PATH,
        "drift_cache": DRIFT_CACHE_PATH,
        "drift_report": DRIFT_REPORT_PATH,
    }
    report = {
        "status": "completed",
        "accepted": False,
        "pde_validated": False,
        "scale_recursion_established": False,
        "constraints_maintained": False,
        "assembled_feasible": bool(selected_feasible),
        "assembled_constraints_maintained": bool(selected_feasible),
        "source": Path(__file__).name,
        "source_sha256": _sha256(Path(__file__)),
        "sources": {name: _source_entry(path) for name, path in source_paths.items()},
        "scope": (
            "Reference-time 1800-profile drift continuation on the frozen 44,400-point grid. "
            "The enriched 264 controls are fixed; the two local 180-control blocks warm-start "
            "from the frozen two-patch candidate. Explicit drift caps continue to the parent "
            "target while four moments, 81 cones, endpoint shape inequalities, and a peak cap "
            "relative to the first constrained patch are checked. No mean replay, holdout, "
            "finite-time, PDE, or recursive-scale acceptance."
        ),
        "inputs": {
            "point_count": int(len(points)),
            "candidate_control_count": 264,
            "local_control_count_per_patch": 180,
            "joint_local_control_count": 360,
            "first_patch_center": first_center.tolist(),
            "first_patch_widths": first_widths.tolist(),
            "second_patch_center": second_center.tolist(),
            "second_patch_widths": second_widths.tolist(),
            "patch_degree": 2,
            "patch_modes": [0, 1, 2],
            "old_controls_held_fixed": True,
            "initialization": warm_start_source,
            "prior_parent_report_sha256": prior_parent_hash,
            "column_normalization": "weighted L2 norm of each patch momentum-response column",
            "variable_normalization": "joint normalized variables divided by frozen two-patch weighted L2",
            "objective": "normalized weighted residual L2 using frozen enriched base and analytic two-patch responses",
            "peak_cap": float(peak_cap),
            "peak_cap_relative_to_first_constrained": PEAK_RELATIVE_CAP,
            "peak_cutting_plane_round_limit": int(max_rounds),
            "parent_drift_cap_reference": float(drift_cap_parent),
            "target_drift_cap": TARGET_DRIFT_CAP,
            "drift_solver_safety_absolute": DRIFT_SOLVER_SAFETY,
            "continuation_caps": [float(value) for value in cap_schedule],
            "tolerances": {
                "moment_max_abs": MOMENT_TOLERANCE,
                "cone_min_margin": CONE_TOLERANCE,
                "cone_optimizer_safety": CONE_OPTIMIZER_SAFETY,
                "endpoint_min_margin": ENDPOINT_TOLERANCE,
                "peak_feasibility_slack": PEAK_FEASIBILITY_SLACK,
            },
            "endpoint_requested_signed_fractional_change": endpoint_requested.tolist(),
        },
        "warm_start_state": {
            "gate_feasible": bool(gate_feasible(initial_state)),
            "parent_cap_feasible": bool(feasible(initial_state)),
            "moment_max_abs": float(np.max(np.abs(initial_state["moment_error"]))),
            "cone_min_margin": float(np.min(initial_state["cone_margin"])),
            "cone_pass_count": int(np.sum(initial_state["cone_margin"] >= -CONE_TOLERANCE)),
            "endpoint_margin": initial_state["endpoint_margin"].tolist(),
            "endpoint_min_margin": float(np.min(initial_state["endpoint_margin"])),
            "momentum": initial_state["metrics"],
            "peak_cap": float(peak_cap),
            "drift_square": initial_state["drift_square"],
            "drift": float(np.sqrt(max(initial_state["drift_square"], 0.0))) if initial_state["drift_square"] is not None else None,
        },
        "selected": {
            "reason": best_reason,
            "feasible": bool(selected_feasible),
            "patch1_coefficients": selected["patch1_controls"].tolist(),
            "patch2_coefficients": selected["patch2_controls"].tolist(),
            "joint_local_coefficients": selected["controls"].tolist(),
            "patch1_coefficient_norm": float(np.linalg.norm(selected["patch1_controls"])),
            "patch2_coefficient_norm": float(np.linalg.norm(selected["patch2_controls"])),
            "joint_coefficient_norm": float(np.linalg.norm(selected["controls"])),
            "momentum": selected["metrics"],
            "peak_detail": _max_detail(selected["residual"], points, weights),
            "momentum_l2_improvement_vs_first": float(first_metric["momentum_volume_L2"] - selected["metrics"]["momentum_volume_L2"]),
            "momentum_max_change_vs_first": float(selected["metrics"]["momentum_max"] - first_metric["momentum_max"]),
            "moment_max_abs": float(np.max(np.abs(selected["moment_error"]))),
            "cone_min_margin": float(np.min(selected["cone_margin"])),
            "cone_pass_count": int(np.sum(selected["cone_margin"] >= -CONE_TOLERANCE)),
            "endpoint_values": selected["endpoint_values"].tolist(),
            "endpoint_margin": selected["endpoint_margin"].tolist(),
            "endpoint_min_margin": float(np.min(selected["endpoint_margin"])),
            "drift_square": selected["drift_square"],
            "drift": float(np.sqrt(max(selected["drift_square"], 0.0))) if selected["drift_square"] is not None else None,
            "old_controls_unchanged": True,
            "active_constraints": {
                "moment_rows_within_1e-8_of_tolerance": [
                    int(index) for index, value in enumerate(selected["moment_error"])
                    if abs(value) >= MOMENT_TOLERANCE - 1.0e-8
                ],
                "cone_rows_within_1e-8_of_tolerance": [
                    int(index) for index, value in enumerate(selected["cone_margin"])
                    if value <= -CONE_TOLERANCE + 1.0e-8
                ],
                "endpoint_rows_within_1e-8": [
                    int(index) for index, value in enumerate(selected["endpoint_margin"])
                    if value <= 1.0e-8
                ],
                "peak_cap_fractional_slack": float(peak_cap / selected["metrics"]["momentum_max"] - 1.0),
            },
        },
        "first_constrained_before": {
            "momentum": first_metric,
            "peak_detail": first_peak,
            "moment_max_abs": float(np.max(np.abs(moment_initial_error))),
            "cone_min_margin": float(np.min(cone_initial_margin)),
            "cone_pass_count": int(np.sum(cone_initial_margin >= -CONE_TOLERANCE)),
        },
        "drift_constraint": {
            "enabled": True,
            "cap": TARGET_DRIFT_CAP,
            "selected_drift": float(np.sqrt(max(selected["drift_square"], 0.0))) if selected["drift_square"] is not None else None,
            "selected_slack": float(TARGET_DRIFT_CAP - np.sqrt(max(selected["drift_square"], 0.0))),
            "parent_cap_reference": float(drift_cap_parent),
            "warm_start_drift": warm_drift,
            "target_feasible_improving_candidate_found": bool(target_best is not None),
            "last_trial_target_cap_feasible": bool(feasible(last_trial)),
        },
        "cutting_plane_rounds": rounds,
        "continuation_stage_results": stage_results,
        "row_diagnostics": {
            "first_moment_shape": list(first_moment_rows.shape),
            "second_moment_shape": list(second_moment_rows.shape),
            "first_cone_shape": list(first_cone_rows.shape),
            "second_cone_shape": list(second_cone_rows.shape),
            "first_moment_max_abs": float(np.max(np.abs(first_moment_rows))),
            "second_moment_max_abs": float(np.max(np.abs(second_moment_rows))),
            "first_cone_max_abs": float(np.max(np.abs(first_cone_rows))),
            "second_cone_max_abs": float(np.max(np.abs(second_cone_rows))),
            "first_moment_rank": int(np.linalg.matrix_rank(first_moment_rows)),
            "second_moment_rank": int(np.linalg.matrix_rank(second_moment_rows)),
            "joint_moment_responsive_rows": int(np.sum(moment_responsive)),
            "joint_cone_responsive_rows": int(np.sum(cone_responsive)),
            "second_moment_geometry": second_moment_geometry,
            "second_cone_geometry": second_cone_geometry,
        },
        "matrix_diagnostics": {
            "patch_design_shape": [133200, 180],
            "joint_design_shape": [133200, 360],
            "patch1_column_scale_min": float(np.min(column_scales_1)),
            "patch1_column_scale_max": float(np.max(column_scales_1)),
            "patch2_column_scale_min": float(np.min(column_scales_2)),
            "patch2_column_scale_max": float(np.max(column_scales_2)),
            "variable_scale_warm_start_weighted_L2": float(variable_scale),
            "objective_gradient_check": objective_gradient_check,
            "endpoint_response_shape": list(endpoint_v_response.shape),
            "endpoint_gradient_response_shape": list(endpoint_j_response.shape),
        },
        "elapsed_seconds": float(time.perf_counter() - started),
    }
    output_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": report["status"],
        "selected_feasible": selected_feasible,
        "best_reason": best_reason,
        "first_l2": first_metric["momentum_volume_L2"],
        "selected_l2": selected["metrics"]["momentum_volume_L2"],
        "first_max": first_metric["momentum_max"],
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
