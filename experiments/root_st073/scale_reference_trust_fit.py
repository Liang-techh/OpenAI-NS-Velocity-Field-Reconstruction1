"""Trust-region least-squares fit for the frozen reference-velocity step.

This is a bounded numerical follow-up to ``scale_reference_velocity_step``.
It rebuilds that module's 360 columns (270 velocity and 90 pressure), but
solves the weighted linear least-squares problem inside the one-percent
velocity trust region instead of scaling an unconstrained direction.  The
pressure variables are eliminated with an SVD, and the velocity trust metric
is reduced with QR/SVD factors.  The selected result is then checked with the
exact quadratic convection remainder and a short nonlinear backtrack.

The report is a spatial diagnostic.  It does not establish support, moments,
cones, shape, finite-time evolution, PDE validity, or scale recursion.
"""

from __future__ import annotations

import argparse
import gc
import hashlib
import json
import os
import sys
import time
from pathlib import Path

for _key in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS"):
    os.environ[_key] = "1"

import numpy as np

ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import scale_reference_velocity_step as previous  # noqa: E402


OUTPUT_PATH = ROOT / "scale_reference_trust_fit.json"
PREVIOUS_REPORT_PATH = ROOT / "scale_reference_velocity_step.json"
RIDGE_RELATIVE = previous.RIDGE_RELATIVE
DEGREE = previous.DEGREE
MODES = previous.MODES
VELOCITY_COLUMNS_PER_PATCH = previous.VELOCITY_COLUMNS_PER_PATCH
PRESSURE_COLUMNS_PER_PATCH = previous.PRESSURE_COLUMNS_PER_PATCH
COLUMNS_PER_PATCH = previous.COLUMNS_PER_PATCH
PATCH_COUNT = 2
JOINT_COLUMNS = PATCH_COUNT * COLUMNS_PER_PATCH
VELOCITY_COLUMNS = PATCH_COUNT * VELOCITY_COLUMNS_PER_PATCH
PRESSURE_COLUMNS = PATCH_COUNT * PRESSURE_COLUMNS_PER_PATCH
TRUST_FRACTION = previous.VELOCITY_TRUST_FRACTION


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _array_sha256(*arrays) -> str:
    digest = hashlib.sha256()
    for array in arrays:
        value = np.ascontiguousarray(np.asarray(array))
        digest.update(str(value.dtype).encode("ascii"))
        digest.update(np.asarray(value.shape, dtype=np.int64).tobytes())
        digest.update(value.tobytes())
    return digest.hexdigest()


def _save(path: Path, report: dict) -> None:
    path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")


def _metric(values, weights, points):
    values = np.asarray(values, dtype=float)
    weights = np.asarray(weights, dtype=float)
    points = np.asarray(points, dtype=float)
    norms = np.linalg.norm(values, axis=1)
    weighted_square = float(np.sum(weights * np.sum(values * values, axis=1)))
    index = int(np.argmax(norms))
    return {
        "point_count": int(len(values)),
        "volume_L2": float(np.sqrt(max(weighted_square, 0.0))),
        "volume_RMS": float(np.sqrt(max(weighted_square, 0.0) / max(float(np.sum(weights)), 1.0e-300))),
        "max_norm": float(norms[index]),
        "max_point": points[index].tolist(),
        "max_index": index,
        "max_vector": values[index].tolist(),
        "physical_volume": float(np.sum(weights)),
    }


def _weighted_column_norms(array):
    return np.sqrt(np.sum(np.asarray(array, dtype=float) ** 2, axis=0))


def _rank_from_singular_values(singular_values):
    singular_values = np.asarray(singular_values, dtype=float)
    if singular_values.size == 0 or singular_values[0] <= 0.0:
        return 0
    cutoff = RIDGE_RELATIVE * max(float(singular_values[0]), 1.0e-300)
    return int(np.sum(singular_values > cutoff))


def _spectral_solution(u, singular_values, vh, projected_rhs, multiplier):
    """Return y for min ||A y-b||^2 + multiplier ||y||^2 via SVD."""

    singular_values = np.asarray(singular_values, dtype=float)
    projected_rhs = np.asarray(projected_rhs, dtype=float)
    # The caller supplies U^T b.  Zero singular directions are discarded
    # rather than divided by a small number.
    cutoff = RIDGE_RELATIVE * max(float(singular_values[0]), 1.0e-300) if singular_values.size else 0.0
    keep = singular_values > cutoff
    if not np.any(keep):
        return np.zeros(vh.shape[1], dtype=float)
    coefficients = np.zeros_like(singular_values)
    coefficients[keep] = singular_values[keep] * projected_rhs[keep] / (
        singular_values[keep] ** 2 + float(multiplier)
    )
    return vh.T @ coefficients


def _solve_multiplier(a_tilde, rhs, trust_bound):
    """Solve the Euclidean trust-region problem with QR followed by SVD."""

    # QR keeps the large-row problem out of normal equations.  Only its thin
    # R factor is needed for the SVD, while Q is used once to form Q^T rhs.
    q, r = np.linalg.qr(a_tilde, mode="reduced")
    u_small, singular_values, vh = np.linalg.svd(r, full_matrices=False)
    projected_rhs = u_small.T @ (q.T @ rhs)
    del q, r, u_small
    unconstrained = _spectral_solution(None, singular_values, vh, projected_rhs, 0.0)
    unconstrained_norm = float(np.linalg.norm(unconstrained))
    trust_bound = float(trust_bound)
    if unconstrained_norm <= trust_bound * (1.0 + 1.0e-12):
        multiplier = 0.0
        solution = unconstrained
        iterations = 0
    else:
        low = 0.0
        high = max(float(singular_values[0] ** 2), 1.0)

        def at(value):
            trial = _spectral_solution(None, singular_values, vh, projected_rhs, value)
            return trial, float(np.linalg.norm(trial))

        solution, norm_high = at(high)
        while norm_high > trust_bound:
            high *= 4.0
            solution, norm_high = at(high)
            if high > 1.0e300:
                raise FloatingPointError("trust-region multiplier could not be bracketed")
        iterations = 0
        for iterations in range(1, 81):
            middle = 0.5 * (low + high)
            trial, norm_middle = at(middle)
            if norm_middle > trust_bound:
                low = middle
            else:
                high = middle
                solution = trial
            if abs(norm_middle - trust_bound) <= max(1.0e-12 * trust_bound, 1.0e-14):
                solution = trial
                high = middle
                break
        multiplier = high
        solution, _ = at(multiplier)
    return {
        "coordinates": solution,
        "multiplier": float(multiplier),
        "unconstrained_norm": unconstrained_norm,
        "constrained_norm": float(np.linalg.norm(solution)),
        "iterations": int(iterations),
        "singular_values": singular_values,
        "vh": vh,
        "projected_rhs": projected_rhs,
    }


def _layout_and_patch_data(data, pressure_design):
    points = data["points"]
    drift_inputs = data["drift_report"]["inputs"]
    pressure_meta = {
        "first_center": np.asarray(drift_inputs["first_patch_center"], dtype=float),
        "first_widths": np.asarray(drift_inputs["first_patch_widths"], dtype=float),
        "second_center": np.asarray(drift_inputs["second_patch_center"], dtype=float),
        "second_widths": np.asarray(drift_inputs["second_patch_widths"], dtype=float),
    }
    patch_data = []
    joint_columns = []
    joint_layout = []
    centers_widths = (
        (pressure_meta["first_center"], pressure_meta["first_widths"]),
        (pressure_meta["second_center"], pressure_meta["second_widths"]),
    )
    for patch_index, (center, widths) in enumerate(centers_widths, 1):
        values_and_jets = previous._mode_velocity_columns(
            points,
            center,
            widths,
            data["carrier"],
            data["nu"],
            data["h"],
            data["velocity"],
            data["jacobian"],
            data["laplacian"],
            data["tau0"],
        )
        velocity_response, values, jacobians, laplacians, mean_swirl, velocity_layout = values_and_jets
        pstart = (patch_index - 1) * PRESSURE_COLUMNS_PER_PATCH
        pend = pstart + PRESSURE_COLUMNS_PER_PATCH
        pressure_response = pressure_design[:, pstart:pend]
        joint_columns.extend((velocity_response, pressure_response))
        joint_layout.extend({"patch": patch_index, **entry} for entry in velocity_layout)
        joint_layout.extend({"patch": patch_index, **entry} for entry in data["pressure_report"]["pressure_layout"]["per_patch"])
        patch_data.append((velocity_response, values, jacobians, laplacians, mean_swirl))
        print(json.dumps({"stage": "patch_basis", "patch": patch_index, "velocity_columns": int(velocity_response.shape[1]), "pressure_columns": int(pressure_response.shape[1])}), flush=True)
    return patch_data, joint_columns, joint_layout, pressure_meta


def _velocity_value_matrix(patch_data):
    return np.column_stack(
        [patch_data[0][1][:, :, j].reshape(-1) for j in range(VELOCITY_COLUMNS_PER_PATCH)]
        + [patch_data[1][1][:, :, j].reshape(-1) for j in range(VELOCITY_COLUMNS_PER_PATCH)]
    )


def _coefficients_by_patch(coefficients):
    coefficients = np.asarray(coefficients, dtype=float)
    return (
        coefficients[:VELOCITY_COLUMNS_PER_PATCH],
        coefficients[COLUMNS_PER_PATCH:COLUMNS_PER_PATCH + VELOCITY_COLUMNS_PER_PATCH],
    )


def _trial(alpha, baseline, linear, quadratic, patch_data, coefficients, weights, points):
    alpha = float(alpha)
    trial_coefficients = alpha * np.asarray(coefficients, dtype=float)
    velocity, jacobian = previous._joint_velocity_from_coefficients(
        patch_data, _coefficients_by_patch(trial_coefficients)
    )
    residual = baseline + alpha * linear + alpha * alpha * quadratic
    return {
        "alpha": alpha,
        "coefficients": trial_coefficients,
        "velocity": velocity,
        "jacobian": jacobian,
        "remainder": alpha * alpha * quadratic,
        "residual": residual,
        "residual_metric": _metric(residual, weights, points),
        "velocity_metric": _metric(velocity, weights, points),
        "remainder_metric": _metric(alpha * alpha * quadratic, weights, points),
    }


def run(output_path=OUTPUT_PATH):
    started = time.perf_counter()
    output_path = Path(output_path)
    data = previous._load_inputs()
    points = data["points"]
    weights = data["weights"]
    baseline_generator = data["generator_residual"]
    baseline, pressure_correction, pressure_design, pressure_meta = previous._pressure_baseline(data)
    report = {
        "status": "initializing",
        "accepted": False,
        "pde_validated": False,
        "scale_recursion_established": False,
        "constraints_maintained": False,
        "scope": (
            "Bounded weighted least-squares reference-velocity trust-region diagnostic. "
            "The same two degree-2 compact patches, modes 0/1/2, 270 velocity columns, "
            "and 90 pressure-gradient columns as scale_reference_velocity_step are rebuilt. "
            "Pressure is eliminated by SVD and the velocity change is constrained to one "
            "percent of the frozen reference velocity norm. Exact quadratic convection and "
            "a nonlinear backtrack are reported. This is not a PDE, trajectory, support, "
            "moment/cone, or scale-recursion acceptance test."
        ),
        "sources": {
            "previous_script": {"path": previous.__file__, "sha256": _sha256(Path(previous.__file__))},
            "previous_report": {"path": PREVIOUS_REPORT_PATH.name, "sha256": _sha256(PREVIOUS_REPORT_PATH)},
            "generator_report": {"path": previous.GENERATOR_REPORT_PATH.name, "sha256": _sha256(previous.GENERATOR_REPORT_PATH)},
            "generator_cache": {"path": previous.GENERATOR_CACHE_PATH.name, "sha256": _sha256(previous.GENERATOR_CACHE_PATH)},
            "pressure_report": {"path": previous.PRESSURE_REPORT_PATH.name, "sha256": _sha256(previous.PRESSURE_REPORT_PATH)},
            "drift_report": {"path": previous.DRIFT_REPORT_PATH.name, "sha256": _sha256(previous.DRIFT_REPORT_PATH)},
            "frozen_geometry": {"path": previous.FROZEN_PATH.name, "sha256": _sha256(previous.FROZEN_PATH)},
            "trust_fit_script": {"path": Path(__file__).name, "sha256": _sha256(Path(__file__))},
        },
        "inputs": {
            "point_count": int(len(points)),
            "physical_volume": float(np.sum(weights)),
            "grid_array_sha256": _array_sha256(points, weights),
            "tau0": data["tau0"],
            "similarity_exponent_h": data["h"],
            "viscosity": data["nu"],
            "degree": DEGREE,
            "modes": list(MODES),
            "patch_count": PATCH_COUNT,
            "velocity_columns_per_patch": VELOCITY_COLUMNS_PER_PATCH,
            "pressure_columns_per_patch": PRESSURE_COLUMNS_PER_PATCH,
            "joint_column_count": JOINT_COLUMNS,
            "velocity_column_count": VELOCITY_COLUMNS,
            "pressure_column_count": PRESSURE_COLUMNS,
            "velocity_trust_fraction": TRUST_FRACTION,
            "ridge_relative_cutoff": RIDGE_RELATIVE,
            "carrier": data["carrier"].tolist(),
            "correction_kind": "initial reference curl-potential velocity; time-tangent coefficients are not reused",
            "optional_wave_amplitude_column": {"included": False, "reason": "primary same-360-column fit requested; optional direction deferred"},
        },
        "patches": {
            "first": {"center": pressure_meta["first_center"].tolist(), "widths": pressure_meta["first_widths"].tolist()},
            "second": {"center": pressure_meta["second_center"].tolist(), "widths": pressure_meta["second_widths"].tolist()},
            "joint_order": "patch 1 [135 velocity, 45 pressure], then patch 2 [135 velocity, 45 pressure]",
        },
        "pressure_baseline": {
            "before": previous._metric(baseline_generator, weights, points),
            "after": previous._metric(baseline, weights, points),
            "correction": previous._metric(pressure_correction, weights, points),
        },
    }
    _save(output_path, report)
    print(json.dumps({"stage": "sources_loaded", "point_count": len(points)}), flush=True)

    patch_data, joint_columns, joint_layout, _ = _layout_and_patch_data(data, pressure_design)
    design = np.column_stack(joint_columns)
    if design.shape != (3 * len(points), JOINT_COLUMNS):
        raise ValueError(f"Unexpected joint response shape {design.shape}")
    row_weight = np.repeat(np.sqrt(weights), 3)
    design *= row_weight[:, None]
    rhs = -baseline.reshape(-1) * row_weight
    velocity_indices = np.r_[np.arange(VELOCITY_COLUMNS_PER_PATCH), np.arange(COLUMNS_PER_PATCH, COLUMNS_PER_PATCH + VELOCITY_COLUMNS_PER_PATCH)]
    pressure_indices = np.r_[np.arange(VELOCITY_COLUMNS_PER_PATCH, COLUMNS_PER_PATCH), np.arange(COLUMNS_PER_PATCH + VELOCITY_COLUMNS_PER_PATCH, JOINT_COLUMNS)]
    av = design[:, velocity_indices]
    ap = design[:, pressure_indices]
    value_matrix = _velocity_value_matrix(patch_data)
    weighted_values = value_matrix * row_weight[:, None]
    velocity_scales = np.maximum(_weighted_column_norms(av), 1.0e-30)
    normalized_velocity_metric = weighted_values / velocity_scales[None, :]
    reference_velocity_norm = float(np.sqrt(np.sum(weights * np.sum(data["velocity"] * data["velocity"], axis=1))))
    trust_bound = TRUST_FRACTION * reference_velocity_norm

    # Eliminate pressure by an orthogonal projection. This avoids forming any
    # normal equations for the nearly dependent compact pressure columns.
    pressure_u, pressure_s, pressure_vh = np.linalg.svd(ap, full_matrices=False)
    pressure_rank = _rank_from_singular_values(pressure_s)
    if pressure_rank == 0:
        raise np.linalg.LinAlgError("pressure design has zero numerical rank")
    pressure_u = pressure_u[:, :pressure_rank]
    pressure_s = pressure_s[:pressure_rank]
    pressure_vh = pressure_vh[:pressure_rank, :]
    pressure_rhs_projection = pressure_u.T @ rhs
    projected_rhs = rhs - pressure_u @ pressure_rhs_projection
    av_projected = av - pressure_u @ (pressure_u.T @ av)
    trust_r = np.linalg.qr(normalized_velocity_metric, mode="r")
    trust_u, trust_s, trust_vh = np.linalg.svd(trust_r, full_matrices=False)
    trust_rank = _rank_from_singular_values(trust_s)
    if trust_rank == 0:
        raise np.linalg.LinAlgError("velocity trust metric has zero numerical rank")
    trust_v = trust_vh[:trust_rank, :].T
    trust_s_kept = trust_s[:trust_rank]
    # z = velocity_scales * velocity_coefficients and
    # z = trust_v @ (trust_coordinates / trust_s_kept).
    av_scaled = av_projected / velocity_scales[None, :]
    a_tilde = (av_scaled @ trust_v) / trust_s_kept[None, :]
    del trust_r, trust_u, trust_vh, normalized_velocity_metric, weighted_values, value_matrix, av_scaled
    gc.collect()
    solver = _solve_multiplier(a_tilde, projected_rhs, trust_bound)
    y = solver["coordinates"]
    z = trust_v @ (y / trust_s_kept)
    velocity_coefficients = z / velocity_scales
    pressure_rhs = pressure_rhs_projection - (pressure_u.T @ av) @ velocity_coefficients
    pressure_coefficients = pressure_vh.T @ (pressure_rhs / pressure_s)
    coefficients = np.zeros(JOINT_COLUMNS, dtype=float)
    coefficients[velocity_indices] = velocity_coefficients
    coefficients[pressure_indices] = pressure_coefficients
    linear_weighted = design @ coefficients
    linear = (linear_weighted / row_weight).reshape(-1, 3)
    objective_weighted = linear_weighted - rhs
    # KKT residuals are evaluated in the same orthogonal coordinates used by
    # the solve.  In the trust coordinates, the reduced gradient is diagonal
    # in the SVD basis and therefore does not require a normal equation.
    trust_spectral_coordinates = solver["vh"] @ y
    trust_spectral_gradient = solver["singular_values"] * (
        solver["singular_values"] * trust_spectral_coordinates - solver["projected_rhs"]
    )
    reduced_kkt_vector = solver["vh"].T @ trust_spectral_gradient + solver["multiplier"] * y
    reduced_kkt_relative = float(
        np.linalg.norm(reduced_kkt_vector)
        / max(np.linalg.norm(solver["projected_rhs"]), 1.0e-300)
    )
    # Dimensionally consistent trust KKT normalization.  In the retained
    # SVD coordinates, ||A^T b|| = ||s * (U^T b)|| and the multiplier term is
    # ||mu*y||, so the denominator has the same units as the gradient norm.
    trust_atb_norm = float(np.linalg.norm(solver["singular_values"] * solver["projected_rhs"]))
    trust_multiplier_y_norm = float(solver["multiplier"] * np.linalg.norm(y))
    reduced_kkt_consistent_relative = float(
        np.linalg.norm(reduced_kkt_vector)
        / max(trust_atb_norm + trust_multiplier_y_norm, 1.0e-300)
    )
    pressure_stationarity_vector = pressure_vh.T @ (
        pressure_s * (pressure_u.T @ objective_weighted)
    )
    pressure_stationarity_relative = float(
        np.linalg.norm(pressure_stationarity_vector)
        / max(np.linalg.norm(objective_weighted), 1.0e-300)
    )
    del a_tilde, av_projected, pressure_u, pressure_s, pressure_vh, pressure_rhs_projection, projected_rhs, av, ap, design
    gc.collect()
    print(json.dumps({"stage": "trust_solution", "pressure_rank": pressure_rank, "trust_rank": trust_rank, "multiplier": solver["multiplier"], "trust_norm": solver["constrained_norm"], "trust_bound": trust_bound}), flush=True)

    raw_velocity, raw_jacobian = previous._joint_velocity_from_coefficients(patch_data, _coefficients_by_patch(coefficients))
    quadratic = np.einsum("nij,nj->ni", raw_jacobian, raw_velocity)
    baseline_metric = _metric(baseline, weights, points)
    linear_metric = _metric(baseline + linear, weights, points)
    raw_velocity_metric = _metric(raw_velocity, weights, points)
    nonlinear_trials = []
    selected = None
    for backtrack in range(13):
        alpha = 0.5 ** backtrack
        trial = _trial(alpha, baseline, linear, quadratic, patch_data, coefficients, weights, points)
        row = {
            "backtrack": backtrack,
            "alpha": alpha,
            "residual": trial["residual_metric"],
            "velocity_change": trial["velocity_metric"],
            "quadratic_remainder": trial["remainder_metric"],
            "trust_region_pass": bool(trial["velocity_metric"]["volume_L2"] <= trust_bound * (1.0 + 1.0e-10)),
            "l2_relative_to_baseline": float((trial["residual_metric"]["volume_L2"] - baseline_metric["volume_L2"]) / max(baseline_metric["volume_L2"], 1.0e-300)),
            "max_relative_to_baseline": float((trial["residual_metric"]["max_norm"] - baseline_metric["max_norm"]) / max(baseline_metric["max_norm"], 1.0e-300)),
        }
        nonlinear_trials.append(row)
        if row["trust_region_pass"] and row["residual"]["volume_L2"] < baseline_metric["volume_L2"] and selected is None:
            selected = trial | {"backtrack": backtrack}
    if selected is None:
        selected = _trial(0.0, baseline, linear, quadratic, patch_data, coefficients, weights, points) | {"backtrack": None}
    selected_coefficients = selected["coefficients"]
    objective_norm = float(np.linalg.norm(objective_weighted))

    design_info = {
        "shape": [3 * len(points), JOINT_COLUMNS],
        "weighted_shape": [3 * len(points), JOINT_COLUMNS],
        "pressure_rank": int(pressure_rank),
        "velocity_metric_rank": int(trust_rank),
        "velocity_metric_singular_min": float(np.min(trust_s_kept)),
        "velocity_metric_singular_max": float(np.max(trust_s_kept)),
        "velocity_column_scale_min": float(np.min(velocity_scales)),
        "velocity_column_scale_max": float(np.max(velocity_scales)),
        "joint_layout": joint_layout,
        "solver": "pressure SVD projection; velocity-metric QR/SVD; scalar Lagrange-multiplier bisection",
        "unconstrained_velocity_change_norm": solver["unconstrained_norm"],
        "constrained_velocity_change_norm": solver["constrained_norm"],
        "reference_velocity_norm": reference_velocity_norm,
        "trust_bound": trust_bound,
        "trust_ratio": float(solver["constrained_norm"] / max(trust_bound, 1.0e-300)),
        "lagrange_multiplier": solver["multiplier"],
        "multiplier_iterations": solver["iterations"],
        "linear_objective_weighted_norm": objective_norm,
        "kkt": {
            "constraint_abs_error": float(solver["constrained_norm"] - trust_bound) if solver["multiplier"] > 0.0 else 0.0,
            "pressure_stationarity_relative": pressure_stationarity_relative,
            "projected_trust_stationarity_relative": reduced_kkt_relative,
            "projected_trust_stationarity_relative_consistent": reduced_kkt_consistent_relative,
            "projected_trust_gradient_norm": float(np.linalg.norm(reduced_kkt_vector)),
            "projected_trust_atb_norm": trust_atb_norm,
            "projected_trust_multiplier_y_norm": trust_multiplier_y_norm,
            "note": "The trust and pressure solves use orthogonal QR/SVD factors; no normal-equation solve is used.",
        },
    }
    old_report = json.loads(PREVIOUS_REPORT_PATH.read_text(encoding="utf-8")) if PREVIOUS_REPORT_PATH.exists() else {}
    report["design"] = design_info
    report["linear_solution"] = {
        "coefficients": coefficients.tolist(),
        "baseline": baseline_metric,
        "linearized_residual": linear_metric,
        "velocity_change": raw_velocity_metric,
        "quadratic_remainder": _metric(quadratic, weights, points),
    }
    report["nonlinear_backtrack"] = nonlinear_trials
    report["selected"] = {
        "feasible": bool(selected["alpha"] > 0.0),
        "reason": "first exact-quadratic trial improving volume L2 inside the one-percent trust region" if selected["alpha"] > 0.0 else "no improving trust-region trial; zero correction retained",
        "alpha": float(selected["alpha"]),
        "backtrack": selected["backtrack"],
        "coefficients": selected_coefficients.tolist(),
        "velocity_coefficients_patch1": selected_coefficients[:COLUMNS_PER_PATCH][:VELOCITY_COLUMNS_PER_PATCH].tolist(),
        "pressure_coefficients_patch1": selected_coefficients[:COLUMNS_PER_PATCH][VELOCITY_COLUMNS_PER_PATCH:].tolist(),
        "velocity_coefficients_patch2": selected_coefficients[COLUMNS_PER_PATCH:][:VELOCITY_COLUMNS_PER_PATCH].tolist(),
        "pressure_coefficients_patch2": selected_coefficients[COLUMNS_PER_PATCH:][VELOCITY_COLUMNS_PER_PATCH:].tolist(),
        "baseline": baseline_metric,
        "selected_residual": selected["residual_metric"],
        "selected_velocity_change": selected["velocity_metric"],
        "selected_quadratic_remainder": selected["remainder_metric"],
        "selected_l2_change": float(selected["residual_metric"]["volume_L2"] - baseline_metric["volume_L2"]),
        "selected_max_change": float(selected["residual_metric"]["max_norm"] - baseline_metric["max_norm"]),
        "peak_improved": bool(selected["residual_metric"]["max_norm"] < baseline_metric["max_norm"]),
    }
    report["comparison"] = {
        "previous_report": PREVIOUS_REPORT_PATH.name,
        "previous_selected": old_report.get("selected", {}),
        "same_grid_baseline": baseline_metric,
        "same_grid_trust_selected": selected["residual_metric"],
    }
    report["status"] = "completed"
    report["elapsed_seconds"] = float(time.perf_counter() - started)
    _save(output_path, report)
    print(json.dumps({"stage": "completed", "baseline_l2": baseline_metric["volume_L2"], "selected_l2": selected["residual_metric"]["volume_L2"], "selected_max": selected["residual_metric"]["max_norm"], "alpha": selected["alpha"], "elapsed_seconds": report["elapsed_seconds"]}), flush=True)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT_PATH)
    args = parser.parse_args()
    run(args.output)
