"""One bounded reference-velocity correction for the scale-generator defect.

The previous diagnostic replaced only ``u_t`` and found a large momentum
defect.  This module changes the reference velocity itself by one compact,
exact-curl patch step.  It uses analytic degree-2 basis jets for modes 0, 1,
and 2 on the same two local supports as the drift fit.  The frozen pressure
projection is the starting pressure; its 90 pressure-gradient degrees of
freedom remain available in this joint 360-column linearized step.

The trust region is a weighted velocity norm: the proposed reference-velocity
change is limited to one percent of the current reference velocity norm.  The
accepted trial is checked with the exact quadratic convection remainder
``(delta J) @ delta u``.  This is a spatial diagnostic only; support,
moments, cones, shape, time integration, and PDE admission are not checked.
"""

from __future__ import annotations

import argparse
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

from localized_residual_enrichment import _mode_block_columns  # noqa: E402
from scale_generator_linearization import velocity_correction_response  # noqa: E402
from scale_generator_pressure_projection import _pressure_patch_design  # noqa: E402
from supported_fourier_analytic_jets import basis_jets  # noqa: E402


GENERATOR_REPORT_PATH = ROOT / "scale_generator_momentum_defect.json"
GENERATOR_CACHE_PATH = ROOT / "scale_generator_momentum_defect.npz"
PRESSURE_REPORT_PATH = ROOT / "scale_generator_pressure_projection.json"
PRESSURE_SOURCE_PATH = ROOT / "scale_generator_pressure_projection.py"
DRIFT_REPORT_PATH = ROOT / "localized_drift_1800_fit.json"
FROZEN_PATH = ROOT / "full_wave_frozen_cache.json"
LINEARIZATION_PATH = ROOT / "scale_generator_linearization.py"
BASIS_PATH = ROOT / "supported_fourier_analytic_jets.py"
OUTPUT_PATH = ROOT / "scale_reference_velocity_step.json"

DEGREE = 2
MODES = (0, 1, 2)
VELOCITY_COLUMNS_PER_PATCH = 135
PRESSURE_COLUMNS_PER_PATCH = 45
COLUMNS_PER_PATCH = VELOCITY_COLUMNS_PER_PATCH + PRESSURE_COLUMNS_PER_PATCH
RIDGE_RELATIVE = 1.0e-10
VELOCITY_TRUST_FRACTION = 1.0e-2


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


def _load_inputs():
    generator_report = json.loads(GENERATOR_REPORT_PATH.read_text(encoding="utf-8"))
    pressure_report = json.loads(PRESSURE_REPORT_PATH.read_text(encoding="utf-8"))
    drift_report = json.loads(DRIFT_REPORT_PATH.read_text(encoding="utf-8"))
    frozen_report = json.loads(FROZEN_PATH.read_text(encoding="utf-8"))
    if generator_report.get("status") != "completed":
        raise ValueError("Generator defect report is not completed")
    if pressure_report.get("status") != "completed":
        raise ValueError("Pressure projection is not completed")
    if drift_report.get("status") != "completed" or not drift_report.get("selected", {}).get("feasible", False):
        raise ValueError("Latest localized drift field is not completed and feasible")
    with np.load(GENERATOR_CACHE_PATH, allow_pickle=False) as loaded:
        required = {
            "points", "weights", "velocity", "gradient", "laplacian", "generator_residual",
        }
        missing = sorted(required.difference(loaded.files))
        if missing:
            raise ValueError(f"Generator cache is missing {missing}")
        points = np.asarray(loaded["points"], dtype=float)
        weights = np.asarray(loaded["weights"], dtype=float)
        velocity = np.asarray(loaded["velocity"], dtype=float)
        jacobian = np.asarray(loaded["gradient"], dtype=float)
        laplacian = np.asarray(loaded["laplacian"], dtype=float)
        generator_residual = np.asarray(loaded["generator_residual"], dtype=float)
    count = len(points)
    if points.shape != (count, 3) or weights.shape != (count,):
        raise ValueError("Invalid generator cache grid shapes")
    if velocity.shape != (count, 3) or jacobian.shape != (count, 3, 3) or laplacian.shape != (count, 3):
        raise ValueError("Invalid generator cache jet shapes")
    if generator_residual.shape != (count, 3):
        raise ValueError("Invalid generator residual shape")
    if not np.all(np.isfinite(points)) or not np.all(np.isfinite(weights)) or np.any(weights <= 0.0):
        raise ValueError("Generator grid is invalid")
    if any(not np.all(np.isfinite(value)) for value in (velocity, jacobian, laplacian, generator_residual)):
        raise ValueError("Generator cache contains nonfinite arrays")
    wave = frozen_report["inputs"]["wave"]
    carrier = np.asarray(wave["carrier"], dtype=float)
    tau0 = float(generator_report["inputs"]["tau0"])
    h = float(generator_report["inputs"]["similarity_exponent_h"])
    nu = float(generator_report["inputs"]["viscosity"])
    return {
        "generator_report": generator_report,
        "pressure_report": pressure_report,
        "drift_report": drift_report,
        "frozen_report": frozen_report,
        "points": points,
        "weights": weights,
        "velocity": velocity,
        "jacobian": jacobian,
        "laplacian": laplacian,
        "generator_residual": generator_residual,
        "carrier": carrier,
        "tau0": tau0,
        "h": h,
        "nu": nu,
    }


def _pressure_baseline(data):
    points = data["points"]
    drift = data["drift_report"]
    carrier = data["carrier"]
    inputs = drift["inputs"]
    first_center = np.asarray(inputs["first_patch_center"], dtype=float)
    first_widths = np.asarray(inputs["first_patch_widths"], dtype=float)
    second_center = np.asarray(inputs["second_patch_center"], dtype=float)
    second_widths = np.asarray(inputs["second_patch_widths"], dtype=float)
    first_design = _pressure_patch_design(points, first_center, first_widths, carrier)
    second_design = _pressure_patch_design(points, second_center, second_widths, carrier)
    design = np.column_stack((first_design, second_design))
    coefficients = np.asarray(data["pressure_report"]["fit"]["coefficients"], dtype=float)
    if design.shape != (3 * len(points), 90) or coefficients.shape != (90,):
        raise ValueError(f"Unexpected pressure design/coefficient shapes {design.shape} {coefficients.shape}")
    correction = (design @ coefficients).reshape(-1, 3)
    baseline = data["generator_residual"] + correction
    return baseline, correction, design, {
        "first_center": first_center,
        "first_widths": first_widths,
        "second_center": second_center,
        "second_widths": second_widths,
        "coefficients": coefficients,
    }


def _column_component(array, component):
    array = np.asarray(array)
    return array.real if component == "real" else -array.imag


def _mode_velocity_columns(points, center, widths, carrier, nu, h, base_velocity, base_jacobian, base_laplacian, tau0):
    """Build analytic linearized momentum columns and replay metadata."""

    q = (DEGREE + 1) ** 2
    velocity_columns = []
    velocity_values = []
    velocity_jacobians = []
    velocity_laplacians = []
    mean_swirl_columns = []
    layout = []
    theta = np.arctan2(points[:, 1], points[:, 0])
    radius = np.linalg.norm(points[:, :2], axis=1)
    for mode in MODES:
        mode_carrier = np.zeros(2, dtype=float) if mode == 0 else (carrier if mode == 1 else 2.0 * carrier)
        values, gradients, viscous, _, _ = basis_jets(
            points, center, widths, mode, DEGREE, mode_carrier, nu
        )
        laplacian = -viscous / nu
        components = ("real",) if mode == 0 else ("real", "imag")
        for index in range(3 * q):
            for component in components:
                value = _column_component(values[:, :, index], component)
                jacobian = _column_component(gradients[:, :, :, index], component)
                lap = _column_component(laplacian[:, :, index], component)
                if mode == 0:
                    mean_swirl = -np.sin(theta) * value[:, 0] + np.cos(theta) * value[:, 1]
                    mean_swirl = np.where(radius > 0.0, mean_swirl, 0.0)
                else:
                    # Fixed five-angle means of m=1 and m=2 real/imaginary
                    # columns vanish exactly, so they have no swirl correction.
                    mean_swirl = np.zeros(len(points), dtype=float)
                response, _ = velocity_correction_response(
                    points,
                    base_velocity,
                    base_jacobian,
                    value,
                    jacobian,
                    lap,
                    tau0,
                    h,
                    mean_swirl,
                    nu,
                )
                velocity_columns.append(response.reshape(-1))
                velocity_values.append(value)
                velocity_jacobians.append(jacobian)
                velocity_laplacians.append(lap)
                mean_swirl_columns.append(mean_swirl)
                layout.append({"kind": "velocity", "mode": mode, "component": component, "index": index})
    if len(velocity_columns) != VELOCITY_COLUMNS_PER_PATCH:
        raise ValueError(f"Unexpected velocity column count {len(velocity_columns)}")
    return (
        np.column_stack(velocity_columns),
        np.stack(velocity_values, axis=2),
        np.stack(velocity_jacobians, axis=3),
        np.stack(velocity_laplacians, axis=2),
        np.stack(mean_swirl_columns, axis=1),
        layout,
    )


def _velocity_from_coefficients(values, coefficients):
    coefficients = np.asarray(coefficients, dtype=float)
    return np.einsum("niq,q->ni", values, coefficients)


def _joint_velocity_from_coefficients(patch_data, coefficients):
    value_parts = []
    jacobian_parts = []
    for data, control in zip(patch_data, coefficients):
        value_parts.append(_velocity_from_coefficients(data[1], control))
        jacobian_parts.append(np.einsum("nijq,q->nij", data[2], control))
    return sum(value_parts), sum(jacobian_parts)


def _candidate_metrics(alpha, baseline, linear, quadratic, values, jacobians, weights, points):
    velocity_change, jacobian_change = _joint_velocity_from_coefficients(values, alpha)
    residual = baseline + alpha[0] * linear + (alpha[0] ** 2) * quadratic
    return residual, velocity_change, jacobian_change


def run(output_path=OUTPUT_PATH):
    started = time.perf_counter()
    output_path = Path(output_path)
    data = _load_inputs()
    points = data["points"]
    weights = data["weights"]
    base_velocity = data["velocity"]
    base_jacobian = data["jacobian"]
    base_laplacian = data["laplacian"]
    baseline_generator = data["generator_residual"]
    baseline, pressure_correction, pressure_design, pressure_meta = _pressure_baseline(data)
    drift_inputs = data["drift_report"]["inputs"]
    first_center = pressure_meta["first_center"]
    first_widths = pressure_meta["first_widths"]
    second_center = pressure_meta["second_center"]
    second_widths = pressure_meta["second_widths"]
    pressure_before_metric = _metric(baseline_generator, weights, points)
    pressure_after_metric = _metric(baseline, weights, points)
    pressure_report_after = data["pressure_report"]["after"]
    pressure_replay_error = {
        "volume_L2_abs": float(pressure_after_metric["volume_L2"] - pressure_report_after["volume_L2"]),
        "max_abs": float(pressure_after_metric["max_norm"] - pressure_report_after["max_norm"]),
    }
    report = {
        "status": "initializing",
        "accepted": False,
        "pde_validated": False,
        "scale_recursion_established": False,
        "constraints_maintained": False,
        "scope": (
            "One bounded analytic reference-velocity correction for the frozen scale-generator "
            "momentum defect. Two compact degree-2 exact-curl patches use modes 0, 1, and 2; "
            "the existing 90 pressure-gradient columns remain in the joint fit. The one-percent "
            "velocity trust region and exact quadratic convection remainder are checked, while "
            "support, moments, cones, shape, trajectory, PDE, and recursive-scale acceptance are not."
        ),
        "sources": {
            "generator_report": {"path": GENERATOR_REPORT_PATH.name, "sha256": _sha256(GENERATOR_REPORT_PATH)},
            "generator_cache": {"path": GENERATOR_CACHE_PATH.name, "sha256": _sha256(GENERATOR_CACHE_PATH)},
            "pressure_report": {"path": PRESSURE_REPORT_PATH.name, "sha256": _sha256(PRESSURE_REPORT_PATH)},
            "pressure_source": {"path": PRESSURE_SOURCE_PATH.name, "sha256": _sha256(PRESSURE_SOURCE_PATH)},
            "drift_report": {"path": DRIFT_REPORT_PATH.name, "sha256": _sha256(DRIFT_REPORT_PATH)},
            "frozen_geometry": {"path": FROZEN_PATH.name, "sha256": _sha256(FROZEN_PATH)},
            "linearization": {"path": LINEARIZATION_PATH.name, "sha256": _sha256(LINEARIZATION_PATH)},
            "analytic_basis_jets": {"path": BASIS_PATH.name, "sha256": _sha256(BASIS_PATH)},
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
            "patch_count": 2,
            "velocity_columns_per_patch": VELOCITY_COLUMNS_PER_PATCH,
            "pressure_columns_per_patch": PRESSURE_COLUMNS_PER_PATCH,
            "joint_column_count": 2 * COLUMNS_PER_PATCH,
            "velocity_trust_fraction": VELOCITY_TRUST_FRACTION,
            "ridge_relative_cutoff": RIDGE_RELATIVE,
            "carrier": data["carrier"].tolist(),
            "correction_kind": "initial reference curl-potential velocity; time-tangent coefficients are not reused as velocity increments",
        },
        "pressure_baseline": {
            "before": pressure_before_metric,
            "after": pressure_after_metric,
            "reported_after": pressure_report_after,
            "replay_error": pressure_replay_error,
            "correction_metric": _metric(pressure_correction, weights, points),
            "pressure_columns": 90,
            "pressure_coefficients_source": PRESSURE_REPORT_PATH.name,
        },
        "patches": {
            "first": {"center": first_center.tolist(), "widths": first_widths.tolist()},
            "second": {"center": second_center.tolist(), "widths": second_widths.tolist()},
            "joint_order": "patch 1 [135 velocity, 45 pressure], then patch 2 [135 velocity, 45 pressure]",
        },
    }
    _save(output_path, report)
    print(json.dumps({"stage": "sources_loaded", "point_count": len(points), "pressure_replay_error": pressure_replay_error}), flush=True)

    patch_data = []
    joint_columns = []
    joint_layout = []
    for patch_index, (center, widths) in enumerate(((first_center, first_widths), (second_center, second_widths)), 1):
        velocity_response, values, jacobians, laplacians, mean_swirl, velocity_layout = _mode_velocity_columns(
            points,
            center,
            widths,
            data["carrier"],
            data["nu"],
            data["h"],
            base_velocity,
            base_jacobian,
            base_laplacian,
            data["tau0"],
        )
        pressure_start = (patch_index - 1) * PRESSURE_COLUMNS_PER_PATCH
        pressure_end = pressure_start + PRESSURE_COLUMNS_PER_PATCH
        pressure_response = pressure_design[:, pressure_start:pressure_end]
        joint_columns.extend((velocity_response, pressure_response))
        for entry in velocity_layout:
            joint_layout.append({"patch": patch_index, **entry})
        for entry in data["pressure_report"]["pressure_layout"]["per_patch"]:
            joint_layout.append({"patch": patch_index, **entry})
        patch_data.append((velocity_response, values, jacobians, laplacians, mean_swirl))
        print(json.dumps({"stage": "patch_basis", "patch": patch_index, "velocity_columns": velocity_response.shape[1], "pressure_columns": pressure_response.shape[1]}), flush=True)

    design = np.column_stack(joint_columns)
    if design.shape != (3 * len(points), 2 * COLUMNS_PER_PATCH):
        raise ValueError(f"Unexpected joint response shape {design.shape}")
    row_weight = np.repeat(np.sqrt(weights), 3)
    weighted_design = design * row_weight[:, None]
    weighted_rhs = -baseline.reshape(-1) * row_weight
    column_scales = np.maximum(np.linalg.norm(weighted_design, axis=0), 1.0e-30)
    normalized_design = weighted_design / column_scales[None, :]
    singular_values = np.linalg.svd(normalized_design, compute_uv=False)
    cutoff = RIDGE_RELATIVE * max(float(singular_values[0]), 1.0e-300)
    rank = int(np.sum(singular_values > cutoff))
    normalized_solution, _, lstsq_rank, _ = np.linalg.lstsq(
        normalized_design, weighted_rhs, rcond=RIDGE_RELATIVE
    )
    raw_coefficients = normalized_solution / column_scales
    linear_raw = (normalized_design @ normalized_solution / row_weight).reshape(-1, 3)
    raw_velocity, raw_jacobian = _joint_velocity_from_coefficients(
        patch_data,
        (raw_coefficients[:COLUMNS_PER_PATCH][0:VELOCITY_COLUMNS_PER_PATCH],
         raw_coefficients[COLUMNS_PER_PATCH:COLUMNS_PER_PATCH + VELOCITY_COLUMNS_PER_PATCH][0:VELOCITY_COLUMNS_PER_PATCH]),
    )
    # Recompute the exact nonlinear convection remainder for the raw trial.
    quadratic_raw = np.einsum("nij,nj->ni", raw_jacobian, raw_velocity)
    reference_velocity_norm = float(np.sqrt(np.sum(weights * np.sum(base_velocity * base_velocity, axis=1))))
    raw_velocity_norm = float(np.sqrt(np.sum(weights * np.sum(raw_velocity * raw_velocity, axis=1))))
    trust_bound = VELOCITY_TRUST_FRACTION * reference_velocity_norm
    alpha_cap = min(1.0, trust_bound / max(raw_velocity_norm, 1.0e-300))
    if alpha_cap < 1.0:
        alpha_cap *= 0.999999
    baseline_metric = _metric(baseline, weights, points)
    trial_rows = []
    selected = None
    for backtrack in range(12):
        alpha = alpha_cap * (0.5 ** backtrack)
        trial_residual = baseline + alpha * linear_raw + (alpha ** 2) * quadratic_raw
        trial_velocity = alpha * raw_velocity
        trial_jacobian = alpha * raw_jacobian
        trial_metric = _metric(trial_residual, weights, points)
        velocity_metric = _metric(trial_velocity, weights, points)
        remainder = (alpha ** 2) * quadratic_raw
        row = {
            "backtrack": backtrack,
            "alpha": float(alpha),
            "velocity_change": velocity_metric,
            "quadratic_remainder": _metric(remainder, weights, points),
            "residual": trial_metric,
            "l2_relative_to_baseline": float((trial_metric["volume_L2"] - baseline_metric["volume_L2"]) / max(baseline_metric["volume_L2"], 1.0e-300)),
            "trust_region_pass": bool(velocity_metric["volume_L2"] <= trust_bound * (1.0 + 1.0e-12)),
        }
        trial_rows.append(row)
        if selected is None and row["trust_region_pass"] and trial_metric["volume_L2"] < baseline_metric["volume_L2"]:
            selected = {
                "alpha": alpha,
                "residual": trial_residual,
                "velocity": trial_velocity,
                "jacobian": trial_jacobian,
                "remainder": remainder,
                "metric": trial_metric,
                "velocity_metric": velocity_metric,
                "backtrack": backtrack,
            }
    if selected is None:
        selected = {
            "alpha": 0.0,
            "residual": baseline.copy(),
            "velocity": np.zeros_like(base_velocity),
            "jacobian": np.zeros_like(base_jacobian),
            "remainder": np.zeros_like(base_velocity),
            "metric": baseline_metric,
            "velocity_metric": _metric(np.zeros_like(base_velocity), weights, points),
            "backtrack": None,
        }
    selected_coefficients = raw_coefficients * float(selected["alpha"])
    report["design"] = {
        "shape": list(design.shape),
        "weighted_shape": list(weighted_design.shape),
        "normalized_rank": rank,
        "lstsq_rank": int(lstsq_rank),
        "normalized_singular_values": singular_values.tolist(),
        "column_scale_min": float(np.min(column_scales)),
        "column_scale_max": float(np.max(column_scales)),
        "raw_coefficient_norm": float(np.linalg.norm(raw_coefficients)),
        "raw_velocity_norm": raw_velocity_norm,
        "reference_velocity_norm": reference_velocity_norm,
        "trust_bound": trust_bound,
        "alpha_cap": float(alpha_cap),
        "joint_layout": joint_layout,
    }
    report["trust_region_trials"] = trial_rows
    report["selected"] = {
        "feasible": bool(selected["alpha"] > 0.0),
        "reason": "first backtracked exact-quadratic trial improving L2 inside one-percent velocity trust region" if selected["alpha"] > 0.0 else "no improving trust-region trial; zero correction retained",
        "alpha": float(selected["alpha"]),
        "backtrack": selected["backtrack"],
        "coefficients": selected_coefficients.tolist(),
        "velocity_coefficients_patch1": selected_coefficients[:COLUMNS_PER_PATCH][:VELOCITY_COLUMNS_PER_PATCH].tolist(),
        "pressure_coefficients_patch1": selected_coefficients[:COLUMNS_PER_PATCH][VELOCITY_COLUMNS_PER_PATCH:].tolist(),
        "velocity_coefficients_patch2": selected_coefficients[COLUMNS_PER_PATCH:][:VELOCITY_COLUMNS_PER_PATCH].tolist(),
        "pressure_coefficients_patch2": selected_coefficients[COLUMNS_PER_PATCH:][VELOCITY_COLUMNS_PER_PATCH:].tolist(),
        "baseline": baseline_metric,
        "selected_residual": selected["metric"],
        "selected_velocity_change": selected["velocity_metric"],
        "selected_quadratic_remainder": _metric(selected["remainder"], weights, points),
        "selected_l2_change": float(selected["metric"]["volume_L2"] - baseline_metric["volume_L2"]),
        "selected_max_change": float(selected["metric"]["max_norm"] - baseline_metric["max_norm"]),
    }
    report["status"] = "completed"
    report["elapsed_seconds"] = float(time.perf_counter() - started)
    _save(output_path, report)
    print(json.dumps({
        "stage": "completed",
        "baseline_l2": baseline_metric["volume_L2"],
        "selected_l2": selected["metric"]["volume_L2"],
        "selected_max": selected["metric"]["max_norm"],
        "alpha": selected["alpha"],
        "elapsed_seconds": report["elapsed_seconds"],
    }), flush=True)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT_PATH)
    args = parser.parse_args()
    run(args.output)
