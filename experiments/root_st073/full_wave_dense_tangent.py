"""Dense spatial-stability check for the locked full-wave tangent fit.

The saved order-9, 12-angle grid from ``full_wave_tangent.json`` becomes the
training set.  Two prescribed SVD cutoffs are compared, then one new order-13
split grid is evaluated.  Only the training-selected fit with the smallest
training residual receives a direct finite-difference field replay on that
new grid.  This is a spatial diagnostic, not a trajectory or PDE result.
"""

from __future__ import annotations

import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from affine_momentum import jets, momentum  # noqa: E402
from broad_meridional_constrained import load_saved_field as load_constrained_mean  # noqa: E402
from broad_wave_mean_fit import norms  # noqa: E402
from grouped_joined_field import install_in_field  # noqa: E402
from wave_residual_harmonics import budget as harmonic_budget  # noqa: E402
from full_wave_tangent import (  # noqa: E402
    LocalPotentialField,
    _basis_columns,
    _decode_complex,
    _grid,
    _unpack_full,
)


SOURCE_PATH = ROOT / "full_wave_tangent.json"
OUTPUT_PATH = ROOT / "full_wave_dense_tangent.json"

NEW_ORDER = 13
NEW_ANGLES = 12
NEW_SHIFT = 0.47
RCOND_VALUES = (1.0e-10, 0.05)


def _save(report):
    OUTPUT_PATH.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")


def _sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _metric(residual, weights):
    result = norms(np.asarray(residual), np.asarray(weights))
    result["point_count"] = int(len(residual))
    return result


def _weighted_fit(design, residual, weights, rcond):
    sqrt_weights = np.sqrt(np.asarray(weights, dtype=float))
    component_weights = np.repeat(sqrt_weights, 3)
    weighted_design = design * component_weights[:, None]
    weighted_rhs = (-np.asarray(residual).reshape(-1)) * component_weights
    scales = np.maximum(np.linalg.norm(weighted_design, axis=0), 1.0e-30)
    normalized = weighted_design / scales[None, :]
    scaled, _, rank, singular = np.linalg.lstsq(
        normalized, weighted_rhs, rcond=float(rcond)
    )
    coefficients = scaled / scales
    return coefficients, {
        "rcond": float(rcond),
        "rank": int(rank),
        "column_count": int(design.shape[1]),
        "singular_values": singular.tolist(),
        "column_scales": scales.tolist(),
    }


def _pack_full(derivatives, pressures, q):
    pieces = [np.asarray(derivatives[0]).real, np.asarray(pressures[0]).real]
    for mode in (1, 2):
        d = np.asarray(derivatives[mode])
        p = np.asarray(pressures[mode])
        pieces.append(np.column_stack((d.real, d.imag)).reshape(-1))
        pieces.append(np.column_stack((p.real, p.imag)).reshape(-1))
    return np.concatenate(pieces)


def _energy_data(points, weights, center, widths, carrier, coefficient, derivative):
    from supported_fourier_basis import basis_data

    velocity = basis_data(points, center, widths, 1, 2, carrier)[0]
    flat = velocity.reshape(-1, velocity.shape[-1])
    matrix = flat.conj().T @ (np.repeat(np.asarray(weights), 3)[:, None] * flat)
    energy = 0.25 * float(np.real(np.vdot(coefficient, matrix @ coefficient)))
    rate = 0.5 * float(np.real(np.vdot(coefficient, matrix @ derivative)))
    return {
        "physical_mode1_energy": energy,
        "physical_mode1_energy_rate": rate,
        "initial_complex_norm": float(np.linalg.norm(coefficient)),
        "derivative_complex_norm": float(np.linalg.norm(derivative)),
    }


def _fit_row(name, coefficients, svd, design, frozen_residual, weights,
             points, center, widths, degree, carriers, initial, q):
    predicted = frozen_residual.reshape(-1) + design @ coefficients
    derivatives, pressures = _unpack_full(coefficients, q)
    packed = _pack_full(derivatives, pressures, q)
    return {
        "name": name,
        "svd": svd,
        "coefficients": coefficients.tolist(),
        "training_predicted": _metric(predicted.reshape(-1, 3), weights),
        "unpack_repack_max_abs_error": float(np.max(np.abs(packed - coefficients))),
        "derivative_norms": [float(np.linalg.norm(value)) for value in derivatives],
        "pressure_norms": [float(np.linalg.norm(value)) for value in pressures],
        "mode1_energy": _energy_data(
            points, weights, center, widths, carriers[1], initial[1], derivatives[1]
        ),
    }


def run():
    started = time.perf_counter()
    source = json.loads(SOURCE_PATH.read_text(encoding="utf-8"))
    cache = source.get("frozen_cache")
    if not cache:
        raise ValueError("full_wave_tangent.json has no frozen cache")
    source_inputs = source["inputs"]
    center = tuple(source_inputs["wave"]["center"])
    widths = tuple(source_inputs["wave"]["widths"])
    degree = int(source_inputs["wave"]["degree"])
    carrier = np.asarray(source_inputs["wave"]["carrier"], dtype=float)
    carriers = {0: np.zeros(2), 1: carrier, 2: 2.0 * carrier}
    q = int(source["basis_q"])
    initial_mode1 = _decode_complex(source_inputs["wave"]["coefficients_original"])
    initial = (np.zeros(3 * q, complex), initial_mode1, np.zeros(3 * q, complex))
    zero_derivatives = tuple(np.zeros(3 * q, complex) for _ in (0, 1, 2))
    zero_pressures = tuple(np.zeros(q, complex) for _ in (0, 1, 2))
    mean, mean_report = load_constrained_mean()
    replacements = install_in_field(mean)
    tau0 = float(source["initial_tau"])
    hspace = float(source["timesteps"]["hspace"])
    htime = float(source["timesteps"]["htime"])
    breaks = list(source_inputs["geometry"]["radial_breaks"])

    report = {
        "status": "initializing",
        "accepted": False,
        "pde_validated": False,
        "wave_integrated": False,
        "constraints_maintained": False,
        "scale_recursion_established": False,
        "scope": (
            "Dense spatial stability correction. The saved order-9, 12-angle "
            "grid is promoted to training; a new order-13, 12-angle split grid "
            "is evaluated. No time integration, moment/cone constraints, PDE, "
            "or scale-recursion acceptance is claimed."
        ),
        "source_report": SOURCE_PATH.name,
        "source_sha256": _sha256(SOURCE_PATH),
        "source_inputs": source_inputs,
        "mean_loader_report": mean_report.get("status", "loaded"),
        "grouped_backend_replacements": int(replacements),
        "center": list(center),
        "widths": list(widths),
        "degree": degree,
        "mode_list": [0, 1, 2],
        "basis_q": q,
        "initial_tau": tau0,
        "viscosity": float(mean.nu),
        "timesteps": {"hspace": hspace, "htime": htime},
        "radial_breaks": breaks,
        "training_definition": "full_wave_tangent.json frozen_cache holdout: order9, 12 angles, physical split weights",
        "new_grid_definition": {
            "order": NEW_ORDER,
            "angles": NEW_ANGLES,
            "angle_shift": NEW_SHIFT,
            "radial_breaks": breaks,
        },
    }
    _save(report)
    print(json.dumps({"stage": "initializing", "grouped_backend_replacements": replacements}), flush=True)

    train_points = np.asarray(cache["holdout_points"], dtype=float)
    train_weights = np.asarray(cache["holdout_weights"], dtype=float)
    train_residual = np.asarray(cache["holdout_residual"], dtype=float)
    train_design, train_layout = _basis_columns(
        train_points, center, widths, carriers, degree
    )
    report["training_grid"] = {
        "point_count": int(len(train_points)),
        "angles": 12,
        "order": 9,
        "weights_sum": float(np.sum(train_weights)),
        "layout": train_layout,
    }
    report["training_harmonic_budget_frozen"] = harmonic_budget(
        train_points, train_weights, train_residual, 12
    )
    report["status"] = "training_design_assembled"
    _save(report)
    print(json.dumps({"stage": "training_design_assembled", "points": len(train_points)}), flush=True)

    fit_rows = {}
    for rcond in RCOND_VALUES:
        label = "rcond_1e-10" if rcond == 1.0e-10 else "rcond_0.05"
        coefficients, svd = _weighted_fit(
            train_design, train_residual, train_weights, rcond
        )
        fit_rows[label] = _fit_row(
            label, coefficients, svd, train_design, train_residual, train_weights,
            train_points, center, widths, degree, carriers, initial, q,
        )
        report["fits"] = fit_rows
        report["status"] = "training_fits_in_progress"
        _save(report)
        print(json.dumps({"stage": "training_fit_done", "name": label,
                          "rank": svd["rank"],
                          "training_l2": fit_rows[label]["training_predicted"]["volume_L2"]}), flush=True)

    selected_name = min(
        fit_rows,
        key=lambda name: fit_rows[name]["training_predicted"]["volume_L2"],
    )
    report["training_selected"] = {
        "name": selected_name,
        "criterion": "smallest predicted physical training volume L2",
    }
    report["status"] = "training_fits_complete"
    _save(report)
    print(json.dumps({"stage": "training_fits_complete", "selected": selected_name}), flush=True)

    print(json.dumps({"stage": "new_grid_frozen_started", "order": NEW_ORDER,
                      "angles": NEW_ANGLES}), flush=True)
    new_points, new_weights, new_geometry = _grid(
        mean, center, widths, tau0, breaks, NEW_ORDER, NEW_ANGLES, NEW_SHIFT
    )
    frozen = LocalPotentialField(
        mean, center, widths, degree, carriers, initial,
        zero_derivatives, zero_pressures, tau0,
    )
    new_frozen_residual = momentum(jets(frozen, new_points, tau0, hspace, htime))
    report["new_grid"] = new_geometry
    report["new_grid"]["weights_sum"] = float(np.sum(new_weights))
    report["new_frozen_cache"] = {
        "points": new_points.tolist(),
        "weights": new_weights.tolist(),
        "residual": new_frozen_residual.tolist(),
        "weighting": "physical cylindrical volume weights repeated over angles",
    }
    report["new_harmonic_budget_frozen"] = harmonic_budget(
        new_points, new_weights, new_frozen_residual, NEW_ANGLES
    )
    report["status"] = "new_grid_frozen_complete"
    _save(report)
    print(json.dumps({"stage": "new_grid_frozen_complete", "points": len(new_points),
                      "frozen_l2": report["new_harmonic_budget_frozen"]["momentum_volume_L2"]}), flush=True)

    new_design, _ = _basis_columns(new_points, center, widths, carriers, degree)
    new_predicted = {}
    for name, row in fit_rows.items():
        coefficients = np.asarray(row["coefficients"], dtype=float)
        predicted = new_frozen_residual.reshape(-1) + new_design @ coefficients
        new_predicted[name] = {
            "metric": _metric(predicted.reshape(-1, 3), new_weights),
            "harmonic_budget": harmonic_budget(
                new_points, new_weights, predicted.reshape(-1, 3), NEW_ANGLES
            ),
        }
    report["new_grid_predicted"] = new_predicted
    report["status"] = "new_grid_predicted"
    _save(report)
    print(json.dumps({"stage": "new_grid_predicted", "selected": selected_name,
                      "selected_predicted_l2": new_predicted[selected_name]["metric"]["volume_L2"]}), flush=True)

    selected_coefficients = np.asarray(fit_rows[selected_name]["coefficients"], dtype=float)
    selected_derivatives, selected_pressures = _unpack_full(selected_coefficients, q)
    selected_field = LocalPotentialField(
        mean, center, widths, degree, carriers, initial,
        selected_derivatives, selected_pressures, tau0,
    )
    actual_selected = momentum(jets(selected_field, new_points, tau0, hspace, htime))
    report["new_grid_actual_selected"] = {
        "fit_name": selected_name,
        "metric": _metric(actual_selected, new_weights),
        "harmonic_budget": harmonic_budget(
            new_points, new_weights, actual_selected, NEW_ANGLES
        ),
        "difference_from_predicted": _metric(
            actual_selected - (
                new_frozen_residual.reshape(-1)
                + new_design @ selected_coefficients
            ).reshape(-1, 3),
            new_weights,
        ),
    }
    report["status"] = "completed"
    report["elapsed_seconds"] = float(time.perf_counter() - started)
    _save(report)
    print(json.dumps({
        "stage": "completed",
        "status": report["status"],
        "selected": selected_name,
        "new_predicted_l2": report["new_grid_predicted"][selected_name]["metric"]["volume_L2"],
        "new_actual_l2": report["new_grid_actual_selected"]["metric"]["volume_L2"],
        "elapsed_seconds": report["elapsed_seconds"],
    }), flush=True)
    return report


if __name__ == "__main__":
    run()
