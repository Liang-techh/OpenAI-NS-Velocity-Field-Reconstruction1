"""Split-panel mode-1 energy screen on the fixed broad-shear patch.

The earlier ``evolved_wave_energy.json`` used an unsplit 16 by 16 tensor
quadrature.  This bounded follow-up evaluates only mode 1 at the initial mean
and the saved constrained ``Delta k = 1e-4`` endpoint.  Physical radial panels
are split at every saved constrained cutoff plus the missing meridional ``0.62``
cutoff through :func:`broad_wave_mean_fit.patch_nodes`.
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

from broad_shear_dynamic_control import load_saved_field  # noqa: E402
from broad_wave_mean_fit import patch_nodes  # noqa: E402
from evolved_wave_energy import (  # noqa: E402
    _dominant_mode1_coefficients,
    _row_summary,
    _stateful_velocity_mean,
    _support,
)
from fourier_patch_energy_budget import screen_mode, spatial_gradient  # noqa: E402
from grouped_joined_field import install_in_field  # noqa: E402
from meridional_state_cache import value_and_pressure_modes  # noqa: E402


GROWTH_PATH = ROOT / "broad_shear_growth.json"
CONSTRAINED_PATH = ROOT / "meridional_constrained_evolution.json"
UNSPLIT_PATH = ROOT / "evolved_wave_energy.json"
OUTPUT_PATH = ROOT / "evolved_wave_energy_split.json"


def _save(report):
    OUTPUT_PATH.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")


def _endpoint(name, dynamic, values, pressures, constrained):
    if name == "initial":
        fit = constrained["initial_fit"]
        state = np.asarray(fit["state"], dtype=float)
        slope = np.asarray(fit["slope"], dtype=float)
        k = float(fit["k"])
    else:
        state = np.asarray(constrained["endpoint_state"], dtype=float)
        slope = np.asarray(constrained["endpoint_fit"]["slope"], dtype=float)
        k = float(constrained["endpoint_k"])
    return _stateful_velocity_mean(dynamic, values, pressures, state, slope, k), k


def run():
    started = time.perf_counter()
    growth = json.loads(GROWTH_PATH.read_text(encoding="utf-8"))
    constrained = json.loads(CONSTRAINED_PATH.read_text(encoding="utf-8"))
    unsplit = json.loads(UNSPLIT_PATH.read_text(encoding="utf-8"))
    center = tuple(float(value) for value in growth["center"])
    widths = tuple(float(value) for value in growth["widths"])
    degree = int(growth["degree"])
    carrier = np.asarray(growth["carriers"]["1"], dtype=float)
    dynamic, source = load_saved_field()
    replacements = install_in_field(dynamic)
    _, values, pressures, _, _ = value_and_pressure_modes(dynamic, source)
    if len(values) != 25 or len(pressures) != 19:
        raise ValueError("Expected 25 velocity and 19 pressure directions")

    saved_breaks = constrained.get("fit_quadrature", {}).get("radial_breaks")
    if saved_breaks is None:
        saved_breaks = constrained.get("constraint_builder", {}).get("radial_breaks")
    if saved_breaks is None:
        raise ValueError("Constrained report has no radial split breaks")
    breaks = sorted(set(float(value) for value in saved_breaks) | {0.62})
    endpoints = {}
    for name in ("initial", "constrained_endpoint"):
        field, k = _endpoint(
            "initial" if name == "initial" else "endpoint",
            dynamic,
            values,
            pressures,
            constrained,
        )
        endpoints[name] = (field, k)

    report = {
        "status": "initializing",
        "accepted": False,
        "pde_validated": False,
        "wave_integrated": False,
        "source_growth_report": GROWTH_PATH.name,
        "source_constrained_report": CONSTRAINED_PATH.name,
        "unsplit_reference_report": UNSPLIT_PATH.name,
        "center": list(center),
        "widths": list(widths),
        "width_factor": 8,
        "degree": degree,
        "mode": 1,
        "carrier": carrier.tolist(),
        "viscosity": float(dynamic.nu),
        "velocity_state_count": len(values),
        "pressure_direction_count": len(pressures),
        "pressure_coefficients_zero": True,
        "grouped_backend_replacements": int(replacements),
        "radial_breaks": breaks,
        "quadrature_definitions": {
            "moderate": {"z_order": 12, "panel_order": 6},
            "higher": {"z_order": 16, "panel_order": 9},
            "weighting": "actual physical cylindrical volume weights from patch_nodes",
            "same_physical_patch": True,
        },
        "matrix_cache": {
            "saved": False,
            "note": "Only spectra and higher-order dominant coefficients were retained; full complex M, G, and whitening matrices were not saved in this bounded run.",
        },
        "rows": {},
        "scope": (
            "Mode-1 diffusion/strain energy screen only. Radial cutoff panels are "
            "split at all saved constrained breaks and 0.62. No wave integration, "
            "nonlinear wave transfer, moment/cone refit, PDE, or acceptance claim."
        ),
    }
    _save(report)
    print(json.dumps({
        "stage": "split_screen_started",
        "break_count": len(breaks),
        "initial_k": endpoints["initial"][1],
        "endpoint_k": endpoints["constrained_endpoint"][1],
        "grouped_backend_replacements": replacements,
    }), flush=True)

    for endpoint_name, (mean, k) in endpoints.items():
        tau = 0.5 * 2.0 ** (-k)
        h = 5.0e-4 * np.sqrt(mean.nu * tau)
        # Check the same physical support including all finite-difference
        # stencil points before calling the evolved mean.
        probe_points, _probe_weights, _probe_geometry = patch_nodes(
            mean, center, widths, tau, 12, 6, breaks
        )
        support = _support(mean.inner, probe_points, tau, h, widths[0])
        endpoint_row = {
            "endpoint": endpoint_name,
            "k": float(k),
            "tau": float(tau),
            "physical_time": float(-tau),
            "support": support,
            "configurations": {},
        }
        if not support["supported"]:
            endpoint_row["status"] = "skipped_unsupported"
            report["rows"][endpoint_name] = endpoint_row
            _save(report)
            continue
        endpoint_row["status"] = "running"
        report["rows"][endpoint_name] = endpoint_row
        _save(report)
        for label, nz, nr in (("moderate", 12, 6), ("higher", 16, 9)):
            points, weights, geometry = patch_nodes(
                mean, center, widths, tau, nz, nr, breaks
            )
            background_gradient = spatial_gradient(mean, points, tau, h)
            mode_row = _row_summary(
                screen_mode(
                    mean,
                    points,
                    weights,
                    center,
                    widths,
                    degree,
                    carrier,
                    1,
                    background_gradient,
                    h,
                )
            )
            mode_row["geometry"] = geometry
            mode_row["weight_sum"] = float(np.sum(weights))
            mode_row["point_count"] = int(len(points))
            if label == "higher":
                mode_row["dominant_coefficients"] = _dominant_mode1_coefficients(
                    mean,
                    points,
                    weights,
                    center,
                    widths,
                    degree,
                    carrier,
                    background_gradient,
                    h,
                )
            endpoint_row["configurations"][label] = mode_row
            _save(report)
            print(json.dumps({
                "stage": "configuration_done",
                "endpoint": endpoint_name,
                "configuration": label,
                "point_count": len(points),
                "combined_lambda": mode_row["combined"]["lambda_max"],
                "energy_rate_2lambda": mode_row["combined"]["energy_rate_2lambda_max"],
                "elapsed": time.perf_counter() - started,
            }), flush=True)
        endpoint_row["status"] = "completed"
        endpoint_row["moderate_to_higher_relative_rate_change"] = float(
            abs(
                endpoint_row["configurations"]["higher"]["combined"]["energy_rate_2lambda_max"]
                - endpoint_row["configurations"]["moderate"]["combined"]["energy_rate_2lambda_max"]
            )
            / max(abs(endpoint_row["configurations"]["higher"]["combined"]["energy_rate_2lambda_max"]), 1.0)
        )
        _save(report)

    # Read the previous unsplit rates for direct comparison; no old case is
    # recomputed in this split-panel run.
    for endpoint_name, unsplit_name in (
        ("initial", "initial"),
        ("constrained_endpoint", "constrained_delta_0.0001"),
    ):
        old_row = unsplit["endpoint_rows"].get(unsplit_name, {})
        row = report["rows"].get(endpoint_name, {})
        if row.get("status") != "completed":
            continue
        row["unsplit_reference"] = {
            "row": unsplit_name,
            "mode1_combined_energy_rate_2lambda": old_row.get("mode1_growth", {}).get("combined_energy_rate_2lambda_max"),
            "mode1_with_divergence_energy_rate_2lambda": old_row.get("mode1_growth", {}).get("combined_with_divergence_energy_rate_2lambda_max"),
        }
    report["status"] = "completed"
    report["elapsed_seconds"] = float(time.perf_counter() - started)
    _save(report)
    print(json.dumps({
        "stage": "completed",
        "elapsed_seconds": report["elapsed_seconds"],
        "initial_higher_rate": report["rows"]["initial"]["configurations"]["higher"]["combined"]["energy_rate_2lambda_max"],
        "constrained_higher_rate": report["rows"]["constrained_endpoint"]["configurations"]["higher"]["combined"]["energy_rate_2lambda_max"],
    }), flush=True)
    return report


if __name__ == "__main__":
    run()
