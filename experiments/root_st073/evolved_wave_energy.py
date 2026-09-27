"""Instantaneous Fourier-patch energy screen along the evolved mean states.

The broad-shear growth candidate supplied the fixed physical factor-8 patch
and carriers.  This diagnostic evaluates that same patch against the initial
mean and the two saved ``Delta k = 0.01`` refreshed endpoints: the dense
two-step state evolution and the four-step sensitivity run.  Each row uses
the weighted velocity mass SVD and reports both ``-nu K - S`` and the sampled
``+0.5 div(U)`` correction.  No wave amplitude is integrated.

Pressure is deliberately discarded.  The mean is represented by
``StatefulMean(dynamic, values, pressures, state, slope, zero_pressure, k)``;
the dynamic field supplies the saved broad-shear state and the nineteen
pressure directions have zero coefficient, while the velocity directions are
the original 25 basis columns.  The compact pressure primitive has zero
velocity and is therefore omitted from this velocity-only screen.
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
from fourier_patch_energy_budget import (  # noqa: E402
    phase_metrics,
    quadrature,
    screen_mode,
    spatial_gradient,
    support_evidence,
)
from fourier_patch_evolution import basis_jets  # noqa: E402
from grouped_joined_field import install_in_field  # noqa: E402
from meridional_state_cache import (  # noqa: E402
    StatefulMean,
    value_and_pressure_modes,
)
from supported_fourier_basis import basis_layout  # noqa: E402


GROWTH_PATH = ROOT / "broad_shear_growth.json"
DENSE_MEAN_PATH = ROOT / "meridional_state_evolution.json"
REFINED_MEAN_PATH = ROOT / "meridional_state_evolution_refined.json"
CONSTRAINED_MEAN_PATH = ROOT / "meridional_constrained_evolution.json"
OUTPUT_PATH = ROOT / "evolved_wave_energy.json"
QUADRATURE_ORDER = 16
EVOLVED_WINDOW_BREAKS = tuple(sorted({
    0.01, 0.10, 0.12, 0.38, 0.40, 0.58, 0.60, 0.62, 0.72,
    0.88, 0.92, 0.93, 0.99,
}))


def _save(report):
    OUTPUT_PATH.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")


def _add_quadrature_note(report):
    report["quadrature_note"] = {
        "order_per_axis": QUADRATURE_ORDER,
        "point_count": QUADRATURE_ORDER ** 2,
        "weighting": "normalized cylindrical r weight on the fixed physical patch",
        "window_split_status": "unsplit_tensor_gauss",
        "evolved_velocity_window_breaks_in_similarity": list(EVOLVED_WINDOW_BREAKS),
        "unresolved_breaks_in_original_radial_breaks": [0.62],
        "note": "The fixed-patch screen records the missing 0.62 meridional cutoff explicitly; endpoint rates are not asserted as cutoff-converged production estimates.",
    }


def _stateful_velocity_mean(dynamic, values, pressures, state, slope, k):
    """Create a state endpoint with all pressure-direction coefficients zero."""

    state = np.asarray(state, dtype=float)
    slope = np.asarray(slope, dtype=float)
    zero_pressure = np.zeros(len(pressures), dtype=float)
    return StatefulMean(
        dynamic,
        values,
        pressures,
        state,
        slope,
        zero_pressure,
        float(k),
    )


def _row_summary(row):
    """Retain the mode matrix diagnostics needed for comparison."""

    keys = (
        "rank",
        "full_dimension",
        "svd_cutoff",
        "singular_max",
        "singular_min_retained",
        "mass_whitened_identity_error",
        "background_divergence_max_abs",
        "background_divergence_rms",
        "background_strain_norm_max",
        "background_strain_norm_rms",
        "divergence_to_strain_rms_ratio",
        "divergence_matrix_term_to_strain_matrix_ratio",
        "diffusion_only",
        "strain_only",
        "combined",
        "combined_with_divergence_energy_term",
    )
    return {key: row[key] for key in keys}


def _support(inner, points, tau, h, radial_halfwidth):
    evidence = support_evidence(inner, points, tau, h, radial_halfwidth)
    evidence["radial_support_condition"] = bool(radial_halfwidth < 1.0)
    evidence["supported"] = bool(
        evidence["radial_support_condition"]
        and evidence["radial_off_axis_supported"]
        and evidence["axial_slab_supported"]
        and evidence["time_slab_supported"]
    )
    return evidence


def _hermitian(matrix):
    return 0.5 * (matrix + matrix.conj().T)


def _dominant_mode1_coefficients(
    mean, points, weights, center, widths, degree, carrier, background_gradient, h
):
    """Export mass-whitened dominant mode-1 coefficients for wave coupling."""

    velocity, gradient, _, _, _ = basis_jets(
        points, center, widths, 1, degree, carrier, mean.nu, h
    )
    strain = 0.5 * (background_gradient + np.swapaxes(background_gradient, 1, 2))
    divergence = np.trace(background_gradient, axis1=1, axis2=2)
    q = velocity.shape[-1]
    velocity_flat = velocity.reshape(-1, q)
    velocity_weight = np.repeat(weights, 3)
    gradient_flat = gradient.reshape(-1, q)
    gradient_weight = np.repeat(weights, 9)
    mass = velocity_flat.conj().T @ (velocity_weight[:, None] * velocity_flat)
    stiffness = gradient_flat.conj().T @ (gradient_weight[:, None] * gradient_flat)
    strain_matrix = np.einsum(
        "nia,nij,njb,n->ab", velocity.conj(), strain, velocity, weights
    )
    divergence_matrix = velocity_flat.conj().T @ (
        np.repeat(weights * divergence, 3)[:, None] * velocity_flat
    )
    weighted_velocity = velocity_flat * np.sqrt(velocity_weight)[:, None]
    _, singular, right = np.linalg.svd(weighted_velocity, full_matrices=False)
    cutoff = max(float(singular[0]), np.finfo(float).tiny) * 1.0e-10
    keep = singular > cutoff
    whitening = right[keep].conj().T / singular[keep][None, :]
    diffusion = -mean.nu * whitening.conj().T @ stiffness @ whitening
    strain_white = whitening.conj().T @ strain_matrix @ whitening
    divergence_white = 0.5 * whitening.conj().T @ divergence_matrix @ whitening
    combined = _hermitian(diffusion - strain_white)
    combined_with_divergence = _hermitian(combined + divergence_white)
    values, vectors = np.linalg.eigh(combined)
    values_div, vectors_div = np.linalg.eigh(combined_with_divergence)
    coefficient = whitening @ vectors[:, -1]
    coefficient_div = whitening @ vectors_div[:, -1]
    pack = lambda value: np.stack((value.real, value.imag), axis=-1).tolist()
    mass_white = whitening.conj().T @ mass @ whitening
    return {
        "basis_layout": basis_layout(degree),
        "operator": "combined = -nu K - S in the mass-whitened basis",
        "operator_with_divergence": "combined + 0.5 V^* W div(U) V",
        "rank": int(np.sum(keep)),
        "full_dimension": int(q),
        "svd_cutoff": float(cutoff),
        "singular_max": float(singular[0]),
        "singular_min_retained": float(singular[keep][-1]),
        "mass_whitened_identity_error": float(
            np.linalg.norm(mass_white - np.eye(int(np.sum(keep))))
        ),
        "combined_lambda_max": float(values[-1]),
        "combined_energy_rate_2lambda_max": float(2.0 * values[-1]),
        "combined_with_divergence_lambda_max": float(values_div[-1]),
        "combined_with_divergence_energy_rate_2lambda_max": float(2.0 * values_div[-1]),
        "potential_coefficients_combined": pack(coefficient),
        "potential_coefficients_with_divergence": pack(coefficient_div),
        "coefficient_convention": "u = -curl(A); coefficients are radial, theta, axial blocks with a-major then b polynomial order",
    }


def _endpoint_report(name, mean, k, inner, center, widths, degree, carriers, modes):
    tau = 0.5 * 2.0 ** (-float(k))
    h = 5.0e-4 * np.sqrt(mean.nu * tau)
    points, weights = quadrature(center, widths, QUADRATURE_ORDER)
    support = _support(inner, points, tau, h, widths[0])
    row = {
        "name": name,
        "k": float(k),
        "tau": float(tau),
        "physical_time": float(-tau),
        "center": list(center),
        "widths": list(widths),
        "quadrature_order_per_axis": QUADRATURE_ORDER,
        "point_count": int(len(points)),
        "support": support,
        "carrier_phase_halfwidths": phase_metrics(center, widths, carriers, modes),
        "status": "support_checked",
    }
    if not support["supported"]:
        row["status"] = "skipped_unsupported"
        row["skipped_reason"] = "Axis, registered axial slab, or registered time slab check failed before mean fields calls."
        return row

    started = time.perf_counter()
    background_gradient = spatial_gradient(mean, points, tau, h)
    row["background_gradient_seconds"] = float(time.perf_counter() - started)
    mode_rows = {}
    for mode in modes:
        mode_rows[str(mode)] = _row_summary(
            screen_mode(
                mean,
                points,
                weights,
                center,
                widths,
                degree,
                carriers[mode],
                mode,
                background_gradient,
                h,
            )
        )
        if mode == 1:
            row["mode1_dominant_coefficients"] = _dominant_mode1_coefficients(
                mean,
                points,
                weights,
                center,
                widths,
                degree,
                carriers[mode],
                background_gradient,
                h,
            )
    row["modes"] = mode_rows
    row["status"] = "completed"
    row["mode1_growth"] = {
        "combined_lambda_max": mode_rows["1"]["combined"]["lambda_max"],
        "combined_energy_rate_2lambda_max": mode_rows["1"]["combined"]["energy_rate_2lambda_max"],
        "combined_with_divergence_lambda_max": mode_rows["1"]["combined_with_divergence_energy_term"]["lambda_max"],
        "combined_with_divergence_energy_rate_2lambda_max": mode_rows["1"]["combined_with_divergence_energy_term"]["energy_rate_2lambda_max"],
        "positive_direction_count_combined": mode_rows["1"]["combined"]["positive_growth_direction_count"],
        "positive_direction_count_with_divergence": mode_rows["1"]["combined_with_divergence_energy_term"]["positive_growth_direction_count"],
    }
    row["largest_combined_energy_rate_mode"] = max(
        modes,
        key=lambda mode: mode_rows[str(mode)]["combined"]["energy_rate_2lambda_max"],
    )
    row["largest_combined_energy_rate_2lambda"] = float(
        mode_rows[str(row["largest_combined_energy_rate_mode"])]
        ["combined"]["energy_rate_2lambda_max"]
    )
    return row


def _endpoint_from_report(report_path, initial_fit, values, pressures):
    report = json.loads(Path(report_path).read_text(encoding="utf-8"))
    if report_path == DENSE_MEAN_PATH:
        run = report["delta_runs"]["0.01"]
        name = "dense_two_step_delta_0.01"
    else:
        run = report["delta_runs"]["0.01"]
        name = "refined_four_step_delta_0.01"
    state = np.asarray(run["refreshed_state"], dtype=float)
    slope = np.asarray(run["refreshed_slope"], dtype=float)
    k = float(run["endpoint_k"])
    return name, k, state, slope


def _load_context():
    """Load the shared fixed patch and velocity basis for append-only checks."""

    growth = json.loads(GROWTH_PATH.read_text(encoding="utf-8"))
    center = tuple(float(value) for value in growth["center"])
    widths = tuple(float(value) for value in growth["widths"])
    degree = int(growth["degree"])
    modes = [int(mode) for mode in growth["modes"] if int(mode) >= 1]
    carriers = {
        int(mode): np.asarray(value, dtype=float)
        for mode, value in growth["carriers"].items()
        if int(mode) >= 1
    }
    dynamic, source = load_saved_field()
    replacements = install_in_field(dynamic)
    _, values, pressures, value_names, pressure_names = value_and_pressure_modes(
        dynamic, source
    )
    if len(values) != 25 or len(pressures) != 19:
        raise ValueError("Expected 25 velocity and 19 pressure directions")
    return {
        "growth": growth,
        "center": center,
        "widths": widths,
        "degree": degree,
        "modes": modes,
        "carriers": carriers,
        "dynamic": dynamic,
        "source": source,
        "values": values,
        "pressures": pressures,
        "value_names": value_names,
        "pressure_names": pressure_names,
        "replacements": replacements,
    }


def _saved_endpoint(path, values):
    report = json.loads(Path(path).read_text(encoding="utf-8"))
    if path == DENSE_MEAN_PATH:
        run = report["delta_runs"]["0.01"]
        name = "dense_two_step_delta_0.01"
        initial_slope = np.asarray(report["initial_fit"]["slope"], dtype=float)
        k0 = float(report["k0"])
        return report, name, k0, np.zeros(len(values)), initial_slope, run
    if path == REFINED_MEAN_PATH:
        run = report["delta_runs"]["0.01"]
        name = "refined_four_step_delta_0.01"
        return report, name, float(run["endpoint_k"]), np.asarray(run["refreshed_state"], dtype=float), np.asarray(run["refreshed_slope"], dtype=float), run
    if path == CONSTRAINED_MEAN_PATH:
        name = "constrained_delta_0.0001"
        return report, name, float(report["endpoint_k"]), np.asarray(report["endpoint_state"], dtype=float), np.asarray(report["endpoint_fit"]["slope"], dtype=float), report
    raise ValueError(f"Unsupported evolution report: {path}")


def run():
    started = time.perf_counter()
    growth = json.loads(GROWTH_PATH.read_text(encoding="utf-8"))
    center = tuple(float(value) for value in growth["center"])
    widths = tuple(float(value) for value in growth["widths"])
    degree = int(growth["degree"])
    modes = [int(mode) for mode in growth["modes"] if int(mode) >= 1]
    carriers = {
        int(mode): np.asarray(value, dtype=float)
        for mode, value in growth["carriers"].items()
        if int(mode) >= 1
    }
    dynamic, source = load_saved_field()
    replacements = install_in_field(dynamic)
    # value_and_pressure_modes preserves the original 25 basis directions and
    # source knots.  The compact pressure object is used only to enumerate the
    # 19 pressure modes; all of their endpoint coefficients are zero here.
    _, values, pressures, value_names, pressure_names = value_and_pressure_modes(
        dynamic, source
    )
    if len(values) != 25 or len(pressures) != 19:
        raise ValueError("Expected 25 velocity and 19 pressure directions")

    dense = json.loads(DENSE_MEAN_PATH.read_text(encoding="utf-8"))
    state0 = np.zeros(len(values), dtype=float)
    slope0 = np.asarray(dense["initial_fit"]["slope"], dtype=float)
    k0 = float(dense["k0"])
    initial_mean = _stateful_velocity_mean(
        dynamic, values, pressures, state0, slope0, k0
    )
    refined_name, refined_k, refined_state, refined_slope = _endpoint_from_report(
        REFINED_MEAN_PATH, dense["initial_fit"], values, pressures
    )
    dense_name, dense_k, dense_state, dense_slope = _endpoint_from_report(
        DENSE_MEAN_PATH, dense["initial_fit"], values, pressures
    )
    dense_mean = _stateful_velocity_mean(
        dynamic, values, pressures, dense_state, dense_slope, dense_k
    )
    refined_mean = _stateful_velocity_mean(
        dynamic, values, pressures, refined_state, refined_slope, refined_k
    )

    report = {
        "status": "initializing",
        "accepted": False,
        "pde_validated": False,
        "wave_integrated": False,
        "scale_recursion_established": False,
        "source_growth_report": GROWTH_PATH.name,
        "dense_mean_report": DENSE_MEAN_PATH.name,
        "refined_mean_report": REFINED_MEAN_PATH.name,
        "initial_k": k0,
        "initial_tau": float(0.5 * 2.0 ** (-k0)),
        "center": list(center),
        "widths": list(widths),
        "width_factor": 8,
        "degree": degree,
        "modes": modes,
        "carriers": {str(mode): carriers[mode].tolist() for mode in modes},
        "viscosity": float(dynamic.nu),
        "velocity_state_count": len(values),
        "pressure_direction_count": len(pressures),
        "pressure_coefficients_zero": True,
        "grouped_backend_replacements": int(replacements),
        "basis_value_names": value_names,
        "basis_pressure_names": pressure_names,
        "same_physical_patch": True,
        "endpoint_rows": {},
        "scope": (
            "Instantaneous homogeneous Fourier-patch energy screen for the broad-shear "
            "seed along two saved mean-state replays. Mass whitening, diffusion, "
            "strain, and sampled divergence correction are reported. No perturbation "
            "wave time integration, nonlinear transfer, total energy, moment/cone "
            "closure, or Navier-Stokes acceptance is claimed."
        ),
    }
    _add_quadrature_note(report)
    _save(report)
    print(json.dumps({
        "stage": "mean_ready",
        "initial_k": k0,
        "dense_endpoint_k": dense_k,
        "refined_endpoint_k": refined_k,
        "point_patch": 256,
        "grouped_backend_replacements": replacements,
    }), flush=True)

    endpoints = (
        ("initial", initial_mean, k0),
        (dense_name, dense_mean, dense_k),
        (refined_name, refined_mean, refined_k),
    )
    for name, mean, k in endpoints:
        row = _endpoint_report(
            name, mean, k, dynamic.inner, center, widths, degree, carriers, modes
        )
        report["endpoint_rows"][name] = row
        _save(report)
        print(json.dumps({
            "stage": "endpoint_done",
            "name": name,
            "k": k,
            "supported": row["support"]["supported"],
            "mode1_energy_rate_2lambda": row.get("mode1_growth", {}).get(
                "combined_with_divergence_energy_rate_2lambda_max"
            ),
            "elapsed": time.perf_counter() - started,
        }), flush=True)

    initial_row = report["endpoint_rows"]["initial"]
    comparison = {}
    for name in (dense_name, refined_name):
        row = report["endpoint_rows"][name]
        if row.get("status") != "completed":
            continue
        comparison[name] = {
            "mode1_combined_energy_rate_2lambda_initial": initial_row["mode1_growth"]["combined_energy_rate_2lambda_max"],
            "mode1_combined_energy_rate_2lambda_endpoint": row["mode1_growth"]["combined_energy_rate_2lambda_max"],
            "mode1_with_divergence_energy_rate_2lambda_initial": initial_row["mode1_growth"]["combined_with_divergence_energy_rate_2lambda_max"],
            "mode1_with_divergence_energy_rate_2lambda_endpoint": row["mode1_growth"]["combined_with_divergence_energy_rate_2lambda_max"],
            "all_mode_combined_rates_endpoint": {
                mode: row["modes"][str(mode)]["combined"]["energy_rate_2lambda_max"]
                for mode in modes
            },
            "all_mode_with_divergence_rates_endpoint": {
                mode: row["modes"][str(mode)]["combined_with_divergence_energy_term"]["energy_rate_2lambda_max"]
                for mode in modes
            },
            "mass_whitened_identity_error_mode1": row["modes"]["1"]["mass_whitened_identity_error"],
            "rank_mode1": row["modes"]["1"]["rank"],
            "divergence_to_strain_rms_ratio_mode1": row["modes"]["1"]["divergence_to_strain_rms_ratio"],
        }
    report["comparison"] = comparison
    report["status"] = "completed"
    report["elapsed_seconds"] = float(time.perf_counter() - started)
    _save(report)
    print(json.dumps({
        "stage": "completed",
        "elapsed_seconds": report["elapsed_seconds"],
        "mode1_initial_with_divergence_rate": initial_row["mode1_growth"]["combined_with_divergence_energy_rate_2lambda_max"],
    }), flush=True)
    return report


def _add_seed_reference(report, growth):
    """Record the broad-shear lambda versus energy-rate convention explicitly."""

    selected = growth.get("selected_candidate") or {}
    higher = growth.get("higher_order_fixed_amplitude") or {}
    refined_lambda = higher.get("fixed_selected_mode_growth_lambda")
    report["broad_shear_seed_reference"] = {
        "selected_mode": int(selected.get("mode", 1)),
        "coarse_lambda": selected.get("growth_lambda"),
        "coarse_energy_rate_2lambda": selected.get("energy_rate_2lambda"),
        "higher_panel_lambda": refined_lambda,
        "higher_panel_energy_rate_2lambda": (
            None if refined_lambda is None else 2.0 * float(refined_lambda)
        ),
        "note": "The energy screen reports 2 lambda as the amplitude-squared rate; broad_shear_growth records both lambda and 2 lambda.",
    }


def augment_existing():
    """Append dominant mode-1 coefficients without rerunning all mode screens."""

    report = json.loads(OUTPUT_PATH.read_text(encoding="utf-8"))
    context = _load_context()
    dense = json.loads(DENSE_MEAN_PATH.read_text(encoding="utf-8"))
    refined = json.loads(REFINED_MEAN_PATH.read_text(encoding="utf-8"))
    endpoints = (
        (
            "initial",
            float(dense["k0"]),
            np.zeros(len(context["values"])),
            np.asarray(dense["initial_fit"]["slope"], dtype=float),
        ),
        (
            "dense_two_step_delta_0.01",
            float(dense["delta_runs"]["0.01"]["endpoint_k"]),
            np.asarray(dense["delta_runs"]["0.01"]["refreshed_state"], dtype=float),
            np.asarray(dense["delta_runs"]["0.01"]["refreshed_slope"], dtype=float),
        ),
        (
            "refined_four_step_delta_0.01",
            float(refined["delta_runs"]["0.01"]["endpoint_k"]),
            np.asarray(refined["delta_runs"]["0.01"]["refreshed_state"], dtype=float),
            np.asarray(refined["delta_runs"]["0.01"]["refreshed_slope"], dtype=float),
        ),
    )
    for name, k, state, slope in endpoints:
        row = report["endpoint_rows"][name]
        tau = 0.5 * 2.0 ** (-k)
        h = 5.0e-4 * np.sqrt(context["dynamic"].nu * tau)
        points, weights = quadrature(context["center"], context["widths"], QUADRATURE_ORDER)
        support = _support(context["dynamic"].inner, points, tau, h, context["widths"][0])
        if not support["supported"]:
            row["mode1_dominant_coefficients"] = {"status": "unsupported", "support": support}
            continue
        mean = _stateful_velocity_mean(
            context["dynamic"], context["values"], context["pressures"],
            state, slope, k,
        )
        background_gradient = spatial_gradient(mean, points, tau, h)
        row["mode1_dominant_coefficients"] = _dominant_mode1_coefficients(
            mean,
            points,
            weights,
            context["center"],
            context["widths"],
            context["degree"],
            context["carriers"][1],
            background_gradient,
            h,
        )
    _add_seed_reference(report, context["growth"])
    _add_quadrature_note(report)
    report["dominant_coefficients_augmented"] = True
    _save(report)
    print(json.dumps({"stage": "dominant_coefficients_augmented"}), flush=True)
    return report


def append_constrained():
    """Add the saved constrained k0 + 1e-4 endpoint to the fixed patch."""

    report = json.loads(OUTPUT_PATH.read_text(encoding="utf-8"))
    context = _load_context()
    constrained = json.loads(CONSTRAINED_MEAN_PATH.read_text(encoding="utf-8"))
    state = np.asarray(constrained["endpoint_state"], dtype=float)
    slope = np.asarray(constrained["endpoint_fit"]["slope"], dtype=float)
    k = float(constrained["endpoint_k"])
    mean = _stateful_velocity_mean(
        context["dynamic"], context["values"], context["pressures"], state, slope, k
    )
    row = _endpoint_report(
        "constrained_delta_0.0001",
        mean,
        k,
        context["dynamic"].inner,
        context["center"],
        context["widths"],
        context["degree"],
        context["carriers"],
        context["modes"],
    )
    row["source_constrained_report"] = CONSTRAINED_MEAN_PATH.name
    row["constraints_maintained_source"] = bool(constrained.get("constraints_maintained"))
    row["source_initial_replay"] = constrained.get("initial_replay", {})
    row["source_endpoint_replay"] = constrained.get("endpoint_replay", {})
    report["endpoint_rows"]["constrained_delta_0.0001"] = row
    initial = report["endpoint_rows"]["initial"]
    report.setdefault("comparison", {})["constrained_delta_0.0001"] = {
        "initial_reference": "initial",
        "mode1_combined_energy_rate_2lambda_initial": initial["mode1_growth"]["combined_energy_rate_2lambda_max"],
        "mode1_combined_energy_rate_2lambda_endpoint": row.get("mode1_growth", {}).get("combined_energy_rate_2lambda_max"),
        "mode1_with_divergence_energy_rate_2lambda_initial": initial["mode1_growth"]["combined_with_divergence_energy_rate_2lambda_max"],
        "mode1_with_divergence_energy_rate_2lambda_endpoint": row.get("mode1_growth", {}).get("combined_with_divergence_energy_rate_2lambda_max"),
        "all_mode_combined_rates_endpoint": {
            mode: row["modes"][str(mode)]["combined"]["energy_rate_2lambda_max"]
            for mode in context["modes"]
        } if row.get("status") == "completed" else {},
        "constraints_maintained_source": bool(constrained.get("constraints_maintained")),
        "source_endpoint_k": k,
    }
    report["constrained_reference"] = {
        "path": CONSTRAINED_MEAN_PATH.name,
        "delta_k": float(constrained.get("delta_k", 1.0e-4)),
        "constraints_maintained": bool(constrained.get("constraints_maintained")),
        "initial_cone_pass_count": constrained.get("initial_replay", {}).get("cone_pass_count"),
        "endpoint_cone_pass_count": constrained.get("endpoint_replay", {}).get("cone_pass_count"),
        "initial_moment_max_abs": constrained.get("initial_replay", {}).get("moment_max_abs"),
        "endpoint_moment_max_abs": constrained.get("endpoint_replay", {}).get("moment_max_abs"),
        "note": "Reuses the saved constrained state and tangent; this energy screen does not re-solve or certify the wave.",
    }
    _add_seed_reference(report, context["growth"])
    _add_quadrature_note(report)
    report["status"] = "completed"
    _save(report)
    print(json.dumps({
        "stage": "constrained_endpoint_appended",
        "k": k,
        "supported": row["support"]["supported"],
        "mode1_energy_rate_2lambda": row.get("mode1_growth", {}).get("combined_with_divergence_energy_rate_2lambda_max"),
    }), flush=True)
    return report


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--augment",
        action="store_true",
        help="append dominant mode-1 coefficients to an existing energy report",
    )
    parser.add_argument(
        "--append-constrained",
        action="store_true",
        help="append the saved constrained k0 + 1e-4 endpoint",
    )
    args = parser.parse_args()
    if args.augment and args.append_constrained:
        raise SystemExit("Choose at most one append operation")
    if args.augment:
        augment_existing()
    elif args.append_constrained:
        append_constrained()
    else:
        run()
