"""Bounded growth and wave-stress co-design for the fixed mode-1 patch.

This diagnostic uses the constrained initial-mean five-node target from
``broad_wave_covariance.json``.  The mode-1 degree-2 potential coefficients
are represented in the mass-whitened exact-curl basis on split physical
panels.  A bounded SLSQP search fits the ten real radial-flux components while
retaining the homogeneous instantaneous energy condition

    y* G y >= 1000 y* y,

where ``G`` is ``-nu K - S`` after mass whitening.  Full problem matrices are
checkpointed before optimization for later joint construction.  This remains
an instantaneous finite-patch diagnostic; it does not integrate a wave or
claim Navier--Stokes acceptance.
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.optimize import minimize


ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from broad_shear_dynamic_control import load_saved_field  # noqa: E402
from broad_wave_mean_fit import patch_nodes  # noqa: E402
from evolved_wave_energy import (  # noqa: E402
    _stateful_velocity_mean,
    _support,
)
from fourier_patch_evolution import basis_jets  # noqa: E402
from fourier_patch_energy_budget import spatial_gradient  # noqa: E402
from grouped_joined_field import install_in_field  # noqa: E402
from meridional_state_cache import value_and_pressure_modes  # noqa: E402
from supported_fourier_basis import basis_layout, basis_data  # noqa: E402


GROWTH_PATH = ROOT / "broad_shear_growth.json"
CONSTRAINED_PATH = ROOT / "meridional_constrained_evolution.json"
TARGET_PATH = ROOT / "broad_wave_covariance.json"
SPLIT_ENERGY_PATH = ROOT / "evolved_wave_energy_split.json"
OUTPUT_PATH = ROOT / "wave_stress_growth_codesign.json"

MODE = 1
ENERGY_THRESHOLD = 1.0e3
GAUGE_THRESHOLD = 1.0e-10
Z_ORDER = 16
PANEL_ORDER = 9
ANGULAR_CHECK_COUNT = 64
COEFFICIENT_BOUND_MULTIPLIER = 10.0


def _save(report):
    OUTPUT_PATH.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")


def _pack_complex(value):
    value = np.asarray(value)
    return np.stack((value.real, value.imag), axis=-1).tolist()


def _decode_complex(value):
    value = np.asarray(value, dtype=float)
    return value[..., 0] + 1j * value[..., 1]


def _hermitian(value):
    return 0.5 * (value + value.conj().T)


def _flux_form(a, b):
    """Hermitian H with c*Hc = .5 Re((a c) conj(b c))."""

    return 0.25 * (
        np.outer(np.conjugate(a), b) + np.outer(np.conjugate(b), a)
    )


def _pack_real_matrix(value):
    return _pack_complex(np.asarray(value, dtype=complex))


def _real_vector(y):
    y = np.asarray(y, dtype=complex)
    return np.concatenate((y.real, y.imag))


def _complex_vector(x, rank):
    x = np.asarray(x, dtype=float)
    return x[:rank] + 1j * x[rank:]


def _quadratic(y, matrix):
    return float(np.real(np.vdot(y, matrix @ y)))


def _quadratic_gradient(y, matrix):
    value = matrix @ y
    return 2.0 * np.concatenate((value.real, value.imag))


def _assemble_problem(mean, points, weights, center, widths, degree, carrier,
                      target_nodes):
    """Assemble physical mass, growth, whitening, and five-node flux forms."""

    tau = float(target_nodes["tau"])
    h = 5.0e-4 * np.sqrt(mean.nu * tau)
    background_gradient = target_nodes["background_gradient"]
    velocity, gradient, _, _, _ = basis_jets(
        points, center, widths, MODE, degree, carrier, mean.nu, h
    )
    q = velocity.shape[-1]
    velocity_flat = velocity.reshape(-1, q)
    gradient_flat = gradient.reshape(-1, q)
    velocity_weight = np.repeat(weights, 3)
    gradient_weight = np.repeat(weights, 9)
    mass = velocity_flat.conj().T @ (velocity_weight[:, None] * velocity_flat)
    stiffness = gradient_flat.conj().T @ (gradient_weight[:, None] * gradient_flat)
    strain = 0.5 * (
        background_gradient + np.swapaxes(background_gradient, 1, 2)
    )
    strain_matrix = np.einsum(
        "nia,nij,njb,n->ab", velocity.conj(), strain, velocity, weights
    )
    growth = _hermitian(-mean.nu * stiffness - strain_matrix)
    weighted_velocity = velocity_flat * np.sqrt(velocity_weight)[:, None]
    _, singular, right = np.linalg.svd(weighted_velocity, full_matrices=False)
    cutoff = max(float(singular[0]), np.finfo(float).tiny) * GAUGE_THRESHOLD
    keep = singular > cutoff
    whitening = right[keep].conj().T / singular[keep][None, :]
    mass_white = _hermitian(whitening.conj().T @ mass @ whitening)
    growth_white = _hermitian(whitening.conj().T @ growth @ whitening)

    flux_original = []
    flux_white = []
    for node in target_nodes["nodes"]:
        node_point = np.asarray(node["point"], dtype=float).reshape(1, 3)
        node_velocity = basis_data(
            node_point, center, widths, MODE, degree, carrier
        )[0][0]
        # All saved physical wave nodes use theta = 0, so Cartesian x/y are
        # cylindrical radial/theta components at the node.
        radial = node_velocity[0]
        theta = node_velocity[1]
        axial = node_velocity[2]
        forms = (_flux_form(radial, theta), _flux_form(radial, axial))
        flux_original.append(forms)
        flux_white.append(tuple(
            _hermitian(whitening.conj().T @ form @ whitening)
            for form in forms
        ))
    return {
        "mass": mass,
        "stiffness": stiffness,
        "strain": strain_matrix,
        "growth": growth,
        "mass_white": mass_white,
        "growth_white": growth_white,
        "whitening": whitening,
        "singular_values": singular,
        "keep": keep,
        "flux_original": flux_original,
        "flux_white": flux_white,
        "rank": int(np.sum(keep)),
        "full_dimension": int(q),
        "gauge_cutoff": float(cutoff),
        "mass_identity_error": float(
            np.linalg.norm(mass_white - np.eye(int(np.sum(keep))))
        ),
        "background_gradient": background_gradient,
        "basis_gradient": gradient,
        "basis_velocity": velocity,
        "h": float(h),
    }


def _flux_values(y, forms):
    return np.asarray([
        [_quadratic(y, forms[node][0]), _quadratic(y, forms[node][1])]
        for node in range(len(forms))
    ])


def _objective_and_gradient(y, forms, target, scales):
    predicted = _flux_values(y, forms)
    residual = (predicted - target) / scales
    value = float(np.sum(residual * residual))
    # SLSQP optimizes the packed real vector [Re y, Im y].
    gradient = np.zeros(2 * len(y), dtype=float)
    for node, forms_node in enumerate(forms):
        for component, matrix in enumerate(forms_node):
            gradient += (
                2.0 * residual[node, component] / scales[node, component]
                * _quadratic_gradient(y, matrix)
            )
    return value, gradient


def _growth_margin(y, growth, threshold):
    return _quadratic(y, growth) - float(threshold) * float(np.vdot(y, y).real)


def _growth_margin_gradient(y, growth, threshold):
    return _quadratic_gradient(y, growth - float(threshold) * np.eye(len(y)))


def _trial_summary(x, problem, target, scales, bound):
    rank = problem["rank"]
    y = _complex_vector(x, rank)
    predicted = _flux_values(y, problem["flux_white"])
    residual = predicted - target
    mass_norm = _quadratic(y, problem["mass_white"])
    growth_value = _quadratic(y, problem["growth_white"])
    return {
        "success": False,
        "coefficient_norm": float(np.linalg.norm(y)),
        "coefficient_bound": float(bound),
        "mass_norm": mass_norm,
        "growth_quadratic": growth_value,
        "growth_lambda": float(growth_value / max(mass_norm, 1.0e-300)),
        "growth_margin": float(growth_value - ENERGY_THRESHOLD * mass_norm),
        "flux_predicted": predicted.tolist(),
        "flux_residual": residual.tolist(),
        "normalized_objective": float(np.sum((residual / scales) ** 2)),
        "relative_flux_l2": float(
            np.linalg.norm(residual) / max(np.linalg.norm(target), 1.0e-300)
        ),
        "coefficients_whitened": _pack_complex(y),
        "coefficients_original": _pack_complex(problem["whitening"] @ y),
    }


def _optimize(problem, target, seed, bound, trial_index):
    rank = problem["rank"]
    scales = np.maximum(np.abs(target), 1.0)
    growth = problem["growth_white"]
    forms = problem["flux_white"]
    eye = np.eye(rank)

    def objective(x):
        y = _complex_vector(x, rank)
        return _objective_and_gradient(y, forms, target, scales)

    def growth_constraint(x):
        y = _complex_vector(x, rank)
        return _growth_margin(y, growth, ENERGY_THRESHOLD)

    def growth_jacobian(x):
        y = _complex_vector(x, rank)
        return _growth_margin_gradient(y, growth, ENERGY_THRESHOLD)

    def norm_constraint(x):
        return float(bound * bound - np.dot(x, x))

    def norm_jacobian(x):
        return -2.0 * np.asarray(x, dtype=float)

    x0 = _real_vector(seed)
    result = minimize(
        lambda x: objective(x),
        x0,
        jac=True,
        method="SLSQP",
        constraints=(
            {"type": "ineq", "fun": growth_constraint, "jac": growth_jacobian},
            {"type": "ineq", "fun": norm_constraint, "jac": norm_jacobian},
        ),
        options={"maxiter": 600, "ftol": 1.0e-11, "disp": False},
    )
    summary = _trial_summary(result.x, problem, target, scales, bound)
    summary.update({
        "trial_index": int(trial_index),
        "success": bool(result.success),
        "status": int(result.status),
        "message": str(result.message),
        "optimizer_objective": float(result.fun),
        "growth_constraint_valid": bool(summary["growth_margin"] >= -1.0e-6),
        "norm_constraint_valid": bool(summary["coefficient_norm"] <= bound + 1.0e-8),
    })
    return summary


def _angular_covariance(coefficient, nodes, center, widths, degree, carrier):
    angles = 2.0 * np.pi * np.arange(ANGULAR_CHECK_COUNT) / ANGULAR_CHECK_COUNT
    rows = []
    for node in nodes:
        radius = float(np.hypot(node["point"][0], node["point"][1]))
        z = float(node["point"][2])
        points = np.column_stack((
            radius * np.cos(angles),
            radius * np.sin(angles),
            np.full(len(angles), z),
        ))
        velocity = basis_data(points, center, widths, MODE, degree, carrier)[0]
        complex_velocity = np.einsum("niq,q->ni", velocity, coefficient)
        cylindrical = np.column_stack((
            np.cos(angles) * complex_velocity[:, 0]
            + np.sin(angles) * complex_velocity[:, 1],
            -np.sin(angles) * complex_velocity[:, 0]
            + np.cos(angles) * complex_velocity[:, 1],
            complex_velocity[:, 2],
        ))
        real_wave = np.real(cylindrical)
        covariance = real_wave.T @ real_wave / float(len(angles))
        rows.append({
            "eta": node["eta"],
            "radial_fraction": node["radial_fraction"],
            "covariance": covariance.tolist(),
            "flux_pair": [float(covariance[0, 1]), float(covariance[0, 2])],
        })
    return rows


def _energy_check(mean, points, weights, center, widths, degree, carrier,
                  coefficient, background_gradient, h):
    velocity, gradient, _, _, _ = basis_jets(
        points, center, widths, MODE, degree, carrier, mean.nu, h
    )
    q = velocity.shape[-1]
    vflat = velocity.reshape(-1, q)
    gflat = gradient.reshape(-1, q)
    vw = np.repeat(weights, 3)
    gw = np.repeat(weights, 9)
    mass = vflat.conj().T @ (vw[:, None] * vflat)
    stiffness = gflat.conj().T @ (gw[:, None] * gflat)
    strain = 0.5 * (background_gradient + np.swapaxes(background_gradient, 1, 2))
    strain_matrix = np.einsum(
        "nia,nij,njb,n->ab", velocity.conj(), strain, velocity, weights
    )
    growth = _hermitian(-mean.nu * stiffness - strain_matrix)
    mass_value = float(np.real(np.vdot(coefficient, mass @ coefficient)))
    growth_value = float(np.real(np.vdot(coefficient, growth @ coefficient)))
    return {
        "point_count": int(len(points)),
        "weight_sum": float(np.sum(weights)),
        "mass_norm": mass_value,
        "growth_quadratic": growth_value,
        "growth_lambda": float(growth_value / max(mass_value, 1.0e-300)),
        "energy_rate_2lambda": float(2.0 * growth_value / max(mass_value, 1.0e-300)),
    }


def _velocity_rms_check(mean, points, weights, tau, center, widths, degree,
                        carrier, coefficient):
    """Report physical real-mode RMS against the local mean RMS.

    The stored complex mode is the amplitude of one nonzero Fourier harmonic;
    angular averaging of its real part contributes one half of the complex
    squared amplitude.
    """
    mode_velocity = basis_data(points, center, widths, MODE, degree, carrier)[0]
    complex_wave = np.einsum("niq,q->ni", mode_velocity, coefficient)
    mean_velocity = mean.fields(points, tau)[0]
    volume = float(np.sum(weights))
    complex_wave_rms = np.sqrt(
        float(np.sum(weights * np.sum(np.abs(complex_wave) ** 2, axis=1)))
        / max(volume, 1.0e-300)
    )
    physical_wave_rms = complex_wave_rms / np.sqrt(2.0)
    mean_rms = np.sqrt(
        float(np.sum(weights * np.sum(np.abs(mean_velocity) ** 2, axis=1)))
        / max(volume, 1.0e-300)
    )
    return {
        "volume": volume,
        "complex_mode_rms": float(complex_wave_rms),
        "physical_real_mode_rms": float(physical_wave_rms),
        "mean_velocity_rms": float(mean_rms),
        "physical_to_mean_rms": float(physical_wave_rms / max(mean_rms, 1.0e-300)),
    }


def run():
    started = time.perf_counter()
    growth_report = json.loads(GROWTH_PATH.read_text(encoding="utf-8"))
    constrained_report = json.loads(CONSTRAINED_PATH.read_text(encoding="utf-8"))
    target_report = json.loads(TARGET_PATH.read_text(encoding="utf-8"))
    center = tuple(float(value) for value in growth_report["center"])
    widths = tuple(float(value) for value in growth_report["widths"])
    degree = int(growth_report["degree"])
    carrier = np.asarray(growth_report["carriers"][str(MODE)], dtype=float)
    dynamic, source = load_saved_field()
    replacements = install_in_field(dynamic)
    _, values, pressures, _, _ = value_and_pressure_modes(dynamic, source)
    if len(values) != 25 or len(pressures) != 19:
        raise ValueError("Expected 25 velocity and 19 pressure directions")
    fit = constrained_report["initial_fit"]
    mean = _stateful_velocity_mean(
        dynamic,
        values,
        pressures,
        np.asarray(fit["state"], dtype=float),
        np.asarray(fit["slope"], dtype=float),
        float(fit["k"]),
    )
    wave_nodes = target_report["wave_nodes"]
    targets = np.asarray(
        target_report["target_comparison"]["constrained_initial_mean"]["targets"],
        dtype=float,
    )
    if targets.shape != (5, 2):
        raise ValueError(f"Expected five constrained two-component targets, got {targets.shape}")
    saved_breaks = constrained_report["fit_quadrature"]["radial_breaks"]
    breaks = sorted(set(float(value) for value in saved_breaks) | {0.62})
    tau = 0.5 * 2.0 ** (-float(fit["k"]))
    h = 5.0e-4 * np.sqrt(mean.nu * tau)
    support_points, support_weights, _ = patch_nodes(
        mean, center, widths, tau, Z_ORDER, PANEL_ORDER, breaks
    )
    support = _support(mean.inner, support_points, tau, h, widths[0])
    if not support["supported"]:
        raise ValueError(f"Fixed patch support is not valid: {support}")
    background_gradient = None
    report = {
        "status": "initializing",
        "accepted": False,
        "pde_validated": False,
        "wave_integrated": False,
        "scale_recursion_established": False,
        "source_target_report": TARGET_PATH.name,
        "source_growth_report": GROWTH_PATH.name,
        "source_constrained_report": CONSTRAINED_PATH.name,
        "split_energy_reference": SPLIT_ENERGY_PATH.name,
        "center": list(center),
        "widths": list(widths),
        "width_factor": 8,
        "degree": degree,
        "mode": MODE,
        "carrier": carrier.tolist(),
        "viscosity": float(mean.nu),
        "initial_k": float(fit["k"]),
        "initial_tau": tau,
        "wave_node_count": len(wave_nodes),
        "target_flux": targets.tolist(),
        "target_definition": target_report["target_sign_convention"],
        "radial_breaks": breaks,
        "quadrature": {"z_order": Z_ORDER, "panel_order": PANEL_ORDER,
                        "point_count": int(len(support_points)),
                        "weighting": "actual physical cylindrical volume weights"},
        "support": support,
        "energy_threshold_lambda": ENERGY_THRESHOLD,
        "coefficient_bound_multiplier": COEFFICIENT_BOUND_MULTIPLIER,
        "grouped_backend_replacements": int(replacements),
        "scope": (
            "Bounded instantaneous mode-1 stress/growth co-design. The target is "
            "the saved constrained initial-mean five-node flux pair. This finite "
            "patch optimization does not integrate waves, refit moments/cones, or "
            "establish PDE, continuum, or Navier-Stokes acceptance."
        ),
    }
    _save(report)
    print(json.dumps({
        "stage": "problem_assembly_started",
        "point_count": len(support_points),
        "target_norm": float(np.linalg.norm(targets)),
        "grouped_backend_replacements": replacements,
    }), flush=True)

    background_gradient = spatial_gradient(mean, support_points, tau, h)
    problem = _assemble_problem(
        mean,
        support_points,
        # patch_nodes returns actual physical volume weights.
        support_weights,
        center,
        widths,
        degree,
        carrier,
        {"tau": tau, "nodes": wave_nodes, "background_gradient": background_gradient},
    )
    eigenvalues, eigenvectors = np.linalg.eigh(problem["growth_white"])
    dominant = eigenvectors[:, -1]
    seed_flux = _flux_values(dominant, problem["flux_white"])
    alpha_numerator = float(np.sum(seed_flux * targets))
    alpha_denominator = float(np.sum(seed_flux * seed_flux))
    alpha = max(0.0, alpha_numerator / max(alpha_denominator, 1.0e-300))
    if alpha <= 0.0:
        alpha = 1.0
    seed = np.sqrt(alpha) * dominant
    coefficient_bound = COEFFICIENT_BOUND_MULTIPLIER * float(np.linalg.norm(seed))
    report["problem"] = {
        "rank": problem["rank"],
        "full_dimension": problem["full_dimension"],
        "gauge_cutoff": problem["gauge_cutoff"],
        "singular_values": problem["singular_values"].tolist(),
        "mass_identity_error": problem["mass_identity_error"],
        "mass_original": _pack_real_matrix(problem["mass"]),
        "growth_original": _pack_real_matrix(problem["growth"]),
        "stiffness_original": _pack_real_matrix(problem["stiffness"]),
        "strain_original": _pack_real_matrix(problem["strain"]),
        "mass_whitened": _pack_real_matrix(problem["mass_white"]),
        "growth_whitened": _pack_real_matrix(problem["growth_white"]),
        "whitening": _pack_real_matrix(problem["whitening"]),
        "growth_eigenvalues": eigenvalues.tolist(),
        "growth_eigenvectors": _pack_real_matrix(eigenvectors),
        "flux_forms_original": [
            [_pack_real_matrix(form) for form in node_forms]
            for node_forms in problem["flux_original"]
        ],
        "flux_forms_whitened": [
            [_pack_real_matrix(form) for form in node_forms]
            for node_forms in problem["flux_white"]
        ],
        "target": targets.tolist(),
        "seed_dominant_lambda": float(eigenvalues[-1]),
        "seed_flux": seed_flux.tolist(),
        "seed_flux_scale_alpha": alpha,
        "seed_coefficients_whitened": _pack_complex(seed),
        "seed_coefficients_original": _pack_complex(problem["whitening"] @ seed),
        "coefficient_bound": coefficient_bound,
        "coefficient_bound_definition": "10 times the mass-whitened norm of the dominant-growth seed scaled by its nonnegative flux least-squares amplitude; bound is a search limit, not an infeasibility proof",
        "basis_layout": basis_layout(degree),
    }
    report["status"] = "problem_assembled"
    _save(report)
    print(json.dumps({
        "stage": "problem_assembled",
        "rank": problem["rank"],
        "mass_identity_error": problem["mass_identity_error"],
        "dominant_lambda": float(eigenvalues[-1]),
        "seed_norm": float(np.linalg.norm(seed)),
        "coefficient_bound": coefficient_bound,
    }), flush=True)

    # Three deterministic starts: dominant direction, its phase rotation, and
    # a small admixture of the next growth eigenvector.
    phase_start = np.exp(1j * np.pi / 4.0) * seed
    second = eigenvectors[:, -2]
    mixed = np.sqrt(alpha) * (np.sqrt(0.96) * dominant + np.sqrt(0.04) * second)
    starts = (seed, phase_start, mixed)
    trials = []
    for index, start in enumerate(starts):
        trial = _optimize(problem, targets, start, coefficient_bound, index)
        trials.append(trial)
        report["trials"] = trials
        report["status"] = "optimization_in_progress"
        _save(report)
        print(json.dumps({
            "stage": "trial_done",
            "trial": index,
            "success": trial["success"],
            "objective": trial["normalized_objective"],
            "growth_lambda": trial["growth_lambda"],
            "growth_margin": trial["growth_margin"],
        }), flush=True)

    feasible = [
        trial for trial in trials
        if trial["success"] and trial["growth_constraint_valid"] and trial["norm_constraint_valid"]
    ]
    selected = min(feasible, key=lambda trial: trial["normalized_objective"]) if feasible else None
    report["selected"] = selected
    report["status"] = "selected" if selected is not None else "optimization_failed_or_no_feasible_trial"
    _save(report)
    if selected is not None:
        coefficient = _decode_complex(selected["coefficients_original"])
        covariance = _angular_covariance(coefficient, wave_nodes, center, widths, degree, carrier)
        high_points, high_weights, high_geometry = patch_nodes(
            mean, center, widths, tau, 20, 12, breaks
        )
        high_gradient = spatial_gradient(mean, high_points, tau, h)
        high_energy = _energy_check(
            mean, high_points, high_weights, center, widths, degree, carrier,
            coefficient, high_gradient, h,
        )
        velocity_metrics = _velocity_rms_check(
            mean, high_points, high_weights, tau, center, widths, degree,
            carrier, coefficient,
        )
        selected["angular_covariance_check"] = covariance
        selected["higher_energy_check"] = {"geometry": high_geometry, **high_energy}
        selected["physical_velocity_metrics"] = velocity_metrics
        selected["target_flux"] = targets.tolist()
        selected["target_flux_mismatch"] = (
            np.asarray(selected["flux_predicted"]) - targets
        ).tolist()
        report["selected"] = selected
        report["independent_checks"] = {
            "split_growth_lambda": float(selected["growth_lambda"]),
            "higher_growth_lambda": float(high_energy["growth_lambda"]),
            "higher_energy_rate_2lambda": float(high_energy["energy_rate_2lambda"]),
            "relative_flux_l2": float(selected["relative_flux_l2"]),
            "physical_to_mean_rms": float(velocity_metrics["physical_to_mean_rms"]),
        }
        report["status"] = "completed"
        _save(report)
        print(json.dumps({
            "stage": "independent_checks_done",
            "growth_lambda_split": selected["growth_lambda"],
            "growth_lambda_higher": high_energy["growth_lambda"],
            "relative_flux_l2": selected["relative_flux_l2"],
        }), flush=True)
    report["elapsed_seconds"] = float(time.perf_counter() - started)
    _save(report)
    print(json.dumps({"stage": "completed", "status": report["status"],
                      "elapsed_seconds": report["elapsed_seconds"]}), flush=True)
    return report


if __name__ == "__main__":
    run()
