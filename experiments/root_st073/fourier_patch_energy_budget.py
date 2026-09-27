"""Bounded instantaneous energy-production screen for the degree 2 patch.

The screen projects diffusion and the symmetric background strain onto the
compact exact-curl velocity basis at the initial saved time.  It is a local
homogeneous perturbation diagnostic: mean forcing, nonlinear perturbation
transfer, total energy, and Navier--Stokes acceptance are outside its scope.
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss


ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from fourier_patch_evolution import basis_jets  # noqa: E402
from joined_field import coordinates  # noqa: E402
from outer_feedback_evolution import Trajectory, build_current  # noqa: E402


DEFAULT_INPUT = ROOT / "fourier_patch_implicit.json"
DEFAULT_MEAN = ROOT / "outer_feedback_evolution.json"
DEFAULT_OUTPUT = ROOT / "fourier_patch_energy_budget.json"


def spatial_gradient(field, points, tau, h):
    """Return Cartesian ``dU_i/dx_j`` using spatial fourth-order FD only."""

    points = np.asarray(points, dtype=float)
    gradient = np.empty((len(points), 3, 3), dtype=float)
    for axis in range(3):
        direction = np.zeros(3)
        direction[axis] = h
        minus2 = field.fields(points - 2.0 * direction, tau)[0]
        minus1 = field.fields(points - direction, tau)[0]
        plus1 = field.fields(points + direction, tau)[0]
        plus2 = field.fields(points + 2.0 * direction, tau)[0]
        gradient[:, :, axis] = (minus2 - 8.0 * minus1 + 8.0 * plus1 - plus2) / (12.0 * h)
    return gradient


def support_evidence(inner, points, tau, h, radial_halfwidth):
    """Check axis and registered axial slab before any mean-field call."""

    points = np.asarray(points, dtype=float)
    stencil = [points]
    for axis in range(3):
        direction = np.zeros(3)
        direction[axis] = 2.0 * h
        stencil.extend((points - direction, points + direction))
    stencil = np.concatenate(stencil, axis=0)
    radii = np.hypot(stencil[:, 0], stencil[:, 1])
    scale = np.sqrt(inner.nu)
    coord = coordinates(radii / scale, stencil[:, 2] / scale, tau, inner.h)
    eta = np.asarray(coord["eta"])
    tau_min = 0.5 * 2.0 ** (-inner.p.k_max)
    return {
        "radial_halfwidth": float(radial_halfwidth),
        "minimum_stencil_radius": float(np.min(radii)),
        "eta_abs_max_including_fd_stencil": float(np.max(np.abs(eta))),
        "eta_max_registered": float(inner.p.eta_max),
        "tau": float(tau),
        "tau_min_registered": float(tau_min),
        "radial_off_axis_supported": bool(np.min(radii) > 0.0),
        "axial_slab_supported": bool(np.max(np.abs(eta)) <= inner.p.eta_max),
        "time_slab_supported": bool(tau_min <= tau <= 0.5),
    }


def phase_metrics(center, widths, carriers, modes):
    """Record carrier phase half-widths and angular carrier ranges."""

    r0 = float(center[0])
    dr, dz = widths
    rmin, rmax = r0 - dr, r0 + dr
    result = {}
    for mode in modes:
        kr, kz = carriers[mode]
        result[str(mode)] = {
            "abs_kr_times_radial_halfwidth": float(abs(kr) * dr),
            "abs_kz_times_axial_halfwidth": float(abs(kz) * dz),
            "abs_mode_over_r_min": float(abs(mode) / rmin),
            "abs_mode_over_r_max": float(abs(mode) / rmax),
        }
    return result


def spectrum(matrix):
    """Hermitian eigenvalue extrema and twice-eigenvalue energy rates."""

    matrix = (matrix + matrix.conj().T) * 0.5
    values = np.linalg.eigvalsh(matrix)
    return {
        "lambda_min": float(values[0]),
        "lambda_max": float(values[-1]),
        "energy_rate_2lambda_min": float(2.0 * values[0]),
        "energy_rate_2lambda_max": float(2.0 * values[-1]),
        "positive_growth_direction_count": int(np.sum(values > 0.0)),
        "dimension": int(len(values)),
    }


def screen_mode(mean, points, weights, center, widths, degree, carrier, mode,
                background_gradient, h):
    """Assemble one mode's mass, gradient, strain, and diffusion screens."""

    velocity, gradient, _, _, _ = basis_jets(
        points, center, widths, mode, degree, carrier, mean.nu, h
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
    _, singular_values, right = np.linalg.svd(weighted_velocity, full_matrices=False)
    cutoff = max(float(singular_values[0]), np.finfo(float).tiny) * 1.0e-10
    keep = singular_values > cutoff
    # c = Q y gives c^* M c = y^* y in the retained weighted basis.
    whitening = right[keep].conj().T / singular_values[keep][None, :]
    rank = int(np.sum(keep))
    mass_white = whitening.conj().T @ mass @ whitening
    stiffness_white = whitening.conj().T @ stiffness @ whitening
    strain_white = whitening.conj().T @ strain_matrix @ whitening
    divergence_white = whitening.conj().T @ divergence_matrix @ whitening
    diffusion_operator = -mean.nu * stiffness_white
    strain_operator = -strain_white
    combined_operator = diffusion_operator + strain_operator
    # For a sampled non-solenoidal background, the advective energy identity
    # contributes +1/2 div(U) to the amplitude operator.  Keep both forms so
    # the requested -nu K-S screen remains directly visible.
    divergence_operator = 0.5 * divergence_white
    combined_with_divergence = combined_operator + divergence_operator

    strain_norm = np.linalg.norm(strain, axis=(1, 2))
    div_rms = float(np.sqrt(np.sum(weights * divergence**2)))
    strain_rms = float(np.sqrt(np.sum(weights * strain_norm**2)))
    matrix_strain_norm = float(np.linalg.norm(strain_white))
    matrix_divergence_norm = float(np.linalg.norm(divergence_operator))
    return {
        "rank": rank,
        "full_dimension": int(q),
        "svd_cutoff": float(cutoff),
        "singular_max": float(singular_values[0]),
        "singular_min_retained": float(singular_values[keep][-1]),
        "mass_whitened_identity_error": float(
            np.linalg.norm(mass_white - np.eye(rank))
        ),
        "background_divergence_max_abs": float(np.max(np.abs(divergence))),
        "background_divergence_rms": div_rms,
        "background_strain_norm_max": float(np.max(strain_norm)),
        "background_strain_norm_rms": strain_rms,
        "divergence_to_strain_rms_ratio": float(
            div_rms / max(strain_rms, np.finfo(float).tiny)
        ),
        "divergence_matrix_term_to_strain_matrix_ratio": float(
            matrix_divergence_norm / max(matrix_strain_norm, np.finfo(float).tiny)
        ),
        "diffusion_only": spectrum(diffusion_operator),
        "strain_only": spectrum(strain_operator),
        "combined": spectrum(combined_operator),
        "combined_with_divergence_energy_term": spectrum(combined_with_divergence),
    }


def quadrature(center, widths, order):
    nodes, weights = leggauss(order)
    points = np.array(
        [
            [center[0] + x * widths[0], 0.0, center[1] + z * widths[1]]
            for x in nodes
            for z in nodes
        ]
    )
    point_weights = np.outer(weights, weights).ravel() * points[:, 0]
    point_weights /= np.sum(point_weights)
    return points, point_weights


def run(
    input_path=DEFAULT_INPUT,
    mean_path=DEFAULT_MEAN,
    output_path=DEFAULT_OUTPUT,
    include_order24=True,
):
    report = json.loads(Path(input_path).read_text(encoding="utf-8"))
    if not report.get("state_data"):
        print("state_data absent; no energy screen written")
        return None
    state_data = report["state_data"]
    times = np.asarray(state_data["physical_times"], dtype=float)
    tau = float(-times[0])
    center = tuple(float(v) for v in report["center"])
    base_widths = tuple(float(v) for v in report["widths"])
    degree = int(report["degree"])
    modes = [int(m) for m in report["modes"] if int(m) >= 1]
    carriers = {
        int(m): np.asarray(c, dtype=float) for m, c in report["carriers"].items()
    }

    # Construct the saved mean trajectory once.  No mean jets are evaluated;
    # spatial_gradient below calls only mean.fields at the FD stencil points.
    _, _, current = build_current()
    saved_mean = json.loads(Path(mean_path).read_text(encoding="utf-8"))
    mean = Trajectory(current, saved_mean["nodes"])
    h = 0.0005 * np.sqrt(mean.nu * tau)
    print(json.dumps({"stage": "mean_ready", "tau": tau, "h": h}), flush=True)

    factors = (1.0, 2.0, 4.0, 8.0)
    baseline = {}
    started = time.perf_counter()
    for factor in factors:
        widths = (factor * base_widths[0], factor * base_widths[1])
        points, weights = quadrature(center, widths, 16)
        evidence = support_evidence(inner=current.inner, points=points, tau=tau, h=h,
                                    radial_halfwidth=widths[0])
        evidence["radial_support_condition"] = bool(
            center[0] - widths[0] > 0.0
        )
        evidence["supported"] = bool(
            evidence["radial_support_condition"]
            and evidence["radial_off_axis_supported"]
            and evidence["axial_slab_supported"]
            and evidence["time_slab_supported"]
        )
        row = {
            "factor": factor,
            "widths": list(widths),
            "quadrature_order_per_axis": 16,
            "point_count": int(len(points)),
            "support": evidence,
            "carrier_phase_halfwidths": phase_metrics(center, widths, carriers, [0] + modes),
        }
        if not evidence["supported"]:
            row["skipped_reason"] = "Registered axis or axial/time slab would be exceeded; no mean-field call was made."
            baseline[str(factor)] = row
            print(json.dumps({"stage": "factor_skipped", "factor": factor, "support": evidence}), flush=True)
            continue
        background_gradient = spatial_gradient(mean, points, tau, h)
        mode_rows = {}
        for mode in modes:
            mode_rows[str(mode)] = screen_mode(
                mean, points, weights, center, widths, degree, carriers[mode], mode,
                background_gradient, h
            )
        row["modes"] = mode_rows
        baseline[str(factor)] = row
        max_growth = max(mode_rows[str(m)]["combined"]["lambda_max"] for m in modes)
        print(json.dumps({"stage": "factor_complete", "factor": factor,
                          "max_combined_lambda": max_growth,
                          "elapsed": time.perf_counter() - started}), flush=True)

    supported_factors = [
        factor for factor in factors if baseline[str(factor)].get("support", {}).get("supported")
        and "modes" in baseline[str(factor)]
    ]
    promising_factor = None
    if supported_factors:
        promising_factor = max(
            supported_factors,
            key=lambda factor: max(
                baseline[str(factor)]["modes"][str(mode)]["combined"]["lambda_max"]
                for mode in modes
            ),
        )

    order24 = {}
    order24_factors = []
    if include_order24 and supported_factors:
        order24_factors = [1.0] if promising_factor == 1.0 else [1.0, promising_factor]
    for factor in order24_factors:
        widths = (factor * base_widths[0], factor * base_widths[1])
        points, weights = quadrature(center, widths, 24)
        evidence = support_evidence(current.inner, points, tau, h, widths[0])
        evidence["radial_support_condition"] = bool(
            center[0] - widths[0] > 0.0
        )
        evidence["supported"] = bool(
            evidence["radial_support_condition"]
            and evidence["radial_off_axis_supported"]
            and evidence["axial_slab_supported"]
            and evidence["time_slab_supported"]
        )
        row = {
            "factor": factor,
            "widths": list(widths),
            "quadrature_order_per_axis": 24,
            "point_count": int(len(points)),
            "support": evidence,
            "carrier_phase_halfwidths": phase_metrics(center, widths, carriers, [0] + modes),
        }
        if evidence["supported"]:
            background_gradient = spatial_gradient(mean, points, tau, h)
            row["modes"] = {
                str(mode): screen_mode(
                    mean, points, weights, center, widths, degree, carriers[mode], mode,
                    background_gradient, h
                )
                for mode in modes
            }
        else:
            row["skipped_reason"] = "Registered axis or axial/time slab would be exceeded; no mean-field call was made."
        order24[str(factor)] = row
        print(json.dumps({"stage": "order24_complete", "factor": factor,
                          "elapsed": time.perf_counter() - started}), flush=True)

    result = {
        "accepted": False,
        "pde_validated": False,
        "scope": "Initial-time homogeneous compact-patch energy-production screen only; no mean forcing, nonlinear perturbation transfer, total finite energy, stress matching, or NS acceptance.",
        "source_report": "fourier_patch_implicit.json",
        "mean_report": "outer_feedback_evolution.json",
        "degree": degree,
        "initial_physical_time": float(times[0]),
        "tau": tau,
        "center": list(center),
        "base_widths": list(base_widths),
        "viscosity": float(mean.nu),
        "basis_fd_step": float(h),
        "modes": modes,
        "carriers": {str(mode): carriers[mode].tolist() for mode in sorted(carriers)},
        "matrix_definition": {
            "mass": "M = V^* W V",
            "stiffness": "K = grad(V)^* W grad(V)",
            "strain": "S = V^* W sym(grad(U)) V",
            "diffusion_operator": "-nu K in the M-whitened basis",
            "strain_operator": "-S in the M-whitened basis",
            "combined_operator": "-nu K - S in the M-whitened basis",
            "divergence_correction": "+0.5 V^* W div(U) V, reported separately for sampled non-solenoidal background",
            "energy_rate": "2 lambda",
            "whitening": "SVD of sqrt(W)V, retaining singular values above 1e-10 of the maximum",
        },
        "registered_geometry": {
            "inner_axial_eta_max": float(current.inner.p.eta_max),
            "inner_k_max": float(current.inner.p.k_max),
            "note": "Unsupported widened factors are recorded and skipped before mean.fields calls.",
        },
        "promising_factor_by_max_combined_lambda": promising_factor,
        "quadrature16": baseline,
        "quadrature24": order24,
        "elapsed_seconds": float(time.perf_counter() - started),
    }
    output_path = Path(output_path)
    output_path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"stage": "saved", "output": str(output_path),
                      "promising_factor": promising_factor,
                      "elapsed": result["elapsed_seconds"]}), flush=True)
    return result


if __name__ == "__main__":
    run()
