"""Diagnose growth-matrix versus angular-cache cancellation for two wave states."""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from broad_meridional_constrained import load_saved_field  # noqa: E402
from broad_wave_mean_fit import patch_nodes  # noqa: E402
from fourier_patch_energy_budget import spatial_gradient  # noqa: E402
from fourier_patch_evolution import basis_jets  # noqa: E402
from grouped_joined_field import install_in_field  # noqa: E402


SOURCE_PATH = ROOT / "wave_stress_growth_codesign.json"
DYNAMICS_PATH = ROOT / "wave_dynamics_codesign.json"
CACHE_PATH = ROOT / "wave_momentum_projection.npz"
OUTPUT_PATH = ROOT / "wave_growth_consistency.json"


def _pack_decode(value):
    raw = np.asarray(value, dtype=float)
    return raw[..., 0] + 1j * raw[..., 1]


def _hermitian(value):
    return 0.5 * (value + value.conj().T)


def _metric_from_matrices(mass, stiffness, strain, coefficient, viscosity):
    mass_value = float(np.real(np.vdot(coefficient, mass @ coefficient)))
    stiffness_value = float(np.real(np.vdot(coefficient, stiffness @ coefficient)))
    strain_value = float(np.real(np.vdot(coefficient, strain @ coefficient)))
    # Real physical wave = Re(V c): angular averaging contributes one half.
    energy = 0.25 * mass_value
    production = -0.5 * strain_value
    dissipation = 0.5 * float(viscosity) * stiffness_value
    return {
        "mass_complex_quadratic": mass_value,
        "physical_energy": energy,
        "strain_quadratic": strain_value,
        "shear_production": production,
        "stiffness_complex_quadratic": stiffness_value,
        "viscous_dissipation": dissipation,
        "growth_numerator": production - dissipation,
        "growth_lambda": (production - dissipation) / max(2.0 * energy, 1.0e-300),
        "energy_rate_2lambda": (production - dissipation) / max(energy, 1.0e-300),
        "convention": "E=.5 integral |Re(V c)|^2=.25 c*Mc; P=-integral Re(Vc)^T S Re(Vc); D=nu integral |grad Re(Vc)|^2",
    }


def _assemble_axis(mean, points, weights, center, widths, carrier, viscosity, h,
                   background_gradient=None):
    velocity, gradient, _, _, _ = basis_jets(
        points, center, widths, 1, 2, carrier, viscosity, h
    )
    if background_gradient is None:
        background_gradient = spatial_gradient(
            mean, points, float(mean._consistency_tau), h
        )
    strain = 0.5 * (background_gradient + np.swapaxes(background_gradient, 1, 2))
    q = velocity.shape[-1]
    velocity_flat = velocity.reshape(-1, q)
    gradient_flat = gradient.reshape(-1, q)
    mass = velocity_flat.conj().T @ (
        np.repeat(weights, 3)[:, None] * velocity_flat
    )
    stiffness = gradient_flat.conj().T @ (
        np.repeat(weights, 9)[:, None] * gradient_flat
    )
    strain_matrix = np.einsum(
        "nia,nij,njb,n->ab", velocity.conj(), strain, velocity, weights
    )
    return {
        "mass": _hermitian(mass),
        "stiffness": _hermitian(stiffness),
        "strain": _hermitian(strain_matrix),
        "point_count": int(len(points)),
        "weight_sum": float(np.sum(weights)),
    }


def _cache_metric(cache, coefficient, viscosity):
    weights = np.asarray(cache["weights"], dtype=float)
    x = np.concatenate((coefficient.real, coefficient.imag))
    velocity = np.einsum("nci,i->nc", cache["A"], x)
    gradient = np.einsum("ncdi,i->ncd", cache["B"], x)
    mean_gradient = np.asarray(cache["mean_gradient"], dtype=float)
    strain = 0.5 * (mean_gradient + np.swapaxes(mean_gradient, 1, 2))
    energy = 0.5 * float(np.sum(weights[:, None] * velocity * velocity))
    production = -float(np.einsum("n,ni,nij,nj->", weights, velocity, strain, velocity))
    dissipation = float(viscosity) * float(
        np.sum(weights[:, None, None] * gradient * gradient)
    )
    return {
        "point_count": int(len(weights)),
        "angles": 12,
        "weight_sum": float(np.sum(weights)),
        "physical_energy": energy,
        "shear_production": production,
        "viscous_dissipation": dissipation,
        "growth_numerator": production - dissipation,
        "growth_lambda": (production - dissipation) / max(2.0 * energy, 1.0e-300),
        "energy_rate_2lambda": (production - dissipation) / max(energy, 1.0e-300),
        "convention": "direct real mode on saved 12-angle rings; E=.5 sum w |w|^2",
    }


def _relative_matrix_error(actual, reference):
    return float(np.linalg.norm(actual - reference) / max(np.linalg.norm(reference), 1.0e-300))


def run():
    started = time.perf_counter()
    source = json.loads(SOURCE_PATH.read_text(encoding="utf-8"))
    dynamics = json.loads(DYNAMICS_PATH.read_text(encoding="utf-8"))
    source_problem = source["problem"]
    center = tuple(float(v) for v in source["center"])
    widths = tuple(float(v) for v in source["widths"])
    carrier = np.asarray(source["carrier"], dtype=float)
    viscosity = float(source["viscosity"])
    tau = float(source["initial_tau"])
    h = 5.0e-4 * np.sqrt(viscosity * tau)
    breaks = list(source["radial_breaks"])
    original = _pack_decode(source["selected"]["coefficients_original"])
    dynamic = _pack_decode(dynamics["selected"]["coefficients_original"])
    mean, mean_report = load_saved_field()
    install_in_field(mean)
    mean._consistency_tau = tau

    packed = lambda value: np.asarray(value, dtype=float)[..., 0] + 1j * np.asarray(value, dtype=float)[..., 1]
    source_mass = packed(source_problem["mass_original"])
    source_stiffness = packed(source_problem["stiffness_original"])
    source_strain = packed(source_problem["strain_original"])
    source_matrix_metrics = {
        "original_selected": _metric_from_matrices(
            source_mass, source_stiffness, source_strain, original, viscosity
        ),
        "wave_dynamics_selected": _metric_from_matrices(
            source_mass, source_stiffness, source_strain, dynamic, viscosity
        ),
    }
    report = {
        "status": "initializing",
        "accepted": False,
        "pde_validated": False,
        "wave_integrated": False,
        "source_reports": {
            "growth": SOURCE_PATH.name,
            "dynamics": DYNAMICS_PATH.name,
            "projection": CACHE_PATH.name,
        },
        "inputs": {
            "center": list(center),
            "widths": list(widths),
            "carrier": carrier.tolist(),
            "degree": 2,
            "mode": 1,
            "initial_k": float(source["initial_k"]),
            "initial_tau": tau,
            "viscosity": viscosity,
            "radial_breaks": breaks,
            "mean_source_report": mean_report.get("status", "loaded"),
            "original_coefficients": source["selected"]["coefficients_original"],
            "wave_dynamics_coefficients": dynamics["selected"]["coefficients_original"],
        },
        "source_matrix_metrics": source_matrix_metrics,
        "source_matrix_reference": {
            "mass": source_problem["mass_original"],
            "stiffness": source_problem["stiffness_original"],
            "strain": source_problem["strain_original"],
        },
        "scope": (
            "Numerical consistency audit of the source complex growth matrix versus "
            "the saved full-angle real-wave projection. It compares split physical "
            "quadratures and two fixed coefficient snapshots; no optimization, PDE, "
            "trajectory, or acceptance claim."
        ),
    }
    OUTPUT_PATH.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"stage": "source_matrices_loaded", "cache": CACHE_PATH.name}), flush=True)

    cache = np.load(CACHE_PATH, allow_pickle=False)
    cache_metrics = {
        "original_selected": _cache_metric(cache, original, viscosity),
        "wave_dynamics_selected": _cache_metric(cache, dynamic, viscosity),
    }
    report["angular_cache"] = {
        "arrays": {key: list(cache[key].shape) for key in cache.files},
        "metrics": cache_metrics,
        "source_cache_convention": "wave_momentum_projection._assemble_wave_columns uses A=Re(V),-Im(V), B=Re(grad V),-Im(grad V) over the 12-angle cached rings",
    }
    report["status"] = "angular_cache_compared"
    OUTPUT_PATH.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"stage": "angular_cache_compared",
                      "original_lambda": cache_metrics["original_selected"]["growth_lambda"],
                      "dynamics_lambda": cache_metrics["wave_dynamics_selected"]["growth_lambda"]}), flush=True)

    # Build the cache-order axis grid, source grid, and higher split matrices
    # in one grouped Cartesian-gradient call.  The cache uses z=9/r=9,
    # the source report uses z=16/r=9, and the check uses z=20/r=11.
    grids = []
    for label, z_order, radial_order in (
        ("cache_order_axis_z9_r9", 9, 9),
        ("source_z16_r9", 16, 9),
        ("higher_z20_r11", 20, 11),
    ):
        points, weights, geometry = patch_nodes(
            mean, np.asarray(center), np.asarray(widths), tau,
            z_order, radial_order, breaks,
        )
        grids.append((label, points, weights, geometry))
    all_points = np.concatenate([row[1] for row in grids], axis=0)
    all_grad = spatial_gradient(mean, all_points, tau, h)
    split_results = {}
    offset = 0
    for label, points, weights, geometry in grids:
        count = len(points)
        mean._consistency_tau = tau
        result = _assemble_axis(
            mean, points, weights, center, widths, carrier, viscosity, h,
            background_gradient=all_grad[offset:offset + count],
        )
        result["source_matrix_relative_errors"] = {
            "mass": _relative_matrix_error(result["mass"], source_mass),
            "stiffness": _relative_matrix_error(result["stiffness"], source_stiffness),
            "strain": _relative_matrix_error(result["strain"], source_strain),
        }
        result["metrics"] = {
            "original_selected": _metric_from_matrices(result["mass"], result["stiffness"], result["strain"], original, viscosity),
            "wave_dynamics_selected": _metric_from_matrices(result["mass"], result["stiffness"], result["strain"], dynamic, viscosity),
        }
        result["geometry"] = geometry
        # Large matrices/gradients are intentionally summarized; the source
        # matrix plus these metrics are the reproducible comparison payload.
        result.pop("mass", None)
        result.pop("stiffness", None)
        result.pop("strain", None)
        split_results[label] = result
        offset += count
        report["split_quadrature"] = split_results
        report["status"] = "split_quadrature_in_progress"
        OUTPUT_PATH.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
        print(json.dumps({"stage": "split_grid_done", "label": label,
                          "point_count": count,
                          "original_lambda": result["metrics"]["original_selected"]["growth_lambda"],
                          "dynamics_lambda": result["metrics"]["wave_dynamics_selected"]["growth_lambda"]}), flush=True)

    # Preserve the source matrix metrics separately; split results intentionally
    # retain only summaries to keep this diagnostic compact.
    report["status"] = "completed"
    report["elapsed_seconds"] = float(time.perf_counter() - started)
    OUTPUT_PATH.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"stage": "completed", "elapsed_seconds": report["elapsed_seconds"]}), flush=True)
    return report


if __name__ == "__main__":
    run()
