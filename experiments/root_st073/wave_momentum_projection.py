"""Cached instantaneous momentum projection for the co-designed mode-1 wave.

The immutable ``full_wave_frozen_cache.json`` holdout grid is intentionally
used as the training grid for this next bounded optimizer.  The mean field is
evaluated once with the grouped backend, while the compact mode-1 potential is
represented by its 27 complex coefficients as 54 real variables::

    x = [Re(c_0), ..., Re(c_26), Im(c_0), ..., Im(c_26)].

At the reference time the wave coefficient is spatially fixed, so its
complete Cartesian momentum is the exact finite-difference-compatible
quadratic form

    R(x) = R0 + L x + (B x) (A x).

Here ``A`` is the real velocity basis, ``B`` is the real gradient basis, and
``L`` contains the mean-gradient/wave-velocity cross term, wave-gradient/mean
velocity cross term, and wave viscosity.  This module does not optimize or
claim a PDE, time integration, cone closure, or scale recursion.
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
from broad_meridional_constrained import (  # noqa: E402
    load_saved_field as load_constrained_mean,
)
from full_wave_tangent import _basis_columns  # noqa: E402
from fourier_patch_evolution import basis_jets  # noqa: E402
from grouped_joined_field import install_in_field  # noqa: E402


FROZEN_PATH = ROOT / "full_wave_frozen_cache.json"
MEAN_PATH = ROOT / "broad_meridional_constrained.json"
WAVE_PATH = ROOT / "wave_stress_growth_codesign.json"
NPZ_PATH = ROOT / "wave_momentum_projection.npz"
OUTPUT_PATH = ROOT / "wave_momentum_projection.json"


def _sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _save(report):
    OUTPUT_PATH.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")


def _pack_complex(value):
    value = np.asarray(value, dtype=complex)
    return np.stack((value.real, value.imag), axis=-1).tolist()


def _decode_complex(value):
    value = np.asarray(value, dtype=float)
    if value.shape != (27, 2):
        raise ValueError(f"Expected packed 27 complex coefficients, got {value.shape}")
    return value[:, 0] + 1j * value[:, 1]


def _wave_inputs(frozen, wave_report):
    frozen_wave = frozen["inputs"]["wave"]
    selected = wave_report.get("selected")
    if not selected or "coefficients_original" not in selected:
        raise ValueError("co-designed wave report has no selected coefficients")
    coefficient = _decode_complex(selected["coefficients_original"])
    center = np.asarray(wave_report["center"], dtype=float)
    widths = np.asarray(wave_report["widths"], dtype=float)
    carrier = np.asarray(wave_report["carrier"], dtype=float)
    mode = int(wave_report["mode"])
    degree = int(wave_report["degree"])
    if mode != 1 or degree != 2:
        raise ValueError("This cache requires the degree-2 mode-1 wave")
    if not np.allclose(center, frozen_wave["center"]):
        raise ValueError("wave center differs from immutable frozen cache")
    if not np.allclose(widths, frozen_wave["widths"]):
        raise ValueError("wave widths differ from immutable frozen cache")
    if not np.allclose(carrier, frozen_wave["carrier"]):
        raise ValueError("wave carrier differs from immutable frozen cache")
    frozen_coefficient = _decode_complex(frozen_wave["coefficients_original"])
    coefficient_difference = float(np.max(np.abs(coefficient - frozen_coefficient)))
    return {
        "center": center,
        "widths": widths,
        "carrier": carrier,
        "mode": mode,
        "degree": degree,
        "coefficient": coefficient,
        "frozen_coefficient": frozen_coefficient,
        "coefficient_difference_from_frozen": coefficient_difference,
    }


def _real_complex_columns(value):
    """Convert complex q-columns to the original all-real/all-imag layout."""

    value = np.asarray(value, dtype=complex)
    return np.concatenate((value.real, -value.imag), axis=-1)


def _assemble_wave_columns(points, mean_u, mean_grad, center, widths, carrier,
                           degree, nu, hspace):
    """Build A, B and L for the complete 27-complex-coefficient mode."""

    velocity, wave_grad, wave_diffusion, _, _ = basis_jets(
        points, center, widths, 1, degree, carrier, nu, hspace
    )
    # L_complex[n, component, q] is the coefficient of c_q in the linear
    # cross and viscosity terms.  The mean gradient acts on wave velocity;
    # wave gradient acts on the mean velocity.
    mean_grad_wave = np.einsum("ncd,ndq->ncq", mean_grad, velocity)
    wave_grad_mean = np.einsum("ncdq,nd->ncq", wave_grad, mean_u)
    linear_complex = mean_grad_wave + wave_grad_mean + wave_diffusion
    A = _real_complex_columns(velocity)
    B = _real_complex_columns(wave_grad)
    L = _real_complex_columns(linear_complex)
    return A, B, L, {
        "velocity_complex_shape": list(velocity.shape),
        "gradient_complex_shape": list(wave_grad.shape),
        "diffusion_complex_shape": list(wave_diffusion.shape),
        "real_coefficient_layout": "[Re(c_0..c_26), Im(c_0..c_26)]",
        "gradient_layout": "B[n, velocity_component, spatial_derivative, coefficient]",
        "linear_definition": "Jmean*V + Jwave*Umean - nu*lap(V)",
    }


def _assemble_tangent_design(points, center, widths, carrier, degree, weights):
    """Build and SVD-scale the full mode-0/1/2 180-column tangent matrix."""

    carriers = {0: np.zeros(2), 1: np.asarray(carrier, dtype=float),
                2: 2.0 * np.asarray(carrier, dtype=float)}
    design, layout = _basis_columns(
        points, center, widths, carriers, int(degree)
    )
    row_weight = np.repeat(np.sqrt(np.asarray(weights, dtype=float)), 3)
    weighted = design * row_weight[:, None]
    scales = np.maximum(np.linalg.norm(weighted, axis=0), 1.0e-30)
    normalized = weighted / scales[None, :]
    U, singular, Vh = np.linalg.svd(normalized, full_matrices=False)
    return design, U, singular, Vh, scales, layout


def _array(cache, name):
    if isinstance(cache, dict):
        return np.asarray(cache[name])
    return np.asarray(cache[name])


def residual_and_jacobian(cache, x):
    """Evaluate cached momentum and its analytic 54-column Jacobian.

    Parameters
    ----------
    cache : mapping or ``numpy.lib.npyio.NpzFile``
        Must provide ``R0``, ``A``, ``B`` and ``L`` arrays from the saved NPZ.
    x : array_like, shape (54,)
        All-real/all-imag coefficient vector.

    Returns
    -------
    residual : ndarray, shape (N, 3)
    jacobian : ndarray, shape (N, 3, 54)
    """

    x = np.asarray(x, dtype=float)
    if x.shape != (54,):
        raise ValueError(f"Expected shape (54,), got {x.shape}")
    R0 = _array(cache, "R0")
    A = _array(cache, "A")
    B = _array(cache, "B")
    L = _array(cache, "L")
    wave_velocity = np.einsum("nci,i->nc", A, x)
    wave_gradient = np.einsum("ncdi,i->ncd", B, x)
    linear = np.einsum("nci,i->nc", L, x)
    residual = (
        R0 + linear
        + np.einsum("ncd,nd->nc", wave_gradient, wave_velocity)
    )
    # d[(B x)(A x)]/dx_j = B_j (A x) + (B x) A_j.
    jacobian = (
        L.copy()
        + np.einsum("ncdi,nd->nci", B, wave_velocity)
        + np.einsum("ncd,ndi->nci", wave_gradient, A)
    )
    return residual, jacobian


def cached_energy_norm(cache, residual):
    """Return the physical-volume weighted L2 norm of a residual."""

    residual = np.asarray(residual, dtype=float)
    weights = _array(cache, "weights")
    if residual.shape != (len(weights), 3):
        raise ValueError("residual shape does not match cached points")
    return float(np.sqrt(np.sum(weights[:, None] * residual * residual)))


def weighted_volume_l2(cache, residual):
    """Alias used by downstream projection code."""

    return cached_energy_norm(cache, residual)


def _metric(residual, weights):
    residual = np.asarray(residual, dtype=float)
    weights = np.asarray(weights, dtype=float)
    point_norm = np.linalg.norm(residual, axis=1)
    energy = float(np.sqrt(np.sum(weights[:, None] * residual * residual)))
    return {
        "point_count": int(len(residual)),
        "momentum_max": float(np.max(point_norm)),
        "momentum_volume_L2": energy,
        "momentum_volume_RMS": float(energy / np.sqrt(np.sum(weights))),
        "physical_volume": float(np.sum(weights)),
    }


def _directional_check(cache, x, direction, step, active_indices):
    direction = np.asarray(direction, dtype=float)
    plus, _ = residual_and_jacobian(cache, x + step * direction)
    minus, _ = residual_and_jacobian(cache, x - step * direction)
    _, jacobian = residual_and_jacobian(cache, x)
    finite = (plus - minus) / (2.0 * step)
    analytic = np.einsum("nci,i->nc", jacobian, direction)
    finite_active = finite[active_indices]
    analytic_active = analytic[active_indices]
    difference = finite_active - analytic_active
    scale = max(float(np.linalg.norm(analytic_active)), 1.0e-30)
    return {
        "step": float(step),
        "direction_norm": float(np.linalg.norm(direction)),
        "active_point_count": int(len(active_indices)),
        "active_indices": np.asarray(active_indices, dtype=int).tolist(),
        "finite_difference_max_abs": float(np.max(np.linalg.norm(difference, axis=1))),
        "finite_difference_l2": float(np.linalg.norm(difference)),
        "analytic_l2": float(np.linalg.norm(analytic_active)),
        "relative_l2_error": float(np.linalg.norm(difference) / scale),
    }


def _load_frozen_training():
    frozen = json.loads(FROZEN_PATH.read_text(encoding="utf-8"))
    cache = frozen["frozen_cache"]
    points = np.asarray(cache["holdout_points"], dtype=float)
    weights = np.asarray(cache["holdout_weights"], dtype=float)
    target = np.asarray(cache["holdout_residual"], dtype=float)
    if points.shape != (9072, 3) or weights.shape != (9072,):
        raise ValueError(f"unexpected immutable holdout shapes: {points.shape}, {weights.shape}")
    if target.shape != (9072, 3):
        raise ValueError(f"unexpected immutable residual shape: {target.shape}")
    return frozen, points, weights, target


def run():
    started = time.perf_counter()
    frozen, points, weights, frozen_residual = _load_frozen_training()
    mean, mean_report = load_constrained_mean(MEAN_PATH)
    grouped_replacements = int(install_in_field(mean))
    k = float(mean_report["k"])
    tau = float(mean_report["tau"])
    nu = float(mean_report["nu"])
    hspace = float(frozen["timesteps"]["hspace"])
    htime = float(frozen["timesteps"]["htime"])
    if not np.isclose(tau, frozen["inputs"]["mean"]["tau"]):
        raise ValueError("mean tau differs from immutable frozen cache")
    if not np.isclose(nu, frozen["inputs"]["mean"]["nu"]):
        raise ValueError("mean viscosity differs from immutable frozen cache")
    wave_report = json.loads(WAVE_PATH.read_text(encoding="utf-8"))
    wave = _wave_inputs(frozen, wave_report)
    report = {
        "status": "initializing",
        "accepted": False,
        "pde_validated": False,
        "scale_recursion_established": False,
        "wave_integrated": False,
        "inputs": {
            "frozen_cache": {
                "path": FROZEN_PATH.name,
                "sha256": _sha256(FROZEN_PATH),
                "training_source": "immutable holdout_points/holdout_weights/holdout_residual",
                "holdout_is_reused_as_training": True,
                "point_count": int(len(points)),
                "point_shape": list(points.shape),
                "weight_shape": list(weights.shape),
                "quadrature_source": frozen["inputs"]["geometry"],
            },
            "mean": {
                "path": MEAN_PATH.name,
                "sha256": _sha256(MEAN_PATH),
                "k": k,
                "tau": tau,
                "physical_time": "t=-tau",
                "nu": nu,
                "coefficients": mean_report.get("coefficients", []),
            },
            "wave": {
                "path": WAVE_PATH.name,
                "sha256": _sha256(WAVE_PATH),
                "mode": wave["mode"],
                "degree": wave["degree"],
                "center": wave["center"].tolist(),
                "widths": wave["widths"].tolist(),
                "carrier": wave["carrier"].tolist(),
                "coefficients_original": _pack_complex(wave["coefficient"]),
                "coefficient_layout": "packed [[Re(c_j), Im(c_j)]], j=0..26",
                "coefficient_difference_from_frozen": wave["coefficient_difference_from_frozen"],
            },
        },
        "timesteps": {
            "hspace": hspace,
            "htime": htime,
            "source": "full_wave_frozen_cache.json",
        },
        "grouped_backend_replacements": grouped_replacements,
        "scope": (
            "Complete instantaneous Cartesian momentum projection as a cached "
            "quadratic function of the 27 complex mode-1 potential coefficients. "
            "The frozen order-9 holdout is reused as training; no independent "
            "validation, optimizer, finite-time evolution, PDE, cone, moment, "
            "or scale-recursion acceptance is claimed."
        ),
    }
    _save(report)
    print(json.dumps({"stage": "initializing", "point_count": len(points),
                      "grouped_backend_replacements": grouped_replacements}), flush=True)

    # One mean FD jet evaluation.  R0 includes the complete mean advection;
    # U and J are reused in every wave column below.
    mean_jet = jets(mean, points, tau, hspace, htime)
    mean_u, mean_grad = mean_jet[0], mean_jet[1]
    R0 = momentum(mean_jet)
    A, B, L, wave_layout = _assemble_wave_columns(
        points, mean_u, mean_grad, wave["center"], wave["widths"],
        wave["carrier"], wave["degree"], nu, hspace,
    )
    tangent_design, tangent_U, tangent_s, tangent_Vh, tangent_scales, tangent_layout = _assemble_tangent_design(
        points, wave["center"], wave["widths"], wave["carrier"],
        wave["degree"], weights,
    )
    arrays = {
        "points": points,
        "weights": weights,
        "mean_velocity": mean_u,
        "mean_gradient": mean_grad,
        "R0": R0,
        "A": A,
        "B": B,
        "L": L,
        "tangent_design": tangent_design,
        "tangent_U": tangent_U,
        "tangent_s": tangent_s,
        "tangent_Vh": tangent_Vh,
        "tangent_scales": tangent_scales,
    }
    np.savez_compressed(NPZ_PATH, **arrays)
    report["status"] = "assembled"
    report["arrays"] = {
        key: {"shape": list(value.shape), "dtype": str(value.dtype)}
        for key, value in arrays.items()
    }
    report["mean_snapshot"] = {
        "R0_metric": _metric(R0, weights),
        "mean_velocity_shape": list(mean_u.shape),
        "mean_gradient_shape": list(mean_grad.shape),
        "mean_jet_definition": "affine_momentum.jets(mean, points, tau, hspace, htime)",
    }
    report["wave_snapshot"] = wave_layout
    report["tangent_snapshot"] = {
        "layout": tangent_layout,
        "column_count": int(tangent_design.shape[1]),
        "rank": int(np.sum(tangent_s > 1.0e-10 * tangent_s[0])),
        "singular_values": tangent_s.tolist(),
        "singular_value_threshold": float(1.0e-10 * tangent_s[0]),
        "weighted_normalization": "rows multiplied by sqrt(weights repeated over 3 Cartesian components); columns divided by tangent_scales",
        "design_definition": "full_wave_tangent._basis_columns(mode 0, 1, 2 velocity and pressure-gradient response)",
    }
    report["elapsed_assembly_seconds"] = float(time.perf_counter() - started)
    _save(report)
    print(json.dumps({"stage": "assembled", "A": list(A.shape), "B": list(B.shape),
                      "L": list(L.shape), "tangent": list(tangent_design.shape),
                      "elapsed_seconds": report["elapsed_assembly_seconds"]}), flush=True)

    final_coefficients = wave["coefficient"]
    x = np.concatenate((final_coefficients.real, final_coefficients.imag))
    predicted, jacobian = residual_and_jacobian(arrays, x)
    difference = predicted - frozen_residual
    report["final_coefficient_vector"] = {
        "x": x.tolist(),
        "layout": "[Re(c_0..c_26), Im(c_0..c_26)]",
        "complex_packed": _pack_complex(final_coefficients),
        "l2_norm": float(np.linalg.norm(x)),
    }
    report["frozen_final_replay"] = {
        "predicted_metric": _metric(predicted, weights),
        "immutable_target_metric": _metric(frozen_residual, weights),
        "difference_max_abs": float(np.max(np.linalg.norm(difference, axis=1))),
        "difference_volume_L2": cached_energy_norm(arrays, difference),
        "difference_relative_to_target_L2": float(
            cached_energy_norm(arrays, difference)
            / max(cached_energy_norm(arrays, frozen_residual), np.finfo(float).tiny)
        ),
        "comparison": "quadratic cached reconstruction versus immutable frozen FD residual",
    }
    residual_norm = np.linalg.norm(predicted, axis=1)
    active_count = min(64, len(points))
    active_indices = np.argsort(residual_norm)[-active_count:]
    direction_one = x / max(float(np.linalg.norm(x)), 1.0)
    rng = np.random.default_rng(73073)
    direction_two = rng.normal(size=54)
    direction_two /= max(float(np.linalg.norm(direction_two)), 1.0e-30)
    step = 1.0e-5 * max(float(np.linalg.norm(x)), 1.0)
    report["directional_derivative_checks"] = {
        "active_point_rule": "64 points with largest predicted residual norm",
        "direction_one": _directional_check(
            arrays, x, direction_one, step, active_indices
        ),
        "direction_two": _directional_check(
            arrays, x, direction_two, step, active_indices
        ),
    }
    report["status"] = "completed"
    report["elapsed_seconds"] = float(time.perf_counter() - started)
    report["cache_path"] = NPZ_PATH.name
    _save(report)
    print(json.dumps({
        "stage": "completed",
        "point_count": len(points),
        "predicted_l2": report["frozen_final_replay"]["predicted_metric"]["momentum_volume_L2"],
        "target_l2": report["frozen_final_replay"]["immutable_target_metric"]["momentum_volume_L2"],
        "difference_max": report["frozen_final_replay"]["difference_max_abs"],
        "difference_l2": report["frozen_final_replay"]["difference_volume_L2"],
        "directional_relative_errors": [
            report["directional_derivative_checks"][name]["relative_l2_error"]
            for name in ("direction_one", "direction_two")
        ],
        "elapsed_seconds": report["elapsed_seconds"],
    }), flush=True)
    return report


if __name__ == "__main__":
    run()
