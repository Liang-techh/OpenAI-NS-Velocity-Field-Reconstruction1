"""Finite-matrix mean-swirl co-design screen for the supported Fourier patch.

The construction searches exact axisymmetric swirl increments in a cached
factor-8 degree-2 patch.  It optimizes an instantaneous homogeneous energy
matrix for one oscillatory Fourier mode; it does not solve the mean moment
constraints, evolve the perturbation, or establish a Navier--Stokes field.
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import minimize


ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from fourier_patch_evolution import basis_jets  # noqa: E402
from midplane_resolved_feasibility import ZeroBackground  # noqa: E402
from outer_feedback_evolution import SwirlValue, Trajectory, build_current  # noqa: E402
from supported_fourier_basis import basis_data  # noqa: E402


INPUT_PATH = ROOT / "fourier_patch_implicit.json"
MEAN_PATH = ROOT / "outer_feedback_evolution.json"
OUTPUT_PATH = ROOT / "fourier_shear_codesign.json"
TARGET_GROWTH = 1.0e4
OPTIMIZATION_TARGET = TARGET_GROWTH * 1.0001
WIDTH_FACTOR = 8.0


def spatial_gradient(field, points, tau, h):
    """Fourth-order Cartesian spatial gradient; no time or Laplacian calls."""

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


def hermitian(matrix):
    return 0.5 * (matrix + matrix.conj().T)


def pack_complex(value):
    value = np.asarray(value)
    return np.stack((value.real, value.imag), axis=-1).tolist()


def spectrum(matrix):
    values, vectors = np.linalg.eigh(hermitian(matrix))
    return values, vectors


def assemble_mode(mean, points, weights, center, widths, degree, carrier, mode,
                  baseline_gradient, unit_gradients, nu, h):
    """Return whitened H0 and linear swirl matrices for one Fourier mode."""

    velocity, gradient, _, _, _ = basis_jets(
        points, center, widths, mode, degree, carrier, nu, h
    )
    strain = 0.5 * (baseline_gradient + np.swapaxes(baseline_gradient, 1, 2))
    q = velocity.shape[-1]
    velocity_flat = velocity.reshape(-1, q)
    velocity_weight = np.repeat(weights, 3)
    gradient_flat = gradient.reshape(-1, q)
    gradient_weight = np.repeat(weights, 9)
    stiffness = gradient_flat.conj().T @ (gradient_weight[:, None] * gradient_flat)
    strain_matrix = np.einsum(
        "nia,nij,njb,n->ab", velocity.conj(), strain, velocity, weights
    )
    unit_strain = []
    for unit_gradient in unit_gradients:
        unit_strain_field = 0.5 * (
            unit_gradient + np.swapaxes(unit_gradient, 1, 2)
        )
        unit_strain.append(
            np.einsum(
                "nia,nij,njb,n->ab",
                velocity.conj(), unit_strain_field, velocity, weights,
            )
        )

    weighted_velocity = velocity_flat * np.sqrt(velocity_weight)[:, None]
    _, singular_values, right = np.linalg.svd(weighted_velocity, full_matrices=False)
    cutoff = max(float(singular_values[0]), np.finfo(float).tiny) * 1.0e-10
    keep = singular_values > cutoff
    whitening = right[keep].conj().T / singular_values[keep][None, :]
    diffusion = -nu * (whitening.conj().T @ stiffness @ whitening)
    base_strain = -(whitening.conj().T @ strain_matrix @ whitening)
    increment_strain = np.array([
        -(whitening.conj().T @ item @ whitening) for item in unit_strain
    ])
    H0 = hermitian(diffusion + base_strain)
    H = np.array([hermitian(item) for item in increment_strain])
    return {
        "H0": H0,
        "H_increments": H,
        "whitening": whitening,
        "rank": int(np.sum(keep)),
        "dimension": int(q),
        "singular_max": float(singular_values[0]),
        "singular_min_retained": float(singular_values[keep][-1]),
        "diffusion_only": spectrum(diffusion)[0],
        "strain_only": spectrum(base_strain)[0],
        "combined_zero": spectrum(H0)[0],
    }


def mode_candidate(matrix_data, coefficient, coefficient_scale=None):
    H0 = matrix_data["H0"]
    H = matrix_data["H_increments"]
    if coefficient_scale is None:
        coefficient_scale = np.ones(len(H))
    actual = np.asarray(coefficient) * coefficient_scale
    total = H0 + np.einsum("j,jab->ab", actual, H)
    values, vectors = spectrum(total)
    return float(values[-1]), vectors[:, -1], total, actual


def node_points(inner, tau):
    locations = list(dict.fromkeys(
        [(eta, y) for eta in (-0.2, 0.0, 0.2) for y in (0.5, 0.75)]
        + [(eta + de, 0.75 + dy) for eta in (-0.2, 0.2)
           for de in (-0.005, 0.0, 0.005)
           for dy in (-0.005, 0.0, 0.005)]
    ))
    points = []
    labels = []
    for eta, y in locations:
            X = inner.p.X_max * (1.0 + 15.0 * y) ** 2
            points.append(inner.from_similarity([X], [eta], tau)[0])
            labels.append({"eta": eta, "y": y})
    return np.asarray(points), labels


def cone_lambda2(points, velocity, gradient):
    """Return the sampled lambda-squared cone diagnostic at theta = 0."""

    values = []
    for point, value, jacobian in zip(points, velocity, gradient):
        radius = float(np.hypot(point[0], point[1]))
        safe = max(radius, np.finfo(float).tiny)
        # Node points have theta zero; retain the general cylindrical
        # conversion for clarity and to make the cached diagnostic reusable.
        ca, sa = point[0] / safe, point[1] / safe
        u_theta = -sa * value[0] + ca * value[1]
        j_theta_r = -sa * jacobian[0, 0] + ca * jacobian[1, 0]
        F = u_theta / safe
        shear = np.array([j_theta_r - F, jacobian[2, 0]], dtype=float)
        shear_norm = float(np.linalg.norm(shear))
        if shear_norm == 0.0:
            values.append(float("nan"))
            continue
        normal = shear / shear_norm
        values.append(float(-2.0 * F * normal[0] * (2.0 * F * normal[0] + shear_norm)))
    return np.asarray(values, dtype=float)


def candidate_node_data(node_cache, coefficient):
    velocity = node_cache["baseline_velocity"].copy()
    gradient = node_cache["baseline_gradient"].copy()
    velocity += np.einsum("ncj,j->nc", node_cache["unit_velocity"], coefficient)
    gradient += np.einsum("ncij,j->nci", node_cache["unit_gradient"], coefficient)
    return cone_lambda2(node_cache["points"], velocity, gradient)


def covariance_direction(center, widths, degree, carrier, mode, coefficients):
    angles = np.arange(64, dtype=float) * 2.0 * np.pi / 64.0
    points = np.column_stack(
        (center[0] * np.cos(angles), center[0] * np.sin(angles),
         np.full(len(angles), center[1]))
    )
    velocity, _, _ = basis_data(
        points, center=center, widths=widths, mode=mode, degree=degree,
        carrier=carrier
    )
    cylindrical_cart = np.einsum("niq,q->ni", velocity, coefficients).real
    ca, sa = np.cos(angles), np.sin(angles)
    cylindrical = np.column_stack(
        (ca * cylindrical_cart[:, 0] + sa * cylindrical_cart[:, 1],
         -sa * cylindrical_cart[:, 0] + ca * cylindrical_cart[:, 1],
         cylindrical_cart[:, 2])
    )
    covariance = np.mean(cylindrical[:, 0, None] * cylindrical[:, 1:], axis=0)
    return covariance


def run(input_path=INPUT_PATH, mean_path=MEAN_PATH, output_path=OUTPUT_PATH,
        include_order24=True):
    started = time.perf_counter()
    report = json.loads(Path(input_path).read_text(encoding="utf-8"))
    if not report.get("state_data"):
        print("state_data absent; no co-design report written")
        return None
    data = report["state_data"]
    tau = float(-data["physical_times"][0])
    center = tuple(float(v) for v in report["center"])
    base_widths = tuple(float(v) for v in report["widths"])
    widths = tuple(WIDTH_FACTOR * v for v in base_widths)
    degree = int(report["degree"])
    modes = [int(m) for m in report["modes"] if int(m) >= 1]
    carriers = {int(m): np.asarray(c, dtype=float) for m, c in report["carriers"].items()}
    _, base, current = build_current()
    saved_mean = json.loads(Path(mean_path).read_text(encoding="utf-8"))
    mean = Trajectory(current, saved_mean["nodes"])
    h = 0.0005 * np.sqrt(mean.nu * tau)
    points, weights = quadrature(center, widths, 16)
    if center[0] - widths[0] <= 0.0:
        raise ValueError("factor-8 radial support reaches the cylindrical axis")
    print(json.dumps({"stage": "mean_ready", "tau": tau, "point_count": len(points)}), flush=True)

    baseline_velocity, _ = mean.fields(points, tau)
    baseline_gradient = spatial_gradient(mean, points, tau, h)
    zero = ZeroBackground(base)
    unit_fields = [SwirlValue(zero, np.eye(9)[j]) for j in range(9)]
    unit_velocity = np.stack([field.fields(points, tau)[0] for field in unit_fields], axis=2)
    unit_gradient = np.stack([spatial_gradient(field, points, tau, h) for field in unit_fields], axis=3)
    print(json.dumps({"stage": "swirl_cache_ready", "units": 9}), flush=True)

    # The solve_control locations are cached once for optional normalized
    # constraints and always reported as a post-optimization geometry check.
    node_pts, node_labels = node_points(current.inner, tau)
    node_cache = {
        "points": node_pts,
        "labels": node_labels,
        "baseline_velocity": mean.fields(node_pts, tau)[0],
        "baseline_gradient": spatial_gradient(mean, node_pts, tau, h),
        "unit_velocity": np.stack([field.fields(node_pts, tau)[0] for field in unit_fields], axis=2),
        "unit_gradient": np.stack([spatial_gradient(field, node_pts, tau, h) for field in unit_fields], axis=3),
    }
    baseline_node_lambda2 = cone_lambda2(
        node_pts, node_cache["baseline_velocity"], node_cache["baseline_gradient"]
    )
    print(json.dumps({"stage": "node_cache_ready", "baseline_lambda2": baseline_node_lambda2.tolist()}), flush=True)

    gram = np.einsum("nci,ncj,n->ij", unit_velocity, unit_velocity, weights)
    gram = 0.5 * (gram + gram.T)
    gram_eigenvalues = np.linalg.eigvalsh(gram)
    gram_regularization = max(float(np.max(gram_eigenvalues)), 1.0) * 1.0e-10
    gram_regularized = gram + gram_regularization * np.eye(9)
    coefficient_scale = 1.0 / np.sqrt(np.maximum(np.diag(gram_regularized), np.finfo(float).tiny))
    gram_scaled = coefficient_scale[:, None] * gram_regularized * coefficient_scale[None, :]

    matrix_data = {}
    for mode in modes:
        # Assemble once per mode while reusing the already-built mean gradient
        # and cached swirl gradients.
        matrix_data[str(mode)] = assemble_mode(
            mean, points, weights, center, widths, degree, carriers[mode], mode,
            baseline_gradient, [unit_gradient[:, :, :, j] for j in range(9)],
            mean.nu, h
        )
        print(json.dumps({"stage": "mode_matrix_ready", "mode": mode}), flush=True)

    # Transform the coefficient vector so the objective has an O(1) Gram
    # metric.  This is only an optimization-coordinate preconditioner; all
    # matrices and reported coefficients remain in physical units.
    optimization_rows = {}
    best_attempt = None
    diagnostic_best = None
    rng = np.random.default_rng(73013)
    for mode in modes:
        item = matrix_data[str(mode)]
        H0, Hincrements = item["H0"], item["H_increments"]
        Hscaled = Hincrements * coefficient_scale[:, None, None]
        node_scale = np.maximum(np.abs(baseline_node_lambda2), 1.0)
        values0, vector0 = spectrum(H0)
        growth_scale = max(abs(float(values0[-1])), OPTIMIZATION_TARGET)
        cone_floor = np.maximum(0.01 * baseline_node_lambda2, 1.0e-8)

        def evaluate(b):
            Htotal = H0 + np.einsum("j,jab->ab", b, Hscaled)
            values, vectors = spectrum(Htotal)
            return float(values[-1]), vectors[:, -1], Htotal

        def objective(b):
            return float(0.5 * b @ gram_scaled @ b + 1.0e-8 * 0.5 * b @ b)

        def objective_jac(b):
            return gram_scaled @ b + 1.0e-8 * b

        def growth_constraint(b):
            return (evaluate(b)[0] - OPTIMIZATION_TARGET) / growth_scale

        def growth_jac(b):
            _, vector, _ = evaluate(b)
            return np.real(np.einsum("a,jab,b->j", vector.conj(), Hscaled, vector)) / growth_scale

        def cone_constraints(b):
            actual = b * coefficient_scale
            return (candidate_node_data(node_cache, actual) - cone_floor) / node_scale

        seeds = [np.zeros(9)]
        first_gradient = np.real(np.einsum("a,jab,b->j", vector0[:, -1].conj(), Hscaled, vector0[:, -1]))
        if np.linalg.norm(first_gradient) > 0.0:
            direction = first_gradient / np.linalg.norm(first_gradient)
            step = OPTIMIZATION_TARGET / max(np.linalg.norm(first_gradient), 1.0)
            seeds.extend([direction * step * multiplier for multiplier in (0.1, 1.0, 10.0)])
        for _ in range(4):
            direction = rng.normal(size=9)
            direction /= max(np.linalg.norm(direction), np.finfo(float).tiny)
            seeds.append(direction * OPTIMIZATION_TARGET / max(np.linalg.norm(first_gradient), 1.0))

        attempts = []
        for seed in seeds:
            constraints = [
                {"type": "ineq", "fun": growth_constraint, "jac": growth_jac},
                {"type": "ineq", "fun": cone_constraints},
            ]
            result = minimize(
                objective,
                seed,
                jac=objective_jac,
                constraints=constraints,
                method="SLSQP",
                options={"maxiter": 300, "ftol": 1.0e-10, "disp": False},
            )
            lam, vector, total = evaluate(result.x)
            feasible = bool(
                lam >= TARGET_GROWTH * (1.0 - 1.0e-7)
                and np.min(cone_constraints(result.x)) >= -1.0e-7
            )
            attempt = {
                "constraint_mode": "growth_and_cone",
                "success": bool(result.success),
                "message": str(result.message),
                "iterations": int(result.nit),
                "objective": float(result.fun),
                "growth_lambda": lam,
                "feasible": feasible,
                "coefficient": (result.x * coefficient_scale).tolist(),
                "node_lambda2": candidate_node_data(node_cache, result.x * coefficient_scale).tolist(),
            }
            attempt["cone_min_margin"] = float(np.min(cone_constraints(result.x)))
            attempts.append(attempt)
            if diagnostic_best is None or attempt["cone_min_margin"] > diagnostic_best["cone_min_margin"]:
                diagnostic_best = dict(attempt, mode=mode, vector=vector, total=total)
            if feasible and (best_attempt is None or attempt["objective"] < best_attempt["objective"]):
                best_attempt = dict(attempt, mode=mode, vector=vector, total=total)
        optimization_rows[str(mode)] = {
            "attempts": attempts,
            "baseline_combined_lambda_max": float(values0[-1]),
        }
        print(json.dumps({"stage": "mode_optimized", "mode": mode,
                          "feasible_attempts": sum(a["feasible"] for a in attempts),
                          "best_lambda": max(a["growth_lambda"] for a in attempts)}), flush=True)

    if best_attempt is None and diagnostic_best is not None:
        best_attempt = dict(
            diagnostic_best,
            rejected=True,
            rejection_reason="No cone-constrained trial met both the growth target and the 1 percent baseline lambda2 floor; this is an optimizer result, not an infeasibility claim.",
            selection_reason="Largest minimum normalized cone margin among recorded constrained trials.",
        )

    selected = None
    if best_attempt is not None and best_attempt.get("feasible", False):
        mode = int(best_attempt["mode"])
        item = matrix_data[str(mode)]
        actual = np.asarray(best_attempt["coefficient"], dtype=float)
        lam, vector, total = mode_candidate(item, actual)
        potential_coefficients = item["whitening"] @ vector
        added_velocity = np.einsum("ncj,j->nc", unit_velocity, actual)
        baseline_rms = float(np.sqrt(np.sum(weights * np.sum(baseline_velocity**2, axis=1))))
        baseline_peak = float(np.max(np.linalg.norm(baseline_velocity, axis=1)))
        added_rms = float(np.sqrt(np.sum(weights * np.sum(added_velocity**2, axis=1))))
        added_peak = float(np.max(np.linalg.norm(added_velocity, axis=1)))
        selected = {
            "delta_state": actual.tolist(),
            "mode": mode,
            "growth_lambda": lam,
            "energy_rate_2lambda": 2.0 * lam,
            "added_swirl_rms": added_rms,
            "added_swirl_peak": added_peak,
            "baseline_velocity_rms": baseline_rms,
            "baseline_velocity_peak": baseline_peak,
            "potential_coefficients": pack_complex(potential_coefficients),
            "potential_coefficient_dimension": int(len(potential_coefficients)),
            "node_lambda2": candidate_node_data(node_cache, actual).tolist(),
            "growing_eigenvector_center_covariance": covariance_direction(
                center, widths, degree, carriers[mode], mode, potential_coefficients
            ).tolist(),
        }

    order24 = None
    if include_order24 and selected is not None:
        mode = int(selected["mode"])
        widths24 = widths
        points24, weights24 = quadrature(center, widths24, 24)
        baseline_gradient24 = spatial_gradient(mean, points24, tau, h)
        item24 = assemble_mode(
            mean, points24, weights24, center, widths24, degree, carriers[mode], mode,
            baseline_gradient24,
            [spatial_gradient(field, points24, tau, h) for field in unit_fields],
            mean.nu, h
        )
        actual = np.asarray(selected["delta_state"], dtype=float)
        lam24, vector24, _ = mode_candidate(item24, actual)
        selected24_coefficients = item24["whitening"] @ vector24
        velocity24 = np.stack([field.fields(points24, tau)[0] for field in unit_fields], axis=2)
        added24 = np.einsum("ncj,j->nc", velocity24, actual)
        selected["order24_replay"] = {
            "quadrature_order_per_axis": 24,
            "point_count": int(len(points24)),
            "growth_lambda": lam24,
            "energy_rate_2lambda": 2.0 * lam24,
            "added_swirl_rms": float(np.sqrt(np.sum(weights24 * np.sum(added24**2, axis=1)))),
            "potential_coefficients": pack_complex(selected24_coefficients),
        }

    matrix_cache = {
        "quadrature_order_per_axis": 16,
        "points": points.tolist(),
        "weights": weights.tolist(),
        "swirl_gram": gram.tolist(),
        "swirl_gram_regularization": gram_regularization,
        "coefficient_scale": coefficient_scale.tolist(),
        "modes": {
            mode: {
                "rank": item["rank"],
                "dimension": item["dimension"],
                "singular_max": item["singular_max"],
                "singular_min_retained": item["singular_min_retained"],
                "whitening": pack_complex(item["whitening"]),
                "H0": pack_complex(item["H0"]),
                "H_increments": pack_complex(item["H_increments"]),
            }
            for mode, item in matrix_data.items()
        },
    }
    result = {
        "accepted": False,
        "pde_validated": False,
        "scope": "Finite-matrix factor-8 mean-swirl co-design candidate only; no moment feasibility, full stress cone, time integration, nonlinear transfer, or Navier--Stokes acceptance.",
        "source_report": "fourier_patch_implicit.json",
        "mean_report": "outer_feedback_evolution.json",
        "initial_k": float(-np.log2(2.0 * tau)),
        "initial_physical_time": -tau,
        "center": list(center),
        "widths": list(widths),
        "base_widths": list(base_widths),
        "width_factor": WIDTH_FACTOR,
        "degree": degree,
        "modes": modes,
        "carriers": {str(m): carriers[m].tolist() for m in sorted(carriers)},
        "target_growth_lambda": TARGET_GROWTH,
        "optimization_target_lambda": OPTIMIZATION_TARGET,
        "optimization": {
            "coordinate": "actual_delta_state = coefficient_scale * optimization_coordinate",
            "seeds": "zero, first-order growth directions at multipliers 0.1/1/10, and four deterministic Gaussian directions",
            "cone_constraints": "22 solve_control eta/y locations with lambda2 >= 1 percent of the baseline lambda2 included in SLSQP and reported",
            "regularization": gram_regularization,
        },
        "baseline_diagnostic": {
            mode: {
                "rank": item["rank"],
                "diffusion_only": {
                    "lambda_min": float(item["diffusion_only"][0]),
                    "lambda_max": float(item["diffusion_only"][-1]),
                },
                "strain_only": {
                    "lambda_min": float(item["strain_only"][0]),
                    "lambda_max": float(item["strain_only"][-1]),
                },
                "combined_zero": {
                    "lambda_min": float(item["combined_zero"][0]),
                    "lambda_max": float(item["combined_zero"][-1]),
                },
            }
            for mode, item in matrix_data.items()
        },
        "baseline_node_lambda2": baseline_node_lambda2.tolist(),
        "node_labels": node_labels,
        "best_attempt": None if best_attempt is None else {
            key: value for key, value in best_attempt.items()
            if key not in ("vector", "total")
        },
        "selected_candidate": selected,
        "optimization_rows": optimization_rows,
        "matrix_cache": matrix_cache,
        "elapsed_seconds": float(time.perf_counter() - started),
    }
    Path(output_path).write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"stage": "saved", "output": str(output_path),
                      "selected": selected is not None}), flush=True)
    return result


if __name__ == "__main__":
    run()
