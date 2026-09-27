"""Bounded refined-grid momentum fit with a sampled peak-norm cutting plane.

The reduced variables preserve the four mean moments exactly.  Every iterate
is checked against all 81 cone rows and the nonlinear endpoint shape oracle.
The pointwise momentum bound starts with the 20 largest seed points and adds
worst violators for at most five rounds, keeping the best actually feasible
iterate rather than trusting an optimizer status flag.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import time
from pathlib import Path

for _name in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS"):
    os.environ[_name] = "1"

import numpy as np
from scipy.optimize import LinearConstraint, minimize

from constrained_tangent_projection import quadratic_moment_target
from enriched_shape_tangent import _weighted_constraint_setup


ROOT = Path(__file__).resolve().parent
DEFAULT_CACHE = ROOT / "refined_wave_momentum_cache.npz"
DEFAULT_OUTPUT = ROOT / "balanced_refined_tangent.json"


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _metric(residual, weights):
    residual = np.asarray(residual, float)
    weights = np.asarray(weights, float)
    magnitudes = np.linalg.norm(residual, axis=1)
    return dict(
        point_count=int(len(residual)),
        momentum_max=float(np.max(magnitudes)),
        momentum_volume_L2=float(np.linalg.norm(residual.reshape(-1) * np.repeat(np.sqrt(weights), 3))),
        momentum_volume_RMS=float(np.sqrt(np.sum(weights * np.sum(residual * residual, axis=1)) / np.sum(weights))),
        physical_volume=float(np.sum(weights)),
    )


def _replay_candidate(control, residual, design, weights, E, target, C, cone_base,
                      cone_lower, oracle, reference, signs, ref_scale, requested):
    control = np.asarray(control, float)
    point_residual = residual + (design @ control).reshape(-1, 3)
    metrics = _metric(point_residual, weights)
    moment_error = E @ control - target
    cone_margins = cone_base + C @ control - cone_lower
    values, _ = oracle.evaluate(control)
    endpoint_margin = signs * (values - reference) / ref_scale - requested
    return dict(
        control=control,
        residual=point_residual,
        metrics=metrics,
        moment_max_abs=float(np.max(np.abs(moment_error))),
        cone_min_margin=float(np.min(cone_margins)),
        cone_pass_count=int(np.sum(cone_margins >= -1.0e-7)),
        endpoint_observables=values,
        endpoint_margin=endpoint_margin,
        endpoint_min_margin=float(np.min(endpoint_margin)),
    )


def run(cache_path=DEFAULT_CACHE, output_path=DEFAULT_OUTPUT, max_rounds=5):
    started = time.perf_counter()
    cache_path = Path(cache_path)
    output_path = Path(output_path)
    paths = {
        "seed": ROOT / "enriched_mean_endpoint_tangent.json",
        "mean_rows": ROOT / "enriched_mean_constraint_rows.json",
        "moment": ROOT / "wave_dynamics_mean_compatibility.json",
        "cone": ROOT / "wave_mean_cone_projection.json",
        "endpoint": ROOT / "enriched_endpoint_shape_cache.json",
        "geometry": ROOT / "full_wave_frozen_cache.json",
    }
    raw = {name: path.read_bytes() for name, path in paths.items()}
    data = {name: json.loads(value) for name, value in raw.items()}
    if data["seed"].get("status") != "completed":
        raise ValueError("Seed report must be completed")
    if any(data[name].get("status") != "completed" for name in ("mean_rows", "moment", "cone", "endpoint")):
        raise ValueError("Constraint and endpoint reports must be completed")
    with np.load(cache_path, allow_pickle=False) as loaded:
        cache = {key: loaded[key] for key in loaded.files}
    points = np.asarray(cache["points"], dtype=float)
    weights = np.asarray(cache["weights"], dtype=float)
    residual = np.asarray(cache["residual"], dtype=float)
    design = np.asarray(cache["tangent_design"], dtype=float)
    if points.shape != (44400, 3) or design.shape != (133200, 264):
        raise ValueError(f"Unexpected refined cache shapes {points.shape} {design.shape}")

    from enriched_endpoint_shape_cache import EnrichedEndpointShape

    oracle = EnrichedEndpointShape()
    seed = data["seed"]["selected"]
    seed_control = np.asarray(seed["tangent_coefficients"], dtype=float)
    if seed_control.shape != (264,):
        raise ValueError("Expected 264-control enriched seed")
    if not np.array_equal(np.asarray(oracle.coefficients_original), np.asarray(seed["coefficients_original"])):
        raise ValueError("Endpoint initial wave differs from enriched seed")
    if "coefficients_original" in cache and not np.array_equal(cache["coefficients_original"], seed["coefficients_original"]):
        raise ValueError("Refined cache initial wave differs from enriched seed")

    rows = data["mean_rows"]
    moment = data["moment"]["reusable_moment_linearization"]
    cone = data["cone"].get("reusable_cone_linearization", data["cone"])
    coeff = np.asarray(seed["coefficients_original"], dtype=float)
    wave_x = np.r_[coeff[:, 0], coeff[:, 1]]
    E = np.column_stack((np.asarray(rows["moment_rows"], dtype=float), np.zeros((4, 200))))
    C = np.column_stack((np.asarray(rows["cone_control_rows"], dtype=float), np.zeros((81, 200))))
    if E.shape != (4, 264) or C.shape != (81, 264):
        raise ValueError("Unexpected enriched moment/cone row shapes")
    target, _ = quadratic_moment_target(moment["baseline_moments"], moment["wave_moment_forms"], wave_x)
    cone_base = np.asarray(cone["cone_baseline"], float) + np.einsum(
        "i,kij,j->k", wave_x, np.asarray(cone["cone_wave_forms"], float), wave_x
    )
    cone_lower = np.asarray(cone["cone_lower"], float)
    setup = _weighted_constraint_setup(design, E, target, residual, weights)
    particular = setup["particular"]
    mapping = setup["mapping"] * setup["scale"]
    weighted_design = setup["weighted_design"]
    sqrtw = setup["sqrt_weights"]
    shifted = residual.reshape(-1) * sqrtw + weighted_design @ particular
    # Rank-whitened seed recovery is deliberate: direct physical lstsq is
    # ill-conditioned for the 252-dimensional reduced tangent basis.
    initial_z = setup["Q"].T @ (weighted_design @ (seed_control - particular)) / setup["scale"]
    mapped_seed = particular + mapping @ initial_z
    weighted_mapping = (weighted_design @ mapping)
    # The setup retains the same 133200 by 264 array.  Drop that reference now
    # that the reduced objective matrix has been formed, keeping the bounded
    # solve within the memory needed by the refined cache and design.
    setup["weighted_design"] = None
    del weighted_design

    reference = np.asarray(oracle.reference, dtype=float)
    signs = np.array((-1.0, 1.0, np.sign(reference[2])), dtype=float)
    ref_scale = np.maximum(np.abs(reference), 1.0e-30)
    requested = np.array((1.0e-6, 1.0e-6, 1.0e-3), dtype=float)
    endpoint_tolerance = 1.0e-9
    moment_tolerance = 1.0e-5
    cone_tolerance = 1.0e-7
    seed_replay = _replay_candidate(
        mapped_seed, residual, design, weights, E, target, C, cone_base,
        cone_lower, oracle, reference, signs, ref_scale, requested
    )
    actual_seed_replay = _replay_candidate(
        seed_control, residual, design, weights, E, target, C, cone_base,
        cone_lower, oracle, reference, signs, ref_scale, requested
    )
    peak_bound = actual_seed_replay["metrics"]["momentum_max"] * (1.0 + 1.0e-6)
    objective_scale = max(seed_replay["metrics"]["momentum_volume_L2"] ** 2, 1.0)
    # The all-point seed peak is the bound; its 20 largest points initialize
    # the cutting plane pool.
    seed_magnitudes = np.linalg.norm(seed_replay["residual"], axis=1)
    pool = set(np.argsort(seed_magnitudes)[-20:].tolist())

    # Normalize the 81 cone rows before passing them to SLSQP.  Their physical
    # row norms differ by two orders of magnitude, while the constraint itself
    # remains the exact frozen lower bound.  The seed was constructed at about
    # 1e-4 above that lower bound, so adding another arbitrary 1e-4 margin would
    # incorrectly discard the known feasible starting point.
    cone_matrix = C @ mapping
    cone_rhs = cone_lower - cone_base - C @ particular
    cone_row_norm = np.linalg.norm(cone_matrix, axis=1)
    cone_responsive = cone_row_norm > 1.0e-14
    cone_fixed_violation = cone_rhs[~cone_responsive]
    if np.any(cone_fixed_violation > cone_tolerance):
        raise ValueError("A frozen cone row is infeasible for every reduced control")
    cone_matrix_scaled = cone_matrix[cone_responsive] / cone_row_norm[cone_responsive, None]
    cone_rhs_scaled = cone_rhs[cone_responsive] / cone_row_norm[cone_responsive]

    def control_from_z(z):
        return particular + mapping @ np.asarray(z, float)

    def weighted_objective(z):
        rw = shifted + weighted_mapping @ z
        return 0.5 * float(rw @ rw) / objective_scale

    def weighted_gradient(z):
        rw = shifted + weighted_mapping @ z
        return (weighted_mapping.T @ rw) / objective_scale

    def endpoint_constraint(z):
        values, jac = oracle.evaluate(control_from_z(z))
        return signs * (values - reference) / ref_scale - requested, (
            signs[:, None] * (jac @ mapping) / ref_scale[:, None]
        )

    def cone_constraint(z):
        return C @ control_from_z(z) + cone_base - cone_lower

    def pool_data(z, indices):
        control = control_from_z(z)
        point_residual = residual + (design @ control).reshape(-1, 3)
        indices = np.asarray(sorted(indices), dtype=int)
        rows = np.concatenate([np.arange(3 * i, 3 * i + 3) for i in indices])
        point_design = design[rows] @ mapping
        values = point_residual[indices]
        margins = peak_bound ** 2 - np.sum(values * values, axis=1)
        jac = np.empty((len(indices), mapping.shape[1]), dtype=float)
        for row, (value, block) in enumerate(zip(values, point_design.reshape(len(indices), 3, -1))):
            jac[row] = -2.0 * value @ block
        return margins, jac

    def peak_eval(z):
        control = control_from_z(z)
        point_residual = residual + (design @ control).reshape(-1, 3)
        magnitudes = np.linalg.norm(point_residual, axis=1)
        return magnitudes, point_residual

    def feasible(replay):
        return (
            replay["moment_max_abs"] <= moment_tolerance
            and replay["cone_min_margin"] >= -cone_tolerance
            and replay["endpoint_min_margin"] >= -endpoint_tolerance
            and replay["metrics"]["momentum_max"] <= peak_bound * (1.0 + 1.0e-10)
        )

    source_name = data["seed"].get("source", paths["seed"].name)
    source_path = ROOT / source_name
    if not source_path.exists():
        raise FileNotFoundError(f"Seed source artifact is missing: {source_path}")
    report = dict(
        status="running", accepted=False, pde_validated=False,
        scale_recursion_established=False, constraints_maintained=False,
        source=source_name,
        source_sha256=_sha(source_path),
        sources={name: {"path": path.name, "sha256": _sha(path)} for name, path in paths.items()},
        frozen_provenance={
            "endpoint_source": data["endpoint"].get("source"),
            "endpoint_source_sha256": data["endpoint"].get("source_sha256"),
            "mean_source": data["endpoint"].get("mean_source"),
            "mean_source_sha256": data["endpoint"].get("mean_source_sha256"),
            "geometry_source": data["endpoint"].get("geometry_source"),
            "geometry_sha256": data["endpoint"].get("geometry_sha256"),
            "endpoint_cache_sha256": data["endpoint"].get("cache_sha256"),
            "endpoint_grid_hash": data["endpoint"].get("grid_hash"),
        },
        momentum_cache={"path": cache_path.name, "sha256": _sha(cache_path),
                        "point_count": int(len(points)), "design_shape": list(design.shape)},
        control_count=264, reduced_rank=int(setup["keep_rank"]),
        cone_constraint_rows={"total": int(len(C)), "responsive": int(np.sum(cone_responsive)),
                              "fixed": int(np.sum(~cone_responsive)),
                              "row_norm_min": float(np.min(cone_row_norm)),
                              "row_norm_max": float(np.max(cone_row_norm)),
                              "enforced_lower_bound": "cone_lower"},
        peak_bound=float(peak_bound), peak_bound_relative_to_seed=1.0e-6,
        tolerances={"moment_max_abs": moment_tolerance, "cone_min_margin": cone_tolerance,
                    "endpoint_min_margin": endpoint_tolerance},
        requested_signed_fractional_endpoint_change=requested.tolist(),
        seed_metrics={"actual_seed": {
            "momentum": actual_seed_replay["metrics"],
            "moment_max": actual_seed_replay["metrics"]["momentum_max"],
            "moment_max_bound": peak_bound,
        }, "rank_whitened_mapped_seed": {
            "momentum": seed_replay["metrics"],
            "mapped_control_max_error": float(np.max(np.abs(mapped_seed - seed_control))),
            "moment_max": seed_replay["metrics"]["momentum_max"],
        }},
        initial_constraints={
            "moment_max_abs": seed_replay["moment_max_abs"],
            "cone_min_margin": seed_replay["cone_min_margin"],
            "endpoint_min_margin": seed_replay["endpoint_min_margin"],
        },
        cutting_plane_rounds=[],
        scope=(
            "Refined 44400-point frozen-grid tangent fit with four moment equalities, "
            "81 cones, nonlinear endpoint geometry, and a finite cutting-plane pool "
            "for the global momentum peak bound. No independent PDE or trajectory replay."
        ),
    )

    best = seed_replay
    best_z = initial_z.copy()
    best_reason = "rank_whitened_seed"
    # Keep a second tracker in case a candidate lowers the peak while preserving
    # all constraints but does not lower L2.
    best_peak = seed_replay
    best_peak_z = initial_z.copy()

    for round_index in range(int(max_rounds)):
        indices = sorted(pool)

        def pool_fun(z):
            return pool_data(z, indices)[0]

        def pool_jac(z):
            return pool_data(z, indices)[1]

        constraints = [
            LinearConstraint(cone_matrix_scaled, cone_rhs_scaled, np.inf),
            {"type": "ineq", "fun": lambda z: endpoint_constraint(z)[0],
             "jac": lambda z: endpoint_constraint(z)[1]},
            {"type": "ineq", "fun": pool_fun, "jac": pool_jac},
        ]
        fit = minimize(
            weighted_objective, best_z, jac=weighted_gradient,
            method="SLSQP", constraints=constraints,
            options=dict(maxiter=80, ftol=1.0e-10, disp=False),
        )
        candidate_z = np.asarray(fit.x if fit.x is not None else best_z, dtype=float)
        candidate = _replay_candidate(
            control_from_z(candidate_z), residual, design, weights, E, target, C,
            cone_base, cone_lower, oracle, reference, signs, ref_scale, requested
        )
        magnitudes, _ = peak_eval(candidate_z)
        new_violations = np.flatnonzero(magnitudes > peak_bound * (1.0 + 1.0e-10))
        if feasible(candidate):
            if candidate["metrics"]["momentum_volume_L2"] < best["metrics"]["momentum_volume_L2"]:
                best = candidate
                best_z = candidate_z.copy()
                best_reason = "feasible_lower_L2"
            if candidate["metrics"]["momentum_max"] < best_peak["metrics"]["momentum_max"]:
                best_peak = candidate
                best_peak_z = candidate_z.copy()
        else:
            # A non-feasible optimizer output is retained only for diagnostics;
            # it cannot replace the frozen feasible seed.
            pass
        worst = np.argsort(magnitudes)[-10:][::-1]
        added = [int(index) for index in worst if int(index) not in pool]
        for index in added:
            pool.add(index)
        report["cutting_plane_rounds"].append(dict(
            round=round_index + 1, pool_size=len(pool),
            optimizer_success=bool(fit.success), optimizer_status=int(fit.status),
            optimizer_message=str(fit.message), optimizer_iterations=int(getattr(fit, "nit", -1)),
            candidate_feasible=bool(feasible(candidate)),
            candidate_metrics=candidate["metrics"],
            candidate_moment_max_abs=candidate["moment_max_abs"],
            candidate_cone_min_margin=candidate["cone_min_margin"],
            candidate_endpoint_min_margin=candidate["endpoint_min_margin"],
            candidate_global_peak=float(np.max(magnitudes)),
            peak_bound=float(peak_bound),
            newly_added_worst_points=added,
            pool_worst_point=int(np.argmax(magnitudes)),
        ))
        # Stop once the candidate is feasible and no unrepresented global peak
        # violations remain; otherwise continue cutting planes up to the bound.
        if feasible(candidate) and len(new_violations) == 0:
            break
        if not added and len(new_violations) == 0:
            break
        # Warm start the next pool round from the best feasible point.
        if feasible(candidate) and candidate["metrics"]["momentum_volume_L2"] <= best["metrics"]["momentum_volume_L2"]:
            best_z = candidate_z.copy()

    # Explicit final replay of the selected feasible iterate.
    selected = best
    selected_z = best_z
    if not feasible(selected):
        selected = seed_replay
        selected_z = initial_z
        best_reason = "seed_retained_no_feasible_improvement"
    report["selected"] = dict(
        reason=best_reason,
        tangent_coefficients=selected["control"].tolist(),
        training_momentum=selected["metrics"],
        momentum_l2_improvement=float(seed_replay["metrics"]["momentum_volume_L2"] - selected["metrics"]["momentum_volume_L2"]),
        momentum_max_change=float(selected["metrics"]["momentum_max"] - seed_replay["metrics"]["momentum_max"]),
        momentum_max_bound=float(peak_bound),
        assembled_moment_max_abs=selected["moment_max_abs"],
        assembled_cone_min_margin=selected["cone_min_margin"],
        assembled_cone_pass_count=selected["cone_pass_count"],
        endpoint_observables=selected["endpoint_observables"].tolist(),
        endpoint_margin=selected["endpoint_margin"].tolist(),
        endpoint_min_margin=selected["endpoint_min_margin"],
        coefficients_original=seed["coefficients_original"],
        coefficients_whitened=seed["coefficients_whitened"],
    )
    selected_assembled_feasible = bool(feasible(selected))
    report["selected"]["assembled_feasible"] = selected_assembled_feasible
    report["assembled_feasible"] = selected_assembled_feasible
    report["assembled_constraints_maintained"] = selected_assembled_feasible
    report["constraints_maintained_scope"] = (
        "Frozen-grid assembled replay only (moments, cones, endpoint geometry, and "
        "all 44400 peak values); no independent trajectory validation."
    )
    report["independent_constraints_validated"] = False
    report["best_peak_feasible"] = dict(
        momentum=best_peak["metrics"],
        momentum_max_change=float(best_peak["metrics"]["momentum_max"] - seed_replay["metrics"]["momentum_max"]),
        moment_max_abs=best_peak["moment_max_abs"],
        cone_min_margin=best_peak["cone_min_margin"],
        endpoint_min_margin=best_peak["endpoint_min_margin"],
    )
    report["status"] = "completed"
    report["optimizer_scope"] = dict(max_rounds=int(max_rounds), rounds_run=len(report["cutting_plane_rounds"]),
                                      pool_final=len(pool), selected_z_norm=float(np.linalg.norm(selected_z)))
    report["elapsed_seconds"] = float(time.perf_counter() - started)
    output_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": report["status"], "reason": best_reason,
        "seed_l2": seed_replay["metrics"]["momentum_volume_L2"],
        "selected_l2": selected["metrics"]["momentum_volume_L2"],
        "seed_max": seed_replay["metrics"]["momentum_max"],
        "selected_max": selected["metrics"]["momentum_max"],
        "rounds": len(report["cutting_plane_rounds"]),
        "elapsed_seconds": report["elapsed_seconds"],
    }), flush=True)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--refined-cache", type=Path, default=DEFAULT_CACHE)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--max-rounds", type=int, default=5)
    args = parser.parse_args()
    run(args.refined_cache, args.output, args.max_rounds)
