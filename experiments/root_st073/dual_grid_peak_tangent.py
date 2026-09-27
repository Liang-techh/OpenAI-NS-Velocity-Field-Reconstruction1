"""Bounded dual-grid peak-constrained tangent fit.

This diagnostic reuses the refined-grid L2 objective while adding the prior
order-13, 12-angle grid as a second peak constraint.  The holdout affine
model is anchored at the previously feasible candidate by reconstructing its
base residual from an actual five-point finite-difference replay.  A finite
union of worst-point peak constraints is expanded for at most five rounds;
every candidate is then checked against both complete sampled grids before it
can replace the old feasible seed.

The output is an assembled spatial diagnostic.  It does not certify a new
disjoint grid, a time trajectory, the PDE, or scale recursion.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import time
from pathlib import Path

for _name in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS"):
    os.environ[_name] = "1"

import numpy as np
from scipy.optimize import LinearConstraint, minimize


ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from affine_momentum import jets, momentum  # noqa: E402
from constrained_tangent_projection import quadratic_moment_target  # noqa: E402
from enriched_shape_replay import build_field  # noqa: E402
from enriched_shape_tangent import _weighted_constraint_setup  # noqa: E402
from wave_higher_harmonic_tangent import _mode_block_columns  # noqa: E402


DEFAULT_DENSE = ROOT / "full_wave_dense_tangent.json"
DEFAULT_REFINED_CACHE = ROOT / "refined_wave_momentum_cache.npz"
DEFAULT_REFINED_REPORT = ROOT / "refined_wave_momentum_cache.json"
DEFAULT_OUTPUT = ROOT / "dual_grid_peak_tangent.json"
DEFAULT_NPZ = ROOT / "dual_grid_peak_cache.npz"

OLD_CANDIDATE = ROOT / "enriched_mean_endpoint_tangent.json"
BALANCED_CANDIDATE = ROOT / "balanced_refined_tangent.json"
OLD_REPLAY = ROOT / "enriched_mean_momentum_replay.json"
BALANCED_REPLAY = ROOT / "balanced_refined_momentum_replay.json"
MEAN_ROWS = ROOT / "enriched_mean_constraint_rows.json"
MOMENT_REPORT = ROOT / "wave_dynamics_mean_compatibility.json"
CONE_REPORT = ROOT / "wave_mean_cone_projection.json"
ENDPOINT_REPORT = ROOT / "enriched_endpoint_shape_cache.json"
FROZEN_GEOMETRY = ROOT / "full_wave_frozen_cache.json"

HOLDOUT_ORDER = 13
HOLDOUT_ANGLES = 12
CONTROL_COUNT = 264
MAX_ROUNDS = 5


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _metric(residual, weights):
    residual = np.asarray(residual, dtype=float)
    weights = np.asarray(weights, dtype=float)
    magnitudes = np.linalg.norm(residual, axis=1)
    return {
        "point_count": int(len(residual)),
        "momentum_max": float(np.max(magnitudes)),
        "momentum_volume_L2": float(
            np.linalg.norm(residual.reshape(-1) * np.repeat(np.sqrt(weights), 3))
        ),
        "momentum_volume_RMS": float(
            np.sqrt(np.sum(weights * np.sum(residual * residual, axis=1)) / np.sum(weights))
        ),
        "physical_volume": float(np.sum(weights)),
    }


def _mixed_tangent_design(points, geometry):
    """Build the 264 analytic velocity/pressure-gradient response columns."""

    center = geometry["center"]
    widths = geometry["widths"]
    carrier = np.asarray(geometry["carrier"], dtype=float)
    mode0_full = _mode_block_columns(points, center, widths, 0, 3, np.zeros(2))
    mode0 = mode0_full[:, ::2]
    mode1 = _mode_block_columns(points, center, widths, 1, 2, carrier)
    mode2 = _mode_block_columns(points, center, widths, 2, 3, 2.0 * carrier)
    design = np.column_stack((mode0, mode1, mode2))
    if design.shape != (3 * len(points), CONTROL_COUNT):
        raise ValueError(f"Unexpected analytic holdout design shape {design.shape}")
    return design


def _actual_replay(candidate, points, snapshot):
    field, replay_snapshot = build_field(candidate)
    if replay_snapshot["inputs"]["mean"]["coefficients"] != snapshot["inputs"]["mean"]["coefficients"]:
        raise ValueError("Replay mean differs from frozen geometry mean")
    tau = float(replay_snapshot["inputs"]["mean"]["tau"])
    hspace = float(replay_snapshot["timesteps"]["hspace"])
    htime = float(replay_snapshot["timesteps"]["htime"])
    return momentum(jets(field, points, tau, hspace, htime)), {
        "tau": tau,
        "hspace": hspace,
        "htime": htime,
    }


def _load_json(path):
    return json.loads(path.read_text(encoding="utf-8"))


def run(
    dense_path=DEFAULT_DENSE,
    refined_cache_path=DEFAULT_REFINED_CACHE,
    output_path=DEFAULT_OUTPUT,
    npz_path=DEFAULT_NPZ,
    max_rounds=MAX_ROUNDS,
):
    started = time.perf_counter()
    dense_path = Path(dense_path)
    refined_cache_path = Path(refined_cache_path)
    output_path = Path(output_path)
    npz_path = Path(npz_path)
    candidate_old = _load_json(OLD_CANDIDATE)
    candidate_balanced = _load_json(BALANCED_CANDIDATE)
    old_replay_report = _load_json(OLD_REPLAY)
    balanced_replay_report = _load_json(BALANCED_REPLAY)
    dense = _load_json(dense_path)
    rows = _load_json(MEAN_ROWS)
    moment = _load_json(MOMENT_REPORT)["reusable_moment_linearization"]
    cone = _load_json(CONE_REPORT).get("reusable_cone_linearization", _load_json(CONE_REPORT))
    endpoint_report = _load_json(ENDPOINT_REPORT)
    frozen = _load_json(FROZEN_GEOMETRY)

    if candidate_old.get("status") != "completed" or not candidate_old.get("assembled_feasible"):
        raise ValueError("Old candidate must be completed and assembled feasible")
    if candidate_balanced.get("status") != "completed" or not candidate_balanced.get("assembled_feasible"):
        raise ValueError("Balanced candidate must be completed and assembled feasible")
    if old_replay_report.get("status") != "completed" or balanced_replay_report.get("status") != "completed":
        raise ValueError("Independent replay reports must be completed")
    if not dense.get("new_frozen_cache"):
        raise ValueError("Dense report has no new_frozen_cache")

    old_control = np.asarray(candidate_old["selected"]["tangent_coefficients"], dtype=float)
    balanced_control = np.asarray(candidate_balanced["selected"]["tangent_coefficients"], dtype=float)
    if old_control.shape != (CONTROL_COUNT,) or balanced_control.shape != (CONTROL_COUNT,):
        raise ValueError("Expected 264-control old and balanced candidates")
    old_wave = np.asarray(candidate_old["selected"]["coefficients_original"], dtype=float)
    balanced_wave = np.asarray(candidate_balanced["selected"]["coefficients_original"], dtype=float)
    if not np.array_equal(old_wave, balanced_wave):
        raise ValueError("Old and balanced candidates use different locked waves")
    dense_cache = dense["new_frozen_cache"]
    holdout_points = np.asarray(dense_cache["points"], dtype=float)
    holdout_weights = np.asarray(dense_cache["weights"], dtype=float)
    if holdout_points.shape != (18720, 3) or holdout_weights.shape != (18720,):
        raise ValueError(f"Unexpected holdout grid shapes {holdout_points.shape} {holdout_weights.shape}")

    with np.load(refined_cache_path, allow_pickle=False) as loaded:
        refined = {key: loaded[key] for key in loaded.files}
    refined_points = np.asarray(refined["points"], dtype=float)
    refined_weights = np.asarray(refined["weights"], dtype=float)
    refined_residual = np.asarray(refined["residual"], dtype=float)
    refined_design = np.asarray(refined["tangent_design"], dtype=float)
    if refined_points.shape != (44400, 3) or refined_design.shape != (133200, CONTROL_COUNT):
        raise ValueError(f"Unexpected refined cache shapes {refined_points.shape} {refined_design.shape}")
    if not np.array_equal(np.asarray(refined["coefficients_original"], dtype=float), old_wave):
        raise ValueError("Refined cache locked wave differs from candidates")

    report = {
        "status": "initializing",
        "accepted": False,
        "pde_validated": False,
        "scale_recursion_established": False,
        "constraints_maintained": False,
        "assembled_feasible": False,
        "scope": (
            "Dual-grid spatial tangent fit. The prior order-13, 12-angle grid "
            "is added as a finite peak-constraint training grid while the order-20 "
            "refined grid supplies the L2 objective. A new disjoint verification "
            "grid, trajectory, PDE, and scale recursion remain outside scope."
        ),
        "sources": {
            "old_candidate": {"path": OLD_CANDIDATE.name, "sha256": _sha(OLD_CANDIDATE)},
            "balanced_candidate": {"path": BALANCED_CANDIDATE.name, "sha256": _sha(BALANCED_CANDIDATE)},
            "old_replay": {"path": OLD_REPLAY.name, "sha256": _sha(OLD_REPLAY)},
            "balanced_replay": {"path": BALANCED_REPLAY.name, "sha256": _sha(BALANCED_REPLAY)},
            "dense_report": {"path": dense_path.name, "sha256": _sha(dense_path)},
            "refined_cache": {"path": refined_cache_path.name, "sha256": _sha(refined_cache_path)},
            "refined_cache_report": {"path": DEFAULT_REFINED_REPORT.name, "sha256": _sha(DEFAULT_REFINED_REPORT)},
            "mean_rows": {"path": MEAN_ROWS.name, "sha256": _sha(MEAN_ROWS)},
            "moment": {"path": MOMENT_REPORT.name, "sha256": _sha(MOMENT_REPORT)},
            "cone": {"path": CONE_REPORT.name, "sha256": _sha(CONE_REPORT)},
            "endpoint": {"path": ENDPOINT_REPORT.name, "sha256": _sha(ENDPOINT_REPORT)},
            "frozen_geometry": {"path": FROZEN_GEOMETRY.name, "sha256": _sha(FROZEN_GEOMETRY)},
        },
        "grid_definitions": {
            "refined": {"point_count": int(len(refined_points)), "design_shape": list(refined_design.shape),
                         "definition": "refined_wave_momentum_cache order-20, 12-angle split grid"},
            "holdout_promoted_to_training": {
                "point_count": int(len(holdout_points)), "angles": HOLDOUT_ANGLES,
                "order": HOLDOUT_ORDER,
                "definition": "full_wave_dense_tangent.new_frozen_cache order-13, 12-angle split grid",
            },
        },
        "analytic_design": {
            "method": "wave_higher_harmonic_tangent._mode_block_columns via supported_fourier_basis.basis_data",
            "layout": "mode0 degree-3 real [64], mode1 degree-2 complex [72], mode2 degree-3 complex [128]",
            "control_count": CONTROL_COUNT,
        },
    }
    output_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"stage": "sources_loaded", "holdout_points": len(holdout_points), "refined_points": len(refined_points)}), flush=True)

    # The holdout base is anchored at the previously feasible old candidate so
    # the fallback is an actual replay, not merely a linearized assertion.  We
    # also replay the balanced candidate to record the known peak regression.
    geometry = frozen["inputs"]["wave"]
    print(json.dumps({"stage": "holdout_actual_replays_started", "point_count": len(holdout_points)}), flush=True)
    actual_old, replay_steps = _actual_replay(candidate_old, holdout_points, frozen)
    actual_balanced, _ = _actual_replay(candidate_balanced, holdout_points, frozen)
    holdout_design = _mixed_tangent_design(holdout_points, geometry)
    holdout_base = actual_old - (holdout_design @ old_control).reshape(-1, 3)
    balanced_affine_holdout = holdout_base + (holdout_design @ balanced_control).reshape(-1, 3)
    old_replay_error = holdout_base + (holdout_design @ old_control).reshape(-1, 3) - actual_old
    report["actual_replays"] = {
        "old_candidate": {
            "metric": _metric(actual_old, holdout_weights),
            "reported_metric": old_replay_report.get("momentum"),
            "affine_anchor_error_max": float(np.max(np.abs(old_replay_error))),
        },
        "balanced_candidate": {
            "metric": _metric(actual_balanced, holdout_weights),
            "reported_metric": balanced_replay_report.get("momentum"),
            "affine_old_anchor_model_metric": _metric(balanced_affine_holdout, holdout_weights),
        },
        "replay_steps": replay_steps,
        "base_reconstruction": "actual old-candidate FD residual minus analytic D_holdout @ old_control",
    }
    report["status"] = "holdout_replays_complete"
    output_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"stage": "holdout_actual_replays_complete",
                      "old_l2": report["actual_replays"]["old_candidate"]["metric"]["momentum_volume_L2"],
                      "old_max": report["actual_replays"]["old_candidate"]["metric"]["momentum_max"],
                      "balanced_l2": report["actual_replays"]["balanced_candidate"]["metric"]["momentum_volume_L2"],
                      "balanced_max": report["actual_replays"]["balanced_candidate"]["metric"]["momentum_max"]}), flush=True)

    # The caps are the old candidate peaks on each grid, with the requested
    # 1e-6 numerical allowance.  They are read from the frozen reports rather
    # than inferred from an optimizer objective.
    old_peak_base = float(old_replay_report["momentum"]["momentum_max"])
    refined_peak_base = float(
        candidate_balanced["seed_metrics"]["actual_seed"]["momentum"]["momentum_max"]
    )
    old_peak_cap = old_peak_base * (1.0 + 1.0e-6)
    refined_peak_cap = refined_peak_base * (1.0 + 1.0e-6)
    report["peak_caps"] = {
        "old_holdout": {"base": old_peak_base, "relative_allowance": 1.0e-6, "cap": old_peak_cap},
        "refined": {"base": refined_peak_base, "relative_allowance": 1.0e-6, "cap": refined_peak_cap},
        "cap_definition": "old-candidate sampled maxima from frozen reports times 1.000001",
    }

    coeff = old_wave
    wave_x = np.r_[coeff[:, 0], coeff[:, 1]]
    E = np.column_stack((np.asarray(rows["moment_rows"], dtype=float), np.zeros((4, 200))))
    C = np.column_stack((np.asarray(rows["cone_control_rows"], dtype=float), np.zeros((81, 200))))
    if E.shape != (4, CONTROL_COUNT) or C.shape != (81, CONTROL_COUNT):
        raise ValueError("Unexpected moment/cone row shapes")
    target, _ = quadratic_moment_target(moment["baseline_moments"], moment["wave_moment_forms"], wave_x)
    cone_base = np.asarray(cone["cone_baseline"], float) + np.einsum(
        "i,kij,j->k", wave_x, np.asarray(cone["cone_wave_forms"], float), wave_x
    )
    cone_lower = np.asarray(cone["cone_lower"], float)
    setup = _weighted_constraint_setup(refined_design, E, target, refined_residual, refined_weights)
    particular = setup["particular"]
    mapping = setup["mapping"] * setup["scale"]
    weighted_design = setup["weighted_design"]
    sqrtw = setup["sqrt_weights"]
    shifted = refined_residual.reshape(-1) * sqrtw + weighted_design @ particular
    initial_z = setup["Q"].T @ (weighted_design @ (old_control - particular)) / setup["scale"]
    mapped_old = particular + mapping @ initial_z
    weighted_mapping = weighted_design @ mapping
    setup["weighted_design"] = None
    del weighted_design

    # Cone rows are normalized for SLSQP conditioning, without introducing an
    # arbitrary extra margin that could reject the known old seed.
    cone_matrix = C @ mapping
    cone_rhs = cone_lower - cone_base - C @ particular
    cone_norms = np.linalg.norm(cone_matrix, axis=1)
    cone_responsive = cone_norms > 1.0e-14
    if np.any(cone_rhs[~cone_responsive] > 1.0e-7):
        raise ValueError("A fixed cone row is infeasible for every reduced control")
    cone_matrix_scaled = cone_matrix[cone_responsive] / cone_norms[cone_responsive, None]
    cone_rhs_scaled = cone_rhs[cone_responsive] / cone_norms[cone_responsive]

    from enriched_endpoint_shape_cache import EnrichedEndpointShape

    oracle = EnrichedEndpointShape()
    if not np.array_equal(np.asarray(oracle.coefficients_original), old_wave):
        raise ValueError("Endpoint oracle locked wave differs from candidate")
    reference = np.asarray(oracle.reference, dtype=float)
    signs = np.array((-1.0, 1.0, np.sign(reference[2])), dtype=float)
    ref_scale = np.maximum(np.abs(reference), 1.0e-30)
    requested = np.array((1.0e-6, 1.0e-6, 1.0e-3), dtype=float)
    endpoint_tolerance = 1.0e-9
    moment_tolerance = 1.0e-5
    cone_tolerance = 1.0e-7
    objective_scale = max(_metric(refined_residual + (refined_design @ old_control).reshape(-1, 3), refined_weights)["momentum_volume_L2"] ** 2, 1.0)

    def control_from_z(z):
        return particular + mapping @ np.asarray(z, dtype=float)

    def endpoint_constraint(z):
        values, jac = oracle.evaluate(control_from_z(z))
        return signs * (values - reference) / ref_scale - requested, signs[:, None] * (jac @ mapping) / ref_scale[:, None]

    def evaluate(control):
        control = np.asarray(control, dtype=float)
        refined_res = refined_residual + (refined_design @ control).reshape(-1, 3)
        holdout_res = holdout_base + (holdout_design @ control).reshape(-1, 3)
        values, _ = oracle.evaluate(control)
        moment_error = E @ control - target
        cone_margins = cone_base + C @ control - cone_lower
        endpoint_margin = signs * (values - reference) / ref_scale - requested
        return {
            "control": control,
            "refined_residual": refined_res,
            "holdout_residual": holdout_res,
            "refined_metric": _metric(refined_res, refined_weights),
            "holdout_metric": _metric(holdout_res, holdout_weights),
            "moment_max_abs": float(np.max(np.abs(moment_error))),
            "cone_min_margin": float(np.min(cone_margins)),
            "cone_pass_count": int(np.sum(cone_margins >= -cone_tolerance)),
            "endpoint_observables": values,
            "endpoint_margin": endpoint_margin,
            "endpoint_min_margin": float(np.min(endpoint_margin)),
        }

    def feasible(item):
        return (
            item["moment_max_abs"] <= moment_tolerance
            and item["cone_min_margin"] >= -cone_tolerance
            and item["endpoint_min_margin"] >= -endpoint_tolerance
            and item["refined_metric"]["momentum_max"] <= refined_peak_cap * (1.0 + 1.0e-10)
            and item["holdout_metric"]["momentum_max"] <= old_peak_cap * (1.0 + 1.0e-10)
        )

    def weighted_objective(z):
        rw = shifted + weighted_mapping @ z
        return 0.5 * float(rw @ rw) / objective_scale

    def weighted_gradient(z):
        rw = shifted + weighted_mapping @ z
        return weighted_mapping.T @ rw / objective_scale

    def point_pool_data(z, grid_name, indices):
        control = control_from_z(z)
        if grid_name == "refined":
            base, design, cap = refined_residual, refined_design, refined_peak_cap
        else:
            base, design, cap = holdout_base, holdout_design, old_peak_cap
        indices = np.asarray(sorted(indices), dtype=int)
        rows_flat = np.concatenate([np.arange(3 * i, 3 * i + 3) for i in indices])
        point_design = design[rows_flat] @ mapping
        values = (base + (design @ control).reshape(-1, 3))[indices]
        margins = cap * cap - np.sum(values * values, axis=1)
        jac = np.empty((len(indices), mapping.shape[1]), dtype=float)
        for row, (value, block) in enumerate(zip(values, point_design.reshape(len(indices), 3, -1))):
            jac[row] = -2.0 * value @ block
        return margins, jac

    def grid_peak(item, grid_name):
        return item["refined_metric"]["momentum_max"] if grid_name == "refined" else item["holdout_metric"]["momentum_max"]

    seed = evaluate(old_control)
    mapped_seed = evaluate(mapped_old)
    pool_refined = set(np.argsort(np.linalg.norm(seed["refined_residual"], axis=1))[-20:].tolist())
    pool_holdout = set(np.argsort(np.linalg.norm(seed["holdout_residual"], axis=1))[-20:].tolist())
    report["reduced_rank"] = int(setup["keep_rank"])
    report["rank_whitened_seed"] = {
        "mapped_control_max_error": float(np.max(np.abs(mapped_old - old_control))),
        "old_seed": {
            "refined": seed["refined_metric"], "holdout": seed["holdout_metric"],
            "moment_max_abs": seed["moment_max_abs"], "cone_min_margin": seed["cone_min_margin"],
            "endpoint_min_margin": seed["endpoint_min_margin"], "feasible": bool(feasible(seed)),
        },
        "mapped_old": {
            "refined": mapped_seed["refined_metric"], "holdout": mapped_seed["holdout_metric"],
            "moment_max_abs": mapped_seed["moment_max_abs"], "cone_min_margin": mapped_seed["cone_min_margin"],
            "endpoint_min_margin": mapped_seed["endpoint_min_margin"], "feasible": bool(feasible(mapped_seed)),
        },
    }
    report["status"] = "analytic_dual_grid_setup_complete"
    output_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    np.savez_compressed(
        npz_path,
        points_holdout=holdout_points,
        weights_holdout=holdout_weights,
        base_residual_holdout=holdout_base,
        tangent_design_holdout=holdout_design,
        actual_old_holdout=actual_old,
        actual_balanced_holdout=actual_balanced,
        old_control=old_control,
        balanced_control=balanced_control,
        source_hash_dense=np.asarray([_sha(dense_path)]),
        source_hash_old=np.asarray([_sha(OLD_CANDIDATE)]),
        source_hash_balanced=np.asarray([_sha(BALANCED_CANDIDATE)]),
    )
    report["dual_grid_cache"] = {
        "path": npz_path.name,
        "sha256": _sha(npz_path),
        "size_bytes": int(npz_path.stat().st_size),
        "keys": {
            "points_holdout": list(holdout_points.shape),
            "weights_holdout": list(holdout_weights.shape),
            "base_residual_holdout": list(holdout_base.shape),
            "tangent_design_holdout": list(holdout_design.shape),
        },
    }
    report["cutting_plane_rounds"] = []
    report["status"] = "optimization_started"
    output_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")

    best = seed if feasible(seed) else mapped_seed
    best_z = initial_z.copy()
    best_reason = "old_feasible_seed"
    best_peak = best
    for round_index in range(int(max_rounds)):
        def pool_fun_refined(z):
            return point_pool_data(z, "refined", pool_refined)[0]

        def pool_jac_refined(z):
            return point_pool_data(z, "refined", pool_refined)[1]

        def pool_fun_holdout(z):
            return point_pool_data(z, "holdout", pool_holdout)[0]

        def pool_jac_holdout(z):
            return point_pool_data(z, "holdout", pool_holdout)[1]

        constraints = [
            LinearConstraint(cone_matrix_scaled, cone_rhs_scaled, np.inf),
            {"type": "ineq", "fun": lambda z: endpoint_constraint(z)[0],
             "jac": lambda z: endpoint_constraint(z)[1]},
            {"type": "ineq", "fun": pool_fun_refined, "jac": pool_jac_refined},
            {"type": "ineq", "fun": pool_fun_holdout, "jac": pool_jac_holdout},
        ]
        fit = minimize(
            weighted_objective,
            best_z,
            jac=weighted_gradient,
            method="SLSQP",
            constraints=constraints,
            options=dict(maxiter=80, ftol=1.0e-10, disp=False),
        )
        candidate_z = np.asarray(fit.x if fit.x is not None else best_z, dtype=float)
        candidate = evaluate(control_from_z(candidate_z))
        mags_refined = np.linalg.norm(candidate["refined_residual"], axis=1)
        mags_holdout = np.linalg.norm(candidate["holdout_residual"], axis=1)
        violating_refined = np.flatnonzero(mags_refined > refined_peak_cap * (1.0 + 1.0e-10))
        violating_holdout = np.flatnonzero(mags_holdout > old_peak_cap * (1.0 + 1.0e-10))
        added_refined = [int(i) for i in violating_refined[np.argsort(mags_refined[violating_refined])[-10:][::-1]] if int(i) not in pool_refined]
        added_holdout = [int(i) for i in violating_holdout[np.argsort(mags_holdout[violating_holdout])[-10:][::-1]] if int(i) not in pool_holdout]
        pool_refined.update(added_refined)
        pool_holdout.update(added_holdout)
        candidate_is_feasible = bool(feasible(candidate))
        if candidate_is_feasible and candidate["refined_metric"]["momentum_volume_L2"] < best["refined_metric"]["momentum_volume_L2"]:
            best = candidate
            best_z = candidate_z.copy()
            best_reason = "feasible_lower_refined_L2"
        if candidate_is_feasible and candidate["holdout_metric"]["momentum_max"] < best_peak["holdout_metric"]["momentum_max"]:
            best_peak = candidate
        report["cutting_plane_rounds"].append({
            "round": round_index + 1,
            "pool_refined": len(pool_refined),
            "pool_holdout": len(pool_holdout),
            "optimizer_success": bool(fit.success),
            "optimizer_status": int(fit.status),
            "optimizer_message": str(fit.message),
            "optimizer_iterations": int(getattr(fit, "nit", -1)),
            "candidate_feasible": candidate_is_feasible,
            "candidate_refined_metric": candidate["refined_metric"],
            "candidate_holdout_metric": candidate["holdout_metric"],
            "candidate_moment_max_abs": candidate["moment_max_abs"],
            "candidate_cone_min_margin": candidate["cone_min_margin"],
            "candidate_endpoint_min_margin": candidate["endpoint_min_margin"],
            "candidate_refined_global_peak": float(np.max(mags_refined)),
            "candidate_holdout_global_peak": float(np.max(mags_holdout)),
            "newly_added_refined_points": added_refined,
            "newly_added_holdout_points": added_holdout,
        })
        output_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
        print(json.dumps({"stage": "cutting_plane_round", "round": round_index + 1,
                          "success": bool(fit.success), "feasible": candidate_is_feasible,
                          "refined_l2": candidate["refined_metric"]["momentum_volume_L2"],
                          "refined_max": candidate["refined_metric"]["momentum_max"],
                          "holdout_max": candidate["holdout_metric"]["momentum_max"],
                          "added_refined": len(added_refined), "added_holdout": len(added_holdout)}), flush=True)
        if candidate_is_feasible and not added_refined and not added_holdout:
            break
        if not added_refined and not added_holdout and not candidate_is_feasible:
            break
        if candidate_is_feasible and candidate["refined_metric"]["momentum_volume_L2"] <= best["refined_metric"]["momentum_volume_L2"]:
            best_z = candidate_z.copy()

    selected = best if feasible(best) else seed
    selected_reason = best_reason if feasible(best) else "old_feasible_seed_retained"
    report["selected"] = {
        "reason": selected_reason,
        "tangent_coefficients": selected["control"].tolist(),
        "training_refined_momentum": selected["refined_metric"],
        "training_holdout_momentum": selected["holdout_metric"],
        "refined_l2_improvement_vs_old_seed": float(seed["refined_metric"]["momentum_volume_L2"] - selected["refined_metric"]["momentum_volume_L2"]),
        "refined_max_change_vs_old_seed": float(selected["refined_metric"]["momentum_max"] - seed["refined_metric"]["momentum_max"]),
        "holdout_max_change_vs_old_seed": float(selected["holdout_metric"]["momentum_max"] - seed["holdout_metric"]["momentum_max"]),
        "momentum_max_caps": {"refined": refined_peak_cap, "holdout": old_peak_cap},
        "assembled_moment_max_abs": selected["moment_max_abs"],
        "assembled_cone_min_margin": selected["cone_min_margin"],
        "assembled_cone_pass_count": selected["cone_pass_count"],
        "endpoint_observables": selected["endpoint_observables"].tolist(),
        "endpoint_margin": selected["endpoint_margin"].tolist(),
        "endpoint_min_margin": selected["endpoint_min_margin"],
        "coefficients_original": candidate_old["selected"]["coefficients_original"],
        "coefficients_whitened": candidate_old["selected"]["coefficients_whitened"],
        "assembled_feasible": bool(feasible(selected)),
    }
    report["assembled_feasible"] = bool(feasible(selected))
    report["assembled_constraints_maintained"] = bool(feasible(selected))
    report["constraints_maintained_scope"] = (
        "Both complete sampled grids plus moments, cones, and endpoint geometry "
        "were checked for the selected affine assembled iterate; no independent "
        "new-grid or trajectory validation was performed."
    )
    report["independent_validation_pending"] = True
    report["status"] = "completed"
    report["optimizer_scope"] = {
        "max_rounds": int(max_rounds),
        "rounds_run": len(report["cutting_plane_rounds"]),
        "pool_refined_final": len(pool_refined),
        "pool_holdout_final": len(pool_holdout),
        "objective": "refined-grid weighted physical volume L2",
    }
    report["elapsed_seconds"] = float(time.perf_counter() - started)
    output_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": report["status"],
        "reason": selected_reason,
        "old_refined_l2": seed["refined_metric"]["momentum_volume_L2"],
        "selected_refined_l2": selected["refined_metric"]["momentum_volume_L2"],
        "old_refined_max": seed["refined_metric"]["momentum_max"],
        "selected_refined_max": selected["refined_metric"]["momentum_max"],
        "old_holdout_max": seed["holdout_metric"]["momentum_max"],
        "selected_holdout_max": selected["holdout_metric"]["momentum_max"],
        "rounds": len(report["cutting_plane_rounds"]),
        "elapsed_seconds": report["elapsed_seconds"],
    }), flush=True)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--dense", type=Path, default=DEFAULT_DENSE)
    parser.add_argument("--refined-cache", type=Path, default=DEFAULT_REFINED_CACHE)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--npz", type=Path, default=DEFAULT_NPZ)
    parser.add_argument("--max-rounds", type=int, default=MAX_ROUNDS)
    args = parser.parse_args()
    run(args.dense, args.refined_cache, args.output, args.npz, args.max_rounds)
