"""Audit the local Pareto cone spanned by two already-screened ST052-M directions.

Preregistered in issue #641 before evaluation.  This round adds no new velocity
basis and constructs no combined child.  It reruns the exact Agent-7 compact
pure-swirl screen (#621) and shoulder-relative timing screen (#634), forms the
centered local Jacobian in the four previously frozen desirability coordinates,
and asks by deterministic linear programming whether their bounded local span
contains a target-free Pareto direction.

This is first-order expression-capacity evidence only.  It is not a visual
correspondence test or Navier--Stokes validation.
"""
from __future__ import annotations

import argparse
import json
import tempfile
from pathlib import Path

import numpy as np
from scipy.optimize import linprog

import agent7_st052m_compact_swirl_screen as swirl
import agent7_st052m_shoulder_temporal_phase as shoulder

TASK_ID = "CR003-ST052M-TWO-DIRECTION-PARETO-CONE-091"
PREREG_ISSUE = 641
SOURCE_TEMPORAL_PR = 587
SOURCE_TEMPORAL_HEAD = "0b93819095f6c8576a7571bdc2d2fbef4154944d"
SOURCE_SWIRL_PR = 621
SOURCE_SWIRL_HEAD = "d01ea33e4e18c0fe8f3a40baaf22b233a66fe96a"
SOURCE_SHOULDER_PR = 634
SOURCE_SHOULDER_HEAD = "9c628326cad1b8637b13016ffefb945901c4dd53"
SOURCE_RENDER_PR = 601
SOURCE_RENDER_HEAD = "02fb726e03268d17f3ecd1a7d503282a225d50eb"

PARETO_TOL = 5.0e-5
CONDITION_MAX = 25.0
METRICS = (
    "aspect_gain",
    "tip_thinning",
    "turns_fidelity",
    "axial_pair_fidelity",
)

TRUTH = {
    "canonical_velocity_changed": False,
    "saved_velocity_changed": False,
    "production_candidate_selected": False,
    "new_spatial_basis_added": False,
    "new_temporal_basis_added": False,
    "combined_velocity_child_constructed": False,
    "parameter_grid_scan_performed": False,
    "pressure_or_force_changed": False,
    "held_out_pde_residual_evaluated": False,
    "parent_pde_receipt_transferred": False,
    "public_image_numeric_target_used": False,
    "visual_acceptance_threshold_defined": False,
    "visualization_ready": False,
    "visual_correspondence_verified": False,
    "pde_validated": False,
    "source_correspondence_verified": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


def _swirl_column(report: dict) -> np.ndarray:
    """Centered derivative for normalized swirl coordinate a in [-1,1]."""
    o = report["oriented_incremental_desirabilities"]
    plus = np.asarray(
        [
            o["plus"]["aspect_gain_increment"],
            o["plus"]["tip_thinning_increment"],
            o["plus"]["turns_fidelity_increment"],
            o["plus"]["axial_pair_fidelity_increment"],
        ],
        dtype=float,
    )
    minus = np.asarray(
        [
            o["minus"]["aspect_gain_increment"],
            o["minus"]["tip_thinning_increment"],
            o["minus"]["turns_fidelity_increment"],
            o["minus"]["axial_pair_fidelity_increment"],
        ],
        dtype=float,
    )
    return 0.5 * (plus - minus)


def _shoulder_column(report: dict) -> np.ndarray:
    """Derivative for normalized shoulder coordinate b=lambda/0.02 in [-1,1]."""
    o = report["oriented_desirability_sensitivities"]
    per_lambda = np.asarray(
        [
            o["aspect_desirability"],
            o["tip_thinning_desirability"],
            o["path_turns_desirability"],
            o["path_pair_axial_desirability"],
        ],
        dtype=float,
    )
    return shoulder.EPSILON * per_lambda


def _response_diagnostic(jacobian: np.ndarray) -> dict:
    j = np.asarray(jacobian, dtype=float)
    if j.shape != (4, 2):
        raise ValueError("jacobian must have shape (4,2)")
    norms = np.linalg.norm(j, axis=0)
    if np.any(norms <= np.finfo(float).tiny):
        return {
            "rank": int(np.linalg.matrix_rank(j, tol=1e-12)),
            "column_norms": norms.tolist(),
            "normalized_singular_values": [0.0, 0.0],
            "normalized_condition_number": float("inf"),
            "normalized_column_cosine": float("nan"),
        }
    normalized = j / norms[None, :]
    singular = np.linalg.svd(normalized, compute_uv=False)
    return {
        "rank": int(np.linalg.matrix_rank(normalized, tol=1e-12)),
        "column_norms": norms.tolist(),
        "normalized_singular_values": singular.tolist(),
        "normalized_condition_number": float(singular[0] / singular[-1]),
        "normalized_column_cosine": float(np.dot(normalized[:, 0], normalized[:, 1])),
    }


def solve_pareto_cone(jacobian: np.ndarray) -> dict:
    """Exact LP feasibility check for a bounded first-order Pareto direction.

    For each desirability row k, maximize (Jx)_k subject to all rows satisfying
    Jx >= -PARETO_TOL and x in [-1,1]^2.  If any optimum improves its row by
    more than PARETO_TOL, a nonzero local Pareto witness exists.
    """
    j = np.asarray(jacobian, dtype=float)
    if j.shape != (4, 2):
        raise ValueError("jacobian must have shape (4,2)")
    attempts = []
    witnesses = []
    for k, metric in enumerate(METRICS):
        result = linprog(
            c=-j[k],
            A_ub=-j,
            b_ub=np.full(4, PARETO_TOL, dtype=float),
            bounds=[(-1.0, 1.0), (-1.0, 1.0)],
            method="highs",
        )
        item = {
            "metric": metric,
            "solver_success": bool(result.success),
            "solver_status": int(result.status),
            "solver_message": str(result.message),
        }
        if result.success:
            x = np.asarray(result.x, dtype=float)
            d = j @ x
            improvement = float(d[k])
            passes = bool(np.all(d >= -PARETO_TOL - 1e-12) and improvement > PARETO_TOL)
            item.update(
                {
                    "witness": {"swirl_a": float(x[0]), "shoulder_b": float(x[1])},
                    "predicted_desirability_delta": {
                        name: float(value) for name, value in zip(METRICS, d)
                    },
                    "optimized_metric_improvement": improvement,
                    "passes_pareto_rule": passes,
                }
            )
            if passes:
                witnesses.append(item)
        attempts.append(item)
    if witnesses:
        best = max(witnesses, key=lambda row: row["optimized_metric_improvement"])
        direction = best["witness"]
    else:
        direction = None
    return {
        "tolerance": PARETO_TOL,
        "coordinate_bounds": {"swirl_a": [-1.0, 1.0], "shoulder_b": [-1.0, 1.0]},
        "parameter_grid_scan_performed": False,
        "lp_attempts": attempts,
        "target_free_local_pareto_direction": direction,
    }


def run(out: Path) -> dict:
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="agent7-pareto-cone-") as td:
        td_path = Path(td)
        swirl_report = swirl.run(td_path / "swirl.json")
        shoulder_report = shoulder.run(td_path / "shoulder.json")

    if swirl_report["prereg_issue"] != 620:
        raise RuntimeError("unexpected compact-swirl source receipt")
    if shoulder_report["prereg_issue"] != 633:
        raise RuntimeError("unexpected shoulder-timing source receipt")
    if swirl_report["compact_swirl_basis_growth_justified"]:
        raise RuntimeError("upstream #621 verdict drifted")
    if shoulder_report["shoulder_relative_temporal_growth_justified"]:
        raise RuntimeError("upstream #634 verdict drifted")

    swirl_col = _swirl_column(swirl_report)
    shoulder_col = _shoulder_column(shoulder_report)
    jacobian = np.column_stack((swirl_col, shoulder_col))
    response = _response_diagnostic(jacobian)
    pareto = solve_pareto_cone(jacobian)

    response_independent = bool(
        response["rank"] == 2
        and np.isfinite(response["normalized_condition_number"])
        and response["normalized_condition_number"] <= CONDITION_MAX
    )
    combination_justified = bool(
        response_independent and pareto["target_free_local_pareto_direction"] is not None
    )

    report = {
        "task_id": TASK_ID,
        "prereg_issue": PREREG_ISSUE,
        "sources": {
            "linear_temporal": {"pr": SOURCE_TEMPORAL_PR, "head": SOURCE_TEMPORAL_HEAD},
            "compact_swirl": {"pr": SOURCE_SWIRL_PR, "head": SOURCE_SWIRL_HEAD},
            "shoulder_temporal": {"pr": SOURCE_SHOULDER_PR, "head": SOURCE_SHOULDER_HEAD},
            "fixed_render": {"pr": SOURCE_RENDER_PR, "head": SOURCE_RENDER_HEAD},
        },
        "frozen_coordinates": {
            "swirl_a": "a=+/-1 is exact +/-epsilon from PR #621",
            "shoulder_b": "b=+/-1 is lambda=+/-0.02 from PR #634",
        },
        "metrics": list(METRICS),
        "jacobian": {
            "rows": list(METRICS),
            "columns": ["compact_swirl_a", "shoulder_timing_b"],
            "values": jacobian.tolist(),
            "compact_swirl_column": swirl_col.tolist(),
            "shoulder_timing_column": shoulder_col.tolist(),
        },
        "response_diagnostic": response,
        "pareto_cone": pareto,
        "basis_growth_combination_justified": combination_justified,
        "routing": (
            "evaluate exactly one combined child at the LP witness before any basis promotion"
            if combination_justified
            else "keep the minimal #587 family; do not combine these two screened degrees"
        ),
        "upstream_truth": {
            "compact_swirl_growth_justified": swirl_report["compact_swirl_basis_growth_justified"],
            "shoulder_relative_temporal_growth_justified": shoulder_report[
                "shoulder_relative_temporal_growth_justified"
            ],
        },
        **TRUTH,
    }
    out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    report = run(args.out)
    print(
        json.dumps(
            {
                "basis_growth_combination_justified": report[
                    "basis_growth_combination_justified"
                ],
                "response_diagnostic": report["response_diagnostic"],
                "target_free_local_pareto_direction": report["pareto_cone"][
                    "target_free_local_pareto_direction"
                ],
                "jacobian": report["jacobian"],
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
