"""Restricted-force cross-check for the supported Phi(1,0) temporal visualization trial.

This increment deliberately combines two already-owned pieces of work instead of
introducing another PDE operator or another forcing family:

* the slope-1.4 support-connected affine ``Phi(1,0)`` visualization trial selected
  by the target-free nonlinear morphology replay; and
* the existing support-connected vorticity-equation / bounded
  ``RestrictedForce(a,c)`` projection machinery.

The slope is frozen before any PDE calculation.  Only the preregistered force
coefficients ``0 <= a,c <= 10`` are fitted on the existing stratified training
probes and then frozen on the disjoint supported-child holdout probes.  The static
supported baseline is projected with the existing production helper, while the
new temporal field reuses exactly the same probe sets, force columns, finite-
difference operator, bounds and solver settings.

This report is a trade-off diagnostic for candidate routing.  It is not a full
momentum validation, does not refit pressure, does not fit OpenAI imagery and can
never promote visual similarity or solver convergence to PDE validity.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Iterable

import numpy as np
from scipy.optimize import lsq_linear

from .constrained_eq45_supported_force_projection import (
    FIT_MAX_ITER,
    FIT_SPATIAL_STEP,
    FIT_TOL,
    FIXED_TIME_STEP,
    FORCE_BOUNDS,
    _assert_disjoint,
    _curl_force,
    _fit_probe_set,
    _force_columns,
    _region_metrics,
    _vector_metrics,
    project_supported_restricted_force,
)
from .constrained_eq45_supported_phi10_temporal_replay import build_phi10_temporal_trial
from .constrained_eq45_supported_vorticity import _probe_set, vorticity_equation_terms
from .constrained_force import RestrictedForce

FROZEN_SLOPE = 1.4


def _project_temporal_trial(
    field,
    *,
    nu: float,
    derivative_steps: Iterable[float],
    fit_spatial_step: float = FIT_SPATIAL_STEP,
    time_step: float = FIXED_TIME_STEP,
) -> dict:
    """Project the same preregistered force family on one frozen temporal field.

    This is intentionally a thin comparison harness.  The derivative operator,
    probe sets, force-curl columns, bounds and metric helpers are imported from the
    already checked supported-child projection lane.
    """
    if not hasattr(field, "at_points") or not hasattr(field, "sha256"):
        raise TypeError("field must expose public at_points and sha256")
    if not np.isfinite(nu) or nu < 0.0:
        raise ValueError("nu must be finite and nonnegative")
    fit_spatial_step = float(fit_spatial_step)
    time_step = float(time_step)
    if not np.isfinite(fit_spatial_step) or fit_spatial_step <= 0.0:
        raise ValueError("fit_spatial_step must be positive and finite")
    if not np.isfinite(time_step) or time_step <= 0.0:
        raise ValueError("time_step must be positive and finite")
    steps = tuple(float(step) for step in derivative_steps)
    if len(steps) < 3 or any((not np.isfinite(step)) or step <= 0.0 for step in steps):
        raise ValueError("derivative_steps must contain >=3 positive finite levels")
    if any(right >= left for left, right in zip(steps, steps[1:])):
        raise ValueError("derivative_steps must be strictly decreasing")

    fit_points, fit_times, fit_labels = _fit_probe_set()
    holdout_points, holdout_times, holdout_labels = _probe_set()
    _assert_disjoint(fit_points, fit_times, holdout_points, holdout_times)

    fit_terms = vorticity_equation_terms(
        field,
        fit_points,
        fit_times,
        spatial_step=fit_spatial_step,
        time_step=time_step,
        nu=nu,
    )
    fit_target = np.asarray(fit_terms["residual"], dtype=float)
    column_a, column_c = _force_columns(fit_points, fit_times, fit_spatial_step)
    design = np.column_stack((column_a.reshape(-1), column_c.reshape(-1)))
    target = fit_target.reshape(-1)
    singular_values = np.linalg.svd(design, compute_uv=False)
    rank = int(np.linalg.matrix_rank(design))
    condition = float(
        np.inf if singular_values[-1] == 0.0 else singular_values[0] / singular_values[-1]
    )

    solution = lsq_linear(
        design,
        target,
        bounds=FORCE_BOUNDS,
        tol=FIT_TOL,
        lsmr_tol="auto",
        max_iter=FIT_MAX_ITER,
    )
    if not solution.success or not np.all(np.isfinite(solution.x)):
        raise RuntimeError(f"temporal-trial restricted-force solve failed: {solution.message}")
    a, c = float(solution.x[0]), float(solution.x[1])
    fitted_force = RestrictedForce(a=a, c=c)
    fit_after = fit_target - (a * column_a + c * column_c)
    fit_before_metrics = _vector_metrics(fit_target)
    fit_after_metrics = _vector_metrics(fit_after)

    fit_record = {
        "a": a,
        "c": c,
        "bounds": list(FORCE_BOUNDS),
        "spatial_step": fit_spatial_step,
        "time_step": time_step,
        "probe_count": int(len(fit_points)),
        "probe_contract": "fixed_stratified_4x4_plateau_radial_axial_corner",
        "design_rank": rank,
        "design_condition": condition,
        "singular_values": [float(value) for value in singular_values],
        "solver_status": int(solution.status),
        "solver_iterations": int(solution.nit),
        "max_iterations": FIT_MAX_ITER,
        "solver_tolerance": FIT_TOL,
        "before": fit_before_metrics,
        "after": fit_after_metrics,
        "rms_reduction_fraction": float(
            1.0 - fit_after_metrics["rms"] / max(fit_before_metrics["rms"], np.finfo(float).tiny)
        ),
        "by_region": _region_metrics(fit_target, fit_after, fit_labels),
    }

    holdout_rows = []
    for step in steps:
        terms = vorticity_equation_terms(
            field,
            holdout_points,
            holdout_times,
            spatial_step=step,
            time_step=time_step,
            nu=nu,
        )
        before = np.asarray(terms["residual"], dtype=float)
        force_curl = _curl_force(fitted_force, holdout_points, holdout_times, step)
        after = before - force_curl
        before_metrics = _vector_metrics(before)
        after_metrics = _vector_metrics(after)
        holdout_rows.append(
            {
                "spatial_step": step,
                "time_step": time_step,
                "before": before_metrics,
                "after": after_metrics,
                "rms_reduction_fraction": float(
                    1.0 - after_metrics["rms"] / max(before_metrics["rms"], np.finfo(float).tiny)
                ),
                "by_region": _region_metrics(before, after, holdout_labels),
            }
        )

    return {
        "candidate_sha256": field.sha256,
        "fit": fit_record,
        "holdout_rows": holdout_rows,
        "holdout_probe_count": int(len(holdout_points)),
    }


def compare_supported_phi10_temporal_force(
    *,
    slope: float = FROZEN_SLOPE,
    constraints_path: str | Path | None = None,
) -> dict:
    """Compare static and slope-frozen supported candidates under the same force family."""
    slope = float(slope)
    if not np.isfinite(slope) or slope <= 0.0:
        raise ValueError("slope must be positive and finite")
    if not np.isclose(slope, FROZEN_SLOPE, rtol=0.0, atol=0.0):
        raise ValueError("this cross-check is frozen to the morphology-selected slope 1.4")

    root = Path(__file__).resolve().parents[2]
    constraints_path = Path(constraints_path or root / "configs/constraints.json")
    constraints = json.loads(constraints_path.read_text(encoding="utf-8"))
    if constraints["forcing"]["mode"] != "restricted_two_parameter_family":
        raise AssertionError("active forcing convention drifted from the preregistered family")
    nu = float(constraints["nu"])
    steps = tuple(float(step) for step in constraints["validation"]["derivative_steps"])

    trial = build_phi10_temporal_trial(slope=slope)
    static = trial.base
    static_report = project_supported_restricted_force(
        static,
        nu=nu,
        derivative_steps=steps,
    )
    temporal_report = _project_temporal_trial(
        trial,
        nu=nu,
        derivative_steps=steps,
    )

    static_finest = static_report["holdout_rows"][-1]
    temporal_finest = temporal_report["holdout_rows"][-1]
    static_zero = float(static_finest["before"]["rms"])
    temporal_zero = float(temporal_finest["before"]["rms"])
    static_forced = float(static_finest["after"]["rms"])
    temporal_forced = float(temporal_finest["after"]["rms"])

    midpoint = 0.5 * (trial.time_start + trial.time_end)
    midpoint_points = np.asarray(
        [[0.37, 0.11, -0.22], [1.25, 0.0, 0.35], [1.72, 0.2, 0.4]],
        dtype=float,
    )
    midpoint_exact = bool(
        np.array_equal(
            trial.at_points(midpoint_points, midpoint),
            static.at_points(midpoint_points, midpoint),
        )
    )
    if not midpoint_exact:
        raise AssertionError("temporal trial must recover the static supported midpoint exactly")

    return {
        "schema": "eq45_supported_phi10_temporal_force_crosscheck_v1",
        "claim_scope": "bounded_restricted_force_tradeoff_for_frozen_visualization_trial",
        "slope_frozen_before_pde_fit": True,
        "slope": slope,
        "mode": {"family": "phi", "index": [1, 0]},
        "static_supported_sha256": static.sha256,
        "temporal_trial_sha256": trial.sha256,
        "midpoint_exact": midpoint_exact,
        "nu": nu,
        "derivative_steps": list(steps),
        "force_family": "preregistered_restricted_two_parameter_family",
        "force_bounds": list(FORCE_BOUNDS),
        "pressure_fitted": False,
        "new_force_direction_added": False,
        "residual_defined_force_allowed": False,
        "static_projection": {
            "fit": static_report["fit"],
            "holdout_rows": static_report["holdout_rows"],
        },
        "temporal_projection": temporal_report,
        "finest_holdout_comparison": {
            "static_zero_force_rms": static_zero,
            "temporal_zero_force_rms": temporal_zero,
            "zero_force_fractional_change": (temporal_zero - static_zero)
            / max(static_zero, np.finfo(float).tiny),
            "static_projected_force_rms": static_forced,
            "temporal_projected_force_rms": temporal_forced,
            "projected_force_fractional_change": (temporal_forced - static_forced)
            / max(static_forced, np.finfo(float).tiny),
            "static_force_reduction_fraction": float(static_finest["rms_reduction_fraction"]),
            "temporal_force_reduction_fraction": float(temporal_finest["rms_reduction_fraction"]),
        },
        "formal_pde_gate_assessed": False,
        "reason_formal_pde_gate_unassessed": (
            "fixed-probe pressure-free vorticity RMS is not the preregistered volume-weighted "
            "momentum L2/global maximum and pressure is not bound in this report"
        ),
        "visual_morphology_source": "PR143_target_free_slope_1p4_trial",
        "visual_correspondence_verified": False,
        "visualization_ready": False,
        "pde_validated": False,
        "paper_exact": False,
        "openai_field_identified": False,
        "blowup_proved": False,
    }


def main() -> None:
    print(json.dumps(compare_supported_phi10_temporal_force(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
