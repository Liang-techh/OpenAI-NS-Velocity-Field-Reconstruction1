"""Bounded restricted-force projection for the support-connected Eq45 candidate.

This module is intentionally stacked on the independent supported-child vorticity
operator in :mod:`constrained_eq45_supported_vorticity`.  It changes neither the
velocity nor the support transform.  It asks one narrow CR004 question: after the
physical-support composition changes ``[u,v,w]`` in the collar, how much of the
pressure-free vorticity-equation defect can the *already preregistered* two-parameter
``RestrictedForce(a,c)`` family remove?

Only ``0 <= a,c <= 10`` are fitted, by bounded linear variable projection on a fixed,
stratified training probe set.  The coefficients are then frozen on the disjoint
supported-child holdout set and the registered spatial derivative ladder.  No
residual-defined force direction, pressure basis, velocity coefficient, taper
parameter, or acceptance threshold is introduced here.

The resulting sampled RMS/max values are a necessary-condition diagnostic.  They are
not the preregistered volume-weighted spatial L2/global maximum, so this module never
promotes ``pde_validated``.
"""
from __future__ import annotations

import json
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Iterable

import numpy as np
from scipy.optimize import lsq_linear

from .constrained_eq45_candidate import Eq45VelocityCandidate
from .constrained_eq45_supported_candidate import Eq45SupportedVelocityCandidate
from .constrained_eq45_supported_vorticity import _probe_set, vorticity_equation_terms
from .constrained_force import RestrictedForce

FORCE_BOUNDS = (0.0, 10.0)
FIT_SPATIAL_STEP = 0.005
FIXED_TIME_STEP = 0.0025
FIT_MAX_ITER = 200
FIT_TOL = 1.0e-12


def _fit_probe_set() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Fixed probes balanced across the untouched plateau and three taper collars."""
    points = np.asarray(
        [
            [0.31, -0.22, 0.20],
            [-0.37, -0.18, -0.27],
            [0.47, 0.29, 0.35],
            [-0.33, 0.41, -0.38],
            [1.72, -0.18, 0.20],
            [-1.68, -0.35, -0.30],
            [1.28, -1.22, 0.42],
            [-1.35, 1.10, -0.48],
            [0.38, -0.12, 1.70],
            [-0.32, -0.28, -1.68],
            [0.58, 0.23, 1.78],
            [-0.22, 0.50, -1.76],
            [1.42, 0.82, 1.70],
            [-1.45, 0.70, -1.72],
            [1.10, -1.30, 1.76],
            [-1.18, -1.22, -1.74],
        ],
        dtype=float,
    )
    labels = np.asarray(
        ["plateau"] * 4
        + ["radial_collar"] * 4
        + ["axial_collar"] * 4
        + ["corner_collar"] * 4,
        dtype=object,
    )
    times = np.asarray([0.37, 0.47, 0.57, 0.63] * 4, dtype=float)
    return points, times, labels


def _vector_metrics(values: np.ndarray) -> dict[str, float]:
    value = np.asarray(values, dtype=float)
    if value.ndim != 2 or value.shape[1] != 3 or len(value) == 0:
        raise ValueError("residual values must have shape (n,3) with n>0")
    if not np.all(np.isfinite(value)):
        raise FloatingPointError("residual values must be finite")
    norms = np.linalg.norm(value, axis=1)
    return {
        "max": float(np.max(norms)),
        "rms": float(np.sqrt(np.mean(norms * norms))),
    }


def _curl_force(
    force: RestrictedForce,
    points: np.ndarray,
    times: np.ndarray,
    step: float,
) -> np.ndarray:
    """Centered Cartesian curl of the preregistered force family."""
    h = float(step)
    if not np.isfinite(h) or h <= 0.0:
        raise ValueError("force-curl step must be positive and finite")
    points = np.asarray(points, dtype=float)
    times = np.asarray(times, dtype=float)
    if points.ndim != 2 or points.shape[1] != 3 or len(points) == 0:
        raise ValueError("points must have shape (n,3) with n>0")
    times = np.broadcast_to(times, (len(points),)).astype(float, copy=False)
    if not np.all(np.isfinite(points)) or not np.all(np.isfinite(times)):
        raise ValueError("force-curl inputs must be finite")

    derivatives: list[np.ndarray] = []
    for axis in range(3):
        shift = np.zeros_like(points)
        shift[:, axis] = h
        plus = np.asarray(force(points + shift, times), dtype=float)
        minus = np.asarray(force(points - shift, times), dtype=float)
        if plus.shape != points.shape or minus.shape != points.shape:
            raise ValueError("RestrictedForce returned malformed values")
        if not np.all(np.isfinite(plus)) or not np.all(np.isfinite(minus)):
            raise FloatingPointError("RestrictedForce returned nonfinite values")
        derivatives.append((plus - minus) / (2.0 * h))

    d_dx, d_dy, d_dz = derivatives
    return np.stack(
        (
            d_dy[:, 2] - d_dz[:, 1],
            d_dz[:, 0] - d_dx[:, 2],
            d_dx[:, 1] - d_dy[:, 0],
        ),
        axis=-1,
    )


def _force_columns(
    points: np.ndarray,
    times: np.ndarray,
    step: float,
) -> tuple[np.ndarray, np.ndarray]:
    return (
        _curl_force(RestrictedForce(a=1.0, c=0.0), points, times, step),
        _curl_force(RestrictedForce(a=0.0, c=1.0), points, times, step),
    )


def _region_metrics(
    before: np.ndarray,
    after: np.ndarray,
    labels: np.ndarray,
) -> dict[str, dict[str, dict[str, float] | float]]:
    result: dict[str, dict[str, dict[str, float] | float]] = {}
    for region in tuple(dict.fromkeys(labels.tolist())):
        mask = labels == region
        before_metrics = _vector_metrics(before[mask])
        after_metrics = _vector_metrics(after[mask])
        result[region] = {
            "before": before_metrics,
            "after": after_metrics,
            "rms_reduction_fraction": float(
                1.0 - after_metrics["rms"] / max(before_metrics["rms"], np.finfo(float).tiny)
            ),
        }
    return result


def _assert_disjoint(
    fit_points: np.ndarray,
    fit_times: np.ndarray,
    holdout_points: np.ndarray,
    holdout_times: np.ndarray,
) -> None:
    for point, time in zip(fit_points, fit_times):
        same_point = np.all(holdout_points == point, axis=1)
        same_time = holdout_times == time
        if np.any(same_point & same_time):
            raise AssertionError("fit and holdout probe tuples must be disjoint")


def project_supported_restricted_force(
    child: Eq45SupportedVelocityCandidate,
    *,
    nu: float,
    derivative_steps: Iterable[float],
    fit_spatial_step: float = FIT_SPATIAL_STEP,
    time_step: float = FIXED_TIME_STEP,
) -> dict:
    """Fit only ``a,c`` on fixed training probes and freeze them on holdout probes."""
    if not isinstance(child, Eq45SupportedVelocityCandidate):
        raise TypeError("child must be an Eq45SupportedVelocityCandidate")
    if not np.isfinite(nu) or nu < 0.0:
        raise ValueError("nu must be finite and nonnegative")
    fit_spatial_step = float(fit_spatial_step)
    time_step = float(time_step)
    if not np.isfinite(fit_spatial_step) or fit_spatial_step <= 0.0:
        raise ValueError("fit_spatial_step must be positive and finite")
    if not np.isfinite(time_step) or time_step <= 0.0:
        raise ValueError("time_step must be positive and finite")
    steps = tuple(float(step) for step in derivative_steps)
    if len(steps) < 3 or any(not np.isfinite(step) or step <= 0.0 for step in steps):
        raise ValueError("derivative_steps must contain >=3 positive finite levels")
    if any(right >= left for left, right in zip(steps, steps[1:])):
        raise ValueError("derivative_steps must be strictly decreasing")

    fit_points, fit_times, fit_labels = _fit_probe_set()
    holdout_points, holdout_times, holdout_labels = _probe_set()
    _assert_disjoint(fit_points, fit_times, holdout_points, holdout_times)

    fit_terms = vorticity_equation_terms(
        child,
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
        raise RuntimeError(f"supported-child restricted-force solve failed: {solution.message}")
    a, c = (float(solution.x[0]), float(solution.x[1]))
    fitted_force = RestrictedForce(a=a, c=c)
    fit_force_curl = a * column_a + c * column_c
    fit_after = fit_target - fit_force_curl
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
            child,
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
                    1.0
                    - after_metrics["rms"] / max(before_metrics["rms"], np.finfo(float).tiny)
                ),
                "by_region": _region_metrics(before, after, holdout_labels),
            }
        )

    return {
        "schema": "eq45_supported_restricted_force_projection_v1",
        "claim_scope": "bounded_same_family_force_projection_on_supported_child_vorticity_obstruction",
        "supported_child_sha256": child.sha256,
        "parent_sha256": child.parent.sha256,
        "nu": float(nu),
        "force_family": "preregistered_restricted_two_parameter_family",
        "force_parameters_fitted": ["a", "c"],
        "force_bounds": list(FORCE_BOUNDS),
        "fit_and_holdout_separate": True,
        "holdout_force_refit": False,
        "untapered_force_transferred": False,
        "pressure_fitted": False,
        "new_force_direction_added": False,
        "residual_defined_force_allowed": False,
        "fit": fit_record,
        "holdout_probe_count": int(len(holdout_points)),
        "holdout_probe_contract": "agent3_supported_vorticity_fixed_16_region_probes",
        "holdout_rows": holdout_rows,
        "formal_pde_gate_assessed": False,
        "reason_formal_pde_gate_unassessed": (
            "fixed-probe vorticity RMS/max are not the preregistered volume-weighted spatial L2/global maximum; "
            "full supported-child momentum pressure binding is also not established here"
        ),
        "velocity_changed": False,
        "support_transform_changed": False,
        "visualization_ready": False,
        "visual_correspondence_verified": False,
        "pde_validated": False,
        "paper_exact": False,
        "openai_field_identified": False,
        "blowup_proved": False,
    }


def audit_supported_eq45_restricted_force(
    *,
    parent_path: str | Path | None = None,
    constraints_path: str | Path | None = None,
) -> dict:
    """Run the governed production supported-child restricted-force projection."""
    root = Path(__file__).resolve().parents[2]
    parent_path = Path(parent_path or root / "artifacts/constrained/eq45_velocity_candidate_seed.json")
    constraints_path = Path(constraints_path or root / "configs/constraints.json")

    parent = Eq45VelocityCandidate.load_json(parent_path)
    constructed_child = Eq45SupportedVelocityCandidate(parent=parent)
    with TemporaryDirectory() as tmp:
        child_path = Path(tmp) / "supported.json"
        constructed_child.save_json(child_path)
        child = Eq45SupportedVelocityCandidate.load_json(child_path)
    if child.sha256 != constructed_child.sha256:
        raise AssertionError("supported-child serialization changed candidate identity")

    constraints = json.loads(constraints_path.read_text(encoding="utf-8"))
    if constraints["forcing"]["mode"] != "restricted_two_parameter_family":
        raise AssertionError("active forcing convention drifted from the preregistered restricted family")
    nu = float(constraints["nu"])
    steps = tuple(float(step) for step in constraints["validation"]["derivative_steps"])
    report = project_supported_restricted_force(
        child,
        nu=nu,
        derivative_steps=steps,
    )
    report["registered_pde_thresholds"] = {
        "max": float(constraints["validation"]["thresholds"]["pde_residual_max"]),
        "L2": float(constraints["validation"]["thresholds"]["pde_residual_L2"]),
    }
    return report


def main() -> None:
    print(json.dumps(audit_supported_eq45_restricted_force(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
