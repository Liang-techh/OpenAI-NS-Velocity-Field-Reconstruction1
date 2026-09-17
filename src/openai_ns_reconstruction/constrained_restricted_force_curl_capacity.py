"""Bounded restricted-force capacity for the pressure-curl obstruction.

For a frozen velocity define the pressure-free momentum term

    M = u_t + (u . grad)u - nu * laplacian(u).

A scalar pressure can only add a curl-free gradient, so compatibility of the
preregistered forced equation requires

    curl(M - f(a,c)) = 0.

This module fits only the two coefficients of the already-preregistered
``RestrictedForce(a,c)`` family, with the original hard bounds ``0 <= a,c <= 10``.
It never constructs a residual-dependent forcing basis and never changes the
velocity.  The bounded fit is a training-side capacity diagnostic; frozen
held-out measurements remain separate and are not PDE validation.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any, Iterable

import numpy as np
from scipy.optimize import lsq_linear

from .constrained_force import RestrictedForce
from .constrained_pressure_compatibility import curl_pressure_free_momentum

Array = np.ndarray
VelocityLike = Any
REGISTERED_NU = 0.01
FORCE_BOUNDS = (0.0, 10.0)
DEFAULT_FIT_SEED = 20260916
DEFAULT_HOLDOUT_SEED = 914027
DEFAULT_SAMPLE_COUNT = 16
DEFAULT_PROBE_HALF_WIDTH = 0.25
DEFAULT_PROBE_TIME_RANGE = (0.38, 0.62)
DEFAULT_FIT_STEP = 0.005
DEFAULT_DERIVATIVE_STEPS = (0.02, 0.01, 0.005)


@dataclass(frozen=True)
class RestrictedForceCurlFit:
    a: float
    c: float
    fit_step: float
    design_rank: int
    design_condition: float
    singular_values: tuple[float, ...]
    active_mask: tuple[int, int]
    solver_status: int
    solver_iterations: int
    nonlinear_function_evaluations: int
    curl_rms_before: float
    curl_max_before: float
    curl_rms_after: float
    curl_max_after: float
    recoverable_rms_fraction: float


@dataclass(frozen=True)
class RestrictedForceCurlLevel:
    spatial_step: float
    time_step: float
    curl_rms_before: float
    curl_max_before: float
    curl_rms_after: float
    curl_max_after: float
    recoverable_rms_fraction: float


@dataclass(frozen=True)
class RestrictedForceCurlCapacityReport:
    nu: float
    force_family: str
    force_bounds: tuple[float, float]
    pressure_role: str
    fit_seed: int | None
    holdout_seed: int | None
    fit_point_count: int
    holdout_point_count: int
    fit: RestrictedForceCurlFit
    holdout_levels: tuple[RestrictedForceCurlLevel, ...]
    fit_and_holdout_separate: bool = True
    velocity_changed: bool = False
    pressure_fitted: bool = False
    forcing_fitted: bool = True
    pde_validated: bool = False

    def to_dict(self) -> dict[str, Any]:
        payload = asdict(self)
        payload["fit"] = asdict(self.fit)
        payload["holdout_levels"] = [asdict(level) for level in self.holdout_levels]
        payload["interpretation"] = (
            "bounded same-family force fit plus frozen held-out curl(M-f) audit; "
            "residual reduction is capacity evidence only, not PDE validation"
        )
        return payload


def deterministic_probe_cloud(
    seed: int,
    count: int = DEFAULT_SAMPLE_COUNT,
    *,
    half_width: float = DEFAULT_PROBE_HALF_WIDTH,
    time_range: tuple[float, float] = DEFAULT_PROBE_TIME_RANGE,
) -> tuple[Array, Array]:
    """Return deterministic interior probes for the fit/holdout diagnostic."""
    if not isinstance(seed, (int, np.integer)):
        raise ValueError("seed must be an integer")
    if not isinstance(count, (int, np.integer)) or count <= 0:
        raise ValueError("count must be a positive integer")
    half_width = float(half_width)
    if not np.isfinite(half_width) or half_width <= 0.0:
        raise ValueError("half_width must be positive and finite")
    if len(time_range) != 2:
        raise ValueError("time_range must contain exactly two values")
    t0, t1 = (float(time_range[0]), float(time_range[1]))
    if not np.isfinite(t0) or not np.isfinite(t1) or not t0 < t1:
        raise ValueError("time_range must be finite and strictly increasing")

    rng = np.random.default_rng(int(seed))
    points = rng.uniform(-half_width, half_width, size=(int(count), 3))
    times = rng.uniform(t0, t1, size=int(count))
    return points, times


def _as_samples(points: Array, times: Array, name: str) -> tuple[Array, Array]:
    pts = np.asarray(points, dtype=float)
    ts = np.asarray(times, dtype=float)
    if pts.ndim != 2 or pts.shape[1] != 3 or pts.shape[0] == 0:
        raise ValueError(f"{name}_points must have shape (N, 3) with N > 0")
    try:
        ts = np.broadcast_to(ts, (pts.shape[0],)).astype(float, copy=False)
    except ValueError as exc:
        raise ValueError(f"{name}_times must broadcast to ({pts.shape[0]},)") from exc
    if not np.all(np.isfinite(pts)) or not np.all(np.isfinite(ts)):
        raise ValueError(f"{name} samples must be finite")
    return pts, ts


def _vector_rms(values: Array) -> float:
    return float(np.sqrt(np.mean(np.sum(values * values, axis=-1))))


def _vector_max(values: Array) -> float:
    return float(np.max(np.linalg.norm(values, axis=-1)))


def _curl_force(force: RestrictedForce, points: Array, times: Array, step: float) -> Array:
    """Centered Cartesian curl of the preregistered force on supplied probes."""
    h = float(step)
    if not np.isfinite(h) or h <= 0.0:
        raise ValueError("force curl step must be positive and finite")
    derivatives: list[Array] = []
    for axis in range(3):
        shift = np.zeros(3)
        shift[axis] = h
        plus = np.asarray(force(points + shift, times), dtype=float)
        minus = np.asarray(force(points - shift, times), dtype=float)
        if plus.shape != points.shape or minus.shape != points.shape:
            raise ValueError("restricted force returned malformed values")
        if not np.all(np.isfinite(plus)) or not np.all(np.isfinite(minus)):
            raise ValueError("restricted force returned nonfinite values")
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


def _force_curl_columns(points: Array, times: Array, step: float) -> tuple[Array, Array]:
    column_a = _curl_force(RestrictedForce(a=1.0, c=0.0), points, times, step)
    column_c = _curl_force(RestrictedForce(a=0.0, c=1.0), points, times, step)
    return column_a, column_c


def audit_restricted_force_curl_capacity(
    velocity: VelocityLike,
    fit_points: Array,
    fit_times: Array,
    holdout_points: Array,
    holdout_times: Array,
    *,
    fit_seed: int | None = None,
    holdout_seed: int | None = None,
    nu: float = REGISTERED_NU,
    fit_step: float = DEFAULT_FIT_STEP,
    derivative_steps: Iterable[float] = DEFAULT_DERIVATIVE_STEPS,
) -> RestrictedForceCurlCapacityReport:
    """Fit only ``a,c`` on training probes and freeze them on held-out probes."""
    if not np.isfinite(nu) or not np.isclose(float(nu), REGISTERED_NU, rtol=0.0, atol=0.0):
        raise ValueError(f"nu is preregistered and must remain {REGISTERED_NU}")
    fit_pts, fit_ts = _as_samples(fit_points, fit_times, "fit")
    hold_pts, hold_ts = _as_samples(holdout_points, holdout_times, "holdout")
    if fit_pts.shape == hold_pts.shape and np.array_equal(fit_pts, hold_pts) and np.array_equal(fit_ts, hold_ts):
        raise ValueError("fit and holdout samples must be disjoint")
    if fit_seed is not None and holdout_seed is not None and int(fit_seed) == int(holdout_seed):
        raise ValueError("fit_seed and holdout_seed must be distinct")

    fit_step = float(fit_step)
    if not np.isfinite(fit_step) or fit_step <= 0.0:
        raise ValueError("fit_step must be positive and finite")
    steps = tuple(float(step) for step in derivative_steps)
    if len(steps) < 3 or any(not np.isfinite(step) or step <= 0.0 for step in steps):
        raise ValueError("derivative_steps must contain at least three positive finite levels")
    if any(next_step >= step for step, next_step in zip(steps, steps[1:])):
        raise ValueError("derivative_steps must be strictly decreasing")

    _, fit_target = curl_pressure_free_momentum(
        velocity,
        fit_pts,
        fit_ts,
        nu=REGISTERED_NU,
        spatial_step=fit_step,
        time_step=fit_step,
    )
    column_a, column_c = _force_curl_columns(fit_pts, fit_ts, fit_step)
    design = np.column_stack((column_a.reshape(-1), column_c.reshape(-1)))
    target = fit_target.reshape(-1)
    singular_values = np.linalg.svd(design, compute_uv=False)
    rank = int(np.linalg.matrix_rank(design))
    condition = float(np.inf if singular_values[-1] == 0.0 else singular_values[0] / singular_values[-1])

    solution = lsq_linear(
        design,
        target,
        bounds=FORCE_BOUNDS,
        tol=1e-12,
        lsmr_tol="auto",
        max_iter=200,
    )
    if not solution.success or not np.all(np.isfinite(solution.x)):
        raise RuntimeError(f"bounded restricted-force solve failed: {solution.message}")
    a, c = (float(solution.x[0]), float(solution.x[1]))
    fitted_force = RestrictedForce(a=a, c=c)
    fit_force_curl = a * column_a + c * column_c
    fit_after = fit_target - fit_force_curl
    fit_before_rms = _vector_rms(fit_target)
    fit_after_rms = _vector_rms(fit_after)

    fit = RestrictedForceCurlFit(
        a=a,
        c=c,
        fit_step=fit_step,
        design_rank=rank,
        design_condition=condition,
        singular_values=tuple(float(value) for value in singular_values),
        active_mask=tuple(int(value) for value in solution.active_mask),
        solver_status=int(solution.status),
        solver_iterations=int(solution.nit),
        nonlinear_function_evaluations=0,
        curl_rms_before=fit_before_rms,
        curl_max_before=_vector_max(fit_target),
        curl_rms_after=fit_after_rms,
        curl_max_after=_vector_max(fit_after),
        recoverable_rms_fraction=float(1.0 - fit_after_rms / max(fit_before_rms, 1e-15)),
    )

    levels: list[RestrictedForceCurlLevel] = []
    for step in steps:
        _, hold_target = curl_pressure_free_momentum(
            velocity,
            hold_pts,
            hold_ts,
            nu=REGISTERED_NU,
            spatial_step=step,
            time_step=step,
        )
        hold_force_curl = _curl_force(fitted_force, hold_pts, hold_ts, step)
        hold_after = hold_target - hold_force_curl
        before_rms = _vector_rms(hold_target)
        after_rms = _vector_rms(hold_after)
        levels.append(
            RestrictedForceCurlLevel(
                spatial_step=step,
                time_step=step,
                curl_rms_before=before_rms,
                curl_max_before=_vector_max(hold_target),
                curl_rms_after=after_rms,
                curl_max_after=_vector_max(hold_after),
                recoverable_rms_fraction=float(1.0 - after_rms / max(before_rms, 1e-15)),
            )
        )

    return RestrictedForceCurlCapacityReport(
        nu=REGISTERED_NU,
        force_family="preregistered_restricted_two_parameter_family",
        force_bounds=FORCE_BOUNDS,
        pressure_role="scalar_gradient_cannot_change_curl_obstruction",
        fit_seed=None if fit_seed is None else int(fit_seed),
        holdout_seed=None if holdout_seed is None else int(holdout_seed),
        fit_point_count=fit_pts.shape[0],
        holdout_point_count=hold_pts.shape[0],
        fit=fit,
        holdout_levels=tuple(levels),
    )


def audit_default_eq45_restricted_force_capacity(velocity: VelocityLike) -> RestrictedForceCurlCapacityReport:
    """Run the checked deterministic Eq45 diagnostic probe contract."""
    fit_points, fit_times = deterministic_probe_cloud(DEFAULT_FIT_SEED)
    holdout_points, holdout_times = deterministic_probe_cloud(DEFAULT_HOLDOUT_SEED)
    return audit_restricted_force_curl_capacity(
        velocity,
        fit_points,
        fit_times,
        holdout_points,
        holdout_times,
        fit_seed=DEFAULT_FIT_SEED,
        holdout_seed=DEFAULT_HOLDOUT_SEED,
    )
