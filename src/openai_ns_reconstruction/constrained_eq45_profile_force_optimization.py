"""Bounded Eq45 profile optimization with restricted-force variable projection.

This training-side diagnostic varies only the existing autonomous ``(0,2)``
coefficients of ``Phi`` and ``F``.  At every nonlinear profile trial it eliminates
only the preregistered two-parameter ``RestrictedForce(a,c)`` subproblem by bounded
linear least squares.  The force basis is fixed before seeing the residual; no
residual-defined forcing direction is allowed.

The returned profile and force parameters are capacity evidence only.  Solver
convergence is not Navier--Stokes validation, and this module never promotes the
candidate to paper-exact, visualization-verified, or PDE-valid status.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, replace
from typing import Any, Iterable

import numpy as np
from scipy.optimize import lsq_linear, minimize

from .constrained_eq45_candidate import Eq45VelocityCandidate
from .constrained_pressure_compatibility import curl_pressure_free_momentum
from .constrained_restricted_force_curl_capacity import (
    FORCE_BOUNDS,
    REGISTERED_NU,
    _curl_force,
    _vector_max,
    _vector_rms,
    deterministic_probe_cloud,
)
from .constrained_force import RestrictedForce

Array = np.ndarray
PROFILE_BOUNDS = (-4.0, 4.0)
DEFAULT_FIT_SEED = 20260917
DEFAULT_HOLDOUT_SEED = 914117
DEFAULT_SAMPLE_COUNT = 16
DEFAULT_FIT_STEP = 0.005
DEFAULT_DERIVATIVE_STEPS = (0.02, 0.01, 0.005)
DEFAULT_MAX_FUNCTION_EVALUATIONS = 80
DEFAULT_FTOL = 1e-9
DEFAULT_GTOL = 1e-5


@dataclass(frozen=True)
class ProfileForceLevel:
    spatial_step: float
    seed_rms_after_force: float
    optimized_rms_before_force: float
    optimized_rms_after_force: float
    optimized_max_after_force: float
    improvement_vs_seed: float


@dataclass(frozen=True)
class Eq45ProfileForceOptimizationReport:
    fit_seed: int
    holdout_seed: int
    sample_count: int
    fit_step: float
    derivative_steps: tuple[float, ...]
    profile_bounds: tuple[float, float]
    force_bounds: tuple[float, float]
    optimizer: str
    max_function_evaluations: int
    optimizer_success: bool
    optimizer_message: str
    optimizer_function_evaluations: int
    optimizer_iterations: int
    phi_02: float
    F_02: float
    force_a: float
    force_c: float
    seed_phi_02: float
    seed_F_02: float
    seed_force_a: float
    seed_force_c: float
    seed_training_rms_after_force: float
    optimized_training_rms_before_force: float
    optimized_training_rms_after_force: float
    optimized_training_max_after_force: float
    holdout_levels: tuple[ProfileForceLevel, ...]
    velocity_seed_rms: float
    velocity_optimized_rms: float
    velocity_delta_rms: float
    velocity_relative_delta_rms: float
    velocity_delta_max: float
    optimized_candidate_sha256: str
    optimized_candidate_roundtrip_equal: bool
    velocity_changed: bool = True
    pressure_fitted: bool = False
    forcing_family_changed: bool = False
    pde_validated: bool = False
    visualization_ready_promoted: bool = False
    paper_exact: bool = False
    openai_field_identified: bool = False

    def to_dict(self) -> dict[str, Any]:
        payload = asdict(self)
        payload["holdout_levels"] = [asdict(level) for level in self.holdout_levels]
        payload["interpretation"] = (
            "bounded existing-profile optimization with preregistered restricted-force "
            "variable projection; training convergence and holdout reduction are capacity "
            "evidence only, not PDE or visualization validation"
        )
        return payload


def candidate_with_profile_pair(
    candidate: Eq45VelocityCandidate,
    phi_02: float,
    F_02: float,
) -> Eq45VelocityCandidate:
    """Return ``candidate`` with only the existing ``(0,2)`` Phi/F modes changed."""
    if not isinstance(candidate, Eq45VelocityCandidate):
        raise TypeError("candidate must be Eq45VelocityCandidate")
    bounds = (-candidate.profile_basis.coefficient_limit, candidate.profile_basis.coefficient_limit)
    values = (float(phi_02), float(F_02))
    if not all(np.isfinite(value) and bounds[0] <= value <= bounds[1] for value in values):
        raise ValueError(f"profile coefficients must lie in [{bounds[0]}, {bounds[1]}]")
    try:
        index = candidate.profile_basis.mode_indices.index((0, 2))
    except ValueError as exc:
        raise ValueError("candidate basis does not contain the required (0,2) mode") from exc

    phi = list(candidate.profile_basis.phi_coefficients)
    swirl = list(candidate.profile_basis.swirl_coefficients)
    phi[index] = values[0]
    swirl[index] = values[1]
    basis = replace(
        candidate.profile_basis,
        phi_coefficients=tuple(phi),
        swirl_coefficients=tuple(swirl),
    )
    return replace(candidate, profile_basis=basis)


def _force_design(points: Array, times: Array, step: float) -> Array:
    column_a = _curl_force(RestrictedForce(a=1.0, c=0.0), points, times, step)
    column_c = _curl_force(RestrictedForce(a=0.0, c=1.0), points, times, step)
    return np.column_stack((column_a.reshape(-1), column_c.reshape(-1)))


def _project_force(target: Array, design: Array):
    solution = lsq_linear(
        design,
        target.reshape(-1),
        bounds=FORCE_BOUNDS,
        tol=1e-12,
        lsmr_tol="auto",
        max_iter=200,
    )
    if not solution.success or not np.all(np.isfinite(solution.x)):
        raise RuntimeError(f"bounded restricted-force solve failed: {solution.message}")
    after = target - (design @ solution.x).reshape(target.shape)
    return solution, after


def optimize_eq45_profile_force(
    candidate: Eq45VelocityCandidate | None = None,
    *,
    fit_seed: int = DEFAULT_FIT_SEED,
    holdout_seed: int = DEFAULT_HOLDOUT_SEED,
    sample_count: int = DEFAULT_SAMPLE_COUNT,
    fit_step: float = DEFAULT_FIT_STEP,
    derivative_steps: Iterable[float] = DEFAULT_DERIVATIVE_STEPS,
    max_function_evaluations: int = DEFAULT_MAX_FUNCTION_EVALUATIONS,
) -> tuple[Eq45VelocityCandidate, Eq45ProfileForceOptimizationReport]:
    """Optimize two existing profile coefficients and freeze them on held-out probes."""
    base = Eq45VelocityCandidate.seed() if candidate is None else candidate
    if not isinstance(base, Eq45VelocityCandidate):
        raise TypeError("candidate must be Eq45VelocityCandidate")
    if fit_seed == holdout_seed:
        raise ValueError("fit and holdout seeds must differ")
    if not isinstance(max_function_evaluations, int) or max_function_evaluations <= 0:
        raise ValueError("max_function_evaluations must be a positive integer")
    steps = tuple(float(step) for step in derivative_steps)
    if len(steps) < 3 or any(not np.isfinite(step) or step <= 0.0 for step in steps):
        raise ValueError("derivative_steps must contain at least three positive finite levels")
    if any(next_step >= step for step, next_step in zip(steps, steps[1:])):
        raise ValueError("derivative_steps must be strictly decreasing")

    fit_points, fit_times = deterministic_probe_cloud(fit_seed, sample_count)
    holdout_points, holdout_times = deterministic_probe_cloud(holdout_seed, sample_count)
    design = _force_design(fit_points, fit_times, fit_step)
    mode_index = base.profile_basis.mode_indices.index((0, 2))
    x0 = np.array(
        [
            base.profile_basis.phi_coefficients[mode_index],
            base.profile_basis.swirl_coefficients[mode_index],
        ],
        dtype=float,
    )

    def objective(parameters: Array) -> float:
        trial = candidate_with_profile_pair(base, parameters[0], parameters[1])
        _, target = curl_pressure_free_momentum(
            trial,
            fit_points,
            fit_times,
            nu=REGISTERED_NU,
            spatial_step=fit_step,
            time_step=fit_step,
        )
        _, after = _project_force(target, design)
        return _vector_rms(after)

    seed_target = curl_pressure_free_momentum(
        base,
        fit_points,
        fit_times,
        nu=REGISTERED_NU,
        spatial_step=fit_step,
        time_step=fit_step,
    )[1]
    seed_force_solution, seed_after = _project_force(seed_target, design)

    nonlinear = minimize(
        objective,
        x0,
        method="L-BFGS-B",
        bounds=(PROFILE_BOUNDS, PROFILE_BOUNDS),
        options={
            "maxfun": max_function_evaluations,
            "maxiter": max_function_evaluations,
            "ftol": DEFAULT_FTOL,
            "gtol": DEFAULT_GTOL,
            "maxls": 20,
        },
    )
    if not np.all(np.isfinite(nonlinear.x)):
        raise RuntimeError("profile optimizer returned nonfinite parameters")

    optimized = candidate_with_profile_pair(base, nonlinear.x[0], nonlinear.x[1])
    optimized_target = curl_pressure_free_momentum(
        optimized,
        fit_points,
        fit_times,
        nu=REGISTERED_NU,
        spatial_step=fit_step,
        time_step=fit_step,
    )[1]
    force_solution, optimized_after = _project_force(optimized_target, design)
    fitted_force = RestrictedForce(a=float(force_solution.x[0]), c=float(force_solution.x[1]))

    levels: list[ProfileForceLevel] = []
    seed_force = RestrictedForce(
        a=float(seed_force_solution.x[0]), c=float(seed_force_solution.x[1])
    )
    for step in steps:
        seed_holdout = curl_pressure_free_momentum(
            base,
            holdout_points,
            holdout_times,
            nu=REGISTERED_NU,
            spatial_step=step,
            time_step=step,
        )[1]
        optimized_holdout = curl_pressure_free_momentum(
            optimized,
            holdout_points,
            holdout_times,
            nu=REGISTERED_NU,
            spatial_step=step,
            time_step=step,
        )[1]
        seed_after_holdout = seed_holdout - _curl_force(
            seed_force, holdout_points, holdout_times, step
        )
        optimized_after_holdout = optimized_holdout - _curl_force(
            fitted_force, holdout_points, holdout_times, step
        )
        seed_rms = _vector_rms(seed_after_holdout)
        optimized_rms = _vector_rms(optimized_after_holdout)
        levels.append(
            ProfileForceLevel(
                spatial_step=step,
                seed_rms_after_force=seed_rms,
                optimized_rms_before_force=_vector_rms(optimized_holdout),
                optimized_rms_after_force=optimized_rms,
                optimized_max_after_force=_vector_max(optimized_after_holdout),
                improvement_vs_seed=float(1.0 - optimized_rms / max(seed_rms, 1e-15)),
            )
        )

    seed_velocity = np.asarray(base.at_points(holdout_points, holdout_times), dtype=float)
    optimized_velocity = np.asarray(
        optimized.at_points(holdout_points, holdout_times), dtype=float
    )
    velocity_delta = optimized_velocity - seed_velocity
    seed_velocity_rms = _vector_rms(seed_velocity)
    roundtrip = Eq45VelocityCandidate.from_dict(optimized.to_dict())

    report = Eq45ProfileForceOptimizationReport(
        fit_seed=int(fit_seed),
        holdout_seed=int(holdout_seed),
        sample_count=int(sample_count),
        fit_step=float(fit_step),
        derivative_steps=steps,
        profile_bounds=PROFILE_BOUNDS,
        force_bounds=FORCE_BOUNDS,
        optimizer="scipy.optimize.minimize(method='L-BFGS-B') with nested lsq_linear variable projection",
        max_function_evaluations=int(max_function_evaluations),
        optimizer_success=bool(nonlinear.success),
        optimizer_message=str(nonlinear.message),
        optimizer_function_evaluations=int(nonlinear.nfev),
        optimizer_iterations=int(nonlinear.nit),
        phi_02=float(nonlinear.x[0]),
        F_02=float(nonlinear.x[1]),
        force_a=float(force_solution.x[0]),
        force_c=float(force_solution.x[1]),
        seed_phi_02=float(x0[0]),
        seed_F_02=float(x0[1]),
        seed_force_a=float(seed_force_solution.x[0]),
        seed_force_c=float(seed_force_solution.x[1]),
        seed_training_rms_after_force=_vector_rms(seed_after),
        optimized_training_rms_before_force=_vector_rms(optimized_target),
        optimized_training_rms_after_force=_vector_rms(optimized_after),
        optimized_training_max_after_force=_vector_max(optimized_after),
        holdout_levels=tuple(levels),
        velocity_seed_rms=seed_velocity_rms,
        velocity_optimized_rms=_vector_rms(optimized_velocity),
        velocity_delta_rms=_vector_rms(velocity_delta),
        velocity_relative_delta_rms=float(
            _vector_rms(velocity_delta) / max(seed_velocity_rms, 1e-15)
        ),
        velocity_delta_max=_vector_max(velocity_delta),
        optimized_candidate_sha256=optimized.sha256,
        optimized_candidate_roundtrip_equal=(roundtrip.sha256 == optimized.sha256),
    )
    return optimized, report
