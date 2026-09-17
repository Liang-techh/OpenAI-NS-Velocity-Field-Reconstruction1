"""Minimal bounded Eq45 profile extension after the (0,2) fit is exhausted.

This training-side diagnostic starts from the checked CR005 Eq45 profile/force
receipt, freezes its fitted ``Phi(0,2)`` and ``F(0,2)`` coefficients, and varies
only the already-existing ``Phi(1,0)`` coefficient.  At every trial it eliminates
only the preregistered two-parameter ``RestrictedForce(a,c)`` subproblem with the
same bounded variable projection used by the parent fit.

No new basis direction, pressure model, forcing direction, validation threshold,
or canonical candidate is introduced.  The resulting field is a candidate for
subsequent geometry/visual-fingerprint review; optimizer convergence and residual
reduction are not Navier--Stokes validation or OpenAI-field identification.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, replace
from typing import Any, Iterable

import numpy as np
from scipy.optimize import minimize

from .constrained_eq45_candidate import Eq45VelocityCandidate
from .constrained_eq45_profile_boundary_diagnostic import load_and_validate_receipt
from .constrained_eq45_profile_force_optimization import (
    DEFAULT_DERIVATIVE_STEPS,
    DEFAULT_FIT_STEP,
    DEFAULT_FTOL,
    DEFAULT_GTOL,
    PROFILE_BOUNDS,
    _force_design,
    _project_force,
    candidate_with_profile_pair,
)
from .constrained_force import RestrictedForce
from .constrained_pressure_compatibility import curl_pressure_free_momentum
from .constrained_restricted_force_curl_capacity import (
    FORCE_BOUNDS,
    REGISTERED_NU,
    _curl_force,
    _vector_max,
    _vector_rms,
    deterministic_probe_cloud,
)

Array = np.ndarray
TASK_ID = "CR005-EQ45-PHI10-PROFILE-FORCE-OPT-017"
MODE = (1, 0)
DEFAULT_MAX_FUNCTION_EVALUATIONS = 48


@dataclass(frozen=True)
class Phi10HoldoutLevel:
    spatial_step: float
    parent_rms_after_force: float
    optimized_rms_before_force: float
    optimized_rms_after_force: float
    optimized_max_after_force: float
    improvement_vs_parent: float


@dataclass(frozen=True)
class Eq45Phi10ProfileForceReport:
    task_id: str
    parent_candidate_sha256: str
    fit_seed: int
    holdout_seed: int
    sample_count: int
    fit_step: float
    derivative_steps: tuple[float, ...]
    profile_bounds: tuple[float, float]
    force_bounds: tuple[float, float]
    varied_profile_parameter: str
    frozen_phi_02: float
    frozen_F_02: float
    initial_phi_10: float
    optimized_phi_10: float
    optimizer: str
    max_function_evaluations: int
    optimizer_success: bool
    optimizer_message: str
    optimizer_function_evaluations: int
    optimizer_iterations: int
    parent_training_rms_after_force: float
    optimized_training_rms_before_force: float
    optimized_training_rms_after_force: float
    optimized_training_max_after_force: float
    training_improvement_vs_parent: float
    force_a: float
    force_c: float
    holdout_levels: tuple[Phi10HoldoutLevel, ...]
    velocity_parent_rms: float
    velocity_optimized_rms: float
    velocity_delta_rms: float
    velocity_relative_delta_rms: float
    velocity_delta_max: float
    optimized_candidate_sha256: str
    optimized_candidate_roundtrip_equal: bool
    velocity_changed: bool = True
    new_basis_added: bool = False
    pressure_fitted: bool = False
    forcing_family_changed: bool = False
    holdout_force_refit: bool = False
    canonical_candidate_promoted: bool = False
    pde_validated: bool = False
    visualization_ready_promoted: bool = False
    visual_correspondence_verified: bool = False
    paper_exact: bool = False
    openai_field_identified: bool = False

    def to_dict(self) -> dict[str, Any]:
        payload = asdict(self)
        payload["holdout_levels"] = [asdict(level) for level in self.holdout_levels]
        payload["interpretation"] = (
            "one existing bounded Phi(1,0) channel optimized after the checked "
            "Phi(0,2)/F(0,2) pair; the preregistered restricted force is projected "
            "only on training probes and frozen on holdout probes. Residual reduction "
            "is optimization/representation evidence only, not PDE or visual validation"
        )
        return payload


def candidate_with_phi10(
    candidate: Eq45VelocityCandidate,
    phi_10: float,
) -> Eq45VelocityCandidate:
    """Return ``candidate`` with only the existing ``Phi(1,0)`` coefficient changed."""
    if not isinstance(candidate, Eq45VelocityCandidate):
        raise TypeError("candidate must be Eq45VelocityCandidate")
    value = float(phi_10)
    limit = float(candidate.profile_basis.coefficient_limit)
    if not np.isfinite(value) or not (-limit <= value <= limit):
        raise ValueError(f"Phi(1,0) must lie in [{-limit}, {limit}]")
    try:
        index = candidate.profile_basis.mode_indices.index(MODE)
    except ValueError as exc:
        raise ValueError("candidate basis does not contain Phi(1,0)") from exc

    phi = list(candidate.profile_basis.phi_coefficients)
    phi[index] = value
    basis = replace(candidate.profile_basis, phi_coefficients=tuple(phi))
    return replace(candidate, profile_basis=basis)


def _validate_steps(steps: Iterable[float]) -> tuple[float, ...]:
    result = tuple(float(step) for step in steps)
    if len(result) < 3 or any(not np.isfinite(step) or step <= 0.0 for step in result):
        raise ValueError("derivative_steps must contain at least three positive finite levels")
    if any(next_step >= step for step, next_step in zip(result, result[1:])):
        raise ValueError("derivative_steps must be strictly decreasing")
    return result


def optimize_eq45_phi10_profile_force(
    *,
    receipt_path=None,
    derivative_steps: Iterable[float] = DEFAULT_DERIVATIVE_STEPS,
    max_function_evaluations: int = DEFAULT_MAX_FUNCTION_EVALUATIONS,
) -> tuple[Eq45VelocityCandidate, Eq45Phi10ProfileForceReport]:
    """Optimize only Phi(1,0), with the registered force projected on fit probes."""
    receipt = load_and_validate_receipt(receipt_path)
    contract = receipt["contract"]
    fit = receipt["fit"]
    expected_parent_sha = receipt["optimized_candidate"]["sha256"]
    steps = _validate_steps(derivative_steps)
    if not isinstance(max_function_evaluations, int) or max_function_evaluations <= 0:
        raise ValueError("max_function_evaluations must be a positive integer")

    seed = Eq45VelocityCandidate.seed()
    parent = candidate_with_profile_pair(seed, float(fit["Phi_02"]), float(fit["F_02"]))
    if parent.sha256 != expected_parent_sha:
        raise ValueError("parent candidate identity no longer matches the checked CR005 receipt")

    fit_seed = int(contract["fit_seed"])
    holdout_seed = int(contract["holdout_seed"])
    sample_count = int(contract["sample_count"])
    fit_step = float(contract["fit_step"])
    if not np.isclose(fit_step, DEFAULT_FIT_STEP, rtol=0.0, atol=0.0):
        raise ValueError("fit-step contract drift")
    fit_points, fit_times = deterministic_probe_cloud(fit_seed, sample_count)
    holdout_points, holdout_times = deterministic_probe_cloud(holdout_seed, sample_count)
    design = _force_design(fit_points, fit_times, fit_step)

    mode_index = parent.profile_basis.mode_indices.index(MODE)
    initial_phi10 = float(parent.profile_basis.phi_coefficients[mode_index])
    cache: dict[float, tuple[Eq45VelocityCandidate, Array, Any, Array]] = {}

    def evaluate(value: float):
        key = float(value)
        if key not in cache:
            trial = candidate_with_phi10(parent, key)
            target = curl_pressure_free_momentum(
                trial,
                fit_points,
                fit_times,
                nu=REGISTERED_NU,
                spatial_step=fit_step,
                time_step=fit_step,
            )[1]
            solution, after = _project_force(target, design)
            cache[key] = (trial, target, solution, after)
        return cache[key]

    parent_trial, _, parent_force_solution, parent_after = evaluate(initial_phi10)
    parent_rms = _vector_rms(parent_after)
    recorded_parent_rms = float(fit["curl_rms_after_force"])
    if not np.isclose(parent_rms, recorded_parent_rms, rtol=5e-10, atol=5e-10):
        raise ValueError("replayed parent objective no longer matches the CR005 receipt")
    recorded_parent_force = np.array([float(fit["force_a"]), float(fit["force_c"])])
    replay_parent_force = np.asarray(parent_force_solution.x, dtype=float)
    if not np.allclose(replay_parent_force, recorded_parent_force, rtol=5e-9, atol=5e-9):
        raise ValueError("replayed parent restricted force no longer matches the CR005 receipt")

    def objective(parameters: Array) -> float:
        value = float(np.asarray(parameters, dtype=float).reshape(-1)[0])
        return _vector_rms(evaluate(value)[3])

    nonlinear = minimize(
        objective,
        np.array([initial_phi10], dtype=float),
        method="L-BFGS-B",
        bounds=(PROFILE_BOUNDS,),
        options={
            "maxfun": int(max_function_evaluations),
            "maxiter": int(max_function_evaluations),
            "ftol": DEFAULT_FTOL,
            "gtol": DEFAULT_GTOL,
            "maxls": 20,
        },
    )
    if not np.all(np.isfinite(nonlinear.x)) or not np.isfinite(nonlinear.fun):
        raise RuntimeError("Phi(1,0) optimizer returned nonfinite state")

    optimized_phi10 = float(nonlinear.x[0])
    optimized, optimized_target, force_solution, optimized_after = evaluate(optimized_phi10)
    fitted_force = RestrictedForce(a=float(force_solution.x[0]), c=float(force_solution.x[1]))
    parent_force = RestrictedForce(
        a=float(parent_force_solution.x[0]), c=float(parent_force_solution.x[1])
    )

    levels: list[Phi10HoldoutLevel] = []
    for step in steps:
        parent_target = curl_pressure_free_momentum(
            parent_trial,
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
        parent_after_holdout = parent_target - _curl_force(
            parent_force, holdout_points, holdout_times, step
        )
        optimized_after_holdout = optimized_holdout - _curl_force(
            fitted_force, holdout_points, holdout_times, step
        )
        parent_holdout_rms = _vector_rms(parent_after_holdout)
        optimized_holdout_rms = _vector_rms(optimized_after_holdout)
        levels.append(
            Phi10HoldoutLevel(
                spatial_step=float(step),
                parent_rms_after_force=parent_holdout_rms,
                optimized_rms_before_force=_vector_rms(optimized_holdout),
                optimized_rms_after_force=optimized_holdout_rms,
                optimized_max_after_force=_vector_max(optimized_after_holdout),
                improvement_vs_parent=float(
                    1.0 - optimized_holdout_rms / max(parent_holdout_rms, 1e-15)
                ),
            )
        )

    parent_velocity = np.asarray(parent.at_points(holdout_points, holdout_times), dtype=float)
    optimized_velocity = np.asarray(optimized.at_points(holdout_points, holdout_times), dtype=float)
    velocity_delta = optimized_velocity - parent_velocity
    parent_velocity_rms = _vector_rms(parent_velocity)
    roundtrip = Eq45VelocityCandidate.from_dict(optimized.to_dict())

    report = Eq45Phi10ProfileForceReport(
        task_id=TASK_ID,
        parent_candidate_sha256=parent.sha256,
        fit_seed=fit_seed,
        holdout_seed=holdout_seed,
        sample_count=sample_count,
        fit_step=fit_step,
        derivative_steps=steps,
        profile_bounds=PROFILE_BOUNDS,
        force_bounds=FORCE_BOUNDS,
        varied_profile_parameter="Phi(1,0)",
        frozen_phi_02=float(fit["Phi_02"]),
        frozen_F_02=float(fit["F_02"]),
        initial_phi_10=initial_phi10,
        optimized_phi_10=optimized_phi10,
        optimizer="scipy.optimize.minimize(method='L-BFGS-B') with nested lsq_linear variable projection",
        max_function_evaluations=int(max_function_evaluations),
        optimizer_success=bool(nonlinear.success),
        optimizer_message=str(nonlinear.message),
        optimizer_function_evaluations=int(nonlinear.nfev),
        optimizer_iterations=int(nonlinear.nit),
        parent_training_rms_after_force=parent_rms,
        optimized_training_rms_before_force=_vector_rms(optimized_target),
        optimized_training_rms_after_force=_vector_rms(optimized_after),
        optimized_training_max_after_force=_vector_max(optimized_after),
        training_improvement_vs_parent=float(
            1.0 - _vector_rms(optimized_after) / max(parent_rms, 1e-15)
        ),
        force_a=float(force_solution.x[0]),
        force_c=float(force_solution.x[1]),
        holdout_levels=tuple(levels),
        velocity_parent_rms=parent_velocity_rms,
        velocity_optimized_rms=_vector_rms(optimized_velocity),
        velocity_delta_rms=_vector_rms(velocity_delta),
        velocity_relative_delta_rms=float(
            _vector_rms(velocity_delta) / max(parent_velocity_rms, 1e-15)
        ),
        velocity_delta_max=_vector_max(velocity_delta),
        optimized_candidate_sha256=optimized.sha256,
        optimized_candidate_roundtrip_equal=(roundtrip.sha256 == optimized.sha256),
    )
    return optimized, report
