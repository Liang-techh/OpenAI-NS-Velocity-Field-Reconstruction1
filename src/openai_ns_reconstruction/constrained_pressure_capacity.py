"""Pressure-capacity diagnostic separating frozen generalization from basis ceiling.

This module deliberately refits the *held-out* residual only as a diagnostic
capacity ceiling.  That refit is data leakage by construction and therefore
must never be reported as independent validation or used to upgrade
``pde_validated``.
"""
from dataclasses import asdict, dataclass, replace

import numpy as np
from scipy.optimize import lsq_linear

from .constrained_pressure_holdout import _independent_momentum
from .constrained_pressure_projection import (
    apply_projected_pressure,
    project_pressure_coefficients,
)


@dataclass(frozen=True)
class PressureCapacityDiagnostic:
    training_coefficients: tuple[float, ...]
    training_rank: int
    training_condition_number: float
    training_active_lower: int
    training_active_upper: int
    training_rms_after: float
    frozen_holdout_rms: float
    capacity_coefficients: tuple[float, ...]
    holdout_rank: int
    holdout_condition_number: float
    holdout_active_lower: int
    holdout_active_upper: int
    capacity_holdout_rms: float
    capacity_fraction_of_frozen: float | None
    recoverable_fraction_of_frozen: float | None
    velocity_max_change: float
    training_seed: int
    holdout_seed: int
    validation_step: float
    holdout_points: int
    truth_boundary: str

    def to_dict(self):
        return asdict(self)


def _pressure_gradient_matrix_fourth_order(candidate, points, times, *, step):
    """Differentiate only the preregistered pressure basis with validation stencil."""
    x = np.asarray(points, dtype=float)
    t = np.asarray(times, dtype=float)
    if (x.ndim != 2 or x.shape[1] != 3 or len(x) == 0
            or not np.isfinite(x).all()):
        raise ValueError("points must be a nonempty finite (N,3) array")
    if t.shape not in ((), (len(x),)) or not np.isfinite(t).all():
        raise ValueError("times must be a finite scalar or length-N array")
    if not np.isfinite(step) or step <= 0:
        raise ValueError("step must be finite and positive")

    probe = np.asarray(candidate.pressure_basis(x, t), dtype=float)
    if (probe.ndim != 2 or probe.shape[0] != len(x)
            or probe.shape[1] == 0 or not np.isfinite(probe).all()):
        raise ValueError("pressure_basis must return a finite nonempty (N,K) array")

    gradient = np.empty((len(x), 3, probe.shape[1]), dtype=float)
    for j in range(3):
        d = np.eye(3)[j] * step
        pm2, pm1, pp1, pp2 = [
            np.asarray(candidate.pressure_basis(x + k * d, t), dtype=float)
            for k in (-2.0, -1.0, 1.0, 2.0)
        ]
        if any(value.shape != probe.shape for value in (pm2, pm1, pp1, pp2)):
            raise ValueError("pressure_basis shape changed under validation stencil")
        gradient[:, j, :] = (pm2 - 8 * pm1 + 8 * pp1 - pp2) / (12 * step)
    if not np.isfinite(gradient).all():
        raise ValueError("nonfinite pressure-gradient design matrix")
    return gradient


def _matrix_diagnostics(matrix):
    singular = np.linalg.svd(matrix, compute_uv=False)
    rank = int(np.linalg.matrix_rank(matrix))
    condition = (
        float("inf")
        if singular.size == 0 or singular[-1] == 0
        else float(singular[0] / singular[-1])
    )
    return rank, condition


def diagnose_pressure_capacity(
    candidate,
    force,
    training_points,
    training_times,
    holdout_points,
    holdout_times,
    nu,
    *,
    training_seed,
    holdout_seed,
    training_step=0.001,
    validation_step=0.005,
    time_bounds=(0.25, 0.75),
    lower=-1.0,
    upper=1.0,
    tol=1e-10,
    max_iter=200,
):
    """Compare frozen train-fit pressure with a diagnostic held-out refit ceiling.

    The first solve is a normal bounded training projection.  The second solve
    intentionally sees held-out residuals and exists only to estimate how much
    of the frozen held-out floor the *same pressure basis* could remove if it
    were allowed to overfit those samples.  It is not validation.
    """
    if int(training_seed) == int(holdout_seed):
        raise ValueError("training and holdout seeds must differ")
    if (not np.isfinite([nu, training_step, validation_step, lower, upper, tol,
                         *time_bounds]).all()
            or nu <= 0 or training_step <= 0 or validation_step <= 0
            or lower >= upper or tol <= 0 or max_iter <= 0):
        raise ValueError("invalid viscosity, steps, bounds, tolerance, or budget")
    if time_bounds[0] >= time_bounds[1]:
        raise ValueError("invalid time bounds")

    current = np.asarray(candidate.pressure_coefficients, dtype=float)
    if current.ndim != 1 or current.size == 0 or not np.isfinite(current).all():
        raise ValueError("candidate must expose finite 1-D pressure_coefficients")

    training = project_pressure_coefficients(
        candidate,
        force,
        training_points,
        training_times,
        nu,
        step=training_step,
        lower=lower,
        upper=upper,
        tol=tol,
        max_iter=max_iter,
    )
    frozen = apply_projected_pressure(candidate, training)

    frozen_momentum = _independent_momentum(
        frozen,
        force,
        holdout_points,
        holdout_times,
        nu=nu,
        step=validation_step,
        time_bounds=time_bounds,
    )
    zero = replace(candidate, pressure_coefficients=tuple(np.zeros_like(current)))
    base = _independent_momentum(
        zero,
        force,
        holdout_points,
        holdout_times,
        nu=nu,
        step=validation_step,
        time_bounds=time_bounds,
    )
    design = _pressure_gradient_matrix_fourth_order(
        zero, holdout_points, holdout_times, step=validation_step
    )
    if base.shape != design.shape[:2]:
        raise ValueError("independent residual and pressure basis shapes are incompatible")

    matrix = design.reshape(-1, current.size)
    rhs = -base.reshape(-1)
    fit = lsq_linear(
        matrix, rhs, bounds=(lower, upper), tol=tol, max_iter=max_iter
    )
    if not fit.success:
        raise RuntimeError("bounded held-out capacity refit failed: " + fit.message)

    capacity = base + np.einsum("nck,k->nc", design, fit.x)
    frozen_norm = np.linalg.norm(frozen_momentum, axis=1)
    capacity_norm = np.linalg.norm(capacity, axis=1)
    frozen_rms = float(np.sqrt(np.mean(frozen_norm ** 2)))
    capacity_rms = float(np.sqrt(np.mean(capacity_norm ** 2)))
    ratio_floor = 100.0 * np.finfo(float).eps
    if frozen_rms > ratio_floor:
        fraction = float(capacity_rms / frozen_rms)
        recoverable = float(max(0.0, 1.0 - fraction))
    else:
        fraction = None
        recoverable = None

    holdout_rank, holdout_condition = _matrix_diagnostics(matrix)
    scale = max(1.0, abs(lower), abs(upper))
    active_tol = 1e-8 * scale
    active_lower = int(np.count_nonzero(np.abs(fit.x - lower) <= active_tol))
    active_upper = int(np.count_nonzero(np.abs(fit.x - upper) <= active_tol))

    x = np.asarray(holdout_points, dtype=float)
    t = np.asarray(holdout_times, dtype=float)
    velocity_change = float(np.max(np.abs(
        candidate.velocity(x, t) - frozen.velocity(x, t)
    )))

    return PressureCapacityDiagnostic(
        training_coefficients=training.coefficients,
        training_rank=training.rank,
        training_condition_number=training.condition_number,
        training_active_lower=training.active_lower,
        training_active_upper=training.active_upper,
        training_rms_after=training.residual_rms_after,
        frozen_holdout_rms=frozen_rms,
        capacity_coefficients=tuple(float(v) for v in fit.x),
        holdout_rank=holdout_rank,
        holdout_condition_number=holdout_condition,
        holdout_active_lower=active_lower,
        holdout_active_upper=active_upper,
        capacity_holdout_rms=capacity_rms,
        capacity_fraction_of_frozen=fraction,
        recoverable_fraction_of_frozen=recoverable,
        velocity_max_change=velocity_change,
        training_seed=int(training_seed),
        holdout_seed=int(holdout_seed),
        validation_step=float(validation_step),
        holdout_points=len(x),
        truth_boundary=(
            "held-out refit is a diagnostic capacity ceiling with deliberate "
            "data leakage; it is not independent validation, PDE acceptance, "
            "visual correspondence, or evidence for an exact OpenAI field"
        ),
    )
