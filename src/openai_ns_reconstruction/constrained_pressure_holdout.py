"""Freeze a bounded pressure projection, then audit it on disjoint data.

The pressure fit uses the training operator from constrained_pressure_projection.
The held-out check uses the independent fourth-order validation operator and never
refits pressure on held-out points. Velocity and force are immutable throughout.
"""
from dataclasses import dataclass

import numpy as np

from .constrained_pressure_projection import (
    apply_projected_pressure,
    project_pressure_coefficients,
)
from .constrained_validation import residual


@dataclass(frozen=True)
class PressureHoldoutAudit:
    training_seed: int
    holdout_seed: int
    projected_coefficients: tuple[float, ...]
    training_rms_before: float
    training_rms_after: float
    heldout_rms_before: float
    heldout_rms_after: float
    heldout_max_before: float
    heldout_max_after: float
    velocity_max_change: float
    validation_step: float
    validation_points: int

    @property
    def heldout_rms_ratio(self):
        return (self.heldout_rms_after / self.heldout_rms_before
                if self.heldout_rms_before else 0.0)


def _independent_momentum(candidate, force, points, times, *, nu, step,
                          time_bounds):
    x = np.asarray(points, dtype=float)
    t = np.asarray(times, dtype=float)
    if (x.ndim != 2 or x.shape[1] != 3 or len(x) == 0
            or not np.isfinite(x).all()):
        raise ValueError('points must be a nonempty finite (N,3) array')
    if t.shape != (len(x),) or not np.isfinite(t).all():
        raise ValueError('times must be a finite length-N array')

    rows = []
    for value in np.unique(t):
        mask = t == value
        result = residual(
            candidate.velocity,
            candidate.pressure,
            force,
            x[mask],
            float(value),
            nu=nu,
            step=step,
            time_bounds=time_bounds,
        )
        momentum = np.asarray(result['momentum'], dtype=float)
        if (momentum.shape != (int(np.count_nonzero(mask)), 3)
                or not np.isfinite(momentum).all()):
            raise ValueError('independent residual returned invalid momentum')
        rows.append(momentum)
    return np.concatenate(rows, axis=0)


def audit_projected_pressure(candidate, force, training_points, training_times,
                             holdout_points, holdout_times, nu, *,
                             training_seed, holdout_seed,
                             training_step=0.001, validation_step=0.005,
                             time_bounds=(0.25, 0.75),
                             lower=-1.0, upper=1.0):
    """Fit pressure once on training data and evaluate frozen pressure held out.

    This diagnoses pressure-block capacity only. It cannot establish PDE validity:
    the velocity and restricted force are held fixed, and acceptance thresholds are
    not changed or consulted here.
    """
    if int(training_seed) == int(holdout_seed):
        raise ValueError('training and holdout seeds must differ')
    if (not np.isfinite([nu, training_step, validation_step, *time_bounds]).all()
            or nu <= 0 or training_step <= 0 or validation_step <= 0):
        raise ValueError('invalid viscosity, steps, or time bounds')
    if time_bounds[0] >= time_bounds[1]:
        raise ValueError('invalid time bounds')

    projection = project_pressure_coefficients(
        candidate,
        force,
        training_points,
        training_times,
        nu,
        step=training_step,
        lower=lower,
        upper=upper,
    )
    projected = apply_projected_pressure(candidate, projection)

    before = _independent_momentum(
        candidate,
        force,
        holdout_points,
        holdout_times,
        nu=nu,
        step=validation_step,
        time_bounds=time_bounds,
    )
    after = _independent_momentum(
        projected,
        force,
        holdout_points,
        holdout_times,
        nu=nu,
        step=validation_step,
        time_bounds=time_bounds,
    )
    before_norm = np.linalg.norm(before, axis=1)
    after_norm = np.linalg.norm(after, axis=1)

    x = np.asarray(holdout_points, dtype=float)
    t = np.asarray(holdout_times, dtype=float)
    velocity_change = float(np.max(np.abs(
        candidate.velocity(x, t) - projected.velocity(x, t)
    )))

    return PressureHoldoutAudit(
        training_seed=int(training_seed),
        holdout_seed=int(holdout_seed),
        projected_coefficients=projection.coefficients,
        training_rms_before=projection.residual_rms_before,
        training_rms_after=projection.residual_rms_after,
        heldout_rms_before=float(np.sqrt(np.mean(before_norm ** 2))),
        heldout_rms_after=float(np.sqrt(np.mean(after_norm ** 2))),
        heldout_max_before=float(before_norm.max()),
        heldout_max_after=float(after_norm.max()),
        velocity_max_change=velocity_change,
        validation_step=float(validation_step),
        validation_points=len(x),
    )
