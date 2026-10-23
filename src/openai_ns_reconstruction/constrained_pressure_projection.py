"""Bounded pressure projection for a fixed constrained velocity/force pair.

This module eliminates only the predeclared linear pressure block.  It never
changes the velocity or force, so it cannot manufacture a PDE pass via a
residual-defined forcing term.  The projection uses the training derivative
operator; independent validation must still be run separately.
"""
from dataclasses import dataclass, replace
import numpy as np
from scipy.optimize import lsq_linear


@dataclass(frozen=True)
class PressureProjectionResult:
    coefficients: tuple[float, ...]
    rank: int
    condition_number: float
    active_lower: int
    active_upper: int
    residual_rms_before: float
    residual_rms_after: float
    residual_max_before: float
    residual_max_after: float
    solver_status: int
    solver_message: str

    @property
    def rms_ratio(self):
        return self.residual_rms_after / self.residual_rms_before if self.residual_rms_before else 0.0


def pressure_gradient_matrix(candidate, points, times, *, step=0.001):
    """Return d/dx of the candidate's linear pressure basis as (N,3,K)."""
    x = np.asarray(points, dtype=float)
    t = np.asarray(times, dtype=float)
    if x.ndim != 2 or x.shape[1] != 3 or not np.isfinite(x).all():
        raise ValueError('points must be a finite (N,3) array')
    if t.shape not in ((), (len(x),)) or not np.isfinite(t).all():
        raise ValueError('times must be a finite scalar or length-N array')
    if not np.isfinite(step) or step <= 0:
        raise ValueError('step must be finite and positive')
    probe = np.asarray(candidate.pressure_basis(x, t), dtype=float)
    if probe.ndim != 2 or probe.shape[0] != len(x) or not np.isfinite(probe).all():
        raise ValueError('pressure_basis must return a finite (N,K) array')
    k = probe.shape[1]
    gradient = np.empty((len(x), 3, k), dtype=float)
    for j in range(3):
        d = np.eye(3)[j] * step
        plus = np.asarray(candidate.pressure_basis(x + d, t), dtype=float)
        minus = np.asarray(candidate.pressure_basis(x - d, t), dtype=float)
        if plus.shape != probe.shape or minus.shape != probe.shape:
            raise ValueError('pressure_basis shape changed under spatial stencil')
        gradient[:, j, :] = (plus - minus) / (2 * step)
    if not np.isfinite(gradient).all():
        raise ValueError('nonfinite pressure gradient basis')
    return gradient


def project_pressure_coefficients(candidate, force, points, times, nu, *, step=0.001,
                                  lower=-1.0, upper=1.0, tol=1e-10,
                                  max_iter=200):
    """Solve the bounded linear pressure subproblem for fixed velocity/force.

    The objective is the second-order training-operator L2 momentum residual.
    This is an exact bounded least-squares elimination for that objective, not
    an independent PDE validation and not the quartic outer training loss.
    """
    from .constrained_optimize import training_residual

    current = np.asarray(candidate.pressure_coefficients, dtype=float)
    if current.ndim != 1 or current.size == 0 or not np.isfinite(current).all():
        raise ValueError('candidate must expose finite 1-D pressure_coefficients')
    if not np.isfinite([nu, lower, upper, tol]).all() or nu <= 0 or lower >= upper or tol <= 0:
        raise ValueError('invalid viscosity, bounds, or tolerance')
    if max_iter <= 0:
        raise ValueError('max_iter must be positive')
    if np.any(current < lower) or np.any(current > upper):
        raise ValueError('current pressure coefficients are outside projection bounds')

    zero = replace(candidate, pressure_coefficients=tuple(np.zeros_like(current)))
    base = np.asarray(training_residual(zero, force, points, times, nu, step), dtype=float)
    design = pressure_gradient_matrix(zero, points, times, step=step)
    if base.shape != design.shape[:2]:
        raise ValueError('training residual and pressure basis shapes are incompatible')

    matrix = design.reshape(-1, current.size)
    rhs = -base.reshape(-1)
    fit = lsq_linear(matrix, rhs, bounds=(lower, upper), tol=tol, max_iter=max_iter)
    if not fit.success:
        raise RuntimeError('bounded pressure projection failed: ' + fit.message)

    before = base + np.einsum('nck,k->nc', design, current)
    after = base + np.einsum('nck,k->nc', design, fit.x)
    before_norm = np.linalg.norm(before, axis=1)
    after_norm = np.linalg.norm(after, axis=1)
    singular = np.linalg.svd(matrix, compute_uv=False)
    rank = int(np.linalg.matrix_rank(matrix))
    condition = float('inf') if singular.size == 0 or singular[-1] == 0 else float(singular[0] / singular[-1])
    scale = max(1.0, abs(lower), abs(upper))
    active_tol = 1e-8 * scale
    return PressureProjectionResult(
        coefficients=tuple(float(v) for v in fit.x),
        rank=rank,
        condition_number=condition,
        active_lower=int(np.count_nonzero(np.abs(fit.x - lower) <= active_tol)),
        active_upper=int(np.count_nonzero(np.abs(fit.x - upper) <= active_tol)),
        residual_rms_before=float(np.sqrt(np.mean(before_norm ** 2))),
        residual_rms_after=float(np.sqrt(np.mean(after_norm ** 2))),
        residual_max_before=float(before_norm.max()),
        residual_max_after=float(after_norm.max()),
        solver_status=int(fit.status),
        solver_message=str(fit.message),
    )


def apply_projected_pressure(candidate, result):
    """Return a copy with projected pressure coefficients; velocity/force stay fixed."""
    if len(result.coefficients) != len(candidate.pressure_coefficients):
        raise ValueError('projection/candidate pressure dimension mismatch')
    return replace(candidate, pressure_coefficients=result.coefficients)
