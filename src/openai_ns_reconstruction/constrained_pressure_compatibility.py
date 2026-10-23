"""Independent pressure-Poisson compatibility diagnostics for CR004.

For an incompressible velocity field, taking divergence of

    u_t + (u·grad)u + grad(p) - nu*laplacian(u) = f

gives the pressure compatibility condition

    laplacian(p) + sum_ij (partial_j u_i)(partial_i u_j) - div(f) = 0,

provided div(u)=0.  This module evaluates that scalar defect with a spatial
finite-difference operator that is separate from the optimizer loss.  It does
not alter pressure or forcing and therefore cannot erase the momentum residual.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Callable

import numpy as np


ArrayField = Callable[[np.ndarray, np.ndarray], np.ndarray]
ScalarField = Callable[[np.ndarray, np.ndarray], np.ndarray]


@dataclass(frozen=True)
class PressureCompatibilityReport:
    sample_count: int
    step: float
    poisson_defect_max_abs: float
    poisson_defect_rms: float
    pressure_laplacian_rms: float
    velocity_quadratic_rms: float
    force_divergence_rms: float
    velocity_divergence_max_abs: float
    exterior_pressure_max_abs: float | None

    def to_dict(self) -> dict[str, float | int | None]:
        return asdict(self)


def _samples(points, times) -> tuple[np.ndarray, np.ndarray]:
    p = np.asarray(points, dtype=float)
    if p.ndim != 2 or p.shape[1] != 3 or len(p) == 0:
        raise ValueError("points must have shape (n,3) with n>0")
    if not np.all(np.isfinite(p)):
        raise ValueError("points must be finite")
    t = np.asarray(times, dtype=float)
    if t.ndim == 0:
        t = np.full(len(p), float(t))
    elif t.shape != (len(p),):
        raise ValueError("times must be scalar or shape (n,)")
    if not np.all(np.isfinite(t)):
        raise ValueError("times must be finite")
    return p, t


def _vector_value(field: ArrayField, points: np.ndarray, times: np.ndarray, name: str) -> np.ndarray:
    value = np.asarray(field(points, times), dtype=float)
    if value.shape != points.shape or not np.all(np.isfinite(value)):
        raise ValueError(f"{name} must return finite shape (n,3)")
    return value


def _scalar_value(field: ScalarField, points: np.ndarray, times: np.ndarray, name: str) -> np.ndarray:
    value = np.asarray(field(points, times), dtype=float)
    if value.shape != (len(points),) or not np.all(np.isfinite(value)):
        raise ValueError(f"{name} must return finite shape (n,)")
    return value


def pressure_poisson_terms(
    velocity: ArrayField,
    pressure: ScalarField,
    force: ArrayField,
    points,
    times,
    *,
    step: float = 1e-3,
) -> dict[str, np.ndarray]:
    """Evaluate the pressure-Poisson terms on fixed samples.

    The returned ``defect`` is zero for an exactly incompressible solution with
    a compatible pressure and prescribed force.  A small defect is not by
    itself a Navier--Stokes acceptance result because momentum balance and
    divergence must be checked independently.
    """
    p, t = _samples(points, times)
    if not np.isfinite(step) or step <= 0:
        raise ValueError("step must be positive and finite")
    h = float(step)
    p0 = _scalar_value(pressure, p, t, "pressure")
    grad_u = np.empty((len(p), 3, 3), dtype=float)
    lap_p = np.zeros(len(p), dtype=float)
    div_f = np.zeros(len(p), dtype=float)
    for axis in range(3):
        delta = np.zeros(3, dtype=float)
        delta[axis] = h
        plus = p + delta
        minus = p - delta
        u_plus = _vector_value(velocity, plus, t, "velocity")
        u_minus = _vector_value(velocity, minus, t, "velocity")
        grad_u[:, :, axis] = (u_plus - u_minus) / (2.0 * h)
        p_plus = _scalar_value(pressure, plus, t, "pressure")
        p_minus = _scalar_value(pressure, minus, t, "pressure")
        lap_p += (p_plus - 2.0 * p0 + p_minus) / (h * h)
        f_plus = _vector_value(force, plus, t, "force")
        f_minus = _vector_value(force, minus, t, "force")
        div_f += (f_plus[:, axis] - f_minus[:, axis]) / (2.0 * h)
    quadratic = np.einsum("nij,nji->n", grad_u, grad_u)
    divergence = np.trace(grad_u, axis1=1, axis2=2)
    defect = lap_p + quadratic - div_f
    return {
        "pressure_laplacian": lap_p,
        "velocity_quadratic": quadratic,
        "force_divergence": div_f,
        "velocity_divergence": divergence,
        "defect": defect,
    }


def pressure_compatibility_report(
    velocity: ArrayField,
    pressure: ScalarField,
    force: ArrayField,
    points,
    times,
    *,
    step: float = 1e-3,
    exterior_points=None,
    exterior_times=None,
) -> PressureCompatibilityReport:
    """Summarize pressure compatibility without fitting any coefficient."""
    p, t = _samples(points, times)
    terms = pressure_poisson_terms(velocity, pressure, force, p, t, step=step)

    exterior_max: float | None = None
    if exterior_points is not None:
        ep, et = _samples(exterior_points, t[0] if exterior_times is None else exterior_times)
        exterior_max = float(np.max(np.abs(_scalar_value(pressure, ep, et, "pressure"))))

    def rms(value: np.ndarray) -> float:
        return float(np.sqrt(np.mean(np.square(value))))

    return PressureCompatibilityReport(
        sample_count=len(p),
        step=float(step),
        poisson_defect_max_abs=float(np.max(np.abs(terms["defect"]))),
        poisson_defect_rms=rms(terms["defect"]),
        pressure_laplacian_rms=rms(terms["pressure_laplacian"]),
        velocity_quadratic_rms=rms(terms["velocity_quadratic"]),
        force_divergence_rms=rms(terms["force_divergence"]),
        velocity_divergence_max_abs=float(np.max(np.abs(terms["velocity_divergence"]))),
        exterior_pressure_max_abs=exterior_max,
    )
