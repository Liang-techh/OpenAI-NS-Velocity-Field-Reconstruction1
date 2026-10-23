"""Pressure-only capacity diagnostic for the registered local pressure family.

This diagnostic freezes velocity, forcing, the base pressure, sampling points,
and derivative stencil. It varies only ``candidate.pressure_coefficients``,
whose preregistered bounds are [-1, 1], and uses the candidate's declared
``pressure_basis``. The least-squares fit is a same-sample family-capacity
ceiling, not independent validation and not evidence of PDE acceptance.
"""
from __future__ import annotations

from dataclasses import dataclass, asdict
from typing import Any

import numpy as np
from scipy.optimize import lsq_linear

from .constrained_validation import residual


PRESSURE_COEFFICIENT_LOWER = -1.0
PRESSURE_COEFFICIENT_UPPER = 1.0


@dataclass(frozen=True)
class PressureCapacityReport:
    coefficient_count: int
    design_rank: int
    design_nullity: int
    design_condition: float | None
    singular_values: tuple[float, ...]
    current_coefficients: tuple[float, ...]
    capacity_coefficients: tuple[float, ...]
    active_lower_bounds: tuple[int, ...]
    active_upper_bounds: tuple[int, ...]
    current_momentum_rms: float
    current_momentum_max: float
    capacity_momentum_rms: float
    capacity_momentum_max: float
    recoverable_fraction_rms: float
    same_sample_capacity: bool = True
    pde_validated: bool = False
    paper_exact: bool = False
    truth_boundary: str = (
        "pressure-only same-sample capacity within the preregistered [-1,1] "
        "local pressure coefficients; velocity, forcing, base pressure, points "
        "and derivative stencil are frozen; not independent PDE validation"
    )

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def _validate_inputs(candidate, points, time, step):
    x = np.asarray(points, dtype=float)
    if x.ndim != 2 or x.shape[1] != 3 or len(x) == 0 or not np.all(np.isfinite(x)):
        raise ValueError("points must be a nonempty finite (N,3) array")
    if not np.isfinite(time):
        raise ValueError("time must be finite")
    if not np.isfinite(step) or step <= 0:
        raise ValueError("step must be positive and finite")
    if not hasattr(candidate, "pressure_basis") or not hasattr(candidate, "pressure_coefficients"):
        raise TypeError("candidate must expose pressure_basis and pressure_coefficients")
    coeff = np.asarray(candidate.pressure_coefficients, dtype=float)
    if coeff.ndim != 1 or coeff.size == 0 or not np.all(np.isfinite(coeff)):
        raise ValueError("pressure_coefficients must be a nonempty finite vector")
    if np.any(coeff < PRESSURE_COEFFICIENT_LOWER) or np.any(coeff > PRESSURE_COEFFICIENT_UPPER):
        raise ValueError("pressure_coefficients must respect preregistered [-1,1] bounds")
    basis = np.asarray(candidate.pressure_basis(x, time), dtype=float)
    if basis.shape != (len(x), len(coeff)) or not np.all(np.isfinite(basis)):
        raise ValueError("pressure_basis must return finite shape (N,K)")
    return x, coeff


def _pressure_gradient_design(candidate, points, time, step, coefficient_count):
    design = np.empty((len(points), 3, coefficient_count), dtype=float)
    for axis in range(3):
        d = np.eye(3)[axis] * step
        bm2, bm1, bp1, bp2 = [
            np.asarray(candidate.pressure_basis(points + k*d, time), dtype=float)
            for k in (-2, -1, 1, 2)
        ]
        for value in (bm2, bm1, bp1, bp2):
            if value.shape != (len(points), coefficient_count) or not np.all(np.isfinite(value)):
                raise ValueError("pressure_basis stencil evaluation is invalid")
        design[:, axis, :] = (bm2 - 8*bm1 + 8*bp1 - bp2) / (12*step)
    return design


def diagnose_pressure_family_capacity(
    candidate,
    force,
    points,
    time,
    *,
    nu=0.01,
    step=0.005,
    time_bounds=(0.25, 0.75),
    rank_rtol=1e-10,
) -> PressureCapacityReport:
    """Compute the best same-sample momentum fit from pressure coefficients only."""
    x, current = _validate_inputs(candidate, points, time, step)
    if not np.isfinite(rank_rtol) or not 0 < rank_rtol < 1:
        raise ValueError("rank_rtol must lie in (0,1)")

    base = residual(
        candidate.velocity,
        candidate.pressure,
        force,
        x,
        time,
        nu=nu,
        step=step,
        time_bounds=time_bounds,
    )
    momentum = np.asarray(base["momentum"], dtype=float)
    if momentum.shape != (len(x), 3) or not np.all(np.isfinite(momentum)):
        raise ValueError("independent residual returned invalid momentum")

    design3 = _pressure_gradient_design(candidate, x, time, step, len(current))
    design = design3.reshape(-1, len(current))
    target = -momentum.reshape(-1)

    lower_delta = PRESSURE_COEFFICIENT_LOWER - current
    upper_delta = PRESSURE_COEFFICIENT_UPPER - current
    fit = lsq_linear(
        design,
        target,
        bounds=(lower_delta, upper_delta),
        method="trf",
        lsmr_tol="auto",
    )
    if not fit.success or not np.all(np.isfinite(fit.x)):
        raise RuntimeError("bounded pressure capacity solve failed")

    delta = np.asarray(fit.x, dtype=float)
    capacity_coeff = current + delta
    remaining = momentum + (design @ delta).reshape(len(x), 3)

    singular = np.linalg.svd(design, compute_uv=False)
    if singular.size and singular[0] > 0:
        threshold = rank_rtol * singular[0]
        rank = int(np.count_nonzero(singular > threshold))
        condition = float(singular[0] / singular[rank-1]) if rank else None
    else:
        rank = 0
        condition = None

    tol = 256*np.finfo(float).eps
    lower_active = tuple(int(i) for i,v in enumerate(capacity_coeff)
                         if abs(v-PRESSURE_COEFFICIENT_LOWER) <= tol)
    upper_active = tuple(int(i) for i,v in enumerate(capacity_coeff)
                         if abs(v-PRESSURE_COEFFICIENT_UPPER) <= tol)

    current_norms = np.linalg.norm(momentum, axis=1)
    capacity_norms = np.linalg.norm(remaining, axis=1)
    current_rms = float(np.sqrt(np.mean(current_norms**2)))
    capacity_rms = float(np.sqrt(np.mean(capacity_norms**2)))
    recoverable = 0.0 if current_rms == 0 else float(1.0-capacity_rms/current_rms)

    return PressureCapacityReport(
        coefficient_count=len(current),
        design_rank=rank,
        design_nullity=len(current)-rank,
        design_condition=condition,
        singular_values=tuple(float(v) for v in singular),
        current_coefficients=tuple(float(v) for v in current),
        capacity_coefficients=tuple(float(v) for v in capacity_coeff),
        active_lower_bounds=lower_active,
        active_upper_bounds=upper_active,
        current_momentum_rms=current_rms,
        current_momentum_max=float(current_norms.max()),
        capacity_momentum_rms=capacity_rms,
        capacity_momentum_max=float(capacity_norms.max()),
        recoverable_fraction_rms=recoverable,
    )
