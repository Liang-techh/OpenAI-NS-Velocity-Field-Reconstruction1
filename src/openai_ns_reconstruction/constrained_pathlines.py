"""Truth-bounded time-dependent pathlines for a frozen velocity callable.

This module integrates particle paths

    dX/dt = velocity(X, t)

using SciPy's public ``solve_ivp`` API. It is a visualization/kinematics
utility only: stable or visually plausible trajectories are not Navier--Stokes
validation and do not identify any hidden OpenAI field.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np
from scipy.integrate import solve_ivp

VelocityCallable = Callable[[np.ndarray, float], np.ndarray]


@dataclass(frozen=True)
class PathlineBundle:
    """Sampled trajectories from a time-dependent velocity field."""

    times: np.ndarray
    positions: np.ndarray
    speeds: np.ndarray
    seed_positions: np.ndarray
    rtol: float
    atol: float
    max_step: float
    solver_method: str
    claim_scope: str = "visualization_kinematics_only"
    pde_validated: bool = False
    paper_exact: bool = False
    openai_field_identified: bool = False
    blowup_proved: bool = False


def _readonly(a: np.ndarray) -> np.ndarray:
    out = np.asarray(a, dtype=float).copy()
    out.setflags(write=False)
    return out


def _validate_seeds(seed_positions: np.ndarray) -> np.ndarray:
    seeds = np.asarray(seed_positions, dtype=float)
    if seeds.ndim != 2 or seeds.shape[1] != 3 or seeds.shape[0] == 0:
        raise ValueError("seed_positions must have shape (n, 3) with n >= 1")
    if not np.all(np.isfinite(seeds)):
        raise ValueError("seed_positions must be finite")
    return seeds


def _validate_times(times: np.ndarray) -> np.ndarray:
    t = np.asarray(times, dtype=float)
    if t.ndim != 1 or t.size < 2:
        raise ValueError("times must be a 1D array with at least two entries")
    if not np.all(np.isfinite(t)):
        raise ValueError("times must be finite")
    if not np.all(np.diff(t) > 0):
        raise ValueError("times must be strictly increasing")
    return t


def integrate_pathlines(
    velocity: VelocityCallable,
    seed_positions: np.ndarray,
    times: np.ndarray,
    *,
    allowed_time_interval: tuple[float, float] | None = None,
    rtol: float = 1e-7,
    atol: float = 1e-9,
    max_step: float | None = None,
    method: str = "RK45",
) -> PathlineBundle:
    """Integrate time-dependent particle paths through ``velocity``.

    ``velocity(points, t)`` must accept points of shape ``(n, 3)`` and return
    an array with the same shape. One ODE system containing all seeds is
    integrated together, so returned trajectories share the same sample times.

    ``allowed_time_interval`` is a governance guard, not a solver domain.
    When supplied, every requested time must lie inside the declared interval.

    The result is visualization/kinematics evidence only. It does not evaluate
    PDE residuals, divergence, forcing, or acceptance thresholds.
    """

    if not callable(velocity):
        raise TypeError("velocity must be callable")

    seeds = _validate_seeds(seed_positions)
    t_eval = _validate_times(times)

    if allowed_time_interval is not None:
        if len(allowed_time_interval) != 2:
            raise ValueError("allowed_time_interval must be (start, end)")
        lo, hi = map(float, allowed_time_interval)
        if not (np.isfinite(lo) and np.isfinite(hi) and lo < hi):
            raise ValueError("allowed_time_interval must be finite and increasing")
        if t_eval[0] < lo or t_eval[-1] > hi:
            raise ValueError("requested times exceed allowed_time_interval")

    rtol = float(rtol)
    atol = float(atol)
    if not (np.isfinite(rtol) and rtol > 0):
        raise ValueError("rtol must be finite and positive")
    if not (np.isfinite(atol) and atol > 0):
        raise ValueError("atol must be finite and positive")

    total_span = float(t_eval[-1] - t_eval[0])
    max_step_value = total_span / 50.0 if max_step is None else float(max_step)
    if not (np.isfinite(max_step_value) and max_step_value > 0):
        raise ValueError("max_step must be finite and positive")

    nseed = seeds.shape[0]

    def rhs(time: float, flat_state: np.ndarray) -> np.ndarray:
        pts = np.asarray(flat_state, dtype=float).reshape(nseed, 3)
        vel = np.asarray(velocity(pts, float(time)), dtype=float)
        if vel.shape != pts.shape:
            raise ValueError(f"velocity must return shape {pts.shape}, got {vel.shape}")
        if not np.all(np.isfinite(vel)):
            raise ValueError("velocity returned non-finite values")
        return vel.reshape(-1)

    sol = solve_ivp(
        rhs,
        (float(t_eval[0]), float(t_eval[-1])),
        seeds.reshape(-1),
        method=method,
        t_eval=t_eval,
        rtol=rtol,
        atol=atol,
        max_step=max_step_value,
        vectorized=False,
    )
    if not sol.success:
        raise RuntimeError(f"pathline integration failed: {sol.message}")
    if sol.t.shape != t_eval.shape or not np.array_equal(sol.t, t_eval):
        raise RuntimeError("solver did not return every requested sample time")

    positions = sol.y.T.reshape(t_eval.size, nseed, 3)
    if not np.all(np.isfinite(positions)):
        raise RuntimeError("solver returned non-finite positions")

    speeds = np.empty((t_eval.size, nseed), dtype=float)
    for i, time in enumerate(t_eval):
        vel = np.asarray(velocity(positions[i], float(time)), dtype=float)
        if vel.shape != (nseed, 3):
            raise ValueError(f"velocity must return shape {(nseed, 3)}, got {vel.shape}")
        if not np.all(np.isfinite(vel)):
            raise ValueError("velocity returned non-finite values")
        speeds[i] = np.linalg.norm(vel, axis=1)

    return PathlineBundle(
        times=_readonly(t_eval),
        positions=_readonly(positions),
        speeds=_readonly(speeds),
        seed_positions=_readonly(seeds),
        rtol=rtol,
        atol=atol,
        max_step=max_step_value,
        solver_method=str(method),
    )
