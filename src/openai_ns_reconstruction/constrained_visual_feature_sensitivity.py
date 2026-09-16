"""Local sensitivity of visualization fingerprints to bounded candidate parameters.

This module is candidate-agnostic.  It never evaluates a PDE residual and does not
decide visual correspondence.  Its purpose is to tell a basis-growth experiment
which parameter directions can locally move declared visualization features.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Callable, Sequence

import numpy as np


@dataclass(frozen=True)
class VisualFeatureSensitivity:
    parameter_labels: tuple[str, ...]
    feature_labels: tuple[str, ...]
    baseline_features: tuple[float, ...]
    raw_jacobian: tuple[tuple[float, ...], ...]
    scaled_jacobian: tuple[tuple[float, ...], ...]
    singular_values: tuple[float, ...]
    numerical_rank: int
    condition_number: float | None
    parameter_response_norms: tuple[float, ...]
    feature_reachability_norms: tuple[float, ...]
    dead_parameter_labels: tuple[str, ...]
    stencil_by_parameter: tuple[str, ...]
    step_by_parameter: tuple[float, ...]
    truth_boundary: str = (
        "local visualization-feature sensitivity only; not a PDE validation, "
        "not proof of visual correspondence, and not identification of an OpenAI field"
    )

    def to_dict(self) -> dict:
        return asdict(self)


def _as_finite_vector(name: str, value, *, length: int | None = None) -> np.ndarray:
    arr = np.asarray(value, dtype=float)
    if arr.ndim != 1 or arr.size == 0:
        raise ValueError(f"{name} must be a nonempty 1D vector")
    if length is not None and arr.size != length:
        raise ValueError(f"{name} must have length {length}")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} must be finite")
    return arr


def _labels(name: str, labels: Sequence[str] | None, count: int) -> tuple[str, ...]:
    if labels is None:
        return tuple(f"{name}_{i}" for i in range(count))
    result = tuple(str(item) for item in labels)
    if len(result) != count or len(set(result)) != count or any(not item for item in result):
        raise ValueError(f"{name} labels must be unique, nonempty, and length-matched")
    return result


def diagnose_visual_feature_sensitivity(
    feature_fn: Callable[[np.ndarray], np.ndarray],
    parameters,
    *,
    steps,
    lower_bounds,
    upper_bounds,
    feature_scales,
    parameter_labels: Sequence[str] | None = None,
    feature_labels: Sequence[str] | None = None,
    rank_rtol: float = 1e-8,
    dead_parameter_rtol: float = 1e-8,
) -> VisualFeatureSensitivity:
    """Finite-difference a visualization fingerprint with bounded perturbations.

    `feature_fn(theta)` must return a finite 1D feature vector.  `feature_scales`
    are declared positive normalization scales (for example, visualization
    tolerances or target magnitudes).  Each Jacobian column is multiplied by its
    actual finite-difference step and divided by those feature scales before SVD,
    making the rank diagnostic dimensionless without inventing a visual target.

    Central differences are preferred.  A one-sided stencil is used only when a
    requested central perturbation would cross an explicit parameter bound.
    """
    theta = _as_finite_vector("parameters", parameters)
    n = theta.size
    step = _as_finite_vector("steps", steps, length=n)
    lower = _as_finite_vector("lower_bounds", lower_bounds, length=n)
    upper = _as_finite_vector("upper_bounds", upper_bounds, length=n)
    if np.any(step <= 0):
        raise ValueError("steps must be strictly positive")
    if np.any(lower >= upper):
        raise ValueError("every lower bound must be below its upper bound")
    if np.any(theta < lower) or np.any(theta > upper):
        raise ValueError("parameters must lie within bounds")
    if not np.isfinite(rank_rtol) or not 0 < rank_rtol < 1:
        raise ValueError("rank_rtol must lie in (0,1)")
    if not np.isfinite(dead_parameter_rtol) or not 0 <= dead_parameter_rtol < 1:
        raise ValueError("dead_parameter_rtol must lie in [0,1)")

    baseline = _as_finite_vector("feature_fn(parameters)", feature_fn(theta.copy()))
    m = baseline.size
    scales = _as_finite_vector("feature_scales", feature_scales, length=m)
    if np.any(scales <= 0):
        raise ValueError("feature_scales must be strictly positive")

    p_labels = _labels("parameter", parameter_labels, n)
    f_labels = _labels("feature", feature_labels, m)

    jac = np.empty((m, n), dtype=float)
    actual_steps = np.empty(n, dtype=float)
    stencils: list[str] = []

    for j in range(n):
        requested = float(step[j])
        can_minus = theta[j] - requested >= lower[j]
        can_plus = theta[j] + requested <= upper[j]

        if can_minus and can_plus:
            minus = theta.copy()
            plus = theta.copy()
            minus[j] -= requested
            plus[j] += requested
            fm = _as_finite_vector("feature_fn(theta-step)", feature_fn(minus), length=m)
            fp = _as_finite_vector("feature_fn(theta+step)", feature_fn(plus), length=m)
            jac[:, j] = (fp - fm) / (2.0 * requested)
            actual_steps[j] = requested
            stencils.append("central")
            continue

        if can_plus:
            plus = theta.copy()
            plus[j] += requested
            fp = _as_finite_vector("feature_fn(theta+step)", feature_fn(plus), length=m)
            jac[:, j] = (fp - baseline) / requested
            actual_steps[j] = requested
            stencils.append("forward")
            continue

        if can_minus:
            minus = theta.copy()
            minus[j] -= requested
            fm = _as_finite_vector("feature_fn(theta-step)", feature_fn(minus), length=m)
            jac[:, j] = (baseline - fm) / requested
            actual_steps[j] = requested
            stencils.append("backward")
            continue

        plus_room = float(upper[j] - theta[j])
        minus_room = float(theta[j] - lower[j])
        direction = "forward" if plus_room >= minus_room else "backward"
        room = max(plus_room, minus_room)
        if room <= 0:
            raise ValueError(f"parameter {p_labels[j]} has no perturbation room")
        used = min(requested, room)
        if used <= np.finfo(float).eps * max(1.0, abs(theta[j])):
            raise ValueError(f"parameter {p_labels[j]} has numerically zero perturbation room")
        probe = theta.copy()
        if direction == "forward":
            probe[j] += used
            fp = _as_finite_vector("feature_fn(bounded forward)", feature_fn(probe), length=m)
            jac[:, j] = (fp - baseline) / used
        else:
            probe[j] -= used
            fm = _as_finite_vector("feature_fn(bounded backward)", feature_fn(probe), length=m)
            jac[:, j] = (baseline - fm) / used
        actual_steps[j] = used
        stencils.append(direction)

    if not np.all(np.isfinite(jac)):
        raise ValueError("feature Jacobian is non-finite")

    scaled = jac * actual_steps[np.newaxis, :] / scales[:, np.newaxis]
    singular = np.linalg.svd(scaled, compute_uv=False)
    if singular.size == 0 or singular[0] <= 0:
        rank = 0
        condition = None
    else:
        threshold = rank_rtol * singular[0]
        rank = int(np.count_nonzero(singular > threshold))
        condition = None if rank == 0 else float(singular[0] / singular[rank - 1])

    param_norms = np.linalg.norm(scaled, axis=0)
    feature_norms = np.linalg.norm(scaled, axis=1)
    max_param_norm = float(np.max(param_norms)) if param_norms.size else 0.0
    dead_cut = dead_parameter_rtol * max_param_norm
    dead = tuple(
        label for label, norm in zip(p_labels, param_norms)
        if norm <= dead_cut
    )

    return VisualFeatureSensitivity(
        parameter_labels=p_labels,
        feature_labels=f_labels,
        baseline_features=tuple(float(v) for v in baseline),
        raw_jacobian=tuple(tuple(float(v) for v in row) for row in jac),
        scaled_jacobian=tuple(tuple(float(v) for v in row) for row in scaled),
        singular_values=tuple(float(v) for v in singular),
        numerical_rank=rank,
        condition_number=condition,
        parameter_response_norms=tuple(float(v) for v in param_norms),
        feature_reachability_norms=tuple(float(v) for v in feature_norms),
        dead_parameter_labels=dead,
        stencil_by_parameter=tuple(stencils),
        step_by_parameter=tuple(float(v) for v in actual_steps),
    )
