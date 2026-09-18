"""Truth-bounded local swirl-axis orientation diagnostics for public velocity fields.

This module is a visualization/morphology diagnostic only.  It never changes a
candidate velocity, pressure, forcing, acceptance threshold, or validation set.
"""
from __future__ import annotations

import math
from typing import Callable, Iterable, Sequence

import numpy as np

_EVAL_MIN = -2.0
_EVAL_MAX = 2.0
_TIME_MIN = 0.25
_TIME_MAX = 0.75
_EPS = np.finfo(float).eps

_FALSE_TRUTH_STATES = {
    "visualization_ready": False,
    "visual_correspondence_verified": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


def _validate_sha256(value: str) -> str:
    value = str(value).lower()
    if len(value) != 64 or any(ch not in "0123456789abcdef" for ch in value):
        raise ValueError("candidate_sha256 must be 64 lowercase hexadecimal characters")
    return value


def _validate_bounds(bounds: Sequence[Sequence[float]]) -> tuple[tuple[float, float], ...]:
    if len(bounds) != 3:
        raise ValueError("bounds must contain x/y/z intervals")
    out: list[tuple[float, float]] = []
    for pair in bounds:
        if len(pair) != 2:
            raise ValueError("each bound must be a (low, high) pair")
        lo, hi = float(pair[0]), float(pair[1])
        if not (math.isfinite(lo) and math.isfinite(hi) and lo < hi):
            raise ValueError("bounds must be finite and strictly increasing")
        if lo < _EVAL_MIN or hi > _EVAL_MAX:
            raise ValueError("bounds must stay inside the governed [-2,2]^3 evaluation box")
        out.append((lo, hi))
    return tuple(out)


def _validate_axis(axis: Sequence[float]) -> np.ndarray:
    arr = np.asarray(axis, dtype=float)
    if arr.shape != (3,) or not np.isfinite(arr).all():
        raise ValueError("target_axis must be one finite 3-vector")
    norm = float(np.linalg.norm(arr))
    if norm <= 0.0:
        raise ValueError("target_axis must be nonzero")
    return arr / norm


def _weighted_quantile(values: np.ndarray, weights: np.ndarray, quantiles: Iterable[float]) -> list[float]:
    order = np.argsort(values)
    v = values[order]
    w = weights[order]
    cumulative = np.cumsum(w)
    total = float(cumulative[-1])
    if not math.isfinite(total) or total <= 0.0:
        raise ValueError("swirl weights must have positive finite sum")
    return [float(v[np.searchsorted(cumulative, q * total, side="left")]) for q in quantiles]


def _real_eigenvector(vec: np.ndarray) -> np.ndarray | None:
    """Phase-align an eigenvector of a real eigenvalue and return a real unit axis."""
    pivot = int(np.argmax(np.abs(vec)))
    magnitude = float(abs(vec[pivot]))
    if magnitude <= 128.0 * _EPS:
        return None
    aligned = vec * np.exp(-1j * np.angle(vec[pivot]))
    real = np.real(aligned)
    imag_norm = float(np.linalg.norm(np.imag(aligned)))
    real_norm = float(np.linalg.norm(real))
    if real_norm <= 128.0 * _EPS:
        return None
    # A real matrix with a real eigenvalue admits a real eigenvector.  If the
    # returned vector is not numerically close to that situation, skip the
    # point rather than manufacturing an orientation.
    if imag_norm > 4096.0 * _EPS * max(1.0, real_norm):
        return None
    return real / real_norm


def diagnose_swirl_axis_orientation(
    velocity: Callable[[np.ndarray, float], np.ndarray],
    *,
    time: float,
    bounds: Sequence[Sequence[float]] = ((-1.0, 1.0), (-1.0, 1.0), (-1.0, 1.0)),
    grid_size: int = 13,
    target_axis: Sequence[float] = (0.0, 0.0, 1.0),
    candidate_sha256: str,
    provenance: str,
) -> dict[str, object]:
    """Measure local swirl-axis orientation from ``velocity(points,time)``.

    ``lambda_ci`` is taken from the imaginary part of the complex-conjugate
    eigenvalue pair of the independently finite-differenced velocity-gradient
    tensor.  At swirl-active points, the eigenvector belonging to the remaining
    real eigenvalue is used as an *unsigned* local axis.  Unsigned orientation
    avoids laundering arbitrary eigenvector signs into a handedness claim.

    The machine-scale complex-eigenvalue guard is numerical noise suppression,
    not a candidate acceptance threshold.  No visual pass/fail threshold is
    defined here.
    """
    if not callable(velocity):
        raise TypeError("velocity must be callable as velocity(points,time)")
    t = float(time)
    if not math.isfinite(t) or not (_TIME_MIN <= t <= _TIME_MAX):
        raise ValueError("time must lie in the governed [0.25,0.75] interval")
    if not isinstance(grid_size, (int, np.integer)) or int(grid_size) < 7 or int(grid_size) % 2 == 0:
        raise ValueError("grid_size must be an odd integer >= 7")
    n = int(grid_size)
    box = _validate_bounds(bounds)
    axis = _validate_axis(target_axis)
    sha = _validate_sha256(candidate_sha256)
    provenance = str(provenance).strip()
    if not provenance:
        raise ValueError("provenance must be a nonempty source/identity string")

    coords = [np.linspace(lo, hi, n, dtype=float) for lo, hi in box]
    xx, yy, zz = np.meshgrid(*coords, indexing="ij")
    points = np.column_stack((xx.ravel(), yy.ravel(), zz.ravel()))
    sampled = np.asarray(velocity(points, t), dtype=float)
    if sampled.shape != (points.shape[0], 3):
        raise ValueError("velocity must return finite values with shape (N,3)")
    if not np.isfinite(sampled).all():
        raise ValueError("velocity returned NaN or Inf")
    field = sampled.reshape((n, n, n, 3))
    interior = field[1:-1, 1:-1, 1:-1, :]
    velocity_rms = float(np.sqrt(np.mean(np.sum(interior * interior, axis=-1))))
    if velocity_rms == 0.0:
        raise ValueError("exact-zero interior velocity is not a valid visualization diagnostic")

    dx, dy, dz = (float(c[1] - c[0]) for c in coords)
    deriv_x = (field[2:, 1:-1, 1:-1, :] - field[:-2, 1:-1, 1:-1, :]) / (2.0 * dx)
    deriv_y = (field[1:-1, 2:, 1:-1, :] - field[1:-1, :-2, 1:-1, :]) / (2.0 * dy)
    deriv_z = (field[1:-1, 1:-1, 2:, :] - field[1:-1, 1:-1, :-2, :]) / (2.0 * dz)
    grad = np.stack((deriv_x, deriv_y, deriv_z), axis=-1)  # (..., component, coordinate)
    flat_grad = grad.reshape((-1, 3, 3))
    divergence = np.trace(flat_grad, axis1=1, axis2=2)
    divergence_rms = float(np.sqrt(np.mean(divergence * divergence)))

    eigenvalues, eigenvectors = np.linalg.eig(flat_grad)
    grad_scale = np.linalg.norm(flat_grad, axis=(1, 2))
    imag_guard = 64.0 * _EPS * np.maximum(1.0, grad_scale)
    lambda_ci_all = np.max(np.abs(np.imag(eigenvalues)), axis=1)

    weights: list[float] = []
    axes: list[np.ndarray] = []
    for i in range(flat_grad.shape[0]):
        if not (lambda_ci_all[i] > imag_guard[i]):
            continue
        complex_mask = np.abs(np.imag(eigenvalues[i])) > imag_guard[i]
        if int(np.count_nonzero(complex_mask)) != 2:
            continue
        real_indices = np.flatnonzero(~complex_mask)
        if real_indices.size != 1:
            continue
        local_axis = _real_eigenvector(eigenvectors[i, :, int(real_indices[0])])
        if local_axis is None:
            continue
        weights.append(float(lambda_ci_all[i]))
        axes.append(local_axis)

    interior_count = int(flat_grad.shape[0])
    active_count = len(weights)
    active_fraction = float(active_count / interior_count)
    lambda_ci_rms = float(np.sqrt(np.mean(lambda_ci_all * lambda_ci_all)))
    lambda_ci_max = float(np.max(lambda_ci_all))

    base: dict[str, object] = {
        "schema": "cr_a9_swirl_axis_orientation_v1",
        "time": t,
        "grid_size": n,
        "bounds": box,
        "target_axis": tuple(float(x) for x in axis),
        "candidate_sha256": sha,
        "provenance": provenance,
        "velocity_rms": velocity_rms,
        "divergence_rms_same_operator": divergence_rms,
        "lambda_ci_rms": lambda_ci_rms,
        "lambda_ci_max": lambda_ci_max,
        "swirl_active_points": active_count,
        "interior_points": interior_count,
        "swirl_active_fraction": active_fraction,
        "complex_eigenvalue_guard": "64*eps*max(1,||grad_u||_F); numerical noise guard only",
        "truth_states": dict(_FALSE_TRUTH_STATES),
    }

    if active_count == 0:
        base.update(
            {
                "swirl_axis_detected": False,
                "weighted_mean_abs_axis_alignment": None,
                "weighted_mean_tilt_degrees": None,
                "tilt_degrees_q10_q50_q90": None,
                "principal_unsigned_axis": None,
                "principal_axis_abs_alignment": None,
                "orientation_concentration": None,
            }
        )
        return base

    w = np.asarray(weights, dtype=float)
    a = np.asarray(axes, dtype=float)
    alignment = np.clip(np.abs(a @ axis), 0.0, 1.0)
    tilt = np.degrees(np.arccos(alignment))
    total_weight = float(np.sum(w))
    weighted_alignment = float(np.sum(w * alignment) / total_weight)
    weighted_tilt = float(np.sum(w * tilt) / total_weight)
    tilt_quantiles = _weighted_quantile(tilt, w, (0.10, 0.50, 0.90))

    tensor = np.einsum("n,ni,nj->ij", w, a, a) / total_weight
    evals, evecs = np.linalg.eigh(tensor)
    principal = evecs[:, int(np.argmax(evals))]
    if float(np.dot(principal, axis)) < 0.0:
        principal = -principal
    principal_alignment = float(np.clip(abs(np.dot(principal, axis)), 0.0, 1.0))
    concentration = float(np.max(evals))  # 1/3 isotropic axes; 1 perfectly aligned axes.

    base.update(
        {
            "swirl_axis_detected": True,
            "weighted_mean_abs_axis_alignment": weighted_alignment,
            "weighted_mean_tilt_degrees": weighted_tilt,
            "tilt_degrees_q10_q50_q90": tuple(tilt_quantiles),
            "principal_unsigned_axis": tuple(float(x) for x in principal),
            "principal_axis_abs_alignment": principal_alignment,
            "orientation_concentration": concentration,
        }
    )
    return base
