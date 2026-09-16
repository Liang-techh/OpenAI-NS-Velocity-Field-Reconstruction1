"""Temporal coefficient-lift capacity diagnostic for constrained basis growth.

This module compares two local linearized parameterizations across multiple
reference times:

1. static coefficients shared by every time, and
2. the same static coefficients plus one declared affine-in-time lift
   direction ``tau(t) * slope``.

It is intended to answer a narrow representation question before growing the
candidate family: does a single time-dependent basis channel materially improve
both residual-side and visualization-side reachability compared with a static
basis?  It does not construct or optimize a velocity field and it does not
evaluate a PDE operator itself.

All inputs are precomputed local response matrices. Results are capacity
diagnostics only, not independent PDE validation or proof of visual
correspondence.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Sequence

import numpy as np
from scipy.optimize import lsq_linear


@dataclass(frozen=True)
class TemporalLiftLevel:
    model: str
    parameter_count: int
    coefficients: tuple[float, ...]
    residual_remaining_ratio: float
    visual_remaining_ratio: float
    joint_remaining_ratio: float
    residual_recovered_fraction: float
    visual_recovered_fraction: float
    numerical_rank: int
    condition_number: float | None
    active_lower_bounds: tuple[str, ...]
    active_upper_bounds: tuple[str, ...]


@dataclass(frozen=True)
class TemporalLiftCapacity:
    times: tuple[float, ...]
    tau: tuple[float, ...]
    base_parameter_labels: tuple[str, ...]
    lifted_parameter_label: str
    static: TemporalLiftLevel
    affine_lift: TemporalLiftLevel
    residual_gain_from_lift: float
    visual_gain_from_lift: float
    joint_gain_from_lift: float
    lift_materially_improves_both: bool
    claim_scope: str = "local_expression_capacity_only"
    pde_validated: bool = False
    visual_correspondence_verified: bool = False
    paper_exact: bool = False
    openai_field_identified: bool = False

    def to_dict(self) -> dict:
        return asdict(self)


def _vector(name: str, value, *, length: int | None = None) -> np.ndarray:
    arr = np.asarray(value, dtype=float)
    if arr.ndim != 1 or arr.size == 0:
        raise ValueError(f"{name} must be a nonempty 1D vector")
    if length is not None and arr.size != length:
        raise ValueError(f"{name} must have length {length}")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} must be finite")
    return arr


def _cube(name: str, value, *, time_count: int | None = None) -> np.ndarray:
    arr = np.asarray(value, dtype=float)
    if arr.ndim != 3 or min(arr.shape) == 0:
        raise ValueError(f"{name} must have shape (time, rows, parameters)")
    if time_count is not None and arr.shape[0] != time_count:
        raise ValueError(f"{name} has the wrong time dimension")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} must be finite")
    return arr


def _targets(name: str, value, *, time_count: int, rows: int) -> np.ndarray:
    arr = np.asarray(value, dtype=float)
    if arr.shape != (time_count, rows):
        raise ValueError(f"{name} must have shape {(time_count, rows)}")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} must be finite")
    return arr


def _labels(labels: Sequence[str] | None, count: int) -> tuple[str, ...]:
    if labels is None:
        return tuple(f"basis_{i}" for i in range(count))
    out = tuple(str(x) for x in labels)
    if len(out) != count or len(set(out)) != count or any(not x for x in out):
        raise ValueError("parameter_labels must be unique, nonempty, and length matched")
    return out


def _rank_condition(matrix: np.ndarray, rtol: float) -> tuple[int, float | None]:
    s = np.linalg.svd(matrix, compute_uv=False)
    if s.size == 0 or s[0] <= 0:
        return 0, None
    rank = int(np.count_nonzero(s > rtol * s[0]))
    if rank == 0:
        return 0, None
    return rank, float(s[0] / s[rank - 1])


def _solve_level(
    A_r: np.ndarray,
    b_r: np.ndarray,
    A_v: np.ndarray,
    b_v: np.ndarray,
    *,
    lower: np.ndarray,
    upper: np.ndarray,
    labels: tuple[str, ...],
    residual_weight: float,
    visual_weight: float,
    residual_norm: float,
    visual_norm: float,
    joint_norm: float,
    rank_rtol: float,
    model: str,
) -> TemporalLiftLevel:
    A = np.vstack((np.sqrt(residual_weight) * A_r, np.sqrt(visual_weight) * A_v))
    b = np.concatenate((np.sqrt(residual_weight) * b_r, np.sqrt(visual_weight) * b_v))
    result = lsq_linear(A, b, bounds=(lower, upper), method="trf", lsmr_tol="auto")
    if not result.success or not np.all(np.isfinite(result.x)):
        raise RuntimeError(f"{model} bounded solve failed")
    coeff = np.asarray(result.x, dtype=float)

    rem_r = b_r - A_r @ coeff
    rem_v = b_v - A_v @ coeff
    rr = float(np.linalg.norm(rem_r) / residual_norm)
    vr = float(np.linalg.norm(rem_v) / visual_norm)
    jr = float(
        np.sqrt(
            residual_weight * np.dot(rem_r, rem_r)
            + visual_weight * np.dot(rem_v, rem_v)
        )
        / joint_norm
    )
    rank, cond = _rank_condition(A, rank_rtol)

    scale = np.maximum(1.0, np.maximum(np.abs(lower), np.abs(upper)))
    tol = 128.0 * np.finfo(float).eps * scale
    active_lo = tuple(labels[i] for i in range(len(labels)) if abs(coeff[i] - lower[i]) <= tol[i])
    active_hi = tuple(labels[i] for i in range(len(labels)) if abs(coeff[i] - upper[i]) <= tol[i])

    return TemporalLiftLevel(
        model=model,
        parameter_count=len(labels),
        coefficients=tuple(float(x) for x in coeff),
        residual_remaining_ratio=rr,
        visual_remaining_ratio=vr,
        joint_remaining_ratio=jr,
        residual_recovered_fraction=float(1.0 - rr),
        visual_recovered_fraction=float(1.0 - vr),
        numerical_rank=rank,
        condition_number=cond,
        active_lower_bounds=active_lo,
        active_upper_bounds=active_hi,
    )


def diagnose_affine_temporal_lift(
    times,
    residual_jacobians,
    residual_correction_targets,
    visual_jacobians,
    visual_correction_targets,
    *,
    residual_scales,
    visual_scales,
    lower_bounds,
    upper_bounds,
    lifted_parameter_index: int,
    slope_bound: float,
    parameter_labels: Sequence[str] | None = None,
    residual_weight: float = 1.0,
    visual_weight: float = 1.0,
    rank_rtol: float = 1e-8,
    material_gain: float = 0.05,
) -> TemporalLiftCapacity:
    """Compare a static basis against one affine-in-time coefficient lift.

    The base Jacobians must have shape ``(T, rows, P)`` and use the same
    parameter ordering at every time. The affine model augments the shared base
    update ``delta`` with one scalar ``s`` on a declared parameter direction:

        delta(t) = delta + tau(t) * s * e_k,

    where ``tau`` is centered and scaled from the supplied times to [-1,1].
    This tests the smallest possible time-dependent growth of an existing basis
    direction without introducing a new spatial shape.

    Bounds for the shared base coefficients are caller-declared. ``slope_bound``
    is symmetric and must be positive. The routine is a local response/capacity
    diagnostic and does not assert that the resulting update is globally valid.
    """
    t = _vector("times", times)
    if t.size < 3 or np.any(np.diff(t) <= 0.0):
        raise ValueError("times must contain at least three strictly increasing values")
    center = 0.5 * (t[0] + t[-1])
    halfspan = 0.5 * (t[-1] - t[0])
    if halfspan <= 0.0:
        raise ValueError("time interval must have positive width")
    tau = (t - center) / halfspan

    jr = _cube("residual_jacobians", residual_jacobians, time_count=t.size)
    jv = _cube("visual_jacobians", visual_jacobians, time_count=t.size)
    if jr.shape[2] != jv.shape[2]:
        raise ValueError("residual and visual Jacobians must share parameter columns")
    p = jr.shape[2]
    rt = _targets(
        "residual_correction_targets",
        residual_correction_targets,
        time_count=t.size,
        rows=jr.shape[1],
    )
    vt = _targets(
        "visual_correction_targets",
        visual_correction_targets,
        time_count=t.size,
        rows=jv.shape[1],
    )

    rs = _vector("residual_scales", residual_scales, length=jr.shape[1])
    vs = _vector("visual_scales", visual_scales, length=jv.shape[1])
    lo = _vector("lower_bounds", lower_bounds, length=p)
    hi = _vector("upper_bounds", upper_bounds, length=p)
    if np.any(rs <= 0.0) or np.any(vs <= 0.0):
        raise ValueError("row scales must be strictly positive")
    if np.any(lo >= hi):
        raise ValueError("every lower bound must be below its upper bound")
    if not isinstance(lifted_parameter_index, (int, np.integer)) or not (0 <= int(lifted_parameter_index) < p):
        raise ValueError("lifted_parameter_index is out of range")
    k = int(lifted_parameter_index)
    if not np.isfinite(slope_bound) or slope_bound <= 0.0:
        raise ValueError("slope_bound must be positive and finite")
    if not np.isfinite(residual_weight) or residual_weight <= 0.0:
        raise ValueError("residual_weight must be positive and finite")
    if not np.isfinite(visual_weight) or visual_weight <= 0.0:
        raise ValueError("visual_weight must be positive and finite")
    if not np.isfinite(rank_rtol) or not 0.0 < rank_rtol < 1.0:
        raise ValueError("rank_rtol must lie in (0,1)")
    if not np.isfinite(material_gain) or not 0.0 <= material_gain < 1.0:
        raise ValueError("material_gain must lie in [0,1)")

    labels = _labels(parameter_labels, p)
    lift_label = f"{labels[k]}:time_slope"

    jr_s = jr / rs[np.newaxis, :, np.newaxis]
    jv_s = jv / vs[np.newaxis, :, np.newaxis]
    rt_s = rt / rs[np.newaxis, :]
    vt_s = vt / vs[np.newaxis, :]

    A_r_static = jr_s.reshape(-1, p)
    A_v_static = jv_s.reshape(-1, p)
    b_r = rt_s.reshape(-1)
    b_v = vt_s.reshape(-1)

    residual_norm = float(np.linalg.norm(b_r))
    visual_norm = float(np.linalg.norm(b_v))
    if residual_norm <= np.finfo(float).eps or visual_norm <= np.finfo(float).eps:
        raise ValueError("scaled correction targets must both be nonzero")
    joint_norm = float(
        np.sqrt(residual_weight * residual_norm**2 + visual_weight * visual_norm**2)
    )

    static = _solve_level(
        A_r_static,
        b_r,
        A_v_static,
        b_v,
        lower=lo,
        upper=hi,
        labels=labels,
        residual_weight=residual_weight,
        visual_weight=visual_weight,
        residual_norm=residual_norm,
        visual_norm=visual_norm,
        joint_norm=joint_norm,
        rank_rtol=rank_rtol,
        model="static_shared_coefficients",
    )

    temporal_r = (jr_s[:, :, k] * tau[:, np.newaxis]).reshape(-1, 1)
    temporal_v = (jv_s[:, :, k] * tau[:, np.newaxis]).reshape(-1, 1)
    A_r_affine = np.hstack((A_r_static, temporal_r))
    A_v_affine = np.hstack((A_v_static, temporal_v))
    lo_affine = np.concatenate((lo, [-float(slope_bound)]))
    hi_affine = np.concatenate((hi, [float(slope_bound)]))
    labels_affine = labels + (lift_label,)

    affine = _solve_level(
        A_r_affine,
        b_r,
        A_v_affine,
        b_v,
        lower=lo_affine,
        upper=hi_affine,
        labels=labels_affine,
        residual_weight=residual_weight,
        visual_weight=visual_weight,
        residual_norm=residual_norm,
        visual_norm=visual_norm,
        joint_norm=joint_norm,
        rank_rtol=rank_rtol,
        model="single_affine_temporal_lift",
    )

    residual_gain = float(static.residual_remaining_ratio - affine.residual_remaining_ratio)
    visual_gain = float(static.visual_remaining_ratio - affine.visual_remaining_ratio)
    joint_gain = float(static.joint_remaining_ratio - affine.joint_remaining_ratio)
    material = residual_gain >= material_gain and visual_gain >= material_gain and joint_gain > 0.0

    return TemporalLiftCapacity(
        times=tuple(float(x) for x in t),
        tau=tuple(float(x) for x in tau),
        base_parameter_labels=labels,
        lifted_parameter_label=lift_label,
        static=static,
        affine_lift=affine,
        residual_gain_from_lift=residual_gain,
        visual_gain_from_lift=visual_gain,
        joint_gain_from_lift=joint_gain,
        lift_materially_improves_both=material,
    )
