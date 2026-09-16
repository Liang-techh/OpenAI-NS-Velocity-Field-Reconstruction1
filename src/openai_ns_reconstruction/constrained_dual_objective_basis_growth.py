"""Bounded dual-objective diagnostics for nested basis growth.

This module does not construct a velocity field and does not evaluate a PDE
operator.  It consumes *declared local response matrices* for two separate
objectives:

1. an optimizer/PDE-side residual correction target, and
2. a visualization-feature correction target.

For each nested basis dimension it solves one bounded linearized coefficient
update against both objectives simultaneously.  Using one shared coefficient
vector is important: a basis direction that looks useful when residual and
visual objectives are fitted independently may be mutually incompatible in the
actual candidate parameterization.

All results are local capacity diagnostics.  They are not optimizer convergence,
PDE validation, visual correspondence, or identification of an OpenAI field.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Sequence

import numpy as np
from scipy.optimize import lsq_linear


@dataclass(frozen=True)
class NestedBasisLevel:
    dimension: int
    parameter_labels: tuple[str, ...]
    coefficients: tuple[float, ...]
    residual_rank: int
    visual_rank: int
    joint_rank: int
    residual_condition: float | None
    visual_condition: float | None
    joint_condition: float | None
    residual_remaining_ratio: float
    visual_remaining_ratio: float
    joint_remaining_ratio: float
    residual_recovered_fraction: float
    visual_recovered_fraction: float
    active_lower_bounds: tuple[str, ...]
    active_upper_bounds: tuple[str, ...]
    jointly_improves_both: bool


@dataclass(frozen=True)
class DualObjectiveBasisGrowth:
    parameter_labels: tuple[str, ...]
    nested_dimensions: tuple[int, ...]
    residual_target_norm: float
    visual_target_norm: float
    residual_weight: float
    visual_weight: float
    levels: tuple[NestedBasisLevel, ...]
    smallest_jointly_improving_dimension: int | None
    truth_boundary: str = (
        "bounded local basis-capacity diagnostic only; not optimizer convergence, "
        "not independent PDE validation, not proof of visual correspondence, "
        "and not identification of an OpenAI field"
    )

    def to_dict(self) -> dict:
        return asdict(self)


def _finite_vector(name: str, value, *, length: int | None = None) -> np.ndarray:
    arr = np.asarray(value, dtype=float)
    if arr.ndim != 1 or arr.size == 0:
        raise ValueError(f"{name} must be a nonempty 1D vector")
    if length is not None and arr.size != length:
        raise ValueError(f"{name} must have length {length}")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} must be finite")
    return arr


def _finite_matrix(name: str, value) -> np.ndarray:
    arr = np.asarray(value, dtype=float)
    if arr.ndim != 2 or min(arr.shape) == 0:
        raise ValueError(f"{name} must be a nonempty 2D matrix")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} must be finite")
    return arr


def _labels(labels: Sequence[str] | None, count: int) -> tuple[str, ...]:
    if labels is None:
        return tuple(f"basis_{i}" for i in range(count))
    result = tuple(str(item) for item in labels)
    if len(result) != count or len(set(result)) != count or any(not item for item in result):
        raise ValueError("parameter_labels must be unique, nonempty, and length-matched")
    return result


def _rank_condition(matrix: np.ndarray, rank_rtol: float) -> tuple[int, float | None]:
    singular = np.linalg.svd(matrix, compute_uv=False)
    if singular.size == 0 or singular[0] <= 0.0:
        return 0, None
    threshold = rank_rtol * singular[0]
    rank = int(np.count_nonzero(singular > threshold))
    if rank == 0:
        return 0, None
    return rank, float(singular[0] / singular[rank - 1])


def diagnose_nested_dual_objective_basis_growth(
    residual_jacobian,
    residual_correction_target,
    visual_jacobian,
    visual_correction_target,
    *,
    residual_scales,
    visual_scales,
    lower_bounds,
    upper_bounds,
    nested_dimensions: Sequence[int],
    parameter_labels: Sequence[str] | None = None,
    residual_weight: float = 1.0,
    visual_weight: float = 1.0,
    rank_rtol: float = 1e-8,
    improvement_rtol: float = 1e-8,
) -> DualObjectiveBasisGrowth:
    """Audit nested basis dimensions with one shared bounded local update.

    The Jacobians may come from independent diagnostics, but both must use the
    same parameter ordering.  Targets are desired *corrections*: for example, a
    residual-side caller may pass ``-current_residual`` and a visualization-side
    caller may pass ``target_features-current_features``.

    ``residual_scales`` and ``visual_scales`` are explicit positive row scales.
    They keep heterogeneous units from silently deciding the joint least-squares
    result.  No default acceptance threshold is inferred from these scales.

    Bounds apply to the local coefficient update itself.  They must therefore be
    chosen from the candidate's declared local perturbation budget rather than
    invented after seeing the result.
    """
    jr = _finite_matrix("residual_jacobian", residual_jacobian)
    jv = _finite_matrix("visual_jacobian", visual_jacobian)
    if jr.shape[1] != jv.shape[1]:
        raise ValueError("residual and visual Jacobians must share parameter columns")
    parameter_count = jr.shape[1]

    rt = _finite_vector(
        "residual_correction_target", residual_correction_target, length=jr.shape[0]
    )
    vt = _finite_vector(
        "visual_correction_target", visual_correction_target, length=jv.shape[0]
    )
    rs = _finite_vector("residual_scales", residual_scales, length=jr.shape[0])
    vs = _finite_vector("visual_scales", visual_scales, length=jv.shape[0])
    lo = _finite_vector("lower_bounds", lower_bounds, length=parameter_count)
    hi = _finite_vector("upper_bounds", upper_bounds, length=parameter_count)

    if np.any(rs <= 0.0) or np.any(vs <= 0.0):
        raise ValueError("row scales must be strictly positive")
    if np.any(lo >= hi):
        raise ValueError("every lower bound must be below its upper bound")
    if not np.isfinite(residual_weight) or residual_weight <= 0.0:
        raise ValueError("residual_weight must be positive and finite")
    if not np.isfinite(visual_weight) or visual_weight <= 0.0:
        raise ValueError("visual_weight must be positive and finite")
    if not np.isfinite(rank_rtol) or not 0.0 < rank_rtol < 1.0:
        raise ValueError("rank_rtol must lie in (0,1)")
    if not np.isfinite(improvement_rtol) or not 0.0 <= improvement_rtol < 1.0:
        raise ValueError("improvement_rtol must lie in [0,1)")

    dims = tuple(int(value) for value in nested_dimensions)
    if not dims:
        raise ValueError("nested_dimensions must be nonempty")
    if any(value <= 0 or value > parameter_count for value in dims):
        raise ValueError("nested dimensions must lie in [1, parameter_count]")
    if any(a >= b for a, b in zip(dims, dims[1:])):
        raise ValueError("nested dimensions must be strictly increasing")

    labels = _labels(parameter_labels, parameter_count)

    jr_scaled = jr / rs[:, np.newaxis]
    jv_scaled = jv / vs[:, np.newaxis]
    rt_scaled = rt / rs
    vt_scaled = vt / vs

    residual_target_norm = float(np.linalg.norm(rt_scaled))
    visual_target_norm = float(np.linalg.norm(vt_scaled))
    if residual_target_norm <= np.finfo(float).eps:
        raise ValueError("scaled residual correction target must be nonzero")
    if visual_target_norm <= np.finfo(float).eps:
        raise ValueError("scaled visual correction target must be nonzero")

    base_joint_sq = (
        residual_weight * residual_target_norm**2
        + visual_weight * visual_target_norm**2
    )
    base_joint_norm = float(np.sqrt(base_joint_sq))

    levels: list[NestedBasisLevel] = []
    first_joint: int | None = None

    for dimension in dims:
        ar = jr_scaled[:, :dimension]
        av = jv_scaled[:, :dimension]
        joint_a = np.vstack(
            (
                np.sqrt(residual_weight) * ar,
                np.sqrt(visual_weight) * av,
            )
        )
        joint_b = np.concatenate(
            (
                np.sqrt(residual_weight) * rt_scaled,
                np.sqrt(visual_weight) * vt_scaled,
            )
        )

        result = lsq_linear(
            joint_a,
            joint_b,
            bounds=(lo[:dimension], hi[:dimension]),
            method="trf",
            lsmr_tol="auto",
        )
        if not result.success or not np.all(np.isfinite(result.x)):
            raise RuntimeError(f"bounded local solve failed at dimension {dimension}")

        coeff = np.asarray(result.x, dtype=float)
        residual_remaining = rt_scaled - ar @ coeff
        visual_remaining = vt_scaled - av @ coeff

        residual_ratio = float(np.linalg.norm(residual_remaining) / residual_target_norm)
        visual_ratio = float(np.linalg.norm(visual_remaining) / visual_target_norm)
        joint_remaining = float(
            np.sqrt(
                residual_weight * np.dot(residual_remaining, residual_remaining)
                + visual_weight * np.dot(visual_remaining, visual_remaining)
            )
            / base_joint_norm
        )

        residual_rank, residual_condition = _rank_condition(ar, rank_rtol)
        visual_rank, visual_condition = _rank_condition(av, rank_rtol)
        joint_rank, joint_condition = _rank_condition(joint_a, rank_rtol)

        scale = np.maximum(1.0, np.maximum(np.abs(lo[:dimension]), np.abs(hi[:dimension])))
        bound_tol = 128.0 * np.finfo(float).eps * scale
        lower_active = tuple(
            labels[index]
            for index in range(dimension)
            if abs(coeff[index] - lo[index]) <= bound_tol[index]
        )
        upper_active = tuple(
            labels[index]
            for index in range(dimension)
            if abs(coeff[index] - hi[index]) <= bound_tol[index]
        )

        jointly_improves = (
            residual_ratio < 1.0 - improvement_rtol
            and visual_ratio < 1.0 - improvement_rtol
        )
        if jointly_improves and first_joint is None:
            first_joint = dimension

        levels.append(
            NestedBasisLevel(
                dimension=dimension,
                parameter_labels=labels[:dimension],
                coefficients=tuple(float(value) for value in coeff),
                residual_rank=residual_rank,
                visual_rank=visual_rank,
                joint_rank=joint_rank,
                residual_condition=residual_condition,
                visual_condition=visual_condition,
                joint_condition=joint_condition,
                residual_remaining_ratio=residual_ratio,
                visual_remaining_ratio=visual_ratio,
                joint_remaining_ratio=joint_remaining,
                residual_recovered_fraction=float(1.0 - residual_ratio),
                visual_recovered_fraction=float(1.0 - visual_ratio),
                active_lower_bounds=lower_active,
                active_upper_bounds=upper_active,
                jointly_improves_both=jointly_improves,
            )
        )

    return DualObjectiveBasisGrowth(
        parameter_labels=labels,
        nested_dimensions=dims,
        residual_target_norm=residual_target_norm,
        visual_target_norm=visual_target_norm,
        residual_weight=float(residual_weight),
        visual_weight=float(visual_weight),
        levels=tuple(levels),
        smallest_jointly_improving_dimension=first_joint,
    )
