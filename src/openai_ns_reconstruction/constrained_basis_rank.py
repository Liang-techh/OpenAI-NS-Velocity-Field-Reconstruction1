"""Rank-revealing screening for linearized constrained-reconstruction bases.

This module is a *capacity diagnostic*. It does not optimize the Navier--Stokes
candidate, evaluate an independent PDE acceptance metric, alter forcing, or
change preregistered thresholds. In particular, the projected residual reported
here is an unbounded first-order least-squares diagnostic; it is not a feasible
nonlinear update and must not be reported as a validated PDE residual.

External method provenance
--------------------------
The column selection uses the public ``scipy.linalg.qr(..., pivoting=True)`` API
(rank-revealing QR / QR with column pivoting). No SciPy implementation code is
copied. Source screened: scipy/scipy commit
``eff78058ed0d4cb8e1f0b5c5585d62d712e948d2`` (BSD-3-Clause).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np
from scipy.linalg import qr


@dataclass(frozen=True)
class LinearizedBasisScreen:
    """Immutable summary of one linearized basis-capacity screen."""

    residual_size: int
    basis_size: int
    numerical_rank: int
    independent_indices: tuple[int, ...]
    dependent_indices: tuple[int, ...]
    pivot_order: tuple[int, ...]
    independent_labels: tuple[str, ...]
    dependent_labels: tuple[str, ...]
    column_norms: tuple[float, ...]
    qr_diagonal: tuple[float, ...]
    rank_threshold: float
    normalized_condition_number: float
    baseline_l2: float
    projected_l2: float
    projected_reduction_fraction: float


def _as_finite_vector(name: str, value: np.ndarray | Sequence[float]) -> np.ndarray:
    array = np.asarray(value, dtype=float)
    if array.ndim != 1:
        raise ValueError(f"{name} must be a one-dimensional array")
    if array.size == 0:
        raise ValueError(f"{name} must not be empty")
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} contains non-finite values")
    return array


def _as_finite_matrix(name: str, value: np.ndarray | Sequence[Sequence[float]]) -> np.ndarray:
    array = np.asarray(value, dtype=float)
    if array.ndim != 2:
        raise ValueError(f"{name} must be a two-dimensional array")
    if array.shape[0] == 0 or array.shape[1] == 0:
        raise ValueError(f"{name} must have nonzero row and column counts")
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} contains non-finite values")
    return array


def screen_linearized_residual_capacity(
    residual: np.ndarray | Sequence[float],
    jacobian: np.ndarray | Sequence[Sequence[float]],
    *,
    labels: Sequence[str] | None = None,
    relative_tolerance: float | None = None,
    absolute_tolerance: float = 0.0,
) -> LinearizedBasisScreen:
    """Screen effective basis rank and the first-order residual subspace.

    ``jacobian[:, j]`` is the residual derivative associated with basis or
    parameter direction ``j``. Columns are normalized before pivoted QR so
    selection diagnoses geometric redundancy rather than merely preferring a
    direction because its coefficient uses a larger numerical scale.

    The returned ``projected_l2`` is the residual norm after an *unbounded*
    linear least-squares projection onto the selected column space. It is only
    a local capacity indicator. Bounds, nonlinear terms, nontriviality,
    boundary/support constraints, forcing restrictions and held-out validation
    still have to be enforced by their owning components.
    """

    r = _as_finite_vector("residual", residual)
    j = _as_finite_matrix("jacobian", jacobian)

    if j.shape[0] != r.size:
        raise ValueError("jacobian row count must equal residual length")

    if absolute_tolerance < 0 or not np.isfinite(absolute_tolerance):
        raise ValueError("absolute_tolerance must be finite and nonnegative")

    if relative_tolerance is None:
        relative_tolerance = np.finfo(float).eps * max(j.shape)
    if relative_tolerance < 0 or not np.isfinite(relative_tolerance):
        raise ValueError("relative_tolerance must be finite and nonnegative")

    n_basis = j.shape[1]
    if labels is None:
        label_tuple = tuple(f"basis_{index}" for index in range(n_basis))
    else:
        label_tuple = tuple(str(label) for label in labels)
        if len(label_tuple) != n_basis:
            raise ValueError("labels length must equal jacobian column count")
        if len(set(label_tuple)) != len(label_tuple):
            raise ValueError("labels must be unique")

    column_norms = np.linalg.norm(j, axis=0)
    normalized = np.zeros_like(j)
    nonzero = column_norms > absolute_tolerance
    normalized[:, nonzero] = j[:, nonzero] / column_norms[nonzero]

    _, r_factor, pivots = qr(
        normalized,
        mode="economic",
        pivoting=True,
        check_finite=False,
    )
    diagonal = np.abs(np.diag(r_factor))

    leading = float(diagonal[0]) if diagonal.size else 0.0
    rank_threshold = max(float(absolute_tolerance), float(relative_tolerance) * leading)
    rank = int(np.count_nonzero(diagonal > rank_threshold))

    independent = tuple(int(index) for index in pivots[:rank])
    independent_set = set(independent)
    dependent = tuple(index for index in range(n_basis) if index not in independent_set)

    baseline_l2 = float(np.linalg.norm(r))

    if rank == 0:
        projected_l2 = baseline_l2
        normalized_condition_number = float("inf")
    else:
        selected = j[:, independent]
        step, *_ = np.linalg.lstsq(selected, -r, rcond=None)
        projected = r + selected @ step
        projected_l2 = float(np.linalg.norm(projected))
        normalized_condition_number = float(np.linalg.cond(normalized[:, independent]))

    if baseline_l2 == 0.0:
        reduction_fraction = 0.0
    else:
        reduction_fraction = (baseline_l2 - projected_l2) / baseline_l2
        reduction_fraction = float(np.clip(reduction_fraction, 0.0, 1.0))

    return LinearizedBasisScreen(
        residual_size=r.size,
        basis_size=n_basis,
        numerical_rank=rank,
        independent_indices=independent,
        dependent_indices=dependent,
        pivot_order=tuple(int(index) for index in pivots),
        independent_labels=tuple(label_tuple[index] for index in independent),
        dependent_labels=tuple(label_tuple[index] for index in dependent),
        column_norms=tuple(float(value) for value in column_norms),
        qr_diagonal=tuple(float(value) for value in diagonal),
        rank_threshold=float(rank_threshold),
        normalized_condition_number=normalized_condition_number,
        baseline_l2=baseline_l2,
        projected_l2=projected_l2,
        projected_reduction_fraction=reduction_fraction,
    )
