"""Diagnostic pure-swirl directions that remove the current frozen-collar blind spot.

This module is deliberately not a production candidate family.  It supplies two
small axisymmetric pure-swirl directions for a CR003/CR005 capacity probe.  The
modes vanish on an open core neighborhood, have compact outer support, and can
be orthogonalized against the existing outer-torque basis without changing
their value in the currently frozen inner annulus.

No force is fitted or inferred from a residual here, and no residual reduction
is claimed merely from obtaining nonzero span in the blind region.
"""
from __future__ import annotations

from math import sqrt
from numbers import Integral

import numpy as np

from .cutoffs import standard_cutoff
from .constrained_momentum_budget import angular_moment

INNER_ZERO_RADIUS = 0.1
INNER_FULL_RADIUS = sqrt(0.02)
OLD_CORRECTION_ZERO_RADIUS = sqrt(0.125)
MODE_COUNT = 2

# Deterministic inner-annulus points, all outside the proposed guard and inside
# the exact zero region of the current parent.outer_basis corrections.
REACH_POINTS = np.array(
    [
        [0.16, 0.0, 0.55],
        [0.19, 0.0, 0.75],
        [0.22, 0.0, 0.90],
        [0.26, 0.0, 1.05],
        [0.30, 0.0, 0.70],
        [0.32, 0.0, 0.8857142857142857],
        [0.34, 0.0, 1.10],
    ],
    dtype=float,
)
FROZEN_COLLAR_POINT = np.array([[0.32, 0.0, 0.8857142857142857]], dtype=float)


def _validated_points(points) -> np.ndarray:
    value = np.asarray(points, dtype=float)
    if value.ndim < 1 or value.shape[-1] != 3 or not np.all(np.isfinite(value)):
        raise ValueError("points must be finite with trailing dimension 3")
    return value


def _cutoff_array(values) -> np.ndarray:
    array = np.asarray(values, dtype=float)
    if not np.all(np.isfinite(array)):
        raise ValueError("cutoff arguments must be finite")
    flat = np.fromiter(
        (standard_cutoff(float(entry)) for entry in array.ravel()),
        dtype=float,
        count=array.size,
    )
    return flat.reshape(array.shape)


def guarded_swirl_basis(points) -> np.ndarray:
    """Return two compact pure-swirl directions with shape ``(...,3,2)``.

    ``1-standard_cutoff(r^2/.02)`` is exactly zero for ``r<=.1`` and exactly
    one for ``r>=sqrt(.02)``.  Thus the new directions stay away from every
    current core probe while reaching the frozen annulus around ``r=.32``.
    The second direction is an axial/radial localizer centered on that observed
    collar defect; the first is a broad guarded direction.
    """
    value = _validated_points(points)
    x, y, z = np.moveaxis(value, -1, 0)
    radius_sq = x * x + y * y
    axial_sq = z * z

    inner_guard = 1.0 - _cutoff_array(radius_sq / 0.02)
    outer = _cutoff_array(radius_sq / 4.0) * _cutoff_array(axial_sq / 4.0)
    common = inner_guard * outer
    localizer = np.exp(-((radius_sq - 0.105) / 0.055) ** 2 - ((axial_sq - 0.80) / 0.32) ** 2)
    scalars = np.stack((common, common * localizer), axis=-1)

    swirl = np.stack((-y, x, np.zeros_like(x)), axis=-1)
    return swirl[..., :, None] * scalars[..., None, :]


def moment_ratios(parent, order: int = 48) -> tuple[float, float]:
    """Ratios that cancel each diagnostic mode's global z-angular moment.

    The same quadrature operator is used in numerator and denominator.  This is
    a numerical cancellation at the selected order, not an interval proof.
    """
    if isinstance(order, bool) or not isinstance(order, Integral) or int(order) < 4:
        raise ValueError("order must be an integer >= 4")
    if not hasattr(parent, "outer_basis"):
        raise TypeError("parent must provide outer_basis(points,time)")
    order = int(order)
    denominator = angular_moment(parent.outer_basis, 0.5, order)
    if not np.isfinite(denominator) or denominator == 0.0:
        raise ValueError("parent outer basis must have nonzero finite angular moment")

    ratios = []
    for mode in range(MODE_COUNT):
        def field(points, time, selected=mode):
            del time
            return guarded_swirl_basis(points)[..., :, selected]

        numerator = angular_moment(field, 0.5, order)
        if not np.isfinite(numerator):
            raise ValueError("diagnostic mode angular moment is nonfinite")
        ratios.append(float(numerator / denominator))
    return tuple(ratios)


def moment_cancelled_basis(parent, points, time, *, ratios=None, order: int = 48) -> np.ndarray:
    """Subtract the existing outer basis so each new direction has zero moment.

    The subtraction does not erase the new inner-annulus reach because the
    existing ``parent.outer_basis`` is exactly zero throughout ``r^2<=.125``.
    """
    value = _validated_points(points)
    if ratios is None:
        ratios = moment_ratios(parent, order)
    ratios = np.asarray(ratios, dtype=float)
    if ratios.shape != (MODE_COUNT,) or not np.all(np.isfinite(ratios)):
        raise ValueError(f"ratios must be {MODE_COUNT} finite values")
    outer = np.asarray(parent.outer_basis(value, time), dtype=float)
    if outer.shape != value.shape or not np.all(np.isfinite(outer)):
        raise ValueError("parent outer_basis returned malformed values")
    return guarded_swirl_basis(value) - outer[..., :, None] * ratios


def reach_rank(points=REACH_POINTS, mode_count: int = MODE_COUNT) -> dict[str, float | int]:
    """Return rank/conditioning of the flattened spatial evaluation matrix."""
    if isinstance(mode_count, bool) or not isinstance(mode_count, Integral):
        raise TypeError("mode_count must be an integer")
    mode_count = int(mode_count)
    if not 1 <= mode_count <= MODE_COUNT:
        raise ValueError(f"mode_count must be in [1,{MODE_COUNT}]")
    basis = guarded_swirl_basis(points)[..., :mode_count]
    matrix = basis.reshape(-1, mode_count)
    norms = np.linalg.norm(matrix, axis=0)
    if np.any(norms == 0.0):
        return {"rank": int(np.linalg.matrix_rank(matrix)), "condition_number": float("inf")}
    normalized = matrix / norms
    return {
        "rank": int(np.linalg.matrix_rank(matrix)),
        "condition_number": float(np.linalg.cond(normalized)),
    }


def default_summary() -> dict:
    """Deterministic support-rank summary used by the checked-in artifact."""
    one = reach_rank(mode_count=1)
    two = reach_rank(mode_count=2)
    frozen = guarded_swirl_basis(FROZEN_COLLAR_POINT)[0]
    return {
        "inner_zero_radius": INNER_ZERO_RADIUS,
        "inner_full_radius": INNER_FULL_RADIUS,
        "old_correction_zero_radius": OLD_CORRECTION_ZERO_RADIUS,
        "reach_points": REACH_POINTS.tolist(),
        "rank_by_basis_dimension": {"1": one, "2": two},
        "frozen_point": FROZEN_COLLAR_POINT[0].tolist(),
        "frozen_point_y_components": frozen[1, :].tolist(),
    }


__all__ = [
    "FROZEN_COLLAR_POINT",
    "INNER_FULL_RADIUS",
    "INNER_ZERO_RADIUS",
    "MODE_COUNT",
    "OLD_CORRECTION_ZERO_RADIUS",
    "REACH_POINTS",
    "default_summary",
    "guarded_swirl_basis",
    "moment_cancelled_basis",
    "moment_ratios",
    "reach_rank",
]
