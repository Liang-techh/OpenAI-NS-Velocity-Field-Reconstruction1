"""Generic finite-interval axial quadratic repair building block.

For supplied scalar ``F(R,Z)``, ``base_U(R,Z)``, ``Rmax`` and
``quadratic_target(Z)``, this module constructs

``U = base_U + sum_i c_i bump_i``

with four compact C-infinity radial bumps by default.  It first enforces
``integral U dR = 0`` and ``integral 2 R F U dR = 0`` linearly.  It then
solves
``integral (U^2 - R F^2) dR = quadratic_target``
inside the nullspace of those two constraints.  Among the quadratic roots it
selects the one with the smallest ``L2`` correction ``integral (U-base_U)^2``.

The implementation is a finite numerical repair primitive.  It does not
construct the Lei--Ren exterior/collar, certify a PDE, prove admissibility,
or claim completion of the five-moment problem.  Source context is Z. Lei
and X. Ren, arXiv:2609.35406v1, Section 2.5.
"""

from __future__ import annotations

import math
import operator
from dataclasses import dataclass, field
from typing import Any, Callable, Iterable, Sequence

import numpy as np
from numpy.polynomial.legendre import leggauss


SOURCE = "https://arxiv.org/html/2609.35406v1"
SOURCE_VERSION = "2609.35406v1"
SOURCE_SECTION = "Section 2.5"

ScalarProfile = Callable[[float, float], float]


def _finite(value: Any, name: str) -> float:
    """Convert one input to a finite real scalar."""

    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be finite") from exc
    if not math.isfinite(result):
        raise ValueError(f"{name} must be finite")
    return result


def _positive_order(value: Any, name: str = "quadrature_order") -> int:
    """Validate a positive integer quadrature order."""

    if isinstance(value, (bool, np.bool_)):
        raise TypeError(f"{name} must be a positive integer")
    try:
        result = operator.index(value)
    except TypeError as exc:
        raise TypeError(f"{name} must be a positive integer") from exc
    if result <= 0:
        raise ValueError(f"{name} must be a positive integer")
    return int(result)


def default_bump_intervals(Rmax: Any) -> tuple[tuple[float, float], ...]:
    """Return four disjoint supports away from both radial endpoints."""

    radius = _finite(Rmax, "Rmax")
    if radius <= 0.0:
        raise ValueError("Rmax must be positive")
    fractions = ((0.08, 0.20), (0.30, 0.42), (0.54, 0.66), (0.78, 0.90))
    return tuple((radius * left, radius * right) for left, right in fractions)


def flat_bump(R: Any, left: Any, right: Any) -> float:
    """Evaluate a unit-height compact C-infinity bump on ``(left,right)``."""

    radius = _finite(R, "R")
    a = _finite(left, "left")
    b = _finite(right, "right")
    if not a < b:
        raise ValueError("require left < right")
    if radius <= a or radius >= b:
        return 0.0
    s = (radius - a) / (b - a)
    # The +4 normalization makes the midpoint equal to one without changing
    # compact support or flatness at either endpoint.
    return float(math.exp(4.0 - 1.0 / s - 1.0 / (1.0 - s)))


def _validate_intervals(
    Rmax: float, bump_intervals: Sequence[tuple[float, float]] | None
) -> tuple[tuple[float, float], ...]:
    """Validate disjoint compact supports strictly inside ``(0,Rmax)``."""

    intervals = (
        default_bump_intervals(Rmax)
        if bump_intervals is None
        else tuple(tuple(pair) for pair in bump_intervals)
    )
    if len(intervals) < 3:
        raise ValueError("at least three radial bumps are required")
    normalized: list[tuple[float, float]] = []
    for index, pair in enumerate(intervals):
        if len(pair) != 2:
            raise ValueError(f"bump interval {index} must have two endpoints")
        left = _finite(pair[0], f"bump_intervals[{index}].left")
        right = _finite(pair[1], f"bump_intervals[{index}].right")
        if not 0.0 < left < right < Rmax:
            raise ValueError("all bump supports must lie strictly inside (0,Rmax)")
        normalized.append((left, right))
    normalized.sort()
    for previous, current in zip(normalized[:-1], normalized[1:]):
        if previous[1] >= current[0]:
            raise ValueError("bump supports must be disjoint")
    return tuple(normalized)


def _quadrature_grid(
    Rmax: float,
    intervals: Sequence[tuple[float, float]],
    order: int,
    breakpoints: Sequence[float] = (),
) -> tuple[np.ndarray, np.ndarray]:
    """Build a Gauss grid split at bump edges and supplied geometry breaks."""

    checked_breakpoints = []
    for index, breakpoint in enumerate(breakpoints):
        value = _finite(breakpoint, f"quadrature_breakpoints[{index}]")
        if not 0.0 < value < Rmax:
            raise ValueError(
                "quadrature breakpoints must lie strictly inside (0,Rmax)"
            )
        checked_breakpoints.append(value)
    edges = sorted(
        {
            0.0,
            Rmax,
            *checked_breakpoints,
            *(edge for pair in intervals for edge in pair),
        }
    )
    nodes, weights = leggauss(order)
    points: list[float] = []
    weights_all: list[float] = []
    for left, right in zip(edges[:-1], edges[1:]):
        width = right - left
        points.extend((left + width * (nodes + 1.0) / 2.0).tolist())
        weights_all.extend((width * weights / 2.0).tolist())
    return np.asarray(points, dtype=float), np.asarray(weights_all, dtype=float)


def _sample_profiles(
    F: ScalarProfile,
    base_U: ScalarProfile,
    Z: float,
    radii: np.ndarray,
    intervals: Sequence[tuple[float, float]],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Sample supplied profiles and all bump columns on one radial grid."""

    f_values = np.empty(radii.size, dtype=float)
    base_values = np.empty(radii.size, dtype=float)
    bumps = np.empty((radii.size, len(intervals)), dtype=float)
    for row, radius in enumerate(radii):
        try:
            f_value = float(F(float(radius), Z))
            base_value = float(base_U(float(radius), Z))
        except (TypeError, ValueError, IndexError) as exc:
            raise ValueError("F and base_U must return finite scalar values") from exc
        if not math.isfinite(f_value) or not math.isfinite(base_value):
            raise ValueError("F and base_U must return finite scalar values")
        f_values[row] = f_value
        base_values[row] = base_value
        for column, (left, right) in enumerate(intervals):
            bumps[row, column] = flat_bump(float(radius), left, right)
    return f_values, base_values, bumps


def _target_value(quadratic_target: Any, Z: float) -> float:
    """Evaluate a scalar or callable quadratic target."""

    value = quadratic_target(Z) if callable(quadratic_target) else quadratic_target
    return _finite(value, "quadratic_target(Z)")


@dataclass
class AxialQuadraticRepair:
    """Finite repaired axial profile and its linear/quadratic diagnostics."""

    status: str
    solvable: bool
    Z: float
    Rmax: float
    quadratic_target: float
    bump_intervals: tuple[tuple[float, float], ...]
    quadrature_breakpoints: tuple[float, ...]
    quadrature_order: int
    linear_rank: int
    linear_condition_number: float
    nullspace_dimension: int
    nullspace_condition_number: float
    particular_coefficients: np.ndarray | None
    nullspace_basis: np.ndarray | None
    coefficients: np.ndarray | None
    unconstrained_quadratic_minimum: float | None
    nullspace_level: float | None
    l2_change: float | None
    _F: ScalarProfile | None = field(default=None, repr=False, compare=False)
    _base_U: ScalarProfile | None = field(default=None, repr=False, compare=False)

    def correction(self, R: Any) -> float:
        """Evaluate the compact radial correction for the solved profile."""

        if not self.solvable or self.coefficients is None:
            raise ValueError("repair has no solvable coefficient vector")
        radius = _finite(R, "R")
        if not 0.0 <= radius <= self.Rmax:
            raise ValueError("R must lie in [0,Rmax]")
        return float(
            sum(
                coefficient * flat_bump(radius, left, right)
                for coefficient, (left, right) in zip(
                    self.coefficients, self.bump_intervals
                )
            )
        )

    def U(self, R: Any) -> float:
        """Evaluate ``base_U + correction`` at the repair's fixed ``Z``."""

        if self._base_U is None:
            raise ValueError("repair does not retain base_U")
        radius = _finite(R, "R")
        if not 0.0 <= radius <= self.Rmax:
            raise ValueError("R must lie in [0,Rmax]")
        return float(self._base_U(radius, self.Z) + self.correction(radius))

    def metadata(self) -> dict[str, Any]:
        """Return JSON-ready solver metadata without callable fields."""

        return {
            "status": self.status,
            "solvable": bool(self.solvable),
            "Z": self.Z,
            "Rmax": self.Rmax,
            "quadratic_target": self.quadratic_target,
            "bump_intervals": [list(pair) for pair in self.bump_intervals],
            "quadrature_breakpoints": list(self.quadrature_breakpoints),
            "quadrature_order": self.quadrature_order,
            "linear_rank": self.linear_rank,
            "linear_condition_number": self.linear_condition_number,
            "nullspace_dimension": self.nullspace_dimension,
            "nullspace_condition_number": self.nullspace_condition_number,
            "particular_coefficients": (
                None
                if self.particular_coefficients is None
                else self.particular_coefficients.tolist()
            ),
            "coefficients": (
                None if self.coefficients is None else self.coefficients.tolist()
            ),
            "unconstrained_quadratic_minimum": self.unconstrained_quadratic_minimum,
            "nullspace_level": self.nullspace_level,
            "l2_change": self.l2_change,
        }


def solve_axial_quadratic_repair(
    F: ScalarProfile,
    base_U: ScalarProfile,
    Rmax: Any,
    Z: Any,
    quadratic_target: Any,
    *,
    bump_intervals: Sequence[tuple[float, float]] | None = None,
    quadrature_breakpoints: Sequence[float] | None = None,
    quadrature_order: Any = 64,
) -> AxialQuadraticRepair:
    """Solve the finite two-linear-plus-one-quadratic axial subproblem.

    The two linear constraints are imposed first with an SVD minimum-norm
    particular coefficient vector.  The quadratic target is then solved on
    the linear nullspace.  The selected point minimizes the weighted radial
    ``L2`` norm of the correction among the roots represented by that
    nullspace.
    """

    if not callable(F) or not callable(base_U):
        raise TypeError("F and base_U must be callable as F(R,Z) and base_U(R,Z)")
    radius_max = _finite(Rmax, "Rmax")
    if radius_max <= 0.0:
        raise ValueError("Rmax must be positive")
    z = _finite(Z, "Z")
    order = _positive_order(quadrature_order)
    intervals = _validate_intervals(radius_max, bump_intervals)
    raw_breakpoints = (
        () if quadrature_breakpoints is None else tuple(quadrature_breakpoints)
    )
    checked_breakpoints = tuple(
        sorted(
            {
                _finite(value, f"quadrature_breakpoints[{index}]")
                for index, value in enumerate(raw_breakpoints)
            }
        )
    )
    for value in checked_breakpoints:
        if not 0.0 < value < radius_max:
            raise ValueError(
                "quadrature breakpoints must lie strictly inside (0,Rmax)"
            )
    target = _target_value(quadratic_target, z)
    radii, weights = _quadrature_grid(
        radius_max, intervals, order, checked_breakpoints
    )
    f_values, base_values, bumps = _sample_profiles(
        F, base_U, z, radii, intervals
    )

    weighted = weights[:, None]
    linear_matrix = np.vstack(
        (
            weights @ bumps,
            (weights * 2.0 * radii * f_values) @ bumps,
        )
    )
    linear_base = np.array(
        [weights @ base_values, weights @ (2.0 * radii * f_values * base_values)]
    )
    rhs = -linear_base
    singular_values = np.linalg.svd(linear_matrix, compute_uv=False)
    scale = singular_values[0] if singular_values.size else 0.0
    rank_tolerance = max(linear_matrix.shape) * np.finfo(float).eps * max(scale, 1.0)
    rank = int(np.sum(singular_values > rank_tolerance))
    linear_condition = (
        float(singular_values[0] / singular_values[rank - 1])
        if rank > 0
        else float("inf")
    )
    nullspace_dimension = len(intervals) - rank
    if rank < 2 or nullspace_dimension < 1:
        return AxialQuadraticRepair(
            status="linear_constraints_rank_deficient",
            solvable=False,
            Z=z,
            Rmax=radius_max,
            quadratic_target=target,
            bump_intervals=intervals,
            quadrature_breakpoints=checked_breakpoints,
            quadrature_order=order,
            linear_rank=rank,
            linear_condition_number=linear_condition,
            nullspace_dimension=nullspace_dimension,
            nullspace_condition_number=float("inf"),
            particular_coefficients=None,
            nullspace_basis=None,
            coefficients=None,
            unconstrained_quadratic_minimum=None,
            nullspace_level=None,
            l2_change=None,
            _F=F,
            _base_U=base_U,
        )

    # SVD minimum-norm particular solution and an orthonormal nullspace.
    U_svd, _, Vt_svd = np.linalg.svd(linear_matrix, full_matrices=True)
    particular = linear_matrix.T @ np.linalg.solve(
        linear_matrix @ linear_matrix.T, rhs
    )
    nullspace = Vt_svd[rank:, :].T
    gram = bumps.T @ (weights[:, None] * bumps)
    corrected_base = base_values + bumps @ particular
    cross = bumps.T @ (weights * corrected_base)
    energy_F = float(weights @ (radii * f_values * f_values))
    base_quadratic = float(weights @ (corrected_base * corrected_base)) - energy_F
    H = nullspace.T @ gram @ nullspace
    H_singular = np.linalg.svd(H, compute_uv=False)
    H_scale = H_singular[0] if H_singular.size else 0.0
    H_tol = max(H.shape) * np.finfo(float).eps * max(H_scale, 1.0)
    H_rank = int(np.sum(H_singular > H_tol))
    H_condition = (
        float(H_singular[0] / H_singular[H_rank - 1])
        if H_rank > 0
        else float("inf")
    )
    if H_rank < nullspace_dimension:
        return AxialQuadraticRepair(
            status="nullspace_quadratic_rank_deficient",
            solvable=False,
            Z=z,
            Rmax=radius_max,
            quadratic_target=target,
            bump_intervals=intervals,
            quadrature_breakpoints=checked_breakpoints,
            quadrature_order=order,
            linear_rank=rank,
            linear_condition_number=linear_condition,
            nullspace_dimension=nullspace_dimension,
            nullspace_condition_number=H_condition,
            particular_coefficients=particular,
            nullspace_basis=nullspace,
            coefficients=None,
            unconstrained_quadratic_minimum=None,
            nullspace_level=None,
            l2_change=None,
            _F=F,
            _base_U=base_U,
        )

    # q(y)=base_quadratic+2*g*y+y^T H y, with y in the linear nullspace.
    g = nullspace.T @ cross
    center_quadratic = -np.linalg.solve(H, g)
    q_minimum = base_quadratic - float(g @ np.linalg.solve(H, g))
    level = target - q_minimum
    level_tolerance = 2.0e-10 * max(1.0, abs(target), abs(q_minimum))
    if level < -level_tolerance:
        return AxialQuadraticRepair(
            status="quadratic_target_unreachable",
            solvable=False,
            Z=z,
            Rmax=radius_max,
            quadratic_target=target,
            bump_intervals=intervals,
            quadrature_breakpoints=checked_breakpoints,
            quadrature_order=order,
            linear_rank=rank,
            linear_condition_number=linear_condition,
            nullspace_dimension=nullspace_dimension,
            nullspace_condition_number=H_condition,
            particular_coefficients=particular,
            nullspace_basis=nullspace,
            coefficients=None,
            unconstrained_quadratic_minimum=q_minimum,
            nullspace_level=level,
            l2_change=None,
            _F=F,
            _base_U=base_U,
        )
    level = max(0.0, level)

    # L2(correction) has the same H metric in nullspace coordinates.  Choose
    # the point on the H sphere nearest its own objective center.
    objective_cross = nullspace.T @ gram @ particular
    center_objective = -np.linalg.solve(H, objective_cross)
    direction = center_objective - center_quadratic
    direction_norm = math.sqrt(max(0.0, float(direction @ H @ direction)))
    if level == 0.0:
        null_coordinate = center_quadratic
    elif direction_norm <= 1.0e-14:
        # The two quadratic centers coincide.  Every point on the level
        # ellipse has the same correction L2, so choose a deterministic unit
        # direction in the nullspace rather than returning its center, which
        # would miss the requested positive level.
        raw_direction = np.zeros(nullspace_dimension, dtype=float)
        raw_direction[0] = 1.0
        raw_norm = math.sqrt(max(0.0, float(raw_direction @ H @ raw_direction)))
        if raw_norm <= 0.0 or not math.isfinite(raw_norm):
            raise ArithmeticError("nullspace level direction is degenerate")
        null_coordinate = center_quadratic + math.sqrt(level) * raw_direction / raw_norm
    else:
        null_coordinate = center_quadratic + math.sqrt(level) * direction / direction_norm
    coefficients = particular + nullspace @ null_coordinate
    correction_values = bumps @ coefficients
    l2_change = float(weights @ (correction_values * correction_values))
    if not np.all(np.isfinite(coefficients)) or not math.isfinite(l2_change):
        raise ArithmeticError("axial quadratic repair produced non-finite coefficients")
    return AxialQuadraticRepair(
        status="solved",
        solvable=True,
        Z=z,
        Rmax=radius_max,
        quadratic_target=target,
        bump_intervals=intervals,
        quadrature_breakpoints=checked_breakpoints,
        quadrature_order=order,
        linear_rank=rank,
        linear_condition_number=linear_condition,
        nullspace_dimension=nullspace_dimension,
        nullspace_condition_number=H_condition,
        particular_coefficients=particular,
        nullspace_basis=nullspace,
        coefficients=coefficients,
        unconstrained_quadratic_minimum=q_minimum,
        nullspace_level=level,
        l2_change=l2_change,
        _F=F,
        _base_U=base_U,
    )


def evaluate_repair(
    repair: AxialQuadraticRepair,
    *,
    F: ScalarProfile | None = None,
    base_U: ScalarProfile | None = None,
    quadrature_order: Any | None = None,
) -> dict[str, Any]:
    """Reintegrate a solved repair independently at a requested order."""

    if not repair.solvable or repair.coefficients is None:
        raise ValueError("cannot evaluate an unsolved repair")
    profile_F = repair._F if F is None else F
    profile_base = repair._base_U if base_U is None else base_U
    if profile_F is None or profile_base is None:
        raise ValueError("F and base_U are required for evaluation")
    order = repair.quadrature_order if quadrature_order is None else _positive_order(
        quadrature_order
    )
    radii, weights = _quadrature_grid(
        repair.Rmax,
        repair.bump_intervals,
        order,
        repair.quadrature_breakpoints,
    )
    f_values, base_values, bumps = _sample_profiles(
        profile_F, profile_base, repair.Z, radii, repair.bump_intervals
    )
    correction = bumps @ repair.coefficients
    values = base_values + correction
    mass = float(weights @ values)
    mixed = float(weights @ (2.0 * radii * f_values * values))
    quadratic = float(weights @ (values * values - radii * f_values * f_values))
    l2_change = float(weights @ (correction * correction))
    return {
        "quadrature_order": order,
        "mass_moment": mass,
        "mixed_moment": mixed,
        "quadratic_value": quadratic,
        "quadratic_residual": quadratic - repair.quadratic_target,
        "l2_change": l2_change,
        "finite": bool(
            all(
                math.isfinite(value)
                for value in (mass, mixed, quadratic, l2_change)
            )
        ),
    }


__all__ = [
    "AxialQuadraticRepair",
    "SOURCE",
    "SOURCE_SECTION",
    "SOURCE_VERSION",
    "default_bump_intervals",
    "evaluate_repair",
    "flat_bump",
    "solve_axial_quadratic_repair",
]
