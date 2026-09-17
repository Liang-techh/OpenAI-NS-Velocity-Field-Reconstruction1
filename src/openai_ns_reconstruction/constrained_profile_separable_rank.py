"""Low-rank expression-capacity audit for the selected paper-core profiles.

This diagnostic asks a narrow basis-growth question: how many separable angular
profile directions are needed to reproduce both a profile value and its radial
X derivative over a declared trusted inner domain? It is intentionally
independent of the optimizer, pressure, forcing, and PDE acceptance logic.

The stacked value/derivative SVD is useful for visualization-oriented basis
selection because matching profile values alone can hide derivative errors that
show up in streamline pitch and vorticity. A low-rank result is a compression /
expression-capacity statement only. It is not a Navier--Stokes validation or a
claim that the compressed modes are the hidden OpenAI profiles.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Iterable

import numpy as np


@dataclass(frozen=True)
class SeparableRankLevel:
    rank: int
    value_relative_error: float
    radial_derivative_relative_error: float
    worst_relative_error: float
    condition_number: float | None


@dataclass(frozen=True)
class SeparableRankReport:
    shape: tuple[int, int]
    value_frobenius_norm: float
    radial_derivative_frobenius_norm: float
    singular_values: tuple[float, ...]
    numerical_rank: int
    rank_rtol: float
    diagnostic_tolerance: float
    levels: tuple[SeparableRankLevel, ...]
    smallest_rank_below_diagnostic_tolerance: int | None
    claim_scope: str = "local_profile_expression_capacity_only"
    pde_validated: bool = False
    visual_correspondence_verified: bool = False
    paper_exact: bool = False
    openai_field_identified: bool = False

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass(frozen=True)
class SelectedPaperCoreRankAudit:
    inner_seed: dict
    X_range: tuple[float, float]
    eta_range: tuple[float, float]
    X_count: int
    eta_count: int
    requested_ranks: tuple[int, ...]
    rank_rtol: float
    diagnostic_tolerance: float
    channels: dict[str, SeparableRankReport]
    production_inner_seed_evaluated: bool = True
    velocity_changed: bool = False
    claim_scope: str = "selected_inner_profile_expression_capacity_only"
    visualization_ready: bool = False
    pde_validated: bool = False
    visual_correspondence_verified: bool = False
    paper_exact: bool = False
    openai_field_identified: bool = False

    def to_dict(self) -> dict:
        result = asdict(self)
        result["channels"] = {
            name: report.to_dict() for name, report in self.channels.items()
        }
        return result


def _matrix(name: str, value) -> np.ndarray:
    array = np.asarray(value, dtype=float)
    if array.ndim != 2 or min(array.shape) < 2:
        raise ValueError(
            f"{name} must be a finite 2D matrix with both dimensions >= 2"
        )
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must contain only finite values")
    return array


def _positive_finite(name: str, value: float) -> float:
    result = float(value)
    if not np.isfinite(result) or result <= 0.0:
        raise ValueError(f"{name} must be positive and finite")
    return result


def _ranks(values: Iterable[int], maximum: int) -> tuple[int, ...]:
    result = tuple(int(value) for value in values)
    if not result or any(value < 1 or value > maximum for value in result):
        raise ValueError(f"ranks must lie in [1, {maximum}]")
    if tuple(sorted(set(result))) != result:
        raise ValueError("ranks must be strictly increasing and unique")
    return result


def diagnose_separable_profile_rank(
    values,
    radial_derivatives,
    *,
    ranks: Iterable[int] = (1, 2, 3, 4, 6, 8),
    rank_rtol: float = 1e-3,
    diagnostic_tolerance: float = 1e-3,
) -> SeparableRankReport:
    """Audit a profile matrix and its X derivative with one shared SVD basis.

    ``values`` and ``radial_derivatives`` have shape ``(X, eta)``. Each block
    is normalized by its own Frobenius norm before stacking, so a large value
    scale cannot make the derivative block irrelevant. The same truncated SVD
    reconstructs both blocks, and errors are then reported in the original
    units as relative Frobenius errors.

    ``diagnostic_tolerance`` is only a basis-compression reporting threshold.
    It is not a PDE or visual-correspondence acceptance threshold.
    """

    values_array = _matrix("values", values)
    derivatives_array = _matrix("radial_derivatives", radial_derivatives)
    if values_array.shape != derivatives_array.shape:
        raise ValueError("values and radial_derivatives must have the same shape")

    value_norm = float(np.linalg.norm(values_array))
    derivative_norm = float(np.linalg.norm(derivatives_array))
    activity_floor = 100.0 * np.finfo(float).eps
    if value_norm <= activity_floor or derivative_norm <= activity_floor:
        raise ValueError("value and derivative blocks must both be numerically active")

    rank_rtol = _positive_finite("rank_rtol", rank_rtol)
    diagnostic_tolerance = _positive_finite(
        "diagnostic_tolerance", diagnostic_tolerance
    )

    combined = np.vstack(
        (values_array / value_norm, derivatives_array / derivative_norm)
    )
    requested_ranks = _ranks(ranks, min(combined.shape))
    left, singular_values, right_t = np.linalg.svd(combined, full_matrices=False)
    if singular_values.size == 0 or singular_values[0] <= 0.0:
        raise ValueError("combined profile matrix is numerically inactive")

    numerical_rank = int(
        np.count_nonzero(singular_values > rank_rtol * singular_values[0])
    )
    levels: list[SeparableRankLevel] = []
    X_count = values_array.shape[0]
    for rank in requested_ranks:
        reconstruction = (
            left[:, :rank] * singular_values[:rank]
        ) @ right_t[:rank, :]
        value_reconstruction = reconstruction[:X_count] * value_norm
        derivative_reconstruction = reconstruction[X_count:] * derivative_norm
        value_error = float(
            np.linalg.norm(values_array - value_reconstruction) / value_norm
        )
        derivative_error = float(
            np.linalg.norm(derivatives_array - derivative_reconstruction)
            / derivative_norm
        )
        denominator = float(singular_values[rank - 1])
        condition = (
            None
            if denominator <= np.finfo(float).tiny
            else float(singular_values[0] / denominator)
        )
        levels.append(
            SeparableRankLevel(
                rank=rank,
                value_relative_error=value_error,
                radial_derivative_relative_error=derivative_error,
                worst_relative_error=max(value_error, derivative_error),
                condition_number=condition,
            )
        )

    smallest = next(
        (
            level.rank
            for level in levels
            if level.worst_relative_error <= diagnostic_tolerance
        ),
        None,
    )
    return SeparableRankReport(
        shape=tuple(int(value) for value in values_array.shape),
        value_frobenius_norm=value_norm,
        radial_derivative_frobenius_norm=derivative_norm,
        singular_values=tuple(float(value) for value in singular_values),
        numerical_rank=numerical_rank,
        rank_rtol=rank_rtol,
        diagnostic_tolerance=diagnostic_tolerance,
        levels=tuple(levels),
        smallest_rank_below_diagnostic_tolerance=smallest,
    )


def _grid(name: str, values, *, lower: float, upper: float) -> np.ndarray:
    array = np.asarray(values, dtype=float)
    if array.ndim != 1 or array.size < 2 or not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must be a finite 1D grid with at least two points")
    if np.any(np.diff(array) <= 0.0):
        raise ValueError(f"{name} must be strictly increasing")
    if array[0] < lower or array[-1] > upper:
        raise ValueError(f"{name} must stay inside [{lower}, {upper}]")
    return array


def audit_selected_paper_core_separable_rank(
    *,
    X_values=None,
    eta_values=None,
    ranks: Iterable[int] = (1, 2, 3, 4, 6, 8),
    rank_rtol: float = 1e-3,
    diagnostic_tolerance: float = 1e-3,
) -> SelectedPaperCoreRankAudit:
    """Run the audit on the currently selected nonlinear inner seed.

    The selected seed matches the active function-first checkpoint:
    ``PaperCoreSeries(PaperCoreReference(sigma=.5), maxdegree=14,
    eta_nodes=257)``. The default sampled rectangle remains strictly inside
    the finite series domain ``Lambda*X <= 4.1`` and avoids eta endpoint
    dominance. These finite grids and compression tolerances are autonomous
    diagnostic choices, not public OpenAI facts.
    """

    from .paper_core_reference import PaperCoreReference
    from .paper_core_series import PaperCoreSeries

    if X_values is None:
        X_values = np.linspace(0.001, 0.409, 25)
    if eta_values is None:
        eta_values = np.linspace(-0.9, 0.9, 65)
    X_grid = _grid("X_values", X_values, lower=0.0, upper=0.409)
    eta_grid = _grid("eta_values", eta_values, lower=-0.9, upper=0.9)

    series = PaperCoreSeries(
        PaperCoreReference(sigma=0.5), maxdegree=14, eta_nodes=257
    )
    X = X_grid[:, None]
    eta = eta_grid[None, :]
    channels = {
        "phi": (
            series.phi_value(X, eta),
            series.phi_radial_derivative(X, eta),
        ),
        "F": (
            series.F(X, eta),
            series.F_radial_derivative(X, eta),
        ),
        "U": (
            series.U(X, eta),
            series.U_radial_derivative(X, eta),
        ),
    }
    reports = {
        name: diagnose_separable_profile_rank(
            values,
            derivatives,
            ranks=ranks,
            rank_rtol=rank_rtol,
            diagnostic_tolerance=diagnostic_tolerance,
        )
        for name, (values, derivatives) in channels.items()
    }
    requested_levels = next(iter(reports.values())).levels
    return SelectedPaperCoreRankAudit(
        inner_seed={
            "family": "PaperCoreSeries",
            "sigma": 0.5,
            "maxdegree": 14,
            "eta_nodes": 257,
            "h": series.reference.h,
            "j": series.reference.j,
            "Lambda": series.reference.Lambda,
            "C": series.reference.C,
            "status": "selected nonlinear inner seed; no exterior/full-NS acceptance",
        },
        X_range=(float(X_grid[0]), float(X_grid[-1])),
        eta_range=(float(eta_grid[0]), float(eta_grid[-1])),
        X_count=int(X_grid.size),
        eta_count=int(eta_grid.size),
        requested_ranks=tuple(level.rank for level in requested_levels),
        rank_rtol=float(rank_rtol),
        diagnostic_tolerance=float(diagnostic_tolerance),
        channels=reports,
    )


def main() -> None:
    import json

    print(json.dumps(audit_selected_paper_core_separable_rank().to_dict(), indent=2))


if __name__ == "__main__":
    main()
