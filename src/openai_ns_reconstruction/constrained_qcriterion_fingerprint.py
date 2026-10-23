"""Truth-bounded Q-criterion morphology diagnostics for a callable velocity field.

This module is a visualization diagnostic only.  It samples a frozen
``velocity(points, time) -> [..., 3]`` callable on a Cartesian grid and
computes the VTK-style Q criterion from a second-order finite-difference
velocity gradient.  It does not alter the candidate, fit any parameter,
compute pressure/forcing, or establish Navier--Stokes validity.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Callable, Iterable

import numpy as np

VelocityCallable = Callable[[np.ndarray, float], np.ndarray]


@dataclass(frozen=True)
class QCriterionFingerprint:
    time: float
    grid_size: int
    half_width: float
    q_peak: float
    q_rms: float
    positive_q_volume_fraction: float
    positive_q_weight: float
    positive_q_r50: float | None
    positive_q_r90: float | None
    positive_q_abs_z50: float | None
    positive_q_abs_z90: float | None
    sampled_divergence_rms: float
    positive_core_present: bool
    claim_scope: str = "visualization_morphology_only"
    pde_validated: bool = False
    paper_exact: bool = False
    openai_field_identified: bool = False
    blowup_proved: bool = False

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class QCriterionResolutionAudit:
    fingerprints: tuple[QCriterionFingerprint, ...]
    relative_changes_to_finest: tuple[dict[str, float | None], ...]
    claim_scope: str = "visualization_resolution_only"
    pde_validated: bool = False
    paper_exact: bool = False
    openai_field_identified: bool = False

    def to_dict(self) -> dict[str, object]:
        return {
            "fingerprints": [item.to_dict() for item in self.fingerprints],
            "relative_changes_to_finest": list(self.relative_changes_to_finest),
            "claim_scope": self.claim_scope,
            "pde_validated": self.pde_validated,
            "paper_exact": self.paper_exact,
            "openai_field_identified": self.openai_field_identified,
        }


def _validate_geometry(time: float, half_width: float, grid_size: int) -> None:
    if not np.isfinite(time):
        raise ValueError("time must be finite")
    if not np.isfinite(half_width) or half_width <= 0.0:
        raise ValueError("half_width must be finite and positive")
    if not isinstance(grid_size, (int, np.integer)) or int(grid_size) < 5:
        raise ValueError("grid_size must be an integer >= 5")
    if int(grid_size) % 2 == 0:
        raise ValueError("grid_size must be odd so the symmetry axis is sampled")


def _weighted_quantile(values: np.ndarray, weights: np.ndarray, quantile: float) -> float:
    order = np.argsort(values, kind="mergesort")
    values = values[order]
    weights = weights[order]
    cumulative = np.cumsum(weights)
    cutoff = quantile * cumulative[-1]
    index = int(np.searchsorted(cumulative, cutoff, side="left"))
    return float(values[min(index, values.size - 1)])


def _relative_change(value: float, reference: float) -> float:
    scale = max(abs(reference), np.finfo(float).tiny)
    return float(abs(value - reference) / scale)


def diagnose_qcriterion_fingerprint(
    velocity: VelocityCallable,
    *,
    time: float,
    half_width: float = 1.5,
    grid_size: int = 33,
) -> QCriterionFingerprint:
    """Compute a finite-grid Q-criterion morphology fingerprint.

    The Q value matches the incompressible form used by VTK's
    ``vtkGradientFilter``:

    ``Q = -0.5 * trace((grad u)^2)``.

    Positive Q marks regions where local rotation dominates strain under the
    incompressible interpretation.  Because the gradient is sampled
    numerically, the returned divergence RMS is reported alongside Q and is
    not converted into a PDE pass/fail claim.
    """

    _validate_geometry(time, half_width, grid_size)
    n = int(grid_size)
    axis = np.linspace(-half_width, half_width, n, dtype=float)
    dx = float(axis[1] - axis[0])
    x, y, z = np.meshgrid(axis, axis, axis, indexing="ij")
    points = np.stack((x, y, z), axis=-1)

    values = np.asarray(velocity(points, float(time)), dtype=float)
    if values.shape != points.shape:
        raise ValueError("velocity must return an array with shape (..., 3)")
    if not np.all(np.isfinite(values)):
        raise ValueError("velocity returned non-finite values")

    u = values[..., 0]
    v = values[..., 1]
    w = values[..., 2]

    # np.gradient uses centered second-order differences in the interior.
    du_dx, du_dy, du_dz = np.gradient(u, dx, dx, dx, edge_order=2)
    dv_dx, dv_dy, dv_dz = np.gradient(v, dx, dx, dx, edge_order=2)
    dw_dx, dw_dy, dw_dz = np.gradient(w, dx, dx, dx, edge_order=2)

    # Restrict morphology statistics to the centered-stencil interior so
    # one-sided boundary differences do not dominate the core diagnostic.
    sl = np.s_[1:-1, 1:-1, 1:-1]
    du_dx = du_dx[sl]
    du_dy = du_dy[sl]
    du_dz = du_dz[sl]
    dv_dx = dv_dx[sl]
    dv_dy = dv_dy[sl]
    dv_dz = dv_dz[sl]
    dw_dx = dw_dx[sl]
    dw_dy = dw_dy[sl]
    dw_dz = dw_dz[sl]

    q = -0.5 * (du_dx**2 + dv_dy**2 + dw_dz**2) - (
        du_dy * dv_dx + du_dz * dw_dx + dv_dz * dw_dy
    )
    divergence = du_dx + dv_dy + dw_dz

    xi = x[sl]
    yi = y[sl]
    zi = z[sl]
    radius = np.sqrt(xi**2 + yi**2)

    if not np.all(np.isfinite(q)) or not np.all(np.isfinite(divergence)):
        raise ValueError("non-finite derivative diagnostic")

    q_positive = np.maximum(q, 0.0)
    mask = q > 0.0
    positive_weight = float(np.sum(q_positive) * dx**3)
    positive_fraction = float(np.mean(mask))
    core_present = bool(np.any(mask) and positive_weight > 0.0)

    r50: float | None = None
    r90: float | None = None
    z50: float | None = None
    z90: float | None = None
    if core_present:
        weights = q_positive[mask].ravel()
        rvals = radius[mask].ravel()
        zvals = np.abs(zi[mask]).ravel()
        r50 = _weighted_quantile(rvals, weights, 0.50)
        r90 = _weighted_quantile(rvals, weights, 0.90)
        z50 = _weighted_quantile(zvals, weights, 0.50)
        z90 = _weighted_quantile(zvals, weights, 0.90)

    return QCriterionFingerprint(
        time=float(time),
        grid_size=n,
        half_width=float(half_width),
        q_peak=float(np.max(q)),
        q_rms=float(np.sqrt(np.mean(q**2))),
        positive_q_volume_fraction=positive_fraction,
        positive_q_weight=positive_weight,
        positive_q_r50=r50,
        positive_q_r90=r90,
        positive_q_abs_z50=z50,
        positive_q_abs_z90=z90,
        sampled_divergence_rms=float(np.sqrt(np.mean(divergence**2))),
        positive_core_present=core_present,
    )


def audit_qcriterion_resolution(
    velocity: VelocityCallable,
    *,
    time: float,
    half_width: float,
    grid_sizes: Iterable[int],
) -> QCriterionResolutionAudit:
    """Compare Q morphology across three or more increasing odd grid sizes."""

    sizes = tuple(int(size) for size in grid_sizes)
    if len(sizes) < 3:
        raise ValueError("at least three grid sizes are required")
    if tuple(sorted(set(sizes))) != sizes:
        raise ValueError("grid_sizes must be strictly increasing and unique")

    fingerprints = tuple(
        diagnose_qcriterion_fingerprint(
            velocity,
            time=time,
            half_width=half_width,
            grid_size=size,
        )
        for size in sizes
    )
    finest = fingerprints[-1]
    rows: list[dict[str, float | None]] = []

    for item in fingerprints:
        row: dict[str, float | None] = {
            "grid_size": float(item.grid_size),
            "q_peak": _relative_change(item.q_peak, finest.q_peak),
            "q_rms": _relative_change(item.q_rms, finest.q_rms),
            "positive_q_volume_fraction": _relative_change(
                item.positive_q_volume_fraction,
                finest.positive_q_volume_fraction,
            )
            if finest.positive_q_volume_fraction > 0.0
            else None,
            "positive_q_weight": _relative_change(item.positive_q_weight, finest.positive_q_weight)
            if finest.positive_q_weight > 0.0
            else None,
            "sampled_divergence_rms": _relative_change(
                item.sampled_divergence_rms,
                finest.sampled_divergence_rms,
            )
            if finest.sampled_divergence_rms > 0.0
            else 0.0,
        }
        for name in (
            "positive_q_r50",
            "positive_q_r90",
            "positive_q_abs_z50",
            "positive_q_abs_z90",
        ):
            value = getattr(item, name)
            reference = getattr(finest, name)
            row[name] = (
                _relative_change(float(value), float(reference))
                if value is not None and reference is not None and reference != 0.0
                else (0.0 if value == reference else None)
            )
        rows.append(row)

    return QCriterionResolutionAudit(
        fingerprints=fingerprints,
        relative_changes_to_finest=tuple(rows),
    )
