"""Independent visualization-grid consistency audit for callable velocity fields.

This module never reads optimizer state or training loss. It compares values obtained
by trilinearly interpolating a sampled visualization grid with direct off-grid calls to
the same public ``velocity``/``VelocityField`` interface. The result is a rendering /
sampling diagnostic only; it is not PDE validation and cannot promote a candidate to
``pde_validated``.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Iterable, Protocol
import json

import numpy as np


class VelocityLike(Protocol):
    def at_points(self, points, time): ...
    def grid(self, x, y, z, times): ...


@dataclass(frozen=True)
class GridConsistencyRow:
    time: float
    grid_size: int
    samples: int
    max_vector_error: float
    rms_vector_error: float
    relative_rms_error: float

    def as_dict(self):
        return asdict(self)


def _validate_box(box):
    bounds = np.asarray(box, dtype=float)
    if bounds.shape != (3, 2) or not np.all(np.isfinite(bounds)):
        raise ValueError("box must be finite shape (3,2)")
    if np.any(bounds[:, 1] <= bounds[:, 0]):
        raise ValueError("each box upper bound must exceed its lower bound")
    return bounds


def _trilinear(values, points, bounds):
    """Interpolate a uniform Cartesian vector grid at interior points."""
    values = np.asarray(values, dtype=float)
    points = np.asarray(points, dtype=float)
    bounds = _validate_box(bounds)
    if values.ndim != 4 or values.shape[-1] != 3:
        raise ValueError("values must have shape (nx,ny,nz,3)")
    if len(set(values.shape[:3])) != 1 or values.shape[0] < 2:
        raise ValueError("audit requires a cubic uniform grid with size >=2")
    if points.ndim != 2 or points.shape[1] != 3 or not np.all(np.isfinite(points)):
        raise ValueError("points must be finite shape (n,3)")
    if np.any(points < bounds[:, 0]) or np.any(points > bounds[:, 1]):
        raise ValueError("off-grid points must lie inside the sampled box")

    n = values.shape[0]
    scale = (n - 1) / (bounds[:, 1] - bounds[:, 0])
    s = (points - bounds[:, 0]) * scale
    lo = np.floor(s).astype(int)
    lo = np.clip(lo, 0, n - 2)
    frac = np.clip(s - lo, 0.0, 1.0)

    out = np.zeros((len(points), 3), dtype=float)
    for dx in (0, 1):
        wx = (1.0 - frac[:, 0]) if dx == 0 else frac[:, 0]
        ix = lo[:, 0] + dx
        for dy in (0, 1):
            wy = (1.0 - frac[:, 1]) if dy == 0 else frac[:, 1]
            iy = lo[:, 1] + dy
            for dz in (0, 1):
                wz = (1.0 - frac[:, 2]) if dz == 0 else frac[:, 2]
                iz = lo[:, 2] + dz
                out += (wx * wy * wz)[:, None] * values[ix, iy, iz]
    return out


def audit_velocity_grid_consistency(
    field: VelocityLike,
    *,
    grid_sizes: Iterable[int] = (9, 17, 33),
    times: Iterable[float] = (0.35, 0.50, 0.65),
    box=((-1.0, 1.0), (-1.0, 1.0), (-1.0, 1.0)),
    sample_count: int = 512,
    seed: int = 350009,
):
    """Compare visualization-grid interpolation against direct off-grid velocity calls.

    The same held-out physical points are reused across grid levels at each time. This
    makes refinement ratios interpretable and prevents a changing sample set from
    masquerading as convergence. The audit deliberately uses only ``field.grid`` and
    ``field.at_points``; it never imports or inspects optimizer internals.
    """
    bounds = _validate_box(box)
    sizes = tuple(int(n) for n in grid_sizes)
    if len(sizes) < 3 or any(n < 2 for n in sizes) or any(b <= a for a, b in zip(sizes, sizes[1:])):
        raise ValueError("grid_sizes must contain at least three strictly increasing sizes >=2")
    times = tuple(float(t) for t in times)
    if not times or not np.all(np.isfinite(times)):
        raise ValueError("times must be a nonempty finite sequence")
    if not isinstance(sample_count, int) or sample_count < 8:
        raise ValueError("sample_count must be an integer >=8")

    rng = np.random.default_rng(seed)
    rows = []
    for time in times:
        points = rng.uniform(bounds[:, 0], bounds[:, 1], size=(sample_count, 3))
        direct = np.asarray(field.at_points(points, time), dtype=float)
        if direct.shape != (sample_count, 3) or not np.all(np.isfinite(direct)):
            raise FloatingPointError("direct velocity evaluator returned invalid values")
        direct_rms = float(np.sqrt(np.mean(np.sum(direct * direct, axis=1))))

        for n in sizes:
            axes = [np.linspace(lo, hi, n) for lo, hi in bounds]
            sampled = np.asarray(field.grid(*axes, [time]), dtype=float)
            if sampled.shape != (1, n, n, n, 3) or not np.all(np.isfinite(sampled)):
                raise FloatingPointError("grid evaluator returned invalid values or shape")
            interp = _trilinear(sampled[0], points, bounds)
            delta = interp - direct
            norms = np.linalg.norm(delta, axis=1)
            rms = float(np.sqrt(np.mean(norms * norms)))
            rows.append(GridConsistencyRow(
                time=time,
                grid_size=n,
                samples=sample_count,
                max_vector_error=float(np.max(norms)),
                rms_vector_error=rms,
                relative_rms_error=rms / max(direct_rms, np.finfo(float).tiny),
            ))
    return rows


def refinement_summary(rows):
    """Return per-time coarse/fine RMS ratios without declaring physical convergence."""
    grouped = {}
    for row in rows:
        grouped.setdefault(row.time, []).append(row)
    result = {}
    for time, group in grouped.items():
        ordered = sorted(group, key=lambda row: row.grid_size)
        if len(ordered) < 3:
            raise ValueError("at least three grid levels are required per time")
        result[time] = {
            "grid_sizes": [row.grid_size for row in ordered],
            "rms_errors": [row.rms_vector_error for row in ordered],
            "coarse_to_fine_ratio": ordered[-1].rms_vector_error / max(ordered[0].rms_vector_error, np.finfo(float).tiny),
            "finest_relative_rms": ordered[-1].relative_rms_error,
            "interpretation": "visualization-grid sampling only; not PDE evidence",
        }
    return result


def audit_default_velocity_field(**kwargs):
    """Load the packaged candidate through the user-facing API and audit its grid."""
    from .velocity_components import VelocityField
    field = VelocityField()
    rows = audit_velocity_grid_consistency(field, **kwargs)
    return {
        "candidate_sha256": field.sha256,
        "rows": [row.as_dict() for row in rows],
        "refinement": refinement_summary(rows),
        "claim_boundary": "visualization sampling audit only; not PDE validation",
    }


def main():
    report = audit_default_velocity_field(
        grid_sizes=(9, 17, 25), times=(0.35, 0.50, 0.65),
        sample_count=128, seed=350009,
    )
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
