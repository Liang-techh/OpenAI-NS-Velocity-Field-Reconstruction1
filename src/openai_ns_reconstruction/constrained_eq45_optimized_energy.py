"""Finite-window kinetic-energy audit for the bounded Eq45 profile update.

This module compares the frozen canonical ``Eq45VelocityCandidate`` against the
bounded two-profile update recorded by CR005.  Every velocity value is obtained
through the public ``candidate.at_points(...)->[u,v,w]`` interface.  The audit
uses deterministic midpoint Cartesian quadrature and decomposes kinetic energy
into cylindrical radial, swirl, and axial channels.

The reported energy is a finite-window visualization/research diagnostic.  The
physical-support connection is still unresolved, so these numbers must not be
used to claim the preregistered CR001 physical-energy normalization, PDE
validity, visual correspondence, paper exactness, or recovery of an OpenAI
hidden field.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Iterable

import numpy as np

from .constrained_eq45_candidate import Eq45VelocityCandidate
from .constrained_eq45_profile_force_optimization import candidate_with_profile_pair


DEFAULT_TIMES = (0.25, 0.50, 0.75)
DEFAULT_HALF_WIDTH = 2.0
DEFAULT_RESOLUTIONS = (24, 32, 48)
# Same cell width, dx=0.125, so this ladder isolates finite-window truncation.
DEFAULT_DOMAIN_LEVELS = ((2.0, 32), (2.5, 40), (3.0, 48))
EXPECTED_CANONICAL_SHA256 = (
    "48f1845fcfd71ec95495c98bb7aac3fca4653e748a8856a5421234e9601525a7"
)


def _validated_times(times: Iterable[float]) -> tuple[float, ...]:
    values = tuple(float(value) for value in times)
    if not values or not np.all(np.isfinite(values)):
        raise ValueError("times must be a nonempty finite sequence")
    if any(next_value <= value for value, next_value in zip(values, values[1:])):
        raise ValueError("times must be strictly increasing")
    return values


def _validated_resolutions(resolutions: Iterable[int]) -> tuple[int, ...]:
    values = tuple(resolutions)
    if len(values) < 3:
        raise ValueError("at least three resolution levels are required")
    if any(not isinstance(value, (int, np.integer)) or int(value) < 4 for value in values):
        raise ValueError("resolutions must be integers >= 4")
    result = tuple(int(value) for value in values)
    if any(next_value <= value for value, next_value in zip(result, result[1:])):
        raise ValueError("resolutions must be strictly increasing")
    return result


def _validated_domain_levels(
    levels: Iterable[tuple[float, int]],
) -> tuple[tuple[float, int], ...]:
    values = tuple(levels)
    if len(values) < 3:
        raise ValueError("at least three domain levels are required")
    checked: list[tuple[float, int]] = []
    cell_widths: list[float] = []
    for half_width, resolution in values:
        half_width = float(half_width)
        if not np.isfinite(half_width) or half_width <= 0.0:
            raise ValueError("domain half-widths must be positive and finite")
        if not isinstance(resolution, (int, np.integer)) or int(resolution) < 4:
            raise ValueError("domain resolutions must be integers >= 4")
        resolution = int(resolution)
        checked.append((half_width, resolution))
        cell_widths.append(2.0 * half_width / resolution)
    if any(b[0] <= a[0] for a, b in zip(checked, checked[1:])):
        raise ValueError("domain half-widths must be strictly increasing")
    reference = cell_widths[0]
    if any(not np.isclose(width, reference, rtol=0.0, atol=1e-14) for width in cell_widths[1:]):
        raise ValueError("domain ladder must keep Cartesian cell width fixed")
    return tuple(checked)


def _midpoint_axis(half_width: float, resolution: int) -> np.ndarray:
    dx = 2.0 * half_width / resolution
    return -half_width + (np.arange(resolution, dtype=float) + 0.5) * dx


def finite_window_component_energy(
    candidate: Eq45VelocityCandidate,
    time: float,
    *,
    half_width: float = DEFAULT_HALF_WIDTH,
    resolution: int = 32,
) -> dict[str, float]:
    """Integrate Cartesian kinetic energy and cylindrical channel contributions."""
    if not isinstance(candidate, Eq45VelocityCandidate):
        raise TypeError("candidate must be Eq45VelocityCandidate")
    time = float(time)
    half_width = float(half_width)
    if not np.isfinite(time):
        raise ValueError("time must be finite")
    if not np.isfinite(half_width) or half_width <= 0.0:
        raise ValueError("half_width must be positive and finite")
    if not isinstance(resolution, (int, np.integer)) or int(resolution) < 4:
        raise ValueError("resolution must be an integer >=4")
    resolution = int(resolution)

    axis = _midpoint_axis(half_width, resolution)
    dx = 2.0 * half_width / resolution
    radial_sum = 0.0
    swirl_sum = 0.0
    axial_sum = 0.0

    # Process one z-slab at a time.  This keeps peak memory bounded even when a
    # downstream caller raises the resolution for a visualization study.
    xx, yy = np.meshgrid(axis, axis, indexing="ij")
    radius = np.sqrt(xx * xx + yy * yy)
    nonaxis = radius > 0.0
    for z_value in axis:
        zz = np.full_like(xx, z_value)
        points = np.stack((xx, yy, zz), axis=-1).reshape(-1, 3)
        values = np.asarray(candidate.at_points(points, time), dtype=float).reshape(
            resolution, resolution, 3
        )
        if not np.all(np.isfinite(values)):
            raise RuntimeError("candidate returned nonfinite velocity during energy audit")
        u = values[..., 0]
        v = values[..., 1]
        w = values[..., 2]
        radial = np.zeros_like(u)
        swirl = np.zeros_like(u)
        radial[nonaxis] = (xx[nonaxis] * u[nonaxis] + yy[nonaxis] * v[nonaxis]) / radius[nonaxis]
        swirl[nonaxis] = (-yy[nonaxis] * u[nonaxis] + xx[nonaxis] * v[nonaxis]) / radius[nonaxis]
        # At the exact axis a regular axisymmetric field has zero transverse
        # velocity.  Midpoint grids used by the production audit are even and do
        # not hit r=0, but this convention keeps the helper robust for odd N.
        radial_sum += float(np.sum(radial * radial))
        swirl_sum += float(np.sum(swirl * swirl))
        axial_sum += float(np.sum(w * w))

    volume = dx**3
    radial_energy = 0.5 * radial_sum * volume
    swirl_energy = 0.5 * swirl_sum * volume
    axial_energy = 0.5 * axial_sum * volume
    total = radial_energy + swirl_energy + axial_energy
    if not np.isfinite(total) or total <= 0.0:
        raise RuntimeError("finite-window kinetic energy must remain positive and finite")
    return {
        "half_width": half_width,
        "resolution": resolution,
        "cell_width": dx,
        "radial_energy": radial_energy,
        "swirl_energy": swirl_energy,
        "axial_energy": axial_energy,
        "total_energy": total,
        "radial_fraction": radial_energy / total,
        "swirl_fraction": swirl_energy / total,
        "axial_fraction": axial_energy / total,
    }


def _relative_change(after: float, before: float) -> float:
    return float((after - before) / max(abs(before), np.finfo(float).tiny))


def _candidate_pair(
    candidate_path: str | Path | None,
    optimization_path: str | Path | None,
) -> tuple[Eq45VelocityCandidate, Eq45VelocityCandidate, dict]:
    root = Path(__file__).resolve().parents[2]
    candidate_path = Path(
        candidate_path or root / "artifacts/constrained/eq45_velocity_candidate_seed.json"
    )
    optimization_path = Path(
        optimization_path or root / "artifacts/constrained/eq45_profile_force_optimization.json"
    )
    baseline = Eq45VelocityCandidate.load_json(candidate_path)
    if baseline.sha256 != EXPECTED_CANONICAL_SHA256:
        raise RuntimeError("canonical Eq45 candidate identity drifted")
    optimization = json.loads(optimization_path.read_text(encoding="utf-8"))
    if optimization.get("task_id") != "CR005-EQ45-PROFILE-FORCE-OPT-015":
        raise RuntimeError("unexpected profile-optimization artifact")
    fit = optimization["fit"]
    optimized = candidate_with_profile_pair(
        baseline,
        float(fit["Phi_02"]),
        float(fit["F_02"]),
    )
    expected_sha = str(optimization["optimized_candidate"]["sha256"])
    if optimized.sha256 != expected_sha:
        raise RuntimeError("optimized candidate identity does not match checked artifact")
    roundtrip = Eq45VelocityCandidate.from_dict(optimized.to_dict())
    if roundtrip.sha256 != expected_sha:
        raise RuntimeError("optimized candidate failed governed serialization round-trip")
    return baseline, roundtrip, optimization


def audit_optimized_eq45_energy(
    *,
    times: Iterable[float] = DEFAULT_TIMES,
    half_width: float = DEFAULT_HALF_WIDTH,
    resolutions: Iterable[int] = DEFAULT_RESOLUTIONS,
    domain_levels: Iterable[tuple[float, int]] = DEFAULT_DOMAIN_LEVELS,
    candidate_path: str | Path | None = None,
    optimization_path: str | Path | None = None,
) -> dict:
    """Compare canonical and residual-optimized Eq45 finite-window energy balance."""
    checked_times = _validated_times(times)
    checked_resolutions = _validated_resolutions(resolutions)
    checked_domains = _validated_domain_levels(domain_levels)
    half_width = float(half_width)
    if not np.isfinite(half_width) or half_width <= 0.0:
        raise ValueError("half_width must be positive and finite")

    baseline, optimized, optimization = _candidate_pair(candidate_path, optimization_path)

    resolution_reports: list[dict] = []
    for time in checked_times:
        rows: list[dict] = []
        for resolution in checked_resolutions:
            canonical = finite_window_component_energy(
                baseline, time, half_width=half_width, resolution=resolution
            )
            changed = finite_window_component_energy(
                optimized, time, half_width=half_width, resolution=resolution
            )
            rows.append(
                {
                    "resolution": resolution,
                    "canonical": canonical,
                    "optimized": changed,
                    "total_energy_relative_change": _relative_change(
                        changed["total_energy"], canonical["total_energy"]
                    ),
                }
            )
        resolution_reports.append({"time": time, "levels": rows})

    domain_reports: list[dict] = []
    for time in checked_times:
        rows = []
        for domain_half_width, resolution in checked_domains:
            canonical = finite_window_component_energy(
                baseline,
                time,
                half_width=domain_half_width,
                resolution=resolution,
            )
            changed = finite_window_component_energy(
                optimized,
                time,
                half_width=domain_half_width,
                resolution=resolution,
            )
            rows.append(
                {
                    "half_width": domain_half_width,
                    "resolution": resolution,
                    "canonical": canonical,
                    "optimized": changed,
                }
            )
        domain_reports.append({"time": time, "levels": rows})

    finest_comparison: list[dict] = []
    for report in resolution_reports:
        finest = report["levels"][-1]
        canonical = finest["canonical"]
        changed = finest["optimized"]
        finest_comparison.append(
            {
                "time": report["time"],
                "resolution": finest["resolution"],
                "canonical_total_energy": canonical["total_energy"],
                "optimized_total_energy": changed["total_energy"],
                "total_energy_relative_change": _relative_change(
                    changed["total_energy"], canonical["total_energy"]
                ),
                "radial_energy_relative_change": _relative_change(
                    changed["radial_energy"], canonical["radial_energy"]
                ),
                "swirl_energy_relative_change": _relative_change(
                    changed["swirl_energy"], canonical["swirl_energy"]
                ),
                "axial_energy_relative_change": _relative_change(
                    changed["axial_energy"], canonical["axial_energy"]
                ),
                "canonical_fractions": {
                    "radial": canonical["radial_fraction"],
                    "swirl": canonical["swirl_fraction"],
                    "axial": canonical["axial_fraction"],
                },
                "optimized_fractions": {
                    "radial": changed["radial_fraction"],
                    "swirl": changed["swirl_fraction"],
                    "axial": changed["axial_fraction"],
                },
            }
        )

    resolution_stability = []
    for report in resolution_reports:
        previous = report["levels"][-2]
        finest = report["levels"][-1]
        resolution_stability.append(
            {
                "time": report["time"],
                "canonical_finest_relative_change": abs(
                    _relative_change(
                        finest["canonical"]["total_energy"],
                        previous["canonical"]["total_energy"],
                    )
                ),
                "optimized_finest_relative_change": abs(
                    _relative_change(
                        finest["optimized"]["total_energy"],
                        previous["optimized"]["total_energy"],
                    )
                ),
            }
        )

    domain_sensitivity = []
    for report in domain_reports:
        smallest = report["levels"][0]
        largest = report["levels"][-1]
        domain_sensitivity.append(
            {
                "time": report["time"],
                "canonical_L2_vs_largest_relative_deficit": _relative_change(
                    largest["canonical"]["total_energy"],
                    smallest["canonical"]["total_energy"],
                ),
                "optimized_L2_vs_largest_relative_deficit": _relative_change(
                    largest["optimized"]["total_energy"],
                    smallest["optimized"]["total_energy"],
                ),
            }
        )

    return {
        "schema": "eq45_optimized_energy_audit_v1",
        "task_id": "CR007-EQ45-OPTIMIZED-ENERGY-AUDIT-015",
        "canonical_candidate_sha256": baseline.sha256,
        "optimized_candidate_sha256": optimized.sha256,
        "source_optimization_task_id": optimization["task_id"],
        "source_profile_parameters": {
            "Phi_02": float(optimization["fit"]["Phi_02"]),
            "F_02": float(optimization["fit"]["F_02"]),
        },
        "quadrature": "cell-centered Cartesian midpoint rule",
        "energy_definition": "0.5 * integral_box |u|^2 dx",
        "component_definition": "cylindrical radial/swirl/axial kinetic-energy split",
        "times": list(checked_times),
        "resolution_half_width": half_width,
        "resolutions": list(checked_resolutions),
        "domain_levels": [list(level) for level in checked_domains],
        "resolution_reports": resolution_reports,
        "domain_reports": domain_reports,
        "finest_comparison": finest_comparison,
        "resolution_stability": resolution_stability,
        "domain_sensitivity": domain_sensitivity,
        "interpretation": (
            "finite-window energy/component-balance review of the bounded residual-"
            "optimized profile; this is visualization/research evidence only because "
            "the registered physical-support connection remains unresolved"
        ),
        "velocity_changed_by_source_optimization": True,
        "velocity_changed_by_this_audit": False,
        "physical_support_validated": False,
        "cr001_energy_normalization_status": "pending_unknown",
        "cr001_energy_gate_assessed": False,
        "visualization_ready": False,
        "visual_correspondence_verified": False,
        "pde_validated": False,
        "paper_exact": False,
        "openai_field_identified": False,
        "blowup_proved": False,
    }


def main() -> None:
    print(json.dumps(audit_optimized_eq45_energy(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
