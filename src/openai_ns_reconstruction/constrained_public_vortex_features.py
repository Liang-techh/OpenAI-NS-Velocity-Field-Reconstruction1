"""Public-observable vortex features for bounded basis-growth experiments.

The OpenAI public visualization describes a spinning vortex with trajectories
that spiral inward and undergo axial stretching.  This module converts only
those public qualitative observations into continuous, dimensionless probe
features that can be passed to local basis-sensitivity tools.

It deliberately defines no visual acceptance threshold and does not compare
against hidden parameters.  A favorable feature value is visualization-side
experiment evidence only; it is not Navier--Stokes validation or paper-exact
field identification.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Callable, Sequence

import numpy as np


OPENAI_PUBLIC_SOURCE = "https://openai.com/index/navier-stokes-solution/"
PUBLIC_OBSERVATIONS = (
    "spinning swirl",
    "inward spiraling trajectories",
    "axial stretching",
)


@dataclass(frozen=True)
class PublicVortexFeatures:
    time: float
    radii: tuple[float, ...]
    z_half_levels: tuple[float, ...]
    angular_samples: int
    active_fraction: float
    inward_radial_margin: float
    circulation_fraction: float
    inward_spiral_score: float
    axial_stretch_margin: float
    source_url: str = OPENAI_PUBLIC_SOURCE
    source_classification: str = "public_source_fact"
    claim_scope: str = "visualization_basis_target_only"
    pde_validated: bool = False
    paper_exact: bool = False
    openai_field_identified: bool = False

    def feature_vector(self) -> np.ndarray:
        return np.asarray(
            (
                self.inward_radial_margin,
                self.circulation_fraction,
                self.inward_spiral_score,
                self.axial_stretch_margin,
            ),
            dtype=float,
        )

    def to_dict(self) -> dict:
        payload = asdict(self)
        payload["public_observations"] = list(PUBLIC_OBSERVATIONS)
        return payload


@dataclass(frozen=True)
class PublicVortexResolutionAudit:
    reports: tuple[PublicVortexFeatures, ...]
    angular_samples: tuple[int, ...]
    absolute_deltas_to_finest: tuple[tuple[float, ...], ...]
    claim_scope: str = "visualization_resolution_sensitivity_only"
    pde_validated: bool = False
    paper_exact: bool = False


def _positive_vector(name: str, values: Sequence[float]) -> np.ndarray:
    arr = np.asarray(tuple(values), dtype=float)
    if arr.ndim != 1 or arr.size == 0:
        raise ValueError(f"{name} must be a nonempty 1D sequence")
    if not np.all(np.isfinite(arr)) or np.any(arr <= 0.0):
        raise ValueError(f"{name} must contain positive finite values")
    if np.any(np.diff(arr) <= 0.0):
        raise ValueError(f"{name} must be strictly increasing")
    return arr


def _probe_points(
    radii: np.ndarray, z_half_levels: np.ndarray, angular_samples: int
) -> tuple[np.ndarray, np.ndarray]:
    if not isinstance(angular_samples, int) or angular_samples < 8:
        raise ValueError("angular_samples must be an integer >= 8")
    angles = 2.0 * np.pi * np.arange(angular_samples, dtype=float) / angular_samples
    rr, zz, aa = np.meshgrid(
        radii,
        np.concatenate((-z_half_levels[::-1], z_half_levels)),
        angles,
        indexing="ij",
    )
    x = rr * np.cos(aa)
    y = rr * np.sin(aa)
    points = np.stack((x, y, zz), axis=-1).reshape(-1, 3)
    return points, zz.reshape(-1)


def diagnose_public_vortex_features(
    velocity: Callable[[np.ndarray, float], np.ndarray],
    time: float,
    *,
    radii: Sequence[float],
    z_half_levels: Sequence[float],
    angular_samples: int = 32,
    speed_floor: float = 1e-12,
) -> PublicVortexFeatures:
    """Measure public-observable vortex morphology on symmetric ring probes.

    Positive ``inward_radial_margin`` means radial motion is inward on average.
    Positive ``axial_stretch_margin`` means motion points away from the midplane
    on average. ``inward_spiral_score`` is positive only when inward radial
    motion coexists with circulation.  No numeric pass threshold is imposed.
    """
    if not np.isfinite(time):
        raise ValueError("time must be finite")
    if not np.isfinite(speed_floor) or speed_floor <= 0.0:
        raise ValueError("speed_floor must be positive and finite")

    r_values = _positive_vector("radii", radii)
    z_values = _positive_vector("z_half_levels", z_half_levels)
    points, z = _probe_points(r_values, z_values, angular_samples)

    values = np.asarray(velocity(points, float(time)), dtype=float)
    if values.shape != points.shape:
        raise ValueError("velocity must return an array with shape (N,3)")
    if not np.all(np.isfinite(values)):
        raise ValueError("velocity returned non-finite values")

    x = points[:, 0]
    y = points[:, 1]
    radius = np.hypot(x, y)
    ur = (x * values[:, 0] + y * values[:, 1]) / radius
    utheta = (-y * values[:, 0] + x * values[:, 1]) / radius
    w = values[:, 2]
    speed = np.linalg.norm(values, axis=1)
    active = speed > speed_floor
    if not np.any(active):
        raise ValueError("velocity is numerically inactive on all public-vortex probes")

    ur = ur[active]
    utheta = utheta[active]
    w = w[active]
    z = z[active]
    speed = speed[active]

    swirl_smooth = np.sqrt(utheta * utheta + speed_floor * speed_floor)
    transverse_sq = ur * ur + utheta * utheta + speed_floor * speed_floor

    inward_margin = float(np.mean(-ur / speed))
    circulation_fraction = float(np.mean(swirl_smooth / speed))
    spiral_score = float(np.mean((-ur * swirl_smooth) / transverse_sq))
    stretch_margin = float(np.mean(np.sign(z) * w / speed))

    features = np.asarray(
        (inward_margin, circulation_fraction, spiral_score, stretch_margin),
        dtype=float,
    )
    if not np.all(np.isfinite(features)):
        raise ValueError("public-vortex feature calculation produced non-finite values")

    return PublicVortexFeatures(
        time=float(time),
        radii=tuple(float(v) for v in r_values),
        z_half_levels=tuple(float(v) for v in z_values),
        angular_samples=angular_samples,
        active_fraction=float(np.mean(active)),
        inward_radial_margin=inward_margin,
        circulation_fraction=circulation_fraction,
        inward_spiral_score=spiral_score,
        axial_stretch_margin=stretch_margin,
    )


def audit_public_vortex_feature_resolution(
    velocity: Callable[[np.ndarray, float], np.ndarray],
    time: float,
    *,
    radii: Sequence[float],
    z_half_levels: Sequence[float],
    angular_samples: Sequence[int] = (16, 32, 64),
    speed_floor: float = 1e-12,
) -> PublicVortexResolutionAudit:
    levels = tuple(int(v) for v in angular_samples)
    if len(levels) < 3:
        raise ValueError("at least three angular resolutions are required")
    if any(v < 8 for v in levels) or any(a >= b for a, b in zip(levels, levels[1:])):
        raise ValueError("angular resolutions must be strictly increasing integers >= 8")

    reports = tuple(
        diagnose_public_vortex_features(
            velocity,
            time,
            radii=radii,
            z_half_levels=z_half_levels,
            angular_samples=level,
            speed_floor=speed_floor,
        )
        for level in levels
    )
    finest = reports[-1].feature_vector()
    deltas = tuple(
        tuple(float(v) for v in np.abs(report.feature_vector() - finest))
        for report in reports
    )
    return PublicVortexResolutionAudit(
        reports=reports,
        angular_samples=levels,
        absolute_deltas_to_finest=deltas,
    )
