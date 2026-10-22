"""Portable transform kernel for the exact ST052-M local-swirl morphology child.

This module migrates only the deterministic post-processing transform used by
Agent-7 PR #559.  It deliberately does not bundle or identify the ST052-M parent
field, pressure, forcing, PDE receipt, or any OpenAI-private numerical data.

The kernel is useful as a delivery building block: once a checksum-bound
ST052-M base evaluator is supplied, :class:`St052LocalSwirlAdapter` exposes a
vectorized ``velocity(x, y, z, t)`` interface while preserving the exact frozen
redistribution, localized radial-Piola taper, and shoulder swirl compensation.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path
from typing import Callable

import numpy as np

Array = np.ndarray
BaseVelocity = Callable[[Array, float], Array]


@dataclass(frozen=True)
class St052LocalSwirlTransformSpec:
    """Frozen identity for the exact static PR #559 transform."""

    schema_version: int = 1
    transform_id: str = "st052m_local_swirl_energy_static_v1"
    parent_candidate_id: str = "ST052-M"
    parent_source_head: str = "b3b8bfdbe1077f9ec967d158602951997d81e17d"
    redistribution_source_head: str = "779ffca71066e2864496d37de55a7aafc45d6f57"
    source_child_pr: int = 559
    source_child_head: str = "39b106ad8cb8df2064cabead3a12682089575e74"
    redistribution_gain: float = 0.05
    redistribution_alpha: float = 2.520520814687742
    redistribution_scale: float = 1.0032534663681094
    inner_window: tuple[float, float] = (0.30, 1.05)
    outer_window: tuple[float, float] = (0.95, 1.85)
    taper_tau: float = 0.05
    taper_window_abs_z_over_2: tuple[float, float] = (0.50, 0.82)
    axial_support: float = 2.0
    shoulder_beta: float = 0.08837490297155456
    shoulder_window_abs_z_over_2: tuple[float, float] = (0.36, 0.49)
    time_interval: tuple[float, float] = (0.25, 0.75)
    post_transform_common_scale: float = 1.0

    def payload(self) -> dict:
        return asdict(self)

    def sha256(self) -> str:
        raw = json.dumps(
            self.payload(), sort_keys=True, separators=(",", ":"), allow_nan=False
        ).encode("utf-8")
        return hashlib.sha256(raw).hexdigest()

    def save(self, path: str | Path) -> None:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(
            json.dumps(
                {"spec": self.payload(), "spec_sha256": self.sha256()},
                indent=2,
                sort_keys=True,
                allow_nan=False,
            )
            + "\n",
            encoding="utf-8",
        )

    @classmethod
    def load(cls, path: str | Path) -> "St052LocalSwirlTransformSpec":
        obj = json.loads(Path(path).read_text(encoding="utf-8"))
        payload = dict(obj["spec"])
        for key in (
            "inner_window",
            "outer_window",
            "taper_window_abs_z_over_2",
            "shoulder_window_abs_z_over_2",
            "time_interval",
        ):
            payload[key] = tuple(payload[key])
        spec = cls(**payload)
        if obj.get("spec_sha256") != spec.sha256():
            raise ValueError("ST052 transform spec checksum mismatch")
        return spec


DEFAULT_SPEC = St052LocalSwirlTransformSpec()


def _compact_bump(value: Array, lo: float, hi: float) -> Array:
    value = np.asarray(value, dtype=float)
    mid = 0.5 * (lo + hi)
    half = 0.5 * (hi - lo)
    xi = (value - mid) / half
    out = np.zeros_like(value)
    mask = np.abs(xi) < 1.0
    if np.any(mask):
        xm = xi[mask]
        out[mask] = np.exp(1.0 - 1.0 / (1.0 - xm * xm))
    return out


def _symmetric_bump_and_derivative(
    z: Array, *, axial_support: float, window: tuple[float, float]
) -> tuple[Array, Array]:
    """C-infinity bump q(|z|/L) and dq/dz on one frozen axial window."""
    z = np.asarray(z, dtype=float)
    if not np.isfinite(z).all():
        raise ValueError("z must be finite")
    scaled = np.abs(z) / float(axial_support)
    lo, hi = window
    mid = 0.5 * (lo + hi)
    half = 0.5 * (hi - lo)
    xi = (scaled - mid) / half
    q = np.zeros_like(scaled)
    dq = np.zeros_like(scaled)
    mask = np.abs(xi) < 1.0
    if np.any(mask):
        xm = xi[mask]
        den = 1.0 - xm * xm
        qm = np.exp(1.0 - 1.0 / den)
        q[mask] = qm
        dq_dxi = qm * (-2.0 * xm) / (den * den)
        dq[mask] = dq_dxi * (1.0 / half) * np.sign(z[mask]) / float(axial_support)
    return q, dq


def _validate_points(points: Array) -> Array:
    points = np.asarray(points, dtype=float)
    if points.ndim != 2 or points.shape[1] != 3:
        raise ValueError("points must have shape (n,3)")
    if not np.isfinite(points).all():
        raise ValueError("points must be finite")
    return points


def _validate_time(time: float, spec: St052LocalSwirlTransformSpec) -> float:
    time = float(time)
    lo, hi = spec.time_interval
    if not np.isfinite(time) or not (lo - 1.0e-12 <= time <= hi + 1.0e-12):
        raise ValueError("time outside frozen ST052 interval")
    return time


def redistributed_control(
    base_velocity: BaseVelocity,
    points: Array,
    time: float,
    spec: St052LocalSwirlTransformSpec = DEFAULT_SPEC,
) -> Array:
    """Apply the frozen kappa=.05 radial swirl redistribution to ST052-M."""
    points = _validate_points(points)
    time = _validate_time(time, spec)
    u = np.asarray(base_velocity(points, time), dtype=float)
    if u.shape != points.shape or not np.isfinite(u).all():
        raise ValueError("base velocity must return finite shape-(n,3) values")
    u = u.copy()
    x, y = points[:, 0], points[:, 1]
    radius = np.hypot(x, y)
    mask = radius > 1.0e-14
    if np.any(mask):
        rx = x[mask] / radius[mask]
        ry = y[mask] / radius[mask]
        ux = u[mask, 0].copy()
        uy = u[mask, 1].copy()
        ur = rx * ux + ry * uy
        utheta = -ry * ux + rx * uy
        h = _compact_bump(radius[mask], *spec.inner_window) - spec.redistribution_alpha * _compact_bump(
            radius[mask], *spec.outer_window
        )
        utheta *= 1.0 + spec.redistribution_gain * h
        u[mask, 0] = rx * ur - ry * utheta
        u[mask, 1] = ry * ur + rx * utheta
    return spec.redistribution_scale * u


def _taper_map(
    points: Array, spec: St052LocalSwirlTransformSpec
) -> tuple[Array, Array, Array]:
    q, dq = _symmetric_bump_and_derivative(
        points[:, 2],
        axial_support=spec.axial_support,
        window=spec.taper_window_abs_z_over_2,
    )
    scale = 1.0 + spec.taper_tau * q
    derivative = spec.taper_tau * dq
    if np.any(scale <= 0.0):
        raise RuntimeError("localized taper lost orientation")
    mapped = points.copy()
    mapped[:, 0] *= scale
    mapped[:, 1] *= scale
    return mapped, scale, derivative


def _apply_shoulder_swirl(
    points: Array,
    velocity: Array,
    spec: St052LocalSwirlTransformSpec,
) -> Array:
    q, _ = _symmetric_bump_and_derivative(
        points[:, 2],
        axial_support=spec.axial_support,
        window=spec.shoulder_window_abs_z_over_2,
    )
    x, y = points[:, 0], points[:, 1]
    radius = np.hypot(x, y)
    out = velocity.copy()
    mask = radius > 1.0e-14
    if np.any(mask):
        xm, ym, rm = x[mask], y[mask], radius[mask]
        source = velocity[mask]
        utheta = (-ym * source[:, 0] + xm * source[:, 1]) / rm
        dtheta = -spec.shoulder_beta * q[mask] * utheta
        out[mask, 0] += dtheta * (-ym / rm)
        out[mask, 1] += dtheta * (xm / rm)
    return out


def transformed_velocity_points(
    base_velocity: BaseVelocity,
    points: Array,
    time: float,
    spec: St052LocalSwirlTransformSpec = DEFAULT_SPEC,
) -> Array:
    """Evaluate the exact frozen static #559 transform on ``(n,3)`` points."""
    points = _validate_points(points)
    time = _validate_time(time, spec)
    mapped, axial_scale, axial_scale_derivative = _taper_map(points, spec)
    u = redistributed_control(base_velocity, mapped, time, spec)
    x, y = points[:, 0], points[:, 1]
    uz = u[:, 2].copy()
    out = np.empty_like(u)
    out[:, 0] = axial_scale * u[:, 0] - axial_scale * axial_scale_derivative * x * uz
    out[:, 1] = axial_scale * u[:, 1] - axial_scale * axial_scale_derivative * y * uz
    out[:, 2] = axial_scale * axial_scale * uz
    out = _apply_shoulder_swirl(points, out, spec)
    return spec.post_transform_common_scale * out


class St052LocalSwirlAdapter:
    """Callable x/y/z/t adapter around a checksum-bound ST052-M base evaluator.

    This class materializes only the migrated transform kernel.  The caller must
    still provide the exact parent evaluator.  Consequently this class by itself
    is not a complete saved candidate and does not imply ``visualization_ready``
    or ``pde_validated``.
    """

    def __init__(
        self,
        base_velocity: BaseVelocity,
        *,
        parent_candidate_id: str,
        parent_source_head: str,
        spec: St052LocalSwirlTransformSpec = DEFAULT_SPEC,
    ) -> None:
        if parent_candidate_id != spec.parent_candidate_id:
            raise ValueError("parent candidate identity mismatch")
        if parent_source_head != spec.parent_source_head:
            raise ValueError("parent source head mismatch")
        self._base_velocity = base_velocity
        self.spec = spec

    def points(self, points: Array, time: float) -> Array:
        return transformed_velocity_points(self._base_velocity, points, time, self.spec)

    def velocity(self, x, y, z, t) -> Array:
        """Broadcast ``x,y,z,t`` and return velocity with final dimension 3."""
        xb, yb, zb, tb = np.broadcast_arrays(
            np.asarray(x, dtype=float),
            np.asarray(y, dtype=float),
            np.asarray(z, dtype=float),
            np.asarray(t, dtype=float),
        )
        if not all(np.isfinite(a).all() for a in (xb, yb, zb, tb)):
            raise ValueError("x, y, z, t must be finite")
        shape = xb.shape
        points = np.column_stack((xb.ravel(), yb.ravel(), zb.ravel()))
        times = tb.ravel()
        out = np.empty((len(points), 3), dtype=float)
        for time in np.unique(times):
            mask = times == time
            out[mask] = self.points(points[mask], float(time))
        return out.reshape(shape + (3,))


TRUTH_BOUNDARY = {
    "full_st052_parent_materialized_here": False,
    "complete_candidate_save_load_ready": False,
    "production_candidate_selected": False,
    "pressure_or_force_changed": False,
    "held_out_pde_residual_evaluated": False,
    "parent_pde_receipt_transferred": False,
    "public_image_numeric_target_used": False,
    "visualization_ready": False,
    "visual_correspondence_verified": False,
    "pde_validated": False,
    "source_correspondence_verified": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}
