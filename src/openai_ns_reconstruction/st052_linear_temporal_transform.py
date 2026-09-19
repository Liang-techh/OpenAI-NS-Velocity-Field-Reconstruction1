"""Portable linear-time activation for the ST052-M local-swirl transform.

This module migrates exactly the single temporal degree screened in Agent-7
PR #587 onto the package-level static transform kernel already integrated on
the constrained delivery ancestry.  It does not materialize the ST052-M parent
field and therefore is not a complete candidate capsule by itself.

The frozen law is

    g(t) = 2 * (t - 0.25)
    tau(t) = 0.05 * g(t)
    beta(t) = 0.08837490297155456 * g(t)

on t in [0.25, 0.75].  Thus the temporal child is exactly the redistributed
control at t=0.25 and exactly the static PR #559 child at t=0.75.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, replace
import hashlib
import json
from pathlib import Path

import numpy as np

from .st052_local_swirl_transform import (
    BaseVelocity,
    DEFAULT_SPEC as STATIC_SPEC,
    St052LocalSwirlTransformSpec,
    redistributed_control,
    transformed_velocity_points,
)

Array = np.ndarray

# Immutable identity of the already-integrated static PR #559 transform spec.
STATIC_TRANSFORM_SPEC_SHA256 = (
    "470103b68fd5ccece5d43d054711c894e40e4cb2f7ce890d52d04c4623bf7c25"
)


@dataclass(frozen=True)
class St052LinearTemporalTransformSpec:
    """Frozen identity for the exact Agent-7 PR #587 temporal activation."""

    schema_version: int = 1
    transform_id: str = "st052m_local_swirl_energy_linear_temporal_v1"
    parent_candidate_id: str = "ST052-M"
    parent_source_head: str = "b3b8bfdbe1077f9ec967d158602951997d81e17d"
    static_transform_spec_sha256: str = STATIC_TRANSFORM_SPEC_SHA256
    static_child_pr: int = 559
    static_child_head: str = "39b106ad8cb8df2064cabead3a12682089575e74"
    source_temporal_pr: int = 587
    source_temporal_head: str = "0b93819095f6c8576a7571bdc2d2fbef4154944d"
    activation_origin: float = 0.25
    activation_slope: float = 2.0
    time_interval: tuple[float, float] = (0.25, 0.75)
    static_taper_tau: float = 0.05
    static_shoulder_beta: float = 0.08837490297155456

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
    def load(cls, path: str | Path) -> "St052LinearTemporalTransformSpec":
        obj = json.loads(Path(path).read_text(encoding="utf-8"))
        payload = dict(obj["spec"])
        payload["time_interval"] = tuple(payload["time_interval"])
        spec = cls(**payload)
        if obj.get("spec_sha256") != spec.sha256():
            raise ValueError("ST052 temporal transform spec checksum mismatch")
        spec.validate_static_binding()
        return spec

    def validate_static_binding(self) -> None:
        if self.static_transform_spec_sha256 != STATIC_TRANSFORM_SPEC_SHA256:
            raise ValueError("temporal spec static-transform identity mismatch")
        if STATIC_SPEC.sha256() != STATIC_TRANSFORM_SPEC_SHA256:
            raise RuntimeError("integrated static ST052 transform spec drifted")
        if self.parent_candidate_id != STATIC_SPEC.parent_candidate_id:
            raise ValueError("temporal/static parent candidate mismatch")
        if self.parent_source_head != STATIC_SPEC.parent_source_head:
            raise ValueError("temporal/static parent source mismatch")
        if self.static_child_pr != STATIC_SPEC.source_child_pr:
            raise ValueError("temporal/static child PR mismatch")
        if self.static_child_head != STATIC_SPEC.source_child_head:
            raise ValueError("temporal/static child head mismatch")
        if self.static_taper_tau != STATIC_SPEC.taper_tau:
            raise ValueError("temporal/static taper endpoint mismatch")
        if self.static_shoulder_beta != STATIC_SPEC.shoulder_beta:
            raise ValueError("temporal/static shoulder endpoint mismatch")
        if self.time_interval != STATIC_SPEC.time_interval:
            raise ValueError("temporal/static time interval mismatch")


DEFAULT_TEMPORAL_SPEC = St052LinearTemporalTransformSpec()


def activation(
    time: float, spec: St052LinearTemporalTransformSpec = DEFAULT_TEMPORAL_SPEC
) -> float:
    """Return the exact frozen linear activation g(t)=2*(t-.25)."""
    time = float(time)
    lo, hi = spec.time_interval
    if not np.isfinite(time) or not (lo - 1.0e-12 <= time <= hi + 1.0e-12):
        raise ValueError("time outside frozen ST052 interval")
    return float(spec.activation_slope * (time - spec.activation_origin))


def effective_static_spec(
    time: float, spec: St052LinearTemporalTransformSpec = DEFAULT_TEMPORAL_SPEC
) -> St052LocalSwirlTransformSpec:
    """Return the static-kernel parameters representing the temporal child at t."""
    spec.validate_static_binding()
    gain = activation(time, spec)
    return replace(
        STATIC_SPEC,
        taper_tau=spec.static_taper_tau * gain,
        shoulder_beta=spec.static_shoulder_beta * gain,
    )


def temporal_transformed_velocity_points(
    base_velocity: BaseVelocity,
    points: Array,
    time: float,
    spec: St052LinearTemporalTransformSpec = DEFAULT_TEMPORAL_SPEC,
) -> Array:
    """Evaluate the exact PR #587 linear temporal transform on ``(n,3)`` points."""
    return transformed_velocity_points(
        base_velocity,
        points,
        float(time),
        effective_static_spec(float(time), spec),
    )


class St052LinearTemporalAdapter:
    """Broadcastable x/y/z/t adapter for the frozen temporal transform kernel.

    The exact ST052-M parent evaluator is still injected externally.  Therefore
    this adapter materializes the temporal transform identity only; it does not
    by itself establish a complete candidate ID/SHA, export readiness, visual
    correspondence, or PDE validity.
    """

    def __init__(
        self,
        base_velocity: BaseVelocity,
        *,
        parent_candidate_id: str,
        parent_source_head: str,
        spec: St052LinearTemporalTransformSpec = DEFAULT_TEMPORAL_SPEC,
    ) -> None:
        spec.validate_static_binding()
        if parent_candidate_id != spec.parent_candidate_id:
            raise ValueError("parent candidate identity mismatch")
        if parent_source_head != spec.parent_source_head:
            raise ValueError("parent source head mismatch")
        self._base_velocity = base_velocity
        self.spec = spec

    def points(self, points: Array, time: float) -> Array:
        return temporal_transformed_velocity_points(
            self._base_velocity, points, time, self.spec
        )

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
    "temporal_transform_kernel_materialized": True,
    "full_st052_parent_materialized_here": False,
    "complete_candidate_save_load_ready": False,
    "production_candidate_selected": False,
    "pressure_or_force_changed": False,
    "held_out_pde_residual_evaluated": False,
    "parent_or_static_pde_receipt_transferred": False,
    "public_image_numeric_target_used": False,
    "visualization_ready": False,
    "visual_correspondence_verified": False,
    "pde_validated": False,
    "source_correspondence_verified": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}
