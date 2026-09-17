"""Callable early-localized temporal Phi(1,0) supported Eq45 candidate.

Agent-7's target-free capacity screen established a minimal quadratic temporal
shape that keeps the successful early Phi(1,0) endpoint from the affine trial,
returns exactly to the static supported field at the midpoint and late endpoint,
and adds no spatial basis term.  This module materializes that already-screened
shape as an ordinary bounded, serializable velocity candidate.

The schedule is

    delta Phi10(tau) = 0.5 * early_delta * tau * (tau - 1),

with tau=-1,0,+1 at the delivery start/mid/end.  The default
``early_delta=-1.4`` therefore gives coefficient displacements ``-1.4, 0, 0``.
Every scalar-time snapshot delegates to the existing support-connected Eq45
candidate, so this module does not duplicate the Eq. (4.5) velocity formulas or
the physical-support transform.

This is a visualization-candidate representation only.  It does not fit a
public image, identify hidden OpenAI profiles, validate Navier--Stokes, establish
paper exactness, or prove blow-up.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from .constrained_eq45_supported_candidate import Eq45SupportedVelocityCandidate
from .constrained_eq45_supported_phi10_quadratic_temporal_capacity import (
    DEFAULT_EARLY_DELTA,
    quadratic_phi10_delta,
    quadratic_phi10_snapshot,
)


SCHEMA = "eq45_supported_phi10_early_localized_temporal_candidate_v1"
_MODE = (1, 0)
_SCHEDULE = "0.5*early_delta*tau*(tau-1)"
_TAU_DEFINITION = "delivery_interval_endpoints_map_to_minus1_plus1"

_CLASSIFICATION = {
    "velocity_interface": "user_requirement",
    "eq45_backbone": "public_source_fact",
    "base_supported_candidate": "parent_provenance",
    "phi10_mode": "autonomous_design",
    "quadratic_temporal_shape": "autonomous_design",
    "early_delta": "autonomous_design",
    "openai_hidden_profiles": "pending_unknown",
}

_TRUTH_BOUNDARY = {
    "callable_serializable": True,
    "velocity_export_ready": True,
    "visualization_candidate_only": True,
    "canonical_velocity_changed": False,
    "screened_temporal_shape_materialized": True,
    "production_temporal_shape_promoted": False,
    "physical_support_connection_implemented": True,
    "physical_support_validated": False,
    "spatial_basis_grown": False,
    "visualization_ready": False,
    "visual_correspondence_verified": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


@dataclass(frozen=True)
class Eq45SupportedPhi10EarlyLocalizedTemporalCandidate:
    """Support-connected Eq45 field with the screened quadratic Phi(1,0) shape."""

    base: Eq45SupportedVelocityCandidate
    early_delta: float = DEFAULT_EARLY_DELTA

    def __post_init__(self) -> None:
        if not isinstance(self.base, Eq45SupportedVelocityCandidate):
            raise TypeError("base must be Eq45SupportedVelocityCandidate")
        if _MODE not in self.profile_basis.mode_indices:
            raise ValueError("base profile basis must contain Phi(1,0)")
        if not np.isscalar(self.early_delta) or not np.isfinite(self.early_delta):
            raise ValueError("early_delta must be a finite scalar")
        if float(self.early_delta) == 0.0:
            raise ValueError("early_delta must be nonzero; keep the static base candidate otherwise")

        # A quadratic on tau in [-1,1] can have an interior extremum at tau=0.5.
        # Checking both endpoints, the midpoint, and that vertex therefore checks
        # the entire registered interval against the existing profile bound.
        check_tau = np.array([-1.0, 0.0, 0.5, 1.0], dtype=float)
        coefficients = np.asarray(self.coefficient_at_tau(check_tau), dtype=float)
        limit = float(self.profile_basis.coefficient_limit)
        if np.min(coefficients) < -limit or np.max(coefficients) > limit:
            raise ValueError(
                "quadratic Phi(1,0) schedule exceeds the existing profile coefficient bound "
                f"[-{limit}, {limit}] over the registered time interval"
            )

    @property
    def profile_basis(self):
        return self.base.parent.profile_basis

    @property
    def mode_index(self) -> int:
        return self.profile_basis.mode_indices.index(_MODE)

    @property
    def midpoint_coefficient(self) -> float:
        return float(self.profile_basis.phi_coefficients[self.mode_index])

    @property
    def time_start(self) -> float:
        return float(self.base.time_start)

    @property
    def time_end(self) -> float:
        return float(self.base.time_end)

    @property
    def time_midpoint(self) -> float:
        return 0.5 * (self.time_start + self.time_end)

    @property
    def base_sha256(self) -> str:
        return self.base.sha256

    def tau(self, time):
        """Return normalized delivery time in [-1,1]."""
        time_arr = np.asarray(time, dtype=float)
        if not np.all(np.isfinite(time_arr)):
            raise ValueError("time must be finite")
        if np.any((time_arr < self.time_start) | (time_arr > self.time_end)):
            raise ValueError(
                f"time must lie in declared interval [{self.time_start}, {self.time_end}]"
            )
        half_width = 0.5 * (self.time_end - self.time_start)
        return (time_arr - self.time_midpoint) / half_width

    def coefficient_at_tau(self, tau):
        """Return the bounded Phi(1,0) coefficient at normalized time ``tau``."""
        return self.midpoint_coefficient + quadratic_phi10_delta(
            tau, early_delta=float(self.early_delta)
        )

    def coefficient_at(self, time):
        """Return the bounded Phi(1,0) coefficient at physical delivery time."""
        return self.coefficient_at_tau(self.tau(time))

    def snapshot(self, time: float) -> Eq45SupportedVelocityCandidate:
        """Materialize one ordinary supported Eq45 snapshot at scalar ``time``."""
        time_arr = np.asarray(time, dtype=float)
        if time_arr.ndim != 0:
            raise ValueError("snapshot time must be a scalar")
        # Reuse the exact screened snapshot construction from PR #150 rather than
        # implementing a second profile/velocity path here.
        return quadratic_phi10_snapshot(
            self.base,
            float(time_arr),
            early_delta=float(self.early_delta),
        )

    def _validated_inputs(self, points, time) -> tuple[np.ndarray, np.ndarray]:
        points_arr = np.asarray(points, dtype=float)
        if points_arr.ndim == 0 or points_arr.shape[-1] != 3:
            raise ValueError("points must have shape (..., 3)")
        if not np.all(np.isfinite(points_arr)):
            raise ValueError("points must be finite")

        time_arr = np.asarray(time, dtype=float)
        if not np.all(np.isfinite(time_arr)):
            raise ValueError("time must be finite")
        try:
            time_arr = np.broadcast_to(time_arr, points_arr.shape[:-1])
        except ValueError as exc:
            raise ValueError("time must broadcast to points.shape[:-1]") from exc
        self.tau(time_arr)
        return points_arr, time_arr

    def velocity(self, points, time) -> np.ndarray:
        """Evaluate Cartesian ``[...,3] = [u,v,w]`` from the screened temporal field."""
        points_arr, time_arr = self._validated_inputs(points, time)
        flat_points = points_arr.reshape((-1, 3))
        flat_times = np.asarray(time_arr, dtype=float).reshape(-1)
        flat_out = np.empty_like(flat_points)

        for scalar_time in np.unique(flat_times):
            mask = flat_times == scalar_time
            flat_out[mask] = self.snapshot(float(scalar_time)).at_points(
                flat_points[mask], float(scalar_time)
            )

        out = flat_out.reshape(points_arr.shape)
        if not np.all(np.isfinite(out)):
            raise RuntimeError("early-localized temporal Eq45 velocity became nonfinite")
        return out

    def at_points(self, points, time) -> np.ndarray:
        return self.velocity(points, time)

    __call__ = velocity

    def velocity_xyz(self, x, y, z, time) -> np.ndarray:
        arrays = [np.asarray(value, dtype=float) for value in (x, y, z, time)]
        if not all(np.all(np.isfinite(value)) for value in arrays):
            raise ValueError("x, y, z and time must be finite")
        try:
            x_arr, y_arr, z_arr, time_arr = np.broadcast_arrays(*arrays)
        except ValueError as exc:
            raise ValueError("x, y, z and time must be broadcast-compatible") from exc
        points = np.stack((x_arr, y_arr, z_arr), axis=-1)
        return self.velocity(points, time_arr)

    def grid(self, x, y, z, times) -> np.ndarray:
        """Return regular samples in ``(time,x,y,z,component)`` order."""
        axes = [np.asarray(value, dtype=float) for value in (times, x, y, z)]
        if any(axis.ndim != 1 or axis.size == 0 for axis in axes):
            raise ValueError("grid axes and times must be nonempty 1D arrays")
        if not all(np.all(np.isfinite(axis)) for axis in axes):
            raise ValueError("grid axes and times must be finite")
        tt, xx, yy, zz = np.meshgrid(*axes, indexing="ij")
        return self.velocity_xyz(xx, yy, zz, tt)

    @property
    def sha256(self) -> str:
        payload = json.dumps(
            self.to_dict(), sort_keys=True, separators=(",", ":"), ensure_ascii=False
        ).encode("utf-8")
        return hashlib.sha256(payload).hexdigest()

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "classification": dict(_CLASSIFICATION),
            "base_sha256": self.base.sha256,
            "base_supported_candidate": self.base.to_dict(),
            "temporal_mode": {
                "family": "phi",
                "mode": list(_MODE),
                "early_delta": float(self.early_delta),
                "schedule": _SCHEDULE,
                "tau_definition": _TAU_DEFINITION,
            },
            "truth_boundary": dict(_TRUTH_BOUNDARY),
        }

    @classmethod
    def from_dict(
        cls, data: Mapping[str, Any]
    ) -> "Eq45SupportedPhi10EarlyLocalizedTemporalCandidate":
        if not isinstance(data, Mapping):
            raise ValueError("early-localized temporal candidate payload must be an object")
        if data.get("schema") != SCHEMA:
            raise ValueError(f"schema must be {SCHEMA!r}")
        if data.get("classification") != _CLASSIFICATION:
            raise ValueError("candidate source classification was modified")
        if data.get("truth_boundary") != _TRUTH_BOUNDARY:
            raise ValueError("candidate truth-boundary metadata was modified")
        if "base_supported_candidate" not in data or "temporal_mode" not in data:
            raise ValueError("candidate payload is missing a required component")

        base = Eq45SupportedVelocityCandidate.from_dict(data["base_supported_candidate"])
        if data.get("base_sha256") != base.sha256:
            raise ValueError("base_sha256 does not match the embedded supported candidate")

        temporal = data["temporal_mode"]
        if not isinstance(temporal, Mapping):
            raise ValueError("temporal_mode must be an object")
        if temporal.get("family") != "phi" or temporal.get("mode") != list(_MODE):
            raise ValueError("temporal mode must remain exactly Phi(1,0)")
        if temporal.get("schedule") != _SCHEDULE:
            raise ValueError("temporal schedule was modified")
        if temporal.get("tau_definition") != _TAU_DEFINITION:
            raise ValueError("temporal tau definition was modified")
        if "early_delta" not in temporal:
            raise ValueError("temporal_mode is missing early_delta")

        return cls(base=base, early_delta=float(temporal["early_delta"]))

    def save_json(self, path: str | Path) -> None:
        Path(path).write_text(
            json.dumps(self.to_dict(), indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )

    @classmethod
    def load_json(
        cls, path: str | Path
    ) -> "Eq45SupportedPhi10EarlyLocalizedTemporalCandidate":
        return cls.from_dict(json.loads(Path(path).read_text(encoding="utf-8")))
