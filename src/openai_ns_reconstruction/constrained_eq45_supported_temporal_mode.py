"""One bounded affine temporal coefficient on a support-connected Eq45 candidate.

This module adds exactly one time-dependent coefficient to an *existing* Phi/F
mode.  It deliberately does not grow the spatial basis.  At every requested
time it materializes an ordinary :class:`Eq45VelocityCandidate`, then evaluates
that snapshot through the already-governed :class:`Eq45SupportedVelocityCandidate`.
Thus the public Eq. (4.5) evaluator and physical-support transform remain the
single numerical implementation of ``[u,v,w]``.

The temporal slope and target mode are autonomous visualization-design choices.
They do not identify hidden OpenAI profile data and do not establish visual
correspondence, Navier--Stokes validity, paper exactness, or blow-up.
"""
from __future__ import annotations

from dataclasses import dataclass, replace
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from .constrained_eq45_profile_basis import Eq45CompactProfileBasis
from .constrained_eq45_supported_candidate import Eq45SupportedVelocityCandidate


SCHEMA = "eq45_supported_affine_temporal_mode_v1"

_CLASSIFICATION = {
    "velocity_interface": "user_requirement",
    "eq45_backbone": "public_source_fact",
    "base_supported_candidate": "parent_provenance",
    "temporal_mode_family": "autonomous_design",
    "temporal_mode_index": "autonomous_design",
    "temporal_slope": "autonomous_design",
    "openai_hidden_profiles": "pending_unknown",
}

_TRUTH_BOUNDARY = {
    "callable_serializable": True,
    "velocity_export_ready": True,
    "physical_support_connection_implemented": True,
    "physical_support_validated": False,
    "visualization_ready": False,
    "visual_correspondence_verified": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
    "spatial_basis_grown": False,
    "affine_temporal_mode_active": True,
}


@dataclass(frozen=True)
class Eq45SupportedAffineTemporalModeCandidate:
    """Support-connected Eq45 field with one affine-in-time existing mode.

    The selected coefficient is

        c(t) = c_mid + slope * tau(t),

    where ``tau=-1`` at ``time_start``, ``0`` at the interval midpoint, and
    ``+1`` at ``time_end``.  The existing profile coefficient limit is enforced
    at both endpoints, which controls the whole affine interval.
    """

    base: Eq45SupportedVelocityCandidate
    family: str
    mode_i: int
    mode_j: int
    slope: float

    def __post_init__(self) -> None:
        if not isinstance(self.base, Eq45SupportedVelocityCandidate):
            raise TypeError("base must be Eq45SupportedVelocityCandidate")
        if self.family not in {"phi", "swirl"}:
            raise ValueError("family must be 'phi' or 'swirl'")
        if not isinstance(self.mode_i, int) or not isinstance(self.mode_j, int):
            raise ValueError("mode indices must be integers")
        if self.mode_i < 0 or self.mode_j < 0:
            raise ValueError("mode indices must be nonnegative")
        if (self.mode_i, self.mode_j) not in self.profile_basis.mode_indices:
            raise ValueError("target mode must already exist in the base profile basis")
        if not np.isscalar(self.slope) or not np.isfinite(self.slope):
            raise ValueError("slope must be a finite scalar")
        if float(self.slope) == 0.0:
            raise ValueError("slope must be nonzero; keep the base candidate for a static field")

        lo = self.midpoint_coefficient - abs(float(self.slope))
        hi = self.midpoint_coefficient + abs(float(self.slope))
        limit = float(self.profile_basis.coefficient_limit)
        if lo < -limit or hi > limit:
            raise ValueError(
                "affine coefficient exceeds the existing profile coefficient bound "
                f"[-{limit}, {limit}] over the registered time interval"
            )

    @property
    def profile_basis(self) -> Eq45CompactProfileBasis:
        return self.base.parent.profile_basis

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
    def mode(self) -> tuple[int, int]:
        return (self.mode_i, self.mode_j)

    @property
    def mode_index(self) -> int:
        return self.profile_basis.mode_indices.index(self.mode)

    @property
    def midpoint_coefficient(self) -> float:
        coefficients = (
            self.profile_basis.phi_coefficients
            if self.family == "phi"
            else self.profile_basis.swirl_coefficients
        )
        return float(coefficients[self.mode_index])

    @property
    def base_sha256(self) -> str:
        return self.base.sha256

    def tau(self, time):
        """Return normalized affine time in [-1, 1] on the registered interval."""
        time_arr = np.asarray(time, dtype=float)
        if not np.all(np.isfinite(time_arr)):
            raise ValueError("time must be finite")
        if np.any((time_arr < self.time_start) | (time_arr > self.time_end)):
            raise ValueError(
                f"time must lie in declared interval [{self.time_start}, {self.time_end}]"
            )
        half_width = 0.5 * (self.time_end - self.time_start)
        return (time_arr - self.time_midpoint) / half_width

    def coefficient_at(self, time):
        """Return the selected bounded coefficient at ``time``."""
        return self.midpoint_coefficient + float(self.slope) * self.tau(time)

    def snapshot(self, time: float) -> Eq45SupportedVelocityCandidate:
        """Materialize an ordinary support-connected candidate at one scalar time."""
        time_arr = np.asarray(time, dtype=float)
        if time_arr.ndim != 0:
            raise ValueError("snapshot time must be a scalar")
        coefficient = float(self.coefficient_at(float(time_arr)))

        phi = list(self.profile_basis.phi_coefficients)
        swirl = list(self.profile_basis.swirl_coefficients)
        target = phi if self.family == "phi" else swirl
        target[self.mode_index] = coefficient

        basis = replace(
            self.profile_basis,
            phi_coefficients=tuple(phi),
            swirl_coefficients=tuple(swirl),
        )
        parent = replace(self.base.parent, profile_basis=basis)
        return Eq45SupportedVelocityCandidate(parent=parent, taper=self.base.taper)

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
        """Evaluate the temporal support-connected field as Cartesian ``[...,3]``."""
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
            raise RuntimeError("temporal support-connected Eq45 velocity became nonfinite")
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
        """Return ``(time,x,y,z,component)`` samples for visualization/export."""
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
                "family": self.family,
                "mode": [self.mode_i, self.mode_j],
                "slope": float(self.slope),
                "tau_definition": "interval_endpoints_map_to_minus1_plus1",
            },
            "truth_boundary": dict(_TRUTH_BOUNDARY),
        }

    @classmethod
    def from_dict(
        cls, data: Mapping[str, Any]
    ) -> "Eq45SupportedAffineTemporalModeCandidate":
        if not isinstance(data, Mapping):
            raise ValueError("temporal candidate payload must be an object")
        if data.get("schema") != SCHEMA:
            raise ValueError(f"schema must be {SCHEMA!r}")
        if data.get("classification") != _CLASSIFICATION:
            raise ValueError("temporal candidate source classification was modified")
        if data.get("truth_boundary") != _TRUTH_BOUNDARY:
            raise ValueError("temporal candidate truth-boundary metadata was modified")
        if "base_supported_candidate" not in data or "temporal_mode" not in data:
            raise ValueError("temporal candidate payload is missing a required component")

        base = Eq45SupportedVelocityCandidate.from_dict(data["base_supported_candidate"])
        if data.get("base_sha256") != base.sha256:
            raise ValueError("base_sha256 does not match the embedded supported candidate")

        temporal = data["temporal_mode"]
        if not isinstance(temporal, Mapping):
            raise ValueError("temporal_mode must be an object")
        if temporal.get("tau_definition") != "interval_endpoints_map_to_minus1_plus1":
            raise ValueError("temporal tau definition was modified")
        mode = temporal.get("mode")
        if not isinstance(mode, (list, tuple)) or len(mode) != 2:
            raise ValueError("temporal mode must contain exactly two indices")

        return cls(
            base=base,
            family=str(temporal["family"]),
            mode_i=int(mode[0]),
            mode_j=int(mode[1]),
            slope=float(temporal["slope"]),
        )

    def save_json(self, path: str | Path) -> None:
        Path(path).write_text(
            json.dumps(self.to_dict(), indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )

    @classmethod
    def load_json(
        cls, path: str | Path
    ) -> "Eq45SupportedAffineTemporalModeCandidate":
        return cls.from_dict(json.loads(Path(path).read_text(encoding="utf-8")))
