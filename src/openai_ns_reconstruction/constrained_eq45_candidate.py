"""Serializable callable candidate that composes the integrated Eq. (4.5) pieces.

This module is intentionally narrow: it connects the already-integrated
source-correct q solver / Cartesian Eq. (4.5) backbone, the bounded autonomous
Phi/F basis, and the divergence-preserving Phi -> (v0,U) adapter into one
reproducible ``velocity(points, time) -> [u,v,w]`` object.

The resulting field is a visualization/research candidate.  The selected
profile coefficients, h value and finite delivery interval are autonomous
project choices.  The object does not claim the hidden OpenAI profiles, the
complete paper field, physical compact-support validation, Navier--Stokes
validation, or blow-up.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from .constrained_eq45_poloidal_streamfunction import (
    eq45_streamfunction_poloidal_jet,
)
from .constrained_eq45_profile_basis import Eq45CompactProfileBasis
from .constrained_eq45_velocity import Eq45VelocityBackbone


SCHEMA = "eq45_velocity_candidate_v1"
DEFAULT_H = 0.005
DEFAULT_TIME_INTERVAL = (0.25, 0.75)
MIN_COEFFICIENT_NORM = 1e-12

_CLASSIFICATION = {
    "velocity_interface": "user_requirement",
    "eq45_backbone": "public_source_fact",
    "profile_family": "autonomous_design",
    "profile_coefficients": "autonomous_design",
    "h": "autonomous_design",
    "time_interval": "autonomous_design",
    "openai_hidden_profiles": "pending_unknown",
}

_TRUTH_BOUNDARY = {
    "callable_serializable": True,
    "velocity_export_ready": True,
    "physical_support_validated": False,
    "visualization_ready": False,
    "visual_correspondence_verified": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
    "requires_physical_support_connection": True,
}


@dataclass(frozen=True)
class Eq45VelocityCandidate:
    """Bounded, serializable Eq. (4.5) leading-field candidate.

    ``Eq45CompactProfileBasis`` supplies one scalar streamfunction profile
    ``Phi`` and one independent swirl profile ``F``.  The streamfunction adapter
    derives the coupled ``v0`` and ``U`` profiles before the public Eq. (4.5)
    Cartesian backbone evaluates the final velocity.

    The current object deliberately stops before the open physical-support taper
    lane.  It is therefore callable/exportable but must not be described as
    satisfying the registered physical compact-support contract.
    """

    profile_basis: Eq45CompactProfileBasis
    h: float = DEFAULT_H
    time_start: float = DEFAULT_TIME_INTERVAL[0]
    time_end: float = DEFAULT_TIME_INTERVAL[1]

    def __post_init__(self) -> None:
        if not isinstance(self.profile_basis, Eq45CompactProfileBasis):
            raise TypeError("profile_basis must be Eq45CompactProfileBasis")
        if not np.isscalar(self.h) or not np.isfinite(self.h) or not (0.0 < float(self.h) < 0.5):
            raise ValueError("h must be finite with 0 < h < 1/2")
        if not np.isfinite(self.time_start) or not np.isfinite(self.time_end):
            raise ValueError("time interval must be finite")
        if not (self.time_start < self.time_end < 1.0):
            raise ValueError("time interval must satisfy start < end < 1")

        coefficients = np.asarray(
            self.profile_basis.phi_coefficients + self.profile_basis.swirl_coefficients,
            dtype=float,
        )
        if float(np.linalg.norm(coefficients)) <= MIN_COEFFICIENT_NORM:
            raise ValueError("candidate profile coefficients must be nontrivial")

    @classmethod
    def seed(cls) -> "Eq45VelocityCandidate":
        """Return the deterministic autonomous nonzero visualization seed."""
        return cls(profile_basis=Eq45CompactProfileBasis.seed())

    def profile_values(self, X, eta) -> np.ndarray:
        """Return vectorized ``[..., (v0,F,U)]`` values for the backbone."""
        jets = self.profile_basis.evaluate(X, eta)
        poloidal = eq45_streamfunction_poloidal_jet(
            X,
            eta,
            self.h,
            jets.phi,
            jets.phi_x,
            jets.phi_eta,
            jets.phi_xx,
            jets.phi_xeta,
        )
        values = np.stack((poloidal.v0, jets.swirl, poloidal.U), axis=-1)
        if not np.all(np.isfinite(values)):
            raise RuntimeError("Eq. (4.5) profile composition became nonfinite")
        return values

    def _validated_time(self, points, time) -> np.ndarray:
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
        if np.any((time_arr < self.time_start) | (time_arr > self.time_end)):
            raise ValueError(
                f"time must lie in declared interval [{self.time_start}, {self.time_end}]"
            )
        return time_arr

    def velocity(self, points, time) -> np.ndarray:
        """Evaluate the composed Eq. (4.5) candidate at arbitrary point batches."""
        time_arr = self._validated_time(points, time)
        backbone = Eq45VelocityBackbone(self.profile_values, h=self.h)
        return backbone.velocity(points, time_arr)

    def at_points(self, points, time) -> np.ndarray:
        """Public-style alias used by downstream validators/visualizers."""
        return self.velocity(points, time)

    def velocity_xyz(self, x, y, z, time) -> np.ndarray:
        """Broadcast x/y/z/time and return the final ``[...,3]`` velocity."""
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
        """Return direct velocity samples in ``(time,x,y,z,component)`` order.

        This mirrors the repository's existing public visualization-grid layout
        while evaluating this candidate directly, so saved/reloaded Eq45
        candidates can feed MATLAB/Python grid consumers without accessing
        optimizer internals or a different velocity implementation.
        """
        axes = [np.asarray(value, dtype=float) for value in (times, x, y, z)]
        if any(axis.ndim != 1 or axis.size == 0 for axis in axes):
            raise ValueError("grid axes and times must be nonempty 1D arrays")
        if not all(np.all(np.isfinite(axis)) for axis in axes):
            raise ValueError("grid axes and times must be finite")
        tt, xx, yy, zz = np.meshgrid(*axes, indexing="ij")
        return self.velocity_xyz(xx, yy, zz, tt)

    __call__ = velocity

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
            "h": float(self.h),
            "time_interval": [float(self.time_start), float(self.time_end)],
            "profile_basis": self.profile_basis.to_dict(),
            "truth_boundary": dict(_TRUTH_BOUNDARY),
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> "Eq45VelocityCandidate":
        if not isinstance(data, Mapping):
            raise ValueError("candidate payload must be an object")
        if data.get("schema") != SCHEMA:
            raise ValueError(f"schema must be {SCHEMA!r}")
        if data.get("classification") != _CLASSIFICATION:
            raise ValueError("candidate source classification does not match the governed contract")
        truth = data.get("truth_boundary")
        if truth != _TRUTH_BOUNDARY:
            raise ValueError("candidate truth-boundary metadata was modified")
        interval = data.get("time_interval")
        if not isinstance(interval, (list, tuple)) or len(interval) != 2:
            raise ValueError("time_interval must contain exactly two values")
        if "profile_basis" not in data:
            raise ValueError("candidate payload is missing profile_basis")
        return cls(
            profile_basis=Eq45CompactProfileBasis.from_dict(data["profile_basis"]),
            h=float(data["h"]),
            time_start=float(interval[0]),
            time_end=float(interval[1]),
        )

    def save_json(self, path: str | Path) -> None:
        target = Path(path)
        target.write_text(
            json.dumps(self.to_dict(), indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )

    @classmethod
    def load_json(cls, path: str | Path) -> "Eq45VelocityCandidate":
        target = Path(path)
        return cls.from_dict(json.loads(target.read_text(encoding="utf-8")))
