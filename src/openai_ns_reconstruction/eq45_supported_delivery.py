"""Named public delivery wrapper for the support-connected Eq. (4.5) field.

This module exposes the post-taper child candidate through the literal engineering
interface ``velocity(x,y,z,t)->[...,3]``.  It does not promote visualization,
PDE, paper-exact, or OpenAI-field claims; those states remain whatever the
serialized child candidate declares.
"""
from __future__ import annotations

from functools import lru_cache
from pathlib import Path
from typing import Any

import numpy as np

from .constrained_eq45_candidate import Eq45VelocityCandidate
from .constrained_eq45_supported_candidate import Eq45SupportedVelocityCandidate


class Eq45SupportedDeliveryField:
    """User-facing wrapper around one support-connected Eq45 child candidate."""

    def __init__(
        self,
        candidate: Eq45SupportedVelocityCandidate | None = None,
        *,
        candidate_path: str | Path | None = None,
    ) -> None:
        if candidate is not None and candidate_path is not None:
            raise ValueError("provide candidate or candidate_path, not both")
        if candidate_path is not None:
            candidate = Eq45SupportedVelocityCandidate.load_json(candidate_path)
        if candidate is None:
            candidate = Eq45SupportedVelocityCandidate(parent=Eq45VelocityCandidate.seed())
        if not isinstance(candidate, Eq45SupportedVelocityCandidate):
            raise TypeError("candidate must be Eq45SupportedVelocityCandidate")
        self._candidate = candidate

    @property
    def candidate(self) -> Eq45SupportedVelocityCandidate:
        return self._candidate

    @property
    def sha256(self) -> str:
        return self._candidate.sha256

    def velocity(self, x, y, z, t) -> np.ndarray:
        """Return broadcast Cartesian velocity in ``(...,3)`` component order."""
        return self._candidate.velocity_xyz(x, y, z, t)

    __call__ = velocity

    def components(self, x, y, z, t):
        """Return ``(u,v,w)`` with the same broadcast shape as the inputs."""
        values = self.velocity(x, y, z, t)
        parts = tuple(values[..., i] for i in range(3))
        if values.ndim == 1:
            return tuple(float(value) for value in parts)
        return parts

    def u(self, x, y, z, t):
        return self.components(x, y, z, t)[0]

    def v(self, x, y, z, t):
        return self.components(x, y, z, t)[1]

    def w(self, x, y, z, t):
        return self.components(x, y, z, t)[2]

    def at_points(self, points, time) -> np.ndarray:
        """Return point-batch velocity through the same support-connected evaluator."""
        return self._candidate.at_points(points, time)

    def grid(self, x, y, z, times) -> np.ndarray:
        """Return ``(time,x,y,z,component)`` samples for visualization/export."""
        return self._candidate.grid(x, y, z, times)

    def save_candidate(self, path: str | Path) -> None:
        self._candidate.save_json(path)

    @classmethod
    def load_candidate(cls, path: str | Path) -> "Eq45SupportedDeliveryField":
        return cls(candidate_path=path)

    def metadata(self) -> dict[str, Any]:
        payload = self._candidate.to_dict()
        return {
            "family": payload["schema"],
            "candidate_sha256": self.sha256,
            "parent_sha256": self._candidate.parent_sha256,
            "entrypoint": "openai_ns_reconstruction.eq45_supported_delivery:velocity",
            "components": ["u", "v", "w"],
            "coordinates": "right-handed Cartesian x,y,z",
            "units": "dimensionless; no physical unit or OpenAI scene calibration claimed",
            "time_interval": [self._candidate.time_start, self._candidate.time_end],
            "grid_layout": ["time", "x", "y", "z", "component"],
            "classification": dict(payload["classification"]),
            "truth_boundary": dict(payload["truth_boundary"]),
        }


@lru_cache(maxsize=1)
def default_field() -> Eq45SupportedDeliveryField:
    """Return the deterministic default support-connected Eq45 delivery field."""
    return Eq45SupportedDeliveryField()


def velocity(x, y, z, t) -> np.ndarray:
    """Literal supported Eq45 interface ``velocity(x,y,z,t)->[...,3]``."""
    return default_field().velocity(x, y, z, t)


def components(x, y, z, t):
    return default_field().components(x, y, z, t)


def u(x, y, z, t):
    return default_field().u(x, y, z, t)


def v(x, y, z, t):
    return default_field().v(x, y, z, t)


def w(x, y, z, t):
    return default_field().w(x, y, z, t)
