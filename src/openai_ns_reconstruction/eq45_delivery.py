"""Named public delivery wrapper for the frozen Eq. (4.5) candidate.

This module deliberately does not replace ``velocity_components.velocity``,
which still exposes the older independently reconstructed candidate.  It gives
callers an explicit Eq. (4.5) entry point while delegating every numerical
velocity evaluation to :class:`Eq45VelocityCandidate`.

The default candidate is the deterministic checked seed.  It is callable and
serializable, but its physical-support connection, visual correspondence and
Navier--Stokes validation remain unverified.
"""
from __future__ import annotations

from functools import lru_cache
from pathlib import Path
from typing import Any

import numpy as np

from .constrained_eq45_candidate import Eq45VelocityCandidate


class Eq45DeliveryField:
    """Small user-facing wrapper around one frozen ``Eq45VelocityCandidate``."""

    def __init__(
        self,
        candidate: Eq45VelocityCandidate | None = None,
        *,
        candidate_path: str | Path | None = None,
    ) -> None:
        if candidate is not None and candidate_path is not None:
            raise ValueError("provide candidate or candidate_path, not both")
        if candidate_path is not None:
            candidate = Eq45VelocityCandidate.load_json(candidate_path)
        if candidate is None:
            candidate = Eq45VelocityCandidate.seed()
        if not isinstance(candidate, Eq45VelocityCandidate):
            raise TypeError("candidate must be Eq45VelocityCandidate")
        self._candidate = candidate

    @property
    def candidate(self) -> Eq45VelocityCandidate:
        return self._candidate

    @property
    def sha256(self) -> str:
        return self._candidate.sha256

    def velocity(self, x, y, z, t) -> np.ndarray:
        """Return broadcast Cartesian velocity in ``(..., 3)`` component order."""
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
        """Return public point-batch velocity in ``(..., 3)`` order."""
        return self._candidate.at_points(points, time)

    def grid(self, x, y, z, times) -> np.ndarray:
        """Return ``(time,x,y,z,component)`` samples from the same evaluator."""
        return self._candidate.grid(x, y, z, times)

    def save_candidate(self, path: str | Path) -> None:
        """Persist exactly the wrapped candidate JSON."""
        self._candidate.save_json(path)

    @classmethod
    def load_candidate(cls, path: str | Path) -> "Eq45DeliveryField":
        return cls(candidate_path=path)

    def metadata(self) -> dict[str, Any]:
        payload = self._candidate.to_dict()
        return {
            "family": "eq45_velocity_candidate_v1",
            "candidate_sha256": self.sha256,
            "entrypoint": "openai_ns_reconstruction.eq45_delivery:velocity",
            "components": ["u", "v", "w"],
            "coordinates": "right-handed Cartesian x,y,z",
            "units": "dimensionless; no physical unit or OpenAI scene calibration claimed",
            "time_interval": list(payload["time_interval"]),
            "grid_layout": ["time", "x", "y", "z", "component"],
            "classification": dict(payload["classification"]),
            "truth_boundary": dict(payload["truth_boundary"]),
        }


@lru_cache(maxsize=1)
def default_field() -> Eq45DeliveryField:
    """Return the deterministic named Eq. (4.5) delivery field."""
    return Eq45DeliveryField()


def velocity(x, y, z, t) -> np.ndarray:
    """Literal named Eq. (4.5) interface ``velocity(x,y,z,t)->[...,3]``."""
    return default_field().velocity(x, y, z, t)


def components(x, y, z, t):
    return default_field().components(x, y, z, t)


def u(x, y, z, t):
    return default_field().u(x, y, z, t)


def v(x, y, z, t):
    return default_field().v(x, y, z, t)


def w(x, y, z, t):
    return default_field().w(x, y, z, t)
