"""Vectorized evaluator for the public-source Eq. (4.5) leading velocity backbone.

The source-governed similarity relation used here is

    q - z**2 * q**(2*h) = 1 - t,

not the stale negative-exponent transcription. This module evaluates the
leading Cartesian mixing once a profile provider supplies ``(v0, F, U)`` at
``(X, eta)``. Profile identification/parameterization is deliberately kept
separate so that later Phi/F streamfunction coupling can replace the provider
without changing the public velocity evaluator.

This is a leading-field visualization backbone. It is not a claim of the full
paper field, OpenAI's hidden numerical profile data, PDE validation, or blow-up.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np

from .constrained_eq45_q import solve_eq45_q


ProfileProvider = Callable[[np.ndarray, np.ndarray], np.ndarray]


@dataclass(frozen=True)
class Eq45Coordinates:
    """Broadcast Eq. (4.5) similarity coordinates for one point/time batch."""

    q: np.ndarray
    X: np.ndarray
    eta: np.ndarray


class Eq45VelocityBackbone:
    """Evaluate Eq. (4.5) from a supplied two-dimensional profile provider.

    The provider receives ``(X, eta)`` arrays and must return an array with
    shape ``X.shape + (3,)`` ordered as ``[..., 0]=v0``, ``[..., 1]=F``, and
    ``[..., 2]=U``. The profile provider is an autonomous numerical interface,
    not evidence that the profiles are the paper's final numerical data.
    """

    def __init__(self, profile_provider: ProfileProvider, *, h: float = 0.005):
        if not callable(profile_provider):
            raise TypeError("profile_provider must be callable")
        if not np.isscalar(h) or not np.isfinite(h) or not (0.0 < float(h) < 0.5):
            raise ValueError("h must be a finite scalar with 0 < h < 1/2")
        self.profile_provider = profile_provider
        self.h = float(h)

    @staticmethod
    def _prepare(points, time):
        points = np.asarray(points, dtype=float)
        if points.ndim == 0 or points.shape[-1] != 3:
            raise ValueError("points must have shape (..., 3)")
        if not np.all(np.isfinite(points)):
            raise ValueError("points must be finite")

        lead_shape = points.shape[:-1]
        time = np.asarray(time, dtype=float)
        if not np.all(np.isfinite(time)):
            raise ValueError("time must be finite")
        try:
            time = np.broadcast_to(time, lead_shape)
        except ValueError as exc:
            raise ValueError("time must broadcast to points.shape[:-1]") from exc
        if np.any(time >= 1.0):
            raise ValueError("Eq. (4.5) backbone requires time < 1")
        return points, time

    def coordinates(self, points, time) -> Eq45Coordinates:
        """Return vectorized ``q, X, eta`` for a point/time batch."""

        points, time = self._prepare(points, time)
        x = points[..., 0]
        y = points[..., 1]
        z = points[..., 2]

        q = solve_eq45_q(z, time, self.h)
        D = 0.5 - self.h
        X = (x * x + y * y) / (2.0 * q)
        eta = z / np.power(q, D)

        if not (np.all(np.isfinite(X)) and np.all(np.isfinite(eta))):
            raise RuntimeError("Eq. (4.5) coordinates became nonfinite")
        if np.any(X < 0.0):
            raise RuntimeError("Eq. (4.5) radial coordinate X became negative")
        return Eq45Coordinates(q=q, X=X, eta=eta)

    def velocity(self, points, time) -> np.ndarray:
        """Return Eq. (4.5) leading velocity with shape ``points.shape``."""

        points, time = self._prepare(points, time)
        coordinates = self.coordinates(points, time)
        values = np.asarray(
            self.profile_provider(coordinates.X, coordinates.eta), dtype=float
        )
        expected_shape = coordinates.X.shape + (3,)
        if values.shape != expected_shape:
            raise ValueError(
                "profile_provider must return shape "
                f"{expected_shape}, got {values.shape}"
            )
        if not np.all(np.isfinite(values)):
            raise ValueError("profile_provider returned nonfinite values")

        v0 = values[..., 0]
        F = values[..., 1]
        U = values[..., 2]
        x = points[..., 0]
        y = points[..., 1]
        q = coordinates.q
        h = self.h

        swirl_scale = np.power(q, -1.0 - h)
        axial_scale = np.power(q, -0.5 - h)
        u = x * v0 / (2.0 * q) - y * swirl_scale * F
        v = y * v0 / (2.0 * q) + x * swirl_scale * F
        w = axial_scale * U
        out = np.stack((u, v, w), axis=-1)
        if not np.all(np.isfinite(out)):
            raise RuntimeError("Eq. (4.5) velocity became nonfinite")
        return out

    __call__ = velocity
