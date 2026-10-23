"""Bounded executable complete-curl oscillatory velocity correction.

Provenance boundary
-------------------
KokunoYumeto's corrected 2026-09-09 reconstruction writes oscillatory velocity
fields as complete curls of localized vector potentials and explicitly retains
the amplitude-derivative remainder.  It also records that applying a cutoff to
the potential produces the extra ``grad(cutoff) x A`` velocity term.

This module preserves those two structural requirements -- exact curl and
potential-level localization -- but it does *not* implement the source's full
cylindrical phase/frame/annular hierarchy.  The Cartesian phase, compact box
window and numerical parameter bounds below are autonomous reconstruction
choices.  Consequently fields produced here are divergence-free surrogate
corrections, not paper-exact Kokuno/OpenAI fields.

Source: KokunoYumeto, "Forced Navier-Stokes blowup: reconstruction and
validation reader (corrected 208-page edition)", version
2026.09.09-consolidated, DOI 10.5281/zenodo.22678406, especially the complete
curl identity in the oscillatory-field section and cutoff identity (A14).
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from math import pi
from typing import Any, Callable

import numpy as np


ArrayLike = Any
VelocityCallable = Callable[[ArrayLike, ArrayLike, ArrayLike, ArrayLike], np.ndarray]


def _triple(values: tuple[float, float, float], name: str) -> np.ndarray:
    array = np.asarray(values, dtype=float)
    if array.shape != (3,) or not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must contain exactly three finite values")
    return array


def _compact_c4_window(s: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return w(s)=(1-s^2)^5_+ and dw/ds.

    The zero extension is C4 because the polynomial has a fifth-order zero at
    |s|=1.  C4 regularity gives margin beyond the derivatives needed by the
    first velocity/vorticity diagnostics while keeping support exact.
    """
    inside = np.abs(s) < 1.0
    base = np.where(inside, 1.0 - s * s, 0.0)
    value = base**5
    derivative = np.where(inside, -10.0 * s * base**4, 0.0)
    return value, derivative


@dataclass(frozen=True)
class KokunoCompleteCurlCorrection:
    """One low-dimensional localized real wave written as an exact curl.

    The autonomous real vector potential is

        A = -amplitude * B(x) * q * sin(psi),
        q = (n x t) / |n|^2,
        psi = n . (x-center) + omega*t + phase,

    where ``t`` is the normalized projection of ``polarization`` onto the plane
    perpendicular to ``n=wave_vector``.  Its complete curl is evaluated
    analytically as

        u_osc = amplitude * [B*t*cos(psi) - (grad B x q)*sin(psi)].

    The second term is deliberately retained: dropping it would destroy the
    complete-curl/localization contract emphasized by the source.
    """

    amplitude: float = 0.25
    wave_vector: tuple[float, float, float] = (5.0, -3.0, 4.0)
    polarization: tuple[float, float, float] = (1.0, 2.0, 1.0)
    center: tuple[float, float, float] = (0.0, 0.0, 0.0)
    half_widths: tuple[float, float, float] = (0.9, 0.8, 0.7)
    omega: float = 2.0
    phase: float = 0.0

    def __post_init__(self) -> None:
        if not np.isfinite(self.amplitude) or abs(self.amplitude) > 2.0:
            raise ValueError("amplitude must be finite and lie in [-2, 2]")
        if not np.isfinite(self.omega) or abs(self.omega) > 16.0 * pi:
            raise ValueError("omega must be finite and lie in [-16*pi, 16*pi]")
        if not np.isfinite(self.phase) or abs(self.phase) > pi:
            raise ValueError("phase must be finite and lie in [-pi, pi]")

        n = _triple(self.wave_vector, "wave_vector")
        if np.max(np.abs(n)) > 32.0 or float(np.dot(n, n)) < 1.0e-12:
            raise ValueError("wave_vector components must lie in [-32, 32] and be nonzero")

        p = _triple(self.polarization, "polarization")
        transverse = p - n * (float(np.dot(p, n)) / float(np.dot(n, n)))
        if float(np.linalg.norm(transverse)) < 1.0e-10:
            raise ValueError("polarization must have a component transverse to wave_vector")

        _triple(self.center, "center")
        widths = _triple(self.half_widths, "half_widths")
        if np.any(widths < 0.05) or np.any(widths > 4.0):
            raise ValueError("half_widths must each lie in [0.05, 4]")

    @property
    def transverse_polarization(self) -> np.ndarray:
        n = _triple(self.wave_vector, "wave_vector")
        p = _triple(self.polarization, "polarization")
        p = p - n * (float(np.dot(p, n)) / float(np.dot(n, n)))
        return p / np.linalg.norm(p)

    @property
    def potential_direction(self) -> np.ndarray:
        n = _triple(self.wave_vector, "wave_vector")
        t = self.transverse_polarization
        return np.cross(n, t) / float(np.dot(n, n))

    def _broadcast(self, x, y, z, time):
        return np.broadcast_arrays(
            np.asarray(x, dtype=float),
            np.asarray(y, dtype=float),
            np.asarray(z, dtype=float),
            np.asarray(time, dtype=float),
        )

    def _envelope(self, x, y, z, time) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        x, y, z, time = self._broadcast(x, y, z, time)
        coords = np.stack((x, y, z), axis=-1)
        center = _triple(self.center, "center")
        widths = _triple(self.half_widths, "half_widths")
        shifted = coords - center

        values: list[np.ndarray] = []
        derivatives: list[np.ndarray] = []
        for axis in range(3):
            value, derivative_s = _compact_c4_window(shifted[..., axis] / widths[axis])
            values.append(value)
            derivatives.append(derivative_s / widths[axis])

        envelope = values[0] * values[1] * values[2]
        gradient = np.stack(
            (
                derivatives[0] * values[1] * values[2],
                values[0] * derivatives[1] * values[2],
                values[0] * values[1] * derivatives[2],
            ),
            axis=-1,
        )
        return shifted, envelope, gradient

    def _phase(self, shifted: np.ndarray, time: np.ndarray) -> np.ndarray:
        n = _triple(self.wave_vector, "wave_vector")
        return np.einsum("...i,i->...", shifted, n) + self.omega * time + self.phase

    def vector_potential(self, x, y, z, time) -> np.ndarray:
        """Evaluate the localized vector potential in ``(...,3)`` layout."""
        _, _, _, time = self._broadcast(x, y, z, time)
        shifted, envelope, _ = self._envelope(x, y, z, time)
        psi = self._phase(shifted, time)
        q = self.potential_direction
        return -self.amplitude * envelope[..., None] * np.sin(psi)[..., None] * q

    def velocity(self, x, y, z, time) -> np.ndarray:
        """Evaluate the exact analytic complete curl in ``(...,3)`` layout."""
        _, _, _, time = self._broadcast(x, y, z, time)
        shifted, envelope, grad_envelope = self._envelope(x, y, z, time)
        psi = self._phase(shifted, time)
        t = self.transverse_polarization
        q = self.potential_direction

        leading = envelope[..., None] * np.cos(psi)[..., None] * t
        cutoff_remainder = -np.cross(grad_envelope, q) * np.sin(psi)[..., None]
        return self.amplitude * (leading + cutoff_remainder)

    __call__ = velocity

    def at_points(self, points: ArrayLike, time: ArrayLike) -> np.ndarray:
        points = np.asarray(points, dtype=float)
        if points.ndim < 1 or points.shape[-1] != 3:
            raise ValueError("points must have final dimension 3")
        return self.velocity(points[..., 0], points[..., 1], points[..., 2], time)

    def metadata(self) -> dict[str, Any]:
        return {
            "family": "kokuno_complete_curl_cartesian_surrogate_v1",
            "source": {
                "author": "KokunoYumeto",
                "edition": "2026.09.09-consolidated",
                "doi": "10.5281/zenodo.22678406",
                "source_structure_used": [
                    "velocity represented as complete curl of vector potential",
                    "potential-level localization retains cutoff-gradient curl term",
                ],
            },
            "autonomous_choices": [
                "Cartesian affine phase",
                "separable compact C4 box window",
                "bounded amplitude/frequency/orientation/width parameterization",
            ],
            "paper_exact": False,
            "openai_field_identified": False,
            "parameters": asdict(self),
        }


@dataclass(frozen=True)
class CompositeVelocityField:
    """Compose an existing public velocity callable with one curl correction."""

    base_velocity: VelocityCallable
    correction: KokunoCompleteCurlCorrection

    def velocity(self, x, y, z, time) -> np.ndarray:
        base = np.asarray(self.base_velocity(x, y, z, time), dtype=float)
        delta = self.correction.velocity(x, y, z, time)
        if base.shape != delta.shape or base.shape[-1] != 3:
            raise ValueError("base velocity and correction must both use broadcast (...,3) layout")
        return base + delta

    __call__ = velocity


def compose_velocity(
    base_velocity: VelocityCallable,
    correction: KokunoCompleteCurlCorrection,
) -> CompositeVelocityField:
    """Return ``u_base + u_osc`` without changing either serialized candidate."""
    return CompositeVelocityField(base_velocity=base_velocity, correction=correction)
