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


def _compact_c4_window_derivatives(
    s: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return ``w=(1-s^2)^5_+`` and its first three ``s`` derivatives.

    The zero extension is C4 because the polynomial has a fifth-order zero at
    ``|s|=1``.  The third derivative is needed only for the analytic Laplacian
    of the complete curl; no derivative of an indicator is introduced because
    all returned derivatives vanish continuously at the support faces.
    """
    inside = np.abs(s) < 1.0
    base = np.where(inside, 1.0 - s * s, 0.0)
    value = base**5
    first = np.where(inside, -10.0 * s * base**4, 0.0)
    second = np.where(inside, 10.0 * base**3 * (9.0 * s * s - 1.0), 0.0)
    third = np.where(inside, 240.0 * s * base**2 * (1.0 - 3.0 * s * s), 0.0)
    return value, first, second, third


def _compact_c4_window(s: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return the C4 compact window and its first derivative."""
    value, first, _, _ = _compact_c4_window_derivatives(s)
    return value, first


@dataclass(frozen=True)
class KokunoCompleteCurlCorrection:
    """One low-dimensional localized real wave written as an exact curl.

    The autonomous real vector potential is

        A = -amplitude * B(x) * q * sin(psi),
        q = (n x t) / |n|^2,
        psi = n . (x-center) + omega*time + phase,

    where ``t`` is the normalized projection of ``polarization`` onto the plane
    perpendicular to ``n=wave_vector``.  Its complete curl is evaluated
    analytically as

        u_osc = amplitude * [B*t*cos(psi) - (grad B x q)*sin(psi)].

    The second term is deliberately retained: dropping it would destroy the
    complete-curl/localization contract emphasized by the source.

    The analytic differential methods in this class differentiate this same
    autonomous surrogate exactly; they do not add source-derived phase data.
    They are exposed so later residual code can consume the correction without
    numerically differentiating the oscillatory field itself.
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

    def _envelope_differential_data(
        self, x, y, z, time
    ) -> tuple[
        np.ndarray,
        np.ndarray,
        np.ndarray,
        np.ndarray,
        np.ndarray,
        np.ndarray,
    ]:
        """Return shifted coordinates and B, grad B, Hess B, Lap B, grad Lap B."""
        x, y, z, time = self._broadcast(x, y, z, time)
        coords = np.stack((x, y, z), axis=-1)
        center = _triple(self.center, "center")
        widths = _triple(self.half_widths, "half_widths")
        shifted = coords - center

        derivatives: list[tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]] = []
        for axis in range(3):
            value, first, second, third = _compact_c4_window_derivatives(
                shifted[..., axis] / widths[axis]
            )
            derivatives.append(
                (
                    value,
                    first / widths[axis],
                    second / widths[axis] ** 2,
                    third / widths[axis] ** 3,
                )
            )

        fx, fx1, fx2, fx3 = derivatives[0]
        fy, fy1, fy2, fy3 = derivatives[1]
        fz, fz1, fz2, fz3 = derivatives[2]

        envelope = fx * fy * fz
        gradient = np.stack(
            (fx1 * fy * fz, fx * fy1 * fz, fx * fy * fz1), axis=-1
        )

        hxx = fx2 * fy * fz
        hyy = fx * fy2 * fz
        hzz = fx * fy * fz2
        hxy = fx1 * fy1 * fz
        hxz = fx1 * fy * fz1
        hyz = fx * fy1 * fz1
        hessian = np.stack(
            (
                np.stack((hxx, hxy, hxz), axis=-1),
                np.stack((hxy, hyy, hyz), axis=-1),
                np.stack((hxz, hyz, hzz), axis=-1),
            ),
            axis=-2,
        )

        laplacian = hxx + hyy + hzz
        grad_laplacian = np.stack(
            (
                fx3 * fy * fz + fx1 * fy2 * fz + fx1 * fy * fz2,
                fx2 * fy1 * fz + fx * fy3 * fz + fx * fy1 * fz2,
                fx2 * fy * fz1 + fx * fy2 * fz1 + fx * fy * fz3,
            ),
            axis=-1,
        )
        return shifted, envelope, gradient, hessian, laplacian, grad_laplacian

    def _envelope(self, x, y, z, time) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        shifted, envelope, gradient, _, _, _ = self._envelope_differential_data(
            x, y, z, time
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

    def time_derivative(self, x, y, z, time) -> np.ndarray:
        """Return analytic ``partial_t u_osc`` in ``(...,3)`` layout."""
        _, _, _, time = self._broadcast(x, y, z, time)
        shifted, envelope, grad_envelope = self._envelope(x, y, z, time)
        psi = self._phase(shifted, time)
        t = self.transverse_polarization
        q = self.potential_direction
        grad_cross_q = np.cross(grad_envelope, q)
        return self.amplitude * self.omega * (
            -envelope[..., None] * np.sin(psi)[..., None] * t
            - grad_cross_q * np.cos(psi)[..., None]
        )

    def spatial_jacobian(self, x, y, z, time) -> np.ndarray:
        """Return analytic ``du_i/dx_j`` with shape ``(...,3,3)``."""
        _, _, _, time = self._broadcast(x, y, z, time)
        shifted, envelope, grad_envelope, hessian, _, _ = self._envelope_differential_data(
            x, y, z, time
        )
        psi = self._phase(shifted, time)
        sine = np.sin(psi)
        cosine = np.cos(psi)
        n = _triple(self.wave_vector, "wave_vector")
        t = self.transverse_polarization
        q = self.potential_direction
        grad_cross_q = np.cross(grad_envelope, q)

        columns = []
        for axis in range(3):
            d_grad_cross_q = np.cross(hessian[..., :, axis], q)
            column = (
                grad_envelope[..., axis, None] * cosine[..., None] * t
                - envelope[..., None] * n[axis] * sine[..., None] * t
                - d_grad_cross_q * sine[..., None]
                - grad_cross_q * n[axis] * cosine[..., None]
            )
            columns.append(self.amplitude * column)
        return np.stack(columns, axis=-1)

    def divergence(self, x, y, z, time) -> np.ndarray:
        """Return the analytic trace of the spatial Jacobian."""
        jacobian = self.spatial_jacobian(x, y, z, time)
        return np.trace(jacobian, axis1=-2, axis2=-1)

    def vorticity(self, x, y, z, time) -> np.ndarray:
        """Return analytic ``curl u_osc`` in Cartesian component order."""
        jacobian = self.spatial_jacobian(x, y, z, time)
        return np.stack(
            (
                jacobian[..., 2, 1] - jacobian[..., 1, 2],
                jacobian[..., 0, 2] - jacobian[..., 2, 0],
                jacobian[..., 1, 0] - jacobian[..., 0, 1],
            ),
            axis=-1,
        )

    def laplacian(self, x, y, z, time) -> np.ndarray:
        """Return analytic componentwise ``Delta u_osc`` in ``(...,3)`` layout."""
        _, _, _, time = self._broadcast(x, y, z, time)
        (
            shifted,
            envelope,
            grad_envelope,
            hessian,
            lap_envelope,
            grad_lap_envelope,
        ) = self._envelope_differential_data(x, y, z, time)
        psi = self._phase(shifted, time)
        sine = np.sin(psi)
        cosine = np.cos(psi)
        n = _triple(self.wave_vector, "wave_vector")
        n2 = float(np.dot(n, n))
        t = self.transverse_polarization
        q = self.potential_direction

        grad_dot_n = np.einsum("...i,i->...", grad_envelope, n)
        hessian_n = np.einsum("...ij,j->...i", hessian, n)
        grad_cross_q = np.cross(grad_envelope, q)
        grad_lap_cross_q = np.cross(grad_lap_envelope, q)
        hessian_n_cross_q = np.cross(hessian_n, q)

        leading_scalar = (
            (lap_envelope - n2 * envelope) * cosine
            - 2.0 * grad_dot_n * sine
        )
        remainder = (
            (grad_lap_cross_q - n2 * grad_cross_q) * sine[..., None]
            + 2.0 * hessian_n_cross_q * cosine[..., None]
        )
        return self.amplitude * (leading_scalar[..., None] * t - remainder)

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
            "analytic_surrogate_derivatives": [
                "partial_t velocity",
                "Cartesian spatial Jacobian",
                "vorticity",
                "componentwise Laplacian",
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
