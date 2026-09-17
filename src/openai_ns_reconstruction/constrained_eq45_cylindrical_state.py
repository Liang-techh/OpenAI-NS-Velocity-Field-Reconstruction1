"""Axis-regular cylindrical state adapter for the production Eq. (4.5) candidate.

The integrated :class:`Eq45VelocityCandidate` already evaluates the public
Eq. (4.5) Cartesian velocity.  The open physical-support taper lane, however,
needs the *same* candidate expressed as the associated axisymmetric
streamfunction and cylindrical components.  This module exposes that narrow
representation seam without changing profile coefficients, forcing, pressure,
or any readiness claim.

For the candidate streamfunction

    psi = q**(1/2-h) * X * Phi(X, eta),

the returned components are

    u_r     = r/(2q) * v0,
    u_theta = r*q**(-1-h) * F,
    u_z     = q**(-1/2-h) * U.

These are exactly the cylindrical form of the existing Cartesian backbone and
are suitable inputs for a divergence-preserving streamfunction-level exterior
connection.  They are not a new velocity model and do not establish physical
support, visual correspondence, PDE validity, or paper exactness.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .constrained_eq45_candidate import Eq45VelocityCandidate
from .constrained_eq45_poloidal_streamfunction import eq45_streamfunction_poloidal_jet
from .constrained_eq45_velocity import Eq45VelocityBackbone


@dataclass(frozen=True)
class Eq45CylindricalState:
    """Broadcast similarity coordinates plus streamfunction/cylindrical velocity."""

    q: np.ndarray
    X: np.ndarray
    eta: np.ndarray
    psi: np.ndarray
    u_r: np.ndarray
    u_theta: np.ndarray
    u_z: np.ndarray

    def cartesian(self, points, *, axis_tolerance: float = 1e-12) -> np.ndarray:
        """Reconstruct ``[...,3]`` Cartesian velocity from this state."""
        points = np.asarray(points, dtype=float)
        if points.ndim == 0 or points.shape[-1] != 3:
            raise ValueError("points must have shape (..., 3)")
        if not np.all(np.isfinite(points)):
            raise ValueError("points must be finite")
        if not np.isfinite(axis_tolerance) or axis_tolerance < 0.0:
            raise ValueError("axis_tolerance must be finite and nonnegative")

        lead = points.shape[:-1]
        for name in ("q", "X", "eta", "psi", "u_r", "u_theta", "u_z"):
            value = np.asarray(getattr(self, name), dtype=float)
            if value.shape != lead:
                raise ValueError(f"state {name} must have shape {lead}")
            if not np.all(np.isfinite(value)):
                raise ValueError(f"state {name} must be finite")

        x = points[..., 0]
        y = points[..., 1]
        r = np.hypot(x, y)
        axis = r <= axis_tolerance
        cos_theta = np.divide(x, r, out=np.ones_like(r), where=~axis)
        sin_theta = np.divide(y, r, out=np.zeros_like(r), where=~axis)

        u = self.u_r * cos_theta - self.u_theta * sin_theta
        v = self.u_r * sin_theta + self.u_theta * cos_theta
        u = np.where(axis, 0.0, u)
        v = np.where(axis, 0.0, v)
        out = np.stack((u, v, self.u_z), axis=-1)
        if not np.all(np.isfinite(out)):
            raise RuntimeError("cylindrical-to-Cartesian reconstruction became nonfinite")
        return out

    def taper_inputs(self) -> dict[str, np.ndarray]:
        """Return the four named arrays consumed by the physical taper primitive."""
        return {
            "psi": self.psi,
            "u_r": self.u_r,
            "u_theta": self.u_theta,
            "u_z": self.u_z,
        }


def eq45_candidate_cylindrical_state(
    candidate: Eq45VelocityCandidate,
    points,
    time,
) -> Eq45CylindricalState:
    """Evaluate the production candidate as ``(psi,u_r,u_theta,u_z)``.

    Time validation and the similarity-coordinate solver are reused from the
    production candidate/backbone.  No independent profile fit or alternate
    velocity formula is introduced.
    """
    if not isinstance(candidate, Eq45VelocityCandidate):
        raise TypeError("candidate must be Eq45VelocityCandidate")

    points_arr = np.asarray(points, dtype=float)
    time_arr = candidate._validated_time(points_arr, time)
    backbone = Eq45VelocityBackbone(candidate.profile_values, h=candidate.h)
    coordinates = backbone.coordinates(points_arr, time_arr)

    jets = candidate.profile_basis.evaluate(coordinates.X, coordinates.eta)
    poloidal = eq45_streamfunction_poloidal_jet(
        coordinates.X,
        coordinates.eta,
        candidate.h,
        jets.phi,
        jets.phi_x,
        jets.phi_eta,
        jets.phi_xx,
        jets.phi_xeta,
    )

    x = points_arr[..., 0]
    y = points_arr[..., 1]
    r = np.hypot(x, y)
    q = coordinates.q
    h = candidate.h

    psi = np.power(q, 0.5 - h) * coordinates.X * jets.phi
    u_r = r * poloidal.v0 / (2.0 * q)
    u_theta = r * np.power(q, -1.0 - h) * jets.swirl
    u_z = np.power(q, -0.5 - h) * poloidal.U

    arrays = (q, coordinates.X, coordinates.eta, psi, u_r, u_theta, u_z)
    if not all(np.all(np.isfinite(value)) for value in arrays):
        raise RuntimeError("Eq45 cylindrical state became nonfinite")

    return Eq45CylindricalState(
        q=q,
        X=coordinates.X,
        eta=coordinates.eta,
        psi=psi,
        u_r=u_r,
        u_theta=u_theta,
        u_z=u_z,
    )
