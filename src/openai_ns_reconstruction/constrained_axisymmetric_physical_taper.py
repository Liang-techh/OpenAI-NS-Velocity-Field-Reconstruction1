"""Divergence-preserving physical-space taper for axisymmetric candidates.

The Eq. (4.5) similarity-profile lane is naturally compact in similarity
coordinates only. This module provides a small autonomous exterior connection
that can impose the preregistered physical cylinder ``r < R, |z| < H`` without
multiplying Cartesian velocity by a scalar (which would generally create
divergence).

The poloidal part is tapered at the streamfunction level. If a base
axisymmetric field is represented by cylindrical values ``(psi, u_r, u_theta,
u_z)`` with

    u_r = -(1/r) d_z psi,
    u_z =  (1/r) d_r psi,

then for ``psi_tilde = T(r^2,z) * psi``

    u_r_tilde = T*u_r - (psi/r) * d_z T,
    u_z_tilde = T*u_z + 2*psi * d_(r^2) T.

The axisymmetric swirl channel may be multiplied by the same scalar taper
without changing divergence. The taper is exactly one on an inner plateau,
exactly zero at/outside the physical support, and is C4 across the plateau and
support interfaces.

This is an autonomous representation device. It is not paper-sourced profile
data, not PDE validation, and not evidence of correspondence to the OpenAI
visualization.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass

import numpy as np


def _cutoff_c4(q: np.ndarray, plateau: float) -> tuple[np.ndarray, np.ndarray]:
    """Return C4 cutoff and derivative with respect to ``q``."""

    q = np.asarray(q, dtype=float)
    value = np.ones_like(q)
    derivative = np.zeros_like(q)

    transition = (q > plateau) & (q < 1.0)
    x = (q[transition] - plateau) / (1.0 - plateau)

    smooth = (
        126.0 * x**5
        - 420.0 * x**6
        + 540.0 * x**7
        - 315.0 * x**8
        + 70.0 * x**9
    )
    dsmooth = (
        630.0 * x**4
        - 2520.0 * x**5
        + 3780.0 * x**6
        - 2520.0 * x**7
        + 630.0 * x**8
    )

    value[transition] = 1.0 - smooth
    derivative[transition] = -dsmooth / (1.0 - plateau)
    value[q >= 1.0] = 0.0
    return value, derivative


@dataclass(frozen=True)
class AxisymmetricPhysicalTaper:
    """Physical compact-support contract for an axisymmetric streamfunction field."""

    radial_support: float = 2.0
    axial_half_height: float = 2.0
    radial_plateau_q: float = 0.64
    axial_plateau_q: float = 0.64
    axis_tolerance: float = 1e-12
    schema: str = "axisymmetric_physical_taper_v1"

    def __post_init__(self) -> None:
        if not np.isfinite(self.radial_support) or self.radial_support <= 0.0:
            raise ValueError("radial_support must be positive and finite")
        if not np.isfinite(self.axial_half_height) or self.axial_half_height <= 0.0:
            raise ValueError("axial_half_height must be positive and finite")
        for name in ("radial_plateau_q", "axial_plateau_q"):
            value = float(getattr(self, name))
            if not np.isfinite(value) or not (0.0 < value < 1.0):
                raise ValueError(f"{name} must lie strictly between 0 and 1")
        if not np.isfinite(self.axis_tolerance) or self.axis_tolerance < 0.0:
            raise ValueError("axis_tolerance must be finite and nonnegative")
        if self.schema != "axisymmetric_physical_taper_v1":
            raise ValueError("unsupported taper schema")

    def to_dict(self) -> dict:
        data = asdict(self)
        data["classification"] = "autonomous_design"
        data["claim_scope"] = "physical_support_representation_only"
        data["pde_validated"] = False
        data["paper_exact"] = False
        data["openai_field_identified"] = False
        return data

    @classmethod
    def from_dict(cls, data: dict) -> "AxisymmetricPhysicalTaper":
        if not isinstance(data, dict):
            raise ValueError("taper payload must be a dict")
        if data.get("schema") != "axisymmetric_physical_taper_v1":
            raise ValueError("unsupported taper schema")
        return cls(
            radial_support=float(data["radial_support"]),
            axial_half_height=float(data["axial_half_height"]),
            radial_plateau_q=float(data["radial_plateau_q"]),
            axial_plateau_q=float(data["axial_plateau_q"]),
            axis_tolerance=float(data.get("axis_tolerance", 1e-12)),
            schema=str(data["schema"]),
        )

    def factors(self, points) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Return ``T``, ``dT/d(r^2)``, and ``dT/dz`` on Cartesian points."""

        points = np.asarray(points, dtype=float)
        if points.ndim == 0 or points.shape[-1] != 3:
            raise ValueError("points must have shape (...,3)")
        if not np.all(np.isfinite(points)):
            raise ValueError("points must be finite")

        x = points[..., 0]
        y = points[..., 1]
        z = points[..., 2]
        s = x * x + y * y

        radial_q = s / (self.radial_support * self.radial_support)
        axial_q = (z * z) / (self.axial_half_height * self.axial_half_height)

        radial, radial_dq = _cutoff_c4(radial_q, self.radial_plateau_q)
        axial, axial_dq = _cutoff_c4(axial_q, self.axial_plateau_q)

        taper = radial * axial
        taper_s = radial_dq / (self.radial_support * self.radial_support) * axial
        taper_z = (
            radial
            * axial_dq
            * (2.0 * z / (self.axial_half_height * self.axial_half_height))
        )
        return taper, taper_s, taper_z

    def apply_cylindrical(
        self,
        points,
        *,
        psi,
        u_r,
        u_theta,
        u_z,
    ) -> np.ndarray:
        """Apply the exterior connection and return Cartesian velocity.

        ``psi`` must be the streamfunction associated with ``u_r`` and ``u_z``.
        This routine does not replace independent divergence/PDE validation.
        """

        points = np.asarray(points, dtype=float)
        if points.ndim == 0 or points.shape[-1] != 3:
            raise ValueError("points must have shape (...,3)")
        if not np.all(np.isfinite(points)):
            raise ValueError("points must be finite")

        lead = points.shape[:-1]
        arrays = {}
        for name, value in {
            "psi": psi,
            "u_r": u_r,
            "u_theta": u_theta,
            "u_z": u_z,
        }.items():
            arr = np.asarray(value, dtype=float)
            try:
                arr = np.broadcast_to(arr, lead)
            except ValueError as exc:
                raise ValueError(f"{name} must broadcast to points.shape[:-1]") from exc
            if not np.all(np.isfinite(arr)):
                raise ValueError(f"{name} must be finite")
            arrays[name] = arr

        x = points[..., 0]
        y = points[..., 1]
        r = np.hypot(x, y)
        axis = r <= self.axis_tolerance
        if np.any(np.abs(arrays["psi"][axis]) > self.axis_tolerance):
            raise ValueError("axis-regular streamfunction must vanish on the axis")

        taper, taper_s, taper_z = self.factors(points)
        psi_over_r = np.divide(
            arrays["psi"],
            r,
            out=np.zeros_like(r),
            where=~axis,
        )

        radial = taper * arrays["u_r"] - psi_over_r * taper_z
        axial = taper * arrays["u_z"] + 2.0 * arrays["psi"] * taper_s
        swirl = taper * arrays["u_theta"]

        cos_theta = np.divide(x, r, out=np.ones_like(r), where=~axis)
        sin_theta = np.divide(y, r, out=np.zeros_like(r), where=~axis)
        u = radial * cos_theta - swirl * sin_theta
        v = radial * sin_theta + swirl * cos_theta
        u = np.where(axis, 0.0, u)
        v = np.where(axis, 0.0, v)

        out = np.stack((u, v, axial), axis=-1)
        if not np.all(np.isfinite(out)):
            raise RuntimeError("tapered velocity became nonfinite")
        return out
