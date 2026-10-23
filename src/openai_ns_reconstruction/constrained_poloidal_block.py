"""Small axis-regular divergence-free poloidal correction block for CR003.

The block is deliberately candidate-agnostic. It supplies a compact velocity
correction generated as curl((-y q, x q, 0)) with q=q(s,z,t), s=x**2+y**2.
This removes cylindrical 1/r factors and gives an exact divergence identity.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np


@dataclass(frozen=True)
class PoloidalKinematics:
    velocity: np.ndarray
    spatial_jacobian: np.ndarray
    time_derivative: np.ndarray
    divergence: np.ndarray


@dataclass(frozen=True)
class CompactPoloidalBlock:
    """Four-parameter compact poloidal correction on a fixed time window.

    The scalar potential is

        q(s,z,t) = R(s) Z(z) [c0 + c1 z + c2 tau + c3 z tau],

    where ``tau`` maps ``time_interval`` to [-1, 1],
    ``R(s)=(1-s/radius**2)^5_+`` and
    ``Z(z)=(1-(z/half_height)**2)^5_+``.

    The induced Cartesian correction is

        (-x q_z, -y q_z, 2 q + 2 s q_s),

    which is ``curl((-y q, x q, 0))`` and is therefore divergence-free for
    the smooth interior formula. The fifth-power cutoff makes q C^4 across
    the declared support boundary, enough for the current velocity/Jacobian
    interface and future second-velocity-derivative work.
    """

    coefficients: tuple[float, float, float, float]
    coefficient_bound: float = 1.0
    support_radius: float = 2.0
    support_half_height: float = 2.0
    time_interval: tuple[float, float] = (0.25, 0.75)

    def __post_init__(self) -> None:
        coeff = np.asarray(self.coefficients, dtype=float)
        if coeff.shape != (4,) or not np.all(np.isfinite(coeff)):
            raise ValueError("coefficients must contain four finite values")
        if not np.isfinite(self.coefficient_bound) or self.coefficient_bound <= 0.0:
            raise ValueError("coefficient_bound must be finite and positive")
        if np.any(np.abs(coeff) > self.coefficient_bound):
            raise ValueError("coefficient exceeds coefficient_bound")
        if not np.isfinite(self.support_radius) or self.support_radius <= 0.0:
            raise ValueError("support_radius must be finite and positive")
        if not np.isfinite(self.support_half_height) or self.support_half_height <= 0.0:
            raise ValueError("support_half_height must be finite and positive")
        t0, t1 = self.time_interval
        if not np.isfinite(t0) or not np.isfinite(t1) or not t0 < t1:
            raise ValueError("time_interval must be finite and increasing")

    @property
    def parameter_count(self) -> int:
        return 4

    @property
    def is_nontrivial(self) -> bool:
        return bool(np.linalg.norm(np.asarray(self.coefficients, dtype=float)) > 0.0)

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema": "compact_poloidal_block_v1",
            "coefficients": [float(value) for value in self.coefficients],
            "coefficient_bound": float(self.coefficient_bound),
            "support_radius": float(self.support_radius),
            "support_half_height": float(self.support_half_height),
            "time_interval": [float(value) for value in self.time_interval],
        }

    @classmethod
    def from_dict(cls, payload: dict[str, Any]) -> "CompactPoloidalBlock":
        if payload.get("schema") != "compact_poloidal_block_v1":
            raise ValueError("unsupported compact poloidal block schema")
        return cls(
            coefficients=tuple(payload["coefficients"]),
            coefficient_bound=float(payload["coefficient_bound"]),
            support_radius=float(payload["support_radius"]),
            support_half_height=float(payload["support_half_height"]),
            time_interval=tuple(payload["time_interval"]),
        )

    def velocity(self, points: np.ndarray, time: float) -> np.ndarray:
        return self.kinematics(points, time).velocity

    def kinematics(self, points: np.ndarray, time: float) -> PoloidalKinematics:
        pts = np.asarray(points, dtype=float)
        if pts.shape[-1:] != (3,) or not np.all(np.isfinite(pts)):
            raise ValueError("points must have shape (..., 3) and be finite")
        if not np.isfinite(time):
            raise ValueError("time must be finite")
        t0, t1 = self.time_interval
        if not (t0 <= time <= t1):
            raise ValueError("time is outside the declared time_interval")

        original_shape = pts.shape[:-1]
        flat = pts.reshape((-1, 3))
        x, y, z = flat[:, 0], flat[:, 1], flat[:, 2]
        s = x * x + y * y

        q = self._potential_derivatives(s, z, float(time))
        q0, qs, qz, qss, qsz, qzz, qt, qst, qzt = q

        velocity = np.column_stack(
            (
                -x * qz,
                -y * qz,
                2.0 * q0 + 2.0 * s * qs,
            )
        )

        jac = np.empty((flat.shape[0], 3, 3), dtype=float)
        jac[:, 0, 0] = -qz - 2.0 * x * x * qsz
        jac[:, 0, 1] = -2.0 * x * y * qsz
        jac[:, 0, 2] = -x * qzz
        jac[:, 1, 0] = -2.0 * x * y * qsz
        jac[:, 1, 1] = -qz - 2.0 * y * y * qsz
        jac[:, 1, 2] = -y * qzz
        jac[:, 2, 0] = 8.0 * x * qs + 4.0 * x * s * qss
        jac[:, 2, 1] = 8.0 * y * qs + 4.0 * y * s * qss
        jac[:, 2, 2] = 2.0 * qz + 2.0 * s * qsz

        dt = np.column_stack(
            (
                -x * qzt,
                -y * qzt,
                2.0 * qt + 2.0 * s * qst,
            )
        )
        divergence = np.trace(jac, axis1=1, axis2=2)

        return PoloidalKinematics(
            velocity=velocity.reshape(original_shape + (3,)),
            spatial_jacobian=jac.reshape(original_shape + (3, 3)),
            time_derivative=dt.reshape(original_shape + (3,)),
            divergence=divergence.reshape(original_shape),
        )

    def _potential_derivatives(
        self, s: np.ndarray, z: np.ndarray, time: float
    ) -> tuple[np.ndarray, ...]:
        radius_sq = self.support_radius**2
        radial_u = 1.0 - s / radius_sq
        radial_mask = radial_u > 0.0
        ru = np.where(radial_mask, radial_u, 0.0)
        r0 = ru**5
        rs = np.where(radial_mask, -5.0 * ru**4 / radius_sq, 0.0)
        rss = np.where(radial_mask, 20.0 * ru**3 / radius_sq**2, 0.0)

        height_sq = self.support_half_height**2
        axial_u = 1.0 - z * z / height_sq
        axial_mask = axial_u > 0.0
        zu = np.where(axial_mask, axial_u, 0.0)
        z0 = zu**5
        zz = np.where(axial_mask, -10.0 * z * zu**4 / height_sq, 0.0)
        zzz = np.where(
            axial_mask,
            -10.0 * zu**4 / height_sq + 80.0 * z * z * zu**3 / height_sq**2,
            0.0,
        )

        t0, t1 = self.time_interval
        midpoint = 0.5 * (t0 + t1)
        half_span = 0.5 * (t1 - t0)
        tau = (time - midpoint) / half_span
        c0, c1, c2, c3 = (float(value) for value in self.coefficients)
        f = c0 + c1 * z + c2 * tau + c3 * z * tau
        fz = c1 + c3 * tau
        ft = (c2 + c3 * z) / half_span
        fzt = np.full_like(z, c3 / half_span)

        axial_f = z0 * f
        axial_f_z = zz * f + z0 * fz
        axial_f_zz = zzz * f + 2.0 * zz * fz
        axial_f_t = z0 * ft
        axial_f_zt = zz * ft + z0 * fzt

        q0 = r0 * axial_f
        qs = rs * axial_f
        qz = r0 * axial_f_z
        qss = rss * axial_f
        qsz = rs * axial_f_z
        qzz = r0 * axial_f_zz
        qt = r0 * axial_f_t
        qst = rs * axial_f_t
        qzt = r0 * axial_f_zt
        return q0, qs, qz, qss, qsz, qzz, qt, qst, qzt
