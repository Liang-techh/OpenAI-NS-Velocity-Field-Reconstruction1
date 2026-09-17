"""Preregistered two-parameter force; independent of any candidate residual."""
from dataclasses import dataclass

import numpy as np


def compact_bump(s):
    """B(s) and B'(s), with the smooth zero extension at s=1."""
    s = np.asarray(s, dtype=float)
    if not np.all(np.isfinite(s)):
        raise ValueError("bump arguments must be finite")
    # Clamp only the inactive branch to avoid division by zero/overflow.
    d = np.where(s < 1, 1 - s, 1.0)
    value = np.where(s < 1, np.exp(1 - 1 / d), 0.0)
    derivative = -value / d**2
    return value, derivative


@dataclass(frozen=True)
class RestrictedForce:
    """f=g curl[b(-a*z*y, a*z*x, -c*r^2/2)], 0 <= a,c <= 10.

    Points have shape (...,3); time broadcasts against their leading shape.
    Length/time support and bounds are fixed by constraints.json v1.
    """

    a: float = 1.0
    c: float = 1.0

    def __post_init__(self):
        if not all(np.isfinite(v) and 0 <= v <= 10 for v in (self.a, self.c)):
            raise ValueError("force coefficients must be finite and in [0,10]")

    def _terms(self, points, time):
        points = np.asarray(points, dtype=float)
        time = np.asarray(time, dtype=float)
        if points.ndim < 1 or points.shape[-1] != 3:
            raise ValueError("points must have shape (...,3)")
        if not np.all(np.isfinite(points)) or not np.all(np.isfinite(time)):
            raise ValueError("points and time must be finite")
        x, y, z, t = np.broadcast_arrays(*np.moveaxis(points, -1, 0), time)
        r2 = x*x + y*y
        radial, radial_d = compact_bump(r2 / 4)
        axial, axial_d = compact_bump(z*z / 4)
        temporal, _ = compact_bump((2*t-1)**2)
        b = radial * axial
        # b_r/r and b_z avoid any singular division by r on the axis.
        br_over_r = radial_d * axial / 2
        bz = radial * axial_d * z / 2
        return x, y, z, r2, b, br_over_r, bz, temporal

    def potential(self, points, time):
        x, y, z, r2, b, _, _, g = self._terms(points, time)
        return g[..., None] * np.stack(
            (-self.a*z*y*b, self.a*z*x*b, -self.c*r2*b/2), axis=-1
        )

    def __call__(self, points, time):
        x, y, z, r2, b, br, bz, g = self._terms(points, time)
        radial_factor = -self.a * (b + z*bz)
        swirl_factor = self.c * (b + r2*br/2)
        axial = self.a*z*(2*b + r2*br)
        return g[..., None] * np.stack(
            (x*radial_factor-y*swirl_factor,
             y*radial_factor+x*swirl_factor, axial), axis=-1
        )
