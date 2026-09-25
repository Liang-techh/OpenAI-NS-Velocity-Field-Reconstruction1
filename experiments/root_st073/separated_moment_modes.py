"""Axisymmetric divergence-free bridge bumps separated in radius and scale.

Exploratory physical moment directions inspired by Appendix A; the smooth
time interpolation and finite mode set are not the paper's recursion.
"""

import numpy as np

from width_field import coordinates


RADIAL_WINDOWS = ((0.12, 0.38), (0.62, 0.88))
RADIAL_WINDOWS_THREE = ((0.12, 0.38), (0.40, 0.60), (0.62, 0.88))


def flat_transition(s):
    if s <= 0.0:
        return 0.0
    if s >= 1.0:
        return 1.0
    left = np.exp(-1.0 / s)
    right = np.exp(-1.0 / (1.0 - s))
    return float(left / (left + right))


def scale_partition(k, knots):
    """C-infinity partition with flat jets at each adjacent scale knot."""
    weights = np.zeros(len(knots))
    if k <= knots[0]:
        weights[0] = 1.0
    elif k >= knots[-1]:
        weights[-1] = 1.0
    else:
        left_index = int(np.searchsorted(knots, k, side="right") - 1)
        s = (k - knots[left_index]) / (knots[left_index + 1] - knots[left_index])
        right = flat_transition(s)
        weights[left_index] = 1.0 - right
        weights[left_index + 1] = right
    return weights


def scale_weight(k):
    """Backward-compatible right-knot weight for the two-knot model."""
    return float(scale_partition(k, (11.0, 19.0))[1])


def flat_bump(y, lo, hi):
    value = np.zeros_like(y)
    derivative = np.zeros_like(y)
    inside = (y > lo) & (y < hi)
    s = (y[inside] - lo) / (hi - lo)
    base = np.exp(4.0 - 1.0 / (s * (1.0 - s)))
    value[inside] = base
    derivative[inside] = (base * (1.0 - 2.0 * s)
                          / (s**2 * (1.0 - s)**2 * (hi - lo)))
    return value, derivative


class SeparatedMomentModes:
    """Scale knots × swirl/poloidal × radial windows × even/odd eta."""

    def __init__(self, base, amplitudes, windows=RADIAL_WINDOWS,
                 knots=(11.0, 19.0)):
        self.base = base
        self.inner = base.inner
        self.nu = base.nu
        self.join_X = base.join_X
        self.ratio = base.ratio
        self.windows = tuple(windows)
        self.knots = tuple(float(k) for k in knots)
        if len(self.knots) < 2 or any(b <= a for a, b in zip(self.knots, self.knots[1:])):
            raise ValueError("Need increasing scale knots")
        self.a = np.asarray(amplitudes, float).reshape(
            len(self.knots), 2, len(self.windows), 2)

    def fields(self, points, tau):
        points = np.asarray(points, float)
        times = np.broadcast_to(np.asarray(tau, float), (len(points),))
        if not np.all(times == times[0]):
            raise ValueError("Separated modes require one remaining time")
        u, p = self.base.fields(points, times)
        radius = np.hypot(points[:, 0], points[:, 1])
        safe = np.where(radius > 0, radius, 1.0)
        sn = np.sqrt(self.nu)
        co = coordinates(radius / sn, points[:, 2] / sn,
                         times, self.inner.h)
        q, eta = np.asarray(co["q"]), np.asarray(co["eta"])
        qz = np.asarray(co["q_z"]) / sn
        etaz = np.asarray(co["eta_z"]) / sn
        ri = np.sqrt(2.0 * self.nu * q * self.join_X)
        width = (self.ratio - 1.0) * ri
        y = (radius - ri) / width
        yz = -(1.0 + (self.ratio - 1.0) * y) * qz / (
            2.0 * q * (self.ratio - 1.0))
        k = -np.log2(2.0 * float(times[0]))
        time_weights = scale_partition(k, self.knots)
        psi_scale = self.nu**1.5 * (2.0 * self.join_X) * q**(1.0 - self.inner.A)
        swirl_scale = sn * q**(-self.inner.A)
        psi_r = np.zeros_like(radius)
        psi_z = np.zeros_like(radius)
        swirl = np.zeros_like(radius)
        for radial_index, (lo, hi) in enumerate(self.windows):
            bump, bump_y = flat_bump(y, lo, hi)
            for parity in (0, 1):
                axial = np.ones_like(eta) if parity == 0 else eta / 0.3
                axial_z = np.zeros_like(eta) if parity == 0 else etaz / 0.3
                for time_index, time_weight in enumerate(time_weights):
                    angular_amp = self.a[time_index, 0, radial_index, parity] * time_weight
                    poloidal_amp = self.a[time_index, 1, radial_index, parity] * time_weight
                    swirl += angular_amp * swirl_scale * bump * axial
                    psi_r += poloidal_amp * psi_scale * bump_y * axial / width
                    psi_z += poloidal_amp * psi_scale * (
                        (1.0 - self.inner.A) * qz / q * bump * axial
                        + bump_y * yz * axial + bump * axial_z)
        ur = -psi_z / safe
        uz = psi_r / safe
        u[:, 0] += ur * points[:, 0] / safe - swirl * points[:, 1] / safe
        u[:, 1] += ur * points[:, 1] / safe + swirl * points[:, 0] / safe
        u[:, 2] += uz
        return u, p
