"""Exact Part I inward heat collar with ODE-integrated stress.

This module implements the finite collar in Lei--Ren, arXiv:2609.35406v1,
Section 5.  The heat exterior is the repository's checked ``HeatExterior``
with ``delta = 2*h``.  On

    R_b * exp(-ell) <= R <= R_b,

the angular profile is multiplied by the flat factor
``f(y) = 1 - epsilon * exp(-4/y**2)``, where ``y = log(R_b/R)``.  The
inertial stress is obtained by integrating the Section 3 source equations
backwards from the exact heat boundary values.  It is not fitted to sampled
stress data.

The implementation is a finite numerical collar receipt.  It does not claim
global moment matching, a complete outer construction, the five moments,
the PDE, or the full admissible background.  Near ``y = 0`` all stress
quantities are retained in the flat-factor normalization used in (5.23),
because their unnormalized values underflow in binary64.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
import operator
from functools import lru_cache
from pathlib import Path
import sys
from typing import Any, Callable

import numpy as np
from numpy.polynomial.laguerre import laggauss
from numpy.polynomial.legendre import leggauss


_SRC = Path(__file__).resolve().parents[2] / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

from openai_ns_reconstruction.heat_exterior import HeatExterior, heat_factor  # noqa: E402


SOURCE = "https://arxiv.org/html/2609.35406v1"
SOURCE_VERSION = "2609.35406v1"
SOURCE_SECTIONS = "Sections 3.2, 3.4 and 5, equations (5.12)-(5.24)"

DEFAULT_H = 0.001
DEFAULT_C_INF = 0.1
DEFAULT_R_B = 2048.0
DEFAULT_ELL = 0.5
DEFAULT_LAGUERRE_ORDER = 64
TAYLOR_ORDER = 12
TAYLOR_XI_MAX = 0.003


def _finite(value: Any, name: str) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be finite") from exc
    if not math.isfinite(result):
        raise ValueError(f"{name} must be finite")
    return result


def _integer(value: Any, name: str, *, allowed: tuple[int, ...] | None = None) -> int:
    if isinstance(value, (bool, np.bool_)):
        raise TypeError(f"{name} must be an integer")
    try:
        result = operator.index(value)
    except TypeError as exc:
        raise TypeError(f"{name} must be an integer") from exc
    result = int(result)
    if allowed is not None and result not in allowed:
        raise ValueError(f"{name} must be one of {allowed}")
    return result


def _rising(value: float, count: int) -> float:
    result = 1.0
    for index in range(count):
        result *= value + index
    return result


@lru_cache(maxsize=8)
def _laguerre_rule(order: int) -> tuple[np.ndarray, np.ndarray]:
    nodes, weights = laggauss(order)
    nodes = np.asarray(nodes, dtype=float)
    weights = np.asarray(weights, dtype=float)
    nodes.setflags(write=False)
    weights.setflags(write=False)
    return nodes, weights


@lru_cache(maxsize=8)
def _unit_legendre_rule(order: int) -> tuple[np.ndarray, np.ndarray]:
    nodes, weights = leggauss(order)
    nodes = (np.asarray(nodes, dtype=float) + 1.0) / 2.0
    weights = np.asarray(weights, dtype=float) / 2.0
    nodes.setflags(write=False)
    weights.setflags(write=False)
    return nodes, weights


@dataclass(frozen=True)
class HeatJet:
    """Finite heat-factor evaluation and its explicit Taylor error bounds."""

    value: float
    derivative: float
    value_error_bound: float | None
    derivative_error_bound: float | None
    method: str


def _heat_jet(xi: float, h: float, *, order: int = TAYLOR_ORDER) -> HeatJet:
    """Evaluate H and H' using a finite jet when xi is small.

    The coefficient identity is
    H^(m)(0) = (-1)^m (h)_m (1+h)_m.  Complete monotonicity gives the
    displayed next-term bounds for the finite Taylor approximations; no
    convergent infinite-series claim is made.
    """

    xi = _finite(xi, "xi")
    h = _finite(h, "h")
    order = _integer(order, "order")
    if xi < 0.0:
        raise ValueError("xi must be nonnegative")
    if order < 1:
        raise ValueError("order must be positive")
    if xi > TAYLOR_XI_MAX:
        return HeatJet(
            value=float(heat_factor(xi, h)),
            derivative=float(heat_factor(xi, h, derivative=1)),
            value_error_bound=None,
            derivative_error_bound=None,
            method="checked_heat_factor",
        )

    value = 1.0
    derivative = 0.0
    power = 1.0
    coefficient = 1.0
    for degree in range(1, order + 1):
        coefficient *= (
            -(h + degree - 1.0) * (1.0 + h + degree - 1.0) / degree
        )
        power *= xi
        value += coefficient * power
        derivative += degree * coefficient * xi ** (degree - 1)

    next_derivative_at_zero = _rising(h, order + 1) * _rising(1.0 + h, order + 1)
    value_error = (
        next_derivative_at_zero
        * xi ** (order + 1)
        / math.factorial(order + 1)
    )
    derivative_error = (
        next_derivative_at_zero * xi**order / math.factorial(order)
    )
    return HeatJet(
        value=float(value),
        derivative=float(derivative),
        value_error_bound=float(value_error),
        derivative_error_bound=float(derivative_error),
        method=f"finite_taylor_{order}",
    )


class HeatCollar:
    """Finite exact-heat inward collar with Section 3 stress integration."""

    def __init__(
        self,
        *,
        h: float = DEFAULT_H,
        c_inf: float = DEFAULT_C_INF,
        R_b: float = DEFAULT_R_B,
        ell: float = DEFAULT_ELL,
        epsilon: float | None = None,
        laguerre_order: int = DEFAULT_LAGUERRE_ORDER,
        heat: HeatExterior | None = None,
    ) -> None:
        h = _finite(h, "h")
        c_inf = _finite(c_inf, "c_inf")
        R_b = _finite(R_b, "R_b")
        ell = _finite(ell, "ell")
        if not 0.0 < h <= 0.005:
            raise ValueError("h must satisfy 0 < h <= .005 for the Part I collar")
        if c_inf <= 0.0:
            raise ValueError("c_inf must be positive")
        if R_b <= 0.0:
            raise ValueError("R_b must be positive")
        if not 0.0 < ell <= 1.0:
            raise ValueError("ell must satisfy 0 < ell <= 1")
        laguerre_order = _integer(
            laguerre_order,
            "laguerre_order",
            allowed=(32, 64),
        )
        delta = 2.0 * h
        if epsilon is None:
            epsilon = delta / 100.0
        epsilon = _finite(epsilon, "epsilon")
        if not 0.0 < epsilon <= 0.5:
            raise ValueError("epsilon must satisfy 0 < epsilon <= .5")
        if heat is None:
            heat = HeatExterior(h=h, c_inf=c_inf)
        if not isinstance(heat, HeatExterior):
            raise TypeError("heat must be a HeatExterior instance")
        if abs(float(heat.h) - h) > 1e-15 or abs(float(heat.c_inf) - c_inf) > 1e-15:
            raise ValueError("heat parameters must agree with h and c_inf")

        self.h = h
        self.delta = delta
        self.c_inf = c_inf
        self.R_b = R_b
        self.ell = ell
        self.epsilon = epsilon
        self.heat = heat
        self.laguerre_order = laguerre_order
        self.taylor_order = TAYLOR_ORDER
        self.taylor_xi_max = TAYLOR_XI_MAX
        self._pressure_legendre_order = 64

    @property
    def collar_inner_radius(self) -> float:
        return self.R_b * math.exp(-self.ell)

    def _coordinates(self, R: Any, Z: Any) -> tuple[float, float, float]:
        radius = _finite(R, "R")
        z = _finite(Z, "Z")
        if radius <= 0.0:
            raise ValueError("collar requires R>0")
        if abs(z) > 1.0:
            raise ValueError("collar requires |Z|<=1")
        y = math.log(self.R_b / radius)
        if y < -1e-13 or y > self.ell + 1e-13:
            raise ValueError(
                "collar coordinate must satisfy R_b*exp(-ell) <= R <= R_b"
            )
        y = min(max(y, 0.0), self.ell)
        return radius, z, y

    def y(self, R: Any) -> float:
        radius = _finite(R, "R")
        if radius <= 0.0:
            raise ValueError("R must be positive")
        result = math.log(self.R_b / radius)
        if result < -1e-13 or result > self.ell + 1e-13:
            raise ValueError("R is outside the collar")
        return float(min(max(result, 0.0), self.ell))

    def R_at_y(self, y: Any) -> float:
        value = _finite(y, "y")
        if not 0.0 <= value <= self.ell:
            raise ValueError("y must lie in [0, ell]")
        return float(self.R_b * math.exp(-value))

    @staticmethod
    def _flat(y: float) -> float:
        if y <= 0.0:
            return 0.0
        exponent = -4.0 / (y * y)
        if exponent < math.log(np.finfo(float).tiny):
            return 0.0
        return float(math.exp(exponent))

    def flat_factor(self, y: Any) -> float:
        value = _finite(y, "y")
        if not 0.0 <= value <= self.ell:
            raise ValueError("y must lie in [0, ell]")
        return float(self._flat(value))

    def f(self, y: Any) -> float:
        value = _finite(y, "y")
        if not 0.0 <= value <= self.ell:
            raise ValueError("y must lie in [0, ell]")
        return float(1.0 - self.epsilon * self._flat(value))

    def f_y(self, y: Any) -> float:
        value = _finite(y, "y")
        if not 0.0 <= value <= self.ell:
            raise ValueError("y must lie in [0, ell]")
        if value == 0.0 or self._flat(value) == 0.0:
            return 0.0
        return float(-8.0 * self.epsilon * self._flat(value) / value**3)

    def heat_jet(self, xi: Any) -> dict[str, float | str | None]:
        jet = _heat_jet(_finite(xi, "xi"), self.h, order=self.taylor_order)
        return {
            "value": jet.value,
            "derivative": jet.derivative,
            "value_error_bound": jet.value_error_bound,
            "derivative_error_bound": jet.derivative_error_bound,
            "method": jet.method,
        }

    def _heat_values(self, radius: float, z: float) -> tuple[float, float, float]:
        xi = 2.0 * (1.0 - z * z) / radius
        jet = _heat_jet(xi, self.h, order=self.taylor_order)
        prefactor = self.c_inf / math.sqrt(2.0) * radius ** (-1.0 - self.h)
        F = prefactor * jet.value
        F_z = prefactor * jet.derivative * (-4.0 * z / radius)
        F_R = prefactor / radius * (-(1.0 + self.h) * jet.value - xi * jet.derivative)
        return float(F), float(F_z), float(F_R)

    def F_heat(self, R: Any, Z: Any) -> float:
        radius = _finite(R, "R")
        z = _finite(Z, "Z")
        if radius <= 0.0 or abs(z) > 1.0:
            raise ValueError("require R>0 and |Z|<=1")
        return self._heat_values(radius, z)[0]

    def F_heat_Z(self, R: Any, Z: Any) -> float:
        radius = _finite(R, "R")
        z = _finite(Z, "Z")
        if radius <= 0.0 or abs(z) > 1.0:
            raise ValueError("require R>0 and |Z|<=1")
        return self._heat_values(radius, z)[1]

    def F_heat_R(self, R: Any, Z: Any) -> float:
        radius = _finite(R, "R")
        z = _finite(Z, "Z")
        if radius <= 0.0 or abs(z) > 1.0:
            raise ValueError("require R>0 and |Z|<=1")
        return self._heat_values(radius, z)[2]

    def _heat_pressure_Z(self, radius: float, z: float) -> float:
        nodes, weights = _unit_legendre_rule(self._pressure_legendre_order)
        exponent = 2.0 * (0.5 + self.h) - 1.0
        xi = 2.0 * (1.0 - z * z) / radius
        H_values = np.empty_like(nodes)
        Hp_values = np.empty_like(nodes)
        for index, point in enumerate(nodes):
            jet = _heat_jet(xi * float(point), self.h, order=self.taylor_order)
            H_values[index] = jet.value
            Hp_values[index] = jet.derivative
        integral = float(
            np.dot(weights * nodes**exponent * nodes, H_values * Hp_values)
        )
        return float(self.c_inf**2 * radius ** (-(1.0 + 2.0 * self.h)) * (4.0 * z / radius) * integral)

    def P_heat(self, R: Any, Z: Any) -> float:
        radius = _finite(R, "R")
        z = _finite(Z, "Z")
        if radius <= 0.0 or abs(z) > 1.0:
            raise ValueError("require R>0 and |Z|<=1")
        return float(
            self.heat.pressure_from_tau(math.sqrt(2.0 * radius), 1.0 - z * z)
        )

    def P_heat_Z(self, R: Any, Z: Any) -> float:
        radius = _finite(R, "R")
        z = _finite(Z, "Z")
        if radius <= 0.0 or abs(z) > 1.0:
            raise ValueError("require R>0 and |Z|<=1")
        return self._heat_pressure_Z(radius, z)

    def _scaled_flat_integral(
        self,
        y: float,
        kernel: Callable[[float], float],
        *,
        radial_exponent: float = 1.0,
    ) -> float:
        """Return integral_0^y exp(-4/s^2) k(s) ds / exp(-4/y^2)."""

        if y <= 0.0:
            return 0.0
        nodes, weights = _laguerre_rule(self.laguerre_order)
        total = 0.0
        y2_over_4 = y * y / 4.0
        for node, weight in zip(nodes, weights):
            s = y / math.sqrt(1.0 + float(node) * y2_over_4)
            total += float(weight) * math.exp(radial_exponent * (y - s)) * float(
                kernel(s)
            ) * s**3 / 8.0
        return float(total)

    def _scaled_flat_integral_zstress(
        self,
        y: float,
        kernel: Callable[[float], float],
    ) -> float:
        """Return I_z/(flat*y^3) from the scaled Section 3 source."""

        if y <= 0.0:
            return 0.0
        nodes, weights = _laguerre_rule(self.laguerre_order)
        total = 0.0
        y2_over_4 = y * y / 4.0
        for node, weight in zip(nodes, weights):
            s = y / math.sqrt(1.0 + float(node) * y2_over_4)
            total += float(weight) * math.exp(0.5 * (y - s)) * float(
                kernel(s)
            ) * (s / y) ** 3 / 8.0
        return float(-total)

    def _heat_U_values(self, radius: float, z: float) -> tuple[float, float, float]:
        F, F_z, F_R = self._heat_values(radius, z)
        scale = math.sqrt(2.0 * radius)
        U = scale * F
        U_z = scale * F_z
        U_R = scale * (F_R + F / (2.0 * radius))
        return float(U), float(U_z), float(U_R)

    def _Ntheta_heat(self, radius: float, z: float) -> float:
        U, U_z, U_R = self._heat_U_values(radius, z)
        L = 1.0 - self.delta * z * z
        bracket = (
            (1.0 + self.delta) * U / 2.0
            + (1.0 - self.delta) * z * U_z / 2.0
            + radius * U_R
        )
        return float(-math.sqrt(radius / 2.0) * bracket / L)

    def _delta_Ntheta_over_flat(self, y: float, z: float) -> float:
        if y <= 0.0:
            raise ValueError("scaled angular source requires y>0")
        radius = self.R_b * math.exp(-y)
        U_heat, _U_z, _U_R = self._heat_U_values(radius, z)
        L = 1.0 - self.delta * z * z
        N_heat = self._Ntheta_heat(radius, z)
        # (f-1)/exp(-4/y^2) = -epsilon and f_y/exp(-4/y^2) = -8 epsilon/y^3.
        return float(
            -self.epsilon * N_heat
            - math.sqrt(radius / 2.0) * U_heat * 8.0 * self.epsilon / (L * y**3)
        )

    def _delta_P_over_flat(self, y: float, z: float) -> float:
        if y <= 0.0:
            return 0.0

        def kernel(s: float) -> float:
            radius = self.R_b * math.exp(-s)
            flat = self._flat(s)
            F_heat = self._heat_values(radius, z)[0]
            return self.epsilon * (2.0 - self.epsilon * flat) * radius * F_heat**2

        # The radial Jacobian ``R_b exp(-s)`` is already part of ``kernel``;
        # do not apply the angular-stress exp(y-s) factor a second time.
        return self._scaled_flat_integral(y, kernel, radial_exponent=0.0)

    def _delta_P_Z_over_flat(self, y: float, z: float) -> float:
        if y <= 0.0:
            return 0.0

        def kernel(s: float) -> float:
            radius = self.R_b * math.exp(-s)
            flat = self._flat(s)
            F_heat, F_heat_Z, _F_heat_R = self._heat_values(radius, z)
            return (
                self.epsilon
                * (2.0 - self.epsilon * flat)
                * radius
                * 2.0
                * F_heat
                * F_heat_Z
            )

        # As above, the radial Jacobian is included in ``kernel``.
        return self._scaled_flat_integral(y, kernel, radial_exponent=0.0)

    def delta_P_R_over_flat(self, R: Any, Z: Any) -> float:
        radius, z, y = self._coordinates(R, Z)
        flat = self._flat(y)
        F_heat = self._heat_values(radius, z)[0]
        return float(-self.epsilon * (2.0 - self.epsilon * flat) * F_heat**2)

    def delta_P_over_flat(self, R: Any, Z: Any) -> float:
        radius, z, y = self._coordinates(R, Z)
        del radius
        return self._delta_P_over_flat(y, z)

    def delta_P_Z_over_flat(self, R: Any, Z: Any) -> float:
        radius, z, y = self._coordinates(R, Z)
        del radius
        return self._delta_P_Z_over_flat(y, z)

    def delta_P(self, R: Any, Z: Any) -> float:
        radius, z, y = self._coordinates(R, Z)
        return float(self._flat(y) * self._delta_P_over_flat(y, z))

    def delta_P_Z(self, R: Any, Z: Any) -> float:
        radius, z, y = self._coordinates(R, Z)
        return float(self._flat(y) * self._delta_P_Z_over_flat(y, z))

    def F(self, R: Any, Z: Any) -> float:
        radius, z, y = self._coordinates(R, Z)
        return float((1.0 - self.epsilon * self._flat(y)) * self._heat_values(radius, z)[0])

    def F_Z(self, R: Any, Z: Any) -> float:
        radius, z, y = self._coordinates(R, Z)
        return float((1.0 - self.epsilon * self._flat(y)) * self._heat_values(radius, z)[1])

    def F_R(self, R: Any, Z: Any) -> float:
        radius, z, y = self._coordinates(R, Z)
        flat = self._flat(y)
        f_y = 0.0 if y == 0.0 else -8.0 * self.epsilon * flat / y**3
        f_R = -f_y / radius
        F_heat, _F_heat_Z, F_heat_R = self._heat_values(radius, z)
        return float((1.0 - self.epsilon * flat) * F_heat_R + f_R * F_heat)

    def P(self, R: Any, Z: Any) -> float:
        radius, z, y = self._coordinates(R, Z)
        return float(self.P_heat(radius, z) + self._flat(y) * self._delta_P_over_flat(y, z))

    def P_Z(self, R: Any, Z: Any) -> float:
        radius, z, y = self._coordinates(R, Z)
        return float(self.P_heat_Z(radius, z) + self._flat(y) * self._delta_P_Z_over_flat(y, z))

    def P_R(self, R: Any, Z: Any) -> float:
        radius, z, _y = self._coordinates(R, Z)
        del radius
        return float(self.F(R, Z) ** 2)

    def N_theta(self, R: Any, Z: Any) -> float:
        radius, z, y = self._coordinates(R, Z)
        N_heat = self._Ntheta_heat(radius, z)
        if y == 0.0:
            return float(N_heat)
        U_heat = self._heat_U_values(radius, z)[0]
        L = 1.0 - self.delta * z * z
        return float(
            (1.0 - self.epsilon * self._flat(y)) * N_heat
            + math.sqrt(radius / 2.0) * U_heat * self.f_y(y) / L
        )

    def N_z(self, R: Any, Z: Any) -> float:
        radius, z, y = self._coordinates(R, Z)
        if y == 0.0:
            return 0.0
        L = 1.0 - self.delta * z * z
        scaled = self._delta_Nz_over_flat(y, z)
        return float(self._flat(y) * scaled)

    def _delta_Nz_over_flat(self, y: float, z: float) -> float:
        if y <= 0.0:
            return 0.0
        radius = self.R_b * math.exp(-y)
        L = 1.0 - self.delta * z * z
        delta_P_scaled = self._delta_P_over_flat(y, z)
        delta_P_Z_scaled = self._delta_P_Z_over_flat(y, z)
        delta_P_R_scaled = self.delta_P_R_over_flat(radius, z)
        return float(
            math.sqrt(radius / 2.0)
            / L
            * (
                2.0 * (1.0 + self.delta) * z * delta_P_scaled
                - (1.0 - z * z) * delta_P_Z_scaled
                + 2.0 * z * radius * delta_P_R_scaled
            )
        )

    def S_theta(self, R: Any, Z: Any) -> float:
        radius, z, _y = self._coordinates(R, Z)
        return float(2.0 * radius * self.F_R(radius, z))

    def S_z(self, R: Any, Z: Any) -> float:
        self._coordinates(R, Z)
        return 0.0

    def _delta_Itheta_over_flat(self, y: float, z: float) -> float:
        if y <= 0.0:
            return 0.0
        return -self._scaled_flat_integral(
            y,
            lambda s: self._delta_Ntheta_over_flat(s, z),
            radial_exponent=1.0,
        )

    def _delta_Iz_over_flat_y3(self, y: float, z: float) -> float:
        if y <= 0.0:
            return 0.0
        return self._scaled_flat_integral_zstress(
            y,
            lambda s: self._delta_Nz_over_flat(s, z),
        )

    def I_theta_heat(self, R: Any, Z: Any) -> float:
        return float(-2.0 * self._coordinates(R, Z)[0] * self._heat_values(
            self._coordinates(R, Z)[0], self._coordinates(R, Z)[1]
        )[2])

    def I_theta(self, R: Any, Z: Any) -> float:
        radius, z, y = self._coordinates(R, Z)
        return float(
            self.I_theta_heat(radius, z)
            + self._flat(y) * self._delta_Itheta_over_flat(y, z)
        )

    def I_z(self, R: Any, Z: Any) -> float:
        radius, z, y = self._coordinates(R, Z)
        if y == 0.0:
            return 0.0
        return float(
            self._flat(y)
            * y**3
            * self._delta_Iz_over_flat_y3(y, z)
        )

    def T_theta(self, R: Any, Z: Any) -> float:
        radius, z, _y = self._coordinates(R, Z)
        return float(self.I_theta(radius, z) + self.S_theta(radius, z))

    def T_z(self, R: Any, Z: Any) -> float:
        radius, z, _y = self._coordinates(R, Z)
        return float(self.I_z(radius, z))

    def b_theta(self, R: Any, Z: Any) -> float:
        radius, z, y = self._coordinates(R, Z)
        F_heat = self._heat_values(radius, z)[0]
        if y == 0.0:
            return float(16.0 * self.epsilon * F_heat)
        S_heat = 2.0 * radius * self._heat_values(radius, z)[2]
        return float(
            y**3 * self._delta_Itheta_over_flat(y, z)
            + self.epsilon * (16.0 * F_heat - y**3 * S_heat)
        )

    def b_z(self, R: Any, Z: Any) -> float:
        radius, z, y = self._coordinates(R, Z)
        F_heat = self._heat_values(radius, z)[0]
        L = 1.0 - self.delta * z * z
        if y == 0.0:
            return float(
                self.epsilon * z * radius**1.5 * F_heat**2 / (2.0 * math.sqrt(2.0) * L)
            )
        return float(self._delta_Iz_over_flat_y3(y, z))

    def normalized(self, R: Any, Z: Any) -> dict[str, float | bool]:
        radius, z, y = self._coordinates(R, Z)
        F = self.F(radius, z)
        S_theta = self.S_theta(radius, z)
        kappa = float(-S_theta / F)
        b_theta = self.b_theta(radius, z)
        b_z = self.b_z(radius, z)
        ratio = 0.0 if y == 0.0 else float(y**6 * b_z / b_theta)
        margin = float(2.0 - (kappa - 2.0) * ratio**2)
        return {
            "R": radius,
            "Z": z,
            "y": y,
            "F": F,
            "S_theta": S_theta,
            "kappa": kappa,
            "btheta": b_theta,
            "bz": b_z,
            "Tz_over_Ttheta": ratio,
            "directional_margin": margin,
            "Ttheta_positive_normalized": bool(b_theta > 0.0),
            "Stheta_negative": bool(S_theta < 0.0),
            "kappa_margin": float(kappa - (2.0 + self.delta / 2.0)),
        }

    def cone_data(self, R: Any, Z: Any) -> dict[str, float | bool]:
        return self.normalized(R, Z)

    def inertial(self, R: Any, Z: Any) -> tuple[float, float]:
        return self.I_theta(R, Z), self.I_z(R, Z)

    def shear(self, R: Any, Z: Any) -> tuple[float, float]:
        return self.S_theta(R, Z), self.S_z(R, Z)

    def stress(self, R: Any, Z: Any) -> tuple[float, float]:
        return self.T_theta(R, Z), self.T_z(R, Z)

    def metadata(self) -> dict[str, Any]:
        return {
            "source": SOURCE,
            "source_version": SOURCE_VERSION,
            "source_sections": SOURCE_SECTIONS,
            "status": "finite exact-heat inward collar with ODE-integrated stress",
            "h": self.h,
            "delta": self.delta,
            "c_inf": self.c_inf,
            "R_b": self.R_b,
            "ell": self.ell,
            "collar_inner_radius": self.collar_inner_radius,
            "epsilon": self.epsilon,
            "laguerre_order": self.laguerre_order,
            "taylor_order": self.taylor_order,
            "taylor_xi_max": self.taylor_xi_max,
            "heat_mapping": "delta=2h; F_heat=HeatExterior.profile_E/sqrt(2R)",
            "stress_source": "Section 3.2 Ntheta,Nz and radial ODE (3.13)",
            "global_moment_matching": False,
            "full_background": False,
            "pde_certified": False,
        }


__all__ = [
    "DEFAULT_C_INF",
    "DEFAULT_ELL",
    "DEFAULT_H",
    "DEFAULT_LAGUERRE_ORDER",
    "DEFAULT_R_B",
    "HeatCollar",
    "HeatJet",
    "SOURCE",
    "SOURCE_SECTIONS",
    "SOURCE_VERSION",
]
