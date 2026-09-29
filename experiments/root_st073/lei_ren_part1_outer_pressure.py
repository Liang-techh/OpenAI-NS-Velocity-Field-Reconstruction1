"""Finite core-to-heat pressure binding for the Lei--Ren Part I route.

The source formula is pinned to Z. Lei and X. Ren, arXiv:2609.35406v1,
Section 2.5, equation (2.22): ``M^p(R,Z)=integral_0^R F(rho,Z)^2 d rho``.
The exterior is the existing checked ``HeatExterior`` implementation.  This
module supplies a generic scalar core ``F_core`` and a flat C-infinity radial
blend to that heat profile, then uses one shared ``F`` to define the finite
interior pressure datum and pressure.

The implementation is deliberately local.  It does not construct ``Uz``, the
other four moments, the paper's O(1)--O(8) exterior/collar, moment closure,
admissible stress, PDE solution, global energy bound, or scale recursion.
"""

from __future__ import annotations

import math
import operator
import sys
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Any, Callable

import numpy as np
from numpy.polynomial.legendre import leggauss

_SRC = Path(__file__).resolve().parents[2] / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

try:
    from openai_ns_reconstruction.heat_exterior import HeatExterior
except ModuleNotFoundError:  # pragma: no cover - source package is required
    raise


SOURCE = "https://arxiv.org/html/2609.35406v1"
SOURCE_VERSION = "2609.35406v1"
SOURCE_SECTION = "Section 2.5, equation (2.22)"


@lru_cache(maxsize=32768)
def _cached_heat_F(heat, R, Z):
    # The heat tail is unchanged across nonlinear core-pressure iterations.
    return heat.profile_E(R, Z) / math.sqrt(2.0 * R)


@lru_cache(maxsize=4096)
def _cached_heat_pressure(heat, R, Z):
    return heat.pressure_from_tau(math.sqrt(2.0 * R), 1.0-Z*Z)


def _finite_scalar(value: Any, name: str) -> float:
    """Convert one input to a finite real scalar."""

    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be a finite real scalar") from exc
    if not math.isfinite(result):
        raise ValueError(f"{name} must be finite")
    return result


def _order(n: Any) -> int:
    """Validate a positive Gauss--Legendre order."""

    if isinstance(n, (bool, np.bool_)):
        raise TypeError("n must be a positive integer")
    try:
        result = operator.index(n)
    except TypeError as exc:
        raise TypeError("n must be a positive integer") from exc
    if result <= 0:
        raise ValueError("n must be a positive integer")
    return int(result)


def finite_pressure_moment(
    F: Callable[[float, float], float], Rmax: Any, Z: Any, *, n: Any = 64
) -> float:
    """Compute ``integral_0^Rmax F(rho,Z)^2 d rho`` on a finite interval.

    The scalar profile is sampled at an n-point Gauss--Legendre rule.  Inputs
    and returned profile values must be finite, and ``Rmax`` must be
    nonnegative.  This is precisely the finite ``M^p`` increment in Part I;
    it does not evaluate an infinite tail.

    Source: arXiv:2609.35406v1, Section 2.5, equation (2.22).
    """

    order = _order(n)
    radius = _finite_scalar(Rmax, "Rmax")
    z = _finite_scalar(Z, "Z")
    if radius < 0:
        raise ValueError("Rmax must be nonnegative")
    if not callable(F):
        raise TypeError("F must be callable as F(R, Z)")
    if radius == 0:
        return 0.0
    nodes, weights = leggauss(order)
    rho = radius * (nodes + 1.0) / 2.0
    values = np.empty(order, dtype=float)
    for index, point in enumerate(rho):
        try:
            value = float(F(float(point), z))
        except (TypeError, ValueError, IndexError) as exc:
            raise ValueError("F(R, Z) must return finite scalar values") from exc
        if not math.isfinite(value):
            raise ValueError("F(R, Z) returned a non-finite value")
        values[index] = value
    result = float(np.dot(radius * weights / 2.0, values * values))
    if not math.isfinite(result):
        raise ArithmeticError("finite pressure quadrature returned non-finite value")
    return result


def _flat_switch(s: float) -> float:
    """C-infinity switch from 0 to 1 on the unit interval."""

    if s <= 0.0:
        return 0.0
    if s >= 1.0:
        return 1.0
    left = math.exp(-1.0 / s)
    right = math.exp(-1.0 / (1.0 - s))
    return left / (left + right)


@dataclass(frozen=True)
class JoinedOuterPressure:
    """Join a supplied scalar core swirl to the checked heat exterior.

    ``F_core`` is used through ``R_core``.  Between ``R_core`` and
    ``R_join`` the profile is ``(1-switch) F_core + switch F_heat`` with the
    flat switch from :func:`_flat_switch`.  At and beyond ``R_join`` the
    profile is exactly ``F_heat = HeatExterior.profile_E(R,Z)/sqrt(2R)``.
    The class accepts any scalar finite core callable, so the
    ``PaperCoreReference`` seed is one executable instance rather than a
    hidden requirement of the interface.
    """

    F_core: Callable[[float, float], float]
    heat: HeatExterior
    R_core: float = 0.05
    R_join: float = 0.20

    def __post_init__(self) -> None:
        if not callable(self.F_core):
            raise TypeError("F_core must be callable as F_core(R, Z)")
        if not isinstance(self.heat, HeatExterior):
            raise TypeError("heat must be a HeatExterior instance")
        core = _finite_scalar(self.R_core, "R_core")
        join = _finite_scalar(self.R_join, "R_join")
        if core < 0 or not join > core:
            raise ValueError("require 0 <= R_core < R_join")

    @property
    def A(self) -> float:
        """The exterior exponent ``A=1/2+h``."""

        return self.heat.A

    def _coordinates(self, R: Any, Z: Any) -> tuple[float, float]:
        radius = _finite_scalar(R, "R")
        z = _finite_scalar(Z, "Z")
        if radius < 0 or abs(z) > 1.0:
            raise ValueError("require finite R>=0 and |Z|<=1")
        return radius, z

    def F_heat(self, R: Any, Z: Any) -> float:
        """Evaluate the exact heat regular profile at positive ``R``."""

        radius, z = self._coordinates(R, Z)
        if radius <= 0:
            raise ValueError("F_heat is undefined at R=0; use the core F")
        value = _cached_heat_F(self.heat, radius, z)
        if not math.isfinite(value):
            raise ValueError("F_heat returned a non-finite value")
        return float(value)

    def F(self, R: Any, Z: Any) -> float:
        """Evaluate the joined regular swirl profile ``F(R,Z)``.

        Heat evaluation is never attempted at ``R=0``.  The core is used on
        ``0 <= R <= R_core``; the flat blend reaches exact heat values at
        ``R_join``; and the heat profile is used for ``R >= R_join``.
        """

        radius, z = self._coordinates(R, Z)
        if radius <= self.R_core:
            value = self.F_core(radius, z)
        elif radius >= self.R_join:
            value = self.F_heat(radius, z)
        else:
            core_value = self.F_core(radius, z)
            heat_value = self.F_heat(radius, z)
            switch = _flat_switch(
                (radius - self.R_core) / (self.R_join - self.R_core)
            )
            value = (1.0 - switch) * core_value + switch * heat_value
        try:
            value = float(value)
        except (TypeError, ValueError) as exc:
            raise ValueError("F_core must return finite scalar values") from exc
        if not math.isfinite(value):
            raise ValueError("joined F returned a non-finite value")
        return value

    def heat_pressure(self, R: Any, Z: Any) -> float:
        """Return exact checked heat pressure for any finite ``R>0``.

        This delegates to ``HeatExterior.pressure_from_tau`` with
        ``r=sqrt(2R)`` and ``tau=1-Z^2``.  At ``R=R_join`` this is the exact
        terminal pressure used to set ``P0``.
        """

        radius, z = self._coordinates(R, Z)
        if radius <= 0:
            raise ValueError("heat pressure requires R>0")
        value = _cached_heat_pressure(self.heat, radius, z)
        if not math.isfinite(value):
            raise ValueError("heat pressure returned a non-finite value")
        return float(value)

    @lru_cache(maxsize=1024)
    def P0(self, Z: float, *, n: int = 64) -> float:
        """Return the joined axis pressure datum at finite quadrature order.

        ``P0(Z) = P_heat(R_join,Z) - integral_0^R_join F(rho,Z)^2 d rho``.
        The terminal pressure is the existing checked heat quadrature, while
        the subtraction is a finite Gauss--Legendre integral of this same
        joined profile.
        """

        z = _finite_scalar(Z, "Z")
        order = _order(n)
        terminal = self.heat_pressure(self.R_join, z)
        increment = finite_pressure_moment(self.F, self.R_join, z, n=order)
        return float(terminal - increment)

    def pressure(self, R: Any, Z: Any, *, n: Any = 64) -> float:
        """Evaluate joined pressure on finite interior and exact exterior.

        For ``0 <= R <= R_join`` this returns ``P0(Z) + integral_0^R F^2``.
        For ``R > R_join`` it delegates to exact ``HeatExterior`` pressure.
        The two branches agree at the join up to the finite interior
        quadrature error recorded by the caller.
        """

        radius, z = self._coordinates(R, Z)
        order = _order(n)
        if radius <= self.R_join:
            return self.P0(z, n=order) + finite_pressure_moment(
                self.F, radius, z, n=order
            )
        return self.heat_pressure(radius, z)


def build_reference_joined_profile(
    *, c_inf: float = 0.1, h: float = 0.001
) -> tuple[JoinedOuterPressure, Any]:
    """Build the requested ``PaperCoreReference``-to-heat executable seed.

    The returned tuple is ``(joined_profile, core_reference)`` so callers can
    report the autonomous core pressure and parameters without changing the
    generic joined-profile API.  The reference is local and finite only.
    """

    from openai_ns_reconstruction.paper_core_reference import PaperCoreReference

    reference = PaperCoreReference(
        h=h, j=0.02, sigma=0.3, Lambda=10.0, C=2.0, pressure_scale=1.0
    )
    heat = HeatExterior(h=h, c_inf=c_inf)
    return JoinedOuterPressure(reference.F, heat), reference


__all__ = [
    "JoinedOuterPressure",
    "SOURCE",
    "SOURCE_SECTION",
    "SOURCE_VERSION",
    "_flat_switch",
    "build_reference_joined_profile",
    "finite_pressure_moment",
]
