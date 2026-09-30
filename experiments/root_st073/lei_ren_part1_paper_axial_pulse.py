"""The fixed Lei--Ren Part I axial pulse from Section 7.5.

This module materializes only the source pulse shape in (7.28).  It does not
construct the incoming moment right-hand sides in (7.31), the corrected
angular tail in (7.34), or the full axial closure.  The stage-local API uses
``xi = mu*t`` directly and never constructs the enormous absolute radius or
the collapsed coordinate ``13/mu``.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from functools import lru_cache
from typing import Any

import numpy as np
from scipy.integrate import quad


SOURCE = "https://arxiv.org/html/2609.35406v1"
SOURCE_VERSION = "2609.35406v1"
SOURCE_SECTION = "Part I Section 7.5, equation (7.28); definitions from (4.1)"
SOURCE_EQUATION = "g_p(xi)=[1-sigma(xi-10)]*int_0^xi sigma(50v)dv for xi>=0, else 0"
SOURCE_KP_BOUND = (0.24, 0.246)


def _finite(value: Any, name: str) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be finite") from exc
    if not math.isfinite(result):
        raise ValueError(f"{name} must be finite")
    return result


def flat_function(t: Any) -> float:
    """The source ``mathfrak f(t)`` from (4.1)."""

    value = _finite(t, "t")
    if value <= 0.0:
        return 0.0
    return math.exp(-1.0 / (value * value))


def smooth_step(t: Any) -> float:
    """The source ``sigma(t)`` from (4.1)."""

    value = _finite(t, "t")
    if value <= 0.0:
        return 0.0
    if value >= 1.0:
        return 1.0
    left = flat_function(value)
    right = flat_function(1.0 - value)
    denominator = left + right
    if denominator == 0.0:
        return 0.0 if value < 0.5 else 1.0
    return left / denominator


def smooth_step_derivative(t: Any) -> float:
    """Derivative of the source smooth step on its transition interval."""

    value = _finite(t, "t")
    if value <= 0.0 or value >= 1.0:
        return 0.0
    left = flat_function(value)
    right = flat_function(1.0 - value)
    denominator = left + right
    if denominator == 0.0:
        return 0.0
    left_prime = 2.0 * left / (value**3)
    right_prime_in_value = 2.0 * right / ((1.0 - value) ** 3)
    return (left_prime * right + left * right_prime_in_value) / (denominator * denominator)


def _primitive_cutoff() -> float:
    value, error = quad(lambda v: smooth_step(50.0 * v), 0.0, 1.0 / 50.0, epsabs=2.0e-14, epsrel=2.0e-14, limit=200)
    if not math.isfinite(value) or error > 1.0e-12:
        raise RuntimeError(f"source pulse primitive quadrature failed: value={value}, error={error}")
    return float(value)


_J50_CUTOFF = _primitive_cutoff()


@lru_cache(maxsize=1024)
def pulse_primitive(xi: float) -> float:
    """Return ``J_50(xi)=int_0^xi sigma(50v)dv`` from (7.28)."""

    value = _finite(xi, "xi")
    if value <= 0.0:
        return 0.0
    cutoff = 1.0 / 50.0
    if value >= cutoff:
        return float(_J50_CUTOFF + value - cutoff)
    result, error = quad(lambda v: smooth_step(50.0 * v), 0.0, value, epsabs=2.0e-14, epsrel=2.0e-14, limit=200)
    if not math.isfinite(result) or error > 1.0e-12:
        raise RuntimeError(f"source pulse primitive quadrature failed: value={result}, error={error}")
    return float(result)


def pulse_primitive_derivative(xi: Any) -> float:
    """Derivative of the source primitive, namely ``sigma(50*xi)``."""

    value = _finite(xi, "xi")
    if value <= 0.0:
        return 0.0
    return smooth_step(50.0 * value)


def pulse_value(xi: Any) -> float:
    """Evaluate the source ``g_p(xi)`` with support in ``[0,11]``."""

    value = _finite(xi, "xi")
    if value <= 0.0 or value >= 11.0:
        return 0.0
    return (1.0 - smooth_step(value - 10.0)) * pulse_primitive(value)


def pulse_derivative(xi: Any) -> float:
    """Evaluate the analytic first derivative of ``g_p``."""

    value = _finite(xi, "xi")
    if value <= 0.0 or value >= 11.0:
        return 0.0
    tail_step = smooth_step(value - 10.0)
    return (
        -smooth_step_derivative(value - 10.0) * pulse_primitive(value)
        + (1.0 - tail_step) * pulse_primitive_derivative(value)
    )


@dataclass(frozen=True)
class PaperAxialPulse:
    """Callable source pulse with a stage-local ``mu*t`` view."""

    source: str = SOURCE
    source_version: str = SOURCE_VERSION

    def value(self, xi: Any) -> float:
        return pulse_value(xi)

    def derivative(self, xi: Any) -> float:
        return pulse_derivative(xi)

    def primitive(self, xi: Any) -> float:
        return pulse_primitive(_finite(xi, "xi"))

    def value_at_stage(self, t: Any, mu: Any) -> float:
        """Evaluate ``g_p(mu*t)`` without constructing an absolute radius."""

        t_value = _finite(t, "t")
        mu_value = _finite(mu, "mu")
        if mu_value <= 0.0:
            raise ValueError("mu must be positive")
        return self.value(mu_value * t_value)

    def derivative_at_stage(self, t: Any, mu: Any) -> float:
        """Return ``d/dt g_p(mu*t)=mu*g_p'(mu*t)``."""

        t_value = _finite(t, "t")
        mu_value = _finite(mu, "mu")
        if mu_value <= 0.0:
            raise ValueError("mu must be positive")
        return mu_value * self.derivative(mu_value * t_value)

    def K_p(self, *, quadrature_order: int = 128) -> float:
        return pulse_K_p(quadrature_order=quadrature_order)


def pulse_K_p(*, quadrature_order: int = 128) -> float:
    """Quadrature for ``K_p=int_0^13 exp(-2xi) g_p(xi)^2 dxi``."""

    try:
        order = int(quadrature_order)
    except (TypeError, ValueError) as exc:
        raise ValueError("quadrature_order must be an integer at least 16") from exc
    if order < 16 or order != quadrature_order:
        raise ValueError("quadrature_order must be an integer at least 16")
    # Split at both flat transition locations and include the zero tail.
    intervals = ((0.0, 0.02), (0.02, 10.0), (10.0, 10.02), (10.02, 11.0), (11.0, 13.0))
    nodes, weights = np.polynomial.legendre.leggauss(order)
    total = 0.0
    for left, right in intervals:
        midpoint = 0.5 * (left + right)
        half_width = 0.5 * (right - left)
        xi = midpoint + half_width * nodes
        values = np.asarray([pulse_value(float(value)) for value in xi], dtype=float)
        total += half_width * float(np.sum(weights * np.exp(-2.0 * xi) * values * values))
    return float(total)


__all__ = [
    "SOURCE",
    "SOURCE_VERSION",
    "SOURCE_SECTION",
    "SOURCE_EQUATION",
    "SOURCE_KP_BOUND",
    "flat_function",
    "smooth_step",
    "smooth_step_derivative",
    "pulse_primitive",
    "pulse_primitive_derivative",
    "pulse_value",
    "pulse_derivative",
    "PaperAxialPulse",
    "pulse_K_p",
]
