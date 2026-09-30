"""Exact heat-exterior tail moments for the Lei--Ren Part I bookkeeping.

Source pinned to Z. Lei and X. Ren, arXiv:2609.35406v1, Sections 2.5,
4.23 and 5.1.  For the exact heat swirl, with ``d=1-Z^2``,
``k=2d/R``, ``A=1/2+h`` and ``H=heat_factor``, this module evaluates the
finite-radius tail data

* ``M^p_tail = -P_heat(R,Z)``;
* ``M^z_tail = M^theta_z_tail = 0``;
* ``M^ztheta_tail = -c_inf^2 R^(-2h)/2 * integral_0^1
  w^(2h-1) H(k w)^2 dw``; and
* the renormalized angular target from equation (4.23), using the equivalent
  tail-difference integral.

The pressure tail delegates to the repository's checked ``HeatExterior``
quadrature.  Small-argument remainder integrands use finite Taylor jets of
``H`` only where ``xi <= 0.01`` (order 8) or ``xi <= 0.02`` (order 12), which
avoids subtracting nearly equal numbers.  These are numerical diagnostics for
the exact supplied heat formula; they are not a global certificate, a core
matching result, a PDE validation, or scale-recursion evidence.
"""

from __future__ import annotations

import math
import operator
import sys
from functools import lru_cache
from pathlib import Path
from typing import Any

import numpy as np
from scipy.integrate import quad


_SRC = Path(__file__).resolve().parents[2] / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

from openai_ns_reconstruction.heat_exterior import HeatExterior, heat_factor


SOURCE = "https://arxiv.org/html/2609.35406v1"
SOURCE_VERSION = "2609.35406v1"
SOURCE_SECTIONS = "Sections 2.5, 4.23 and 5.1"


def _finite_scalar(value: Any, name: str) -> float:
    """Convert one input to a finite real scalar."""

    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be a finite real scalar") from exc
    if not math.isfinite(result):
        raise ValueError(f"{name} must be finite")
    return result


def _positive_order(value: Any, name: str = "series_order") -> int:
    """Validate a positive integer quadrature/Taylor order."""

    if isinstance(value, (bool, np.bool_)):
        raise TypeError(f"{name} must be a positive integer")
    try:
        result = operator.index(value)
    except TypeError as exc:
        raise TypeError(f"{name} must be a positive integer") from exc
    if result <= 0:
        raise ValueError(f"{name} must be a positive integer")
    return int(result)


def _validate_inputs(R: Any, Z: Any, heat: HeatExterior) -> tuple[float, float]:
    """Validate tail coordinates and the existing heat object."""

    if not isinstance(heat, HeatExterior):
        raise TypeError("heat must be an existing HeatExterior instance")
    radius = _finite_scalar(R, "R")
    z = _finite_scalar(Z, "Z")
    if radius <= 0:
        raise ValueError("tail moments require finite R>0")
    if abs(z) > 1.0:
        raise ValueError("tail moments require |Z|<=1")
    return radius, z


@lru_cache(maxsize=65536)
def _cached_heat_factor(xi: float, h: float) -> float:
    """Cache scalar heat factors during the nested tail quadratures."""

    return float(heat_factor(float(xi), float(h)))


def _heat_series(xi: float, h: float, order: int) -> float:
    """Evaluate the requested finite smooth jet of H at xi=0."""

    # H(xi) = sum_m (-1)^m (h)_m (1+h)_m xi^m / m!.
    total = 1.0
    coefficient = 1.0
    power = 1.0
    for m in range(1, order + 1):
        coefficient *= -(h + m - 1.0) * (1.0 + h + m - 1.0) / m
        power *= xi
        total += coefficient * power
    return float(total)


def _heat_series_remainder(
    xi: float, h: float, order: int, first_power: int
) -> float:
    """Sum selected Taylor terms directly, without subtractive cancellation."""

    if xi == 0.0 or first_power > order:
        return 0.0
    total = 0.0
    coefficient = 1.0
    power = 1.0
    for m in range(1, order + 1):
        coefficient *= -(h + m - 1.0) * (1.0 + h + m - 1.0) / m
        power *= xi
        if m >= first_power:
            total += coefficient * power
    return float(total)


def _small_xi_cutoff(order: int) -> float:
    """Use the requested xi thresholds for order 8 and order 12 jets."""

    return 0.01 if order <= 8 else 0.02


def _h_minus_one(xi: float, h: float, order: int) -> float:
    """Stable H(xi)-1 using the finite jet at sufficiently small xi."""

    if xi <= _small_xi_cutoff(order):
        return _heat_series_remainder(xi, h, order, 1)
    return _cached_heat_factor(xi, h) - 1.0


def _quadratic_remainder(xi: float, h: float, order: int) -> float:
    """Stable H(xi)^2-1 for the pressure/energy tail correction."""

    if xi <= _small_xi_cutoff(order):
        # Compute H-1 directly, then factor H^2-1=(H-1)(H+1).  This avoids
        # subtracting two unit-scale numbers at the singular endpoint.
        h_minus_one = _heat_series_remainder(xi, h, order, 1)
        return h_minus_one * (2.0 + h_minus_one)
    value = _cached_heat_factor(xi, h)
    return (value - 1.0) * (value + 1.0)


def _checked_unit_integral(function, *, rtol: float = 2e-10,
                           atol: float = 2e-12,
                           split: float | None = None,
                           endpoint_power: float | None = None) -> float:
    """Integrate an improper unit-interval function with endpoint mapping.

    The remainder integrands behave like ``w**endpoint_power`` at zero.  A
    power substitution removes the nearly singular derivative that otherwise
    makes adaptive quadrature report a pessimistic roundoff estimate for
    ``h=.001``.
    """

    def integrate_interval(lower: float, upper: float) -> tuple[float, float]:
        result = quad(
            function,
            lower,
            upper,
            epsabs=atol,
            epsrel=rtol,
            limit=300,
        )
        return float(result[0]), float(result[1])

    if split is not None and 0.0 < split < 1.0:
        if endpoint_power is None:
            lower_value, lower_error = integrate_interval(0.0, split)
        else:
            exponent = 1.0 / (endpoint_power + 1.0)

            def transformed(parameter: float) -> float:
                if parameter <= 0.0:
                    return 0.0
                w = split * parameter**exponent
                jacobian = split * exponent * parameter ** (exponent - 1.0)
                return function(w) * jacobian

            lower_result = quad(
                transformed,
                0.0,
                1.0,
                epsabs=atol,
                epsrel=rtol,
                limit=300,
            )
            lower_value, lower_error = float(lower_result[0]), float(lower_result[1])
        upper_value, upper_error = integrate_interval(split, 1.0)
        value = lower_value + upper_value
        error = lower_error + upper_error
    elif endpoint_power is not None:
        exponent = 1.0 / (endpoint_power + 1.0)

        def transformed(parameter: float) -> float:
            if parameter <= 0.0:
                return 0.0
            w = parameter**exponent
            jacobian = exponent * parameter ** (exponent - 1.0)
            return function(w) * jacobian

        result = quad(
            transformed,
            0.0,
            1.0,
            epsabs=atol,
            epsrel=rtol,
            limit=300,
        )
        value, error = float(result[0]), float(result[1])
    else:
        value, error = integrate_interval(0.0, 1.0)
    if not (math.isfinite(value) and math.isfinite(error)):
        raise ArithmeticError("heat-tail quadrature returned non-finite output")
    if error > 20.0 * max(atol, rtol * abs(value)):
        raise ArithmeticError(
            f"heat-tail quadrature error estimate too large: {error:g}"
        )
    return value


def heat_tail_moments(
    R: Any,
    Z: Any,
    heat: HeatExterior,
    *,
    series_order: Any = 12,
) -> dict[str, float]:
    """Return exact heat-tail moments at finite ``(R,Z)``.

    Parameters
    ----------
    R, Z:
        Similarity radius and axial coordinate with ``R>0`` and ``|Z|<=1``.
    heat:
        Existing :class:`HeatExterior`; its ``h`` and ``c_inf`` are used
        directly, including the checked pressure quadrature.
    series_order:
        Finite order of the small-``xi`` smooth jet.  Orders 8 and 12 are the
        registered comparison choices; lower/higher positive orders are
        accepted for diagnostics.

    Returns
    -------
    dict
        ``angular`` is the renormalized target from (4.23), ``axial`` and
        ``mixed`` are exactly zero for the pure heat swirl, ``quadratic`` is
        ``M^{ztheta}_tail``, and ``pressure`` is ``M^p_tail=-P_heat``.
        ``angular_tail_difference`` is the scaled tail-difference term
        subtracted from the leading angular term, while
        ``angular_tail_difference_bracket`` and ``quadratic_integral`` expose
        the finite unit-interval quantities used to form those values.

    This evaluates the exact supplied heat profile numerically on a finite
    lower radius.  It does not assert that a core is matched to the tail or
    that any global Part I moment condition, PDE, energy, or recursion holds.
    """

    radius, z = _validate_inputs(R, Z, heat)
    order = _positive_order(series_order)
    h = _finite_scalar(heat.h, "heat.h")
    c_inf = _finite_scalar(heat.c_inf, "heat.c_inf")
    d = 1.0 - z * z
    k = 2.0 * d / radius
    b = h * (1.0 + h)

    cutoff = _small_xi_cutoff(order)
    split = cutoff / k if k > 0.0 else None

    def angular_integrand(w: float) -> float:
        if w <= 0.0:
            return 0.0
        xi = k * w
        remainder = _h_minus_one(xi, h, order) + b * xi
        return w ** (h - 2.0) * remainder

    angular_tail_difference_bracket = (
        0.0
        if k == 0.0
        else _checked_unit_integral(
            angular_integrand, split=split, endpoint_power=h
        )
    )
    angular_tail_difference_bracket += -b * k / h
    angular_tail_prefactor = (
        math.sqrt(2.0) * c_inf * radius ** (1.0 - h)
    )
    angular_tail_difference = (
        angular_tail_prefactor * angular_tail_difference_bracket
    )
    angular_target = (
        math.sqrt(2.0) * c_inf / (1.0 - h) * radius ** (1.0 - h)
        - angular_tail_difference
    )

    def quadratic_integrand(w: float) -> float:
        if w <= 0.0:
            return 0.0
        xi = k * w
        remainder = _quadratic_remainder(xi, h, order)
        return w ** (2.0 * h - 1.0) * remainder

    quadratic_correction = (
        0.0
        if k == 0.0
        else _checked_unit_integral(
            quadratic_integrand, split=split, endpoint_power=2.0 * h
        )
    )
    quadratic_integral = 1.0 / (2.0 * h) + quadratic_correction
    quadratic_tail = -0.5 * c_inf**2 * radius ** (-2.0 * h) * quadratic_integral

    # Existing HeatExterior.pressure_from_tau uses the checked quadrature in
    # the exact heat pressure formula.  The pressure moment is its negative.
    pressure_tail = -float(
        heat.pressure_from_tau(math.sqrt(2.0 * radius), 1.0 - z * z)
    )
    if not all(
        math.isfinite(value)
        for value in (
            angular_target,
            angular_tail_difference,
            angular_tail_difference_bracket,
            quadratic_tail,
            quadratic_integral,
            pressure_tail,
        )
    ):
        raise ArithmeticError("heat-tail moments returned non-finite output")

    return {
        "angular": float(angular_target),
        "angular_target": float(angular_target),
        "axial": 0.0,
        "mixed": 0.0,
        "quadratic": float(quadratic_tail),
        "quadratic_tail": float(quadratic_tail),
        "pressure": float(pressure_tail),
        "pressure_tail": float(pressure_tail),
        "angular_tail_difference": float(angular_tail_difference),
        "taildiff": float(angular_tail_difference),
        "angular_taildiff": float(angular_tail_difference),
        "angular_tail_difference_bracket": float(
            angular_tail_difference_bracket
        ),
        "quadratic_integral": float(quadratic_integral),
        "series_order": order,
        "R": radius,
        "Z": z,
    }


def direct_heat_tail_integrals(
    R: Any, Z: Any, heat: HeatExterior
) -> dict[str, float]:
    """Independently direct-integrate the two heat tail unit integrals.

    This diagnostic intentionally evaluates ``H`` directly in the integrands,
    without the small-``xi`` remainder factorization.  It is intended for the
    larger-h comparison in the checks, where cancellation is less severe; it
    is not a replacement for the stable primary evaluator.
    """

    radius, z = _validate_inputs(R, Z, heat)
    h = _finite_scalar(heat.h, "heat.h")
    d = 1.0 - z * z
    k = 2.0 * d / radius
    b = h * (1.0 + h)

    def direct_angular_regular(w: float) -> float:
        """The regular factor in w^(h-1) * ((H-1)/w)."""

        if w <= 0.0:
            return -b * k
        xi = k * w
        # Direct H evaluation is retained away from the endpoint.  At the
        # final machine-scale endpoint use its analytic first derivative so
        # this independent quadrature does not divide roundoff by w.
        if xi < 1.0e-7:
            return -b * k
        return (_cached_heat_factor(xi, h) - 1.0) / w

    if k == 0.0:
        direct_angular_difference_bracket = 0.0
    else:
        direct_result = quad(
            direct_angular_regular,
            0.0,
            1.0,
            epsabs=2e-11,
            epsrel=2e-10,
            limit=300,
            weight="alg",
            wvar=(h - 1.0, 0.0),
        )
        direct_angular_difference_bracket, direct_error = map(float, direct_result[:2])
        if not math.isfinite(direct_angular_difference_bracket) or not math.isfinite(direct_error):
            raise ArithmeticError("direct angular tail quadrature failed")
    direct_angular_prefactor = (
        math.sqrt(2.0) * heat.c_inf * radius ** (1.0 - h)
    )
    direct_angular_difference = (
        direct_angular_prefactor
        * direct_angular_difference_bracket
    )
    direct_angular_target = (
        math.sqrt(2.0) * heat.c_inf / (1.0 - h) * radius ** (1.0 - h)
        - direct_angular_difference
    )

    if k == 0.0:
        direct_quadratic_integral = 1.0 / (2.0 * h)
    else:
        direct_result = quad(
            lambda w: _cached_heat_factor(k * w, h) ** 2,
            0.0,
            1.0,
            epsabs=2e-11,
            epsrel=2e-10,
            limit=300,
            weight="alg",
            wvar=(2.0 * h - 1.0, 0.0),
        )
        direct_quadratic_integral, direct_error = map(float, direct_result[:2])
        if not math.isfinite(direct_quadratic_integral) or not math.isfinite(direct_error):
            raise ArithmeticError("direct quadratic tail quadrature failed")
    direct_quadratic_tail = (
        -0.5 * heat.c_inf**2 * radius ** (-2.0 * h) * direct_quadratic_integral
    )
    return {
        "angular": float(direct_angular_target),
        "angular_tail_difference": float(direct_angular_difference),
        "taildiff": float(direct_angular_difference),
        "angular_taildiff": float(direct_angular_difference),
        "angular_tail_difference_bracket": float(
            direct_angular_difference_bracket
        ),
        "quadratic": float(direct_quadratic_tail),
        "quadratic_integral": float(direct_quadratic_integral),
    }


__all__ = [
    "SOURCE",
    "SOURCE_SECTIONS",
    "SOURCE_VERSION",
    "direct_heat_tail_integrals",
    "heat_tail_moments",
]
