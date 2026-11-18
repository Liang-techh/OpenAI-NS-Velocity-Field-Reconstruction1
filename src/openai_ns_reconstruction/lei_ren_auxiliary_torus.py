"""Finite Fourier and radial operators for Lei–Ren Part II, §§2.5 and 4.2.

The coefficients are supplied at fixed slow coordinates. These operators do
not infer them from physical samples or certify the background assumptions.
"""
from __future__ import annotations

from collections.abc import Callable, Mapping
from dataclasses import dataclass
import math
from numbers import Integral

from scipy.integrate import quad_vec

from .charts import J_G, LAMBDA_G, T_G, V_R

Mode = tuple[int, int]
Coefficients = Mapping[Mode, complex]


def _level(level: int) -> int:
    if isinstance(level, bool) or not isinstance(level, Integral) or level < 0:
        raise ValueError("covering levels must be nonnegative integers")
    return int(level)


def _mode(mode: Mode) -> Mode:
    if len(mode) != 2 or any(isinstance(n, bool) or not isinstance(n, Integral) for n in mode):
        raise ValueError("expected a two-dimensional integer Fourier mode")
    return int(mode[0]), int(mode[1])


def pullback(coefficients: Coefficients, *, level: int, common_level: int) -> dict[Mode, complex]:
    """Pull back by E**(level-common_level) via exact integer mode indices."""
    gap = _level(level) - _level(common_level)
    if gap < 0:
        raise ValueError("the common level cannot exceed the source level")
    a, b = J_G.T[0]
    c, d = J_G.T[1]
    result = {}
    for mode, coefficient in coefficients.items():
        n, m = _mode(mode)
        for _ in range(gap):
            n, m = a*n + b*m, c*n + d*m
        result[n, m] = coefficient
    return result


def _direction_dot(n: int, m: int, direction: str) -> float:
    # v·(n,m) = a+b sqrt(2). Preserve the small divisor when the terms
    # nearly cancel by using the exact integer norm a²-2b².
    a, b = (n + m, -m) if direction == "radial" else (m - n, n)
    root_term = b * math.sqrt(2.0)
    if a * b < 0:
        return (a*a - 2*b*b) / (a - root_term)
    return a + root_term


def inverse_directional_derivative(
    coefficients: Coefficients, *, level: int, direction: str = "time"
) -> dict[Mode, complex]:
    """Invert Λ**level v·∇_Y on finite zero-mean Fourier data (Lemma 2.5).

    To cancel the oscillatory part of a tangential source, pass its negative
    after taking out the zero mode. The solver does not discard a nonzero mean.
    """
    level = _level(level)
    if direction not in {"radial", "time"}:
        raise ValueError("direction must be radial or time")
    eigenvalue = LAMBDA_G if direction == "radial" else T_G
    inverse_scale = eigenvalue ** (-level)
    if inverse_scale == 0:
        raise ArithmeticError("inverse covering scale underflowed")
    result = {}
    for mode, coefficient in coefficients.items():
        n, m = _mode(mode)
        if (n, m) == (0, 0):
            if coefficient != 0:
                raise ValueError("directional inverse requires zero mean")
        else:
            result[n, m] = (coefficient / (2j * math.pi * _direction_dot(n, m, direction))) * inverse_scale
    return result


@dataclass(frozen=True)
class RadialPrimitive:
    value: complex
    physical_radial_derivative: complex
    quadrature_error_estimate: float
    derivative_error_estimate: float


def transported_radial_primitive(
    source: Callable[[float, tuple[float, float]], complex],
    *,
    radius: float,
    auxiliary: tuple[float, float],
    annulus: tuple[float, float],
    cutoff: Callable[[float], tuple[float, float]],
    level: int,
    radial_exponent: float,
    atol: float = 1e-11,
    rtol: float = 1e-11,
) -> RadialPrimitive:
    """Evaluate I₀h and D_r I₀h in Eq. (4.33), at fixed slow t,z.

    The caller supplies h supported in the annulus and a smooth cutoff giving
    (χ, ∂rχ), equal to (0,0) near the inner edge and (1,0) near the outer edge.
    For Eq. (4.34) use h=r B_aux and Ψ=I₀h/r. The derivative uses the
    transported FTC, so it differentiates neither a numerical quadrature nor
    the torus phase. Returned quadrature estimates are numerical diagnostics.
    """
    level = _level(level)
    lo, hi = annulus
    y0, y1 = auxiliary
    if (
        not all(math.isfinite(v) for v in (lo, hi, radius, y0, y1, radial_exponent, atol, rtol))
        or not 0 < lo < hi or radius < 0 or min(radial_exponent, atol, rtol) <= 0
    ):
        raise ValueError("expected finite positive scales and a positive annulus")
    if radius <= lo or radius >= hi:
        return RadialPrimitive(0j, 0j, 0.0, 0.0)

    weight, weight_r = cutoff(radius)
    if not 0 <= weight <= 1 or not math.isfinite(weight_r):
        raise ValueError("cutoff must return a weight in [0,1] and finite radial derivative")
    radial_scale = math.exp(level * math.log(LAMBDA_G) + radial_exponent * math.log(radius))
    y0, y1 = y0 % 1, y1 % 1

    def transported(s: float) -> complex:
        offset = (s - radius) / radius
        log_ratio = math.log1p(offset) if abs(offset) <= 0.5 else math.log(s) - math.log(radius)
        shift = radial_scale * math.expm1(radial_exponent * log_ratio)
        torus = ((y0 + shift*V_R[0]) % 1, (y1 + shift*V_R[1]) % 1)
        return complex(source(s, torus))

    def integrate(a: float, b: float) -> tuple[complex, float]:
        value, error, info = quad_vec(transported, a, b, epsabs=atol, epsrel=rtol, full_output=True)
        if not info.success:
            raise ArithmeticError(f"transported radial quadrature failed (status {info.status})")
        return complex(value), float(error)

    inner, inner_error = integrate(lo, radius)
    outer, outer_error = integrate(radius, hi)
    return RadialPrimitive(
        (1 - weight)*inner - weight*outer,
        complex(source(radius, (y0, y1))) - weight_r*(inner + outer),
        (1 - weight)*inner_error + weight*outer_error,
        abs(weight_r)*(inner_error + outer_error),
    )
