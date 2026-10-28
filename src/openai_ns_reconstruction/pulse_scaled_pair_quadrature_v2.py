"""Arbitrary-precision, normalized signed pulse moment on one directed chart.

The caller supplies same-coordinate frame and signed-difference ports; this
module does not bind or authenticate their sources. Radii are normalized by
exp(log_radius_scale), which is never materialized. The integral and boundary
have different scales and are returned separately, including ordinary Z1.
The quadrature estimate is diagnostic, not an enclosure or source error bound.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Any
import mpmath as mp


@dataclass(frozen=True)
class ScaledPairSample:
    """Values, not bounds. dr is independently supplied as rP-rF."""
    rP: Any
    rF: Any
    dr: Any
    S_F: Any
    S_F_Z: Any
    dS: Any
    dS_Z: Any
    Qy_F: Any
    Qy_F_Z: Any
    dQy: Any
    dQy_Z: Any
    measure: Any


@dataclass(frozen=True)
class ScaledPairResult:
    lower: Any
    upper: Any
    precision_dps: int
    integral_log_scale: Any
    boundary_log_scale: Any
    integral: tuple[Any, Any]
    lower_boundary: tuple[Any, Any]
    upper_boundary: tuple[Any, Any]
    quadrature_error_estimate: tuple[Any, Any]
    quadrature_error_certified: bool = False
    source_identity_proved: bool = False


def _number(value: Any, name: str) -> Any:
    if isinstance(value, bool):
        raise ValueError(f"{name} must be finite and real")
    try:
        result = mp.mpf(value)
    except (ValueError, TypeError, OverflowError) as exc:
        raise ValueError(f"{name} must be finite and real") from exc
    if not mp.isfinite(result):
        raise ValueError(f"{name} must be finite and real")
    return result


def _sample(provider: Callable, x: Any) -> dict[str, Any]:
    raw = provider(x)
    if not isinstance(raw, ScaledPairSample):
        raise TypeError("provider must return ScaledPairSample")
    v = {name: _number(getattr(raw, name), name)
         for name in ScaledPairSample.__dataclass_fields__}
    if v['rP'] <= 0 or v['rF'] <= 0:
        raise ValueError("normalized radii must be positive")
    if v['measure'] <= 0:
        raise ValueError("physical dy/dx measure must be positive")
    return v


def integrate_scaled_pair(provider: Callable, lower: Any, upper: Any, *,
                          log_radius_scale: Any, difference_scale: Any, precision_dps: int,
                          radius_and_endpoints_Z_independent: bool,
                          measure_Z_independent: bool) -> ScaledPairResult:
    """Integrate one chart, retaining both oriented finite endpoint terms.

    The density ports supply S_F and dS=S_P-S_F, plus ordinary Z1. The
    radius ports supply rP,rF and dr=(R_P-R_F)/exp(log_radius_scale), with
    dr independently computed from source-consistent data. For dy=measure dx,
    the output integral is half the integral of
    rP**2*dS + dr*(rP+rF)*S_F, at log scale 2*log_radius_scale.
    Each endpoint is rP*dQy + dr*Qy_F, at log scale log_radius_scale.
    Both mantissas are divided by the positive Z-independent difference_scale;
    its logarithm is added to both returned scales. This keeps tiny supplied
    differences at order one without evaluating exp(log_radius_scale).
    Thus m15_P-m15_F is the scaled integral minus (upper minus lower).
    No exp(log_radius_scale), source identity, source errors or seam
    cancellation is inferred. Endpoints, radii, scale and measure must be
    ordinary-Z-independent. The provider must be deterministic. The returned
    mp.quad estimate does not certify truncation, source or rounding errors.
    """
    if not callable(provider):
        raise TypeError("provider must be callable")
    if isinstance(precision_dps, bool) or not isinstance(precision_dps, int) or precision_dps < 15:
        raise ValueError("precision_dps must be an integer >= 15")
    if radius_and_endpoints_Z_independent is not True or measure_Z_independent is not True:
        raise ValueError("Z-independent radii, scale, endpoints and measure are required")
    with mp.workdps(precision_dps):
        a = _number(lower, 'lower')
        b = _number(upper, 'upper')
        scale = _number(log_radius_scale, 'log_radius_scale')
        difference = _number(difference_scale, 'difference_scale')
        if difference <= 0:
            raise ValueError('difference_scale must be positive')
        difference_log = mp.log(difference)
        def component(x: Any, order: int) -> Any:
            v = _sample(provider, x)
            sf = v['S_F'] if order == 0 else v['S_F_Z']
            ds = v['dS'] if order == 0 else v['dS_Z']
            return v['measure'] * (v['rP']**2 * ds + v['dr'] * (v['rP'] + v['rF']) * sf) / (2 * difference)
        def endpoint(x: Any) -> tuple[Any, Any]:
            v = _sample(provider, x)
            return ((v['rP'] * v['dQy'] + v['dr'] * v['Qy_F']) / difference,
                    (v['rP'] * v['dQy_Z'] + v['dr'] * v['Qy_F_Z']) / difference)
        lo, hi = endpoint(a), endpoint(b)
        if a == b:
            values = (mp.mpf(0), mp.mpf(0))
            errors = (mp.mpf(0), mp.mpf(0))
        else:
            first = mp.quad(lambda x: component(x, 0), [a, b], error=True)
            second = mp.quad(lambda x: component(x, 1), [a, b], error=True)
            values = (_number(first[0], 'integral Z0'), _number(second[0], 'integral Z1'))
            errors = (_number(first[1], 'diagnostic Z0'), _number(second[1], 'diagnostic Z1'))
        return ScaledPairResult(a, b, precision_dps, 2 * scale + difference_log, scale + difference_log, values, lo, hi, errors)
