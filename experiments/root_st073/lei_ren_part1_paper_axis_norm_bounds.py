"""Conservative real interval bounds for the axis profile derivatives.

This is only the real ``Z in [-1, 1]`` contribution to an axis norm.  It
does not bound mixed radial derivatives, a complex ``A_Omega`` norm, the
pressure amplitude, or the nonlinear core.  For an injected
``RegularCoreAxisJets`` object, ``axis_norm_bounds(axis)`` bounds
``G, G', G'', G'''`` using coefficient absolute sums for ``H`` and the
positive denominator ``H^2 + sigma0^2``.

Writing ``f(H)=H/(H^2+sigma0^2)``, the quotient derivatives used here are

    f'  = (s^2-H^2)/Q^2,
    f'' = 2 H (H^2-3s^2)/Q^3,
    f'''= -6 (H^4-6H^2s^2+s^4)/Q^4,

where ``Q=H^2+s^2``.  The bounds use the elementary real inequalities

    |f| <= 1/(2s), |f'| <= 1/s^2,
    |f''| <= 3/(2s^3), |f'''| <= 48/s^4.

The last constant is intentionally loose and keeps the estimate finite and
transparent.  No complex-domain or PDE certification is implied.
"""

from __future__ import annotations

from collections.abc import Mapping
import json
from pathlib import Path
import sys
from typing import Any

import mpmath as mp


HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from lei_ren_part1_paper_axis_jets import RegularCoreAxisJets  # noqa: E402


SOURCE = "https://arxiv.org/html/2609.35406v1"
SOURCE_VERSION = "2609.35406v1"


def _mp(value: Any) -> mp.mpf:
    if isinstance(value, mp.mpf):
        return mp.mpf(value)
    return mp.mpf(str(value))


def _nstr(value: Any, digits: int = 50) -> str:
    return mp.nstr(_mp(value), digits)


def _poly_derivative_abs_sum(coefficients: Mapping[int, Any], order: int) -> mp.mpf:
    """Coefficient absolute-sum bound on [-1, 1] for one polynomial derivative."""

    total = mp.mpf(0)
    for degree, coefficient in coefficients.items():
        if degree < order:
            continue
        falling = mp.mpf(1)
        for step in range(order):
            falling *= degree - step
        total += abs(_mp(coefficient)) * falling
    return total


def _h_coefficients(axis: Any) -> dict[int, mp.mpf]:
    """Exact coefficients of H0=(1-delta)Z/2+(1-Z^2)(4Z+j)."""

    # H0 = -4 Z^3 - j Z^2 + (9/2-delta/2) Z + j.
    return {
        0: axis.j,
        1: mp.mpf(9) / 2 - axis.delta / 2,
        2: -axis.j,
        3: mp.mpf(-4),
    }


def _l_coefficients(axis: Any) -> dict[int, mp.mpf]:
    return {0: mp.mpf(1), 2: -axis.delta}


def _real_polynomial_bounds(axis: Any) -> dict[str, list[mp.mpf]]:
    h_coefficients = _h_coefficients(axis)
    l_coefficients = _l_coefficients(axis)
    h = [_poly_derivative_abs_sum(h_coefficients, order) for order in range(4)]
    l = [_poly_derivative_abs_sum(l_coefficients, order) for order in range(4)]
    return {"H_derivative_abs_sum": h, "L_derivative_abs_sum": l}


def _quotient_derivative_bounds(sigma0: mp.mpf) -> list[mp.mpf]:
    """Bounds for f(H), f'(H), f''(H), f'''(H) over real H."""

    s = abs(sigma0)
    if s == 0:
        raise ValueError("sigma0 must be nonzero for a real denominator bound")
    return [
        1 / (2 * s),
        1 / s ** 2,
        mp.mpf(3) / (2 * s ** 3),
        mp.mpf(48) / s ** 4,
    ]


def _g_derivative_bounds(axis: Any) -> list[mp.mpf]:
    """Leibniz bounds for G' = L f(H) through its third derivative."""

    polynomial = _real_polynomial_bounds(axis)
    h = polynomial["H_derivative_abs_sum"]
    l = polynomial["L_derivative_abs_sum"]
    f = _quotient_derivative_bounds(axis.sigma0)
    return [
        l[0] * f[0],
        l[1] * f[0] + l[0] * h[1] * f[1],
        l[2] * f[0] + 2 * l[1] * h[1] * f[1] + l[0] * (h[2] * f[1] + h[1] ** 2 * f[2]),
        l[3] * f[0]
        + 3 * l[2] * h[1] * f[1]
        + 3 * l[1] * (h[2] * f[1] + h[1] ** 2 * f[2])
        + l[0] * (h[3] * f[1] + 3 * h[1] * h[2] * f[2] + h[1] ** 3 * f[3]),
    ]


def _u0_c3_bound(axis: Any) -> mp.mpf:
    """C^3 sum bound for U0z=4Z+j on [-1,1]."""

    return (mp.mpf(4) + axis.j) + mp.mpf(4)


def axis_norm_bounds(axis: Any) -> dict[str, Any]:
    """Return conservative real-axis MP bounds for ``G`` and its derivatives."""

    precision = int(getattr(axis, "precision", 160))
    with mp.workdps(precision):
        if not hasattr(axis, "sigma0") or not hasattr(axis, "H0"):
            raise TypeError("axis must expose sigma0 and H0 from RegularCoreAxisJets")
        polynomial = _real_polynomial_bounds(axis)
        quotient = _quotient_derivative_bounds(axis.sigma0)
        g_bounds = _g_derivative_bounds(axis)
        g_value_upper = 2 * g_bounds[0]
        u0_c3 = _u0_c3_bound(axis)
        A_axis_upper = mp.mpf(10) + axis.Lambda * sum(g_bounds) + u0_c3
        return {
            "source": SOURCE,
            "source_version": SOURCE_VERSION,
            "interval": (mp.mpf(-1), mp.mpf(1)),
            "j": axis.j,
            "Lambda": axis.Lambda,
            "delta": axis.delta,
            "sigma0": axis.sigma0,
            "Z0": axis.Z0,
            "H_derivative_abs_sum": polynomial["H_derivative_abs_sum"],
            "L_derivative_abs_sum": polynomial["L_derivative_abs_sum"],
            "quotient_derivative_bounds": quotient,
            "G_derivative_bounds": g_bounds,
            "G_value_upper": g_value_upper,
            "G_actual_endpoint_values": {"-1": axis.G(-1), "1": axis.G(1)},
            "U0z_C3_sum_bound": u0_c3,
            "A_axis_upper": A_axis_upper,
            "scope": "real interval axis contribution only; no mixed radial, complex A_Omega, K, cone, PDE, or full-core certificate",
        }


def _direct_gradient(axis: Any, z: mp.mpf) -> mp.mpf:
    return axis.L(z) * axis.H0(z) / (axis.H0(z) ** 2 + axis.sigma0 ** 2)


def _derivative_checks(axis: Any) -> list[dict[str, Any]]:
    bounds = _g_derivative_bounds(axis)
    points = [mp.mpf("-.8"), mp.mpf("-.2"), axis.Z0, mp.mpf(".2"), mp.mpf(".8")]
    rows = []
    with mp.workdps(axis.precision):
        for point in points:
            values = []
            for order in range(4):
                actual = mp.diff(lambda value: _direct_gradient(axis, value), point, order)
                values.append(
                    {
                        "order": order,
                        "value": _nstr(actual),
                        "abs_value": _nstr(abs(actual)),
                        "bound": _nstr(bounds[order]),
                        "bound_covers": abs(actual) <= bounds[order],
                    }
                )
            rows.append({"Z": _nstr(point), "derivatives": values})
    return rows


def run_sample() -> dict[str, Any]:
    with mp.workdps(160):
        axis = RegularCoreAxisJets(
            j="1e-14",
            Lambda="1e36",
            delta="1e-200",
            precision=160,
        )
        result = axis_norm_bounds(axis)
        result["derivative_checks"] = _derivative_checks(axis)
        result["all_sample_derivative_checks_covered"] = all(
            item["bound_covers"]
            for row in result["derivative_checks"]
            for item in row["derivatives"]
        )
        return result


def _json_value(value: Any) -> Any:
    if isinstance(value, mp.mpf):
        return _nstr(value)
    if isinstance(value, Mapping):
        return {str(key): _json_value(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_json_value(item) for item in value]
    return value


if __name__ == "__main__":
    receipt = run_sample()
    output = Path(__file__).with_suffix(".json")
    output.write_text(json.dumps(_json_value(receipt), indent=2) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "A_axis_upper": _nstr(receipt["A_axis_upper"]),
                "G_value_upper": _nstr(receipt["G_value_upper"]),
                "all_sample_derivative_checks_covered": receipt["all_sample_derivative_checks_covered"],
            }
        ),
        flush=True,
    )


__all__ = ["axis_norm_bounds", "run_sample"]
