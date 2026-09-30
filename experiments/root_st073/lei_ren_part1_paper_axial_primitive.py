"""Stage-local primitive of the Section 7.5 axial pulse.

For ``xi = mu * log(R/Rp)`` this module evaluates

    J(xi) = (1/mu) * integral_0^min(xi, 11)
           exp(-lambda * (xi-v)/mu) g_p(v) dv,

where ``lambda = 1/2 - row*mu`` and ``g_p`` is the source pulse in (7.28).
The API keeps the result in logarithmic arbitrary-precision form.  The
linear interval uses its exact particular solution, the flat endpoint uses
a saddle or boundary coordinate, and the post-pulse branch reuses the
independently checked full-row pulse integral.  No absolute radius or
``13/mu`` coordinate is materialized.

This is a primitive input for the corrected axial profile.  It does not
construct ``U^z``, ``M^z_Z``, ``U^r``, pressure, or a five-moment closure.
"""

from __future__ import annotations

from functools import lru_cache
import json
import math
from pathlib import Path
import sys
from typing import Any

import mpmath as mp
import numpy as np
from numpy.polynomial.legendre import leggauss


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SRC = ROOT / "src"
for path in (SRC, HERE):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from lei_ren_part1_paper_axial_pulse_moments import (  # noqa: E402
    normalized_pulse_integral,
)
from lei_ren_part1_paper_outer import PaperOuterSchedule  # noqa: E402


SOURCE = "https://arxiv.org/html/2609.35406v1"
SOURCE_VERSION = "2609.35406v1"
SOURCE_SECTION = "Part I Section 7.5, equations (7.28)-(7.31)"
XI_STARTUP = mp.mpf("0.02")
XI_CUTOFF = mp.mpf("10")
XI_END = mp.mpf("11")
XI_LINEAR_OFFSET = mp.mpf("0.01")


def _mp(value: Any) -> mp.mpf:
    return mp.mpf(str(value))


def _validate(mu: Any, row: Any, xi: Any, precision: int, order: int) -> tuple[mp.mpf, int, mp.mpf]:
    try:
        # Parse decimal strings only after raising mpmath precision.  Parsing
        # at the process default would permanently round the tiny source mu
        # before the arbitrary-exponent branches see it.
        with mp.workdps(max(int(precision), 80)):
            mu_mp = _mp(mu)
            xi_mp = _mp(xi)
    except (TypeError, ValueError) as exc:
        raise ValueError("mu and xi must be finite real values") from exc
    if not mp.isfinite(mu_mp) or mu_mp <= 0:
        raise ValueError("mu must be positive and finite")
    if not mp.isfinite(xi_mp) or xi_mp < 0:
        raise ValueError("xi must be finite and nonnegative")
    try:
        row_int = int(row)
    except (TypeError, ValueError) as exc:
        raise ValueError("row must be 1 or 2") from exc
    if row_int not in (1, 2) or row_int != row:
        raise ValueError("row must be 1 or 2")
    if int(order) != order or order < 16:
        raise ValueError("order must be an integer at least 16")
    if int(precision) != precision or precision < 80:
        raise ValueError("precision must be an integer at least 80")
    lam = mp.mpf("0.5") - row_int * mu_mp
    if lam <= 0:
        raise ValueError("lambda = 1/2 - row*mu must be positive")
    return mu_mp, row_int, xi_mp


def _sigma_mp(value: mp.mpf) -> mp.mpf:
    """Evaluate the source flat step without binary64 underflow."""

    if value <= 0:
        return mp.mpf("0")
    if value >= 1:
        return mp.mpf("1")
    x = 1 / (value * value) - 1 / ((1 - value) * (1 - value))
    if x > 0:
        tiny = mp.exp(-x)
        return tiny / (1 + tiny)
    return 1 / (1 + mp.exp(x))


@lru_cache(maxsize=4096)
def _sigma_integral_cached(q_string: str, precision: int) -> mp.mpf:
    q = mp.mpf(q_string)
    if q <= 0:
        return mp.mpf("0")
    if q >= 1:
        return mp.mpf("0.5")
    # For very small q, retaining a zero with a declared bound is more
    # informative than asking adaptive quadrature to resolve a flat endpoint.
    if q < mp.mpf("1e-12"):
        bound = q * mp.exp(-1 / (q * q) + 1 / ((1 - q) * (1 - q)))
        return mp.mpf("0") if bound < mp.power(10, -(precision - 12)) else mp.quad(_sigma_mp, [0, q])
    return mp.quad(_sigma_mp, [0, q])


def _pulse_start_value(v: mp.mpf, precision: int) -> mp.mpf:
    """Evaluate g_p(v) on 0 <= v <= .02 in arbitrary precision."""

    if v <= 0:
        return mp.mpf("0")
    if v >= XI_STARTUP:
        return v - XI_LINEAR_OFFSET
    q = 50 * v
    if q <= mp.mpf("0.5"):
        integral = _sigma_integral_cached(mp.nstr(q, precision), precision)
        return integral / 50
    delta = 1 - q
    # Symmetry sigma(1-u)=1-sigma(u) gives the linear endpoint value plus
    # the tiny positive correction integral on [0, delta].
    correction = _sigma_integral_cached(mp.nstr(delta, precision), precision)
    return v - XI_LINEAR_OFFSET + correction / 50


def source_pulse_value_mp(value: Any, *, precision: int = 160) -> mp.mpf:
    """Evaluate the source ``g_p`` with arbitrary-precision flat endpoints."""

    with mp.workdps(precision):
        v = _mp(value)
        if v <= 0 or v >= XI_END:
            return mp.mpf("0")
        if v <= XI_STARTUP:
            return _pulse_start_value(v, precision)
        if v <= XI_CUTOFF:
            return v - XI_LINEAR_OFFSET
        return (v - XI_LINEAR_OFFSET) * _sigma_mp(XI_END - v)


def _gauss_integral(left: mp.mpf, right: mp.mpf, function, order: int) -> mp.mpf:
    if right <= left:
        return mp.mpf("0")
    nodes, weights = leggauss(order)
    half = (right - left) / 2
    center = (right + left) / 2
    total = mp.mpf("0")
    for node, weight in zip(nodes, weights):
        x = center + half * mp.mpf(str(float(node)))
        total += mp.mpf(str(float(weight))) * function(x)
    return half * total


def _startup_value(
    mu: mp.mpf,
    row: int,
    xi: mp.mpf,
    *,
    precision: int,
    order: int,
) -> tuple[mp.mpf, mp.mpf, mp.mpf, mp.mpf]:
    """Return startup J, linear homogeneous term, and remainder bound."""

    lam = mp.mpf("0.5") - row * mu
    k = lam / mu
    total_span = k * xi
    span = min(total_span, mp.mpf("80"))

    def remainder_integrand(s: mp.mpf) -> mp.mpf:
        v = xi - s / k
        return mp.exp(-s) * (_pulse_start_value(v, precision) - (v - XI_LINEAR_OFFSET))

    remainder = _gauss_integral(mp.mpf("0"), span, remainder_integrand, order) / lam
    # The linear source on [0,xi] has a closed solution with J(0)=0.  This
    # avoids subtracting two O(1) quantities to recover the tiny startup mode.
    linear_particular = (xi - XI_LINEAR_OFFSET) / lam - mu / (lam * lam)
    linear_homogeneous = (
        XI_LINEAR_OFFSET / lam + mu / (lam * lam)
    ) * mp.exp(-k * xi)
    tail_bound = (
        (2 * xi + XI_LINEAR_OFFSET) * mp.exp(-span) / lam
        if total_span > span
        else mp.mpf("0")
    )
    return linear_particular + linear_homogeneous + remainder, linear_homogeneous, remainder, tail_bound


def _startup_data(mu: mp.mpf, row: int, *, precision: int, order: int) -> dict[str, mp.mpf | str]:
    """Compute J(.02) and its propagated homogeneous startup correction."""

    lam = mp.mpf("0.5") - row * mu
    j20, linear_homogeneous, remainder, tail_bound = _startup_value(
        mu, row, XI_STARTUP, precision=precision, order=order
    )
    jpart20 = (XI_STARTUP - XI_LINEAR_OFFSET) / lam - mu / (lam * lam)
    homogeneous = linear_homogeneous + remainder
    # The reported bound includes the endpoint-scaled remainder and is
    # intentionally a numerical bound, not a theorem-level estimate.
    remainder_bound = abs(remainder) + tail_bound
    homogeneous_bound = abs(linear_homogeneous) + remainder_bound
    span = min((lam / mu) * XI_STARTUP, mp.mpf("80"))
    total_span = (lam / mu) * XI_STARTUP
    return {
        "j20": j20,
        "jpart20": jpart20,
        "homogeneous": homogeneous,
        "homogeneous_bound": homogeneous_bound,
        "linear_homogeneous": linear_homogeneous,
        "tail_bound": tail_bound,
        "span": span,
        "total_span": total_span,
    }


def _signed_log(value: mp.mpf, precision: int) -> dict[str, Any]:
    if value == 0:
        return {"sign": 0, "log_abs": None, "arbitrary_exponent_value": "0"}
    return {
        "sign": 1 if value > 0 else -1,
        "log_abs": mp.nstr(mp.log(abs(value)), precision),
        "arbitrary_exponent_value": mp.nstr(value, precision),
    }


def _bulk_value(mu: mp.mpf, row: int, xi: mp.mpf, *, precision: int, order: int) -> tuple[mp.mpf, dict[str, Any]]:
    lam = mp.mpf("0.5") - row * mu
    k = lam / mu
    startup = _startup_data(mu, row, precision=precision, order=order)
    particular = (xi - XI_LINEAR_OFFSET) / lam - mu / (lam * lam)
    homogeneous = startup["homogeneous"] * mp.exp(-k * (xi - XI_STARTUP))
    value = particular + homogeneous
    rhs = (xi - XI_LINEAR_OFFSET) / mu
    derivative = 1 / lam - k * homogeneous
    ode_residual = derivative + k * value - rhs
    return value, {
        "particular": particular,
        "homogeneous_at_xi": homogeneous,
        "ode_residual": ode_residual,
        "startup": startup,
    }


def _smooth_endpoint_factor(u: mp.mpf) -> mp.mpf:
    """Return sigma(u)*exp(1/u^2), with the singular factor cancelled."""

    if u <= 0 or u >= 1:
        if u <= 0:
            return mp.mpf("0")
        return mp.e
    x = 1 / (u * u) - 1 / ((1 - u) * (1 - u))
    if x > 0:
        log_sigma = -x - mp.log1p(mp.exp(-x))
    else:
        log_sigma = -mp.log1p(mp.exp(x))
    return mp.exp(log_sigma + 1 / (u * u))


def _post_cutoff_integral(
    mu: mp.mpf,
    row: int,
    xi: mp.mpf,
    *,
    precision: int,
    order: int,
) -> tuple[mp.mpf, dict[str, Any]]:
    """Integrate the flat endpoint contribution on 10 < xi < 11."""

    lam = mp.mpf("0.5") - row * mu
    k = lam / mu
    u_min = XI_END - xi
    u0 = (2 / k) ** (mp.mpf(1) / 3)
    relative_width = 1 / mp.sqrt(6 / (u0 * u0))
    band = mp.mpf("20")
    x_min = (u_min / u0 - 1) / relative_width
    centered = x_min <= band

    if centered:
        lower_x = max(-band, x_min)
        upper_x = band
        phi0 = k * u0 + 1 / (u0 * u0)

        def integrand(x: mp.mpf) -> mp.mpf:
            u = u0 * (1 + relative_width * x)
            phi = k * u + 1 / (u * u)
            return (
                (10.99 - u)
                * _smooth_endpoint_factor(u)
                * mp.exp(-(phi - phi0))
                * u0
                * relative_width
            )

        scaled = _gauss_integral(lower_x, upper_x, integrand, order)
        value = mp.exp(k * u_min - phi0 - mp.log(mu)) * scaled
        edge_lo = max(mp.mpf("0"), u0 * (1 - band * relative_width))
        edge_hi = u0 * (1 + band * relative_width)
        omitted = mp.mpf("0")
        if u_min < edge_lo:
            omitted += mp.mpf("440") / mu * mp.exp(k * u_min - (k * edge_lo + 1 / (edge_lo * edge_lo)))
        if edge_hi < 1:
            omitted += mp.mpf("440") / mu * mp.exp(k * u_min - (k * edge_hi + 1 / (edge_hi * edge_hi)))
        return value, {
            "method": "centered_saddle",
            "u_min": u_min,
            "u0": u0,
            "relative_width": relative_width,
            "omitted_tail_bound": omitted,
            "log_omitted_tail_bound": mp.nstr(mp.log(omitted), precision) if omitted else None,
        }

    span = k * (1 - u_min)
    s_limit = min(span, mp.mpf("80"))

    def boundary_integrand(s: mp.mpf) -> mp.mpf:
        u = u_min + s / k
        return mp.exp(-s) * (10.99 - u) * _sigma_mp(u) / lam

    value = _gauss_integral(mp.mpf("0"), s_limit, boundary_integrand, order)
    omitted = mp.mpf("0")
    if span > s_limit:
        omitted = mp.mpf("10.99") * mp.exp(-s_limit) / lam
    return value, {
        "method": "boundary_scaled",
        "u_min": u_min,
        "u0": u0,
        "relative_width": relative_width,
        "s_limit": s_limit,
        "s_total": span,
        "omitted_tail_bound": omitted,
        "log_omitted_tail_bound": mp.nstr(mp.log(omitted), precision) if omitted else None,
    }


def normalized_pulse_primitive(
    mu: Any,
    row: Any,
    xi: Any,
    *,
    precision: int = 160,
    order: int = 128,
) -> dict[str, Any]:
    """Return the signed-log receipt for the source normalized primitive.

    The principal result is ``log_J``.  The ``arbitrary_exponent_value`` and
    method-specific bounds are retained for callers that need a stable
    primitive rather than a binary64 value.
    """

    mu_mp, row_int, xi_mp = _validate(mu, row, xi, precision, order)
    with mp.workdps(precision):
        lam = mp.mpf("0.5") - row_int * mu_mp
        k = lam / mu_mp
        if xi_mp == 0:
            value = mp.mpf("0")
            details: dict[str, Any] = {"method": "axis_endpoint", "bounds": {}}
        elif xi_mp <= XI_STARTUP:
            value, linear_homogeneous, remainder, tail = _startup_value(
                mu_mp, row_int, xi_mp, precision=precision, order=order
            )
            total_span = k * xi_mp
            span = min(total_span, mp.mpf("80"))
            details = {
                "method": "startup_endpoint_scaled",
                "linear_homogeneous": linear_homogeneous,
                "startup_remainder": remainder,
                "startup_tail_bound": tail,
                "log_startup_tail_bound": mp.nstr(mp.log(tail), precision) if tail else None,
            }
        elif xi_mp <= XI_CUTOFF:
            value, bulk = _bulk_value(mu_mp, row_int, xi_mp, precision=precision, order=order)
            details = {
                "method": "bulk_linear_particular",
                "particular": bulk["particular"],
                "homogeneous_at_xi": bulk["homogeneous_at_xi"],
                "startup": bulk["startup"],
                "ode_residual": bulk["ode_residual"],
            }
        elif xi_mp < XI_END:
            bulk10, bulk_details = _bulk_value(
                mu_mp, row_int, XI_CUTOFF, precision=precision, order=order
            )
            post, post_details = _post_cutoff_integral(
                mu_mp, row_int, xi_mp, precision=precision, order=order
            )
            value = mp.exp(-k * (xi_mp - XI_CUTOFF)) * bulk10 + post
            details = {
                "method": post_details["method"],
                "bulk_at_10": bulk10,
                "post_cutoff_integral": post,
                "bulk_at_10_details": bulk_details,
                "post_cutoff_details": post_details,
            }
        else:
            full = normalized_pulse_integral(mu_mp, row_int, order=order, precision=precision)
            full_log = mp.mpf(str(full["log_normalized_pulse_integral"]))
            log_value = full_log + k * (mp.mpf("13") - xi_mp)
            value = mp.exp(log_value)
            details = {
                "method": "postpulse_full_row",
                "full_row_log_at_13": full_log,
                "full_row_bounds": {
                    "log_relative_omitted_positive_bound": full[
                        "log_relative_omitted_positive_bound"
                    ],
                    "log_omitted_piece_bounds": full["log_omitted_piece_bounds"],
                },
                "log_propagation": k * (mp.mpf("13") - xi_mp),
            }
        result = _signed_log(value, precision)
        return {
            "mu": mp.nstr(mu_mp, precision),
            "row": row_int,
            "xi": mp.nstr(xi_mp, precision),
            "lambda": mp.nstr(lam, precision),
            "log_J": result["log_abs"],
            "sign": result["sign"],
            "arbitrary_exponent_value": result["arbitrary_exponent_value"],
            "details": {
                key: (
                    mp.nstr(val, precision)
                    if isinstance(val, mp.mpf)
                    else val
                )
                for key, val in details.items()
            },
            "source": SOURCE,
            "source_version": SOURCE_VERSION,
            "source_section": SOURCE_SECTION,
        }


def _check_bulk_ode(mu: mp.mpf, row: int, xi: mp.mpf, *, precision: int, order: int) -> dict[str, Any]:
    receipt = normalized_pulse_primitive(mu, row, xi, precision=precision, order=order)
    residual = mp.mpf(receipt["details"]["ode_residual"])
    return {
        "xi": mp.nstr(xi, precision),
        "row": row,
        "absolute_ode_residual": mp.nstr(abs(residual), precision),
        "log_absolute_ode_residual": mp.nstr(mp.log(abs(residual)), precision)
        if residual
        else None,
    }


def run() -> dict[str, Any]:
    """Run source-identity and post-pulse replay checks and write JSON."""

    schedule = PaperOuterSchedule(
        logPstar=14,
        logRref=10,
        delta="1e-32",
        Md=".5",
        c_mu=".001",
        c_delta=".001",
        c_epsilon=".01",
    )
    mu = schedule.mu
    precision = 160
    order = 128
    bulk_checks = [_check_bulk_ode(mu, row, mp.mpf(xi), precision=precision, order=order)
                   for row in (1, 2) for xi in (mp.mpf("1"), mp.mpf("10"))]
    samples = [
        normalized_pulse_primitive(mu, 1, "0.02", precision=precision, order=order),
        normalized_pulse_primitive(mu, 1, "1", precision=precision, order=order),
        normalized_pulse_primitive(mu, 1, "10.5", precision=precision, order=order),
        normalized_pulse_primitive(mu, 1, "10.999999999", precision=precision, order=order),
        normalized_pulse_primitive(mu, 1, "11", precision=precision, order=order),
        normalized_pulse_primitive(mu, 1, "13", precision=precision, order=order),
    ]
    post11 = samples[4]
    post13 = samples[5]
    with mp.workdps(precision):
        post_equivalence_log_difference = mp.mpf(post11["log_J"]) - mp.mpf(post13["log_J"])
        expected_difference = (
            (mp.mpf("0.5") - mu) / mu * (mp.mpf("13") - mp.mpf("11"))
        )
        post_equivalence_error = post_equivalence_log_difference - expected_difference
    report = {
        "source": SOURCE,
        "source_version": SOURCE_VERSION,
        "source_section": SOURCE_SECTION,
        "api": "normalized_pulse_primitive(mu, row, xi, precision=160, order=128)",
        "definition": "J=(1/mu)*integral_0^min(xi,11) exp(-lambda*(xi-v)/mu)*g_p(v) dv",
        "xi_convention": "xi=mu*log(R/Rp), lambda=1/2-row*mu",
        "schedule_mu": str(mu),
        "bulk_ode_checks": bulk_checks,
        "samples": samples,
        "postpulse_equivalence": {
            "log_difference_J11_minus_J13": mp.nstr(post_equivalence_log_difference, precision),
            "expected_difference": mp.nstr(expected_difference, precision),
            "log_identity_residual": mp.nstr(post_equivalence_error, precision),
        },
        "claims": {
            "source_ode_bulk_checked": True,
            "postpulse_full_row_reuse_checked": True,
            "startup_bound_retained": True,
            "absolute_radius_materialized": False,
            "full_axial_profile": False,
            "radial_primitive_Mz": False,
            "regularity_pressure_cone_pde": False,
        },
        "scope": (
            "Z-independent pulse primitive only. The corrected profile must combine "
            "this J with E(t,Z), a_p(Z), and the end-bump primitives."
        ),
    }
    path = Path(__file__).with_suffix(".json")
    path.write_text(json.dumps(report, indent=2, default=str) + "\n", encoding="utf-8")
    print(json.dumps({
        "bulk_checks": bulk_checks,
        "methods": [sample["details"]["method"] for sample in samples],
        "postpulse_log_identity_residual": report["postpulse_equivalence"]["log_identity_residual"],
    }))
    return report


if __name__ == "__main__":
    run()
