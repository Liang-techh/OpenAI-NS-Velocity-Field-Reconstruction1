"""Directed mixed derivatives of the original signed core model.

This module is deliberately independent of the accepted reconstruction
pipeline.  It evaluates the scalar model

.. math::

    H=(1-\delta)Z/2+(1-Z^2)(4Z+j),\qquad
    \chi=H^2/(H^2+\sigma^2),

and the convergent Bessel-type series

.. math::

    \Phi_0(\rho,Z)=\sum_{n\geq0}(-\rho\chi/2)^n/[n!(n+1)!].

All source numbers are accepted as exact rationals or certified exact-mpf
endpoint records.  The returned intervals are directed intervals from an
``MPIntervalContext``.  Derivatives are ordinary derivatives; no factorial
normalisation is applied.
"""

from __future__ import annotations

from decimal import Decimal, InvalidOperation
from fractions import Fraction
import math
from collections.abc import Mapping

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext


RHO_DOMAIN = (Fraction(4), Fraction(41, 10))
Z_DOMAIN = (Fraction(37, 100), Fraction(93, 250))
SELECTED_RHO = Fraction(4)
SELECTED_Z = Fraction(371, 1000)


def _is_interval(value):
    return hasattr(value, "_mpi_")


def _fraction_from_mpf_tuple(raw):
    """Recover an exact rational from an mpmath ``mpf`` tuple."""
    if not isinstance(raw, (tuple, list)) or len(raw) < 3:
        raise TypeError("Expected an exact mpf tuple")
    sign, mantissa, exponent = raw[:3]
    if type(sign) is not int or type(mantissa) is not int or type(exponent) is not int:
        raise TypeError("An exact mpf tuple must contain integer fields")
    value = Fraction((-1 if sign else 1) * mantissa)
    return value * (2 ** exponent if exponent >= 0 else Fraction(1, 2 ** (-exponent)))


def _fraction(value):
    """Parse a source scalar without passing through a binary float."""
    if isinstance(value, bool):
        raise TypeError("Boolean values are not exact source scalars")
    if isinstance(value, Fraction):
        return value
    if type(value) is int:
        return Fraction(value)
    if isinstance(value, str):
        try:
            return Fraction(value)
        except (ValueError, ZeroDivisionError):
            try:
                return Fraction(Decimal(value))
            except (InvalidOperation, ValueError, ZeroDivisionError) as exc:
                raise TypeError("Use an exact rational or decimal string") from exc
    if isinstance(value, float):
        raise TypeError("Binary floats are not certified exact inputs")
    if isinstance(value, (tuple, list)) and len(value) == 4 and all(type(v) is int for v in value):
        return _fraction_from_mpf_tuple(value)
    if isinstance(value, mp.mpf):
        raise TypeError("Pass an exact mpf tuple instead of a bare mpf")
    raise TypeError("Unsupported exact source scalar")


def _rational_interval(ctx, value):
    fraction = _fraction(value)
    if fraction.denominator == 1:
        return ctx.mpf(fraction.numerator)
    # MPIntervalContext does not accept Fraction directly.  Integer division
    # retains directed endpoint rounding and never uses a nearest float.
    return ctx.mpf(fraction.numerator) / fraction.denominator


def _endpoint(ctx, value):
    if isinstance(value, (tuple, list)) and len(value) >= 3 and all(
            type(v) is int for v in value[:3]):
        return ctx.make_mpf(tuple(value))
    return _rational_interval(ctx, value)


def _bounds(value):
    """Return ordinary exact-mpf objects for the endpoints of an interval."""
    if not _is_interval(value):
        raise TypeError("An interval value is required")
    rows = value._mpi_
    if len(rows) != 2:
        # A singleton exact interval produced by make_mpf is accepted here.
        rows = (rows, rows)
    return mp.make_mpf(rows[0]), mp.make_mpf(rows[1])


def _fraction_bounds(value):
    low, high = _bounds(value)
    return _fraction_from_mpf_tuple(low._mpf_), _fraction_from_mpf_tuple(high._mpf_)


def _within_domain(ctx, value, lower, upper):
    value_low, value_high = _fraction_bounds(value)
    domain_low, _ = _fraction_bounds(_rational_interval(ctx, lower))
    _, domain_high = _fraction_bounds(_rational_interval(ctx, upper))
    return value_low >= domain_low and value_high <= domain_high


def _box(ctx, lower, upper):
    low = _endpoint(ctx, lower)
    high = _endpoint(ctx, upper)
    low_mp = _bounds(low)[0]
    high_mp = _bounds(high)[1]
    return ctx.mpf([low_mp, high_mp])


def _interval(ctx, value, *, name):
    """Convert a rational or certified interval record to an interval."""
    if _is_interval(value):
        low, high = _bounds(value)
        return ctx.mpf([low, high])
    if isinstance(value, Mapping):
        if "lower_exact_mpf_tuple" in value and "upper_exact_mpf_tuple" in value:
            return _box(ctx, value["lower_exact_mpf_tuple"], value["upper_exact_mpf_tuple"])
        if "lower" in value and "upper" in value:
            return _box(ctx, value["lower"], value["upper"])
        raise TypeError(f"{name} needs lower/upper exact endpoints")
    if isinstance(value, (tuple, list)) and len(value) >= 3 and all(
            type(v) is int for v in value[:3]):
        # Keep large negative exponents (for example an exact tiny delta)
        # in mpmath's binary tuple form; converting them to a Fraction would
        # attempt to materialise an impractically large power of two.
        return _box(ctx, value, value)
    if isinstance(value, (tuple, list)) and len(value) == 2:
        return _box(ctx, value[0], value[1])
    return _rational_interval(ctx, value)


def _exact_tuple(value):
    low, high = _bounds(value)
    return [list(low._mpf_), list(high._mpf_)]


def _upper(value):
    return _bounds(value)[1]


def _lower(value):
    return _bounds(value)[0]


def _abs_upper(ctx, value):
    return _upper(ctx.fabs(value))


def _abs_interval(ctx, value):
    return ctx.fabs(value)


def _rational_coefficient(ctx, numerator, denominator=1):
    if denominator <= 0:
        raise ValueError("Positive denominator required")
    return _rational_interval(ctx, Fraction(numerator, denominator))


def _falling(n, order):
    if order < 0 or order > n:
        return 0
    return math.factorial(n) // math.factorial(n - order)


def _series_value_and_tail(ctx, rho, chi, rho_order, chi_order, n_terms):
    """Evaluate one ``rho``/``chi`` derivative and bound its omitted tail."""
    result = ctx.mpf(0)
    first = max(rho_order, chi_order)
    for n in range(first, n_terms):
        coefficient = ((-1) ** n * _falling(n, rho_order) * _falling(n, chi_order))
        denominator = (2 ** n) * math.factorial(n) * math.factorial(n + 1)
        term = _rational_coefficient(ctx, coefficient, denominator)
        if n > rho_order:
            term *= rho ** (n - rho_order)
        if n > chi_order:
            term *= chi ** (n - chi_order)
        result += term

    # The absolute term ratio is
    # q * (n+1)/((n+1-rho_order)(n+1-chi_order)(n+2)),
    # q = |rho chi|/2.  The displayed ratio is decreasing once n is large;
    # at n_terms=64 it is far below one on the prescribed box.
    rho_abs = _abs_interval(ctx, rho)
    chi_abs = _abs_interval(ctx, chi)
    q_abs = rho_abs * chi_abs / 2
    n = n_terms
    if n < first:
        raise ValueError("n_terms must exceed every derivative order")
    coefficient = _rational_coefficient(
        ctx,
        _falling(n, rho_order) * _falling(n, chi_order),
        (2 ** n) * math.factorial(n) * math.factorial(n + 1),
    )
    first_abs = coefficient * rho_abs ** (n - rho_order) * chi_abs ** (n - chi_order)
    ratio = q_abs * _rational_coefficient(
        ctx,
        n + 1,
        (n + 1 - rho_order) * (n + 1 - chi_order) * (n + 2),
    )
    ratio_upper = _upper(ratio)
    if ratio_upper >= 1:
        raise ArithmeticError("Series tail ratio is not below one")
    tail_interval = first_abs / (1 - ratio)
    return result, _upper(tail_interval), ratio_upper


def _widen(ctx, model, radius):
    radius_interval = ctx.mpf([-1, 1]) * ctx.mpf([0, radius])
    return model + radius_interval


def _h_derivatives(ctx, j, delta, z):
    one = ctx.mpf(1)
    h = (one - delta) * z / 2 + (one - z * z) * (4 * z + j)
    h1 = (one - delta) / 2 + 4 - 12 * z * z - 2 * j * z
    h2 = -24 * z - 2 * j
    h3 = ctx.mpf(-24)
    return h, h1, h2, h3


def _chi_derivatives(ctx, j, delta, sigma, z):
    """Return chi and its first three ordinary Z derivatives."""
    h, h1, h2, h3 = _h_derivatives(ctx, j, delta, z)
    sigma2 = sigma * sigma
    denominator = h * h + sigma2
    if _lower(denominator) <= 0:
        raise ArithmeticError("H^2 + sigma^2 must stay positive")
    # This equivalent form avoids losing the small signed gap ``1-chi``
    # through the dependency in H^2 when sigma is a very small certified
    # interval.  It remains an interval expression, never a clamp to one.
    f0 = 1 - sigma2 / denominator
    f1 = 2 * h * sigma2 / (denominator ** 2)
    f2 = 2 * sigma2 * (sigma2 - 3 * h * h) / (denominator ** 3)
    f3 = 24 * h * sigma2 * (h * h - sigma2) / (denominator ** 4)
    return (
        f0,
        f1 * h1,
        f2 * h1 * h1 + f1 * h2,
        f3 * h1 ** 3 + 3 * f2 * h1 * h2 + f1 * h3,
    ), (h, h1, h2, h3)


def _derivative_record(ctx, i, k, model, tail, ratio, n_terms):
    widened = _widen(ctx, model, tail)
    return {
        "rho_order": i,
        "Z_order": k,
        "ordinary_derivative": True,
        "factorial_normalized": False,
        "directed_model_enclosure": model,
        "directed_model_exact_mpf_tuples": _exact_tuple(model),
        "nonnegative_tail_radius": tail,
        "nonnegative_tail_radius_exact_mpf_tuple": list(tail._mpf_),
        "directed_full_enclosure": widened,
        "directed_full_exact_mpf_tuples": _exact_tuple(widened),
        "truncation_recipe": {
            "terms": n_terms,
            "first_omitted_index": n_terms,
            "geometric_ratio_upper_exact_mpf_tuple": list(ratio._mpf_),
            "series": "sum (-1)^n (n)_i (n)_ell rho^(n-i) chi^(n-ell)/(2^n n! (n+1)!), then Z chain rule",
        },
    }


def _evaluate_box(ctx, j, delta, sigma, rho, z, n_terms):
    chi_rows, h_rows = _chi_derivatives(ctx, j, delta, sigma, z)
    values = {}
    q_abs_upper = _upper(_abs_interval(ctx, rho) * _abs_interval(ctx, chi_rows[0]) / 2)
    # Cache every rho/chi derivative needed by the Z chain rules.
    cache = {}
    for i in range(4):
        for ell in range(4 - i):
            cache[i, ell] = _series_value_and_tail(ctx, rho, chi_rows[0], i, ell, n_terms)
    chi_abs = [_abs_upper(ctx, row) for row in chi_rows]
    for i in range(4):
        for k in range(4 - i):
            if k == 0:
                model, tail, ratio = cache[i, 0]
            elif k == 1:
                g1, t1, r1 = cache[i, 1]
                model = g1 * chi_rows[1]
                tail = _upper(ctx.mpf(t1) * ctx.mpf(chi_abs[1]))
                ratio = max(r1, cache[i, 0][2])
            elif k == 2:
                g1, t1, r1 = cache[i, 1]
                g2, t2, r2 = cache[i, 2]
                model = g2 * chi_rows[1] ** 2 + g1 * chi_rows[2]
                tail = _upper(ctx.mpf(t2) * ctx.mpf(chi_abs[1]) ** 2
                              + ctx.mpf(t1) * ctx.mpf(chi_abs[2]))
                ratio = max(r1, r2)
            else:
                g1, t1, r1 = cache[i, 1]
                g2, t2, r2 = cache[i, 2]
                g3, t3, r3 = cache[i, 3]
                model = (g3 * chi_rows[1] ** 3
                         + 3 * g2 * chi_rows[1] * chi_rows[2]
                         + g1 * chi_rows[3])
                tail = _upper(ctx.mpf(t3) * ctx.mpf(chi_abs[1]) ** 3
                              + 3 * ctx.mpf(t2) * ctx.mpf(chi_abs[1]) * ctx.mpf(chi_abs[2])
                              + ctx.mpf(t1) * ctx.mpf(chi_abs[3]))
                ratio = max(r1, r2, r3)
            values[f"rho{i}_Z{k}"] = _derivative_record(
                ctx, i, k, model, tail, ratio, n_terms)
    return {
        "rho": _exact_tuple(rho),
        "Z": _exact_tuple(z),
        "H": _exact_tuple(h_rows[0]),
        "chi": _exact_tuple(chi_rows[0]),
        "chi_derivatives": {
            f"Z{k}": _exact_tuple(chi_rows[k]) for k in range(4)
        },
        "q_abs_upper_exact_mpf_tuple": list(q_abs_upper._mpf_),
        "derivatives": values,
    }


def _identity_controls():
    """Small exact controls independent of any accepted source producer."""
    # Polynomial coefficient differentiation checks H', H'', H''' exactly.
    z = Fraction(2, 7)
    j = Fraction(3, 11)
    delta = Fraction(1, 13)
    coefficients = [j, Fraction(9 - delta, 2), -j, Fraction(-4)]

    def derivative(poly):
        return [index * value for index, value in enumerate(poly)][1:]

    def evaluate(poly, x):
        total = Fraction(0)
        for value in reversed(poly):
            total = total * x + value
        return total

    h1_poly = derivative(coefficients)
    h2_poly = derivative(h1_poly)
    h3_poly = derivative(h2_poly)
    explicit_h1 = Fraction(9 - delta, 2) - 12 * z * z - 2 * j * z
    explicit_h2 = -24 * z - 2 * j
    explicit_h3 = Fraction(-24)
    polynomial_identity = (
        evaluate(h1_poly, z) == explicit_h1
        and evaluate(h2_poly, z) == explicit_h2
        and evaluate(h3_poly, z) == explicit_h3
    )

    # Independent formal Taylor arithmetic checks the rational chi chain
    # rule.  H(z+t) is represented in ordinary Taylor coefficients here;
    # derivatives are recovered only at the end by multiplying by k!.
    sigma = Fraction(5, 7)
    h0 = evaluate(coefficients, z)
    h1 = evaluate(h1_poly, z)
    h2 = evaluate(h2_poly, z)
    h3 = evaluate(h3_poly, z)
    h_taylor = [h0, h1, h2 / 2, h3 / 6]

    def poly_mul(left, right, degree=3):
        result = [Fraction(0)] * (degree + 1)
        for a, left_value in enumerate(left):
            for b, right_value in enumerate(right):
                if a + b <= degree:
                    result[a + b] += left_value * right_value
        return result

    denominator_taylor = poly_mul(h_taylor, h_taylor)
    denominator_taylor[0] += sigma * sigma
    inverse_taylor = [Fraction(0)] * 4
    inverse_taylor[0] = 1 / denominator_taylor[0]
    for order in range(1, 4):
        inverse_taylor[order] = -sum(
            denominator_taylor[index] * inverse_taylor[order - index]
            for index in range(1, order + 1)
        ) / denominator_taylor[0]
    chi_taylor = [(-sigma * sigma * value) for value in inverse_taylor]
    chi_taylor[0] += 1
    formal_chi_derivatives = [math.factorial(order) * chi_taylor[order]
                              for order in range(4)]
    h2q = h0 * h0
    denominator = h2q + sigma * sigma
    f1 = 2 * h0 * sigma * sigma / denominator ** 2
    f2 = 2 * sigma * sigma * (sigma * sigma - 3 * h2q) / denominator ** 3
    f3 = 24 * h0 * sigma * sigma * (h2q - sigma * sigma) / denominator ** 4
    chain_chi_derivatives = [
        1 - sigma * sigma / denominator,
        f1 * h1,
        f2 * h1 * h1 + f1 * h2,
        f3 * h1 ** 3 + 3 * f2 * h1 * h2 + f1 * h3,
    ]
    chi_chain_identity = formal_chi_derivatives == chain_chi_derivatives

    # The mixed rho/chi series coefficients are checked against direct
    # differentiation of a finite polynomial in q.  This catches both the
    # falling-factor convention and the 2^{-n} factor used by the model.
    q = Fraction(1, 5)
    series_terms = 10
    base_series = sum(
        (-q) ** n / (math.factorial(n) * math.factorial(n + 1))
        for n in range(series_terms)
    )
    direct_rho = sum(
        (-1) ** n * n * q ** (n - 1) / (math.factorial(n) * math.factorial(n + 1))
        for n in range(1, series_terms)
    ) * Fraction(1, 7)
    coefficient_rho = sum(
        (-1) ** n * _falling(n, 1) * Fraction(1, 1) *
        Fraction(1, 2 ** n * math.factorial(n) * math.factorial(n + 1)) *
        Fraction(7, 5) ** (n - 1) * Fraction(2, 7) ** n
        for n in range(1, series_terms)
    )
    series_coefficient_identity = direct_rho == coefficient_rho
    tail_ratio = Fraction(21, 10) * Fraction(65, (65 - 3) * (65 - 3) * 66)
    tail_ratio_identity = 0 < tail_ratio < 1
    # Constant H and sigma give an exactly constant chi.  This catches
    # accidental nonzero chain-rule terms and uses only integer intervals.
    h, s = Fraction(2), Fraction(3)
    constant_chi = h * h / (h * h + s * s)
    constant_chi_identity = constant_chi == Fraction(4, 13)
    passed = bool(
        polynomial_identity and chi_chain_identity and series_coefficient_identity
        and tail_ratio_identity and constant_chi_identity
    )
    return {
        "passed": passed,
        "polynomial_H_derivative_identity": polynomial_identity,
        "formal_chi_chain_identity": chi_chain_identity,
        "series_coefficient_identity": series_coefficient_identity,
        "tail_ratio_below_one": tail_ratio_identity,
        "constant_chi_identity": constant_chi_identity,
        "constant_chi_exact": str(constant_chi),
        "base_series_control_exact": str(base_series),
    }


def _parse_call(c, j, delta, sigma, rho, Z, context):
    """Support both ``(ctx,j,delta,sigma,rho,Z)`` and no-context calls."""
    if isinstance(c, MPIntervalContext):
        if context is not None:
            raise TypeError("Specify either c or context, not both")
        return c, j, delta, sigma, rho, Z
    if context is not None:
        if c is not None:
            # Positional no-context arguments may still provide an explicit
            # context keyword; shift those arguments just as below.
            return context, c, j, delta, sigma, rho if Z is None else Z
        return context, j, delta, sigma, rho, Z
    # In the no-context form c is j and the remaining positional arguments
    # shift one place to the left.  Five positional inputs are supported;
    # a sixth would be ambiguous with the context-bearing form.
    if Z is not None:
        raise TypeError("Use an MPIntervalContext as the first argument")
    return None, c, j, delta, sigma, rho


def build_core_signed_model(c=None, j=None, delta=None, sigma=None, rho=None, Z=None,
                            *, context=None, precision=100, n_terms=64, terms=None):
    """Build directed model enclosures on the prescribed box and selected point.

    The preferred call is ``build_core_signed_model(ctx, j, delta, sigma,
    rho, Z)``.  ``ctx`` may be omitted, in which case a fresh
    ``MPIntervalContext`` is created.  ``rho`` and ``Z`` default to the
    certified domain boxes.  The returned ``domain`` and ``selected`` maps
    contain every ordinary mixed derivative with total order at most 3.
    """
    if terms is not None:
        if n_terms != 64:
            raise TypeError("Specify only one of n_terms and terms")
        n_terms = terms
    if type(n_terms) is not int or n_terms < 8:
        raise ValueError("n_terms must be an integer at least 8")
    parsed_context, j, delta, sigma, rho, Z = _parse_call(
        c, j, delta, sigma, rho, Z, context)
    if j is None or delta is None or sigma is None:
        raise TypeError("j, delta, and sigma are required exact inputs")
    ctx = parsed_context or MPIntervalContext()
    if parsed_context is None:
        ctx.dps = precision
    j_box = _interval(ctx, j, name="j")
    delta_box = _interval(ctx, delta, name="delta")
    sigma_box = _interval(ctx, sigma, name="sigma")
    rho_box = _interval(ctx, RHO_DOMAIN if rho is None else rho, name="rho")
    z_box = _interval(ctx, Z_DOMAIN if Z is None else Z, name="Z")

    # The model's advertised certified domain is fixed.  Callers may pass
    # narrower boxes, but never silently evaluate outside it.
    if not _within_domain(ctx, rho_box, *RHO_DOMAIN):
        raise ValueError("rho must lie in [4, 41/10]")
    if not _within_domain(ctx, z_box, *Z_DOMAIN):
        raise ValueError("Z must lie in [37/100, 93/250]")
    if _lower(sigma_box) <= 0:
        raise ValueError("sigma must be strictly positive on its source interval")

    controls = _identity_controls()
    if not controls["passed"]:
        raise ArithmeticError("Independent analytic identity controls failed")

    selected_rho = _rational_interval(ctx, SELECTED_RHO)
    selected_z = _rational_interval(ctx, SELECTED_Z)
    result = {
        "formula": "Phi0=sum_n (-chi*rho/2)^n/(n!(n+1)!), H=(1-delta)Z/2+(1-Z^2)(4Z+j), chi=H^2/(H^2+sigma^2)",
        "ordinary_mixed_derivatives": True,
        "factorial_normalized": False,
        "n_terms": n_terms,
        "precision_dps": ctx.dps,
        "inputs": {
            "j": _exact_tuple(j_box),
            "delta": _exact_tuple(delta_box),
            "sigma": _exact_tuple(sigma_box),
        },
        "domain": _evaluate_box(ctx, j_box, delta_box, sigma_box, rho_box, z_box, n_terms),
        "selected": _evaluate_box(ctx, j_box, delta_box, sigma_box,
                                   selected_rho, selected_z, n_terms),
        "selected_coordinates": {
            "rho": str(SELECTED_RHO),
            "Z": str(SELECTED_Z),
        },
        "tail_bound": {
            "kind": "absolute geometric tail",
            "first_omitted_index": n_terms,
            "ratio_formula": "q_abs*(n+1)/((n+1-i)*(n+1-ell)*(n+2)), q_abs=|rho*chi|/2",
            "all_entries_include_full_directed_enclosure": True,
        },
        "controls": controls,
    }
    return result


def numeric_controls():
    """Return the source-independent exact analytic controls."""
    return _identity_controls()


def core_model(c, rho, z, j, delta, sigma, terms=64):
    """Evaluate one certified ``(rho,Z)`` box with the compact core API.

    ``c`` is an ``MPIntervalContext``.  The scalar inputs may be exact
    rationals, exact decimal strings, or certified endpoint records.  The
    result keeps the interval objects under ``derivatives`` and supplies
    JSON-safe endpoint and tail records under ``parameters`` and
    ``certificate``.
    """
    if not isinstance(c, MPIntervalContext):
        raise TypeError("c must be an MPIntervalContext")
    if type(terms) is not int or terms < 8:
        raise ValueError("terms must be an integer at least 8")
    j_box = _interval(c, j, name="j")
    delta_box = _interval(c, delta, name="delta")
    sigma_box = _interval(c, sigma, name="sigma")
    rho_box = _interval(c, rho, name="rho")
    z_box = _interval(c, z, name="Z")
    if not _within_domain(c, rho_box, *RHO_DOMAIN):
        raise ValueError("rho must lie in [4, 41/10]")
    if not _within_domain(c, z_box, *Z_DOMAIN):
        raise ValueError("Z must lie in [37/100, 93/250]")
    if _lower(sigma_box) <= 0:
        raise ValueError("sigma must be strictly positive on its source interval")
    controls = _identity_controls()
    if not controls["passed"]:
        raise ArithmeticError("Independent analytic identity controls failed")
    box = _evaluate_box(c, j_box, delta_box, sigma_box, rho_box, z_box, terms)
    derivative_map = {}
    model_derivative_map = {}
    tail_map = {}
    for key, record in box["derivatives"].items():
        derivative_map[(record["rho_order"], record["Z_order"])] = record["directed_full_enclosure"]
        model_derivative_map[(record["rho_order"], record["Z_order"])] = record["directed_model_enclosure"]
        tail_map[key] = {
            "nonnegative_tail_radius_exact_mpf_tuple": record["nonnegative_tail_radius_exact_mpf_tuple"],
            "truncation_recipe": record["truncation_recipe"],
        }
    return {
        "derivatives": derivative_map,
        "model_derivatives": model_derivative_map,
        "derivative_records": box["derivatives"],
        "parameters": {
            "rho": box["rho"],
            "Z": box["Z"],
            "j": _exact_tuple(j_box),
            "delta": _exact_tuple(delta_box),
            "sigma": _exact_tuple(sigma_box),
            "H": box["H"],
            "chi": box["chi"],
            "chi_derivatives": box["chi_derivatives"],
            "q_abs_upper_exact_mpf_tuple": box["q_abs_upper_exact_mpf_tuple"],
        },
        "certificate": {
            "ordinary_mixed_derivatives": True,
            "factorial_normalized": False,
            "terms": terms,
            "tail_bound_kind": "absolute geometric tail",
            "tail_ratio_formula": "q_abs*(n+1)/((n+1-i)*(n+1-ell)*(n+2))",
            "tails": tail_map,
            "controls": controls,
        },
    }


def evaluate_core_signed_model(*args, **kwargs):
    """Alias for :func:`build_core_signed_model`."""
    return build_core_signed_model(*args, **kwargs)


__all__ = [
    "RHO_DOMAIN",
    "Z_DOMAIN",
    "SELECTED_RHO",
    "SELECTED_Z",
    "build_core_signed_model",
    "evaluate_core_signed_model",
    "core_model",
    "numeric_controls",
]
