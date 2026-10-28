"""Numerical port of the compact radial correction profile.

This module is deliberately a small caller-owned numerical adapter.  It does
not import the historical compact implementation and it does not manufacture
an amplitude, a control, an original bump integral, a finite ``N`` or a
bound.  The caller supplies the six ordinary-``Z`` rows of ``A`` and
``h2``, ``h3`` and ``h4``; the adapter only evaluates the displayed compact
profile formula at one fixed ``x``.

The ordinary jets use ordinary derivatives, rather than Taylor coefficients:

``(uv)[n] = sum(k=0..n) binom(n,k) u[k] v[n-k]``.

The numerical backend is a private mpmath context.  This matters for a
selected integer ``N`` with thousands of digits: it is never converted to a
binary64 float, and ``1/N`` is not allowed to become a rounded zero.  Decimal
and Fraction caller values are converted through mpmath's exact-value
converter.  The result is numerical evidence only; it is not a source
certificate or an error bound.
"""

from __future__ import annotations

from decimal import Decimal
from fractions import Fraction
from math import comb
from numbers import Integral, Rational, Real
from typing import Any, Sequence

import mpmath


COMPONENT = "compact_profile_function_port_v1"
SCHEMA = COMPONENT
ORDERS = tuple(range(6))
DEFAULT_PRECISION_DPS = 100
MAX_ORDER = 5

# These are source descriptions, not imports.  They document the exact
# formula ported by this leaf while keeping the old source modules out of the
# runtime dependency graph.
SOURCE_FORMULA_PATH = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_true_compact_C5_functions.py"
)
SOURCE_CONTROL_PATH = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_native_Rc_functional_controls.py"
)
RAW_BETA_DEFINITION = "exp(-1/(1-t^2)) for abs(t)<1, zero otherwise"
# The original native Rc route admits selected source candidates only for
# N >= 160.  This leaf intentionally keeps the displayed algebra available
# for every positive integer, but never treats that as source admission.
SOURCE_MIN_N = 160


class CompactProfileInputError(ValueError):
    """Raised when a caller supplies an invalid numerical profile input."""


class CompactProfileNumericalOverflow(ArithmeticError):
    """Raised only when a caller's data cannot be represented numerically."""


def _exact_fraction(value: Any) -> Fraction | None:
    """Return an exact real value when the caller representation exposes one.

    This is used only for the radius domain and endpoint decision.  It must
    run before conversion into the requested mpmath precision, otherwise a
    value such as ``2 + 10**-100`` can be rounded to the endpoint ``2``.
    """

    if isinstance(value, bool):
        return None
    if isinstance(value, Fraction):
        return value
    if isinstance(value, Decimal):
        try:
            return Fraction(value) if value.is_finite() else None
        except (ValueError, ZeroDivisionError, OverflowError):
            return None
    if isinstance(value, Integral):
        return Fraction(int(value), 1)
    if isinstance(value, float):
        try:
            return Fraction.from_float(value) if value == value and abs(value) != float("inf") else None
        except (ValueError, OverflowError):
            return None
    if isinstance(value, (str, bytes)):
        try:
            text = value.decode() if isinstance(value, bytes) else value
            return Fraction(text.strip())
        except (UnicodeDecodeError, ValueError, ZeroDivisionError, OverflowError):
            return None
    raw = getattr(value, "_mpf_", None)
    if raw is not None and isinstance(raw, tuple) and len(raw) == 4:
        sign, mantissa, exponent, bit_count = raw
        # mpmath uses a zero mantissa with negative bit_count sentinels for
        # NaN and infinities; only the all-zero tuple is finite zero.
        if mantissa == 0:
            if exponent == 0 and bit_count == 0:
                return Fraction(0, 1)
            return None
        try:
            magnitude = Fraction(int(mantissa), 1)
            if int(exponent) >= 0:
                magnitude *= 1 << int(exponent)
            else:
                magnitude /= 1 << (-int(exponent))
            return -magnitude if sign else magnitude
        except (TypeError, ValueError, OverflowError):
            return None
    raw_rational = getattr(value, "_mpq_", None)
    if raw_rational is not None and isinstance(raw_rational, tuple) and len(raw_rational) == 2:
        try:
            return Fraction(int(raw_rational[0]), int(raw_rational[1]))
        except (TypeError, ValueError, ZeroDivisionError, OverflowError):
            return None
    raw_interval = getattr(value, "_mpi_", None)
    if raw_interval is not None and isinstance(raw_interval, tuple) and len(raw_interval) == 2:
        lower = _exact_fraction(raw_interval[0])
        upper = _exact_fraction(raw_interval[1])
        return lower if lower is not None and lower == upper else None
    if isinstance(value, Rational):
        try:
            return Fraction(int(value.numerator), int(value.denominator))
        except (TypeError, ValueError, ZeroDivisionError, OverflowError):
            return None
    if hasattr(value, "as_integer_ratio"):
        try:
            numerator, denominator = value.as_integer_ratio()
            return Fraction(int(numerator), int(denominator))
        except (AttributeError, TypeError, ValueError, ZeroDivisionError, OverflowError):
            return None
    if isinstance(value, Real):
        try:
            numerator, denominator = value.as_integer_ratio()
            return Fraction(int(numerator), int(denominator))
        except (AttributeError, TypeError, ValueError, ZeroDivisionError, OverflowError):
            return None
    return None


def _validate_radius(context: Any, value: Any) -> tuple[Any, bool]:
    """Convert a radius only after exact pre-rounding domain validation."""

    exact = _exact_fraction(value)
    if exact is None:
        # Preserve the ordinary finite/real diagnostic for NaN, infinity and
        # unsupported complex values before reporting the stricter x contract.
        _as_mp(context, value, "x")
        raise CompactProfileInputError(
            "x must be an exactly comparable real value before precision rounding"
        )
    lower = Fraction(1, 1)
    upper = Fraction(2, 1)
    if exact < lower or exact > upper:
        raise CompactProfileInputError("x must lie in [1, 2] before precision rounding")
    radius = _as_mp(context, value, "x")
    return radius, exact == lower or exact == upper


def _new_context(precision_dps: int):
    """Return an isolated mpmath context with an explicit decimal precision."""

    if isinstance(precision_dps, bool) or not isinstance(precision_dps, Integral):
        raise CompactProfileInputError("precision_dps must be an integer")
    if precision_dps < 50:
        raise CompactProfileInputError("precision_dps must be at least 50")
    context = mpmath.mp.clone()
    context.dps = int(precision_dps)
    return context


def _as_mp(context: Any, value: Any, name: str):
    """Convert a caller number without routing it through binary64.

    ``context.convert`` is intentional here: unlike the raw ``mpf``
    constructor it accepts exact ``Fraction`` and ``Decimal`` values.
    """

    if isinstance(value, bool):
        raise CompactProfileInputError(f"{name} must be numeric, not bool")
    try:
        converted = context.convert(value)
        if hasattr(converted, "_mpi_"):
            lower, upper = converted._mpi_
            if lower != upper:
                raise TypeError("a nonzero-width interval is not a scalar real")
            converted = context.make_mpf(lower)
        if not hasattr(converted, "_mpf_"):
            raise TypeError("a complex value is not a real scalar")
    except (TypeError, ValueError, OverflowError, ArithmeticError) as exc:
        raise CompactProfileInputError(f"{name} must be a finite real number") from exc
    if not context.isfinite(converted):
        raise CompactProfileInputError(f"{name} must be finite")
    return converted


def _jet(context: Any, value: Any, name: str) -> tuple[Any, ...]:
    """Validate and convert one complete ordinary jet of order five."""

    if isinstance(value, (str, bytes)):
        raise CompactProfileInputError(f"{name} must contain six ordinary derivatives")
    try:
        values = tuple(value)
    except TypeError as exc:
        raise CompactProfileInputError(
            f"{name} must contain six ordinary derivatives"
        ) from exc
    if len(values) != 6:
        raise CompactProfileInputError(f"{name} must contain exactly six entries")
    return tuple(_as_mp(context, item, f"{name}[{index}]") for index, item in enumerate(values))


def _three_h_jets(
    context: Any,
    h2: Any,
    h3: Any,
    h4: Any,
    h: Any,
) -> tuple[tuple[Any, ...], tuple[Any, ...], tuple[Any, ...]]:
    """Accept either explicit h2/h3/h4 arguments or one ``h`` triplet."""

    if h is not None:
        if any(item is not None for item in (h2, h3, h4)):
            raise CompactProfileInputError("use h or h2/h3/h4, not both")
        try:
            supplied = tuple(h)
        except TypeError as exc:
            raise CompactProfileInputError("h must contain h2, h3 and h4 jets") from exc
        if len(supplied) != 3:
            raise CompactProfileInputError("h must contain exactly h2, h3 and h4")
        h2, h3, h4 = supplied
    elif h3 is None and h4 is None and h2 is not None:
        # A convenient positional form is
        # ``(..., h=(h2,h3,h4), J0=..., N=..., x=...)``.  A nested triplet
        # passed as the second positional argument is accepted too.
        try:
            supplied = tuple(h2)
        except TypeError:
            supplied = ()
        if len(supplied) == 3 and all(
            not isinstance(item, (str, bytes)) for item in supplied
        ):
            h2, h3, h4 = supplied

    if h2 is None or h3 is None or h4 is None:
        raise CompactProfileInputError("h2, h3 and h4 ordinary jets are required")
    return (
        _jet(context, h2, "h2"),
        _jet(context, h3, "h3"),
        _jet(context, h4, "h4"),
    )


def ordinary_product_jets(
    left: Sequence[Any],
    right: Sequence[Any],
    *,
    precision_dps: int = DEFAULT_PRECISION_DPS,
) -> tuple[Any, ...]:
    """Return ordinary derivative jets of a product through order five.

    The helper is public because the binomial ordinary-product convention is
    part of this port's caller-facing contract.  It performs no source lookup
    and accepts caller-provided values only.
    """

    context = _new_context(precision_dps)
    u = _jet(context, left, "left")
    v = _jet(context, right, "right")
    result = []
    for order in ORDERS:
        total = context.zero
        for k in range(order + 1):
            total += context.mpf(comb(order, k)) * u[k] * v[order - k]
        if not context.isfinite(total):
            raise CompactProfileNumericalOverflow(
                f"ordinary product jet overflow at order {order}"
            )
        result.append(total)
    return tuple(result)


def _support_edge_tolerance(context: Any):
    """Conservative working-precision margin around ``abs(t) == 1``."""

    # No interval certificate is claimed.  The margin only prevents a value
    # whose support side is unresolved at the selected precision from being
    # mislabeled as a source-exact support zero.
    return context.mpf(10) ** (-max(8, int(context.dps) - 8))


def _raw_beta(context: Any, t: Any) -> tuple[Any, bool, str]:
    """Evaluate raw beta and retain working-precision support diagnostics."""

    absolute_t = abs(t)
    edge_uncertain = abs(absolute_t - context.one) <= _support_edge_tolerance(context)
    if absolute_t >= context.one:
        return (
            context.zero,
            False,
            "rounded_support_edge" if edge_uncertain else "outside_open_support",
        )
    denominator = context.one - t * t
    # With a finite-precision context, a value mathematically inside the
    # interval should still have a positive denominator.  Failing closed here
    # avoids an accidental division by zero at a rounded support edge.
    if denominator <= context.zero:
        return context.zero, False, "rounded_support_edge"
    exponent = -context.one / denominator
    beta = context.exp(exponent)
    if not context.isfinite(beta):
        raise CompactProfileNumericalOverflow("raw beta could not be represented")
    return beta, True, "near_support_edge_inside" if edge_uncertain else "inside_open_support"


def _bump_rows(
    context: Any,
    *,
    J0: Any,
    x: Any,
) -> tuple[tuple[Any, ...], tuple[dict[str, Any], ...]]:
    """Return three scalar normalized bumps and their support diagnostics."""

    log_two = context.log(context.mpf(2))
    ell = log_two / context.mpf(40)
    centers = (
        log_two / context.mpf(5),
        log_two / context.mpf(2),
        (context.mpf(4) * log_two) / context.mpf(5),
    )
    log_x = context.log(x)
    bumps = []
    diagnostics = []
    for index, center in enumerate(centers):
        t = (log_x - center) / ell
        beta, in_support, reason = _raw_beta(context, t)
        support_edge_uncertain = reason in {
            "rounded_support_edge",
            "near_support_edge_inside",
        }
        if beta == context.zero:
            # Far outside support this is an exact zero for the working
            # argument; a rounded edge is kept as a separate diagnostic and
            # is never promoted to a source-exact support claim.
            bump = context.zero
            if reason == "outside_open_support":
                zero_kind = "support_zero"
            elif reason == "rounded_support_edge":
                zero_kind = "rounded_support_edge_zero"
            else:
                zero_kind = "numerically_rounded_zero"
        else:
            bump = beta / (ell * J0 * x)
            if not context.isfinite(bump):
                raise CompactProfileNumericalOverflow(
                    f"normalized bump {index} could not be represented"
                )
            zero_kind = "nonzero_numerical"
        bumps.append(bump)
        diagnostics.append(
            {
                "index": index,
                "center_fraction_of_log2": ("1/5", "1/2", "4/5")[index],
                "t": t,
                "raw_beta": beta,
                "in_open_support": in_support,
                "raw_beta_zero_reason": reason if beta == context.zero else None,
                "zero_kind": zero_kind,
                "zero_is_support_zero_at_working_precision": zero_kind == "support_zero",
                "zero_is_numerically_rounded": zero_kind == "numerically_rounded_zero",
                "support_classification_uncertain_at_working_precision": support_edge_uncertain,
                "source_support_zero_claimed": False,
            }
        )
    return tuple(bumps), tuple(diagnostics)


def _ordinary_product(context: Any, left: Sequence[Any], right: Sequence[Any]) -> tuple[Any, ...]:
    """Internal product helper using an already isolated context."""

    result = []
    for order in ORDERS:
        total = context.zero
        for k in range(order + 1):
            total += context.mpf(comb(order, k)) * left[k] * right[order - k]
        if not context.isfinite(total):
            raise CompactProfileNumericalOverflow(
                f"ordinary product jet overflow at order {order}"
            )
        result.append(total)
    return tuple(result)


def compact_radial_correction_delta_E(
    A: Any,
    h2: Any = None,
    h3: Any = None,
    h4: Any = None,
    J0: Any = None,
    N: Any = None,
    x: Any = None,
    *,
    h: Any = None,
    precision_dps: int = DEFAULT_PRECISION_DPS,
) -> dict[str, Any]:
    """Evaluate compact ``F`` and ``delta_E = A * F`` ordinary jets.

    Parameters
    ----------
    A, h2, h3, h4:
        Six-entry ordinary derivative jets supplied by the caller.  ``h`` may
        instead be used as ``(h2, h3, h4)``.
    J0:
        The caller-supplied positive original raw-bump integral.  It is never
        computed or replaced here.
    N:
        The caller-supplied positive integer.  Arbitrarily large Python
        integers are accepted and converted directly to the local mpmath
        context, never to ``float``.
    x:
        A fixed radius in ``[1, 2]``.  The original caller representation is
        checked before conversion to working precision, and exact endpoint
        status is retained separately from the rounded numerical radius.
        ``x`` is not differentiated with respect to ``Z`` by this function.
    precision_dps:
        Decimal precision of the isolated numerical context.

    Returns
    -------
    dict
        ``F`` and ``delta_E`` are six-entry mpmath tuples.  Endpoint support
        diagnostics distinguish support classification at working precision
        from numerical values; near ``abs(t)==1`` they refuse a source-exact
        support-zero claim.  The result is numerical and uncertified.
    """

    context = _new_context(precision_dps)
    amplitude = _jet(context, A, "A")
    h_jets = _three_h_jets(context, h2, h3, h4, h)
    if J0 is None:
        raise CompactProfileInputError("positive caller-supplied J0 is required")
    original_J0 = _as_mp(context, J0, "J0")
    if original_J0 <= context.zero:
        raise CompactProfileInputError("J0 must be strictly positive")
    if isinstance(N, bool) or not isinstance(N, Integral) or N <= 0:
        raise CompactProfileInputError("N must be a positive integer")
    if x is None:
        raise CompactProfileInputError("x is required")
    radius, exact_endpoint = _validate_radius(context, x)

    bumps, bump_diagnostics = _bump_rows(context, J0=original_J0, x=radius)
    inverse_N = context.one / context.mpf(N)
    # F_n = (h2_n*bump_0 + h3_n*bump_1 + h4_n*bump_2)/N.
    F = []
    for order in ORDERS:
        value = inverse_N * (
            h_jets[0][order] * bumps[0]
            + h_jets[1][order] * bumps[1]
            + h_jets[2][order] * bumps[2]
        )
        if not context.isfinite(value):
            raise CompactProfileNumericalOverflow(f"F jet overflow at order {order}")
        F.append(value)
    F = tuple(F)
    delta_E = _ordinary_product(context, amplitude, F)
    endpoint_condition = (
        "All-Z endpoint support is concluded only when x(Z) is identically 1 or "
        "2 on the whole Z domain (so every fixed-x evaluation is outside all "
        "three open bump supports).  A single point x(Z0) equal to an endpoint "
        "does not provide x derivatives and therefore does not establish an "
        "all-Z jet or history conclusion."
    )
    return {
        "schema": SCHEMA,
        "status": "numerical_returned",
        "F": F,
        "F_jet": F,
        "F_value": F[0],
        "delta_E": delta_E,
        "delta_E_jet": delta_E,
        "delta_E_jets": delta_E,
        "delta_E_value": delta_E[0],
        "bumps": bumps,
        "bump_diagnostics": bump_diagnostics,
        "A_jet": amplitude,
        "h_jets": h_jets,
        "J0": original_J0,
        "N": N,
        "x": radius,
        "x_exactly_classified_before_rounding": True,
        "x_exact_endpoint": exact_endpoint,
        "ordinary_orders": ORDERS,
        "precision_dps": int(precision_dps),
        "backend": "isolated_mpmath_context",
        "numeric_values_computed": True,
        "numerical": True,
        "numeric_evaluation": True,
        "source_inputs_provenance": "caller-supplied and not verified by this port",
        "certified": False,
        "source_admissibility": {
            "source_N_lower_bound": SOURCE_MIN_N,
            "N": N,
            "N_below_source_lower_bound": N < SOURCE_MIN_N,
            "source_admissible": None,
            "source_admissibility_claimed": False,
            "reason": (
                "This caller-owned port does not verify the original selected-N "
                "controls, source inputs or J0; N<160 is algebraically accepted "
                "but is not source-admissible."
            ),
        },
        "all_Z_endpoint_support": {
            "condition": endpoint_condition,
            "fixed_x_is_endpoint": exact_endpoint,
            "pointwise_endpoint_jet_is_zero": exact_endpoint,
            "all_Z_endpoint_support_concluded": False,
            "x_was_differentiated": False,
        },
        "h_history_conclusion": (
            "None: h2, h3 and h4 are caller-provided pointwise ordinary jets; "
            "this port infers no h history, transport, source identity or seam."
        ),
        "h_history_inferred": False,
        "all_Z_endpoint_support_condition": endpoint_condition,
    }


# Readable aliases for callers that prefer a shorter verb.  All aliases point
# to the same implementation and therefore carry the same fail-closed input
# validation and numerical-only status.
compact_radial_correction = compact_radial_correction_delta_E
evaluate_compact_radial_correction = compact_radial_correction_delta_E
evaluate_compact_profile = compact_radial_correction_delta_E
compute_compact_radial_correction = compact_radial_correction_delta_E
compute_delta_E = compact_radial_correction_delta_E
compact_profile_delta_E = compact_radial_correction_delta_E


def _fixture_jet(context: Any, coefficients: Sequence[int]) -> tuple[Any, ...]:
    """Ordinary derivatives at Z=0 of a known polynomial fixture."""

    return tuple(context.factorial(order) * context.mpf(coefficients[order]) for order in ORDERS)


def _fixture_expected(
    context: Any,
    amplitude: Sequence[Any],
    h_jets: Sequence[Sequence[Any]],
    J0: Any,
    N: int,
    x: Any,
) -> tuple[tuple[Any, ...], tuple[Any, ...]]:
    """Independent fixture arithmetic, kept separate from production helper."""

    log_two = context.log(context.mpf(2))
    ell = log_two / context.mpf(40)
    centers = (log_two / 5, log_two / 2, 4 * log_two / 5)
    log_x = context.log(x)
    bump_values = []
    for center in centers:
        t = (log_x - center) / ell
        if abs(t) >= 1:
            bump_values.append(context.zero)
        else:
            bump_values.append(context.exp(-1 / (1 - t * t)) / (ell * J0 * x))
    f = tuple(
        (sum(h_jets[row][order] * bump_values[row] for row in range(3)) / context.mpf(N))
        for order in ORDERS
    )
    delta = tuple(
        sum(context.mpf(comb(order, k)) * amplitude[k] * f[order - k] for k in range(order + 1))
        for order in ORDERS
    )
    return f, delta


def _close(context: Any, actual: Any, expected: Any) -> bool:
    """High-precision comparison for the synthetic verification fixture."""

    difference = abs(actual - expected)
    scale = max(context.one, abs(expected), abs(actual))
    return difference <= context.mpf("1e-70") * scale


def verify_synthetic_fixture() -> dict[str, Any]:
    """Verify the production port on explicit synthetic polynomials.

    The fixture is intentionally independent of the historical source and
    uses distinct coefficient rows, an interior point, and both radial
    endpoints.  It is the one function the staged runner invokes.
    """

    precision_dps = 100
    context = _new_context(precision_dps)
    amplitude = _fixture_jet(context, (2, -3, 5, -7, 11, -13))
    h2 = _fixture_jet(context, (17, -19, 23, -29, 31, -37))
    h3 = _fixture_jet(context, (41, -43, 47, -53, 59, -61))
    h4 = _fixture_jet(context, (67, -71, 73, -79, 83, -89))
    h_jets = (h2, h3, h4)
    J0 = context.mpf(11) / 13
    N = 37
    # log(x) is c_0 + ell/4, so this is strictly interior to the first bump.
    log_two = context.log(context.mpf(2))
    interior_x = context.exp(log_two / 5 + log_two / 160)
    expected_F, expected_delta = _fixture_expected(
        context, amplitude, h_jets, J0, N, interior_x
    )
    interior = compact_radial_correction_delta_E(
        amplitude,
        h2,
        h3,
        h4,
        J0,
        N,
        interior_x,
        precision_dps=precision_dps,
    )
    checks: dict[str, bool] = {
        "interior_distinct_h_coefficients": len({h[0] for h in h_jets}) == 3,
        "interior_F_binomial_formula": all(
            _close(context, actual, expected)
            for actual, expected in zip(interior["F"], expected_F)
        ),
        "interior_delta_E_binomial_formula": all(
            _close(context, actual, expected)
            for actual, expected in zip(interior["delta_E"], expected_delta)
        ),
        "interior_has_one_open_bump": sum(
            item["in_open_support"] for item in interior["bump_diagnostics"]
        )
        == 1,
        "interior_values_are_not_support_zero": all(
            item["zero_kind"] == "nonzero_numerical"
            for item in interior["bump_diagnostics"]
            if item["in_open_support"]
        ),
    }

    # Independently check each bump at its own center, where beta(0)=exp(-1).
    # Only that control may contribute because the source bumps are disjoint.
    ell = log_two / 40
    for idx, fraction in enumerate((context.mpf(1)/5, context.mpf(1)/2, context.mpf(4)/5)):
        center_x = context.exp(log_two * fraction)
        center_case = compact_radial_correction_delta_E(amplitude, h2, h3, h4, J0, N, center_x, precision_dps=precision_dps)
        expected0 = amplitude[0] * h_jets[idx][0] * context.exp(-1) / (ell * J0 * center_x * N)
        checks[f"center_{idx}_direct_value"] = _close(context, center_case["delta_E"][0], expected0)
        checks[f"center_{idx}_only_own_bump"] = sum(item["in_open_support"] for item in center_case["bump_diagnostics"]) == 1 and center_case["bump_diagnostics"][idx]["in_open_support"]
    huge_N = 10 ** 400
    large_case = compact_radial_correction_delta_E(amplitude, h2, h3, h4, J0, huge_N, interior_x, precision_dps=precision_dps)
    checks["large_N_remains_nonzero"] = large_case["delta_E"][0] != 0
    checks["large_N_ratio"] = _close(context, large_case["delta_E"][0] * context.mpf(huge_N) / (interior["delta_E"][0] * N), context.one)

    endpoint_cases = {}
    for label, endpoint_x in (("left_endpoint", context.one), ("right_endpoint", context.mpf(2))):
        endpoint = compact_radial_correction_delta_E(
            amplitude,
            h2,
            h3,
            h4,
            J0,
            N,
            endpoint_x,
            precision_dps=precision_dps,
        )
        endpoint_cases[label] = endpoint
        checks[f"{label}_all_bumps_support_zero"] = all(
            item["zero_kind"] == "support_zero"
            for item in endpoint["bump_diagnostics"]
        )
        checks[f"{label}_F_exact_zero"] = all(value == context.zero for value in endpoint["F"])
        checks[f"{label}_delta_E_exact_zero"] = all(
            value == context.zero for value in endpoint["delta_E"]
        )

    if not all(checks.values()):
        failed = ", ".join(name for name, passed in checks.items() if not passed)
        raise AssertionError(f"synthetic compact profile fixture failed: {failed}")
    return {
        "schema": SCHEMA,
        "status": "synthetic_fixture_returned",
        "fixture_tag": (
            "synthetic_explicit_interior_distinct_coefficients_known_polynomial_derivatives_"
            "and_both_endpoints_v1"
        ),
        "checks": checks,
        "checks_passed": True,
        "interior_x": interior_x,
        "endpoint_cases": {
            label: {
                "x": case["x"],
                "pointwise_endpoint_jet_is_zero": case["all_Z_endpoint_support"][
                    "pointwise_endpoint_jet_is_zero"
                ],
            }
            for label, case in endpoint_cases.items()
        },
        "source_value_computation": False,
        "source_values_computed": False,
        "synthetic_value_computation": True,
        "numeric_values_computed": True,
        "numerical": True,
        "certified": False,
        "source_admissibility": {
            "source_N_lower_bound": SOURCE_MIN_N,
            "N": N,
            "N_below_source_lower_bound": N < SOURCE_MIN_N,
            "source_admissible": False,
            "source_admissibility_claimed": False,
            "reason": "N=37 is a synthetic algebra fixture, not a selected source candidate.",
        },
        "precision_dps": precision_dps,
        "backend": "isolated_mpmath_context",
        "h_history_conclusion": (
            "No h history conclusion: fixture checks only caller-provided h2/h3/h4 jets."
        ),
    }


def verify_input_guards_fixture() -> dict[str, Any]:
    """Check exact radius guards and caller numeric representations.

    This fixture is independent of the frozen v1 fixture.  It deliberately
    keeps ``N=37`` as synthetic algebra data and does not claim source input,
    source normalization or source admissibility.
    """

    precision_dps = 50
    N = 37
    amplitude = tuple(Fraction(order + 1, 7) for order in ORDERS)
    h2 = tuple(Fraction(order + 2, 11) for order in ORDERS)
    h3 = tuple(Fraction(order + 3, 13) for order in ORDERS)
    h4 = tuple(Fraction(order + 4, 17) for order in ORDERS)
    J0 = Fraction(11, 13)
    decimal_inside = Decimal("1.15")
    # Decimal constructors preserve all written digits without applying the
    # ambient Decimal arithmetic context.
    x_above = Decimal("2." + "0" * 99 + "1")
    x_below = Decimal("1." + "9" * 100)

    checks: dict[str, bool] = {}
    try:
        compact_radial_correction_delta_E(
            amplitude,
            h2,
            h3,
            h4,
            J0,
            N,
            x_above,
            precision_dps=precision_dps,
        )
    except CompactProfileInputError:
        checks["x_2_plus_10_minus_100_rejected_before_rounding"] = True
    else:
        checks["x_2_plus_10_minus_100_rejected_before_rounding"] = False

    inside = compact_radial_correction_delta_E(
        amplitude,
        h2,
        h3,
        h4,
        J0,
        N,
        decimal_inside,
        precision_dps=precision_dps,
    )
    below = compact_radial_correction_delta_E(
        amplitude,
        h2,
        h3,
        h4,
        J0,
        N,
        x_below,
        precision_dps=precision_dps,
    )
    checks["decimal_and_fraction_inputs_supported"] = (
        inside["x_exactly_classified_before_rounding"] is True
        and inside["J0"] > 0
        and below["x_exactly_classified_before_rounding"] is True
    )
    checks["2_minus_10_minus_100_is_not_exact_endpoint"] = (
        below["x_exact_endpoint"] is False
        and below["all_Z_endpoint_support"]["fixed_x_is_endpoint"] is False
        and below["all_Z_endpoint_support"]["pointwise_endpoint_jet_is_zero"] is False
    )

    left = compact_radial_correction_delta_E(
        amplitude,
        h2,
        h3,
        h4,
        J0,
        N,
        Fraction(1, 1),
        precision_dps=precision_dps,
    )
    right = compact_radial_correction_delta_E(
        amplitude,
        h2,
        h3,
        h4,
        J0,
        N,
        Decimal("2"),
        precision_dps=precision_dps,
    )
    checks["exact_endpoints_are_identified"] = (
        left["x_exact_endpoint"] is True
        and right["x_exact_endpoint"] is True
        and left["all_Z_endpoint_support"]["pointwise_endpoint_jet_is_zero"] is True
        and right["all_Z_endpoint_support"]["pointwise_endpoint_jet_is_zero"] is True
    )
    checks["N37_is_explicitly_synthetic_not_source_admissible"] = (
        inside["source_admissibility"]["N"] == N
        and inside["source_admissibility"]["N_below_source_lower_bound"] is True
        and inside["source_admissibility"]["source_admissible"] is None
        and inside["source_admissibility"]["source_admissibility_claimed"] is False
        and inside["certified"] is False
    )

    edge_context = _new_context(precision_dps)
    near_edge_t = edge_context.one - edge_context.mpf("1e-100")
    edge_beta, edge_in_support, edge_reason = _raw_beta(edge_context, near_edge_t)
    checks["near_t_edge_is_not_source_support_zero"] = (
        edge_reason in {"rounded_support_edge", "near_support_edge_inside"}
        and not (edge_reason == "outside_open_support" and edge_beta == edge_context.zero)
    )
    checks["caller_source_provenance_not_verified"] = (
        inside["source_inputs_provenance"] == "caller-supplied and not verified by this port"
    )
    if not all(checks.values()):
        failed = ", ".join(name for name, passed in checks.items() if not passed)
        raise AssertionError(f"compact profile input guard fixture failed: {failed}")
    return {
        "schema": "compact_profile_input_guards_v1",
        "status": "input_guards_fixture_returned",
        "fixture_tag": "exact_radius_pre_rounding_decimal_fraction_support_and_near_edge_v1",
        "checks": checks,
        "checks_passed": True,
        "precision_dps": precision_dps,
        "synthetic_N": N,
        "source_admissibility": inside["source_admissibility"],
        "near_support_edge": {
            "t": near_edge_t,
            "raw_beta": edge_beta,
            "in_open_support_at_working_precision": edge_in_support,
            "reason": edge_reason,
            "source_support_zero_claimed": False,
        },
        "source_value_computation": False,
        "source_values_computed": False,
        "numeric_values_computed": True,
        "numerical": True,
        "certified": False,
        "backend": "isolated_mpmath_context",
    }


__all__ = [
    "COMPONENT",
    "SCHEMA",
    "CompactProfileInputError",
    "CompactProfileNumericalOverflow",
    "ordinary_product_jets",
    "compact_radial_correction_delta_E",
    "compact_radial_correction",
    "evaluate_compact_radial_correction",
    "evaluate_compact_profile",
    "compute_compact_radial_correction",
    "compute_delta_E",
    "compact_profile_delta_E",
    "SOURCE_MIN_N",
    "verify_synthetic_fixture",
    "verify_input_guards_fixture",
]
