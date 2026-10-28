"""Independent numerical compact profiles with an unmaterialized ``N``.

This leaf is a caller-owned numerical adapter for the compact profile formulas
used by the staged reconstruction.  It deliberately does not import a source
leaf, choose a value of ``N``, compute ``J0``, or infer any source/function
identity.  The profile entry point returns the scaled numerators

``X = N * delta_E = A * F`` and ``Y = N * delta_V = A * G``

at a fixed ``x``.  The factor ``1/N`` is retained as a formal descriptor, so a
caller can carry a selected (possibly non-materializable) ``N`` without this
module inventing one.  All input jets are ordinary derivative jets through
order five.  The radial rows are derivatives with respect to ``y = log(x)``
of orders zero, one and two.

The compact bump is evaluated independently in a private mpmath context.  At
and outside a support boundary the branch returns exact working-precision
zeros before evaluating any singular-looking derivative expression.  This is
only a working-precision support classification; it is not an exact boundary
certificate.  Numerical values and synthetic fixture checks do not certify a
source function, an error bound, a physical field, or a completed
reconstruction layer.
"""

from __future__ import annotations

from decimal import Decimal
from fractions import Fraction
from math import comb, isfinite
from numbers import Integral
from typing import Any, Mapping, Sequence

import mpmath


COMPONENT = "compact_scaled_profiles_v1"
SCHEMA = COMPONENT
ORDERS = tuple(range(6))
RADIAL_ORDERS = tuple(range(3))
MAX_ORDER = 5
DEFAULT_PRECISION_DPS = 100

SOURCE_FORMULA_PATH = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_true_compact_C5_functions.py"
)


class CompactScaledProfilesInputError(ValueError):
    """Raised when caller-supplied numerical data fail closed validation."""


class CompactScaledProfilesNumericalOverflow(ArithmeticError):
    """Raised when a finite caller input cannot produce a finite mpmath value."""


def _new_context(precision_dps: int):
    """Return an isolated mpmath context with explicit decimal precision."""

    if isinstance(precision_dps, bool) or not isinstance(precision_dps, Integral):
        raise CompactScaledProfilesInputError("precision_dps must be an integer")
    if precision_dps < 50:
        raise CompactScaledProfilesInputError("precision_dps must be at least 50")
    context = mpmath.mp.clone()
    context.dps = int(precision_dps)
    return context


def _as_mp(context: Any, value: Any, name: str):
    """Convert a caller value directly in the private context, never via float."""

    if isinstance(value, bool):
        raise CompactScaledProfilesInputError(f"{name} must be numeric, not bool")
    try:
        if isinstance(value, Fraction):
            # Do the rational conversion in the selected context without
            # routing the numerator or denominator through binary64.
            converted = context.mpf(value.numerator) / context.mpf(value.denominator)
        elif isinstance(value, Decimal):
            converted = context.mpf(str(value))
        else:
            converted = context.mpf(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise CompactScaledProfilesInputError(
            f"{name} must be a finite real number"
        ) from exc
    if not context.isfinite(converted):
        raise CompactScaledProfilesInputError(f"{name} must be finite")
    return converted


def _strict_interval_check(value: Any, name: str) -> None:
    """Check x's original representation before any mpmath rounding."""

    lower = 1
    upper = 2
    try:
        if isinstance(value, Fraction):
            exact = value
            if exact < lower or exact > upper:
                raise CompactScaledProfilesInputError(f"{name} must lie in [1, 2]")
            return
        if isinstance(value, Decimal):
            if not value.is_finite():
                raise CompactScaledProfilesInputError(f"{name} must be finite")
            if value < Decimal(1) or value > Decimal(2):
                raise CompactScaledProfilesInputError(f"{name} must lie in [1, 2]")
            return
        if isinstance(value, Integral) and not isinstance(value, bool):
            if value < lower or value > upper:
                raise CompactScaledProfilesInputError(f"{name} must lie in [1, 2]")
            return
        if isinstance(value, float):
            if not isfinite(value):
                raise CompactScaledProfilesInputError(f"{name} must be finite")
            if value < lower or value > upper:
                raise CompactScaledProfilesInputError(f"{name} must lie in [1, 2]")
            return
        if isinstance(value, str):
            # Decimal handles ordinary decimal and scientific notation exactly
            # enough to reject a just-outside endpoint before mpmath rounding.
            exact = Decimal(value.strip())
            if not exact.is_finite():
                raise CompactScaledProfilesInputError(f"{name} must be finite")
            if exact < Decimal(1) or exact > Decimal(2):
                raise CompactScaledProfilesInputError(f"{name} must lie in [1, 2]")
            return
        if hasattr(value, "_mpf_"):
            # An mpf already carries its own working-precision value; compare
            # that value before copying it into the private context.
            if value != value:
                raise CompactScaledProfilesInputError(f"{name} must be finite")
            if value < lower or value > upper:
                raise CompactScaledProfilesInputError(f"{name} must lie in [1, 2]")
            return
    except (ValueError, TypeError, ArithmeticError) as exc:
        if isinstance(exc, CompactScaledProfilesInputError):
            raise
        raise CompactScaledProfilesInputError(
            f"{name} must be a supported finite real number"
        ) from exc
    raise CompactScaledProfilesInputError(
        f"{name} must be a supported real number: use int, float, str, "
        "Decimal, Fraction, or mpmath mpf"
    )


def _jet(context: Any, value: Any, name: str) -> tuple[Any, ...]:
    """Validate one complete ordinary-derivative jet of order zero through five."""

    if isinstance(value, (str, bytes)):
        raise CompactScaledProfilesInputError(
            f"{name} must contain six ordinary derivatives"
        )
    try:
        values = tuple(value)
    except TypeError as exc:
        raise CompactScaledProfilesInputError(
            f"{name} must contain six ordinary derivatives"
        ) from exc
    if len(values) != len(ORDERS):
        raise CompactScaledProfilesInputError(
            f"{name} must contain exactly six entries"
        )
    return tuple(
        _as_mp(context, item, f"{name}[{order}]")
        for order, item in enumerate(values)
    )


def _coerce_h_jets(
    context: Any,
    h0: Any,
    h1: Any,
    h2: Any,
    h3: Any,
    h4: Any,
    h: Any,
) -> tuple[tuple[Any, ...], ...]:
    """Accept explicit ``h0`` through ``h4`` or one five-jet sequence."""

    explicit = (h0, h1, h2, h3, h4)
    if h is not None:
        if any(item is not None for item in explicit):
            raise CompactScaledProfilesInputError("use h or h0/h1/h2/h3/h4, not both")
        if isinstance(h, (str, bytes)):
            raise CompactScaledProfilesInputError("h must contain h0 through h4 jets")
        try:
            supplied = tuple(h)
        except TypeError as exc:
            raise CompactScaledProfilesInputError(
                "h must contain h0 through h4 jets"
            ) from exc
        if len(supplied) != 5:
            raise CompactScaledProfilesInputError(
                "h must contain exactly five ordinary jets"
            )
        explicit = supplied
    elif any(item is None for item in explicit):
        raise CompactScaledProfilesInputError(
            "h0, h1, h2, h3 and h4 ordinary jets are required"
        )
    return tuple(
        _jet(context, value, f"h{index}")
        for index, value in enumerate(explicit)
    )


def _ordinary_product(
    context: Any, left: Sequence[Any], right: Sequence[Any]
) -> tuple[Any, ...]:
    """Multiply two ordinary jets using binomial derivative coefficients."""

    result = []
    for order in ORDERS:
        total = context.zero
        for k in range(order + 1):
            total += context.mpf(comb(order, k)) * left[k] * right[order - k]
        if not context.isfinite(total):
            raise CompactScaledProfilesNumericalOverflow(
                f"ordinary product jet overflow at order {order}"
            )
        result.append(total)
    return tuple(result)


def ordinary_product_jets(
    left: Sequence[Any],
    right: Sequence[Any],
    *,
    precision_dps: int = DEFAULT_PRECISION_DPS,
) -> tuple[Any, ...]:
    """Return the ordinary-Z product jet of two caller-supplied jets."""

    context = _new_context(precision_dps)
    return _ordinary_product(
        context,
        _jet(context, left, "left"),
        _jet(context, right, "right"),
    )


def _add_jets(
    context: Any, *jets: Sequence[Any]
) -> tuple[Any, ...]:
    """Add ordinary jets without using a binary floating-point accumulator."""

    result = []
    for order in ORDERS:
        total = context.zero
        for jet in jets:
            total += jet[order]
        if not context.isfinite(total):
            raise CompactScaledProfilesNumericalOverflow(
                f"ordinary jet sum overflow at order {order}"
            )
        result.append(total)
    return tuple(result)


def _scale_jet(context: Any, jet: Sequence[Any], scalar: Any) -> tuple[Any, ...]:
    """Scale an ordinary jet by a fixed scalar."""

    result = tuple(scalar * value for value in jet)
    if not all(context.isfinite(value) for value in result):
        raise CompactScaledProfilesNumericalOverflow("ordinary jet scale overflow")
    return result


def _raw_bump(context: Any, t: Any) -> dict[str, Any]:
    """Evaluate beta and its first two t-derivatives with lazy exterior support."""

    one = context.one
    absolute_t = abs(t)
    if absolute_t >= one:
        # Do not form 1/(1-t**2) here.  This branch covers exact and rounded
        # support boundaries and is intentionally lazy with respect to all
        # singular-looking derivative factors.
        reason = (
            "support_boundary_at_working_precision"
            if absolute_t == one
            else "outside_open_support"
        )
        return {
            "t": t,
            "raw_beta": context.zero,
            "raw_beta_prime": context.zero,
            "raw_beta_second": context.zero,
            "in_open_support": False,
            "support_reason": reason,
            "lazy_exterior": True,
        }

    denominator = one - t * t
    if denominator <= context.zero:
        # This guard is defensive for a value that rounded across the open
        # interval; it retains the same lazy-zero semantics.
        return {
            "t": t,
            "raw_beta": context.zero,
            "raw_beta_prime": context.zero,
            "raw_beta_second": context.zero,
            "in_open_support": False,
            "support_reason": "rounded_support_edge",
            "lazy_exterior": True,
        }

    inverse_denominator = one / denominator
    beta = context.exp(-inverse_denominator)
    two = context.mpf(2)
    beta_prime = beta * (-two * t / (denominator * denominator))
    beta_second = beta * (
        -two / (denominator * denominator)
        - 8 * t * t / (denominator * denominator * denominator)
        + 4 * t * t
        / (denominator * denominator * denominator * denominator)
    )
    values = (beta, beta_prime, beta_second)
    if not all(context.isfinite(value) for value in values):
        raise CompactScaledProfilesNumericalOverflow(
            "compact bump derivative could not be represented"
        )
    return {
        "t": t,
        "raw_beta": beta,
        "raw_beta_prime": beta_prime,
        "raw_beta_second": beta_second,
        "in_open_support": True,
        "support_reason": "inside_open_support",
        "lazy_exterior": False,
    }


def _bump_rows(
    context: Any,
    *,
    J0: Any,
    x: Any,
) -> tuple[tuple[dict[str, Any], ...], Any, tuple[Any, ...]]:
    """Return three ``(b, b_y, b_yy)`` rows and compact geometry metadata."""

    log_two = context.log(context.mpf(2))
    ell = log_two / context.mpf(40)
    centers = (
        log_two / context.mpf(5),
        log_two / context.mpf(2),
        context.mpf(4) * log_two / context.mpf(5),
    )
    log_x = context.log(x)
    scale = context.one / (ell * J0 * x)
    rows = []
    for index, center in enumerate(centers):
        raw = _raw_bump(context, (log_x - center) / ell)
        if not raw["in_open_support"]:
            b = by = byy = context.zero
        else:
            beta = raw["raw_beta"]
            beta_prime = raw["raw_beta_prime"]
            beta_second = raw["raw_beta_second"]
            b = scale * beta
            by = scale * (beta_prime / ell - beta)
            byy = scale * (
                beta_second / (ell * ell)
                - 2 * beta_prime / ell
                + beta
            )
            if not all(context.isfinite(value) for value in (b, by, byy)):
                raise CompactScaledProfilesNumericalOverflow(
                    f"normalized bump {index} could not be represented"
                )
        rows.append(
            {
                "index": index,
                "center_fraction_of_log2": ("1/5", "1/2", "4/5")[index],
                "center": center,
                "b": b,
                "by": by,
                "byy": byy,
                "b_y": by,
                "b_yy": byy,
                **raw,
                "support_classification_scope": "working_precision_only",
                "support_boundary_near_rounding_unresolved": True,
            }
        )
    return tuple(rows), ell, centers


def _linear_bump_jets(
    context: Any,
    h_jets: Sequence[Sequence[Any]],
    bump_rows: Sequence[Mapping[str, Any]],
    *,
    field: str,
) -> tuple[Any, ...]:
    """Form one fixed-x ordinary-Z jet from selected h rows and bump scalars."""

    pieces = [
        _scale_jet(context, h_jets[index], bump_rows[index][field])
        for index in range(3)
    ]
    return _add_jets(context, *pieces)


def _validate_N_expression(expression: Any) -> str:
    """Keep a caller's symbolic N expression textual and never evaluate it."""

    if expression is None:
        return "N"
    if isinstance(expression, bool) or not isinstance(expression, str):
        raise CompactScaledProfilesInputError(
            "N_expression must be a non-empty symbolic string"
        )
    stripped = expression.strip()
    if not stripped:
        raise CompactScaledProfilesInputError(
            "N_expression must be a non-empty symbolic string"
        )
    return stripped


def _formal_inverse_N(
    power: int = 1,
    *,
    expression: str = "N",
) -> dict[str, Any]:
    """Describe an inverse-N factor without accepting or materializing N."""

    if power < 1:
        raise ValueError("formal inverse-N power must be positive")
    return {
        "symbol": "N",
        "source_expression": expression,
        "power": power,
        "formal_factor": "1/N" if power == 1 else f"1/N^{power}",
        "materialized": False,
        "selected_value": None,
    }


def _formal_unscaled_delta(label: str, expression: str) -> dict[str, Any]:
    """Describe the deferred conversion from a scaled numerator to delta."""

    return {
        "numerator": label,
        "formal_expression": f"{label}/({expression})",
        "inverse_N_factor": _formal_inverse_N(1, expression=expression),
        "numeric_jet": None,
    }


def compact_scaled_profiles(
    A: Any,
    h0: Any = None,
    h1: Any = None,
    h2: Any = None,
    h3: Any = None,
    h4: Any = None,
    J0: Any = None,
    x: Any = None,
    *,
    h: Any = None,
    precision_dps: int = DEFAULT_PRECISION_DPS,
    N_expression: Any = None,
    source_N_expression: Any = None,
) -> dict[str, Any]:
    """Return fixed-x ordinary-Z jets for ``N*delta_E`` and ``N*delta_V``.

    ``A`` and ``h0`` through ``h4`` are six-entry ordinary derivative jets.
    ``J0`` is a caller-supplied positive number and ``x`` is a fixed radius in
    ``[1, 2]``.  No numerical ``N`` argument is accepted: the returned scaled
    rows are the numerators ``X = N*delta_E`` and ``Y = N*delta_V`` and each
    deferred ``delta`` conversion carries an explicit, unmaterialized ``1/N``
    record.  ``N_expression`` (or its spelling ``source_N_expression``) may
    carry a symbolic source expression as text; it is recorded verbatim and
    never evaluated.
    """

    if N_expression is not None and source_N_expression is not None:
        raise CompactScaledProfilesInputError(
            "use N_expression or source_N_expression, not both"
        )
    formal_N_expression = _validate_N_expression(
        N_expression if N_expression is not None else source_N_expression
    )
    context = _new_context(precision_dps)
    amplitude = _jet(context, A, "A")
    h_jets = _coerce_h_jets(context, h0, h1, h2, h3, h4, h)
    if J0 is None:
        raise CompactScaledProfilesInputError("positive caller-supplied J0 is required")
    original_J0 = _as_mp(context, J0, "J0")
    if original_J0 <= context.zero:
        raise CompactScaledProfilesInputError("J0 must be strictly positive")
    if x is None:
        raise CompactScaledProfilesInputError("x is required")
    _strict_interval_check(x, "x")
    radius = _as_mp(context, x, "x")
    if radius < context.one or radius > context.mpf(2):
        raise CompactScaledProfilesInputError("x must lie in [1, 2]")

    bumps, ell, centers = _bump_rows(context, J0=original_J0, x=radius)
    radial_rows = []
    f_rows = []
    g_rows = []
    x_rows = []
    y_rows = []
    for radial_order, field in enumerate(("b", "by", "byy")):
        # F = h2*b0 + h3*b1 + h4*b2; G = h0*b0 + h1*b2.
        f_jet = _linear_bump_jets(
            context,
            h_jets[2:5],
            bumps,
            field=field,
        )
        g_jet = _add_jets(
            context,
            _scale_jet(context, h_jets[0], bumps[0][field]),
            _scale_jet(context, h_jets[1], bumps[2][field]),
        )
        scaled_delta_E = _ordinary_product(context, amplitude, f_jet)
        scaled_delta_V = _ordinary_product(context, amplitude, g_jet)
        f_rows.append(f_jet)
        g_rows.append(g_jet)
        x_rows.append(scaled_delta_E)
        y_rows.append(scaled_delta_V)
        radial_rows.append(
            {
                "radial_order": radial_order,
                "y_derivative_order": radial_order,
                "variable": "y=log(x)",
                "F": f_jet,
                "G": g_jet,
                "N_delta_E": scaled_delta_E,
                "N_delta_V": scaled_delta_V,
                "X": scaled_delta_E,
                "Y": scaled_delta_V,
                "scaled_delta_E": scaled_delta_E,
                "scaled_delta_V": scaled_delta_V,
                "delta_E_formal": _formal_unscaled_delta(
                    "N_delta_E", formal_N_expression
                ),
                "delta_V_formal": _formal_unscaled_delta(
                    "N_delta_V", formal_N_expression
                ),
                "inverse_N_factor": _formal_inverse_N(
                    1, expression=formal_N_expression
                ),
            }
        )

    formal_scaling = {
        "symbol": "N",
        "source_expression": formal_N_expression,
        "materialized": False,
        "selected_value": None,
        "source_admissibility": "not_verified",
        "source_N_condition": "N>=160 not verified; no numerical N supplied",
        "inverse_factor": "1/N",
        "scaled_numerators": {
            "N_delta_E": "A*F",
            "N_delta_V": "A*G",
        },
        "deferred_conversions": {
            "delta_E": "(N_delta_E)/N",
            "delta_V": "(N_delta_V)/N",
        },
    }
    result = {
        "schema": SCHEMA,
        "component": COMPONENT,
        "status": "numerical_returned",
        "A_jet": amplitude,
        "h_jets": h_jets,
        "J0": original_J0,
        "x": radius,
        "ell": ell,
        "centers": centers,
        "bumps": bumps,
        "bump_rows": bumps,
        "bump_diagnostics": bumps,
        "bump_values": tuple(
            (row["b"], row["by"], row["byy"]) for row in bumps
        ),
        "radial_y_rows": tuple(radial_rows),
        "radial_rows": tuple(radial_rows),
        "F_rows": tuple(f_rows),
        "G_rows": tuple(g_rows),
        "F": tuple(f_rows),
        "G": tuple(g_rows),
        "N_delta_E_rows": tuple(x_rows),
        "N_delta_V_rows": tuple(y_rows),
        "N_delta_E": tuple(x_rows),
        "N_delta_V": tuple(y_rows),
        "X": tuple(x_rows),
        "Y": tuple(y_rows),
        "X_rows": tuple(x_rows),
        "Y_rows": tuple(y_rows),
        "scaled_delta_E_rows": tuple(x_rows),
        "scaled_delta_V_rows": tuple(y_rows),
        "ordinary_Z_orders": ORDERS,
        "radial_y_orders": RADIAL_ORDERS,
        "fixed_x": True,
        "x_was_differentiated": False,
        "support_classification": "working_precision_only_noncertified",
        "formal_scaling": formal_scaling,
        "formal_inverse_N": _formal_inverse_N(
            1, expression=formal_N_expression
        ),
        "source_N_expression": formal_N_expression,
        "numeric_N_materialized": False,
        "numeric_values_computed": True,
        "numerical": True,
        "source_inputs_provenance": "caller-supplied and not verified by this port",
        "source_admissibility": "not_verified",
        "source_function_equality": "not_assessed",
        "certified": False,
        "certified_error": False,
        "physical_field_completed": False,
        "limitations": (
            "This is an independent numerical formula port only.  It does not "
            "authenticate caller source data, prove source-function equality, "
            "provide a certified error, or complete a physical field."
        ),
    }
    return result


def correction_density_decomposition(
    E: Any,
    X: Any,
    Y: Any,
    *,
    precision_dps: int = DEFAULT_PRECISION_DPS,
    N_expression: Any = None,
    source_N_expression: Any = None,
) -> dict[str, Any]:
    """Expose exact inverse-``N`` powers for the five correction densities.

    ``E`` is a caller-provided source-background ordinary jet, while ``X`` and
    ``Y`` are the scaled profile jets ``N*delta_E`` and ``N*delta_V``.  The
    function intentionally has no ``N`` argument.  Every numerator term is
    returned separately, including the quadratic terms, and no combined
    numerical density is produced.
    """

    if N_expression is not None and source_N_expression is not None:
        raise CompactScaledProfilesInputError(
            "use N_expression or source_N_expression, not both"
        )
    formal_N_expression = _validate_N_expression(
        N_expression if N_expression is not None else source_N_expression
    )
    context = _new_context(precision_dps)
    e_jet = _jet(context, E, "E")
    x_jet = _jet(context, X, "X")
    y_jet = _jet(context, Y, "Y")
    ex = _ordinary_product(context, e_jet, x_jet)
    ey = _ordinary_product(context, e_jet, y_jet)
    xy = _ordinary_product(context, x_jet, y_jet)
    xx = _ordinary_product(context, x_jet, x_jet)
    yy = _ordinary_product(context, y_jet, y_jet)
    half = context.one / context.mpf(2)
    neg_ex = _scale_jet(context, ex, -context.one)
    half_xx = _scale_jet(context, xx, half)
    neg_half_xx = _scale_jet(context, xx, -half)

    def term(
        label: str,
        inverse_power: int,
        numerator: str,
        jet: Sequence[Any],
    ) -> dict[str, Any]:
        return {
            "label": label,
            "inverse_N_power": inverse_power,
            "formal_factor": (
                "1/N" if inverse_power == 1 else f"1/N^{inverse_power}"
            ),
            "numerator": numerator,
            "jet": tuple(jet),
        }

    density_terms = {
        "m": [term("m_linear", 1, "Y", y_jet)],
        "h": [term("h_linear", 1, "X", x_jet)],
        "k": [
            term("k_linear", 1, "E*Y", ey),
            term("k_quadratic", 2, "X*Y", xy),
        ],
        "e": [
            term("e_linear", 1, "-E*X", neg_ex),
            term("e_quadratic_Y", 2, "Y^2", yy),
            term("e_quadratic_X", 2, "-X^2/2", neg_half_xx),
        ],
        "p": [
            term("p_linear", 1, "E*X", ex),
            term("p_quadratic", 2, "X^2/2", half_xx),
        ],
    }
    expressions = {
        "m": "Y/N",
        "h": "X/N",
        "k": "E*Y/N + X*Y/N^2",
        "e": "-E*X/N + (Y^2 - X^2/2)/N^2",
        "p": "E*X/N + X^2/(2*N^2)",
    }
    densities: dict[str, Any] = {}
    for name, terms in density_terms.items():
        by_power: dict[str, tuple[Any, ...]] = {}
        for power in (1, 2):
            selected = [
                item["jet"]
                for item in terms
                if item["inverse_N_power"] == power
            ]
            if selected:
                by_power[str(power)] = _add_jets(context, *selected)
        densities[name] = {
            "name": name,
            "formal_expression": expressions[name],
            "terms": tuple(terms),
            "coefficient_jets_by_inverse_N_power": by_power,
            "ordinary_Z_jets_with_N_unmaterialized": None,
            "numeric_N_materialized": False,
        }
    result = {
        "schema": SCHEMA,
        "component": COMPONENT,
        "status": "formal_N_power_decomposition_returned",
        "formal_scaling": {
            "symbol": "N",
            "source_expression": formal_N_expression,
            "materialized": False,
            "selected_value": None,
            "source_admissibility": "not_verified",
            "source_N_condition": "N>=160 not verified; no numerical N supplied",
            "allowed_inverse_powers": (1, 2),
        },
        "input_jets": {"E": e_jet, "X": x_jet, "Y": y_jet},
        "densities": densities,
        "density_names": ("m", "h", "k", "e", "p"),
        "numeric_N_materialized": False,
        "source_inputs_provenance": "caller-supplied and not verified by this port",
        "source_admissibility": "not_verified",
        "source_function_equality": "not_assessed",
        "certified": False,
        "certified_error": False,
        "physical_field_completed": False,
        "limitations": (
            "Only exact formal inverse-N powers of caller-supplied ordinary jets "
            "are exposed; no selected N or combined numerical density is made."
        ),
    }
    # Keep the structured ``densities`` namespace while making each named
    # density directly addressable for small callers.
    result.update(densities)
    return result


def _fixture_jet(context: Any, coefficients: Sequence[int]) -> tuple[Any, ...]:
    """Ordinary derivatives at zero of a known polynomial fixture."""

    return tuple(
        context.factorial(order) * context.mpf(coefficients[order])
        for order in ORDERS
    )


def _fixture_product(
    context: Any, left: Sequence[Any], right: Sequence[Any]
) -> tuple[Any, ...]:
    """Independent product arithmetic kept separate from production helpers."""

    return tuple(
        sum(
            context.mpf(comb(order, k)) * left[k] * right[order - k]
            for k in range(order + 1)
        )
        for order in ORDERS
    )


def _close(context: Any, actual: Any, expected: Any) -> bool:
    difference = abs(actual - expected)
    scale = max(context.one, abs(actual), abs(expected))
    return difference <= context.mpf("1e-70") * scale


def verify_synthetic_fixture() -> dict[str, Any]:
    """Run one self-contained formula fixture without loading a source leaf.

    The fixture checks all three bumps at their own centers, both fixed-radius
    endpoints, the closed-form center radial values, ordinary-Z product
    coefficients, and the explicit inverse-N density terms.  It uses no
    selected ``N`` and computes no source value.
    """

    precision_dps = 100
    context = _new_context(precision_dps)
    A = _fixture_jet(context, (2, -3, 5, -7, 11, -13))
    h_jets = (
        _fixture_jet(context, (17, -19, 23, -29, 31, -37)),
        _fixture_jet(context, (41, -43, 47, -53, 59, -61)),
        _fixture_jet(context, (67, -71, 73, -79, 83, -89)),
        _fixture_jet(context, (97, -101, 103, -107, 109, -113)),
        _fixture_jet(context, (127, -131, 137, -139, 149, -151)),
    )
    J0 = context.mpf(11) / context.mpf(13)
    log_two = context.log(context.mpf(2))
    ell = log_two / context.mpf(40)
    center_fractions = (
        context.mpf(1) / context.mpf(5),
        context.mpf(1) / context.mpf(2),
        context.mpf(4) / context.mpf(5),
    )
    centers_x = tuple(context.exp(log_two * fraction) for fraction in center_fractions)
    checks: dict[str, bool] = {}

    center_cases = []
    for index, center_x in enumerate(centers_x):
        case = compact_scaled_profiles(
            A,
            *h_jets,
            J0,
            center_x,
            precision_dps=precision_dps,
        )
        center_cases.append(case)
        checks[f"center_{index}_one_open_bump"] = sum(
            bool(row["in_open_support"]) for row in case["bump_rows"]
        ) == 1
        checks[f"center_{index}_own_bump_is_open"] = bool(
            case["bump_rows"][index]["in_open_support"]
        )
        beta0 = context.exp(-context.one)
        expected_b = beta0 / (ell * J0 * center_x)
        expected_by = -expected_b
        expected_byy = expected_b * (
            context.one - context.mpf(2) / (ell * ell)
        )
        center_row = case["bump_rows"][index]
        checks[f"center_{index}_beta_e_minus_one"] = _close(
            context, center_row["raw_beta"], beta0
        )
        checks[f"center_{index}_b_closed_form"] = _close(
            context, center_row["b"], expected_b
        )
        checks[f"center_{index}_by_closed_form"] = _close(
            context, center_row["by"], expected_by
        )
        checks[f"center_{index}_byy_closed_form"] = _close(
            context, center_row["byy"], expected_byy
        )

        # Validate all fixed-x F/G rows and A*F/A*G ordinary product
        # coefficients using independently assembled scalar rows.
        for radial_order, field in enumerate(("b", "by", "byy")):
            scalar_bumps = tuple(row[field] for row in case["bump_rows"])
            expected_f = tuple(
                sum(
                    h_jets[2 + bump_index][order] * scalar_bumps[bump_index]
                    for bump_index in range(3)
                )
                for order in ORDERS
            )
            expected_g = tuple(
                h_jets[0][order] * scalar_bumps[0]
                + h_jets[1][order] * scalar_bumps[2]
                for order in ORDERS
            )
            expected_x = _fixture_product(context, A, expected_f)
            expected_y = _fixture_product(context, A, expected_g)
            checks[f"center_{index}_radial_{radial_order}_F_jet"] = all(
                _close(context, actual, expected)
                for actual, expected in zip(case["F_rows"][radial_order], expected_f)
            )
            checks[f"center_{index}_radial_{radial_order}_G_jet"] = all(
                _close(context, actual, expected)
                for actual, expected in zip(case["G_rows"][radial_order], expected_g)
            )
            checks[f"center_{index}_radial_{radial_order}_A_product"] = all(
                _close(context, actual, expected)
                for actual, expected in zip(
                    case["N_delta_E_rows"][radial_order], expected_x
                )
            )
            checks[f"center_{index}_radial_{radial_order}_V_product"] = all(
                _close(context, actual, expected)
                for actual, expected in zip(
                    case["N_delta_V_rows"][radial_order], expected_y
                )
            )
        checks[f"center_{index}_N_unmaterialized"] = (
            case["formal_scaling"]["materialized"] is False
            and case["formal_scaling"]["selected_value"] is None
            and case["numeric_N_materialized"] is False
        )

    for label, endpoint_x in (
        ("left_endpoint", context.one),
        ("right_endpoint", context.mpf(2)),
    ):
        endpoint = compact_scaled_profiles(
            A,
            *h_jets,
            J0,
            endpoint_x,
            precision_dps=precision_dps,
        )
        checks[f"{label}_all_bumps_lazy_support_zero"] = all(
            row["lazy_exterior"] is True
            and row["in_open_support"] is False
            and row["b"] == context.zero
            and row["by"] == context.zero
            and row["byy"] == context.zero
            for row in endpoint["bump_rows"]
        )
        checks[f"{label}_all_F_rows_zero"] = all(
            value == context.zero
            for radial_row in endpoint["F_rows"]
            for value in radial_row
        )
        checks[f"{label}_all_G_rows_zero"] = all(
            value == context.zero
            for radial_row in endpoint["G_rows"]
            for value in radial_row
        )
        checks[f"{label}_all_scaled_rows_zero"] = all(
            value == context.zero
            for rows in (endpoint["N_delta_E_rows"], endpoint["N_delta_V_rows"])
            for radial_row in rows
            for value in radial_row
        )

    # Check that quadratic density terms remain visible as formal inverse-N
    # coefficients, without selecting a numerical N.
    density_E = _fixture_jet(context, (3, 5, -7, 11, -13, 17))
    density = correction_density_decomposition(
        density_E,
        center_cases[0]["N_delta_E_rows"][0],
        center_cases[0]["N_delta_V_rows"][0],
        precision_dps=precision_dps,
    )
    checks["density_names_are_complete"] = density["density_names"] == (
        "m",
        "h",
        "k",
        "e",
        "p",
    )
    checks["density_N_is_unmaterialized"] = (
        density["numeric_N_materialized"] is False
        and density["formal_scaling"]["materialized"] is False
    )
    checks["density_k_keeps_two_powers"] = [
        term["inverse_N_power"] for term in density["densities"]["k"]["terms"]
    ] == [1, 2]
    checks["density_e_keeps_linear_and_two_quadratics"] = [
        term["numerator"] for term in density["densities"]["e"]["terms"]
    ] == ["-E*X", "Y^2", "-X^2/2"]
    checks["density_p_keeps_quadratic"] = density["densities"]["p"]["terms"][1][
        "numerator"
    ] == "X^2/2"

    # Independent value and first-derivative checks for every formal density
    # coefficient; retain both inverse powers instead of adding a numerical N.
    X = center_cases[0]["N_delta_E_rows"][0]
    Y = center_cases[0]["N_delta_V_rows"][0]
    E0,E1 = density_E[0],density_E[1]
    X0,X1,Y0,Y1 = X[0],X[1],Y[0],Y[1]
    expected = {
        ("m", "1"): (Y0,Y1),
        ("h", "1"): (X0,X1),
        ("k", "1"): (E0*Y0,E1*Y0+E0*Y1),
        ("k", "2"): (X0*Y0,X1*Y0+X0*Y1),
        ("e", "1"): (-E0*X0,-E1*X0-E0*X1),
        ("e", "2"): (Y0*Y0-X0*X0/2,2*Y0*Y1-X0*X1),
        ("p", "1"): (E0*X0,E1*X0+E0*X1),
        ("p", "2"): (X0*X0/2,X0*X1),
    }
    for (name,power), target in expected.items():
        actual = density["densities"][name]["coefficient_jets_by_inverse_N_power"][power]
        checks[f"density_{name}_Nminus{power}_Z0_Z1"] = all(_close(context, actual[i], target[i]) for i in (0,1))
    # Input interval is tested before conversion to working precision.
    try:
        compact_scaled_profiles(A, h=h_jets, J0=Fraction(1,2), x=Decimal("2." + "0"*99 + "1"), precision_dps=50)
    except CompactScaledProfilesInputError:
        checks["outside_x_rejected_before_rounding"] = True
    else:
        checks["outside_x_rejected_before_rounding"] = False

    if not all(checks.values()):
        failed = ", ".join(name for name, passed in checks.items() if not passed)
        raise AssertionError(f"synthetic compact scaled-profile fixture failed: {failed}")
    return {
        "schema": SCHEMA,
        "component": COMPONENT,
        "status": "synthetic_fixture_returned",
        "fixture_tag": (
            "independent_scaled_profile_centers_endpoints_ordinary_Z_products_"
            "and_formal_density_powers_v1"
        ),
        "checks": checks,
        "checks_passed": True,
        "center_x": centers_x,
        "endpoint_x": (context.one, context.mpf(2)),
        "source_value_computation": False,
        "source_values_computed": False,
        "numeric_N_materialized": False,
        "numeric_values_computed": True,
        "numerical": True,
        "certified": False,
        "certified_error": False,
        "physical_field_completed": False,
        "precision_dps": precision_dps,
        "backend": "isolated_mpmath_context",
        "source_inputs_provenance": "synthetic caller data only; not authenticated",
        "source_admissibility": "not_verified",
        "source_N_condition": "N>=160 not verified; no numerical N supplied",
        "source_function_equality": "not_assessed",
    }


# Readable aliases for callers that use a verb-oriented name.
evaluate_compact_scaled_profiles = compact_scaled_profiles
compute_compact_scaled_profiles = compact_scaled_profiles
compact_profile_scaled = compact_scaled_profiles
decompose_correction_densities = correction_density_decomposition
formal_correction_density_decomposition = correction_density_decomposition
scaled_correction_density_decomposition = correction_density_decomposition


__all__ = [
    "COMPONENT",
    "SCHEMA",
    "ORDERS",
    "RADIAL_ORDERS",
    "CompactScaledProfilesInputError",
    "CompactScaledProfilesNumericalOverflow",
    "ordinary_product_jets",
    "compact_scaled_profiles",
    "evaluate_compact_scaled_profiles",
    "compute_compact_scaled_profiles",
    "compact_profile_scaled",
    "correction_density_decomposition",
    "decompose_correction_densities",
    "formal_correction_density_decomposition",
    "scaled_correction_density_decomposition",
    "verify_synthetic_fixture",
]
