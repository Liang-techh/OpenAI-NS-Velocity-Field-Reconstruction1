"""Closed special-function forms for the genuine compact bump normalizer.

The raw normalizer in this leaf is the exact symbolic quantity

    J0 = integral[-1, 1] exp(-1/(1 - t**2)) dt.

There are two independent reductions of that same integral.  The first uses
``x = t**2`` followed by ``s = x/(1-x)``.  Since

    dt on the two symmetric halves gives x**(-1/2) dx,
    1/(1-x) = 1+s, and dx = ds/(1+s)**2,

it gives

    J0 = exp(-1) integral[0, infinity]
         exp(-s) s**(-1/2) (1+s)**(-3/2) ds
       = exp(-1) sqrt(pi) U(1/2, 0, 1).

For an independent Bessel reduction, put ``t = tanh(u)``.  Then
``1/(1-t**2) = cosh(u)**2`` and ``dt = sech(u)**2 du``.  With ``x = 2u``
and ``z = 1/2`` this becomes

    J0 = exp(-1/2) integral[-infinity, infinity]
         exp(-z cosh(x))/(1+cosh(x)) dx.

For ``f(x) = sinh(x)/(1+cosh(x))``, ``f'(x) = 1/(1+cosh(x))``.  The
boundary term in integration by parts vanishes, and therefore

    integral[0, infinity] exp(-z cosh(x))/(1+cosh(x)) dx
      = z integral[0, infinity] exp(-z cosh(x)) (cosh(x)-1) dx
      = z (K_1(z) - K_0(z)).

The even extension contributes a factor two.  At the present ``z=1/2``
that factor is ``2z=1``, so

    J0 = exp(-1/2) (K_1(1/2) - K_0(1/2)).

The endpoint behavior is harmless: if ``delta = 1-|t|`` then
``1-t**2 ~ 2 delta`` and the integrand is bounded by a constant multiple of
``exp(-1/(2 delta))``.  The transformed Tricomi integrand behaves like
``s**(-1/2)`` at zero and like ``exp(-s) s**(-2)`` at infinity.  Thus both
changes of variables preserve an absolutely convergent integral.

This module evaluates only the two displayed special-function expressions.
It never calls a numerical integrator, imports an ancestor normalizer, or
traverses a source-provider graph.  The returned values are ordinary,
uncertified numerical evidence; the expression strings in the result are
the exact symbolic statements and are kept separate from those values.
Each expression is evaluated in its own cloned mpmath context, so the
caller's requested decimal precision does not mutate the process-global
mpmath context.
"""

from __future__ import annotations

import ast
from pathlib import Path
from numbers import Integral
from typing import Any

import mpmath


COMPONENT = "compact_bump_normalizer_closed_form_v1"
SCHEMA = COMPONENT
DEFAULT_PRECISION_DPS = 80
MIN_PRECISION_DPS = 20
FIXTURE_PRECISIONS = (50, 90)

EXACT_SYMBOLIC_EXPRESSIONS = {
    "integral_definition": (
        "J0 = integral_{-1}^{1} exp(-1/(1-t^2)) dt"
    ),
    "bessel": "exp(-1/2)*(K_1(1/2) - K_0(1/2))",
    "tricomi_u": "exp(-1)*sqrt(pi)*U(1/2, 0, 1)",
}


def _validate_precision(precision_dps: int) -> int:
    """Validate and normalize a caller-controlled decimal precision."""

    if isinstance(precision_dps, bool) or not isinstance(precision_dps, Integral):
        raise ValueError("precision_dps must be an integer")
    precision = int(precision_dps)
    if precision < MIN_PRECISION_DPS:
        raise ValueError(
            f"precision_dps must be at least {MIN_PRECISION_DPS} decimal digits"
        )
    return precision


def _private_context(precision_dps: int):
    """Return a cloned context without changing ``mpmath.mp.dps``."""

    precision = _validate_precision(precision_dps)
    context = mpmath.mp.clone()
    context.dps = precision
    return context


def _bessel_expression(context: Any):
    """Evaluate only the Bessel-K representation in ``context``."""

    half = context.mpf("0.5")
    return context.exp(-half) * (context.besselk(1, half) - context.besselk(0, half))


def _tricomi_u_expression(context: Any):
    """Evaluate only the Tricomi-U representation in ``context``."""

    half = context.mpf("0.5")
    return context.exp(-context.one) * context.sqrt(context.pi) * context.hyperu(
        half, context.zero, context.one
    )


def _stringify(context: Any, value: Any) -> str:
    """Render an mpmath value without routing it through binary64."""

    return context.nstr(value, n=int(context.dps), strip_zeros=False)


def evaluate_closed_forms(
    precision_dps: int = DEFAULT_PRECISION_DPS,
) -> dict[str, Any]:
    """Evaluate the two exact forms at caller-selected decimal precision.

    The Bessel and Tricomi values are calculated in separate private cloned
    contexts.  ``comparison`` contains only discrepancy diagnostics between
    those two numerical evaluations; it is not a certified error bound.
    No numerical quadrature is used.
    """

    precision = _validate_precision(precision_dps)
    bessel_context = _private_context(precision)
    tricomi_context = _private_context(precision)
    bessel_value = _bessel_expression(bessel_context)
    tricomi_value = _tricomi_u_expression(tricomi_context)

    # Compare in a third private context.  This keeps the two expression
    # evaluations independent while avoiding an accidental global-mpmath
    # operation when forming the diagnostics.
    comparison_context = _private_context(precision + 8)
    bessel_for_comparison = comparison_context.mpf(bessel_value)
    tricomi_for_comparison = comparison_context.mpf(tricomi_value)
    absolute_discrepancy = comparison_context.fabs(
        bessel_for_comparison - tricomi_for_comparison
    )
    scale = max(
        comparison_context.fabs(bessel_for_comparison),
        comparison_context.fabs(tricomi_for_comparison),
    )
    relative_discrepancy = absolute_discrepancy / scale

    return {
        "schema": SCHEMA,
        "status": "special_function_numerical_returned",
        "precision_dps": precision,
        "backend": "two_independent_private_mpmath_contexts",
        "exact_symbolic_expression": dict(EXACT_SYMBOLIC_EXPRESSIONS),
        "symbolic_identity_status": "exact_closed_form_derivation",
        "bessel_expression": {
            "exact": EXACT_SYMBOLIC_EXPRESSIONS["bessel"],
            "value": bessel_value,
            "decimal": _stringify(bessel_context, bessel_value),
        },
        "tricomi_u_expression": {
            "exact": EXACT_SYMBOLIC_EXPRESSIONS["tricomi_u"],
            "value": tricomi_value,
            "decimal": _stringify(tricomi_context, tricomi_value),
        },
        "comparison": {
            "absolute_discrepancy": absolute_discrepancy,
            "absolute_discrepancy_decimal": _stringify(
                comparison_context, absolute_discrepancy
            ),
            "relative_discrepancy": relative_discrepancy,
            "relative_discrepancy_decimal": _stringify(
                comparison_context, relative_discrepancy
            ),
            "interpretation": (
                "Agreement diagnostic between two independently evaluated "
                "special-function expressions; not a certified error bound."
            ),
        },
        "numeric_status": "uncertified_special_function_approximation",
        "numeric_values_computed": True,
        "numerical": True,
        "numeric_evaluation": True,
        "certified": False,
        "source_values_computed": False,
        "source_value_computation": False,
        "source_provider_graph_same_proof": False,
        "source_provider_graph_same_proof_claimed": False,
        "source_provenance": (
            "No source-provider module, source graph traversal, or ancestor "
            "normalizer was imported or evaluated by this leaf."
        ),
    }


# Explicit aliases make the callable easy to discover without changing its
# one implementation or introducing a second numerical route.
evaluate_j0_closed_form = evaluate_closed_forms
compute_j0_closed_form = evaluate_closed_forms


def _source_definition_check(root: Path) -> dict[str, bool]:
    """Check source syntax only; do not evaluate producer or its graph."""
    path = root / "experiments/root_st073/lei_ren_part1_paper_compliant_current_native_Rc_functional_controls.py"
    tree = ast.parse(path.read_text(encoding="utf-8"))
    functions = [n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == "exact_weights"]
    if len(functions) != 1:
        raise ValueError("source exact_weights function not unique")
    expected = {
        "raw_beta_definition": "beta=g.node('compact_raw_beta',argument=t.node,definition='exp(-1/(1-t^2)) on abs(t)<1, zero otherwise',defining_module=PREFIX+'outer_pulse_map.py',defining_module_sha256=sha(PREFIX+'outer_pulse_map.py'))",
        "integral_coordinate_and_bounds": "integrate=lambda value:integral(g,value,variable,g.constant(-1),g.one,measure='normalized original log-bump coordinate dt')",
        "normal_is_beta_integral": "normal=integrate(beta)",
    }
    actual = {ast.dump(n,include_attributes=False) for n in ast.walk(functions[0]) if isinstance(n, ast.Assign)}
    return {key: ast.dump(ast.parse(value).body[0],include_attributes=False) in actual for key,value in expected.items()}


def verify_special_function_fixture(root: Path) -> dict[str, Any]:
    """Compare both independent forms once at two precision levels.

    The refinement and cross-expression discrepancies below are diagnostics
    only.  They deliberately do not claim a rigorous numerical enclosure.
    """

    source_checks = _source_definition_check(root)
    rows = [evaluate_closed_forms(precision) for precision in FIXTURE_PRECISIONS]
    checks: dict[str, bool] = {
        **{f"source_{k}": v for k,v in source_checks.items()},
        "exact_bessel_expression_is_frozen": all(
            row["bessel_expression"]["exact"] == EXACT_SYMBOLIC_EXPRESSIONS["bessel"]
            for row in rows
        ),
        "exact_tricomi_expression_is_frozen": all(
            row["tricomi_u_expression"]["exact"]
            == EXACT_SYMBOLIC_EXPRESSIONS["tricomi_u"]
            for row in rows
        ),
        "two_distinct_special_function_routes": (
            EXACT_SYMBOLIC_EXPRESSIONS["bessel"]
            != EXACT_SYMBOLIC_EXPRESSIONS["tricomi_u"]
        ),
    }

    # A working-precision agreement threshold is useful as a diagnostic, but
    # it is intentionally labeled and returned as an uncertified indicator.
    for row in rows:
        precision = int(row["precision_dps"])
        context = _private_context(precision + 8)
        discrepancy = context.mpf(row["comparison"]["absolute_discrepancy"])
        threshold = context.power(10, -max(8, precision - 8))
        checks[f"dps_{precision}_routes_agree_at_working_precision"] = (
            discrepancy < threshold
        )

    low, high = rows
    refinement_context = _private_context(int(high["precision_dps"]) + 8)
    low_value = refinement_context.mpf(low["bessel_expression"]["value"])
    high_value = refinement_context.mpf(high["bessel_expression"]["value"])
    refinement_delta = refinement_context.fabs(high_value - low_value)
    refinement_threshold = refinement_context.power(
        10, -max(8, int(low["precision_dps"]) - 8)
    )
    checks["precision_refinement_is_finite"] = refinement_context.isfinite(
        refinement_delta
    )
    checks["precision_refinement_is_small_working_indicator"] = (
        refinement_delta < refinement_threshold
    )
    checks["no_source_provider_same_proof_claim"] = all(
        row["source_provider_graph_same_proof"] is False for row in rows
    )
    checks["all_values_are_explicitly_uncertified"] = all(
        row["certified"] is False
        and row["numeric_status"] == "uncertified_special_function_approximation"
        for row in rows
    )

    if not all(checks.values()):
        failed = ", ".join(name for name, passed in checks.items() if not passed)
        raise AssertionError(f"closed-form special-function fixture failed: {failed}")

    return {
        "schema": SCHEMA,
        "status": "special_function_fixture_returned",
        "fixture_tag": (
            "bessel_K_difference_vs_tricomi_U_two_private_contexts_two_precisions_v1"
        ),
        "checks": checks,
        "checks_passed": True,
        "precision_levels_dps": list(FIXTURE_PRECISIONS),
        "rows": rows,
        "precision_refinement": {
            "bessel_absolute_delta": refinement_delta,
            "bessel_absolute_delta_decimal": _stringify(
                refinement_context, refinement_delta
            ),
            "interpretation": (
                "Difference between two requested precisions; an approximate "
                "stability indicator, not a certified error bound."
            ),
        },
        "exact_symbolic_expression": dict(EXACT_SYMBOLIC_EXPRESSIONS),
        "numeric_status": "uncertified_special_function_approximation",
        "numeric_values_computed": True,
        "numerical": True,
        "numeric_evaluation": True,
        "certified": False,
        "source_values_computed": False,
        "source_value_computation": False,
        "source_provider_graph_same_proof": False,
        "source_provider_graph_same_proof_claimed": False,
        "source_provenance": (
            "This independent leaf has no source-provider graph identity proof; "
            "only the exact symbolic derivation and special-function diagnostics "
            "are reported."
        ),
        "backend": "two_independent_private_mpmath_contexts",
    }


__all__ = [
    "COMPONENT",
    "SCHEMA",
    "DEFAULT_PRECISION_DPS",
    "MIN_PRECISION_DPS",
    "FIXTURE_PRECISIONS",
    "EXACT_SYMBOLIC_EXPRESSIONS",
    "evaluate_closed_forms",
    "evaluate_j0_closed_form",
    "compute_j0_closed_form",
    "verify_special_function_fixture",
]
