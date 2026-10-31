"""Synthetic fixture for resolving one original Q4 endpoint source port.

Only small JSON-like dictionaries are constructed here.  The fixture does not
decode an archive, inspect a real graph, or evaluate a field.  Every case
builds a new graph, source port, and environment before calling the supplied
resolver.
"""
from __future__ import annotations

from collections.abc import Callable, Mapping
from fractions import Fraction
from typing import Any


SCHEMA = "original_q4_endpoint_resolver_fixture_v1"
ARCHIVE_SHA = "a" * 64
PARAMETER_EXACT_MPF_TUPLE = [0, 3, -6, 2]
COORDINATE_EXACT_MPF_TUPLE = [0, 3, -7, 2]
ENVIRONMENT_Z0 = (0, 1, -3, 1)
EXPECTED_POSITIVE_VALUE = Fraction(19, 64)


def _add(nodes: list[dict[str, Any]], row: dict[str, Any]) -> int:
    nodes.append(row)
    return len(nodes) - 1


def _base_case(*, active_invalid: bool) -> tuple[
    dict[str, Any], dict[str, Any], dict[int, tuple[int, int, int, int]], dict[str, int]
]:
    """Build one fresh graph with the supported synthetic operations."""

    # The first two rows are deliberately the exact zero and exact one nodes.
    nodes: list[dict[str, Any]] = [
        {"operation": "exact_rational", "numerator": 0, "denominator": 1},
        {"operation": "exact_rational", "numerator": 1, "denominator": 1},
    ]
    parameter = _add(
        nodes,
        {
            "operation": "current_original_source_parameter",
            "name": "sc",
            "Z_independent": True,
            "source_family": {"fixture": "synthetic_q4_endpoint_resolver_v1"},
        },
    )
    coordinate = _add(
        nodes,
        {
            "operation": "bound_variable",
            "name": "first_micro_coordinate",
            "Z_independent": True,
        },
    )
    z0 = _add(
        nodes,
        {
            "operation": "symbol",
            "name": "Z",
            "order": 0,
            "ordinary_Z_order": 0,
        },
    )
    half = _add(
        nodes,
        {"operation": "exact_rational", "numerator": 1, "denominator": 2},
    )
    parameter_half = _add(
        nodes,
        {"operation": "product", "arguments": [parameter, half]},
    )
    upper = _add(
        nodes,
        {"operation": "exact_rational", "numerator": 1, "denominator": 4},
    )
    comparison = _add(
        nodes,
        {
            "operation": "exact_real_comparison",
            "operator": "<=",
            "left": coordinate,
            "right": upper,
        },
    )
    predicate = _add(
        nodes,
        {"operation": "logical_or", "arguments": [comparison, 0]},
    )
    active_body = _add(
        nodes,
        (
            {
                "operation": "synthetic_unsupported_active_body",
                "label": "must_not_be_read_when_flat",
            }
            if active_invalid
            else {"operation": "exact_rational", "numerator": 99, "denominator": 1}
        ),
    )
    lazy = _add(
        nodes,
        {
            "operation": "original_lazy_flat_branch",
            "flat_predicate": predicate,
            "flat_value": 0,
            "active_body": active_body,
            "ordinary_slow_Z_derivative_order": 0,
            "active_body_not_evaluated_on_flat_branch": True,
        },
    )
    reciprocal = _add(
        nodes,
        {"operation": "positive_quotient", "numerator": 1, "denominator": half},
    )
    subtotal = _add(
        nodes,
        {"operation": "sum", "arguments": [parameter_half, z0, lazy]},
    )
    root = _add(nodes, {"operation": "product", "arguments": [subtotal, reciprocal]})

    graph = {"exact_graph_nodes": nodes}
    port = {
        "schema": "q4_physical_Uy_endpoint_source_port_v1",
        "archive_sha256": ARCHIVE_SHA,
        "expression_node": root,
        "ordinary_Z_order": 0,
        "parameter_node": parameter,
        "parameter_exact_mpf_tuple": list(PARAMETER_EXACT_MPF_TUPLE),
        "coordinate_node": coordinate,
        "coordinate_exact_mpf_tuple": list(COORDINATE_EXACT_MPF_TUPLE),
    }
    # Environment bindings are node-id bindings, never values for arithmetic
    # nodes.  Z is a symbolic node whose exact dyadic is supplied here.
    environment = {z0: ENVIRONMENT_Z0}
    refs = {
        "parameter": parameter,
        "coordinate": coordinate,
        "z0": z0,
        "comparison": comparison,
        "lazy": lazy,
        "root": root,
    }
    return graph, port, environment, refs


def _case(
    name: str,
) -> tuple[dict[str, Any], dict[str, Any], dict[int, tuple[int, int, int, int]]]:
    """Construct a named case without sharing mutable input dictionaries."""

    graph, port, environment, refs = _base_case(
        active_invalid=name == "lazy_flat_active_invalid"
    )
    nodes = graph["exact_graph_nodes"]

    if name == "unknown_operation":
        port["expression_node"] = _add(
            nodes, {"operation": "synthetic_unknown_operation"}
        )
    elif name == "missing_argument":
        port["expression_node"] = _add(
            nodes, {"operation": "positive_quotient", "numerator": 1}
        )
    elif name == "missing_environment":
        environment = {}
    elif name == "env_root_override":
        # A resolver must reject an environment attempting to overwrite the
        # arithmetic root rather than treating it as a field value.
        environment = {refs["root"]: ENVIRONMENT_Z0}
    elif name == "unproven_denominator":
        denominator = _add(
            nodes,
            {
                "operation": "current_original_source_parameter",
                "name": "unproved_denominator",
                "Z_independent": False,
            },
        )
        port["expression_node"] = _add(
            nodes,
            {
                "operation": "positive_quotient",
                "numerator": 1,
                "denominator": denominator,
            },
        )
    elif name == "unsupported_condition":
        nodes[refs["comparison"]]["operator"] = "=="
    elif name == "z_order_limit":
        z3 = _add(
            nodes,
            {
                "operation": "symbol",
                "name": "Z",
                "order": 3,
                "ordinary_Z_order": 3,
            },
        )
        environment = {z3: ENVIRONMENT_Z0}
        port["ordinary_Z_order"] = 3
        port["expression_node"] = z3
    elif name == "bad_sc":
        port["parameter_exact_mpf_tuple"] = [0, 6, -7, 3]
    elif name == "bad_coordinate":
        port["coordinate_exact_mpf_tuple"] = list(PARAMETER_EXACT_MPF_TUPLE)
    elif name not in ("positive", "lazy_flat_active_invalid"):
        raise ValueError("unknown fixture case: " + name)

    return graph, port, environment


def _scope_flag_matches(result: Any, expected: str) -> bool:
    scope = result.get("scope") if isinstance(result, Mapping) else None
    return (
        isinstance(scope, Mapping)
        and scope.get("numeric_field_value_computed") is (expected == "resolved_exact")
    )


def _exact_expression_shape(result: Any) -> bool:
    """Check only the exact-expression contract, without evaluating it."""

    if not isinstance(result, Mapping):
        return False
    expression = result.get("expression")
    if not isinstance(expression, Mapping):
        return False
    if "approx" in expression or "bound" in expression:
        return False
    if expression.get("kind") == "exact_binary_rational":
        return (
            type(expression.get("numerator")) is int
            and type(expression.get("denominator")) is int
            and type(expression.get("exponent")) is int
        )
    # Symbolic exact expressions are allowed; they still need a named kind.
    return isinstance(expression.get("kind"), str) and bool(expression.get("kind"))


def _exact_expression_equals(result: Any, expected: Fraction) -> bool:
    """Check a small synthetic exact value without evaluating any field."""

    if not isinstance(result, Mapping):
        return False
    expression = result.get("expression")
    if not isinstance(expression, Mapping):
        return False
    if expression.get("kind") != "exact_binary_rational":
        return False
    numerator = expression.get("numerator")
    denominator = expression.get("denominator")
    exponent = expression.get("exponent")
    if (
        type(numerator) is not int
        or type(denominator) is not int
        or type(exponent) is not int
        or denominator == 0
    ):
        return False
    try:
        value = Fraction(numerator, denominator) * (Fraction(2, 1) ** exponent)
    except (TypeError, ValueError, ZeroDivisionError):
        return False
    return value == expected


def _diagnostic(case: str, result: Any, expected: str) -> dict[str, Any]:
    if not isinstance(result, Mapping):
        return {"case": case, "expected": expected, "status": None, "passed": False}
    return {
        "case": case,
        "expected": expected,
        "status": result.get("status"),
        "passed": result.get("status") == expected,
    }


def run_fixture(resolve: Callable[..., Any]) -> dict[str, Any]:
    """Run the synthetic exact endpoint resolver admission matrix."""

    checks: dict[str, bool] = {}
    diagnostics: list[dict[str, Any]] = []
    if not callable(resolve):
        return {
            "schema": SCHEMA,
            "passed": False,
            "checks": {"resolve_callable": False},
            "diagnostics": [{"case": "input", "passed": False}],
            "synthetic_only": True,
            "scope": {"numeric_field_value_computed": False},
        }
    checks["resolve_callable"] = True

    cases = (
        ("positive", "resolved_exact"),
        ("lazy_flat_active_invalid", "resolved_exact"),
        ("unknown_operation", "unresolved"),
        ("missing_argument", "mismatch"),
        ("missing_environment", "unresolved"),
        ("env_root_override", "mismatch"),
        ("unproven_denominator", "unresolved"),
        ("unsupported_condition", "unresolved"),
        ("z_order_limit", "mismatch"),
        ("bad_sc", "mismatch"),
        ("bad_coordinate", "mismatch"),
    )

    for case, expected in cases:
        graph, port, environment = _case(case)
        try:
            result = resolve(graph, port, environment)
            checks[case + "_status"] = (
                isinstance(result, Mapping) and result.get("status") == expected
            )
            checks[case + "_scope"] = _scope_flag_matches(result, expected)
            if expected == "resolved_exact":
                checks[case + "_exact_expression"] = _exact_expression_shape(result)
                checks[case + "_exact_value_19_over_64"] = _exact_expression_equals(
                    result, EXPECTED_POSITIVE_VALUE
                )
            diagnostics.append(_diagnostic(case, result, expected))
        except Exception as exc:
            checks[case + "_status"] = False
            checks[case + "_scope"] = False
            if expected == "resolved_exact":
                checks[case + "_exact_expression"] = False
                checks[case + "_exact_value_19_over_64"] = False
            diagnostics.append(
                {
                    "case": case,
                    "expected": expected,
                    "status": None,
                    "passed": False,
                    "exception": type(exc).__name__,
                }
            )

    return {
        "schema": SCHEMA,
        "passed": bool(checks) and all(value is True for value in checks.values()),
        "checks": checks,
        "diagnostics": diagnostics,
        "test_count": len(checks),
        "synthetic_only": True,
        "source_code_executed": False,
        "source_graph_walked": False,
        "source_values_evaluated": False,
        "scope": {"numeric_field_value_computed": False},
    }


__all__ = ["SCHEMA", "run_fixture"]
