"""Synthetic fail-closed fixture for the Q4 lexical-closure step.

The preceding expansion step leaves a single ``expanded_source_node`` for a
``function_substitution``.  This fixture supplies a fresh, deliberately tiny
graph and checks a caller-provided ``close(graph, expression, env)``
implementation.  It does not decode an archive, walk a source graph, or
evaluate a source value.  The only exact tuple below is an inert environment
binding used to make identity and capture checks visible.

The candidate data shape is node-identity based: its requested substitution
binder is graph node ``2`` (the synthetic bound variable ``x``), and the name
is not used to replace an environment entry.  The candidate replacement is a
source thunk for node ``3`` with a copy of the *outer* environment, before the
new lexical binding is appended to the closure environment.  These are only
synthetic closure-data rules: the fixture does not prove that a source
occurrence has the claimed lexical identity, prove capture avoidance, or
perform the replacement.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping
from copy import deepcopy
from typing import Any


SCHEMA = "original_q4_lexical_closure_fixture_v1"
_EXACT_DYADIC = [0, 3, -6, 2]
_VALID_STATUSES = frozenset({"closure_created", "unresolved", "mismatch"})
_SCOPE_FLAGS = (
    "numeric_field_value_computed",
    "source_value_evaluated",
    "old_errors_removed",
    "occurrence_binding_identity_verified",
    "capture_avoidance_proved",
    "replacement_value_evaluated",
)


def _exact_binding(variable_node: int) -> dict[str, Any]:
    """Return one fresh inert environment binding."""

    return {
        "variable_node": variable_node,
        "binding": {"kind": "exact_dyadic", "tuple": list(_EXACT_DYADIC)},
    }


def _fresh_graph() -> dict[str, Any]:
    """Build the independent synthetic graph used by every case."""

    source_row = {
        "operation": "function_substitution",
        "expression": 1,
        "variable": 2,
        "value": 3,
        "Z_independent_substitution": True,
    }
    # Node 2 is intentionally exact: negative cases replace the row only in
    # their own fresh graph, never in the normal fixture graph.
    nodes = [
        source_row,
        {
            "operation": "synthetic_source_expression",
            "label": "outer-expression-node",
        },
        {"operation": "bound_variable", "name": "x"},
        {
            "operation": "synthetic_replacement_source",
            "label": "replacement-thunk-node",
        },
        {"operation": "bound_variable", "name": "x"},
        {"operation": "bound_variable", "name": "y"},
    ]
    return {"exact_graph_nodes": nodes}


def _expanded_expression(graph: Mapping[str, Any]) -> dict[str, Any]:
    """Return exactly one previous-round expanded source node."""

    row = graph["exact_graph_nodes"][0]
    return {
        "kind": "expanded_source_node",
        "node": 0,
        "operation": "function_substitution",
        "source_row": deepcopy(row),
        "source_children": {
            "expression": {"kind": "source_graph_ref", "node": 1},
            "variable": {"kind": "source_graph_ref", "node": 2},
            "value": {"kind": "source_graph_ref", "node": 3},
        },
        "substitution_applied": False,
    }


def _base_case() -> tuple[dict[str, Any], dict[str, Any], list[dict[str, Any]]]:
    graph = _fresh_graph()
    expression = _expanded_expression(graph)
    # Node 5 is unrelated to the substitution binder and makes preservation of
    # the caller's environment observable without using a numeric value.
    environment = [_exact_binding(5)]
    return graph, expression, environment


def _case(
    name: str,
) -> tuple[dict[str, Any], dict[str, Any], list[dict[str, Any]], str]:
    """Construct one fresh positive or malformed lexical-closure case."""

    graph, expression, environment = _base_case()
    expected = "closure_created"

    if name == "normal":
        pass
    elif name == "same_id_shadow":
        # The old exact binding remains in the outer environment.  The new
        # thunk must be appended and therefore be the nearest node-2 binding.
        environment = [_exact_binding(5), _exact_binding(2)]
    elif name == "same_name_different_id":
        # Node 4 is another binder named x.  Name-based replacement must not
        # treat it as the node-2 substitution variable.
        environment = [_exact_binding(4)]
    elif name == "parameter_environment_unchanged":
        environment = [_exact_binding(5), _exact_binding(4)]
    elif name == "untrusted_source_row_tamper":
        expression["source_row"]["value"] = 4
        expected = "fail_closed"
    elif name == "non_bound_variable":
        graph["exact_graph_nodes"][2] = {
            "operation": "synthetic_not_a_binder",
            "name": "x",
        }
        expected = "fail_closed"
    elif name == "missing_source_field":
        del expression["source_row"]["value"]
        expected = "fail_closed"
    elif name == "missing_graph_source_field":
        del graph["exact_graph_nodes"][0]["value"]
        expected = "fail_closed"
    elif name == "missing_source_child":
        del expression["source_children"]["value"]
        expected = "fail_closed"
    elif name == "out_of_bounds_source_child":
        expression["source_children"]["value"]["node"] = 999
        expected = "fail_closed"
    elif name == "out_of_bounds_graph_row":
        graph["exact_graph_nodes"][0]["value"] = 999
        expected = "fail_closed"
    else:
        raise ValueError("unknown lexical-closure fixture case: " + name)

    return graph, expression, environment, expected


def _scope_is_non_numeric(result: Any) -> bool:
    if not isinstance(result, Mapping):
        return False
    scope = result.get("scope")
    return isinstance(scope, Mapping) and all(scope.get(key) is False for key in _SCOPE_FLAGS)


def _result_status(result: Any) -> str | None:
    if not isinstance(result, Mapping):
        return None
    status = result.get("status")
    return status if isinstance(status, str) else None


def _closure_parts(result: Any) -> tuple[dict[str, Any] | None, list[Any] | None]:
    if not isinstance(result, Mapping):
        return None, None
    expression = result.get("expression")
    if not isinstance(expression, Mapping):
        return None, None
    environment = expression.get("environment")
    if not isinstance(environment, list):
        return None, None
    return dict(expression), environment


def _expected_new_binding(outer_environment: list[dict[str, Any]]) -> dict[str, Any]:
    return {
        "variable_node": 2,
        "binding": {
            "kind": "source_thunk",
            "source_node": 3,
            "environment": deepcopy(outer_environment),
        },
    }


def _valid_closure_shape(
    result: Any,
    outer_environment: list[dict[str, Any]],
) -> tuple[bool, bool, bool]:
    """Return shape and synthetic environment-snapshot checks."""

    expression, environment = _closure_parts(result)
    if expression is None or environment is None:
        return False, False, False
    shape = (
        expression.get("kind") == "source_closure"
        and expression.get("expression_node") == 1
    )
    expected_binding = _expected_new_binding(outer_environment)
    binding = (
        bool(environment)
        and environment[-1] == expected_binding
        and environment[:-1] == outer_environment
        and environment[-1].get("variable_node") == 2
    )
    thunk_environment = None
    if binding:
        thunk_environment = environment[-1]["binding"].get("environment")
    outer_environment_snapshot = (
        binding
        and thunk_environment == outer_environment
        and thunk_environment is not outer_environment
        and all(
            captured is not original
            for captured, original in zip(thunk_environment, outer_environment)
        )
    )
    return shape, binding, outer_environment_snapshot


def _nearest(environment: list[Any], variable_node: int) -> Any:
    """Perform only the synthetic last-binding lookup used by this fixture."""

    for item in reversed(environment):
        if isinstance(item, Mapping) and item.get("variable_node") == variable_node:
            return item.get("binding")
    return None


def _diagnostic(
    name: str,
    expected: str,
    result: Any,
    checks: Mapping[str, bool],
    exception: str | None = None,
) -> dict[str, Any]:
    item: dict[str, Any] = {
        "case": name,
        "expected": expected,
        "status": _result_status(result),
        "passed": all(value is True for value in checks.values()),
        "failed_checks": [key for key, value in checks.items() if value is not True],
    }
    if exception is not None:
        item["exception"] = exception
    return item


def run_fixture(close: Callable[..., Any]) -> dict[str, Any]:
    """Run synthetic lexical-closure and fail-closed admission checks."""

    if not callable(close):
        return {
            "schema": SCHEMA,
            "passed": False,
            "checks": {"close_callable": False},
            "diagnostics": [
                {"case": "input", "expected": "callable", "passed": False}
            ],
            "case_count": 0,
            "test_count": 1,
            "synthetic_only": True,
            "source_code_executed": False,
            "source_graph_walked": False,
            "source_values_evaluated": False,
            "scope": {
                "numeric_field_value_computed": False,
                "source_value_evaluated": False,
                "old_errors_removed": False,
                "occurrence_binding_identity_verified": False,
                "capture_avoidance_proved": False,
                "replacement_value_evaluated": False,
            },
        }

    checks: dict[str, bool] = {"close_callable": True}
    diagnostics: list[dict[str, Any]] = []
    cases = (
        "normal",
        "same_id_shadow",
        "same_name_different_id",
        "parameter_environment_unchanged",
        "untrusted_source_row_tamper",
        "non_bound_variable",
        "missing_source_field",
        "missing_graph_source_field",
        "missing_source_child",
        "out_of_bounds_source_child",
        "out_of_bounds_graph_row",
    )

    for name in cases:
        graph, expression, environment, expected = _case(name)
        original_environment = deepcopy(environment)
        result: Any = None
        exception: str | None = None
        try:
            # The fixture intentionally passes the three public arguments
            # separately; no helper is allowed to infer them from a graph.
            result = close(graph, expression, environment)
        except Exception as exc:  # A raised error is not a returned contract.
            exception = type(exc).__name__

        prefix = name + "_"
        checks[prefix + "input_environment_unchanged"] = (
            environment == original_environment
        )
        checks[prefix + "scope_non_numeric"] = _scope_is_non_numeric(result)
        status = _result_status(result)
        checks[prefix + "status_enum"] = status in _VALID_STATUSES

        if expected == "closure_created":
            shape, binding, outer_snapshot = _valid_closure_shape(
                result, original_environment
            )
            checks[prefix + "status"] = status == "closure_created"
            checks[prefix + "shape"] = shape
            checks[prefix + "new_binding"] = binding
            checks[prefix + "synthetic_outer_environment_snapshot"] = outer_snapshot

            expression_out, environment_out = _closure_parts(result)
            checks[prefix + "synthetic_nearest_variable_id"] = (
                expression_out is not None
                and environment_out is not None
                and _nearest(environment_out, 2)
                == _expected_new_binding(original_environment)["binding"]
            )
            if name == "same_name_different_id":
                checks[prefix + "synthetic_same_name_id_preserved"] = (
                    expression_out is not None
                    and environment_out is not None
                    and _nearest(environment_out, 4)
                    == {"kind": "exact_dyadic", "tuple": list(_EXACT_DYADIC)}
                )
        else:
            checks[prefix + "fail_closed"] = status in {"mismatch", "unresolved"}
            checks[prefix + "no_closure_created"] = status != "closure_created"

        diagnostics.append(_diagnostic(name, expected, result, {
            key: value for key, value in checks.items() if key.startswith(prefix)
        }, exception))

    return {
        "schema": SCHEMA,
        "passed": bool(checks) and all(value is True for value in checks.values()),
        "checks": checks,
        "diagnostics": diagnostics,
        "case_count": len(cases),
        "test_count": len(checks),
        "synthetic_only": True,
        "source_code_executed": False,
        "source_graph_walked": False,
        "source_values_evaluated": False,
        "scope": {
            "numeric_field_value_computed": False,
            "source_value_evaluated": False,
            "old_errors_removed": False,
            "occurrence_binding_identity_verified": False,
            "capture_avoidance_proved": False,
            "replacement_value_evaluated": False,
        },
        "occurrence_binding_identity_verified": False,
        "capture_avoidance_proved": False,
        "replacement_value_evaluated": False,
    }


__all__ = ["SCHEMA", "run_fixture"]
