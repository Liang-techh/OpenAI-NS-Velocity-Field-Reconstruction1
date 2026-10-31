"""Fresh synthetic contract fixture for one-layer Q4 endpoint expansion.

The fixture intentionally has no connection to a saved graph packet.  A caller
passes an implementation of ``expand(graph, prior_result)`` to
``run_fixture``.  The implementation is expected to replace only an existing
``opaque_source_ref`` whose graph row is one of ``negative``,
``analytic_unary`` or ``function_substitution``.

The exact source shape used here is deliberately small and typed:

* the result is ``kind == "expanded_source_node"``;
* ``source_ref`` is the complete old opaque reference and ``source_row`` is
  the complete graph row, so source ids, operation names, assumptions and
  extra provenance fields remain available;
* direct children are in ``source_children`` as
  ``{kind: "source_graph_ref", node: ...}`` leaves;
* ``substitution_applied`` is always ``False``.

Graph references are intentionally leaves.  The fixture therefore catches an
implementation that recursively expands a newly referenced graph node, does
numeric evaluation, or performs substitution.  The ``flat_branches`` field is
also checked byte-for-byte at the Python object level; its inactive body is
not part of this v1 expansion.
"""

from __future__ import annotations

from copy import deepcopy
from typing import Any, Callable, Mapping


Expand = Callable[[dict[str, Any], dict[str, Any]], dict[str, Any]]
_EXPANDABLE = frozenset({"negative", "analytic_unary", "function_substitution"})


def _opaque(node: int, operation: str, reason: str, **extra: Any) -> dict[str, Any]:
    """Make a v1 opaque AST ref, with deliberately visible provenance."""

    value: dict[str, Any] = {
        "kind": "opaque_source_ref",
        "node": node,
        "operation": operation,
        "reason": reason,
    }
    value.update(extra)
    return value


def _graph_ref(node: int) -> dict[str, Any]:
    return {"kind": "source_graph_ref", "node": node}


def _exact_source_expression(row: Mapping[str, Any], opaque: Mapping[str, Any]) -> dict[str, Any]:
    """Build the fixture's one-layer expected source expression.

    This only records direct integer references.  It deliberately does not
    inspect any referenced row.
    """

    operation = row["operation"]
    result: dict[str, Any] = {
        "kind": "expanded_source_node",
        "operation": operation,
        "node": opaque["node"],
        "source_ref": deepcopy(dict(opaque)),
        "source_row": deepcopy(dict(row)),
        "source_children": {},
        "substitution_applied": False,
    }
    if "source_reference_id" in row:
        result["source_reference_id"] = deepcopy(row["source_reference_id"])
    if "assumptions" in row:
        result["assumptions"] = deepcopy(row["assumptions"])

    if operation == "negative":
        result["source_children"] = {"argument": _graph_ref(row["argument"])}
    elif operation == "analytic_unary":
        result["source_children"] = {"argument": _graph_ref(row["argument"])}
    elif operation == "function_substitution":
        result["Z_independent_substitution"] = row["Z_independent_substitution"]
        result["source_children"] = {
            "expression": _graph_ref(row["expression"]),
            "variable": _graph_ref(row["variable"]),
            "value": _graph_ref(row["value"]),
        }
    return result


def _scope() -> dict[str, Any]:
    return {
        "status": "partial",
        "owner": "q4-endpoint-v1",
        "parameter_owner": "outer-scoped-owner",
        "coordinate_owner": "outer-scoped-owner",
        "error": {"kind": "unresolved", "code": "opaque-source"},
    }


def _outer_parameter() -> dict[str, Any]:
    return {
        "kind": "scoped_parameter",
        "name": "lambda_parameter",
        "owner": "outer-scoped-owner",
        "scope": "outer",
        "status": "symbolic",
    }


def _outer_coordinate() -> dict[str, Any]:
    return {
        "kind": "scoped_coordinate",
        "name": "sc",
        "owner": "outer-scoped-owner",
        "scope": "outer",
        "status": "symbolic",
    }


def _flat_branches() -> list[dict[str, Any]]:
    return [
        {
            "label": "inactive-flat-body",
            "active": False,
            "body": _opaque(
                80,
                "inactive_flat_fixture_body",
                "v1 flat branch is inactive and must not be expanded",
                branch_owner="flat-v1",
            ),
        }
    ]


def _prior(expression: dict[str, Any]) -> dict[str, Any]:
    return {
        "status": "unresolved",
        "expression": expression,
        "flat_branches": _flat_branches(),
    }


def _negative_case() -> dict[str, Any]:
    graph = {
        "exact_graph_nodes": [
            {
                "operation": "negative",
                "argument": 1,
                "source_reference_id": "neg-source-101",
                "assumptions": {"sign": "formal", "scope": "endpoint"},
            },
            {
                "operation": "analytic_unary",
                "name": "exp",
                "argument": 2,
                "source_reference_id": "child-unary-102",
                "assumptions": {"child": "must-remain-a-ref"},
            },
            {
                "operation": "scoped_coordinate_leaf",
                "name": "sc",
                "source_reference_id": "coordinate-103",
            },
        ]
    }
    opaque = _opaque(
        0,
        "negative",
        "negative source was opaque in endpoint v1",
        source_reference_id="opaque-neg-ref-101",
        assumptions={"sign": "formal", "scope": "endpoint"},
    )
    expression = {
        "kind": "sum",
        "terms": [
            _outer_parameter(),
            opaque,
        ],
        "scoped_parameter": _outer_parameter(),
        "scoped_coordinate": _outer_coordinate(),
        "scope": _scope(),
    }
    expected = deepcopy(expression)
    expected["terms"][1] = _exact_source_expression(graph["exact_graph_nodes"][0], opaque)
    return {
        "name": "negative_one_layer",
        "graph": graph,
        "prior_result": _prior(expression),
        "expected_status": "expanded",
        "expected_expression": expected,
        "one_layer_path": ("terms", 1, "source_children", "argument"),
        "one_layer_node": 1,
    }


def _unary_case() -> dict[str, Any]:
    graph = {
        "exact_graph_nodes": [
            {
                "operation": "analytic_unary",
                "name": "exp",
                "argument": 1,
                "source_reference_id": "unary-source-201",
                "assumptions": {"analytic_name": "exp", "domain": "formal"},
            },
            {
                "operation": "negative",
                "argument": 2,
                "source_reference_id": "nested-negative-202",
                "assumptions": {"nested": "opaque-ref-only"},
            },
            {
                "operation": "scoped_parameter_leaf",
                "name": "lambda_parameter",
                "source_reference_id": "parameter-203",
            },
        ]
    }
    opaque = _opaque(
        0,
        "analytic_unary",
        "analytic unary source was opaque in endpoint v1",
        source_reference_id="opaque-unary-ref-201",
        assumptions={"analytic_name": "exp", "domain": "formal"},
    )
    expression = {
        "kind": "product",
        "factors": [opaque, _outer_coordinate()],
        "scoped_parameter": _outer_parameter(),
        "scoped_coordinate": _outer_coordinate(),
        "scope": _scope(),
    }
    expected = deepcopy(expression)
    expected["factors"][0] = _exact_source_expression(graph["exact_graph_nodes"][0], opaque)
    return {
        "name": "analytic_unary_one_layer",
        "graph": graph,
        "prior_result": _prior(expression),
        "expected_status": "expanded",
        "expected_expression": expected,
        "one_layer_path": ("factors", 0, "source_children", "argument"),
        "one_layer_node": 1,
    }


def _substitution_case() -> dict[str, Any]:
    graph = {
        "exact_graph_nodes": [
            {
                "operation": "function_substitution",
                "expression": 1,
                "variable": 2,
                "value": 3,
                "Z_independent_substitution": True,
                "source_reference_id": "substitution-source-301",
                "assumptions": {
                    "binder": "lambda_x",
                    "coordinate": "sc",
                    "Z_independent": True,
                },
            },
            {
                "operation": "analytic_unary",
                "name": "exp",
                "argument": 4,
                "source_reference_id": "substitution-expression-302",
            },
            {
                "operation": "bound_variable",
                "name": "lambda_x",
                "binder": True,
                "source_reference_id": "substitution-binder-303",
            },
            {
                "operation": "scoped_parameter_leaf",
                "name": "lambda_parameter",
                "source_reference_id": "substitution-value-304",
            },
            {
                "operation": "scoped_coordinate_leaf",
                "name": "sc",
                "source_reference_id": "substitution-coordinate-305",
            },
        ]
    }
    opaque = _opaque(
        0,
        "function_substitution",
        "function substitution source was opaque in endpoint v1",
        source_reference_id="opaque-substitution-ref-301",
        assumptions={
            "binder": "lambda_x",
            "coordinate": "sc",
            "Z_independent": True,
        },
    )
    binder = {
        "kind": "lambda_binder",
        "node": 2,
        "name": "lambda_x",
        "owner": "lambda-local-owner",
        "scope": "lambda",
    }
    expression = {
        "kind": "lambda",
        "binder": binder,
        "body": opaque,
        "scoped_parameter": _outer_parameter(),
        "scoped_coordinate": _outer_coordinate(),
        "scope": _scope(),
    }
    expected = deepcopy(expression)
    expected["body"] = _exact_source_expression(graph["exact_graph_nodes"][0], opaque)
    return {
        "name": "function_substitution_preserves_lambda_binder",
        "graph": graph,
        "prior_result": _prior(expression),
        "expected_status": "expanded",
        "expected_expression": expected,
        "binder_path": ("binder",),
        "one_layer_path": ("body", "source_children", "expression"),
        "one_layer_node": 1,
    }


def _unknown_case() -> dict[str, Any]:
    graph = {
        "exact_graph_nodes": [
            {
                "operation": "original_unknown_fixture_operation",
                "source_reference_id": "unknown-source-401",
                "assumptions": {"must": "remain-opaque"},
            }
        ]
    }
    opaque = _opaque(
        0,
        "original_unknown_fixture_operation",
        "unknown operation is intentionally unresolved",
        source_reference_id="opaque-unknown-ref-401",
    )
    expression = {
        "kind": "sum",
        "terms": [opaque, _outer_parameter()],
        "scoped_parameter": _outer_parameter(),
        "scoped_coordinate": _outer_coordinate(),
        "scope": _scope(),
    }
    return {
        "name": "unknown_operation_stays_opaque",
        "graph": graph,
        "prior_result": _prior(expression),
        "expected_status": "unresolved",
        "expected_expression": expression,
    }


def _nonexpandable_case() -> dict[str, Any]:
    graph = {
        "exact_graph_nodes": [
            {
                "operation": "sum",
                "arguments": [1, 2],
                "source_reference_id": "sum-source-501",
            },
            {"operation": "scoped_parameter_leaf", "name": "lambda_parameter"},
            {"operation": "scoped_coordinate_leaf", "name": "sc"},
        ]
    }
    opaque = _opaque(
        0,
        "sum",
        "sum is a v1 opaque operation and is not an expansion target",
    )
    expression = {
        "kind": "product",
        "factors": [opaque, _outer_coordinate()],
        "scoped_parameter": _outer_parameter(),
        "scoped_coordinate": _outer_coordinate(),
        "scope": _scope(),
    }
    return {
        "name": "nonexpandable_sum_stays_unresolved",
        "graph": graph,
        "prior_result": _prior(expression),
        "expected_status": "unresolved",
        "expected_expression": expression,
    }


def _mismatched_operation_case() -> dict[str, Any]:
    graph = {
        "exact_graph_nodes": [
            {
                "operation": "analytic_unary",
                "name": "exp",
                "argument": 1,
                "source_reference_id": "tampered-source-601",
            },
            {"operation": "scoped_parameter_leaf", "name": "lambda_parameter"},
        ]
    }
    expression = {
        "kind": "opaque",
        "value": _opaque(0, "negative", "operation was tampered in prior AST"),
        "scoped_parameter": _outer_parameter(),
        "scoped_coordinate": _outer_coordinate(),
        "scope": _scope(),
    }
    return {
        "name": "tampered_prior_operation_is_mismatch",
        "graph": graph,
        "prior_result": _prior(expression),
        "expected_status": "mismatch",
        "expected_expression": expression,
    }


def _missing_field_case() -> dict[str, Any]:
    graph = {
        "exact_graph_nodes": [
            {
                "operation": "negative",
                "source_reference_id": "missing-argument-source-701",
            }
        ]
    }
    expression = {
        "kind": "opaque",
        "value": _opaque(0, "negative", "required argument was removed"),
        "scoped_parameter": _outer_parameter(),
        "scoped_coordinate": _outer_coordinate(),
        "scope": _scope(),
    }
    return {
        "name": "missing_negative_argument_is_mismatch",
        "graph": graph,
        "prior_result": _prior(expression),
        "expected_status": "mismatch",
        "expected_expression": expression,
    }


def _wrong_graph_ref_case() -> dict[str, Any]:
    graph = {
        "exact_graph_nodes": [
            {
                "operation": "negative",
                "argument": 1,
                "source_reference_id": "wrong-ref-source-801",
            },
            {"operation": "scoped_coordinate_leaf", "name": "sc"},
        ]
    }
    expression = {
        "kind": "opaque",
        "value": _opaque(99, "negative", "graph reference is outside exact node list"),
        "scoped_parameter": _outer_parameter(),
        "scoped_coordinate": _outer_coordinate(),
        "scope": _scope(),
    }
    return {
        "name": "wrong_graph_reference_is_mismatch",
        "graph": graph,
        "prior_result": _prior(expression),
        "expected_status": "mismatch",
        "expected_expression": expression,
    }


def _cases() -> tuple[dict[str, Any], ...]:
    """Build every case independently; no graph list is shared between cases."""

    return (
        _negative_case(),
        _unary_case(),
        _substitution_case(),
        _unknown_case(),
        _nonexpandable_case(),
        _mismatched_operation_case(),
        _missing_field_case(),
        _wrong_graph_ref_case(),
    )


def _at_path(value: Any, path: tuple[Any, ...]) -> Any:
    current = value
    for part in path:
        if isinstance(current, Mapping):
            current = current.get(part)
        elif isinstance(current, (list, tuple)) and isinstance(part, int) and not isinstance(part, bool):
            if part < 0 or part >= len(current):
                return None
            current = current[part]
        else:
            return None
    return current


def _scope_receipt(before: Any, after: Any) -> dict[str, Any]:
    before_scope = before.get("scope") if isinstance(before, Mapping) else None
    after_scope = after.get("scope") if isinstance(after, Mapping) else None
    before_is_map = isinstance(before_scope, Mapping)
    after_is_map = isinstance(after_scope, Mapping)
    if not before_is_map or not after_is_map:
        return {
            "passed": False,
            "status_preserved": False,
            "no_numerical_field": False,
            "owners_preserved": False,
            "error_preserved": False,
        }
    owner_keys = ("owner", "parameter_owner", "coordinate_owner")
    owners_preserved = all(after_scope.get(key) == before_scope.get(key) for key in owner_keys)
    no_numerical = all("numerical" not in str(key).lower() for key in after_scope)
    status_preserved = after_scope.get("status") == before_scope.get("status") and after_scope.get("status") != "numerical"
    error_preserved = after_scope.get("error") == before_scope.get("error")
    return {
        "passed": bool(status_preserved and no_numerical and owners_preserved and error_preserved),
        "status_preserved": status_preserved,
        "no_numerical_field": no_numerical,
        "owners_preserved": owners_preserved,
        "error_preserved": error_preserved,
    }


def _case_receipt(case: Mapping[str, Any], result: Any, original_graph: Any, original_prior: Any, graph: Any, prior: Any) -> dict[str, Any]:
    expected_status = case["expected_status"]
    is_map = isinstance(result, Mapping)
    actual_status = result.get("status") if is_map else None
    actual_expression = result.get("expression") if is_map else None
    expected_expression = case["expected_expression"]
    expression_ok = is_map and actual_expression == expected_expression
    status_ok = actual_status == expected_status
    flat_ok = is_map and result.get("flat_branches") == original_prior.get("flat_branches")
    graph_untouched = graph == original_graph
    prior_untouched = prior == original_prior
    scope = _scope_receipt(original_prior.get("expression"), actual_expression)

    one_layer = {"checked": False, "passed": True}
    if "one_layer_path" in case:
        one_layer["checked"] = True
        leaf = _at_path(actual_expression, case["one_layer_path"])
        one_layer["passed"] = (
            isinstance(leaf, Mapping)
            and leaf.get("kind") == "source_graph_ref"
            and leaf.get("node") == case["one_layer_node"]
        )

    binder = {"checked": False, "passed": True}
    if "binder_path" in case:
        binder["checked"] = True
        binder["passed"] = _at_path(actual_expression, case["binder_path"]) == _at_path(original_prior["expression"], case["binder_path"])

    outer_layer = {
        "parameter_preserved": _at_path(actual_expression, ("scoped_parameter",)) == _at_path(original_prior["expression"], ("scoped_parameter",)),
        "coordinate_preserved": _at_path(actual_expression, ("scoped_coordinate",)) == _at_path(original_prior["expression"], ("scoped_coordinate",)),
    }
    outer_layer["passed"] = bool(outer_layer["parameter_preserved"] and outer_layer["coordinate_preserved"])

    passed = bool(
        status_ok
        and expression_ok
        and flat_ok
        and graph_untouched
        and prior_untouched
        and scope["passed"]
        and one_layer["passed"]
        and binder["passed"]
        and outer_layer["passed"]
    )
    return {
        "name": case["name"],
        "expected_status": expected_status,
        "actual_status": actual_status,
        "passed": passed,
        "expression_matches_expected": expression_ok,
        "flat_branches_unchanged": flat_ok,
        "graph_input_unchanged": graph_untouched,
        "prior_input_unchanged": prior_untouched,
        "scope": scope,
        "outer_scoped_inputs": outer_layer,
        "one_layer_graph_refs": one_layer,
        "lambda_binder": binder,
    }


def run_fixture(expand: Expand) -> dict[str, Any]:
    """Run all fresh synthetic contract cases against ``expand``.

    The return value is a receipt rather than an assertion exception so a
    parent proof can report which contract edge failed.  No case uses a real
    saved graph, and no result is interpreted numerically.
    """

    if not callable(expand):
        raise TypeError("run_fixture expects a callable expand(graph, prior_result)")

    receipts: list[dict[str, Any]] = []
    for case in _cases():
        graph = case["graph"]
        prior = case["prior_result"]
        original_graph = deepcopy(graph)
        original_prior = deepcopy(prior)
        try:
            result = expand(graph, prior)
        except Exception as exc:  # A fixture receipt should identify the bad case.
            receipts.append(
                {
                    "name": case["name"],
                    "expected_status": case["expected_status"],
                    "actual_status": "raised",
                    "passed": False,
                    "error_type": type(exc).__name__,
                    "error": str(exc),
                }
            )
            continue
        receipts.append(_case_receipt(case, result, original_graph, original_prior, graph, prior))

    passed = bool(receipts) and all(receipt["passed"] for receipt in receipts)
    return {
        "fixture": "original_q4_endpoint_expand_fixture_v1",
        "api": "expand(graph, prior_result)",
        "passed": passed,
        "all_passed": passed,
        "case_count": len(receipts),
        "cases": receipts,
        "contract": {
            "expandable_operations": sorted(_EXPANDABLE),
            "one_layer_only": True,
            "numeric_evaluation": False,
            "substitution_application": False,
            "flat_branches_recomputed": False,
            "unknown_operation_status": "unresolved",
            "graph_reference_mismatch_status": "mismatch",
            "typed_result_kind": "expanded_source_node",
        },
    }


__all__ = ["run_fixture"]
