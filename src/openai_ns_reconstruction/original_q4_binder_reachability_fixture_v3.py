"""Fresh v3 synthetic checks for Q4 binder reachability.

Every case builds an independent ``exact_graph_nodes`` table with five
unreachable v3 dummy rows before its live handles.  Binder names and invalid
references are v3-specific, so this fixture does not reuse the v1/v2 input
graphs.  It never imports or runs an analyzer, evaluates a value, or walks a
real graph.

The occurrence report is checked by target identity, not by the total number
of bound-variable leaves: a shared-context body contains an unrelated beta
leaf in addition to its shadowed and outer alpha leaves, and the replacement
body contains an unrelated beta leaf in addition to its outer alpha leaf.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping
from typing import Any


SCHEMA = "original_q4_binder_reachability_fixture_v3"
_VALID_STATUSES = frozenset({"complete", "unresolved", "mismatch"})


def _add(nodes: list[dict[str, Any]], operation: str, **fields: Any) -> int:
    row: dict[str, Any] = {"operation": operation}
    row.update(fields)
    nodes.append(row)
    return len(nodes) - 1


def _substitution(
    nodes: list[dict[str, Any]], expression: Any, variable: Any, value: Any
) -> int:
    return _add(
        nodes,
        "function_substitution",
        expression=expression,
        variable=variable,
        value=value,
    )


def _base() -> tuple[list[dict[str, Any]], dict[str, int]]:
    """Build five unreachable v3 rows before all live node handles."""

    nodes: list[dict[str, Any]] = [
        {
            "operation": "fresh_v3_unreachable_dummy",
            "label": "dummy-v3-delta",
        },
        {
            "operation": "fresh_v3_unreachable_dummy",
            "label": "dummy-v3-epsilon",
        },
        {
            "operation": "fresh_v3_unreachable_dummy",
            "label": "dummy-v3-zeta",
        },
        {
            "operation": "fresh_v3_unreachable_dummy",
            "label": "dummy-v3-eta",
        },
        {
            "operation": "fresh_v3_unreachable_dummy",
            "label": "dummy-v3-theta",
        },
    ]
    binder_alpha = _add(nodes, "bound_variable", name="fresh_v3_alpha")
    binder_alpha_same_name = _add(
        nodes, "bound_variable", name="fresh_v3_alpha"
    )
    binder_beta = _add(nodes, "bound_variable", name="fresh_v3_beta")
    return nodes, {
        "alpha": binder_alpha,
        "same_name_alpha": binder_alpha_same_name,
        "beta": binder_beta,
    }


def _case(name: str) -> tuple[dict[str, Any], int, str | None, frozenset[str]]:
    """Build a fresh v3 graph and expected status contract."""

    nodes, ids = _base()
    alpha = ids["alpha"]
    beta = ids["beta"]

    if name == "direct_same_binder_id":
        substitution = _substitution(nodes, alpha, alpha, beta)
        expected = "complete"
        allowed = frozenset({expected})
    elif name == "nested_same_binder_body_shadow_value_outer":
        nested = _substitution(nodes, alpha, alpha, alpha)
        substitution = _substitution(nodes, nested, alpha, beta)
        expected = "complete"
        allowed = frozenset({expected})
    elif name == "integral_string_binder_bounds_outer":
        integral = _add(
            nodes,
            "definite_integral",
            variable="fresh_v3_alpha",
            integrand=alpha,
            lower=alpha,
            upper=alpha,
        )
        substitution = _substitution(nodes, integral, alpha, beta)
        expected = "unresolved"
        allowed = frozenset({expected})
    elif name == "same_name_different_id":
        substitution = _substitution(nodes, ids["same_name_alpha"], alpha, beta)
        expected = None
        allowed = frozenset({"unresolved", "mismatch"})
    elif name == "opaque_unknown_active_body_not_read":
        unknown = _add(
            nodes,
            "opaque_unknown",
            active_body=30_000,
            reason="fresh v3 opaque source",
        )
        substitution = _substitution(nodes, unknown, alpha, beta)
        expected = "unresolved"
        allowed = frozenset({expected})
    elif name == "lazy_branch_active_body_not_read":
        lazy = _add(
            nodes,
            "lazy_branch",
            condition=beta,
            active_body=30_001,
            inactive_body=beta,
        )
        substitution = _substitution(nodes, lazy, alpha, beta)
        expected = "unresolved"
        allowed = frozenset({expected})
    elif name == "cycle_fail_closed":
        cycle = _add(nodes, "sum", arguments=[])
        nodes[cycle]["arguments"] = [cycle]
        substitution = _substitution(nodes, cycle, alpha, beta)
        expected = None
        allowed = frozenset({"unresolved", "mismatch"})
    elif name == "invalid_expression_id_fail_closed":
        substitution = _substitution(nodes, 30_002, alpha, beta)
        expected = None
        allowed = frozenset({"unresolved", "mismatch"})
    elif name == "invalid_substitution_id_fail_closed":
        substitution = 30_003
        expected = None
        allowed = frozenset({"unresolved", "mismatch"})
    elif name == "malformed_binder_fail_closed":
        malformed = _add(nodes, "bound_variable")
        substitution = _substitution(nodes, malformed, malformed, beta)
        expected = None
        allowed = frozenset({"unresolved", "mismatch"})
    elif name == "shared_node_distinct_shadow_contexts":
        # Body leaves are shadowed alpha, unrelated beta, and outer alpha.
        shared = alpha
        shadowed_branch = _substitution(nodes, shared, alpha, beta)
        branches = _add(
            nodes,
            "sum",
            arguments=[shadowed_branch, shared],
        )
        substitution = _substitution(nodes, branches, alpha, beta)
        expected = "complete"
        allowed = frozenset({expected})
    elif name == "replacement_occurrence_outer_environment":
        # Body leaves are unrelated beta and the nested replacement alpha in
        # the old outer environment.
        nested_body = _add(nodes, "sum", arguments=[beta])
        nested = _substitution(nodes, nested_body, alpha, alpha)
        substitution = _substitution(nodes, nested, alpha, beta)
        expected = "complete"
        allowed = frozenset({expected})
    else:
        raise ValueError("unknown v3 binder reachability case: " + name)

    return {"exact_graph_nodes": nodes}, substitution, expected, allowed


def _scope(result: Any) -> Mapping[str, Any] | None:
    if not isinstance(result, Mapping):
        return None
    value = result.get("scope")
    return value if isinstance(value, Mapping) else None


def _status(result: Any) -> str | None:
    if not isinstance(result, Mapping):
        return None
    value = result.get("status")
    return value if isinstance(value, str) else None


def _scope_is_non_numeric(result: Any) -> bool:
    scope = _scope(result)
    return (
        scope is not None
        and scope.get("numeric_field_value_computed") is False
        and scope.get("old_errors_removed") is False
        and scope.get("complete_identity_proved") is False
    )


def _counts_match(result: Any, expected: Mapping[str, int]) -> bool:
    if not isinstance(result, Mapping):
        return False
    counts = result.get("counts")
    return isinstance(counts, Mapping) and all(
        counts.get(key) == value for key, value in expected.items()
    )


def _occurrence_flags(result: Any, side: str) -> list[tuple[Any, Any, Any]]:
    """Return flags only for target-identity records; ignore beta leaves."""

    if not isinstance(result, Mapping):
        return []
    occurrences = result.get("occurrences")
    if not isinstance(occurrences, Mapping):
        return []
    rows = occurrences.get(side)
    if not isinstance(rows, list):
        return []
    return [
        (
            row.get("is_target_id"),
            row.get("shadowed"),
            row.get("target_id_unshadowed"),
        )
        for row in rows
        if isinstance(row, Mapping) and row.get("is_target_id") is True
    ]


def _shadow_record_for_expression(result: Any) -> bool:
    if not isinstance(result, Mapping):
        return False
    records = result.get("shadow_records")
    if not isinstance(records, Mapping):
        return False
    body = records.get("body")
    return isinstance(body, list) and any(
        isinstance(item, Mapping)
        and item.get("scope") == "expression"
        and item.get("shadows_target") is True
        for item in body
    )


def _diagnostic(
    name: str,
    result: Any,
    checks: Mapping[str, bool],
    exception: str | None,
) -> dict[str, Any]:
    item: dict[str, Any] = {
        "case": name,
        "status": _status(result),
        "passed": all(value is True for value in checks.values()),
        "failed_checks": [key for key, value in checks.items() if value is not True],
    }
    if exception is not None:
        item["exception"] = exception
    return item


def run_fixture(analyze: Callable[..., Any]) -> dict[str, Any]:
    """Run fresh v3 binder reachability and fail-closed checks."""

    if not callable(analyze):
        return {
            "schema": SCHEMA,
            "passed": False,
            "checks": {"analyze_callable": False},
            "diagnostics": [
                {"case": "input", "passed": False, "expected": "callable"}
            ],
            "case_count": 0,
            "test_count": 1,
            "synthetic_only": True,
            "source_code_executed": False,
            "source_graph_walked": False,
            "scope": {
                "numeric_field_value_computed": False,
                "old_errors_removed": False,
                "complete_identity_proved": False,
            },
        }

    checks: dict[str, bool] = {"analyze_callable": True}
    diagnostics: list[dict[str, Any]] = []
    cases = (
        "direct_same_binder_id",
        "nested_same_binder_body_shadow_value_outer",
        "integral_string_binder_bounds_outer",
        "same_name_different_id",
        "opaque_unknown_active_body_not_read",
        "lazy_branch_active_body_not_read",
        "cycle_fail_closed",
        "invalid_expression_id_fail_closed",
        "invalid_substitution_id_fail_closed",
        "malformed_binder_fail_closed",
        "shared_node_distinct_shadow_contexts",
        "replacement_occurrence_outer_environment",
    )

    for name in cases:
        graph, substitution_node, expected, allowed = _case(name)
        result: Any = None
        exception: str | None = None
        try:
            result = analyze(graph, substitution_node)
        except Exception as exc:  # A raised error is not a returned contract.
            exception = type(exc).__name__

        prefix = name + "_"
        status = _status(result)
        scope = _scope(result)
        checks[prefix + "status_enum"] = status in _VALID_STATUSES
        checks[prefix + "scope_shape"] = _scope_is_non_numeric(result)
        checks[prefix + "numeric_field_not_computed"] = (
            scope is not None
            and scope.get("numeric_field_value_computed") is False
        )
        checks[prefix + "old_errors_not_removed"] = (
            scope is not None and scope.get("old_errors_removed") is False
        )
        checks[prefix + "status_allowed"] = status in allowed

        if expected is not None:
            checks[prefix + "status"] = status == expected
            checks[prefix + "identity_not_proved"] = (
                scope is not None
                and scope.get("complete_identity_proved") is False
            )
        else:
            checks[prefix + "not_complete"] = status != "complete"
            checks[prefix + "identity_not_proved"] = (
                scope is not None
                and scope.get("complete_identity_proved") is False
            )

        if name == "direct_same_binder_id":
            checks[prefix + "body_target_occurrence_counts"] = _counts_match(
                result,
                {
                    "body_occurrences": 1,
                    "value_occurrences": 1,
                    "body_target_id_occurrences": 1,
                    "body_target_id_unshadowed_occurrences": 1,
                    "body_target_id_shadowed_occurrences": 0,
                    "value_target_id_occurrences": 0,
                },
            )
            checks[prefix + "body_occurrence_is_unshadowed"] = (
                _occurrence_flags(result, "body") == [(True, False, True)]
            )
        elif name == "nested_same_binder_body_shadow_value_outer":
            checks[prefix + "body_shadow_and_outer_counts"] = _counts_match(
                result,
                {
                    "body_occurrences": 2,
                    "value_occurrences": 1,
                    "body_target_id_occurrences": 2,
                    "body_target_id_unshadowed_occurrences": 1,
                    "body_target_id_shadowed_occurrences": 1,
                    "value_target_id_occurrences": 0,
                    "body_shadow_records": 1,
                },
            )
            checks[prefix + "body_has_shadowed_and_outer_records"] = (
                sorted(
                    _occurrence_flags(result, "body"),
                    key=lambda item: (bool(item[1]), bool(item[2])),
                )
                == [(True, False, True), (True, True, False)]
                and _shadow_record_for_expression(result)
            )
        elif name == "integral_string_binder_bounds_outer":
            checks[prefix + "bounds_scanned_integrand_not_scanned"] = _counts_match(
                result,
                {
                    "body_occurrences": 2,
                    "value_occurrences": 1,
                    "body_target_id_occurrences": 2,
                    "body_target_id_unshadowed_occurrences": 2,
                    "body_target_id_shadowed_occurrences": 0,
                    "body_shadow_records": 1,
                    "value_target_id_occurrences": 0,
                },
            )
        elif name == "shared_node_distinct_shadow_contexts":
            checks[prefix + "shared_context_occurrence_counts"] = _counts_match(
                result,
                {
                    "body_occurrences": 3,
                    "value_occurrences": 1,
                    "body_target_id_occurrences": 2,
                    "body_target_id_unshadowed_occurrences": 1,
                    "body_target_id_shadowed_occurrences": 1,
                    "body_shadow_records": 1,
                    "value_target_id_occurrences": 0,
                },
            )
            checks[prefix + "shared_node_has_both_context_records"] = (
                sorted(
                    _occurrence_flags(result, "body"),
                    key=lambda item: (bool(item[1]), bool(item[2])),
                )
                == [(True, False, True), (True, True, False)]
            )
        elif name == "replacement_occurrence_outer_environment":
            checks[prefix + "replacement_outer_occurrence_counts"] = _counts_match(
                result,
                {
                    "body_occurrences": 2,
                    "value_occurrences": 1,
                    "body_target_id_occurrences": 1,
                    "body_target_id_unshadowed_occurrences": 1,
                    "body_target_id_shadowed_occurrences": 0,
                    "body_shadow_records": 1,
                    "value_target_id_occurrences": 0,
                },
            )
            checks[prefix + "replacement_occurrence_is_outer"] = (
                _occurrence_flags(result, "body") == [(True, False, True)]
            )

        diagnostics.append(
            _diagnostic(
                name,
                result,
                {
                    key: value
                    for key, value in checks.items()
                    if key.startswith(prefix)
                },
                exception,
            )
        )

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
        "numeric_values_computed": False,
        "scope": {
            "numeric_field_value_computed": False,
            "old_errors_removed": False,
            "complete_identity_proved": False,
        },
    }


__all__ = ["SCHEMA", "run_fixture"]
