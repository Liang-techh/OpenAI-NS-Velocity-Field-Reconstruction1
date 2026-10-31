"""Fresh synthetic receipts for the Q4 source-only consumers.

The fixtures accept caller-provided consumers and never import, run, or
evaluate a source implementation.  Every graph is rebuilt per case.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping
from typing import Any


SCHEMA = "original_q4_integrand_binding_fixture_v1"
_VALID_STATUSES = frozenset({"observed", "mismatch"})
_INTEGRAND_SCOPE = (
    "integral_binder_identity_proved",
    "numeric_field_value_computed",
    "old_errors_removed",
)


def _add(nodes: list[dict[str, Any]], operation: str, **fields: Any) -> int:
    row = {"operation": operation}
    row.update(fields)
    nodes.append(row)
    return len(nodes) - 1


def _base(tag: str) -> tuple[list[dict[str, Any]], int, int]:
    return [
        {"operation": f"{tag}_unreachable_dummy", "active_body": 90_000},
        {"operation": "exact_rational", "numerator": 0, "denominator": 1},
        {"operation": "exact_rational", "numerator": 1, "denominator": 1},
    ], 1, 2


def _binder(nodes: list[dict[str, Any]], label: str) -> int:
    return _add(nodes, "bound_variable", name="fresh_moment_t", record=label)


def _frontier(integral: int, body: int, **extra: Any) -> dict[str, Any]:
    value = {
        "node_id": integral,
        "integrand_node_id": body,
        "binder_name": "fresh_moment_t",
        "operation": "definite_integral",
    }
    value.update(extra)
    return value


def _integral(
    nodes: list[dict[str, Any]],
    body: int,
    lower: int,
    upper: int,
    *,
    operation: str = "definite_integral",
    variable: str = "fresh_moment_t",
    **extra: Any,
) -> int:
    return _add(
        nodes,
        operation,
        variable=variable,
        integrand=body,
        lower=lower,
        upper=upper,
        **extra,
    )


def _observed(
    nodes: list[dict[str, Any]],
    frontier: dict[str, Any],
    *,
    unshadowed: list[int],
    shadowed: list[int],
    conflict: bool,
    visited: int,
    complete: bool,
    opaque: tuple[tuple[int, str], ...] = (),
) -> dict[str, Any]:
    return {
        "graph": {"exact_graph_nodes": nodes},
        "frontier": frontier,
        "expected_status": "observed",
        "expected_evidence": {
            "integral_node": frontier["node_id"],
            "integrand_node": frontier["integrand_node_id"],
            "binder_name": frontier["binder_name"],
            "unshadowed_binder_nodes": sorted(unshadowed),
            "shadowed_binder_nodes": sorted(shadowed),
            "same_name_node_conflict": conflict,
            "unit_interval_bounds_verified": True,
            "visited_nodes": visited,
        },
        "expected_scan_complete": complete,
        "expected_opaque": list(opaque),
    }


def _mismatch(nodes: list[dict[str, Any]], frontier: dict[str, Any]) -> dict[str, Any]:
    return {
        "graph": {"exact_graph_nodes": nodes},
        "frontier": frontier,
        "expected_status": "mismatch",
        "expected_evidence": None,
        "expected_scan_complete": False,
        "expected_opaque": None,
    }


def _integrand_case(name: str) -> dict[str, Any]:
    if name == "product_single_binder":
        nodes, zero, one = _base("product-single")
        binder = _binder(nodes, "product-single")
        coefficient = _add(nodes, "exact_rational", numerator=3, denominator=2)
        body = _add(nodes, "product", arguments=[binder, coefficient])
        integral = _integral(nodes, body, zero, one)
        return _observed(nodes, _frontier(integral, body), unshadowed=[binder],
                         shadowed=[], conflict=False, visited=3, complete=True)

    if name == "source_record_unit_interval":
        nodes, zero, one = _base("source-record")
        binder = _binder(nodes, "source-record")
        body = _add(nodes, "sum", arguments=[binder])
        integral = _integral(
            nodes, body, zero, one,
            source_record={"name": "fresh_moment_t", "binder_node": binder},
        )
        return _observed(nodes, _frontier(integral, body), unshadowed=[binder],
                         shadowed=[], conflict=False, visited=2, complete=True)

    if name == "same_name_handles_conflict":
        nodes, zero, one = _base("same-name")
        first, second = _binder(nodes, "first"), _binder(nodes, "second")
        body = _add(nodes, "product", arguments=[first, second])
        integral = _integral(nodes, body, zero, one)
        return _observed(nodes, _frontier(integral, body),
                         unshadowed=[first, second], shadowed=[], conflict=True,
                         visited=3, complete=False)

    if name == "frontier_name_mismatch":
        nodes, zero, one = _base("name-mismatch")
        body = _add(nodes, "product", arguments=[_binder(nodes, "name")])
        integral = _integral(nodes, body, zero, one)
        frontier = _frontier(integral, body, binder_name="stale_moment_t")
        return _mismatch(nodes, frontier)

    if name == "frontier_reference_mismatch":
        nodes, zero, one = _base("reference-mismatch")
        body = _add(nodes, "product", arguments=[_binder(nodes, "live")])
        decoy = _binder(nodes, "decoy")
        integral = _integral(nodes, body, zero, one)
        return _mismatch(nodes, _frontier(integral, decoy))

    if name == "wrong_integral_operation":
        nodes, zero, one = _base("wrong-integral")
        body = _add(nodes, "product", arguments=[_binder(nodes, "wrong")])
        integral = _integral(nodes, body, zero, one, operation="indefinite_integral")
        return _mismatch(nodes, _frontier(integral, body))

    if name == "invalid_bound_reference":
        nodes, _zero, one = _base("invalid-bound")
        body = _add(nodes, "product", arguments=[_binder(nodes, "bound")])
        integral = _integral(nodes, body, 99_101, one)
        return _mismatch(nodes, _frontier(integral, body))

    if name == "nested_substitution_shadow_value_outer":
        nodes, zero, one = _base("nested-shadow")
        outer, inner = _binder(nodes, "outer"), _binder(nodes, "inner")
        nested = _add(nodes, "function_substitution", expression=inner,
                      variable=inner, value=outer)
        body = _add(nodes, "product", arguments=[nested])
        integral = _integral(nodes, body, zero, one)
        return _observed(nodes, _frontier(integral, body), unshadowed=[outer],
                         shadowed=[inner], conflict=True, visited=4, complete=False)

    if name == "opaque_unknown_lazy_nested_integral":
        nodes, zero, one = _base("opaque-frontier")
        unknown = _add(nodes, "fresh_unknown_source_operation",
                       active_body=99_201, inactive_body=99_202)
        lazy = _add(nodes, "original_lazy_flat_branch", active_body=99_203,
                    flat_predicate=99_204, flat_value=99_205)
        nested_integral = _add(nodes, "definite_integral", variable="fresh_moment_t",
                               integrand=99_206, lower=99_207, upper=99_208)
        body = _add(nodes, "product", arguments=[unknown, lazy, nested_integral])
        integral = _integral(nodes, body, zero, one)
        return _observed(
            nodes, _frontier(integral, body), unshadowed=[], shadowed=[],
            conflict=False, visited=4, complete=False,
            opaque=((unknown, "fresh_unknown_source_operation"),
                    (lazy, "original_lazy_flat_branch"),
                    (nested_integral, "definite_integral")),
        )

    if name == "shadow_binding_not_miscounted":
        nodes, zero, one = _base("shadow-count")
        outer, inner = _binder(nodes, "outer"), _binder(nodes, "inner")
        shadowed_body = _add(nodes, "product", arguments=[inner, inner])
        nested = _add(nodes, "function_substitution", expression=shadowed_body,
                      variable=inner, value=outer)
        body = _add(nodes, "product", arguments=[nested, outer])
        integral = _integral(nodes, body, zero, one)
        return _observed(nodes, _frontier(integral, body), unshadowed=[outer],
                         shadowed=[inner], conflict=True, visited=5, complete=False)

    raise ValueError("unknown integrand fixture case: " + name)


def _scope_ok(result: Any, complete: bool, extra: tuple[str, ...] = ()) -> bool:
    if not isinstance(result, Mapping) or not isinstance(result.get("scope"), Mapping):
        return False
    scope = result["scope"]
    return all(scope.get(key) is False for key in _INTEGRAND_SCOPE + extra) and (
        scope.get("integrand_scan_complete") is complete
    )


def _evidence_ok(result: Any, expected: Mapping[str, Any]) -> bool:
    if not isinstance(result, Mapping) or not isinstance(result.get("evidence"), Mapping):
        return False
    evidence = result["evidence"]
    keys = (
        "integral_node", "integrand_node", "binder_name",
        "unshadowed_binder_nodes", "shadowed_binder_nodes",
        "same_name_node_conflict", "unit_interval_bounds_verified", "visited_nodes",
    )
    return all(evidence.get(key) == expected.get(key) for key in keys)


def _frontier_ok(result: Any, expected: list[tuple[int, str]] | None) -> bool:
    if expected is None:
        return True
    actual = result.get("frontier") if isinstance(result, Mapping) else None
    if not isinstance(actual, list) or len(actual) != len(expected):
        return False
    return all(
        isinstance(item, Mapping)
        and item.get("node") == node
        and item.get("operation") == operation
        and isinstance(item.get("reason"), str)
        and bool(item.get("reason"))
        for item, (node, operation) in zip(actual, expected)
    )


def _receipt(case: Mapping[str, Any], result: Any, exc: str | None) -> dict[str, Any]:
    status = result.get("status") if isinstance(result, Mapping) else None
    checks = {
        "status_enum": status in _VALID_STATUSES,
        "status": status == case["expected_status"],
        "scope": _scope_ok(result, case["expected_scan_complete"]),
    }
    if case["expected_status"] == "observed":
        checks["evidence"] = _evidence_ok(result, case["expected_evidence"])
        checks["frontier"] = _frontier_ok(result, case["expected_opaque"])
    else:
        checks["no_positive_evidence"] = not (
            isinstance(result, Mapping)
            and (
                isinstance(result.get("evidence"), Mapping)
                or isinstance(result.get("source_ref"), Mapping)
            )
        )
    item: dict[str, Any] = {
        "case": case["name"],
        "expected_status": case["expected_status"],
        "status": status,
        "passed": all(checks.values()),
        "checks": checks,
    }
    if exc is not None:
        item["exception"] = exc
    return item


def run_fixture(inspect: Callable[..., Any]) -> dict[str, Any]:
    """Run ten fresh ``inspect(graph, frontier)`` cases."""

    names = (
        "product_single_binder", "source_record_unit_interval",
        "same_name_handles_conflict", "frontier_name_mismatch",
        "frontier_reference_mismatch", "wrong_integral_operation",
        "invalid_bound_reference", "nested_substitution_shadow_value_outer",
        "opaque_unknown_lazy_nested_integral", "shadow_binding_not_miscounted",
    )
    if not callable(inspect):
        return {"schema": SCHEMA, "api": "inspect(graph, frontier)", "passed": False,
                "case_count": 0, "test_count": 1, "checks": {"callable": False},
                "diagnostics": [{"case": "input", "passed": False}],
                "synthetic_only": True, "source_code_executed": False,
                "source_graph_walked": False,
                "scope": {"integral_binder_identity_proved": False,
                          "integrand_scan_complete": False,
                          "numeric_field_value_computed": False,
                          "old_errors_removed": False}}

    diagnostics: list[dict[str, Any]] = []
    checks: dict[str, bool] = {"inspect_callable": True}
    for name in names:
        case = _integrand_case(name)
        result: Any = None
        exc: str | None = None
        try:
            result = inspect(case["graph"], case["frontier"])
        except Exception as error:
            exc = type(error).__name__
        item = _receipt({**case, "name": name}, result, exc)
        diagnostics.append(item)
        checks.update({f"{name}_{key}": value for key, value in item["checks"].items()})
    return {
        "schema": SCHEMA, "api": "inspect(graph, frontier)",
        "passed": bool(checks) and all(checks.values()), "checks": checks,
        "diagnostics": diagnostics, "case_count": len(names),
        "test_count": len(checks), "synthetic_only": True,
        "source_code_executed": False, "source_graph_walked": False,
        "numeric_values_computed": False,
        "scope": {"integral_binder_identity_proved": False,
                  "integrand_scan_complete": False,
                  "numeric_field_value_computed": False,
                  "old_errors_removed": False},
    }


# ---------------------------------------------------------------------------
# Independent synthetic fixture for the adjacent analytic-core source ref.

_CORE_SCOPE_FLAGS = ("numeric_field_value_computed", "core_solution_admission_proved")
_CORE_QUANTITIES = frozenset({"Phi", "Psi", "mean_Psi"})


def _core_graph(case: str, *, quantity: Any = "Phi", tamper: str | None = None) -> tuple[dict[str, Any], dict[str, Any], dict[str, Any]]:
    nodes: list[dict[str, Any]] = [
        {"operation": f"{case}_unreachable_dummy"},
        {"operation": "bound_variable", "name": "fresh_core_rho"},
        {"operation": "bound_variable", "name": "fresh_core_Z"},
    ]
    row: dict[str, Any] = {
        "operation": "current_original_unique_analytic_core_component",
        "quantity": quantity,
        "radial_order": 2,
        "ordinary_slow_Z_derivative_order": 1,
        "rho": 1,
        "Z_variable": 2,
        "actual_twenty_term_operator": True,
        "unique_analytic_Banach_solution": True,
        "admitted_solution_not_a_model_or_enclosure_center": True,
        "source_family": "fresh_original_core",
        "component_definition": "fresh analytic core component definition",
        "source_provider": "fresh-core-source-provider",
        "radial_recurrence_provider": "fresh-core-radial-provider",
        "actual_source_provider": "fresh-core-actual-source-provider",
        "actual_green_inverse_provider": "fresh-core-green-provider",
    }
    if tamper == "flag":
        row["actual_twenty_term_operator"] = False
    elif tamper == "rho":
        row["rho"] = 99_301
    elif tamper == "Z":
        row["Z_variable"] = 99_302
    elif tamper == "order":
        row["radial_order"] = -1
    nodes.append(row)
    source = len(nodes) - 1
    frontier = {
        "node_id": source,
        "operation": "current_original_unique_analytic_core_component",
        "expected_binder_node": 1,
    }
    return {"exact_graph_nodes": nodes}, frontier, {"source": source, "rho": 1, "Z_variable": 2, "quantity": quantity}


def _core_scope_ok(result: Any) -> bool:
    if not isinstance(result, Mapping) or not isinstance(result.get("scope"), Mapping):
        return False
    scope = result["scope"]
    return all(scope.get(key) is False for key in _CORE_SCOPE_FLAGS)


def run_core_fixture(inspect_core: Callable[..., Any]) -> dict[str, Any]:
    """Run core source-reference admission checks without numeric evaluation."""

    cases = (
        ("positive_Phi", "Phi", None),
        ("positive_Psi", "Psi", None),
        ("positive_mean_Psi", "mean_Psi", None),
        ("wrong_quantity", "not_a_core_quantity", None),
        ("tampered_flag", "Phi", "flag"),
        ("invalid_rho", "Phi", "rho"),
        ("invalid_Z", "Phi", "Z"),
        ("illegal_order", "Phi", "order"),
    )
    if not callable(inspect_core):
        return {"schema": SCHEMA, "api": "inspect_core(graph, frontier)", "passed": False,
                "case_count": 0, "test_count": 1, "checks": {"callable": False},
                "diagnostics": [{"case": "input", "passed": False}]}

    checks: dict[str, bool] = {"inspect_core_callable": True}
    diagnostics: list[dict[str, Any]] = []
    for name, quantity, tamper in cases:
        graph, frontier, expected = _core_graph(name, quantity=quantity, tamper=tamper)
        result: Any = None
        error: str | None = None
        try:
            result = inspect_core(graph, frontier)
        except Exception as exc:
            error = type(exc).__name__
        status = result.get("status") if isinstance(result, Mapping) else None
        positive = tamper is None and quantity in _CORE_QUANTITIES
        local: dict[str, bool] = {
            "status": status == ("source_ref" if positive else "mismatch"),
            "scope": _core_scope_ok(result),
        }
        if positive:
            source_ref = result.get("source_ref") if isinstance(result, Mapping) else None
            evidence = result.get("evidence") if isinstance(result, Mapping) else None
            retained = source_ref if isinstance(source_ref, Mapping) else evidence
            rho_ref = retained.get("rho") if isinstance(retained, Mapping) else None
            if rho_ref is None and isinstance(retained, Mapping):
                rho_ref = retained.get("rho_node")
            z_ref = retained.get("Z_variable") if isinstance(retained, Mapping) else None
            if z_ref is None and isinstance(retained, Mapping):
                z_ref = retained.get("Z_node")
            local["source_ref_shape"] = (
                isinstance(retained, Mapping)
                and retained.get("quantity") == quantity
                and rho_ref == expected["rho"]
                and z_ref == expected["Z_variable"]
            )
        else:
            local["no_positive_ref"] = not (
                isinstance(result, Mapping)
                and (
                    isinstance(result.get("source_ref"), Mapping)
                    or isinstance(result.get("evidence"), Mapping)
                )
            )
        item: dict[str, Any] = {"case": name, "status": status,
                                "passed": all(local.values()), "checks": local}
        if error is not None:
            item["exception"] = error
        diagnostics.append(item)
        checks.update({f"{name}_{key}": value for key, value in local.items()})
    return {
        "schema": SCHEMA, "api": "inspect_core(graph, frontier)",
        "passed": bool(checks) and all(checks.values()), "checks": checks,
        "diagnostics": diagnostics, "case_count": len(cases),
        "test_count": len(checks), "synthetic_only": True,
        "source_code_executed": False, "source_graph_walked": False,
        "scope": {"numeric_field_value_computed": False,
                  "core_solution_admission_proved": False},
    }


__all__ = ["SCHEMA", "run_fixture", "run_core_fixture"]
