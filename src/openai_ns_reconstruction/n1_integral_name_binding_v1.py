"""Conditional source-semantic name binding for the saved N1 integrals.

This leaf consumes two caller-owned values only: a decoded scoped-catalog
result and a mapping of raw Python source texts.  It parses those texts with
``ast`` and checks a small, fixed set of defining statements.  It never reads
the graph, walks graph edges, imports a producer, evaluates a source, or
installs a binder mapping.

The result deliberately says ``source_name_binding_candidate``.  A source
string is not a graph-node identity witness, so even a passing candidate keeps
actual owner identity, callback certification, and all numerical/physical
claims false.  Saved occurrence counts are copied as receipt metadata only;
their paths are not re-walked here.
"""

from __future__ import annotations

import ast
from collections.abc import Mapping
from typing import Any


SCHEMA = "n1_integral_name_binding_v1"
RESULT_SCHEMA = SCHEMA + "_result"
FIXTURE_SCHEMA = "n1_integral_name_binding_fixture_v1"
EXPECTED_COUNT = 7

CORE_VARIABLE = "original_true_core_gauge_v"
LOWER_VARIABLE = "original_lower_collar_ratio_v"
DERIVATIVE_MARKERS = frozenset(
    {
        "derivative_of_same_defining_integral",
        "ordinary_Z_derivative_order",
        "differentiation_under_integral",
    }
)


def _clone(value: Any) -> Any:
    """Copy JSON-like receipt data without invoking arbitrary objects."""

    if isinstance(value, Mapping):
        return {str(key): _clone(item) for key, item in value.items()}
    if isinstance(value, list):
        return [_clone(item) for item in value]
    if isinstance(value, tuple):
        return [_clone(item) for item in value]
    return value


def _is_int(value: Any) -> bool:
    return type(value) is int


def _loc(node: ast.AST | None) -> dict[str, int] | None:
    if node is None or not hasattr(node, "lineno"):
        return None
    return {
        "line": int(getattr(node, "lineno")),
        "column": int(getattr(node, "col_offset", 0)),
    }


def _name(node: ast.AST | None, wanted: str | None = None) -> bool:
    return isinstance(node, ast.Name) and (wanted is None or node.id == wanted)


def _const(node: ast.AST | None, value: Any) -> bool:
    return isinstance(node, ast.Constant) and node.value == value


def _attribute(node: ast.AST | None, owner: str, attribute: str) -> bool:
    return (
        isinstance(node, ast.Attribute)
        and node.attr == attribute
        and _name(node.value, owner)
    )


def _call_attribute(
    node: ast.AST | None,
    owner: str,
    attribute: str,
) -> ast.Call | None:
    if not isinstance(node, ast.Call) or not _attribute(node.func, owner, attribute):
        return None
    return node


def _functions(tree: ast.AST, name: str | None = None) -> list[ast.FunctionDef | ast.AsyncFunctionDef]:
    rows = [
        item
        for item in ast.walk(tree)
        if isinstance(item, (ast.FunctionDef, ast.AsyncFunctionDef))
    ]
    return [item for item in rows if name is None or item.name == name]


def _classes(tree: ast.AST, name: str) -> list[ast.ClassDef]:
    return [item for item in ast.walk(tree) if isinstance(item, ast.ClassDef) and item.name == name]


def _method(class_node: ast.ClassDef, name: str) -> list[ast.FunctionDef | ast.AsyncFunctionDef]:
    return [
        item
        for item in class_node.body
        if isinstance(item, (ast.FunctionDef, ast.AsyncFunctionDef)) and item.name == name
    ]


def _argument_names(function: ast.FunctionDef | ast.AsyncFunctionDef) -> set[str]:
    arguments = function.args
    result = {item.arg for item in arguments.posonlyargs}
    result.update(item.arg for item in arguments.args)
    result.update(item.arg for item in arguments.kwonlyargs)
    if arguments.vararg is not None:
        result.add(arguments.vararg.arg)
    if arguments.kwarg is not None:
        result.add(arguments.kwarg.arg)
    return result


def _slice_constant(node: ast.Subscript) -> Any:
    # Python 3.8 exposed ast.Index while current Python exposes the value
    # directly.  Keeping this helper local makes the AST contract explicit.
    value = node.slice
    if isinstance(value, ast.Index):  # pragma: no cover - old Python shape
        value = value.value
    return value.value if isinstance(value, ast.Constant) else None


def _g_nodes_name(node: ast.AST | None, variable: str) -> bool:
    """Match g.nodes[var.node]['name'] exactly, including the handle name."""

    if not isinstance(node, ast.Subscript) or _slice_constant(node) != "name":
        return False
    outer = node.value
    if not isinstance(outer, ast.Subscript):
        return False
    inner = outer.value
    if not _attribute(inner, "g", "nodes"):
        return False
    return _attribute(outer.slice, variable, "node")


def _symbol_call(
    node: ast.AST | None,
    expected: str | None = None,
    *,
    owner: str = "g",
) -> ast.Call | None:
    call = _call_attribute(node, owner, "symbol")
    if call is None or len(call.args) != 1 or not isinstance(call.args[0], ast.Constant):
        return None
    if not isinstance(call.args[0].value, str):
        return None
    if expected is not None and call.args[0].value != expected:
        return None
    return call


def _assignment_name(node: ast.AST | None) -> str | None:
    if not isinstance(node, ast.Assign) or len(node.targets) != 1:
        return None
    target = node.targets[0]
    return target.id if isinstance(target, ast.Name) else None


def _result(
    key: str,
    *,
    passed: bool,
    evidence: Mapping[str, Any] | None = None,
    error: str | None = None,
    status: str | None = None,
) -> dict[str, Any]:
    item: dict[str, Any] = {
        "name": key,
        "passed": bool(passed),
        "status": status or ("ast_pattern_observed" if passed else "ast_pattern_mismatch"),
        "evidence_scope": "ast_fragment_pattern_only",
        "data_flow_proved": False,
        "source_call_binding_proved": False,
        "graph_node_identity_proved": False,
        "evidence": _clone(evidence or {}),
    }
    if error is not None:
        item["error"] = error
        if status is None:
            item["status"] = "unresolved"
    return item


def _parse_source(source_text: Any, label: str) -> tuple[ast.Module | None, dict[str, Any]]:
    if not isinstance(source_text, str) or not source_text.strip():
        return None, _result(label, passed=False, error="required raw source text is missing", status="unresolved")
    try:
        tree = ast.parse(source_text, filename=label, mode="exec")
    except (SyntaxError, TypeError, ValueError) as exc:
        return None, _result(label, passed=False, error=type(exc).__name__ + ": " + str(exc), status="unresolved")
    return tree, _result(label, passed=True, evidence={"ast": "Module", "source_text_consumed": True})


def _check_graph(source_text: Any) -> dict[str, Any]:
    tree, parsed = _parse_source(source_text, "graph")
    if tree is None:
        return parsed
    classes = _classes(tree, "FunctionTransportGraph")
    if len(classes) != 1:
        return _result("graph", passed=False, error="FunctionTransportGraph class is missing or ambiguous", status="unresolved")
    graph = classes[0]
    node_methods = _method(graph, "node")
    symbol_methods = _method(graph, "symbol")
    if len(node_methods) != 1 or len(symbol_methods) != 1:
        return _result("graph", passed=False, error="node/symbol method is missing or ambiguous", status="unresolved")
    node_method, symbol_method = node_methods[0], symbol_methods[0]

    symbol_hits: list[ast.Call] = []
    for item in ast.walk(symbol_method):
        if isinstance(item, ast.Return):
            call = _symbol_call(item.value, "bound_variable", owner="self")
            if call is not None and any(
                keyword.arg == "name" and _name(keyword.value, "name")
                for keyword in call.keywords
            ):
                symbol_hits.append(call)
    symbol_ok = len(symbol_hits) == 1 and "name" in _argument_names(symbol_method)

    row_ok = False
    key_ok = False
    intern_ok = False
    return_ok = False
    row_line: dict[str, int] | None = None
    key_line: dict[str, int] | None = None
    intern_line: dict[str, int] | None = None
    return_line: dict[str, int] | None = None
    for item in ast.walk(node_method):
        if not isinstance(item, ast.Assign) or len(item.targets) != 1:
            continue
        target = item.targets[0]
        if not _name(target, "row") or not isinstance(item.value, ast.Call) or not _name(item.value.func, "dict"):
            continue
        operation_kw = next((kw for kw in item.value.keywords if kw.arg == "operation"), None)
        attributes_kw = next((kw for kw in item.value.keywords if kw.arg is None), None)
        row_ok = _name(operation_kw.value, "operation") if operation_kw is not None else False
        row_ok = row_ok and _name(attributes_kw.value, "attributes") if attributes_kw is not None else False
        if row_ok:
            row_line = _loc(item)
    for item in ast.walk(node_method):
        if not isinstance(item, ast.Assign) or len(item.targets) != 1 or not _name(item.targets[0], "key"):
            continue
        call = item.value
        if not _call_attribute(call, "json", "dumps"):
            continue
        assert isinstance(call, ast.Call)
        sort_kw = next((kw for kw in call.keywords if kw.arg == "sort_keys"), None)
        separators_kw = next((kw for kw in call.keywords if kw.arg == "separators"), None)
        tuple_value = separators_kw.value if separators_kw is not None else None
        tuple_ok = (
            isinstance(tuple_value, (ast.Tuple, ast.List))
            and len(tuple_value.elts) == 2
            and _const(tuple_value.elts[0], ",")
            and _const(tuple_value.elts[1], ":")
        )
        key_ok = (
            len(call.args) == 1
            and _name(call.args[0], "row")
            and sort_kw is not None
            and _const(sort_kw.value, True)
            and separators_kw is not None
            and tuple_ok
        )
        if key_ok:
            key_line = _loc(item)

    for item in ast.walk(node_method):
        if not isinstance(item, ast.If):
            continue
        test = item.test
        compare_ok = (
            isinstance(test, ast.Compare)
            and len(test.ops) == 1
            and isinstance(test.ops[0], ast.NotIn)
            and _name(test.left, "key")
            and len(test.comparators) == 1
            and _attribute(test.comparators[0], "self", "keys")
        )
        if not compare_ok:
            continue
        has_index = False
        has_append = False
        for child in ast.walk(ast.Module(body=item.body, type_ignores=[])):
            if isinstance(child, ast.Assign) and len(child.targets) == 1:
                target = child.targets[0]
                if isinstance(target, ast.Subscript) and _attribute(target.value, "self", "keys") and _slice_constant(target) == "key":
                    value = child.value
                    has_index = (
                        isinstance(value, ast.Call)
                        and _name(value.func, "len")
                        and len(value.args) == 1
                        and _attribute(value.args[0], "self", "nodes")
                    )
            if isinstance(child, ast.Expr) and isinstance(child.value, ast.Call):
                call = child.value
                if (
                    isinstance(call.func, ast.Attribute)
                    and call.func.attr == "append"
                    and _attribute(call.func.value, "self", "nodes")
                ):
                    has_append = len(call.args) == 1 and _name(call.args[0], "row")
        if has_index and has_append:
            intern_ok = True
            intern_line = _loc(item)

    for item in ast.walk(node_method):
        if not isinstance(item, ast.Return) or not isinstance(item.value, ast.Call):
            continue
        call = item.value
        if not _name(call.func, "FunctionRef") or len(call.args) != 2:
            continue
        if not _name(call.args[0], "self"):
            continue
        handle = call.args[1]
        location = _loc(item)
        if (
            isinstance(handle, ast.Subscript)
            and _attribute(handle.value, "self", "keys")
            and _slice_constant(handle) == "key"
            and intern_line is not None
            and location is not None
            and location["line"] > intern_line["line"]
        ):
            return_ok = True
            return_line = location

    passed = symbol_ok and row_ok and key_ok and intern_ok and return_ok
    return _result(
        "graph",
        passed=passed,
        evidence={
            "class": "FunctionTransportGraph",
            "symbol": {
                "calls_bound_variable": symbol_ok,
                "location": _loc(symbol_hits[0]) if symbol_hits else None,
            },
            "node": {
                "row_dict": row_ok,
                "canonical_json_key": key_ok,
                "keys_interning": intern_ok,
                "returns_interned_function_ref": return_ok,
                "row_location": row_line,
                "key_location": key_line,
                "intern_location": intern_line,
                "return_location": return_line,
            },
        },
        error=None if passed else "fixed FunctionTransportGraph symbol/node/interned-handle contract is not fully present",
    )


def _integral_calls(function: ast.AST) -> list[ast.Call]:
    return [
        item
        for item in ast.walk(function)
        if isinstance(item, ast.Call)
        and _call_attribute(item, "g", "node") is not None
        and item.args
        and _const(item.args[0], "definite_integral")
    ]


def _check_controls(source_text: Any) -> dict[str, Any]:
    tree, parsed = _parse_source(source_text, "controls")
    if tree is None:
        return parsed
    functions = _functions(tree, "integral")
    if len(functions) != 1:
        return _result("controls_integral", passed=False, error="controls.integral function is missing or ambiguous", status="unresolved")
    function = functions[0]
    names = _argument_names(function)
    hits: list[ast.Call] = []
    direct_variable: list[ast.Call] = []
    for call in _integral_calls(function):
        if "variable" not in names:
            continue
        hits.append(call)
        keyword = next((kw for kw in call.keywords if kw.arg == "variable"), None)
        if keyword is not None and _name(keyword.value, "variable"):
            direct_variable.append(call)
    passed = "variable" in names and len(hits) >= 1 and len(direct_variable) >= 1
    return _result(
        "controls_integral",
        passed=passed,
        evidence={
            "function": function.name,
            "variable_parameter": "variable" in names,
            "definite_integral_calls": len(hits),
            "variable_saved_as_string_argument": bool(direct_variable),
            "location": _loc(direct_variable[0]) if direct_variable else None,
        },
        error=None if passed else "controls.integral does not pass variable directly to the definite-integral row",
    )


def _controls_integral_call(node: ast.AST) -> bool:
    return (
        isinstance(node, ast.Call)
        and _call_attribute(node, "controls", "integral") is not None
    )


def _symbol_assignments(function: ast.AST, expected: str) -> list[tuple[str, ast.Assign, ast.Call]]:
    hits: list[tuple[str, ast.Assign, ast.Call]] = []
    for item in ast.walk(function):
        variable = _assignment_name(item)
        if variable is None:
            continue
        call = _symbol_call(item.value, expected)
        if call is not None:
            hits.append((variable, item, call))
    return hits


def _check_core(source_text: Any) -> dict[str, Any]:
    tree, parsed = _parse_source(source_text, "core")
    if tree is None:
        return parsed
    evidence: list[dict[str, Any]] = []
    for function in _functions(tree):
        assignments = _symbol_assignments(function, CORE_VARIABLE)
        for handle, assignment, call in assignments:
            calls = [item for item in ast.walk(function) if _controls_integral_call(item)]
            lookup_calls = [
                item
                for item in calls
                if hasattr(item, "lineno")
                and hasattr(assignment, "lineno")
                and item.lineno > assignment.lineno
                if any(_g_nodes_name(argument, handle) for argument in item.args)
                or any(_g_nodes_name(keyword.value, handle) for keyword in item.keywords)
            ]
            for item in lookup_calls:
                evidence.append(
                    {
                        "function": function.name,
                        "handle": handle,
                        "symbol_name": CORE_VARIABLE,
                        "symbol_location": _loc(call),
                        "integral_location": _loc(item),
                        "g_nodes_name_lookup": True,
                        "same_function_lexical_order": True,
                    }
                )
    passed = len(evidence) == 1
    return _result(
        "core_constructor",
        passed=passed,
        evidence={
            "target_name": CORE_VARIABLE if evidence else None,
            "symbol_to_controls_integral": evidence,
            "constructor_scope_found": bool(evidence),
        },
        error=None if passed else "core gauge symbol and controls.integral name lookup are not uniquely connected",
    )


def _function_call_name(call: ast.Call) -> str | None:
    if isinstance(call.func, ast.Name):
        return call.func.id
    if isinstance(call.func, ast.Attribute):
        return call.func.attr
    return None


def _check_lower(source_text: Any) -> dict[str, Any]:
    tree, parsed = _parse_source(source_text, "lower_collar")
    if tree is None:
        return parsed
    helper_evidence: list[dict[str, Any]] = []
    for helper in _functions(tree, "integral"):
        names = _argument_names(helper)
        variable_names = [name for name in names if name not in {"g", "jet", "lower", "upper", "metadata"}]
        for variable in variable_names:
            for call in ast.walk(helper):
                if not _controls_integral_call(call):
                    continue
                assert isinstance(call, ast.Call)
                if any(_g_nodes_name(argument, variable) for argument in call.args) or any(
                    _g_nodes_name(keyword.value, variable) for keyword in call.keywords
                ):
                    helper_evidence.append(
                        {
                            "function": helper.name,
                            "variable_parameter": variable,
                            "location": _loc(call),
                        }
                    )

    closure_evidence: list[dict[str, Any]] = []
    for function in _functions(tree):
        assignments = _symbol_assignments(function, LOWER_VARIABLE)
        for handle, assignment, call in assignments:
            helper_calls = [
                item
                for item in ast.walk(function)
                if isinstance(item, ast.Call)
                and _function_call_name(item) == "integral"
                and hasattr(item, "lineno")
                and hasattr(assignment, "lineno")
                and item.lineno > assignment.lineno
                and any(_name(argument, handle) for argument in item.args)
            ]
            for helper_call in helper_calls:
                closure_evidence.append(
                    {
                        "function": function.name,
                        "handle": handle,
                        "symbol_name": LOWER_VARIABLE,
                        "symbol_location": _loc(call),
                        "helper_call_location": _loc(helper_call),
                        "same_function_lexical_order": True,
                    }
                )
    passed = len(helper_evidence) >= 1 and len(closure_evidence) == 1
    return _result(
        "lower_collar_closure",
        passed=passed,
        evidence={
            "target_name": LOWER_VARIABLE if closure_evidence else None,
            "helper_converts_handle_to_name": helper_evidence,
            "closure_passes_symbol_handle": closure_evidence,
        },
        error=None if passed else "lower-collar symbol/helper conversion is not fully present",
    )


def _definite_branch(test: ast.AST) -> bool:
    return (
        isinstance(test, ast.Compare)
        and len(test.ops) == 1
        and isinstance(test.ops[0], ast.Eq)
        and _name(test.left, "op")
        and len(test.comparators) == 1
        and _const(test.comparators[0], "definite_integral")
    )


def _copy_entire_row(value: ast.AST) -> bool:
    if isinstance(value, ast.Call):
        if _attribute(value.func, "copy", "deepcopy"):
            return len(value.args) == 1 and _name(value.args[0], "row")
        if _name(value.func, "dict"):
            return len(value.args) == 1 and _name(value.args[0], "row")
        if _attribute(value.func, "row", "copy"):
            return len(value.args) == 0
    return False


def _subscript_key(node: ast.AST | None, owner: str, key: str) -> bool:
    return isinstance(node, ast.Subscript) and _name(node.value, owner) and _slice_constant(node) == key


def _check_c3(source_text: Any) -> dict[str, Any]:
    tree, parsed = _parse_source(source_text, "c3_derivative")
    if tree is None:
        return parsed
    classes = _classes(tree, "OrdinaryZDerivative")
    functions = [item for item in _functions(tree, "rule") if not classes or any(item in cls.body for cls in classes)]
    branches: list[ast.If] = []
    for function in functions:
        branches.extend(item for item in ast.walk(function) if isinstance(item, ast.If) and _definite_branch(item.test))
    if len(branches) != 1:
        return _result("ordinary_z_derivative", passed=False, error="OrdinaryZDerivative definite_integral rule is missing or ambiguous", status="unresolved")
    branch = branches[0]
    copied = [
        item
        for item in ast.walk(ast.Module(body=branch.body, type_ignores=[]))
        if isinstance(item, ast.Assign)
        and len(item.targets) == 1
        and _name(item.targets[0], "attrs")
        and _copy_entire_row(item.value)
    ]
    changed = [
        item
        for item in ast.walk(ast.Module(body=branch.body, type_ignores=[]))
        if isinstance(item, ast.Assign)
        and len(item.targets) == 1
        and _subscript_key(item.targets[0], "attrs", "integrand")
    ]
    returned = []
    for item in ast.walk(ast.Module(body=branch.body, type_ignores=[])):
        if not isinstance(item, ast.Return) or not isinstance(item.value, ast.Call):
            continue
        call = item.value
        if _call_attribute(call, "g", "node") is None or not call.args or not _const(call.args[0], "definite_integral"):
            continue
        if any(keyword.arg is None and _name(keyword.value, "attrs") for keyword in call.keywords):
            returned.append(item)
    variable_writes = []
    variable_pops = []
    for item in ast.walk(ast.Module(body=branch.body, type_ignores=[])):
        if isinstance(item, ast.Assign) and len(item.targets) == 1 and _subscript_key(item.targets[0], "attrs", "variable"):
            variable_writes.append(item)
        if isinstance(item, ast.Expr) and isinstance(item.value, ast.Call) and _attribute(item.value.func, "attrs", "pop"):
            if item.value.args and _const(item.value.args[0], "variable"):
                variable_pops.append(item)
    lexical_order = (
        bool(copied)
        and bool(changed)
        and bool(returned)
        and copied[0].lineno < changed[0].lineno < returned[0].lineno
    )
    passed = lexical_order and not variable_writes and not variable_pops
    return _result(
        "ordinary_z_derivative",
        passed=passed,
        evidence={
            "rule_location": _loc(branch),
            "copied_entire_attrs": bool(copied),
            "changed_integrand": bool(changed),
            "returned_definite_integral_with_attrs": bool(returned),
            "lexical_order_preserved": lexical_order,
            "variable_preserved_by_copy": bool(copied) and not variable_writes and not variable_pops,
            "variable_overwrites": len(variable_writes) + len(variable_pops),
        },
        error=None if passed else "ordinary-Z definite-integral derivative does not copy attrs/change only integrand/preserve variable",
    )


def _source_value(source_texts: Mapping[str, Any], aliases: tuple[str, ...]) -> Any:
    for key in aliases:
        if key in source_texts:
            value = source_texts[key]
            if isinstance(value, Mapping) and isinstance(value.get("source_text"), str):
                return value["source_text"]
            return value
    return None


def _sources(source_texts: Any) -> dict[str, Any]:
    if not isinstance(source_texts, Mapping):
        return {}
    nested = source_texts.get("source_texts")
    if isinstance(nested, Mapping):
        source_texts = nested
    return {
        "graph": _source_value(source_texts, ("graph", "graph_source", "function_transport_graph", "function_transport_graph_source", "graph_function_transport", "FunctionTransportGraph")),
        "controls": _source_value(source_texts, ("controls", "controls_source", "controls_integral", "integral_controls")),
        "core": _source_value(source_texts, ("core", "core_source", "core_constructor", "core_source_constructor")),
        "lower_collar": _source_value(source_texts, ("lower_collar", "lower", "lower_source", "lower_collar_source")),
        "c3_derivative": _source_value(source_texts, ("c3", "c3_source", "c3_derivative_source", "derivative", "derivative_source", "ordinary_z_derivative")),
    }


def _catalog(saved_result: Any) -> tuple[list[Any], str]:
    if isinstance(saved_result, list):
        return list(saved_result), "direct_catalog_list"
    if not isinstance(saved_result, Mapping):
        return [], "missing"
    for key in ("integral_candidates", "integral_records", "scoped_catalog", "catalog"):
        value = saved_result.get(key)
        if isinstance(value, list):
            return list(value), key
    progress = saved_result.get("progress")
    if isinstance(progress, Mapping):
        for key in ("integral_records", "integral_candidates", "scoped_catalog", "catalog"):
            value = progress.get(key)
            if isinstance(value, list):
                return list(value), "progress." + key
    return [], "missing"


def _occurrence_records(saved_result: Any) -> dict[int, Mapping[str, Any]]:
    holders: list[Any] = []
    if isinstance(saved_result, Mapping):
        holders.extend(
            [
                saved_result.get("occurrence_output"),
                saved_result.get("saved_occurrence_output"),
                saved_result.get("leaf_raw_output"),
                saved_result,
            ]
        )
        progress = saved_result.get("progress")
        if isinstance(progress, Mapping):
            holders.extend([progress.get("leaf_raw_output"), progress.get("leaf_output"), progress])
    index: dict[int, Mapping[str, Any]] = {}
    for holder in holders:
        if not isinstance(holder, Mapping):
            continue
        rows = holder.get("records")
        if not isinstance(rows, list):
            continue
        for row in rows:
            if isinstance(row, Mapping) and _is_int(row.get("integral_node_id")):
                index.setdefault(row["integral_node_id"], row)
    return index


def _saved_occurrence_summary(row: Mapping[str, Any] | None) -> dict[str, Any]:
    if row is None:
        return {
            "supplied": False,
            "witness_count": None,
            "opaque_frontier_count": None,
            "paths_rewalked": False,
        }
    occurrence_rows = row.get("unshadowed_occurrences")
    if not isinstance(occurrence_rows, list):
        occurrence_rows = row.get("occurrences")
    frontier_rows = row.get("opaque_frontier")
    if not isinstance(frontier_rows, list):
        frontier_rows = row.get("frontier")
    summary: dict[str, Any] = {
        "supplied": True,
        "witness_count": len(occurrence_rows) if isinstance(occurrence_rows, list) else None,
        "opaque_frontier_count": len(frontier_rows) if isinstance(frontier_rows, list) else None,
        "paths_rewalked": False,
        "saved_status": row.get("status"),
    }
    if isinstance(row.get("errors"), list):
        summary["saved_errors"] = _clone(row.get("errors"))
    return summary


def _row_record(
    row: Any,
    *,
    source_checks: Mapping[str, Mapping[str, Any]],
    occurrence: Mapping[str, Any] | None,
) -> dict[str, Any]:
    record: dict[str, Any] = {
        "catalog_record": _clone(row),
        "source_name_binding_candidate": False,
        "candidate_basis": "ast_fragment_pattern_only_not_source_call_proof",
        "candidate_actual_integral_binder": False,
        "candidate_actual_integral_binder_proved": False,
        "actual_graph_owner_identity": False,
        "binder_identity_proved": False,
        "parent_identity_proved": False,
        "parent_identity_verified": False,
        "binder_mapping_installed": False,
        "source_callback_certified": False,
        "source_function_equality_proved": False,
        "source_call_binding_proved": False,
        "source_data_flow_proved": False,
        "ast_fragment_pattern_only": True,
        "m15_claimed": False,
        "parent_identity_unproved": True,
        "saved_occurrence": _saved_occurrence_summary(occurrence),
        "unresolved_reasons": [],
    }
    if not isinstance(row, Mapping):
        record["unresolved_reasons"].append("catalog record is not an object")
        record["status"] = "unresolved"
        return record
    integral_id = row.get("integral_node_id")
    candidate_id = row.get("unique_node_id")
    variable = row.get("variable_name")
    original = row.get("original_record")
    record.update(
        {
            "integral_node_id": integral_id,
            "unique_node_id": candidate_id,
            "variable_name": variable,
            "candidate_node_ids": _clone(row.get("candidate_node_ids")),
            "original_record": _clone(original),
        }
    )
    valid = True
    if not _is_int(integral_id) or integral_id < 0:
        valid = False
        record["unresolved_reasons"].append("integral_node_id is not a nonnegative integer")
    if not _is_int(candidate_id) or candidate_id < 0:
        valid = False
        record["unresolved_reasons"].append("unique_node_id is not a nonnegative integer")
    candidate_ids = row.get("candidate_node_ids")
    if not isinstance(candidate_ids, list) or len(candidate_ids) != 1 or candidate_ids[0] != candidate_id:
        valid = False
        record["unresolved_reasons"].append("candidate_node_ids is not the unique candidate")
    if row.get("binding_status") != "unique_candidate":
        valid = False
        record["unresolved_reasons"].append("binding_status is not unique_candidate")
    if not isinstance(variable, str) or not variable:
        valid = False
        record["unresolved_reasons"].append("variable_name is missing")
    if not isinstance(original, Mapping) or original.get("operation") != "definite_integral":
        valid = False
        record["unresolved_reasons"].append("original_record is not a definite_integral record")
    else:
        if original.get("variable") != variable:
            valid = False
            record["unresolved_reasons"].append("catalog/original variable mismatch")
        if not _is_int(original.get("integrand")):
            valid = False
            record["unresolved_reasons"].append("original integrand is not a node id")
        explicit_parent = original.get("parent_node_id", row.get("parent_node_id"))
        if _is_int(explicit_parent) and explicit_parent >= 0:
            record["parent_node_id"] = explicit_parent
            record["parent_node_id_unverified_reference"] = True
        # In particular, upper/lower metadata is never promoted to a parent
        # identity.  No source path is independently checked by this leaf, so
        # parent_identity_unproved stays true even when a caller supplies a
        # parent_node_id-looking integer.  The three derivative rows have no
        # parent node id in the saved catalog, but the same conservative rule
        # applies to every row.
        record["parent_identity_unproved"] = True

    if not valid:
        record["status"] = "unresolved"
        return record

    graph = source_checks.get("graph", {})
    graph_return = source_checks.get("graph_return_handle", {})
    controls = source_checks.get("controls_integral", {})
    required = ["graph", "graph_return_handle", "controls_integral"]
    core = variable == CORE_VARIABLE
    lower = variable == LOWER_VARIABLE
    derivative = isinstance(original, Mapping) and any(marker in original for marker in DERIVATIVE_MARKERS)
    if core:
        required.append("core_constructor")
    elif lower:
        required.append("lower_collar_closure")
        if derivative:
            required.append("ordinary_z_derivative")
    else:
        record["unresolved_reasons"].append("source variable is outside the fixed core/lower contracts")
    if not graph.get("passed"):
        record["unresolved_reasons"].append("graph symbol/node AST evidence is unresolved")
    if not graph_return.get("passed"):
        record["unresolved_reasons"].append("graph node return-handle AST evidence is unresolved")
    if not controls.get("passed"):
        record["unresolved_reasons"].append("controls.integral variable-string AST evidence is unresolved")
    for check in required[3:]:
        if not source_checks.get(check, {}).get("passed"):
            record["unresolved_reasons"].append(check + " AST evidence is unresolved")
    candidate = valid and (core or lower) and all(source_checks.get(check, {}).get("passed") is True for check in required)
    record["source_name_binding_candidate"] = bool(candidate)
    record["required_source_checks"] = required
    record["source_role"] = "core" if core else "lower" if lower else "unknown"
    record["derivative_rule_required"] = bool(derivative and lower)
    record["status"] = "conditional_candidate" if candidate else "unresolved"
    return record


def _looks_like_sources(value: Any) -> bool:
    if not isinstance(value, Mapping):
        return False
    keys = {str(key) for key in value}
    return bool(
        keys
        & {
            "graph",
            "graph_source",
            "function_transport_graph",
            "function_transport_graph_source",
            "controls",
            "controls_integral",
            "core",
            "core_source",
            "lower_collar",
            "lower_collar_source",
            "c3",
            "c3_source",
            "c3_derivative_source",
            "source_texts",
        }
    ) and not bool(
        keys & {"integral_candidates", "integral_records", "progress", "original_record"}
    )


def analyze(saved_result: Any, source_texts: Any = None) -> dict[str, Any]:
    """Return conditional source-name candidates for a saved seven-row catalog."""

    # Accept either natural positional order without changing the public
    # contract: analyze(saved_result, source_texts).  This also helps a caller
    # that naturally writes analyze(source_texts, saved_result).
    if _looks_like_sources(saved_result) and not _looks_like_sources(source_texts):
        saved_result, source_texts = source_texts, saved_result
    if source_texts is None and isinstance(saved_result, Mapping):
        source_texts = saved_result.get("source_texts")

    sources = _sources(source_texts)
    graph_check = _check_graph(sources.get("graph"))
    graph_return_ok = bool(
        graph_check.get("evidence", {})
        .get("node", {})
        .get("returns_interned_function_ref")
        is True
    )
    source_checks = {
        "graph": graph_check,
        "graph_return_handle": _result(
            "graph_return_handle",
            passed=graph_return_ok,
            evidence={
                "returns_interned_function_ref": graph_return_ok,
                "source_check": "graph.evidence.node.returns_interned_function_ref",
            },
            error=None if graph_return_ok else "node() return does not preserve self.keys[key] in FunctionRef",
            status=None if graph_return_ok else "unresolved",
        ),
        "controls_integral": _check_controls(sources.get("controls")),
        "core_constructor": _check_core(sources.get("core")),
        "lower_collar_closure": _check_lower(sources.get("lower_collar")),
        "ordinary_z_derivative": _check_c3(sources.get("c3_derivative")),
    }
    catalog, catalog_path = _catalog(saved_result)
    occurrences = _occurrence_records(saved_result)
    records = [
        _row_record(
            row,
            source_checks=source_checks,
            occurrence=occurrences.get(row.get("integral_node_id")) if isinstance(row, Mapping) else None,
        )
        for row in catalog
    ]
    if len(records) != EXPECTED_COUNT:
        missing = max(0, EXPECTED_COUNT - len(records))
        for index in range(missing):
            records.append(
                {
                    "status": "unresolved",
                    "source_name_binding_candidate": False,
                    "candidate_basis": "ast_fragment_pattern_only_not_source_call_proof",
                    "candidate_actual_integral_binder": False,
                    "candidate_actual_integral_binder_proved": False,
                    "actual_graph_owner_identity": False,
                    "binder_identity_proved": False,
                    "parent_identity_proved": False,
                    "parent_identity_verified": False,
                    "binder_mapping_installed": False,
                    "source_callback_certified": False,
                    "source_function_equality_proved": False,
                    "source_call_binding_proved": False,
                    "source_data_flow_proved": False,
                    "ast_fragment_pattern_only": True,
                    "m15_claimed": False,
                    "parent_identity_unproved": True,
                    "saved_occurrence": _saved_occurrence_summary(None),
                    "unresolved_reasons": ["catalog has fewer than seven records; missing row " + str(index)],
                }
            )
    candidate_count = sum(item.get("source_name_binding_candidate") is True for item in records)
    if catalog_path == "missing":
        status = "unresolved"
    elif candidate_count == EXPECTED_COUNT and len(catalog) == EXPECTED_COUNT:
        status = "conditional_candidates"
    elif candidate_count:
        status = "partial_conditional_candidates"
    else:
        status = "unresolved"
    prior_error_fields: dict[str, Any] = {}
    if isinstance(saved_result, Mapping):
        for key in ("error", "traceback", "saved_error"):
            if key in saved_result:
                prior_error_fields["saved_result_" + key] = _clone(saved_result.get(key))
    result = {
        "schema": RESULT_SCHEMA,
        "status": status,
        "catalog_path": catalog_path,
        "catalog_record_count": len(catalog),
        "records": records,
        "source_checks": source_checks,
        "source_name_binding_candidate": bool(candidate_count == EXPECTED_COUNT and len(catalog) == EXPECTED_COUNT),
        "candidate_basis": "ast_fragment_pattern_only_not_source_call_proof",
        "conditional_candidate_count": candidate_count,
        "saved_occurrence_counts_preserved": bool(occurrences),
        "prior_errors_preserved": True,
        "graph_walk_performed": False,
        "graph_read": False,
        "graph_traversed": False,
        "source_graph_imported": False,
        "source_value_evaluated": False,
        "integrals_evaluated": False,
        "numeric_field_value_computed": False,
        "actual_graph_owner_identity": False,
        "candidate_actual_integral_binder": False,
        "candidate_actual_integral_binder_proved": False,
        "binder_identity_proved": False,
        "parent_identity_proved": False,
        "parent_identity_verified": False,
        "binder_mapping_installed": False,
        "source_callback_certified": False,
        "source_function_equality_proved": False,
        "source_call_binding_proved": False,
        "source_data_flow_proved": False,
        "ast_fragment_pattern_only": True,
        "m15_claimed": False,
        "scope": {
            "source_name_binding_candidate_only": True,
            "actual_graph_owner_identity": False,
            "binder_identity_proved": False,
            "parent_identity_proved": False,
            "parent_identity_verified": False,
            "binder_mapping_installed": False,
            "source_callback_certified": False,
            "source_function_equality_proved": False,
            "source_call_binding_proved": False,
            "source_data_flow_proved": False,
            "ast_fragment_pattern_only": True,
            "parent_identity_unproved_for_missing_parent_ids": True,
            "occurrence_paths_rewalked": False,
            "numeric_values": False,
            "m15": False,
        },
    }
    result.update(prior_error_fields)
    return result


def audit_source_name_binding(saved_result: Any, source_texts: Any = None) -> dict[str, Any]:
    """Descriptive alias for callers that prefer an audit-named API."""

    return analyze(saved_result, source_texts)


def prove(saved_result: Any, source_texts: Any = None) -> dict[str, Any]:
    """Compatibility alias; the result remains conditional, never a proof."""

    return analyze(saved_result, source_texts)


def check_source(saved_result: Any, source_texts: Any = None) -> dict[str, Any]:
    """Short alias for a runner that names its source consumer ``check_source``."""

    return analyze(saved_result, source_texts)


__all__ = [
    "SCHEMA",
    "RESULT_SCHEMA",
    "FIXTURE_SCHEMA",
    "EXPECTED_COUNT",
    "CORE_VARIABLE",
    "LOWER_VARIABLE",
    "analyze",
    "audit_source_name_binding",
    "prove",
    "check_source",
]
