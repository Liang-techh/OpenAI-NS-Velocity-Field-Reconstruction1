"""Independent synthetic fixtures for ``original_q4_moment_lowering_v1``.

The fixture deliberately builds graph rows instead of borrowing an archived
graph or an evaluator.  ``run_fixture`` accepts the public ``build`` callable
so its checks exercise observable IR shape and fail-closed behavior rather than
duplicating the compiler's private helpers.
"""

_FAMILY = {
    "name": "original_q4",
    "actual_five_defect_family_sha256": "a" * 64,
    "implicit_source_sha256": "b" * 64,
    "datum_enclosure_sha256": "c" * 64,
}
_P1_FAMILY = {
    "name": "p1",
    "actual_five_defect_family_sha256": "d" * 64,
    "implicit_source_sha256": "e" * 64,
    "datum_enclosure_sha256": "f" * 64,
}
_OUTER = "original_n1_true_core_rho"
_T = "original_core_inlet_moment_t"


def _clone(value):
    if isinstance(value, dict):
        return {key: _clone(item) for key, item in value.items()}
    if isinstance(value, list):
        return [_clone(item) for item in value]
    return value


def _exact(value):
    return {
        "op": "exact_rational",
        "numerator": value,
        "denominator": 1,
    }


def _base_graph():
    """Build a non-degenerate graph containing every supported standard op."""

    nodes = {
        "entry": {
            "op": "function_substitution",
            "variable": "rho_ref",
            "value": "four",
            "expression": "moment",
        },
        "four": _exact(4),
        "zero": _exact(0),
        "one": _exact(1),
        "two": _exact(2),
        "three": _exact(3),
        "rho_ref": {"op": "bound_variable", "name": _OUTER},
        "t_ref": {"op": "bound_variable", "name": _T},
        "outer_binder": {"op": "bound_variable", "name": _OUTER},
        "z_ref": {"op": "bound_variable", "name": "Z"},
        "z_parameter": {
            "op": "current_original_source_parameter",
            "source_family": _FAMILY,
            "Z_independent": True,
            "name": "Z",
            "archive_value": 101,
        },
        "component": {
            "op": "current_original_unique_analytic_core_component",
            "source_family": _FAMILY,
            "rho": "rho_ref",
            "Z_variable": "z_ref",
            "quantity": "Phi",
            "radial_order": 1,
            "ordinary_slow_Z_derivative_order": 0,
            "provider_phi_ref": "provider-phi",
            "provider_radial_ref": "provider-radial",
            "provider_slow_ref": "provider-slow",
            "provider_source_ref": "provider-source",
        },
        "phi_at_sample": {
            "op": "function_substitution",
            "variable": "rho_ref",
            "value": "rho_times_t",
            "expression": "component",
        },
        "rho_times_t": {
            "op": "product",
            "arguments": ["rho_ref", "t_ref"],
        },
        "negative": {
            "op": "negative",
            "argument": "three",
        },
        "quotient": {
            "op": "positive_quotient",
            "numerator": "three",
            "denominator": "two",
        },
        "unary": {
            "op": "analytic_unary",
            "function": "exp",
            "argument": "t_ref",
        },
        "shadow_product": {
            "op": "product",
            "arguments": ["rho_ref", "t_ref"],
        },
        "inner_shadow": {
            "op": "function_substitution",
            "variable": "rho_ref",
            "value": "two",
            "expression": "shadow_product",
        },
        "nested_integral": {
            "op": "definite_integral",
            "lower": "zero",
            "upper": "one",
            "variable": _T,
            "integrand": "t_ref",
        },
        "integrand_sum": {
            "op": "sum",
            "arguments": [
                "phi_at_sample",
                "z_parameter",
                "negative",
                "quotient",
                "unary",
                "inner_shadow",
                "nested_integral",
            ],
        },
        "moment": {
            "op": "definite_integral",
            "lower": "zero",
            "upper": "one",
            "variable": _T,
            "integrand": "integrand_sum",
        },
    }
    handles = {handle: index + 150 for index, handle in enumerate(nodes)}
    exact_graph_nodes = [{"operation": "exact_rational", "numerator": 11, "denominator": 1} for _ in range(150)]
    for row in nodes.values():
        item = _clone(row)
        item["operation"] = item.pop("op")
        exact_graph_nodes.append(item)
    exact_graph_nodes = [_rewrite_refs(item, handles) for item in exact_graph_nodes]
    return {
        "exact_graph_nodes": exact_graph_nodes,
        "_fixture_handles": handles,
        "source_family": _FAMILY,
        "bounds": {"lower": 0, "upper": 1, "unit": True},
        "old_errors": ["archived-error-preserved"],
        "numeric_field": {"value": "not-evaluated"},
    }


def fixture_graph():
    """Return a fresh copy of the independent happy-path synthetic graph."""

    return _clone(_base_graph())


def _row(graph, node_id):
    handles = graph.get("_fixture_handles", {})
    index = handles.get(node_id)
    if isinstance(index, int) and index < len(graph.get("exact_graph_nodes", [])):
        return graph["exact_graph_nodes"][index]
    return None


def _rewrite_refs(value, handles):
    if isinstance(value, str) and value in handles:
        return handles[value]
    if isinstance(value, dict):
        # Operator names and other metadata are literals, never node handles.
        protected = {"operation", "name", "function", "quantity", "source_family"}
        return {key: _clone(item) if key in protected or key.startswith("provider")
                else _rewrite_refs(item, handles) for key, item in value.items()}
    if isinstance(value, list):
        return [_rewrite_refs(item, handles) for item in value]
    return value


def _reindex_fixture(graph):
    handles = graph.get("_fixture_handles", {})
    graph["exact_graph_nodes"] = [_rewrite_refs(row, handles) for row in graph.get("exact_graph_nodes", [])]
    return graph


def _append_row(graph, node_id, operation, **fields):
    row = {"operation": operation}
    row.update(fields)
    graph.setdefault("_fixture_handles", {})[node_id] = len(graph.setdefault("exact_graph_nodes", []))
    graph.setdefault("exact_graph_nodes", []).append(row)


def _walk(value, found):
    if isinstance(value, dict):
        found.append(value)
        for child in value.values():
            _walk(child, found)
    elif isinstance(value, list):
        for child in value:
            _walk(child, found)
    return found


def _nodes_with(result, key, expected=None):
    values = []
    for item in _walk(result, []):
        if key in item and (expected is None or item.get(key) == expected):
            values.append(item)
    return values


def _case(build, name, graph, inlet=150, family=_FAMILY, check=None):
    try:
        result = build(graph, inlet, family)
    except Exception as exc:
        return {
            "status": "failed",
            "name": name,
            "error": {"kind": "builder_exception", "text": str(exc)},
        }
    try:
        passed = bool(check(result)) if check is not None else True
    except Exception as exc:
        return {
            "status": "failed",
            "name": name,
            "error": {"kind": "check_exception", "text": str(exc)},
        }
    return {
        "status": "passed" if passed else "failed",
        "name": name,
        "builder_status": result.get("status") if isinstance(result, dict) else None,
        "observed": {
            "scope": result.get("scope") if isinstance(result, dict) else None,
            "errors": result.get("errors") if isinstance(result, dict) else None,
        },
    }


def _is_four(node):
    return (
        isinstance(node, dict)
        and node.get("op") == "exact_rational"
        and node.get("numerator") == 4
        and node.get("denominator") == 1
    )


def _is_fresh(node):
    return isinstance(node, dict) and node.get("kind") == "fresh_binder" and isinstance(node.get("token"), str)


def _happy_check(result):
    if not isinstance(result, dict) or result.get("status") != "lowered":
        return False
    if result.get("source_family") != _FAMILY:
        return False
    scope = result.get("scope", {})
    if scope.get("capture_safe_for_this_IR") is not True:
        return False
    if scope.get("graph_runtime_semantics_proved") is not False:
        return False
    if scope.get("numeric_field_value_computed") is not False:
        return False
    report = result.get("report", {}).get("inner_replacement", {})
    if report.get("is_exact_rational_4_times_integral_fresh_binder") is not True:
        return False
    atoms = result.get("source_atoms", [])
    atom_ops = {item.get("op") for item in atoms if isinstance(item, dict)}
    if "current_original_source_parameter" not in atom_ops:
        return False
    if "current_original_unique_analytic_core_component" not in atom_ops:
        return False
    # Inspect the lowered component's arguments independently of source row
    # construction.  This checks Phi(4*t,Z) structurally and not numerically.
    components = _nodes_with(result.get("ir"), "kind", "analytic_core_component")
    if not components:
        return False
    component = components[0]
    rho_arg = component.get("arguments", {}).get("rho")
    if not isinstance(rho_arg, dict) or rho_arg.get("op") != "product":
        return False
    factors = rho_arg.get("factors", [])
    if len(factors) != 2 or not any(_is_four(item) for item in factors) or not any(_is_fresh(item) for item in factors):
        return False
    z_arg = component.get("arguments", {}).get("Z_variable")
    if not isinstance(z_arg, dict) or z_arg.get("op") != "bound_variable" or z_arg.get("name") != "Z":
        return False
    if component.get("provider_refs", {}).get("provider_phi_ref") != "provider-phi":
        return False
    if component.get("radial_order") != 1 or component.get("ordinary_slow_Z_derivative_order") != 0:
        return False
    component_reports = result.get("report", {}).get("core_components", [])
    if not component_reports or component_reports[0].get("rho_is_exact_rational_4_times_fresh_binder") is not True:
        return False
    return True


def _shadow_check(result):
    if not isinstance(result, dict) or result.get("status") != "lowered":
        return False
    # The inner substitution must shadow outer rho = 4 with exact 2.
    for item in _nodes_with(result.get("ir"), "substitution_variable", _OUTER):
        body = item.get("factors")
        if not isinstance(body, list) or len(body) != 2:
            continue
        if any(_is_four(factor) for factor in body):
            continue
        if any(
            isinstance(factor, dict)
            and factor.get("op") == "exact_rational"
            and factor.get("numerator") == 2
            for factor in body
        ) and any(_is_fresh(factor) for factor in body):
            return True
    return False


def _fresh_binder_check(result):
    if not isinstance(result, dict) or result.get("status") != "lowered":
        return False
    binders = [item.get("binder") for item in _nodes_with(result.get("ir"), "op", "definite_integral")]
    tokens = [item.get("token") for item in binders if _is_fresh(item)]
    return len(tokens) >= 2 and len(set(tokens)) == len(tokens) and all(token != _T for token in tokens)


def _unresolved(result):
    return isinstance(result, dict) and result.get("status") == "unresolved" and result.get("scope", {}).get("capture_safe_for_this_IR") is False


def _free_variable_graph():
    graph = _base_graph()
    _row(graph, "entry")["expression"] = "free_wrapper"
    _append_row(
        graph,
        "free_wrapper",
        "function_substitution",
        variable="free_x_binder",
        value="free_t",
        expression="free_product",
    )
    _append_row(graph, "free_x_binder", "bound_variable", name="free_x")
    _append_row(graph, "free_t", "bound_variable", name=_T)
    _append_row(
        graph,
        "free_product",
        "product",
        arguments=["free_x_binder", "moment"],
    )
    return _reindex_fixture(graph)


def _unknown_graph():
    graph = _base_graph()
    _row(graph, "component")["operation"] = "future_unknown_operator"
    return _reindex_fixture(graph)


def _cycle_graph():
    graph = _base_graph()
    _append_row(graph, "cycle", "product", arguments=["cycle"])
    _row(graph, "integrand_sum")["arguments"].append("cycle")
    return _reindex_fixture(graph)


def _wrong_endpoint_graph():
    graph = _base_graph()
    _row(graph, "moment")["upper"] = "two"
    return _reindex_fixture(graph)


def _unknown_unary_graph():
    graph = _base_graph()
    _row(graph, "unary")["function"] = "unissued_operation"
    return _reindex_fixture(graph)


def _missing_component_family_graph():
    graph = _base_graph()
    del _row(graph, "component")["source_family"]
    return _reindex_fixture(graph)


def _invalid_binder_graph():
    graph = _base_graph()
    graph["exact_graph_nodes"][graph["_fixture_handles"]["phi_at_sample"]]["variable"] = graph["_fixture_handles"]["component"]
    return graph


def _same_name_handle_conflict_graph():
    graph = _base_graph()
    rho_times_t = _row(graph, "rho_times_t")
    rho_times_t["arguments"][0] = graph["_fixture_handles"]["outer_binder"]
    return graph


def run_fixture(build):
    """Run independent happy-path and negative checks against ``build``.

    The returned object is JSON-ready and reports only IR observations.  It
    intentionally does not claim evaluator semantics, global handle identity,
    numeric computation, or removal of archived errors.
    """

    cases = {}
    cases["happy_path"] = _case(build, "happy_path", fixture_graph(), check=_happy_check)
    cases["alpha_shadow_replacement"] = _case(
        build, "alpha_shadow_replacement", fixture_graph(), check=_shadow_check
    )
    cases["same_name_binders_and_nested_shadow"] = _case(
        build, "same_name_binders_and_nested_shadow", fixture_graph(), check=_fresh_binder_check
    )
    cases["free_variable_under_replacement"] = _case(
        build,
        "free_variable_under_replacement",
        _free_variable_graph(),
        check=lambda result: (
            isinstance(result, dict)
            and result.get("status") == "lowered"
            and bool(_nodes_with(result.get("ir"), "name", _T))
            and not any(
                item.get("name") == _T and item.get("kind") == "fresh_binder"
                for item in _nodes_with(result.get("ir"), "name", _T)
            )
        ),
    )
    cases["component_provider_provenance"] = _case(
        build, "component_provider_provenance", fixture_graph(), check=_happy_check
    )
    cases["unknown_op"] = _case(build, "unknown_op", _unknown_graph(), check=_unresolved)
    cases["malformed_cycle"] = _case(build, "malformed_cycle", _cycle_graph(), check=_unresolved)
    cases["family_mismatch"] = _case(
        build, "family_mismatch", fixture_graph(), family=_P1_FAMILY, check=_unresolved
    )
    same_name_bad_hash = _clone(_FAMILY)
    same_name_bad_hash['datum_enclosure_sha256'] = 'f' * 64
    cases["same_name_different_family_hash"] = _case(
        build, "same_name_different_family_hash", fixture_graph(), family=same_name_bad_hash, check=_unresolved
    )
    cases["unknown_unary"] = _case(build, "unknown_unary", _unknown_unary_graph(), check=_unresolved)
    cases["missing_component_family"] = _case(
        build, "missing_component_family", _missing_component_family_graph(), check=_unresolved
    )
    cases["wrong_endpoints"] = _case(build, "wrong_endpoints", _wrong_endpoint_graph(), check=_unresolved)
    cases["invalid_substitution_binder_row"] = _case(
        build, "invalid_substitution_binder_row", _invalid_binder_graph(), check=_unresolved
    )
    cases["same_name_different_handle"] = _case(
        build, "same_name_different_handle", _same_name_handle_conflict_graph(), check=_unresolved
    )

    passed = all(item.get("status") == "passed" for item in cases.values())
    return {
        "status": "passed" if passed else "failed",
        "passed": passed,
        "source_family": _FAMILY,
        "cases": cases,
        "checks": {name: item.get("status") == "passed" for name, item in cases.items()},
        "scope": {
            "capture_safe_for_this_IR": passed,
            "graph_runtime_semantics_proved": False,
            "global_string_to_node_identity_proved": False,
            "numeric_field_value_computed": False,
            "old_errors_removed": False,
            "terminal_closure_proved": False,
        },
    }


__all__ = ["fixture_graph", "run_fixture"]
