"""Fresh synthetic checks for :mod:`n1_integral_name_binding_v1`.

The fixture passes a new seven-row catalog and short source fragments to the
supplied module.  It never imports a producer, reads a graph archive, or
reuses an earlier checker.  Cases cover a changed symbol, a deleted name
lookup, same-name candidate forgery, a missing derivative branch, and missing
required source text.
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any


SCHEMA = "n1_integral_name_binding_fixture_v1"


GRAPH = '''
import json
class FunctionRef: pass
class FunctionTransportGraph:
 def __init__(self):
  self.nodes=[]; self.keys={}
 def node(self, operation, **attributes):
  row = dict(operation=operation, **attributes)
  key = json.dumps(row, sort_keys=True, separators=(',', ':'))
  if key not in self.keys:
   self.keys[key] = len(self.nodes)
   self.nodes.append(row)
  return FunctionRef(self, self.keys[key])
 def symbol(self, name):
  return self.node('bound_variable', name=name)
'''

CONTROLS = '''
def integral(g, integrand, variable, lower, upper, **metadata):
 g.ids([integrand, lower, upper])
 return g.node('definite_integral', integrand=integrand.node,
  variable=variable, lower=lower.node, upper=upper.node,
  exact_function_integral=True, numerical_value_not_installed=True,
  **metadata)
'''

CORE = '''
class CurrentOriginalN1CoreFixedPointFunctions:
 def __init__(self):
  v = g.symbol('original_true_core_gauge_v')
  gv = g.node('function_substitution', expression=gradient.value.node,
      variable=z.node, value=v.node)
  G = controls.integral(g, gv, g.nodes[v.node]['name'], anchor, z,
      original_H_root_anchor=True)
'''

LOWER = '''
def integral(g, jet, var, lower, upper, **metadata):
 return controls.integral(g, jet, g.nodes[var.node]['name'], lower, upper,
     endpoint_independent_of_Z=True, **metadata)

def collar_fields(g):
 v = g.symbol('original_lower_collar_ratio_v')
 body = g.mul(v, v)
 return integral(g, body, v, g.zero, s, original_ratio_ODE=True)
'''

C3 = '''
import copy
class OrdinaryZDerivative:
 def rule(self, i):
  row = self.g.nodes[i]
  op = row['operation']
  if op == 'definite_integral':
   body = self.derivative(row['integrand'])
   attrs = copy.deepcopy(row)
   attrs.pop('operation')
   attrs['integrand'] = body.node
   attrs['ordinary_Z_derivative_order'] = 1
   return self.g.node('definite_integral', **attrs)
'''


def _catalog() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    specs = [
        (9100, 9101, 9000, "original_true_core_gauge_v", {"actual_anchored_gauge_integral": True}),
        (9200, 9201, 9001, "original_lower_collar_ratio_v", {"ordinary_slow_Z_derivative_order": 0}),
        (9202, 9203, 9001, "original_lower_collar_ratio_v", {"ordinary_slow_Z_derivative_order": 1}),
        (9204, 9205, 9001, "original_lower_collar_ratio_v", {"ordinary_slow_Z_derivative_order": 2}),
        (9300, 9301, 9001, "original_lower_collar_ratio_v", {
            "derivative_of_same_defining_integral": True,
            "ordinary_Z_derivative_order": 1,
        }),
        (9302, 9303, 9001, "original_lower_collar_ratio_v", {
            "derivative_of_same_defining_integral": True,
            "ordinary_Z_derivative_order": 1,
        }),
        (9304, 9305, 9001, "original_lower_collar_ratio_v", {
            "derivative_of_same_defining_integral": True,
            "ordinary_Z_derivative_order": 1,
        }),
    ]
    for integral_id, integrand, candidate, variable, extra in specs:
        original = {
            "operation": "definite_integral",
            "integrand": integrand,
            "lower": 0,
            "upper": 1,
            "variable": variable,
        }
        original.update(extra)
        rows.append(
            {
                "integral_node_id": integral_id,
                "unique_node_id": candidate,
                "variable_name": variable,
                "binding_status": "unique_candidate",
                "candidate_node_ids": [candidate],
                "original_record": original,
            }
        )
    return rows


def _occurrences(catalog: list[Mapping[str, Any]]) -> list[dict[str, Any]]:
    result = []
    for index, row in enumerate(catalog):
        result.append(
            {
                "integral_node_id": row["integral_node_id"],
                "status": "found_but_unresolved",
                "unshadowed_occurrences": [{"saved": index}],
                "opaque_frontier": [{"reason": "synthetic_frontier"} for _ in range(index + 1)],
                "errors": [{"code": "synthetic_saved_error"}],
            }
        )
    return result


def _saved() -> dict[str, Any]:
    catalog = _catalog()
    return {
        "schema": "synthetic_runner_decoded_scoped_catalog_v1",
        "authenticated_and_decoded": True,
        "integral_candidates": catalog,
        "leaf_raw_output": {"records": _occurrences(catalog)},
    }


def _sources() -> dict[str, str]:
    return {"graph": GRAPH, "controls": CONTROLS, "core": CORE, "lower_collar": LOWER, "c3": C3}


def _records(result: Any) -> list[Mapping[str, Any]]:
    if not isinstance(result, Mapping) or not isinstance(result.get("records"), list):
        return []
    return [row for row in result["records"] if isinstance(row, Mapping)]


def _by_variable(result: Any, variable: str) -> list[Mapping[str, Any]]:
    return [row for row in _records(result) if row.get("variable_name") == variable]


def _identity_parent_flags_false(rows: list[Mapping[str, Any]]) -> bool:
    return bool(rows) and all(
        row.get("candidate_actual_integral_binder") is False
        and row.get("candidate_actual_integral_binder_proved") is False
        and row.get("actual_graph_owner_identity") is False
        and row.get("binder_identity_proved") is False
        and row.get("parent_identity_proved") is False
        and row.get("parent_identity_verified") is False
        and row.get("parent_identity_unproved") is True
        and row.get("binder_mapping_installed") is False
        and row.get("source_callback_certified") is False
        and row.get("source_call_binding_proved") is False
        and row.get("source_data_flow_proved") is False
        and row.get("source_function_equality_proved") is False
        and row.get("m15_claimed") is False
        for row in rows
    )


def _run_case(module: Any, saved: Any, sources: Any) -> tuple[Any, str | None]:
    function = getattr(module, "analyze", None)
    if not callable(function):
        return None, "analyze_not_callable"
    try:
        return function(saved, sources), None
    except Exception as exc:
        return None, type(exc).__name__


def _diagnostic(name: str, result: Any, checks: Mapping[str, bool], error: str | None) -> dict[str, Any]:
    item: dict[str, Any] = {
        "case": name,
        "passed": bool(checks) and all(value is True for value in checks.values()),
        "failed_checks": [key for key, value in checks.items() if value is not True],
        "status": result.get("status") if isinstance(result, Mapping) else None,
    }
    if error is not None:
        item["exception"] = error
    return item


def run_fixture(module: Any) -> dict[str, Any]:
    """Run all source-only cases against the supplied new module."""

    analyze = getattr(module, "analyze", None)
    if not callable(analyze):
        return {
            "schema": SCHEMA,
            "passed": False,
            "checks": {"analyze_callable": False},
            "diagnostics": [{"case": "input", "passed": False}],
            "case_count": 0,
            "test_count": 1,
            "synthetic_only": True,
            "graph_walk_performed": False,
        }

    base_saved = _saved()
    base_sources = _sources()
    checks: dict[str, bool] = {"analyze_callable": True}
    diagnostics: list[dict[str, Any]] = []

    result, error = _run_case(module, base_saved, base_sources)
    rows = _records(result)
    checks["valid_result_shape"] = isinstance(result, Mapping) and len(rows) == 7
    checks["valid_all_conditional"] = bool(rows) and all(row.get("source_name_binding_candidate") is True for row in rows)
    checks["valid_identity_false"] = isinstance(result, Mapping) and result.get("binder_identity_proved") is False
    checks["valid_mapping_false"] = isinstance(result, Mapping) and result.get("binder_mapping_installed") is False
    checks["valid_callback_false"] = isinstance(result, Mapping) and result.get("source_callback_certified") is False
    checks["valid_all_identity_parent_flags_false"] = _identity_parent_flags_false(rows)
    checks["valid_source_checks_pattern_only"] = isinstance(result, Mapping) and all(
        isinstance(check, Mapping)
        and check.get("evidence_scope") == "ast_fragment_pattern_only"
        and check.get("data_flow_proved") is False
        and check.get("source_call_binding_proved") is False
        and check.get("graph_node_identity_proved") is False
        for check in result.get("source_checks", {}).values()
    )
    checks["valid_saved_counts_preserved"] = bool(rows) and all(
        isinstance(row.get("saved_occurrence"), Mapping)
        and row["saved_occurrence"].get("witness_count") == 1
        for row in rows
    )
    diagnostics.append(_diagnostic("valid", result, {key: value for key, value in checks.items() if key.startswith("valid_")}, error))

    changed_graph = dict(base_sources)
    changed_graph["graph"] = GRAPH.replace("'bound_variable'", "'changed_name'")
    result, error = _run_case(module, base_saved, changed_graph)
    rows = _records(result)
    checks["modified_symbol_rejected"] = bool(rows) and all(row.get("source_name_binding_candidate") is False for row in rows)
    diagnostics.append(_diagnostic("modified_symbol", result, {"rejected": checks["modified_symbol_rejected"]}, error))

    changed_return = dict(base_sources)
    changed_return["graph"] = GRAPH.replace(
        "return FunctionRef(self, self.keys[key])",
        "return FunctionRef(self, self.keys[0])",
    )
    result, error = _run_case(module, base_saved, changed_return)
    rows = _records(result)
    checks["changed_return_handle_rejected"] = bool(rows) and all(
        row.get("source_name_binding_candidate") is False for row in rows
    )
    diagnostics.append(_diagnostic("changed_return_handle", result, {"rejected": checks["changed_return_handle_rejected"]}, error))

    deleted_lookup = dict(base_sources)
    deleted_lookup["core"] = CORE.replace("g.nodes[v.node]['name']", "v.node")
    result, error = _run_case(module, base_saved, deleted_lookup)
    core_rows = _by_variable(result, "original_true_core_gauge_v")
    lower_rows = _by_variable(result, "original_lower_collar_ratio_v")
    checks["deleted_variable_lookup_rejected"] = bool(core_rows) and all(
        row.get("source_name_binding_candidate") is False for row in core_rows
    ) and bool(lower_rows) and any(row.get("source_name_binding_candidate") is True for row in lower_rows)
    diagnostics.append(_diagnostic("deleted_variable_lookup", result, {"rejected": checks["deleted_variable_lookup_rejected"]}, error))

    forged_saved = _saved()
    forged_saved["integral_candidates"][0]["unique_node_id"] = 99999
    forged_saved["integral_candidates"][0]["candidate_node_ids"] = [99999]
    result, error = _run_case(module, forged_saved, base_sources)
    forged_rows = _by_variable(result, "original_true_core_gauge_v")
    checks["same_name_forgery_not_identity"] = bool(forged_rows) and forged_rows[0].get("source_name_binding_candidate") is True and forged_rows[0].get("binder_identity_proved") is False
    checks["same_name_forgery_owner_false"] = isinstance(result, Mapping) and result.get("actual_graph_owner_identity") is False
    diagnostics.append(_diagnostic("same_name_forgery", result, {key: value for key, value in checks.items() if key.startswith("same_name_forgery_")}, error))

    forged_parent_saved = _saved()
    forged_parent_saved["integral_candidates"][4]["parent_node_id"] = 777777
    forged_parent_saved["integral_candidates"][4]["original_record"]["parent_node_id"] = 777777
    result, error = _run_case(module, forged_parent_saved, base_sources)
    parent_rows = [
        row for row in _records(result) if row.get("integral_node_id") == 9300
    ]
    checks["forged_parent_not_identity"] = bool(parent_rows) and parent_rows[0].get("parent_node_id") == 777777 and parent_rows[0].get("parent_identity_unproved") is True
    checks["forged_parent_reference_marked_unverified"] = bool(parent_rows) and parent_rows[0].get("parent_node_id_unverified_reference") is True
    checks["forged_parent_all_identity_flags_false"] = _identity_parent_flags_false(parent_rows)
    diagnostics.append(_diagnostic("forged_parent_id", result, {key: value for key, value in checks.items() if key.startswith("forged_parent_")}, error))

    missing_derivative = dict(base_sources)
    missing_derivative["c3"] = C3.replace("if op == 'definite_integral':", "if op == 'not_an_integral':")
    result, error = _run_case(module, base_saved, missing_derivative)
    derivative_rows = [
        row for row in _by_variable(result, "original_lower_collar_ratio_v")
        if row.get("derivative_rule_required") is True
    ]
    base_lower_rows = [
        row for row in _by_variable(result, "original_lower_collar_ratio_v")
        if row.get("derivative_rule_required") is False
    ]
    checks["missing_derivative_rejected"] = bool(derivative_rows) and all(
        row.get("source_name_binding_candidate") is False for row in derivative_rows
    ) and bool(base_lower_rows) and any(row.get("source_name_binding_candidate") is True for row in base_lower_rows)
    diagnostics.append(_diagnostic("missing_derivative_rule", result, {"rejected": checks["missing_derivative_rejected"]}, error))

    missing_required = dict(base_sources)
    del missing_required["lower_collar"]
    result, error = _run_case(module, base_saved, missing_required)
    core_rows = _by_variable(result, "original_true_core_gauge_v")
    lower_rows = _by_variable(result, "original_lower_collar_ratio_v")
    checks["missing_required_source_unresolved"] = bool(core_rows) and all(
        row.get("source_name_binding_candidate") is True for row in core_rows
    ) and bool(lower_rows) and all(
        row.get("source_name_binding_candidate") is False for row in lower_rows
    )
    diagnostics.append(_diagnostic("missing_lower_source", result, {"unresolved": checks["missing_required_source_unresolved"]}, error))

    return {
        "schema": SCHEMA,
        "passed": all(checks.values()),
        "checks": checks,
        "diagnostics": diagnostics,
        "case_count": len(diagnostics),
        "test_count": len(checks),
        "synthetic_only": True,
        "graph_walk_performed": False,
        "graph_read": False,
        "source_graph_imported": False,
        "source_value_evaluated": False,
        "integrals_evaluated": False,
        "numeric_values_computed": False,
        "binder_identity_proved": False,
        "binder_mapping_installed": False,
        "source_callback_certified": False,
    }


__all__ = ["SCHEMA", "run_fixture"]
