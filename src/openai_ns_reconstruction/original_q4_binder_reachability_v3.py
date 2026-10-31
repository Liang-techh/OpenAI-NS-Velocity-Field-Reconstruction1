"""Fail-closed lexical reachability for one saved Q4 substitution (v2).

Only the supplied exact graph table is read.  This is structural evidence: no
numeric field, source value, integral, or replacement is evaluated.  The body
and value roots are scanned independently, with substitution values remaining
in the outer scope.  A same-name nested substitution binder shadows its body
even when the binder row reuses the target node id; a same-name different-id
binder is retained as an explicit unresolved frontier.
"""

from __future__ import annotations


SCHEMA = "original_q4_binder_reachability_v3"
RESULT_SCHEMA = SCHEMA
MAX_VISITS = 3000
MAX_DEPTH = 128
MAX_RECORDS = 512

_SAFE_CHILD_FIELDS = {
    "sum": ("arguments",),
    "product": ("arguments",),
    "negative": ("argument",),
    "positive_quotient": ("numerator", "denominator"),
    "analytic_unary": ("argument",),
    "function_substitution": ("expression", "variable", "value"),
    "definite_integral": ("integrand", "variable", "lower", "upper"),
    "exact_real_comparison": ("left", "right"),
    "logical_or": ("arguments",),
}
_KNOWN_OPERATIONS = {"exact_rational", "bound_variable"} | set(_SAFE_CHILD_FIELDS)


def _empty_summary():
    return {
        "target_id_observed": False,
        "target_id_unshadowed": False,
        "same_name_different_id_conflict": False,
        "unsupported_binder": False,
        "opaque": False,
        "budget_truncated": False,
        "cycle": False,
        "format_error": False,
        "evidence_truncated": False,
        "occurrences": [],
        "shadow_records": [],
        "frontier": [],
        "errors": [],
    }


def _copy_summary(summary):
    copied = {
        key: bool(summary[key])
        for key in (
            "target_id_observed",
            "target_id_unshadowed",
            "same_name_different_id_conflict",
            "unsupported_binder",
            "opaque",
            "budget_truncated",
            "cycle",
            "format_error",
            "evidence_truncated",
        )
    }
    for key in ("occurrences", "shadow_records", "frontier", "errors"):
        copied[key] = [dict(item) for item in summary[key]]
    return copied


def _merge_summary(into, child):
    for key in (
        "target_id_observed",
        "target_id_unshadowed",
        "same_name_different_id_conflict",
        "unsupported_binder",
        "opaque",
        "budget_truncated",
        "cycle",
        "format_error",
        "evidence_truncated",
    ):
        into[key] = bool(into[key] or child[key])
    if child["evidence_truncated"]:
        into["opaque"] = True
    for key in ("occurrences", "shadow_records", "frontier", "errors"):
        for item in child[key]:
            if len(into[key]) >= MAX_RECORDS:
                into["evidence_truncated"] = True
                into["opaque"] = True
                break
            into[key].append(dict(item))


def _is_handle(value, nodes):
    return type(value) is int and 0 <= value < len(nodes) and isinstance(nodes[value], dict)


def _state_shadowed(state, target_name):
    return any(token[2] == target_name for token in state)


def _state_json(state):
    return [
        {"kind": token[0], "binder_node_id": token[1], "name": token[2]}
        for token in state
    ]


def _frontier(node_id, operation, reason):
    return {"node_id": node_id, "operation": operation, "reason": reason}


class _Scan:
    def __init__(self, nodes, target_id, target_name, root):
        self.nodes = nodes
        self.target_id = target_id
        self.target_name = target_name
        self.root = root
        self.memo = {}
        self.active = set()
        self.visits = 0
        self.memo_hits = 0
        self.unique_nodes = 0
        self.max_depth_observed = 0

    def _append(self, summary, key, item):
        if len(summary[key]) >= MAX_RECORDS:
            summary["evidence_truncated"] = True
            summary["opaque"] = True
            return
        summary[key].append(item)

    def _error(self, summary, code, node_id=None, operation=None, field=None, detail=None):
        item = {"code": code}
        if node_id is not None:
            item["node_id"] = node_id
        if operation is not None:
            item["operation"] = operation
        if field is not None:
            item["field"] = field
        if detail is not None:
            item["detail"] = detail
        summary["format_error"] = True
        self._append(summary, "errors", item)

    def _opaque(self, summary, node_id, operation, reason, **extra):
        summary["opaque"] = True
        item = _frontier(node_id, operation, reason)
        item.update(extra)
        self._append(summary, "frontier", item)

    def _child(self, summary, node_id, operation, field, value):
        if field == "arguments":
            if not isinstance(value, list):
                self._error(summary, "operand_list_shape", node_id, operation, field, "arguments must be a list")
                return None
            children = []
            for index, child in enumerate(value):
                if not _is_handle(child, self.nodes):
                    self._error(summary, "operand_handle_shape", node_id, operation, field, "invalid child at index %d" % index)
                    return None
                children.append(child)
            return children
        if not _is_handle(value, self.nodes):
            self._error(summary, "operand_handle_shape", node_id, operation, field, "operand is not a node handle")
            return None
        return value

    def _shadow(self, summary, node_id, binder_node_id, name, state, shadows_target, authenticated=True):
        self._append(
            summary,
            "shadow_records",
            {
                "node_id": node_id,
                "kind": "function_substitution",
                "binder_node_id": binder_node_id,
                "binder_name": name,
                "scope": "expression",
                "shadows_target": shadows_target,
                # A matching local handle is a syntactic observation, not a
                # proof of lexical source identity or authenticated binding.
                "identity_authenticated": False,
                "binder_same_target_id": binder_node_id == self.target_id,
                "outer_shadow_state": _state_json(state),
            },
        )

    def _visit(self, node_id, state, depth):
        summary = _empty_summary()
        if self.visits >= MAX_VISITS:
            summary["budget_truncated"] = True
            self._opaque(summary, node_id, self.nodes[node_id].get("operation"), "visit_budget")
            return summary
        self.visits += 1
        self.max_depth_observed = max(self.max_depth_observed, depth)
        if depth > MAX_DEPTH:
            summary["budget_truncated"] = True
            self._opaque(summary, node_id, self.nodes[node_id].get("operation"), "depth_budget")
            return summary

        key = (node_id, tuple(state))
        if key in self.active:
            summary["cycle"] = True
            self._opaque(summary, node_id, self.nodes[node_id].get("operation"), "cycle")
            return summary
        cached = self.memo.get(key)
        if cached is not None:
            self.memo_hits += 1
            return _copy_summary(cached)

        self.unique_nodes += 1
        row = self.nodes[node_id]
        operation = row.get("operation")
        if not isinstance(operation, str):
            self._error(summary, "operation_shape", node_id, None, "operation", "operation must be a string")
            self._opaque(summary, node_id, operation, "invalid_operation")
            return summary
        if operation not in _KNOWN_OPERATIONS:
            # This includes original_lazy_flat_branch; no child or active /
            # inactive branch metadata is inspected.
            self._opaque(summary, node_id, operation, "unsupported_operation")
            return summary

        self.active.add(key)
        try:
            if operation == "exact_rational":
                return summary

            if operation == "bound_variable":
                name = row.get("name")
                if not isinstance(name, str) or not name:
                    self._error(summary, "bound_variable_name_shape", node_id, operation, "name", "name must be a nonempty string")
                    return summary
                shadowed = _state_shadowed(state, self.target_name)
                same_name = name == self.target_name
                target_id = node_id == self.target_id
                conflict = same_name and not target_id and not shadowed
                summary["target_id_observed"] |= target_id
                summary["target_id_unshadowed"] |= target_id and not shadowed
                summary["same_name_different_id_conflict"] |= conflict
                self._append(
                    summary,
                    "occurrences",
                    {
                        "node_id": node_id,
                        "name": name,
                        "same_name_as_target": same_name,
                        "is_target_id": target_id,
                        "shadowed": shadowed,
                        "target_id_observed": target_id,
                        "target_id_unshadowed": target_id and not shadowed,
                        "same_name_different_id_conflict": conflict,
                        "shadow_state": _state_json(state),
                    },
                )
                return summary

            if operation == "function_substitution":
                expression = self._child(summary, node_id, operation, "expression", row.get("expression"))
                variable = self._child(summary, node_id, operation, "variable", row.get("variable"))
                value = self._child(summary, node_id, operation, "value", row.get("value"))
                if expression is None or variable is None or value is None:
                    return summary
                binder_row = self.nodes[variable]
                binder_operation = binder_row.get("operation")
                binder_name = binder_row.get("name")
                binder_valid = binder_operation == "bound_variable" and isinstance(binder_name, str) and bool(binder_name)
                if not binder_valid:
                    summary["unsupported_binder"] = True
                    self._opaque(
                        summary,
                        node_id,
                        operation,
                        "unsupported_binder_semantics",
                        binder_node_id=variable,
                        binder_operation=binder_operation,
                        body_node_id=expression,
                    )
                    _merge_summary(summary, self._visit(value, state, depth + 1))
                    return summary

                same_name = binder_name == self.target_name
                different_id_same_name = same_name and variable != self.target_id
                # Every same-name nested binder shadows the body, including a
                # row that reuses target_id.  A different-id declaration is
                # additionally unresolved because source identity is not
                # authenticated by this structural consumer.
                self._shadow(summary, node_id, variable, binder_name, state, same_name, not different_id_same_name)
                if different_id_same_name:
                    summary["unsupported_binder"] = True
                    self._opaque(
                        summary,
                        node_id,
                        operation,
                        "same_name_nested_binder_identity_unresolved",
                        binder_node_id=variable,
                        binder_name=binder_name,
                        body_node_id=expression,
                    )
                body_state = state + (("function_substitution", variable, binder_name),) if same_name else state
                _merge_summary(summary, self._visit(expression, body_state, depth + 1))
                # Replacement values are always in the outer state.
                _merge_summary(summary, self._visit(value, state, depth + 1))
                return summary

            if operation == "definite_integral":
                # Do not dereference the integrand.  It remains an opaque
                # frontier until a string-to-binder identity is authenticated.
                raw_integrand = row.get("integrand")
                integrand = raw_integrand if type(raw_integrand) is int and raw_integrand >= 0 else None
                if integrand is None:
                    self._error(summary, "integrand_handle_shape", node_id, operation, "integrand", "integrand must be a nonnegative integer handle")
                lower = self._child(summary, node_id, operation, "lower", row.get("lower"))
                upper = self._child(summary, node_id, operation, "upper", row.get("upper"))
                name = row.get("variable")
                if not isinstance(name, str) or not name:
                    summary["unsupported_binder"] = True
                    self._error(summary, "integral_variable_shape", node_id, operation, "variable", "variable must be a nonempty string")
                    return summary
                if integrand is None or lower is None or upper is None:
                    return summary
                summary["unsupported_binder"] = True
                self._opaque(
                    summary,
                    node_id,
                    operation,
                    "integral_string_binder_unresolved",
                    binder_name=name,
                    integrand_node_id=integrand,
                )
                self._append(
                    summary,
                    "shadow_records",
                    {
                        "node_id": node_id,
                        "kind": "definite_integral",
                        "binder_node_id": None,
                        "binder_name": name,
                        "scope": "integrand",
                        "name_matches_target": name == self.target_name,
                        "shadows_target": None,
                        "identity_authenticated": False,
                        "outer_shadow_state": _state_json(state),
                    },
                )
                _merge_summary(summary, self._visit(lower, state, depth + 1))
                _merge_summary(summary, self._visit(upper, state, depth + 1))
                return summary

            for field in _SAFE_CHILD_FIELDS[operation]:
                child = self._child(summary, node_id, operation, field, row.get(field))
                if child is None:
                    continue
                if field == "arguments":
                    for item in child:
                        _merge_summary(summary, self._visit(item, state, depth + 1))
                else:
                    _merge_summary(summary, self._visit(child, state, depth + 1))
            return summary
        finally:
            self.active.discard(key)
            if not summary["budget_truncated"] and not summary["cycle"] and not summary["format_error"]:
                self.memo[key] = _copy_summary(summary)

    def run(self, root_id):
        summary = self._visit(root_id, (), 0)
        occurrences = summary["occurrences"]
        target_id_occurrences = sum(1 for item in occurrences if item.get("is_target_id"))
        target_id_unshadowed_occurrences = sum(1 for item in occurrences if item.get("target_id_unshadowed"))
        target_id_shadowed_occurrences = sum(1 for item in occurrences if item.get("is_target_id") and item.get("shadowed"))
        same_name_different_id_occurrences = sum(1 for item in occurrences if item.get("same_name_as_target") and not item.get("is_target_id"))
        same_name_different_id_conflict_occurrences = sum(1 for item in occurrences if item.get("same_name_different_id_conflict"))
        output = {
            "root_node_id": root_id,
            "root": self.root,
            "visits": self.visits,
            "unique_nodes": self.unique_nodes,
            "memo_hits": self.memo_hits,
            "max_depth_observed": self.max_depth_observed,
            "max_visits": MAX_VISITS,
            "max_depth": MAX_DEPTH,
            "max_records": MAX_RECORDS,
            "budget_truncated": bool(summary["budget_truncated"]),
            "evidence_truncated": bool(summary["evidence_truncated"]),
            "cycle": bool(summary["cycle"]),
            "target_id_observed": bool(summary["target_id_observed"]),
            "target_id_unshadowed": bool(summary["target_id_unshadowed"]),
            "target_id_occurrences": target_id_occurrences,
            "target_id_unshadowed_occurrences": target_id_unshadowed_occurrences,
            "target_id_shadowed_occurrences": target_id_shadowed_occurrences,
            "same_name_different_id_occurrences": same_name_different_id_occurrences,
            "same_name_different_id_conflict_occurrences": same_name_different_id_conflict_occurrences,
            "same_name_different_id_conflict": bool(summary["same_name_different_id_conflict"]),
            "unsupported_binder": bool(summary["unsupported_binder"]),
            "opaque": bool(summary["opaque"]),
            "format_error": bool(summary["format_error"]),
            "occurrences": [dict(item) for item in occurrences],
            "shadow_records": [dict(item) for item in summary["shadow_records"]],
            "frontier": [dict(item) for item in summary["frontier"]],
            "errors": [dict(item) for item in summary["errors"]],
        }
        for collection in (output["occurrences"], output["shadow_records"], output["frontier"], output["errors"]):
            for item in collection:
                item["root"] = self.root
        return output


def _empty_scan():
    return {
        "root_node_id": None,
        "root": None,
        "visits": 0,
        "unique_nodes": 0,
        "memo_hits": 0,
        "max_depth_observed": 0,
        "max_visits": MAX_VISITS,
        "max_depth": MAX_DEPTH,
        "max_records": MAX_RECORDS,
        "budget_truncated": False,
        "evidence_truncated": False,
        "cycle": False,
        "target_id_observed": False,
        "target_id_unshadowed": False,
        "target_id_occurrences": 0,
        "target_id_unshadowed_occurrences": 0,
        "target_id_shadowed_occurrences": 0,
        "same_name_different_id_occurrences": 0,
        "same_name_different_id_conflict_occurrences": 0,
        "same_name_different_id_conflict": False,
        "unsupported_binder": False,
        "opaque": False,
        "format_error": True,
        "occurrences": [],
        "shadow_records": [],
        "frontier": [],
        "errors": [],
    }


def _invalid_result(error, substitution_node=None):
    empty = _empty_scan()
    empty["errors"] = [dict(error)]
    return {
        "schema": SCHEMA,
        "status": "mismatch",
        "substitution_node": substitution_node,
        "target": None,
        "body": dict(empty),
        "value": dict(empty),
        "scans": {"body": dict(empty), "value": dict(empty)},
        "occurrences": {"body": [], "value": []},
        "shadow_records": {"body": [], "value": []},
        "frontier": [],
        "errors": [dict(error)],
        "counts": {"body_visits": 0, "value_visits": 0, "frontier": 0, "shadow_records": 0, "budget_truncated": 0, "evidence_truncated": 0},
        "scope": {
            "conditional_on_source_operator_semantics": True,
            "numeric_field_value_computed": False,
            "old_errors_removed": False,
            "complete_identity_proved": False,
            "syntactic_scan_complete": False,
        },
    }


def _normalise_substitution_node(value):
    if isinstance(value, dict) and set(value) == {"node"}:
        value = value.get("node")
    return value if type(value) is int else None


def analyze(graph, substitution_node):
    """Analyse one function-substitution row without evaluating any value."""
    requested = substitution_node
    if not isinstance(graph, dict):
        return _invalid_result({"code": "graph_shape", "detail": "graph must be an object"}, requested)
    nodes = graph.get("exact_graph_nodes")
    if not isinstance(nodes, list):
        return _invalid_result({"code": "exact_graph_nodes_shape", "detail": "exact_graph_nodes must be a list"}, requested)

    target_id = _normalise_substitution_node(substitution_node)
    if target_id is None or not _is_handle(target_id, nodes):
        return _invalid_result({"code": "substitution_node_shape", "detail": "substitution_node is not a valid node handle"}, requested)
    row = nodes[target_id]
    if row.get("operation") != "function_substitution":
        return _invalid_result({"code": "substitution_operation", "node_id": target_id, "detail": "target row is not function_substitution"}, requested)

    fields = {}
    for field in ("expression", "variable", "value"):
        value = row.get(field)
        if not _is_handle(value, nodes):
            return _invalid_result({"code": "substitution_child_shape", "node_id": target_id, "field": field, "detail": "child is not a valid node handle"}, requested)
        fields[field] = value
    binder = nodes[fields["variable"]]
    if binder.get("operation") != "bound_variable" or not isinstance(binder.get("name"), str) or not binder.get("name"):
        return _invalid_result({"code": "substitution_binder_shape", "node_id": fields["variable"], "detail": "binder must be a bound_variable with a nonempty name"}, requested)

    target_name = binder["name"]
    body_scan = _Scan(nodes, fields["variable"], target_name, "body").run(fields["expression"])
    value_scan = _Scan(nodes, fields["variable"], target_name, "value").run(fields["value"])
    body_occurrences = [dict(item) for item in body_scan["occurrences"]]
    value_occurrences = [dict(item) for item in value_scan["occurrences"]]
    body_shadow = [dict(item) for item in body_scan["shadow_records"]]
    value_shadow = [dict(item) for item in value_scan["shadow_records"]]
    all_frontier = [dict(item) for item in body_scan["frontier"] + value_scan["frontier"]]
    all_errors = [dict(item) for item in body_scan["errors"] + value_scan["errors"]]
    conflict = bool(body_scan["same_name_different_id_conflict"] or value_scan["same_name_different_id_conflict"])
    budget = bool(body_scan["budget_truncated"] or value_scan["budget_truncated"])
    evidence_truncated = bool(body_scan["evidence_truncated"] or value_scan["evidence_truncated"])
    opaque = bool(body_scan["opaque"] or value_scan["opaque"])
    unsupported_binder = bool(body_scan["unsupported_binder"] or value_scan["unsupported_binder"])
    format_error = bool(body_scan["format_error"] or value_scan["format_error"])
    target_observed = bool(body_scan["target_id_observed"] or value_scan["target_id_observed"])
    target_unshadowed = bool(body_scan["target_id_unshadowed"] or value_scan["target_id_unshadowed"])
    syntactic_scan_complete = bool(
        not format_error
        and not opaque
        and not conflict
        and not budget
        and not evidence_truncated
        and not unsupported_binder
        and not body_scan["cycle"]
        and not value_scan["cycle"]
    )
    status = "mismatch" if format_error else ("complete" if syntactic_scan_complete else "unresolved")

    return {
        "schema": SCHEMA,
        "status": status,
        "substitution_node": target_id,
        "target": {
            "substitution_node_id": target_id,
            "binder_node_id": fields["variable"],
            "binder_name": target_name,
            "body_node_id": fields["expression"],
            "value_node_id": fields["value"],
        },
        "body": body_scan,
        "value": value_scan,
        "scans": {"body": body_scan, "value": value_scan},
        "occurrences": {"body": body_occurrences, "value": value_occurrences},
        "shadow_records": {"body": body_shadow, "value": value_shadow},
        "frontier": all_frontier,
        "errors": all_errors,
        "counts": {
            "body_visits": body_scan["visits"],
            "value_visits": value_scan["visits"],
            "body_unique_nodes": body_scan["unique_nodes"],
            "value_unique_nodes": value_scan["unique_nodes"],
            "body_memo_hits": body_scan["memo_hits"],
            "value_memo_hits": value_scan["memo_hits"],
            "body_occurrences": len(body_occurrences),
            "value_occurrences": len(value_occurrences),
            "body_target_id_occurrences": body_scan["target_id_occurrences"],
            "value_target_id_occurrences": value_scan["target_id_occurrences"],
            "body_target_id_unshadowed_occurrences": body_scan["target_id_unshadowed_occurrences"],
            "value_target_id_unshadowed_occurrences": value_scan["target_id_unshadowed_occurrences"],
            "body_target_id_shadowed_occurrences": body_scan["target_id_shadowed_occurrences"],
            "value_target_id_shadowed_occurrences": value_scan["target_id_shadowed_occurrences"],
            "body_same_name_different_id_occurrences": body_scan["same_name_different_id_occurrences"],
            "value_same_name_different_id_occurrences": value_scan["same_name_different_id_occurrences"],
            "body_shadow_records": len(body_shadow),
            "value_shadow_records": len(value_shadow),
            "frontier": len(all_frontier),
            "shadow_records": len(body_shadow) + len(value_shadow),
            "budget_truncated": int(budget),
            "evidence_truncated": int(evidence_truncated),
            "format_errors": len(all_errors),
            "target_id_observed": int(target_observed),
            "target_id_unshadowed": int(target_unshadowed),
            "same_name_different_id_conflict": int(conflict),
        },
        "scope": {
            "conditional_on_source_operator_semantics": True,
            "numeric_field_value_computed": False,
            "old_errors_removed": False,
            "complete_identity_proved": False,
            "syntactic_scan_complete": syntactic_scan_complete,
            "target_id_observed": target_observed,
            "target_id_unshadowed": target_unshadowed,
            "same_name_different_id_conflict": conflict,
            "opaque_frontier_present": opaque,
            "budget_truncated": budget,
            "evidence_truncated": evidence_truncated,
            "unsupported_binder_semantics": unsupported_binder,
            "source_value_evaluated": False,
            "replacement_value_evaluated": False,
        },
    }


__all__ = ["SCHEMA", "RESULT_SCHEMA", "MAX_VISITS", "MAX_DEPTH", "MAX_RECORDS", "analyze"]

