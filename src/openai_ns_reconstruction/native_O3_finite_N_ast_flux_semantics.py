"""Independent source-semantic finite-``N`` O3 flux kernel.

The packet consumed by this module is an authenticated *data* packet.  The
finite-N producer is never imported or called.  Its small source function is
parsed with :mod:`ast` into a typed term environment, and the resulting term
definitions are joined to the selected native graph handles.  The returned
ordinary-Z rows are formal Taylor-coefficient expressions; they do not carry
source values.

This module intentionally keeps the source leaf separate from the previous
families.  In particular, a native graph identifier is accepted only when it
is present in the supplied literal packet.  Internal typed terms have no
native identifiers.
"""
from __future__ import annotations

import ast
import hashlib
from dataclasses import dataclass, field
from typing import Any, Callable, Mapping


COMPONENT = "native_O3_finite_N_ast_flux_semantics"
INPUT_SCHEMA = "native_O3_finite_N_ast_flux_semantics_input_v1"
RESULT_SCHEMA = "native_O3_finite_N_ast_flux_semantics_result_v1"
RAW_SCHEMA = "native_O3_density_flux_raw_literals_v1"
ORDERS = (0, 1, 2, 3, 4)
FIELDS = ("m", "h", "k", "e", "p")
BRANCHES = ("background", "modified")
FACT = (1, 1, 2, 6, 24)

# These are wrapper pins, not executable source paths.  They are checked
# against the packet's declared bytes and are never used to open a file.
WRAPPER_PINS = {
    "saved_history_result": "0521dc68caec6e046c0de1b252c518869d20250c824c8e461de1a149d3955eab",
    "source_literals": "ae3a6afbdc8de9c557dff18b3552a6b1c6c1d902decca5298b0b6adbcc216b59",
    "pointwise_lemma": "9612240a921b20b006589a7def78cc7900e17c2b6507e7cddb233bfcedba0639",
}


class _ContractError(ValueError):
    """A malformed data contract, kept private so callers receive a clear error."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise _ContractError(message)


def _clone(value: Any) -> Any:
    """Small JSON-shaped clone; no external object or producer is touched."""
    if isinstance(value, dict):
        return {key: _clone(item) for key, item in value.items()}
    if isinstance(value, list):
        return [_clone(item) for item in value]
    if isinstance(value, tuple):
        return tuple(_clone(item) for item in value)
    return value


@dataclass
class _Term:
    """Typed source term used by the restricted AST compiler."""

    kind: str
    op: str = ""
    children: tuple["_Term", ...] = ()
    data: dict[str, Any] = field(default_factory=dict)


def _term(kind: str, op: str = "", children: tuple[_Term, ...] = (), **data: Any) -> _Term:
    return _Term(kind=kind, op=op, children=children, data=data)


def _term_dict(value: _Term) -> dict[str, Any]:
    out: dict[str, Any] = {"kind": value.kind}
    if value.op:
        out["operation"] = value.op
    if value.children:
        out["children"] = [_term_dict(child) for child in value.children]
    if value.data:
        data: dict[str, Any] = {}
        for key, item in value.data.items():
            data[key] = _term_dict(item) if isinstance(item, _Term) else _clone(item)
        out["data"] = data
    return out


def _same_term(left: _Term, right: _Term) -> bool:
    if left.kind != right.kind or left.op != right.op or len(left.children) != len(right.children):
        return False
    if left.data.keys() != right.data.keys():
        return False
    for key in left.data:
        lv, rv = left.data[key], right.data[key]
        if isinstance(lv, _Term) and isinstance(rv, _Term):
            if not _same_term(lv, rv):
                return False
        elif lv != rv:
            return False
    return all(_same_term(a, b) for a, b in zip(left.children, right.children))


def _literal(node: ast.AST) -> Any:
    if isinstance(node, ast.Constant) and isinstance(node.value, (str, int, float, bool, type(None))):
        return node.value
    raise _ContractError("AST literal is not a restricted scalar")


def _subscript_key(node: ast.AST) -> str:
    part = node.slice
    if isinstance(part, ast.Index):  # pragma: no cover - Python 3.8 compatibility
        part = part.value
    value = _literal(part)
    _require(isinstance(value, str), "source dictionary key is not a string literal")
    return value


class _SourceCompiler:
    """Compile only the finite-N function's Assign/Return def-use fragment."""

    _JET_METHODS = {"add", "mul", "sub", "scale", "exp"}
    _GRAPH_METHODS = {"one", "quotient", "unary"}

    def __init__(self, source_text: str) -> None:
        _require(isinstance(source_text, str) and source_text.strip(), "finite_N source text is missing")
        try:
            self.tree = ast.parse(source_text, mode="exec")
        except SyntaxError as exc:  # make malformed source a contract error
            raise _ContractError(f"finite_N source AST failed: {exc.msg}") from exc
        self.env: dict[str, _Term] = {}
        self.operations: list[str] = []
        self.return_fields: dict[str, _Term] = {}
        functions = [node for node in self.tree.body if isinstance(node, ast.FunctionDef)]
        _require(len(functions) == 1 and functions[0].name == "finite_N_functions",
                 "finite_N source must contain exactly finite_N_functions")
        self.function = functions[0]
        _require([arg.arg for arg in self.function.args.args] == ["g", "roots", "loop", "N"],
                 "finite_N source argument contract changed")
        self.env.update({
            "g": _term("namespace", data_name="graph"),
            "roots": _term("namespace", data_name="roots"),
            "loop": _term("namespace", data_name="loop"),
            "N": _term("source_parameter", data_name="N", materialized=False),
        })

    def compile(self) -> dict[str, Any]:
        for statement in self.function.body:
            if isinstance(statement, ast.Expr) and isinstance(statement.value, ast.Constant) and isinstance(statement.value.value, str):
                continue
            if isinstance(statement, ast.Assign):
                self._assign(statement)
            elif isinstance(statement, ast.Return):
                self._return(statement)
            else:
                raise _ContractError(f"unsupported finite_N AST statement: {type(statement).__name__}")
        _require(self.return_fields, "finite_N function has no return mapping")
        for key in ("E_N", "V_N", "delta_E", "delta_V", "A_over_N", "five_signed_density_increments"):
            _require(key in self.return_fields, f"finite_N return field {key} is missing")
        increments = self.return_fields["five_signed_density_increments"]
        _require(increments.kind == "mapping", "five density increments are not a typed mapping")
        for key in FIELDS:
            _require(key in increments.data, f"finite_N increment {key} is missing")
        defs = {name: _term_dict(term) for name, term in self.env.items() if name not in {"g", "roots", "loop", "N"}}
        return {
            "function_name": self.function.name,
            "typed_definitions": defs,
            "return_fields": {key: _term_dict(value) for key, value in self.return_fields.items()},
            "validated_operations": list(self.operations),
            "source_ast_sha256": hashlib.sha256(ast.dump(self.tree, include_attributes=False).encode()).hexdigest(),
            "ordinary_jet_order": 5,
        }

    def term(self) -> dict[str, _Term]:
        """Return the already compiled typed fields without a second walk."""
        _require(self.return_fields, "typed finite_N fields requested before compile")
        return self.return_fields

    def _assign(self, statement: ast.Assign) -> None:
        _require(len(statement.targets) == 1, "chained finite_N assignment is not admitted")
        value = self._expr(statement.value)
        target = statement.targets[0]
        if isinstance(target, ast.Name):
            self.env[target.id] = value
            return
        if isinstance(target, (ast.Tuple, ast.List)):
            _require(value.kind == "tuple" and len(value.children) == len(target.elts),
                     "finite_N tuple assignment shape changed")
            for child, element in zip(value.children, target.elts):
                _require(isinstance(element, ast.Name), "tuple assignment target is not a name")
                self.env[element.id] = child
            return
        raise _ContractError("finite_N assignment target is not a name or tuple")

    def _return(self, statement: ast.Return) -> None:
        value = self._expr(statement.value)
        _require(value.kind == "mapping", "finite_N return is not a typed mapping")
        self.return_fields = dict(value.data)

    def _expr(self, node: ast.AST) -> _Term:
        if isinstance(node, ast.Name):
            _require(node.id in self.env, f"unbound finite_N name {node.id}")
            return self.env[node.id]
        if isinstance(node, ast.Constant):
            return _term("literal", value=_literal(node))
        if isinstance(node, (ast.Tuple, ast.List)):
            return _term("tuple", children=tuple(self._expr(item) for item in node.elts))
        if isinstance(node, ast.Subscript):
            base = self._expr(node.value)
            if base.kind == "rows":
                sl = node.slice
                _require(isinstance(sl, ast.Slice) and _literal(sl.lower) == 1 and sl.upper is None,
                         "only rows(jet)[1:] splice is admitted")
                return _term("rows_tail", children=base.children)
            key = _subscript_key(node)
            if base.kind == "namespace":
                return _term("source_ref", namespace=base.data.get("data_name"), key=key)
            raise _ContractError("subscript base is not a typed source namespace")
        if isinstance(node, ast.Attribute):
            base = self._expr(node.value)
            if base.kind == "namespace" and base.data.get("data_name") == "graph" and node.attr == "one":
                return _term("graph_method", op="one")
            _require(node.attr == "value" and base.kind in {"jet", "jet_blend"},
                     "only FunctionJet.value is admitted")
            return _term("jet_value", children=(base,))
        if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Add):
            left, right = self._expr(node.left), self._expr(node.right)
            if left.kind == "tuple" and right.kind == "rows_tail":
                _require(len(left.children) == 1, "FunctionJet value splice must have one value term")
                return _term("jet_blend", children=(left.children[0], right), value_splice=True)
            raise _ContractError("unrecognized finite_N addition expression")
        if isinstance(node, ast.Call):
            return self._call(node)
        raise _ContractError(f"unsupported finite_N AST expression: {type(node).__name__}")

    def _call(self, node: ast.Call) -> _Term:
        _require(not node.args or all(not isinstance(arg, ast.Starred) for arg in node.args),
                 "starred finite_N arguments are not admitted")
        if isinstance(node.func, ast.Name):
            name = node.func.id
            if name == "LoopAlgebra":
                _require(len(node.args) == 2 and not node.keywords, "LoopAlgebra constructor shape changed")
                graph, order = self._expr(node.args[0]), self._expr(node.args[1])
                _require(graph.kind == "namespace" and graph.data.get("data_name") == "graph" and
                         order.kind == "literal" and order.data.get("value") == 5,
                         "LoopAlgebra must be a typed order-five namespace")
                self.operations.append("typed_namespace")
                return _term("namespace", data_name="jet", order=5)
            if name == "FunctionJet":
                _require(len(node.args) == 1 and not node.keywords, "FunctionJet constructor shape changed")
                value = self._expr(node.args[0])
                _require(value.kind == "jet_blend", "FunctionJet must splice exp higher rows")
                return value
            if name == "rows":
                _require(len(node.args) == 1 and not node.keywords, "rows helper shape changed")
                return _term("rows", children=(self._expr(node.args[0]),))
            if name == "dict":
                _require(not node.args, "density mapping must use keyword fields")
                data: dict[str, _Term] = {}
                for keyword in node.keywords:
                    _require(keyword.arg is not None, "mapping unpack is not admitted")
                    data[keyword.arg] = self._expr(keyword.value)
                return _term("mapping", **data)
            raise _ContractError(f"unknown finite_N call {name}")
        _require(isinstance(node.func, ast.Attribute), "finite_N call target is not typed")
        base = self._expr(node.func.value)
        method = node.func.attr
        if base.kind != "namespace":
            raise _ContractError(f"method {method} is not called on a typed namespace")
        if base.data.get("data_name") == "jet":
            _require(method in self._JET_METHODS, f"unknown typed jet operation {method}")
            args = tuple(self._expr(arg) for arg in node.args)
            _require(not node.keywords, f"jet operation {method} has keyword arguments")
            if method in {"add", "mul"}:
                _require(len(args) >= 1, f"jet operation {method} is empty")
            elif method == "exp":
                _require(len(args) == 1, "jet operation exp arity changed")
            else:
                _require(len(args) == 2, f"jet operation {method} arity changed")
            self.operations.append(method)
            return _term("jet", op=method, children=args)
        _require(base.data.get("data_name") == "graph" and method in self._GRAPH_METHODS,
                 f"unknown typed graph operation {method}")
        args = tuple(self._expr(arg) for arg in node.args)
        if method == "one":
            _require(len(args) == 0 and not node.keywords, "g.one has arguments")
            self.operations.append("one")
            return _term("fixed_one")
        if method == "quotient":
            _require(len(args) >= 2 and not node.keywords, "g.quotient shape changed")
            self.operations.append("quotient")
            return _term("quotient", op="quotient", children=args)
        _require(method == "unary" and len(args) == 2 and not node.keywords,
                 "g.unary shape changed")
        _require(args[0].kind == "literal" and args[0].data.get("value") in {"exp", "expm1"},
                 "analytic unary name is not exp or expm1")
        self.operations.append(str(args[0].data["value"]))
        return _term("analytic_unary", op=str(args[0].data["value"]), children=(args[1],))


def _native_ref(node: int, source_path: str, order: int | None = None, **extra: Any) -> dict[str, Any]:
    _require(type(node) is int and node >= 0, "native graph node is not an integer")
    out: dict[str, Any] = {"kind": "native_graph_ref", "node": node, "source_path": source_path}
    if order is not None:
        out["ordinary_Z_order"] = order
    out.update(extra)
    return out


def _row_refs(literal: Mapping[str, Any], group: str, name: str) -> list[dict[str, Any]]:
    groups = literal.get("literal_rows")
    _require(isinstance(groups, Mapping) and isinstance(groups.get(group), Mapping),
             f"literal row group {group} is missing")
    entry = groups[group].get(name)
    _require(isinstance(entry, Mapping), f"literal row {group}.{name} is missing")
    ids = entry.get("selected_node_ids")
    orders = entry.get("selected_orders")
    _require(isinstance(ids, list) and isinstance(orders, list) and len(ids) == len(ORDERS) and orders == list(ORDERS),
             f"literal row {group}.{name} is not a genuine C4 selection")
    rows = entry.get("rows")
    _require(isinstance(rows, list) and len(rows) >= len(ORDERS), f"literal row {group}.{name} rows are missing")
    for index in ORDERS:
        _require(rows[index].get("ordinary_Z_order") == index and rows[index].get("node") == ids[index],
                 f"literal row {group}.{name} order metadata changed")
    path = str(entry.get("archive_json_path", f"{group}.{name}"))
    return [_native_ref(int(node), f"{path}[{order}]", order) for order, node in enumerate(ids)]


def _head(graph: Mapping[str, Any], node: int) -> Mapping[str, Any]:
    heads = graph.get("heads")
    _require(isinstance(heads, Mapping), "source graph heads are missing")
    value = heads.get(str(node), heads.get(node))
    _require(isinstance(value, Mapping), f"selected source head {node} is absent")
    return value


def _wrapper_map(spec: Mapping[str, Any], payload: Mapping[str, Any]) -> dict[str, Mapping[str, Any]]:
    wrappers = spec.get("source_wrappers", payload.get("source_wrappers"))
    if not isinstance(wrappers, Mapping):
        direct = {name: spec.get(name) for name in ("saved_history_result", "source_literals", "pointwise_lemma")}
        if all(isinstance(value, Mapping) for value in direct.values()):
            wrappers = direct
    _require(isinstance(wrappers, Mapping), "exactly three source wrappers are required")
    aliases = {
        "saved_history_result": ("saved_history_result",),
        "source_literals": ("source_literals", "density_flux_literals"),
        "pointwise_lemma": ("pointwise_lemma", "cached_density_scope"),
    }
    result: dict[str, Mapping[str, Any]] = {}
    for canonical, choices in aliases.items():
        found = next((wrappers[name] for name in choices if name in wrappers), None)
        _require(isinstance(found, Mapping), f"wrapper {canonical} is missing")
        _require(set(found.keys()) == {"sha256", "literal"}, f"wrapper {canonical} shape changed")
        _require(found.get("sha256") == WRAPPER_PINS[canonical], f"wrapper {canonical} digest changed")
        result[canonical] = found
    _require(len(wrappers) == 3, "source wrapper count changed")
    return result


def _validate_packets(wrappers: Mapping[str, Mapping[str, Any]]) -> tuple[dict[str, Any], dict[str, Any], dict[str, Any], dict[str, Any]]:
    saved = wrappers["saved_history_result"]["literal"]
    source = wrappers["source_literals"]["literal"]
    lemma = wrappers["pointwise_lemma"]["literal"]
    _require(isinstance(saved, Mapping) and isinstance(source, Mapping) and isinstance(lemma, Mapping),
             "source wrapper literals are not mappings")
    _require(source.get("schema") == RAW_SCHEMA and source.get("packet_revision") == 2,
             "raw source packet schema or revision changed")
    _require(source.get("status") == "static_literal_extraction_only", "raw source packet status changed")
    _require(source.get("scope", {}).get("chart") == "O3_transition" and
             source.get("scope", {}).get("ordinary_Z_orders") == list(ORDERS),
             "raw source chart or ordinary-Z scope changed")
    _require(lemma.get("schema") == "native_O3_generic_density_certificate_scope_v1",
             "pointwise lemma schema changed")
    abstract = lemma.get("abstract_density_schema")
    _require(isinstance(abstract, Mapping), "pointwise polynomial schema is missing")
    _require(isinstance(abstract.get("base_density_functions"), Mapping) and
             set(abstract["base_density_functions"]) == set(FIELDS),
             "pointwise base density schema changed")
    _require(isinstance(saved, Mapping), "saved history result is missing")
    saved_result = saved.get("result") if isinstance(saved.get("result"), Mapping) else saved
    _require(isinstance(saved_result, Mapping) and saved_result.get("passed") is True and
             saved_result.get("source_only") is True, "saved history result is not accepted source data")
    parent_payload = saved_result.get("payload")
    _require(isinstance(parent_payload, Mapping), "saved history payload is missing")
    return dict(source), dict(lemma), dict(saved_result), dict(parent_payload)


def _rate_map(parent_payload: Mapping[str, Any]) -> dict[str, str]:
    rows = parent_payload.get("duhamel_source_rows")
    _require(isinstance(rows, list), "saved Duhamel rows are missing")
    rates: dict[str, str] = {}
    for row in rows:
        if row.get("branch") == "background" and row.get("ordinary_Z_order") == 0 and row.get("history_key") in FIELDS:
            rates[str(row["history_key"])] = str(row.get("rate_lambda"))
    _require(set(rates) == set(FIELDS), "saved Duhamel rate map is incomplete")
    return rates


def _selector_error(branch: Any, field_name: Any, order: Any, radial_order: Any, phase_mode: Any) -> str:
    if branch not in BRANCHES:
        return "branch is not one of background or modified"
    if field_name not in FIELDS:
        return "field is not one of the five density fields"
    if isinstance(order, bool) or not isinstance(order, int) or order not in ORDERS:
        return "ordinary_Z_order must be integer Z0 through Z4"
    if radial_order != 1:
        return "only first radial order is admitted"
    if phase_mode != "ordinary_Z":
        return "held phase is not admitted"
    return "radius selector contract rejected"


def select_radius_binding(branch: str, field: str, ordinary_Z_order: int,
                          radial_order: int = 1, phase_mode: str = "ordinary_Z") -> dict[str, Any]:
    """Select a typed first-radial, ordinary-Z source binding.

    This public selector has no graph access and therefore returns a typed
    selector token.  ``derive`` enriches it with authenticated native refs.
    """
    error = _selector_error(branch, field, ordinary_Z_order, radial_order, phase_mode)
    _require(error == "radius selector contract rejected", error)
    return {
        "kind": "typed_radius_binding_selector",
        "branch": branch,
        "field": field,
        "ordinary_Z_order": ordinary_Z_order,
        "radial_order": radial_order,
        "phase_mode": phase_mode,
        "operator": "D_R",
        "numeric_evaluation": False,
    }


def _metadata_counts_from_journal(journal: list[dict[str, Any]]) -> dict[str, int]:
    """Recover only the durable callback prefix after a mid-family failure."""
    names = (
        "source_family_calls", "new_source_identity_family_calls",
        "new_ast_semantics_family_calls", "new_algebra_family_calls",
        "new_mathematical_proofs", "ast_binding_rows", "functional_alias_rows",
        "analytic_alias_rows", "density_function_identity_rows", "density_flux_rows",
        "jet_lift_rows", "derivative_rewrite_rows", "my_compatibility_rows",
        "saved_m_y_alias_rows", "final_binding_rows", "physical_binding_rows",
        "scope_binding_rows", "rewritten_physical_bindings", "retained_physical_bindings",
        "actual_negative_api_refusals", "old_math_calls", "old_provider_calls",
        "numeric_evaluation_calls", "N_materialization_calls", "source_value_reducer_calls",
        "cached_algebra_replay_calls", "failed_family_calls", "old_jet_recurrence_calls",
    )
    counts = {name: 0 for name in names}
    for record in journal:
        event, row = record.get("event"), record.get("row", {})
        if event == "family_entered":
            counts["source_family_calls"] = 1
            counts["new_source_identity_family_calls"] = 1
            counts["new_ast_semantics_family_calls"] = 1
            counts["new_algebra_family_calls"] = 1
        elif event == "ast_binding":
            counts["ast_binding_rows"] += 1
            counts["scope_binding_rows"] += 1
        elif event == "functional_alias_row":
            counts["functional_alias_rows"] += 1
            counts["analytic_alias_rows"] += 1
        elif event == "density_function_identity":
            counts["density_function_identity_rows"] += 1
        elif event == "jet_lift_row":
            counts["jet_lift_rows"] += 1
            counts["density_flux_rows"] += 1
        elif event == "derivative_rewrite_row":
            counts["derivative_rewrite_rows"] += 1
        elif event == "my_compatibility_row":
            counts["my_compatibility_rows"] += 1
            counts["saved_m_y_alias_rows"] += 1
        elif event == "final_binding_row":
            counts["final_binding_rows"] += 1
            counts["physical_binding_rows"] += 1
            rewritten = (
                row.get("branch") == "modified" and
                row.get("field") in {"m", "h", "k", "e", "p", "pressure"}
            )
            counts["rewritten_physical_bindings" if rewritten else "retained_physical_bindings"] += 1
        elif event == "api_refusal":
            counts["actual_negative_api_refusals"] += 1
    # A family that does not reach its completed event has no mathematical
    # proof count, even when a prefix contained all earlier proof rows.
    counts["new_mathematical_proofs"] = 0
    counts["failed_family_calls"] = 1
    return counts


def _failure_result(parent_commit: Any, reason: str, journal: list[dict[str, Any]], counts: Mapping[str, int] | None = None) -> dict[str, Any]:
    false_flags = {
        "conditional_native_O3_finite_N_definition_semantics_proved": False,
        "conditional_native_O3_density_flux_functions_proved": False,
        "conditional_native_O3_my_compatibility_proved": False,
        "native_O3_flux_radial_interface_defined": False,
        "native_O3_exp_expm1_functional_alias_proved": False,
        "native_O3_density_flux_source_identities_proved": False,
        "native_O3_saved_m_y_compatibility_proved": False,
        "unconditional_Z4_source_admission": False,
        "unconditional_ordinary_Z4_source_admission": False,
        "native_field_derivative_identity_proved": False,
        "actual_incoming_x_independence_proved": False,
        "global_phase_seam_regularity_proved": False,
        "global_across_phase_seams_proved": False,
        "oracle_source_value_identity": False,
        "physical_layer_completed": False,
        "full_goal_complete": False,
        "numeric_closure": False,
        "cap_deletion": False,
        "temporal_recursion": False,
    }
    return {
        "schema": RESULT_SCHEMA,
        "component": COMPONENT,
        "passed": False,
        "source_only": True,
        "metadata_only": False,
        "runtime_model_metadata": None,
        "parent_commit": parent_commit,
        "failure": reason,
        "proof_count": 0,
        "new_flags": false_flags,
        "actual_counts": _metadata_counts_from_journal(journal),
        "callback_event_count": len(journal),
        "callback_events": journal,
    }


def _product_coefficient(left: list[dict[str, Any]], right: list[dict[str, Any]], order: int) -> dict[str, Any]:
    summands: list[dict[str, Any]] = []
    for split in range(order + 1):
        summands.append({
            "kind": "ordinary_Z_coefficient_product",
            "left": _clone(left[split]),
            "right": _clone(right[order - split]),
            "left_factorial": FACT[split],
            "right_factorial": FACT[order - split],
        })
    return {
        "kind": "ordinary_Z_product_coefficient",
        "ordinary_Z_order": order,
        "coefficient_summands": summands,
        "coefficient_factorial": FACT[order],
        "ordinary_derivative": {"kind": "factorial_times_coefficient", "factorial": FACT[order]},
    }


def _jet_expression(field_name: str, order: int, en_rows: list[dict[str, Any]], vn_rows: list[dict[str, Any]]) -> dict[str, Any]:
    if field_name == "m":
        return {"kind": "native_graph_ref", **_clone(vn_rows[order]), "meaning": "D_Z^j V_N"}
    if field_name == "h":
        return {"kind": "native_graph_ref", **_clone(en_rows[order]), "meaning": "D_Z^j E_N"}
    product = _product_coefficient(en_rows if field_name in {"k", "p", "e"} else vn_rows,
                                   vn_rows if field_name == "k" else en_rows if field_name in {"p", "e"} else vn_rows,
                                   order)
    if field_name == "p":
        return {"kind": "ordinary_Z_scale", "factor": "1/2", "child": product}
    if field_name == "e":
        vv = _product_coefficient(vn_rows, vn_rows, order)
        return {"kind": "ordinary_Z_difference", "positive": vv, "negative": {"kind": "ordinary_Z_scale", "factor": "1/2", "child": product}}
    return product


def derive(payload: Mapping[str, Any], emit: Callable[[str, Mapping[str, Any]], None] | None = None) -> dict[str, Any]:
    """Consume the three literal wrappers and emit one durable row at a time."""
    if not isinstance(payload, Mapping):
        return _failure_result(None, "payload is not a mapping", [], {})
    spec = payload.get("input") if isinstance(payload.get("input"), Mapping) else payload
    parent_context = spec.get("frozen_context", payload.get("frozen_context", {}))
    parent_commit = parent_context.get("expected_parent_commit") if isinstance(parent_context, Mapping) else None
    if parent_commit is None:
        parent_commit = spec.get("expected_parent_commit", payload.get("expected_parent_commit"))
    if parent_commit is None:
        parent_commit = spec.get("parent_commit", payload.get("parent_commit"))
    if parent_commit is None:
        driver = spec.get("request_driver", payload.get("request_driver", {}))
        driver_context = driver.get("frozen_context", {}) if isinstance(driver, Mapping) else {}
        if isinstance(driver_context, Mapping):
            parent_commit = driver_context.get("expected_parent_commit", driver_context.get("parent_commit"))
    journal: list[dict[str, Any]] = []
    event_counts: dict[str, int] = {}

    def send(event: str, row: Mapping[str, Any]) -> None:
        event_counts[event] = event_counts.get(event, 0) + 1
        record = {"event": event, "row": _clone(dict(row))}
        journal.append(record)
        if emit is not None:
            emit(event, dict(row))

    # This is deliberately the first callback, before wrapper or source checks.
    send("family_entered", {
        "family": COMPONENT,
        "schema": INPUT_SCHEMA,
        "chart": "O3_transition",
        "ordinary_Z_orders": list(ORDERS),
        "phase_mode": "ordinary_Z",
        "radial_order": 1,
        "source_only": True,
        "numeric_evaluation": False,
        "proof_count": 0,
        "flags_all_false_until_success": True,
        "count_names": {
            "source_family_calls": "source_family_calls",
            "new_source_identity_family_calls": "new_source_identity_family_calls",
            "new_ast_semantics_family_calls": "new_ast_semantics_family_calls",
            "new_mathematical_proofs": "new_mathematical_proofs",
            "ast_binding_rows": "ast_binding_rows",
            "functional_alias_rows": "functional_alias_rows",
            "analytic_alias_rows": "analytic_alias_rows",
            "density_function_identity_rows": "density_function_identity_rows",
            "density_flux_rows": "density_flux_rows",
            "jet_lift_rows": "jet_lift_rows",
            "derivative_rewrite_rows": "derivative_rewrite_rows",
            "my_compatibility_rows": "my_compatibility_rows",
            "saved_m_y_alias_rows": "saved_m_y_alias_rows",
            "final_binding_rows": "final_binding_rows",
            "physical_binding_rows": "physical_binding_rows",
            "scope_binding_rows": "scope_binding_rows",
            "rewritten_physical_bindings": "rewritten_physical_bindings",
            "retained_physical_bindings": "retained_physical_bindings",
            "actual_negative_api_refusals": "actual_negative_api_refusals",
            "old_math_calls": "old_math_calls",
            "old_provider_calls": "old_provider_calls",
            "numeric_evaluation_calls": "numeric_evaluation_calls",
            "N_materialization_calls": "N_materialization_calls",
            "source_value_reducer_calls": "source_value_reducer_calls",
            "failed_family_calls": "failed_family_calls",
            "old_jet_recurrence_calls": "old_jet_recurrence_calls",
        },
    })

    counts: dict[str, int] = {
        "source_family_calls": 0,
        "new_source_identity_family_calls": 0,
        "new_ast_semantics_family_calls": 0,
        "new_algebra_family_calls": 0,
        "new_mathematical_proofs": 0,
        "ast_binding_rows": 0,
        "functional_alias_rows": 0,
        "analytic_alias_rows": 0,
        "density_function_identity_rows": 0,
        "density_flux_rows": 0,
        "jet_lift_rows": 0,
        "derivative_rewrite_rows": 0,
        "my_compatibility_rows": 0,
        "saved_m_y_alias_rows": 0,
        "final_binding_rows": 0,
        "physical_binding_rows": 0,
        "scope_binding_rows": 0,
        "rewritten_physical_bindings": 0,
        "retained_physical_bindings": 0,
        "actual_negative_api_refusals": 0,
        "old_math_calls": 0,
        "old_provider_calls": 0,
        "numeric_evaluation_calls": 0,
        "N_materialization_calls": 0,
        "source_value_reducer_calls": 0,
        "cached_algebra_replay_calls": 0,
        "failed_family_calls": 0,
        "old_jet_recurrence_calls": 0,
    }
    try:
        _require(spec.get("schema") == INPUT_SCHEMA, "input schema changed")
        _require(spec.get("component", COMPONENT) == COMPONENT, "component changed")
        wrappers = _wrapper_map(spec, payload)
        source, lemma, saved, parent_payload = _validate_packets(wrappers)
        counts["source_family_calls"] = 1
        counts["new_source_identity_family_calls"] = 1
        counts["new_ast_semantics_family_calls"] = 1
        counts["new_algebra_family_calls"] = 1

        opaque = source.get("opaque_source_parameters")
        n_data = opaque.get("N") if isinstance(opaque, Mapping) else None
        inv_data = opaque.get("inverse_N_graph_port") if isinstance(opaque, Mapping) else None
        _require(isinstance(n_data, Mapping) and n_data.get("node") == 18 and n_data.get("materialized") is False,
                 "opaque N certificate is not retained")
        _require(isinstance(inv_data, Mapping) and inv_data.get("node") == 99411 and inv_data.get("materialized") is False,
                 "inverse-N graph port changed")
        target_n = source.get("target_certificates_data_only", {}).get("N_root_18")
        _require(isinstance(target_n, Mapping) and target_n.get("certificate", {}).get("root") == 18,
                 "N_root_18 data certificate is missing")

        source_ast = source.get("source_ast")
        _require(isinstance(source_ast, Mapping) and isinstance(source_ast.get("finite_N_functions"), Mapping),
                 "finite_N source AST witness is missing")
        _require(set(source_ast) == {
            "BridgeAlgebra.exp", "JetAlgebra.add", "JetAlgebra.mul", "JetAlgebra.scale",
            "background_density", "finite_N_functions",
        }, "source AST witness keys changed")
        finite_source = source_ast["finite_N_functions"].get("source_text")
        compiler = _SourceCompiler(finite_source)
        typed_ir = compiler.compile()  # deterministic AST data, not execution
        fields = compiler.return_fields
        send("ast_binding", {
            "family": COMPONENT,
            "function": typed_ir["function_name"],
            "typed_function_ir": typed_ir,
            "ordinary_Z_orders": list(ORDERS),
            "source_ast_only": True,
            "numeric_evaluation": False,
        })
        counts["ast_binding_rows"] = 1
        counts["scope_binding_rows"] = 1

        increments = fields["five_signed_density_increments"]
        x_term = fields["A_over_N"]
        en_term, vn_term = fields["E_N"], fields["V_N"]
        de_term, dv_term = fields["delta_E"], fields["delta_V"]
        _require(x_term.kind == "jet" and x_term.op == "scale" and len(x_term.children) == 2,
                 "A_over_N typed definition changed")
        _require(en_term.kind == "jet" and en_term.op == "mul" and len(en_term.children) == 2,
                 "E_N typed definition changed")
        _require(vn_term.kind == "jet" and vn_term.op == "add" and len(vn_term.children) == 2,
                 "V_N typed definition changed")
        _require(de_term.kind == "jet" and de_term.op == "mul" and len(de_term.children) == 2,
                 "delta_E typed definition changed")
        _require(dv_term.kind == "jet" and dv_term.op == "scale" and len(dv_term.children) == 2,
                 "delta_V typed definition changed")
        _require(de_term.children[0].kind == "source_ref" and de_term.children[0].data.get("key") == "E",
                 "delta_E root binding changed")
        _require(dv_term.children[0].kind == "source_ref" and dv_term.children[0].data.get("key") == "B_over_Pstar",
                 "delta_V loop binding changed")
        _require(en_term.children[0].kind == "source_ref" and en_term.children[0].data.get("key") == "E",
                 "E_N root binding changed")
        _require(vn_term.children[0].kind == "source_ref" and vn_term.children[0].data.get("key") == "V",
                 "V_N root binding changed")
        _require(en_term.children[1].kind == "jet" and en_term.children[1].op == "exp",
                 "E_N exponential term changed")
        exp_term = en_term.children[1]
        _require(de_term.children[1].kind == "jet_blend" and de_term.children[1].children[0].kind == "analytic_unary" and
                 de_term.children[1].children[0].op == "expm1",
                 "finite_N value/higher expm1 splice changed")
        expm1_value = de_term.children[1].children[0]
        _require(len(exp_term.children) == 1 and len(expm1_value.children) == 1 and
                 expm1_value.children[0].kind == "jet_value" and
                 len(expm1_value.children[0].children) == 1 and
                 _same_term(exp_term.children[0], expm1_value.children[0].children[0]),
                 "exp and expm1 do not share the same typed argument")
        _require(_same_term(vn_term.children[1], dv_term), "V_N does not use the same s term")

        # Structural term-definition matching is intentionally local: it
        # checks the five source definitions without constructing a generic
        # ring or reordering terms.
        roots_e = _term("source_ref", namespace="roots", key="E")
        roots_v = _term("source_ref", namespace="roots", key="V")
        loop_b = _term("source_ref", namespace="loop", key="B_over_Pstar")
        expected_k = _term("jet", op="add", children=(
            _term("jet", op="mul", children=(roots_v, de_term)),
            _term("jet", op="mul", children=(roots_e, dv_term)),
            _term("jet", op="mul", children=(de_term, dv_term)),
        ))
        expected_e = _term("jet", op="sub", children=(
            _term("jet", op="add", children=(
                _term("jet", op="scale", children=(
                    _term("jet", op="mul", children=(roots_v, dv_term)), _term("literal", value=2)
                )),
                _term("jet", op="mul", children=(dv_term, dv_term)),
            )),
            _term("jet", op="add", children=(
                _term("jet", op="mul", children=(roots_e, de_term)),
                _term("jet", op="scale", children=(
                    _term("jet", op="mul", children=(de_term, de_term)), _term("literal", value="1/2")
                )),
            )),
        ))
        expected_p = _term("jet", op="add", children=(
            _term("jet", op="mul", children=(roots_e, de_term)),
            _term("jet", op="scale", children=(
                _term("jet", op="mul", children=(de_term, de_term)), _term("literal", value="1/2")
            )),
        ))
        _require(_same_term(increments.data["m"], dv_term) and _same_term(increments.data["h"], de_term),
                 "linear density increment definitions changed")
        _require(_same_term(increments.data["k"], expected_k) and
                 _same_term(increments.data["e"], expected_e) and
                 _same_term(increments.data["p"], expected_p),
                 "nonlinear density increment definition changed")
        _require(dv_term.children[0].kind == "source_ref" and
                 _same_term(dv_term.children[0], loop_b), "delta_V loop source changed")

        constructor_exp = _row_refs(source, "finite_N_constructor_rows", "exp_A_over_N")
        constructor_expm1 = _row_refs(source, "finite_N_constructor_rows", "expm1_A_over_N")
        x_rows = _row_refs(source, "finite_N", "A_over_N")
        en_rows = _row_refs(source, "finite_N", "E_N")
        vn_rows = _row_refs(source, "finite_N", "V_N")
        de_rows = _row_refs(source, "finite_N", "delta_E")
        dv_rows = _row_refs(source, "finite_N", "delta_V")
        density_rows = {name: _row_refs(source, "five_signed_density_increments", name) for name in FIELDS}
        background_e = _row_refs(source, "background_fields", "E")
        background_v = _row_refs(source, "background_fields", "V")
        loop_a = _row_refs(source, "loop", "A")
        loop_b = _row_refs(source, "loop", "B_over_Pstar")
        modified_rows = {name: _row_refs(source, "modified_fields", name) for name in ("E", "V", "m")}
        for order in ORDERS:
            _require(modified_rows["E"][order]["node"] == en_rows[order]["node"],
                     f"modified E source port changed at ordinary Z{order}")
            _require(modified_rows["V"][order]["node"] == vn_rows[order]["node"],
                     f"modified V source port changed at ordinary Z{order}")
        saved_my_rows = _row_refs(source, "saved_total_radial", "m_y")
        selected_row_count = 0
        for group in source.get("literal_rows", {}).values():
            if isinstance(group, Mapping):
                for entry in group.values():
                    if isinstance(entry, Mapping):
                        selected_row_count += int(entry.get("selected_row_count", 0))
        _require(selected_row_count == 120, "raw packet selected-row count changed")
        graph = source.get("source_graph")
        _require(isinstance(graph, Mapping) and graph.get("head_count") == 84 and
                 isinstance(graph.get("heads"), Mapping) and len(graph["heads"]) == 84,
                 "source graph head contract changed")
        constructor = graph.get("constructor_skeleton")
        _require(isinstance(constructor, Mapping) and isinstance(constructor.get("records"), Mapping) and
                 len(constructor["records"]) == 189, "constructor skeleton record count changed")

        alias_rows = [
            {
                "alias": "expM1_value",
                "equation": "expm1(x)=exp(x)-1",
                "argument": _clone(x_rows[0]),
                "typed_argument": _term_dict(x_term),
                "loop_A_ref": _clone(loop_a[0]),
                "exp_ref": _clone(constructor_exp[0]),
                "expm1_ref": _clone(constructor_expm1[0]),
                "higher_exp_refs": _clone(constructor_exp[1:]),
                "higher_expm1_refs": _clone(constructor_expm1[1:]),
                "higher_rows_share_exp_derivative_sources": [
                    constructor_exp[index]["node"] == constructor_expm1[index]["node"] for index in range(1, 5)
                ],
                "expm1_value_row_distinct": constructor_exp[0]["node"] != constructor_expm1[0]["node"],
                "alias_has_no_A_division": True,
                "A_zero_branch": "x=0 gives exp(0)-1=0 without division by A",
                "functional": True,
            },
            {
                "alias": "E_N_equals_E_plus_r",
                "equation": "E_N=E+r",
                "E_ref": _clone(background_e[0]),
                "E_N_ref": _clone(en_rows[0]),
                "r_ref": _clone(de_rows[0]),
                "exp_ref": _clone(constructor_exp[0]),
                "expm1_ref": _clone(constructor_expm1[0]),
                "functional": True,
            },
            {
                "alias": "V_N_equals_V_plus_s",
                "equation": "V_N=V+s",
                "V_ref": _clone(background_v[0]),
                "V_N_ref": _clone(vn_rows[0]),
                "s_ref": _clone(dv_rows[0]),
                "loop_B_over_Pstar_ref": _clone(loop_b[0]),
                "functional": True,
            },
        ]
        for row in alias_rows:
            send("functional_alias_row", row)
        counts["functional_alias_rows"] = len(alias_rows)
        counts["analytic_alias_rows"] = len(alias_rows)

        schema = lemma["abstract_density_schema"]
        base_schema = schema["base_density_functions"]
        increment_schema = schema["flux_or_density_increment"]
        expected_increment_terms = {
            "m": "s", "h": "r", "k": "V*r+E*s+r*s",
            "e": "2*V*s+s^2-E*r-r^2/2", "p": "E*r+r^2/2",
        }
        density_identity_rows: list[dict[str, Any]] = []
        for field_name in FIELDS:
            term = increments.data[field_name]
            _require(term.kind == "jet", f"increment {field_name} is not a typed jet")
            identity = {
                "field": field_name,
                "base_function": base_schema[field_name],
                "increment_function": increment_schema[field_name],
                "expected_structural_increment": expected_increment_terms[field_name],
                "source_increment_refs": _clone(density_rows[field_name]),
                "identity": f"F_{field_name}(E_N,V_N)-F_{field_name}(E,V)={increment_schema[field_name]}",
                "ordinary_Z_orders": list(ORDERS),
                "structural_term_definition_match": True,
                "generic_pointwise_lemma_consumed_as_data": True,
                "numeric_evaluation": False,
            }
            density_identity_rows.append(identity)
            send("density_function_identity", identity)
        counts["density_function_identity_rows"] = len(density_identity_rows)

        jet_rows: list[dict[str, Any]] = []
        for field_name in FIELDS:
            for order in ORDERS:
                expression = _jet_expression(field_name, order, en_rows, vn_rows)
                row = {
                    "field": field_name,
                    "ordinary_Z_order": order,
                    "phase_mode": "ordinary_Z",
                    "radial_order": 1,
                    "expression": expression,
                    "factorial_table_0_through_4": list(FACT),
                    "coefficient_rule": "jet_coefficient_j=f_j/j!; product=sum_l f_l/l!*g_(j-l)/(j-l)!",
                    "ordinary_derivative_rule": "D_Z^j=j!*coefficient_j",
                    "conditional_source_jet_semantics": True,
                    "unconditional_Z4_source_admission": False,
                    "numeric_evaluation": False,
                }
                jet_rows.append(row)
                send("jet_lift_row", row)
        counts["jet_lift_rows"] = len(jet_rows)
        counts["density_flux_rows"] = len(jet_rows)
        full_flux_expressions = {
            (row["field"], row["ordinary_Z_order"]): row["expression"] for row in jet_rows
        }
        _require(len(full_flux_expressions) == 25, "fresh full density jet lookup is incomplete")

        rates = _rate_map(parent_payload)
        parent_modified = parent_payload.get("modified_history_rows")
        _require(isinstance(parent_modified, list) and len(parent_modified) == 25,
                 "saved modified history rows changed")
        rewrite_rows: list[dict[str, Any]] = []
        rewrite_by_key: dict[tuple[str, int], dict[str, Any]] = {}
        for parent_index, parent_row in enumerate(parent_modified):
            key = str(parent_row.get("history_key"))
            order = parent_row.get("ordinary_Z_order")
            _require(key in FIELDS and isinstance(order, int) and order in ORDERS,
                     "saved modified history row key changed")
            same_h_node = parent_row.get("modified_node")
            _require(type(same_h_node) is int and same_h_node >= 0, "saved nonzero incoming history node missing")
            source_increment_ref = density_rows[key][order]
            full_flux_expression = full_flux_expressions[(key, order)]
            same_history_ref = _native_ref(
                same_h_node, f"payload.modified_history_rows[{parent_index}].modified_node", order
            )
            rewrite = {
                "branch": "modified",
                "history_key": key,
                "ordinary_Z_order": order,
                "rate_lambda": rates[key],
                "new_flux_ref": _clone(full_flux_expression),
                "source_increment_ref": _clone(source_increment_ref),
                "same_history_ref": _clone(same_history_ref),
                "source_expression": {
                    "kind": "typed_source_derivative_rewrite",
                    "operator": "D_y",
                    "equation": "newflux_j-lambda*sameH_j",
                    "newflux": _clone(full_flux_expression),
                    "sameH": _clone(same_history_ref),
                    "rate_lambda": rates[key],
                    "evaluated": False,
                },
                "incoming_nonzero_histories_retained": bool(parent_row.get("premises", {}).get("nonzero_incoming_histories_retained", True)),
                "parent_data_pointer": f"payload.modified_history_rows[{parent_index}]",
                "phase_seam_regularity_proved": False,
                "numeric_evaluation": False,
            }
            _require(rewrite["incoming_nonzero_histories_retained"], "nonzero incoming history premise was dropped")
            rewrite_rows.append(rewrite)
            rewrite_by_key[(key, order)] = rewrite
            send("derivative_rewrite_row", rewrite)
        counts["derivative_rewrite_rows"] = len(rewrite_rows)

        # Prove the five saved m_y compatibility rows only from actual sum and
        # negative-child records in the packet.
        compatibility_rows: list[dict[str, Any]] = []
        for order, sum_ref in enumerate(saved_my_rows):
            sum_head = _head(graph, sum_ref["node"])
            children = sum_head.get("ordered_children")
            _require(sum_head.get("operation") == "sum" and isinstance(children, list) and len(children) == 2,
                     "saved m_y head is not a binary sum")
            v_ref = vn_rows[order]
            _require(children[0] == v_ref["node"], "saved m_y first child is not V_N")
            negative_node = children[1]
            negative_head = _head(graph, negative_node)
            m_ref = modified_rows["m"][order]
            _require(negative_head.get("operation") == "negative" and
                     negative_head.get("ordered_children") == [m_ref["node"]],
                     "saved m_y negative child does not point to modified m_j")
            raw_record = negative_head.get("raw_record")
            _require(isinstance(raw_record, Mapping) and raw_record.get("operation") == "negative" and
                     raw_record.get("argument") == m_ref["node"],
                     "saved m_y negative raw record changed")
            row = {
                "ordinary_Z_order": order,
                "equation": "m_y=V_N_j+(-modified_m_j)",
                "sum_head": _native_ref(sum_ref["node"], sum_ref["source_path"], order),
                "V_N_ref": _clone(v_ref),
                "negative_child": _native_ref(negative_node, f"saved_total_radial.m_y[{order}].second_child", order),
                "negative_argument": _clone(m_ref),
                "negative_head_operation": "negative",
                "raw_record_argument_matches": True,
                "strict_actual_native_head": True,
                "numeric_evaluation": False,
            }
            compatibility_rows.append(row)
            send("my_compatibility_row", row)
        counts["my_compatibility_rows"] = len(compatibility_rows)
        counts["saved_m_y_alias_rows"] = len(compatibility_rows)

        parent_physical = parent_payload.get("physical_radius_rows")
        _require(isinstance(parent_physical, list) and len(parent_physical) == 90,
                 "saved physical radius rows changed")
        all_keys: set[tuple[str, str, int]] = set()
        for row in parent_physical:
            key = (row.get("branch"), row.get("field"), row.get("ordinary_Z_order"))
            _require(key not in all_keys, "saved physical radius row keys are not unique")
            all_keys.add(key)
        native_radius_cert = parent_payload.get("cached_data", {}).get("nativeR_certificate", {})
        native_radius_node = native_radius_cert.get("root")
        _require(native_radius_node == 241787, "native radius certificate root changed")
        frame_rows = parent_payload.get("cached_data", {}).get("frame_specialization_rows", [])
        _require(isinstance(frame_rows, list) and any(
            isinstance(item, Mapping) and isinstance(item.get("branch_operator_manifest"), Mapping) and
            item["branch_operator_manifest"].get("frame", {}).get("coordinate_node") == 241536
            for item in frame_rows
        ), "native chart coordinate x241536 is not retained as source data")

        final_rows: list[dict[str, Any]] = []
        retained_rows: list[dict[str, Any]] = []
        fresh_rows: list[dict[str, Any]] = []
        for parent_row in parent_physical:
            branch, field_name, order = parent_row.get("branch"), parent_row.get("field"), parent_row.get("ordinary_Z_order")
            _require(branch in BRANCHES and field_name in ("E", "V", "m", "h", "k", "e", "p", "P0", "pressure") and order in ORDERS,
                     "saved physical binding key changed")
            rewritten = branch == "modified" and field_name in {"m", "h", "k", "e", "p", "pressure"}
            if not rewritten:
                out_row = _clone(parent_row)
                retained_rows.append(out_row)
                final_rows.append(out_row)
                send("final_binding_row", out_row)
                continue
            if field_name == "pressure":
                density_field = "p"
                rewrite = rewrite_by_key[(density_field, order)]
                numerator = _clone(rewrite["source_expression"])
                source_expression: dict[str, Any] = {
                    "kind": "typed_pressure_radial_source",
                    "operator": "D_R",
                    "equation": "D_R(P0)=0 and D_R(p)=(newflux_j-lambda*sameH_j)/R_native",
                    "P0_radial_zero_only": True,
                    "ordinary_Z_zero_not_implied": True,
                    "numerator": _clone(numerator),
                    "denominator": _native_ref(native_radius_node, "payload.cached_data.nativeR_certificate", None, name="R_native"),
                    "P0_radial_zero_term": {
                        "kind": "declared_zero_radial_source",
                        "field": "P0",
                        "radial_derivative": "0",
                        "pointer": f"payload.source_proof_rows.background.P0.{order}.source_expression",
                    },
                    "parent_P0_pointer": f"payload.source_proof_rows.background.P0.{order}.source_expression",
                    "fresh_source_proof_namespace": f"payload.fresh_source_proof_rows.modified.pressure.{order}",
                }
            else:
                rewrite = rewrite_by_key[(field_name, order)]
                numerator = _clone(rewrite["source_expression"])
                source_expression = {
                    "kind": "typed_first_radial_source",
                    "operator": "D_R",
                    "equation": "(newflux_j-lambda*sameH_j)/R_native",
                    "numerator": _clone(numerator),
                    "denominator": _native_ref(native_radius_node, "payload.cached_data.nativeR_certificate", None, name="R_native"),
                    "fresh_source_proof_namespace": f"payload.fresh_source_proof_rows.modified.{field_name}.{order}",
                }
            out_row = {
                "branch": branch,
                "field": field_name,
                "ordinary_Z_order": order,
                "phase_mode": "ordinary_Z",
                "radial_order": 1,
                "operator": "D_R",
                "native_radius": _native_ref(native_radius_node, "payload.cached_data.nativeR_certificate", None, name="R_native"),
                "source_expression": source_expression,
                "source_numerator": _clone(numerator),
                "parent_binding_kind": "replaced_modified_history_or_pressure",
                "parent_relative_pointer_namespace": "payload.cached_data/source_proof_rows",
                "actual_field_derivative_identity_proved": False,
                "phase_seam_regularity_proved": False,
                "numeric_evaluation": False,
            }
            fresh_rows.append(out_row)
            final_rows.append(out_row)
            send("final_binding_row", out_row)
        _require(len(final_rows) == 90 and len(retained_rows) == 60 and len(fresh_rows) == 30,
                 "final physical binding multiplicity changed")
        counts["retained_physical_bindings"] = len(retained_rows)
        counts["rewritten_physical_bindings"] = len(fresh_rows)
        counts["final_binding_rows"] = len(final_rows)
        counts["physical_binding_rows"] = len(final_rows)

        # Four negative controls are exercised exactly once, after all positive
        # source rows have been emitted.
        negative_controls = [
            ("held_phase", "background", "m", 0, 1, "held_phase"),
            ("second_radial", "background", "m", 0, 2, "ordinary_Z"),
            ("wrong_branch", "comparison", "m", 0, 1, "ordinary_Z"),
            ("Z5", "background", "m", 5, 1, "ordinary_Z"),
        ]
        for case, branch, field_name, order, radial_order, phase_mode in negative_controls:
            try:
                select_radius_binding(branch, field_name, order, radial_order, phase_mode)
            except ValueError as exc:
                counts["actual_negative_api_refusals"] += 1
                send("api_refusal", {"case": case, "refused": True, "reason": str(exc), "numeric_evaluation": False})
            else:
                raise _ContractError(f"negative selector control {case} was admitted")
        _require(counts["actual_negative_api_refusals"] == 4, "negative selector count changed")

        flags = {
            "conditional_native_O3_finite_N_definition_semantics_proved": True,
            "conditional_native_O3_density_flux_functions_proved": True,
            "conditional_native_O3_my_compatibility_proved": True,
            "native_O3_flux_radial_interface_defined": True,
            "native_O3_exp_expm1_functional_alias_proved": True,
            "native_O3_density_flux_source_identities_proved": True,
            "native_O3_saved_m_y_compatibility_proved": True,
            "unconditional_Z4_source_admission": False,
            "unconditional_ordinary_Z4_source_admission": False,
            "native_field_derivative_identity_proved": False,
            "actual_incoming_x_independence_proved": False,
            "global_phase_seam_regularity_proved": False,
            "global_across_phase_seams_proved": False,
            "oracle_source_value_identity": False,
            "physical_layer_completed": False,
            "full_goal_complete": False,
            "numeric_closure": False,
            "cap_deletion": False,
            "temporal_recursion": False,
        }
        completed_counts = dict(counts)
        completed_counts["new_mathematical_proofs"] = 1
        completed_row = {
            "family": COMPONENT,
            "proof_count": 1,
            "ordinary_Z_orders": list(ORDERS),
            "conditional_source_semantics": True,
            "unconditional_Z4_source_admission": False,
            "flags": _clone(flags),
            "actual_counts": completed_counts,
            "source_only": True,
            "numeric_evaluation": False,
        }
        send("family_completed", completed_row)
        counts["new_mathematical_proofs"] = 1
        output_payload = _clone(dict(parent_payload))
        # The parent-relative namespaces remain at their original payload keys;
        # all fresh proof rows are in a separate namespace.
        output_payload.update({
            "typed_function_ir": typed_ir,
            "functional_alias_rows": alias_rows,
            "density_function_identity_rows": density_identity_rows,
            "jet_lift_rows": jet_rows,
            "derivative_rewrite_rows": rewrite_rows,
            "my_compatibility_rows": compatibility_rows,
            "physical_radius_rows": final_rows,
            "physical_radius_bindings": {
                branch: {field_name: {str(row["ordinary_Z_order"]): _clone(row) for row in final_rows
                                       if row.get("branch") == branch and row.get("field") == field_name}
                         for field_name in ("E", "V", "m", "h", "k", "e", "p", "P0", "pressure")}
                for branch in BRANCHES
            },
            "fresh_source_proof_rows": {
                "modified": {
                    field_name: {str(row["ordinary_Z_order"]): {"proof_kind": "fresh_D_R_binding", "ordinary_Z_order": row["ordinary_Z_order"],
                                              "source_expression": _clone(row["source_expression"])}
                                  for row in fresh_rows if row.get("field") == field_name and row.get("ordinary_Z_order") in ORDERS}
                    for field_name in ("m", "h", "k", "e", "p", "pressure")
                }
            },
            "retained_parent_physical_bindings": retained_rows,
            "fresh_rewritten_physical_bindings": fresh_rows,
            "retained_parent_proof_namespace": "payload.cached_data/source_proof_rows",
            "fresh_proof_namespace": "payload.fresh_source_proof_rows",
            "source_scope": {
                "chart": "O3_transition",
                "fields": ["E", "V", "m", "h", "k", "e", "p", "P0", "pressure"],
                "branches": list(BRANCHES),
                "ordinary_Z_orders": list(ORDERS),
                "radial_order": 1,
                "phase_mode": "ordinary_Z",
                "conditional_source_parameter_semantics": True,
                "native_radius_node": native_radius_node,
                "opaque_N_node": 18,
                "inverse_N_node": 99411,
            },
        })
        result = {
            "schema": RESULT_SCHEMA,
            "component": COMPONENT,
            "passed": True,
            "source_only": True,
            "metadata_only": False,
            "conditional": True,
            "conditional_on_declared_source_atoms": True,
            "runtime_model_metadata": None,
            "parent_commit": parent_commit,
            "source_wrapper_pins": {name: {"sha256": value["sha256"]} for name, value in wrappers.items()},
            "source_archive_sha256": source.get("allowed_inputs", {}).get("native_archive", {}).get("raw_sha256"),
            "proof_count": 1,
            "new_flags": flags,
            "actual_counts": counts,
            "callback_event_count": len(journal),
            "callback_event_counts": dict(event_counts),
            "callback_events": journal,
            "payload": output_payload,
            "failure": None,
        }
        _require(len(journal) == 160, "callback row count changed")
        return result
    except Exception as exc:
        counts["failed_family_calls"] += 1
        return _failure_result(parent_commit, str(exc), journal, counts)


__all__ = ["COMPONENT", "INPUT_SCHEMA", "RESULT_SCHEMA", "derive", "select_radius_binding"]
