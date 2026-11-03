"""AST-and-graph candidate witness for the original core gauge integral.

This module is deliberately a small, non-executing boundary.  It consumes an
explicit serialized node list and raw source text supplied by its caller.  It
does not import a producer, evaluate source, look up a node by name, install a
binder sidecar, or inspect a live graph object.  A passing result therefore
means only that the two supplied artifacts support one candidate binding.

The source contract is intentionally narrow.  In one lexical scope it must
contain the ordered chain::

    v = g.symbol("<binder name>")
    gv = g.node("function_substitution", value=v.node, ...)
    G = controls.integral(g, gv, g.nodes[v.node]["name"], self.parameters["anchor"], z, ...)

The graph contract then checks the explicit integral row, its explicit
``function_substitution`` integrand, and the exact ``value`` edge to the
caller-supplied binder node.  This is a source-plus-serialized-graph
candidate witness; it is not a live-runtime owner proof and it does not prove
that an archived owner is the same object.
"""

from __future__ import annotations

import ast
import hashlib
import json
from collections.abc import Iterable, Mapping, Sequence
from typing import Any


SCHEMA = "n1_core_gauge_binding_witness_v2"
RESULT_SCHEMA = SCHEMA + "_result"

_REQUIRED_INTEGRAL_FLAGS = (
    "exact_function_integral",
    "numerical_value_not_installed",
    "original_H_root_anchor",
    "anchor_is_exact_root_not_disk_center",
    "actual_anchored_gauge_integral",
)

_NESTED_SCOPE_TYPES = (
    ast.FunctionDef,
    ast.AsyncFunctionDef,
    ast.ClassDef,
    ast.Lambda,
    ast.ListComp,
    ast.SetComp,
    ast.DictComp,
    ast.GeneratorExp,
)


class GaugeBindingWitnessError(ValueError):
    """Raised by the strict builder when the supplied evidence is rejected."""

    def __init__(self, message: str, result: Mapping[str, Any] | None = None) -> None:
        super().__init__(message)
        self.result = dict(result) if isinstance(result, Mapping) else None


# A descriptive alias is useful to callers that use the word validation.
WitnessValidationError = GaugeBindingWitnessError


def _is_int(value: Any) -> bool:
    return type(value) is int


def _is_true(value: Any) -> bool:
    return value is True


def _copy_json(value: Any) -> Any:
    """Copy JSON-shaped data without invoking arbitrary user methods."""

    if isinstance(value, Mapping):
        return {str(key): _copy_json(item) for key, item in value.items()}
    if isinstance(value, list):
        return [_copy_json(item) for item in value]
    if isinstance(value, tuple):
        return [_copy_json(item) for item in value]
    return value


def _loc(node: ast.AST | None) -> dict[str, int] | None:
    if node is None or not hasattr(node, "lineno"):
        return None
    return {
        "line": int(getattr(node, "lineno")),
        "column": int(getattr(node, "col_offset", 0)),
    }


def _position(node: ast.AST) -> tuple[int, int]:
    return (
        int(getattr(node, "lineno", 0)),
        int(getattr(node, "col_offset", 0)),
    )


def _name(node: ast.AST | None, wanted: str | None = None) -> bool:
    return isinstance(node, ast.Name) and (wanted is None or node.id == wanted)


def _attribute(node: ast.AST | None, owner: str, attribute: str) -> bool:
    return (
        isinstance(node, ast.Attribute)
        and node.attr == attribute
        and _name(node.value, owner)
    )


def _constant(node: ast.AST | None, value: Any) -> bool:
    return isinstance(node, ast.Constant) and node.value == value


def _subscript_slice(node: ast.Subscript) -> ast.AST | None:
    # ``ast.Index`` was used by old Python versions; current versions expose
    # the slice directly.  Keeping this compatibility local makes the exact
    # graph-row expression below independent of the interpreter version.
    value = node.slice
    if isinstance(value, ast.Index):  # pragma: no cover - old Python shape
        value = value.value
    return value


def _g_nodes_name(node: ast.AST | None, binder_local: str) -> bool:
    """Match ``g.nodes[v.node]['name']`` with both handles fixed."""

    if not isinstance(node, ast.Subscript):
        return False
    if not _constant(_subscript_slice(node), "name"):
        return False
    outer = node.value
    if not isinstance(outer, ast.Subscript):
        return False
    if not _attribute(outer.value, "g", "nodes"):
        return False
    return _attribute(_subscript_slice(outer), binder_local, "node")


def _self_anchor(node: ast.AST | None) -> bool:
    """Match the source constructor's explicit ``self.parameters['anchor']``."""

    if not isinstance(node, ast.Subscript):
        return False
    if not _constant(_subscript_slice(node), "anchor"):
        return False
    return _attribute(node.value, "self", "parameters")


def _node_reference(node: ast.AST | None, *, local: str | None = None) -> bool:
    if not isinstance(node, ast.Attribute) or node.attr != "node":
        return False
    if local is not None:
        return _name(node.value, local)
    # A source expression such as ``gradient.value.node`` is a pure attribute
    # chain.  Do not accept calls, subscripts, or other computed receivers.
    receiver: ast.AST = node.value
    while isinstance(receiver, ast.Attribute):
        receiver = receiver.value
    return isinstance(receiver, ast.Name)


def _simple_node_reference(node: ast.AST | None) -> bool:
    """Match a local variable handle such as ``z.node`` exactly."""

    return isinstance(node, ast.Attribute) and node.attr == "node" and isinstance(node.value, ast.Name)


def _call(node: ast.AST | None, owner: str, attribute: str) -> ast.Call | None:
    if not isinstance(node, ast.Call) or not _attribute(node.func, owner, attribute):
        return None
    return node


def _keyword_values(call: ast.Call) -> dict[str, ast.AST] | None:
    """Return keyword values, rejecting ``**mapping`` and duplicate names."""

    result: dict[str, ast.AST] = {}
    for keyword in call.keywords:
        if keyword.arg is None or keyword.arg in result:
            return None
        result[keyword.arg] = keyword.value
    return result


def _assignment_name(node: ast.AST) -> str | None:
    if not isinstance(node, ast.Assign) or len(node.targets) != 1:
        return None
    target = node.targets[0]
    return target.id if isinstance(target, ast.Name) else None


def _target_names(node: ast.AST | None) -> set[str]:
    """Collect names introduced by one assignment-like target."""

    if isinstance(node, ast.Name):
        return {node.id}
    if isinstance(node, (ast.Tuple, ast.List)):
        result: set[str] = set()
        for item in node.elts:
            result.update(_target_names(item))
        return result
    # Attribute/subscript stores do not rebind a lexical local name.
    return set()


def _binding_events(scope: ast.AST) -> dict[str, list[ast.AST]]:
    """Collect local rebinding events while excluding nested lexical scopes."""

    events: dict[str, list[ast.AST]] = {}
    for node in _scope_nodes(scope):
        targets: Iterable[ast.AST] = ()
        if isinstance(node, ast.Assign):
            targets = node.targets
        elif isinstance(node, (ast.AnnAssign, ast.AugAssign, ast.NamedExpr)):
            targets = (node.target,)
        elif isinstance(node, (ast.For, ast.AsyncFor)):
            targets = (node.target,)
        elif isinstance(node, (ast.With, ast.AsyncWith)):
            targets = tuple(
                item.optional_vars
                for item in node.items
                if item.optional_vars is not None
            )
        elif isinstance(node, ast.ExceptHandler) and node.name:
            # ``ExceptHandler.name`` is a string, not an AST target.
            events.setdefault(node.name, []).append(node)
            continue
        for target in targets:
            for name in _target_names(target):
                events.setdefault(name, []).append(node)
    return events


def _scope_nodes(scope: ast.AST) -> Iterable[ast.AST]:
    """Yield one scope's nodes, stopping at nested function/class scopes."""

    pending: list[tuple[ast.AST, bool]] = [(scope, True)]
    while pending:
        node, is_root = pending.pop()
        yield node
        # The root is the scope being inspected.  Any nested function or class
        # is a separate lexical scope and none of its body belongs here.
        if not is_root and isinstance(node, _NESTED_SCOPE_TYPES):
            continue
        for child in ast.iter_child_nodes(node):
            if isinstance(child, _NESTED_SCOPE_TYPES):
                continue
            pending.append((child, False))


def _function_arguments(node: ast.FunctionDef | ast.AsyncFunctionDef) -> set[str]:
    args = node.args
    result = {item.arg for item in args.posonlyargs}
    result.update(item.arg for item in args.args)
    result.update(item.arg for item in args.kwonlyargs)
    if args.vararg is not None:
        result.add(args.vararg.arg)
    if args.kwarg is not None:
        result.add(args.kwarg.arg)
    return result


def _scopes(tree: ast.Module) -> list[tuple[str, ast.AST]]:
    """Return module and every function as separate lexical scopes."""

    result: list[tuple[str, ast.AST]] = [("<module>", tree)]
    for node in ast.walk(tree):
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            result.append((node.name, node))
    return result


def _exact_symbol_assignment(
    node: ast.AST,
    *,
    binder_local: str,
    binder_name: str,
) -> ast.Assign | None:
    if _assignment_name(node) != binder_local:
        return None
    assert isinstance(node, ast.Assign)
    call = _call(node.value, "g", "symbol")
    if call is None or len(call.args) != 1 or call.keywords:
        return None
    return node if _constant(call.args[0], binder_name) else None


def _exact_substitution_assignment(
    node: ast.AST,
    *,
    binder_local: str,
    substitution_local: str,
) -> tuple[ast.Assign, str] | None:
    if _assignment_name(node) != substitution_local:
        return None
    assert isinstance(node, ast.Assign)
    call = _call(node.value, "g", "node")
    if call is None or len(call.args) != 1 or not _constant(call.args[0], "function_substitution"):
        return None
    keywords = _keyword_values(call)
    if keywords is None:
        return None
    required = {"expression", "variable", "value", "Z_independent_substitution"}
    if not required.issubset(keywords):
        return None
    if not _node_reference(keywords["value"], local=binder_local):
        return None
    if not _simple_node_reference(keywords["variable"]):
        return None
    # The expression is intentionally not tied to a mathematical producer;
    # it must nevertheless be a graph-node reference, as in gradient.value.node.
    if not _node_reference(keywords["expression"]):
        return None
    if not _constant(keywords["Z_independent_substitution"], True):
        return None
    variable_local = keywords["variable"].value.id  # type: ignore[union-attr]
    assert isinstance(variable_local, str)
    return node, variable_local


def _exact_integral_call(
    node: ast.AST,
    *,
    binder_local: str,
    substitution_local: str,
    variable_local: str,
) -> ast.Call | None:
    call = node if isinstance(node, ast.Call) else None
    if call is None:
        return None
    if not _attribute(call.func, "controls", "integral"):
        return None
    if len(call.args) < 5 or any(isinstance(arg, ast.Starred) for arg in call.args):
        return None
    if not _name(call.args[0], "g") or not _name(call.args[1], substitution_local):
        return None
    if not _g_nodes_name(call.args[2], binder_local):
        return None
    if not _self_anchor(call.args[3]) or not _name(call.args[4], variable_local):
        return None
    keywords = _keyword_values(call)
    if keywords is None:
        return None
    for flag in (
        "original_H_root_anchor",
        "anchor_is_exact_root_not_disk_center",
        "actual_anchored_gauge_integral",
    ):
        if flag not in keywords or not _constant(keywords[flag], True):
            return None
    return call


def _assignment_for_call(scope: ast.AST, call: ast.Call) -> str | None:
    for node in _scope_nodes(scope):
        if isinstance(node, ast.Assign) and node.value is call:
            return _assignment_name(node)
    return None


def _source_chain(
    source_text: str,
    *,
    binder_name: str,
    binder_local: str,
    substitution_local: str,
    integral_target: str | None,
) -> tuple[dict[str, Any] | None, list[str]]:
    errors: list[str] = []
    try:
        tree = ast.parse(source_text, filename="<authenticated-source>", mode="exec")
    except (SyntaxError, TypeError, ValueError) as exc:
        return None, ["source AST parse failed: " + type(exc).__name__ + ": " + str(exc)]

    candidates: list[dict[str, Any]] = []
    for scope_name, scope in _scopes(tree):
        if isinstance(scope, (ast.FunctionDef, ast.AsyncFunctionDef)):
            arguments = _function_arguments(scope)
            if binder_local in arguments or substitution_local in arguments:
                continue
        events = _binding_events(scope)
        # A local rebind of the two handles is ambiguous even if one of the
        # assignments happens to have the desired spelling.
        if len(events.get(binder_local, [])) != 1 or len(events.get(substitution_local, [])) != 1:
            continue
        nodes = sorted(_scope_nodes(scope), key=_position)
        binder_assignments = [
            item
            for item in nodes
            if _exact_symbol_assignment(
                item,
                binder_local=binder_local,
                binder_name=binder_name,
            )
        ]
        substitutions: list[tuple[ast.Assign, str]] = []
        for item in nodes:
            match = _exact_substitution_assignment(
                item,
                binder_local=binder_local,
                substitution_local=substitution_local,
            )
            if match is not None:
                substitutions.append(match)
        if len(binder_assignments) != 1 or len(substitutions) != 1:
            continue
        binder_assignment = binder_assignments[0]
        substitution_assignment, variable_local = substitutions[0]
        # The actual constructor may establish ``g`` in the same function
        # (for example ``self.g = g = self.pulse.g``).  Such a definition is
        # admissible only before the witness chain and only once; a later
        # rebind would make the graph owner ambiguous.  ``controls`` follows
        # the same conservative rule for source fragments that localize it.
        if len(events.get("g", [])) > 1 or len(events.get("controls", [])) > 1:
            continue
        chain_start = _position(binder_assignment)
        if any(_position(item) > chain_start for item in events.get("g", [])):
            continue
        if any(_position(item) > chain_start for item in events.get("controls", [])):
            continue
        calls = [
            item
            for item in nodes
            if isinstance(item, ast.Call)
            and _exact_integral_call(
                item,
                binder_local=binder_local,
                substitution_local=substitution_local,
                variable_local=variable_local,
            )
        ]
        if len(calls) != 1:
            continue
        integral_call = calls[0]
        if not (
            _position(binder_assignment)
            < _position(substitution_assignment)
            < _position(integral_call)
        ):
            continue
        target = _assignment_for_call(scope, integral_call)
        if integral_target is not None and target != integral_target:
            continue
        candidates.append(
            {
                "scope": scope_name,
                "binder_local": binder_local,
                "substitution_local": substitution_local,
                "variable_local": variable_local,
                "integral_target": target,
                "locations": {
                    "binder_definition": _loc(binder_assignment),
                    "substitution_definition": _loc(substitution_assignment),
                    "integral_constructor": _loc(integral_call),
                },
                "chain": [
                    "symbol -> substitution.value -> integral.integrand",
                    "binder.node -> g.nodes[binder.node]['name']",
                    "substitution.variable.node -> integral upper variable",
                ],
            }
        )
    if not candidates:
        errors.append("no single lexical scope contains the required source constructor chain")
        return None, errors
    if len(candidates) != 1:
        errors.append("source constructor chain is ambiguous across lexical scopes")
        return None, errors
    return candidates[0], errors


def _canonical_nodes(graph_nodes: Any) -> tuple[list[dict[str, Any]] | None, str | None, list[str]]:
    errors: list[str] = []
    if isinstance(graph_nodes, Mapping):
        # Accept only explicit container keys; no name-based node lookup is
        # performed.  This keeps callers from accidentally passing a wrapper.
        present = [key for key in ("graph_nodes", "exact_graph_nodes", "nodes") if key in graph_nodes]
        if len(present) != 1:
            return None, None, ["graph_nodes must be an explicit list, not an ambiguous mapping"]
        graph_nodes = graph_nodes[present[0]]
    if isinstance(graph_nodes, (str, bytes, bytearray)) or not isinstance(graph_nodes, Sequence):
        return None, None, ["graph_nodes must be an explicit sequence of serialized rows"]
    rows: list[dict[str, Any]] = []
    for index, row in enumerate(graph_nodes):
        if not isinstance(row, Mapping):
            errors.append(f"graph row {index} is not a mapping")
            continue
        rows.append({str(key): _copy_json(value) for key, value in row.items()})
    if errors:
        return None, None, errors
    try:
        serialized = json.dumps(
            rows,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=False,
            allow_nan=False,
        )
    except (TypeError, ValueError) as exc:
        return None, None, ["graph rows are not JSON-serializable: " + str(exc)]
    digest = hashlib.sha256(serialized.encode("utf-8")).hexdigest()
    # Round-trip through JSON so the returned support is exactly the structure
    # whose digest was calculated, not a caller-owned mutable object.
    try:
        canonical_rows = json.loads(serialized)
    except (TypeError, ValueError) as exc:  # pragma: no cover - defensive
        return None, None, ["canonical graph serialization failed: " + str(exc)]
    return canonical_rows, digest, []


def _node_at(rows: Sequence[Mapping[str, Any]], node_id: Any, label: str, errors: list[str]) -> Mapping[str, Any] | None:
    if not _is_int(node_id) or not 0 <= node_id < len(rows):
        errors.append(f"{label} must be an in-range integer node id")
        return None
    return rows[node_id]


def _validate_graph(
    rows: Sequence[Mapping[str, Any]],
    *,
    integral_node_id: int,
    binder_node_id: int,
) -> tuple[dict[str, Any] | None, list[str]]:
    errors: list[str] = []
    integral = _node_at(rows, integral_node_id, "integral_node_id", errors)
    binder = _node_at(rows, binder_node_id, "binder_node_id", errors)
    if integral is None or binder is None:
        return None, errors
    if binder.get("operation") != "bound_variable":
        errors.append("binder node is not a bound_variable row")
    binder_name = binder.get("name")
    if type(binder_name) is not str or not binder_name:
        errors.append("binder row name must be a non-empty string")
    if integral.get("operation") != "definite_integral":
        errors.append("integral node is not a definite_integral row")
    for field in ("lower", "upper"):
        _node_at(rows, integral.get(field), f"integral.{field}", errors)
    for flag in _REQUIRED_INTEGRAL_FLAGS:
        if not _is_true(integral.get(flag)):
            errors.append(f"integral row flag {flag!r} is not exactly True")
    integrand_id = integral.get("integrand")
    substitution = _node_at(rows, integrand_id, "integral.integrand", errors)
    if substitution is None:
        return None, errors
    if substitution.get("operation") != "function_substitution":
        errors.append("integral.integrand is not a function_substitution row")
    substitution_value = substitution.get("value")
    if substitution_value != binder_node_id or type(substitution_value) is not int:
        errors.append("function_substitution.value does not exactly reference binder_node_id")
    if _node_at(rows, substitution_value, "function_substitution.value", errors) is None:
        pass
    substitution_variable_id = substitution.get("variable")
    if integral.get("upper") != substitution_variable_id:
        errors.append("integral.upper is not exactly function_substitution.variable")
    variable_node = _node_at(rows, substitution.get("variable"), "function_substitution.variable", errors)
    if variable_node is None or variable_node.get("operation") != "bound_variable":
        errors.append("function_substitution.variable is not a bound_variable node")
    expression_node = _node_at(rows, substitution.get("expression"), "function_substitution.expression", errors)
    if expression_node is None:
        pass
    if not _is_true(substitution.get("Z_independent_substitution")):
        errors.append("function_substitution lacks exact Z_independent_substitution=True")
    lower_node = _node_at(rows, integral.get("lower"), "integral.lower", errors)
    if lower_node is not None and (
        lower_node.get("operation") != "current_original_source_parameter"
        or lower_node.get("name") != "true_core_anchor"
        or lower_node.get("Z_independent") is not True
        or lower_node.get("exact_fixed_source_parameter") is not True
    ):
        errors.append("integral.lower is not the explicit true_core_anchor source parameter")
    if type(integral.get("variable")) is not str:
        errors.append("integral.variable must be a string")
    elif binder_name and integral.get("variable") != binder_name:
        errors.append("integral.variable is not exactly binder_row['name']")
    if errors:
        return None, errors
    assert isinstance(binder_name, str)
    assert isinstance(integrand_id, int)
    return (
        {
            "integral_node_id": integral_node_id,
            "binder_node_id": binder_node_id,
            "substitution_node_id": integrand_id,
            "binder_name": binder_name,
            "integral_variable": integral.get("variable"),
            "substitution_variable_node_id": substitution.get("variable"),
            "substitution_expression_node_id": substitution.get("expression"),
            "flags": {
                key: integral.get(key)
                for key in _REQUIRED_INTEGRAL_FLAGS
            }
            | {"Z_independent_substitution": substitution.get("Z_independent_substitution")},
        },
        errors,
    )


def _base_result(*, passed: bool, status: str, errors: Sequence[str]) -> dict[str, Any]:
    return {
        "schema": RESULT_SCHEMA,
        "passed": bool(passed),
        "status": status,
        "errors": [str(item) for item in errors],
        "claims": {
            "candidate_witness_supported": bool(passed),
            "source_ast_chain_proved": False,
            "serialized_graph_edge_chain_proved": False,
            "global_unique_name_inference_used": False,
            "binder_mapping_installed": False,
            "live_runtime_owner_proved": False,
            "archived_owner_same_proved": False,
            "source_identity_verified": False,
            "runtime_control_flow_proved": False,
            "shadowing_control_flow_proved": False,
            "numerical_or_physical_claim": False,
        },
        "limitations": [
            "This result does not identify or install a live runtime owner.",
            "This result does not establish that an archived owner is the same object.",
            "The source digest is an integrity check only; source provenance remains caller-supplied.",
            "AST lexical presence does not prove branch execution, lambda/comprehension execution, or exception-path control flow.",
            "The lexical rebind checks do not authenticate runtime alias ownership.",
        ],
    }


def _validate(
    graph_nodes: Any,
    source_text: Any,
    *,
    integral_node_id: Any,
    binder_node_id: Any,
    source_sha256: str | None,
    expected_source_sha256: str | None,
    binder_local: str,
    substitution_local: str,
    integral_target: str | None,
) -> dict[str, Any]:
    rows, graph_digest, errors = _canonical_nodes(graph_nodes)
    if errors or rows is None or graph_digest is None:
        return _base_result(passed=False, status="rejected", errors=errors)
    if not _is_int(integral_node_id) or not _is_int(binder_node_id):
        return _base_result(
            passed=False,
            status="rejected",
            errors=["integral_node_id and binder_node_id must be explicit integer ids"],
        )
    if not isinstance(source_text, str) or not source_text.strip():
        return _base_result(passed=False, status="rejected", errors=["source_text is required raw text"])
    if not isinstance(binder_local, str) or not binder_local.isidentifier():
        return _base_result(passed=False, status="rejected", errors=["binder_local must be an identifier"])
    if not isinstance(substitution_local, str) or not substitution_local.isidentifier():
        return _base_result(passed=False, status="rejected", errors=["substitution_local must be an identifier"])
    if integral_target is not None and (
        not isinstance(integral_target, str) or not integral_target.isidentifier()
    ):
        return _base_result(passed=False, status="rejected", errors=["integral_target must be an identifier"])

    source_digest = hashlib.sha256(source_text.encode("utf-8")).hexdigest()
    expected_digests = [item for item in (source_sha256, expected_source_sha256) if item is not None]
    if any(not isinstance(item, str) or item.lower() != source_digest for item in expected_digests):
        return _base_result(passed=False, status="rejected", errors=["source SHA-256 does not match source_text"])
    if len(expected_digests) == 2 and expected_digests[0].lower() != expected_digests[1].lower():
        return _base_result(passed=False, status="rejected", errors=["source digest assertions disagree"])

    graph_evidence, graph_errors = _validate_graph(
        rows,
        integral_node_id=integral_node_id,
        binder_node_id=binder_node_id,
    )
    if graph_errors or graph_evidence is None:
        result = _base_result(passed=False, status="rejected", errors=graph_errors)
        result["graph_support"] = {
            "node_count": len(rows),
            "serialization_sha256": graph_digest,
        }
        return result
    source_evidence, source_errors = _source_chain(
        source_text,
        binder_name=graph_evidence["binder_name"],
        binder_local=binder_local,
        substitution_local=substitution_local,
        integral_target=integral_target,
    )
    if source_errors or source_evidence is None:
        result = _base_result(passed=False, status="rejected", errors=source_errors)
        result["graph_support"] = {
            "node_count": len(rows),
            "serialization_sha256": graph_digest,
            "validated_ids": {
                "integral_node_id": integral_node_id,
                "binder_node_id": binder_node_id,
                "substitution_node_id": graph_evidence["substitution_node_id"],
            },
        }
        result["source_support"] = {
            "sha256": source_digest,
            "digest_checked": bool(expected_digests),
        }
        return result

    result = _base_result(passed=True, status="candidate_witness", errors=[])
    result["candidate_witness"] = {
        "integral_node_id": integral_node_id,
        "binder_node_id": binder_node_id,
        "substitution_node_id": graph_evidence["substitution_node_id"],
        "variable_name": graph_evidence["binder_name"],
        "evidence_kind": "source_ast_plus_same_serialized_graph",
    }
    result["graph_support"] = {
        "node_count": len(rows),
        "serialization": "json(sort_keys=True,separators=(',',':'),ensure_ascii=False,allow_nan=False)",
        "serialization_sha256": graph_digest,
        "same_serialized_graph": True,
        "integral_row": _copy_json(rows[integral_node_id]),
        "substitution_row": _copy_json(rows[graph_evidence["substitution_node_id"]]),
        "binder_row": _copy_json(rows[binder_node_id]),
        "edge_chain": {
            "integral_integrand_is_substitution": True,
            "substitution_value_is_exact_binder": True,
            "integral_upper_is_exact_substitution_variable": True,
            "integral_lower_is_true_core_anchor_parameter": True,
            "integral_variable_is_exact_binder_name": True,
        },
    }
    result["source_support"] = {
        "sha256": source_digest,
        "digest_checked": bool(expected_digests),
        "lexical_scope": source_evidence["scope"],
        "locations": source_evidence["locations"],
        "definition_chain": source_evidence["chain"],
        "integral_target": source_evidence["integral_target"],
        "required_source_flags": [
            "Z_independent_substitution=True",
            "original_H_root_anchor=True",
            "anchor_is_exact_root_not_disk_center=True",
            "actual_anchored_gauge_integral=True",
        ],
    }
    result["claims"].update(
        {
            "source_ast_chain_proved": True,
            "serialized_graph_edge_chain_proved": True,
            "candidate_witness_supported": True,
        }
    )
    return result


def validate_core_gauge_binding_witness(
    graph_nodes: Any,
    source_text: Any,
    *,
    integral_node_id: Any,
    binder_node_id: Any,
    source_sha256: str | None = None,
    expected_source_sha256: str | None = None,
    binder_local: str = "v",
    substitution_local: str = "gv",
    integral_target: str | None = None,
) -> dict[str, Any]:
    """Return a fail-closed validation result without executing source code."""

    return _validate(
        graph_nodes,
        source_text,
        integral_node_id=integral_node_id,
        binder_node_id=binder_node_id,
        source_sha256=source_sha256,
        expected_source_sha256=expected_source_sha256,
        binder_local=binder_local,
        substitution_local=substitution_local,
        integral_target=integral_target,
    )


def build_core_gauge_binding_witness(
    graph_nodes: Any,
    source_text: Any,
    *,
    integral_node_id: Any,
    binder_node_id: Any,
    source_sha256: str | None = None,
    expected_source_sha256: str | None = None,
    binder_local: str = "v",
    substitution_local: str = "gv",
    integral_target: str | None = None,
) -> dict[str, Any]:
    """Build a candidate witness or raise :class:`GaugeBindingWitnessError`.

    The strict entry point is useful at installation boundaries: a rejected
    artifact cannot accidentally be mistaken for a witness.  Call
    :func:`validate_core_gauge_binding_witness` when a structured rejection is
    preferred.
    """

    result = validate_core_gauge_binding_witness(
        graph_nodes,
        source_text,
        integral_node_id=integral_node_id,
        binder_node_id=binder_node_id,
        source_sha256=source_sha256,
        expected_source_sha256=expected_source_sha256,
        binder_local=binder_local,
        substitution_local=substitution_local,
        integral_target=integral_target,
    )
    if not result.get("passed"):
        errors = result.get("errors")
        message = "; ".join(str(item) for item in errors) if isinstance(errors, list) else "witness rejected"
        raise GaugeBindingWitnessError(message, result)
    return result


# Short aliases make the boundary convenient without adding alternate logic.
build_witness = build_core_gauge_binding_witness
validate_witness = validate_core_gauge_binding_witness


__all__ = [
    "SCHEMA",
    "RESULT_SCHEMA",
    "GaugeBindingWitnessError",
    "WitnessValidationError",
    "validate_core_gauge_binding_witness",
    "build_core_gauge_binding_witness",
    "build_witness",
    "validate_witness",
]
