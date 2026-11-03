"""A source-preserving, executable program for saved raw g5 DAG nodes.

This module is deliberately independent of the producers that wrote an
archive.  ``compile_program`` consumes an already decoded graph (or its
``exact_graph_nodes`` list), and constructs a flat program by walking only
the explicitly supported arithmetic node schemas.  A source parameter,
bound variable, integral, substitution, analytic component, or unknown
record is an opaque boundary.  It is retained verbatim and can only be
evaluated through the caller's oracle.

The program is executable, but it is not a proof object.  Even an exact
backend and a successful oracle are caller-supplied components, so every
successful execution is reported as conditional and ``certified`` is always
false.  No function equality, physical m15, Kchi, or Navier--Stokes claim is
made here.  Integral/substitution nodes never execute in this version;
``bound_variable`` is only an explicit ``variable`` oracle boundary and does
not implement lexical binder capture.  ``graph_hash`` is caller-declared
metadata and is not authenticated by this module.
"""

from __future__ import annotations

from collections.abc import Mapping
import copy
import inspect
import math
import numbers
from typing import Any


__all__ = [
    "SCHEMA",
    "EXECUTION_SCHEMA",
    "CompileError",
    "GraphValidationError",
    "CycleError",
    "ExecutionError",
    "OpaqueOracleRequired",
    "compile_program",
    "execute",
]


SCHEMA = "n1_raw_g5_dag_program_v1"
EXECUTION_SCHEMA = "n1_raw_g5_dag_program_execution_v1"

_MISSING = object()
_MAX_NODES = 1_000_000

# These are the deterministic source schemas that can be represented by the
# flat program.  Everything else remains an opaque leaf.  The aliases cover
# the two records used by the saved source graphs (``operation``/``op`` and
# ``arguments``/``terms``/``factors``) without importing either producer.
_SUPPORTED_OPS = frozenset(
    {
        "exact_rational",
        "sum",
        "product",
        "negative",
        "positive_quotient",
        "analytic_unary",
    }
)

_OPAQUE_OPS = frozenset(
    {
        "bound_variable",
        "source_parameter",
        "current_original_source_parameter",
        "original_source_parameter",
        "current_original_physical_source_parameter",
        "component",
        "source_component",
        "current_original_source_component",
        "current_original_component",
        "current_original_unique_regular_n1_component",
        "analytic_component",
        "analytic_core_component",
        "current_original_unique_analytic_core_component",
        "current_original_analytic_core_component",
        "definite_integral",
        "function_substitution",
    }
)

# The q4 runtime and the saved n1 producers use these analytic names.  An
# unrecognised function name is an opaque boundary; it must not silently turn
# into an arbitrary backend call.
_SUPPORTED_ANALYTIC_FUNCTIONS = frozenset(
    {"abs", "exp", "log", "log1p", "sqrt", "positive_sqrt"}
)


class CompileError(ValueError):
    """Base class for fail-closed graph compilation errors."""


class GraphValidationError(CompileError):
    """The graph, root list, or a supported node schema is malformed."""


class CycleError(GraphValidationError):
    """A cycle was found in the traversed supported DAG portion."""


class ExecutionError(RuntimeError):
    """A backend or executable program operation failed closed."""

    def __init__(self, message: str, *, code: str = "execution_error", node_id: Any = None):
        self.code = code
        self.node_id = node_id
        super().__init__(message)


class OpaqueOracleRequired(ExecutionError):
    """An opaque node was reached without an explicit caller oracle."""

    def __init__(self, node_id: int, operation: str = "unknown"):
        super().__init__(
            f"opaque {operation} node {node_id} requires an explicit oracle",
            code="oracle_required",
            node_id=node_id,
        )


def _copy(value: Any, *, label: str) -> Any:
    """Copy a source record without allowing compiler-side mutation."""

    try:
        return copy.deepcopy(value)
    except Exception as exc:  # pragma: no cover - defensive boundary
        raise GraphValidationError(f"{label} is not copyable") from exc


def _is_node_id(value: Any) -> bool:
    """Accept only a real integer; bool is intentionally rejected."""

    return type(value) is int


def _node_id(value: Any, *, label: str, count: int) -> int:
    if not _is_node_id(value):
        raise GraphValidationError(f"{label}: node id must be an int, got {value!r}")
    if value < 0 or value >= count:
        raise GraphValidationError(
            f"{label}: node id {value!r} is outside exact_graph_nodes[0:{count}]"
        )
    return value


def _graph_nodes(graph: Any) -> list[Any]:
    """Extract an already decoded indexed graph without importing a producer."""

    candidate: Any = graph
    if isinstance(graph, Mapping):
        candidate = _MISSING
        for key in ("exact_graph_nodes", "graph_nodes", "nodes"):
            if key in graph:
                candidate = graph[key]
                break
        if candidate is _MISSING:
            raise GraphValidationError(
                "graph must contain exact_graph_nodes (or graph_nodes/nodes)"
            )
    if not isinstance(candidate, (list, tuple)):
        raise GraphValidationError("graph nodes must be an indexed list or tuple")
    if len(candidate) > _MAX_NODES:
        raise GraphValidationError("graph node limit exceeded")
    if not candidate:
        raise GraphValidationError("graph nodes must be non-empty")
    return list(candidate)


def _roots(roots: Any, *, count: int) -> list[int]:
    # A root list is intentionally required.  This keeps the multiple-root
    # boundary explicit and avoids treating a single integer as a hidden list.
    if isinstance(roots, (str, bytes, bytearray, Mapping)) or not isinstance(
        roots, (list, tuple)
    ):
        raise GraphValidationError("roots must be an explicit list or tuple of node ids")
    if not roots:
        raise GraphValidationError("at least one root node id is required")
    return [
        _node_id(root, label=f"roots[{index}]", count=count)
        for index, root in enumerate(roots)
    ]


def _operation(record: Any, *, node_id: int) -> str | None:
    if not isinstance(record, Mapping):
        return None
    operation = record.get("operation", _MISSING)
    shorthand = record.get("op", _MISSING)
    if operation is not _MISSING and shorthand is not _MISSING and operation != shorthand:
        raise GraphValidationError(
            f"node {node_id}: operation and op disagree ({operation!r}, {shorthand!r})"
        )
    value = operation if operation is not _MISSING else shorthand
    if value is _MISSING:
        return None
    if not isinstance(value, str) or not value:
        # A non-string operation is unknown, hence an opaque boundary.  It is
        # preserved in the original record rather than guessed.
        return None
    return value


def _field(record: Mapping[str, Any], names: tuple[str, ...], *, node_id: int) -> Any:
    present = [(name, record[name]) for name in names if name in record]
    if not present:
        raise GraphValidationError(f"node {node_id}: missing {'/'.join(names)}")
    # If both aliases are present, only accept the same value.  This avoids
    # compiling a source row whose bracket/operand type is ambiguous.
    first_name, first = present[0]
    for name, value in present[1:]:
        if value != first:
            raise GraphValidationError(
                f"node {node_id}: conflicting fields {first_name} and {name}"
            )
    return first


def _children_for_supported(
    record: Mapping[str, Any], op: str, *, node_id: int, count: int
) -> tuple[list[int], dict[str, Any]]:
    """Validate a supported record and return child IDs plus canonical fields."""

    if op == "exact_rational":
        has_numerator = "numerator" in record
        has_value = "value" in record
        if has_numerator and has_value and record["numerator"] != record["value"]:
            raise GraphValidationError(
                f"node {node_id}: numerator and value disagree for exact rational"
            )
        numerator = (
            record["numerator"]
            if has_numerator
            else record.get("value", _MISSING)
        )
        denominator = record.get("denominator", 1)
        if type(numerator) is not int or type(denominator) is not int:
            raise GraphValidationError(
                f"node {node_id}: exact rational needs integer numerator/denominator"
            )
        if denominator == 0:
            raise GraphValidationError(f"node {node_id}: exact rational denominator is zero")
        return [], {"numerator": numerator, "denominator": denominator}

    if op in ("sum", "product"):
        names = ("arguments", "terms") if op == "sum" else ("arguments", "factors")
        raw_children = _field(record, names, node_id=node_id)
        if not isinstance(raw_children, (list, tuple)) or not raw_children:
            raise GraphValidationError(f"node {node_id}: {op} needs a non-empty list/tuple")
        children = [
            _node_id(item, label=f"node {node_id} {op}[{index}]", count=count)
            for index, item in enumerate(raw_children)
        ]
        return children, {"children": children}

    if op == "negative":
        child = _node_id(
            _field(record, ("argument",), node_id=node_id),
            label=f"node {node_id}.argument",
            count=count,
        )
        return [child], {"argument": child}

    if op == "positive_quotient":
        numerator = _node_id(
            _field(record, ("numerator",), node_id=node_id),
            label=f"node {node_id}.numerator",
            count=count,
        )
        denominator = _node_id(
            _field(record, ("denominator",), node_id=node_id),
            label=f"node {node_id}.denominator",
            count=count,
        )
        return [numerator, denominator], {
            "numerator": numerator,
            "denominator": denominator,
        }

    if op == "analytic_unary":
        function = _field(record, ("function", "name"), node_id=node_id)
        if not isinstance(function, str) or not function.strip():
            raise GraphValidationError(
                f"node {node_id}: analytic function must be a non-empty string"
            )
        function = function.strip().lower()
        if function not in _SUPPORTED_ANALYTIC_FUNCTIONS:
            # Unknown analytic semantics are an opaque boundary.  Returning a
            # special marker lets the caller retain the row without walking
            # its argument.
            return [], {"opaque_reason": "unknown_analytic_function", "function": function}
        child = _node_id(
            _field(record, ("argument",), node_id=node_id),
            label=f"node {node_id}.argument",
            count=count,
        )
        return [child], {"argument": child, "function": function}

    raise GraphValidationError(f"node {node_id}: unsupported compiler operation {op!r}")


def _opaque_kind(operation: str | None) -> str:
    if operation == "bound_variable":
        return "variable"
    if operation in {
        "source_parameter",
        "current_original_source_parameter",
        "original_source_parameter",
        "current_original_physical_source_parameter",
    }:
        return "parameter"
    if operation in {
        "component",
        "source_component",
        "current_original_source_component",
        "current_original_component",
        "current_original_unique_regular_n1_component",
        "analytic_component",
    }:
        return "regular_component"
    if operation in {
        "analytic_core_component",
        "current_original_unique_analytic_core_component",
        "current_original_analytic_core_component",
    }:
        return "component"
    if operation == "definite_integral":
        return "integral"
    if operation == "function_substitution":
        return "substitution"
    return "unknown"


def _classify(record: Any, *, node_id: int, count: int) -> dict[str, Any]:
    """Classify one record without traversing opaque children."""

    op = _operation(record, node_id=node_id)
    if op not in _SUPPORTED_OPS:
        reason = "unknown_operation" if op is not None else "unknown_record"
        if op in _OPAQUE_OPS:
            reason = f"opaque_{op}"
        return {
            "supported": False,
            "op": "opaque",
            "opaque_operation": op,
            "opaque_kind": _opaque_kind(op),
            "opaque_reason": reason,
            "children": [],
        }

    if not isinstance(record, Mapping):  # currently unreachable, kept explicit
        return {
            "supported": False,
            "op": "opaque",
            "opaque_operation": op,
            "opaque_kind": "unknown",
            "opaque_reason": "non_mapping_record",
            "children": [],
        }

    children, fields = _children_for_supported(record, op, node_id=node_id, count=count)
    if fields.get("opaque_reason"):
        return {
            "supported": False,
            "op": "opaque",
            "opaque_operation": op,
            "opaque_kind": "unknown",
            "opaque_reason": fields["opaque_reason"],
            "function": fields.get("function"),
            "children": [],
        }
    return {"supported": True, "op": op, "children": children, **fields}


def _claims() -> dict[str, bool]:
    return {
        "function_equality": False,
        "physical_m15": False,
        "Kchi": False,
        "Navier_Stokes": False,
        "certified": False,
    }


def compile_program(graph: Any, roots: Any, graph_hash: Any) -> dict[str, Any]:
    """Compile explicit root IDs into a flat, executable source DAG program.

    The traversal is iterative.  Supported operands are indexed and checked;
    opaque records are emitted as leaves, so an integral binder or an unknown
    nested payload is never interpreted or traversed.  The returned mapping
    keeps each original row alongside its canonical instruction and retains
    one instruction per shared source node.
    """

    nodes = _graph_nodes(graph)
    root_ids = _roots(roots, count=len(nodes))
    if graph_hash is None or type(graph_hash) is bool:
        raise GraphValidationError("graph_hash must be a non-null non-boolean value")

    # A three-state iterative DFS gives dependency-first topological order.
    # State 1 is on the active stack, state 2 is emitted.  Opaque nodes are
    # emitted immediately and therefore form deliberate traversal boundaries.
    state: dict[int, int] = {}
    classified: dict[int, dict[str, Any]] = {}
    topological: list[int] = []
    for root in root_ids:
        if state.get(root) == 2:
            continue
        stack: list[tuple[int, bool]] = [(root, False)]
        while stack:
            node_id, leaving = stack.pop()
            current = state.get(node_id, 0)
            if leaving:
                if current != 1:
                    # A repeated dependency was already emitted; no second
                    # instruction is allowed for a shared source node.
                    if current == 2:
                        continue
                    raise GraphValidationError(f"node {node_id}: invalid traversal state")
                state[node_id] = 2
                topological.append(node_id)
                continue
            if current == 2:
                continue
            if current == 1:
                raise CycleError(f"cycle detected through source node {node_id}")

            state[node_id] = 1
            info = classified.get(node_id)
            if info is None:
                info = _classify(nodes[node_id], node_id=node_id, count=len(nodes))
                classified[node_id] = info
            children = info["children"]
            if not info["supported"]:
                # Opaque nodes intentionally do not inspect any payload child.
                state[node_id] = 2
                topological.append(node_id)
                continue

            stack.append((node_id, True))
            for child in reversed(children):
                if state.get(child, 0) == 1:
                    raise CycleError(
                        f"cycle detected: node {node_id} references active node {child}"
                    )
                if state.get(child, 0) != 2:
                    stack.append((child, False))

    instructions: list[dict[str, Any]] = []
    by_id: dict[int, dict[str, Any]] = {}
    opaque_ids: list[int] = []
    opaque_reachable: set[int] = set()
    records: dict[int, Any] = {}

    for node_id in topological:
        info = classified[node_id]
        original = _copy(nodes[node_id], label=f"node {node_id} original_record")
        records[node_id] = _copy(original, label=f"node {node_id} record copy")
        instruction: dict[str, Any] = {
            "node_id": node_id,
            "record": original,
            "original_record": _copy(original, label=f"node {node_id} original_record copy"),
            "status": "supported" if info["supported"] else "opaque",
            "op": info["op"],
        }
        if info["supported"]:
            instruction["children"] = list(info["children"])
            for field in ("numerator", "denominator", "argument", "function"):
                if field in info:
                    instruction[field] = info[field]
        else:
            instruction.update(
                {
                    "opaque_operation": info.get("opaque_operation"),
                    "opaque_kind": info.get("opaque_kind", "unknown"),
                    "opaque_reason": info.get("opaque_reason", "opaque_boundary"),
                }
            )
            if info.get("function") is not None:
                instruction["function"] = info["function"]
            opaque_ids.append(node_id)
            opaque_reachable.add(node_id)
        instructions.append(instruction)
        by_id[node_id] = instruction

    # Mark supported nodes whose result depends on an opaque leaf.  This is a
    # status annotation only; it never changes the executable instruction.
    for instruction in instructions:
        if instruction["status"] != "supported":
            continue
        if any(child in opaque_reachable for child in instruction["children"]):
            opaque_reachable.add(instruction["node_id"])

    root_entries: list[dict[str, Any]] = []
    for root in root_ids:
        root_instruction = by_id[root]
        if root in opaque_reachable:
            root_status = (
                "opaque"
                if root_instruction["status"] == "opaque"
                else "conditional_opaque"
            )
        else:
            root_status = "supported"
        root_entries.append(
            {
                "root_node_id": root,
                "program_node_id": root,
                "status": root_status,
                "original_record": _copy(
                    nodes[root], label=f"root {root} original_record"
                ),
            }
        )

    # The returned mapping is intentionally redundant at the top level: both
    # an instruction list (stable topo order) and an ID index are useful to
    # callers, while records gives a direct source-preserving audit view.
    return {
        "schema": SCHEMA,
        "status": "compiled",
        "compile_status": "compiled",
        "graph_hash": _copy(graph_hash, label="graph_hash"),
        "node_count": len(nodes),
        "roots": list(root_ids),
        "root_entries": root_entries,
        "root_statuses": [entry["status"] for entry in root_entries],
        "topological_order": list(topological),
        "instructions": instructions,
        "nodes": by_id,
        "records": records,
        "missing_oracle_nodes": list(opaque_ids),
        "opaque_nodes": list(opaque_ids),
        "conditional": True,
        "certified": False,
        "claims": _claims(),
        "provenance": {
            "source_records_preserved": True,
            "source_identity_verified": False,
            "backend_verified": False,
            "oracle_verified": False,
            "graph_hash_caller_declared": True,
            "graph_hash_authenticated": False,
            "lexical_binder_supported": False,
            "integral_runtime_supported": False,
            "substitution_runtime_supported": False,
            "claim": "conditional executable arithmetic only",
        },
    }


def _member(owner: Any, names: tuple[str, ...]) -> Any:
    if isinstance(owner, Mapping):
        for name in names:
            if name in owner:
                return owner[name]
        return _MISSING
    for name in names:
        try:
            return getattr(owner, name)
        except Exception:
            continue
    return _MISSING


def _backend_call(
    backend: Any,
    names: tuple[str, ...],
    args: tuple[Any, ...],
    *,
    node_id: int,
    code: str,
) -> Any:
    method = _member(backend, names)
    if not callable(method):
        raise ExecutionError(
            f"backend method {'/'.join(names)} is required for node {node_id}",
            code=code,
            node_id=node_id,
        )
    try:
        value = method(*args)
    except ExecutionError:
        raise
    except Exception as exc:
        raise ExecutionError(
            f"backend operation failed at node {node_id}: {type(exc).__name__}",
            code="backend_error",
            node_id=node_id,
        ) from exc
    if value is None or type(value) is bool:
        raise ExecutionError(
            f"backend returned an invalid value at node {node_id}",
            code="backend_invalid_value",
            node_id=node_id,
        )
    _validate_scalar(value, backend, node_id=node_id, source="backend")
    return value


def _validate_scalar(value: Any, backend: Any, *, node_id: int, source: str) -> None:
    """Reject metadata-shaped or non-finite values at a scalar boundary.

    The caller still supplies the semantic meaning and authenticity of the
    scalar.  This check only prevents a record/cap container, a boolean, or a
    non-finite numeric value from being smuggled into arithmetic.  If a
    backend advertises a scalar validator, it is authoritative for its own
    scalar type and must return exactly ``True``.
    """

    if value is None or type(value) is bool:
        raise ExecutionError(
            f"{source} returned a non-scalar value at node {node_id}",
            code=f"{source}_invalid_scalar",
            node_id=node_id,
        )
    if isinstance(value, (Mapping, list, tuple, set, frozenset, str, bytes, bytearray)):
        raise ExecutionError(
            f"{source} returned a container/non-numeric value at node {node_id}",
            code=f"{source}_invalid_scalar",
            node_id=node_id,
        )
    if isinstance(value, float) and not math.isfinite(value):
        raise ExecutionError(
            f"{source} returned a non-finite float at node {node_id}",
            code=f"{source}_non_finite",
            node_id=node_id,
        )
    if isinstance(value, complex):
        raise ExecutionError(
            f"{source} returned an unsupported complex value at node {node_id}",
            code=f"{source}_invalid_scalar",
            node_id=node_id,
        )
    finite_member = _member(value, ("is_finite",))
    if finite_member is not _MISSING:
        try:
            finite = finite_member() if callable(finite_member) else finite_member
        except Exception as exc:
            raise ExecutionError(
                f"{source} scalar finiteness check failed at node {node_id}",
                code=f"{source}_non_finite",
                node_id=node_id,
            ) from exc
        if type(finite) is not bool or not finite:
            raise ExecutionError(
                f"{source} returned a non-finite scalar at node {node_id}",
                code=f"{source}_non_finite",
                node_id=node_id,
            )

    validator = _member(
        backend,
        ("validate_scalar", "scalar_validator", "is_scalar", "is_numeric_scalar"),
    )
    if not isinstance(value, numbers.Real):
        # Fraction is a numbers.Real and remains accepted without an adapter.
        # Decimal, interval, and other custom scalar classes must be admitted
        # explicitly by the backend's validate_scalar contract.
        explicit_validator = _member(backend, ("validate_scalar",))
        if not callable(explicit_validator):
            raise ExecutionError(
                f"{source} returned a non-Real scalar without validate_scalar at node {node_id}",
                code=f"{source}_scalar_validator_required",
                node_id=node_id,
            )
        validator = explicit_validator
    if validator is not _MISSING:
        if not callable(validator):
            raise ExecutionError(
                "backend scalar validator is not callable",
                code="backend_scalar_validator_error",
                node_id=node_id,
            )
        try:
            valid = validator(value)
        except Exception as exc:
            raise ExecutionError(
                f"backend scalar validator failed at node {node_id}",
                code="backend_scalar_validator_error",
                node_id=node_id,
            ) from exc
        if type(valid) is not bool or not valid:
            raise ExecutionError(
                f"backend scalar validator rejected value at node {node_id}",
                code="backend_scalar_rejected",
                node_id=node_id,
            )


def _positive(backend: Any, value: Any, *, node_id: int) -> bool:
    method = _member(
        backend,
        ("is_positive", "strictly_positive", "is_strictly_positive", "positive_denominator"),
    )
    if not callable(method):
        raise ExecutionError(
            f"backend strict-positive validation is required at node {node_id}",
            code="backend_positive_required",
            node_id=node_id,
        )
    try:
        result = method(value)
    except Exception as exc:
        raise ExecutionError(
            f"backend positivity check failed at node {node_id}",
            code="backend_positive_error",
            node_id=node_id,
        ) from exc
    if type(result) is not bool:
        raise ExecutionError(
            f"backend positivity check was not bool at node {node_id}",
            code="backend_positive_not_bool",
            node_id=node_id,
        )
    return result


def _unary(backend: Any, function: str, value: Any, *, node_id: int) -> Any:
    generic = _member(backend, ("analytic_unary", "unary", "apply_unary"))
    if callable(generic):
        return _backend_call(
            backend,
            ("analytic_unary", "unary", "apply_unary"),
            (function, value),
            node_id=node_id,
            code="backend_unary_required",
        )
    names = {
        "abs": ("abs", "absolute"),
        "exp": ("exp",),
        "log": ("log",),
        "log1p": ("log1p",),
        "sqrt": ("sqrt",),
        "positive_sqrt": ("positive_sqrt", "sqrt"),
    }.get(function, ())
    if not names:
        raise ExecutionError(
            f"unsupported analytic function {function!r} at node {node_id}",
            code="backend_unary_required",
            node_id=node_id,
        )
    return _backend_call(backend, names, (value,), node_id=node_id, code="backend_unary_required")


def _call_oracle(
    method: Any, node_id: int, record: Any, operation: str, backend: Any
) -> Any:
    """Call a category-aware oracle without guessing scalar values.

    Preferred callable signature is ``oracle(node_id, original_record,
    operation)``.  A two-argument callable remains accepted for compatibility
    with the public boundary in the task statement, but the record still
    retains its operation and no binder/endpoint payload is interpreted.
    Category-named object methods (``variable``, ``parameter``,
    ``regular_component``, ``analytic_core_component``, ``integral``,
    ``substitution``, ``unknown``) are also accepted.  There is no lexical
    binder scope or implicit ``resolve`` callback.
    """

    if not callable(method):
        raise OpaqueOracleRequired(node_id, operation)
    try:
        try:
            signature = inspect.signature(method)
        except (TypeError, ValueError):
            signature = None
        if signature is not None:
            parameters = list(signature.parameters.values())
            accepts_varargs = any(
                item.kind is inspect.Parameter.VAR_POSITIONAL for item in parameters
            )
            positional = [
                item
                for item in parameters
                if item.kind
                in (inspect.Parameter.POSITIONAL_ONLY, inspect.Parameter.POSITIONAL_OR_KEYWORD)
            ]
            if accepts_varargs or len(positional) >= 3:
                value = method(node_id, record, operation)
            else:
                value = method(node_id, record)
        else:
            value = method(node_id, record, operation)
    except OpaqueOracleRequired:
        raise
    except Exception as exc:
        raise ExecutionError(
            f"oracle failed for opaque {operation} node {node_id}: {type(exc).__name__}",
            code="oracle_error",
            node_id=node_id,
        ) from exc
    # Returning the index itself is a common accidental placeholder and is
    # never accepted as an oracle result.  A caller that genuinely needs that
    # integer must return it through a typed/exact scalar backend value.
    if type(value) is int and value == node_id:
        raise ExecutionError(
            f"oracle returned node id {node_id} as its own value",
            code="oracle_node_id_value",
            node_id=node_id,
        )
    _validate_scalar(value, backend, node_id=node_id, source="oracle")
    return value


def _oracle_method(oracle: Any, operation: str, node_id: int) -> Any:
    if callable(oracle):
        return oracle
    if isinstance(oracle, Mapping):
        method = oracle.get(operation, _MISSING)
        if callable(method):
            return method
    method = _member(oracle, (operation,))
    if callable(method):
        return method
    raise OpaqueOracleRequired(node_id, operation)


def _validate_program(program: Any) -> tuple[list[dict[str, Any]], list[int], Any]:
    if not isinstance(program, Mapping) or program.get("schema") != SCHEMA:
        raise ExecutionError(
            "a compiled n1_raw_g5_dag_program_v1 is required", code="invalid_program"
        )
    instructions = program.get("instructions")
    roots = program.get("roots")
    if not isinstance(instructions, list) or not isinstance(roots, list) or not roots:
        raise ExecutionError(
            "compiled program has invalid instructions or roots", code="invalid_program"
        )
    by_id: dict[int, dict[str, Any]] = {}
    for index, instruction in enumerate(instructions):
        if not isinstance(instruction, Mapping):
            raise ExecutionError(
                f"instruction {index} is not a mapping", code="invalid_program"
            )
        node_id = instruction.get("node_id")
        if not _is_node_id(node_id) or node_id in by_id:
            raise ExecutionError(
                f"instruction {index} has an invalid node id", code="invalid_program"
            )
        by_id[node_id] = dict(instruction)
    for index, root in enumerate(roots):
        if not _is_node_id(root) or root not in by_id:
            raise ExecutionError(f"root {index} is absent from program", code="invalid_program")
    # Preserve the compiler's order and reject an instruction that references
    # a node emitted later; this keeps execution genuinely topological.
    seen: set[int] = set()
    for instruction in instructions:
        node_id = instruction["node_id"]
        if instruction.get("status") == "supported":
            children = instruction.get("children")
            if not isinstance(children, list) or not children:
                # exact_rational has no children; all other supported ops do.
                if instruction.get("op") != "exact_rational":
                    raise ExecutionError(
                        f"node {node_id} has invalid children", code="invalid_program"
                    )
            for child in children or []:
                if not _is_node_id(child) or child not in by_id or child not in seen:
                    raise ExecutionError(
                        f"node {node_id} has a non-topological child {child!r}",
                        code="invalid_program",
                    )
        elif instruction.get("status") != "opaque":
            raise ExecutionError(
                f"node {node_id} has invalid status", code="invalid_program"
            )
        seen.add(node_id)
    return [dict(item) for item in instructions], list(roots), program.get("graph_hash")


def execute(program: Any, backend: Any, oracle: Any = None) -> dict[str, Any]:
    """Execute supported arithmetic, resolving opaque leaves only by oracle.

    Exact rationals are passed as integer numerator/denominator pairs to an
    explicit backend method.  No floating conversion or implicit zero is
    available.  Successful results are conditional and uncertified by design.
    """

    instructions, root_ids, graph_hash = _validate_program(program)
    values: dict[int, Any] = {}
    opaque_used: list[int] = []

    for instruction in instructions:
        node_id = instruction["node_id"]
        if instruction["status"] == "opaque":
            operation = instruction.get("opaque_kind", "unknown")
            if operation in {"integral", "substitution"}:
                raise ExecutionError(
                    f"opaque {operation} node {node_id} has no scalar runtime",
                    code="opaque_runtime_not_supported",
                    node_id=node_id,
                )
            if oracle is None:
                raise OpaqueOracleRequired(node_id, operation)
            method = _oracle_method(oracle, operation, node_id)
            value = _call_oracle(
                method,
                node_id,
                instruction.get("original_record", instruction.get("record")),
                operation,
                backend,
            )
            values[node_id] = value
            opaque_used.append(node_id)
            continue

        op = instruction.get("op")
        if op == "exact_rational":
            values[node_id] = _backend_call(
                backend,
                ("rational", "exact_rational", "make_rational", "from_rational"),
                (instruction.get("numerator"), instruction.get("denominator", 1)),
                node_id=node_id,
                code="backend_rational_required",
            )
        elif op == "sum":
            children = instruction["children"]
            value = values[children[0]]
            for child in children[1:]:
                value = _backend_call(
                    backend,
                    ("add", "plus"),
                    (value, values[child]),
                    node_id=node_id,
                    code="backend_add_required",
                )
            values[node_id] = value
        elif op == "product":
            children = instruction["children"]
            value = values[children[0]]
            for child in children[1:]:
                value = _backend_call(
                    backend,
                    ("multiply", "mul", "times"),
                    (value, values[child]),
                    node_id=node_id,
                    code="backend_multiply_required",
                )
            values[node_id] = value
        elif op == "negative":
            values[node_id] = _backend_call(
                backend,
                ("negative", "negate", "neg"),
                (values[instruction["argument"]],),
                node_id=node_id,
                code="backend_negative_required",
            )
        elif op == "positive_quotient":
            numerator = values[instruction["numerator"]]
            denominator = values[instruction["denominator"]]
            if not _positive(backend, denominator, node_id=node_id):
                raise ExecutionError(
                    f"quotient denominator is not strictly positive at node {node_id}",
                    code="quotient_denominator_not_positive",
                    node_id=node_id,
                )
            positive_div = _member(
                backend, ("positive_div", "positive_divide", "positive_quotient")
            )
            if callable(positive_div):
                values[node_id] = _backend_call(
                    backend,
                    ("positive_div", "positive_divide", "positive_quotient"),
                    (numerator, denominator),
                    node_id=node_id,
                    code="backend_positive_divide_required",
                )
            else:
                values[node_id] = _backend_call(
                    backend,
                    ("divide", "div", "quotient"),
                    (numerator, denominator),
                    node_id=node_id,
                    code="backend_divide_required",
                )
        elif op == "analytic_unary":
            values[node_id] = _unary(
                backend,
                instruction["function"],
                values[instruction["argument"]],
                node_id=node_id,
            )
        else:  # pragma: no cover - compiler never emits another supported op
            raise ExecutionError(
                f"unknown executable operation {op!r}",
                code="unknown_operation",
                node_id=node_id,
            )

        if values[node_id] is None or type(values[node_id]) is bool:
            raise ExecutionError(
                f"backend returned an invalid value at node {node_id}",
                code="backend_invalid_value",
                node_id=node_id,
            )

    root_values = [values[root] for root in root_ids]
    root_entries = []
    compiled_entries = program.get("root_entries", [])
    for index, (root, value) in enumerate(zip(root_ids, root_values)):
        compile_status = (
            compiled_entries[index].get("status")
            if isinstance(compiled_entries, list)
            and index < len(compiled_entries)
            and isinstance(compiled_entries[index], Mapping)
            else "supported"
        )
        root_entries.append(
            {
                "root_node_id": root,
                "status": "conditional_opaque" if compile_status != "supported" else "conditional",
                "value": value,
                "actual_value_selected": False,
                "certified": False,
            }
        )

    return {
        "schema": EXECUTION_SCHEMA,
        "status": "conditional",
        "ok": True,
        "conditional": True,
        "certified": False,
        "graph_hash": _copy(graph_hash, label="program graph_hash"),
        "roots": list(root_ids),
        "root_results": root_entries,
        "results": root_entries,
        "root_values": root_values,
        "value": root_values[0] if len(root_values) == 1 else None,
        "values": dict(values),
        "opaque_nodes_used": opaque_used,
        "missing_oracle_nodes": [],
        "claims": _claims(),
        "provenance": {
            "source_identity_verified": False,
            "backend_verified": False,
            "oracle_verified": False,
            "actual_value_selected": False,
            "graph_hash_caller_declared": True,
            "graph_hash_authenticated": False,
            "lexical_binder_supported": False,
            "integral_runtime_supported": False,
            "substitution_runtime_supported": False,
            "claim": "conditional executable arithmetic only",
        },
    }

