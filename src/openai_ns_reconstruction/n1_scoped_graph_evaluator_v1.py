"""A small lexical evaluator for the real raw n1 graph row schema.

This module is deliberately separate from ``n1_raw_g5_dag_program_v1``.
The latter preserves integral and substitution rows as opaque leaves.  This
evaluator consumes the already authenticated ``exact_graph_nodes`` list and
adds only the missing scoped semantics; it does not import, rebuild, or
renumber a source graph.

The raw schema has two different binder encodings:

* ``function_substitution`` stores ``variable`` as a bound-variable node id.
* ``definite_integral`` stores ``variable`` as a string.  A caller must pass
  ``integral_binders[integral_node_id]`` explicitly.  The sidecar is checked
  against the bound-variable row's name; no name-based guessing is allowed.

All values are conditional backend results.  The evaluator makes no source
identity, mathematical, physical, Kchi, Navier--Stokes, or certification
claim.  Unknown rows remain opaque and require a category-specific oracle.
"""

from __future__ import annotations

from collections.abc import Mapping
import copy
from dataclasses import dataclass
import math
import numbers
from typing import Any


SCHEMA = "n1_scoped_graph_evaluator_v1"
RESULT_SCHEMA = "n1_scoped_graph_evaluator_result_v1"

_MISSING = object()
_MAX_NODES = 1_000_000
_MAX_DEPTH = 100_000

_SUPPORTED_ANALYTIC_FUNCTIONS = frozenset(
    {"abs", "exp", "log", "log1p", "sqrt", "positive_sqrt", "exprel"}
)

_OPAQUE_CATEGORIES = {
    "source_parameter": "parameter",
    "current_original_source_parameter": "parameter",
    "original_source_parameter": "parameter",
    "current_original_physical_source_parameter": "parameter",
    "component": "regular_component",
    "source_component": "regular_component",
    "current_original_source_component": "regular_component",
    "current_original_component": "regular_component",
    "current_original_unique_regular_n1_component": "regular_component",
    "analytic_component": "regular_component",
    "analytic_core_component": "component",
    "current_original_unique_analytic_core_component": "component",
    "current_original_analytic_core_component": "component",
}

_CLAIMS_FALSE = {
    "function_equality": False,
    "physical_m15": False,
    "Kchi": False,
    "Navier_Stokes": False,
    "certified": False,
}


class ScopedGraphError(ValueError):
    """Base class for fail-closed graph and scope errors."""


class GraphValidationError(ScopedGraphError):
    """The graph row, root, reference, or sidecar is malformed."""


class ScopeResolutionError(ScopedGraphError):
    """A lexical binder or free variable could not be resolved."""


class CycleError(ScopedGraphError):
    """A graph cycle or an evaluator recursion limit was reached."""


class ExecutionError(RuntimeError):
    """A backend or opaque oracle failed closed."""

    def __init__(self, message: str, *, code: str = "execution_error", node_id: Any = None):
        self.code = code
        self.node_id = node_id
        super().__init__(message)


class OpaqueOracleRequired(ExecutionError):
    """An opaque category was reached without its explicit oracle."""

    def __init__(self, node_id: int, category: str):
        super().__init__(
            f"opaque {category} node {node_id} requires an explicit oracle",
            code="oracle_required",
            node_id=node_id,
        )


@dataclass(frozen=True)
class ScopeEnv:
    """Persistent lexical frames keyed by bound-variable node id.

    Frames are prepended and never mutated.  Node ids, rather than variable
    names, are the identity key, so same-name binders remain distinct.
    """

    frames: tuple[tuple[int, Any], ...] = ()

    def bind(self, node_id: int, value: Any) -> "ScopeEnv":
        if type(node_id) is not int or node_id < 0:
            raise GraphValidationError("scope binding key must be a non-negative int")
        return ScopeEnv(((node_id, value),) + self.frames)

    def lookup(self, node_id: int) -> Any:
        for key, value in self.frames:
            if key == node_id:
                return value
        return _MISSING


def _copy(value: Any, *, label: str) -> Any:
    try:
        return copy.deepcopy(value)
    except Exception as exc:  # pragma: no cover - defensive boundary
        raise GraphValidationError(f"{label} is not copyable") from exc


def _is_node_id(value: Any) -> bool:
    return type(value) is int and value >= 0


def _node_id(value: Any, *, label: str, count: int) -> int:
    if not _is_node_id(value) or value >= count:
        raise GraphValidationError(
            f"{label} must be a node id in [0, {count}), got {value!r}"
        )
    return value


def _graph_nodes(graph: Any) -> tuple[list[Any], Any]:
    candidate = graph
    graph_hash = None
    if isinstance(graph, Mapping):
        graph_hash = graph.get("graph_hash")
        candidate = _MISSING
        for key in ("exact_graph_nodes", "graph_nodes", "nodes"):
            if key in graph:
                candidate = graph[key]
                break
        if candidate is _MISSING:
            raise GraphValidationError(
                "graph must contain exact_graph_nodes (or graph_nodes/nodes)"
            )
    if isinstance(candidate, (str, bytes, bytearray)) or not isinstance(
        candidate, (list, tuple)
    ):
        raise GraphValidationError("graph nodes must be an indexed list or tuple")
    if not candidate:
        raise GraphValidationError("graph nodes must be non-empty")
    if len(candidate) > _MAX_NODES:
        raise GraphValidationError("graph node limit exceeded")
    return list(candidate), _copy(graph_hash, label="graph_hash")


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
        return None
    return value


def _field(record: Mapping[str, Any], names: tuple[str, ...], *, node_id: int) -> Any:
    present = [(name, record[name]) for name in names if name in record]
    if not present:
        raise GraphValidationError(f"node {node_id}: missing {'/'.join(names)}")
    first_name, first = present[0]
    for name, value in present[1:]:
        if value != first:
            raise GraphValidationError(
                f"node {node_id}: conflicting fields {first_name} and {name}"
            )
    return first


def _ref(record: Mapping[str, Any], field: str, *, node_id: int, count: int) -> int:
    return _node_id(record.get(field, _MISSING), label=f"node {node_id}.{field}", count=count)


def _member(owner: Any, names: tuple[str, ...]) -> Any:
    if isinstance(owner, Mapping):
        for name in names:
            value = owner.get(name, _MISSING)
            if value is not _MISSING:
                return value
        return _MISSING
    for name in names:
        try:
            return getattr(owner, name)
        except Exception:
            continue
    return _MISSING


def _validate_scalar(value: Any, backend: Any, *, node_id: int, source: str) -> None:
    """Reject metadata-shaped, boolean, complex, and non-finite values."""

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
    if isinstance(value, complex):
        raise ExecutionError(
            f"{source} returned an unsupported complex value at node {node_id}",
            code=f"{source}_invalid_scalar",
            node_id=node_id,
        )
    if isinstance(value, float) and not math.isfinite(value):
        raise ExecutionError(
            f"{source} returned a non-finite float at node {node_id}",
            code=f"{source}_non_finite",
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

    validator = _member(backend, ("validate_scalar",))
    if callable(validator):
        try:
            accepted = validator(value)
        except Exception as exc:
            raise ExecutionError(
                f"backend scalar validation failed at node {node_id}",
                code="backend_invalid_scalar",
                node_id=node_id,
            ) from exc
        if accepted is True:
            return
    elif isinstance(value, numbers.Number):
        # A backend without a domain validator may use the standard numeric
        # scalar boundary.  When it does provide a validator, the validator
        # above is authoritative even for numbers.Real/Fraction values.
        return
    raise ExecutionError(
        f"{source} returned an unsupported scalar at node {node_id}",
        code=f"{source}_invalid_scalar",
        node_id=node_id,
    )


def _backend_method(backend: Any, names: tuple[str, ...], *, node_id: int, code: str) -> Any:
    method = _member(backend, names)
    if not callable(method):
        raise ExecutionError(
            f"backend method {'/'.join(names)} is required for node {node_id}",
            code=code,
            node_id=node_id,
        )
    return method


def _backend_call(
    backend: Any,
    names: tuple[str, ...],
    args: tuple[Any, ...],
    *,
    node_id: int,
    code: str,
) -> Any:
    method = _backend_method(backend, names, node_id=node_id, code=code)
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
    _validate_scalar(value, backend, node_id=node_id, source="backend")
    return value


def _backend_positive(backend: Any, value: Any, *, node_id: int) -> bool:
    method = _backend_method(
        backend, ("is_positive", "positive"), node_id=node_id, code="backend_positive_required"
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
            f"backend positivity check did not return bool at node {node_id}",
            code="backend_positive_invalid",
            node_id=node_id,
        )
    return result


def _opaque_category(operation: str | None) -> str:
    if operation == "bound_variable":
        return "variable"
    if operation in {"definite_integral", "function_substitution"}:
        # These are handled before this helper is called.
        return operation
    if operation in _OPAQUE_CATEGORIES:
        return _OPAQUE_CATEGORIES[operation]
    return "unknown"


def _oracle_method(oracle: Any, category: str, *, node_id: int) -> Any:
    if callable(oracle):
        return oracle
    if isinstance(oracle, Mapping):
        method = oracle.get(category, _MISSING)
        if callable(method):
            return method
    method = _member(oracle, (category,))
    if callable(method):
        return method
    raise OpaqueOracleRequired(node_id, category)


def _call_oracle(
    oracle: Any,
    node_id: int,
    record: Any,
    category: str,
    env: ScopeEnv,
    backend: Any,
) -> Any:
    method = _oracle_method(oracle, category, node_id=node_id)
    try:
        value = method(node_id, _copy(record, label=f"node {node_id} oracle record"), category, env)
    except OpaqueOracleRequired:
        raise
    except ExecutionError:
        raise
    except Exception as exc:
        raise ExecutionError(
            f"oracle failed for {category} node {node_id}: {type(exc).__name__}",
            code="oracle_error",
            node_id=node_id,
        ) from exc
    if type(value) is int and value == node_id:
        raise ExecutionError(
            f"oracle returned node id {node_id} as its own value",
            code="oracle_node_id_value",
            node_id=node_id,
        )
    _validate_scalar(value, backend, node_id=node_id, source="oracle")
    return value


class ScopedEvaluator:
    """Evaluate one root with explicit graph-index and lexical scope inputs."""

    def __init__(
        self,
        graph: Any,
        backend: Any,
        *,
        free_variables: Mapping[int, Any] | None = None,
        integral_binders: Mapping[int, int] | None = None,
        oracle: Any = None,
        max_depth: int = _MAX_DEPTH,
    ) -> None:
        self.nodes, self.graph_hash = _graph_nodes(graph)
        if type(max_depth) is not int or max_depth <= 0:
            raise GraphValidationError("max_depth must be a positive int")
        self.backend = backend
        self.oracle = oracle
        self.max_depth = max_depth
        if free_variables is not None and not isinstance(free_variables, Mapping):
            raise GraphValidationError("free_variables must be a node-id mapping")
        if integral_binders is not None and not isinstance(integral_binders, Mapping):
            raise GraphValidationError("integral_binders must be a node-id mapping")
        self.free_variables = dict(free_variables or {})
        self.integral_binders = dict(integral_binders or {})
        self._active: set[int] = set()
        self._opaque_nodes_used: list[int] = []
        self._initial_env = self._make_initial_env()

    def _make_initial_env(self) -> ScopeEnv:
        env = ScopeEnv()
        for key, value in self.free_variables.items():
            node_id = _node_id(key, label="free_variables key", count=len(self.nodes))
            row = self.nodes[node_id]
            if _operation(row, node_id=node_id) != "bound_variable":
                raise GraphValidationError(
                    f"free_variables key {node_id} is not a bound_variable row"
                )
            _validate_scalar(value, self.backend, node_id=node_id, source="free_variable")
            env = env.bind(node_id, value)
        return env

    def evaluate(self, root: Any) -> dict[str, Any]:
        root_id = _node_id(root, label="root", count=len(self.nodes))
        self._active.clear()
        self._opaque_nodes_used.clear()
        value = self._eval(root_id, self._initial_env, 0)
        return {
            "schema": RESULT_SCHEMA,
            "status": "conditional",
            "ok": True,
            "conditional": True,
            "certified": False,
            "graph_hash": _copy(self.graph_hash, label="graph_hash"),
            "root": root_id,
            "value": value,
            "opaque_nodes_used": list(self._opaque_nodes_used),
            "claims": dict(_CLAIMS_FALSE),
            "provenance": {
                "source_identity_verified": False,
                "backend_verified": False,
                "oracle_verified": False,
                "graph_hash_authenticated": False,
                "lexical_binder_supported": True,
                "integral_runtime_supported": True,
                "substitution_runtime_supported": True,
                "claim": "conditional scoped backend result only",
            },
        }

    def _eval(self, node_id: int, env: ScopeEnv, depth: int) -> Any:
        if depth > self.max_depth:
            raise CycleError(f"evaluation depth exceeded at node {node_id}")
        if node_id in self._active:
            raise CycleError(f"cycle detected through source node {node_id}")
        self._active.add(node_id)
        try:
            record = self.nodes[node_id]
            operation = _operation(record, node_id=node_id)
            if operation == "bound_variable":
                return self._eval_bound_variable(node_id, record, env)
            if operation == "exact_rational":
                return self._eval_exact_rational(node_id, record)
            if operation in {"sum", "product"}:
                return self._eval_nary(node_id, record, operation, env, depth)
            if operation == "negative":
                child = _ref(record, "argument", node_id=node_id, count=len(self.nodes))
                value = self._eval(child, env, depth + 1)
                minus_one = _backend_call(
                    self.backend,
                    ("rational",),
                    (-1, 1),
                    node_id=node_id,
                    code="backend_rational_required",
                )
                return _backend_call(
                    self.backend,
                    ("mul", "multiply", "times"),
                    (minus_one, value),
                    node_id=node_id,
                    code="backend_multiply_required",
                )
            if operation == "positive_quotient":
                return self._eval_quotient(node_id, record, env, depth)
            if operation == "analytic_unary":
                return self._eval_unary(node_id, record, env, depth)
            if operation == "function_substitution":
                return self._eval_substitution(node_id, record, env, depth)
            if operation == "definite_integral":
                return self._eval_integral(node_id, record, env, depth)
            return self._eval_opaque(node_id, record, _opaque_category(operation), env)
        finally:
            self._active.remove(node_id)

    def _eval_bound_variable(self, node_id: int, record: Any, env: ScopeEnv) -> Any:
        value = env.lookup(node_id)
        if value is not _MISSING:
            _validate_scalar(value, self.backend, node_id=node_id, source="scope")
            return value
        return self._eval_opaque(node_id, record, "variable", env)

    def _eval_exact_rational(self, node_id: int, record: Any) -> Any:
        if not isinstance(record, Mapping):
            raise GraphValidationError(f"node {node_id}: exact_rational row is not a mapping")
        has_numerator = "numerator" in record
        has_value = "value" in record
        if has_numerator and has_value and record["numerator"] != record["value"]:
            raise GraphValidationError(f"node {node_id}: numerator and value disagree")
        numerator = record.get("numerator", record.get("value", _MISSING))
        denominator = record.get("denominator", 1)
        if type(numerator) is not int or type(denominator) is not int:
            raise GraphValidationError(f"node {node_id}: exact rational needs integer fields")
        if denominator == 0:
            raise GraphValidationError(f"node {node_id}: exact rational denominator is zero")
        return _backend_call(
            self.backend,
            ("rational",),
            (numerator, denominator),
            node_id=node_id,
            code="backend_rational_required",
        )

    def _eval_nary(
        self, node_id: int, record: Any, operation: str, env: ScopeEnv, depth: int
    ) -> Any:
        if not isinstance(record, Mapping):
            raise GraphValidationError(f"node {node_id}: {operation} row is not a mapping")
        names = ("arguments", "terms") if operation == "sum" else ("arguments", "factors")
        raw_children = _field(record, names, node_id=node_id)
        if not isinstance(raw_children, (list, tuple)) or not raw_children:
            raise GraphValidationError(f"node {node_id}: {operation} needs children")
        children = [
            _node_id(item, label=f"node {node_id} {operation}[{i}]", count=len(self.nodes))
            for i, item in enumerate(raw_children)
        ]
        value = self._eval(children[0], env, depth + 1)
        method_names = ("add", "plus") if operation == "sum" else ("mul", "multiply", "times")
        code = "backend_add_required" if operation == "sum" else "backend_multiply_required"
        for child in children[1:]:
            value = _backend_call(
                self.backend,
                method_names,
                (value, self._eval(child, env, depth + 1)),
                node_id=node_id,
                code=code,
            )
        return value

    def _eval_quotient(self, node_id: int, record: Any, env: ScopeEnv, depth: int) -> Any:
        if not isinstance(record, Mapping):
            raise GraphValidationError(f"node {node_id}: positive_quotient row is not a mapping")
        numerator_id = _ref(record, "numerator", node_id=node_id, count=len(self.nodes))
        denominator_id = _ref(record, "denominator", node_id=node_id, count=len(self.nodes))
        numerator = self._eval(numerator_id, env, depth + 1)
        denominator = self._eval(denominator_id, env, depth + 1)
        if not _backend_positive(self.backend, denominator, node_id=node_id):
            raise ExecutionError(
                f"quotient denominator is not strictly positive at node {node_id}",
                code="quotient_denominator_not_positive",
                node_id=node_id,
            )
        return _backend_call(
            self.backend,
            ("positive_div", "positive_divide", "divide", "div", "quotient"),
            (numerator, denominator),
            node_id=node_id,
            code="backend_division_required",
        )

    def _eval_unary(self, node_id: int, record: Any, env: ScopeEnv, depth: int) -> Any:
        if not isinstance(record, Mapping):
            raise GraphValidationError(f"node {node_id}: analytic_unary row is not a mapping")
        function = _field(record, ("function", "name"), node_id=node_id)
        if not isinstance(function, str) or not function.strip():
            raise GraphValidationError(f"node {node_id}: unary function must be a string")
        function = function.strip().lower()
        if function not in _SUPPORTED_ANALYTIC_FUNCTIONS:
            return self._eval_opaque(node_id, record, "unknown", env)
        argument_id = _ref(record, "argument", node_id=node_id, count=len(self.nodes))
        value = self._eval(argument_id, env, depth + 1)
        return _backend_call(
            self.backend,
            ("unary",),
            (function, value),
            node_id=node_id,
            code="backend_unary_required",
        )

    def _eval_substitution(
        self, node_id: int, record: Any, env: ScopeEnv, depth: int
    ) -> Any:
        if not isinstance(record, Mapping):
            raise GraphValidationError(
                f"node {node_id}: function_substitution row is not a mapping"
            )
        expression_id = _ref(record, "expression", node_id=node_id, count=len(self.nodes))
        variable_id = _ref(record, "variable", node_id=node_id, count=len(self.nodes))
        value_id = _ref(record, "value", node_id=node_id, count=len(self.nodes))
        variable_row = self.nodes[variable_id]
        if _operation(variable_row, node_id=variable_id) != "bound_variable":
            raise ScopeResolutionError(
                f"substitution node {node_id}: variable {variable_id} is not bound_variable"
            )
        value = self._eval(value_id, env, depth + 1)
        return self._eval(expression_id, env.bind(variable_id, value), depth + 1)

    def _eval_integral(self, node_id: int, record: Any, env: ScopeEnv, depth: int) -> Any:
        if not isinstance(record, Mapping):
            raise GraphValidationError(f"node {node_id}: definite_integral row is not a mapping")
        integrand_id = _ref(record, "integrand", node_id=node_id, count=len(self.nodes))
        lower_id = _ref(record, "lower", node_id=node_id, count=len(self.nodes))
        upper_id = _ref(record, "upper", node_id=node_id, count=len(self.nodes))
        variable_name = record.get("variable", _MISSING)
        if not isinstance(variable_name, str) or not variable_name:
            raise ScopeResolutionError(
                f"integral node {node_id}: variable must be a non-empty string"
            )
        if node_id not in self.integral_binders:
            raise ScopeResolutionError(
                f"integral node {node_id}: explicit binder sidecar entry is required"
            )
        binder_id = _node_id(
            self.integral_binders[node_id],
            label=f"integral_binders[{node_id}]",
            count=len(self.nodes),
        )
        binder_row = self.nodes[binder_id]
        if _operation(binder_row, node_id=binder_id) != "bound_variable":
            raise ScopeResolutionError(
                f"integral node {node_id}: sidecar binder {binder_id} is not bound_variable"
            )
        binder_name = binder_row.get("name") if isinstance(binder_row, Mapping) else _MISSING
        if binder_name != variable_name:
            raise ScopeResolutionError(
                f"integral node {node_id}: sidecar binder name does not match variable"
            )

        # The raw integral row carries only a string.  Even with a sidecar,
        # accepting one of several same-name bound-variable ids could let the
        # integrand resolve a different id (including an outer free binding)
        # without an error.  Scan the complete authenticated node list and
        # fail closed unless this name has exactly one bound-variable identity.
        same_name_ids = [
            candidate_id
            for candidate_id, candidate_row in enumerate(self.nodes)
            if _operation(candidate_row, node_id=candidate_id) == "bound_variable"
            and isinstance(candidate_row, Mapping)
            and candidate_row.get("name", _MISSING) == variable_name
        ]
        if len(same_name_ids) != 1 or same_name_ids[0] != binder_id:
            raise ScopeResolutionError(
                f"integral node {node_id}: variable name is not a unique sidecar identity"
            )

        # Both endpoints are evaluated in the outer environment.  The
        # backend receives scalar endpoints and a closure, never endpoint ids.
        lower = self._eval(lower_id, env, depth + 1)
        upper = self._eval(upper_id, env, depth + 1)
        metadata = {
            key: _copy(value, label=f"integral node {node_id} metadata")
            for key, value in record.items()
            if key not in {"operation", "op", "integrand", "lower", "upper", "variable"}
        }

        def body(argument: Any) -> Any:
            _validate_scalar(
                argument,
                self.backend,
                node_id=node_id,
                source="integral_argument",
            )
            return self._eval(integrand_id, env.bind(binder_id, argument), depth + 1)

        return _backend_call(
            self.backend,
            ("integral",),
            (lower, upper, body, metadata),
            node_id=node_id,
            code="backend_integral_required",
        )

    def _eval_opaque(self, node_id: int, record: Any, category: str, env: ScopeEnv) -> Any:
        if self.oracle is None:
            raise OpaqueOracleRequired(node_id, category)
        value = _call_oracle(self.oracle, node_id, record, category, env, self.backend)
        self._opaque_nodes_used.append(node_id)
        return value


def evaluate(
    graph: Any,
    root: Any,
    backend: Any,
    *,
    free_variables: Mapping[int, Any] | None = None,
    integral_binders: Mapping[int, int] | None = None,
    oracle: Any = None,
    max_depth: int = _MAX_DEPTH,
) -> dict[str, Any]:
    """Evaluate one graph root with explicit lexical and opaque boundaries."""

    evaluator = ScopedEvaluator(
        graph,
        backend,
        free_variables=free_variables,
        integral_binders=integral_binders,
        oracle=oracle,
        max_depth=max_depth,
    )
    return evaluator.evaluate(root)


def execute(
    graph: Any,
    root: Any,
    backend: Any,
    *,
    free_variables: Mapping[int, Any] | None = None,
    integral_binders: Mapping[int, int] | None = None,
    oracle: Any = None,
    max_depth: int = _MAX_DEPTH,
) -> dict[str, Any]:
    """Compatibility spelling for callers that use an execute entry point."""

    return evaluate(
        graph,
        root,
        backend,
        free_variables=free_variables,
        integral_binders=integral_binders,
        oracle=oracle,
        max_depth=max_depth,
    )


__all__ = [
    "SCHEMA",
    "RESULT_SCHEMA",
    "ScopedGraphError",
    "GraphValidationError",
    "ScopeResolutionError",
    "CycleError",
    "ExecutionError",
    "OpaqueOracleRequired",
    "ScopeEnv",
    "ScopedEvaluator",
    "evaluate",
    "execute",
]
