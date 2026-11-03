"""Restore the saved reserved C1 graph without invoking its old producer.

The caller supplies JSON documents whose *raw gzip bytes* were authenticated
by the runner before decoding.  This consumer does not hash or decode raw
bytes; it only requires the runner's authentication marker and checks its
declared archive pin.  It checks the saved graph relationship, and restores
opaque rows into a small same-graph transport object.  It interprets only the
dependency closures of ``mu``, the reserved matrix, and its saved inverse.
Every other full-graph row is retained but intentionally opaque.

The restored source owner is the newly created graph object.  It is not the
historical live source owner, and the saved inverse remains a caller/archive
premise; no inverse identity is proved here.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass
from fractions import Fraction
import json
from types import MappingProxyType
from typing import Any, Final

from .reserved_matrix_c1_action_v1 import MatrixC1Action, inverse_c1_action


SCHEMA: Final[str] = "reserved_matrix_actual_graph_v1"
RESERVED_ARCHIVE_SHA256: Final[str] = (
    "98fd45a48d72656ec91f19a5564577baed0dbb43fe85ce8c842336856dc17a17"
)
FULL_ARCHIVE_SHA256: Final[str] = (
    "d6b193bf527d4e8cd4afe415b1e941c928f32288b63362783074381ef7480fc1"
)
RESERVED_GRAPH_NODE_COUNT: Final[int] = 25_199
FULL_GRAPH_NODE_COUNT: Final[int] = 53_398
INTEGRAL_VARIABLE: Final[str] = "n1_reserved_normalized_bump_t"

M_NODE_IDS: Final[tuple[tuple[int, int, int], ...]] = (
    (25044, 25046, 25048),
    (25057, 25060, 25063),
    (25072, 25075, 25078),
)
INVERSE_NODE_IDS: Final[tuple[tuple[int, int, int], ...]] = (
    (25106, 25112, 25117),
    (25123, 25128, 25134),
    (25139, 25145, 25150),
)
Y_NODE_IDS: Final[tuple[int, int, int]] = (53204, 53214, 53222)
Y_Z_NODE_IDS: Final[tuple[int, int, int]] = (53205, 53215, 53223)
MOMENT_NODE_IDS: Final[Mapping[str, tuple[int, int]]] = MappingProxyType(
    {
        "m11": (53145, 53146),
        "m12": (53149, 53150),
        "m13": (53169, 53170),
        "m14": (53157, 53158),
        "m15": (53175, 53176),
    }
)
COEFFICIENT_NODE_IDS: Final[Mapping[str, tuple[int, int]]] = MappingProxyType(
    {
        "alpha1": (53228, 53250),
        "alpha2": (53233, 53255),
        "beta1": (53237, 53259),
        "beta2": (53241, 53263),
        "beta3": (53245, 53267),
    }
)


class ActualGraphError(ValueError):
    """Base error for unauthenticated or malformed saved graph input."""


class ArchiveAuthenticationError(ActualGraphError):
    """Raised when runner-owned raw-SHA authentication metadata is missing."""


class GraphRestoreError(ActualGraphError):
    """Raised when the saved rows cannot preserve transport-graph semantics."""


class DependencyError(ActualGraphError):
    """Raised when a checked dependency is unknown or is not Z-independent."""


@dataclass(frozen=True, slots=True)
class AuthenticatedDecodedArchive:
    """A decoded archive plus the external raw-byte authentication result.

    ``raw_sha256_verified`` is a runner-owned assertion made only after the
    complete raw gzip bytes have been checked.  This consumer never hashes or
    decompresses anything itself.
    """

    document: Mapping[str, Any]
    raw_sha256: str
    raw_sha256_verified: bool


def authenticated_archive(document: Mapping[str, Any], raw_sha256: str) -> AuthenticatedDecodedArchive:
    """Wrap a document after the caller has checked its complete raw bytes."""

    return AuthenticatedDecodedArchive(document, raw_sha256, True)


@dataclass(frozen=True, slots=True)
class RestoredFunctionRef:
    """Opaque same-graph function handle used by the C1 consumer."""

    graph: "FunctionTransportGraphAdapter"
    node: int


@dataclass(frozen=True, slots=True)
class RestoredC1Function:
    value: RestoredFunctionRef
    Z: RestoredFunctionRef


class FunctionTransportGraphAdapter:
    """Minimal saved-graph transport matching the native interning rules."""

    def __init__(self, nodes: Sequence[Mapping[str, Any]]) -> None:
        rows: list[dict[str, Any]] = []
        for index, row in enumerate(nodes):
            if not isinstance(row, Mapping):
                raise GraphRestoreError(f"graph row {index} is not a JSON object")
            rows.append(dict(row))
        if len(rows) < 2:
            raise GraphRestoreError("restored graph must contain zero and one rows")
        if rows[0] != {"operation": "exact_rational", "numerator": 0, "denominator": 1}:
            raise GraphRestoreError("graph row zero is not the exact zero handle")
        if rows[1] != {"operation": "exact_rational", "numerator": 1, "denominator": 1}:
            raise GraphRestoreError("graph row one is not the exact one handle")
        self.nodes = rows
        self.keys: dict[str, int] = {}
        for index, row in enumerate(rows):
            key = self._key(row)
            if key in self.keys:
                raise GraphRestoreError(f"duplicate interned graph row at {index}")
            self.keys[key] = index
        self.zero = RestoredFunctionRef(self, 0)
        self.one = RestoredFunctionRef(self, 1)
        self.restored_node_count = len(rows)

    @staticmethod
    def _key(row: Mapping[str, Any]) -> str:
        try:
            return json.dumps(dict(row), sort_keys=True, separators=(",", ":"))
        except (TypeError, ValueError) as exc:
            raise GraphRestoreError("graph row is not JSON-internable") from exc

    def _ids(self, arguments: Sequence[RestoredFunctionRef]) -> list[int]:
        if any(type(value) is not RestoredFunctionRef or value.graph is not self for value in arguments):
            raise TypeError("same-restored-graph function handles required")
        return [value.node for value in arguments]

    def node(self, operation: str, **attributes: Any) -> RestoredFunctionRef:
        if type(operation) is not str or not operation:
            raise GraphRestoreError("operation name must be a nonempty string")
        row = dict(operation=operation, **attributes)
        key = self._key(row)
        existing = self.keys.get(key)
        if existing is not None:
            return RestoredFunctionRef(self, existing)
        node = len(self.nodes)
        self.nodes.append(row)
        self.keys[key] = node
        return RestoredFunctionRef(self, node)

    def constant(self, value: Any) -> RestoredFunctionRef:
        fraction = Fraction(value)
        return self.node("exact_rational", numerator=fraction.numerator, denominator=fraction.denominator)

    def add(self, *arguments: RestoredFunctionRef) -> RestoredFunctionRef:
        ids = self._ids(arguments)
        ids = [value for value in ids if value != self.zero.node]
        if not ids:
            return self.zero
        if len(ids) == 1:
            return RestoredFunctionRef(self, ids[0])
        return self.node("sum", arguments=ids)

    def mul(self, *arguments: RestoredFunctionRef) -> RestoredFunctionRef:
        ids = self._ids(arguments)
        if self.zero.node in ids:
            return self.zero
        ids = [value for value in ids if value != self.one.node]
        if not ids:
            return self.one
        if len(ids) == 1:
            return RestoredFunctionRef(self, ids[0])
        return self.node("product", arguments=ids)

    def neg(self, value: RestoredFunctionRef) -> RestoredFunctionRef:
        self._ids((value,))
        if value.node == self.zero.node:
            return self.zero
        return self.node("negative", argument=value.node)


@dataclass(frozen=True, slots=True)
class DependencyProof:
    roots: tuple[int, ...]
    node_count: int
    operations: tuple[str, ...]
    Z_independent: bool


@dataclass(frozen=True, slots=True)
class SymbolicCoefficientMatch:
    checked: Mapping[str, bool | None]
    methods: Mapping[str, str]
    complete: bool
    status: str


@dataclass(frozen=True, slots=True)
class ActualReservedGraph:
    graph: FunctionTransportGraphAdapter
    source_family: Mapping[str, Any]
    M: tuple[tuple[RestoredFunctionRef, ...], ...]
    M_Z: tuple[tuple[RestoredFunctionRef, ...], ...]
    M_inverse: tuple[tuple[RestoredFunctionRef, ...], ...]
    y: tuple[RestoredFunctionRef, ...]
    y_Z: tuple[RestoredFunctionRef, ...]
    beta: tuple[RestoredFunctionRef, ...]
    beta_Z: tuple[RestoredFunctionRef, ...]
    beta_pairs: tuple[Any, ...]
    coefficients: Mapping[str, RestoredC1Function]
    action: MatrixC1Action
    coefficient_match: SymbolicCoefficientMatch
    mu_dependency: DependencyProof
    matrix_dependency: DependencyProof
    inverse_dependency: DependencyProof
    reserved_prefix_verified: bool
    source_owner: str
    actual_source_owner_authenticated: bool
    inverse_identity_verified: bool


def _unwrap_archive(
    archive: Mapping[str, Any] | AuthenticatedDecodedArchive,
    supplied_sha256: str | None,
    raw_sha256_verified: bool,
    expected_sha256: str,
    label: str,
) -> Mapping[str, Any]:
    if isinstance(archive, AuthenticatedDecodedArchive):
        document = archive.document
        digest = archive.raw_sha256
        verified = archive.raw_sha256_verified
        if supplied_sha256 is not None and supplied_sha256.lower() != digest.lower():
            raise ArchiveAuthenticationError(f"{label} supplied SHA differs from archive wrapper")
    else:
        document = archive
        digest = supplied_sha256
        verified = raw_sha256_verified
    if not verified or not isinstance(digest, str) or digest.lower() != expected_sha256:
        raise ArchiveAuthenticationError(f"{label} must carry externally verified raw SHA {expected_sha256}")
    if not isinstance(document, Mapping):
        raise GraphRestoreError(f"{label} decoded document is not an object")
    return document


def _nodes(document: Mapping[str, Any], label: str) -> list[Mapping[str, Any]]:
    rows = document.get("exact_graph_nodes")
    if not isinstance(rows, list):
        raise GraphRestoreError(f"{label}.exact_graph_nodes must be a list")
    return rows


def _require_node(rows: Sequence[Mapping[str, Any]], node: int, label: str) -> Mapping[str, Any]:
    if type(node) is not int or node < 0 or node >= len(rows):
        raise GraphRestoreError(f"{label} node id is outside restored graph")
    row = rows[node]
    if not isinstance(row, Mapping):
        raise GraphRestoreError(f"{label} row is not an object")
    return row


def _require_field(row: Mapping[str, Any], key: str, expected: Any, label: str) -> None:
    if row.get(key) != expected:
        raise GraphRestoreError(f"{label}.{key} has an unexpected saved semantic reference")


def _semantic_edges(
    row: Mapping[str, Any], label: str, scope: frozenset[str]
) -> tuple[tuple[int, frozenset[str]], ...]:
    operation = row.get("operation")
    if operation == "exact_rational":
        if type(row.get("numerator")) is not int or type(row.get("denominator")) is not int:
            raise DependencyError(f"{label} exact rational has non-integer literals")
        if row["denominator"] == 0:
            raise DependencyError(f"{label} exact rational has zero denominator")
        return ()
    if operation in {"sum", "product"}:
        values = row.get("arguments")
        if not isinstance(values, list) or any(type(value) is not int for value in values):
            raise DependencyError(f"{label} {operation} has malformed semantic references")
        return tuple((value, scope) for value in values)
    if operation in {"negative", "analytic_unary", "compact_raw_beta"}:
        value = row.get("argument")
        if type(value) is not int:
            raise DependencyError(f"{label} {operation} has no integer argument reference")
        if operation == "analytic_unary" and row.get("name") not in {"exp", "exprel", "log"}:
            raise DependencyError(f"{label} analytic unary is not a known reserved operation")
        return ((value, scope),)
    if operation == "positive_quotient":
        numerator, denominator = row.get("numerator"), row.get("denominator")
        if type(numerator) is not int or type(denominator) is not int:
            raise DependencyError(f"{label} quotient has malformed semantic references")
        return ((numerator, scope), (denominator, scope))
    if operation == "definite_integral":
        integrand, lower, upper = (row.get(key) for key in ("integrand", "lower", "upper"))
        variable = row.get("variable")
        if (
            any(type(value) is not int for value in (integrand, lower, upper))
            or variable != INTEGRAL_VARIABLE
        ):
            raise DependencyError(f"{label} integral has malformed semantic references")
        # Bounds are evaluated in the outer lexical scope.  Only the
        # integrand sees the newly introduced integration variable.
        bound_scope = scope | {variable}
        return ((lower, scope), (upper, scope), (integrand, bound_scope))
    if operation == "bound_variable":
        name = row.get("name")
        if not isinstance(name, str) or name not in scope:
            raise DependencyError(f"{label} contains a free or unknown bound variable")
        return ()
    raise DependencyError(f"{label} contains opaque or unknown operation {operation!r}")


def _dependency_proof(
    rows: Sequence[Mapping[str, Any]], roots: Sequence[int], label: str
) -> DependencyProof:
    seen: set[tuple[int, frozenset[str]]] = set()
    pending = [(node, frozenset()) for node in roots]
    operations: set[str] = set()
    while pending:
        node, scope = pending.pop()
        state = (node, scope)
        if state in seen:
            continue
        row = _require_node(rows, node, f"{label}[{node}]")
        if row.get("Z_independent") is not None and row.get("Z_independent") is not True:
            raise DependencyError(f"{label}[{node}] is not marked Z-independent")
        seen.add(state)
        operation = row.get("operation")
        if not isinstance(operation, str):
            raise DependencyError(f"{label}[{node}] has no operation name")
        operations.add(operation)
        pending.extend(_semantic_edges(row, f"{label}[{node}]", scope))
    return DependencyProof(tuple(roots), len(seen), tuple(sorted(operations)), True)


def _outer_sum_atoms(graph: FunctionTransportGraphAdapter, node: int) -> tuple[int, ...]:
    """Flatten only outer ``sum`` nodes; every other node is an ID atom.

    This deliberately does not normalize products or recursively inspect
    integrals, source recipes, or moment DAGs.  Those expressions are already
    interned in one restored graph, so equal product/source expressions have
    the same node ID.  The only permitted non-identity comparison is the
    associativity of an outer sum.
    """

    atoms: list[int] = []
    pending = [node]
    while pending:
        current = pending.pop()
        row = _require_node(graph.nodes, current, f"comparison[{current}]")
        if row.get("operation") != "sum":
            atoms.append(current)
            continue
        arguments = row.get("arguments")
        if not isinstance(arguments, list) or any(type(value) is not int for value in arguments):
            raise DependencyError(f"comparison[{current}] has malformed sum references")
        pending.extend(reversed(arguments))
    return tuple(atoms)


def _symbolic_match(
    graph: FunctionTransportGraphAdapter, actual: RestoredFunctionRef, expected: int
) -> tuple[bool, str]:
    if actual.graph is not graph:
        return False, "different_graph"
    if actual.node == expected:
        return True, "saved_node_identity"
    actual_row = _require_node(graph.nodes, actual.node, f"comparison[{actual.node}]")
    expected_row = _require_node(graph.nodes, expected, f"comparison[{expected}]")
    if actual_row.get("operation") != "sum" or expected_row.get("operation") != "sum":
        return False, "node_id_atom"
    try:
        return _outer_sum_atoms(graph, actual.node) == _outer_sum_atoms(graph, expected), "outer_sum_associativity"
    except (ActualGraphError, KeyError, TypeError):
        return False, "unresolved_sum_expression"


def _validate_y_nodes(rows: Sequence[Mapping[str, Any]], moments: Mapping[str, Any]) -> None:
    expected_moments = {key: list(value) for key, value in MOMENT_NODE_IDS.items()}
    if dict(moments) != expected_moments:
        raise GraphRestoreError("full archive moment handles are not the saved signed family")
    descriptors = {
        53204: {"operation": "product", "arguments": [53203, 53149]},
        53205: {"operation": "product", "arguments": [53203, 53150]},
        53214: {"operation": "product", "arguments": [53213, 53206]},
        53215: {"operation": "product", "arguments": [53213, 53211]},
        53222: {"operation": "product", "arguments": [53189, 53216]},
        53223: {"operation": "product", "arguments": [53189, 53221]},
    }
    for node, descriptor in descriptors.items():
        row = _require_node(rows, node, f"signed y[{node}]")
        for key, value in descriptor.items():
            _require_field(row, key, value, f"signed y[{node}]")


def restore_actual_graph(
    reserved_archive: Mapping[str, Any] | AuthenticatedDecodedArchive,
    full_archive: Mapping[str, Any] | AuthenticatedDecodedArchive,
    *,
    reserved_raw_sha256: str | None = None,
    full_raw_sha256: str | None = None,
    raw_sha256_verified: bool = False,
    action: Callable[..., MatrixC1Action] = inverse_c1_action,
) -> ActualReservedGraph:
    """Restore and bind the real saved reserved matrix C1 graph.

    Inputs are already decoded JSON documents.  Either wrap each document in
    :class:`AuthenticatedDecodedArchive`, or pass the runner's externally
    verified raw-SHA label through the corresponding keyword and set
    ``raw_sha256_verified``.  The runner, not this consumer, hashes and
    authenticates raw bytes; no file is opened or decompressed here.
    """

    reserved = _unwrap_archive(
        reserved_archive, reserved_raw_sha256, raw_sha256_verified, RESERVED_ARCHIVE_SHA256, "reserved archive"
    )
    full = _unwrap_archive(
        full_archive, full_raw_sha256, raw_sha256_verified, FULL_ARCHIVE_SHA256, "full archive"
    )
    reserved_nodes = _nodes(reserved, "reserved archive")
    full_nodes = _nodes(full, "full archive")
    if len(reserved_nodes) != RESERVED_GRAPH_NODE_COUNT or len(full_nodes) != FULL_GRAPH_NODE_COUNT:
        raise GraphRestoreError("saved graph node counts do not match the pinned archives")
    if full_nodes[:RESERVED_GRAPH_NODE_COUNT] != reserved_nodes:
        raise GraphRestoreError("reserved graph is not the exact full-defect graph prefix")
    reserved_family = reserved.get("source_family")
    full_family = full.get("source_family")
    if not isinstance(reserved_family, Mapping) or reserved_family != full_family:
        raise GraphRestoreError("reserved and full archives have different source families")

    basis = reserved.get("exact_normalized_basis_and_matrix_functions")
    if not isinstance(basis, Mapping):
        raise GraphRestoreError("reserved normalized matrix binding is missing")
    if basis.get("mu") != 17 or basis.get("angular_rows") != [list(row) for row in M_NODE_IDS]:
        raise GraphRestoreError("reserved M/mu semantic references changed")
    if basis.get("inverse_angular") != [list(row) for row in INVERSE_NODE_IDS]:
        raise GraphRestoreError("reserved inverse semantic references changed")

    bindings = full.get("full_defect_linear_bindings")
    if not isinstance(bindings, Mapping):
        raise GraphRestoreError("full signed-defect bindings are missing")
    moments = bindings.get("moments")
    coefficients = bindings.get("coefficients")
    if not isinstance(moments, Mapping) or not isinstance(coefficients, Mapping):
        raise GraphRestoreError("full moments or coefficients are missing")
    _validate_y_nodes(full_nodes, moments)
    if dict(coefficients) != {key: list(value) for key, value in COEFFICIENT_NODE_IDS.items()}:
        raise GraphRestoreError("full coefficient handles are not the saved archive bindings")

    graph = FunctionTransportGraphAdapter(full_nodes)
    matrix = tuple(tuple(RestoredFunctionRef(graph, node) for node in row) for row in M_NODE_IDS)
    matrix_inverse = tuple(tuple(RestoredFunctionRef(graph, node) for node in row) for row in INVERSE_NODE_IDS)
    y = tuple(RestoredFunctionRef(graph, node) for node in Y_NODE_IDS)
    y_Z = tuple(RestoredFunctionRef(graph, node) for node in Y_Z_NODE_IDS)
    matrix_Z = tuple(tuple(graph.zero for _ in range(3)) for _ in range(3))

    mu_dependency = _dependency_proof(full_nodes, (17,), "mu")
    matrix_dependency = _dependency_proof(full_nodes, tuple(node for row in M_NODE_IDS for node in row), "M")
    inverse_dependency = _dependency_proof(
        full_nodes, tuple(node for row in INVERSE_NODE_IDS for node in row), "M_inverse"
    )

    if not callable(action):
        raise TypeError("action must be the reserved_matrix_c1_action_v1 inverse action callable")
    try:
        applied = action(
            matrix,
            matrix_inverse,
            y,
            y_Z,
            matrix_Z=matrix_Z,
            add=graph.add,
            mul=graph.mul,
            neg=graph.neg,
        )
    except Exception as exc:
        raise ActualGraphError("reserved matrix C1 action rejected restored handles") from exc
    if not isinstance(applied, MatrixC1Action):
        raise ActualGraphError("reserved matrix action did not return MatrixC1Action")

    archive_coefficients = {
        key: RestoredC1Function(RestoredFunctionRef(graph, pair[0]), RestoredFunctionRef(graph, pair[1]))
        for key, pair in COEFFICIENT_NODE_IDS.items()
    }
    checked: dict[str, bool | None] = {}
    methods: dict[str, str] = {}
    for index, key in enumerate(("beta1", "beta2", "beta3")):
        value_ok, value_method = _symbolic_match(graph, applied.beta[index], COEFFICIENT_NODE_IDS[key][0])
        z_ok, z_method = _symbolic_match(graph, applied.beta_Z[index], COEFFICIENT_NODE_IDS[key][1])
        checked[key] = value_ok and z_ok
        methods[key] = value_method if value_method == z_method else value_method + "+" + z_method
    checked["alpha1"] = None
    checked["alpha2"] = None
    methods["alpha1"] = "not_in_three_channel_action"
    methods["alpha2"] = "not_in_three_channel_action"
    angular_complete = all(checked[key] is True for key in ("beta1", "beta2", "beta3"))
    match = SymbolicCoefficientMatch(
        MappingProxyType(checked),
        MappingProxyType(methods),
        False,
        "angular_partial_match" if angular_complete else "angular_match_unresolved_or_failed",
    )
    return ActualReservedGraph(
        graph=graph,
        source_family=MappingProxyType(dict(reserved_family)),
        M=matrix,
        M_Z=matrix_Z,
        M_inverse=matrix_inverse,
        y=y,
        y_Z=y_Z,
        beta=tuple(applied.beta),
        beta_Z=tuple(applied.beta_Z),
        beta_pairs=tuple(applied.beta_pairs),
        coefficients=MappingProxyType(archive_coefficients),
        action=applied,
        coefficient_match=match,
        mu_dependency=mu_dependency,
        matrix_dependency=matrix_dependency,
        inverse_dependency=inverse_dependency,
        reserved_prefix_verified=True,
        source_owner="restored_archive",
        actual_source_owner_authenticated=False,
        inverse_identity_verified=False,
    )


restore_reserved_actual_graph = restore_actual_graph
load_reserved_actual_graph = restore_actual_graph


__all__ = [
    "ActualGraphError",
    "ActualReservedGraph",
    "ArchiveAuthenticationError",
    "AuthenticatedDecodedArchive",
    "COEFFICIENT_NODE_IDS",
    "DependencyError",
    "DependencyProof",
    "FULL_ARCHIVE_SHA256",
    "FULL_GRAPH_NODE_COUNT",
    "FunctionTransportGraphAdapter",
    "GraphRestoreError",
    "INVERSE_NODE_IDS",
    "M_NODE_IDS",
    "MOMENT_NODE_IDS",
    "RESERVED_ARCHIVE_SHA256",
    "RESERVED_GRAPH_NODE_COUNT",
    "RestoredC1Function",
    "RestoredFunctionRef",
    "SCHEMA",
    "SymbolicCoefficientMatch",
    "Y_NODE_IDS",
    "Y_Z_NODE_IDS",
    "authenticated_archive",
    "load_reserved_actual_graph",
    "restore_actual_graph",
    "restore_reserved_actual_graph",
]
