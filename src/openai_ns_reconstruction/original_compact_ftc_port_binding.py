"""Typed, bounded source-port codec for compact FTC derivations.

The codec consumes an already authenticated decoded ``exact_graph_nodes``
array.  It only checks the fixed physical rows and the closed arithmetic
slices named by a caller.  It never imports a producer, replays an archived
consumer, evaluates a profile or integral, or performs proof normalization.
The caller owns all residual and transport mathematics.
"""

from __future__ import annotations

from collections import deque
import hashlib
import heapq
import json
from typing import Any, Mapping, Sequence


CODEC_SCHEMA = "compact_ftc_port_codec_v1"
ALLOWED_OPERATIONS = (
    "exact_rational",
    "negative",
    "positive_quotient",
    "product",
    "sum",
)
ALLOWED_OPERATION_SET = frozenset(ALLOWED_OPERATIONS)
DEFAULT_SLICE_CAP = 1024

NATIVE_ARCHIVE_SHA256 = "b6a7227da5e2427561cd22fcff614003cc14f3182537e1dbe278a6de84322510"
OLD_ARCHIVE_SHA256 = "348df68bb426a2dfda6012e89f678b1e9c4ceb68b12f65f36853d5579977bc64"


class PortCodecError(ValueError):
    """A malformed typed port, arithmetic slice, or slice bound."""


class PortCodecUnsupportedOperation(PortCodecError):
    """A selected row uses an operation outside the five arithmetic ports."""


class PortCodecCycle(PortCodecError):
    """A selected arithmetic slice is cyclic."""


STATIC_NEGATIVE_FIXTURE_PLAN = (
    {
        "name": "reversed_physical_factors",
        "kind": "physical_acceptance",
        "expected": "accept_and_preserve_raw_argument_order",
        "execute_once_inside_future_derive": True,
        "executed": False,
    },
    {
        "name": "missing_physical_factor",
        "kind": "physical_rejection",
        "expected": "raise PortCodecError",
        "execute_once_inside_future_derive": True,
        "executed": False,
    },
    {
        "name": "duplicate_pstar",
        "kind": "physical_rejection",
        "expected": "raise PortCodecError",
        "execute_once_inside_future_derive": True,
        "executed": False,
    },
    {
        "name": "unsupported_operation",
        "kind": "compile_rejection",
        "expected": "raise PortCodecUnsupportedOperation",
        "execute_once_inside_future_derive": True,
        "executed": False,
    },
)


# These mappings are the fixed physical rows checked against the authenticated
# native and old decoded archives before this file was written.  Native ``my``
# is derived as U-M in its archive and therefore has no physical-product row.
AUTHENTICATED_PHYSICAL_PORTS = {
    "native": {
        "decoded_sha256": NATIVE_ARCHIVE_SHA256,
        "pstar": 33572,
        "U": {"rows": [262269, 262270, 262271, 262272, 262273, 262274], "fields": [260375, 260378, 260382, 260387, 260393, 260400]},
        "V": {"rows": [262269, 262270, 262271, 262272, 262273, 262274], "fields": [260375, 260378, 260382, 260387, 260393, 260400]},
        "U_y": {"rows": [262275, 262276, 262277, 262278, 262279, 262280], "fields": [260494, 260497, 260501, 260506, 260512, 260519]},
        "V_y": {"rows": [262275, 262276, 262277, 262278, 262279, 262280], "fields": [260494, 260497, 260501, 260506, 260512, 260519]},
        "M": {"rows": [262281, 262282, 262283, 262284, 262285, 262286], "fields": [261716, 261717, 261718, 261719, 261720, 261721]},
        "m": {"rows": [262281, 262282, 262283, 262284, 262285, 262286], "fields": [261716, 261717, 261718, 261719, 261720, 261721]},
        "my": {"rows": [262293, 262294, 262295, 262296, 262297, 262298], "derived": "U-M", "physical_product": False, "negative_parents": [262287, 262288, 262289, 262290, 262291, 262292]},
    },
    "old": {
        "decoded_sha256": OLD_ARCHIVE_SHA256,
        "pstar": 25203,
        "U": {"rows": [30709, 30710, 30711, 30712], "fields": [13182, 13185, 13310, 24596]},
        "V": {"rows": [30709, 30710, 30711, 30712], "fields": [13182, 13185, 13310, 24596]},
        "U_y": {"rows": [30713, 30714, 30715, 30716], "fields": [13231, 13232, 13644, 24837]},
        "V_y": {"rows": [30713, 30714, 30715, 30716], "fields": [13231, 13232, 13644, 24837]},
        "M": {"rows": [30717, 30718, 30719, 30720], "fields": [13227, 13230, 13344, 24627]},
        "m": {"rows": [30717, 30718, 30719, 30720], "fields": [13227, 13230, 13344, 24627]},
        "M_y": {"rows": [30721, 30722, 30723, 30724], "fields": [13247, 13248, 13658, 24850]},
        "my": {"rows": [30721, 30722, 30723, 30724], "fields": [13247, 13248, 13658, 24850]},
    },
}


def _canonical_hash(value: Any) -> str:
    encoded = json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


class PortCodec:
    """Compile fixed arithmetic slices into caller-supplied symbolic objects."""

    allowed_operations = ALLOWED_OPERATIONS

    def __init__(self, nodes: Sequence[Any], archive_sha: str, chart: str, sp: Any):
        if not isinstance(nodes, list):
            raise PortCodecError("nodes must be the decoded exact_graph_nodes list")
        if not isinstance(archive_sha, str) or len(archive_sha) != 64:
            raise PortCodecError("archive_sha must be a hexadecimal SHA-256 string")
        try:
            int(archive_sha, 16)
        except ValueError as exc:
            raise PortCodecError("archive_sha must be hexadecimal") from exc
        if not isinstance(chart, str) or not chart:
            raise PortCodecError("chart must be a non-empty string")
        if sp is None:
            raise PortCodecError("sp must be the caller-supplied symbolic arithmetic module")
        expected = AUTHENTICATED_PHYSICAL_PORTS.get(chart)
        if expected is not None and expected["decoded_sha256"] != archive_sha:
            raise PortCodecError(f"{chart}: decoded archive SHA-256 does not match fixed port mapping")
        self.nodes = nodes
        self.archive_sha = archive_sha
        self.chart = chart
        self.sp = sp
        self.audit: dict[str, Any] = {
            "physical_calls": 0,
            "compile_calls": 0,
            "rows_read": 0,
            "physical_rows_checked": 0,
            "slice_rows_checked": 0,
            "actual_expansions": 0,
            "unsupported_operations_rejected": 0,
            "invalid_ids_rejected": 0,
            "missing_fields_rejected": 0,
            "cycles_rejected": 0,
            "caps_rejected": 0,
            "profile_evaluations": 0,
            "integral_evaluations": 0,
            "source_evaluations": 0,
            "whole_archive_walks": 0,
            "proof_normalizations": 0,
        }
        self.counters = self.audit
        self.slices: list[dict[str, Any]] = []
        self.physical_receipts: list[dict[str, Any]] = []

    def _row(self, node_id: int, label: str) -> dict[str, Any]:
        if isinstance(node_id, bool) or not isinstance(node_id, int) or node_id < 0 or node_id >= len(self.nodes):
            self.audit["invalid_ids_rejected"] += 1
            raise PortCodecError(f"{label}: invalid node id {node_id!r}")
        row = self.nodes[node_id]
        self.audit["rows_read"] += 1
        if not isinstance(row, dict):
            raise PortCodecError(f"{label}: node {node_id} is not an object")
        return row

    def _children(self, node_id: int, row: Mapping[str, Any], label: str) -> list[int]:
        operation = row.get("operation")
        if operation not in ALLOWED_OPERATION_SET:
            self.audit["unsupported_operations_rejected"] += 1
            raise PortCodecUnsupportedOperation(
                f"{label}: node {node_id} operation {operation!r} is outside {list(ALLOWED_OPERATIONS)!r}"
            )
        if operation == "exact_rational":
            numerator = row.get("numerator")
            denominator = row.get("denominator")
            if isinstance(numerator, bool) or not isinstance(numerator, int) or isinstance(denominator, bool) or not isinstance(denominator, int):
                self.audit["missing_fields_rejected"] += 1
                raise PortCodecError(f"{label}: exact_rational node {node_id} lacks integer numerator/denominator")
            if denominator == 0:
                raise PortCodecError(f"{label}: exact_rational node {node_id} has zero denominator")
            return []
        if operation in {"sum", "product"}:
            arguments = row.get("arguments")
            if not isinstance(arguments, list):
                self.audit["missing_fields_rejected"] += 1
                raise PortCodecError(f"{label}: {operation} node {node_id} lacks an arguments list")
            children = list(arguments)
        elif operation == "negative":
            if "argument" not in row:
                self.audit["missing_fields_rejected"] += 1
                raise PortCodecError(f"{label}: negative node {node_id} lacks argument")
            children = [row["argument"]]
        else:
            if "numerator" not in row or "denominator" not in row:
                self.audit["missing_fields_rejected"] += 1
                raise PortCodecError(f"{label}: positive_quotient node {node_id} lacks numerator/denominator")
            children = [row["numerator"], row["denominator"]]
        for child in children:
            if isinstance(child, bool) or not isinstance(child, int):
                self.audit["invalid_ids_rejected"] += 1
                raise PortCodecError(f"{label}: node {node_id} has non-integer child {child!r}")
            self._row(child, label)
        return children

    def physical(
        self,
        rows: Sequence[int],
        pstar: int,
        expected_fields: Sequence[int],
        label: str,
    ) -> list[int]:
        """Validate Pstar times one underlying field and return field IDs."""
        if not isinstance(rows, (list, tuple)) or not isinstance(expected_fields, (list, tuple)):
            raise PortCodecError(f"{label}: rows and expected_fields must be sequences")
        if len(rows) != len(expected_fields):
            raise PortCodecError(f"{label}: rows and expected_fields have different lengths")
        if isinstance(pstar, bool) or not isinstance(pstar, int):
            raise PortCodecError(f"{label}: pstar must be an integer node id")
        if pstar in expected_fields:
            raise PortCodecError(f"{label}: Pstar node {pstar} cannot also be an expected underlying field")
        for order, expected in enumerate(expected_fields):
            if isinstance(expected, bool) or not isinstance(expected, int):
                raise PortCodecError(f"{label}[{order}]: expected underlying field must be an integer node id")
        self._row(pstar, label)
        self.audit["physical_calls"] += 1
        underlying: list[int] = []
        raw_rows: list[dict[str, Any]] = []
        for order, (row_id, expected) in enumerate(zip(rows, expected_fields)):
            row = self._row(row_id, label)
            if row.get("operation") != "product":
                raise PortCodecError(f"{label}[{order}]: node {row_id} is not a product")
            arguments = row.get("arguments")
            if not isinstance(arguments, list) or len(arguments) != 2:
                raise PortCodecError(f"{label}[{order}]: node {row_id} does not have exactly two factors")
            if any(isinstance(value, bool) or not isinstance(value, int) for value in arguments):
                raise PortCodecError(f"{label}[{order}]: node {row_id} contains a non-integer literal factor")
            pstar_positions = [index for index, value in enumerate(arguments) if value == pstar]
            if len(pstar_positions) != 1:
                raise PortCodecError(f"{label}[{order}]: node {row_id} does not contain exactly one Pstar factor")
            other_positions = [index for index in range(2) if index not in pstar_positions]
            other = arguments[other_positions[0]]
            self._row(other, label)
            if other != expected:
                raise PortCodecError(f"{label}[{order}]: node {row_id} underlying field {other!r}, expected {expected!r}")
            underlying.append(other)
            raw_rows.append({"row_id": row_id, "operation": row.get("operation"), "arguments": list(arguments)})
            self.audit["physical_rows_checked"] += 1
        receipt = {
            "kind": "physical",
            "label": label,
            "archive_decoded_sha256": self.archive_sha,
            "pstar": pstar,
            "row_ids": list(rows),
            "expected_fields": list(expected_fields),
            "underlying_field_ids": list(underlying),
            "raw_rows_in_input_order": raw_rows,
            "pstar_positions": [[index for index, value in enumerate(raw["arguments"]) if value == pstar][0] for raw in raw_rows],
            "underlying_positions": [[index for index, value in enumerate(raw["arguments"]) if value != pstar][0] for raw in raw_rows],
            "raw_row_hash_sha256": _canonical_hash(raw_rows),
        }
        self.physical_receipts.append(receipt)
        return underlying

    def compile(
        self,
        roots: Sequence[int],
        seeds: Mapping[int, Any],
        label: str,
        cap: int = DEFAULT_SLICE_CAP,
    ) -> list[Any]:
        """Compile one bounded closed arithmetic slice using a fresh cache."""
        if not isinstance(roots, (list, tuple)):
            raise PortCodecError(f"{label}: roots must be a sequence")
        if not isinstance(seeds, Mapping):
            raise PortCodecError(f"{label}: seeds must be a mapping")
        if isinstance(cap, bool) or not isinstance(cap, int) or cap <= 0:
            raise PortCodecError(f"{label}: cap must be a positive integer")
        root_ids = list(roots)
        seed_ids = list(seeds.keys())
        for node_id in root_ids + seed_ids:
            self._row(node_id, label)
        self.audit["compile_calls"] += 1
        self.audit["slice_rows_checked"] += len(root_ids) + len(seed_ids)

        # Breadth-first discovery stops exactly at caller-owned symbolic seeds.
        queue: deque[int] = deque(root_ids)
        discovered: set[int] = set()
        dependencies: dict[int, list[int]] = {}
        while queue:
            node_id = queue.popleft()
            if node_id in seeds or node_id in discovered:
                continue
            if len(discovered) >= cap:
                self.audit["caps_rejected"] += 1
                raise PortCodecError(f"{label}: closed slice exceeds cap {cap}")
            row = self._row(node_id, label)
            children = self._children(node_id, row, label)
            discovered.add(node_id)
            dependencies[node_id] = list(children)
            for child in children:
                if child not in seeds and child not in discovered:
                    queue.append(child)

        # Kahn's algorithm gives a deterministic topological order: the
        # smallest ready node is compiled first, while raw argument order is
        # retained in each operation's expression construction.
        nonseed = set(discovered)
        indegree = {node_id: len({child for child in dependencies[node_id] if child in nonseed}) for node_id in nonseed}
        dependents: dict[int, set[int]] = {node_id: set() for node_id in nonseed}
        for node_id, children in dependencies.items():
            for child in set(children):
                if child in nonseed:
                    dependents[child].add(node_id)
        ready = [node_id for node_id, degree in indegree.items() if degree == 0]
        heapq.heapify(ready)
        topological: list[int] = []
        while ready:
            node_id = heapq.heappop(ready)
            topological.append(node_id)
            for parent in sorted(dependents[node_id]):
                indegree[parent] -= 1
                if indegree[parent] == 0:
                    heapq.heappush(ready, parent)
        if len(topological) != len(nonseed):
            self.audit["cycles_rejected"] += 1
            raise PortCodecCycle(f"{label}: closed slice contains a cycle")

        cache = dict(seeds)
        opcounts = {operation: 0 for operation in ALLOWED_OPERATIONS}
        for node_id in topological:
            row = self.nodes[node_id]
            operation = row["operation"]
            opcounts[operation] += 1
            if operation == "exact_rational":
                value = self.sp.Rational(row["numerator"], row["denominator"])
            elif operation == "sum":
                value = self.sp.Add(*(cache[child] for child in row["arguments"]))
            elif operation == "product":
                value = self.sp.Mul(*(cache[child] for child in row["arguments"]))
            elif operation == "negative":
                value = -cache[row["argument"]]
            else:
                value = cache[row["numerator"]] / cache[row["denominator"]]
            cache[node_id] = value
            self.audit["actual_expansions"] += 1

        raw_rows = []
        for node_id in sorted(discovered):
            row = self.nodes[node_id]
            raw_rows.append({"row_id": node_id, **dict(row)})
        slice_receipt = {
            "kind": "compile",
            "label": label,
            "archive_decoded_sha256": self.archive_sha,
            "roots_in_input_order": root_ids,
            "seed_ids": sorted(seed_ids),
            "cap": cap,
            "allowed_operations": list(ALLOWED_OPERATIONS),
            "slice_ids_ascending": sorted(discovered),
            "topological_ids_ascending_ready_order": topological,
            "verified_dependencies_raw_order": {str(node_id): list(dependencies[node_id]) for node_id in sorted(dependencies)},
            "raw_rows_in_ascending_id_order": raw_rows,
            "raw_row_hash_sha256": _canonical_hash(raw_rows),
            "opcounts": opcounts,
            "actual_expansions": len(topological),
            "cache_entries": len(cache),
            "proof_normalization_performed": False,
            "profile_or_integral_evaluation_performed": False,
        }
        self.slices.append(slice_receipt)
        return [cache[node_id] for node_id in root_ids]


"""Fresh source-specific compact FTC port recipes.

This module is intentionally a caller-owned symbolic layer.  It accepts an
already authenticated native record, an already authenticated old record, a
symbolic module, and the separately supplied :class:`PortCodec`.  The codec
owns bounded graph slices; this file owns the source equations and C1 jet
checks.  It never imports a producer, a previous proof, a checker, a
normalizer, or SymPy, and it never evaluates a profile or an integral.

The universal FTC algebra is deliberately not repeated here.  The output is
the fresh source-row proof that the root composition may combine with its
cached formal FTC result.
"""

RECIPE_SCHEMA = "compact_ftc_port_recipes_v1"
NATIVE_RAW_SHA256 = "3593dd1ba3e510081090dd315d8581c1253afd5192233a284079c0646a556269"
NATIVE_DECODED_SHA256 = "b6a7227da5e2427561cd22fcff614003cc14f3182537e1dbe278a6de84322510"
OLD_RAW_SHA256 = "3921b728c3782d42802fef9e5a6f01bc244f65441811756d89dc98754aac0cce"
OLD_DECODED_SHA256 = "348df68bb426a2dfda6012e89f678b1e9c4ceb68b12f65f36853d5579977bc64"
SLICE_CAP = 1024

NATIVE_ROWS = {
    "U": [262269, 262270, 262271, 262272, 262273, 262274],
    "U_y": [262275, 262276, 262277, 262278, 262279, 262280],
    "M": [262281, 262282, 262283, 262284, 262285, 262286],
    "M_y": [262293, 262294, 262295, 262296, 262297, 262298],
    "Q": [262364, 262368, 262374],
    "Q_y": [262452, 262456],
    "S": [262603, 262604],
    "R_Q_y": [262611, 262612],
    "R_Q_y_x1": [263155, 263156],
    "R_Q_y_x2": [262914, 262915],
    "S1_weight": [262633],
    "S1_weighted_products_C1": [262634, 262635],
    "S1_pointwise_C1": [262638, 262639],
    "S1_integral_C1": [262642, 262643],
    "S1_full_history_C1": [262646, 262647],
}

OLD_ROWS = {
    "U": [30709, 30710, 30711, 30712],
    "U_y": [30713, 30714, 30715, 30716],
    "M": [30717, 30718, 30719, 30720],
    "M_y": [30721, 30722, 30723, 30724],
    "Q": [30758, 30762, 30768],
    "Q_y": [30802, 30806],
    "I1_q": [30848, 30850],
    "I1_q_Z": [30851, 30853],
    "I1_q_squared": [30855, 30857],
    "I1_U_q_Z": [30859, 30861],
    "I1_U_squared": [30863, 30865],
    "I1_U_U_Z": [30867, 30869],
    "S1": [31135, 31136],
    "m15": [31155, 31156],
}

NATIVE_FIELD_IDS = {
    "U": [260375, 260378, 260382, 260387, 260393, 260400],
    "U_y": [260494, 260497, 260501, 260506, 260512, 260519],
    "M": [261716, 261717, 261718, 261719, 261720, 261721],
}

OLD_FIELD_IDS = {
    "U": [13182, 13185, 13310, 24596],
    "U_y": [13231, 13232, 13644, 24837],
    "M": [13227, 13230, 13344, 24627],
    "M_y": [13247, 13248, 13658, 24850],
}

OLD_ENDPOINT_ROWS = {
    "lower": {
        "R_q": [30876, 30877],
        "R_squared_q": [30879, 30880],
        "R_q_squared": [30881, 30884],
        "R_squared_q_squared": [30885, 30888],
        "R_U_squared": [30891, 30894],
        "R_squared_U_squared": [30897, 30900],
        "q_plus_q_y": [30901, 30902],
        "R_q_y": [30903, 30904],
    },
    "upper": {
        "R_q": [30912, 30913],
        "R_squared_q": [30915, 30916],
        "R_q_squared": [30917, 30920],
        "R_squared_q_squared": [30921, 30924],
        "R_U_squared": [30927, 30930],
        "R_squared_U_squared": [30933, 30936],
        "q_plus_q_y": [30937, 30938],
        "R_q_y": [30939, 30940],
    },
}

OLD_I1_NAMES = (
    "I1_q",
    "I1_q_Z",
    "I1_q_squared",
    "I1_U_q_Z",
    "I1_U_squared",
    "I1_U_U_Z",
)


class RecipeFailure(ValueError):
    """A selected source row, recipe, or symbolic port is inconsistent."""


def _json_safe(value):
    """Convert receipts without importing a serializer or exposing objects."""
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    if isinstance(value, dict):
        return {str(key): _json_safe(item) for key, item in value.items()}
    return str(value)


def _same(actual, expected, label):
    if actual != expected:
        raise RecipeFailure(f"{label}: got {actual!r}, expected {expected!r}")


def _nodes(document, chart):
    if not isinstance(document, dict):
        raise RecipeFailure(f"{chart}: authenticated document must be a mapping")
    nodes = document.get("exact_graph_nodes")
    if not isinstance(nodes, list):
        raise RecipeFailure(f"{chart}: exact_graph_nodes is not a list")
    return nodes


def _row(nodes, node_id, label):
    if isinstance(node_id, bool) or not isinstance(node_id, int) or node_id < 0 or node_id >= len(nodes):
        raise RecipeFailure(f"{label}: invalid node id {node_id!r}")
    value = nodes[node_id]
    if not isinstance(value, dict):
        raise RecipeFailure(f"{label}: node {node_id} is not an object")
    return value


def _expect(nodes, node_id, label, operation, arguments=None, fields=None):
    value = _row(nodes, node_id, label)
    _same(value.get("operation"), operation, f"{label} node {node_id} operation")
    if arguments is not None:
        _same(value.get("arguments"), list(arguments), f"{label} node {node_id} arguments")
    if fields:
        for key, expected in fields.items():
            _same(value.get(key), expected, f"{label} node {node_id} field {key}")
    return value


def _expect_exact(nodes, node_id, label, numerator, denominator=1):
    return _expect(
        nodes,
        node_id,
        label,
        "exact_rational",
        fields={"numerator": numerator, "denominator": denominator},
    )


def _expect_physical(nodes, rows, pstar, fields, label):
    if len(rows) != len(fields):
        raise RecipeFailure(f"{label}: rows and fields have different lengths")
    for order, (node_id, field_id) in enumerate(zip(rows, fields)):
        value = _expect(nodes, node_id, f"{label}[{order}]", "product")
        arguments = value.get("arguments")
        if not isinstance(arguments, list) or len(arguments) != 2 or arguments.count(pstar) != 1:
            raise RecipeFailure(f"{label}[{order}] does not have exactly one Pstar factor")
        other = arguments[1] if arguments[0] == pstar else arguments[0]
        _same(other, field_id, f"{label}[{order}] primitive field")


def _derive_old_fields(nodes, rows, pstar, label):
    result = []
    for order, node_id in enumerate(rows):
        value = _expect(nodes, node_id, f"{label}[{order}]", "product")
        arguments = value.get("arguments")
        if not isinstance(arguments, list) or len(arguments) != 2 or arguments.count(pstar) != 1:
            raise RecipeFailure(f"{label}[{order}] does not have exactly one Pstar factor")
        result.append(arguments[1] if arguments[0] == pstar else arguments[0])
    return result


def _sym(sp, name, nonzero=False):
    try:
        return sp.Symbol(name, nonzero=nonzero)
    except TypeError:
        return sp.Symbol(name)


def _record(state, sp, label, actual, expected):
    """Record one actual symbolic residual with its complete text."""
    state["fresh_source_normalizations"] += 1
    state["residual_checks"] += 1
    residue = sp.cancel(sp.together(actual - expected))
    state["cancel_count"] += 1
    text = str(residue)
    try:
        passed = bool(residue == 0)
    except Exception:
        passed = text == "0"
    result = {
        "label": label,
        "passed": passed,
        "residual": "0" if passed else text,
        "residual_full": text,
    }
    state["partial_residuals"].append(result)
    return result


def _proof(state, label, residuals, **metadata):
    result = {
        "label": label,
        "passed": bool(residuals) and all(item.get("passed") is True for item in residuals),
        "residuals": list(residuals),
    }
    result.update(metadata)
    state["proofs"].append(result)
    return result


def _emit(state, emit, stage, payload):
    event = {"stage": stage, "payload": _json_safe(payload)}
    state["diagnostics"].append(event)
    if emit is not None:
        emit(stage, event["payload"])


def _D1(sp, expression, Z, jet_map):
    """Ordinary Z derivative for independent C1 value/dZ symbols."""
    result = sp.diff(expression, Z)
    for value, derivative in jet_map.items():
        result += sp.diff(expression, value) * derivative
    return result


def _source_symbols(sp, chart):
    Z = _sym(sp, f"{chart}_Z")
    delta = _sym(sp, f"{chart}_delta")
    Pstar = _sym(sp, f"{chart}_Pstar", nonzero=True)
    Rc = _sym(sp, f"{chart}_Rc")
    x = _sym(sp, f"{chart}_x")
    V = sp.Function(f"{chart}_V")(Z)
    Vy = sp.Function(f"{chart}_Vy")(Z)
    M = sp.Function(f"{chart}_m")(Z)
    My = sp.Function(f"{chart}_my")(Z)
    U = Pstar * V
    Uy = Pstar * Vy
    Mass = Pstar * M
    Mass_y = Pstar * My
    L = 1 - delta * Z**2
    d = 1 - Z**2
    h = 1 - delta
    return {
        "Z": Z,
        "delta": delta,
        "Pstar": Pstar,
        "Rc": Rc,
        "x": x,
        "R": Rc * x,
        "V": V,
        "Vy": Vy,
        "Mfun": M,
        "Myfun": My,
        "U": U,
        "Uy": Uy,
        "M": Mass,
        "My": Mass_y,
        "L": L,
        "d": d,
        "h": h,
    }


def _seed_physical(seeds, rows, fields, physical_values, raw_values):
    """Stop at a physical row with Pstar*raw, but expose raw leaves as raw."""
    for row_id, field_id, physical_value, raw_value in zip(rows, fields, physical_values, raw_values):
        seeds[row_id] = physical_value
        seeds[field_id] = raw_value


def _native_field_expectations(document):
    bindings = document.get("true_compact_C5_bindings")
    if not isinstance(bindings, dict):
        raise RecipeFailure("native: true_compact_C5_bindings is missing")
    physical = bindings.get("physical")
    profiles = bindings.get("profiles", {}).get("radial_y_rows")
    modified = bindings.get("modified_fields")
    if not isinstance(physical, dict) or not isinstance(profiles, list) or len(profiles) < 2 or not isinstance(modified, dict):
        raise RecipeFailure("native: primitive field metadata is incomplete")
    fields = {
        "U": list(profiles[0].get("V", [])),
        "U_y": list(profiles[1].get("V", [])),
        "M": list(modified.get("m", [])),
    }
    for name in ("U", "U_y", "M"):
        _same(fields[name], NATIVE_FIELD_IDS[name], f"native authenticated primitive fields {name}")
        physical_name = {"U": "U", "U_y": "U_y", "M": "m"}[name]
        _same(list(physical.get(physical_name, [])), NATIVE_ROWS[name], f"native authenticated physical rows {name} ({physical_name})")
    return fields


def _guard_native(document):
    nodes = _nodes(document, "native")
    fields = _native_field_expectations(document)
    _expect_exact(nodes, 1, "native", 1)
    _expect_exact(nodes, 48, "native", 2)
    _expect(nodes, 20, "native Z", "bound_variable", fields={"name": "Z"})
    _expect(nodes, 25220, "native delta", "analytic_unary", fields={"name": "exp", "argument": 25219})
    _expect(nodes, 251066, "native Rc", "function_substitution", fields={"expression": 250457, "variable": 250215, "value": 48, "Z_independent_substitution": True})
    _expect(nodes, 260181, "native x", "bound_variable", fields={"name": "true_compact_repair_x"})
    _expect(nodes, 260862, "native R", "product", arguments=[251066, 260181])
    _expect(nodes, 260183, "native 1/x", "positive_quotient", fields={"numerator": 1, "denominator": 260181})
    _expect(nodes, 25249, "native Z^2", "product", arguments=[20, 20])
    _expect(nodes, 25270, "native delta Z^2", "product", arguments=[25220, 25249])
    _expect(nodes, 25271, "native negative delta Z^2", "negative", fields={"argument": 25270})
    _expect(nodes, 25272, "native L", "sum", arguments=[1, 25271])
    _expect(nodes, 33572, "native Pstar", "analytic_unary", fields={"name": "exp", "argument": 5})
    for name in ("U", "U_y", "M"):
        _expect_physical(nodes, NATIVE_ROWS[name], 33572, fields[name], f"native {name}")
    for order, (row_id, u_id, m_id, neg_id) in enumerate(zip(NATIVE_ROWS["M_y"], NATIVE_ROWS["U"], NATIVE_ROWS["M"], range(262287, 262293))):
        _expect(nodes, neg_id, f"native M_y[{order}] negative", "negative", fields={"argument": m_id})
        _expect(nodes, row_id, f"native M_y[{order}] U-M", "sum", arguments=[u_id, neg_id])
    for row_id, numerator in zip(NATIVE_ROWS["Q"], (262359, 262367, 262373)):
        _expect(nodes, row_id, "native Q", "positive_quotient", fields={"numerator": numerator, "denominator": 25272})
    for row_id, numerator in zip(NATIVE_ROWS["Q_y"], (262447, 262455)):
        _expect(nodes, row_id, "native Q_y", "positive_quotient", fields={"numerator": numerator, "denominator": 25272})
    _expect(nodes, 262603, "native S0", "sum", arguments=[262496, 262521, 262586])
    _expect(nodes, 262604, "native S1", "sum", arguments=[262500, 262524, 262590])
    _expect(nodes, 262633, "native S1 weight", "product", arguments=[260183, 260862, 260862])
    _expect(nodes, 262634, "native S1 value product", "product", arguments=[262603, 262633])
    _expect(nodes, 262635, "native S1 dZ product", "product", arguments=[262604, 262633])
    _expect(nodes, 262615, "native S1 dummy", "bound_variable", fields={"name": "true_outer_pressure_true_compact"})
    for sub_id, point_id, integral_id, order in zip((262638, 262639), (262634, 262635), (262642, 262643), (0, 1)):
        _expect(nodes, sub_id, "native S1 substituted integrand", "function_substitution", fields={"expression": point_id, "variable": 260181, "value": 262615, "Z_independent_substitution": True})
        _expect(nodes, integral_id, "native S1 integral", "definite_integral", fields={"integrand": sub_id, "variable": "true_outer_pressure_true_compact", "lower": 1, "upper": 260181, "exact_function_integral": True, "ordinary_slow_Z_derivative_order": order, "true_core_inlet_complete_history": True, "numerical_cap_not_defining_integrand": True, "original_signed_reduced_pressure": True, "original_chart": "true_compact", "physical_radius_power": 2, "physical_Jacobian_applied_once": True, "endpoint_independent_of_Z": True, "nonzero_preceding_S0_S1_retained": True, "independent_P0_retained_in_original_source": True, "numerical_value_not_installed": True})
    _expect(nodes, 262646, "native full-history S1[0]", "sum", arguments=[258837, 262642])
    _expect(nodes, 262647, "native full-history S1[1]", "sum", arguments=[258838, 262643])
    for row_id, expression, value in ((263155, 262611, 1), (263156, 262612, 1), (262914, 262611, 48), (262915, 262612, 48)):
        _expect(nodes, row_id, "native RQy endpoint", "function_substitution", fields={"expression": expression, "variable": 260181, "value": value, "Z_independent_substitution": True})
    _expect(nodes, 262611, "native RQy[0]", "product", arguments=[262452, 260862])
    _expect(nodes, 262612, "native RQy[1]", "product", arguments=[262456, 260862])
    return {
        "nodes": nodes,
        "fields": fields,
        "rows": {key: list(value) for key, value in NATIVE_ROWS.items()},
        "local_S1_guarded": True,
        "local_S1_guard_rows": {
            "weight": list(NATIVE_ROWS["S1_weight"]),
            "pointwise_products": list(NATIVE_ROWS["S1_weighted_products_C1"]),
            "substitutions": list(NATIVE_ROWS["S1_pointwise_C1"]),
            "integrals": list(NATIVE_ROWS["S1_integral_C1"]),
            "full_history": list(NATIVE_ROWS["S1_full_history_C1"]),
            "dummy_bound_variable": 262615,
            "upper_coordinate": 260181,
            "arbitrary_upper_x_then_first_signed_graph_subx2": {
                "upper_coordinate_node": 260181,
                "first_pointwise_product": 262634,
                "first_substitution": 262638,
                "first_integral": 262642,
                "ordinary_Z_orders": [0, 1],
            },
        },
        "geometry_chain_guarded": True,
        "local_S1_measure_guarded": True,
        "incoming_full_history_retained": True,
    }


def _old_expected_fields(nodes):
    result = {}
    for name, rows in OLD_ROWS.items():
        if name in ("U", "U_y", "M", "M_y"):
            result[name] = _derive_old_fields(nodes, rows, 25203, f"old {name}")
            _same(result[name], OLD_FIELD_IDS[name], f"old authenticated primitive fields {name}")
    return result


def _guard_old(document):
    nodes = _nodes(document, "old")
    fields = _old_expected_fields(nodes)
    _expect_exact(nodes, 1, "old", 1)
    _expect_exact(nodes, 28, "old", 1, 2)
    _expect_exact(nodes, 48, "old", 2)
    _expect_exact(nodes, 220, "old", -1)
    _expect_exact(nodes, 639, "old", -2)
    _expect_exact(nodes, 13766, "old", 3)
    _expect(nodes, 20, "old Z", "bound_variable", fields={"name": "Z"})
    _expect(nodes, 25220, "old delta", "analytic_unary", fields={"name": "exp", "argument": 25219})
    _expect(nodes, 5556, "old x", "bound_variable", fields={"name": "repair_x"})
    _expect(nodes, 30828, "old Rc", "analytic_unary", fields={"name": "exp", "argument": 24956})
    _expect(nodes, 30829, "old R", "product", arguments=[30828, 5556])
    _expect(nodes, 25249, "old Z^2", "product", arguments=[20, 20])
    _expect(nodes, 25270, "old delta Z^2", "product", arguments=[25220, 25249])
    _expect(nodes, 25271, "old negative delta Z^2", "negative", fields={"argument": 25270})
    _expect(nodes, 25272, "old L", "sum", arguments=[1, 25271])
    _expect(nodes, 30377, "old source delta Z", "product", arguments=[25220, 20])
    _expect(nodes, 30378, "old source delta Z^2", "product", arguments=[30377, 20])
    _expect(nodes, 30380, "old source negative delta Z^2", "product", arguments=[220, 30378])
    _expect(nodes, 30382, "old source L", "sum", arguments=[1, 30380])
    _expect(nodes, 30409, "old source L squared", "product", arguments=[30382, 30382])
    _expect(nodes, 30538, "old source L fourth power", "product", arguments=[30409, 30409])
    _expect(nodes, 25203, "old Pstar", "current_original_physical_source_scale", fields={"name": "Pstar", "Z_independent": True, "physical_common_unit": True})
    for name in ("U", "U_y", "M", "M_y"):
        _expect_physical(nodes, OLD_ROWS[name], 25203, fields[name], f"old {name}")
    for row_id, numerator in zip(OLD_ROWS["Q"], (30755, 30761, 30767)):
        _expect(nodes, row_id, "old Q", "positive_quotient", fields={"numerator": numerator, "denominator": 25272})
    for row_id, numerator in zip(OLD_ROWS["Q_y"], (30799, 30805)):
        _expect(nodes, row_id, "old Q_y", "positive_quotient", fields={"numerator": numerator, "denominator": 25272})
    for name in OLD_I1_NAMES:
        rows = list(OLD_ROWS[name])
        for order, row_id in enumerate(rows):
            value = _expect(nodes, row_id, f"old {name}[{order}]", "definite_integral")
            _same(value.get("variable"), "repair_x", f"old {name}[{order}] variable")
            _same(value.get("lower"), 1, f"old {name}[{order}] lower")
            _same(value.get("upper"), 48, f"old {name}[{order}] upper")
            _same(value.get("ordinary_slow_Z_derivative_order"), order, f"old {name}[{order}] Z order")
            _same(value.get("measure"), "physical R*dR=Rc^2*x*dx", f"old {name}[{order}] measure")
            _same(value.get("numerical_value_not_installed"), True, f"old {name}[{order}] numeric guard")
    q_d_row = _row(nodes, OLD_ROWS["I1_q"][1], "old I1_q ordinary-Z row")
    q_z_row = _row(nodes, OLD_ROWS["I1_q_Z"][0], "old I1_q_Z value row")
    _same(q_d_row.get("integrand"), 30849, "old I1_q ordinary-Z integrand")
    _same(q_z_row.get("integrand"), 30849, "old I1_q_Z value integrand alias")
    for key in ("variable", "lower", "upper", "measure", "endpoint_independent_of_Z", "numerical_value_not_installed"):
        _same(q_d_row.get(key), q_z_row.get(key), f"old I1_q/I1_q_Z shared {key}")
    endpoint_shapes = {
        30876: ("product", [30828, 30872]), 30877: ("product", [30828, 30873]),
        30879: ("product", [30878, 30872]), 30880: ("product", [30878, 30873]),
        30881: ("product", [30876, 30872]), 30884: ("sum", [30882, 30883]),
        30885: ("product", [30879, 30872]), 30888: ("sum", [30886, 30887]),
        30891: ("product", [30889, 30870]), 30894: ("sum", [30892, 30893]),
        30897: ("product", [30895, 30870]), 30900: ("sum", [30898, 30899]),
        30901: ("sum", [30872, 30874]), 30902: ("sum", [30873, 30875]),
        30903: ("product", [30828, 30874]), 30904: ("product", [30828, 30875]),
        30912: ("product", [30911, 30907]), 30913: ("product", [30911, 30908]),
        30915: ("product", [30914, 30907]), 30916: ("product", [30914, 30908]),
        30917: ("product", [30912, 30907]), 30920: ("sum", [30918, 30919]),
        30921: ("product", [30915, 30907]), 30924: ("sum", [30922, 30923]),
        30925: ("product", [30911, 30905]), 30926: ("product", [30911, 30906]),
        30927: ("product", [30925, 30905]), 30928: ("product", [30926, 30905]),
        30929: ("product", [30925, 30906]), 30930: ("sum", [30928, 30929]),
        30931: ("product", [30914, 30905]), 30932: ("product", [30914, 30906]),
        30933: ("product", [30931, 30905]), 30934: ("product", [30932, 30905]),
        30935: ("product", [30931, 30906]), 30936: ("sum", [30934, 30935]),
        30937: ("sum", [30907, 30909]), 30938: ("sum", [30908, 30910]),
        30939: ("product", [30911, 30909]), 30940: ("product", [30911, 30910]),
    }
    for row_id, (operation, arguments) in endpoint_shapes.items():
        _expect(nodes, row_id, "old endpoint", operation, arguments=arguments)
    substitutions = {
        30870: (30709, 1), 30871: (30710, 1), 30872: (30758, 1), 30873: (30762, 1), 30874: (30802, 1), 30875: (30806, 1),
        30905: (30709, 48), 30906: (30710, 48), 30907: (30758, 48), 30908: (30762, 48), 30909: (30802, 48), 30910: (30806, 48),
    }
    for row_id, (expression, value) in substitutions.items():
        _expect(nodes, row_id, "old endpoint substitution", "function_substitution", fields={"expression": expression, "variable": 5556, "value": value, "Z_independent_substitution": True})
    _expect(nodes, 30911, "old upper radius", "product", arguments=[30828, 48])
    for row_id, operation, arguments in ((31133, "sum", [31131, 31107]), (31134, "sum", [31132, 31108]), (31135, "sum", [31133, 31119]), (31136, "sum", [31134, 31120]), (31147, "product", [28, 31135]), (31148, "product", [28, 31136]), (31149, "product", [220, 30903]), (31150, "product", [220, 30904]), (31151, "sum", [30939, 31149]), (31152, "sum", [30940, 31150]), (31153, "product", [220, 31151]), (31154, "product", [220, 31152]), (31155, "sum", [31147, 31153]), (31156, "sum", [31148, 31154])):
        _expect(nodes, row_id, "old S1/m15 parent", operation, arguments=arguments)
    return {
        "nodes": nodes,
        "fields": fields,
        "rows": {key: list(value) for key, value in OLD_ROWS.items()},
        "endpoints": {edge: {key: list(value) for key, value in group.items()} for edge, group in OLD_ENDPOINT_ROWS.items()},
        "q_Z_alias_guarded": True,
        "physical_I1_port_count": 6,
        "independent_C1_channel_count": 5,
    }


def _native_group(document, sp, Codec, state, emit):
    guard = _guard_native(document)
    nodes = guard["nodes"]
    symbols = _source_symbols(sp, "native")
    Z, delta, P = symbols["Z"], symbols["delta"], symbols["Pstar"]
    U, Uy, M = symbols["U"], symbols["Uy"], symbols["M"]
    L, d, h = symbols["L"], symbols["d"], symbols["h"]
    codec = Codec(nodes, NATIVE_DECODED_SHA256, "native", sp)
    state["codecs"]["native"] = codec
    geometry_actual = codec.compile([25272], {1: sp.Integer(1), 20: Z, 25220: delta}, "native L geometry chain", cap=SLICE_CAP)[0]
    geometry_proof = _proof(state, "native L=1-delta*Z^2 source geometry", [_record(state, sp, "native L geometry", geometry_actual, L)], rows=[25249, 25270, 25271, 25272], source_formula="1-delta*Z^2", geometry_chain_guarded=True)
    if not geometry_proof["passed"]:
        raise RecipeFailure("native L geometry proof failed before source seeds")
    seeds = {1: sp.Integer(1), 20: Z, 48: sp.Integer(2), 25220: delta, 25272: geometry_actual, 33572: P}
    U_raw_values = [sp.diff(symbols["V"], Z, order) for order in range(6)]
    Uy_raw_values = [sp.diff(symbols["Vy"], Z, order) for order in range(6)]
    M_raw_values = [sp.diff(symbols["Mfun"], Z, order) for order in range(6)]
    U_values = [P * value for value in U_raw_values]
    Uy_values = [P * value for value in Uy_raw_values]
    M_values = [P * value for value in M_raw_values]
    fields = guard["fields"]
    physical_receipts = {}
    for name, values in (("U", U_values), ("U_y", Uy_values), ("M", M_values)):
        physical_receipts[name] = codec.physical(NATIVE_ROWS[name], 33572, fields[name], f"native {name}")
        raw_values = {"U": U_raw_values, "U_y": Uy_raw_values, "M": M_raw_values}[name]
        _seed_physical(seeds, NATIVE_ROWS[name], fields[name], values, raw_values)
    my_actual = codec.compile(NATIVE_ROWS["M_y"][:3], seeds, "native M_y=U-M through C2", cap=SLICE_CAP)
    my_expected = [U - M, sp.diff(U - M, Z), sp.diff(U - M, Z, 2)]
    my_proof = _proof(state, "native M_y=U-M through Z2", [_record(state, sp, f"native M_y Z{order}", actual, expected) for order, (actual, expected) in enumerate(zip(my_actual, my_expected))], rows=NATIVE_ROWS["M_y"][:3])
    _emit(state, emit, "native_physical_ports", {"physical_receipts": physical_receipts, "M_y_rows": NATIVE_ROWS["M_y"][:3]})

    q_actual = codec.compile(NATIVE_ROWS["Q"], seeds, "native Q source quotient through C2", cap=SLICE_CAP)
    q_base = (2 * Z * U - h * Z * M - d * sp.diff(M, Z)) / L
    q_proof = _proof(state, "native Q=N(U,M)/L through Z2", [_record(state, sp, f"native Q Z{order}", actual, expected) for order, (actual, expected) in enumerate(zip(q_actual, (q_base, sp.diff(q_base, Z), sp.diff(q_base, Z, 2))))], rows=NATIVE_ROWS["Q"])
    qy_actual = codec.compile(NATIVE_ROWS["Q_y"], seeds, "native Q_y source quotient through C1", cap=SLICE_CAP)
    qy_base = (2 * Z * Uy - h * Z * (U - M) - d * sp.diff(U - M, Z)) / L
    qy_proof = _proof(state, "native Q_y=N(U_y,U-M)/L through C1", [_record(state, sp, f"native Q_y Z{order}", actual, expected) for order, (actual, expected) in enumerate(zip(qy_actual, (qy_base, sp.diff(qy_base, Z))))], rows=NATIVE_ROWS["Q_y"])
    rho_rhs = ((1 + delta) * Z * U + 2 * Z * Uy - d * sp.diff(U, Z)) / L
    rho_proof = _proof(state, "native rho and ordinary-Z derivative", [_record(state, sp, "native rho value", q_actual[0] + qy_actual[0], rho_rhs), _record(state, sp, "native rho ordinary-Z derivative", q_actual[1] + qy_actual[1], sp.diff(rho_rhs, Z))], rows=NATIVE_ROWS["Q"] + NATIVE_ROWS["Q_y"])
    s_actual = codec.compile(NATIVE_ROWS["S"], seeds, "native displayed S source rows", cap=SLICE_CAP)
    s_formula = (q_actual[0] + h * Z * q_actual[1] / 2 + qy_actual[0]) / L + q_actual[0] * (q_actual[0] / 2 + qy_actual[0]) + U * (-2 * Z * q_actual[0] + d * q_actual[1] - 2 * Z * qy_actual[0]) / L
    s_proof = _proof(state, "native displayed S formula and D(S)", [_record(state, sp, "native S value", s_actual[0], s_formula), _record(state, sp, "native S ordinary-Z derivative", s_actual[1], sp.diff(s_formula, Z))], rows=NATIVE_ROWS["S"])

    weight_seeds = {1: sp.Integer(1), 260183: 1 / symbols["x"], 260862: symbols["R"]}
    weight_actual = codec.compile(NATIVE_ROWS["S1_weight"], weight_seeds, "native S1 physical weight", cap=SLICE_CAP)[0]
    weight_proof = _proof(state, "native local S1 physical weight Rc^2*x", [_record(state, sp, "native S1 weight", weight_actual, symbols["R"] ** 2 / symbols["x"])], rows=NATIVE_ROWS["S1_weight"], integral_rows=NATIVE_ROWS["S1_integral_C1"], full_history_rows=NATIVE_ROWS["S1_full_history_C1"], integral_evaluations=0)
    weighted_seeds = {262633: weight_actual, 262603: s_actual[0], 262604: s_actual[1]}
    weighted_actual = codec.compile(NATIVE_ROWS["S1_weighted_products_C1"], weighted_seeds, "native weighted S1 pointwise C1 products", cap=SLICE_CAP)
    weighted_expected = [weight_actual * s_actual[0], weight_actual * s_actual[1]]
    weighted_proof = _proof(state, "native weighted S1 pointwise C1 products", [_record(state, sp, f"native weighted S1 Z{order}", actual, expected) for order, (actual, expected) in enumerate(zip(weighted_actual, weighted_expected))], rows=NATIVE_ROWS["S1_weighted_products_C1"], weight_row=262633, S_rows=NATIVE_ROWS["S"], arbitrary_upper_coordinate=260181, source_measure="R^2/x with physical Jacobian exactly once", integral_evaluations=0, full_history_rows=NATIVE_ROWS["S1_full_history_C1"])
    rqy_actual = codec.compile(NATIVE_ROWS["R_Q_y"], {**seeds, 260862: symbols["R"], 262452: qy_actual[0], 262456: qy_actual[1]}, "native RQy local boundary", cap=SLICE_CAP)
    rqy_proof = _proof(state, "native R Q_y one-radius boundary rows", [_record(state, sp, f"native RQy Z{order}", actual, symbols["R"] * expected) for order, (actual, expected) in enumerate(zip(rqy_actual, qy_actual))], rows=NATIVE_ROWS["R_Q_y"], endpoint_rows=NATIVE_ROWS["R_Q_y_x1"] + NATIVE_ROWS["R_Q_y_x2"])
    proofs = [geometry_proof, my_proof, q_proof, qy_proof, rho_proof, s_proof, weight_proof, weighted_proof, rqy_proof]
    result = {
        "chart": "native",
        "passed": all(item["passed"] for item in proofs),
        "proofs": proofs,
        "source_ports": guard["rows"],
        "physical_receipts": physical_receipts,
        "codec": codec,
        "native_local_S1_binding": bool(guard["local_S1_guarded"] and guard["local_S1_measure_guarded"] and guard["incoming_full_history_retained"] and weight_proof["passed"] and weighted_proof["passed"]),
        "native_C1_binding": all(item["passed"] for item in (my_proof, q_proof, qy_proof, rho_proof, s_proof, rqy_proof)),
    }
    _emit(state, emit, "native_source_proofs", {"passed": result["passed"], "proof_count": len(proofs), "residual_checks": state["residual_checks"]})
    return result


def _old_group(document, sp, Codec, state, emit):
    guard = _guard_old(document)
    nodes = guard["nodes"]
    symbols = _source_symbols(sp, "old")
    Z, delta, P = symbols["Z"], symbols["delta"], symbols["Pstar"]
    U, Uy, M, My = symbols["U"], symbols["Uy"], symbols["M"], symbols["My"]
    L, d, h = symbols["L"], symbols["d"], symbols["h"]
    codec = Codec(nodes, OLD_DECODED_SHA256, "old", sp)
    state["codecs"]["old"] = codec
    geometry_actual = codec.compile([25272, 30382, 30409, 30538], {1: sp.Integer(1), 20: Z, 25220: delta}, "old L source geometry chain", cap=SLICE_CAP)
    geometry_expected = (L, L, L**2, L**4)
    geometry_proof = _proof(state, "old L source geometry powers", [_record(state, sp, f"old geometry row {row_id}", actual, expected) for row_id, (actual, expected) in zip((25272, 30382, 30409, 30538), zip(geometry_actual, geometry_expected))], rows=[25272, 30382, 30409, 30538], source_formula="30382=L, 30409=L^2, 30538=L^4", logical_L_row=25272, source_L_row=30382, source_L_squared_row=30409, source_L_fourth_row=30538)
    if not geometry_proof["passed"]:
        raise RecipeFailure("old L source geometry proof failed before source seeds")
    old_L, old_source_L, old_source_L_squared, old_source_L_fourth = geometry_actual
    seeds = {1: sp.Integer(1), 20: Z, 28: sp.Rational(1, 2), 48: sp.Integer(2), 220: sp.Integer(-1), 25220: delta, 25203: P, 25272: old_L}
    U_raw_values = [sp.diff(symbols["V"], Z, order) for order in range(4)]
    Uy_raw_values = [sp.diff(symbols["Vy"], Z, order) for order in range(4)]
    M_raw_values = [sp.diff(symbols["Mfun"], Z, order) for order in range(4)]
    My_raw_values = [sp.diff(symbols["Myfun"], Z, order) for order in range(4)]
    U_values = [P * value for value in U_raw_values]
    Uy_values = [P * value for value in Uy_raw_values]
    M_values = [P * value for value in M_raw_values]
    My_values = [P * value for value in My_raw_values]
    physical_receipts = {}
    for name, values in (("U", U_values), ("U_y", Uy_values), ("M", M_values), ("M_y", My_values)):
        physical_receipts[name] = codec.physical(OLD_ROWS[name], 25203, guard["fields"][name], f"old {name}")
        raw_values = {"U": U_raw_values, "U_y": Uy_raw_values, "M": M_raw_values, "M_y": My_raw_values}[name]
        _seed_physical(seeds, OLD_ROWS[name], guard["fields"][name], values, raw_values)
    _emit(state, emit, "old_physical_ports", {"physical_receipts": physical_receipts})
    q_actual = codec.compile(OLD_ROWS["Q"], seeds, "old Q source quotient through C2", cap=SLICE_CAP)
    q_base = (2 * Z * U - h * Z * M - d * sp.diff(M, Z)) / L
    q_proof = _proof(state, "old Q=N(U,M)/L through Z2", [_record(state, sp, f"old Q Z{order}", actual, expected) for order, (actual, expected) in enumerate(zip(q_actual, (q_base, sp.diff(q_base, Z), sp.diff(q_base, Z, 2))))], rows=OLD_ROWS["Q"])
    qy_actual = codec.compile(OLD_ROWS["Q_y"], seeds, "old Q_y source quotient through C1", cap=SLICE_CAP)
    qy_base = (2 * Z * Uy - h * Z * My - d * sp.diff(My, Z)) / L
    qy_proof = _proof(state, "old Q_y=N(U_y,M_y)/L through C1", [_record(state, sp, f"old Q_y Z{order}", actual, expected) for order, (actual, expected) in enumerate(zip(qy_actual, (qy_base, sp.diff(qy_base, Z))))], rows=OLD_ROWS["Q_y"])

    x, Rc, R = symbols["x"], symbols["Rc"], symbols["R"]
    integrand_seeds = {1: sp.Integer(1), 20: Z, 48: sp.Integer(2), 25220: delta, 25203: P, 30828: Rc, 30829: R, 5556: x, 25272: old_L}
    integrand_seeds.update({30709: U_values[0], 30710: U_values[1], 30711: U_values[2], 30758: q_actual[0], 30762: q_actual[1], 30768: q_actual[2]})
    measure = Rc**2 * x
    integrand_targets = {
        "I1_q": measure * q_actual[0],
        "I1_q_Z": measure * q_actual[1],
        "I1_q_squared": measure * q_actual[0] ** 2,
        "I1_U_q_Z": measure * U_values[0] * q_actual[1],
        "I1_U_squared": measure * U_values[0] ** 2,
        "I1_U_U_Z": measure * U_values[0] * U_values[1],
    }
    integrand_proof_rows = {}
    integrand_residuals = []
    for name in OLD_I1_NAMES:
        integral_rows = OLD_ROWS[name]
        integrand_rows = [_row(nodes, row_id, f"old {name}").get("integrand") for row_id in integral_rows]
        integrand_proof_rows[name] = list(integrand_rows)
        actual = codec.compile(integrand_rows, integrand_seeds, f"old {name} integrand", cap=SLICE_CAP)
        expected = [integrand_targets[name], sp.diff(integrand_targets[name], Z)]
        for order, (a, e) in enumerate(zip(actual, expected)):
            integrand_residuals.append(_record(state, sp, f"old {name} integrand Z{order}", a, e))
    i1_proof = _proof(state, "old six I1 physical integrand pairs", integrand_residuals, rows=integrand_proof_rows, measure="Rc^2*x", integral_evaluations=0, q_Z_alias_guarded=guard["q_Z_alias_guarded"], independent_C1_channel_count=guard["independent_C1_channel_count"], physical_I1_port_count=guard["physical_I1_port_count"])
    _emit(state, emit, "old_i1_integrands", {"passed": i1_proof["passed"], "pairs": len(OLD_I1_NAMES), "residuals": len(integrand_residuals)})

    # The assembled old S1 is C1 in ordinary Z.  There are six physical I1
    # ports but five independent C1 channels: I1_q's dZ row is the same
    # integrand as the I1_q_Z value row, so the source alias is A -> B -> C.
    i1_values = {}
    jet_map = {}
    s1_seeds = {1: sp.Integer(1), 20: Z, 28: sp.Rational(1, 2), 48: sp.Integer(2), 220: sp.Integer(-1), 25220: delta, 30409: old_source_L_squared, 30538: old_source_L_fourth}
    iq_value = _sym(sp, "old_I1_q_value")
    iq_dz = _sym(sp, "old_I1_q_dZ")
    iqz_dz = _sym(sp, "old_I1_q_Z_dZ")
    i1_values["I1_q"] = iq_value
    i1_values["I1_q_Z"] = iq_dz
    jet_map[iq_value] = iq_dz
    jet_map[iq_dz] = iqz_dz
    s1_seeds[OLD_ROWS["I1_q"][0]] = iq_value
    s1_seeds[OLD_ROWS["I1_q"][1]] = iq_dz
    s1_seeds[OLD_ROWS["I1_q_Z"][0]] = iq_dz
    s1_seeds[OLD_ROWS["I1_q_Z"][1]] = iqz_dz
    for name in OLD_I1_NAMES:
        if name in ("I1_q", "I1_q_Z"):
            continue
        value = _sym(sp, f"old_{name}_value")
        derivative = _sym(sp, f"old_{name}_dZ")
        i1_values[name] = value
        jet_map[value] = derivative
        s1_seeds[OLD_ROWS[name][0]] = value
        s1_seeds[OLD_ROWS[name][1]] = derivative
    endpoint_values = {}
    endpoint_substitution_rows = {
        "lower": {"U": 30870, "U_Z": 30871, "Q": 30872, "Q_Z": 30873, "Q_y": 30874, "Q_y_Z": 30875},
        "upper": {"U": 30905, "U_Z": 30906, "Q": 30907, "Q_Z": 30908, "Q_y": 30909, "Q_y_Z": 30910},
    }
    for edge in ("lower", "upper"):
        endpoint_values[edge] = {}
        for field in ("U", "Q", "Q_y"):
            value = _sym(sp, f"old_{edge}_{field}")
            derivative = _sym(sp, f"old_{edge}_{field}_dZ")
            endpoint_values[edge][field] = value
            jet_map[value] = derivative
            sub_value = endpoint_substitution_rows[edge][field]
            sub_derivative = endpoint_substitution_rows[edge][field + "_Z"]
            s1_seeds[sub_value] = value
            s1_seeds[sub_derivative] = derivative
    endpoint_products = {}
    for edge in ("lower", "upper"):
        for name in ("R_squared_q", "R_squared_q_squared", "R_squared_U_squared", "R_q_y"):
            rows = OLD_ENDPOINT_ROWS[edge][name]
            if name == "R_q_y":
                base = _sym(sp, f"old_{edge}_R_q_y")
            else:
                base = _sym(sp, f"old_{edge}_{name}")
            derivative = _sym(sp, f"old_{edge}_{name}_dZ")
            endpoint_products[(edge, name)] = base
            jet_map[base] = derivative
            s1_seeds[rows[0]] = base
            s1_seeds[rows[1]] = derivative
    old_s1_actual = codec.compile(OLD_ROWS["S1"], s1_seeds, "old assembled S1 C1", cap=SLICE_CAP)
    boundary = lambda name: endpoint_products[("upper", name)] - endpoint_products[("lower", name)]
    s1_formula = boundary("R_squared_q") / L - i1_values["I1_q"] / L + h * Z * i1_values["I1_q_Z"] / (2 * L) + boundary("R_squared_q_squared") / 2 - i1_values["I1_q_squared"] / 2 + d * i1_values["I1_U_q_Z"] / L - 2 * Z**2 * boundary("R_squared_U_squared") / L**2 + 2 * h * Z**2 * i1_values["I1_U_squared"] / L**2 + 2 * Z * d * i1_values["I1_U_U_Z"] / L**2
    s1_proof = _proof(state, "old assembled S1 nine-term IBP formula and C1 derivative", [_record(state, sp, "old S1 value", old_s1_actual[0], s1_formula), _record(state, sp, "old S1 ordinary-Z derivative", old_s1_actual[1], _D1(sp, s1_formula, Z, jet_map))], rows=OLD_ROWS["S1"], old_pointwise_S_assumed=False, q_Z_alias_guarded=guard["q_Z_alias_guarded"], independent_C1_channel_count=guard["independent_C1_channel_count"], physical_I1_port_count=guard["physical_I1_port_count"])
    _emit(state, emit, "old_S1_ibp", {"passed": s1_proof["passed"], "residuals": len(s1_proof["residuals"])})

    endpoint_seeds = {1: sp.Integer(1), 20: Z, 28: sp.Rational(1, 2), 48: sp.Integer(2), 220: sp.Integer(-1), 25220: delta, 30828: Rc, 30911: 2 * Rc, 5556: x}
    endpoint_symbols = {}
    for edge in ("lower", "upper"):
        endpoint_symbols[edge] = {}
        for field in ("U", "Q", "Q_y"):
            value = _sym(sp, f"old_endpoint_{edge}_{field}_products")
            derivative = _sym(sp, f"old_endpoint_{edge}_{field}_products_dZ")
            endpoint_symbols[edge][field] = value
            endpoint_seeds[endpoint_substitution_rows[edge][field]] = value
            endpoint_seeds[endpoint_substitution_rows[edge][field + "_Z"]] = derivative
            jet_map[value] = derivative
    endpoint_residuals = []
    endpoint_product_names = ("R_squared_q", "R_squared_q_squared", "R_squared_U_squared")
    endpoint_product_rows = {edge: {} for edge in ("lower", "upper")}
    for edge in ("lower", "upper"):
        radius = Rc if edge == "lower" else 2 * Rc
        for name in endpoint_product_names:
            rows = OLD_ENDPOINT_ROWS[edge][name]
            endpoint_product_rows[edge][name] = list(rows)
            if name == "R_squared_q":
                target = radius**2 * endpoint_symbols[edge]["Q"]
            elif name == "R_squared_q_squared":
                target = radius**2 * endpoint_symbols[edge]["Q"]**2
            else:
                target = radius**2 * endpoint_symbols[edge]["U"]**2
            actual = codec.compile(rows, endpoint_seeds, f"old {edge} {name} endpoint", cap=SLICE_CAP)
            endpoint_residuals.append(_record(state, sp, f"old {edge} {name} endpoint value", actual[0], target))
            endpoint_residuals.append(_record(state, sp, f"old {edge} {name} endpoint dZ", actual[1], _D1(sp, target, Z, jet_map)))
    for edge in ("lower", "upper"):
        rows = OLD_ENDPOINT_ROWS[edge]["R_q_y"]
        target = (Rc if edge == "lower" else 2 * Rc) * endpoint_symbols[edge]["Q_y"]
        actual = codec.compile(rows, endpoint_seeds, f"old {edge} RQy endpoint", cap=SLICE_CAP)
        endpoint_residuals.append(_record(state, sp, f"old {edge} RQy endpoint value", actual[0], target))
        endpoint_residuals.append(_record(state, sp, f"old {edge} RQy endpoint dZ", actual[1], _D1(sp, target, Z, jet_map)))
    endpoint_proof = _proof(state, "old 16 finite endpoint products C1", endpoint_residuals, rows=endpoint_product_rows, R_Q_y_rows={edge: list(OLD_ENDPOINT_ROWS[edge]["R_q_y"]) for edge in ("lower", "upper")}, product_count=16)
    _emit(state, emit, "old_endpoint_products", {"passed": endpoint_proof["passed"], "product_count": 16, "residuals": len(endpoint_residuals)})

    m15_s1_value = _sym(sp, "old_m15_S1")
    m15_s1_dz = _sym(sp, "old_m15_S1_dZ")
    m15_lower_rqy_dz = _sym(sp, "old_lower_R_q_y_dZ")
    m15_upper_rqy_dz = _sym(sp, "old_upper_R_q_y_dZ")
    m15_seeds = {
        1: sp.Integer(1),
        28: sp.Rational(1, 2),
        220: sp.Integer(-1),
        31135: m15_s1_value,
        31136: m15_s1_dz,
        30903: endpoint_products[("lower", "R_q_y")],
        30904: m15_lower_rqy_dz,
        30939: endpoint_products[("upper", "R_q_y")],
        30940: m15_upper_rqy_dz,
    }
    # The C1 endpoint seeds used above are explicit independent pairs.  Keep
    # their derivative rows available without aliasing the value channel.
    m15_jet_map = {m15_s1_value: m15_s1_dz, endpoint_products[("lower", "R_q_y")]: m15_lower_rqy_dz, endpoint_products[("upper", "R_q_y")]: m15_upper_rqy_dz}
    m15_actual = codec.compile(OLD_ROWS["m15"], m15_seeds, "old m15 assembled C1", cap=SLICE_CAP)
    m15_formula = m15_s1_value / 2 - (endpoint_products[("upper", "R_q_y")] - endpoint_products[("lower", "R_q_y")])
    m15_proof = _proof(state, "old m15 half S1 minus RQy boundary", [_record(state, sp, "old m15 value", m15_actual[0], m15_formula), _record(state, sp, "old m15 ordinary-Z derivative", m15_actual[1], _D1(sp, m15_formula, Z, m15_jet_map))], rows=OLD_ROWS["m15"], lower_RQy=OLD_ENDPOINT_ROWS["lower"]["R_q_y"], upper_RQy=OLD_ENDPOINT_ROWS["upper"]["R_q_y"])
    _emit(state, emit, "old_m15", {"passed": m15_proof["passed"], "residuals": len(m15_proof["residuals"])})

    proofs = [geometry_proof, q_proof, qy_proof, i1_proof, s1_proof, endpoint_proof, m15_proof]
    result = {
        "chart": "old",
        "passed": all(item["passed"] for item in proofs),
        "proofs": proofs,
        "source_ports": guard["rows"],
        "physical_receipts": physical_receipts,
        "codec": codec,
        "old_C1_IBP_binding": bool(guard["q_Z_alias_guarded"] and s1_proof["passed"]),
        "endpoint_and_m15_C1_binding": endpoint_proof["passed"] and m15_proof["passed"],
        "old_mass_closure_proved": False,
        "old_pointwise_S_assumed": False,
        "q_Z_alias_guarded": guard["q_Z_alias_guarded"],
        "physical_I1_port_count": guard["physical_I1_port_count"],
        "independent_C1_channel_count": guard["independent_C1_channel_count"],
        "geometry_source_rows": {"logical_L": 25272, "source_L": 30382, "source_L_squared": 30409, "source_L_fourth": 30538},
    }
    _emit(state, emit, "old_source_proofs", {"passed": result["passed"], "proof_count": len(proofs), "residual_checks": state["residual_checks"]})
    return result


def _codec_receipt(codec):
    return {
        "audit": _json_safe(getattr(codec, "audit", {})),
        "physical_receipts": _json_safe(getattr(codec, "physical_receipts", [])),
        "slices": _json_safe(getattr(codec, "slices", [])),
    }


def evaluate_port_recipes(native_doc, old_doc, sp, Codec, emit=None):
    """Evaluate the fresh native/old source recipes once on supplied records.

    ``sp`` and ``Codec`` are dependency-injected so this source file has no
    import-time symbolic action.  Every selected slice is bounded by the
    codec's fixed roots and ``SLICE_CAP``.  A failed group keeps all proofs
    and residuals already produced, records the exception, and leaves the
    aggregate ``passed`` false.
    """
    state = {
        "fresh_source_normalizations": 0,
        "cancel_count": 0,
        "residual_checks": 0,
        "diagnostics": [],
        "proofs": [],
        "partial_residuals": [],
        "codecs": {},
    }
    groups = {}
    errors = []
    for chart, document, function in (("native", native_doc, _native_group), ("old", old_doc, _old_group)):
        try:
            groups[chart] = function(document, sp, Codec, state, emit)
        except Exception as error:
            failure = {"type": type(error).__name__, "message": str(error), "partial_proofs": list(state["proofs"]), "partial_residuals": list(state["partial_residuals"])}
            groups[chart] = {"chart": chart, "passed": False, "proofs": list(state["proofs"]), "failure": failure, "source_ports": {}, "physical_receipts": {}, "codec": state["codecs"].get(chart)}
            errors.append({"chart": chart, **failure})
            _emit(state, emit, f"{chart}_source_failure", failure)
    native_result = groups.get("native", {})
    old_result = groups.get("old", {})
    native_codec = native_result.get("codec")
    old_codec = old_result.get("codec")
    native_gate = bool(native_result.get("native_C1_binding"))
    local_gate = bool(native_result.get("native_local_S1_binding"))
    old_gate = bool(old_result.get("old_C1_IBP_binding"))
    endpoint_gate = bool(old_result.get("endpoint_and_m15_C1_binding"))
    passed = bool(native_result.get("passed") and old_result.get("passed") and native_gate and local_gate and old_gate and endpoint_gate and not errors)
    return {
        "schema": RECIPE_SCHEMA,
        "passed": passed,
        "gates": {
            "native_C1_binding": native_gate,
            "old_C1_IBP_binding": old_gate,
            "endpoint_and_m15_C1_binding": endpoint_gate,
            "native_local_S1_binding": local_gate,
        },
        "proofs": list(state["proofs"]),
        "source_ports": {"native": native_result.get("source_ports", {}), "old": old_result.get("source_ports", {})},
        "groups": {"native": {key: value for key, value in native_result.items() if key not in ("codec",)}, "old": {key: value for key, value in old_result.items() if key not in ("codec",)}},
        "codec_audits": {"native": _codec_receipt(native_codec) if native_codec is not None else {}, "old": _codec_receipt(old_codec) if old_codec is not None else {}},
        "physical_receipts": {"native": _json_safe(native_result.get("physical_receipts", {})), "old": _json_safe(old_result.get("physical_receipts", {}))},
        "slices": {"native": _json_safe(getattr(native_codec, "slices", [])) if native_codec is not None else [], "old": _json_safe(getattr(old_codec, "slices", [])) if old_codec is not None else []},
        "fresh_source_normalizations": state["fresh_source_normalizations"],
        "cancel_count": state["cancel_count"],
        "residual_checks": state["residual_checks"],
        "diagnostics": list(state["diagnostics"]),
        "failures": errors,
        "runtime_symbolic_version": getattr(sp, "__version__", None),
        "universal_ftc_reproved": False,
        "cached_universal_ftc_required": True,
        "old_mass_closure_proved": False,
        "old_pointwise_S_assumed": False,
        "physical_admission": False,
        "original_source_evaluations": 0,
        "profile_evaluations": 0,
        "integral_evaluations": 0,
        "quadrature_calls": 0,
        "recurrence_calls": 0,
    }


# Root composition for a new source-port compiler; cached universal proof only.
import base64
import traceback

PORT_BINDING_SCHEMA = "original_compact_ftc_port_binding"
PORT_BINDING_COMPONENT = "compact_ftc_port_binding"
PARENT_RECORD = "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_compact_ftc_transport.json"
SIGNED_RECORD = "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_signed_pullback_difference.json"
NATIVE_RECORD = "lei_ren_part1_paper_compliant_current_original_true_compact_C5_functions.json.gz"
OLD_RECORD = "lei_ren_part1_paper_compliant_current_original_n1_preRp_moment_functions.json.gz"
PORT_RECORD_HASHES = {
    PARENT_RECORD: ("c9857568a2e553c1963a94f670be4f6687a4b9c296dd2d10746b323aba64e033", "c9857568a2e553c1963a94f670be4f6687a4b9c296dd2d10746b323aba64e033"),
    SIGNED_RECORD: ("4a44b291afb4a0de8f82c652b0e99b4cd6e02eece8c86f44e2e94c76c5e4b9f7", "4a44b291afb4a0de8f82c652b0e99b4cd6e02eece8c86f44e2e94c76c5e4b9f7"),
    NATIVE_RECORD: ("3593dd1ba3e510081090dd315d8581c1253afd5192233a284079c0646a556269", "b6a7227da5e2427561cd22fcff614003cc14f3182537e1dbe278a6de84322510"),
    OLD_RECORD: ("3921b728c3782d42802fef9e5a6f01bc244f65441811756d89dc98754aac0cce", "348df68bb426a2dfda6012e89f678b1e9c4ceb68b12f65f36853d5579977bc64"),
}
CACHED_RESULT_SHA256 = "ff1242aac805000a11a38bc88d0a49c45a1b2695dfd49e991cea9962c927e6d5"
CACHED_FORMAL_SHA256 = "7737108e6b64c64825a36574c438b17b0e9d837c72ced9749b0e60003bd22690"
PORT_FALSE_FLAGS = (
    "native_old_density_identity", "native_old_density_closure", "source_ftc_zero_identity",
    "old212_adapter_identity", "old212_quantclosure", "old_errors_removed", "old_error_terms_removed",
    "quantitative_packet_error_reduction_proved", "numeric_values_computed", "numeric_containment_claimed",
    "source_C5_installation_claimed", "new_loop_domain_admission_claimed", "admissible_cone_admission_claimed",
    "global_five_moment_closure", "physical_layer_completed", "time_order_recursion_claimed", "full_goal_complete",
    "old_operators_replayed", "old_pointwise_S_recovered", "old_mass_relation_proved",
)


def _port_require(condition, message):
    if not condition:
        raise ValueError(message)


def _port_json_hash(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")).hexdigest()


def _port_authenticate(records):
    _port_require(isinstance(records, dict) and set(records) == set(PORT_RECORD_HASHES), "Four selected source records required")
    docs = {}
    evidence = {}
    for name, (raw_sha, decoded_sha) in PORT_RECORD_HASHES.items():
        item = records[name]
        _port_require(isinstance(item, dict), "Malformed record payload: " + name)
        _port_require(item.get("raw_sha256") == raw_sha and item.get("decoded_sha256") == decoded_sha and item.get("raw_authenticated_before_decode") is True, "Raw-first authentication evidence mismatch: " + name)
        _port_require(item.get("is_gzip") is name.endswith(".gz"), "Compression metadata mismatch: " + name)
        source = item.get("source_text")
        _port_require(isinstance(source, str) and hashlib.sha256(source.encode("utf-8")).hexdigest() == decoded_sha, "Decoded source hash mismatch: " + name)
        docs[name] = json.loads(source)
        evidence[name] = {"raw_sha256": raw_sha, "decoded_sha256": decoded_sha, "raw_authenticated_before_decode": True}
    return docs, evidence


def _port_cached_formal(parent):
    _port_require(parent.get("component") == "compact_ftc_transport" and parent.get("attempt_status") == "failed", "Cached parent must be the recorded failed source binding")
    _port_require(parent.get("source_ftc_transport_defined") is False and parent.get("all_16_native_B_accepted_from_parent") is True, "Parent scope changed")
    artifacts = [row for row in parent["lossless_artifacts"] if row.get("role") == "result_checkpoint"]
    _port_require(len(artifacts) == 1, "Cached raw result artifact missing")
    artifact = artifacts[0]
    saved = base64.b64decode(artifact["base64_bytes"], validate=True)
    _port_require(artifact["sha256"] == CACHED_RESULT_SHA256 and hashlib.sha256(saved).hexdigest() == CACHED_RESULT_SHA256, "Cached raw result hash changed")
    full = json.loads(saved.decode("utf-8"))
    _port_require(full == parent["component_result"] == parent["full_result"], "Cached result projection is not exact")
    formal = full["formal_algebra"]
    _port_require(_port_json_hash(formal) == CACHED_FORMAL_SHA256, "Cached formal object changed")
    _port_require(formal["passed"] is True and formal["runtime"]["sympy_version"] == "1.14.0", "Cached formal proof is not accepted")
    _port_require(formal["audit"]["new_symbolic_normalizations"]["sympy_cancel_together"] == 32, "Cached formal work count changed")
    _port_require(formal["formal_proof"]["formal_identity_passed"] is True and formal["source_jet_proof"]["source_function_identity_passed"] is True, "Cached arbitrary-coordinate proof missing")
    _port_require(formal["old_S1_IBP"]["old_C1_IBP_coefficient_identity_proved"] is True and formal["pullback_product_rule"]["passed"] is True, "Cached IBP or measure identity missing")
    negatives = formal["negative_controls"]["controls"]
    _port_require(len(negatives) == 9 and all(row["witness"]["is_zero"] is False for row in negatives.values()), "Cached negative controls changed")
    return {"passed": True, "parent_report_sha256": PORT_RECORD_HASHES[PARENT_RECORD][0], "raw_result_sha256": CACHED_RESULT_SHA256, "formal_object_sha256": CACHED_FORMAL_SHA256, "selector": "$.component_result.formal_algebra", "fresh_normalizations_in_cached_attempt": 32, "negative_controls_in_cached_attempt": 9, "cached_formal_proof_calls": 0, "rerun": False}


def _port_saved_signature(signed):
    _port_require(signed.get("component") == "signed_pullback_difference" and signed.get("accepted_control_stage") is True, "Signed function graph parent changed")
    view = signed["component_result"]["source_interface"]
    expected = {
        "native": {"U": [262269,262270,262271], "Q": [262364,262368,262374], "Q_y": [262452,262456], "pointwise_S": [262603,262604,262605,262606], "S1_integral_export": [262642,262643], "Rc_node": 251066, "R_node": 260862, "coordinate_node": 260181, "delta_node": 25220, "Z_node": 20, "endpoint_R_Q_y": {"x1":[263155,263156],"x2":[262914,262915]}},
        "old": {"U": [30709,30710,30711], "Q": [30758,30762,30768], "Q_y": [30802,30806], "S1_integral_export": [31135,31136], "Rc_node": 30828, "R_node": 30829, "coordinate_node": 5556, "delta_node": 25220, "Z_node": 20, "endpoint_R_q_y_already_includes_R": {"x1":[30903,30904],"x2":[30939,30940]}},
    }
    compared = []
    for chart, fields in expected.items():
        for key, value in fields.items():
            _port_require(view[chart].get(key) == value, "Saved signed source port changed: " + chart + "." + key)
            compared.append(chart + "." + key)
    _port_require(view["pullback"]["domain"] == ["1","2"], "Saved common integration interval changed")
    _port_require(view["source_formulas"]["transport"] == "integral_1^2 weighted_Delta_v dx + [Delta(R^2*F-R*Q_y)]_x1^x2", "Saved signed transport expression changed")
    return {"passed": True, "compared_ports": compared, "source_ports": expected, "source_report_sha256": PORT_RECORD_HASHES[SIGNED_RECORD][0], "whole_graph_replayed": False}


def _port_codec_controls(sp):
    # Synthetic typed rows only; these are not original profile values.
    nodes = [
        {"operation":"fixture"}, {"operation":"fixture_scale"}, {"operation":"fixture_field"},
        {"operation":"product","arguments":[2,1]},
        {"operation":"product","arguments":[1,2]},
        {"operation":"product","arguments":[1,1]},
        {"operation":"product","arguments":[2,2]},
        {"operation":"product","arguments":[1,8]},
        {"operation":"sum","arguments":[1,2]},
        {"operation":"product","arguments":[True,2]},
        {"operation":"product","arguments":[1,2,2]},
        {"operation":"not_permitted"},
        {"operation":"sum","arguments":[13,1]},
        {"operation":"product","arguments":[12,2]},
        {"operation":"sum","arguments":[1,2]},
        {"operation":"product","arguments":[14,2]},
        {"operation":"product","arguments":[1,999]},
    ]
    codec = PortCodec(nodes, "0"*64, "synthetic_fixture_only", sp)
    controls = []
    for name, node in [("field_first_accepted",3),("scale_first_accepted",4)]:
        value = codec.physical([node],1,[2],name)
        controls.append({"name":name,"passed":value==[2],"raw_arguments":nodes[node]["arguments"]})
    for name, node, expected in [("duplicate_scale_rejected",5,2),("missing_scale_rejected",6,2),("wrong_field_rejected",7,2),("wrong_operation_rejected",8,2),("boolean_factor_rejected",9,2),("three_factors_rejected",10,2),("out_of_range_field_rejected",16,999)]:
        try:
            codec.physical([node],1,[expected],name)
        except PortCodecError as error:
            controls.append({"name":name,"passed":True,"exception_type":type(error).__name__,"message":str(error)})
        else:
            controls.append({"name":name,"passed":False,"unexpected_acceptance":True})
    seeds = {1: sp.Symbol("fixture_scale"), 2: sp.Symbol("fixture_field")}
    for name, roots, cap, expected_error in [("unsupported_compile_rejected",[11],1024,PortCodecUnsupportedOperation),("cycle_rejected",[12],1024,PortCodecCycle),("slice_cap_rejected",[15],1,PortCodecError)]:
        try:
            codec.compile(roots,seeds,name,cap=cap)
        except expected_error as error:
            controls.append({"name":name,"passed":True,"exception_type":type(error).__name__,"message":str(error)})
        else:
            controls.append({"name":name,"passed":False,"unexpected_acceptance":True})
    return {"passed":all(row["passed"] for row in controls),"controls":controls,"count":len(controls),"audit":dict(codec.audit),"source_evaluations":0,"formal_proof_calls":0,"chart_scope":"synthetic fixtures only; production recipe codecs must use native/old fixed decoded archive hashes"}



def _port_recipe_link(recipe, signature):
    """Join newly compiled source ports to the previously saved graph ports."""
    _port_require(recipe.get("passed") is True, "Fresh source recipe did not pass")
    proof_rows = recipe.get("proofs", [])
    _port_require(bool(proof_rows) and all(row.get("passed") is True for row in proof_rows), "A source residual group is unproved")
    n = recipe.get("fresh_source_normalizations")
    _port_require(isinstance(n, int) and not isinstance(n, bool) and n > 0 and n == recipe.get("cancel_count") == recipe.get("residual_checks"), "Actual source normalization counters disagree")
    _port_require(n == sum(len(row["residuals"]) for row in proof_rows), "Source residual inventory differs from work counter")
    for chart, expected_sha, expected_pstar in [("native", NATIVE_ARCHIVE_SHA256, 33572),("old",OLD_ARCHIVE_SHA256,25203)]:
        receipt = recipe["codec_audits"][chart]
        slices = receipt["slices"]
        physical = receipt["physical_receipts"]
        _port_require(bool(slices) and bool(physical), "Missing production source context: " + chart)
        _port_require(all(row["archive_decoded_sha256"] == expected_sha and 0 < row["cap"] <= 1024 and len(row["slice_ids_ascending"]) <= row["cap"] for row in slices), "Production source chart or slice cap changed: " + chart)
        _port_require(all(row["archive_decoded_sha256"] == expected_sha and row["pstar"] == expected_pstar for row in physical), "Production physical scale namespace changed: " + chart)
        _port_require(receipt["audit"]["actual_expansions"] == sum(row["actual_expansions"] for row in slices), "Source expansion counter changed: " + chart)
    new = recipe["source_ports"]
    saved = signature["source_ports"]
    comparisons = []
    selected = {
        "native": [("U","U",3),("Q","Q",3),("Q_y","Q_y",2),("S","pointwise_S",2),("S1_integral_C1","S1_integral_export",2)],
        "old": [("U","U",3),("Q","Q",3),("Q_y","Q_y",2),("S1","S1_integral_export",2)],
    }
    for chart, ports in selected.items():
        for fresh_key, saved_key, count in ports:
            _port_require(new[chart][fresh_key][:count] == saved[chart][saved_key][:count], "Fresh/saved source port mismatch: " + chart + "." + fresh_key)
            comparisons.append(chart + "." + fresh_key)
    _port_require(new["native"]["R_Q_y_x1"] == saved["native"]["endpoint_R_Q_y"]["x1"] and new["native"]["R_Q_y_x2"] == saved["native"]["endpoint_R_Q_y"]["x2"], "Native endpoint ports changed")
    _port_require(new["old"]["m15"] == [31155,31156], "Old m15 source ports changed")
    return {"passed": True, "matched_ports": comparisons, "production_chart_names": ["native","old"], "source_proof_groups":len(proof_rows),"actual_source_normalizations":n,"cached_formal_proof_calls":0}


def derive(records, diagnostic_callback=None):
    """Compile current source ports once, using the saved universal proof."""
    def emit(stage, payload):
        if diagnostic_callback is not None:
            diagnostic_callback(stage, payload)

    result = {"schema":PORT_BINDING_SCHEMA,"component":PORT_BINDING_COMPONENT,"passed":False,"metadata_only":False,"source_only":True,"source_ftc_transport_defined":False,"source_bound_native_C1_FTC_proved":False,"source_bound_old_C1_IBP_proved":False,"source_bound_signed_transport_matches_saved_expression_C1_proved":False,"cached_formal_proof_calls":0,"original_source_evaluations":0,"original_numerical_field_evaluations":0,"source_refs_opaque":False,"primitive_profile_leaves_opaque":True,"source_record_dag_replay":False,"whole_archive_DAG_walks":0,"runtime_model_metadata":None,"source_recipe_result":None,"codec_controls":None,"cached_formal_algebra":None,"signed_graph_source_signature":None,"record_count":4}
    result.update({key:False for key in PORT_FALSE_FLAGS})
    result["formal_preconditions"] = {"conditional_source_identity":True,"native_radial_jet_contracts_required":["D_y U=U_y","D_y M=M_y","D_y Q=Q_y"],"mixed_derivatives_commute_required":True,"L_native_and_L_old_nonzero_required":True,"ordinary_Z_contract":"Z_Z=1; delta_Z=Rc_Z=R_Z=mu_Z=0","full_history_ODE_reproved_here":False,"physical_domain_admission":False}
    result["ordinary_Z_orders"]=[0,1]
    result["proof_coordinate_scope"]="arbitrary ordinary Z within the common source formula domain; no sampled-Z acceptance"
    result["owned_Z_cell"]=["0","1/2"]
    try:
        docs, authentication = _port_authenticate(records)
        result["authenticated_records"] = authentication
        emit("authenticated_four_records", {"count":4,"raw_first":True})
        parent = docs[PARENT_RECORD]
        result["cached_formal_algebra"] = _port_cached_formal(parent)
        emit("cached_formal_result_read",result["cached_formal_algebra"])
        result["signed_graph_source_signature"] = _port_saved_signature(docs[SIGNED_RECORD])
        for key in ["retained_construction_data","retained_integral_term_inventory"]:
            result[key] = json.loads(json.dumps(parent[key]))
        ledger = result["retained_construction_data"]
        _port_require(ledger["cumulative_terms"]==212 and ledger["raw_signed_terms"]==213 and ledger["reduced_terms"]==183 and ledger["reserved_contribution_registered_once"] is True and ledger["all_old_error_data_retained"] is True,"Retained signed ledger changed")
        result["accepted_native_B_nodes"] = list(range(520,536))
        result["remaining_native_B_nodes"] = []
        result["all_16_native_B_accepted_from_parent"] = True
        import sympy as sp
        result["runtime"]={"sympy_version":sp.__version__,"recorded_after_entry":True}
        result["codec_controls"]=_port_codec_controls(sp)
        emit("new_typed_codec_controls",result["codec_controls"])
        recipe = evaluate_port_recipes(docs[NATIVE_RECORD], docs[OLD_RECORD], sp, PortCodec, emit)
        result["source_recipe_result"] = recipe
        link = _port_recipe_link(recipe, result["signed_graph_source_signature"])
        result["source_recipe_graph_link"] = link
        gates=recipe.get("gates",{})
        native=bool(recipe.get("passed") is True and gates.get("native_C1_binding") is True and gates.get("native_local_S1_binding") is True and result["codec_controls"]["passed"] and result["cached_formal_algebra"]["passed"])
        old=bool(recipe.get("passed") is True and gates.get("old_C1_IBP_binding") is True and gates.get("endpoint_and_m15_C1_binding") is True and result["codec_controls"]["passed"] and result["cached_formal_algebra"]["passed"])
        matched=bool(native and old and result["codec_controls"]["passed"] and result["cached_formal_algebra"]["passed"] and result["signed_graph_source_signature"]["passed"])
        result["source_bound_native_C1_FTC_proved"]=native
        result["source_bound_old_C1_IBP_proved"]=old
        result["source_bound_signed_transport_matches_saved_expression_C1_proved"]=matched
        result["source_ftc_transport_defined"]=matched
        result["passed"]=matched
        result["arbitrary_Z_functional_identity"]=matched
    except Exception as error:
        result["failure"]={"type":type(error).__name__,"message":str(error),"traceback":traceback.format_exc()}
        emit("composition_exception",result["failure"])
    recipe=result.get("source_recipe_result") or {}
    result["fresh_source_normalizations"]=recipe.get("fresh_source_normalizations")
    result["codec_audits"]=recipe.get("codec_audits",{})
    result["source_slices"]=recipe.get("slices",{})
    result["physical_port_receipts"]=recipe.get("physical_receipts",{})
    result["source_inspection_audit"]={"cached_formal_proof_calls":0,"old_source_consumer_calls":0,"old_proof_functions_called":0,"old_operator_calls":0,"old_source_normalizer_calls":0,"original_source_evaluations":0,"original_numerical_field_evaluations":0,"profile_evaluations":0,"integral_evaluations":0,"quadratures":0,"recurrence_calls":0,"whole_archive_DAG_walks":0,"fresh_source_normalizations":recipe.get("fresh_source_normalizations"),"codec_audits":recipe.get("codec_audits",{}),"source_refs_opaque":False,"primitive_profile_leaves_opaque":True,"runtime_model_metadata":None}
    emit("new_port_binding_complete",{"passed":result["passed"],"source_ftc_transport_defined":result["source_ftc_transport_defined"],"fresh_source_normalizations":result["fresh_source_normalizations"],"cached_formal_proof_calls":0})
    return result

__all__ = ["derive", "PortCodec", "evaluate_port_recipes"]
