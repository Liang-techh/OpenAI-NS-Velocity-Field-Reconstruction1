"""Bounded reachability over the already-saved F gap fresh DAG.

This module is deliberately a consumer of the published fresh-DAG receipt.  It
does not construct a graph, import a project constructor, read the old source
graph, or follow a saved/external namespace.  The only edges considered are
operand fields from the frozen source-row schema.  A saved reference to node
``33982`` is a terminal opaque boundary, not a node whose operator row may be
read.

The result is structural reachability evidence only.  In particular, a hit is
not an F-radius owner proof, a function-identity proof, a capture-free proof,
an ordinary-Z transport proof, or an m15/physical equality proof.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any, Iterable, Mapping


SCHEMA = "f_gap_saved_radius_reachability_v1"
COMPONENT = SCHEMA

PRIOR_RESULT = Path("tmp/frame_gap_fresh_substitution_v1/result.json")
PRIOR_RESULT_SHA256 = (
    "a37fc375ff06c06d39b4250ffa3d36c6c315182a81e5d783f4bfb1f981ba0234"
)
PRIOR_AUDIT_RESULT = Path("tmp/frame_gap_fresh_substitution_audit_v1/result.json")
PRIOR_AUDIT_RESULT_SHA256 = (
    "a3de519abfdc42fe18a0fa6830105c9c2c39a7bb4188d27693123117123087a3"
)
PRIOR_AUDIT_SCHEMA = "frame_gap_fresh_substitution_audit_v1"
PRIOR_SCHEMA = "frame_gap_fresh_substitution_v1"
PRIOR_PACKET_STATUS = "returned"

FRESH_NAMESPACE = "fresh_frame_gap_substitution_v1"
SAVED_NAMESPACE = "saved_Rp_native_C5_pulse_source"
TARGET_SAVED_NODE = 33982
EXPECTED_FRESH_NODE_COUNT = 1131

# The previous constructor's operand schema is frozen locally.  No source-row
# metadata key is an edge unless it appears in this table.
FROZEN_OPERAND_SCHEMA: dict[str, tuple[str, ...]] = {
    "mathematical_positive_infinity": (),
    "exact_original_source_heat_phi": ("argument",),
    "exact_rational": (),
    "bound_variable": (),
    "mathematical_pi": (),
    "current_original_source_parameter": (),
    "original_source_parameter": (),
    "exact_positive_integer_expression": (),
    "shared_positive_integer": (),
    "sum": ("arguments",),
    "product": ("arguments",),
    "negative": ("argument",),
    "analytic_unary": ("argument",),
    "positive_quotient": ("numerator", "denominator"),
    "exact_real_comparison": ("left", "right"),
    "logical_or": ("arguments",),
    "compact_raw_beta": ("argument",),
    "exact_original_source_flat_sigma": ("argument",),
    "current_original_leading_function_recipe": ("coordinate", "Z_variable"),
    "current_original_physical_source_scale": (),
    "current_original_physical_source_quantity": ("coordinate", "Z_variable"),
    "current_original_physical_radius": ("log_radius", "free_query"),
    "substitute_original_native_coordinate": ("body", "query_coordinate"),
    "substitute_free_query_at_native_coordinate": ("body", "free_query"),
    "substitute_original_inverse_angle": ("body", "inverse_angle"),
    "function_substitution": ("expression", "variable", "value"),
    "definite_integral": ("integrand", "lower", "upper"),
    "original_lazy_flat_branch": (
        "active_body",
        "Delta",
        "eta",
        "flat_predicate",
        "original_collar_coordinate",
        "original_collar_upper",
    ),
    "original_monotone_phase_inverse": (
        "phase_function",
        "target_phase",
        "angle_lower",
        "angle_upper",
        "Delta",
        "eta",
    ),
}

# The previous fresh receipt has exactly these roots.  The field-radius root is
# intentionally excluded: it is an authenticated boundary record, not a path
# query root.
ROOT_SPEC: tuple[tuple[str, str, int], ...] = (
    ("gap_m15", "gap_m15", 0),
    ("gap_m15", "gap_m15", 1),
    ("gap_m15", "gap_m15", 2),
    ("gap_m15", "gap_m15", 3),
    ("lower", "R_q_y", 0),
    ("lower", "R_q_y", 1),
    ("lower", "R_q_y", 2),
    ("lower", "R_q_y", 3),
    ("lower", "q_plus_q_y", 0),
    ("lower", "q_plus_q_y", 1),
    ("lower", "q_plus_q_y", 2),
    ("lower", "q_plus_q_y", 3),
    ("upper", "R_q_y", 0),
    ("upper", "R_q_y", 1),
    ("upper", "R_q_y", 2),
    ("upper", "R_q_y", 3),
    ("upper", "q_plus_q_y", 0),
    ("upper", "q_plus_q_y", 1),
    ("upper", "q_plus_q_y", 2),
    ("upper", "q_plus_q_y", 3),
)

PROTECTED_NODE_IDS = frozenset(
    {
        0,
        1,
        17,
        33572,
        25220,
        263344,
        263348,
        263349,
        262918,
        205147,
        263346,
    }
)

REQUIRED_INPUTS = {
    PRIOR_RESULT.as_posix(): PRIOR_RESULT_SHA256,
    PRIOR_AUDIT_RESULT.as_posix(): PRIOR_AUDIT_RESULT_SHA256,
}


class FGapSavedRadiusReachabilityInputError(ValueError):
    """A frozen manifest, receipt, or bounded DAG scope is malformed."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise FGapSavedRadiusReachabilityInputError(message)


def _sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def _path(value: Any) -> str:
    _require(isinstance(value, str) and value, "manifest path missing")
    return value.replace("\\", "/")


def _digest(value: Any, label: str) -> str:
    _require(isinstance(value, str) and len(value) == 64, label + " SHA malformed")
    value = value.lower()
    _require(
        all(char in "0123456789abcdef" for char in value),
        label + " SHA malformed",
    )
    return value


def _opaque_node_id(value: Any) -> bool:
    """Accept only a non-bool integer as a saved/source node handle."""

    return type(value) is int and value >= 0


def _safe(root: Path, name: str) -> Path:
    relative = Path(name)
    _require(
        not relative.is_absolute() and ".." not in relative.parts,
        "manifest path escapes root",
    )
    result = (root / relative).resolve()
    _require(result.is_relative_to(root.resolve()), "manifest path escapes root")
    return result


def _manifest_hashes(
    source: Mapping[str, Any], hashes: Mapping[str, str] | None
) -> dict[str, str]:
    entries: dict[str, str] = {}
    files = source.get("files")
    if isinstance(files, list):
        for index, item in enumerate(files):
            _require(isinstance(item, Mapping), "manifest files entry malformed")
            name = _path(item.get("path"))
            value = _digest(item.get("sha256"), f"manifest.files[{index}]")
            _require(name not in entries, "manifest duplicate path: " + name)
            entries[name] = value
    for field in ("sha256", "input_hashes", "hashes"):
        table = source.get(field)
        if isinstance(table, Mapping):
            for name, value in table.items():
                key = _path(str(name))
                digest = _digest(value, "manifest." + field + "." + key)
                _require(
                    key not in entries or entries[key] == digest,
                    "manifest hash conflict: " + key,
                )
                entries[key] = digest
    if isinstance(hashes, Mapping):
        for name, value in hashes.items():
            key = _path(str(name))
            digest = _digest(value, "caller hashes." + key)
            _require(
                key not in entries or entries[key] == digest,
                "caller hash conflict: " + key,
            )
            entries[key] = digest
    _require(entries, "caller frozen manifest has no hashes")
    return entries


def _find_hash(entries: Mapping[str, str], name: str) -> str:
    if name in entries:
        return entries[name]
    normalized = name.replace("\\", "/")
    candidates = [
        value
        for key, value in entries.items()
        if key.replace("\\", "/").endswith("/" + normalized)
    ]
    _require(len(candidates) == 1, "required manifest hash missing: " + name)
    return candidates[0]


def _authenticate(
    root: Path, source: Mapping[str, Any], hashes: Mapping[str, str] | None
) -> tuple[dict[str, bytes], dict[str, str], int]:
    """Hash every manifest file before any prior-result JSON is decoded."""

    entries = _manifest_hashes(source, hashes)
    for name, expected in REQUIRED_INPUTS.items():
        _require(
            _find_hash(entries, name) == expected,
            "required SHA changed: " + name,
        )

    raw_by_name: dict[str, bytes] = {}
    actual: dict[str, str] = {}
    # This loop is intentionally raw-byte-only.  JSON decoding occurs later.
    for raw_name, expected in entries.items():
        name = raw_name.replace("\\", "/")
        if name.startswith("/") or ":/" in name:
            aliases = [item for item in REQUIRED_INPUTS if name.endswith("/" + item)]
            _require(len(aliases) == 1, "absolute manifest alias is not required input")
            name = aliases[0]
        data = _safe(root, name).read_bytes()
        observed = _sha(data)
        _require(observed == expected, "raw SHA changed: " + name)
        raw_by_name[name] = data
        actual[name] = observed

    # A required input may have been supplied only through a suffix alias or a
    # manifest with a noncanonical path.  Still authenticate it before parsing.
    for name, expected in REQUIRED_INPUTS.items():
        if name not in raw_by_name:
            data = _safe(root, name).read_bytes()
            observed = _sha(data)
            _require(observed == expected, "required raw SHA changed: " + name)
            raw_by_name[name] = data
            actual[name] = observed
    return raw_by_name, actual, len(entries)


def _decode_prior(raw: bytes) -> Mapping[str, Any]:
    try:
        packet = json.loads(raw.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise FGapSavedRadiusReachabilityInputError(
            "prior fresh-DAG result is not valid JSON"
        ) from exc
    _require(isinstance(packet, Mapping), "prior fresh-DAG packet is not an object")
    _require(packet.get("status") == PRIOR_PACKET_STATUS, "prior packet status changed")
    body = packet.get("result")
    _require(isinstance(body, Mapping), "prior packet result is not an object")
    _require(body.get("schema") == PRIOR_SCHEMA, "prior result schema changed")
    return body


def _decode_audit(raw: bytes) -> Mapping[str, Any]:
    try:
        packet = json.loads(raw.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise FGapSavedRadiusReachabilityInputError(
            "prior fresh-DAG audit result is not valid JSON"
        ) from exc
    _require(isinstance(packet, Mapping), "prior audit packet is not an object")
    _require(packet.get("status") == PRIOR_PACKET_STATUS, "prior audit packet status changed")
    body = packet.get("result")
    _require(isinstance(body, Mapping), "prior audit result is not an object")
    _require(body.get("schema") == PRIOR_AUDIT_SCHEMA, "prior audit schema changed")
    _require(body.get("status") == "passed", "prior fresh-DAG audit did not pass")
    # The independently-authenticated audit must refer to this exact fresh
    # result.  The receipt exposes both the observed and expected digest names;
    # accepting neither is safer than accepting an unrelated passed audit.
    _require(
        body.get("result_sha256") == PRIOR_RESULT_SHA256
        or body.get("expected_result_sha256") == PRIOR_RESULT_SHA256,
        "prior audit does not correspond to the frozen fresh result",
    )
    return body


def _mapping(value: Any, label: str) -> Mapping[str, Any]:
    _require(isinstance(value, Mapping), label + " must be an object")
    return value


def _fresh_ref(value: Any, label: str) -> dict[str, Any]:
    ref = _mapping(value, label)
    _require(ref.get("namespace") == FRESH_NAMESPACE, label + " is not fresh")
    node = ref.get("node")
    _require(isinstance(node, str) and node, label + " fresh node is malformed")
    return {"namespace": FRESH_NAMESPACE, "node": node}


def _saved_ref(value: Any, label: str) -> tuple[str, Any]:
    ref = _mapping(value, label)
    namespace = ref.get("namespace")
    _require(isinstance(namespace, str) and namespace, label + " namespace missing")
    node = ref.get("node")
    if namespace == SAVED_NAMESPACE:
        _require(_opaque_node_id(node), label + " saved node is malformed")
    return namespace, node


def _root_entries(body: Mapping[str, Any]) -> list[dict[str, Any]]:
    roots = _mapping(body.get("roots"), "prior roots")
    _require(set(roots) >= {"gap_m15", "endpoints", "field_radius"}, "prior roots are incomplete")
    # Read the radius boundary only as metadata.  It is not appended to the
    # root list and its target is never dereferenced.
    boundary = _mapping(roots.get("field_radius"), "field_radius boundary")
    _require(
        _opaque_node_id(boundary.get("original_ref"))
        and boundary.get("original_ref") == TARGET_SAVED_NODE,
        "radius boundary original_ref changed",
    )
    boundary_ref = _mapping(boundary.get("node_ref"), "field_radius boundary node_ref")
    _require(boundary_ref.get("namespace") == SAVED_NAMESPACE, "radius boundary namespace changed")
    _require(
        _opaque_node_id(boundary_ref.get("node"))
        and boundary_ref.get("node") == TARGET_SAVED_NODE,
        "radius boundary node changed",
    )
    radius_boundary = body.get("radius_boundary")
    radius_boundary = _mapping(radius_boundary, "radius_boundary")
    _require(
        _opaque_node_id(radius_boundary.get("frame_ref"))
        and radius_boundary.get("frame_ref") == TARGET_SAVED_NODE,
        "radius boundary frame_ref changed",
    )
    _require(
        _opaque_node_id(radius_boundary.get("original_ref"))
        and radius_boundary.get("original_ref") == TARGET_SAVED_NODE,
        "radius boundary original_ref metadata changed",
    )
    _require(radius_boundary.get("same_opaque_id") is True, "radius boundary identity flag changed")

    output: list[dict[str, Any]] = []
    gap_rows = roots.get("gap_m15")
    _require(isinstance(gap_rows, list) and len(gap_rows) == 4, "gap_m15 roots are not four rows")
    for order, entry in enumerate(gap_rows):
        item = _mapping(entry, "gap_m15 root")
        _require(item.get("row_order") == order, "gap_m15 root order changed")
        output.append(
            {
                "role": "gap_m15",
                "field": "gap_m15",
                "row_order": order,
                "root_ref": item.get("node_ref"),
                "original_ref": item.get("original_ref"),
            }
        )

    endpoints = _mapping(roots.get("endpoints"), "endpoint roots")
    for side in ("lower", "upper"):
        side_rows = _mapping(endpoints.get(side), side + " endpoint roots")
        for field in ("R_q_y", "q_plus_q_y"):
            rows = side_rows.get(field)
            _require(
                isinstance(rows, list) and len(rows) == 4,
                side + "." + field + " roots are not four rows",
            )
            for order, entry in enumerate(rows):
                item = _mapping(entry, side + "." + field + " root")
                _require(
                    item.get("row_order") == order,
                    side + "." + field + " root order changed",
                )
                output.append(
                    {
                        "role": side + "." + field,
                        "field": field,
                        "side": side,
                        "row_order": order,
                        "root_ref": item.get("node_ref"),
                        "original_ref": item.get("original_ref"),
                    }
                )
    _require(len(output) == len(ROOT_SPEC), "bounded root scope is not exactly twenty roots")
    actual_scope = [
        (item["role"], item["field"], item["row_order"]) for item in output
    ]
    _require(actual_scope == list(ROOT_SPEC), "bounded root role/order scope changed")
    return output


def _index_nodes(body: Mapping[str, Any]) -> tuple[dict[str, Mapping[str, Any]], list[dict[str, Any]]]:
    nodes = body.get("nodes")
    _require(isinstance(nodes, list), "prior fresh nodes are not a list")
    _require(len(nodes) == EXPECTED_FRESH_NODE_COUNT, "prior fresh node count changed")
    _require(
        body.get("fresh_node_count") == EXPECTED_FRESH_NODE_COUNT,
        "declared fresh node count changed",
    )
    by_id: dict[str, Mapping[str, Any]] = {}
    malformed: list[dict[str, Any]] = []
    for position, raw_node in enumerate(nodes):
        if not isinstance(raw_node, Mapping):
            malformed.append({"code": "node_not_object", "position": position})
            continue
        node_id = raw_node.get("node_id")
        if not isinstance(node_id, str) or not node_id:
            malformed.append({"code": "node_id_malformed", "position": position})
            continue
        if node_id in by_id:
            malformed.append({"code": "node_id_duplicate", "node": node_id})
            continue
        by_id[node_id] = raw_node
        if raw_node.get("namespace") != FRESH_NAMESPACE:
            malformed.append({"code": "node_namespace_mismatch", "node": node_id})
    return by_id, malformed


def _operand_edges(
    node: Mapping[str, Any],
) -> tuple[list[tuple[str, Any]], list[dict[str, Any]]]:
    operator_row = node.get("operator_row")
    if not isinstance(operator_row, Mapping):
        return [], [{"code": "operator_row_missing"}]
    operation = operator_row.get("operation")
    if not isinstance(operation, str) or operation not in FROZEN_OPERAND_SCHEMA:
        return [], [{"code": "operation_not_frozen", "operation": operation}]
    if node.get("operand_schema") != list(FROZEN_OPERAND_SCHEMA[operation]):
        return [], [
            {
                "code": "operand_schema_mismatch",
                "operation": operation,
                "expected": list(FROZEN_OPERAND_SCHEMA[operation]),
                "actual": node.get("operand_schema"),
            }
        ]
    if operator_row.get("operation") != operation:
        return [], [{"code": "operation_mismatch", "operation": operation}]

    edges: list[tuple[str, Any]] = []
    unresolved: list[dict[str, Any]] = []
    for field in FROZEN_OPERAND_SCHEMA[operation]:
        if field not in operator_row or operator_row.get(field) is None:
            continue
        value = operator_row.get(field)
        values: Iterable[Any]
        labels: Iterable[str]
        if isinstance(value, list):
            values = value
            labels = (field + "[" + str(index) + "]" for index in range(len(value)))
        else:
            values = (value,)
            labels = (field,)
        for label, operand in zip(labels, values):
            if not isinstance(operand, Mapping):
                unresolved.append(
                    {
                        "code": "operand_ref_malformed",
                        "operand": label,
                        "value_kind": type(operand).__name__,
                    }
                )
                continue
            if not isinstance(operand.get("namespace"), str) or not operand.get("namespace"):
                unresolved.append(
                    {"code": "operand_namespace_missing", "operand": label}
                )
                continue
            if "node" not in operand:
                unresolved.append(
                    {"code": "operand_node_missing", "operand": label}
                )
                continue
            edges.append((label, operand))
    return edges, unresolved


def _node_step(
    node_id: str, operation: Any, operand_label: str | None = None
) -> dict[str, Any]:
    return {
        "node_ref": {"namespace": FRESH_NAMESPACE, "node": node_id},
        "operand_label": operand_label,
        "operator": operation if isinstance(operation, str) else None,
        "dereferenced": True,
    }


def _terminal_step(
    namespace: str, node: Any, operand_label: str, status: str
) -> dict[str, Any]:
    return {
        "node_ref": {"namespace": namespace, "node": node},
        "operand_label": operand_label,
        "operator": None,
        "terminal": status,
        "dereferenced": False,
    }


def _search_root(
    root: Mapping[str, Any], nodes: Mapping[str, Mapping[str, Any]]
) -> dict[str, Any]:
    root_ref = root.get("root_ref")
    record: dict[str, Any] = {
        "role": root.get("role"),
        "field": root.get("field"),
        "side": root.get("side"),
        "row_order": root.get("row_order"),
        "original_ref": root.get("original_ref"),
        "root_ref": root_ref,
        "hit": False,
        "path": [],
        "visited_node_count": 0,
        "visited_node_ids": [],
        "saved_boundary_terminations": 0,
        "saved_boundary_target_counts": {},
        "saved_edges": [],
        "external_namespace_terminations": 0,
        "protected_terminations": 0,
        "unresolved": [],
        "no_path_is_not_zero_claim": True,
    }

    if not isinstance(root_ref, Mapping):
        record["unresolved"].append({"code": "root_ref_malformed"})
        return record
    if root_ref.get("namespace") != FRESH_NAMESPACE:
        namespace = root_ref.get("namespace")
        node = root_ref.get("node")
        valid_saved_node = _opaque_node_id(node)
        protected = valid_saved_node and node in PROTECTED_NODE_IDS
        if namespace == SAVED_NAMESPACE:
            record["saved_boundary_terminations"] = 1
            record["saved_boundary_target_counts"] = {str(node): 1}
            target_is_radius = valid_saved_node and node == TARGET_SAVED_NODE
            record["saved_edges"] = [
                {
                    "from_node_ref": None,
                    "operand_label": "<root>",
                    "target_ref": {"namespace": namespace, "node": node},
                    "status": "saved_radius_boundary"
                    if target_is_radius
                    else "saved_boundary_malformed"
                    if not valid_saved_node
                    else "protected_saved_boundary"
                    if protected
                    else "saved_boundary_non_target",
                    "target_is_radius": target_is_radius,
                    "protected": protected,
                    "dereferenced": False,
                }
            ]
            if target_is_radius:
                record["hit"] = True
                record["path"] = [
                    _terminal_step(namespace, node, "<root>", "saved_radius_boundary")
                ]
            elif not valid_saved_node:
                record["unresolved"].append(
                    {"code": "saved_node_malformed", "node": node, "dereferenced": False}
                )
            else:
                record["unresolved"].append(
                    {"code": "root_saved_boundary_non_target", "node": node}
                )
        else:
            record["external_namespace_terminations"] = 1
            record["unresolved"].append(
                {"code": "root_external_namespace", "namespace": namespace, "node": node}
            )
        if protected:
            record["protected_terminations"] = 1
        return record

    root_id = root_ref.get("node")
    if not isinstance(root_id, str) or not root_id:
        record["unresolved"].append({"code": "root_fresh_node_malformed"})
        return record

    # Each root has an independent traversal.  ``prefix`` contains the path to
    # the current node, while ``via_operand`` labels the edge entering it.
    # ``active`` is only for cycle detection; saved/external refs never enter
    # the stack.
    stack: list[tuple[str, list[dict[str, Any]], frozenset[str], str | None]] = [
        (root_id, [], frozenset(), None)
    ]
    visited: set[str] = set()
    while stack:
        node_id, prefix, active, via_operand = stack.pop()
        if node_id in active:
            record["unresolved"].append(
                {"code": "fresh_cycle", "node": node_id}
            )
            continue
        if node_id in visited:
            continue
        node = nodes.get(node_id)
        if node is None:
            record["unresolved"].append(
                {"code": "fresh_child_missing", "node": node_id}
            )
            continue
        if node.get("namespace") != FRESH_NAMESPACE:
            record["unresolved"].append(
                {
                    "code": "fresh_node_namespace_mismatch",
                    "node": node_id,
                    "namespace": node.get("namespace"),
                }
            )
            continue
        visited.add(node_id)
        operation = (
            node.get("operator_row", {}).get("operation")
            if isinstance(node.get("operator_row"), Mapping)
            else None
        )
        current_path = prefix + [_node_step(node_id, operation, via_operand)]
        edges, edge_errors = _operand_edges(node)
        for error in edge_errors:
            item = dict(error)
            item["node"] = node_id
            record["unresolved"].append(item)
        if edge_errors:
            continue
        active_next = active | {node_id}
        # Reverse push gives a stable source-order depth-first witness.
        for operand_label, operand in reversed(edges):
            namespace = operand.get("namespace")
            target = operand.get("node")
            if namespace == FRESH_NAMESPACE:
                if not isinstance(target, str) or not target:
                    record["unresolved"].append(
                        {
                            "code": "fresh_child_malformed",
                            "node": node_id,
                            "operand": operand_label,
                        }
                    )
                    continue
                if target in active_next:
                    record["unresolved"].append(
                        {
                            "code": "fresh_cycle",
                            "node": node_id,
                            "operand": operand_label,
                            "target": target,
                        }
                    )
                    continue
                if target not in nodes:
                    record["unresolved"].append(
                        {
                            "code": "fresh_child_missing",
                            "node": node_id,
                            "operand": operand_label,
                            "target": target,
                        }
                    )
                    continue
                stack.append((target, current_path, active_next, operand_label))
            elif namespace == SAVED_NAMESPACE:
                record["saved_boundary_terminations"] += 1
                key = str(target)
                counts = record["saved_boundary_target_counts"]
                counts[key] = counts.get(key, 0) + 1
                valid_saved_node = _opaque_node_id(target)
                protected = valid_saved_node and target in PROTECTED_NODE_IDS
                target_is_radius = valid_saved_node and target == TARGET_SAVED_NODE
                edge_status = (
                    "saved_radius_boundary"
                    if target_is_radius
                    else "saved_boundary_malformed"
                    if not valid_saved_node
                    else "protected_saved_boundary"
                    if protected
                    else "saved_boundary_non_target"
                )
                saved_edge = {
                    "from_node_ref": {
                        "namespace": FRESH_NAMESPACE,
                        "node": node_id,
                    },
                    "operand_label": operand_label,
                    "target_ref": {
                        "namespace": SAVED_NAMESPACE,
                        "node": target,
                    },
                    "status": edge_status,
                    "target_is_radius": target_is_radius,
                    "protected": protected,
                    "dereferenced": False,
                }
                if target_is_radius:
                    saved_edge["path"] = current_path + [
                        _terminal_step(
                            SAVED_NAMESPACE,
                            target,
                            operand_label,
                            "saved_radius_boundary",
                        )
                    ]
                else:
                    saved_edge["path_to_edge"] = current_path + [
                        _terminal_step(
                            SAVED_NAMESPACE,
                            target,
                            operand_label,
                            edge_status,
                        )
                    ]
                record["saved_edges"].append(saved_edge)
                if protected:
                    record["protected_terminations"] += 1
                if target_is_radius:
                    if not record["hit"]:
                        record["hit"] = True
                        record["path"] = current_path + [
                            _terminal_step(
                                SAVED_NAMESPACE,
                                target,
                                operand_label,
                                "saved_radius_boundary",
                            )
                        ]
                elif not valid_saved_node:
                    record["unresolved"].append(
                        {
                            "code": "saved_node_malformed",
                            "node": node_id,
                            "operand": operand_label,
                            "target": target,
                            "dereferenced": False,
                        }
                    )
                else:
                    record["unresolved"].append(
                        {
                            "code": "saved_boundary_non_target",
                            "node": node_id,
                            "operand": operand_label,
                            "target": target,
                            "dereferenced": False,
                        }
                    )
            else:
                record["external_namespace_terminations"] += 1
                if _opaque_node_id(target) and target in PROTECTED_NODE_IDS:
                    record["protected_terminations"] += 1
                record["unresolved"].append(
                    {
                        "code": "external_namespace_terminal",
                        "node": node_id,
                        "operand": operand_label,
                        "namespace": namespace,
                        "target": target,
                        "dereferenced": False,
                    }
                )

    record["visited_node_count"] = len(visited)
    record["visited_node_ids"] = sorted(visited)
    return record


def _query_body(body: Mapping[str, Any]) -> dict[str, Any]:
    roots = _root_entries(body)
    nodes, malformed = _index_nodes(body)
    root_results = [_search_root(root, nodes) for root in roots]
    unresolved = list(malformed)
    unresolved.extend(
        {
            "root": item["role"] + ".Z" + str(item["row_order"]),
            "code": error.get("code"),
            "detail": error,
        }
        for item in root_results
        for error in item["unresolved"]
    )
    hit_count = sum(bool(item["hit"]) for item in root_results)
    return {
        "schema": SCHEMA,
        "component": COMPONENT,
        "status": "bounded_structural_reachability_complete",
        "source_only": True,
        "fresh_saved_DAG_traversed": True,
        "old_source_graph_traversed": False,
        "prior_schema": PRIOR_SCHEMA,
        "prior_fresh_node_count": EXPECTED_FRESH_NODE_COUNT,
        "fresh_namespace": FRESH_NAMESPACE,
        "saved_namespace": SAVED_NAMESPACE,
        "target_saved_node": TARGET_SAVED_NODE,
        "protected_node_ids": sorted(PROTECTED_NODE_IDS),
        "root_count": len(root_results),
        "hit_count": hit_count,
        "roots": root_results,
        "field_radius_boundary": {
            "original_ref": TARGET_SAVED_NODE,
            "node_ref": {
                "namespace": SAVED_NAMESPACE,
                "node": TARGET_SAVED_NODE,
            },
            "included_as_root_path": False,
            "metadata_only": True,
            "dereferenced": False,
        },
        "saved_boundary_terminations": sum(
            item["saved_boundary_terminations"] for item in root_results
        ),
        "external_namespace_terminations": sum(
            item["external_namespace_terminations"] for item in root_results
        ),
        "protected_terminations": sum(
            item["protected_terminations"] for item in root_results
        ),
        "unresolved_count": len(unresolved),
        "unresolved": unresolved,
        "scope": {
            "only_saved_fresh_nodes_read": True,
            "fresh_saved_DAG_traversed": True,
            "only_fresh_namespace_operand_edges_followed": True,
            "saved_or_external_namespace_terminates": True,
            "saved_boundaries_never_dereferenced": True,
            "protected_refs_never_dereferenced": True,
            "old_source_graph_indexed": False,
            "old_source_graph_traversed": False,
            "old_graph_rebuilt": False,
            "ancestor_replayed": False,
            "source_values_computed": False,
            "numeric_values_computed": False,
            "integrals_evaluated": False,
            "field_radius_boundary_counted_as_root": False,
            "metadata_values_used_as_edges": False,
            "paths_are_structural_only": True,
        },
        "source_gates": {
            "F_field_radius_owner_bound": False,
            "same_source_function_proved": False,
            "function_identity_proved": False,
            "capture_free_substitution_proved": False,
            "ordinary_Z_transport_proved": False,
            "m15_closure_proved": False,
            "local_to_global_m15_identity_proved": False,
            "P_equals_F_proved": False,
            "G_equals_F_proved": False,
            "old_error_terms_removed": False,
            "error_transport_certified": False,
            "numeric_use_admissible": False,
            "full_goal_complete": False,
        },
    }


def build_saved_radius_reachability(
    root: str | Path = ".",
    source: Mapping[str, Any] | None = None,
    hashes: Mapping[str, str] | None = None,
) -> dict[str, Any]:
    """Authenticate a frozen manifest, then query only the saved fresh DAG."""

    _require(isinstance(source, Mapping), "caller frozen manifest is required")
    repo = Path(root).resolve()
    raw, actual, manifest_entry_count = _authenticate(repo, source, hashes)
    prior_raw = raw.get(PRIOR_RESULT.as_posix())
    _require(prior_raw is not None, "authenticated prior result is missing")
    audit_raw = raw.get(PRIOR_AUDIT_RESULT.as_posix())
    _require(audit_raw is not None, "authenticated prior audit result is missing")
    body = _decode_prior(prior_raw)
    audit = _decode_audit(audit_raw)
    result = _query_body(body)
    result["raw_hash_authentication"] = {
        "all_manifest_raw_hashes_checked_before_prior_json_read": True,
        "manifest_entry_count": manifest_entry_count,
        "required_inputs": {
            PRIOR_RESULT.as_posix(): {
                "expected_sha256": PRIOR_RESULT_SHA256,
                "actual_sha256": actual[PRIOR_RESULT.as_posix()],
            },
            PRIOR_AUDIT_RESULT.as_posix(): {
                "expected_sha256": PRIOR_AUDIT_RESULT_SHA256,
                "actual_sha256": actual[PRIOR_AUDIT_RESULT.as_posix()],
            },
        },
    }
    result["prior_result"] = {
        "path": PRIOR_RESULT.as_posix(),
        "sha256": PRIOR_RESULT_SHA256,
        "json_decoded_after_raw_authentication": True,
    }
    result["prior_audit"] = {
        "path": PRIOR_AUDIT_RESULT.as_posix(),
        "sha256": PRIOR_AUDIT_RESULT_SHA256,
        "status": audit.get("status"),
        "schema": audit.get("schema"),
        "passed": True,
        "json_decoded_after_raw_authentication": True,
    }
    result["prior_audit_passed"] = True
    return result


def build_f_gap_saved_radius_reachability(
    root: str | Path = ".",
    source: Mapping[str, Any] | None = None,
    hashes: Mapping[str, str] | None = None,
) -> dict[str, Any]:
    """Explicit F-gap-named alias for :func:`build_saved_radius_reachability`."""

    return build_saved_radius_reachability(root, source=source, hashes=hashes)


def build_radius_reachability(
    root: str | Path = ".",
    source: Mapping[str, Any] | None = None,
    hashes: Mapping[str, str] | None = None,
) -> dict[str, Any]:
    """Runner-compatible alias for the bounded saved-radius query."""

    return build_saved_radius_reachability(root, source=source, hashes=hashes)


def _synthetic_body(nodes: list[dict[str, Any]], root_ref: dict[str, Any]) -> dict[str, Any]:
    roots = {
        "gap_m15": [
            {"row_order": order, "original_ref": 1000 + order, "node_ref": root_ref}
            for order in range(4)
        ],
        "endpoints": {
            side: {
                field: [
                    {
                        "row_order": order,
                        "original_ref": 2000 + order,
                        "node_ref": root_ref,
                    }
                    for order in range(4)
                ]
                for field in ("R_q_y", "q_plus_q_y")
            }
            for side in ("lower", "upper")
        },
        "field_radius": {
            "original_ref": TARGET_SAVED_NODE,
            "node_ref": {"namespace": SAVED_NAMESPACE, "node": TARGET_SAVED_NODE},
        },
    }
    return {
        "schema": PRIOR_SCHEMA,
        "fresh_node_count": EXPECTED_FRESH_NODE_COUNT,
        "nodes": nodes + [
            {
                "namespace": FRESH_NAMESPACE,
                "node_id": "padding_%04d" % index,
                "original_ref": 900000 + index,
                "root_context": {"role": "padding", "row_order": 0},
                "operand_schema": [],
                "operator_row": {"operation": "exact_rational"},
            }
            for index in range(EXPECTED_FRESH_NODE_COUNT - len(nodes))
        ],
        "roots": roots,
        "radius_boundary": {
            "frame_ref": TARGET_SAVED_NODE,
            "original_ref": TARGET_SAVED_NODE,
            "same_opaque_id": True,
        },
    }


def _synthetic_packet(body: Mapping[str, Any]) -> Mapping[str, Any]:
    return {"status": PRIOR_PACKET_STATUS, "result": body}


def verify_synthetic_fixture() -> dict[str, Any]:
    """Return in-memory path and negative controls without repository input."""

    base = {
        "namespace": FRESH_NAMESPACE,
        "node_id": "n00000000",
        "original_ref": 500,
        "root_context": {"role": "fixture", "row_order": 0},
        "operand_schema": ["argument"],
        "operator_row": {
            "operation": "analytic_unary",
            "argument": {
                "namespace": SAVED_NAMESPACE,
                "node": TARGET_SAVED_NODE,
            },
        },
    }
    matching = _query_body(
        _synthetic_body(
            [base],
            {"namespace": FRESH_NAMESPACE, "node": "n00000000"},
        )
    )

    no_path_node = dict(base)
    no_path_node["operator_row"] = {
        "operation": "analytic_unary",
        "argument": {"namespace": SAVED_NAMESPACE, "node": 12345},
    }
    no_path = _query_body(
        _synthetic_body(
            [no_path_node],
            {"namespace": FRESH_NAMESPACE, "node": "n00000000"},
        )
    )

    same_name_external_node = dict(base)
    same_name_external_node["operator_row"] = {
        "operation": "analytic_unary",
        "argument": {"namespace": "external_namespace", "node": TARGET_SAVED_NODE},
    }
    same_name_external = _query_body(
        _synthetic_body(
            [same_name_external_node],
            {"namespace": FRESH_NAMESPACE, "node": "n00000000"},
        )
    )

    protected_node = dict(base)
    protected_node["operator_row"] = {
        "operation": "analytic_unary",
        "argument": {"namespace": SAVED_NAMESPACE, "node": 17},
    }
    protected = _query_body(
        _synthetic_body(
            [protected_node],
            {"namespace": FRESH_NAMESPACE, "node": "n00000000"},
        )
    )

    float_target_node = dict(base)
    float_target_node["operator_row"] = {
        "operation": "analytic_unary",
        "argument": {
            "namespace": SAVED_NAMESPACE,
            "node": float(TARGET_SAVED_NODE),
        },
    }
    float_target = _query_body(
        _synthetic_body(
            [float_target_node],
            {"namespace": FRESH_NAMESPACE, "node": "n00000000"},
        )
    )

    cycle_a = dict(base)
    cycle_a["node_id"] = "cycle_a"
    cycle_a["operator_row"] = {
        "operation": "analytic_unary",
        "argument": {"namespace": FRESH_NAMESPACE, "node": "cycle_b"},
    }
    cycle_b = dict(base)
    cycle_b["node_id"] = "cycle_b"
    cycle_b["operator_row"] = {
        "operation": "analytic_unary",
        "argument": {"namespace": FRESH_NAMESPACE, "node": "cycle_a"},
    }
    cycle = _query_body(
        _synthetic_body(
            [cycle_a, cycle_b],
            {"namespace": FRESH_NAMESPACE, "node": "cycle_a"},
        )
    )

    missing_node = dict(base)
    missing_node["operator_row"] = {
        "operation": "analytic_unary",
        "argument": {"namespace": FRESH_NAMESPACE, "node": "missing"},
    }
    missing = _query_body(
        _synthetic_body(
            [missing_node],
            {"namespace": FRESH_NAMESPACE, "node": "n00000000"},
        )
    )

    checks = {
        "matching_path_hits": matching["hit_count"] == 20
        and all(item["hit"] for item in matching["roots"]),
        "no_path_is_not_zero": no_path["hit_count"] == 0
        and all(item["no_path_is_not_zero_claim"] for item in no_path["roots"]),
        "same_name_different_namespace_terminates": same_name_external["hit_count"] == 0
        and same_name_external["external_namespace_terminations"] == 20,
        "protected_ref_terminates_without_deref": protected["hit_count"] == 0
        and protected["protected_terminations"] == 20,
        "float_saved_target_rejected": float_target["hit_count"] == 0
        and float_target["saved_boundary_terminations"] == 20
        and any(
            item.get("code") == "saved_node_malformed"
            for item in float_target["unresolved"]
        ),
        "cycle_is_reported": cycle["hit_count"] == 0
        and cycle["unresolved_count"] > 0
        and any(item.get("code") == "fresh_cycle" for item in cycle["unresolved"]),
        "missing_child_is_reported": missing["hit_count"] == 0
        and missing["unresolved_count"] > 0
        and any(item.get("code") == "fresh_child_missing" for item in missing["unresolved"]),
    }
    return {
        "schema": SCHEMA + "_synthetic_fixture_v1",
        "component": COMPONENT,
        "status": "synthetic_negative_controls_passed"
        if all(checks.values())
        else "synthetic_negative_controls_failed",
        "synthetic_only": True,
        "actual_source_read": False,
        "fresh_saved_DAG_traversed": True,
        "old_source_graph_traversed": False,
        "source_graph_indexed": False,
        "checks": checks,
        "matching": matching,
        "no_path": no_path,
        "same_name_external": same_name_external,
        "protected": protected,
        "float_target": float_target,
        "cycle": cycle,
        "missing": missing,
        "source_gates": {
            "F_field_radius_owner_bound": False,
            "same_source_function_proved": False,
            "function_identity_proved": False,
            "capture_free_substitution_proved": False,
            "ordinary_Z_transport_proved": False,
            "m15_closure_proved": False,
            "old_error_terms_removed": False,
            "numeric_use_admissible": False,
            "full_goal_complete": False,
        },
    }


__all__ = [
    "COMPONENT",
    "FRESH_NAMESPACE",
    "FROZEN_OPERAND_SCHEMA",
    "FGapSavedRadiusReachabilityInputError",
    "PRIOR_AUDIT_RESULT",
    "PRIOR_AUDIT_RESULT_SHA256",
    "PRIOR_RESULT",
    "PRIOR_RESULT_SHA256",
    "PROTECTED_NODE_IDS",
    "SAVED_NAMESPACE",
    "SCHEMA",
    "TARGET_SAVED_NODE",
    "build_f_gap_saved_radius_reachability",
    "build_radius_reachability",
    "build_saved_radius_reachability",
    "verify_synthetic_fixture",
]
