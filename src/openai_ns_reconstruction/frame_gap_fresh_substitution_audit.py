"""Independent structural audit for the frozen frame-gap substitution receipt.

This module intentionally contains no project imports and has no import-time or
module-entry execution.  ``recover`` reads only the saved receipt and the
authenticated Rp archive when a caller explicitly invokes it.  Passing this
audit does not authenticate semantic or physical equalities.
"""

from __future__ import annotations

import base64
import gzip
import hashlib
import json
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Set, Tuple


RESULT_REPORT = Path("experiments/root_st073/lei_ren_part1_paper_compliant_current_frame_gap_fresh_substitution.json")
RESULT_REPORT_SHA256 = "780a667e601422a1dfb3b13864f8ffa2e8ea12f8c2e875183f025e145d95ff95"
ARCHIVE_RELATIVE = Path(
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_Rp_native_C5_pulse_source.json.gz"
)

RESULT_SHA256 = "a37fc375ff06c06d39b4250ffa3d36c6c315182a81e5d783f4bfb1f981ba0234"
ARCHIVE_RAW_SHA256 = "442df70e5ab4051066b241e943c8d0531067d057e774a04dcc341dadbd703ebe"

FRESH_NAMESPACE = "fresh_frame_gap_substitution_v1"
SAVED_NAMESPACE = "saved_Rp_native_C5_pulse_source"
AUDIT_SCHEMA = "frame_gap_fresh_substitution_audit_v1"
MAX_FAILURE_DETAILS = 100


# This is the frozen source-row operand schema.  It is deliberately repeated
# here rather than imported from the constructor or any earlier incident code.
FROZEN_OPERAND_SCHEMA: Dict[str, Tuple[str, ...]] = {
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


# The five root roles and their four Z-row source references are frozen receipt
# literals.  Keeping the role order explicit prevents an accidental endpoint
# permutation from being accepted as a successful audit.
EXPECTED_ROOT_ORIGINAL_REFS: Dict[str, Tuple[int, ...]] = {
    "gap_m15": (35007, 35008, 61112, 61466),
    "lower.R_q_y": (34762, 34763, 61108, 61462),
    "lower.q_plus_q_y": (34760, 34761, 60604, 60912),
    "upper.R_q_y": (34791, 34792, 61105, 61456),
    "upper.q_plus_q_y": (34789, 34790, 60602, 60910),
}


def _is_int(value: Any) -> bool:
    return isinstance(value, int) and not isinstance(value, bool)


def _read_json(path: Path) -> Tuple[Optional[Any], Optional[str], Optional[str]]:
    """Read JSON and return (value, raw_sha256, error) without project imports."""

    try:
        raw = path.read_bytes()
    except Exception as exc:  # pragma: no cover - error reporting path
        return None, None, "read %s failed: %s" % (path, exc)
    digest = hashlib.sha256(raw).hexdigest()
    try:
        return json.loads(raw.decode("utf-8")), digest, None
    except Exception as exc:  # pragma: no cover - error reporting path
        return None, digest, "parse %s failed: %s" % (path, exc)


class _Audit:
    def __init__(self) -> None:
        self.failures: List[Dict[str, Any]] = []
        self.failure_count = 0
        self.counts: Dict[str, Any] = {
            "source_node_count": 0,
            "saved_mapping_count": 0,
            "fresh_node_count": 0,
            "expected_root_count": 20,
            "roots_checked": 0,
            "fresh_nodes_checked": 0,
            "root_context_checks": 0,
            "root_membership_checks": 0,
            "metadata_checks": 0,
            "operand_fields_checked": 0,
            "operand_refs_checked": 0,
            "saved_operand_checks": 0,
            "fresh_operand_checks": 0,
            "postorder_checks": 0,
        }

    def fail(self, code: str, detail: str, **context: Any) -> None:
        self.failure_count += 1
        if len(self.failures) >= MAX_FAILURE_DETAILS:
            return
        item: Dict[str, Any] = {"code": code, "detail": detail}
        if context:
            item["context"] = context
        self.failures.append(item)

    def report(self, *, result_digest: Optional[str] = None, archive_digest: Optional[str] = None) -> Dict[str, Any]:
        flags = {
            "G_equals_F_proved": False,
            "P_equals_F_proved": False,
            "semantic_equivalence_proved": False,
            "physical_equivalence_proved": False,
            "physical_equalities_authenticated": False,
            "numeric_values_computed": False,
            "source_values_computed": False,
            "integrals_evaluated": False,
            "old_builders_called": False,
            "old_constructors_called": False,
            "old_operators_replayed": False,
            "extra_ledger_terms_added": False,
            "old_error_terms_removed": False,
            "full_goal_complete": False,
            "admissible_for_numeric_use": False,
        }
        return {
            "schema": AUDIT_SCHEMA,
            "status": "passed" if self.failure_count == 0 else "failed",
            "result_sha256": result_digest,
            "expected_result_sha256": RESULT_SHA256,
            "archive_raw_sha256": archive_digest,
            "expected_archive_raw_sha256": ARCHIVE_RAW_SHA256,
            "counts": dict(self.counts),
            "failure_count": self.failure_count,
            "failure_details": list(self.failures),
            "failure_details_truncated": self.failure_count > len(self.failures),
            "semantic_flags": dict(flags),
            # Keep the receipt's safety claims explicit at top level as well.
            **flags,
        }


def _node_at(nodes: Any, index: int) -> Optional[Mapping[str, Any]]:
    if not _is_int(index) or index < 0:
        return None
    item: Any = None
    if isinstance(nodes, list):
        if index >= len(nodes):
            return None
        item = nodes[index]
    elif isinstance(nodes, Mapping):
        item = nodes.get(str(index))
        if item is None:
            item = nodes.get(index)
    if isinstance(item, Mapping):
        return item
    return None


def _root_entries(roots: Any, role: str) -> Optional[Any]:
    if not isinstance(roots, Mapping):
        return None
    if role == "gap_m15":
        return roots.get("gap_m15")
    side, name = role.split(".", 1)
    endpoints = roots.get("endpoints")
    if not isinstance(endpoints, Mapping):
        return None
    side_rows = endpoints.get(side)
    if not isinstance(side_rows, Mapping):
        return None
    return side_rows.get(name)


def _iter_operand_values(value: Any) -> Iterable[Any]:
    if isinstance(value, list):
        return value
    return (value,)


def _load_archive(root: Path, audit: _Audit) -> Tuple[Optional[Any], Optional[str]]:
    path = root / ARCHIVE_RELATIVE
    try:
        raw = path.read_bytes()
    except Exception as exc:  # pragma: no cover - error reporting path
        audit.fail("archive_read", "read archive failed: %s" % exc)
        return None, None
    digest = hashlib.sha256(raw).hexdigest()
    if digest != ARCHIVE_RAW_SHA256:
        audit.fail(
            "archive_raw_sha256_mismatch",
            "archive raw SHA-256 does not match the frozen receipt",
            actual=digest,
            expected=ARCHIVE_RAW_SHA256,
        )
        return None, digest
    try:
        decoded = gzip.decompress(raw)
        archive = json.loads(decoded.decode("utf-8"))
    except Exception as exc:  # pragma: no cover - error reporting path
        audit.fail("archive_decode", "verified archive could not be decoded: %s" % exc)
        return None, digest
    return archive, digest


def _validate_saved_ref(
    audit: _Audit,
    ref: Any,
    original_operand: Any,
    mapping: Mapping[str, Any],
    source_count: int,
    node_id: str,
    field: str,
) -> None:
    audit.counts["saved_operand_checks"] += 1
    if not _is_int(original_operand):
        audit.fail("source_operand_not_index", "source operand is not an integer archive index", node=node_id, field=field)
        return
    if not isinstance(ref, Mapping) or ref.get("namespace") != SAVED_NAMESPACE:
        audit.fail("saved_ref_shape", "saved operand reference has the wrong shape", node=node_id, field=field)
        return
    target = ref.get("node")
    if not _is_int(target):
        audit.fail("saved_target_not_index", "saved target is not an integer archive index", node=node_id, field=field)
        return
    if target < 0 or target >= source_count:
        audit.fail("saved_target_out_of_range", "saved target is outside the source graph", node=node_id, field=field)
    key = str(original_operand)
    if key not in mapping:
        audit.fail("saved_mapping_missing", "saved mapping has no entry for the source operand", node=node_id, field=field, original_ref=original_operand)
    elif mapping.get(key) != target:
        audit.fail(
            "saved_mapping_mismatch",
            "saved operand target differs from rewritten_nodes",
            node=node_id,
            field=field,
            original_ref=original_operand,
            target=target,
            mapped_target=mapping.get(key),
        )


def _validate_fresh_ref(
    audit: _Audit,
    ref: Any,
    original_operand: Any,
    parent: Mapping[str, Any],
    parent_id: str,
    field: str,
    fresh_by_id: Mapping[str, Mapping[str, Any]],
    position_by_id: Mapping[str, int],
) -> None:
    audit.counts["fresh_operand_checks"] += 1
    if not _is_int(original_operand):
        audit.fail("source_operand_not_index", "source operand is not an integer archive index", node=parent_id, field=field)
        return
    if not isinstance(ref, Mapping) or ref.get("namespace") != FRESH_NAMESPACE:
        audit.fail("fresh_ref_shape", "fresh operand reference has the wrong shape", node=parent_id, field=field)
        return
    child_id = ref.get("node")
    child = fresh_by_id.get(child_id) if isinstance(child_id, str) else None
    if child is None:
        audit.fail("fresh_target_missing", "fresh operand references an unknown new node", node=parent_id, field=field, target=child_id)
        return
    if child.get("original_ref") != original_operand:
        audit.fail(
            "fresh_original_ref_mismatch",
            "fresh child original_ref differs from the source operand",
            node=parent_id,
            field=field,
            target=child_id,
            expected=original_operand,
            actual=child.get("original_ref"),
        )
    if child.get("root_context") != parent.get("root_context"):
        audit.fail(
            "fresh_root_context_mismatch",
            "fresh child does not stay in the parent's root_context",
            node=parent_id,
            field=field,
            target=child_id,
        )
    parent_position = position_by_id.get(parent_id)
    child_position = position_by_id.get(child_id)
    audit.counts["postorder_checks"] += 1
    if parent_position is None or child_position is None or child_position >= parent_position:
        audit.fail(
            "fresh_postorder_violation",
            "fresh child is not strictly before its parent",
            node=parent_id,
            field=field,
            target=child_id,
        )


def _validate_node(
    audit: _Audit,
    node: Mapping[str, Any],
    position: int,
    source_nodes: Any,
    source_count: int,
    mapping: Mapping[str, Any],
    fresh_by_id: Mapping[str, Mapping[str, Any]],
    position_by_id: Mapping[str, int],
) -> None:
    node_id = node.get("node_id")
    node_label = node_id if isinstance(node_id, str) else "<malformed>"
    audit.counts["fresh_nodes_checked"] += 1
    original_ref = node.get("original_ref")
    source_row = _node_at(source_nodes, original_ref)
    if source_row is None:
        audit.fail("source_row_missing", "fresh node original_ref has no source row", node=node_label, original_ref=original_ref)
        return
    operation = source_row.get("operation")
    schema = FROZEN_OPERAND_SCHEMA.get(operation)
    if schema is None:
        audit.fail("operation_not_frozen", "source operation is absent from the frozen operand schema", node=node_label, operation=operation)
        return

    declared_schema = node.get("operand_schema")
    if declared_schema != list(schema):
        audit.fail(
            "operand_schema_mismatch",
            "fresh node operand_schema differs from the frozen schema",
            node=node_label,
            expected=list(schema),
            actual=declared_schema,
        )
    operator_row = node.get("operator_row")
    if not isinstance(operator_row, Mapping):
        audit.fail("operator_row_shape", "fresh node operator_row is not an object", node=node_label)
        return
    if operator_row.get("operation") != operation:
        audit.fail("operation_mismatch", "fresh operation differs from the source operation", node=node_label, expected=operation, actual=operator_row.get("operation"))

    source_keys = set(source_row.keys())
    fresh_keys = set(operator_row.keys())
    if source_keys != fresh_keys:
        audit.fail("operator_key_shape", "fresh operator row does not preserve source key shape", node=node_label)
    for key in source_keys - set(schema) - {"operation"}:
        audit.counts["metadata_checks"] += 1
        if operator_row.get(key) != source_row.get(key):
            audit.fail("nonoperand_metadata_mismatch", "non-operand metadata differs from the source row", node=node_label, field=key)

    parent_context = node.get("root_context")
    audit.counts["root_context_checks"] += 1
    if not isinstance(parent_context, Mapping):
        audit.fail("root_context_shape", "fresh node root_context is not an object", node=node_label)

    for field in schema:
        audit.counts["operand_fields_checked"] += 1
        source_has = field in source_row
        fresh_has = field in operator_row
        if source_has != fresh_has:
            audit.fail("operand_field_presence", "source and fresh operand field presence differs", node=node_label, field=field)
            continue
        if not source_has:
            continue
        source_value = source_row.get(field)
        fresh_value = operator_row.get(field)
        source_is_list = isinstance(source_value, list)
        fresh_is_list = isinstance(fresh_value, list)
        if source_is_list != fresh_is_list:
            audit.fail("operand_list_shape", "operand scalar/list shape changed", node=node_label, field=field)
            continue
        if source_is_list and len(source_value) != len(fresh_value):
            audit.fail("operand_list_length", "operand list length changed", node=node_label, field=field)
            continue
        fresh_values = list(_iter_operand_values(fresh_value))
        source_values = list(_iter_operand_values(source_value))
        if len(source_values) != len(fresh_values):
            audit.fail("operand_value_shape", "operand value cardinality changed", node=node_label, field=field)
            continue
        for source_operand, fresh_operand in zip(source_values, fresh_values):
            audit.counts["operand_refs_checked"] += 1
            if source_operand is None:
                if fresh_operand is not None:
                    audit.fail("operand_null_shape", "null operand shape changed", node=node_label, field=field)
                continue
            if not _is_int(source_operand):
                audit.fail("source_operand_not_index", "source operand is not an integer archive index", node=node_label, field=field)
                continue
            if isinstance(fresh_operand, Mapping) and fresh_operand.get("namespace") == SAVED_NAMESPACE:
                _validate_saved_ref(audit, fresh_operand, source_operand, mapping, source_count, node_label, field)
            elif isinstance(fresh_operand, Mapping) and fresh_operand.get("namespace") == FRESH_NAMESPACE:
                _validate_fresh_ref(
                    audit,
                    fresh_operand,
                    source_operand,
                    node,
                    node_label,
                    field,
                    fresh_by_id,
                    position_by_id,
                )
            else:
                audit.fail("operand_namespace", "operand reference namespace is not saved or fresh", node=node_label, field=field)


def _fresh_children(node: Mapping[str, Any], source_nodes: Any) -> Iterable[str]:
    operation = None
    operator_row = node.get("operator_row")
    if isinstance(operator_row, Mapping):
        operation = operator_row.get("operation")
    schema = FROZEN_OPERAND_SCHEMA.get(operation, ())
    if not isinstance(operator_row, Mapping):
        return ()
    children: List[str] = []
    for field in schema:
        for ref in _iter_operand_values(operator_row.get(field)):
            if isinstance(ref, Mapping) and ref.get("namespace") == FRESH_NAMESPACE and isinstance(ref.get("node"), str):
                children.append(ref["node"])
    return children


def _validate_root_structure(
    audit: _Audit,
    result: Mapping[str, Any],
    fresh_by_id: Mapping[str, Mapping[str, Any]],
    mapping: Mapping[str, Any],
    source_count: int,
) -> Dict[str, Set[str]]:
    roots = result.get("roots")
    owners: Dict[str, Set[str]] = {node_id: set() for node_id in fresh_by_id}
    if not isinstance(roots, Mapping):
        audit.fail("roots_shape", "result roots is not an object")
        return owners

    expected_root_keys = {"gap_m15", "endpoints", "field_radius"}
    if set(roots.keys()) != expected_root_keys:
        audit.fail("roots_key_shape", "result roots does not have the frozen key shape")

    for role, expected_refs in EXPECTED_ROOT_ORIGINAL_REFS.items():
        entries = _root_entries(roots, role)
        if not isinstance(entries, list) or len(entries) != len(expected_refs):
            audit.fail("root_row_shape", "root role does not contain exactly four Z rows", role=role)
            continue
        for row_order, expected_original_ref in enumerate(expected_refs):
            audit.counts["roots_checked"] += 1
            entry = entries[row_order]
            if not isinstance(entry, Mapping):
                audit.fail("root_entry_shape", "root entry is not an object", role=role, row_order=row_order)
                continue
            if entry.get("row_order") != row_order:
                audit.fail("root_row_order", "root row_order is not the frozen Z order", role=role, row_order=row_order)
            if entry.get("original_ref") != expected_original_ref:
                audit.fail(
                    "root_original_ref",
                    "root original_ref differs from the frozen twenty-root receipt",
                    role=role,
                    row_order=row_order,
                    expected=expected_original_ref,
                    actual=entry.get("original_ref"),
                )
            node_ref = entry.get("node_ref")
            node_id = node_ref.get("node") if isinstance(node_ref, Mapping) else None
            if not isinstance(node_ref, Mapping) or node_ref.get("namespace") != FRESH_NAMESPACE or not isinstance(node_id, str):
                audit.fail("root_node_ref", "root does not reference a fresh node", role=role, row_order=row_order)
                continue
            node = fresh_by_id.get(node_id)
            if node is None:
                audit.fail("root_node_missing", "root references an unknown fresh node", role=role, row_order=row_order, target=node_id)
                continue
            if node.get("original_ref") != expected_original_ref:
                audit.fail("root_node_original_ref", "root node original_ref differs from root entry", role=role, row_order=row_order, target=node_id)
            expected_context = {"role": role, "row_order": row_order}
            if node.get("root_context") != expected_context:
                audit.fail("root_context", "root node root_context is not the frozen role and Z row", role=role, row_order=row_order, target=node_id)
            root_key = "%s#%d" % (role, row_order)
            stack = [node_id]
            visited: Set[str] = set()
            while stack:
                current = stack.pop()
                if current in visited:
                    continue
                visited.add(current)
                if current not in fresh_by_id:
                    audit.fail("root_reachability", "root traversal reached an unknown fresh node", role=role, row_order=row_order, target=current)
                    continue
                owners.setdefault(current, set()).add(root_key)
                audit.counts["root_membership_checks"] += 1
                if fresh_by_id[current].get("root_context") != expected_context:
                    audit.fail("reachable_root_context", "reachable node does not stay in the root_context", role=role, row_order=row_order, target=current)
                stack.extend(_fresh_children(fresh_by_id[current], result.get("nodes")))

    field_radius = roots.get("field_radius")
    if not isinstance(field_radius, Mapping):
        audit.fail("radius_root_shape", "field_radius root is not an object")
    else:
        if field_radius.get("original_ref") != 33982:
            audit.fail("radius_original_ref", "field_radius original_ref is not the frozen radius node")
        radius_ref = field_radius.get("node_ref")
        if not isinstance(radius_ref, Mapping) or radius_ref.get("namespace") != SAVED_NAMESPACE or radius_ref.get("node") != 33982:
            audit.fail("radius_ref", "field_radius does not retain the saved opaque radius reference")
        _validate_saved_ref(audit, radius_ref, 33982, mapping, source_count, "field_radius", "field_radius")
    radius_boundary = result.get("radius_boundary")
    expected_boundary = {"frame_ref": 33982, "original_ref": 33982, "same_opaque_id": True}
    if radius_boundary != expected_boundary:
        audit.fail("radius_boundary", "radius_boundary is not the frozen opaque boundary")
    return owners


def recover(root: Path) -> Dict[str, Any]:
    """Audit the saved rewrite receipt without evaluating or constructing a DAG."""

    audit = _Audit()
    report_path = root / RESULT_REPORT
    try:
        report_raw = report_path.read_bytes()
        if hashlib.sha256(report_raw).hexdigest() != RESULT_REPORT_SHA256:
            raise ValueError("published result report hash mismatch")
        report = json.loads(report_raw)
        raw = base64.b64decode(report["artifacts"]["tmp/frame_gap_fresh_substitution_v1/result.json"]["base64"], validate=True)
        result_digest = hashlib.sha256(raw).hexdigest()
        result = json.loads(raw)
    except Exception as exc:
        audit.fail("result_input", "published result decode failed: " + str(exc))
        return audit.report()
    if result_digest != RESULT_SHA256:
        audit.fail(
            "result_sha256_mismatch",
            "result raw SHA-256 does not match the frozen receipt",
            actual=result_digest,
            expected=RESULT_SHA256,
        )
        return audit.report(result_digest=result_digest)
    if not isinstance(result, Mapping) or not isinstance(result.get("result"), Mapping):
        audit.fail("result_shape", "receipt does not contain an object-valued result")
        return audit.report(result_digest=result_digest)
    result_body = result["result"]

    archive, archive_digest = _load_archive(root, audit)
    if archive is None:
        return audit.report(result_digest=result_digest, archive_digest=archive_digest)
    if result_body.get("frame_archive_raw_sha256") != ARCHIVE_RAW_SHA256:
        audit.fail("receipt_archive_sha256", "receipt archive hash differs from the frozen archive hash")
    if result_body.get("schema") != "frame_gap_fresh_substitution_v1":
        audit.fail("receipt_schema", "receipt schema is not frame_gap_fresh_substitution_v1")
    if result_body.get("namespace") != FRESH_NAMESPACE:
        audit.fail("fresh_namespace", "receipt fresh namespace is not frozen")
    if result_body.get("existing_frame_namespace") != SAVED_NAMESPACE:
        audit.fail("saved_namespace", "receipt saved namespace is not frozen")

    for flag in (
        "source_values_computed",
        "integrals_evaluated",
        "old_builders_called",
        "G_equals_F_proved",
        "P_equals_F_proved",
        "admissible_for_numeric_use",
        "full_goal_complete",
        "extra_ledger_terms_added",
        "old_error_terms_removed",
    ):
        if result_body.get(flag) is not False:
            audit.fail("receipt_safety_flag", "receipt safety flag is not explicitly false", flag=flag)

    fresh_nodes = result_body.get("nodes")
    if not isinstance(fresh_nodes, list):
        audit.fail("fresh_nodes_shape", "receipt nodes is not a list")
        return audit.report(result_digest=result_digest, archive_digest=archive_digest)
    audit.counts["fresh_node_count"] = len(fresh_nodes)
    if result_body.get("fresh_node_count") != len(fresh_nodes):
        audit.fail("fresh_node_count", "declared fresh_node_count differs from the saved DAG")

    if not isinstance(archive, Mapping):
        audit.fail("archive_shape", "decoded archive is not an object")
        return audit.report(result_digest=result_digest, archive_digest=archive_digest)
    source_nodes = archive.get("exact_graph_nodes")
    frame = archive.get("Rp_native_C5_pulse_source_bindings")
    if not isinstance(source_nodes, (list, Mapping)):
        audit.fail("source_nodes_shape", "archive exact_graph_nodes is not an indexable collection")
        return audit.report(result_digest=result_digest, archive_digest=archive_digest)
    source_count = len(source_nodes)
    audit.counts["source_node_count"] = source_count
    if not isinstance(frame, Mapping) or not isinstance(frame.get("rewritten_nodes"), Mapping):
        audit.fail("saved_mapping_shape", "archive saved rewritten_nodes is not an object")
        return audit.report(result_digest=result_digest, archive_digest=archive_digest)
    mapping = frame["rewritten_nodes"]
    audit.counts["saved_mapping_count"] = len(mapping)
    for key, value in mapping.items():
        if not isinstance(key, str) or not key.isdigit() or not _is_int(value) or value < 0 or value >= source_count:
            audit.fail("saved_mapping_entry", "saved rewritten_nodes contains a malformed entry", key=key)

    fresh_by_id: Dict[str, Mapping[str, Any]] = {}
    position_by_id: Dict[str, int] = {}
    for position, raw_node in enumerate(fresh_nodes):
        if not isinstance(raw_node, Mapping):
            audit.fail("fresh_node_shape", "fresh DAG entry is not an object", position=position)
            continue
        node_id = raw_node.get("node_id")
        if not isinstance(node_id, str):
            audit.fail("fresh_node_id", "fresh DAG entry node_id is not a string", position=position)
            continue
        if node_id in fresh_by_id:
            audit.fail("fresh_node_id_duplicate", "fresh DAG node_id is duplicated", node=node_id)
            continue
        fresh_by_id[node_id] = raw_node
        position_by_id[node_id] = position

    for position, raw_node in enumerate(fresh_nodes):
        if isinstance(raw_node, Mapping) and isinstance(raw_node.get("node_id"), str):
            _validate_node(
                audit,
                raw_node,
                position,
                source_nodes,
                source_count,
                mapping,
                fresh_by_id,
                position_by_id,
            )
            if raw_node.get("namespace") != FRESH_NAMESPACE:
                audit.fail("fresh_node_namespace", "fresh node namespace is not frozen", node=raw_node.get("node_id"))
            if not _is_int(raw_node.get("original_ref")) or _node_at(source_nodes, raw_node.get("original_ref")) is None:
                audit.fail("fresh_original_ref", "fresh node original_ref is not a valid source index", node=raw_node.get("node_id"))

    owners = _validate_root_structure(audit, result_body, fresh_by_id, mapping, source_count)
    for node_id in fresh_by_id:
        root_owners = owners.get(node_id, set())
        if not root_owners:
            audit.fail("unowned_fresh_node", "fresh node is not reachable from any frozen root", node=node_id)
        elif len(root_owners) != 1:
            audit.fail("multi_owned_fresh_node", "fresh node is reachable from more than one root", node=node_id)

    audit.counts["unowned_fresh_nodes"] = sum(1 for node_id in fresh_by_id if not owners.get(node_id))
    return audit.report(result_digest=result_digest, archive_digest=archive_digest)
