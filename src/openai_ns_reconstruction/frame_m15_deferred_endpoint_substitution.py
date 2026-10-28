"""Unevaluated independent source DAG for deferred m15 endpoint rows.

The published mismatch result selects the forty-four rows.  This module reads
that result, the independently audited fresh-gap result, and the authenticated
Rp archive, then builds a new namespace with one memoized DAG per occurrence.
Saved ``rewritten_nodes`` and prior fresh nodes are opaque boundary caches.
No project constructor, checker, old operator, value, integral, P/F comparison,
or endpoint identity is invoked here.
"""

from __future__ import annotations

import base64
import copy
import gzip
import hashlib
import json
from collections.abc import Mapping
from pathlib import Path
from typing import Any


MISMATCH_REPORT_RELATIVE = Path(
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_frame_m15_endpoint_mismatch.json"
)
MISMATCH_REPORT_SHA256 = (
    "cb8d8aa0cf5a877b03bcda35f7cae32c6761171d9b739e4035ca2302f938d884"
)
MISMATCH_RESULT_ARTIFACT = "tmp/frame_m15_endpoint_mismatch_v1/result.json"
MISMATCH_RESULT_SHA256 = (
    "6adc3c1a272b2fcc6f6259c8ad844b8560a9893ac4efdcd8dc689e7bcd177d93"
)

ARCHIVE_RELATIVE = Path(
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_Rp_native_C5_pulse_source.json.gz"
)
ARCHIVE_RAW_SHA256 = (
    "442df70e5ab4051066b241e943c8d0531067d057e774a04dcc341dadbd703ebe"
)

FRESH_REPORT_RELATIVE = Path(
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_frame_gap_fresh_substitution.json"
)
FRESH_REPORT_SHA256 = (
    "780a667e601422a1dfb3b13864f8ffa2e8ea12f8c2e875183f025e145d95ff95"
)
FRESH_RESULT_ARTIFACT = "tmp/frame_gap_fresh_substitution_v1/result.json"
FRESH_RESULT_SHA256 = (
    "a37fc375ff06c06d39b4250ffa3d36c6c315182a81e5d783f4bfb1f981ba0234"
)

AUDIT_REPORT_RELATIVE = Path(
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_frame_gap_fresh_substitution_audit.json"
)
AUDIT_REPORT_SHA256 = (
    "de6de458582dab9b77eee5ed81ef6d1669bff245416aa118536d3fadb68f5432"
)
AUDIT_RESULT_ARTIFACT = "tmp/frame_gap_fresh_substitution_audit_v1/result.json"
AUDIT_RESULT_SHA256 = (
    "a3de519abfdc42fe18a0fa6830105c9c2c39a7bb4188d27693123117123087a3"
)

C3_ARCHIVE_NAME = (
    "lei_ren_part1_paper_compliant_current_original_n1_pulse_pressure_C3_functions.json.gz"
)
C3_ARCHIVE_SHA256 = (
    "bdfcc70f57d8cea981f508074690a7e940f884382c2452805f126a3f96bdb59d"
)

SAVED_NAMESPACE = "saved_Rp_native_C5_pulse_source"
PRIOR_FRESH_NAMESPACE = "fresh_frame_gap_substitution_v1"
FRESH_NAMESPACE = "fresh_frame_m15_deferred_endpoint_substitution_v1"
DEFERRED_NAMESPACE = "unexecuted_saved_inlet_substitution_v1"

M15_ENDPOINT_FIELDS = (
    "R_squared_q",
    "R_squared_q_squared",
    "R_squared_U_squared",
    "R_q_y",
)
ORDINARY_Z_ROWS = (0, 1)
ENDPOINT_LAYOUT = (
    ("main.lower", "chart", "main", 0),
    ("main.upper", "chart", "main", 1),
    ("gap.lower", "chart", "gap", 0),
    ("gap.upper", "chart", "gap", 1),
    ("end.lower", "chart", "end", 0),
    ("end.upper", "chart", "end", 1),
    ("global.lower", "global", "global", 0),
    ("global.upper", "global", "global", 1),
)
EXPECTED_GAP_R_Q_Y = {"lower": (34762, 34763), "upper": (34791, 34792)}

# Frozen and local: only these fields are graph operands.  All other source
# row attributes are copied without interpretation.
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
        "active_body", "Delta", "eta", "flat_predicate",
        "original_collar_coordinate", "original_collar_upper",
    ),
    "original_monotone_phase_inverse": (
        "phase_function", "target_phase", "angle_lower", "angle_upper",
        "Delta", "eta",
    ),
}
REFERENCE_KEYS = FROZEN_OPERAND_SCHEMA
LIST_OPERAND_FIELDS = frozenset({"arguments"})


def _read_hashed_bytes(root: Path, relative: Path, expected: str) -> bytes:
    raw = (root / relative).read_bytes()
    if hashlib.sha256(raw).hexdigest() != expected:
        raise ValueError("authenticated input hash mismatch: " + str(relative))
    return raw


def _mapping(value: Any, label: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise ValueError(label + " must be an object")
    return value


def _list(value: Any, label: str) -> list[Any]:
    if not isinstance(value, list):
        raise ValueError(label + " must be a list")
    return value


def _index(value: Any, label: str, count: int) -> int:
    if type(value) is not int or value < 0 or value >= count:
        raise ValueError(label + " is outside authenticated archive graph")
    return value


def _json_result_from_report(
    root: Path,
    relative: Path,
    report_sha256: str,
    artifact_name: str,
    result_sha256: str,
    label: str,
) -> tuple[dict[str, Any], bytes]:
    report = _mapping(
        json.loads(_read_hashed_bytes(root, relative, report_sha256).decode("utf-8")),
        label + " report",
    )
    artifact = _mapping(
        _mapping(report.get("artifacts"), label + " artifacts").get(artifact_name),
        label + " artifact",
    )
    encoded = artifact.get("base64")
    if not isinstance(encoded, str):
        raise ValueError(label + " artifact has no base64 payload")
    try:
        raw = base64.b64decode(encoded, validate=True)
    except Exception as exc:
        raise ValueError(label + " artifact base64 is invalid") from exc
    if hashlib.sha256(raw).hexdigest() != result_sha256:
        raise ValueError(label + " result hash mismatch")
    if artifact.get("sha256") not in (None, result_sha256):
        raise ValueError(label + " artifact metadata hash mismatch")
    packet = _mapping(json.loads(raw.decode("utf-8")), label + " packet")
    if packet.get("status") != "returned":
        raise ValueError(label + " packet is not returned")
    return dict(_mapping(packet.get("result"), label + " result")), raw


def _load_archive(root: Path) -> tuple[dict[str, Any], bytes]:
    raw = _read_hashed_bytes(root, ARCHIVE_RELATIVE, ARCHIVE_RAW_SHA256)
    try:
        archive = json.loads(gzip.decompress(raw).decode("utf-8"))
    except Exception as exc:
        raise ValueError("authenticated Rp archive cannot be decoded") from exc
    return dict(_mapping(archive, "Rp archive")), raw


def _canonical_mapping(value: Any, label: str, source_count: int) -> dict[str, int]:
    mapping = _mapping(value, label)
    if not mapping:
        raise ValueError(label + " must be nonempty")
    out: dict[str, int] = {}
    for key, target in mapping.items():
        if not isinstance(key, str) or not key.isdecimal() or str(int(key)) != key:
            raise ValueError(label + " contains a noncanonical source index")
        if key in out:
            raise ValueError(label + " contains a duplicate source index")
        out[key] = _index(target, label + "[" + key + "]", source_count)
    return out


def _substitution_digest(substitution: Mapping[str, Any]) -> str:
    payload = json.dumps(
        substitution, sort_keys=True, separators=(",", ":"), ensure_ascii=False
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _endpoint_from_c3(c3: Mapping[str, Any], label: str) -> Mapping[str, Any]:
    for candidate, kind, container, index in ENDPOINT_LAYOUT:
        if candidate != label:
            continue
        if kind == "global":
            rows = c3.get("global_endpoints")
            if not isinstance(rows, list) or len(rows) != 2:
                raise ValueError("saved C3 global_endpoints must contain lower and upper")
            return _mapping(rows[index], label + " archive endpoint")
        charts = c3.get("endpoints")
        if isinstance(charts, list):
            if len(charts) != 3:
                raise ValueError("saved C3 endpoints list has wrong chart count")
            chart = charts[("main", "gap", "end").index(container)]
        elif isinstance(charts, Mapping):
            chart = charts.get(container)
        else:
            raise ValueError("saved C3 endpoints are missing")
        if isinstance(chart, list):
            if len(chart) != 2:
                raise ValueError(label + " archive endpoint list has wrong length")
            return _mapping(chart[index], label + " archive endpoint")
        if isinstance(chart, Mapping):
            return _mapping(chart.get("lower" if index == 0 else "upper"), label + " archive endpoint")
        raise ValueError(label + " archive chart endpoint is missing")
    raise ValueError("unknown endpoint layout label: " + label)


def _select_report_rows(
    root: Path,
    archive: Mapping[str, Any],
    source_count: int,
    rewritten_nodes: Mapping[str, int],
    substitution_digest: str,
) -> tuple[dict[str, Any], list[dict[str, Any]], list[dict[str, Any]]]:
    """Authenticate the prior result and select rows by its real status."""

    result, _ = _json_result_from_report(
        root,
        MISMATCH_REPORT_RELATIVE,
        MISMATCH_REPORT_SHA256,
        MISMATCH_RESULT_ARTIFACT,
        MISMATCH_RESULT_SHA256,
        "frame m15 mismatch",
    )
    if result.get("schema") != "frame_m15_endpoint_mismatch_v1":
        raise ValueError("frame m15 mismatch result schema mismatch")
    if result.get("status") != "unexecuted_source_bound_endpoint_mismatch":
        raise ValueError("frame m15 mismatch result status mismatch")
    if result.get("source_references_are_opaque") is not True:
        raise ValueError("frame m15 mismatch result is not source-opaque")
    if result.get("endpoint_count") != 8 or result.get("source_row_count") != 64 or result.get("ordinary_Z_rows") != [0, 1]:
        raise ValueError("frame m15 top-level row scope disagrees")
    provenance = _mapping(result.get("provenance"), "frame m15 provenance")
    for key, expected in (
        ("expected_archive_raw_sha256", ARCHIVE_RAW_SHA256),
        ("archive_raw_sha256", ARCHIVE_RAW_SHA256),
        ("C3_archive_name", C3_ARCHIVE_NAME),
        ("C3_archive_sha256", C3_ARCHIVE_SHA256),
        ("saved_input_substitution_sha256", substitution_digest),
    ):
        if provenance.get(key) != expected:
            raise ValueError("frame m15 provenance disagrees: " + key)
    counts = _mapping(result.get("binding_status_counts"), "frame m15 binding counts")
    if dict(counts) != {
        "missing_saved_reference": 44,
        "saved_frame_reference": 16,
        "reused_committed_gap_R_q_y_reference": 4,
    }:
        raise ValueError("frame m15 result does not contain the published 44/16/4 scope")

    frame = _mapping(
        archive.get("Rp_native_C5_pulse_source_bindings"), "Rp frame bindings"
    )
    c3 = _mapping(
        _mapping(frame.get("inputs"), "Rp frame inputs").get(C3_ARCHIVE_NAME),
        "saved C3 pulse input",
    )
    endpoint_map = _mapping(result.get("endpoints"), "frame m15 endpoints")
    if set(endpoint_map) != {item[0] for item in ENDPOINT_LAYOUT}:
        raise ValueError("frame m15 endpoint labels disagree")

    missing: list[dict[str, Any]] = []
    retained: list[dict[str, Any]] = []
    for label, _, _, _ in ENDPOINT_LAYOUT:
        endpoint = _mapping(endpoint_map.get(label), label + " result")
        if endpoint.get("endpoint") != label:
            raise ValueError(label + " endpoint label disagrees")
        fields = _mapping(endpoint.get("fields"), label + ".fields")
        if set(fields) != set(M15_ENDPOINT_FIELDS):
            raise ValueError(label + " has unexpected endpoint field roles")
        archive_endpoint = _endpoint_from_c3(c3, label)
        for field in M15_ENDPOINT_FIELDS:
            selected = _mapping(fields.get(field), label + "." + field)
            originals = _list(selected.get("original_refs_Z0_Z1"), label + "." + field + ".refs")
            if len(originals) != 2:
                raise ValueError(label + "." + field + " does not have Z0/Z1")
            rows = _mapping(selected.get("rows"), label + "." + field + ".rows")
            archive_rows = _list(archive_endpoint.get(field), label + " archive " + field)
            if len(archive_rows) < 2:
                raise ValueError(label + " archive " + field + " lacks Z0/Z1")
            for order in ORDINARY_Z_ROWS:
                z_label = "Z" + str(order)
                original = _index(originals[order], label + "." + field + "." + z_label, source_count)
                row = _mapping(rows.get(z_label), label + "." + field + "." + z_label)
                if _index(row.get("original_ref"), "report original_ref", source_count) != original:
                    raise ValueError(label + "." + field + "." + z_label + " ref disagrees")
                if _index(archive_rows[order], label + " archive row", source_count) != original:
                    raise ValueError("source frame input disagrees with published mismatch row")
                occurrence = {
                    "endpoint": label,
                    "field": field,
                    "z_order": order,
                    "z_label": z_label,
                    "original_ref": original,
                    "report_status": row.get("status"),
                    "report_binding": copy.deepcopy(dict(row)),
                }
                status = row.get("status")
                if status == "missing_saved_reference":
                    if str(original) in rewritten_nodes:
                        raise ValueError("missing report row is present in saved rewritten_nodes")
                    deferred = _mapping(row.get("saved_inlet_substitution"), "missing substitution")
                    if deferred.get("executed") is not False:
                        raise ValueError("missing substitution is marked executed")
                    if deferred.get("namespace") != DEFERRED_NAMESPACE:
                        raise ValueError("missing substitution namespace disagrees")
                    if deferred.get("mapping_sha256") != substitution_digest:
                        raise ValueError("missing substitution digest disagrees")
                    if deferred.get("operation") != "simultaneous_saved_inlet_substitution" or deferred.get("mapping_pointer") != "$.Rp_native_C5_pulse_source_bindings.input_substitution":
                        raise ValueError("missing substitution operation or owner disagrees")
                    missing.append(occurrence)
                elif status == "saved_frame_reference":
                    target = _index(row.get("node"), "saved report target", source_count)
                    if rewritten_nodes.get(str(original)) != target:
                        raise ValueError("saved report target disagrees with rewritten_nodes")
                    retained.append(occurrence)
                elif status == "reused_committed_gap_R_q_y_reference":
                    if label not in ("gap.lower", "gap.upper") or field != "R_q_y":
                        raise ValueError("prior fresh status appears on an unexpected field")
                    retained.append(occurrence)
                else:
                    raise ValueError("unknown frame m15 endpoint binding status")
    if len(missing) != 44 or len(retained) != 20:
        raise ValueError("published frame m15 occurrence scope is not 44 plus 20")
    return result, missing, retained


def _load_prior_boundaries(
    root: Path, source_count: int, rewritten_nodes: Mapping[str, int]
) -> tuple[dict[int, list[str]], dict[str, Any], set[int], dict[str, Any]]:
    """Index the previously audited DAG by original reference as opaque boundaries."""

    fresh, _ = _json_result_from_report(
        root,
        FRESH_REPORT_RELATIVE,
        FRESH_REPORT_SHA256,
        FRESH_RESULT_ARTIFACT,
        FRESH_RESULT_SHA256,
        "prior fresh substitution",
    )
    if fresh.get("schema") != "frame_gap_fresh_substitution_v1":
        raise ValueError("prior fresh result schema mismatch")
    if fresh.get("status") != "unevaluated_local_source_dag":
        raise ValueError("prior fresh result status mismatch")
    if fresh.get("fresh_node_count") != 1131 or fresh.get("missing_root_count") != 20:
        raise ValueError("prior fresh count metadata mismatch")
    if fresh.get("frame_archive_raw_sha256") != ARCHIVE_RAW_SHA256:
        raise ValueError("prior fresh result is pinned to another archive")
    nodes = _list(fresh.get("nodes"), "prior fresh nodes")
    if len(nodes) != 1131:
        raise ValueError("prior fresh node list count mismatch")
    node_ids: set[str] = set()
    for position, raw_node in enumerate(nodes):
        node = _mapping(raw_node, "prior fresh node " + str(position))
        if node.get("namespace") != PRIOR_FRESH_NAMESPACE:
            raise ValueError("prior fresh node namespace mismatch")
        node_id = node.get("node_id")
        if not isinstance(node_id, str) or not node_id or node_id in node_ids:
            raise ValueError("prior fresh node id is not unique")
        node_ids.add(node_id)
        _index(node.get("original_ref"), "prior fresh original_ref", source_count)
        if not isinstance(node.get("root_context"), Mapping):
            raise ValueError("prior fresh root context is missing")
    # The independently saved prior audit is authenticated below; do not replay it.
    prior_by_original: dict[int, list[str]] = {}
    for raw_node in nodes:
        node = _mapping(raw_node, "prior fresh node")
        original = _index(node.get("original_ref"), "prior fresh original_ref", source_count)
        prior_by_original.setdefault(original, []).append(str(node["node_id"]))

    audit, _ = _json_result_from_report(
        root,
        AUDIT_REPORT_RELATIVE,
        AUDIT_REPORT_SHA256,
        AUDIT_RESULT_ARTIFACT,
        AUDIT_RESULT_SHA256,
        "prior fresh structural audit",
    )
    if audit.get("schema") != "frame_gap_fresh_substitution_audit_v1" or audit.get("status") != "passed":
        raise ValueError("prior fresh structural audit is not passed")
    if audit.get("archive_raw_sha256") != ARCHIVE_RAW_SHA256 or audit.get("expected_archive_raw_sha256") != ARCHIVE_RAW_SHA256:
        raise ValueError("prior fresh audit archive hash mismatch")
    if audit.get("expected_result_sha256") != FRESH_RESULT_SHA256:
        raise ValueError("prior fresh audit result hash mismatch")
    audit_counts = _mapping(audit.get("counts"), "prior fresh audit counts")
    if audit_counts.get("fresh_node_count") != 1131 or audit_counts.get("fresh_nodes_checked") != 1131:
        raise ValueError("prior fresh audit node count mismatch")
    if audit_counts.get("expected_root_count") != 20 or audit_counts.get("roots_checked") != 20:
        raise ValueError("prior fresh audit root count mismatch")
    if audit.get("failure_count") != 0 or audit.get("failure_details") != []:
        raise ValueError("prior fresh audit contains failures")
    semantic_flags = _mapping(audit.get("semantic_flags"), "prior fresh audit flags")
    if any(value is not False for value in semantic_flags.values()):
        raise ValueError("prior fresh audit semantic flag is not false")
    for key in (
        "source_values_computed", "integrals_evaluated", "numeric_values_computed",
        "old_builders_called", "old_constructors_called", "old_operators_replayed",
        "old_error_terms_removed", "extra_ledger_terms_added", "full_goal_complete",
    ):
        if audit.get(key) is not False:
            raise ValueError("prior fresh audit flag is not false: " + key)

    roots = _mapping(fresh.get("roots"), "prior fresh roots")
    prior_root_originals: set[int] = set()

    def check_roots(value: Any, label: str) -> None:
        rows = _list(value, label)
        if len(rows) != 4:
            raise ValueError(label + " must contain four rows")
        for position, raw_row in enumerate(rows):
            row = _mapping(raw_row, label + " row")
            original = _index(row.get("original_ref"), label + ".original_ref", source_count)
            node_ref = _mapping(row.get("node_ref"), label + ".node_ref")
            if node_ref.get("namespace") != PRIOR_FRESH_NAMESPACE:
                raise ValueError(label + " root is not prior fresh")
            node_id = node_ref.get("node")
            if not isinstance(node_id, str) or node_id not in node_ids:
                raise ValueError(label + " root node is unknown")
            if node_id not in prior_by_original.get(original, []):
                raise ValueError(label + " root original/node pair disagrees")
            if row.get("row_order") != position:
                raise ValueError(label + " row order disagrees")
            prior_root_originals.add(original)

    check_roots(roots.get("gap_m15"), "prior fresh gap_m15")
    endpoint_roots = _mapping(roots.get("endpoints"), "prior fresh endpoint roots")
    for side in ("lower", "upper"):
        side_roots = _mapping(endpoint_roots.get(side), "prior fresh " + side + " roots")
        check_roots(side_roots.get("R_q_y"), "prior fresh " + side + " R_q_y")
        check_roots(side_roots.get("q_plus_q_y"), "prior fresh " + side + " q_plus_q_y")
    radius = _mapping(roots.get("field_radius"), "prior fresh radius")
    radius_original = _index(radius.get("original_ref"), "prior fresh radius original_ref", source_count)
    radius_ref = _mapping(radius.get("node_ref"), "prior fresh radius node_ref")
    if radius_ref.get("namespace") != SAVED_NAMESPACE or rewritten_nodes.get(str(radius_original)) != radius_ref.get("node"):
        raise ValueError("prior fresh radius boundary disagrees with saved map")
    receipt = {
        "report_sha256": AUDIT_REPORT_SHA256,
        "result_sha256": AUDIT_RESULT_SHA256,
        "schema": audit.get("schema"),
        "status": audit.get("status"),
        "fresh_node_count": len(nodes),
        "fresh_root_count": len(prior_root_originals),
        "semantic_flags": dict(semantic_flags),
    }
    return prior_by_original, receipt, prior_root_originals, fresh


def _boundary_ref(
    original: int,
    rewritten_nodes: Mapping[str, int],
    prior_by_original: Mapping[int, list[str]],
    source_count: int,
) -> tuple[dict[str, Any] | None, dict[str, Any] | None]:
    _index(original, "boundary original_ref", source_count)
    saved = rewritten_nodes.get(str(original))
    if saved is not None:
        return (
            {"namespace": SAVED_NAMESPACE, "node": saved},
            {"kind": "saved_rewritten_nodes", "namespace": SAVED_NAMESPACE,
             "original_ref": original, "node": saved},
        )
    candidates = prior_by_original.get(original)
    if candidates:
        chosen = candidates[0]
        return (
            {"namespace": PRIOR_FRESH_NAMESPACE, "node": chosen},
            {"kind": "prior_fresh_boundary", "namespace": PRIOR_FRESH_NAMESPACE,
             "original_ref": original, "node": chosen,
             "candidate_nodes": list(candidates),
             "selection": "deterministic_first_prior_result_order",
             "duplicate_original_ref": len(candidates) > 1},
        )
    return None, None


def _operand_children(
    row: Mapping[str, Any], operation: str, source_count: int, label: str
) -> dict[str, list[int]]:
    if operation not in FROZEN_OPERAND_SCHEMA:
        raise ValueError("unknown archive operation: " + operation)
    children: dict[str, list[int]] = {}
    for field in FROZEN_OPERAND_SCHEMA[operation]:
        if field not in row or row[field] is None:
            continue
        value = row[field]
        values = value if isinstance(value, list) else [value]
        children[field] = [
            _index(operand, label + "." + field + "[" + str(position) + "]", source_count)
            for position, operand in enumerate(values)
        ]
    return children


def _build_occurrence_dag(
    occurrence: Mapping[str, Any],
    source_nodes: list[Any],
    rewritten_nodes: Mapping[str, int],
    prior_by_original: Mapping[int, list[str]],
    source_count: int,
    occurrence_index: int,
    boundary_collisions: set[int],
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    root_original = _index(occurrence.get("original_ref"), "occurrence original_ref", source_count)
    local_nodes: list[dict[str, Any]] = []
    local_by_original: dict[int, str] = {}
    visiting: set[int] = set()
    boundary_events: list[dict[str, Any]] = []

    def boundary(original: int) -> dict[str, Any] | None:
        if rewritten_nodes.get(str(original)) is not None and prior_by_original.get(original):
            boundary_collisions.add(original)
        return _boundary_ref(original, rewritten_nodes, prior_by_original, source_count)[0]

    def record_boundary(original: int) -> dict[str, Any] | None:
        ref, event = _boundary_ref(original, rewritten_nodes, prior_by_original, source_count)
        if event is not None:
            boundary_events.append(event)
        return ref

    root_boundary = record_boundary(root_original)
    if root_boundary is not None:
        root_status = (
            "saved_rewritten_boundary"
            if root_boundary["namespace"] == SAVED_NAMESPACE
            else "prior_fresh_boundary"
        )
        return (
            {
                "root_ref": root_boundary,
                "root_status": root_status,
                "nodes": [],
                "boundary_events": boundary_events,
                "boundary_refs": {
                    "saved": [event for event in boundary_events if event["kind"] == "saved_rewritten_nodes"],
                    "prior_fresh": [event for event in boundary_events if event["kind"] == "prior_fresh_boundary"],
                },
            },
            local_nodes,
        )

    stack: list[tuple[int, bool]] = [(root_original, False)]
    while stack:
        node_index, ready = stack.pop()
        _index(node_index, "DAG node", source_count)
        if node_index in local_by_original or boundary(node_index) is not None:
            continue
        raw_row = _mapping(source_nodes[node_index], "archive source node")
        operation = raw_row.get("operation")
        if not isinstance(operation, str) or operation not in FROZEN_OPERAND_SCHEMA:
            raise ValueError("unknown archive operation: " + repr(operation))
        operand_lists = _operand_children(raw_row, operation, source_count, "archive node " + str(node_index))
        if ready:
            if node_index not in visiting:
                raise ValueError("DAG postorder state is inconsistent")
            rewritten_row = copy.deepcopy(dict(raw_row))
            for field, children in operand_lists.items():
                refs: list[dict[str, Any]] = []
                for child in children:
                    ref = record_boundary(child)
                    if ref is None:
                        child_node = local_by_original.get(child)
                        if child_node is None:
                            raise ValueError("unresolved child in source DAG postorder")
                        ref = {"namespace": FRESH_NAMESPACE, "node": child_node}
                    refs.append(ref)
                rewritten_row[field] = refs if isinstance(raw_row[field], list) else refs[0]
            node_id = "occ%04d_n%08d" % (occurrence_index, len(local_nodes))
            local_nodes.append({
                "namespace": FRESH_NAMESPACE,
                "node_id": node_id,
                "original_ref": node_index,
                "root_context": {
                    "occurrence_id": "occ%04d" % occurrence_index,
                    "endpoint": occurrence.get("endpoint"),
                    "field": occurrence.get("field"),
                    "z_order": occurrence.get("z_order"),
                    "z_label": occurrence.get("z_label"),
                },
                "operand_schema": list(FROZEN_OPERAND_SCHEMA[operation]),
                "operator_row": rewritten_row,
            })
            local_by_original[node_index] = node_id
            visiting.remove(node_index)
            continue
        if node_index in visiting:
            raise ValueError("cycle in deferred source DAG at archive ref " + str(node_index))
        visiting.add(node_index)
        stack.append((node_index, True))
        children = [child for values in operand_lists.values() for child in values]
        for child in reversed(children):
            if boundary(child) is None and child not in local_by_original:
                stack.append((child, False))

    root_node = local_by_original.get(root_original)
    if root_node is None:
        raise ValueError("deferred source root was not constructed")
    return (
        {
            "root_ref": {"namespace": FRESH_NAMESPACE, "node": root_node},
            "root_status": "new_fresh_dag",
            "nodes": [{"namespace": FRESH_NAMESPACE, "node": node["node_id"]} for node in local_nodes],
            "boundary_events": boundary_events,
            "boundary_refs": {
                "saved": [event for event in boundary_events if event["kind"] == "saved_rewritten_nodes"],
                "prior_fresh": [event for event in boundary_events if event["kind"] == "prior_fresh_boundary"],
            },
        },
        local_nodes,
    )


def _false_flags() -> dict[str, bool]:
    return {
        "semantic_identity_proved": False,
        "ordinary_Z1_commutation_proved": False,
        "Z1_commutation_proved": False,
        "seam_identity_proved": False,
        "endpoint_equalities_proved": False,
        "P_equals_F_proved": False,
        "G_equals_F_proved": False,
        "source_values_computed": False,
        "numeric_values_computed": False,
        "integrals_evaluated": False,
        "physical_equivalence_proved": False,
        "admissible_for_numeric_use": False,
        "old_builders_called": False,
        "old_constructors_called": False,
        "old_operators_replayed": False,
        "old_error_terms_removed": False,
        "extra_ledger_terms_added": False,
        "full_goal_complete": False,
    }


def recover(root: Path) -> dict[str, Any]:
    """Return the independent, unevaluated forty-four-row source DAG packet."""

    root = Path(root)
    archive, archive_raw = _load_archive(root)
    if _mapping(archive.get("input_hashes"), "Rp archive input hashes").get(C3_ARCHIVE_NAME) != C3_ARCHIVE_SHA256:
        raise ValueError("Rp archive does not pin the saved C3 input")
    source_nodes = _list(archive.get("exact_graph_nodes"), "Rp exact_graph_nodes")
    if not source_nodes:
        raise ValueError("Rp exact_graph_nodes must be nonempty")
    source_count = len(source_nodes)
    frame = _mapping(archive.get("Rp_native_C5_pulse_source_bindings"), "Rp frame bindings")
    substitution = _canonical_mapping(frame.get("input_substitution"), "saved input_substitution", source_count)
    rewritten_nodes = _canonical_mapping(frame.get("rewritten_nodes"), "saved rewritten_nodes", source_count)
    if len(rewritten_nodes) < len(substitution):
        raise ValueError("saved rewritten_nodes is smaller than input_substitution")
    for key, target in substitution.items():
        if key in ("0", "1"):
            raise ValueError("saved input_substitution illegally rewrites reserved root")
        if rewritten_nodes.get(key) != target:
            raise ValueError("saved input_substitution disagrees with rewritten_nodes")
    substitution_digest = _substitution_digest(frame["input_substitution"])

    mismatch_result, missing, retained = _select_report_rows(
        root, archive, source_count, rewritten_nodes, substitution_digest
    )
    prior_by_original, audit_receipt, prior_root_originals, fresh_result = _load_prior_boundaries(
        root, source_count, rewritten_nodes
    )

    if any(str(original) in rewritten_nodes for original in prior_by_original):
        raise ValueError("prior fresh node unexpectedly overlaps saved rewrite boundary")
    for item in retained:
        if item["report_status"] != "reused_committed_gap_R_q_y_reference":
            continue
        side = item["endpoint"].split(".", 1)[1]
        if item["original_ref"] != EXPECTED_GAP_R_Q_Y[side][item["z_order"]]:
            raise ValueError("retained prior fresh gap root differs from published state")
        binding = _mapping(item["report_binding"], "retained prior fresh binding")
        if binding.get("namespace") != PRIOR_FRESH_NAMESPACE or binding.get("node") not in prior_by_original.get(item["original_ref"], []):
            raise ValueError("retained prior fresh root ref differs from published state")

    boundary_collisions: set[int] = set()
    occurrences: list[dict[str, Any]] = []
    fresh_nodes: list[dict[str, Any]] = []
    prior_root_reuses = 0
    for occurrence_index, source_occurrence in enumerate(missing):
        dag_packet, local_nodes = _build_occurrence_dag(
            source_occurrence,
            source_nodes,
            rewritten_nodes,
            prior_by_original,
            source_count,
            occurrence_index,
            boundary_collisions,
        )
        if dag_packet["root_status"] == "prior_fresh_boundary":
            prior_root_reuses += 1
        fresh_nodes.extend(local_nodes)
        occurrence_id = "occ%04d" % occurrence_index
        occurrences.append({
            "occurrence_id": occurrence_id,
            "endpoint": source_occurrence["endpoint"],
            "field": source_occurrence["field"],
            "z_order": source_occurrence["z_order"],
            "z_label": source_occurrence["z_label"],
            "role": source_occurrence["endpoint"] + "." + source_occurrence["field"] + "." + source_occurrence["z_label"],
            "original_ref": source_occurrence["original_ref"],
            "report_status": source_occurrence["report_status"],
            "source_ref": dict(dag_packet["root_ref"]),
            "root_ref": dict(dag_packet["root_ref"]),
            "root_status": dag_packet["root_status"],
            "dag_namespace": FRESH_NAMESPACE,
            "dag_node_refs": list(dag_packet["nodes"]),
            "dag_node_count": len(local_nodes),
            "boundary_refs": copy.deepcopy(dag_packet["boundary_refs"]),
            "boundary_events": copy.deepcopy(dag_packet["boundary_events"]),
            "source_values_computed": False,
            "integrals_evaluated": False,
            "numeric_values_computed": False,
            "semantic_identity_proved": False,
            "ordinary_Z1_commutation_proved": False,
            "seam_identity_proved": False,
        })
    if len(occurrences) != 44:
        raise ValueError("deferred occurrence output is not exactly forty-four rows")
    if len(fresh_nodes) != len({node["node_id"] for node in fresh_nodes}):
        raise ValueError("new deferred DAG node IDs are not unique")

    retained_summary = []
    for item in retained:
        binding = _mapping(item["report_binding"], "retained report binding")
        retained_summary.append({
            "endpoint": item["endpoint"],
            "field": item["field"],
            "z_order": item["z_order"],
            "z_label": item["z_label"],
            "original_ref": item["original_ref"],
            "status": item["report_status"],
            "source_ref": {key: binding[key] for key in ("namespace", "node") if key in binding},
        })
    duplicate_map = {
        str(original): {
            "candidate_nodes": list(candidates),
            "chosen_node": candidates[0],
            "selection": "deterministic_first_prior_result_order",
        }
        for original, candidates in prior_by_original.items()
        if len(candidates) > 1
    }
    flags = _false_flags()
    return {
        "schema": "frame_m15_deferred_endpoint_substitution_v1",
        "status": "unevaluated_deferred_source_dag",
        "namespace": FRESH_NAMESPACE,
        "endpoint_fields": list(M15_ENDPOINT_FIELDS),
        "ordinary_Z_rows": list(ORDINARY_Z_ROWS),
        "deferred_occurrence_count": len(occurrences),
        "retained_occurrence_counts": {
            "missing_saved_reference_deferred": len(occurrences),
            "saved_frame_reference_retained": sum(item["report_status"] == "saved_frame_reference" for item in retained),
            "prior_fresh_reference_retained": sum(item["report_status"] == "reused_committed_gap_R_q_y_reference" for item in retained),
        },
        "occurrences": occurrences,
        "retained_occurrences": retained_summary,
        "nodes": fresh_nodes,
        "fresh_node_count": len(fresh_nodes),
        "prior_fresh_root_reuses": prior_root_reuses,
        "prior_fresh_boundary_original_count": len(prior_by_original),
        "prior_fresh_duplicate_originals": duplicate_map,
        "boundary_reuse_counts": {
            "saved_rewritten_nodes": sum(len(item["boundary_refs"]["saved"]) for item in occurrences),
            "prior_fresh_nodes": sum(len(item["boundary_refs"]["prior_fresh"]) for item in occurrences),
            "prior_fresh_unique_originals": len({
                event["original_ref"]
                for item in occurrences
                for event in item["boundary_refs"]["prior_fresh"]
            }),
            "saved_prior_namespace_collisions": len(boundary_collisions),
        },
        "boundary_cache": {
            "saved_namespace": SAVED_NAMESPACE,
            "prior_fresh_namespace": PRIOR_FRESH_NAMESPACE,
            "new_namespace": FRESH_NAMESPACE,
            "saved_owner": "Rp_native_C5_pulse_source_bindings.rewritten_nodes",
            "prior_fresh_owner": "frame_gap_fresh_substitution_v1.nodes",
            "saved_semantics": "opaque authenticated archive boundary; not followed",
            "prior_fresh_semantics": "opaque audited structural boundary; not re-expanded",
            "prior_fresh_root_originals": sorted(prior_root_originals),
        },
        "provenance": {
            "mismatch_report_relative": str(MISMATCH_REPORT_RELATIVE),
            "mismatch_report_sha256": MISMATCH_REPORT_SHA256,
            "mismatch_result_sha256": MISMATCH_RESULT_SHA256,
            "archive_relative": str(ARCHIVE_RELATIVE),
            "archive_raw_sha256": hashlib.sha256(archive_raw).hexdigest(),
            "expected_archive_raw_sha256": ARCHIVE_RAW_SHA256,
            "C3_archive_name": C3_ARCHIVE_NAME,
            "C3_archive_sha256": C3_ARCHIVE_SHA256,
            "fresh_report_relative": str(FRESH_REPORT_RELATIVE),
            "fresh_report_sha256": FRESH_REPORT_SHA256,
            "fresh_result_sha256": FRESH_RESULT_SHA256,
            "fresh_structural_audit_report_sha256": AUDIT_REPORT_SHA256,
            "fresh_structural_audit_result_sha256": AUDIT_RESULT_SHA256,
            "saved_rewritten_nodes_count": len(rewritten_nodes),
            "saved_input_substitution_sha256": substitution_digest,
            "prior_audit_receipt": audit_receipt,
            "prior_fresh_schema": fresh_result.get("schema"),
        },
        "input_selection": {
            "source_frame_inputs_checked": True,
            "source_frame_endpoint_rows_checked": 64,
            "missing_status_selected_from_published_result": True,
            "report_schema": mismatch_result.get("schema"),
            "report_status": mismatch_result.get("status"),
        },
        "frozen_operand_schema": {operation: list(fields) for operation, fields in FROZEN_OPERAND_SCHEMA.items()},
        "semantic_flags": flags,
        "source_references_are_opaque": True,
        "separate_root_memo_scope": True,
        "old_twenty_roots_replayed": False,
        "old_builders_called": False,
        "old_constructors_called": False,
        "old_operators_replayed": False,
        "old_error_terms_removed": False,
        "extra_ledger_terms_added": False,
        "source_values_computed": False,
        "integrals_evaluated": False,
        "numeric_values_computed": False,
        "admissible_for_numeric_use": False,
        "full_goal_complete": False,
        "limitations": [
            "The forty-four rows remain unevaluated source references and operator metadata.",
            "Saved rewritten_nodes and prior fresh nodes are opaque boundaries and are not followed.",
            "No P/F comparison, endpoint seam identity, ordinary-Z1 commutation, or numerical claim is made.",
            "Previously saved sixteen rows and prior fresh four rows are retained.",
            "Historical error terms are not removed or replaced by this packet.",
        ],
    }



