"""Narrow, source-backed closure audit for the Q4 frontier rows.

This leaf is intentionally smaller than a graph walker.  The caller supplies
an explicit mapping whose only admissible keys are observed one at a time for
the ten frontier ids below.  A caller supplied adapter may translate an
unknown source-row schema into the small canonical record described by
``audit``.  The adapter is not a source authenticator: it is only a typed
view of the row that the caller has already selected.

The leaf records row metadata, declared edges, and a source-binding identity.
It never dereferences an edge, even when that edge names another selected
frontier.  Protected and other external ids are copied as opaque unresolved
edges.  Provider-reference or source-text equality is deliberately ignored.

``conditionally_provable`` means that the supplied row claims complete,
scope-aware, authenticated occurrence/free-variable evidence and complete
declared edges.  The claim remains conditional on the caller's source
authentication; this module never proves source bytes, a live provider, or a
whole-graph binder theorem.  Missing evidence is ``unresolved``.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from typing import Any


SCHEMA = "original_q4_frontier_closure_audit_v1"
RESULT_SCHEMA = SCHEMA

# Keep this order stable: it is the order used by the source-backed handoff.
TARGET_NODE_IDS = (
    35559,
    35674,
    35675,
    35676,
    35677,
    35678,
    35679,
    35680,
    35558,
    20,
)
TARGET_IDS = TARGET_NODE_IDS

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
        33982,
    }
)
PROTECTED_IDS = PROTECTED_NODE_IDS

_MISSING = object()
_IGNORED_IDENTITY_KEYS = frozenset(
    {"provider_refs", "provider_ref", "text", "source_text"}
)


def _shallow_mapping(value: Any) -> dict[Any, Any] | None:
    """Copy only top-level mapping entries; never recurse into row data."""

    if not isinstance(value, Mapping):
        return None
    try:
        return {key: item for key, item in value.items()}
    except Exception:
        return None


def _exact_true(value: Any) -> bool:
    return value is True


def _field(mapping: Mapping[Any, Any] | None, key: str, default: Any = _MISSING) -> Any:
    if mapping is None:
        return default
    try:
        return mapping[key] if key in mapping else default
    except Exception:
        return default


def _error(code: str, target_id: int | None = None, detail: str | None = None) -> dict[str, Any]:
    item: dict[str, Any] = {"code": code}
    if target_id is not None:
        item["node_id"] = target_id
    if detail is not None:
        item["detail"] = detail
    return item


def _empty_scope() -> dict[str, Any]:
    return {
        "row_mapping_only": True,
        "adapter_required_for_unknown_schema": True,
        "source_rows_loaded": False,
        "source_graph_traversed": False,
        "declared_edges_traversed": False,
        "source_values_evaluated": False,
        "numeric_field_value_computed": False,
        "protected_ids_indexed": False,
        "non_target_ids_indexed": False,
        "provider_refs_or_text_used_as_proof": False,
        "source_authentication_proved": False,
        "conditional_on_caller_source_authentication": True,
        "closure_provable": False,
        "old_errors_removed": False,
    }


def _base_report() -> dict[str, Any]:
    return {
        "schema": SCHEMA,
        "status": "unresolved",
        "passed": False,
        "target_ids": list(TARGET_NODE_IDS),
        "rows": [],
        "nodes": [],
        "frontier": [],
        "errors": [],
        "access": {
            "requested_ids": list(TARGET_NODE_IDS),
            "indexed_ids": [],
            "non_target_ids_indexed": [],
            "protected_ids_indexed": [],
        },
        "counts": {
            "target_ids": len(TARGET_NODE_IDS),
            "rows_present": 0,
            "rows_missing": len(TARGET_NODE_IDS),
            "rows_conditionally_provable": 0,
            "rows_unresolved": 0,
            "rows_mismatch": 0,
            "declared_edges": 0,
            "opaque_external_edges": 0,
            "protected_edges": 0,
        },
        "scope": _empty_scope(),
    }


def _canonical_row(target_id: int, value: Mapping[Any, Any]) -> dict[str, Any]:
    """Read canonical fields without interpreting arbitrary integer metadata."""

    node_id = _field(value, "node_id", target_id)
    metadata_value = _field(value, "metadata", _MISSING)
    metadata = _shallow_mapping(metadata_value)
    if metadata is None:
        # A direct raw source row is retained as metadata.  No raw integer
        # field is promoted to an edge unless the explicit canonical
        # ``declared_edges`` field is present.
        metadata = _shallow_mapping(value) or {}

    declared_edges = _field(value, "declared_edges", _MISSING)
    binding = _field(value, "source_binding_identity", _MISSING)
    occurrence = _field(value, "occurrence_evidence", _MISSING)
    free_variables = _field(value, "free_variable_evidence", _MISSING)
    source_auth = _field(value, "source_authentication", _MISSING)

    return {
        "node_id": node_id,
        "metadata": metadata,
        "declared_edges": declared_edges,
        "source_binding_identity": _shallow_mapping(binding),
        "occurrence_evidence": _shallow_mapping(occurrence),
        "free_variable_evidence": _shallow_mapping(free_variables),
        "source_authentication": _shallow_mapping(source_auth),
        "declared_edges_present": declared_edges is not _MISSING,
    }


def _edge_records(
    target_id: int, raw_edges: Any
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]], bool]:
    """Normalize explicit edge descriptors without looking up any child row."""

    records: list[dict[str, Any]] = []
    frontiers: list[dict[str, Any]] = []
    errors: list[dict[str, Any]] = []
    malformed = False

    if raw_edges is _MISSING:
        return records, frontiers, [_error("declared_edges_missing", target_id)], False
    if not isinstance(raw_edges, Sequence) or isinstance(raw_edges, (str, bytes, bytearray)):
        return records, frontiers, [
            _error("declared_edges_shape", target_id, "declared_edges must be a sequence")
        ], True

    for index, raw_edge in enumerate(raw_edges):
        edge = _shallow_mapping(raw_edge)
        if edge is None and type(raw_edge) is int:
            edge = {"node_id": raw_edge}
        if edge is None:
            malformed = True
            error = _error("declared_edge_shape", target_id, f"edge {index} is not a mapping or integer")
            errors.append(error)
            frontiers.append(
                {
                    "owner_node_id": target_id,
                    "edge_index": index,
                    "edge": raw_edge,
                    "resolution": "unresolved",
                    "traversed": False,
                    "reason": "opaque_edge_shape",
                }
            )
            continue

        edge_id = _field(edge, "node_id", _MISSING)
        item: dict[str, Any] = {
            "owner_node_id": target_id,
            "edge_index": index,
            "node_id": edge_id if edge_id is not _MISSING else None,
            "edge": edge,
            "traversed": False,
        }
        if type(edge_id) is not int:
            malformed = True
            item.update({"resolution": "unresolved", "reason": "edge_node_id_missing_or_invalid"})
            errors.append(_error("declared_edge_node_id_shape", target_id, f"edge {index} has no integer node_id"))
            frontiers.append(item)
        elif edge_id in PROTECTED_NODE_IDS:
            item.update({"resolution": "unresolved", "protected": True, "reason": "protected_opaque_edge"})
            frontiers.append(item)
        elif edge_id not in TARGET_NODE_IDS:
            item.update({"resolution": "unresolved", "external": True, "reason": "external_opaque_edge"})
            frontiers.append(item)
        else:
            # The selected child is still not dereferenced.  It is merely a
            # declared member of this audit's explicit target set.
            item.update({"resolution": "declared_in_scope", "in_scope": True})
        records.append(item)

    return records, frontiers, errors, malformed


def _row_audit(target_id: int, record: Mapping[Any, Any], adapter_used: bool) -> dict[str, Any]:
    """Audit one already selected row; no child row is ever read."""

    item: dict[str, Any] = {
        "node_id": target_id,
        "adapter_used": adapter_used,
        "status": "unresolved",
        "closure_provable": False,
        "closure_proof_is_conditional": True,
        "raw_row": _shallow_mapping(record) or {},
        "metadata": {},
        "declared_edges": [],
        "frontier": [],
        "errors": [],
        "source_binding_identity": None,
        "occurrence_evidence": None,
        "free_variable_evidence": None,
        "source_authentication": None,
        "checks": {},
    }

    canonical = _canonical_row(target_id, record)
    item["metadata"] = canonical["metadata"]
    item["source_binding_identity"] = canonical["source_binding_identity"]
    item["occurrence_evidence"] = canonical["occurrence_evidence"]
    item["free_variable_evidence"] = canonical["free_variable_evidence"]
    item["source_authentication"] = canonical["source_authentication"]

    node_id = canonical["node_id"]
    node_id_ok = type(node_id) is int and node_id == target_id
    item["checks"]["node_id_matches_requested_target"] = node_id_ok
    if not node_id_ok:
        item["errors"].append(_error("row_node_id_mismatch", target_id))

    edges, frontiers, edge_errors, edge_malformed = _edge_records(
        target_id, canonical["declared_edges"]
    )
    item["declared_edges"] = edges
    item["frontier"].extend(frontiers)
    item["errors"].extend(edge_errors)
    item["checks"]["declared_edges_present"] = canonical["declared_edges_present"]
    item["checks"]["declared_edges_complete"] = (
        _field(record, "declared_edges_complete") is True
    )
    item["checks"]["no_external_or_protected_edges"] = not any(
        edge.get("resolution") == "unresolved" for edge in frontiers
    )

    binding = canonical["source_binding_identity"]
    binding_id = _field(binding, "binder_node_id")
    binding_relation = _field(binding, "binding_relation")
    binding_identity_claim = _exact_true(_field(binding, "identity_proved"))
    binding_node_in_scope = (
        type(binding_id) is int
        and binding_id in TARGET_NODE_IDS
        and binding_id not in PROTECTED_NODE_IDS
    )
    binding_shape = (
        binding is not None
        and binding_node_in_scope
        and isinstance(binding_relation, str)
        and bool(binding_relation)
        and binding_identity_claim
    )
    item["checks"]["source_binding_node_is_in_scope"] = binding_node_in_scope
    item["checks"]["source_binding_identity_explicit"] = binding_shape
    item["checks"]["source_binding_identity_is_not_provider_or_text_equality"] = bool(
        binding is not None
        and not (
            set(binding).issubset(_IGNORED_IDENTITY_KEYS)
            and bool(binding)
        )
    )

    occurrence = canonical["occurrence_evidence"]
    occurrence_complete = (
        occurrence is not None
        and _exact_true(_field(occurrence, "complete"))
        and _exact_true(_field(occurrence, "scope_aware"))
        and _exact_true(_field(occurrence, "authenticated"))
        and _exact_true(_field(occurrence, "binding_identity_proved"))
    )
    item["checks"]["occurrence_evidence_complete_authenticated"] = occurrence_complete

    free_variables = canonical["free_variable_evidence"]
    free_complete = (
        free_variables is not None
        and _exact_true(_field(free_variables, "complete"))
        and _exact_true(_field(free_variables, "scope_aware"))
        and _exact_true(_field(free_variables, "authenticated"))
        and _exact_true(_field(free_variables, "free_variable_set_complete"))
        and _exact_true(_field(free_variables, "target_binder_absent"))
    )
    item["checks"]["free_variable_evidence_complete_authenticated"] = free_complete

    source_auth = canonical["source_authentication"]
    source_auth_claim = (
        source_auth is not None
        and _exact_true(_field(source_auth, "authenticated"))
        and isinstance(_field(source_auth, "basis"), str)
        and bool(_field(source_auth, "basis"))
    )
    item["checks"]["caller_source_authentication_claim_present"] = source_auth_claim

    all_claims = (
        node_id_ok
        and canonical["declared_edges_present"]
        and item["checks"]["declared_edges_complete"]
        and item["checks"]["no_external_or_protected_edges"]
        and _exact_true(_field(record, "scope_aware"))
        and binding_shape
        and occurrence_complete
        and free_complete
        and source_auth_claim
        and not edge_malformed
    )
    item["checks"]["row_scope_aware"] = _exact_true(_field(record, "scope_aware"))
    item["checks"]["opaque_edges_not_traversed"] = all(
        edge.get("traversed") is False for edge in edges
    )
    item["closure_provable"] = bool(all_claims)
    item["status"] = "conditionally_provable" if all_claims else "unresolved"
    if any(error.get("code") in {"row_node_id_mismatch", "declared_edges_shape", "declared_edge_node_id_shape"} for error in item["errors"]):
        item["status"] = "mismatch"
    return item


def audit(
    rows: Mapping[int, Any],
    adapter: Callable[[int, Any], Mapping[Any, Any]] | None = None,
    *,
    row_adapter: Callable[[int, Any], Mapping[Any, Any]] | None = None,
) -> dict[str, Any]:
    """Audit exactly the supplied Q4 frontier rows.

    ``rows`` is a caller-owned mapping.  Only ``rows[target_id]`` for each
    member of :data:`TARGET_NODE_IDS` is requested; mapping keys outside that
    set are never enumerated or indexed.  If the source schema is not already
    canonical, pass ``adapter(target_id, raw_row)``.  The adapter should return
    a mapping with ``node_id``, ``metadata``, ``declared_edges``,
    ``declared_edges_complete``, ``scope_aware``,
    ``source_binding_identity``, ``occurrence_evidence``,
    ``free_variable_evidence``, and ``source_authentication``.  Every evidence
    flag is treated as a caller claim and yields only a conditional result.
    """

    report = _base_report()
    if not isinstance(rows, Mapping):
        report["status"] = "mismatch"
        report["errors"] = [_error("row_mapping_shape", detail="rows must be a mapping")]
        report["scope"]["source_rows_loaded"] = False
        return report

    if adapter is not None and row_adapter is not None:
        report["status"] = "mismatch"
        report["errors"] = [_error("multiple_adapters_supplied")]
        return report
    if row_adapter is not None:
        adapter = row_adapter
    if adapter is not None and not callable(adapter):
        report["status"] = "mismatch"
        report["errors"] = [_error("adapter_shape", detail="adapter must be callable")]
        return report

    report["scope"]["source_rows_loaded"] = True
    report["scope"]["adapter_used"] = adapter is not None
    row_reports: list[dict[str, Any]] = []
    aggregate_frontier: list[dict[str, Any]] = []
    aggregate_errors: list[dict[str, Any]] = []
    mismatch_count = 0

    # This is deliberately a fixed loop.  Do not replace it with a walk over
    # rows.keys(), and do not ask the adapter to resolve any edge child.
    for target_id in TARGET_NODE_IDS:
        try:
            raw_row = rows[target_id]
        except KeyError:
            aggregate_errors.append(_error("target_row_missing", target_id))
            continue
        except Exception as exc:
            mismatch_count += 1
            aggregate_errors.append(
                _error("target_row_read_error", target_id, type(exc).__name__)
            )
            continue

        report["access"]["indexed_ids"].append(target_id)
        report["counts"]["rows_present"] += 1
        if adapter is None:
            selected = raw_row
        else:
            try:
                selected = adapter(target_id, raw_row)
            except Exception as exc:
                mismatch_count += 1
                aggregate_errors.append(
                    _error("adapter_error", target_id, type(exc).__name__)
                )
                continue
        if not isinstance(selected, Mapping):
            mismatch_count += 1
            aggregate_errors.append(
                _error("adapter_result_shape", target_id, "adapter must return a mapping")
            )
            continue

        row_report = _row_audit(target_id, selected, adapter is not None)
        row_reports.append(row_report)
        aggregate_frontier.extend(row_report["frontier"])
        aggregate_errors.extend(row_report["errors"])
        if row_report["status"] == "mismatch":
            mismatch_count += 1

        report["counts"]["declared_edges"] += len(row_report["declared_edges"])
        report["counts"]["opaque_external_edges"] += sum(
            edge.get("resolution") == "unresolved"
            and not edge.get("protected", False)
            for edge in row_report["frontier"]
        )
        report["counts"]["protected_edges"] += sum(
            edge.get("protected", False) for edge in row_report["frontier"]
        )

    report["rows"] = row_reports
    report["nodes"] = row_reports
    report["frontier"] = aggregate_frontier
    report["errors"] = aggregate_errors
    report["counts"]["rows_missing"] = len(TARGET_NODE_IDS) - report["counts"]["rows_present"]
    conditional = sum(row["status"] == "conditionally_provable" for row in row_reports)
    unresolved = sum(row["status"] == "unresolved" for row in row_reports)
    report["counts"]["rows_conditionally_provable"] = conditional
    report["counts"]["rows_unresolved"] = unresolved + report["counts"]["rows_missing"]
    report["counts"]["rows_mismatch"] = mismatch_count

    complete = (
        report["counts"]["rows_present"] == len(TARGET_NODE_IDS)
        and conditional == len(TARGET_NODE_IDS)
        and mismatch_count == 0
    )
    if mismatch_count:
        report["status"] = "mismatch"
    elif complete:
        report["status"] = "conditionally_provable"
        report["passed"] = True
    else:
        report["status"] = "unresolved"

    # The fixed target loop is the only indexing performed by this module.
    report["scope"]["indexed_target_ids"] = list(report["access"]["indexed_ids"])
    report["scope"]["non_target_ids_indexed"] = False
    report["scope"]["protected_ids_indexed"] = False
    report["scope"]["opaque_external_edges_unresolved"] = any(
        edge.get("resolution") == "unresolved" for edge in aggregate_frontier
    )
    report["scope"]["closure_provable"] = complete
    return report


# The alias keeps the leaf easy to call from generic source-audit runners.
analyze = audit


__all__ = [
    "SCHEMA",
    "RESULT_SCHEMA",
    "TARGET_NODE_IDS",
    "TARGET_IDS",
    "PROTECTED_NODE_IDS",
    "PROTECTED_IDS",
    "audit",
    "analyze",
]
