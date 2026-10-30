"""Bounded, source-only evidence for the archived F gap Q rows.

The runner freezes and authenticates its manifest before calling
``build_q_source_rows(root, source, hashes)``.  This leaf authenticates the
same raw bytes again before decoding either JSON packet, then indexes exactly
the two top-level rows 277042 and 277055.  Row fields are copied verbatim;
integer arguments are recorded as opaque typed edge references and are never
used to index a child.  Opaque references, directly observed rows, and
semantic conclusions are separate output layers.

The synthetic public fixture is completely in-memory.  This module does not
import a producer, an earlier leaf, or a source graph helper, and it has no
executable entry point.
"""

from __future__ import annotations

import copy
import gzip
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping


SCHEMA = "f_gap_q_source_rows_v1"
F_GRAPH_NAMESPACE = "F_exact_graph"
PRIOR_COMPONENT = "f_gap_qsquare_transport_v1"
PRIOR_WRAPPER_SCHEMA = PRIOR_COMPONENT + "_result_v1"
PRIOR_INNER_STATUSES = frozenset(
    {
        "bounded_F_qsquare_structural_match_conditional",
        "bounded_F_qsquare_structural_unresolved",
    }
)
F_ARCHIVE = Path(
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_Rp_native_C5_pulse_source.json.gz"
)
F_ARCHIVE_SHA256 = "442df70e5ab4051066b241e943c8d0531067d057e774a04dcc341dadbd703ebe"
PRIOR_RESULT = Path("tmp/f_gap_qsquare_transport_v1/result.json")
PRIOR_RESULT_SHA256 = "b516e5a72290062d788114f25237ab5019abc217838a68d125de24caa1323d95"
Q_NODE_IDS = {0: 277042, 1: 277055}
PROTECTED_NODE_IDS = frozenset(
    {0, 1, 17, 33572, 25220, 263344, 263348, 263349, 262918, 205147, 263346, 33982}
)
MAX_TYPED_EDGE_REFS = 32
# Frozen operand schema copied as data from the source importer.  Only fields
# registered here can become an edge descriptor.  Other integer attributes
# remain raw row metadata and are never guessed to be node IDs.
REFERENCE_KEYS: dict[str, tuple[str, ...]] = {
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
REQUIRED_INPUTS = {
    F_ARCHIVE.as_posix(): F_ARCHIVE_SHA256,
    PRIOR_RESULT.as_posix(): PRIOR_RESULT_SHA256,
}


class FGapQSourceRowsInputError(ValueError):
    """A frozen input or bounded source-row contract is malformed."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise FGapQSourceRowsInputError(message)


def _sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def _canon(value: Any) -> bytes:
    return json.dumps(
        value, ensure_ascii=False, sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode("utf-8")


def _digest(value: Any, label: str) -> str:
    _require(isinstance(value, str) and len(value) == 64, label + " SHA malformed")
    value = value.lower()
    _require(all(ch in "0123456789abcdef" for ch in value), label + " SHA malformed")
    return value


def _manifest_hashes(source: Mapping[str, Any], hashes: Mapping[str, Any] | None) -> dict[str, str]:
    entries: dict[str, str] = {}

    def add(name: Any, digest: Any, label: str) -> None:
        _require(isinstance(name, str) and name, "manifest path missing")
        key = name.replace("\\", "/")
        value = _digest(digest, label)
        _require(key not in entries or entries[key] == value, "manifest hash conflict: " + key)
        entries[key] = value

    files = source.get("files")
    if isinstance(files, list):
        for index, item in enumerate(files):
            _require(isinstance(item, Mapping), "manifest files entry malformed")
            add(item.get("path"), item.get("sha256"), f"manifest.files[{index}]")
    for field in ("sha256", "input_hashes", "hashes"):
        table = source.get(field)
        if isinstance(table, Mapping):
            for name, digest in table.items():
                add(name, digest, f"manifest.{field}.{name}")
    if isinstance(hashes, Mapping):
        for name, digest in hashes.items():
            add(name, digest, f"caller hashes.{name}")
    _require(entries, "caller frozen manifest has no hashes")

    for required, expected in REQUIRED_INPUTS.items():
        if required in entries:
            observed = entries[required]
        else:
            aliases = [
                value for name, value in entries.items()
                if name.endswith("/" + required)
            ]
            _require(len(aliases) == 1, "required manifest hash missing: " + required)
            observed = aliases[0]
        _require(observed == expected, "required SHA changed: " + required)
    return entries


def _manifest_name(name: str) -> str:
    normalized = name.replace("\\", "/")
    absolute = normalized.startswith("/") or (
        len(normalized) >= 3 and normalized[1] == ":" and normalized[2] == "/"
    )
    if not absolute:
        return normalized
    aliases = [
        required for required in REQUIRED_INPUTS
        if normalized == required or normalized.endswith("/" + required)
    ]
    _require(len(aliases) == 1, "absolute manifest alias is outside required inputs")
    return aliases[0]


def _safe(root: Path, name: str) -> Path:
    relative = Path(name)
    _require(not relative.is_absolute() and ".." not in relative.parts, "manifest path escapes root")
    result = (root / relative).resolve()
    _require(result.is_relative_to(root.resolve()), "manifest path escapes root")
    return result


def _authenticate(
    root: Path, source: Mapping[str, Any], hashes: Mapping[str, Any] | None
) -> tuple[dict[str, bytes], dict[str, str], int]:
    """Hash every supplied byte before any JSON parse or gzip decompression."""

    entries = _manifest_hashes(source, hashes)
    raw: dict[str, bytes] = {}
    actual: dict[str, str] = {}
    for raw_name, expected in entries.items():
        name = _manifest_name(raw_name)
        data = _safe(root, name).read_bytes()
        observed = _sha(data)
        _require(observed == expected, "raw SHA changed: " + name)
        _require(name not in actual or actual[name] == observed, "manifest alias conflict: " + name)
        raw[name], actual[name] = data, observed
    for required, expected in REQUIRED_INPUTS.items():
        if required in raw:
            continue
        data = _safe(root, required).read_bytes()
        observed = _sha(data)
        _require(observed == expected, "required raw SHA changed: " + required)
        raw[required], actual[required] = data, observed
    return raw, actual, len(entries)


def _json(raw: bytes, label: str) -> Mapping[str, Any]:
    try:
        value = json.loads(raw.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise FGapQSourceRowsInputError("invalid JSON: " + label) from exc
    _require(isinstance(value, Mapping), label + " is not an object")
    return value


def _gzip_json(raw: bytes, label: str) -> Mapping[str, Any]:
    try:
        decoded = gzip.decompress(raw)
    except (OSError, EOFError) as exc:
        raise FGapQSourceRowsInputError("invalid gzip: " + label) from exc
    return _json(decoded, label)


def _opaque_ref(node: Any, role: str) -> dict[str, Any]:
    valid = type(node) is int and node >= 0
    return {
        "status": "opaque_source_ref" if valid else "unresolved_source_ref",
        "namespace": F_GRAPH_NAMESPACE,
        "node": node if valid else None,
        "source_node_id": node if valid else None,
        "source_role": role,
        "node_id_is_not_a_numeric_value": True,
        "source_value_selected": False,
        "indexed": False,
        "protected_node_indexed": False,
    }


def _node_from_ref(value: Any) -> int | None:
    """Inspect one ref object only; never follow its node ID."""

    if type(value) is int and value >= 0:
        return value
    if not isinstance(value, Mapping):
        return None
    for key in ("node", "source_node_id", "node_id"):
        item = value.get(key)
        if type(item) is int and item >= 0:
            return item
    nested = value.get("ref")
    if isinstance(nested, Mapping):
        for key in ("node", "source_node_id", "node_id"):
            item = nested.get(key)
            if type(item) is int and item >= 0:
                return item
    return None


def _edge(field: str, node: int, role: str) -> dict[str, Any]:
    ref = _opaque_ref(node, role + "." + field)
    return {
        "field": field,
        "typed_edge": True,
        "ref": ref,
        "child_indexed": False,
        "edge_followed": False,
        "protected_target": node in PROTECTED_NODE_IDS,
    }


def _typed_edge_refs(row: Mapping[str, Any], node: int) -> dict[str, Any]:
    """Extract registered operand fields only, without child access.

    An operation absent from ``REFERENCE_KEYS`` is explicitly unknown.  Its
    integer attributes, including an ``arguments`` field, remain raw metadata
    and are not guessed to be graph references.  For known operations, only
    the registered fields can produce opaque edge descriptors; this keeps
    literal fields such as ``exact_rational.value`` out of the edge set.
    """

    operation = row.get("operation")
    fields = REFERENCE_KEYS.get(operation) if isinstance(operation, str) else None
    unknown = fields is None
    edges: list[dict[str, Any]] = []
    skipped = 0
    unresolved: list[str] = []
    if fields is not None:
        for field in fields:
            if field not in row:
                continue
            value = row.get(field)
            if field == "arguments":
                if not isinstance(value, list):
                    unresolved.append(field)
                    continue
                values = ((f"{field}[{index}]", item) for index, item in enumerate(value))
            else:
                values = ((field, value),)
            for edge_field, candidate in values:
                if len(edges) >= MAX_TYPED_EDGE_REFS:
                    skipped += 1
                    continue
                child = _node_from_ref(candidate)
                if child is None:
                    unresolved.append(edge_field)
                    continue
                edges.append(_edge(edge_field, child, f"F.Q[{node}].row"))

    registered = set(fields or ())
    ignored_integer_attributes = sorted(
        str(key)
        for key, value in row.items()
        if key != "operation"
        and key not in registered
        and (
            type(value) is int
            or isinstance(value, list) and any(type(item) is int for item in value)
        )
    )
    return {
        "operation": operation,
        "operation_known": not unknown,
        "unknown_operation": unknown,
        "registered_reference_fields": list(fields or ()),
        "refs": edges,
        "count": len(edges),
        "max": MAX_TYPED_EDGE_REFS,
        "skipped_or_unresolved": skipped + len(unresolved),
        "unresolved_registered_fields": unresolved,
        "ignored_integer_attributes": ignored_integer_attributes,
        "children_indexed": False,
        "edges_followed": False,
    }


def _semantic_false_flags() -> dict[str, bool]:
    return {
        "source_identity_proved": False,
        "same_source_function_proved": False,
        "ordinary_Z_derivative_proved": False,
        "ordinary_Z_transport_proved": False,
        "Q1_is_partial_Z_Q0_proved": False,
        "capture_free_substitution_proved": False,
        "endpoint_transport_proved": False,
        "source_values_selected": False,
        "source_values_computed": False,
        "integrals_evaluated": False,
        "numeric_use_admissible": False,
        "error_transport_proved": False,
        "error_removal": False,
        "closure_proved": False,
    }


def _observed_row(row: Mapping[str, Any], node: int, role: str) -> dict[str, Any]:
    copied = copy.deepcopy(dict(row))
    arguments = copied.get("arguments")
    edge_packet = _typed_edge_refs(copied, node)
    return {
        "status": "direct_observed_row",
        "role": role,
        "node": node,
        "namespace": F_GRAPH_NAMESPACE,
        "operation": copied.get("operation"),
        "field_count": len(copied),
        "argument_count": len(arguments) if isinstance(arguments, list) else None,
        "typed_edge_ref_count": edge_packet["count"],
        "typed_edge_refs": edge_packet["refs"],
        "typed_edge_scope": edge_packet,
        "row": copied,
        "raw_fields": copy.deepcopy(copied),
        "row_sha256": _sha(_canon(copied)),
        "source_value_selected": False,
        "children_indexed": False,
        "child_rows_read": 0,
        "semantic_identity_proved": False,
    }


def _direct_row(nodes: Any, node: int, order: int) -> dict[str, Any]:
    _require(node not in PROTECTED_NODE_IDS, "requested Q node is protected")
    _require(isinstance(nodes, list), "F exact_graph_nodes missing")
    _require(node < len(nodes), "requested Q node outside archive")
    row = nodes[node]
    _require(isinstance(row, Mapping), "requested Q row malformed")
    return {
        "opaque_ref": _opaque_ref(node, f"F_archive.leading_fields.gap.Q[{order}]"),
        "direct_observed_row": _observed_row(
            row, node, f"F_archive.leading_fields.gap.Q[{order}]"
        ),
        "semantic_proof": _semantic_false_flags(),
    }


def _ref_at(value: Any, role: str) -> dict[str, Any] | None:
    if not isinstance(value, Mapping):
        return None
    node = _node_from_ref(value)
    if node is None:
        return None
    return {"prior_ref": copy.deepcopy(dict(value)), "opaque_ref": _opaque_ref(node, role)}


def _prior_qsquare(packet: Mapping[str, Any]) -> Mapping[str, Any]:
    """Accept only the successful, fixed prior runner envelope."""

    _require(packet.get("schema") == PRIOR_WRAPPER_SCHEMA, "prior wrapper schema changed")
    _require(packet.get("component") == PRIOR_COMPONENT, "prior wrapper component changed")
    _require(packet.get("status") == "returned", "prior wrapper status is not returned")
    _require("error" in packet and packet.get("error") is None, "prior wrapper contains an error")
    result = packet.get("result")
    _require(isinstance(result, Mapping), "prior wrapper result missing")
    candidate = result.get("qsquare_transport")
    _require(isinstance(candidate, Mapping), "prior wrapper qsquare_transport missing")
    _require(candidate.get("schema") == PRIOR_COMPONENT, "prior inner schema changed")
    _require(candidate.get("component") == PRIOR_COMPONENT, "prior inner component changed")
    _require(
        candidate.get("status") in PRIOR_INNER_STATUSES,
        "prior inner status is not a recognized structural result",
    )
    return candidate


def _prior_binding(packet: Mapping[str, Any]) -> dict[str, Any]:
    qsquare = _prior_qsquare(packet)
    issues: list[str] = []
    prior_scope = qsquare.get("prior_scope")
    leading_prior = prior_scope.get("leading_Q_prior") if isinstance(prior_scope, Mapping) else None
    snapshot = qsquare.get("F_snapshot")
    leading_snapshot = snapshot.get("leading_Q_roles") if isinstance(snapshot, Mapping) else None
    prior_refs: dict[str, Any] = {}
    snapshot_refs: dict[str, Any] = {}
    for order, expected in Q_NODE_IDS.items():
        key = "Q" + str(order)
        prior_ref = _ref_at(
            leading_prior.get(key) if isinstance(leading_prior, Mapping) else None,
            f"prior.leading_Q_prior.{key}",
        )
        snapshot_ref = _ref_at(
            leading_snapshot.get(key) if isinstance(leading_snapshot, Mapping) else None,
            f"prior.F_snapshot.leading_Q_roles.{key}",
        )
        if prior_ref is None:
            issues.append(key + "_prior_ref_missing")
        elif prior_ref["opaque_ref"]["node"] != expected:
            issues.append(key + "_prior_ref_mismatch")
        if snapshot_ref is None:
            issues.append(key + "_snapshot_ref_missing")
        elif snapshot_ref["opaque_ref"]["node"] != expected:
            issues.append(key + "_snapshot_ref_mismatch")
        prior_refs[key], snapshot_refs[key] = prior_ref, snapshot_ref
    namespace = snapshot.get("namespace") if isinstance(snapshot, Mapping) else None
    if namespace != F_GRAPH_NAMESPACE:
        issues.append("F_namespace_mismatch")
    if qsquare.get("same_F_snapshot_role_ref_correspondence") is not True:
        issues.append("prior_F_snapshot_correspondence_not_true")
    # Preserve the prior operation-shape conclusion as a scalar.  Do not copy
    # the old Q-squared shape/direct-row tables into this packet: the runner's
    # recursive scope audit requires every direct-row table in this component
    # to name only the two newly observed Q rows.
    shape = {"actual_row_shape_matched": qsquare.get("actual_row_shape_matched") is True}
    if not shape["actual_row_shape_matched"]:
        issues.append("prior_operation_structure_unmatched")
    return {
        "status": "matched" if not issues else "unresolved",
        "path": PRIOR_RESULT.as_posix(),
        "sha256": PRIOR_RESULT_SHA256,
        "schema": packet.get("schema"),
        "component": packet.get("component"),
        "wrapper_status": packet.get("status"),
        "wrapper_error": None,
        "inner_schema": qsquare.get("schema"),
        "inner_component": qsquare.get("component"),
        "inner_status": qsquare.get("status"),
        "Q0": prior_refs["Q0"],
        "Q1": prior_refs["Q1"],
        "F_snapshot_Q0": snapshot_refs["Q0"],
        "F_snapshot_Q1": snapshot_refs["Q1"],
        "namespace": namespace,
        "issues": sorted(set(issues)),
        "operation_structure_evidence": {
            "status": "observed_in_prior_result"
            if shape["actual_row_shape_matched"]
            else "unresolved_in_prior_result",
            **shape,
            "semantic_identity_proved": False,
        },
        "source_values_selected": False,
    }


def _row_summary(records: Mapping[str, Any]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key in ("Q0", "Q1"):
        packet = records.get(key)
        row = packet.get("direct_observed_row") if isinstance(packet, Mapping) else None
        result[key] = {
            "operation": row.get("operation") if isinstance(row, Mapping) else None,
            "field_count": row.get("field_count") if isinstance(row, Mapping) else None,
            "argument_count": row.get("argument_count") if isinstance(row, Mapping) else None,
            "typed_edge_ref_count": row.get("typed_edge_ref_count") if isinstance(row, Mapping) else None,
            "row_sha256": row.get("row_sha256") if isinstance(row, Mapping) else None,
        }
    return result


def verify_synthetic_fixture() -> dict[str, Any]:
    """Check the evidence layers using only caller-owned in-memory rows."""

    q0 = _observed_row(
        {"operation": "product", "arguments": [41, 41], "literal": 99},
        401,
        "synthetic.Q0",
    )
    q1 = _observed_row(
        {"operation": "sum", "arguments": [402, 403], "typed_edges": []},
        404,
        "synthetic.Q1",
    )
    quotient = _observed_row(
        {"operation": "positive_quotient", "numerator": 17, "denominator": 19, "literal": 23},
        405,
        "synthetic.quotient",
    )
    unknown = _observed_row(
        {"operation": "future_unknown", "arguments": [501], "literal": 502},
        406,
        "synthetic.unknown",
    )
    synthetic_inner = {
        "schema": PRIOR_COMPONENT,
        "component": PRIOR_COMPONENT,
        "status": "bounded_F_qsquare_structural_match_conditional",
    }
    synthetic_wrapper = {
        "schema": PRIOR_WRAPPER_SCHEMA,
        "component": PRIOR_COMPONENT,
        "status": "returned",
        "error": None,
        "result": {"qsquare_transport": synthetic_inner},
    }
    wrapper_success = False
    try:
        wrapper_success = _prior_qsquare(synthetic_wrapper) is synthetic_inner
    except FGapQSourceRowsInputError:
        wrapper_success = False
    synthetic_failure = copy.deepcopy(synthetic_wrapper)
    synthetic_failure["status"] = "failed"
    synthetic_failure["error"] = {"type": "SyntheticFailure"}
    wrapper_failure_rejected = False
    try:
        _prior_qsquare(synthetic_failure)
    except FGapQSourceRowsInputError:
        wrapper_failure_rejected = True
    checks = {
        "opaque_ref_keeps_namespace_and_node": (
            _opaque_ref(401, "synthetic")["namespace"] == F_GRAPH_NAMESPACE
            and _opaque_ref(401, "synthetic")["node"] == 401
            and _opaque_ref(401, "synthetic")["source_value_selected"] is False
        ),
        "raw_fields_are_preserved": q0["raw_fields"] == {
            "operation": "product",
            "arguments": [41, 41],
            "literal": 99,
        },
        "typed_refs_are_finite_and_unfollowed": (
            q1["typed_edge_ref_count"] == 2
            and q1["typed_edge_scope"]["edges_followed"] is False
            and q1["typed_edge_scope"]["children_indexed"] is False
        ),
        "quotient_registered_fields_are_extracted": (
            quotient["typed_edge_ref_count"] == 2
            and quotient["typed_edge_scope"]["registered_reference_fields"]
            == ["numerator", "denominator"]
            and quotient["typed_edge_scope"]["unknown_operation"] is False
        ),
        "unknown_operation_is_not_guessed": (
            unknown["typed_edge_ref_count"] == 0
            and unknown["typed_edge_scope"]["unknown_operation"] is True
            and unknown["typed_edge_scope"]["operation_known"] is False
        ),
        "literal_integer_attributes_are_ignored": (
            q0["typed_edge_scope"]["ignored_integer_attributes"] == ["literal"]
            and quotient["typed_edge_scope"]["ignored_integer_attributes"] == ["literal"]
            and "arguments" in unknown["typed_edge_scope"]["ignored_integer_attributes"]
        ),
        "successful_prior_wrapper_is_accepted": wrapper_success,
        "failed_prior_wrapper_is_rejected": wrapper_failure_rejected,
        "row_sha_is_structural": len(q0["row_sha256"]) == 64 and q0["semantic_identity_proved"] is False,
        "semantic_flags_stay_false": all(not value for value in _semantic_false_flags().values()),
    }
    return {
        "schema": SCHEMA + "_synthetic_fixture_v1",
        "component": SCHEMA,
        "status": "synthetic_checks_passed" if all(checks.values()) else "synthetic_checks_failed",
        "synthetic_only": True,
        "checks": checks,
        "check_count": len(checks),
        "fixtures": {"Q0": q0, "Q1": q1, "quotient": quotient, "unknown": unknown},
        "semantic_proof": _semantic_false_flags(),
        "real_inputs_read": False,
        "archive_decoded": False,
        "source_values_selected": False,
    }


def build_q_source_rows(
    root: str | Path = ".",
    source: Mapping[str, Any] | None = None,
    hashes: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Read exactly Q0 and Q1 from the authenticated F archive."""

    _require(isinstance(source, Mapping), "caller frozen manifest is required")
    raw, actual, manifest_count = _authenticate(Path(root).resolve(), source, hashes)
    prior_packet = _json(raw[PRIOR_RESULT.as_posix()], PRIOR_RESULT.as_posix())
    prior = _prior_binding(prior_packet)
    archive = _gzip_json(raw[F_ARCHIVE.as_posix()], F_ARCHIVE.as_posix())
    nodes = archive.get("exact_graph_nodes")
    _require(isinstance(nodes, list), "F exact_graph_nodes missing")
    records = {
        "Q0": _direct_row(nodes, Q_NODE_IDS[0], 0),
        "Q1": _direct_row(nodes, Q_NODE_IDS[1], 1),
    }
    indexed = [Q_NODE_IDS[0], Q_NODE_IDS[1]]
    summary = _row_summary(records)
    semantic = _semantic_false_flags()
    return {
        "schema": SCHEMA,
        "component": SCHEMA,
        "status": "bounded_F_Q_source_rows_observed"
        if prior["status"] == "matched"
        else "bounded_F_Q_source_rows_observed_prior_unresolved",
        "source_only": True,
        "runner_interface": {
            "entry": "build_q_source_rows",
            "arguments": ["root", "source", "hashes"],
            "one_shot": True,
            "coordination": "/root/q_source_runner",
        },
        "raw_hash_authentication": {
            "all_manifest_raw_hashes_checked_before_decode": True,
            "manifest_entry_count": manifest_count,
            "required_inputs": {
                name: {"path": name, "expected_sha256": digest, "actual_sha256": actual[name]}
                for name, digest in REQUIRED_INPUTS.items()
            },
        },
        "prior_result_binding": prior,
        "F_snapshot": {
            "archive_path": F_ARCHIVE.as_posix(),
            "namespace": F_GRAPH_NAMESPACE,
            "Q0": _opaque_ref(Q_NODE_IDS[0], "F_archive.leading_fields.gap.Q[0]"),
            "Q1": _opaque_ref(Q_NODE_IDS[1], "F_archive.leading_fields.gap.Q[1]"),
            "whole_json_decoded": True,
            "source_identity_proved": False,
        },
        "direct_graph_scope": {
            "requested_Q_node_ids": indexed,
            "indexed_rows": indexed,
            "indexed_row_count": 2,
            "maximum_indexed_rows": 2,
            "Q_rows_are_direct_observations": True,
            "Q_references_in_edges_are_opaque": True,
            "children_not_indexed": True,
            "child_rows_read": 0,
            "graph_edges_followed": False,
            "source_graph_traversed": False,
            "ancestor_replay": False,
            "numeric_evaluation": False,
            "protected_node_ids": sorted(PROTECTED_NODE_IDS),
            "protected_nodes_indexed": False,
            "forbidden_nodes_indexed": False,
        },
        "rows": records,
        "row_operation_summary": summary,
        "indexed_row_count": 2,
        "operation_structure_evidence": {
            "prior_result": copy.deepcopy(prior["operation_structure_evidence"]),
            "direct_rows": {"Q0": summary["Q0"], "Q1": summary["Q1"], "semantic_identity_proved": False},
        },
        "semantic_proof": semantic,
        "source_gates": {
            "prior_operation_structure_evidence": bool(
                prior["operation_structure_evidence"]["actual_row_shape_matched"]
            ),
            **semantic,
        },
        "scope_statement": (
            "Exactly top-level F graph rows 277042 and 277055 were directly observed "
            "after raw hash authentication. Raw fields, operation metadata, finite "
            "typed edge references, and canonical row SHA-256 values are retained. "
            "Edge references are opaque; no child, ancestor, source value, derivative, "
            "integral, error, or closure semantic is established."
        ),
    }


build_f_gap_q_source_rows = build_q_source_rows
build_q_source_row_evidence = build_q_source_rows

__all__ = [
    "FGapQSourceRowsInputError",
    "build_q_source_rows",
    "build_f_gap_q_source_rows",
    "build_q_source_row_evidence",
    "verify_synthetic_fixture",
]
