"""Bounded source-only evidence for the historical F gap-radius spine.

The consumer authenticates the caller's frozen manifest before decoding the F
and old C3 gzip archives.  It directly indexes only rows 33982, 33981, and
33830 on each side.  Every child reference is an opaque frontier; protected
nodes are never explicitly indexed or dereferenced.  No graph edge, old leaf,
numerical value, integral, or result reducer is called here.  Whole-archive
JSON decoding may still contain protected-node bytes.

Equal handles, equal row hashes, equal source-family metadata, and finite
ordinary-Z rows are metadata facts only.  They never promote the F owner,
capture freedom, ordinary-Z transport, or function equality gates.
"""

from __future__ import annotations

import copy
import gzip
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping


SCHEMA = "f_gap_radius_spine_binding_v1"
COMPONENT = SCHEMA
F_ARCHIVE = Path(
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_Rp_native_C5_pulse_source.json.gz"
)
OLD_ARCHIVE = Path(
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_pulse_pressure_C3_functions.json.gz"
)
DEFERRED_RESULT = Path("tmp/historical_frame_gap_deferred_roles_v1/result.json")
LOOKUP_RESULT = Path("tmp/historical_frame_gap_saved_rewrite_lookup_v1/result.json")
F_ARCHIVE_SHA256 = "442df70e5ab4051066b241e943c8d0531067d057e774a04dcc341dadbd703ebe"
OLD_ARCHIVE_SHA256 = "bdfcc70f57d8cea981f508074690a7e940f884382c2452805f126a3f96bdb59d"
DEFERRED_RESULT_SHA256 = "3ba2e27cb66459ca75592de04bf6359fa3295733cf87aa1b02f4efcf62d30ec6"
LOOKUP_RESULT_SHA256 = "44da94d04c52e020799692049e8a6c58f1585a1ea4f794a0146ba933c71cfc4d"
RADIUS_NODE_IDS = (33982, 33981, 33830)
GAP_QUERY_FRONTIER = 33829
PROTECTED_NODE_IDS = frozenset(
    {0, 1, 17, 33572, 25220, 263344, 263348, 263349, 262918, 205147, 263346}
)
REQUIRED_INPUTS = {
    F_ARCHIVE.as_posix(): F_ARCHIVE_SHA256,
    OLD_ARCHIVE.as_posix(): OLD_ARCHIVE_SHA256,
    DEFERRED_RESULT.as_posix(): DEFERRED_RESULT_SHA256,
    LOOKUP_RESULT.as_posix(): LOOKUP_RESULT_SHA256,
}


class FGapRadiusSpineInputError(ValueError):
    """A frozen spine input is missing, malformed, or changed."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise FGapRadiusSpineInputError(message)


def _sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def _canon(value: Any) -> bytes:
    return json.dumps(
        value, ensure_ascii=False, sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode("utf-8")


def _canonical_sha(value: Any) -> str:
    return _sha(_canon(value))


def _path(value: Any) -> str:
    _require(isinstance(value, str) and value, "manifest path missing")
    return value.replace("\\", "/")


def _safe(root: Path, name: str) -> Path:
    relative = Path(name)
    _require(not relative.is_absolute() and ".." not in relative.parts, "manifest path escapes root")
    result = (root / relative).resolve()
    _require(result.is_relative_to(root.resolve()), "manifest path escapes root")
    return result


def _digest(value: Any, label: str) -> str:
    _require(isinstance(value, str) and len(value) == 64, label + " SHA malformed")
    value = value.lower()
    _require(all(char in "0123456789abcdef" for char in value), label + " SHA malformed")
    return value


def _manifest_hashes(source: Mapping[str, Any], hashes: Mapping[str, str] | None) -> dict[str, str]:
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
                _require(key not in entries or entries[key] == digest, "manifest hash conflict: " + key)
                entries[key] = digest
    if isinstance(hashes, Mapping):
        for name, value in hashes.items():
            key = _path(str(name))
            digest = _digest(value, "caller hashes." + key)
            _require(key not in entries or entries[key] == digest, "caller hash conflict: " + key)
            entries[key] = digest
    _require(entries, "caller frozen manifest has no hashes")
    return entries


def _find_hash(entries: Mapping[str, str], name: str) -> str:
    if name in entries:
        return entries[name]
    normalized = name.replace("\\", "/")
    candidates = [
        value for key, value in entries.items()
        if key.replace("\\", "/").endswith("/" + normalized)
    ]
    _require(len(candidates) == 1, "required manifest hash missing: " + name)
    return candidates[0]


def _authenticate(
    root: Path, source: Mapping[str, Any], hashes: Mapping[str, str] | None
) -> tuple[dict[str, bytes], dict[str, str]]:
    entries = _manifest_hashes(source, hashes)
    for name, expected in REQUIRED_INPUTS.items():
        _require(_find_hash(entries, name) == expected, "required SHA changed: " + name)
    raw_by_name: dict[str, bytes] = {}
    actual: dict[str, str] = {}
    # Hash every manifest entry before parsing JSON or decompressing gzip.
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
    for name, expected in REQUIRED_INPUTS.items():
        if name not in raw_by_name:
            data = _safe(root, name).read_bytes()
            observed = _sha(data)
            _require(observed == expected, "required raw SHA changed: " + name)
            raw_by_name[name] = data
            actual[name] = observed
    return raw_by_name, actual


def _json(data: bytes, label: str) -> Mapping[str, Any]:
    try:
        value = json.loads(data.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise FGapRadiusSpineInputError("invalid JSON: " + label) from exc
    _require(isinstance(value, Mapping), label + " is not an object")
    return value


def _archive(data: bytes, label: str) -> Mapping[str, Any]:
    try:
        decoded = gzip.decompress(data)
    except (OSError, EOFError) as exc:
        raise FGapRadiusSpineInputError("invalid gzip: " + label) from exc
    return _json(decoded, label)


def _prior(packet: Mapping[str, Any], schema: str, label: str) -> Mapping[str, Any]:
    _require(packet.get("status") == "returned", label + " wrapper status changed")
    result = packet.get("result")
    _require(isinstance(result, Mapping), label + " result missing")
    _require(result.get("schema") == schema, label + " schema changed")
    return result


def _opaque(value: Any, role: str) -> dict[str, Any]:
    return {
        "status": "opaque_source_ref" if type(value) is int else "unresolved_source_ref",
        "source_role": role,
        "source_node_id": value if type(value) is int else None,
        "node_id_is_not_a_numeric_value": True,
        "source_value_selected": False,
    }


def _frontier(value: Any, prefix: str = "") -> list[dict[str, Any]]:
    found: list[dict[str, Any]] = []
    if isinstance(value, Mapping):
        for key, item in value.items():
            found.extend(_frontier(item, f"{prefix}.{key}" if prefix else str(key)))
    elif isinstance(value, list):
        for index, item in enumerate(value):
            found.extend(_frontier(item, f"{prefix}[{index}]"))
    elif type(value) is int:
        found.append({
            "path": prefix,
            "node_id": value,
            "node_id_is_not_a_numeric_value": True,
            "child_row_indexed": False,
            "protected_root_candidate": value in PROTECTED_NODE_IDS,
            "protected_root_not_indexed": True,
            "owner": "unresolved_opaque_frontier",
        })
    return found


def _row(nodes: Any, node_id: int, side: str) -> dict[str, Any]:
    _require(isinstance(nodes, list), side + " exact_graph_nodes missing")
    _require(0 <= node_id < len(nodes), side + " fixed spine row outside archive")
    value = nodes[node_id]
    _require(isinstance(value, Mapping), side + " fixed spine row malformed")
    row = copy.deepcopy(dict(value))
    return {
        "side": side,
        "node_id": node_id,
        "node_id_is_not_a_numeric_value": True,
        "operation": row.get("operation") if isinstance(row.get("operation"), str) else None,
        "attributes": {str(key): copy.deepcopy(item) for key, item in row.items() if key != "operation"},
        "attributes_complete": True,
        "row_sha256": _canonical_sha(row),
        "opaque_frontier": _frontier({key: item for key, item in row.items() if key != "operation"}),
        "children_not_indexed": True,
        "protected_roots_not_indexed": True,
        "source_values_selected": False,
    }


def _tokens(row: Mapping[str, Any]) -> set[str]:
    values: set[str] = set()
    operation = row.get("operation")
    if isinstance(operation, str):
        values.add(operation.lower())
    attrs = row.get("attributes")
    if isinstance(attrs, Mapping):
        for key in ("name", "operator", "kind", "definition", "formula", "formula_text"):
            if isinstance(attrs.get(key), str):
                values.add(attrs[key].lower())
    return values


def _contains_id(value: Any, target: int) -> bool:
    if type(value) is int:
        return value == target
    if isinstance(value, Mapping):
        return any(_contains_id(item, target) for item in value.values())
    if isinstance(value, list):
        return any(_contains_id(item, target) for item in value)
    return False


def _formula(rows: Mapping[int, Mapping[str, Any]], side: str) -> dict[str, Any]:
    root, log_radius, scale = (rows.get(item) for item in RADIUS_NODE_IDS)
    if not all(isinstance(item, Mapping) for item in (root, log_radius, scale)):
        return {"status": "unresolved_missing_row", "side": side, "full_expression_proved": False}
    root_attrs = root.get("attributes")
    log_attrs = log_radius.get("attributes")
    scale_attrs = scale.get("attributes")
    if not isinstance(root_attrs, Mapping):
        root_attrs = {}
    if not isinstance(log_attrs, Mapping):
        log_attrs = {}
    if not isinstance(scale_attrs, Mapping):
        scale_attrs = {}

    # Match the FunctionTransportGraph schema exactly. Text containing exp or
    # an arbitrary occurrence of an ID is not operation evidence.
    root_argument = root_attrs.get("argument")
    root_shape = (
        root.get("operation") == "analytic_unary"
        and root_attrs.get("name") == "exp"
        and type(root_argument) is int
        and root_argument == 33981
    )

    log_arguments = log_attrs.get("arguments")
    log_arguments_are_two = isinstance(log_arguments, list) and len(log_arguments) == 2
    log_radius_refs = (
        [item for item in log_arguments if type(item) is int and item == 33830]
        if log_arguments_are_two
        else []
    )
    log_shape = (
        log_radius.get("operation") == "sum"
        and log_arguments_are_two
        and len(log_radius_refs) == 1
    )
    log_other = (
        next(
            (
                item for item in log_arguments
                if not (type(item) is int and item == 33830)
            ),
            None,
        )
        if log_arguments_are_two
        else None
    )

    scale_denominator_present = "denominator" in scale_attrs
    scale_numerator = scale_attrs.get("numerator")
    scale_shape = (
        scale.get("operation") == "positive_quotient"
        and type(scale_numerator) is int
        and scale_numerator == GAP_QUERY_FRONTIER
        and scale_denominator_present
    )

    reasons = []
    if not root_shape:
        reasons.append("radius root is not analytic_unary exp(argument=33981)")
    if not log_shape:
        reasons.append("log-radius row is not sum(arguments=[33830, opaque])")
    if not scale_shape:
        reasons.append("33830 is not positive_quotient(numerator=33829, denominator=opaque)")
    reasons.append("logRp sibling, quotient denominator, children, and binder/source-literal semantics are opaque")
    return {
        "status": "conditional_direct_shape" if root_shape and log_shape and scale_shape else "unresolved_direct_shape",
        "side": side,
        "expected": "R=exp(logRp+xi/mu)",
        "direct_root_exp_evidence": root_shape,
        "direct_root_log_radius_evidence": root_shape,
        "direct_log_radius_add_evidence": log_shape,
        "direct_scale_arithmetic_evidence": scale_shape,
        "root_operation": root.get("operation"),
        "root_name": root_attrs.get("name"),
        "root_argument": _opaque(root_argument, side + ".radius.argument"),
        "log_radius_operation": log_radius.get("operation"),
        "log_radius_arguments_are_exactly_two": log_arguments_are_two,
        "log_radius_contains_scale_row_33830": len(log_radius_refs) == 1,
        "log_radius_other_child": _opaque(log_other, side + ".log_radius.other_child"),
        "scale_operation": scale.get("operation"),
        "scale_numerator_is_old_gap_query_33829": (
            type(scale_numerator) is int and scale_numerator == GAP_QUERY_FRONTIER
        ),
        "scale_denominator": _opaque(
            scale_attrs.get("denominator"),
            side + ".positive_quotient.denominator",
        ) if scale_denominator_present else None,
        "query_frontier": _opaque(GAP_QUERY_FRONTIER, side + ".gap.query"),
        "full_expression_proved": False,
        "reasons": reasons,
    }


def _direct_three_operations_nonbinding(rows: Mapping[int, Mapping[str, Any]]) -> bool:
    """Recognize three operation labels without claiming binder semantics."""

    root, log_radius, scale = (rows.get(item) for item in RADIUS_NODE_IDS)
    if not all(isinstance(item, Mapping) for item in (root, log_radius, scale)):
        return False
    return (
        root.get("operation") == "analytic_unary"
        and log_radius.get("operation") == "sum"
        and scale.get("operation") == "positive_quotient"
    )


def _rewrites(bindings: Mapping[str, Any]) -> dict[str, Any]:
    mapping = bindings.get("rewritten_nodes")
    _require(isinstance(mapping, Mapping), "F rewritten_nodes missing")
    records = []
    for node_id in RADIUS_NODE_IDS:
        key = str(node_id)
        if key not in mapping:
            records.append({
                "source_ref": _opaque(node_id, "F.rewritten_nodes." + key),
                "status": "unresolved_missing_key",
                "rewritten_ref": None,
                "queried_key_only": True,
            })
            continue
        target = mapping[key]
        records.append({
            "source_ref": _opaque(node_id, "F.rewritten_nodes." + key),
            "status": "direct_key_present" if type(target) is int else "unresolved_target",
            "rewritten_ref": _opaque(target, "F.rewritten_nodes." + key),
            "same_opaque_id": type(target) is int and target == node_id,
            "target_row_read": False,
            "queried_key_only": True,
        })
    return {
        "map_owner": "Rp_native_C5_pulse_source_bindings.rewritten_nodes",
        "queried_keys": [str(item) for item in RADIUS_NODE_IDS],
        "records": records,
        "other_keys_not_read": True,
        "same_opaque_id_is_not_function_identity": True,
    }


def _membership(bindings: Mapping[str, Any], deferred: Mapping[str, Any]) -> dict[str, Any]:
    current = bindings.get("input_substitution")
    prior = deferred.get("archive_input_substitution")
    if not isinstance(current, Mapping) or not isinstance(prior, Mapping):
        return {
            "status": "unresolved_missing_substitution",
            "target_values_read": False,
            "lower": {"status": "unresolved", "target_values_read": False},
            "upper": {"status": "unresolved", "target_values_read": False},
        }
    current_keys = sorted(str(item) for item in current)
    prior_keys = sorted(str(item) for item in prior)
    return {
        "status": "membership_only",
        "F_keys": current_keys,
        "deferred_keys": prior_keys,
        "key_membership_equal": current_keys == prior_keys,
        "target_values_read": False,
        "lower": {"status": "unresolved_lower_owner", "membership_only": True},
        "upper": {"status": "unresolved_upper_owner", "membership_only": True},
    }


def _radius_binding_membership(bindings: Mapping[str, Any], radius_ref: Any) -> dict[str, Any]:
    substitutions = bindings.get("input_substitution")
    if type(radius_ref) is not int or not isinstance(substitutions, Mapping):
        return {
            "status": "unresolved_missing_input_substitution",
            "radius_key_membership": False,
            "target_value_read": False,
            "seeded_input_entry_is_not_recursive_proof": True,
        }
    key = str(radius_ref)
    present = key in substitutions
    target = substitutions.get(key)
    return {
        "status": "key_present" if present else "key_absent",
        "radius_key": _opaque(radius_ref, "F.input_substitution.radius"),
        "radius_key_membership": present,
        "target_ref": _opaque(target, "F.input_substitution.radius.target") if present else None,
        "target_value_read": False,
        "seeded_input_entry_is_not_recursive_proof": True,
        "function_identity_proved": False,
    }


def _bindings(archive: Mapping[str, Any], side: str) -> Mapping[str, Any]:
    key = "Rp_native_C5_pulse_source_bindings" if side == "F" else "C5_pulse_C3_pressure_bindings"
    value = archive.get(key)
    _require(isinstance(value, Mapping), side + " source bindings missing")
    return value


def _family(archive: Mapping[str, Any]) -> dict[str, Any]:
    value = archive.get("source_family")
    return {
        "present": value is not None,
        "metadata_sha256": _canonical_sha(value) if value is not None else None,
        "source_family_is_not_function_identity": True,
    }


def _prefix(old: Mapping[str, Any], frame: Mapping[str, Any]) -> dict[str, Any]:
    names = ("graph_prefix", "graph_namespace", "graph_owner", "graph_id")
    old_values = {name: old.get(name) for name in names if old.get(name) is not None}
    frame_values = {name: frame.get(name) for name in names if frame.get(name) is not None}
    common = sorted(set(old_values) & set(frame_values))
    equal = bool(common) and all(old_values[name] == frame_values[name] for name in common)
    return {
        "old": old_values,
        "F": frame_values,
        "common_keys": common,
        "same_graph_prefix_gate": equal,
        "status": "explicit_equal_metadata" if equal else "unresolved_no_explicit_common_prefix",
        "metadata_is_not_function_identity": True,
    }


def build_radius_spine_binding(
    root: str | Path = ".",
    source: Mapping[str, Any] | None = None,
    hashes: Mapping[str, str] | None = None,
) -> dict[str, Any]:
    """Build JSON-safe bounded source evidence; semantic gaps stay false."""

    _require(isinstance(source, Mapping), "caller frozen manifest is required")
    repo = Path(root).resolve()
    raw, actual = _authenticate(repo, source, hashes)
    deferred_packet = _json(raw[DEFERRED_RESULT.as_posix()], DEFERRED_RESULT.as_posix())
    lookup_packet = _json(raw[LOOKUP_RESULT.as_posix()], LOOKUP_RESULT.as_posix())
    deferred = _prior(deferred_packet, "historical_frame_gap_deferred_roles_v1", "deferred")
    lookup = _prior(lookup_packet, "historical_frame_gap_saved_rewrite_lookup_v1", "lookup")
    template = deferred.get("original_own_gap_field_template")
    roles = lookup.get("roles")
    _require(isinstance(template, Mapping), "deferred template missing")
    _require(isinstance(roles, Mapping) and isinstance(roles.get("field_radius"), Mapping), "lookup radius role missing")
    radius_role = roles["field_radius"]

    # This is the only archive decoding point, after the complete hash pass.
    old_archive = _archive(raw[OLD_ARCHIVE.as_posix()], OLD_ARCHIVE.as_posix())
    f_archive = _archive(raw[F_ARCHIVE.as_posix()], F_ARCHIVE.as_posix())
    old_bindings = _bindings(old_archive, "old")
    f_bindings = _bindings(f_archive, "F")
    old_nodes, f_nodes = old_archive.get("exact_graph_nodes"), f_archive.get("exact_graph_nodes")
    _require(isinstance(old_nodes, list) and isinstance(f_nodes, list), "exact_graph_nodes missing")
    old_rows = {item: _row(old_nodes, item, "old") for item in RADIUS_NODE_IDS}
    f_rows = {item: _row(f_nodes, item, "F") for item in RADIUS_NODE_IDS}
    old_ref, lookup_ref = template.get("radius"), radius_role.get("original_ref")
    radius_boundary_matches = type(old_ref) is int and type(lookup_ref) is int and old_ref == lookup_ref
    radius_same_handle = radius_boundary_matches and type(radius_role.get("rewritten_ref")) is int and radius_role.get("rewritten_ref") == lookup_ref
    # Check membership before interpreting any seeded rewritten-node entry;
    # membership is not a recursive-call or function-identity proof.
    radius_binding = _radius_binding_membership(f_bindings, old_ref)
    rewrite = _rewrites(f_bindings)
    row_records = []
    for index, item in enumerate(RADIUS_NODE_IDS):
        old_row, f_row = old_rows[item], f_rows[item]
        row_records.append({
            "node_id": item,
            "node_id_is_not_a_numeric_value": True,
            "old_row_sha256": old_row["row_sha256"],
            "F_row_sha256": f_row["row_sha256"],
            "direct_row_sha_equal": old_row["row_sha256"] == f_row["row_sha256"],
            "rewrite_record": rewrite["records"][index],
            "rewrite_target_row_not_indexed": True,
            "function_identity_proved": False,
        })
    direct_equal = all(item["direct_row_sha_equal"] for item in row_records)
    old_family, f_family = _family(old_archive), _family(f_archive)
    family_equal = bool(
        old_family["present"] and f_family["present"]
        and old_family["metadata_sha256"] == f_family["metadata_sha256"]
    )
    old_direct_three_operations_nonbinding = _direct_three_operations_nonbinding(old_rows)
    f_direct_three_operations_nonbinding = _direct_three_operations_nonbinding(f_rows)
    direct_three_operations_nonbinding = (
        old_direct_three_operations_nonbinding and f_direct_three_operations_nonbinding
    )
    return {
        "schema": SCHEMA,
        "component": COMPONENT,
        "status": "bounded_source_spine_unresolved",
        "source_only": True,
        "source_values_selected": False,
        "source_values_computed": False,
        "source_graph_traversed": False,
        "integrals_evaluated": False,
        "numeric_use_admissible": False,
        "certified": False,
        "duplicate_graph_computation": False,
        "old_leaf_imported": False,
        "distinct_archive_indices": 3,
        "archive_row_lookups": 6,
        "gzip_archives_parsed": 2,
        "direct_three_operations_nonbinding": direct_three_operations_nonbinding,
        "row_attributes_binder_semantics_proved": False,
        "frontier_binder_free_proved": False,
        "scope_flags": {
            "source_values_computed": False,
            "integrals_evaluated": False,
            "numeric_use_admissible": False,
            "certified": False,
            "source_graph_traversed": False,
            "duplicate_graph_computation": False,
            "old_leaf_imported": False,
            "old_producer_called": False,
            "ancestor_math_called": False,
            "function_identity_proved": False,
            "protected_root_indexed": False,
            "full_graph_traversal": False,
        },
        "raw_hash_authentication": {
            "all_manifest_raw_hashes_checked_before_decompression": True,
            "archives_decompressed_after_all_raw_hashes": True,
            "required_inputs": {
                name: {"path": name, "expected_sha256": digest, "actual_sha256": actual[name]}
                for name, digest in REQUIRED_INPUTS.items()
            },
        },
        "archive_scope": {
            "archive_json_decoded_not_graph_traversed": True,
            "fixed_rows_only": list(RADIUS_NODE_IDS),
            "protected_root_indexed": False,
            "protected_roots_not_indexed": True,
            "protected_root_ids_kept_opaque": sorted(PROTECTED_NODE_IDS),
            "gap_query_frontier": _opaque(GAP_QUERY_FRONTIER, "old/F.gap.query"),
        },
        "prior_scope": {
            "deferred_result_sha256": DEFERRED_RESULT_SHA256,
            "lookup_result_sha256": LOOKUP_RESULT_SHA256,
            "deferred_radius": {
                "radius": _opaque(old_ref, "deferred.original_own_gap_field_template.radius"),
                "query": _opaque(template.get("query"), "deferred.original_own_gap_field_template.query"),
                "log_radius": _opaque(template.get("log_radius"), "deferred.original_own_gap_field_template.log_radius"),
            },
            "lookup_radius": {
                "original_ref": _opaque(radius_role.get("original_ref"), "lookup.roles.field_radius.original_ref"),
                "rewritten_ref": _opaque(radius_role.get("rewritten_ref"), "lookup.roles.field_radius.rewritten_ref"),
                "status": radius_role.get("status"),
                "same_opaque_id": radius_same_handle,
                "same_opaque_id_is_not_function_identity": True,
            },
            "radius_boundary_record_matches_prior": radius_boundary_matches,
        },
        "source_bindings": {
            "old_archive": OLD_ARCHIVE.as_posix(),
            "F_archive": F_ARCHIVE.as_posix(),
            "archive_paths_distinct": True,
            "old_source_family": old_family,
            "F_source_family": f_family,
            "source_family_metadata_equal": family_equal,
            "same_graph_prefix": _prefix(old_archive, f_archive),
            "input_substitution_membership": _membership(f_bindings, deferred),
            "radius_input_substitution_membership": radius_binding,
            "source_family_or_prefix_is_not_function_identity": True,
        },
        "spine": {
            "old_rows": [old_rows[item] for item in RADIUS_NODE_IDS],
            "F_rows": [f_rows[item] for item in RADIUS_NODE_IDS],
            "rewritten_nodes": rewrite,
            "row_comparison": {
                "rows": row_records,
                "all_direct_row_sha_equal": direct_equal,
                "local_spine_rows_equal": {
                    "value": direct_equal,
                    "status": "conditional_direct_row_sha_fact",
                    "function_identity_proved": False,
                    "target_rows_and_children_not_read": True,
                },
            },
            "direct_row_sha_is_not_function_identity": True,
        },
        "formula_structure": {
            "expected": "R=exp(logRp+xi/mu)",
            "old": _formula(old_rows, "old"),
            "F": _formula(f_rows, "F"),
            "full_expression_proved": False,
            "all_children_opaque": True,
        },
        "source_literal_frontier_owners": {
            "old": [item for row in old_rows.values() for item in row["opaque_frontier"]],
            "F": [item for row in f_rows.values() for item in row["opaque_frontier"]],
            "all_owners_bound": False,
            "source_literals_are_not_values": True,
        },
        "conditional_facts": {
            "local_spine_rows_equal": direct_equal,
            "direct_three_operations_nonbinding": {
                "value": direct_three_operations_nonbinding,
                "old": old_direct_three_operations_nonbinding,
                "F": f_direct_three_operations_nonbinding,
                "status": "operation_labels_only",
                "rows_only": list(RADIUS_NODE_IDS),
                "row_attributes_binder_semantics_proved": False,
                "frontier_binder_free_proved": False,
                "capture_free_entire_expression": False,
                "subgraph_semantics_proved": False,
            },
            "local_spine_contains_only_binder_free_arithmetic": {
                "value": False,
                "status": "unresolved_opaque_frontiers",
                "direct_rows_only": True,
                "frontier_binder_free_proved": False,
            },
        },
        "source_gates": {
            "F_field_radius_owner_bound": False,
            "F_radius_input_substitution_membership_checked": True,
            "F_radius_input_substitution_recursive_proof": False,
            "capture_free_entire_expression": False,
            "ordinary_Z_transport_proved": False,
            "numeric_use_admissible": False,
            "old_errors_removed": False,
            "same_source_function_proved": False,
            "input_substitution_lower_closed": False,
            "input_substitution_upper_closed": False,
        },
        "retained_legacy_scope": {
            "numeric_use_admissible": False,
            "old_errors_removed": False,
            "integrals_evaluated": False,
            "source_values_selected": False,
        },
    }


def _synthetic_gate(
    left: Mapping[str, Any],
    right: Mapping[str, Any],
    left_frontier: Mapping[str, str],
    right_frontier: Mapping[str, str],
    capture_free: bool,
) -> bool:
    return _canonical_sha(left) == _canonical_sha(right) and dict(left_frontier) == dict(right_frontier) and capture_free


def verify_synthetic_fixture() -> dict[str, Any]:
    """Independent negative controls; this function reads no repository input."""

    row = {"operation": "exp", "argument": "n33981"}
    different_child = {"operation": "exp", "argument": "n33981_F"}
    old_frontier = {"logRp": "old.logRp", "xi": "old.xi", "mu": "old.mu"}
    f_frontier = {"logRp": "F.logRp", "xi": "F.xi", "mu": "F.mu"}
    checks = {
        "same_root_id_capture_rejected": not _synthetic_gate(row, row, old_frontier, f_frontier, False),
        "same_root_id_frontier_mismatch_rejected": not _synthetic_gate(row, row, old_frontier, f_frontier, True),
        "different_child_rejected": not _synthetic_gate(row, different_child, old_frontier, old_frontier, True),
        "equal_row_frontier_is_only_conditional": _synthetic_gate(row, row, old_frontier, old_frontier, True),
    }
    return {
        "schema": "f_gap_radius_spine_binding_synthetic_fixture_v1",
        "component": COMPONENT,
        "status": "synthetic_negative_controls_passed" if all(checks.values()) else "synthetic_negative_controls_failed",
        "synthetic_only": True,
        "actual_source_read": False,
        "source_graph_traversed": False,
        "source_values_computed": False,
        "source_values_selected": False,
        "integrals_evaluated": False,
        "numeric_use_admissible": False,
        "certified": False,
        "duplicate_graph_computation": False,
        "old_leaf_imported": False,
        "scope_flags": {
            "source_values_computed": False,
            "integrals_evaluated": False,
            "numeric_use_admissible": False,
            "certified": False,
            "source_graph_traversed": False,
            "duplicate_graph_computation": False,
            "old_leaf_imported": False,
            "old_producer_called": False,
            "ancestor_math_called": False,
            "function_identity_proved": False,
            "protected_root_indexed": False,
            "full_graph_traversal": False,
        },
        "checks": checks,
        "checks_passed": bool(all(checks.values())),
        "same_handle_is_not_function_identity": True,
        "finite_rows_are_not_function_values": True,
    }


__all__ = [
    "COMPONENT",
    "DEFERRED_RESULT",
    "DEFERRED_RESULT_SHA256",
    "F_ARCHIVE",
    "F_ARCHIVE_SHA256",
    "FGapRadiusSpineInputError",
    "LOOKUP_RESULT",
    "LOOKUP_RESULT_SHA256",
    "OLD_ARCHIVE",
    "OLD_ARCHIVE_SHA256",
    "RADIUS_NODE_IDS",
    "SCHEMA",
    "build_radius_spine_binding",
    "verify_synthetic_fixture",
]
