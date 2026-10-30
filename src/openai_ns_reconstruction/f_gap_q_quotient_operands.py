"""Bounded, source-only observation of the three F-Q quotient operands.

This leaf is the post-freeze consumer for the returned
``f_gap_q_source_rows_v1`` packet.  It authenticates every byte named by the
caller manifest before decoding either JSON packet, checks the two root
quotient references in that returned packet, and then indexes exactly
``277041``, ``277054``, and ``25272`` in the F archive.  Operand references
are opaque and no referenced child is followed.

The operation schema is copied as data from the importer contract.  It is
deliberately not imported: an operation not present in ``REFERENCE_KEYS`` is
unknown and none of its integer fields is guessed to be a graph edge.
"""

from __future__ import annotations

import copy
import gzip
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping


SCHEMA = "f_gap_q_quotient_operands_v1"
F_GRAPH_NAMESPACE = "F_exact_graph"

F_ARCHIVE = Path(
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_Rp_native_C5_pulse_source.json.gz"
)
F_ARCHIVE_SHA256 = "442df70e5ab4051066b241e943c8d0531067d057e774a04dcc341dadbd703ebe"

PRIOR_COMPONENT = "f_gap_q_source_rows_v1"
PRIOR_WRAPPER_SCHEMA = PRIOR_COMPONENT + "_result_v1"
PRIOR_RESULT = Path("tmp/f_gap_q_source_rows_v1/result.json")
PRIOR_RESULT_SHA256 = "27c9e86eabd988f5494086c142dfe03669561558389a880b4a6891f3af350f36"
PRIOR_INNER_STATUSES = frozenset(
    {
        "bounded_F_Q_source_rows_observed",
        "bounded_F_Q_source_rows_observed_prior_unresolved",
    }
)

# These are copied from the static Q-source-row definition.  The quotient
# root rows are prior observations; this leaf never indexes them again.
Q_ROOTS: dict[str, dict[str, int]] = {
    "Q0": {"root": 277042, "numerator": 277041, "denominator": 25272},
    "Q1": {"root": 277055, "numerator": 277054, "denominator": 25272},
}
OPERAND_NODE_IDS = (277041, 277054, 25272)

# Frozen operand schema copied as data from the source importer.  Only fields
# registered here can become an edge descriptor.  All other row fields stay
# in the verbatim row packet and are not interpreted.
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

MAX_TYPED_EDGE_REFS = 32


class FGapQQuotientOperandsInputError(ValueError):
    """A frozen input or bounded operand contract is malformed."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise FGapQQuotientOperandsInputError(message)


def _sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def _canon(value: Any) -> bytes:
    return json.dumps(
        value,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    ).encode("utf-8")


def _digest(value: Any, label: str) -> str:
    _require(isinstance(value, str) and len(value) == 64, label + " SHA malformed")
    value = value.lower()
    _require(all(ch in "0123456789abcdef" for ch in value), label + " SHA malformed")
    return value


def _manifest_hashes(
    source: Mapping[str, Any], hashes: Mapping[str, Any] | None
) -> dict[str, str]:
    """Collect and validate every caller-supplied manifest hash."""

    entries: dict[str, str] = {}

    def add(name: Any, digest: Any, label: str) -> None:
        _require(isinstance(name, str) and name, "manifest path missing")
        key = name.replace("\\", "/")
        value = _digest(digest, label)
        _require(
            key not in entries or entries[key] == value,
            "manifest hash conflict: " + key,
        )
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
                value
                for name, value in entries.items()
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
        required
        for required in REQUIRED_INPUTS
        if normalized == required or normalized.endswith("/" + required)
    ]
    _require(len(aliases) == 1, "absolute manifest alias is outside required inputs")
    return aliases[0]


def _safe(root: Path, name: str) -> Path:
    relative = Path(name)
    _require(
        not relative.is_absolute() and ".." not in relative.parts,
        "manifest path escapes root",
    )
    result = (root / relative).resolve()
    _require(result.is_relative_to(root.resolve()), "manifest path escapes root")
    return result


def _authenticate(
    root: Path, source: Mapping[str, Any], hashes: Mapping[str, Any] | None
) -> tuple[dict[str, bytes], dict[str, str], int]:
    """Hash every supplied byte before parsing JSON or decompressing gzip."""

    entries = _manifest_hashes(source, hashes)
    raw: dict[str, bytes] = {}
    actual: dict[str, str] = {}
    for raw_name, expected in entries.items():
        name = _manifest_name(raw_name)
        data = _safe(root, name).read_bytes()
        observed = _sha(data)
        _require(observed == expected, "raw SHA changed: " + name)
        _require(
            name not in actual or actual[name] == observed,
            "manifest alias conflict: " + name,
        )
        raw[name], actual[name] = data, observed

    # A caller may provide a manifest that records only aliases.  The two
    # required bytes are still read and checked before either is decoded.
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
        raise FGapQQuotientOperandsInputError("invalid JSON: " + label) from exc
    _require(isinstance(value, Mapping), label + " is not an object")
    return value


def _gzip_json(raw: bytes, label: str) -> Mapping[str, Any]:
    try:
        decoded = gzip.decompress(raw)
    except (OSError, EOFError) as exc:
        raise FGapQQuotientOperandsInputError("invalid gzip: " + label) from exc
    return _json(decoded, label)


def _node_from_ref(value: Any) -> int | None:
    """Inspect one opaque reference object only; never follow its node ID."""

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


def _opaque_ref(node: Any, source_role: str) -> dict[str, Any]:
    valid = type(node) is int and node >= 0
    return {
        "status": "opaque_source_ref" if valid else "unresolved_source_ref",
        "namespace": F_GRAPH_NAMESPACE,
        "node": node if valid else None,
        "source_node_id": node if valid else None,
        "source_role": source_role,
        "node_id_is_not_a_numeric_value": True,
        "source_value_selected": False,
        "indexed": False,
        "protected_node_indexed": False,
    }


def _edge(field: str, node: int, source_role: str) -> dict[str, Any]:
    return {
        "field": field,
        "typed_edge": True,
        "ref": _opaque_ref(node, source_role + "." + field),
        "child_indexed": False,
        "edge_followed": False,
        "protected_target": False,
    }


def _typed_edge_refs(row: Mapping[str, Any], node: int) -> dict[str, Any]:
    """Extract only fields registered for the row's operation."""

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
                edges.append(
                    _edge(
                        edge_field,
                        child,
                        f"F_archive.exact_graph_nodes[{node}]",
                    )
                )

    registered = set(fields or ())
    ignored_integer_attributes = sorted(
        str(key)
        for key, value in row.items()
        if key != "operation"
        and key not in registered
        and (
            type(value) is int
            or isinstance(value, list)
            and any(type(item) is int for item in value)
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
        "function_identity_proved": False,
        "ordinary_Z_derivative_proved": False,
        "derivative_proved": False,
        "Q1_is_partial_Z_Q0_proved": False,
        "positivity_proved": False,
        "denominator_positive_proved": False,
        "source_values_selected": False,
        "source_values_computed": False,
        "numeric_evaluation": False,
        "integrals_evaluated": False,
        "numeric_use_admissible": False,
        "error_transport_proved": False,
        "error_removal": False,
        "closure_proved": False,
    }


def _observed_row(row: Mapping[str, Any], node: int, record_key: str) -> dict[str, Any]:
    copied = copy.deepcopy(dict(row))
    edge_packet = _typed_edge_refs(copied, node)
    return {
        "status": "direct_observed_row",
        "record_key": record_key,
        "node": node,
        "source_node_id": node,
        "namespace": F_GRAPH_NAMESPACE,
        "operation": copied.get("operation"),
        "field_count": len(copied),
        "argument_count": (
            len(copied.get("arguments"))
            if isinstance(copied.get("arguments"), list)
            else None
        ),
        "typed_edge_ref_count": edge_packet["count"],
        "typed_edge_refs": edge_packet["refs"],
        "typed_edge_scope": edge_packet,
        "row": copied,
        "raw_fields": copy.deepcopy(copied),
        "row_sha256": _sha(_canon(copied)),
        "canonical_row_sha256": _sha(_canon(copied)),
        "source_value_selected": False,
        "children_indexed": False,
        "child_rows_read": 0,
        "graph_edges_followed": False,
        "semantic_identity_proved": False,
    }


def _prior_q_source(packet: Mapping[str, Any]) -> Mapping[str, Any]:
    """Accept only the returned, fixed prior source-row wrapper."""

    _require(packet.get("schema") == PRIOR_WRAPPER_SCHEMA, "prior wrapper schema changed")
    _require(packet.get("component") == PRIOR_COMPONENT, "prior wrapper component changed")
    _require(packet.get("status") == "returned", "prior wrapper status is not returned")
    _require("error" in packet and packet.get("error") is None, "prior wrapper contains an error")
    result = packet.get("result")
    _require(isinstance(result, Mapping), "prior wrapper result missing")
    candidate = result.get("q_source_rows")
    _require(isinstance(candidate, Mapping), "prior wrapper q_source_rows missing")
    _require(candidate.get("schema") == PRIOR_COMPONENT, "prior inner schema changed")
    _require(candidate.get("component") == PRIOR_COMPONENT, "prior inner component changed")
    _require(
        candidate.get("status") in PRIOR_INNER_STATUSES,
        "prior inner status is not a recognized source-row result",
    )
    _require(candidate.get("source_only") is True, "prior inner packet is not source-only")
    return candidate


def _require_opaque_ref(value: Any, node: int, label: str) -> dict[str, Any]:
    _require(isinstance(value, Mapping), label + " missing")
    _require(value.get("status") == "opaque_source_ref", label + " status changed")
    _require(value.get("node") == node, label + " node mismatch")
    _require(value.get("source_node_id") == node, label + " source node mismatch")
    _require(value.get("node_id_is_not_a_numeric_value") is True, label + " is numeric")
    _require(value.get("namespace") == F_GRAPH_NAMESPACE, label + " namespace mismatch")
    _require(value.get("source_value_selected") is False, label + " selected a source value")
    _require(value.get("indexed") is False, label + " was indexed")
    return copy.deepcopy(dict(value))


def _prior_root_binding(prior: Mapping[str, Any]) -> dict[str, Any]:
    """Check root quotient refs without using them to select archive rows."""

    rows = prior.get("rows")
    _require(isinstance(rows, Mapping), "prior source rows missing")
    direct_scope = prior.get("direct_graph_scope")
    _require(isinstance(direct_scope, Mapping), "prior direct graph scope missing")
    prior_indexed = direct_scope.get("indexed_rows")
    _require(
        isinstance(prior_indexed, list)
        and prior_indexed == [Q_ROOTS["Q0"]["root"], Q_ROOTS["Q1"]["root"]],
        "prior root row scope changed",
    )
    for key in ("children_not_indexed", "graph_edges_followed", "source_graph_traversed"):
        if key in direct_scope:
            _require(direct_scope[key] is (True if key == "children_not_indexed" else False),
                     "prior direct scope marker changed: " + key)

    bound: dict[str, Any] = {}
    for key, expected in Q_ROOTS.items():
        packet = rows.get(key)
        _require(isinstance(packet, Mapping), "prior root packet missing: " + key)
        root_ref = _require_opaque_ref(packet.get("opaque_ref"), expected["root"], key + ".opaque_ref")
        direct = packet.get("direct_observed_row")
        _require(isinstance(direct, Mapping), "prior direct root row missing: " + key)
        _require(direct.get("node") == expected["root"], key + " root node changed")
        _require(direct.get("namespace") == F_GRAPH_NAMESPACE, key + " root namespace changed")
        _require(direct.get("children_indexed") is False, key + " child scope changed")
        _require(direct.get("child_rows_read") == 0, key + " child row count changed")

        row = direct.get("row")
        _require(isinstance(row, Mapping), key + " root raw row missing")
        _require(row.get("operation") == "positive_quotient", key + " root operation changed")
        _require(row.get("numerator") == expected["numerator"], key + " numerator changed")
        _require(row.get("denominator") == expected["denominator"], key + " denominator changed")

        scope = direct.get("typed_edge_scope")
        _require(isinstance(scope, Mapping), key + " typed edge scope missing")
        _require(scope.get("operation") == "positive_quotient", key + " typed operation changed")
        _require(scope.get("operation_known") is True, key + " operation is not known")
        _require(scope.get("unknown_operation") is False, key + " operation is unknown")
        _require(
            scope.get("registered_reference_fields") == ["numerator", "denominator"],
            key + " quotient schema changed",
        )
        _require(scope.get("children_indexed") is False, key + " typed child scope changed")
        _require(scope.get("edges_followed") is False, key + " edge-follow scope changed")

        typed = direct.get("typed_edge_refs")
        _require(isinstance(typed, list), key + " typed refs missing")
        _require(len(typed) == 2, key + " typed ref count changed")
        by_field: dict[str, Mapping[str, Any]] = {}
        for edge in typed:
            _require(isinstance(edge, Mapping), key + " typed edge malformed")
            field = edge.get("field")
            _require(field in ("numerator", "denominator"), key + " typed field changed")
            _require(field not in by_field, key + " duplicate typed field")
            _require(edge.get("typed_edge") is True, key + " edge is not typed")
            _require(edge.get("child_indexed") is False, key + " edge indexed a child")
            _require(edge.get("edge_followed") is False, key + " edge was followed")
            by_field[field] = edge
        _require(set(by_field) == {"numerator", "denominator"}, key + " typed fields missing")
        refs: dict[str, Any] = {}
        for field in ("numerator", "denominator"):
            expected_node = expected[field]
            edge = by_field[field]
            ref = _require_opaque_ref(edge.get("ref"), expected_node, key + "." + field + ".ref")
            _require(_node_from_ref(row.get(field)) == expected_node, key + " row/ref mismatch")
            refs[field] = ref
        bound[key] = {
            "root_ref": root_ref,
            "root_node": expected["root"],
            "operation": "positive_quotient",
            "typed_refs": refs,
            "row_sha256": direct.get("row_sha256"),
            "canonical_row_sha256": direct.get("canonical_row_sha256", direct.get("row_sha256")),
            "source_value_selected": False,
            "children_indexed": False,
            "child_rows_read": 0,
        }

    snapshot = prior.get("F_snapshot")
    _require(isinstance(snapshot, Mapping), "prior F snapshot missing")
    _require(snapshot.get("namespace") == F_GRAPH_NAMESPACE, "prior F snapshot namespace changed")
    for key, expected in Q_ROOTS.items():
        _require_opaque_ref(snapshot.get(key), expected["root"], "prior F snapshot " + key)
    _require(snapshot.get("whole_json_decoded") is True, "prior F snapshot decode marker missing")
    return {
        "status": "matched",
        "schema": prior.get("schema"),
        "component": prior.get("component"),
        "inner_status": prior.get("status"),
        "namespace": F_GRAPH_NAMESPACE,
        "root_quotients": bound,
        "root_quotient_typed_refs_bound": True,
        "source_value_selected": False,
        "semantic_identity_proved": False,
    }


def _prior_binding(packet: Mapping[str, Any]) -> dict[str, Any]:
    prior = _prior_q_source(packet)
    binding = _prior_root_binding(prior)
    return {
        "path": PRIOR_RESULT.as_posix(),
        "sha256": PRIOR_RESULT_SHA256,
        "wrapper_schema": packet.get("schema"),
        "wrapper_component": packet.get("component"),
        "wrapper_status": packet.get("status"),
        "wrapper_error": None,
        **binding,
    }


def _direct_row(nodes: Any, node: int, record_key: str) -> dict[str, Any]:
    _require(isinstance(nodes, list), "F exact_graph_nodes missing")
    _require(node < len(nodes), "requested F operand node outside archive")
    row = nodes[node]
    _require(isinstance(row, Mapping), "requested F operand row malformed")
    return _observed_row(row, node, record_key)


def _row_summary(records: Mapping[str, Any]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, packet in records.items():
        row = packet if isinstance(packet, Mapping) else {}
        result[key] = {
            "operation": row.get("operation"),
            "field_count": row.get("field_count"),
            "argument_count": row.get("argument_count"),
            "typed_edge_ref_count": row.get("typed_edge_ref_count"),
            "row_sha256": row.get("row_sha256"),
        }
    return result


def verify_synthetic_fixture() -> dict[str, Any]:
    """Exercise schema handling and failure gates using in-memory rows only."""

    quotient = _observed_row(
        {
            "operation": "positive_quotient",
            "numerator": 11,
            "denominator": 12,
            "literal": 13,
        },
        100,
        "synthetic.quotient",
    )
    unknown = _observed_row(
        {"operation": "future_unknown", "arguments": [201], "literal": 202},
        101,
        "synthetic.unknown",
    )

    bad_prior = {
        "schema": PRIOR_WRAPPER_SCHEMA,
        "component": PRIOR_COMPONENT,
        "status": "failed",
        "error": {"type": "SyntheticFailure"},
        "result": {},
    }
    bad_prior_rejected = False
    try:
        _prior_q_source(bad_prior)
    except FGapQQuotientOperandsInputError:
        bad_prior_rejected = True

    class GuardedNodes(list):
        def __init__(self, values: list[Any]) -> None:
            super().__init__(values)
            self.accessed: list[int] = []

        def __getitem__(self, index: Any) -> Any:
            if type(index) is int:
                self.accessed.append(index)
            return super().__getitem__(index)

    guarded = GuardedNodes(
        [
            {"operation": "exact_rational", "value": 0},
            {"operation": "sum", "arguments": [999]},
            {"operation": "product", "arguments": [998]},
            {"operation": "bound_variable", "name": "synthetic"},
        ]
    )
    direct = _direct_row(guarded, 3, "synthetic.direct")

    checks = {
        "quotient_schema_fields_are_typed": (
            quotient["typed_edge_scope"]["operation_known"] is True
            and quotient["typed_edge_scope"]["registered_reference_fields"]
            == ["numerator", "denominator"]
            and quotient["typed_edge_ref_count"] == 2
        ),
        "unknown_operation_is_not_guessed": (
            unknown["typed_edge_ref_count"] == 0
            and unknown["typed_edge_scope"]["unknown_operation"] is True
            and unknown["typed_edge_scope"]["operation_known"] is False
            and "arguments" in unknown["typed_edge_scope"]["ignored_integer_attributes"]
        ),
        "literal_fields_are_retained_but_not_edges": (
            quotient["raw_fields"]["literal"] == 13
            and quotient["typed_edge_scope"]["ignored_integer_attributes"] == ["literal"]
        ),
        "bad_prior_status_and_schema_fail_closed": bad_prior_rejected,
        "direct_read_does_not_read_children": (
            guarded.accessed == [3]
            and direct["children_indexed"] is False
            and direct["child_rows_read"] == 0
            and direct["typed_edge_scope"]["edges_followed"] is False
        ),
        "row_sha_is_canonical_and_structural": (
            len(direct["row_sha256"]) == 64
            and direct["row_sha256"] == direct["canonical_row_sha256"]
            and direct["semantic_identity_proved"] is False
        ),
        "semantic_flags_stay_false": all(
            not value for value in _semantic_false_flags().values()
        ),
    }
    return {
        "schema": SCHEMA + "_synthetic_fixture_v1",
        "component": SCHEMA,
        "status": "synthetic_checks_passed" if all(checks.values()) else "synthetic_checks_failed",
        "passed": all(checks.values()),
        "checks_passed": all(checks.values()),
        "synthetic_only": True,
        "checks": checks,
        "check_count": len(checks),
        "fixtures": {"quotient": quotient, "unknown": unknown, "direct": direct},
        "semantic_proof": _semantic_false_flags(),
        "real_inputs_read": False,
        "archive_decoded": False,
        "source_values_selected": False,
        "children_indexed": False,
        "child_rows_read": 0,
    }


def build_quotient_operands(
    root: str | Path = ".",
    source: Mapping[str, Any] | None = None,
    hashes: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Authenticate frozen inputs, then observe exactly three F graph rows."""

    _require(isinstance(source, Mapping), "caller frozen manifest is required")
    raw, actual, manifest_count = _authenticate(Path(root).resolve(), source, hashes)

    # No packet is decoded until every manifest byte, including both required
    # bytes, has been hashed above.
    prior_packet = _json(raw[PRIOR_RESULT.as_posix()], PRIOR_RESULT.as_posix())
    prior = _prior_binding(prior_packet)
    archive = _gzip_json(raw[F_ARCHIVE.as_posix()], F_ARCHIVE.as_posix())
    nodes = archive.get("exact_graph_nodes")
    _require(isinstance(nodes, list), "F exact_graph_nodes missing")

    records = {
        "Q0_numerator": _direct_row(nodes, Q_ROOTS["Q0"]["numerator"], "Q0.numerator"),
        "Q1_numerator": _direct_row(nodes, Q_ROOTS["Q1"]["numerator"], "Q1.numerator"),
        "shared_denominator": _direct_row(
            nodes, Q_ROOTS["Q0"]["denominator"], "shared.denominator"
        ),
    }
    indexed = list(OPERAND_NODE_IDS)
    summary = _row_summary(records)
    semantic = _semantic_false_flags()
    return {
        "schema": SCHEMA,
        "component": SCHEMA,
        "status": "bounded_F_Q_quotient_operands_observed"
        if prior["inner_status"] == "bounded_F_Q_source_rows_observed"
        else "bounded_F_Q_quotient_operands_observed_prior_unresolved",
        "source_only": True,
        "source_rows_only": True,
        "runner_interface": {
            "entry": "build_quotient_operands",
            "arguments": ["root", "source", "hashes"],
            "one_shot": True,
            "coordination": "/root/operands_runner",
        },
        "raw_hash_authentication": {
            "all_manifest_raw_hashes_checked_before_decode": True,
            "manifest_entry_count": manifest_count,
            "required_inputs": {
                name: {
                    "path": name,
                    "expected_sha256": digest,
                    "actual_sha256": actual[name],
                }
                for name, digest in REQUIRED_INPUTS.items()
            },
        },
        "prior_result_binding": prior,
        "F_snapshot": {
            "archive_path": F_ARCHIVE.as_posix(),
            "namespace": F_GRAPH_NAMESPACE,
            "operand_node_ids": indexed,
            "whole_json_decoded": True,
            "source_identity_proved": False,
        },
        "direct_graph_scope": {
            "requested_node_ids": indexed,
            "direct_rows": indexed,
            "indexed_rows": indexed,
            "indexed_row_count": 3,
            "maximum_indexed_rows": 3,
            "operand_rows_are_direct_observations": True,
            "operand_references_are_opaque": True,
            "children_not_indexed": True,
            "child_rows_read": 0,
            "graph_edges_followed": False,
            "source_graph_traversed": False,
            "ancestor_replay": False,
            "numeric_evaluation": False,
            "protected_nodes_indexed": False,
            "forbidden_nodes_indexed": False,
        },
        "rows": records,
        "row_operation_summary": summary,
        "indexed_row_count": 3,
        "semantic_proof": semantic,
        "source_gates": {
            "prior_wrapper_status_schema_verified": True,
            "root_quotient_typed_refs_bound": True,
            **semantic,
        },
        "scope_statement": (
            "Exactly F graph rows 277041, 277054, and 25272 were directly "
            "observed after raw hash authentication.  The prior returned "
            "root quotient packet was checked for its operation-registered "
            "numerator and denominator opaque refs.  Operand row fields and "
            "canonical row SHA-256 values are retained; no child, ancestor, "
            "source value, formula identity, derivative, positivity, numeric "
            "value, error, or closure semantic is established."
        ),
    }


build_f_gap_q_quotient_operands = build_quotient_operands
build_q_quotient_operands = build_quotient_operands


__all__ = [
    "FGapQQuotientOperandsInputError",
    "build_quotient_operands",
    "build_f_gap_q_quotient_operands",
    "build_q_quotient_operands",
    "verify_synthetic_fixture",
]
