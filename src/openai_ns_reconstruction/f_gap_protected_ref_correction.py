"""Correct only the stale protection labels in the returned Q-operands packet.

This consumer authenticates the returned JSON bytes before decoding them.  It
does not import the producer, open the F archive, follow an edge, or interpret
an opaque node identifier as a scalar.  The producer's reference-field table
is copied below as data so a row field becomes an edge only when the frozen
producer schema says that it is one.
"""

from __future__ import annotations

import copy
import hashlib
import json
from collections.abc import Mapping
from pathlib import Path
from typing import Any


SCHEMA = "f_gap_protected_ref_correction_v1"
PRIOR_COMPONENT = "f_gap_q_quotient_operands_v1"
PRIOR_WRAPPER_SCHEMA = PRIOR_COMPONENT + "_result_v1"
PRIOR_RESULT = Path("tmp/f_gap_q_quotient_operands_v1/result.json")
PRIOR_RESULT_SHA256 = "0ae2076688bd01072504f282bf82a2c3346555abdbb7ec4a8d2828a4ab2a80cc"
F_GRAPH_NAMESPACE = "F_exact_graph"

PROTECTED_NODE_IDS = frozenset(
    {0, 1, 17, 33572, 25220, 263344, 263348, 263349, 262918, 205147, 263346, 33982}
)

# Static snapshot of the producer's operation schema.  Literal integer fields
# not listed here are metadata, never graph-edge candidates.
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

EXPECTED_ROWS: dict[str, dict[str, Any]] = {
    "Q0_numerator": {"node": 277041, "arguments": [277038, 277040]},
    "Q1_numerator": {"node": 277054, "arguments": [277051, 277053]},
    "shared_denominator": {"node": 25272, "arguments": [1, 25271]},
}

SEMANTIC_FALSE_KEYS = (
    "source_identity_proved",
    "same_source_function_proved",
    "function_identity_proved",
    "ordinary_Z_derivative_proved",
    "derivative_proved",
    "Q1_is_partial_Z_Q0_proved",
    "positivity_proved",
    "denominator_positive_proved",
    "source_values_selected",
    "source_values_computed",
    "numeric_evaluation",
    "integrals_evaluated",
    "numeric_use_admissible",
    "error_transport_proved",
    "error_removal",
    "closure_proved",
)


class FGapProtectedRefCorrectionInputError(ValueError):
    """A frozen correction input or direct-edge contract is malformed."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise FGapProtectedRefCorrectionInputError(message)


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
    _require(all(char in "0123456789abcdef" for char in value), label + " SHA malformed")
    return value


def _manifest_entries(
    source: Mapping[str, Any], hashes: Mapping[str, Any] | None
) -> list[tuple[Any, Any, str]]:
    entries: list[tuple[Any, Any, str]] = []

    if "files" in source:
        files = source.get("files")
        _require(isinstance(files, list), "manifest files malformed")
        for index, item in enumerate(files):
            _require(isinstance(item, Mapping), "manifest files entry malformed")
            entries.append((item.get("path"), item.get("sha256"), f"manifest.files[{index}]"))

    for field in ("sha256", "input_hashes", "hashes", "input_sha256", "raw_input_sha256"):
        if field not in source:
            continue
        table = source.get(field)
        _require(isinstance(table, Mapping), "manifest." + field + " malformed")
        entries.extend((name, digest, "manifest." + field) for name, digest in table.items())

    if hashes is not None:
        _require(isinstance(hashes, Mapping), "caller hashes malformed")
        entries.extend((name, digest, "caller hashes") for name, digest in hashes.items())

    _require(entries, "caller frozen manifest has no hashes")
    return entries


def _manifest_name(name: Any) -> str:
    _require(isinstance(name, str) and name, "manifest path missing")
    normalized = name.replace("\\", "/")
    absolute = normalized.startswith("/") or (
        len(normalized) >= 3 and normalized[1] == ":" and normalized[2] == "/"
    )
    if not absolute:
        return normalized
    _require(
        normalized == PRIOR_RESULT.as_posix()
        or normalized.endswith("/" + PRIOR_RESULT.as_posix()),
        "absolute manifest alias is outside required input",
    )
    return PRIOR_RESULT.as_posix()


def _safe(root: Path, name: str) -> Path:
    relative = Path(name)
    _require(not relative.is_absolute() and ".." not in relative.parts, "manifest path escapes root")
    candidate = (root / relative).resolve()
    _require(candidate.is_relative_to(root.resolve()), "manifest path escapes root")
    return candidate


def _authenticate(
    root: Path, source: Mapping[str, Any], hashes: Mapping[str, Any] | None
) -> tuple[dict[str, bytes], dict[str, str], int]:
    """Hash every named byte before decoding the required JSON packet."""

    expected: dict[str, str] = {}
    for raw_name, digest, label in _manifest_entries(source, hashes):
        name = _manifest_name(raw_name)
        value = _digest(digest, label + "." + name)
        _require(name not in expected or expected[name] == value, "manifest hash conflict: " + name)
        expected[name] = value

    _require(
        expected.get(PRIOR_RESULT.as_posix()) == PRIOR_RESULT_SHA256,
        "required SHA changed: " + PRIOR_RESULT.as_posix(),
    )

    raw: dict[str, bytes] = {}
    actual: dict[str, str] = {}
    for name, digest in expected.items():
        data = _safe(root, name).read_bytes()
        observed = _sha(data)
        _require(observed == digest, "raw SHA changed: " + name)
        raw[name] = data
        actual[name] = observed

    # This also supports an absolute alias manifest while keeping the required
    # byte in the post-hash decode set.
    if PRIOR_RESULT.as_posix() not in raw:
        data = _safe(root, PRIOR_RESULT.as_posix()).read_bytes()
        observed = _sha(data)
        _require(observed == PRIOR_RESULT_SHA256, "required raw SHA changed")
        raw[PRIOR_RESULT.as_posix()] = data
        actual[PRIOR_RESULT.as_posix()] = observed
    return raw, actual, len(expected)


def _json(raw: bytes, label: str) -> Mapping[str, Any]:
    try:
        value = json.loads(raw.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise FGapProtectedRefCorrectionInputError("invalid JSON: " + label) from exc
    _require(isinstance(value, Mapping), label + " is not an object")
    return value


def _false_flags() -> dict[str, bool]:
    return {name: False for name in SEMANTIC_FALSE_KEYS}


def _node_from_ref(ref: Mapping[str, Any], label: str) -> int:
    _require(ref.get("status") == "opaque_source_ref", label + " status changed")
    _require(ref.get("namespace") == F_GRAPH_NAMESPACE, label + " namespace changed")
    node = ref.get("node")
    source_node = ref.get("source_node_id")
    _require(type(node) is int and node >= 0, label + " node malformed")
    _require(source_node == node, label + " source node mismatch")
    _require(ref.get("node_id_is_not_a_numeric_value") is True, label + " became numeric")
    _require(ref.get("source_value_selected") is False, label + " selected a source value")
    _require(ref.get("indexed") is False, label + " was indexed")
    _require(ref.get("protected_node_indexed") is False, label + " protected node was indexed")
    return node


def _edge_review(
    edge: Mapping[str, Any], expected_field: str, expected_node: int, label: str
) -> tuple[dict[str, Any], dict[str, Any]]:
    _require(edge.get("field") == expected_field, label + " field changed")
    _require(edge.get("typed_edge") is True, label + " is not typed")
    ref = edge.get("ref")
    _require(isinstance(ref, Mapping), label + " ref missing")
    node = _node_from_ref(ref, label + ".ref")
    _require(node == expected_node, label + " node changed")
    _require(edge.get("child_indexed") is False, label + " child indexed")
    _require(edge.get("edge_followed") is False, label + " edge followed")
    previous = edge.get("protected_target")
    _require(type(previous) is bool, label + " previous protected flag malformed")
    corrected = node in PROTECTED_NODE_IDS
    review = {
        "field": expected_field,
        "node": node,
        "source_node_id": node,
        "previous_protected_target": previous,
        "protected_target": corrected,
        "corrected_protected_target": corrected,
        "changed": previous != corrected,
        "previous_false_new_true": previous is False and corrected is True,
        "direct_dependency": True,
        "target_is_opaque": True,
        "indexed": False,
        "edge_followed": False,
    }
    fixed = copy.deepcopy(dict(edge))
    fixed["protected_target"] = corrected
    return review, fixed


def _row(
    key: str, packet: Mapping[str, Any], expected: Mapping[str, Any]
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    node = expected["node"]
    _require(packet.get("status") == "direct_observed_row", key + " row status changed")
    _require(packet.get("node") == node, key + " row node changed")
    _require(packet.get("source_node_id") == node, key + " source node changed")
    _require(packet.get("namespace") == F_GRAPH_NAMESPACE, key + " namespace changed")
    _require(packet.get("children_indexed") is False, key + " children indexed")
    _require(packet.get("child_rows_read") == 0, key + " child row read")
    _require(packet.get("graph_edges_followed") is False, key + " edge followed")
    _require(packet.get("source_value_selected") is False, key + " source value selected")

    raw_row = packet.get("row")
    _require(isinstance(raw_row, Mapping), key + " raw row missing")
    _require(raw_row.get("operation") == "sum", key + " operation changed")
    arguments = raw_row.get("arguments")
    _require(arguments == expected["arguments"], key + " arguments changed")
    _require(all(type(item) is int and item >= 0 for item in arguments), key + " arguments malformed")

    scope = packet.get("typed_edge_scope")
    typed = packet.get("typed_edge_refs")
    _require(isinstance(scope, Mapping), key + " typed scope missing")
    _require(isinstance(typed, list), key + " typed refs missing")
    _require(scope.get("operation") == "sum", key + " typed operation changed")
    _require(scope.get("operation_known") is True, key + " operation unknown")
    _require(scope.get("unknown_operation") is False, key + " operation unknown")
    _require(
        scope.get("registered_reference_fields") == list(REFERENCE_KEYS["sum"]),
        key + " registered fields changed",
    )
    _require(scope.get("refs") == typed, key + " duplicate typed-ref views differ")
    _require(scope.get("count") == len(typed) == len(arguments), key + " typed count changed")
    _require(scope.get("unresolved_registered_fields") == [], key + " unresolved typed field")
    _require(scope.get("skipped_or_unresolved") == 0, key + " skipped typed field")
    _require(scope.get("ignored_integer_attributes") == [], key + " guessed integer field")
    _require(scope.get("children_indexed") is False, key + " typed scope indexed children")
    _require(scope.get("edges_followed") is False, key + " typed scope followed edge")

    reviews: list[dict[str, Any]] = []
    corrected: list[dict[str, Any]] = []
    for index, (edge, target) in enumerate(zip(typed, arguments)):
        _require(isinstance(edge, Mapping), key + " edge malformed")
        review, fixed = _edge_review(edge, f"arguments[{index}]", target, f"{key}.edge[{index}]")
        review["row"] = key
        review["edge_index"] = index
        reviews.append(review)
        corrected.append(fixed)

    corrected_scope = copy.deepcopy(dict(scope))
    corrected_scope["refs"] = corrected
    corrected_scope["count"] = len(corrected)

    return {
        "node": node,
        "operation": "sum",
        "row": copy.deepcopy(dict(raw_row)),
        "prior_typed_edge_refs": copy.deepcopy(typed),
        "prior_typed_edge_scope": copy.deepcopy(dict(scope)),
        "typed_edge_refs": corrected,
        "corrected_typed_edge_refs": corrected,
        "typed_edge_scope": corrected_scope,
        "edges_reviewed": reviews,
        "typed_edge_ref_count": len(reviews),
        "children_indexed": False,
        "child_rows_read": 0,
        "graph_edges_followed": False,
        "source_value_selected": False,
    }, reviews


def _validate_q(packet: Mapping[str, Any]) -> tuple[Mapping[str, Any], dict[str, Any], list[dict[str, Any]]]:
    _require(packet.get("schema") == PRIOR_COMPONENT, "inner schema changed")
    _require(packet.get("component") == PRIOR_COMPONENT, "inner component changed")
    _require(
        packet.get("status")
        in {
            "bounded_F_Q_quotient_operands_observed",
            "bounded_F_Q_quotient_operands_observed_prior_unresolved",
        },
        "inner status changed",
    )
    _require(packet.get("source_only") is True, "inner packet is not source-only")
    _require(packet.get("source_rows_only") is True, "inner packet is not row-only")
    scope = packet.get("direct_graph_scope")
    _require(isinstance(scope, Mapping), "direct graph scope missing")
    _require(scope.get("indexed_rows") == [277041, 277054, 25272], "indexed rows changed")
    _require(scope.get("direct_rows") == [277041, 277054, 25272], "direct rows changed")
    for name, expected in (
        ("children_not_indexed", True),
        ("child_rows_read", 0),
        ("graph_edges_followed", False),
        ("source_graph_traversed", False),
        ("protected_nodes_indexed", False),
        ("forbidden_nodes_indexed", False),
    ):
        actual = scope.get(name)
        matches = actual is expected if type(expected) is bool else actual == expected
        _require(matches, "direct scope marker changed: " + name)

    rows = packet.get("rows")
    _require(isinstance(rows, Mapping), "inner rows missing")
    _require(set(rows) == set(EXPECTED_ROWS), "inner row set changed")
    corrected_rows: dict[str, Any] = {}
    all_reviews: list[dict[str, Any]] = []
    for key, expected in EXPECTED_ROWS.items():
        corrected, reviews = _row(key, rows[key], expected)
        corrected_rows[key] = corrected
        all_reviews.extend(reviews)
    _require(len(all_reviews) == 6, "direct edge count changed")
    node_one = [item for item in all_reviews if item["node"] == 1]
    _require(len(node_one) == 1, "opaque node 1 edge count changed")
    _require(node_one[0]["protected_target"] is True, "node 1 was not protected")
    return packet, corrected_rows, all_reviews


def _outer(packet: Mapping[str, Any]) -> Mapping[str, Any]:
    _require(packet.get("schema") == PRIOR_WRAPPER_SCHEMA, "outer schema changed")
    _require(packet.get("component") == PRIOR_COMPONENT, "outer component changed")
    _require(packet.get("status") == "returned", "outer status is not returned")
    _require("error" in packet, "outer error flag missing")
    result = packet.get("result")
    _require(isinstance(result, Mapping), "outer result missing")
    inner = result.get("q_quotient_operands")
    _require(isinstance(inner, Mapping), "outer q_quotient_operands missing")
    return inner


def build_protected_ref_correction(
    root: str | Path = ".",
    source: Mapping[str, Any] | None = None,
    hashes: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Authenticate and correct labels on the six direct operand edges only."""

    _require(isinstance(source, Mapping), "caller frozen manifest is required")
    raw, actual, manifest_count = _authenticate(Path(root).resolve(), source, hashes)
    # No result byte is decoded before every supplied byte has been hashed.
    original = _json(raw[PRIOR_RESULT.as_posix()], PRIOR_RESULT.as_posix())
    inner = _outer(original)
    _inner, corrected_rows, reviews = _validate_q(inner)

    false_to_true = [item for item in reviews if item["previous_false_new_true"]]
    other_edges = [item for item in reviews if not item["previous_false_new_true"]]
    corrected_refs = {
        key: value["typed_edge_refs"] for key, value in corrected_rows.items()
    }
    original_error = copy.deepcopy(original.get("error"))
    return {
        "schema": SCHEMA,
        "component": SCHEMA,
        "status": "protected_ref_correction_returned",
        "source_only": True,
        "direct_dependency_edges_only": True,
        "runner_interface": {
            "entry": "build_protected_ref_correction",
            "arguments": ["root", "source", "hashes"],
            "one_shot": True,
            "coordination": "/root/correction_runner",
        },
        "raw_hash_authentication": {
            "all_manifest_raw_hashes_checked_before_decode": True,
            "manifest_entry_count": manifest_count,
            "required_input": {
                "path": PRIOR_RESULT.as_posix(),
                "expected_sha256": PRIOR_RESULT_SHA256,
                "actual_sha256": actual[PRIOR_RESULT.as_posix()],
            },
        },
        "prior_result": {
            "path": PRIOR_RESULT.as_posix(),
            "sha256": actual[PRIOR_RESULT.as_posix()],
            "wrapper_schema": original.get("schema"),
            "wrapper_component": original.get("component"),
            "wrapper_status": original.get("status"),
            "original_error": original_error,
            "error_present": original_error is not None,
        },
        # The returned producer packet, including its stale flags and error,
        # is retained byte-for-byte as decoded data and is never overwritten.
        "original_result": copy.deepcopy(dict(original)),
        "original_error": original_error,
        "error_preserved": original_error == original.get("error"),
        "protected_node_ids": sorted(PROTECTED_NODE_IDS),
        "relabel_rule": "protected_target = (target node in frozen protected set)",
        "rows": corrected_rows,
        "corrected_typed_edge_refs": corrected_refs,
        "edges_reviewed": reviews,
        "all_direct_edges_reviewed": True,
        "previous_false_new_true": false_to_true,
        "previously_false_now_true": false_to_true,
        "other_edges_reviewed": other_edges,
        "edge_summary": {
            "direct_edge_count": len(reviews),
            "previous_false_new_true_count": len(false_to_true),
            "other_edge_count": len(other_edges),
            "changed_count": sum(item["changed"] for item in reviews),
        },
        "node_1_handling": {
            "opaque_reference_verified": True,
            "indexed": False,
            "edge_followed": False,
            "source_value_selected": False,
            "treated_as_constant": False,
        },
        "archive_read": False,
        "archive_decompressed": False,
        "source_graph_traversed": False,
        "children_indexed": False,
        "graph_edges_followed": False,
        "semantic_proof": _false_flags(),
        "scope_statement": (
            "Only the six typed direct dependency edges in the three returned "
            "sum rows were relabeled. Node 1 remains an opaque reference; no "
            "archive row, child, source value, identity, derivative, positivity, "
            "error removal, or closure was established."
        ),
    }


def _synthetic_edge(field: str, node: int, previous: bool) -> dict[str, Any]:
    return {
        "field": field,
        "typed_edge": True,
        "ref": {
            "status": "opaque_source_ref",
            "namespace": F_GRAPH_NAMESPACE,
            "node": node,
            "source_node_id": node,
            "source_role": "synthetic.direct." + field,
            "node_id_is_not_a_numeric_value": True,
            "source_value_selected": False,
            "indexed": False,
            "protected_node_indexed": False,
        },
        "child_indexed": False,
        "edge_followed": False,
        "protected_target": previous,
    }


def verify_synthetic_fixture() -> dict[str, Any]:
    """Exercise both relabel outcomes without reading any real input."""

    protected, protected_fixed = _edge_review(
        _synthetic_edge("arguments[0]", 1, False),
        "arguments[0]",
        1,
        "synthetic.protected",
    )
    unprotected, unprotected_fixed = _edge_review(
        _synthetic_edge("arguments[1]", 99, False),
        "arguments[1]",
        99,
        "synthetic.unprotected",
    )
    checks = {
        "protected_target_changes_false_to_true": (
            protected["previous_protected_target"] is False
            and protected["protected_target"] is True
            and protected["previous_false_new_true"] is True
        ),
        "unprotected_target_remains_false": (
            unprotected["previous_protected_target"] is False
            and unprotected["protected_target"] is False
            and unprotected["changed"] is False
        ),
        "node_1_stays_opaque_and_unindexed": (
            protected["node"] == 1
            and protected_fixed["ref"]["status"] == "opaque_source_ref"
            and protected_fixed["ref"]["indexed"] is False
            and protected_fixed["child_indexed"] is False
            and protected_fixed["edge_followed"] is False
        ),
        "unprotected_edge_stays_unindexed": (
            unprotected_fixed["protected_target"] is False
            and unprotected_fixed["ref"]["indexed"] is False
            and unprotected_fixed["child_indexed"] is False
        ),
        "node_1_is_not_treated_as_constant": (
            protected["node"] == protected["source_node_id"] == 1
            and protected_fixed["ref"]["source_value_selected"] is False
        ),
    }
    passed = all(checks.values())
    return {
        "schema": SCHEMA + "_synthetic_fixture_v1",
        "component": SCHEMA,
        "status": "synthetic_checks_passed" if passed else "synthetic_checks_failed",
        "passed": passed,
        "checks_passed": passed,
        "synthetic_only": True,
        "checks": checks,
        "fixtures": {
            "protected_edge": {"review": protected, "corrected": protected_fixed},
            "unprotected_edge": {"review": unprotected, "corrected": unprotected_fixed},
        },
        "protected_node_ids": sorted(PROTECTED_NODE_IDS),
        "scope": {
            "synthetic_only": True,
            "archive_read": False,
            "archive_decompressed": False,
            "children_indexed": False,
            "graph_edges_followed": False,
            "node_1_indexed": False,
        },
        "semantic_proof": _false_flags(),
        "real_inputs_read": False,
    }


__all__ = [
    "FGapProtectedRefCorrectionInputError",
    "build_protected_ref_correction",
    "verify_synthetic_fixture",
]
