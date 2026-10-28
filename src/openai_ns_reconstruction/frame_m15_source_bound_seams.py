"""Reference-level signed seam assembly for the authenticated m15 rows.

This module is deliberately source-only.  It authenticates exactly three
submitted report envelopes, decodes only their base64 result artifacts, and
assembles references already present in those results.  It does not read the
compressed archive named by the prior reports, import a project source module,
follow a source graph, evaluate a source expression, or compute a value.

The forty-four deferred occurrences are selected from the published mismatch
result and are taken from the new deferred-result roots.  The remaining twenty
occurrences are retained from the deferred result (sixteen saved roots and four
prior-fresh roots).  Seam cancellation is intentionally strict: a positive and
negative term can pair only when its source namespace and node, field/z role,
and unmodified symbolic coefficient are identical.  In particular, node ``0``
is an ordinary opaque source reference and is never treated as numeric zero.
"""

from __future__ import annotations

import base64
import copy
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
MISMATCH_ARTIFACT = "tmp/frame_m15_endpoint_mismatch_v1/result.json"
MISMATCH_RESULT_SHA256 = (
    "6adc3c1a272b2fcc6f6259c8ad844b8560a9893ac4efdcd8dc689e7bcd177d93"
)

DEFERRED_REPORT_RELATIVE = Path(
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_frame_m15_deferred_endpoint_substitution.json"
)
DEFERRED_REPORT_SHA256 = (
    "65be948ec715056f310e17c827b3132f37f06b6439d932743cfb3b57e2c50d39"
)
DEFERRED_ARTIFACT = "tmp/frame_m15_deferred_endpoint_substitution_v1/result.json"
DEFERRED_RESULT_SHA256 = (
    "bf07cdce2ba58cd9beae53fc576a3ff4d023f56e4a027de495f194cd0047d5f9"
)

AUDIT_REPORT_RELATIVE = Path(
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_frame_m15_deferred_endpoint_substitution_audit.json"
)
AUDIT_REPORT_SHA256 = (
    "bbbf20f63cb0021c6c110a6b173f74e966fb3d719a0e8f6d1e08a1d34848caf9"
)
AUDIT_ARTIFACT = "tmp/frame_m15_deferred_endpoint_substitution_audit_v1/result.json"
AUDIT_RESULT_SHA256 = (
    "3d714d2100ad05b6f37fcca6a3f9e81a6a7da967a2d1bc58998eb1169c2f1566"
)


ENDPOINTS = (
    "main.lower",
    "main.upper",
    "gap.lower",
    "gap.upper",
    "end.lower",
    "end.upper",
    "global.lower",
    "global.upper",
)
FIELDS = (
    "R_squared_q",
    "R_squared_q_squared",
    "R_squared_U_squared",
    "R_q_y",
)
Z_LABELS = ("Z0", "Z1")
FIELD_ALIAS = {
    "R_squared_q": "A",
    "R_squared_q_squared": "B",
    "R_squared_U_squared": "C",
    "R_q_y": "D",
}

SAVED_NAMESPACE = "saved_Rp_native_C5_pulse_source"
PRIOR_FRESH_NAMESPACE = "fresh_frame_gap_substitution_v1"
NEW_NAMESPACE = "fresh_frame_m15_deferred_endpoint_substitution_v1"
DEFERRED_PLACEHOLDER_NAMESPACE = "unexecuted_saved_inlet_substitution_v1"

K_EXPRESSION = "A0/(2*L) + B0/4 - Z^2*C0/L^2 - D0"
K_Z_EXPRESSION = (
    "A1/(2*L) - A0*L_Z/(2*L^2) + B1/4 - Z^2*C1/L^2 "
    "+ (-2*Z/L^2 + 2*Z^2*L_Z/L^3)*C0 - D1"
)
REPORT_K_Z_EXPRESSION = (
    "A1/(2*L) - A0*L_Z/(2*L^2) + B1/4 - "
    "(2*Z*C0 + Z^2*C1)/L^2 + 2*Z^2*C0*L_Z/L^3 - D1"
)
L_EXPRESSION = "1-delta*Z^2"
L_Z_EXPRESSION = "-2*delta*Z"

Z0_TERM_SPECS = (
    ("A0", "R_squared_q", 0, "1/(2*L)"),
    ("B0", "R_squared_q_squared", 0, "1/4"),
    ("C0", "R_squared_U_squared", 0, "-Z^2/L^2"),
    ("D0", "R_q_y", 0, "-1"),
)
Z1_TERM_SPECS = (
    ("A1", "R_squared_q", 1, "1/(2*L)"),
    ("A0", "R_squared_q", 0, "-L_Z/(2*L^2)"),
    ("B1", "R_squared_q_squared", 1, "1/4"),
    ("C1", "R_squared_U_squared", 1, "-Z^2/L^2"),
    ("C0", "R_squared_U_squared", 0, "(-2*Z/L^2 + 2*Z^2*L_Z/L^3)"),
    ("D1", "R_q_y", 1, "-1"),
)

# The conditional identity's signs are kept as endpoint weights.  They are
# not inferred from endpoint names or from source-node numbers.
SEAM_LAYOUT = (
    ("global.lower-minus-main.lower", "global.lower", "main.lower"),
    ("main.upper-minus-gap.lower", "main.upper", "gap.lower"),
    ("gap.upper-minus-end.lower", "gap.upper", "end.lower"),
    ("end.upper-minus-global.upper", "end.upper", "global.upper"),
)
SUM_Z0_EXPRESSION = (
    "+K(global.lower) -K(main.lower) +K(main.upper) -K(gap.lower) "
    "+K(gap.upper) -K(end.lower) +K(end.upper) -K(global.upper)"
)
REPORT_Z1_TERMS_EXPRESSION = (
    "+K_Z(global.lower) -K_Z(main.lower) +K_Z(main.upper) -K_Z(gap.lower) "
    "+K_Z(gap.upper) -K_Z(end.lower) +K_Z(end.upper) -K_Z(global.upper)"
)
SUM_Z1_EXPRESSION = (
    "K_Z(global.lower)-K_Z(main.lower)+K_Z(main.upper)-K_Z(gap.lower) "
    "+K_Z(gap.upper)-K_Z(end.lower)+K_Z(end.upper)-K_Z(global.upper)"
)


def _mapping(value: Any, label: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise ValueError(label + " must be an object")
    return value


def _list(value: Any, label: str) -> list[Any]:
    if not isinstance(value, list):
        raise ValueError(label + " must be a list")
    return value


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _property(value: Any, name: str, label: str) -> Any:
    mapping = _mapping(value, label)
    if name not in mapping:
        raise ValueError(label + " is missing " + name)
    return mapping[name]


def _node_ref(value: Any, label: str) -> dict[str, Any]:
    mapping = _mapping(value, label)
    namespace = mapping.get("namespace")
    node = mapping.get("node")
    _require(isinstance(namespace, str) and bool(namespace), label + " namespace is invalid")
    _require(
        (type(node) is int and node >= 0) or (isinstance(node, str) and bool(node)),
        label + " node is invalid",
    )
    # Keep the node's type and value exactly.  This prevents integer 0 and the
    # string "0" (or two nodes in distinct namespaces) from being conflated.
    return {"namespace": namespace, "node": node}


def _node_key(value: Mapping[str, Any]) -> tuple[str, str, Any]:
    node = value["node"]
    return (str(value["namespace"]), type(node).__name__, node)


def _role_key(endpoint: str, field: str, z_order: int, original_ref: int) -> tuple[str, str, int, int]:
    return endpoint, field, z_order, original_ref


def _read_authenticated_result(
    root: Path,
    report_relative: Path,
    report_sha256: str,
    artifact_name: str,
    result_sha256: str,
    label: str,
) -> tuple[dict[str, Any], dict[str, Any]]:
    """Read one report envelope and authenticate its base64 result artifact."""

    report_raw = (root / report_relative).read_bytes()
    _require(
        hashlib.sha256(report_raw).hexdigest() == report_sha256,
        label + " report SHA-256 mismatch",
    )
    report = _mapping(json.loads(report_raw.decode("utf-8")), label + " report")
    artifacts = _mapping(report.get("artifacts"), label + " artifacts")
    artifact = _mapping(artifacts.get(artifact_name), label + " result artifact")
    encoded = artifact.get("base64")
    _require(isinstance(encoded, str), label + " result artifact has no base64 payload")
    try:
        result_raw = base64.b64decode(encoded, validate=True)
    except Exception as exc:
        raise ValueError(label + " result artifact base64 is invalid") from exc
    _require(
        hashlib.sha256(result_raw).hexdigest() == result_sha256,
        label + " result SHA-256 mismatch",
    )
    metadata_hash = artifact.get("sha256")
    _require(
        metadata_hash in (None, result_sha256),
        label + " result artifact metadata SHA-256 mismatch",
    )
    packet = _mapping(json.loads(result_raw.decode("utf-8")), label + " result packet")
    _require(packet.get("status") == "returned", label + " result packet is not returned")
    result = dict(_mapping(packet.get("result"), label + " result"))
    citation = {
        "report_relative": report_relative.as_posix(),
        "report_sha256": report_sha256,
        "artifact": artifact_name,
        "artifact_sha256": metadata_hash,
        "result_sha256": result_sha256,
        "report_status": report.get("status"),
        "report_component": report.get("component"),
        "report_full_goal_complete": report.get("full_goal_complete"),
        "packet_status": packet.get("status"),
    }
    return result, citation


def _source_row_ref(row: Mapping[str, Any], label: str) -> dict[str, Any] | None:
    status = row.get("status")
    if status == "saved_frame_reference":
        return _node_ref(row, label + " saved source")
    if status == "reused_committed_gap_R_q_y_reference":
        return _node_ref(row, label + " prior fresh source")
    if status == "missing_saved_reference":
        _require(row.get("namespace") == DEFERRED_PLACEHOLDER_NAMESPACE, label + " missing namespace")
        _require("node" not in row, label + " missing row unexpectedly has a node")
        return None
    raise ValueError(label + " has an unknown status")


def _validate_mismatch(result: Mapping[str, Any]) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    _require(result.get("schema") == "frame_m15_endpoint_mismatch_v1", "mismatch schema mismatch")
    _require(
        result.get("status") == "unexecuted_source_bound_endpoint_mismatch",
        "mismatch status mismatch",
    )
    _require(result.get("source_references_are_opaque") is True, "mismatch references are not opaque")
    _require(result.get("endpoint_count") == 8, "mismatch endpoint count mismatch")
    _require(result.get("source_row_count") == 64, "mismatch source row count mismatch")
    _require(result.get("endpoint_layout") == list(ENDPOINTS), "mismatch endpoint layout mismatch")
    _require(result.get("endpoint_fields") == list(FIELDS), "mismatch endpoint fields mismatch")
    _require(result.get("ordinary_Z_rows") == [0, 1], "mismatch Z rows mismatch")
    counts = _mapping(result.get("binding_status_counts"), "mismatch binding status counts")
    _require(counts.get("missing_saved_reference") == 44, "mismatch missing count mismatch")
    _require(counts.get("saved_frame_reference") == 16, "mismatch saved count mismatch")
    _require(
        counts.get("reused_committed_gap_R_q_y_reference") == 4,
        "mismatch prior fresh count mismatch",
    )
    _require(result.get("old_builders_called") is False, "mismatch old builders were called")
    _require(result.get("old_operators_replayed") is False, "mismatch old operators were replayed")
    _require(result.get("source_values_computed") is False, "mismatch source values were computed")
    _require(result.get("numeric_values_computed") is False, "mismatch numeric values were computed")
    _require(result.get("integrals_evaluated") is False, "mismatch integrals were evaluated")

    conditional = _mapping(result.get("conditional_identity"), "mismatch conditional identity")
    z0_terms = _mapping(conditional.get("Z0_terms"), "mismatch Z0 terms")
    z1_terms = _mapping(conditional.get("Z1_terms"), "mismatch Z1 terms")
    _require(z0_terms.get("rational_form") == SUM_Z0_EXPRESSION, "mismatch Z0 expression mismatch")
    _require(z1_terms.get("rational_form") == REPORT_Z1_TERMS_EXPRESSION, "mismatch Z1 terms expression mismatch")
    _require(conditional.get("Z1_rational_form") == SUM_Z1_EXPRESSION, "mismatch Z1 expression mismatch")
    _require(conditional.get("endpoint_equalities_proved") is False, "mismatch endpoint equality claimed")
    _require(conditional.get("residual_zero_proved") is False, "mismatch residual zero claimed")
    _require(result.get("unvalued_K_forms") is not None, "mismatch K forms are missing")

    endpoints = _mapping(result.get("endpoints"), "mismatch endpoints")
    rows: list[dict[str, Any]] = []
    seen_roles: set[tuple[str, str, int, int]] = set()
    for endpoint in ENDPOINTS:
        endpoint_record = _mapping(endpoints.get(endpoint), endpoint + " mismatch endpoint")
        _require(endpoint_record.get("endpoint") == endpoint, endpoint + " endpoint label mismatch")
        fields = _mapping(endpoint_record.get("fields"), endpoint + " fields")
        _require(set(fields) == set(FIELDS), endpoint + " field set mismatch")
        for field in FIELDS:
            field_record = _mapping(fields.get(field), endpoint + "." + field)
            original_refs = _list(
                field_record.get("original_refs_Z0_Z1"), endpoint + "." + field + " original refs"
            )
            _require(len(original_refs) == 2, endpoint + "." + field + " original ref count mismatch")
            rows_map = _mapping(field_record.get("rows"), endpoint + "." + field + " rows")
            for z_order, z_label in enumerate(Z_LABELS):
                row = dict(_mapping(rows_map.get(z_label), endpoint + "." + field + "." + z_label))
                original_ref = row.get("original_ref")
                _require(type(original_ref) is int and original_ref >= 0, "mismatch original_ref is invalid")
                _require(original_ref == original_refs[z_order], "mismatch original_ref list disagrees")
                role = endpoint + "." + field + "." + z_label
                key = _role_key(endpoint, field, z_order, original_ref)
                _require(key not in seen_roles, "duplicate mismatch source role: " + role)
                seen_roles.add(key)
                _source_row_ref(row, role)
                rows.append(
                    {
                        "endpoint": endpoint,
                        "field": field,
                        "z_label": z_label,
                        "z_order": z_order,
                        "role": role,
                        "original_ref": original_ref,
                        "status": row.get("status"),
                        "report_row": copy.deepcopy(row),
                    }
                )
    _require(len(rows) == 64, "mismatch rows do not contain exactly 64 occurrences")

    # Check the report's A/B/C/D aliases against the same 64 rows.  This is a
    # reference-level consistency check; it never follows an alias graph.
    unvalued = _mapping(result.get("unvalued_K_forms"), "mismatch unvalued K forms")
    by_role = {row["role"]: row for row in rows}
    for endpoint in ENDPOINTS:
        form = _mapping(unvalued.get(endpoint), endpoint + " K form")
        _require(form.get("formal_symbol") == "Z", endpoint + " formal symbol mismatch")
        _require(form.get("unvalued") is True, endpoint + " K form is valued")
        weights = _mapping(form.get("weights"), endpoint + " K weights")
        _require(weights.get("L") == L_EXPRESSION, endpoint + " L expression mismatch")
        _require(weights.get("L_Z") == L_Z_EXPRESSION, endpoint + " L_Z expression mismatch")
        _require(_mapping(form.get("delta_source_ref"), endpoint + " delta source") is not None, endpoint + " delta source missing")
        z0_form = _mapping(form.get("Z0"), endpoint + " Z0 K form")
        z1_form = _mapping(form.get("Z1"), endpoint + " Z1 K form")
        _require(z0_form.get("K") == K_EXPRESSION, endpoint + " K expression mismatch")
        _require(z1_form.get("K_Z") == REPORT_K_Z_EXPRESSION, endpoint + " K_Z expression mismatch")
        for z_form, suffix in ((z0_form, "0"), (z1_form, "1")):
            aliases = _mapping(z_form.get("aliases"), endpoint + " aliases " + suffix)
            for alias, field in (("A", "R_squared_q"), ("B", "R_squared_q_squared"), ("C", "R_squared_U_squared"), ("D", "R_q_y")):
                role = endpoint + "." + field + ".Z" + suffix
                alias_record = _mapping(aliases.get(alias + suffix), role + " alias")
                row = by_role[role]
                _require(alias_record.get("original_ref") == row["original_ref"], role + " alias original_ref mismatch")
                _require(alias_record.get("status") == row["status"], role + " alias status mismatch")
                if row["status"] == "missing_saved_reference":
                    _require(alias_record.get("namespace") == DEFERRED_PLACEHOLDER_NAMESPACE, role + " alias namespace mismatch")
                else:
                    _require(
                        _node_key(_node_ref(alias_record, role + " alias source"))
                        == _node_key(_node_ref(row["report_row"], role + " row source")),
                        role + " alias source mismatch",
                    )

    delta_sources = []
    for endpoint in ENDPOINTS:
        form = _mapping(unvalued[endpoint], endpoint + " K form")
        delta_sources.append(_node_ref(form.get("delta_source_ref"), endpoint + " delta source"))
    first_delta = delta_sources[0]
    _require(all(_node_key(item) == _node_key(first_delta) for item in delta_sources), "delta source refs disagree")
    return rows, {"delta_source_ref": copy.deepcopy(first_delta)}


def _validate_deferred(
    result: Mapping[str, Any],
) -> tuple[dict[tuple[str, str, int, int], dict[str, Any]], dict[tuple[str, str, int, int], dict[str, Any]]]:
    _require(result.get("schema") == "frame_m15_deferred_endpoint_substitution_v1", "deferred schema mismatch")
    _require(result.get("status") == "unevaluated_deferred_source_dag", "deferred status mismatch")
    _require(result.get("source_references_are_opaque") is True, "deferred references are not opaque")
    _require(result.get("separate_root_memo_scope") is True, "deferred root scope is not separate")
    _require(result.get("deferred_occurrence_count") == 44, "deferred occurrence count mismatch")
    _require(result.get("fresh_node_count") == 168, "deferred new-node count mismatch")
    _require(result.get("ordinary_Z_rows") == [0, 1], "deferred Z rows mismatch")
    _require(result.get("endpoint_fields") == list(FIELDS), "deferred endpoint fields mismatch")
    _require(result.get("old_builders_called") is False, "deferred old builders were called")
    _require(result.get("old_constructors_called") is False, "deferred old constructors were called")
    _require(result.get("old_operators_replayed") is False, "deferred old operators were replayed")
    _require(result.get("old_twenty_roots_replayed") is False, "deferred old roots were replayed")
    _require(result.get("source_values_computed") is False, "deferred source values were computed")
    _require(result.get("numeric_values_computed") is False, "deferred numeric values were computed")
    _require(result.get("integrals_evaluated") is False, "deferred integrals were evaluated")
    _require(result.get("old_error_terms_removed") is False, "deferred old errors were removed")
    _require(result.get("extra_ledger_terms_added") is False, "deferred extra terms were added")

    boundary = _mapping(result.get("boundary_cache"), "deferred boundary cache")
    _require(boundary.get("saved_namespace") == SAVED_NAMESPACE, "deferred saved namespace mismatch")
    _require(boundary.get("prior_fresh_namespace") == PRIOR_FRESH_NAMESPACE, "deferred prior namespace mismatch")
    _require(boundary.get("new_namespace") == NEW_NAMESPACE, "deferred new namespace mismatch")
    _require(boundary.get("saved_semantics") == "opaque authenticated archive boundary; not followed", "deferred saved boundary semantics changed")
    _require(boundary.get("prior_fresh_semantics") == "opaque audited structural boundary; not re-expanded", "deferred prior boundary semantics changed")

    occurrences = _list(result.get("occurrences"), "deferred occurrences")
    _require(len(occurrences) == 44, "deferred occurrences list count mismatch")
    deferred_by_role: dict[tuple[str, str, int, int], dict[str, Any]] = {}
    occurrence_ids: set[str] = set()
    for raw_value in occurrences:
        raw = dict(_mapping(raw_value, "deferred occurrence"))
        endpoint = raw.get("endpoint")
        field = raw.get("field")
        z_order = raw.get("z_order")
        z_label = raw.get("z_label")
        original_ref = raw.get("original_ref")
        _require(endpoint in ENDPOINTS, "deferred occurrence endpoint is unknown")
        _require(field in FIELDS, "deferred occurrence field is unknown")
        _require(type(z_order) is int and z_order in (0, 1), "deferred occurrence z_order is invalid")
        _require(z_label == Z_LABELS[z_order], "deferred occurrence z_label mismatch")
        _require(type(original_ref) is int and original_ref >= 0, "deferred occurrence original_ref is invalid")
        role = endpoint + "." + field + "." + z_label
        _require(raw.get("role") == role, "deferred occurrence role mismatch")
        key = _role_key(endpoint, field, z_order, original_ref)
        _require(key not in deferred_by_role, "duplicate deferred occurrence role: " + role)
        occurrence_id = raw.get("occurrence_id")
        _require(isinstance(occurrence_id, str) and occurrence_id, "deferred occurrence id is invalid")
        _require(occurrence_id not in occurrence_ids, "duplicate deferred occurrence id")
        occurrence_ids.add(occurrence_id)
        source_ref = _node_ref(raw.get("source_ref"), role + " source_ref")
        root_ref = _node_ref(raw.get("root_ref"), role + " root_ref")
        _require(_node_key(source_ref) == _node_key(root_ref), role + " source/root mismatch")
        root_status = raw.get("root_status")
        if root_status == "new_fresh_dag":
            _require(source_ref["namespace"] == NEW_NAMESPACE, role + " new namespace mismatch")
        elif root_status == "prior_fresh_boundary":
            _require(source_ref["namespace"] == PRIOR_FRESH_NAMESPACE, role + " prior boundary namespace mismatch")
        else:
            raise ValueError(role + " unexpected root status")
        _require(raw.get("dag_namespace") == NEW_NAMESPACE, role + " DAG namespace mismatch")
        _require(raw.get("report_status") == "missing_saved_reference", role + " deferred status mismatch")
        _require(raw.get("semantic_identity_proved") is False, role + " semantic identity was claimed")
        _require(raw.get("seam_identity_proved") is False, role + " seam identity was claimed")
        deferred_by_role[key] = {
            "record": copy.deepcopy(raw),
            "source_ref": source_ref,
            "role": role,
            "endpoint": endpoint,
            "field": field,
            "z_label": z_label,
            "z_order": z_order,
            "original_ref": original_ref,
        }

    root_counts = {"new_fresh_dag": 0, "prior_fresh_boundary": 0}
    for entry in deferred_by_role.values():
        root_counts[entry["record"]["root_status"]] += 1
    _require(root_counts == {"new_fresh_dag": 32, "prior_fresh_boundary": 12}, "deferred root counts mismatch")

    retained = _list(result.get("retained_occurrences"), "deferred retained occurrences")
    _require(len(retained) == 20, "deferred retained occurrence count mismatch")
    retained_by_role: dict[tuple[str, str, int, int], dict[str, Any]] = {}
    retained_status_counts = {"saved_frame_reference": 0, "reused_committed_gap_R_q_y_reference": 0}
    for raw_value in retained:
        raw = dict(_mapping(raw_value, "retained occurrence"))
        endpoint = raw.get("endpoint")
        field = raw.get("field")
        z_order = raw.get("z_order")
        z_label = raw.get("z_label")
        original_ref = raw.get("original_ref")
        _require(endpoint in ENDPOINTS, "retained occurrence endpoint is unknown")
        _require(field in FIELDS, "retained occurrence field is unknown")
        _require(type(z_order) is int and z_order in (0, 1), "retained occurrence z_order is invalid")
        _require(z_label == Z_LABELS[z_order], "retained occurrence z_label mismatch")
        _require(type(original_ref) is int and original_ref >= 0, "retained occurrence original_ref is invalid")
        role = endpoint + "." + field + "." + z_label
        _require(raw.get("role") in (None, role), "retained occurrence role disagrees")
        key = _role_key(endpoint, field, z_order, original_ref)
        _require(key not in retained_by_role, "duplicate retained occurrence role: " + role)
        status = raw.get("status")
        _require(status in retained_status_counts, "retained occurrence status is unknown")
        retained_status_counts[status] += 1
        source_ref = _node_ref(raw.get("source_ref"), role + " retained source_ref")
        if status == "saved_frame_reference":
            _require(source_ref["namespace"] == SAVED_NAMESPACE, role + " saved namespace mismatch")
        else:
            _require(source_ref["namespace"] == PRIOR_FRESH_NAMESPACE, role + " prior fresh namespace mismatch")
            _require(field == "R_q_y", role + " prior fresh field mismatch")
        retained_by_role[key] = {
            "record": copy.deepcopy(raw),
            "source_ref": source_ref,
            "role": role,
            "endpoint": endpoint,
            "field": field,
            "z_label": z_label,
            "z_order": z_order,
            "original_ref": original_ref,
        }
    _require(retained_status_counts["saved_frame_reference"] == 16, "retained saved count mismatch")
    _require(retained_status_counts["reused_committed_gap_R_q_y_reference"] == 4, "retained prior count mismatch")

    provenance = _mapping(result.get("provenance"), "deferred provenance")
    _require(provenance.get("mismatch_report_sha256") == MISMATCH_REPORT_SHA256, "deferred mismatch report citation changed")
    _require(provenance.get("mismatch_result_sha256") == MISMATCH_RESULT_SHA256, "deferred mismatch result citation changed")
    _require(provenance.get("saved_rewritten_nodes_count") == 7752, "deferred saved history count changed")
    return deferred_by_role, retained_by_role


def _validate_audit(result: Mapping[str, Any]) -> None:
    _require(result.get("schema") == "frame_m15_deferred_endpoint_substitution_audit_v1", "audit schema mismatch")
    _require(result.get("status") == "passed", "audit status is not passed")
    _require(result.get("failure_count") == 0, "audit reports failures")
    _require(result.get("failure_details") == [], "audit failure details are not empty")
    _require(result.get("failure_details_truncated") is False, "audit failure details are truncated")
    _require(result.get("expected_result_sha256") == DEFERRED_RESULT_SHA256, "audit expected result citation changed")
    _require(result.get("result_sha256") == DEFERRED_RESULT_SHA256, "audit result citation changed")
    counts = _mapping(result.get("counts"), "audit counts")
    for key, expected in (("deferred_occurrences_checked", 44), ("retained_occurrences_checked", 20), ("new_nodes_checked", 168), ("new_operand_refs_checked", 408)):
        _require(counts.get(key) == expected, "audit coverage changed: " + key)
    flags = _mapping(result.get("semantic_flags"), "audit semantic flags")
    for key, value in flags.items():
        _require(value is False, "audit semantic flag is true: " + str(key))
    for key in (
        "physical_equalities_authenticated",
        "physical_equivalence_proved",
        "seam_identity_proved",
        "semantic_identity_proved",
        "Z1_commutation_proved",
        "ordinary_Z1_commutation_proved",
        "G_equals_F_proved",
        "P_equals_F_proved",
    ):
        _require(result.get(key) is False, "audit physical/semantic claim is true: " + key)


def _assemble_occurrences(
    rows: list[dict[str, Any]],
    deferred_by_role: Mapping[tuple[str, str, int, int], Mapping[str, Any]],
    retained_by_role: Mapping[tuple[str, str, int, int], Mapping[str, Any]],
    mismatch_citation: Mapping[str, Any],
    deferred_citation: Mapping[str, Any],
) -> list[dict[str, Any]]:
    assembled: list[dict[str, Any]] = []
    selected_keys: set[tuple[str, str, int, int]] = set()
    expected_missing = 44
    expected_saved = 16
    expected_prior = 4
    counts = {"missing_saved_reference": 0, "saved_frame_reference": 0, "reused_committed_gap_R_q_y_reference": 0}
    for row in rows:
        key = _role_key(row["endpoint"], row["field"], row["z_order"], row["original_ref"])
        status = row["status"]
        if status == "missing_saved_reference":
            selected = deferred_by_role.get(key)
            selection = "deferred_result.occurrences"
        else:
            selected = retained_by_role.get(key)
            selection = "deferred_result.retained_occurrences"
        _require(selected is not None, "no authenticated replacement for source role: " + row["role"])
        _require(key not in selected_keys, "replacement source role selected twice: " + row["role"])
        selected_keys.add(key)
        selected_status = selected["record"].get("report_status", selected["record"].get("status"))
        _require(selected_status == status, row["role"] + " replacement status mismatch")
        source_ref = _node_ref(selected["source_ref"], row["role"] + " selected source")
        counts[status] += 1
        assembled.append(
            {
                "role": row["role"],
                "endpoint": row["endpoint"],
                "field": row["field"],
                "z_label": row["z_label"],
                "z_order": row["z_order"],
                "original_ref": row["original_ref"],
                "binding_status": status,
                "source_ref": source_ref,
                "source_namespace": source_ref["namespace"],
                "source_node": source_ref["node"],
                "source_selection": selection,
                "source_provenance": {
                    "mismatch_report": {
                        "report_relative": mismatch_citation["report_relative"],
                        "report_sha256": mismatch_citation["report_sha256"],
                        "result_sha256": mismatch_citation["result_sha256"],
                    },
                    "replacement_report": {
                        "report_relative": deferred_citation["report_relative"],
                        "report_sha256": deferred_citation["report_sha256"],
                        "result_sha256": deferred_citation["result_sha256"],
                    },
                    "replacement_record": selection,
                },
                # Both records are retained so that a later audit can recover
                # the exact source role and its submitted historical origin.
                "mismatch_row": copy.deepcopy(row["report_row"]),
                "replacement_record": copy.deepcopy(selected["record"]),
            }
        )
    _require(len(assembled) == 64, "assembled occurrence count is not 64")
    _require(counts["missing_saved_reference"] == expected_missing, "assembled missing count mismatch")
    _require(counts["saved_frame_reference"] == expected_saved, "assembled saved count mismatch")
    _require(counts["reused_committed_gap_R_q_y_reference"] == expected_prior, "assembled prior count mismatch")
    _require(set(selected_keys) == set(deferred_by_role) | set(retained_by_role), "replacement history was not preserved exactly")
    return assembled


def _occurrence_index(occurrences: list[dict[str, Any]]) -> dict[tuple[str, str, int], dict[str, Any]]:
    index: dict[tuple[str, str, int], dict[str, Any]] = {}
    for occurrence in occurrences:
        key = (occurrence["endpoint"], occurrence["field"], occurrence["z_order"])
        _require(key not in index, "assembled occurrence role is not unique: " + occurrence["role"])
        index[key] = occurrence
    return index


def _signed_weighted_coefficient(weight: int, kernel_coefficient: str) -> str:
    sign = "+1" if weight == 1 else "-1"
    return sign + "*" + "(" + kernel_coefficient + ")"


def _build_term(
    seam_id: str,
    endpoint_weight: int,
    endpoint: str,
    z_label: str,
    component: str,
    field: str,
    z_order: int,
    kernel_coefficient: str,
    occurrence_index: Mapping[tuple[str, str, int], Mapping[str, Any]],
) -> dict[str, Any]:
    occurrence = occurrence_index[(endpoint, field, z_order)]
    source_ref = _node_ref(occurrence["source_ref"], component + " source")
    return {
        "seam_id": seam_id,
        "endpoint": endpoint,
        "endpoint_weight": endpoint_weight,
        "z_label": z_label,
        "component": component,
        "field": field,
        "z_order": z_order,
        "source_role": occurrence["role"],
        "original_ref": occurrence["original_ref"],
        "source_ref": source_ref,
        "source_namespace": source_ref["namespace"],
        "source_node": source_ref["node"],
        "binding_status": occurrence["binding_status"],
        "kernel_coefficient": kernel_coefficient,
        "weighted_coefficient": _signed_weighted_coefficient(endpoint_weight, kernel_coefficient),
    }


def _pair_terms(positive: list[dict[str, Any]], negative: list[dict[str, Any]]) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    used_negative: set[int] = set()
    cancellations: list[dict[str, Any]] = []
    for positive_index, positive_term in enumerate(positive):
        positive_key = (
            positive_term["field"],
            positive_term["z_order"],
            positive_term["kernel_coefficient"],
            _node_key(positive_term["source_ref"]),
        )
        for negative_index, negative_term in enumerate(negative):
            if negative_index in used_negative:
                continue
            negative_key = (
                negative_term["field"],
                negative_term["z_order"],
                negative_term["kernel_coefficient"],
                _node_key(negative_term["source_ref"]),
            )
            if positive_key != negative_key:
                continue
            _require(positive_term["endpoint_weight"] == 1, "positive seam term has wrong weight")
            _require(negative_term["endpoint_weight"] == -1, "negative seam term has wrong weight")
            used_negative.add(negative_index)
            cancellations.append(
                {
                    "field": positive_term["field"],
                    "z_label": positive_term["z_label"],
                    "z_order": positive_term["z_order"],
                    "component": positive_term["component"],
                    "kernel_coefficient": positive_term["kernel_coefficient"],
                    "namespace": positive_term["source_ref"]["namespace"],
                    "node": positive_term["source_ref"]["node"],
                    "positive_source_role": positive_term["source_role"],
                    "negative_source_role": negative_term["source_role"],
                    "positive_original_ref": positive_term["original_ref"],
                    "negative_original_ref": negative_term["original_ref"],
                    "positive_weighted_coefficient": positive_term["weighted_coefficient"],
                    "negative_weighted_coefficient": negative_term["weighted_coefficient"],
                    "match_rule": "exact namespace/node + exact field/Z order + exact kernel coefficient",
                }
            )
            break
    # Use stable identities rather than source-node keys when removing terms;
    # repeated node 0 terms must remain distinct until an exact pair is found.
    paired_positive_roles = {item["positive_source_role"] for item in cancellations}
    paired_negative_roles = {item["negative_source_role"] for item in cancellations}
    remaining = [
        copy.deepcopy(term)
        for term in positive
        if term["source_role"] not in paired_positive_roles
    ] + [
        copy.deepcopy(term)
        for index, term in enumerate(negative)
        if index not in used_negative and term["source_role"] not in paired_negative_roles
    ]
    return cancellations, remaining


def _build_seam(
    seam_id: str,
    positive_endpoint: str,
    negative_endpoint: str,
    occurrence_index: Mapping[tuple[str, str, int], Mapping[str, Any]],
) -> dict[str, Any]:
    output: dict[str, Any] = {
        "seam_id": seam_id,
        "positive_endpoint": positive_endpoint,
        "negative_endpoint": negative_endpoint,
        "positive_weight": 1,
        "negative_weight": -1,
        "delta_source_ref": None,
    }
    for z_label, specs, kernel, expression in (
        ("Z0", Z0_TERM_SPECS, "K", "+K(" + positive_endpoint + ") -K(" + negative_endpoint + ")"),
        ("Z1", Z1_TERM_SPECS, "K_Z", "K_Z(" + positive_endpoint + ")-K_Z(" + negative_endpoint + ")"),
    ):
        positive_terms = [
            _build_term(
                seam_id,
                1,
                positive_endpoint,
                z_label,
                component,
                field,
                z_order,
                coefficient,
                occurrence_index,
            )
            for component, field, z_order, coefficient in specs
        ]
        negative_terms = [
            _build_term(
                seam_id,
                -1,
                negative_endpoint,
                z_label,
                component,
                field,
                z_order,
                coefficient,
                occurrence_index,
            )
            for component, field, z_order, coefficient in specs
        ]
        cancellations, remaining = _pair_terms(positive_terms, negative_terms)
        all_terms = copy.deepcopy(positive_terms) + copy.deepcopy(negative_terms)
        output[z_label] = {
            "kernel": kernel,
            "expression": expression,
            "expression_terms": all_terms,
            "expression_total_terms": len(all_terms),
            "cancellations": cancellations,
            "cancellation_count": len(cancellations),
            "remaining_terms": remaining,
            "remaining_term_count": len(remaining),
            "pairing_is_source_reference_only": True,
        }
    return output


def _aggregate_sum(seams: list[dict[str, Any]], z_label: str, expression: str) -> dict[str, Any]:
    terms: list[dict[str, Any]] = []
    cancellations: list[dict[str, Any]] = []
    remaining: list[dict[str, Any]] = []
    for seam in seams:
        part = _mapping(seam[z_label], seam["seam_id"] + " " + z_label)
        terms.extend(copy.deepcopy(_list(part["expression_terms"], "seam expression terms")))
        cancellations.extend(copy.deepcopy(_list(part["cancellations"], "seam cancellations")))
        remaining.extend(copy.deepcopy(_list(part["remaining_terms"], "seam remaining terms")))
    return {
        "expression": expression,
        "kernel": "K" if z_label == "Z0" else "K_Z",
        "expression_terms": terms,
        "expression_total_terms": len(terms),
        "cancellations": cancellations,
        "cancellation_count": len(cancellations),
        "remaining_terms": remaining,
        "remaining_term_count": len(remaining),
        "pairing_is_source_reference_only": True,
    }


def recover(root: Path) -> dict[str, Any]:
    """Return authenticated source-bound seam references without execution."""

    root = Path(root)
    mismatch_result, mismatch_citation = _read_authenticated_result(
        root,
        MISMATCH_REPORT_RELATIVE,
        MISMATCH_REPORT_SHA256,
        MISMATCH_ARTIFACT,
        MISMATCH_RESULT_SHA256,
        "mismatch",
    )
    deferred_result, deferred_citation = _read_authenticated_result(
        root,
        DEFERRED_REPORT_RELATIVE,
        DEFERRED_REPORT_SHA256,
        DEFERRED_ARTIFACT,
        DEFERRED_RESULT_SHA256,
        "deferred substitution",
    )
    audit_result, audit_citation = _read_authenticated_result(
        root,
        AUDIT_REPORT_RELATIVE,
        AUDIT_REPORT_SHA256,
        AUDIT_ARTIFACT,
        AUDIT_RESULT_SHA256,
        "deferred audit",
    )

    rows, mismatch_context = _validate_mismatch(mismatch_result)
    deferred_by_role, retained_by_role = _validate_deferred(deferred_result)
    _validate_audit(audit_result)
    occurrences = _assemble_occurrences(
        rows,
        deferred_by_role,
        retained_by_role,
        mismatch_citation,
        deferred_citation,
    )
    occurrence_index = _occurrence_index(occurrences)
    seams = [
        _build_seam(seam_id, positive_endpoint, negative_endpoint, occurrence_index)
        for seam_id, positive_endpoint, negative_endpoint in SEAM_LAYOUT
    ]
    sum_local_global = {
        "Z0": _aggregate_sum(seams, "Z0", SUM_Z0_EXPRESSION),
        "Z1": _aggregate_sum(seams, "Z1", SUM_Z1_EXPRESSION),
    }
    for seam in seams:
        seam["delta_source_ref"] = copy.deepcopy(mismatch_context["delta_source_ref"])

    report_citations = {
        "mismatch": copy.deepcopy(mismatch_citation),
        "deferred_substitution": copy.deepcopy(deferred_citation),
        "deferred_audit": copy.deepcopy(audit_citation),
    }
    return {
        "schema": "frame_m15_source_bound_seams_v1",
        "status": "unevaluated_source_bound_signed_seams",
        "source_only": True,
        "source_executed": False,
        "source_references_are_opaque": True,
        "graph_traversal_performed": False,
        "source_values_computed": False,
        "numeric_values_computed": False,
        "integrals_evaluated": False,
        "old_builders_called": False,
        "old_constructors_called": False,
        "old_operators_replayed": False,
        "old_error_terms_removed": False,
        "extra_ledger_terms_added": False,
        "ordinary_Z_rows": [0, 1],
        "occurrence_count": len(occurrences),
        "binding_counts": {
            "resolved_deferred_occurrences": 44,
            "new_fresh_dag_roots": 32,
            "prior_fresh_boundary_roots": 12,
            "saved_frame_reference_retained": 16,
            "prior_fresh_reference_retained": 4,
        },
        "field_aliases": copy.deepcopy(FIELD_ALIAS),
        "algebra": {
            "K": K_EXPRESSION,
            "K_Z": K_Z_EXPRESSION,
            "report_K_Z": REPORT_K_Z_EXPRESSION,
            "L": L_EXPRESSION,
            "L_Z": L_Z_EXPRESSION,
            "delta_source_ref": copy.deepcopy(mismatch_context["delta_source_ref"]),
            "term_coefficients": {
                "Z0": [item[3] for item in Z0_TERM_SPECS],
                "Z1": [item[3] for item in Z1_TERM_SPECS],
            },
        },
        "authenticated_reports": report_citations,
        # Keep the authenticated result snapshots as history.  No result is
        # silently replaced by the newly assembled references.
        "historical_inputs": {
            "mismatch_result": copy.deepcopy(mismatch_result),
            "deferred_substitution_result": copy.deepcopy(deferred_result),
            "deferred_audit_result": copy.deepcopy(audit_result),
        },
        "occurrences": occurrences,
        "seams": {seam["seam_id"]: seam for seam in seams},
        "sum_local_global": sum_local_global,
        "conclusion": {
            "source_role_exact_match_only": True,
            "physical_closure_claimed": False,
            "source_semantic_equivalence_proved": False,
            "same_physical_domain_proved": False,
            "Z1_exchange_proved": False,
            "old_error_deletion_claimed": False,
            "numeric_error_reduction_claimed": False,
            "statement": (
                "The reported cancellations and remaining terms are reference-level "
                "signed algebra only; they do not establish physical closure, source "
                "semantic equivalence, Z1 exchange, or numerical error reduction."
            ),
        },
        "limitations": [
            "Only the three authenticated report envelopes are read.",
            "The compressed archive named by prior reports is not read.",
            "No source module is imported and no source graph is traversed.",
            "No source expression, integral, endpoint equality, or numeric value is evaluated.",
            "Cancellation requires exact namespace/node, field/Z role, and symbolic coefficient.",
        ],
    }


__all__ = ["recover"]
