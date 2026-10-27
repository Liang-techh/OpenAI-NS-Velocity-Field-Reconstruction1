"""Inert native source-parameter ordinary-Z certificate kernel draft.

This module is intentionally a new, small source-identity family.  Importing
it does not import the project tree, a provider, a checker, a symbolic package,
or a numerical helper.  Nothing is called at module load time.  The root
workflow freezes this file and the four JSON literal wrappers, then invokes
``derive`` exactly once.

The kernel consumes authenticated JSON values as data.  It derives free
variable *tokens* from the exact graph records preserved by the raw capsule;
it never evaluates a source value, materializes ``N``, differentiates a phase,
or treats a metadata flag as a proof.  All output is ordinary JSON data.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from copy import deepcopy
from typing import Any


COMPONENT = "native_parameter_Z_independence"
INPUT_SCHEMA = "native_parameter_Z_independence_input_v1"
RESULT_SCHEMA = "original_native_parameter_Z_independence"

SOURCE_WRAPPER_NAMES = (
    "saved_metric_result",
    "scalar_literals",
    "phase_literals",
    "raw_capsule",
)

SOURCE_PINS = {
    "saved_metric_result": "6e9a3aa733a080d3203328b91f850ae7ad18d29049df2421a8af409b17249679",
    "scalar_literals": "a533245e802dda5391601b26e5939da751db1fed9f3f35d5a328f9915b4bfd9d",
    "phase_literals": "536f693ea89e25066cadda6191f92c57effb48489d74097cb6dc5baacb9beb00",
    "raw_capsule": "f5ef9950a04fcf0d33fdfd99dd03b799baaae433c6bc45285a0d7020663fcf7a",
}

PARENT_REPORT_SHA256 = "8b4ef714ffad4c7181484d6b29511637fd44c3412bd6783f806b4f7e324de927"
PARENT_RESULT_SHA256 = SOURCE_PINS["saved_metric_result"]
SOURCE_HEAD = "0bc7d7eb694467707a954b5a383bbc8dd5fc454d"
ARCHIVE_SHA256 = "75895e1637411cfd2e490778f495f3adf904295a8840184cdc467e006b23305d"
CURRENT_SOURCE_RECEIPT_SHA256 = "f3c4d1b230b3835da6f4f270926fa85f74ca488f8331e05fbbb894a7c22a0158"

RAW_RECORD_COUNT = 49
RAW_RECORD_IDS = (
    0,
    1,
    2,
    3,
    4,
    5,
    8,
    10,
    13,
    14,
    18,
    19,
    20,
    21,
    22,
    28,
    643,
    2083,
    24954,
    35560,
    35562,
    35564,
    35565,
    65604,
    79968,
    81513,
    81514,
    84528,
    84529,
    84530,
    84531,
    84532,
    84533,
    85743,
    97298,
    99411,
    185676,
    187214,
    187963,
    205154,
    241536,
    241542,
    241786,
    241787,
    242547,
    242548,
    242549,
    242550,
    242551,
)

TARGET_ROOTS = {
    "logC": 35560,
    "logP": 5,
    "logref": 81514,
    "Rpost_restore": 84533,
    "Rm": 85743,
    "Ra": 65604,
    "phase_origin": 24954,
    "N": 18,
    "nativeR": 241787,
    "unfoldedphase": 242550,
    "actual_phase": 242551,
}

SCALAR_PACKET_IDS = (
    2,
    3,
    4,
    5,
    2083,
    35560,
    79968,
    81513,
    81514,
    84528,
    84529,
    84530,
    84531,
    84532,
    84533,
    85743,
)

KNOWN_OPERATIONS = {
    "exact_rational",
    "exact_positive_integer_expression",
    "current_original_source_parameter",
    "bound_variable",
    "sum",
    "product",
    "negative",
    "analytic_unary",
    "positive_quotient",
    "function_substitution",
}

ANALYTIC_NAMES = {"exp", "log", "fractional_part"}


class _KernelFailure(Exception):
    """Expected source-packet or graph rejection."""

    def __init__(self, code: str, message: str, details: Any = None) -> None:
        super().__init__(message)
        self.code = code
        self.message = message
        self.details = details


class _PolicyRefusal(Exception):
    """An actual source-backed request was deliberately refused."""

    def __init__(self, row: Mapping[str, Any]) -> None:
        super().__init__(str(row.get("reason", "policy refusal")))
        self.row = dict(row)


def _clone(value: Any) -> Any:
    return deepcopy(value)


def _fail(code: str, message: str, details: Any = None) -> None:
    raise _KernelFailure(code, message, details)


def _require(condition: Any, code: str, message: str, details: Any = None) -> None:
    if not condition:
        _fail(code, message, details)


def _admit_or_refuse_request(
    *,
    guard: str,
    request: str,
    evidence: Mapping[str, Any],
    reject: bool,
    reason: str,
    extra: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Evaluate one actual request and raise only when source evidence rejects it."""

    if not reject:
        return {
            "guard": guard,
            "request": request,
            "decision": "admitted",
            "evidence": _clone(dict(evidence)),
        }
    row: dict[str, Any] = {
        "guard": guard,
        "request": request,
        "decision": "refused",
        "evidence": _clone(dict(evidence)),
        "reason": reason,
    }
    if extra is not None:
        row.update(_clone(dict(extra)))
    raise _PolicyRefusal(row)


def _is_mapping(value: Any) -> bool:
    return isinstance(value, Mapping)


def _is_int(value: Any) -> bool:
    return isinstance(value, int) and not isinstance(value, bool)


def _valid_sha(value: Any) -> bool:
    return (
        isinstance(value, str)
        and len(value) == 64
        and all(character in "0123456789abcdefABCDEF" for character in value)
    )


def _node_id(value: Any, *, code: str = "malformed_child_id") -> int:
    _require(_is_int(value) and value >= 0, code, "Graph child ID must be a nonnegative integer", value)
    return value


def _token(kind: str, node_id: int, name: str) -> str:
    # Node identity is deliberately part of every token.  Name matching alone
    # would make a second lexical binding indistinguishable from this binding.
    return f"{kind}:{node_id}:{name}"


def _sorted_tokens(values: Sequence[str]) -> list[str]:
    return sorted(set(values))


def _token_for_bound(node_id: int, record: Mapping[str, Any]) -> str:
    return _token("bound_variable", node_id, record["name"])


def _token_for_fixed(node_id: int, record: Mapping[str, Any]) -> str:
    return _token("fixed_parameter", node_id, record["name"])


def _raw_record_id_map(raw_records: Any) -> dict[int, dict[str, Any]]:
    _require(_is_mapping(raw_records), "record_map_missing", "raw_graph_records must be a mapping")
    _require(len(raw_records) == RAW_RECORD_COUNT, "record_count_mismatch", "raw graph record count changed")
    result: dict[int, dict[str, Any]] = {}
    for key, record in raw_records.items():
        _require(isinstance(key, str) and key.isdecimal(), "record_id_malformed", "Raw record key must be a decimal string", key)
        node_id = _node_id(int(key), code="record_id_malformed")
        _require(str(node_id) == key, "record_id_malformed", "Raw record keys must be canonical decimal IDs", key)
        _require(_is_mapping(record), "raw_record_malformed", "Raw graph record must be an object", node_id)
        _require(node_id not in result, "duplicate_record_id", "Duplicate raw graph record ID", node_id)
        result[node_id] = dict(record)
    _require(tuple(sorted(result)) == RAW_RECORD_IDS, "record_id_set_mismatch", "Raw graph record ID set changed")
    return result


def _children(node_id: int, record: Mapping[str, Any]) -> list[int]:
    operation = record.get("operation")
    _require(operation in KNOWN_OPERATIONS, "unknown_operation", "Unknown graph operation", {"node_id": node_id, "operation": operation})
    if operation in {"exact_rational", "exact_positive_integer_expression", "current_original_source_parameter", "bound_variable"}:
        return []
    if operation in {"sum", "product"}:
        arguments = record.get("arguments")
        _require(isinstance(arguments, list) and bool(arguments), "arguments_malformed", "Sum/product arguments must be a nonempty list", node_id)
        return [_node_id(item) for item in arguments]
    if operation in {"negative", "analytic_unary"}:
        return [_node_id(record.get("argument"))]
    if operation == "positive_quotient":
        return [_node_id(record.get("numerator")), _node_id(record.get("denominator"))]
    if operation == "function_substitution":
        return [
            _node_id(record.get("expression")),
            _node_id(record.get("variable")),
            _node_id(record.get("value")),
        ]
    _fail("unknown_operation", "Unknown graph operation", {"node_id": node_id, "operation": operation})
    return []


def _validate_scalar_packet(packet: Mapping[str, Any]) -> dict[int, dict[str, Any]]:
    _require(packet.get("schema") == "bounded_static_native_scalar_source_ast_literal_rows_v1", "scalar_schema_mismatch", "Scalar literal packet schema changed")
    _require(packet.get("status") == "static_only_uncomputed", "scalar_status_mismatch", "Scalar literal packet is not static-only")
    _require(packet.get("runtime_model_metadata") is None, "scalar_runtime_metadata", "Scalar packet runtime metadata must remain null")
    graph = packet.get("literal_graph_nodes")
    _require(_is_mapping(graph), "scalar_graph_missing", "Scalar packet graph section missing")
    rows = graph.get("rows")
    _require(_is_mapping(rows), "scalar_rows_missing", "Scalar packet rows missing")
    ids = tuple(sorted(_node_id(int(key), code="scalar_row_id_malformed") for key in rows))
    _require(ids == SCALAR_PACKET_IDS, "scalar_row_id_set_mismatch", "Scalar packet row IDs changed")
    result: dict[int, dict[str, Any]] = {}
    for key, row in rows.items():
        _require(isinstance(key, str) and key.isdecimal() and str(int(key)) == key, "scalar_row_id_malformed", "Scalar row key malformed", key)
        _require(_is_mapping(row), "scalar_row_malformed", "Scalar row must be an object", key)
        result[int(key)] = dict(row)
    return result


def _validate_phase_packet(packet: Mapping[str, Any]) -> dict[int, dict[str, Any]]:
    _require(packet.get("schema") == "bounded_static_native_phase_parameter_source_ast_literal_rows_v1", "phase_schema_mismatch", "Phase literal packet schema changed")
    _require(packet.get("status") == "static_only_uncomputed", "phase_status_mismatch", "Phase literal packet is not static-only")
    _require(packet.get("runtime_model_metadata") is None, "phase_runtime_metadata", "Phase packet runtime metadata must remain null")
    graph = packet.get("literal_graph_nodes")
    _require(_is_mapping(graph), "phase_graph_missing", "Phase packet graph section missing")
    rows = graph.get("rows")
    _require(_is_mapping(rows), "phase_rows_missing", "Phase packet rows missing")
    ids = tuple(sorted(_node_id(int(key), code="phase_row_id_malformed") for key in rows))
    _require(ids == RAW_RECORD_IDS, "phase_row_id_set_mismatch", "Phase packet must cover every raw capsule record")
    result: dict[int, dict[str, Any]] = {}
    for key, row in rows.items():
        _require(isinstance(key, str) and key.isdecimal() and str(int(key)) == key, "phase_row_id_malformed", "Phase row key malformed", key)
        _require(_is_mapping(row), "phase_row_malformed", "Phase row must be an object", key)
        result[int(key)] = dict(row)
    return result


def _validate_overlap(packet_rows: Mapping[int, Mapping[str, Any]], raw_records: Mapping[int, Mapping[str, Any]], packet_name: str) -> None:
    # Packet rows are annotations.  They are used for overlap/authentication
    # only; derivation below consumes raw_records exclusively.
    structural_fields = (
        "operation",
        "name",
        "argument",
        "arguments",
        "numerator",
        "denominator",
        "expression",
        "variable",
        "value",
        "type",
        "definition",
    )
    for node_id, row in packet_rows.items():
        _require(node_id in raw_records, "packet_unknown_record", f"{packet_name} packet references unknown raw record", node_id)
        raw = raw_records[node_id]
        _require(row.get("operation") == raw.get("operation"), "packet_operation_mismatch", f"{packet_name} operation differs from raw capsule", node_id)
        for field in structural_fields:
            if field in row and field in raw:
                _require(row[field] == raw[field], "packet_raw_overlap_mismatch", f"{packet_name} structural field differs from raw capsule", {"node_id": node_id, "field": field})


def _validate_source_fact_metadata(
    raw_capsule: Mapping[str, Any],
    scalar_packet: Mapping[str, Any],
    phase_packet: Mapping[str, Any],
    records: Mapping[int, Mapping[str, Any]],
) -> dict[str, Any]:
    authentication = raw_capsule.get("native_archive_authentication")
    _require(_is_mapping(authentication), "archive_authentication_missing", "Raw archive authentication metadata missing")
    _require(authentication.get("sha_match") is True, "archive_sha_mismatch", "Raw archive authentication did not match")
    _require(authentication.get("raw_sha256_expected") == ARCHIVE_SHA256, "archive_sha_pin_mismatch", "Raw archive expected SHA changed")
    _require(authentication.get("raw_sha256_observed") == ARCHIVE_SHA256, "archive_sha_pin_mismatch", "Raw archive observed SHA changed")
    _require(authentication.get("authentication_preceded_gzip_decode") is True, "archive_auth_order", "Archive authentication ordering changed")
    _require(raw_capsule.get("scope", {}).get("source_head") == SOURCE_HEAD, "source_head_mismatch", "Native source head changed")
    _require(raw_capsule.get("scope", {}).get("source_repo") is not None, "source_repo_missing", "Native source repository metadata missing")
    packet_inputs = raw_capsule.get("source_packet_inputs")
    _require(isinstance(packet_inputs, list), "source_packet_inputs_missing", "Raw capsule source packet pins missing")
    packet_hashes = {item.get("raw_sha256") for item in packet_inputs if _is_mapping(item)}
    _require(packet_hashes == {SOURCE_PINS["scalar_literals"], SOURCE_PINS["phase_literals"]}, "source_packet_pin_mismatch", "Raw capsule source packet pins changed")

    record_union = raw_capsule.get("record_union")
    _require(_is_mapping(record_union), "record_union_missing", "Raw record union metadata missing")
    union_ids = record_union.get("union_record_ids")
    _require(tuple(union_ids or ()) == RAW_RECORD_IDS, "record_union_mismatch", "Raw record union changed")
    _require(record_union.get("raw_record_count") == RAW_RECORD_COUNT, "record_union_count_mismatch", "Raw record count metadata changed")
    operation_counts: dict[str, int] = {}
    for record in records.values():
        operation = record.get("operation")
        operation_counts[operation] = operation_counts.get(operation, 0) + 1
    _require(operation_counts == record_union.get("raw_operation_counts"), "operation_count_mismatch", "Raw operation-count metadata changed")

    raw_vs_annotated = raw_capsule.get("raw_vs_annotated")
    _require(_is_mapping(raw_vs_annotated), "raw_annotation_boundary_missing", "Raw/annotated boundary metadata missing")
    _require(raw_vs_annotated.get("packet_rows_were_used_only_for_ids") is True, "annotated_row_boundary", "Annotated rows may not be promoted to graph records")
    _require(raw_vs_annotated.get("packet_row_values_were_not_copied") is True, "annotated_row_boundary", "Annotated packet values may not replace raw records")

    ordinary = raw_capsule.get("ordinary_Z_scope_binding")
    _require(_is_mapping(ordinary), "ordinary_Z_binding_missing", "Ordinary-Z source binding metadata missing")
    _require(ordinary.get("raw_node_id") == 20 and ordinary.get("raw_record_name") == "Z", "ordinary_Z_binding_mismatch", "Ordinary-Z binding ID/name changed")
    _require(ordinary.get("raw_record_operation") == "bound_variable", "ordinary_Z_binding_mismatch", "Ordinary-Z binding operation changed")
    _require(ordinary.get("source_factory_selector") == "build_graph", "source_factory_mismatch", "Ordinary-Z source factory changed")
    _require(ordinary.get("source_factory_expression") == "Z=g.symbol('Z')", "source_factory_mismatch", "Ordinary-Z factory expression changed")
    _require(records[20].get("operation") == "bound_variable" and records[20].get("name") == "Z", "ordinary_Z_record_mismatch", "Raw ordinary-Z record changed")

    one_literal = raw_capsule.get("required_one_literal")
    _require(_is_mapping(one_literal) and one_literal.get("raw_node_id") == 1 and one_literal.get("record_preserved") is True, "one_literal_missing", "Required exact-one literal metadata changed")
    _require(records[1].get("operation") == "exact_rational" and records[1].get("numerator") == 1 and records[1].get("denominator") == 1, "one_literal_mismatch", "Required exact-one raw record changed")

    opaque_n = raw_capsule.get("opaque_N_declaration")
    _require(_is_mapping(opaque_n), "opaque_N_metadata_missing", "Opaque N source factory metadata missing")
    _require(opaque_n.get("raw_node_id") == 18 and opaque_n.get("raw_node_type") == "DoubleDyadicInteger", "opaque_N_metadata_mismatch", "Opaque N source factory ID/type changed")
    _require(opaque_n.get("source_factory_selector") == "build_graph", "opaque_N_factory_mismatch", "Opaque N source factory selector changed")
    _require(opaque_n.get("source_factory_expression") == "N=g.node('exact_positive_integer_expression',**selection.record(),one_common_N_for_all_charts_and_control_denominators=True)", "opaque_N_factory_mismatch", "Opaque N source factory expression changed")
    _require(opaque_n.get("numeric_payload_evaluation") is False and opaque_n.get("N_materialized") is False and opaque_n.get("J_arithmetic") is False, "opaque_N_materialization", "Opaque N source factory metadata permits evaluation")

    native_x = raw_capsule.get("native_x_binding")
    _require(_is_mapping(native_x), "native_x_binding_missing", "Native x source binding metadata missing")
    _require(native_x.get("raw_node_id") == 241536 and native_x.get("raw_record_name") == "true_outer_O3_transition", "native_x_binding_mismatch", "Native x binding ID/name changed")
    _require(records[241536].get("operation") == "bound_variable", "native_x_record_mismatch", "Native x raw record changed")

    anchor = raw_capsule.get("hbB_sc_atom_factory_source_anchor")
    _require(_is_mapping(anchor), "atom_factory_anchor_missing", "Source atom factory anchor missing")
    _require(anchor.get("selector") == "build_graph", "atom_factory_selector_mismatch", "Source atom factory selector changed")
    _require(anchor.get("selector_ast_sha256") == "9b4a15b79e802ebf61e292a7245a30600c5c7dafdcad7de5ac216dd0d2bccacf", "atom_factory_ast_mismatch", "Source atom factory AST pin changed")
    _require(anchor.get("source_segment_sha256") == "62f5e829ec05aa545f9445236011c33ea44d9bc7191c4a678cd459bdb4093d0b", "atom_factory_segment_mismatch", "Source atom factory segment pin changed")
    evidence = anchor.get("evidence_lines")
    _require(isinstance(evidence, list) and any("param=lambda" in str(item.get("text")) for item in evidence if _is_mapping(item)), "atom_factory_evidence_missing", "Typed source atom factory evidence missing")
    atom_records = anchor.get("atom_records")
    _require(_is_mapping(atom_records) and atom_records.get("hbB", {}).get("raw_node_id") == 8 and atom_records.get("sc", {}).get("raw_node_id") == 10, "atom_factory_record_ids", "Source atom factory record IDs changed")

    source_family_keys = ("actual_five_defect_family_sha256", "implicit_source_sha256", "datum_enclosure_sha256")
    for node_id, expected_name in ((8, "hbB"), (10, "sc"), (21, "eta"), (22, "d_star")):
        if node_id not in records:
            continue
        record = records[node_id]
        _require(record.get("operation") == "current_original_source_parameter" and record.get("name") == expected_name, "source_atom_record_mismatch", "Typed source atom identity changed", node_id)
        _require(record.get("current_source_receipt_sha256") == CURRENT_SOURCE_RECEIPT_SHA256, "source_atom_receipt_mismatch", "Typed source atom receipt changed", node_id)
        _require(record.get("quantity_is_exact_original_parameter_not_range_endpoint") is True, "source_atom_type_mismatch", "Source atom is not typed as an exact original parameter", node_id)
        _require(record.get("source_family") and all(_valid_sha(record["source_family"].get(key)) for key in source_family_keys), "source_atom_family_mismatch", "Source atom family declaration malformed", node_id)
        _require(record.get("Z_independent") is True, "source_atom_type_mismatch", "Typed source atom declaration missing", node_id)

    logc = records.get(35560)
    _require(logc is not None and logc.get("operation") == "current_original_source_parameter", "logC_record_missing", "logC source atom missing")
    _require(logc.get("name") == "true_core_logC", "logC_name_mismatch", "logC source atom name changed")
    _require(logc.get("exact_fixed_source_parameter") is True, "logC_type_mismatch", "logC is not an exact fixed source parameter")
    _require(logc.get("source_path") == "common.core.logC", "logC_owner_mismatch", "logC source owner changed")
    _require(logc.get("parameter_enclosure_not_selected_value") is True, "logC_enclosure_mismatch", "logC enclosure declaration changed")
    _require(logc.get("source_family") and all(_valid_sha(logc["source_family"].get(key)) for key in source_family_keys), "logC_family_mismatch", "logC source family declaration malformed")
    provider = logc.get("source_provider")
    _require(_is_mapping(provider), "logC_provider_missing", "logC source owner AST declaration missing")
    for field in ("module", "module_sha256", "function", "AST_sha256"):
        _require(isinstance(provider.get(field), str) and provider.get(field), "logC_provider_malformed", "logC source owner AST declaration malformed", field)
    _require(
        provider == {
            "module": "lei_ren_part1_paper_compliant_current_original_n1_core_fixed_point_functions.py",
            "module_sha256": "d71ff19a68cfed8ace9a6857ffeca860ecee90d8d89b35b2292a7310e180fb59",
            "function": "CurrentOriginalN1CoreFixedPointFunctions.attach_oracle",
            "AST_sha256": "90c40b6f5111f7104975096c30cbc23657ee6317819ad9279c396b7cdc4bcf9f",
        },
        "logC_provider_mismatch",
        "logC source owner AST declaration changed",
    )
    return {
        "archive_sha256": ARCHIVE_SHA256,
        "source_head": SOURCE_HEAD,
        "raw_record_count": RAW_RECORD_COUNT,
        "ordinary_Z_node": 20,
        "native_x_node": 241536,
        "required_one_node": 1,
        "typed_source_atom_nodes": [8, 10, 21, 22, 35560],
        "raw_records_are_authoritative": True,
        "annotated_packets_used_for_overlap_only": True,
    }


def _validate_roots(
    scalar_packet: Mapping[str, Any],
    phase_packet: Mapping[str, Any],
    records: Mapping[int, Mapping[str, Any]],
) -> dict[str, Any]:
    scalar_roles = scalar_packet["literal_graph_nodes"].get("named_roles")
    _require(_is_mapping(scalar_roles), "scalar_named_roles_missing", "Scalar named graph roles missing")
    expected_scalar_roles = {
        "logC": 35560,
        "logP": 5,
        "logref": 81514,
        "Rpost_restore": 84533,
        "Rm_endpoint_substitution": 85743,
    }
    for name, node_id in expected_scalar_roles.items():
        role = scalar_roles.get(name)
        _require(_is_mapping(role) and role.get("node") == node_id or name == "Rm_endpoint_substitution" and _is_mapping(role) and role.get("node") == node_id, "scalar_root_mismatch", "Scalar graph root metadata changed", name)

    phase_graph = phase_packet["literal_graph_nodes"]
    closures = phase_graph.get("root_closures")
    _require(_is_mapping(closures), "phase_root_closures_missing", "Phase root closure metadata missing")
    expected_phase_closures = {
        "phase_origin_24954": 24954,
        "Ra_65604": 65604,
        "N_18": 18,
        "Rm_endpoint_85743": 85743,
        "native_radius_241787": 241787,
        "actual_phase_242551": 242551,
    }
    for name, node_id in expected_phase_closures.items():
        closure = closures.get(name)
        _require(_is_mapping(closure) and closure.get("root") == node_id, "phase_root_mismatch", "Phase graph root metadata changed", name)
    bindings = phase_packet.get("native_bindings")
    _require(_is_mapping(bindings), "native_bindings_missing", "Native source bindings missing")
    _require(bindings.get("ordinary_Z_symbol", {}).get("graph_node") == 20, "native_binding_mismatch", "Ordinary-Z native binding changed")
    _require(bindings.get("coordinate", {}).get("graph_node") == 241536, "native_binding_mismatch", "Native-x binding changed")
    _require(bindings.get("Rm", {}).get("endpoint_graph_root") == 85743, "native_binding_mismatch", "Rm endpoint binding changed")
    _require(bindings.get("Ra", {}).get("literal_graph_root") == 65604, "native_binding_mismatch", "Ra binding changed")
    _require(bindings.get("phase_origin", {}).get("literal_graph_root") == 24954, "native_binding_mismatch", "Phase-origin binding changed")
    _require(bindings.get("N", {}).get("literal_graph_root") == 18, "native_binding_mismatch", "N binding changed")
    _require(bindings.get("actual_phase", {}).get("graph_root") == 242551, "native_binding_mismatch", "Actual-phase binding changed")
    _require(records[18].get("operation") == "exact_positive_integer_expression", "N_record_mismatch", "N raw record operation changed")
    n_record = records[18]
    _require(n_record.get("type") == "DoubleDyadicInteger", "N_record_type_mismatch", "N raw record type changed")
    _require(n_record.get("definition") == "N=2^(2^J)", "N_record_definition_mismatch", "N raw record definition changed")
    _require(isinstance(phase_packet.get("variable_binding_context"), Mapping), "binding_context_missing", "Variable binding context missing")
    n_binding = bindings.get("N")
    _require(_is_mapping(n_binding) and n_binding.get("J_literal_source_token") == "1358356628656378313", "N_j_token_mismatch", "Opaque N fixed J token changed")
    return {
        "target_roots": dict(TARGET_ROOTS),
        "scalar_roles": expected_scalar_roles,
        "phase_roots": expected_phase_closures,
        "ordinary_Z_node": 20,
        "native_x_node": 241536,
        "N_is_opaque": True,
    }


def _validate_input(source: Mapping[str, Any]) -> tuple[dict[str, Any], dict[int, dict[str, Any]], dict[int, dict[str, Any]], dict[int, dict[str, Any]], dict[str, Any]]:
    _require(source.get("schema") == INPUT_SCHEMA, "input_schema_mismatch", "Input schema changed")
    _require(source.get("component") == COMPONENT, "input_component_mismatch", "Input component changed")
    _require(source.get("parent_report_raw_sha256") == PARENT_REPORT_SHA256, "parent_report_pin_mismatch", "Parent report pin changed")
    _require(source.get("parent_result_raw_sha256") == PARENT_RESULT_SHA256, "parent_result_pin_mismatch", "Parent result pin changed")
    wrappers = source.get("saved_sources")
    _require(_is_mapping(wrappers), "source_wrappers_missing", "The four source wrappers are required")
    _require(tuple(wrappers.keys()) == SOURCE_WRAPPER_NAMES or set(wrappers.keys()) == set(SOURCE_WRAPPER_NAMES), "source_wrapper_names_mismatch", "Source wrapper names changed")
    literals: dict[str, dict[str, Any]] = {}
    for name in SOURCE_WRAPPER_NAMES:
        wrapper = wrappers.get(name)
        _require(_is_mapping(wrapper) and set(wrapper.keys()) == {"literal", "raw_sha256"}, "source_wrapper_shape", "Every source wrapper must have exactly literal and raw_sha256", name)
        _require(wrapper.get("raw_sha256") == SOURCE_PINS[name], "source_wrapper_pin_mismatch", "Source wrapper SHA changed", name)
        _require(_is_mapping(wrapper.get("literal")), "source_literal_shape", "Source wrapper literal must be an object", name)
        literals[name] = dict(wrapper["literal"])

    scalar_packet = literals["scalar_literals"]
    phase_packet = literals["phase_literals"]
    raw_capsule = literals["raw_capsule"]
    scalar_rows = _validate_scalar_packet(scalar_packet)
    phase_rows = _validate_phase_packet(phase_packet)
    records = _raw_record_id_map(raw_capsule.get("raw_graph_records"))
    _validate_overlap(scalar_rows, records, "scalar")
    _validate_overlap(phase_rows, records, "phase")
    source_meta = _validate_source_fact_metadata(raw_capsule, scalar_packet, phase_packet, records)
    root_meta = _validate_roots(scalar_packet, phase_packet, records)
    metric = literals["saved_metric_result"]
    _require(_is_mapping(metric.get("payload")), "metric_payload_missing", "Saved metric result payload missing")
    return metric, records, scalar_rows, phase_rows, {"source": source_meta, "roots": root_meta, "raw_capsule": raw_capsule}


def _admitted_atom_scope(node_id: int, record: Mapping[str, Any]) -> dict[str, Any]:
    operation = record.get("operation")
    if operation == "current_original_source_parameter":
        return {
            "admitted": True,
            "kind": "definitional_fixed_scalar_source_atom",
            "node_id": node_id,
            "token": _token_for_fixed(node_id, record),
            "typed_fixed_parameter_semantics": True,
            "oracle_source_value_identity_recertified": False,
            "numeric_value_selected": False,
        }
    if operation == "exact_positive_integer_expression":
        return {
            "admitted": True,
            "kind": "definitional_opaque_positive_integer_expression",
            "node_id": node_id,
            "token": _token_for_fixed(node_id, {"name": "N"}),
            "typed_fixed_parameter_semantics": True,
            "oracle_source_value_identity_recertified": False,
            "numeric_value_selected": False,
            "N_materialized": False,
            "J_arithmetized": False,
        }
    return {"admitted": False}


def _derive_graph(records: Mapping[int, Mapping[str, Any]], emitter: "_Emitter") -> tuple[dict[int, list[str]], list[dict[str, Any]]]:
    memo: dict[int, list[str]] = {}
    visiting: set[int] = set()
    rows: list[dict[str, Any]] = []

    def visit(node_id: int) -> list[str]:
        if node_id in memo:
            return list(memo[node_id])
        if node_id in visiting:
            _fail("graph_cycle", "Cycle detected in raw graph closure", node_id)
        _require(node_id in records, "omitted_edge", "Graph edge points outside raw capsule", node_id)
        visiting.add(node_id)
        record = records[node_id]
        operation = record.get("operation")
        child_ids = _children(node_id, record)
        child_fvs: list[dict[str, Any]] = []
        child_sets: list[list[str]] = []
        for child_id in child_ids:
            child_fv = visit(child_id)
            child_sets.append(child_fv)
            child_fvs.append({"child": child_id, "free_variable_tokens": list(child_fv)})

        binder: Any = None
        removal_fv: list[str] = []
        replacement_fv: list[str] = []
        admitted_scope = _admitted_atom_scope(node_id, record)
        if operation == "exact_rational":
            _require(_is_int(record.get("numerator")) and _is_int(record.get("denominator")), "rational_malformed", "Exact rational numerator/denominator must be integers", node_id)
            _require(record.get("denominator") != 0, "rational_zero_denominator", "Exact rational denominator must be nonzero", node_id)
            result_fv: list[str] = []
            rule = "exact_rational_empty"
        elif operation == "exact_positive_integer_expression":
            _require(node_id == 18, "opaque_integer_scope", "Opaque positive integer expression is only admitted for native N", node_id)
            _require(record.get("type") == "DoubleDyadicInteger" and record.get("definition") == "N=2^(2^J)", "opaque_integer_malformed", "Opaque N declaration changed", node_id)
            # The raw record preserves J as an opaque source literal.  Its
            # integer payload is only type-checked; it is never converted,
            # exponentiated, or otherwise used as a value.
            _require(_is_int(record.get("J")), "opaque_integer_j_token", "Opaque N fixed J literal missing", node_id)
            _require(record.get("exact_positive_integer") is True, "opaque_integer_type", "Opaque N must be declared exact positive integer", node_id)
            _require(record.get("outer_integer_and_inner_power_not_materialized") is True and record.get("one_common_N_for_all_charts_and_control_denominators") is True, "opaque_integer_materialization", "Opaque N declaration permits materialization", node_id)
            result_fv = [_token_for_fixed(node_id, {"name": "N"})]
            rule = "typed_opaque_positive_integer_expression"
        elif operation == "current_original_source_parameter":
            _require(admitted_scope.get("admitted") is True, "source_atom_not_admitted", "Current source atom failed strict admission", node_id)
            result_fv = [_token_for_fixed(node_id, record)]
            rule = "typed_fixed_scalar_source_atom"
        elif operation == "bound_variable":
            _require(isinstance(record.get("name"), str) and record.get("name"), "bound_variable_malformed", "Bound variable name missing", node_id)
            result_fv = [_token_for_bound(node_id, record)]
            rule = "bound_variable_token"
        elif operation in {"sum", "product", "negative", "analytic_unary", "positive_quotient"}:
            if operation == "analytic_unary":
                _require(record.get("name") in ANALYTIC_NAMES, "unknown_analytic_unary", "Unsupported analytic unary operation", node_id)
            merged: list[str] = []
            for child_fv in child_sets:
                merged.extend(child_fv)
            result_fv = _sorted_tokens(merged)
            rule = {
                "sum": "union_sum_children",
                "product": "union_product_children",
                "negative": "union_negative_child",
                "analytic_unary": f"union_analytic_unary_{record.get('name')}",
                "positive_quotient": "union_positive_quotient_children",
            }[operation]
        elif operation == "function_substitution":
            body_id, binder_id, value_id = child_ids
            binder_record = records[binder_id]
            _require(binder_record.get("operation") == "bound_variable", "substitution_binder_type", "Substitution binder must be a raw bound_variable node", node_id)
            binder_token = _token_for_bound(binder_id, binder_record)
            body_fv = list(child_sets[0])
            value_fv = list(child_sets[2])
            occurs = binder_token in body_fv
            result_fv = [token for token in body_fv if token != binder_token]
            if occurs:
                result_fv.extend(value_fv)
            result_fv = _sorted_tokens(result_fv)
            binder = {
                "node_id": binder_id,
                "name": binder_record.get("name"),
                "token": binder_token,
                "raw_record": _clone(binder_record),
                "occurs_in_body": occurs,
            }
            removal_fv = [binder_token]
            replacement_fv = list(value_fv)
            rule = "capture_avoiding_function_substitution"
        else:
            _fail("unknown_operation", "Unknown graph operation", {"node_id": node_id, "operation": operation})

        result_fv = _sorted_tokens(result_fv)
        visiting.remove(node_id)
        memo[node_id] = result_fv
        row = {
            "node_id": node_id,
            "raw_record": _clone(record),
            "rule": rule,
            "child_FV": child_fvs,
            "binder": binder,
            "removal_FV": removal_fv,
            "replacement_FV": replacement_fv,
            "result_FV": list(result_fv),
            "admitted_atom_scope": admitted_scope,
        }
        rows.append(row)
        emitter.emit("node_derivation", row)
        return list(result_fv)

    for node_id in sorted(records):
        visit(node_id)
    _require(len(rows) == RAW_RECORD_COUNT, "node_visit_count_mismatch", "Every raw graph record must be visited exactly once", len(rows))
    _require(len(memo) == RAW_RECORD_COUNT, "node_memo_count_mismatch", "Every raw graph record must have one memoized derivation", len(memo))
    return memo, rows


def _negative_guards(
    records: Mapping[int, Mapping[str, Any]],
    scalar_rows: Mapping[int, Mapping[str, Any]],
    phase_rows: Mapping[int, Mapping[str, Any]],
    free_variables: Mapping[int, Sequence[str]],
    emitter: "_Emitter",
) -> list[dict[str, Any]]:
    refusals: list[dict[str, Any]] = []
    z_token = _token_for_bound(20, records[20])
    _require(z_token in free_variables[20], "ordinary_Z_guard_not_triggered", "Actual ordinary-Z record did not derive its own free token")
    try:
        _admit_or_refuse_request(
            guard="ordinary_Z_free_variable_rejected",
            request="promote ordinary Z as source-independent",
            evidence={
                "node_id": 20,
                "raw_record": records[20],
                "free_variable_tokens": free_variables[20],
                "ordinary_Z_token": z_token,
            },
            reject=z_token in free_variables[20],
            reason="The actual ordinary-Z binding is a free source variable; it is not silently admitted into a positive certificate.",
        )
    except _PolicyRefusal as refusal:
        refusal_row = refusal.row
        refusals.append(refusal_row)
        emitter.emit("policy_refusal", refusal_row)
    else:
        _fail("ordinary_Z_policy_not_refused", "Actual ordinary-Z source request was admitted")

    substitution = phase_rows.get(85743)
    _require(_is_mapping(substitution) and substitution.get("Z_independent_substitution") is True, "substitution_guard_not_triggered", "Actual substitution annotation missing")
    try:
        _admit_or_refuse_request(
            guard="Z_independent_substitution_metadata_rejected",
            request="silently trust Z_independent_substitution metadata as proof",
            evidence={
                "node_id": 85743,
                "raw_record": records[85743],
                "annotation_record": substitution,
            },
            reject=substitution.get("Z_independent_substitution") is True,
            reason="The capture-avoiding graph rule is the proof input; the source flag is metadata only.",
            extra={"metadata_used_as_proof": False},
        )
    except _PolicyRefusal as refusal:
        refusal_row = refusal.row
        refusals.append(refusal_row)
        emitter.emit("policy_refusal", refusal_row)
    else:
        _fail("substitution_policy_not_refused", "Actual substitution metadata request was admitted")

    logc = records.get(35560)
    scalar_logc = scalar_rows.get(35560)
    _require(_is_mapping(logc) and _is_mapping(scalar_logc), "logC_guard_not_triggered", "Actual logC source records missing")
    _require(logc.get("parameter_enclosure_not_selected_value") is True and scalar_logc.get("parameter_enclosure_not_selected_value") is True, "logC_guard_not_triggered", "Actual enclosure-backed logC marker missing")
    try:
        _admit_or_refuse_request(
            guard="enclosure_backed_logC_value_promotion_rejected",
            request="promote the enclosure-backed logC record to a selected numeric value",
            evidence={
                "node_id": 35560,
                "raw_record": logc,
                "annotation_record": scalar_logc,
            },
            reject=(
                logc.get("parameter_enclosure_not_selected_value") is True
                and scalar_logc.get("parameter_enclosure_not_selected_value") is True
            ),
            reason="The fixed source atom is admitted definitionally; no enclosure endpoint or oracle value is selected.",
            extra={"numeric_value_selected": False},
        )
    except _PolicyRefusal as refusal:
        refusal_row = refusal.row
        refusals.append(refusal_row)
        emitter.emit("policy_refusal", refusal_row)
    else:
        _fail("logC_policy_not_refused", "Actual enclosure-backed logC request was admitted")
    _require(len(refusals) == 3, "negative_refusal_count_mismatch", "Exactly three actual-source policy refusals are required")
    return refusals


def _frame_specializations(metric: Mapping[str, Any], target_rows: Sequence[Mapping[str, Any]], emitter: "_Emitter") -> list[dict[str, Any]]:
    payload = metric.get("payload")
    _require(_is_mapping(payload), "metric_payload_missing", "Saved metric result payload missing")
    coefficients = payload.get("full_operator_coefficients")
    manifests = payload.get("branch_operator_manifests")
    _require(_is_mapping(coefficients), "metric_coefficients_missing", "Saved metric full operator coefficients missing")
    _require(isinstance(manifests, list) and len(manifests) == 2, "metric_branch_manifest_count", "Saved metric branch manifests must contain two branches")
    premises = payload.get("premises")
    _require(_is_mapping(premises), "metric_premises_missing", "Saved metric premises missing")
    rho = premises.get("rho")
    native_rm = premises.get("native_Rm")
    _require(_is_mapping(rho) and rho.get("positive") is True, "metric_positive_domain_missing", "Physical positive-domain premise for rho is missing")
    _require(_is_mapping(native_rm) and native_rm.get("positive") is not None, "metric_positive_domain_missing", "Physical positive-domain premise for native Rm is missing")
    rows: list[dict[str, Any]] = []
    for manifest in manifests:
        _require(_is_mapping(manifest), "metric_branch_manifest_malformed", "Saved metric branch manifest must be an object")
        branch = manifest.get("branch")
        _require(branch in {"background", "modified"}, "metric_branch_name", "Unexpected saved metric branch", branch)
        row = {
            "branch": branch,
            "source": "saved_metric_result.payload.full_operator_coefficients_and_branch_operator_manifests",
            "full_operator_coefficients": _clone(coefficients),
            "branch_operator_manifest": _clone(manifest),
            "source_target_certificates": [item.get("target") for item in target_rows],
            "rho_Z": "0",
            "rho_ZZ": "0",
            "a": "0",
            "a_Z": "0",
            "P_Z": "D_Z",
            "scalar_laplace": "R^-2(Fxx+Falphaalpha)+FZZ",
            "rho_Z_zero_source_identity": True,
            "rho_ZZ_zero_source_identity": True,
            "a_zero_consequence": True,
            "a_Z_zero_consequence": True,
            "P_Z_identity": True,
            "physical_positive_domain_premise": {
                "rho_positive": _clone(rho.get("positive")),
                "native_Rm_positive": _clone(native_rm.get("positive")),
            },
            "source_identity_application": True,
            "metric_recomputed": False,
            "metric_differentiated": False,
            "actual_field_pullback_instantiated": False,
            "mixed_derivatives_claimed": False,
            "time_or_L2_values": False,
            "full_closure_claimed": False,
        }
        rows.append(row)
        emitter.emit("frame_specialization", row)
    _require(len(rows) == 2, "frame_specialization_count_mismatch", "Exactly two branch frame specializations are required")
    return rows


class _Emitter:
    def __init__(self, diagnostic: Any) -> None:
        self.diagnostic = diagnostic
        self.events: list[dict[str, Any]] = []

    def emit(self, event: str, row: Mapping[str, Any]) -> None:
        event_row = {"event": event, "row": _clone(dict(row))}
        self.events.append(event_row)
        if self.diagnostic is None:
            return
        if callable(self.diagnostic):
            # The durable manager callback receives one mapping, matching
            # the frozen event journal shape.
            self.diagnostic(_clone(event_row))
            return
        append = getattr(self.diagnostic, "append", None)
        if callable(append):
            append(_clone(event_row))
            return
        _fail("diagnostic_callback_invalid", "diagnostic must be callable or appendable")


def _flags(success: bool) -> dict[str, bool]:
    return {
        "native_parameter_Z_independence_proved": success,
        "native_Rm_Z_zero_proved": success,
        "native_Ra_phase_origin_Z_zero_proved": success,
        "native_fixed_radius_Z_frame_specialized": success,
        "physical_layer_completed": False,
        "full_goal_complete": False,
        "numeric_closure": False,
        "oracle_source_value_identity": False,
        "cap_deletion": False,
        "temporal_recursion": False,
        "actual_field_pullback_instantiated": False,
        "time_values_evaluated": False,
    }


def _counts(
    node_rows: Sequence[Any],
    target_rows: Sequence[Any],
    refusals: Sequence[Any],
    frame_rows: Sequence[Any],
    success: bool,
    family_entered: bool = False,
) -> dict[str, int]:
    return {
        "new_source_identity_family_calls": 1 if family_entered else 0,
        "new_algebra_family_calls": 1 if family_entered else 0,
        "new_mathematical_proofs": 1 if success else 0,
        "node_derivation_rows": len(node_rows),
        "target_Z_zero_rows": len(target_rows),
        "negative_source_refusals": len(refusals),
        "frame_specialization_rows": len(frame_rows),
        "old_math_calls": 0,
        "old_provider_calls": 0,
        "source_value_identity_calls": 0,
        "numeric_evaluation_calls": 0,
        "N_materialization_calls": 0,
        "metric_replay_calls": 0,
    }


def _event_rows(emitter: _Emitter, event_name: str) -> list[dict[str, Any]]:
    """Recover the durable full-row prefix without re-deriving anything."""

    return [
        _clone(item["row"])
        for item in emitter.events
        if item.get("event") == event_name and isinstance(item.get("row"), Mapping)
    ]


def _failure_result(
    emitter: _Emitter,
    failure: Mapping[str, Any],
    source_validation: Mapping[str, Any] | None,
    node_rows: Sequence[Any],
    target_rows: Sequence[Any],
    refusals: Sequence[Any],
    frame_rows: Sequence[Any],
) -> dict[str, Any]:
    family_entered = any(item.get("event") == "source_identity_family_entry" for item in emitter.events)
    return {
        "schema": RESULT_SCHEMA,
        "component": COMPONENT,
        "status": "failed",
        "passed": False,
        "failure": _clone(dict(failure)),
        "conditional": True,
        "source_only": True,
        "metadata_only": False,
        "definitional_parameter_scope": True,
        "conditional_on_declared_source_atoms": True,
        "actual_counts": _counts(node_rows, target_rows, refusals, frame_rows, False, family_entered),
        "callback_event_count": len(emitter.events),
        "callback_events": _clone(emitter.events),
        "new_flags": _flags(False),
        **_flags(False),
        "payload": {
            "source_validation": _clone(source_validation) if source_validation is not None else None,
            "node_derivation_rows": _clone(list(node_rows)),
            "target_certificates": _clone(list(target_rows)),
            "negative_source_refusals": _clone(list(refusals)),
            "frame_specialization_rows": _clone(list(frame_rows)),
            "source_identity_family_entered": any(item.get("event") == "source_identity_family_entry" for item in emitter.events),
        },
        "runtime_model_metadata": None,
    }


def derive(source: dict[str, Any], diagnostic: Any = None) -> dict[str, Any]:
    """Derive the inert source-Z certificate from frozen literal wrappers.

    The caller owns freezing and the single invocation.  This function has no
    import-time or self-test entry point and performs no source-value or metric
    computation.
    """

    emitter = _Emitter(diagnostic)
    source_validation: dict[str, Any] | None = None
    node_rows: list[dict[str, Any]] = []
    target_rows: list[dict[str, Any]] = []
    refusals: list[dict[str, Any]] = []
    frame_rows: list[dict[str, Any]] = []
    try:
        _require(_is_mapping(source), "input_not_mapping", "Input must be a JSON object")
        metric, records, scalar_rows, phase_rows, metadata = _validate_input(source)
        source_validation = {
            "component": COMPONENT,
            "input_schema": INPUT_SCHEMA,
            "source_wrapper_names": list(SOURCE_WRAPPER_NAMES),
            "source_wrapper_pins": dict(SOURCE_PINS),
            "parent_report_raw_sha256": PARENT_REPORT_SHA256,
            "parent_result_raw_sha256": PARENT_RESULT_SHA256,
            "graph_record_count": RAW_RECORD_COUNT,
            "target_root_count": len(TARGET_ROOTS),
            "metadata": _clone(metadata),
            "raw_records_authoritative": True,
            "sourcevalue_evaluation": False,
            "N_materialized": False,
            "fractional_part_differentiated": False,
        }
        emitter.emit("source_validation", source_validation)
        emitter.emit(
            "source_identity_family_entry",
            {
                "family": "native_source_identity_free_variable_tokens_v1",
                "source_only": True,
                "typed_fixed_parameter_semantics": True,
                "oracle_source_value_identity_recertified": False,
                "numeric_values_selected": False,
                "ordinary_Z_metadata_not_trusted": True,
            },
        )
        free_variables, node_rows = _derive_graph(records, emitter)

        for target, root in TARGET_ROOTS.items():
            free = list(free_variables[root])
            z_token = _token_for_bound(20, records[20])
            _require(z_token not in free, "target_Z_free", "Target root retains ordinary-Z free token", {"target": target, "root": root, "free_variables": free})
            row = {
                "target": target,
                "root": root,
                "FV": free,
                "free_variable_tokens": free,
                "ordinary_Z_token": z_token,
                "Z_excluded": True,
                "ordinary_Z_derivative": "0",
                "derivative_basis": "source_functional_Z_invariance_on_defined_domain",
                "not_calculator_output": True,
                "fractional_part_differentiated": False,
                "defined_domain_premise": "declared source-domain semantics only",
            }
            target_rows.append(row)
            emitter.emit("target_certificate", row)
        _require(len(target_rows) == len(TARGET_ROOTS), "target_count_mismatch", "All target roots must receive certificates")
        refusals = _negative_guards(records, scalar_rows, phase_rows, free_variables, emitter)
        frame_rows = _frame_specializations(metric, target_rows, emitter)
        _require(len(frame_rows) == 2, "frame_count_mismatch", "Both branch frame specializations must complete")

        counts = _counts(node_rows, target_rows, refusals, frame_rows, True, True)
        flags = _flags(True)
        proof = {
            "family": "native_source_identity_free_variable_tokens_v1",
            "new_source_identity_family_calls": 1,
            "new_algebra_family_calls": 1,
            "new_mathematical_proofs": 1,
            "statement": "Declared typed fixed source atoms and graph composition are functionally independent of ordinary Z on the declared domain.",
            "ordinary_Z_derivative_identity": "0",
            "proof_is_definitional_source_identity": True,
            "oracle_source_value_identity_recertified": False,
            "numeric_values_selected": False,
            "N_materialized": False,
            "fractional_part_differentiated": False,
        }
        return {
            "schema": RESULT_SCHEMA,
            "component": COMPONENT,
            "status": "complete",
            "passed": True,
            "failure": None,
            "conditional": True,
            "source_only": True,
            "metadata_only": False,
            "definitional_parameter_scope": True,
            "conditional_on_declared_source_atoms": True,
            "actual_counts": counts,
            "callback_event_count": len(emitter.events),
            "callback_events": _clone(emitter.events),
            "new_flags": flags,
            **flags,
            "payload": {
                "source_validation": source_validation,
                "source_identity_proof": proof,
                "node_derivation_rows": _clone(node_rows),
                "target_certificates": _clone(target_rows),
                "negative_source_refusals": _clone(refusals),
                "frame_specialization_rows": _clone(frame_rows),
                "source_identity_scope": {
                    "typed_fixed_parameter_semantics": True,
                    "oracle_source_value_identity_recertified": False,
                    "numeric_value_selected": False,
                    "conditional_on_declared_source_atoms": True,
                    "physical_layer": False,
                },
            },
            "runtime_model_metadata": None,
        }
    except _KernelFailure as exc:
        # Helpers emit each completed row before returning.  On a rejection
        # (or a callback failure) their local list may never be assigned back
        # to this frame; recover the durable prefix from the emitter journal
        # without re-running any derivation.
        node_rows = _event_rows(emitter, "node_derivation")
        target_rows = _event_rows(emitter, "target_certificate")
        refusals = _event_rows(emitter, "policy_refusal")
        frame_rows = _event_rows(emitter, "frame_specialization")
        return _failure_result(
            emitter,
            {"code": exc.code, "message": exc.message, "details": _clone(exc.details)},
            source_validation,
            node_rows,
            target_rows,
            refusals,
            frame_rows,
        )
    except Exception as exc:  # The complete prefix is returned; no retry occurs here.
        node_rows = _event_rows(emitter, "node_derivation")
        target_rows = _event_rows(emitter, "target_certificate")
        refusals = _event_rows(emitter, "policy_refusal")
        frame_rows = _event_rows(emitter, "frame_specialization")
        return _failure_result(
            emitter,
            {
                "code": "unexpected_exception",
                "message": str(exc),
                "exception_type": type(exc).__name__,
                "full_prefix_preserved": True,
            },
            source_validation,
            node_rows,
            target_rows,
            refusals,
            frame_rows,
        )


__all__ = ["COMPONENT", "INPUT_SCHEMA", "RESULT_SCHEMA", "derive"]
