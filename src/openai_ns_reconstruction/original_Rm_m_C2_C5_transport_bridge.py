"""Source-only bridge between the finite-N C2 and prescribed C5 m histories.

The C2 archive stores the N-scaled density and its native partitioned
transport, while the C5 archive stores the corresponding rate-one loop
transport.  This module authenticates those two saved source projections and
checks their incoming and non-reset endpoint recurrence through restoration.
It does not run either producer, evaluate an operator, or infer the density
identity from graph handles.  Admission of the relation is deliberately
conditional on an explicit density-provider certificate supplied in
``source_records['same_density_premise']``.

``source_records`` is a mapping from each name in :data:`RECORD_HASHES` to a
record with the exact ``experiments/root_st073/<name>`` path, the raw
``sha256``, and UTF-8 ``source_text``.  For a gzip record the runner must also
provide ``raw_bytes_base64`` (or ``raw_bytes`` / ``gzip_bytes``) and the raw
hash is checked before decompression.  The optional density premise is a
mapping, not a file record.  To admit it, the mapping must contain an exact
``certificate_sha256`` equal to the canonical digest of the mapping with that
field removed, the expected source graph digest and count, all
:data:`RECORD_HASHES`, and ``per_chart_bindings`` for each first-nine chart.
Each binding names m, orders 0, 1, 2, relation ``C2=N*C5``, densities
``B`` and ``B_over_Pstar/N``, and the exact C2/C5 handles returned by this
module.  A bare ``passed`` or ``admitted`` flag is never sufficient.
"""

from __future__ import annotations

import ast
import base64
import binascii
import copy
from collections.abc import Mapping, Sequence
import gzip
import hashlib
import json
from typing import Any


SCHEMA = "original_Rm_m_C2_C5_transport_bridge_v1"
PREFIX = "lei_ren_part1_paper_compliant_current_original_"
EXPECTED_NODES_DIGEST = "a2dcc782ba3033ff6ca5f54a66632ac6ef007bb0c5b9960f10567d44ef2059a3"
EXPECTED_NODE_COUNT = 205300
PARENT_REPORT_SHA256 = "8903b59267bae01d5182bf2dcc84ed3c02c44eb8c7d3dd9d7b14af67ac092746"
POSITIVE_REPORT_SHA256 = "4c3b730e57966c0c4a8c2481564a6bf6ec19b1c1136ab01f7e7e2ade1246f675"

C2_PHASE_PY = PREFIX + "C2_phase_primitives.py"
C2_PHASE_CHECK = PREFIX + "C2_phase_primitives_check.json"
C2_PHASE_GZ = PREFIX + "C2_phase_primitives.json.gz"
C2_DENSITY_PY = PREFIX + "C2_density_transport.py"
C2_DENSITY_CHECK = PREFIX + "C2_density_transport_check.json"
C2_DENSITY_GZ = PREFIX + "C2_density_transport.json.gz"
C5_PY = PREFIX + "prescribed_loop_C5_functions.py"
C5_CHECK = PREFIX + "prescribed_loop_C5_functions_check.json"
C5_GZ = PREFIX + "prescribed_loop_C5_functions.json.gz"

# These are the accepted parent input pins.  The check.py files are present in
# the workspace but are intentionally not listed: the accepted parent report
# does not pin them as inputs, and the integration runner rejects new pins.
RECORD_HASHES = {
    C2_PHASE_PY: "84f5da3eb64096f82df88dd828b22b09bce18192a03925d91cfee0b193902d2c",
    C2_PHASE_CHECK: "76ed1af50ec799c49081ab67881b189b42e0d3137bcca770a67766d4111dff8b",
    C2_PHASE_GZ: "80e8f6df15edd94be6cec6aaf24528ab1eb1f06c7f7c5b61da47479e733a6bb0",
    C2_DENSITY_PY: "8bfa32c87c73c5bedf3056c7171707ae934dafd950bd0df874dbd3b3450379d1",
    C2_DENSITY_CHECK: "c4d2c62dd67e9dd7e171fabbdb9ba18308e420c5e0c98d01147c708840b8f583",
    C2_DENSITY_GZ: "0f347824bea798ea692eca7967d7acadba49fe4336989bea1a5c86c922260877",
    C5_PY: "69d6391120000488472b7961b4bff1901bff872def95d183c4ee056b0948e17d",
    C5_CHECK: "2c77e80df04fae615f8ea909fc7b606bfccfa27659cb67744ae8d84020b1a9e7",
    C5_GZ: "8419fdc338dc746099c7cad119ea61a8e928652319d5caf5314f8a730c533caa",
}

CHARTS = (
    "first_micro",
    "second_micro",
    "frozen_macro",
    "first_switch",
    "second_switch",
    "post_power",
    "long_reshape",
    "reference",
    "restoration",
)
POSTRESTORE = "postrestore"
ORDERS = (0, 1, 2)
M_COMPONENT = "m"


def canonical_digest(value: Any) -> str:
    raw = json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
        allow_nan=False,
    ).encode("utf-8")
    return hashlib.sha256(raw).hexdigest()


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _row(nodes: Sequence[Any], index: int) -> Mapping[str, Any]:
    _require(type(index) is int and 0 <= index < len(nodes), f"graph handle out of range: {index}")
    value = nodes[index]
    _require(isinstance(value, Mapping), f"graph handle is not a mapping: {index}")
    return value


def _bytes_from_record(record: Mapping[str, Any]) -> bytes | None:
    for key in ("raw_bytes", "source_bytes", "gzip_bytes"):
        value = record.get(key)
        if isinstance(value, bytes):
            return value
        if isinstance(value, bytearray):
            return bytes(value)
    encoded = record.get("raw_bytes_base64")
    if isinstance(encoded, str):
        try:
            return base64.b64decode(encoded, validate=True)
        except (ValueError, binascii.Error):
            return None
    return None


def _record(source_records: Mapping[str, Any], name: str, gates: list[dict[str, Any]]) -> dict[str, Any] | None:
    expected = RECORD_HASHES[name]
    record = source_records.get(name) if isinstance(source_records, Mapping) else None
    if not isinstance(record, Mapping):
        gates.append({"gate": "source_record_missing", "record": name})
        return None
    expected_path = "experiments/root_st073/" + name
    if record.get("path") != expected_path:
        gates.append({"gate": "source_record_path_changed", "record": name, "expected": expected_path})
        return None
    if record.get("sha256") != expected:
        gates.append({"gate": "source_record_hash_changed", "record": name, "expected_sha256": expected})
        return None

    is_gzip = name.endswith(".gz")
    raw: bytes | None = None
    if is_gzip:
        # A compressed input must be authenticated from the compressed bytes
        # before any decompression is accepted.
        raw = _bytes_from_record(record)
        if raw is None:
            gates.append({"gate": "gzip_raw_bytes_missing", "record": name})
            return None
        raw_sha = hashlib.sha256(raw).hexdigest()
        if raw_sha != expected:
            gates.append({"gate": "gzip_raw_hash_changed", "record": name, "expected_sha256": expected})
            return None
        try:
            decoded = gzip.decompress(raw)
        except (OSError, EOFError) as exc:
            gates.append({"gate": "gzip_decode_failed", "record": name, "error": str(exc)})
            return None
        supplied_decoded = record.get("source_text")
        try:
            text = decoded.decode("utf-8")
        except UnicodeDecodeError as exc:
            gates.append({"gate": "gzip_utf8_decode_failed", "record": name, "error": str(exc)})
            return None
        if isinstance(supplied_decoded, str) and supplied_decoded != text:
            gates.append({"gate": "gzip_decoded_text_changed", "record": name})
            return None
        decoded_sha = record.get("decoded_sha256")
        if decoded_sha is not None and decoded_sha != hashlib.sha256(decoded).hexdigest():
            gates.append({"gate": "gzip_decoded_hash_changed", "record": name})
            return None
    else:
        text = record.get("source_text")
        if not isinstance(text, str):
            gates.append({"gate": "source_text_missing", "record": name})
            return None
        if hashlib.sha256(text.encode("utf-8")).hexdigest() != expected:
            gates.append({"gate": "source_text_hash_changed", "record": name, "expected_sha256": expected})
            return None
        decoded = text.encode("utf-8")
    return {"name": name, "sha256": expected, "text": text, "decoded": decoded, "raw_sha256": hashlib.sha256(raw).hexdigest() if raw is not None else None}


def _definitions(text: str, name: str, wanted: tuple[str, ...], gates: list[dict[str, Any]]) -> dict[str, dict[str, int]]:
    try:
        tree = ast.parse(text, filename=name)
    except SyntaxError as exc:
        gates.append({"gate": "source_ast_invalid", "record": name, "error": str(exc)})
        return {}
    found: dict[str, dict[str, int]] = {}
    for item in ast.walk(tree):
        if isinstance(item, (ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)) and item.name in wanted:
            found[item.name] = {"line": item.lineno, "end_line": item.end_lineno or item.lineno}
    missing = [name for name in wanted if name not in found]
    if missing:
        gates.append({"gate": "source_declaration_missing", "record": name, "declarations": missing})
    return found


def _tokens(text: str, name: str, wanted: tuple[str, ...], gates: list[dict[str, Any]]) -> bool:
    missing = [token for token in wanted if token not in text]
    if missing:
        gates.append({"gate": "source_constructor_token_missing", "record": name, "tokens": missing})
        return False
    return True


def _json_record(record: dict[str, Any], name: str, gates: list[dict[str, Any]]) -> Mapping[str, Any]:
    try:
        value = json.loads(record["text"])
    except (KeyError, TypeError, json.JSONDecodeError) as exc:
        gates.append({"gate": "accepted_source_receipt_invalid", "record": name, "error": str(exc)})
        return {}
    if not isinstance(value, Mapping):
        gates.append({"gate": "accepted_source_receipt_not_mapping", "record": name})
        return {}
    return value


def _flag(value: Mapping[str, Any], path: tuple[str, ...]) -> Any:
    current: Any = value
    for key in path:
        current = current.get(key) if isinstance(current, Mapping) else None
    return current


def _receipt(record: dict[str, Any], name: str, required: tuple[tuple[str, ...], ...], gates: list[dict[str, Any]]) -> tuple[bool, Mapping[str, Any]]:
    value = _json_record(record, name, gates)
    ok = value.get("all_passed") is True
    if not ok:
        gates.append({"gate": "accepted_source_receipt_header_missing", "record": name})
    for path in required:
        if _flag(value, path) is not True:
            gates.append({"gate": "accepted_source_receipt_flag_missing", "record": name, "path": ".".join(path)})
            ok = False
    return ok, value


def _source_witnesses(source_records: Mapping[str, Any], gates: list[dict[str, Any]]) -> tuple[bool, dict[str, Any], dict[str, Any]]:
    records = {name: _record(source_records, name, gates) for name in RECORD_HASHES}
    if any(value is None for value in records.values()):
        return False, {"records": dict(RECORD_HASHES)}, {}
    assert all(isinstance(value, dict) for value in records.values())
    declarations = {
        C2_PHASE_PY: _definitions(records[C2_PHASE_PY]["text"], C2_PHASE_PY, ("C2Algebra", "primitive_C2", "source_bindings", "CurrentC2PhasePrimitives"), gates),
        C2_DENSITY_PY: _definitions(records[C2_DENSITY_PY]["text"], C2_DENSITY_PY, ("coefficient_C2", "build", "CurrentC2DensityTransport"), gates),
        C5_PY: _definitions(records[C5_PY]["text"], C5_PY, ("LoopAlgebra", "loop_functions", "finite_N_functions", "transport", "CurrentOriginalPrescribedLoopC5Functions"), gates),
    }
    token_ok = all((
        _tokens(records[C2_PHASE_PY]["text"], C2_PHASE_PY, ("source_bindings", "ordinary", "Z2", "inverse"), gates),
        _tokens(records[C2_DENSITY_PY]["text"], C2_DENSITY_PY, ("g.unary('exp'", "g.sub(ref(old['right_offset']),ref(old['offset']))", "g.mul(kernel,q,ref(old['Jacobian']))", "N_scaled_density=True", "ordinary_slow_Z_derivative_order=2", "original dlogR: Jacobian exactly once", "integration_endpoints_kernel_radius_and_global_phase_Z_independent=True"), gates),
        _tokens(records[C5_PY]["text"], C5_PY, ("g.unary('exp'", "g.mul(g.constant(-rate),jacobian,g.sub(coordinate,t))", "physical_Jacobian_applied_once=True", "zero_correction_inlet_only_at_original_rminus=True", "predecessor_memory_never_reset=True"), gates),
    ))
    phase_ok, phase_receipt = _receipt(records[C2_PHASE_CHECK], C2_PHASE_CHECK, (
        ("current_original_17_chart_source_Z2_rows_installed",),
        ("current_original_same_inverse_C2_phase_primitive_functions_installed",),
        ("ordinary_source_second_derivatives_use_actual_Taylor_rows",),
    ), gates)
    density_ok, density_receipt = _receipt(records[C2_DENSITY_CHECK], C2_DENSITY_CHECK, (
        ("current_original_17_chart_C2_signed_density_functions_installed",),
        ("current_original_17_chart_C2_own_rate_transport_and_five_target_functions_installed",),
        ("actual_five_target_ZZ_function_representations_installed",),
        ("current_C2_target_functions_are_enclosure_free_definitions",),
    ), gates)
    c5_ok, c5_receipt = _receipt(records[C5_CHECK], C5_CHECK, (
        ("checks", "actual_source", "passed"),
        ("checks", "transport", "passed"),
        ("checks", "transport", "all_five_original_rates_and_single_physical_measure_checked"),
        ("current_original_true_prescribed_ten_chart_signed_loop_C5_functions_installed",),
        ("current_original_true_prescribed_ten_chart_finite_N_density_C5_functions_installed",),
        ("current_original_true_prescribed_rminus_to_Rm_correction_history_C5_functions_installed",),
    ), gates)
    parsed = {}
    for name in (C2_PHASE_GZ, C2_DENSITY_GZ, C5_GZ):
        try:
            parsed[name] = json.loads(records[name]["text"])
        except json.JSONDecodeError as exc:
            gates.append({"gate": "source_archive_json_invalid", "record": name, "error": str(exc)})
    ok = token_ok and phase_ok and density_ok and c5_ok and all(bool(value) for value in declarations.values()) and len(parsed) == 3
    witness = {
        "records": dict(RECORD_HASHES),
        "raw_gzip_hashes_verified_before_decode": {name: records[name]["raw_sha256"] for name in (C2_PHASE_GZ, C2_DENSITY_GZ, C5_GZ)},
        "declarations": declarations,
        "accepted_receipts": {
            C2_PHASE_CHECK: {"sha256": RECORD_HASHES[C2_PHASE_CHECK], "passed": phase_ok},
            C2_DENSITY_CHECK: {"sha256": RECORD_HASHES[C2_DENSITY_CHECK], "passed": density_ok},
            C5_CHECK: {"sha256": RECORD_HASHES[C5_CHECK], "passed": c5_ok},
        },
        "parent_report_sha256": PARENT_REPORT_SHA256,
        "accepted_a_positivity_report_sha256": POSITIVE_REPORT_SHA256,
        "source_chain_structural": ok,
    }
    return ok, witness, parsed


def _copy_rows(nodes: Sequence[Any], reference: Sequence[Any], handles: set[int], label: str, gates: list[dict[str, Any]]) -> bool:
    ok = True
    for index in sorted(handles):
        try:
            actual = _row(nodes, index)
            expected = _row(reference, index)
        except ValueError as exc:
            gates.append({"gate": "saved_source_graph_handle_invalid", "label": label, "error": str(exc), "handle": index})
            ok = False
            continue
        if actual != expected:
            gates.append({"gate": "current_graph_source_row_changed", "label": label, "handle": index})
            ok = False
    return ok


def _as_int_list(value: Any, length: int, label: str, gates: list[dict[str, Any]]) -> list[int] | None:
    if not isinstance(value, list) or len(value) < length or any(type(item) is not int for item in value[:length]):
        gates.append({"gate": "source_handle_vector_changed", "label": label})
        return None
    return value[:length]


def _c2_witness(c2_archive: Mapping[str, Any], nodes: Sequence[Any], gates: list[dict[str, Any]]) -> tuple[bool, dict[str, Any]]:
    target = c2_archive.get("actual_C2_density_transport_and_target_functions")
    graph = c2_archive.get("exact_graph_nodes")
    if not isinstance(target, Mapping) or not isinstance(graph, list) or c2_archive.get("candidate_C2_transport_constructed") is not True:
        gates.append({"gate": "c2_archive_shape_changed"})
        return False, {}
    archive_flags_ok = True
    for flag in ("current_original_17_chart_C2_signed_density_functions_installed", "current_original_17_chart_C2_own_rate_transport_and_five_target_functions_installed"):
        if c2_archive.get(flag) is not True:
            gates.append({"gate": "c2_archive_acceptance_flag_missing", "flag": flag})
            archive_flags_ok = False
    windows = target.get("windows")
    if not isinstance(windows, list) or [w.get("chart") for w in windows[:len(CHARTS)] if isinstance(w, Mapping)] != list(CHARTS):
        gates.append({"gate": "c2_first_nine_chart_order_changed"})
        return False, {}
    rows: list[dict[str, Any]] = []
    handles: set[int] = set()
    ok = True
    previous_outgoing: list[int] | None = None
    previous_partitions: list[Mapping[str, Any]] | None = None
    partition_summary: dict[str, Any] = {}
    for position, chart in enumerate(CHARTS):
        window = windows[position]
        if not isinstance(window, Mapping):
            ok = False
            gates.append({"gate": "c2_window_not_mapping", "chart": chart})
            continue
        incoming = _as_int_list(window.get("incoming", {}).get("m"), 3, f"C2 {chart} incoming m", gates)
        outgoing = _as_int_list(window.get("outgoing", {}).get("m"), 3, f"C2 {chart} outgoing m", gates)
        contribution = _as_int_list(window.get("contributions", {}).get("m"), 3, f"C2 {chart} contribution m", gates)
        memory = window.get("original_own_rate_memory", {}).get("m")
        if not isinstance(memory, int):
            gates.append({"gate": "c2_memory_handle_changed", "chart": chart})
            ok = False
        if incoming is None or outgoing is None or contribution is None or not isinstance(memory, int):
            ok = False
            continue
        if window.get("predecessor_ZZ_memory_retained") is not True or window.get("quiet_local_zero_does_not_reset_incoming") is not True:
            gates.append({"gate": "c2_memory_or_zero_inlet_guard_changed", "chart": chart})
            ok = False
        if position == 0:
            first_rows = [_row(graph, item) for item in incoming]
            if any(row.get("operation") != "exact_rational" or row.get("numerator") != 0 or row.get("denominator") != 1 for row in first_rows):
                gates.append({"gate": "c2_zero_inlet_changed", "chart": chart})
                ok = False
        elif incoming != previous_outgoing:
            gates.append({"gate": "c2_incoming_recurrence_changed", "chart": chart})
            ok = False
        previous_outgoing = outgoing
        pieces = window.get("original_physical_measure_and_partitions")
        expected_count = 4 if chart in ("first_micro", "second_micro") else 2
        if not isinstance(pieces, list) or len(pieces) != expected_count:
            gates.append({"gate": "c2_native_partition_count_changed", "chart": chart, "expected": expected_count})
            ok = False
            pieces = []
        for left, right in zip(pieces, pieces[1:]):
            if not isinstance(left, Mapping) or not isinstance(right, Mapping) or left.get("right_offset") != right.get("left_offset") or left.get("exact_right") != right.get("exact_left"):
                gates.append({"gate": "c2_native_partition_endpoint_additivity_changed", "chart": chart})
                ok = False
        if previous_partitions is not None and pieces and previous_partitions:
            if previous_partitions[-1].get("right_offset") != pieces[0].get("left_offset"):
                gates.append({"gate": "c2_cross_chart_offset_additivity_changed", "chart": chart})
                ok = False
        previous_partitions = pieces
        partition_summary[chart] = copy.deepcopy(pieces)

        # Authenticate the source archive's m recurrence and the current graph
        # rows.  The actual outgoing root is sum(memory * incoming, integral).
        handles.update(incoming + outgoing + contribution + [memory])
        for order in ORDERS:
            try:
                out_row = _row(graph, outgoing[order])
                if chart == CHARTS[0]:
                    # The zero inlet is simplified away by the source graph;
                    # the outgoing root is exactly the native contribution.
                    if outgoing[order] != contribution[order] or out_row.get("operation") != "sum":
                        raise ValueError("C2 zero-inlet outgoing topology changed")
                    handles.update(child for child in out_row.get("arguments", []) if isinstance(child, int))
                else:
                    if out_row.get("operation") != "sum" or not isinstance(out_row.get("arguments"), list) or len(out_row["arguments"]) != 2:
                        raise ValueError("C2 outgoing is not a memory-plus-integral sum")
                    product = _row(graph, out_row["arguments"][0])
                    integral_sum = _row(graph, out_row["arguments"][1])
                    if product.get("operation") != "product" or product.get("arguments", [None, None])[1] != incoming[order] or integral_sum.get("operation") != "sum":
                        raise ValueError("C2 outgoing topology changed")
                    if contribution[order] not in integral_sum.get("arguments", []):
                        raise ValueError("C2 contribution root changed")
                    handles.update(child for child in out_row.get("arguments", []) if isinstance(child, int))
                    handles.update(child for child in product.get("arguments", []) if isinstance(child, int))
                    handles.update(child for child in integral_sum.get("arguments", []) if isinstance(child, int))
                _row(graph, contribution[order])
                jacobian_handles: list[int] = []
                kernel_handles: list[int] = []
                for part in pieces:
                    integral_handle = _find_integral(graph, contribution[order], chart, part)
                    if integral_handle is None:
                        raise ValueError("C2 native partition integral missing")
                    handles.add(integral_handle)
                    integral_row = _row(graph, integral_handle)
                    if integral_row.get("lower") != part.get("lower") or integral_row.get("upper") != part.get("upper"):
                        raise ValueError("C2 native partition domain endpoint changed")
                    integrand_handle = integral_row.get("integrand")
                    integrand_row = _row(graph, integrand_handle)
                    args = integrand_row.get("arguments")
                    density_value = window.get("density", {}).get("m", [None, None, None])[order]
                    if integrand_row.get("operation") != "product" or not isinstance(args, list) or len(args) not in (2, 3) or (len(args) == 2 and chart != "restoration") or args[1] != density_value:
                        raise ValueError("C2 native kernel/density/Jacobian product changed")
                    kernel_handle = args[0]
                    jacobian_handle = args[2] if len(args) == 3 else None
                    kernel_row = _row(graph, kernel_handle)
                    kernel_arg = _row(graph, kernel_row.get("argument")) if kernel_row.get("operation") == "analytic_unary" else {}
                    offset_difference = _row(graph, kernel_arg.get("argument")) if kernel_arg.get("operation") == "negative" else {}
                    if kernel_row.get("operation") != "analytic_unary" or kernel_row.get("name") != "exp" or kernel_arg.get("operation") != "negative" or offset_difference.get("operation") != "sum":
                        raise ValueError("C2 native rate-one exponential kernel changed")
                    kernel_handles.append(kernel_handle)
                    if jacobian_handle is not None:
                        jacobian_handles.append(jacobian_handle)
                    handles.update([integrand_handle, kernel_handle, kernel_row.get("argument"), kernel_arg.get("argument")])
                    handles.update(child for child in offset_difference.get("arguments", []) if isinstance(child, int))
                    if jacobian_handle is not None:
                        handles.add(jacobian_handle)
                if kernel_handles and len(set(kernel_handles)) != 1:
                    raise ValueError("C2 partition kernel was not shared")
                if jacobian_handles and len(set(jacobian_handles)) != 1:
                    raise ValueError("C2 physical Jacobian was not applied once and shared")
            except (TypeError, ValueError, IndexError) as exc:
                gates.append({"gate": "c2_saved_transport_topology_changed", "chart": chart, "order": order, "error": str(exc)})
                ok = False
            rows.append({"chart": chart, "ordinary_order": order, "incoming": incoming[order], "outgoing": outgoing[order], "contribution": contribution[order], "memory": memory})
    graph_ok = _copy_rows(nodes, graph, handles, "C2 source prefix", gates)
    ok = ok and archive_flags_ok and graph_ok
    return ok, {"rows": rows, "native_partitions": partition_summary, "archive_graph_length": len(graph), "current_rows_authenticated": graph_ok}


def _find_integral(graph: Sequence[Any], root: int, chart: str, part: Mapping[str, Any]) -> int | None:
    row = _row(graph, root)
    if row.get("operation") == "definite_integral" and row.get("actual_current_chart") == chart and row.get("exact_native_left") == part.get("exact_left") and row.get("exact_native_right") == part.get("exact_right"):
        return root
    if row.get("operation") != "sum":
        return None
    for child in row.get("arguments", []):
        if isinstance(child, int):
            found = _find_integral(graph, child, chart, part)
            if found is not None:
                return found
    return None


def _c5_witness(c5_archive: Mapping[str, Any], nodes: Sequence[Any], gates: list[dict[str, Any]]) -> tuple[bool, dict[str, Any]]:
    graph = c5_archive.get("exact_graph_nodes")
    bindings = c5_archive.get("prescribed_loop_C5_bindings")
    charts = bindings.get("charts") if isinstance(bindings, Mapping) else None
    endpoints = bindings.get("endpoints") if isinstance(bindings, Mapping) else None
    if not isinstance(graph, list) or not isinstance(charts, Mapping) or not isinstance(endpoints, Mapping):
        gates.append({"gate": "c5_archive_shape_changed"})
        return False, {}
    archive_flags_ok = True
    for flag in ("current_original_true_prescribed_ten_chart_signed_loop_C5_functions_installed", "current_original_true_prescribed_ten_chart_finite_N_density_C5_functions_installed", "current_original_true_prescribed_rminus_to_Rm_correction_history_C5_functions_installed"):
        if c5_archive.get(flag) is not True:
            gates.append({"gate": "c5_archive_acceptance_flag_missing", "flag": flag})
            archive_flags_ok = False
    expected_order = list(CHARTS) + [POSTRESTORE]
    if list(charts)[:len(expected_order)] != expected_order:
        gates.append({"gate": "c5_chart_order_changed"})
        return False, {}
    rows: list[dict[str, Any]] = []
    handles: set[int] = set()
    ok = True
    previous_endpoint: list[int] | None = None
    offsets: dict[str, Any] = {}
    for position, chart in enumerate(expected_order):
        value = charts.get(chart)
        if not isinstance(value, Mapping):
            gates.append({"gate": "c5_chart_not_mapping", "chart": chart})
            ok = False
            continue
        if chart == POSTRESTORE:
            # Scope stops at the restoration endpoint.  Authenticate only the
            # incoming map used by the next chart; do not replay or inspect
            # the postrestore history itself.
            incoming = _as_int_list(value.get("incoming", {}).get("m"), 3, "C5 postrestore incoming m", gates)
            if incoming is None or incoming != previous_endpoint:
                gates.append({"gate": "c5_restoration_to_postrestore_endpoint_changed"})
                ok = False
            if incoming is not None:
                handles.update(incoming)
            continue
        incoming = _as_int_list(value.get("incoming", {}).get("m"), 3, f"C5 {chart} incoming m", gates)
        history = _as_int_list(value.get("history", {}).get("m"), 3, f"C5 {chart} history m", gates)
        contribution = _as_int_list(value.get("contributions", {}).get("m"), 3, f"C5 {chart} contribution m", gates)
        endpoint = _as_int_list(endpoints.get(chart, {}).get("m"), 3, f"C5 {chart} endpoint m", gates)
        memory = value.get("memory", {}).get("m")
        coordinate = value.get("coordinate")
        physical_measure = value.get("physical_measure")
        global_offset = value.get("global_offset")
        offsets[chart] = {"coordinate": coordinate, "physical_measure": physical_measure, "global_offset": global_offset}
        if not isinstance(memory, int) or not isinstance(coordinate, int) or not isinstance(physical_measure, int) or not isinstance(global_offset, int):
            gates.append({"gate": "c5_chart_geometry_changed", "chart": chart})
            ok = False
        if incoming is None or history is None or contribution is None or endpoint is None or not isinstance(memory, int):
            ok = False
            continue
        if position == 0:
            zero_rows = [_row(graph, item) for item in incoming]
            if any(row.get("operation") != "exact_rational" or row.get("numerator") != 0 or row.get("denominator") != 1 for row in zero_rows):
                gates.append({"gate": "c5_zero_inlet_changed", "chart": chart})
                ok = False
        elif incoming != previous_endpoint:
            gates.append({"gate": "c5_incoming_endpoint_recurrence_changed", "chart": chart})
            ok = False
        previous_endpoint = endpoint
        if endpoint != _as_int_list(endpoints.get(chart, {}).get("m"), 3, f"C5 endpoint {chart}", gates):
            ok = False
        for order in ORDERS:
            try:
                wrapper = _row(graph, endpoint[order])
                upper = value.get("domain", [None, None])[1]
                if wrapper.get("operation") != "function_substitution" or wrapper.get("expression") != history[order] or wrapper.get("variable") != coordinate or wrapper.get("value") != upper or wrapper.get("Z_independent_substitution") is not True:
                    raise ValueError("C5 endpoint substitution changed")
                hrow = _row(graph, history[order])
                if chart == CHARTS[0]:
                    # At the original rminus inlet the zero incoming term is
                    # simplified away; history and contribution share the
                    # exact definite-integral root.
                    if history[order] != contribution[order] or hrow.get("operation") != "definite_integral":
                        raise ValueError("C5 zero-inlet history topology changed")
                    integral = hrow
                    product = None
                else:
                    if hrow.get("operation") != "sum" or not isinstance(hrow.get("arguments"), list) or len(hrow["arguments"]) != 2:
                        raise ValueError("C5 history is not a sum")
                    product = _row(graph, hrow["arguments"][0])
                    integral = _row(graph, hrow["arguments"][1])
                    if product.get("operation") != "product" or product.get("arguments") != [incoming[order], memory] or integral.get("operation") != "definite_integral" or integral.get("original_chart") != chart or integral.get("original_rate") != "1":
                        raise ValueError("C5 history recurrence changed")
                if integral.get("operation") != "definite_integral" or integral.get("original_chart") != chart or integral.get("original_rate") != "1":
                    raise ValueError("C5 transport integral changed")
                if integral.get("physical_Jacobian_applied_once") is not True or integral.get("zero_correction_inlet_only_at_original_rminus") is not True or integral.get("predecessor_memory_never_reset") is not True or integral.get("endpoint_independent_of_Z") is not True:
                    raise ValueError("C5 transport guard changed")
                domain = value.get("domain")
                if not isinstance(domain, list) or len(domain) != 2 or integral.get("lower") != domain[0] or integral.get("upper") != coordinate or integral.get("variable") != "true_loop_C5_transport_" + chart:
                    raise ValueError("C5 native endpoint or variable changed")
                integrand = integral.get("integrand")
                if not isinstance(integrand, int):
                    raise ValueError("C5 integrand handle missing")
                integrand_row = _row(graph, integrand)
                integrand_args = integrand_row.get("arguments")
                if integrand_row.get("operation") != "product" or not isinstance(integrand_args, list) or len(integrand_args) != 2:
                    raise ValueError("C5 source/Jacobian/kernel product changed")
                physical_measure = value.get("physical_measure")
                body_row = _row(graph, integrand_args[1])
                if body_row.get("operation") != "product" or not isinstance(body_row.get("arguments"), list) or body_row["arguments"].count(physical_measure) != 1:
                    raise ValueError("C5 physical measure was not applied exactly once")
                kernel_candidates = [item for item in body_row["arguments"] if item != physical_measure]
                if len(kernel_candidates) != 1:
                    raise ValueError("C5 kernel/Jacobian product changed")
                kernel_row = _row(graph, kernel_candidates[0])
                kernel_arg = _row(graph, kernel_row.get("argument")) if kernel_row.get("operation") == "analytic_unary" else {}
                kernel_arg_children = kernel_arg.get("arguments") if isinstance(kernel_arg, Mapping) else None
                if kernel_row.get("operation") != "analytic_unary" or kernel_row.get("name") != "exp" or kernel_arg.get("operation") != "product" or not isinstance(kernel_arg_children, list) or physical_measure not in kernel_arg_children or not any(_row(graph, child).get("operation") == "sum" for child in kernel_arg_children if isinstance(child, int)):
                    raise ValueError("C5 rate-one exp[-J*(x-t)] kernel changed")
                handles.add(integrand)
                handles.update([endpoint[order], history[order], contribution[order], incoming[order], memory, integral.get("integrand"), integrand_args[0], integrand_args[1], kernel_candidates[0], kernel_row.get("argument")])
                handles.update(child for child in hrow.get("arguments", []) if isinstance(child, int))
                if product is not None:
                    handles.update(child for child in product.get("arguments", []) if isinstance(child, int))
                handles.update(child for child in body_row.get("arguments", []) if isinstance(child, int))
                handles.update(child for child in kernel_arg_children if isinstance(child, int))
            except (TypeError, ValueError, IndexError) as exc:
                gates.append({"gate": "c5_saved_transport_topology_changed", "chart": chart, "order": order, "error": str(exc)})
                ok = False
            rows.append({"chart": chart, "ordinary_order": order, "incoming": incoming[order], "history": history[order], "contribution": contribution[order], "endpoint": endpoint[order], "memory": memory})
    # The C5 saved chart offsets are the only authenticated source of the
    # physical d(offset)/dx measure used here; retain them without evaluating.
    geometry_ok = True
    for chart in CHARTS:
        value = charts[chart]
        if not (value.get("total_Z_phase_derivative_exact_zero") is True and value.get("physical_measure") is not None and value.get("global_offset") is not None):
            geometry_ok = False
    if not geometry_ok:
        gates.append({"gate": "c5_shared_phase_or_offset_geometry_changed"})
    graph_ok = _copy_rows(nodes, graph, handles, "C5 current source prefix", gates)
    ok = ok and archive_flags_ok and geometry_ok and graph_ok
    return ok, {"rows": rows, "chart_geometry": offsets, "archive_graph_length": len(graph), "current_rows_authenticated": graph_ok}


def _density_premise(source_records: Mapping[str, Any], c2_witness: Mapping[str, Any], c5_witness: Mapping[str, Any], gates: list[dict[str, Any]]) -> tuple[bool, dict[str, Any]]:
    value = None
    if isinstance(source_records, Mapping):
        for key in ("same_density_premise", "density_premise", "density_certificate"):
            candidate = source_records.get(key)
            if isinstance(candidate, Mapping):
                value = candidate
                break
    if not isinstance(value, Mapping):
        gates.append({"gate": "same_density_premise_missing", "required": "explicit provider certificate C2=N*C5"})
        return False, {"present": False, "admitted": False}
    relation = value.get("relation")
    accepted_relation = relation == "C2=N*C5"
    c2_density = value.get("rho_C2")
    c5_density = value.get("rho_C5")
    fields = value.get("fields")
    orders = value.get("ordinary_orders")
    source_graph_ok = value.get("source_graph_digest") == EXPECTED_NODES_DIGEST and value.get("source_graph_node_count") == EXPECTED_NODE_COUNT
    source_hashes = value.get("source_record_hashes")
    hashes_ok = isinstance(source_hashes, Mapping) and all(source_hashes.get(name) == digest for name, digest in RECORD_HASHES.items())
    certificate_digest = value.get("certificate_sha256", value.get("certificate_digest"))
    digest_body = copy.deepcopy(dict(value))
    digest_body.pop("certificate_sha256", None)
    digest_body.pop("certificate_digest", None)
    certificate_digest_ok = isinstance(certificate_digest, str) and certificate_digest == canonical_digest(digest_body)
    field_ok = fields == ["m"] or fields == "m"
    order_ok = orders == list(ORDERS)
    bindings = value.get("per_chart_bindings")
    expected_c2 = {(row.get("chart"), row.get("ordinary_order")): row for row in c2_witness.get("rows", []) if isinstance(row, Mapping)}
    expected_c5 = {(row.get("chart"), row.get("ordinary_order")): row for row in c5_witness.get("rows", []) if isinstance(row, Mapping)}
    bindings_ok = isinstance(bindings, Mapping) and set(bindings) == set(CHARTS)
    if bindings_ok:
        for chart in CHARTS:
            binding = bindings.get(chart)
            if not isinstance(binding, Mapping):
                bindings_ok = False
                break
            if binding.get("field") != "m" or binding.get("ordinary_orders") != list(ORDERS) or binding.get("relation") != "C2=N*C5" or binding.get("rho_C2") != "B" or binding.get("rho_C5") != "B_over_Pstar/N":
                bindings_ok = False
                break
            order_bindings = binding.get("orders") if isinstance(binding.get("orders"), Mapping) else binding
            for order in ORDERS:
                left, right = expected_c2.get((chart, order), {}), expected_c5.get((chart, order), {})
                expected_binding = {
                    "c2_incoming": left.get("incoming"), "c2_outgoing": left.get("outgoing"),
                    "c5_incoming": right.get("incoming"), "c5_history": right.get("history"), "c5_endpoint": right.get("endpoint"),
                }
                if order_bindings.get(str(order), order_bindings.get(order)) != expected_binding:
                    bindings_ok = False
                    break
            if not bindings_ok:
                break
    structurally_complete = bool(accepted_relation and c2_density == "B" and c5_density == "B_over_Pstar/N" and field_ok and order_ok and source_graph_ok and hashes_ok and certificate_digest_ok and bindings_ok)
    # A self-digest authenticates only the supplied object, not a theorem.
    # No independently pinned actual provider-identity receipt exists here.
    admitted = False
    gates.append({"gate": "independent_density_provider_receipt_not_pinned",
                  "self_described_certificate_complete": structurally_complete})
    if not admitted:
        gates.append({"gate": "same_density_premise_not_admitted", "required": {"relation": "C2=N*C5", "rho_C2": "B", "rho_C5": "B_over_Pstar/N", "field": ["m"], "ordinary_orders": list(ORDERS), "source_graph_digest": EXPECTED_NODES_DIGEST, "full_per_chart_bindings": True, "certificate_sha256": "canonical digest of certificate body excluding certificate_sha256"}})
    return admitted, {"present": True, "admitted": admitted, "relation": relation, "rho_C2": c2_density, "rho_C5": c5_density, "fields": copy.deepcopy(fields), "ordinary_orders": copy.deepcopy(orders), "source_graph_digest": value.get("source_graph_digest"), "source_graph_node_count": value.get("source_graph_node_count"), "source_record_hashes": copy.deepcopy(source_hashes), "certificate_sha256": certificate_digest, "certificate_digest_matches": certificate_digest_ok, "full_per_chart_bindings": bindings_ok, "provider": value.get("provider")}


def _relation_rows(c2: Mapping[str, Any], c5: Mapping[str, Any]) -> list[dict[str, Any]]:
    c2_by = {(row["chart"], row["ordinary_order"]): row for row in c2.get("rows", [])}
    c5_by = {(row["chart"], row["ordinary_order"]): row for row in c5.get("rows", [])}
    rows = []
    for chart in CHARTS:
        for order in ORDERS:
            left, right = c2_by.get((chart, order), {}), c5_by.get((chart, order), {})
            rows.append({
                "chart": chart,
                "field": M_COMPONENT,
                "ordinary_order": order,
                "c2_incoming": left.get("incoming"),
                "c2_outgoing": left.get("outgoing"),
                "c5_incoming": right.get("incoming"),
                "c5_history": right.get("history"),
                "c5_endpoint": right.get("endpoint"),
                "same_density_relation": "conditional C2=N*C5",
                "history_relation_formula": "H_C2=N*H_C5",
                "relation_proved": False,
            })
    return rows


def derive(nodes: Sequence[Any], source_records: Mapping[str, Any], *, expected_nodes_digest: str) -> dict[str, Any]:
    """Return a fail-closed source certificate for the first nine m charts."""

    _require(expected_nodes_digest == EXPECTED_NODES_DIGEST, "unexpected source graph digest assertion")
    _require(len(nodes) == EXPECTED_NODE_COUNT and canonical_digest(nodes) == EXPECTED_NODES_DIGEST, "saved source graph changed")
    gates: list[dict[str, Any]] = []
    source_ok, source_witness, parsed = _source_witnesses(source_records, gates)
    c2_ok = c5_ok = False
    c2_witness: dict[str, Any] = {}
    c5_witness: dict[str, Any] = {}
    if parsed:
        c2_ok, c2_witness = _c2_witness(parsed.get(C2_DENSITY_GZ, {}), nodes, gates)
        c5_ok, c5_witness = _c5_witness(parsed.get(C5_GZ, {}), nodes, gates)
    density_ok, density_witness = _density_premise(source_records, c2_witness, c5_witness, gates)
    structural_ok = bool(source_ok and c2_ok and c5_ok)
    gates.extend([{"gate": "transport_offset_kernel_equality_unbound",
        "required": "exact rate-one operands, native coordinate/integration variable and cross-archive offset/Jacobian identity"},
        {"gate": "native_partition_additivity_unbound",
        "required": "each contribution root equals the sum of exactly the adjacent partition integrals with shared chart endpoints"}])
    passed = False
    c2_rows = c2_witness.get("rows", [])
    c5_rows = c5_witness.get("rows", [])
    theorem = {
        "premise": "explicit same linear density rho_C2=B and rho_C5=B_over_Pstar/N, hence C2=N*C5",
        "conclusion": "for m and ordinary Z orders 0,1,2, corresponding native chart histories satisfy H_C2=N*H_C5 through restoration",
        "C2_kernel": "exp[-(right_offset-offset(t))]",
        "C5_kernel": "exp[-J*(x-t)]",
        "shared_phase_offset_rule": "d(offset)/dx=physical_measure; phase is held fixed in ordinary Z rows",
        "physical_Jacobian_applied_once": True,
        "zero_inlet_only_at_original_rminus": True,
        "predecessor_memory_never_reset": True,
        "native_partition_endpoint_additivity_required": True,
        "proved": False,
        "numeric_values_computed": False,
    }
    return {
        "schema": SCHEMA,
        "source_only": True,
        "nodes_digest": EXPECTED_NODES_DIGEST,
        "nodes_count": EXPECTED_NODE_COUNT,
        "scope": {"field": M_COMPONENT, "ordinary_orders": list(ORDERS), "charts_through_restoration": list(CHARTS), "postrestore_endpoint_checked_on_C5": True, "full_history_C2_claimed": False, "Qy1_claimed": False},
        "passed": passed,
        "structural_source_rows_checked": structural_ok,
        "transport_calculus_checked": False,
        "conditional_transport_calculus_supported": False,
        "conditional_transport_theorem": theorem,
        "conditional_transport_calculus": theorem,
        "density_premise": density_witness,
        "rows": _relation_rows(c2_witness, c5_witness),
        "restoration_edges": {
            "C2_reference_to_restoration": [2323, 2324, 11099],
            "C2_restoration_to_post": [2607, 2608, 11202],
            "C5_reference_to_restoration_incoming": [116764, 116765, 116766],
            "C5_restoration_history": [119228, 119229, 119230],
            "C5_restoration_to_post_incoming": [119382, 119383, 119384],
            "relation_required": "C2=N*C5; handles are distinct source maps",
            "relation_proved": False,
        },
        "witnesses": {"producer_sources": source_witness, "C2_saved_transport": c2_witness, "C5_saved_transport": c5_witness, "source_chain_structural": source_ok, "C2_graph_structural": c2_ok, "C5_graph_structural": c5_ok},
        "unresolved_gates": gates,
        "numeric_values_computed": False,
        "old_operators_replayed": False,
        "operators_reexecuted": False,
        "full_goal_complete": False,
        "full_Qy1_or_full_goal_claimed": False,
        "source_record_hashes": dict(RECORD_HASHES),
    }
