"""Source-only exact native transport checks for the first nine m charts.

The checker authenticates the saved C2 and C5 source projections, then checks
their graph topology without importing, executing, or replaying any producer.
The C2 contribution roots are required to be the exact flattened sums of their
connected native partition integrals.  The C5 histories are checked as the
non-reset recurrence with the original-rminus zero inlet.  Kernel operands are
bound by graph handle and operation, including the integration variable and
the physical Jacobian.

The saved source files are inputs only.  A gzip record is authenticated from
its compressed bytes before decoding.  Producer archive flags are retained as
metadata and deliberately do not control acceptance: the separate pinned
check JSON records are the source acceptance receipts.  The module reports a
conditional transport theorem, but never admits the C2 = N C5 history
relation from a self-described density object or a bare ``passed`` flag.
"""

from __future__ import annotations

import ast
import base64
import binascii
import copy
from collections import Counter
from collections.abc import Mapping, Sequence
import gzip
import hashlib
import json
from typing import Any


SCHEMA = "original_Rm_m_C2_C5_native_transport_v1"
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

# These are the accepted parent input pins.  In particular, the check JSON
# receipts are pinned separately from the candidate archive flags.
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
    """Return the pinned JSON digest without evaluating any source object."""

    raw = json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
        allow_nan=False,
    ).encode("utf-8")
    return hashlib.sha256(raw).hexdigest()


def _gate(gates: list[dict[str, Any]], name: str, **details: Any) -> None:
    item: dict[str, Any] = {"gate": name}
    item.update(details)
    gates.append(item)


def _mapping(value: Any) -> Mapping[str, Any]:
    """Adapt a descriptor field only after checking that it is a mapping."""

    return value if isinstance(value, Mapping) else {}


def _descriptor(value: Any, label: str, gates: list[dict[str, Any]]) -> Mapping[str, Any]:
    """Require a named archive descriptor to be a mapping before field access."""

    if isinstance(value, Mapping):
        return value
    _gate(gates, "descriptor_field_not_mapping", field=label)
    return {}


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _row(nodes: Sequence[Any], index: int) -> Mapping[str, Any]:
    _require(type(index) is int and 0 <= index < len(nodes), f"graph handle out of range: {index}")
    value = nodes[index]
    _require(isinstance(value, Mapping), f"graph handle is not a mapping: {index}")
    return value


def _bytes_from_record(record: Mapping[str, Any]) -> bytes | None:
    """Read only the runner's raw-byte fields; no path is opened here."""

    encoded = record.get("raw_bytes_base64")
    if isinstance(encoded, str):
        try:
            return base64.b64decode(encoded, validate=True)
        except (ValueError, binascii.Error):
            return None
    for key in ("raw_bytes", "source_bytes", "gzip_bytes"):
        value = record.get(key)
        if isinstance(value, bytes):
            return value
        if isinstance(value, bytearray):
            return bytes(value)
    return None


def _record(source_records: Mapping[str, Any], name: str, gates: list[dict[str, Any]]) -> dict[str, Any] | None:
    """Authenticate one source text record, hashing gzip bytes before decode."""

    expected = RECORD_HASHES[name]
    record = source_records.get(name) if isinstance(source_records, Mapping) else None
    if not isinstance(record, Mapping):
        _gate(gates, "source_record_missing", record=name)
        return None
    if not isinstance(record.get("path"), str) or not isinstance(record.get("sha256"), str):
        _gate(gates, "source_descriptor_field_type_changed", record=name, fields=("path", "sha256"))
        return None
    expected_path = "experiments/root_st073/" + name
    if record.get("path") != expected_path:
        _gate(gates, "source_record_path_changed", record=name, expected=expected_path)
        return None
    if record.get("sha256") != expected:
        _gate(gates, "source_record_hash_changed", record=name, expected_sha256=expected)
        return None

    if name.endswith(".gz"):
        if not isinstance(record.get("raw_bytes_base64"), str) and not any(
            isinstance(record.get(key), (bytes, bytearray)) for key in ("raw_bytes", "source_bytes", "gzip_bytes")
        ):
            _gate(gates, "gzip_raw_descriptor_field_type_changed", record=name)
            return None
        if not isinstance(record.get("source_text"), str) or not isinstance(record.get("decoded_sha256"), str):
            _gate(gates, "gzip_decoded_descriptor_field_type_changed", record=name, fields=("source_text", "decoded_sha256"))
            return None
        raw = _bytes_from_record(record)
        if raw is None:
            _gate(gates, "gzip_raw_bytes_missing", record=name)
            return None
        raw_sha = hashlib.sha256(raw).hexdigest()
        if raw_sha != expected:
            _gate(gates, "gzip_raw_hash_changed", record=name, expected_sha256=expected)
            return None
        try:
            decoded = gzip.decompress(raw)
        except (OSError, EOFError) as exc:
            _gate(gates, "gzip_decode_failed", record=name, error=str(exc))
            return None
        try:
            text = decoded.decode("utf-8")
        except UnicodeDecodeError as exc:
            _gate(gates, "gzip_utf8_decode_failed", record=name, error=str(exc))
            return None
        if record.get("source_text") != text:
            _gate(gates, "gzip_decoded_text_changed", record=name)
            return None
        if record.get("decoded_sha256") != hashlib.sha256(decoded).hexdigest():
            _gate(gates, "gzip_decoded_hash_changed", record=name)
            return None
        return {
            "name": name,
            "sha256": expected,
            "text": text,
            "decoded": decoded,
            "raw_sha256": raw_sha,
        }

    text = record.get("source_text")
    if not isinstance(text, str):
        _gate(gates, "source_text_missing", record=name)
        return None
    encoded = text.encode("utf-8")
    if hashlib.sha256(encoded).hexdigest() != expected:
        _gate(gates, "source_text_hash_changed", record=name, expected_sha256=expected)
        return None
    return {"name": name, "sha256": expected, "text": text, "decoded": encoded, "raw_sha256": None}


def _definitions(text: str, name: str, wanted: tuple[str, ...], gates: list[dict[str, Any]]) -> dict[str, dict[str, int]]:
    try:
        tree = ast.parse(text, filename=name)
    except SyntaxError as exc:
        _gate(gates, "source_ast_invalid", record=name, error=str(exc))
        return {}
    found: dict[str, dict[str, int]] = {}
    for item in ast.walk(tree):
        if isinstance(item, (ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)) and item.name in wanted:
            found[item.name] = {"line": item.lineno, "end_line": item.end_lineno or item.lineno}
    missing = [item for item in wanted if item not in found]
    if missing:
        _gate(gates, "source_declaration_missing", record=name, declarations=missing)
    return found


def _tokens(text: str, name: str, wanted: tuple[str, ...], gates: list[dict[str, Any]]) -> bool:
    missing = [token for token in wanted if token not in text]
    if missing:
        _gate(gates, "source_constructor_token_missing", record=name, tokens=missing)
        return False
    return True


def _call_arities(text: str, name: str, requirements: Mapping[str, tuple[int, ...]], gates: list[dict[str, Any]]) -> bool:
    """Check static direct source calls without invoking the source AST."""

    try:
        tree = ast.parse(text, filename=name)
    except SyntaxError:
        return False
    seen: dict[str, set[int]] = {key: set() for key in requirements}
    for item in ast.walk(tree):
        if not isinstance(item, ast.Call) or not isinstance(item.func, ast.Attribute):
            continue
        owner = item.func.value
        if not isinstance(owner, ast.Name) or owner.id != "g":
            continue
        key = item.func.attr
        if key in seen:
            if any(isinstance(keyword, ast.Starred) for keyword in item.args):
                _gate(gates, "source_dynamic_call_arguments", record=name, call=f"g.{key}")
                continue
            seen[key].add(len(item.args))
    ok = True
    for key, expected in requirements.items():
        if not seen[key].intersection(expected):
            _gate(gates, "source_static_call_arity_missing", record=name, call=f"g.{key}", expected=list(expected))
            ok = False
    return ok


def _json_record(record: Mapping[str, Any], name: str, gates: list[dict[str, Any]]) -> Mapping[str, Any]:
    try:
        value = json.loads(record["text"])
    except (KeyError, TypeError, json.JSONDecodeError) as exc:
        _gate(gates, "accepted_source_receipt_invalid", record=name, error=str(exc))
        return {}
    if not isinstance(value, Mapping):
        _gate(gates, "accepted_source_receipt_not_mapping", record=name)
        return {}
    return value


def _flag(value: Mapping[str, Any], path: tuple[str, ...]) -> Any:
    current: Any = value
    for key in path:
        current = current.get(key) if isinstance(current, Mapping) else None
    return current


def _receipt(record: Mapping[str, Any], name: str, required: tuple[tuple[str, ...], ...], gates: list[dict[str, Any]]) -> tuple[bool, Mapping[str, Any]]:
    value = _json_record(record, name, gates)
    ok = value.get("all_passed") is True
    if not ok:
        _gate(gates, "accepted_source_receipt_header_missing", record=name)
    for path in required:
        if _flag(value, path) is not True:
            _gate(gates, "accepted_source_receipt_flag_missing", record=name, path=".".join(path))
            ok = False
    return ok, value


def _source_witnesses(source_records: Mapping[str, Any], gates: list[dict[str, Any]]) -> tuple[bool, dict[str, Any], dict[str, Any]]:
    records = {name: _record(source_records, name, gates) for name in RECORD_HASHES}
    if any(value is None for value in records.values()):
        return False, {"records": dict(RECORD_HASHES)}, {}
    if any(not isinstance(value, dict) for value in records.values()):
        _gate(gates, "source_record_descriptor_invalid")
        return False, {"records": dict(RECORD_HASHES)}, {}

    declarations = {
        C2_PHASE_PY: _definitions(records[C2_PHASE_PY]["text"], C2_PHASE_PY, ("C2Algebra", "primitive_C2", "source_bindings", "CurrentC2PhasePrimitives"), gates),
        C2_DENSITY_PY: _definitions(records[C2_DENSITY_PY]["text"], C2_DENSITY_PY, ("coefficient_C2", "build", "CurrentC2DensityTransport"), gates),
        C5_PY: _definitions(records[C5_PY]["text"], C5_PY, ("LoopAlgebra", "loop_functions", "finite_N_functions", "transport", "CurrentOriginalPrescribedLoopC5Functions"), gates),
    }
    token_ok = all(
        (
            _tokens(records[C2_PHASE_PY]["text"], C2_PHASE_PY, ("source_bindings", "ordinary", "Z2", "inverse"), gates),
            _tokens(records[C2_DENSITY_PY]["text"], C2_DENSITY_PY, ("g.unary('exp'", "g.sub(ref(old['right_offset']),ref(old['offset']))", "g.mul(kernel,q,ref(old['Jacobian']))", "N_scaled_density=True", "ordinary_slow_Z_derivative_order=2", "original dlogR: Jacobian exactly once", "integration_endpoints_kernel_radius_and_global_phase_Z_independent=True"), gates),
            _tokens(records[C5_PY]["text"], C5_PY, ("g.unary('exp'", "g.mul(g.constant(-rate),jacobian,g.sub(coordinate,t))", "physical_Jacobian_applied_once=True", "zero_correction_inlet_only_at_original_rminus=True", "predecessor_memory_never_reset=True"), gates),
        )
    )
    api_ok = all(
        (
            _call_arities(records[C2_DENSITY_PY]["text"], C2_DENSITY_PY, {"sub": (2,), "mul": (3,), "unary": (2,)}, gates),
            _call_arities(records[C5_PY]["text"], C5_PY, {"sub": (2,), "mul": (3,), "unary": (2,)}, gates),
        )
    )
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
    # A quiet original-power branch is admissible only when the pinned source
    # receipt itself carries the theorem.  Archive candidate flags and
    # arbitrary caller booleans never authorize a zero increment.
    zero_source_admitted = bool(
        _flag(phase_receipt, ("original_power_zero_increment_source_admitted",)) is True
        or _flag(density_receipt, ("original_power_zero_increment_source_admitted",)) is True
    )

    parsed: dict[str, Any] = {}
    for name in (C2_PHASE_GZ, C2_DENSITY_GZ, C5_GZ):
        try:
            value = json.loads(records[name]["text"])
        except json.JSONDecodeError as exc:
            _gate(gates, "source_archive_json_invalid", record=name, error=str(exc))
            continue
        if not isinstance(value, Mapping):
            _gate(gates, "source_archive_json_not_mapping", record=name)
            continue
        parsed[name] = value

    structural_ok = bool(
        api_ok
        and token_ok
        and phase_ok
        and density_ok
        and c5_ok
        and all(bool(value) for value in declarations.values())
        and len(parsed) == 3
    )
    witness = {
        "records": dict(RECORD_HASHES),
        "raw_gzip_hashes_verified_before_decode": {name: records[name]["raw_sha256"] for name in (C2_PHASE_GZ, C2_DENSITY_GZ, C5_GZ)},
        "declarations": declarations,
        "static_direct_api_arity_checked": api_ok,
        "dynamic_mapping_checks": True,
        "accepted_receipts": {
            C2_PHASE_CHECK: {"sha256": RECORD_HASHES[C2_PHASE_CHECK], "passed": phase_ok},
            C2_DENSITY_CHECK: {"sha256": RECORD_HASHES[C2_DENSITY_CHECK], "passed": density_ok},
            C5_CHECK: {"sha256": RECORD_HASHES[C5_CHECK], "passed": c5_ok},
        },
        "candidate_archive_flags_ignored": True,
        "original_power_zero_increment_source_admitted": zero_source_admitted,
        "parent_report_sha256": PARENT_REPORT_SHA256,
        "accepted_a_positivity_report_sha256": POSITIVE_REPORT_SHA256,
        "source_chain_structural": structural_ok,
    }
    return structural_ok, witness, parsed


_REFERENCE_KEYS = (
    "arguments",
    "argument",
    "integrand",
    "expression",
    "variable",
    "value",
    "lower",
    "upper",
    "active_body",
    "Delta",
    "eta",
    "flat_predicate",
    "original_collar_coordinate",
    "original_collar_upper",
)


def _collect_graph_refs(graph: Sequence[Any], index: int, handles: set[int], seen: set[int] | None = None) -> None:
    if type(index) is not int or index in handles:
        return
    if index < 0 or index >= len(graph) or not isinstance(graph[index], Mapping):
        return
    handles.add(index)
    seen = seen or set()
    if index in seen:
        return
    seen.add(index)
    row = graph[index]
    for key in _REFERENCE_KEYS:
        value = row.get(key)
        if type(value) is int:
            _collect_graph_refs(graph, value, handles, seen)
        elif isinstance(value, list):
            for item in value:
                if type(item) is int:
                    _collect_graph_refs(graph, item, handles, seen)


def _copy_rows(nodes: Sequence[Any], reference: Sequence[Any], handles: set[int], label: str, gates: list[dict[str, Any]]) -> bool:
    ok = True
    for index in sorted(handles):
        try:
            actual = _row(nodes, index)
            expected = _row(reference, index)
        except ValueError as exc:
            _gate(gates, "saved_graph_handle_invalid", label=label, handle=index, error=str(exc))
            ok = False
            continue
        if actual != expected:
            _gate(gates, "current_graph_source_row_changed", label=label, handle=index)
            ok = False
    return ok


def _as_int_list(value: Any, length: int, label: str, gates: list[dict[str, Any]]) -> list[int] | None:
    if not isinstance(value, list) or len(value) != length or any(type(item) is not int for item in value):
        _gate(gates, "source_handle_vector_changed", label=label, expected_length=length)
        return None
    return value


def _is_zero_row(row: Mapping[str, Any]) -> bool:
    return row.get("operation") == "exact_rational" and row.get("numerator") == 0 and row.get("denominator") == 1


def _all_zero_handles(graph: Sequence[Any], handles: Sequence[int], label: str, gates: list[dict[str, Any]]) -> bool:
    try:
        return all(_is_zero_row(_row(graph, item)) for item in handles)
    except (TypeError, ValueError, IndexError) as exc:
        _gate(gates, "zero_inlet_handle_invalid", label=label, error=str(exc))
        return False


def _flatten_add(graph: Sequence[Any], root: int, gates: list[dict[str, Any]], label: str, seen: set[int] | None = None) -> list[int] | None:
    """Flatten graph ``sum`` nodes and erase exact zeros, preserving order."""

    if type(root) is not int or root < 0 or root >= len(graph):
        _gate(gates, "additive_graph_handle_invalid", label=label, handle=root)
        return None
    seen = seen or set()
    if root in seen:
        _gate(gates, "additive_graph_cycle", label=label, handle=root)
        return None
    seen.add(root)
    row = graph[root]
    if not isinstance(row, Mapping):
        _gate(gates, "additive_graph_row_invalid", label=label, handle=root)
        return None
    if row.get("operation") == "sum":
        args = row.get("arguments")
        if not isinstance(args, list) or not args:
            _gate(gates, "additive_sum_arguments_invalid", label=label, handle=root)
            return None
        flattened: list[int] = []
        for child in args:
            if type(child) is not int:
                _gate(gates, "additive_sum_child_invalid", label=label, handle=root)
                return None
            child_row = graph[child] if 0 <= child < len(graph) else {}
            if isinstance(child_row, Mapping) and _is_zero_row(child_row):
                continue
            child_terms = _flatten_add(graph, child, gates, label, seen.copy())
            if child_terms is None:
                return None
            flattened.extend(child_terms)
        return flattened
    if _is_zero_row(row):
        return []
    return [root]


def _flatten_product(graph: Sequence[Any], root: int, gates: list[dict[str, Any]], label: str, seen: set[int] | None = None) -> list[int] | None:
    if type(root) is not int or root < 0 or root >= len(graph):
        _gate(gates, "multiplicative_graph_handle_invalid", label=label, handle=root)
        return None
    seen = seen or set()
    if root in seen:
        _gate(gates, "multiplicative_graph_cycle", label=label, handle=root)
        return None
    seen.add(root)
    row = graph[root]
    if not isinstance(row, Mapping):
        _gate(gates, "multiplicative_graph_row_invalid", label=label, handle=root)
        return None
    if row.get("operation") != "product":
        return [root]
    args = row.get("arguments")
    if not isinstance(args, list) or not args:
        _gate(gates, "multiplicative_product_arguments_invalid", label=label, handle=root)
        return None
    result: list[int] = []
    for child in args:
        if type(child) is not int:
            _gate(gates, "multiplicative_product_child_invalid", label=label, handle=root)
            return None
        terms = _flatten_product(graph, child, gates, label, seen.copy())
        if terms is None:
            return None
        result.extend(terms)
    return result


def _is_exact_zero_expression(graph: Sequence[Any], root: int, gates: list[dict[str, Any]], label: str) -> bool:
    """Recognize literal zero and sums that flatten to only literal zeros."""

    return _flatten_add(graph, root, gates, label) == []


def _find_bound_variable(graph: Sequence[Any], name: str) -> int | None:
    for index, value in enumerate(graph):
        if isinstance(value, Mapping) and value.get("operation") == "bound_variable" and value.get("name") == name:
            return index
    return None


def _contains_handle(graph: Sequence[Any], root: int, wanted: int, seen: set[int] | None = None) -> bool:
    if type(root) is not int or root < 0 or root >= len(graph):
        return False
    if root == wanted:
        return True
    seen = seen or set()
    if root in seen:
        return False
    seen.add(root)
    row = graph[root]
    if not isinstance(row, Mapping):
        return False
    for key in _REFERENCE_KEYS:
        value = row.get(key)
        if type(value) is int and _contains_handle(graph, value, wanted, seen):
            return True
        if isinstance(value, list) and any(type(item) is int and _contains_handle(graph, item, wanted, seen) for item in value):
            return True
    return False


def _c2_kernel_binding(
    graph: Sequence[Any],
    integral: int,
    density_handle: int,
    parts: Sequence[Mapping[str, Any]],
    order: int,
    chart: str,
    gates: list[dict[str, Any]],
) -> tuple[bool, dict[str, Any]]:
    """Check one C2 integral's exact exp[-(right_offset-offset(t))] operands."""

    info: dict[str, Any] = {"integral": integral, "ordinary_order": order, "chart": chart}
    try:
        irow = _row(graph, integral)
        integrand_handle = irow.get("integrand")
        integrand = _row(graph, integrand_handle)
        args = integrand.get("arguments")
        if integrand.get("operation") != "product" or not isinstance(args, list):
            raise ValueError("C2 integral integrand is not a product")
        kernels = [item for item in args if type(item) is int and isinstance(_row(graph, item), Mapping) and _row(graph, item).get("operation") == "analytic_unary" and _row(graph, item).get("name") == "exp"]
        if len(kernels) != 1 or density_handle not in args:
            raise ValueError("C2 kernel or density operand changed")
        if args.count(density_handle) != 1 or len(args) not in (2, 3):
            raise ValueError("C2 density/Jacobian factor count changed")
        kernel_handle = kernels[0]
        kernel = _row(graph, kernel_handle)
        negative_handle = kernel.get("argument")
        negative = _row(graph, negative_handle)
        difference_handle = negative.get("argument") if negative.get("operation") == "negative" else None
        difference = _row(graph, difference_handle)
        diff_args = difference.get("arguments")
        if negative.get("operation") != "negative" or difference.get("operation") != "sum" or not isinstance(diff_args, list) or len(diff_args) != 2:
            raise ValueError("C2 exponential operand is not a signed offset difference")
        terminal = parts[-1].get("right_offset") if parts else None
        negative_terms = [item for item in diff_args if type(item) is int and _row(graph, item).get("operation") == "negative"]
        positive_terms = [item for item in diff_args if item != negative_terms[0] and type(item) is int] if negative_terms else []
        if len(negative_terms) != 1 or len(positive_terms) != 1 or positive_terms[0] != terminal:
            raise ValueError("C2 right_offset operand is not the terminal native offset")
        offset_handle = _row(graph, negative_terms[0]).get("argument")
        variable_name = irow.get("variable")
        variable_handle = _find_bound_variable(graph, variable_name) if isinstance(variable_name, str) else None
        if type(offset_handle) is not int or variable_handle is None or not _contains_handle(graph, offset_handle, variable_handle):
            raise ValueError("C2 offset(t) is not bound to the integration variable")
        jacobians = [item for item in args if item not in (density_handle, kernel_handle)]
        if len(jacobians) > 1:
            raise ValueError("C2 physical Jacobian appears more than once")
        jacobian = jacobians[0] if jacobians else None
        derivative_direct = jacobian is not None and _contains_handle(graph, offset_handle, jacobian)
        info.update(
            {
                "kernel": kernel_handle,
                "negative_exponent": negative_handle,
                "difference": difference_handle,
                "right_offset": terminal,
                "offset_of_t": offset_handle,
                "integration_variable": variable_handle,
                "jacobian": jacobian,
                "jacobian_factor_count": len(jacobians),
                "unit_jacobian_absorbed": jacobian is None,
                "offset_derivative_directly_bound_to_jacobian": derivative_direct,
            }
        )
        return True, info
    except (TypeError, ValueError, IndexError) as exc:
        _gate(gates, "c2_kernel_operands_changed", chart=chart, ordinary_order=order, error=str(exc))
        info["error"] = str(exc)
        return False, info


def _c2_witness(
    c2_archive: Mapping[str, Any],
    nodes: Sequence[Any],
    gates: list[dict[str, Any]],
    *,
    zero_source_admitted: bool = False,
) -> tuple[bool, dict[str, Any]]:
    target = c2_archive.get("actual_C2_density_transport_and_target_functions")
    graph = c2_archive.get("exact_graph_nodes")
    # Candidate/archive flags are descriptive only.  Separate check JSON
    # receipts in _source_witnesses control source acceptance.
    if not isinstance(target, Mapping) or not isinstance(graph, list):
        _gate(gates, "c2_archive_shape_changed")
        return False, {}
    windows = target.get("windows")
    if not isinstance(windows, list) or [item.get("chart") for item in windows[: len(CHARTS)] if isinstance(item, Mapping)] != list(CHARTS):
        _gate(gates, "c2_first_nine_chart_order_changed")
        return False, {}

    rows: list[dict[str, Any]] = []
    handles: set[int] = set()
    ok = True
    partition_ok = True
    kernel_ok = True
    derivative_direct = True
    previous_outgoing: list[int] | None = None
    previous_parts: Sequence[Mapping[str, Any]] | None = None
    cross_chart_offset_handles: list[dict[str, Any]] = []
    partition_summary: dict[str, Any] = {}
    kernel_summary: dict[str, Any] = {}
    zero_branches: dict[str, bool] = {}

    for position, chart in enumerate(CHARTS):
        window = windows[position]
        if not isinstance(window, Mapping):
            _gate(gates, "c2_window_not_mapping", chart=chart)
            ok = False
            partition_ok = False
            kernel_ok = False
            continue
        incoming = _as_int_list(_descriptor(window.get("incoming"), f"C2 {chart}.incoming", gates).get(M_COMPONENT), 3, f"C2 {chart} incoming m", gates)
        outgoing = _as_int_list(_descriptor(window.get("outgoing"), f"C2 {chart}.outgoing", gates).get(M_COMPONENT), 3, f"C2 {chart} outgoing m", gates)
        contribution = _as_int_list(_descriptor(window.get("contributions"), f"C2 {chart}.contributions", gates).get(M_COMPONENT), 3, f"C2 {chart} contribution m", gates)
        density = _as_int_list(_descriptor(window.get("density"), f"C2 {chart}.density", gates).get(M_COMPONENT), 3, f"C2 {chart} density m", gates)
        memory = _descriptor(window.get("original_own_rate_memory"), f"C2 {chart}.original_own_rate_memory", gates).get(M_COMPONENT)
        if not isinstance(memory, int):
            _gate(gates, "c2_memory_handle_changed", chart=chart)
            ok = False
        if incoming is None or outgoing is None or contribution is None or density is None or not isinstance(memory, int):
            ok = False
            partition_ok = False
            kernel_ok = False
            continue
        if window.get("predecessor_ZZ_memory_retained") is not True or window.get("quiet_local_zero_does_not_reset_incoming") is not True:
            _gate(gates, "c2_memory_or_zero_inlet_guard_changed", chart=chart)
            ok = False
        if position == 0:
            if not _all_zero_handles(graph, incoming, f"C2 {chart} incoming", gates):
                _gate(gates, "c2_zero_inlet_changed", chart=chart)
                ok = False
        elif incoming != previous_outgoing:
            _gate(gates, "c2_incoming_recurrence_changed", chart=chart)
            ok = False
        previous_outgoing = outgoing

        parts = window.get("original_physical_measure_and_partitions")
        expected_count = 4 if chart in ("first_micro", "second_micro") else 2
        if not isinstance(parts, list) or len(parts) != expected_count or any(not isinstance(item, Mapping) for item in parts):
            _gate(gates, "c2_native_partition_count_changed", chart=chart, expected=expected_count)
            ok = False
            partition_ok = False
            parts = []
        for part in parts:
            exact_left = part.get("exact_left")
            exact_right = part.get("exact_right")
            if (
                not all(type(part.get(key)) is int for key in ("lower", "upper", "left_offset", "right_offset"))
                or not isinstance(exact_left, list)
                or not isinstance(exact_right, list)
                or not exact_left
                or not exact_right
                or any(not isinstance(item, (int, str)) for item in exact_left + exact_right)
            ):
                _gate(gates, "c2_partition_handle_changed", chart=chart)
                ok = False
                partition_ok = False
        for left, right in zip(parts, parts[1:]):
            if left.get("right_offset") != right.get("left_offset") or left.get("exact_right") != right.get("exact_left"):
                _gate(gates, "c2_native_partition_endpoint_additivity_changed", chart=chart)
                ok = False
                partition_ok = False
        if previous_parts is not None and parts and previous_parts:
            if previous_parts[-1].get("right_offset") != parts[0].get("left_offset"):
                _gate(gates, "c2_cross_chart_offset_additivity_changed", chart=chart)
                ok = False
                partition_ok = False
            else:
                cross_chart_offset_handles.append(
                    {
                        "from_chart": CHARTS[position - 1],
                        "to_chart": chart,
                        "shared_offset_handle": previous_parts[-1].get("right_offset"),
                    }
                )
        previous_parts = parts
        partition_summary[chart] = copy.deepcopy(parts)
        loop_none = "loop" in window and window.get("loop") is None
        quiet_local_zero = window.get("quiet_local_source_zero") is True or loop_none
        zero_increment = chart == "post_power" and quiet_local_zero
        zero_branches[chart] = bool(zero_increment and zero_source_admitted)
        if zero_increment and not zero_source_admitted:
            _gate(
                gates,
                "c2_original_power_zero_increment_receipt_missing",
                chart=chart,
                required="pinned original_power_zero_increment_source_admitted receipt",
            )
            ok = False
        if zero_increment and zero_source_admitted:
            if any(
                not _is_exact_zero_expression(graph, item, gates, f"C2 {chart} density zero")
                for item in density
            ) or any(
                not _is_exact_zero_expression(graph, item, gates, f"C2 {chart} contribution zero")
                for item in contribution
            ):
                _gate(gates, "c2_quiet_zero_branch_not_literal_zero", chart=chart)
                ok = False

        for order in ORDERS:
            try:
                root = contribution[order]
                terms = _flatten_add(graph, root, gates, f"C2 {chart} order {order} contribution")
                expected_integrals = [] if zero_increment and zero_source_admitted else list(terms or [])
                if zero_increment and zero_source_admitted:
                    if terms:
                        raise ValueError("quiet original-power contribution is not exact zero")
                else:
                    if len(expected_integrals) != len(parts) or any(graph[item].get("operation") != "definite_integral" for item in expected_integrals):
                        raise ValueError("contribution is not the exact connected partition sum")
                    for part, integral in zip(parts, expected_integrals):
                        integral_row = _row(graph, integral)
                        if integral_row.get("actual_current_chart") != chart or integral_row.get("exact_native_left") != part.get("exact_left") or integral_row.get("exact_native_right") != part.get("exact_right"):
                            raise ValueError("connected partition integral order or endpoints changed")
                    integrand_roots = {graph[item].get("integrand") for item in expected_integrals}
                    if len(integrand_roots) != 1:
                        raise ValueError("partition integrands are not connected")
                    for part, integral in zip(parts, expected_integrals):
                        irow = _row(graph, integral)
                        if irow.get("lower") != part.get("lower") or irow.get("upper") != part.get("upper"):
                            raise ValueError("partition domain endpoint changed")
                        good_kernel, kernel_info = _c2_kernel_binding(graph, integral, density[order], parts, order, chart, gates)
                        kernel_ok = kernel_ok and good_kernel
                        if not good_kernel:
                            ok = False
                        derivative_direct = derivative_direct and bool(kernel_info.get("offset_derivative_directly_bound_to_jacobian"))
                for key in (root, incoming[order], outgoing[order], density[order], memory):
                    _collect_graph_refs(graph, key, handles)
                for integral in expected_integrals:
                    _collect_graph_refs(graph, integral, handles)
                if chart not in kernel_summary:
                    kernel_summary[chart] = []
                kernel_summary[chart].append({"ordinary_order": order, "integrals": expected_integrals, "literal_zero": bool(zero_increment and zero_source_admitted)})
                if position == 0:
                    out_terms = _flatten_add(graph, outgoing[order], gates, f"C2 {chart} order {order} outgoing")
                    if outgoing[order] != contribution[order] or out_terms != expected_integrals:
                        raise ValueError("zero-inlet outgoing is not the exact contribution sum")
                else:
                    out_terms = _flatten_add(graph, outgoing[order], gates, f"C2 {chart} order {order} outgoing")
                    product_terms = [item for item in out_terms or [] if graph[item].get("operation") == "product"]
                    integral_terms = [item for item in out_terms or [] if graph[item].get("operation") == "definite_integral"]
                    if (
                        len(product_terms) != 1
                        or len(out_terms or []) != 1 + len(expected_integrals)
                        or Counter(out_terms or []) != Counter([product_terms[0], *expected_integrals])
                    ):
                        raise ValueError("outgoing is not memory plus the exact contribution sum")
                    product_args = graph[product_terms[0]].get("arguments")
                    if not isinstance(product_args, list) or Counter(product_args) != Counter([incoming[order], memory]):
                        raise ValueError("C2 predecessor-memory product changed")
                    _collect_graph_refs(graph, product_terms[0], handles)
                rows.append({"chart": chart, "field": M_COMPONENT, "ordinary_order": order, "incoming": incoming[order], "outgoing": outgoing[order], "contribution": contribution[order], "memory": memory, "density": density[order]})
            except (TypeError, ValueError, IndexError) as exc:
                _gate(gates, "c2_saved_transport_topology_changed", chart=chart, ordinary_order=order, error=str(exc))
                ok = False
                partition_ok = False
    graph_ok = _copy_rows(nodes, graph, handles, "C2 exact native transport", gates)
    ok = ok and graph_ok
    if not derivative_direct:
        _gate(gates, "c2_offset_derivative_equals_jacobian_unbound", required="saved source offset function identity or parameter derivative binding")
    return ok, {
        "rows": rows,
        "native_partitions": partition_summary,
        "cross_chart_shared_offset_handles": cross_chart_offset_handles,
        "cross_chart_comparison_uses_shared_offset_handles": True,
        "kernels": kernel_summary,
        "literal_zero_branches": zero_branches,
        "archive_graph_length": len(graph),
        "current_rows_authenticated": graph_ok,
        "exact_connected_partition_sums": partition_ok,
        "kernel_operand_bindings": kernel_ok,
        "offset_derivative_equals_jacobian": derivative_direct,
        "original_power_zero_increment_source_admitted": zero_source_admitted,
        "candidate_archive_flags": {key: c2_archive.get(key) for key in c2_archive if isinstance(key, str) and (key.endswith("installed") or key.startswith("candidate_"))},
        "candidate_archive_flags_ignored": True,
    }


def _is_exact_one_row(row: Mapping[str, Any]) -> bool:
    return row.get("operation") == "exact_rational" and row.get("numerator") == 1 and row.get("denominator") == 1


def _c5_memory_binding(
    graph: Sequence[Any],
    memory: int,
    coordinate: int,
    lower: int,
    measure: int,
    chart: str,
    gates: list[dict[str, Any]],
) -> tuple[bool, dict[str, Any]]:
    """Bind the C5 m-memory exp(-J*(coordinate-lower)) operands."""

    info: dict[str, Any] = {
        "memory": memory,
        "chart": chart,
        "coordinate": coordinate,
        "lower": lower,
        "physical_measure": measure,
        "rate": -1,
    }
    try:
        memory_row = _row(graph, memory)
        if memory_row.get("operation") != "analytic_unary" or memory_row.get("name") != "exp":
            raise ValueError("C5 memory is not an exponential")
        if _row(graph, coordinate).get("operation") != "bound_variable":
            raise ValueError("C5 memory coordinate is not a bound variable")
        exponent_handle = memory_row.get("argument")
        exponent = _row(graph, exponent_handle)
        factors = _flatten_product(graph, exponent_handle, gates, f"C5 {chart} memory exponent")
        if exponent.get("operation") != "product" or factors is None:
            raise ValueError("C5 memory exponent is not a product")
        rate_rows = [
            item
            for item in factors
            if _row(graph, item).get("operation") == "exact_rational"
            and _row(graph, item).get("numerator") == -1
            and _row(graph, item).get("denominator") == 1
        ]
        if len(rate_rows) != 1:
            raise ValueError("C5 memory rate operand changed")
        unit_measure = measure == 1 and _is_exact_one_row(_row(graph, measure))
        measure_count = factors.count(measure)
        if unit_measure:
            if measure_count not in (0, 1):
                raise ValueError("C5 unit memory Jacobian multiplicity changed")
        elif measure_count != 1:
            raise ValueError("C5 memory physical measure is not present exactly once")

        coordinate_difference: list[int] = []
        lower_is_zero = _is_zero_row(_row(graph, lower))
        for item in factors:
            if item in rate_rows or item == measure:
                continue
            if item == coordinate and lower_is_zero:
                coordinate_difference.append(item)
                continue
            candidate = _row(graph, item)
            candidate_args = candidate.get("arguments")
            if candidate.get("operation") != "sum" or not isinstance(candidate_args, list) or len(candidate_args) != 2:
                continue
            negative_terms = [
                child
                for child in candidate_args
                if type(child) is int
                and _row(graph, child).get("operation") == "negative"
                and _row(graph, child).get("argument") == lower
            ]
            positive_terms = [child for child in candidate_args if child == coordinate]
            if len(negative_terms) == 1 and len(positive_terms) == 1:
                coordinate_difference.append(item)
        if len(coordinate_difference) != 1:
            raise ValueError("C5 memory coordinate-minus-lower operand changed")
        expected_factor_count = 2 + (0 if unit_measure and measure_count == 0 else 1)
        if len(factors) != expected_factor_count:
            raise ValueError("C5 memory exponent has an extra factor")
        info.update(
            {
                "coordinate_difference": coordinate_difference[0],
                "physical_measure_factor_count": measure_count,
                "unit_jacobian_absorbed": unit_measure and measure_count == 0,
                "physical_measure_once_verified": True,
            }
        )
        return True, info
    except (TypeError, ValueError, IndexError) as exc:
        _gate(gates, "c5_memory_operands_changed", chart=chart, error=str(exc))
        info["error"] = str(exc)
        info["physical_measure_once_verified"] = False
        return False, info


def _c5_kernel_binding(
    graph: Sequence[Any],
    integral: int,
    density_handle: int,
    coordinate: int,
    measure: int,
    variable_handle: int,
    chart: str,
    order: int,
    gates: list[dict[str, Any]],
) -> tuple[bool, dict[str, Any]]:
    info: dict[str, Any] = {"integral": integral, "chart": chart, "ordinary_order": order, "physical_measure": measure}
    physical_measure_verified = False
    try:
        irow = _row(graph, integral)
        integrand_handle = irow.get("integrand")
        integrand = _row(graph, integrand_handle)
        outer_args = integrand.get("arguments")
        if integrand.get("operation") != "product" or not isinstance(outer_args, list) or len(outer_args) != 2:
            raise ValueError("C5 integral density/kernel product changed")
        density_rows = [item for item in outer_args if type(item) is int and _row(graph, item).get("operation") == "function_substitution" and _row(graph, item).get("expression") == density_handle]
        if len(density_rows) != 1:
            raise ValueError("C5 finite-N density substitution changed")
        body_candidates = [item for item in outer_args if item != density_rows[0]]
        if len(body_candidates) != 1:
            raise ValueError("C5 transport body changed")
        body = body_candidates[0]
        unit_measure = (
            measure == 1
            and _is_exact_one_row(_row(graph, measure))
        )
        body_factors = _flatten_product(graph, body, gates, f"C5 {chart} order {order} body")
        if body_factors is None:
            raise ValueError("C5 body product cannot be flattened")
        measure_count = body_factors.count(measure)
        kernels = [item for item in body_factors if _row(graph, item).get("operation") == "analytic_unary" and _row(graph, item).get("name") == "exp"]
        if unit_measure:
            if measure_count not in (0, 1) or len(kernels) != 1 or any(item not in (measure, kernels[0]) for item in body_factors):
                raise ValueError("unit C5 Jacobian absorption changed")
        elif measure_count != 1 or len(kernels) != 1 or len(body_factors) != 2:
            raise ValueError("C5 physical measure was not applied exactly once")
        kernel_handle = kernels[0]
        kernel = _row(graph, kernel_handle)
        exponent_handle = kernel.get("argument")
        exponent = _row(graph, exponent_handle)
        exponent_factors = _flatten_product(graph, exponent_handle, gates, f"C5 {chart} order {order} exponent")
        if exponent.get("operation") != "product" or exponent_factors is None:
            raise ValueError("C5 kernel exponent is not a product")
        rate_rows = [item for item in exponent_factors if _row(graph, item).get("operation") == "exact_rational" and _row(graph, item).get("numerator") == -1 and _row(graph, item).get("denominator") == 1]
        # A composite physical-Jacobian handle may itself be a sum.  Select
        # the coordinate-minus-integration-variable sum by its operands,
        # rather than treating every sum factor as the kernel difference.
        difference_rows: list[int] = []
        for item in exponent_factors:
            candidate = _row(graph, item)
            if candidate.get("operation") != "sum":
                continue
            candidate_args = candidate.get("arguments")
            if not isinstance(candidate_args, list) or len(candidate_args) != 2:
                continue
            candidate_negative = [
                child
                for child in candidate_args
                if type(child) is int
                and _row(graph, child).get("operation") == "negative"
                and _row(graph, child).get("argument") == variable_handle
            ]
            candidate_positive = [child for child in candidate_args if child == coordinate]
            if len(candidate_negative) == 1 and len(candidate_positive) == 1:
                difference_rows.append(item)
        if len(rate_rows) != 1 or len(difference_rows) != 1:
            raise ValueError("C5 rate-one kernel operands changed")
        if unit_measure:
            allowed = {rate_rows[0], difference_rows[0], measure}
            measure_factor_count = exponent_factors.count(measure)
            if (measure_factor_count not in (0, 1) or len(exponent_factors) != 2 + measure_factor_count
                    or any(item not in allowed for item in exponent_factors)
                    or exponent_factors.count(difference_rows[0]) != 1):
                raise ValueError("C5 unit-Jacobian exponent operands changed")
        else:
            if (exponent_factors.count(measure) != 1 or len(exponent_factors) != 3
                    or any(item not in (rate_rows[0], measure, difference_rows[0]) for item in exponent_factors)):
                raise ValueError("C5 physical Jacobian is not the exact kernel factor")
        difference = _row(graph, difference_rows[0])
        diff_args = difference.get("arguments")
        negative_terms = [
            item
            for item in diff_args
            if type(item) is int
            and _row(graph, item).get("operation") == "negative"
            and _row(graph, item).get("argument") == variable_handle
        ] if isinstance(diff_args, list) else []
        positive_terms = [item for item in diff_args if item == coordinate] if isinstance(diff_args, list) else []
        if not isinstance(diff_args, list) or len(diff_args) != 2 or len(negative_terms) != 1 or len(positive_terms) != 1:
            raise ValueError("C5 coordinate-minus-integration-variable operands changed")
        if _row(graph, coordinate).get("operation") != "bound_variable":
            raise ValueError("C5 coordinate is not a bound variable")
        physical_measure_verified = True
        info.update({
            "kernel": kernel_handle,
            "rate": rate_rows[0],
            "coordinate": coordinate,
            "integration_variable": variable_handle,
            "difference": difference_rows[0],
            "physical_measure_factor_count": measure_count,
            "unit_jacobian_absorbed": unit_measure and measure_count == 0,
            "kernel_jacobian_factor_count": exponent_factors.count(measure),
            "unit_measure_exact_one": unit_measure,
            "physical_measure_once_verified": physical_measure_verified,
        })
        return True, info
    except (TypeError, ValueError, IndexError) as exc:
        _gate(gates, "c5_kernel_operands_changed", chart=chart, ordinary_order=order, error=str(exc))
        info["error"] = str(exc)
        info["physical_measure_once_verified"] = physical_measure_verified
        return False, info


def _c5_witness(c5_archive: Mapping[str, Any], nodes: Sequence[Any], gates: list[dict[str, Any]]) -> tuple[bool, dict[str, Any]]:
    graph = c5_archive.get("exact_graph_nodes")
    bindings = c5_archive.get("prescribed_loop_C5_bindings")
    charts = bindings.get("charts") if isinstance(bindings, Mapping) else None
    endpoints = bindings.get("endpoints") if isinstance(bindings, Mapping) else None
    if not isinstance(graph, list) or not isinstance(charts, Mapping) or not isinstance(endpoints, Mapping):
        _gate(gates, "c5_archive_shape_changed")
        return False, {}
    expected_order = list(CHARTS) + [POSTRESTORE]
    if list(charts)[: len(expected_order)] != expected_order:
        _gate(gates, "c5_chart_order_changed")
        return False, {}
    bound_variables = {value.get("name"): index for index, value in enumerate(graph) if isinstance(value, Mapping) and value.get("operation") == "bound_variable" and isinstance(value.get("name"), str)}
    rows: list[dict[str, Any]] = []
    handles: set[int] = set()
    ok = True
    kernel_ok = True
    memory_ok = True
    measure_ok = True
    geometry_ok = True
    previous_endpoint: list[int] | None = None
    offsets: dict[str, Any] = {}
    memory_bindings: dict[str, Any] = {}
    measure_rows: dict[str, bool] = {}
    unit_measure_rows: dict[str, bool] = {}

    for position, chart in enumerate(expected_order):
        value = charts.get(chart)
        if not isinstance(value, Mapping):
            _gate(gates, "c5_chart_not_mapping", chart=chart)
            ok = False
            memory_ok = False
            continue
        if chart == POSTRESTORE:
            incoming = _as_int_list(_descriptor(value.get("incoming"), "C5 postrestore.incoming", gates).get(M_COMPONENT), 6, "C5 postrestore incoming m", gates)
            if incoming is None or incoming != previous_endpoint:
                _gate(gates, "c5_restoration_to_postrestore_endpoint_changed")
                ok = False
                memory_ok = False
            if incoming is not None:
                handles.update(incoming)
            continue
        incoming = _as_int_list(_descriptor(value.get("incoming"), f"C5 {chart}.incoming", gates).get(M_COMPONENT), 6, f"C5 {chart} incoming m", gates)
        history = _as_int_list(_descriptor(value.get("history"), f"C5 {chart}.history", gates).get(M_COMPONENT), 6, f"C5 {chart} history m", gates)
        contribution = _as_int_list(_descriptor(value.get("contributions"), f"C5 {chart}.contributions", gates).get(M_COMPONENT), 6, f"C5 {chart} contribution m", gates)
        endpoint = _as_int_list(_descriptor(_descriptor(endpoints, "C5.endpoints", gates).get(chart), f"C5 {chart}.endpoint", gates).get(M_COMPONENT), 6, f"C5 {chart} endpoint m", gates)
        memory = _descriptor(value.get("memory"), f"C5 {chart}.memory", gates).get(M_COMPONENT)
        coordinate = value.get("coordinate")
        physical_measure = value.get("physical_measure")
        domain = value.get("domain")
        offsets[chart] = {"coordinate": coordinate, "physical_measure": physical_measure, "global_offset": value.get("global_offset")}
        if not all(type(item) is int for item in (memory, coordinate, physical_measure, value.get("global_offset"))):
            _gate(gates, "c5_chart_geometry_changed", chart=chart)
            ok = False
            geometry_ok = False
        if incoming is None or history is None or contribution is None or endpoint is None or not isinstance(memory, int):
            ok = False
            memory_ok = False
            continue
        memory_good, memory_info = _c5_memory_binding(
            graph,
            memory,
            coordinate,
            domain[0] if isinstance(domain, list) and len(domain) == 2 and type(domain[0]) is int else -1,
            physical_measure,
            chart,
            gates,
        )
        memory_bindings[chart] = memory_info
        if not memory_good:
            ok = False
            memory_ok = False
        if position == 0:
            if not _all_zero_handles(graph, incoming, f"C5 {chart} incoming", gates):
                _gate(gates, "c5_zero_inlet_changed", chart=chart)
                ok = False
        elif incoming != previous_endpoint:
            _gate(gates, "c5_incoming_endpoint_recurrence_changed", chart=chart)
            ok = False
        previous_endpoint = endpoint
        if endpoint != _as_int_list(_descriptor(_descriptor(endpoints, "C5.endpoints", gates).get(chart), f"C5 {chart}.endpoint", gates).get(M_COMPONENT), 6, f"C5 endpoint {chart}", gates):
            ok = False

        finite_n = value.get("finite_N")
        density_vector = _descriptor(_descriptor(finite_n, f"C5 {chart}.finite_N", gates).get("five_signed_density_increments"), f"C5 {chart}.five_signed_density_increments", gates).get(M_COMPONENT)
        if not isinstance(density_vector, list) or len(density_vector) != 6 or any(type(item) is not int for item in density_vector):
            _gate(gates, "c5_density_source_map_changed", chart=chart)
            ok = False
            density_vector = [None] * 6
        variable_name = "true_loop_C5_transport_" + chart
        variable_handle = bound_variables.get(variable_name)
        if variable_handle is None:
            _gate(gates, "c5_integration_variable_missing", chart=chart, variable=variable_name)
            ok = False
            variable_handle = -1
        for order in ORDERS:
            try:
                endpoint_row = _row(graph, endpoint[order])
                if endpoint_row.get("operation") != "function_substitution" or endpoint_row.get("expression") != history[order] or endpoint_row.get("variable") != coordinate or not isinstance(domain, list) or len(domain) != 2 or endpoint_row.get("value") != domain[1] or endpoint_row.get("Z_independent_substitution") is not True:
                    raise ValueError("C5 endpoint substitution changed")
                history_row = _row(graph, history[order])
                if chart == CHARTS[0]:
                    if history[order] != contribution[order] or history_row.get("operation") != "definite_integral":
                        raise ValueError("C5 original-rminus zero-inlet history changed")
                else:
                    history_terms = _flatten_add(graph, history[order], gates, f"C5 {chart} order {order} history")
                    if history_terms is None:
                        raise ValueError("C5 history contribution term changed")
                    product_terms = [item for item in history_terms if graph[item].get("operation") == "product"]
                    if (
                        len(product_terms) != 1
                        or len(history_terms) != 2
                        or Counter(history_terms) != Counter([product_terms[0], contribution[order]])
                    ):
                        raise ValueError("C5 history is not predecessor memory plus contribution")
                    product_args = graph[product_terms[0]].get("arguments")
                    if not isinstance(product_args, list) or Counter(product_args) != Counter([incoming[order], memory]):
                        raise ValueError("C5 predecessor-memory product changed")
                integral = contribution[order]
                integral_row = _row(graph, integral)
                if integral_row.get("operation") != "definite_integral" or integral_row.get("original_chart") != chart or integral_row.get("original_rate") != "1":
                    raise ValueError("C5 transport integral metadata changed")
                if integral_row.get("physical_Jacobian_applied_once") is not True or integral_row.get("zero_correction_inlet_only_at_original_rminus") is not True or integral_row.get("predecessor_memory_never_reset") is not True or integral_row.get("endpoint_independent_of_Z") is not True:
                    raise ValueError("C5 transport guard changed")
                if not isinstance(domain, list) or len(domain) != 2 or integral_row.get("lower") != domain[0] or integral_row.get("upper") != coordinate or integral_row.get("variable") != variable_name:
                    raise ValueError("C5 native endpoint or integration variable changed")
                good_kernel, kernel_info = _c5_kernel_binding(graph, integral, density_vector[order], coordinate, physical_measure, variable_handle, chart, order, gates)
                kernel_ok = kernel_ok and good_kernel
                row_key = f"{chart}:{order}"
                row_measure_ok = bool(good_kernel and kernel_info.get("physical_measure_once_verified") is True)
                measure_rows[row_key] = row_measure_ok
                measure_ok = measure_ok and row_measure_ok
                if physical_measure == 1:
                    unit_row_ok = bool(
                        good_kernel
                        and kernel_info.get("physical_measure_once_verified") is True
                        and kernel_info.get("unit_measure_exact_one") is True
                        and (
                            kernel_info.get("unit_jacobian_absorbed") is True
                            or kernel_info.get("kernel_jacobian_factor_count") == 1
                        )
                    )
                    unit_measure_rows[row_key] = unit_row_ok
                if not good_kernel:
                    ok = False
                for key in (endpoint[order], history[order], contribution[order], incoming[order], memory, coordinate, physical_measure):
                    if type(key) is int and key >= 0:
                        _collect_graph_refs(graph, key, handles)
                rows.append({"chart": chart, "field": M_COMPONENT, "ordinary_order": order, "incoming": incoming[order], "history": history[order], "contribution": contribution[order], "endpoint": endpoint[order], "memory": memory, "density": density_vector[order]})
            except (TypeError, ValueError, IndexError) as exc:
                _gate(gates, "c5_saved_transport_topology_changed", chart=chart, ordinary_order=order, error=str(exc))
                ok = False
                memory_ok = False
    for chart in CHARTS:
        value = charts.get(chart)
        if not (isinstance(value, Mapping) and value.get("total_Z_phase_derivative_exact_zero") is True and value.get("physical_measure") is not None and value.get("global_offset") is not None):
            _gate(gates, "c5_shared_phase_or_offset_geometry_changed", chart=chart)
            geometry_ok = False
    graph_ok = _copy_rows(nodes, graph, handles, "C5 exact native transport", gates)
    ok = ok and graph_ok
    expected_unit_rows = {f"restoration:{order}" for order in ORDERS}
    return ok, {
        "rows": rows,
        "chart_geometry": offsets,
        "memory_operand_bindings": memory_bindings,
        "physical_measure_rows": measure_rows,
        "unit_measure_rows": unit_measure_rows,
        "archive_graph_length": len(graph),
        "current_rows_authenticated": graph_ok,
        "kernel_operand_bindings": kernel_ok,
        "nonreset_memory_recurrence": memory_ok,
        "physical_measure_once_semantic": bool(measure_ok and len(measure_rows) == len(CHARTS) * len(ORDERS)),
        "unit_jacobian_one_aggregate": bool(
            set(unit_measure_rows) == expected_unit_rows
            and all(unit_measure_rows.get(key) is True for key in expected_unit_rows)
        ),
        "geometry": geometry_ok,
        "candidate_archive_flags": {key: c5_archive.get(key) for key in c5_archive if isinstance(key, str) and key.endswith("installed")},
        "candidate_archive_flags_ignored": True,
    }


def _density_premise(source_records: Mapping[str, Any], c2: Mapping[str, Any], c5: Mapping[str, Any], gates: list[dict[str, Any]]) -> tuple[bool, dict[str, Any]]:
    """Validate an optional certificate, but never admit it as a theorem."""

    candidate = None
    if isinstance(source_records, Mapping):
        for key in ("same_density_premise", "density_premise", "density_certificate"):
            if isinstance(source_records.get(key), Mapping):
                candidate = source_records[key]
                break
    if not isinstance(candidate, Mapping):
        _gate(gates, "same_density_certificate_missing", required="independently pinned genuine C2=N*C5 density provider identity")
        return False, {"present": False, "admitted": False, "actual_relation_admitted": False}

    relation = candidate.get("relation")
    c2_density = candidate.get("rho_C2")
    c5_density = candidate.get("rho_C5")
    source_graph_ok = candidate.get("source_graph_digest") == EXPECTED_NODES_DIGEST and candidate.get("source_graph_node_count") == EXPECTED_NODE_COUNT
    source_hashes = candidate.get("source_record_hashes")
    hashes_ok = isinstance(source_hashes, Mapping) and all(source_hashes.get(name) == digest for name, digest in RECORD_HASHES.items())
    cert_digest = candidate.get("certificate_sha256", candidate.get("certificate_digest"))
    digest_body = copy.deepcopy(dict(candidate))
    digest_body.pop("certificate_sha256", None)
    digest_body.pop("certificate_digest", None)
    cert_digest_ok = isinstance(cert_digest, str) and cert_digest == canonical_digest(digest_body)
    fields_ok = candidate.get("fields") in (["m"], "m")
    orders_ok = candidate.get("ordinary_orders") == list(ORDERS)
    bindings = candidate.get("per_chart_bindings")
    expected_c2 = {(row.get("chart"), row.get("ordinary_order")): row for row in c2.get("rows", []) if isinstance(row, Mapping)}
    expected_c5 = {(row.get("chart"), row.get("ordinary_order")): row for row in c5.get("rows", []) if isinstance(row, Mapping)}
    bindings_ok = isinstance(bindings, Mapping) and set(bindings) == set(CHARTS)
    if bindings_ok:
        for chart in CHARTS:
            binding = bindings.get(chart)
            if not isinstance(binding, Mapping) or binding.get("field") != M_COMPONENT or binding.get("ordinary_orders") != list(ORDERS) or binding.get("relation") != "C2=N*C5" or binding.get("rho_C2") != "B" or binding.get("rho_C5") != "B_over_Pstar/N":
                bindings_ok = False
                break
            order_bindings = binding.get("orders") if isinstance(binding.get("orders"), Mapping) else binding
            for order in ORDERS:
                left, right = expected_c2.get((chart, order), {}), expected_c5.get((chart, order), {})
                expected = {"c2_incoming": left.get("incoming"), "c2_outgoing": left.get("outgoing"), "c5_incoming": right.get("incoming"), "c5_history": right.get("history"), "c5_endpoint": right.get("endpoint")}
                if order_bindings.get(str(order), order_bindings.get(order)) != expected:
                    bindings_ok = False
                    break
            if not bindings_ok:
                break
    structural_complete = bool(relation == "C2=N*C5" and c2_density == "B" and c5_density == "B_over_Pstar/N" and fields_ok and orders_ok and source_graph_ok and hashes_ok and cert_digest_ok and bindings_ok)
    _gate(gates, "independent_density_provider_receipt_not_pinned", self_described_certificate_complete=structural_complete)
    _gate(gates, "same_density_relation_not_admitted", required="independently pinned genuine density provider certificate")
    return False, {
        "present": True,
        "admitted": False,
        "actual_relation_admitted": False,
        "relation": relation,
        "rho_C2": c2_density,
        "rho_C5": c5_density,
        "fields": copy.deepcopy(candidate.get("fields")),
        "ordinary_orders": copy.deepcopy(candidate.get("ordinary_orders")),
        "source_graph_digest": candidate.get("source_graph_digest"),
        "source_graph_node_count": candidate.get("source_graph_node_count"),
        "source_record_hashes": copy.deepcopy(source_hashes),
        "certificate_sha256": cert_digest,
        "certificate_digest_matches": cert_digest_ok,
        "full_per_chart_bindings": bindings_ok,
        "independent_provider_receipt": candidate.get("independent_provider_receipt"),
    }


def _cross_archive_kernel_summary(c2: Mapping[str, Any], c5: Mapping[str, Any], source_witness: Mapping[str, Any], gates: list[dict[str, Any]]) -> dict[str, Any]:
    """Keep the missing source offset identity explicit and fail closed."""

    formula_tokens = bool(source_witness.get("static_direct_api_arity_checked") and source_witness.get("source_chain_structural"))
    _gate(gates, "source_offset_function_identity_unbound", required="an independently pinned source identity for d(offset)/dcoordinate = J")
    _gate(gates, "cross_archive_kernel_equality_unproved", required="C2 exp[-(right_offset-offset(t))] equals C5 exp[-rate*J*(coordinate-t)]")
    return {
        "c2_kernel_formula": "exp[-(right_offset-offset(t))]",
        "c5_kernel_formula": "exp[-rate*J*(coordinate-t)]",
        "formula_tokens_authenticated": formula_tokens,
        "actual_offset_derivative_equals_J": False,
        "cross_archive_kernel_equality": False,
        "source_offset_function_identity_available": False,
        "structural_directed_checks_independent": bool(c2.get("kernel_operand_bindings") and c5.get("kernel_operand_bindings")),
        "proved": False,
    }


def _relation_rows(c2: Mapping[str, Any], c5: Mapping[str, Any]) -> list[dict[str, Any]]:
    c2_by = {(row.get("chart"), row.get("ordinary_order")): row for row in c2.get("rows", []) if isinstance(row, Mapping)}
    c5_by = {(row.get("chart"), row.get("ordinary_order")): row for row in c5.get("rows", []) if isinstance(row, Mapping)}
    rows: list[dict[str, Any]] = []
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
    """Return a fail-closed exact native transport source certificate."""

    # The root graph checkpoint is deliberately first: no archive decoding or
    # source-record inspection is reached for an unexpected graph.
    _require(expected_nodes_digest == EXPECTED_NODES_DIGEST, "unexpected source graph digest assertion")
    _require(len(nodes) == EXPECTED_NODE_COUNT and canonical_digest(nodes) == EXPECTED_NODES_DIGEST, "saved source graph changed")

    gates: list[dict[str, Any]] = []
    source_ok, source_witness, parsed = _source_witnesses(source_records, gates)
    c2_ok = False
    c5_ok = False
    c2_witness: dict[str, Any] = {}
    c5_witness: dict[str, Any] = {}
    if parsed:
        c2_ok, c2_witness = _c2_witness(
            parsed.get(C2_DENSITY_GZ, {}),
            nodes,
            gates,
            zero_source_admitted=bool(source_witness.get("original_power_zero_increment_source_admitted")),
        )
        c5_ok, c5_witness = _c5_witness(parsed.get(C5_GZ, {}), nodes, gates)
    density_ok, density_witness = _density_premise(source_records, c2_witness, c5_witness, gates)
    cross_kernel = _cross_archive_kernel_summary(c2_witness, c5_witness, source_witness, gates)

    structural_ok = bool(source_ok and c2_ok and c5_ok)
    transport_shape_ok = bool(
        source_ok
        and
        c2_ok
        and c5_ok
        and c2_witness.get("exact_connected_partition_sums")
        and c2_witness.get("kernel_operand_bindings")
        and c5_witness.get("kernel_operand_bindings")
        and c5_witness.get("physical_measure_once_semantic")
        and c5_witness.get("unit_jacobian_one_aggregate")
        and c5_witness.get("nonreset_memory_recurrence")
    )
    theorem = {
        "premises": {
            "source_records_authenticated": source_ok,
            "exact_connected_C2_partition_sums": bool(c2_witness.get("exact_connected_partition_sums")),
            "C2_kernel_operands_bound": bool(c2_witness.get("kernel_operand_bindings")),
            "C5_kernel_operands_bound": bool(c5_witness.get("kernel_operand_bindings")),
            "C5_unit_Jacobian_one_absorbed": bool(c5_witness.get("unit_jacobian_one_aggregate")),
            "physical_measure_once_semantically": bool(c5_witness.get("physical_measure_once_semantic")),
            "zero_inlet_only_at_original_rminus": bool(c5_witness.get("nonreset_memory_recurrence")),
            "predecessor_memory_never_reset": bool(c5_witness.get("nonreset_memory_recurrence")),
            "actual_offset_derivative_equals_J": False,
            "cross_archive_kernel_equality": False,
            "same_density_C2_equals_N_C5": bool(density_ok),
        },
        "conclusion": "for m and ordinary Z orders 0,1,2, corresponding native chart histories satisfy H_C2=N*H_C5 through restoration",
        "C2_kernel": "exp[-(right_offset-offset(t))]",
        "C5_kernel": "exp[-rate*J*(coordinate-t)]",
        "shared_phase_offset_rule": "d(offset)/dcoordinate=J; phase is held fixed in ordinary Z rows",
        "native_partition_endpoint_additivity_required": True,
        "proved": False,
        "actual_relation_admitted": False,
        "numeric_values_computed": False,
    }

    return {
        "schema": SCHEMA,
        "source_only": True,
        "nodes_digest": EXPECTED_NODES_DIGEST,
        "nodes_count": EXPECTED_NODE_COUNT,
        "scope": {
            "field": M_COMPONENT,
            "ordinary_orders": list(ORDERS),
            "charts_through_restoration": list(CHARTS),
            "postrestore_endpoint_checked_on_C5": bool(c5_ok),
            "full_history_C2_claimed": False,
            "Qy1_claimed": False,
        },
        "passed": False,
        "structural_source_rows_checked": structural_ok,
        "native_transport_definitions_checked": transport_shape_ok,
        "transport_calculus_checked": False,
        "conditional_transport_calculus_supported": False,
        "conditional_transport_theorem": theorem,
        "conditional_transport_calculus": theorem,
        "density_premise": density_witness,
        "cross_archive_kernel": cross_kernel,
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
        "strongest_proved_scope": {
            "saved_json_and_source_ast_hashes": source_ok,
            "saved_graph_rows_authenticated": structural_ok,
            "C2_exact_connected_partition_integrals": bool(c2_witness.get("exact_connected_partition_sums")),
            "C2_kernel_exp_right_offset_minus_offset_t_operands": bool(c2_witness.get("kernel_operand_bindings")),
            "C5_kernel_exp_rate_J_coordinate_minus_t_operands": bool(c5_witness.get("kernel_operand_bindings")),
            "C5_unit_Jacobian_one_absorption": bool(c5_witness.get("unit_jacobian_one_aggregate")),
            "C5_zero_inlet_and_nonreset_memory": bool(c5_witness.get("nonreset_memory_recurrence")),
            "same_density_identity": False,
            "cross_archive_offset_identity": False,
            "history_relation_H_C2_equals_N_H_C5": False,
        },
        "witnesses": {
            "producer_sources": source_witness,
            "C2_saved_transport": c2_witness,
            "C5_saved_transport": c5_witness,
            "source_chain_structural": source_ok,
            "C2_graph_structural": c2_ok,
            "C5_graph_structural": c5_ok,
        },
        "unresolved_gates": gates,
        "numeric_values_computed": False,
        "old_operators_replayed": False,
        "operators_reexecuted": False,
        "full_goal_complete": False,
        "full_Qy1_or_full_goal_claimed": False,
        "source_record_hashes": dict(RECORD_HASHES),
    }


__all__ = [
    "SCHEMA",
    "RECORD_HASHES",
    "EXPECTED_NODES_DIGEST",
    "EXPECTED_NODE_COUNT",
    "canonical_digest",
    "derive",
]
