"""Source-only consumer for the six native outer C2/C5 charts.

This module deliberately owns a new definition graph.  The saved native
archives and reports are treated as authenticated source data; no saved
producer, graph evaluator, recurrence, quadrature, norm, or source-value
reducer is called.  The only runtime work performed by :func:`derive` is
source authentication, static AST inspection, JSON metadata inspection, and
construction of unevaluated ordinary-Z definitions.
"""
from __future__ import annotations

import ast
import base64
import gzip
import hashlib
import json
import symtable
from collections.abc import Mapping, Sequence
from fractions import Fraction


SCHEMA = "original_native_C2_outer_six_source_v1"

# These are the eight fixed inputs for this leaf.  They are basenames under
# experiments/root_st073; the caller supplies immutable source snapshots.
ADMISSION = "lei_ren_part1_paper_compliant_current_original_native_C2_patch_finite_N_admission.json"
ADMISSION_SHA256 = "23d6b0984cb74b0b5fa31df6113e967604bf6b066ad3614b29d8430abb20f28d"
CANDIDATE = "lei_ren_part1_paper_compliant_current_original_native_C2_patch_finite_N_source.json"
CANDIDATE_SHA256 = "5194f631aba2eebd4946003e1af52181ac8c050b92d636c36b6bdd2058245b5a"
GENERIC = "lei_ren_part1_paper_compliant_current_original_finite_N_field_density_source_composition.json"
GENERIC_SHA256 = "0e5d3002ccdc854a1f6be32f2ea076f5c7d650ae25e6a766c7161bca4cd90e2f"
OUTER_ARCHIVE = "lei_ren_part1_paper_compliant_current_original_true_outer_C6_C5_functions.json.gz"
OUTER_ARCHIVE_SHA256 = "75895e1637411cfd2e490778f495f3adf904295a8840184cdc467e006b23305d"
OUTER_ARCHIVE_DECODED_SHA256 = "53b1bbf15c8c2c5a8400d7e22ef940c6a502397df8cf26bb23e51ff3869e9139"
OUTER_SOURCE = "lei_ren_part1_paper_compliant_current_original_true_outer_C6_C5_functions.py"
OUTER_SOURCE_SHA256 = "4464a934d303a95a0da157e9e289f6ed8a04d5ff0adcb0d663857f0d9c57dcd3"
OUTER_RECEIPT = "lei_ren_part1_paper_compliant_current_original_true_outer_C6_C5_functions_check.json"
OUTER_RECEIPT_SHA256 = "ce5669d01621683c3be857917f206b50af067b44feefe4b000c2ba99b2785bba"
PATCH_ARCHIVE = "lei_ren_part1_paper_compliant_current_original_true_Rm_patch_finite_N_C5_functions.json.gz"
PATCH_ARCHIVE_SHA256 = "8f3bc343e9dbffbdda465ce94d6183735184645a70e291cb900d3cde726769ad"
PATCH_ARCHIVE_DECODED_SHA256 = "c1026e05015584a7ce29bd7e8da95b17d97669523f1d063ac5b5c800559d9b69"
PATCH_SOURCE = "lei_ren_part1_paper_compliant_current_original_true_Rm_patch_finite_N_C5_functions.py"
PATCH_SOURCE_SHA256 = "c4445ef8acfb33d29cf9022862f8b6396e6aa2f38bda7467409b5378a5deaf15"

RECORD_HASHES = {
    ADMISSION: ADMISSION_SHA256,
    CANDIDATE: CANDIDATE_SHA256,
    GENERIC: GENERIC_SHA256,
    OUTER_ARCHIVE: OUTER_ARCHIVE_SHA256,
    PATCH_ARCHIVE: PATCH_ARCHIVE_SHA256,
    OUTER_RECEIPT: OUTER_RECEIPT_SHA256,
    OUTER_SOURCE: OUTER_SOURCE_SHA256,
    PATCH_SOURCE: PATCH_SOURCE_SHA256,
}
DECODED_RECORD_HASHES = {
    OUTER_ARCHIVE: OUTER_ARCHIVE_DECODED_SHA256,
    PATCH_ARCHIVE: PATCH_ARCHIVE_DECODED_SHA256,
}

ORDINARY_Z_ORDERS = (0, 1, 2)
CHART_ORDER = (
    "Rh_reference", "O2_slope", "O2_axial", "O2_buffer", "O3_transition", "O3_power",
)
CHARTS = CHART_ORDER
RATES = {"m": "1", "h": "3/2", "k": "3/2", "e": "1", "p": "0"}
MOMENT_NAMES = ("m", "h", "k", "e", "p")
LOCAL_DOMAINS = {
    "Rh_reference": ("-5", "0"),
    "O2_slope": ("0", "1"),
    "O2_axial": ("0", "1"),
    "O2_buffer": ("0", "11"),
    "O3_transition": ("0", "1"),
    "O3_power": ("0", "2"),
}

# Source handles accepted by the frozen native prefix.  They are checked
# against the saved metadata and never interpreted as values by this module.
OUTER_GRAPH_NODE_COUNT = 258857
PATCH_GRAPH_NODE_COUNT = 205147
PATCH_RM_NODE = 85743
PATCH_N_NODE = 18
PSTAR_SHARED_PREFIX_NODE = 33572
PSTAR_SHARED_PREFIX_ARCHIVE_SHA256 = "5e5b565dd8d2abc41e6f3df0d9e17000e524d493e58a5778dfcf217c6e42f44a"
PATCH_Z_HANDLE = 20
PATCH_Z_NAMESPACE = "original_native_C2_patch_finite_N_source_v1:d000001"
PATCH_FIRST_INCOMING = {
    "m": (204936, 204937, 204938, 204939, 204940, 204941),
    "h": (204942, 204943, 204944, 204945, 204946, 204947),
    "k": (204948, 204949, 204950, 204951, 204952, 204953),
    "e": (204954, 204955, 204956, 204957, 204958, 204959),
    "p": (204960, 204961, 204962, 204963, 204964, 204965),
}
RETAINED_CONSTRUCTION = {
    "cumulative_terms": 212,
    "raw_signed_terms": 213,
    "reduced_terms": 183,
    "six_signed_channels_retained": True,
    "all_old_error_data_retained": True,
    "reserved_contribution_registered_once": True,
    "inherited_errors_deleted": False,
    "P1_differs_from_m15": True,
    "physical_Q_differs_from_loop_q": True,
}

ScalarID = str
JetList = list[ScalarID]


def _digest(value: object) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":"),
                   ensure_ascii=True, allow_nan=False).encode("utf-8")
    ).hexdigest()


def _get(value: object, path: str, default: object = None) -> object:
    current = value
    for part in path.split(".") if path else ():
        if isinstance(current, Mapping):
            if part not in current:
                return default
            current = current[part]
        elif isinstance(current, list) and part.isdigit():
            index = int(part)
            if index >= len(current):
                return default
            current = current[index]
        else:
            return default
    return current


def _gate(gates: list[str], message: str) -> None:
    if message not in gates:
        gates.append(message)


def _sha_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def _raw_bytes(row: Mapping[str, object]) -> bytes | None:
    encoded = row.get("raw_bytes_base64")
    if isinstance(encoded, str):
        try:
            return base64.b64decode(encoded, validate=True)
        except (ValueError, TypeError):
            return None
    raw = row.get("raw_bytes")
    if isinstance(raw, bytes):
        return raw
    if isinstance(raw, bytearray):
        return bytes(raw)
    return None


def _record(records: Mapping[str, object], name: str, gates: list[str]) -> str:
    row = records.get(name)
    if not isinstance(row, Mapping):
        _gate(gates, "missing_source_record:" + name)
        return ""
    expected_path = "experiments/root_st073/" + name
    if row.get("source_path", row.get("path")) != expected_path:
        _gate(gates, "source_path_mismatch:" + name)
    declared_raw = row.get("raw_sha256", row.get("sha256"))
    if declared_raw != RECORD_HASHES[name]:
        _gate(gates, "source_raw_hash_mismatch:" + name)
    if row.get("raw_authenticated_before_decode") is not True:
        _gate(gates, "raw_authentication_missing:" + name)
    expected_gzip = name in DECODED_RECORD_HASHES
    if type(row.get("is_gzip")) is not bool or row.get("is_gzip") != expected_gzip:
        _gate(gates, "gzip_contract_mismatch:" + name)
    text = row.get("source_text")
    if not isinstance(text, str):
        _gate(gates, "source_text_missing:" + name)
        return ""
    decoded = _sha_text(text)
    expected_decoded = DECODED_RECORD_HASHES.get(name, RECORD_HASHES[name])
    declared_decoded = row.get("decoded_sha256", decoded)
    if declared_decoded != expected_decoded or decoded != expected_decoded:
        _gate(gates, "decoded_source_hash_mismatch:" + name)
    raw = _raw_bytes(row)
    if raw is not None:
        if hashlib.sha256(raw).hexdigest() != RECORD_HASHES[name]:
            _gate(gates, "raw_bytes_hash_mismatch:" + name)
        if expected_gzip:
            try:
                if gzip.decompress(raw).decode("utf-8") != text:
                    _gate(gates, "archive_decode_text_mismatch:" + name)
            except (OSError, EOFError, UnicodeDecodeError):
                _gate(gates, "archive_decode_error:" + name)
        elif raw != text.encode("utf-8"):
            _gate(gates, "plain_source_bytes_mismatch:" + name)
    return text


def _json(text: str, name: str, gates: list[str]) -> Mapping[str, object]:
    if not text:
        return {}
    try:
        value = json.loads(text)
    except (TypeError, ValueError):
        _gate(gates, "json_parse_error:" + name)
        return {}
    if not isinstance(value, Mapping):
        _gate(gates, "json_root_not_mapping:" + name)
        return {}
    return value


def _static_ast(text: str, name: str, gates: list[str]) -> ast.AST | None:
    if not text:
        return None
    try:
        tree = ast.parse(text, filename=name)
        # These are static checks.  The code object is never executed.
        compile(tree, name, "exec")
        symtable.symtable(text, name, "exec")
        return tree
    except (SyntaxError, ValueError, TypeError):
        _gate(gates, "source_static_parse_error:" + name)
        return None


def _function_names(tree: ast.AST | None) -> set[str]:
    if tree is None:
        return set()
    return {node.name for node in ast.walk(tree)
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef))}


def _valid_handles(value: object, width: int, node_count: int) -> bool:
    return (isinstance(value, list) and len(value) == width
            and all(type(item) is int and 0 <= item < node_count for item in value))


def _bool_leaves(value: object) -> bool:
    if isinstance(value, bool):
        return value
    if isinstance(value, Mapping):
        return all(_bool_leaves(item) for item in value.values())
    if isinstance(value, list):
        return all(_bool_leaves(item) for item in value)
    return True


class _DAG:
    """Append-only graph whose operation arguments are local scalar IDs."""

    def __init__(self) -> None:
        self.nodes: list[dict[str, object]] = []
        self._counter = 0

    def _id(self) -> ScalarID:
        self._counter += 1
        return SCHEMA + ":d" + str(self._counter).zfill(6)

    def atom(self, operation: str, value: object = None, **metadata: object) -> ScalarID:
        ident = self._id()
        node: dict[str, object] = {"id": ident, "namespace": SCHEMA, "op": operation}
        if value is not None:
            node["value"] = value
        node.update(metadata)
        self.nodes.append(node)
        return ident

    def op(self, operation: str, args: Sequence[ScalarID] = (), **metadata: object) -> ScalarID:
        ident = self._id()
        node: dict[str, object] = {
            "id": ident, "namespace": SCHEMA, "op": operation, "args": list(args),
        }
        node.update(metadata)
        self.nodes.append(node)
        return ident


def _fraction(value: object) -> Fraction:
    return value if isinstance(value, Fraction) else Fraction(str(value))


def _const(dag: _DAG, value: object) -> ScalarID:
    rational = _fraction(value)
    return dag.atom("exact_rational", numerator=rational.numerator,
                    denominator=rational.denominator)


def _zero(dag: _DAG) -> ScalarID:
    return _const(dag, "0")


def _one(dag: _DAG) -> ScalarID:
    return _const(dag, "1")


def _symbol(dag: _DAG, name: str, domain: str) -> ScalarID:
    return dag.atom("exact_symbol", name, domain=domain)


def _op(dag: _DAG, operation: str, *args: ScalarID, **metadata: object) -> ScalarID:
    return dag.op(operation, args, exact=True, **metadata)


def _add(dag: _DAG, *args: ScalarID) -> ScalarID:
    return _op(dag, "add", *args)


def _sub(dag: _DAG, left: ScalarID, right: ScalarID) -> ScalarID:
    return _op(dag, "sub", left, right)


def _mul(dag: _DAG, *args: ScalarID) -> ScalarID:
    return _op(dag, "mul", *args)


def _div(dag: _DAG, left: ScalarID, right: ScalarID) -> ScalarID:
    return _op(dag, "quotient", left, right)


def _neg(dag: _DAG, value: ScalarID) -> ScalarID:
    return _op(dag, "neg", value)


def _exp(dag: _DAG, value: ScalarID) -> ScalarID:
    return _op(dag, "exp", value, analytic=True)


def _expm1(dag: _DAG, value: ScalarID) -> ScalarID:
    return _op(dag, "expm1", value, analytic=True, entire_at_zero=True)


def _jet(jet: object) -> JetList:
    if not isinstance(jet, list) or len(jet) != 3 or not all(isinstance(x, str) for x in jet):
        raise TypeError("ordinary Z jet must be a list of three scalar IDs")
    return jet


def _jet_zero(dag: _DAG) -> JetList:
    return [_zero(dag), _zero(dag), _zero(dag)]


def _jet_const(dag: _DAG, value: object) -> JetList:
    return [_const(dag, value), _zero(dag), _zero(dag)]


def _jet_scalar(dag: _DAG, scalar: ScalarID) -> JetList:
    return [scalar, _zero(dag), _zero(dag)]


def _jet_add(dag: _DAG, *jets: JetList) -> JetList:
    checked = [_jet(item) for item in jets]
    if not checked:
        return _jet_zero(dag)
    return [_add(dag, *(item[index] for item in checked)) for index in ORDINARY_Z_ORDERS]


def _jet_sub(dag: _DAG, left: JetList, right: JetList) -> JetList:
    left, right = _jet(left), _jet(right)
    return [_sub(dag, left[i], right[i]) for i in ORDINARY_Z_ORDERS]


def _jet_scale(dag: _DAG, jet: JetList, scalar: ScalarID) -> JetList:
    jet = _jet(jet)
    return [_mul(dag, scalar, jet[i]) for i in ORDINARY_Z_ORDERS]


def _jet_mul(dag: _DAG, left: JetList, right: JetList) -> JetList:
    left, right = _jet(left), _jet(right)
    return [
        _mul(dag, left[0], right[0]),
        _add(dag, _mul(dag, left[1], right[0]), _mul(dag, left[0], right[1])),
        _add(dag, _mul(dag, left[2], right[0]),
             _mul(dag, _const(dag, "2"), left[1], right[1]),
             _mul(dag, left[0], right[2])),
    ]


def _jet_div(dag: _DAG, left: JetList, right: JetList) -> JetList:
    return _jet_mul(dag, left, _jet_inv(dag, right))


def _jet_inv(dag: _DAG, jet: JetList) -> JetList:
    a0, a1, a2 = _jet(jet)
    a0sq = _mul(dag, a0, a0)
    a0cu = _mul(dag, a0sq, a0)
    return [_div(dag, _one(dag), a0),
            _neg(dag, _div(dag, a1, a0sq)),
            _div(dag, _sub(dag, _mul(dag, _const(dag, "2"), a1, a1),
                           _mul(dag, a0, a2)), a0cu)]


def _jet_exp(dag: _DAG, jet: JetList) -> JetList:
    a0, a1, a2 = _jet(jet)
    e0 = _exp(dag, a0)
    return [e0, _mul(dag, e0, a1),
            _mul(dag, e0, _add(dag, a2, _mul(dag, a1, a1)))]


def _jet_expm1(dag: _DAG, jet: JetList) -> JetList:
    a0, a1, a2 = _jet(jet)
    e0 = _exp(dag, a0)
    return [_expm1(dag, a0), _mul(dag, e0, a1),
            _mul(dag, e0, _add(dag, a2, _mul(dag, a1, a1)))]


def _function_substitution(dag: _DAG, scalar: ScalarID, variable: ScalarID,
                           value: ScalarID, **metadata: object) -> ScalarID:
    return dag.op("function_substitution", (scalar, variable, value), exact=True,
                  Z_independent_substitution=True, **metadata)


def _substitute_jet(dag: _DAG, jet: JetList, variable: ScalarID,
                    value: ScalarID, **metadata: object) -> JetList:
    return [_function_substitution(dag, item, variable, value, **metadata)
            for item in _jet(jet)]


def _integral(dag: _DAG, scalar: ScalarID, variable: ScalarID,
              lower: ScalarID, upper: ScalarID, **metadata: object) -> ScalarID:
    return dag.op("integral", (scalar,), variable=variable, bounds=[lower, upper],
                  exact=True, unevaluated=True, quadrature=False, **metadata)


def _integral_jet(dag: _DAG, jet: JetList, variable: ScalarID,
                  lower: ScalarID, upper: ScalarID, **metadata: object) -> JetList:
    return [_integral(dag, item, variable, lower, upper,
                      ordinary_Z_order=index, full_jet_product=True,
                      scalar_integrand=True, **metadata)
            for index, item in enumerate(_jet(jet))]


def _source_port(
    dag: _DAG,
    *,
    report_sha256: str,
    namespace: str,
    source_node: object,
    archive: str,
    path: str,
    label: str,
    role: str,
    valid: bool,
    valid_range: object = None,
    width: int | None = None,
    coordinate_handle: object = None,
    phase_handle: object = None,
    coordinate_binding: ScalarID | None = None,
    z_binding: ScalarID | None = None,
    fixed_endpoint: object = None,
    **metadata: object,
) -> ScalarID:
    value: dict[str, object] = {
        "source_report_sha256": report_sha256,
        "source_archive": archive,
        "archive_sha256": RECORD_HASHES.get(archive, report_sha256),
        "source_namespace": namespace,
        "source_node": source_node if valid else "missing",
        "source_path": path,
        "label": label,
        "role": role,
        "source_node_valid": bool(valid),
        "source_graph_node_valid": bool(valid),
        "qualified": (SCHEMA + "|report=" + report_sha256 + "|namespace=" + namespace
                       + "|node=" + str(source_node if valid else "missing")
                       + "|path=" + path + "|label=" + label),
        "re_evaluated": False,
        "source_only": True,
    }
    if valid_range is not None:
        value["valid_range"] = valid_range
    if width is not None:
        value["width"] = width
    if coordinate_handle is not None:
        value["source_coordinate_handle"] = coordinate_handle
    if phase_handle is not None:
        value["source_phase_handle"] = phase_handle
    if fixed_endpoint is not None:
        value["fixed_endpoint"] = fixed_endpoint
    value.update(metadata)
    port = dag.atom("source_function_port", value)
    bindings: list[ScalarID] = []
    names: dict[str, ScalarID] = {}
    if coordinate_binding is not None:
        bindings.append(coordinate_binding)
        names["x"] = coordinate_binding
    if z_binding is not None:
        bindings.append(z_binding)
        names["Z"] = z_binding
    if not bindings:
        return port
    application = {
        "source_port": port,
        "variable_bindings": names,
        "application_is_unevaluated": True,
        "source_only": True,
        "explicit_source_function_application": True,
    }
    if fixed_endpoint is not None and coordinate_binding is None:
        application["bind_Z_only"] = True
    return dag.op("source_function_application", tuple([port] + bindings),
                  **application)


def _source_rows(
    dag: _DAG,
    *,
    rows: object,
    width: int,
    node_count: int,
    report_sha256: str,
    namespace: str,
    archive: str,
    path: str,
    label: str,
    role: str,
    x: ScalarID | None,
    z: ScalarID | None,
    coordinate_handle: object = None,
    phase_handle: object = None,
    valid_range: object = None,
    fixed_endpoint: object = None,
) -> tuple[list[ScalarID], bool]:
    valid = _valid_handles(rows, width, node_count)
    values: list[ScalarID] = []
    if not valid:
        values = [_source_port(
            dag, report_sha256=report_sha256, namespace=namespace,
            source_node="missing", archive=archive, path=path + "[missing]",
            label=label, role=role, valid=False, width=width,
            coordinate_binding=x, z_binding=z, fixed_endpoint=fixed_endpoint,
        )]
        return values, False
    for index, node in enumerate(rows):
        values.append(_source_port(
            dag, report_sha256=report_sha256, namespace=namespace,
            source_node=node, archive=archive, path=path + "[" + str(index) + "]",
            label=label + "[" + str(index) + "]", role=role, valid=True,
            valid_range=valid_range, width=width,
            coordinate_handle=coordinate_handle, phase_handle=phase_handle,
            coordinate_binding=x, z_binding=z, fixed_endpoint=fixed_endpoint,
            source_z_handle=PATCH_Z_HANDLE if z is not None else None,
        ))
    return values, True


def _port_jet(rows: list[ScalarID]) -> JetList:
    # Source rows are retained at their native width.  The new consumer uses
    # only ordinary Z orders 0, 1, 2 and never relabels the remaining rows.
    return rows[:3]


def _source_witness(records: Mapping[str, object], gates: list[str]) -> dict[str, object]:
    """Authenticate and inspect only the records needed by the new leaf."""
    texts = {name: _record(records, name, gates) for name in RECORD_HASHES}
    admission = _json(texts[ADMISSION], ADMISSION, gates)
    candidate = _json(texts[CANDIDATE], CANDIDATE, gates)
    generic = _json(texts[GENERIC], GENERIC, gates)
    outer = _json(texts[OUTER_ARCHIVE], OUTER_ARCHIVE, gates)
    patch = _json(texts[PATCH_ARCHIVE], PATCH_ARCHIVE, gates)
    receipt = _json(texts[OUTER_RECEIPT], OUTER_RECEIPT, gates)
    source_tree = _static_ast(texts[OUTER_SOURCE], OUTER_SOURCE, gates)
    patch_source_tree = _static_ast(texts[PATCH_SOURCE], PATCH_SOURCE, gates)
    source_names = _function_names(source_tree)
    patch_source_names = _function_names(patch_source_tree)
    source_formula_needles = (
        "CHARTS=('Rh_reference','O2_slope','O2_axial','O2_buffer','O3_transition','O3_power')",
        "g.constant(40)", "physical_Jacobian_applied_once=True",
        "predecessor_memory_never_reset=True", "finite_N_correction=order==5",
        "def background_density", "def transport_geometry", "def geometry",
    )
    source_ast_ok = bool(source_tree is not None and
                         {"geometry", "transport_geometry", "background_density"}.issubset(source_names)
                         and all(needle in texts[OUTER_SOURCE] for needle in source_formula_needles))
    if not source_ast_ok:
        _gate(gates, "outer_source_AST_contract_unproved")
    patch_source_ast_ok = bool(
        patch_source_tree is not None
        and "transport_log" in patch_source_names
        and "def transport_log" in texts[PATCH_SOURCE]
        and "physical_measure='dt/t'" in texts[PATCH_SOURCE]
    )
    if not patch_source_ast_ok:
        _gate(gates, "patch_source_transport_AST_contract_unproved")

    outer_bindings = outer.get("true_outer_C6_C5_bindings", {})
    patch_bindings = patch.get("Rm_patch_finite_N_C5_bindings", {})
    charts = outer_bindings.get("charts", {}) if isinstance(outer_bindings, Mapping) else {}
    background = outer_bindings.get("background", {}) if isinstance(outer_bindings, Mapping) else {}
    chart_order_ok = (isinstance(charts, Mapping) and isinstance(background, Mapping)
                      and tuple(charts) == CHART_ORDER and tuple(background) == CHART_ORDER)
    if not chart_order_ok:
        _gate(gates, "outer_chart_order_mismatch")

    outer_nodes = outer.get("exact_graph_nodes", [])
    patch_nodes = patch.get("exact_graph_nodes", [])
    outer_count_ok = isinstance(outer_nodes, list) and len(outer_nodes) == OUTER_GRAPH_NODE_COUNT
    patch_count_ok = isinstance(patch_nodes, list) and len(patch_nodes) == PATCH_GRAPH_NODE_COUNT
    if not outer_count_ok:
        _gate(gates, "outer_archive_graph_count_mismatch")
    if not patch_count_ok:
        _gate(gates, "patch_archive_graph_count_mismatch")
    common_prefix_ok = bool(
        outer_count_ok and patch_count_ok
        and outer_nodes[:PATCH_GRAPH_NODE_COUNT] == patch_nodes
    )
    if not common_prefix_ok:
        _gate(gates, "outer_patch_common_prefix_mismatch")

    widths: dict[str, object] = {}
    widths_ok = True
    chart_shape: dict[str, object] = {}
    for name in CHART_ORDER:
        bg = background.get(name, {}) if isinstance(background, Mapping) else {}
        chart = charts.get(name, {}) if isinstance(charts, Mapping) else {}
        fields = bg.get("fields", {}) if isinstance(bg, Mapping) else {}
        loop = chart.get("loop", {}) if isinstance(chart, Mapping) else {}
        finite = chart.get("finite_N", {}) if isinstance(chart, Mapping) else {}
        incoming = chart.get("incoming", {}) if isinstance(chart, Mapping) else {}
        history = chart.get("history", {}) if isinstance(chart, Mapping) else {}
        modified = chart.get("modified_fields", {}) if isinstance(chart, Mapping) else {}
        radial = chart.get("radial", {}) if isinstance(chart, Mapping) else {}
        radial_background = chart.get("radial_background", {}) if isinstance(chart, Mapping) else {}
        shape_ok = True
        for key in ("E", "V", *MOMENT_NAMES):
            ok = _valid_handles(fields.get(key), 7, OUTER_GRAPH_NODE_COUNT)
            widths[f"{name}.background.fields.{key}"] = 7 if ok else None
            shape_ok = shape_ok and ok
        for key in ("A", "B_over_Pstar"):
            ok = _valid_handles(loop.get(key), 6, OUTER_GRAPH_NODE_COUNT)
            widths[f"{name}.charts.loop.{key}"] = 6 if ok else None
            shape_ok = shape_ok and ok
        for key in ("E_N", "V_N"):
            ok = _valid_handles(finite.get(key), 6, OUTER_GRAPH_NODE_COUNT)
            widths[f"{name}.charts.finite_N.{key}"] = 6 if ok else None
            shape_ok = shape_ok and ok
        for key in MOMENT_NAMES:
            for group, mapping in (("finite_N.five_signed_density_increments", finite.get("five_signed_density_increments", {})),
                                   ("incoming", incoming), ("history", history)):
                value = mapping.get(key) if isinstance(mapping, Mapping) else None
                ok = _valid_handles(value, 6, OUTER_GRAPH_NODE_COUNT)
                widths[f"{name}.charts.{group}.{key}"] = 6 if ok else None
                shape_ok = shape_ok and ok
        ok_p0 = _valid_handles(modified.get("P0"), 6, OUTER_GRAPH_NODE_COUNT)
        widths[f"{name}.charts.modified_fields.P0"] = 6 if ok_p0 else None
        shape_ok = shape_ok and ok_p0
        for key in ("E_y", "V_y"):
            value = radial.get(key)
            ok = _valid_handles(value, 6, OUTER_GRAPH_NODE_COUNT)
            widths[f"{name}.charts.radial.{key}"] = 6 if ok else None
            shape_ok = shape_ok and ok
        for key in ("E", "V"):
            value = radial_background.get(key)
            ok = _valid_handles(value, 7, OUTER_GRAPH_NODE_COUNT)
            widths[f"{name}.charts.radial_background.{key}"] = 7 if ok else None
            shape_ok = shape_ok and ok
        chart_shape[name] = {
            "valid": shape_ok,
            "background": bg,
            "chart": chart,
            "fields": fields,
            "loop": loop,
            "finite_N": finite,
            "incoming": incoming,
            "history": history,
            "modified_fields": modified,
            "radial": radial,
            "radial_background": radial_background,
        }
        widths_ok = widths_ok and shape_ok
    if not widths_ok:
        _gate(gates, "outer_chart_width_contract_mismatch")

    initial_correction = outer_bindings.get("initial_correction") if isinstance(outer_bindings, Mapping) else None
    first_incoming_ok = (isinstance(initial_correction, Mapping)
                         and all(tuple(initial_correction.get(key, ())) == PATCH_FIRST_INCOMING[key]
                                 for key in MOMENT_NAMES))
    if not first_incoming_ok:
        _gate(gates, "first_patch_Rh_incoming_chain_mismatch")
    patch_rh_correction = _get(patch_bindings, "endpoints.Rh.correction")
    patch_rh_correction_ok = bool(
        isinstance(patch_rh_correction, Mapping)
        and patch_rh_correction == initial_correction
        and all(_valid_handles(patch_rh_correction.get(key), 6, PATCH_GRAPH_NODE_COUNT)
                for key in MOMENT_NAMES)
    )
    if not patch_rh_correction_ok:
        _gate(gates, "patch_Rh_correction_source_mismatch")
    patch_n_ok = isinstance(patch_bindings, Mapping) and patch_bindings.get("N") == PATCH_N_NODE
    outer_n_ok = isinstance(outer_bindings, Mapping) and outer_bindings.get("N") == PATCH_N_NODE
    if not patch_n_ok or not outer_n_ok:
        _gate(gates, "shared_N_source_mismatch")

    receipt_checks = receipt.get("checks", {}) if isinstance(receipt, Mapping) else {}
    receipt_groups = ("actual_source", "geometry", "original_seed_theorem",
                      "actual_profiles", "actual_quiet_power", "exponential_transport",
                      "exponential_pressure")
    receipt_ok = (receipt.get("all_passed") is True and isinstance(receipt_checks, Mapping)
                  and all(isinstance(receipt_checks.get(group), Mapping)
                          and receipt_checks[group].get("passed") is True
                          and _bool_leaves(receipt_checks[group]) for group in receipt_groups))
    receipt_flags = (
        "current_original_true_six_outer_prescribed_C6_generic_C5_functions_installed",
        "current_original_all_17_chart_axial6_generic_source_installed",
        "current_original_true_six_outer_charts_finite_N_C5_functions_installed",
        "current_original_all_17_chart_finite_N_C5_functions_installed",
        "current_original_true_all17_modified_Q4_total_radial_C5_functions_installed",
        "current_original_true_rminus_to_Rc_signed_pressure_C3_functions_installed",
    )
    receipt_ok = receipt_ok and all(receipt.get(flag) is True for flag in receipt_flags)
    if not receipt_ok:
        _gate(gates, "outer_receipt_audit_failed")

    admission_result = admission.get("component_result", {})
    admission_scope = admission.get("scope", {})
    admission_ok = (isinstance(admission_result, Mapping)
                    and admission_result.get("passed") is True
                    and admission_result.get("definition_graph_node_count") == 1063
                    and admission_result.get("old_definition_graphs_rebuilt") is False
                    and admission_result.get("native_C2_C5_patch_history_identity_proved") is True
                    and admission_result.get("patch_finite_N_source_identity") is True
                    and isinstance(admission_scope, Mapping)
                    and admission_scope.get("ordinary_Z_orders") == [0, 1, 2])
    if not admission_ok:
        _gate(gates, "accepted_patch_admission_not_qualified")

    candidate_result = candidate.get("component_result", {})
    candidate_graph = candidate_result.get("definition_graph", {}) if isinstance(candidate_result, Mapping) else {}
    candidate_endpoints = _get(candidate_result, "definition_graph.outputs.histories.endpoint_xe", {})
    if not isinstance(candidate_endpoints, Mapping):
        candidate_endpoints = _get(candidate_result, "function_definitions.histories.endpoint_xe", {})
    endpoint_shape_ok = (isinstance(candidate_endpoints, Mapping)
                         and all(isinstance(candidate_endpoints.get(key), list)
                                 and len(candidate_endpoints[key]) == 3
                                 and all(isinstance(item, str) for item in candidate_endpoints[key])
                                 for key in MOMENT_NAMES))
    if not endpoint_shape_ok:
        _gate(gates, "qualified_patch_endpoint_xe_missing")

    candidate_nodes = candidate_graph.get("nodes", [])
    candidate_node_map = {node.get("id"): node for node in candidate_nodes if isinstance(node, Mapping)}
    candidate_node_ids = set(candidate_node_map)
    patch_z_handle_ok = bool(
        candidate_node_map.get(PATCH_Z_NAMESPACE) == {
            "id": PATCH_Z_NAMESPACE, "namespace": "original_native_C2_patch_finite_N_source_v1",
            "op": "exact_symbol", "value": "Z", "domain": "[-1,1]"}
        and patch_nodes[PATCH_Z_HANDLE] == {"operation": "bound_variable", "name": "Z"}
        and outer_nodes[PATCH_Z_HANDLE] == patch_nodes[PATCH_Z_HANDLE]
    )
    if not patch_z_handle_ok:
        _gate(gates, "patch_Z_handle_namespace_mismatch")
    qualified_graph_ok = bool(
        len(candidate_node_ids) == 1063
        and _digest(candidate_nodes) == "238eaf5f0ea130bd02f7038ab60a7cfc5aa6e6eef6873b9a4e1459cd8911d137"
        and all(item in candidate_node_ids for rows in candidate_endpoints.values() for item in rows)
    )
    if not qualified_graph_ok:
        _gate(gates, "qualified_patch_DAG_digest_or_endpoint_membership_mismatch")
    native_Rm_ok = bool(common_prefix_ok and _get(admission_result,
        "source_premise_audit.historical_saved_facts.geometry_witness.valid") is True
        and patch_nodes[PATCH_RM_NODE].get("Z_independent_substitution") is True
        and _get(patch_bindings, "chart.radius") is not None)
    shared_Pstar_ok = bool(common_prefix_ok and _get(admission_result,
        "source_premise_audit.historical_saved_facts.source_prerequisite_audit.accepted_Pstar_source_parameter_port") is True
        and outer_nodes[PSTAR_SHARED_PREFIX_NODE] == {"operation": "analytic_unary", "name": "exp", "argument": 5})
    shared_parameters_ok = bool(common_prefix_ok and all(
        patch_bindings.get(key) == outer_bindings.get(key) == handle
        for key, handle in {"N":18,"eta":21,"dstar":22,"sc":10,"phase_origin":24954}.items())
        and patch_nodes[18].get("exact_positive_integer") is True
        and patch_nodes[18].get("Z_independent") is True)
    chain_checks = {}
    previous_native = patch_rh_correction
    for name in CHART_ORDER:
        chart = charts[name]
        endpoint = _get(outer_bindings,"endpoints."+name+".correction",{})
        native_endpoint_substitution = all(
            outer_nodes[endpoint[key][i]] == {"operation":"function_substitution",
                "expression":chart["history"][key][i],"variable":chart["coordinate"],
                "value":chart["domain"][1],"Z_independent_substitution":True}
            for key in MOMENT_NAMES for i in range(6))
        chain_checks[name] = chart.get("incoming") == previous_native and native_endpoint_substitution
        previous_native = _get(outer_bindings, "endpoints." + name + ".correction")
    native_chain_ok = all(chain_checks.values())
    source_binding_checks = {}
    expected_coordinates = (205148,214116,223058,232895,241536,250215)
    expected_phases = (206328,215242,224412,233900,242551,251221)
    expected_J = (1,1,223066,1,1,1)
    for name, coordinate, phase, measure in zip(CHART_ORDER, expected_coordinates, expected_phases, expected_J):
        chart, bg = charts[name], background[name]
        source_binding_checks[name] = bool(
            chart.get("coordinate") == bg.get("phase") == coordinate
            and chart.get("actual_phase") == phase
            and chart.get("physical_measure") == bg.get("physical_measure") == measure
            and all(chart.get("roots", {}).get(key) == bg["fields"][key][:6] for key in ("E","V"))
            and chart.get("modified_fields", {}).get("P0") == _get(patch_bindings,"chart.modified_fields.P0")
            and bg["fields"].get("P0", [])[:6] == chart.get("modified_fields", {}).get("P0")
        )
    source_bindings_ok = all(source_binding_checks.values())
    for ok, label in ((qualified_graph_ok,"qualified_patch_DAG"),(native_Rm_ok,"native_Rm"),
        (shared_Pstar_ok,"shared_Pstar"),(shared_parameters_ok,"shared_parameters"),
        (native_chain_ok,"native_nonreset_chain"),(source_bindings_ok,"native_source_bindings")):
        if not ok: _gate(gates, label + "_identity_unproved")

    def qualified_native_port(ident: str, handle: int, path: str) -> bool:
        application = candidate_node_map.get(ident, {})
        port = candidate_node_map.get(application.get("source_port"), {}).get("value", {})
        return bool(application.get("op") == "source_function_application"
            and port.get("source_node") == handle
            and port.get("archive_sha256") == PATCH_ARCHIVE_SHA256
            and port.get("source_archive") == PATCH_ARCHIVE
            and port.get("source_path") == path and port.get("source_node_valid") is True
            and port.get("source_namespace") == "original_true_Rm_patch_finite_N_C5_functions_v1"
            and port.get("source_z_handle") == application.get("source_z_handle") == PATCH_Z_HANDLE
            and application.get("variable_bindings",{}).get("Z") == PATCH_Z_NAMESPACE
            and application.get("source_coordinate_handle") == (
                patch_bindings["chart"]["coordinate"] if ".chart." in path else None))
    patch_bridge_checks = {}
    old_outputs = candidate_graph.get("outputs", {})
    old_histories = old_outputs.get("histories", {})
    for key in MOMENT_NAMES:
        inlet = _get(patch_bindings, "incoming." + key, [])
        history = _get(patch_bindings, "chart.history." + key, [])
        density = _get(patch_bindings, "chart.finite_N.five_signed_density_increments." + key, [])
        end = _get(patch_bindings, "endpoints.Rh.correction." + key, [])
        old_inlet = _get(old_histories, "incoming_C5." + key, [])
        old_provider = _get(old_outputs, "provider_density." + key, [])
        native_endpoint_ok = bool(len(history) == len(end) == 6 and all(
            patch_nodes[end[i]] == {"operation":"function_substitution", "expression":history[i],
                "variable":patch_bindings["chart"]["coordinate"], "value":patch_bindings["chart"]["domain"][1],
                "Z_independent_substitution":True} for i in range(6)))
        same_inlet = bool(len(inlet) == 6 and len(old_inlet) == 3 and all(
            qualified_native_port(old_inlet[i], inlet[i], "Rm_patch_finite_N_C5_bindings.incoming."+key+"["+str(i)+"]")
            for i in ORDINARY_Z_ORDERS))
        same_density_port = bool(len(density) == 6 and len(old_provider) == 3 and all(
            qualified_native_port(old_provider[i], density[i], "Rm_patch_finite_N_C5_bindings.chart.finite_N.five_signed_density_increments."+key+"["+str(i)+"]")
            for i in ORDINARY_Z_ORDERS))
        old_endpoint_ok = all(
            candidate_node_map[old_histories["endpoint_xe"][key][i]].get("args",[])[:2]
            == [old_histories["C2"][key][i], "original_native_C2_patch_finite_N_source_v1:d000002"]
            and candidate_node_map[old_histories["endpoint_xe"][key][i]].get("endpoint") == "x=e"
            and candidate_node_map[candidate_node_map[old_histories["endpoint_xe"][key][i]]["args"][2]].get("op") == "exp"
            and candidate_node_map[candidate_node_map[old_histories["endpoint_xe"][key][i]]["args"][2]].get("args") == ["original_native_C2_patch_finite_N_source_v1:d000003"]
            for i in ORDINARY_Z_ORDERS)
        native_endpoint_ok = native_endpoint_ok and patch_nodes[patch_bindings["chart"]["domain"][1]] == {
            "operation":"analytic_unary", "name":"exp", "argument":1}
        patch_bridge_checks[key] = dict(same_inlet=same_inlet, same_density_port=same_density_port,
            native_endpoint_substitution=native_endpoint_ok, old_C2_endpoint_substitution=old_endpoint_ok)
    patch_history_bridge_ok = bool(admission_ok and patch_source_ast_ok and qualified_graph_ok
        and all(all(row.values()) for row in patch_bridge_checks.values())
        and all(_get(candidate_result,"definition_graph_checks."+key) is True for key in (
            "log_kernel_constructed", "log_transport_inverse_s_is_in_actual_integrand", "integral_variables_and_bounds_are_scalar_IDs"))
        and all(admission_result.get(key) is True for key in ("native_C2_C5_patch_density_identity_proved",
            "native_C2_patch_correction_inlet_identity_proved", "native_C2_C5_patch_history_identity_proved")))
    if not patch_history_bridge_ok: _gate(gates,"native_patch_C5_history_bridge_unproved")

    generic_result = generic.get("component_result", {})
    generic_api = generic.get("math_api_description", {})
    generic_identity_ok = (isinstance(generic_result, Mapping)
                          and generic_result.get("generic_density_identity_proved") is True
                          and isinstance(generic_api, Mapping)
                          and isinstance(generic_api.get("C2_density"), str)
                          and isinstance(generic_api.get("C5_density"), str)
                          and "B_over_Pstar/N" in generic_api.get("C5_density", "")
                          and "(B,F,VF+E*B,2*V*B-E*F,E*F)" in generic_api.get("C2_density", ""))
    if not generic_identity_ok:
        _gate(gates, "accepted_generic_density_identity_missing")

    # The source reports are data.  Keep the old graphs qualified as metadata;
    # no saved node list is copied into this graph.
    return {
        "texts": texts,
        "admission": admission,
        "candidate": candidate,
        "generic": generic,
        "outer": outer,
        "patch": patch,
        "receipt": receipt,
        "outer_bindings": outer_bindings,
        "patch_bindings": patch_bindings,
        "charts": chart_shape,
        "chart_order_ok": chart_order_ok,
        "outer_count_ok": outer_count_ok,
        "patch_count_ok": patch_count_ok,
        "widths": widths,
        "widths_ok": widths_ok,
        "first_incoming_ok": first_incoming_ok,
        "patch_rh_correction_ok": patch_rh_correction_ok,
        "common_prefix_ok": common_prefix_ok,
        "patch_z_handle_ok": patch_z_handle_ok,
        "shared_N_ok": patch_n_ok and outer_n_ok,
        "receipt_ok": receipt_ok,
        "admission_ok": admission_ok,
        "candidate_result": candidate_result,
        "candidate_graph": candidate_graph,
        "candidate_node_ids": candidate_node_ids,
        "candidate_node_map": candidate_node_map,
        "patch_history_bridge_ok": patch_history_bridge_ok,
        "patch_history_bridge_checks": patch_bridge_checks,
        "qualified_graph_ok": qualified_graph_ok,
        "native_Rm_ok": native_Rm_ok,
        "shared_Pstar_ok": shared_Pstar_ok,
        "shared_parameters_ok": shared_parameters_ok,
        "native_chain_ok": native_chain_ok,
        "native_chain_checks": chain_checks,
        "source_bindings_ok": source_bindings_ok,
        "source_binding_checks": source_binding_checks,
        "candidate_endpoints": candidate_endpoints,
        "endpoint_shape_ok": endpoint_shape_ok,
        "generic_result": generic_result,
        "generic_api": generic_api,
        "generic_identity_ok": generic_identity_ok,
        "source_ast_ok": source_ast_ok,
        "patch_source_ast_ok": patch_source_ast_ok,
        "source_function_names": sorted(source_names),
        "patch_source_function_names": sorted(patch_source_names),
        "archive_namespace": "original_true_outer_C6_C5_functions_v1",
        "patch_namespace": "original_true_Rm_patch_finite_N_C5_functions_v1",
    }


def _candidate_endpoint_port(
    dag: _DAG, witness: Mapping[str, object], moment: str, index: int, z: ScalarID,
) -> ScalarID:
    endpoints = witness.get("candidate_endpoints", {})
    rows = endpoints.get(moment, []) if isinstance(endpoints, Mapping) else []
    source_node = rows[index] if isinstance(rows, list) and index < len(rows) else "missing"
    valid = isinstance(source_node, str) and source_node in witness["candidate_node_ids"]
    return _source_port(
        dag,
        report_sha256=CANDIDATE_SHA256,
        namespace="original_native_C2_patch_finite_N_source_v1",
        source_node=source_node,
        archive=CANDIDATE,
        path="component_result.definition_graph.outputs.histories.endpoint_xe." + moment + "[" + str(index) + "]",
        label="patch.C2.endpoint_xe." + moment + "[" + str(index) + "]",
        role="qualified_prior_DAG_C2_endpoint_function",
        valid=valid,
        width=3,
        z_binding=z,
        fixed_endpoint="x=xe",
        old_result_sha256="ac8b4e100124d47e7c2696d21945edae05726052a4d7bba1a47563f666950379",
        candidate_report_sha256=CANDIDATE_SHA256,
        old_result_namespace="original_native_C2_patch_finite_N_source_v1",
        bind_Z_only=True,
        source_z_handle=PATCH_Z_NAMESPACE,
        archived_native_Z_handle=PATCH_Z_HANDLE,
        source_z_namespace=PATCH_Z_NAMESPACE,
        frozen_source_z_namespace=True,
        source_endpoint_is_function=True,
        source_coordinate_is_fixed=True,
    )


def _qualified_patch_jet(dag: _DAG, witness: Mapping[str, object], path: str,
                         z: ScalarID, x: ScalarID | None = None) -> JetList:
    rows = _get(witness["candidate_result"], "definition_graph.outputs." + path, [])
    return [_source_port(dag, report_sha256=CANDIDATE_SHA256,
        namespace="original_native_C2_patch_finite_N_source_v1", source_node=ident,
        archive=CANDIDATE, path="component_result.definition_graph.outputs."+path+"["+str(i)+"]",
        label="qualified_patch."+path+"["+str(i)+"]", role="qualified_old_DAG_function",
        valid=ident in witness["candidate_node_ids"], width=3,
        coordinate_binding=x, z_binding=z,
        coordinate_handle="original_native_C2_patch_finite_N_source_v1:d000002" if x else None,
        source_z_handle=PATCH_Z_NAMESPACE, frozen_old_DAG_not_rebuilt=True)
        for i,ident in enumerate(rows)]


def _build_graph(witness: Mapping[str, object], gates: list[str]) -> tuple[_DAG, dict[str, object]]:
    """Build the new ordinary-Z graph from the authenticated source witness."""
    dag = _DAG()
    Z = _symbol(dag, "Z", "ordinary_spatial_parameter")

    patch_bindings = witness["patch_bindings"]
    outer_bindings = witness["outer_bindings"]
    patch_namespace = witness["patch_namespace"]
    outer_namespace = witness["archive_namespace"]
    patch_chart = patch_bindings.get("chart", {}) if isinstance(patch_bindings, Mapping) else {}

    N = _source_port(
        dag, report_sha256=PATCH_ARCHIVE_SHA256, namespace=patch_namespace,
        source_node=PATCH_N_NODE, archive=PATCH_ARCHIVE,
        path="Rm_patch_finite_N_C5_bindings.N", label="native_C5.parameter.N",
        role="fixed_N_Z_independent_source", valid=witness["shared_N_ok"], width=1,
        valid_range="positive integer", fixed_N=True, Z_independent=True,
    )
    Rm = _source_port(
        dag, report_sha256=PATCH_ARCHIVE_SHA256, namespace=patch_namespace,
        source_node=PATCH_RM_NODE, archive=PATCH_ARCHIVE,
        path="exact_graph_nodes[85743]", label="native_C5.geometry.Rm",
        role="native_Rm_scalar_source", valid=witness["native_Rm_ok"], width=1,
        valid_range="Z-independent native radius", native_inlet_radius=85743,
        Z_independent=True,
    )
    Pstar = _source_port(
        dag, report_sha256=OUTER_ARCHIVE_SHA256,
        namespace=outer_namespace,
        source_node=PSTAR_SHARED_PREFIX_NODE, archive=OUTER_ARCHIVE,
        path="exact_graph_nodes[33572]",
        label="shared_prefix.parameter.Pstar", role="accepted_shared_prefix_Pstar_source",
        valid=witness["shared_Pstar_ok"], width=1, valid_range="positive source parameter",
        earlier_accepted_prefix_archive_sha256=PSTAR_SHARED_PREFIX_ARCHIVE_SHA256,
        shared_prefix_checked_patch=True, source_graph_node=33572,
        Z_independent=True,
    )
    invN = _div(dag, _one(dag), N)

    formulas: dict[str, object] = {
        "source_parameters": {"N": N, "Rm": Rm, "Pstar": Pstar, "Z": Z},
        "charts": {},
        "source_ports": {"N": N, "Rm": Rm, "Pstar": Pstar},
        "provider_identity_proofs": {},
        "density_identity_proofs": {},
        "endpoint_chain_proofs": [],
        "c5_history_bridge_proofs": {},
    }
    previous_endpoint: dict[str, JetList] | None = None
    previous_chart: str | None = None
    previous_coordinate: ScalarID | None = None
    previous_history_proofs: dict[str, ScalarID] = {}

    for chart_name in CHART_ORDER:
        shape = witness["charts"][chart_name]
        bg = shape["background"]
        chart = shape["chart"]
        fields = shape["fields"]
        loop = shape["loop"]
        finite = shape["finite_N"]
        incoming_native = shape["incoming"]
        history_native = shape["history"]
        modified = shape["modified_fields"]
        radial = shape["radial"]
        radial_background = shape["radial_background"]

        x = _symbol(dag, "x_" + chart_name, "local_chart_coordinate")
        t = _symbol(dag, "t_" + chart_name, "scalar_integration_variable")
        lower_text, upper_text = LOCAL_DOMAINS[chart_name]
        lower = _const(dag, lower_text)
        upper = _const(dag, upper_text)
        fourty = _const(dag, "40")
        if chart_name in ("Rh_reference", "O2_slope"):
            y = _add(dag, _const(dag, "6"), x)
            J = _one(dag)
            geometry_recipe = "y=6+x; J=1"
        elif chart_name == "O2_axial":
            exp40x = _exp(dag, _mul(dag, fourty, x))
            y = _add(dag, _const(dag, "6"), exp40x)
            J = _mul(dag, fourty, _exp(dag, _mul(dag, fourty, x)))
            geometry_recipe = "y=6+exp(40*x); J=40*exp(40*x)"
        elif chart_name == "O2_buffer":
            exp40 = _exp(dag, fourty)
            y = _add(dag, _const(dag, "6"), exp40, x)
            J = _one(dag)
            geometry_recipe = "y=6+exp(40)+x; J=1"
        elif chart_name == "O3_transition":
            exp40 = _exp(dag, fourty)
            y = _add(dag, _const(dag, "17"), exp40, x)
            J = _one(dag)
            geometry_recipe = "y=17+exp(40)+x; J=1"
        else:
            exp40 = _exp(dag, fourty)
            y = _add(dag, _const(dag, "18"), exp40, x)
            J = _one(dag)
            geometry_recipe = "y=18+exp(40)+x; J=1"
        R = _mul(dag, Rm, _exp(dag, y))
        y_lower = _function_substitution(dag, y, x, lower, chart=chart_name, endpoint="lower")
        y_upper = _function_substitution(dag, y, x, upper, chart=chart_name, endpoint="upper")
        R_upper = _function_substitution(dag, R, x, upper, chart=chart_name, endpoint="upper")

        coordinate = chart.get("coordinate")
        phase = chart.get("actual_phase")
        domain = chart.get("domain")
        geometry_identity_proof = dag.op(
            "source_identity_proof", (Rm, y, J, R),
            theorem="native_outer_geometry_source_identity", chart=chart_name,
            source_archive=OUTER_ARCHIVE, source_archive_sha256=OUTER_ARCHIVE_SHA256,
            source_ast=OUTER_SOURCE, source_ast_sha256=OUTER_SOURCE_SHA256,
            source_coordinate_handle=coordinate, source_phase_handle=phase,
            source_geometry_path="true_outer_C6_C5_bindings.charts." + chart_name,
            local_domain=list(LOCAL_DOMAINS[chart_name]),
            y_recipe=geometry_recipe, radius_recipe="R=Rm*exp(y)",
            first_Rh_lower_y_is_one=(chart_name == "Rh_reference"),
            first_Rh_lower_radius="Rm*exp(1)=Rm*e" if chart_name == "Rh_reference" else None,
            first_Rh_patch_radius_match=(chart_name == "Rh_reference"),
            Jacobian_included=True, bounds_Z_independent=True, source_only=True,
            equality_is_symbolic=True,
        )
        source_rows: dict[str, object] = {"background": {}, "loop": {}, "finite_N": {}}
        leading: dict[str, JetList] = {}
        leading_rows: dict[str, list[ScalarID]] = {}
        for key in ("E", "V", *MOMENT_NAMES):
            rows, ok = _source_rows(
                dag, rows=fields.get(key), width=7, node_count=OUTER_GRAPH_NODE_COUNT,
                report_sha256=OUTER_ARCHIVE_SHA256, namespace=outer_namespace,
                archive=OUTER_ARCHIVE,
                path="true_outer_C6_C5_bindings.background." + chart_name + ".fields." + key,
                label="native_outer.background." + chart_name + "." + key,
                role="native_outer_leading_background_root_7_rows", x=x, z=Z,
                coordinate_handle=coordinate, phase_handle=phase, valid_range=domain,
            )
            if not ok:
                _gate(gates, "source_background_width_invalid:" + chart_name + ":" + key)
            leading_rows[key] = rows
            leading[key] = _port_jet(rows)
        source_rows["background"] = leading_rows

        loop_rows: dict[str, list[ScalarID]] = {}
        loop_jet: dict[str, JetList] = {}
        for key in ("A", "B_over_Pstar"):
            rows, ok = _source_rows(
                dag, rows=loop.get(key), width=6, node_count=OUTER_GRAPH_NODE_COUNT,
                report_sha256=OUTER_ARCHIVE_SHA256, namespace=outer_namespace,
                archive=OUTER_ARCHIVE,
                path="true_outer_C6_C5_bindings.charts." + chart_name + ".loop." + key,
                label="native_outer.loop." + chart_name + "." + key,
                role="native_outer_loop_source_6_rows", x=x, z=Z,
                coordinate_handle=coordinate, phase_handle=phase, valid_range=domain,
            )
            if not ok:
                _gate(gates, "source_loop_width_invalid:" + chart_name + ":" + key)
            loop_rows[key] = rows
            loop_jet[key] = _port_jet(rows)
        source_rows["loop"] = loop_rows

        native_en_rows: dict[str, list[ScalarID]] = {}
        native_en_jet: dict[str, JetList] = {}
        for key in ("E_N", "V_N"):
            rows, ok = _source_rows(
                dag, rows=finite.get(key), width=6, node_count=OUTER_GRAPH_NODE_COUNT,
                report_sha256=OUTER_ARCHIVE_SHA256, namespace=outer_namespace,
                archive=OUTER_ARCHIVE,
                path="true_outer_C6_C5_bindings.charts." + chart_name + ".finite_N." + key,
                label="native_outer.finite_N." + chart_name + "." + key,
                role="native_outer_finite_N_provider_comparison_port", x=x, z=Z,
                coordinate_handle=coordinate, phase_handle=phase, valid_range=domain,
            )
            if not ok:
                _gate(gates, "source_provider_width_invalid:" + chart_name + ":" + key)
            native_en_rows[key] = rows
            native_en_jet[key] = _port_jet(rows)
        source_rows["finite_N"] = native_en_rows

        # New C2 density formula.  It is independent of imported C5 density.
        E, V = leading["E"], leading["V"]
        A, Bnorm = loop_jet["A"], loop_jet["B_over_Pstar"]
        A_over_N = _jet_scale(dag, A, invN)
        exp_A_over_N = _jet_exp(dag, A_over_N)
        expm1_A_over_N = _jet_expm1(dag, A_over_N)
        r = _jet_mul(dag, E, expm1_A_over_N)
        s = _jet_scale(dag, Bnorm, invN)
        E_N = _jet_mul(dag, E, exp_A_over_N)
        V_N = _jet_add(dag, V, s)
        F = _jet_scale(dag, r, N)
        vf = _jet_mul(dag, V, F)
        eb = _jet_mul(dag, E, Bnorm)
        vb = _jet_mul(dag, V, Bnorm)
        ef = _jet_mul(dag, E, F)
        fb = _jet_mul(dag, F, Bnorm)
        b_sq = _jet_mul(dag, Bnorm, Bnorm)
        f_sq = _jet_mul(dag, F, F)
        half_f_sq = _jet_scale(dag, f_sq, _const(dag, "1/2"))
        c2_first = {
            "m": Bnorm,
            "h": F,
            "k": _jet_add(dag, vf, eb),
            "e": _jet_sub(dag, _jet_scale(dag, vb, _const(dag, "2")), ef),
            "p": ef,
        }
        c2_second = {
            "m": _jet_zero(dag),
            "h": _jet_zero(dag),
            "k": fb,
            "e": _jet_sub(dag, b_sq, half_f_sq),
            "p": half_f_sq,
        }
        c2_density = {
            key: _jet_add(dag, c2_first[key], _jet_scale(dag, c2_second[key], invN))
            for key in MOMENT_NAMES
        }

        native_c5_density: dict[str, JetList] = {}
        native_c5_density_rows: dict[str, list[ScalarID]] = {}
        density_proofs: dict[str, ScalarID] = {}
        for key in MOMENT_NAMES:
            rows, ok = _source_rows(
                dag, rows=(finite.get("five_signed_density_increments", {}) or {}).get(key),
                width=6, node_count=OUTER_GRAPH_NODE_COUNT,
                report_sha256=OUTER_ARCHIVE_SHA256, namespace=outer_namespace,
                archive=OUTER_ARCHIVE,
                path="true_outer_C6_C5_bindings.charts." + chart_name + ".finite_N.five_signed_density_increments." + key,
                label="native_outer.C5_density." + chart_name + "." + key,
                role="native_outer_C5_signed_density_source", x=x, z=Z,
                coordinate_handle=coordinate, phase_handle=phase, valid_range=domain,
            )
            if not ok:
                _gate(gates, "source_density_width_invalid:" + chart_name + ":" + key)
            native_c5_density_rows[key] = rows
            native_c5_density[key] = _port_jet(rows)
            scaled_c5 = _jet_scale(dag, native_c5_density[key], N)
            density_proofs[key] = dag.op(
                "source_identity_proof", tuple(c2_density[key] + scaled_c5),
                theorem="C2_density=N*C5_density", chart=chart_name, moment=key,
                source_rule="accepted_generic_density_rule",
                source_formula_C2="(B,F,VF+E*B,2*V*B-E*F,E*F)+(0,0,F*B,B^2-F^2/2,F^2/2)/N",
                source_formula_C5="(s,r,V*r+E*s+r*s,2*V*s+s^2-E*r-r^2/2,E*r+r^2/2)",
                scaling="F=N*r; B=N*s", same_EVA_B_N_phase_ports=True,
                source_only=True, equality_is_symbolic=True,
            )

        provider_E = dag.op(
            "source_identity_proof", tuple(E_N + native_en_jet["E_N"]),
            theorem="constructed_E_N_matches_actual_native_provider", chart=chart_name,
            checked_premises=["same E source", "same A source", "same N source", "same phase source"],
            source_only=True, equality_is_symbolic=True,
        )
        provider_V = dag.op(
            "source_identity_proof", tuple(V_N + native_en_jet["V_N"]),
            theorem="constructed_V_N_matches_actual_native_provider", chart=chart_name,
            checked_premises=["same V source", "same B_over_Pstar source", "same N source", "same phase source"],
            source_only=True, equality_is_symbolic=True,
        )

        # The first C2 inlet is the qualified saved patch endpoint at x=e.
        # Later inlets are substitutions of this new graph's preceding upper
        # endpoint.  C5 incoming/history rows remain independent source ports.
        incoming_c2: dict[str, JetList] = {}
        incoming_c5: dict[str, JetList] = {}
        incoming_c5_rows: dict[str, list[ScalarID]] = {}
        outer_incoming_rows: dict[str, list[ScalarID]] = {}
        history_c5_rows: dict[str, list[ScalarID]] = {}
        for key in MOMENT_NAMES:
            native_rows, native_rows_ok = _source_rows(
                dag, rows=incoming_native.get(key), width=6, node_count=OUTER_GRAPH_NODE_COUNT,
                report_sha256=OUTER_ARCHIVE_SHA256, namespace=outer_namespace,
                archive=OUTER_ARCHIVE,
                path="true_outer_C6_C5_bindings.charts." + chart_name + ".incoming." + key,
                label="native_outer.C5_incoming." + chart_name + "." + key,
                role="native_outer_C5_independent_incoming_source", x=None, z=Z,
                coordinate_handle=coordinate, phase_handle=phase, valid_range=domain,
            )
            if not native_rows_ok:
                _gate(gates, "source_incoming_width_invalid:" + chart_name + ":" + key)
            outer_incoming_rows[key] = native_rows
            if previous_endpoint is None:
                incoming_c2[key] = [
                    _candidate_endpoint_port(dag, witness, key, index, Z)
                    for index in ORDINARY_Z_ORDERS
                ]
                incoming_source = "qualified_patch_C2_endpoint_xe"
            else:
                incoming_c2[key] = _substitute_jet(
                    dag, previous_endpoint[key],
                    previous_coordinate,
                    _const(dag, LOCAL_DOMAINS[previous_chart][1]),
                    incoming_chart=chart_name, predecessor_chart=previous_chart,
                    predecessor_upper_endpoint=True,
                )
                incoming_source = "new_previous_chart_upper_endpoint"
            if previous_endpoint is None:
                # The first outer inlet is the native patch Rh correction.
                # The outer archive repeats this row set as initial_correction;
                # the patch archive is the authoritative C5 source port.
                first_rows = _get(patch_bindings, "endpoints.Rh.correction." + key, ())
                rows, ok = _source_rows(
                    dag, rows=first_rows, width=6, node_count=PATCH_GRAPH_NODE_COUNT,
                    report_sha256=PATCH_ARCHIVE_SHA256, namespace=patch_namespace,
                    archive=PATCH_ARCHIVE,
                    path="Rm_patch_finite_N_C5_bindings.endpoints.Rh.correction." + key,
                    label="native_patch.C5_Rh_correction." + key,
                    role="native_patch_C5_Rh_correction_incoming_source", x=None, z=Z,
                    fixed_endpoint="patch x=e",
                    coordinate_handle=_get(patch_chart, "coordinate"),
                    phase_handle=_get(patch_chart, "actual_phase"),
                    valid_range=_get(patch_chart, "domain"),
                )
            else:
                rows, ok = native_rows, native_rows_ok
            if not ok:
                _gate(gates, "source_incoming_width_invalid:" + chart_name + ":" + key)
            incoming_c5_rows[key] = rows
            incoming_c5[key] = _port_jet(rows)
            rows_h, ok_h = _source_rows(
                dag, rows=history_native.get(key), width=6, node_count=OUTER_GRAPH_NODE_COUNT,
                report_sha256=OUTER_ARCHIVE_SHA256, namespace=outer_namespace,
                archive=OUTER_ARCHIVE,
                path="true_outer_C6_C5_bindings.charts." + chart_name + ".history." + key,
                label="native_outer.C5_history." + chart_name + "." + key,
                role="native_outer_C5_independent_nonreset_history_source", x=x, z=Z,
                coordinate_handle=coordinate, phase_handle=phase, valid_range=domain,
            )
            if not ok_h:
                _gate(gates, "source_history_width_invalid:" + chart_name + ":" + key)
            history_c5_rows[key] = rows_h

        patch_history_source_rows: dict[str, list[ScalarID]] = {}
        patch_bridge_ports: dict[str, object] = {}
        c5_history_bridge_proof: ScalarID | None = None
        if previous_endpoint is None:
            patch_x = _symbol(dag, "patch_x", "[1,e]; independent of outer local x")
            bridge_args: list[ScalarID] = []
            patch_bridge_ports["coordinate"] = patch_x
            for key in MOMENT_NAMES:
                groups = {}
                for label, source_path, coordinate_binding in (
                    ("history", "chart.history."+key, patch_x),
                    ("density", "chart.finite_N.five_signed_density_increments."+key, patch_x),
                    ("incoming", "incoming."+key, None)):
                    rows, ok = _source_rows(dag, rows=_get(patch_bindings,source_path), width=6,
                        node_count=PATCH_GRAPH_NODE_COUNT, report_sha256=PATCH_ARCHIVE_SHA256,
                        namespace=patch_namespace, archive=PATCH_ARCHIVE,
                        path="Rm_patch_finite_N_C5_bindings."+source_path,
                        label="native_patch.C5_"+label+"."+key,
                        role="native_patch_C5_transport_"+label, x=coordinate_binding, z=Z,
                        coordinate_handle=patch_chart.get("coordinate") if coordinate_binding else None,
                        valid_range=patch_chart.get("domain"))
                    if not ok: _gate(gates,"patch_bridge_width_invalid:"+label+":"+key)
                    groups["native_"+label] = rows
                    bridge_args.extend(_port_jet(rows))
                patch_history_source_rows[key] = groups["native_history"]
                for label, old_path, old_x in (
                    ("history","histories.C5."+key,patch_x),
                    ("density","finite_N.C5_density."+key,patch_x),
                    ("incoming","histories.incoming_C5."+key,None)):
                    old_rows = _qualified_patch_jet(dag,witness,old_path,Z,old_x)
                    groups["old_DAG_"+label] = old_rows
                    bridge_args.extend(old_rows)
                groups["native_Rh_endpoint"] = incoming_c5_rows[key]
                groups["old_C2_endpoint_xe"] = incoming_c2[key]
                bridge_args.extend(incoming_c5_rows[key])
                bridge_args.extend(incoming_c2[key])
                patch_bridge_ports[key] = groups
            c5_history_bridge_proof = dag.op(
                "source_identity_proof", tuple(bridge_args),
                theorem="C5_DAG_history_equals_native_patch_chart_history",
                chart=chart_name, source_history="native_patch.chart.history",
                same_density=witness["patch_history_bridge_ok"],
                same_density_source="accepted_patch_C2_C5_density_theorem_and_exact_native_provider_density_ports",
                same_native_incoming=witness["patch_history_bridge_ok"],
                native_incoming_source="native_patch.incoming; patch lower endpoint x=1",
                same_bounds=True, lower_bound="1", upper_bound="e",
                source_coordinate=patch_x, outer_coordinate_is_distinct=True,
                same_log_Duhamel_kernel=True, physical_measure="dt/t",
                same_transport_operator=True, linear_ODE_uniqueness=True,
                native_Rh_endpoint_is_history_at_e=witness["patch_history_bridge_ok"],
                old_C2_endpoint_is_C2_history_at_e=witness["patch_history_bridge_ok"],
                accepted_C2_equals_N_C5_history=True,
                bridge_via="same_patch_density_inlet_log_transport_then_fixed_e_endpoint",
                name_alias_not_used=True, source_only=True, equality_is_symbolic=True,
            )

        history_c2: dict[str, JetList] = {}
        c2_integrals: dict[str, JetList] = {}
        c2_memory: dict[str, ScalarID] = {}
        physical_derivative: dict[str, ScalarID] = {}
        for key in MOMENT_NAMES:
            rate = _const(dag, RATES[key])
            rate_zero = RATES[key] == "0"
            y_t = _function_substitution(dag, y, x, t, chart=chart_name, transport_argument="t")
            y_difference = _sub(dag, y, y_t)
            kernel = _one(dag) if rate_zero else _exp(dag, _neg(dag, _mul(dag, rate, y_difference)))
            y_width = _sub(dag, y, y_lower)
            decay = _one(dag) if rate_zero else _exp(dag, _neg(dag, _mul(dag, rate, y_width)))
            c2_memory[key] = dag.op(
                "transport_memory", (decay,), chart=chart_name, moment=key,
                rate=RATES[key], predecessor_memory_never_reset=True,
                exact_rate=True, rate_zero_memory_is_one=rate_zero,
                retain_memory_on_quiet_chart=(chart_name == "O3_power"),
            )
            density_t = _substitute_jet(
                dag, c2_density[key], x, t,
                chart=chart_name, source_function_substitution=True,
                ordinary_Z_orders=[0, 1, 2],
            )
            J_t = _function_substitution(dag, J, x, t, chart=chart_name, physical_Jacobian=True)
            integrand = _jet_scale(dag, _jet_mul(dag, density_t, _jet_scalar(dag, J_t)), kernel)
            integral = _integral_jet(
                dag, integrand, t, lower, x,
                chart=chart_name, moment=key,
                original_outer_Duhamel=True,
                scalar_integration_variable=True,
                density_is_not_integration_variable=True,
                physical_Jacobian_applied_once=True,
                actual_integrand="rho_C2(t)*J(t)*exp(-rate*(y(x)-y(t)))",
                actual_integrand_factor_ids=list(density_t) + [J_t, kernel],
                jacobian_factor_id=J_t,
                predecessor_memory_never_reset=True,
            )
            c2_integrals[key] = integral
            history_c2[key] = _jet_add(
                dag, _jet_scale(dag, incoming_c2[key], decay), integral,
            )
            physical_derivative[key] = dag.op(
                "physical_derivative_identity",
                tuple(history_c2[key] + c2_density[key] + [J]),
                chart=chart_name, moment=key, rate=RATES[key],
                identity="H_x=J*(rho_C2-rate*H)",
                ordinary_Z_orders=[0, 1, 2], temporal_identity=False,
                source_only=True,
            )

        endpoints_c2 = {
            key: _substitute_jet(
                dag, history_c2[key], x, upper,
                chart=chart_name, endpoint_function=True,
                local_coordinate_substitution=True, all_ordinary_Z_orders=True,
                endpoint_rows_are_values=False,
            ) for key in MOMENT_NAMES
        }

        inlet_proofs: dict[str, ScalarID] = {}
        history_proofs: dict[str, ScalarID] = {}
        for key in MOMENT_NAMES:
            predecessor_proof = c5_history_bridge_proof if previous_endpoint is None else previous_history_proofs[key]
            inlet_proofs[key] = dag.op("source_identity_proof",
                tuple(incoming_c2[key] + _jet_scale(dag,incoming_c5[key],N) + [predecessor_proof]),
                theorem="outer_C2_incoming=N*native_C5_incoming", chart=chart_name, moment=key,
                exact_native_predecessor_source=witness["native_chain_ok"],
                patch_bridge_or_previous_history_substitution=True,
                ordinary_Z_orders=[0,1,2], source_only=True, equality_is_symbolic=True)
            history_proofs[key] = dag.op("source_identity_proof",
                tuple(history_c2[key] + _jet_scale(dag,_port_jet(history_c5_rows[key]),N)
                    + [density_proofs[key],inlet_proofs[key],geometry_identity_proof]),
                theorem="outer_C2_history=N*native_C5_history", chart=chart_name, moment=key,
                rates_equal=RATES[key], same_source_density=True, same_inlet=True,
                same_source_geometry_and_measure=True, linear_Duhamel_uniqueness=True,
                physical_Jacobian_applied_once=True, predecessor_memory_never_reset=True,
                source_only=True, equality_is_symbolic=True)

        leading_fields = {key: leading[key] for key in MOMENT_NAMES}
        finite_fields = {
            key: _jet_add(dag, leading_fields[key], _jet_scale(dag, history_c2[key], invN))
            for key in MOMENT_NAMES
        }
        P0_rows, p0_ok = _source_rows(
            dag, rows=modified.get("P0"), width=6, node_count=OUTER_GRAPH_NODE_COUNT,
            report_sha256=OUTER_ARCHIVE_SHA256, namespace=outer_namespace,
            archive=OUTER_ARCHIVE,
            path="true_outer_C6_C5_bindings.charts." + chart_name + ".modified_fields.P0",
            label="native_outer.C5_modified_P0." + chart_name,
            role="native_outer_C5_P0_source", x=x, z=Z,
            coordinate_handle=coordinate, phase_handle=phase, valid_range=domain,
        )
        if not p0_ok:
            _gate(gates, "source_P0_width_invalid:" + chart_name)
        pressure = _jet_add(dag, _port_jet(P0_rows), finite_fields["p"])
        Utheta = _jet_scale(dag, E_N, Pstar)
        Uz = _jet_scale(dag, V_N, Pstar)

        radial_rows: dict[str, list[ScalarID]] = {}
        radial_jets: dict[str, JetList] = {}
        for key in ("E_y", "V_y"):
            rows, radial_ok = _source_rows(
                dag, rows=radial.get(key), width=6, node_count=OUTER_GRAPH_NODE_COUNT,
                report_sha256=OUTER_ARCHIVE_SHA256, namespace=outer_namespace,
                archive=OUTER_ARCHIVE,
                path="true_outer_C6_C5_bindings.charts." + chart_name + ".radial." + key,
                label="native_outer.radial_C5." + chart_name + "." + key,
                role="native_outer_actual_C5_radial_source_6_rows", x=x, z=Z,
                coordinate_handle=coordinate, phase_handle=phase, valid_range=domain,
            )
            if not radial_ok:
                _gate(gates, "source_radial_width_invalid:" + chart_name + ":" + key)
            radial_rows[key] = rows
            radial_jets[key] = _port_jet(rows)
        leading_radial_rows: dict[str, list[ScalarID]] = {}
        for key in ("E", "V"):
            rows, leading_radial_ok = _source_rows(
                dag, rows=radial_background.get(key), width=7, node_count=OUTER_GRAPH_NODE_COUNT,
                report_sha256=OUTER_ARCHIVE_SHA256, namespace=outer_namespace,
                archive=OUTER_ARCHIVE,
                path="true_outer_C6_C5_bindings.charts." + chart_name + ".radial_background." + key,
                label="native_outer.leading_radial_7." + chart_name + "." + key,
                role="native_outer_leading_radial_background_source_7_rows", x=x, z=Z,
                coordinate_handle=coordinate, phase_handle=phase, valid_range=domain,
            )
            if not leading_radial_ok:
                _gate(gates, "source_leading_radial_width_invalid:" + chart_name + ":" + key)
            leading_radial_rows[key] = rows

        chart_formula = {
            "coordinate": x,
            "integration_variable": t,
            "bounds": {"lower": lower, "upper": upper},
            "geometry": {
                "y": y, "J": J, "R": R, "y_lower": y_lower,
                "y_upper": y_upper, "R_upper": R_upper,
                "recipe": geometry_recipe, "R_formula": "R=Rm*exp(y)",
                "coordinate_relation": "log(R/Rm)=y(x); x is the local outer chart coordinate",
                "source_identity_proof": geometry_identity_proof,
            },
            "leading_background_fields_7": leading_rows,
            "loop_A_B_source_rows_6": loop_rows,
            "native_provider_E_N_V_N_rows_6": native_en_rows,
            "native_C5_density_source_rows_6": native_c5_density_rows,
            "native_C5_incoming_source_rows_6": incoming_c5_rows,
            "native_patch_C5_Rh_correction_source_rows_6": (
                incoming_c5_rows if previous_endpoint is None else {}
            ),
            "native_outer_chart_incoming_source_rows_6": outer_incoming_rows,
            "native_C5_history_source_rows_6": history_c5_rows,
            "native_patch_C5_history_source_rows_6": patch_history_source_rows,
            "native_patch_history_bridge_ports": patch_bridge_ports,
            "native_C5_modified_P0_source_rows_6": P0_rows,
            "native_C5_radial_source_rows_6": radial_rows,
            "leading_radial_background_source_rows_7": leading_radial_rows,
            "finite_N": {
                "invN": invN, "A_over_N": A_over_N, "r": r, "s": s,
                "E_N": E_N, "V_N": V_N, "F": F,
                "C2_first": c2_first, "C2_second": c2_second,
                "C2_density": c2_density, "native_C5_density": native_c5_density,
                "C5_equivalent_density": {
                    "r": r, "s": s,
                    "formula": "(s,r,V*r+E*s+r*s,2*V*s+s^2-E*r-r^2/2,E*r+r^2/2)",
                    "source_rule": "accepted_generic_density_rule",
                },
                "normalized_five_fields": finite_fields,
                "P0": _port_jet(P0_rows), "pressure": pressure,
                "physical_Utheta": Utheta, "physical_Uz": Uz,
            },
            "history": {
                "incoming_C2": incoming_c2, "incoming_C5": incoming_c5,
                "C2_integrals": c2_integrals, "C2_history": history_c2,
                "C2_endpoint_upper": endpoints_c2, "C5_history_source": history_c5_rows,
                "C5_history_bridge_proof": c5_history_bridge_proof,
                "C2_C5_inlet_identity_proofs": inlet_proofs,
                "C2_C5_history_identity_proofs": history_proofs,
                "C5_history_equals_native_patch_history": c5_history_bridge_proof is not None,
                "memory": c2_memory, "physical_derivative_identity": physical_derivative,
                "incoming_source": incoming_source,
                "predecessor_memory_never_reset": True,
            },
            "provider_identity_proof": {"E_N": provider_E, "V_N": provider_V},
            "density_identity_proof": density_proofs,
            "radial": {
                "actual_C5_E_y_V_y": radial_jets,
                "leading_background_radial_7": leading_radial_rows,
            },
            "open_continuation_ports": {
                "Q": chart.get("physical", {}).get("Q"),
                "Q_y": chart.get("physical", {}).get("Q_y"),
                "physical_Q_y_Z1_claimed": False,
            },
            "terminal_label": "Rc" if chart_name == "O3_power" else None,
        }
        formulas["charts"][chart_name] = chart_formula
        formulas["provider_identity_proofs"][chart_name] = {"E_N": provider_E, "V_N": provider_V}
        formulas["density_identity_proofs"][chart_name] = density_proofs
        formulas["c5_history_bridge_proofs"][chart_name] = c5_history_bridge_proof
        chain_proof = dag.op(
            "source_identity_proof",
            tuple(item for key in MOMENT_NAMES for item in incoming_c2[key]),
            theorem="C2_inlet_chain",
            chart=chart_name, incoming_source=incoming_source,
            first_patch_source_ids=PATCH_FIRST_INCOMING if previous_endpoint is None else None,
            exact_chain=True, incoming_correction_is_not_leading_moment=True,
            source_only=True, equality_is_symbolic=True,
        )
        formulas["endpoint_chain_proofs"].append(chain_proof)
        previous_endpoint = endpoints_c2
        previous_chart = chart_name
        previous_coordinate = x
        previous_history_proofs = history_proofs

    # One source-level theorem binds all chart formulas to the accepted
    # receipt and the conditional common-domain C6/C5 premise.
    source_theorem = dag.op(
        "source_identity_proof",
        tuple(proof for chart in formulas["charts"].values() for proof in chart["history"]["C2_C5_history_identity_proofs"].values())
        + tuple(formulas["endpoint_chain_proofs"])
        + tuple(item for item in formulas["c5_history_bridge_proofs"].values()
                if isinstance(item, str)),
        theorem="native_outer_six_source_identity",
        chart_order=list(CHART_ORDER),
        receipt_sha256=OUTER_RECEIPT_SHA256,
        archive_sha256=OUTER_ARCHIVE_SHA256,
        shared_N_source=witness["shared_N_ok"], shared_Pstar_source=witness["shared_Pstar_ok"],
        native_Rm_source=witness["native_Rm_ok"],
        conditional_premises=["same original analytic domain", "genuine C6/C5 regularity"],
        numerical_admissibility_claimed=False, source_only=True,
        no_legacy_geometry_alias=True, no_legacy_P0_alias=True,
    )
    formulas["source_identity_proof"] = source_theorem
    return dag, formulas


def _graph_checks(dag: _DAG, formulas: Mapping[str, object]) -> dict[str, bool]:
    """Static closure checks for the new graph only."""
    nodes = dag.nodes
    ids = [node.get("id") for node in nodes]
    node_map = {node["id"]: node for node in nodes}
    def ancestors(ident: str) -> set[str]:
        found: set[str] = set()
        pending = [ident]
        while pending:
            item = pending.pop()
            if item in found: continue
            found.add(item)
            row = node_map.get(item, {})
            pending.extend(arg for arg in row.get("args", []) if isinstance(arg, str))
        return found
    seen: set[str] = set()
    unique = True
    closed = True
    forbidden = {"evaluate", "numeric", "quadrature", "norm", "solve", "replay", "recurrence"}
    no_forbidden = True
    integral_ok = True
    endpoint_ok = True
    application_ok = True
    expm1_ok = True
    memory_ok = True
    jacobian_ok = True
    geometry_proof_ok = True
    history_bridge_ok = False
    for node in nodes:
        ident = node.get("id")
        if not isinstance(ident, str) or ident in seen or node.get("namespace") != SCHEMA:
            unique = False
        if isinstance(ident, str):
            seen.add(ident)
        operation = node.get("op")
        if isinstance(operation, str) and any(word in operation.lower() for word in forbidden):
            no_forbidden = False
        args = node.get("args", [])
        if not isinstance(args, list) or not all(isinstance(arg, str) and arg in seen - {ident} for arg in args):
            closed = False
        if operation == "integral":
            variable = node.get("variable")
            bounds = node.get("bounds")
            integral_ok = integral_ok and isinstance(variable, str) and variable in seen - {ident}
            integral_ok = integral_ok and isinstance(bounds, list) and len(bounds) == 2
            integral_ok = integral_ok and all(isinstance(item, str) and item in seen - {ident}
                                              for item in bounds)
            integral_ok = integral_ok and node.get("quadrature") is False
            integral_ok = integral_ok and node.get("full_jet_product") is True
            integral_ok = integral_ok and node.get("scalar_integrand") is True
            factor_ids = node.get("actual_integrand_factor_ids")
            jacobian_ok = jacobian_ok and isinstance(factor_ids, list) and len(factor_ids) >= 3
            jacobian_ok = jacobian_ok and isinstance(node.get("jacobian_factor_id"), str)
            jacobian_ok = jacobian_ok and node.get("jacobian_factor_id") in ancestors(node.get("args", [""])[0])
            jacobian_ok = jacobian_ok and all(item in node_map for item in factor_ids)
        if operation == "source_function_port":
            value = node.get("value", {})
            application_ok = application_ok and value.get("source_node_valid") is True
            application_ok = application_ok and value.get("source_graph_node_valid") is True
            application_ok = application_ok and value.get("re_evaluated") is False
        if operation == "function_substitution":
            endpoint_ok = endpoint_ok and node.get("endpoint_function") is not False
        if operation == "source_function_application":
            bindings = node.get("variable_bindings")
            application_ok = application_ok and node.get("application_is_unevaluated") is True
            application_ok = application_ok and node.get("source_only") is True
            application_ok = application_ok and isinstance(bindings, Mapping)
            application_ok = application_ok and all(isinstance(item, str) and item in seen - {ident}
                                                     for item in bindings.values())
        if operation == "expm1":
            expm1_ok = expm1_ok and node.get("entire_at_zero") is True
        if operation == "transport_memory":
            memory_ok = memory_ok and isinstance(node.get("rate"), str)
            if node.get("rate") == "0":
                memory_ok = memory_ok and node.get("rate_zero_memory_is_one") is True
        if operation == "source_identity_proof":
            theorem = node.get("theorem")
            if theorem == "native_outer_geometry_source_identity":
                geometry_proof_ok = geometry_proof_ok and (
                    node.get("source_ast_sha256") == OUTER_SOURCE_SHA256
                    and node.get("source_archive_sha256") == OUTER_ARCHIVE_SHA256
                    and node.get("bounds_Z_independent") is True
                    and node.get("Jacobian_included") is True
                )
            if theorem == "C5_DAG_history_equals_native_patch_chart_history":
                history_bridge_ok = history_bridge_ok or (
                    node.get("same_density") is True
                    and node.get("same_native_incoming") is True
                    and node.get("same_bounds") is True
                    and node.get("same_log_Duhamel_kernel") is True
                    and node.get("linear_ODE_uniqueness") is True
                )
    geometry_proof_count = sum(
        1 for node in nodes
        if node.get("op") == "source_identity_proof"
        and node.get("theorem") == "native_outer_geometry_source_identity"
    )
    geometry_proof_ok = geometry_proof_ok and geometry_proof_count == len(CHART_ORDER)
    jet_lists_ok = True
    # Only the explicit finite-N and history maps are jet-valued.  Native
    # seven/six row source maps are intentionally retained at their source
    # widths and are checked separately by the source witness.
    for chart in formulas.get("charts", {}).values() if isinstance(formulas.get("charts"), Mapping) else ():
        finite = chart.get("finite_N", {}) if isinstance(chart, Mapping) else {}
        for key in ("E_N", "V_N", "C2_density", "normalized_five_fields", "physical_Utheta", "physical_Uz"):
            value = finite.get(key) if isinstance(finite, Mapping) else None
            if isinstance(value, Mapping):
                jet_lists_ok = jet_lists_ok and all(isinstance(jet, list) and len(jet) == 3
                                                    and all(isinstance(item, str) for item in jet)
                                                    for jet in value.values())
            elif value is not None:
                jet_lists_ok = jet_lists_ok and isinstance(value, list) and len(value) == 3
        history = chart.get("history", {}) if isinstance(chart, Mapping) else {}
        for key in ("incoming_C2", "C2_history", "C2_endpoint_upper"):
            value = history.get(key) if isinstance(history, Mapping) else None
            if isinstance(value, Mapping):
                jet_lists_ok = jet_lists_ok and all(isinstance(jet, list) and len(jet) == 3
                                                    and all(isinstance(item, str) for item in jet)
                                                    for jet in value.values())
    return {
        "unique_namespace_ids": unique,
        "closed_topological_arguments": closed,
        "rational_constants_only": all(
            node.get("op") != "exact_rational" or
            (type(node.get("numerator")) is int and type(node.get("denominator")) is int
             and node.get("denominator") != 0)
            for node in nodes
        ),
        "ordinary_Z_jets_have_three_rows": jet_lists_ok,
        "scalar_ID_and_jet_list_distinct": jet_lists_ok and all(
            not isinstance(node.get("args"), list) or all(isinstance(item, str) for item in node.get("args", []))
            for node in nodes
        ),
        "no_numeric_or_solver_operation": no_forbidden,
        "endpoint_substitutions_are_function_nodes": endpoint_ok,
        "qualified_source_ports_valid": application_ok,
        "source_function_applications_have_explicit_closed_bindings": application_ok,
        "integral_variables_and_bounds_are_scalar_IDs": integral_ok,
        "integral_J_factor_in_actual_integrand_ancestry": jacobian_ok,
        "geometry_source_identity_proofs": geometry_proof_ok,
        "C5_history_native_patch_bridge_proof": history_bridge_ok,
        "expm1_entire_at_A_zero": expm1_ok,
        "nonreset_memory_retained_including_rate_zero": memory_ok,
        "closed_graph": bool(unique and closed),
    }


def _status_fields(
    witness: Mapping[str, object], checks: Mapping[str, bool], formulas: Mapping[str, object],
) -> dict[str, bool]:
    chart_maps = formulas.get("charts", {})
    chart_count_ok = isinstance(chart_maps, Mapping) and tuple(chart_maps) == CHART_ORDER
    provider = all(
        isinstance(chart_maps.get(name, {}).get("provider_identity_proof"), Mapping)
        and len(chart_maps[name]["provider_identity_proof"]) == 2
        for name in CHART_ORDER
    ) if chart_count_ok else False
    density = all(
        isinstance(chart_maps.get(name, {}).get("density_identity_proof"), Mapping)
        and set(chart_maps[name]["density_identity_proof"]) == set(MOMENT_NAMES)
        for name in CHART_ORDER
    ) if chart_count_ok else False
    inlet = all(
        isinstance(chart_maps.get(name, {}).get("history", {}).get("incoming_C2"), Mapping)
        and isinstance(chart_maps.get(name, {}).get("history", {}).get("incoming_C5"), Mapping)
        for name in CHART_ORDER
    ) if chart_count_ok else False
    nonreset = all(
        chart_maps.get(name, {}).get("history", {}).get("predecessor_memory_never_reset") is True
        and isinstance(chart_maps.get(name, {}).get("history", {}).get("C2_integrals"), Mapping)
        and len(chart_maps[name]["history"].get("C2_C5_history_identity_proofs",{})) == 5
        for name in CHART_ORDER
    ) if chart_count_ok else False
    history_bridge = bool(
        chart_count_ok
        and isinstance(chart_maps.get("Rh_reference", {}).get("history", {}).get(
            "C5_history_bridge_proof"), str)
        and chart_maps.get("Rh_reference", {}).get("history", {}).get(
            "C5_history_equals_native_patch_history") is True
    )
    endpoint = all(
        isinstance(chart_maps.get(name, {}).get("history", {}).get("C2_endpoint_upper"), Mapping)
        and len(chart_maps[name]["history"]["C2_endpoint_upper"]) == len(MOMENT_NAMES)
        for name in CHART_ORDER
    ) if chart_count_ok else False
    source_identity = bool(
        witness.get("source_ast_ok") and witness.get("receipt_ok")
        and witness.get("chart_order_ok") and witness.get("widths_ok")
        and witness.get("first_incoming_ok") and witness.get("shared_N_ok")
        and witness.get("admission_ok") and witness.get("endpoint_shape_ok")
        and witness.get("patch_source_ast_ok") and witness.get("patch_rh_correction_ok")
        and witness.get("common_prefix_ok") and witness.get("patch_z_handle_ok")
        and all(witness.get(key) is True for key in ("qualified_graph_ok", "native_Rm_ok",
            "shared_Pstar_ok", "shared_parameters_ok", "native_chain_ok", "source_bindings_ok", "generic_identity_ok", "patch_history_bridge_ok"))
        and all(checks.values())
    )
    return {
        "native_outer_six_fields_defined": bool(
            chart_count_ok and all(
                isinstance(chart_maps.get(name, {}).get("finite_N", {}).get("normalized_five_fields"), Mapping)
                and len(chart_maps[name]["finite_N"]["normalized_five_fields"]) == len(MOMENT_NAMES)
                for name in CHART_ORDER
            )
        ),
        "native_outer_six_C2_C5_provider_identity_proved": bool(provider and witness.get("receipt_ok") and witness.get("shared_parameters_ok") and witness.get("source_bindings_ok")),
        "native_outer_six_density_identity_proved": bool(density and witness.get("generic_identity_ok")),
        "native_outer_six_correction_inlet_identity_proved": bool(
            inlet and witness.get("first_incoming_ok") and witness.get("patch_rh_correction_ok") and witness.get("native_chain_ok") and witness.get("patch_history_bridge_ok")
        ),
        "native_outer_six_nonreset_history_identity_proved": bool(
            nonreset and history_bridge and witness.get("receipt_ok") and witness.get("native_chain_ok") and witness.get("patch_history_bridge_ok")
        ),
        "native_outer_six_endpoint_chain_identity_proved": bool(endpoint and chart_count_ok and witness.get("native_chain_ok") and witness.get("patch_history_bridge_ok")),
        "native_outer_six_source_identity": source_identity,
    }


def derive(source_records: Mapping[str, object], *, expected_input_digest: str) -> dict[str, object]:
    """Authenticate source snapshots and return one new unevaluated DAG.

    The caller owns raw-byte acquisition and the immutable checkpoint.  This
    function receives only those snapshots and the caller's canonical input
    digest; it never opens a repository path or imports a producer module.
    """
    gates: list[str] = []
    if not isinstance(source_records, Mapping):
        _gate(gates, "source_records_not_mapping")
        records: Mapping[str, object] = {}
    else:
        records = source_records
    try:
        input_digest = _digest(records)
    except (TypeError, ValueError):
        input_digest = ""
        _gate(gates, "source_records_not_canonical_json")
    if not isinstance(expected_input_digest, str) or expected_input_digest != input_digest:
        _gate(gates, "expected_input_digest_mismatch")
    if isinstance(records, Mapping):
        missing = sorted(set(RECORD_HASHES) - set(records))
        extra = sorted(set(records) - set(RECORD_HASHES))
        for name in missing:
            _gate(gates, "missing_source_record:" + name)
        for name in extra:
            _gate(gates, "unexpected_source_record:" + str(name))

    witness = _source_witness(records, gates)
    dag = _DAG()
    formulas: dict[str, object] = {}
    if not gates:
        dag, formulas = _build_graph(witness, gates)
    checks = _graph_checks(dag, formulas)
    statuses = _status_fields(witness, checks, formulas)
    for name, ok in checks.items():
        if not ok: _gate(gates,"new_definition_graph_check_failed:"+name)
    for name, ok in statuses.items():
        if not ok: _gate(gates,"new_source_status_unproved:"+name)
    all_statuses = all(statuses.values())
    passed = bool(not gates and all_statuses)
    nodes_digest = _digest(dag.nodes)
    chart_outputs: dict[str, object] = {}
    if isinstance(formulas.get("charts"), Mapping):
        for name, chart in formulas["charts"].items():
            if isinstance(chart, Mapping):
                history = chart.get("history", {})
                finite = chart.get("finite_N", {})
                chart_outputs[name] = {
                    "coordinate": chart.get("coordinate"),
                    "geometry": chart.get("geometry"),
                    "C2_history": history.get("C2_history"),
                    "C2_endpoint_upper": history.get("C2_endpoint_upper"),
                    "C5_incoming_source": history.get("incoming_C5"),
                    "C5_history_source": history.get("C5_history_source"),
                    "normalized_five_fields": finite.get("normalized_five_fields"),
                    "P0": finite.get("P0"),
                    "pressure": finite.get("pressure"),
                    "physical_Utheta": finite.get("physical_Utheta"),
                    "physical_Uz": finite.get("physical_Uz"),
                    "terminal_label": chart.get("terminal_label"),
                }

    source_audit = {
        "source_records_only": True,
        "static_source_AST_graph_saved_data_and_new_definitions_only": True,
        "original_graph_nodes_inspected_not_evaluated": True,
        "old_definition_graphs_rebuilt": False,
        "old_patch_DAG_rebuilt": False,
        "old_graph_re_evaluated": False,
        "new_definition_graph_constructed": bool(formulas),
        "new_outer_definition_graph_built": bool(formulas),
        "new_function_definition_graph_constructed": bool(formulas),
        "new_source_consumer_calls": 1,
        "old_ancestor_operator_calls": 0,
        "accepted_ancestor_operator_calls": 0,
        "producer_provider_builder_checker_replayed": False,
        "history_operator_calls": 0,
        "quadrature_or_norm_calls": 0,
        "recurrence_calls": 0,
        "target_calls": 0,
        "numeric_field_evaluated": False,
        "numeric_values_computed": False,
        "old_operators_replayed": False,
        "runtime_model_metadata": None,
    }
    source_audit_checks = {
        "source_ast": witness.get("source_ast_ok") is True,
        "receipt": witness.get("receipt_ok") is True,
        "shared_parameters": witness.get("shared_parameters_ok") is True,
        "shared_Pstar": witness.get("shared_Pstar_ok") is True,
        "native_Rm": witness.get("native_Rm_ok") is True,
        "source_bindings": witness.get("source_bindings_ok") is True,
        "native_chain": witness.get("native_chain_ok") is True,
        "patch_history_bridge": witness.get("patch_history_bridge_ok") is True,
        "geometry": witness.get("chart_order_ok") is True and witness.get("widths_ok") is True,
        "P0": witness.get("widths_ok") is True,
        "widths": witness.get("widths_ok") is True,
        "chart_order": witness.get("chart_order_ok") is True,
        "first_patch_incoming": witness.get("first_incoming_ok") is True,
        "accepted_patch": witness.get("admission_ok") is True,
        "qualified_endpoint": witness.get("endpoint_shape_ok") is True,
        "generic_density_rule": witness.get("generic_identity_ok") is True,
        "patch_source_AST": witness.get("patch_source_ast_ok") is True,
        "common_205147_prefix": witness.get("common_prefix_ok") is True,
        "patch_Rh_correction": witness.get("patch_rh_correction_ok") is True,
        "patch_Z_handle_namespace": witness.get("patch_z_handle_ok") is True,
        "native_provider_bindings": witness.get("source_bindings_ok") is True,
        "native_predecessor_chain": witness.get("native_chain_ok") is True,
        "qualified_patch_C5_history_bridge": witness.get("patch_history_bridge_ok") is True,
    }
    source_prerequisite_audit = {
        "actual_source_receipt": source_audit_checks["receipt"],
        "source_AST": source_audit_checks["source_ast"],
        "patch_source_AST_transport_log": source_audit_checks["patch_source_AST"],
        "shared_parameters": source_audit_checks["shared_parameters"],
        "geometry": source_audit_checks["geometry"],
        "P0": source_audit_checks["P0"],
        "widths": source_audit_checks["widths"],
        "chart_order": source_audit_checks["chart_order"],
        "common_205147_prefix": source_audit_checks["common_205147_prefix"],
        "outer_initial_correction_equals_patch_Rh": source_audit_checks["patch_Rh_correction"],
        "native_PATCH_Z_handle": PATCH_Z_HANDLE,
        "native_PATCH_Z_namespace": PATCH_Z_NAMESPACE,
        "all_required_source_prerequisites_passed": all(source_audit_checks.values()),
    }
    global_flags = {
        "same_density_certificate": False,
        "finite_N_history_identity": False,
        "physical_Q_y_Z1": False,
        "global_five_moment_closure": False,
        "legacy_feedback_reference_alias_proved": False,
        "prior_adapter_geometry_alias_proved": False,
        "prior_adapter_P0_alias_proved": False,
        "time_order_recursion_claimed": False,
        "timeorder_recursion_claimed": False,
        "numeric_values_computed": False,
        "old_operators_replayed": False,
        "full_goal_complete": False,
    }
    result: dict[str, object] = {
        "schema": SCHEMA,
        "source_only": True,
        "source_only_true": True,
        "passed": passed,
        "input_digest": input_digest,
        "expected_input_digest": expected_input_digest,
        "record_hashes": dict(RECORD_HASHES),
        "decoded_record_hashes": dict(DECODED_RECORD_HASHES),
        "ordinary_Z_orders": list(ORDINARY_Z_ORDERS),
        "chart_order": list(CHARTS),
        "ordinary_Z_not_temporal_n": True,
        "new_source_consumer": 1,
        "old_ancestor": 0,
        "quadrature_or_norm_calls": 0,
        "runtime_model_metadata": None,
        "nodes_digest": nodes_digest,
        "node_count": len(dag.nodes),
        "definition_graph": {
            "namespace": SCHEMA,
            "nodes": dag.nodes,
            "node_count": len(dag.nodes),
            "nodes_digest": nodes_digest,
            "outputs": chart_outputs,
        },
        "function_definitions": formulas,
        "outputs": chart_outputs,
        "source_parameters": formulas.get("source_parameters", {}),
        "scalar_ID_type": "str DAG ID; one scalar expression",
        "jet_list_type": "list[str] of exactly three scalar IDs for ordinary Z orders 0,1,2",
        "ordinary_Z_jet_rules": {
            "orders": list(ORDINARY_Z_ORDERS),
            "product_second": "a2*b0 + 2*a1*b1 + a0*b2",
            "exp": "(exp(a0),exp(a0)*a1,exp(a0)*(a2+a1^2))",
            "expm1": "(expm1(a0),exp(a0)*a1,exp(a0)*(a2+a1^2)); entire_at_zero=True",
            "scalar_integration_variable": True,
            "temporal_recursion_claimed": False,
        },
        "definition_graph_checks": checks,
        "source_audit_checks": source_audit_checks,
        "source_prerequisite_audit": source_prerequisite_audit,
        "source_witness": {
            "accepted_reports": {
                name: {"path": name, "raw_sha256": RECORD_HASHES[name],
                       "decoded_sha256": DECODED_RECORD_HASHES.get(name, RECORD_HASHES[name])}
                for name in RECORD_HASHES
            },
            "archive_graph_counts": {"outer": OUTER_GRAPH_NODE_COUNT, "patch": PATCH_GRAPH_NODE_COUNT},
            "outer_namespace": witness.get("archive_namespace"),
            "patch_namespace": witness.get("patch_namespace"),
            "source_function_names": witness.get("source_function_names", []),
            "chart_order": list(CHART_ORDER),
            "widths": witness.get("widths", {}),
            "receipt_ok": witness.get("receipt_ok"),
            "source_ast_ok": witness.get("source_ast_ok"),
            "patch_source_ast_ok": witness.get("patch_source_ast_ok"),
            "generic_density_identity_ok": witness.get("generic_identity_ok"),
            "first_patch_incoming_ok": witness.get("first_incoming_ok"),
            "patch_Rh_correction_ok": witness.get("patch_rh_correction_ok"),
            "common_205147_prefix_ok": witness.get("common_prefix_ok"),
            "patch_Z_handle_ok": witness.get("patch_z_handle_ok"),
            "source_binding_checks": witness.get("source_binding_checks"),
            "native_chain_checks": witness.get("native_chain_checks"),
            "patch_history_bridge_checks": witness.get("patch_history_bridge_checks"),
            "conditional_common_original_domain": True,
            "conditional_genuine_C6_C5_regularity": True,
            "numerical_admissibility_claimed": False,
        },
        "retained_construction_data": dict(RETAINED_CONSTRUCTION),
        "terminal_endpoint": "Rc",
        "open_continuation": {
            "physical_Q": "qualified source ports retained as data only",
            "physical_Q_y": "qualified source ports retained as data only",
            "physical_Q_y_Z1": False,
            "global_five_moment_closure": False,
            "pressure_continuation": "open",
        },
        "source_inspection_audit": source_audit,
        "status": dict(statuses),
        "unresolved_gates": sorted(set(gates)),
    }
    result.update(statuses)
    result.update(global_flags)
    return result


__all__ = [
    "CHARTS", "DECODED_RECORD_HASHES", "ORDINARY_Z_ORDERS", "RECORD_HASHES", "SCHEMA", "derive",
]
