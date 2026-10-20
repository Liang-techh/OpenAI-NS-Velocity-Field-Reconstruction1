"""Source-only finite-``N`` consumer for the native C2 Rm patch.

The module authenticates saved source records, inspects their AST and saved
binding data, and builds a new unevaluated ordinary-Z definition graph.  It
does not import a project producer, execute a saved graph, solve a control
system, or evaluate a field.  Every source function is represented by a
qualified source port; the finite-N formulas below are new DAG nodes.
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


SCHEMA = "original_native_C2_patch_finite_N_source_v1"
EXPECTED_NODE_COUNT = 205300
EXPECTED_NODES_DIGEST = "a2dcc782ba3033ff6ca5f54a66632ac6ef007bb0c5b9960f10567d44ef2059a3"

LATEST_REPORT = "lei_ren_part1_paper_compliant_current_original_native_C2_patch_composition.json"
LATEST_REPORT_SHA256 = "f426f396b3b9f55aea1bf1b6c07cd7678addae7644093b97babb110a47619136"
PROVIDER_REPORT = "lei_ren_part1_paper_compliant_current_original_finite_N_provider_source_binding.json"
PROVIDER_REPORT_SHA256 = "1f23b6188abffb6d99cf00fb1274a58a2fba1100a326334a1a798d0b920f6c41"
MEMORY_REPORT = "lei_ren_part1_paper_compliant_current_original_finite_N_nonreset_memory_source_binding.json"
MEMORY_REPORT_SHA256 = "41defb67ef05207a079075328f3f7afd3980193c4cec7adb9574a27d41f33335"
GENERIC_REPORT = "lei_ren_part1_paper_compliant_current_original_finite_N_field_density_source_composition.json"
GENERIC_REPORT_SHA256 = "0e5d3002ccdc854a1f6be32f2ea076f5c7d650ae25e6a766c7161bca4cd90e2f"
C5_SOURCE = "lei_ren_part1_paper_compliant_current_original_true_Rm_patch_finite_N_C5_functions.py"
C5_SOURCE_SHA256 = "c4445ef8acfb33d29cf9022862f8b6396e6aa2f38bda7467409b5378a5deaf15"
C5_ARCHIVE = "lei_ren_part1_paper_compliant_current_original_true_Rm_patch_finite_N_C5_functions.json.gz"
C5_ARCHIVE_SHA256 = "8f3bc343e9dbffbdda465ce94d6183735184645a70e291cb900d3cde726769ad"
C5_ARCHIVE_DECODED_SHA256 = "c1026e05015584a7ce29bd7e8da95b17d97669523f1d063ac5b5c800559d9b69"
LOOP_SOURCE = "lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions.py"
LOOP_SOURCE_SHA256 = "69d6391120000488472b7961b4bff1901bff872def95d183c4ee056b0948e17d"
LOOP_ARCHIVE = "lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions.json.gz"
LOOP_ARCHIVE_SHA256 = "8419fdc338dc746099c7cad119ea61a8e928652319d5caf5314f8a730c533caa"
LOOP_ARCHIVE_DECODED_SHA256 = "f4825741f568c001f0f4658dfc3ce8eec3dc0466db976c29f2f4b609b8ceadc2"
LOOP_ARCHIVE_PREFIX = 121834
LOOP_ARCHIVE_PREFIX_SHA256 = "bc66540d371bcaf875acd3dcd3b926fd3dfaa3172368e0b3e65b93eaa3772525"
ADAPTER_REPORT = "lei_ren_part1_paper_compliant_current_original_native_C2_reference_adapter.json"
ADAPTER_REPORT_SHA256 = "3fd6359a456345f15f3323c594558c9edb89169a21114837c160be3ef86e3db7"
R100_C6_ARCHIVE = "lei_ren_part1_paper_compliant_current_original_prescribed_R100_Rm_C6_functions.json.gz"
R100_C6_ARCHIVE_SHA256 = "5e5b565dd8d2abc41e6f3df0d9e17000e524d493e58a5778dfcf217c6e42f44a"
R100_C6_ARCHIVE_DECODED_SHA256 = "18d06627b935a46de5eca87c286a9a51c5a2d5ce14f11ae0c1c06a7c6308221a"
R100_C6_GRAPH_COUNT = 97295
TRUE_C6_ARCHIVE = "lei_ren_part1_paper_compliant_current_original_true_Rm_patch_C6_functions.json.gz"
TRUE_C6_ARCHIVE_SHA256 = "4352ebd8862b5c92e3e4026728183014063baf556802a88421cdbd517ff1346a"
TRUE_C6_ARCHIVE_DECODED_SHA256 = "b5671cc0c60821db9a5c591ea899c9616c3cc949cec4e9214cc57a31ae6cdc10"
TRUE_C6_PREFIX = 196042
TRUE_C6_PREFIX_SHA256 = "eb5b2d5ee9a0332b3e6dce597b3971f39c50413b701f4c58773ae84ba1e75276"
PSTAR_SOURCE_HANDLE = 33572

# The runner must pass exactly these authenticated records.  They are source
# snapshots, not paths to be opened by this consumer.
RECORD_HASHES = {
    LATEST_REPORT: LATEST_REPORT_SHA256,
    PROVIDER_REPORT: PROVIDER_REPORT_SHA256,
    MEMORY_REPORT: MEMORY_REPORT_SHA256,
    GENERIC_REPORT: GENERIC_REPORT_SHA256,
    C5_SOURCE: C5_SOURCE_SHA256,
    C5_ARCHIVE: C5_ARCHIVE_SHA256,
    LOOP_SOURCE: LOOP_SOURCE_SHA256,
    LOOP_ARCHIVE: LOOP_ARCHIVE_SHA256,
    ADAPTER_REPORT: ADAPTER_REPORT_SHA256,
    R100_C6_ARCHIVE: R100_C6_ARCHIVE_SHA256,
    TRUE_C6_ARCHIVE: TRUE_C6_ARCHIVE_SHA256,
}
DECODED_RECORD_HASHES = {
    C5_ARCHIVE: C5_ARCHIVE_DECODED_SHA256,
    LOOP_ARCHIVE: LOOP_ARCHIVE_DECODED_SHA256,
    R100_C6_ARCHIVE: R100_C6_ARCHIVE_DECODED_SHA256,
    TRUE_C6_ARCHIVE: TRUE_C6_ARCHIVE_DECODED_SHA256,
}

RATES = {"m": "1", "h": "3/2", "k": "3/2", "e": "1", "p": "0"}
MOMENT_NAMES = tuple(RATES)
ORDINARY_Z_ORDERS = (0, 1, 2)
PATCH_ARCHIVE_PREFIX = 205147
PATCH_COORDINATE_HANDLE = 185676
PATCH_Z_HANDLE = 20
PATCH_PHASE_HANDLE = 196046
LATEST_X_NODE = "original_native_C2_patch_composition_v1:d000002"
LATEST_Z_NODE = "original_native_C2_patch_composition_v1:d000001"
PATCH_LOOP_A = (197890, 197891, 197892, 197893, 197894, 197895)
PATCH_LOOP_B = (197928, 197929, 197930, 197931, 197932, 197933)
PATCH_C2_INCOMING = {
    "m": (2891, 2892, 11305),
    "h": (2909, 2910, 11311),
    "k": (2921, 2922, 11317),
    "e": (2933, 2934, 11323),
    "p": (2943, 2944, 11328),
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

# Type aliases make the scalar-vs-jet distinction explicit.  A scalar ID is a
# single DAG ID; a JetList is exactly the three ordinary Z rows.
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


def _record(records: Mapping[str, object], name: str, gates: list[str]) -> str:
    row = records.get(name)
    if not isinstance(row, Mapping):
        _gate(gates, "missing_source_record:" + name)
        return ""
    expected_path = "experiments/root_st073/" + name
    if row.get("path") != expected_path:
        _gate(gates, "source_path_mismatch:" + name)
    if row.get("sha256") != RECORD_HASHES[name]:
        _gate(gates, "source_declared_hash_mismatch:" + name)
    source = row.get("source_text")
    if not isinstance(source, str):
        _gate(gates, "source_text_missing:" + name)
        return ""
    if name in DECODED_RECORD_HASHES:
        if row.get("encoding") != "gzip":
            _gate(gates, "encoded_archive_contract_mismatch:" + name)
        encoded = row.get("raw_bytes_base64")
        try:
            # Authenticate the raw bytes before any decompression or JSON parse.
            raw = base64.b64decode(encoded if isinstance(encoded, str) else "", validate=True)
        except Exception as error:  # malformed input remains an unresolved gate
            _gate(gates, "archive_raw_base64_error:" + name + ":" + type(error).__name__)
            return ""
        if hashlib.sha256(raw).hexdigest() != RECORD_HASHES[name]:
            _gate(gates, "archive_raw_hash_mismatch:" + name)
            return ""
        try:
            decoded = gzip.decompress(raw)
        except Exception as error:
            _gate(gates, "archive_decode_error:" + name + ":" + type(error).__name__)
            return ""
        decoded_hash = hashlib.sha256(decoded).hexdigest()
        if row.get("decoded_sha256") != DECODED_RECORD_HASHES[name] or decoded_hash != DECODED_RECORD_HASHES[name]:
            _gate(gates, "archive_decoded_hash_mismatch:" + name)
        if decoded != source.encode("utf-8"):
            _gate(gates, "archive_source_text_mismatch:" + name)
        return source
    if row.get("encoding") is not None or row.get("raw_bytes_base64") is not None:
        _gate(gates, "unexpected_encoded_record:" + name)
    if hashlib.sha256(source.encode("utf-8")).hexdigest() != RECORD_HASHES[name]:
        _gate(gates, "source_hash_mismatch:" + name)
        return ""
    return source


def _json(text: str, name: str, gates: list[str]) -> Mapping[str, object]:
    if not text:
        return {}
    try:
        value = json.loads(text)
    except Exception as error:
        _gate(gates, "json_parse_error:" + name + ":" + type(error).__name__)
        return {}
    if not isinstance(value, Mapping):
        _gate(gates, "json_root_not_mapping:" + name)
        return {}
    return value


def _ast(text: str, name: str, gates: list[str]) -> ast.AST | None:
    if not text:
        return None
    try:
        tree = ast.parse(text, filename=name)
        # Compile and symbol-table construction are static checks only.  Never
        # execute the resulting code object.
        compile(tree, name, "exec")
        symtable.symtable(text, name, "exec")
        return tree
    except (SyntaxError, ValueError, TypeError) as error:
        _gate(gates, "source_static_parse_error:" + name + ":" + type(error).__name__)
        return None


def _function_names(tree: ast.AST | None) -> set[str]:
    if tree is None:
        return set()
    return {node.name for node in ast.walk(tree)
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef))}


def _has(text: str, *needles: str) -> bool:
    return all(needle in text for needle in needles)


def _valid_handle(value: object, count: int) -> bool:
    return type(value) is int and 0 <= value < count


def _valid_row(value: object, count: int, width: int | None = None) -> bool:
    if not isinstance(value, list):
        return False
    if width is not None and len(value) != width:
        return False
    return all(_valid_handle(item, count) for item in value)


def _valid_jet(value: object, ids: set[str]) -> bool:
    return isinstance(value, list) and len(value) == 3 and all(item in ids for item in value)


class _DAG:
    """Small append-only graph whose operation arguments are local IDs."""

    def __init__(self) -> None:
        self.nodes: list[dict[str, object]] = []
        self._counter = 0

    def _id(self) -> str:
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
    if isinstance(value, Fraction):
        return value
    return Fraction(str(value))


def _const(dag: _DAG, value: object) -> ScalarID:
    rational = _fraction(value)
    return dag.atom("exact_rational", numerator=rational.numerator,
                    denominator=rational.denominator)


def _zero(dag: _DAG) -> ScalarID:
    return _const(dag, 0)


def _one(dag: _DAG) -> ScalarID:
    return _const(dag, 1)


def _symbol(dag: _DAG, name: str, domain: str) -> ScalarID:
    return dag.atom("exact_symbol", name, domain=domain)


def _source_port(
    dag: _DAG,
    *,
    report_sha256: str,
    namespace: str,
    source_node: object,
    archive: str,
    label: str,
    role: str,
    source_path: str,
    valid: bool,
    coordinate_handle: object = None,
    phase_handle: object = None,
    bound_variable: str | None = None,
    coordinate_binding: ScalarID | None = None,
    z_binding: ScalarID | None = None,
    source_coordinate_handle: object = None,
    source_z_handle: object = None,
    source_coordinate_name: str | None = None,
    source_z_name: str | None = None,
    **metadata: object,
) -> ScalarID:
    value: dict[str, object] = {
        "source_report_sha256": report_sha256,
        "source_namespace": namespace,
        "source_node": source_node if valid else "missing",
        "source_archive": archive,
        "archive_sha256": RECORD_HASHES.get(archive),
        "label": label,
        "role": role,
        "source_path": source_path,
        "source_node_valid": valid,
        "source_graph_node_valid": valid,
        "qualified": (SCHEMA + "|report=" + report_sha256 + "|namespace=" + namespace
                       + "|node=" + str(source_node if valid else "missing")
                       + "|path=" + source_path + "|label=" + label),
        "re_evaluated": False,
    }
    if coordinate_handle is not None:
        value["coordinate_handle"] = coordinate_handle
    if phase_handle is not None:
        value["phase_handle"] = phase_handle
    if bound_variable is not None:
        value["bound_variable"] = bound_variable
    if source_coordinate_handle is not None:
        value["source_coordinate_handle"] = source_coordinate_handle
    if source_z_handle is not None:
        value["source_z_handle"] = source_z_handle
    value.update(metadata)
    port = dag.atom("source_function_port", value)
    bindings: list[tuple[str, ScalarID]] = []
    if coordinate_binding is not None:
        bindings.append((source_coordinate_name or bound_variable or "x", coordinate_binding))
    if z_binding is not None:
        bindings.append((source_z_name or "Z", z_binding))
    if not bindings:
        return port
    return dag.op(
        "source_function_application", tuple([port] + [node for _, node in bindings]),
        exact=True, source_port=port,
        variable_bindings={name: node for name, node in bindings},
        source_coordinate_handle=source_coordinate_handle,
        source_z_handle=source_z_handle,
        application_is_unevaluated=True,
        source_only=True,
    )


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


def _log(dag: _DAG, value: ScalarID) -> ScalarID:
    return _op(dag, "log", value, analytic=True, positive_argument=True)


def _power(dag: _DAG, value: ScalarID, exponent: object) -> ScalarID:
    return _op(dag, "power", value, _const(dag, exponent),
               exponent=str(_fraction(exponent)), positive_base=True)


def _jet_check(jet: object) -> JetList:
    if not isinstance(jet, list) or len(jet) != 3 or not all(isinstance(x, str) for x in jet):
        raise TypeError("ordinary Z jet must be a list of three scalar IDs")
    return jet


def _jet_source(source: ScalarID, dag: _DAG) -> JetList:
    return [source, _zero(dag), _zero(dag)]


def _jet_const(dag: _DAG, value: object) -> JetList:
    return [_const(dag, value), _zero(dag), _zero(dag)]


def _jet_scalar(dag: _DAG, scalar: ScalarID) -> JetList:
    return [scalar, _zero(dag), _zero(dag)]


def _jet_add(dag: _DAG, *jets: JetList) -> JetList:
    checked = [_jet_check(jet) for jet in jets]
    if not checked:
        return _jet_const(dag, 0)
    return [_add(dag, *(jet[index] for jet in checked)) for index in ORDINARY_Z_ORDERS]


def _jet_sub(dag: _DAG, left: JetList, right: JetList) -> JetList:
    left, right = _jet_check(left), _jet_check(right)
    return [_sub(dag, left[i], right[i]) for i in ORDINARY_Z_ORDERS]


def _jet_scale(dag: _DAG, jet: JetList, scalar: ScalarID | JetList) -> JetList:
    jet = _jet_check(jet)
    if isinstance(scalar, list):
        return _jet_mul(dag, jet, scalar)
    if not isinstance(scalar, str):
        raise TypeError("jet scale requires scalar ID or jet list")
    return [_mul(dag, scalar, jet[i]) for i in ORDINARY_Z_ORDERS]


def _jet_mul(dag: _DAG, left: JetList, right: JetList) -> JetList:
    left, right = _jet_check(left), _jet_check(right)
    return [
        _mul(dag, left[0], right[0]),
        _add(dag, _mul(dag, left[1], right[0]), _mul(dag, left[0], right[1])),
        _add(dag, _mul(dag, left[2], right[0]),
             _mul(dag, _const(dag, 2), left[1], right[1]),
             _mul(dag, left[0], right[2])),
    ]


def _jet_inv(dag: _DAG, jet: JetList) -> JetList:
    jet = _jet_check(jet)
    a0, a1, a2 = jet
    a0sq = _mul(dag, a0, a0)
    a0cu = _mul(dag, a0sq, a0)
    return [_div(dag, _one(dag), a0),
            _neg(dag, _div(dag, a1, a0sq)),
            _div(dag, _sub(dag, _mul(dag, _const(dag, 2), a1, a1),
                           _mul(dag, a0, a2)), a0cu)]


def _jet_div(dag: _DAG, left: JetList, right: JetList) -> JetList:
    return _jet_mul(dag, left, _jet_inv(dag, right))


def _jet_exp(dag: _DAG, jet: JetList) -> JetList:
    jet = _jet_check(jet)
    e0 = _exp(dag, jet[0])
    return [e0, _mul(dag, e0, jet[1]),
            _mul(dag, e0, _add(dag, jet[2], _mul(dag, jet[1], jet[1])))]


def _jet_expm1(dag: _DAG, jet: JetList) -> JetList:
    """Ordinary jet of expm1 with no division by its argument."""
    jet = _jet_check(jet)
    e0 = _exp(dag, jet[0])
    return [_expm1(dag, jet[0]), _mul(dag, e0, jet[1]),
            _mul(dag, e0, _add(dag, jet[2], _mul(dag, jet[1], jet[1])))]


def _integral(dag: _DAG, scalar: ScalarID, variable: ScalarID,
              lower: ScalarID, upper: ScalarID, **metadata: object) -> ScalarID:
    return dag.op("integral", (scalar,), variable=variable, bounds=[lower, upper],
                  exact=True, unevaluated=True, quadrature=False, **metadata)


def _function_substitution(dag: _DAG, scalar: ScalarID, variable: ScalarID,
                           value: ScalarID, **metadata: object) -> ScalarID:
    return dag.op("function_substitution", (scalar, variable, value), exact=True,
                  **metadata)


def _substitute_jet(dag: _DAG, jet: JetList, variable: ScalarID,
                    value: ScalarID, **metadata: object) -> JetList:
    jet = _jet_check(jet)
    return [_function_substitution(dag, row, variable, value,
                                   Z_independent_substitution=True,
                                   **metadata) for row in jet]


def _integral_jet(dag: _DAG, jet: JetList, variable: ScalarID,
                  lower: ScalarID, upper: ScalarID, **metadata: object) -> JetList:
    """Integrate a complete product jet row-by-row after multiplication."""
    jet = _jet_check(jet)
    return [_integral(dag, jet[index], variable, lower, upper,
                      jet_component=index, full_jet_product=True,
                      ordinary_Z_order=index, **metadata)
            for index in ORDINARY_Z_ORDERS]


def _safe_source_ports(
    dag: _DAG,
    rows: object,
    *,
    report_sha256: str,
    namespace: str,
    archive: str,
    label: str,
    path: str,
    valid_count: int,
    gates: list[str],
    role: str,
    coordinate_handle: object = None,
    phase_handle: object = None,
    bound_variable: str | None = None,
    **metadata: object,
) -> JetList:
    valid = _valid_row(rows, valid_count, 3)
    if not valid:
        _gate(gates, "invalid_source_jet:" + path)
    values = list(rows)[:3] if isinstance(rows, list) else ["missing"] * 3
    while len(values) < 3:
        values.append("missing")
    return [
        _source_port(dag, report_sha256=report_sha256, namespace=namespace,
                     source_node=value, archive=archive,
                     label=label + "[" + str(index) + "]", role=role,
                     source_path=path + "[" + str(index) + "]",
                     valid=valid and _valid_handle(value, valid_count),
                     coordinate_handle=coordinate_handle, phase_handle=phase_handle,
                     bound_variable=bound_variable, **metadata)
        for index, value in enumerate(values)
    ]


def _latest_graph_ids(latest: Mapping[str, object]) -> set[str]:
    graph = latest.get("component_result", {}).get("definition_graph", {}) if isinstance(latest.get("component_result"), Mapping) else {}
    nodes = graph.get("nodes", []) if isinstance(graph, Mapping) else []
    return {node.get("id") for node in nodes if isinstance(node, Mapping) and isinstance(node.get("id"), str)}


def _saved_latest_jet(
    dag: _DAG, latest: Mapping[str, object], name: str, *, gates: list[str], path: str,
    section: str = "normalized_fields", coordinate_binding: ScalarID | None = None,
    z_binding: ScalarID | None = None,
) -> JetList:
    component = latest.get("component_result", {})
    fields = component.get("function_definitions", {}) if isinstance(component, Mapping) else {}
    group = fields.get(section, {}) if isinstance(fields, Mapping) else None
    row = group.get(name) if isinstance(group, Mapping) else None
    ids = _latest_graph_ids(latest)
    valid = _valid_jet(row, ids)
    if not valid:
        _gate(gates, "invalid_previous_definition_jet:" + path)
    values = list(row)[:3] if isinstance(row, list) else ["missing"] * 3
    while len(values) < 3:
        values.append("missing")
    return [_source_port(
        dag, report_sha256=LATEST_REPORT_SHA256,
        namespace="original_native_C2_patch_composition_v1",
        source_node=value, archive=LATEST_REPORT, label="previous." + name + "[" + str(index) + "]",
        role="qualified_previous_DAG_reference", source_path=path + "[" + str(index) + "]",
        valid=valid and value in ids, source_graph_digest=EXPECTED_NODES_DIGEST,
        source_graph_node_valid=valid and value in ids,
        coordinate_binding=coordinate_binding, z_binding=z_binding,
        source_coordinate_handle=LATEST_X_NODE, source_z_handle=LATEST_Z_NODE,
        source_coordinate_name="x", source_z_name="Z",
    ) for index, value in enumerate(values)]


def _patch_archive_witness(archive: Mapping[str, object], gates: list[str]) -> dict[str, object]:
    graph = archive.get("exact_graph_nodes", [])
    bindings = archive.get("Rm_patch_finite_N_C5_bindings", {})
    chart = bindings.get("chart", {}) if isinstance(bindings, Mapping) else {}
    finite = chart.get("finite_N", {}) if isinstance(chart, Mapping) else {}
    loop = chart.get("loop", {}) if isinstance(chart, Mapping) else {}
    incoming = bindings.get("incoming", {}) if isinstance(bindings, Mapping) else {}
    history = chart.get("history", {}) if isinstance(chart, Mapping) else {}
    shape = {
        "graph_count": len(graph) if isinstance(graph, list) else None,
        "graph_count_exact": isinstance(graph, list) and len(graph) == PATCH_ARCHIVE_PREFIX,
        "coordinate": chart.get("coordinate"),
        "phase": chart.get("actual_phase"),
        "radius": chart.get("radius"),
        "physical_measure": chart.get("physical_measure"),
        "loop_A": loop.get("A"),
        "loop_B_over_Pstar": loop.get("B_over_Pstar"),
        "finite_N": finite,
        "incoming": incoming,
        "history": history,
        "source_family": archive.get("source_family"),
        "definition": bindings.get("definition", {}),
    }
    if not shape["graph_count_exact"]:
        _gate(gates, "patch_archive_prefix_mismatch")
    if chart.get("coordinate") != PATCH_COORDINATE_HANDLE:
        _gate(gates, "patch_coordinate_handle_mismatch")
    if chart.get("actual_phase") != PATCH_PHASE_HANDLE:
        _gate(gates, "patch_phase_handle_mismatch")
    if loop.get("A") != list(PATCH_LOOP_A):
        _gate(gates, "patch_loop_A_handle_mismatch")
    if loop.get("B_over_Pstar") != list(PATCH_LOOP_B):
        _gate(gates, "patch_loop_B_handle_mismatch")
    for key, row in (("incoming", incoming), ("history", history)):
        if not isinstance(row, Mapping) or set(row) != set(MOMENT_NAMES):
            _gate(gates, "patch_" + key + "_moment_rows_missing")
        else:
            for moment in MOMENT_NAMES:
                if not _valid_row(row[moment], PATCH_ARCHIVE_PREFIX, 6):
                    _gate(gates, "patch_" + key + "_row_invalid:" + moment)
    if not isinstance(finite, Mapping) or not all(
        _valid_row(finite.get(key), PATCH_ARCHIVE_PREFIX, 6)
        for key in ("A_over_N", "E_N", "V_N", "delta_E", "delta_V")
    ):
        _gate(gates, "patch_finite_N_binding_rows_invalid")
    densities = finite.get("five_signed_density_increments", {}) if isinstance(finite, Mapping) else {}
    if not isinstance(densities, Mapping) or set(densities) != set(MOMENT_NAMES) or not all(
        _valid_row(densities.get(key), PATCH_ARCHIVE_PREFIX, 6) for key in MOMENT_NAMES
    ):
        _gate(gates, "patch_finite_N_density_rows_invalid")
    if not isinstance(bindings.get("N"), int) or not _valid_handle(bindings.get("N"), PATCH_ARCHIVE_PREFIX):
        _gate(gates, "patch_N_source_parameter_invalid")
    if not isinstance(bindings.get("phase_origin"), int) or not _valid_handle(bindings.get("phase_origin"), PATCH_ARCHIVE_PREFIX):
        _gate(gates, "patch_phase_origin_source_parameter_invalid")
    selected = list(PATCH_LOOP_A) + list(PATCH_LOOP_B)
    for handle in selected:
        if _valid_handle(handle, len(graph)):
            node = graph[handle]
            if not isinstance(node, Mapping) or node.get("operation") != "original_lazy_flat_branch":
                _gate(gates, "patch_loop_source_operation_mismatch:" + str(handle))
        else:
            _gate(gates, "patch_loop_source_handle_invalid:" + str(handle))
    if isinstance(graph, list) and _valid_handle(PATCH_COORDINATE_HANDLE, len(graph)):
        coordinate = graph[PATCH_COORDINATE_HANDLE]
        if not (isinstance(coordinate, Mapping) and coordinate.get("operation") == "bound_variable"
                and coordinate.get("name") == "true_Rm_patch_x"):
            _gate(gates, "patch_coordinate_source_operation_mismatch")
    if isinstance(graph, list) and _valid_handle(PATCH_Z_HANDLE, len(graph)):
        z_source = graph[PATCH_Z_HANDLE]
        if not (isinstance(z_source, Mapping) and z_source.get("operation") == "bound_variable"
                and z_source.get("name") == "Z"):
            _gate(gates, "patch_Z_source_operation_mismatch")
    if isinstance(graph, list) and _valid_handle(PATCH_PHASE_HANDLE, len(graph)):
        phase = graph[PATCH_PHASE_HANDLE]
        if not (isinstance(phase, Mapping) and phase.get("operation") == "analytic_unary"
                and phase.get("name") == "fractional_part"):
            _gate(gates, "patch_phase_source_operation_mismatch")
    return shape


def _native_geometry_witness(graph: object, chart: object, inlet: object) -> dict[str, object]:
    """Inspect only the bounded saved geometry subgraph, never evaluate it."""
    if not isinstance(graph, list) or not isinstance(chart, Mapping) or not isinstance(inlet, Mapping):
        return {"valid": False, "closure": []}
    checks = {
        "native_Rm_source": inlet.get("radius") == 85743,
        "native_chart_ports": chart.get("radius") == 187214 and chart.get("physical_measure") == 187215,
        "native_radius_product": graph[187214] == {"operation": "product", "arguments": [85743, PATCH_COORDINATE_HANDLE]},
        "native_measure_quotient": graph[187215].get("operation") == "positive_quotient"
            and graph[187215].get("numerator") == 1
            and graph[187215].get("denominator") == PATCH_COORDINATE_HANDLE,
        "native_one": graph[1] == {"operation": "exact_rational", "numerator": 1, "denominator": 1},
    }
    pending = [85743, 187214, 187215]
    closure: set[int] = set()
    closure_valid = True
    while pending:
        ident = pending.pop()
        if type(ident) is not int or not 0 <= ident < len(graph):
            closure_valid = False
            continue
        if ident in closure:
            continue
        closure.add(ident)
        node = graph[ident]
        if not isinstance(node, Mapping):
            closure_valid = False
            continue
        operation = node.get("operation")
        refs: list[object] = []
        if operation in ("product", "sum"):
            arguments = node.get("arguments")
            if not isinstance(arguments, list):
                closure_valid = False
            else:
                refs = arguments
        elif operation == "positive_quotient":
            refs = [node.get("numerator"), node.get("denominator")]
        elif operation == "analytic_unary":
            closure_valid &= node.get("name") == "exp"
            refs = [node.get("argument")]
        elif operation == "function_substitution":
            closure_valid &= (node.get("Z_independent_substitution") is True
                              and type(node.get("value")) is int and node.get("value") == 1)
            refs = [node.get("expression"), node.get("variable")]
        elif operation == "bound_variable":
            closure_valid &= node.get("name") in ("true_prescribed_post_restore_fraction", "true_Rm_patch_x")
        elif operation == "exact_rational":
            closure_valid &= (type(node.get("numerator")) is int
                              and type(node.get("denominator")) is int and node["denominator"] != 0)
        elif operation == "current_original_source_parameter":
            closure_valid &= node.get("Z_independent") is True
        else:
            closure_valid = False
        if any(type(ref) is not int or not 0 <= ref < ident for ref in refs):
            closure_valid = False
        pending.extend(refs)
    checks["complete_known_operations_Z_free_closure"] = bool(closure_valid)
    return {"valid": all(checks.values()), "checks": checks, "closure": sorted(closure),
            "source_only": True, "numeric_evaluation": False,
            "adapter_Rrestore_Rm_alias_assumed": False}


def _source_witness(records: Mapping[str, object], gates: list[str]) -> dict[str, object]:
    texts = {name: _record(records, name, gates) for name in RECORD_HASHES}
    latest = _json(texts[LATEST_REPORT], LATEST_REPORT, gates)
    provider = _json(texts[PROVIDER_REPORT], PROVIDER_REPORT, gates)
    memory = _json(texts[MEMORY_REPORT], MEMORY_REPORT, gates)
    generic = _json(texts[GENERIC_REPORT], GENERIC_REPORT, gates)
    archive = _json(texts[C5_ARCHIVE], C5_ARCHIVE, gates)
    loop_archive = _json(texts[LOOP_ARCHIVE], LOOP_ARCHIVE, gates)
    adapter = _json(texts[ADAPTER_REPORT], ADAPTER_REPORT, gates)
    r100_archive = _json(texts[R100_C6_ARCHIVE], R100_C6_ARCHIVE, gates)
    true_c6_archive = _json(texts[TRUE_C6_ARCHIVE], TRUE_C6_ARCHIVE, gates)
    c5_tree = _ast(texts[C5_SOURCE], C5_SOURCE, gates)
    loop_tree = _ast(texts[LOOP_SOURCE], LOOP_SOURCE, gates)
    c5_names = _function_names(c5_tree)
    loop_names = _function_names(loop_tree)
    component = latest.get("component_result", {})
    latest_interface = latest.get("next_finite_N_interface", {})
    latest_scope = latest.get("scope", {})
    latest_result_flags = latest.get("result_flags", {})
    provider_component = provider.get("component_result", {})
    provider_scope = provider.get("scope", {})
    generic_component = generic.get("component_result", {})
    memory_component = memory.get("component_result", {})
    patch_provider = memory_component.get("patch_provider", {}) if isinstance(memory_component, Mapping) else {}
    patch_checks = patch_provider.get("checks", {}) if isinstance(patch_provider, Mapping) else {}

    latest_ok = (
        isinstance(component, Mapping)
        and component.get("passed") is True
        and component.get("nodes_digest") == EXPECTED_NODES_DIGEST
        and component.get("node_count") == EXPECTED_NODE_COUNT
        and component.get("native_C2_leading_patch_defined") is True
        and component.get("native_C2_C6_patch_source_identity_proved") is True
        and component.get("same_unique_leading_control_map_proved") is True
        and component.get("legacy_feedback_reference_alias_proved") is False
        and component.get("ordinary_Z_orders") == list(ORDINARY_Z_ORDERS)
        and latest_scope.get("first_ten_nonreset_history_identity_proved") is True
        and all(latest_scope.get(flag) is False for flag in (
            "same_density_certificate", "finite_N_history_identity", "physical_Q_y_Z1",
            "global_five_moment_closure", "full_goal_complete"))
        and latest_result_flags.get("source_only") is True
        and latest_result_flags.get("numeric_values_computed") is False
        and latest_result_flags.get("old_operators_replayed") is False
    )
    if not latest_ok:
        _gate(gates, "latest_native_C2_composition_prerequisite_unproved")
    interface_ok = (
        isinstance(latest_interface, Mapping)
        and latest_interface.get("schema") == "native_C2_patch_next_finite_N_interface_v1"
        and latest_interface.get("source_scan_only") is True
        and latest_interface.get("source_hashes_verified") is True
        and latest_interface.get("patch_archive_prefix") == PATCH_ARCHIVE_PREFIX
        and latest_interface.get("phase_handle") == PATCH_PHASE_HANDLE
        and latest_interface.get("coordinate_handle") == PATCH_COORDINATE_HANDLE
        and latest_interface.get("patch_loop_A") == list(PATCH_LOOP_A[:3])
        and latest_interface.get("patch_loop_B_over_Pstar") == list(PATCH_LOOP_B[:3])
        and _get(latest_interface, "patch_transport.coordinate") == "x=R/Rm"
        and _get(latest_interface, "patch_transport.measure") == "dx/x"
        and _get(latest_interface, "patch_transport.p_memory") == "exactly one"
        and _get(latest_interface, "phase.origin") == "hb*sc/2"
    )
    if not interface_ok:
        _gate(gates, "next_finite_N_interface_prerequisite_unproved")
    provider_ok = (
        isinstance(provider_component, Mapping)
        and provider_component.get("passed") is True
        and provider_component.get("source_only") is True
        and provider_component.get("nodes_digest") == EXPECTED_NODES_DIGEST
        and provider_component.get("node_count") == EXPECTED_NODE_COUNT
        and provider_component.get("first_ten_finite_N_field_source_composition_proved") is True
        and provider_component.get("first_ten_density_C2_NC5_identity_proved") is True
        and provider_component.get("generic_density_identity_proved") is True
        and provider_component.get("ordinary_Z_orders") == list(ORDINARY_Z_ORDERS)
        and provider_scope.get("numeric_values_computed") is False
        and provider_scope.get("old_operators_replayed") is False
    )
    if not provider_ok:
        _gate(gates, "accepted_finite_N_provider_prerequisite_unproved")
    source_equality = provider_component.get("actual_source_equality_theorem", {}) if isinstance(provider_component, Mapping) else {}
    source_equality_checks = source_equality.get("checks", {}) if isinstance(source_equality, Mapping) else {}
    source_domain = source_equality.get("domain_contract", {}) if isinstance(source_equality, Mapping) else {}
    source_generic = provider_component.get("strong_generic_recovery_theorem", {}) if isinstance(provider_component, Mapping) else {}
    source_generic_checks = source_generic.get("checks", {}) if isinstance(source_generic, Mapping) else {}
    source_smoothness_ok = (
        isinstance(source_equality, Mapping)
        and source_equality.get("passed") is True
        and isinstance(source_equality_checks, Mapping)
        and source_equality_checks.get("ordinary_derivative_conversion_and_smoothness") is True
        and source_equality.get("same_actual_analytic_source_domain") is True
        and isinstance(source_domain, Mapping)
        and source_domain.get("identity_is_conditional_on_common_original_domain") is True
        and isinstance(source_domain.get("fixed_parameters"), list)
        and isinstance(source_domain.get("active_branch"), list)
        and isinstance(source_domain.get("underlying_regularities"), str)
        and "genuine per-side C3" in source_domain.get("underlying_regularities", "")
        and isinstance(source_generic, Mapping)
        and source_generic.get("passed") is True
        and isinstance(source_generic_checks, Mapping)
        and source_generic_checks.get("C6_recovery") is True
        and "genuine C6 definition supplies C3 regularity" in source_generic.get("regularity_lift", "")
        and source_equality.get("numerical_admissibility_or_conditioning_proved") is False
        and source_equality.get("numerical_cone_or_admissibility_margins_claimed") is False
    )
    if not source_smoothness_ok:
        _gate(gates, "accepted_C6_smoothness_common_domain_premise_unproved")
    generic_check = generic_component.get("generic_density_identity_check", {}) if isinstance(generic_component, Mapping) else {}
    generic_ok = (
        isinstance(generic_component, Mapping)
        and generic_component.get("generic_density_identity_proved") is True
        and isinstance(generic_check, Mapping)
        and generic_check.get("proved") is True
        and generic_check.get("signed_densities") == {name: True for name in MOMENT_NAMES}
        and generic_check.get("ordinary_Z_orders") == list(ORDINARY_Z_ORDERS)
        and _get(generic_check, "public_api.checks.F_equals_N_r") is True
        and _get(generic_check, "public_api.checks.A_zero_entire_extension") is True
        and _get(generic_check, "public_api.checks.constant_inverse_N_rows") is True
    )
    if not generic_ok:
        _gate(gates, "accepted_generic_C2_C5_density_identity_unproved")
    memory_ok = (
        isinstance(memory_component, Mapping)
        and memory_component.get("passed") is True
        and memory_component.get("source_only") is True
        and memory_component.get("nodes_digest") == EXPECTED_NODES_DIGEST
        and memory_component.get("first_ten_nonreset_history_identity_proved") is True
        and memory_component.get("first_ten_history_endpoint_identity_proved") is True
        and memory_component.get("ordinary_Z_orders") == list(ORDINARY_Z_ORDERS)
        and memory_component.get("numeric_values_computed") is False
        and memory_component.get("old_operators_replayed") is False
    )
    if not memory_ok:
        _gate(gates, "accepted_nonreset_memory_prerequisite_unproved")
    patch_receipt_witness = bool(
        isinstance(patch_provider, Mapping)
        and all(patch_checks.get(key) is True for key in (
            "patch_binding_present", "patch_chart_present", "patch_nonzero_incoming_port_present",
            "patch_definition_retains_nonzero_Rm_inlet", "patch_definition_retains_log_transport",
            "postrestore_C2_outgoing_retained", "incoming_is_saved_C5_postrestore_endpoint"))
        and patch_provider.get("patch_identity_stays_false") is True
    )
    if not patch_receipt_witness:
        _gate(gates, "accepted_patch_inlet_binding_witness_unproved")
    source_formula_ok = (
        _has(texts[C5_SOURCE], "def transport_log", "def native_derivative", "def pressure_continuation",
             "finite=engine.finite_N_functions(g,roots,loop,self.N)",
             "self.phase_origin=self.loop.phase_origin", "self.incoming=self.loop.Rm_correction",
             "history,parts,memory=transport_log", "physical_measure='dt/t'",
             "predecessor_memory_never_reset=True")
        and _has(texts[LOOP_SOURCE], "def finite_N_functions", "def transport", "def loop_functions",
                 "phase_origin=g.mul(background.bridge.hb,self.sc,g.constant('1/2'))",
                 "phase=g.unary('fractional_part',g.mul(self.N,offset))",
                 "A=lazy(a.scale(a.mul(roots['a'],chi),'1/2'))",
                 "B=lazy(a.scale(a.mul(roots['E'],M),'1/2'))",
                 "exm1=FunctionJet", "five_signed_density_increments=inc",
                 "def __init__(self,background,require_checked=True)")
        and {"transport_log", "native_derivative", "pressure_continuation",
             "CurrentOriginalTrueRmPatchFiniteNC5Functions"}.issubset(c5_names)
        and {"finite_N_functions", "transport", "loop_functions", "CurrentOriginalPrescribedLoopC5Functions"}.issubset(loop_names)
    )
    if not source_formula_ok:
        _gate(gates, "finite_N_patch_source_formula_AST_unproved")
    archive_shape = _patch_archive_witness(archive, gates)
    retained = component.get("retained_construction_data", {}) if isinstance(component, Mapping) else {}
    retained_ok = all(retained.get(key) == value for key, value in RETAINED_CONSTRUCTION.items())
    if not retained_ok:
        _gate(gates, "retained_construction_data_changed")
    # Admit the correction inlet only from authenticated snapshots.  The
    # prescribed-loop archive must share its exact canonical graph prefix with
    # the actual finite-N patch, including the postrestore C5 rows, while the
    # accepted memory report must bind the C2 outgoing endpoint map.
    prescribed_graph = loop_archive.get("exact_graph_nodes", []) if isinstance(loop_archive, Mapping) else []
    patch_graph = archive.get("exact_graph_nodes", []) if isinstance(archive, Mapping) else []
    prefix_equal = (
        isinstance(prescribed_graph, list) and isinstance(patch_graph, list)
        and len(prescribed_graph) == LOOP_ARCHIVE_PREFIX
        and len(patch_graph) == PATCH_ARCHIVE_PREFIX
        and prescribed_graph == patch_graph[:LOOP_ARCHIVE_PREFIX]
        and _digest(prescribed_graph) == LOOP_ARCHIVE_PREFIX_SHA256
    )
    if not prefix_equal:
        _gate(gates, "prescribed_C5_patch_common_prefix_unproved")
    true_c6_graph = true_c6_archive.get("exact_graph_nodes", []) if isinstance(true_c6_archive, Mapping) else []
    true_c6_bindings = true_c6_archive.get("true_Rm_patch_C6_bindings", {}) if isinstance(true_c6_archive, Mapping) else {}
    true_c6_chart = true_c6_bindings.get("chart", {}) if isinstance(true_c6_bindings, Mapping) else {}
    c6_prefix_equal = (
        isinstance(true_c6_graph, list) and isinstance(patch_graph, list)
        and len(true_c6_graph) == TRUE_C6_PREFIX
        and len(patch_graph) == PATCH_ARCHIVE_PREFIX
        and true_c6_graph == patch_graph[:TRUE_C6_PREFIX]
        and _digest(true_c6_graph) == TRUE_C6_PREFIX_SHA256
    )
    c6_fields = true_c6_chart.get("fields", {}) if isinstance(true_c6_chart, Mapping) else {}
    c6_radial = true_c6_chart.get("radial", {}) if isinstance(true_c6_chart, Mapping) else {}
    patch_chart = archive.get("Rm_patch_finite_N_C5_bindings", {}).get("chart", {}) if isinstance(archive, Mapping) else {}
    patch_roots = patch_chart.get("roots", {}) if isinstance(patch_chart, Mapping) else {}
    patch_background_radial = archive.get("Rm_patch_finite_N_C5_bindings", {}).get("radial_background", {}) if isinstance(archive, Mapping) else {}
    c6_chart_ports_equal = (
        isinstance(c6_fields, Mapping) and isinstance(c6_radial, Mapping)
        and isinstance(patch_roots, Mapping) and isinstance(patch_background_radial, Mapping)
        and all(c6_fields.get(name, [])[:6] == patch_roots.get(name) for name in ("E", "V"))
        and c6_radial.get("E_y", [])[:6] == patch_background_radial.get("E", [])[:6]
        and c6_radial.get("V_y", [])[:6] == patch_background_radial.get("V", [])[:6]
        and true_c6_chart.get("radius") == patch_chart.get("radius")
        and true_c6_chart.get("physical_measure") == patch_chart.get("physical_measure")
    )
    c6_parameter_prefix_equal = (
        isinstance(true_c6_graph, list) and isinstance(patch_graph, list)
        and all(true_c6_graph[index] == patch_graph[index]
                for index in (10, 18, 21, 22, 24954))
        and _get(archive, "Rm_patch_finite_N_C5_bindings.N") == 18
        and _get(archive, "Rm_patch_finite_N_C5_bindings.eta") == 21
        and _get(archive, "Rm_patch_finite_N_C5_bindings.dstar") == 22
        and _get(archive, "Rm_patch_finite_N_C5_bindings.sc") == 10
        and _get(archive, "Rm_patch_finite_N_C5_bindings.phase_origin") == 24954
    )
    c6_definition = true_c6_bindings.get("definition", {}) if isinstance(true_c6_bindings, Mapping) else {}
    c6_provider_identity = (
        c6_prefix_equal and c6_chart_ports_equal and c6_parameter_prefix_equal
        and _get(component, "qualified_source_ports.C6_archive_sha256") == TRUE_C6_ARCHIVE_SHA256
        and _get(component, "qualified_source_ports.C6_prefix") == TRUE_C6_PREFIX
        and _get(component, "qualified_source_ports.C6_system") == 193470
        and source_smoothness_ok
        and isinstance(c6_definition, Mapping)
        and c6_definition.get("actual_E_V_m_h_k_e_p_P0_pressure_and_Ey_Vy_C6") is True
        and c6_definition.get("full_signed_generic_C5") is True
        and c6_definition.get("leading_control_not_finite_N_correction_and_no_history_reset") is True
    )
    if not c6_provider_identity:
        _gate(gates, "native_C6_to_C5_patch_provider_identity_unproved")
    native_geometry = _native_geometry_witness(patch_graph, patch_chart, true_c6_bindings.get("inlet", {}))
    native_geometry_source_ok = bool(c6_prefix_equal and native_geometry.get("valid"))
    if not native_geometry_source_ok:
        _gate(gates, "native_C5_geometry_Z_independence_source_unproved")
    prescribed_bindings = loop_archive.get("prescribed_loop_C5_bindings", {}) if isinstance(loop_archive, Mapping) else {}
    prescribed_incoming = prescribed_bindings.get("Rm_correction", {}) if isinstance(prescribed_bindings, Mapping) else {}
    patch_incoming = _get(archive, "Rm_patch_finite_N_C5_bindings.incoming", {})
    bridge_rows_equal = (
        isinstance(prescribed_incoming, Mapping) and isinstance(patch_incoming, Mapping)
        and prescribed_incoming == patch_incoming
        and set(prescribed_incoming) == set(MOMENT_NAMES)
        and all(prescribed_incoming.get(name) == [
            121756 + 6 * index + offset for offset in range(6)
        ] for index, name in enumerate(MOMENT_NAMES))
    )
    if not bridge_rows_equal:
        _gate(gates, "prescribed_C5_patch_incoming_snapshot_mismatch")
    adapter_component = adapter.get("component_result", {}) if isinstance(adapter, Mapping) else {}
    adapter_ports = _get(adapter, "component_result.qualified_endpoint_ports.parameter", {})
    r100_graph = r100_archive.get("exact_graph_nodes", []) if isinstance(r100_archive, Mapping) else []
    pstar_binding_ok = (
        isinstance(adapter_component, Mapping)
        and adapter_component.get("passed") is True
        and adapter_component.get("source_only") is True
        and adapter_component.get("nodes_digest") == EXPECTED_NODES_DIGEST
        and adapter_component.get("legacy_feedback_reference_alias_proved") is False
        and isinstance(adapter_ports, Mapping)
        and adapter_ports.get("Pstar") == PSTAR_SOURCE_HANDLE
        and isinstance(r100_graph, list)
        and len(r100_graph) == R100_C6_GRAPH_COUNT
        and isinstance(r100_graph[PSTAR_SOURCE_HANDLE], Mapping)
        and r100_graph[PSTAR_SOURCE_HANDLE].get("operation") == "analytic_unary"
        and r100_graph[PSTAR_SOURCE_HANDLE].get("name") == "exp"
    )
    if not pstar_binding_ok:
        _gate(gates, "accepted_Pstar_source_parameter_port_unproved")
    memory_outgoing = _get(memory_component, "patch_provider.postrestore_outgoing_Z012", {})
    memory_outgoing_witness = (
        isinstance(memory_component, Mapping)
        and memory_component.get("nodes_digest") == EXPECTED_NODES_DIGEST
        and isinstance(memory_outgoing, Mapping)
        and memory_outgoing == {name: list(row) for name, row in PATCH_C2_INCOMING.items()}
        and _get(memory_component, "patch_provider.postrestore_outgoing_source")
        == "C2_density_transport.windows[postrestore].outgoing"
        and _get(memory_component, "patch_provider.patch_identity_stays_false") is True
    )
    if not memory_outgoing_witness:
        _gate(gates, "accepted_memory_C2_outgoing_endpoint_unproved")
    bridge_identity = prefix_equal and bridge_rows_equal and patch_receipt_witness and memory_outgoing_witness
    incoming_identity = bridge_identity
    if not incoming_identity:
        _gate(gates, "patch_correction_inlet_identity_unproved")
    return {
        "texts": texts,
        "latest": latest,
        "provider": provider,
        "memory": memory,
        "generic": generic,
        "archive": archive,
        "loop_archive": loop_archive,
        "adapter": adapter,
        "r100_archive": r100_archive,
        "archive_shape": archive_shape,
        "latest_ok": latest_ok,
        "interface_ok": interface_ok,
        "provider_ok": provider_ok,
        "source_smoothness_ok": source_smoothness_ok,
        "generic_ok": generic_ok,
        "memory_ok": memory_ok,
        "patch_receipt_witness": patch_receipt_witness,
        "source_formula_ok": source_formula_ok,
        "incoming_identity": incoming_identity,
        "bridge_prefix_witness": prefix_equal,
        "bridge_rows_equal": bridge_rows_equal,
        "memory_outgoing_witness": memory_outgoing_witness,
        "pstar_binding_ok": pstar_binding_ok,
        "c6_provider_identity": c6_provider_identity,
        "c6_prefix_equal": c6_prefix_equal,
        "c6_chart_ports_equal": c6_chart_ports_equal,
        "c6_parameter_prefix_equal": c6_parameter_prefix_equal,
        "native_geometry_source_ok": native_geometry_source_ok,
        "native_geometry_source_witness": native_geometry,
        "retained_ok": retained_ok,
        "c5_names": sorted(c5_names),
        "loop_names": sorted(loop_names),
    }


def _qualified_archive_jet(
    dag: _DAG, rows: object, *, label: str, path: str, role: str,
    witness: Mapping[str, object], gates: list[str], metadata: Mapping[str, object] | None = None,
    coordinate_binding: ScalarID | None = None, z_binding: ScalarID | None = None,
    source_width: int = 6,
) -> JetList:
    values = list(rows)[:3] if isinstance(rows, list) else ["missing"] * 3
    valid_count = PATCH_ARCHIVE_PREFIX
    valid = _valid_row(rows, valid_count, source_width)
    if not valid:
        _gate(gates, "invalid_patch_archive_jet:" + path)
    while len(values) < 3:
        values.append("missing")
    extra = dict(metadata or {})
    return [_source_port(
        dag, report_sha256=C5_ARCHIVE_SHA256,
        namespace="original_true_Rm_patch_finite_N_C5_functions_v1",
        source_node=value, archive=C5_ARCHIVE, label=label + "[" + str(index) + "]",
        role=role, source_path=path + "[" + str(index) + "]",
        valid=valid and _valid_handle(value, valid_count),
        coordinate_handle=PATCH_COORDINATE_HANDLE if coordinate_binding is not None else None,
        phase_handle=PATCH_PHASE_HANDLE if coordinate_binding is not None else None,
        bound_variable="true_Rm_patch_x" if coordinate_binding is not None else None,
        coordinate_binding=coordinate_binding, z_binding=z_binding,
        source_coordinate_handle=PATCH_COORDINATE_HANDLE if coordinate_binding is not None else None,
        source_z_handle=PATCH_Z_HANDLE if z_binding is not None else None,
        source_coordinate_name="true_Rm_patch_x" if coordinate_binding is not None else None,
        source_z_name="Z" if z_binding is not None else None,
        source_row_width=source_width, ordinary_rows_imported=[0, 1, 2], **extra)
        for index, value in enumerate(values)]


def _latest_inlet_jet(
    dag: _DAG, name: str, gates: list[str], memory: Mapping[str, object],
    *, z_binding: ScalarID | None = None,
) -> JetList:
    # These are accepted first-ten C2 correction ports from the memory report's
    # postrestore endpoint map.  The next-finite-N scan is descriptive only;
    # it is never used as the source owner for these integers.
    component = memory.get("component_result", {}) if isinstance(memory, Mapping) else {}
    patch_provider = component.get("patch_provider", {}) if isinstance(component, Mapping) else {}
    rows = patch_provider.get("postrestore_outgoing_Z012", {}).get(name) if isinstance(patch_provider, Mapping) else None
    valid = rows == list(PATCH_C2_INCOMING[name])
    if not valid:
        _gate(gates, "memory_C2_incoming_handle_mismatch:" + name)
    source_valid = (
        valid
        and component.get("nodes_digest") == EXPECTED_NODES_DIGEST
        and patch_provider.get("postrestore_outgoing_source") == "C2_density_transport.windows[postrestore].outgoing"
        and patch_provider.get("patch_identity_stays_false") is True
    ) if isinstance(component, Mapping) and isinstance(patch_provider, Mapping) else False
    if not source_valid:
        _gate(gates, "memory_C2_incoming_source_identity_unproved:" + name)
    values = list(rows) if isinstance(rows, list) else list(PATCH_C2_INCOMING[name])
    return [_source_port(
        dag, report_sha256=MEMORY_REPORT_SHA256,
        namespace="original_finite_N_nonreset_memory_source_binding_v1",
        source_node=value, archive=MEMORY_REPORT,
        label="accepted_C2_patch_outgoing." + name + "[" + str(index) + "]",
        role="accepted_first_ten_C2_outgoing_correction_port",
        source_path="component_result.patch_provider.postrestore_outgoing_Z012." + name + "[" + str(index) + "]",
        valid=source_valid and type(value) is int and 0 <= value < EXPECTED_NODE_COUNT,
        incoming_is_correction_port=True, leading_moment_port=False,
        z_binding=z_binding, source_z_handle=PATCH_Z_HANDLE, source_z_name="Z",
        source_graph_digest=EXPECTED_NODES_DIGEST,
        source_graph_node_valid=source_valid,
        source_owner="memory_component.patch_provider.postrestore_outgoing_Z012",
    ) for index, value in enumerate(values)]


def _build_dag(nodes: Sequence[object], witness: Mapping[str, object], gates: list[str]) -> tuple[_DAG, dict[str, object]]:
    dag = _DAG()
    latest = witness["latest"]
    archive = witness["archive"]
    shape = witness["archive_shape"]
    graph = archive.get("exact_graph_nodes", []) if isinstance(archive, Mapping) else []
    # Native symbols and exact source parameters.
    Z = _symbol(dag, "Z", "[-1,1]")
    x = _symbol(dag, "x", "[1,e]")
    one = _one(dag)
    zero = _zero(dag)
    N = _source_port(dag, report_sha256=C5_ARCHIVE_SHA256,
                     namespace="original_true_Rm_patch_finite_N_C5_functions_v1",
                     source_node=_get(archive, "Rm_patch_finite_N_C5_bindings.N"),
                     archive=C5_ARCHIVE, label="patch.parameter.N",
                     role="source_parameter", source_path="Rm_patch_finite_N_C5_bindings.N",
                     valid=_valid_handle(_get(archive, "Rm_patch_finite_N_C5_bindings.N"), PATCH_ARCHIVE_PREFIX),
                     bound_variable="N_Z_independent", source_z_independent=True)
    Pstar = _source_port(
        dag, report_sha256=ADAPTER_REPORT_SHA256,
        namespace="original_native_C2_reference_adapter_v1",
        source_node=PSTAR_SOURCE_HANDLE, archive=R100_C6_ARCHIVE,
        label="adapter.parameter.Pstar", role="accepted_source_parameter",
        source_path="component_result.qualified_endpoint_ports.parameter.Pstar",
        valid=bool(witness.get("pstar_binding_ok")),
        bound_variable="Pstar_Z_independent", source_graph_node_valid=bool(witness.get("pstar_binding_ok")),
        source_graph_node=PSTAR_SOURCE_HANDLE, source_graph_count=R100_C6_GRAPH_COUNT,
        source_binding_report=ADAPTER_REPORT, source_binding_report_sha256=ADAPTER_REPORT_SHA256,
        source_archive_owner=R100_C6_ARCHIVE, source_archive_sha256=R100_C6_ARCHIVE_SHA256)
    invN = _jet_inv(dag, _jet_scalar(dag, N))
    njet = _jet_scalar(dag, N)
    xjet = _jet_scalar(dag, x)

    # Previous accepted graph fields and the native C5 patch provider ports.
    E = _saved_latest_jet(dag, latest, "E", gates=gates, path="previous.function_definitions.normalized_fields.E",
                           coordinate_binding=x, z_binding=Z)
    V = _saved_latest_jet(dag, latest, "V", gates=gates, path="previous.function_definitions.normalized_fields.V",
                           coordinate_binding=x, z_binding=Z)
    previous_P0 = _saved_latest_jet(dag, latest, "P0", gates=gates,
                                   path="previous.function_definitions.normalized_fields.P0",
                                   coordinate_binding=x, z_binding=Z)
    P0 = _qualified_archive_jet(
        dag, _get(archive, "Rm_patch_finite_N_C5_bindings.chart.modified_fields.P0"),
        label="native_C5.modified_fields.P0",
        path="Rm_patch_finite_N_C5_bindings.chart.modified_fields.P0",
        role="native_C5_independent_P0_source", witness=witness, gates=gates,
        coordinate_binding=x, z_binding=Z)
    radial_E = _get(latest, "component_result.function_definitions.radial_fields.E_y")
    radial_V = _get(latest, "component_result.function_definitions.radial_fields.V_y")
    radial_E_y = _saved_latest_jet(
        dag, latest, "E_y", section="radial_fields", gates=gates,
        path="previous.function_definitions.radial_fields.E_y", coordinate_binding=x, z_binding=Z)
    radial_V_y = _saved_latest_jet(
        dag, latest, "V_y", section="radial_fields", gates=gates,
        path="previous.function_definitions.radial_fields.V_y", coordinate_binding=x, z_binding=Z)
    # Keep the source IDs of the radial rows explicit even though the helper
    # above reads only normalized fields; these are qualified imported refs.
    if not isinstance(radial_E, list) or not isinstance(radial_V, list):
        _gate(gates, "previous_radial_source_fields_missing")
    chart = _get(archive, "Rm_patch_finite_N_C5_bindings.chart", {})
    roots = chart.get("roots", {}) if isinstance(chart, Mapping) else {}
    loop = chart.get("loop", {}) if isinstance(chart, Mapping) else {}
    finite = chart.get("finite_N", {}) if isinstance(chart, Mapping) else {}
    native_E = _qualified_archive_jet(dag, roots.get("E"), label="native_C5.roots.E",
                                      path="Rm_patch_finite_N_C5_bindings.chart.roots.E", role="native_C5_graph_ref",
                                      witness=witness, gates=gates, coordinate_binding=x, z_binding=Z)
    native_V = _qualified_archive_jet(dag, roots.get("V"), label="native_C5.roots.V",
                                      path="Rm_patch_finite_N_C5_bindings.chart.roots.V", role="native_C5_graph_ref",
                                      witness=witness, gates=gates, coordinate_binding=x, z_binding=Z)
    previous_E, previous_V = E, V
    leading_source_identity = {
        "E": dag.op("source_identity", tuple(previous_E + native_E), exact=True,
                    relation="accepted native C2 leading E equals same C6/C5 root E",
                    source_premise_report=LATEST_REPORT_SHA256,
                    source_premise="native_C2_C6_patch_source_identity_proved",
                    common_prefix_and_ports_checked=bool(witness.get("c6_provider_identity")),
                    conditional_common_original_domain=True, source_only=True),
        "V": dag.op("source_identity", tuple(previous_V + native_V), exact=True,
                    relation="accepted native C2 leading V equals same C6/C5 root V",
                    source_premise_report=LATEST_REPORT_SHA256,
                    source_premise="native_C2_C6_patch_source_identity_proved",
                    common_prefix_and_ports_checked=bool(witness.get("c6_provider_identity")),
                    conditional_common_original_domain=True, source_only=True),
    }
    # Evaluate no source: use the native provider roots themselves in every new
    # finite-N formula, retaining the prior-owner equality as a separate theorem.
    E, V = native_E, native_V
    native_A = _qualified_archive_jet(dag, loop.get("A"), label="native_C5.loop.A",
                                      path="Rm_patch_finite_N_C5_bindings.chart.loop.A", role="native_C5_graph_ref",
                                      witness=witness, gates=gates, coordinate_binding=x, z_binding=Z)
    native_B = _qualified_archive_jet(dag, loop.get("B_over_Pstar"), label="native_C5.loop.B_over_Pstar",
                                      path="Rm_patch_finite_N_C5_bindings.chart.loop.B_over_Pstar", role="native_C5_graph_ref",
                                      witness=witness, gates=gates, coordinate_binding=x, z_binding=Z,
                                      metadata={"normalized_units": "B_over_Pstar"})

    A_over_N = _jet_mul(dag, native_A, invN)
    exp_A_over_N = _jet_exp(dag, A_over_N)
    expm1_A_over_N = _jet_expm1(dag, A_over_N)
    r = _jet_mul(dag, E, expm1_A_over_N)
    s = _jet_mul(dag, native_B, invN)
    E_N = _jet_mul(dag, E, exp_A_over_N)
    V_N = _jet_add(dag, V, s)
    c5 = {
        "m": s,
        "h": r,
        "k": _jet_add(dag, _jet_mul(dag, V, r), _jet_mul(dag, E, s), _jet_mul(dag, r, s)),
        "e": _jet_sub(dag, _jet_add(dag, _jet_scale(dag, _jet_mul(dag, V, s), _const(dag, 2)),
                                  _jet_mul(dag, s, s)),
                       _jet_add(dag, _jet_mul(dag, E, r),
                                _jet_scale(dag, _jet_mul(dag, r, r), _const(dag, Fraction(1, 2))))),
        "p": _jet_add(dag, _jet_mul(dag, E, r),
                      _jet_scale(dag, _jet_mul(dag, r, r), _const(dag, Fraction(1, 2)))),
    }
    G = _jet_mul(dag, njet, expm1_A_over_N)
    F = _jet_mul(dag, E, G)
    c2_first = {
        "m": native_B,
        "h": F,
        "k": _jet_add(dag, _jet_mul(dag, V, F), _jet_mul(dag, E, native_B)),
        "e": _jet_sub(dag, _jet_scale(dag, _jet_mul(dag, V, native_B), _const(dag, 2)),
                       _jet_mul(dag, E, F)),
        "p": _jet_mul(dag, E, F),
    }
    c2_second = {
        "m": _jet_const(dag, 0),
        "h": _jet_const(dag, 0),
        "k": _jet_mul(dag, F, native_B),
        "e": _jet_sub(dag, _jet_mul(dag, native_B, native_B),
                       _jet_scale(dag, _jet_mul(dag, F, F), _const(dag, Fraction(1, 2)))),
        "p": _jet_scale(dag, _jet_mul(dag, F, F), _const(dag, Fraction(1, 2))),
    }
    c2 = {name: _jet_add(dag, c2_first[name], _jet_scale(dag, c2_second[name], invN))
          for name in MOMENT_NAMES}
    density_relations = {
        name: dag.op("identity_claim", tuple(c2[name] + c5[name]), exact=True,
                     relation="C2=N*C5", coefficient=name,
                     proof_scope="formal ordinary-Z coefficientwise identity; source bindings required",
                     source_only=True)
        for name in MOMENT_NAMES
    }

    # Provider rows are retained as qualified native ports and compared by a
    # source identity node.  They are never treated as numeric field values.
    native_finite = finite if isinstance(finite, Mapping) else {}
    provider_fields = {
        "E_N": _qualified_archive_jet(dag, native_finite.get("E_N"), label="native_C5.finite_N.E_N",
                                       path="Rm_patch_finite_N_C5_bindings.chart.finite_N.E_N", role="native_C5_graph_ref",
                                       witness=witness, gates=gates, coordinate_binding=x, z_binding=Z),
        "V_N": _qualified_archive_jet(dag, native_finite.get("V_N"), label="native_C5.finite_N.V_N",
                                       path="Rm_patch_finite_N_C5_bindings.chart.finite_N.V_N", role="native_C5_graph_ref",
                                       witness=witness, gates=gates, coordinate_binding=x, z_binding=Z),
    }
    provider_density = {
        name: _qualified_archive_jet(
            dag, _get(native_finite, "five_signed_density_increments." + name),
            label="native_C5.finite_N.density." + name,
            path="Rm_patch_finite_N_C5_bindings.chart.finite_N.five_signed_density_increments." + name,
            role="native_C5_graph_ref", witness=witness, gates=gates, coordinate_binding=x, z_binding=Z)
        for name in MOMENT_NAMES
    }
    provider_identity = {
        "E_N": dag.op("source_identity", tuple(E_N + provider_fields["E_N"]), exact=True,
                       relation="same native C5 finite-N E_N function", source_only=True),
        "V_N": dag.op("source_identity", tuple(V_N + provider_fields["V_N"]), exact=True,
                       relation="same native C5 finite-N V_N function", source_only=True),
    }
    provider_density_identity = {
        name: dag.op("source_identity", tuple(c5[name] + provider_density[name]), exact=True,
                     relation="same native C5 signed density function", coefficient=name,
                     source_only=True)
        for name in MOMENT_NAMES
    }

    # Section 11 correction histories.  The incoming C2 ports are the accepted
    # first-ten outgoing correction rows, never the leading moment rows.
    s_variable = _symbol(dag, "s", "[1,x]")
    histories: dict[str, JetList] = {}
    c5_histories: dict[str, JetList] = {}
    incoming: dict[str, JetList] = {}
    c5_incoming: dict[str, JetList] = {}
    kernels: dict[str, ScalarID] = {}
    endpoint_x1: dict[str, JetList] = {}
    endpoint_xe: dict[str, JetList] = {}
    derivative_identities: dict[str, ScalarID] = {}
    incoming_source = _get(archive, "Rm_patch_finite_N_C5_bindings.incoming", {})
    for name, rate in RATES.items():
        incoming[name] = _latest_inlet_jet(dag, name, gates, witness["memory"], z_binding=Z)
        c5_incoming[name] = _qualified_archive_jet(
            dag, incoming_source.get(name) if isinstance(incoming_source, Mapping) else None,
            label="native_C5.incoming." + name,
            path="Rm_patch_finite_N_C5_bindings.incoming." + name,
            role="native_C5_incoming_correction_port", witness=witness, gates=gates,
            z_binding=Z, metadata={"predecessor_memory_retained": True})
        rate_fraction = _fraction(rate)
        if rate_fraction == 0:
            kernel = one
            decay = one
        else:
            log_s = _log(dag, s_variable)
            log_x = _log(dag, x)
            log_kernel = _mul(dag, _const(dag, rate_fraction), _sub(dag, log_s, log_x))
            kernel = _exp(dag, log_kernel)
            decay = _exp(dag, _mul(dag, _const(dag, -rate_fraction), log_x))
        kernels[name] = kernel
        rho_at_s = _substitute_jet(
            dag, c2[name], x, s_variable,
            source_coordinate_substitution=True, source_coordinate="true_Rm_patch_x",
            substituted_coordinate="s", full_jet_substitution=True)
        weighted_kernel = _div(dag, kernel, s_variable)
        integrand = _jet_mul(dag, _jet_scalar(dag, weighted_kernel), rho_at_s)
        contribution = _integral_jet(
            dag, integrand, s_variable, one, x,
            radial_measure="ds/s", log_kernel=True,
            inverse_radial_factor_in_integrand=True,
            kernel_with_measure=weighted_kernel,
            kernel_formula="exp(rate*(log(s)-log(x)))",
            original_rate=rate, original_chart="Rm_patch",
            physical_jacobian_applied_once=True,
            predecessor_memory_never_reset=True,
            full_jet_product_before_component_integral=True)
        histories[name] = _jet_add(dag, _jet_mul(dag, _jet_scalar(dag, decay), incoming[name]), contribution)
        c5_rho_at_s = _substitute_jet(
            dag, c5[name], x, s_variable,
            source_coordinate_substitution=True, source_coordinate="true_Rm_patch_x",
            substituted_coordinate="s", full_jet_substitution=True)
        c5_integrand = _jet_mul(dag, _jet_scalar(dag, weighted_kernel), c5_rho_at_s)
        c5_contribution = _integral_jet(
            dag, c5_integrand, s_variable, one, x,
            radial_measure="ds/s", log_kernel=True,
            inverse_radial_factor_in_integrand=True,
            kernel_with_measure=weighted_kernel,
            kernel_formula="exp(rate*(log(s)-log(x)))",
            original_rate=rate, original_chart="Rm_patch", physical_jacobian_applied_once=True,
            predecessor_memory_never_reset=True, full_jet_product_before_component_integral=True)
        c5_histories[name] = _jet_add(dag, _jet_mul(dag, _jet_scalar(dag, decay), c5_incoming[name]), c5_contribution)
        endpoint_x1[name] = _substitute_jet(dag, histories[name], x, one,
                                            endpoint="x=1", function_substitution=True)
        x_exp_one = _exp(dag, one)
        endpoint_xe[name] = _substitute_jet(dag, histories[name], x, x_exp_one,
                                            endpoint="x=e", function_substitution=True)
        derivative_identities[name] = dag.op(
            "source_identity", (x, histories[name][0], c2[name][0]), exact=True,
            identity="x*H_x=rho-rate*H", rate=rate,
            coordinate="x", source_density="C2rho", source_only=True)

    history_relations = {
        name: dag.op("identity_claim", tuple(histories[name] + c5_histories[name]), exact=True,
                     relation="C2_history=N*C5_history", coefficient=name,
                     incoming_source="qualified C2 outgoing and native C5 postrestore inlet",
                     source_only=True)
        for name in MOMENT_NAMES
    }

    # Geometry and radial source fields are explicit.  Build radius = Rm*x and
    # physical measure = 1/x once; no physical Jacobian is hidden in kernels.
    previous_geometry = _get(latest, "component_result.function_definitions.geometry", {})
    previous_rm = previous_geometry.get("Rm_radius") if isinstance(previous_geometry, Mapping) else None
    if not _valid_jet(previous_rm, _latest_graph_ids(latest)):
        _gate(gates, "invalid_previous_geometry_Rm_radius")
    previous_rm_radius = [
        _source_port(dag, report_sha256=LATEST_REPORT_SHA256,
                     namespace="original_native_C2_patch_composition_v1",
                     source_node=value, archive=LATEST_REPORT,
                     label="previous.geometry.Rm_radius[" + str(index) + "]",
                     role="qualified_previous_DAG_reference",
                     source_path="component_result.function_definitions.geometry.Rm_radius[" + str(index) + "]",
                     valid=value in _latest_graph_ids(latest), source_graph_digest=EXPECTED_NODES_DIGEST,
                     source_graph_node_valid=value in _latest_graph_ids(latest),
                     coordinate_binding=x, z_binding=Z,
                     source_coordinate_handle=LATEST_X_NODE, source_z_handle=LATEST_Z_NODE,
                     source_coordinate_name="x", source_z_name="Z")
        for index, value in enumerate(list(previous_rm)[:3] if isinstance(previous_rm, list) else ["missing"] * 3)
    ]
    native_rm_scalar = _source_port(
        dag, report_sha256=C5_ARCHIVE_SHA256,
        namespace="original_true_Rm_patch_finite_N_C5_functions_v1",
        source_node=85743, archive=C5_ARCHIVE,
        label="native_C5.geometry.Rm_radius", role="native_C5_Rm_radius_source",
        source_path="exact_graph_nodes[85743]",
        valid=bool(witness.get("native_geometry_source_ok")),
        native_C6_inlet_source_path="true_Rm_patch_C6_bindings.inlet.radius",
        ordinary_Z_independent=bool(witness.get("native_geometry_source_ok")),
        source_ancestor_closure_Z_free=bool(witness.get("native_geometry_source_ok")))
    rm_radius = _jet_scalar(dag, native_rm_scalar)
    radius = _jet_mul(dag, rm_radius, xjet)
    physical_measure = _jet_inv(dag, xjet)
    native_radius = _source_port(
        dag, report_sha256=C5_ARCHIVE_SHA256,
        namespace="original_true_Rm_patch_finite_N_C5_functions_v1",
        source_node=chart.get("radius"), archive=C5_ARCHIVE,
        label="native_C5.chart.radius", role="native_C5_graph_ref",
        source_path="Rm_patch_finite_N_C5_bindings.chart.radius",
        valid=_valid_handle(chart.get("radius"), PATCH_ARCHIVE_PREFIX),
        coordinate_binding=x, z_binding=Z,
        source_coordinate_handle=PATCH_COORDINATE_HANDLE, source_z_handle=PATCH_Z_HANDLE,
        source_coordinate_name="true_Rm_patch_x", source_z_name="Z")
    geometry_identity = dag.op(
        "source_identity", tuple(radius + _jet_scalar(dag, native_radius)), exact=True,
        relation="native chart radius=Rm*x; all ordinary Z derivatives zero",
        native_graph_product_and_Z_free_ancestor_closure_checked=bool(witness.get("native_geometry_source_ok")),
        prior_adapter_Rrestore_to_Rm_alias_assumed=False,
        conditional_common_original_domain=True, source_only=True)
    leading_radial_native = {
        name + "_y": _qualified_archive_jet(
            dag, _get(archive, "Rm_patch_finite_N_C5_bindings.radial_background." + name),
            label="native_C5.leading_radial." + name,
            path="Rm_patch_finite_N_C5_bindings.radial_background." + name,
            role="native_C5_leading_radial_ref", witness=witness, gates=gates,
            coordinate_binding=x, z_binding=Z, source_width=7)
        for name in ("E", "V")
    }
    native_radial = {
        name: _qualified_archive_jet(
            dag, _get(archive, "Rm_patch_finite_N_C5_bindings.radial." + name),
            label="native_C5.finite_N_radial." + name,
            path="Rm_patch_finite_N_C5_bindings.radial." + name,
            role="native_C5_finite_N_radial_ref", witness=witness, gates=gates,
            coordinate_binding=x, z_binding=Z)
        for name in ("E_y", "V_y")
    }
    radial_identities = {
        name: dag.op("source_identity", tuple(previous + leading_radial_native[name]),
                     exact=True, relation="same native leading radial derivative source",
                     source_premise_report=LATEST_REPORT_SHA256,
                     source_premise="native_C2_C6_patch_source_identity_proved",
                     common_prefix_and_ports_checked=bool(witness.get("c6_provider_identity")),
                     conditional_common_original_domain=True, source_only=True)
        for name, previous in (("E_y", radial_E_y), ("V_y", radial_V_y))
    }
    leading_moments = {
        name: _saved_latest_jet(
            dag, latest, name, gates=gates,
            path="previous.function_definitions.normalized_fields." + name,
            coordinate_binding=x, z_binding=Z)
        for name in MOMENT_NAMES
    }
    finite_N_moments = {name: _jet_add(dag, leading_moments[name], c5_histories[name])
                        for name in MOMENT_NAMES}
    finite_N_pressure = _jet_add(dag, P0, finite_N_moments["p"])
    canonical_amplitude = _exp(dag, _const(dag, Fraction(-3, 5)))
    source_derivative_rule = dag.op(
        "source_identity", tuple(radial_E_y + radial_V_y), exact=True,
        identity="radial physical rows are source functions; no numeric evaluation",
        physical_jacobian_once=True, source_only=True,
        higher_derivative_source_premise="accepted C6 source equality and genuine C6 regularity on the common original domain",
        ordinary_Z012_is_not_higher_order_certificate=True)
    outputs = {
        "normalized_fields": {
            "E": E_N, "V": V_N, **finite_N_moments,
            "P0": P0, "pressure": finite_N_pressure,
        },
        "physical_fields": {"U": _jet_scale(dag, V_N, Pstar),
                            "pressure": finite_N_pressure,
                            "P0": P0, "Utheta": _jet_scale(dag, E_N, Pstar),
                            "Uz": _jet_scale(dag, V_N, Pstar)},
        "radial_fields": native_radial,
        "leading_fields": {"E": E, "V": V, **leading_moments},
        "leading_radial_fields": {"E_y": radial_E_y, "V_y": radial_V_y},
        "retained_unaliased_adapter_data": {"P0": previous_P0,
                                              "Rm_radius_label_from_prior_owner": previous_rm_radius,
                                              "geometry_alias_proved": False,
                                              "P0_alias_proved": False},
        "geometry": {"Rm_radius": rm_radius, "radius": radius,
                      "physical_measure": physical_measure, "native_x": xjet},
        "finite_N": {"A_over_N": A_over_N, "r": r, "s": s,
                     "E_N": E_N, "V_N": V_N, "C5_density": c5,
                     "C2_G": G, "C2_F": F, "C2_first": c2_first,
                     "C2_second": c2_second, "C2_density": c2},
        "provider_fields": provider_fields,
        "provider_density": provider_density,
        "source_identity_nodes": {"fields": provider_identity,
                                   "leading_fields": leading_source_identity,
                                   "leading_radial_fields": radial_identities,
                                   "geometry": geometry_identity,
                                   "density": provider_density_identity,
                                   "C2_C5_density": density_relations,
                                   "histories": history_relations,
                                   "radial": source_derivative_rule},
        "histories": {"C2": histories, "C5": c5_histories,
                      "incoming_C2": incoming, "incoming_C5": c5_incoming,
                      "log_kernels": kernels, "endpoint_x1": endpoint_x1,
                      "endpoint_xe": endpoint_xe,
                      "derivative_identity": derivative_identities},
        "source_parameters": {"N": N, "Pstar": Pstar,
                               "phase": _source_port(
                                   dag, report_sha256=C5_ARCHIVE_SHA256,
                                   namespace="original_true_Rm_patch_finite_N_C5_functions_v1",
                                   source_node=_get(archive, "Rm_patch_finite_N_C5_bindings.chart.actual_phase"),
                                   archive=C5_ARCHIVE, label="native_C5.phase",
                                   role="source_parameter", source_path="Rm_patch_finite_N_C5_bindings.chart.actual_phase",
                                   valid=_valid_handle(_get(archive, "Rm_patch_finite_N_C5_bindings.chart.actual_phase"), PATCH_ARCHIVE_PREFIX),
                                   coordinate_handle=PATCH_COORDINATE_HANDLE, phase_handle=PATCH_PHASE_HANDLE,
                                   bound_variable="native_phase", coordinate_binding=x, z_binding=Z,
                                   source_coordinate_handle=PATCH_COORDINATE_HANDLE, source_z_handle=PATCH_Z_HANDLE,
                                   source_coordinate_name="true_Rm_patch_x", source_z_name="Z"),
                               "phase_origin": _source_port(
                                   dag, report_sha256=C5_ARCHIVE_SHA256,
                                   namespace="original_true_Rm_patch_finite_N_C5_functions_v1",
                                   source_node=_get(archive, "Rm_patch_finite_N_C5_bindings.phase_origin"),
                                   archive=C5_ARCHIVE, label="native_C5.phase_origin",
                                   role="source_parameter", source_path="Rm_patch_finite_N_C5_bindings.phase_origin",
                                   valid=_valid_handle(_get(archive, "Rm_patch_finite_N_C5_bindings.phase_origin"), PATCH_ARCHIVE_PREFIX),
                                   coordinate_handle=PATCH_COORDINATE_HANDLE, phase_handle=PATCH_PHASE_HANDLE,
                                   bound_variable="hb*sc/2")},
        "canonical_amplitude": canonical_amplitude,
        "coordinate_bindings": {
            "native_x": "true_Rm_patch_x",
            "coordinate_handle": PATCH_COORDINATE_HANDLE,
            "phase_handle": PATCH_PHASE_HANDLE,
            "domain": "x in [1,e]",
            "radius_formula": "radius=Rm*x",
            "physical_measure_formula": "1/x",
        },
    }
    return dag, outputs


def _graph_checks(dag: _DAG, outputs: Mapping[str, object]) -> dict[str, bool]:
    ids = [node.get("id") for node in dag.nodes]
    id_set = set(ids)
    unique_ids = len(ids) == len(id_set) and all(isinstance(item, str) for item in ids)
    topological = True
    for node in dag.nodes:
        args = node.get("args", [])
        if not isinstance(args, list):
            topological = False
            continue
        for arg in args:
            if not isinstance(arg, str) or arg not in id_set:
                topological = False
            elif ids.index(arg) >= ids.index(node.get("id")):
                topological = False
    rational = all(node.get("op") != "exact_constant" for node in dag.nodes)
    rational = rational and all(
        node.get("op") != "string_constant" and
        (node.get("op") != "exact_rational" or
         type(node.get("numerator")) is int and type(node.get("denominator")) is int and node.get("denominator") != 0)
        for node in dag.nodes
    )
    no_eval = all(node.get("op") not in {"numeric_evaluation", "quadrature", "solver", "norm", "recurrence"}
                  for node in dag.nodes)
    jets = []
    def scan(value: object) -> None:
        if isinstance(value, list):
            if len(value) == 3 and all(isinstance(item, str) and item in id_set for item in value):
                jets.append(value)
            for item in value:
                scan(item)
        elif isinstance(value, Mapping):
            for item in value.values():
                scan(item)
    scan(outputs)
    ordinary = bool(jets) and all(len(jet) == 3 for jet in jets)
    function_substitution = any(node.get("op") == "function_substitution" for node in dag.nodes)
    log_kernel = any(node.get("op") == "exp" and node.get("analytic") is True for node in dag.nodes)
    source_ports = [node for node in dag.nodes if node.get("op") == "source_function_port"]
    source_ports_valid = bool(source_ports) and all(
        isinstance(node.get("value"), Mapping)
        and node["value"].get("source_node_valid") is True
        and node["value"].get("source_graph_node_valid") is True
        and node["value"].get("source_node") != "missing"
        for node in source_ports
    )
    by_id = {node["id"]: node for node in dag.nodes}
    applications = [node for node in dag.nodes if node.get("op") == "source_function_application"]
    applications_valid = bool(applications) and all(
        isinstance(node.get("source_port"), str)
        and node["source_port"] in by_id
        and by_id[node["source_port"]].get("op") == "source_function_port"
        and isinstance(node.get("variable_bindings"), Mapping)
        and bool(node["variable_bindings"])
        and node.get("args") == [node["source_port"]] + list(node["variable_bindings"].values())
        and all(isinstance(value, str) and value in id_set
                for value in node["variable_bindings"].values())
        and node.get("application_is_unevaluated") is True
        for node in applications
    )
    integrals = [node for node in dag.nodes if node.get("op") == "integral"]
    integral_bindings_valid = bool(integrals) and all(
        isinstance(node.get("variable"), str) and node["variable"] in id_set
        and by_id[node["variable"]].get("op") == "exact_symbol"
        and isinstance(node.get("bounds"), list) and len(node["bounds"]) == 2
        and all(isinstance(value, str) and value in id_set for value in node["bounds"])
        and node.get("exact") is True and node.get("unevaluated") is True
        and node.get("quadrature") is False
        for node in integrals
    )
    def ancestors(ident: str) -> set[str]:
        found: set[str] = set()
        pending = [ident]
        while pending:
            current = pending.pop()
            if current in found or current not in by_id:
                continue
            found.add(current)
            pending.extend(by_id[current].get("args", []))
        return found
    integral_measure_valid = bool(integrals) and all(
        isinstance(node.get("kernel_with_measure"), str)
        and node["kernel_with_measure"] in id_set
        and by_id[node["kernel_with_measure"]].get("op") == "quotient"
        and by_id[node["kernel_with_measure"]].get("args", [None])[-1] == node.get("variable")
        and node["kernel_with_measure"] in ancestors(node.get("args", [""])[0])
        and node.get("inverse_radial_factor_in_integrand") is True
        for node in integrals
    )
    return {
        "unique_namespace_ids": unique_ids and all(node.get("namespace") == SCHEMA for node in dag.nodes),
        "closed_topological_arguments": topological,
        "rational_constants_only": rational,
        "ordinary_Z_jets_have_three_rows": ordinary,
        "no_numeric_or_solver_operation": no_eval,
        "endpoint_substitutions_are_function_nodes": function_substitution,
        "log_kernel_constructed": log_kernel,
        "scalar_ID_and_jet_list_distinct": topological and integral_bindings_valid and applications_valid,
        "expm1_entire_at_A_zero": any(node.get("op") == "expm1" and node.get("entire_at_zero") is True for node in dag.nodes),
        "qualified_source_ports_valid": source_ports_valid,
        "source_function_applications_have_explicit_closed_bindings": applications_valid,
        "integral_variables_and_bounds_are_scalar_IDs": integral_bindings_valid,
        "log_transport_inverse_s_is_in_actual_integrand": integral_measure_valid,
    }


def _status_fields(witness: Mapping[str, object], checks: Mapping[str, bool]) -> dict[str, bool]:
    fields_defined = bool(
        all(checks.values()) and witness.get("latest_ok") and witness.get("interface_ok") and witness.get("source_formula_ok")
        and witness.get("pstar_binding_ok")
        and witness.get("native_geometry_source_ok")
        and checks.get("closed_topological_arguments") and checks.get("rational_constants_only")
        and checks.get("ordinary_Z_jets_have_three_rows") and checks.get("expm1_entire_at_A_zero")
        and checks.get("qualified_source_ports_valid")
    )
    provider_identity = bool(fields_defined and witness.get("patch_receipt_witness")
                            and witness.get("archive_shape", {}).get("graph_count_exact")
                            and witness.get("c6_provider_identity"))
    density_identity = bool(provider_identity and witness.get("provider_ok") and witness.get("generic_ok"))
    inlet_identity = bool(witness.get("incoming_identity"))
    history_identity = bool(density_identity and inlet_identity and witness.get("memory_ok")
                           and checks.get("endpoint_substitutions_are_function_nodes")
                           and checks.get("log_kernel_constructed"))
    return {
        "native_C2_patch_finite_N_fields_defined": fields_defined,
        "native_C2_C5_patch_provider_identity_proved": provider_identity,
        "native_C2_C5_patch_density_identity_proved": density_identity,
        "native_C2_patch_correction_inlet_identity_proved": inlet_identity,
        "native_C2_C5_patch_history_identity_proved": history_identity,
    }


def derive(nodes: Sequence[object], source_records: Mapping[str, object], *, expected_nodes_digest: str) -> dict[str, object]:
    """Authenticate source records and return a new unevaluated finite-N DAG."""
    gates: list[str] = []
    try:
        actual_digest = _digest(nodes)
    except Exception as error:
        actual_digest = ""
        _gate(gates, "nodes_digest_error:" + type(error).__name__)
    if expected_nodes_digest != EXPECTED_NODES_DIGEST:
        _gate(gates, "expected_nodes_digest_mismatch")
    if actual_digest != EXPECTED_NODES_DIGEST:
        _gate(gates, "nodes_digest_mismatch")
    if not isinstance(nodes, Sequence) or isinstance(nodes, (str, bytes, bytearray)):
        _gate(gates, "nodes_not_sequence")
    elif len(nodes) != EXPECTED_NODE_COUNT:
        _gate(gates, "node_count_mismatch")
    if not isinstance(source_records, Mapping):
        _gate(gates, "source_records_not_mapping")
        records: Mapping[str, object] = {}
    else:
        records = source_records
    try:
        witness = _source_witness(records, gates)
    except Exception as error:
        _gate(gates, "source_inspection_error:" + type(error).__name__)
        witness = {
            "latest": {}, "provider": {}, "memory": {}, "generic": {}, "archive": {},
            "archive_shape": {}, "latest_ok": False, "interface_ok": False,
            "provider_ok": False, "generic_ok": False, "memory_ok": False,
            "patch_receipt_witness": False, "source_formula_ok": False,
            "incoming_identity": False, "retained_ok": False,
        }
    try:
        dag, dag_outputs = _build_dag(nodes, witness, gates)
    except Exception as error:
        _gate(gates, "definition_graph_construction_error:" + type(error).__name__)
        dag = _DAG()
        dag_outputs = {"error": dag.atom("source_definition_error", str(type(error).__name__))}
    dag_nodes = dag.nodes
    checks = _graph_checks(dag, dag_outputs)
    for name, passed in checks.items():
        if not passed:
            _gate(gates, "definition_graph_check_failed:" + name)
    statuses = _status_fields(witness, checks)
    status_values = dict(statuses)
    patch_identity = bool(all(status_values.values()) and not gates)
    if not patch_identity:
        _gate(gates, "patch_finite_N_source_identity_unproved")
    gates = list(dict.fromkeys(gates))
    result = {
        "schema": SCHEMA,
        "source_only": True,
        "source_only_true": True,
        "ordinary_Z_orders": list(ORDINARY_Z_ORDERS),
        "passed": bool(patch_identity and not gates),
        "nodes_digest": actual_digest,
        "expected_nodes_digest": EXPECTED_NODES_DIGEST,
        "node_count": len(nodes) if isinstance(nodes, Sequence) and not isinstance(nodes, (str, bytes, bytearray)) else None,
        **status_values,
        "patch_finite_N_source_identity": patch_identity,
        "definition_graph": {
            "namespace": SCHEMA,
            "active_graph_digest": EXPECTED_NODES_DIGEST,
            "nodes": dag_nodes,
            "node_count": len(dag_nodes),
            "outputs": dag_outputs,
        },
        "function_definitions": {
            "normalized_fields": dag_outputs.get("normalized_fields") if isinstance(dag_outputs, Mapping) else {},
            "physical_fields": dag_outputs.get("physical_fields") if isinstance(dag_outputs, Mapping) else {},
            "radial_fields": dag_outputs.get("radial_fields") if isinstance(dag_outputs, Mapping) else {},
            "geometry": dag_outputs.get("geometry") if isinstance(dag_outputs, Mapping) else {},
            "finite_N": dag_outputs.get("finite_N") if isinstance(dag_outputs, Mapping) else {},
            "histories": dag_outputs.get("histories") if isinstance(dag_outputs, Mapping) else {},
            "source_parameters": dag_outputs.get("source_parameters") if isinstance(dag_outputs, Mapping) else {},
        },
        "finite_N_source_recipe": {
            "ordinary_Z_orders": list(ORDINARY_Z_ORDERS),
            "fields": {"x": "A/N", "r": "E*expm1(A/N)", "s": "B_over_Pstar/N",
                       "E_N": "E*exp(A/N)", "V_N": "V+s"},
            "C5_density": {"m": "s", "h": "r", "k": "V*r+E*s+r*s",
                           "e": "2*V*s+s^2-E*r-r^2/2", "p": "E*r+r^2/2"},
            "C2_density": {"G": "N*expm1(A/N)", "F": "E*G",
                           "first": {"m": "B", "h": "F", "k": "V*F+E*B",
                                      "e": "2*V*B-E*F", "p": "E*F"},
                           "second": {"m": "0", "h": "0", "k": "F*B",
                                       "e": "B^2-F^2/2", "p": "F^2/2"},
                           "density": "first+second/N"},
            "identity": {"relation": "C2=N*C5", "coefficientwise_fields": list(MOMENT_NAMES),
                         "F_equals_N_r_entire_at_A_zero": True,
                         "same_B_provider_required": True,
                         "same_Pstar_units_required": True},
        },
        "section11_correction_recipe": {
            "moments": list(MOMENT_NAMES), "rates": dict(RATES),
            "formula": "H=x^-rate H_in+integral_1^x (s/x)^rate*C2rho(s)*ds/s",
            "kernel": "exp(rate*(log(s)-log(x)))",
            "coordinate": "x=R/Rm", "domain": "[1,e]", "p_memory": "exactly one",
            "reset": "only original Rminus; retain predecessor inlet",
            "derivative_identity": "x*H_x=rho-rate*H",
        },
        "qualified_source_ports": {
            "latest_report_sha256": LATEST_REPORT_SHA256,
            "C5_archive_sha256": C5_ARCHIVE_SHA256,
            "coordinate_handle": PATCH_COORDINATE_HANDLE,
            "phase_handle": PATCH_PHASE_HANDLE,
            "loop_A_handles": list(PATCH_LOOP_A),
            "loop_B_over_Pstar_handles": list(PATCH_LOOP_B),
            "C2_incoming_handles": {name: list(row) for name, row in PATCH_C2_INCOMING.items()},
            "incoming_is_correction_port": True,
            "leading_moment_port": False,
        },
        "ordinary_Z_jet_rules": {
            "orders": list(ORDINARY_Z_ORDERS),
            "product_second": "a2*b0 + 2*a1*b1 + a0*b2",
            "reciprocal": "(1/a0,-a1/a0^2,(2*a1^2-a0*a2)/a0^3)",
            "exp": "(exp(a0),exp(a0)*a1,exp(a0)*(a2+a1^2))",
            "expm1": "(expm1(a0),exp(a0)*a1,exp(a0)*(a2+a1^2))",
            "spatial_variable": "Z", "radial_variable": "x in [1,e]",
            "temporal_recursion_claimed": False,
        },
        "definition_graph_checks": checks,
        "source_prerequisite_audit": {
            "latest_native_C2_composition": witness.get("latest_ok"),
            "next_finite_N_interface": witness.get("interface_ok"),
            "accepted_provider_report": witness.get("provider_ok"),
            "accepted_generic_density_identity": witness.get("generic_ok"),
            "accepted_nonreset_memory_report": witness.get("memory_ok"),
            "patch_C5_source_formula_AST": witness.get("source_formula_ok"),
            "patch_C5_archive_shape": witness.get("archive_shape", {}).get("graph_count_exact"),
            "native_C6_patch_prefix_identity": witness.get("c6_provider_identity"),
            "native_C5_geometry_and_Z_independence": witness.get("native_geometry_source_ok"),
            "native_C6_smoothness_common_domain_premise": witness.get("source_smoothness_ok"),
            "native_C6_C5_chart_ports": witness.get("c6_chart_ports_equal"),
            "native_C6_shared_parameter_prefix": witness.get("c6_parameter_prefix_equal"),
            "accepted_Pstar_source_parameter_port": witness.get("pstar_binding_ok"),
            "patch_inlet_witness": witness.get("patch_receipt_witness"),
            "patch_correction_inlet_identity": witness.get("incoming_identity"),
            "retained_construction_data": witness.get("retained_ok"),
        },
        "native_geometry_source_witness": witness.get("native_geometry_source_witness", {}),
        "source_identity_argument": {
            "new_owner_only": True,
            "qualified_previous_DAG_refs_and_native_C5_graph_refs": True,
            "same_N_eta_delta_Pstar_phase_origin_required": True,
            "phase_origin": "hb*sc/2",
            "B_normalization": "B_over_Pstar; physical B=2*V_y remains distinct",
            "higher_derivative_source_premise": "accepted C6 source equality and genuine C6 regularity on the common original domain",
            "ordinary_Z012_is_not_higher_order_certificate": True,
            "incoming_data": "accepted C2 outgoing correction ports plus qualified native C5 postrestore inlet",
            "numeric_equality_or_interval_overlap_used": False,
            "bridge_identity_status": bool(witness.get("incoming_identity")),
            "P0_owner": "native C5 chart.modified_fields.P0",
            "geometry_owner": "native C5 chart radius/measure and native C6 inlet radius 85743",
            "prior_adapter_P0_and_geometry_alias_assumed": False,
        },
        "retained_construction_data": dict(RETAINED_CONSTRUCTION),
        "legacy_feedback_reference_alias_proved": False,
        "prior_adapter_geometry_alias_proved": False,
        "prior_adapter_P0_alias_proved": False,
        "same_density_certificate": False,
        "finite_N_history_identity": False,
        "physical_Q_y_Z1": False,
        "global_five_moment_closure": False,
        "time_order_recursion_claimed": False,
        "numeric_values_computed": False,
        "old_operators_replayed": False,
        "full_goal_complete": False,
        "runtime_model_metadata": None,
        "record_hashes": dict(RECORD_HASHES),
        "decoded_record_hashes": dict(DECODED_RECORD_HASHES),
        "source_inspection_audit": {
            "source_records_only": True,
            "source_AST_graph_saved_data_and_new_definitions_only": True,
            "nodes_hashed_only": False,
            "original_graph_nodes_inspected_not_evaluated": True,
            "new_definition_graph_constructed": bool(dag_nodes),
            "new_function_definition_graph_constructed": bool(dag_nodes),
            "new_source_consumer_calls": 1,
            "accepted_ancestor_operator_calls": 0,
            "producer_provider_builder_checker_replayed": False,
            "quadrature_or_norm_calls": 0,
            "numeric_field_evaluated": False,
            "numeric_values_computed": False,
            "old_operators_replayed": False,
            "runtime_model_metadata": None,
        },
        "scope_limits": {
            "outer_six_unresolved": ["Rh_reference", "O2_slope", "O2_axial", "O2_buffer", "O3_transition", "O3_power"],
            "global_same_density_certificate": False,
            "global_finite_N_history_identity": False,
            "endpoint_rows_are_values": False,
        },
        "unresolved_gates": gates,
    }
    return result


__all__ = [
    "EXPECTED_NODE_COUNT", "EXPECTED_NODES_DIGEST", "RECORD_HASHES",
    "DECODED_RECORD_HASHES", "SCHEMA", "derive",
]
