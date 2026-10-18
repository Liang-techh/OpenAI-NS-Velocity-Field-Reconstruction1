"""Source-only first-micro C2 versus prescribed-C6 identity gate.

The first-micro source has two live owners.  The current whole-Z owner emits
the actual ``Phi``, ``U``, six own histories, pressure datum, and ``F0`` that
the C2 projection reads through ``source_query``.  The prescribed C6 owner
constructs its own ``micro`` packet from the checked core-inlet source.  The
two source trees have the same named family and compatible equations, but
those facts do not bind the functions to the same values.  This module
therefore records the exact common algebra and fails closed until a signed
per-quantity cross-owner binding is supplied.

No generated source module is imported.  ``derive`` authenticates the supplied
graph, source text, and receipts, parses Python with :mod:`ast`, and inspects
saved JSON only.  A gzip record is required to provide its raw bytes; the raw
SHA-256 is checked before decompression.  The record contract is deliberately
small so callers can provide compact source/AST records instead of the large
saved graph archives.
"""

from __future__ import annotations

import ast
import base64
from collections.abc import Mapping, Sequence
import gzip
import hashlib
import json
from typing import Any


SCHEMA = "original_first_micro_C2_C6_source_identity_v1"
EXPECTED_NODES_DIGEST = "a2dcc782ba3033ff6ca5f54a66632ac6ef007bb0c5b9960f10567d44ef2059a3"
EXPECTED_NODE_COUNT = 205300

SOURCE_FAMILY = {
    "actual_five_defect_family_sha256": "3983d0ddb33fa85e6ab152ef7e29960f8b95aca3e1e86f1bda0882b39d825894",
    "implicit_source_sha256": "5aec111986d745459eb2e2fece291f1dc3fa7bf986529494df802c2aa2daceae",
    "datum_enclosure_sha256": "dd4040ee3ded75a56f65407fdd1cebee6aa6f559cce9d1bb57e2146a8a436010",
}

PREFIX = "lei_ren_part1_paper_compliant_current_original_"

# These are the exact parent pins used by the current original source report.
# Gzip entries are accepted as raw bytes and are never trusted through a
# decoded-text hash alone.
RECORD_HASHES = {
    PREFIX + "C2_phase_primitives.py": "84f5da3eb64096f82df88dd828b22b09bce18192a03925d91cfee0b193902d2c",
    PREFIX + "C2_phase_primitives_check.json": "76ed1af50ec799c49081ab67881b189b42e0d3137bcca770a67766d4111dff8b",
    PREFIX + "C2_phase_primitives.json.gz": "80e8f6df15edd94be6cec6aaf24528ab1eb1f06c7f7c5b61da47479e733a6bb0",
    PREFIX + "C2_density_transport.py": "8bfa32c87c73c5bedf3056c7171707ae934dafd950bd0df874dbd3b3450379d1",
    PREFIX + "C2_density_transport_check.json": "c4d2c62dd67e9dd7e171fabbdb9ba18308e420c5e0c98d01147c708840b8f583",
    PREFIX + "C2_density_transport.json.gz": "0f347824bea798ea692eca7967d7acadba49fe4336989bea1a5c86c922260877",
    PREFIX + "prescribed_loop_C5_functions.py": "69d6391120000488472b7961b4bff1901bff872def95d183c4ee056b0948e17d",
    PREFIX + "prescribed_loop_C5_functions_check.json": "2c77e80df04fae615f8ea909fc7b606bfccfa27659cb67744ae8d84020b1a9e7",
    PREFIX + "prescribed_loop_C5_functions.json.gz": "8419fdc338dc746099c7cad119ea61a8e928652319d5caf5314f8a730c533caa",
    PREFIX + "whole_Z_bridge_source.py": "df8b204d023c1d8762da16d07190601ebd1e11447077b01a734bba8baa595a6a",
    PREFIX + "whole_Z_bridge_source_check.json": "91d9ba8119d9dba85322ea33791ac3126c22ad85dfd94bf9a961c49fc18268a3",
    PREFIX + "whole_Z_bridge_source.json.gz": "91f3f68b5e6ef98e851134e680777cb1d29e2410a8afc6f9ae4380e70b129c3b",
    PREFIX + "whole_Z_bridge_correlated_functions.py": "bb29100e671b6a1927f9897385056645a17d859108e537eec760a3904d66f832",
    PREFIX + "whole_Z_bridge_correlated_functions_check.json": "bbf33386a8a591eef946ca5f3f262dd92ac6dd1490bc8bfc7421ee1c6eca49da",
    PREFIX + "whole_Z_bridge_correlated_functions.json.gz": "5c0b01b829b3dd65c4886709ecade889f2b21c38de21247f9cb36e4cb8200060",
    PREFIX + "prescribed_bridge_C6_functions.py": "8ab2c26092db06453ed5a6047d8905e37fdad17eaf9c4947048498551425c605",
    PREFIX + "prescribed_bridge_C6_functions_check.json": "9a0969a8cd4a5989dbd6272d6dd223533400113824f9352b4e57a83196bff61e",
    PREFIX + "core_inlet_axial6_generic_source.py": "29738d46c6eb5897f0fbeefaec741d0eda4b97baed57b43ce9a7890a7dc03890",
    PREFIX + "core_inlet_axial6_generic_source_check.json": "9d8b746ef78502956cbb2b53f96d452be636eecbf645d021aaca3f9c2e21885d",
    PREFIX + "prescribed_bridge_C6_functions.json.gz": "5df06c4312efc01073516023fd886c071f439c11de30d3411b2e358bf76abd92",
    PREFIX + "Rm_generic_inputs.py": "62df6a6513b6e9c521cb6cf7cff25710dec17bf04f7b23f4bc0f15dbc1a57a88",
    PREFIX + "Rm_generic_inputs_check.json": "ba4c448e316cffd08232b41cbe6783e4b7cd36705c657a27117d2e9fd6b8629e",
    PREFIX + "Rm_generic_inputs.json.gz": "0d84df4c85f2bba91b33c27c567b0c62777a181df731a0c87212ba2a6d7e3e9f",
    PREFIX + "reference_restore_finite_N.py": "f338527b4139f549f2f76faaaa65f94923c423bfbc0d133151a7aca0dfa3a47d",
    PREFIX + "reference_restore_finite_N_check.json": "1a0a41eac31099b52fc76ce57754aadfc43459c7b7d0cd9a29e83dc35a236c4a",
    PREFIX + "reference_restore_finite_N.json.gz": "319b285ac9a4c2b9cc5dd7e2896a3ca5dadcc2aade3ca8768838b2f528c47a24",
    "lei_ren_part1_paper_compliant_core_coefficient_rebuild.py": "d37bfa6fe89b8955e17b38a00e1c8ae53b6725bb5a73c4648067737ce3be990f",
    "lei_ren_part1_paper_compliant_core_coefficient_rebuild_check.json": "02790017aeddc309624dcff46f2d4a2cd9d557075cc0e88096ff0c63cd9b65d0",
    "lei_ren_part1_paper_compliant_current_core_common_fixed_point.py": "dd44110f48fb40fee715b358801dd6b8944bda5e6d5c0054eefa63c66f82255f",
    "lei_ren_part1_paper_compliant_current_core_common_fixed_point_check.json": "a1dee74791b2f9c4a9cafcad62e573baf788c7b6f6b39cef7cde90d23a2ff672",
    "lei_ren_part1_paper_compliant_current_core_nonlinear_operator.py": "d09f1bd0e1a55b14fd8a4006fd38f0991ea7787b88bc55ebb136c576736abf9f",
    "lei_ren_part1_paper_compliant_current_core_nonlinear_operator_check.json": "a1f9a3632a7765f6bc8d38ac1bd25c28850e128a8b25adaa0d42f769f3704d62",
    PREFIX + "n1_analytic_core_lower_pressure_C3_functions.py": "14178f73828e278d9368721b93b1efa1687dd190f837ca92a0d78a10ecaeb5c0",
    PREFIX + "n1_analytic_core_lower_pressure_C3_functions_check.json": "e3b52bd2e0331b3d8956bad328dc413b276cd8cfee251c198f18ba9f5fb703f5",
    "lei_ren_part1_paper_compliant_current_original_Rm_m_density_transport_source_interfaces.json": "6f1be1e1f86667ad24fd58139289ca10bf77e25ea0eb6d2f713e9c11e602715b",
    "lei_ren_part1_paper_compliant_current_original_Rm_physical_Qy_Z_calculus.json": "8903b59267bae01d5182bf2dcc84ed3c02c44eb8c7d3dd9d7b14af67ac092746",
}

DECODED_RECORD_HASHES = {
    PREFIX + "C2_phase_primitives.json.gz": "e614b601f5ee489bb713ba4e3ae9cd297ed318f23aa3ac2f9f0c34124fa1fff2",
    PREFIX + "C2_density_transport.json.gz": "2a9e4a55b53110944beb2883ed6ea11368ee3f695af9185c6ca72ffe8f966dff",
    PREFIX + "prescribed_loop_C5_functions.json.gz": "f4825741f568c001f0f4658dfc3ce8eec3dc0466db976c29f2f4b609b8ceadc2",
    PREFIX + "whole_Z_bridge_source.json.gz": "cbbf060c08dad8502ac3c5d837c6fa811dafb3d400ba5534f1849dd0303b650a",
    PREFIX + "whole_Z_bridge_correlated_functions.json.gz": "1eb004c0d14ce861fd8f9a72fa7cc7c03014f9562464e3ce359892924ad96285",
    PREFIX + "prescribed_bridge_C6_functions.json.gz": "f361ab39e7599f5ed770df630bccadc4056b3e7ad6c0e1de32c975963f070298",
    PREFIX + "Rm_generic_inputs.json.gz": "512f7d9d94e734f3904b16ba24bc0ef5c5df7f10c637bbd35e00ec826ae86400",
    PREFIX + "reference_restore_finite_N.json.gz": "94ba3cbbe09ddaa69ca96f9ecc551871b6884b20724e2015b086d6db8abdfd07",
}

ARCHIVE_RECORDS = frozenset(name for name in RECORD_HASHES if name.endswith(".json.gz"))
REPORT_RECORDS = frozenset({
    "lei_ren_part1_paper_compliant_current_original_Rm_m_density_transport_source_interfaces.json",
    "lei_ren_part1_paper_compliant_current_original_Rm_physical_Qy_Z_calculus.json",
})

_FOUNDATION_RECORDS = frozenset({
    PREFIX + "core_inlet_axial6_generic_source.py",
    PREFIX + "core_inlet_axial6_generic_source_check.json",
})

_CORE_PROJECTION_RECORDS = frozenset({
    "lei_ren_part1_paper_compliant_core_coefficient_rebuild.py",
    "lei_ren_part1_paper_compliant_core_coefficient_rebuild_check.json",
    "lei_ren_part1_paper_compliant_current_core_common_fixed_point.py",
    "lei_ren_part1_paper_compliant_current_core_common_fixed_point_check.json",
    "lei_ren_part1_paper_compliant_current_core_nonlinear_operator.py",
    "lei_ren_part1_paper_compliant_current_core_nonlinear_operator_check.json",
    PREFIX + "n1_analytic_core_lower_pressure_C3_functions.py",
    PREFIX + "n1_analytic_core_lower_pressure_C3_functions_check.json",
})

_C2_PHASE = PREFIX + "C2_phase_primitives.py"
_C2_PHASE_CHECK = PREFIX + "C2_phase_primitives_check.json"
_C2_DENSITY = PREFIX + "C2_density_transport.py"
_C2_DENSITY_CHECK = PREFIX + "C2_density_transport_check.json"
_C5 = PREFIX + "prescribed_loop_C5_functions.py"
_C5_CHECK = PREFIX + "prescribed_loop_C5_functions_check.json"
_BRIDGE = PREFIX + "whole_Z_bridge_source.py"
_BRIDGE_CHECK = PREFIX + "whole_Z_bridge_source_check.json"
_CORRELATED = PREFIX + "whole_Z_bridge_correlated_functions.py"
_CORRELATED_CHECK = PREFIX + "whole_Z_bridge_correlated_functions_check.json"
_C6 = PREFIX + "prescribed_bridge_C6_functions.py"
_C6_CHECK = PREFIX + "prescribed_bridge_C6_functions_check.json"
_FOUNDATION = PREFIX + "core_inlet_axial6_generic_source.py"
_FOUNDATION_CHECK = PREFIX + "core_inlet_axial6_generic_source_check.json"
_CORE_REBUILD = "lei_ren_part1_paper_compliant_core_coefficient_rebuild.py"
_CORE_REBUILD_CHECK = "lei_ren_part1_paper_compliant_core_coefficient_rebuild_check.json"
_CORE_COMMON = "lei_ren_part1_paper_compliant_current_core_common_fixed_point.py"
_CORE_COMMON_CHECK = "lei_ren_part1_paper_compliant_current_core_common_fixed_point_check.json"
_CORE_OPERATOR = "lei_ren_part1_paper_compliant_current_core_nonlinear_operator.py"
_CORE_OPERATOR_CHECK = "lei_ren_part1_paper_compliant_current_core_nonlinear_operator_check.json"
_N1_ANALYTIC = PREFIX + "n1_analytic_core_lower_pressure_C3_functions.py"
_N1_ANALYTIC_CHECK = PREFIX + "n1_analytic_core_lower_pressure_C3_functions_check.json"
_N1_LOWER = PREFIX + "n1_lower_collar_functions.py"
_N1_LOWER_CHECK = PREFIX + "n1_lower_collar_functions_check.json"
_N1_CORE = PREFIX + "n1_core_fixed_point_functions.py"
_N1_CORE_CHECK = PREFIX + "n1_core_fixed_point_functions_check.json"
_N1_ANALYTIC_EXPECTED_INPUT_PINS = {
    _N1_CORE: "d71ff19a68cfed8ace9a6857ffeca860ecee90d8d89b35b2292a7310e180fb59",
    _N1_CORE_CHECK: "45809af92d5bfff3f518c5623f6b53f9175dc2321cf654008671d4b9bd832ead",
    _N1_ANALYTIC: "14178f73828e278d9368721b93b1efa1687dd190f837ca92a0d78a10ecaeb5c0",
    _N1_ANALYTIC_CHECK: "e3b52bd2e0331b3d8956bad328dc413b276cd8cfee251c198f18ba9f5fb703f5",
    _N1_LOWER: "d00dfe285735844eabc985f99d94190caa1cfd187f668eb14a2686055cf175fa",
    _N1_LOWER_CHECK: "1a7b6b370b5ecf7d6b7e97001e0311e9525698ecf3389fc3da608af308a13ba9",
}
_GENERIC = PREFIX + "Rm_generic_inputs.py"
_GENERIC_CHECK = PREFIX + "Rm_generic_inputs_check.json"
_RESTORE = PREFIX + "reference_restore_finite_N.py"
_RESTORE_CHECK = PREFIX + "reference_restore_finite_N_check.json"
_LATEST_REPORT = "lei_ren_part1_paper_compliant_current_original_Rm_m_density_transport_source_interfaces.json"
_PARENT_REPORT = "lei_ren_part1_paper_compliant_current_original_Rm_physical_Qy_Z_calculus.json"

_REQUIRED_RECORDS = frozenset(set(RECORD_HASHES) - ARCHIVE_RECORDS - _FOUNDATION_RECORDS - _CORE_PROJECTION_RECORDS)
FIRST_MICRO_DOMAIN = ("sc/2", "1")


def canonical_digest(value: Any) -> str:
    """Hash JSON-like graph data without evaluating any operation node."""

    raw = json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True, allow_nan=False).encode("utf-8")
    return hashlib.sha256(raw).hexdigest()


def _gate(gates: list[dict[str, Any]], name: str, **details: Any) -> None:
    row = {"gate": name, **details}
    if row not in gates:
        gates.append(row)


def _lookup(source_records: Mapping[str, Any], name: str) -> Any:
    if not isinstance(source_records, Mapping):
        return None
    return source_records.get(name, source_records.get("experiments/root_st073/" + name))


def _raw_bytes(value: Any, encoded: Any = None) -> bytes | None:
    if isinstance(value, bytes):
        return value
    if isinstance(value, (bytearray, memoryview)):
        return bytes(value)
    if isinstance(encoded, str):
        try:
            return base64.b64decode(encoded, validate=True)
        except (ValueError, TypeError):
            return None
    return None


def _record(
    source_records: Mapping[str, Any],
    name: str,
    gates: list[dict[str, Any]],
    *,
    optional: bool = False,
) -> tuple[str | None, dict[str, Any]]:
    """Authenticate one record and return decoded text plus static metadata.

    A ``.json.gz`` record must carry raw bytes directly or as
    ``raw_bytes_base64``.  The digest check is made on those bytes before
    ``gzip.decompress`` is called.  Plain records may use UTF-8
    ``source_text`` or the equivalent raw bytes.
    """

    record = _lookup(source_records, name)
    if not isinstance(record, Mapping):
        if not optional:
            _gate(gates, "source_record_missing", record=name)
        return None, {"present": False, "archive": name.endswith(".json.gz")}
    expected_path = "experiments/root_st073/" + name
    if record.get("path") != expected_path:
        _gate(gates, "source_record_path_changed", record=name, expected=expected_path)
        return None, {"present": False, "archive": name.endswith(".json.gz")}
    expected = RECORD_HASHES[name]
    archive = name.endswith(".json.gz")
    encoded_raw = record.get("raw_bytes_base64")
    raw = _raw_bytes(record.get("raw_bytes"), encoded_raw)
    if archive:
        if raw is None:
            if encoded_raw is not None:
                _gate(gates, "gzip_raw_bytes_base64_invalid", record=name)
            else:
                _gate(gates, "gzip_raw_bytes_required", record=name)
            return None, {"present": False, "archive": True}
        actual = hashlib.sha256(raw).hexdigest()
        if record.get("sha256") != expected or actual != expected:
            _gate(gates, "source_record_hash_changed", record=name, expected_sha256=expected)
            return None, {"present": False, "archive": True}
        try:
            decoded = gzip.decompress(raw)
            decoded_hash = hashlib.sha256(decoded).hexdigest()
            text = decoded.decode("utf-8")
        except (OSError, UnicodeDecodeError) as exc:
            _gate(gates, "gzip_decode_failed_after_raw_hash", record=name, error=str(exc))
            return None, {"present": False, "archive": True}
        supplied_text = record.get("source_text")
        if supplied_text is not None and supplied_text != text:
            _gate(gates, "gzip_decoded_source_text_mismatch", record=name)
            return None, {"present": False, "archive": True}
        expected_decoded = DECODED_RECORD_HASHES[name]
        if decoded_hash != expected_decoded or record.get("decoded_sha256", expected_decoded) != expected_decoded:
            _gate(gates, "decoded_source_hash_changed", record=name, expected_sha256=expected_decoded)
            return None, {"present": False, "archive": True}
        return text, {"present": True, "archive": True, "raw_sha256": actual, "decoded_sha256": decoded_hash}
    if raw is None:
        text_value = record.get("source_text")
        if not isinstance(text_value, str):
            _gate(gates, "plain_source_text_required", record=name)
            return None, {"present": False, "archive": False}
        raw = text_value.encode("utf-8")
    else:
        try:
            text_value = raw.decode("utf-8")
        except UnicodeDecodeError as exc:
            _gate(gates, "plain_source_utf8_invalid", record=name, error=str(exc))
            return None, {"present": False, "archive": False}
        supplied = record.get("source_text")
        if supplied is not None and supplied != text_value:
            _gate(gates, "plain_source_text_mismatch", record=name)
            return None, {"present": False, "archive": False}
    actual = hashlib.sha256(raw).hexdigest()
    if record.get("sha256") != expected or actual != expected:
        _gate(gates, "source_record_hash_changed", record=name, expected_sha256=expected)
        return None, {"present": False, "archive": False}
    return text_value, {"present": True, "archive": False, "raw_sha256": actual}


def _json(text: str | None, name: str, gates: list[dict[str, Any]]) -> Any:
    if text is None:
        return None
    try:
        return json.loads(text, parse_int=str)
    except (TypeError, json.JSONDecodeError) as exc:
        _gate(gates, "saved_json_invalid", record=name, error=str(exc))
        return None


def _ast(text: str | None, name: str, gates: list[dict[str, Any]]) -> ast.AST | None:
    if text is None:
        return None
    try:
        return ast.parse(text, filename=name)
    except SyntaxError as exc:
        _gate(gates, "source_ast_invalid", record=name, error=str(exc))
        return None


def _dotted(node: ast.AST | None) -> str | None:
    if isinstance(node, ast.Name):
        return node.id
    if isinstance(node, ast.Attribute):
        left = _dotted(node.value)
        return None if left is None else left + "." + node.attr
    return None


def _call(node: ast.AST | None, dotted: str) -> ast.Call | None:
    return node if isinstance(node, ast.Call) and _dotted(node.func) == dotted else None


def _constant(node: ast.AST | None, value: Any) -> bool:
    return isinstance(node, ast.Constant) and node.value == value


def _subscript(node: ast.AST | None, base: str, key: str) -> bool:
    if not isinstance(node, ast.Subscript) or _dotted(node.value) != base:
        return False
    return isinstance(node.slice, ast.Constant) and node.slice.value == key


def _method(tree: ast.AST | None, class_name: str, method_name: str) -> ast.AST | None:
    if tree is None:
        return None
    for node in ast.walk(tree):
        if isinstance(node, ast.ClassDef) and node.name == class_name:
            for child in node.body:
                if isinstance(child, (ast.FunctionDef, ast.AsyncFunctionDef)) and child.name == method_name:
                    return child
    return None


def _function(tree: ast.AST | None, name: str) -> ast.AST | None:
    if tree is None:
        return None
    for node in ast.walk(tree):
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and node.name == name:
            return node
    return None


def _assignment(node: ast.AST | None, name: str) -> list[ast.AST]:
    values: list[ast.AST] = []
    if node is None:
        return values
    for child in ast.walk(node):
        if isinstance(child, ast.Assign):
            for target in child.targets:
                if isinstance(target, ast.Name) and target.id == name:
                    values.append(child.value)
        elif isinstance(child, ast.AnnAssign) and isinstance(child.target, ast.Name) and child.target.id == name:
            values.append(child.value)
    return values


def _tuple_assignment_values(
    node: ast.AST | None,
    names: tuple[str, ...],
) -> tuple[ast.AST, ...] | None:
    """Read one exact tuple-target assignment with the requested names."""

    if node is None:
        return None
    for child in ast.walk(node):
        if not isinstance(child, ast.Assign) or len(child.targets) != 1:
            continue
        target = child.targets[0]
        value = child.value
        if not isinstance(target, ast.Tuple) or not isinstance(value, ast.Tuple):
            continue
        if len(target.elts) != len(names) or len(value.elts) != len(names):
            continue
        if not all(
            isinstance(item, ast.Name) and item.id == expected
            for item, expected in zip(target.elts, names)
        ):
            continue
        return tuple(value.elts)
    return None


def _returned_dict_keyword(node: ast.AST | None, key: str) -> ast.AST | None:
    """Read one keyword from a direct ``return dict(...)`` expression only."""

    if node is None:
        return None
    for child in ast.walk(node):
        if not isinstance(child, ast.Return) or not isinstance(child.value, ast.Call):
            continue
        if _dotted(child.value.func) != "dict":
            continue
        for keyword in child.value.keywords:
            if keyword.arg == key:
                return keyword.value
    return None


def _attribute_assignment(node: ast.AST | None, attribute: str) -> list[ast.AST]:
    values: list[ast.AST] = []
    if node is None:
        return values
    for child in ast.walk(node):
        if isinstance(child, ast.Assign):
            if any(isinstance(target, ast.Attribute) and target.attr == attribute for target in child.targets):
                values.append(child.value)
    return values


def _dict_keyword(node: ast.AST | None, variable: str, key: str) -> ast.AST | None:
    for value in _assignment(node, variable):
        if isinstance(value, ast.Call) and _dotted(value.func) == "dict":
            for keyword in value.keywords:
                if keyword.arg == key:
                    return keyword.value
        if isinstance(value, ast.Dict):
            for k, v in zip(value.keys, value.values):
                if _constant(k, key):
                    return v
    return None


def _attribute_dict_keyword(node: ast.AST | None, attribute: str, key: str) -> ast.AST | None:
    for value in _attribute_assignment(node, attribute):
        if isinstance(value, ast.Call) and _dotted(value.func) == "dict":
            for keyword in value.keywords:
                if keyword.arg == key:
                    return keyword.value
        if isinstance(value, ast.Dict):
            for k, v in zip(value.keys, value.values):
                if _constant(k, key):
                    return v
    return None


def _call_keyword(node: ast.AST | None, call_name: str, key: str) -> ast.AST | None:
    """Read one keyword from a direct call, including a returned ``dict``."""

    if node is None:
        return None
    for child in ast.walk(node):
        if isinstance(child, ast.Call) and _dotted(child.func) == call_name:
            for keyword in child.keywords:
                if keyword.arg == key:
                    return keyword.value
    return None


def _constructor_keyword(
    node: ast.AST | None,
    variable: str,
    constructor: str,
    key: str,
) -> ast.AST | None:
    """Read one keyword from one named constructor assignment only."""

    for value in _assignment(node, variable):
        if isinstance(value, ast.Call) and _dotted(value.func) == constructor:
            for keyword in value.keywords:
                if keyword.arg == key:
                    return keyword.value
    return None


def _contains_call(node: ast.AST | None, dotted: str) -> bool:
    return any(_call(child, dotted) is not None for child in ast.walk(node)) if node is not None else False


def _contains_subscript(node: ast.AST | None, base: str, key: str) -> bool:
    return any(_subscript(child, base, key) for child in ast.walk(node)) if node is not None else False


def _contains_name(node: ast.AST | None, name: str) -> bool:
    return any(isinstance(child, ast.Name) and child.id == name for child in ast.walk(node)) if node is not None else False


def _contains_constant(node: ast.AST | None, value: Any) -> bool:
    return any(_constant(child, value) for child in ast.walk(node)) if node is not None else False


def _contains_attribute(node: ast.AST | None, dotted: str) -> bool:
    return any(_dotted(child) == dotted for child in ast.walk(node)) if node is not None else False


def _contains_is_not(node: ast.AST | None, left: str, right: str) -> bool:
    if node is None:
        return False
    return any(
        isinstance(child, ast.Compare)
        and len(child.ops) == 1
        and isinstance(child.ops[0], ast.IsNot)
        and _dotted(child.left) == left
        and len(child.comparators) == 1
        and _dotted(child.comparators[0]) == right
        for child in ast.walk(node)
    )


def _contains_not_equal(node: ast.AST | None, left: str, right: str) -> bool:
    if node is None:
        return False
    return any(
        isinstance(child, ast.Compare)
        and len(child.ops) == 1
        and isinstance(child.ops[0], ast.NotEq)
        and _dotted(child.left) == left
        and len(child.comparators) == 1
        and _dotted(child.comparators[0]) == right
        for child in ast.walk(node)
    )


def _contains_subscript_is_not(node: ast.AST | None, base: str, key: str, right: str) -> bool:
    if node is None:
        return False
    return any(
        isinstance(child, ast.Compare)
        and len(child.ops) == 1
        and isinstance(child.ops[0], ast.IsNot)
        and _subscript(child.left, base, key)
        and len(child.comparators) == 1
        and _dotted(child.comparators[0]) == right
        for child in ast.walk(node)
    )


def _contains_subscript_call(node: ast.AST | None, dotted: str, key: str) -> bool:
    if node is None:
        return False
    return any(
        isinstance(child, ast.Subscript)
        and isinstance(child.slice, ast.Constant)
        and child.slice.value == key
        and _call(child.value, dotted) is not None
        for child in ast.walk(node)
    )


def _contains_subscript_assignment(node: ast.AST | None, base: str, key: str) -> bool:
    if node is None:
        return False
    for child in ast.walk(node):
        targets: Sequence[ast.AST] = ()
        if isinstance(child, ast.Assign):
            targets = child.targets
        elif isinstance(child, ast.AnnAssign):
            targets = (child.target,)
        if any(_subscript(target, base, key) for target in targets):
            return True
    return False


def _contains_type_identity(node: ast.AST | None, argument: str, qualified_type: str) -> bool:
    """Require the source's exact ``type(x) is module.Type`` guard."""

    if node is None:
        return False
    for child in ast.walk(node):
        if not isinstance(child, ast.Compare) or len(child.ops) != 1 or not isinstance(child.ops[0], (ast.Is, ast.IsNot)):
            continue
        if not isinstance(child.left, ast.Call) or _dotted(child.left.func) != "type" or len(child.left.args) != 1:
            continue
        if child.left.keywords or _dotted(child.left.args[0]) != argument or len(child.comparators) != 1:
            continue
        if _dotted(child.comparators[0]) == qualified_type:
            return True
    return False


def _contains_exact_call_text(node: ast.AST | None, qualified: str, expression: str) -> bool:
    """Match one source call by its unparsed AST, including every argument."""

    if node is None:
        return False
    for child in ast.walk(node):
        if isinstance(child, ast.Call) and _dotted(child.func) == qualified and ast.unparse(child) == expression:
            return True
    return False


def _expr_dump(node: ast.AST | None) -> str | None:
    return None if node is None else ast.dump(node, include_attributes=False)


def _expr_equal(node: ast.AST | None, expression: str) -> bool:
    """Compare literal expression ASTs; the RHS is parsed in eval mode only."""

    if node is None:
        return False
    try:
        expected = ast.parse(expression, mode="eval").body
    except SyntaxError:
        return False
    return _expr_dump(node) == _expr_dump(expected)


def _require(condition: bool, gates: list[dict[str, Any]], name: str, **details: Any) -> bool:
    if not condition:
        _gate(gates, name, **details)
    return condition


def _receipt(
    text: str | None,
    name: str,
    required_paths: Sequence[Sequence[str]],
    gates: list[dict[str, Any]],
) -> tuple[bool, Any]:
    value = _json(text, name, gates)
    if not isinstance(value, Mapping):
        _gate(gates, "accepted_source_receipt_missing", record=name)
        return False, value
    ok = value.get("all_passed") is True
    if not ok:
        _gate(gates, "accepted_source_receipt_header_missing", record=name)
    family = value.get("source_family")
    if isinstance(family, Mapping):
        for field, expected in SOURCE_FAMILY.items():
            if family.get(field) != expected:
                _gate(gates, "source_family_changed", record=name, field=field, expected=expected)
                ok = False
    elif family is not None and family != SOURCE_FAMILY.get("actual_five_defect_family_sha256"):
        _gate(gates, "source_family_changed", record=name, expected=SOURCE_FAMILY)
        ok = False
    for path in required_paths:
        current: Any = value
        for key in path:
            current = current.get(key) if isinstance(current, Mapping) else None
        if current is not True:
            _gate(gates, "accepted_source_receipt_flag_missing", record=name, path=".".join(path))
            ok = False
    return ok, value


def _saved_path(value: Any, path: Sequence[str]) -> Any:
    """Read a saved JSON path without assuming an untrusted container shape."""

    current = value
    for key in path:
        if not isinstance(current, Mapping):
            return None
        current = current.get(key)
    return current


def _c2_witness(tree: ast.AST | None, gates: list[dict[str, Any]]) -> dict[str, Any]:
    leading = _method(tree, "CurrentC2PhasePrimitives", "leading_source")
    packet = _method(tree, "CurrentC2PhasePrimitives", "source_packet")
    leading_ok = _require(
        leading is not None and _contains_call(leading, "self.outer.owner.source_query"),
        gates,
        "c2_leading_source_does_not_delegate_whole_z_source_query",
    )
    raw_ok = _require(
        packet is not None
        and _contains_subscript(packet, "packet", "original_generic_source")
        and _contains_subscript(packet, "packet", "original_full_source_quotients"),
        gates,
        "c2_projection_missing_actual_generic_or_quotient_source",
    )
    rows = _assignment(packet, "rows")
    rows_value = rows[-1] if rows else None
    e_node = _dict_keyword(packet, "rows", "E")
    v_node = _dict_keyword(packet, "rows", "V")
    p0_node = _dict_keyword(packet, "result", "exact_common_P0_axial5")
    e_ok = _subscript(e_node, "raw", "common_velocity_E_axial5")
    v_ok = _subscript(v_node, "raw", "common_velocity_V_axial5")
    p0_ok = _dotted(p0_node) == "op.P0"
    _require(e_ok, gates, "c2_E_root_not_common_velocity_E_axial5")
    _require(v_ok, gates, "c2_V_root_not_common_velocity_V_axial5")
    _require(p0_ok, gates, "c2_P0_not_exact_live_op_P0")
    projection_ok = leading_ok and raw_ok and e_ok and v_ok and p0_ok and rows_value is not None
    density_ok = _function(tree, "build") is not None
    return {
        "passed": bool(projection_ok),
        "source_owner": "CurrentC2PhasePrimitives.source_packet -> WholeZCorrelatedBridgeFunctions.source_query",
        "actual_paths": {
            "Phi": "whole-Z actual prepared data phi0 / live actual source state",
            "U": "whole-Z actual prepared data V0 / live actual source state",
            "own_moments": "original_full_source_quotients and C2 source transport",
            "E": "original_generic_source.common_velocity_E_axial5",
            "V": "original_generic_source.common_velocity_V_axial5",
            "P0": "exact_common_P0_axial5 = op.P0",
        },
        "density_transport_ast_present": density_ok,
        "rows_descriptor_checked": rows_value is not None,
        "operator_evaluation": False,
    }


def _whole_z_witness(tree: ast.AST | None, gates: list[dict[str, Any]]) -> dict[str, Any]:
    owner = _method(tree, "OriginalWholeZBridgeSource", "owner")
    inputs = _method(tree, "OriginalWholeZBridgeSource", "inputs")
    actual_core = _require(
        owner is not None
        and bool(_assignment(owner, "prepared"))
        and bool(_assignment(owner, "core"))
        and bool(_assignment(owner, "flow"))
        and _contains_call(owner, "self.actual.prepare")
        and _contains_call(owner, "self.actual.packet")
        and _contains_call(owner, "macro.MacroFlow"),
        gates,
        "whole_z_actual_core_or_flow_owner_not_bound",
    )
    source_inputs = _require(
        inputs is not None and _contains_call(inputs, "self.admission.pressure.normalized_jets"),
        gates,
        "whole_z_pressure_normalized_jets_not_source_bound",
    )
    proof = _dict_keyword(owner, "proof", "fresh_actual_core")
    g_actual = _dict_keyword(owner, "proof", "actual_G_not_Gbar_used_in_flow")
    pressure = _dict_keyword(owner, "proof", "original_pressure_Z0_through_Z6")
    proof_ok = proof is not None and _constant(g_actual, True) and pressure is not None
    _require(proof_ok, gates, "whole_z_actual_core_pressure_or_G_proof_missing")
    f0_ok = inputs is not None and _dict_keyword(inputs, "result", "F0_ratios") is not None and _dict_keyword(inputs, "result", "F0_squared_ratios") is not None
    _require(f0_ok, gates, "whole_z_F0_gradient_source_not_bound")
    return {
        "passed": bool(actual_core and source_inputs and proof_ok and f0_ok),
        "owner": "OriginalWholeZBridgeSource.owner",
        "actual_core_packet": "self.actual.packet(Z, 0, 'first')",
        "actual_state_preparation": "self.actual.prepare(Z)",
        "common_pressure_source": "self.admission.pressure.normalized_jets -> p0",
        "common_F0_source": "amplitude.gradient_jets -> F0 ratios and F0 squared ratios",
        "actual_G_used": True,
        "whole_Z_pressure_rows": "original_pressure_Z0_through_Z6",
        "numeric_values": False,
    }


def _correlated_witness(tree: ast.AST | None, gates: list[dict[str, Any]]) -> dict[str, Any]:
    query = _method(tree, "WholeZCorrelatedBridgeFunctions", "source_query")
    packet = _assignment(query, "packet")
    n_node = _assignment(query, "n")
    bind = _contains_subscript_assignment(query, "packet", "original_C_B_E_C0_Z")
    ok = _require(
        query is not None
        and bool(packet)
        and bool(n_node)
        and _contains_subscript(query, "packet", "original_generic_source")
        and _contains_subscript(query, "packet", "original_C_B_E_C0_Z")
        and _contains_subscript(query, "packet", "original_eta_log"),
        gates,
        "correlated_source_query_does_not_expose_actual_source_packet",
    )
    return {
        "passed": bool(ok),
        "owner": "WholeZCorrelatedBridgeFunctions.source_query",
        "entry": "source_query",
        "actual_generic_source": "packet['original_generic_source']",
        "actual_C_B_E_C0_Z": "C/B/E copied from actual_generic_source_numerators",
        "eta_log": "self.switch_source.eta_log",
        "descriptor_assignment_present": bool(bind),
    }


def _c6_actual_witness(tree: ast.AST | None, gates: list[dict[str, Any]]) -> dict[str, Any]:
    actual = _function(tree, "actual_fields")
    histories = _function(tree, "histories")
    generic = _function(tree, "generic_fields")
    init = _method(tree, "CurrentOriginalPrescribedBridgeC6Functions", "__init__")
    actual_present = all(x is not None for x in (actual, histories, generic, init))
    _require(actual_present, gates, "c6_actual_fields_histories_or_generic_fields_missing")
    angular = _assignment(actual, "angular")
    phi = _assignment(actual, "phi")
    drive = _assignment(actual, "full_drive")
    vy = _assignment(actual, "Vy")
    u = _assignment(actual, "U")
    phiy = _assignment(actual, "phiy")
    # The semantic checks below inspect the actual call tree, including the
    # negative half and the double negative in U, rather than relying on a
    # field-name or argument-count heuristic.
    angular_ok = (
        bool(angular)
        and _contains_subscript(angular[-1], "known", "D_over_R")
        and _contains_name(angular[-1], "chi")
        and _contains_name(angular[-1], "radius")
        and _contains_name(angular[-1], "jacobian")
        and _contains_constant(angular[-1], "-1/2")
    )
    phi_ok = bool(phi) and _contains_subscript(phi[-1], "initial", "phi") and _contains_call(phi[-1], "a.exp")
    drive_ok = (
        bool(drive)
        and _call(drive[-1], "a.add") is not None
        and _contains_name(drive[-1], "radius")
        and _contains_name(drive[-1], "Pstar")
        and _contains_name(drive[-1], "F0")
        and _contains_name(drive[-1], "known")
    )
    vy_ok = bool(vy) and _contains_call(vy[-1], "g.neg") and _contains_name(vy[-1], "chi") and _contains_name(vy[-1], "ratio") and _contains_name(vy[-1], "full_drive")
    u_ok = bool(u) and _call(u[-1], "a.sub") is not None and _contains_call(u[-1], "a.integral") and _contains_call(u[-1], "g.neg") and _contains_name(u[-1], "jacobian")
    phiy_ok = bool(phiy) and _contains_subscript(phiy[-1], "known", "D_over_R") and _contains_name(phiy[-1], "phi") and _contains_name(phiy[-1], "chi") and _contains_name(phiy[-1], "radius") and _contains_constant(phiy[-1], "-1/2")
    _require(angular_ok, gates, "c6_phi_angular_equation_missing")
    _require(phi_ok, gates, "c6_phi_actual_initial_data_or_exponential_missing")
    _require(drive_ok, gates, "c6_full_signed_drive_missing_hydro_pressure_swirl_terms")
    _require(vy_ok, gates, "c6_Vy_missing_full_axial_minus_sign")
    _require(u_ok, gates, "c6_U_missing_double_negative_jacobian_integral")
    _require(phiy_ok, gates, "c6_phi_y_missing_negative_half_radius_factor")

    fields = {name: _dict_keyword(generic, "fields", name) for name in ("E", "V", "m", "h", "k", "e", "p", "P0", "pressure")}
    radial = {name: _dict_keyword(generic, "radial", name) for name in ("E_y", "V_y")}
    inverse_nodes = _assignment(generic, "inverse")
    inverse_ok = bool(inverse_nodes) and _call(inverse_nodes[-1], "g.quotient") is not None and _contains_name(inverse_nodes[-1], "Pstar")
    w_nodes = _assignment(generic, "W")
    ff_nodes = _assignment(generic, "FF")
    w_ok = bool(w_nodes) and _contains_name(w_nodes[-1], "F0") and _contains_name(w_nodes[-1], "radius") and _contains_name(w_nodes[-1], "inverse") and _contains_call(w_nodes[-1], "g.unary") and _contains_constant(w_nodes[-1], "sqrt")
    ff_ok = bool(ff_nodes) and _contains_name(ff_nodes[-1], "F0") and _contains_name(ff_nodes[-1], "radius") and _contains_name(ff_nodes[-1], "inverse")
    normalization_ok = (
        inverse_ok
        and w_ok
        and ff_ok
        and _contains_subscript(fields["E"], "actual", "phi")
        and _contains_subscript(fields["V"], "actual", "U")
        and _contains_subscript(radial["E_y"], "actual", "phi_y")
        and _contains_subscript(radial["E_y"], "actual", "phi")
        and _contains_subscript(radial["V_y"], "actual", "U_y")
        and _contains_call(fields["E"], "a.mul")
        and _contains_call(fields["V"], "a.scale")
    )
    _require(normalization_ok, gates, "c6_normalized_E_V_Ey_Vy_fields_not_bound")
    f0_ok = bool(_assignment(init, "F0")) and any(_contains_subscript(value, "self.core7", "F0") for value in _assignment(init, "F0"))
    p0_init = _attribute_assignment(init, "P0")
    p0_init_ok = bool(p0_init) and _contains_call(p0_init[-1], "self.pressure7")
    initial_ok = bool(_assignment(init, "phi0")) and bool(_assignment(init, "U0")) and bool(_assignment(init, "actualJ"))
    _require(f0_ok, gates, "c6_F0_not_taken_from_checked_core7")
    _require(p0_init_ok, gates, "c6_P0_not_taken_from_pressure7")
    _require(initial_ok, gates, "c6_common_initial_phi_U_or_six_history_rows_missing")
    micro_actual = _attribute_dict_keyword(init, "micro", "actual")
    micro_ok = bool(micro_actual) and _dotted(micro_actual) == "actual"
    _require(micro_ok, gates, "c6_micro_packet_not_bound_to_actual_fields")
    actual_core_ok = _constant(
        _attribute_dict_keyword(init, "definition", "actual_twenty_term_core_and_own_six_histories_retained"),
        True,
    )
    comparison_separate_ok = _constant(
        _attribute_dict_keyword(init, "definition", "known_direction_not_replaced_by_actual_moments"),
        True,
    )
    _require(actual_core_ok, gates, "c6_actual_twenty_term_core_binding_missing")
    _require(comparison_separate_ok, gates, "c6_comparison_profile_separation_missing")
    return {
        "passed": bool(actual_present and angular_ok and phi_ok and drive_ok and vy_ok and u_ok and phiy_ok and normalization_ok and f0_ok and p0_init_ok and initial_ok and micro_ok and actual_core_ok and comparison_separate_ok),
        "owner": "CurrentOriginalPrescribedBridgeC6Functions.micro",
        "actual_fields": {
            "Phi": "initial['phi'] * exp(integral(-chi * R * D_over_R / 2))",
            "U_y": "-chi * (Phi / barPhi) * (R*hydro + R*Pstar^2*pressure + R^2*F0^2*swirl)",
            "U": "initial['U'] - integral((-jacobian) * U_y) = initial['U'] + integral(jacobian * U_y)",
            "Phi_y": "-chi * R * D_over_R * Phi / 2",
        },
        "generic_fields": {
            "W": "F0 * sqrt(2*R) / Pstar",
            "FF": "F0^2 * R / Pstar^2",
            "E": "W * Phi",
            "V": "U / Pstar",
            "E_y": "W * (Phi/2 + Phi_y)",
            "V_y": "U_y / Pstar",
            "pressure": "P0 + p",
        },
        "own_moments": "actualJ = histories(actual Phi, actual U, initial7, phase, jacobian, actual_micro)",
        "P0": "self.P0 = self.pressure7()",
        "F0": "F0 = trim(self.core7['F0'])",
        "common_initial_data": "phi0 = core7 Phi at rho=4; U0 = core7 U0 + epsilon * core7 Psi at rho=4",
        "actual_micro_packet": "self.micro.actual, self.micro.moments, self.micro.fields, self.micro.radial, self.micro.generic",
        "actual_core_binding": "definition.actual_twenty_term_core_and_own_six_histories_retained=True",
        "comparison_profile_separation": "definition.known_direction_not_replaced_by_actual_moments=True",
        "normalization_factors_checked": normalization_ok,
        "physical_measure": "dy = hb * ds; radius R = Ra * exp(hb*s)",
        "numeric_values": False,
    }


def _pressure_witness(
    texts: Mapping[str, str],
    trees: Mapping[str, ast.AST | None],
    receipts: Mapping[str, bool],
    gates: list[dict[str, Any]],
) -> dict[str, Any]:
    """Bind the shared independent normalized ``P0`` source only.

    This is a proper positive sublemma: C2's ``op.P0`` is traced into the
    whole-Z pressure jet provider, while the C6 inlet source checks source,
    datum, input-hash, and ``normalized_jets.__func__`` identity before its
    ordinary Z-row projection.  It does not identify the Phi/U histories or
    F0, which remain separate gates below.
    """

    c2 = trees.get(_C2_PHASE)
    whole = trees.get(_BRIDGE)
    c6 = trees.get(_C6)
    foundation = trees.get(_FOUNDATION)
    if texts.get(_FOUNDATION) is None or texts.get(_FOUNDATION_CHECK) is None:
        _gate(gates, "pressure_foundation_records_not_supplied",
              records=sorted(_FOUNDATION_RECORDS))
        return {
            "passed": False,
            "quantity": "independent normalized P0 and ordinary Z rows",
            "status": "foundation source/check records are optional and were not supplied",
            "numeric_values": False,
        }
    c2_packet = _method(c2, "CurrentC2PhasePrimitives", "source_packet")
    c2_leading = _method(c2, "CurrentC2PhasePrimitives", "leading_source")
    c2_op = _assignment(c2_packet, "op")
    c2_owner_ok = any(_call(value, "self.outer.owner.owner") is not None for value in c2_op)
    c2_result_p0 = _dict_keyword(c2_packet, "result", "exact_common_P0_axial5")
    c2_p0_ok = (
        c2_owner_ok
        and _dotted(c2_result_p0) == "op.P0"
        and _contains_subscript_is_not(c2_packet, "packet", "exact_common_P0_axial5", "op.P0")
        and _contains_call(c2_leading, "self.outer.owner.source_query")
    )
    _require(c2_p0_ok, gates, "pressure_c2_op_P0_source_binding_missing")

    whole_owner = _method(whole, "OriginalWholeZBridgeSource", "__init__")
    whole_inputs = _method(whole, "OriginalWholeZBridgeSource", "inputs")
    whole_core_datum = _constructor_keyword(whole_owner, "core", "SimpleNamespace", "datum")
    whole_p0_ok = (
        whole_inputs is not None
        and _contains_call(whole_inputs, "self.admission.pressure.normalized_jets")
        and _contains_subscript(whole_inputs, "pressure", "normalized_pressure_coefficients")
        and _dict_keyword(whole_inputs, "result", "p0") is not None
        and _dotted(whole_core_datum) == "admission.pressure"
    )
    _require(whole_p0_ok, gates, "pressure_whole_z_normalized_jets_source_binding_missing")

    foundation_init = _method(foundation, "CurrentOriginalCoreInletAxial6GenericSource", "__init__")
    foundation_p0 = _method(foundation, "CurrentOriginalCoreInletAxial6GenericSource", "P0_functions")
    foundation_rows = _method(foundation, "CurrentOriginalCoreInletAxial6GenericSource", "P0_rows")
    foundation_common_ok = (
        _dotted(next(iter(_attribute_assignment(foundation_init, "common")), None)) == "self.core._oracle"
        and any(
            _call(value, "self.core.parameter_identity") is not None
            and isinstance(value, ast.Call)
            and len(value.args) == 1
            and not value.keywords
            and _dotted(value.args[0]) == "self.common"
            for value in _attribute_assignment(foundation_init, "parameter_identity")
        )
        and any(_call(value, "CoreSixMoments") is not None for value in _attribute_assignment(foundation_init, "moments"))
    )
    provider_values = _assignment(foundation_p0, "provider")
    datum_values = _assignment(foundation_p0, "datum")
    foundation_provider_ok = (
        any(_dotted(value) == "self.hierarchy.source.admission.pressure" for value in provider_values)
        and any(_dotted(value) == "self.common.core.datum" for value in datum_values)
        and _contains_not_equal(foundation_p0, "provider.source_sha", "datum.source_sha")
        and _contains_not_equal(foundation_p0, "provider.datum_sha", "datum.datum_sha")
        and _contains_not_equal(foundation_p0, "provider.input_hashes", "datum.input_hashes")
        and _contains_is_not(foundation_p0, "provider.normalized_jets.__func__", "datum.normalized_jets.__func__")
    )
    foundation_rows_ok = (
        _contains_subscript_call(foundation_rows, "self.hierarchy.source.inputs", "p0")
        and _contains_attribute(foundation_rows, "self.hierarchy.source")
    )
    _require(foundation_common_ok, gates, "pressure_core_parameter_identity_trace_missing")
    _require(foundation_provider_ok, gates, "pressure_core_provider_datum_identity_checks_missing")
    _require(foundation_rows_ok, gates, "pressure_core_P0_rows_source_projection_missing")

    c6_init = _method(c6, "CurrentOriginalPrescribedBridgeC6Functions", "__init__")
    c6_pressure7 = _method(c6, "CurrentOriginalPrescribedBridgeC6Functions", "pressure7")
    c6_rows = _method(c6, "CurrentOriginalPrescribedBridgeC6Functions", "P0_rows")
    c6_inlet_ok = (
        any(_dotted(value) == "inlet_owner" for value in _attribute_assignment(c6_init, "inlet"))
        and any(_call(value, "self.pressure7") is not None for value in _attribute_assignment(c6_init, "P0"))
        and _contains_attribute(c6_pressure7, "self.inlet.P0")
    )
    c6_rows_ok = (
        any(_dotted(value) == "self.inlet.hierarchy.source.admission.pressure" for value in _assignment(c6_rows, "provider"))
        and _contains_call(c6_rows, "provider.normalized_jets")
        and _contains_subscript(c6_rows, "packet", "all_14_pressure_atoms_included")
        and _contains_subscript(c6_rows, "packet", "normalized_pressure_coefficients")
    )
    _require(c6_inlet_ok, gates, "pressure_c6_inlet_P0_projection_binding_missing")
    _require(c6_rows_ok, gates, "pressure_c6_normalized_jets_rows_binding_missing")
    receipt_ok = receipts.get(_FOUNDATION_CHECK) is True
    _require(receipt_ok, gates, "pressure_foundation_receipt_not_accepted")
    passed = bool(c2_p0_ok and whole_p0_ok and foundation_common_ok and foundation_provider_ok and foundation_rows_ok and c6_inlet_ok and c6_rows_ok and receipt_ok)
    return {
        "passed": passed,
        "quantity": "independent normalized P0 and ordinary Z rows",
        "c2_binding": "source_packet exact_common_P0_axial5 is op.P0; op=self.outer.owner.owner(ends)",
        "whole_z_binding": "owner.inputs p0 = admission.pressure.normalized_jets(..., 6)",
        "core_branch_binding": "core.parameter_identity(common) and CoreSixMoments(common) are retained in the checked inlet source",
        "provider_identity": "provider.source_sha, datum_sha, input_hashes and normalized_jets.__func__ are compared",
        "c6_binding": "C6 pressure7 extends inlet.P0 and P0_rows reads inlet.hierarchy.source.inputs(Z)['p0']",
        "ordinary_rows": "C6 order-7 and core-inlet order-6 projections remain source rows, not range values",
        "numeric_values": False,
        "receipt": receipt_ok,
    }


def _core_projection_witness(
    texts: Mapping[str, str],
    trees: Mapping[str, ast.AST | None],
    receipts: Mapping[str, bool],
    gates: list[dict[str, Any]],
) -> dict[str, Any]:
    """Admit the common analytic core projection as a bounded sublemma.

    The accepted common fixed-point receipt closes the finite rows and the
    analytic tail to one source-defined solution.  The AST chain below then
    binds that receipt to its exact operator/rebuild constructors and checks
    that the C6 inlet consumes the checked N1 owner.  This is deliberately a
    core projection statement: it does not identify the first-micro Phi/U
    histories or F0 values owned by the two downstream providers.
    """

    required = (
        _CORE_REBUILD,
        _CORE_REBUILD_CHECK,
        _CORE_COMMON,
        _CORE_COMMON_CHECK,
        _CORE_OPERATOR,
        _CORE_OPERATOR_CHECK,
        _N1_ANALYTIC,
        _N1_ANALYTIC_CHECK,
        _FOUNDATION,
        _FOUNDATION_CHECK,
        _C6,
        _C6_CHECK,
    )
    missing = [name for name in required if texts.get(name) is None]
    if missing:
        _gate(gates, "common_core_projection_records_not_supplied", records=sorted(missing))
        return {
            "passed": False,
            "full_micro_source_identity_proved": False,
            "status": "bounded common-core source records are optional and were not supplied",
            "records": sorted(required),
            "numeric_values": False,
        }

    common_init = _method(trees.get(_CORE_COMMON), "CurrentCoreCommonFixedPoint", "__init__")
    common_eval = _method(trees.get(_CORE_COMMON), "CurrentCoreCommonFixedPoint", "evaluate")
    common_consumer = _function(trees.get(_CORE_COMMON), "consumer_source_bindings")
    common_operator = _attribute_assignment(common_init, "operator")
    common_default_ok = any(
        _expr_equal(value, "operator if operator is not None else CurrentCoreNonlinearOperator()")
        for value in common_operator
    )
    common_chain_ok = (
        any(_dotted(value) == "self.operator.source" for value in _attribute_assignment(common_init, "source"))
        and any(_dotted(value) == "self.source.core" for value in _attribute_assignment(common_init, "core"))
        and any(_dotted(value) == "self.source.rebuild" for value in _attribute_assignment(common_init, "rebuild"))
        and any(_dotted(value) == "self.core.ctx" for value in _attribute_assignment(common_init, "ctx"))
        and _contains_attribute(common_init, "self.operator.acceptance_loaded")
    )
    common_consumer_ok = common_consumer is not None and all(
        _contains_constant(common_consumer, value)
        for value in (
            "advance_scaled_one(c, fixed, rows, n, c.mpf(Z), self.core.delta, self.core.epsilon)",
        )
    ) and all(
        _constant(_call_keyword(common_consumer, "dict", key), True)
        for key in (
            "actual_all_order_map_and_integrals_AST_bound",
            "original_factorial_model_tail_consumer_AST_bound",
            "original_model_plus_nonlinear_norm_tail_consumers_AST_bound",
            "rooted_nonlinear_tail_uses_same_actual_map_correction",
            "same_pressure_primitive_coefficient_tail_consumer_AST_bound",
        )
    )
    common_eval_ok = (
        common_eval is not None
        and _contains_call(common_eval, "self.rebuild.rebuild")
        and _contains_call(common_eval, "self.rebuild.profile")
    )
    _require(common_default_ok, gates, "common_fixed_point_operator_constructor_not_exactly_bound")
    _require(common_chain_ok, gates, "common_fixed_point_source_core_rebuild_chain_not_bound")
    _require(common_consumer_ok, gates, "common_fixed_point_old_rebuild_consumer_AST_not_bound")
    _require(common_eval_ok, gates, "common_fixed_point_evaluation_source_not_bound")

    rebuild_init = _method(trees.get(_CORE_REBUILD), "CompliantCoreCoefficientRebuild", "__init__")
    rebuild_method = _method(trees.get(_CORE_REBUILD), "CompliantCoreCoefficientRebuild", "rebuild")
    rebuild_calls = [
        child
        for child in ast.walk(rebuild_method) if isinstance(child, ast.Call)
        and _dotted(child.func) == "advance_scaled_one"
        and ast.unparse(child) == "advance_scaled_one(c, fixed, rows, n, c.mpf(Z), self.core.delta, self.core.epsilon)"
    ] if rebuild_method is not None else []
    rebuild_ast_ok = (
        rebuild_init is not None
        and _contains_call(rebuild_init, "CompliantCorePhysicalField")
        and len(rebuild_calls) == 1
        and _contains_call(rebuild_method, "self.seed")
    )
    _require(rebuild_ast_ok, gates, "core_coefficient_rebuild_exact_advance_scaled_one_call_missing")

    operator_init = _method(trees.get(_CORE_OPERATOR), "CurrentCoreNonlinearOperator", "__init__")
    operator_default_ok = any(
        _expr_equal(value, "source if source is not None else CurrentCoreScaledSwirlSource()")
        for value in _attribute_assignment(operator_init, "source")
    )
    operator_chain_ok = (
        operator_init is not None
        and _contains_attribute(operator_init, "self.source.acceptance_loaded")
        and any(_dotted(value) == "self.source.ctx" for value in _attribute_assignment(operator_init, "ctx"))
    )
    _require(operator_default_ok, gates, "core_nonlinear_operator_source_constructor_not_exactly_bound")
    _require(operator_chain_ok, gates, "core_nonlinear_operator_acceptance_source_chain_not_bound")

    n1_init = _method(trees.get(_N1_ANALYTIC), "CurrentOriginalN1AnalyticCoreLowerPressureC3Functions", "__init__")
    n1_derivative_init = _method(trees.get(_N1_ANALYTIC), "HigherCoreZDerivative", "__init__")
    n1_guard_ok = (
        n1_init is not None
        and _contains_type_identity(n1_init, "lower_owner", "lower.CurrentOriginalN1LowerCollarFunctions")
        and _contains_call(n1_init, "lower_owner.assert_graph")
        and any(_dotted(value) == "lower_owner.core" for value in _attribute_assignment(n1_init, "core"))
    )
    n1_common_ok = (
        n1_derivative_init is not None
        and any(_dotted(value) == "lower_owner.core" for value in _attribute_assignment(n1_derivative_init, "core"))
        and _constant(_attribute_dict_keyword(n1_init, "definition", "same_unique_analytic_nonlinear_core_components"), True)
        and _constant(_attribute_dict_keyword(n1_init, "definition", "actual_twenty_term_source_and_Green_inverse_not_changed"), True)
    )
    _require(n1_guard_ok, gates, "n1_analytic_core_lower_owner_guard_or_core_binding_missing")
    _require(n1_common_ok, gates, "n1_analytic_core_common_component_definition_missing")

    foundation_init = _method(trees.get(_FOUNDATION), "CurrentOriginalCoreInletAxial6GenericSource", "__init__")
    foundation_n1_ok = (
        foundation_init is not None
        and _contains_type_identity(
            foundation_init,
            "core_lower_owner",
            "previous.CurrentOriginalN1AnalyticCoreLowerPressureC3Functions",
        )
        and _contains_call(foundation_init, "core_lower_owner.assert_graph")
        and any(_dotted(value) == "core_lower_owner.core" for value in _attribute_assignment(foundation_init, "core"))
        and any(_dotted(value) == "core_lower_owner" for value in _attribute_assignment(foundation_init, "previous"))
        and _contains_attribute(foundation_init, "self.core._oracle")
        and _contains_call(foundation_init, "self.core.parameter_identity")
        and _contains_call(foundation_init, "CoreSixMoments")
    )
    c6_init = _method(trees.get(_C6), "CurrentOriginalPrescribedBridgeC6Functions", "__init__")
    c6_inlet_ok = (
        c6_init is not None
        and _contains_type_identity(c6_init, "inlet_owner", "inlet.CurrentOriginalCoreInletAxial6GenericSource")
        and _contains_call(c6_init, "inlet_owner.assert_graph")
        and any(_dotted(value) == "inlet_owner" for value in _attribute_assignment(c6_init, "inlet"))
        and any(_dotted(value) == "inlet_owner.core" for value in _attribute_assignment(c6_init, "core"))
    )
    _require(foundation_n1_ok, gates, "core_inlet_does_not_consume_checked_n1_source")
    _require(c6_inlet_ok, gates, "c6_does_not_consume_checked_core_inlet_source")

    common_receipt_ok = receipts.get(_CORE_COMMON_CHECK) is True
    operator_receipt_ok = receipts.get(_CORE_OPERATOR_CHECK) is True
    n1_receipt = _json(texts.get(_N1_ANALYTIC_CHECK), _N1_ANALYTIC_CHECK, gates)
    parent_receipt = _json(texts.get(_PARENT_REPORT), _PARENT_REPORT, gates)
    n1_input_hashes = n1_receipt.get("input_hashes") if isinstance(n1_receipt, Mapping) else None
    parent_input_hashes = parent_receipt.get("input_hashes") if isinstance(parent_receipt, Mapping) else None
    # The C3 receipt cannot self-pin its own check JSON.  Its input_hashes
    # therefore bind the actual N1 source and lower-collar source/check, while
    # the authenticated immutable parent binds all six paths, including the
    # C3 check JSON itself.
    n1_provider_input_pins_ok = (
        isinstance(n1_input_hashes, Mapping)
        and all(
            n1_input_hashes.get(name) == digest
            for name, digest in _N1_ANALYTIC_EXPECTED_INPUT_PINS.items()
            if name != _N1_ANALYTIC_CHECK
        )
    )
    parent_input_pins_ok = (
        isinstance(parent_input_hashes, Mapping)
        and all(parent_input_hashes.get(name) == digest for name, digest in _N1_ANALYTIC_EXPECTED_INPUT_PINS.items())
    )
    n1_receipt_nested_ok = receipts.get(_N1_ANALYTIC_CHECK) is True
    n1_receipt_ok = bool(n1_receipt_nested_ok and n1_provider_input_pins_ok and parent_input_pins_ok)
    foundation_receipt_ok = receipts.get(_FOUNDATION_CHECK) is True
    c6_receipt_ok = receipts.get(_C6_CHECK) is True
    rebuild_receipt = _json(texts.get(_CORE_REBUILD_CHECK), _CORE_REBUILD_CHECK, gates)
    rebuild_receipt_ok = (
        isinstance(rebuild_receipt, Mapping)
        and rebuild_receipt.get("all_passed") is True
        and rebuild_receipt.get("old_finite_coefficient_rows_read") is False
        and rebuild_receipt.get("same_selected_Cstar_and_compliant_pressure_source") is True
        and rebuild_receipt.get("positive_swirl_source_not_dropped") is True
        and rebuild_receipt.get("infinite_radial_tails_bound") is True
    )
    _require(common_receipt_ok, gates, "accepted_common_fixed_point_projection_receipt_missing")
    _require(operator_receipt_ok, gates, "accepted_core_nonlinear_operator_receipt_missing")
    _require(rebuild_receipt_ok, gates, "accepted_core_rebuild_projection_receipt_missing")
    _require(n1_provider_input_pins_ok, gates, "n1_c3_receipt_input_hashes_not_exact_parent_pins")
    _require(parent_input_pins_ok, gates, "n1_c3_parent_input_hashes_not_exact_authenticated_pins")
    _require(n1_receipt_ok, gates, "accepted_n1_analytic_core_receipt_missing")
    _require(foundation_receipt_ok, gates, "accepted_core_inlet_projection_receipt_missing")
    _require(c6_receipt_ok, gates, "accepted_c6_inlet_projection_receipt_missing")

    passed = bool(
        common_default_ok and common_chain_ok and common_consumer_ok and common_eval_ok
        and rebuild_ast_ok and operator_default_ok and operator_chain_ok
        and n1_guard_ok and n1_common_ok and foundation_n1_ok and c6_inlet_ok
        and common_receipt_ok and operator_receipt_ok and rebuild_receipt_ok
        and n1_receipt_ok and foundation_receipt_ok and c6_receipt_ok
    )
    return {
        "passed": passed,
        "full_micro_source_identity_proved": False,
        "status": "bounded common analytic core projection admitted" if passed else "bounded core projection remains gated",
        "accepted_projection_flags": {
            "common_fixed_point_receipt": common_receipt_ok,
            "operator_receipt": operator_receipt_ok,
            "coefficient_rebuild_receipt": rebuild_receipt_ok,
            "n1_analytic_core_receipt": n1_receipt_ok,
            "n1_c3_nested_theorem_flags": n1_receipt_nested_ok,
            "n1_c3_receipt_input_hashes": n1_provider_input_pins_ok,
            "n1_c3_parent_authenticated_input_hashes": parent_input_pins_ok,
            "core_inlet_receipt": foundation_receipt_ok,
            "c6_receipt": c6_receipt_ok,
        },
        "premise_checks": {
            "c3_nested_analytic_source_passed": _saved_path(n1_receipt, ("checks", "analytic_source", "passed")) is True,
            "c3_nested_same_unique_analytic_twenty_term_core_and_radial_definitions": _saved_path(n1_receipt, ("checks", "analytic_source", "same_unique_analytic_twenty_term_core_and_radial_definitions")) is True,
            "c3_nested_lower_source_passed": _saved_path(n1_receipt, ("checks", "lower_source", "passed")) is True,
            "c3_receipt_input_hashes": {
                name: isinstance(n1_input_hashes, Mapping) and n1_input_hashes.get(name) == digest
                for name, digest in _N1_ANALYTIC_EXPECTED_INPUT_PINS.items()
                if name != _N1_ANALYTIC_CHECK
            },
            "full_parent_authenticated_input_hashes": {
                name: isinstance(parent_input_hashes, Mapping) and parent_input_hashes.get(name) == digest
                for name, digest in _N1_ANALYTIC_EXPECTED_INPUT_PINS.items()
            },
        },
        "common_source_binding": {
            "operator_default": "CurrentCoreNonlinearOperator()",
            "source_core_rebuild_ctx": "self.operator.source -> source.core/rebuild -> core.ctx",
            "finite_tail_consumer": "consumer_source_bindings exact advance_scaled_one call and admitted tail consumers",
            "common_runtime_evaluation": "evaluate -> rebuild.rebuild and rebuild.profile",
        },
        "coefficient_rebuild_binding": {
            "exact_call": "advance_scaled_one(c, fixed, rows, n, c.mpf(Z), self.core.delta, self.core.epsilon)",
            "call_count": len(rebuild_calls),
            "old_finite_rows_read": rebuild_receipt.get("old_finite_coefficient_rows_read") is False if isinstance(rebuild_receipt, Mapping) else False,
        },
        "n1_common_core_consumption": {
            "lower_owner_guard": n1_guard_ok,
            "same_unique_analytic_core_component": n1_common_ok,
            "core_binding": "self.core = lower_owner.core and HigherCoreZDerivative.core = lower_owner.core",
            "accepted_nested_theorem_and_input_pins": n1_receipt_ok,
            "authenticated_input_paths": sorted(_N1_ANALYTIC_EXPECTED_INPUT_PINS),
        },
        "c6_inlet_n1_consumption": {
            "core_inlet_type_guard": foundation_n1_ok,
            "c6_inlet_type_guard": c6_inlet_ok,
            "core_inlet_common_oracle": "self.common = self.core._oracle; parameter_identity(common); CoreSixMoments(common)",
        },
        "analytic_scope": {
            "finite_rows_and_tails_same_unique_solution": common_receipt_ok,
            "ordinary_factorials_and_radial_units_preserved": common_receipt_ok,
            "source_functions_not_receipt_ids_only": common_receipt_ok,
        },
        "numeric_values": False,
    }


def _root_witness(texts: Mapping[str, str], trees: Mapping[str, ast.AST | None], gates: list[dict[str, Any]]) -> dict[str, Any]:
    c5 = _function(trees.get(_C5), "generic_roots")
    generic = _function(trees.get(_GENERIC), "recover_inputs")
    restore = _function(trees.get(_RESTORE), "quotient_cell")
    loop = _function(trees.get(_C5), "loop_functions")
    c_node = _assignment(c5, "C")
    b_node = _assignment(c5, "B")
    aa_bb_nodes = _tuple_assignment_values(c5, ("aa", "bb"))
    aa_node = [aa_bb_nodes[0]] if aa_bb_nodes is not None else []
    bb_node = [aa_bb_nodes[1]] if aa_bb_nodes is not None else []
    t0_node = _assignment(c5, "t0")
    p2_node = _assignment(c5, "p2")
    delta_node = _returned_dict_keyword(c5, "Delta")
    aa_bb_ok = (
        aa_bb_nodes is not None
        and all(_call(value, "divide") is not None for value in aa_bb_nodes)
        and len(aa_bb_nodes[0].args) == 1
        and len(aa_bb_nodes[1].args) == 1
        and _contains_name(aa_bb_nodes[0].args[0], "C")
        and _contains_name(aa_bb_nodes[1].args[0], "B")
    )
    delta_call = _call(delta_node, "a.sub")
    delta_fixed = (
        _call(delta_call.args[1], "a.fixed")
        if delta_call is not None and len(delta_call.args) == 2
        else None
    )
    delta_ok = (
        delta_call is not None
        and len(delta_call.args) == 2
        and _contains_name(delta_call.args[0], "kap")
        and delta_fixed is not None
        and len(delta_fixed.args) == 1
        and _constant(delta_fixed.args[0], 2)
    )
    roots_ok = (
        bool(c_node) and bool(b_node) and aa_bb_ok
        and bool(t0_node) and bool(p2_node) and delta_ok
        and _contains_subscript(c_node[-1], "radial", "E_y")
        and _contains_subscript(b_node[-1], "radial", "V_y")
        and _contains_name(t0_node[-1], "B")
        and _contains_name(t0_node[-1], "C")
        and _subscript(p2_node[-1], "generic", "p2")
    )
    # Exact source equations are reported independently of the structural
    # booleans so a reviewer can see the signed factors and ownership.
    _require(roots_ok, gates, "prescribed_C6_generic_root_recipe_not_structurally_bound")
    generic_ok = (
        generic is not None
        and all(bool(_assignment(generic, name)) for name in ("E", "Ey", "V", "Vy", "C", "B", "R", "It", "Iz"))
        and _contains_subscript(generic, "packet", "original_P0_normalized_axial5")
        and _contains_subscript(generic, "packet", "geometry")
    )
    _require(generic_ok, gates, "actual_Rm_generic_inputs_not_structurally_bound")
    p2_ok = restore is not None and bool(_assignment(restore, "p2")) and _contains_attribute(restore, "proxy.source_radius")
    _require(p2_ok, gates, "full_signed_p2_radius_factor_not_structurally_bound")
    loop_distinct = loop is not None and _contains_subscript(loop, "roots", "p2") and _contains_subscript(loop, "roots", "a")
    _require(loop_distinct, gates, "prescribed_loop_root_inputs_missing")
    return {
        "passed": bool(roots_ok and generic_ok and p2_ok and loop_distinct),
        "normalized_defining_equations": {
            "C": "E - 2*E_y",
            "B": "2*V_y = 2*U_y/Pstar (actual source B; distinct from loop B=E*M/2)",
            "a": "C/E = (E - 2*E_y)/E",
            "b": "B/E = 2*V_y/E",
            "t0": "-B/C = -b/a",
            "kappa": "a + b^2/a",
            "Delta": "kappa - 2",
            "p2": "R * (axial_linear + Pstar*axial_quadratic) / E, with full signed sectors retained",
        },
        "root_paths": {
            "E": "common_velocity_E_axial5",
            "V": "common_velocity_V_axial5",
            "a": "actual_correlated_a_axial5",
            "b": "actual_nonzero_b_axial5",
            "t0": "actual_t0_axial5",
            "p2": "full_signed_p2_axial4",
            "Delta": "actual_Delta_axial5",
        },
        "source_p2_factors": {
            "source_radius": "proxy.source_radius",
            "Pstar": "one full Pstar multiplier on axial_quadratic",
            "signed_sectors": "axial_linear + Pstar*axial_quadratic",
        },
        "loop_B_distinction": "loop_functions B_over_Pstar = E*M/2; it is not actual generic source B = 2*V_y",
    }


def _report_witness(reports: Mapping[str, Any], gates: list[dict[str, Any]]) -> dict[str, Any]:
    latest = reports.get(_LATEST_REPORT)
    parent = reports.get(_PARENT_REPORT)
    latest_ok = isinstance(latest, Mapping) and latest.get("parent_report") == _PARENT_REPORT and latest.get("parent_report_sha256") == RECORD_HASHES[_PARENT_REPORT]
    _require(latest_ok, gates, "latest_source_report_parent_pin_missing")
    parent_ok = isinstance(parent, Mapping)
    _require(parent_ok, gates, "immutable_parent_report_missing")
    family_ok = True
    for name, report in (("latest", latest), ("parent", parent)):
        if not isinstance(report, Mapping):
            continue
        family = report.get("source_family")
        if family is not None and family != SOURCE_FAMILY:
            family_ok = False
            _gate(gates, "report_source_family_changed", report=name)
    return {"latest_parent_pin": latest_ok, "parent_present": parent_ok, "report_source_family": family_ok}


def _source_witness(source_records: Mapping[str, Any], gates: list[dict[str, Any]]) -> dict[str, Any]:
    texts: dict[str, str] = {}
    metadata: dict[str, Any] = {}
    for name in RECORD_HASHES:
        text, meta = _record(
            source_records,
            name,
            gates,
            optional=name in ARCHIVE_RECORDS or name in _FOUNDATION_RECORDS or name in _CORE_PROJECTION_RECORDS,
        )
        metadata[name] = meta
        if text is not None:
            texts[name] = text

    missing = sorted(_REQUIRED_RECORDS - set(texts))
    if missing:
        _gate(gates, "required_source_records_missing", records=missing)
        return {
            "records_present": sorted(texts),
            "records_metadata": metadata,
            "missing_required": missing,
            "receipts": {},
            "reports": {},
            "c2": {},
            "whole_z": {},
            "correlated": {},
            "c6": {},
            "roots": {},
            "actual_math_coverage": {},
        }

    trees = {name: _ast(text, name, gates) for name, text in texts.items() if name.endswith(".py")}
    receipts = {
        _C2_PHASE_CHECK: _receipt(texts.get(_C2_PHASE_CHECK), _C2_PHASE_CHECK, (
            ("current_original_17_chart_source_Z2_rows_installed",),
            ("current_original_same_inverse_C2_phase_primitive_functions_installed",),
            ("upstream_six_chart_C2_provider_installed",),
            ("same_independent_P0_and_source_basis_ledger_checked",),
        ), gates)[0],
        _C2_DENSITY_CHECK: _receipt(texts.get(_C2_DENSITY_CHECK), _C2_DENSITY_CHECK, (
            ("current_original_17_chart_C2_signed_density_functions_installed",),
            ("current_original_17_chart_C2_own_rate_transport_and_five_target_functions_installed",),
            ("actual_five_target_ZZ_function_representations_installed",),
            ("same_actual_terminal_amplitude_P0_source_and_global_phase_retained",),
        ), gates)[0],
        _BRIDGE_CHECK: _receipt(texts.get(_BRIDGE_CHECK), _BRIDGE_CHECK, (
            ("current_original_whole_Z_live_leading_bridge_source_functions_installed",),
            ("symbolic_source", "original_all_six_own_ODEs_and_core_inlet_identities"),
        ), gates)[0],
        _CORRELATED_CHECK: _receipt(texts.get(_CORRELATED_CHECK), _CORRELATED_CHECK, (
            ("current_original_whole_Z_correlated_bridge_C1_functions_and_Rc_retransport_installed",),
        ), gates)[0],
        _C6_CHECK: _receipt(texts.get(_C6_CHECK), _C6_CHECK, (
            ("current_original_true_comparison_bridge_C7_functions_installed",),
            ("current_original_true_prescribed_micro_macro_C6_functions_installed",),
            ("current_original_true_prescribed_bridge_Q_inertial_p2_C5_functions_installed",),
            ("checks", "actual_source", "passed"),
            ("checks", "actual_source", "original_P0_C6_preserved_and_executable_Z7_projection"),
            ("checks", "original_direction", "physical_swirl_F0_squared_dressing_and_external_scale_once"),
            ("checks", "signed_transport", "full_axial_minus_sign_and_pressure_swirl_physical_scales_once"),
            ("checks", "own_six_histories", "all_nonzero_incoming_histories_and_sixth_derivatives_retained"),
        ), gates)[0],
        _FOUNDATION_CHECK: False,
        _CORE_COMMON_CHECK: _receipt(texts.get(_CORE_COMMON_CHECK), _CORE_COMMON_CHECK, (
            ("actual_Banach_and_common_coefficient_proof", "actual_source_operator_defined_on_same_complete_Xh_pair"),
            ("actual_Banach_and_common_coefficient_proof", "unique_actual_analytic_fixed_point_for_each_fixed_implicit_source"),
            ("actual_Banach_and_common_coefficient_proof", "all_original_finite_rows_enclose_coefficients_of_this_unique_solution"),
            ("actual_Banach_and_common_coefficient_proof", "no_interval_overlap_or_fitted_source_equality_argument"),
            ("actual_Banach_and_common_coefficient_proof", "same_source_analytic_functions_not_just_shared_receipt_ids"),
            ("actual_Banach_and_common_coefficient_proof", "passed"),
            ("same_solution_all_order_tail_proof", "same_actual_fixed_point_coefficient_norm_used"),
            ("same_solution_all_order_tail_proof", "finite_and_infinite_parts_of_same_analytic_solution"),
            ("same_solution_all_order_tail_proof", "ordinary_derivative_factorials_and_radial_units_preserved"),
            ("same_solution_all_order_tail_proof", "passed"),
            ("actual_all_order_pressure_source_proof", "original_all_order_normalized_jets_enclose_this_same_function"),
            ("actual_all_order_pressure_source_proof", "passed"),
            ("actual_twenty_majorant_rows_replayed",),
            ("finite_rows_and_tails_bound_to_same_nonlinear_fixed_point",),
            ("current_core_actual_operator_Banach_fixed_point_admitted",),
            ("current_core_actual_operator_contraction_and_tail_certified",),
        ), gates)[0] if texts.get(_CORE_COMMON_CHECK) is not None else False,
        _CORE_OPERATOR_CHECK: _receipt(texts.get(_CORE_OPERATOR_CHECK), _CORE_OPERATOR_CHECK, (
            ("source_bindings_and_incidence", "passed"),
            ("actual_pressure_primitive_and_its_Z_derivative_use_same_source",),
            ("formal_operator_does_not_replace_exact_source_with_Cauchy_box",),
            ("current_core_actual_twenty_term_operator_and_radial_extraction_certified",),
            ("current_core_actual_twenty_term_norm_incidence_certified",),
        ), gates)[0] if texts.get(_CORE_OPERATOR_CHECK) is not None else False,
        _CORE_REBUILD_CHECK: _receipt(texts.get(_CORE_REBUILD_CHECK), _CORE_REBUILD_CHECK, (
            ("independent_original_equations", "passed"),
            ("independent_infinite_tail_checks", "passed"),
            ("same_selected_Cstar_and_compliant_pressure_source",),
            ("positive_swirl_source_not_dropped",),
            ("infinite_radial_tails_bound",),
        ), gates)[0] if texts.get(_CORE_REBUILD_CHECK) is not None else False,
        _N1_ANALYTIC_CHECK: _receipt(texts.get(_N1_ANALYTIC_CHECK), _N1_ANALYTIC_CHECK, (
            ("current_original_n1_actual_analytic_core_required_C5_functions_installed",),
            ("current_original_n1_actual_lower_U_mean_C5_Q_C4_functions_installed",),
            ("current_original_n1_true_lower_pressure_C3_functions_installed",),
            ("checks", "analytic_source", "passed"),
            ("checks", "analytic_source", "same_unique_analytic_twenty_term_core_and_radial_definitions"),
            ("checks", "lower_source", "passed"),
        ), gates)[0] if texts.get(_N1_ANALYTIC_CHECK) is not None else False,
        _C5_CHECK: _receipt(texts.get(_C5_CHECK), _C5_CHECK, (
            ("current_original_true_prescribed_ten_chart_signed_loop_C5_functions_installed",),
            ("current_original_true_prescribed_ten_chart_finite_N_density_C5_functions_installed",),
            ("checks", "actual_source", "genuine_C6_inputs_and_C5_generic_roots"),
            ("checks", "actual_source", "full_signed_p2_preserved_by_identity"),
            ("checks", "actual_source", "same_original_N_eta_dstar_and_global_phase"),
        ), gates)[0],
        _GENERIC_CHECK: _receipt(texts.get(_GENERIC_CHECK), _GENERIC_CHECK, (
            ("original_actual_two_frame_Rm_full_signed_generic_inputs_installed",),
        ), gates)[0],
        _RESTORE_CHECK: _receipt(texts.get(_RESTORE_CHECK), _RESTORE_CHECK, (
            ("current_original_Rsh_Rm_full_reference_restoration_source_and_finite_N_local_drivers_installed",),
        ), gates)[0],
    }
    if texts.get(_FOUNDATION_CHECK) is not None:
        receipts[_FOUNDATION_CHECK] = _receipt(texts.get(_FOUNDATION_CHECK), _FOUNDATION_CHECK, (
            ("current_original_true_core_inlet_C6_generic_source_functions_installed",),
            ("current_original_true_core_inlet_axial6_numeric_source_installed",),
            ("checks", "source_rows", "passed"),
            ("checks", "source_rows", "original_P0_C3_preserved_and_higher_rows_are_explicit_function_derivatives"),
            ("checks", "source_rows", "executable_order6_pressure_projection_bound"),
            ("checks", "actual_numeric_source", "same_live_formal_bases_and_directed_ledger"),
        ), gates)[0]
    reports = {}
    for name in REPORT_RECORDS:
        reports[name] = _json(texts.get(name), name, gates)
    report_witness = _report_witness(reports, gates)
    c2 = _c2_witness(trees.get(_C2_PHASE), gates)
    whole_z = _whole_z_witness(trees.get(_BRIDGE), gates)
    correlated = _correlated_witness(trees.get(_CORRELATED), gates)
    c6 = _c6_actual_witness(trees.get(_C6), gates)
    pressure = _pressure_witness(texts, trees, receipts, gates)
    core_projection = _core_projection_witness(texts, trees, receipts, gates)
    roots = _root_witness(texts, trees, gates)
    for name, passed in receipts.items():
        if not passed and name != _FOUNDATION_CHECK and name not in _CORE_PROJECTION_RECORDS:
            _gate(gates, "source_receipt_not_accepted", record=name)
    conditional = bool(c2.get("passed") and whole_z.get("passed") and correlated.get("passed") and c6.get("passed") and roots.get("passed") and all(receipts.values()) and report_witness.get("latest_parent_pin") and report_witness.get("parent_present"))
    # The pressure subledger is independently bound above.  The remaining
    # actual fields still lack a signed per-quantity cross-owner identity.
    _gate(gates, "full_micro_evolution_initial_data_parameter_moments_F0_correspondence_unbound",
          required="same full first-micro evolution equations, initial Phi/U data, live parameters, all six incoming histories, own moments, F0, and source basis")
    if not pressure.get("passed"):
        _gate(gates, "cross_owner_common_pressure_P0_identity_unbound",
              required="checked core-inlet source record and provider/datum normalized_jets identity")
    _gate(gates, "cross_owner_normalized_root_value_identity_unbound",
          required="same C2 source-query actual fields and C6 micro fields before E/V/Ey/Vy normalization")
    return {
        "records_present": sorted(texts),
        "records_metadata": metadata,
        "missing_required": missing,
        "receipts": receipts,
        "reports": report_witness,
        "c2": c2,
        "whole_z": whole_z,
        "correlated": correlated,
        "c6": c6,
        "pressure": pressure,
        "common_core_projection_sublemma": core_projection,
        "roots": roots,
        "conditional_formula_witness": conditional,
        "actual_math_coverage": {
            "whole_Z_correlated_source_query": correlated.get("passed") is True,
            "whole_Z_actual_Phi_U_own_moments": whole_z.get("passed") is True,
            "prescribed_C6_actual_micro_Phi_U_own_moments": c6.get("passed") is True,
            "common_initial_data": False,
            "common_pressure_P0": pressure.get("passed") is True,
            "common_analytic_core_projection": core_projection.get("passed") is True,
            "common_amplitude_F0": False,
            "normalized_E_V_Ey_Vy": roots.get("passed") is True and c6.get("normalization_factors_checked") is True,
            "signed_p2_with_radius_and_Pstar": roots.get("passed") is True,
            "actual_same_value_identity": False,
        },
    }


def _row(witness: Mapping[str, Any], gates: list[dict[str, Any]]) -> dict[str, Any]:
    return {
        "chart": "first_micro",
        "domain": list(FIRST_MICRO_DOMAIN),
        "status": "unresolved",
        "exact_identity": None,
        "actual_C2_provider": "WholeZCorrelatedBridgeFunctions.source_query -> CurrentC2PhasePrimitives.source_packet",
        "prescribed_C6_provider": "CurrentOriginalPrescribedBridgeC6Functions.micro",
        "actual_fields": ["Phi", "U", "own_moments", "P0", "F0"],
        "normalized_fields": ["E", "V", "E_y", "V_y"],
        "root_paths": witness.get("roots", {}).get("root_paths", {}),
        "conditional_equations": witness.get("roots", {}).get("normalized_defining_equations", {}),
        "positive_sublemmas": {
            "common_pressure_P0": witness.get("pressure", {}).get("passed") is True,
            "common_analytic_core_projection": witness.get("common_core_projection_sublemma", {}).get("passed") is True,
        },
        "same_actual_source_identity": False,
        "unresolved_gates": [dict(gate) for gate in gates],
        "source_only": True,
        "numeric_values_computed": False,
        "operator_calls": 0,
    }


def derive(nodes: Sequence[Any], source_records: Mapping[str, Any], *, expected_nodes_digest: str) -> dict[str, Any]:
    """Return a fail-closed first-micro source identity packet.

    ``source_records`` is a mapping keyed by the basename in ``RECORD_HASHES``
    (the corresponding ``experiments/root_st073/`` prefix is also accepted).
    Plain records carry ``source_text`` or UTF-8 ``raw_bytes``.  Gzip records
    carry raw ``raw_bytes`` or ``raw_bytes_base64`` and must also carry
    ``path`` and ``sha256``.  No
    record path is read from disk by this function.
    """

    if expected_nodes_digest != EXPECTED_NODES_DIGEST:
        raise ValueError("unexpected source graph digest assertion")
    try:
        node_count = len(nodes)
        actual_digest = canonical_digest(nodes)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError("nodes are not canonical JSON graph data") from exc
    if node_count != EXPECTED_NODE_COUNT or actual_digest != EXPECTED_NODES_DIGEST:
        raise ValueError("nodes canonical digest or count does not match expected source graph")

    gates: list[dict[str, Any]] = []
    witness = _source_witness(source_records, gates)
    row = _row(witness, gates)
    return {
        "schema": SCHEMA,
        "source_only": True,
        "passed": False,
        "full_goal_complete": False,
        "nodes_digest": actual_digest,
        "expected_nodes_digest": expected_nodes_digest,
        "node_count": node_count,
        "graph": {
            "canonical_graph_authenticated": True,
            "operator_nodes_evaluated": False,
            "node_count": node_count,
            "sha256": actual_digest,
        },
        "rows": [row],
        "admitted_rows": [],
        "unresolved_rows": [row],
        "unresolved_gates": gates,
        "component_proofs": {
            "actual_math_coverage": witness.get("actual_math_coverage", {}),
            "conditional_formula_witness": witness.get("conditional_formula_witness", False),
            "c2_actual_source": witness.get("c2", {}),
            "whole_Z_actual_owner": witness.get("whole_z", {}),
            "correlated_source_query": witness.get("correlated", {}),
            "prescribed_C6_actual_micro": witness.get("c6", {}),
            "common_pressure_P0_sublemma": witness.get("pressure", {}),
            "common_core_projection_sublemma": witness.get("common_core_projection_sublemma", {}),
            "normalized_roots": witness.get("roots", {}),
            "source_receipts": witness.get("receipts", {}),
            "report_pins": witness.get("reports", {}),
        },
        "source_root_paths": witness.get("roots", {}).get("root_paths", {}),
        "preserved_source_handle_assertions": {
            "N": {"handle": 18, "is_parameter_value": False},
            "eta": {"handle": 21, "is_parameter_value": False},
            "d_star": {"handle": 22, "is_parameter_value": False},
            "handles_are_not_values": True,
            "cross_owner_constant_equality_proved": False,
            "phase_origin": "original first_micro phase; no new phase origin introduced",
        },
        "physical_factor_assertions": {
            "measure": "dy = hb * ds",
            "radius": "R = Ra * exp(hb*s)",
            "W": "F0 * sqrt(2*R) / Pstar",
            "FF": "F0^2 * R / Pstar^2",
            "actual_V_sign": "U = U0 + integral(hb * U_y) after preserving C6 double negative",
            "pressure": "pressure = P0 + (F0^2 * R / Pstar^2) * C",
        },
        "numeric_values_computed": False,
        "mathematical_operator_replayed": False,
        "old_operators_replayed": False,
        "operator_calls": 0,
        "producer_provider_builder_checker_replayed": False,
        "runtime_metadata": None,
        "record_hashes": dict(RECORD_HASHES),
        "decoded_record_hashes": dict(DECODED_RECORD_HASHES),
        "record_contract": {
            "plain": "source_text or raw_bytes, exact UTF-8 SHA-256",
            "gzip": "raw_bytes required; raw SHA-256 checked before decompression; decoded SHA checked after",
            "path": "experiments/root_st073/<basename>",
            "required_records": sorted(_REQUIRED_RECORDS),
            "optional_source_records": sorted(_FOUNDATION_RECORDS | _CORE_PROJECTION_RECORDS),
            "optional_saved_json_archives": sorted(ARCHIVE_RECORDS),
        },
        "proof_limitations": [
            "The bounded common analytic core projection uses the accepted finite-row/tail fixed-point receipt and exact rebuild AST consumer, then checks the N1 -> core-inlet -> C6 constructor chain; it does not prove first-micro Phi/U/own-history/F0 equality.",
            "The current whole-Z and prescribed-C6 providers have separate actual Phi/U/own-history/F0 paths; the independent normalized P0 subledger is admitted separately when its checked foundation records are supplied.",
            "Matching source family, names, class chain, receipt headers, or capacity flags is not a same-value function identity.",
            "The conditional normalized equations are source algebra only; they do not identify C2 values with C6 values.",
            "The full first-micro evolution, initial data, live parameter, own-moment/history, F0, and source-basis correspondence remains unresolved until a signed per-quantity shared branch is supplied.",
            "The actual generic B is 2*V_y; the prescribed loop B_over_Pstar = E*M/2 is a different quantity.",
            "No numerical roots, bounds, quadrature, recurrence, norm, producer, provider, builder, or checker was executed.",
        ],
    }


__all__ = [
    "DECODED_RECORD_HASHES",
    "EXPECTED_NODES_DIGEST",
    "EXPECTED_NODE_COUNT",
    "RECORD_HASHES",
    "SCHEMA",
    "SOURCE_FAMILY",
    "canonical_digest",
    "derive",
]
